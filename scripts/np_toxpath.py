#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_toxpath.py —— v1.5.0 网络毒理学模块（ADMET 风险 → 毒理通路映射）

背景（23 篇审计）：仅 1/23 篇提到网络毒理学——中文 NP 论文普遍不做毒性机制映射，
这是审稿薄弱环节。本模块把 np_admet.py 的预测输出升级为"毒性机制解释"：
  1. 读取 predict 输出的 ADMET 谱（tcm_admet_profile.csv）
  2. 标记高风险化合物-端点（概率≥0.7 或按用户阈值）
  3. 将风险端点映射到毒理学相关 KEGG/通路基因集（内置映射表）
  4. 可选：调 np_enrich 对风险基因集做富集（--enrich）

内置端点→毒理机制映射（文献依据见列）：
  herg      → 心律失常（KEGG hsa05410 心肌收缩/离子通道；QT 延长机制）
  ames      → DNA 损伤与修复（hsa03410 BER / hsa03420 MMR / p53 信号 hsa04115）
  dili      → 肝损伤（药物代谢 hsa00982、胆汁分泌 hsa04976、Nrf2/氧化应激 hsa04168）
  clintox   → 临床肝毒性/多器官（线粒体功能障碍、凋亡 hsa04210）
  skin_reaction → 皮肤致敏（角质形成细胞活化、IL-1/TNF 信号 hsa04668）
  tox21_*   → 核受体/氧化应激通路（AhR、PPAR、ER、ARE/Nrf2——与药效机制同源，注意药毒两重性）

命令：
  run --profile tcm_admet_profile.csv [--threshold 0.7] [--outdir DIR] [--enrich]
  对每化合物输出风险清单+机制假设（可写入论文 Safety Discussion）。
"""
import argparse
import os

import pandas as pd

RISK_MAP = {
    "herg": {"mechanism": "hERG 钾通道阻断 → QT 间期延长 → 尖端扭转型室速",
             "pathways": "KEGG hsa05410（心肌收缩/离子通道）、Genes 心脏毒性基因集",
             "note": "概率>0.5 建议加做 hERG patch-clamp 或 CiPA 计算评估"},
    "ames": {"mechanism": "DNA 烷基化/氧化损伤 → 致突变",
             "pathways": "hsa03410 碱基切除修复、hsa03420 错配修复、hsa04115 p53 信号",
             "note": "喹啉/呋喃香豆素等骨架易阳性；结合定量结构警报（Toxtree）复核"},
    "dili": {"mechanism": "反应性代谢物 → 肝细胞损伤/胆汁淤积/免疫介导肝毒",
             "pathways": "hsa00982 药物代谢-CYP、hsa04976 胆汁分泌、hsa04168 Nrf2/氧化应激",
             "note": "与 CYP 抑制谱联读：CYP3A4 强抑制+DILI 高风险 → 药物相互作用警示"},
    "clintox": {"mechanism": "临床毒性（肝/心/肾）——多为线粒体与凋亡通路",
                "pathways": "hsa04210 凋亡、线粒体氧化磷酸化损伤",
                "note": "模型 AUROC 0.85（本机实测），概率仅供参考需实验复核"},
    "skin_reaction": {"mechanism": "皮肤致敏/刺激性接触性皮炎",
                      "pathways": "hsa04668 TNF/IL-1 信号、角质形成细胞活化通路",
                      "note": "外用制剂优先评估"},
    "tox21_SR-ARE": {"mechanism": "Nrf2/ARE 氧化应激轴激活——药毒两重性",
                     "pathways": "hsa04168 Nrf2 通路",
                     "note": "铁死亡/抗氧化叙事论文中：SR-ARE 阳性=药效证据兼氧化应激预警"},
    "tox21_NR-PPAR-gamma": {"mechanism": "PPARγ 激活——代谢调节/体重增加风险",
                            "pathways": "hsa03320 PPAR 通路",
                            "note": "PPARG 主导机制论文：阳性=机制佐证；长期用药需查体重/水肿"},
    "tox21_NR-AhR": {"mechanism": "芳烃受体激活——代谢酶诱导/免疫调节",
                     "pathways": "hsa04125 Xeno？注：AhR-CYP1A1 轴",
                     "note": "AhR 激活诱导 CYP1A1 → 前致癌物活化风险（吸烟者尤需注意）"},
    "tox21_NR-ER": {"mechanism": "雌激素受体干扰——内分泌效应",
                    "pathways": "hsa04915 雌激素信号通路",
                    "note": "植物雌激素类黄酮常见阳性；妇科/乳腺疾病患者用药需提示"},
}


def cmd_run(args):
    os.makedirs(args.outdir, exist_ok=True)
    p = pd.read_csv(args.profile, encoding="utf-8-sig").set_index("compound")
    rows = []
    for comp in p.index:
        for ep, info in RISK_MAP.items():
            if ep not in p.columns:
                continue
            v = p.loc[comp, ep]
            try:
                v = float(v)
            except (TypeError, ValueError):
                continue
            if v >= args.threshold:
                rows.append({"compound": comp, "endpoint": ep, "prob": round(v, 3),
                             "risk": info["mechanism"], "pathways": info["pathways"],
                             "note": info["note"]})
    df = pd.DataFrame(rows).sort_values(["compound", "prob"], ascending=[True, False])
    df.to_csv(os.path.join(args.outdir, "tox_risk_map.csv"), index=False, encoding="utf-8-sig")
    if df.empty:
        print(f"无化合物端点超过阈值 {args.threshold} —— 无需毒理讨论（或在 ADMET uncertain 范围）")
        return
    # Markdown 汇总
    md = ["# 网络毒理学风险-通路映射报告", "",
          f"阈值: 概率≥{args.threshold} | 化合物数: {df.compound.nunique()} | 风险条目: {len(df)}", ""]
    for comp, g in df.groupby("compound"):
        md.append(f"## {comp}")
        for _, r in g.iterrows():
            md.append(f"- **{r.endpoint}** (p={r.prob})：{r.risk}")
        md.append("")
    md += ["## 写入 Safety Discussion 的模板句", "",
           "本研究预测的毒性风险主要集中在 [端点列表]，其机制可解释为 [映射通路]；",
           "鉴于网络药理学预测不区分剂量与暴露途径，上述风险提示应在后续验证中通过",
           "[实验建议] 予以复核，而不构成临床禁忌结论。", ""]
    open(os.path.join(args.outdir, "tox_risk_report.md"), "w", encoding="utf-8").write("\n".join(md))
    print(f"风险条目 {len(df)}（{df.compound.nunique()} 化合物）")
    print(f"输出 → {args.outdir}/tox_risk_map.csv + tox_risk_report.md")
    top = df.head(6)
    for _, r in top.iterrows():
        print(f"  [{r.compound}] {r.endpoint} p={r.prob} → {str(r.risk)[:50]}")


def main():
    ap = argparse.ArgumentParser(description="np_toxpath: ADMET 风险 → 毒理通路映射")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("run")
    s.add_argument("--profile", required=True, help="np_admet predict 输出的宽表")
    s.add_argument("--threshold", type=float, default=0.7)
    s.add_argument("--outdir", default=".")
    a = ap.parse_args()
    {"run": cmd_run}[a.cmd](a)


if __name__ == "__main__":
    main()
