#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_lab.py —— v1.5.0 反向湿实验指导模块（Hub 基因 → 实验方案骨架）

背景（23 篇审计）：18/23 篇用 WB/qPCR 验证、6/23 用动物实验——湿实验是中文 NP 论文
标配收尾。本模块对核心 hub 基因生成结构化实验方案骨架：
  - 基因功能摘要（UniProt REST 实时拉取，不编造）
  - qPCR 引物设计参数（exon-spanning/产物长度）+ PrimerBank 检索链接（不给臆造序列）
  - 抗体检索链接（CiteAb，聚合文献引用数——不编造 catalog 号）
  - 疾病匹配的细胞系/动物模型推荐（内置 20 病种映射表，来自文献常用实践）

诚实边界：不生成具体抗体货号与引物序列（需实时数据库+校验，编造=幻觉风险），
全部给出权威检索入口；细胞系/模型为"文献常用实践"而非唯一推荐。

命令：
  plan --genes STAT3,PTGS2 --disease 胃癌 --outdir DIR
"""
import argparse
import json
import os

import pandas as pd

# 病种 → 细胞系/动物模型（文献常用实践）
DISEASE_MODELS = {
 "胃癌": {"lines": ["MKN-45", "AGS", "HGC-27", "BGC-823"], "animals": ["C57BL/6 原位移植瘤", "BALB/c-nu 荷瘤鼠", "MNNG 诱导胃癌模型"]},
 "肝癌": {"lines": ["HepG2", "Huh7", "HCCLM3", "SMMC-7721"], "animals": ["DEN 诱导 HCC 小鼠", "C57BL/6 原位肝接种"]},
 "结直肠癌": {"lines": ["HCT116", "SW480", "LoVo", "DLD-1"], "animals": ["AOM/DSS 结肠炎相关癌变模型", "CT-26 荷瘤 BALB/c"]},
 "糖尿病肾病": {"lines": ["HK-2", "SV40 MES13", "NRK-52E"], "animals": ["db/db 小鼠", "STZ 诱导 DKD 大鼠", "5/6 肾切除"]},
 "慢性肾病": {"lines": ["HK-2", "NRK-52E", "RPTEC"], "animals": ["腺嘌呤诱导 CKD 大鼠", "5/6 肾切除", "UUO 纤维化模型"]},
 "阿尔茨海默病": {"lines": ["SH-SY5Y", "PC12", "APP/PSEN1 转 Aβ 细胞模型"], "animals": ["APP/PS1 双转基因小鼠", "SAMP8 衰老加速鼠", "D-半乳糖+AlCl3 模型"]},
 "椎间盘退变": {"lines": ["大鼠/人 NP 细胞（原代）", "MSC 分化 NP 样细胞"], "animals": ["针刺尾椎 IDD 大鼠", "自然衰老沙鼠"]},
 "抑郁": {"lines": ["CORT 处理 SH-SY5Y", "原代海马神经元"], "animals": ["CUMS 慢性不可预知应激", "CSDS 社会挫败", "OLC 颊黏膜刺激"]},
 "骨质疏松": {"lines": ["MC3T3-E1", "RAW264.7 破骨分化"], "animals": ["OVX 去卵巢小鼠/大鼠", "糖皮质激素性 OP"]},
 "动脉粥样硬化": {"lines": ["THP-1 巨噬细胞", "HUVEC", "VSMC"], "animals": ["ApoE-/- 高脂小鼠", "LDLR-/- 小鼠"]},
 "肺纤维化": {"lines": ["A549", "MRC-5", "原代成纤维细胞"], "animals": ["BLM 博来霉素诱导 PF 小鼠"]},
 "银屑病": {"lines": ["HaCaT", "NHEK"], "animals": ["IMQ 咪喹莫特模型"]},
 "类风湿关节炎": {"lines": ["RA-FLS", "RAW264.7", "MH7A"], "animals": ["CIA 胶原诱导关节炎", "AIA 佐剂关节炎"]},
 "心肌缺血": {"lines": ["H9c2", "AC16"], "animals": ["LAD 结扎 MI 小鼠", "I/R 缺血再灌注模型"]},
 "非酒精性脂肪肝": {"lines": ["HepG2 + FFA", "AML-12"], "animals": ["HFD 高脂饮食", "MCD 饮食", "C57BL/6 Lep ob/ob"]},
 "肺腺癌": {"lines": ["A549", "PC-9", "H1975"], "animals": ["LLC 皮下/原位", "PDX 模型"]},
 "乳腺癌": {"lines": ["MCF-7", "MDA-MB-231", "4T1"], "animals": ["4T1 荷瘤 BALB/c", "MDA-MB-231 骨转移模型"]},
 "前列腺癌": {"lines": ["PC-3", "LNCaP", "DU145"], "animals": ["PC-3 荷瘤裸鼠", "TRAMP 转基因"]},
 "胰岛素抵抗": {"lines": ["C2C12 + PA", "HepG2 + PA", "3T3-L1 脂肪细胞"], "animals": ["HFD+STZ", "db/db"]},
 "肠道菌群-肠病": {"lines": ["Caco-2", "NCM460", "HT-29"], "animals": ["DSS 结肠炎", "无菌鼠粪菌移植（FMT）"]},
}


def _uniprot(gene):
    import requests
    try:
        r = requests.get("https://rest.uniprot.org/uniprotkb/search",
                         params={"query": f'gene_exact:"{gene}" AND organism_id:9606 AND reviewed:true',
                                 "format": "json", "size": 1, "fields": "accession,function_name,protein_name"},
                         timeout=60)
        res = r.json().get("results") or []
        if res:
            f = (res[0].get("comments") or [{}])
            fn = next((c["texts"][0]["value"] for c in f if c.get("commentType") == "FUNCTION"), "")
            return {"accession": res[0]["primaryAccession"], "protein": res[0].get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""), "function": (fn or "")[:220]}
    except Exception:
        pass
    return None


def cmd_plan(args):
    os.makedirs(args.outdir, exist_ok=True)
    genes = [g.strip().upper() for g in args.genes.split(",") if g.strip()]
    dis = args.disease
    models = DISEASE_MODELS.get(dis, {})
    if not models:
        print(f"[提示] 病种库无 '{dis}'，可用病种: {', '.join(list(DISEASE_MODELS)[:10])}...")
    rows = []
    for g in genes:
        u = _uniprot(g)
        rows.append({
            "gene": g,
            "uniprot": u["accession"] if u else "",
            "protein": u["protein"] if u else "",
            "function_摘要": u["function"] if u else "（UniProt 查询失败，请核查基因名）",
            "qPCR_引物": f"PrimerBank 检索 https://pga.hs.harvard.edu/primerbank/ （gene={g}）；"
                        f"设计原则：跨外显子（exon-exon junction）、产物 80-150 bp、Tm 58-62℃、GC 40-60%",
            "抗体": f"CiteAb 检索 https://citeab.com/search?q={g}+antibody （按引用数排序选 WB/IF 验证级）；"
                   f"另参考 Human Protein Atlas https://www.proteinatlas.org/search/{g} （组织定位+抗体验证级别）",
            "WB_对照": "β-actin/GAPDH 内参；如靶点为磷酸化通路（STAT3/AKT1）加 total 蛋白对照与刺激/抑制对照",
        })
    pd.DataFrame(rows).to_csv(os.path.join(args.outdir, "hub_gene_lab_plan.csv"),
                              index=False, encoding="utf-8-sig")
    md = ["# 湿实验方案骨架（反向验证指导）", "",
          f"病种: **{dis}** | Hub 基因: **{', '.join(genes)}**", ""]
    if models:
        md += ["## 细胞模型（文献常用）", "", "\n".join(f"- {l}" for l in models["lines"]), "",
               "## 动物模型（文献常用）", "", "\n".join(f"- {a}" for a in models["animals"]), ""]
    md += ["## 推荐验证组合（依 23 篇文献套路频次）", "",
           "1. **WB/qPCR**（18/23 篇采用）：核心 hub 蛋白表达+磷酸化水平",
           "2. **免疫荧光/组化**：定位验证（如 Nrf2 核转位、p-AKT 膜定位）",
           "3. **功能回复实验**：siRNA/过表达反转（对应 NP 预测的抑制/激活方向）",
           "4. **动物实验**（6/23 篇采用）：模型见上表；给药剂量按临床等效剂量换算（体表面积法）",
           "5. **加分项**（仅 antrodia 一篇用过）：TCGA 生存分析+GEO 验证 → 可用 np_survival.py 直接产出",
           "",
           "## 诚实边界",
           "",
           "- 抗体货号/引物序列不自动生成（实时校验数据库+交叉验证成本高，臆造=幻觉）；",
           "- 上述细胞系/模型为文献常用实践汇总，不构成唯一推荐；实验前请核对细胞系 STR 鉴定与文献最近复现情况。",
           ""]
    open(os.path.join(args.outdir, "lab_plan.md"), "w", encoding="utf-8").write("\n".join(md))
    print(f"方案骨架 → {args.outdir}/lab_plan.md + hub_gene_lab_plan.csv")


def main():
    ap = argparse.ArgumentParser(description="np_lab: 湿实验方案生成器")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("plan")
    s.add_argument("--genes", required=True)
    s.add_argument("--disease", required=True)
    s.add_argument("--outdir", default=".")
    a = ap.parse_args()
    {"plan": cmd_plan}[a.cmd](a)


if __name__ == "__main__":
    main()
