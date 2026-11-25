#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_bench.py — 复现基准集管理、归因审计、置信度自评、阈值自适应推断、学术差异解释生成器

子命令
  score      给每篇复现打置信度分（0-1），列出 <0.6 的链路断裂点
  audit      失败模式聚类统计（成分漏检/靶点源/阈值/未披露/版本/功效）
  infer      阈值自适应推断：从原文参数文本建议 OB/DL/STRING/score 阈值
  explain    为不可复现案例生成可直接写入 Discussion 的学术差异解释
  full       一次性输出全部（JSON + 人类可读）

零第三方依赖。数据来自 data/benchmark_23.json（23 篇真实复现日志聚类）。
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BENCH = os.path.join(ROOT, "data", "benchmark_23.json")


def load_bench(path=BENCH):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def confidence(case, weights=None):
    """按 5 维加权求置信度；返回 (score, 最弱环节)。"""
    w = weights or {"compound": .15, "intersection": .15, "hub": .25, "pathway": .30, "docking": .15}
    dims = ["compound", "intersection", "hub", "pathway", "docking"]
    s = sum(case[d] * w[d] for d in dims)
    weakest = min(dims, key=lambda d: case[d])
    return round(s, 3), weakest


def cmd_score(args):
    b = load_bench()
    rows = []
    for c in b["cases"]:
        s, weak = confidence(c, b.get("weights"))
        rows.append((s, c["id"], c["key"], c["paper"], weak, c.get("breakpoint") or "—"))
    rows.sort()
    print("=== 23 篇复现置信度评分（0-1）===")
    print(f"{'分':>5} {'id':>3} {'key':<14} {'最弱环节':<14} 论文")
    for s, i, k, paper, weak, bp in rows:
        print(f"{s:5.2f} {i:3d} {k:<14} {weak:<14} {paper[:34]}")
    low = [r for r in rows if r[0] < 0.6]
    print(f"\n=== 低于 0.6 的案例：{len(low)}/{len(rows)} ===")
    for s, i, k, paper, weak, bp in low:
        print(f"[{k}] {s:.2f} 链路断裂于「{weak}」：{bp}")
    if args.json:
        out = [{"id": r[1], "key": r[2], "score": r[0], "weakest": r[4], "breakpoint": r[5]} for r in rows]
        json.dump(out, open(args.json, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print(f"\n已写 {args.json}")


def cmd_audit(args):
    b = load_bench()
    from collections import Counter, defaultdict
    cnt = Counter()
    by_mode = defaultdict(list)
    for c in b["cases"]:
        for fm in c.get("failure", []):
            cnt[fm] += 1
            by_mode[fm].append(c["key"])
    n = len(b["cases"])
    print("=== 23 篇失败模式聚类归因（可重叠）===")
    for fm, k in cnt.most_common():
        print(f"{k:2d}/{n} ({k/n:5.1%})  {fm}  {b['failure_modes'][fm]}")
        print(f"        案例: {', '.join(by_mode[fm])}")
    clean = [c['key'] for c in b['cases'] if not c.get('failure')]
    print(f"\n无显著失败（高保真）：{len(clean)}/{n}: {', '.join(clean)}")
    # 链路层级归因：把 failure 映射到阶段
    stage = Counter()
    mp = {"FM1_compound_gap": "02 成分层", "FM5_version_drift": "02 成分层",
          "FM2_target_source": "01/03 靶点层", "FM3_threshold": "参数层",
          "FM4_undisclosed": "原文报告层", "FM6_low_power": "03/06 统计功效层"}
    for c in b["cases"]:
        for fm in c.get("failure", []):
            stage[mp.get(fm, fm)] += 1
    print("\n=== 按流水线阶段归因 ===")
    for k, v in stage.most_common():
        print(f"  {k:<18} {v} 次")


# ---------------- 阈值自适应推断 ----------------
THRESH_PATTERNS = [
    ("ob",   [r'OB\s*[≥>=]+\s*(\d+)', r'生物利用度[^0-9]{0,6}(\d+)\s*%', r'oral bioavailability[^\d]{0,20}(\d+)'],
     "OB 阈值", "compounds.csv 的 ob 列"),
    ("dl",   [r'DL\s*[≥>=]+\s*(0?\.\d+)', r' drug[- ]?likeness[^\d]{0,20}(0?\.\d+)'],
     "DL 阈值", "compounds.csv 的 dl 列"),
    ("string",[r'STRING[^.]{0,40}?(0?\.\d{1,2})', r'置信度[^\d]{0,8}(0?\.\d{1,2})',
              r'confidence[^\d]{0,20}(0?\.\d{1,2})', r'minimum required interaction score[^\d]{0,20}(0?\.\d{1,2})'],
     "STRING 置信度", "string_edges.tsv 的 score 列"),
    ("genecards", [r'GeneCards[^.]{0,40}?(?:score|相关性)[^\d]{0,8}([≥>=]*\s*\d+(?:\.\d+)?)',
                   r'relevance score\s*[≥>=]+\s*(\d+(?:\.\d+)?)'],
     "GeneCards score", "疾病靶点表（本机用 Open Targets，需说明差异）"),
    ("pvalue", [r'p\s*<\s*(0?\.\d+)', r'P\s*值\s*<\s*(0?\.\d+)'],
     "富集 p 阈值", "enrich_significant.csv 的 pvalue/qvalue"),
]
DEFAULTS = {"ob": 30.0, "dl": 0.18, "string": 0.4, "genecards": None, "pvalue": 0.05}


def infer_thresholds(text):
    found = {}
    for key, pats, label, applies in THRESH_PATTERNS:
        vals = []
        for pat in pats:
            for m in re.finditer(pat, text, re.I):
                try:
                    vals.append(float(re.sub(r'[≥>=]|\s', '', m.group(1))))
                except ValueError:
                    pass
        found[key] = {"label": label, "values": sorted(set(vals)), "applies_to": applies,
                      "default": DEFAULTS[key],
                      "suggest": sorted(set(vals))[-1] if vals else DEFAULTS[key],
                      "source": "原文显式" if vals else "回退默认（原文未披露）"}
    return found


def cmd_infer(args):
    text = open(args.input, encoding="utf-8", errors="replace").read() if args.input else sys.stdin.read()
    r = infer_thresholds(text)
    print("=== 阈值自适应推断 ===")
    for k, v in r.items():
        vs = ", ".join(str(x) for x in v["values"]) or "（未见显式值）"
        print(f"{v['label']:<12} 原文值={vs:<18} 建议={v['suggest']!s:<8} 来源={v['source']}")
        print(f"{'':<12} 作用于 {v['applies_to']}")
    if args.json:
        json.dump(r, open(args.json, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print(f"已写 {args.json}")


# ---------------- 学术差异解释生成 ----------------
def academic_explanation(case, bench):
    c = case
    modes = [bench["failure_modes"][m] for m in c.get("failure", [])]
    weak = confidence(c, bench.get("weights"))[1]
    reasons = []
    if "FM1_compound_gap" in c.get("failure", []):
        reasons.append("本研究采用的成分-靶点数据库（TCMSP）对动物药、矿物药及高极性糖苷类成分的覆盖有限，"
                       "而原文化学成分谱来源于 UPLC-MS/MS 血清药物化学鉴定，二者在成分集合上存在系统性差异")
    if "FM2_target_source" in c.get("failure", []):
        reasons.append("原文靶点主要来自反向药效团/结构相似性预测服务器（如 SwissTargetPrediction、PharmMapper、Super-PRED），"
                       "而本研究采用可机器核查的实测活性数据库（ChEMBL/TCMSP）；预测型靶点集合通常较实测集合扩大 5–100 倍，"
                       "其 hub 排名更多反映预测算法对小分子骨架的偏好，而非真实的结合证据")
    if "FM5_version_drift" in c.get("failure", []):
        reasons.append("数据库版本更迭导致部分成分在现版数据库中已无法检索，属可溯源的版本漂移而非算法错误")
    if "FM3_threshold" in c.get("failure", []):
        reasons.append("原文部分筛选阈值（如 GeneCards 相关性中位数截断、靶点去重规则）未完整披露或无法复算，"
                       "本研究依可核查的默认值执行，可能对交集规模产生影响")
    if "FM4_undisclosed" in c.get("failure", []):
        reasons.append("原文分子对接未报告完整的 PDB 编号、结合口袋定义与结合能数值，"
                       "本研究对缺乏可用共晶配体的受体执行严格重对接门禁（RMSD<2 Å），因此不出分而非结果不一致")
    if "FM6_low_power" in c.get("failure", []):
        reasons.append("在交集样本量较小的情况下，超几何富集与网络拓扑分析的统计功效受限，"
                       "本研究对低于最小功效阈值的分析作拒绝处理，以避免过度解读")
    if not reasons:
        reasons.append("本研究在成分、交集、hub 与通路层均与原文高度一致，剩余差异主要来自疾病靶点库覆盖范围")
    body = "；".join(reasons) + "。"
    return (f"【{c['paper']}｜复现置信度 {confidence(c, bench.get('weights'))[0]:.2f}，"
            f"最弱环节 {weak}】\n"
            f"本研究与原文的差异可归因于以下客观因素：{body}\n"
            f"需强调，上述差异属于数据源性质与算法原理层面的差异，不构成对原文生物学结论的否定；"
            f"原文提供的实验验证（L3/L4 级证据）仍是机制判断的主要依据，而网络药理学结果应定位为假设生成。")


def cmd_explain(args):
    b = load_bench()
    if args.key:
        cases = [c for c in b["cases"] if c["key"] == args.key]
        if not cases:
            sys.exit(f"未找到 key={args.key}")
    else:
        cases = [c for c in b["cases"] if confidence(c, b.get("weights"))[0] < args.below]
    for c in cases:
        print(academic_explanation(c, b))
        print("\n" + "=" * 78 + "\n")


def cmd_full(args):
    b = load_bench()
    out = {"n": len(b["cases"]), "scores": [], "audit": {}}
    for c in b["cases"]:
        s, weak = confidence(c, b.get("weights"))
        out["scores"].append({"id": c["id"], "key": c["key"], "score": s, "weakest": weak,
                              "failure": c.get("failure", []), "explain": academic_explanation(c, b)})
    from collections import Counter
    cnt = Counter(fm for c in b["cases"] for fm in c.get("failure", []))
    out["audit"] = {fm: {"n": n, "label": b["failure_modes"][fm]} for fm, n in cnt.most_common()}
    out["mean_confidence"] = round(sum(x["score"] for x in out["scores"]) / len(out["scores"]), 3)
    path = args.json or os.path.join(ROOT, "data", "benchmark_report.json")
    json.dump(out, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"基准集报告 → {path}")
    print(f"平均置信度 {out['mean_confidence']}")
    print("失败模式：")
    for fm, v in out["audit"].items():
        print(f"  {v['n']:2d}/23  {fm}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("score"); s.add_argument("--json"); s.set_defaults(fn=cmd_score)
    a = sub.add_parser("audit"); a.set_defaults(fn=cmd_audit)
    i = sub.add_parser("infer"); i.add_argument("--input"); i.add_argument("--json"); i.set_defaults(fn=cmd_infer)
    e = sub.add_parser("explain"); e.add_argument("--key"); e.add_argument("--below", type=float, default=.6); e.set_defaults(fn=cmd_explain)
    f = sub.add_parser("full"); f.add_argument("--json"); f.set_defaults(fn=cmd_full)
    args = p.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
