#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_guard.py — 复现合理性边界检查与异常告警

当复现出的数量与原文偏差超阈值时自动告警并提示可能原因：
  - 成分-靶点对数量偏差 >50%
  - PPI 网络节点数 < 原文的 70%
  - 交集数偏差 >50%
  - 靶点映射率 <80%
  - 成分数偏差 >50%
  - 无任何核心 hub 命中（命中率 0）
  - 交集 <5（统计功效不足）

用法
  python np_guard.py --local local.json --paper paper.json [--outdir .]
  python np_guard.py --demo     # 用内置 23 篇基准的阈值演示
"""
from __future__ import annotations

import argparse, json, os, sys

RULES = [
    ("n_compounds", 0.5, "成分数偏差>50%", "成分漏检：database coverage；非植物药/矿物药/糖苷；OB/DL 过严"),
    ("n_ct_pairs", 0.5, "成分-靶点对偏差>50%", "靶点源差异（预测vs实测）；靶点映射丢失；成分未匹配"),
    ("n_intersection", 0.5, "交集数偏差>50%", "疾病靶点库规模差异；GeneCards score 阈值；物种限定"),
    ("n_ppi_nodes", 0.3, "PPI节点数<原文70%", "STRING 置信度过高；交集基因过多孤立；阈值过严"),
    ("mapping_rate", 0.2, "靶点映射率<80%", "UniProt 词典版本；同义词未覆盖；物种ID错误"),
    ("hub_hit_rate", 0.5, "hub命中率<50%", "靶点源差异导致 hub 重排；原文 hub 为预测工具偏好"),
]


def check(local, paper):
    alerts = []
    for key, tol, label, hint in RULES:
        if key not in local or key not in paper:
            continue
        l, p = float(local[key]), float(paper[key])
        if p <= 0:
            continue
        if key == "n_ppi_nodes":
            bad = l < 0.7 * p
            detail = f"本机 {l:.0f} vs 原文 {p:.0f}（阈值 70%）"
        elif key == "mapping_rate":
            bad = l < p - tol if p < 1 else l < 0.8
            detail = f"本机 {l:.3f}"
        else:
            bad = abs(l - p) / p > tol
            detail = f"本机 {l:.0f} vs 原文 {p:.0f}（偏差 {abs(l-p)/p:.0%}）"
        if bad:
            alerts.append({"rule": key, "level": "WARN", "label": label,
                           "detail": detail, "possible_causes": hint})
    if local.get("n_intersection", 999) < 5:
        alerts.append({"rule": "low_power", "level": "BLOCK", "label": "交集<5，统计功效不足",
                       "detail": f"交集={local['n_intersection']}",
                       "possible_causes": "改走通路/文献路线，不做 hub 分析"})
    return alerts


def cmd_run(args):
    local = json.load(open(args.local, encoding="utf-8"))
    paper = json.load(open(args.paper, encoding="utf-8"))
    alerts = check(local, paper)
    print(f"=== 复现合理性检查: {len(alerts)} 条告警 ===")
    for a in alerts:
        print(f"[{a['level']}] {a['label']} — {a['detail']}")
        print(f"        可能原因：{a['possible_causes']}")
    if not alerts:
        print("未触发告警；各数量级与原文一致。")
    if args.output:
        json.dump(alerts, open(args.output, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print(f"已写 {args.output}")
    return 1 if any(a["level"] == "BLOCK" for a in alerts) else 0


def cmd_demo(args):
    print("=== 阈值示例（用 23 篇基准的真实案例）===")
    demo = [
        ({"n_compounds": 1, "n_ct_pairs": 27, "n_intersection": 17, "n_ppi_nodes": 15, "hub_hit_rate": 0.0},
         {"n_compounds": 1, "n_ct_pairs": 124, "n_intersection": 89, "n_ppi_nodes": 90, "hub_hit_rate": 1.0}, "celastrol"),
        ({"n_compounds": 69, "n_ct_pairs": 584, "n_intersection": 132, "n_ppi_nodes": 122, "hub_hit_rate": 0.7},
         {"n_compounds": 69, "n_ct_pairs": 600, "n_intersection": 180, "n_ppi_nodes": 180, "hub_hit_rate": 1.0}, "liuwei"),
        ({"n_compounds": 1, "n_ct_pairs": 1, "n_intersection": 1, "n_ppi_nodes": 0},
         {"n_compounds": 1, "n_ct_pairs": 300, "n_intersection": 100, "n_ppi_nodes": 90}, "crocetin"),
    ]
    for l, p, name in demo:
        al = check(l, p)
        print(f"\n[{name}] {len(al)} 条告警")
        for a in al:
            print(f"  [{a['level']}] {a['label']}: {a['detail']}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("--local", required=True); r.add_argument("--paper", required=True)
    r.add_argument("--output"); r.set_defaults(fn=cmd_run)
    d = sub.add_parser("demo"); d.set_defaults(fn=cmd_demo)
    args = p.parse_args(); sys.exit(args.fn(args) or 0)


if __name__ == "__main__":
    main()
