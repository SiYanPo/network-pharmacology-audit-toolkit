#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_qa.py —— v1.5.0 一键复现质量报告（QA Report）

背景（23 篇痛点）：复现收尾时需人工拼装"数据版本+参数+相似度+薄弱环节"——
本模块一键扫描复现输出目录，生成可附论文投稿材料包的 QA 报告。

命令：
  report --dir 复现输出目录 [--paper 原文参数.json] [--outdir DIR]
扫描内容：
  - 产物清单与完整性（按 np_export 五表/富集/对接产物分类）
  - 参数与数据源（environment.json / run_manifest.json / ligands / gmt 文件时间戳）
  - 与原文数量对比（--paper 时：交集基因数/PPI 节点边数/对接对数）
  - 审稿人薄弱环节清单（基于 np_bench 置信度规则 + np_guard 告警规则内联实现）
"""
import argparse
import hashlib
import json
import os
from datetime import datetime

import pandas as pd

EXPECTED = {
    "herbs": ["herbs.csv", "herbs.json"],
    "compounds": ["compounds.csv", "compounds_kept.csv", "compounds_standard.csv"],
    "targets": ["targets_mapped.csv", "intersection_genes.txt", "disease_targets.tsv"],
    "ppi": ["intersect_edges_0.4.tsv", "ppi_edges_clean.tsv", "hubs.csv", "topology.csv", "network.svg"],
    "enrich": ["enrich_significant.csv", "enrich_all.csv", "enrich_bubble.svg"],
    "docking": ["docking_scores.csv", "docking_validation.csv", "poses"],
    "quality": ["provenance.json", "audit_report.md", "admet_metrics.csv"],
}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def cmd_report(args):
    os.makedirs(args.outdir, exist_ok=True)
    d = args.dir
    files = []
    for root, _, fns in os.walk(d):
        for fn in fns:
            p = os.path.join(root, fn)
            rel = os.path.relpath(p, d)
            try:
                sz = os.path.getsize(p)
                mt = datetime.fromtimestamp(os.path.getmtime(p)).strftime("%Y-%m-%d %H:%M")
                files.append({"file": rel, "size_kb": round(sz / 1024, 1), "mtime": mt,
                              "sha256_16": sha256(p) if sz < 5 * 1024 * 1024 else "(>5MB)"})
            except OSError:
                pass
    fl = pd.DataFrame(files)
    # 完整性
    found = {f for f in fl.file}
    completeness = {}
    for cat, cands in EXPECTED.items():
        hit = [c for c in cands if any(f.endswith(c) for f in found)]
        completeness[cat] = ("OK" if hit else "MISSING")
    # environment/manifest
    env = {}
    for cand in ("environment.json", "run_manifest.json", "provenance.json"):
        p = os.path.join(d, cand)
        if os.path.exists(p):
            try:
                env[cand] = json.load(open(p, encoding="utf-8"))
            except Exception:
                env[cand] = {}
    # 数量对比
    compare = []
    paper = {}
    if args.paper and os.path.exists(args.paper):
        paper = json.load(open(args.paper, encoding="utf-8-sig"))
    def count_of(fn, col=None):
        p = os.path.join(d, fn)
        if not os.path.exists(p):
            return None
        try:
            df = pd.read_csv(p, encoding="utf-8-sig")
            return len(df)
        except Exception:
            try:
                return len(open(p, encoding="utf-8-sig").read().split())
            except Exception:
                return None
    for key, fn in [("intersection_genes", "intersection_genes.txt"),
                    ("compounds", "compounds.csv"), ("compounds_kept", "compounds_kept.csv"),
                    ("ppi_edges", "ppi_edges_clean.tsv"), ("enrich_significant", "enrich_significant.csv")]:
        local = count_of(fn)
        paper_v = paper.get(key)
        if local is not None or paper_v is not None:
            ratio = f"{local/paper_v:.0%}" if local and paper_v else "—"
            compare.append({"item": key, "local": local, "paper": paper_v, "ratio": ratio,
                            "flag": "OK" if (paper_v and local and 0.7 <= local / paper_v <= 1.3) else
                                    ("WARN" if local and paper_v else "—")})
    # 薄弱环节（规则内联自 np_guard/np_bench）
    weak = []
    if completeness["herbs"] == "MISSING":
        weak.append("药材清单缺失 → 数据溯源断链（审稿人会问数据来自哪个库哪个日期）")
    if completeness["quality"] == "MISSING":
        weak.append("无 provenance/快照 → 无法精确归因版本漂移")
    for c in compare:
        if c["flag"] == "WARN":
            weak.append(f"{c['item']} 本机 {c['local']} vs 原文 {c['paper']}（偏差>30%）→ 需在 Discussion 归因")
    if env.get("environment.json", {}).get("python_version"):
        pass
    if not weak:
        weak.append("未检出结构性薄弱环节——仍建议人工核查对接受体选择与富集库版本")
    # 报告
    md = ["# 复现质量报告（QA Report）", "",
          f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M')} | 目录: `{os.path.abspath(d)}`",
          f"产物文件: {len(files)} 个 | SHA256 指纹已记录（≤5MB 产物）", "",
          "## 产物完整性", "", "| 类别 | 状态 |", "|---|---|"]
    md += [f"| {k} | {v} |" for k, v in completeness.items()]
    if compare:
        md += ["", "## 与原文数量对比", "", "| 项目 | 本机 | 原文 | 比例 | 判定 |", "|---|---|---|---|---|"]
        md += [f"| {c['item']} | {c['local']} | {c['paper']} | {c['ratio']} | {c['flag']} |" for c in compare]
    if env:
        md += ["", "## 参数与环境快照", "", "```json",
               json.dumps(env, ensure_ascii=False, indent=1)[:1200], "```"]
    md += ["", "## 审稿人薄弱环节与应对", ""]
    md += [f"- {w}" for w in weak]
    md += ["", "## 产物 SHA256（抽查前 12 个）", "", "| 文件 | SHA256(16位) |", "|---|---|"]
    md += [f"| {r.file} | {r.sha256_16} |" for _, r in fl.head(12).iterrows()]
    out = os.path.join(args.outdir, "qa_report.md")
    open(out, "w", encoding="utf-8").write("\n".join(md))
    fl.to_csv(os.path.join(args.outdir, "qa_files.csv"), index=False, encoding="utf-8-sig")
    print(f"产物 {len(files)} 文件 | 完整性: {sum(1 for v in completeness.values() if v=='OK')}/{len(completed_ok(completeness))} 类命中")
    print(f"薄弱环节 {len(weak)} 条")
    print(f"报告 → {out}")


def completed_ok(c):
    return [k for k, v in c.items()]


def main():
    ap = argparse.ArgumentParser(description="np_qa: 一键复现质量报告")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("report")
    s.add_argument("--dir", required=True)
    s.add_argument("--paper", default=None)
    s.add_argument("--outdir", default=".")
    a = ap.parse_args()
    {"report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
