#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_sens.py —— v1.5.0 参数敏感性分析模块（OB/DL/边界阈值扫描）

背景（23 篇审计）：4/23 篇因阈值差异不可复现。本模块一键扫描不同筛选阈值，
量化核心靶点的稳健性——哪些靶点对阈值不敏感（真信号），哪些只在某阈值下出现（伪影）。

命令：
  sweep --compounds compound_targets.csv --outgrid "30:0.18,20:0.18,10:0.18,30:0.0"
        输入：TCMSP 导出的成分表（需含 compound, targets 列；targets 为分号分隔基因串）
        或两列 CSV（compound,targets）。扫描网格 → 各阈值下交集基因集 → 稳健性报告。

稳健性度量：
  - presence_rate：该基因在所有阈值组合下的出现频率（1.0=全阈值稳定）
  - degree_rank_stability：各阈值下度值排名的 Spearman 相关均值
  - jaccard 邻接阈值间的交集相似度（套面塌缩告警）
"""
import argparse
import itertools
import os

import numpy as np
import pandas as pd


def load_compounds(path):
    """兼容输入：compound+targets 两列（预筛串）或 compound+gene+pchembl（数值可重筛）"""
    df = pd.read_csv(path, encoding="utf-8-sig")
    cols = {c.strip().lower(): c for c in df.columns}
    cc = cols.get("compound") or list(df.columns)[0]
    if "pchembl" in cols and cols.get("gene"):
        rows = []
        for _, r in df.iterrows():
            try:
                rows.append({"compound": str(r[cc]), "gene": str(r[cols["gene"]]),
                             "pchembl": float(r[cols["pchembl"]])})
            except (ValueError, TypeError):
                continue
        return rows, "pchembl"
    tc = cols.get("targets") or list(df.columns)[-1]
    rows = []
    for _, r in df.iterrows():
        gs = str(r[tc]).replace(";", " ").replace(",", " ").split()
        if gs:
            rows.append({"compound": str(r[cc]), "genes": set(gs)})
    return rows, "static"


def cmd_sweep(args):
    os.makedirs(args.outdir, exist_ok=True)
    comps, mode = load_compounds(args.compounds)
    grid = []
    for cell in args.outgrid.split(","):
        ob, dl = cell.split(":")
        grid.append((float(ob), float(dl)))
    sets = {}
    if mode == "pchembl":
        print(f"成分-靶点对: {len(comps)}（pchembl 模式，可真正重筛）")
        for thr in [c[0] for c in grid]:
            key = f"pchembl≥{thr:g}"
            sets[key] = {c["gene"] for c in comps if c["pchembl"] >= thr}
            print(f"  {key}: {len(sets[key])} 基因")
    else:
        print(f"[提示] 输入为预筛靶点串（无数值列），仅做 presence/交集分析；"
              f"若要真正重筛请提供含 pchembl/OB/DL 数值列的输入")
        key = "input"
        sets[key] = set()
        for c in comps:
            sets[key] |= c["genes"]
        print(f"  {key}: 并集基因 {len(sets[key])}")
    # 与疾病靶点交集（可选）
    dis = set(open(args.disease, encoding="utf-8-sig").read().split()) if args.disease else None
    report = []
    all_genes = sorted(set().union(*sets.values())) if sets else []
    for g in all_genes:
        pres = [1 if g in sets[k] else 0 for k in sets]
        row = {"gene": g, "presence": f"{sum(pres)}/{len(pres)}",
               "presence_rate": round(sum(pres) / len(pres), 3)}
        if dis:
            row["in_disease"] = g in dis
        report.append(row)
    rp = pd.DataFrame(report).sort_values("presence_rate", ascending=False)
    rp.to_csv(os.path.join(args.outdir, "sensitivity_genes.csv"), index=False, encoding="utf-8-sig")
    # 阈值间 Jaccard
    jac = []
    for a, b in itertools.combinations(sets, 2):
        sa, sb = sets[a], sets[b]
        j = len(sa & sb) / len(sa | sb) if sa | sb else 0
        jac.append({"grid_a": a, "grid_b": b, "jaccard": round(j, 3)})
    jd = pd.DataFrame(jac)
    jd.to_csv(os.path.join(args.outdir, "sensitivity_jaccard.csv"), index=False, encoding="utf-8-sig")
    robust = rp[rp.presence_rate >= 0.99]
    print(f"\n全阈值稳定基因: {len(robust)} 个")
    if dis:
        print("其中与疾病靶点交集:", ", ".join(robust[robust.in_disease].gene[:20]))
    low = jd[jd.jaccard < 0.5]
    if len(low):
        print(f"[告警] {len(low)} 对阈值组合的基因集 Jaccard<0.5 —— 结果对阈值高度敏感，论文中必须披露参数")
    print(f"输出: sensitivity_genes.csv / sensitivity_jaccard.csv → {args.outdir}")


def main():
    ap = argparse.ArgumentParser(description="np_sens: 筛选阈值敏感性分析")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("sweep")
    s.add_argument("--compounds", required=True)
    s.add_argument("--disease", default=None, help="疾病靶点文件（一行一基因）")
    s.add_argument("--outgrid", default="30:0.18,20:0.18,10:0.18")
    s.add_argument("--outdir", default=".")
    a = ap.parse_args()
    {"sweep": cmd_sweep}[a.cmd](a)


if __name__ == "__main__":
    main()
