#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""np_wetlab_result.py — 湿实验结果回灌：把实验结果变成可参与证据分级的记录。

背景（SCM0061-2021 缺口 G-3）
------------------------------
标准 6.1.3 / 6.3.3 明确把"实验验证"列为结果验证的一类，并要求给出**验证模型**与**检测指标**。
本 skill 此前只有实验**计划**（np_wetlab.py 的引物特异性、亚细胞定位），没有结果回灌的入口，
导致湿实验证据无法自动进入证据账本——这是"标准条款已实现"表里唯一的空白之一。

本模块做三件事：
1. 按 `12_wetlab/wetlab_result.csv` 契约校验表结构（交给 np_contracts）；
2. 按预注册规则给每条结果定级：
     L3（实验支持）要求 n_replicates >= 3 且 p_value < alpha（默认 0.05）且给出 readout；
     n=1 或 p 值缺失 → 降为 L2 提示（"有实验但不足以支撑"）；
3. 输出 `wetlab_evidence.json`，字段与 np_evidence_ledger 的 L3 条目兼容。

用法：
  python scripts/np_wetlab_result.py --csv 12_wetlab/wetlab_result.csv [--alpha 0.05]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from np_contracts import CONTRACTS, check_file
except Exception:  # 允许独立运行
    CONTRACTS, check_file = {}, None

REQUIRED = ["experiment", "model", "readout", "value", "n_replicates"]


def _f(row: dict, *names: str) -> float | None:
    for n in names:
        if n in row and str(row[n]).strip() not in ("", "NA", "na", "-"):
            try:
                return float(str(row[n]).strip())
            except ValueError:
                continue
    return None


def _get(row: dict, *names: str) -> str:
    for n in names:
        if n in row and str(row[n]).strip():
            return str(row[n]).strip()
    return ""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="np_wetlab_result: 湿实验结果回灌与定级")
    ap.add_argument("--csv", required=True)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args(argv)

    if not os.path.exists(a.csv):
        raise SystemExit(f"[FAIL] 找不到 {a.csv}")
    spec = CONTRACTS.get("12_wetlab/wetlab_result.csv")
    if check_file and spec:
        ok, probs = check_file(a.csv, spec)
        if not ok:
            print("[FAIL] 契约不通过：")
            for p in probs:
                print("   -", p)
            return 2

    rows = list(csv.DictReader(open(a.csv, encoding="utf-8-sig")))
    recs, n_l3, n_l2 = [], 0, 0
    for r in rows:
        n = _f(r, "n_replicates", "n", "replicates", "重复数")
        p = _f(r, "p_value", "p", "pval", "P")
        readout = _get(r, "readout", "indicator", "metric", "检测指标")
        model = _get(r, "model", "cell_line", "system", "模型")
        rec = {
            "experiment": _get(r, "experiment", "assay", "实验"),
            "model": model, "readout": readout,
            "value": _f(r, "value"), "n_replicates": int(n) if n else None,
            "p_value": p, "direction": _get(r, "direction", "effect", "方向"),
        }
        reasons = []
        if not model:
            reasons.append("缺验证模型（SCM0061 6.1.3 要求）")
        if not readout:
            reasons.append("缺检测指标（SCM0061 6.3.3 要求）")
        if n is None or n < a.min_n:
            reasons.append(f"重复数 {int(n) if n else 0} < {a.min_n}，不足以定 L3")
        if p is None:
            reasons.append("缺 p 值，无法判定显著性")
        elif p >= a.alpha:
            reasons.append(f"p={p:.3g} ≥ α={a.alpha}，不显著")

        if not reasons:
            rec["evidence_level"] = "L3"
            rec["rationale"] = f"n={int(n)}≥{a.min_n} 且 p={p:.3g}<{a.alpha}，模型与指标齐备"
            n_l3 += 1
        else:
            rec["evidence_level"] = "L2"
            rec["rationale"] = "；".join(reasons)
            n_l2 += 1
        recs.append(rec)

    out = {
        "source_csv": os.path.abspath(a.csv),
        "rules": {"min_replicates": a.min_n, "alpha": a.alpha,
                  "note": "n<min 或 p 缺失/不显著 → 只能作 L2 提示，不得计入 L3 实验支持"},
        "n_records": len(recs), "n_L3": n_l3, "n_L2": n_l2,
        "records": recs,
        "scm0061_mapping": "6.1.3 结果验证（方法/模型/终点）、6.3.3 结果验证（模型/检测指标）",
    }
    outdir = a.outdir or os.path.dirname(os.path.abspath(a.csv))
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, "wetlab_evidence.json")
    json.dump(out, open(p, "w", encoding="utf-8"), ensure_ascii=False, indent=1)

    print(f"湿实验证据：{len(recs)} 条 → L3 {n_l3} 条 / L2 {n_l2} 条")
    for r in recs:
        print(f"  [{r['evidence_level']}] {r['experiment'][:28]:<28} {r['readout'][:18]:<18} {r['rationale'][:60]}")
    print(f"→ {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
