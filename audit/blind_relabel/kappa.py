#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""kappa.py — 用第二位标注者的盲法结果计算评分者间信度。

用法：
    python kappa.py blind_sheet.csv ground_truth.csv [--min-conf 1]

输入
    blind_sheet.csv   标注者填好的表（需含 pmid, rater_label, rater_confidence_1to3）
    ground_truth.csv  作者的标签（pmid, author_label, verdict_targets）——标注者不应查看

输出
    混淆矩阵、percent agreement、Cohen's κ + 95%CI（bootstrap）、按置信度分层的 κ

为什么用 bootstrap 算 CI 而不是解析式：五个类别里 `thr` 只有 3 例，解析式方差在稀有类别上
不稳定；bootstrap 至少让区间反映实际抽样分布。n=25 时区间会很宽——这是事实，不该被藏起来。
"""
from __future__ import annotations

import argparse
import csv
import random
import sys

LABELS = ["prov", "both", "none", "herb", "thr"]


def load(path: str) -> dict:
    rows = list(csv.DictReader(open(path, encoding="utf-8-sig")))
    if not rows:
        raise SystemExit(f"{path} 为空")
    key = "pmid"
    out = {}
    for r in rows:
        p = (r.get(key) or "").strip()
        if p:
            out[p] = r
    return out


def cohen_kappa(pairs: list[tuple[str, str]]) -> float:
    """标准 Cohen's κ。pairs = [(rater1, rater2)]。"""
    n = len(pairs)
    if n == 0:
        return float("nan")
    labels = sorted({a for a, _ in pairs} | {b for _, b in pairs})
    idx = {l: i for i, l in enumerate(labels)}
    k = len(labels)
    m = [[0] * k for _ in range(k)]
    for a, b in pairs:
        m[idx[a]][idx[b]] += 1
    po = sum(m[i][i] for i in range(k)) / n
    row = [sum(m[i]) / n for i in range(k)]
    col = [sum(m[i][j] for i in range(k)) / n for j in range(k)]
    pe = sum(row[i] * col[i] for i in range(k))
    if abs(1 - pe) < 1e-12:
        return 1.0
    return (po - pe) / (1 - pe)


def boot_ci(pairs: list[tuple[str, str]], n_iter: int = 5000, seed: int = 20261125):
    rng = random.Random(seed)
    n = len(pairs)
    vals = []
    for _ in range(n_iter):
        s = [pairs[rng.randrange(n)] for _ in range(n)]
        vals.append(cohen_kappa(s))
    vals.sort()
    return vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals))]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("sheet")
    ap.add_argument("truth")
    ap.add_argument("--min-conf", type=int, default=1,
                    help="只用置信度 >= 该值的条目（默认 1，即全部）")
    a = ap.parse_args()

    sheet, truth = load(a.sheet), load(a.truth)
    pairs, skipped_blank, skipped_conf = [], 0, 0
    for pmid, t in truth.items():
        s = sheet.get(pmid)
        if not s:
            continue
        lab = (s.get("rater_label") or "").strip().lower()
        if not lab:
            skipped_blank += 1
            continue
        try:
            conf = int((s.get("rater_confidence_1to3") or "3").strip() or 3)
        except ValueError:
            conf = 3
        if conf < a.min_conf:
            skipped_conf += 1
            continue
        pairs.append((lab, t["author_label"].strip().lower()))

    if not pairs:
        print("没有可比较的条目：请确认 rater_label 已填、且 pmid 与真值表一致")
        return 2

    n = len(pairs)
    agree = sum(1 for x, y in pairs if x == y)
    kappa = cohen_kappa(pairs)
    lo, hi = boot_ci(pairs)

    print(f"条目 n = {n}（空标签 {skipped_blank}，被置信度过滤 {skipped_conf}）")
    print(f"percent agreement = {agree}/{n} = {100*agree/n:.1f}%")
    print(f"Cohen's κ = {kappa:.3f}   95%CI [{lo:.3f}, {hi:.3f}]  （bootstrap, 5000 次）")
    print()
    print("混淆矩阵（行 = 第二位标注者，列 = 作者）：")
    cols = LABELS
    print("            " + "".join(f"{c:>8}" for c in cols))
    for r in cols:
        line = f"{r:<12}"
        for c in cols:
            line += f"{sum(1 for x, y in pairs if x == r and y == c):>8}"
        print(line)
    print()
    per = {l: (sum(1 for x, y in pairs if x == l == y),
               sum(1 for x, y in pairs if x == l))
           for l in sorted({x for x, _ in pairs})}
    print("按类别一致率：")
    for l, (ok, tot) in per.items():
        print(f"  {l:<6} {ok}/{tot}" + ("  (n<5，κ 不稳)" if tot < 5 else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
