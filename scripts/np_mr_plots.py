#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_mr_plots.py —— MR / COLOC 的出版级标准图（森林图 / 散点图 / LocusZoom 风格区域图）

为什么必须有图
    MR 的数值表格审稿人不会细看，但森林图看一眼就知道"方向一致不一致、有没有离群点"；
    LocusZoom 区域图能直接看出来"表达信号与疾病信号是不是同一个峰"。
    只有数值没有图的因果分析，在方法学审稿人眼里是不完整的。

三张图（都用 matplotlib，300 dpi PNG + 矢量 SVG；图内英文，避免字体缺失）
    1 forest  : 逐 SNP Wald 比（OR 与 95% CI）+ IVW / MR-Egger / 加权中位数汇总
    2 scatter : β_y vs β_x 散点 + 各方法拟合线（斜率即因果估计）
    3 locus   : 按染色体位置画 -log10(p) 双向轨道（上=暴露，下=结局），同一区域对齐

输入
    `np_mr.py analyze` 产出的 mr_result.json（含 snps 数组；旧结果无 snps 时只画汇总）
    COLOC 的 LocusZoom 需要 --exposure/--outcome 两份原始 TSV（含 pos 列）

用法
    python np_mr_plots.py --mr np_out/mr/mr_result.json -o np_out/mr/figs
    python np_mr_plots.py locus --exposure exp.tsv --outcome out.tsv -o figs/
    python np_mr_plots.py selftest
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    })
    return plt


def _save(fig, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, stem + ".png")
    svg = os.path.join(outdir, stem + ".svg")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    return png, svg


def _wald(s):
    """单 SNP Wald 比与 delta 法 CI（OR 尺度）。"""
    bx, by = s.get("beta_x"), s.get("beta_y")
    sx, sy = s.get("se_x") or 0, s.get("se_y") or 0
    if not bx or abs(bx) < 1e-12 or by is None:
        return None
    ratio = by / bx
    se = math.sqrt((sy ** 2) / (bx ** 2) + (by ** 2) * (sx ** 2) / (bx ** 4)) if bx else None
    if not se:
        return None
    return {"snp": s.get("snp", ""), "or": math.exp(ratio),
            "lo": math.exp(ratio - 1.96 * se), "hi": math.exp(ratio + 1.96 * se)}


def forest(mr, outdir, title="Mendelian randomization"):
    plt = _mpl()
    snps = mr.get("snps") or []
    rows = [w for w in (_wald(s) for s in snps) if w]
    pooled = []
    for key, lab in (("IVW", "IVW"), ("MR_Egger", "MR-Egger"), ("WeightedMedian", "Weighted median")):
        m = mr.get(key) or {}
        b, se = m.get("beta"), m.get("se") or m.get("se_approx")
        if b is None:
            continue
        if se:
            lo, hi = math.exp(b - 1.96 * se), math.exp(b + 1.96 * se)
        else:
            lo = hi = math.exp(b)
        pooled.append({"snp": lab, "or": math.exp(b), "lo": lo, "hi": hi})
    if not pooled:
        return None, "无汇总估计"
    labels = [r["snp"] for r in rows] + [r["snp"] for r in pooled]
    ors = [r["or"] for r in rows] + [r["or"] for r in pooled]
    los = [r["lo"] for r in rows] + [r["lo"] for r in pooled]
    his = [r["hi"] for r in rows] + [r["hi"] for r in pooled]
    n = len(labels)
    fig, ax = plt.subplots(figsize=(5.4, max(2.4, 0.25 * n + 1.2)), dpi=300)
    y = list(range(n))[::-1]
    for i, (yy, o, lo, hi) in enumerate(zip(y, ors, los, his)):
        is_pool = i >= len(rows)
        ax.plot([lo, hi], [yy, yy], color="#007B7A" if is_pool else "#6b7280",
                lw=1.8 if is_pool else 1.0)
        ax.plot(o, yy, marker="D" if is_pool else "o", ms=5 if is_pool else 3.5,
                color="#007B7A" if is_pool else "#374151", zorder=3)
    ax.axvline(1.0, ls="--", lw=0.8, color="#b91c1c")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.5 if n > 12 else 7.5)
    ax.set_xscale("log")
    # 对数轴的主/次刻度：默认会生成 1.3×10^0 这类科学计数并互相压字，这里改为固定小数刻度
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    lo_v = min(los + [1.0])
    hi_v = max(his + [1.0])
    span = hi_v - lo_v
    if span <= 0:
        span = 0.2
    step = 10 ** math.floor(math.log10(span / 4.0)) if span > 0 else 0.1
    for mult in (1, 2, 2.5, 5, 10):
        if span / (step * mult) <= 7:
            step *= mult
            break
    ticks = [lo_v + i * step for i in range(0, int(span / step) + 2)]
    ticks = [t for t in ticks if t > 0 and t <= hi_v + step]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: ("%.2f" % v).rstrip("0").rstrip(".")))
    ax.set_xlabel("Odds ratio per genetically predicted expression (95% CI)", fontsize=8)
    trans = ax.get_yaxis_transform()          # x 用轴比例、y 用数据坐标 → 不会互相压字
    ax.text(1.03, -0.75, "OR [95% CI]", transform=trans, fontsize=6.5,
            color="#6b7280", va="center")
    for i, r in enumerate(pooled):
        ax.text(1.03, y[len(rows) + i], "%.2f [%.2f, %.2f]" % (r["or"], r["lo"], r["hi"]),
                transform=trans, fontsize=6.4, color="#007B7A", va="center")
    ax.set_ylim(-1.25, n - 0.2)
    ax.set_title(title, fontsize=9)
    if rows:
        ax.axhline(y[len(rows) - 1] - 0.5, color="#e5e7eb", lw=0.8)
    fig.tight_layout()
    png, svg = _save(fig, outdir, "mr_forest")
    plt.close(fig)
    return png, ""


def scatter(mr, outdir, title="SNP effect: exposure vs outcome"):
    plt = _mpl()
    snps = mr.get("snps") or []
    if len(snps) < 3:
        return None, "SNP<3，无法画散点"
    x = [s.get("beta_x") or 0 for s in snps]
    y = [s.get("beta_y") or 0 for s in snps]
    fig, ax = plt.subplots(figsize=(3.8, 3.6), dpi=300)
    ax.scatter(x, y, s=14, color="#374151", alpha=0.75, edgecolor="white", linewidth=0.4, zorder=3)
    xr = [min(x), max(x)]
    styles = [("IVW", "#007B7A", "-"), ("MR_Egger", "#c0504d", "--"), ("WeightedMedian", "#4f8a5b", ":")]
    for key, color, ls in styles:
        b = (mr.get(key) or {}).get("beta")
        if b is None:
            continue
        ax.plot(xr, [b * v for v in xr], ls=ls, lw=1.1, color=color, label="%s (β=%.3f)" % (key, b))
    ax.axhline(0, color="#e5e7eb", lw=0.8)
    ax.axvline(0, color="#e5e7eb", lw=0.8)
    ax.set_xlabel("SNP effect on exposure (β)", fontsize=8)
    ax.set_ylabel("SNP effect on outcome (β)", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=6.3, frameon=False, loc="best")
    fig.tight_layout()
    png, svg = _save(fig, outdir, "mr_scatter")
    plt.close(fig)
    return png, ""


def _read_tsv(path):
    import csv
    out = []
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        for r in csv.DictReader(f, delimiter=delim):
            try:
                out.append({"snp": r.get("SNP") or r.get("rsid") or "",
                            "pos": float(r.get("pos") or r.get("position") or 0),
                            "p": float(r.get("p") or 1),
                            "beta": float(r.get("beta") or 0)})
            except (TypeError, ValueError):
                continue
    return out


def locus(exposure_tsv, outcome_tsv, outdir, title="Regional association (LocusZoom style)"):
    """上下双轨：上=暴露(eQTL)，下=结局(GWAS)，共享 x 轴位置。"""
    plt = _mpl()
    ex = [r for r in _read_tsv(exposure_tsv) if r["pos"] > 0]
    ou = [r for r in _read_tsv(outcome_tsv) if r["pos"] > 0]
    if not ex or not ou:
        return None, "缺 pos 列或数据为空"
    lo = min(min(r["pos"] for r in ex), min(r["pos"] for r in ou))
    hi = max(max(r["pos"] for r in ex), max(r["pos"] for r in ou))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.0, 3.8), dpi=300, sharex=True)
    for ax, data, color, lab in ((ax1, ex, "#007B7A", "Exposure (eQTL)"),
                                 (ax2, ou, "#c0504d", "Outcome (GWAS)")):
        xs = [r["pos"] / 1e6 for r in data]
        ys = [-math.log10(max(r["p"], 1e-300)) for r in data]
        ax.scatter(xs, ys, s=8, color=color, alpha=0.65, edgecolor="none")
        ax.set_ylabel("-log10 P", fontsize=7.5)
        ax.text(0.01, 0.92, lab, transform=ax.transAxes, fontsize=7.5, color=color)
        ax.axhline(-math.log10(5e-8), ls="--", lw=0.7, color="#9aa3ad")
    ax2.set_xlabel("Chromosomal position (Mb)", fontsize=8)
    ax1.set_title(title, fontsize=9)
    ax1.set_xlim(lo / 1e6, hi / 1e6)
    fig.tight_layout()
    png, svg = _save(fig, outdir, "coloc_locuszoom")
    plt.close(fig)
    return png, ""


def selftest():
    import tempfile
    plt = _mpl()
    d = tempfile.mkdtemp(prefix="np_mrplots_")
    snps = [{"snp": "rs%d" % i, "beta_x": 0.1 + 0.02 * i, "se_x": 0.02,
             "beta_y": 0.3 * (0.1 + 0.02 * i) + (0.002 if i % 2 else -0.002),
             "se_y": 0.01, "pos": 1000000 + i * 50000} for i in range(8)]
    mr = {"status": "OK", "snps": snps,
          "IVW": {"beta": 0.30, "se": 0.02}, "MR_Egger": {"beta": 0.31, "se": 0.05},
          "WeightedMedian": {"beta": 0.295, "se_approx": 0.03}}
    p1, e1 = forest(mr, d)
    p2, e2 = scatter(mr, d)
    e = os.path.join(d, "e.tsv")
    o = os.path.join(d, "o.tsv")
    open(e, "w", encoding="utf-8").write("SNP\tpos\tp\tbeta\n" + "\n".join(
        "rs%d\t%d\t%.2g\t0.2" % (i, 1000000 + i * 50000, 10 ** (-6 - i)) for i in range(8)))
    open(o, "w", encoding="utf-8").write("SNP\tpos\tp\tbeta\n" + "\n".join(
        "rs%d\t%d\t%.2g\t0.05" % (i, 1000000 + i * 50000, 10 ** (-3 - i)) for i in range(8)))
    p3, e3 = locus(e, o, d)
    ok = all(p and os.path.getsize(p) > 3000 for p in (p1, p2, p3))
    print("selftest: %s（forest %s / scatter %s / locus %s）" % (
        "PASS" if ok else "FAIL",
        "%.0fKB" % (os.path.getsize(p1) / 1024) if p1 else e1,
        "%.0fKB" % (os.path.getsize(p2) / 1024) if p2 else e2,
        "%.0fKB" % (os.path.getsize(p3) / 1024) if p3 else e3))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = argparse.ArgumentParser(description="np_mr_plots: MR/COLOC 出版级图")
    sub = ap.add_subparsers(dest="cmd")
    p = sub.add_parser("all")
    p.add_argument("--mr", required=True)
    p.add_argument("-o", "--outdir", default="np_out/mr/figs")
    p.add_argument("--exposure", default=None)
    p.add_argument("--outcome", default=None)
    l = sub.add_parser("locus")
    l.add_argument("--exposure", required=True)
    l.add_argument("--outcome", required=True)
    l.add_argument("-o", "--outdir", default="np_out/mr/figs")
    a = ap.parse_args(argv)

    def _run_all(mr_path, outdir, exp=None, out=None):
        mr = json.load(open(mr_path, encoding="utf-8"))
        if mr.get("status") != "OK":
            print("[mrplots] MR 状态=%s，跳过" % mr.get("status"))
            return 2
        made = []
        for fn, args in ((forest, (mr, outdir)), (scatter, (mr, outdir))):
            p, e = fn(*args)
            if p:
                made.append(p)
            else:
                print("[mrplots] %s 跳过：%s" % (fn.__name__, e))
        if exp and out:
            p, e = locus(exp, out, outdir)
            if p:
                made.append(p)
            else:
                print("[mrplots] locus 跳过：%s" % e)
        for p in made:
            print("[mrplots] → %s" % p)
        return 0 if made else 2

    if a.cmd == "locus":
        p, e = locus(a.exposure, a.outcome, a.outdir)
        print("[mrplots] → %s" % p if p else "[mrplots] 失败：%s" % e)
        return 0 if p else 2
    if a.cmd == "all":
        return _run_all(a.mr, a.outdir, a.exposure, a.outcome)
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
