#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_benchmark.py —— 一键跑完 4 个经典方剂的金标准基准，并出论文级柱状图/ROC

为什么需要它
    上一轮把基准引擎写好了，但报告里只留了一句"请自行跑 evaluate"。
    审稿人不接受"你自己去跑"，合作者也不会。这个脚本把"准备输入 → 逐方剂评测 → 汇总表 → 出图"
    串成一条命令，缺数据时**明确列出缺什么**，绝不编造数字。

输入怎么给（三种，按优先级）
    1. --map CDDP=cddp_edges.csv,HLJDT=hljdt_edges.csv ...   直接指定每个方剂的成分-靶点边表
    2. --edges-dir DIR                                       目录下按 <ID>*.csv/tsv 自动找
    3. --pipeline-dir DIR                                    目录下有 inputs_<ID>.json 时，先跑流水线再评测
    都没有的方剂 → 记为 SKIP 并写明原因（缺该方剂的成分-靶点表），summary 里不出现伪造数值。

输出（默认 benchmarks/）
    benchmark_summary.csv          每个方剂一行：AUC / AUPRC / Recall@K / Precision@K / EF / coverage
    benchmark_summary.json         同上（含 skip 原因）
    fig_benchmark_bars.png/svg     分组柱状图：Precision@20、Recall@20、Recall@50 + 随机基线
    fig_benchmark_roc.png/svg      多曲线 ROC（每个方剂一条）
    benchmark_report.md            由 np_benchmark 生成的完整报告（含 Methods 段）

用法
    python run_benchmark.py --demo                       # 构造数据自检（AUC≈0.5，仅验证管线）
    python run_benchmark.py --edges-dir data/formula_edges/
    python run_benchmark.py --map CDDP=a.csv,HLJDT=b.csv --outdir benchmarks
    python run_benchmark.py selftest
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "scripts")
sys.path.insert(0, SCRIPTS)

import np_benchmark as NB  # noqa: E402

FORMULAS = ["CDDP", "HLJDT", "LHQW", "XFZYT"]
LABEL_EN = {"CDDP": "CDDP (Danshen)", "HLJDT": "HLJDT", "LHQW": "LHQW", "XFZYT": "XFZYT"}
COLORS = ["#007B7A", "#c0504d", "#4f8a5b", "#7a6fb0"]


def resolve_inputs(args):
    """返回 {formula_id: (edges_path, note)}，缺的也保留并写 note。"""
    out = {}
    mapping = {}
    if args.map:
        for kv in args.map.split(","):
            if "=" in kv:
                k, v = kv.split("=", 1)
                mapping[k.strip().upper()] = v.strip()
    for fid in FORMULAS:
        if fid in mapping:
            p = mapping[fid]
            out[fid] = (p if os.path.exists(p) else None,
                        "" if os.path.exists(p) else "指定的边表不存在：%s" % p)
            continue
        if args.edges_dir:
            hits = []
            for ext in ("csv", "tsv", "txt"):
                hits += glob.glob(os.path.join(args.edges_dir, "%s*%s" % (fid, ext)))
                hits += glob.glob(os.path.join(args.edges_dir, "*%s*%s" % (fid.lower(), ext)))
            if hits:
                out[fid] = (sorted(hits)[0], "")
                continue
        out[fid] = (None, "缺该方剂的成分-靶点边表（用 --map %s=文件 指定）" % fid)
    return out


def run_pipeline_for(fid, pipeline_dir, outroot):
    """若提供了 inputs_<ID>.json，就先跑流水线，再返回其标准化边表。"""
    inp = os.path.join(pipeline_dir, "inputs_%s.json" % fid)
    if not os.path.exists(inp):
        return None, "无 inputs_%s.json" % fid
    od = os.path.join(outroot, "np_out_%s" % fid)
    cmd = [sys.executable, os.path.join(SCRIPTS, "np_pipeline.py"), "--inputs", inp, "-o", od]
    print("[benchmark] 跑流水线 %s：%s" % (fid, " ".join(os.path.basename(c) for c in cmd[:3])))
    rc = subprocess.call(cmd)
    if rc != 0:
        return None, "流水线失败（rc=%d），见 %s" % (rc, od)
    std = os.path.join(od, "02_standardize", "targets_standardized.csv")
    return (std if os.path.exists(std) else None), ("" if os.path.exists(std) else "标准化产物缺失")


def evaluate_one(fid, edges, gold, tier="lenient", topk=(10, 20, 50), boot=300):
    ranked = NB.rank_from_edges(edges, weight_mode="penalized")
    formula = next(f for f in gold["formulas"] if f["id"] == fid)
    pos, amb, ginfo = NB.gold_sets(formula, tier)
    r = NB.evaluate(ranked, pos, amb, topk=topk, n_boot=boot)
    r["label"] = formula["name_zh"]
    r["formula_id"] = fid
    r["tier_label"] = ginfo["tier"]
    r["gold_meta"] = dict(ginfo)
    r["edges_file"] = os.path.abspath(edges)
    r["n_edges"] = sum(1 for _ in open(edges, encoding="utf-8-sig", errors="replace")) - 1
    return r, formula


def plot_bars(summary, outdir):
    """分组柱状图：Precision@20 / Recall@20 / Recall@50，附随机基线。图内文字用英文（期刊安全）。"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None, "matplotlib 不可用"
    rows = [r for r in summary if r.get("status") == "OK"]
    if not rows:
        return None, "无可用结果"
    import numpy as np
    metrics = [("precision", "20", "Precision@20"), ("recall", "20", "Recall@20"),
               ("recall", "50", "Recall@50")]
    x = np.arange(len(rows))
    w = 0.26
    fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=300)
    for i, (m, k, lab) in enumerate(metrics):
        vals = [r["at_k"][k][m] for r in rows]
        bars = ax.bar(x + (i - 1) * w, vals, w, label=lab,
                      color=["#007B7A", "#4f8a5b", "#7a6fb0"][i], edgecolor="white", linewidth=0.6)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.008, "%.2f" % v,
                    ha="center", va="bottom", fontsize=6.5)
    base = [r["at_k"]["20"]["expected_hits"] / 20.0 for r in rows]
    ax.plot(x, base, "k--", lw=1, marker="o", ms=3, label="Random baseline (P@20)")
    ax.set_xticks(x)
    ax.set_xticklabels([LABEL_EN.get(r["formula_id"], r["formula_id"]) for r in rows], fontsize=8)
    ax.set_ylabel("Value", fontsize=8.5)
    ax.set_title("Gold-standard benchmark (core+extended targets)", fontsize=9.5)
    ax.legend(fontsize=7, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.14))
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7.5)
    fig.tight_layout()
    png = os.path.join(outdir, "fig_benchmark_bars.png")
    svg = os.path.join(outdir, "fig_benchmark_bars.svg")
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    return png, ""


def plot_roc(summary, outdir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None, "matplotlib 不可用"
    rows = [r for r in summary if r.get("status") == "OK" and r.get("roc_points")]
    if not rows:
        return None, "无 ROC 数据"
    fig, ax = plt.subplots(figsize=(3.9, 3.9), dpi=300)
    for i, r in enumerate(rows):
        pts = r["roc_points"]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], lw=1.4,
                color=COLORS[i % len(COLORS)],
                label="%s (AUC %.2f)" % (LABEL_EN.get(r["formula_id"], r["formula_id"]), r["AUC"]))
    ax.plot([0, 1], [0, 1], ls="--", lw=0.8, color="#9aa3ad")
    ax.set_xlabel("1 - Specificity", fontsize=8)
    ax.set_ylabel("Sensitivity", fontsize=8)
    ax.set_title("ROC — formula target recovery", fontsize=9)
    ax.legend(fontsize=6.5, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    fig.tight_layout()
    png = os.path.join(outdir, "fig_benchmark_roc.png")
    svg = os.path.join(outdir, "fig_benchmark_roc.svg")
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    return png, ""


def write_summary(summary, outdir, skips):
    cols = ["formula_id", "label", "tier_label", "status", "n_edges", "n_evaluated", "n_pos_total",
            "n_positive_found", "AUC", "AUC_95CI", "AUPRC", "coverage",
            "recall@10", "recall@20", "recall@50", "precision@20", "EF@20", "edges_file"]
    with open(os.path.join(outdir, "benchmark_summary.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in summary:
            row = dict(r)
            row["AUC_95CI"] = "–".join(str(x) for x in (r.get("AUC_95CI") or [])) or "NA"
            for k in ("10", "20", "50"):
                if r.get("status") == "OK":
                    row["recall@" + k] = r["at_k"][k]["recall"]
                    if k == "20":
                        row["precision@20"] = r["at_k"][k]["precision"]
                        row["EF@20"] = r["at_k"][k]["enrichment_factor"]
            w.writerow(row)
    json.dump({"results": [{k: v for k, v in r.items() if k != "roc_points"} for r in summary],
               "skipped": skips,
               "note": "AUC/Recall 只对 status=OK 的方剂有意义；skipped 表示缺该方剂的成分-靶点边表，未做任何估算"},
              open(os.path.join(outdir, "benchmark_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)


def selftest():
    """用临时构造边表跑通评测+出图（验证管线，不代表任何生物学结论）。"""
    import tempfile
    d = tempfile.mkdtemp(prefix="np_runbench_")
    gold = NB.load_gold()
    f = next(x for x in gold["formulas"] if x["id"] == "CDDP")
    core, fakes = list(f["core"]), ["FAKEGENE1", "FAKEGENE2", "FAKEGENE3"]
    rows = ["compound,gene"]
    for g in core:                      # 两个成分覆盖全部金标准靶点
        rows += ["c1,%s" % g, "c2,%s" % g]
    for g in fakes:                     # 第三个成分只覆盖假靶点 → 金标准应排在前面
        rows.append("c3,%s" % g)
    p = os.path.join(d, "cddp.csv")
    open(p, "w", encoding="utf-8").write("\n".join(rows) + "\n")
    r, _ = evaluate_one("CDDP", p, gold, tier="lenient", topk=(10, 20, 50), boot=30)
    png, err = plot_bars([r], d)
    ok = (r["status"] == "OK" and r["AUC"] > 0.8 and png and os.path.getsize(png) > 5000)
    print("selftest: %s（AUC=%.3f，Recall@20=%.2f，图=%.0f KB）" % (
        "PASS" if ok else "FAIL", r.get("AUC") or 0,
        r["at_k"]["20"]["recall"] if r.get("status") == "OK" else 0,
        (os.path.getsize(png) / 1024) if png else 0))
    if not ok:
        print("       详情：status=%s err=%s" % (r.get("status"), err))
    shutil.rmtree(d, ignore_errors=True)
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = argparse.ArgumentParser(description="run_benchmark: 一键跑 4 方剂金标准基准并出图")
    ap.add_argument("--map", default=None, help="CDDP=a.csv,HLJDT=b.csv")
    ap.add_argument("--edges-dir", default=None)
    ap.add_argument("--pipeline-dir", default=None, help="含 inputs_<ID>.json 的目录")
    ap.add_argument("--demo", action="store_true", help="用 examples/ 构造数据（仅验证管线）")
    ap.add_argument("--tier", default="lenient", choices=["strict", "lenient", "strict-l3"])
    ap.add_argument("--topk", default="10,20,50")
    ap.add_argument("--bootstrap", type=int, default=300)
    ap.add_argument("--outdir", default=os.path.join(HERE, "benchmarks"))
    ap.add_argument("--no-fig", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    gold = NB.load_gold()
    topk = tuple(int(x) for x in a.topk.split(",") if x.strip())

    if a.demo:
        # 示例边表是"蛋白全名"，必须先标准化成基因符号，否则与金标准无交集
        src0 = os.path.join(HERE, "examples", "compound_targets.csv")
        std_dir = os.path.join(a.outdir, "_demo_standardized")
        os.makedirs(std_dir, exist_ok=True)
        cmd = [sys.executable, os.path.join(SCRIPTS, "np_standardize.py"), "-i", src0, "-o", std_dir, "--quiet"]
        rc = subprocess.call(cmd)
        src = os.path.join(std_dir, "targets_standardized.csv")
        if rc != 0 or not os.path.exists(src):
            print("[benchmark] DEMO 标准化失败（rc=%d）；请检查 examples/ 与 np_standardize" % rc)
            return 2
        inputs = {fid: (src, "构造数据（synthetic），仅验证管线，不得作生物学解读") for fid in FORMULAS}
        print("[benchmark] DEMO 模式：examples/ 构造数据 → 标准化后评测，AUC≈0.5 属预期")
    else:
        inputs = resolve_inputs(a)

    summary, skips = [], []
    for fid in FORMULAS:
        edge, note = inputs.get(fid, (None, "未提供"))
        if not edge and a.pipeline_dir:
            edge, note = run_pipeline_for(fid, a.pipeline_dir, a.outdir)
        if not edge:
            skips.append({"formula_id": fid, "reason": note})
            print("[benchmark] %-6s SKIP：%s" % (fid, note))
            continue
        try:
            r, formula = evaluate_one(fid, edge, gold, a.tier, topk, a.bootstrap)
        except SystemExit as e:
            skips.append({"formula_id": fid, "reason": str(e)})
            print("[benchmark] %-6s SKIP：%s" % (fid, e))
            continue
        if r.get("status") != "OK":
            skips.append({"formula_id": fid, "reason": r.get("reason", r.get("status"))})
            print("[benchmark] %-6s 无交集：%s" % (fid, r.get("reason")))
            continue
        # ROC 取点（供多曲线图）
        genes = [g for g, _ in NB.rank_from_edges(edge, weight_mode="penalized")]
        pos, amb, _ = NB.gold_sets(formula, a.tier)
        uni = [g for g in genes if g not in amb]
        sc = dict(NB.rank_from_edges(edge, weight_mode="penalized"))
        r["roc_points"] = NB.roc_points([1 if g in pos else 0 for g in uni], [sc[g] for g in uni])
        summary.append(r)
        print("[benchmark] %-6s AUC=%.3f  Recall@20=%.2f  EF@20=%.2fx（%d 成分边）" % (
            fid, r["AUC"], r["at_k"]["20"]["recall"], r["at_k"]["20"]["enrichment_factor"], r["n_edges"]))

    if not summary:
        print("\n[benchmark] 没有可评测的方剂。请用 --map 指定边表，例如：")
        print("  python run_benchmark.py --map CDDP=path/to/cddp_edges.csv,HLJDT=path/to/hljdt_edges.csv")
        write_summary([], a.outdir, skips)
        return 2

    write_summary(summary, a.outdir, skips)
    figs = []
    if not a.no_fig:
        p1, e1 = plot_bars(summary, a.outdir)
        p2, e2 = plot_roc(summary, a.outdir)
        figs = [p for p in (p1, p2) if p]
        if e1:
            print("[benchmark] 柱状图未生成：%s" % e1)
    # 复用 np_benchmark 的报告生成（含 Methods 段与 Limitations）
    for r in summary:
        r["top_missed"] = []
    rep = NB.write_report(a.outdir, summary,
                          note="由 run_benchmark.py 自动生成；tier=%s，bootstrap=%d。"
                               "未评测的方剂见 benchmark_summary.json 的 skipped 段。" % (a.tier, a.bootstrap))
    print("\n[benchmark] 汇总 → %s" % os.path.join(a.outdir, "benchmark_summary.csv"))
    print("[benchmark] 报告 → %s" % rep)
    for p in figs:
        print("[benchmark] 图 → %s" % p)
    if skips:
        print("[benchmark] 跳过 %d 个方剂：%s" % (len(skips), ", ".join(s["formula_id"] for s in skips)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
