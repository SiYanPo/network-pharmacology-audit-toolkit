#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""md_replicates_merge.py — 把 N 条独立重复 MD 轨迹汇成"可发表级"统计：重复间均值、t 区间、
套合方差分解（ICC）、异质性检验（Cochran Q / I²），并给出"再加一条重复能买到多少精度"的定量预测。

为什么需要它：
单条轨迹的逐帧 SEM 把帧当独立样本，会低估 2–3 倍（自相关未扣）。只有**独立重复**才给得出
真正的重复间 SD，才能在论文里写 "mean ± SD (n=3 independent replicates)"。
n=2 时 σ_b 只有 1 个自由度、t(0.975,1)=12.71，区间宽到没有信息量；本工具会直接量化这一点。

输入：--reps 目录下每个子目录含 rep_summary.json 与 rep_perframe.csv
      （rep_perframe.csv 列：t_ns,rmsd_backbone_A,ligand_disp_A,rg_A,used）
输出：reps_merged.csv / reps_summary.json / reps_table.md / reps_forest.png

用法：
  python tools/md_replicates_merge.py --reps 全流程实测/md_reuse/replicates \
      --outdir 全流程实测/md_reuse/replicate_stats --common-start 0.5
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

import numpy as np

METRICS = ["rmsd_backbone_A", "ligand_disp_A", "rg_A"]
LABEL = {"rmsd_backbone_A": "Backbone RMSD (Å)",
         "ligand_disp_A": "Ligand displacement (Å)",
         "rg_A": "Radius of gyration (Å)"}


# ---------- 统计小工具（不依赖 scipy） ----------
def t_crit(df: int, alpha: float = 0.05) -> float:
    """双尾 t 临界值，df<=30 用查表，更大用正态近似 + 修正。"""
    table = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
             8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160,
             14: 2.145, 15: 2.131, 16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093,
             20: 2.086, 21: 2.080, 22: 2.074, 23: 2.069, 24: 2.064, 25: 2.060,
             26: 2.056, 27: 2.052, 28: 2.048, 29: 2.045, 30: 2.042}
    if df <= 0:
        return float("inf")
    if df in table:
        return table[df]
    if df < 30:
        # 表已覆盖 1–30，此分支仅防未来改动；用相邻两点线性内插
        lo = max(k for k in table if k < df)
        hi = min(k for k in table if k > df)
        return table[lo] + (table[hi] - table[lo]) * (df - lo) / (hi - lo)
    return 1.960 + 2.4 / df                                 # df>30 粗略外推


def chi2_sf(x: float, df: int) -> float:
    """上尾概率，df 小用级数精确，大用 Wilson–Hilferty。"""
    if df <= 0 or x <= 0:
        return 1.0
    if df <= 40:
        # P(X>x) = exp(-x/2) * sum_{k<df/2} (x/2)^k / k!   （df 为偶数精确；奇数时近似）
        h = x / 2.0
        term, s = 1.0, 1.0
        for k in range(1, int(df / 2) + 1):
            term *= h / k
            s += term
        return min(1.0, math.exp(-h) * s)
    z = ((x / df) ** (1 / 3) - (1 - 2 / (9 * df))) / math.sqrt(2 / (9 * df))
    return 0.5 * math.erfc(z / math.sqrt(2))


def block_sem(x: np.ndarray, block: int) -> tuple[float, int]:
    """块平均 SEM：把序列切成 block 长的块，块均值当独立样本。"""
    nb = len(x) // block
    if nb < 2:
        return float("nan"), 0
    bm = x[:nb * block].reshape(nb, block).mean(1)
    return float(bm.std(ddof=1) / math.sqrt(nb)), nb


# ---------- 数据装载 ----------
def load_rep(d: str) -> dict:
    per = os.path.join(d, "rep_perframe.csv")
    summ = os.path.join(d, "rep_summary.json")
    if not os.path.exists(per):
        return {}
    tt, cols = [], {m: [] for m in METRICS}
    with open(per, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            if r.get("used", "1") not in ("1", "True", "true"):
                continue
            tt.append(float(r["t_ns"]))
            for m in METRICS:
                cols[m].append(float(r[m]))
    meta = {}
    if os.path.exists(summ):
        meta = json.load(open(summ, encoding="utf-8"))
    return {"name": os.path.basename(d.rstrip("/\\")),
            "t": np.array(tt), **{m: np.array(cols[m]) for m in METRICS},
            "frames": len(tt), "time_ns": float(tt[-1] - tt[0]) if tt else 0.0,
            "meta": meta}


def main() -> int:
    ap = argparse.ArgumentParser(description="多重复 MD 轨迹统计合并")
    ap.add_argument("--reps", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--common-start", type=float, default=0.5,
                    help="丢弃前若干 ns 的弛豫段（默认 0.5 ns）")
    ap.add_argument("--block-ns", type=float, default=1.0)
    ap.add_argument("--frame-ns", type=float, default=0.01)
    ap.add_argument("--project-n", default="2,3,4,5",
                    help="预测再加到 N 条重复时的 95%% 区间半宽")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    reps = []
    for name in sorted(os.listdir(a.reps)):
        d = os.path.join(a.reps, name)
        if os.path.isdir(d):
            r = load_rep(d)
            if r:
                reps.append(r)
    if not reps:
        sys.exit("未找到任何 rep_perframe.csv")

    # 公共时间窗：所有重复都覆盖的区间（保证可比）
    t_end = min(r["time_ns"] for r in reps)
    t_end = math.floor(t_end / 0.1) * 0.1
    if t_end <= a.common_start + 0.2:
        sys.exit(f"公共时间窗过短（{a.common_start}-{t_end} ns），无法比较")
    print(f"[窗口] {a.common_start}–{t_end} ns（受最短重复限制，共 {len(reps)} 条重复）")

    block = max(2, int(round(a.block_ns / a.frame_ns)))
    per_rep, summary = {}, {}
    for m in METRICS:
        means, sds, ns, sem_blocks = [], [], [], []
        for r in reps:
            mask = (r["t"] >= a.common_start) & (r["t"] <= t_end)
            x = r[m][mask]
            means.append(float(x.mean()))
            sds.append(float(x.std(ddof=1)))
            ns.append(int(x.size))
            sb, _ = block_sem(x, block)
            sem_blocks.append(sb)
            per_rep.setdefault(r["name"], {})[m] = {
                "mean": round(float(x.mean()), 3), "sd": round(float(x.std(ddof=1)), 3),
                "n_frames": int(x.size), "block_sem": None if math.isnan(sb) else round(sb, 3)}
        means_a = np.array(means)
        N = len(means_a)
        grand = float(means_a.mean())
        sd_b = float(means_a.std(ddof=1)) if N > 1 else float("nan")
        se_b = sd_b / math.sqrt(N) if N > 1 else float("nan")
        tc = t_crit(N - 1)
        ci = (grand - tc * se_b, grand + tc * se_b) if N > 1 else (float("nan"), float("nan"))

        # 套合方差：重复间 vs 重复内（重复内用块 SEM 的平方，已扣自相关）
        within = [s for s in sem_blocks if s is not None and not math.isnan(s)]
        within_ms = float(np.mean([s ** 2 for s in within])) if within else float("nan")
        between_ms = sd_b ** 2 if N > 1 else float("nan")
        icc = (between_ms / (between_ms + within_ms)
               if between_ms == between_ms and within_ms == within_ms and
               (between_ms + within_ms) > 0 else float("nan"))

        # 异质性：Cochran Q（用重复均值的块 SEM 作方差）
        Q = p = I2 = float("nan")
        if N > 2 and within:
            w = np.array([1.0 / s ** 2 for s in within])
            mu = float((w * means_a).sum() / w.sum())
            Q = float((w * (means_a - mu) ** 2).sum())
            p = chi2_sf(Q, N - 1)
            I2 = max(0.0, (Q - (N - 1)) / Q) * 100 if Q > 0 else 0.0

        summary[m] = {
            "per_replicate": {r["name"]: per_rep[r["name"]][m]["mean"] for r in reps},
            "grand_mean": round(grand, 3),
            "sd_between_replicates": None if sd_b != sd_b else round(sd_b, 3),
            "sem_between": None if se_b != se_b else round(se_b, 3),
            "ci95": None if ci[0] != ci[0] else [round(ci[0], 3), round(ci[1], 3)],
            "ci95_halfwidth": None if ci[0] != ci[0] else round((ci[1] - ci[0]) / 2, 3),
            "t_crit": None if tc == float("inf") else tc,
            "df": N - 1,
            "pooled_within_block_sem": None if within_ms != within_ms else round(math.sqrt(within_ms), 3),
            "icc_replicate": None if icc != icc else round(icc, 3),
            "cochran_Q": None if Q != Q else round(Q, 3),
            "Q_p": None if p != p else round(p, 4),
            "I2_percent": None if I2 != I2 else round(I2, 1),
            "n_replicates": N,
        }

    # ---- 加一条重复能买到什么：给定 σ_b，CI 半宽 ∝ t(N-1)/√N ----
    proj = {}
    for m in METRICS:
        sb = summary[m]["sd_between_replicates"]
        if not sb:
            continue
        ref = summary[m]["ci95_halfwidth"]
        row = {}
        for n in [int(x) for x in a.project_n.split(",") if x.strip()]:
            row[n] = round(t_crit(n - 1) * sb / math.sqrt(n), 3)
        proj[m] = {"sd_b_used": sb, "current_halfwidth": ref,
                   "halfwidth_by_n": row,
                   "shrink_vs_current": ({n: (round(ref / v, 2) if ref and v else None)
                                          for n, v in row.items()} if ref else {})}

    # ---- 落盘 ----
    rows = []
    for r in reps:
        for m in METRICS:
            d = per_rep[r["name"]][m]
            rows.append({"replicate": r["name"], "metric": m, "mean": d["mean"],
                         "sd_frames": d["sd"], "block_sem": d["block_sem"],
                         "n_frames_common": d["n_frames"], "window_ns": t_end - a.common_start})
    with open(os.path.join(a.outdir, "reps_merged.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    out = {"window_ns": [a.common_start, t_end], "n_replicates": len(reps),
           "replicates": [r["name"] for r in reps],
           "drop_first_ns": a.common_start, "block_ns": a.block_ns,
           "metrics": summary, "projection_add_replicate": proj,
           "caveat": ("n=2 时重复间 SD 仅 1 个自由度；t(0.975,1)=12.71，区间没有信息量，"
                      "不能作为可发表的误差估计。ICC 与 Cochran Q 需 n>=3 才有意义。")}
    json.dump(out, open(os.path.join(a.outdir, "reps_summary.json"), "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)

    # ---- 森林图 ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, len(METRICS), figsize=(12, 3.2), dpi=160)
        for ax, m in zip(np.atleast_1d(axes), METRICS):
            vals = [per_rep[r["name"]][m]["mean"] for r in reps]
            errs = [per_rep[r["name"]][m]["block_sem"] or 0 for r in reps]
            y = np.arange(len(reps))[::-1]
            ax.errorbar(vals, y, xerr=errs, fmt="o", color="#c0392b", capsize=3, ms=6)
            for yi, v in zip(y, vals):
                ax.annotate(f"{v:.2f}", (v, yi), textcoords="offset points",
                            xytext=(0, 8), ha="center", fontsize=8)
            ax.axvline(summary[m]["grand_mean"], color="#2c3e50", ls="--", lw=1,
                       label=f"grand mean {summary[m]['grand_mean']:.2f}")
            ax.set_yticks(y); ax.set_yticklabels([r["name"].replace("replicate_", "rep") for r in reps])
            ax.set_title(f"{LABEL[m]}\n(±block SEM)", fontsize=9)
            ax.grid(alpha=.25, lw=.5, axis="x"); ax.legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, "reps_forest.png"))
    except Exception as exc:      # 图失败不影响统计
        print(f"[warn] 森林图未生成: {exc}")

    # ---- 控制台 ----
    for m in METRICS:
        s = summary[m]
        print(f"\n{LABEL[m]}")
        print(f"  per-rep: " + ", ".join(f"{k.replace('replicate_','rep')}={v:.2f}"
                                         for k, v in s["per_replicate"].items()))
        print(f"  grand mean = {s['grand_mean']:.2f} ± {s['sd_between_replicates']} (SD, n={s['n_replicates']})"
              f"   95%CI {s['ci95']}  halfwidth={s['ci95_halfwidth']}")
        print(f"  pooled within block SEM = {s['pooled_within_block_sem']}   ICC(rep) = {s['icc_replicate']}")
        if s["cochran_Q"] is not None:
            print(f"  heterogeneity: Q={s['cochran_Q']} df={s['df']} p={s['Q_p']} I²={s['I2_percent']}%")
        if m in proj:
            print(f"  再加重复的 95%CI 半宽: " +
                  ", ".join(f"n={n}:{v}" for n, v in proj[m]["halfwidth_by_n"].items()))
    print(f"\n→ {os.path.join(a.outdir, 'reps_summary.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
