#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_fdr.py — 全局多重性控制与"多重性声明"（v3.0）

为什么需要（评审报告 §4，七寸之四）
------------------------------------
此前只有 `np_enrich` 做了 BH；但全流程有 20+ 模块、几十个靶点、多个终点
（MR p、富集 p、基准超几何 p）。每个模块各自用 0.05，等于在 20 个地方各开一扇门。
target-wide MR 里"p<0.05 就叫因果证据"是硬伤——检测 30 个靶点时，期望假阳性 1.5 个。

本模块做三件事
--------------
1. **收口**：把所有（终点 × 假设）的 p 值汇总成一张表，附来源文件哈希。
2. **校正**：族内 BH（同一终点内部）+ 全局 BH（跨全部终点），输出 q 值。
3. **声明**：生成可直接写进 Methods 的"多重性声明"段落，含检测数、方法、校正前后显著数。

为什么同时报族内与全局
----------------------
族内 BH 回答"在 MR 这一族里，哪些靶点还站得住"；全局 BH 回答"整个研究结论里，
有多少是被多重比较撑起来的"。只报前者会被质疑"你在 20 个族里各挑一个显著的"；
只报后者会过度保守（不同终点的先验不同）。两者并列，读者可自行判断。

用法
----
    python np_fdr.py run --dir np_out -o np_out/09_evidence
    python np_fdr.py mr  --dir np_out          # 只看 MR 一族
    python np_fdr.py statement --dir np_out    # 只打印 Methods 段落
    python np_fdr.py selftest
"""

import argparse
import csv
import glob
import hashlib
import json
import os
import sys
import time

VERSION = "3.0"


# ============================ 统计 ============================
def bh(pvals):
    """Benjamini–Hochberg。返回 q 值（与输入同序），并保证单调。"""
    n = len(pvals)
    if n == 0:
        return []
    idx = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = idx[rank - 1]
        val = pvals[i] * n / rank
        prev = min(prev, val)
        q[i] = min(1.0, prev)
    return q


def bonferroni(pvals):
    n = max(1, len(pvals))
    return [min(1.0, p * n) for p in pvals]


# ============================ 采集 ============================
def _sha16(path):
    try:
        return hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]
    except Exception:
        return ""


def _read_csv(path):
    if not path or not os.path.exists(path):
        return [], []
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        rows = list(csv.DictReader(f, delimiter=delim))
        hdrs = list(rows[0].keys()) if rows else []
    return rows, hdrs


def _f(x, default=None):
    try:
        v = float(x)
        return v if v == v else default
    except (TypeError, ValueError):
        return default


def gather(outdir):
    """汇总所有可校正的 p 值。返回 rows（每行一个假设）。"""
    rows = []

    # 1) MR：每个靶点一个 IVW p
    for p in glob.glob(os.path.join(outdir, "**", "mr_result.json"), recursive=True):
        try:
            j = json.load(open(p, encoding="utf-8"))
        except Exception:
            continue
        if j.get("status") != "OK":
            continue
        ivw = j.get("IVW") or {}
        pv = _f(ivw.get("p"))
        if pv is None:
            continue
        gene = j.get("gene") or os.path.basename(os.path.dirname(p))
        rows.append({"family": "MR_IVW", "unit": gene, "test": "causal_effect",
                     "stat": _f(ivw.get("beta")), "p": pv,
                     "source": os.path.relpath(p, outdir), "sha256_16": _sha16(p)})
        col = j.get("coloc") or {}
        if col.get("status") == "OK" and _f(col.get("PP.H4")) is not None:
            # 共定位用 PP.H4 而非 p；不进 BH，只做记录（避免口径混淆）
            rows.append({"family": "COLOC_PP4", "unit": gene, "test": "pp4",
                         "stat": _f(col.get("PP.H4")), "p": None,
                         "source": os.path.relpath(p, outdir), "sha256_16": _sha16(p)})

    # 2) 富集：每个 term 一个超几何 p
    for p in glob.glob(os.path.join(outdir, "**", "enrich_all.csv"), recursive=True):
        r, h = _read_csv(p)
        if not r:
            continue
        low = {x.lower(): x for x in h}
        tc = low.get("term") or low.get("term_name")
        pc = low.get("pvalue") or low.get("p_value") or low.get("p")
        if not (tc and pc):
            continue
        for x in r:
            pv = _f(x.get(pc))
            if pv is None:
                continue
            rows.append({"family": "ENRICH", "unit": x.get(tc, ""), "test": "hypergeom",
                         "stat": _f(x.get("n_hit")), "p": pv,
                         "source": os.path.relpath(p, outdir), "sha256_16": _sha16(p)})

    # 3) 基准：超几何 p（每个 K 一个检验）
    for p in glob.glob(os.path.join(outdir, "**", "benchmark_*.json"), recursive=True):
        try:
            j = json.load(open(p, encoding="utf-8"))
        except Exception:
            continue
        if j.get("status") != "OK" or j.get("kind") == "concordance":
            continue
        for K, v in (j.get("at_k") or {}).items():
            pv = _f(v.get("hypergeom_p"))
            if pv is None:
                continue
            rows.append({"family": "BENCHMARK_AT_K", "unit": "%s@%s" % (j.get("label", "?"), K),
                         "test": "hypergeom_vs_random", "stat": v.get("hits"), "p": pv,
                         "source": os.path.relpath(p, outdir), "sha256_16": _sha16(p)})

    # 4) 亚网/零模型（经验 p）
    for p in glob.glob(os.path.join(outdir, "**", "nullmodel*.json"), recursive=True) + \
            glob.glob(os.path.join(outdir, "**", "null_*.json"), recursive=True):
        try:
            j = json.load(open(p, encoding="utf-8"))
        except Exception:
            continue
        pv = _f(j.get("p_value") or j.get("p"))
        if pv is None:
            continue
        rows.append({"family": "NULL_MODEL", "unit": os.path.basename(p), "test": "empirical",
                     "stat": _f(j.get("observed")), "p": pv,
                     "source": os.path.relpath(p, outdir), "sha256_16": _sha16(p)})
    return rows


def apply(rows, alpha=0.05):
    """族内 + 全局 BH。返回带 q 的新行。"""
    testable = [r for r in rows if r.get("p") is not None]
    if not testable:
        return [], {"n_hypotheses": 0}
    # 全局
    gq = bh([r["p"] for r in testable])
    for r, q in zip(testable, gq):
        r["q_global"] = q
    # 族内
    fams = {}
    for i, r in enumerate(testable):
        fams.setdefault(r["family"], []).append(i)
    for fam, idx in fams.items():
        qs = bh([testable[i]["p"] for i in idx])
        for i, q in zip(idx, qs):
            testable[i]["q_family"] = q
    for r in testable:
        r["significant_raw"] = int(r["p"] < alpha)
        r["significant_family"] = int(r.get("q_family", 1) < alpha)
        r["significant_global"] = int(r.get("q_global", 1) < alpha)
    summary = {
        "n_hypotheses": len(testable),
        "n_families": len({r["family"] for r in testable}),
        "alpha": alpha,
        "n_raw": sum(r["significant_raw"] for r in testable),
        "n_family": sum(r["significant_family"] for r in testable),
        "n_global": sum(r["significant_global"] for r in testable),
        "per_family": {f: {"n": len(idx),
                           "n_raw": sum(testable[i]["significant_raw"] for i in idx),
                           "n_family": sum(testable[i]["significant_family"] for i in idx)}
                       for f, idx in fams.items()},
        "method": "Benjamini–Hochberg（族内与全局各一次）",
        "n_non_p_endpoints": len(rows) - len(testable),
    }
    return testable, summary


def statement(summary):
    """生成 Methods 可直接引用的多重性声明（中英双语）。"""
    s = summary
    if not s.get("n_hypotheses"):
        return {"zh": "本流程无可校正的 p 值终点（未运行 MR/富集/基准）。",
                "en": "No p-value endpoints were available for multiplicity control."}
    zh = ("本流程共检测 %d 个统计假设（%d 类终点）。族内与全局均采用 Benjamini–Hochberg（BH）校正："
          "未校正时 %d 个达到 p<%.2g；族内校正后 %d 个（q_family<%.2g）；"
          "全局校正后 %d 个（q_global<%.2g）。**q_global 是唯一可写入结论的显著性口径**；"
          "仅 p<0.05 而未通过 q 的结果只能表述为\"提示性关联\"。"
          % (s["n_hypotheses"], s["n_families"], s["n_raw"], s["alpha"],
             s["n_family"], s["alpha"], s["n_global"], s["alpha"]))
    en = ("A total of %d statistical hypotheses across %d endpoint families were tested. "
          "Multiplicity was controlled with Benjamini–Hochberg correction both within family and "
          "globally: %d hypotheses reached p<%.2g before correction, %d after within-family "
          "correction (q_family<%.2g), and %d after global correction (q_global<%.2g). "
          "Only q_global<%.2g is reported as significant; results passing nominal p alone are "
          "described as suggestive."
          % (s["n_hypotheses"], s["n_families"], s["n_raw"], s["alpha"],
             s["n_family"], s["alpha"], s["n_global"], s["alpha"], s["alpha"]))
    return {"zh": zh, "en": en}


def write(outdir, testable, summary, out_sub="09_evidence"):
    d = os.path.join(outdir, out_sub)
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, "fdr_table.csv")
    cols = ["family", "unit", "test", "stat", "p", "q_family", "q_global",
            "significant_raw", "significant_family", "significant_global", "source", "sha256_16"]
    with open(p, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(testable, key=lambda x: (x["family"], x.get("q_global", 1))):
            w.writerow(r)
    st = statement(summary)
    summary["statement_zh"] = st["zh"]
    summary["statement_en"] = st["en"]
    summary["generated_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    json.dump(summary, open(os.path.join(d, "multiplicity_statement.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    return p, summary


# ============================ CLI ============================
def selftest():
    checks = []
    # 1) BH 与已知结果一致（R p.adjust("BH")）
    p = [0.01, 0.02, 0.03, 0.04, 0.05]
    q = [round(x, 6) for x in bh(p)]
    checks.append(("bh_known", q == [0.05] * 5, str(q)))
    # 2) 单值不变
    checks.append(("bh_single", [round(x, 6) for x in bh([0.03])] == [0.03], ""))
    # 3) 单调性：q 随 p 非降
    p2 = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]
    q2 = bh(p2)
    order = sorted(range(len(p2)), key=lambda i: p2[i])
    mono = all(q2[order[i]] <= q2[order[i + 1]] + 1e-12 for i in range(len(p2) - 1))
    checks.append(("bh_monotone", mono, ""))
    # 4) 全域 p=1 → q=1
    checks.append(("bh_ones", all(abs(x - 1) < 1e-12 for x in bh([1.0] * 5)), ""))
    # 5) Bonferroni ≥ BH
    qb, qbf = bh(p2), bonferroni(p2)
    checks.append(("bonf_ge_bh", all(a >= b - 1e-12 for a, b in zip(qbf, qb)), ""))
    # 6) 族内/全局分离：一族全显著、另一族全不显著
    rows = [{"family": "A", "unit": "a%d" % i, "test": "t", "stat": 1, "p": 1e-6,
             "source": "x", "sha256_16": "0"} for i in range(5)]
    rows += [{"family": "B", "unit": "b%d" % i, "test": "t", "stat": 1, "p": 0.4,
              "source": "x", "sha256_16": "0"} for i in range(5)]
    t, s = apply(rows)
    checks.append(("family_vs_global", s["n_family"] == 5 and s["n_global"] == 5,
                   "%d/%d" % (s["n_family"], s["n_global"])))
    # 7) 无 p 的终点不进校正
    rows2 = rows[:2] + [{"family": "C", "unit": "c", "test": "pp4", "stat": 0.9, "p": None,
                         "source": "x", "sha256_16": "0"}]
    t2, s2 = apply(rows2)
    checks.append(("pp4_not_in_bh", s2["n_hypotheses"] == 2 and s2["n_non_p_endpoints"] == 1,
                   "%d/%d" % (s2["n_hypotheses"], s2["n_non_p_endpoints"])))
    # 8) 声明含关键量
    st = statement(s)
    checks.append(("statement_numbers", "10" in st["zh"] and "q_global" in st["zh"], ""))
    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for n, g, d in checks:
        print("  [%s] %-22s %s" % ("✓" if g else "✗", n, d))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = argparse.ArgumentParser(description="np_fdr: 全局多重性控制与声明")
    sub = ap.add_subparsers(dest="fn")
    r = sub.add_parser("run")
    r.add_argument("--dir", required=True)
    r.add_argument("--alpha", type=float, default=0.05)
    r.add_argument("-o", "--outdir", default=None)
    m = sub.add_parser("mr")
    m.add_argument("--dir", required=True)
    m.add_argument("--alpha", type=float, default=0.05)
    s = sub.add_parser("statement")
    s.add_argument("--dir", required=True)
    a = ap.parse_args(argv)

    if a.fn == "run":
        rows = gather(a.dir)
        testable, summary = apply(rows, a.alpha)
        p, summary = write(a.dir, testable, summary, a.outdir or "09_evidence")
        print("[fdr] %d 个假设 / %d 族：p<%.2g 有 %d，族内 q 有 %d，全局 q 有 %d → %s" % (
            summary["n_hypotheses"], summary["n_families"], a.alpha, summary["n_raw"],
            summary["n_family"], summary["n_global"], p))
        return 0
    if a.fn == "mr":
        rows = [r for r in gather(a.dir) if r["family"] == "MR_IVW"]
        testable, summary = apply(rows, a.alpha)
        print("[fdr:mr] %d 个靶点：p<%.2g 有 %d；族内 q 有 %d" % (
            summary.get("n_hypotheses", 0), a.alpha, summary.get("n_raw", 0),
            summary.get("n_family", 0)))
        for r in sorted(testable, key=lambda x: x.get("q_family", 1)):
            print("  %-12s p=%.3g  q=%.3g  %s" % (
                r["unit"], r["p"], r.get("q_family", 1),
                "显著" if r["significant_family"] else "不显著"))
        return 0
    if a.fn == "statement":
        rows = gather(a.dir)
        _, summary = apply(rows, 0.05)
        st = statement(summary)
        print(st["zh"])
        print()
        print(st["en"])
        return 0
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
