#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_coloc_compare.py — 共定位近似 vs R `coloc` 的一致性核验（v3.0）

为什么需要（评审报告 §8、§9.3）
--------------------------------
`np_mr.coloc_abf()` 是 **ABF 近似**（W=0.15², p1=p2=1e-4, p12=1e-5），
它**不等于** `coloc::coloc.abf()`：官方实现用的是 Wakefield 近似贝叶斯因子 +
完整先验组合，且对输入做 MAF/SD 换算。在没对齐之前，Methods 只能写
"an ABF-based colocalisation approximation"，不能写 "we performed colocalisation"。

本模块做的是**核验**而不是"再实现一遍"：
1. 生成/读取同一批位点的暴露-结局汇总统计；
2. 用 np_mr 的 ABF 近似算 PP.H0–H4；
3. 调用 R `coloc::coloc.abf()` 算官方值（若 R 与 coloc 可用）；
4. 输出逐位点对照表 + 相关系数 + 偏差分布 + 判定（PASS/WARN/FAIL）。

R 不可用时**不伪造**：状态写 `R_UNAVAILABLE`，并把可直接粘贴运行的 R 脚本
写到输出目录，让用户在装有 R 的机器上补跑。

用法
----
    python np_coloc_compare.py make  --n 20 -o np_out/coloc_cmp       # 生成合成位点
    python np_coloc_compare.py run   --loci np_out/coloc_cmp/loci.tsv -o np_out/coloc_cmp
    python np_coloc_compare.py run   --exposure e.tsv --outcome o.tsv -o out/
    python np_coloc_compare.py selftest
"""

import argparse
import csv
import glob
import json
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import np_mr as MR  # noqa: E402

VERSION = "3.0"

R_SNIPPET = r'''# 由 np_coloc_compare.py 生成：在装有 R + coloc 的机器上运行
# install.packages("coloc")
library(coloc)
d <- read.delim("%(loci)s", stringsAsFactors = FALSE)
res <- lapply(split(d, d$locus), function(x) {
  e <- list(snp = x$snp, beta = x$beta_x, varbeta = x$se_x^2, type = "quant",
            N = %(n_x)d, sdY = 1)
  o <- list(snp = x$snp, beta = x$beta_y, varbeta = x$se_y^2, type = "quant",
            N = %(n_y)d, sdY = 1)
  tryCatch(coloc.abf(e, o)$summary, error = function(err) NULL)
})
out <- do.call(rbind, lapply(names(res), function(k) {
  if (is.null(res[[k]])) return(data.frame(locus = k, PP.H4 = NA))
  data.frame(locus = k, PP.H0 = res[[k]][["PP.H0.abf"]], PP.H1 = res[[k]][["PP.H1.abf"]],
             PP.H2 = res[[k]][["PP.H2.abf"]], PP.H3 = res[[k]][["PP.H3.abf"]],
             PP.H4 = res[[k]][["PP.H4.abf"]])
}))
write.csv(out, "%(out_csv)s", row.names = FALSE)
cat("R coloc done ->", "%(out_csv)s", "\n")
'''


def _csv(path):
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        return list(csv.DictReader(f, delimiter=delim))


def make_loci(outdir, n_loci=20, n_snp=40, seed=7, signal=True):
    """生成 n 个位点：一半有共享因果变异（H4），一半两个独立变异（H3）。"""
    rnd = random.Random(seed)
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, "loci.tsv")
    rows = []
    for li in range(n_loci):
        shared = (li % 2 == 0) and signal
        causal = rnd.sample(range(n_snp), 1 if shared else 2)
        for s in range(n_snp):
            b_x = rnd.gauss(0, 0.02)
            b_y = rnd.gauss(0, 0.02)
            if shared and s == causal[0]:
                b_x, b_y = 0.35, 0.28
            if (not shared) and len(causal) == 2:
                if s == causal[0]:
                    b_x = 0.35
                if s == causal[1]:
                    b_y = 0.30
            rows.append({"locus": "L%02d" % li, "snp": "rs%d_%d" % (li, s),
                         "beta_x": round(b_x, 5), "se_x": 0.02,
                         "beta_y": round(b_y, 5), "se_y": 0.02,
                         "p_x": 1e-6 if (shared and s == causal[0]) else 0.5})
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["locus", "snp", "beta_x", "se_x", "beta_y", "se_y", "p_x"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return p


def find_rscript():
    try:
        import np_tools
        for cand in ("Rscript", "Rscript.exe"):
            w = shutil.which(cand)
            if w:
                return w
    except Exception:
        pass
    for pat in (r"C:/Program Files/R/R-*/bin/Rscript.exe", "/usr/bin/Rscript",
                "/usr/local/bin/Rscript"):
        g = sorted(glob.glob(pat))
        if g:
            return g[-1]
    return None


def has_coloc(rscript, lib=None):
    code = ('cat(requireNamespace("coloc", quietly=TRUE))')
    env = dict(os.environ)
    if lib:
        env["R_LIBS_USER"] = lib
    try:
        r = subprocess.run([rscript, "-e", code], capture_output=True, text=True,
                           timeout=120, env=env)
        return (r.stdout or "").strip().endswith("TRUE")
    except Exception:
        return False


def our_pp(loci, n_x=1000, n_y=5000):
    """用 np_mr 的 ABF 近似逐位点算 PP。"""
    out = {}
    by = {}
    for r in loci:
        by.setdefault(r.get("locus") or "all", []).append(r)
    for loc, rows in by.items():
        snps = []
        for r in rows:
            try:
                snps.append({"snp": r["snp"], "beta_x": float(r["beta_x"]),
                             "se_x": float(r["se_x"]), "beta_y": float(r["beta_y"]),
                             "se_y": float(r["se_y"])})
            except (KeyError, TypeError, ValueError):
                continue
        res = MR.coloc_abf(snps)
        out[loc] = res
    return out


def run(loci_path, outdir, lib=None, n_x=1000, n_y=5000, s_prop=0.3):
    loci = _csv(loci_path)
    if not loci:
        return {"status": "NO_LOCI"}
    ours = our_pp(loci, n_x, n_y)
    os.makedirs(outdir, exist_ok=True)
    rscript = find_rscript()
    r_csv = os.path.join(outdir, "r_coloc_pp.csv")
    info = {"rscript": rscript, "coloc_available": False}
    if rscript:
        info["coloc_available"] = has_coloc(rscript, lib)
    status = "OK"
    if not (rscript and info["coloc_available"]):
        status = "R_UNAVAILABLE"
        script = os.path.join(outdir, "run_coloc.R")
        open(script, "w", encoding="utf-8").write(
            R_SNIPPET % {"loci": loci_path.replace("\\", "/"), "out_csv": r_csv.replace("\\", "/"),
                         "n_x": n_x, "n_y": n_y, "s_prop": s_prop})
        info["r_script"] = script
        info["hint"] = ("R 或 coloc 包不可用。两步：① 安装 R 与 coloc（install.packages('coloc')）；"
                        "② Rscript %s。之后重跑本命令即可得到对照表。" % os.path.basename(script))
    if status == "OK":
        # R 在非 UTF-8 locale 下按本地编码解析脚本：路径含中文会被截断成乱码。
        # 解决方式不是"让用户改路径"，而是把输入/输出/脚本都放到纯 ASCII 临时目录，
        # 跑完再把结果搬回 outdir。
        work = tempfile.mkdtemp(prefix="np_coloc_run_")
        try:
            loci_ascii = os.path.join(work, "loci.tsv")
            shutil.copyfile(loci_path, loci_ascii) if os.path.abspath(loci_path) != os.path.abspath(loci_ascii) else None
            r_csv_ascii = os.path.join(work, "r_coloc_pp.csv")
            script = os.path.join(work, "run_coloc.R")
            open(script, "w", encoding="utf-8").write(
                R_SNIPPET % {"loci": loci_ascii.replace("\\", "/"),
                             "out_csv": r_csv_ascii.replace("\\", "/"),
                             "n_x": n_x, "n_y": n_y, "s_prop": s_prop})
            env = dict(os.environ)
            if lib:
                env["R_LIBS_USER"] = lib
                env["R_LIBS"] = lib
            try:
                r = subprocess.run([rscript, script], capture_output=True, text=True,
                                   timeout=600, env=env, cwd=work)
                info["r_stdout"] = (r.stdout or "")[-500:]
                info["r_stderr"] = (r.stderr or "")[-500:]
                if r.returncode == 0 and os.path.exists(r_csv_ascii):
                    shutil.copyfile(r_csv_ascii, r_csv)
                else:
                    status = "R_FAILED"
            except Exception as e:
                status = "R_FAILED"
                info["error"] = "%s: %s" % (type(e).__name__, str(e)[:120])
            info["ascii_workdir"] = work
            # 把脚本也留一份在 outdir，便于用户在本地复跑
            shutil.copyfile(script, os.path.join(outdir, "run_coloc.R")) if os.path.exists(script) else None
        finally:
            shutil.rmtree(work, ignore_errors=True)
    pairs = []
    if status == "OK":
        for r in _csv(r_csv):
            loc = r.get("locus")
            o = ours.get(loc) or {}
            try:
                pairs.append({"locus": loc, "ours_pp4": float(o.get("PP.H4", 0)),
                              "r_pp4": float(r.get("PP.H4") or 0),
                              "ours_pp3": float(o.get("PP.H3", 0)),
                              "r_pp3": float(r.get("PP.H3") or 0)})
            except (TypeError, ValueError):
                continue
    verdict = "NOT_ASSESSED"
    stats = {}
    if pairs:
        xs = [p["r_pp4"] for p in pairs]
        ys = [p["ours_pp4"] for p in pairs]
        if len(xs) >= 3:
            mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
            num = sum((xs[i] - mx) * (ys[i] - my) for i in range(len(xs)))
            dx = math.sqrt(sum((v - mx) ** 2 for v in xs))
            dy = math.sqrt(sum((v - my) ** 2 for v in ys))
            rho = num / (dx * dy) if dx * dy else 0.0
            mae = sum(abs(xs[i] - ys[i]) for i in range(len(xs))) / len(xs)
            stats = {"n": len(xs), "pearson_r": round(rho, 4), "mae_pp4": round(mae, 4),
                     "max_abs_diff": round(max(abs(xs[i] - ys[i]) for i in range(len(xs))), 4)}
            # 判定：一致性只作为"能不能写 colocalisation"的门槛
            if rho >= 0.90 and mae <= 0.10:
                verdict = "PASS"
            elif rho >= 0.70:
                verdict = "WARN"
            else:
                verdict = "FAIL"
    res = {"status": status, "n_loci": len(ours), "comparison": pairs, "stats": stats,
           "verdict": verdict,
           "method_ours": "ABF 近似（W=0.15², p1=p2=1e-4, p12=1e-5）",
           "method_r": "coloc::coloc.abf（Wakefield ABF + 完整先验组合）",
           "wording_rule": ("PASS → Methods 可写 'colocalisation analysis (coloc.abf-equivalent)'；"
                            "WARN/FAIL/NOT_ASSESSED → 只能写 'an ABF-based colocalisation approximation'"),
           "info": info, "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")}
    json.dump(res, open(os.path.join(outdir, "coloc_comparison.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    if pairs:
        with open(os.path.join(outdir, "coloc_comparison.csv"), "w", encoding="utf-8-sig",
                  newline="") as f:
            w = csv.DictWriter(f, fieldnames=["locus", "ours_pp4", "r_pp4", "ours_pp3", "r_pp3"])
            w.writeheader()
            w.writerows(pairs)
    return res


def selftest():
    checks = []
    d = tempfile.mkdtemp(prefix="np_coloc_")
    p = make_loci(d, n_loci=6, n_snp=30)
    rows = _csv(p)
    checks.append(("make_loci", len(rows) == 6 * 30, "n=%d" % len(rows)))
    ours = our_pp(rows)
    checks.append(("our_pp_covers_loci", len(ours) == 6, "n=%d" % len(ours)))
    # 共享因果变异的位点，PP.H4 应显著高于独立变异位点（偶数 L00/L02/L04）
    h4_shared = [ours["L%02d" % i].get("PP.H4", 0) for i in (0, 2, 4)]
    h4_indep = [ours["L%02d" % i].get("PP.H4", 0) for i in (1, 3, 5)]
    checks.append(("pp4_direction", sum(h4_shared) / 3 > sum(h4_indep) / 3,
                   "%.3f vs %.3f" % (sum(h4_shared) / 3, sum(h4_indep) / 3)))
    checks.append(("pp_sums_to_one", all(abs(sum(ours[k][q] for q in
                   ("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")) - 1) < 1e-6 for k in ours), ""))
    r = run(p, d)
    checks.append(("status_honest", r["status"] in ("OK", "R_UNAVAILABLE", "R_FAILED"), r["status"]))
    checks.append(("script_written", os.path.exists(os.path.join(d, "run_coloc.R")), ""))
    checks.append(("wording_rule_present", "approximation" in r["wording_rule"], ""))
    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for n, g, dd in checks:
        print("  [%s] %-22s %s" % ("✓" if g else "✗", n, str(dd)[:70]))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = argparse.ArgumentParser(description="np_coloc_compare: ABF 近似 vs R coloc 一致性")
    sub = ap.add_subparsers(dest="fn")
    m = sub.add_parser("make")
    m.add_argument("--n", type=int, default=20)
    m.add_argument("--snp", type=int, default=40)
    m.add_argument("-o", "--outdir", default="np_out/coloc_cmp")
    r = sub.add_parser("run")
    r.add_argument("--loci", default=None)
    r.add_argument("--exposure", default=None)
    r.add_argument("--outcome", default=None)
    r.add_argument("--r-lib", default=None)
    r.add_argument("--n-exposure", type=int, default=1000)
    r.add_argument("--n-outcome", type=int, default=5000)
    r.add_argument("-o", "--outdir", default="np_out/coloc_cmp")
    a = ap.parse_args(argv)
    if a.fn == "make":
        p = make_loci(a.outdir, a.n, a.snp)
        print("[coloc] 合成位点 → %s" % p)
        return 0
    if a.fn == "run":
        loci = a.loci
        if not loci and a.exposure and a.outcome:
            ex = {r["SNP"]: r for r in _csv(a.exposure)}
            loci_o = []
            for r in _csv(a.outcome):
                s = r.get("SNP")
                if s in ex:
                    loci_o.append({"locus": "locus1", "snp": s, "beta_x": ex[s].get("beta"),
                                   "se_x": ex[s].get("se"), "beta_y": r.get("beta"),
                                   "se_y": r.get("se")})
            os.makedirs(a.outdir, exist_ok=True)
            loci = os.path.join(a.outdir, "loci.tsv")
            with open(loci, "w", encoding="utf-8", newline="") as f:
                w = csv.DictWriter(f, fieldnames=["locus", "snp", "beta_x", "se_x", "beta_y", "se_y"],
                                   delimiter="\t")
                w.writeheader()
                w.writerows(loci_o)
        if not loci:
            print("[coloc] 需要 --loci 或 --exposure/--outcome")
            return 2
        res = run(loci, a.outdir, a.r_lib, a.n_exposure, a.n_outcome)
        print("[coloc] %s | verdict=%s | %s" % (
            res["status"], res["verdict"], json.dumps(res.get("stats"), ensure_ascii=False)))
        if res["status"] == "R_UNAVAILABLE":
            print("  " + res["info"].get("hint", ""))
        return 0 if res["status"] in ("OK", "R_UNAVAILABLE") else 2
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
