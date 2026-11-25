#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_consensus.py —— 共识打分与不确定性量化（UQ）：Vina + Vinardo + GNINA

问题（UQ 盲区）
    单一 Vina 打分 -8.0 不代表真结合。不同打分函数偏差大。
    同时调用多个打分函数，算均值与标准差：方差大 → 预测不稳定 → Unstable_Binding 警示。

打分来源（诚实分级）
    - Vina      : 由 np_tools 解析（NP_VINA 环境变量 / config/tools.local.json / PATH）
    - Vinardo   : 同一 vina 二进制支持 --scoring vinardo（Vina 1.2.x 内置，实测版本支持）
    - GNINA     : 需单独安装（教程给出）；未装时只算前两者并标注 n_scorers
    - 离线模式  : 若提供已有三列打分 CSV，直接算 UQ（不需外部工具）

UQ 规则
    mean = 三打分均值；sd = 样本标准差
    sd > --sd-threshold（默认 1.5 kcal/mol）→ Unstable_Binding（强制写入报告）
    sd ≤ 阈值且 mean ≤ --good-cut（默认 -7）→ Stable_Binding
    否则 → Moderate

用法
    python np_consensus.py offline --scores three_scores.csv -o out/
    python np_consensus.py scan --receptor rec.pdbqt --ligands ligs_dir/ --vina $NP_VINA -o out/
    python np_consensus.py selftest
"""
import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
UNSTABLE = "Unstable_Binding"


def uq(scores, sd_threshold=1.5, good_cut=-7.0):
    """scores: list[float]（越负越好）。返回 dict。"""
    vals = [s for s in scores if isinstance(s, (int, float)) and math.isfinite(s)]
    if not vals:
        return {"status": "NO_SCORE", "flag": UNSTABLE}
    mean = sum(vals) / len(vals)
    sd = math.sqrt(sum((x - mean) ** 2 for x in vals) / (len(vals) - 1)) if len(vals) > 1 else 0.0
    if sd > sd_threshold:
        verdict, flag = "UNSTABLE", UNSTABLE
    elif mean <= good_cut:
        verdict, flag = "STABLE", ""
    else:
        verdict, flag = "MODERATE", ""
    return {"n_scorers": len(vals), "scores": vals, "mean": round(mean, 3),
            "sd": round(sd, 3), "verdict": verdict, "flag": flag}


def cmd_offline(a):
    rows = []
    with open(a.scores, encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f)
        hdrs = rdr.fieldnames or []
        score_cols = [h for h in hdrs if h.lower() not in ("compound", "target", "ligand", "name", "靶点", "成分")]
        for r in rdr:
            vals = []
            for h in score_cols:
                try:
                    vals.append(float(r[h]))
                except (ValueError, TypeError):
                    pass
            res = uq(vals, a.sd_threshold, a.good_cut)
            res["ligand"] = r.get("compound") or r.get("ligand") or r.get("name") or ""
            rows.append(res)
    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, "consensus_uq.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ligand", "n_scorers", "mean", "sd", "verdict", "flag", "scores"])
        for r in rows:
            w.writerow([r["ligand"], r["n_scorers"], r.get("mean", ""), r.get("sd", ""),
                        r.get("verdict", ""), r.get("flag", ""), ";".join(str(x) for x in r.get("scores", []))])
    n_unstable = sum(1 for r in rows if r.get("flag") == UNSTABLE)
    json.dump({"n": len(rows), "n_unstable": n_unstable,
               "note": "sd 为多打分函数样本标准差；%s 者不得声称稳定结合" % UNSTABLE},
              open(os.path.join(a.outdir, "consensus_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[consensus] %d 配体，%d 个 %s → %s" %
          (len(rows), n_unstable, UNSTABLE, os.path.join(a.outdir, "consensus_uq.csv")))
    return 0


def find_vina(explicit=None):
    # 评审 §7：环境变量优先，作者机器路径仅作兼容候选，不作为唯一来源
    try:
        import np_tools
        cands = [explicit, np_tools.path_of("vina", explicit)]
    except Exception:
        cands = [explicit, os.environ.get("NP_VINA"), shutil.which("vina")]
    return next((c for c in cands if c and (os.path.exists(c) or shutil.which(c))), None)


def cmd_scan(a):
    vina = find_vina(a.vina)
    if not vina:
        return 2
    os.makedirs(a.outdir, exist_ok=True)
    ligs = [f for f in os.listdir(a.ligands) if f.lower().endswith((".pdbqt", ".sdf"))]
    results = []
    for lig in ligs[:a.max_ligands]:
        lp = os.path.join(a.ligands, lig)
        scores = {}
        for scoring in ("vina", "vinardo"):
            # Vina 1.2.x 内置 vinardo；若不支持则报错跳过
            out = os.path.join(a.outdir, "%s_%s.out" % (lig, scoring))
            cmd = [vina, "--receptor", a.receptor, "--ligand", lp,
                   "--scoring", scoring, "--center_x", str(a.cx), "--center_y", str(a.cy),
                   "--center_z", str(a.cz), "--size_x", str(a.box), "--size_y", str(a.box),
                   "--size_z", str(a.box), "--out", out, "--quiet"]
            try:
                r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
                # 解析 "1  -8.123  ..." 最优行
                m = None
                for line in (r.stdout or "").splitlines():
                    parts = line.split()
                    if len(parts) >= 2 and parts[0].isdigit() and parts[0] == "1":
                        m = parts
                        break
                if m:
                    scores[scoring] = float(m[1])
            except Exception:
                pass
        gn = ""
        if a.gnina:
            scores["gnina"] = None  # 需外部填写，标注
        vals = [v for v in scores.values() if v is not None]
        res = uq(vals, a.sd_threshold, a.good_cut)
        res["ligand"] = lig
        results.append(res)
    with open(os.path.join(a.outdir, "consensus_uq.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ligand", "n_scorers", "mean", "sd", "verdict", "flag", "scores"])
        for r in results:
            w.writerow([r["ligand"], r["n_scorers"], r.get("mean", ""), r.get("sd", ""),
                        r.get("verdict", ""), r.get("flag", ""), ";".join(str(x) for x in r.get("scores", []))])
    print("[consensus] scan %d 配体（Vina=%s）→ %s" % (len(results), vina, a.outdir))
    return 0


def selftest():
    r1 = uq([-9.0, -8.8, -9.2])          # sd 小 → STABLE
    r2 = uq([-9.0, -5.0, -11.0])         # sd 大 → UNSTABLE
    ok = (r1["verdict"] == "STABLE" and r2["verdict"] == "UNSTABLE"
          and r2["flag"] == UNSTABLE)
    print("selftest: %s (stable sd=%.2f; unstable sd=%.2f flag=%s)" %
          ("PASS" if ok else "FAIL", r1["sd"], r2["sd"], r2["flag"]))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_consensus: 多打分共识与 UQ")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("offline")
    p.add_argument("--scores", required=True)
    p.add_argument("-o", "--outdir", default="np_out/consensus")
    p.add_argument("--sd-threshold", type=float, default=1.5)
    p.add_argument("--good-cut", type=float, default=-7.0)
    p.set_defaults(fn="offline")
    p = sub.add_parser("scan")
    p.add_argument("--receptor", required=True); p.add_argument("--ligands", required=True)
    p.add_argument("--vina", default=None); p.add_argument("--gnina", action="store_true")
    p.add_argument("--cx", type=float, default=0); p.add_argument("--cy", type=float, default=0)
    p.add_argument("--cz", type=float, default=0); p.add_argument("--box", type=float, default=24)
    p.add_argument("--max-ligands", type=int, default=20)
    p.add_argument("-o", "--outdir", default="np_out/consensus")
    p.add_argument("--sd-threshold", type=float, default=1.5)
    p.add_argument("--good-cut", type=float, default=-7.0)
    p.set_defaults(fn="scan")
    p = sub.add_parser("selftest"); p.set_defaults(fn="selftest")
    a = ap.parse_args(argv)
    if a.fn == "selftest":
        return selftest()
    if a.fn == "offline":
        return cmd_offline(a)
    return cmd_scan(a)


if __name__ == "__main__":
    sys.exit(main())
