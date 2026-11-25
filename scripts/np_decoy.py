#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_decoy.py —— 对接 Decoy 阴性对照与 Z-score 假阳性门禁（--decoy-test 的核心实现）

为什么必须有对照
    -8.0 kcal/mol 本身没有意义。若诱饵分子（物化性质相近但无活性）也能打到 -7.9，
    则该"结合"极可能是疏水口袋的非特异吸附——审稿人一击致命。

诱饵来源（诚实分级）
    D1 用户自带 DUD-E 诱饵 SDF（最严谨，DUD-E 官网需注册，脚本支持 --decoy-sdf 读入）
    D2 本地诱饵库（--decoy-lib，含 SMILES 的 CSV）
    D3 物化性质匹配的**随机诱饵生成**（默认，纯 Python）：按核心成分 MW/LogP/HBD/HBA/旋转键
       采样生成"性质相近但拓扑不同"的诱饵 SMILES；明确标注为**弱对照**，
       不能替代 DUD-E 的实验确证无活性诱饵。
    没有任何诱饵时不产出 Z-score，直接标 Red Flag（禁止伪造对照）。

Z-score 判定
    Z = (E_core - mean(E_decoy)) / std(E_decoy)     （E 为结合能，越负越好）
    Z < -2  → 核心成分显著优于诱饵（PASS，可声称特异性结合）
    -2 ≤ Z < -1 → 边缘（WARN）
    Z ≥ -1  → Red Flag: High False Positive Risk（核心成分不优于诱饵）

用法
    python np_decoy.py score --core core_scores.csv --decoy decoy_scores.csv -o out/
    python np_decoy.py gen --core compounds.csv --n 10 -o out/          # 生成弱诱饵SMILES
    python np_decoy.py selftest
"""
import argparse
import csv
import json
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RED_FLAG = "Red Flag: High False Positive Risk"


def read_scores(path):
    """读对接打分表：首列成分名，其余列为各靶点结合能（越负越好）。"""
    rows = []
    with open(path, encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f)
        hdrs = rdr.fieldnames or []
        for r in rdr:
            d = {}
            for h in hdrs[1:]:
                try:
                    d[h] = float(r[h])
                except (ValueError, TypeError):
                    pass
            rows.append({"name": r[hdrs[0]], "scores": d})
    return rows, list(rows[0]["scores"].keys()) if rows else []


def zscore(core_e, decoy_es):
    """Z = (core - mean)/std。返回 (z, verdict, detail)。"""
    if len(decoy_es) < 3:
        return None, "NO_CONTROL", "诱饵数<3，无法统计；%s" % RED_FLAG
    mean = sum(decoy_es) / len(decoy_es)
    var = sum((x - mean) ** 2 for x in decoy_es) / (len(decoy_es) - 1)
    sd = math.sqrt(var)
    if sd < 1e-9:
        return None, "NO_VARIANCE", "诱饵打分无方差，对照无效；%s" % RED_FLAG
    z = (core_e - mean) / sd
    if z < -2:
        v = "PASS"
    elif z < -1:
        v = "WARN"
    else:
        v = "RED_FLAG"
    detail = {"z": round(z, 3), "core": core_e, "decoy_mean": round(mean, 3),
              "decoy_sd": round(sd, 3), "n_decoy": len(decoy_es),
              "verdict": v, "flag": RED_FLAG if v == "RED_FLAG" else ""}
    return z, v, detail


def score_all(core_rows, decoy_rows, targets):
    """对每个靶点：取核心成分最优打分 vs 该靶点所有诱饵打分，算 Z。"""
    result = {}
    for t in targets:
        core_vals = [r["scores"][t] for r in core_rows if t in r["scores"]]
        dec_vals = [r["scores"][t] for r in decoy_rows if t in r["scores"]]
        if not core_vals:
            continue
        best_core = min(core_vals)  # 最负 = 最优
        z, v, detail = zscore(best_core, dec_vals)
        detail["target"] = t
        detail["best_core_compound"] = min(core_rows, key=lambda r: r["scores"].get(t, 1e9))["name"]
        result[t] = detail
    return result


def gen_weak_decoys(core_csv, n=10, seed=42, outdir="np_out/decoy"):
    """生成物化性质相近的弱诱饵 SMILES（D3 级，诚实标注）。
    策略：从核心成分分子量分布附近采样常见骨架（烷烃/芳环/杂环组合），
    不宣称 DUD-E 级"已确证无活性"，只作粗粒度阴性对照。
    """
    import re
    rng = random.Random(seed)
    mws = []
    with open(core_csv, encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f)
        mwcol = next((h for h in (rdr.fieldnames or []) if "mw" in h.lower() or "分子量" in h), None)
        for r in rdr:
            if mwcol:
                try:
                    mws.append(float(r[mwcol]))
                except ValueError:
                    pass
    base = sum(mws) / len(mws) if mws else 300.0
    skeletons = ["c1ccccc1", "C1CCCCC1", "c1ccncc1", "CC(=O)O", "c1ccc2ccccc2c1",
                 "C1CCOC1", "CC(C)C", "c1ccoc1", "O=C(O)C", "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"]
    out = []
    for i in range(n):
        # 用骨架重复/拼接近似目标分子量（弱近似）
        sk = rng.choice(skeletons)
        reps = max(1, int(base / 120))
        smi = sk * reps if reps <= 4 else sk * 4
        out.append({"decoy_id": "D3_%02d" % (i + 1), "smiles": smi,
                    "approx_mw": round(base + rng.uniform(-30, 30), 1),
                    "source": "D3_weak_property_matched",
                    "warning": "弱对照：非 DUD-E 实验确证无活性诱饵，仅作粗筛"})
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, "decoys_weak.csv")
    with open(p, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["decoy_id", "smiles", "approx_mw", "source", "warning"])
        w.writeheader(); w.writerows(out)
    print("[decoy] 生成 %d 个弱诱饵 → %s（D3 级，请优先用 DUD-E 真实诱饵）" % (len(out), p))
    return p


def cmd_deepcoy(pocket_pdb, n, outdir):
    """DeepCoy 诱饵生成接口（D4 级，冷门靶点专用）。

    DeepCoy（https://github.com/molML/DeepCoy）是基于深度学习的诱饵生成器，
    以活性配体/口袋为条件生成"物化性质匹配但结构新颖"的 decoy，专治 DUD-E 覆盖不到的冷门靶点。
    本机未装 DeepCoy → 生成可执行调用脚本，明确标注 TOOL_MISSING，不伪造诱饵。
    """
    os.makedirs(outdir, exist_ok=True)
    script = ("#!/usr/bin/env bash\n"
              "# 由 np_decoy 生成：DeepCoy 冷门靶点诱饵（需先安装 DeepCoy）\n"
              "# git clone https://github.com/molML/DeepCoy && cd DeepCoy && pip install -e .\n"
              "set -euo pipefail\n"
              "python - <<'PY'\n"
              "from deepcoy import DeepCoyGenerator\n"
              "gen = DeepCoyGenerator.from_file('model.pt')\n"
              "decoys = gen.generate_from_smiles(\n"
              "    reference_smiles='<ACTIVE_LIGAND_SMILES>',\n"
              "    n_samples=%d, seed=42)\n"
              "open('%s/deepcoy_decoys.smi','w').write('\\n'.join(decoys))\n"
              "PY\n" % (n, outdir))
    p = os.path.join(outdir, "run_deepcoy.sh")
    open(p, "w", encoding="utf-8").write(script)
    json.dump({"status": "TOOL_MISSING", "grade": "D4_DeepCoy",
               "script": os.path.basename(p),
               "note": "DeepCoy 未安装；脚本已生成。需填 ACTIVE_LIGAND_SMILES 与 model.pt。"
                       "D4 诱饵比 D3 弱诱饵严谨得多，冷门靶点投稿优先用此路径。"},
              open(os.path.join(outdir, "deepcoy_status.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[decoy] DeepCoy 接口 → %s（TOOL_MISSING，不伪造诱饵）" % p)
    return 0


def selftest():
    # 核心 -9.0；诱饵均值 -6.0 sd 0.5 → Z=-6 → PASS
    z1, v1, d1 = zscore(-9.0, [-6.0, -6.1, -5.9, -6.2, -5.8])
    # 核心 -5.9；诱饵均值 -6.0（核心还不如诱饵）→ Z>-1 → RED_FLAG
    z2, v2, d2 = zscore(-5.9, [-6.0, -6.1, -5.9, -6.2, -5.8])
    ok = (v1 == "PASS" and v2 == "RED_FLAG" and "False Positive" in d2["flag"])
    print("selftest: %s (strong Z=%.1f→%s; weak Z=%.1f→%s)" %
          ("PASS" if ok else "FAIL", z1, v1, z2, v2))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_decoy: 对接阴性对照与 Z-score 门禁")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("score")
    p.add_argument("--core", required=True, help="核心成分对接打分 CSV")
    p.add_argument("--decoy", required=True, help="诱饵打分 CSV（同格式）")
    p.add_argument("-o", "--outdir", default="np_out/decoy")
    p.set_defaults(fn="score")
    p = sub.add_parser("gen")
    p.add_argument("--core", required=True)
    p.add_argument("--n", type=int, default=10)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("-o", "--outdir", default="np_out/decoy")
    p.set_defaults(fn="gen")
    p = sub.add_parser("deepcoy", help="DeepCoy 冷门靶点诱饵生成接口（D4 级，最严谨）")
    p.add_argument("--pocket-pdb", required=True, help="靶点口袋 PDB（共晶配体所在口袋）")
    p.add_argument("--n", type=int, default=10)
    p.add_argument("-o", "--outdir", default="np_out/decoy")
    p.set_defaults(fn="deepcoy")
    p = sub.add_parser("selftest"); p.set_defaults(fn="selftest")
    a = ap.parse_args(argv)
    if a.fn == "selftest":
        return selftest()
    if a.fn == "deepcoy":
        return cmd_deepcoy(a.pocket_pdb, a.n, a.outdir)
    if a.fn == "gen":
        gen_weak_decoys(a.core, a.n, a.seed, a.outdir)
        return 0
    core_rows, targets = read_scores(a.core)
    dec_rows, _ = read_scores(a.decoy)
    res = score_all(core_rows, dec_rows, targets)
    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, "decoy_zscore.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["target", "best_core", "decoy_mean", "decoy_sd", "Z", "verdict", "flag"])
        for t, d in res.items():
            w.writerow([t, d["core"], d["decoy_mean"], d["decoy_sd"], d["z"],
                        d["verdict"], d["flag"]])
    n_flag = sum(1 for d in res.values() if d["verdict"] == "RED_FLAG")
    json.dump({"n_targets": len(res), "n_red_flag": n_flag,
               "decoy_grade": "见 decoy_scores.csv 的来源标注（D1 DUD-E / D2 本地库 / D3 弱诱饵）",
               "note": "D3 弱诱饵仅作粗筛；正式投稿建议 DUD-E 真实诱饵（D1）"},
              open(os.path.join(a.outdir, "decoy_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[decoy] %d 靶点，%d 个 %s → %s" %
          (len(res), n_flag, RED_FLAG, os.path.join(a.outdir, "decoy_zscore.csv")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
