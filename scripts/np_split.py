#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_split.py — 基准切分（随机 / 复合物-家族 / 留一数据库 LODO）与外部工具对照（v3.0）

为什么需要（评审报告 §2，七寸之二）
------------------------------------
此前的基准只有"一个随机靶点池"。文献明确指出：网络传播/靶点预测任务的**普通随机切分**
会因蛋白复合物与家族成员高度相关而系统性乐观（PLOS Comput Biol 2019, PMID 31479437；
Pahikkala 2014, PMID 24723570）。审稿人只要问一句"你的切分是怎么做的"，随机切分就撑不住。

本模块给三种切分，并**明确各自的解释边界**：

1. `random`    —— 阳性随机分折。**乐观**，只作参照，不得单独作为性能证据。
2. `complex`   —— 按"复合物/家族"整组留出（leave-one-group-out）。
                  没有 CORUM 时用**基因家族前缀**作代理，输出里标 `proxy=gene_family`，
                  绝不冒充真实复合物。
3. `lodo`      —— 留一数据库：用库 A 建的排序去预测库 B 的金标准（反之亦然），
                  前提是 A、B 不是子集关系（会先做包含检验并报警）。

注意（必须写进 Methods）
------------------------
对一个**固定排序**（不是训练出来的模型）做切分，测的不是"泛化误差"，而是
"性能是否集中在少数家族/数据库上"。因此本模块输出的是**每折指标的分布**与
`concentration`（前 20% 家族贡献的命中占比），而不是一个漂亮的单一 AUC。

用法
----
    python np_split.py run --edges pred.csv --gold-csv gold.csv --split complex -o out/
    python np_split.py run --edges pred.csv --gold-csv gold.csv --split random  -o out/
    python np_split.py lodo --db-edges all_edges.csv --db-col database -o out/
    python np_split.py compare --gold-csv gold.csv --pred A=a.csv --pred B=b.csv -o out/
    python np_split.py selftest
"""

import argparse
import csv
import glob
import json
import math
import os
import random
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import np_benchmark as B  # noqa: E402

VERSION = "3.0"

# 基因家族前缀（复合物代理）。只做"同家族→同组"的粗分组，不做生物学断言。
FAMILY_PREFIX = [
    "MAPK", "CASP", "BCL", "BAX", "BAK", "AKT", "PIK3", "PRK", "PLCG", "STAT",
    "JAK", "SRC", "FYN", "LYN", "YES", "ABL", "EGFR", "ERBB", "FGFR", "IGF1R",
    "TNFRSF", "TRAF", "IRAK", "MYD", "NFKB", "RELA", "RELB", "NFKBIA", "IKBK",
    "CCL", "CXCL", "CXCR", "CCR", "IL", "IFN", "TLR", "NLR", "MMP", "ADAM",
    "COL", "ITGA", "ITGB", "CDH", "CTNN", "WNT", "FZD", "SMAD", "TGFB", "BMP",
    "GATA", "FOX", "SOX", "HOX", "RUNX", "ETS", "MYC", "JUN", "FOS", "ATF",
    "NR", "PPAR", "ESR", "AR", "VDR", "RXR", "RAR", "THR", "HIF", "EPAS",
    "SLC", "ABC", "CYP", "UGT", "SULT", "GST", "NQO", "HMOX", "NOS", "PTGS",
    "ALOX", "PLA2", "LTA", "LTB", "HRH", "ADRB", "CHRM", "DRD", "HTR", "OPR",
    "GABR", "GRIN", "GRIA", "SCN", "KCN", "CACNA", "RYR", "ATP", "VDAC",
    "HSP", "DNAJ", "CALR", "CANX", "PDIA", "UBB", "UBC", "PSM", "USP", "OTU",
    "E2F", "CDK", "CCNA", "CCNB", "CCND", "CCNE", "CDC", "CHEK", "ATM", "ATR",
    "BRCA", "RAD", "XRCC", "MLH", "MSH", "POLD", "POLE", "TOP", "H2A", "H2B",
    "H3", "H4", "HDAC", "SIRT", "KAT", "KDM", "DNMT", "TET",
]


def family_of(gene):
    """基因符号 → 家族代理。返回 (family, method)。"""
    g = (gene or "").strip().upper()
    for p in sorted(FAMILY_PREFIX, key=len, reverse=True):
        if g.startswith(p):
            return p + " family", "prefix"
    m = re.match(r"^([A-Z]+)[0-9]+$", g)
    if m and len(m.group(1)) >= 3:
        return m.group(1) + " family", "prefix"
    return g, "singleton"


def load_groups(map_file=None):
    """可选外部复合物表：gene,complex。缺省返回空（用家族代理）。"""
    out = {}
    if map_file and os.path.exists(map_file):
        with open(map_file, encoding="utf-8-sig", errors="replace") as f:
            first = f.readline()
            delim = "\t" if first.count("\t") > first.count(",") else ","
            f.seek(0)
            for r in csv.DictReader(f, delimiter=delim):
                g = (r.get("gene") or r.get("symbol") or "").strip().upper()
                c = (r.get("complex") or r.get("group") or "").strip()
                if g and c:
                    out[g] = c
    return out


def _folds_random(genes, y, k, seed):
    pos = [i for i, v in enumerate(y) if v == 1]
    neg = [i for i, v in enumerate(y) if v == 0]
    rnd = random.Random(seed)
    rnd.shuffle(pos)
    rnd.shuffle(neg)
    folds = [[] for _ in range(k)]
    for arr in (pos, neg):
        for j, i in enumerate(arr):
            folds[j % k].append(i)
    return folds


def _folds_group(genes, y, k, group_of, seed):
    """按组整组留出：同组阳性必在同一折（这才是防泄漏的关键）。"""
    groups = {}
    for i, g in enumerate(genes):
        if y[i] == 1:
            groups.setdefault(group_of.get(g, g), []).append(i)
    keys = sorted(groups)
    rnd = random.Random(seed)
    rnd.shuffle(keys)
    folds = [[] for _ in range(k)]
    load = [0] * k
    for key in keys:
        j = min(range(k), key=lambda x: load[x])
        folds[j].extend(groups[key])
        load[j] += len(groups[key])
    neg = [i for i, v in enumerate(y) if v == 0]
    rnd.shuffle(neg)
    for j, i in enumerate(neg):
        folds[j % k].append(i)
    return folds


def evaluate_folds(ranked, pos, amb, folds, topk=(20,)):
    """每折：阳性只取该折内的（其余为阴性）；返回逐折指标 + 汇总。"""
    genes = [g for g, _ in ranked]
    sc = dict(ranked)
    per = []
    for fi, idx in enumerate(folds):
        fold_pos = {genes[i] for i in idx if genes[i] in pos}
        if not fold_pos:
            continue
        r = B.evaluate(ranked, fold_pos, amb, topk=topk, n_boot=0)
        if r.get("status") != "OK":
            per.append({"fold": fi, "status": r["status"], "n_pos": len(fold_pos)})
            continue
        per.append({"fold": fi, "n_pos": len(fold_pos), "auc": r["AUC"],
                    "ap": r["AUPRC"],
                    "topk_hits": {k: v["hits"] for k, v in r["at_k"].items()},
                    "topk_precision": {k: v["precision"] for k, v in r["at_k"].items()}})
    ok = [p for p in per if p.get("auc") is not None]
    if not ok:
        return {"status": "NO_VALID_FOLD", "folds": per}
    aucs = [p["auc"] for p in ok]
    mu = sum(aucs) / len(aucs)
    sd = math.sqrt(sum((a - mu) ** 2 for a in aucs) / max(1, len(aucs) - 1)) if len(aucs) > 1 else 0.0
    return {"status": "OK", "n_folds": len(ok), "auc_mean": round(mu, 4),
            "auc_sd": round(sd, 4), "auc_min": round(min(aucs), 4),
            "auc_max": round(max(aucs), 4), "folds": per}


def concentration(ranked, pos, top_frac=0.2):
    """命中集中度：按家族汇总命中数，前 top_frac 家族贡献的命中占比。

    高集中度 = "表面性能由少数家族撑起来"，这比一个 AUC 更能说明问题。
    """
    hits = [g for g, _ in ranked if g in pos]
    fam = {}
    for g in hits:
        f, _ = family_of(g)
        fam[f] = fam.get(f, 0) + 1
    if not fam:
        return {"n_hits": 0}
    vals = sorted(fam.values(), reverse=True)
    k = max(1, int(math.ceil(len(vals) * top_frac)))
    return {"n_hits": len(hits), "n_families": len(fam),
            "top_families": sorted(fam.items(), key=lambda x: -x[1])[:8],
            "top_frac": top_frac, "concentration": round(sum(vals[:k]) / len(hits), 4)}


def run(edges=None, pred=None, gold_csv=None, split="random", k=5, seed=42,
        topk=(20,), map_file=None, outdir="np_out/benchmark_split", label=None):
    ranked = B.load_ranking(pred) if pred else B.rank_from_edges(edges)
    if not ranked:
        return {"status": "EMPTY_RANKING"}
    rows, hdrs = B.read_table(gold_csv)
    low = {h.lower(): h for h in hdrs}
    gc = next((low[c] for c in B.GENE_COLS if c in low), hdrs[0])
    tc = next((low[c] for c in ("tier", "evidence", "level", "class") if c in low), None)
    pos, amb = set(), set()
    for r in rows:
        g = (r.get(gc) or "").strip().upper()
        if not g:
            continue
        t = (r.get(tc) or "core").strip().lower() if tc else "core"
        (amb if t.startswith("ext") else pos).add(g)
    genes = [g for g, _ in ranked]
    y = [1 if g in pos else 0 for g in genes]
    ext_map = load_groups(map_file)
    group_of, proxy = {}, "gene_family"
    for g in genes:
        if g in ext_map:
            group_of[g] = ext_map[g]
            proxy = "provided_complex_map"
        else:
            group_of[g], _ = family_of(g)
    if split == "random":
        folds = _folds_random(genes, y, k, seed)
    elif split == "complex":
        folds = _folds_group(genes, y, k, group_of, seed)
    else:
        return {"status": "BAD_SPLIT", "split": split}
    res = evaluate_folds(ranked, pos, amb, folds, topk=topk)
    res.update({"split": split, "proxy": proxy if split == "complex" else "none",
                "k": k, "seed": seed, "label": label or os.path.basename(gold_csv or ""),
                "n_ranked": len(ranked), "n_pos": len(pos),
                "concentration": concentration(ranked, pos),
                "interpretation": (
                    "对固定排序做留组评估：测的是性能是否集中在少数家族，不是训练泛化误差"
                    if split == "complex" else
                    "随机切分：阳性高度相关（复合物/家族），指标偏乐观，仅作参照"),
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")})
    os.makedirs(outdir, exist_ok=True)
    tag = "%s_%s" % (split, (label or "gold").replace(" ", "_")[:24])
    json.dump(res, open(os.path.join(outdir, "split_%s.json" % tag), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    return res


def _db_edges(path, db_col="database"):
    """读带数据库标签的边表：返回 {db: set((compound,gene))}。"""
    rows, hdrs = B.read_table(path)
    low = {h.lower(): h for h in hdrs}
    gc = next((low[c] for c in B.GENE_COLS if c in low), None)
    cc = next((low[c] for c in ("compound", "mol_id", "molecule", "成分") if c in low), None)
    dc = low.get(db_col.lower())
    if not (gc and dc):
        return {}, "缺 gene 列或 %s 列" % db_col
    out = {}
    for r in rows:
        g = (r.get(gc) or "").strip().upper()
        d = (r.get(dc) or "").strip() or "unknown"
        c = (r.get(cc) or "_").strip().lower() if cc else "_"
        if g:
            out.setdefault(d, set()).add((c, g))
    return out, ""


def lodo(path, db_col="database", outdir="np_out/benchmark_split"):
    """留一数据库：库 d 的边作金标准，其余库的边建排序，互相预测。"""
    dbs, err = _db_edges(path, db_col)
    if err:
        return {"status": "BAD_INPUT", "reason": err}
    if len(dbs) < 2:
        return {"status": "TOO_FEW_DBS", "dbs": list(dbs)}
    results, subset = [], {}
    names = sorted(dbs)
    for a in names:
        for b in names:
            if a == b:
                continue
            A, Bs = dbs[a], dbs[b]
            inter = len(A & Bs)
            subset["%s⊂%s" % (a, b)] = round(inter / max(1, len(A)), 3)
    for held in names:
        train = set().union(*[dbs[x] for x in names if x != held])
        sup = {}
        for c, g in train:
            sup.setdefault(g, set()).add(c)
        ranked = sorted([(g, len(cs)) for g, cs in sup.items()], key=lambda x: (-x[1], x[0]))
        gold = {g for _, g in dbs[held]}
        r = B.evaluate(ranked, gold, set(), topk=(10, 20, 50), n_boot=0)
        results.append({"held_out_db": held, "n_train_edges": len(train),
                        "n_holdout_edges": len(dbs[held]), **{k: r.get(k) for k in
                        ("status", "AUC", "AUPRC", "n_evaluated", "n_positive_found", "coverage")},
                        "at_k": r.get("at_k")})
    suspicious = [k for k, v in subset.items() if v > 0.95]
    res = {"status": "OK", "split": "lodo", "databases": names, "n_edges": {k: len(v) for k, v in dbs.items()},
           "pairwise_containment": subset, "results": results,
           "subset_warning": (("疑似子集关系（包含率>0.95）：%s → 该对照**不是独立验证**，"
                               "须先查清两库谱系" % "；".join(suspicious)) if suspicious else ""),
           "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")}
    os.makedirs(outdir, exist_ok=True)
    json.dump(res, open(os.path.join(outdir, "split_lodo.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    return res


def compare(gold_csv, preds, topk=(20,), outdir="np_out/benchmark_split"):
    """外部工具对照：同一金标准、同一候选池、同一指标，逐个工具出数。"""
    rows, hdrs = B.read_table(gold_csv)
    low = {h.lower(): h for h in hdrs}
    gc = next((low[c] for c in B.GENE_COLS if c in low), hdrs[0])
    tc = next((low[c] for c in ("tier", "evidence", "level", "class") if c in low), None)
    pos, amb = set(), set()
    for r in rows:
        g = (r.get(gc) or "").strip().upper()
        if not g:
            continue
        t = (r.get(tc) or "core").strip().lower() if tc else "core"
        (amb if t.startswith("ext") else pos).add(g)
    out = []
    for name, path in preds.items():
        ranked = B.load_ranking(path)
        if not ranked:
            out.append({"tool": name, "status": "EMPTY"})
            continue
        r = B.evaluate(ranked, pos, amb, topk=topk, n_boot=200)
        out.append({"tool": name, "status": r.get("status"),
                    "n_ranked": r.get("n_ranked"), "AUC": r.get("AUC"), "AUPRC": r.get("AUPRC"),
                    "coverage": r.get("coverage"),
                    "topk_hits": {k: v["hits"] for k, v in (r.get("at_k") or {}).items()},
                    "ef": {k: v["enrichment_factor"] for k, v in (r.get("at_k") or {}).items()}})
    res = {"status": "OK", "gold": os.path.basename(gold_csv), "n_pos": len(pos),
           "tools": out,
           "note": ("所有工具必须用同一基因命名、同一候选池；覆盖率不同会导致 Recall 不可比，"
                    "因此同时报 top-K 命中数与 EF。"),
           "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")}
    os.makedirs(outdir, exist_ok=True)
    json.dump(res, open(os.path.join(outdir, "external_comparison.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    with open(os.path.join(outdir, "external_comparison.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tool", "status", "n_ranked", "AUC", "AUPRC", "coverage",
                    "hits@20", "EF@20"])
        for t in out:
            w.writerow([t.get("tool"), t.get("status"), t.get("n_ranked"), t.get("AUC"),
                        t.get("AUPRC"), t.get("coverage"),
                        (t.get("topk_hits") or {}).get("20", ""),
                        (t.get("ef") or {}).get("20", "")])
    return res


def selftest():
    checks = []
    # 1) 家族代理：同家族必须同组，数字后缀不影响
    f1, m1 = family_of("MAPK14")
    f2, _ = family_of("MAPK1")
    checks.append(("family_mapk", f1 == f2 and m1 == "prefix", "%s|%s" % (f1, f2)))
    f3, _ = family_of("TP53")
    checks.append(("family_fallback", f3 is not None, f3))

    # 2) 分组切分：同组阳性不得跨折
    genes = ["MAPK1", "MAPK3", "AKT1", "AKT2", "IL6", "TNF", "BCL2", "BAX"]
    y = [1, 1, 1, 1, 1, 1, 1, 1]
    group_of = {g: family_of(g)[0] for g in genes}
    folds = _folds_group(genes, y, 3, group_of, 42)
    fold_of = {genes[i]: fi for fi, idx in enumerate(folds) for i in idx}
    same = all(fold_of["MAPK1"] == fold_of["MAPK3"] for _ in (0,))
    checks.append(("group_leakage_blocked", same, "MAPK1=%s MAPK3=%s" % (
        fold_of["MAPK1"], fold_of["MAPK3"])))
    checks.append(("all_assigned", len({i for f in folds for i in f}) == len(genes), ""))

    # 3) 集中度：全命中在单一家族 → 1.0
    ranked = [("MAPK1", 3.0), ("MAPK3", 2.0), ("TP53", 1.0)]
    c = concentration(ranked, {"MAPK1", "MAPK3"})
    checks.append(("concentration_max", c["concentration"] == 1.0, str(c["concentration"])))

    # 4) 逐折评估：可分排序 AUC 高
    import tempfile
    d = tempfile.mkdtemp(prefix="np_split_")
    gpath = os.path.join(d, "gold.csv")
    open(gpath, "w", encoding="utf-8").write("gene,tier\nMAPK1,core\nAKT1,core\nIL6,core\n")
    epath = os.path.join(d, "edges.csv")
    open(epath, "w", encoding="utf-8").write(
        "compound,gene\nc1,MAPK1\nc2,MAPK1\nc1,AKT1\nc2,AKT1\nc1,IL6\nc1,X1\nc1,X2\n")
    r = run(edges=epath, gold_csv=gpath, split="complex", k=3, outdir=d)
    checks.append(("run_ok", r.get("status") == "OK", r.get("status")))
    checks.append(("auc_reported", r.get("auc_mean") is not None and r["auc_mean"] > 0.5,
                   str(r.get("auc_mean"))))
    checks.append(("proxy_declared", r.get("proxy") in ("gene_family", "provided_complex_map"),
                   str(r.get("proxy"))))

    # 5) 子集关系检测
    p = os.path.join(d, "db_edges.csv")
    open(p, "w", encoding="utf-8").write(
        "compound,gene,database\n"
        "c1,MAPK1,TCMSP\nc1,AKT1,TCMSP\nc1,IL6,TCMSP\n"
        "c1,MAPK1,ETCM2\nc1,AKT1,ETCM2\n")
    lr = lodo(p, outdir=d)
    checks.append(("lodo_ok", lr.get("status") == "OK", lr.get("status")))
    checks.append(("subset_detected", bool(lr.get("subset_warning")), lr.get("subset_warning", "")[:40]))

    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for n, g_, d_ in checks:
        print("  [%s] %-24s %s" % ("✓" if g_ else "✗", n, d_))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = argparse.ArgumentParser(description="np_split: 基准切分与外部工具对照")
    sub = ap.add_subparsers(dest="fn")

    r = sub.add_parser("run")
    r.add_argument("--edges", default=None)
    r.add_argument("--pred", default=None)
    r.add_argument("--gold-csv", required=True)
    r.add_argument("--split", default="complex", choices=["random", "complex"])
    r.add_argument("--k", type=int, default=5)
    r.add_argument("--seed", type=int, default=42)
    r.add_argument("--topk", default="20")
    r.add_argument("--complex-map", default=None, help="gene,complex 两列；缺省用基因家族代理")
    r.add_argument("--label", default=None)
    r.add_argument("-o", "--outdir", default="np_out/benchmark_split")

    l = sub.add_parser("lodo")
    l.add_argument("--db-edges", required=True)
    l.add_argument("--db-col", default="database")
    l.add_argument("-o", "--outdir", default="np_out/benchmark_split")

    c = sub.add_parser("compare")
    c.add_argument("--gold-csv", required=True)
    c.add_argument("--pred", action="append", default=[], help="name=path，可重复")
    c.add_argument("--topk", default="20")
    c.add_argument("-o", "--outdir", default="np_out/benchmark_split")

    a = ap.parse_args(argv)
    if a.fn == "run":
        res = run(a.edges, a.pred, a.gold_csv, a.split, a.k, a.seed,
                  tuple(int(x) for x in a.topk.split(",")), a.complex_map, a.outdir, a.label)
        if res.get("status") != "OK":
            print("[split] %s %s" % (res.get("status"), res.get("reason", "")))
            return 2
        print("[split] %s(k=%d, proxy=%s)：AUC %.3f±%.3f（%.3f–%.3f），命中集中度 %.2f" % (
            res["split"], res["k"], res["proxy"], res["auc_mean"], res["auc_sd"],
            res["auc_min"], res["auc_max"], res["concentration"]["concentration"]))
        return 0
    if a.fn == "lodo":
        res = lodo(a.db_edges, a.db_col, a.outdir)
        if res.get("status") != "OK":
            print("[lodo] %s" % json.dumps(res, ensure_ascii=False)[:300])
            return 2
        for r in res["results"]:
            print("[lodo] 留出 %-8s AUC=%-6s 覆盖=%s" % (r["held_out_db"], r.get("AUC"), r.get("coverage")))
        if res["subset_warning"]:
            print("  ⚠️ " + res["subset_warning"])
        return 0
    if a.fn == "compare":
        preds = {}
        for x in a.pred:
            if "=" in x:
                k, v = x.split("=", 1)
                preds[k] = v
        res = compare(a.gold_csv, preds, tuple(int(x) for x in a.topk.split(",")), a.outdir)
        for t in res["tools"]:
            print("[compare] %-18s AUC=%-6s hits@20=%-3s EF@20=%s" % (
                t.get("tool"), t.get("AUC"), (t.get("topk_hits") or {}).get("20", ""),
                (t.get("ef") or {}).get("20", "")))
        return 0
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
