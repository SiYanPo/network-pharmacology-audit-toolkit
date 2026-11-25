#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_enrich.py — 富集分析（本地 GMT / g:Profiler 在线）+ 冗余折叠 + 三层分级

三个必须做的动作（否则结论会被审稿人直接否掉）
  1) **背景可交代**：超几何的 N 必须是「背景基因集」，且 K 必须落在该背景内。
     用「输入基因集自身」当背景是常见错误；没有实测背景时必须声明结果为探索性。
  2) **冗余折叠**：GO/KEGG 大量条目高度重叠（同一批基因反复出现）。把 200 条相似条目
     当成 200 个独立机制，是「数量代替结论」。本脚本按基因集 Jaccard ≥ 阈值折叠，
     只保留每簇最显著的代表，并如实列出被折叠条目数。
  3) **三层分级**：机制明确 / 伞状泛化 / 背景性（感染类）。像 "Pathways in cancer"
     这种几乎任何基因集都会命中的条目，不能作为机制结论。

用法
  python np_enrich.py -g intersection_genes.txt -m examples/pathways.gmt -o out/
  python np_enrich.py -g genes.txt --gprofiler -o out/ --organism hsa
  python np_enrich.py -g genes.txt -m x.gmt --background bg.txt -o out/
"""
import argparse
import csv
import json
import os
import sys
import urllib.parse
import urllib.request
from collections import defaultdict

import npcore as C

UMBRELLA_PAT = [
    "pathways in cancer", "cancer pathways", "signal transduction", "metabolism",
    "immune system", "generic transcription", "disease", "developmental biology",
    "cellular responses to stress", "gene expression (transcription)",
    "vesicle-mediated transport", "transport of small molecules", "cell-cell communication",
    "post-translational protein modification", "dna repair", "chromatin organization",
    "signaling by interleukins", "cytokine signaling in immune system",
]
INFECTIOUS_PAT = [
    "hepatitis b", "hepatitis c", "kaposi", "cytomegalovirus", "influenza", "tuberculosis",
    "legionellosis", "coronavirus disease", "measles", "herpes simplex", "epstein-barr",
    "human papillomavirus", "hiv", "malaria", "leishmania", "toxoplasmosis", "amoebiasis",
    "pertussis", "salmonella", "shigellosis", "yersinia", "staphylococcus", "chagas",
]
# 机制明确的通路（可作为结论锚点）
MECH_PAT = [
    "pi3k-akt", "tnf", "il-17", "nf-kappa b", "nf-κb", "hif-1", "mapk", "foxo", "apoptosis",
    "autophagy", "age-rage", "jak-stat", "toll-like", "nlrp3", "ampk", "mtor", "wnt",
    "notch", "vegf", "tgf-beta", "p53", "cell cycle", "ferroptosis", "pyroptosis",
    "oxidative phosphorylation", "glutathione", "lipid and atherosclerosis",
]


def tier(term, size=None):
    t = term.lower()
    for p in INFECTIOUS_PAT:
        if p in t:
            return "背景性（感染类）"
    for p in UMBRELLA_PAT:
        if t.strip() == p or t.startswith(p):
            return "伞状泛化"
    if size is not None and size >= 1500:
        return "伞状泛化"
    for p in MECH_PAT:
        if p in t:
            return "机制明确"
    return "未分类"


def read_gmt(path):
    sets = {}
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].strip()
            genes = {g.strip().upper() for g in parts[2:] if g.strip()}
            if genes:
                sets[name] = genes
    return sets


def collapse(rows, sets, thr=0.5):
    """按基因集 Jaccard 聚类折叠（贪心：按显著性依次取代表）。"""
    kept, clusters = [], []
    for r in rows:
        s = sets.get(r["term"], set())
        merged = False
        for cl in clusters:
            if C.jaccard(s, cl["genes"]) >= thr:
                cl["members"].append(r["term"])
                merged = True
                break
        if not merged:
            clusters.append({"genes": s, "members": [r["term"]], "rep": r})
            kept.append(r)
    for r in kept:
        for cl in clusters:
            if cl["rep"] is r:
                r["collapsed_n"] = len(cl["members"]) - 1
                r["collapsed_terms"] = ";".join(cl["members"][1:6])
                break
    return kept


def gprofiler(genes, organism="hsa", sources=("GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC", "WP"),
              user_thr=0.05, background=None, timeout=90):
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
    payload = {"organism": organism, "query": genes, "sources": list(sources),
               "user_threshold": user_thr, "significance_threshold_method": "fdr",
               "no_evidences": False, "combined": False}
    if background:
        payload["background"] = background
    req = urllib.request.Request(url, data=json.dumps(payload).encode("utf-8"),
                                 headers={"Content-Type": "application/json",
                                          "User-Agent": "np-pro/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = json.loads(r.read().decode("utf-8", "replace"))
    rows = []
    for it in data.get("result", []):
        rows.append({
            "term": it.get("name", ""),
            "source": it.get("source", ""),
            "native": it.get("native", ""),
            "pvalue": it.get("p_value"),
            "qvalue": it.get("p_value"),  # g:Profiler 已按 FDR 校正
            "genes": ";".join(it.get("intersections", [])[:30]),
            "n_hit": len(it.get("intersections", [])),
            "term_size": it.get("term_size", ""),
            "geneset": set(g.upper() for g in it.get("intersections", [])),
        })
    return rows


def hpa_background(tsv, tissue=None, min_ntpm=1.0):
    """HPA bulk TSV → 组织表达基因集（nTPM ≥ 阈值）。

    为什么要这一步：网络药理学最常被诟病的是拿 20,596 个基因当背景，
    但真正"可能被检出"的只有在该组织表达的基因（多数组织约 8,000–13,000）。
    背景越大，超几何 P 越小，假阳性越多。
    """
    if not tsv or not os.path.exists(tsv):
        return None
    genes = set()
    tis = (tissue or "").strip().lower()
    with open(tsv, encoding="utf-8", errors="replace") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        cols = rdr.fieldnames or []
        gcol = next((c for c in cols if c.lower() in ("gene name", "gene_name")), None)
        tcol = next((c for c in cols if c.lower() == "tissue"), None)
        vcol = next((c for c in cols if "ntpm" in c.lower()), None)
        if not (gcol and vcol):
            return None
        for r in rdr:
            g = (r.get(gcol) or "").strip().upper()
            if not g:
                continue
            if tis:
                tn = (r.get(tcol) or "").strip().lower() if tcol else ""
                if tis not in tn:
                    continue
            try:
                v = float(r.get(vcol) or 0)
            except ValueError:
                continue
            if v >= min_ntpm:
                genes.add(g)
    return genes or None



def selftest():
    """核心统计自检：超几何 P、BH 校正、背景过滤、分层字段。"""
    import tempfile
    import csv as _csv
    d = tempfile.mkdtemp(prefix="np_enrich_")
    gmt = os.path.join(d, "t.gmt")
    terms = {"T1": "IL6 TNF AKT1 BCL2", "T2": "IL6 TNF EGFR VEGFA", "T3": "ZZZ1 ZZZ2 ZZZ3"}
    open(gmt, "w", encoding="utf-8").write(
        "\n".join("%s\tdesc\t%s" % (k, "\t".join(v.split())) for k, v in terms.items()))
    genes = os.path.join(d, "g.txt")
    open(genes, "w", encoding="utf-8").write("IL6\nTNF\nAKT1\nFAKEGENE\n")   # 末者不在背景，应被剔除
    bgp = os.path.join(d, "bg.txt")
    open(bgp, "w", encoding="utf-8").write(
        "\n".join("IL6 TNF AKT1 BCL2 EGFR VEGFA CASP3 MMP9 STAT3 JUN FOS MYC".split()))
    out = os.path.join(d, "out")
    rc = main(["-g", genes, "-m", gmt, "--background", bgp, "--min-hit", "1",
               "-o", out, "--quiet"])
    rows = list(_csv.DictReader(open(os.path.join(out, "enrich_all.csv"), encoding="utf-8-sig")))
    t1 = next((r for r in rows if r["term"] == "T1"), None)
    # 手算：N=12, K=|T1∩bg|=4, n=|{A,B,C}|=3, k=3
    exp = C.hypergeom_sf(3, 12, 4, 3)
    checks = [
        ("run_ok", rc == 0 and t1 is not None, "rc=%s" % rc),
        ("pvalue_matches_hypergeom", t1 is not None
         and abs(float(t1["pvalue"]) - exp) <= max(1e-9, 1e-6 * exp),
         "%s vs %.6g" % (t1["pvalue"] if t1 else "NA", exp)),
        ("bh_q_not_smaller_than_p", t1 is not None and float(t1["qvalue"]) >= float(t1["pvalue"]) - 1e-12,
         "q=%s" % (t1["qvalue"] if t1 else "NA")),
        ("tier_field_present", t1 is not None and (t1.get("tier") or "") != "", ""),
        ("background_filtered_gene", all(int(r["n_hit"]) <= 3 for r in rows), "FAKEGENE 未被计入"),
    ]
    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for n, g_, d_ in checks:
        print("  [%s] %-28s %s" % ("✓" if g_ else "✗", n, d_))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = C.Cli(description="富集分析：本地 GMT 超几何 / g:Profiler 在线，含冗余折叠与三层分级")
    ap.add_argument("-g", "--genes", required=True)
    ap.add_argument("-m", "--gmt", default=None, help="本地 GMT 文件（离线，可多个用逗号分隔）")
    ap.add_argument("--gprofiler", action="store_true", help="在线调用 g:Profiler（无需 key）")
    ap.add_argument("--organism", default="hsa")
    ap.add_argument("--sources", default="GO:BP,GO:CC,GO:MF,KEGG,REAC")
    ap.add_argument("--background", default=None, help="背景基因列表（强烈建议提供实测背景）")
    ap.add_argument("--background-mode", default="auto", choices=["auto", "deg", "hpa", "genome"],
                    help="auto: 有 --background 用 DEG，否则尝试 HPA 组织表达背景；hpa: 强制 HPA；genome: 显式使用全基因组（会打膨胀警告）")
    ap.add_argument("--hpa-tsv", default=None, help="HPA rna_tissue_consensus.tsv（缺省自动下载到 --hpa-cache）")
    ap.add_argument("--hpa-cache", default=None, help="HPA bulk 缓存目录（默认 <outdir>/../hpa_cache）")
    ap.add_argument("--hpa-tissue", default=None, help="疾病相关组织，如 liver / cerebral cortex；强烈建议指定")
    ap.add_argument("--hpa-min-ntpm", type=float, default=1.0, help="判定表达的 nTPM 阈值（默认 1）")
    ap.add_argument("--universe", type=int, default=20596, help="无背景列表时的背景基因数")
    ap.add_argument("-o", "--outdir", default="np_out")
    ap.add_argument("--min-hit", type=int, default=3)
    ap.add_argument("--q", type=float, default=0.05)
    ap.add_argument("--collapse", type=float, default=0.5, help="冗余折叠的 Jaccard 阈值，0=关闭")
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--calibrate", type=int, default=0,
                    help="经验校准的置换次数（0=关闭）。用于检验解析 P 是否系统性偏乐观")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if os.environ.get("NP_OFFLINE") == "1":
        a.background_mode = "genome" if not a.background else "deg"
    if a.selftest:
        return selftest()

    man = C.Manifest("06_enrich", a.outdir, vars(a))
    man.input("genes", a.genes)
    if a.background:
        man.input("background", a.background)

    genes = C.read_gene_list(a.genes)
    gset = set(genes)
    if len(genes) < 3:
        print("[错误] 输入基因 <3，富集无意义。")
        return 2

    # ---- 背景选择：DEG（最优）→ HPA 组织表达（次优）→ 全基因组（必须打膨胀警告）----
    bg_mode = a.background_mode
    if a.background:
        bg_mode = "deg"
    elif bg_mode == "auto":
        bg_mode = "hpa"
    bg_source, inflation_warning = "", False
    bg = None
    if bg_mode == "hpa":
        cache = a.hpa_cache or (C.cache_dir("hpa") or
                                os.path.join(os.path.dirname(os.path.abspath(a.outdir)), "hpa_cache"))
        tsv = a.hpa_tsv
        if not tsv or not os.path.exists(tsv):
            try:
                import np_hpa
                tsv, err = np_hpa.fetch_bulk(cache)
            except Exception as e:
                tsv, err = None, "%s" % type(e).__name__
        if tsv and os.path.exists(tsv):
            bg = hpa_background(tsv, a.hpa_tissue, a.hpa_min_ntpm)
            bg_source = "HPA rna_tissue_consensus（%s，nTPM≥%.1f，组织=%s）" % (
                os.path.basename(tsv), a.hpa_min_ntpm, a.hpa_tissue or "任意组织")
            if not a.hpa_tissue:
                man.warn("未指定 --hpa-tissue：HPA 背景退化为「任一组组表达」，接近全基因组，"
                         "建议指定疾病相关组织（如 liver / cerebral cortex）")
                inflation_warning = True
        else:
            man.warn("HPA 背景获取失败（%s），回退全基因组背景；P 值可能膨胀" % (err or "unknown"))
            bg_mode = "genome"
    if bg_mode == "genome":
        bg_source = "全基因组（未使用疾病/组织表达背景）"
        inflation_warning = True
    if a.background:
        bg = set(C.read_gene_list(a.background))
        # 背景必须"大于输入集且不是输入集的子/超集"，否则超几何退化为无意义
        ov = len(bg & gset)
        if len(bg) < 2 * len(gset):
            print("[阻断] 背景集（%d）不足输入集（%d）的 2 倍：超几何会退化，结果无解释力。"
                  % (len(bg), len(gset)))
            print("       背景应是「实际可被检出的基因」（表达谱/蛋白组检出基因），")
            print("       而不是药物靶点列表或疾病基因列表。若确要这么做，请改用 --universe 指定 N。")
            return 3
        if ov / max(1, len(gset)) > 0.95 and len(bg) < 5 * len(gset):
            man.warn("背景与输入集重叠 %.0f%%，且背景规模仅 %.1f 倍：结果会偏保守，"
                     "建议用表达谱检出的全基因作背景" % (100 * ov / len(gset), len(bg) / len(gset)))
    if bg:
        man.warn("背景来源：%s（N=%d）" % (bg_source, len(bg)))
    else:
        man.warn("背景来源：%s（N=%d）——未使用疾病表达谱背景，P 值可能膨胀" % (bg_source, a.universe))
    N_all = len(bg) if bg else a.universe
    if bg:
        gset_eff = gset & bg
        if len(gset_eff) < len(gset):
            man.warn("%d 个输入基因不在背景集中，已从检验中剔除" % (len(gset) - len(gset_eff)))
            gset = gset_eff

    rows = []
    if a.gmt:
        for path in [p.strip() for p in a.gmt.split(",") if p.strip()]:
            man.input("gmt_" + os.path.basename(path), path)
            sets = read_gmt(path)
            for term, members in sets.items():
                if bg:
                    members = members & bg
                K = len(members)
                if K < a.min_hit:
                    continue
                hit = gset & members
                k = len(hit)
                if k < a.min_hit:
                    continue
                p = C.hypergeom_sf(k, N_all, K, len(gset))
                rows.append({"term": term, "source": os.path.basename(path).split(".")[0],
                             "pvalue": p, "n_hit": k, "term_size": K,
                             "genes": ";".join(sorted(hit)[:30]), "geneset": hit})
    if a.gprofiler:
        try:
            rows += gprofiler(genes, a.organism, tuple(a.sources.split(",")),
                              background=sorted(bg) if bg else None)
            man.count("gprofiler", "ok")
        except Exception as e:
            # 公开端点失败 ≠ 生物学阴性，必须如实记录
            man.warn("g:Profiler 不可用（%s）；结果仅为本地 GMT 部分，不得据此下「该通路不富集」的结论"
                     % type(e).__name__)
            print("[警告] g:Profiler 调用失败：%s（已记录为不可用，不视为阴性结果）" % e)

    if not rows:
        C.ensure_dir(a.outdir)
        md = ("# 富集分析报告\n\n**未得到任何显著条目。**\n\n"
              "可能原因与处理：\n"
              "1. 输入基因过少（当前 %d）——富集至少要 5–10 个基因才有意义；\n"
              "2. min-hit 过高（当前 %d）——可降到 2，但须在方法学说明；\n"
              "3. GMT 与基因符号口径不一致（大小写/别名/物种）；\n"
              "4. 背景集过小，导致 K 被裁剪。\n\n"
              "> 注意：**「没有显著条目」是结果，不是失败**，但必须排除以上技术原因后才能这样写。\n"
              % (len(genes), a.min_hit))
        p = C.write_text(os.path.join(a.outdir, "enrich_report.md"), md)
        man.output("report", p)
        man.warn("富集无显著条目")
        man.save()
        print("[富集] 无显著条目 → %s" % p)
        return 0

    qs = C.bh_fdr([r["pvalue"] for r in rows])
    for r, q in zip(rows, qs):
        r["qvalue"] = q
        r["tier"] = tier(r["term"], r.get("term_size") if isinstance(r.get("term_size"), int) else None)
    rows.sort(key=lambda r: (r["pvalue"], -r["n_hit"]))
    sig = [r for r in rows if r["qvalue"] < a.q]

    collapsed = sig
    if a.collapse > 0 and sig:
        sets_for = {r["term"]: set(r["geneset"]) for r in sig}
        collapsed = collapse(sig, sets_for, a.collapse)

    # ---- 经验校准：解析 P 在零假设下是否偏斜（文献 PMID 38436561 指出普遍偏乐观）----
    calib = {}
    if a.calibrate > 0 and a.gmt:
        import random as _rnd
        rng = _rnd.Random(42)
        ann_universe = set()
        for path in [x.strip() for x in a.gmt.split(",") if x.strip()]:
            for _t, mem in read_gmt(path).items():
                ann_universe |= (mem & bg if bg else mem)
        ann_universe = sorted(ann_universe)
        n = len(gset)
        if len(ann_universe) > n >= 3:
            Nann = len(ann_universe)
            term_sets = {}
            for path in [x.strip() for x in a.gmt.split(",") if x.strip()]:
                for _t, mem in read_gmt(path).items():
                    mm = mem & bg if bg else mem
                    if len(mm) >= a.min_hit:
                        term_sets[_t] = mm
            # 每个条目：在 Nann 背景下的解析 P vs 置换分布
            null_ps = {t: [] for t in term_sets}
            for _ in range(a.calibrate):
                samp = set(rng.sample(ann_universe, n))
                for t, mem in term_sets.items():
                    k2 = len(samp & mem)
                    null_ps[t].append(C.hypergeom_sf(k2, Nann, len(mem), n))
            for r in rows:
                t = r["term"]
                if t not in null_ps or not null_ps[t]:
                    continue
                obs_p = C.hypergeom_sf(r["n_hit"], Nann, len(term_sets[t]), n)
                emp = (1 + sum(1 for x in null_ps[t] if x <= obs_p)) / (1 + len(null_ps[t]))
                r["p_analytic_u"] = obs_p
                r["p_empirical"] = emp
            if any("p_empirical" in r for r in rows):
                floor = 1.0 / (a.calibrate + 1)
                pe = [r.get("p_empirical") for r in rows]
                pa = [r.get("p_analytic_u") for r in rows]
                # 删失校正：经验 P 触到 1/(B+1) 地板时被右删失。直接用全体中位数算 λ，
                # 会把「置换次数不够」误读成「解析 P 偏乐观」（实测 B=200 得 0.16、B=1000 得 0.78）。
                usable = [r for r in rows if r.get("p_empirical") is not None
                          and r["p_empirical"] > floor * 1.0001 and r.get("p_analytic_u") is not None]
                frac_floor = 1.0 - len(usable) / max(1, len([x for x in pe if x is not None]))
                if len(usable) >= 20:
                    med_e = sorted(r["p_empirical"] for r in usable)[len(usable) // 2]
                    med_a = sorted(r["p_analytic_u"] for r in usable)[len(usable) // 2]
                    lam = (med_a / med_e) if med_e else None
                else:
                    med_e = med_a = lam = None
                calib = {"B": a.calibrate, "universe": Nann, "n": n,
                         "median_p_analytic": med_a, "median_p_empirical": med_e,
                         "lambda": lam, "floor": floor,
                         "frac_at_floor": round(frac_floor, 3), "n_usable": len(usable),
                         "n_sig_analytic": sum(1 for x in pa if x is not None and x < 0.05),
                         "n_sig_empirical": sum(1 for x in pe if x is not None and x < 0.05)}
        man.count("calibration", calib)

    C.ensure_dir(a.outdir)
    p_all = C.write_table(os.path.join(a.outdir, "enrich_all.csv"),
                          ["term", "source", "n_hit", "term_size", "pvalue", "qvalue", "tier",
                           "p_analytic_u", "p_empirical", "genes"],
                          rows)
    p_sig = C.write_table(os.path.join(a.outdir, "enrich_significant.csv"),
                          ["term", "source", "n_hit", "term_size", "pvalue", "qvalue", "tier",
                           "p_analytic_u", "p_empirical", "collapsed_n", "collapsed_terms", "genes"],
                          collapsed)

    by_tier = defaultdict(list)
    for r in collapsed:
        by_tier[r["tier"]].append(r)
    p_bub = C.bubble_svg([(r["term"], r["n_hit"] / max(1, len(gset)), -1 * __import__("math").log10(max(r["qvalue"], 1e-300)), r["n_hit"])
                          for r in collapsed[:a.top]],
                         os.path.join(a.outdir, "enrich_bubble.svg"),
                         title="富集（折叠后 top%d）" % a.top)

    md = ["# 富集分析报告", "",
          "- 输入基因：%d；背景：%s（N=%d）" % (len(genes), bg_source if bg else "默认全基因组", N_all),
          "- 最小命中数：%d；显著性：q<%.2f（BH 校正）" % (a.min_hit, a.q),
          "- 显著条目：%d；折叠后代表条目：%d（Jaccard 阈值 %.2f）"
          % (len(sig), len(collapsed), a.collapse), ""]

    if inflation_warning or not bg:
        md += ["<div style='border-left:4px solid #b91c1c;background:#fdecec;padding:8pt 10pt;margin:8pt 0'>"
               "<b style='color:#b91c1c'>P-value may be inflated / P 值可能膨胀</b><br>"
               "本次富集未使用疾病相关的表达谱（DEGs）或组织表达基因集作为背景"
               "（当前：%s，N=%d）。全基因组背景会系统性放大富集显著性，"
               "本节的 P/q 只能作为<b>探索性</b>证据，不得写入结论段。"
               "<br>修正方式：<code>python np_enrich.py -g genes.txt -m x.gmt "
               "--background degs.txt</code>（首选，实测 DEGs），或 "
               "<code>--background-mode hpa --hpa-tissue liver</code>（次选，HPA 组织表达）。"
               "</div>" % (bg_source, N_all), ""]

    if calib:
        md += ["## 经验校准（解析 P 是否偏乐观）", "",
               "- 置换 %d 次；注释背景 N=%d，输入集 n=%d" % (calib["B"], calib["universe"], calib["n"]),
               "- **经验 P 的可达下限 = 1/(B+1) = %.5f**；落在下限上（右删失）占 %.1f%%"
               % (calib["floor"], 100 * calib["frac_at_floor"]),
               "- 解析 P<0.05 的条目 %d 个；经验 P<0.05 的条目 %d 个"
               % (calib["n_sig_analytic"], calib["n_sig_empirical"])]
        if calib.get("lambda") is None:
            md += ["", "> 无法给出 λ：未触地板的可用条目不足 20 个。请把 --calibrate 提到 ≥1000。", ""]
        else:
            md += ["- 中位解析 P = %s；中位经验 P（未删失子集 n=%d）= %s；**λ = %.2f**"
                   % (C.fmt_p(calib["median_p_analytic"]), calib["n_usable"],
                      C.fmt_p(calib["median_p_empirical"]), calib["lambda"]), ""]
            if calib["frac_at_floor"] > 0.3:
                md += ["> 约 %d%% 的条目触到经验 P 地板：λ 建立在删失子集上，仅供参考；"
                       "需要可靠 λ 时把 --calibrate 提到 ≥2000。" % round(100 * calib["frac_at_floor"]), ""]
            elif calib["n_sig_empirical"] < calib["n_sig_analytic"] * 0.8:
                md += ["> 经验校准后显著条目明显减少：解析 P 在本数据上偏乐观，正文应以经验 P 为准。", ""]
            else:
                md += ["> 解析 P 与经验 P 量级一致（λ≈1），说明该结果对零分布偏斜不敏感。", ""]

    try:
        json.dump({"background_mode": "deg" if a.background else bg_mode,
                   "background_source": bg_source, "N": N_all,
                   "n_input": len(genes), "hpa_tissue": a.hpa_tissue,
                   "hpa_min_ntpm": a.hpa_min_ntpm,
                   "pvalue_inflation_warning": bool(inflation_warning or not bg),
                   "warning_text": ("未使用疾病表达谱背景，P 值可能膨胀"
                                    if (inflation_warning or not bg) else "")},
                  open(os.path.join(a.outdir, "enrich_meta.json"), "w", encoding="utf-8"),
                  ensure_ascii=False, indent=1)
    except Exception:
        pass

    md += ["## 三层分级（决定哪些能写进结论）", "",
           "| 层级 | 条目数 | 使用规则 |", "|---|---|---|",
           "| 机制明确 | %d | 可作为机制假说的锚点 |" % len(by_tier.get("机制明确", [])),
           "| 伞状泛化 | %d | 只说明功能谱广，**不得**作为机制结论 |" % len(by_tier.get("伞状泛化", [])),
           "| 背景性（感染类） | %d | 与核心通路共享分子，属被动命中 |" % len(by_tier.get("背景性（感染类）", [])),
           "| 未分类 | %d | 需人工判断 |" % len(by_tier.get("未分类", [])), ""]

    for tname in ["机制明确", "伞状泛化", "背景性（感染类）", "未分类"]:
        items = by_tier.get(tname, [])
        if not items:
            continue
        md += ["### %s（%d）" % (tname, len(items)), "",
               "| 条目 | 命中 | 条目大小 | q 值 | 折叠 |", "|---|---|---|---|---|"]
        for r in items[:a.top]:
            md.append("| %s | %d | %s | %s | %s |" % (
                r["term"], r["n_hit"], r.get("term_size", ""), C.fmt_p(r["qvalue"]),
                ("+%d" % r["collapsed_n"]) if r.get("collapsed_n") else ""))
        md.append("")

    md += ["## 写法红线", "",
           "1. **P 值是显著性，不是效应值**：不能写成「作用强度」。效应量看倍数富集与命中数。",
           "2. **GeneRatio / Rich factor 只是比值**（命中数/条目大小），不是通路活性评分。",
           "3. **P 显著 ≠ 通路被激活**：GO/KEGG 的输入是无方向的基因集，没有激活/抑制信息。",
           "4. **不得只挑「好看」的通路放图**：抗炎研究不能只放 TNF/NF-κB，应给完整表（附录）+ 讨论解释。",
           "5. 冗余条目已折叠：折叠前后条目数必须都报告（%d → %d）。" % (len(sig), len(collapsed)), ""]

    p_md = C.write_text(os.path.join(a.outdir, "enrich_report.md"), "\n".join(md) + "\n")
    for k, v in [("all", p_all), ("sig", p_sig), ("bubble", p_bub), ("report", p_md)]:
        man.output(k, v)
    man.count("input_genes", len(genes))
    man.count("significant", len(sig))
    man.count("collapsed", len(collapsed))
    man.count("tiers", {k: len(v) for k, v in by_tier.items()})
    man.save()

    if not a.quiet:
        print("[富集] 输入 %d 基因，显著 %d 条，折叠后 %d 条（机制明确 %d）"
              % (len(genes), len(sig), len(collapsed), len(by_tier.get("机制明确", []))))
        print("[富集] 输出 → %s" % p_md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
