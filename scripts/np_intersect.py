#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_intersect.py — 靶点交集与"成分级"富集（C-T-D 的核心一步）

比"取交集画韦恩图"多做三件事
  1) 每个交集靶点带**证据强度**：疾病侧打分、来源数、来源类型（遗传证据 / 表达证据 / 文本挖掘）
  2) 每个成分做**精确超几何检验**：它的靶点撞上疾病基因，是巧合还是确实偏多
     （回答"谁的靶点最集中在疾病上"，而不是只看靶点数多少）
  3) 零交集/低交集时的**恢复序列**（这是最常见的卡点，不能只说"结果不显著"）

输入
  --drug       成分-靶点表（compound,gene 或 compound,target_name）+ --disease 疾病靶点表
  或 --sets    A.txt,B.txt[,...] 任意多组基因集
输出
  intersection.csv / intersection_genes.txt / venn.svg / compound_enrichment.csv / intersect_report.md
"""
import argparse
import os
import sys
from collections import defaultdict

import npcore as C

SOURCE_TIER = {
    "omim": ("遗传证据", 3), "clinvar": ("遗传证据", 3), "gwas": ("遗传证据", 3),
    "disgenet": (" curated 关联", 2), "ctd": (" curated 关联", 2), "ttd": (" curated 关联", 2),
    "drugbank": (" curated 关联", 2), "open targets": (" curated 关联", 2),
    "genecards": ("整合/文本挖掘", 1), "malacards": ("整合/文本挖掘", 1),
    "geo": ("表达证据", 2), "deg": ("表达证据", 2), "wgcna": ("表达证据", 2),
}


def tier_of(src):
    s = (src or "").lower()
    for k, v in SOURCE_TIER.items():
        if k in s:
            return v
    return ("未标注来源", 1)



def selftest():
    """交集逻辑自检：集合运算正确、非交集基因不混入、证据分层字段非空。"""
    import tempfile
    import csv as _csv
    d = tempfile.mkdtemp(prefix="np_intersect_")
    drug = os.path.join(d, "drug.csv")
    open(drug, "w", encoding="utf-8").write(
        "compound,gene\nc1,A\nc1,B\nc2,B\nc2,C\n")
    dis = os.path.join(d, "dis.tsv")
    open(dis, "w", encoding="utf-8").write(
        "gene\tscore\tsource\nB\t5\tGeneCards\nC\t1\tCTD\nZ\t9\tOMIM\n")
    out = os.path.join(d, "out")
    rc = main(["--drug", drug, "--disease", dis, "-o", out, "--quiet"])
    rows = list(_csv.DictReader(open(os.path.join(out, "intersection.csv"), encoding="utf-8-sig")))
    genes = {r["gene"] for r in rows}
    b = next((r for r in rows if r["gene"] == "B"), {})
    checks = [
        ("run_ok", rc == 0 and bool(rows), "rc=%s n=%d" % (rc, len(rows))),
        ("intersection_exact", genes == {"B", "C"}, ",".join(sorted(genes))),
        ("excludes_non_shared", "A" not in genes and "Z" not in genes, ""),
        ("multi_compound_support", str(b.get("n_compounds")) == "2", str(b.get("n_compounds"))),
        ("evidence_tier_present", all((r.get("evidence_tier") or "") for r in rows), ""),
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
    ap = C.Cli(description="靶点交集 + 成分级超几何富集 + 证据强度分层")
    ap.add_argument("--drug", help="成分-靶点表（long format）")
    ap.add_argument("--drug-col", default=None, help="成分列名（默认自动识别）")
    ap.add_argument("--gene-col", default=None, help="基因/靶点列名（默认自动识别）")
    ap.add_argument("--disease", help="疾病靶点表（gene,score,source）")
    ap.add_argument("--sets", help="逗号分隔的基因列表文件，做 N 组交集（替代 --drug/--disease）")
    ap.add_argument("-o", "--outdir", default="np_out")
    ap.add_argument("--min-score", type=float, default=None, help="疾病靶点最低打分")
    ap.add_argument("--universe", default="genome",
                    help="富集背景：genome=人类蛋白编码基因 20596（默认）；tested=全部成分靶点并集"
                         "（即可被检出的空间，更保守）；或直接给整数")
    ap.add_argument("--label", default="药物-疾病", help="报告标题用的标签")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)

    man = C.Manifest("03_intersect", a.outdir, vars(a))
    sets, comp_targets = {}, {}
    disease_meta = {}

    if a.sets:
        for path in [p.strip() for p in a.sets.split(",") if p.strip()]:
            name = os.path.splitext(os.path.basename(path))[0]
            sets[name] = set(C.read_gene_list(path))
            man.input("set_" + name, path)
    else:
        if not a.drug or not a.disease:
            print("[错误] 需要 --drug + --disease，或用 --sets 指定多个基因列表。")
            return 2
        h, rows = C.read_table(a.drug)
        cd = a.drug_col or C.resolve_col(h, "compound", required=True)
        # 优先用已标准化的 symbol 列
        cg = a.gene_col or ("symbol" if "symbol" in h else None) or C.resolve_col(h, "gene", required=True)
        for r in rows:
            comp = (r.get(cd) or "").strip()
            gene = (r.get(cg) or "").strip().upper()
            if comp and gene and gene not in ("—", "-", "NA"):
                comp_targets.setdefault(comp, set()).add(gene)
        sets["药物靶点"] = set().union(*comp_targets.values()) if comp_targets else set()
        man.input("drug", a.drug)

        h2, rows2 = C.read_table(a.disease)
        cg2 = C.resolve_col(h2, "gene", required=True)
        cs2 = C.resolve_col(h2, "score")
        csrc = C.resolve_col(h2, "source")
        for r in rows2:
            g = (r.get(cg2) or "").strip().upper()
            if not g:
                continue
            sc = C.to_float(r.get(cs2), 0.0) if cs2 else 0.0
            if a.min_score is not None and sc < a.min_score:
                continue
            src = (r.get(csrc) or "") if csrc else ""
            prev = disease_meta.get(g)
            if prev is None or sc > prev["score"]:
                disease_meta[g] = {"score": sc, "source": src,
                                   "sources": {src} if src else set()}
            else:
                prev["sources"].add(src)
        sets["疾病靶点"] = set(disease_meta)
        man.input("disease", a.disease)

    names = list(sets)
    inter = set.intersection(*sets.values()) if sets else set()

    C.ensure_dir(a.outdir)
    p_venn = C.venn_svg(sets, os.path.join(a.outdir, "venn.svg"),
                        title="%s 靶点交集" % a.label)

    # ---- 交集表 ----
    rows_out = []
    for g in sorted(inter):
        comps = sorted([c for c, tg in comp_targets.items() if g in tg]) if comp_targets else []
        dm = disease_meta.get(g, {})
        tname, tval = tier_of(dm.get("source", ""))
        rows_out.append({
            "gene": g,
            "n_compounds": len(comps),
            "compounds": ";".join(comps[:12]) + ("…" if len(comps) > 12 else ""),
            "disease_score": dm.get("score", ""),
            "disease_sources": ";".join(sorted(dm.get("sources", set()))),
            "evidence_tier": tname,
            "evidence_weight": tval,
        })
    rows_out.sort(key=lambda r: (-(r["n_compounds"]), -(C.to_float(r["disease_score"], 0) or 0)))
    p_inter = C.write_table(os.path.join(a.outdir, "intersection.csv"),
                            ["gene", "n_compounds", "compounds", "disease_score",
                             "disease_sources", "evidence_tier", "evidence_weight"], rows_out)
    p_genes = C.write_text(os.path.join(a.outdir, "intersection_genes.txt"),
                           "\n".join(sorted(inter)) + "\n")

    # ---- 成分级超几何 ----
    p_enr = None
    enr_rows = []
    if comp_targets and inter:
        K = max(1, len(sets.get("疾病靶点", set())))
        # 两套背景：全基因组 / "可检出空间"（所有成分靶点并集）
        # 前者回答"是否高于随机基因"，后者回答"是否高于随机可成药基因"——后者更保守、更贴近机制问题
        universe_all = set().union(*comp_targets.values()) or set()
        dis_set = sets.get("疾病靶点", set())
        if a.universe == "genome":
            N, mode = 20596, "genome"
        elif a.universe == "tested":
            N, mode = (len(universe_all) or 20596), "tested"
        else:
            N, mode = int(a.universe), "custom"
        K = max(1, len(dis_set))
        # 保守背景：只在"可检出空间"内比较——K2 只数该空间里的疾病基因，避免 K>N 的非法情形
        N2 = len(universe_all) or N
        K2 = max(1, len(dis_set & universe_all))
        for comp, tg in comp_targets.items():
            n = len(tg)
            k = len(tg & inter)
            exp = n * K / N
            exp2 = n * K2 / N2
            row = {"compound": comp, "n_targets": n, "n_shared": k,
                   "expected": round(exp, 2),
                   "fold": round(k / max(1e-9, exp), 2) if n else 0,
                   "pvalue": C.hypergeom_sf(k, N, K, n),
                   "expected_tested": round(exp2, 2),
                   "fold_tested": round(k / max(1e-9, exp2), 2) if n else 0,
                   "pvalue_tested": C.hypergeom_sf(k, N2, K2, n)}
            enr_rows.append(row)
        for r, q in zip(enr_rows, C.bh_fdr([r["pvalue"] for r in enr_rows])):
            r["qvalue"] = q
        for r, q in zip(enr_rows, C.bh_fdr([r["pvalue_tested"] for r in enr_rows])):
            r["qvalue_tested"] = q
        n_sig = sum(1 for r in enr_rows if r["qvalue"] < 0.05)
        n_sig2 = sum(1 for r in enr_rows if r["qvalue_tested"] < 0.05)
        if n_sig == len(enr_rows) and len(enr_rows) >= 5:
            man_warn_power = ("成分级富集在「全基因组」背景下全部显著（%d/%d）：该检验在此数据下已无区分度。"
                              "应以「可检出空间」背景（N=%d, K=%d，显著 %d/%d）与倍数富集（fold）为准，"
                              "并在报告中写明背景口径。" % (n_sig, len(enr_rows), N2, K2, n_sig2, len(enr_rows)))
        else:
            man_warn_power = ""
        enr_rows.sort(key=lambda r: (r["n_shared"], r["fold"]), reverse=True)
        cols = ["compound", "n_targets", "n_shared", "expected", "fold", "pvalue", "qvalue",
                "expected_tested", "fold_tested", "pvalue_tested", "qvalue_tested"]
        p_enr = C.write_table(os.path.join(a.outdir, "compound_enrichment.csv"), cols, enr_rows)
        if n_sig == len(enr_rows) and len(enr_rows) >= 5:
            man_warn_power = ("成分级富集「全部显著」（%d/%d）：说明该数据下超几何检验已无区分度，"
                              "结论不能靠 P 值，应改用倍数富集（fold）排序并报告背景口径。" % (n_sig, len(enr_rows)))
        else:
            man_warn_power = ""

    # ---- 报告 ----
    md = ["# 交集与成分富集报告（%s）" % a.label, "",
          "| 集合 | 规模 |", "|---|---|"]
    for n in names:
        md.append("| %s | %d |" % (n, len(sets[n])))
    md += ["| **交集** | **%d** |" % len(inter), ""]

    if not inter:
        md += ["## 零交集：按顺序排查（不要直接放弃）", "",
               "1. **ID 空间**：药侧是蛋白名还是基因符号？未标准化的蛋白名与符号取交集必然为零。"
               "先跑 `np_standardize.py`。",
               "2. **物种**：药侧是否混入大鼠/小鼠/菌源靶点？STRING 与 TCMSP 默认人类。",
               "3. **疾病集过窄**：GeneCards 只取高分区（如 score≥10）时可能只剩几十个基因；"
               "先放宽到 score≥1 或取中位数以上。",
               "4. **成分过少**：OB/DL 把方剂筛到 <5 个成分时，靶点数不足以覆盖疾病集；"
               "改用含量/入血成分表，或补充例外纳入。",
               "5. **疾病名不规范**：用英文标准病名 + 同义词（如 diabetic nephropathy / diabetic kidney disease）"
               "重取一次并取并集。",
               "6. 以上都排除后仍为零：这是一个**有信息量的阴性结果**——说明在所选证据口径下，"
               "该方与该病没有共享靶点，不应通过换阈值硬凑。", ""]
    else:
        if len(inter) < 5:
            md += ["## ⚠ 交集仅 %d 个基因：统计功效不足" % len(inter), "",
                   "富集分析在 <5 个基因时几乎不可能得到稳定结果，且随机对照检验无意义。",
                   "处理：扩充疾病靶点来源（多库并集）、增加药味、放宽成分阈值（并声明），"
                   "或改为「以疾病为中心」的反向网络药理学（`np_optimize.py`）。", ""]

        tiers = defaultdict(int)
        for r in rows_out:
            tiers[r["evidence_tier"]] += 1
        md += ["## 交集靶点的证据构成", "", "| 证据层级 | 靶点数 | 解读 |", "|---|---|---|"]
        TIER_NOTE = {"遗传证据": "最强（OMIM/ClinVar/GWAS），可作为机制锚点",
                     " curated 关联": "较强（人工审编数据库）",
                     "表达证据": "情境依赖（特定样本/组织）",
                     "整合/文本挖掘": "最弱（GeneCards/MalaCards 汇总，非因果）",
                     "未标注来源": "需补来源"}
        for k, v in sorted(tiers.items(), key=lambda x: -x[1]):
            md.append("| %s | %d | %s |" % (k, v, TIER_NOTE.get(k, "")))
        md += ["", "> 写法提醒：GeneCards 的 score 是**整合相关性**，不是因果证据；"
               "GEO 差异表达是「特定样本 vs 对照」的变化，不等于「疾病靶点」。"
               "报告中应按上表分层陈述，不能把三者混为「疾病相关基因」。", ""]

        top = rows_out[:15]
        md += ["## 交集靶点（按共享成分数排序，前 15）", "",
               "| 基因 | 共享成分数 | 疾病打分 | 来源 | 证据层级 |", "|---|---|---|---|---|"]
        for r in top:
            md.append("| %s | %d | %s | %s | %s |" % (
                r["gene"], r["n_compounds"],
                ("%.1f" % r["disease_score"]) if isinstance(r["disease_score"], float) else r["disease_score"],
                r["disease_sources"] or "—", r["evidence_tier"]))
        md.append("")

        if enr_rows:
            sig = [r for r in enr_rows if r["qvalue"] < 0.05]
            md += ["## 成分级富集（谁的靶点最集中在疾病上）", "",
                   "- 检验：精确超几何；背景 N=%s（%s），疾病集 K=%s" % (N, mode, K),
                   "- 显著成分（q<0.05）：%d / %d" % (len(sig), len(enr_rows)), "",
                   "| 成分 | 靶点数 | 共享 | 期望 | 倍数 | q 值 |", "|---|---|---|---|---|---|"]
            for r in enr_rows[:12]:
                md.append("| %s | %d | %d | %.1f | %.2f | %s |" % (
                    r["compound"], r["n_targets"], r["n_shared"], r["expected"], r["fold"],
                    C.fmt_p(r["qvalue"])))
            if man_warn_power:
                md += ["", "> ⚠ " + man_warn_power]
            md += ["", "> 注意：这是「靶点是否富集于疾病基因集」的统计，**不是**「该成分有疗效」。"
                   "靶点数多的成分（常为泛成分）天然更容易显著，解释时须结合含量与实测证据。", ""]

    p_md = C.write_text(os.path.join(a.outdir, "intersect_report.md"), "\n".join(md) + "\n")

    for k, v in [("intersection", p_inter), ("genes", p_genes), ("venn", p_venn),
                 ("enrichment", p_enr), ("report", p_md)]:
        if v:
            man.output(k, v)
    man.count("sets", {n: len(sets[n]) for n in names})
    man.count("intersection", len(inter))
    if "man_warn_power" in dir() and man_warn_power:
        man.warn(man_warn_power)
    if not inter:
        man.warn("交集为零，已输出恢复序列")
    elif len(inter) < 5:
        man.warn("交集 <5，统计功效不足")
    man.save()

    if not a.quiet:
        print("[交集] " + " ∩ ".join("%s(%d)" % (n, len(sets[n])) for n in names)
              + " = %d 个基因" % len(inter))
        if enr_rows:
            nsig = sum(1 for r in enr_rows if r["qvalue"] < 0.05)
            print("[交集] 成分级富集显著：%d/%d（q<0.05）" % (nsig, len(enr_rows)))
        print("[交集] 输出 → %s" % p_inter)
    return 0


if __name__ == "__main__":
    sys.exit(main())
