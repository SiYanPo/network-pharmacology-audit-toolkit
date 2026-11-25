#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_goldtest.py —— 外部金标准盲测：给 np_fetch / np_enrich 算 Recall 与 FPR（防御"自己测自己"）

为什么需要（SCI 审稿人视角）
    benchmark_23.json 是"自己复现自己"，无法证明方法对未知数据也准。
    本模块引入**独立于本 skill 的外部金标准**：
      * CTD (Comparative Toxicogenomics Database) 的 disease–gene 关联（真实文献策展）
      * DrugBank 已上市 PD 药物靶点（药品-靶点是监管事实，非预测）
    用它们当 Ground Truth 盲测"疾病靶点获取"与"富集通路"的召回与假阳性。

盲测设计（严格防止数据窥探）
    1. 金标准集合 G 从外部文件读入（CTD TSV / DrugBank CSV），skill 内部方法不得参与构建 G。
    2. 预测集合 P = np_fetch 产出的 disease_targets.tsv 的基因集。
    3. Recall = |P ∩ G| / |G|          —— 金标准里有多少被我们的方法覆盖。
    4. FPR 的负集 N：从"与 PD 无关的其他疾病特异基因"构建（如 CTD 中仅关联
       非神经疾病的基因），FPR = |P ∩ N| / |N| —— 预测集里有多大比例是无关噪音。
       注意：NP 的 disease_targets 是"相关基因集"而非二分类器，FPR 定义为
       "预测集中落在已知无关基因集的比例"，是集合层面的假阳性代理指标，论文需如此声明。
    5. 输出 goldtest_report.md 含混淆矩阵、点估计 + 95% Wilson 区间、盲测声明。

用法
    python np_goldtest.py --gold-ctd ctd_pd_genes.tsv --pred disease_targets.tsv \
                          --decoy ctd_decoy_genes.tsv -o out/
    python np_goldtest.py --gold-drugbank drugbank_pd_targets.csv --pred ... -o out/
    python np_goldtest.py --selftest
"""
import argparse
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def read_genes(path, cols=("gene", "Gene Symbol", "symbol", "gene_symbol", "Gene", "靶点")):
    """从 TSV/CSV/txt 读基因集，自动识别列。"""
    out = set()
    if not path or not os.path.exists(path):
        return out
    if path.lower().endswith((".txt", ".list")):
        for line in open(path, encoding="utf-8-sig"):
            g = line.strip()
            if g:
                out.add(g.upper())
        return out
    with open(path, encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f, delimiter="\t" if path.lower().endswith(".tsv") else ",")
        hdrs = rdr.fieldnames or []
        col = None
        for c in cols:
            if c in hdrs:
                col = c
                break
        if col is None:
            for c in hdrs:
                if any(k in c.lower() for k in ("gene", "symbol", "靶点")):
                    col = c
                    break
        if col is None and hdrs:
            col = hdrs[0]
        for r in rdr:
            g = (r.get(col) or "").strip().upper()
            if g and g not in ("GENE", "SYMBOL", "NA"):
                out.add(g)
    return out


def wilson(k, n, z=1.96):
    """比例的 95% Wilson 置信区间。"""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    m = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (c - m) / d), min(1.0, (c + m) / d))


def metrics(pred, gold, decoy):
    """集合层面 Recall / FPR + Wilson 区间。"""
    tp = len(pred & gold)
    fn = len(gold - pred)
    fp = len(pred & decoy)
    tn = max(0, len(decoy) - fp)
    recall = tp / max(1, tp + fn)
    fpr = fp / max(1, fp + tn)
    return {
        "gold_n": len(gold), "pred_n": len(pred), "decoy_n": len(decoy),
        "tp": tp, "fn": fn, "fp": fp, "tn": tn,
        "recall": round(recall, 4), "recall_ci95": [round(x, 4) for x in wilson(tp, tp + fn)],
        "fpr": round(fpr, 4), "fpr_ci95": [round(x, 4) for x in wilson(fp, fp + tn)],
        "overlap_genes": sorted(pred & gold)[:30],
        "false_positive_genes": sorted(pred & decoy)[:30],
    }


def enrichment_goldtest(pred_terms, gold_terms):
    """富集通路盲测：预测富集通路 vs 金标准通路（如 CTD/KEGG PD 通路集）的召回。"""
    P = {str(x).strip().lower() for x in pred_terms if str(x).strip()}
    G = {str(x).strip().lower() for x in gold_terms if str(x).strip()}
    tp = len(P & G)
    return {
        "pred_terms": len(P), "gold_terms": len(G), "matched": tp,
        "recall": round(tp / max(1, len(G)), 4),
        "matched_names": sorted(P & G)[:30],
        "missed": sorted(G - P)[:20],
    }


def report(m, eg, args, out):
    os.makedirs(out, exist_ok=True)
    lines = [
        "# 外部金标准盲测报告", "",
        "## 盲测声明", "",
        "金标准来自独立于本 skill 的外部策展数据（CTD / DrugBank），构建过程未使用本 skill 的任何",
        "方法、阈值或中间产物。预测集来自 `np_fetch` 产出的 disease_targets。两者仅在评测阶段相遇。",
        "", "**指标定义（集合层面，非二分类）**：",
        "- Recall = |预测∩金标准| / |金标准|：我们的方法覆盖了已知 PD 基因的比例。",
        "- FPR = |预测∩负集| / |负集|：预测集落入『已知与其他疾病特异相关、与 PD 无关』基因集的比例，",
        "  是集合层面假阳性代理指标；论文方法学须照此声明，不可表述为『分类器 FPR』。",
        "", "## 混淆矩阵与点估计", "",
        "| 指标 | 值 | 95% CI (Wilson) |", "|---|---|---|",
        "| Recall | %.3f | [%.3f, %.3f] |" % (m["recall"], *m["recall_ci95"]),
        "| FPR | %.3f | [%.3f, %.3f] |" % (m["fpr"], *m["fpr_ci95"]),
        "", "| | 金标准+ | 金标准− |", "|---|---|---|",
        "| 预测+ | TP=%d | FP=%d |" % (m["tp"], m["fp"]),
        "| 预测− | FN=%d | TN=%d |" % (m["fn"], m["tn"]),
        "", "规模：金标准 %d，预测集 %d，负集 %d。" % (m["gold_n"], m["pred_n"], m["decoy_n"]),
        "", "### 命中金标准的基因（前 30）", "",
        ", ".join(m["overlap_genes"]) or "无",
        "", "### 落入负集的假阳性基因（前 30）", "",
        ", ".join(m["false_positive_genes"]) or "无",
    ]
    if eg:
        lines += ["", "## 富集通路盲测", "",
                  "预测富集通路 %d，金标准通路 %d，命中 %d，Recall=%.3f。" %
                  (eg["pred_terms"], eg["gold_terms"], eg["matched"], eg["recall"]),
                  "", "命中：" + ", ".join(eg["matched_names"]), "",
                  "未命中：" + ", ".join(eg["missed"])]
    lines += ["", "## 复核建议", "",
              "- Recall 偏低：检查 np_fetch 的疾病靶点阈值（如 OpenTargets score）是否过严，或金标准版本更新。",
              "- FPR 偏高：检查负集构建是否混入神经相关基因；高 FPR 时不得声称方法特异性。",
              "- 每次上游数据库版本变更后须重跑本盲测，并把本报告纳入 run_manifest。"]
    rp = os.path.join(out, "goldtest_report.md")
    open(rp, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    jp = os.path.join(out, "goldtest_result.json")
    json.dump({"target_metrics": m, "enrichment_metrics": eg,
               "inputs": {k: v for k, v in vars(args).items() if k not in ("fn",)}},
              open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=str)
    print("[goldtest] Recall=%.3f %s | FPR=%.3f %s" %
          (m["recall"], m["recall_ci95"], m["fpr"], m["fpr_ci95"]))
    print("[goldtest] → %s" % rp)
    return m


def selftest():
    pred = {"SNCA", "LRRK2", "PARK7", "TP53", "TNF", "RANDOM1", "RANDOM2"}
    gold = {"SNCA", "LRRK2", "PARK7", "GBA1"}
    decoy = {"RANDOM1", "RANDOM2", "HBB", "ALB"}
    m = metrics(pred, gold, decoy)
    ok = (m["tp"] == 3 and m["fp"] == 2 and abs(m["recall"] - 0.75) < 1e-6)
    print("selftest: %s (tp=%d fp=%d recall=%.3f fpr=%.3f)" %
          ("PASS" if ok else "FAIL", m["tp"], m["fp"], m["recall"], m["fpr"]))
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="np_goldtest: 外部金标准盲测 Recall/FPR")
    ap.add_argument("--gold-ctd", default=None, help="CTD PD 疾病-基因 TSV")
    ap.add_argument("--gold-drugbank", default=None, help="DrugBank PD 药物靶点 CSV")
    ap.add_argument("--pred", default=None, help="本 skill 预测的 disease_targets.tsv")
    ap.add_argument("--decoy", default=None, help="负集：其他疾病特异基因 TSV/txt")
    ap.add_argument("--pred-terms", default=None, help="预测富集通路文件（每行/每列一个通路名）")
    ap.add_argument("--gold-terms", default=None, help="金标准通路文件")
    ap.add_argument("-o", "--outdir", default="np_out/goldtest")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if not a.pred or not (a.gold_ctd or a.gold_drugbank):
        ap.error("需 --pred 且至少一种金标准（--gold-ctd / --gold-drugbank）")
    pred = read_genes(a.pred)
    gold = read_genes(a.gold_ctd or a.gold_drugbank)
    decoy = read_genes(a.decoy) if a.decoy else set()
    if not decoy:
        print("[WARN] 未提供 --decoy，FPR 无法计算（将按 0 处理且不具特异性意义）", file=sys.stderr)
    m = metrics(pred, gold, decoy)
    eg = None
    if a.pred_terms and a.gold_terms:
        pt = [l.strip() for l in open(a.pred_terms, encoding="utf-8-sig") if l.strip()]
        gt = [l.strip() for l in open(a.gold_terms, encoding="utf-8-sig") if l.strip()]
        eg = enrichment_goldtest(pt, gt)
    report(m, eg, a, a.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
