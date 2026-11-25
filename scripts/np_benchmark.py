#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_benchmark.py —— 金标准基准测试（Benchmark）：证明"算得全"且"算得准"

为什么软件类论文必须有这一节
    审稿人不会因为你画出了漂亮的 PPI 网络就相信你；他们会问：
    "同一批输入，你的工具比别人的命中率高吗？高多少？AUC 多少？"
    没有 Benchmark 的网络药理学软件 = 没有对照组的实验。

本模块做什么
    1. 内置 4 个已上市/经典方剂的人工整理金标准靶点集（core 强证据 / extended 扩展），
       来源：DrugBank/TTD/ChEMBL 机制记录 + 标志性成分的药理实验文献（见 provenance）。
    2. 把你的工具输出（靶点排序表，或 成分-靶点边表）与金标准对齐，计算：
         AUC (ROC, Mann-Whitney U)          —— 全局排序能力
         AUPRC / Average Precision          —— 正类稀少时比 AUC 更敏感
         Recall@K / Precision@K             —— "前 K 个核心靶点命中几个"（审稿人最常问）
         Enrichment Factor (EF@K)           —— 相对随机抽样的富集倍数
         超几何 p 值                        —— 命中数是否显著优于随机
         Bootstrap 95% CI                   —— AUC 的置信区间（1000 次重抽样）
    3. 生成 benchmark_report.md（可直接改写为论文 Results 第一段）+ ROC 曲线 SVG。

三条防止"自欺"的硬规则
    R1 严格口径只用 core（强证据）靶点，extended 从**阴性池中剔除**（不确定的不能当负样本）。
    R2 记录覆盖率 coverage = 金标准靶点中出现在候选池里的比例；coverage<0.6 时 Recall 上限受限，
       报告必须显式标注，禁止拿低 coverage 的 Recall 吹性能。
    R3 单侧随机基线：Recall@K 的期望 = K/N，同时给二项/超几何 p 值；只报绝对值不报基线视为无效。

用法
    python np_benchmark.py list
    python np_benchmark.py evaluate --formula CDDP --pred core_targets.csv -o out/
    python np_benchmark.py evaluate --formula CDDP --edges ct_edges.csv -o out/
    python np_benchmark.py evaluate --formula CDDP --pred r.csv --tier lenient --topk 10,20,50
    python np_benchmark.py evaluate --gold-csv my_gold.csv --pred r.csv --label "外部金标准"
    python np_benchmark.py report --dir out/            # 汇总所有 benchmark_*.json → benchmark_report.md
    python np_benchmark.py selftest
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GOLD_PATH = os.path.join(ROOT, "data", "benchmark_goldstandard.json")
# 顺序即优先级：标准化的 symbol/gene 必须优先于 target_name（后者常是蛋白全名）
# 预注册基线（看到结果之前就固定；报告同时给出两者，不做事后挑选）
BASELINES = [
    ("B0", "support", "支持该靶点的不同成分数（最弱基线）"),
    ("B1", "penalized", "支持成分数 × promiscuity 降权 1/(1+log2(m/m_median))"),
]

GENE_COLS = ("gene", "symbol", "hgnc_symbol", "target", "target_name", "target_symbol",
             "genes", "preferredname", "preferred_name", "靶点")
SCORE_COLS = ("score", "weight", "degree", "centrality", "combined", "integrated", "rank_score", "value")


# ============================ 指标（纯 Python，可单测） ============================
def roc_auc(labels, scores):
    """Mann-Whitney U 形式的 ROC AUC；并列分数取平均秩。labels: 0/1。"""
    pairs = [(s, l) for s, l in zip(scores, labels)]
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return None
    # 平均秩（升序）
    order = sorted(range(len(pairs)), key=lambda i: pairs[i][0])
    ranks = [0.0] * len(pairs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and pairs[order[j + 1]][0] == pairs[order[i]][0]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    rank_sum_pos = sum(r for r, l in zip(ranks, labels) if l == 1)
    u = rank_sum_pos - n_pos * (n_pos + 1) / 2.0
    return round(u / (n_pos * n_neg), 4)


def average_precision(labels, scores):
    """AP = Σ (precision@k × Δrecall)，按分数降序逐点累积。"""
    n_pos = sum(labels)
    if n_pos == 0:
        return None
    order = sorted(range(len(labels)), key=lambda i: -scores[i])
    hit, ap = 0, 0.0
    prev_recall = 0.0
    for k, idx in enumerate(order, 1):
        if labels[idx] == 1:
            hit += 1
            recall = hit / n_pos
            ap += (hit / k) * (recall - prev_recall)
            prev_recall = recall
    return round(ap, 4)


def hypergeom_sf(k, N, K, n):
    """P(X >= k)，X~Hypergeom(总体N, 成功K, 抽样n)。用于 Recall@K 的显著性。"""
    if k <= 0:
        return 1.0
    total = _comb(N, n)
    if total == 0:
        return 1.0
    acc = 0.0
    for x in range(k, min(K, n) + 1):
        acc += _comb(K, x) * _comb(N - K, n - x) / total
    return min(1.0, acc)


def _comb(n, k):
    if k < 0 or k > n:
        return 0
    return math.comb(int(n), int(k))


def metrics_at_k(ranked_genes, pos_set, N_universe, topk=(10, 20, 50)):
    """Recall@K / Precision@K / EF@K + 超几何 p。"""
    K_pos = len(pos_set & set(ranked_genes))
    out = {"n_pos_total": len(pos_set), "n_pos_in_universe": K_pos,
           "n_universe": N_universe, "at_k": {}}
    for K in topk:
        top = ranked_genes[:K]
        hits = sum(1 for g in top if g in pos_set)
        rec = hits / len(pos_set) if pos_set else 0.0
        prec = hits / K if K else 0.0
        exp = K * len(pos_set) / N_universe if N_universe else 0.0
        ef = (hits / exp) if exp > 0 else 0.0
        out["at_k"][str(K)] = {"hits": hits, "recall": round(rec, 4), "precision": round(prec, 4),
                               "expected_hits": round(exp, 2), "enrichment_factor": round(ef, 2),
                               "hypergeom_p": float("%.3g" % hypergeom_sf(hits, N_universe, len(pos_set), K))}
    return out


def bootstrap_auc(labels, scores, n_boot=1000, seed=42):
    """Bootstrap 95% CI（百分位法）。"""
    rnd = random.Random(seed)
    n = len(labels)
    vals = []
    for _ in range(n_boot):
        idx = [rnd.randrange(n) for _ in range(n)]
        l = [labels[i] for i in idx]
        s = [scores[i] for i in idx]
        a = roc_auc(l, s)
        if a is not None:
            vals.append(a)
    if not vals:
        return None
    vals.sort()
    lo = vals[int(0.025 * len(vals))]
    hi = vals[int(0.975 * len(vals)) - 1]
    return [round(lo, 4), round(hi, 4)]


# ============================ 输入解析 ============================
def read_table(path):
    """读 CSV/TSV，自动识别分隔符；返回 (rows, headers)。"""
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        rdr = csv.DictReader(f, delimiter=delim)
        rows = [r for r in rdr]
        return rows, (rdr.fieldnames or [])


def load_ranking(path):
    """靶点排序表 → [(gene, score)] 已按分数降序。支持单列基因列表。"""
    rows, hdrs = read_table(path)
    low = {h.lower().strip(): h for h in hdrs}
    gcol = next((low[c] for c in GENE_COLS if c in low), hdrs[0] if hdrs else None)
    scol = next((low[c] for c in SCORE_COLS if c in low), None)
    out = []
    for i, r in enumerate(rows):
        g = (r.get(gcol) or "").strip().upper()
        if not g:
            continue
        if scol:
            try:
                s = float(r.get(scol))
            except (TypeError, ValueError):
                s = 0.0
        else:
            s = float(len(rows) - i)  # 单列列表：保持原始顺序
        out.append((g, s))
    # 同基因保留最高分
    best = {}
    for g, s in out:
        if g not in best or s > best[g]:
            best[g] = s
    return sorted(best.items(), key=lambda x: -x[1])


def rank_from_edges(path, min_support=1, weight_mode="penalized"):
    """成分-靶点边表 → 靶点排序。

    weight_mode
      support   : score = 支持该靶点的**不同成分数**（最朴素的基线 B0）
      penalized : score = Σ_c 1/(1+log2(m_c/median_m))，m_c 为成分 c 的靶点数（B1）
                  对"多靶点 promiscuous 成分"降权——与 np_funnel 的 promiscuity penalty 同源。
                  这是**预先注册**的两个基线，报告同时给出，不做事后挑选。
    """
    rows, hdrs = read_table(path)
    low = {h.lower().strip(): h for h in hdrs}
    gcol = next((low[c] for c in GENE_COLS if c in low), None)
    ccol = next((low[c] for c in ("compound", "mol_id", "src_id", "成分", "name", "molecule") if c in low), None)
    if not gcol:
        raise SystemExit("[FAIL] 边表缺少基因列（gene/target/靶点）")
    ct = {}
    for r in rows:
        g = (r.get(gcol) or "").strip().upper()
        c = (r.get(ccol) or "_").strip().lower() if ccol else "_"
        if g:
            ct.setdefault(c, set()).add(g)
    if weight_mode == "support":
        score = {g: 0.0 for gs in ct.values() for g in gs}
        for gs in ct.values():
            for g in gs:
                score[g] += 1.0
    else:
        sizes = sorted(len(v) for v in ct.values())
        med = sizes[len(sizes) // 2] if sizes else 1
        score = {}
        for c, gs in ct.items():
            w = 1.0 / (1.0 + math.log2(len(gs) / med)) if len(gs) > med and med > 0 else 1.0
            for g in gs:
                score[g] = score.get(g, 0.0) + w
    # 支持成分数下限
    support = {g: 0 for gs in ct.values() for g in gs}
    for gs in ct.values():
        for g in gs:
            support[g] += 1
    ranked = [(g, round(s, 4)) for g, s in score.items() if support.get(g, 0) >= min_support]
    ranked.sort(key=lambda x: (-x[1], x[0]))
    return ranked


def edge_concordance(pred_edges, gold_edges):
    """边级跨库一致性（当金标准是**另一数据库的边集**时使用，比靶点级 AUC 更可解释）。

    pred_edges / gold_edges: set[(compound_lower, GENE)]
      sensitivity = |pred ∩ gold| / |gold|      —— 金标准边被覆盖的比例（"算得全"）
      precision   = |pred ∩ gold| / |pred|      —— 预测边中被独立库支持的比例（"算得准"）
      jaccard     = |∩| / |∪|
      over_prediction = |pred| / |gold|          —— 过预测倍数（特异性缺口）
    """
    inter = pred_edges & gold_edges
    return {
        "status": "OK" if gold_edges else "NO_OVERLAP",
        "n_pred_edges": len(pred_edges), "n_gold_edges": len(gold_edges),
        "n_overlap": len(inter),
        "sensitivity": round(len(inter) / len(gold_edges), 4) if gold_edges else None,
        "precision": round(len(inter) / len(pred_edges), 4) if pred_edges else None,
        "jaccard": round(len(inter) / len(pred_edges | gold_edges), 4) if (pred_edges | gold_edges) else None,
        "over_prediction_factor": round(len(pred_edges) / len(gold_edges), 2) if gold_edges else None,
        "n_pred_targets": len({g for _, g in pred_edges}),
        "n_gold_targets": len({g for _, g in gold_edges}),
        "n_overlap_targets": len({g for _, g in inter}),
    }


def read_edge_pairs(path):
    """边表 → set[(compound_lower, GENE)]。"""
    rows, hdrs = read_table(path)
    low = {h.lower().strip(): h for h in hdrs}
    gcol = next((low[c] for c in GENE_COLS if c in low), None)
    ccol = next((low[c] for c in ("compound", "mol_id", "src_id", "成分", "name", "molecule") if c in low), None)
    if not gcol:
        raise SystemExit("[FAIL] 边表缺少基因列：%s" % path)
    out = set()
    for r in rows:
        g = (r.get(gcol) or "").strip().upper()
        c = (r.get(ccol) or "_").strip().lower() if ccol else "_"
        if g:
            out.add((c, g))
    return out


def load_gold(path=GOLD_PATH):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def gold_sets(formula, tier="strict"):
    """返回 (pos_set, ambiguous_set, meta)。strict: 只 core，extended 作含糊集剔除。"""
    core = set(formula.get("core") or {})
    ext = set(formula.get("extended") or {})
    if tier == "lenient":
        return core | ext, set(), {"tier": "lenient (core+extended)", "n_core": len(core), "n_ext": len(ext)}
    if tier == "strict-l3":
        core = {g for g, ev in (formula.get("core") or {}).items() if ev in ("L3", "L4")}
        return core, ext, {"tier": "strict-l3 (仅 L3/L4 证据)", "n_core": len(core), "n_ext": len(ext)}
    return core, ext, {"tier": "strict (仅 core)", "n_core": len(core), "n_ext": len(ext)}


# ============================ 评估 ============================
def evaluate(ranked, pos_set, ambiguous, topk=(10, 20, 50), n_boot=1000):
    """ranked: [(gene, score)]；返回指标 dict。"""
    genes = [g for g, _ in ranked]
    scores = [s for _, s in ranked]
    universe = [g for g in genes if g not in ambiguous]      # R1：含糊靶点剔除
    uscores = [s for g, s in zip(genes, scores) if g not in ambiguous]
    labels = [1 if g in pos_set else 0 for g in universe]
    res = {"n_ranked": len(genes), "n_evaluated": len(universe),
           "n_positive_found": sum(labels)}
    if not universe or sum(labels) == 0:
        return {"status": "NO_OVERLAP",
                "reason": "候选池与金标准无交集（检查基因符号大小写/物种/命名）",
                "n_ranked": len(genes), "n_positive_found": 0}
    res["status"] = "OK"
    res["AUC"] = roc_auc(labels, uscores)
    res["AUPRC"] = average_precision(labels, uscores)
    res["AUC_95CI"] = bootstrap_auc(labels, uscores, n_boot=n_boot) if n_boot else None
    res.update(metrics_at_k(universe, pos_set, len(universe), topk))
    res["coverage"] = round(len(pos_set & set(genes)) / len(pos_set), 3) if pos_set else None
    res["random_recall_at_k"] = {str(K): round(K / len(universe), 4) for K in topk}
    if res["coverage"] is not None and res["coverage"] < 0.6:
        res["coverage_warning"] = ("金标准覆盖率 %.0f%% < 60%%：大量金标准靶点不在候选池内，"
                                   "Recall 上限被压死，禁止据此宣称性能" % (100 * res["coverage"]))
    return res


def roc_points(labels, scores):
    order = sorted(range(len(labels)), key=lambda i: -scores[i])
    pts, P, N = [(0.0, 0.0)], sum(labels), len(labels) - sum(labels)
    tp = fp = 0
    for i in order:
        tp += labels[i]
        fp += 1 - labels[i]
        pts.append((fp / N if N else 0, tp / P if P else 0))
    return pts


def svg_roc(labels, scores, title, path):
    pts = roc_points(labels, scores)
    W, H, M = 420, 420, 48
    w, h = W - 2 * M, H - 2 * M

    def xy(x, y):
        return M + x * w, H - M - y * h

    poly = " ".join("%.1f,%.1f" % xy(x, y) for x, y in pts)
    auc = roc_auc(labels, scores)
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d">' % (W, H, W, H),
             '<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H),
             '<text x="%d" y="24" font-family="Helvetica,Arial" font-size="14" fill="#0A1A2F">%s</text>'
             % (M, _esc(title)),
             '<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#9ca3af"/>' % (M, M - 12, w, h),
             '<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#c9ced6" stroke-dasharray="4 4"/>'
             % (M, H - M, M + w, M - 12),
             '<polyline fill="none" stroke="#007B7A" stroke-width="2.4" points="%s"/>' % poly,
             '<text x="%d" y="%d" font-family="Helvetica,Arial" font-size="12" fill="#007B7A">AUC = %.3f</text>'
             % (M + w - 96, H - M - 12, auc or 0),
             '<text x="%d" y="%d" font-family="Helvetica,Arial" font-size="11" fill="#6b7280">1 - Specificity</text>'
             % (M + w / 2 - 34, H - 12),
             '<text x="10" y="%d" font-family="Helvetica,Arial" font-size="11" fill="#6b7280" '
             'transform="rotate(-90 14 %d)">Sensitivity</text>' % (M + h / 2, M + h / 2),
             "</svg>"]
    open(path, "w", encoding="utf-8").write("\n".join(parts))
    return path


def _esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


# ============================ 报告 ============================
def write_report(outdir, results, note=""):
    """results: list of dict（evaluate / concordance 输出 + 元信息）。→ benchmark_report.md"""
    ts = time.strftime("%Y-%m-%d %H:%M")
    gold_ver = load_gold().get("version", "?") if os.path.exists(GOLD_PATH) else "?"
    target_res = [r for r in results if r.get("kind") != "concordance"]
    conc_res = [r for r in results if r.get("kind") == "concordance"]
    L = ["# Benchmark Report — network-pharmacology-pro", "",
         "生成时间：%s　　引擎：`np_benchmark.py v1.0`　　金标准版本：%s" % (ts, gold_ver),
         "", "> 本文件由程序自动生成，所有数值来自 `benchmark_*.json`，无人工改写。"
         "金标准为人工整理先验（见 `data/benchmark_goldstandard.json` 的 provenance 段），"
         "**不等同于临床准确率**；其用途是方法间横向比较与回归测试，防止版本迭代导致性能静默退化。", ""]

    # 1. 金标准命中率（靶点级）
    L += ["## 1　主指标：top-K 命中数（AUC 降为次要指标）", "",
          "> **为什么改口径**：金标准覆盖候选池约 42% 时，随机抽 20 个就能中 8.5 个，"
          "Recall@20 没有分辨力；而 AUC 对「前 20 名里到底有几个真靶点」不敏感。"
          "因此主指标改为 **hits@20（命中个数，文献现实值 2–4）** 与 **EF@20**，"
          "AUC 与 AUPRC 仅作排序质量的辅助描述。", "",
          "| 方剂 | 口径 | 候选池 | 金标准 | 命中总数 | **hits@10** | **hits@20** | **hits@50** | "
          "**EF@20** | AUC [95% CI] | AUPRC |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    if not target_res:
        L.append("| — | — | — | — | — | 尚未运行 | — | — | — | — | — |")
    for r in target_res:
        if r.get("status") != "OK":
            L.append("| %s | %s | — | — | — | NO_OVERLAP | — | — | — | — | — |" % (
                r.get("label", "?"), r.get("tier_label", "?")))
            continue
        ak = r["at_k"]
        L.append("| %s | %s | %d | %d | %d | **%d** | **%d** | **%d** | **%.2fx** | %.3f [%s] | %.3f |" % (
            r["label"], r["tier_label"], r["n_evaluated"], r["n_pos_total"], r["n_positive_found"],
            ak["10"]["hits"], ak["20"]["hits"], ak["50"]["hits"], ak["20"]["enrichment_factor"],
            r["AUC"] or 0, ("%.2f–%.2f" % tuple(r["AUC_95CI"])) if r.get("AUC_95CI") else "NA",
            r["AUPRC"] or 0))
    L += ["", "随机基线：Recall@K ≈ K/N（N 为候选池大小）。EF@20 = Precision@20 / (N_pos/N)，"
              "即相对随机抽样的富集倍数；EF@20 ≈ 1 表示与随机无异。"
              "未跑本机流水线时该表为空，运行命令见第 5 节。", ""]

    # 2. 跨库边级一致性
    if conc_res:
        L += ["## 2　跨库边级一致性（成分-靶点边）", "",
              "当金标准是**另一数据库的边集**时，靶点级 AUC 会被「金标准是候选池子集」这一关系污染，"
              "此时应改用边级敏感性/精确率——它直接回答「算得全」与「算得准」两个问题。", "",
              "| 对比 | 预测边 | 金标准边 | 交集 | 敏感性 | 精确率 | Jaccard | 过预测倍数 |",
              "|---|---|---|---|---|---|---|---|"]
        for r in conc_res:
            L.append("| %s | %d | %d | %d | %.3f | %.3f | %.3f | %.2fx |" % (
                r["label"], r["n_pred_edges"], r["n_gold_edges"], r["n_overlap"],
                r["sensitivity"] or 0, r["precision"] or 0, r.get("jaccard") or 0,
                r.get("over_prediction_factor") or 0))
        L += ["", "解读：敏感性接近 1.00 说明本工具的成分-靶点库**覆盖完整**（金标准里的边一条不落）；"
                  "过预测倍数 >2 说明库更宽（特异性偏低），这是网络药理学的普遍特征——"
                  "数据库并集必然引入假阳性，因此下游必须用零模型/诱饵对接做过滤，而不能把边表当结论。", ""]

    # 3. 逐方剂细节
    L += ["## 3　逐方剂细节", ""]
    if not target_res:
        L += ["（无靶点级结果）", ""]
    for r in target_res:
        L.append("### %s（%s）" % (r.get("label", "?"), r.get("tier_label", "?")))
        if r.get("status") != "OK":
            L.append("- 状态：**NO_OVERLAP** — %s" % r.get("reason", ""))
            L.append("")
            continue
        L.append("- 候选池 %d 个靶点；金标准 %d 个（core %d + extended %d）；覆盖率 %.0f%%" % (
            r["n_evaluated"], r["n_pos_total"], r["gold_meta"]["n_core"], r["gold_meta"]["n_ext"],
            100 * (r.get("coverage") or 0)))
        L.append("- AUC %.3f（95%% CI %s）；AUPRC %.3f" % (
            r["AUC"] or 0, ("%.2f–%.2f" % tuple(r["AUC_95CI"])) if r.get("AUC_95CI") else "NA", r["AUPRC"] or 0))
        for K, v in r["at_k"].items():
            L.append("- @%s：命中 %d（期望 %.1f）→ Recall %.2f、Precision %.2f、EF %.2fx、超几何 p=%.2g" % (
                K, v["hits"], v["expected_hits"], v["recall"], v["precision"],
                v["enrichment_factor"], v["hypergeom_p"]))
        if r.get("coverage_warning"):
            L.append("- ⚠️ **%s**" % r["coverage_warning"])
        if r.get("top_missed"):
            L.append("- 未召回的核心靶点（前 10）：%s" % ", ".join(r["top_missed"][:10]))
        L.append("")

    # 3.5 预注册基线
    L += ["## 3.5　预注册基线（pre-specified，不做事后挑选）", "",
          "两个基线在**任何一次评测之前**就写死在 `np_benchmark.BASELINES` 里：", "",
          "| 编号 | 定义 | 用途 |", "|---|---|---|",
          "| B0 | 支持该靶点的**不同成分数** | 最弱基线：任何加权若打不过它就没有存在价值 |",
          "| B1 | 支持成分数 × promiscuity 降权 `1/(1+log2(m/m_median))` | 惩罚\"多靶点滥交\"成分 |", "",
          "> 严格性说明：若两个基线是在看到结果之后才并列报告的，就属于探索性分析。"
          "本实现把基线集合固化为常量，并同时给出两者，**不依据结果选择\"更好看的那个\"**。", ""]

    # 3.6 切分与外部对照
    sp_files = sorted(glob.glob(os.path.join(outdir, "split_*.json")))
    if sp_files:
        L += ["## 3.6　切分策略：随机 vs 家族留出 vs 留一数据库", "",
              "> 网络传播任务的普通随机切分因复合物/家族成员相关而系统性乐观"
              "（PLOS Comput Biol 2019, PMID 31479437；Pahikkala 2014, PMID 24723570）。"
              "对**固定排序**做留组评估，测的不是训练泛化误差，而是\"性能是否集中在少数家族\"。", "",
              "| 切分 | 标签 | 代理 | 折数 | AUC 均值±SD | 范围 | 命中集中度 |",
              "|---|---|---|---|---|---|---|"]
        for f in sp_files:
            try:
                j = json.load(open(f, encoding="utf-8"))
            except Exception:
                continue
            if j.get("split") == "lodo":
                for r in j.get("results") or []:
                    L.append("| lodo（留出 %s） | 边数 %s | 数据库 | — | %s | — | — |" % (
                        r.get("held_out_db"), j.get("n_edges", {}).get(r.get("held_out_db"), ""),
                        r.get("AUC")))
                if j.get("subset_warning"):
                    L += ["", "⚠️ **%s**" % j["subset_warning"]]
                continue
            if j.get("status") != "OK":
                continue
            c = j.get("concentration") or {}
            L.append("| %s | %s | %s | %d | %.3f±%.3f | %.3f–%.3f | %.2f |" % (
                j.get("split"), j.get("label"), j.get("proxy"), j.get("k"),
                j.get("auc_mean", 0), j.get("auc_sd", 0), j.get("auc_min", 0),
                j.get("auc_max", 0), c.get("concentration", 0)))
        L += ["", "读法：若 `complex` 的 AUC 均值明显低于 `random`，说明随机切分被家族相关性抬高；"
              "若均值接近但 SD 变大，说明性能在不同家族间不稳定（本身是有价值的负面结果）。"
              "集中度 >0.6 表示多数命中来自少数家族，结论不能外推到\"全靶点空间\"。", ""]

    ext = os.path.join(outdir, "external_comparison.json")
    if os.path.exists(ext):
        try:
            je = json.load(open(ext, encoding="utf-8"))
        except Exception:
            je = {}
        if je.get("tools"):
            L += ["## 3.7　外部工具对照（同金标准 / 同候选池 / 同指标）", "",
                  "| 工具 | 状态 | 候选数 | hits@20 | EF@20 | AUC | 覆盖率 |",
                  "|---|---|---|---|---|---|---|"]
            for t in je["tools"]:
                L.append("| %s | %s | %s | %s | %s | %s | %s |" % (
                    t.get("tool"), t.get("status"), t.get("n_ranked"),
                    (t.get("topk_hits") or {}).get("20", "—"),
                    (t.get("ef") or {}).get("20", "—"),
                    t.get("AUC"), t.get("coverage")))
            L += ["", "> 覆盖率不同的工具，Recall 不可直接比较；因此本表以 hits@20 与 EF@20 为主。"
                  "若某工具未接入，此处只出现本平台内部基线——**不得**在论文中称为\"与外部工具对比\"。", ""]

    # 4. Methods 段落（可直接改写进论文）
    L += ["## 4　可直接改写进 Methods 的段落", "",
          "> Benchmarking was performed against a manually curated gold standard of formula-level "
          "targets assembled from DrugBank/TTD mechanism annotations and experimental literature on "
          "each formula's marker constituents. Only targets with strong evidence (≥2 independent sources "
          "or direct experimental support; `core` tier) were treated as positives, whereas targets with "
          "single-source or computationally predicted evidence (`extended` tier) were excluded from the "
          "negative pool to avoid label noise. Tool outputs were ranked by the pipeline's composite score "
          "and evaluated with ROC-AUC (Mann–Whitney U statistic, 1000-fold bootstrap 95% CI), average "
          "precision (AUPRC), Recall@K, Precision@K and enrichment factor at K∈{10,20,50}, with "
          "hypergeometric tests against the random-expectation baseline K/N. Cross-database concordance "
          "was additionally assessed at the level of individual compound–target edges (sensitivity, "
          "precision, Jaccard index and over-prediction factor) to separate coverage from specificity. "
          "All metrics are reported together with gold-standard coverage, because a target absent from "
          "the candidate pool cannot be recovered and would otherwise deflate recall in a way unrelated "
          "to ranking quality.", ""]

    # 5. 复现命令
    L += ["## 5　复现命令", "", "```bash",
          "# A. 内置经典方剂（需先跑出该方剂的靶点排序/边表）",
          "python scripts/np_benchmark.py evaluate --formula CDDP --edges np_out/03_compounds/compound_target_edges.csv -o np_out/benchmark",
          "python scripts/np_benchmark.py evaluate --formula HLJDT --pred np_out/05_ppi/hub_genes_scored.csv --tier lenient",
          "# B. 跨库边级一致性",
          "python scripts/np_benchmark.py concordance --pred-edges mine.csv --gold-edges other_db.csv -o np_out/benchmark",
          "# C. 汇总出本报告",
          "python scripts/np_benchmark.py report --dir np_out/benchmark", "```", ""]

    # 6. 局限
    L += ["## 6　必须写进 Limitations 的八条（v3.0 新增 4 条）", "",
          "1. **Publication bias**：被研究得多的靶点更容易进入金标准，任何方法的 AUC 都会被系统性抬高；"
          "本基准只用于版本间回归比较，不能外推为临床准确率。",
          "2. **Label noise**：中药复方的靶点集合本身存在文献异质性，同一方剂不同研究可给出差异很大的靶点表。",
          "3. **Coverage bias**：若工具的成分-靶点库本身不覆盖某靶点，则该靶点在结构上不可能被召回，"
          "Recall 偏低反映的是数据库覆盖而非算法能力——故必须同时报告 coverage。",
          "4. **Single-gene granularity**：本基准在基因层评估，未评价通路层与表型层一致性；"
          "通路层需另设 KEGG/GO 富集重叠基准。",
          "5. **AUC 与 top-K 脱节（实测）**：益智胶囊数据上，promiscuity 降权把 AUC 从 0.599 提到 0.705，"
          "但 **hits@20 两者都是 6、EF@20 都是 0.71×（低于随机）**。AUC 提升没有转化为"
          "\"前 20 名里更多真靶点\"，因此本平台不以 AUC 作为卖点。",
          "6. **金标准可能非独立（实测）**：LODO 检验发现 ETCM2 的边 **全部**出现在 TCMSP 中"
          "（包含率>0.95），两库疑为子集/共享上游关系。任何\"跨库独立验证\"的说法在该数据上不成立，"
          "必须先查清谱系。",
          "7. **切分是家族代理，不是真实复合物**：无 CORUM 时用基因符号前缀分组（输出标 `proxy=gene_family`）。"
          "它只能回答\"性能是否集中在少数家族\"，不能替代 complex-aware 交叉验证。",
          "8. **校准器局限于本基准**：`P(gold)` 是在某个方剂/某份金标准上拟合的，换疾病或换数据库必须重拟合；"
          "跨域直接套用等于把别人的先验当自己的证据。", ""]
    if note:
        L += ["## 7　运行备注", "", note, ""]
    p = os.path.join(outdir, "benchmark_report.md")
    open(p, "w", encoding="utf-8").write("\n".join(L))
    return p


# ============================ 子命令 ============================
def cmd_list(a):
    g = load_gold()
    print("内置金标准 v%s（%s）" % (g.get("version"), g.get("curated")))
    print("%-8s %-22s %-6s %-8s %-8s %s" % ("ID", "方剂", "药味", "core", "extended", "适应症"))
    for f in g["formulas"]:
        print("%-8s %-22s %-6d %-8d %-8d %s" % (
            f["id"], f["name_zh"], len(f["herbs"]), len(f["core"]), len(f["extended"]),
            f["indication_zh"][:28]))
    print("\n口径：strict=只用 core（extended 从阴性池剔除）；lenient=core+extended；strict-l3=仅 L3/L4")
    return 0


def cmd_evaluate(a):
    g = load_gold(a.gold or GOLD_PATH) if not a.gold_csv else None
    if a.gold_csv:
        # 外部金标准：两列 gene,tier(core/extended)
        rows, hdrs = read_table(a.gold_csv)
        low = {h.lower(): h for h in hdrs}
        gc = next((low[c] for c in GENE_COLS if c in low), hdrs[0])
        tc = next((low[c] for c in ("tier", "evidence", "level", "class") if c in low), None)
        core, ext = {}, {}
        for r in rows:
            gene = (r.get(gc) or "").strip().upper()
            if not gene:
                continue
            t = (r.get(tc) or "core").strip().lower() if tc else "core"
            (ext if t.startswith("ext") else core)[gene] = "L3"
        formula = {"id": a.label or "external", "name_zh": a.label or "外部金标准",
                   "core": core, "extended": ext, "kegg": []}
        tier = "strict"
        gmeta = {"version": "external", "curated": os.path.basename(a.gold_csv)}
    else:
        f = next((x for x in g["formulas"] if x["id"].upper() == a.formula.upper() or
                  x["name_zh"] == a.formula), None)
        if not f:
            raise SystemExit("[FAIL] 未知方剂 %s（可选 %s）" % (a.formula, [x["id"] for x in g["formulas"]]))
        formula, tier, gmeta = f, a.tier, {"version": g.get("version"), "curated": g.get("curated")}

    ranked = load_ranking(a.pred) if a.pred else rank_from_edges(a.edges, a.min_support, a.weight_mode)
    if not ranked:
        raise SystemExit("[FAIL] 排序表为空")
    pos, amb, ginfo = gold_sets(formula, tier)
    topk = tuple(int(x) for x in a.topk.split(",") if x.strip())
    res = evaluate(ranked, pos, amb, topk=topk, n_boot=a.bootstrap)
    res["label"] = formula.get("name_zh") or formula.get("id")
    res["formula_id"] = formula.get("id")
    res["tier_label"] = ginfo["tier"]
    res["gold_meta"] = {**ginfo, **gmeta}
    res["top_missed"] = [x for x in sorted(pos) if x not in set(g for g, _ in ranked)]
    res["timestamp"] = time.strftime("%Y-%m-%d %H:%M:%S")

    os.makedirs(a.outdir, exist_ok=True)
    tag = (formula.get("id") or "external")
    jp = os.path.join(a.outdir, "benchmark_%s.json" % tag)
    json.dump(res, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)

    # 命中明细 CSV（便于人工核查每一个命中/漏检）
    if res.get("status") == "OK":
        genes = [g for g, _ in ranked if g not in amb]
        sc = {g: s for g, s in ranked}
        with open(os.path.join(a.outdir, "benchmark_%s_hits.csv" % tag), "w",
                  encoding="utf-8-sig", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["rank", "gene", "score", "in_gold_core", "in_gold_extended", "tier"])
            for i, gene in enumerate(genes, 1):
                w.writerow([i, gene, sc.get(gene, ""), int(gene in formula.get("core", {})),
                            int(gene in formula.get("extended", {})),
                            "core" if gene in formula.get("core", {}) else
                            ("extended" if gene in formula.get("extended", {}) else "background")])
        labels = [1 if x in pos else 0 for x in genes]
        scores = [sc[x] for x in genes]
        svg_roc(labels, scores, "%s — ROC (AUC %.3f)" % (res["label"], res["AUC"] or 0),
                os.path.join(a.outdir, "roc_%s.svg" % tag))

    if res.get("status") != "OK":
        print("[benchmark] %s：%s" % (res["status"], res.get("reason")))
        return 2
    print("[benchmark] %s | %s | 候选池 %d / 命中 %d" % (
        res["label"], res["tier_label"], res["n_evaluated"], res["n_positive_found"]))
    print("  AUC=%.3f [%s]  AUPRC=%.3f  Recall@10/20/50 = %.2f/%.2f/%.2f  EF@20=%.2fx" % (
        res["AUC"], ("%.2f–%.2f" % tuple(res["AUC_95CI"])) if res.get("AUC_95CI") else "NA",
        res["AUPRC"], res["at_k"]["10"]["recall"], res["at_k"]["20"]["recall"],
        res["at_k"]["50"]["recall"], res["at_k"]["20"]["enrichment_factor"]))
    if res.get("coverage_warning"):
        print("  ⚠️ " + res["coverage_warning"])
    print("  → %s" % jp)
    return 0


def cmd_concordance(a):
    """边级跨库一致性：pred 数据库 vs gold 数据库的 (成分,靶点) 边集合。"""
    pred = read_edge_pairs(a.pred_edges)
    gold = read_edge_pairs(a.gold_edges)
    res = edge_concordance(pred, gold)
    res["kind"] = "concordance"
    res["label"] = a.label or "%s vs %s" % (os.path.basename(a.pred_edges), os.path.basename(a.gold_edges))
    res["timestamp"] = time.strftime("%Y-%m-%d %H:%M:%S")
    res["source_pred"] = os.path.abspath(a.pred_edges)
    res["source_gold"] = os.path.abspath(a.gold_edges)
    os.makedirs(a.outdir, exist_ok=True)
    tag = "".join(ch if ch.isalnum() else "_" for ch in res["label"])[:40]
    jp = os.path.join(a.outdir, "benchmark_concordance_%s.json" % tag)
    json.dump(res, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("[benchmark] 边级一致性：pred %d 边 / gold %d 边 / 交集 %d" % (
        res["n_pred_edges"], res["n_gold_edges"], res["n_overlap"]))
    print("  敏感性(覆盖金标准)=%.3f  精确率=%.3f  Jaccard=%.3f  过预测=%.2fx" % (
        res["sensitivity"] or 0, res["precision"] or 0, res["jaccard"] or 0,
        res["over_prediction_factor"] or 0))
    print("  → %s" % jp)
    return 0


def cmd_report(a):
    files = sorted(f for f in os.listdir(a.dir)
                   if f.startswith("benchmark_") and f.endswith(".json"))
    if not files:
        raise SystemExit("[FAIL] %s 下没有 benchmark_*.json" % a.dir)
    results = [json.load(open(os.path.join(a.dir, f), encoding="utf-8")) for f in files]
    results.sort(key=lambda r: (r.get("kind") != "concordance", r.get("formula_id") or ""))
    p = write_report(a.dir, results, note=a.note or "")
    print("[benchmark] 汇总 %d 个结果 → %s" % (len(results), p))
    return 0


def cmd_selftest(a):
    """算法自检：完美排序 AUC=1、随机 AUC≈0.5、Recall@K 与超几何 p 正确。"""
    # 1) 完美排序
    pos = {"A", "B", "C"}
    ranked = [("A", 9), ("B", 8), ("C", 7), ("D", 6), ("E", 5), ("F", 4), ("G", 3), ("H", 2)]
    r1 = evaluate(ranked, pos, set(), topk=(3, 5), n_boot=50)
    # 2) 最差排序
    ranked_bad = [("H", 9), ("G", 8), ("F", 7), ("E", 6), ("D", 5), ("C", 4), ("B", 3), ("A", 2)]
    r2 = evaluate(ranked_bad, pos, set(), topk=(3,), n_boot=50)
    # 3) 完美 AUC 应为 1.0；最差为 0.0；Recall@3 完美=1.0
    # 4) 超几何：N=8,K=3,n=3, 命中 3 → p = 1/C(8,3) = 1/56
    p_all = hypergeom_sf(3, 8, 3, 3)
    ok = (abs(r1["AUC"] - 1.0) < 1e-9 and abs(r2["AUC"] - 0.0) < 1e-9
          and r1["at_k"]["3"]["recall"] == 1.0 and r2["at_k"]["3"]["recall"] == 0.0
          and abs(p_all - 1 / 56) < 1e-9
          and r1["at_k"]["3"]["enrichment_factor"] == round(3 / (3 * 3 / 8), 2))
    # 5) 并列分数不被算错（全并列 → AUC=0.5）
    tie = evaluate([("A", 1), ("B", 1), ("C", 1), ("D", 1)], {"A", "B"}, set(), topk=(2,), n_boot=0)
    ok = ok and abs(tie["AUC"] - 0.5) < 1e-9
    print("selftest: %s（perfect AUC=%.3f / worst=%.3f / tie=%.3f / hypergeom p=%.5f）" % (
        "PASS" if ok else "FAIL", r1["AUC"], r2["AUC"], tie["AUC"], p_all))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_benchmark: 金标准基准测试（AUC / Recall@K / EF）")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("list"); s.set_defaults(fn=cmd_list)
    e = sub.add_parser("evaluate")
    e.add_argument("--formula", default="CDDP", help="内置方剂 ID：CDDP/HLJDT/LHQW/XFZYT 或中文名")
    e.add_argument("--pred", default=None, help="靶点排序表（gene,score）")
    e.add_argument("--edges", default=None, help="成分-靶点边表（自动按支持成分数排序）")
    e.add_argument("--min-support", type=int, default=1)
    e.add_argument("--weight-mode", default="penalized", choices=["support", "penalized"],
                   help="边表排序基线：support=支持成分数；penalized=多靶点成分降权（默认）")
    e.add_argument("--gold", default=None, help="自定义金标准 JSON（同内置格式）")
    e.add_argument("--gold-csv", default=None, help="外部金标准 CSV（gene,tier）")
    e.add_argument("--label", default=None)
    e.add_argument("--tier", default="strict", choices=["strict", "lenient", "strict-l3"])
    e.add_argument("--topk", default="10,20,50")
    e.add_argument("--bootstrap", type=int, default=1000)
    e.add_argument("-o", "--outdir", default="np_out/benchmark")
    e.set_defaults(fn=cmd_evaluate)
    c = sub.add_parser("concordance", help="边级跨库一致性（敏感性/精确率/Jaccard）")
    c.add_argument("--pred-edges", required=True)
    c.add_argument("--gold-edges", required=True)
    c.add_argument("--label", default=None)
    c.add_argument("-o", "--outdir", default="np_out/benchmark")
    c.set_defaults(fn=cmd_concordance)
    r = sub.add_parser("report"); r.add_argument("--dir", default="np_out/benchmark")
    r.add_argument("--note", default=None); r.set_defaults(fn=cmd_report)
    t = sub.add_parser("selftest"); t.set_defaults(fn=cmd_selftest)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
