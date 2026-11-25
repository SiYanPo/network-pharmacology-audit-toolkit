#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_evidence_ledger.py — 证据账本与校准打分（v3.0）

为什么需要这个模块（评审报告 §1，七寸之一）
------------------------------------------------
此前流程产出：交集表（按支持成分数排）、hub 表（按拓扑共识排）、MR（按 β 排）、
COLOC、协同、对接……每一块各有指标，但**没有任何一步把它们合成"该靶点值得做实验"的判断**。
方法学审稿人的第一个问题必然是："你的最终排序函数是什么？"

本模块给出唯一入口：

    score_targets(outdir) -> [{gene, p_gold, ci_low, ci_high, n_evidence_types, decision, ...}]

三个不可妥协的设计
------------------
1. **账本先于打分**：每一行是"（靶点 × 证据类型）"，带 level/direction/effect_size/uncertainty/
   provenance_hash。打分只能读账本，不能绕过它去读原始文件（否则又变成各模块自说自话）。
2. **概率必须校准，且必须有 CI**：输出 P(gold) 与 bootstrap 95% CI；同时给 Brier/ECE/校准曲线/DCA。
   没校准过的分数只能叫"排序分"，不能叫"概率"——本模块在缺模型时**拒绝**输出概率（返回 UNCALIBRATED）。
3. **决策规则写死**：P≥0.70 且 ≥2 类独立证据 → `wetlab_candidate`；P≥0.70 但单一证据 → `hypothesis`；
   其余 → `insufficient`。规则随模型一起落盘，报告只引用，不在模板里另写一套。

证据分级（L0–L4）
-----------------
    L0 文本挖掘/共现（最低，仅提示）
    L1 计算预测（对接、网络邻近、同源推断）
    L2 数据库注释（TCMSP/ETCM2/ChEMBL 记录、通路成员、组织表达）
    L3 实验支持（结合/活性测定、共表达+实验 PPI 证据）
    L4 遗传因果（MR、共定位）

注意：级别是"证据强度先验"，不是"真值"。模型会学到级别的实际权重；若某级别权重为负，
说明该级别在本基准上是噪声——报告必须如实展示系数。

用法
----
    python np_evidence_ledger.py build  --dir np_out -o np_out/09_evidence
    python np_evidence_ledger.py fit    --ledger np_out/09_evidence/evidence_ledger.csv \
                                        --gold np_out/benchmark/benchmark_CDDP_hits.csv \
                                        -o np_out/09_evidence
    python np_evidence_ledger.py score  --dir np_out --model np_out/09_evidence/calibration_model.json \
                                        -o np_out/09_evidence
    python np_evidence_ledger.py calib  --model ... --ledger ... --gold ...
    python np_evidence_ledger.py selftest
"""

import argparse
import csv
import glob
import hashlib
import json
import math
import os
import random
import sys
import time

VERSION = "3.1"

# 校准器最低可用线：OOF AUC 低于此值时，概率不具区分度，必须拒绝输出 p_gold。
# 2026-11-22 实测：真实 24 靶点 / 13 特征 / 10 正例，OOF AUC=0.471、系数全 0；
# 旧实现仍标 status=OK，score 会给出"概率"——这是把噪声包装成概率，必须拦。
MIN_CV_AUC = 0.60

# 证据级别 → 先验强度（仅用于展示与缺模型时的兜底排序，不参与校准概率）
LEVELS = {
    "L0": ("textmining", 0.10),
    "L1": ("computational", 0.30),
    "L2": ("database", 0.55),
    "L3": ("experimental", 0.80),
    "L4": ("genetic_causal", 1.00),
}
LEVEL_ORDER = ["L0", "L1", "L2", "L3", "L4"]

# 决策规则（写死在此，报告与门禁只引用，不另写）
DECISION_RULES = {
    "wetlab_candidate": {"p_min": 0.70, "min_evidence_types": 2,
                         "desc": "P(gold)≥0.70 且 ≥2 类独立证据 → 可进入湿实验候选"},
    "hypothesis": {"p_min": 0.70, "min_evidence_types": 1,
                   "desc": "P(gold)≥0.70 但仅单一证据类 → 假设生成，需补第二类证据"},
    "insufficient": {"p_min": 0.0, "min_evidence_types": 0,
                     "desc": "P(gold)<0.70 → 证据不足，不得写入结论"},
}

FEATURES = [
    "log_n_compounds",      # 支持该靶点的成分数（log1p）
    "degree_pct",           # PPI 度分位（0–1）
    "betweenness_pct",      # 介数分位
    "hub_robust",           # 拓扑共识判定为"稳健"
    "mr_significant",       # MR p<0.05
    "mr_direction",         # MR β 符号（+1/-1/0）
    "coloc_pp4",            # 最大 PP.H4
    "tissue_expressed",     # HPA 组织表达
    "enrich_hit",           # 落在显著通路
    "ortholog_ok",          # 有直系同源（湿实验可行性）
    "dock_norm",            # 对接分数归一化
    "n_evidence_types",     # 独立证据类计数
    "has_L4",               # 是否存在遗传因果证据
]
FEATURE_ZH = {
    "log_n_compounds": "支持成分数(log)", "degree_pct": "PPI度分位",
    "betweenness_pct": "介数分位", "hub_robust": "拓扑共识稳健",
    "mr_significant": "MR显著", "mr_direction": "MR方向", "coloc_pp4": "共定位PP.H4",
    "tissue_expressed": "组织表达", "enrich_hit": "通路命中", "ortholog_ok": "直系同源",
    "dock_norm": "对接分数", "n_evidence_types": "证据类数", "has_L4": "有遗传因果",
}


# ============================ 基础 IO ============================
def _read_csv(path):
    if not path or not os.path.exists(path):
        return [], []
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        rows = list(csv.DictReader(f, delimiter=delim))
        hdrs = list(rows[0].keys()) if rows else [h.strip() for h in first.strip().split(delim)]
    return rows, hdrs


def _sha16(path):
    try:
        return hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]
    except Exception:
        return ""


def _f(x, default=0.0):
    try:
        v = float(x)
        return v if v == v else default          # 过滤 NaN
    except (TypeError, ValueError):
        return default


def _pick(hdrs, cands):
    low = {h.lower().strip(): h for h in hdrs}
    for c in cands:
        if c in low:
            return low[c]
    return None


def _pct_rank(vals):
    """值 → 分位（0–1，同值取平均秩）。空输入返回空 dict。"""
    n = len(vals)
    if not n:
        return {}
    order = sorted(range(n), key=lambda i: vals[i])
    out, i = {}, 0
    while i < n:
        j = i
        while j + 1 < n and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 / max(1, n - 1) if n > 1 else 0.5
        for k in range(i, j + 1):
            out[order[k]] = avg
        i = j + 1
    return out


# ============================ 账本构建 ============================
def _add(rows, gene, source, etype, level, direction="", effect=None, se=None,
         detail="", prov=""):
    gene = (gene or "").strip().upper()
    if not gene or level not in LEVELS:
        return
    rows.append({
        "gene": gene, "source": source, "evidence_type": etype,
        "level": level, "level_name": LEVELS[level][0], "direction": direction,
        "effect_size": "" if effect is None else round(float(effect), 6),
        "uncertainty": "" if se is None else round(float(se), 6),
        "detail": detail, "provenance_hash": prov,
    })


def collect(outdir):
    """从磁盘产物汇总账本。缺失的模块不补行——账本必须反映真实存在的证据。"""
    rows, sources = [], {}

    # 1) 交集（支持成分数 = L2 数据库注释；disease_score = 疾病侧证据）
    p = os.path.join(outdir, "04_intersect", "intersection.csv")
    r, h = _read_csv(p)
    if r:
        gc = _pick(h, ["gene", "symbol", "target"])
        nc = _pick(h, ["n_compounds"])
        dc = _pick(h, ["disease_score"])
        prov = _sha16(p)
        sources["intersection"] = {"file": p, "rows": len(r), "sha256_16": prov}
        for x in r:
            g = x.get(gc, "")
            _add(rows, g, "compound_target_network", "compound_support", "L2",
                 "support", _f(x.get(nc)), None,
                 "%s 个成分支持" % x.get(nc, "?"), prov)
            if dc and _f(x.get(dc)) > 0:
                _add(rows, g, "disease_targets", "disease_association", "L2",
                     "association", _f(x.get(dc)), None,
                     "疾病侧得分 %s" % x.get(dc), prov)

    # 2) PPI 拓扑（L2：数据库注释的网络位置；不当作实验证据）
    p = os.path.join(outdir, "05_ppi", "hubs.csv")
    r, h = _read_csv(p)
    if r:
        nc = _pick(h, ["node", "gene", "symbol"])
        dg = _pick(h, ["degree"])
        bt = _pick(h, ["betweenness"])
        vd = _pick(h, ["verdict"])
        prov = _sha16(p)
        sources["hubs"] = {"file": p, "rows": len(r), "sha256_16": prov}
        for x in r:
            _add(rows, x.get(nc, ""), "string_ppi", "network_topology", "L2",
                 "", _f(x.get(dg)), None,
                 "degree=%s verdict=%s" % (x.get(dg, "?"), x.get(vd, "")), prov)

    # 3) 富集命中（L2：通路成员关系）
    eg = glob.glob(os.path.join(outdir, "**", "enrich_all.csv"), recursive=True)
    if eg:
        r, h = _read_csv(eg[0])
        tc = _pick(h, ["term", "term_name", "name"])
        qc = _pick(h, ["qvalue", "q_value", "padj"])
        # 富集表通常不含基因列；若有 genes 列则展开
        gc = _pick(h, ["genes", "gene", "intersection"])
        prov = _sha16(eg[0])
        n_sig = 0
        for x in r:
            if qc and _f(x.get(qc), 1) <= 0.05:
                n_sig += 1
                for g in str(x.get(gc) or "").replace(";", ",").split(","):
                    _add(rows, g, "pathway_enrichment", "pathway_membership", "L2",
                         "", _f(x.get(qc)), None, "term=%s" % x.get(tc, ""), prov)
        sources["enrich"] = {"file": eg[0], "n_significant_terms": n_sig, "sha256_16": prov}

    # 4) MR（L4）+ COLOC（L4）
    for mp in glob.glob(os.path.join(outdir, "**", "mr_result.json"), recursive=True):
        try:
            j = json.load(open(mp, encoding="utf-8"))
        except Exception:
            continue
        if j.get("status") != "OK":
            continue
        gene = j.get("gene") or os.path.basename(os.path.dirname(mp))
        ivw = j.get("IVW") or {}
        prov = _sha16(mp)
        _add(rows, gene, "mendelian_randomization", "causal_estimate", "L4",
             "risk_increasing" if _f(ivw.get("beta")) > 0 else "risk_decreasing",
             _f(ivw.get("beta")), _f(ivw.get("se")),
             "IVW p=%.3g F=%s" % (_f(ivw.get("p"), 1), j.get("F_statistic", "?")), prov)
        col = j.get("coloc") or {}
        if col.get("status") == "OK":
            _add(rows, gene, "coloc_abf", "colocalisation", "L4", "",
                 _f(col.get("PP.H4")), None, "ABF 近似 PP.H4", prov)

    # 5) 组织表达（L1：上下文，不作因果）
    em = os.path.join(outdir, "07_enrich", "enrich_meta.json")
    if os.path.exists(em):
        try:
            jm = json.load(open(em, encoding="utf-8"))
        except Exception:
            jm = {}
        if jm.get("background_mode") in ("hpa", "deg"):
            sources["expression_background"] = {
                "file": em, "mode": jm.get("background_mode"),
                "N": jm.get("N"), "sha256_16": _sha16(em)}
            # 逐基因组织表达：否则 tissue_expressed 永远是 0，等于白留一个特征
            tissue = jm.get("hpa_tissue") or ""
            tsv = _find_hpa_tsv(outdir)
            if tsv and tissue:
                try:
                    import np_enrich as EN
                    bg = EN.hpa_background(tsv, tissue, float(jm.get("hpa_min_ntpm") or 1.0))
                except Exception:
                    bg = None
                if bg:
                    prov = _sha16(tsv)
                    sources["expression_background"]["tsv"] = tsv
                    sources["expression_background"]["n_expressed"] = len(bg)
                    # 只标注**候选宇宙内**的基因：HPA 背景有 1.3 万个基因，
                    # 若全量写入，账本会被"非候选靶点"淹没（实测 361 行 → 13900 行）。
                    universe = {r["gene"] for r in rows}
                    for g in bg:
                        if g in universe:
                            _add(rows, g, "hpa_tissue_expression", "expression_context", "L2",
                                 "expressed", 1.0, None, "tissue=%s" % tissue, prov)

    # 6) 同源（L1：湿实验可行性，不是靶点证据）
    for op in glob.glob(os.path.join(outdir, "**", "*_orthologs.csv"), recursive=True):
        r, h = _read_csv(op)
        gcl = _pick(h, ["human"])
        st = _pick(h, ["status"])
        mt = _pick(h, ["mapping_type"])
        prov = _sha16(op)
        for x in r:
            ok = str(x.get(st) or "").startswith("resolved") and (x.get(mt) or "ortholog") == "ortholog"
            _add(rows, x.get(gcl, ""), "ortholog_mapping", "wetlab_feasibility", "L1",
                 "direct" if ok else "analog", 1.0 if ok else 0.0, None,
                 "status=%s type=%s" % (x.get(st, ""), x.get(mt, "")), prov)

    # 7) 对接（L1：计算预测）
    for dp in glob.glob(os.path.join(outdir, "**", "dock_scores.csv"), recursive=True) + \
            glob.glob(os.path.join(outdir, "**", "docking_results.csv"), recursive=True):
        r, h = _read_csv(dp)
        gc = _pick(h, ["gene", "target", "receptor"])
        ac = _pick(h, ["affinity", "score", "best_affinity"])
        prov = _sha16(dp)
        for x in r:
            _add(rows, x.get(gc, ""), "molecular_docking", "binding_prediction", "L1",
                 "binding", _f(x.get(ac)), None, "docking score", prov)

    # 8) 干预方向一致性（L4 的反向证据：矛盾时记负向）
    ip = os.path.join(outdir, "08_intervention", "intervention_report.json")
    if not os.path.exists(ip):
        cand = glob.glob(os.path.join(outdir, "**", "intervention*.json"), recursive=True)
        ip = cand[0] if cand else ""
    if ip and os.path.exists(ip):
        try:
            j = json.load(open(ip, encoding="utf-8"))
        except Exception:
            j = {}
        for r in (j.get("rows") or []):
            if r.get("verdict") in ("Target_Intervention_Mismatch", "Directional_Conflict"):
                _add(rows, r.get("gene", ""), "intervention_check", "direction_conflict", "L4",
                     "conflict", -1.0, None, r.get("message", "")[:120], _sha16(ip))

    return rows, sources



def _find_hpa_tsv(outdir):
    """在常见位置找 HPA bulk TSV：缓存目录 / 输出目录上层 / 同目录 hpa_cache。"""
    cands = []
    try:
        import npcore as C
        d = C.cache_dir("hpa")
        if d:
            cands.append(os.path.join(d, "rna_tissue_consensus.tsv"))
    except Exception:
        pass
    cands += [os.path.join(outdir, "..", "hpa_cache", "rna_tissue_consensus.tsv"),
              os.path.join(outdir, "hpa_cache", "rna_tissue_consensus.tsv"),
              os.path.join(os.path.expanduser("~"), ".np_cache", "hpa",
                           "rna_tissue_consensus.tsv")]
    for c in cands:
        if c and os.path.exists(c):
            return os.path.abspath(c)
    return None


def ledger_path(outdir):
    return os.path.join(outdir, "09_evidence", "evidence_ledger.csv")


def write_ledger(rows, outdir):
    d = os.path.join(outdir, "09_evidence")
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, "evidence_ledger.csv")
    cols = ["gene", "source", "evidence_type", "level", "level_name", "direction",
            "effect_size", "uncertainty", "detail", "provenance_hash"]
    with open(p, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    return p


# ============================ 特征工程 ============================
def build_features(ledger_rows, outdir=None):
    """账本 → 每基因特征字典。**只从账本推导**（除拓扑分位需要原始 degree 分布）。"""
    by_gene = {}
    for r in ledger_rows:
        g = r["gene"]
        d = by_gene.setdefault(g, {"levels": set(), "types": set(), "n_compounds": 0,
                                   "degree": 0.0, "betw": 0.0, "robust": 0,
                                   "mr_p": None, "mr_beta": None, "pp4": 0.0,
                                   "tissue": 0, "enrich": 0, "ortholog": 0, "dock": []})
        d["levels"].add(r["level"])
        d["types"].add(r["evidence_type"])
        et = r["evidence_type"]
        if et == "compound_support":
            d["n_compounds"] = max(d["n_compounds"], _f(r["effect_size"]))
        elif et == "network_topology":
            d["degree"] = max(d["degree"], _f(r["effect_size"]))
            if "稳健" in (r.get("detail") or ""):
                d["robust"] = 1
            m = (r.get("detail") or "")
            if "verdict=" in m:
                pass
        elif et == "causal_estimate":
            d["mr_beta"] = _f(r["effect_size"])
            pv = None
            det = r.get("detail") or ""
            if "p=" in det:
                try:
                    pv = float(det.split("p=")[1].split()[0])
                except Exception:
                    pv = None
            d["mr_p"] = pv
        elif et == "colocalisation":
            d["pp4"] = max(d["pp4"], _f(r["effect_size"]))
        elif et == "pathway_membership":
            d["enrich"] = 1
        elif et == "expression_context":
            d["tissue"] = 1
        elif et == "wetlab_feasibility":
            d["ortholog"] = 1 if _f(r["effect_size"]) >= 1 else 0
        elif et == "binding_prediction":
            d["dock"].append(_f(r["effect_size"]))

    # 介数分位需要全表分布
    betw = {}
    if outdir:
        r, h = _read_csv(os.path.join(outdir, "05_ppi", "hubs.csv"))
        nc = _pick(h, ["node", "gene"])
        bc = _pick(h, ["betweenness"])
        if nc and bc:
            for x in r:
                betw[(x.get(nc) or "").strip().upper()] = _f(x.get(bc))

    genes = sorted(by_gene)
    degs = [by_gene[g]["degree"] for g in genes]
    deg_pct = {g: v for g, v in zip(genes, [_pct_rank(degs)[i] for i in range(len(genes))])}
    bl = [betw.get(g, 0.0) for g in genes]
    brank = _pct_rank(bl) if bl else {}
    dock_vals = [max(by_gene[g]["dock"]) if by_gene[g]["dock"] else 0.0 for g in genes]
    dmin, dmax = (min(dock_vals), max(dock_vals)) if dock_vals else (0, 1)

    out = {}
    for i, g in enumerate(genes):
        d = by_gene[g]
        dv = max(d["dock"]) if d["dock"] else 0.0
        out[g] = {
            "log_n_compounds": math.log1p(max(0.0, d["n_compounds"])),
            "degree_pct": round(deg_pct.get(g, 0.0), 6),
            "betweenness_pct": round(brank.get(i, 0.0), 6) if brank else 0.0,
            "hub_robust": float(d["robust"]),
            "mr_significant": 1.0 if (d["mr_p"] is not None and d["mr_p"] < 0.05) else 0.0,
            "mr_direction": (1.0 if (d["mr_beta"] or 0) > 0 else (-1.0 if (d["mr_beta"] or 0) < 0 else 0.0)),
            "coloc_pp4": round(d["pp4"], 6),
            "tissue_expressed": float(d["tissue"]),
            "enrich_hit": float(d["enrich"]),
            "ortholog_ok": float(d["ortholog"]),
            "dock_norm": round((dv - dmin) / (dmax - dmin), 6) if dmax > dmin else 0.0,
            "n_evidence_types": float(len(d["types"])),
            "has_L4": 1.0 if "L4" in d["levels"] else 0.0,
        }
    return out, {"n_genes": len(out), "betweenness_available": bool(betw and any(betw.values()))}


# ============================ 纯 Python 逻辑回归（IRLS） ============================
def _solve(A, b):
    """高斯消元（带部分主元）解 A x = b。A 为 n×n 列表。"""
    n = len(A)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-12:
            return None
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        for j in range(c, n + 1):
            M[c][j] /= pv
        for r in range(n):
            if r != c and abs(M[r][c]) > 0:
                f = M[r][c]
                for j in range(c, n + 1):
                    M[r][j] -= f * M[c][j]
    return [M[i][n] for i in range(n)]


def _sigmoid(z):
    if z >= 0:
        return 1.0 / (1.0 + math.exp(-z))
    e = math.exp(z)
    return e / (1.0 + e)


def fit_logistic(X, y, l2=1.0, iters=60, eps=1e-8):
    """IRLS 拟合 L2 正则逻辑回归。X 已含截距列。返回 (coef, se, info)。"""
    n, p = len(X), len(X[0])
    beta = [0.0] * p
    for it in range(iters):
        eta = [sum(beta[j] * X[i][j] for j in range(p)) for i in range(n)]
        mu = [_sigmoid(z) for z in eta]
        w = [max(1e-6, mu[i] * (1 - mu[i])) for i in range(n)]
        # 梯度与 Hessian（不含截距正则）
        grad = [sum(X[i][j] * (y[i] - mu[i]) for i in range(n)) - (l2 * beta[j] if j else 0.0)
                for j in range(p)]
        H = [[sum(w[i] * X[i][j] * X[i][k] for i in range(n)) + (l2 if (j == k and j) else 0.0)
              for k in range(p)] for j in range(p)]
        step = _solve(H, grad)
        if step is None:
            return None, None, {"status": "SINGULAR"}
        beta = [beta[j] + step[j] for j in range(p)]
        if max(abs(s) for s in step) < eps:
            break
    # 标准差：H^-1 对角
    H = [[sum(max(1e-6, _sigmoid(sum(beta[t] * X[i][t] for t in range(p))) *
                    (1 - _sigmoid(sum(beta[t] * X[i][t] for t in range(p))))) *
           X[i][j] * X[i][k] for i in range(n)) + (l2 if (j == k and j) else 0.0)
          for k in range(p)] for j in range(p)]
    inv_diag = []
    for j in range(p):
        e = [1.0 if k == j else 0.0 for k in range(p)]
        sol = _solve(H, e)
        inv_diag.append(math.sqrt(max(0.0, sol[j])) if sol else float("nan"))
    return beta, inv_diag, {"status": "OK", "iterations": it + 1}


def _design(feat_dict, genes, scaler=None):
    """特征 → 设计矩阵（含截距）。scaler 为 {mu, sd}，缺省时按本批计算并返回。"""
    cols = FEATURES
    raw = [[float(feat_dict[g].get(c, 0.0)) for c in cols] for g in genes]
    if scaler is None:
        mu = [sum(r[j] for r in raw) / max(1, len(raw)) for j in range(len(cols))]
        sd = []
        for j in range(len(cols)):
            v = sum((r[j] - mu[j]) ** 2 for r in raw) / max(1, len(raw) - 1)
            sd.append(math.sqrt(v) if v > 1e-12 else 1.0)
        scaler = {"mu": mu, "sd": sd, "features": cols}
    X = [[1.0] + [(r[j] - scaler["mu"][j]) / scaler["sd"][j] for j in range(len(cols))] for r in raw]
    return X, scaler


# ============================ 校准指标 ============================
def brier(y, p):
    return sum((pi - yi) ** 2 for pi, yi in zip(p, y)) / max(1, len(y))


def ece(y, p, bins=10):
    """期望校准误差：等宽分箱，Σ(n_b/N)·|acc_b − conf_b|。"""
    n = len(y)
    if not n:
        return 0.0, []
    buckets = [[] for _ in range(bins)]
    for yi, pi in zip(y, p):
        k = min(bins - 1, max(0, int(pi * bins)))
        buckets[k].append((pi, yi))
    tot, curve = 0.0, []
    for k, b in enumerate(buckets):
        if not b:
            continue
        conf = sum(x[0] for x in b) / len(b)
        acc = sum(x[1] for x in b) / len(b)
        tot += len(b) / n * abs(acc - conf)
        curve.append({"bin": k, "n": len(b), "conf": round(conf, 4), "acc": round(acc, 4)})
    return tot, curve


def net_benefit(y, p, thresholds=None):
    """决策曲线（DCA）：NB = TP/N − FP/N × pt/(1−pt)。"""
    thresholds = thresholds or [round(0.05 * i, 2) for i in range(1, 20)]
    n = max(1, len(y))
    out = []
    prev = sum(y) / n
    for t in thresholds:
        tp = sum(1 for yi, pi in zip(y, p) if pi >= t and yi == 1)
        fp = sum(1 for yi, pi in zip(y, p) if pi >= t and yi == 0)
        nb = tp / n - fp / n * (t / (1 - t))
        nb_all = prev - (1 - prev) * (t / (1 - t))
        out.append({"threshold": t, "net_benefit": round(nb, 4),
                    "treat_all": round(nb_all, 4), "treat_none": 0.0})
    return out


def auc(y, p):
    """Mann–Whitney U。"""
    pos = [pi for yi, pi in zip(y, p) if yi == 1]
    neg = [pi for yi, pi in zip(y, p) if yi == 0]
    if not pos or not neg:
        return None
    s = 0.0
    for a in pos:
        for b in neg:
            s += 1.0 if a > b else (0.5 if a == b else 0.0)
    return s / (len(pos) * len(neg))


# ============================ 主流程 ============================
def fit(ledger_csv, gold_genes, outdir, l2=1.0, seed=42, n_boot=300):
    """在基准标签上拟合校准模型。gold_genes 为金标准正类集合。"""
    rows, _ = _read_csv(ledger_csv)
    if not rows:
        return {"status": "NO_LEDGER"}
    d = os.path.dirname(os.path.abspath(ledger_csv))
    feats, meta = build_features(rows, os.path.dirname(d) if os.path.basename(d) == "09_evidence" else d)
    gold = {(g or "").strip().upper() for g in gold_genes if g}
    genes = [g for g in feats if g in gold] + [g for g in feats if g not in gold]
    if len(genes) < 10:
        return {"status": "TOO_FEW_GENES", "n": len(genes)}
    y = [1.0 if g in gold else 0.0 for g in genes]
    if sum(y) < 2 or sum(y) > len(y) - 2:
        return {"status": "DEGENERATE_LABELS", "n_pos": int(sum(y)), "n": len(y)}
    X, scaler = _design(feats, genes)
    beta, se, info = fit_logistic(X, y, l2=l2)
    if beta is None:
        return {"status": "FIT_FAILED"}
    p = [_sigmoid(sum(beta[j] * X[i][j] for j in range(len(beta)))) for i in range(len(X))]
    b, ec, curve = brier(y, p), None, None
    ec, curve = ece(y, p)
    model = {
        "version": VERSION, "status": "OK", "features": FEATURES,
        "feature_zh": FEATURE_ZH, "scaler": scaler, "coef": beta, "se": se,
        "l2": l2, "seed": seed, "n_train": len(genes), "n_pos": int(sum(y)),
        "train_metrics": {"brier": round(b, 4), "ece": round(ec, 4),
                          "auc": None if auc(y, p) is None else round(auc(y, p), 4)},
        "reliability_curve": curve, "dca": net_benefit(y, p),
        "decision_rules": DECISION_RULES,
        "fitted_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "caveat": ("这是**在基准标签上拟合的校准器**，其数值只在该基准的定义域内有效；"
                   "换疾病/换数据库后必须先重拟合，否则概率不可解释。"),
    }
    # bootstrap 系数 → 预测区间
    rnd = random.Random(seed)
    idx_all = list(range(len(genes)))
    boots = []
    for _ in range(n_boot):
        idx = [rnd.choice(idx_all) for _ in idx_all]
        if len({y[i] for i in idx}) < 2:
            continue
        bb, _, _ = fit_logistic([X[i] for i in idx], [y[i] for i in idx], l2=l2, iters=30)
        if bb:
            boots.append(bb)
    cv = cv_calibration(feats, gold, k=5, l2=l2, seed=seed)
    model["cv"] = cv
    model["n_bootstrap"] = len(boots)
    model["bootstrap"] = boots[:200]        # 控制文件体积

    # 校准可用性闸门：OOF AUC 太低或系数全 0 → 概率不可解释，拒绝标 OK。
    cv_auc = cv.get("auc") if isinstance(cv, dict) else None
    weak_reasons = []
    if not isinstance(cv, dict) or cv.get("status") != "OK":
        weak_reasons.append("CV_NOT_OK:%s" % (cv.get("status") if isinstance(cv, dict) else "NA"))
    elif cv_auc is None or float(cv_auc) < MIN_CV_AUC:
        weak_reasons.append("OOF_AUC_TOO_LOW:%s<%.2f" % (cv_auc, MIN_CV_AUC))
    if all(abs(float(b)) < 1e-9 for b in beta[1:]):
        weak_reasons.append("ALL_COEFFICIENTS_ZERO")
    if weak_reasons:
        model["status"] = "WEAK"
        model["weak_reasons"] = weak_reasons
        model["caveat"] = ("校准器未通过可用性闸门（%s）：只能作为诊断工件，"
                           "score 将返回 UNCALIBRATED，不得把 p_gold 写成概率。"
                           % "; ".join(weak_reasons))
    os.makedirs(outdir, exist_ok=True)
    mp = os.path.join(outdir, "calibration_model.json")
    json.dump(model, open(mp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return {"status": "OK" if model["status"] == "OK" else "WEAK_CALIBRATION",
            "model": mp, "metrics_train": model["train_metrics"],
            "metrics_cv": {k: v for k, v in cv.items() if k in ("brier", "ece", "auc", "n_oof")},
            "n_train": len(genes), "n_pos": int(sum(y)),
            "weak_reasons": weak_reasons,
            "coefficients": dict(zip(FEATURES, beta[1:]))}


def _predict(model, feats, genes):
    sc = model["scaler"]
    X, _ = _design(feats, genes, scaler=sc)
    beta = model["coef"]
    p = [_sigmoid(sum(beta[j] * X[i][j] for j in range(len(beta)))) for i in range(len(X))]
    return p, X


def cv_calibration(feats, gold, k=5, l2=1.0, seed=42):
    """分层 k 折**样本外**校准指标。

    为什么必须有：样本内 Brier/ECE 只会告诉你"模型能拟合训练集"。审稿人要的是
    out-of-fold 校准；若 OOF ECE 明显大于样本内，说明模型在过拟合，概率不可用。
    """
    genes = [g for g in feats if g in gold] + [g for g in feats if g not in gold]
    y = [1.0 if g in gold else 0.0 for g in genes]
    pos = [i for i, v in enumerate(y) if v == 1]
    neg = [i for i, v in enumerate(y) if v == 0]
    rnd = random.Random(seed)
    rnd.shuffle(pos)
    rnd.shuffle(neg)
    folds = [[] for _ in range(k)]
    for arr in (pos, neg):
        for j, i in enumerate(arr):
            folds[j % k].append(i)
    oof = [None] * len(genes)
    for fi in range(k):
        te = folds[fi]
        tr = [i for i in range(len(genes)) if i not in set(te)]
        if len({y[i] for i in tr}) < 2 or not te:
            continue
        Xtr, scaler = _design(feats, [genes[i] for i in tr])
        beta, _, _ = fit_logistic(Xtr, [y[i] for i in tr], l2=l2)
        if beta is None:
            continue
        Xte, _ = _design(feats, [genes[i] for i in te], scaler=scaler)
        for row, i in enumerate(te):
            oof[i] = _sigmoid(sum(beta[j] * Xte[row][j] for j in range(len(beta))))
    pairs = [(y[i], oof[i]) for i in range(len(genes)) if oof[i] is not None]
    if len(pairs) < 10:
        return {"status": "TOO_FEW_OOF", "n": len(pairs)}
    yy = [a for a, _ in pairs]
    pp = [b for _, b in pairs]
    e_val, curve = ece(yy, pp)
    return {"status": "OK", "k": k, "n_oof": len(pairs),
            "brier": round(brier(yy, pp), 4), "ece": round(e_val, 4),
            "auc": None if auc(yy, pp) is None else round(auc(yy, pp), 4),
            "reliability_curve": curve,
            "note": "out-of-fold 指标；与样本内指标的差值即过拟合幅度"}


def score(outdir, model_path=None, ledger=None):
    """唯一打分入口。无模型时**拒绝**输出概率（返回 UNCALIBRATED + 先验分）。"""
    lc = ledger or ledger_path(outdir)
    rows, _ = _read_csv(lc)
    if not rows:
        return {"status": "NO_LEDGER", "targets": []}
    feats, fmeta = build_features(rows, outdir)
    genes = list(feats)
    model = None
    if model_path and os.path.exists(model_path):
        try:
            model = json.load(open(model_path, encoding="utf-8"))
        except Exception:
            model = None
    if not model or model.get("status") != "OK":
        # 兜底：级别先验加权，明确标 UNCALIBRATED
        prior = {}
        for r in rows:
            g = r["gene"]
            w = LEVELS.get(r["level"], ("", 0.0))[1]
            prior[g] = max(prior.get(g, 0.0), w)
        targets = [{"gene": g, "score_prior": round(prior.get(g, 0.0), 4),
                    "p_gold": None, "ci_low": None, "ci_high": None,
                    "n_evidence_types": int(feats[g]["n_evidence_types"]),
                    "decision": "uncalibrated",
                    "decision_zh": "未校准：仅可按级别先验排序，不得写作概率",
                    **{k: round(v, 4) for k, v in feats[g].items()}} for g in genes]
        targets.sort(key=lambda x: (-x["score_prior"], x["gene"]))
        reason = "缺 calibration_model.json；先运行 fit 再打分"
        if model and model.get("status") == "WEAK":
            reason = "校准器未通过可用性闸门：%s" % "; ".join(model.get("weak_reasons") or [])
        return {"status": "UNCALIBRATED", "targets": targets, "feature_meta": fmeta,
                "reason": reason}

    p, X = _predict(model, feats, genes)
    boots = model.get("bootstrap") or []
    ci = []
    for i in range(len(genes)):
        vals = []
        for b in boots:
            vals.append(_sigmoid(sum(b[j] * X[i][j] for j in range(len(b)))))
        if len(vals) >= 20:
            vals.sort()
            ci.append((vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals)) - 1]))
        else:
            ci.append((None, None))
    out = []
    for i, g in enumerate(genes):
        lo, hi = ci[i]
        ntypes = int(feats[g]["n_evidence_types"])
        if p[i] >= 0.70 and ntypes >= 2:
            dec, zh = "wetlab_candidate", "可进入湿实验候选"
        elif p[i] >= 0.70:
            dec, zh = "hypothesis", "单一证据类，需补第二类"
        else:
            dec, zh = "insufficient", "证据不足，不得写入结论"
        out.append({"gene": g, "p_gold": round(p[i], 4),
                    "ci_low": None if lo is None else round(lo, 4),
                    "ci_high": None if hi is None else round(hi, 4),
                    "n_evidence_types": ntypes, "decision": dec, "decision_zh": zh,
                    **{k: round(v, 4) for k, v in feats[g].items()}})
    out.sort(key=lambda x: (-x["p_gold"], x["gene"]))
    return {"status": "OK", "targets": out, "feature_meta": fmeta,
            "model": model_path, "n_bootstrap": model.get("n_bootstrap")}


def score_targets(outdir, model_path=None, ledger=None):
    """对外稳定名（DoD 要求）：等价于 score()，便于论文/API 只引用一个入口。"""
    return score(outdir, model_path=model_path, ledger=ledger)


def write_scores(res, outdir):
    d = os.path.join(outdir, "09_evidence")
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, "target_scores.csv")
    if not res.get("targets"):
        return None
    cols = list(dict.fromkeys(
        ["gene", "p_gold", "ci_low", "ci_high", "n_evidence_types", "decision",
         "decision_zh"] + FEATURES))
    with open(p, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in res["targets"]:
            w.writerow(r)
    json.dump({k: v for k, v in res.items() if k != "targets"},
              open(os.path.join(d, "score_manifest.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    json.dump({"decision_rules": DECISION_RULES, "features": FEATURES,
               "feature_zh": FEATURE_ZH, "levels": {k: v[0] for k, v in LEVELS.items()},
               "applied_at": time.strftime("%Y-%m-%d %H:%M:%S"),
               "n_targets": len(res["targets"])},
              open(os.path.join(d, "decision_rules.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    return p


# ============================ 报告片段 ============================
def report_section(res, model_path=None):
    L = ["## 证据聚合与校准（单一排序函数）", ""]
    if res.get("status") == "UNCALIBRATED":
        L += ["> **未校准**：本报告未提供 `calibration_model.json`，因此只给出按证据级别先验的排序，"
              "**不输出概率**。概率必须在基准标签上拟合校准器后才能报告（见 Methods）。",
              "", "| 靶点 | 先验分 | 证据类数 |", "|---|---|---|"]
        for r in res["targets"][:15]:
            L.append("| %s | %.2f | %d |" % (r["gene"], r["score_prior"], r["n_evidence_types"]))
        return "\n".join(L)
    L += ["| 靶点 | P(gold) | 95% CI | 证据类数 | 决策 |", "|---|---|---|---|---|"]
    for r in res["targets"][:20]:
        ci = ("%.2f–%.2f" % (r["ci_low"], r["ci_high"])) if r["ci_low"] is not None else "NA"
        L.append("| %s | %.3f | %s | %d | %s |" % (
            r["gene"], r["p_gold"], ci, r["n_evidence_types"], r["decision_zh"]))
    L += ["", "决策规则（与模型一同落盘，不在模板里另写）：",
          "- `P≥0.70 且 ≥2 类独立证据` → 湿实验候选；",
          "- `P≥0.70 但单一证据类` → 假设生成；",
          "- `P<0.70` → 证据不足，不得写入结论。", ""]
    return "\n".join(L)


# ============================ CLI ============================
def selftest():
    checks = []

    # 1) 账本构建：缺失产物不补行
    import tempfile
    d = tempfile.mkdtemp(prefix="np_ledger_")
    os.makedirs(os.path.join(d, "04_intersect"), exist_ok=True)
    os.makedirs(os.path.join(d, "05_ppi"), exist_ok=True)
    open(os.path.join(d, "04_intersect", "intersection.csv"), "w", encoding="utf-8").write(
        "gene,n_compounds,disease_score\nIL6,10,6.5\nTNF,8,5.0\nFAKE1,1,0.1\n")
    open(os.path.join(d, "05_ppi", "hubs.csv"), "w", encoding="utf-8").write(
        "rank,node,degree,betweenness,verdict\n1,IL6,30,0.13,稳健\n2,TNF,20,0.08,稳健\n3,FAKE1,2,0.0,不稳定\n")
    rows, src = collect(d)
    checks.append(("ledger_rows", len(rows) >= 4, "n=%d" % len(rows)))
    checks.append(("levels_valid", all(r["level"] in LEVELS for r in rows), ""))
    checks.append(("provenance_present", all(r["provenance_hash"] for r in rows), ""))

    # 2) 特征：缺失 dock/ortholog 时为零，不报错
    feats, fmeta = build_features(rows, d)
    checks.append(("features_complete", all(set(FEATURES) <= set(v) for v in feats.values()),
                   "genes=%d" % len(feats)))
    checks.append(("degree_pct_range", all(0.0 <= v["degree_pct"] <= 1.0 for v in feats.values()), ""))

    # 3) 校准指标：完美预测 → Brier 0、ECE 0
    y = [1, 1, 0, 0]
    checks.append(("brier_perfect", abs(brier(y, [1.0, 1.0, 0.0, 0.0])) < 1e-12, ""))
    e_val, curve = ece(y, [1.0, 1.0, 0.0, 0.0])
    checks.append(("ece_perfect", abs(e_val) < 1e-12, "%.3g" % e_val))
    checks.append(("auc_perfect", abs(auc(y, [0.9, 0.8, 0.2, 0.1]) - 1.0) < 1e-12, ""))

    # 4) 逻辑回归：可分数据应收敛到高 AUC
    rnd = random.Random(7)
    X, yy = [], []
    for _ in range(120):
        a = rnd.gauss(0, 1)
        X.append([1.0, a, rnd.gauss(0, 0.3)])
        yy.append(1.0 if a + rnd.gauss(0, 0.4) > 0 else 0.0)
    beta, se, info = fit_logistic(X, yy, l2=0.5)
    pp = [_sigmoid(sum(beta[j] * X[i][j] for j in range(len(beta)))) for i in range(len(X))]
    a = auc(yy, pp)
    checks.append(("logistic_fit", beta is not None and info["status"] == "OK", info.get("status", "")))
    checks.append(("logistic_auc", a is not None and a > 0.85, "AUC=%.3f" % (a or 0)))
    checks.append(("prob_in_range", all(0.0 <= v <= 1.0 for v in pp), ""))

    # 5) 无模型时拒绝输出概率
    write_ledger(rows, d)
    res = score(d, model_path=None)
    checks.append(("refuse_uncalibrated",
                   res["status"] == "UNCALIBRATED" and all(t["p_gold"] is None for t in res["targets"]),
                   res["status"]))

    # 6) 决策规则阈值与文档一致
    checks.append(("decision_rules", DECISION_RULES["wetlab_candidate"]["p_min"] == 0.70
                   and DECISION_RULES["wetlab_candidate"]["min_evidence_types"] == 2, ""))

    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for n, g, d_ in checks:
        print("  [%s] %-26s %s" % ("✓" if g else "✗", n, d_))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = argparse.ArgumentParser(description="np_evidence_ledger: 证据账本 + 校准打分")
    sub = ap.add_subparsers(dest="fn")

    b = sub.add_parser("build", help="从产物汇总账本")
    b.add_argument("--dir", required=True)

    f = sub.add_parser("fit", help="在基准标签上拟合校准器")
    f.add_argument("--ledger", required=True)
    f.add_argument("--gold", required=True, help="金标准：基因列表文件或 *_hits.csv（label=1）")
    f.add_argument("-o", "--outdir", required=True)
    f.add_argument("--l2", type=float, default=1.0)
    f.add_argument("--bootstrap", type=int, default=300)

    s = sub.add_parser("score", help="打分（唯一入口）")
    s.add_argument("--dir", required=True)
    s.add_argument("--model", default=None)
    s.add_argument("--ledger", default=None)
    s.add_argument("-o", "--outdir", default=None)

    au = sub.add_parser("auto", help="build + score（有模型则校准打分，否则 UNCALIBRATED）")
    au.add_argument("--dir", required=True)
    au.add_argument("--model", default=None)

    c = sub.add_parser("calib", help="只算校准指标")
    c.add_argument("--ledger", required=True)
    c.add_argument("--gold", required=True)
    c.add_argument("--model", required=True)

    a = ap.parse_args(argv)
    if a.fn == "build":
        rows, src = collect(a.dir)
        p = write_ledger(rows, a.dir)
        json.dump({"n_rows": len(rows), "sources": src, "built_at": time.strftime("%Y-%m-%d %H:%M:%S")},
                  open(os.path.join(a.dir, "09_evidence", "ledger_manifest.json"), "w", encoding="utf-8"),
                  ensure_ascii=False, indent=1)
        print("[ledger] %d 行 → %s" % (len(rows), p))
        return 0
    if a.fn == "fit":
        g = read_gold(a.gold)
        r = fit(a.ledger, g, a.outdir, l2=a.l2, n_boot=a.bootstrap)
        print("[fit] %s" % json.dumps(r, ensure_ascii=False)[:600])
        return 0 if r["status"] == "OK" else 2
    if a.fn == "score":
        res = score(a.dir, a.model or os.path.join(a.dir, "09_evidence", "calibration_model.json"),
                    ledger=a.ledger)
        p = write_scores(res, a.outdir or a.dir)
        print("[score] %s：%d 个靶点 → %s" % (res["status"], len(res.get("targets") or []), p))
        return 0 if res["status"] == "OK" else 3
    if a.fn == "auto":
        rows, src = collect(a.dir)
        lp = write_ledger(rows, a.dir)
        model = a.model or os.path.join(a.dir, "09_evidence", "calibration_model.json")
        res = score(a.dir, model_path=model, ledger=lp)
        sp = write_scores(res, a.dir)
        print("[evidence] 账本 %d 行 → %s；打分 %s（%s，%d 靶点）" % (
            len(rows), os.path.basename(lp), res["status"], os.path.basename(sp or ""),
            len(res.get("targets") or [])))
        return 0
    if a.fn == "calib":
        rows, _ = _read_csv(a.ledger)
        model = json.load(open(a.model, encoding="utf-8"))
        g = read_gold(a.gold)
        feats, _ = build_features(rows)
        genes = [x for x in feats if x in g] + [x for x in feats if x not in g]
        y = [1.0 if x in g else 0.0 for x in genes]
        p, _ = _predict(model, feats, genes)
        e_val, curve = ece(y, p)
        print("[calib] Brier=%.4f ECE=%.4f AUC=%s n=%d" % (
            brier(y, p), e_val, auc(y, p), len(y)))
        return 0
    ap.print_help()
    return 0


def read_gold(path):
    """金标准读取：纯列表文件，或 np_benchmark 的 *_hits.csv（label=1 记为阳性）。"""
    out = set()
    if not path or not os.path.exists(path):
        return out
    rows, h = _read_csv(path)
    if rows and _pick(h, ["label", "y", "is_positive"]):
        gc = _pick(h, ["gene", "target", "node"])
        lc = _pick(h, ["label", "y", "is_positive"])
        for r in rows:
            if _f(r.get(lc)) >= 1:
                out.add((r.get(gc) or "").strip().upper())
        return out
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            t = line.strip().split("\t")[0].split(",")[0].strip().upper()
            if t and t not in ("GENE", "TARGET"):
                out.add(t)
    return out


if __name__ == "__main__":
    sys.exit(main())
