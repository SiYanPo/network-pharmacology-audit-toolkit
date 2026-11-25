#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_funnel.py —— 大复方"先筛选后对接"漏斗（解决 50+ 味药、上千成分的 Vina 计算灾难）

问题
    大复方可能筛出 500–1000 个活性成分，若全部对接 4–5 个靶点 = 数千次 Vina，
    单机需数天。且大量成分度值低、暴露量不可达，对接意义有限。

漏斗（四级，逐级收窄；每级记录淘汰原因，可审计、可回溯）
    Stage 0  全量成分 N0
    Stage 1  拓扑筛选：保留 degree ≥ top_degree_ratio 分位（默认 top 30%）
    Stage 2  暴露量筛选：np_dose 估算的游离 Cmax ≥ 阈值（默认 0.01 µM），剔除不可达成分
    Stage 3  ADMET 筛选：BBB/毒性过滤（若提供 ADMET 表），优先血脑屏障可透过者（PD 中枢）
    Stage 4  去重骨架：按 Murcko/简化骨架去重，每骨架保留度值最高者
    → 输出 final_docking_list.csv（Top N），传给 np_dockprep.py

诚实边界
    * 漏斗是"计算预算下的优先级排序"，不是否定被淘汰成分的生物活性。
    * 每级淘汰写入 funnel_trace.csv，论文须报告各层级剩余数与阈值。

用法
    python np_funnel.py --edges edges.tsv --admet admet.csv --dose dose.csv -o out/
    python np_funnel.py --edges edges.tsv --top-ratio 0.2 --max-final 40 -o out/
    python np_funnel.py --selftest
"""
import argparse
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def read_csv(path, delim=","):
    with open(path, encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f, delimiter=delim)
        return list(rdr), (rdr.fieldnames or [])


def compound_degree(edges_path):
    """从成分-靶点边表算每个成分的靶点数（degree）。"""
    rows, hdrs = read_csv(edges_path, "\t" if edges_path.lower().endswith(".tsv") else ",")
    sc = next((h for h in hdrs if h.lower() in ("compound", "src_id", "成分", "name")), hdrs[0] if hdrs else None)
    dc = next((h for h in hdrs if h.lower() in ("gene", "dst_id", "target", "靶点")), hdrs[1] if len(hdrs) > 1 else None)
    deg = {}
    for r in rows:
        c = (r.get(sc) or "").strip()
        t = (r.get(dc) or "").strip().upper()
        if c and t:
            deg.setdefault(c.lower(), {"name": c, "targets": set()})["targets"].add(t)
    return {k: {"name": v["name"], "degree": len(v["targets"])} for k, v in deg.items()}


def target_promiscuity(edges_path):
    """反向统计：每个靶点连接的成分数（multiplicity）。
    这是 Degree 伪指标的核心——槲皮素等 promiscuous 成分会把靶点度撑高，
    但这些靶点多为非特异结合（假阳性 hub）。"""
    rows, hdrs = read_csv(edges_path, "\t" if edges_path.lower().endswith(".tsv") else ",")
    sc = next((h for h in hdrs if h.lower() in ("compound", "src_id", "成分", "name")), hdrs[0] if hdrs else None)
    dc = next((h for h in hdrs if h.lower() in ("gene", "dst_id", "target", "靶点")), hdrs[1] if len(hdrs) > 1 else None)
    mult = {}
    for r in rows:
        c = (r.get(sc) or "").strip().lower()
        t = (r.get(dc) or "").strip().upper()
        if c and t:
            mult.setdefault(t, set()).add(c)
    return {t: len(cs) for t, cs in mult.items()}


def load_disease_context(path):
    """读 CTD 疾病-基因关联（np_goldfetch 产出 gold_pd_genes.tsv），返回 {gene: score}。
    CTD 的 score/relevance 列若存在则用之，否则用出现次数作代理（均一化为 0~1）。
    """
    if not path or not os.path.exists(path):
        return {}
    raw = {}
    with open(path, encoding="utf-8-sig") as f:
        first = f.readline()
        delim = "\t" if "\t" in first else ","
        f.seek(0)
        rdr = csv.DictReader(f, delimiter=delim)
        hdrs = rdr.fieldnames or []
        gcol = next((h for h in hdrs if "gene" in h.lower() or "symbol" in h.lower()), hdrs[0])
        scol = next((h for h in hdrs if "score" in h.lower() or "relevance" in h.lower()), None)
        for r in rdr:
            g = (r.get(gcol) or "").strip().upper()
            if not g:
                continue
            if scol:
                try:
                    v = float(r[scol])
                    raw[g] = max(raw.get(g, 0), v)
                    continue
                except ValueError:
                    pass
            raw[g] = raw.get(g, 0) + 1
    if not raw:
        return {}
    mx = max(raw.values())
    return {g: v / mx for g, v in raw.items()}  # 归一 0~1


def promiscuity_penalty(mult, threshold=200, pct=None, penalize=True, disease_context=None,
                         exemption_cut=0.5):
    """多效性惩罚 + 疾病语境豁免（Disease Context Exemption）。

    疾病语境豁免（避免误杀 PTGS2 这类"高连接但 disease-relevant"靶点）：
        d_t    = 该靶点在目标疾病库（CTD）中的归一化关联强度（0~1）
        exempt = max(0, d_t - cut) / (1 - cut)         cut 默认 0.5
        w_t    = w_base + (1 - w_base) * exempt
      → d_t 越高，权重越接近 1.0（完全豁免）；d_t ≤ cut 则保留原惩罚。
      豁免的靶点在 flag 中标注 disease_exempt，不掩盖其 nonspecific 属性。
    """
    if not mult:
        return {}, {}, {}
    vals = sorted(mult.values())
    n = len(vals)
    median = vals[n // 2]
    qcut = None
    if pct is not None and 0 < pct < 1:
        qcut = vals[min(n - 1, int(round(pct * (n - 1))))]
    dc = disease_context or {}
    flags, weights = {}, {}
    for t, m in mult.items():
        is_ns = m > threshold or (qcut is not None and m >= qcut)
        factor = 1.0 / (1.0 + math.log2(max(1.0, m / max(1, median)))) if (penalize and m > median) else 1.0
        d_t = dc.get(t, 0.0)
        exempt = max(0.0, d_t - exemption_cut) / max(1e-9, 1 - exemption_cut) if dc else 0.0
        w_final = factor + (1.0 - factor) * exempt
        if is_ns:
            flags[t] = "disease_exempt" if exempt > 0 else "nonspecific_hub"
        weights[t] = round(w_final, 4)
    stats = {"n_targets": n, "median_multiplicity": median, "max_multiplicity": max(vals),
             "threshold": threshold, "pct": pct, "qcut": qcut,
             "n_nonspecific": sum(1 for v in flags.values() if v == "nonspecific_hub"),
             "n_disease_exempt": sum(1 for v in flags.values() if v == "disease_exempt"),
             "disease_context_enabled": bool(dc), "exemption_cut": exemption_cut}
    return weights, flags, stats


def promiscuity_penalty_legacy(mult, threshold=200, pct=None, penalize=True):
    """多效性惩罚（Promiscuity Penalty）。

    算法（透明、可写进 Methods）：
        对靶点 t，multiplicity m_t = 连接的不同成分数。
        判定 nonspecific_hub：
          (a) 绝对阈值：m_t > --promiscuity-threshold（默认 200）
          (b) 相对分位：m_t >= 分布的 --promiscuity-pct 分位（如 0.95），取两者并集
        惩罚权重（对数衰减，保留靶点但降级）：
          w_t = 1 / (1 + log2(m_t / m_median))  当 m_t > m_median，否则 w_t = 1.0
        即：连接成分数是中位数的 2 倍→权重 1/2，4 倍→1/3，8 倍→1/4……
        输出 (penalty_weights, flags, stats)，不删除任何靶点（诚实：降权≠造假）。
    """
    if not mult:
        return {}, {}, {}
    vals = sorted(mult.values())
    n = len(vals)
    median = vals[n // 2]
    qcut = None
    if pct is not None and 0 < pct < 1:
        idx = min(n - 1, int(round(pct * (n - 1))))
        qcut = vals[idx]
    flags, weights = {}, {}
    for t, m in mult.items():
        is_ns = m > threshold or (qcut is not None and m >= qcut)
        if is_ns:
            flags[t] = "nonspecific_hub"
        factor = 1.0 / (1.0 + math.log2(max(1.0, m / max(1, median)))) if (penalize and m > median) else 1.0
        weights[t] = round(factor, 4)
    stats = {"n_targets": n, "median_multiplicity": median, "max_multiplicity": max(vals),
             "threshold": threshold, "pct": pct, "qcut": qcut,
             "n_nonspecific": len(flags)}
    return weights, flags, stats


def load_admet(path):
    if not path or not os.path.exists(path):
        return {}
    rows, hdrs = read_csv(path, ",")
    d = {}
    bbb = next((h for h in hdrs if "bbb" in h.lower() or "血脑" in h), None)
    tox = next((h for h in hdrs if "tox" in h.lower() or "毒性" in h), None)
    for r in rows:
        c = (r.get("compound") or r.get("name") or r.get("成分") or "").strip().lower()
        if not c:
            continue
        d[c] = {"bbb": (r.get(bbb) or "") if bbb else "",
                "tox": (r.get(tox) or "") if tox else ""}
    return d


def barrier_filter(cur, admet, cns_disease=False, target_tissue=None, whitelist=None):
    """生理屏障过滤（PK/PD 空间桥接）——决定工具能否落地的关键。

    CNS 疾病：BBB=low/no/false/0/低 的成分 → 强制剔除（物理不可达脑内靶点），
              除非在白名单中（文献报道的转运体介导入脑成分可豁免）。
    其他组织：仅提示，不剔除（如肝/肾分布通常不受 BBB 限制）。

    返回 (kept, trace, unreachable)
    """
    whitelist = whitelist or {}
    BBB_NEG = {"no", "false", "0", "low", "低", "n", "negative", "bbb-", "-"}
    trace, unreachable, kept = [], [], {}
    for k, v in cur.items():
        a = admet.get(k, {})
        bbb = str(a.get("bbb", "")).strip().lower()
        is_low = bbb in BBB_NEG
        if cns_disease and is_low and k not in whitelist:
            unreachable.append({"compound": v["name"], "reason": "BBB=%s，物理不可达脑内靶点" % bbb,
                                "label": "Physically Unreachable"})
            continue
        kept[k] = v
    if cns_disease:
        trace.append({"stage": "barrier_BBB", "n": len(kept),
                      "rule": "CNS 疾病：强制剔除 BBB 阴性成分（白名单豁免 %d 个）" %
                              sum(1 for k in whitelist if k in cur)})
    elif target_tissue:
        trace.append({"stage": "barrier_hint", "n": len(kept),
                      "rule": "目标组织「%s」：非 BBB 限制，仅提示膜通透性需核查" % target_tissue})
    return kept, trace, unreachable


def simple_scaffold(name):
    """极简骨架近似：去取代基/糖苷/盐，取母核字符串（够用于粗去重）。"""
    s = name.lower()
    for token in ("-o-beta-d-glucoside", "-glucoside", "-glycoside", "sulfate", "chloride", "acetate"):
        s = s.replace(token, "")
    s = "".join(ch for ch in s if ch.isalnum())
    return s[:12]  # 母核前 12 字符作近似键


def load_whitelist(path):
    """白名单：文献报道的已知活性成分，即使低 degree/低暴露也强制进入对接。"""
    if not path or not os.path.exists(path):
        return {}
    out = {}
    with open(path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            c = (r.get("compound") or r.get("name") or r.get("成分") or "").strip()
            if c:
                out[c.lower()] = {"reason": r.get("reason") or r.get("文献") or "whitelist",
                                  "pmid": r.get("pmid") or ""}
    return out


def run_funnel(degree, admet, dose, top_ratio, min_exposure, max_final, keep_bbb_only, whitelist=None,
               cns_disease=False, target_tissue=None):
    whitelist = whitelist or {}
    trace = []
    # Stage 0
    cur = dict(degree)
    trace.append({"stage": "0_all", "n": len(cur), "rule": "全量成分"})
    # Stage 1 拓扑
    if cur:
        vals = sorted((v["degree"] for v in cur.values()), reverse=True)
        k = max(1, int(len(vals) * top_ratio))
        thr = vals[min(k - 1, len(vals) - 1)]
        cur = {k2: v for k2, v in cur.items() if v["degree"] >= thr}
        trace.append({"stage": "1_topology", "n": len(cur), "rule": "degree≥%d (top %.0f%%)" % (thr, top_ratio * 100)})
    # Stage 2 暴露量（若提供剂量）
    if dose and cur:
        kept = {}
        for k2, v in cur.items():
            pd = dose.get(k2)
            if pd is None:
                kept[k2] = v  # 无剂量信息不淘汰，但标注
                continue
            intake = pd.get("daily_g", 0) * 1000 * pd.get("content_pct", 0) / 100
            cmax = intake * (pd.get("f_pct", 30) / 100) / max(1e-9, pd.get("mw", 300)) * pd.get("fu", 0.1)
            if cmax >= min_exposure:
                kept[k2] = v
        cur = kept
        trace.append({"stage": "2_exposure", "n": len(cur), "rule": "free Cmax≥%.4f µM 或无剂量数据" % min_exposure})
    # Stage 3 ADMET
    if admet and cur:
        kept = {}
        for k2, v in cur.items():
            a = admet.get(k2, {})
            tox = a.get("tox", "").lower()
            if any(x in tox for x in ("high", "致命", "strong", "danger")):
                continue
            if keep_bbb_only and a.get("bbb", "").lower() in ("no", "false", "0", "低", "low"):
                continue
            kept[k2] = v
        cur = kept
        trace.append({"stage": "3_admet", "n": len(cur), "rule": "低毒 + BBB=%s" % ("要求可透过" if keep_bbb_only else "不限")})
    # 生理屏障过滤（CNS 疾病强制 BBB，其他组织仅提示）
    cur, barrier_trace, unreachable = barrier_filter(cur, admet, cns_disease, target_tissue, whitelist)
    trace.extend(barrier_trace)
    if unreachable:
        trace.append({"stage": "physically_unreachable", "n": len(unreachable),
                      "rule": "BBB 阴性成分（物理不可达）：%s" % ",".join(x["compound"] for x in unreachable[:8]),
                      "details": unreachable})
    # 白名单豁免：过滤后强制补回（保护"微量但高活性"成分不被误杀）
    if whitelist:
        rescued = {k: v for k, v in degree.items() if k in whitelist and k not in cur}
        cur.update(rescued)
        if rescued:
            trace.append({"stage": "whitelist_rescue", "n": len(rescued),
                          "rule": "白名单豁免补回：%s" % ",".join(sorted(rescued)[:8])})
    # Stage 4 骨架去重 + Top N（白名单优先，确保不被 max_final 截断）
    seen, dedup = set(), {}
    for k2, v in sorted(cur.items(), key=lambda x: (0 if x[0] in whitelist else 1, -x[1]["degree"])):
        sc = simple_scaffold(v["name"])
        if sc in seen and k2 not in whitelist:
            continue
        seen.add(sc)
        dedup[k2] = v
    trace.append({"stage": "4_scaffold", "n": len(dedup), "rule": "骨架去重（每骨架留度值最高；白名单豁免）"})
    final = sorted(dedup.values(),
                   key=lambda x: (0 if x["name"].lower() in whitelist else 1, -x["degree"]))[:max_final]
    n_wl_final = sum(1 for v in final if v["name"].lower() in whitelist)
    trace.append({"stage": "5_final", "n": len(final),
                  "rule": "Top %d 进入对接（含 %d 白名单豁免）" % (max_final, n_wl_final)})
    return final, trace


def selftest():
    deg = {"c%d" % i: {"name": "compound_%d" % i, "degree": 100 - i} for i in range(20)}
    final, trace = run_funnel(deg, {}, {}, 0.3, 0.01, 5, False)
    ok = len(final) <= 5 and trace[0]["n"] == 20 and trace[-1]["n"] <= 5
    # 多效性惩罚自测：中位=10，一个靶点连 200 成分应被标记且权重<1
    mult = {"T%d" % i: 10 for i in range(20)}
    mult["SUPERHUB"] = 200
    pw, flags, st = promiscuity_penalty(mult, threshold=150)
    ok2 = ("SUPERHUB" in flags and pw["SUPERHUB"] < 1.0 and st["n_nonspecific"] == 1)
    # 疾病语境豁免：同一 SUPERHUB 若强关联疾病(d=1.0)，权重应恢复到接近 1 且标记 disease_exempt
    pw2, fl2, st2 = promiscuity_penalty(mult, threshold=150,
                                        disease_context={"SUPERHUB": 1.0}, exemption_cut=0.5)
    ok3 = (fl2.get("SUPERHUB") == "disease_exempt" and pw2["SUPERHUB"] > pw["SUPERHUB"]
           and st2["n_disease_exempt"] == 1)
    # BBB 屏障：CNS 模式剔除 low，非 CNS 保留
    cur = {"q1": {"name": "q1", "degree": 10}, "q2": {"name": "q2", "degree": 8}}
    adm = {"q1": {"bbb": "low"}, "q2": {"bbb": "high"}}
    kept_cns, _, un = barrier_filter(cur, adm, cns_disease=True)
    kept_non, _, _ = barrier_filter(cur, adm, cns_disease=False)
    ok4 = (len(kept_cns) == 1 and "q2" in kept_cns and len(un) == 1
           and len(kept_non) == 2 and un[0]["label"] == "Physically Unreachable")
    print("selftest: %s (all=%d→%d; promiscuity %s→exempt %s; CNS BBB kept %d, unreachable %d)" %
          ("PASS" if ok and ok2 and ok3 and ok4 else "FAIL", trace[0]["n"], len(final),
           flags.get("SUPERHUB", "-"), fl2.get("SUPERHUB", "-"), len(kept_cns), len(un)))
    return 0 if ok and ok2 and ok3 and ok4 else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="np_funnel: 先筛选后对接漏斗")
    ap.add_argument("--edges", help="成分-靶点边表")
    ap.add_argument("--admet", default=None, help="ADMET 表 CSV")
    ap.add_argument("--dose", default=None, help="剂量表 CSV")
    ap.add_argument("--top-ratio", type=float, default=0.3)
    ap.add_argument("--min-exposure", type=float, default=0.01)
    ap.add_argument("--max-final", type=int, default=40)
    ap.add_argument("--whitelist", default=None,
                    help="白名单 CSV（compound,reason,pmid）——文献报道的微量高活性成分，强制进入对接")
    ap.add_argument("--promiscuity-threshold", type=int, default=200,
                    help="靶点连接成分数 > 此值标记为 nonspecific_hub（默认 200）")
    ap.add_argument("--promiscuity-pct", type=float, default=None,
                    help="可选：超过此分位也标记（如 0.95）")
    ap.add_argument("--no-promiscuity-penalty", action="store_true",
                    help="关闭对数降权（仅标记不降权）")
    ap.add_argument("--disease-context", default=None,
                    help="CTD 疾病关联文件（np_goldfetch 产出）——疾病语境豁免，强相关 hub 减免惩罚")
    ap.add_argument("--exemption-cut", type=float, default=0.5,
                    help="疾病关联强度豁免阈值（默认 0.5，越高越难豁免）")
    ap.add_argument("--cns-disease", action="store_true",
                    help="中枢神经系统疾病：强制剔除 BBB 阴性成分（物理不可达脑内靶点）")
    ap.add_argument("--target-tissue", default=None,
                    help="目标组织（如 Liver/Lung）——非 CNS 时仅提示膜通透性，不强制过滤")
    ap.add_argument("--keep-bbb-only", action="store_true", help="仅保留 BBB 可透过（中枢疾病）")
    ap.add_argument("-o", "--outdir", default="np_out/funnel")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if not a.edges:
        ap.error("需 --edges 或 --selftest")
    degree = compound_degree(a.edges)
    admet = load_admet(a.admet)
    dose = {}
    if a.dose and os.path.exists(a.dose):
        rows, _ = read_csv(a.dose, ",")
        for r in rows:
            c = (r.get("compound") or "").strip().lower()
            if c:
                dose[c] = {"daily_g": float(r.get("daily_g") or 0),
                           "content_pct": float(r.get("content_pct") or 0),
                           "f_pct": float(r.get("f_pct") or 30),
                           "fu": float(r.get("fu") or 0.1),
                           "mw": float(r.get("mw") or 300)}
    os.makedirs(a.outdir, exist_ok=True)
    whitelist = load_whitelist(a.whitelist)
    if whitelist:
        print("[funnel] 白名单 %d 味（强制保留）：%s" % (len(whitelist), list(whitelist)[:6]))
    final, trace = run_funnel(degree, admet, dose, a.top_ratio, a.min_exposure,
                              a.max_final, a.keep_bbb_only, whitelist,
                              cns_disease=a.cns_disease, target_tissue=a.target_tissue)
    # 写白名单命中记录
    if whitelist:
        wl_hit = [v["name"] for v in final if v["name"].lower() in whitelist]
        with open(os.path.join(a.outdir, "whitelist_trace.json"), "w", encoding="utf-8") as f:
            json.dump({"whitelist": whitelist, "rescued_to_final": wl_hit}, f,
                      ensure_ascii=False, indent=1)
    with open(os.path.join(a.outdir, "final_docking_list.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f); w.writerow(["compound", "degree"])
        for v in final:
            w.writerow([v["name"], v["degree"]])
    json.dump(trace, open(os.path.join(a.outdir, "funnel_trace.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[funnel] " + " → ".join("%s:%d" % (t["stage"], t["n"]) for t in trace))
    print("[funnel] → %s" % os.path.join(a.outdir, "final_docking_list.csv"))

    # 多效性惩罚：靶点优先级表（纠偏 Degree 伪指标）
    mult = target_promiscuity(a.edges)
    dc = load_disease_context(a.disease_context)
    if dc:
        print("[funnel] 疾病语境豁免：载入 %d 个疾病关联靶点（cut=%.2f）" % (len(dc), a.exemption_cut))
    pw, flags, pstats = promiscuity_penalty(mult, a.promiscuity_threshold, a.promiscuity_pct,
                                            penalize=not a.no_promiscuity_penalty,
                                            disease_context=dc, exemption_cut=a.exemption_cut)
    prio = sorted(mult.items(), key=lambda x: -x[1])
    with open(os.path.join(a.outdir, "target_priority.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["target", "n_compounds", "penalty_weight", "flag", "priority_rank"])
        for i, (t, m) in enumerate(prio, 1):
            w.writerow([t, m, pw.get(t, 1.0), flags.get(t, ""), i])
    pstats["note"] = ("多效性惩罚 w=1/(1+log2(m/m_median))；nonspecific_hub 为高连接假阳性候选，"
                      "对接优先级须结合 penalty_weight 排序，不得仅按 multiplicity")
    json.dump(pstats, open(os.path.join(a.outdir, "promiscuity_stats.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[funnel] 多效性：%d 靶点，中位连接 %d，标记 nonspecific_hub %d 个 → target_priority.csv" %
          (pstats["n_targets"], pstats["median_multiplicity"], pstats["n_nonspecific"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
