#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_synergy.py —— 多成分协同/拮抗分析（揭示中药配伍的现代科学内涵）

问题（中药复方科学性盲区）
    单体预测回答不了：A+B  together 是 1+1>2（协同）、1+1<2（拮抗）还是无关？
    网络拓扑图再漂亮也无法回答"这两味药加一起是好是坏"。

算法（透明、可写进 Methods）
    成分 i 的靶点集 T_i；两成分 i,j：
    1) Jaccard 重叠  J = |Ti∩Tj| / |Ti∪Tj|
    2) 网络邻近度    N = 1/(1+平均最短路径)  —— 用 PPI 上 Ti 与 Tj 的距离衡量互补/趋同
    3) 互补性        C = |Ti Δ Tj| / |Ti∪Tj|  （对称差比例，高=靶点互补）
    4) 协同评分      S = α·C + β·N + γ·(共享通路数归一)
       高 S 且共享通路≥1 → Synergy pair（不同成分打不同靶点但汇聚同一通路）
    5) 拮抗判定      若提供 --direction 表（靶点 up/down 注释）：
       共享靶点中方向相反比例 > 0.5 → Antagonism pair（A 激活 B 抑制同一靶点）
       无方向表时，仅基于"高 Jaccard + 已知 antagonist 库"标注 potential_antagonism（诚实弱化）

输出
    synergy_pairs.csv / antagonism_pairs.csv / synergy_subnetwork.svg / synergy_manifest.json

用法
    python np_synergy.py analyze --edges compound_targets.tsv --ppi string_edges.tsv \
        --pathways enrich.csv -o out/
    python np_synergy.py analyze --edges ct.tsv --direction direction.csv -o out/
    python np_synergy.py selftest
"""
import argparse
import csv
import itertools
import json
import math
import os
import sys

try:
    import networkx as nx
except ImportError:
    nx = None

HERE = os.path.dirname(os.path.abspath(__file__))


def read_edges(path, sc=None, dc=None):
    rows = []
    with open(path, encoding="utf-8-sig") as f:
        first = f.readline()
        delim = "\t" if "\t" in first else ","
        f.seek(0)
        rdr = csv.DictReader(f, delimiter=delim)
        hdrs = rdr.fieldnames or []
        # v3.0.2: 显式指定的列名若在表头不存在，自动回退到常见别名。
        # 典型场景：np_ppi 输出 node1/node2，而调用方按旧约定传 source/target。
        def _pick(preferred, aliases, fallback):
            if preferred and preferred in hdrs:
                return preferred
            for a in aliases:
                if a in hdrs:
                    return a
            return fallback
        sc = _pick(sc, ("source", "node1", "from", "protein1", "compound", "src_id", "成分", "name"), hdrs[0] if hdrs else None)
        dc = _pick(dc, ("target", "node2", "to", "protein2", "gene", "dst_id", "target_name",
                        "symbol", "靶点"), hdrs[1] if len(hdrs) > 1 else None)
        for r in rdr:
            c = (r.get(sc) or "").strip()
            g = (r.get(dc) or "").strip().upper()
            if c and g:
                rows.append((c.lower(), g))
    return rows


def compound_targets(rows):
    d = {}
    for c, g in rows:
        d.setdefault(c, set()).add(g)
    return d


def load_pathways(path):
    """富集表 genes 列 → {term: set(genes)}。"""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            genes = set(g.strip().upper() for g in (r.get("genes") or "").replace(";", ",").split(",") if g.strip())
            if genes:
                out[r.get("term", "")] = genes
    return out


def load_direction(path):
    """靶点方向表 {gene: 'up'/'down'}。"""
    d = {}
    if not path or not os.path.exists(path):
        return d
    with open(path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            g = (r.get("gene") or r.get("target") or "").strip().upper()
            v = (r.get("direction") or r.get("change") or "").strip().lower()
            if g and v in ("up", "down", "上调", "下调"):
                d[g] = "up" if v in ("up", "上调") else "down"
    return d


_DIST_CACHE = {}


def ppi_distance(G, ta, tb, max_pairs=100):
    """两靶点集间的平均最短路径（采样+全局缓存），返回邻近度 1/(1+dist)。"""
    if G is None or not ta or not tb:
        return 0.0
    common = ta & tb
    if common:
        return 1.0
    dists = []
    # 最多取每边 8 个代表靶点，配对上限 100；缓存单对距离
    for a in list(ta)[:8]:
        for b in list(tb)[:8]:
            key = (a, b) if a <= b else (b, a)
            if key in _DIST_CACHE:
                d = _DIST_CACHE[key]
            else:
                try:
                    d = nx.shortest_path_length(G, a, b)
                except (nx.NetworkXNoPath, nx.NodeNotFound):
                    d = None
                _DIST_CACHE[key] = d
            if d is not None and d <= 6:
                dists.append(d)
            if len(dists) >= max_pairs:
                break
        if len(dists) >= max_pairs:
            break
    return 1.0 / (1.0 + sum(dists) / len(dists)) if dists else 0.0


def analyze(ct, pathways=None, direction=None, G=None, alpha=0.4, beta=0.35, gamma=0.25,
            max_compounds=60, score_cut=0.55):
    """主算法。返回 (synergy_pairs, antagonism_pairs, stats)。"""
    comps = sorted(ct)
    if len(comps) > max_compounds:
        comps = comps[:max_compounds]  # 大复方控制 O(n²) 预算
    synergy, antagonism = [], []
    for a, b in itertools.combinations(comps, 2):
        ta, tb = ct[a], ct[b]
        if not ta or not tb:
            continue
        inter = ta & tb
        union = ta | tb
        jac = len(inter) / len(union) if union else 0
        comp = len(ta ^ tb) / len(union) if union else 0
        prox = ppi_distance(G, ta, tb)
        # 共享通路：同时命中两成分靶点的通路
        shared_pw = [t for t, genes in pathways.items()
                     if (genes & ta) and (genes & tb)]
        pw_score = min(1.0, len(shared_pw) / 5.0)
        S = alpha * comp + beta * prox + gamma * pw_score
        if S >= score_cut and shared_pw:
            synergy.append({"compound_a": a, "compound_b": b, "score": round(S, 3),
                            "jaccard": round(jac, 3), "complementarity": round(comp, 3),
                            "network_proximity": round(prox, 3),
                            "shared_pathways": ";".join(shared_pw[:4]),
                            "label": "Synergy"})
        # 拮抗：共享靶点中方向表标注了相反方向（up 与 down 同时出现）。
        # v3.0.3: 移除旧的 _hint()/opp 死代码（恒返回空、结果未被使用）。
        # 基于成分自身方向的精确推断请用 --mode mechanistic + --actions；
        # 此处 topological 分支仅按共享靶点方向表的方向多样性作粗粒度提示。
        if direction:
            rev = sum(1 for g in inter if g in direction) and len(
                {direction.get(g) for g in inter if g in direction}) >= 2
            if len(inter) >= 2 and rev:
                antagonism.append({"compound_a": a, "compound_b": b,
                                   "shared_targets": ";".join(sorted(inter)[:8]),
                                   "n_shared": len(inter),
                                   "label": "Antagonism (opposite directions on shared targets)"})
    stats = {"n_compounds": len(comps), "n_pairs_evaluated": len(comps) * (len(comps) - 1) // 2,
             "n_synergy": len(synergy), "n_antagonism": len(antagonism)}
    return synergy, antagonism, stats


def draw_subnetwork(synergy, ct, out):
    """画协同成分对-共享通路子网络。"""
    if not synergy or nx is None:
        return None
    G = nx.Graph()
    # 只画 top 30 高分边，避免 spring_layout 在大图上卡住
    for p in sorted(synergy, key=lambda x: -x["score"])[:30]:
        G.add_edge(p["compound_a"], p["compound_b"], weight=p["score"])
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        pos = nx.circular_layout(G)  # 确定性布局，避免 spring_layout 在大图上耗时
        fig, ax = plt.subplots(figsize=(8, 6))
        nx.draw_networkx(G, pos, ax=ax, node_size=300, node_color="#4F8A5B",
                         font_size=6, edge_color="#999999")
        ax.set_title("Synergy pairs network")
        ax.set_axis_off()
        p = os.path.join(out, "synergy_subnetwork.svg")
        fig.savefig(p, bbox_inches="tight")
        plt.close(fig)
        return p
    except Exception:
        return None


# ==================== 机制定向协同（Mechanistic Synergy） ====================
# 为什么拓扑法不够：两成分共享同一靶点，多半是抢同一个口袋（竞争性拮抗），
# 不是协同。真正的协同要求"同一条通路的不同层级 + 干预方向不冲突"。
# 判定优先级（硬规则，写进 Methods 可复现）：
#   1) 共享靶点 → Competitive_Antagonism（同一蛋白即同一口袋；若有残基数据则按 Jaccard≥0.5 判"同一口袋"）
#   2) 同通路不同层级 + 方向冲突（一抑一激）→ Directional_Antagonism
#   3) 同通路不同层级 + 方向兼容（同为抑制/同为激活/方向未知）→ Mechanistic_Synergy
#   4) 仅拓扑相邻、无通路层级证据 → Topological_Only（不得写作协同）
RECEPTORS = {
    "EGFR", "ERBB2", "ERBB3", "ERBB4", "INSR", "IGF1R", "KDR", "FLT1", "FLT4", "FGFR1", "FGFR2",
    "PDGFRA", "PDGFRB", "MET", "ALK", "RET", "KIT", "NTRK1", "TGFBR1", "TGFBR2", "ACVR1", "BMPR1A",
    "ADRB1", "ADRB2", "ADRA1A", "AGTR1", "AVPR2", "CNR1", "CNR2", "DRD2", "HTR2A", "OPRM1",
    "CHRM1", "CHRM2", "GABBR1", "GRM5", "SSTR2", "PTGER2", "PTGER4", "CXCR4", "CCR5", "IL6R", "IL1R1",
    "TNFRSF1A", "TNFRSF1B", "TLR4", "IFNGR1", "LEPR", "PRLR", "TSHR", "FSHR", "LHCGR",
    "ESR1", "ESR2", "AR", "PGR", "NR3C1", "PPARG", "PPARA", "NR1H4", "VDR", "RARA", "THRA",
    "SLC6A4", "SLC6A2", "SCN5A", "KCNH2", "CFTR", "ABCG2",
}


def load_actions(path):
    """成分级动作表 → {compound_lower: {GENE: 'inhibit'|'activate'}}。"""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, encoding="utf-8-sig") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        rdr = csv.DictReader(f, delimiter=delim)
        low = {(h or "").lower().strip(): h for h in (rdr.fieldnames or [])}
        cc = next((low[c] for c in ("compound", "mol_id", "成分", "name") if c in low), None)
        gc = next((low[c] for c in ("gene", "target", "靶点") if c in low), None)
        ac = next((low[c] for c in ("action", "effect", "direction", "change") if c in low), None)
        if not (cc and gc and ac):
            return out
        for r in rdr:
            c = (r.get(cc) or "").strip().lower()
            g = (r.get(gc) or "").strip().upper()
            v = (r.get(ac) or "").strip().lower()
            if not (c and g):
                continue
            if v in ("inhibit", "down", "antagonist", "suppress", "抑制", "下调"):
                out.setdefault(c, {})[g] = "inhibit"
            elif v in ("activate", "up", "agonist", "induce", "激活", "上调"):
                out.setdefault(c, {})[g] = "activate"
    return out


def load_pockets(path):
    """结合位点残基表 → {GENE: set(残基)}。支持 JSON（gene→list/逗号串）或 CSV(gene,residues)。"""
    if not path or not os.path.exists(path):
        return {}
    out = {}
    if path.lower().endswith(".json"):
        try:
            j = json.load(open(path, encoding="utf-8"))
        except Exception:
            return {}
        for g, v in (j.items() if isinstance(j, dict) else []):
            if isinstance(v, (list, tuple)):
                out[str(g).upper()] = {str(x) for x in v}
            elif isinstance(v, str):
                out[str(g).upper()] = {x.strip() for x in v.replace(";", ",").split(",") if x.strip()}
        return out
    with open(path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            g = (r.get("gene") or r.get("target") or "").strip().upper()
            v = r.get("residues") or r.get("residue") or ""
            if g and v:
                out[g] = {x.strip() for x in v.replace(";", ",").split(",") if x.strip()}
    return out


def _level_map(sub, sources):
    """通路子图内各基因到"受体源集"的 BFS 距离（层级）。"""
    if nx is None or sub is None or sub.number_of_nodes() == 0:
        return {}
    src = [s for s in sources if s in sub]
    if not src:
        return {}
    return nx.multi_source_dijkstra_path_length(sub, src)


def pathway_layers(pathways, G, receptors=None):
    """每条通路 → {gene: level}。G 为 PPI 图；无 G 时按通路基因顺序给不出层级，返回 {}。"""
    src = set(receptors or RECEPTORS)
    out = {}
    for term, genes in (pathways or {}).items():
        if G is None:
            continue
        nodes = [g for g in genes if g in G]
        if len(nodes) < 3:
            continue
        sub = G.subgraph(nodes)
        lv = _level_map(sub, src)
        if lv:
            out[term] = lv
    return out


def _same_pocket(g1, g2, pockets):
    """同一基因必定同口袋；有残基数据时按 Jaccard≥0.5 判同一口袋。"""
    if g1 == g2:
        a, b = pockets.get(g1), pockets.get(g2)
        if a and b:
            return len(a & b) / max(1, len(a | b)) >= 0.5
        return True
    a, b = pockets.get(g1), pockets.get(g2)
    if a and b:
        return len(a & b) / max(1, len(a | b)) >= 0.5
    return False


def classify_pair(ta, tb, layers=None, actions=None, pockets=None, min_gap=1):
    """两成分靶点集的机制关系判定（纯函数，可单测）。

    ta/tb   : set(gene)
    layers  : {term: {gene: level}}，来自 pathway_layers()
    actions : (dict gene→action for A, dict for B)，未知则传 None
    pockets : {gene: set(residues)}
    返回 dict: label / reason / term / upstream / downstream / level_gap / direction
    """
    pockets = pockets or {}
    common = set(ta) & set(tb)
    if common:
        # 同一靶点 = 同一结合口袋（除非有残基数据证明不同位点）
        same = [g for g in common if _same_pocket(g, g, pockets)]
        if same:
            return {"label": "Competitive_Antagonism", "reason": "共享靶点 %s：争夺同一结合口袋，"
                    "属竞争性拮抗而非协同" % ",".join(sorted(same)[:4]),
                    "term": "", "upstream": "", "downstream": "", "level_gap": 0, "direction": ""}
    best = None
    for term, lv in (layers or {}).items():
        ga = [g for g in ta if g in lv]
        gb = [g for g in tb if g in lv]
        for x in ga:
            for y in gb:
                gap = abs(lv[x] - lv[y])
                if gap < min_gap:
                    continue
                up, dn = (x, y) if lv[x] < lv[y] else (y, x)
                # 取「上游节点所属成分」与「下游节点所属成分」的动作
                if lv[x] < lv[y]:
                    act_up, act_dn = (actions or ({}, {}))[0].get(x), (actions or ({}, {}))[1].get(y)
                else:
                    act_up, act_dn = (actions or ({}, {}))[1].get(y), (actions or ({}, {}))[0].get(x)
                conflict = bool(act_up and act_dn and act_up != act_dn)
                cand = {"term": term, "upstream": up, "downstream": dn, "level_gap": gap,
                        "direction": "%s/%s" % (act_up or "?", act_dn or "?"),
                        "conflict": conflict}
                if conflict:
                    cand["label"] = "Directional_Antagonism"
                    cand["reason"] = ("同一通路 %s 的不同层级（%s→%s）但干预方向相反（%s/%s）："
                                      "一个抑制上游、一个激活下游，净通量互相抵消" %
                                      (term, up, dn, act_up, act_dn))
                    return cand
                cand["label"] = "Mechanistic_Synergy"
                cand["reason"] = ("同一通路 %s 的上下游节点（%s[L%d] → %s[L%d]，间隔 %d），"
                                  "方向兼容（%s/%s）：构成串联阻断" %
                                  (term, up, lv[up], dn, lv[dn], gap, act_up or "未标注", act_dn or "未标注"))
                if best is None or gap > best["level_gap"]:
                    best = cand
    if best:
        return best
    # 无共享靶点、无通路层级证据
    return {"label": "Topological_Only", "reason": "无共享靶点，且未在任一通路内形成上下层级关系；"
            "只能作拓扑相邻，不得写作协同", "term": "", "upstream": "", "downstream": "",
            "level_gap": 0, "direction": ""}


def analyze_mechanistic(ct, pathways=None, G=None, actions=None, pockets=None,
                        receptors=None, max_compounds=60, min_gap=1):
    """机制定向协同主流程。返回 (mechanistic, antagonism, stats)。"""
    comps = sorted(ct)
    if len(comps) > max_compounds:
        comps = comps[:max_compounds]
    layers = pathway_layers(pathways or {}, G, receptors)
    mech, ant = [], []
    for a, b in itertools.combinations(comps, 2):
        ta, tb = ct[a], ct[b]
        if not ta or not tb:
            continue
        acts = (actions.get(a, {}) if actions else {}, actions.get(b, {}) if actions else {})
        r = classify_pair(ta, tb, layers, acts, pockets, min_gap)
        if r["label"] in ("Topological_Only",):
            continue
        row = {"compound_a": a, "compound_b": b, "label": r["label"], "term": r["term"],
               "upstream": r["upstream"], "downstream": r["downstream"],
               "level_gap": r["level_gap"], "direction": r["direction"],
               "jaccard": round(len(ta & tb) / len(ta | tb), 3) if (ta | tb) else 0,
               "network_proximity": round(ppi_distance(G, ta, tb), 3),
               "score": round(min(1.0, 0.45 * min(1.0, r["level_gap"] / 2.0) +
                                  0.35 * (1.0 if r["label"] == "Mechanistic_Synergy" else 0.0) +
                                  0.20 * (1.0 / (1.0 + r["level_gap"]))), 3),
               "reason": r["reason"]}
        (mech if r["label"] == "Mechanistic_Synergy" else ant).append(row)
    mech.sort(key=lambda x: -x["score"])
    stats = {"mode": "mechanistic", "n_compounds": len(comps),
             "n_pairs_evaluated": len(comps) * (len(comps) - 1) // 2,
             "n_mechanistic_synergy": len(mech), "n_antagonism": len(ant),
             "n_pathways_layered": len(layers),
             "n_receptors": len([g for g in (receptors or RECEPTORS)]),
             "caveat": ("同一靶点默认判为竞争性拮抗；若有变构位点证据（不同残基集合）应提供 --pockets 覆盖。"
                        "无 --actions 时方向视为未知，按兼容处理并在 direction 列标 ?/?" )}
    return mech, ant, stats


def selftest():
    checks = []
    # --- A. 旧拓扑算法回归（不能因为升级把它弄坏） ---
    ct = {"a": {"TP53", "CASP3", "BAX"}, "b": {"IL6", "TNF", "TP53"},
          "c": {"EGFR", "MAPK1", "AKT1"}, "d": {"TP53", "BAX"}}
    pw = {"apoptosis": {"TP53", "CASP3", "BAX", "IL6"}, "inflammation": {"IL6", "TNF"}}
    syn, ant, st = analyze(ct, pathways=pw, max_compounds=10, score_cut=0.3)
    checks.append(("topological_regression",
                   st["n_pairs_evaluated"] == 6 and len(syn) >= 1
                   and all(p["score"] >= 0.3 for p in syn), "pairs=%d syn=%d" % (st["n_pairs_evaluated"], len(syn))))

    # --- B. 机制定向协同 ---
    if nx is None:
        print("selftest: SKIP（未安装 networkx，无法构建通路层级）")
        return 0
    G = nx.Graph()
    for e in [("EGFR", "PIK3CA"), ("PIK3CA", "AKT1"), ("AKT1", "MTOR"), ("TNF", "RELA")]:
        G.add_edge(*e)
    paths = {"PI3K-Akt": {"EGFR", "PIK3CA", "AKT1", "MTOR"}}
    layers = pathway_layers(paths, G)
    checks.append(("pathway_layers", layers.get("PI3K-Akt", {}).get("EGFR") == 0
                   and layers["PI3K-Akt"].get("AKT1") == 2, str(layers.get("PI3K-Akt"))))

    r1 = classify_pair({"EGFR"}, {"AKT1"}, layers)          # 受体 vs 下游激酶 → 协同
    checks.append(("mech_synergy_upstream_downstream", r1["label"] == "Mechanistic_Synergy",
                   "%s/%s" % (r1["label"], r1.get("level_gap"))))

    r2 = classify_pair({"AKT1"}, {"AKT1", "MTOR"}, layers)  # 共享 AKT1 → 竞争性拮抗
    checks.append(("shared_target_competitive", r2["label"] == "Competitive_Antagonism", r2["label"]))

    r3 = classify_pair({"EGFR"}, {"AKT1"}, layers,
                       actions=({"EGFR": "inhibit"}, {"AKT1": "activate"}))   # 方向冲突
    checks.append(("direction_conflict", r3["label"] == "Directional_Antagonism", r3["label"]))

    r4 = classify_pair({"EGFR"}, {"AKT1"}, layers,
                       actions=({"EGFR": "inhibit"}, {"AKT1": "inhibit"}))    # 串联阻断
    checks.append(("series_block", r4["label"] == "Mechanistic_Synergy", r4["label"]))

    r5 = classify_pair({"TNF"}, {"MTOR"}, layers)           # 不同通路，无层级 → 仅拓扑
    checks.append(("no_pathway_no_synergy", r5["label"] == "Topological_Only", r5["label"]))

    # --- C. 主流程：两个成分在 PI3K-Akt 上下游 ---
    ct2 = {"tanshinone": {"EGFR"}, "ginsenoside": {"AKT1"}, "quercetin": {"AKT1", "MTOR"}}
    mech, antagon, st2 = analyze_mechanistic(ct2, paths, G, max_compounds=10)
    checks.append(("analyze_mechanistic_end_to_end",
                   st2["n_mechanistic_synergy"] >= 1 and st2["n_antagonism"] >= 1,
                   "mech=%d ant=%d" % (st2["n_mechanistic_synergy"], st2["n_antagonism"])))

    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for name, good, d in checks:
        print("  [%s] %-34s %s" % ("✓" if good else "✗", name, d))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["selftest"]:        # 兼容统一写法：selftest 与 --selftest 等价
        argv = ["--selftest"]
    ap = argparse.ArgumentParser(description="np_synergy: 多成分协同/拮抗分析")
    ap.add_argument("--edges", default=None)
    ap.add_argument("--ppi", default=None)
    ap.add_argument("--pathways", default=None)
    ap.add_argument("--direction", default=None)
    ap.add_argument("--max-compounds", type=int, default=60)
    ap.add_argument("-o", "--outdir", default="np_out/synergy")
    ap.add_argument("--draw", action="store_true", help="绘制协同子网络图（大复方建议关闭以提速）")
    ap.add_argument("--score-cut", type=float, default=0.55, help="协同判定阈值（默认 0.55）")
    ap.add_argument("--mode", default="auto", choices=["auto", "mechanistic", "topological"],
                    help="mechanistic=通路上下游+方向兼容；topological=旧版拓扑评分")
    ap.add_argument("--actions", default=None, help="成分级动作表 compound,gene,action（mechanistic 模式推荐）")
    ap.add_argument("--pockets", default=None, help="结合位点残基表 gene,residues（判定同一口袋）")
    ap.add_argument("--receptors", default=None, help="受体基因列表文件（每行一个；默认内置常见受体集）")
    ap.add_argument("--min-gap", type=int, default=1, help="上下游最小层级差（默认 1）")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if not a.edges:
        ap.error("需 --edges 或 --selftest")
    rows = read_edges(a.edges)
    ct = compound_targets(rows)
    pathways = load_pathways(a.pathways)
    direction = load_direction(a.direction)
    G = None
    if a.ppi and nx and os.path.exists(a.ppi):
        G = nx.Graph()
        for r in read_edges(a.ppi, "source", "target"):
            G.add_edge(*r)
        print("[synergy] PPI: %d 节点 %d 边" % (G.number_of_nodes(), G.number_of_edges()))
    actions = load_actions(a.actions)
    pockets = load_pockets(a.pockets)
    receptors = None
    if a.receptors and os.path.exists(a.receptors):
        receptors = {l.strip().upper() for l in open(a.receptors, encoding="utf-8-sig") if l.strip()}
    mode = a.mode
    if mode == "auto":
        mode = "mechanistic" if (G is not None and pathways) else "topological"
    if mode == "mechanistic":
        syn, ant, st = analyze_mechanistic(ct, pathways, G, actions, pockets, receptors,
                                          max_compounds=a.max_compounds, min_gap=a.min_gap)
        fn_syn, cols_syn = "mechanistic_synergy.csv", ["compound_a", "compound_b", "label", "score",
                                                       "term", "upstream", "downstream", "level_gap",
                                                       "direction", "jaccard", "network_proximity", "reason"]
        fn_ant, cols_ant = "antagonism_pairs.csv", ["compound_a", "compound_b", "label", "score", "term",
                                                    "upstream", "downstream", "direction", "jaccard", "reason"]
    else:
        syn, ant, st = analyze(ct, pathways, direction, G,
                               max_compounds=a.max_compounds, score_cut=a.score_cut)
        fn_syn, cols_syn = "synergy_pairs.csv", ["compound_a", "compound_b", "score", "jaccard",
                                                 "complementarity", "network_proximity", "shared_pathways", "label"]
        fn_ant, cols_ant = "antagonism_pairs.csv", ["compound_a", "compound_b", "shared_targets",
                                                    "n_shared", "label"]
    os.makedirs(a.outdir, exist_ok=True)
    for fn, data, cols in [(fn_syn, syn, cols_syn), (fn_ant, ant, cols_ant)]:
        with open(os.path.join(a.outdir, fn), "w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader(); w.writerows(data)
    fig = draw_subnetwork(syn, ct, a.outdir) if a.draw else None
    st["figure"] = os.path.basename(fig) if fig else "not_drawn (use --draw)"
    st["mode"] = mode
    st["note"] = ("mechanistic：同一通路上下游 + 方向兼容才判协同；共享靶点强制 Competitive_Antagonism"
                  if mode == "mechanistic" else
                  "topological：互补靶点汇聚同一通路；共享靶点方向相反判拮抗")
    json.dump(st, open(os.path.join(a.outdir, "synergy_manifest.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    if mode == "mechanistic":
        print("[synergy] %s：%d 成分 → %d 对：机制协同 %d / 拮抗 %d（通路分层 %d 条）→ %s" %
              (mode, st["n_compounds"], st["n_pairs_evaluated"], st["n_mechanistic_synergy"],
               st["n_antagonism"], st["n_pathways_layered"], a.outdir))
    else:
        print("[synergy] %s：%d 成分 → %d 对：协同 %d / 拮抗 %d → %s" %
              (mode, st["n_compounds"], st["n_pairs_evaluated"], st["n_synergy"], st["n_antagonism"], a.outdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
