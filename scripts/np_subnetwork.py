#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_subnetwork.py —— 巨图拆解：Louvain 聚类 → 生物学子网络 → 分图出图（解决 Hairball）

问题
    PPI/成分-靶点网络超过 500 节点就是一团乱麻（hairball），无法用于论文。
    直接整图出 SVG 审稿人看不清，必须拆成有生物学意义的子网络。

方案
    Louvain 社区发现（networkx.algorithms.community.louvain_communities）把大图拆成模块；
    每个模块 ≥ min_size 节点者单独出一张 SVG（含节点标签、度值着色）；
    同时输出一张"总览图"（社区用颜色区分，不标节点名，只看模块结构）；
    附 enrichment-per-community：对每个社区跑一次轻量富集（若提供 GMT），给社区起生物学名字。

诚实边界
    * Louvain 有随机性 → 固定 seed（默认 42），meta 记录 seed 与 version。
    * 社区数由算法决定，不强行凑 3-4 个；若结果过多，用 --top-communities 取最大的几个。
    * 子网络是"拓扑模块"，不等同于"生物学通路"，需富集注释后才可命名。

用法
    python np_subnetwork.py --edges string_edges.tsv --outdir figs/subnetworks --seed 42
    python np_subnetwork.py --edges edges.tsv --top-communities 4 --min-size 10
    python np_subnetwork.py --selftest
"""
import argparse
import json
import os
import sys

try:
    import networkx as nx
except ImportError:
    nx = None

HERE = os.path.dirname(os.path.abspath(__file__))


def load_edges(path, min_score=0.0):
    rows = []
    with open(path, encoding="utf-8-sig") as f:
        first = f.readline()
        delim = "\t" if "\t" in first else ","
        f.seek(0)
        import csv
        rdr = csv.DictReader(f, delimiter=delim)
        hdrs = rdr.fieldnames or []
        a = next((h for h in hdrs if h.lower() in ("preferredname_a", "node1", "source", "src_id")), hdrs[0])
        b = next((h for h in hdrs if h.lower() in ("preferredname_b", "node2", "target", "dst_id")), hdrs[1] if len(hdrs) > 1 else None)
        sc = next((h for h in hdrs if h.lower() in ("score", "combined_score", "weight")), None)
        for r in rdr:
            u, v = (r.get(a) or "").strip(), (r.get(b) or "").strip()
            if not u or not v or u == v:
                continue
            s = 1.0
            if sc:
                try:
                    s = float(r.get(sc) or 1)
                except ValueError:
                    s = 1.0
            if s >= min_score:
                rows.append((u, v, s))
    return rows


def load_significant_pathways(enrich_path, adjp_cut=0.05):
    """载入显著富集通路 {term: set(genes)}（来自 np_enrich 的富集_*_三库.csv）。"""
    import csv
    out = {}
    if not enrich_path or not os.path.exists(enrich_path):
        return out
    with open(enrich_path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            try:
                adjp = float(r.get("adjp") or r.get("adjP") or 1.0)
            except ValueError:
                continue
            if adjp > adjp_cut:
                continue
            genes = set(g.strip().upper() for g in (r.get("genes") or "").replace(";", ",").split(",") if g.strip())
            if genes:
                out[r.get("term", "")] = genes
    return out


def hypergeom_p(k, M, N, n):
    """社区 gene set 与通路 gene set 重叠的超几何尾概率 P(X>=k)。纯 math，零依赖。"""
    import math
    if k <= 0:
        return 1.0
    total = math.comb(M, n)
    s = 0.0
    for i in range(k, min(n, N) + 1):
        if N - i >= 0 and M - N >= n - i:
            s += math.comb(N, i) * math.comb(M - N, n - i)
    return min(1.0, s / total) if total else 1.0


def annotate_community(comm_genes, pathways, universe_n, adjp_cut=0.05):
    """对一个社区做富集注释，返回 (best_term, best_p, n_sig)。"""
    comm = set(comm_genes)
    n = len(comm)
    best_term, best_p, n_sig = "", 1.0, 0
    for term, pgenes in pathways.items():
        k = len(comm & pgenes)
        if k == 0:
            continue
        p = hypergeom_p(k, universe_n, len(pgenes), n)
        if p < adjp_cut:
            n_sig += 1
            if p < best_p:
                best_p, best_term = p, term
    return best_term, best_p, n_sig


def communities(G, seed=42, resolution=1.0):
    if hasattr(nx.algorithms.community, "louvain_communities"):
        return nx.algorithms.community.louvain_communities(G, seed=seed, resolution=resolution)
    # 回退：greedy modularity
    from networkx.algorithms.community import greedy_modularity_communities
    return greedy_modularity_communities(G)


def draw_subgraph(H, path, title, color):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    deg = dict(H.degree())
    pos = nx.spring_layout(H, seed=42, k=0.5, iterations=200)
    fig = plt.figure(figsize=(7, 5.2))
    ax = fig.add_axes([0.05, 0.08, 0.9, 0.82]); ax.set_axis_off()
    nx.draw_networkx_edges(H, pos, ax=ax, alpha=0.3, edge_color="#B0B8C0", width=0.6)
    sizes = [60 + 12 * deg.get(n, 0) for n in H]
    nx.draw_networkx_nodes(H, pos, ax=ax, node_color=color, node_size=sizes,
                           alpha=0.85, edgecolors="white", linewidths=0.5)
    labels = {n: n for n in sorted(H, key=lambda x: -deg.get(x, 0))[:25]}
    nx.draw_networkx_labels(H, pos, labels, ax=ax, font_size=6)
    ax.set_title(title, fontsize=9, fontweight="bold", loc="left")
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description="np_subnetwork: Louvain 巨图拆解与子网络出图")
    ap.add_argument("--edges", help="边表 TSV/CSV")
    ap.add_argument("--outdir", default="figs/subnetworks")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--min-size", type=int, default=5)
    ap.add_argument("--top-communities", type=int, default=None)
    ap.add_argument("--min-score", type=float, default=0.0)
    ap.add_argument("--enrich", default=None,
                    help="富集结果 CSV（含 genes/adjp 列）；启用 KEGG 显著性过滤，无显著通路的社区降级合并")
    ap.add_argument("--adjp-cut", type=float, default=0.05, help="社区富集显著性阈值")
    ap.add_argument("--keep-orphan", action="store_true",
                    help="无显著通路的社区不合并，单独作 Unknown module 出图（保留原拓扑）")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        if nx is None:
            print("selftest: SKIP (需 networkx)"); return 0
        G = nx.karate_club_graph()
        cs = communities(G, seed=42)
        ok = len(cs) >= 2
        print("selftest: %s (karate → %d 社区)" % ("PASS" if ok else "FAIL", len(cs)))
        return 0 if ok else 1
    if nx is None:
        raise SystemExit("需 networkx：pip install networkx")
    if not a.edges:
        ap.error("需 --edges 或 --selftest")
    rows = load_edges(a.edges, a.min_score)
    G = nx.Graph()
    for u, v, s in rows:
        G.add_edge(u, v, weight=s)
    cs = communities(G, seed=a.seed, resolution=a.resolution)
    cs = sorted(cs, key=len, reverse=True)
    os.makedirs(a.outdir, exist_ok=True)

    # ---- KEGG 显著性过滤（堵住"强行找关联"）----
    pathways = load_significant_pathways(a.enrich, a.adjp_cut)
    annot = {}
    if pathways:
        print("[subnetwork] 载入 %d 条显著通路（adjP<%.3f）" % (len(pathways), a.adjp_cut))
        for idx, comm in enumerate(cs):
            term, p, nsig = annotate_community(comm, pathways, G.number_of_nodes(), a.adjp_cut)
            annot[idx] = {"term": term, "p": p, "n_sig": nsig, "pass": nsig > 0}
        n_pass = sum(1 for v in annot.values() if v["pass"])
        print("[subnetwork] 社区富集过滤：%d/%d 有显著通路" % (n_pass, len(cs)))
        # 降级合并：无显著通路的小社区并入有显著通路且重叠边最多的社区
        keep_idx = [i for i in range(len(cs)) if annot[i]["pass"]]
        if not keep_idx:
            print("[warn] 所有社区均无显著富集通路；按拓扑保留前 %d 个，但已在 manifest 标注 biological_support=false"
                  % min(len(cs), a.top_communities or 4))
            keep_idx = list(range(min(len(cs), a.top_communities or 4)))
        drop_idx = [i for i in range(len(cs)) if i not in keep_idx]
        merged, orphans = [], []
        if a.keep_orphan:
            # 保留孤立模块：不合并，标注为 Unknown module 单独出图
            for di in drop_idx:
                annot[di]["orphan"] = True
                orphans.append(di)
            keep_idx = sorted(keep_idx + orphans)
            print("[subnetwork] --keep-orphan：%d 个无注释社区保留为 Unknown module（不合并）" % len(orphans))
        else:
            for di in drop_idx:
                # 找与 di 连边最多的一个保留社区
                best, best_n = None, -1
                for ki in keep_idx:
                    n_edges = sum(1 for u, v in G.edges() if (u in cs[di] and v in cs[ki]) or (v in cs[di] and u in cs[ki]))
                    if n_edges > best_n:
                        best, best_n = ki, n_edges
                if best is not None:
                    cs[best] = cs[best] | cs[di]  # 合并节点集
                    merged.append({"community": di, "merged_into": best, "n_nodes": len(cs[di])})
            print("[subnetwork] 降级合并 %d 个无显著通路社区" % len(merged))
    else:
        for idx in range(len(cs)):
            annot[idx] = {"term": "", "p": 1.0, "n_sig": 0, "pass": False}
        merged = []
        if a.enrich:
            print("[warn] 指定的富集文件不可用，跳过显著性过滤（社区仅拓扑依据）")

    palette = ["#2E5E8C", "#C0504D", "#4F8A5B", "#8064A2", "#E3A857", "#4BACC6", "#7F7F7F"]
    kept, manifest = [], []
    out_i = 0
    for i, comm in enumerate(cs):
        if len(comm) < a.min_size:
            continue
        # 合并模式下：被合并掉的社区不再单独出图（其内容已并入保留社区）
        if pathways and not a.keep_orphan and i in drop_idx:
            continue
        if a.top_communities and out_i >= a.top_communities:
            break
        H = G.subgraph(comm).copy()
        # 重新取该社区的富集注释（合并不影响 index 的用原值）
        ann = annot.get(i, {"term": "", "p": 1.0, "n_sig": 0, "pass": False})
        fn = os.path.join(a.outdir, "subnetwork_%d.svg" % (out_i + 1))
        if ann.get("orphan"):
            title = "Unknown module %d (%d nodes, no significant pathway)" % (out_i + 1, len(H))
        else:
            title = "Subnetwork %d (%d nodes)%s" % (
                out_i + 1, len(H),
                " | %s (P=%.1e)" % (ann["term"], ann["p"]) if ann["term"] else " | no pathway annotation")
        draw_subgraph(H, fn, title, palette[out_i % len(palette)])
        top = sorted(H.degree(), key=lambda x: -x[1])[:5]
        manifest.append({"id": out_i + 1, "n_nodes": len(H), "n_edges": H.number_of_edges(),
                         "file": os.path.basename(fn), "top_degree": [n for n, _ in top],
                         "enriched_term": ann["term"], "enrich_p": round(ann["p"], 6),
                         "n_sig_pathways": ann["n_sig"], "biological_support": ann["pass"]})
        kept.append(fn)
        out_i += 1
    # 总览图（社区着色，不标节点名）
    overview = os.path.join(a.outdir, "overview.svg")
    try:
        node2c = {}
        for i, comm in enumerate(cs):
            for n in comm:
                node2c[n] = i
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        pos = nx.spring_layout(G, seed=42, iterations=200)
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_axes([0.03, 0.03, 0.94, 0.94]); ax.set_axis_off()
        nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.15, edge_color="#CCCCCC", width=0.4)
        cols = [palette[node2c.get(n, 0) % len(palette)] for n in G]
        nx.draw_networkx_nodes(G, pos, ax=ax, node_color=cols, node_size=18, alpha=0.8)
        ax.set_title("Community overview (%d nodes, %d communities, Louvain seed=%d)" %
                     (G.number_of_nodes(), len(cs), a.seed), fontsize=10)
        fig.savefig(overview, format="svg", bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        print("[warn] 总览图失败：%s" % e)
    json.dump({"seed": a.seed, "resolution": a.resolution, "n_communities": len(cs),
               "n_nodes": G.number_of_nodes(), "n_edges": G.number_of_edges(),
               "enrich_filter": bool(pathways), "adjp_cut": a.adjp_cut if pathways else None,
               "merged_low_support": merged,
               "orphan_modules": orphans if a.keep_orphan else [],
               "keep_orphan": bool(a.keep_orphan),
               "note": "biological_support=false 的社区无显著富集通路；keep_orphan=0 时合并，=1 时单独作 Unknown module；勿据此声称功能模块",
               "communities": manifest, "overview": "overview.svg"},
              open(os.path.join(a.outdir, "subnetwork_manifest.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[subnetwork] 图 %d 节点 %d 边 → %d 个社区（seed=%d）" %
          (G.number_of_nodes(), G.number_of_edges(), len(cs), a.seed))
    for m in manifest:
        print("  %s: %d 节点 %d 边, top=%s" % (m["file"], m["n_nodes"], m["n_edges"], m["top_degree"][:3]))
    print("[subnetwork] → %s" % a.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
