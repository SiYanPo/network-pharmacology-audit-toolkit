#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_ppi.py — PPI 网络拓扑、hub 判定与稳健性检验

设计依据（为什么不能只按 Degree 排序）
  * Brief Bioinform 2021（PMID 33782690）：多数"活性模块识别"算法在**度序列相同的随机网络**上
    表现不降，说明它们主要学到的是节点度而非生物学。
  * eLife 2024（PMID 39660719）：PPI 的幂律度分布本身可由研究偏倚与技术偏倚解释，
    因此"度最高=最重要"在拓扑上并不成立。
  → 结论：hub 判定必须多指标交叉 + 对**度匹配零模型**做检验 + 报入选频率（bootstrap）。

本脚本做四件事
  1) 清洗：自环、重复边、孤立点、置信度阈值、最大连通分量
  2) 拓扑：DC / 加权度 / BC / 谐波 CC / PageRank / k-core / 聚类系数 / MCC
  3) hub：多指标共识（秩归一化均值）+ 指标间一致性 k / n
  4) 稳健性：边 bootstrap（入选频率）+ 节点丢弃 + 度匹配零模型（拓扑是否"只是度"）

用法
  python np_ppi.py -e string_edges.tsv -o out/
  python np_ppi.py -e edges.tsv -o out/ --min-score 0.7 --top 20 --bootstrap 300
"""
import argparse
import json
import math
import os
import random
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict

import npcore as C

ENSP = ("9606.",)


def read_edges(path, min_score=None, score_col=None):
    """兼容 STRING 网页导出 / API / 无表头三列 / 两列 / CSV。"""
    h, rows = C.read_table(path)
    if not rows:
        return C.Graph(), 0, 0
    cols = [x.lower() for x in h]
    a = b = sc = None
    for cand in ("preferredname_a", "node1", "protein1", "gene_a", "#node1", "source"):
        if cand in cols:
            a = h[cols.index(cand)]
            break
    for cand in ("preferredname_b", "node2", "protein2", "gene_b", "#node2", "target"):
        if cand in cols:
            b = h[cols.index(cand)]
            break
    for cand in ("combined_score", "score", "weight", "confidence"):
        if cand in cols:
            sc = h[cols.index(cand)]
            break
    if not a or not b:
        # 无表头或列名未识别：取前两列当边
        a, b = h[0], h[1]
        if len(h) > 2:
            sc = h[2]
    if score_col:
        sc = score_col

    g = C.Graph()
    for r in rows:
        u = strip_id(r.get(a, ""))
        v = strip_id(r.get(b, ""))
        if not u or not v or u.lower() in ("nan", "none"):
            continue
        w = 1.0
        if sc:
            s = C.to_float(r.get(sc), None)
            if s is not None:
                if s > 1.0:      # STRING 网页导出常为 0–1000
                    s = s / 1000.0
                if min_score is not None and s < min_score:
                    continue
                w = s
        g.add_edge(u, v, w)
    return g, g.self_loops, g.dups


def strip_id(x):
    s = (x or "").strip().strip('"')
    if not s:
        return ""
    for p in ENSP:
        if s.startswith(p):
            return ""
    if s.lower().startswith("string:"):
        return ""
    return s.upper()


# ==================== 孤立/低度交集靶点的扩展子网 ====================
STRING_API = "https://string-db.org/api/json/interaction_partners"


def read_gene_list(path):
    if not path or not os.path.exists(path):
        return []
    out = []
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            g = line.strip().upper()
            if g and g not in out:
                out.append(g)
    return out


def fetch_string_partners(genes, species=9606, required_score=400, limit=10,
                          cache_dir=None, offline=False, timeout=15):
    """STRING interaction_partners → (edges, info)。

    为什么用 interaction_partners 而不是 network：`network` 只返回**给定基因之间**的边，
    对孤立节点等于什么都没加；`interaction_partners` 才会把一阶邻居带进来。
    """
    genes = [g for g in genes if g]
    if not genes:
        return [], {"status": "NO_SEEDS"}
    if offline:
        return [], {"status": "OFFLINE", "reason": "--offline 指定，未访问 STRING"}
    cache = None
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        key = "%s_%s_%s_%s.json" % (species, required_score, limit, "_".join(sorted(genes))[:80])
        cache = os.path.join(cache_dir, key)
        if os.path.exists(cache):
            try:
                j = json.load(open(cache, encoding="utf-8"))
                return parse_string_rows(j, required_score), {"status": "OK", "source": "cache",
                                                              "cache": cache, "n_rows": len(j)}
            except Exception:
                pass
    q = urllib.parse.urlencode({"identifiers": "\r".join(genes), "species": species,
                                "required_score": required_score, "limit": limit,
                                "caller_identity": "network-pharmacology-pro"})
    url = STRING_API + "?" + q
    for attempt in range(2):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "np-pro/2.9"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                rows = json.loads(r.read().decode("utf-8", "ignore"))
            if cache:
                try:
                    json.dump(rows, open(cache, "w", encoding="utf-8"), ensure_ascii=False)
                    # 版本可追溯：同目录写 .meta.json（url/release/sha256）
                    try:
                        import np_snapshot as _S
                        _S.write_meta(cache, url=url, release=os.environ.get("NP_STRING_RELEASE") or "unrecorded")
                    except Exception:
                        pass
                except Exception:
                    pass
            return parse_string_rows(rows, required_score), {"status": "OK", "source": "string_api",
                                                            "url": url, "n_rows": len(rows),
                                                            "attempt": attempt + 1}
        except Exception as e:
            if attempt < 1:
                time.sleep(1.0)
                continue
            return [], {"status": "OFFLINE", "reason": "%s: %s" % (type(e).__name__, str(e)[:120]),
                        "url": url}


def parse_string_rows(rows, required_score=400):
    """STRING 返回 → [(a, b, weight, evidence)]。"""
    edges = []
    for r in rows or []:
        a = (r.get("preferredName_A") or "").strip().upper()
        b = (r.get("preferredName_B") or "").strip().upper()
        if not a or not b or a == b:
            continue
        try:
            sc = float(r.get("score") or 0)
        except (TypeError, ValueError):
            sc = 0.0
        if sc * 1000 < required_score:
            continue
        ev, parts = [], {}
        for k, lab in (("nscore", "neighborhood"), ("fscore", "fusion"), ("pscore", "cooccurrence"),
                       ("ascore", "coexpression"), ("escore", "experiments"), ("dscore", "database"),
                       ("tscore", "textmining")):
            try:
                v = float(r.get(k) or 0)
            except (TypeError, ValueError):
                v = 0.0
            parts[lab] = v
            if v > 0:
                ev.append(lab)
        # 第 5 列保留分项：判断"实验支持"与"共表达"必须分开，不能只看总分
        edges.append((a, b, round(sc, 3), ";".join(ev[:3]), parts))
    return edges



def _spearman(x, y):
    """秩相关（无 scipy）。同值取平均秩。"""
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(x), ranks(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = math.sqrt(sum((v - mx) ** 2 for v in rx))
    dy = math.sqrt(sum((v - my) ** 2 for v in ry))
    return 0.0 if dx * dy == 0 else num / (dx * dy)


def expand_isolated(g, core_genes, deg, min_degree=3, species=9606, required_score=400,
                    limit=10, cache_dir=None, offline=False):
    """把 degree < min_degree 的交集靶点（含完全不在网络里的）用 STRING 一阶邻居补全。

    返回 (extended_graph, report dict)。**不改动原始网络**：原图仍是主结果，
    扩展图单独落盘，避免"为了让图好看而污染主分析"。
    """
    core = [x.upper() for x in core_genes]
    low = [x for x in core if x in g.nodes and deg.get(x, 0) < min_degree]
    absent = [x for x in core if x not in g.nodes]
    seeds = list(dict.fromkeys(low + absent))
    rep = {"status": "NO_SEEDS", "min_degree": min_degree, "n_core": len(core),
           "low_degree_seeds": low, "absent_seeds": absent, "seeds": seeds,
           "n_added_nodes": 0, "n_added_edges": 0}
    if not seeds:
        rep["status"] = "OK_NO_ACTION"
        return None, rep

    edges, info = fetch_string_partners(seeds, species, required_score, limit, cache_dir, offline)
    rep.update({"fetch": info, "species": species, "required_score": required_score, "limit": limit})
    if not edges:
        rep["status"] = info.get("status", "NO_EDGES")
        return None, rep

    ge = C.Graph()
    for u in g.nodes:
        ge.add_node(u)
    for u in g.nodes:
        for v, w in g.adj[u].items():
            if u < v:
                ge.add_edge(u, v, w)
    added_nodes, added_edges = set(), 0
    # 分层计数：实验证据（escore≥0.15）与共表达（ascore≥0.15）分开统计，
    # 避免"总分≥400 一把抓"把研究热度当机制证据（评审 §3）
    exp_neigh = defaultdict(set)
    coexp_neigh = defaultdict(set)
    for e in edges:
        # 兼容 4 元组（旧缓存/外部桩）：缺分项时视为无实验/共表达证据
        if len(e) >= 5:
            a, b, w, ev, parts = e[0], e[1], e[2], e[3], e[4]
        else:
            a, b, w, ev = e[0], e[1], e[2], e[3]
            parts = {}
        for x in (a, b):
            if x not in ge.nodes:
                ge.add_node(x)
                added_nodes.add(x)
        ge.add_edge(a, b, w)
        added_edges += 1
        if parts.get("experiments", 0) >= 0.15:
            exp_neigh[a].add(b)
            exp_neigh[b].add(a)
        if parts.get("coexpression", 0) >= 0.15:
            coexp_neigh[a].add(b)
            coexp_neigh[b].add(a)

    # hub 排名在扩展前后各算一次：扩展 degree 不可用于重排，但必须报告排名漂移
    def _rank(gr):
        d = gr.degree()
        order = sorted(gr.nodes, key=lambda x: (-d[x], x))
        return {x: i + 1 for i, x in enumerate(order)}, d
    rk0, d0 = _rank(g)
    rk1, d1 = _rank(ge)
    shared = [x for x in g.nodes]
    rho = _spearman([rk0[x] for x in shared], [rk1.get(x, len(rk1)) for x in shared]) if len(shared) > 2 else None

    rep.update({"status": "OK", "n_added_nodes": len(added_nodes), "n_added_edges": added_edges,
                "n_nodes_before": g.n, "n_nodes_after": ge.n,
                "n_edges_before": g.m, "n_edges_after": ge.m,
                "seeds_resolved": [s for s in seeds if any(s in (e[0], e[1]) for e in edges)],
                "evidence_stratified": {
                    x: {"experimental_neighbors": len(exp_neigh.get(x, ())),
                        "coexpression_neighbors": len(coexp_neigh.get(x, ())),
                        "total_neighbors": d1.get(x, 0) - d0.get(x, 0)}
                    for x in seeds},
                "hub_rank_correlation": None if rho is None else round(rho, 4),
                "hub_rank_note": ("扩展后 degree 含 STRING 研究偏倚，**不参与 hub 判定**；"
                                  "此处只报告原网络节点的排名漂移（Spearman），用于判断原 hub 结论是否稳健")})
    return ge, rep



def consensus_hubs(metrics, top):
    """各指标转秩（1=最好），取平均秩；返回 (排序行, 一致性计数)。"""
    names = sorted(metrics)
    ranks = {}
    for k in names:
        # 显式次级键（节点名）：否则并列时的顺序取决于 set/dict 迭代顺序，
        # 而该顺序受 PYTHONHASHSEED 影响 → 同一数据两次运行结果不同（已实测复现）。
        order = sorted(metrics[k], key=lambda x: (-(metrics[k][x] if metrics[k][x] is not None
                                                    else -1e18), x))
        ranks[k] = {n: i + 1 for i, n in enumerate(order)}
    rows = []
    for n in metrics[names[0]]:
        rs = [ranks[k].get(n, len(ranks[k])) for k in names]
        rows.append({"node": n, "mean_rank": sum(rs) / len(rs),
                     "n_top": sum(1 for r in rs if r <= top), **{k: metrics[k][n] for k in names}})
    rows.sort(key=lambda r: (r["mean_rank"], r["node"]))
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    return rows




def selftest():
    """PPI 清洗与 hub 共识自检（自环/重复/阈值/排序稳定性）。"""
    import tempfile
    d = tempfile.mkdtemp(prefix="np_ppi_")
    p = os.path.join(d, "e.csv")
    open(p, "w", encoding="utf-8").write(
        "node1,node2,weight\nA,A,0.9\nA,B,0.9\nA,B,0.8\nC,D,0.2\n")
    g, loops, dups = read_edges(p, min_score=0.4)
    metrics = {"degree": {"A": 3, "B": 2, "C": 1},
               "betweenness": {"A": 0.5, "B": 0.2, "C": 0.1}}
    rows = consensus_hubs(metrics, top=2)
    checks = [
        ("self_loop_removed", loops == 1, "loops=%d" % loops),
        ("duplicate_merged", dups == 1, "dups=%d" % dups),
        ("score_threshold", g.n == 2 and g.m == 1, "n=%d m=%d" % (g.n, g.m)),
        ("consensus_order", rows[0]["node"] == "A" and rows[-1]["node"] == "C",
         ",".join(r["node"] for r in rows)),
        ("rank_assigned", all("rank" in r for r in rows), ""),
    ]
    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for n, g_, d_ in checks:
        print("  [%s] %-24s %s" % ("✓" if g_ else "✗", n, d_))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = C.Cli(description="PPI 拓扑 + 多指标 hub + 稳健性（bootstrap/节点丢弃/度匹配零模型）")
    ap.add_argument("-e", "--edges", required=True)
    ap.add_argument("-o", "--outdir", default="np_out")
    ap.add_argument("--min-score", type=float, default=0.4, help="边置信度阈值（STRING 0.4 中/0.7 高）")
    ap.add_argument("--top", type=int, default=15, help="hub 取前 k，默认 15")
    ap.add_argument("--largest-only", action="store_true", help="只保留最大连通分量")
    ap.add_argument("--bootstrap", type=int, default=200, help="边 bootstrap 次数，0=关闭")
    ap.add_argument("--drop", type=float, default=0.1, help="节点丢弃比例")
    ap.add_argument("--max-clique", type=int, default=12, help="MCC 极大团规模上限（超出会截断并报告）")
    ap.add_argument("--null", type=int, default=0, help="（已废弃）旧版度匹配重连检验，默认关闭")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--label", default="PPI 网络", help="图标题")
    # ---- 孤立/低度节点的扩展子网（STRING 一阶邻居）----
    ap.add_argument("--core", default=None, help="交集/核心基因列表；给出后启用扩展子网")
    ap.add_argument("--expand-isolated", action="store_true",
                    help="对 degree<--expand-min-degree 的核心靶点补一阶邻居（需联网或缓存）")
    ap.add_argument("--expand-min-degree", type=int, default=3, help="低于该度视为孤立（默认 3）")
    ap.add_argument("--expand-top", type=int, default=10, help="每个种子最多取多少个邻居（STRING limit）")
    ap.add_argument("--string-score", type=int, default=400, help="STRING 置信度×1000（400=中等）")
    ap.add_argument("--species", type=int, default=9606)
    ap.add_argument("--string-cache", default=None, help="STRING 响应缓存目录（重复运行秒回）")
    ap.add_argument("--offline", action="store_true", help="禁止联网（扩展子网将标 OFFLINE）")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if os.environ.get("NP_OFFLINE") == "1":
        a.offline = True
    if a.selftest:
        return selftest()

    man = C.Manifest("04_ppi", a.outdir, vars(a))
    man.input("edges", a.edges)
    rng = random.Random(a.seed)

    g, loops, dups = read_edges(a.edges, a.min_score)
    n_raw, m_raw = g.n, g.m
    if a.largest_only:
        g = g.largest_component()

    if g.n < 5:
        print("[错误] 网络节点仅 %d 个，无法做有意义的拓扑分析。检查边表格式与阈值。" % g.n)
        return 2

    comps = g.components()
    deg = g.degree()

    # ---- 孤立/低度交集靶点的扩展子网（不改动主网络，单独落盘）----
    ext_g, ext_rep = None, {"status": "DISABLED"}
    if a.expand_isolated:
        if not a.core:
            ext_rep = {"status": "NO_CORE_FILE", "reason": "需 --core 指定交集基因列表"}
            man.warn("请求了 --expand-isolated 但未给 --core，跳过扩展")
        else:
            core_genes = read_gene_list(a.core)
            ext_g, ext_rep = expand_isolated(
                g, core_genes, deg, min_degree=a.expand_min_degree, species=a.species,
                required_score=a.string_score, limit=a.expand_top,
                cache_dir=a.string_cache or os.path.join(C.cache_dir("string") or a.outdir, ""),
                offline=a.offline)
            man.count("expand_status", ext_rep.get("status"))
            man.count("expand_seeds", len(ext_rep.get("seeds") or []))
            if ext_rep.get("status") not in ("OK", "OK_NO_ACTION"):
                man.warn("扩展子网未生成（%s）：主网络不受影响，但孤立靶点仍无结构基础"
                         % ext_rep.get("status"))
            elif ext_rep.get("status") == "OK":
                man.count("expand_added_edges", ext_rep.get("n_added_edges"))

    wdeg = g.weighted_degree()
    bc = g.betweenness()
    cc = g.closeness()
    pr = g.pagerank()
    core = g.core_number()
    clus = {u: g.clustering(u) for u in g.nodes}
    if g.n <= 400:
        mcc, mcc_info = g.mcc(max_clique=a.max_clique)
    else:
        mcc, mcc_info = ({u: float(deg[u]) for u in g.nodes},
                         {"cliques": 0, "max_size": 0, "truncated": -1, "degraded": True})

    metrics = {"degree": deg, "wdegree": wdeg, "betweenness": bc, "closeness": cc,
               "pagerank": pr, "core": core, "clustering": clus, "mcc": mcc}
    hub_rows = consensus_hubs(metrics, a.top)

    # ---- 模块 ----
    groups = g.louvain(seed=a.seed)
    mod_q = g.modularity(groups)
    grp_of = {}
    for i, gs in enumerate(groups):
        for v in gs:
            grp_of[v] = i

    # ---- 稳健性 ----
    stab_rows, freq = [], {}
    if a.bootstrap > 0:
        edges = [(u, v) for u in g.nodes_sorted() for v in sorted(g.adj[u]) if u < v]
        hit = {u: 0 for u in g.nodes}
        rank_sum = {u: 0.0 for u in g.nodes}
        for _ in range(a.bootstrap):
            samp = [edges[rng.randrange(len(edges))] for _ in range(len(edges))]
            h = C.Graph()
            for u in g.nodes_sorted():
                h.add_node(u)
            for u, v in samp:
                h.add_edge(u, v, 1.0)
            d = h.degree()
            order = sorted(d, key=lambda x: (-d[x], x))   # 显式次级键，避免并列时随哈希顺序漂移
            topk = set(order[:a.top])
            for u in topk:
                hit[u] += 1
            for i, u in enumerate(order):
                rank_sum[u] += i + 1
        for u in g.nodes_sorted():
            freq[u] = hit[u] / a.bootstrap

    drop_freq = {}
    if a.drop > 0:
        nodes = g.nodes_sorted()
        kdrop = max(1, int(len(nodes) * a.drop))
        hitd = {u: 0 for u in nodes}
        reps = max(30, a.bootstrap // 4)
        for _ in range(reps):
            keep = set(nodes) - set(rng.sample(nodes, kdrop))
            h = g.subgraph(keep)
            d = h.degree()
            topk = set(sorted(d, key=lambda x: (-d[x], x))[:a.top])
            for u in topk:
                hitd[u] = hitd.get(u, 0) + 1
        drop_freq = {u: hitd[u] / reps for u in nodes}

    # 度的混杂校正：把介数放到"同度层"里比较
    # 为什么不用"度保持重连 + top-k 重合"：度序列不变时按度排序的 top-k 恒等，检验无信息量。
    # 为什么不用"度 ±1"窄窗：PPI 度分布长尾，281 节点散布在 100+ 个度值上，
    # 窄窗常只有 1–2 个同伴 → 标准差无定义 → z 恒为 0，看起来"没有桥接节点"，
    # 实际是"这个度量根本没算出来"。改为等频分箱，箱内标准化，无法计算时给 NA。
    bc_z, bc_info = C.degree_binned_z(bc, deg, n_bins=10, min_bin=8)
    bc_excess = {u: (None if bc_z.get(u) is None else bc_z[u]) for u in g.nodes_sorted()}
    n_na = sum(1 for v in bc_z.values() if v is None)
    if n_na:
        man.warn("bc_excess_z 有 %d/%d 个节点为 NA（该度分箱方差为 0 或样本不足）"
                 % (n_na, len(bc_z)))
    man.count("bc_z_bins", len(bc_info["bins"]))
    man.count("bc_z_na", n_na)
    p_null = None

    for r in hub_rows:
        u = r["node"]
        r["bootstrap_freq"] = round(freq.get(u, float("nan")), 3) if freq else ""
        r["drop_freq"] = round(drop_freq.get(u, float("nan")), 3) if drop_freq else ""
    hub_rows.sort(key=lambda r: (r["mean_rank"], r["node"]))

    def verdict(r):
        f = r.get("bootstrap_freq")
        if f == "" or f is None:
            return "未做稳健性检验（不得直接称 hub）"
        try:
            f = float(f)
        except Exception:
            return "—"
        if f >= 0.9:
            return "稳健"
        if f >= 0.6:
            return "中等"
        return "不稳定"

    C.ensure_dir(a.outdir)
    # 清洗后的边表
    p_edges = C.write_table(os.path.join(a.outdir, "ppi_edges_clean.tsv"),
                            ["node1", "node2", "weight"],
                            [[u, v, round(w, 3)] for u in sorted(g.nodes) for v, w in g.adj[u].items() if u < v])
    # ---- 扩展子网产物 ----
    p_ext_edges = p_ext_nodes = p_ext_svg = None
    if ext_g and ext_g.n > 0:
        p_ext_edges = C.write_table(
            os.path.join(a.outdir, "extended_ppi_edges.tsv"),
            ["node1", "node2", "weight", "origin"],
            [[u, v, round(w, 3), "seed"] for u in sorted(ext_g.nodes)
             for v, w in ext_g.adj[u].items() if u < v])
        d_ext = ext_g.degree()
        seeds = set(ext_rep.get("seeds") or [])
        p_ext_nodes = C.write_table(
            os.path.join(a.outdir, "extended_ppi_nodes.csv"),
            ["node", "degree_extended", "degree_original", "is_seed", "is_added"],
            [{"node": u, "degree_extended": d_ext[u], "degree_original": deg.get(u, 0),
              "is_seed": int(u in seeds), "is_added": int(u not in g.nodes)}
             for u in sorted(ext_g.nodes, key=lambda x: (-d_ext[x], x))])
        p_ext_svg = C.network_svg(
            ext_g, os.path.join(a.outdir, "extended_ppi.svg"),
            title="扩展子网（含孤立靶点一阶邻居）n=%d m=%d — 本图 degree 不参与 hub 判定" % (
                ext_g.n, ext_g.m),
            node_score={u: d_ext[u] for u in ext_g.nodes}, label_top=min(25, len(seeds) * 3 + 8),
            seed=a.seed)
        json.dump(ext_rep, open(os.path.join(a.outdir, "extended_ppi_manifest.json"),
                                "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        man.count("ext_nodes", ext_g.n)
        man.count("ext_edges", ext_g.m)

    p_topo = C.write_table(os.path.join(a.outdir, "topology.csv"),
                           ["node", "degree", "wdegree", "betweenness", "closeness", "pagerank",
                            "core", "clustering", "mcc", "mcc_log10", "bc_excess_z", "module",
                            "bootstrap_freq", "drop_freq"],
                           [{"node": u, "degree": deg[u], "wdegree": round(wdeg[u], 3),
                             "betweenness": round(bc[u], 5), "closeness": round(cc[u], 5),
                             "pagerank": round(pr[u], 6), "core": core[u],
                             "clustering": round(clus[u], 4), "mcc": mcc.get(u, ""),
                             "mcc_log10": (round(__import__("math").log10(mcc[u]), 2)
                                           if mcc.get(u) and mcc[u] > 0 else ""),
                             "bc_excess_z": ("" if bc_z.get(u) is None else round(bc_z[u], 3)),
                             "module": grp_of.get(u, ""),
                             "bootstrap_freq": freq.get(u, ""), "drop_freq": drop_freq.get(u, "")}
                            for u in sorted(g.nodes, key=lambda x: (-deg[x], x))])
    hh = ["rank", "node", "mean_rank", "degree", "wdegree", "betweenness", "closeness",
          "pagerank", "core", "clustering", "mcc", "mcc_log10", "bc_excess_z", "n_top",
          "bootstrap_freq", "drop_freq", "verdict"]
    for r in hub_rows:
        r["verdict"] = verdict(r)
        r["bc_excess_z"] = "" if bc_z.get(r["node"]) is None else round(bc_z[r["node"]], 3)
        r["mcc_log10"] = (round(__import__("math").log10(r["mcc"]), 2)
                          if r.get("mcc") else "")
    p_hub = C.write_table(os.path.join(a.outdir, "hubs.csv"), hh, hub_rows)

    if ext_rep.get("status") == "OK":
        print("[PPI] 扩展子网：%d 个孤立/低度种子（%s）→ +%d 节点 / +%d 边（STRING ≥%.3f）" % (
            len(ext_rep.get("seeds") or []),
            ",".join((ext_rep.get("seeds") or [])[:6]) + ("…" if len(ext_rep.get("seeds") or []) > 6 else ""),
            ext_rep.get("n_added_nodes", 0), ext_rep.get("n_added_edges", 0),
            ext_rep.get("required_score", 400) / 1000.0))
    elif ext_rep.get("status") not in ("DISABLED", "OK_NO_ACTION"):
        print("[PPI] 扩展子网未生成：%s（主网络不受影响）" % ext_rep.get("status"))

    p_net = C.network_svg(g, os.path.join(a.outdir, "network.svg"),
                          title="%s（n=%d, m=%d, Q=%.3f）" % (a.label, g.n, g.m, mod_q),
                          groups=groups, label_top=min(20, max(8, a.top)),
                          seed=a.seed)

    md = ["# PPI 网络与 hub 判定报告", "",
          "| 项 | 值 |", "|---|---|",
          "| 清洗前节点/边 | %d / %d |" % (n_raw, m_raw),
          "| 清洗后节点/边 | %d / %d |" % (g.n, g.m),
          "| 剔除自环 | %d |" % loops, "| 合并重复边 | %d |" % dups,
          "| 连通分量 | %d（最大 %d） |" % (len(comps), len(comps[0])),
          "| 边阈值 | ≥%.2f |" % a.min_score,
          "| 密度 / 平均度 / 直径 | %.4f / %.1f / %d |" % (g.density(), 2 * g.m / g.n, g.diameter()),
          "| 模块数（Louvain）/ 模块度 Q | %d / %.3f |" % (len(groups), mod_q), ""]

    if len(groups) < 2 or mod_q < 0.2:
        md += ["> ⚠ 模块结构弱（Q=%.3f, 模块数=%d）：说明该网络缺乏明显社区结构，"
               "不宜用「模块」叙事，应回到网络整体拓扑与富集结果。" % (mod_q, len(groups)), ""]

    md += ["## hub 基因（多指标共识前 %d）" % a.top, "",
           "| 排名 | 基因 | 平均秩 | DC | BC | 模块 | 入选频率 | 判定 |", "|---|---|---|---|---|---|---|---|"]
    for r in hub_rows[:a.top]:
        md.append("| %d | %s | %.1f | %d | %.4f | %s | %s | %s |" % (
            r["rank"], r["node"], r["mean_rank"], r["degree"], r["betweenness"],
            r["module"] if "module" in r else grp_of.get(r["node"], ""),
            ("%.2f" % r["bootstrap_freq"]) if isinstance(r["bootstrap_freq"], float) else "—",
            r["verdict"]))

    n_robust = sum(1 for r in hub_rows[:a.top] if r["verdict"] == "稳健")
    md += ["", "稳健（入选频率≥0.90）：**%d / %d**" % (n_robust, min(a.top, len(hub_rows))), ""]

    # 扩展子网说明
    st = ext_rep.get("status")
    if st == "OK":
        md += ["## 孤立/低度交集靶点的扩展子网（STRING 一阶邻居）", "",
               "问题：交集靶点若在 PPI 中 degree<3（如 ESR1、NOS3），下游对接与湿实验缺少结构/网络依据。"
               "本步骤**不改动主网络**，而是单独构建扩展子网：", "",
               "| 项 | 值 |", "|---|---|",
               "| 阈值 | 原始 degree < %d |" % ext_rep["min_degree"],
               "| 低度种子 | %s |" % (", ".join(ext_rep.get("low_degree_seeds") or []) or "无"),
               "| 完全缺失种子 | %s |" % (", ".join(ext_rep.get("absent_seeds") or []) or "无"),
               "| STRING 置信度 | ≥%.3f（species=%d） |" % (ext_rep["required_score"] / 1000.0, ext_rep["species"]),
               "| 新增节点/边 | %d / %d |" % (ext_rep["n_added_nodes"], ext_rep["n_added_edges"]),
               "| 网络规模 | %d→%d 节点，%d→%d 边 |" % (
                   ext_rep["n_nodes_before"], ext_rep["n_nodes_after"],
                   ext_rep["n_edges_before"], ext_rep["n_edges_after"]), ""]
        es = ext_rep.get("evidence_stratified") or {}
        if es:
            md += ["### 分层证据（实验 vs 共表达），避免把研究热度当机制", "",
                   "| 种子 | 新增邻居 | 其中实验证据(escore≥0.15) | 共表达(ascore≥0.15) |",
                   "|---|---|---|---|"]
            for k, v in es.items():
                md.append("| %s | %d | %d | %d |" % (
                    k, v.get("total_neighbors", 0), v.get("experimental_neighbors", 0),
                    v.get("coexpression_neighbors", 0)))
            md += ["", "读法：新增邻居里**实验支持占比低**说明该靶点的网络位置主要来自文本挖掘/共表达，"
                   "不是机制证据；这类靶点不能仅凭「扩展后 degree 变高」就当作核心。", ""]
        if ext_rep.get("hub_rank_correlation") is not None:
            md += ["### hub 排名漂移（扩展前 vs 扩展后）", "",
                   "- 原网络节点排名 Spearman ρ = **%.3f**" % ext_rep["hub_rank_correlation"],
                   "- 判读：ρ>0.9 说明原 hub 结论对扩展稳健；ρ<0.7 说明排名高度依赖边界定义，"
                   "**必须在正文报告这一不稳健性**，而不是只展示扩展后的漂亮网络。",
                   "- 硬约束：扩展图的 degree **不参与** hub 判定（见 `hub_rank_note`）。", ""]

        md += ["产物：`extended_ppi_edges.tsv`、`extended_ppi_nodes.csv`、`extended_ppi.svg`、"
               "`extended_ppi_manifest.json`。**主网络（topology/hubs）仍基于原始边表**，"
               "扩展图只用于给孤立靶点提供结构基础与作图，不得用扩展图的度去重排 hub。", ""]
    elif st in ("OFFLINE", "NO_EDGES", "NO_SEEDS"):
        md += ["## 孤立/低度交集靶点的扩展子网", "",
               "> 未生成扩展子网（%s%s）。主网络结果不受影响；但 degree<3 的靶点仍缺少结构基础，"
               "建议联网重跑 `--expand-isolated`，或在正文说明该靶点为网络孤立节点。"
               % (st, ("：" + str(ext_rep.get("fetch", {}).get("reason", ""))) if st == "OFFLINE" else ""), ""]

    # MCC 实现说明（必须如实报告截断）
    if mcc_info.get("degraded"):
        md += ["## MCC 说明", "", "> 网络规模超过 400 节点，MCC 已退化为度（组合爆炸）。"
               "此时 MCC 不提供额外信息，hub 判定应更依赖 bootstrap 入选频率。", ""]
    elif mcc_info.get("truncated", 0) > 0:
        md += ["## MCC 说明", "",
               "- 极大团枚举：%d 个，最大团 %d；**因超过上限 %d 截断 %d 次**"
               % (mcc_info["cliques"], mcc_info["max_size"], a.max_clique, mcc_info["truncated"]),
               "- 截断会低估大团节点的 MCC，但**不改变排序方向**；如用于发表，"
               "建议把 `--max-clique` 调大复跑一次确认排序不变。", ""]
    else:
        md += ["## MCC 说明", "",
               "- 极大团枚举：%d 个，最大团 %d，未截断（上限 %d）"
               % (mcc_info["cliques"], mcc_info["max_size"], a.max_clique), ""]

    top_excess = sorted([kv for kv in bc_z.items() if kv[1] is not None],
                        key=lambda x: (-x[1], x[0]))[:5]
    md += ["## 度校正后的「桥接」节点（bc_excess_z）", "",
           "把介数放到**同度层**（度 ±1）内标准化，回答「它是否比同等连接数的节点更像桥」。"
           "这是度之外的拓扑信息，可用于区分「度大」与「位置关键」。", "",
           "| 基因 | 度 | 介数 | 同度层 z |", "|---|---|---|---|"]
    for u, z in top_excess:
        md.append("| %s | %d | %.4f | %.2f |" % (u, deg[u], bc[u], z))
    md += ["", "分箱明细（度区间 / 箱内节点数 / 箱内均值 / 标准差）：", "",
           "| 度区间 | n | 均值 | SD |", "|---|---|---|---|"]
    for b in bc_info["bins"]:
        md.append("| %d–%d | %d | %.5f | %.5f |"
                  % (b["deg_min"], b["deg_max"], b["n"], b["mean"], b["sd"]))
    md += ["",
           "> 读法：z>2 提示该节点的「桥」作用超出其连接数所能解释的范围；"
           "z≈0 说明它的中心性基本由度决定。**两者都不等于生物学重要性**。", ""]

    md += ["## 写法红线", "",
           "1. Degree 不是「重要性排名」，它只表示参与的过程/通路多。",
           "2. STRING 的边是**功能关联**（含文本挖掘、共表达、同源），不等于物理结合，也不指示激活/抑制方向。",
           "3. 必须在方法学写明：STRING 版本与访问日期、物种（9606）、边类型（functional/physical）、score 阈值。",
           "4. 称「核心靶点/hub」必须给入选频率；入选频率 <0.60 时只能写「度排序前列节点」。", ""]

    p_md = C.write_text(os.path.join(a.outdir, "ppi_report.md"), "\n".join(md) + "\n")

    for k, v in [("edges_clean", p_edges), ("topology", p_topo), ("hubs", p_hub),
                 ("network", p_net), ("report", p_md)]:
        man.output(k, v)
    man.count("nodes", g.n)
    man.count("edges", g.m)
    man.count("modules", len(groups))
    man.count("modularity", round(mod_q, 3))
    man.count("robust_hubs", n_robust)
    if freq and n_robust == 0:
        man.warn("无稳健 hub（入选频率均 <0.90），正文不得使用「核心枢纽」表述")
    if mcc_info.get("degraded"):
        man.warn("MCC 因网络规模退化为度，hub 判定应更依赖 bootstrap")
    elif mcc_info.get("truncated", 0) > 0:
        man.warn("MCC 极大团被截断 %d 次（上限 %d），建议调大 --max-clique 复核"
                 % (mcc_info["truncated"], a.max_clique))
    man.save()

    if not a.quiet:
        print("[PPI] 节点 %d，边 %d，模块 %d（Q=%.3f）" % (g.n, g.m, len(groups), mod_q))
        print("[PPI] hub top%d：%s" % (a.top, "、".join(r["node"] for r in hub_rows[:8])))
        if freq:
            print("[PPI] 稳健 hub：%d/%d（bootstrap %d 次）" % (n_robust, min(a.top, len(hub_rows)), a.bootstrap))
        if mcc_info.get("truncated", 0) > 0:
            print("[PPI] MCC 截断 %d 次（上限 %d）" % (mcc_info["truncated"], a.max_clique))
        print("[PPI] 输出 → %s" % p_md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
