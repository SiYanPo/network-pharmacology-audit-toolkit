#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
npcore.py — 网络药理学研究平台核心库（零第三方依赖）

设计原则
1. 只依赖 Python 3.8+ 标准库；numpy 等仅作可选加速，缺失自动降级。
2. 所有统计量给出可复算的定义（不藏参数），随机过程必须可设种子。
3. 图表直接输出 SVG（不依赖 matplotlib / 中文字体），可在任何浏览器与 Word 中打开。
4. 每个阶段写 manifest，记录输入哈希、参数、随机种子、行数，保证可复现。

模块
    io       表格读写（编码/分隔符/列名别名自适应）
    num      精确超几何、BH-FDR、经验 P、Z 分数、log-sum-exp
    graph    邻接表、中心性、k-core、社区发现、布局
    svg      SVG 绘图原语（柱状、气泡、网络、韦恩、折线）
    prov     溯源清单
    cli      统一参数与中文错误提示
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io as _io
import json
import math
import os
import random
import re
import shutil
import sys
import time
from collections import Counter, defaultdict, deque

VERSION = "1.0.0"

# ----------------------------------------------------------------------------
# io
# ----------------------------------------------------------------------------

_ENCODINGS = ("utf-8-sig", "utf-8", "gb18030", "gbk", "latin-1")
_DELIMS = ("\t", ",", ";", "|")


def _read_text(path):
    if str(path).lower().endswith(".gz"):
        import gzip
        return gzip.open(path, "rt", encoding="utf-8", errors="replace").read(), "utf-8"
    raw = open(path, "rb").read()
    for enc in _ENCODINGS:
        try:
            return raw.decode(enc), enc
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", "replace"), "utf-8?"


def sniff_delim(text):
    """按表头行统计候选分隔符出现次数，取最稳定的一个。"""
    head = []
    for line in text.splitlines():
        if line.strip():
            head.append(line)
        if len(head) >= 5:
            break
    if not head:
        return "\t"
    best, best_score = "\t", -1.0
    for d in _DELIMS:
        counts = [ln.count(d) for ln in head]
        if max(counts) == 0:
            continue
        # 出现次数多且各行一致 → 更可能是真分隔符
        mode = Counter(counts).most_common(1)[0][0]
        consistency = counts.count(mode) / len(counts)
        score = mode * consistency
        if score > best_score:
            best, best_score = d, score
    return best


def read_table(path, has_header=True):
    """读表格 → (headers, rows[list[dict]])。支持 csv/tsv/txt/xlsx。"""
    path = str(path)
    if path.lower().endswith((".xlsx", ".xlsm")):
        return _read_xlsx(path, has_header)
    text, enc = _read_text(path)
    if not text.strip():
        return [], []
    delim = sniff_delim(text)
    rdr = csv.reader(_io.StringIO(text), delimiter=delim)
    rows_raw = [r for r in rdr if any(c.strip() for c in r)]
    if not rows_raw:
        return [], []
    if not has_header:
        width = max(len(r) for r in rows_raw)
        headers = ["c%d" % (i + 1) for i in range(width)]
        body = rows_raw
    else:
        headers = [h.strip().lstrip("\ufeff") for h in rows_raw[0]]
        body = rows_raw[1:]
    headers = _dedupe_headers(headers)
    rows = []
    for r in body:
        if len(r) < len(headers):
            r = list(r) + [""] * (len(headers) - len(r))
        rows.append({headers[i]: (r[i].strip() if i < len(r) else "") for i in range(len(headers))})
    return headers, rows


def _dedupe_headers(headers):
    seen, out = Counter(), []
    for h in headers:
        h = h or "col"
        if seen[h]:
            seen[h] += 1
            out.append("%s_%d" % (h, seen[h]))
        else:
            seen[h] += 1
            out.append(h)
    return out


def _read_xlsx(path, has_header=True):
    try:
        import openpyxl
    except ImportError:
        raise SystemExit("读取 xlsx 需要 openpyxl：python -m pip install openpyxl（或先另存为 csv）")
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    ws = wb[wb.sheetnames[0]]
    grid = [["" if c is None else str(c).strip() for c in row] for row in ws.iter_rows(values_only=True)]
    grid = [r for r in grid if any(r)]
    if not grid:
        return [], []
    if has_header:
        headers, body = _dedupe_headers(grid[0]), grid[1:]
    else:
        headers, body = ["c%d" % (i + 1) for i in range(len(grid[0]))], grid
    rows = []
    for r in body:
        r = list(r) + [""] * (len(headers) - len(r))
        rows.append({headers[i]: r[i] for i in range(len(headers))})
    return headers, rows


def write_table(path, headers, rows, encoding="utf-8-sig"):
    """rows 可为 list[dict] 或 list[list]。"""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding=encoding, newline="") as f:
        w = csv.writer(f)
        w.writerow(headers)
        for r in rows:
            if isinstance(r, dict):
                w.writerow([_fmt_cell(r.get(h, "")) for h in headers])
            else:
                w.writerow([_fmt_cell(x) for x in r])
    return path


def _fmt_cell(x):
    if x is None:
        return ""
    if isinstance(x, float):
        if x != x:            # NaN
            return ""
        if x in (float("inf"), float("-inf")):   # 无穷：CSV 里如实写，别崩
            return "inf" if x > 0 else "-inf"
        if x == int(x) and abs(x) < 1e15:
            return str(int(x))
        return ("%.6g" % x)
    return str(x)


def write_text(path, text, encoding="utf-8"):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding=encoding, newline="\n") as f:
        f.write(text)
    return path


# 列名别名（中英文混排场景）
ALIASES = {
    "gene": ["gene", "symbol", "gene_symbol", "genesymbol", "hgnc", "target", "靶点", "基因",
             "基因名", "基因符号", "target_gene", "hugo_symbol", "genename", "gene_name",
             "preferredname_a", "preferredname_b", "node1", "node2", "protein"],
    "compound": ["compound", "molecule", "chemical", "ingredient", "成分", "化合物", "分子",
                 "molecule_name", "name", "compound_name", "mol_name", "chem"],
    "mol_id": ["mol_id", "molid", "mol", "cid", "pubchem_cid", "id", "编号"],
    "herb": ["herb", "drug", "tcm", "药材", "中药", "药味", "herb_name", "formula"],
    "ob": ["ob", "ob(%)", "ob%", "oral_bioavailability", "口服生物利用度", "ob值"],
    "dl": ["dl", "dl(%)", "drug_likeness", "类药性", "dl值"],
    "hl": ["hl", "half_life", "半衰期", "hl(h)"],
    "content": ["content", "含量", "content_pct", "含量(%)", "百分含量"],
    "score": ["score", "combined_score", "confidence", "关联度", "score_gda", "relevance_score",
              "evidence", "weight", "prob", "probability", "得分"],
    "disease": ["disease", "疾病", "disease_name", "phenotype", "表型"],
    "source": ["source", "来源", "database", "db", "数据来源"],
    "pvalue": ["pvalue", "p_value", "p", "p-value", "pval", "adj_p", "padj", "qvalue", "fdr"],
    "term": ["term", "pathway", "description", "name", "通路", "条目", "term_name", "go", "kegg"],
    "target_name": ["target_name", "protein_name", "靶点名称", "target protein", "protein"],
    "mw": ["mw", "molecular_weight", "分子量", "mol_weight"],
    "smiles": ["smiles", "canonical_smiles", "isomeric_smiles"],
    "herb_count": ["n_herb", "herb_count", "herbs", "药材数"],
}


def resolve_col(headers, key, required=False, extra=None):
    """在表头里找 key 对应的列名；返回 None 表示没找到。"""
    # 注意：key 自身也必须作为候选。否则 resolve_col(h, "organism") 在 ALIASES 里
    # 没有该键时会返回 None —— 即使表头就叫 organism（这是一个静默失效的坑）。
    cands = [key] + list(ALIASES.get(key, [])) + list(extra or [])
    low = {h.lower().strip(): h for h in headers}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    for c in cands:
        for h in headers:
            if c.lower() in h.lower():
                return h
    if required:
        raise SystemExit("未找到【%s】列。现有列：%s" % (key, "、".join(headers[:12])))
    return None


def to_float(x, default=None):
    if x is None:
        return default
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip().replace("%", "").replace(",", "")
    if not s or s.lower() in ("na", "n/a", "nan", "none", "-", "—", "null"):
        return default
    try:
        return float(s)
    except ValueError:
        m = re.search(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", s)
        return float(m.group(0)) if m else default


# ----------------------------------------------------------------------------
# num
# ----------------------------------------------------------------------------

def _logsumexp(vals):
    vals = [v for v in vals if v != -math.inf]
    if not vals:
        return -math.inf
    m = max(vals)
    return m + math.log(sum(math.exp(v - m) for v in vals))


def hypergeom_sf(k, N, K, n):
    """P(X >= k)，X ~ Hypergeom(N, K, n)。log 空间求和，数值稳定。"""
    if N <= 0 or K <= 0 or n <= 0:
        return 1.0
    lo = max(0, n - (N - K))
    hi = min(K, n)
    k = max(k, lo)
    if k > hi:
        return 0.0
    lgn = math.lgamma(N + 1)
    terms = []
    for i in range(k, hi + 1):
        lp = (math.lgamma(K + 1) - math.lgamma(i + 1) - math.lgamma(K - i + 1)
              + math.lgamma(N - K + 1) - math.lgamma(n - i + 1) - math.lgamma(N - K - n + i + 1)
              - lgn + math.lgamma(n + 1) + math.lgamma(N - n + 1))
        terms.append(lp)
    val = math.exp(_logsumexp(terms))
    return min(1.0, max(val, 1e-300))


def bh_fdr(pvals):
    """Benjamini-Hochberg，返回与输入同序的 q 值。"""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: (pvals[i], i))
    q = [1.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        val = pvals[i] * n / rank
        prev = min(prev, val)
        q[i] = min(1.0, prev)
    return q


def empirical_p(obs, null, tail="greater"):
    """经验 P：(1 + #{null >= obs}) / (1 + B)，永不为 0。"""
    b = len(null)
    if b == 0:
        return 1.0, 0.0
    if tail == "greater":
        c = sum(1 for v in null if v >= obs)
    else:
        c = sum(1 for v in null if v <= obs)
    return (1.0 + c) / (1.0 + b), c / b


def zscore(obs, null):
    b = len(null)
    if b < 2:
        return float("nan")
    mu = sum(null) / b
    sd = math.sqrt(sum((v - mu) ** 2 for v in null) / (b - 1))
    if sd == 0:
        return float("inf") if obs > mu else (0.0 if obs == mu else float("-inf"))
    return (obs - mu) / sd


def jaccard(a, b):
    a, b = set(a), set(b)
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def ks_distance(a, b):
    """两样本 Kolmogorov–Smirnov 统计量（经验分布函数的最大差）。

    注意：必须在**并集的所有取值点**上比较 F_a(x) 与 F_b(x)（含等号），
    否则相同样本集会得到非零距离（这是本函数早期版本的 bug）。
    """
    a = sorted(float(x) for x in a)
    b = sorted(float(x) for x in b)
    if not a or not b:
        return float("nan")
    i = j = 0
    d = 0.0
    for v in sorted(set(a) | set(b)):
        while i < len(a) and a[i] <= v:
            i += 1
        while j < len(b) and b[j] <= v:
            j += 1
        d = max(d, abs(i / len(a) - j / len(b)))
    return d


def degree_matched_sample(obs_degrees, pool, degree, rng):
    """按**观察集的度多重集**逐一匹配抽样（而非按总体分箱比例）。

    为什么不能用"按总体比例分箱"：那复现的是总体度分布，而观察集（交集基因）偏向高度节点，
    会得到系统性偏稀疏的零模型 → 反保守（anti-conservative）。

    策略：
      1) 对观察集每个度 d，从 pool 中随机取一个度恰为 d 的节点（不放回）；
      2) 该度已抽空时，取"最近可用度"的节点，并计入失配；
    返回 (sample, quality)。quality 含平均绝对度偏差、KS 距离、精确匹配比例。
    """
    from collections import defaultdict
    buckets = defaultdict(list)
    for n in pool:
        buckets[degree.get(n, 0)].append(n)
    for d in buckets:
        rng.shuffle(buckets[d])

    used = set()
    sample = []
    diffs = []
    exact = 0
    for d in obs_degrees:
        pick = None
        if buckets.get(d):
            while buckets[d]:
                c = buckets[d].pop()
                if c not in used:
                    pick, exact = c, exact + 1
                    break
        if pick is None:
            # 最近可用度
            best = None
            for delta in range(1, 200):
                for dd in (d - delta, d + delta):
                    if buckets.get(dd):
                        while buckets[dd]:
                            c = buckets[dd].pop()
                            if c not in used:
                                best = c
                                break
                        if best:
                            break
                if best:
                    break
            if best is None:
                continue
            pick = best
        used.add(pick)
        sample.append(pick)
        diffs.append(abs(degree.get(pick, 0) - d))

    k = len(obs_degrees)
    quality = {
        "n": len(sample),
        "mean_abs_ddeg": (sum(diffs) / len(diffs)) if diffs else float("nan"),
        "exact_frac": (exact / k) if k else 0.0,
        "ks": ks_distance(obs_degrees, [degree.get(x, 0) for x in sample]),
    }
    return sample, quality


def fmt_p(p):
    if p is None or (isinstance(p, float) and p != p):
        return "NA"
    if p == 0:
        return "<1e-300"
    if p < 1e-4:
        return "%.2e" % p
    return "%.4f" % p




# ----------------------------------------------------------------------------
# graph
# ----------------------------------------------------------------------------

class Graph:
    """无向简单图（自环剔除、重复边合并、权重取最大）。"""

    def __init__(self):
        self.adj = defaultdict(dict)   # u -> {v: w}
        self.nodes = set()
        self.self_loops = 0
        self.dups = 0

    def add_node(self, n):
        self.nodes.add(n)
        self.adj[n]

    def add_edge(self, u, v, w=1.0):
        if not u or not v:
            return
        if u == v:
            self.self_loops += 1
            return
        self.add_node(u)
        self.add_node(v)
        if v in self.adj[u]:
            self.dups += 1
        self.adj[u][v] = max(self.adj[u].get(v, 0.0), float(w))
        self.adj[v][u] = max(self.adj[v].get(u, 0.0), float(w))

    @property
    def n(self):
        return len(self.nodes)

    @property
    def m(self):
        return sum(len(v) for v in self.adj.values()) // 2

    def degree(self):
        return {u: len(v) for u, v in self.adj.items()}

    def weighted_degree(self):
        return {u: sum(v.values()) for u, v in self.adj.items()}


    def nodes_sorted(self):
        """稳定节点顺序：集合迭代顺序受 PYTHONHASHSEED 影响，凡涉及随机过程必须用这个。"""
        return sorted(self.nodes)

    def components(self):
        seen, comps = set(), []
        for s in self.nodes_sorted():
            if s in seen:
                continue
            q, comp = deque([s]), []
            seen.add(s)
            while q:
                u = q.popleft()
                comp.append(u)
                for v in self.adj[u]:
                    if v not in seen:
                        seen.add(v)
                        q.append(v)
            comps.append(comp)
        comps.sort(key=len, reverse=True)
        return comps

    def largest_component(self):
        comps = self.components()
        if not comps:
            return Graph()
        keep = set(comps[0])
        g = Graph()
        for u in keep:
            g.add_node(u)
        for u in keep:
            for v, w in self.adj[u].items():
                if v in keep and u < v:
                    g.add_edge(u, v, w)
        return g

    def subgraph(self, keep):
        keep = set(keep)
        g = Graph()
        for u in keep:
            g.add_node(u)
        for u in keep:
            for v, w in self.adj[u].items():
                if v in keep and u < v:
                    g.add_edge(u, v, w)
        return g

    # ---- 中心性 ----
    def betweenness(self):
        """Brandes 算法（无权）。O(V*E)，<2000 节点可用。"""
        bc = {v: 0.0 for v in self.nodes}
        for s in self.nodes_sorted():
            S, P = [], {v: [] for v in self.nodes}
            sigma = dict.fromkeys(self.nodes, 0.0)
            sigma[s] = 1.0
            d = dict.fromkeys(self.nodes, -1)
            d[s] = 0
            q = deque([s])
            while q:
                v = q.popleft()
                S.append(v)
                for w in self.adj[v]:
                    if d[w] < 0:
                        d[w] = d[v] + 1
                        q.append(w)
                    if d[w] == d[v] + 1:
                        sigma[w] += sigma[v]
                        P[w].append(v)
            delta = dict.fromkeys(self.nodes, 0.0)
            while S:
                w = S.pop()
                for v in P[w]:
                    if sigma[w] > 0:
                        delta[v] += sigma[v] / sigma[w] * (1.0 + delta[w])
                if w != s:
                    bc[w] += delta[w]
        if len(self.nodes) > 2:
            scale = 1.0 / ((len(self.nodes) - 1) * (len(self.nodes) - 2))
            bc = {k: v * scale for k, v in bc.items()}
        return bc

    def closeness(self):
        """调和中心性（对不连通图稳健）：sum(1/d)。"""
        out = {}
        for s in self.nodes_sorted():
            d = {s: 0}
            q = deque([s])
            while q:
                u = q.popleft()
                for v in self.adj[u]:
                    if v not in d:
                        d[v] = d[u] + 1
                        q.append(v)
            ssum = sum(1.0 / dd for dd in d.values() if dd > 0)
            out[s] = ssum / max(1, (len(self.nodes) - 1))
        return out

    def pagerank(self, alpha=0.85, iters=100, tol=1e-9):
        n = self.n
        if n == 0:
            return {}
        order = self.nodes_sorted()
        pr = {v: 1.0 / n for v in order}
        deg = self.degree()
        for _ in range(iters):
            new = {v: (1 - alpha) / n for v in order}
            dangling = sum(pr[v] for v in order if deg[v] == 0)
            for v in order:
                if deg[v] == 0:
                    continue
                share = alpha * pr[v] / deg[v]
                for w in self.adj[v]:
                    new[w] += share
            if dangling:
                for v in order:
                    new[v] += alpha * dangling / n
            err = sum(abs(new[v] - pr[v]) for v in order)
            pr = new
            if err < tol:
                break
        return pr

    def core_number(self):
        deg = {u: len(v) for u, v in self.adj.items()}
        core = {}
        nodes = sorted(deg, key=lambda x: deg[x])
        cur = dict(deg)
        order = []
        remaining = set(nodes)
        k = 0
        while remaining:
            changed = True
            while changed:
                changed = False
                for u in sorted(remaining):
                    if cur[u] <= k:
                        core[u] = k
                        remaining.discard(u)
                        order.append(u)
                        for v in self.adj[u]:
                            if v in remaining:
                                cur[v] -= 1
                        changed = True
            if remaining:
                k = min(cur[v] for v in remaining)
        return core

    def clustering(self, u):
        nb = list(self.adj[u])
        k = len(nb)
        if k < 2:
            return 0.0
        links = 0
        for i in range(k):
            for j in range(i + 1, k):
                if nb[j] in self.adj[nb[i]]:
                    links += 1
        return 2.0 * links / (k * (k - 1))

    def diameter(self, cap=10000):
        best = 0
        count = 0
        for s in self.nodes_sorted():
            d = {s: 0}
            q = deque([s])
            while q:
                u = q.popleft()
                for v in self.adj[u]:
                    if v not in d:
                        d[v] = d[u] + 1
                        q.append(v)
            if d:
                best = max(best, max(d.values()))
            count += 1
            if count >= cap:
                break
        return best

    def density(self):
        n = self.n
        if n < 2:
            return 0.0
        return 2.0 * self.m / (n * (n - 1))

    def mcc(self, max_clique=12, max_nodes=2000):
        """MCC 中心性（CytoHubba 定义）：对每个节点，累加所有包含它的极大团的 (|C|-1)!。

        实现：Bron–Kerbosch（带枢轴）枚举极大团；团大小超过 max_clique 时按 max_clique 截断，
        并返回截断信息供调用方如实报告（截断会低估大团节点的 MCC）。

        返回 (mcc, info)，info = {"cliques": n, "max_size": k, "truncated": m}
        """
        nodes = sorted(self.nodes)
        if len(nodes) > max_nodes:
            # 超大图：退化为度（并如实标注），避免组合爆炸
            return ({u: float(len(self.adj[u])) for u in nodes},
                    {"cliques": 0, "max_size": 0, "truncated": -1, "degraded": True})

        cliques = []
        truncated = 0
        adj = {u: set(self.adj[u]) for u in nodes}
        order = {u: i for i, u in enumerate(nodes)}

        def bk(R, P, X):
            nonlocal truncated
            if not P and not X:
                if len(R) > max_clique:
                    truncated += 1
                    cliques.append(R[:max_clique])
                else:
                    cliques.append(list(R))
                return
            if len(R) >= max_clique:
                # 达到上限：不再深入，记录一次截断
                if P:
                    truncated += 1
                    cliques.append(list(R))
                return
            pivot = max(sorted(P | X), key=lambda u: len(adj[u] & P)) if (P or X) else None
            for v in list(P - (adj[pivot] if pivot else set())):
                bk(R + [v], P & adj[v], X & adj[v])
                P = P - {v}
                X = X | {v}

        bk([], set(nodes), set())

        mcc = {u: 0.0 for u in nodes}
        for cl in cliques:
            w = 1
            for i in range(2, len(cl) + 1):
                w *= (i - 1)
            for v in cl:
                mcc[v] += float(w)
        info = {"cliques": len(cliques),
                "max_size": max((len(c) for c in cliques), default=0),
                "truncated": truncated, "degraded": False}
        return mcc, info


    def louvain(self, resolution=1.0, seed=42):
        """多层 Louvain（局部移动 + 聚合），返回互斥社区列表（原始节点标签）。

        边权重参与计算；无权图（w=1）等价于经典无权 Louvain。随机顺序固定种子。
        """
        rng = random.Random(seed)
        adj = {u: dict(v) for u, v in self.adj.items()}
        members = {u: [u] for u in adj}
        comm_of = {u: u for u in adj}
        total_w = sum(sum(nb.values()) for nb in adj.values()) / 2.0
        if total_w <= 0:
            return [[u] for u in adj]
        m2 = 2.0 * total_w

        for _level in range(12):
            tot = defaultdict(float)
            for u in adj:
                tot[comm_of[u]] += sum(adj[u].values())
            order = list(adj)
            rng.shuffle(order)
            for _pass in range(60):
                moved = 0
                for u in order:
                    cu = comm_of[u]
                    k_u = sum(adj[u].values())
                    if k_u <= 0:
                        continue
                    nb_comm = defaultdict(float)
                    for v, w in adj[u].items():
                        nb_comm[comm_of[v]] += w
                    tot[cu] -= k_u
                    best_c, best_gain = cu, nb_comm.get(cu, 0.0) - resolution * tot[cu] * k_u / m2
                    for c, w in nb_comm.items():
                        gain = w - resolution * tot[c] * k_u / m2
                        if gain > best_gain + 1e-10:
                            best_c, best_gain = c, gain
                    tot[best_c] += k_u
                    if best_c != cu:
                        comm_of[u] = best_c
                        moved += 1
                if moved == 0:
                    break

            groups = defaultdict(list)
            for u, c in comm_of.items():
                groups[c].append(u)
            if len(groups) == len(adj):
                return sorted((sorted(x for m_ in mem for x in members[m_]) for mem in groups.values()),
                              key=len, reverse=True)

            new_adj, new_members = {}, {}
            for c, mem in groups.items():
                new_members[c] = [x for m_ in mem for x in members[m_]]
                nb = defaultdict(float)
                inner = 0.0
                mem_set = set(mem)
                for u in mem:
                    for v, w in adj[u].items():
                        cv = comm_of[v]
                        if cv == c:
                            if v in mem_set and u <= v:
                                inner += w if u != v else w / 2.0
                        else:
                            nb[cv] += w
                # 自环权重计入聚合节点的度，否则 k 被低估、模块度增益虚高
                nb[c] += 2.0 * inner
                new_adj[c] = dict(nb)
            adj, members = new_adj, new_members
            comm_of = {c: c for c in adj}

        return sorted((sorted(mem) for mem in members.values()), key=len, reverse=True)

    def modularity(self, groups, resolution=1.0):
        total_w = sum(sum(nb.values()) for nb in self.adj.values()) / 2.0
        if total_w <= 0:
            return 0.0
        k = {u: sum(v.values()) for u, v in self.adj.items()}
        q = 0.0
        for g in groups:
            gs = set(g)
            lc = sum(self.adj[u].get(v, 0.0) for u in gs for v in gs) / 2.0
            dc = sum(k.get(u, 0.0) for u in gs)
            q += lc / total_w - resolution * (dc / (2.0 * total_w)) ** 2
        return q

    def layout_fr(self, iters=80, seed=42, width=900, height=640):
        """Fruchterman-Reingold 力导向布局。"""
        rng = random.Random(seed)
        nodes = sorted(self.nodes)
        n = len(nodes)
        if n == 0:
            return {}
        if n == 1:
            return {nodes[0]: (width / 2, height / 2)}
        pos = {v: (rng.uniform(0, width), rng.uniform(0, height)) for v in nodes}
        area = width * height
        k = math.sqrt(area / n)
        t = width / 10.0
        for _ in range(iters):
            disp = {v: [0.0, 0.0] for v in nodes}
            for i in range(n):
                u = nodes[i]
                ux, uy = pos[u]
                for j in range(i + 1, n):
                    v = nodes[j]
                    vx, vy = pos[v]
                    dx, dy = ux - vx, uy - vy
                    dist = math.hypot(dx, dy) or 0.01
                    force = k * k / dist
                    ux2, uy2 = dx / dist * force, dy / dist * force
                    disp[u][0] += ux2
                    disp[u][1] += uy2
                    disp[v][0] -= ux2
                    disp[v][1] -= uy2
            for u in nodes:
                for v in self.adj[u]:
                    if u >= v:
                        continue
                    ux, uy = pos[u]
                    vx, vy = pos[v]
                    dx, dy = ux - vx, uy - vy
                    dist = math.hypot(dx, dy) or 0.01
                    force = dist * dist / k
                    fx, fy = dx / dist * force, dy / dist * force
                    disp[u][0] -= fx
                    disp[u][1] -= fy
                    disp[v][0] += fx
                    disp[v][1] += fy
            for u in nodes:
                dx, dy = disp[u]
                dist = math.hypot(dx, dy) or 0.01
                step = min(dist, t)
                x = pos[u][0] + dx / dist * step
                y = pos[u][1] + dy / dist * step
                pos[u] = (min(width - 30, max(30, x)), min(height - 30, max(30, y)))
            t *= 0.95
        return pos


# ----------------------------------------------------------------------------
# svg
# ----------------------------------------------------------------------------

PALETTE = ["#2F6FB2", "#C0504D", "#4F8A5B", "#8064A2", "#E3A33B", "#3E9BA6",
           "#8C6D46", "#B65A8C", "#5B7C99", "#7E9E4E"]


def _esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;"))


def svg_doc(body, w, h, title=""):
    head = ('<?xml version="1.0" encoding="UTF-8"?>\n'
            '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" '
            'font-family="-apple-system,Segoe UI,Microsoft YaHei,SimHei,sans-serif">\n'
            '<rect width="100%%" height="100%%" fill="#ffffff"/>\n') % (w, h, w, h)
    t = ('<text x="24" y="26" font-size="16" font-weight="600" fill="#1f2937">%s</text>'
         % _esc(title)) if title else ""
    return head + t + "\n" + body + "\n</svg>\n"


def bar_svg(items, path, title="", xlabel="", width=720, value_fmt="%.3g", color="#2F6FB2"):
    """items: [(label, value)]，横向条形图。"""
    items = [(str(k), float(v)) for k, v in items]
    if not items:
        write_text(path, svg_doc('<text x="24" y="60" font-size="13" fill="#6b7280">无数据</text>',
                                 width, 120, title))
        return path
    row_h, pad_l, pad_t = 24, 260, 52
    h = pad_t + row_h * len(items) + 40
    vmax = max(v for _, v in items) or 1.0
    body = []
    plot_w = width - pad_l - 60
    for i, (lab, val) in enumerate(items):
        y = pad_t + i * row_h
        w = max(1.5, plot_w * val / vmax)
        body.append('<rect x="%d" y="%d" width="%.1f" height="%d" fill="%s" rx="3"/>'
                    % (pad_l, y, w, row_h - 8, color))
        body.append('<text x="%d" y="%d" font-size="12" fill="#374151" text-anchor="end">%s</text>'
                    % (pad_l - 8, y + 12, _esc(lab[:38])))
        body.append('<text x="%.1f" y="%d" font-size="11" fill="#4b5563">%s</text>'
                    % (pad_l + w + 6, y + 12, value_fmt % val))
    if xlabel:
        body.append('<text x="%d" y="%d" font-size="12" fill="#6b7280">%s</text>'
                    % (pad_l, h - 12, _esc(xlabel)))
    write_text(path, svg_doc("\n".join(body), width, h, title))
    return path


def bubble_svg(items, path, title="", width=760, top=20):
    """items: [(term, gene_ratio, neglog10q, count)]，气泡图。"""
    items = items[:top]
    if not items:
        write_text(path, svg_doc('<text x="24" y="60" font-size="13" fill="#6b7280">无显著条目</text>',
                                 width, 120, title))
        return path
    pad_l, pad_t, row_h = 330, 56, 26
    h = pad_t + row_h * len(items) + 46
    rmax = max(i[3] for i in items) or 1
    nlmax = max(i[2] for i in items) or 1.0
    grmax = max(i[1] for i in items) or 1.0
    body = []
    plot_w = width - pad_l - 70
    for i, (term, gr, nl, cnt) in enumerate(items):
        y = pad_t + i * row_h + row_h / 2
        x = pad_l + plot_w * (gr / grmax)
        r = 6 + 12 * math.sqrt(cnt / rmax)
        t = nl / nlmax
        col = "#%02x%02x%02x" % (int(214 - 120 * t), int(90 + 60 * t), int(70 + 90 * t))
        body.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" fill-opacity="0.72" stroke="#ffffff"/>'
                    % (x, y, r, col))
        body.append('<text x="%d" y="%.1f" font-size="11.5" fill="#374151" text-anchor="end">%s</text>'
                    % (pad_l - 10, y + 4, _esc(str(term)[:52])))
    body.append('<text x="%d" y="%d" font-size="11" fill="#6b7280">气泡大小=命中基因数；颜色越红=−log10(q) 越大</text>'
                % (pad_l, h - 14))
    body.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#d1d5db"/>'
                % (pad_l, h - 34, pad_l + plot_w, h - 34))
    write_text(path, svg_doc("\n".join(body), width, h, title))
    return path


def line_svg(series, path, title="", xlabel="", ylabel="", width=720, height=340):
    """折线/ECDF 图。series = [(label, [(x, y), ...])]。"""
    pad_l, pad_r, pad_t, pad_b = 64, 24, 52, 52
    body = []
    xs = [x for _, pts in series for x, _ in pts] or [0, 1]
    ys = [y for _, pts in series for _, y in pts] or [0, 1]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if xmax == xmin:
        xmax = xmin + 1
    if ymax == ymin:
        ymax = ymin + 1
    pw, ph = width - pad_l - pad_r, height - pad_t - pad_b

    def sx(x):
        return pad_l + pw * (x - xmin) / (xmax - xmin)

    def sy(y):
        return pad_t + ph * (1 - (y - ymin) / (ymax - ymin))

    body.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#e5e7eb"/>'
                % (pad_l, pad_t, pw, ph))
    for i in range(5):
        gy = pad_t + ph * i / 4.0
        val = ymax - (ymax - ymin) * i / 4.0
        body.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#f3f4f6"/>'
                    % (pad_l, gy, pad_l + pw, gy))
        body.append('<text x="%d" y="%.1f" font-size="10" fill="#9ca3af" text-anchor="end">%.3g</text>'
                    % (pad_l - 6, gy + 3, val))
    for i, (lab, pts) in enumerate(series):
        col = PALETTE[i % len(PALETTE)]
        d = " ".join(("M" if j == 0 else "L") + "%.1f,%.1f" % (sx(x), sy(y))
                     for j, (x, y) in enumerate(pts))
        body.append('<path d="%s" fill="none" stroke="%s" stroke-width="1.8"/>' % (d, col))
        body.append('<rect x="%d" y="%d" width="10" height="10" fill="%s"/>'
                    % (pad_l + pw - 150, pad_t - 22 + 16 * i, col))
        body.append('<text x="%d" y="%d" font-size="11" fill="#374151">%s</text>'
                    % (pad_l + pw - 134, pad_t - 13 + 16 * i, _esc(lab)))
    body.append('<text x="%d" y="%d" font-size="11" fill="#6b7280" text-anchor="middle">%s</text>'
                % (pad_l + pw / 2, height - 16, _esc(xlabel)))
    body.append('<text x="16" y="%d" font-size="11" fill="#6b7280">%s</text>'
                % (pad_t + ph / 2, _esc(ylabel)))
    write_text(path, svg_doc("\n".join(body), width, height, title))
    return path


def edges_from_rows(rows, gene_col_a, gene_col_b, score_col=None, min_score=None):
    """从表格行构造 Graph（统一清洗：去空、去自环、重复边取最大权重）。"""
    g = Graph()
    for r in rows:
        u = (r.get(gene_col_a) or "").strip().upper()
        v = (r.get(gene_col_b) or "").strip().upper()
        if not u or not v:
            continue
        w = 1.0
        if score_col:
            sc = to_float(r.get(score_col), None)
            if sc is not None:
                if sc > 1.0:
                    sc = sc / 1000.0
                if min_score is not None and sc < min_score:
                    continue
                w = sc
        g.add_edge(u, v, w)
    return g


def degree_binned_z(values, degrees, n_bins=10, min_bin=5):
    """把指标按**度分位分箱**后组内标准化，返回 (z, info)。

    为什么不能用"度 ±1"窗口：PPI 的度分布是长尾，281 个节点散布在 100 多个度值上，
    窄窗口里常常只有 1–2 个同伴，标准差无定义，于是 z 恒为 0 —— 看起来"没有桥接节点"，
    实际是"这个度量根本没算出来"（0 与 NA 的语义差别足以误导结论）。
    改为等频分箱（每箱 ≥min_bin 个节点），箱内标准化；箱内方差为 0 时返回 None（NA）。
    """
    vals = {k: float(v) for k, v in values.items()}
    order = sorted(vals, key=lambda k: (degrees.get(k, 0), k))
    n = len(order)
    if n == 0:
        return {}, {"bins": []}
    k = max(1, min(n_bins, n // max(1, min_bin)))
    size = max(1, n // k)
    bins = [order[i:i + size] for i in range(0, n, size)]
    z, info = {}, {"bins": [], "n_bins": len(bins)}
    for b in bins:
        arr = [vals[u] for u in b]
        ds = [degrees.get(u, 0) for u in b]
        mu = sum(arr) / len(arr)
        var = sum((x - mu) ** 2 for x in arr) / (len(arr) - 1) if len(arr) > 1 else 0.0
        sd = var ** 0.5
        info["bins"].append({"n": len(b), "deg_min": min(ds), "deg_max": max(ds),
                             "mean": round(mu, 6), "sd": round(sd, 6)})
        for u in b:
            z[u] = ((vals[u] - mu) / sd) if sd > 1e-12 and len(arr) >= min_bin else None
    return z, info


def network_svg(g, path, title="", node_score=None, groups=None, label_top=15, width=920, height=680,
                seed=42):
    """网络图：节点大小∝度数，颜色按社区，标签只标 top。"""
    if g.n == 0:
        write_text(path, svg_doc('<text x="24" y="60" font-size="13" fill="#6b7280">无节点</text>',
                                 width, 120, title))
        return path
    pos = g.layout_fr(iters=70, seed=seed, width=width - 60, height=height - 70)
    pos = {k: (v[0] + 30, v[1] + 46) for k, v in pos.items()}
    deg = g.degree()
    dmax = max(deg.values()) or 1
    grp = {}
    if groups:
        for i, gset in enumerate(groups):
            for v in gset:
                grp[v] = i
    else:
        for v in g.nodes:
            grp[v] = 0
    body = []
    for u in sorted(g.nodes):
        for v, w in sorted(g.adj[u].items()):
            if u < v:
                x1, y1 = pos[u]
                x2, y2 = pos[v]
                body.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#cbd5e1" '
                            'stroke-width="0.6" stroke-opacity="0.75"/>' % (x1, y1, x2, y2))
    order = sorted(g.nodes, key=lambda x: (deg[x], x))
    for u in order:
        x, y = pos[u]
        r = 3.0 + 11.0 * math.sqrt(deg[u] / dmax)
        col = PALETTE[grp.get(u, 0) % len(PALETTE)]
        body.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" fill-opacity="0.85" stroke="#ffffff" '
                    'stroke-width="0.8"/>' % (x, y, r, col))
    if label_top:
        for u in sorted(g.nodes, key=lambda x: (-deg[x], x))[:label_top]:
            x, y = pos[u]
            body.append('<text x="%.1f" y="%.1f" font-size="10.5" fill="#111827" '
                        'text-anchor="middle">%s</text>' % (x, y - 9, _esc(u)))
    body.append('<text x="24" y="%d" font-size="11" fill="#6b7280">节点数 %d · 边数 %d · 节点大小∝度 · '
                '颜色=模块</text>' % (height - 12, g.n, g.m))
    write_text(path, svg_doc("\n".join(body), width, height, title))
    return path


def venn_svg(sets, path, title="", width=680, height=460):
    """2–4 集合韦恩图（椭圆近似；>4 集合退化为条形）。"""
    names = list(sets.keys())
    k = len(names)
    if k < 2:
        write_text(path, svg_doc('<text x="24" y="60" font-size="13" fill="#6b7280">'
                                 '至少需要 2 组</text>', width, 120, title))
        return path
    if k > 4:
        items = [(n, len(sets[n])) for n in names]
        return bar_svg(items, path, title=title + "（集合规模）", xlabel="基因数")
    body = []
    cx, cy = width / 2, height / 2 + 10
    if k == 2:
        r = 118
        centers = [(cx - 62, cy), (cx + 62, cy)]
    elif k == 3:
        r = 120
        centers = [(cx, cy - 58), (cx - 88, cy + 52), (cx + 88, cy + 52)]
    else:
        r = 112
        centers = [(cx - 60, cy - 52), (cx + 60, cy - 52), (cx - 60, cy + 52), (cx + 60, cy + 52)]
    for i, (x, y) in enumerate(centers):
        body.append('<circle cx="%.1f" cy="%.1f" r="%d" fill="%s" fill-opacity="0.28" stroke="%s" '
                    'stroke-width="1.2"/>' % (x, y, r, PALETTE[i], PALETTE[i]))
    # 计数：每个交集区域
    from itertools import combinations
    region_labels = []
    for size in range(1, k + 1):
        for combo in combinations(range(k), size):
            key = set(names[i] for i in combo)
            inter = set.intersection(*[sets[names[i]] for i in combo]) if combo else set()
            others = [sets[names[i]] for i in range(k) if i not in combo]
            val = len(inter - set().union(*others)) if others else len(inter)
            if size == 1:
                px = centers[combo[0]][0] + (0 if k > 2 else (-30 if combo[0] == 0 else 30))
                py = centers[combo[0]][1] + (0 if k < 3 else (-62 if combo[0] == 0 else 66))
            else:
                px = sum(centers[i][0] for i in combo) / size
                py = sum(centers[i][1] for i in combo) / size
            region_labels.append((px, py, val, size))
    for px, py, val, size in region_labels:
        body.append('<text x="%.1f" y="%.1f" font-size="%d" font-weight="600" fill="#111827" '
                    'text-anchor="middle">%d</text>' % (px, py, 15 if size == 1 else 13, val))
    for i, (x, y) in enumerate(centers):
        ly = y - r - 10 if k > 1 and i == 0 else y + r + 18
        lx = x
        if k == 2:
            ly = y - r - 12
            lx = x + (-10 if i == 0 else 10)
        body.append('<text x="%.1f" y="%.1f" font-size="12" fill="%s" text-anchor="middle" '
                    'font-weight="600">%s (%d)</text>'
                    % (lx, ly, PALETTE[i], _esc(names[i]), len(sets[names[i]])))
    write_text(path, svg_doc("\n".join(body), width, height, title))
    return path


# ----------------------------------------------------------------------------
# prov 溯源
# ----------------------------------------------------------------------------

def sha256_file(path, limit=None):
    """全量哈希（默认）。limit 仅用于显式声明"部分哈希"的场景。"""
    try:
        h = hashlib.sha256()
        read = 0
        with open(path, "rb") as f:
            while True:
                b = f.read(1 << 20)
                if not b:
                    break
                h.update(b)
                read += len(b)
                if limit and read >= limit:
                    break
        return h.hexdigest()[:16]
    except Exception:
        return "NA"


def environment_snapshot():
    """可复现所需的最小环境快照。换机排障时先看这个。"""
    import platform as _p
    import locale
    snap = {
        "python": sys.version.split()[0],
        "implementation": _p.python_implementation(),
        "platform": _p.platform(),
        "machine": _p.machine(),
        "locale": ".".join(str(x) for x in locale.getlocale()),
        "filesystem_encoding": sys.getfilesystemencoding(),
        "stdout_encoding": (getattr(sys.stdout, "encoding", "") or ""),
        "cwd": os.getcwd(),
        "tz": time.strftime("%Z"),
        "packages": {},
    }
    for mod in ("numpy", "scipy", "openpyxl", "networkx", "matplotlib", "rdkit"):
        try:
            m = __import__(mod)
            snap["packages"][mod] = getattr(m, "__version__", "unknown")
        except Exception:
            snap["packages"][mod] = None
    return snap


def build_claims(outdir, label="分析"):
    """从阶段产物生成结论账本 claims.json（结论↔数值↔来源文件↔参数来源）。

    为什么由核心库提供：报告与审计都要用，且必须在审计**之前**生成，
    否则"存在账本"这条结构性规则会永远 FAIL（顺序依赖 bug）。
    """
    man_p = os.path.join(outdir, "run_manifest.json")
    stages = []
    if os.path.exists(man_p):
        try:
            stages = json.load(open(man_p, encoding="utf-8")).get("stages", [])
        except Exception:
            stages = []

    def param(prefix, key):
        for st in stages:
            if str(st.get("stage", "")).startswith(prefix):
                pr = st.get("params") or {}
                if key in pr and pr[key] is not None:
                    return pr[key], st.get("stage")
        return None, None

    counts = {}
    for st in stages:
        counts.update(st.get("counts", {}) or {})

    def nrows(name):
        pth = os.path.join(outdir, name)
        if not os.path.exists(pth):
            return None
        try:
            _h, rows = read_table(pth)
            return len(rows)
        except Exception:
            return None

    claims, missing = [], []

    def need(label, prefix, key):
        v, src = param(prefix, key)
        if v is None:
            missing.append("%s（阶段 %s 未记录 %s）" % (label, prefix, key))
        return v, src

    ob, s2 = need("OB 阈值", "02", "ob")
    dl, _ = need("DL 阈值", "02", "dl")
    route, _ = need("给药途径", "02", "route")
    ms, s4 = need("PPI 置信度阈值", "04", "min_score")
    bs, _ = need("bootstrap 次数", "04", "bootstrap")
    qv, s6 = need("富集 q 阈值", "06", "q")

    kept = nrows("compounds_kept.csv")
    if kept is not None:
        # 结论文字必须与实际门控一致：用类药性门时写 OB≥30 等于伪造方法
        gate = param("02", "gate")[0]
        n_obdl = (counts or {}).get("gate_obdl", 0) or 0
        n_drug = (counts or {}).get("gate_druglike", 0) or 0
        if gate in ("auto", "druglike") and n_drug and not n_obdl:
            desc = "门控=类药性（QED≥%s，Veber 违例≤%s）" % (param("02", "qed_min")[0],
                                                          param("02", "max_veber")[0])
        elif n_obdl and n_drug:
            desc = "门控=混合（OB/DL %d 个、类药性 %d 个）" % (n_obdl, n_drug)
        elif gate == "evidence":
            desc = "门控=证据（按实测活性/入血证据，类药性仅标记）"
        else:
            desc = "门控=OB≥%s, DL≥%s" % (ob, dl)
        claims.append({"id": "C1", "text": "纳入 %d 个成分；%s；route=%s"
                       % (kept, desc, route), "level": "L1",
                       "value": {"kept": kept, "gate": gate, "ob": ob, "dl": dl, "route": route},
                       "source": "compounds_kept.csv", "params_ref": s2})
    nint = nrows("intersection.csv")
    if nint is not None:
        claims.append({"id": "C2", "text": "药物靶点与疾病靶点交集 %d 个" % nint, "level": "L1",
                       "value": {"intersection": nint}, "source": "intersection.csv",
                       "params_ref": "03_intersect"})
    nh = nrows("hubs.csv")
    if nh is not None:
        claims.append({"id": "C3", "text": "hub 候选 %d 个（8 指标秩均值共识，STRING≥%s，bootstrap %s 次）"
                       % (nh, ms, bs), "level": "L1",
                       "value": {"hubs": nh, "min_score": ms, "bootstrap": bs},
                       "source": "hubs.csv", "params_ref": s4})
    ne = nrows("enrich_significant.csv")
    if ne is not None:
        mech = None
        try:
            _h, rows = read_table(os.path.join(outdir, "enrich_significant.csv"))
            mech = sum(1 for r in rows if r.get("tier") == "机制明确")
        except Exception:
            pass
        claims.append({"id": "C4", "text": "显著富集 %d 条（q<%s），其中机制明确类 %s 条"
                       % (ne, qv, mech), "level": "L1",
                       "value": {"sig": ne, "mechanistic": mech, "q": qv},
                       "source": "enrich_significant.csv", "params_ref": s6})
    nd = nrows("null_tests.csv")
    if nd is not None:
        claims.append({"id": "C5", "text": "构建 %d 套零模型" % nd, "level": "L1",
                       "value": {"nulls": nd}, "source": "null_tests.csv",
                       "params_ref": "05_null"})

    data = {"label": label, "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "claims": claims, "missing_params": missing}
    os.makedirs(outdir, exist_ok=True)
    json.dump(data, open(os.path.join(outdir, "claims.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    return data


class Manifest:
    """阶段级溯源清单：把"参数—版本—时间—输入—输出—种子"写进 JSON。"""

    def __init__(self, stage, outdir, params=None):
        self.data = {
            "stage": stage,
            "tool": "network-pharmacology-pro",
            "tool_version": VERSION,
            "python": sys.version.split()[0],
            "platform": sys.platform,
            "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "argv": list(sys.argv),
            "env": environment_snapshot(),
            "params": params or {},
            "inputs": {},
            "outputs": {},
            "warnings": [],
            "counts": {},
        }
        self.outdir = outdir

    def input(self, label, path):
        if path and os.path.exists(str(path)):
            self.data["inputs"][label] = {
                "path": os.path.abspath(str(path)),
                "sha256_16": sha256_file(str(path)),
                "bytes": os.path.getsize(str(path)),
                "hash": "full",
            }
        elif path:
            self.data["inputs"][label] = {"path": str(path), "missing": True}

    def output(self, label, path):
        if path and os.path.exists(str(path)):
            self.data["outputs"][label] = {
                "path": os.path.abspath(str(path)),
                "sha256_16": sha256_file(str(path)),
                "bytes": os.path.getsize(str(path)),
                "hash": "full",
            }

    def count(self, key, val):
        self.data["counts"][key] = val

    def warn(self, msg):
        self.data["warnings"].append(msg)

    def save(self, name="run_manifest.json"):
        self.data["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        path = os.path.join(self.outdir, name)
        os.makedirs(self.outdir, exist_ok=True)
        prev = {}
        if os.path.exists(path):
            try:
                prev = json.load(open(path, encoding="utf-8"))
            except Exception:
                prev = {}
        stages = prev.get("stages", []) if isinstance(prev, dict) else []
        stages = [s for s in stages if s.get("stage") != self.data["stage"]] + [self.data]
        json.dump({"stages": stages}, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
        # 环境快照单独落一份，换机排障时先看它
        envp = os.path.join(self.outdir, "environment.json")
        if not os.path.exists(envp):
            json.dump(environment_snapshot(), open(envp, "w", encoding="utf-8"),
                      ensure_ascii=False, indent=2)
        return path


# ----------------------------------------------------------------------------
# cli
# ----------------------------------------------------------------------------

class Cli(argparse.ArgumentParser):
    def error(self, message):
        self.print_usage(sys.stderr)
        sys.stderr.write("\n[参数错误] %s\n提示：用 --help 查看可用参数。\n" % message)
        raise SystemExit(2)



def tool_path(env_var, *candidates, exe=None):
    """工具解析：环境变量优先 → 显式候选 → PATH；不把作者机器路径写死。

    env_var 如 NP_VINA；candidates 可含历史路径（兼容旧环境），但不作为唯一来源。
    返回绝对路径或 None。
    """
    v = os.environ.get(env_var)
    if v and os.path.exists(v):
        return v
    for c in candidates:
        if c and os.path.exists(c):
            return c
    for name in ([exe] if exe else []):
        w = shutil.which(name)
        if w:
            return w
    return None


def cache_dir(kind, root=None):
    """跨运行共享的缓存目录：$NP_CACHE/<kind> 或 ~/.np_cache/<kind>。

    为什么不能放在输出目录里：网络响应（STRING/Ensembl/HPA）如果跟着每次运行重建，
    重跑就要重新联网，既慢（实测 23s → 423s）又不可复现。缓存集中后，
    第二次运行可完全离线、且结果哈希一致。
    """
    base = (root or os.environ.get("NP_SNAPSHOT") or os.environ.get("NP_CACHE")
            or os.path.join(os.path.expanduser("~"), ".np_cache"))
    p = os.path.join(base, kind)
    try:
        os.makedirs(p, exist_ok=True)
    except Exception:
        return None
    return p


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def gene_upper(x):
    return (x or "").strip().upper()


def split_multi(cell):
    if not cell:
        return []
    parts = re.split(r"[;,/|、\s]+", str(cell))
    return [p.strip() for p in parts if p.strip()]


def read_gene_list(path, warn=None):
    """支持一行一基因、一行多基因、带表头 csv/tsv。

    输入契约（此前未写明，导致"A/B/C 读成空列表"被当成 bug —— 评审 §6）：
      * 接受：`IL6`（一行一基因）、`IL6 TNF AKT1`（一行多基因，空格/逗号/分号分隔）、
        带表头 csv/tsv（任意列值只要形如基因符号都会被收集）。
      * **丢弃**：单字符 token（如 `A`）、纯数字、NA/NAN/NONE/NULL/GENE/SYMBOL 占位符。
        理由：没有 1 字符的人类基因符号，这些几乎总是表头片段或占位符。
      * 需要保留单字符 token（例如自定义 ID）时，请改用 `read_table()` 直接读，别用本函数。
      * `warn=True` 或环境变量 `NP_WARN_INPUT=1` 时，把丢弃项打到 stderr。
    """
    text, _ = _read_text(str(path))
    genes = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "\t" in line or "," in line:
            _h, rows = read_table(str(path))
            for r in rows:
                for v in r.values():
                    if v and re.fullmatch(r"[A-Za-z0-9\-\._@]{1,20}", v.strip()):
                        genes.append(v.strip())
            break
        genes.extend(split_multi(line))
    out, seen, dropped = [], set(), []
    for g in genes:
        gu = gene_upper(g)
        if (gu and gu not in seen and len(gu) >= 2
                and not re.fullmatch(r"\d+", gu)
                and gu not in ("NA", "NAN", "NONE", "NULL", "GENE", "SYMBOL")):
            seen.add(gu)
            out.append(gu)
        elif gu:
            dropped.append(gu)
    if dropped and (warn or os.environ.get("NP_WARN_INPUT") == "1"):
        sys.stderr.write("[input] %s：丢弃 %d 个 token（单字符/纯数字/占位符）：%s%s\n" % (
            os.path.basename(str(path)), len(dropped), ",".join(dropped[:10]),
            "…" if len(dropped) > 10 else ""))
    return out
