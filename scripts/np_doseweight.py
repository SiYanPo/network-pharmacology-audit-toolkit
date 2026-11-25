#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_doseweight.py —— 剂量加权拓扑模块（传统 NP 的创新点：让"10g vs 50g 黄芪"产生不同网络）

问题（审稿人视角）
    传统网络药理学把每个活性成分等权或仅按 degree 加权，完全忽略处方剂量。
    结果是"某成分是否在方中"决定一切，而"用了 10 g 还是 50 g"毫无区别——这在药理学上不成立。

方案（暴露量加权的成分-靶点网络）
    每个成分 i 的剂量权重 w_i：
        w_i = daily_g × content_i(%) × extract(%) × F_i × fu_i   （来自 np_dose.chain 的日暴露量）
    边权 w(i→target) = w_i（成分对靶点的作用强度按暴露量缩放）
    节点度（剂量加权度）= Σ 边权（而非简单计数）
    核心靶点重排：用剂量加权度/介数替代无权度 → 得到"在给定剂量下真正主导网络的靶点"

    两剂量对比模式：给 --dose-a / --dose-b 两个处方剂量表，输出 hub 排名漂移，
    直接回答"剂量变化时核心靶点如何改变"，可作论文 Figures。

诚实边界
    * content/F/fu 多为文献估计或默认值，须在 meta 声明为假设（np_dose 已做）。
    * 权重是"暴露量代理"，不等于真实药效强度；报告中须写明是 dose-exposure weighting。

用法
    python np_doseweight.py --edges edges.tsv --dose dose.csv -o out/
    python np_doseweight.py --edges edges.tsv --dose doseA.csv --dose-b doseB.csv -o out/
    python np_doseweight.py --selftest
    剂量表列: compound, daily_g, content_pct, [extract_pct, f_pct, fu, doses, mw]
"""
import argparse
import csv
import json
import os
import sys

try:
    import networkx as nx
except ImportError:
    nx = None

HERE = os.path.dirname(os.path.abspath(__file__))


def read_csv(path, delim="\t"):
    with open(path, encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f, delimiter=delim)
        return list(rdr), (rdr.fieldnames or [])


def exposure(daily_g, content_pct, extract_pct=50.0, f_pct=30.0, fu=0.1, doses=2, mw=300.0):
    """与 np_dose.chain 一致的暴露量（相对值，用于加权，单位 µM 量级）。"""
    intake = daily_g * 1000.0 * (content_pct / 100.0) * (extract_pct / 100.0)
    per = intake / max(1, doses)
    cmax = (f_pct / 100.0) * per * 1000.0 / max(1e-9, mw * doses)
    return cmax * fu


def load_dose(path):
    rows, hdrs = read_csv(path, ",")
    d = {}
    for r in rows:
        c = (r.get("compound") or r.get("name") or r.get("成分") or "").strip()
        if not c:
            continue
        try:
            daily = float(r.get("daily_g") or r.get("剂量") or 0)
            content = float(r.get("content_pct") or r.get("含量") or 0)
        except ValueError:
            continue
        d[c.lower()] = {
            "daily_g": daily, "content_pct": content,
            "extract_pct": float(r.get("extract_pct") or 50),
            "f_pct": float(r.get("f_pct") or 30),
            "fu": float(r.get("fu") or 0.1),
            "doses": int(float(r.get("doses") or 2)),
            "mw": float(r.get("mw") or 300),
        }
    return d


def build_weighted_graph(edges_path, dose, src_col=None, dst_col=None):
    """成分-靶点边表 → 剂量加权图。返回 (G, unmatched_compounds)。"""
    rows, hdrs = read_csv(edges_path, "\t" if edges_path.lower().endswith(".tsv") else ",")
    sc = src_col or next((h for h in hdrs if h.lower() in ("compound", "src_id", "成分", "source")), hdrs[0] if hdrs else None)
    dc = dst_col or next((h for h in hdrs if h.lower() in ("gene", "dst_id", "target", "靶点")), hdrs[1] if len(hdrs) > 1 else None)
    if nx is None:
        raise RuntimeError("需 networkx：pip install networkx")
    G = nx.Graph()
    unmatched = set()
    for r in rows:
        comp = (r.get(sc) or "").strip()
        tgt = (r.get(dc) or "").strip().upper()
        if not comp or not tgt:
            continue
        pd = dose.get(comp.lower())
        if pd is None:
            unmatched.add(comp)
            w = 0.0  # 无剂量信息的成分不贡献权重（但保留节点，透明）
        else:
            w = exposure(pd["daily_g"], pd["content_pct"], pd["extract_pct"],
                         pd["f_pct"], pd["fu"], pd["doses"], pd["mw"])
        if G.has_edge(comp, tgt):
            G[comp][tgt]["weight"] += w
        else:
            G.add_edge(comp, tgt, weight=w)
    return G, unmatched


def hub_ranking(G, top=15, weighted=True):
    """按剂量加权度给靶点排名（靶点 = 图中度>0 的非成分节点近似：取 dst 侧）。"""
    # 靶点集合：出现在边表 dst 的节点；这里用连通关系近似——权重总和最高的节点即靶点侧
    scored = []
    for n in G.nodes():
        wd = sum(d.get("weight", 0) for _, _, d in G.edges(n, data=True))
        deg = G.degree(n)
        scored.append((n, wd, deg))
    scored.sort(key=lambda x: (-x[1], -x[2]))
    return scored[:top]


def cmd_compare(args):
    doseA = load_dose(args.dose)
    doseB = load_dose(args.dose_b)
    GA, una = build_weighted_graph(args.edges, doseA)
    GB, unb = build_weighted_graph(args.edges, doseB)
    hA = hub_ranking(GA, args.top)
    hB = hub_ranking(GB, args.top)
    rankA = {n: i + 1 for i, (n, _, _) in enumerate(hA)}
    rankB = {n: i + 1 for i, (n, _, _) in enumerate(hB)}
    all_nodes = set(rankA) | set(rankB)
    drift = []
    for n in all_nodes:
        ra, rb = rankA.get(n, 999), rankB.get(n, 999)
        drift.append({"target": n, "rank_A": ra, "rank_B": rb,
                      "delta": (ra - rb) if (ra < 999 and rb < 999) else None})
    drift.sort(key=lambda x: (x["rank_A"] if x["rank_A"] < 999 else 999))
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "dose_hub_drift.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["target", "rank_A", "rank_B", "delta"])
        w.writeheader(); w.writerows(drift)
    meta = {
        "note": "权重=日暴露量代理(daily_g×content×extract×F×fu)，非药效强度；content/F/fu 多为假设",
        "doseA_file": args.dose, "doseB_file": args.dose_b,
        "unmatched_compounds_A": sorted(una), "unmatched_compounds_B": sorted(unb),
        "top_A": [{"target": n, "weight": round(w, 4)} for n, w, _ in hA[:10]],
        "top_B": [{"target": n, "weight": round(w, 4)} for n, w, _ in hB[:10]],
    }
    json.dump(meta, open(os.path.join(args.outdir, "doseweight_meta.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[doseweight] A top5: %s" % [n for n, _, _ in hA[:5]])
    print("[doseweight] B top5: %s" % [n for n, _, _ in hB[:5]])
    print("[doseweight] 排名漂移 → %s" % os.path.join(args.outdir, "dose_hub_drift.csv"))
    return 0


def cmd_single(args):
    dose = load_dose(args.dose)
    G, unmatched = build_weighted_graph(args.edges, dose)
    hubs = hub_ranking(G, args.top)
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "dose_weighted_hubs.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "node", "dose_weighted_degree", "unweighted_degree"])
        for i, (n, wd, deg) in enumerate(hubs, 1):
            w.writerow([i, n, round(wd, 4), deg])
    print("[doseweight] 剂量加权 hubs top10:")
    for i, (n, wd, deg) in enumerate(hubs[:10], 1):
        print("  %2d. %-10s w=%.4f (unw=%d)" % (i, n, wd, deg))
    if unmatched:
        print("[doseweight] %d 个成分无剂量信息，权重=0：%s" % (len(unmatched), sorted(unmatched)[:8]))
    return 0


def selftest():
    d = {"daily_g": 9, "content_pct": 0.5, "extract_pct": 50,
         "f_pct": 30, "fu": 0.1, "doses": 2, "mw": 302}
    e = exposure(**d)
    e10 = exposure(**{**d, "daily_g": 10})
    e50 = exposure(**{**d, "daily_g": 50})
    ok = e > 0 and abs(e50 / e10 - 5.0) < 1e-6
    print("selftest: %s (exposure@10g=%.5f, @50g=%.5f, ratio=%.1f)" %
          ("PASS" if ok else "FAIL", e10, e50, e50 / e10))
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="np_doseweight: 剂量加权拓扑与 hub 漂移")
    ap.add_argument("--edges", help="成分-靶点边表 TSV/CSV")
    ap.add_argument("--dose", help="剂量表 CSV（剂量 A）")
    ap.add_argument("--dose-b", default=None, help="剂量表 CSV（剂量 B，触发对比模式）")
    ap.add_argument("-o", "--outdir", default="np_out/doseweight")
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if not a.edges or not a.dose:
        ap.error("需 --edges 与 --dose，或 --selftest")
    return cmd_compare(a) if a.dose_b else cmd_single(a)


if __name__ == "__main__":
    sys.exit(main())
