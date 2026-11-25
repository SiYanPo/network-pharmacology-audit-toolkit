#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
selftest.py — 离线自检：算法正确性 + IO 鲁棒性 + CLI 完整性 + 端到端

设计原则
  * 断言要能证伪：数值型断言用独立参考实现（Fraction 精确枚举 / 手算值）对照，
    不做"函数能跑就算通过"的假测试。
  * 不依赖网络：全部离线。
  * 失败即报出具体断言与上下文，便于定位。

用法
  python evals/selftest.py            # 全部
  python evals/selftest.py --quick    # 跳过端到端（约 3 秒）
"""
import argparse
import csv
import math
import os
import random
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRIPTS = os.path.join(ROOT, "scripts")
sys.path.insert(0, SCRIPTS)

import npcore as C  # noqa: E402

PASS = FAIL = 0
FAILS = []


def check(cond, name, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
    else:
        FAIL += 1
        FAILS.append((name, detail))
        print("  ✗ %s %s" % (name, ("— " + detail) if detail else ""))


def approx(a, b, tol=1e-9):
    if a is None or b is None:
        return False
    return abs(a - b) <= tol * max(1.0, abs(b))


# ---------------------------------------------------------------- 统计
def ref_hypergeom_sf(k, N, K, n):
    """用 Fraction 精确枚举作参考实现。"""
    lo, hi = max(0, n - (N - K)), min(K, n)
    if k > hi:
        return 0.0
    num = sum(math.comb(K, i) * math.comb(N - K, n - i) for i in range(max(k, lo), hi + 1))
    den = math.comb(N, n)
    return float(Fraction(num, den))


def test_stats():
    print("[统计] 超几何 / BH / 经验P / Z")
    cases = [(2, 100, 10, 10), (0, 50, 5, 5), (5, 20, 8, 6), (3, 200, 40, 20), (0, 20596, 135, 37)]
    for k, N, K, n in cases:
        got = C.hypergeom_sf(k, N, K, n)
        exp = ref_hypergeom_sf(k, N, K, n)
        check(approx(got, exp, 1e-8), "hypergeom_sf(%d,%d,%d,%d)" % (k, N, K, n),
              "got %.3e exp %.3e" % (got, exp))
    check(C.hypergeom_sf(0, 100, 10, 10) == 1.0, "hypergeom 下界=1")
    check(C.hypergeom_sf(999, 100, 10, 10) == 0.0, "hypergeom 越界=0")
    check(C.hypergeom_sf(5, 10, 10, 10) == 1.0, "全命中=1")

    p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]
    q = C.bh_fdr(p)
    check(all(q[i] <= q[i + 1] + 1e-12 for i in range(len(q) - 1)), "BH 单调",
          str([round(x, 4) for x in q]))
    check(all(0 <= x <= 1 for x in q), "BH 值域")
    idx = sorted(range(len(p)), key=lambda i: p[i])
    for rank, i in enumerate(idx, 1):
        if q[i] < 1.0:
            check(q[i] <= p[i] * len(p) / rank + 1e-9, "BH 上界 rank=%d" % rank)

    e, _ = C.empirical_p(10, [1, 2, 3, 4], "greater")
    check(approx(e, 1 / 5), "empirical_p 永不为 0", str(e))
    e2, _ = C.empirical_p(5, [1, 2, 3, 4, 5], "greater")
    check(e2 == 2 / 6, "empirical_p 含等号", str(e2))
    check(C.zscore(5, [1, 2, 3, 4, 5]) > 0, "zscore 方向")
    check(math.isnan(C.zscore(5, [3])), "zscore 单样本=nan")
    check(C.jaccard({1, 2}, {2, 3}) == 1 / 3, "jaccard")
    check(C.fmt_p(0) == "<1e-300", "fmt_p 下溢")


# ---------------------------------------------------------------- 图
def test_graph():
    print("[图] 拓扑 / 中心性 / 社区 / 布局")
    g = C.Graph()
    for a, b in [("A", "B"), ("B", "C"), ("C", "D")]:
        g.add_edge(a, b, 1.0)
    bc = g.betweenness()
    check(approx(bc["B"], 2 / 3, 1e-6), "路径图 P4 介数 B=2/3", str(bc))
    check(approx(bc["A"], 0.0), "路径图端点介数=0")
    check(g.diameter() == 3, "直径=3")
    check(g.components() == [["A", "B", "C", "D"]] or len(g.components()) == 1, "连通分量=1")
    check(approx(g.density(), 3 / 6), "密度=0.5")

    g2 = C.Graph()
    for a, b in [("A", "B"), ("B", "C"), ("C", "A")]:
        g2.add_edge(a, b, 1.0)
    check(g2.core_number()["A"] == 2, "三角形 k-core=2")
    check(approx(g2.clustering("A"), 1.0), "三角形聚类系数=1")
    _m2, _i2 = g2.mcc()
    # CytoHubba 定义：Σ(|C|-1)!；三角形只有一个极大团，(3-1)! = 2
    check(approx(_m2["A"], 2.0), "三角形 MCC=2", str(_m2))
    check(_i2["cliques"] >= 1, "MCC 团枚举有结果")

    # 两个 5-clique 弱连接 → Louvain 应分 2 组
    g3 = C.Graph()
    for base in ("X", "Y"):
        for i in range(5):
            for j in range(i + 1, 5):
                g3.add_edge("%s%d" % (base, i), "%s%d" % (base, j), 1.0)
    g3.add_edge("X0", "Y0", 1.0)
    gs = g3.louvain(seed=7)
    check(len(gs) == 2, "Louvain 两社团", "得到 %d 组" % len(gs))
    check(g3.modularity(gs) > 0.3, "模块度 >0.3", "%.3f" % g3.modularity(gs))
    check(sorted(len(x) for x in gs) == [5, 5], "社团大小 [5,5]", str([len(x) for x in gs]))

    # 自环与重复边
    g4 = C.Graph()
    g4.add_edge("A", "A", 1.0)
    g4.add_edge("A", "B", 0.5)
    g4.add_edge("A", "B", 0.9)
    check(g4.self_loops == 1, "自环计数")
    check(g4.dups == 1, "重复边计数")
    check(approx(g4.adj["A"]["B"], 0.9), "重复边取最大权重")

    # 最大连通分量
    g5 = C.Graph()
    g5.add_edge("A", "B", 1.0)
    g5.add_node("Z")
    check(g5.largest_component().n == 2, "最大连通分量大小")

    pos = g3.layout_fr(iters=10, seed=1, width=400, height=300)
    check(len(pos) == 10 and all(0 <= x <= 400 and 0 <= y <= 300 for x, y in pos.values()),
          "布局坐标在画布内")

    pr = g3.pagerank()
    check(approx(sum(pr.values()), 1.0, 1e-6), "PageRank 归一化")


# ---------------------------------------------------------------- IO
def test_io():
    print("[IO] 编码 / 分隔符 / 中文路径 / 列名别名")
    d = tempfile.mkdtemp(prefix="np 自检 ")
    # BOM + CRLF + 中文路径
    p1 = os.path.join(d, "测试 表.csv")
    with open(p1, "w", encoding="utf-8-sig", newline="") as f:
        f.write("基因,score\r\nAKT1,1.2\r\nTNF,3.4\r\n")
    h, rows = C.read_table(p1)
    check(h == ["基因", "score"], "表头去 BOM", str(h))
    check(len(rows) == 2 and rows[0]["基因"] == "AKT1", "CRLF 行解析")
    check(C.resolve_col(h, "gene") == "基因", "中文列名→gene")

    # 分号分隔 + GB18030
    p2 = os.path.join(d, "semi.tsv")
    open(p2, "wb").write("gene;score;source\nIL6;9.1;CTD\n".encode("gb18030"))
    h2, r2 = C.read_table(p2)
    check(h2 == ["gene", "score", "source"], "GB18030 + 分号", str(h2))

    # 别名
    check(C.resolve_col(["Hugo_Symbol"], "gene") == "Hugo_Symbol", "别名 Hugo_Symbol")
    check(C.resolve_col(["靶点名称"], "target_name") == "靶点名称", "别名 靶点名称")
    check(C.resolve_col(["combined_score"], "score") == "combined_score", "别名 combined_score")

    # 写读往返
    p3 = os.path.join(d, "out.csv")
    C.write_table(p3, ["a", "b"], [{"a": 1.0, "b": "中文"}])
    h3, r3 = C.read_table(p3)
    check(r3[0]["b"] == "中文", "写读往返（含中文）")
    check(r3[0]["a"] == "1", "浮点整值不写成 1.0", r3[0]["a"])

    # 基因列表
    p4 = os.path.join(d, "g.txt")
    open(p4, "w", encoding="utf-8").write("AKT1\nTNF; IL6\n\n# 注释\nTP53\n")
    genes = C.read_gene_list(p4)
    check("AKT1" in genes and "IL6" in genes and "TP53" in genes, "基因列表多分隔符", str(genes))

    check(C.split_multi("A;B|C,D") == ["A", "B", "C", "D"], "split_multi")
    check(C.to_float("46.43%") == 46.43, "to_float 去百分号")
    check(C.to_float("NA") is None, "to_float NA")


# ---------------------------------------------------------------- SVG
def test_svg():
    print("[SVG] 图形生成与 XML 合法性")
    d = tempfile.mkdtemp(prefix="np_svg_")
    cases = []
    p = os.path.join(d, "b.svg")
    C.bar_svg([("甲", 3), ("乙", 5)], p, "标题")
    cases.append(("bar", p))
    p2 = os.path.join(d, "v.svg")
    C.venn_svg({"A": {"x", "y"}, "B": {"y", "z"}}, p2, "venn")
    cases.append(("venn", p2))
    g = C.Graph()
    for a, b in [("A", "B"), ("B", "C")]:
        g.add_edge(a, b, 1.0)
    p3 = os.path.join(d, "n.svg")
    C.network_svg(g, p3, "net")
    cases.append(("network", p3))
    p4 = os.path.join(d, "e.svg")
    C.bubble_svg([("通路1", 0.3, 5.0, 8), ("通路2", 0.5, 3.0, 4)], p4, "bubble")
    cases.append(("bubble", p4))
    p5 = os.path.join(d, "l.svg")
    C.line_svg([("s1", [(0, 1), (1, 2), (2, 3)])], p5, "line")
    cases.append(("line", p5))
    for name, path in cases:
        try:
            ET.parse(path)
            ok = True
        except Exception as e:
            ok = False
            print("    ", e)
        check(ok, "SVG 可解析：%s" % name)
        check(os.path.getsize(path) > 200, "SVG 非空：%s" % name)

    # 空数据不应崩
    p6 = os.path.join(d, "empty.svg")
    C.bar_svg([], p6, "空")
    check(os.path.exists(p6), "空数据降级不崩")
    p7 = os.path.join(d, "one.svg")
    C.venn_svg({"A": {"x"}}, p7, "单集合")
    check(os.path.exists(p7), "单集合降级不崩")


# ---------------------------------------------------------------- 溯源
def test_manifest():
    print("[溯源] manifest 合并与哈希")
    d = tempfile.mkdtemp(prefix="np_man_")
    f = os.path.join(d, "in.txt")
    open(f, "w", encoding="utf-8").write("AKT1\n")
    m1 = C.Manifest("01_test", d, {"x": 1})
    m1.input("raw", f)
    m1.count("n", 1)
    m1.warn("测试警告")
    m1.save()
    m2 = C.Manifest("02_test", d, {"y": 2})
    m2.save()
    import json
    data = json.load(open(os.path.join(d, "run_manifest.json"), encoding="utf-8"))
    check(len(data["stages"]) == 2, "manifest 两阶段合并", str(len(data["stages"])))
    check(data["stages"][0]["warnings"] == ["测试警告"], "manifest 保留警告")
    check(len(data["stages"][0]["inputs"]["raw"]["sha256_16"]) == 16, "SHA256 前缀 16 位")


# ---------------------------------------------------------------- CLI
def test_cli():
    print("[CLI] 每个脚本 --help 可执行")
    for fn in sorted(os.listdir(SCRIPTS)):
        if not fn.startswith("np_") or not fn.endswith(".py"):
            continue
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, fn), "--help"],
                           capture_output=True, text=True, timeout=60)
        check(r.returncode == 0, "%s --help" % fn, (r.stderr or "")[:120])
        check("usage" in (r.stdout or "").lower(), "%s 输出 usage" % fn)


# ---------------------------------------------------------------- 端到端
def test_e2e():
    print("[端到端] 演示数据全流程 + 审计")
    out = tempfile.mkdtemp(prefix="np_e2e_")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "np_run.py"), "--demo",
                        "-o", out, "--quiet"], capture_output=True, text=True, timeout=600)
    check(r.returncode == 0, "np_run --demo 退出码 0", (r.stderr or "")[-300:])
    must = ["targets_standardized.csv", "compounds_kept.csv", "intersection.csv",
            "intersection_genes.txt", "hubs.csv", "topology.csv", "null_tests.csv",
            "enrich_significant.csv", "dose_result.json", "docking_graded.csv",
            "reverse_ranking.csv", "np_report.html", "run_manifest.json",
            "venn.svg", "network.svg", "enrich_bubble.svg"]
    for f in must:
        check(os.path.exists(os.path.join(out, f)), "产物存在：%s" % f)
    au = os.path.join(out, "audit", "audit.json")
    check(os.path.exists(au), "审计 JSON 存在")
    if os.path.exists(au):
        import json
        a = json.load(open(au, encoding="utf-8"))
        st = a.get("structural", a)
        check(st.get("fail", 1) == 0, "审计无 FAIL", "fail=%s" % st.get("fail"))
        check(a["score"] >= 85, "审计评分 ≥85", "score=%.1f" % a["score"])
    html = os.path.join(out, "np_report.html")
    if os.path.exists(html):
        t = open(html, encoding="utf-8").read()
        check("<!DOCTYPE html>" in t, "HTML 有 doctype")
        check("<svg" in t, "HTML 内嵌 SVG")
        check("证据等级" in t, "HTML 含证据等级图例")
        check(t.count("<table") >= 5, "HTML 表格数 ≥5", str(t.count("<table")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    print("=" * 60)
    print(" network-pharmacology-pro · 自检")
    print("=" * 60)
    test_stats()
    test_graph()
    test_io()
    test_svg()
    test_manifest()
    test_cli()
    if not a.quick:
        test_e2e()
    print("-" * 60)
    print(" 通过 %d · 失败 %d" % (PASS, FAIL))
    if FAILS:
        print("\n失败明细：")
        for n, d in FAILS:
            print("  · %s %s" % (n, d))
        return 1
    print(" 全部通过 ✅")
    return 0


if __name__ == "__main__":
    sys.exit(main())
