# -*- coding: utf-8 -*-
"""第八轮新增：数据流一致性(G21) / 交互图 / 同源兜底 / 解读引擎 / PDF 与基准脚本。"""
import json
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))


# ---------------- G21 数据流一致性 ----------------
def _mk_flow(tmp_path, genes, edges):
    d = tmp_path / "np_out"
    (d / "04_intersect").mkdir(parents=True)
    (d / "05_ppi").mkdir(parents=True)
    (d / "04_intersect" / "intersection_genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
    lines = ["node1,node2,weight"] + ["%s,%s,0.9" % e for e in edges]
    (d / "05_ppi" / "ppi_edges_clean.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(d)


def test_g21_pass_when_all_intersection_genes_in_ppi(tmp_path):
    import np_finalcheck as F
    d = _mk_flow(tmp_path, ["A", "B", "C"], [("A", "B"), ("B", "C")])
    r = F.check_dataflow(d)
    assert r["coverage"] == 1.0 and r["verdict"] == "PASS"
    assert r["n_missing"] == 0
    assert os.path.exists(os.path.join(d, "dataflow_check.json"))


def test_g21_flags_missing_gene(tmp_path):
    import np_finalcheck as F
    d = _mk_flow(tmp_path, ["A", "B", "ZZZ"], [("A", "B")])
    r = F.check_dataflow(d)
    assert "ZZZ" in r["missing_genes"]
    assert r["verdict"] in ("FAIL", "WARN")
    gate = F.check_gate("G21", "x", "dataflow_check.json", "hard", d)
    assert gate["pass"] is False and "ZZZ" in gate["reason"]


def test_g21_registered_as_hard_gate():
    import np_finalcheck as F
    g = [x for x in F.GATES if x[0] == "G21"]
    assert g and g[0][3] == "hard"


# ---------------- 交互图 ----------------
def test_interactive_build_embeds_data(tmp_path):
    import np_interactive as I
    e = tmp_path / "e.csv"
    e.write_text("node1,node2,weight\nIL6,AKT1,0.9\nAKT1,MTOR,0.8\n", encoding="utf-8")
    h = tmp_path / "h.csv"
    h.write_text("rank,node,degree,verdict\n1,IL6,2,稳健\n2,AKT1,2,稳健\n3,MTOR,1,不稳定\n", encoding="utf-8")
    c = tmp_path / "c.txt"
    c.write_text("IL6\n", encoding="utf-8")
    out = I.build(str(e), str(h), str(c), title="t")
    assert "cytoscape" in out and '"id": "IL6"' in out
    assert '__DATA__' not in out and "节点属性" in out
    assert out.count('"class": "core"') == 1


# ---------------- 同源兜底 ----------------
def test_ortholog_local_covers_textbook_genes():
    import np_ortholog as O
    local = O.load_local()
    for human, mouse in (("BCL2", "Bcl2"), ("VEGFA", "Vegfa"), ("MAPK14", "Mapk14"), ("CXCL8", None)):
        r = O.convert(human, "mouse", local, online=False)
        if mouse:
            assert r["ortholog"] == mouse, human
            assert r["status"] == "resolved_local"
        else:
            # CXCL8 在啮齿类无直系同源：必须给功能类似物 + 说明，而不是裸 NA
            assert r["status"] == "no_direct_ortholog"
            assert "Cxcl1" in r["ortholog"] and r["note"]


def test_ortholog_offline_unresolved_is_explicit():
    import np_ortholog as O
    r = O.convert("NOTAREALGENE", "mouse", {}, online=False)
    assert r["status"] == "unresolved_offline"
    assert r["ortholog"] == ""


# ---------------- 解读引擎 ----------------
def test_interpret_coloc_levels():
    import np_interpret as NI
    assert NI.interpret_coloc(0.9)["level"] == "strong"
    assert NI.interpret_coloc(0.65)["level"] == "moderate"
    assert NI.interpret_coloc(0.2)["level"] == "weak"
    assert NI.interpret_coloc(0.55)["level"] == "weak"     # 0.5–0.6 不再算"中等"
    # 中等档必须包含"排除多效性"的意思（审稿人关注点）
    assert "多效性" in NI.interpret_coloc(0.65)["zh"]


def test_interpret_mr_downgrades_on_weak_instrument():
    import np_interpret as NI
    strong = NI.interpret_mr(0.3, 1e-5, f=40, direction_consistent=True)
    weak = NI.interpret_mr(0.3, 1e-5, f=4, direction_consistent=True)
    assert strong["level"] == "strong"
    assert weak["level"] == "moderate"
    assert "弱工具" in weak["zh"]
    assert strong["en"] and weak["en"]        # 必须给英文话术


def test_interpret_synergy_antagonism_is_not_synergy():
    import np_interpret as NI
    b = NI.interpret_synergy(0.9, "Competitive_Antagonism")
    assert b["level"] == "antagonism"
    assert "竞争性拮抗" in b["zh"]


# ---------------- PDF 与基准脚本 ----------------
def test_pdf_find_browser_returns_tuple():
    import np_pdf
    exe, src = np_pdf.find_browser()
    assert (exe is None) or os.path.exists(exe)
    assert isinstance(src, str) and src


def test_run_benchmark_evaluate_one(tmp_path):
    import run_benchmark as RB
    import np_benchmark as NB
    gold = NB.load_gold()
    f = next(x for x in gold["formulas"] if x["id"] == "CDDP")
    edge = tmp_path / "e.csv"
    lines = ["compound,gene"]
    for g in f["core"]:
        lines += ["c1,%s" % g, "c2,%s" % g]
    for g in ("FAKE1", "FAKE2", "FAKE3"):      # 必须有负类，否则 AUC 无定义
        lines.append("c3,%s" % g)
    edge.write_text("\n".join(lines) + "\n", encoding="utf-8")
    r, formula = RB.evaluate_one("CDDP", str(edge), gold, tier="lenient", boot=0)
    assert r["status"] == "OK"
    assert r["AUC"] > 0.8          # 金标准全命中应接近完美排序
    assert formula["id"] == "CDDP"
