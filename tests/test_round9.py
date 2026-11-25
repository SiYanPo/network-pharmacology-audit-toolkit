# -*- coding: utf-8 -*-
"""第九轮：孤立节点扩展 / FuncRep 隔离 / MR 方向校验 / 打印水印 / 自动修复。"""
import json
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))


# ---------------- 孤立节点扩展 ----------------
def _graph(pairs):
    import npcore as C
    g = C.Graph()
    for a, b in pairs:
        g.add_node(a)
        g.add_node(b)
        g.add_edge(a, b, 0.9)
    return g


def test_expand_offline_does_not_fabricate():
    import np_ppi as P
    g = _graph([("A", "B"), ("B", "C")])
    ext, rep = P.expand_isolated(g, ["A", "ZZZ"], g.degree(), min_degree=3, offline=True)
    assert ext is None
    assert rep["status"] == "OFFLINE"
    assert "ZZZ" in rep["seeds"]          # 完全缺失的靶点也要被发现


def test_expand_adds_first_order_neighbors(monkeypatch):
    import np_ppi as P
    g = _graph([("A", "B"), ("B", "C")])
    monkeypatch.setattr(P, "fetch_string_partners",
                        lambda *a, **k: ([("A", "D", 0.9, "experiments",
                                           {"experiments": 0.5, "coexpression": 0.0}),
                                          ("ZZZ", "E", 0.8, "database",
                                           {"experiments": 0.0, "coexpression": 0.3})],
                                         {"status": "OK", "source": "stub"}))
    ext, rep = P.expand_isolated(g, ["A", "ZZZ"], g.degree(), min_degree=3)
    assert rep["status"] == "OK"
    assert "D" in ext.nodes and "E" in ext.nodes
    assert rep["n_added_edges"] == 2
    # 主图不被修改
    assert g.n == 3
    # 分层证据：A 有 1 个实验邻居，ZZZ 有 1 个共表达邻居、0 个实验邻居
    es = rep["evidence_stratified"]
    assert es["A"]["experimental_neighbors"] == 1
    assert es["ZZZ"]["experimental_neighbors"] == 0
    assert es["ZZZ"]["coexpression_neighbors"] == 1


def test_expand_no_action_when_all_degrees_ok(monkeypatch):
    import np_ppi as P
    g = _graph([("A", "B"), ("B", "C"), ("C", "A"), ("A", "D")])
    ext, rep = P.expand_isolated(g, ["A"], g.degree(), min_degree=3)
    assert rep["status"] == "OK_NO_ACTION" and ext is None


# ---------------- FuncRep 隔离 ----------------
def test_funcrep_suffix_and_type():
    import np_ortholog as O
    r = O.convert("CXCL8", "mouse", O.load_local(), online=False)
    assert "(FuncRep)" in r["ortholog"]
    assert r["mapping_type"] == "functional_analog"
    assert r["direct_ortholog"] == "no"


def test_direct_ortholog_is_not_marked():
    import np_ortholog as O
    r = O.convert("TP53", "mouse", O.load_local(), online=False)
    assert "FuncRep" not in r["ortholog"]
    assert r["mapping_type"] == "ortholog" and r["direct_ortholog"] == "yes"


def _write_ortholog(d, rows, warn=True):
    os.makedirs(os.path.join(d, "11_ortholog"), exist_ok=True)
    p = os.path.join(d, "11_ortholog", "mouse_orthologs.csv")
    cols = ["human", "ortholog", "mapping_type", "direct_ortholog", "status"]
    with open(p, "w", encoding="utf-8") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(r.get(c, "")) for c in cols) + "\n")
    if warn:
        open(os.path.join(d, "11_ortholog", "func_rep_warning.md"), "w", encoding="utf-8").write("X")
    return p


def test_g22_passes_when_marked(tmp_path):
    import np_finalcheck as F
    d = str(tmp_path)
    _write_ortholog(d, [{"human": "CXCL8", "ortholog": "Cxcl1(FuncRep)", "mapping_type": "functional_analog",
                         "direct_ortholog": "no", "status": "no_direct_ortholog"}])
    g = F.check_gate("G22", "x", "func_rep_warning.md", "hard", d)
    assert g["pass"] is True


def test_g22_fails_when_funcrep_unmarked(tmp_path):
    import np_finalcheck as F
    d = str(tmp_path)
    _write_ortholog(d, [{"human": "CXCL8", "ortholog": "Cxcl1", "mapping_type": "functional_analog",
                         "direct_ortholog": "no", "status": "no_direct_ortholog"}])
    g = F.check_gate("G22", "x", "func_rep_warning.md", "hard", d)
    assert g["pass"] is False
    assert "(FuncRep)" in g["reason"]


def test_g22_fails_when_warning_file_missing(tmp_path):
    import np_finalcheck as F
    d = str(tmp_path)
    _write_ortholog(d, [{"human": "CXCL8", "ortholog": "Cxcl1(FuncRep)",
                         "mapping_type": "functional_analog", "status": "no_direct_ortholog"}], warn=False)
    g = F.check_gate("G22", "x", "func_rep_warning.md", "hard", d)
    assert g["pass"] is False and "警示文件" in g["reason"]


# ---------------- MR 药物方向 ----------------
def _mr_snps(sign=1):
    return [{"snp": "rs%d" % i, "beta_x": 0.1 + 0.02 * i, "se_x": 0.02,
             "beta_y": sign * 0.3 * (0.1 + 0.02 * i), "se_y": 0.01,
             "p_x": 1e-9, "pos": 1000 + i * 100000} for i in range(8)]


def test_mr_direction_mismatch_flag():
    import np_mr as M
    r = M.run_mr(_mr_snps(-1), n_x=1000, n_y=50000, kb=10, drug_direction="inhibit")
    assert r["IVW"]["beta"] < 0
    assert r["intervention"]["required_action"] == "activate"
    assert r["intervention"]["severity"] == "SEVERE"
    assert "TARGET_INTERVENTION_MISMATCH" in r["flags"]
    assert "升高与疾病风险降低" in r["intervention"]["message"]


def test_mr_direction_consistent():
    import np_mr as M
    r = M.run_mr(_mr_snps(1), n_x=1000, n_y=50000, kb=10, drug_direction="inhibit")
    assert r["intervention"]["consistent"] is True
    assert r["intervention"]["severity"] == "OK"


# ---------------- 报告：打印水印与红灯 ----------------
def test_report_hides_watermark_in_print_and_shows_mr_alert(tmp_path):
    import np_uireport as UR
    d = str(tmp_path)
    os.makedirs(os.path.join(d, "mr", "IL6"))
    mr = {"status": "OK", "IVW": {"beta": -0.3, "p": 1e-6}, "F_statistic": 40,
          "MR_Egger": {"pleiotropy_warning": False}, "coloc": {"status": "OK", "PP.H4": 0.9},
          "direction_consistent": True, "flags": ["TARGET_INTERVENTION_MISMATCH"],
          "intervention": {"drug_direction": "inhibit", "required_action": "activate",
                           "consistent": False, "severity": "SEVERE",
                           "message": "方向矛盾：遗传证据提示该靶点表达升高与疾病风险降低相关"}}
    json.dump(mr, open(os.path.join(d, "mr", "IL6", "mr_result.json"), "w", encoding="utf-8"),
              ensure_ascii=False)
    h = UR.build(d, "test", "DRAFT")
    assert ".wm{display:none!important}" in h          # 打印隐藏平铺水印
    assert "print-only print-note" in h
    assert "靶点-干预方向红灯" in h
    assert "方向矛盾" in h


# ---------------- 自动修复 ----------------
def test_auto_fix_dry_run_plan(tmp_path):
    import np_finalcheck as F
    d = str(tmp_path)
    os.makedirs(os.path.join(d, "05_ppi"))
    open(os.path.join(d, "05_ppi", "hubs.csv"), "w").write("rank,node\n1,IL6\n")
    steps = F.auto_fix(d, dry=True, do_pdf=False)
    names = [s["name"] for s in steps]
    assert any("同源" in n for n in names)
    assert all(s["status"] in ("DRY", "SKIP") for s in steps)
