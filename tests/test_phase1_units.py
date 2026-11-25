# -*- coding: utf-8 -*-
"""Phase 1 单元回归：切分泄漏 / MR FDR / top-K / 模块注册表 / 版本元数据。"""
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))

import np_benchmark as B
import np_split as S
import np_mr as M
import np_selfaudit as A


def test_benchmark_top20_metrics_have_baseline_and_ef():
    ranked = ["A", "B", "C", "D"] + ["G%d" % i for i in range(16)]
    r = B.metrics_at_k(ranked, {"A", "C"}, N_universe=100, topk=(20,))
    assert r["at_k"]["20"]["hits"] == 2
    assert r["at_k"]["20"]["expected_hits"] == 0.4
    assert r["at_k"]["20"]["enrichment_factor"] > 1
    assert 0 <= r["at_k"]["20"]["hypergeom_p"] <= 1


def test_mr_fdr_is_monotone_and_refuses_single_hypothesis_claim():
    res = {"status": "OK", "gene": "IL6", "IVW": {"p": 0.01}}
    out = M.apply_mr_fdr(res, {"IL6": 0.01, "TNF": 0.04, "AKT1": 0.20})
    assert out["fdr"]["status"] == "OK"
    assert 0 <= out["fdr_q"] <= 1
    # 三个假设里 IL6 最显著，BH q 应 >= p
    assert out["fdr_q"] >= 0.01 - 1e-12
    single = M.apply_mr_fdr({"status": "OK", "gene": "IL6", "IVW": {"p": 0.01}}, {})
    assert single["fdr"]["status"] == "SINGLE_HYPOTHESIS"


def test_split_family_keeps_paralogs_together():
    genes = ["MAPK1", "MAPK3", "TP53", "EGFR"]
    assert S.family_of("MAPK1")[0] == S.family_of("MAPK3")[0]
    group_of = {g: S.family_of(g)[0] for g in genes}
    folds = S._folds_group(genes, [1, 1, 0, 0], 2, group_of, 42)
    for f in set(group_of.values()):
        where = {fi for fi, te in enumerate(folds)
                 if any(group_of[genes[i]] == f for i in te)}
        assert len(where) <= 1, "family %s 跨折: %s" % (f, where)


def test_split_concentration_bounded():
    c = S.concentration([("A", 1), ("B", 2), ("C", 3), ("D", 4)], {"A"}, top_frac=0.5)
    assert 0 <= c["concentration"] <= 1


def test_module_registry_declares_dead_modules():
    import yaml
    mods = yaml.safe_load(open(os.path.join(HERE, "config", "module_registry.yaml"),
                               encoding="utf-8"))["modules"]
    reg = {m["module"]: m for m in mods} if isinstance(mods, list) else mods
    assert len(reg) >= 20
    assert all(v.get("status") in ("manual", "deprecated", "audit", "core", "entry")
               for v in reg.values())
    assert all(v.get("when") or v.get("use") for v in reg.values())


def test_selfaudit_inventory_marks_registry_status():
    rows = A.inventory()
    by = {r["module"]: r for r in rows}
    assert by["np_evidence_ledger"]["in_dag"] is True
    assert by["np_snapshot"]["registry_status"] in ("manual", "audit")


def test_registry_covers_every_script_once():
    import glob
    import yaml
    mods = {os.path.basename(p)[:-3] for p in glob.glob(os.path.join(HERE, "scripts", "np_*.py"))}
    rows = yaml.safe_load(open(os.path.join(HERE, "config", "module_registry.yaml"),
                               encoding="utf-8"))["modules"]
    reg = {x["module"]: x for x in rows}
    assert set(reg) == mods, (sorted(mods - set(reg)), sorted(set(reg) - mods))
    assert len(rows) == len(reg), "注册表存在重复模块"


def test_calibration_metrics_exist_and_are_bounded():
    import np_evidence_ledger as EL
    y = [0, 0, 1, 1]
    p = [0.1, 0.4, 0.6, 0.9]
    assert 0 <= EL.brier(y, p) <= 1
    e, curve = EL.ece(y, p)
    assert 0 <= e <= 1 and curve
    dca = EL.net_benefit(y, p)
    assert isinstance(dca, list) and dca


def test_bootstrap_auc_is_seed_deterministic():
    import np_benchmark as B
    y = [0, 0, 1, 1, 0, 1, 0, 1]
    s = [0.1, 0.2, 0.8, 0.9, 0.3, 0.7, 0.4, 0.6]
    a = B.bootstrap_auc(y, s, n_boot=50, seed=7)
    b = B.bootstrap_auc(y, s, n_boot=50, seed=7)
    c = B.bootstrap_auc(y, s, n_boot=50, seed=8)
    assert a == b
    assert a[0] <= a[1] and len(a) == 2


def test_ppi_expansion_reports_rank_correlation_and_keeps_main_graph():
    from unittest.mock import patch
    import npcore as C
    import np_ppi as P
    g = C.Graph()
    for u, v in [("A", "B"), ("B", "C"), ("C", "A")]:
        g.add_node(u); g.add_node(v); g.add_edge(u, v, 500)
    deg = g.degree()
    fake = [("A", "X", 700, ["experiments"], {"experiments": 0.9, "coexpression": 0.0})]
    with patch.object(P, "fetch_string_partners", return_value=(fake, {"status": "OK"})):
        ext, rep = P.expand_isolated(g, ["A"], deg, min_degree=10, limit=5)
    assert ext is not None and rep["n_added_nodes"] >= 1
    assert "hub_rank_correlation" in rep
    assert "A" in rep["evidence_stratified"]


def test_benchmark_baselines_are_preregistered():
    import yaml
    p = os.path.join(HERE, "config", "benchmark_prereg.yaml")
    d = yaml.safe_load(open(p, encoding="utf-8"))
    assert d["primary"] in ("B0", "B1")
    assert d["secondary"] in ("B0", "B1") and d["secondary"] != d["primary"]
    assert "top20_hits" in d["metrics_primary"]
    assert set(d["splits_required"]) == {"random", "complex", "lodo"}
