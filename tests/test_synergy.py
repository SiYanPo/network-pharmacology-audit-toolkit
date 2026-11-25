# -*- coding: utf-8 -*-
"""np_synergy：机制定向协同的判定规则（共享靶点=竞争性拮抗，上下游=协同）。"""
import pytest

nx = pytest.importorskip("networkx")
import np_synergy as S


@pytest.fixture
def layers():
    G = nx.Graph()
    for e in [("EGFR", "PIK3CA"), ("PIK3CA", "AKT1"), ("AKT1", "MTOR"), ("TNF", "RELA")]:
        G.add_edge(*e)
    return S.pathway_layers({"PI3K-Akt": {"EGFR", "PIK3CA", "AKT1", "MTOR"}}, G)


def test_layers_are_ordered_from_receptor(layers):
    lv = layers["PI3K-Akt"]
    assert lv["EGFR"] == 0
    assert lv["AKT1"] > lv["EGFR"]


def test_upstream_downstream_is_mechanistic_synergy(layers):
    r = S.classify_pair({"EGFR"}, {"AKT1"}, layers)
    assert r["label"] == "Mechanistic_Synergy"
    assert r["level_gap"] >= 1


def test_shared_target_forces_competitive_antagonism(layers):
    # 共享 AKT1：无论方向如何，先判"抢同一个口袋"
    r = S.classify_pair({"AKT1"}, {"AKT1", "MTOR"}, layers,
                        actions=({"AKT1": "inhibit"}, {"AKT1": "activate"}))
    assert r["label"] == "Competitive_Antagonism"
    assert "口袋" in r["reason"]


def test_opposite_directions_on_series_nodes_is_antagonism(layers):
    r = S.classify_pair({"EGFR"}, {"AKT1"}, layers,
                        actions=({"EGFR": "inhibit"}, {"AKT1": "activate"}))
    assert r["label"] == "Directional_Antagonism"


def test_same_direction_series_block_is_synergy(layers):
    r = S.classify_pair({"EGFR"}, {"AKT1"}, layers,
                        actions=({"EGFR": "inhibit"}, {"AKT1": "inhibit"}))
    assert r["label"] == "Mechanistic_Synergy"


def test_no_shared_pathway_is_only_topological(layers):
    r = S.classify_pair({"TNF"}, {"MTOR"}, layers)
    assert r["label"] == "Topological_Only"


def test_analyze_mechanistic_separates_labels(layers):
    ct = {"a": {"EGFR"}, "b": {"AKT1"}, "c": {"AKT1", "MTOR"}}
    mech, ant, st = S.analyze_mechanistic(ct, {"PI3K-Akt": {"EGFR", "PIK3CA", "AKT1", "MTOR"}},
                                          None)  # 无 PPI 时用不到 layers，靠 analyse 内部重建
    assert st["mode"] == "mechanistic"
    assert st["n_pairs_evaluated"] == 3
    # 无 PPI 图 → 无法分层 → 只可能出拮抗（共享 AKT1）
    assert st["n_antagonism"] >= 1


def test_backward_compatible_topological_mode():
    ct = {"a": {"TP53", "CASP3", "BAX"}, "b": {"IL6", "TNF", "TP53"}}
    syn, ant, st = S.analyze(ct, pathways={"apoptosis": {"TP53", "CASP3", "BAX", "IL6"}},
                             max_compounds=10, score_cut=0.3)
    assert st["n_pairs_evaluated"] == 1
    assert isinstance(syn, list)
