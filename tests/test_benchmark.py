# -*- coding: utf-8 -*-
"""np_benchmark：指标实现必须能被手算验证（AUC/AP/Recall@K/EF/超几何/边级一致性）。"""
import os

import pytest

import np_benchmark as B


def test_roc_auc_perfect_and_worst():
    labels = [1, 1, 0, 0]
    assert B.roc_auc(labels, [4, 3, 2, 1]) == 1.0
    assert B.roc_auc(labels, [1, 2, 3, 4]) == 0.0


def test_roc_auc_all_ties_is_half():
    assert B.roc_auc([1, 1, 0, 0], [1, 1, 1, 1]) == 0.5


def test_average_precision_hand_computed():
    # 排序：1,0,1,0 → precision@1=1, @3=2/3；AP = (1*0.5) + (2/3*0.5)
    ap = B.average_precision([1, 0, 1, 0], [4, 3, 2, 1])
    assert ap == pytest.approx(1 * 0.5 + (2 / 3) * 0.5, abs=1e-4)


def test_recall_and_enrichment():
    ranked = ["A", "B", "C", "D", "E", "F", "G", "H"]
    pos = {"A", "B", "C"}
    m = B.metrics_at_k(ranked, pos, len(ranked), topk=(3, 6))
    assert m["at_k"]["3"]["recall"] == 1.0
    assert m["at_k"]["3"]["precision"] == 1.0
    # 期望命中 = 3 * 3/8 = 1.125 → EF = 3/1.125 = 2.67
    assert m["at_k"]["3"]["enrichment_factor"] == 2.67


def test_hypergeometric_tail():
    # N=8, K=3, n=3, 全中 → p = 1/C(8,3) = 1/56
    assert B.hypergeom_sf(3, 8, 3, 3) == pytest.approx(1 / 56)
    assert B.hypergeom_sf(0, 8, 3, 3) == 1.0


def test_extended_targets_are_excluded_from_negative_pool():
    # strict 口径：extended 不能当负样本（否则 AUC 被系统性压低）
    ranked = [("CORE1", 10), ("EXT1", 9), ("BG1", 8), ("BG2", 7)]
    pos, amb, meta = B.gold_sets({"core": {"CORE1": "L3"}, "extended": {"EXT1": "L1"}}, "strict")
    r = B.evaluate(ranked, pos, amb, topk=(2,), n_boot=0)
    assert r["n_evaluated"] == 3          # EXT1 被剔除
    assert r["n_positive_found"] == 1
    assert meta["tier"].startswith("strict")


def test_rank_from_edges_penalizes_promiscuous_compound(tmp_path):
    p = tmp_path / "edges.csv"
    p.write_text("compound,gene\n"
                 "promiscuous,T1\npromiscuous,T2\npromiscuous,T3\npromiscuous,T4\n"
                 "promiscuous,T5\npromiscuous,T6\npromiscuous,T7\npromiscuous,T8\n"
                 "medium,T9\nmedium,T10\nmedium,T11\n"
                 "selective,T1\nselective,T2\n", encoding="utf-8")
    by_support = dict(B.rank_from_edges(str(p), weight_mode="support"))
    by_pen = dict(B.rank_from_edges(str(p), weight_mode="penalized"))
    # 支持成分数：T1/T2 = 2，其余 = 1
    assert by_support["T1"] == 2
    # 降权后：T1 = 1 + 1 + 0.414（散弹成分贡献被压），仍高于只由散弹成分支持的 T3
    assert by_pen["T1"] > by_pen["T3"]
    assert by_pen["T9"] > by_pen["T3"]      # 3 靶点成分 > 8 靶点成分的单靶贡献


def test_edge_concordance_metrics(tmp_path):
    pred = {("a", "G1"), ("a", "G2"), ("b", "G3")}
    gold = {("a", "G1"), ("a", "G2")}
    r = B.edge_concordance(pred, gold)
    assert r["sensitivity"] == 1.0            # 金标准全覆盖
    assert r["precision"] == pytest.approx(2 / 3, abs=1e-3)
    assert r["over_prediction_factor"] == 1.5


def test_selftest_runs():
    assert B.main(["selftest"]) == 0
