# -*- coding: utf-8 -*-
"""np_alias：中药名四层映射的回归测试（改名/改表最容易静默破坏的地方）。"""
import np_alias as A


def test_exact_alias():
    assert A.match_herb("龟板")[0] == "龟甲"
    assert A.match_herb("龟板")[1] == "L1_exact"


def test_processed_prefix_stripped():
    # "酒龟甲" → 剥离炮制前缀「酒」后命中
    name, layer, _ = A.match_herb("酒龟甲")
    assert name == "龟甲"
    assert layer in ("L2_normalized", "L3_structural")


def test_structural_parentheses():
    # "龟板（醋制）" → 结构拆解出「龟板」→ 龟甲
    name, layer, _ = A.match_herb("龟板（醋制）")
    assert name == "龟甲"
    assert layer == "L3_structural"


def test_unmatched_is_honest():
    name, layer, note = A.match_herb("完全不存在的药材XYZ")
    assert name == ""
    assert layer == "unmatched"
    assert "别名" in note


def test_norm_is_stable_for_case_and_space():
    assert A.norm("  DAN-Shen ") == A.norm("danshen")
