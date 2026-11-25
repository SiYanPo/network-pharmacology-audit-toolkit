# -*- coding: utf-8 -*-
"""np_decoy：对接阴性对照的 Z-score 判定阈值（决定会不会把假阳性当命中）。"""
import math

import pytest

import np_decoy as D


def test_significant_binding_passes():
    # 核心结合能 -9.0，诱饵平均 -6.0、sd≈0.3 → Z 远小于 -2 → PASS
    decoys = [-6.0, -6.2, -5.8, -6.1, -5.9]
    z, verdict, detail = D.zscore(-9.0, decoys)
    assert verdict == "PASS"
    assert z < -2
    assert detail["n_decoy"] == 5
    assert detail["flag"] == ""


def test_borderline_is_red_flag():
    # 核心与诱饵几乎无差别 → Z≈0 → 必须报 Red Flag，而不是"弱阳性"
    decoys = [-7.0, -7.1, -6.9, -7.05, -6.95]
    z, verdict, detail = D.zscore(-7.0, decoys)
    assert verdict == "RED_FLAG"
    assert D.RED_FLAG in detail["flag"]


def test_too_few_decoys_refuses_to_judge():
    z, verdict, _ = D.zscore(-9.0, [-6.0, -6.1])
    assert z is None
    assert verdict == "NO_CONTROL"


def test_zero_variance_decoys_refuses_to_judge():
    z, verdict, _ = D.zscore(-9.0, [-6.0, -6.0, -6.0, -6.0])
    assert z is None
    assert verdict == "NO_VARIANCE"


def test_zscore_matches_manual_formula():
    decoys = [-5.0, -6.0, -7.0]
    z, _, detail = D.zscore(-8.0, decoys)
    mean = -6.0
    sd = math.sqrt(1.0)  # 样本方差 = ((-1)^2+0+1^2)/2 = 1
    assert z == pytest.approx((-8.0 - mean) / sd, abs=1e-6)
    assert detail["decoy_sd"] == pytest.approx(sd, abs=1e-6)
