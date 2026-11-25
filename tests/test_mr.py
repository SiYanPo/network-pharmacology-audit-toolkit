# -*- coding: utf-8 -*-
"""np_mr：IVW / Steiger / LD 去冗余 / COLOC 的核心断言。

这些是"因果结论"的地基：任何一个静默出错，论文的因果段落就失效。
"""
import pytest

import np_mr as M


def _homogeneous(n=8, true_beta=0.3, se_y=0.01):
    bx = [0.10 + 0.02 * i for i in range(n)]
    by = [true_beta * x for x in bx]
    return bx, [0.02] * n, by, [se_y] * n


def test_ivw_recovers_known_effect():
    bx, sx, by, sy = _homogeneous()
    r = M.ivw(bx, sx, by, sy)
    assert r["beta"] == pytest.approx(0.3, abs=1e-3)
    assert r["heterogeneity"] is False       # 同质数据不应报异质性
    assert r["Q_p"] > 0.05


def test_ivw_flags_heterogeneity():
    bx, sx, by, sy = _homogeneous()
    by = list(by)
    by[0] += 0.5                              # 人为制造离群
    r = M.ivw(bx, sx, by, sy)
    assert r["heterogeneity"] is True
    assert r["Q_p"] < 0.05


def test_harmonise_flips_to_positive_exposure():
    snps = [{"beta_x": -0.2, "se_x": 0.01, "beta_y": 0.4, "se_y": 0.01}]
    out, n = M.harmonise(snps)
    assert n == 1
    assert out[0]["beta_x"] > 0 and out[0]["beta_y"] < 0


def test_steiger_removes_reverse_direction_snp():
    snps = [
        {"snp": "rs_ok", "beta_x": 0.40, "se_x": 0.05, "beta_y": 0.02, "se_y": 0.02, "p_x": 1e-9},
        {"snp": "rs_rev", "beta_x": 0.40, "se_x": 0.05, "beta_y": 0.50, "se_y": 0.005, "p_x": 1e-9},
    ]
    kept, info = M.steiger_filter(snps, n_x=5000, n_y=100000)
    assert [s["snp"] for s in kept] == ["rs_ok"]
    assert info["n_removed"] == 1


def test_steiger_without_sample_size_is_skipped_not_faked():
    snps = [{"snp": "rs1", "beta_x": 0.4, "se_x": 0.05, "beta_y": 0.02, "se_y": 0.02}]
    kept, info = M.steiger_filter(snps, None, None)
    assert info["status"] == "STEIGER_SKIPPED"
    assert kept == snps                       # 不静默丢 SNP，但必须打 flag


def test_distance_pruning_keeps_one_per_window():
    cand = [{"snp": "a", "p_x": 1e-9, "pos": 1_000_000},
            {"snp": "b", "p_x": 1e-8, "pos": 1_050_000},
            {"snp": "c", "p_x": 1e-7, "pos": 5_000_000}]
    kept, info = M.clump(cand, kb=100)
    assert [s["snp"] for s in kept] == ["a", "c"]
    assert info["mode"] == "distance_pruning"
    assert "不是真正的 LD" in info["caveat"]


def test_genotype_ld_clumping_uses_real_r2():
    geno = {"a": [0, 0, 1, 1, 2, 2, 1, 1], "b": [0, 0, 1, 1, 2, 2, 1, 1],
            "c": [0, 2, 0, 2, 0, 2, 0, 2]}
    cand = [{"snp": "a", "p_x": 1e-9}, {"snp": "b", "p_x": 1e-8}, {"snp": "c", "p_x": 1e-7}]
    kept, info = M.clump(cand, geno=geno, r2_thresh=0.8)
    assert [s["snp"] for s in kept] == ["a", "c"]   # b 与 a 完全连锁（r²=1）被剔除
    assert info["mode"] == "genotype_ld"


def test_r2_from_dosage_extremes():
    a = [0, 0, 1, 1, 2, 2]
    assert M.r2_from_dosage(a, a) == pytest.approx(1.0)
    assert M.r2_from_dosage(a, [0, 2, 0, 2, 0, 2]) < 0.2


def test_coloc_prefers_shared_causal_variant():
    shared = [{"beta_x": 0.5, "se_x": 0.05, "beta_y": 0.4, "se_y": 0.05},
              {"beta_x": 0.02, "se_x": 0.05, "beta_y": 0.01, "se_y": 0.05}]
    r = M.coloc_abf(shared)
    assert r["PP.H4"] > 0.5
    assert r["PP.H4"] > r["PP.H3"]


def test_run_mr_end_to_end_flags_are_recorded():
    snps = [{"snp": "rs%d" % i, "beta_x": 0.1 + 0.02 * i, "se_x": 0.02,
             "beta_y": 0.3 * (0.1 + 0.02 * i), "se_y": 0.01, "p_x": 1e-9,
             "pos": 1000 + i * 100_000} for i in range(8)]
    r = M.run_mr(snps, n_x=1000, n_y=50000, kb=10)
    assert r["status"] == "OK"
    assert r["IVW"]["beta"] > 0
    assert r["exposure_effect"] == "risk_increasing"
    assert "NOT_TRUE_LD_CLUMPING" in r["flags"]      # 无 LD 面板必须自曝其短
    assert r["coloc"]["status"] == "OK"
    assert r["leave_one_out"]["status"] == "OK"
