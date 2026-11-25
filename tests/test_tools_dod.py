# -*- coding: utf-8 -*-
"""DoD 收尾工具链测试（tools/stp_convert、tools/dod_status）。

重点不是"能跑"，而是"不会把假数据认成已完成"：
- STP 导出的列名有多种写法，解析必须给出明确失败而不是静默产出空表；
- DoD 状态检查必须把 synthetic/demo 之类排除在"外部工具对照"之外。
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for p in (os.path.join(ROOT, "tools"), os.path.join(ROOT, "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import dod_status as DS  # noqa: E402
import stp_convert as SC  # noqa: E402


def _w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w", encoding="utf-8").write(text)
    return path


def test_parse_stp_csv_standard_columns(tmp_path):
    p = _w(str(tmp_path / "STP.csv"),
           "Target,Common name,Uniprot ID,ChEMBL ID,Target Class,Probability\n"
           "PTGS2,COX-2,P35354,CHEMBL230,Enzyme,0.91\n"
           "AKT1,AKT,P31749,CHEMBL321,Kinase,0.42\n")
    pairs, err = SC.parse_one(p)
    assert err == ""
    assert dict(pairs)["PTGS2"] == 0.91
    assert len(pairs) == 2


def test_parse_reports_missing_probability_column(tmp_path):
    p = _w(str(tmp_path / "bad.csv"), "Target,Common name\nPTGS2,COX-2\n")
    pairs, err = SC.parse_one(p)
    assert pairs == []
    assert "no_probability_column" in err


def test_agg_max_vs_mean(tmp_path):
    pairs = [("PTGS2", 0.9), ("PTGS2", 0.5), ("AKT1", 0.4)]
    mx = {r["gene"]: r for r in SC.agg(pairs, "max")}
    mn = {r["gene"]: r for r in SC.agg(pairs, "mean")}
    assert mx["PTGS2"]["score"] == 0.9 and mx["PTGS2"]["n_compounds"] == 2
    assert abs(mn["PTGS2"]["score"] - 0.7) < 1e-9
    # 排序必须按分数降序
    assert SC.agg(pairs, "max")[0]["gene"] == "PTGS2"


def test_semicolon_delimiter_supported(tmp_path):
    p = _w(str(tmp_path / "semi.csv"), "Target;Probability\nTNF;0.8\n")
    pairs, err = SC.parse_one(p)
    assert err == "" and pairs[0][0] == "TNF"


def test_restrict_to_pool(tmp_path):
    pool = _w(str(tmp_path / "pool.txt"), "PTGS2\nVEGFA\n")
    assert SC.read_list(pool) == {"PTGS2", "VEGFA"}


def test_dod_rejects_synthetic_external_tool(tmp_path):
    d = tmp_path / "run"
    d.mkdir()
    (d / "external_comparison.json").write_text(json.dumps(
        {"tools": [{"tool": "STP_synthetic"}, {"tool": "M3_penalized"}]}), encoding="utf-8")
    ok, detail, _ = DS.d4_external(str(d))
    assert not ok, "合成/内部基线不得算作外部工具对照：%s" % detail
    assert "合成" in detail or "内部" in detail


def test_dod_accepts_named_external_tool(tmp_path):
    d = tmp_path / "run"
    d.mkdir()
    (d / "external_comparison.json").write_text(json.dumps(
        {"tools": [{"tool": "SwissTargetPrediction"}, {"tool": "M3_penalized"}]}), encoding="utf-8")
    ok, detail, _ = DS.d4_external(str(d))
    assert ok and "SwissTargetPrediction" in detail


def test_dod_d3_requires_plink_clump_and_f(tmp_path):
    d = tmp_path / "mr" / "IL6"
    d.mkdir(parents=True)
    (d / "mr_result.json").write_text(json.dumps(
        {"status": "OK", "F_statistic": 25, "exploratory_only": False,
         "clump": {"mode": "plink_clump", "status": "OK"}}), encoding="utf-8")
    ok, detail, _ = DS.d3_real_mr(str(tmp_path))
    assert ok and "F=25" in detail
    # 只有距离剪枝 → 不通过
    (d / "mr_result.json").write_text(json.dumps(
        {"status": "OK", "F_statistic": 25, "exploratory_only": True,
         "clump": {"mode": "distance_pruning", "status": "OK"}}), encoding="utf-8")
    ok2, detail2, _ = DS.d3_real_mr(str(tmp_path))
    assert not ok2 and "distance_pruning" in detail2
