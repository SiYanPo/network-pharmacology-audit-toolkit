# -*- coding: utf-8 -*-
"""Phase 1 模块回归：全局 FDR / 切分 / 快照。"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(HERE, "scripts")
sys.path.insert(0, SCRIPTS)


def _selftest(name):
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, name), "selftest"],
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=180)
    assert p.returncode == 0, p.stdout + p.stderr
    assert "PASS" in p.stdout


def test_fdr_bh_known_values_and_family_global():
    import np_fdr as F
    q = F.bh([0.01, 0.02, 0.03, 0.04, 0.05])
    assert abs(q[0] - 0.05) < 1e-9
    rows, summary = F.apply([
        {"family": "A", "p": 0.01}, {"family": "A", "p": 0.04},
        {"family": "B", "p": 0.20}], alpha=0.05)
    assert summary["n_hypotheses"] == 3 and summary["n_families"] == 2
    assert all("q_global" in r and "q_family" in r for r in rows)
    assert "BH" in F.statement(summary)["zh"]


def test_fdr_selftest():
    _selftest("np_fdr.py")


def test_split_and_snapshot_selftest():
    _selftest("np_split.py")
    _selftest("np_snapshot.py")


def test_snapshot_meta_has_required_version_fields(tmp_path):
    import json
    import np_snapshot as S
    p = tmp_path / "resp.json"
    p.write_text('{"ok": true}', encoding="utf-8")
    meta = S.write_meta(str(p), url="https://example.org/x", release="v1",
                        extra={"etag": "abc"})
    d = json.load(open(meta, encoding="utf-8"))
    for k in ("url", "release", "fetched_at", "sha256", "etag"):
        assert k in d, (k, d)
    assert d["release"] == "v1" and len(d["sha256"]) == 64


def test_lodo_detects_subset_and_reports_both_databases(tmp_path):
    import np_split as S
    p = tmp_path / "edges.tsv"
    p.write_text("compound\tgene\tdatabase\nC1\tIL6\tA\nC1\tTNF\tA\n"
                 "C1\tIL6\tB\nC1\tTNF\tB\n", encoding="utf-8")
    r = S.lodo(str(p), outdir=str(tmp_path / "out"))
    assert r["status"] == "OK" and set(r["databases"]) == {"A", "B"}
    assert r["pairwise_containment"]["A⊂B"] == 1.0
    assert "疑似子集" in r["subset_warning"]


def test_evidence_decision_rules_are_explicit():
    import np_evidence_ledger as EL
    assert EL.DECISION_RULES["wetlab_candidate"]["p_min"] == 0.70
    assert EL.DECISION_RULES["wetlab_candidate"]["min_evidence_types"] == 2
    assert EL.MIN_CV_AUC >= 0.5


def test_coloc_abf_matches_r_coloc_source_formula_and_snapshot():
    """锁住 R coloc 对照：20 位点 PASS，且非饱和集存在。"""
    import json
    import os
    base = os.path.join(HERE, "audit")
    p1 = os.path.join(base, "coloc_cmp_20261122_exact", "coloc_comparison.json")
    p2 = os.path.join(base, "coloc_cmp_soft_20261122", "coloc_comparison.json")
    for p in (p1, p2):
        assert os.path.exists(p), p
        d = json.load(open(p, encoding="utf-8"))
        assert d["status"] == "OK" and d["verdict"] == "PASS", d["verdict"]
        assert d["stats"]["n"] == 20
        assert d["stats"]["pearson_r"] >= 0.99
        assert d["stats"]["mae_pp4"] <= 0.01
