# -*- coding: utf-8 -*-
"""np_intervention：靶点-干预方向一致性（防"反向因果"写出相反临床结论）。"""
import np_intervention as I


def _mr(beta=0.42, p=1e-6, f=55, consistent=True, pleio=False, het=False):
    return {"gene": "IL6", "beta": beta, "p": p, "F": f, "direction_consistent": consistent,
            "pleiotropy": pleio, "heterogeneity": het}


def test_required_action_mapping():
    assert I.required_action(+0.3) == "inhibit"     # 靶点升高→风险升高→需抑制
    assert I.required_action(-0.3) == "activate"


def test_grade_a_requires_full_evidence():
    g, cav = I.grade_mr(_mr())
    assert g == "A" and cav == []
    g2, cav2 = I.grade_mr(_mr(f=5))
    assert g2 == "B" and any("弱工具" in c for c in cav2)


def test_no_causal_evidence_blocks_verdict():
    rows, s = I.check([_mr(p=0.3)], I.load_action_from_dict({"IL6": {"inhibit": ["X"], "activate": [], "unknown": []}}))
    assert rows[0]["verdict"] == "No_Causal_Evidence"
    assert s["n_severe"] == 0


def test_opposite_action_raises_severe_mismatch():
    acts = I.load_action_from_dict({"IL6": {"inhibit": [], "activate": ["X"], "unknown": []}})
    rows, s = I.check([_mr()], acts)
    assert rows[0]["verdict"] == I.SEVERE
    assert rows[0]["severity"] == "SEVERE"
    assert s["verdict"] == "FAIL"
    assert I.discussion_snippet(rows)      # 必须能直接写进 Discussion


def test_consistent_action_passes():
    acts = I.load_action_from_dict({"IL6": {"inhibit": ["Tanshinone IIA"], "activate": [], "unknown": []}})
    rows, s = I.check([_mr()], acts)
    assert rows[0]["verdict"] == "Consistent"
    assert s["verdict"] == "PASS"


def test_weak_evidence_downgrades_to_warn():
    acts = I.load_action_from_dict({"IL6": {"inhibit": [], "activate": ["X"], "unknown": []}})
    rows, s = I.check([_mr(f=5)], acts)
    assert rows[0]["severity"] == "WARN"
    assert s["verdict"] == "WARN"


def test_unknown_action_is_not_judged():
    rows, s = I.check([_mr()], {})
    assert rows[0]["verdict"] == "Unknown_Action"
    assert s["n_severe"] == 0


def test_exit_code_three_on_severe(tmp_path):
    act = tmp_path / "action.csv"
    act.write_text("gene,compound,action\nIL6,X,activate\n", encoding="utf-8")
    mrt = tmp_path / "mr.csv"
    mrt.write_text("gene,beta,p,F\nIL6,0.42,1e-6,55\n", encoding="utf-8")
    rc = I.main(["check", "--mr-table", str(mrt), "--action", str(act), "-o", str(tmp_path / "o")])
    assert rc == 3
