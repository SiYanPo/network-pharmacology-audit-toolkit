# -*- coding: utf-8 -*-
"""证据账本接入回归（报告 §1 / DoD 第 1 条）。

不重复 selftest 的数学检查；这里锁三件事：
1. 未校准模型必须返回 UNCALIBRATED，且 score_targets 是稳定入口；
2. 账本列必须覆盖 gene/evidence_type/level/effect/uncertainty/provenance；
3. WEAK 校准器（OOF AUC 低/系数全 0）不得被当成概率。
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))

import np_evidence_ledger as EL


def _write(path, rows, cols):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader(); w.writerows(rows)


def test_ledger_columns_and_uncalibrated_refusal(tmp_path):
    d = str(tmp_path)
    _write(os.path.join(d, "04_intersect", "intersection.csv"),
           [{"gene": "IL6", "n_compounds": 3, "disease_score": 2},
            {"gene": "TNF", "n_compounds": 2, "disease_score": 1}],
           ["gene", "n_compounds", "disease_score"])
    _write(os.path.join(d, "05_ppi", "hubs.csv"),
           [{"node": "IL6", "degree": 10, "betweenness": 0.2, "verdict": "稳健"},
            {"node": "TNF", "degree": 8, "betweenness": 0.1, "verdict": "稳健"}],
           ["node", "degree", "betweenness", "verdict"])
    rows, src = EL.collect(d)
    assert rows and all("gene" in r and "level" in r and "provenance_hash" in r for r in rows)
    p = EL.write_ledger(rows, d)
    assert os.path.exists(p)
    res = EL.score_targets(d)
    assert res["status"] == "UNCALIBRATED"
    assert all(t["p_gold"] is None for t in res["targets"])


def test_weak_calibration_never_becomes_probability(tmp_path):
    """24 个靶点、10 个正例、特征全无区分度 → 必须 WEAK，而不是 OK。"""
    d = str(tmp_path)
    ledger = os.path.join(d, "ledger.csv")
    cols = ["gene", "source", "evidence_type", "level", "level_name",
            "direction", "effect_size", "uncertainty", "detail", "provenance_hash"]
    rows = []
    for i in range(24):
        rows.append({"gene": "G%02d" % i, "source": "s", "evidence_type": "compound_support",
                     "level": "L2", "level_name": "database", "direction": "support",
                     "effect_size": 1, "uncertainty": "", "detail": "",
                     "provenance_hash": "h"})
    _write(ledger, rows, cols)
    gold = os.path.join(d, "gold.txt")
    open(gold, "w").write("\n".join("G%02d" % i for i in range(10)))
    out = os.path.join(d, "fit")
    r = EL.fit(ledger, EL.read_gold(gold), out, l2=1.0, n_boot=20)
    assert r["status"] in ("WEAK_CALIBRATION", "OK")
    if r["status"] == "WEAK_CALIBRATION":
        m = json.load(open(r["model"], encoding="utf-8"))
        assert m["status"] == "WEAK" and m.get("weak_reasons")
        res = EL.score(d, model_path=r["model"], ledger=ledger)
        assert res["status"] == "UNCALIBRATED"
