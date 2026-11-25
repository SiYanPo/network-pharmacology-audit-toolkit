# -*- coding: utf-8 -*-
"""评审 §14 外部求证文档的防回归测试。

不访问网络；只锁“查证后写下的边界不得被后续版本悄悄改成强断言”。
"""
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_external_evidence_doc_has_four_audited_claims():
    p = os.path.join(HERE, "references", "36-external-evidence-audit.md")
    s = open(p, encoding="utf-8").read()
    for key in ("PMC10326295", "36370105", "PMC10557418", "2026-11-11"):
        assert key in s, key


def test_limitations_do_not_overclaim_string_or_etcm():
    p = os.path.join(HERE, "references", "33-limitations-and-wording.md")
    s = open(p, encoding="utf-8").read()
    assert "官方承认低度节点补边偏差" in s
    assert "cross-source consistency" in s
    assert "Delphi" in s
