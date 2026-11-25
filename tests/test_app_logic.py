# -*- coding: utf-8 -*-
"""app.py 的非 UI 逻辑：执行计划与 inputs.json 组装（UI 崩在第一步最常见）。"""
import json
import os
import sys

import pytest


def _import_app():
    pytest.importorskip("streamlit")
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import app
    return app


def test_plan_nodes_follows_dag_order():
    app = _import_app()
    base = app.plan_nodes([])
    assert base[0] == "safety"                  # 安全前置必须第一
    assert "alias" in base and "report" in base
    assert base.index("alias") < base.index("intersect") < base.index("ppi")
    # 插件默认不进计划，开启后追加且仍满足依赖顺序
    with_funnel = app.plan_nodes(["enable-funnel"])
    assert "funnel" in with_funnel
    assert with_funnel.index("compounds") < with_funnel.index("funnel")


def test_build_inputs_writes_resolvable_paths(tmp_path):
    app = _import_app()
    out = str(tmp_path / "run")
    p = app.build_inputs(out, "丹参,三七,冰片", "冠心病")
    cfg = json.load(open(p, encoding="utf-8"))
    for k in ("compounds", "compound_edges", "disease", "edges"):
        assert os.path.exists(cfg["inputs"][k]), k
    assert cfg["inputs"]["herbs"] == "丹参,三七,冰片"
    herbs = open(cfg["inputs"]["herbs_file"], encoding="utf-8").read().split()
    assert herbs == ["丹参", "三七", "冰片"]
    assert cfg["params"]["label"].startswith("丹参,三七,冰片")
