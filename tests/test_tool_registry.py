# -*- coding: utf-8 -*-
"""工具注册表测试（评审 §7）：任何模块都不得再写死机器路径。

这条 lint 规则的价值：换机时"能不能跑"取决于是否只有一处知道 D:/MDtools。
一旦有模块偷偷写回机器路径，这个测试立刻红——而不是等到同事拉代码时才炸。
"""
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")

ALLOW = {"np_tools.py"}          # 唯一允许出现"旧默认根"字面量的文件


def test_no_hardcoded_tool_root_outside_registry():
    offenders = []
    for f in sorted(os.listdir(SCRIPTS)):
        if not f.endswith(".py") or f in ALLOW:
            continue
        txt = open(os.path.join(SCRIPTS, f), encoding="utf-8", errors="replace").read()
        for i, line in enumerate(txt.splitlines(), 1):
            # 只拦"作者机器专有路径"；C:/Program Files 是系统标准位置（浏览器探测用，
            # 且已被 os.environ["ProgramFiles"] 覆盖），不算硬编码。
            if re.search(r"[A-Za-z]:[\\/]MDtools", line):
                offenders.append("%s:%d %s" % (f, i, line.strip()[:70]))
    assert not offenders, "发现硬编码机器路径（应改用 np_tools）：\n" + "\n".join(offenders)


def test_env_vars_documented_in_registry():
    import sys
    sys.path.insert(0, SCRIPTS)
    import np_tools as T
    for name, (env_key, exes, rels) in T.TOOLS.items():
        assert env_key.startswith("NP_"), name
    # 关键工具必须有环境变量（评审 §7 点名的四个）
    for must in ("vina", "obabel", "gmx", "plink"):
        assert must in T.TOOLS


def test_resolve_priority_explicit_then_env(monkeypatch):
    import sys
    sys.path.insert(0, SCRIPTS)
    import np_tools as T
    exe = sys.executable
    r = T.resolve("vina", explicit=exe)
    assert r["status"] == "ok" and r["origin"] == "explicit"
    monkeypatch.setenv("NP_VINA", exe)
    r2 = T.resolve("vina")
    assert r2["origin"] in ("NP_VINA", "explicit")


def test_missing_tool_message_is_actionable():
    import sys
    sys.path.insert(0, SCRIPTS)
    import np_tools as T
    old = os.environ.pop("NP_GNINA", None)
    try:
        with pytest.raises(RuntimeError) as ei:
            T.require("gnina")
        msg = str(ei.value)
        assert "NP_GNINA" in msg and "tools.local.json" in msg and "PATH" in msg
    finally:
        if old:
            os.environ["NP_GNINA"] = old
