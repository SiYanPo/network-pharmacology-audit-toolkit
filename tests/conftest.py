# -*- coding: utf-8 -*-
"""pytest 公共夹具：把 scripts/ 加入 sys.path，并提供临时目录与合成数据。"""
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


@pytest.fixture(scope="session")
def root():
    return ROOT


@pytest.fixture
def tmp_out(tmp_path):
    d = tmp_path / "out"
    d.mkdir()
    return str(d)


def write(path, text):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    return path
