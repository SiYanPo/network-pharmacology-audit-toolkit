# -*- coding: utf-8 -*-
"""DAG / 产物契约回归（评审 §6 / DoD 第 4 条）。

只锁“声明与执行一致”：
1. stage id 唯一、needs 指向存在节点、artifacts 不得绝对路径；
2. 标准五表表头必须与 contracts/standard_tables.md 同步；
3. np_export check 缺文件必须非零退出（不能打印缺失还 rc=0）。
"""
import csv
import json
import os
import re
import subprocess
import sys

import yaml

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = HERE


def _cfg():
    return yaml.safe_load(open(os.path.join(ROOT, "config", "pipeline.yaml"), encoding="utf-8"))


def test_dag_structure_is_resolvable():
    cfg = _cfg()
    stages = cfg["stages"]
    ids = [s["id"] for s in stages]
    assert len(ids) == len(set(ids)), "stage id 重复"
    for s in stages:
        assert s.get("module") and s.get("cli"), s["id"]
        for need in s.get("needs") or []:
            assert need in ids, "needs 指向不存在 stage: %s -> %s" % (s["id"], need)
        for art in s.get("artifacts") or []:
            assert not os.path.isabs(art), "artifact 不得绝对路径: %s" % art
            assert ".." not in art.split("/"), "artifact 不得越界: %s" % art


def test_standard_table_headers_match_contract_doc():
    md = open(os.path.join(ROOT, "contracts", "standard_tables.md"), encoding="utf-8").read()
    blocks = re.findall(r"##\s+([A-Za-z0-9_]+\.(?:csv|sif))\s*\n`([^`]+)`", md)
    assert blocks, "contracts/standard_tables.md 无可解析表头"
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import np_export
    for fn, cols in blocks:
        want = [c.strip() for c in cols.split(",")]
        assert np_export.STANDARD.get(fn) == want, (fn, np_export.STANDARD.get(fn), want)


def test_output_schema_parses_and_requires_status():
    s = json.load(open(os.path.join(ROOT, "contracts", "output.schema.json"), encoding="utf-8"))
    assert s.get("$schema") and "status" in s.get("required", [])
    assert "ok" in s["properties"]["status"]["enum"]


def test_export_check_fails_on_missing_artifacts(tmp_path):
    # 空目录：五张标准表都缺，必须 rc=1，不能静默通过
    r = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "np_export.py"),
                        "check", "--dir", str(tmp_path)],
                       capture_output=True, text=True, encoding="utf-8", timeout=60)
    assert r.returncode == 1, r.stdout + r.stderr
    assert "[缺失]" in r.stdout


def test_export_check_passes_on_minimal_tables(tmp_path):
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import np_export
    for fn, cols in np_export.STANDARD.items():
        with open(tmp_path / fn, "w", encoding="utf-8-sig", newline="") as f:
            csv.writer(f).writerow(cols)
    r = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "np_export.py"),
                        "check", "--dir", str(tmp_path)],
                       capture_output=True, text=True, encoding="utf-8", timeout=60)
    assert r.returncode == 0, r.stdout + r.stderr


# ==================== v3.0：产物契约校验器 + 模块注册表 ====================
# 上面测的是"DAG 声明自洽"；下面测的是"每个产物文件的列名/类型是否仍是下游期待的样子"。
# 前者防配置写错，后者防实现漂移——两者缺一不可。
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import np_contracts as CT  # noqa: E402


def _w(path, text):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    open(path, "w", encoding="utf-8").write(text)
    return path


def test_contracts_cover_declared_artifacts():
    r = CT.coverage_vs_config(os.path.join(ROOT, "config", "pipeline.yaml"))
    assert r["status"] == "OK", r
    assert r["uncovered"] == [], "未覆盖的产物契约：%s" % r["uncovered"]


def test_every_table_contract_declares_required():
    for name, spec in CT.CONTRACTS.items():
        assert spec["kind"] in ("table", "list", "json", "report"), name
        if spec["kind"] == "table":
            assert spec["required"], "%s：table 契约必须声明 required" % name


def test_contract_missing_column_reported(tmp_path):
    p = _w(str(tmp_path / "hubs.csv"), "rank,node\n1,IL6\n")
    ok, probs = CT.check_file(p, CT.CONTRACTS["05_ppi/hubs.csv"])
    assert not ok and any("degree" in x for x in probs)


def test_contract_alias_accepted(tmp_path):
    p = _w(str(tmp_path / "hubs.csv"), "gene,degree\nIL6,30\nTNF,20\n")
    ok, probs = CT.check_file(p, CT.CONTRACTS["05_ppi/hubs.csv"])
    assert ok, probs


def test_contract_type_violation_reported(tmp_path):
    p = _w(str(tmp_path / "hubs.csv"), "node,degree\nIL6,not_a_number\n")
    ok, probs = CT.check_file(p, CT.CONTRACTS["05_ppi/hubs.csv"])
    assert not ok and any("number" in x for x in probs)


def test_contract_header_only_is_empty(tmp_path):
    p = _w(str(tmp_path / "hubs.csv"), "node,degree\n")
    ok, _ = CT.check_file(p, CT.CONTRACTS["05_ppi/hubs.csv"])
    assert not ok, "只有表头的文件必须判为空表"


def test_contract_optional_missing_is_ok(tmp_path):
    ok, probs = CT.check_file(str(tmp_path / "none.csv"),
                              CT.CONTRACTS["09_evidence/target_scores.csv"])
    assert ok and probs == []


def test_contract_list_and_json_kinds(tmp_path):
    p = _w(str(tmp_path / "g.txt"), "IL6\nTNF\n")
    assert CT.check_file(p, CT.CONTRACTS["04_intersect/intersection_genes.txt"])[0]
    bad = _w(str(tmp_path / "x.json"), "{not json")
    spec = {"required": [], "aliases": {}, "types": {}, "min_rows": 1, "optional": False,
            "kind": "json", "json_keys": [], "note": ""}
    assert not CT.check_file(bad, spec)[0]


def test_every_script_registered_in_registry():
    reg = yaml.safe_load(open(os.path.join(ROOT, "config", "module_registry.yaml"),
                              encoding="utf-8"))
    registered = {m["module"] for m in reg["modules"]}
    scripts = {f[:-3] for f in os.listdir(os.path.join(ROOT, "scripts"))
               if f.startswith("np_") and f.endswith(".py")}
    assert not (scripts - registered), "未登记模块：%s" % sorted(scripts - registered)


def test_manual_modules_documented():
    reg = yaml.safe_load(open(os.path.join(ROOT, "config", "module_registry.yaml"),
                              encoding="utf-8"))
    bad = [m["module"] for m in reg["modules"]
           if m["status"] == "manual" and (not m.get("when") or not m.get("min_cmd"))]
    assert not bad, "manual 模块缺 when/min_cmd：%s" % bad


def test_deprecated_modules_not_used_by_pipeline():
    reg = yaml.safe_load(open(os.path.join(ROOT, "config", "module_registry.yaml"),
                              encoding="utf-8"))
    dep = {m["module"] for m in reg["modules"] if m["status"] == "deprecated"}
    cfg = _cfg()
    used = {s.get("module") for s in (cfg.get("stages") or []) + (cfg.get("audit_tools") or [])}
    assert not (dep & used), "已弃用模块仍被引用：%s" % sorted(dep & used)


def test_registry_status_values_valid():
    reg = yaml.safe_load(open(os.path.join(ROOT, "config", "module_registry.yaml"),
                              encoding="utf-8"))
    ok = {"entry", "core", "audit", "manual", "deprecated"}
    assert not [m["module"] for m in reg["modules"] if m["status"] not in ok]


def test_real_run_contracts_if_available():
    d = os.environ.get("NP_OUT")
    if not d or not os.path.isdir(d):
        import pytest
        pytest.skip("未设置 NP_OUT，跳过真实产物校验")
    r = CT.check_dir(d)
    assert r["verdict"] == "PASS", [x for x in r["results"] if not x["ok"]]
