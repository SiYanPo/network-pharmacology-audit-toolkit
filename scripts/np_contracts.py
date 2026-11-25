#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_contracts.py — 产物契约（schema contract）校验（v3.0）

为什么需要（评审报告 §6）
--------------------------
`config/pipeline.yaml` 声明了 20+ 个 artifacts，但此前**没有任何测试**断言
"这些文件的列名/类型是否还是下游期待的样子"。改一个列名，下游静默出错，
直到出图时才炸——这就是"几十个脚本、随便改个参数把下游弄崩"的根因。

契约的定义方式
--------------
    CONTRACTS = {
        "05_ppi/hubs.csv": {
            "required": ["node", "degree"],        # 必须存在（任一别名命中即可）
            "aliases": {"node": ["node", "gene", "symbol"]},
            "types": {"degree": "number"},
            "min_rows": 1,
        },
        ...
    }

`required` 用**逻辑名**，`aliases` 给出可接受的物理列名——这样既约束了语义，
又不把上游锁死在某一种命名上。

用法
----
    python np_contracts.py check --dir np_out          # 校验一次运行的全部产物
    python np_contracts.py list                        # 打印契约表
    python np_contracts.py selftest
"""

import argparse
import csv
import json
import os
import sys

VERSION = "3.0"


def _c(required, aliases=None, types=None, min_rows=1, optional=False, note="", kind="table",
       json_keys=None):
    """kind: table（带表头）/ list（一行一值）/ json（可解析）/ report（非空文本）"""
    return {"required": required, "aliases": aliases or {}, "types": types or {},
            "min_rows": min_rows, "optional": optional, "note": note, "kind": kind,
            "json_keys": json_keys or []}


CONTRACTS = {
    "12_wetlab/wetlab_result.csv": _c(
        ["experiment", "model", "readout", "value", "n_replicates"],
        {"experiment": ["experiment", "assay", "实验"],
         "model": ["model", "cell_line", "system", "模型"],
         "readout": ["readout", "indicator", "metric", "检测指标"],
         "n_replicates": ["n_replicates", "n", "replicates", "重复数"]},
        {"value": "number", "n_replicates": "int", "p_value": "number"},
        optional=True,
        note=("SCM0061-2021 6.1.3/6.3.3 要求「结果验证」给出模型与检测指标；"
              "此前 skill 只产出实验计划（specificity/localization），没有结果回灌的契约，"
              "导致湿实验证据无法自动进入证据账本。n_replicates 为整数是硬约束："
              "n=1 不得生成 L3 证据")),
    "02_standardize/targets_standardized.csv": _c(
        ["gene"], {"gene": ["gene", "symbol", "target", "target_name", "靶点"]},
        note="标准化后的基因列是全部下游的锚点"),
    "04_intersect/intersection.csv": _c(
        ["gene", "n_compounds"], {"gene": ["gene", "symbol", "target"]},
        {"n_compounds": "int"},
        note="n_compounds 必须是整数：它是第 4 节排序键，也是证据账本的 compound_support"),
    "04_intersect/intersection_genes.txt": _c(
        ["gene"], kind="list", note="纯列表（每行一个基因），不走表头校验"),
    "05_ppi/hubs.csv": _c(
        ["node", "degree"], {"node": ["node", "gene", "symbol"]},
        {"degree": "number", "betweenness": "number"},
        note="node 与 degree 是 G21 与拓扑章节的输入"),
    "05_ppi/ppi_edges_clean.tsv": _c(
        ["node1", "node2", "weight"], {"node1": ["node1", "source", "preferredname_a"],
                                       "node2": ["node2", "target", "preferredname_b"]},
        {"weight": "number"},
        note="扩展子网与交互图都读这张表"),
    "05_ppi/extended_ppi_edges.tsv": _c(
        ["node1", "node2", "weight"], optional=True,
        note="可选：仅 --expand-isolated 时生成；degree 不参与 hub 判定"),
    "07_enrich/enrich_all.csv": _c(
        ["term", "pvalue"], {"term": ["term", "term_name", "name"]},
        {"pvalue": "number", "qvalue": "number"},
        note="pvalue 必须可转 float，否则 FDR 汇总会静默漏掉"),
    "09_evidence/evidence_ledger.csv": _c(
        ["gene", "source", "evidence_type", "level", "provenance_hash"],
        types={"level": "str"}, optional=True,
        note="账本是唯一排序函数的输入；provenance_hash 不可为空"),
    "09_evidence/target_scores.csv": _c(
        ["gene", "decision", "n_evidence_types"], optional=True,
        note="decision 取值必须是 wetlab_candidate/hypothesis/insufficient/uncalibrated"),
    "09_evidence/fdr_table.csv": _c(
        ["family", "unit", "p", "q_family", "q_global"],
        {"p": ["p"]}, {"p": "number", "q_family": "number", "q_global": "number"},
        optional=True, note="q_global 是唯一可写入结论的显著性口径"),
    "00_safety/contraindication.json": _c([], kind="json", note="配伍禁忌结论；键名允许 level/verdict 二选一"),
    "01_alias/herb_alias_result.csv": _c(
        ["raw"], {"raw": ["raw", "herb_raw"]}, note="未匹配的药材必须显式列出（不允许静默丢弃）"),
    "03_compounds/compounds_kept.csv": _c(
        ["compound"], types={"OB": "number", "DL": "number"},
        note="OB/DL 必须为数值：下游阈值筛选依赖它们"),
    "05_ppi/ppi_metrics.csv": _c([], kind="list", optional=True,
                                 note="列随实现变化：只要求存在且非空，不锁列名"),
    "05_ppi/topology.csv": _c(["node"], optional=True,
                              note="PPI 8 指标拓扑宽表；node 为基因列，其余度量列随实现变化"),
    "06_null/null_report.md": _c([], kind="report", optional=True),
    "07_enrich/enrich_report.md": _c([], kind="report", optional=True),
    "08_doseweight/dose_weighted_hubs.csv": _c(
        [], kind="list", optional=True, note="仅 --enable-doseweight 时生成"),
    "09_funnel/final_docking_list.csv": _c(
        ["compound"], optional=True, note="仅 --enable-funnel 时生成"),
    "10_subnetworks/subnetwork_manifest.json": _c([], kind="json", optional=True),
    "09_evidence/decision_rules.json": _c(
        ["decision_rules"], kind="json", optional=True,
        note="决策规则随模型落盘，报告只引用不另写"),
    "09_evidence/multiplicity_statement.json": _c(
        ["n_hypotheses", "n_global"], kind="json", optional=True),
    "versions.json": _c([], kind="json", optional=True, note="数据库版本；缺则 Manuscript 只能写 latest"),
    "np_report.html": _c([], kind="report", optional=True),
    "11_ortholog/mouse_orthologs.csv": _c(
        ["human", "ortholog"], {"human": ["human", "gene"]}, optional=True,
        note="FuncRep 必须带 (FuncRep) 后缀（G22 校验）"),
}


def _read_header(path):
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        hdr = [h.strip().lstrip("\ufeff") for h in first.rstrip("\n").split(delim)]
        total = sum(1 for _ in f)          # seek(0) 后含表头
        n = max(0, total - 1)              # 数据行 = 总行数 − 表头
    return hdr, n


def _is_type(val, kind):
    if val is None or str(val).strip() == "":
        return True                     # 空值不算违约（下游自行处理缺失）
    try:
        if kind == "int":
            int(float(val))
        elif kind == "number":
            float(val)
        return True
    except (TypeError, ValueError):
        return False


def _check_nonempty(path, min_rows=1):
    n = 0
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            if line.strip():
                n += 1
    return (n >= min_rows), ([] if n >= min_rows else ["非空行数 %d < %d" % (n, min_rows)])


def check_file(path, spec):
    """返回 (ok, problems[])。"""
    if not os.path.exists(path):
        return (True, []) if spec["optional"] else (False, ["文件不存在"])
    kind = spec.get("kind") or "table"
    if kind == "list":
        return _check_nonempty(path, spec["min_rows"])
    if kind == "report":
        return _check_nonempty(path, 1)
    if kind == "json":
        try:
            j = json.load(open(path, encoding="utf-8"))
        except Exception as e:
            return False, ["JSON 不可解析：%s" % type(e).__name__]
        probs = [("缺键 %s" % k) for k in spec.get("json_keys") or [] if k not in (j or {})]
        return (not probs), probs
    hdr, n = _read_header(path)
    low = {h.lower(): h for h in hdr}
    probs, col_of = [], {}
    for logical in spec["required"]:
        names = [logical] + spec["aliases"].get(logical, [])
        hit = next((low[x.lower()] for x in names if x.lower() in low), None)
        if not hit:
            uniq = list(dict.fromkeys(names))          # 逻辑名常在别名表里重复，去重后再提示
            probs.append("缺列 %s（可接受：%s）" % (logical, "/".join(uniq)))
        else:
            col_of[logical] = hit
    if n < spec["min_rows"]:
        probs.append("行数 %d < 最小要求 %d" % (n, spec["min_rows"]))
    # 类型抽查（最多 200 行）
    if not probs and spec["types"]:
        with open(path, encoding="utf-8-sig", errors="replace") as f:
            first = f.readline()
            delim = "\t" if first.count("\t") > first.count(",") else ","
            f.seek(0)
            for i, r in enumerate(csv.DictReader(f, delimiter=delim)):
                if i >= 200:
                    break
                for logical, kind in spec["types"].items():
                    col = col_of.get(logical)
                    if col and not _is_type(r.get(col), kind):
                        probs.append("列 %s 存在非 %s 值：%r" % (logical, kind, r.get(col)))
                        break
                if probs:
                    break
    return (not probs), probs


def check_dir(outdir, contracts=None):
    contracts = contracts or CONTRACTS
    res = []
    for rel, spec in contracts.items():
        p = os.path.join(outdir, rel.replace("/", os.sep))
        ok, probs = check_file(p, spec)
        res.append({"artifact": rel, "exists": os.path.exists(p), "ok": ok,
                    "optional": spec["optional"], "problems": probs})
    hard_fail = [r for r in res if not r["ok"] and not r["optional"]]
    return {"outdir": outdir, "n_contracts": len(res),
            "n_present": sum(1 for r in res if r["exists"]),
            "n_fail": len(hard_fail), "verdict": "PASS" if not hard_fail else "FAIL",
            "results": res}


def coverage_vs_config(config_path="config/pipeline.yaml", contracts=None):
    """配置声明的 artifacts 是否都有契约（或显式豁免）。"""
    contracts = contracts or CONTRACTS
    try:
        import yaml
        cfg = yaml.safe_load(open(config_path, encoding="utf-8"))
    except Exception as e:
        return {"status": "NO_CONFIG", "error": str(e)}
    declared = []
    for s in (cfg.get("stages") or []) + (cfg.get("audit_tools") or []):
        for art in (s.get("artifacts") or []):
            declared.append(art.replace("\\", "/"))
    missing = [a for a in declared if a not in contracts]
    return {"status": "OK", "n_declared": len(declared),
            "n_covered": len(declared) - len(missing), "uncovered": missing}


def selftest():
    checks = []
    import tempfile
    d = tempfile.mkdtemp(prefix="np_contract_")
    os.makedirs(os.path.join(d, "05_ppi"), exist_ok=True)
    good = os.path.join(d, "05_ppi", "hubs.csv")
    open(good, "w", encoding="utf-8").write("rank,node,degree,betweenness\n1,IL6,30,0.13\n2,TNF,20,0.08\n")
    ok, probs = check_file(good, CONTRACTS["05_ppi/hubs.csv"])
    checks.append(("good_passes", ok and not probs, str(probs)))

    # 别名命中：用 gene 代替 node
    alt = os.path.join(d, "alt.csv")
    open(alt, "w", encoding="utf-8").write("gene,degree\nIL6,30\n")
    ok2, p2 = check_file(alt, CONTRACTS["05_ppi/hubs.csv"])
    checks.append(("alias_accepted", ok2, str(p2)))

    # 缺列
    bad = os.path.join(d, "bad.csv")
    open(bad, "w", encoding="utf-8").write("rank,node\n1,IL6\n")
    ok3, p3 = check_file(bad, CONTRACTS["05_ppi/hubs.csv"])
    checks.append(("missing_column_detected", (not ok3) and any("degree" in x for x in p3), str(p3)))

    # 类型违约
    badt = os.path.join(d, "badt.csv")
    open(badt, "w", encoding="utf-8").write("node,degree\nIL6,abc\n")
    ok4, p4 = check_file(badt, CONTRACTS["05_ppi/hubs.csv"])
    checks.append(("type_violation_detected", not ok4, str(p4)))

    # 空表（只有表头）
    empty = os.path.join(d, "empty.csv")
    open(empty, "w", encoding="utf-8").write("node,degree\n")
    ok5, p5 = check_file(empty, CONTRACTS["05_ppi/hubs.csv"])
    checks.append(("empty_table_detected", not ok5, str(p5)))

    # 可选产物缺失不算失败
    ok6, _ = check_file(os.path.join(d, "nope.csv"), CONTRACTS["09_evidence/target_scores.csv"])
    checks.append(("optional_missing_ok", ok6, ""))

    # 目录级
    r = check_dir(d)
    checks.append(("dir_verdict_fail", r["verdict"] == "FAIL" and r["n_fail"] >= 1,
                   "%s/%d" % (r["verdict"], r["n_fail"])))

    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for n, g, dd in checks:
        print("  [%s] %-26s %s" % ("✓" if g else "✗", n, dd[:80]))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = argparse.ArgumentParser(description="np_contracts: 产物契约校验")
    sub = ap.add_subparsers(dest="fn")
    c = sub.add_parser("check")
    c.add_argument("--dir", required=True)
    c.add_argument("--json", default=None)
    sp = sub.add_parser("spec", help="用外部契约规范校验单个文件（与内置契约同一套定义语言）")
    sp.add_argument("--file", required=True)
    sp.add_argument("--spec", required=True, help="契约 JSON，字段与内置契约一致")
    l = sub.add_parser("list")
    sub.add_parser("coverage")
    a = ap.parse_args(argv)
    if a.fn == "check":
        r = check_dir(a.dir)
        if a.json:
            json.dump(r, open(a.json, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print("[contracts] %s：%d 个契约，存在 %d，硬失败 %d" % (
            r["verdict"], r["n_contracts"], r["n_present"], r["n_fail"]))
        for x in r["results"]:
            if not x["ok"]:
                print("  [%s] %-46s %s" % ("缺" if not x["exists"] else "违约",
                                           x["artifact"], "; ".join(x["problems"])[:110]))
        return 0 if r["verdict"] == "PASS" else 2
    if a.fn == "spec":
        spec = json.load(open(a.spec, encoding="utf-8"))
        spec = _c(**spec) if "required" in spec else spec
        ok, problems = check_file(a.file, spec)
        name = os.path.basename(a.file)
        print("[contract] %s %s" % ("PASS" if ok else "FAIL", name))
        for x in problems:
            print("  - %s" % x)
        return 0 if ok else 2
    if a.fn == "list":
        for k, v in CONTRACTS.items():
            print("%-46s required=%s optional=%s" % (k, ",".join(v["required"]), v["optional"]))
        return 0
    if a.fn == "coverage":
        r = coverage_vs_config()
        print(json.dumps(r, ensure_ascii=False, indent=1))
        return 0
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
