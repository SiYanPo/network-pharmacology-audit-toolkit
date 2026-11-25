#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""NP 技能工程完整性回归（不替代科学结果回归，但能挡住 80% 升级事故）。

检查：输入契约 → 关键脚本可执行 → 数据/参考齐全 → 策略与契约存在 → 输出摘要 JSON
用法：python evals/integrity.py --input herbs.csv --out out_dir
退出码：0 通过；2 失败（stderr 必含 advice）
"""
import argparse
import csv
import json
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PY = sys.executable
KEY_SCRIPTS = ["np_fetch.py", "np_enrich.py", "np_figs.py", "np_qa.py"]


def fail(msg, advice):
    print(json.dumps({"status": "failed", "message": msg, "advice": advice}, ensure_ascii=False),
          file=sys.stderr)
    return 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    # 1) 输入契约
    if not os.path.exists(a.input):
        return fail(f"输入不存在: {a.input}", "提供 herb 单列 CSV，或用 manual_upload 模板")
    try:
        with open(a.input, encoding="utf-8-sig", newline="") as f:
            rows = list(csv.reader(f, strict=True))
    except csv.Error as e:
        return fail(f"CSV 解析失败: {e}", "单列 CSV、逗号分隔、引号成对")
    if not rows or rows[0][:1] != ["herb"]:
        return fail("表头必须是 herb", "首行写 herb，下面每行一个药名")
    body = [r for r in rows[1:] if r and any(x.strip() for x in r)]
    if not body:
        return fail("输入为 0 行（空用例）", "确认药名列表是否为空；矿物药/动物药可能无记录")
    if any(len(r) != 1 for r in body):
        return fail("存在多列行", "每行仅 1 个药名；多药请分行")

    checks = {}
    # 2) 关键脚本 --help
    broken = []
    for s in KEY_SCRIPTS:
        p = os.path.join(ROOT, "scripts", s)
        if not os.path.exists(p):
            broken.append(s); continue
        r = subprocess.run([PY, p, "--help"], capture_output=True, text=True, timeout=60)
        if r.returncode != 0 or "usage" not in (r.stdout or "").lower():
            broken.append(s)
    checks["scripts_ok"] = not broken
    if broken:
        return fail(f"脚本 --help 失败: {broken}", "修复 argparse 或缺失依赖")

    # 3) 数据与策略
    checks["data_dir"] = os.path.isdir(os.path.join(ROOT, "data"))
    checks["references"] = len([f for f in os.listdir(os.path.join(ROOT, "references"))
                                if f.endswith(".md")]) >= 10
    checks["fallback"] = os.path.exists(os.path.join(ROOT, "policies", "fallback.json"))
    checks["hitl"] = os.path.exists(os.path.join(ROOT, "policies", "hitl.json"))
    checks["contracts"] = os.path.exists(os.path.join(ROOT, "contracts", "standard_tables.md"))
    checks["selftest"] = os.path.exists(os.path.join(ROOT, "evals", "selftest.py"))
    failed = [k for k, v in checks.items() if not v]

    summary = {"status": "ok" if not failed else "failed", "input_rows": len(body),
               "checks": checks, "artifacts": []}
    json.dump(summary, open(os.path.join(a.out, "integrity.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    if failed:
        return fail(f"完整性检查未通过: {failed}", "按 SKILL.md 与 references/ 补齐缺失件")
    print(json.dumps(summary, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
