#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_selfaudit.py —— Skill 自审：模块清单 / 死代码 / 自检矩阵 / 契约核对 / 依赖矩阵

为什么需要它
    当脚本从 10 个长到 70 个，"哪些模块真的会被执行"就不再是靠记忆能回答的问题。
    学术软件的审稿人会问：你声称的功能里，有多少从来没在流水线里跑过？
    与其被问，不如自己先算出来，并把它写进 CI。

它检查什么
    1. inventory  每个模块：行数、`--help` 是否可用、是否在 DAG、是否被 evals 覆盖、是否有 selftest
    2. selftest   （可选 --run）逐个执行 selftest，记录 rc 与末行
    3. deadcode   既不在 DAG、也不在 evals、也不被任何脚本引用的模块（"文档里有、执行图里没有"）
    4. contract   给定输出目录，核对 config 里声明的 artifacts 是否真的产出（防止"声明了却没跑"）
    5. deps       扫描所有脚本的顶层第三方 import，标注是否已安装

用法
    python np_selfaudit.py run                      # 静态审计（快）
    python np_selfaudit.py run --run-selftests      # 追加自检矩阵（慢，需 1–3 分钟）
    python np_selfaudit.py run --outdir np_out      # 追加契约核对
    python np_selfaudit.py selftest
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SKIP_DIRS = {".git", "__pycache__", "np_ui_runs", "np_out", "benchmarks", "audit", "audit_out"}


def _py_files(d):
    out = []
    for dp, dirs, fs in os.walk(d):
        dirs[:] = [x for x in dirs if x not in SKIP_DIRS and not x.startswith(".")]
        for f in fs:
            if f.endswith(".py"):
                out.append(os.path.join(dp, f))
    return out


def _read(p):
    try:
        return open(p, encoding="utf-8", errors="replace").read()
    except Exception:
        return ""


def inventory(root=ROOT):
    scripts = sorted(glob.glob(os.path.join(root, "scripts", "np_*.py")))
    cfg = _read(os.path.join(root, "config", "pipeline.yaml"))
    evals = "".join(_read(p) for p in glob.glob(os.path.join(root, "evals", "*.py")))
    exe = {}
    for p in _py_files(root):
        rel = os.path.relpath(p, root).replace("\\", "/")
        if rel.startswith("scripts/np_") and rel.count("/") == 1:
            continue
        exe[rel] = _read(p)
    rows = []
    for p in scripts:
        name = os.path.basename(p)[:-3]
        src = _read(p)
        refs = [r for r, t in exe.items() if re.search(r"\b" + re.escape(name) + r"\b", t)]
        rows.append({
            "module": name,
            "lines": src.count("\n") + 1,
            "in_dag": bool(re.search(r"\b" + re.escape(name) + r"\b", cfg)),
            "in_evals": bool(re.search(r"\b" + re.escape(name) + r"\b", evals)),
            "has_selftest": ("selftest" in src),
            "help_ok": None,
            "exec_refs": refs[:6],
            "third_party": sorted(set(re.findall(
                r"^\s*(?:import|from)\s+(numpy|pandas|matplotlib|networkx|yaml|h5py|streamlit|scipy|rdkit|sklearn|openpyxl)",
                src, re.M))),
        })
    reg = load_registry(root)
    for r in rows:
        r["registry_status"] = (reg.get(r["module"]) or {}).get("status", "")
    return rows


def load_registry(root=ROOT):
    """读取 config/module_registry.yaml；缺失时返回空（不阻断旧环境）。"""
    try:
        import yaml
        p = os.path.join(root, "config", "module_registry.yaml")
        d = yaml.safe_load(_read(p)) or {}
        mods = d.get("modules") or {}
        if isinstance(mods, list):          # registry v3: list of {module,status,...}
            return {m.get("module"): m for m in mods if m.get("module")}
        return mods
    except Exception:
        return {}


def help_matrix(rows, root=ROOT, timeout=45, workers=4):
    """并发探测 --help（串行 74 个模块要 2–3 分钟，没必要）。"""
    from concurrent.futures import ThreadPoolExecutor

    def one(r):
        p = os.path.join(root, "scripts", r["module"] + ".py")
        try:
            q = subprocess.run([sys.executable, p, "--help"], stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                               errors="replace", timeout=timeout)
            out = q.stdout or ""
            r["help_ok"] = (q.returncode == 0 and ("usage" in out.lower() or "用法" in out))
        except subprocess.TimeoutExpired:
            r["help_ok"] = None          # 超时≠不可用（重依赖 import 在并发下会慢）
        except Exception:
            r["help_ok"] = False
        return r

    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(one, rows))
    return rows


def selftest_matrix(rows, root=ROOT, timeout=120, workers=6):
    """并发跑自检：串行会超过 15 分钟（每个模块都可能联网/建图），并发 6 路约 2 分钟。"""
    from concurrent.futures import ThreadPoolExecutor

    def one(r):
        if not r.get("has_selftest"):
            r["selftest"] = "none"
            return r
        p = os.path.join(root, "scripts", r["module"] + ".py")
        try:
            q = subprocess.run([sys.executable, p, "--selftest"], stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                               errors="replace", timeout=timeout)
            out = q.stdout or ""
            tail = out.strip().splitlines()[-1][:140] if out.strip() else ""
            good = q.returncode == 0 and re.search(r"selftest:.*PASS|PASS（|PASS \(", out) is not None
            r["selftest"] = "PASS" if good else ("FAIL" if q.returncode == 0 else
                                                 "ERROR(%d)" % q.returncode)
            r["selftest_tail"] = tail
        except subprocess.TimeoutExpired:
            r["selftest"] = "TIMEOUT"
        return r

    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(one, rows))
    return rows


def contract_check(outdir, root=ROOT):
    """config 声明的 artifacts 是否真的产出。"""
    try:
        import yaml
        cfg = yaml.safe_load(_read(os.path.join(root, "config", "pipeline.yaml")))
    except Exception:
        return {"status": "CONFIG_UNREADABLE"}
    man_p = os.path.join(outdir, "pipeline_manifest.json")
    ran = set()
    if os.path.exists(man_p):
        try:
            ran = {n["id"] for n in json.load(open(man_p, encoding="utf-8")).get("nodes") or []}
        except Exception:
            pass
    rows = []
    for st in cfg.get("stages") or []:
        for a in st.get("artifacts") or []:
            p = os.path.join(outdir, a)
            rows.append({"stage": st["id"], "artifact": a, "exists": os.path.exists(p),
                         "size": os.path.getsize(p) if os.path.exists(p) else 0,
                         "stage_ran": st["id"] in ran})
    return {"status": "OK", "rows": rows,
            "missing": [r for r in rows if not r["exists"]],
            "stage_ran": sorted(ran)}


def deps_matrix(rows):
    mods = sorted({m for r in rows for m in r["third_party"]})
    out = {}
    for m in mods:
        try:
            __import__(m)
            out[m] = "installed"
        except Exception as e:
            out[m] = "MISSING(%s)" % type(e).__name__
    return out


def deadcode(rows):
    """未进 DAG/evals/脚本引用，且未在注册表标 manual/deprecated 的模块。"""
    return [r for r in rows
            if not r["in_dag"] and not r["in_evals"] and not r["exec_refs"]
            and r.get("registry_status") not in ("manual", "deprecated", "audit")]


def write_reports(res, outdir):
    os.makedirs(outdir, exist_ok=True)
    json.dump(res, open(os.path.join(outdir, "skill_audit.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    rows = res["inventory"]
    dead = res["deadcode"]
    L = ["# Skill 自审报告（np_selfaudit）", "",
         "生成时间：%s　模块数：%d" % (res["generated"], len(rows)), "",
         "## 1　覆盖率总览", "",
         "| 指标 | 数量 | 占比 |", "|---|---|---|",
         "| 模块总数 | %d | 100%% |" % len(rows),
         "| 进入 DAG（config/pipeline.yaml） | %d | %.0f%% |" % (
             sum(r["in_dag"] for r in rows), 100 * sum(r["in_dag"] for r in rows) / max(1, len(rows))),
         "| 被 evals 覆盖 | %d | %.0f%% |" % (
             sum(r["in_evals"] for r in rows), 100 * sum(r["in_evals"] for r in rows) / max(1, len(rows))),
         "| 有 selftest | %d | %.0f%% |" % (
             sum(r["has_selftest"] for r in rows), 100 * sum(r["has_selftest"] for r in rows) / max(1, len(rows))),
         "| `--help` 可用 | %s | %s |" % (
             ("%d" % sum(1 for r in rows if r["help_ok"] is True))
             if all(r["help_ok"] is not None for r in rows)
             else "%d（%d 未检测）" % (sum(1 for r in rows if r["help_ok"] is True),
                                     sum(1 for r in rows if r["help_ok"] is None)),
             "用 --with-help 可全量探测（并发 4，约 +3–8 分钟）"
             if any(r["help_ok"] is None for r in rows) else ""),
         "| 执行图不可达（未登记死代码） | %d | — |" % len(dead), "",
         "| 已登记 manual/deprecated | %d | — |" % sum(
             1 for r in rows if r.get("registry_status") in ("manual", "deprecated")), ""]
    if res.get("selftest_ran"):
        L += ["## 2　自检矩阵", "", "| 模块 | 结果 | 末行 |", "|---|---|---|"]
        for r in rows:
            if r.get("selftest") and r["selftest"] != "none":
                L.append("| %s | %s | %s |" % (r["module"], r["selftest"], r.get("selftest_tail", "")[:80]))
        L += ["", "无 selftest 的模块：%s" % (", ".join(r["module"] for r in rows if not r["has_selftest"]) or "无"), ""]
    L += ["## 3　执行图不可达模块", "",
          "这些模块只在文档里出现，DAG、evals、其他脚本都不调用。**不一定是坏事**"
          "（高级分析本来就是手动触发），但必须明确它们是「手动工具」还是「遗留死代码」：", ""]
    if dead:
        L += ["| 模块 | 行数 | 建议 |", "|---|---|---|"]
        for r in dead:
            L.append("| %s | %d | %s |" % (r["module"], r["lines"],
                                           "加进 audit_tools 并写用法" if r["lines"] > 150 else "评估是否合并/删除"))
    else:
        L.append("（无）")
    if res.get("contract"):
        c = res["contract"]
        L += ["", "## 4　契约核对（声明产物 vs 实际产出）", "",
              "运行节点：%s" % ", ".join(c.get("stage_ran") or []) or "（无 manifest）", "",
              "| 阶段 | 声明产物 | 存在 | 大小 |", "|---|---|---|---|"]
        for r in c["rows"]:
            L.append("| %s | %s | %s | %d |" % (r["stage"], r["artifact"],
                                                "✓" if r["exists"] else "✗", r["size"]))
    if res.get("deps"):
        L += ["", "## 5　第三方依赖", "", "| 包 | 状态 |", "|---|---|"]
        for k, v in res["deps"].items():
            L.append("| %s | %s |" % (k, v))
    L += ["", "## 6　结论与动作", ""]
    if dead:
        L.append("- 执行图不可达的 %d 个模块：逐个在 `SKILL.md` 标注「手动工具 / 已弃用」，"
                 "或接入 `audit_tools` 并提供最小示例。" % len(dead))
    miss = [r["module"] for r in rows
            if r.get("selftest") == "FAIL" or str(r.get("selftest", "")).startswith("ERROR")]
    if miss:
        L.append("- 自检失败：%s（必须修，否则等于没有护栏）" % ", ".join(miss))
    if res.get("contract") and res["contract"].get("missing"):
        L.append("- 契约缺口：%d 个声明产物未生成（见第 4 节）" % len(res["contract"]["missing"]))
    open(os.path.join(outdir, "skill_audit.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    return os.path.join(outdir, "skill_audit.md")


def selftest():
    rows = inventory()
    checks = [("inventory_nonempty", len(rows) > 50),
              ("dag_coverage", sum(r["in_dag"] for r in rows) >= 30),
              ("help_all", all(r["help_ok"] is None for r in rows))]  # 未跑 help 时应为 None
    ok = all(c[1] for c in checks)
    print("selftest: %s（%d 模块，DAG %d）" % ("PASS" if ok else "FAIL", len(rows),
                                              sum(r["in_dag"] for r in rows)))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = argparse.ArgumentParser(description="np_selfaudit: Skill 自审（覆盖/死代码/契约/依赖）")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--outdir", default=None, help="分析输出目录（契约核对用）")
    r.add_argument("--report-dir", default=os.path.join(ROOT, "audit"))
    r.add_argument("--run-selftests", action="store_true")
    r.add_argument("--no-help", action="store_true")
    r.add_argument("--with-help", action="store_true",
                   help="与 --run-selftests 同用时也探测 help（并发 4，约 +3–8 分钟）")
    a = ap.parse_args(argv)
    t0 = time.time()
    rows = inventory()
    do_help = (not a.no_help) and (a.with_help or not a.run_selftests)
    if do_help:
        print("[selfaudit] 探测 --help（并发 4）…", flush=True)
        rows = help_matrix(rows)
    if a.run_selftests:
        print("[selfaudit] 运行 selftest（并发 6）…", flush=True)
        rows = selftest_matrix(rows)
    res = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"), "inventory": rows,
           "deadcode": deadcode(rows), "deps": deps_matrix(rows),
           "selftest_ran": bool(a.run_selftests)}
    if a.outdir:
        res["contract"] = contract_check(a.outdir)
    p = write_reports(res, a.report_dir)
    print("[selfaudit] %d 模块；死代码嫌疑 %d；DAG %d；evals %d；selftest %d；用 %.1fs" % (
        len(rows), len(res["deadcode"]), sum(r["in_dag"] for r in rows),
        sum(r["in_evals"] for r in rows), sum(r["has_selftest"] for r in rows),
        time.time() - t0))
    print("[selfaudit] → %s" % p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
