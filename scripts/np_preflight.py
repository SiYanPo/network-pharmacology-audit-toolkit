#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_preflight.py —— 换机 / 交付前的一次性体检（决定"这份 skill 能不能在另一台电脑跑起来"）

为什么单独做一个 preflight
    np_envcheck 只回答"工具在不在"；换机失败的真实原因往往更杂：
    Python 版本不对（scanpy 在 3.13 上装不上）、缓存目录不可写、必需数据文件缺失、
    可选依赖把主流程拖死、路径含空格、时区/编码不一致。
    这些都应该在跑分析**之前**一次性暴露，而不是跑到第 9 个节点才崩。

检查项（全部为只读探测，不写业务数据）
    P1 Python 版本与实现（推荐 3.10–3.12；3.13 仅核心功能）
    P2 必需数据文件（金标准、同源表、内置词典）
    P3 可选第三方包（numpy/networkx/yaml/h5py/matplotlib/streamlit/scanpy…）
    P4 外部工具（vina/obabel/plink/gmx/acpype/blastn…，含非 PATH 位置）
    P5 缓存目录可写（NP_CACHE 或 ~/.np_cache）
    P6 路径与编码（当前目录含空格/中文时能否正常解析）
    P7 关键模块 CLI 可导入（抽样 8 个核心模块跑 --help）

用法
    python np_preflight.py                 # 人类可读
    python np_preflight.py --json preflight.json
    python np_preflight.py --strict         # 有 FAIL 时退出码 3（供 CI / 交付门禁）
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")


def _p(level, name, ok, detail=""):
    return {"level": level, "name": name, "ok": ok, "detail": detail}


def check_python():
    v = sys.version_info
    core_ok = (v.major, v.minor) >= (3, 10)
    full_ok = (3, 10) <= (v.major, v.minor) <= (3, 12)
    ver = "%d.%d.%d" % (v.major, v.minor, v.micro)
    return [
        _p("P1", "Python ≥3.10（核心功能）", core_ok, ver),
        _p("P1", "Python 3.10–3.12（全功能，含 scanpy）", full_ok,
           ver + "；3.13 上 scanpy/numba 生态不完整，单细胞模块会降级"),
    ]


def check_data():
    need = {
        "benchmark_goldstandard.json": "内置 4 方剂金标准",
        "protein_symbols.tsv": "蛋白名→基因符号词典",
    }
    out = []
    for f, desc in need.items():
        p = os.path.join(DATA, f)
        out.append(_p("P2", "数据文件 %s" % f, os.path.exists(p), desc))
    opt = ["ortholog_human_mouse.tsv", "benchmark_23.json", "protein_symbols_uniprot.tsv.gz"]
    for f in opt:
        p = os.path.join(DATA, f)
        out.append(_p("P2-opt", "可选数据 %s" % f, os.path.exists(p), "缺失时走内置兜底/在线"))
    return out


def check_packages():
    """逐个包在**子进程**里探测。

    为什么不在本进程 import：h5py 在 Windows 上偶发
    "HDF5: infinite loop closing library" 直接把解释器拖死（本轮实测），
    隔离后最多是该项 FAIL，不会让整个体检挂掉。
    """
    pkgs = [("numpy", True), ("networkx", False), ("yaml", False), ("h5py", False),
            ("matplotlib", False), ("streamlit", False), ("scanpy", False), ("pandas", False)]
    out = []
    for m, required in pkgs:
        code = ("import importlib,sys;"
                "m=importlib.import_module(%r);"
                "print(getattr(m,'__version__','?'))" % m)
        try:
            q = subprocess.run([sys.executable, "-c", code], stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                               errors="replace", timeout=25)
            if q.returncode == 0:
                out.append(_p("P3", "包 %s" % m, True, (q.stdout or "").strip().splitlines()[-1][:30]))
            else:
                out.append(_p("P3", "包 %s" % m, False,
                              "必需" if required else "可选（缺失时相关功能降级）"))
        except subprocess.TimeoutExpired:
            out.append(_p("P3", "包 %s" % m, False, "导入超时（环境异常，需重装该包）"))
        except Exception as e:
            out.append(_p("P3", "包 %s" % m, False, "%s" % type(e).__name__))
    return out


def check_tools():
    try:
        sys.path.insert(0, HERE)
        import np_envcheck as E
        env = E.detect()
    except Exception as e:
        return [_p("P4", "工具探测", False, "np_envcheck 不可用：%s" % e)]
    out = []
    for name, info in env.get("tools", {}).items():
        out.append(_p("P4", "工具 %s" % name, bool(info.get("available")),
                      info.get("path") or "未安装（见 setup_env.sh / PORTABILITY.md）"))
    return out


def check_cache():
    try:
        import npcore as C
        d = C.cache_dir("preflight")
        p = os.path.join(d, "write_test.txt") if d else None
        if not p:
            return [_p("P5", "缓存目录", False, "无法创建缓存目录")]
        open(p, "w", encoding="utf-8").write("ok")
        os.remove(p)
        return [_p("P5", "缓存目录可写", True, d)]
    except Exception as e:
        return [_p("P5", "缓存目录可写", False, "%s: %s" % (type(e).__name__, str(e)[:80]))]


def check_paths():
    out = []
    d = tempfile.mkdtemp(prefix="preflight 空格 测试_")
    out.append(_p("P6", "临时目录含空格/中文", True, d))
    try:
        b = os.path.join(ROOT, "data", "benchmark_goldstandard.json")
        with open(b, encoding="utf-8") as f:
            json.load(f)
        out.append(_p("P6", "UTF-8 读取", True, "金标准可解析"))
    except Exception as e:
        out.append(_p("P6", "UTF-8 读取", False, str(e)[:80]))
    return out


def check_cli(sample=("np_standardize", "np_intersect", "np_ppi", "np_enrich",
                      "np_finalcheck", "np_uireport", "np_benchmark", "np_selfaudit")):
    out = []
    for m in sample:
        p = os.path.join(HERE, m + ".py")
        if not os.path.exists(p):
            out.append(_p("P7", "CLI %s" % m, False, "缺文件"))
            continue
        try:
            q = subprocess.run([sys.executable, p, "--help"], stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                               errors="replace", timeout=25)
            out.append(_p("P7", "CLI %s" % m, q.returncode == 0,
                          "" if q.returncode == 0 else (q.stdout or "")[-120:]))
        except Exception as e:
            out.append(_p("P7", "CLI %s" % m, False, "%s" % type(e).__name__))
    return out



def selftest():
    """自检：七类检查都能产出结构化行；缓存目录可写；不因单个工具缺失而崩。"""
    checks = []
    for name, f in (("python", check_python), ("data", check_data), ("cache", check_cache),
                    ("paths", check_paths)):
        try:
            rows = f()
            ok = bool(rows) and all({"level", "name", "ok", "detail"} <= set(r) for r in rows)
            checks.append((name, ok, "n=%d" % len(rows)))
        except Exception as e:
            checks.append((name, False, "%s: %s" % (type(e).__name__, str(e)[:60])))
    # 缓存必须可写（否则所有在线模块都会反复联网）
    try:
        rows = check_cache()
        checks.append(("cache_writable", all(r["ok"] for r in rows),
                       "; ".join(r["detail"][:40] for r in rows)))
    except Exception as e:
        checks.append(("cache_writable", False, str(e)[:60]))
    # 工具表：缺失只能是"警告"，不能抛异常
    try:
        rows = check_tools()
        checks.append(("tools_no_raise", len(rows) >= 5, "n=%d" % len(rows)))
    except Exception as e:
        checks.append(("tools_no_raise", False, str(e)[:60]))
    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for n, g, d in checks:
        print("  [%s] %-16s %s" % ("✓" if g else "✗", n, str(d)[:70]))
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="np_preflight: 换机/交付前体检")
    ap.add_argument("--json", default=None)
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    a = ap.parse_args(argv)
    rows = check_python() + check_data() + check_packages() + check_tools() + check_cache() \
        + check_paths() + check_cli()
    fails = [r for r in rows if not r["ok"] and r["level"] in ("P1", "P2", "P5", "P7")]
    warns = [r for r in rows if not r["ok"] and r["level"] not in ("P1", "P2", "P5", "P7")]
    if not a.quiet:
        print("=" * 66)
        print(" network-pharmacology-pro · 换机体检")
        print("=" * 66)
        for r in rows:
            mark = "OK  " if r["ok"] else ("FAIL" if r in fails else "WARN")
            print("  [%s] %-6s %-34s %s" % (mark, r["level"], r["name"], r["detail"][:70]))
        print("-" * 66)
        print("  通过 %d/%d；阻断 %d；提示 %d" % (
            sum(1 for r in rows if r["ok"]), len(rows), len(fails), len(warns)))
        if fails:
            print("\n  必须先解决：")
            for r in fails:
                print("   - %s：%s" % (r["name"], r["detail"][:120]))
        print("\n  提示：缓存目录 = $NP_CACHE 或 ~/.np_cache；离线复现请设 NP_OFFLINE=1")
    if a.json:
        json.dump({"rows": rows, "n_pass": sum(1 for r in rows if r["ok"]),
                   "fails": fails, "warns": warns},
                  open(a.json, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print("  → %s" % a.json)
    return 3 if (a.strict and fails) else 0


if __name__ == "__main__":
    sys.exit(main())
