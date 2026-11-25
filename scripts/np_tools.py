#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_tools.py — 外部可执行文件的**唯一解析入口**（v3.0）

为什么需要（评审报告 §7）
--------------------------
此前 6 个模块各自写着 `D:/MDtools/...`。换一台机器（或同事的 Mac）就全部失效，
而且"到底读了哪个二进制"没有统一记录，Methods 里也写不清楚。

本模块把解析顺序固定为四级，并把"用了哪个路径、来自哪一级"写进返回结构：

    1. 显式参数（--vina 之类）
    2. 环境变量  NP_VINA / NP_OBABEL / NP_GMX / NP_PLINK / NP_FPOCKET / NP_GNINA / NP_MD_PYTHON
    3. 机器本地配置 config/tools.local.json  {"vina": "C:/tools/vina.exe"}（不进版本库）
    4. PATH 查找（shutil.which）
    5. 旧默认根目录（**仅此一处**）：$NP_TOOL_ROOT 或 Windows 的 D:/MDtools

约定：`resolve()` 返回 dict，`status` ∈ {ok, missing}，`origin` 说明来源。
任何模块都不得再写死机器路径——`tests/test_tool_registry.py` 会 grep 检查。

用法
----
    python np_tools.py list
    python np_tools.py resolve vina
    python np_tools.py where          # 打印 config/tools.local.json 的路径
    python np_tools.py selftest
"""

import argparse
import json
import os
import shutil
import sys

VERSION = "3.0"

# 逻辑名 → (环境变量, 可执行文件名候选, 旧默认根下的相对路径)
TOOLS = {
    "vina": ("NP_VINA", ["vina", "vina.exe"], ["vina.exe", "vina"]),
    "obabel": ("NP_OBABEL", ["obabel", "obabel.exe"],
               ["envs/md/Library/bin/obabel.exe", "obabel.exe"]),
    "gmx": ("NP_GMX", ["gmx", "gmx.exe", "gmx_mpi"],
            ["envs/gmx/bin/gmx.exe", "envs/gmx/Library/bin/gmx.exe"]),
    "plink": ("NP_PLINK", ["plink", "plink.exe", "plink2"], ["plink/plink.exe"]),
    "fpocket": ("NP_FPOCKET", ["fpocket", "fpocket.exe"], ["fpocket/fpocket.exe"]),
    "gnina": ("NP_GNINA", ["gnina", "gnina.exe"], ["gnina.exe"]),
    "acpype": ("NP_ACPYPE", ["acpype", "acpype.exe"], []),
    "antechamber": ("NP_ANTECHAMBER", ["antechamber"], ["antechamber"]),
    "python_md": ("NP_MD_PYTHON", [], ["envs/md/python.exe", "envs/md/bin/python"]),
    "pdbfixer": ("NP_PDBFIXER", ["pdbfixer"], []),
}

# 只有这些工具的 "解释器 Scripts/" 条目可信（纯 Python console_scripts）
CONSOLE_ONLY = {"acpype", "pdbfixer"}

LOCAL_CFG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                         "config", "tools.local.json")


def legacy_roots():
    """旧默认根：只在这里出现机器相关路径，其它模块必须走 resolve()。"""
    out = []
    env = os.environ.get("NP_TOOL_ROOT")
    if env:
        out.append(env)
    if os.name == "nt":
        out.append("D:/MDtools")
    out += ["/opt", "/usr/local", os.path.expanduser("~/tools")]
    return out


def load_local():
    if not os.path.exists(LOCAL_CFG):
        return {}
    try:
        return json.load(open(LOCAL_CFG, encoding="utf-8")) or {}
    except Exception:
        return {}


def resolve(name, explicit=None):
    """按四级顺序解析。返回 {name, status, path, origin, tried}。"""
    if name not in TOOLS:
        return {"name": name, "status": "unknown", "path": None, "origin": "n/a", "tried": []}
    env_key, exes, rels = TOOLS[name]
    tried = []

    if explicit:
        tried.append(("explicit", explicit))
        if os.path.exists(explicit):
            return {"name": name, "status": "ok", "path": explicit, "origin": "explicit", "tried": tried}

    env_val = os.environ.get(env_key)
    if env_val:
        tried.append((env_key, env_val))
        # 允许直接给可执行文件，或给目录（在目录里找候选名）
        if os.path.isdir(env_val):
            for e in exes:
                p = os.path.join(env_val, e)
                if os.path.exists(p):
                    return {"name": name, "status": "ok", "path": p, "origin": env_key, "tried": tried}
        elif os.path.exists(env_val):
            return {"name": name, "status": "ok", "path": env_val, "origin": env_key, "tried": tried}

    loc = load_local().get(name)
    if loc:
        tried.append(("tools.local.json", loc))
        if os.path.exists(loc):
            return {"name": name, "status": "ok", "path": loc, "origin": "tools.local.json", "tried": tried}

    for e in exes:
        w = shutil.which(e)
        if w:
            tried.append(("PATH", w))
            return {"name": name, "status": "ok", "path": w, "origin": "PATH", "tried": tried}

    # 4.5) 当前解释器自带的 console_scripts（仅限纯 Python 工具：pip 装的 acpype/pdbfixer
    #      落在 Scripts/ 下而常不在 PATH；对 vina/obabel 这类编译型二进制不做此查找，
    #      否则会命中 pip 的薄包装，行为与真实二进制不同）
    for e in (exes if name in CONSOLE_ONLY else []):
        for sub in ("Scripts", "bin"):
            p_ = os.path.join(sys.exec_prefix, sub, e)
            tried.append(("interpreter_scripts", p_))
            if os.path.exists(p_):
                return {"name": name, "status": "ok", "path": p_, "origin": "interpreter_scripts",
                        "tried": tried}

    for root in legacy_roots():
        for rel in rels:
            p = os.path.join(root, rel)
            tried.append(("legacy_root", p))
            if os.path.exists(p):
                return {"name": name, "status": "ok", "path": p, "origin": "legacy_root", "tried": tried}

    return {"name": name, "status": "missing", "path": None, "origin": "not_found", "tried": tried}


def path_of(name, explicit=None):
    r = resolve(name, explicit)
    return r["path"] if r["status"] == "ok" else None


def require(name, explicit=None):
    """需要二进制但缺失时，给出可执行的修复指令（而不是只抛 FileNotFoundError）。"""
    r = resolve(name, explicit)
    if r["status"] != "ok":
        raise RuntimeError(
            "缺少外部工具 %s。三种修法任选其一：\n"
            "  1) 设环境变量 %s=/绝对路径/%s\n"
            "  2) 写 config/tools.local.json：{\"%s\": \"C:/path/to/%s\"}\n"
            "  3) 放进 PATH\n"
            "（已尝试：%s）" % (name, TOOLS[name][0], name, name, name,
                             ", ".join("%s=%s" % (o, p) for o, p in r["tried"]) or "无"))
    return r["path"]


def inventory():
    return {n: resolve(n) for n in TOOLS}


def selftest():
    checks = []
    # 1) 未知名
    r = resolve("no_such_tool")
    checks.append(("unknown_tool", r["status"] == "unknown", r["status"]))
    # 2) 显式路径
    exe = sys.executable
    r2 = resolve("vina", explicit=exe)
    checks.append(("explicit_wins", r2["status"] == "ok" and r2["origin"] == "explicit", r2["origin"]))
    # 3) 环境变量优先级高于 PATH
    old = os.environ.get("NP_VINA")
    try:
        os.environ["NP_VINA"] = exe
        r3 = resolve("vina")
        checks.append(("env_beats_path", r3["origin"] == "NP_VINA", r3["origin"]))
    finally:
        if old is None:
            os.environ.pop("NP_VINA", None)
        else:
            os.environ["NP_VINA"] = old
    # 4) 缺失时给出修复指令而不是裸异常
    os.environ.pop("NP_VINA", None)
    os.environ.pop("NP_GNINA", None)
    try:
        require("gnina")
        checks.append(("require_raises", False, "未抛异常"))
    except RuntimeError as e:
        checks.append(("require_raises", "NP_GNINA" in str(e) and "tools.local.json" in str(e), ""))
    # 5) inventory 覆盖全部工具且结构完整
    inv = inventory()
    checks.append(("inventory_complete", set(inv) == set(TOOLS)
                   and all("status" in v and "origin" in v for v in inv.values()), ""))
    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for n, g, d in checks:
        print("  [%s] %-22s %s" % ("✓" if g else "✗", n, d))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = argparse.ArgumentParser(description="np_tools: 外部工具解析")
    sub = ap.add_subparsers(dest="fn")
    sub.add_parser("list")
    r = sub.add_parser("resolve")
    r.add_argument("name")
    sub.add_parser("where")
    a = ap.parse_args(argv)
    if a.fn == "list":
        inv = inventory()
        for n, v in inv.items():
            print("%-12s %-8s %-16s %s" % (n, v["status"], v["origin"], v["path"] or "-"))
        return 0
    if a.fn == "resolve":
        print(json.dumps(resolve(a.name), ensure_ascii=False, indent=1))
        return 0
    if a.fn == "where":
        print(LOCAL_CFG)
        return 0
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
