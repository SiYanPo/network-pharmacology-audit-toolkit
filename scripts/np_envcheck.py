#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_envcheck.py —— 真实环境探测（替代"凭 PATH 猜测 TOOL_MISSING"）

背景（2026-11-22 实测纠错）
    之前多个模块用 shutil.which() 判断工具缺失 → 误判：工具常装在非 PATH 位置。
    v3.0 起所有路径解析统一交给 np_tools（唯一允许知道"旧默认根"的地方）：
      1) 显式参数 → 2) NP_VINA/NP_GMX/... 环境变量 → 3) config/tools.local.json
      → 4) PATH → 5) 旧默认根（$NP_TOOL_ROOT 或 Windows 的 <author_machine>/MDtools）
    本模块只**报告**：每个工具的路径、来源（origin）、版本。

用法
    python np_envcheck.py                 # 打印完整环境报告
    python np_envcheck.py --json out.json # 落盘供其他模块读取
    python np_envcheck.py --selftest
"""
import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# 工具解析统一交给 np_tools（唯一允许知道"旧默认根"的地方）。
# 本模块只负责**报告**：每个工具的路径、来源、版本。
try:
    import np_tools as _T
except Exception:                                   # 允许独立运行（无 sys.path）
    sys.path.insert(0, HERE)
    import np_tools as _T

EXTRA_TOOLS = ["blastn", "isPcr", "bowtie2"]        # 暂无注册表条目，按 PATH 探测


def _resolve(name):
    if name in _T.TOOLS:
        return _T.resolve(name)
    p = shutil.which(name)
    return {"name": name, "status": "ok" if p else "missing", "path": p,
            "origin": "PATH" if p else "not_found", "tried": []}


def detect():
    env = {"os": sys.platform, "python": sys.version.split()[0], "tools": {},
           "python_packages": {}, "wsl": {}, "tool_root": _T.legacy_roots()[:2]}
    for name in list(_T.TOOLS) + EXTRA_TOOLS:
        r = _resolve(name)
        env["tools"][name] = {"path": r["path"], "available": r["status"] == "ok",
                              "origin": r["origin"]}
    # Vina 版本实测
    vp = env["tools"].get("vina", {}).get("path")
    if vp:
        try:
            rr = subprocess.run([vp, "--version"], capture_output=True, text=True, timeout=10)
            env["tools"]["vina"]["version"] = (rr.stdout or rr.stderr).strip()[:80]
        except Exception:
            pass
    return _rest(env)


def _rest(env):
    """版本/包/WSL 等"附加"探测（工具路径已在 detect 里解析完）。"""
    # md 环境 Python 包
    md_py = _T.path_of("python_md") or sys.executable
    for py in dict.fromkeys([sys.executable, md_py]):
        if not os.path.exists(py):
            continue
        try:
            r = subprocess.run([py, "-c",
                "import importlib.util,json;print(json.dumps({m:bool(importlib.util.find_spec(m)) "
                "for m in ['rdkit','meeko','openbabel','scanpy','anndata','vina','pandas']}))"],
                capture_output=True, text=True, timeout=20)
            if r.returncode == 0:
                env["python_packages"][py] = json.loads(r.stdout.strip().splitlines()[-1])
        except Exception:
            pass
    # WSL
    try:
        r = subprocess.run(["wsl", "-l", "-q"], capture_output=True, text=True, timeout=10)
        env["wsl"]["distros"] = [x.strip().strip("\x00") for x in r.stdout.splitlines() if x.strip()]
    except Exception:
        env["wsl"]["distros"] = []
    return env


def install_hint(env):
    miss = [k for k, v in env["tools"].items() if not v["available"]]
    out = ["缺失工具及安装方式："]
    if "gmx" in miss:
        out.append("  GROMACS (WSL): sudo apt install -y gromacs   # 或 conda install -c conda-forge gromacs")
    if "acpype" in miss:
        out.append("  acpype (WSL):  pip install acpype  # 或 conda install -c conda-forge acpype")
    if "blastn" in miss:
        out.append("  BLAST+ (WSL):  sudo apt install -y ncbi-blast+")
    if "isPcr" in miss:
        out.append("  isPcr:          conda install -c bioconda ispcr")
    if "scanpy" in str(env.get("python_packages", {})) and "False" in str(env.get("python_packages", {})):
        out.append("  scanpy:         pip install scanpy anndata")
    out.append("  一键安装（推荐）：bash scripts/setup_env.sh")
    return "\n".join(out)


def selftest():
    e = detect()
    ok = isinstance(e.get("tools"), dict) and "vina" in e["tools"]
    print("selftest: %s (Vina=%s, obabel=%s, gmx=%s, antechamber=%s)" %
          ("PASS" if ok else "FAIL",
           e["tools"].get("vina", {}).get("available"),
           e["tools"].get("obabel", {}).get("available"),
           e["tools"].get("gmx", {}).get("available"),
           e["tools"].get("antechamber", {}).get("available")))
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="np_envcheck: 真实工具环境探测")
    ap.add_argument("--json", default=None)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    e = detect()
    print("=" * 60)
    print(" network-pharmacology-pro 环境报告")
    print("=" * 60)
    for k, v in e["tools"].items():
        mark = "OK " if v["available"] else "MISS"
        print("  [%s] %-12s %-14s %s %s" % (mark, k, v.get("origin", ""),
                                            v.get("path", ""), v.get("version", "")))
    print("\nPython 包：")
    for py, pkgs in e["python_packages"].items():
        print("  %s" % py)
        print("    " + ", ".join("%s=%s" % (k, "Y" if v else "N") for k, v in pkgs.items()))
    print("\nWSL: %s" % ", ".join(e["wsl"]["distros"]))
    print("\n" + install_hint(e))
    if a.json:
        json.dump(e, open(a.json, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print("\n→ %s" % a.json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
