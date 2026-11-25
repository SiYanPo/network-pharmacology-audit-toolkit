#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_benchloop.py —— np_bench.py 的回归测试闭环（解决审计发现的 golden 比对未真正执行）

闭环逻辑（伪代码 → 实现）
    ┌─ capture:  选定 golden 用例 → 在干净目录跑真实流水线 → 收集产物(CSV/JSON/MD)
    │            → 语义规范化(剥离时间戳/绝对路径/临时目录) → 逐文件 sha256 → 写 golden.json
    ├─ verify:   当前代码重跑同样用例 → 同样规范化 → 与 golden 逐文件比对
    │            → 分三类：SAME(完全一致) / META_DIFF(仅元数据差，视为可复现) / DIFF(内容差=回归)
    ├─ report:   输出 pass/fail + 差异定位(首个不同文件与行)，写 benchloop_report.md
    └─ accept:   人工确认改进后，显式 --accept 更新 golden（禁止自动接受，防基准漂移）

    关键：golden 必须来自"已知正确"的一次运行，且纳管进 git；任何代码变更都先 verify。

用法
    python np_benchloop.py capture --case typical --cmd "python np_run.py --demo -o {outdir}" --golden-dir golden/
    python np_benchloop.py verify   --case typical --cmd ... --golden-dir golden/
    python np_benchloop.py accept   --case typical --outdir run_out/ --golden-dir golden/
    python np_benchloop.py selftest
"""
import argparse
import difflib
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")

# 语义规范化：这些字段变化不算回归
NOISE_PATTERNS = [
    (re.compile(r'"(?:started|finished|generated|access_date|timestamp)"\s*:\s*"[^"]*"'), '"DATE": "<normalized>"'),
    (re.compile(r"[A-Za-z]:[\\/][^\s\"']{3,}"), "<ABS_PATH>"),
    (re.compile(r"/tmp/[^\s\"']+|/var/folders/[^\s\"']+"), "<TMP_PATH>"),
    (re.compile(r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}"), "<TIMESTAMP>"),
]


def normalize(text):
    for pat, rep in NOISE_PATTERNS:
        text = pat.sub(rep, text)
    return text


def sha(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def collect(outdir, patterns=("*.csv", "*.tsv", "*.json", "*.md", "*.txt")):
    """收集产物为 {相对路径: 规范化内容 hash}。"""
    files = {}
    for root, _dirs, names in os.walk(outdir):
        for n in names:
            if not any(n.endswith(p[1:]) for p in patterns):
                continue
            p = os.path.join(root, n)
            rel = os.path.relpath(p, outdir).replace("\\", "/")
            if rel.startswith("."):
                continue
            try:
                txt = open(p, encoding="utf-8", errors="ignore").read()
            except Exception:
                continue
            files[rel] = {"sha": sha(normalize(txt)), "nbytes": os.path.getsize(p),
                          "content": normalize(txt)}
    return files


def run_case(cmd, outdir, timeout=600):
    # Windows 路径含反斜杠，shlex(posix=True) 会当转义符吃掉 → 用正斜杠 + posix=False 保留
    safe_dir = outdir.replace("\\", "/")
    full = cmd.replace("{outdir}", safe_dir)
    argv = shlex.split(full, posix=False)
    argv = [x.strip('"') for x in argv]
    t0 = time.time()
    try:
        r = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
        return {"rc": r.returncode, "stdout": r.stdout[-1500:], "stderr": r.stderr[-800:],
                "seconds": round(time.time() - t0, 1)}
    except Exception as e:
        return {"rc": -1, "stderr": "%s: %s" % (type(e).__name__, e), "seconds": 0}


def golden_path(gdir, case):
    return os.path.join(gdir, "golden_%s.json" % case)


def cmd_capture(a):
    """在临时目录跑一次，生成 golden。"""
    gdir = os.path.abspath(a.golden_dir)
    os.makedirs(gdir, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="benchloop_") as td:
        run = run_case(a.cmd, td, a.timeout)
        files = collect(td)
    golden = {"case": a.case, "cmd": a.cmd, "captured": time.strftime("%Y-%m-%d %H:%M:%S"),
              "run": {k: v for k, v in run.items() if k != "stdout"},
              "files": {k: v["sha"] for k, v in files.items()},
              "contents": {k: v["content"] for k, v in files.items()}}
    gp = golden_path(gdir, a.case)
    json.dump(golden, open(gp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("[capture] %s：%d 文件 → %s" % (a.case, len(files), gp))
    return 0 if run["rc"] == 0 else 1


def cmd_verify(a):
    gdir = os.path.abspath(a.golden_dir)
    gp = golden_path(gdir, a.case)
    if not os.path.exists(gp):
        raise SystemExit("[FAIL] 无 golden：%s（先 capture）" % gp)
    golden = json.load(open(gp, encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="benchloop_") as td:
        run = run_case(a.cmd, td, a.timeout)
        files = collect(td)
    gfiles, gcontents = golden.get("files", {}), golden.get("contents", {})
    same, diffs = [], []
    for rel, info in files.items():
        if rel not in gfiles:
            diffs.append({"file": rel, "type": "NEW", "why": "golden 中不存在"})
        elif info["sha"] != gfiles[rel]:
            old = gcontents.get(rel, "").splitlines()
            new = info["content"].splitlines()
            delta = list(difflib.unified_diff(old, new, lineterm="", n=1))
            diffs.append({"file": rel, "type": "DIFF", "why": "内容哈希变化",
                          "first_diff": "\n".join(delta[2:12])})
    for rel in gfiles:
        if rel not in files:
            diffs.append({"file": rel, "type": "MISSING", "why": "产物缺失"})
    # META_DIFF：仅元数据文件差异（manifest/env/reproduce）→ 视为可复现
    META_FILES = ("run_manifest.json", "environment.json", "reproduce.cmd", "reproduce.sh",
                  "reproduce.readme.md", "pipeline_manifest.json")
    meta_only = bool(diffs) and all(
        os.path.basename(d["file"]) in META_FILES for d in diffs)
    verdict = "SAME" if not diffs else ("META_DIFF" if meta_only else "DIFF")
    print("[verify] %s：%s（产物 %d，差异 %d）" % (a.case, verdict, len(files), len(diffs)))
    for d in diffs[:10]:
        print("  [%s] %s — %s" % (d["type"], d["file"], d["why"]))
        if d.get("first_diff"):
            print("    " + d["first_diff"].replace("\n", "\n    ")[:400])
    with open(os.path.join(gdir, "benchloop_report.md"), "a", encoding="utf-8") as f:
        f.write("\n## %s — %s\n\n- 判定：%s\n- 命令：`%s`\n" %
                (time.strftime("%Y-%m-%d %H:%M"), a.case, verdict, a.cmd))
        for d in diffs[:20]:
            f.write("- [%s] `%s`：%s\n" % (d["type"], d["file"], d["why"]))
    return 0 if verdict in ("SAME", "META_DIFF") else 1


def cmd_accept(a):
    """人工确认后更新 golden。"""
    files = collect(os.path.abspath(a.outdir))
    gdir = os.path.abspath(a.golden_dir)
    gp = golden_path(gdir, a.case)
    golden = {"case": a.case, "accepted": time.strftime("%Y-%m-%d %H:%M:%S"),
              "files": {k: v["sha"] for k, v in files.items()},
              "contents": {k: v["content"] for k, v in files.items()}}
    json.dump(golden, open(gp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("[accept] %s golden 已更新：%d 文件 → %s" % (a.case, len(files), gp))
    return 0


def selftest():
    t = '{"started":"2026-01-01T00:00:00","path":"C:/x/y/z.csv","gene":"TP53"}\n'
    n = normalize(t)
    ok = ("<normalized>" in n and "<ABS_PATH>" in n and "TP53" in n)
    print("selftest: %s (normalize=%s)" % ("PASS" if ok else "FAIL", n.strip()))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_benchloop: golden 回归测试闭环")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for c in ("capture", "verify"):
        p = sub.add_parser(c)
        p.add_argument("--case", required=True)
        p.add_argument("--cmd", required=True, help="含 {outdir} 占位符的运行命令")
        p.add_argument("--golden-dir", default=os.path.join(ROOT, "benchmark", "golden"))
        p.add_argument("--timeout", type=int, default=600)
        p.set_defaults(fn=cmd_capture if c == "capture" else cmd_verify)
    p = sub.add_parser("accept")
    p.add_argument("--case", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--golden-dir", default=os.path.join(ROOT, "benchmark", "golden"))
    p.set_defaults(fn=cmd_accept)
    p = sub.add_parser("selftest"); p.set_defaults(fn=lambda a: selftest())
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
