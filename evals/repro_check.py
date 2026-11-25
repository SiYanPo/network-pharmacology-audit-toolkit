#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
repro_check.py — 跨机/跨次复现性校验器

用途：把两份运行产物（同一分析、不同机器或不同时间）逐文件比对，
      区分「科学结果的差异」与「仅时间戳的差异」。

判定
  SAME        逐字节一致
  META_DIFF   仅元数据（时间戳/路径）不同 —— 视为可复现
  DIFF        内容不同 —— **不可复现**，需排查（先看 environment.json）

用法
  python evals/repro_check.py --a out_run1 --b out_run2
  python evals/repro_check.py --a out_run1 --b out_run2 --json diff.json
"""
import argparse
import hashlib
import json
import os
import re
import sys

META_PATTERNS = [
    re.compile(r'"generated":\s*"[^"]*"'),
    re.compile(r'"started":\s*"[^"]*"'),
    re.compile(r'"finished":\s*"[^"]*"'),
    re.compile(r'"cwd":\s*"[^"]*"'),
    re.compile(r'"path":\s*"[^"]*"'),
    re.compile(r'生成时间[^<]*'),
    re.compile(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}'),
    re.compile(r'[A-Za-z]:[\\/]{1,2}[^"\s]+'),      # Windows 绝对路径（含正/反斜杠）
    re.compile(r'/(?:home|Users)/[^\s"]+'),          # POSIX 绝对路径
]


# 这些文件含"产物哈希/时间戳"，天然随运行变化，只作元数据看待
META_ONLY = {"run_manifest.json", "environment.json"}


def digest(path, strip_meta=False):
    raw = open(path, "rb").read()
    if strip_meta:
        try:
            t = raw.decode("utf-8")
        except UnicodeDecodeError:
            return hashlib.sha256(raw).hexdigest()[:16], raw
        for p in META_PATTERNS:
            t = p.sub("<META>", t)
        raw = t.encode("utf-8")
    return hashlib.sha256(raw).hexdigest()[:16], raw


def main(argv=None):
    ap = argparse.ArgumentParser(description="比对两次运行的产物，区分内容差异与时间戳差异")
    ap.add_argument("--a", required=True, help="基准产物目录")
    ap.add_argument("--b", required=True, help="待比对产物目录")
    ap.add_argument("--json", default=None)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    if not (os.path.isdir(a.a) and os.path.isdir(a.b)):
        print("[错误] 目录不存在")
        return 2

    def files(root):
        out = {}
        for dp, _d, fs in os.walk(root):
            for f in fs:
                p = os.path.join(dp, f)
                out[os.path.relpath(p, root)] = p
        return out

    fa, fb = files(a.a), files(a.b)
    only_a = sorted(set(fa) - set(fb))
    only_b = sorted(set(fb) - set(fa))
    rows, n_same, n_meta, n_diff = [], 0, 0, 0

    for name in sorted(set(fa) & set(fb)):
        try:
            h1, r1 = digest(fa[name])
            h2, r2 = digest(fb[name])
        except Exception as e:
            rows.append({"file": name, "status": "ERROR", "detail": str(e)})
            n_diff += 1
            continue
        if h1 == h2:
            n_same += 1
            rows.append({"file": name, "status": "SAME", "hash": h1})
            continue
        if os.path.basename(name) in META_ONLY:
            n_meta += 1
            rows.append({"file": name, "status": "META_DIFF", "hash": h1,
                         "note": "含产物哈希/时间戳，天然变化"})
            continue
        h1m, _ = digest(fa[name], strip_meta=True)
        h2m, _ = digest(fb[name], strip_meta=True)
        if h1m == h2m:
            n_meta += 1
            rows.append({"file": name, "status": "META_DIFF", "hash": h1})
        else:
            n_diff += 1
            rows.append({"file": name, "status": "DIFF", "hash_a": h1, "hash_b": h2})

    print("=" * 62)
    print(" 复现性校验")
    print("=" * 62)
    print(" SAME      %3d   逐字节一致" % n_same)
    print(" META_DIFF %3d   仅时间戳/路径不同（可接受）" % n_meta)
    print(" DIFF      %3d   内容不同（不可复现）" % n_diff)
    if only_a:
        print(" 仅存在于 A：%s" % "、".join(only_a[:8]))
    if only_b:
        print(" 仅存在于 B：%s" % "、".join(only_b[:8]))

    if n_diff:
        print("\n内容不同的文件：")
        for r in rows:
            if r["status"] in ("DIFF", "ERROR"):
                print("  ⛔ %s" % r["file"])
        print("\n排查顺序：")
        print("  1) 比对 environment.json（python / 包版本 / locale）")
        print("  2) 检查是否使用了不同的输入文件（看 run_manifest.json 的 sha256_16）")
        print("  3) 检查随机种子（各阶段 params.seed）")
        print("  4) 若为浮点末位差异，属跨平台预期，报告时应写明平台")
    else:
        print("\n✅ 科学产物完全可复现（差异仅限元数据）")

    if a.json:
        json.dump({"same": n_same, "meta_diff": n_meta, "diff": n_diff,
                   "only_a": only_a, "only_b": only_b, "files": rows},
                  open(a.json, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
        if not a.quiet:
            print("JSON → %s" % a.json)

    return 0 if n_diff == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
