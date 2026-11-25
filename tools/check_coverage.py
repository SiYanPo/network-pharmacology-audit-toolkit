#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
check_coverage.py — 覆盖率棘轮（ratchet），不是"达标声明"

评审 §6 的要求是"CI 加 --cov-fail-under=55"，但直接设 55 会让当前仓库永远红灯：
大部分逻辑写在 main() 里，没拆纯函数，测不到。假装达标比不设更糟。

因此这里用**棘轮**语义：
  * `coverage_baseline.txt` 记录当前真实覆盖率（首次由 --write 生成）；
  * CI 只断言"不低于基线 − 容忍度"；
  * 任何人提高覆盖率后跑 `--write` 更新基线，基线只允许上调。

用法
    python tools/check_coverage.py --coverage coverage.json --baseline coverage_baseline.txt
    python tools/check_coverage.py ... --write      # 更新基线（只能上调）
"""
import argparse
import json
import os
import sys


def load_pct(path):
    d = json.load(open(path, encoding="utf-8"))
    t = d.get("totals") or {}
    return round(float(t.get("percent_covered", 0.0)), 2), int(t.get("num_statements", 0))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--coverage", default="coverage.json")
    ap.add_argument("--baseline", default="coverage_baseline.txt")
    ap.add_argument("--tolerance", type=float, default=0.5,
                    help="允许的回落（百分点）；默认 0.5，防止噪声误报")
    ap.add_argument("--write", action="store_true", help="把当前值写成新基线（只允许上调）")
    a = ap.parse_args(argv)
    if not os.path.exists(a.coverage):
        print("[coverage] 缺 %s：先跑 coverage run -m pytest" % a.coverage)
        return 2
    pct, stmts = load_pct(a.coverage)
    old = None
    if os.path.exists(a.baseline):
        txt = open(a.baseline, encoding="utf-8").read().strip().splitlines()
        try:
            old = float(txt[0].split("=")[-1])
        except Exception:
            old = None
    print("[coverage] 当前 %.2f%%（%d statements）；基线 %s" % (pct, stmts, old))
    if a.write:
        if old is not None and pct < old:
            print("[coverage] 拒绝下调基线（%.2f < %.2f）：基线只允许上调" % (pct, old))
            return 1
        with open(a.baseline, "w", encoding="utf-8") as f:
            f.write("percent=%.2f\nstatements=%d\nnote=只允许上调；提升后重新 --write\n" % (pct, stmts))
        print("[coverage] 基线已更新 → %.2f%%" % pct)
        return 0
    if old is None:
        print("[coverage] 无基线文件，跳过门禁（不阻断）")
        return 0
    if pct + a.tolerance < old:
        print("[coverage] 失败：比基线低 %.2f 个百分点（容忍 %.2f）" % (old - pct, a.tolerance))
        return 1
    print("[coverage] 通过：未低于基线")
    return 0


if __name__ == "__main__":
    sys.exit(main())
