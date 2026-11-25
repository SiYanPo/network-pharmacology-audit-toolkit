#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""analyze_titles.py — 从 CNKI 语料题名里统计"方法做到哪一步"。

用途：回答"这批中文网络药理学论文的能力边界在哪"——不是靠主观印象，而是靠题名里的
方法词频。注意边界：只统计**题名**，题名未提不代表正文没做（会低估），
所以结果是**下界**，报告中必须这样表述。

输出：method_profile.csv / method_profile.md / method_profile.json
"""
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict

corpus, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

# 方法阶梯（由浅到深），每个词给正则与归类
LAYERS = [
    ("L1 网络药理学基础（成分-靶点-通路）", [r"网络药理"]),
    ("L2 分子对接", [r"分子对接", r"对接"]),
    ("L3 分子动力学/MM-PBSA", [r"分子动力学", r"动力学模拟", r"MD模拟", r"MM-?PBSA", r"MM/?GBSA"]),
    ("L4 湿实验验证", [r"实验验证", r"体内", r"体外", r"细胞实验", r"动物实验", r"大鼠", r"小鼠",
                    r"细胞", r"机制研究", r"作用机制"]),
    ("L5 多组学/网络毒理/其他", [r"转录组", r"测序", r"网络毒理", r"代谢组", r"蛋白组", r"肠道菌群", r"分子对接和"]),
]

rows = list(csv.DictReader(open(corpus, encoding="utf-8-sig")))
n = len(rows)

hit = Counter()
by_year = defaultdict(Counter)
for r in rows:
    t = r.get("题名", "")
    for name, pats in LAYERS:
        if any(re.search(p, t, re.I) for p in pats):
            hit[name] += 1
            y = (re.match(r"(\d{4})", r.get("发表时间", "") or "") or [None, "未知"])[1]
            by_year[y][name] += 1

# 组合：只到对接 vs 做到动力学 vs 做到湿实验
combo = Counter()
for r in rows:
    t = r.get("题名", "")
    d = bool(re.search(r"分子对接|对接", t))
    m = bool(re.search(r"分子动力学|动力学模拟|MM-?PBSA", t, re.I))
    w = bool(re.search(r"实验验证|体内|体外|细胞|大鼠|小鼠|动物", t))
    key = ("对接" if d else "—") + "+" + ("动力学" if m else "—") + "+" + ("实验" if w else "—")
    combo[key] += 1

with open(os.path.join(outdir, "method_profile.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.writer(fh)
    w.writerow(["方法层级", "题名命中篇数", "占比"])
    for name, _ in LAYERS:
        w.writerow([name, hit[name], f"{hit[name]/n:.1%}"])
    w.writerow([])
    w.writerow(["组合（对接+动力学+实验）", "篇数", "占比"])
    for k, v in combo.most_common():
        w.writerow([k, v, f"{v/n:.1%}"])

with open(os.path.join(outdir, "method_profile.md"), "w", encoding="utf-8") as fh:
    fh.write(f"# 题名方法画像（n={n}，仅题名，结果为下界）\n\n| 方法层级 | 命中 | 占比 |\n|---|---|---|\n")
    for name, _ in LAYERS:
        fh.write(f"| {name} | {hit[name]} | {hit[name]/n:.1%} |\n")
    fh.write("\n## 组合分布\n\n| 对接 | 动力学 | 实验 | 篇数 | 占比 |\n|---|---|---|---|---|\n")
    for k, v in combo.most_common():
        d, m, w = k.split("+")
        fh.write(f"| {d} | {m} | {w} | {v} | {v/n:.1%} |\n")
    fh.write("\n## 年度趋势（题名命中）\n\n| 年 | 篇数 | 含对接 | 含动力学 | 含实验 |\n|---|---|---|---|---|\n")
    for y in sorted(by_year):
        tot = sum(1 for r in rows if (re.match(r"(\d{4})", r.get("发表时间", "") or "") or [None, "未知"])[1] == y)
        fh.write(f"| {y} | {tot} | {by_year[y]['L2 分子对接']} | {by_year[y]['L3 分子动力学/MM-PBSA']} | {by_year[y]['L4 湿实验验证']} |\n")

json.dump({"n": n, "layers": {k: hit[k] for k, _ in LAYERS},
           "combo": dict(combo), "by_year": {y: dict(c) for y, c in by_year.items()}},
          open(os.path.join(outdir, "method_profile.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)

print(f"n={n}")
for name, _ in LAYERS:
    print(f"  {name:<28}{hit[name]:>4}  {hit[name]/n:>6.1%}")
print("组合:", combo.most_common(5))
