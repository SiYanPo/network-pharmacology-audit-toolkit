#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""merge_cnki.py — 合并逐页 CNKI CSV，去重，出语料台账与计量摘要。

输入：目录下 p1.csv … pN.csv（列含 题名/作者/来源/发表时间/数据库/被引/下载/操作）
输出：cnki_corpus.csv（去重后全量）、cnki_summary.md（年度/期刊/被引/作者）、cnki_stats.json

去重键：题名 + 来源 + 发表时间（CNKI 同一篇在不同排序下可能重复出现）。
"""
import csv
import glob
import json
import os
import re
import sys
from collections import Counter

src, outdir, label = sys.argv[1], sys.argv[2], (sys.argv[3] if len(sys.argv) > 3 else "CNKI corpus")
os.makedirs(outdir, exist_ok=True)

FIELDS = ["题名", "作者", "来源", "发表时间", "数据库", "被引", "下载"]


def pick(row, hdr, name):
    for i, h in enumerate(hdr):
        if h.strip() == name:
            return row[i] if i < len(row) else ""
    return ""


recs, per_page = [], {}
for p in sorted(glob.glob(os.path.join(src, "p*.csv")),
                key=lambda x: int(re.search(r"p(\d+)\.csv$", x).group(1))):
    page = int(re.search(r"p(\d+)\.csv$", p).group(1))
    rows = list(csv.reader(open(p, encoding="utf-8-sig")))
    if len(rows) < 2:
        continue
    hdr = rows[0]
    n = 0
    for r in rows[1:]:
        t = pick(r, hdr, "题名").strip()
        if not t:
            continue
        rec = {k: pick(r, hdr, k).strip() for k in FIELDS}
        rec["题名"] = re.sub(r"(网络首发|增强出版|免费)", "", rec["题名"]).strip()
        rec["page"] = page
        recs.append(rec)
        n += 1
    per_page[page] = n

seen, uniq = set(), []
for r in recs:
    k = (r["题名"], r["来源"], r["发表时间"])
    if k in seen:
        continue
    seen.add(k)
    uniq.append(r)

with open(os.path.join(outdir, "cnki_corpus.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS + ["page"])
    w.writeheader()
    w.writerows(uniq)


def year(s):
    m = re.match(r"(\d{4})", s or "")
    return m.group(1) if m else "未知"


def num(s):
    try:
        return int(float(s))
    except Exception:
        return 0


years = Counter(year(r["发表时间"]) for r in uniq)
journals = Counter(r["来源"] for r in uniq if r["来源"])
authors = Counter()
for r in uniq:
    for a in re.split(r"[;；]", r["作者"]):
        a = a.strip()
        if a:
            authors[a] += 1
cited = sorted(uniq, key=lambda r: -num(r["被引"]))[:15]
downs = sorted(uniq, key=lambda r: -num(r["下载"]))[:15]

stats = {
    "label": label,
    "pages": per_page,
    "raw_rows": len(recs),
    "unique_rows": len(uniq),
    "dup_removed": len(recs) - len(uniq),
    "year_distribution": dict(sorted(years.items())),
    "n_journals": len(journals),
    "top_journals": journals.most_common(12),
    "top_authors": authors.most_common(12),
    "cited_sum": sum(num(r["被引"]) for r in uniq),
    "download_sum": sum(num(r["下载"]) for r in uniq),
    "top_cited": [{"题名": r["题名"], "来源": r["来源"], "年": year(r["发表时间"]),
                   "被引": num(r["被引"])} for r in cited],
    "top_downloaded": [{"题名": r["题名"], "来源": r["来源"], "年": year(r["发表时间"]),
                        "下载": num(r["下载"])} for r in downs],
}
json.dump(stats, open(os.path.join(outdir, "cnki_stats.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)

with open(os.path.join(outdir, "cnki_summary.md"), "w", encoding="utf-8") as fh:
    fh.write(f"# {label}\n\n")
    fh.write(f"- 抽取页数：{len(per_page)}（每页 20 条，末页 {list(per_page.values())[-1]} 条）\n")
    fh.write(f"- 原始行 {len(recs)}，去重后 **{len(uniq)}**（移除 {len(recs)-len(uniq)}）\n")
    fh.write(f"- 期刊数 {len(journals)}；总被引 {stats['cited_sum']}；总下载 {stats['download_sum']}\n\n")
    fh.write("## 年度分布\n\n| 年 | 篇数 |\n|---|---|\n")
    for y, c in sorted(years.items()):
        fh.write(f"| {y} | {c} |\n")
    fh.write("\n## 发文最多的期刊\n\n| 期刊 | 篇数 |\n|---|---|\n")
    for j, c in journals.most_common(12):
        fh.write(f"| {j} | {c} |\n")
    fh.write("\n## 被引最高的 15 篇\n\n| 年 | 被引 | 题名 | 来源 |\n|---|---|---|---|\n")
    for r in cited:
        fh.write(f"| {year(r['发表时间'])} | {num(r['被引'])} | {r['题名'][:60]} | {r['来源']} |\n")
    fh.write("\n## 下载最高的 15 篇\n\n| 年 | 下载 | 题名 | 来源 |\n|---|---|---|---|\n")
    for r in downs:
        fh.write(f"| {year(r['发表时间'])} | {num(r['下载'])} | {r['题名'][:60]} | {r['来源']} |\n")

print(f"原始 {len(recs)} 行 → 去重 {len(uniq)} 篇；期刊 {len(journals)}；总被引 {stats['cited_sum']}")
print("年度:", dict(sorted(years.items())))
print("Top 期刊:", journals.most_common(5))
print(f"→ {outdir}/cnki_corpus.csv, cnki_summary.md, cnki_stats.json")
