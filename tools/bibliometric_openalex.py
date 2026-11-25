#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""bibliometric_openalex.py — 用 OpenAlex 做可复现的"竞品 / 空白"计量扫描。

为什么存在：CNKI 的检索接口会返回滑块验证（code -403 + verify/captcha URL），
无法在无浏览器会话下取数。OpenAlex 免 key、可复现、带 DOI，足以回答
"这个架构/这套指标有没有人发过""该主题年产量与期刊分布如何"这两类问题。
中文库（CNKI/万方）覆盖仍是盲区，必须在报告里显式声明。

用法：
  python tools/bibliometric_openalex.py --outdir 全流程实测/bibliometric \
      --mailto you@example.org
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request

BASE = "https://api.openalex.org/works"

# 每条：标签、检索式（title_and_abstract.search）、附加 filter
QUERIES = [
    # OpenAlex 的 search 是 OR 语义，必须写「引号短语 + AND」才得到可解释的计数
    ("np_all", '"network pharmacology"', ""),
    ("np_crc", '"network pharmacology" AND "colorectal cancer"', ""),
    ("np_zh", '"网络药理学"', ""),
    ("np_quality", '"network pharmacology" AND "quality control"', ""),
    ("np_repro", '"network pharmacology" AND "reproducibility"', ""),
    ("np_benchmark", '"network pharmacology" AND "benchmark"', ""),
    ("np_audit", '"network pharmacology" AND "audit"', ""),
    ("np_md", '"network pharmacology" AND "molecular dynamics"', ""),
    ("np_mmpbsa", '"network pharmacology" AND "MM/PBSA"', ""),
    ("np_ml", '"network pharmacology" AND "machine learning"', ""),
]


def get(url: str, tries: int = 3) -> dict:
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "np-bibliometric/1.0"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as exc:
            if i == tries - 1:
                raise
            print(f"  [retry {i+1}] {exc}")
            time.sleep(2 * (i + 1))
    return {}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mailto", default="")
    ap.add_argument("--from-year", type=int, default=2015)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    out = {}
    for tag, terms, extra in QUERIES:
        filt = f"title_and_abstract.search:{terms}"
        if extra:
            filt += f",{extra}"
        q = urllib.parse.quote(filt, safe=":,")
        tail = f"&mailto={a.mailto}" if a.mailto else ""

        try:
            d = get(f"{BASE}?filter={q}&per-page=1{tail}")
        except Exception as exc:
            print(f"{tag:<15} [跳过] {exc}")
            out[tag] = {"query": terms, "total": None, "per_year": {}, "top_sources": [],
                        "top_languages": [], "error": str(exc)}
            continue
        total = d.get("meta", {}).get("count", 0)
        yr = get(f"{BASE}?filter={q}&group_by=publication_year{tail}") if total else {"group_by": []}
        hist = {g["key_display_name"]: g["count"] for g in yr.get("group_by", [])
                if g["key_display_name"].isdigit()}
        recent = {int(k): v for k, v in hist.items() if int(k) >= a.from_year}
        # 期刊名不能直接 group_by（仅 source.id 合法）；改用被引最高的 5 篇反查宿主期刊
        try:
            topw = get(f"{BASE}?filter={q}&sort=cited_by_count:desc&per-page=5{tail}")
            top_src = [(w.get("primary_location", {}).get("source") or {}).get("display_name") or "-"
                       for w in topw.get("results", [])]
        except Exception:
            top_src = []
        try:
            lang = get(f"{BASE}?filter={q}&group_by=language{tail}")
            top_lang = [(g["key_display_name"], g["count"])
                        for g in sorted(lang.get("group_by", []), key=lambda x: -x["count"])[:5]]
        except Exception:
            top_lang = []

        out[tag] = {"query": terms, "total": total,
                    "per_year": dict(sorted(recent.items())),
                    "top_sources": top_src, "top_languages": top_lang}
        peak = max(recent.items(), key=lambda kv: kv[1]) if recent else ("-", 0)
        print(f"{tag:<15} total={total:>6}  峰值年={peak[0]}({peak[1]})  "
              f"{'>'.join(f'{k}:{v}' for k, v in list(sorted(recent.items()))[-4:])}")

    json.dump(out, open(os.path.join(a.outdir, "bibliometric_sweep.json"), "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)

    with open(os.path.join(a.outdir, "bibliometric_sweep.md"), "w", encoding="utf-8") as fh:
        fh.write("| 标签 | 检索式 | 命中总数 | 近五年年产量 | 首位期刊 |\n|---|---|---|---|---|\n")
        for tag, v in out.items():
            recent = v["per_year"]
            last5 = {k: c for k, c in recent.items() if k >= 2021}
            fh.write(f"| {tag} | `{v['query']}` | {v['total']} | "
                     f"{', '.join(f'{k}:{c}' for k, c in sorted(last5.items()))} | "
                     f"{(v['top_sources'][0][:40] if v['top_sources'] else '-')} |\n")
    print(f"\n→ {os.path.join(a.outdir, 'bibliometric_sweep.json')}")
    print("注意：OpenAlex 对 CNKI/万方 中文期刊覆盖不全，中文发文量需另用 CNKI（需浏览器会话 + 滑块验证）。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
