#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""drugbank_targets_to_genes.py — 把 DrugBank 靶点表里的"蛋白全名"映射成基因符号。

为什么需要：登录态 drug card 的 Targets 表字段是 `TARGET / ACTIONS / ORGANISM`
（蛋白全名），基因符号与 UniProt AC 只出现在折叠的摘要 chip 里，抽表抽不到。
本脚本用 UniProt REST 做名称→(AC, 基因符号) 的反查，并把每条映射的**来源与匹配状态**
一并落盘——不允许静默猜。

匹配规则（写死，避免"看似成功"）：
  1) protein_name:"<全名>" AND organism_id:9606，取第一条 reviewed 优先；
  2) 命中多条不同基因 → status=ambiguous，保留第一条但标注候选数；
  3) 无命中 → status=unmatched，保留原名，供人工补；
  4) 表里 TARGET 首字符是链接图标字母（U/A/T 等），按"首字母后紧跟非空格"剥离一次。

用法：
  python tools/drugbank_targets_to_genes.py --csv-dir <目录> --out <输出CSV>
  目录内文件名需形如 q_targets.csv，并用 --map 指定 文件→(化合物,DrugBankID)
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

API = "https://rest.uniprot.org/uniprotkb/search"

# 站点被浏览器翻译成中文时出现过的行，人工确认后写死，避免依赖翻译文本做模糊匹配
OVERRIDE = {
    "酪氨酸蛋白磷酸酶非受体类型1": ("PTPN1", "P18031"),
    "糖原磷酸化酶，肌肉形式": ("PYGM", "P11217"),
    "糖原磷酸化酶,肌肉形式": ("PYGM", "P11217"),
    # DrugBank 用旧式复合体命名，UniProt 官方名不同；AC 与 drug card 摘要 chip 一致，可交叉核对
    "ATP synthase F(1) complex subunit alpha, mitochondrial": ("ATP5F1A", "P25705"),
    "ATP synthase F(1) complex subunit beta, mitochondrial": ("ATP5F1B", "P06576"),
    "ATP synthase F(1) complex subunit gamma, mitochondrial": ("ATP5F1C", "P36542"),
    "DNA (cytosine-5)-methyltransferase 3B": ("DNMT3B", "Q9UBC3"),
}


def clean(name: str) -> str:
    s = (name or "").strip()
    # 图标字母：单个大写字母紧接正文（U酪氨酸… / A5'-AMP… / USignal…）
    m = re.match(r"^([A-Z])(?=[\u4e00-\u9fff0-9A-Z])", s)
    if m and not s[1:2].islower():
        s = s[1:]
    s = re.sub(r"^[\s:：]+", "", s)
    return s.strip()


def _query(name: str, reviewed: bool) -> list[list[str]]:
    q = f'protein_name:"{name}" AND organism_id:9606'
    if reviewed:
        q += ' AND reviewed:true'
    url = API + "?" + urllib.parse.urlencode(
        {"query": q, "fields": "accession,gene_primary,protein_name",
         "format": "tsv", "size": "25"})
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "np-skill/3.0.6"})
            with urllib.request.urlopen(req, timeout=30) as r:
                txt = r.read().decode("utf-8", "ignore")
            return [l.split("\t") for l in txt.strip().split("\n")[1:] if l.strip()]
        except Exception:
            time.sleep(2 * (attempt + 1))
    return []


def _exact(rows: list[list[str]], name: str) -> list[list[str]]:
    """只接受"官方蛋白名 == 查询名"的条目（忽略括号补充/同义后缀）。

    必要性：UniProt 的 protein_name 检索会命中别名，直接取第一条曾把
    Amyloid-beta precursor protein 匹配成 APPBP2、PPARγ 匹配成 PPRC1。
    """
    n = name.strip().lower()
    keep = []
    for r in rows:
        if len(r) < 3:
            continue
        full = re.split(r"\s*\(", r[2])[0].strip().lower()
        if full == n:
            keep.append(r)
    return keep


def uniprot_lookup(name: str, cache: dict) -> tuple[str, str, str, int]:
    if name in OVERRIDE:
        g, ac = OVERRIDE[name]
        return g, ac, "override", 1
    if name in cache:
        return tuple(cache[name])  # type: ignore[return-value]

    hits = _exact(_query(name, reviewed=True), name)
    if not hits:
        hits = _exact(_query(name, reviewed=False), name)
    genes = {r[1] for r in hits if len(r) > 1 and r[1]}
    if not hits:
        res = ("", "", "unmatched", 0)
    elif len(genes) == 1:
        g = next(iter(genes))
        ac = next((r[0] for r in hits if len(r) > 1 and r[1] == g), hits[0][0])
        res = (g, ac, "ok", len(hits))
    else:
        res = ("", "", "ambiguous", len(genes))
    cache[name] = res
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache", default="")
    a = ap.parse_args()

    # 文件 → (化合物, DrugBank ID)
    table = {
        "q_targets.csv": ("Quercetin", "DB04216"),
        "c_targets.csv": ("Curcumin", "DB11672"),
        "b_targets.csv": ("Berberine", "DB04115"),
        "u_targets.csv": ("Ursolic acid", "DB15588"),
    }
    cache: dict = {}
    if a.cache and os.path.exists(a.cache):
        cache = {k: tuple(v) for k, v in json.load(open(a.cache, encoding="utf-8")).items()}

    out_rows = []
    for fn, (cpd, dbid) in table.items():
        p = os.path.join(a.csv_dir, fn)
        if not os.path.exists(p):
            print(f"[skip] {fn} 不存在")
            continue
        rows = list(csv.reader(open(p, encoding="utf-8-sig")))
        hdr = [h.strip().upper() for h in rows[0]]
        i_t = hdr.index("TARGET") if "TARGET" in hdr else 0
        # 表头可能被翻译（目标/行动/有机体）
        i_a = next((i for i, h in enumerate(hdr) if h in ("ACTIONS", "行动", "作用")), 1)
        i_o = next((i for i, h in enumerate(hdr) if h in ("ORGANISM", "有机体", "物种")), 2)
        for r in rows[1:]:
            if len(r) <= i_t:
                continue
            raw = r[i_t]
            name = clean(raw)
            act = r[i_a].strip() if len(r) > i_a else ""
            org = r[i_o].strip() if len(r) > i_o else ""
            gene, ac, status, n = uniprot_lookup(name, cache)
            out_rows.append({
                "compound": cpd, "drugbank_id": dbid, "gene": gene, "uniprot": ac,
                "protein_name": name, "action": act, "organism": org,
                "match_status": status, "n_candidates": n,
            })
        print(f"[{cpd}] {len([x for x in out_rows if x['compound'] == cpd])} 条已映射")

    with open(a.out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    if a.cache:
        json.dump({k: list(v) for k, v in cache.items()},
                  open(a.cache, "w", encoding="utf-8"), ensure_ascii=False, indent=0)

    bad = [r for r in out_rows if r["match_status"] in ("unmatched", "ambiguous")]
    print(f"\n总 {len(out_rows)} 条；ok={len(out_rows)-len(bad)}，需人工={len(bad)}")
    for r in bad:
        print(f"  [{r['match_status']}] {r['compound']}: {r['protein_name'][:60]} ({r['organism']})")
    print(f"→ {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
