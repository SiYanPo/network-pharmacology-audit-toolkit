#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""drugbank_known_target_audit.py — 用 DrugBank 公开药品卡片的已标注靶点做"已知靶点召回"盲审。

背景与边界（必须随结果一起报告）：
- 数据来源是 go.drugbank.com 公开药品卡片（drug card）的 Targets / Transporters 区块，
  每条含基因符号 + UniProt ID + bio_entity 链接，可逐条溯源。
- 公开卡片对靶点列表做截断展示：Curcumin 显示 8 条 / 实际 21 条，Quercetin 8 / 32，
  Berberine 8 / 12（另有 1 条细菌靶点 ftsZ 不计入人源评估）。因此本脚本的 recall
  分母只取"可见子集"，是**下界意义上不完整**的估计：若命中，则必然是真阳性；
  若未命中，可能是预测缺失，也可能是被截断隐藏。结论必须按此口径表述。
- 不做任何 ID 映射推断：基因符号大小写不敏感精确匹配，UniProt 仅作记录。

产出口径：
  recall_visible  = 命中的可见靶点数 / 可见靶点数
  support_rate    = 预测靶点中，能在"本批 5 化合物 DrugBank 标注并集"里找到的比例（精确率代理）
  transporter_hit = 转运体是否在预测集中（用于解释 ADMET 的 P-gp/BCRP/OCT 结果）

用法：
  python tools/drugbank_known_target_audit.py \
      --curated 全流程实测/drugbank_audit/drugbank_curated_2026-11-23.csv \
      --np-out  全流程实测/np_out \
      --outdir  全流程实测/drugbank_audit
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys


def norm(g: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", (g or "").upper())


def read_genes(path: str) -> set[str]:
    """读一行一个基因符号的列表；容忍逗号/制表符分隔与表头。"""
    if not os.path.exists(path):
        return set()
    out: set[str] = set()
    with open(path, encoding="utf-8-sig", errors="ignore") as fh:
        for i, line in enumerate(fh):
            for tok in re.split(r"[,;\t\s]+", line.strip()):
                if not tok:
                    continue
                if i == 0 and norm(tok) in {"GENE", "SYMBOL", "GENES", "TARGET"}:
                    continue
                out.add(norm(tok))
    return out


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def load_curated(path: str) -> dict:
    rows = list(csv.DictReader(open(path, encoding="utf-8-sig")))
    db: dict[str, dict] = {}
    for r in rows:
        c = r["compound"].strip()
        d = db.setdefault(c, {
            "drugbank_id": r.get("drugbank_id", "").strip(),
            "group": r.get("group", "").strip(),
            "n_total": int(r["n_targets_total"]) if r.get("n_targets_total") else 0,
            "n_visible": int(r["n_visible"]) if r.get("n_visible") else 0,
            "targets": [], "transporters": [], "source": r.get("source_url", "").strip(),
            "absent": (r.get("kind", "").strip() == "absent"),
        })
        gene = norm(r.get("gene", ""))
        if not gene:
            continue
        kind = r.get("kind", "").strip()
        if kind == "target":
            d["targets"].append((gene, r.get("uniprot", "").strip()))
        elif kind == "transporter":
            d["transporters"].append((gene, r.get("uniprot", "").strip()))
    return db


def main() -> int:
    ap = argparse.ArgumentParser(description="DrugBank 已知靶点召回盲审")
    ap.add_argument("--curated", required=True)
    ap.add_argument("--np-out", required=True, help="np_out 目录（内含 <Compound>/02_standardize/targets_gene_list.txt）")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--compound-map", default="",
                    help="可选 CSV: drugbank_name,pipeline_dir，用于名称不一致时映射")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    curated = load_curated(a.curated)

    name_map: dict[str, str] = {}
    if a.compound_map and os.path.exists(a.compound_map):
        for r in csv.DictReader(open(a.compound_map, encoding="utf-8-sig")):
            name_map[r["drugbank_name"].strip()] = r["pipeline_dir"].strip()

    # 全部 5 化合物的 DrugBank 标注并集（用于 support_rate 这个精确率代理）
    univ: set[str] = set()
    for d in curated.values():
        univ.update(g for g, _ in d["targets"])

    per_compound: dict[str, dict] = {}
    for cname, d in curated.items():
        pdir = name_map.get(cname, cname)
        base = os.path.join(a.np_out, pdir)
        # 目录名常把空格写成下划线（Ursolic_acid），做一次别名回退
        if not os.path.isdir(base) and os.path.isdir(os.path.join(a.np_out, pdir.replace(" ", "_"))):
            pdir = pdir.replace(" ", "_")
            base = os.path.join(a.np_out, pdir)
        pred = read_genes(os.path.join(base, "02_standardize", "targets_gene_list.txt"))
        inter_path = os.path.join(base, "04_intersect", "intersection_genes.txt")
        inter = read_genes(inter_path)

        rec = {
            "compound": cname,
            "pipeline_dir": pdir,
            "drugbank_id": d["drugbank_id"],
            "drugbank_group": d["group"],
            "drugbank_targets_total": d["n_total"],
            "drugbank_targets_visible": len(d["targets"]),
            "in_drugbank": not d["absent"],
            "predicted_targets": len(pred),
            "disease_intersection_targets": len(inter),
            "source_url": d["source"],
        }
        if d["absent"] or not d["targets"]:
            rec.update({"recovered": [], "missed": [], "recall_visible": None,
                        "recall_ci95": None, "transporter_hit": None,
                        "support_in_union": None,
                        "note": "该化合物在 DrugBank 公开库无实体（或不可检索），无法做已知靶点召回"})
            per_compound[cname] = rec
            continue

        seen, uniq = set(), []
        for g, u in d["targets"]:
            if g not in seen:
                seen.add(g)
                uniq.append((g, u))
        hit = [(g, u) for g, u in uniq if g in pred]
        miss = [(g, u) for g, u in uniq if g not in pred]
        lo, hi = wilson(len(hit), len(uniq))
        rec.update({
            "recovered": [g for g, _ in hit],
            "missed": [g for g, _ in miss],
            "recall_visible": round(len(hit) / len(uniq), 4),
            "recall_ci95": [round(lo, 4), round(hi, 4)],
            "transporter_hit": ([g for g, _ in d["transporters"] if g in pred]
                                if d["transporters"] else []),
            "transporters_visible": [g for g, _ in d["transporters"]],
            "support_in_union": (round(len(pred & univ) / len(pred), 4) if pred else None),
            "note": "",
        })
        per_compound[cname] = rec

    # ---- 汇总 ----
    rows = []
    for c, r in per_compound.items():
        rows.append({
            "compound": c, "drugbank_id": r["drugbank_id"], "group": r["drugbank_group"],
            "db_targets_total": r["drugbank_targets_total"],
            "db_targets_visible": r["drugbank_targets_visible"],
            "predicted": r["predicted_targets"],
            "recovered": len(r["recovered"]) if r["recovered"] is not None else "",
            "recall_visible": "" if r["recall_visible"] is None else f"{r['recall_visible']:.3f}",
            "recall_ci95": "" if r["recall_ci95"] is None else f"{r['recall_ci95'][0]:.2f}-{r['recall_ci95'][1]:.2f}",
            "recovered_genes": "|".join(r["recovered"]) if r["recovered"] else "",
            "missed_genes": "|".join(r["missed"]) if r["missed"] else "",
            "transporter_hits": "|".join(r["transporter_hit"]) if r["transporter_hit"] else "",
        })
    with open(os.path.join(a.outdir, "drugbank_recall.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    ev = [r for r in per_compound.values()
          if r["recall_visible"] is not None and r["drugbank_targets_visible"]]
    K = sum(len(r["recovered"]) for r in ev)
    N = sum(r["drugbank_targets_visible"] for r in ev)
    pooled_lo, pooled_hi = wilson(K, N)

    summary = {
        "generated_from": "go.drugbank.com public drug cards (Targets / Transporters blocks)",
        "caveat": ("公开卡片对靶点列表做截断（Curcumin 8/21、Quercetin 8/32、Berberine 8/12）；"
                   "故 recall 的分母是被展示子集，命中即真阳性，未命中可能是截断隐藏。"),
        "compounds_in_drugbank": sum(1 for r in per_compound.values() if r["in_drugbank"]),
        "pooled_recovered": K,
        "pooled_visible_targets": N,
        "pooled_recall_visible": round(K / N, 4) if N else None,
        "pooled_recall_ci95": [round(pooled_lo, 4), round(pooled_hi, 4)],
        "per_compound": per_compound,
    }
    with open(os.path.join(a.outdir, "drugbank_recall.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False)

    # ---- 控制台表 ----
    print(f"{'compound':<15}{'DB ID':<10}{'vis/tot':>9}{'pred':>6}{'hit':>4}{'recall':>8}{'95%CI':>14}")
    for c, r in per_compound.items():
        if r["recall_visible"] is None:
            print(f"{c:<15}{(r['drugbank_id'] or '-'):<10}{'-':>9}{r['predicted_targets']:>6}"
                  f"{'-':>4}{'n/a':>8}{'not in DrugBank':>14}")
            continue
        ci = r["recall_ci95"]
        print(f"{c:<15}{r['drugbank_id']:<10}"
              f"{r['drugbank_targets_visible']:>4}/{r['drugbank_targets_total']:<4}"
              f"{r['predicted_targets']:>6}{len(r['recovered']):>4}"
              f"{r['recall_visible']:>8.3f}   [{ci[0]:.2f},{ci[1]:.2f}]")
    print(f"\npooled: {K}/{N} = {K / N:.3f}  95%CI [{pooled_lo:.2f},{pooled_hi:.2f}]" if N else "\npooled: n/a")
    print(f"→ {os.path.join(a.outdir, 'drugbank_recall.csv')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
