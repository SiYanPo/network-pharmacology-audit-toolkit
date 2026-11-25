#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
stp_convert.py — SwissTargetPrediction 导出 CSV → `np_split compare` 可用的排序表

背景（DoD 收尾第 2 项）
-----------------------
`np_split compare --pred STP=file.csv` 要的是 `gene,score` 两列（大小写/别名自动识别）。
STP 网页"Export to CSV"给的是它自己的列（Target / Probability / Target Class…），
而且**每个成分一个文件**。本脚本负责：

1. 多个成分文件合并；
2. 同一基因的多个概率按 **max**（默认，取"最强证据"）或 **mean**（跨成分共识）聚合；
3. 可选 `--restrict-to` 把排序限制在**与本地工具相同的候选池**内——
   否则覆盖率不同，hits@K 无法比较（这一点比数值本身更重要）；
4. 输出 `stp_ranking.csv`（gene, score, n_compounds, max_prob）与 manifest。

用法
----
    python tools/stp_convert.py --csv "Quercetin=out/stp/Quercetin.csv" \\
        --csv "Tanshinone IIA=out/stp/Tanshinone.csv" -o out/stp
    python tools/stp_convert.py --csv-dir out/stp/raw -o out/stp \\
        --restrict-to np_out/04_intersect/intersection_genes.txt --agg mean

验收
----
    python scripts/np_split.py compare --gold-csv gold.csv --pred STP=out/stp/stp_ranking.csv
"""

import argparse
import csv
import glob
import json
import os
import sys
import time

VERSION = "1.0"

GENE_COLS = ("target", "gene", "symbol", "gene symbol", "target name", "uniprot gene")
# 注意顺序：先精确 "probability"，再兜底 score
PROB_COLS = ("probability", "probability*", "score", "prob", "p")


def _read(path):
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        sample = f.read(4096)
        f.seek(0)
        delim = ";" if sample.count(";") > sample.count(",") else ","
        rows = list(csv.DictReader(f, delimiter=delim))
        hdrs = list(rows[0].keys()) if rows else []
    return rows, hdrs


def _pick(hdrs, cands):
    low = {h.lower().strip(): h for h in hdrs}
    for c in cands:
        if c in low:
            return low[c]
    # 前缀匹配（STP 有时带 * 或括号注释）
    for c in cands:
        for k, v in low.items():
            if k.startswith(c):
                return v
    return None


def parse_one(path):
    rows, hdrs = _read(path)
    if not rows:
        return [], "empty_file"
    gc = _pick(hdrs, GENE_COLS)
    pc = _pick(hdrs, PROB_COLS)
    if not gc:
        return [], "no_gene_column（表头：%s）" % ",".join(hdrs[:6])
    if not pc:
        return [], "no_probability_column（表头：%s）" % ",".join(hdrs[:6])
    out = []
    for r in rows:
        g = (r.get(gc) or "").strip().upper()
        if not g:
            continue
        try:
            p = float(str(r.get(pc) or "0").replace(",", ".").strip())
        except ValueError:
            continue
        out.append((g, p))
    return out, ""


def agg(pairs, mode="max"):
    d = {}
    for g, p in pairs:
        d.setdefault(g, []).append(p)
    out = []
    for g, ps in d.items():
        v = max(ps) if mode == "max" else sum(ps) / len(ps)
        out.append({"gene": g, "score": round(v, 6), "n_compounds": len(ps),
                    "max_prob": round(max(ps), 6)})
    return sorted(out, key=lambda x: (-x["score"], x["gene"]))


def read_list(path):
    s = set()
    if not path or not os.path.exists(path):
        return s
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            t = line.strip().split("\t")[0].split(",")[0].strip().upper()
            if t and t not in ("GENE", "SYMBOL"):
                s.add(t)
    return s


def main(argv=None):
    ap = argparse.ArgumentParser(description="stp_convert: STP 导出 → gene,score 排序表")
    ap.add_argument("--csv", action="append", default=[], help="name=path，可重复")
    ap.add_argument("--csv-dir", default=None, help="目录下所有 *.csv 都读（文件名作来源名）")
    ap.add_argument("--agg", default="max", choices=["max", "mean"])
    ap.add_argument("--restrict-to", default=None,
                    help="候选池基因列表；给定后排序只保留池内基因（保证 hits@K 可比）")
    ap.add_argument("-o", "--outdir", default="out/stp")
    a = ap.parse_args(argv)

    sources = []
    for x in a.csv:
        if "=" not in x:
            print("[stp] --csv 需要 name=path 形式：%s" % x)
            return 2
        k, v = x.split("=", 1)
        sources.append((k, v))
    if a.csv_dir:
        for p in sorted(glob.glob(os.path.join(a.csv_dir, "*.csv"))):
            sources.append((os.path.splitext(os.path.basename(p))[0], p))
    if not sources:
        print("[stp] 需要 --csv 或 --csv-dir")
        return 2

    all_pairs, detail, bad = [], [], []
    for name, path in sources:
        pairs, err = parse_one(path)
        detail.append({"source": name, "file": path, "n_genes": len(pairs), "error": err})
        if err:
            bad.append("%s: %s" % (name, err))
        all_pairs.extend(pairs)
    if not all_pairs:
        print("[stp] 没有解析出任何基因：\n  " + "\n  ".join(bad))
        return 3

    ranked = agg(all_pairs, a.agg)
    pool = read_list(a.restrict_to)
    n_before = len(ranked)
    if pool:
        ranked = [r for r in ranked if r["gene"] in pool]
    os.makedirs(a.outdir, exist_ok=True)
    out = os.path.join(a.outdir, "stp_ranking.csv")
    with open(out, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "score", "n_compounds", "max_prob"])
        w.writeheader()
        w.writerows(ranked)
    man = {"tool": "SwissTargetPrediction", "agg": a.agg,
           "n_genes_before_restrict": n_before, "n_genes_after_restrict": len(ranked),
           "restrict_to": a.restrict_to or "",
           "n_compounds": len(sources), "sources": detail,
           "output": out, "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
           "note": ("STP 的概率是\"被预测为该靶点的可能性\"，与本地打分不同尺度；"
                    "比较时只能用排序类指标（hits@K、EF@K、AUC），不能直接比分数值。")}
    json.dump(man, open(os.path.join(a.outdir, "stp_convert_manifest.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)

    print("[stp] %d 个成分 → %d 个基因（限池后 %d）→ %s" % (
        len(sources), n_before, len(ranked), out))
    if bad:
        print("  解析失败：\n  " + "\n  ".join(bad))
    print("  下一步：python scripts/np_split.py compare --gold-csv <gold.csv> --pred STP=%s" % out)
    return 0 if ranked else 3


if __name__ == "__main__":
    sys.exit(main())
