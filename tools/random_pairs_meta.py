#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
random_pairs_meta.py —— 随机"成分×疾病"对的阴性对照元实验（方法学证据）

科学问题
    当前网络药理学普遍以"交集靶点显著富集于 PPI 紧密模块 / GO-KEGG 通路"作为
    机制证据。但若随机选取成分与疾病也能频繁得到"显著"，则说明该领域存在
    **背景校准缺失**——大量论文的"核心模块"可能只是 PPI 度分布与通路规模的
    必然产物，而非真实生物学信号。本脚本用 skill 自带 np_null 的四套零模型
    对 N 组随机药-病对做批量检验，统计假阳性率。

流程
    从成分-靶点表随机抽 m 个成分基因集、从疾病表随机抽 d 个疾病基因集，
    求交集 → 调 np_null（B 组疾病内 / C 组度匹配）→ 收集经验 P，
    统计有多少比例 random pairs 在 B 或 C 组 P<0.05（应为假阳性）。

用法
    python tools/random_pairs_meta.py --targets compound_targets.csv \
        --disease disease_targets.tsv --edges string_edges.tsv \
        --n 30 --bgene 20 --dgene 100 -o out/meta/
"""
import argparse
import csv
import json
import os
import random
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))


def read_genes_from_targets(path, field="gene"):
    genes = set()
    with open(path, encoding="utf-8-sig") as f:
        first = f.readline()
        delim = "\t" if "\t" in first else ","
        f.seek(0)
        for r in csv.DictReader(f, delimiter=delim):
            g = (r.get(field) or "").strip().upper()
            if g:
                genes.add(g)
    return sorted(genes)


def parse_p(report_path):
    """从 null_report.md 解析 B/C 组经验 P。"""
    if not os.path.exists(report_path):
        return {}
    txt = open(report_path, encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"\| ([ABC])_[^|]*\|[^|]*\|[^|]*\|[^|]*\|\s*([\d.]+)\s*\|", txt):
        out[m.group(1) + "_p"] = float(m.group(2))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="随机药-病对阴性对照元实验")
    ap.add_argument("--targets", required=True)
    ap.add_argument("--disease", required=True)
    ap.add_argument("--edges", required=True)
    ap.add_argument("--n", type=int, default=30, help="随机对数")
    ap.add_argument("--bgene", type=int, default=20, help="每对抽成分基因数")
    ap.add_argument("--dgene", type=int, default=100, help="每对抽疾病基因数")
    ap.add_argument("--mode", choices=["independent", "pseudo-drug"], default="pseudo-drug",
                    help="independent=两池独立随机抽（交集期望≈0，不适合小基因集）；"
                         "pseudo-drug=从疾病基因池内抽伪成分集（交集=bgene，"
                         "用于检验'从疾病集随机抽基因也显著紧密'的假阳性，推荐）")
    ap.add_argument("--iter", type=int, default=200, help="np_null 迭代数（轻量）")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("-o", "--outdir", default="np_out/random_pairs_meta")
    a = ap.parse_args(argv)

    import np_null
    os.makedirs(a.outdir, exist_ok=True)
    drug_pool = read_genes_from_targets(a.targets, "gene")
    dis_pool = read_genes_from_targets(a.disease, "gene")
    if len(drug_pool) < a.bgene or len(dis_pool) < a.dgene:
        print("[meta] 基因池不足：drug=%d disease=%d" % (len(drug_pool), len(dis_pool)))
        return 2
    rng = random.Random(a.seed)
    rows = []
    for i in range(1, a.n + 1):
        if a.mode == "pseudo-drug":
            # 伪成分集从疾病基因池内抽 → 交集规模=bgene，与真实成分-疾病交集可比。
            # 检验命题：从疾病集随机抽 bgene 个基因，是否也显得"紧密互联"（假阳性）。
            ds = rng.sample(dis_pool, min(a.dgene, len(dis_pool)))
            inter = sorted(rng.sample(ds, min(a.bgene, len(ds))))
        else:
            dg = rng.sample(drug_pool, a.bgene)
            ds = rng.sample(dis_pool, a.dgene)
            inter = sorted(set(dg) & set(ds))
        pair_dir = os.path.join(a.outdir, "pair_%02d" % i)
        os.makedirs(pair_dir, exist_ok=True)
        gf = os.path.join(pair_dir, "genes.txt")
        open(gf, "w", encoding="utf-8").write("\n".join(inter) + "\n")
        if len(inter) < 3:
            rows.append({"pair": i, "n_inter": len(inter), "A_p": "", "B_p": "", "C_p": "",
                         "verdict": "INSUFFICIENT"})
            continue
        try:
            np_null.main(["-g", gf, "-e", a.edges, "-d", a.disease,
                          "-o", pair_dir, "-B", str(a.iter), "--seed", str(a.seed + i), "--quiet"])
        except SystemExit:
            pass
        ps = parse_p(os.path.join(pair_dir, "null_report.md"))
        b, c = ps.get("B_p"), ps.get("C_p")
        sig = any(isinstance(x, float) and x < 0.05 for x in (b, c))
        rows.append({"pair": i, "n_inter": len(inter),
                     "A_p": ps.get("A_p", ""), "B_p": b if b is not None else "",
                     "C_p": c if c is not None else "",
                     "verdict": "SIGNIFICANT" if sig else "NULL"})
        if i % 5 == 0:
            print("  %d/%d done" % (i, a.n), flush=True)

    with open(os.path.join(a.outdir, "random_pairs_results.csv"), "w",
              encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pair", "n_inter", "A_p", "B_p", "C_p", "verdict"])
        w.writeheader(); w.writerows(rows)
    eval_rows = [r for r in rows if r["verdict"] != "INSUFFICIENT"]
    n_sig = sum(1 for r in eval_rows if r["verdict"] == "SIGNIFICANT")
    fpr = n_sig / len(eval_rows) if eval_rows else 0.0
    summary = {"n_pairs": a.n, "n_evaluable": len(eval_rows), "n_significant": n_sig,
               "false_positive_rate": round(fpr, 3),
               "note": ("pseudo-drug 模式：从疾病基因池随机抽 bgene 个基因作伪成分集，"
                        "检验随机子集是否也判为'紧密互联'；若假阳性率远高于 0.05，"
                        "说明网络紧密性判据存在背景校准缺失，不能单独作为机制证据。")
                       if a.mode == "pseudo-drug" else
                       "independent 模式：两池独立随机抽，交集期望≈0（小基因集下多不可评估）。",
               "params": {"mode": a.mode, "bgene": a.bgene, "dgene": a.dgene,
                          "iter": a.iter, "seed": a.seed}}
    json.dump(summary, open(os.path.join(a.outdir, "meta_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[meta] %d/%d 随机对判为显著（假阳性率 %.1f%%）→ %s"
          % (n_sig, len(eval_rows), fpr * 100, a.outdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
