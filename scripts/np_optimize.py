#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_optimize.py — 反向网络药理学：从疾病基因集反推药味/成分，并做组合优化

适用场景
  * 已有一批候选药材（或候选成分），想知道"哪几味对这批疾病基因贡献最大"
  * 交集成百上千、无法用「核心靶点」叙事时，改问"哪味药覆盖了哪些靶点、彼此是否互补"
  * 组方优化：用最少的药味覆盖最多的疾病相关靶点（集合覆盖），并量化协同增量

三个指标（缺一不可，单独用都会被误导）
  coverage    覆盖的疾病基因数（广度）
  precision   coverage / 该药全部靶点数（命中率；低 precision = 广撒网）
  marginal    在已有其他药味的前提下，该药还能新增多少（互补性）——这才是「配伍」的量化

注意
  这里的「优化」是**计算意义上的集合覆盖**，不是药效学最优。它只回答「哪些药在靶点层面互补」，
  不回答剂量、配比、煎煮相互作用。撰写时必须如此限定。

用法
  python np_optimize.py --drug compound_targets.csv --disease disease_targets.tsv --herb-col herb -o out/
  python np_optimize.py --drug x.csv --disease y.tsv --max-herbs 4
"""
import argparse
import os
import sys
from collections import defaultdict

import npcore as C


def main(argv=None):
    ap = C.Cli(description="反向网络药理学：药味/成分覆盖度、互补性与最小覆盖组合")
    ap.add_argument("--drug", required=True, help="成分-靶点表（需含 herb 与 symbol/gene 列）")
    ap.add_argument("--disease", required=True, help="疾病靶点表")
    ap.add_argument("--herb-col", default=None)
    ap.add_argument("--gene-col", default=None)
    ap.add_argument("--unit", choices=["herb", "compound"], default="herb",
                    help="按药味还是按成分做优化，默认 herb")
    ap.add_argument("--max-k", type=int, default=5, help="组合最多包含几个单元")
    ap.add_argument("--min-gain", type=int, default=2, help="边际新增低于此值即停止")
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--kept", default=None,
                    help="上游成分筛选的 kept 名单（一行一个）；给了就只在这些成分里优化，"
                         "避免把已被门控剔除的成分算进「最优组合」")
    ap.add_argument("-o", "--outdir", default="np_out")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    man = C.Manifest("09_optimize", a.outdir, vars(a))
    man.input("drug", a.drug)
    man.input("disease", a.disease)

    h, rows = C.read_table(a.drug)
    cu = a.herb_col or C.resolve_col(h, "herb") or C.resolve_col(h, "compound", required=True)
    cg = a.gene_col or ("symbol" if "symbol" in h else None) or C.resolve_col(h, "gene", required=True)
    if a.unit == "compound":
        cu = a.herb_col or C.resolve_col(h, "compound", required=True)

    kept = None
    if a.kept:
        if not os.path.exists(a.kept):
            print("[错误] --kept 文件不存在：%s" % a.kept)
            return 2
        kept = {x.strip().lower() for x in open(a.kept, encoding="utf-8") if x.strip()}
        man.count("kept_filter", len(kept))

    unit_targets = defaultdict(set)
    n_dropped = 0
    for r in rows:
        u = (r.get(cu) or "").strip()
        g = (r.get(cg) or "").strip().upper()
        if u and g and g not in ("—", "-", "NA"):
            if kept is not None and u.strip().lower() not in kept:
                n_dropped += 1
                continue
            unit_targets[u].add(g)
    if kept is not None:
        man.count("dropped_by_kept", n_dropped)
        if n_dropped and not a.quiet:
            print("[反向] 按 kept 名单剔除 %d 条记录（成分-%d 个单元）"
                  % (n_dropped, len(kept)))
    if not unit_targets:
        print("[错误] 未解析到任何 药味-靶点 关系，检查 --herb-col / --gene-col。")
        return 2

    hd, rd = C.read_table(a.disease)
    cgd = C.resolve_col(hd, "gene", required=True)
    disease = {(r.get(cgd) or "").strip().upper() for r in rd} - {""}

    all_targets = set().union(*unit_targets.values())
    rows_stat = []
    for u, tg in unit_targets.items():
        cov = tg & disease
        rows_stat.append({
            "unit": u, "n_targets": len(tg), "coverage": len(cov),
            "precision": round(len(cov) / max(1, len(tg)), 3),
            "fold": round(len(cov) / max(1e-9, len(tg) * len(disease) / 20596.0), 2),
            "genes": ";".join(sorted(cov)[:20]),
        })
    rows_stat.sort(key=lambda r: (-r["coverage"], -r["precision"]))

    # 边际贡献（leave-one-out）：去掉该单元后，其余并集还能覆盖多少
    union_all = set().union(*unit_targets.values()) & disease
    for r in rows_stat:
        rest = set()
        for u2, tg2 in unit_targets.items():
            if u2 != r["unit"]:
                rest |= (tg2 & disease)
        r["marginal"] = len(union_all - rest)

    # 贪心最小覆盖
    chosen, covered = [], set()
    trace = []
    while len(chosen) < a.max_k:
        best, best_gain, best_prec = None, 0, 0.0
        for r in rows_stat:
            if r["unit"] in chosen:
                continue
            gain = len((unit_targets[r["unit"]] & disease) - covered)
            prec = gain / max(1, len(unit_targets[r["unit"]]))
            if gain > best_gain or (gain == best_gain and prec > best_prec):
                best, best_gain, best_prec = r, gain, prec
        if best is None or best_gain < a.min_gain:
            break
        chosen.append(best["unit"])
        covered |= (unit_targets[best["unit"]] & disease)
        trace.append({"step": len(chosen), "unit": best["unit"], "gain": best_gain,
                      "cumulative": len(covered), "precision": round(best_prec, 3)})

    # 协同增量：组合覆盖 − 各单元单独覆盖的并集（若为 0 则说明无互补）
    union_of_individual = set().union(*[unit_targets[u] & disease for u in chosen]) if chosen else set()
    synergy_gain = len(covered) - len(union_of_individual)

    C.ensure_dir(a.outdir)
    p_rank = C.write_table(os.path.join(a.outdir, "reverse_ranking.csv"),
                           ["unit", "n_targets", "coverage", "precision", "fold", "marginal", "genes"],
                           rows_stat)
    p_comb = C.write_table(os.path.join(a.outdir, "combination_trace.csv"),
                           ["step", "unit", "gain", "cumulative", "precision"], trace)

    if rows_stat:
        p_bar = C.bar_svg([(r["unit"], r["coverage"]) for r in rows_stat[:a.top]],
                          os.path.join(a.outdir, "coverage_bar.svg"),
                          title="各%s覆盖的疾病基因数" % ("药味" if a.unit == "herb" else "成分"),
                          xlabel="覆盖基因数")
    else:
        p_bar = None

    md = ["# 反向网络药理学与组合优化", "",
          "- 疾病基因集：%d 个；候选%s：%d 个；候选靶点并集：%d 个"
          % (len(disease), "药味" if a.unit == "herb" else "成分", len(unit_targets), len(all_targets)),
          "- 全体并集覆盖疾病基因：%d（%.0f%%）" % (len(union_all), 100 * len(union_all) / max(1, len(disease))), ""]

    md += ["## 覆盖度排名（前 %d）" % a.top, "",
           "| %s | 靶点数 | 覆盖 | 命中率 | 富集倍数 | 边际贡献 |" % ("药味" if a.unit == "herb" else "成分"),
           "|---|---|---|---|---|---|"]
    for r in rows_stat[:a.top]:
        md.append("| %s | %d | %d | %.2f | %.2f | %d |" % (
            r["unit"], r["n_targets"], r["coverage"], r["precision"], r["fold"], r["marginal"]))
    md += ["", "> 读法：**coverage 高但 precision 低** = 广撒网（靶点多、命中率低），"
           "这类药味在真实方剂中常是「佐使」而非主药；**marginal 高** = 与其他药味互补性强，"
           "是配伍分析里最有信息量的指标。", ""]

    if chosen:
        md += ["## 贪心最小覆盖组合（计算意义上的最省组合）", "",
               "| 步 | 加入 | 新增覆盖 | 累计覆盖 | 该步命中率 |", "|---|---|---|---|---|"]
        for t in trace:
            md.append("| %d | %s | +%d | %d | %.2f |" % (t["step"], t["unit"], t["gain"],
                                                       t["cumulative"], t["precision"]))
        md += ["", "- 组合：**%s**" % " + ".join(chosen),
               "- 累计覆盖：%d / %d 个疾病基因（%.0f%%）"
               % (len(covered), len(disease), 100 * len(covered) / max(1, len(disease))),
               "- 协同增量（组合覆盖 − 各自单独覆盖之并）：**%d**" % synergy_gain, ""]
        if synergy_gain == 0:
            md += ["> 协同增量为 0：所选单元之间**没有靶点层面的互补**，组合的价值只能来自"
                   "剂量、代谢或通路上下游，不能声称「多靶点协同」。", ""]
        else:
            md += ["> 协同增量为 %d：说明存在「单用覆盖不到、合用才覆盖」的靶点，"
                   "这是配伍协同在靶点层面的可量化证据（仍属预测层）。" % synergy_gain, ""]

    md += ["## 边界（必须写进论文）", "",
           "1. 本优化是**集合覆盖**，不是药效学最优：不涉及剂量、配比、煎煮与代谢相互作用。",
           "2. 覆盖度高不等于疗效好：靶点覆盖只是「可能干预的面」，不等于干预强度。",
           "3. 若所有药味的 precision 都 <0.3，说明靶点预测噪声过大，"
           "应先做 `np_standardize` 与多库交叉，而不是继续优化组合。", ""]

    p_md = C.write_text(os.path.join(a.outdir, "optimize_report.md"), "\n".join(md) + "\n")
    for k, v in [("ranking", p_rank), ("trace", p_comb), ("bar", p_bar), ("report", p_md)]:
        if v:
            man.output(k, v)
    man.count("units", len(unit_targets))
    man.count("disease_genes", len(disease))
    man.count("chosen", chosen)
    man.count("synergy_gain", synergy_gain)
    if synergy_gain == 0 and chosen:
        man.warn("所选组合无靶点层面互补（协同增量 0）")
    man.save()

    if not a.quiet:
        print("[反向] 候选 %d 个，覆盖疾病基因 %d/%d" % (len(unit_targets), len(union_all), len(disease)))
        if rows_stat:
            print("[反向] 覆盖最高：%s（%d）" % (rows_stat[0]["unit"], rows_stat[0]["coverage"]))
        if chosen:
            print("[反向] 最小组合：%s（覆盖 %d，协同增量 %d）"
                  % (" + ".join(chosen), len(covered), synergy_gain))
        print("[反向] 报告 → %s" % p_md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
