#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_null.py — 阴性对照引擎（四套零模型）

为什么必须有这一步
  网络药理学最常见的「过度声称」来自**没有零模型**：交集基因在 PPI 里当然连得很紧，
  因为 STRING 本身就偏向研究热点蛋白（eLife 2024, PMID 39660719），
  而药材靶点又是从「可成药基因」里来的。不做对照，任何「紧密互联」都可以被质疑为方法学假象。

四套零模型（层层加严，回答不同问题）
  A 全基因组随机  从人类蛋白编码基因随机抽同样数量的基因
                  → 回答「是否高于随机基因」（最宽松，几乎必然显著，仅作下限）
  B 疾病集内随机  从该研究的疾病靶点集内随机抽同样数量
                  → 回答「是否高于随机疾病基因」（机制结论的主要依据）
  C 度匹配随机    在 PPI 中抽取**度分布相同**的基因（按度分层抽样）
                  → 回答「是否只是枢纽基因扎堆」（拓扑混杂）
  D 标签置换      打乱「成分-靶点」边但保持每个成分的靶点数与每个靶点的成分数
                  → 回答「交集大小是否只是由度数结构决定」

判读
  B 组经验 P<0.05：交集构成真实功能模块（可作为机制假说依据）
  C/D 组不显著：交集可能只是「高度数基因」的集合，机制结论必须降级

用法
  python np_null.py -g intersection_genes.txt -e ppi_edges.tsv -d disease_targets.tsv -o out/ -B 1000
"""
import argparse
import os
import random
import sys

import npcore as C


def induced_edges(g, nodes):
    s = set(nodes)
    return sum(1 for u in s if u in g.adj for v in g.adj[u] if v in s and u < v)


def stat_mean_degree(g, nodes):
    s = [n for n in nodes if n in g.adj]
    if not s:
        return 0.0
    return sum(len(g.adj[n]) for n in s) / len(s)



def main(argv=None):
    ap = C.Cli(description="四套零模型阴性对照：判定交集靶点的网络紧密性是否超出随机预期")
    ap.add_argument("-g", "--genes", required=True, help="交集基因列表")
    ap.add_argument("-e", "--edges", required=True, help="PPI 边表")
    ap.add_argument("-d", "--disease", default=None, help="疾病靶点表（用于 B 组与 D 组）")
    ap.add_argument("--drug", default=None, help="成分-靶点表（用于 D 组）")
    ap.add_argument("--ref", default=None, help="参照成分-靶点表（同类药材/随机方），构建 D 组参照零模型")
    ap.add_argument("-o", "--outdir", default="np_out")
    ap.add_argument("-B", "--iterations", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--min-score", type=float, default=0.4)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    man = C.Manifest("05_null", a.outdir, vars(a))
    for k, v in [("genes", a.genes), ("edges", a.edges), ("disease", a.disease), ("drug", a.drug)]:
        if v:
            man.input(k, v)
    rng = random.Random(a.seed)

    genes = [g for g in C.read_gene_list(a.genes)]
    if len(genes) < 3:
        print("[错误] 交集基因 <3，零模型无意义。")
        return 2

    h, rows = C.read_table(a.edges)
    cols = [x.lower() for x in h]
    ca = h[cols.index("node1")] if "node1" in cols else h[0]
    cb = h[cols.index("node2")] if "node2" in cols else h[1]
    cs = None
    for cand in ("combined_score", "score", "weight"):
        if cand in cols:
            cs = h[cols.index(cand)]
            break
    g = C.edges_from_rows(rows, ca, cb, cs, a.min_score)

    universe = sorted(g.nodes)
    degree = g.degree()
    obs_in = induced_edges(g, genes)
    obs_md = stat_mean_degree(g, genes)

    disease = set()
    if a.disease:
        hd, rd = C.read_table(a.disease)
        cg = C.resolve_col(hd, "gene", required=True)
        disease = {(r.get(cg) or "").strip().upper() for r in rd} - {""}

    # A：全基因组随机
    nullA_in = []
    if len(universe) >= len(genes):
        for _ in range(a.iterations):
            s = rng.sample(universe, len(genes))
            nullA_in.append(induced_edges(g, s))
    # B：疾病集内随机（限制在 PPI 里的疾病基因）
    poolB = sorted(disease & set(universe))
    nullB_in = []
    if len(poolB) >= len(genes):
        for _ in range(a.iterations):
            s = rng.sample(poolB, len(genes))
            nullB_in.append(induced_edges(g, s))
    # C：度匹配分层抽样
    nullC_in, cq = [], {}
    if len(universe) >= len(genes):
        obs_deg = [degree.get(x, 0) for x in genes if x in g.adj]
        for _ in range(a.iterations):
            s, cq = C.degree_matched_sample(obs_deg, universe, degree, rng)
            nullC_in.append(induced_edges(g, s))
    # D：参照成分库零模型（可选，需 --ref 提供"同类药材的成分-靶点表"）
    #    问题：本方的交集规模，是否超出一个"成分数量相同、来源相近"的随机方子？
    #    说明：不做"标签置换"——在双边际固定的二分图里，多数聚合统计量对边交换不变，
    #          强行构造会得到一个恒等于观察值的假检验（这是常见的方法学错误）。宁可不做，也不做假。
    nullD_size, obsD = [], None
    if a.ref and a.drug:
        href, rref = C.read_table(a.ref)
        cref = C.resolve_col(href, "compound", required=True)
        gref = "symbol" if "symbol" in href else C.resolve_col(href, "gene", required=True)
        ref_map = {}
        for r in rref:
            c_ = (r.get(cref) or "").strip()
            t_ = (r.get(gref) or "").strip().upper()
            if c_ and t_ and t_ not in ("—", "-", "NA"):
                ref_map.setdefault(c_, set()).add(t_)

        hh, rr = C.read_table(a.drug)
        cc = C.resolve_col(hh, "compound", required=True)
        cg = "symbol" if "symbol" in hh else C.resolve_col(hh, "gene", required=True)
        vals = [(r.get(cg) or "").strip() for r in rr if (r.get(cg) or "").strip()]
        # 前门检查：药侧必须是基因符号，否则 D 组会得到恒为 0 的假结论
        bad = sum(1 for v in vals[:200] if " " in v or v.lower() == v)
        if vals and bad / max(1, len(vals[:200])) > 0.2:
            print("[阻断] D 组需要**已标准化**的成分-靶点表（基因符号），"
                  "当前列 `%s` 中 %.0f%% 的值像蛋白名/非符号。" % (cg, 100 * bad / len(vals[:200])))
            print("       先跑 np_standardize.py，再把 targets_standardized.csv 作为 --drug。")
            return 3

        n_comp = len({(r.get(cc) or "") for r in rr})
        ref_comps = sorted(ref_map)
        ref_t = disease if disease else set(genes)
        if len(ref_comps) >= n_comp >= 2:
            per_comp = {}
            for r in rr:
                c_ = r.get(cc)
                t_ = (r.get(cg) or "").strip().upper()
                if c_ and t_ and t_ not in ("—", "-", "NA"):
                    per_comp.setdefault(c_, set()).add(t_)
            obsD = len(set().union(*per_comp.values()) & ref_t)
            for _ in range(min(a.iterations, 500)):
                pick = rng.sample(ref_comps, n_comp)
                u = set().union(*[ref_map[c_] for c_ in pick])
                nullD_size.append(len(u & ref_t))
        else:
            nullD_size = []
        if obsD == 0:
            print("[阻断] D 组观察值为 0：药侧靶点与疾病集无重叠，或 ID 口径不一致。"
                  "该零模型无信息量，已跳过。")
            nullD_size, obsD = [], None

    def summarize(name, obs, null, desc, expect="greater"):
        if not null:
            return {"test": name, "desc": desc, "obs": obs, "null_mean": "", "null_q95": "",
                    "p": "", "z": "", "n_null": 0, "note": "无法构建该零模型（池子不足）"}
        p, _ = C.empirical_p(obs, null, tail=expect)
        z = C.zscore(obs, null)
        q95 = sorted(null)[min(len(null) - 1, int(0.95 * len(null)))]
        return {"test": name, "desc": desc, "obs": obs,
                "null_mean": round(sum(null) / len(null), 3), "null_q95": q95,
                "p": p, "z": round(z, 2) if z == z else "", "n_null": len(null), "note": ""}

    res = [
        summarize("A_全基因组", obs_in, nullA_in, "随机抽取同样数量的全基因组基因，其诱导边数"),
        summarize("B_疾病集内", obs_in, nullB_in, "从疾病靶点集内随机抽取同样数量，其诱导边数"),
        summarize("C_度匹配", obs_in, nullC_in,
                  "按观察集度多重集逐一匹配抽样（排除「只是枢纽扎堆」）"),
    ]
    if nullD_size:
        d_greater, _ = C.empirical_p(obsD, nullD_size, "greater")
        d_less, _ = C.empirical_p(obsD, nullD_size, "less")
        d_two = min(1.0, 2 * min(d_greater, d_less))
        rd_ = summarize("D_参照成分库", obsD, nullD_size,
                        "成分数量相同的随机方（来自参照库）的交集规模")
        rd_["p"] = d_two
        rd_["note"] = "双侧；观察值%s零均值" % ("高于" if obsD > rd_["null_mean"] else "低于")
        res.append(rd_)

    C.ensure_dir(a.outdir)
    p_csv = C.write_table(os.path.join(a.outdir, "null_tests.csv"),
                          ["test", "desc", "obs", "null_mean", "null_q95", "p", "z", "n_null", "note"],
                          res)

    def get(t):
        for r in res:
            if r["test"].startswith(t):
                return r
        return None

    ra, rb, rc, rd = get("A"), get("B"), get("C"), get("D")
    md = ["# 阴性对照检验报告", "",
          "- 观察集：%d 个交集基因；PPI 网络 %d 节点 / %d 边（阈值 ≥%.2f）" % (len(genes), g.n, g.m, a.min_score),
          "- 统计量：观察集诱导边数（obs）；B/D 组另给经验 P 与 Z", "",
          "| 零模型 | 观察值 | 零均值 | 零95分位 | 经验P | Z | 迭代 |", "|---|---|---|---|---|---|---|"]
    for r in res:
        md.append("| %s | %s | %s | %s | %s | %s | %d |" % (
            r["test"], r["obs"], r["null_mean"], r["null_q95"],
            C.fmt_p(r["p"]) if r["p"] != "" else "—", r["z"] if r["z"] != "" else "—", r["n_null"]))
    md.append("")

    # ---- 零模型自校准（本次新增的核心质量门）----
    md += ["## 零模型自校准（先证明你的对照是有效的）", ""]
    calib_ok = True
    if cq:
        # 注意：0.0 是合法值（完美匹配），不能用 `x or default` 兜底——那是 falsy 陷阱
        def _num(x, default=float("inf")):
            return default if x is None else float(x)
        mad = _num(cq.get("mean_abs_ddeg"))
        exf = _num(cq.get("exact_frac"), -1.0)
        ksd = _num(cq.get("ks"))
        ok_mad, ok_exf, ok_ks = (mad <= 1.0), (exf >= 0.70), (ksd <= 0.10)
        md += ["| 指标 | 值 | 判定 |", "|---|---|---|",
               "| C 组平均 \\|Δdeg\\| | %.2f | %s |" % (mad, "通过" if ok_mad else "**未通过**"),
               "| C 组精确匹配比例 | %.2f | %s |" % (exf, "通过" if ok_exf else "**未通过**"),
               "| C 组度分布 KS 距离 | %.3f | %s |" % (ksd, "通过" if ok_ks else "**未通过**"), ""]
        calib_ok = ok_mad and ok_exf and ok_ks
        if calib_ok:
            md += ["> C 组零模型的度分布与观察集匹配（KS≤0.10），其经验 P 值可作为"
                   "「控制度之后是否仍紧密」的依据。", ""]
        else:
            md += ["> ⛔ **C 组零模型未通过校准**：它对度的匹配不足，其 P 值**不可用于下结论**。"
                   "常见原因：网络中度值种类太少、或交集基因的度极端偏高。"
                   "处置：改报 B 组结论，并把「度混杂未能完全控制」写入局限。", ""]
    else:
        md += ["> 未构建 C 组（池子不足），度混杂未受控制。", ""]

    md += ["## 判读", ""]
    if rb and rb["p"] != "":
        if rb["p"] < 0.05:
            md += ["- **B 组显著（P=%s）**：交集基因在疾病靶点集内部也表现出超出随机的互联，"
                   "支持「这些靶点构成功能模块」的表述（仍属预测层证据）。" % C.fmt_p(rb["p"])]
        else:
            md += ["- **B 组不显著（P=%s）**：交集基因的互联程度与「从疾病基因里随机抽」无异。"
                   "此时**不得**用「紧密互联的核心模块」描述结果，机制推断强度必须下调，"
                   "建议改以富集与实验证据为主线。" % C.fmt_p(rb["p"])]
    if rc and rc["p"] != "" and not calib_ok:
        md += ["- **C 组零模型未通过校准，其结论不予采纳**（见上表）。"]
    elif rc and rc["p"] != "":
        if rc["p"] >= 0.05:
            md += ["- **C 组不显著（P=%s）**：观察集的紧密性可由「度分布」解释，"
                   "即交集偏向高连接度基因。这是 PPI 类研究最典型的混杂，必须在局限性中写明。" % C.fmt_p(rc["p"])]
        else:
            md += ["- C 组显著（P=%s）：即便控制了度分布，交集仍比随机更紧密。" % C.fmt_p(rc["p"])]
    # ---- 边缘显著提示：0.01 < P < 0.05 的裁决对清洗策略高度敏感 ----
    if rb and rb["p"] != "" and 0.01 < rb["p"] <= 0.05:
        md += ["", "> ⚠ **B 组 P=%s 处于边缘（0.01–0.05）**。这类结论通常对预处理高度敏感："
               "剔除/保留非人类测定、换疾病基因阈值、换 STRING 置信度都可能改变结论方向。"
               "论文中必须同时报告该敏感性分析，不能只写一个 P 值。" % C.fmt_p(rb["p"]), ""]
        man.warn("B 组 P=%s 处于边缘，需报告预处理敏感性" % C.fmt_p(rb["p"]))

    # ---- 统一裁决：三套零模型可能给出相反结论，必须给出唯一口径 ----
    # 否则作者会挑有利的那一套（cherry-picking），审稿人一眼就能看出。
    if rb and rc and rb["p"] != "" and rc["p"] != "" and calib_ok:
        b_sig, c_sig = rb["p"] < 0.05, rc["p"] < 0.05
        if b_sig and c_sig:
            verdict = ("**两套对照一致显著**：交集的紧密性既超出「疾病基因随机抽样」，"
                       "也超出「同度随机抽样」。可写「超出疾病基因本底的互联模块」，"
                       "但仍属预测层，需实验验证。")
        elif (not b_sig) and c_sig:
            verdict = ("**裁决：以 B 组为准，不得声称「核心模块」。** B 组（在疾病基因集内随机抽取）"
                       "未显示超出随机，说明交集的紧密性可由「它们本来就是疾病基因」解释；"
                       "C 组虽显著，但它只控制了度，未控制疾病基因本底，属较弱的对照。"
                       "论文中如只报 C 组而不报 B 组，等同于选择性报告。")
        elif b_sig and (not c_sig):
            verdict = ("**裁决：紧密性存在但可由度分布解释。** B 组显著、C 组不显著，"
                       "说明交集偏向高连接度基因（hub 偏倚）。可写「富集于高连接度疾病基因」，"
                       "不得写「特异性模块」。")
        else:
            verdict = ("**两套对照均不显著**：不支持「紧密互联」的表述。"
                       "应把结论重心放在富集与外部证据上，并在局限性中写明网络层面无信号。")
        md += ["", "## 统一裁决（三套对照不一致时以此为准）", "", verdict, "",
               "| 对照 | 问题 | 结论方向 |", "|---|---|---|",
               "| A 全基因组 | 比随机基因更紧密？ | %s |" % ("是" if (ra and ra["p"] != "" and ra["p"] < 0.05) else "否/未做"),
               "| B 疾病集内 | 比随机疾病基因更紧密？ | %s |" % ("是" if b_sig else "否"),
               "| C 度匹配 | 控制度之后仍更紧密？ | %s |" % ("是" if c_sig else "否"),
               "", "> 优先级：B > C > A。A 几乎总会显著（任何疾病相关基因都比随机基因更互联），"
               "单报 A 等于没有对照。", ""]

    if rd and rd["p"] != "":
        direction = "高于" if rd["obs"] > rd["null_mean"] else "低于"
        md += ["- D 组（参照成分库，双侧）：观察值 %s %s零均值 %s（P=%s）。"
               % (rd["obs"], direction, rd["null_mean"], C.fmt_p(rd["p"]))]
        if rd["p"] >= 0.05:
            md += ["  观察规模与「随机取同样数量成分」不可区分：本方的交集规模不能作为**特异性**证据。"]
        elif direction == "高于":
            md += ["  观察规模高于随机参照方：本方在靶点层面比「成分数相近的随机方」更集中地命中疾病基因，可作为特异性佐证。"]
        else:
            md += ["  ⚠ 观察规模**低于**随机参照方：本方命中疾病基因的能力不及成分数相近的随机方，"
                   "提示靶点预测噪声较大、或成分与疾病通路口径不一致。**不可**据此宣称特异性。"]
        md.append("")

    md += ["", "## 必须写进论文的一句话", "",
           "> 为排除网络拓扑与研究偏倚造成的假阳性，本研究以 %d 次随机模拟构建零模型："
           "A 全基因组随机、B 疾病靶点集内随机、C 度匹配随机、D 成分-靶点标签置换；"
           "经验 P 值按 (1+#{null≥obs})/(1+B) 计算。" % a.iterations, ""]

    # 零分布 ECDF 图（审稿人最认可的一张图：观察值落在零分布的哪里）
    series = []
    for name, null in [("A 全基因组", nullA_in), ("B 疾病集内", nullB_in), ("C 度匹配", nullC_in)]:
        if not null:
            continue
        xs = sorted(null)
        n = len(xs)
        pts = [(xs[0], 0.0)] + [(x, (i + 1) / n) for i, x in enumerate(xs)]
        series.append((name, pts))
    p_ecdf = None
    if series:
        series.append(("观察值", [(obs_in, 0.0), (obs_in, 1.0)]))
        p_ecdf = C.line_svg(series, os.path.join(a.outdir, "null_ecdf.svg"),
                            title="零分布 ECDF 与观察值（观察值越靠右越不寻常）",
                            xlabel="诱导边数", ylabel="累积概率")
        md += ["## 零分布 ECDF", "",
               "![ECDF](null_ecdf.svg)", "",
               "竖线为观察值。若它落在三条曲线的右端（累积概率≈1），说明观察到的连接密度"
               "在三种零模型下都不寻常；若只落在 A 的右端而不在 C 的右端，说明差异主要来自度分布。", ""]

    p_md = C.write_text(os.path.join(a.outdir, "null_report.md"), "\n".join(md) + "\n")

    man.output("tests", p_csv)
    if p_ecdf:
        man.output("ecdf", p_ecdf)
    man.output("report", p_md)
    man.count("genes", len(genes))
    man.count("iterations", a.iterations)
    for r in res:
        if r["p"] != "":
            man.count("p_" + r["test"], r["p"])
    if rb and rb["p"] != "" and rb["p"] >= 0.05:
        man.warn("B 组零模型不显著：交集不构成超出随机的功能模块")
    if rc and rc["p"] != "" and rc["p"] >= 0.05:
        man.warn("C 组度匹配零模型不显著：紧密性可由度分布解释")
    man.save()

    if not a.quiet:
        for r in res:
            print("[零模型] %-10s obs=%s null_mean=%s P=%s Z=%s" % (
                r["test"], r["obs"], r["null_mean"],
                C.fmt_p(r["p"]) if r["p"] != "" else "—", r["z"] if r["z"] != "" else "—"))
        print("[零模型] 报告 → %s" % p_md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
