#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_dock.py — 分子对接结果批量评估与门禁

三条硬门禁（不过就不许写结论）
  G1 重对接验证：共晶配体回对接 RMSD ≤ 2.0 Å 才说明该靶点的对接协议可复现晶体结合模式。
     文献基准（J Chem Inf Model 2024, PMID 38108627）：全 PDB 自动流程 top-1 成功率仅 30.1%，
     经结构质量过滤后 61.8%。即「能对接」远不等于「对得准」，逐靶点验证不可省。
  G2 打分口径：结合能是**预测值**，不是实测亲和力（不能等同 Kd/IC50）。
     前瞻性基准里 docking 的 auROC 可以接近 0.5（Mol Syst Biol 2022, PMID 36065847）。
  G3 覆盖说明：只对接了 top N 个成分时，必须说明选择依据（否则属选择性报告）。

用法
  python np_dock.py -i docking_scores.csv -o out/
  python np_dock.py -i vina_out/ -o out/          # 目录模式：解析 *.log / *.txt 中的 affinity
"""
import argparse
import os
import re
import sys
from collections import defaultdict

import npcore as C


def grade(aff):
    if aff is None:
        return "未评分"
    if aff <= -9.0:
        return "很强"
    if aff <= -7.0:
        return "较强"
    if aff <= -5.0:
        return "弱结合"
    return "无结合提示"


def parse_vina_dir(root):
    rows = []
    pat = re.compile(r"^\s*(\d+)\s+(-?\d+\.?\d*)\s+", re.M)
    for dirpath, _dirs, files in os.walk(root):
        for fn in files:
            if not fn.lower().endswith((".txt", ".log", ".pdbqt")):
                continue
            p = os.path.join(dirpath, fn)
            try:
                txt = open(p, encoding="utf-8", errors="replace").read()
            except Exception:
                continue
            m = pat.findall(txt)
            if not m:
                continue
            best = min(float(x[1]) for x in m)
            rows.append({"compound": os.path.splitext(fn)[0], "target": os.path.basename(dirpath),
                         "affinity": best, "redock_rmsd": ""})
    return rows


def main(argv=None):
    ap = C.Cli(description="对接结果批量评估：打分分级 + 重对接门禁 + 选择性报告检查")
    ap.add_argument("-i", "--input", required=True, help="CSV 表（compound,target,affinity[,redock_rmsd]）或 Vina 输出目录")
    ap.add_argument("-o", "--outdir", default="np_out")
    ap.add_argument("--rmsd-cut", type=float, default=2.0, help="重对接 RMSD 门限 Å，默认 2.0")
    ap.add_argument("--top", type=int, default=None, help="每个靶点只保留前 N 个（默认全部）")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    man = C.Manifest("08_dock", a.outdir, vars(a))
    man.input("input", a.input)

    if os.path.isdir(a.input):
        rows = parse_vina_dir(a.input)
        if not rows:
            print("[错误] 目录中未找到可解析的 Vina 输出（含 affinity 表的 .txt/.log）。")
            return 2
        heads = ["compound", "target", "affinity", "redock_rmsd"]
    else:
        h, rows = C.read_table(a.input)
        cc = C.resolve_col(h, "compound", required=True)
        ct = C.resolve_col(h, "term", None, extra=["target", "gene", "protein", "receptor", "靶点"])
        ca = C.resolve_col(h, "score", None, extra=["affinity", "binding_energy", "vina", "结合能", "kcal"])
        # 注意：不能用 key="gene" —— ALIASES 里 "gene" 与 "target" 同义，
        # 会把「靶点名」当成 RMSD 列读进来，于是门禁永远拿不到数值（静默失效）。
        cr = C.resolve_col(h, "rmsd", None,
                           extra=["redock_rmsd", "redock", "重对接rmsd", "重对接_rmsd"])
        if not ca:
            print("[错误] 未找到打分列（affinity/binding_energy/score）。现有列：%s" % "、".join(h))
            return 2
        norm = []
        for r in rows:
            norm.append({"compound": (r.get(cc) or "").strip(),
                         "target": (r.get(ct) or "target").strip() if ct else "target",
                         "affinity": C.to_float(r.get(ca)),
                         "redock_rmsd": C.to_float(r.get(cr)) if cr else ""})
        rows = norm

    # 去掉明显是正值的打分（对接软件异常的常见表现）
    bad = [r for r in rows if r["affinity"] is not None and r["affinity"] > 0]
    rows = [r for r in rows if r["affinity"] is not None and r["affinity"] <= 0]

    by_target = defaultdict(list)
    for r in rows:
        r["grade"] = grade(r["affinity"])
        by_target[r["target"]].append(r)

    # G1 重对接门禁
    gate = {}
    for t, lst in by_target.items():
        rms = [r["redock_rmsd"] for r in lst if isinstance(r["redock_rmsd"], float)]
        gate[t] = {"rmsd": max(rms) if rms else None,
                   "passed": (max(rms) <= a.rmsd_cut) if rms else None,
                   "n": len(lst),
                   "has_col": bool(cr) if a.input and not os.path.isdir(a.input) else False}
    if not any(g["rmsd"] is not None for g in gate.values()):
        man.warn("输入表没有可解析的 RMSD 数值：重对接门禁**未生效**，本报告只能作为打分描述")

    for t, lst in by_target.items():
        lst.sort(key=lambda r: r["affinity"])
        if a.top:
            for i, r in enumerate(lst, 1):
                r["rank_in_target"] = i
            by_target[t] = lst[:a.top]
        else:
            for i, r in enumerate(lst, 1):
                r["rank_in_target"] = i
    flat = [r for lst in by_target.values() for r in lst]
    flat.sort(key=lambda r: (r["target"], r["affinity"]))

    C.ensure_dir(a.outdir)
    p_csv = C.write_table(os.path.join(a.outdir, "docking_graded.csv"),
                          ["compound", "target", "affinity", "grade", "rank_in_target", "redock_rmsd",
                           "target_gate"],
                          [{**r, "target_gate": ("通过" if gate[r["target"]]["passed"] else
                                                 ("未验证" if gate[r["target"]]["passed"] is None else "未通过"))}
                           for r in flat])

    if flat:
        best = sorted(flat, key=lambda r: r["affinity"])[:15]
        p_bar = C.bar_svg([("%s–%s" % (r["compound"][:18], r["target"]), -r["affinity"]) for r in best],
                          os.path.join(a.outdir, "docking_bar.svg"),
                          title="结合能（−kcal/mol，越长越强；仅为预测）", xlabel="−ΔG (kcal/mol)")
    else:
        p_bar = None

    n_total = len(rows) + len(bad)
    md = ["# 分子对接结果评估", "",
          "| 项 | 值 |", "|---|---|",
          "| 记录数 | %d |" % n_total,
          "| 剔除（打分 >0，异常） | %d |" % len(bad),
          "| 靶点数 | %d |" % len(by_target),
          "| 重对接验证通过 | %d |" % sum(1 for v in gate.values() if v["passed"] is True),
          "| 未做重对接验证 | %d |" % sum(1 for v in gate.values() if v["passed"] is None),
          "| 重对接未通过 | %d |" % sum(1 for v in gate.values() if v["passed"] is False), ""]

    failed = [t for t, v in gate.items() if v["passed"] is False]
    if failed:
        md += ["## ⛔ 门禁未通过：%s" % "、".join(failed), "",
               "这些靶点的重对接 RMSD > %.1f Å，说明口袋定义、质子化或受体处理有问题。"
               "**该靶点上的打分不得用于机制论证**，应先修正协议再重跑。" % a.rmsd_cut, ""]

    novalidate = [t for t, v in gate.items() if v["passed"] is None]
    if novalidate:
        md += ["## ⚠ 未做重对接验证：%s" % "、".join(novalidate[:10]), "",
               "无共晶配体（apo 结构）时无法做重对接，属客观限制。此时应："
               "① 用同家族共晶结构做替代验证；② 或明确声明该靶点结果为探索性。", ""]

    if flat:
        md += ["## 打分分布", "", "| 分级 | 条数 |", "|---|---|"]
        cnt = defaultdict(int)
        for r in flat:
            cnt[r["grade"]] += 1
        for k in ["很强", "较强", "弱结合", "无结合提示"]:
            if cnt.get(k):
                md.append("| %s | %d |" % (k, cnt[k]))
        md += ["", "## 最强的 12 对（仅表示预测结合可能性）", "",
               "| 成分 | 靶点 | 结合能 | 分级 | 靶点内排名 | 门禁 |", "|---|---|---|---|---|---|"]
        for r in sorted(flat, key=lambda x: x["affinity"])[:12]:
            md.append("| %s | %s | %.2f | %s | %d | %s |" % (
                r["compound"], r["target"], r["affinity"], r["grade"], r["rank_in_target"],
                "通过" if gate[r["target"]]["passed"] else
                ("未验证" if gate[r["target"]]["passed"] is None else "未通过")))
        md.append("")

    md += ["## 判读与写作规范（逐条对照）", "",
           "1. **打分 ≠ 亲和力**：−8.5 kcal/mol 不能写成「结合亲和力强」，只能写「预测结合能较低」。",
           "2. **必须报口袋定义方式**：共晶配体中心 / 指定残基 / 盲对接，三者结论强度不同。",
           "3. **必须报质子化与辅因子处理**：pH 7.4 下的质子化状态会改变氢键网络。",
           "4. **只对接部分成分时必须说明选择依据**：按丰度、按前序分析、还是按文献；"
           "否则属选择性报告。",
           "5. **阳性对照**：同一协议下应有已知活性分子（或共晶配体）作参照，否则打分无可比性。",
           "6. **要主张「结合」需实验**：SPR / MST / ITC；要主张「药效」需 rescue 实验。", ""]

    if n_total and len(by_target) and len(flat) / max(1, len(by_target)) < 2:
        md += ["> ⚠ 平均每靶点不足 2 个配体：无法做靶点内排序，"
               "「最优成分」的说法缺乏比较基础。", ""]

    p_md = C.write_text(os.path.join(a.outdir, "docking_report.md"), "\n".join(md) + "\n")

    man.output("graded", p_csv)
    man.output("report", p_md)
    if p_bar:
        man.output("bar", p_bar)
    man.count("records", n_total)
    man.count("targets", len(by_target))
    if failed:
        man.warn("重对接门禁未通过：%s" % "、".join(failed))
    if novalidate:
        man.warn("未做重对接验证的靶点：%s" % "、".join(novalidate[:8]))
    man.save()

    if not a.quiet:
        print("[对接] %d 条记录，%d 个靶点；很强 %d，较强 %d"
              % (n_total, len(by_target),
                 sum(1 for r in flat if r["grade"] == "很强"),
                 sum(1 for r in flat if r["grade"] == "较强")))
        if failed:
            print("[对接] ⛔ 重对接未通过：%s" % "、".join(failed))
        print("[对接] 报告 → %s" % p_md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
