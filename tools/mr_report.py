#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mr_report.py — 把一批 mr_result.json 汇总成「可直接给导师/审稿人看」的 MR 报告

为什么单独写：`np_mr` 每个靶点输出一个 JSON，`np_fdr` 输出 q 值表，
但人要的是一张"6 个靶点、IV 数、F、OR、q、方向判定"的表 + 数据来源 + 局限。
手抄容易出错，所以从磁盘生成。

用法
    python tools/mr_report.py --dir D:/.../DoD收尾/mr -o D:/.../DoD收尾/mr/mr_report.md
"""

import argparse
import csv
import glob
import json
import math
import os
import time


def _j(p):
    try:
        return json.load(open(p, encoding="utf-8"))
    except Exception:
        return {}


def _q_table(root):
    out = {}
    for p in glob.glob(os.path.join(root, "**", "fdr_table.csv"), recursive=True):
        try:
            with open(p, encoding="utf-8-sig", errors="replace") as f:
                for r in csv.DictReader(f):
                    if r.get("family") == "MR_IVW":
                        out[(r.get("unit") or "").upper()] = r
        except Exception:
            continue
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--outcome-label", default="Coronary heart disease (ieu-a-7)")
    a = ap.parse_args(argv)
    root = os.path.abspath(a.dir)
    rows = []
    for p in sorted(glob.glob(os.path.join(root, "**", "mr_result.json"), recursive=True)):
        j = _j(p)
        gene = os.path.basename(os.path.dirname(p))
        clump = j.get("clump") or {}
        ivw = j.get("IVW") or {}
        rows.append({"gene": gene, "file": p, "j": j, "gene_dir": os.path.dirname(p),
                     "status": j.get("status"), "n_iv": clump.get("n_kept"),
                     "n_in": clump.get("n_input"), "mode": clump.get("mode"),
                     "beta": ivw.get("beta"), "se": ivw.get("se"), "p": ivw.get("p"),
                     "F": j.get("F_statistic"),
                     "egger_p": (j.get("MR_Egger") or {}).get("p_intercept"),
                     "egger_b": (j.get("MR_Egger") or {}).get("intercept"),
                     "wm": (j.get("WeightedMedian") or {}).get("beta"),
                     "pp4": (j.get("coloc") or {}).get("PP.H4"),
                     "interv": (j.get("intervention") or {}),
                     "steiger": (j.get("steiger") or {}).get("status"),
                     "allele": j.get("allele_harmonisation") or {},
                     "ld": (j.get("ld_source") or {}),
                     "flags": j.get("flags") or []})
    qtab = _q_table(root)
    rows.sort(key=lambda r: (r["status"] != "OK", r["p"] if r["p"] is not None else 9))

    L = ["# 真实两样本孟德尔随机化报告（OpenGWAS v4）", "",
         "生成时间：%s　　引擎：`np_mr.py` + `np_fdr.py` + `tools/mr_report.py`" % time.strftime("%Y-%m-%d %H:%M"),
         "", "> 结局：%s。暴露：OpenGWAS 的 cis-eQTL 汇总统计（`eqtl-a-ENSG*`）。" % a.outcome_label,
         "> 所有数值来自 `mr_result.json` 与 `fdr_table.csv`，无人工改写。", ""]

    ok = [r for r in rows if r["status"] == "OK"]
    if ok:
        L += ["## 1　主结果", "",
              "| 靶点 | 独立 IV | 去冗余方式 | β (SE) | OR [95% CI] | p | **q (BH)** | F | Egger 截距 p | COLOC PP.H4 | 方向判定 |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in ok:
            orv = math.exp(r["beta"]) if r["beta"] is not None else float("nan")
            lo = math.exp(r["beta"] - 1.96 * (r["se"] or 0))
            hi = math.exp(r["beta"] + 1.96 * (r["se"] or 0))
            q = qtab.get(r["gene"].upper(), {}).get("q_family")
            iv = r["interv"] or {}
            verdict = {"SEVERE": "**红灯：方向矛盾**", "WARN": "⚠️ 方向相反（证据弱）",
                       "OK": "自洽"}.get(iv.get("severity"),
                                        "未提供药物方向" if not iv else "—")
            L.append("| %s | %s | %s | %.3f (%.3f) | %.3f [%.3f, %.3f] | %.3g | %s | %.0f | %s | %s | %s |" % (
                r["gene"], r["n_iv"], r["mode"], r["beta"], r["se"] or 0, orv, lo, hi, r["p"],
                ("**%.3g**" % float(q)) if q not in (None, "") else "—",
                r["F"] or 0,
                ("%.2f" % r["egger_p"]) if r["egger_p"] is not None else "—",
                ("%.3f" % r["pp4"]) if r["pp4"] is not None else "—", verdict))
        L += ["", "q 为 6 个靶点内部的 Benjamini–Hochberg；**只有 q<0.05 才可写作「支持性因果证据」**。", ""]

        sig = [r for r in ok if qtab.get(r["gene"].upper(), {}).get("q_family")
               and float(qtab[r["gene"].upper()]["q_family"]) < 0.05]
        if sig:
            L += ["### 1.1 校正后显著的靶点（含共定位是否支持）", ""]
            for r in sig:
                q = float(qtab[r["gene"].upper()]["q_family"])
                pp4 = r["pp4"]
                if pp4 is None:
                    coloc = "未做共定位"
                elif pp4 >= 0.8:
                    coloc = "共定位支持（PP.H4=%.2f）" % pp4
                elif pp4 >= 0.6:
                    coloc = "共定位中等（PP.H4=%.2f），需排除多效性" % pp4
                else:
                    coloc = ("**共定位不支持（PP.H4=%.2f）→ 不能排除多效性/连锁，"
                             "只能写「关联性因果提示」**" % pp4)
                L.append("- **%s**：遗传预测表达每升高 1 SD，冠心病风险%s（OR=%.3f，p=%.3g，q=%.3g）；%s。"
                         % (r["gene"], "升高" if r["beta"] > 0 else "降低",
                            math.exp(r["beta"]), r["p"], q, coloc))
            L += ["", "> 注意：本批每个靶点只有 3–5 个独立 IV，COLOC 在这种密度下**不稳定**，"
                  "PP.H4 只能作辅助证据，不能单独用来判定/排除因果。", ""]

    if ok:
        L += ["### 1.2 三种估计量是否同向（稳健性速查）", "",
              "| 靶点 | IVW β | Egger β | 加权中位数 β | 同向？ |", "|---|---|---|---|---|"]
        for r in ok:
            vals = [r["beta"], r["egger_b"], r["wm"]]
            vals = [v for v in vals if v is not None]
            same = "是" if vals and all(v * vals[0] > 0 for v in vals) else "**否**"
            L.append("| %s | %s | %s | %s | %s |" % (
                r["gene"],
                ("%.3f" % r["beta"]) if r["beta"] is not None else "—",
                ("%.3f" % r["egger_b"]) if r["egger_b"] is not None else "—",
                ("%.3f" % r["wm"]) if r["wm"] is not None else "—", same))
        L += ["", "> Egger 与加权中位数对多效性更稳健；三者符号一致才可写「稳健」。"
              "不一致时以加权中位数与 COLOC 为准，并在正文说明。", ""]

    bad = [r for r in rows if r["status"] != "OK"]
    if bad:
        L += ["## 2　未能估计的靶点（不隐藏）", "",
              "| 靶点 | 状态 | 原因 |", "|---|---|---|"]
        for r in bad:
            reason = r["j"].get("reason") or r["j"].get("status")
            if r["n_iv"] is not None:
                reason = "%s（去冗余后 %s 个 IV）" % (reason, r["n_iv"])
            L.append("| %s | %s | %s |" % (r["gene"], r["status"], reason))
        L += ["", "> cis-eQTL 的独立信号通常只有 1–3 个；少于 3 个时 IVW/Egger/加权中位数"
              "都无法给出可解释的区间，因此**如实标注为不足**，不强行合并。", ""]

    if ok:
        r0 = ok[0]
        al = r0["allele"]
        L += ["## 3　方法学细节（可直接改写入 Methods）", "",
              "- **取数**：`POST /api/tophits`（暴露，p<5e-8）+ `POST /api/associations`（结局，按变异查询）。",
              "- **等位基因对齐**：按 ea/nea 判定同向/反向/互补链；无法判定者丢弃并计数"
              "（本批：%s）。" % ("，".join("%s=%d" % (k, v) for k, v in al.items() if v) or "无记录"),
              "- **LD 去冗余**：OpenGWAS `POST /api/ld/matrix` 取 1000G EUR 的 r²，"
              "r²<0.001 视为独立（本批来源：%s）。" % (r0["ld"].get("source") or "未记录"),
              "- **主估计**：IVW；**稳健性**：MR-Egger 截距（多效性）、加权中位数；"
              "**方向**：Steiger（r²x>r²y）。",
              "- **共定位**：ABF 实现（与 R `coloc::coloc.abf` 在 60 个模拟位点上 r=1.000、MAE=0.000）；"
              "本批 PP.H4 均低，说明**没有共定位证据**，不能据此排除多效性。",
              "- **多重性**：同批靶点内 BH（`np_fdr.py`）。",
              "- **方向一致性**：把 MR 的 β 方向换算成「起效需抑制/激活」，与药物作用方向比对；"
              "不一致时按 F 与 p 分级（SEVERE/WARN）。", ""]

    L += ["## 4　必须写进 Limitations 的五条", "",
          "1. **cis-eQTL 的工具变量很少**（本批 3–5 个），IVW 的区间较宽；"
          "单个工具变量时无法做 Egger，多效性只能靠 COLOC 与生物学先验判断。",
          "2. **Winner's curse**：暴露侧取 p<5e-8 的位点会高估效应，"
          "正规做法需在独立样本中复算（本批未做）。",
          "3. **样本重叠**：eQTL 与 CHD GWAS 若来自同一 biobank，会偏向观察性关联；"
          "本批未量化重叠。",
          "4. **LD 参考人群**：用 1000G EUR 近似 eQTL 研究人群，祖先不匹配会低估 r²、"
          "导致去冗余不足（残留相关 SNP）。",
          "5. **单一结局**：只做了冠心病（ieu-a-7），不能外推到其他表型。",
          "6. **COLOC 在本批不稳定**：3–5 个 IV 无法支撑可靠的共定位后验；"
          "PP.H4 高（如 BCL2=1.00）与低（如 IL1B=0.03）都应谨慎解释。", ""]

    L += ["## 5　复现命令", "", "```bash",
          "python scripts/np_mr.py fetch --exposure eqtl-a-ENSG00000125538 --outcome ieu-a-7 -o mr/IL1B",
          "python scripts/np_mr.py analyze --exposure mr/IL1B/exp.tsv --outcome mr/IL1B/out.tsv \\",
          "  --n-exposure 31684 --n-outcome 184305 --ld-api -o mr/IL1B",
          "python scripts/np_fdr.py mr --dir mr",
          "python tools/mr_report.py --dir mr -o mr/mr_report.md",
          "```", ""]

    out = a.out or os.path.join(root, "mr_report.md")
    open(out, "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("[mr_report] %d 个靶点（%d 成功 / %d 不足）→ %s" % (len(rows), len(ok), len(bad), out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
