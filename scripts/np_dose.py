#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_dose.py — 剂量可达性核算（饮片 → 成分摄入 → 游离暴露 → 与体外浓度比较）

为什么必须做
  "某成分在 50 µM 抑制了通路"这类体外结论，若体内游离浓度只有 0.05 µM，
  则整条机制链在暴露层面就不成立。这是中药网络药理学最容易被审稿人击穿的一环，
  也是最容易被作者忽略的一环（因为数据库只给 OB/DL，不给剂量）。

链条（每一环都可替换为实测值）
  饮片日剂量 g → ×含量% → ×提取率% → 日摄入 mg → ÷分服次数 → 单次剂量 mg
  → ×生物利用度 F → ÷分布容积 Vd·体重 → 血浆总 Cmax（µM）→ ×游离分数 fu → 游离 Cmax
  → 与体外有效浓度比值 → 分级

分级（体外浓度 / 体内游离 Cmax）
  ≤1      可达：体内游离浓度不低于体外有效浓度
  1–10    同量级：谨慎，需考虑组织蓄积与主动转运
  10–100  超出：体外浓度高出体内 1–2 个数量级，机制解释需降级
  >100    不可达：该体外实验的浓度在体内不可能达到，不能作为机制证据

用法
  python np_dose.py --daily 9 --content 0.042 --mw 302.24 --invitro 50 -o out/
  python np_dose.py --daily 9 --content 0.042 --mw 302.24 --target 0.5   # 反推所需剂量
"""
import argparse
import json
import os
import sys

import npcore as C


def grade(r):
    if r is None:
        return "未计算"
    if r <= 1:
        return "可达"
    if r <= 10:
        return "同量级"
    if r <= 100:
        return "超出"
    return "不可达"


def chain(daily_g, content_pct, extract_pct, f_pct, bw, vd, fu, doses, mw):
    intake_mg = daily_g * 1000.0 * (content_pct / 100.0) * (extract_pct / 100.0)
    per_dose = intake_mg / max(1, doses)
    cmax_uM = (f_pct / 100.0) * per_dose * 1000.0 / max(1e-9, mw * vd * bw)
    free_uM = cmax_uM * fu
    return {"intake_mg": intake_mg, "per_dose_mg": per_dose,
            "cmax_uM": cmax_uM, "free_uM": free_uM}


def main(argv=None):
    ap = C.Cli(description="剂量可达性核算：饮片→成分→游离暴露→与体外浓度比较")
    ap.add_argument("--herb", default="受试药")
    ap.add_argument("--compound", default="目标成分")
    ap.add_argument("--daily", type=float, required=True, help="饮片日剂量 g/day（药典用量）")
    ap.add_argument("--content", type=float, required=True, help="成分在饮片中的含量 %%（0.05 表示 0.05%%）")
    ap.add_argument("--mw", type=float, required=True, help="分子量 g/mol")
    ap.add_argument("--extract", type=float, default=50.0, help="提取率 %%（水煎默认 50）")
    ap.add_argument("--f", type=float, default=30.0, help="口服生物利用度 F %%（默认 30，属假设）")
    ap.add_argument("--bw", type=float, default=60.0, help="体重 kg")
    ap.add_argument("--vd", type=float, default=1.0, help="分布容积 L/kg")
    ap.add_argument("--fu", type=float, default=0.1, help="血浆游离分数（黄酮/皂苷常 <0.1）")
    ap.add_argument("--doses", type=int, default=2, help="每日分服次数")
    ap.add_argument("--invitro", type=float, default=None, help="体外有效/实验浓度 µM")
    ap.add_argument("--target", type=float, default=None, help="反推：目标游离浓度 µM")
    ap.add_argument("--range-content", default="0.5,2", help="含量敏感性倍数 下限,上限")
    ap.add_argument("--range-extract", default="30,70", help="提取率敏感性 %% 下限,上限")
    ap.add_argument("--range-fu", default="0.02,0.2", help="fu 敏感性 下限,上限")
    ap.add_argument("-o", "--outdir", default="np_out")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    man = C.Manifest("07_dose", a.outdir, vars(a))
    base = chain(a.daily, a.content, a.extract, a.f, a.bw, a.vd, a.fu, a.doses, a.mw)

    ratio = None
    if a.invitro:
        ratio = a.invitro / max(1e-12, base["free_uM"])

    # 敏感性：三参数取极值组合，给出比值区间
    lo_mul, hi_mul = [float(x) for x in a.range_content.split(",")]
    ex_lo, ex_hi = [float(x) for x in a.range_extract.split(",")]
    fu_lo, fu_hi = [float(x) for x in a.range_fu.split(",")]
    scenarios = {}
    for tag, (cm, em, fm) in {
        "最不利（低含量/低提取/低游离）": (lo_mul, ex_lo, fu_lo),
        "基准": (1.0, a.extract, a.fu),
        "最有利（高含量/高提取/高游离）": (hi_mul, ex_hi, fu_hi),
    }.items():
        r = chain(a.daily, a.content * cm, em, a.f, a.bw, a.vd, fm, a.doses, a.mw)
        scenarios[tag] = {"free_uM": r["free_uM"],
                          "ratio": (a.invitro / max(1e-12, r["free_uM"])) if a.invitro else None}

    rev = None
    if a.target:
        need_free = a.target
        need_cmax = need_free / max(1e-9, a.fu)
        # daily = cmax * mw * vd * bw * doses / (F * 1000) / (content% * extract%)
        daily_need = (need_cmax * a.mw * a.vd * a.bw * max(1, a.doses)
                      / ((a.f / 100.0) * 1000.0)
                      / ((a.content / 100.0) * (a.extract / 100.0)) / 1000.0)
        rev = {"target_free_uM": need_free, "required_daily_g": daily_need}

    C.ensure_dir(a.outdir)
    out = {"herb": a.herb, "compound": a.compound, "params": vars(a), "base": base,
           "grade": grade(ratio), "ratio": ratio, "scenarios": scenarios, "reverse": rev}
    p_json = C.write_text(os.path.join(a.outdir, "dose_result.json"),
                          json.dumps(out, ensure_ascii=False, indent=2))

    md = ["# 剂量可达性核算：%s · %s" % (a.herb, a.compound), "",
          "| 环节 | 取值 | 说明 |", "|---|---|---|",
          "| 饮片日剂量 | %.1f g | 须为**饮片**，非提取物 |" % a.daily,
          "| 成分含量 | %.4f%% | 无实测值时标注「估」 |" % a.content,
          "| 提取率 | %.0f%% | 水煎 30–70%%，默认 50 |" % a.extract,
          "| 日摄入量 | %.2f mg | = 日剂量×含量×提取率 |" % base["intake_mg"],
          "| 单次剂量 | %.2f mg | 分 %d 次服 |" % (base["per_dose_mg"], a.doses),
          "| 生物利用度 F | %.0f%% | **假设值**，非实测 |" % a.f,
          "| Vd / 体重 | %.2f L/kg / %.0f kg | 默认 1.0 / 60 |" % (a.vd, a.bw),
          "| 血浆总 Cmax | %.4f µM | 稳态未计，属上限估计 |" % base["cmax_uM"],
          "| 游离分数 fu | %.3f | 黄酮/皂苷常 <0.1 |" % a.fu,
          "| **游离 Cmax** | **%.4f µM** | 与体外对话的应是这个值 |" % base["free_uM"], ""]

    if ratio is not None:
        md += ["## 与体外浓度比较", "",
               "- 体外浓度：%.3g µM" % a.invitro,
               "- 比值（体外/体内游离）：**%.1f 倍** → **%s**" % (ratio, grade(ratio)), ""]
        if ratio > 10:
            md += ["> ⚠ 比值 >10：该体外浓度在体内不可能达到。若论文据此宣称机制，"
                   "必须在局限性中写明，或改用更低浓度复做实验。", ""]
        md += ["### 敏感性区间", "", "| 情景 | 游离 Cmax µM | 比值 | 分级 |", "|---|---|---|---|"]
        for tag, s in scenarios.items():
            md.append("| %s | %.4f | %s | %s |" % (
                tag, s["free_uM"],
                ("%.1f" % s["ratio"]) if s["ratio"] else "—",
                grade(s["ratio"]) if s["ratio"] else "—"))
        md.append("")
        lo = min(s["ratio"] for s in scenarios.values() if s["ratio"])
        hi = max(s["ratio"] for s in scenarios.values() if s["ratio"])
        md += ["> 比值区间 %.1f–%.1f：跨越数量级时，结论只能说「在有利假设下可达」，"
               "不得写成「体内可达有效浓度」。" % (lo, hi), ""]

    if rev:
        md += ["## 反推：要达到该游离浓度需要多少饮片", "",
               "- 目标游离浓度：%.3g µM" % rev["target_free_uM"],
               "- 所需饮片日剂量：**%.1f g/day**" % rev["required_daily_g"], ""]
        if rev["required_daily_g"] > 30:
            md += ["> ⚠ 所需剂量远超临床常用范围（>30 g/day）：该浓度在人用剂量下不可达，"
                   "体外结论不能外推到临床。", ""]

    md += ["## 未计算的部分（必须如实声明）", "",
           "- 未计**首过代谢的部位差异**、肠菌代谢活化、肝肠循环；",
           "- 未计**组织分布**（靶器官浓度可能高于或低于血浆）；",
           "- 未计**蛋白结合竞争**与多成分相互作用；",
           "- Cmax 为单次给药上限估计，非稳态谷/峰；",
           "- 全部参数为文献默认值时，本核算属**量级估计**，用于排除「差 2 个数量级」的情形，"
           "不能用于精确 PK 预测。", ""]

    p_md = C.write_text(os.path.join(a.outdir, "dose_report.md"), "\n".join(md) + "\n")
    man.output("result", p_json)
    man.output("report", p_md)
    man.count("free_cmax_uM", round(base["free_uM"], 6))
    if ratio is not None:
        man.count("ratio", round(ratio, 3))
        man.count("grade", grade(ratio))
        if ratio > 10:
            man.warn("体外浓度超出体内可达范围 %.0f 倍，机制结论须降级" % ratio)
    man.save()

    if not a.quiet:
        print("[剂量] %s·%s：日摄入 %.2f mg → 游离 Cmax %.4f µM（fu=%.3f）"
              % (a.herb, a.compound, base["intake_mg"], base["free_uM"], a.fu))
        if ratio is not None:
            print("[剂量] 体外/体内 = %.1f 倍 → %s" % (ratio, grade(ratio)))
        print("[剂量] 报告 → %s" % p_md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
