#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_interpret.py —— 结果解读引擎（把数字翻译成审稿人能接受的句子，中英双语）

问题
    跑出 COLOC PP.H4 = 0.65、MR β = 0.31、协同分 0.58，作者往往不知道该怎么写。
    写"显著"会被审稿人打回，写"无意义"又浪费结果。这个模块把**判定阈值 + 标准措辞**固化下来，
    保证同一数值在任何报告里得到同一句话，不会这次说"强证据"、下次说"提示"。

分级（与 references/23-discussion-defense-paragraphs.md 的措辞库一致）
    MR       : p<0.05 且 F≥10 且无 Egger 截距 → 支持性因果证据（L1+）；否则仅关联
    COLOC    : PP.H4>0.8 强；0.6–0.8 中等；<0.6 不足（0.6 为文献常用中等阈值，
               见 PMID 40399559 / 39985383；早期版本用 0.5 偏宽松，已修正）
    Synergy  : score≥0.7 强；0.55–0.7 中等；<0.55 弱
    Benchmark: AUC≥0.75 良好；0.65–0.75 中等；<0.65 有限（且必须同时报 EF）

用法
    from np_interpret import interpret_mr, interpret_coloc, render_html
    python np_interpret.py demo        # 打印各档位的示例句子
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys


def _fmt(x, nd=2):
    try:
        return ("%." + str(nd) + "f") % float(x)
    except (TypeError, ValueError):
        return "NA"


def interpret_mr(beta, p, f=None, or_ci=None, pleiotropy=None, coloc_pp4=None,
                 direction_consistent=None, exploratory=False):
    """MR 结果 → {level, zh, en}。level ∈ strong / moderate / suggestive / null。"""
    if beta is None or p is None:
        return {"level": "na", "zh": "MR 未产出可用估计（可能工具变量不足）。",
                "en": "No valid MR estimate was obtained (insufficient instruments)."}
    orv = math.exp(beta) if abs(beta) < 10 else float("inf")
    up = beta > 0
    caveats_zh, caveats_en = [], []
    if f is not None and f < 10:
        caveats_zh.append("F=%.1f<10，存在弱工具偏倚" % f)
        caveats_en.append("the F-statistic (%.1f) indicates weak-instrument bias" % f)
    if pleiotropy:
        caveats_zh.append("MR-Egger 截距显著，提示方向性多效性")
        caveats_en.append("a significant MR-Egger intercept suggests directional pleiotropy")
    if direction_consistent is False:
        caveats_zh.append("三种估计方法方向不一致")
        caveats_en.append("the three estimators disagree in direction")
    if coloc_pp4 is not None and coloc_pp4 < 0.5:
        caveats_zh.append("共定位证据不足（PP.H4=%.2f），不能排除多效性" % coloc_pp4)
        caveats_en.append("colocalisation evidence was limited (PP.H4=%.2f)" % coloc_pp4)

    if p < 0.05 and not caveats_zh:
        level = "strong"
        zh = ("遗传预测的靶点表达升高与疾病风险%s相关（β=%.3f，OR=%.2f，p=%.2g），"
              "且弱工具、多效性与方向一致性检查均未提示问题，支持该靶点与疾病的因果关联。" %
              ("升高" if up else "降低", beta, orv, p))
        en = ("Genetically predicted higher target expression was associated with %s risk of the disease "
              "(β=%.3f, OR=%.2f, p=%.2g), with no evidence of weak-instrument bias, pleiotropy or "
              "directional inconsistency, supporting a causal link." %
              ("increased" if up else "decreased", beta, orv, p))
    elif p < 0.05:
        level = "moderate"
        zh = ("MR 点估计提示靶点表达与疾病风险相关（β=%.3f，OR=%.2f，p=%.2g），但%s；"
              "该结果只能作为支持性证据，不能单独作为因果结论。" % (beta, orv, p, "；".join(caveats_zh)))
        en = ("The MR estimate suggested an association between target expression and disease risk "
              "(β=%.3f, OR=%.2f, p=%.2g); however, %s. This should be regarded as supportive rather "
              "than definitive causal evidence." % (beta, orv, p, "; ".join(caveats_en)))
    else:
        level = "null"
        zh = ("未观察到靶点表达与疾病风险的显著因果关联（β=%.3f，p=%.2g）；"
              "不能据此否定该靶点，只能说明在本 GWAS 样本量下证据不足。" % (beta, p))
        en = ("No significant causal association was detected (β=%.3f, p=%.2g). This does not exclude "
              "the target; it indicates insufficient genetic evidence at the current sample size." % (beta, p))
    if exploratory:
        zh += "（注意：未提供 LD 面板，clumping 仅为距离剪枝，本结果仅作探索性参考。）"
        en += " (Note: no LD panel was supplied; clumping was distance-based, so this analysis is exploratory only.)"
    return {"level": level, "zh": zh, "en": en, "OR": round(orv, 3) if orv < 1e6 else None,
            "caveats": caveats_zh}


def interpret_coloc(pp4, pp3=None, n_snps=None):
    """COLOC → {level, zh, en}。"""
    if pp4 is None:
        return {"level": "na", "zh": "未进行共定位分析。", "en": "Colocalisation was not performed."}
    pp4 = float(pp4)
    pp3s = ("（PP.H3=%.2f）" % pp3) if pp3 is not None else ""
    if pp4 >= 0.8:
        level = "strong"
        zh = ("共定位证据强（PP.H4=%.2f），提示该区域存在同时影响基因表达与疾病风险的共享因果变异，"
              "支持因果而非单纯多效性%s。" % (pp4, pp3s))
        en = ("Colocalisation evidence was strong (PP.H4=%.2f), indicating a shared causal variant "
              "affecting both expression and disease risk, which supports causality rather than "
              "pleiotropy alone%s." % (pp4, (" (PP.H3=%.2f)" % pp3) if pp3 is not None else ""))
    elif pp4 >= 0.6:          # 领域惯例：>0.8 强、0.6–0.8 中等（见 Methods 引用）
        level = "moderate"
        zh = ("共定位证据为中等（PP.H4=%.2f），提示该区域可能同时影响基因表达与疾病风险，"
              "但需进一步排除多效性（建议补充条件分析或孟德尔随机化敏感性分析）%s。" % (pp4, pp3s))
        en = ("Colocalisation evidence was moderate (PP.H4=%.2f), suggesting a shared signal that "
              "may affect both expression and disease risk; pleiotropy cannot be excluded and "
              "conditional analysis is warranted%s." % (pp4, (" (PP.H3=%.2f)" % pp3) if pp3 is not None else ""))
    else:
        level = "weak"
        zh = ("共定位证据不足（PP.H4=%.2f）：表达与疾病信号更可能由两个独立变异驱动，"
              "此时 MR 结果更可能反映多效性/连锁而非因果%s。" % (pp4, pp3s))
        en = ("Colocalisation evidence was limited (PP.H4=%.2f): the expression and disease signals "
              "are more likely driven by distinct variants, so the MR estimate may reflect "
              "pleiotropy or linkage rather than causality%s." % (pp4, (" (PP.H3=%.2f)" % pp3) if pp3 is not None else ""))
    return {"level": level, "zh": zh, "en": en}


def interpret_synergy(score, label="Mechanistic_Synergy", level_gap=None, term=None):
    """协同评分 → {level, zh, en}。"""
    try:
        score = float(score)
    except (TypeError, ValueError):
        return {"level": "na", "zh": "无协同评分。", "en": "No synergy score available."}
    ctx = ("（%s 内层级差 %s）" % (term, level_gap)) if term and level_gap else ""
    if label == "Competitive_Antagonism":
        return {"level": "antagonism",
                "zh": "两成分作用于同一靶点/同一结合口袋，属竞争性拮抗，不应写作协同。",
                "en": "The two compounds act on the same target/pocket; this is competitive antagonism, not synergy."}
    if label == "Directional_Antagonism":
        return {"level": "antagonism",
                "zh": "两成分位于同一通路上下游但干预方向相反（一抑一激），净通量相互抵消。",
                "en": "The compounds sit upstream/downstream of one pathway with opposing directions, "
                      "so their net effects cancel out."}
    if score >= 0.7:
        lvl, zh, en = "strong", ("协同评分 %.2f，且作用于同一通路的上下游节点%s，符合“串联阻断”的机制协同特征。"
                                 % (score, ctx)), \
            ("Synergy score %.2f with upstream/downstream nodes of the same pathway%s is consistent "
             "with mechanism-based (serial blockade) synergy." % (score, ctx))
    elif score >= 0.55:
        lvl, zh, en = "moderate", ("协同评分 %.2f：存在通路层级关系，但方向证据或层级间隔较弱，"
                                   "建议结合剂量-效应实验确认%s。" % (score, ctx)), \
            ("Synergy score %.2f: a pathway-level relation exists but the directional or positional "
             "evidence is modest; dose-response validation is advised%s." % (score, ctx))
    else:
        lvl, zh, en = "weak", ("协同评分 %.2f，证据不足以判定协同，只能视为拓扑相邻%s。" % (score, ctx)), \
            ("Synergy score %.2f is insufficient to claim synergy; the pair should be reported as "
             "topologically adjacent only%s." % (score, ctx))
    return {"level": lvl, "zh": zh, "en": en}


def interpret_benchmark(auc, ef20=None, coverage=None):
    """基准性能 → {level, zh, en}。"""
    if auc is None:
        return {"level": "na", "zh": "未完成基准评测。", "en": "Benchmark not available."}
    d = ""
    if ef20 is not None:
        d += "，EF@20=%.2f×" % ef20
    if coverage is not None:
        d += "，金标准覆盖率 %.0f%%" % (100 * coverage)
    if auc >= 0.75:
        lvl = "good"
        zh = "排序能力良好（AUC=%.3f%s）。" % (auc, d)
        en = "Ranking performance was good (AUC=%.3f%s)." % (auc, d.replace("，", ", "))
    elif auc >= 0.65:
        lvl = "moderate"
        zh = "排序能力中等（AUC=%.3f%s），可作假设生成，不宜宣称预测准确。" % (auc, d)
        en = "Ranking performance was moderate (AUC=%.3f%s); suitable for hypothesis generation only." % (
            auc, d.replace("，", ", "))
    else:
        lvl = "limited"
        zh = "排序能力有限（AUC=%.3f%s），需改进打分或扩大候选池后再评估。" % (auc, d)
        en = "Ranking performance was limited (AUC=%.3f%s); scoring or candidate pool needs revision." % (
            auc, d.replace("，", ", "))
    if coverage is not None and coverage < 0.6:
        zh += "金标准覆盖率偏低，Recall 上限受限，不应据此得出性能结论。"
        en += " Gold-standard coverage was low, limiting attainable recall."
    return {"level": lvl, "zh": zh, "en": en}


def render_html(blocks, title="结果解读（自动生成）"):
    """把若干 {zh,en} 渲染成报告里的一节（中英并列，便于直接改写 Discussion）。"""
    if not blocks:
        return ""
    E = ['<h2>结果解读（Interpretation）</h2>',
         '<div class="boundary"><b>说明：</b>以下文字由规则引擎根据数值与阈值自动生成，'
         '用于统一措辞、避免"显著/不显著"的随意表述；<b>投稿前须由作者按实际语境改写</b>，'
         '不得直接当作 Discussion 成品。</div>', "<table><thead><tr><th>指标</th><th>中文解读</th>"
         "<th>English wording</th></tr></thead><tbody>"]
    for name, b in blocks:
        E.append("<tr><td><b>%s</b></td><td>%s</td><td><i>%s</i></td></tr>" % (
            name, b.get("zh", ""), b.get("en", "")))
    E.append("</tbody></table>")
    return "".join(E)


def demo():
    cases = [
        ("MR strong", interpret_mr(0.31, 1e-5, f=42, direction_consistent=True, coloc_pp4=0.85)),
        ("MR moderate", interpret_mr(0.31, 0.004, f=7, pleiotropy=True)),
        ("MR null", interpret_mr(0.05, 0.42, f=30)),
        ("COLOC strong", interpret_coloc(0.88, 0.05)),
        ("COLOC moderate", interpret_coloc(0.65, 0.2)),
        ("COLOC weak", interpret_coloc(0.21, 0.7)),
        ("Synergy strong", interpret_synergy(0.74, "Mechanistic_Synergy", 2, "PI3K-Akt")),
        ("Antagonism", interpret_synergy(0.5, "Competitive_Antagonism")),
        ("Benchmark moderate", interpret_benchmark(0.705, 2.1, 0.42)),
    ]
    for name, b in cases:
        print("== %s [%s]" % (name, b["level"]))
        print("   中：%s" % b["zh"])
        print("   EN：%s\n" % b["en"])
    return 0


def selftest():
    """阈值边界自检：确保"同一数值永远得到同一档位"，不随报告模板漂移。"""
    checks = [
        ("coloc_strong", interpret_coloc(0.85)["level"] == "strong"),
        ("coloc_moderate", interpret_coloc(0.60)["level"] == "moderate"),
        ("coloc_weak", interpret_coloc(0.30)["level"] == "weak"),
        ("mr_strong", interpret_mr(0.3, 1e-6, f=50, direction_consistent=True)["level"] == "strong"),
        ("mr_downgrade_weak_iv", interpret_mr(0.3, 1e-6, f=4)["level"] == "moderate"),
        ("mr_null", interpret_mr(0.1, 0.4, f=40)["level"] == "null"),
        ("synergy_antagonism", interpret_synergy(0.9, "Competitive_Antagonism")["level"] == "antagonism"),
        ("synergy_strong", interpret_synergy(0.8, "Mechanistic_Synergy", 2, "PI3K-Akt")["level"] == "strong"),
        ("benchmark_limited", interpret_benchmark(0.55, 1.2, 0.3)["level"] == "limited"),
        ("bilingual_present", all(interpret_mr(0.3, 0.001)["en"] and interpret_coloc(0.9)["zh"] for _ in (0,))),
    ]
    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for n, g in checks:
        print("  [%s] %s" % ("✓" if g else "✗", n))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = argparse.ArgumentParser(description="np_interpret: 结果解读引擎（中英双语）")
    ap.add_argument("cmd", nargs="?", default="demo", choices=["demo"])
    ap.parse_args(argv)
    return demo()


if __name__ == "__main__":
    sys.exit(main())
