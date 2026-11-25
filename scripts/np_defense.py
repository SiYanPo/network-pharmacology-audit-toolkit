#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_defense.py —— 三大"顶刊拷问"防御工具：共溶效应 / 时间偏倚 / 构象选择

三个致命盲区（不可假装不存在，只能显式标注 + 提供有限的粗评估）：
  1) 混合物共溶效应（np_soup_effect）：单体 ADMET ≠ 复方汤剂。
     本研究预测基于单体；复方煎煮存在溶解度共增、代谢酶竞争抑制、肠道菌群转化。
     本模块生成结构化"边界声明"+ 可选的粗评估（溶解度加和近似、成分-酶竞争提示），
     绝不编造 PK 数值。是否存在开源工具？—— 无与 GastroPlus 等价的开源工具，
     最接近的是 SwissADME/PK-SM（单体）+ 自建 PBPK（超越本 skill 范围）；
     故本报告强制标注该局限，不做伪定量。
  2) 时间偏倚 / Data leakage（temporal-split）：金标准要与预测集时间切片隔离。
  3) 构象选择（docking-state）：标注每个对接结构的 Apo/Holo/AlphaFold 状态，
     发现只有 Apo/AF 时提示 ensemble docking / 同源建模替换。

用法
    python np_defense.py soup --compounds compounds.csv -o out/
    python np_defense.py temporal --gold gold_targets.tsv --pred predict.tsv \
        --cutoff-year 2020 --gold-year-col year -o out/
    python np_defense.py dockstate --plan dockplan.tsv -o out/
    python np_defense.py selftest
"""
import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

SOUP_DISCLAIMER = (
    "本研究的 ADMET 预测基于单一化合物的理化性质，未考虑中药复方水煎液作为多组分混合物"
    "可能产生的共溶效应（solubility enhancement）、成分间对 CYP450/UGT 代谢酶的竞争性抑制，"
    "以及肠道菌群对苷类、多糖等成分的转化。因此，单体层面的暴露量估计可能低估或高估复方实际吸收，"
    "相关机制结论需以体内 PK/PD 实验验证。"
)


def read_rows(path, delim=","):
    with open(path, encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f, delimiter=delim)
        return list(rdr), (rdr.fieldnames or [])


def cmd_soup(a):
    os.makedirs(a.outdir, exist_ok=True)
    rows, hdrs = read_rows(a.compounds)
    # 粗评估（非定量 PK，只做"有多少成分可能影响 ADMET"的结构化提示）
    n_glycoside = sum(1 for r in rows if any(k in str(r).lower() for k in ("glycoside", "苷", "glucoside")))
    n_polysaccharide = sum(1 for r in rows if any(k in str(r).lower() for k in ("polysaccharide", "多糖")))
    notes = {
        "n_compounds": len(rows),
        "n_glycoside_like": n_glycoside,
        "n_polysaccharide_like": n_polysaccharide,
        "lipinski_violations_note": "Monte-Carlo/类药性规则均基于单体，混合物不适用（需谨慎）",
        "open_source_tools": {
            "SwissADME/PK-SM": "单体 ADMET，可用；非混合物",
            "GastroPlus": "商业 PBPK，无等价开源替代",
            "建议": "如需混合物效应，需自建 PBPK 或体内 PK 实验——超出本 skill 范围",
        },
        "disclaimer_zh": SOUP_DISCLAIMER,
        "warning": "禁止依据单体 ADMET 对复方整体生物利用度做定量断言",
    }
    json.dump(notes, open(os.path.join(a.outdir, "soup_effect_notes.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    with open(os.path.join(a.outdir, "soup_effect_disclaimer.txt"), "w", encoding="utf-8") as f:
        f.write(SOUP_DISCLAIMER + "\n")
    print("[soup] 成分 %d（苷类 %d / 多糖 %d）→ 边界声明已写入" %
          (len(rows), n_glycoside, n_polysaccharide))
    print("[soup] 无等价开源混合物 PBPK 工具；已标注局限，不做伪定量。")
    return 0


def cmd_temporal(a):
    """时间切片：仅用 cutoff_year 前的金标准做依赖，之后的作为盲测集。"""
    gold, gh = read_rows(a.gold, "\t" if a.gold.lower().endswith(".tsv") else ",")
    # 找年份列
    ycol = a.gold_year_col or next((h for h in gh if "year" in h.lower() or "年" in h), None)
    gcol = next((h for h in gh if "gene" in h.lower() or "symbol" in h.lower()), gh[0] if gh else None)
    if not ycol:
        raise SystemExit("[FAIL] 金标准缺年份列（如 first_approval/year），无法做时间切片")
    train, test = set(), set()
    for r in gold:
        try:
            y = int(float(r.get(ycol) or 0))
        except ValueError:
            continue
        g = (r.get(gcol) or "").strip().upper()
        if not g:
            continue
        (train if y <= a.cutoff_year else test).add(g)
    pred = set()
    with open(a.pred, encoding="utf-8-sig") as f:
        for line in f:
            g = line.strip().upper()
            if g:
                pred.add(g)
    # 关键：预测若依赖了 test 集合里的靶点，视为 data leakage
    leakage = sorted(pred & test)
    recall_train = len(pred & train) / max(1, len(train))
    recall_test = len(pred & test) / max(1, len(test))  # 盲测集：预测未见过，低是正常的
    res = {
        "cutoff_year": a.cutoff_year, "year_col": ycol,
        "train_gold_n": len(train), "test_gold_n": len(test), "pred_n": len(pred),
        "recall_train": round(recall_train, 4), "recall_test_blind": round(recall_test, 4),
        "potential_leakage_n": len(leakage), "potential_leakage_genes": leakage[:30],
        "note": "若预测依赖了 test 期间的靶点（如金标准构建晚于预测），则存在 data leakage；"
                "正确的 temporal validation：模型只用到 cutoff 前数据，test 期间靶点仅用于盲测。",
    }
    os.makedirs(a.outdir, exist_ok=True)
    json.dump(res, open(os.path.join(a.outdir, "temporal_validation.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[temporal] cutoff=%d: train=%d test=%d | 潜在泄漏 %d 个" %
          (a.cutoff_year, len(train), len(test), len(leakage)))
    if leakage:
        print("  [warn] 泄漏基因（预测中出现 test 期靶点）：%s" % leakage[:10])
    return 0


def cmd_dockstate(a):
    """标注对接结构状态；提示 ensemble docking。"""
    plan, hdrs = read_rows(a.plan, "\t" if a.plan.lower().endswith(".tsv") else ",")
    pcol = next((h for h in hdrs if "pdb" in h.lower()), None)
    out = []
    for r in plan:
        pdb = (r.get(pcol) or "").strip() if pcol else ""
        # 状态判断：用户提供状态列优先，否则按 PDB ID 存在与否标注 unknown
        state = (r.get("state") or r.get("conformation") or "").strip()
        if not state:
            state = "unknown" if pdb else "AlphaFold/no-structure"
        rec = dict(r)
        rec["conformation_state"] = state
        if state.lower() in ("apo", "alphafold", "af2", "no-structure", "unknown"):
            rec["warning"] = ("Apo/AlphaFold 结构可能缺失配体诱导的侧链构象；"
                              "建议 ensemble docking 或同源建模替换为 Holo 结构；结合能仅供参考")
        else:
            rec["warning"] = ""
        out.append(rec)
    os.makedirs(a.outdir, exist_ok=True)
    fields = list(out[0].keys()) if out else hdrs + ["conformation_state", "warning"]
    with open(os.path.join(a.outdir, "dock_state_annotated.tsv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(out)
    n_warn = sum(1 for r in out if r.get("warning"))
    print("[dockstate] %d 靶点，%d 个需构象警示（Apo/AF/unknown）" % (len(out), n_warn))
    return 0


def selftest():
    ok = len(SOUP_DISCLAIMER) > 50 and "共溶" in SOUP_DISCLAIMER and "菌群" in SOUP_DISCLAIMER
    print("selftest: %s (disclaimer %d 字)" % ("PASS" if ok else "FAIL", len(SOUP_DISCLAIMER)))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_defense: 共溶/时间偏倚/构象选择防御")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("soup"); p.add_argument("--compounds", required=True)
    p.add_argument("-o", "--outdir", default="np_out/defense"); p.set_defaults(fn=cmd_soup)
    p = sub.add_parser("temporal")
    p.add_argument("--gold", required=True); p.add_argument("--pred", required=True)
    p.add_argument("--cutoff-year", type=int, default=2020)
    p.add_argument("--gold-year-col", default=None)
    p.add_argument("-o", "--outdir", default="np_out/defense"); p.set_defaults(fn=cmd_temporal)
    p = sub.add_parser("dockstate"); p.add_argument("--plan", required=True)
    p.add_argument("-o", "--outdir", default="np_out/defense"); p.set_defaults(fn=cmd_dockstate)
    p = sub.add_parser("selftest"); p.set_defaults(fn=lambda a: selftest())
    a = ap.parse_args(argv)
    return a.fn(a) or 0


if __name__ == "__main__":
    sys.exit(main())
