#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_intervention.py —— 靶点-干预方向一致性校验（Target_Intervention_Mismatch）

要解决的"反向因果"陷阱
    MR 估的是"基因表达量 → 疾病风险"。药物做的是"阻断/激活蛋白"。两者方向必须对齐：
        MR β > 0（靶点升高 → 风险升高）→ 起效的干预必须是**抑制**该靶点
        MR β < 0（靶点升高 → 风险降低）→ 起效的干预必须是**激活**该靶点
    如果 MR 说"升高保护"，而你的成分网络预测"下调该靶点"，那么这条链路是**南辕北辙**的：
    成分越有效，病人越糟。不做这一步校验，Discussion 就可能写出与临床相反的结论。

判定分级（避免把弱证据当铁律）
    Grade A : p<0.05 且 F≥10 且 三法方向一致 且 Egger 截距不显著（+ COLOC PP.H4>0.5 更佳）
    Grade B : p<0.05 但存在异质性/多效性/弱工具之一 → 只作提示
    Grade C : p≥0.05（无因果证据）→ 不做方向判定，标 No_Causal_Evidence
    Mismatch 仅在 A/B 下抛出；C 下最多 INFO。SEVERE 只在 A 下。

用法
    # 单靶点：把 mr_result.json 与成分网络的动作表对齐
    python np_intervention.py check --mr mr/IL6/mr_result.json --action action.tsv -o out/
    # 批量：mr_dir/<GENE>/mr_result.json 或 mr_dir/<GENE>.json
    python np_intervention.py check --mr-dir mr/ --action action.tsv --edges ct_edges.csv -o out/
    # 从 MR 汇总表（gene,beta,p,F）
    python np_intervention.py check --mr-table mr_summary.csv --action action.tsv -o out/
    python np_intervention.py selftest

动作表格式（三选一，列名大小写不敏感）
    gene,action                action ∈ inhibit/activate/down/up/antagonist/agonist
    gene,compound,action       按成分逐条给（推荐，可追溯是哪个成分反向）
    gene,compound,direction     direction ∈ up/down（表达方向，down≈inhibit）
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
INHIBIT = {"inhibit", "inhibits", "inhibitor", "antagonist", "down", "decrease", "suppress", "抑制", "下调"}
ACTIVATE = {"activate", "activates", "agonist", "up", "increase", "induce", "激活", "上调"}
SEVERE = "Target_Intervention_Mismatch"


def norm_action(v):
    v = (v or "").strip().lower()
    if v in INHIBIT:
        return "inhibit"
    if v in ACTIVATE:
        return "activate"
    return None


def required_action(mr_beta):
    """MR β>0 → 需抑制；β<0 → 需激活。"""
    if mr_beta is None:
        return None
    return "inhibit" if mr_beta > 0 else "activate"


def grade_mr(rec):
    """MR 证据分级 A/B/C。rec 至少含 beta/p；可选 F/direction_consistent/pleiotropy/coloc。"""
    p = rec.get("p")
    if p is None or p >= 0.05:
        return "C", ["p>=0.05，无因果证据"]
    caveats = []
    f = rec.get("F")
    if f is not None and f < 10:
        caveats.append("弱工具 F<10")
    if rec.get("direction_consistent") is False:
        caveats.append("三法方向不一致")
    if rec.get("pleiotropy"):
        caveats.append("Egger 截距显著")
    if rec.get("heterogeneity"):
        caveats.append("异质性 Q 显著")
    return ("A" if not caveats else "B"), caveats


def load_action(path):
    """动作表 → {gene: {"n": k, "inhibit": [...], "activate": [...], "unknown": [...]}}"""
    out = {}
    with open(path, encoding="utf-8-sig") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        rdr = csv.DictReader(f, delimiter=delim)
        low = {(h or "").lower().strip(): h for h in (rdr.fieldnames or [])}
        gcol = next((low[c] for c in ("gene", "target", "symbol", "靶点") if c in low), None)
        ccol = next((low[c] for c in ("compound", "mol_id", "成分", "name") if c in low), None)
        acol = next((low[c] for c in ("action", "effect", "direction", "change", "调控") if c in low), None)
        if not gcol:
            raise SystemExit("[FAIL] 动作表缺 gene/target 列")
        if not acol:
            raise SystemExit("[FAIL] 动作表缺 action/direction 列")
        for r in rdr:
            g = (r.get(gcol) or "").strip().upper()
            if not g:
                continue
            a = norm_action(r.get(acol))
            c = (r.get(ccol) or "").strip() if ccol else ""
            slot = out.setdefault(g, {"inhibit": [], "activate": [], "unknown": []})
            slot[a if a else "unknown"].append(c or "?")
    return out


def load_mr(args):
    """返回 [{"gene":..., "beta":..., "p":..., "F":..., ...}]"""
    recs = []
    if args.mr_table:
        with open(args.mr_table, encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                g = (r.get("gene") or r.get("target") or "").strip().upper()
                if not g:
                    continue
                recs.append(_rec_from_row(g, r))
    paths = []
    if args.mr:
        paths.append(args.mr)
    if args.mr_dir and os.path.isdir(args.mr_dir):
        for name in sorted(os.listdir(args.mr_dir)):
            p = os.path.join(args.mr_dir, name)
            if os.path.isdir(p):
                p = os.path.join(p, "mr_result.json")
            if p.endswith(".json") and os.path.exists(p):
                gene = os.path.basename(os.path.dirname(p)) if os.path.basename(p) == "mr_result.json" \
                    else os.path.splitext(os.path.basename(p))[0]
                paths.append((p, gene.upper()))
    for item in paths:
        p, gene = (item, None) if isinstance(item, str) else item
        try:
            j = json.load(open(p, encoding="utf-8"))
        except Exception:
            continue
        if j.get("status") != "OK":
            continue
        ivw = j.get("IVW") or {}
        if not gene:
            gene = (j.get("gene") or j.get("exposure_gene") or
                    os.path.basename(os.path.dirname(os.path.abspath(p)))).upper()
        coloc = j.get("coloc") or {}
        recs.append({"gene": gene, "beta": ivw.get("beta"), "p": ivw.get("p"),
                     "F": j.get("F_statistic"),
                     "direction_consistent": j.get("direction_consistent"),
                     "pleiotropy": bool((j.get("MR_Egger") or {}).get("pleiotropy_warning")),
                     "heterogeneity": bool(ivw.get("heterogeneity")),
                     "coloc_pp4": coloc.get("PP.H4"),
                     "source": os.path.abspath(p)})
    return recs


def _rec_from_row(g, r):
    def f(k, cast=float):
        try:
            return cast(r.get(k))
        except (TypeError, ValueError):
            return None

    def b(k):
        v = (r.get(k) or "").strip().lower()
        return True if v in ("1", "true", "yes") else (False if v in ("0", "false", "no") else None)

    return {"gene": g, "beta": f("beta"), "p": f("p"), "F": f("F") or f("F_statistic"),
            "direction_consistent": b("direction_consistent"),
            "pleiotropy": b("pleiotropy") or False, "heterogeneity": b("heterogeneity") or False,
            "coloc_pp4": f("coloc_pp4"), "source": "table"}


def check(recs, actions):
    """逐靶点判定。返回 (rows, summary)。"""
    rows = []
    n_sev = n_warn = 0
    for rec in recs:
        g = rec["gene"]
        act = actions.get(g, {"inhibit": [], "activate": [], "unknown": []})
        grade, caveats = grade_mr(rec)
        need = required_action(rec.get("beta"))
        pred = None
        detail = ""
        if act["inhibit"] and act["activate"]:
            pred = "conflict"
        elif act["inhibit"]:
            pred = "inhibit"
        elif act["activate"]:
            pred = "activate"
        elif act["unknown"]:
            pred = "unknown"

        if grade == "C":
            verdict, severity, msg = "No_Causal_Evidence", "INFO", \
                "MR 无显著因果证据（p≥0.05），不做方向判定；该靶点不能作为因果卖点"
        elif need is None:
            verdict, severity, msg = "MR_Beta_Missing", "INFO", "MR 结果缺 β，无法判定"
        elif pred in (None, "unknown"):
            verdict, severity, msg = "Unknown_Action", "INFO", \
                "成分网络未给出该靶点的调控方向（需 --action 表），仅提示：因果要求%s该靶点" % \
                ("抑制" if need == "inhibit" else "激活")
        elif pred == "conflict":
            verdict, severity, msg = "Action_Conflict", "WARN", \
                "同一靶点上不同成分方向相反（%s vs %s），协同/拮抗需先厘清" % (
                    ",".join(act["inhibit"][:3]), ",".join(act["activate"][:3]))
            n_warn += 1
        elif pred == need:
            verdict, severity, msg = "Consistent", "OK", \
                "MR 要求%s该靶点，成分预测一致 → 方向自洽" % ("抑制" if need == "inhibit" else "激活")
        else:
            severity = "SEVERE" if grade == "A" else "WARN"
            verdict = SEVERE
            msg = ("遗传证据提示该靶点表达升高与疾病风险%s相关（β=%.3f, p=%.3g），"
                   "起效干预应为【%s】；而成分网络预测为【%s】→ 方向相反，临床结论会南辕北辙" % (
                       "升高" if (rec.get("beta") or 0) > 0 else "降低",
                       rec.get("beta") or 0, rec.get("p") or 1,
                       "抑制" if need == "inhibit" else "激活",
                       "抑制" if pred == "inhibit" else "激活"))
            if severity == "SEVERE":
                n_sev += 1
            else:
                n_warn += 1
        rows.append({"gene": g, "mr_beta": rec.get("beta"), "mr_p": rec.get("p"),
                     "mr_F": rec.get("F"), "mr_grade": grade,
                     "grade_caveats": ";".join(caveats),
                     "coloc_pp4": rec.get("coloc_pp4"),
                     "required_action": need or "", "predicted_action": pred or "",
                     "n_compounds_inhibit": len(act["inhibit"]),
                     "n_compounds_activate": len(act["activate"]),
                     "compounds": ";".join((act["inhibit"] + act["activate"])[:6]),
                     "verdict": verdict, "severity": severity, "message": msg})
    order = {"SEVERE": 0, "WARN": 1, "OK": 2, "INFO": 3}
    rows.sort(key=lambda r: (order.get(r["severity"], 9), -(abs(r["mr_beta"] or 0))))
    summary = {"n_targets": len(rows), "n_severe": n_sev, "n_warn": n_warn,
               "n_consistent": sum(1 for r in rows if r["verdict"] == "Consistent"),
               "n_no_evidence": sum(1 for r in rows if r["verdict"] == "No_Causal_Evidence"),
               "verdict": "FAIL" if n_sev else ("WARN" if n_warn else "PASS"),
               "rule": "MR β>0 → 需抑制靶点；MR β<0 → 需激活靶点；相反即 Mismatch",
               "note": "SEVERE 仅在 Grade A 证据下抛出；Grade B 降为 WARN，避免弱 MR 一票否决全部分析"}
    return rows, summary


def discussion_snippet(rows):
    sev = [r for r in rows if r["severity"] == "SEVERE"]
    if not sev:
        return ""
    items = "；".join("%s（MR β=%.2f, p=%.2g）" % (r["gene"], r["mr_beta"], r["mr_p"]) for r in sev[:5])
    return ("Directional consistency check. For %d target(s), genetic evidence indicated that "
            "increased expression was protective, whereas the compound–target network predicted "
            "inhibition: %s. Such target–intervention mismatch argues against a causal reading of "
            "these nodes and they were therefore downgraded to hypothesis-generating only." % (len(sev), items))


def selftest():
    checks = []
    # 1) β>0 → 需抑制
    checks.append(("required_inhibit", required_action(0.3) == "inhibit", required_action(0.3)))
    checks.append(("required_activate", required_action(-0.3) == "activate", required_action(-0.3)))
    # 2) Grade A
    g, cav = grade_mr({"beta": 0.3, "p": 1e-5, "F": 40, "direction_consistent": True,
                       "pleiotropy": False, "heterogeneity": False})
    checks.append(("grade_A", g == "A" and not cav, g))
    # 3) Grade C（无证据）
    checks.append(("grade_C", grade_mr({"beta": 0.3, "p": 0.2})[0] == "C", "C"))
    # 4) β>0 需抑制，但成分预测"激活" → SEVERE
    recs = [{"gene": "IL6", "beta": 0.42, "p": 1e-6, "F": 55, "direction_consistent": True,
             "pleiotropy": False, "heterogeneity": False}]
    acts_bad = load_action_from_dict({"IL6": {"inhibit": [], "activate": ["X"], "unknown": []}})
    rows, s = check(recs, acts_bad)
    checks.append(("mismatch_severe", rows[0]["verdict"] == SEVERE and s["verdict"] == "FAIL",
                   "%s/%s" % (rows[0]["verdict"], rows[0]["severity"])))
    # 5) β>0 需抑制，成分也抑制 → OK
    acts_ok = load_action_from_dict({"IL6": {"inhibit": ["Tanshinone IIA"], "activate": [], "unknown": []}})
    rows2, s2 = check(recs, acts_ok)
    checks.append(("consistent_ok", rows2[0]["verdict"] == "Consistent" and s2["verdict"] == "PASS",
                   rows2[0]["verdict"]))
    # 6) 弱证据(B)下相反 → WARN 而非 SEVERE
    recs_b = [dict(recs[0], F=5)]
    rows3, s3 = check(recs_b, acts_bad)
    checks.append(("gradeB_warn_only", rows3[0]["severity"] == "WARN" and s3["verdict"] == "WARN",
                   rows3[0]["severity"]))
    # 7) 无因果证据 → 不判定
    recs_c = [dict(recs[0], p=0.3)]
    rows4, _ = check(recs_c, acts_ok)
    checks.append(("no_evidence_no_verdict", rows4[0]["verdict"] == "No_Causal_Evidence",
                   rows4[0]["verdict"]))
    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for name, good, d in checks:
        print("  [%s] %-26s %s" % ("✓" if good else "✗", name, d))
    return 0 if ok else 1


def load_action_from_dict(d):
    return d


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_intervention: 靶点-干预方向一致性校验")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("check")
    p.add_argument("--mr", default=None, help="单个 mr_result.json")
    p.add_argument("--mr-dir", default=None, help="目录：<GENE>/mr_result.json 或 <GENE>.json")
    p.add_argument("--mr-table", default=None, help="CSV：gene,beta,p[,F]")
    p.add_argument("--action", required=True, help="成分→靶点动作表（gene,action 或 gene,compound,action）")
    p.add_argument("--edges", default=None, help="成分-靶点边表（仅用于补全成分名，可省略）")
    p.add_argument("-o", "--outdir", default="np_out/intervention")
    p.set_defaults(fn="check")
    t = sub.add_parser("selftest")
    t.set_defaults(fn="selftest")
    a = ap.parse_args(argv)
    if a.fn == "selftest":
        return selftest()
    recs = load_mr(a)
    if not recs:
        raise SystemExit("[FAIL] 未读到任何 MR 结果（--mr / --mr-dir / --mr-table）")
    acts = load_action(a.action)
    rows, summary = check(recs, acts)
    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, "intervention_check.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    snippet = discussion_snippet(rows)
    summary["discussion_snippet"] = snippet
    summary["n_checked"] = len(rows)
    json.dump(summary, open(os.path.join(a.outdir, "intervention_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[intervention] %d 靶点：SEVERE %d / WARN %d / Consistent %d / 无证据 %d → %s" % (
        len(rows), summary["n_severe"], summary["n_warn"], summary["n_consistent"],
        summary["n_no_evidence"], summary["verdict"]))
    for r in rows[:8]:
        if r["severity"] in ("SEVERE", "WARN"):
            print("  [%s] %s：%s" % (r["severity"], r["gene"], r["message"]))
    if snippet:
        print("\n[可直接写进 Discussion]\n%s" % snippet)
    return 3 if summary["n_severe"] else 0


if __name__ == "__main__":
    sys.exit(main())
