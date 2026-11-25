#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_compounds.py — 成分筛选（三道闸：阈值门 / 给药途径 / 泛成分）

本阶段解决的四个真实审稿问题
  1) 阈值不交代、阈值与给药途径不匹配（腹腔注射/静脉给药却用口服生物利用度 OB 筛）
  2) 只报"筛出 N 个成分"，不报剔除理由 → 无法判断是否漏掉关键物质
  3) 槲皮素/山柰酚/木犀草素/β-谷甾醇等"泛成分"出现在几乎所有方剂里，被当成方剂特异性证据
  4) 【v1.2 新增】非中药来源的化合物（已上市药物、合成物、单体）没有 OB/DL，
     旧版会因缺 OB 一律剔除 —— 这不是「严谨」，而是「把无数据的分子当成不合格的分子」。
     现按行判定门控来源，且**每一行都写明它凭什么被纳入或剔除**。

门控模式（--gate）
  auto      逐行判定：有 OB+DL → OB/DL 阈值；无 OB/DL 但有 SMILES → 类药性
            （QED + Veber + PAINS/Brenk 结构警示）；两者都没有 → unassessed（不入下游，单独列出）
  obdl      仅用 OB/DL（TCMSP 体系，要求 ob/dl 列齐全）
  druglike  仅用类药性（要求 smiles 列）
  evidence  以「是否有实测活性/入血证据」为纳入依据（推荐用于已上市药物与单体），
            类药性只作标记、不作门槛。**这是唯一不会把阿托伐他汀、二甲双胍、地高辛
            这类已上市药物剔除的模式**——它们恰恰是 QED/Veber 眼中的"非类药"。
  none      不用任何门控（等同 evidence，但不计算标记），会写警告

输出
  compounds_kept.csv / compounds_exception.csv / compounds_rejected.csv / compounds_unassessed.csv
  每行含 gate_used 与 gate_detail —— 前者说明"用哪把尺子"，后者给出具体数值，可逐行复算。

用法
  python np_compounds.py -i compounds.csv -o out/
  python np_compounds.py -i compounds.csv -o out/ --route injection --force
  python np_compounds.py -i compounds.csv -o out/ --ob 30 --dl 0.18 --hl 4 --logic and
  python np_compounds.py -i compounds.csv -o out/ --gate auto --qed-min 0.30
  python np_compounds.py -i compounds.csv -o out/ --obdl-snapshot tcmsp_export.tsv
"""
import argparse
import os
import sys

import npcore as C

# 泛成分（pan-compound）：在大量方剂/药材中普遍出现，不能作为方剂特异性证据
UBIQUITOUS = {
    "quercetin", "kaempferol", "luteolin", "beta-sitosterol", "β-sitosterol", "stigmasterol",
    "isorhamnetin", "naringenin", "apigenin", "formononetin", "ferulic acid", "caffeic acid",
    "gallic acid", "vanillic acid", "p-coumaric acid", "rutin", "hyperoside", "quercitrin",
    "kaempferide", "chryseriol", "acacetin", "genistein", "daidzein", "hesperetin",
    "beta-carotene", "sitosterol", "campesterol", "oleic acid", "linoleic acid", "palmitic acid",
    "adenosine", "uracil", "nicotinic acid", "succinic acid", "salicylic acid",
}

ROUTE_NOTE = {
    "oral": "口服：OB/DL 阈值适用（TCMSP 体系的原设场景）。",
    "injection": "注射：OB 描述的是口服吸收，与注射暴露无关，**不应**用 OB 筛成分；"
                 "应以制剂中实测含量/入血成分为准，并在方法学中说明。",
    "topical": "外用/局部：OB/DL 不适用（皮肤渗透与局部浓度才是关键），应改用含量与渗透/滞留证据。",
    "unknown": "给药途径未声明：默认按口服处理，但必须在方法学中写明；否则该筛选不可复现。",
}

# Veber 口服生物利用度相关规则（Veber et al., J Med Chem 2002;45:2615）
VEBER = {"mw": 500.0, "clogp": 5.0, "hbd": 5, "hba": 10, "tpsa": 140.0, "rotb": 10}


def rdkit_props(smiles):
    """用 RDKit 计算类药性描述符与结构警示。RDKit 缺失或 SMILES 非法时返回 None。"""
    if not smiles:
        return None
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors
        from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    except Exception:
        return None
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None

    out = {
        "mw": round(Descriptors.MolWt(m), 2),
        "clogp": round(Crippen.MolLogP(m), 2),
        "hbd": int(rdMolDescriptors.CalcNumHBD(m)),
        "hba": int(rdMolDescriptors.CalcNumHBA(m)),
        "tpsa": round(rdMolDescriptors.CalcTPSA(m), 1),
        "rotb": int(rdMolDescriptors.CalcNumRotatableBonds(m)),
        "qed": round(QED.qed(m), 3),
        "heavy": m.GetNumHeavyAtoms(),
    }
    viol, why = 0, []
    for k in ("mw", "clogp", "hbd", "hba", "tpsa", "rotb"):
        if out[k] > VEBER[k]:
            viol += 1
            why.append("%s=%s>%s" % (k, out[k], VEBER[k]))
    out["veber_violations"] = viol
    out["veber_detail"] = "；".join(why)

    for name, params in (("pains", FilterCatalogParams.FilterCatalogs.PAINS),
                         ("brenk", FilterCatalogParams.FilterCatalogs.BRENK)):
        try:
            fp = FilterCatalogParams()
            fp.AddCatalog(params)
            out[name] = int(len(FilterCatalog(fp).GetMatches(m)))
        except Exception:
            out[name] = None
    return out


def load_snapshot(path):
    """读取用户自备的 OB/DL 快照（TCMSP 导出/自建表），用于离线复现 OB/DL 门控。

    要求列：compound, ob, dl（可选 hl, content, source）。
    返回 ({lower(name): row}, sha256)。哈希写入溯源，保证"用了哪一版快照"可核验。
    """
    if not path or not os.path.exists(path):
        return {}, None
    h, rows = C.read_table(path)
    c_n = C.resolve_col(h, "compound", required=True)
    c_ob, c_dl = C.resolve_col(h, "ob"), C.resolve_col(h, "dl")
    c_hl, c_ct = C.resolve_col(h, "hl"), C.resolve_col(h, "content")
    c_src = C.resolve_col(h, "source")
    d = {}
    for r in rows:
        nm = (r.get(c_n) or "").strip()
        if not nm:
            continue
        d[nm.lower()] = {
            "ob": C.to_float(r.get(c_ob)) if c_ob else None,
            "dl": C.to_float(r.get(c_dl)) if c_dl else None,
            "hl": C.to_float(r.get(c_hl)) if c_hl else None,
            "content": C.to_float(r.get(c_ct)) if c_ct else None,
            "source": (r.get(c_src) or "obdl-snapshot") if c_src else "obdl-snapshot",
        }
    return d, C.sha256_file(path)


def main(argv=None):
    ap = C.Cli(description="成分筛选：阈值门（OB/DL 或类药性）+ 给药途径一致性 + 泛成分标记")
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--outdir", default="np_out")
    ap.add_argument("--ob", type=float, default=30.0, help="口服生物利用度阈值%%，默认 30")
    ap.add_argument("--dl", type=float, default=0.18, help="类药性阈值，默认 0.18")
    ap.add_argument("--hl", type=float, default=None, help="半衰期阈值 h（可选，如 4）")
    ap.add_argument("--content", type=float, default=None, help="最低含量%%（可选，如 0.01）")
    ap.add_argument("--logic", choices=["and", "or"], default="and", help="阈值组合逻辑，默认 and")
    ap.add_argument("--route", choices=list(ROUTE_NOTE), default="unknown", help="给药途径")
    ap.add_argument("--force", action="store_true", help="明知途径不匹配仍按阈值筛（会写入警告）")
    ap.add_argument("--exceptions", default=None,
                    help="例外纳入名单（逗号分隔成分名或 CSV 文件），必须有实测/文献依据")
    ap.add_argument("--gate", choices=["auto", "obdl", "druglike", "evidence", "none"],
                    default="auto",
                    help="分子层面门控来源；auto=逐行判定（默认）")
    ap.add_argument("--qed-min", type=float, default=0.30,
                    help="druglike 门的 QED 下限（默认 0.30；QED>0.5 为通常意义的类药）")
    ap.add_argument("--max-veber", type=int, default=1,
                    help="druglike 门允许的 Veber 违例数上限（默认 1）")
    ap.add_argument("--obdl-snapshot", default=None,
                    help="OB/DL 快照表（TCMSP 导出等），用于离线填补 ob/dl 并留痕")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    man = C.Manifest("02_compounds", a.outdir, vars(a))
    man.input("compounds", a.input)

    heads, rows = C.read_table(a.input)
    c_name = C.resolve_col(heads, "compound", required=True)
    c_ob, c_dl = C.resolve_col(heads, "ob"), C.resolve_col(heads, "dl")
    c_hl, c_ct = C.resolve_col(heads, "hl"), C.resolve_col(heads, "content")
    c_herb, c_id = C.resolve_col(heads, "herb"), C.resolve_col(heads, "mol_id")
    c_smi = C.resolve_col(heads, "smiles")
    c_cat = C.resolve_col(heads, "category")
    c_cohort = C.resolve_col(heads, "cohort", extra=["approved", "drug_status", "group"])
    c_coh = C.resolve_col(heads, "cohort", extra=["approved", "drug_status", "group"])

    if a.gate == "obdl" and (not c_ob or not c_dl):
        print("[错误] --gate obdl 需要 ob/dl 列。现有列：%s" % "、".join(heads))
        print("       若该表是已上市药物/合成物（无 OB/DL），请用 --gate auto 或 --gate druglike。")
        return 2
    if a.gate == "druglike" and not c_smi:
        print("[错误] --gate druglike 需要 smiles 列。现有列：%s" % "、".join(heads))
        return 2

    snap, snap_sha = load_snapshot(a.obdl_snapshot)
    if snap:
        man.count("obdl_snapshot_rows", len(snap))
        man.input("obdl_snapshot", a.obdl_snapshot)

    # 途径一致性闸门
    warnings = []
    if a.route in ("injection", "topical") and not a.force:
        print("[阻断] 给药途径=%s，但 OB/DL 是口服场景参数，用其筛选会产生「验了没筛过的分子」式错误。" % a.route)
        print("       %s" % ROUTE_NOTE[a.route])
        print("       处理方式：")
        print("         A. 改用实测含量/入血成分表作为筛选依据（推荐）")
        print("         B. 明确声明这是「成分可成药性预筛」而非「暴露筛选」，并加 --force 继续")
        return 3
    if a.route in ("injection", "topical"):
        warnings.append("给药途径=%s 与 OB/DL 阈值不匹配，已按 --force 继续；方法学必须声明该局限。"
                        % a.route)

    exc = set()
    if a.exceptions:
        if os.path.exists(a.exceptions):
            _h, er = C.read_table(a.exceptions)
            cc = C.resolve_col(_h, "compound") or list(er[0].keys())[0]
            exc = {str(r.get(cc, "")).strip().lower() for r in er}
        else:
            exc = {x.strip().lower() for x in a.exceptions.split(",") if x.strip()}
    man.count("exceptions", len(exc))

    kept, exception, rejected, unassessed, seen = [], [], [], [], {}
    rejected_marketed = []
    herb_presence = {}
    gate_counts = {"obdl": 0, "druglike": 0, "evidence": 0, "none": 0}

    for r in rows:
        nm = (r.get(c_name) or "").strip()
        if not nm:
            continue
        key = nm.lower()
        herb = (r.get(c_herb) or "").strip() if c_herb else ""
        for h in C.split_multi(herb.replace("|", ";")):
            herb_presence.setdefault(key, set()).add(h)

        ob = C.to_float(r.get(c_ob)) if c_ob else None
        dl = C.to_float(r.get(c_dl)) if c_dl else None
        hl = C.to_float(r.get(c_hl)) if c_hl else None
        ct = C.to_float(r.get(c_ct)) if c_ct else None
        smi = (r.get(c_smi) or "").strip() if c_smi else ""
        prov = ""

        # 快照填补（只在原表缺值时，且来源写进 gate_detail）
        if key in snap:
            s = snap[key]
            if ob is None and s.get("ob") is not None:
                ob, prov = s["ob"], "ob:%s" % s.get("source", "snapshot")
            if dl is None and s.get("dl") is not None:
                dl = s["dl"]
                prov = (prov + "; " if prov else "") + "dl:%s" % s.get("source", "snapshot")
            if hl is None and s.get("hl") is not None:
                hl = s["hl"]
            if ct is None and s.get("content") is not None:
                ct = s["content"]

        # ---- 逐行判定门控来源 ----
        gate_used, gate_detail, props = "none", "", None
        if a.gate == "evidence":
            gate_used = "evidence"
            if smi:
                props = rdkit_props(smi)
            gate_detail = "按实测活性/入血证据纳入，类药性仅标记"
            if props:
                gate_detail += "（QED=%.2f，Veber 违例 %d；不作门槛）" % (
                    props["qed"], props["veber_violations"])
        elif a.gate == "none":
            gate_used = "none"
            gate_detail = "未做分子层面门控（用户指定）"
        elif a.gate == "obdl" or (a.gate == "auto" and ob is not None and dl is not None):
            gate_used = "obdl"
        elif a.gate == "druglike" or (a.gate == "auto" and smi):
            gate_used = "druglike"
            props = rdkit_props(smi)
            if props is None:
                gate_used = "unassessed"
                gate_detail = ("RDKit 不可用" if not smi else "SMILES 无法解析") + "，无法计算类药性"
        elif a.gate == "auto":
            gate_used = "unassessed"
            gate_detail = "既无 OB/DL 也无 SMILES，无法判定"

        passed = None
        fails = []
        if gate_used == "obdl":
            checks = []
            if ob is not None:
                checks.append(ob >= a.ob)
            if dl is not None:
                checks.append(dl >= a.dl)
            if a.hl is not None and hl is not None:
                checks.append(hl >= a.hl)
            if a.content is not None and ct is not None:
                checks.append(ct >= a.content)
            passed = (all(checks) if a.logic == "and" else any(checks)) if checks else False
            if ob is None or ob < a.ob:
                fails.append("OB%s" % ("=NA" if ob is None else "=%.2f<%g" % (ob, a.ob)))
            if dl is None or dl < a.dl:
                fails.append("DL%s" % ("=NA" if dl is None else "=%.2f<%g" % (dl, a.dl)))
            if a.hl is not None and (hl is None or hl < a.hl):
                fails.append("HL%s" % ("=NA" if hl is None else "=%.2f<%g" % (hl, a.hl)))
            if a.content is not None and (ct is None or ct < a.content):
                fails.append("含量%s" % ("=NA" if ct is None else "=%.4f%%<%g%%" % (ct, a.content)))
            gate_detail = "OB=%s, DL=%s%s" % (
                "%.2f" % ob if ob is not None else "NA",
                "%.2f" % dl if dl is not None else "NA",
                ("（来源：%s）" % prov) if prov else "")
        elif gate_used == "druglike":
            p = props
            ok_qed = p["qed"] >= a.qed_min
            ok_veber = p["veber_violations"] <= a.max_veber
            passed = bool(ok_qed and ok_veber)
            if not ok_qed:
                fails.append("QED=%.2f<%g" % (p["qed"], a.qed_min))
            if not ok_veber:
                fails.append("Veber 违例 %d>%d（%s）"
                             % (p["veber_violations"], a.max_veber, p["veber_detail"]))
            gate_detail = ("QED=%.2f; Veber 违例 %d/%d; MW=%.0f; cLogP=%.1f; TPSA=%.0f; RotB=%d"
                           % (p["qed"], p["veber_violations"], a.max_veber, p["mw"], p["clogp"],
                              p["tpsa"], p["rotb"]))
        elif gate_used in ("none", "evidence"):
            passed = True            # 证据模式不做门槛判定，纳入依据在 gate_detail 里
            gate_detail = gate_detail or "未门控"

        rec = {
            "compound": nm, "mol_id": r.get(c_id, "") if c_id else "",
            "herb": "|".join(sorted(herb_presence.get(key, []))),
            "OB": ob, "DL": dl, "HL": hl, "content": ct,
            "ubiquitous": "Y" if key in UBIQUITOUS else "",
            "n_herb": len(herb_presence.get(key, [])),
            "category": r.get(c_cat, "") if c_cat else "",
            "gate_used": gate_used, "gate_detail": gate_detail,
        }
        if props:
            rec.update({"qed": props["qed"], "mw": props["mw"], "clogp": props["clogp"],
                        "hbd": props["hbd"], "hba": props["hba"], "tpsa": props["tpsa"],
                        "rotb": props["rotb"], "veber_violations": props["veber_violations"],
                        "pains": props["pains"], "brenk": props["brenk"]})

        if key in exc:
            rec["reason"] = "例外纳入（用户提供实测/文献依据）"
            exception.append(rec)
            gate_counts[gate_used] = gate_counts.get(gate_used, 0) + 1
        elif gate_used == "unassessed":
            rec["reason"] = gate_detail
            unassessed.append(rec)
        elif passed:
            kept.append(rec)
            gate_counts[gate_used] = gate_counts.get(gate_used, 0) + 1
        else:
            rec["reason"] = "；".join(fails) if fails else "缺少数值"
            rejected.append(rec)
            # 类药性门控会把"已上市药物"判为不合格 —— 这不是数据问题，是判据错位。
            coh = (r.get(c_cohort) or "").strip().lower() if c_cohort else ""
            if gate_used in ("druglike", "auto") and coh in ("drug", "approved", "marketed"):
                rejected_marketed.append(nm)
        seen[key] = rec

    # 泛成分在 kept 中的占比 → 一个很实用的"特异性警示指标"
    n_ub = sum(1 for r in kept if r["ubiquitous"] == "Y")
    ub_ratio = n_ub / max(1, len(kept))
    n_pains = sum(1 for r in kept if r.get("pains"))
    mixed_gate = sum(1 for k in ("obdl", "druglike") if gate_counts.get(k))

    C.ensure_dir(a.outdir)
    cols = ["compound", "mol_id", "herb", "category", "OB", "DL", "HL", "content",
            "ubiquitous", "n_herb", "gate_used", "gate_detail",
            "qed", "mw", "clogp", "hbd", "hba", "tpsa", "rotb", "veber_violations",
            "pains", "brenk", "reason"]
    p_keep = C.write_table(os.path.join(a.outdir, "compounds_kept.csv"), cols, kept)
    p_exc = C.write_table(os.path.join(a.outdir, "compounds_exception.csv"), cols, exception)
    p_rej = C.write_table(os.path.join(a.outdir, "compounds_rejected.csv"), cols, rejected)
    p_una = C.write_table(os.path.join(a.outdir, "compounds_unassessed.csv"), cols, unassessed)
    p_list = C.write_text(os.path.join(a.outdir, "compounds_kept_names.txt"),
                          "\n".join(sorted(r["compound"] for r in kept)) + "\n")

    thr_desc = ("OB≥%g%%；DL≥%g；%s；逻辑=%s" % (
        a.ob, a.dl, ("HL≥%gh" % a.hl) if a.hl is not None else "未设 HL", a.logic.upper()))
    md = ["# 成分筛选报告", "",
          "| 项 | 值 |", "|---|---|",
          "| 输入记录 | %d |" % len(rows),
          "| 达标成分（kept） | %d |" % len(kept),
          "| 例外纳入（exception） | %d |" % len(exception),
          "| 剔除（rejected） | %d |" % len(rejected),
          "| 无法判定（unassessed） | %d |" % len(unassessed),
          "| 门控模式 | %s（OB/DL %d 行，类药性 %d 行） |"
          % (a.gate, gate_counts.get("obdl", 0), gate_counts.get("druglike", 0)),
          "| OB/DL 阈值 | %s |" % thr_desc,
          "| 给药途径 | %s |" % a.route,
          "| 泛成分占比 | %.0f%%（%d/%d） |" % (100 * ub_ratio, n_ub, len(kept)),
          "| 结构警示（PAINS） | %d 个达标成分命中 |" % n_pains,
          ""]

    # ---- 方法学表述：按实际用到的门控来源分别写，绝不混写 ----
    md += ["## 方法学可直接引用的表述", ""]
    if gate_counts.get("obdl"):
        md += ["> 成分筛选取 TCMSP 体系常用阈值（OB ≥ %g%% 且 DL ≥ %g），以 %s 逻辑组合；"
               "OB 与 DL 均为**预测值**（OBioavail 模型与 DrugBank Tanimoto 相似度），"
               "不等同于实测生物利用度。该阈值会系统性漏掉苷类、多糖与挥发性成分，"
               "故另设「例外纳入」清单（n=%d）。" % (a.ob, a.dl, a.logic.upper(), len(exception)), ""]
    if gate_counts.get("druglike"):
        md += ["> 对无 OB/DL 数据的化合物（已上市药物、合成单体等），改用可复现的类药性判据："
               "QED ≥ %g（Bickerton 2012, Nat Chem Biol 8:849）且 Veber 规则违例 ≤ %d 项"
               "（MW ≤ 500、cLogP ≤ 5、HBD ≤ 5、HBA ≤ 10、TPSA ≤ 140 Å²、可旋转键 ≤ 10；"
               "Veber 2002, J Med Chem 45:2615）。描述符由 RDKit 从规范 SMILES 计算，"
               "逐行数值见 `compounds_kept.csv`，可原样复算。" % (a.qed_min, a.max_veber), ""]
    if gate_counts.get("obdl") and gate_counts.get("druglike"):
        md += ["> 两类判据**不具可比性**：OB/DL 描述吸收与相似性，QED/Veber 描述理化性质空间。"
               "本表按行标注 `gate_used`，跨队列比较时不得混用同一阈值口径。", ""]
    if a.route != "oral":
        md += ["", "> 给药途径为 %s：%s" % (a.route, ROUTE_NOTE[a.route]), ""]

    if ub_ratio > 0.5:
        md += ["## ⚠ 特异性警示", "",
               "达标成分中 **%.0f%%** 属于在多数方剂中普遍出现的「泛成分」（槲皮素、山柰酚、"
               "木犀草素、β-谷甾醇等）。这类成分不能作为本方的特异性证据。" % (100 * ub_ratio),
               "建议：① 在结果中单列「泛成分」并说明；② 以非泛成分（novel/non-ubiquitous）"
               "重跑一次 PPI 与富集作为敏感性分析；③ 讨论中避免用泛成分解释方剂特有功效。", ""]
        warnings.append("泛成分占比 %.0f%%，特异性存疑" % (100 * ub_ratio))

    if rejected_marketed:
        md += ["## ⛔ 门控与「有效性」的错位（必须处理）", "",
               "本次以类药性（QED/Veber）剔除的成分中包含 **%d 个已上市药物**：%s。"
               % (len(rejected_marketed), "、".join(rejected_marketed)),
               "", "这些药物有确证的人体暴露与疗效，把它们判为「不合格」说明"
               "**类药性描述的是理化性质空间，不是药效**。若研究问题包含这些药物"
               "（阳性对照、机制对照、联合用药），必须改用 `--gate evidence`，"
               "或逐个列入 `--exceptions` 并在方法学说明依据。", "",
               "> 同理，OB/DL 阈值会系统性剔除苷类、皂苷、多糖——这与"
               "「它们没有活性」是两件事。这是网络药理学被诟病的主要方法论问题之一。", ""]
        warnings.append("%d 个已上市药物被类药性门控剔除（Atorvastatin/Metformin 等），"
                        "建议改用 --gate evidence" % len(rejected_marketed))

    if n_pains:
        md += ["## ⚠ 结构警示（PAINS）", "",
               "达标成分中有 **%d** 个命中 PAINS 警示结构。PAINS 是**实验假阳性**风险提示，"
               "不否定其药理活性，但意味着基于比色/荧光读数的体外活性证据需谨慎，"
               "并应优先采信非光学读数的实验（SPR、ITC、细胞表型）。" % n_pains, ""]
        warnings.append("%d 个达标成分命中 PAINS 结构警示" % n_pains)

    if unassessed:
        md += ["## 无法判定（未进入下游）", "",
               "| 成分 | 原因 |", "|---|---|"]
        for r in unassessed:
            md.append("| %s | %s |" % (r["compound"], r["reason"]))
        md += ["", "> 处置：补齐 OB/DL（用 `--obdl-snapshot` 导入 TCMSP 导出表）或补 SMILES。"
               "**不得**因为缺数据就默认纳入或默认剔除。", ""]
        warnings.append("%d 个成分因缺 OB/DL 与 SMILES 无法判定" % len(unassessed))

    if kept:
        top = sorted(kept, key=lambda r: (-(r["n_herb"] or 0),
                                          -(r["OB"] if r["OB"] is not None else 0)))[:15]
        md += ["## 达标成分（前 15）", "",
               "| 成分 | 类别 | 药材数 | 门控 | 依据 | 泛成分 |", "|---|---|---|---|---|---|"]
        for r in top:
            md.append("| %s | %s | %d | %s | %s | %s |" % (
                r["compound"], r.get("category", ""), r["n_herb"], r["gate_used"],
                r["gate_detail"], "是" if r["ubiquitous"] == "Y" else ""))
        md.append("")

    if exception:
        md += ["## 例外纳入（必须在方法学交代依据）", ""]
        for r in exception:
            md.append("- **%s**（门控=%s；%s）：%s" % (
                r["compound"], r["gate_used"], r["gate_detail"], r["reason"]))
        md.append("")

    if rejected:
        md += ["## 剔除理由分布", ""]
        cnt = {}
        for r in rejected:
            k = r["reason"].split("；")[0].split("<")[0].split("=")[0]
            cnt[k] = cnt.get(k, 0) + 1
        for k, v in sorted(cnt.items(), key=lambda x: -x[1]):
            md.append("- %s：%d 条" % (k, v))
        md.append("")

    p_md = C.write_text(os.path.join(a.outdir, "compounds_report.md"), "\n".join(md) + "\n")

    for k, v in [("kept", p_keep), ("exception", p_exc), ("rejected", p_rej),
                 ("unassessed", p_una), ("report", p_md)]:
        man.output(k, v)
    man.count("kept", len(kept))
    man.count("exception", len(exception))
    man.count("rejected", len(rejected))
    man.count("unassessed", len(unassessed))
    man.count("ubiquitous_ratio", round(ub_ratio, 3))
    man.count("pains", n_pains)
    man.count("gate_mixed", mixed_gate)
    man.count("gate_obdl", gate_counts.get("obdl", 0))
    man.count("gate_druglike", gate_counts.get("druglike", 0))
    man.count("gate_evidence", gate_counts.get("evidence", 0))
    man.count("rejected_marketed", len(rejected_marketed))
    for w in warnings:
        man.warn(w)
    man.save()

    if not a.quiet:
        print("[成分筛选] 达标 %d / 例外 %d / 剔除 %d / 无法判定 %d（门控：OB/DL %d、类药性 %d；泛成分 %.0f%%）"
              % (len(kept), len(exception), len(rejected), len(unassessed),
                 gate_counts.get("obdl", 0), gate_counts.get("druglike", 0), 100 * ub_ratio))
        for w in warnings:
            print("[警告]", w)
        print("[成分筛选] 报告 → %s" % p_md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
