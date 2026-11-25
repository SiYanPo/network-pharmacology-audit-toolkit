#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_finalcheck.py —— 投稿前 Final Check 门禁（防"代码跑通 ≠ 能投稿"）

硬性检查（用户指定 + 审计补充 + 终极形态新增）：
    G1 时间切分盲测：存在 temporal_validation.json（np_defense temporal 产出）
    G2 Apo/Holo 状态标注：存在 dock_state_annotated.tsv
    G3 剂量权重：存在 dose_weighted_hubs.csv
    G8 Decoy 阴性对照：decoy_summary.json 且 n_red_flag=0（防对接假阳性）
    G13 共识打分UQ：consensus_summary.json 且 n_unstable=0（防打分不稳定）

审计补充项（不阻断，但写入报告）：
    G4 配伍禁忌质检（contraindication.json 存在且 verdict != ERROR）
    G5 外部金标准盲测（goldtest_result.json 存在）
    G6 图件自检通过（figure_qa.md 含 FAIL=0）
    G7 证据锚点校验（evidence 报告存在且 verdict=PASS）

门禁逻辑
    硬性项全通过 → 输出 Ready_for_Submission 水印（写入 HTML 报告页脚）
    任一未通过     → 只能输出 Draft 报告，并逐项列出缺什么、怎么补

诚实边界：水印只是"本 skill 自查清单通过"的声明，**不代表期刊接受**；
         水印文本明确写明这一点，避免误导。

用法
    python np_finalcheck.py --dir np_out/ -o np_out/
    python np_finalcheck.py --dir np_out/ --html report.html --stamp   # 给 HTML 加水印
    python np_finalcheck.py --selftest
"""
import argparse
import glob
import json
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

GATES = [
    ("G1", "时间切分盲测", "temporal_validation.json", "hard",
     "python np_defense.py temporal --gold ... --pred ... --cutoff-year 2020"),
    ("G2", "Apo/Holo 状态标注", "dock_state_annotated.tsv", "hard",
     "python np_defense.py dockstate --plan ..."),
    ("G3", "剂量权重", "dose_weighted_hubs.csv", "hard",
     "python np_doseweight.py --edges ... --dose ..."),
    ("G4", "配伍禁忌质检", "contraindication.json", "soft",
     "python np_contraindication.py --herbs-file ... -o out/"),
    ("G5", "外部金标准盲测", "goldtest_result.json", "soft",
     "python np_goldtest.py --gold-ctd ... --pred ..."),
    ("G6", "图件自检 FAIL=0", "figure_qa.md", "soft",
     "python np_figcheck.py --dir figs"),
    ("G7", "证据锚点校验", "evidence_report.json", "soft",
     "python np_evidence.py validate --text ... --out evidence_report.json"),
    # ---- 终极形态新增风险点（2026-11-22）----
    ("G8", "Decoy 阴性对照(无 Red Flag)", "decoy_summary.json", "hard",
     "python np_decoy.py score --core ... --decoy ...；n_red_flag 必须为 0"),
    ("G9", "单细胞投影", "singlecell_manifest.json", "soft",
     "python np_singlecell.py project --h5ad ... --targets ...（demo_mode 必须为 false 才能投稿）"),
    ("G10", "亚细胞定位分层", "localization_layers.csv", "soft",
     "python np_wetlab.py localization --targets ..."),
    ("G11", "引物特异性验证", "specificity_plan.json", "soft",
     "python np_wetlab.py specificity --primers ... --genome GRCm39.fa"),
    ("G12", "MD 桥接", "md_bridge_meta.json", "soft",
     "python np_wetlab.py md-bridge --receptor ..."),
    # ---- 第五轮终极形态新增（环境/屏障/UQ/组织）----
    ("G13", "共识打分UQ(无Unstable)", "consensus_summary.json", "hard",
     "python np_consensus.py ...；n_unstable 必须为 0，否则标 Unstable_Binding"),
    ("G14", "生理屏障过滤", "funnel_trace.json", "soft",
     "python np_funnel.py --cns-disease ...；检查 barrier_BBB / Physically Unreachable"),
    ("G15", "组织特异性校验", "hpa_summary.json", "soft",
     "python np_hpa.py check --tissue ...；核心靶点须在病变组织富集"),
    ("G16", "细胞通讯(LIANA)", "liana_status.json", "soft",
     "python np_liana.py run --h5ad ...；TOOL_MISSING 时提供脚本"),
    ("G17", "配体拓扑就绪", "ligand_top_result.json", "soft",
     "python np_ligand_top.py gen --ligand ...；MD 前须生成 GAFF .itp/.gro"),
    # ---- 第六轮因果与协同升维 ----
    ("G18", "MR 因果证据", "mr_result.json", "soft",
     "python np_mr.py analyze --exposure ... --outcome ...；靶点-疾病因果性加分项"),
    ("G19", "协同拮抗分析", "synergy_manifest.json", "soft",
     "python np_synergy.py --edges ...；多成分组合效应，揭示配伍科学内涵"),
    ("G20", "LINCS 重定位", "repurposing_candidates.json", "soft",
     "python np_lincs.py score --signature ... --l1000 ...；老药新用候选"),
    # ---- 第八轮：数据流一致性（防报告自相矛盾）----
    ("G21", "基因数据流一致性(交集⊆PPI)", "dataflow_check.json", "hard",
     "python np_pipeline.py --inputs inputs.json -o np_out --only ppi；"
     "或改用覆盖更全的 PPI 边表（STRING score≥0.4）后重跑 ppi 与 report"),
    ("G23", "证据聚合与校准(ECE<0.10)", "decision_rules.json", "soft",
     "python scripts/np_evidence_ledger.py fit --ledger np_out/09_evidence/evidence_ledger.csv "
     "--gold gold.txt -o np_out/09_evidence；无金标准时只报 UNCALIBRATED，不得写概率"),
    ("G22", "同源类型隔离(FuncRep≠直系同源)", "func_rep_warning.md", "hard",
     "python np_ortholog.py --input hubs.csv --col node -o 11_ortholog；"
     "FuncRep 必须带 (FuncRep) 后缀且不得当作直系同源引物"),
]



def _read_gene_col(path, preferred=("gene", "symbol", "node", "target")):
    """读 CSV/TSV 的基因列，返回大写符号列表（保持顺序去重）。"""
    import csv as _csv
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        rows = list(_csv.DictReader(f, delimiter=delim))
        hdrs = list(rows[0].keys()) if rows else []
        low = {h.lower().strip(): h for h in hdrs}
        col = next((low[c] for c in preferred if c in low), hdrs[0] if hdrs else None)
        if not col:
            return []
        return list(dict.fromkeys([(r.get(col) or "").strip().upper() for r in rows if (r.get(col) or "").strip()]))


def _ppi_nodes(d):
    """从 PPI 清洗边表取节点集合（node1/node2 或 source/target）。"""
    import csv as _csv
    cand = glob.glob(os.path.join(d, "**", "ppi_edges_clean.*"), recursive=True) + \
           glob.glob(os.path.join(d, "**", "string_edges.*"), recursive=True)
    if not cand:
        return set(), None
    p = cand[0]
    with open(p, encoding="utf-8-sig", errors="replace") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        rdr = _csv.DictReader(f, delimiter=delim)
        hdrs = {h.lower(): h for h in (rdr.fieldnames or [])}
        a = hdrs.get("node1") or hdrs.get("source") or hdrs.get("preferredname_a")
        b = hdrs.get("node2") or hdrs.get("target") or hdrs.get("preferredname_b")
        nodes = set()
        if a and b:
            for r in rdr:
                for k in (a, b):
                    v = (r.get(k) or "").strip().upper()
                    if v:
                        nodes.add(v)
    return nodes, p


def check_dataflow(d):
    """G21：交集基因必须真实存在于 PPI 网络节点中。

    为什么需要：报告 S4（交集表，按支持成分数排序）与 S5（PPI hub，按拓扑中心性排序）
    用的是**两套排名**，读者看到"交集第一名不在 hub 前 10"会怀疑数据造假。
    本函数把成员一致性（可判定）与排名顺序（本就不同，不可判定）分开：
      * 成员一致性 → coverage = |交集 ∩ PPI 节点| / |交集|；≥0.90 通过
      * 排名差异   → 记录 rank_map（交集基因在 hub 榜的名次），供报告并排展示，避免误读
    """
    inter_p = glob.glob(os.path.join(d, "**", "intersection_genes.txt"), recursive=True)
    genes = []
    if inter_p:
        genes = list(dict.fromkeys([l.strip().upper() for l in
                                    open(inter_p[0], encoding="utf-8-sig", errors="replace") if l.strip()]))
    else:
        ic = glob.glob(os.path.join(d, "**", "intersection.csv"), recursive=True)
        if ic:
            genes = _read_gene_col(ic[0])
    nodes, edge_file = _ppi_nodes(d)
    missing = [g for g in genes if g not in nodes]
    coverage = round(1 - len(missing) / len(genes), 4) if genes else None

    # 扩展子网（若已生成）：把孤立靶点的一阶邻居算进来，单独统计覆盖率
    ext_file = os.path.join(d, "05_ppi", "extended_ppi_edges.tsv")
    nodes_ext, missing_ext, coverage_ext = set(), list(missing), None
    if os.path.exists(ext_file):
        with open(ext_file, encoding="utf-8-sig", errors="replace") as f:
            first = f.readline()
            delim = "\t" if first.count("\t") > first.count(",") else ","
            f.seek(0)
            import csv as _csv
            for r in _csv.DictReader(f, delimiter=delim):
                for k in ("node1", "node2"):
                    v = (r.get(k) or "").strip().upper()
                    if v:
                        nodes_ext.add(v)
        nodes_ext |= nodes
        missing_ext = [g for g in genes if g not in nodes_ext]
        coverage_ext = round(1 - len(missing_ext) / len(genes), 4) if genes else None

    hubs_p = glob.glob(os.path.join(d, "**", "hubs.csv"), recursive=True)
    rank_map, hub_top = {}, []
    if hubs_p:
        import csv as _csv
        with open(hubs_p[0], encoding="utf-8-sig", errors="replace") as f:
            first = f.readline()
            delim = "\t" if first.count("\t") > first.count(",") else ","
            f.seek(0)
            for i, r in enumerate(_csv.DictReader(f, delimiter=delim), 1):
                g = (r.get("node") or r.get("gene") or "").strip().upper()
                if g:
                    rank_map[g] = {"hub_rank": int(r.get("rank") or i),
                                   "degree": r.get("degree", ""),
                                   "verdict": r.get("verdict", "")}
                    if i <= 10:
                        hub_top.append(g)
    rank_info = [{"gene": g, **rank_map.get(g, {"hub_rank": None})} for g in genes[:15]]
    best_cov = max([c for c in (coverage, coverage_ext) if c is not None] or [0])
    verdict = ("PASS" if best_cov >= 0.90 else ("FAIL" if best_cov < 0.75 else "WARN"))
    out = {"gate": "G21", "n_intersection": len(genes), "n_ppi_nodes": len(nodes),
           "n_missing": len(missing), "missing_genes": missing[:50], "coverage": coverage,
           "n_ppi_nodes_extended": len(nodes_ext) or None,
           "n_missing_extended": len(missing_ext), "missing_genes_extended": missing_ext[:50],
           "coverage_extended": coverage_ext,
           "verdict": verdict, "ppi_edge_file": edge_file,
           "hub_top10": hub_top, "intersection_rank_in_ppi": rank_info,
           "ranking_note": ("S4 交集表按支持成分数(n_compounds)排序，S5 hub 按拓扑中心性(degree/共识)排序，"
                            "两者排序键不同，名单不完全重叠属预期；本门禁校验的是成员一致性，不是排名一致性。"),
           "fix_hint": ("交集基因不在 PPI 中：①该基因在 STRING 无任何边（孤立节点）；"
                        "②PPI 边表被 score 阈值过滤过严。处理：降低 --min-score 或换覆盖更全的边表后重跑 ppi，"
                        "并在正文说明孤立节点数。")}
    try:
        json.dump(out, open(os.path.join(d, "dataflow_check.json"), "w", encoding="utf-8"),
                  ensure_ascii=False, indent=1)
    except Exception:
        pass
    return out


def _gate_g21(d):
    r = check_dataflow(d)
    ok = r["verdict"] == "PASS"
    reason = "%d/%d 个交集基因在 PPI 节点中（主网络 coverage %s）" % (
        r["n_intersection"] - r["n_missing"], r["n_intersection"], r["coverage"])
    if r.get("coverage_extended") is not None:
        reason += "；扩展子网 %s（%d 个节点）" % (r["coverage_extended"], r["n_ppi_nodes_extended"])
    if r["missing_genes"]:
        reason += "；缺失：" + ",".join(r["missing_genes"][:6]) + ("…" if len(r["missing_genes"]) > 6 else "")
    return {"id": "G21", "name": "基因数据流一致性(交集⊆PPI)", "level": "hard", "pass": ok,
            "reason": reason, "file": os.path.join(d, "dataflow_check.json")}



def _gate_g22(d):
    """G22：功能替代（FuncRep）与直系同源必须隔离，否则湿实验会用错引物。

    判定：
      * 无同源结果文件 → 通过（不阻断；只校验"一旦跑了就必须标对"）
      * 有 functional_analog / no_direct_ortholog 行 → ortholog 字段必须含 (FuncRep)
      * 出现 (FuncRep) 但 mapping_type≠functional_analog → 视为污染，失败
      * 有 FuncRep 却没有 func_rep_warning.md → 失败（警示文件缺失）
    """
    import csv as _csv
    csvs = glob.glob(os.path.join(d, "**", "*_orthologs.csv"), recursive=True)
    if not csvs:
        return {"id": "G22", "name": "同源类型隔离(FuncRep≠直系同源)", "level": "hard",
                "pass": True, "reason": "未运行同源转换（若运行则强制校验标记）"}
    bad, n_func = [], 0
    for p in csvs:
        with open(p, encoding="utf-8-sig", errors="replace") as f:
            for r in _csv.DictReader(f):
                orth = " ".join(str(r.get(k) or "") for k in
                                ("ortholog", "mouse", "rat", "zebrafish"))
                mt = (r.get("mapping_type") or "").strip().lower()
                st = (r.get("status") or "").strip().lower()
                has_mark = "(FuncRep)" in orth or "(funcrep)" in orth.lower()
                is_func = mt == "functional_analog" or st == "no_direct_ortholog"
                if is_func:
                    n_func += 1
                    if not has_mark:
                        bad.append("%s 标为功能替代但缺 (FuncRep) 后缀" % r.get("human", "?"))
                elif has_mark:
                    bad.append("%s 带 (FuncRep) 却声明为直系同源（类型污染）" % r.get("human", "?"))
    warn = glob.glob(os.path.join(d, "**", "func_rep_warning.md"), recursive=True)
    if n_func and not warn:
        bad.append("存在 %d 个功能替代但缺少 func_rep_warning.md 警示文件" % n_func)
    ok = not bad
    reason = ("%d 个功能替代均带 (FuncRep) 标记，且警示文件齐备" % n_func) if ok else         "；".join(bad[:4]) + "（严禁将 FuncRep 当作直系同源引物）"
    return {"id": "G22", "name": "同源类型隔离(FuncRep≠直系同源)", "level": "hard",
            "pass": ok, "reason": reason, "file": warn[0] if warn else ""}



def _gate_g23(d):
    """G23：是否存在可辩护的单一排序函数，且样本外校准可接受（ECE<0.10）。

    这是"你的最终排序函数是什么"的机器答案：没有账本 → 未通过；有账本无校准器 →
    概率不可报告；有校准器但 OOF ECE≥0.10 → 概率不可作决策依据。
    """
    dr = os.path.join(d, "09_evidence", "decision_rules.json")
    mp = os.path.join(d, "09_evidence", "calibration_model.json")
    if not os.path.exists(dr):
        return {"id": "G23", "name": "证据聚合与校准(ECE<0.10)", "level": "soft", "pass": False,
                "reason": "未生成证据账本/决策规则；报告只能并列多套排序，"
                          "审稿人会问「最终排序函数是什么」",
                "fix": "python scripts/np_evidence_ledger.py auto --dir <outdir>"}
    if not os.path.exists(mp):
        return {"id": "G23", "name": "证据聚合与校准(ECE<0.10)", "level": "soft", "pass": False,
                "reason": "有账本但无校准器：概率不可报告（UNCALIBRATED），只能给先验排序",
                "fix": ("python scripts/np_evidence_ledger.py fit --ledger "
                        "<outdir>/09_evidence/evidence_ledger.csv --gold gold.txt -o <outdir>/09_evidence")}
    try:
        m = json.load(open(mp, encoding="utf-8"))
    except Exception:
        return {"id": "G23", "name": "证据聚合与校准(ECE<0.10)", "level": "soft", "pass": False,
                "reason": "calibration_model.json 不可解析"}
    cv = m.get("cv") or {}
    e = cv.get("ece")
    if e is None:
        return {"id": "G23", "name": "证据聚合与校准(ECE<0.10)", "level": "soft", "pass": False,
                "reason": "校准器缺样本外（OOF）指标，无法判断校准质量"}
    ok = float(e) < 0.10
    return {"id": "G23", "name": "证据聚合与校准(ECE<0.10)", "level": "soft", "pass": ok,
            "reason": "样本外 ECE=%.4f（阈值 0.10）、Brier=%s、AUC=%s，n_oof=%s" % (
                float(e), cv.get("brier", "?"), cv.get("auc", "?"), cv.get("n_oof", "?")),
            "file": mp}


def find_file(d, name):
    """在目录树内找文件（容忍子目录）。"""
    hits = glob.glob(os.path.join(d, "**", name), recursive=True)
    return hits[0] if hits else None


def check_gate(gid, name, fname, level, d):
    if gid == "G21":                      # 计算型门禁：不看文件是否存在，而是现算数据流
        return _gate_g21(d)
    if gid == "G23":
        return _gate_g23(d)
    if gid == "G22":
        return _gate_g22(d)
    p = find_file(d, fname)
    if not p:
        return {"id": gid, "name": name, "level": level, "pass": False,
                "reason": "缺文件 %s" % fname}
    try:
        if fname.endswith(".json"):
            data = json.load(open(p, encoding="utf-8"))
            # 特殊判定
            if gid == "G4" and data.get("verdict") == "ERROR":
                return {"id": gid, "name": name, "level": level, "pass": False,
                        "reason": "配伍禁忌质检触发 ERROR（方剂含十八反/十九畏）"}
            if gid == "G7" and data.get("report", data).get("verdict") not in (None, "PASS"):
                return {"id": gid, "name": name, "level": level, "pass": False,
                        "reason": "证据锚点校验未通过"}
            if gid == "G8" and data.get("n_red_flag", 0) > 0:
                return {"id": gid, "name": name, "level": level, "pass": False,
                        "reason": "有 %d 个靶点触发 Red Flag: High False Positive Risk" % data["n_red_flag"]}
            if gid == "G9" and data.get("demo_mode", False):
                return {"id": gid, "name": name, "level": level, "pass": False,
                        "reason": "单细胞投影为 DEMO 合成数据，不可投稿；需真实 h5ad"}
            if gid in ("G11", "G12", "G17") and data.get("status") in ("TOOL_MISSING", "FAILED"):
                return {"id": gid, "name": name, "level": level, "pass": False,
                        "reason": "脚本已生成但工具未安装/生成失败（%s）" % data.get("status")}
            if gid == "G13" and data.get("n_unstable", 0) > 0:
                return {"id": gid, "name": name, "level": level, "pass": False,
                        "reason": "有 %d 个配体触发 Unstable_Binding（多打分方差过大）" % data["n_unstable"]}
            if gid == "G16" and data.get("status") == "TOOL_MISSING":
                return {"id": gid, "name": name, "level": level, "pass": False,
                        "reason": "LIANA 未装，已生成脚本但未产出真实通讯结果"}
            return {"id": gid, "name": name, "level": level, "pass": True, "file": p}
        if fname.endswith(".md") and gid == "G6":
            txt = open(p, encoding="utf-8", errors="ignore").read()
            m = re.search(r"FAIL\s*(\d+)", txt)
            fail = int(m.group(1)) if m else -1
            if fail > 0:
                return {"id": gid, "name": name, "level": level, "pass": False,
                        "reason": "图件自检有 %d 个 FAIL" % fail}
            return {"id": gid, "name": name, "level": level, "pass": True, "file": p}
        return {"id": gid, "name": name, "level": level, "pass": True, "file": p}
    except Exception as e:
        return {"id": gid, "name": name, "level": level, "pass": False,
                "reason": "解析异常：%s" % type(e).__name__}


def fix_command(gid, outdir):
    """为未通过项生成**可直接复制执行**的修复命令（智能诊断自愈）。

    返回 (诊断说明, 命令字符串)。命令使用本模块目录的绝对路径，避免 cwd 歧义。
    """
    S = os.path.join(HERE, "..", "scripts").replace("\\", "/")
    py = sys.executable
    cmds = {
        "G1": ("时间切分盲测缺失：需 CTD 金标准 + 预测靶点做 temporal split",
               '"%s" "%s/np_defense.py" temporal --gold gold/gold_pd_genes.tsv '
               '--pred disease_targets.tsv --cutoff-year 2020 -o %s' % (py, S, outdir)),
        "G2": ("缺少 Apo/Holo 构象状态标注：对接结构状态未声明",
               '"%s" "%s/np_defense.py" dockstate --plan docking_plan.tsv -o %s' % (py, S, outdir)),
        "G3": ("缺少剂量权重分析：无法证明成分暴露可达",
               '"%s" "%s/np_doseweight.py" --edges compound_targets.tsv --dose dose.csv -o %s' % (py, S, outdir)),
        "G4": ("配伍禁忌质检未通过：方剂触发禁忌（ERROR）",
               '"%s" "%s/np_contraindication.py" --herbs-file herbs.txt -o %s  # 查看具体禁忌对并调整处方'
               % (py, S, outdir)),
        "G5": ("外部金标准盲测缺失：结果不可证伪",
               '"%s" "%s/np_goldtest.py" --gold-ctd gold/gold_pd_genes.tsv '
               '--pred disease_targets.tsv --decoy gold/decoy.tsv -o %s' % (py, S, outdir)),
        "G6": ("图件自检存在 FAIL：需修复图片字号/重叠/越界",
               '"%s" "%s/np_figcheck.py" --dir figs --out figs/figure_qa.md' % (py, S)),
        "G7": ("Discussion 存在无证据锚点的断言句",
               '"%s" "%s/np_evidence.py" validate --text discussion.txt --datadir %s '
               '--out %s/evidence_report.json  # 按输出逐句补 [Data:]/[PMID:]'
               % (py, S, outdir, outdir)),
        "G8": ("Decoy 对照触发 Red Flag：对接假阳性风险高",
               '"%s" "%s/np_decoy.py" score --core docking_scores.csv --decoy decoy_scores.csv '
               '-o %s  # 检查 Z<-2 才可作为特异结合' % (py, S, outdir)),
        "G9": ("单细胞投影为 DEMO 或缺失：不可投稿",
               '"%s" "%s/np_singlecell.py" project --h5ad data.h5ad --targets core_targets.txt '
               '--cluster-key celltype -o figs/singlecell' % (py, S)),
        "G10": ("亚细胞定位缺失：无法排除核内错配",
                '"%s" "%s/np_wetlab.py" localization --targets core_with_uniprot.csv '
                '-o %s' % (py, S, outdir)),
        "G11": ("引物特异性未验证：可能非特异扩增",
                '"%s" "%s/np_wetlab.py" specificity --primers %s/primer/mouse_primers.csv '
                '--genome GRCm39.fa -o %s' % (py, S, outdir, outdir)),
        "G12": ("MD 桥接缺失：静态对接稳定性未验证",
                '"%s" "%s/np_wetlab.py" md-bridge --receptor dock/receptor.pdbqt -o %s/md' % (py, S, outdir)),
        "G13": ("共识打分不稳定（Unstable_Binding）：单一打分函数偏差大",
                '"%s" "%s/np_consensus.py" offline --scores three_scores.csv -o %s  '
                '# 或补 GNINA 重打分：np_consensus.py scan --gnina ...' % (py, S, outdir)),
        "G14": ("生理屏障过滤缺失：可能包含物理不可达成分",
                '"%s" "%s/np_funnel.py" --edges compound_targets.tsv --admet admet.csv '
                '--cns-disease -o %s/funnel  # 非CNS用 --target-tissue <组织>' % (py, S, outdir)),
        "G15": ("组织特异性校验缺失：无法排除细胞系假阳性",
                '"%s" "%s/np_hpa.py" fetch -o %s/hpa && "%s" "%s/np_hpa.py" check '
                '--targets core_targets.txt --tissue <疾病组织> --bulk %s/hpa/rna_tissue_consensus.tsv -o %s/hpa'
                % (py, S, outdir, py, S, outdir, outdir)),
        "G16": ("细胞通讯分析缺失：LIANA 脚本未跑",
                'bash %s/liana_status.json 中给出的 run_liana.py  # 先 pip install scanpy liana' % outdir),
        "G17": ("配体 GAFF 拓扑缺失：MD 无法启动",
                '"%s" "%s/np_ligand_top.py" gen --ligand dock/best_ligand.sdf -o %s/ligand_top' % (py, S, outdir)),
        "G18": ("缺少 MR 因果证据：靶点-疾病仅为相关性",
                '"%s" "%s/np_mr.py" analyze --exposure exp.tsv --outcome out.tsv -o %s/mr' % (py, S, outdir)),
        "G19": ("缺少协同拮抗分析：多成分组合效应未回答",
                '"%s" "%s/np_synergy.py" --edges compound_targets.tsv --pathways enrich.csv -o %s/synergy'
                % (py, S, outdir)),
        "G20": ("缺少 LINCS 重定位分析：临床转化价值未挖掘",
                '"%s" "%s/np_lincs.py" score --signature signature.json --l1000 profiles.csv -o %s/lincs'
                % (py, S, outdir)),
    }
    return cmds.get(gid, ("未知项，请查阅对应模块 --help", "# 无自动修复命令"))


def stamp_html(html_path, verdict, results):
    """在 </body> 前插入水印。"""
    if not os.path.exists(html_path):
        return None
    txt = open(html_path, encoding="utf-8").read()
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    if verdict == "READY":
        badge = ("<div style='position:fixed;bottom:12px;right:12px;padding:8px 14px;"
                 "background:#1b7a3d;color:#fff;font:13px sans-serif;border-radius:6px;"
                 "box-shadow:0 2px 8px #0004;z-index:9999'>"
                 "✓ Ready_for_Submission · np-pro 自查 %d/%d 通过 · %s</div>" %
                 (sum(1 for r in results if r["pass"]), len(results), ts))
    else:
        missing = [r["id"] for r in results if r["level"] == "hard" and not r["pass"]]
        badge = ("<div style='position:fixed;bottom:12px;right:12px;padding:8px 14px;"
                 "background:#b91c1c;color:#fff;font:13px sans-serif;border-radius:6px;"
                 "box-shadow:0 2px 8px #0004;z-index:9999'>"
                 "DRAFT · 硬性项未通过：%s · %s</div>" % (",".join(missing), ts))
    txt = txt.replace("</body>", badge + "</body>")
    open(html_path, "w", encoding="utf-8").write(txt)
    return html_path



# ==================== 故障自愈：把能自动化的门禁跑通 ====================
def _pipeline_inputs(d):
    mp = os.path.join(d, "pipeline_manifest.json")
    if not os.path.exists(mp):
        return {}
    try:
        return json.load(open(mp, encoding="utf-8")).get("inputs") or {}
    except Exception:
        return {}


def _run(cmd, timeout=600):
    import subprocess
    t0 = time.time()
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           timeout=timeout, text=True, encoding="utf-8", errors="replace")
        return p.returncode, (p.stdout or "")[-600:], round(time.time() - t0, 1)
    except Exception as e:
        return -1, "%s: %s" % (type(e).__name__, str(e)[:150]), round(time.time() - t0, 1)


def auto_fix(d, tissue=None, dry=False, do_pdf=True):
    """按依赖顺序尝试自动修复可自动化项；每一步都记录 rc 与命令，失败不掩盖。

    只做"输入已存在、结果可重算"的事：同源补全、扩展子网、HPA 富集重跑、报告与 PDF。
    绝不伪造：缺输入就 SKIP，并说明缺什么。
    """
    d = os.path.abspath(d)
    py = sys.executable
    S = os.path.dirname(os.path.abspath(__file__))
    inputs = _pipeline_inputs(d)
    steps = []

    def add(name, gate, cmd, ready, why=""):
        steps.append({"name": name, "gate": gate, "cmd": cmd, "ready": ready, "why": why,
                      "status": "PENDING", "rc": None, "out": "", "sec": 0})

    hubs = os.path.join(d, "05_ppi", "hubs.csv")
    edges = os.path.join(d, "05_ppi", "ppi_edges_clean.tsv")
    core = os.path.join(d, "04_intersect", "intersection_genes.txt")
    gmt = inputs.get("gmt") or ""
    add("同源转换（在线补全）", "G22",
        [py, os.path.join(S, "np_ortholog.py"), "--input", hubs, "--col", "node",
         "-o", os.path.join(d, "11_ortholog")], os.path.exists(hubs), "缺 05_ppi/hubs.csv")
    add("孤立靶点扩展子网", "G21",
        [py, os.path.join(S, "np_ppi.py"), "-e", edges, "--core", core, "--expand-isolated",
         "-o", os.path.join(d, "05_ppi_extended")],
        os.path.exists(edges) and os.path.exists(core), "缺 ppi_edges_clean.tsv 或交集基因表")
    enr = [py, os.path.join(S, "np_enrich.py"), "-g", core,
           "-o", os.path.join(d, "07_enrich"), "--background-mode", "hpa"]
    if gmt:
        enr += ["-m", gmt]
    if tissue:
        enr += ["--hpa-tissue", tissue]
    add("富集（HPA 组织背景）", "G7", enr, os.path.exists(core), "缺交集基因表")
    add("重新生成报告", "-",
        [py, os.path.join(S, "np_uireport.py"), "--dir", d, "--label", "自动修复后报告",
         "--verdict", "DRAFT"], True)

    for st in steps:
        if not st["ready"]:
            st["status"] = "SKIP"
            st["out"] = st["why"]
            continue
        if dry:
            st["status"] = "DRY"
            continue
        rc, out, sec = _run(st["cmd"])
        st.update({"status": "OK" if rc == 0 else "FAIL", "rc": rc, "out": out, "sec": sec})

    rep = os.path.join(d, "ui_report.html")
    if do_pdf and os.path.exists(rep) and not dry:
        rc, out, sec = _run([py, os.path.join(S, "np_pdf.py"), "--html", rep])
        steps.append({"name": "导出 PDF（无头浏览器）", "gate": "-", "cmd": [], "ready": True,
                      "status": "OK" if rc == 0 else "SKIP", "rc": rc, "out": out, "sec": sec,
                      "why": "未找到 Chrome/Edge 时会 SKIP"})

    print("\n[auto-fix] 执行结果")
    for st in steps:
        print("  %-6s %-24s rc=%s %5.1fs  %s" % (
            st["status"], st["name"], st["rc"] if st["rc"] is not None else "-", st["sec"],
            (st["out"] or "").splitlines()[-1][:90] if st["out"] else ""))
    json.dump({"steps": [{k: v for k, v in st.items() if k != "cmd"} for st in steps]},
              open(os.path.join(d, "autofix_log.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    return steps


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["selftest"]:        # 兼容统一写法：selftest 与 --selftest 等价
        argv = ["--selftest"]
    ap = argparse.ArgumentParser(description="np_finalcheck: 投稿门禁与水印")
    ap.add_argument("--dir", default="np_out")
    ap.add_argument("--html", default=None, help="待加水印的 HTML 报告")
    ap.add_argument("--stamp", action="store_true", help="写入水印（默认只检查）")
    ap.add_argument("--json", default=None, help="结果 JSON 输出")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--auto-fix", action="store_true", help="按依赖顺序自动跑可自动化的项，再复查门禁")
    ap.add_argument("--tissue", default=None, help="自动修复富集时使用的 HPA 组织")
    ap.add_argument("--dry-run", action="store_true", help="只打印将执行的命令")
    ap.add_argument("--no-pdf", action="store_true", help="自动修复时不导出 PDF")
    a = ap.parse_args(argv)
    if a.selftest:
        ok = all(len(g) == 5 for g in GATES) and any(g[3] == "hard" for g in GATES)
        print("selftest: %s (%d gates)" % ("PASS" if ok else "FAIL", len(GATES)))
        return 0 if ok else 1
    d = os.path.abspath(a.dir)
    if a.auto_fix:
        auto_fix(d, tissue=a.tissue, dry=a.dry_run, do_pdf=not a.no_pdf)
    results = [check_gate(g[0], g[1], g[2], g[3], d) for g in GATES]
    hard = [r for r in results if r["level"] == "hard"]
    soft = [r for r in results if r["level"] == "soft"]
    n_hard_pass = sum(1 for r in hard if r["pass"])
    verdict = "READY" if n_hard_pass == len(hard) else "DRAFT"

    md = ["# 投稿前 Final Check", "", "**判定：%s**" %
          ("Ready_for_Submission" if verdict == "READY" else "DRAFT（不可投稿）"), "",
          "检查时间：%s　目录：%s" % (time.strftime("%Y-%m-%d %H:%M:%S"), d), "",
          "| 项 | 名称 | 级别 | 状态 | 说明 |", "|---|---|---|---|---|"]
    for r in results:
        md.append("| %s | %s | %s | %s | %s |" %
                  (r["id"], r["name"], r["level"], "PASS" if r["pass"] else "FAIL",
                   r.get("reason", r.get("file", ""))[:80]))
    md += ["", "## 硬性项（必须全部通过才能出 Ready 水印）", ""]
    for r in hard:
        g = next(x for x in GATES if x[0] == r["id"])
        md.append("- %s %s：%s%s" % (r["id"], r["name"],
                  "✓ 通过" if r["pass"] else "✗ 缺失", "" if r["pass"] else " → 补：`%s`" % g[4]))
    # 智能诊断与自愈指令（DRAFT 时逐项给出可执行命令）
    failed = [r for r in results if not r["pass"]]
    if failed:
        md += ["", "## 自愈指令（逐项可复制执行）", "",
               "> 以下每条命令可直接复制运行；运行后重新执行本脚本直到判定 READY。", ""]
        for r in failed:
            diag, cmd = fix_command(r["id"], d)
            md += ["### %s %s — %s" % (r["id"], r["name"], r.get("reason", "")), "",
                   "**诊断**：%s" % diag, "", "```bash", cmd, "```", ""]
    md += ["", "## 水印声明", "",
           "> 本水印仅表示 network-pharmacology-pro 的**自查清单**通过，",
           "> 不代表期刊接受或同行评议通过。所有结论仍需实验验证与作者负责任解读。"]

    os.makedirs(d, exist_ok=True)
    rp = os.path.join(d, "finalcheck_report.md")
    open(rp, "w", encoding="utf-8").write("\n".join(md) + "\n")
    if a.json:
        json.dump({"verdict": verdict, "results": results},
                  open(a.json, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("[finalcheck] %s（硬性 %d/%d，软性 %d/%d）→ %s" %
          (verdict, n_hard_pass, len(hard), sum(1 for r in soft if r["pass"]), len(soft), rp))
    if failed:
        print("\n" + "=" * 64)
        print(" 智能诊断与自愈指令（%d 项未通过）" % len(failed))
        print("=" * 64)
        for r in failed:
            diag, cmd = fix_command(r["id"], d)
            print("\n  ✗ %s %s：%s" % (r["id"], r["name"], r.get("reason", "")))
            print("    诊断：%s" % diag)
            print("    修复：%s" % cmd)
    if a.stamp and a.html:
        p = stamp_html(a.html, verdict, results)
        print("[finalcheck] 水印已写入 %s" % p)
    return 0 if verdict == "READY" else 2


if __name__ == "__main__":
    sys.exit(main())
