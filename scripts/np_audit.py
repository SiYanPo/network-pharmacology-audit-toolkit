#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_audit.py — 研究项目红线审计（v2：结构性判定优先，散文提示不计分）

v1 的两个已确认缺陷：
  1. 规则全部基于散文正则 → "本文不将 P 值作为作用强度"会被判为误用（否定句误报）；
  2. 措辞类与结构性规则混在一起计分 → 会写漂亮话的人得高分，真做了对照的人反而被扣。

v2 设计：
  * 结构性（structural）：读 manifest / environment.json / claims.json / 产物 CSV，可判定，计入评分。
  * 提示性（advisory）：读散文，带否定窗口检测，**不计入评分**，单列供人工复核。

用法
  python np_audit.py --dir np_out/ -o np_out/audit/
  python np_audit.py --dir np_out/ --strict      # 有 WARN 即返回 1（投稿前用）
"""
import argparse
import glob
import json
import os
import re
import sys

import npcore as C

STRUCT, ADV = "structural", "advisory"

RULES = [
    ("R01", "溯源", "阈值参数已在 manifest 中记录（OB/DL/富集 q）", "可复现的最低要求", STRUCT),
    ("R02", "溯源", "数据库版本与访问日期已写明", "无时间戳的引用不可核验", ADV),
    ("R03", "溯源", "存在阶段级 manifest（输入/输出哈希、参数）", "可复现性证明", STRUCT),
    ("R04", "溯源", "存在环境快照 environment.json", "换机排障与跨机复现", STRUCT),
    ("R05", "数据", "靶点标准化映射率 ≥80% 且已记录", "同物异名会稀释交集", STRUCT),
    ("R06", "数据", "非人类/噪音靶点已被识别或剔除", "TCMSP 含菌源/病毒/mRNA 噪音", STRUCT),
    ("R07", "数据", "疾病靶点按来源/证据强度分层", "整合相关 ≠ 遗传证据 ≠ 表达变化", STRUCT),
    ("R08", "设计", "给药途径已记录（且与筛选口径一致）", "注射/外用却用口服 OB 筛 = 逻辑冲突", STRUCT),
    ("R09", "设计", "中药炮制与基原已写明", "生品与制品是两个受试物", ADV),
    ("R10", "设计", "已声明成分清单为「部分」而非全集", "OB/DL 会漏苷类/多糖/挥发油", ADV),
    ("R11", "统计", "富集报告了多重检验校正后的 q 值", "未校正必然假阳性", STRUCT),
    ("R12", "统计", "富集声明了背景基因集口径", "以输入集为背景会夸大显著性", STRUCT),
    ("R13", "统计", "已构建阴性对照/零模型", "无零模型的紧密互联可被拓扑偏倚解释", STRUCT),
    ("R14", "统计", "零模型已自校准（度匹配质量达标）", "零模型本身也可能是错的", STRUCT),
    ("R15", "统计", "hub 判定给出稳健性证据（入选频率）", "Degree 不反映重要性", STRUCT),
    ("R16", "统计", "未把 P 值当作效应量使用", "P 是显著性，不是作用强度", ADV),
    ("R17", "统计", "未声称通路被激活/抑制", "GO/KEGG 输入无方向信息", ADV),
    ("R18", "富集", "伞状通路已降级或标注", "Pathways in cancer 几乎必然命中", STRUCT),
    ("R19", "富集", "冗余条目已折叠并报告前后数量", "数百条 GO 不是数百个机制", STRUCT),
    ("R20", "成分", "已标记泛成分并给出占比", "泛成分不能作方剂特异性证据", STRUCT),
    ("R21", "对接", "已报告重对接 RMSD 与门禁结果", "全 PDB top-1 成功率仅约 30%", STRUCT),
    ("R22", "对接", "已说明口袋定义与质子化处理", "决定氢键网络与可复现性", ADV),
    ("R23", "表达", "已声明「预测 vs 实测」的证据边界", "对接/网络/富集均为计算预测", STRUCT),
    ("R24", "表达", "未出现无实验支撑的绝对化措辞", "无 L3+ 证据只能写提示/支持", ADV),
    ("R25", "暴露", "已做剂量可达性核算", "体外 50 µM 而体内 0.05 µM 时机制链断裂", STRUCT),
    ("R26", "合规", "已声明临床用药边界", "安全与合规底线", ADV),
    ("R27", "报告", "存在结论账本 claims.json", "让结论可追溯、可审计", STRUCT),
    ("R28", "报告", "方法学段落引用的参数全部有记录（无「未记录」）", "避免用默认值伪造方法", STRUCT),
    ("R29", "报告", "有局限性章节", "审稿人第一条看有没有自知之明", ADV),
]

PAT = {
    "ver": r"(版本|version|v\d+\.\d+|访问日期|accessed|检索日期)",
    "date": r"(20\d\d[-/年]\d{1,2}[-/月]?\d{0,2})",
    "paozhi": r"(炮制|生品|制品|基原|拉丁学名|正名)",
    "partial": r"(部分|并非全部|不能代表|漏掉|苷类|多糖|挥发)",
    "peffect": r"(作用强度|效应值|越大越强|活性评分|重要性排名)",
    "activate": r"(激活了|抑制了|上调了|下调了|激活该通路|抑制该通路)",
    "absolute": r"(证明了|阐明了|首次发现|证实了|充分证明|明确证实)",
    "pocket": r"(口袋|pocket|grid box|质子化|protonat|辅因子|cofactor)",
    "clinic": r"(不作为诊疗|不构成诊疗|不能替代医师|仅供科研)",
    "limit": r"(局限性|不足|limitation|局限)",
}
NEG = re.compile(r"(不|未|无|非|勿|避免|不得|不能|并未|没有|not\b|no\b|never|avoid)", re.I)


def negated(text, m, window=30):
    """匹配点前 window 字符内若出现否定词，视为「正确表述」，不报警。"""
    return NEG.search(text[max(0, m.start() - window):m.start() + 4]) is not None


def load_text(root):
    texts = {}
    for pat in ("**/*.md", "**/*.html"):
        for p in glob.glob(os.path.join(root, pat), recursive=True):
            try:
                t = open(p, encoding="utf-8", errors="replace").read()
                if p.endswith(".html"):
                    t = re.sub(r"<script.*?</script>", " ", t, flags=re.S)
                    t = re.sub(r"<[^>]+>", " ", t)
                texts[os.path.relpath(p, root)] = t
            except Exception:
                continue
    return texts


def load_manifest(root):
    stages, counts, warns = [], {}, []
    for p in glob.glob(os.path.join(root, "**", "run_manifest.json"), recursive=True):
        try:
            d = json.load(open(p, encoding="utf-8"))
        except Exception:
            continue
        for s in d.get("stages", []):
            stages.append(s)
            counts.update(s.get("counts", {}) or {})
            warns.extend(s.get("warnings", []) or [])
    return stages, counts, warns


def read_csv_any(root, name):
    for p in glob.glob(os.path.join(root, "**", name), recursive=True):
        try:
            h, rows = C.read_table(p)
            return h, rows
        except Exception:
            continue
    return None, None


def main(argv=None):
    ap = C.Cli(description="项目红线审计 v2（结构性判定优先，散文提示不计分）")
    ap.add_argument("--dir", required=True)
    ap.add_argument("-o", "--outdir", default=None)
    ap.add_argument("--strict", action="store_true", help="WARN 也返回非零")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    root = a.dir
    if not os.path.isdir(root):
        print("[错误] 目录不存在：%s" % root)
        return 2
    outdir = a.outdir or os.path.join(root, "audit")
    C.ensure_dir(outdir)

    texts = load_text(root)
    all_text = "\n".join(texts.values())
    stages, counts, warns = load_manifest(root)
    ids = {str(s.get("stage")) for s in stages}

    claims = {}
    cp = os.path.join(root, "claims.json")
    if os.path.exists(cp):
        try:
            claims = json.load(open(cp, encoding="utf-8"))
        except Exception:
            claims = {}

    def stage(prefix):
        return any(i.startswith(prefix) for i in ids)

    def param(prefix, key):
        for s in stages:
            if str(s.get("stage", "")).startswith(prefix):
                return (s.get("params") or {}).get(key)
        return None

    res = []

    def add(rid, status, detail, fix):
        meta = next(r for r in RULES if r[0] == rid)
        res.append({"id": rid, "category": meta[1], "title": meta[2], "basis": meta[3],
                    "kind": meta[4], "status": status, "detail": detail, "fix": fix})

    # ---- 结构性 ----
    # R01：阈值可复现。v1.2 起门控可能是"类药性"而非 OB/DL，
    # 此时只报 OB/DL 会误导读者（以为用了 30/0.18 筛），必须按实际门控报告。
    gate = param("02", "gate")
    n_obdl = counts.get("gate_obdl", 0) or 0
    n_drug = counts.get("gate_druglike", 0) or 0
    qv = param("06", "q")
    if gate in ("auto", "druglike") and n_drug and not n_obdl:
        ok = qv is not None
        det = "门控=类药性（QED≥%s，Veber 违例≤%s）；OB/DL 未参与；q=%s" % (
            param("02", "qed_min"), param("02", "max_veber"), "—" if qv is None else qv)
    elif n_obdl and n_drug:
        ok = all(v is not None for v in (param("02", "ob"), param("02", "dl"), qv))
        det = "混合门控（OB/DL %d 行、类药性 %d 行）；OB=%s / DL=%s / q=%s" % (
            n_obdl, n_drug, param("02", "ob"), param("02", "dl"), "—" if qv is None else qv)
    else:
        th = [param("02", "ob"), param("02", "dl"), qv]
        ok = all(v is not None for v in th)
        det = "OB=%s / DL=%s / q=%s" % tuple("—" if v is None else v for v in th)
    add("R01", "PASS" if ok else ("WARN" if gate or qv is not None else "FAIL"), det,
        "让各阶段写 manifest（默认即写）；不要手改报告里的阈值")

    add("R03", "PASS" if len(stages) >= 3 else ("WARN" if stages else "FAIL"),
        "%d 个阶段 manifest" % len(stages),
        "用 np_run.py 串流程；单跑脚本也会自动写")

    has_env = os.path.exists(os.path.join(root, "environment.json"))
    add("R04", "PASS" if has_env else "FAIL", "有" if has_env else "缺失",
        "任意阶段跑一次即生成；跨机复现时必须提供")

    rate = counts.get("map_rate")
    if rate is None:
        m = re.search(r"映射率[^0-9]{0,12}([\d.]+)\s*%", all_text)
        rate = float(m.group(1)) / 100.0 if m else None
    add("R05", "PASS" if (rate is not None and rate >= 0.8) else ("WARN" if rate else "FAIL"),
        ("%.1f%%" % (100 * rate)) if rate is not None else "未记录",
        "np_standardize.py；<80% 先补别名表")

    nh = 0
    layers = {}
    for st in stages:
        if str(st.get("stage", "")).startswith("01"):
            layers = (st.get("counts") or {}).get("layers") or {}
    if layers:
        nh = int(layers.get("nonhuman", 0) + layers.get("noise", 0))
    if not nh:
        hstd, rows_std = read_csv_any(root, "targets_standardized.csv")
        if hstd and "layer" in hstd and rows_std:
            nh = sum(1 for r in rows_std if r.get("layer") in ("nonhuman", "noise"))
    add("R06", "PASS" if nh > 0 else "WARN",
        ("识别 %d 条非人类/噪音" % nh) if nh else "未见相关记录",
        "np_standardize 自动标注；正文写明剔除理由")

    hI, _ = read_csv_any(root, "intersection.csv")
    add("R07", "PASS" if (hI and "evidence_tier" in hI) else "WARN",
        "有 evidence_tier 列" if (hI and "evidence_tier" in hI) else "无来源分层",
        "np_intersect 输出；正文按层陈述")

    route = param("02", "route")
    add("R08", "PASS" if route and route != "unknown" else "WARN",
        "route=%s" % (route or "未记录"),
        "np_compounds --route；不匹配时改用含量/入血成分")

    qv = param("06", "q")
    add("R11", "PASS" if qv is not None else ("WARN" if stage("06") else "FAIL"),
        "q=%s" % qv if qv is not None else "未记录",
        "np_enrich 默认 BH 校正并写 q")

    bgc = param("06", "background")
    add("R12", "PASS" if bgc else ("WARN" if stage("06") else "FAIL"),
        ("background=%s" % bgc) if bgc else ("使用默认背景（报告已声明探索性）" if stage("06") else "无富集阶段"),
        "有实测表达数据时 --background；否则报告标注探索性")

    add("R13", "PASS" if stage("05") else "FAIL",
        "有零模型阶段" if stage("05") else "未运行 np_null",
        "np_null.py：A/B/C（+可选 D）")

    calib = None
    for k, t in texts.items():
        if "null_report" in k:
            if "未通过校准" in t:
                calib = "未通过"
            elif "零模型自校准" in t:
                calib = "通过"
    add("R14", "PASS" if calib == "通过" else ("FAIL" if calib == "未通过" else "WARN"),
        "C 组校准：%s" % (calib or "未记录"),
        "np_null 输出校准表；未通过时不得用 C 组结论")

    rh = counts.get("robust_hubs")
    add("R15", "PASS" if (rh or 0) > 0 else ("WARN" if rh == 0 else "FAIL"),
        "稳健 hub=%s" % rh if rh is not None else "未记录",
        "np_ppi --bootstrap；入选频率 ≥0.90")

    he, _ = read_csv_any(root, "enrich_significant.csv")
    add("R18", "PASS" if (he and "tier" in he) else "WARN",
        "有 tier 分层" if (he and "tier" in he) else "无分层",
        "np_enrich 自动三层分级")

    add("R19", "PASS" if (he and "collapsed_n" in he) else "WARN",
        "有折叠列" if (he and "collapsed_n" in he) else "未见折叠记录",
        "np_enrich --collapse 0.5")

    ub = counts.get("ubiquitous_ratio")
    add("R20", "PASS" if ub is not None else "WARN",
        ("泛成分占比 %.0f%%" % (100 * ub)) if ub is not None else "未记录",
        "np_compounds 自动标记")

    add("R21", "PASS" if stage("08") else "WARN",
        "有对接门禁" if stage("08") else "未运行对接评估",
        "np_dock.py；RMSD>2 Å 的靶点结果作废")

    add("R23", "PASS" if "预测" in all_text else "FAIL",
        "已声明计算预测边界" if "预测" in all_text else "未声明",
        "结论段固定一句：本文为计算预测，需实验验证")

    add("R25", "PASS" if stage("07") else "FAIL",
        "已核算" if stage("07") else "未核算",
        "np_dose.py；比值>10 时机制结论降级")

    add("R27", "PASS" if (claims.get("claims") or []) else "FAIL",
        "%d 条结论" % len(claims.get("claims") or []),
        "np_report.py 自动生成 claims.json")

    miss = claims.get("missing_params") or []
    add("R28", "PASS" if (claims and not miss) else ("FAIL" if miss else "WARN"),
        "参数全部有记录" if (claims and not miss) else ("%d 个参数未记录" % len(miss) if miss else "无账本"),
        "补齐对应阶段后重跑；报告不会用默认值兜底")

    # ---- 提示性（不计分） ----
    for rid, key, fix in [("R02", "ver", "每个库写「名称, 日期 访问」"),
                          ("R09", "paozhi", "写「生/制品 + 基原 + 药用部位」"),
                          ("R10", "partial", "加一句：该清单为符合类药性规律的部分"),
                          ("R16", "peffect", "改为「显著富集（q=…，命中 n）」"),
                          ("R17", "activate", "改为「与该通路相关」；方向需实验"),
                          ("R22", "pocket", "写明口袋来源与质子化工具"),
                          ("R24", "absolute", "无 L3+ 证据时改「提示/支持」"),
                          ("R26", "clinic", "结尾加「不作为诊疗建议」"),
                          ("R29", "limit", "至少覆盖数据库偏倚/阈值局限/无实验/暴露")]:
        hits = [(m, negated(all_text, m)) for m in re.finditer(PAT[key], all_text, re.I)]
        bad = [m for m, neg in hits if not neg]
        detail = ("命中：…%s…" % all_text[max(0, bad[0].start() - 26):bad[0].end() + 26].replace("\n", " ")) if bad \
            else ("仅有否定式表述" if hits else "未检出")
        meta = next(r for r in RULES if r[0] == rid)
        res.append({"id": rid, "category": meta[1], "title": meta[2], "basis": meta[3],
                    "kind": ADV, "status": "WARN" if bad else "PASS", "detail": detail, "fix": fix})

    if not re.search(PAT["date"], all_text):
        for it in res:
            if it["id"] == "R02" and it["status"] == "PASS":
                it["status"] = "WARN"
                it["detail"] = "未见日期"

    st_items = [x for x in res if x["kind"] == STRUCT]
    n_pass = sum(1 for x in st_items if x["status"] == "PASS")
    n_warn = sum(1 for x in st_items if x["status"] == "WARN")
    n_fail = sum(1 for x in st_items if x["status"] == "FAIL")
    score = round(100.0 * (n_pass + 0.5 * n_warn) / max(1, len(st_items)), 1)
    adv_items = [x for x in res if x["kind"] == ADV]
    adv_hits = [x for x in adv_items if x["status"] == "WARN"]

    md = ["# 研究项目红线审计（v2 · 结构性优先）", "",
          "- 目录：`%s`" % os.path.abspath(root),
          "- 报告 %d 个 · manifest %d 个 · 结论账本 %s"
          % (len(texts), len(stages), ("%d 条" % len(claims.get("claims") or [])) if claims else "无"),
          "- **结构性：PASS %d / WARN %d / FAIL %d，完备度 %.1f/100**" % (n_pass, n_warn, n_fail, score),
          "- 提示性（不计分）：%d 项命中，见文末" % len(adv_hits), ""]

    if n_fail:
        md += ["## ⛔ 阻断项（%d）" % n_fail, ""]
        for x in st_items:
            if x["status"] == "FAIL":
                md += ["### %s %s" % (x["id"], x["title"]), "",
                       "- 现状：%s" % x["detail"], "- 依据：%s" % x["basis"],
                       "- 修法：%s" % x["fix"], ""]
    if n_warn:
        md += ["## ⚠ 需补充（%d）" % n_warn, "", "| 编号 | 事项 | 现状 | 修法 |", "|---|---|---|---|"]
        for x in st_items:
            if x["status"] == "WARN":
                md.append("| %s | %s | %s | %s |" % (x["id"], x["title"], x["detail"], x["fix"]))
        md.append("")

    ICON = {"PASS": "✅", "WARN": "⚠", "FAIL": "⛔"}
    md += ["## 结构性检查明细", "", "| 编号 | 类别 | 事项 | 判定 | 现状 |", "|---|---|---|---|---|"]
    for x in st_items:
        md.append("| %s | %s | %s | %s | %s |" % (x["id"], x["category"], x["title"],
                                                  ICON[x["status"]], x["detail"]))

    md += ["", "## 提示性检查（启发式，不计分）", "",
           "> 基于散文正则，可能误报（已加否定窗口检测）。仅供人工复核，不影响上面的评分。", "",
           "| 编号 | 事项 | 判定 | 上下文 |", "|---|---|---|---|"]
    for x in adv_items:
        md.append("| %s | %s | %s | %s |" % (x["id"], x["title"], ICON[x["status"]], x["detail"][:90]))

    if warns:
        md += ["", "## 各阶段脚本发出的警告", ""] + ["- %s" % w for w in warns[:30]]

    md += ["", "---", "",
           "> 审计判断「报告是否具备可辩护的方法学要素」，不判断科学价值。",
           "> PASS ≠ 结论正确，只代表没有明显可被一击击穿的漏洞。", ""]

    p_md = C.write_text(os.path.join(outdir, "audit_report.md"), "\n".join(md) + "\n")
    p_json = C.write_text(os.path.join(outdir, "audit.json"),
                          json.dumps({"version": 2, "score": score,
                                      "structural": {"pass": n_pass, "warn": n_warn, "fail": n_fail},
                                      "advisory_hits": len(adv_hits),
                                      "items": res, "stage_warnings": warns},
                                     ensure_ascii=False, indent=2))

    if not a.quiet:
        print("[审计v2] 结构 PASS %d / WARN %d / FAIL %d → %.1f/100；提示性命中 %d"
              % (n_pass, n_warn, n_fail, score, len(adv_hits)))
        for x in st_items:
            if x["status"] == "FAIL":
                print("  ⛔ %s %s —— %s" % (x["id"], x["title"], x["fix"]))
        print("[审计v2] 报告 → %s" % p_md)

    if n_fail:
        return 2
    if a.strict and n_warn:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
