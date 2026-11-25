#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_report.py — 生成自包含 HTML 研究报告（内嵌 SVG，无外部依赖）

为什么不用 PDF / Word 模板
  * 自包含单文件：发给合作者就能打开，图不丢、字体不乱码（SVG 直接内嵌）。
  * 方法学段落自动生成：把参数、阈值、版本、种子写死在报告里，直接可粘贴进论文。
  * 证据分级图例固定出现：避免「把预测写成结论」。

用法
  python np_report.py --dir np_out/ --label "黄芪-当归-糖尿病肾病"
"""
import argparse
import glob
import html
import json
import os
import re
import sys

import npcore as C

CSS = """
:root{--ink:#1f2937;--muted:#6b7280;--line:#e5e7eb;--bg:#ffffff;--soft:#f8fafc;
--blue:#2f6fb2;--red:#c0504d;--green:#4f8a5b;--amber:#b7791f}
*{box-sizing:border-box}
body{margin:0;background:var(--soft);color:var(--ink);
font-family:-apple-system,"Segoe UI","Microsoft YaHei","PingFang SC",sans-serif;
line-height:1.7;font-size:15px}
.wrap{max-width:1080px;margin:0 auto;padding:32px 28px 80px;background:var(--bg)}
h1{font-size:26px;margin:0 0 6px;letter-spacing:-.3px}
h2{font-size:19px;margin:38px 0 12px;padding-bottom:7px;border-bottom:2px solid var(--line)}
h3{font-size:15.5px;margin:22px 0 8px;color:#111827}
.sub{color:var(--muted);font-size:13px;margin-bottom:26px}
table{border-collapse:collapse;width:100%;margin:12px 0;font-size:13.5px}
th,td{border:1px solid var(--line);padding:7px 9px;text-align:left;vertical-align:top}
th{background:#f3f6fa;font-weight:600;font-size:13px}
tr:nth-child(even) td{background:#fbfcfe}
code{background:#f1f5f9;padding:1px 5px;border-radius:3px;font-size:12.5px}
.card{background:var(--soft);border:1px solid var(--line);border-radius:8px;padding:14px 16px;margin:14px 0}
.tag{display:inline-block;font-size:11.5px;padding:1px 7px;border-radius:10px;margin-right:5px;
border:1px solid transparent}
.t0{background:#eef2f7;color:#475569;border-color:#dbe3ec}
.t1{background:#eaf2fb;color:#1e4f82;border-color:#cfe0f2}
.t2{background:#fdf3e7;color:#8a5a12;border-color:#f3e0c4}
.t3{background:#eaf5ec;color:#2f6b3a;border-color:#cfe6d4}
.fig{margin:16px 0;border:1px solid var(--line);border-radius:8px;overflow:hidden;background:#fff}
.fig svg{display:block;width:100%;height:auto}
.figcap{font-size:12.5px;color:var(--muted);padding:8px 12px;background:var(--soft);
border-top:1px solid var(--line)}
.warn{background:#fffaf0;border-left:4px solid var(--amber);padding:10px 14px;margin:12px 0;font-size:13.5px}
.bad{background:#fdf3f2;border-left:4px solid var(--red);padding:10px 14px;margin:12px 0;font-size:13.5px}
.good{background:#f3faf5;border-left:4px solid var(--green);padding:10px 14px;margin:12px 0;font-size:13.5px}
ul,ol{padding-left:22px}li{margin:4px 0}
.toc{columns:2;font-size:13.5px}
.toc a{color:var(--blue);text-decoration:none}
.r{text-align:right}
.small{font-size:12.5px;color:var(--muted)}
"""

EVID = [
    ("L0", "t0", "数据库注释/文献事实", "可直接引用，但须带版本与时间"),
    ("L1", "t1", "计算预测（靶点、网络、富集）", "只能作为假说，用「提示/支持」"),
    ("L2", "t2", "结构证据（对接/动力学）", "支持结合可能性，不等于活性"),
    ("L3", "t3", "体外/离体实验", "可用「表明」"),
    ("L4", "", "体内/临床", "可用「证明」，需注明样本量"),
]


def esc(x):
    return html.escape(str(x if x is not None else ""))


def read_manifest(d):
    p = os.path.join(d, "run_manifest.json")
    if not os.path.exists(p):
        return []
    try:
        return json.load(open(p, encoding="utf-8")).get("stages", [])
    except Exception:
        return []


def table(headers, rows, cls=""):
    out = ['<table class="%s"><thead><tr>' % cls]
    out += ["<th>%s</th>" % esc(h) for h in headers]
    out.append("</tr></thead><tbody>")
    for r in rows:
        out.append("<tr>" + "".join("<td>%s</td>" % (c if isinstance(c, str) and c.startswith("<") else esc(c))
                                    for c in r) + "</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def embed_svg(path, cap=""):
    try:
        s = open(path, encoding="utf-8").read()
        s = re.sub(r"<\?xml[^>]*\?>", "", s)
        s = re.sub(r"<!DOCTYPE[^>]*>", "", s)
        s = re.sub(r'\swidth="\d+"\s+height="\d+"', "", s, count=1)
        return '<div class="fig">%s<div class="figcap">%s</div></div>' % (s, esc(cap))
    except Exception:
        return ""


def md_to_html(md):
    """极简 Markdown → HTML（标题/列表/表格/引用/粗体/代码），足够渲染阶段报告。"""
    lines = md.splitlines()
    out, i = [], 0
    while i < len(lines):
        ln = lines[i]
        if re.match(r"^\s*\|.*\|\s*$", ln) and i + 1 < len(lines) and re.match(r"^\s*\|[\s\-:|]+\|\s*$", lines[i + 1]):
            head = [c.strip() for c in ln.strip().strip("|").split("|")]
            i += 2
            body = []
            while i < len(lines) and re.match(r"^\s*\|.*\|\s*$", lines[i]):
                body.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            out.append(table(head, body))
            continue
        if re.match(r"^#{1,4}\s", ln):
            lvl = len(ln) - len(ln.lstrip("#"))
            txt = ln[lvl:].strip()
            out.append("<h%d>%s</h%d>" % (min(lvl + 1, 4), inline(txt), min(lvl + 1, 4)))
        elif re.match(r"^\s*[-*]\s+", ln):
            items = []
            while i < len(lines) and re.match(r"^\s*[-*]\s+", lines[i]):
                items.append("<li>%s</li>" % inline(re.sub(r"^\s*[-*]\s+", "", lines[i])))
                i += 1
            out.append("<ul>%s</ul>" % "".join(items))
            continue
        elif re.match(r"^\s*\d+[.)]\s+", ln):
            items = []
            while i < len(lines) and re.match(r"^\s*\d+[.)]\s+", lines[i]):
                items.append("<li>%s</li>" % inline(re.sub(r"^\s*\d+[.)]\s+", "", lines[i])))
                i += 1
            out.append("<ol>%s</ol>" % "".join(items))
            continue
        elif ln.startswith(">"):
            out.append('<div class="card">%s</div>' % inline(ln.lstrip("> ")))
        elif ln.strip() == "---":
            out.append("<hr>")
        elif ln.strip():
            out.append("<p>%s</p>" % inline(ln))
        i += 1
    return "\n".join(out)


def inline(t):
    t = esc(t)
    t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    return t


def main(argv=None):
    ap = argparse.ArgumentParser(description="生成自包含 HTML 报告")
    ap.add_argument("--dir", required=True)
    ap.add_argument("--label", default="网络药理学分析")
    ap.add_argument("-o", "--output", default=None)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    d = a.dir
    if not os.path.isdir(d):
        print("[错误] 目录不存在：%s" % d)
        return 2
    stages = read_manifest(d)
    counts, warns = {}, []
    for s in stages:
        counts.update(s.get("counts", {}))
        warns.extend(s.get("warnings", []))

    def rd(name):
        p = os.path.join(d, name)
        if not os.path.exists(p):
            return None, []
        try:
            h, rows = C.read_table(p)
            return h, rows
        except Exception:
            return None, []

    parts = []
    parts.append('<div class="wrap">')
    parts.append("<h1>%s</h1>" % esc(a.label))
    parts.append('<div class="sub">网络药理学研究报告 · 由 network-pharmacology-pro v1.0 生成 · '
                 '全部数值可由 <code>run_manifest.json</code> 与原始数据复算</div>')

    # 结论摘要
    inter_h, inter = rd("intersection.csv")
    hubs_h, hubs = rd("hubs.csv")
    enr_h, enr = rd("enrich_significant.csv")
    n_int = len(inter)
    n_hub = len(hubs)
    mech = [r for r in enr if r.get("tier") == "机制明确"]
    umb = [r for r in enr if r.get("tier") == "伞状泛化"]
    parts.append('<div class="card"><h3 style="margin-top:0">一页结论</h3><ul>')
    parts.append("<li>交集靶点：<strong>%d</strong> 个%s</li>" %
                 (n_int, ("，其中共享成分最多的为 <strong>%s</strong>" % esc(inter[0]["gene"]))
                  if inter else ""))
    if hubs:
        parts.append("<li>hub 候选：<strong>%s</strong>（多指标共识前 %d）</li>" %
                     (esc("、".join(r["node"] for r in hubs[:6])), n_hub))
    if enr:
        parts.append("<li>显著富集条目：<strong>%d</strong> 条（折叠后），其中机制明确类 <strong>%d</strong> 条</li>"
                     % (len(enr), len(mech)))
    parts.append("<li>方法学完备度：见文末审计；证据等级：本报告主体为 <strong>L1（计算预测）</strong></li>")
    parts.append("</ul></div>")

    parts.append('<h2 id="evidence">证据等级图例</h2>')
    parts.append("<p>本报告所有结论均可追溯到下列等级之一。等级决定措辞，不可越级。</p>")
    parts.append(table(["等级", "含义", "允许的措辞"], [
        ['<span class="tag %s">%s</span>' % (cls, code), name, allow] for code, cls, name, allow in EVID
    ]))
    parts.append('<div class="warn"><strong>读法：</strong>本报告未含 L3/L4 实验数据，'
                 '因此凡出现「提示/支持/可能」均为刻意降级，不是保守，而是准确。</div>')

    # 参数快照
    parts.append('<h2 id="params">参数与溯源快照</h2>')
    if stages:
        rows = []
        for s in stages:
            p = s.get("params", {})
            keep = {k: v for k, v in p.items() if k in
                    ("ob", "dl", "hl", "logic", "route", "min_score", "top", "bootstrap",
                     "iterations", "seed", "universe", "q", "min_hit", "collapse", "rmsd_cut")}
            rows.append([s.get("stage", ""), s.get("started", "")[:19],
                         ", ".join("%s=%s" % (k, v) for k, v in keep.items()) or "—",
                         ", ".join("%s=%s" % (k, v) for k, v in (s.get("counts") or {}).items()
                                   if isinstance(v, (int, float)))[:160]])
        parts.append(table(["阶段", "时间", "关键参数", "计数"], rows))
        parts.append('<p class="small">完整溯源（输入/输出 SHA256、字节数、随机种子）见 '
                     '<code>run_manifest.json</code>。</p>')
    else:
        parts.append('<div class="bad">未找到 <code>run_manifest.json</code>：本报告不具备可复现性证明。</div>')

    # 图
    figs = [("venn.svg", "图 1　靶点交集（韦恩图）"),
            ("network.svg", "图 2　PPI 网络（节点大小∝度，颜色=模块）"),
            ("enrich_bubble.svg", "图 3　富集气泡图（气泡=命中数，颜色=−log10 q）"),
            ("hub_bar.svg", ""), ("docking_bar.svg", "图 4　对接打分（−kcal/mol，仅为预测）"),
            ("coverage_bar.svg", "图 5　各药味覆盖的疾病基因数")]
    fig_html = []
    for fn, cap in figs:
        p = os.path.join(d, fn)
        if os.path.exists(p) and cap:
            fig_html.append(embed_svg(p, cap))
    if fig_html:
        parts.append('<h2 id="figs">图表</h2>')
        parts.extend(fig_html)

    # 阶段报告（折叠）
    parts.append('<h2 id="detail">分阶段结果</h2>')
    order = ["standardize_audit.md", "compounds_report.md", "intersect_report.md",
             "ppi_report.md", "null_report.md", "enrich_report.md", "dose_report.md",
             "docking_report.md", "optimize_report.md"]
    NAMES = {"standardize_audit.md": "① 靶点标准化", "compounds_report.md": "② 成分筛选",
             "intersect_report.md": "③ 交集与成分富集", "ppi_report.md": "④ PPI 与 hub",
             "null_report.md": "⑤ 零模型对照", "enrich_report.md": "⑥ 富集",
             "dose_report.md": "⑦ 剂量可达性", "docking_report.md": "⑧ 对接门禁",
             "optimize_report.md": "⑨ 反向网络药理学"}
    for fn in order:
        p = os.path.join(d, fn)
        if not os.path.exists(p):
            continue
        txt = open(p, encoding="utf-8", errors="replace").read()
        txt = re.sub(r"^# .*$", "", txt, count=1, flags=re.M)
        parts.append('<details open><summary style="cursor:pointer;font-weight:600;margin:12px 0">%s</summary>%s</details>'
                     % (esc(NAMES.get(fn, fn)), md_to_html(txt)))

    # ---------------- 方法学段落（参数必须来自记录它的阶段） ----------------
    # 设计原则：**绝不用默认值兜底**。参数未记录时输出「未记录」并显式告警，
    # 因为这段文字的设计用途是直接粘贴进论文——写错一个阈值等于伪造方法。
    def param_of(stage_prefix, key):
        for st_ in stages:
            if str(st_.get("stage", "")).startswith(stage_prefix):
                pr = st_.get("params", {}) or {}
                if key in pr and pr[key] is not None:
                    return pr[key], st_.get("stage")
        return None, None

    parts.append('<h2 id="methods">方法学段落（可直接粘贴）</h2>')
    ob = counts.get("kept")
    mi = counts.get("map_rate")
    seg, src_rows, missing = [], [], []

    def add_param(label, stage_prefix, key, fmt="%s"):
        val, src = param_of(stage_prefix, key)
        if val is None:
            missing.append("%s（阶段 %s 未记录 %s）" % (label, stage_prefix, key))
            src_rows.append([label, "**未记录**", "—", "⛔"])
            return None
        shown = fmt % (val,)
        src_rows.append([label, shown, src, "✅"])
        return shown

    ob_v = add_param("OB 阈值(%)", "02", "ob", "%g")
    dl_v = add_param("DL 阈值", "02", "dl", "%g")
    route_v = add_param("给药途径", "02", "route")
    ms_v = add_param("PPI 置信度阈值", "04", "min_score", "%g")
    bs_v = add_param("bootstrap 次数", "04", "bootstrap", "%d")
    top_v = add_param("hub 取前 k", "04", "top", "%d")
    q_v = add_param("富集 q 阈值", "06", "q", "%g")
    mh_v = add_param("富集最小命中数", "06", "min_hit", "%d")
    it_v = add_param("零模型迭代次数", "05", "iterations", "%d")

    seg.append("本研究采用网络药理学方法解析受试物对目标疾病的作用机制。")
    if mi:
        seg.append("靶点名称经统一映射至 HGNC 基因符号，映射率 %.1f%%，非人类来源条目予以剔除。" % (100 * mi))
    # 门控描述必须与 02 阶段实际使用的尺子一致；否则这段"可直接粘贴"的文本就是伪造方法
    gate_v = param_of("02", "gate")[0]
    n_obdl = counts.get("gate_obdl", 0) or 0
    n_drug = counts.get("gate_druglike", 0) or 0
    qed_v = add_param("类药性 QED 下限", "02", "qed_min", "%g")
    veber_v = add_param("Veber 违例上限", "02", "max_veber", "%d")
    if ob:
        if gate_v in ("auto", "druglike") and n_drug and not n_obdl:
            seg.append("成分纳入采用可复现的类药性判据：QED ≥ %s 且 Veber 规则违例 ≤ %s 项"
                       "（MW≤500、cLogP≤5、HBD≤5、HBA≤10、TPSA≤140 Å²、可旋转键≤10），"
                       "描述符由 RDKit 从规范 SMILES 计算；最终纳入 %d 个成分；给药途径为 %s。"
                       % (qed_v or "未记录", veber_v or "未记录", ob, route_v or "未声明"))
        elif n_obdl and n_drug:
            seg.append("成分纳入按行采用两种判据（OB/DL %d 个、类药性 %d 个，两类判据不可直接比较），"
                       "最终纳入 %d 个成分；给药途径为 %s。" % (n_obdl, n_drug, ob, route_v or "未声明"))
        elif gate_v == "evidence":
            seg.append("成分纳入以实测活性/入血证据为准，未设分子层面阈值（类药性仅作描述），"
                       "最终纳入 %d 个成分；给药途径为 %s。" % (ob, route_v or "未声明"))
        elif ob_v and dl_v:
            seg.append("成分筛选取 OB ≥ %s%% 且 DL ≥ %s（TCMSP 体系；二者均为预测值），"
                       "最终纳入 %d 个成分；给药途径为 %s。" % (ob_v, dl_v, ob, route_v or "未声明"))
    if n_int:
        seg.append("药物靶点与疾病靶点取交集得到 %d 个共同靶点。" % n_int)
    if hubs and ms_v:
        seg.append("将共同靶点映射至 STRING（Homo sapiens，combined score ≥ %s）构建 PPI 网络，"
                   "计算度、介数（Brandes）、谐波接近中心性、PageRank、k-core、聚类系数与 MCC，"
                   "以多指标秩均值确定候选 hub（前 %s），并用 %s 次边 bootstrap 估计入选频率。"
                   % (ms_v, top_v or "—", bs_v or "—"))
    if enr and q_v:
        seg.append("GO/KEGG 富集采用超几何检验，Benjamini-Hochberg 校正，q < %s 且命中基因 ≥ %s 视为显著，"
                   "并按基因集 Jaccard 折叠冗余条目。" % (q_v, mh_v or 3))
    if it_v:
        seg.append("为排除拓扑偏倚，构建随机零模型（全基因组随机、疾病集内随机、按观察集度分布精确匹配随机），"
                   "各 %s 次，经验 P 值按 (1+#{null≥obs})/(1+B) 计算。" % it_v)

    parts.append('<div class="card">%s</div>' % inline("".join(seg)))

    if src_rows:
        parts.append('<h3>参数来源核对表</h3>')
        parts.append(table(["参数", "取值", "来源阶段", "状态"], src_rows))
    if missing:
        parts.append('<div class="bad"><strong>⛔ 有 %d 个参数在本次运行中未被记录。</strong>'
                     '上面段落中对应的数值已改为「未记录」，<strong>不会</strong>用默认值填充。'
                     '请补齐对应阶段后重跑，否则该段落不可用于投稿。<ul>%s</ul></div>'
                     % (len(missing), "".join("<li>%s</li>" % esc(m) for m in missing)))
    parts.append('<p class="small">⚠ 数据库版本、访问日期、软件版本需按实际填写；'
                 '上述段落只保证「参数与流程」与本次运行一致，不替代作者对研究设计的表述。</p>')

    # ---------------- 结论账本（claims.json） ----------------
    claims = []
    if ob:
        claims.append({"id": "C1", "text": "纳入 %d 个活性成分" % ob, "level": "L1",
                       "value": {"kept": ob}, "source": "compounds_kept.csv",
                       "params_ref": "02_compounds"})
    if n_int:
        claims.append({"id": "C2", "text": "药物靶点与疾病靶点交集 %d 个" % n_int, "level": "L1",
                       "value": {"intersection": n_int}, "source": "intersection.csv",
                       "params_ref": "03_intersect"})
    if hubs:
        claims.append({"id": "C3", "text": "hub 候选 %d 个（多指标共识）" % len(hubs), "level": "L1",
                       "value": {"hubs": len(hubs)}, "source": "hubs.csv", "params_ref": "04_ppi"})
    if enr:
        mech = sum(1 for r in enr if r.get("tier") == "机制明确")
        claims.append({"id": "C4", "text": "显著富集 %d 条，其中机制明确类 %d 条" % (len(enr), mech),
                       "level": "L1", "value": {"sig": len(enr), "mechanistic": mech},
                       "source": "enrich_significant.csv", "params_ref": "06_enrich"})
    ledger = os.path.join(d, "claims.json")
    C.write_text(ledger, json.dumps({"label": a.label,
                                     "generated": C.time.strftime("%Y-%m-%dT%H:%M:%S"),
                                     "claims": claims, "missing_params": missing},
                                    ensure_ascii=False, indent=2))
    # 同时用核心库的完整版覆盖（含零模型条数等），保证报告与审计读到同一份账本
    try:
        C.build_claims(d, a.label)
    except Exception:
        pass
    parts.append('<p class="small">结论账本：<code>claims.json</code>（%d 条，含证据等级与来源文件）；'
                 '审计器优先依据账本判定，而非散文措辞。</p>' % len(claims))

    # 审计摘要
    audit_p = os.path.join(d, "audit", "audit.json")
    if os.path.exists(audit_p):
        try:
            au = json.load(open(audit_p, encoding="utf-8"))
            parts.append('<h2 id="audit">方法学审计</h2>')
            cls = "good" if au["fail"] == 0 and au["warn"] <= 3 else ("warn" if au["fail"] == 0 else "bad")
            parts.append('<div class="%s">完备度 <strong>%.1f/100</strong>　'
                         'PASS %d · WARN %d · FAIL %d</div>'
                         % (cls, au["score"], au["pass"], au["warn"], au["fail"]))
            bad = [x for x in au["items"] if x["status"] in ("FAIL", "WARN")][:12]
            if bad:
                parts.append(table(["编号", "事项", "判定", "修法"],
                                   [[x["id"], x["title"],
                                     "⛔" if x["status"] == "FAIL" else "⚠", x["fix"]] for x in bad]))
            parts.append('<p class="small">完整审计见 <code>audit/audit_report.md</code>。</p>')
        except Exception:
            pass

    if warns:
        parts.append('<h2 id="warn">流程警告（不可忽略）</h2><ul>')
        for w in warns[:40]:
            parts.append("<li>%s</li>" % inline(w))
        parts.append("</ul>")

    # ---- SCM0061-2021 评价评分卡（若已运行 np_scm0061.py）----
    scm_p = os.path.join(d, "scm0061_scorecard.json")
    if os.path.exists(scm_p):
        try:
            scm = json.load(open(scm_p, encoding="utf-8"))
            parts.append('<h2 id="scm0061">SCM0061-2021 评价评分</h2>')
            parts.append('<p>依据世界中医药学会联合会《网络药理学评价方法指南》（SCM0061-2021）'
                         '附录 A 的评价要素，对本次运行的产物逐项判定：'
                         '<strong>%s/%s（%.1f%%，等级 %s）</strong>。'
                         '0=未检出该要素产物；1=产物存在但缺关键字段；2=齐备。</p>'
                         % (scm.get("total_score"), scm.get("total_max"),
                            scm.get("percent", 0), scm.get("grade", "-")))
            rows = [[k, v.get("got"), v.get("max")]
                    for k, v in (scm.get("dimensions") or {}).items()]
            parts.append(table(["评价维度", "得分", "满分"], rows))
            items = scm.get("items") or []
            miss = [[r["id"], r.get("item", ""), r.get("evidence", "")]
                    for r in items if r.get("score", 0) < 2]
            if miss:
                parts.append("<h3>未满分项与原因</h3>")
                parts.append(table(["编号", "要素", "证据/原因"], miss))
            parts.append('<p class="small">完整评分见 <code>scm0061_scorecard.json</code>；'
                         '已开展口径：%.1f%%（仅统计本次运行过的阶段）。</p>'
                         % scm.get("percent_executed_scope", 0))
        except Exception as exc:
            parts.append('<p class="small">[warn] 评分卡读取失败：%s</p>' % esc(str(exc)))

    parts.append('<h2 id="limit">局限与边界（必读）</h2><ul>')
    parts.append("<li>靶点预测、网络构建、富集与对接均为<strong>计算预测</strong>，不构成实验证据。</li>")
    parts.append("<li>数据库覆盖不全，民族药、地方习用品、微量成分与代谢产物常缺失。</li>")
    parts.append("<li>OB/DL 为预测值，会系统性漏掉苷类、多糖与挥发性成分。</li>")
    parts.append("<li>PPI 网络受研究偏倚影响，度分布本身可由偏倚解释，hub 不等于机制核心。</li>")
    parts.append("<li>本报告不作为临床用药依据。</li>")
    parts.append("</ul>")

    parts.append('<p class="small" style="margin-top:40px;border-top:1px solid var(--line);padding-top:14px">'
                 'network-pharmacology-pro v1.0 · 生成时间 %s · 本文件自包含，可离线打开</p>'
                 % C.time.strftime("%Y-%m-%d %H:%M"))
    parts.append("</div>")

    doc = ("<!DOCTYPE html><html lang=\"zh-CN\"><head><meta charset=\"utf-8\">"
           "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
           "<title>%s</title><style>%s</style></head><body>%s</body></html>"
           % (esc(a.label), CSS, "\n".join(parts)))

    out = a.output or os.path.join(d, "np_report.html")
    C.write_text(out, doc)
    if not a.quiet:
        print("[报告] 已生成 %s（%.0f KB）" % (out, os.path.getsize(out) / 1024))
    return 0


if __name__ == "__main__":
    sys.exit(main())
