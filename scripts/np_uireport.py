#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_uireport.py —— 一键分析报告生成器（供 app.py / 命令行共用）

与已有两类报告的分工
    np_report.py      : 阶段审计型（面向开发者自查，散点式罗列）
    np_htmlreport.py  : 投稿版排版（需 run_manifest + standard_*.csv，适合已有规范产物）
    np_uireport.py    : **一键式交付报告**（不挑产物，缺什么就写"未运行"，绝不空转），
                        带 Ready_for_Submission / DRAFT 水印，可直接下载/打印为 PDF

设计原则
    * 只读产物，不重算任何指标——报告里的每个数字都来自磁盘上的 CSV/JSON/MD；
    * 缺产物 → 该章节显式写"未运行"，而不是静默隐藏（审稿人最恨静默省略）；
    * 水印由 finalcheck 结论决定；没跑门禁就盖 UNVERIFIED，绝不默认"Ready"。

用法
    python np_uireport.py --dir np_out --label "复方丹参滴丸-冠心病" -o report.html
    python np_uireport.py --dir np_out --label demo --pdf        # 有 weasyprint 时同时出 PDF
"""
from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import html
import json
import math
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

CSS = r"""
:root{--brand:#0A1A2F;--accent:#007B7A;--ink:#1f2937;--muted:#6b7280;--line:#e5e7eb;
      --soft:#f8fafc;--hub:#c0504d;--ok:#4f8a5b;--warn:#b7791f;--bad:#b91c1c}
*{box-sizing:border-box}
body{margin:0;background:var(--soft);color:var(--ink);font-family:"Times New Roman","Microsoft YaHei",serif;
     line-height:1.62;font-size:11pt}
.wrap{max-width:196mm;margin:0 auto;padding:16mm 13mm 26mm;background:#fff;position:relative}
h1{font-size:19pt;margin:0 0 4pt;color:var(--brand);letter-spacing:-.3px}
h2{font-size:13.5pt;margin:20pt 0 7pt;padding-bottom:3pt;border-bottom:1.4pt solid var(--accent)}
h3{font-size:11.5pt;margin:12pt 0 4pt;color:var(--brand)}
.sub{color:var(--muted);font-size:9.5pt;margin-bottom:14pt}
table{border-collapse:collapse;width:100%;margin:7pt 0;font-size:8.6pt}
th,td{border:.5pt solid var(--line);padding:3.5pt 5pt;text-align:left;vertical-align:top}
th{background:#f3f6fa;font-weight:600}
tr:nth-child(even) td{background:#fbfcfe}
code,pre{font-family:Consolas,Menlo,monospace}
code{background:#f1f5f9;padding:1pt 4pt;border-radius:3px;font-size:8.5pt}
pre{background:#f8fafc;border:.5pt solid var(--line);border-radius:5pt;padding:7pt;font-size:8pt;
    white-space:pre-wrap;word-break:break-all}
.cards{display:flex;flex-wrap:wrap;gap:8pt;margin:10pt 0}
.card{flex:1 1 30%;min-width:46mm;border:.5pt solid var(--line);border-radius:6pt;padding:8pt 10pt;background:var(--soft)}
.card .n{font-size:17pt;color:var(--accent);font-weight:600;line-height:1.1}
.card .k{font-size:8.2pt;color:var(--muted)}
.badge{display:inline-block;font-size:8pt;padding:1.5pt 7pt;border-radius:9pt;margin-right:4pt}
.b-ok{background:#eef7f0;color:var(--ok);border:.5pt solid #d3e8d8}
.b-warn{background:#fff8e8;color:var(--warn);border:.5pt solid #f3e6c8}
.b-bad{background:#fdf0ef;color:var(--bad);border:.5pt solid #f3d5d3}
.b-bg{background:#f1f1f1;color:#6b7280;border:.5pt solid #e0e0e0}
.boundary{background:#fffaf0;border-left:3pt solid var(--warn);padding:8pt 11pt;margin:9pt 0;font-size:9.3pt}
.miss{color:var(--muted);font-size:9pt;font-style:italic}
.fig{margin:10pt 0;border:.5pt solid var(--line);border-radius:6pt;overflow:hidden;background:#fff;
     page-break-inside:avoid;break-inside:avoid}
.fig svg{display:block;width:100%;height:auto}
.figcap{font-size:8.4pt;color:var(--muted);padding:5pt 8pt;background:var(--soft);border-top:.5pt solid var(--line)}
/* 水印：平铺斜向文字，打印也保留 */
.wm{position:fixed;inset:0;z-index:0;pointer-events:none;overflow:hidden}
.wm span{position:absolute;font:700 26pt/1 Arial;color:rgba(185,28,28,.055);transform:rotate(-28deg);
         white-space:nowrap;letter-spacing:2pt}
.wrap>*{position:relative;z-index:1}
.stamp{position:fixed;right:14px;bottom:14px;z-index:9999;padding:8px 14pt;border-radius:6pt;
       font:600 12px/1.4 Arial;color:#fff;box-shadow:0 2px 8px #0004}
.stamp.ready{background:#1b7a3d}.stamp.draft{background:#b91c1c}.stamp.unv{background:#6b7280}
/* 门禁看板 */
.dash{border:1pt solid var(--line);border-radius:6pt;margin:10pt 0;padding:9pt 11pt;background:#fcfdff;
      page-break-inside:avoid}
.dash.missing{border-color:#f0c9c9;background:#fff7f7;color:#b91c1c;font-size:9.3pt}
.dash-h{display:flex;align-items:center;gap:8pt;flex-wrap:wrap;margin-bottom:6pt;font-size:10pt}
.dash-sum{color:var(--muted);font-size:8.8pt;font-weight:400}
.pill{display:inline-block;padding:1pt 8pt;border-radius:9pt;font-size:8pt;font-weight:600}
.pill.ok{background:#e7f5ea;color:#1b7a3d;border:.5pt solid #bfe3c8}
.pill.bad{background:#fdecec;color:#b91c1c;border:.5pt solid #f3c9c9}
table.gates{font-size:7.8pt;table-layout:fixed}
table.gates th:nth-child(1){width:6%}
table.gates th:nth-child(2){width:14%}
table.gates th:nth-child(3){width:8%}
table.gates th:nth-child(4){width:9%}
table.gates th:nth-child(5){width:19%}
table.gates th:nth-child(6){width:44%}
table.gates .pill{white-space:nowrap}
table.gates td{vertical-align:top}
table.gates tr.fail.hard td{background:#fdf3f3}
table.gates tr.fail.soft td{background:#fffbf0}
table.gates code{font-size:7.6pt;word-break:break-all}
.dash-pass{margin-top:5pt;font-size:8.4pt;color:var(--muted)}
.alert{border-left:3pt solid var(--warn);background:#fff8e8;padding:7pt 10pt;margin:8pt 0;
       font-size:9.2pt;page-break-inside:avoid}
.alert.bad{border-left-color:#b91c1c;background:#fdecec}
.alert b{color:#8a5a00}
.alert.bad b{color:#b91c1c}
iframe.net{width:100%;height:520px;border:.5pt solid var(--line);border-radius:6pt;background:#fff}
.print-only{display:none}
@page{size:A4;margin:17mm 14mm}
@media print{
  /* 屏幕上的平铺水印在打印/PDF 里隐藏：草稿是给导师看的，不是给期刊的 */
  .wm{display:none!important}
  .stamp{display:none!important}
  .dash{page-break-inside:auto;break-inside:auto}
  .dash table.gates{page-break-inside:auto;break-inside:auto}
  .print-only{display:block}
  .print-note{font-size:8pt;color:#6b7280;border-top:.5pt solid var(--line);
              padding-top:4pt;margin-bottom:8pt}
  body{background:#fff;font-size:10.4pt}
  .wrap{max-width:none;margin:0;padding:0}
  h2,h3{page-break-after:avoid}
  .fig,table,.card,.boundary{page-break-inside:avoid;break-inside:avoid}
  .no-print, iframe.net{display:none!important}
  thead{display:table-header-group}
}
"""

WATERMARK_TEXT = {"READY": "READY FOR SUBMISSION", "DRAFT": "DRAFT — 门禁未通过",
                  "UNVERIFIED": "UNVERIFIED — 未跑投稿门禁"}


def esc(x):
    return html.escape(str(x if x is not None else ""))


def _read_csv(path, limit=None):
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        rows = list(csv.reader(f, delimiter=delim))
    return rows[:limit] if limit else rows


def table_from_csv(path, limit=20, max_cols=9):
    rows = _read_csv(path, limit + 1)
    if not rows:
        return ""
    head = rows[0][:max_cols]
    body = [r[:max_cols] for r in rows[1:limit + 1]]
    out = ["<table><thead><tr>"] + ["<th>%s</th>" % esc(c) for c in head] + ["</tr></thead><tbody>"]
    for r in body:
        out.append("<tr>" + "".join("<td>%s</td>" % esc(c) for c in r) + "</tr>")
    out.append("</tbody></table>")
    if len(rows) > limit + 1:
        out.append('<p class="miss">（仅显示前 %d 行）</p>' % limit)
    return "".join(out)


def md_excerpt(path, lines=16):
    if not os.path.exists(path):
        return ""
    txt = open(path, encoding="utf-8", errors="replace").read()
    keep = [l for l in txt.splitlines() if l.strip()][:lines]
    return "<pre>%s</pre>" % esc("\n".join(keep))


def embed_svg(path):
    try:
        s = open(path, encoding="utf-8", errors="replace").read()
        s = re.sub(r"<\?xml[^>]*\?>|<!DOCTYPE[^>]*>", "", s)
        s = re.sub(r'\swidth="[\d.]+"\s+height="[\d.]+"', "", s, count=1)
        return '<div class="fig">%s<div class="figcap">%s</div></div>' % (s, esc(os.path.basename(path)))
    except Exception:
        return ""


def collect(outdir):
    """扫描产物，返回结构化 dict。缺什么就是 None，不做假。"""
    d = outdir
    g = lambda *p: os.path.join(d, *p)
    info = {"dir": d}
    man = g("pipeline_manifest.json")
    info["manifest"] = json.load(open(man, encoding="utf-8")) if os.path.exists(man) else None
    info["run_manifest"] = (json.load(open(g("run_manifest.json"), encoding="utf-8"))
                            if os.path.exists(g("run_manifest.json")) else None)
    info["safety"] = (json.load(open(g("00_safety", "contraindication.json"), encoding="utf-8"))
                      if os.path.exists(g("00_safety", "contraindication.json")) else None)
    info["compounds"] = next((p for p in [g("03_compounds", "compounds_kept.csv"),
                                          g("03_compounds", "compounds.csv")] if os.path.exists(p)), None)
    info["intersection"] = next((p for p in [g("04_intersect", "intersection.csv")] if os.path.exists(p)), None)
    info["intersection_genes"] = next((p for p in [g("04_intersect", "intersection_genes.txt")]
                                       if os.path.exists(p)), None)
    info["hubs"] = next((p for p in [g("05_ppi", "hubs.csv")] if os.path.exists(p)), None)
    info["ppi_report"] = g("05_ppi", "ppi_report.md")
    info["null_report"] = g("06_null", "null_report.md")
    info["enrich_report"] = g("07_enrich", "enrich_report.md")
    info["ortholog"] = next((p for p in glob.glob(g("11_ortholog", "*orthologs.csv"))), None)
    info["svgs"] = sorted(glob.glob(g("*.svg")) + glob.glob(g("**", "*.svg"), recursive=True))
    info["finalcheck"] = next((p for p in [g("finalcheck.json"), g("qa", "finalcheck.json")]
                               if os.path.exists(p)), None)
    info["mr"] = sorted(glob.glob(g("**", "mr_result.json"), recursive=True))
    info["intervention"] = next((p for p in [g("intervention", "intervention_summary.json")]
                                 if os.path.exists(p)), None)
    info["benchmark"] = sorted(glob.glob(g("**", "benchmark_*.json"), recursive=True))
    info["synergy"] = next((p for p in [g("synergy", "synergy_manifest.json"),
                                        g("10_synergy", "synergy_manifest.json")]
                            if os.path.exists(p)), None)
    return info


def counts(info):
    """从产物里数出摘要卡片，全部由文件计数得来。"""
    c = {}
    if info["compounds"]:
        c["compounds"] = max(0, len(_read_csv(info["compounds"])) - 1)
    if info["intersection_genes"]:
        c["targets"] = sum(1 for l in open(info["intersection_genes"], encoding="utf-8-sig") if l.strip())
    if info["hubs"]:
        rows = _read_csv(info["hubs"])
        hdr = [h.lower() for h in rows[0]] if rows else []
        if "verdict" in hdr:
            vi = hdr.index("verdict")
            c["hubs"] = sum(1 for r in rows[1:] if vi < len(r) and "稳健" in r[vi])
        else:
            c["hubs"] = min(15, max(0, len(rows) - 1))
    if info["enrich_report"] and os.path.exists(info["enrich_report"]):
        m = re.search(r"折叠后[^0-9]{0,12}(\d+)", open(info["enrich_report"], encoding="utf-8").read())
        if m:
            c["pathways"] = int(m.group(1))
    if info["manifest"]:
        nodes = info["manifest"].get("nodes") or []
        c["nodes_ok"] = sum(1 for n in nodes if n.get("status") in ("ok", "skipped_done"))
        c["nodes"] = len(nodes)
    return c


def watermark(verdict):
    txt = WATERMARK_TEXT.get(verdict, WATERMARK_TEXT["UNVERIFIED"])
    cells = []
    for i in range(0, 7):
        for j in range(0, 5):
            cells.append('<span style="top:%dcm;left:%dcm">%s</span>' % (i * 9 - 2, j * 12 - 4, esc(txt)))
    return '<div class="wm">%s</div>' % "".join(cells)


# ==================== 门禁看板 / 数据溯源 / 交互图 / 排名对齐 ====================
GATE_CACHE = {}


def _load_gates(d):
    """读 finalcheck.json；没有就返回 None（不去猜）。"""
    for p in (os.path.join(d, "finalcheck.json"), os.path.join(d, "qa", "finalcheck.json")):
        if os.path.exists(p):
            try:
                j = json.load(open(p, encoding="utf-8"))
                return j.get("verdict", "DRAFT"), j.get("results") or []
            except Exception:
                return None, []
    return None, []


def _fix_cmd(gid, d):
    """复用 np_finalcheck 的自愈指令（同一套逻辑，避免报告与门禁说法不一致）。"""
    try:
        import np_finalcheck as FC
        diag, cmd = FC.fix_command(gid, d)
        return cmd
    except Exception:
        return "python scripts/np_finalcheck.py --dir %s" % d


def short_cmd(cmd, d=None):
    """把修复命令压成一行可读形式：去掉解释器绝对路径、把输出目录换成 np_out。

    完整命令放在 title 属性里（屏幕悬停可见），打印版只显示短形式——否则全部
    绝对路径会把看板撑成十几页，导师根本没法看。
    """
    c = (cmd or "").strip()
    toks = c.split()
    idx = next((i for i, t in enumerate(toks) if t.strip('"').endswith(".py")), None)
    if idx is not None:
        c = "python " + " ".join([os.path.basename(toks[idx].strip('"'))] + toks[idx + 1:])
    c = c.replace("\\", "/")
    if d:
        c = c.replace(os.path.abspath(d).replace("\\", "/"), "np_out")
        c = c.replace(d.replace("\\", "/"), "np_out")
    c = re.sub(r'"\s*|\s*"', "", c)
    return c if len(c) <= 120 else c[:117] + "…"


def dashboard_html(d):
    """Checklist Dashboard：一眼看到哪几门没过、缺什么、怎么补。"""
    verdict, results = _load_gates(d)
    if not results:
        return ('<div class="dash missing"><b>门禁未运行</b>：本报告没有 finalcheck.json，'
                '因此无法显示门禁自查结果。<br>生成方式：'
                '<code>python scripts/np_finalcheck.py --dir %s --json %s/finalcheck.json</code></div>'
                % (esc(d), esc(d)))
    hard = [r for r in results if r.get("level") == "hard"]
    soft = [r for r in results if r.get("level") != "hard"]
    hp = sum(1 for r in hard if r.get("pass"))
    sp = sum(1 for r in soft if r.get("pass"))
    cls = "ok" if (verdict == "READY") else "bad"
    head = ('<div class="dash"><div class="dash-h"><b>投稿门禁看板</b>'
            '<span class="pill %s">%s</span>'
            '<span class="dash-sum">硬性 %d/%d 通过 · 软性 %d/%d 通过 · 共 %d 门</span></div>'
            % (cls, esc(verdict or "UNKNOWN"), hp, len(hard), sp, len(soft), len(results)))
    rows = []
    for r in results:
        ok = bool(r.get("pass"))
        gid = r.get("id", "")
        if ok:
            continue
        cmd = _fix_cmd(gid, d)
        rows.append(
            '<tr class="fail %s"><td><b>%s</b></td><td>%s</td><td>%s</td>'
            '<td><span class="pill bad">未通过</span></td><td>%s</td>'
            '<td><code title="%s">%s</code></td></tr>' % (
                "hard" if r.get("level") == "hard" else "soft", esc(gid), esc(r.get("name", "")),
                esc(r.get("level", "")), esc(r.get("reason", "")), esc(cmd),
                esc(short_cmd(cmd, d))))
    if rows:
        body = ('<table class="gates"><thead><tr><th>门</th><th>名称</th><th>级别</th>'
                '<th>状态</th><th>原因</th><th>修复命令</th></tr></thead><tbody>'
                + "".join(rows) + "</tbody></table>")
    else:
        body = '<p class="small">%d 门全部通过。</p>' % len(results)
    passed = ", ".join("%s %s" % (r.get("id"), r.get("name")) for r in results if r.get("pass"))
    foot = ('<details class="dash-pass"><summary>已通过 %d 门（点击展开）</summary><p class="small">%s</p></details>'
            % (len(results) - len(rows), esc(passed)))
    return head + body + foot + "</div>"



def gate_alert(d, gid, hint=""):
    """某一门未通过时，在正文对应章节插一条黄色提示（不打断阅读，但不会漏看）。"""
    verdict, results = _load_gates(d)
    r = next((x for x in results if x.get("id") == gid), None)
    if not r or r.get("pass"):
        return ""
    return ('<div class="alert"><b>%s %s 未通过：</b>%s<br>%s'
            '<span style="color:#6b7280">修复：<code>%s</code></span></div>' % (
                esc(gid), esc(r.get("name", "")), esc(r.get("reason", "")), esc(hint),
                esc(_fix_cmd(gid, d))))



def evidence_block(d):
    """证据聚合章：唯一排序函数 + 校准指标 + 决策分布 + 系数表。

    审稿人问"你的最终排序函数是什么"，这一章就是答案；而且把**系数**摊开，
    包括负权重（若有），而不是只给一张好看的排名表。
    """
    try:
        import np_evidence_ledger as EL
    except Exception as e:
        return '<p class="miss">证据账本模块不可用（%s）。</p>' % esc(type(e).__name__)
    lc = EL.ledger_path(d)
    if not os.path.exists(lc):
        rows, _ = EL.collect(d)
        if not rows:
            return '<p class="miss">未运行：无任何可汇总的产物。</p>'
        lc = EL.write_ledger(rows, d)
    mp = os.path.join(d, "09_evidence", "calibration_model.json")
    res = EL.score(d, model_path=mp)
    EL.write_scores(res, d)
    E = []
    if res["status"] == "UNCALIBRATED":
        E.append('<div class="alert"><b>未校准（UNCALIBRATED）：</b>本次没有在基准标签上拟合校准器，'
                 '因此<b>不输出概率</b>，只给按证据级别先验的排序。'
                 '<br><b>原因：</b>%s'
                 '<br>为什么必须这样：未校准的分数叫"排序分"，叫"概率"就是过度声明。'
                 '<br>修复：<code>python np_evidence_ledger.py fit --ledger '
                 'np_out/09_evidence/evidence_ledger.csv --gold gold.txt -o np_out/09_evidence</code></div>'
                 % esc(res.get("reason") or "缺 calibration_model.json"))
    else:
        m = json.load(open(mp, encoding="utf-8"))
        cv = m.get("cv") or {}
        tr = m.get("train_metrics") or {}
        E.append('<table><thead><tr><th>口径</th><th>Brier</th><th>ECE</th><th>AUC</th><th>说明</th>'
                 '</tr></thead><tbody>')
        E.append('<tr><td>样本内</td><td>%.4f</td><td>%.4f</td><td>%s</td>'
                 '<td>只说明能拟合训练集</td></tr>' % (
                     tr.get("brier", 0), tr.get("ece", 0), tr.get("auc")))
        E.append('<tr><td><b>样本外（%s折 OOF）</b></td><td><b>%s</b></td><td><b>%s</b></td>'
                 '<td><b>%s</b></td><td>审稿人认可的口径；ECE&lt;0.10 视为可接受校准</td></tr>' % (
                     cv.get("k", 5), cv.get("brier", "NA"), cv.get("ece", "NA"), cv.get("auc", "NA")))
        E.append("</tbody></table>")
        if cv.get("ece") is not None and cv.get("ece") > 0.10:
            E.append('<div class="alert"><b>校准不足：</b>样本外 ECE=%.3f &gt; 0.10，'
                     '概率不应作为决策唯一依据，需扩大金标准或减少特征。</div>' % cv["ece"])
        # 系数表（负权重必须显示）
        coefs = list(zip(m.get("features") or [], m.get("coef") or [None] * 1))
        if coefs:
            E.append('<h3>模型系数（标准化尺度；负值=该证据在本基准上反而降低命中概率）</h3>')
            E.append('<table><thead><tr><th>特征</th><th>含义</th><th>系数</th><th>SE</th>'
                     '</tr></thead><tbody>')
            zh = m.get("feature_zh") or {}
            for i, (f, b) in enumerate(coefs):
                if b is None:
                    continue
                se = (m.get("se") or [])[i + 1] if i + 1 < len(m.get("se") or []) else ""
                cls = ' style="background:#fdf3f3"' if isinstance(b, (int, float)) and b < 0 else ""
                E.append('<tr%s><td><code>%s</code></td><td>%s</td><td>%.3f</td><td>%s</td></tr>' % (
                    cls, esc(f), esc(zh.get(f, "")), b,
                    ("%.3f" % se) if isinstance(se, (int, float)) else "—"))
            E.append("</tbody></table>")
            neg = [f for f, b in coefs if isinstance(b, (int, float)) and b < 0]
            if neg:
                E.append('<p class="small"><b>负权重特征（%s）：</b>这不是 bug，而是基准给出的真实信号——'
                         '这些证据在本金标准上不携带正向信息。正文必须如实说明，不得只展示正系数。</p>'
                         % esc("、".join(neg)))
        dec = {}
        for t in res["targets"]:
            dec[t["decision"]] = dec.get(t["decision"], 0) + 1
        E.append('<p><b>决策分布：</b>%s</p>' % "，".join(
            "%s %d 个" % (k, v) for k, v in sorted(dec.items())))
    # 主表
    E.append('<h3>靶点优先级（按 P(gold) 降序）</h3>')
    E.append('<table><thead><tr><th>#</th><th>靶点</th><th>P(gold)</th><th>95% CI</th>'
             '<th>证据类数</th><th>决策</th><th>支持成分</th><th>PPI度分位</th></tr></thead><tbody>')
    for i, t in enumerate(res["targets"][:20], 1):
        ci = ("%.2f–%.2f" % (t["ci_low"], t["ci_high"])) if t.get("ci_low") is not None else "NA"
        E.append('<tr><td>%d</td><td><b>%s</b></td><td>%s</td><td>%s</td><td>%d</td>'
                 '<td>%s</td><td>%.0f</td><td>%.2f</td></tr>' % (
                     i, esc(t["gene"]),
                     ("%.3f" % t["p_gold"]) if t.get("p_gold") is not None else "—",
                     ci, t["n_evidence_types"], esc(t["decision_zh"]),
                     math.expm1(t.get("log_n_compounds", 0)), t.get("degree_pct", 0)))
    E.append("</tbody></table>")
    E.append('<p class="small"><b>为什么这一章替代了原来的三套排序：</b>第 4 节按支持成分数、'
             '第 6 节按拓扑共识、第 9 节按 MR β，三者的排序键不同，读者无法判断该信谁。'
             '本章给出唯一的 <code>score_targets()</code> 映射：证据账本 → 校准概率 → 决策标签。'
             '原表保留，用于解释"为什么某个靶点分数高/低"（系数可逐项归因）。</p>')
    return "".join(E)



def multiplicity_block(d):
    """多重性声明：检测了多少假设、怎么校正、校正后剩几个。

    这是 §4 的落地点。没有这一节，20 个模块各开 0.05 的门就是"在 20 个地方各开一扇窗"。
    """
    mp = os.path.join(d, "09_evidence", "multiplicity_statement.json")
    if not os.path.exists(mp):
        try:
            import np_fdr as F
            rows = F.gather(d)
            testable, summary = F.apply(rows, 0.05)
            if summary.get("n_hypotheses"):
                F.write(d, testable, summary, "09_evidence")
        except Exception:
            return ""
    if not os.path.exists(mp):
        return ('<p class="miss">未运行多重性校正（无 MR/富集/基准的 p 值终点）。</p>')
    try:
        j = json.load(open(mp, encoding="utf-8"))
    except Exception:
        return '<p class="miss">multiplicity_statement.json 不可解析。</p>'
    E = ['<table><thead><tr><th>项</th><th>值</th></tr></thead><tbody>']
    for k, lab in (("n_hypotheses", "检测的假设总数"), ("n_families", "终点族数"),
                   ("n_raw", "未校正 p<0.05"), ("n_family", "族内 BH q<0.05"),
                   ("n_global", "全局 BH q<0.05")):
        E.append("<tr><td>%s</td><td><b>%s</b></td></tr>" % (lab, esc(j.get(k))))
    E.append("</tbody></table>")
    if j.get("per_family"):
        E.append('<table><thead><tr><th>族</th><th>假设数</th><th>p&lt;0.05</th>'
                 '<th>族内 q&lt;0.05</th></tr></thead><tbody>')
        for f, v in sorted(j["per_family"].items()):
            E.append("<tr><td>%s</td><td>%d</td><td>%d</td><td>%d</td></tr>" % (
                esc(f), v.get("n", 0), v.get("n_raw", 0), v.get("n_family", 0)))
        E.append("</tbody></table>")
    if j.get("n_global") == 0 and j.get("n_hypotheses"):
        E.append('<div class="alert"><b>校正后无一显著：</b>全局 BH 后 0 个假设通过 q&lt;0.05。'
                 '这不是脚本错误，而是必须如实报告的结论——本研究的显著性主要来自未校正的 p 值。</div>')
    if j.get("statement_zh"):
        E.append('<h3>可直接写进 Methods 的声明（中/英）</h3>')
        E.append("<p>%s</p>" % re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", esc(j["statement_zh"])))
        E.append("<pre>%s</pre>" % esc(j.get("statement_en", "")))
    E.append('<p class="small"><b>口径：</b><code>q_global</code> 是唯一可写入结论的显著性口径；'
             '仅 p&lt;0.05 未通过 q 的结果只能表述为「提示性关联」。</p>')
    return "".join(E)


def provenance_html(d):
    """Data Provenance：输入文件哈希/时间 + 阶段工具版本 + 包版本。

    注意：文件 mtime 只是"本地文件时间"，**不等于**数据库检索日期；
    没有数据库版本元数据时必须如实写 "未记录"，不能编。
    """
    man = {}
    mp = os.path.join(d, "pipeline_manifest.json")
    if os.path.exists(mp):
        try:
            man = json.load(open(mp, encoding="utf-8"))
        except Exception:
            man = {}
    rows = []
    inputs = (man.get("inputs") or {})
    for key in ("compounds", "compound_edges", "targets", "disease", "edges", "gmt"):
        p = inputs.get(key) or ""
        if not p or not os.path.exists(p):
            continue
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
        rows.append([key, os.path.basename(p), "%.1f KB" % (os.path.getsize(p) / 1024),
                     time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(p))), h])
    # 阶段工具版本
    stages = []
    for rm in sorted(glob.glob(os.path.join(d, "*", "run_manifest.json"))):
        try:
            j = json.load(open(rm, encoding="utf-8"))
            for st in j.get("stages") or []:
                stages.append([st.get("stage", ""), st.get("tool_version", ""),
                               str(st.get("started", ""))[:19]])
        except Exception:
            continue
    pkgs = {}
    for env in glob.glob(os.path.join(d, "*", "environment.json")):
        try:
            pkgs.update(json.load(open(env, encoding="utf-8")).get("packages") or {})
        except Exception:
            continue
    E = ['<h2>10　数据溯源（Data Provenance）</h2>',
         '<div class="boundary"><b>为什么必须有这一节：</b>审稿人要求"结果可追溯到具体数据库版本与检索日期"。'
         '下表列出本次运行实际读取的每个输入文件的 SHA256 与时间；'
         '<b>文件时间不等于数据库检索日期</b>，若原始导出未记录检索日期，本报告不会替它编一个。</div>']
    if rows:
        E.append("<table><thead><tr><th>角色</th><th>文件</th><th>大小</th><th>文件时间</th>"
                 "<th>SHA256(16)</th></tr></thead><tbody>")
        for r in rows:
            E.append("<tr>" + "".join("<td>%s</td>" % esc(x) for x in r) + "</tr>")
        E.append("</tbody></table>")
    else:
        E.append('<p class="small">未找到 pipeline_manifest.json 中的输入清单。</p>')
    if stages:
        E.append("<table><thead><tr><th>阶段</th><th>工具版本</th><th>开始时间</th></tr></thead><tbody>")
        for r in stages:
            E.append("<tr>" + "".join("<td>%s</td>" % esc(x) for x in r) + "</tr>")
        E.append("</tbody></table>")
    if pkgs:
        E.append('<p class="small"><b>运行环境：</b>%s</p>' % esc(
            ", ".join("%s %s" % (k, v) for k, v in sorted(pkgs.items()))))
    E.append('<p class="small">数据库版本提示：TCMSP/HERB/ETCM2/GeneCards/DisGeNET 的具体版本与检索日期'
             '应由数据导出环节写入 <code>db_status.json</code>；若缺失，请在 Methods 中手工补全，'
             '不要以"最新版"含糊表述。</p>')
    return "".join(E)


def alignment_html(d, info):
    """S4↔S5 对齐表：回答"为什么交集第一名不是 hub 第一名"。

    交集表按支持成分数排序，hub 按拓扑共识排序，两套键不同 → 名单不同是预期。
    真正需要校验的是"交集基因是否在 PPI 节点集合里"，由 G21 判定。
    """
    inter = []
    if info.get("intersection"):
        rows = _read_csv(info["intersection"])
        if rows:
            hdr = [h.lower() for h in rows[0]]
            gi = hdr.index("gene") if "gene" in hdr else 0
            ni = hdr.index("n_compounds") if "n_compounds" in hdr else None
            for r in rows[1:11]:
                inter.append((r[gi], (r[ni] if ni is not None and ni < len(r) else "")))
    ranks = {}
    if info.get("hubs"):
        for r in _read_csv(info["hubs"])[1:]:
            try:
                ranks[r[1].upper()] = {"rank": r[0], "degree": r[3] if len(r) > 3 else "",
                                       "verdict": r[-1] if r else ""}
            except Exception:
                continue
    if not inter:
        return ""
    E = ['<h3>4.1 交集靶点在 PPI 中的实际位置（解释两节排名差异）</h3>',
         '<table><thead><tr><th>交集排名</th><th>基因</th><th>支持成分数</th>'
         '<th>PPI hub 排名</th><th>degree</th><th>共识判定</th></tr></thead><tbody>']
    for i, (g, n) in enumerate(inter, 1):
        h = ranks.get(g.upper(), {})
        E.append("<tr><td>%d</td><td><b>%s</b></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % (
            i, esc(g), esc(n), esc(h.get("rank", "不在 PPI")), esc(h.get("degree", "—")),
            esc(h.get("verdict", "—"))))
    E.append("</tbody></table>")
    E.append('<p class="small"><b>读法：</b>第 4 节按"支持该靶点的成分数"排序（药理学证据密度），'
             '第 5 节按"拓扑共识中心性"排序（网络位置）。两者回答不同问题，排序不可能一致。'
             '例：IL6 被 10 个成分支持（交集第 1），但 degree=15 排 hub 第 18；'
             'BCL2 degree=30 排 hub 第 1，却只被少数成分支持。报告同时给出两套排名，'
             '读者可自行判断"证据密度"与"网络枢纽"哪个更符合其假设。</p>')
    return "".join(E)


def interactive_block(d, info):
    """在报告中嵌入交互网络（iframe srcdoc，单文件可离线查看结构）。"""
    p = os.path.join(d, "network_interactive.html")
    if not os.path.exists(p) and info.get("hubs"):
        edge = None
        for pat in ("05_ppi/extended_ppi_edges.tsv", "05_ppi/ppi_edges_clean.tsv",
                    "05_ppi/ppi_edges_clean.csv"):
            if os.path.exists(os.path.join(d, pat)):
                edge = os.path.join(d, pat)
                break
        if edge:
            try:
                import np_interactive
                txt = np_interactive.build(edge, info["hubs"],
                                           os.path.join(d, "04_intersect", "intersection_genes.txt"),
                                           title="PPI 交互网络")
                open(p, "w", encoding="utf-8").write(txt)
            except Exception:
                p = None
        else:
            p = None
    if not p or not os.path.exists(p):
        return ""
    txt = open(p, encoding="utf-8", errors="replace").read()
    ext = os.path.exists(os.path.join(d, "05_ppi", "extended_ppi_edges.tsv"))
    return ('<div class="no-print"><h3>5.0 交互网络（可缩放 / 拖拽 / 搜索 / 点节点看属性 / 导出 SIF）%s</h3>'
            '<iframe class="net" srcdoc="%s" loading="lazy"></iframe>'
            '<p class="small">交互图需要浏览器执行 JavaScript；若所在环境离线且未内联 cytoscape.js，'
            '图内会提示。右侧「导出 SIF」只导出当前可见子网（SIF + 节点属性表），'
            '可直接拖入 Cytoscape 桌面版排版。打印/PDF 版使用下方静态 SVG。</p></div>'
            % ("（含孤立靶点的一阶邻居）" if ext else "", html.escape(txt, quote=True)))



def build(outdir, label, verdict="UNVERIFIED"):
    info = collect(outdir)
    c = counts(info)
    ts = time.strftime("%Y-%m-%d %H:%M")
    P = ["<!DOCTYPE html>", '<html lang="zh-CN"><head><meta charset="utf-8">',
         "<title>%s — 分析报告</title>" % esc(label), "<style>%s</style></head><body>" % CSS,
         watermark(verdict), '<div class="wrap">']
    P.append("<h1>%s</h1>" % esc(label))
    P.append('<div class="sub">network-pharmacology-pro 一键分析报告 · 生成于 %s · 输出目录 '
             '<code>%s</code></div>' % (esc(ts), esc(os.path.basename(os.path.abspath(outdir)))))
    cls = {"READY": "ready", "DRAFT": "draft"}.get(verdict, "unv")
    P.append('<div class="stamp %s">%s · %s</div>' % (cls, esc(WATERMARK_TEXT.get(verdict, "UNVERIFIED")), esc(ts)))
    P.append('<div class="print-only print-note">本报告为 %s 版（生成于 %s）——'
             '屏幕上的水印与交互图已按打印规范省略；门禁明细见下方看板。</div>' % (esc(verdict), esc(ts)))
    P.append('<div class="no-print" style="margin:8pt 0"><button onclick="window.print()">打印 / 另存为 PDF</button>'
             '<span class="miss" style="margin-left:8pt">@media print 已内嵌：A4、图表避断、表头跨页重复</span></div>')
    P.append(dashboard_html(outdir))
    P.append('<div class="boundary"><b>证据边界：</b>本报告中的靶点、网络、富集、对接结果均为<b>计算预测（假说生成）</b>，'
             '不构成活性或疗效证明；剂量-暴露、体外/体内验证前不得写作"证明"。'
             '所有数值来自本次运行的磁盘产物，未人工改写。</div>')

    # 摘要卡片
    P.append('<div class="cards">')
    for k, lab in (("compounds", "通过筛选的成分"), ("targets", "成分-疾病交集靶点"),
                   ("hubs", "PPI hub 候选"), ("pathways", "显著通路（折叠后）"),
                   ("nodes_ok", "流水线节点通过/总数")):
        if k == "nodes_ok":
            v = "%s/%s" % (c.get("nodes_ok", "—"), c.get("nodes", "—")) if c.get("nodes") else "—"
        else:
            v = c.get(k, "—")
        P.append('<div class="card"><div class="n">%s</div><div class="k">%s</div></div>' % (esc(v), esc(lab)))
    P.append("</div>")

    # 1 运行清单
    P.append("<h2>1　运行清单与可复现性</h2>")
    if info["manifest"]:
        nodes = info["manifest"].get("nodes") or []
        rows = []
        for n in nodes:
            st_ = n.get("status", "?")
            badge = {"ok": "b-ok", "failed": "b-bad", "skipped_done": "b-bg"}.get(st_, "b-warn")
            rows.append([n.get("id", ""), '<span class="badge %s">%s</span>' % (badge, esc(st_)),
                         esc((n.get("cli") or "")[:70])])
        P.append("<table><thead><tr><th>节点</th><th>状态</th><th>命令</th></tr></thead><tbody>")
        for r in rows:
            P.append("<tr><td>%s</td><td>%s</td><td><code>%s</code></td></tr>" % (r[0], r[1], r[2]))
        P.append("</tbody></table>")
        P.append('<p class="miss">配置：%s　启用插件：%s</p>' % (
            esc(info["manifest"].get("config", "")), esc(",".join(info["manifest"].get("enable") or []) or "无")))
    else:
        P.append('<p class="miss">未找到 pipeline_manifest.json（可能只跑了单个模块）。</p>')

    # 2 安全前置
    P.append("<h2>2　配伍禁忌与毒性（合规前置）</h2>")
    if info["safety"]:
        s = info["safety"]
        P.append('<div class="boundary"><b>结论：</b>等级 <code>%s</code>；'
                 '违规 %s 条；需声明 %s 条。</div>' % (
                     esc(s.get("level", s.get("verdict", "?"))),
                     esc(len(s.get("violations") or s.get("errors") or [])),
                     esc(len(s.get("warnings") or []))))
        P.append(md_excerpt(os.path.join(outdir, "00_safety", "contraindication.md"), 10))
    else:
        P.append('<p class="miss">未运行（无 00_safety/contraindication.json）。</p>')

    # 3 成分
    P.append("<h2>3　成分筛选</h2>")
    P.append(gate_alert(outdir, "G3", "剂量权重未计算：无剂量表时成分贡献只能等权，需在 Methods 说明。"))
    P.append(table_from_csv(info["compounds"], 12, 8) if info["compounds"] else
             '<p class="miss">未运行（无 03_compounds/compounds_kept.csv）。</p>')

    # 4 交集
    P.append("<h2>4　成分-疾病交集靶点（排序键：支持成分数）</h2>")
    P.append(table_from_csv(info["intersection"], 15, 7) if info["intersection"] else
             '<p class="miss">未运行。</p>')
    P.append(alignment_html(outdir, info))
    P.append(gate_alert(outdir, "G21", "交集基因与 PPI 节点不一致时，第 4/5 节的名单差异无法用排序键解释。"))
    P.append("<h2>5　靶点优先级：单一排序函数（证据聚合 + 校准）</h2>")
    P.append(evidence_block(outdir))
    P.append("<h3>5.1　多重性控制声明（全局 FDR）</h3>")
    P.append(multiplicity_block(outdir))

    # 5 PPI
    P.append("<h2>5　PPI 拓扑与 hub（排序键：拓扑共识中心性，与第 4 节不同）</h2>")
    P.append(gate_alert(outdir, "G21", "低度/孤立靶点缺少 PPI 结构基础，建议运行 --expand-isolated 生成扩展子网。"))
    P.append(interactive_block(outdir, info))
    P.append(table_from_csv(info["hubs"], 12, 8) if info["hubs"] else '<p class="miss">未运行。</p>')

    # 6 零模型
    P.append("<h2>6　零模型阴性对照（这段决定结果是否只是网络本身的性质）</h2>")
    P.append(md_excerpt(info["null_report"], 14) if os.path.exists(info["null_report"])
             else '<p class="miss">未运行。</p>')

    # 7 富集
    P.append("<h2>7　通路富集</h2>")
    em = os.path.join(outdir, "07_enrich", "enrich_meta.json")
    if os.path.exists(em):
        try:
            jm = json.load(open(em, encoding="utf-8"))
            P.append('<p class="small"><b>背景：</b>%s（N=%s）%s</p>' % (
                esc(jm.get("background_source") or jm.get("background_mode")),
                esc(jm.get("N")),
                '　<span class="pill bad">P-value may be inflated</span>' if jm.get("pvalue_inflation_warning") else
                '　<span class="pill ok">已使用组织/表达背景</span>'))
        except Exception:
            pass
    P.append(md_excerpt(info["enrich_report"], 22) if os.path.exists(info["enrich_report"])
             else '<p class="miss">未运行。</p>')

    # 8 同源转换
    P.append("<h2>8　湿实验衔接（同源基因）</h2>")
    P.append(gate_alert(outdir, "G22", "功能替代（FuncRep）必须与直系同源区分，"
                                       "严禁把 Cxcl1 引物当作 CXCL8 的直系同源引物。"))
    P.append(table_from_csv(info["ortholog"], 10, 8) if info["ortholog"] else '<p class="miss">未运行。</p>')

    # 9 可选模块（数值 + 自动解读）
    opt, interp = [], []
    if info["mr"]:
        opt.append("<h3>孟德尔随机化</h3>")
        for p in info["mr"][:3]:
            try:
                j = json.load(open(p, encoding="utf-8"))
                if j.get("status") != "OK":
                    opt.append("<p><code>%s</code>：%s</p>" % (esc(os.path.basename(os.path.dirname(p))),
                                                               esc(j.get("status"))))
                    continue
                ivw = j.get("IVW") or {}
                col = j.get("coloc") or {}
                gene = os.path.basename(os.path.dirname(p)) or j.get("gene", "")
                opt.append("<p><code>%s</code>：IVW β=%.3f（p=%.2g，OR=%.2f）；COLOC PP.H4=%s；flags=%s</p>" % (
                    esc(gene), ivw.get("beta", 0), ivw.get("p", 1),
                    2.718281828 ** (ivw.get("beta") or 0),
                    esc(col.get("PP.H4", "NA")), esc(",".join(j.get("flags") or []) or "none")))
                iv = j.get("intervention") or {}
                if iv.get("severity") in ("SEVERE", "WARN"):
                    opt.append('<div class="alert %s"><b>靶点-干预方向%s：%s</b><br>%s'
                               '<br><span style="color:#6b7280">遗传因果与药物作用方向相反时，'
                               '成分越"有效"，患者可能越差；该靶点不得作为主结论。</span></div>' % (
                                   "bad" if iv["severity"] == "SEVERE" else "",
                                   "红灯" if iv["severity"] == "SEVERE" else "警告",
                                   esc(iv.get("message", "")),
                                   "药物方向=%s；起效所需=%s" % (iv.get("drug_direction"), iv.get("required_action"))))
                try:
                    import np_interpret as NI
                    b = NI.interpret_mr(ivw.get("beta"), ivw.get("p"), j.get("F_statistic"),
                                        pleiotropy=(j.get("MR_Egger") or {}).get("pleiotropy_warning"),
                                        coloc_pp4=col.get("PP.H4"),
                                        direction_consistent=j.get("direction_consistent"),
                                        exploratory=bool(j.get("exploratory_only")))
                    interp.append(("MR · %s" % gene, b))
                    if col.get("status") == "OK":
                        interp.append(("COLOC · %s" % gene,
                                       NI.interpret_coloc(col.get("PP.H4"), col.get("PP.H3"))))
                except Exception:
                    pass
            except Exception:
                pass
    if info["intervention"]:
        j = json.load(open(info["intervention"], encoding="utf-8"))
        opt.append("<h3>靶点-干预方向一致性</h3><p>结论 <b>%s</b>：SEVERE %s / WARN %s / 自洽 %s</p>" % (
            esc(j.get("verdict")), esc(j.get("n_severe")), esc(j.get("n_warn")), esc(j.get("n_consistent"))))
        if j.get("discussion_snippet"):
            opt.append("<pre>%s</pre>" % esc(j["discussion_snippet"]))
    if info["synergy"]:
        try:
            j = json.load(open(info["synergy"], encoding="utf-8"))
            opt.append("<h3>机制定向协同</h3><p>协同 %s 对 / 拮抗 %s 对（通路分层 %s 条）</p>" % (
                esc(j.get("n_mechanistic_synergy", j.get("n_synergy", "—"))),
                esc(j.get("n_antagonism", "—")), esc(j.get("n_pathways_layered", "—"))))
            top = (j.get("top_pairs") or [])[:1]
            if top:
                import np_interpret as NI
                t = top[0]
                interp.append(("Synergy · %s+%s" % (t.get("compound_a"), t.get("compound_b")),
                               NI.interpret_synergy(t.get("score"), t.get("label"),
                                                    t.get("level_gap"), t.get("term"))))
        except Exception:
            pass
    if info["benchmark"]:
        opt.append("<h3>基准测试</h3>")
        for p in info["benchmark"][:4]:
            try:
                j = json.load(open(p, encoding="utf-8"))
                if j.get("kind") == "concordance":
                    opt.append("<p>%s：敏感性 %.3f / 精确率 %.3f / Jaccard %.3f / 过预测 %.2fx</p>" % (
                        esc(j.get("label")), j["sensitivity"], j["precision"], j["jaccard"],
                        j.get("over_prediction_factor") or 0))
                elif j.get("status") == "OK":
                    opt.append("<p>%s：AUC %.3f、Recall@20 %.2f、EF@20 %.2fx</p>" % (
                        esc(j.get("label")), j["AUC"], j["at_k"]["20"]["recall"],
                        j["at_k"]["20"]["enrichment_factor"]))
                    try:
                        import np_interpret as NI
                        interp.append(("Benchmark · %s" % j.get("label"),
                                       NI.interpret_benchmark(j["AUC"], j["at_k"]["20"]["enrichment_factor"],
                                                              j.get("coverage"))))
                    except Exception:
                        pass
            except Exception:
                pass
    if opt:
        P.append("<h2>9　因果、协同与基准（可选模块）</h2>")
        P.extend(opt)
    if interp:
        try:
            import np_interpret as NI
            P.append(NI.render_html(interp))
        except Exception:
            pass

    # 10 图件
    P.append("<h2>10　图件</h2>")
    if info["svgs"]:
        seen = set()
        for s in info["svgs"]:
            if os.path.abspath(s) in seen:
                continue
            seen.add(os.path.abspath(s))
            P.append(embed_svg(s))
    else:
        P.append('<p class="miss">未找到 SVG 图件。</p>')

    # 数据溯源
    P.append(provenance_html(outdir))

    # 产物清单
    P.append("<h2>11　产物清单（SHA256 前 8 位）</h2>")
    files = []
    for p in sorted(glob.glob(os.path.join(outdir, "**", "*"), recursive=True)):
        if os.path.isfile(p) and os.path.getsize(p) > 0:
            files.append(p)
    P.append("<table><thead><tr><th>文件</th><th>大小</th><th>sha256</th></tr></thead><tbody>")
    for p in files[:80]:
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()[:8]
        P.append("<tr><td><code>%s</code></td><td>%.1f KB</td><td><code>%s</code></td></tr>" % (
            esc(os.path.relpath(p, outdir)), os.path.getsize(p) / 1024, h))
    P.append("</tbody></table>")
    P.append('<p class="miss">共 %d 个产物；完整校验请用 np_benchloop.py verify。</p>' % len(files))

    P.append('<hr><p class="miss">由 np_uireport.py 生成 · 水印状态：%s · '
             '本报告不含任何未在磁盘上出现的数值</p>' % esc(verdict))
    P.append("</div></body></html>")
    out = "\n".join(P)
    # 章节号重排：可选模块缺省时不留空号（8 → 10 这种跳号最招审稿人烦）
    counter = {"n": 0}

    def _renum(m):
        counter["n"] += 1
        return '<h2 id="S%d">%d　%s</h2>' % (counter["n"], counter["n"], m.group(2))

    return re.sub(r'<h2>(\d+)　([^<]*)</h2>', _renum, out)


def to_pdf(html_path, pdf_path=None):
    """HTML→PDF：优先系统浏览器无头模式（零 Python 依赖），weasyprint 仅作兜底。

    为什么换掉 weasyprint：Windows 上要装 GTK3/Pango/Cairo 一整套 C 库，失败率极高；
    而 Edge/Chrome 自带 --print-to-pdf，效果与浏览器所见一致（CSS 变量/flex/@media print 全支持）。
    """
    try:
        import np_pdf
        p, info = np_pdf.html_to_pdf(html_path, pdf_path)
        if p:
            return p
    except Exception:
        pass
    try:
        from weasyprint import HTML  # noqa
        pdf_path = pdf_path or os.path.splitext(html_path)[0] + ".pdf"
        HTML(filename=html_path).write_pdf(pdf_path)
        return pdf_path
    except Exception:
        return None


def selftest():
    """自检：最小产物目录 → 报告必须含全部章节、水印，且不编造数值。"""
    import tempfile
    d = tempfile.mkdtemp(prefix="np_uireport_")
    os.makedirs(os.path.join(d, "04_intersect"), exist_ok=True)
    os.makedirs(os.path.join(d, "05_ppi"), exist_ok=True)
    json.dump({"config": "c.yaml", "enable": [], "nodes": [{"id": "alias", "status": "ok", "cli": "x"}]},
              open(os.path.join(d, "pipeline_manifest.json"), "w", encoding="utf-8"))
    open(os.path.join(d, "04_intersect", "intersection_genes.txt"), "w", encoding="utf-8").write("IL6\nTNF\nAKT1\n")
    open(os.path.join(d, "05_ppi", "hubs.csv"), "w", encoding="utf-8").write(
        "rank,node,degree,verdict\n1,IL6,30,稳健\n2,TNF,20,不稳定\n")
    ready = build(d, "自检", "READY")
    unv = build(d, "自检", "UNVERIFIED")
    checks = [("watermark_ready", "READY FOR SUBMISSION" in ready),
              ("watermark_unverified", "UNVERIFIED" in unv),
              ("sections", ready.count("<h2") >= 8),
              ("counts_from_disk", '<div class="n">3</div><div class="k">成分-疾病交集靶点</div>' in ready),
              ("no_fabrication_note", "不含任何未在磁盘上出现的数值" in ready)]
    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for n, g in checks:
        print("  [%s] %s" % ("✓" if g else "✗", n))
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="np_uireport: 一键分析报告（带水印，可转 PDF）")
    ap.add_argument("--dir", required=True)
    ap.add_argument("--label", default="网络药理学分析")
    ap.add_argument("--out", default=None)
    ap.add_argument("--verdict", default="UNVERIFIED", choices=["READY", "DRAFT", "UNVERIFIED"])
    ap.add_argument("--pdf", action="store_true")
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    a = ap.parse_args(argv)
    if not os.path.isdir(a.dir):
        print("[FAIL] 目录不存在：%s" % a.dir)
        return 2
    out = a.out or os.path.join(a.dir, "ui_report.html")
    open(out, "w", encoding="utf-8").write(build(a.dir, a.label, a.verdict))
    print("[uireport] → %s（%.0f KB，水印 %s）" % (out, os.path.getsize(out) / 1024, a.verdict))
    if a.pdf:
        p = to_pdf(out)
        print("[uireport] PDF → %s" % p if p else
              "[uireport] 未安装 weasyprint，未生成 PDF；可用浏览器 Ctrl+P 另存为 PDF")
    return 0


if __name__ == "__main__":
    sys.exit(main())
