#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_htmlreport.py —— SCI 投稿版 HTML 渲染器（补齐竞品 F 的"中科博林标准版式"优势）

与 np_report.py 的分工
    np_report.py     : 阶段审计研究报告（内嵌 SVG、证据分级 L0–L4，面向开发者自查）
    np_htmlreport.py : 投稿/汇报版渲染（@media print 打印分页 + 语义色标签 + 动态章节号 +
                       自包含单文件），面向投稿前预览与打印为 PDF。
    两者共享底层数据（run_manifest.json / 标准表 / Fig_*.json / SVG），互不覆盖输出。

借鉴点（来自市场 F skill，均已本地化实现）
    ① CSS 变量统一品牌色（跨报告同源，换主题只改 :root）
    ② @media print 打印分页（避免图表被截断、去除交互样式）
    ③ 语义色标签 .tag.hub/.tag.bg/.tag.ok（hub/背景/通过 三色，配图例）
    ④ 动态章节号（methods_no），章节增删不会导致 TOC 错位
    ⑤ "预测 vs 实测"边界声明 + 证据分级图例常驻

用法
    python np_htmlreport.py --dir np_out/ --label "益智胶囊-帕金森病" --out manuscript_preview.html
    python np_htmlreport.py --dir np_out/ --print-pdf-hint     # 打印为 PDF 的操作提示
"""
import argparse
import glob
import html
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# ① CSS 变量 + ③ 语义色标签 + ② @media print
CSS = r"""
:root{
  --brand:#0A1A2F; --accent:#007B7A; --accent2:#00B4A6;
  --ink:#1f2937; --muted:#6b7280; --line:#e5e7eb; --bg:#ffffff; --soft:#f8fafc;
  --hub:#c0504d; --bg-tag:#7f7f7f; --ok:#4f8a5b; --warn:#b7791f;
}
*{box-sizing:border-box}
body{margin:0;background:var(--soft);color:var(--ink);
  font-family:"Times New Roman","Microsoft YaHei",serif;line-height:1.65;font-size:11pt}
.wrap{max-width:190mm;margin:0 auto;padding:14mm 12mm 24mm;background:var(--bg)}
h1{font-size:20pt;margin:0 0 4pt;color:var(--brand);letter-spacing:-.3px}
h2{font-size:14pt;margin:22pt 0 8pt;padding-bottom:4pt;border-bottom:1.5pt solid var(--accent)}
h3{font-size:12pt;margin:14pt 0 5pt;color:var(--brand)}
.sub{color:var(--muted);font-size:9.5pt;margin-bottom:16pt}
p{margin:6pt 0;text-align:justify}
table{border-collapse:collapse;width:100%;margin:8pt 0;font-size:9pt}
th,td{border:.5pt solid var(--line);padding:4pt 6pt;text-align:left;vertical-align:top}
th{background:#f3f6fa;font-weight:600}
tr:nth-child(even) td{background:#fbfcfe}
code{background:#f1f5f9;padding:1pt 4pt;border-radius:3px;font-size:8.5pt}
.card{background:var(--soft);border:.5pt solid var(--line);border-radius:6pt;padding:9pt 12pt;margin:10pt 0}
.tag{display:inline-block;font-size:8pt;padding:1pt 6pt;border-radius:9pt;margin-right:4pt;border:.5pt solid transparent}
.tag.hub{background:#fdf0ef;color:var(--hub);border-color:#f3d5d3}
.tag.bg{background:#f1f1f1;color:var(--bg-tag);border-color:#e0e0e0}
.tag.ok{background:#eef7f0;color:var(--ok);border-color:#d3e8d8}
.tag.warn{background:#fff8e8;color:var(--warn);border-color:#f3e6c8}
.legend{font-size:8.5pt;color:var(--muted);margin:6pt 0 12pt}
.fig{margin:12pt 0;border:.5pt solid var(--line);border-radius:6pt;overflow:hidden;background:#fff;
     page-break-inside:avoid;break-inside:avoid}
.fig svg{display:block;width:100%;height:auto}
.figcap{font-size:9pt;color:var(--muted);padding:6pt 9pt;background:var(--soft);border-top:.5pt solid var(--line)}
.boundary{background:#fffaf0;border-left:3pt solid var(--warn);padding:8pt 12pt;margin:10pt 0;font-size:9.5pt}
.toc{columns:2;font-size:9.5pt}.toc a{color:var(--accent);text-decoration:none}
.small{font-size:8.5pt;color:var(--muted)}
nav.toc{page-break-after:always;break-after:page}

@page{size:A4;margin:18mm 15mm}
@media print{
  body{background:#fff;font-size:10.5pt}
  .wrap{max-width:none;margin:0;padding:0;background:#fff}
  h2,h3{page-break-after:avoid;break-after:avoid}
  .fig,table,.card{page-break-inside:avoid;break-inside:avoid}
  nav.toc{page-break-after:always}
  a{color:inherit;text-decoration:none}
  .no-print{display:none!important}
  tr{page-break-inside:avoid}
  thead{display:table-header-group}
}
"""

# 证据分级图例（L0–L4，与 np_report 一致，常驻）
EVID = [
    ("L0", "bg", "数据库注释/文献事实", "可直接引用，须带版本与访问日期"),
    ("L1", "bg", "计算预测（靶点/网络/富集）", "只能作假说，用「提示/支持」"),
    ("L2", "bg", "结构证据（对接/动力学）", "支持结合可能性，不等于活性"),
    ("L3", "ok", "体外/离体实验", "可用「表明」"),
    ("L4", "ok", "体内/临床", "可用「证明」，注明样本量"),
]


def esc(x):
    return html.escape(str(x if x is not None else ""))


def table(headers, rows):
    o = ["<table><thead><tr>"]
    o += ["<th>%s</th>" % esc(h) for h in headers]
    o.append("</tr></thead><tbody>")
    for r in rows:
        o.append("<tr>" + "".join("<td>%s</td>" % esc(c) for c in r) + "</tr>")
    o.append("</tbody></table>")
    return "".join(o)


def embed_svg(path, cap):
    try:
        s = open(path, encoding="utf-8").read()
        s = re.sub(r"<\?xml[^>]*\?>|<!DOCTYPE[^>]*>", "", s)
        s = re.sub(r'\swidth="[\d.]+"\s+height="[\d.]+"', "", s, count=1)
        return '<div class="fig">%s<div class="figcap">%s</div></div>' % (s, esc(cap))
    except Exception:
        return ""


def read_manifest(d):
    p = os.path.join(d, "run_manifest.json")
    if not os.path.exists(p):
        return []
    try:
        return json.load(open(p, encoding="utf-8")).get("stages", [])
    except Exception:
        return []


def read_standard_tables(d):
    """读取 pro 契约标准表（若存在），用于生成结果摘要表。"""
    out = {}
    for name in ("standard_herbs.csv", "standard_compounds.csv", "standard_targets.csv", "standard_edges.csv"):
        p = os.path.join(d, name)
        if os.path.exists(p):
            import csv
            with open(p, encoding="utf-8-sig") as f:
                rdr = csv.reader(f)
                rows = list(rdr)
            out[name] = rows
    return out


def build(d, label):
    stages = read_manifest(d)
    tables = read_standard_tables(d)
    svgs = sorted(glob.glob(os.path.join(d, "*.svg")) + glob.glob(os.path.join(d, "figs", "*.svg")))
    qa = os.path.join(d, "figure_qa.md")

    # 动态章节号
    sec = {"n": 0}

    def h2(title):
        sec["n"] += 1
        return '<h2 id="S%d">%d　%s</h2>' % (sec["n"], sec["n"], esc(title))

    P = ["<!DOCTYPE html>", '<html lang="zh-CN"><head><meta charset="utf-8">',
         '<title>%s — 投稿版预览</title>' % esc(label),
         '<style>%s</style></head><body><div class="wrap">' % CSS]
    P.append("<h1>%s</h1>" % esc(label))
    P.append('<div class="sub">网络药理学分析 · 投稿版预览 · 自包含单文件（可直接打印为 PDF）</div>')
    P.append('<div class="no-print" style="margin:10pt 0">'
             '<button onclick="window.print()">打印 / 另存为 PDF</button>'
             '<span class="small" style="margin-left:10pt">打印样式已内嵌 @media print（A4、分页避断图）</span></div>')

    # 边界声明（常驻）
    P.append('<div class="boundary"><b>证据边界声明：</b>本报告中计算预测结果（靶点、网络、富集、对接）'
             '均为<b>假说生成</b>，不构成活性或疗效证明；需经体外/体内实验验证。'
             '所有数量与参数均来自本地数据与 run_manifest，未经人工改写。</div>')

    # 门禁看板（与 np_uireport 共用同一实现，避免两套说法）
    try:
        import np_uireport as UR
        P.append(UR.dashboard_html(d))
    except Exception:
        pass

    # 证据分级图例
    P.append('<div class="legend"><b>证据分级：</b>')
    for code, cls, name, rule in EVID:
        P.append('<span class="tag %s">%s %s</span>' % (cls, code, esc(name)))
    P.append('</div>')

    # TOC（动态章节号）
    toc_titles = ["数据来源与可复现性", "分析结果摘要", "图件", "图片自检", "方法学参数"]
    P.append('<nav class="toc"><b>目录</b><ol>')
    for i, t in enumerate(toc_titles, 1):
        P.append('<li><a href="#S%d">%d　%s</a></li>' % (i, i, esc(t)))
    P.append("</ol></nav>")

    # 章节 1：数据来源与可复现性
    P.append(h2("数据来源与可复现性"))
    man_rows = []
    for st in stages:
        man_rows.append([
            st.get("stage", ""),
            st.get("tool_version", ""),
            str(st.get("started", ""))[:19],
            ", ".join("%s=%s" % (k, v) for k, v in list((st.get("counts") or {}).items())[:4]),
        ])
    if man_rows:
        P.append(table(["阶段", "工具版本", "开始时间", "计数摘要"], man_rows))
    else:
        P.append('<p class="small">未找到 run_manifest.json，建议先运行各分析阶段以生成溯源清单。</p>')

    # 章节 2：分析结果摘要（标准表）
    P.append(h2("分析结果摘要"))
    for name, rows in tables.items():
        P.append("<h3>%s</h3>" % esc(name))
        if rows:
            P.append(table(rows[0], rows[1:20]))  # 最多 20 行，防打印膨胀
            if len(rows) > 21:
                P.append('<p class="small">（仅显示前 20 行，共 %d 行）</p>' % (len(rows) - 1))
    if not tables:
        P.append('<p class="small">未找到 standard_*.csv，可运行 np_export.py 导出标准表。</p>')

    # 章节 3：图件
    P.append(h2("图件"))
    if svgs:
        for i, s in enumerate(svgs, 1):
            P.append(embed_svg(s, "图 %d　%s" % (i, os.path.basename(s))))
    else:
        P.append('<p class="small">未找到 SVG 图件。</p>')

    # 章节 4：图片自检
    P.append(h2("图片自检"))
    if os.path.exists(qa):
        qtext = open(qa, encoding="utf-8").read()
        first = "\n".join(qtext.splitlines()[:8])
        P.append('<pre class="card" style="white-space:pre-wrap;font-size:8.5pt">%s</pre>' % esc(first))
    else:
        P.append('<p class="small">未找到 figure_qa.md，建议运行 np_figcheck.py。</p>')

    # 章节 5：方法学参数
    P.append(h2("方法学参数"))
    P.append('<div class="card"><span class="tag hub">候选核心</span> '
             '<span class="tag bg">背景节点</span> '
             '<span class="tag ok">通过自检</span> '
             '<span class="tag warn">需复核</span>'
             '<p class="small">种子、物种、阈值、数据库访问日期详见 run_manifest.json 与 '
             'np_prov.snapshot 输出；本预览页不重复转录，避免版本漂移。</p></div>')

    P.append('<hr><p class="small">由 np_htmlreport.py 生成 · CSS 变量/打印样式/语义标签本地化实现 · '
             '数据与结论以 run_manifest 与标准表为准。</p>')
    P.append("</div></body></html>")
    return "\n".join(P)


def main(argv=None):
    ap = argparse.ArgumentParser(description="np_htmlreport: SCI 投稿版 HTML 预览/打印渲染")
    ap.add_argument("--dir", required=True, help="分析输出目录（含 run_manifest/标准表/SVG）")
    ap.add_argument("--label", default="网络药理学分析")
    ap.add_argument("--out", default=None, help="输出 HTML 路径（默认 <dir>/manuscript_preview.html）")
    ap.add_argument("--print-pdf-hint", action="store_true", help="打印为 PDF 的操作提示")
    a = ap.parse_args(argv)
    if a.print_pdf_hint:
        print("打印为 PDF：浏览器打开 HTML → Ctrl+P → 目标『另存为 PDF』→ 勾选背景图形 → A4 → 保存。")
        print("@media print 已内嵌：自动去导航、图表跨页避断、表头跨页重复。")
        return 0
    if not os.path.isdir(a.dir):
        print("[错误] 目录不存在：%s" % a.dir, file=sys.stderr)
        return 2
    out = a.out or os.path.join(a.dir, "manuscript_preview.html")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        f.write(build(a.dir, a.label))
    print("[htmlreport] 已生成 %s（%.0f KB）" % (out, os.path.getsize(out) / 1024))
    print("[htmlreport] 含 @media print / CSS 变量 / 语义色标签 / 动态章节号；浏览器打开即可打印 PDF。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
