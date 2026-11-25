#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_interactive.py —— 交互式网络图（Cytoscape.js，单文件 HTML）

为什么放弃静态 SVG
    几百个节点的 PPI 用 SVG 画出来必然是 hairball：节点重叠、标签互相压、想放大看某个 hub 只能截图。
    软件类文章的审稿人一定会说"图太乱、体验差"。Cytoscape.js 给的是可缩放、可拖拽、可搜索、
    点击高亮邻居的交互图——同一个数据，信息密度和可用性完全不同。

设计
    * 输出**单个自包含 HTML**：数据以 JSON 内嵌，Cytoscape.js 默认走 CDN；
      需要完全离线时用 `--inline-js cytoscape.min.js` 把库内联进文件（约 300 KB）。
    * 节点属性（degree / hub rank / 共识 verdict / 是否核心）来自 hubs.csv 与交集基因表，
      点击节点在右侧面板显示，不靠 tooltip 猜。
    * 视觉编码固定：大小 = degree，颜色 = 核心/背景/不稳定，边宽 = 置信度。
    * 提供搜索、布局切换、只看核心、导出 PNG —— 都是审稿人现场会点的按钮。

用法
    python np_interactive.py --edges ppi_edges_clean.tsv --hubs hubs.csv \
        --core intersection_genes.txt --out network_interactive.html \
        --title "复方丹参滴丸 PPI" [--inline-js cytoscape.min.js]
    python np_interactive.py selftest
"""
from __future__ import annotations

import argparse
import csv
import glob
import html
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CDN = "https://unpkg.com/cytoscape@3.30.2/dist/cytoscape.min.js"


def _read_edges(path):
    """PPI 边表 → [(a, b, w)]；自动识别 node1/node2、source/target、preferredName_A/B。"""
    rows = []
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        rdr = csv.DictReader(f, delimiter=delim)
        low = {(h or "").lower().strip(): h for h in (rdr.fieldnames or [])}
        a = next((low[k] for k in ("node1", "source", "preferredname_a", "from", "gene_a") if k in low), None)
        b = next((low[k] for k in ("node2", "target", "preferredname_b", "to", "gene_b") if k in low), None)
        w = next((low[k] for k in ("weight", "combined_score", "score", "confidence") if k in low), None)
        if not (a and b):
            raise SystemExit("[FAIL] 边表需含 node1/node2（或 source/target）列：%s" % list(low))
        for r in rdr:
            x = (r.get(a) or "").strip().upper()
            y = (r.get(b) or "").strip().upper()
            if not x or not y or x == y:
                continue
            try:
                wt = float(r.get(w)) if w else 1.0
            except (TypeError, ValueError):
                wt = 1.0
            rows.append((x, y, wt))
    return rows


def _read_hubs(path):
    """hubs.csv → {gene: {rank, degree, betweenness, verdict}}。"""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        for r in csv.DictReader(f, delimiter=delim):
            g = (r.get("node") or r.get("gene") or "").strip().upper()
            if not g:
                continue
            out[g] = {k: r.get(k, "") for k in
                      ("rank", "degree", "betweenness", "closeness", "pagerank", "verdict") if k in r}
    return out


def _read_list(path):
    if not path or not os.path.exists(path):
        return set()
    return {l.strip().upper() for l in open(path, encoding="utf-8-sig", errors="replace") if l.strip()}


def build(edges_path, hubs_path=None, core_path=None, title="PPI 交互网络",
          inline_js=None, max_edges=4000):
    """返回完整 HTML 字符串（单文件）。"""
    edges = _read_edges(edges_path)
    truncated = False
    if len(edges) > max_edges:
        edges = sorted(edges, key=lambda x: -x[2])[:max_edges]
        truncated = True
    hubs = _read_hubs(hubs_path)
    core = _read_list(core_path)

    nodes = {}
    for a, b, w in edges:
        for g in (a, b):
            nodes.setdefault(g, 0)
            nodes[g] += 1
    for g in hubs:
        nodes.setdefault(g, int(float(hubs[g].get("degree") or 0) or 0))

    node_data = []
    for g, deg in sorted(nodes.items()):
        h = hubs.get(g, {})
        is_core = g in core
        verdict = h.get("verdict", "")
        node_data.append({
            "id": g, "label": g, "degree": deg,
            "rank": h.get("rank", ""), "verdict": verdict,
            "betweenness": h.get("betweenness", ""),
            "class": "core" if is_core else ("hub" if str(h.get("rank", "999")) in
                                             [str(i) for i in range(1, 11)] else "bg"),
        })
    edge_data = [{"source": a, "target": b, "weight": round(w, 3),
                  "width": round(0.6 + min(3.4, max(0.0, (w - 0.4)) * 3), 2)} for a, b, w in edges]

    n_core = sum(1 for n in node_data if n["class"] == "core")
    lib = ""
    if inline_js and os.path.exists(inline_js):
        lib = "<script>%s</script>" % open(inline_js, encoding="utf-8", errors="replace").read()
    else:
        lib = '<script src="%s"></script>' % CDN

    payload = json.dumps({"nodes": node_data, "edges": edge_data}, ensure_ascii=False)
    tpl = open(os.path.join(HERE, "assets_interactive_tpl.html"), encoding="utf-8").read() \
        if os.path.exists(os.path.join(HERE, "assets_interactive_tpl.html")) else TEMPLATE
    return (tpl.replace("__TITLE__", html.escape(title))
               .replace("__LIB__", lib)
               .replace("__DATA__", payload)
               .replace("__N_NODES__", str(len(node_data)))
               .replace("__N_EDGES__", str(len(edge_data)))
               .replace("__N_CORE__", str(n_core))
               .replace("__TRUNC__", ("边数超过上限，已按置信度保留前 %d 条" % max_edges) if truncated else ""))


TEMPLATE = r"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8">
<title>__TITLE__</title>
<style>
:root{--ink:#1f2937;--muted:#6b7280;--line:#e5e7eb;--accent:#007B7A;--core:#c0504d;--hub:#007B7A;--bg:#9aa3ad}
*{box-sizing:border-box}
html,body{margin:0;height:100%;font-family:"Helvetica Neue",Arial,"Microsoft YaHei",sans-serif;color:var(--ink)}
#app{display:flex;height:100%;background:#fff}
#left{flex:1;position:relative;min-width:0}
#cy{position:absolute;inset:0}
#bar{position:absolute;top:10px;left:10px;right:10px;display:flex;gap:8px;flex-wrap:wrap;align-items:center;
     background:rgba(255,255,255,.94);border:1px solid var(--line);border-radius:8px;padding:8px 10px;z-index:5;
     font-size:12px;box-shadow:0 2px 8px #0001}
#bar input{padding:5px 8px;border:1px solid var(--line);border-radius:6px;font-size:12px;width:170px}
#bar button{padding:5px 10px;border:1px solid var(--line);background:#f8fafc;border-radius:6px;cursor:pointer;font-size:12px}
#bar button:hover{background:#eef2f7}
.hint{color:var(--muted)}
#panel{width:280px;border-left:1px solid var(--line);padding:14px;overflow:auto;font-size:12.5px;background:#fcfdfe}
#panel h3{margin:0 0 6px;font-size:14px;color:var(--accent)}
#panel .kv{display:flex;justify-content:space-between;padding:3px 0;border-bottom:1px dashed #eef1f5}
#panel .kv span:first-child{color:var(--muted)}
.legend{display:flex;gap:10px;align-items:center;font-size:11.5px;color:var(--muted)}
.dot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:4px}
#err{display:none;position:absolute;inset:0;padding:30px;background:#fff;z-index:9;color:#b91c1c;font-size:13px;line-height:1.7}
.badge{display:inline-block;padding:1px 7px;border-radius:8px;font-size:11px;background:#eef7f0;color:#4f8a5b}
</style>
__LIB__
</head><body>
<div id="app">
  <div id="left">
    <div id="bar">
      <input id="q" placeholder="搜索基因，如 IL6">
      <button id="btn-fit">重置视图</button>
      <button id="btn-layout">切换布局</button>
      <button id="btn-core">只看核心</button>
      <button id="btn-png">导出 PNG</button>
      <button id="btn-sif" title="导出当前可见子网的 SIF + 节点属性表，可直接导入 Cytoscape 桌面版">导出 SIF</button>
      <span class="legend"><i class="dot" style="background:var(--core)"></i>核心靶点 __N_CORE__
        <i class="dot" style="background:var(--hub)"></i>Top10 hub <i class="dot" style="background:var(--bg)"></i>背景
        ｜节点 __N_NODES__ · 边 __N_EDGES__</span>
      <span class="hint" id="trunc">__TRUNC__</span>
      <span class="hint" id="sif-tip"></span>
    </div>
    <div id="cy"></div>
    <div id="err">Cytoscape.js 未能加载（离线环境）。<br>请在联网环境打开本文件，或用
      <code>--inline-js cytoscape.min.js</code> 重新生成以完全离线。</div>
  </div>
  <div id="panel">
    <h3>节点属性</h3>
    <div id="info" class="hint">点击任意节点查看 degree / hub 排名 / 共识判定；点击空白处取消高亮。</div>
    <hr style="border:none;border-top:1px solid var(--line);margin:12px 0">
    <div class="hint">操作：滚轮缩放 · 拖拽平移 · 拖动节点 · 点击节点高亮其一阶邻居 · 双击空白适配全图</div>
    <div class="hint" style="margin-top:8px">导出 SIF：只导当前可见节点/边（先用搜索或「只看核心」筛选），
      连同属性表一起拖入 Cytoscape 桌面版即可继续排版。</div>
  </div>
</div>
<script>
(function(){
  var DATA = __DATA__;
  if (typeof cytoscape === "undefined") { document.getElementById("err").style.display = "block"; return; }
  var COLORS = {core:"#c0504d", hub:"#007B7A", bg:"#9aa3ad"};
  var cy = cytoscape({
    container: document.getElementById("cy"),
    elements: { nodes: DATA.nodes.map(function(n){ return {data:n}; }),
                edges: DATA.edges.map(function(e){ return {data:e}; }) },
    style: [
      {selector:"node", style:{
        "background-color": function(ele){ return COLORS[ele.data("class")] || COLORS.bg; },
        "width": function(ele){ return Math.min(46, 12 + Math.sqrt(ele.data("degree")||1) * 6); },
        "height": function(ele){ return Math.min(46, 12 + Math.sqrt(ele.data("degree")||1) * 6); },
        "label":"data(label)", "font-size":8, "color":"#374151",
        "text-valign":"bottom", "text-margin-y":3, "text-outline-width":2, "text-outline-color":"#fff",
        "border-width":1, "border-color":"#ffffff"}},
      {selector:"node:selected", style:{"border-width":3, "border-color":"#111827"}},
      {selector:"edge", style:{"width":"data(width)", "line-color":"#cfd6de",
        "curve-style":"haystack", "opacity":.8}},
      {selector:".faded", style:{"opacity":.12}},
      {selector:".hl", style:{"line-color":"#007B7A", "opacity":1, "width":2.6}},
      {selector:"node.hl", style:{"border-width":3, "border-color":"#007B7A"}}
    ],
    wheelSensitivity: 0.2,
    layout: {name:"cose", animate:false, nodeRepulsion:9000, idealEdgeLength:70, numIter:800}
  });

  var layouts = ["cose","concentric","circle","grid","breadthfirst"], li = 0;
  document.getElementById("btn-layout").onclick = function(){
    li = (li + 1) % layouts.length;
    cy.layout({name: layouts[li], animate:true, animationDuration:400,
               nodeRepulsion:9000, idealEdgeLength:70, numIter:600}).run();
  };
  document.getElementById("btn-fit").onclick = function(){ cy.fit(null, 40); cy.elements().removeClass("faded hl"); };
  function download(name, text, mime){
    var blob = new Blob([text], {type: mime || "text/plain;charset=utf-8"});
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = name;
    document.body.appendChild(a); a.click();
    setTimeout(function(){ URL.revokeObjectURL(a.href); a.remove(); }, 500);
  }
  document.getElementById("btn-sif").onclick = function(){
    // 只导出"当前可见"的元素：用户筛选过的子网才是他想要带走的
    var ns = cy.nodes(":visible"), es = cy.edges(":visible");
    var sif = ["# source\tinteraction\ttarget"];
    es.forEach(function(e){ sif.push(e.source().id() + "\tpp\t" + e.target().id()); });
    var attrs = ["id\tdegree\thub_rank\tverdict\tclass\tbetweenness"];
    ns.forEach(function(n){
      var d = n.data();
      attrs.push([d.id, d.degree, d.rank || "", d.verdict || "", d.class || "", d.betweenness || ""].join("\t"));
    });
    download("network_visible.sif", sif.join("\n"));
    setTimeout(function(){ download("network_visible_attributes.csv", attrs.join("\n"), "text/csv;charset=utf-8"); }, 250);
    var tip = document.getElementById("sif-tip");
    if (tip) { tip.textContent = "已导出 " + es.length + " 条边 / " + ns.length + " 个节点（SIF + 属性表）"; }
  };
  document.getElementById("btn-png").onclick = function(){
    var png = cy.png({full:true, scale:2, bg:"#ffffff"});
    var a = document.createElement("a"); a.href = png; a.download = "network.png"; a.click();
  };
  var coreOnly = false;
  document.getElementById("btn-core").onclick = function(){
    coreOnly = !coreOnly;
    cy.nodes().forEach(function(n){
      if (coreOnly && n.data("class") !== "core") n.style("display","none"); else n.style("display","element");
    });
    cy.fit(null, 40);
  };
  var q = document.getElementById("q");
  q.oninput = function(){
    var v = q.value.trim().toUpperCase();
    cy.elements().removeClass("faded hl");
    if (!v) return;
    var hits = cy.nodes().filter(function(n){ return n.id().indexOf(v) >= 0; });
    if (hits.length) {
      cy.elements().addClass("faded");
      hits.removeClass("faded").addClass("hl");
      hits.connectedEdges().removeClass("faded").addClass("hl");
      hits.connectedNodes().removeClass("faded");
      cy.fit(hits, 120);
    }
  };
  cy.on("tap", "node", function(evt){
    var n = evt.target, d = n.data();
    cy.elements().removeClass("faded hl");
    cy.elements().addClass("faded");
    n.removeClass("faded").addClass("hl");
    n.connectedEdges().removeClass("faded").addClass("hl");
    n.neighborhood("node").removeClass("faded");
    var rows = [["基因", d.id], ["连接度 degree", d.degree], ["hub 排名", d.rank || "未入榜"],
                ["共识判定", d.verdict || "—"], ["中介中心性", d.betweenness || "—"],
                ["类别", d.class === "core" ? "核心靶点" : (d.class === "hub" ? "Top10 hub" : "背景")]];
    document.getElementById("info").innerHTML = rows.map(function(r){
      return '<div class="kv"><span>'+r[0]+'</span><b>'+r[1]+'</b></div>';
    }).join("") + '<div style="margin-top:10px"><span class="badge">一阶邻居 '+n.neighborhood("node").length+' 个</span></div>';
  });
  cy.on("tap", function(evt){ if (evt.target === cy) cy.elements().removeClass("faded hl"); });
})();
</script>
</body></html>
"""


def selftest():
    import tempfile
    d = tempfile.mkdtemp(prefix="np_interactive_")
    e = os.path.join(d, "edges.csv")
    open(e, "w", encoding="utf-8").write("node1,node2,weight\nIL6,AKT1,0.9\nIL6,TNF,0.8\nAKT1,MTOR,0.7\n")
    h = os.path.join(d, "hubs.csv")
    open(h, "w", encoding="utf-8").write("rank,node,degree,verdict\n1,IL6,2,稳健\n2,AKT1,2,稳健\n3,MTOR,1,不稳定\n")
    c = os.path.join(d, "core.txt")
    open(c, "w", encoding="utf-8").write("IL6\nTNF\n")
    out = build(e, h, c, title="test")
    checks = [("has_lib", "cytoscape" in out),
              ("has_nodes", '"id": "IL6"' in out or '"id":"IL6"' in out),
              ("core_count", "__N_CORE__" not in out and '"class": "core"' in out),
              ("has_panel", "节点属性" in out),
              ("no_placeholder", "__DATA__" not in out and "__TITLE__" not in out)]
    ok = all(x[1] for x in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(x[1] for x in checks), len(checks)))
    for n, g in checks:
        print("  [%s] %s" % ("✓" if g else "✗", n))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = argparse.ArgumentParser(description="np_interactive: Cytoscape.js 交互式 PPI")
    ap.add_argument("--edges", required=True)
    ap.add_argument("--hubs", default=None)
    ap.add_argument("--core", default=None, help="核心/交集基因列表")
    ap.add_argument("--title", default="PPI 交互网络")
    ap.add_argument("--inline-js", default=None, help="本地 cytoscape.min.js（离线内联）")
    ap.add_argument("--max-edges", type=int, default=4000)
    ap.add_argument("-o", "--out", default="network_interactive.html")
    a = ap.parse_args(argv)
    html_txt = build(a.edges, a.hubs, a.core, a.title, a.inline_js, a.max_edges)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w", encoding="utf-8").write(html_txt)
    print("[interactive] → %s（%.0f KB）" % (a.out, os.path.getsize(a.out) / 1024))
    if not a.inline_js:
        print("[interactive] 提示：默认从 CDN 加载 cytoscape.js；离线场景请加 --inline-js cytoscape.min.js")
    return 0


if __name__ == "__main__":
    sys.exit(main())
