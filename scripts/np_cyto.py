#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_cyto.py —— v1.5.0 Cytoscape 工程文件与超网络图输出

背景（23 篇复现痛点）：PPI 网络图/富集气泡图已有 SVG（np_ppi network.svg /
np_enrich enrich_bubble.svg），但论文级"成分-靶点-通路"超网络仍需手动 Cytoscape
排版。本模块输出：
  - xgmml：带节点属性（类型/度值/来源库）+ 圆形布局坐标的 Cytoscape 可导入文件
    （File→Import→Network→XGMML；.cys 是 Cytoscape 会话二进制格式，无法脱离
    Cytoscape 生成——诚实边界，XGMML 是官方推荐的可编程替代）
  - sif：简单交互格式（Cytoscape 全版本兼容）
  - supernetwork.svg：成分-靶点-通路三层超网络矢量图（零依赖，直接进论文）

命令：
  export --edges intersect_edges.tsv --nodes hub_nodes.csv(可选) --paths enrich_significant.csv(可选) --outdir DIR
    --edges 需 node1/node2(/score) 列；--paths 需 term 基因串列或 genes 列
"""
import argparse
import csv
import math
import os

import numpy as np
import pandas as pd


def cmd_export(args):
    os.makedirs(args.outdir, exist_ok=True)
    e = pd.read_csv(args.edges, encoding="utf-8-sig")
    ecols = {c.lower(): c for c in e.columns}
    n1, n2 = ecols.get("node1"), ecols.get("node2")
    sc = ecols.get("score")
    edges = [(r[n1], r[n2], float(r[sc]) if sc and sc in r and r[sc] == r[sc] else None) for _, r in e.iterrows()]
    deg = {}
    for a, b, _ in edges:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    # 节点类型： hubs > target; 通路节点
    path_nodes = {}
    if args.paths and os.path.exists(args.paths):
        p = pd.read_csv(args.paths, encoding="utf-8-sig")
        gcol = next((c for c in p.columns if "gene" in c.lower()), None)
        tcol = next((c for c in p.columns if "term" in c.lower()), None)
        if gcol and tcol:
            for _, r in p.iterrows():
                genes = [g for g in str(r[gcol]).replace(";", " ").replace(",", " ").split() if g]
                if 0 < len(genes) <= 400:
                    path_nodes[str(r[tcol])] = set(genes)
    # SIF
    with open(os.path.join(args.outdir, "network.sif"), "w", encoding="utf-8") as f:
        for a, b, _ in edges:
            f.write(f"{a}\tpp\t{b}\n")
        for t, gs in path_nodes.items():
            for g in gs:
                f.write(f"{t}\tinclude\t{g}\n")
    # XGMML（圆形布局）
    nodes = sorted(deg)
    n = len(nodes) or 1
    pos = {}
    for i, nd in enumerate(nodes):
        ang = 2 * math.pi * i / n
        pos[nd] = (300 + 250 * math.cos(ang), 300 + 250 * math.sin(ang))
    hub_thr = np.percentile(list(deg.values()), 95) if deg else 0
    x = ['<?xml version="1.0" encoding="UTF-8"?>',
         '<graph label="NP network" xmlns="http://www.cs.rpi.edu/~gregod/boost/tools/doc/1.0" xmlns:cy="http://www.cytoscape.org" directed="0">']
    x.append(' <att name="documentVersion" value="1.1"/>')
    for nd in nodes:
        typ = "hub" if deg[nd] >= hub_thr else "target"
        x.append(f' <node id="{nd}" label="{nd}">'
                 f'<att name="type" value="{typ}"/><att name="degree" value="{deg[nd]}"/>'
                 f'<graphics type="ELLIPSE" x="{pos[nd][0]:.0f}" y="{pos[nd][1]:.0f}" w="40" h="40" fill="#{ "c0392b" if typ=="hub" else "2471a3"}"/></node>')
    for a, b, s in edges:
        x.append(f' <edge label="{a}-{b}" source="{a}" target="{b}"'
                 + (f'><att name="score" value="{s:.3f}"/></edge>' if s is not None else '/>'))
    for t in path_nodes:
        x.append(f' <node id="{t}" label="{t[:40]}">'
                 f'<att name="type" value="pathway"/>'
                 f'<graphics type="ELLIPSE" x="300" y="300" w="60" h="60" fill="#1e8449"/></node>')
    x.append('</graph>')
    with open(os.path.join(args.outdir, "network.xgmml"), "w", encoding="utf-8") as f:
        f.write("\n".join(x))
    # 超网络 SVG（成分/靶点 hub + 通路三环）
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 680 680" font-family="sans-serif">',
           '<rect width="680" height="680" fill="#ffffff"/>',
           '<text x="340" y="24" text-anchor="middle" font-size="14" font-weight="bold">成分-靶点-通路超网络（Top 度值靶点）</text>']
    top = [k for k, v in sorted(deg.items(), key=lambda kv: -kv[1])[:18]]
    top_set = set(top)
    # 靶点环
    for i, g in enumerate(top):
        ang = 2 * math.pi * i / len(top)
        x_, y_ = 340 + 220 * math.cos(ang), 350 + 220 * math.sin(ang)
        svg.append(f'<circle cx="{x_:.0f}" cy="{y_:.0f}" r="7" fill="#2471a3"/>')
        svg.append(f'<text x="{x_:.0f}" y="{y_ - 12:.0f}" font-size="9" text-anchor="middle">{g}</text>')
    # 通路外环（取 6 条最显著）
    if path_nodes:
        pw = list(path_nodes)[:6]
        for i, t in enumerate(pw):
            ang = 2 * math.pi * i / len(pw)
            px_, py_ = 340 + 300 * math.cos(ang), 350 + 300 * math.sin(ang)
            shared = top_set & path_nodes[t]
            svg.append(f'<circle cx="{px_:.0f}" cy="{py_:.0f}" r="12" fill="#1e8449"/>')
            svg.append(f'<text x="{px_:.0f}" y="{py_ - 16:.0f}" font-size="8" text-anchor="middle">{t[:24]}</text>')
            for g in shared:
                j = top.index(g)
                ga = 2 * math.pi * j / len(top)
                gx, gy = 340 + 220 * math.cos(ga), 350 + 220 * math.sin(ga)
                svg.append(f'<line x1="{px_:.0f}" y1="{py_:.0f}" x2="{gx:.0f}" y2="{gy:.0f}" stroke="#1e8449" stroke-opacity="0.4"/>')
    svg.append('</svg>')
    open(os.path.join(args.outdir, "supernetwork.svg"), "w", encoding="utf-8").write("".join(svg))
    print(f"节点 {len(nodes)}（hub 阈值=度值P95={hub_thr:.0f}）| 边 {len(edges)} | 通路 {len(path_nodes)}")
    print(f"输出 → {args.outdir}/: network.sif, network.xgmml, supernetwork.svg")


def main():
    ap = argparse.ArgumentParser(description="np_cyto: Cytoscape XGMML/SIF/超网络 SVG")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("export")
    s.add_argument("--edges", required=True)
    s.add_argument("--paths", default=None)
    s.add_argument("--outdir", default=".")
    a = ap.parse_args()
    {"export": cmd_export}[a.cmd](a)


if __name__ == "__main__":
    main()
