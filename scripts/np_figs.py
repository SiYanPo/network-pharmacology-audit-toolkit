#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_figs.py —— v2.0.0 投稿级图片自动生成（矢量 + 600 dpi 位图 + 自检元数据）

为什么单独做：500 篇审计显示 heatmap 46%、bubble 33%、PPI 网络常规化，但作者仍手动
在 Cytoscape/Origin 里排版。本模块直接产出符合期刊要求的四件套：
  fig.pdf（矢量，投稿主用） / fig.svg（可编辑） / fig.png@600dpi / fig.tiff@600dpi(LZW)
并写 fig.json 元数据（尺寸/字号/颜色/文本重叠检测/数据源哈希），供 np_figcheck 自检。

期刊尺寸预设（mm）：single=85、1.5col=114、double=174、full=180
字体：Arial（无则 DejaVu Sans）；正文字号 7pt、标题 8pt、面板字母 8pt bold
配色：Okabe-Ito 色盲友好

命令：
  bubble  --enrich enrich_significant.csv [--top 15] [--width double]
  bar     --enrich enrich_significant.csv [--top 15]
  heatmap --long compound_targets.csv（compound,gene[,value]）[--top-compounds 20 --top-genes 25]
  ppi     --edges ppi_edges_clean.tsv [--top 60]
  prisma  --counts prisma_counts.json
  所有子命令支持 --outdir DIR --name fig1 --formats pdf,svg,png,tiff --dpi 600
"""
import argparse
import hashlib
import json
import os
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

WIDTH_MM = {"single": 85, "1.5col": 114, "double": 174, "full": 180}
MM2IN = 1 / 25.4

# Okabe-Ito 色盲友好
PAL = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]


def _font():
    names = {f.name for f in font_manager.fontManager.ttflist}
    for cand in ("Arial", "Helvetica", "DejaVu Sans"):
        if cand in names:
            return cand
    return "DejaVu Sans"


def _setup(width="double", height_mm=None):
    f = _font()
    plt.rcParams.update({
        "font.family": f, "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
        "axes.linewidth": 0.6, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    })
    w = WIDTH_MM.get(width, 174) * MM2IN
    h = (height_mm or WIDTH_MM.get(width, 174) * 0.62) * MM2IN
    return plt.figure(figsize=(w, h)), f


def _meta(fig, ax_list, source, palette=PAL, fonts=None):
    """绘制后检测：文本重叠、越界、字号"""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # 以 tight bbox（实际保存区域）为边界，避免 bbox_inches="tight" 造成误报
    tb = fig.get_tightbbox(renderer)
    dpi = fig.dpi
    x0px, y0px = tb.x0 * dpi, tb.y0 * dpi
    x1px, y1px = tb.x1 * dpi, tb.y1 * dpi
    texts, overlaps, offcanvas, tick_overhang, sizes = [], 0, 0, 0, []
    for ax in ax_list:
        items = list(ax.texts) + [ax.title, ax.xaxis.label, ax.yaxis.label]
        if getattr(ax, "axison", True):
            items += ax.get_xticklabels() + ax.get_yticklabels()
        boxes = []
        for t in items:
            if not t.get_text().strip():
                continue
            try:
                bb = t.get_window_extent(renderer=renderer)
            except Exception:
                continue
            sizes.append(t.get_fontsize())
            boxes.append(bb)
            if bb.x0 < x0px - 2 or bb.y0 < y0px - 2 or bb.x1 > x1px + 2 or bb.y1 > y1px + 2:
                is_tick = t in (list(ax.get_xticklabels()) + list(ax.get_yticklabels())) if getattr(ax, "axison", True) else False
                if is_tick:
                    tick_overhang += 1  # 刻度标签轻微外溢：保存时 tight bbox 会包含，记 WARN 级
                else:
                    offcanvas += 1
                if os.environ.get("NP_FIGS_DEBUG"):
                    print("OFFCANVAS:", repr(t.get_text())[:30], "tick" if is_tick else "label")
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                a, b = boxes[i], boxes[j]
                if a.overlaps(b):
                    iw = min(a.x1, b.x1) - max(a.x0, b.x0)
                    ih = min(a.y1, b.y1) - max(a.y0, b.y0)
                    if iw > 2 and ih > 2:
                        overlaps += 1
                        if os.environ.get("NP_FIGS_DEBUG"):
                            print("OVERLAP:", repr(t.get_text())[:40])
    return {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "source_file": os.path.basename(source) if source else None,
        "source_sha256_16": _sha(source) if source else None,
        "fig_size_inch": [round(x, 3) for x in fig.get_size_inches()],
        "fig_size_mm": [round(x * 25.4, 1) for x in fig.get_size_inches()],
        "font_family": fonts or _font(),
        "min_font_pt": round(min(sizes), 1) if sizes else None,
        "max_font_pt": round(max(sizes), 1) if sizes else None,
        "n_text": len(sizes),
        "text_overlaps": overlaps,
        "text_offcanvas": offcanvas,
        "tick_label_overhang": tick_overhang,
        "palette": palette[:6],
    }


def _sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(65536), b""):
            h.update(c)
    return h.hexdigest()[:16]


def _save(fig, outdir, name, formats, dpi, meta):
    os.makedirs(outdir, exist_ok=True)
    files = []
    for fmt in formats:
        p = os.path.join(outdir, f"{name}.{fmt}")
        kw = {"format": fmt}
        if fmt in ("png", "tif", "tiff"):
            from io import BytesIO
            from PIL import Image as _Img
            buf = BytesIO()
            fig.savefig(buf, format="png", dpi=dpi, facecolor="white")
            buf.seek(0)
            im = _Img.open(buf).convert("RGB")
            if fmt.startswith("ti"):
                im.save(p, format="TIFF", dpi=(dpi, dpi), compression="tiff_lzw")
            else:
                im.save(p, format="PNG", dpi=(dpi, dpi))
            files.append(p)
            continue
        fig.savefig(p, **kw)
        files.append(p)
    meta["files"] = [os.path.basename(f) for f in files]
    meta["dpi"] = dpi if any(f.endswith(("png", "tif", "tiff")) for f in files) else None
    json.dump(meta, open(os.path.join(outdir, f"{name}.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"✓ {name}: " + ", ".join(meta["files"]) + f" | 尺寸 {meta['fig_size_mm']} mm"
          + f" | 文本重叠 {meta['text_overlaps']} 处")


# ---------------- 图型 ----------------

def cmd_bubble(a):
    df = pd.read_csv(a.enrich, encoding="utf-8-sig")
    df = df.sort_values("pvalue").head(a.top).copy()
    df["logp"] = -np.log10(df.pvalue.clip(lower=1e-300))
    df["ratio"] = df.n_hit / df.term_size
    fig, _ = _setup(a.width, height_mm=max(64, 6.0 * a.top))
    ax = fig.add_axes([0.42, 0.12, 0.46, 0.80])
    cax = fig.add_axes([0.91, 0.12, 0.018, 0.80])
    y = np.arange(len(df))[::-1]
    sc = ax.scatter(df.ratio, y, s=df.n_hit * 6, c=df.logp, cmap="viridis",
                    edgecolor="white", linewidth=0.5, vmin=df.logp.min(), vmax=df.logp.max())
    ax.set_yticks(y)
    import textwrap as _tw
    ax.set_yticklabels([_tw.fill(t, 44) for t in df.term], fontsize=6.0)
    ax.set_xlabel("Gene ratio (n hit / term size)")
    # 气泡大小图例（期刊常规要求）
    for q in (df.n_hit.min(), df.n_hit.median(), df.n_hit.max()):
        ax.scatter([], [], s=q * 6, c="#bbb", edgecolor="white", label=f"{int(q)} genes")
    leg = ax.legend(title="Count", loc="lower left", bbox_to_anchor=(0.02, 0.02), frameon=False, fontsize=6,
                    title_fontsize=6, labelspacing=0.9, handletextpad=0.6, borderpad=0.2)
    for t in leg.get_texts():
        t.set_fontsize(6)
    ax.set_title(a.title or f"Enrichment bubble (Top {len(df)})", loc="left", fontweight="bold")
    ax.grid(axis="x", ls=":", lw=0.4, alpha=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("−log10(P)", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    cb.outline.set_linewidth(0.4)
    meta = _meta(fig, [ax], a.enrich, PAL)
    meta["type"] = "bubble"; meta["n_terms"] = len(df)
    _save(fig, a.outdir, a.name, a.formats, a.dpi, meta)
    plt.close(fig)


def cmd_bar(a):
    df = pd.read_csv(a.enrich, encoding="utf-8-sig").sort_values("pvalue").head(a.top).copy()
    df["logp"] = -np.log10(df.pvalue.clip(lower=1e-300))
    fig, _ = _setup(a.width, height_mm=max(52, 5.4 * a.top))
    ax = fig.add_axes([0.44, 0.13, 0.49, 0.78])
    y = np.arange(len(df))[::-1]
    ax.barh(y, df.logp, color=PAL[0], height=0.62, edgecolor="white", lw=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels([(t[:46] + "…") if len(t) > 47 else t for t in df.term], fontsize=6.0)
    ax.set_xlabel("−log10(P)")
    ax.set_xlim(0, float(df.logp.max()) * 1.08)
    ax.set_title(a.title or f"Top {len(df)} enriched terms", loc="left", fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", ls=":", lw=0.4, alpha=0.5)
    meta = _meta(fig, [ax], a.enrich, PAL)
    meta["type"] = "bar"; meta["n_terms"] = len(df)
    _save(fig, a.outdir, a.name, a.formats, a.dpi, meta)
    plt.close(fig)


def cmd_heatmap(a):
    df = pd.read_csv(a.long, encoding="utf-8-sig")
    c0, c1 = df.columns[0], df.columns[1]
    val = df.columns[2] if len(df.columns) > 2 else None
    if val is not None:
        try:
            pd.to_numeric(df[val])
        except (ValueError, TypeError):
            val = None  # 第三列是注释文本（如 target_name），按出现次数计数
    if val is None:
        df = df.copy(); df["_n"] = 1
        mat = df.pivot_table(index=c0, columns=c1, values="_n", aggfunc="sum")
    else:
        mat = df.pivot_table(index=c0, columns=c1, values=val, aggfunc="max")
    mat = mat.fillna(0).astype(float)
    # 选信息量最大的行列
    mat = mat.loc[mat.sum(axis=1).sort_values(ascending=False).index[:a.top_compounds]]
    mat = mat[mat.sum(axis=0).sort_values(ascending=False).index[:a.top_genes]]
    # 简单层次排序（按行均值）
    mat = mat.loc[mat.mean(axis=1).sort_values(ascending=False).index]
    fig, _ = _setup(a.width, height_mm=max(55, 3.2 * len(mat) + 18))
    ax = fig.add_axes([0.34, 0.22, 0.60, 0.70])
    im = ax.imshow(mat.values, aspect="auto", cmap="YlGnBu", interpolation="nearest")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels([str(i)[:26] for i in mat.index], fontsize=6)
    ax.set_xlabel(a.xlabel or "Targets")
    ax.set_ylabel(a.ylabel or "Compounds")
    ax.set_title(a.title or "Compound–target matrix", loc="left", fontweight="bold")
    for s in ax.spines.values():
        s.set_linewidth(0.4)
    cax = fig.add_axes([0.955, 0.22, 0.012, 0.70])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(a.cbar or "Evidence (count / score)", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    cb.outline.set_linewidth(0.4)
    meta = _meta(fig, [ax], a.long, ["#F7FBFF", "#08306B"])
    meta["type"] = "heatmap"; meta["shape"] = list(mat.shape)
    _save(fig, a.outdir, a.name, a.formats, a.dpi, meta)
    plt.close(fig)


def cmd_ppi(a):
    import networkx as nx
    e = pd.read_csv(a.edges, sep=None, engine="python", encoding="utf-8-sig")
    c = {x.lower(): x for x in e.columns}
    n1, n2 = c.get("node1"), c.get("node2")
    w = c.get("weight") or c.get("score")
    G = nx.Graph()
    for _, r in e.iterrows():
        if w and pd.notna(r.get(w)):
            G.add_edge(r[n1], r[n2], weight=float(r[w]))
        else:
            G.add_edge(r[n1], r[n2], weight=1.0)
    # 取度最高的子图，保留连通性最好的分量
    top = {n for n, _ in sorted(dict(G.degree).items(), key=lambda kv: -kv[1])[: a.top]}
    H = G.subgraph(top).copy()
    # 按连通分量分块布局（避免丢弃孤立簇，也避免 spring 把簇拉得过远）
    comps = sorted(nx.connected_components(H), key=len, reverse=True)
    ncomp = len(comps)
    cols = int(np.ceil(np.sqrt(ncomp))) if ncomp else 1
    rows = int(np.ceil(ncomp / cols)) if ncomp else 1
    pos = {}
    for ci, comp in enumerate(comps):
        sub = H.subgraph(comp)
        if len(comp) == 1:
            local = {next(iter(comp)): np.array([0.0, 0.0])}
        else:
            local = nx.spring_layout(sub, k=a.spring_k, iterations=600, seed=42,
                                     weight="weight", scale=1.0)
        xs = np.array([local[n][0] for n in comp]); ys = np.array([local[n][1] for n in comp])
        def _n(v, lo, hi):
            rng = (v.max() - v.min()) or 1.0
            return lo + (hi - lo) * (v - v.min()) / rng
        r_, c_ = divmod(ci, cols)
        cx = 0.12 + 0.76 * ((c_ + 0.5) / cols)
        cy = 0.90 - 0.80 * ((r_ + 0.5) / rows)
        half_w = 0.38 / cols * (0.95 if ncomp > 1 else 1.0)
        half_h = 0.40 / rows * (0.95 if ncomp > 1 else 1.0)
        if len(comp) >= 3:
            px = _n(xs, cx - half_w, cx + half_w); py = _n(ys, cy - half_h, cy + half_h)
        elif len(comp) == 2:
            px = np.array([cx - 0.035, cx + 0.035]); py = np.array([cy, cy])
        else:
            px = np.array([cx]); py = np.array([cy])
        for i, n in enumerate(sorted(comp)):
            pos[n] = np.array([px[i], py[i]])
    fig, _ = _setup(a.width, height_mm=WIDTH_MM[a.width] * 0.80)
    ax = fig.add_axes([0.03, 0.05, 0.94, 0.87]); ax.set_axis_off()
    for u, v, d in H.edges(data=True):
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color="#9aa5b1", lw=0.35 + 0.45 * min(d.get("weight", 1), 1.0), alpha=0.72, zorder=1)
    deg = dict(H.degree)
    ax.scatter([pos[n][0] for n in H], [pos[n][1] for n in H],
               s=[30 + 7 * deg[n] for n in H], c=PAL[0], alpha=0.92,
               edgecolor="white", linewidth=0.6, zorder=2)
    for n in sorted(H.nodes, key=lambda x: -deg[x])[: a.label_top]:
        ax.text(pos[n][0], pos[n][1] + 0.045, n, fontsize=6.2, ha="center", va="bottom", zorder=3)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title(f"{a.title or 'PPI network'} | {H.number_of_nodes()} nodes, {H.number_of_edges()} edges",
                 loc="left", fontweight="bold")
    meta = _meta(fig, [ax], a.edges, PAL)
    meta["type"] = "ppi"; meta["n_nodes"] = H.number_of_nodes(); meta["n_edges"] = H.number_of_edges()
    _save(fig, a.outdir, a.name, a.formats, a.dpi, meta)
    plt.close(fig)

def cmd_prisma(a):
    c = json.load(open(a.counts, encoding="utf-8"))
    steps = [
        ("Records identified (PubMed/CNKI)", c.get("identified")),
        ("Duplicates removed", c.get("duplicates")),
        ("Screened by title/abstract", c.get("screened")),
        ("Excluded (not NP-relevant)", c.get("excluded")),
        ("Full texts assessed", c.get("fulltext")),
        ("Included in reproduction", c.get("included")),
    ]
    fig, _ = _setup("double", height_mm=118)
    ax = fig.add_axes([0.06, 0.04, 0.88, 0.86]); ax.set_axis_off()
    ax.set_xlim(1.0, 6.6); ax.set_ylim(0.4, 9.9)
    y = 9.0
    for i, (lab, val) in enumerate(steps):
        if val is None:
            continue
        ax.add_patch(plt.Rectangle((1.15, y - 0.85), 5.2, 1.15, facecolor="#EAF2FB",
                                   edgecolor=PAL[0], lw=0.8))
        ax.text(3.75, y - 0.28, lab, fontsize=7.0, ha="center", va="center", fontweight="bold")
        ax.text(3.75, y - 0.62, f"n = {val}", fontsize=7.4, ha="center", va="center")
        if i < len(steps) - 1:
            ax.annotate("", xy=(3.75, y - 1.55), xytext=(3.75, y - 0.85),
                        arrowprops=dict(arrowstyle="-|>", lw=0.7, color="#444"))
        y -= 1.55
    ax.set_title(a.title or "Study selection flow (PRISMA-style)", loc="left",
                 fontweight="bold", y=0.98)
    meta = _meta(fig, [ax], a.counts, PAL)
    meta["type"] = "prisma"
    _save(fig, a.outdir, a.name, a.formats, a.dpi, meta)
    plt.close(fig)



def cmd_batch(a):
    """一键生成全部可用图型并自检（全流程自动化入口）"""
    made = []
    if a.enrich:
        ns = argparse.Namespace(enrich=a.enrich, top=a.top, outdir=a.outdir, name="fig_bubble",
                                width=a.width, formats=a.formats, dpi=a.dpi, title=None)
        cmd_bubble(ns); made.append("fig_bubble")
        ns2 = argparse.Namespace(enrich=a.enrich, top=min(a.top, 10), outdir=a.outdir, name="fig_bar",
                                 width=a.width, formats=a.formats, dpi=a.dpi, title=None)
        cmd_bar(ns2); made.append("fig_bar")
    if a.long:
        ns = argparse.Namespace(long=a.long, top_compounds=a.top_compounds, top_genes=a.top_genes,
                                outdir=a.outdir, name="fig_heatmap", width=a.width, formats=a.formats,
                                dpi=a.dpi, title=None, xlabel=None, ylabel=None, cbar=None)
        cmd_heatmap(ns); made.append("fig_heatmap")
    if a.edges:
        ns = argparse.Namespace(edges=a.edges, top=a.top_nodes, label_top=a.label_top, spring_k=0.9,
                                outdir=a.outdir, name="fig_ppi", width=a.width, formats=a.formats,
                                dpi=a.dpi, title=None)
        cmd_ppi(ns); made.append("fig_ppi")
    if a.counts:
        ns = argparse.Namespace(counts=a.counts, outdir=a.outdir, name="fig_flow", width=a.width,
                                formats=a.formats, dpi=a.dpi, title=None)
        cmd_prisma(ns); made.append("fig_flow")
    if not made:
        raise SystemExit("[FAIL] 至少提供 --enrich / --long / --edges / --counts 之一")
    print(f"\n生成 {len(made)} 张图: {', '.join(made)}")
    if a.check:
        import subprocess as _sp
        _sp.run([os.sys.executable, os.path.join(os.path.dirname(__file__), "np_figcheck.py"),
                 "--dir", a.outdir, "--out", os.path.join(a.outdir, "figure_qa.md")])


def main():
    ap = argparse.ArgumentParser(description="np_figs: 投稿级图片生成")
    sub = ap.add_subparsers(dest="cmd", required=True)
    def common(s):
        s.add_argument("--outdir", default="figs")
        s.add_argument("--name", default=None)
        s.add_argument("--width", default="double", choices=list(WIDTH_MM))
        s.add_argument("--formats", default="pdf,svg,png,tiff")
        s.add_argument("--dpi", type=int, default=600)
        s.add_argument("--title", default=None)
    s = sub.add_parser("bubble"); s.add_argument("--enrich", required=True)
    s.add_argument("--top", type=int, default=15); common(s)
    s = sub.add_parser("bar"); s.add_argument("--enrich", required=True)
    s.add_argument("--top", type=int, default=15); common(s)
    s = sub.add_parser("heatmap"); s.add_argument("--long", required=True)
    s.add_argument("--top-compounds", type=int, default=18)
    s.add_argument("--top-genes", type=int, default=24)
    s.add_argument("--xlabel", default=None); s.add_argument("--ylabel", default=None)
    s.add_argument("--cbar", default=None); common(s)
    s = sub.add_parser("ppi"); s.add_argument("--edges", required=True)
    s.add_argument("--top", type=int, default=60); s.add_argument("--label-top", type=int, default=18)
    s.add_argument("--spring-k", type=float, default=0.9)
    common(s)
    s = sub.add_parser("prisma"); s.add_argument("--counts", required=True); common(s)
    s = sub.add_parser("batch")
    s.add_argument("--enrich", default=None); s.add_argument("--long", default=None)
    s.add_argument("--edges", default=None); s.add_argument("--counts", default=None)
    s.add_argument("--top", type=int, default=12); s.add_argument("--top-nodes", type=int, default=40)
    s.add_argument("--top-compounds", type=int, default=14); s.add_argument("--top-genes", type=int, default=20)
    s.add_argument("--label-top", type=int, default=18)
    s.add_argument("--check", action="store_true", default=True)
    common(s)
    a = ap.parse_args()
    a.formats = [f.strip() for f in a.formats.split(",") if f.strip()]
    a.name = getattr(a, "name", None) or {"bubble": "fig_bubble", "bar": "fig_bar", "heatmap": "fig_heatmap",
                        "ppi": "fig_ppi", "prisma": "fig_flow", "batch": "fig"}[a.cmd]
    {"bubble": cmd_bubble, "bar": cmd_bar, "heatmap": cmd_heatmap,
     "ppi": cmd_ppi, "prisma": cmd_prisma, "batch": cmd_batch}[a.cmd](a)


if __name__ == "__main__":
    main()
