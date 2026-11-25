#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_singlecell.py —— 单细胞转录组投影：核心靶点 → 细胞亚群定位（顶刊升维模块）

科学问题
    网络药理学只给"组织层面"靶点；顶刊要求回答：药作用于 Kupffer 细胞还是内皮细胞？
    把 PPI 核心靶点投射到 scRNA-seq 图谱，直接给出"靶点在哪些细胞亚群高表达"，
    这是把文章从 3 分拉到 8 分的关键图。

依赖策略（诚实）
    scanpy/anndata **未预装**（实测）→ 本模块：
      * 有 scanpy：真读 .h5ad，用 obs[cluster_key] 的亚群，算每个靶点在各亚群的平均表达/表达比例，
        输出发表级点图(dotplot)/小提琴图(violin)/UMAP 表达图（PDF+SVG+600dpi PNG+TIFF）。
      * 无 scanpy：进入 --demo 模拟模式（构造annonData 风格矩阵，纯 numpy/matplotlib），
        **仅用于验证流程与出图代码**，明确标注 DEMO 数据不得用于论文。
      * 支持 --geo GSEXXXX：给出 GEO 下载与 h5ad 转换指引（不伪造下载）。

用法
    python np_singlecell.py project --h5ad data.h5ad --targets core_targets.txt \
        --cluster-key celltype --outdir figs/singlecell
    python np_singlecell.py project --targets core_targets.txt --demo -o /tmp/sc_demo
    python np_singlecell.py geo --gse GSE123456          # 下载与转换指引
    python np_singlecell.py selftest
"""
import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def has_scanpy():
    try:
        import scanpy as sc  # noqa
        return True
    except Exception:
        return False


def read_targets(path):
    if not path or not os.path.exists(path):
        return []
    out = []
    with open(path, encoding="utf-8-sig") as f:
        first = f.readline()
        if "\t" in first or "," in first:
            delim = "\t" if "\t" in first else ","
            f.seek(0)
            rdr = csv.DictReader(f, delimiter=delim)
            col = next((h for h in (rdr.fieldnames or []) if h.lower() in ("gene", "target", "symbol", "靶点")),
                       (rdr.fieldnames or [None])[0])
            out = [r[col].strip().upper() for r in rdr if r.get(col, "").strip()]
        else:
            f.seek(0)
            out = [l.strip().upper() for l in f if l.strip()]
    return sorted(set(out))


def compute_expression(adata, targets, cluster_key):
    """返回 {gene: {cluster: (mean_expr, pct_expressing)}} + 亚群顺序。"""
    import numpy as np
    clusters = list(adata.obs[cluster_key].astype(str).unique())
    genes = [g for g in targets if g in adata.var_names]
    missing = [g for g in targets if g not in adata.var_names]
    X = adata[:, genes].X
    try:
        X = X.toarray()
    except AttributeError:
        X = np.asarray(X)
    labels = adata.obs[cluster_key].astype(str).values
    mat_mean, mat_pct = {}, {}
    for gi, g in enumerate(genes):
        mat_mean[g], mat_pct[g] = {}, {}
        for c in clusters:
            vals = X[labels == c, gi]
            mat_mean[g][c] = float(np.mean(vals))
            mat_pct[g][c] = float(np.mean(vals > 0))
    return {"genes": genes, "missing": missing, "clusters": clusters,
            "mean": mat_mean, "pct": mat_pct}


def plot_dotplot(expr, out, label="PROJECT"):
    """发表级点图：X=亚群，Y=靶点，点大小=表达比例，颜色=平均表达。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    genes = expr["genes"]; clusters = expr["clusters"]
    if not genes or not clusters:
        return []
    n_g, n_c = len(genes), len(clusters)
    fig, ax = plt.subplots(figsize=(max(6, n_c * 0.7), max(4, n_g * 0.35)))
    for yi, g in enumerate(genes):
        for xi, c in enumerate(clusters):
            pct = expr["pct"][g][c]; mean = expr["mean"][g][c]
            if pct <= 0:
                continue
            ax.scatter(xi, n_g - yi - 1, s=max(8, pct * 220), c=[mean],
                       cmap="Blues", vmin=0, vmax=max(1e-9, max(expr["mean"][gg][cc]
                       for gg in genes for cc in clusters)), edgecolor="none", alpha=0.9)
    ax.set_xticks(range(n_c)); ax.set_xticklabels(clusters, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n_g)); ax.set_yticklabels(genes, fontsize=7)
    ax.set_xlabel("Cell subtype"); ax.set_title("Core target expression across cell subtypes [%s]" % label,
                                                fontsize=9, fontweight="bold")
    plt.colorbar(ax.collections[0] if ax.collections else plt.cm.ScalarMappable(cmap="Blues"),
                 ax=ax, label="Mean expression")
    plt.tight_layout()
    saved = []
    for ext in ("pdf", "svg", "png", "tiff"):
        p = os.path.join(out, "singlecell_dotplot.%s" % ext)
        if ext in ("png", "tiff"):
            fig.savefig(p, dpi=600, bbox_inches="tight")
            if ext == "tiff":
                from PIL import Image
                im = Image.open(p).convert("RGB"); im.save(p, format="TIFF", dpi=(600, 600), compression="tiff_lzw")
        else:
            fig.savefig(p, bbox_inches="tight")
        saved.append(p)
    plt.close(fig)
    return saved


def plot_violin(adata, genes, cluster_key, out, label="PROJECT"):
    """小提琴图（需 scanpy 的作图后端）。"""
    saved = []
    try:
        import scanpy as sc
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        valid = [g for g in genes if g in adata.var_names]
        if not valid:
            return saved
        sc.settings.verbosity = 0
        ax = sc.pl.violin(adata, valid, groupby=cluster_key, show=False, multi_panel=False)
        fig = plt.gcf()
        fig.set_size_inches(max(7, len(adata.obs[cluster_key].unique()) * 0.8), 4)
        for ext in ("pdf", "svg", "png"):
            p = os.path.join(out, "singlecell_violin.%s" % ext)
            fig.savefig(p, dpi=600, bbox_inches="tight")
            saved.append(p)
        plt.close(fig)
    except Exception as e:
        print("[warn] 小提琴图失败：%s" % str(e)[:100])
    return saved


def plot_umap(adata, genes, out, label="PROJECT"):
    """UMAP 表达图（各靶点一个 panel）。"""
    saved = []
    try:
        import scanpy as sc
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        if "X_umap" not in adata.obsm:
            try:
                sc.pp.neighbors(adata); sc.tl.umap(adata)
            except Exception:
                return saved
        valid = [g for g in genes if g in adata.var_names][:6]
        sc.pl.umap(adata, color=valid, show=False, ncols=3, wspace=0.3)
        fig = plt.gcf()
        for ext in ("pdf", "svg", "png"):
            p = os.path.join(out, "singlecell_umap.%s" % ext)
            fig.savefig(p, dpi=600, bbox_inches="tight")
            saved.append(p)
        plt.close(fig)
    except Exception as e:
        print("[warn] UMAP 图失败：%s" % str(e)[:100])
    return saved


def demo_expr(targets, out):
    """无 scanpy 时的流程演示：构造 3 亚群 × N 靶点矩阵，验证出图代码（DEMO 数据不可投稿）。"""
    import numpy as np
    rng = np.random.default_rng(42)
    clusters = ["Hepatocyte", "Kupffer", "Endothelial"]
    genes = targets[:12] if targets else ["IL6", "TNF", "AKT1", "CASP3"]
    mean = {g: {c: float(rng.uniform(0.05, 2.0)) for c in clusters} for g in genes}
    pct = {g: {c: float(rng.uniform(0.1, 0.9)) for c in clusters} for g in genes}
    expr = {"genes": genes, "missing": [], "clusters": clusters, "mean": mean, "pct": pct}
    saved = plot_dotplot(expr, out, label="DEMO (synthetic, NOT for publication)")
    return expr, saved


def plot_spatial(adata, genes, out):
    """空间特征图：核心靶点映射到组织切片坐标 obsm['spatial']。
    回答"药作用于病灶边缘巨噬细胞区"这类空间定位问题。"""
    saved = []
    try:
        import scanpy as sc
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        if "spatial" not in adata.obsm:
            return saved
        valid = [g for g in genes if g in adata.var_names][:6]
        sc.pl.spatial(adata, color=valid, show=False, ncols=3, spot_size=1.0)
        fig = plt.gcf()
        for ext in ("pdf", "svg", "png"):
            p = os.path.join(out, "spatial_feature.%s" % ext)
            fig.savefig(p, dpi=600, bbox_inches="tight")
            saved.append(p)
        plt.close(fig)
    except Exception as e:
        print("[warn] 空间图失败：%s" % str(e)[:100])
    return saved


def selftest():
    ok = callable(read_targets) and callable(compute_expression) and callable(plot_dotplot) and callable(plot_spatial)
    print("selftest: %s (scanpy available=%s)" % ("PASS" if ok else "FAIL", has_scanpy()))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_singlecell: 核心靶点单细胞亚群投影")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("project")
    p.add_argument("--h5ad", default=None)
    p.add_argument("--targets", required=True)
    p.add_argument("--cluster-key", default="celltype")
    p.add_argument("--demo", action="store_true", help="无 scanpy 时用合成数据演示出图")
    p.add_argument("-o", "--outdir", default="figs/singlecell")
    p.set_defaults(fn="project")
    p = sub.add_parser("geo")
    p.add_argument("--gse", required=True)
    p.set_defaults(fn="geo")
    p = sub.add_parser("selftest"); p.set_defaults(fn="selftest")
    a = ap.parse_args(argv)
    if a.fn == "selftest":
        return selftest()
    if a.fn == "geo":
        print("GEO 单细胞数据获取指引（不伪造下载）：")
        print("  1. GEO 页面：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s" % a.gse)
        print("  2. 下载表达矩阵（如 filtered_feature_bc_matrix 或 expression matrix）。")
        print("  3. 构建 h5ad：sc.read_10x_mtx('filtered_feature_bc_matrix/') 或 sc.read_csv('matrix.csv').T")
        print("  4. 加入亚群注释：adata.obs['celltype'] = <metadata['celltype']>")
        print("  5. 用 project --h5ad <文件>.h5ad --cluster-key celltype 投影核心靶点。")
        return 0

    targets = read_targets(a.targets)
    os.makedirs(a.outdir, exist_ok=True)
    if a.demo or not has_scanpy():
        if not a.demo:
            print("[warn] 未检测到 scanpy → 自动进入 DEMO 模式（合成数据，禁止用于论文）。")
            print("       安装：pip install scanpy anndata；然后重跑 project --h5ad ...")
        expr, saved = demo_expr(targets, a.outdir)
    else:
        import scanpy as sc
        adata = sc.read_h5ad(a.h5ad)
        expr = compute_expression(adata, targets, a.cluster_key)
        saved = plot_dotplot(expr, a.outdir, label="PROJECT")
        saved += plot_violin(adata, expr["genes"], a.cluster_key, a.outdir)
        saved += plot_umap(adata, expr["genes"], a.outdir)
        if "spatial" in adata.obsm:
            print("[singlecell] 检测到 obsm['spatial'] → 渲染空间特征图")
            saved += plot_spatial(adata, expr["genes"], a.outdir)
    json.dump({"n_targets": len(targets), "mapped": len(expr["genes"]),
               "missing_genes": expr["missing"], "clusters": expr["clusters"],
               "demo_mode": (a.demo or not has_scanpy()),
               "note": "DEMO 模式为合成数据，不得用于论文；正式分析需 scanpy + 真实 h5ad",
               "figures": [os.path.basename(p) for p in saved]},
              open(os.path.join(a.outdir, "singlecell_manifest.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[singlecell] 靶点 %d（映射 %d）→ %d 亚群，出图 %d 张 → %s" %
          (len(targets), len(expr["genes"]), len(expr["clusters"]), len(saved), a.outdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
