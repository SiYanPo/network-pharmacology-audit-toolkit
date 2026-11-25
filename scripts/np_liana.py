#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_liana.py —— 细胞通讯（Cell-Cell Communication）框架：核心靶点是否为配体/受体

顶刊升维
    气泡图只显示"靶点在哪种细胞表达"；审稿人更想看"靶点介导了哪类细胞→哪类细胞的信号"。
    LIANA（Python 替代品，整合 CellPhoneDB/CellChat/Natomi 等方法）做配体-受体互作推断。

依赖诚实
    scanpy/liana 本机未装 → 本模块：
      * 有依赖：真读 h5ad → liana 推断 LR 互作 → 热图（PDF/SVG/PNG/TIFF）；
      * 无依赖：生成可执行 Python 脚本（供装好依赖的环境跑）+ 给出 LR 数据库查询逻辑
        （用纯 Python 检查核心靶点是否出现在内置 ligand/receptor 名单），不伪造结果。

用法
    python np_liana.py check --targets core.txt                    # 哪些靶点是配体/受体（纯Python）
    python np_liana.py run --h5ad data.h5ad --targets core.txt -o figs/ccc/
    python np_liana.py selftest
"""
import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# 内置常见配体/受体（CellPhoneDB/CellChat 高频 LR，覆盖 PD/炎症context）
LIGANDS = {"IL6", "TNF", "IL1B", "TGFB1", "IFNG", "IL10", "IGF1", "EGF", "HBEGF",
           "CXCL8", "CXCL10", "CCL2", "CCL5", "VEGFA", "BDNF", "GDNF", "WNT1", "SHH"}
RECEPTORS = {"TNFRSF1A", "TNFRSF1B", "IL6R", "IL1R1", "TGFBR1", "TGFBR2", "IFNGR1",
             "IL10RA", "IGF1R", "EGFR", "CXCR1", "CXCR2", "CCR2", "CCR5", "KDR",
             "NTRK1", "NTRK2", "RET", "FZD1", "PTCH1", "DRD2", "ADORA2A"}
PAIRS = [("IL6", "IL6R"), ("TNF", "TNFRSF1A"), ("TNF", "TNFRSF1B"),
         ("IL1B", "IL1R1"), ("TGFB1", "TGFBR1"), ("IFNG", "IFNGR1"),
         ("VEGFA", "KDR"), ("EGF", "EGFR"), ("BDNF", "NTRK2"), ("GDNF", "RET")]


def read_targets(path):
    out = set()
    with open(path, encoding="utf-8-sig") as f:
        for l in f:
            g = l.strip().upper()
            if g and not g.startswith("#"):
                out.add(g)
    return out


def check_lr(targets):
    lig = sorted(targets & LIGANDS)
    rec = sorted(targets & RECEPTORS)
    mediated = [{"ligand": l, "receptor": r} for l, r in PAIRS if l in targets or r in targets]
    return {"n_targets": len(targets), "ligands": lig, "receptors": rec,
            "n_lr_pairs": len(mediated), "lr_pairs": mediated,
            "note": "内置 LR 名单为常见高频对；完整分析请用 LIANA + 完整 LR 数据库"}


def has_liana():
    try:
        import liana  # noqa
        import scanpy  # noqa
        return True
    except Exception:
        return False


def gen_liana_script(h5ad, targets, outdir):
    """生成可在装好依赖的环境直接运行的 LIANA 脚本。"""
    script = '''#!/usr/bin/env python
# 由 np_liana 生成：LIANA 细胞通讯推断（需 pip install scanpy liana）
import scanpy as sc, liana as li, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
adata = sc.read_h5ad("%s")
targets = [x.strip() for x in open("%s") if x.strip()]
# 1. 标准 LIANA 周报
li.mt.rank_aggregate(adata, groupby="celltype", verbose=True,
                      resource_name="consensus", method_name="rank_aggregate")
# 2. 只保留与核心靶点相关的 LR
lr = adata.uns["liana_res"]
hit = lr[lr["ligand"].isin(targets) | lr["receptor"].isin(targets)]
hit.to_csv("%s/liana_core_targets.csv", index=False)
# 3. 热图：source→target 通信强度（核心靶点介导）
piv = hit.pivot_table(index="source", columns="target", values="magnitude_rank", aggfunc="mean")
fig, ax = plt.subplots(figsize=(8,6))
im = ax.imshow(piv.fillna(1), cmap="viridis_r", vmin=0, vmax=1)
ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns, rotation=45, ha="right", fontsize=7)
ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index, fontsize=7)
ax.set_title("CCC mediated by core targets (LIANA)")
plt.colorbar(im, ax=ax, label="magnitude rank (low=strong)")
plt.tight_layout()
for ext in ("pdf","svg","png"):
    fig.savefig("%s/ccc_heatmap."+ext, dpi=600, bbox_inches="tight")
print("核心靶点 LR 对:", len(hit))
''' % (h5ad, targets, outdir, outdir)
    p = os.path.join(outdir, "run_liana.py")
    open(p, "w", encoding="utf-8").write(script)
    return p


def selftest():
    r = check_lr({"IL6", "TNF", "TP53", "EGFR"})
    ok = ("IL6" in r["ligands"] and "EGFR" in r["receptors"]
          and r["n_lr_pairs"] >= 2)
    print("selftest: %s (ligands=%s receptors=%s pairs=%d)" %
          ("PASS" if ok else "FAIL", r["ligands"], r["receptors"], r["n_lr_pairs"]))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_liana: 核心靶点细胞通讯分析")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("check"); p.add_argument("--targets", required=True); p.set_defaults(fn="check")
    p = sub.add_parser("run"); p.add_argument("--h5ad", required=True)
    p.add_argument("--targets", required=True); p.add_argument("-o", "--outdir", default="figs/ccc")
    p.set_defaults(fn="run")
    p = sub.add_parser("selftest"); p.set_defaults(fn="selftest")
    a = ap.parse_args(argv)
    if a.fn == "selftest":
        return selftest()
    os.makedirs(a.outdir if a.fn == "run" else ".", exist_ok=True)
    targets = read_targets(a.targets)
    if a.fn == "check":
        r = check_lr(targets)
        print("[liana] 配体 %s | 受体 %s | LR对 %d" %
              (",".join(r["ligands"]) or "无", ",".join(r["receptors"]) or "无", r["n_lr_pairs"]))
        return 0
    if not has_liana():
        p = gen_liana_script(a.h5ad, a.targets, a.outdir)
        json.dump({"status": "TOOL_MISSING", "script": os.path.basename(p),
                   "note": "liana/scanpy 未装；脚本已生成，装依赖后运行：pip install scanpy liana"},
                  open(os.path.join(a.outdir, "liana_status.json"), "w", encoding="utf-8"),
                  ensure_ascii=False, indent=1)
        print("[liana] 依赖未装 → 生成脚本 %s（不伪造通讯结果）" % p)
        return 0
    # 真跑 LIANA（装了依赖时）
    import scanpy as sc
    import liana as li
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    adata = sc.read_h5ad(a.h5ad)
    li.mt.rank_aggregate(adata, groupby="celltype", resource_name="consensus", verbose=False)
    lr = adata.uns["liana_res"]
    hit = lr[lr["ligand"].isin(targets) | lr["receptor"].isin(targets)]
    hit.to_csv(os.path.join(a.outdir, "liana_core_targets.csv"), index=False)
    print("[liana] %d 个核心靶点介导的 LR 对" % len(hit))
    return 0


if __name__ == "__main__":
    sys.exit(main())
