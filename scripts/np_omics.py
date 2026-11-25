#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_omics.py —— v1.5.0 多组学交叉映射模块（GEO DEGs × NP 靶点）

背景（23 篇审计）：GEO DEGs 1/23（antrodia 用 4 个 GSE+RRA）、scRNA 1/23、WGCNA 0/23。
高频论文的套路是"NP 预测靶点 ∩ 疾病 DEGs → 交集基因进富集/对接"。

命令：
  deg    --gse GSE112790 [--outdir DIR] [--logfc 1] [--fdr 0.05]
         下载 GEO series matrix → 解析表达矩阵 → 分组 t 检验 → BH-FDR → DEGs.csv
  cross  --gse GSE112790 --targets targets.txt [--outdir DIR]
         DEGs ∩ NP 预测靶点 → 交集清单 + Venn 数据（可直接喂给 np_enrich）
  scan   --gse GSE112790
         只列出该 GSE 的样本分组信息（辅助确认 disease/control 列）

GEO 访问：https://www.ncbi.nlm.nih.gov/geo/download/?acc=<GSE>&format=file&file=<GSE>_series_matrix.txt.gz
（2026-11-20 实测 api_ok；NCBI eutils 亦可查 GSE 元数据）
"""
import argparse
import gzip
import io
import os
import re
import sys

import numpy as np
import pandas as pd

GEO_FILE = "https://ftp.ncbi.nlm.nih.gov/geo/series/{head}nnn/{acc}/matrix/{fn}"


def _gpl_of(acc):
    """NCBI eutils gds 库查 GSE 的平台 GPL 编号"""
    import requests, time
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    r = requests.get(base + "esearch.fcgi", params={"db": "gds", "term": f"{acc}[ACCN]", "retmode": "json"}, timeout=60)
    ids = (r.json().get("esearchresult") or {}).get("idlist") or []
    if not ids:
        return None
    time.sleep(0.4)
    r = requests.get(base + "esummary.fcgi", params={"db": "gds", "id": ids[0], "retmode": "json"}, timeout=60)
    doc = (r.json().get("result") or {}).get(ids[0], {})
    return doc.get("gpl")  # 如 GPL571


def fetch_series_matrix(gse, cachedir):
    """下载（带本地缓存）并解压 series matrix；直接拼文件名（目录列表偶发404不走）"""
    os.makedirs(cachedir, exist_ok=True)
    gz = os.path.join(cachedir, f"{gse}_series_matrix.txt.gz")
    if not os.path.exists(gz):
        import requests, time
        acc = gse.upper()
        digits = acc[3:]
        head = digits[:-3]  # 模板再补 nnn（112790→112, 14520→14）
        cands = [f"{acc}_series_matrix.txt.gz"]
        gpl = _gpl_of(acc)
        if gpl:
            for p in re.split(r"[;,]", str(gpl)):
                p = p.strip().upper()
                if p and not p.startswith("GPL"):
                    p = "GPL" + p
                if re.match(r"^GPL\d+$", p):
                    cands.insert(0, f"{acc}-{p}_series_matrix.txt.gz")
        got = None
        for fn in cands:
            u = GEO_FILE.format(head=head, acc=acc, fn=fn)
            for k in range(3):
                d = requests.get(u, timeout=600, headers={"User-Agent": "np-skill/1.5"})
                if d.status_code == 200 and len(d.content) > 500:
                    got = d.content
                    print(f"下载 {fn}（{len(d.content)//1024} KB）", flush=True)
                    break
                time.sleep(3 * (k + 1))
            if got:
                break
        if not got:
            raise SystemExit(f"[FAIL] {acc} series matrix 下载失败（试过 {cands}）；"
                             f"RNA-seq 类 GSE 无此文件，请改用芯片类 GSE 或手工下载 suppl 文件")
        open(gz, "wb").write(got)
    with gzip.open(gz, "rt", encoding="utf-8", errors="ignore") as f:
        return f.read()


def parse_series_matrix(txt):
    """解析 !Sample_* 元数据 + 值矩阵；返回 expr DataFrame + 特征表"""
    chars, platforms, titles, geos = {}, [], [], []
    expr = None
    for line in txt.splitlines():
        line = line.strip()
        if expr is None and line.startswith("!Sample_characteristics_ch1"):
            parts = line.split("\t")
            chars[len(chars)] = [p.strip().strip('"') for p in parts[1:]]
        elif line.startswith("!Sample_title"):
            titles = [p.strip().strip('"') for p in line.split("\t")[1:]]
        elif line.startswith("!Sample_geo_accession"):
            geos = [p.strip().strip('"') for p in line.split("\t")[1:]]
        elif line.startswith("!series_matrix_table_begin"):
            mode = "table"; continue
        elif line.startswith("!series_matrix_table_end"):
            mode = None; continue
        if line.lstrip('"').startswith("ID_REF") and expr is None:
            hdr = line.replace('"', "").split("\t")
            cols = hdr[1:]
            expr = {"cols": cols, "rows": []}
        elif expr is not None and mode == "table" and not line.startswith("!"):
            parts = line.replace('"', "").split("\t")
            if len(parts) == len(expr["cols"]) + 1:
                expr["rows"].append((parts[0], [p for p in parts[1:]]))
    if expr is None or not expr["rows"]:
        raise SystemExit("[FAIL] series matrix 中无表达矩阵（该 GSE 可能只有元数据或用 suppl 文件）")
    mat = pd.DataFrame([r[1] for r in expr["rows"]], index=[r[0] for r in expr["rows"]],
                       columns=expr["cols"], dtype=float)
    meta = pd.DataFrame(chars, index=geos[:len(chars[0])] if geos else range(len(chars[0])))
    meta.index.name = "sample"
    return mat, meta, titles


def guess_groups(meta):
    """从 characteristics 行猜测 disease/control：找含 disease/normal/control 关键词的行"""
    lab = {}
    for i, row in meta.iterrows():
        lab[i] = [str(v).strip() for v in row.values]
    disease_kw = ["disease", "tumor", "cancer", "carcinoma", "hcc", "patient", "case", "模型", " diabetic", "injur"]
    control_kw = ["normal", "control", "healthy", "adjacent", "non-tumor", "对照", "健康"]
    best = None
    for ridx in range(meta.shape[1]):
        vals = [str(v).lower() for v in meta.iloc[:, ridx].values]
        n_d = sum(any(k in v for k in disease_kw) for v in vals)
        n_c = sum(any(k in v for k in control_kw) for v in vals)
        if n_d >= 2 and n_c >= 2:
            best = ridx
            break
    if best is None:
        return None
    col = meta.iloc[:, best]
    g1, g2 = [], []
    for s, v in col.items():
        lv = str(v).lower()
        # 先查 control 关键词（non-tumor/adjacent 含 tumor 字样，必须先判）
        if any(k in lv for k in control_kw):
            g2.append(s)
        elif any(k in lv for k in disease_kw):
            g1.append(s)
    return g1, g2


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(ranked, 0, 1)
    return out


def cmd_deg(args):
    txt = fetch_series_matrix(args.gse, args.outdir)
    mat, meta, titles = parse_series_matrix(txt)
    print(f"矩阵: {mat.shape[0]} 探针 × {mat.shape[1]} 样本")
    g = guess_groups(meta)
    if not g:
        print("\n自动分组失败。样本特征表：")
        print(meta.to_string()[:2000])
        raise SystemExit("[提示] 请人工确认分组后用 --disease/--control 指定样本 GSM 编号")
    g1, g2 = g
    if args.disease: g1 = args.disease.split(",")
    if args.control: g2 = args.control.split(",")
    print(f"分组: disease={len(g1)} vs control={len(g2)}")
    from scipy import stats
    X1 = mat[g1].values; X2 = mat[g2].values
    t, p = stats.ttest_ind(X1, X2, axis=1)
    fc = X1.mean(axis=1) - X2.mean(axis=1)  # 已是 log2 尺度则即 log2FC
    fdr = bh_fdr(p)
    res = pd.DataFrame({"probe": mat.index, "log2FC": fc, "t": t, "p": p, "fdr": fdr})
    res = res.sort_values("p")
    res.to_csv(os.path.join(args.outdir, f"{args.gse}_degs_all.csv"), index=False, encoding="utf-8-sig")
    sig = res[(res.fdr <= args.fdr) & (res.log2FC.abs() >= args.logfc)]
    sig.to_csv(os.path.join(args.outdir, f"{args.gse}_degs.csv"), index=False, encoding="utf-8-sig")
    print(f"DEGs（|log2FC|≥{args.logfc}, FDR≤{args.fdr}）: {len(sig)} 条 → {args.outdir}/{args.gse}_degs.csv")


GPL_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL{m1}nnn/{gpl}/annot/{gpl}.annot.gz"


def probe2gene(gpl, cachedir):
    """下载 NCBI 官方平台注释 → probe→gene 映射 dict（多基因探针取第一个）"""
    os.makedirs(cachedir, exist_ok=True)
    gpl = gpl.upper().replace("GPL", "GPL")
    n = gpl[3:]
    fn = os.path.join(cachedir, f"{gpl}_probe2gene.csv")
    if os.path.exists(fn):
        d = pd.read_csv(fn, encoding="utf-8-sig", dtype=str)
        return dict(zip(d.probe, d.gene))
    u = GPL_ANNOT.format(m1=n[:3] if len(n) > 3 else "", gpl=gpl)
    import requests, gzip as _gz, time
    got = None
    for k in range(3):
        r = requests.get(u, timeout=300, headers={"User-Agent": "np-skill/1.5"})
        if r.status_code == 200 and len(r.content) > 1000:
            got = r.content; break
        time.sleep(3 * (k + 1))
    if not got:
        raise SystemExit(f"[FAIL] 平台注释下载失败：{u}")
    t = _gz.decompress(got).decode("utf-8", "ignore")
    i = t.find("!platform_table_begin")
    lines = t[i:].splitlines()
    hdr = lines[1].split("\t")
    si = hdr.index("Gene symbol") if "Gene symbol" in hdr else None
    if si is None:
        raise SystemExit("[FAIL] 注释无 Gene symbol 列")
    rows = []
    for l in lines[2:]:
        if l.startswith("!platform_table_end"):
            break
        p = l.split("\t")
        if len(p) > si and p[si].strip():
            gene = p[si].split("///")[0].strip()
            if gene:
                rows.append((p[0], gene))
    pd.DataFrame(rows, columns=["probe", "gene"]).to_csv(fn, index=False, encoding="utf-8-sig")
    print(f"平台注释 {gpl}: {len(rows)} 探针映射 → {fn}")
    return dict(rows)


def platform_of(txt):
    m = re.search(r'!Series_platform_id\s+"?(GPL\d+)', txt)
    return m.group(1) if m else None


def cmd_scan(args):
    txt = fetch_series_matrix(args.gse, args.outdir)
    mat, meta, titles = parse_series_matrix(txt)
    print(f"GSE: {args.gse} | 样本数: {mat.shape[1]}")
    for i in range(meta.shape[1]):
        vals = pd.Series(meta.iloc[:, i].values).value_counts()
        print(f"特征{i}: {dict(list(vals.items())[:6])}")


def cmd_cross(args):
    txt = fetch_series_matrix(args.gse, args.outdir)
    degs = os.path.join(args.outdir, f"{args.gse}_degs.csv")
    if not os.path.exists(degs):
        raise SystemExit("[FAIL] 先运行 deg 子命令")
    d = pd.read_csv(degs, encoding="utf-8-sig")
    targets = set(open(args.targets, encoding="utf-8-sig").read().split())
    targets = {t.upper() for t in targets if t.strip()}
    # 探针→基因：DEGs 表自带 gene 列，或从 NCBI 平台注释自动映射
    if "gene" not in d.columns:
        gpl = platform_of(txt)
        if not gpl:
            raise SystemExit("[FAIL] 无法确定平台号（!Series_platform_id 缺失）")
        print(f"自动注释平台 {gpl} ...")
        p2g = probe2gene(gpl, args.outdir)
        d["gene"] = d.probe.astype(str).map(p2g)
    genes = set(d.gene.dropna().str.upper())
    inter = sorted(genes & targets)
    up = set(d[d.log2FC > 0].gene.str.upper()); down = set(d[d.log2FC < 0].gene.str.upper())
    rows = [{"gene": g, "in_np_targets": True,
             "direction": "up" if g in up else ("down" if g in down else "ambiguous")} for g in inter]
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(args.outdir, f"{args.gse}_np_intersection.csv"), index=False, encoding="utf-8-sig")
    print(f"NP 靶点 {len(targets)} ∩ DEGs {len(genes)} = {len(inter)} 个")
    print(f"交集 → {args.outdir}/{args.gse}_np_intersection.csv")
    print("基因:", ", ".join(inter[:40]))


def main():
    ap = argparse.ArgumentParser(description="np_omics: GEO DEGs × NP 靶点交叉映射")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("deg")
    s.add_argument("--gse", required=True); s.add_argument("--outdir", default=".")
    s.add_argument("--logfc", type=float, default=1.0); s.add_argument("--fdr", type=float, default=0.05)
    s.add_argument("--disease", default=None); s.add_argument("--control", default=None)
    s = sub.add_parser("cross")
    s.add_argument("--gse", required=True); s.add_argument("--targets", required=True)
    s.add_argument("--outdir", default=".")
    s = sub.add_parser("scan")
    s.add_argument("--gse", required=True); s.add_argument("--outdir", default=".")
    a = ap.parse_args()
    {"deg": cmd_deg, "cross": cmd_cross, "scan": cmd_scan}[a.cmd](a)


if __name__ == "__main__":
    main()
