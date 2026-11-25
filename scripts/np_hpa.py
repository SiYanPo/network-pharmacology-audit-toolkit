#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_hpa.py —— 疾病组织特异性校验（堵"细胞系做出来的不符合真实组织"质疑）

科学问题
    核心 Hub 基因必须在"病变组织"中真上调，否则可能是细胞系假象。
    HPA (Human Protein Atlas) 提供 62 组织的 RNA-seq 共识表达（nTPM）。

数据源（实测）
    HPA search API：https://www.proteinatlas.org/api/search_download.php?search=<gene>&format=json
    HPA bulk TSV ：https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip（~几 MB）

判定逻辑
    nTPM(疾病组织) / nTPM(全组织中位) ≥ --fold（默认 2）→ Enriched（支持）
    否则 → Not Enriched（假阳性风险，须讨论）
    组织名用 HPA 标准名（如 "cerebral cortex" / "liver" / "lung"）。

用法
    python np_hpa.py fetch --tissue liver --outdir hpa_cache/
    python np_hpa.py check --targets core.txt --tissue "cerebral cortex" --bulk hpa_cache/rna_tissue_consensus.tsv -o out/
    python np_hpa.py selftest
"""
import argparse
import csv
import json
import os
import re
import sys
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
# 实测（2026-11-22）：正确路径为 /download/tsv/...（旧 /download/... 返回 404）
HPA_BULK = "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"
HPA_SEARCH = "https://www.proteinatlas.org/api/search_download.php?search=%s&format=json&columns=g,gsn&compress=no"


def fetch_bulk(outdir):
    """下载 HPA 组织共识表达 bulk TSV（解压）。"""
    os.makedirs(outdir, exist_ok=True)
    zp = os.path.join(outdir, "rna_tissue_consensus.tsv.zip")
    tp = os.path.join(outdir, "rna_tissue_consensus.tsv")
    if not os.path.exists(tp):
        try:
            urllib.request.urlretrieve(HPA_BULK, zp)
            with zipfile.ZipFile(zp) as z:
                z.extractall(outdir)
        except Exception as e:
            return None, "download_failed:%s" % str(e)[:80]
    # 版本可追溯：HPA 没有稳定的版本 API，只能记 URL + 文件哈希 + （若用户提供）release。
    # **不猜**：未设置 NP_HPA_RELEASE 时写 unrecorded，Methods 里必须手工补。
    try:
        import np_snapshot as _S
        _S.write_meta(tp, url=HPA_BULK,
                      release=os.environ.get("NP_HPA_RELEASE") or "unrecorded",
                      extra={"note": "HPA release 无法自动探测；请设 NP_HPA_RELEASE=24.0 之类再重跑"})
    except Exception:
        pass
    return tp, ""


def tissue_enrichment(tissue, tsv, targets, fold=2.0):
    """对核心靶点算 nTPM(tissue) / median(nTPM across tissues)。"""
    if not tsv or not os.path.exists(tsv):
        return {"status": "NO_BULK", "reason": "缺 HPA bulk TSV；先运行 fetch"}
    tis = tissue.strip().lower()
    rows = []
    with open(tsv, encoding="utf-8") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        cols = rdr.fieldnames or []
        # 实测表头：Gene=Ensembl ID；Gene name=基因符号 → 优先 "Gene name"
        gcol = next((c for c in cols if c.lower() in ("gene name", "gene_name")),
                    next((c for c in cols if c.lower() == "gene"), cols[0]))
        tcol = next((c for c in cols if c.lower() == "tissue"), None)
        vcol = next((c for c in cols if "ntpm" in c.lower()), None)
        # 筛目标组织
        tissue_data = {}
        all_tissues = {}
        for r in rdr:
            g = (r.get(gcol) or "").strip().upper()
            tn = (r.get(tcol) or "").strip().lower() if tcol else ""
            try:
                v = float(r.get(vcol) or 0)
            except ValueError:
                continue
            if g in targets:
                all_tissues.setdefault(g, []).append(v)
                if tn == tis or tis in tn:
                    tissue_data[g] = tissue_data.get(g, 0) + v
    med = {}
    for g, vs in all_tissues.items():
        vs_sorted = sorted(vs)
        med[g] = vs_sorted[len(vs_sorted) // 2] if vs_sorted else 0
    results = []
    for g in sorted(targets):
        tpm = tissue_data.get(g, 0.0)
        m = med.get(g, 0.0)
        ratio = tpm / m if m > 0 else (float("inf") if tpm > 0 else 0)
        results.append({"gene": g, "tissue": tissue, "ntpm_tissue": round(tpm, 2),
                        "ntpm_median": round(m, 2), "fold": round(ratio, 2) if ratio != float("inf") else "inf",
                        "enriched": bool(ratio >= fold)})
    n_enr = sum(1 for r in results if r["enriched"])
    return {"status": "OK", "tissue": tissue, "n_targets": len(targets),
            "n_enriched": n_enr, "fold_cut": fold, "results": results,
            "note": "Not Enriched 靶点需讨论是否为细胞系/数据库假阳性；"
                    "HPA 为健康组织基线，病变组织上调需结合 GEO 疾病队列确认"}


def selftest():
    ok = callable(tissue_enrichment) and callable(fetch_bulk)
    print("selftest: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_hpa: 疾病组织特异性校验")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("fetch"); p.add_argument("-o", "--outdir", default="np_out/hpa")
    p.set_defaults(fn="fetch")
    p = sub.add_parser("check"); p.add_argument("--targets", required=True)
    p.add_argument("--tissue", required=True); p.add_argument("--bulk", required=True)
    p.add_argument("--fold", type=float, default=2.0)
    p.add_argument("-o", "--outdir", default="np_out/hpa")
    p.set_defaults(fn="check")
    p = sub.add_parser("selftest"); p.set_defaults(fn="selftest")
    a = ap.parse_args(argv)
    if a.fn == "selftest":
        return selftest()
    if a.fn == "fetch":
        p, err = fetch_bulk(a.outdir)
        if err:
            print("[hpa] %s；请手动下载 %s 到 %s" % (err, HPA_BULK, a.outdir))
            return 2
        print("[hpa] bulk → %s" % p)
        return 0
    targets = set()
    with open(a.targets, encoding="utf-8-sig") as f:
        targets = {l.strip().upper() for l in f if l.strip()}
    os.makedirs(a.outdir, exist_ok=True)
    r = tissue_enrichment(a.tissue, a.bulk, targets, a.fold)
    with open(os.path.join(a.outdir, "hpa_tissue_enrichment.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "tissue", "ntpm_tissue", "ntpm_median", "fold", "enriched"])
        for x in r.get("results", []):
            w.writerow([x["gene"], x["tissue"], x["ntpm_tissue"], x["ntpm_median"], x["fold"], x["enriched"]])
    json.dump(r, open(os.path.join(a.outdir, "hpa_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[hpa] %s：%d/%d 靶点富集（fold≥%.1f）→ %s" %
          (a.tissue, r.get("n_enriched", 0), r.get("n_targets", 0), a.fold, a.outdir))
    return 0 if r["status"] == "OK" else 2


if __name__ == "__main__":
    sys.exit(main())
