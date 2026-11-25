#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_export.py — 标准化中间产物导出（可追溯的核心痛点）

无论输入是什么格式，统一输出 5 张标准表：
  herbs.csv               药味表（药名/来源库/是否命中/降级策略）
  compounds.csv           成分表（成分/来源库/筛选阈值/OB/DL/是否保留）
  targets_mapped.csv      靶点映射表（原始名→UniProt ID→HGNC 符号→物种→映射层级）
  PPI_network.sif         Cytoscape 可导入的 PPI 网络（node1 interaction node2 score）
  enrichment_results.csv  富集表（term/source/pvalue/qvalue/n_hit/term_size/ID）

子命令
  all       从现有 skill 产物自动归集 5 张表
  check     检查标准产物是否齐全、列是否合规
"""
from __future__ import annotations

import argparse, csv, os, sys

STANDARD = {
    "herbs.csv": ["herb", "source_db", "matched", "fallback_strategy"],
    "compounds.csv": ["compound", "source_db", "filter_rule", "ob", "dl", "kept"],
    "targets_mapped.csv": ["raw", "uniprot_id", "symbol", "organism", "layer"],
    "PPI_network.sif": ["node1", "interaction", "node2", "score"],
    "enrichment_results.csv": ["term", "source", "pvalue", "qvalue", "n_hit", "term_size", "term_id"],
}


def _read(path, delim=None):
    if not os.path.exists(path):
        return []
    d = delim or ("\t" if path.endswith(".tsv") else ",")
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        return list(csv.DictReader(f, delimiter=d))


def _write(path, header, rows):
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow([r.get(h, "") for h in header])


def cmd_all(args):
    d = args.dir
    # 输入主目录与产物子目录都可传：自动向上找 compounds/STRING 原始表
    parent = os.path.dirname(os.path.abspath(d))
    def find_rel(*names):
        for base in (d, parent):
            for n in names:
                p = os.path.join(base, n)
                if os.path.exists(p):
                    return p
        return os.path.join(d, names[0])
    # compounds: 优先当前目录，其次父目录；再 compounds_kept
    comp = (_read(os.path.join(d, "compounds.csv"))
            or _read(os.path.join(parent, "compounds.csv"))
            or _read(find_rel("compounds_kept.csv")))
    comp_rows = []
    for r in comp:
        comp_rows.append({
            "compound": r.get("compound") or r.get("name") or "",
            "source_db": r.get("source_db") or r.get("cohort") or "unknown",
            "filter_rule": r.get("gate_used") or r.get("category") or "未记录",
            "ob": r.get("ob", ""), "dl": r.get("dl", ""),
            "kept": r.get("kept", "yes")})
    # 若输入已是 compounds.csv（父目录），当前目录无同名文件则直接写；否则写标准视图
    src_is_comp = comp and _read(os.path.join(d, "compounds.csv"))
    out_comp = os.path.join(d, "compounds_standard.csv" if src_is_comp else "compounds.csv")
    _write(out_comp, STANDARD["compounds.csv"], comp_rows)

    # targets
    tg = _read(os.path.join(d, "targets_standardized.csv"))
    tg_rows = [{
        "raw": r.get("raw", ""), "uniprot_id": r.get("uniprot_id") or r.get("id", ""),
        "symbol": r.get("symbol", ""),
        "organism": r.get("organism", "Homo sapiens"),
        "layer": r.get("layer", "")} for r in tg]
    _write(os.path.join(d, "targets_mapped.csv"), STANDARD["targets_mapped.csv"], tg_rows)

    # PPI sif
    edges = _read(find_rel("string_edges.tsv")) or _read(find_rel("intersect_edges.tsv"))
    sif = [{"node1": r.get("node1") or r.get("source", ""),
            "interaction": "pp",
            "node2": r.get("node2") or r.get("target", ""),
            "score": r.get("score", "")} for r in edges]
    _write(os.path.join(d, "PPI_network.sif"), STANDARD["PPI_network.sif"], sif)

    # enrichment
    en = _read(os.path.join(d, "enrich_significant.csv")) or _read(os.path.join(d, "enrich_all.csv"))
    en_rows = [{
        "term": r.get("term", ""), "source": r.get("source", ""),
        "pvalue": r.get("pvalue", ""), "qvalue": r.get("qvalue", ""),
        "n_hit": r.get("n_hit", ""), "term_size": r.get("term_size", ""),
        "term_id": r.get("term_id") or r.get("ID", "")} for r in en]
    _write(os.path.join(d, "enrichment_results.csv"), STANDARD["enrichment_results.csv"], en_rows)

    # herbs（从 compounds 的 cohort 推断，无则空表+提示）
    herbs = sorted(set(r["source_db"] for r in comp_rows if r["source_db"] and r["source_db"] != "unknown"))
    herb_rows = [{"herb": h, "source_db": "TCMSP", "matched": "yes",
                  "fallback_strategy": "见 provenance.json"} for h in herbs]
    _write(os.path.join(d, "herbs.csv"), STANDARD["herbs.csv"], herb_rows)

    print(f"标准产物已导出到 {d}")
    for fn in STANDARD:
        p = os.path.join(d, fn)
        n = sum(1 for _ in open(p, encoding="utf-8-sig")) - 1 if os.path.exists(p) else -1
        print(f"  {fn:<28} {n:>6} 行" if n >= 0 else f"  {fn:<28} 缺失")


def cmd_check(args):
    d = args.dir
    ok = True
    print(f"=== 标准产物检查: {d} ===")
    for fn, header in STANDARD.items():
        p = os.path.join(d, fn)
        alt = os.path.join(d, "compounds_standard.csv") if fn == "compounds.csv" else ""
        if not os.path.exists(p) and not (alt and os.path.exists(alt)):
            print(f"  [缺失] {fn}"); ok = False; continue
        if not os.path.exists(p):
            p = alt
        with open(p, encoding="utf-8-sig") as f:
            got = next(csv.reader(f), [])
        missing = [h for h in header if h not in got]
        if missing:
            print(f"  [列不全] {fn} 缺 {missing}"); ok = False
        else:
            print(f"  [OK] {fn} ({len(got)} 列)")
    return 0 if ok else 1


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("all"); a.add_argument("--dir", required=True); a.set_defaults(fn=cmd_all)
    c = sub.add_parser("check"); c.add_argument("--dir", required=True); c.set_defaults(fn=cmd_check)
    args = p.parse_args(); sys.exit(args.fn(args) or 0)


if __name__ == "__main__":
    main()
