#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tcmsp_adapter.py —— TCMSP 原始输出 → pro 标准表（BaseAdapter 的具体示例）

演示三种契约差异的统一：
    TCMSP "成分-靶点" 表常见列：mol_id / molecule_name / OB / DL / Herbs / target_name
    → standard_compounds.csv : mol_id, name, ob, dl, herb, source
    → standard_targets.csv   : gene, protein_name, source, confidence（target_name 作 protein_name，gene 留空待 np_standardize）
    → standard_edges.csv     : src_id(mol_id), dst_id(target_name), relation=compound-target, weight=DL, source, evidence

用法
    python tcmsp_adapter.py --input tcmsp_herb_result.csv --outdir np_out/adapted
    python tcmsp_adapter.py --selftest
"""
import argparse
import os
import sys

from base_adapter import BaseAdapter, STANDARD_HEADERS, read_csv, resolve_col

HERE = os.path.dirname(os.path.abspath(__file__))


class TcmspAdapter(BaseAdapter):
    source_name = "TCMSP"

    COMPOUND_COLS = ["mol_id", "mol id", "molecule_id", "MOL-ID", "分子ID"]
    NAME_COLS = ["molecule_name", "name", "成分名", "成分", "化合物", "molecule name"]
    OB_COLS = ["ob", "oral bioavailability", "口服生物利用度", "OB(%)"]
    DL_COLS = ["dl", "drug likeness", "类药性", "DL"]
    HERB_COLS = ["herb", "herbs", "药味", "中药", "药材"]
    TARGET_COLS = ["target_name", "target", "靶点", "Targets", "蛋白名"]

    def _load(self, src_path):
        rows, headers = read_csv(src_path)
        c = {
            "mol_id": resolve_col(headers, self.COMPOUND_COLS),
            "name": resolve_col(headers, self.NAME_COLS),
            "ob": resolve_col(headers, self.OB_COLS),
            "dl": resolve_col(headers, self.DL_COLS),
            "herb": resolve_col(headers, self.HERB_COLS),
            "target": resolve_col(headers, self.TARGET_COLS),
        }
        self.meta["resolved_columns"] = c
        self.log("列映射：%s" % c)
        if not c["name"] and not c["mol_id"]:
            self.meta["warnings"].append("未找到成分名列/ID列，compounds 将为空")
        return rows, c

    def convert_compounds(self, src_path):
        rows, c = self._load(src_path)
        out, seen = [], set()
        for r in rows:
            name = (r.get(c["name"]) or "") if c["name"] else ""
            mid = (r.get(c["mol_id"]) or "") if c["mol_id"] else ""
            key = (mid or name).strip()
            if not key or key in seen:
                if not key:
                    self.meta["dropped_rows"] += 1
                continue
            seen.add(key)
            out.append(self._row({
                "mol_id": mid, "name": name or mid,
                "ob": (r.get(c["ob"]) if c["ob"] else ""),
                "dl": (r.get(c["dl"]) if c["dl"] else ""),
                "herb": (r.get(c["herb"]) if c["herb"] else ""),
                "source": "TCMSP",
            }))
        return out

    def convert_targets(self, src_path):
        rows, c = self._load(src_path)
        if not c["target"]:
            return []
        out, seen = [], set()
        for r in rows:
            tn = (r.get(c["target"]) or "").strip()
            if not tn or tn in seen:
                continue
            seen.add(tn)
            # gene 留空：由 np_standardize.py 统一映射，避免适配器里做不可审计的猜测
            out.append(self._row({
                "gene": "", "protein_name": tn, "source": "TCMSP", "confidence": "database",
            }))
        return out

    def convert_edges(self, src_path):
        rows, c = self._load(src_path)
        if not c["target"]:
            return []
        out, seen = [], set()
        for r in rows:
            mid = (r.get(c["mol_id"]) or r.get(c["name"]) or "").strip() if (c["mol_id"] or c["name"]) else ""
            tn = (r.get(c["target"]) or "").strip()
            if not mid or not tn:
                self.meta["dropped_rows"] += 1
                continue
            key = (mid, tn)
            if key in seen:
                continue
            seen.add(key)
            out.append(self._row({
                "src_id": mid, "dst_id": tn, "relation": "compound-target",
                "weight": (r.get(c["dl"]) if c["dl"] else ""),
                "source": "TCMSP", "evidence": "TCMSP database",
            }))
        return out

    def convert_herbs(self, src_path):
        rows, c = self._load(src_path)
        if not c["herb"]:
            return []
        out, seen = [], set()
        for r in rows:
            hb = (r.get(c["herb"]) or "").strip()
            if not hb or hb in seen:
                continue
            seen.add(hb)
            out.append(self._row({
                "herb_cn": hb, "herb_pinyin": "", "herb_latin": "", "source": "TCMSP",
                "n_comp": "", "n_active": "",
            }))
        return out


def selftest():
    import csv as _csv
    tmp = os.path.join(HERE, "_selftest_tcmsp.csv")
    with open(tmp, "w", encoding="utf-8-sig", newline="") as f:
        w = _csv.writer(f)
        w.writerow(["MOL-ID", "molecule_name", "OB(%)", "DL", "Herbs", "target_name"])
        w.writerow(["MOL001", "quercetin", "46.4", "0.28", "丹参", "AKT1"])
        w.writerow(["MOL002", "luteolin", "36.2", "0.25", "丹参", "TNF"])
        w.writerow(["MOL002", "luteolin", "36.2", "0.25", "丹参", "TNF"])  # 重复行
        w.writerow(["", "", "", "", "", ""])  # 空行
    a = TcmspAdapter(verbose=False)
    meta = a.run(tmp, os.path.join(HERE, "_selftest_out"))
    ok = (meta["produced"]["compounds"]["n_rows"] == 2
          and meta["produced"]["edges"]["n_rows"] == 2
          and meta["dropped_rows"] >= 1)
    print("selftest: %s (compounds=%d edges=%d dropped=%d)" % (
        "PASS" if ok else "FAIL",
        meta["produced"]["compounds"]["n_rows"],
        meta["produced"]["edges"]["n_rows"],
        meta["dropped_rows"]))
    os.remove(tmp)
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="TCMSP → pro 标准表适配器")
    ap.add_argument("--input", help="TCMSP 原始 CSV")
    ap.add_argument("--outdir", default="np_out/adapted")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if not a.input:
        ap.error("需 --input 或 --selftest")
    TcmspAdapter().run(a.input, a.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
