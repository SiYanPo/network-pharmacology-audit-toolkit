#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
base_adapter.py —— 统一数据契约适配器基类（解决跨 skill 数据契约不统一）

问题（审计发现）
    市场各 skill 的输出列名五花八门：ctd-network 用 compound/disease，network-pharmacology
    用 gene/symbol，TCMSP 用 mol_id/成分名，pro 的 contracts/standard_tables.md 定义五张标准表。
    "缝合怪"跨 skill 调用的最大壁垒就是契约不统一 → 本基类定义统一转换接口。

设计（契约驱动 + 审计内建）
    * 每个适配器把某来源的原始表 → 标准表（standard_herbs/compounds/targets/edges）。
    * 列名映射用"候选列名表"而非硬编码（L1 精确 → L2 归一化 → L3 关键词模糊）。
    * 所有转换记录"丢弃行/缺失字段/命中列"到 meta，写 degraded_fields（延续 fallback 原则：禁止静默丢数据）。
    * 改标准表列名 = 破坏性变更（bump major），见 contracts/standard_tables.md。

标准表契约（来自 contracts/standard_tables.md）
    standard_herbs.csv     : herb_cn, herb_pinyin, herb_latin, source, n_comp, n_active
    standard_compounds.csv : mol_id, name, ob, dl, herb, source
    standard_targets.csv   : gene, protein_name, source, confidence
    standard_edges.csv     : src_id, dst_id, relation, weight, source, evidence
"""
import csv
import json
import os
import re
from abc import ABC, abstractmethod

SCHEMA_VERSION = "1.0"

STANDARD_HEADERS = {
    "herbs": ["herb_cn", "herb_pinyin", "herb_latin", "source", "n_comp", "n_active"],
    "compounds": ["mol_id", "name", "ob", "dl", "herb", "source"],
    "targets": ["gene", "protein_name", "source", "confidence"],
    "edges": ["src_id", "dst_id", "relation", "weight", "source", "evidence"],
}


def norm(s):
    return re.sub(r"[\s\-_/,.、，;；:：]+", "", str(s or "").strip().lower())


def resolve_col(headers, candidates):
    """L1 精确 → L2 归一化 → L3 关键词包含。返回真实列名或 None。"""
    hmap = {norm(h): h for h in headers}
    for c in candidates:
        if c in headers:
            return c
    for c in candidates:
        if norm(c) in hmap:
            return hmap[norm(c)]
    for c in candidates:
        nc = norm(c)
        for k, orig in hmap.items():
            if nc and (nc in k or k in nc):
                return orig
    return None


def read_csv(path):
    with open(path, encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f)
        return list(rdr), (rdr.fieldnames or [])


def write_csv(path, headers, rows):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


class BaseAdapter(ABC):
    """所有来源适配器的基类。子类实现 source_name 与 convert* 方法。"""

    source_name = "base"

    def __init__(self, verbose=True):
        self.verbose = verbose
        self.meta = {
            "adapter": self.__class__.__name__,
            "source": self.source_name,
            "schema_version": SCHEMA_VERSION,
            "dropped_rows": 0,
            "missing_fields": {},
            "resolved_columns": {},
            "warnings": [],
        }

    def log(self, msg):
        if self.verbose:
            print("[%s] %s" % (self.source_name, msg))

    def _row(self, mapping):
        """按标准表头补齐缺失字段为空串，并记录 degraded_fields。"""
        out = {}
        for k, v in mapping.items():
            out[k] = "" if v is None else str(v).strip()
            if v in (None, "", "NA", "nan", "None"):
                self.meta["missing_fields"][k] = self.meta["missing_fields"].get(k, 0) + 1
        return out

    @abstractmethod
    def convert_compounds(self, src_path):
        """→ list[dict] 符合 standard_compounds.csv。"""

    def convert_targets(self, src_path):
        return []

    def convert_edges(self, src_path):
        return []

    def convert_herbs(self, src_path):
        return []

    def run(self, src_path, outdir, tables=("compounds", "targets", "edges", "herbs")):
        """执行转换并写出标准表 + meta.json。返回 meta。"""
        os.makedirs(outdir, exist_ok=True)
        converters = {
            "compounds": self.convert_compounds,
            "targets": self.convert_targets,
            "edges": self.convert_edges,
            "herbs": self.convert_herbs,
        }
        produced = {}
        for t in tables:
            rows = converters[t](src_path) or []
            headers = STANDARD_HEADERS[t]
            outp = os.path.join(outdir, "standard_%s.csv" % t)
            write_csv(outp, headers, rows)
            produced[t] = {"file": os.path.basename(outp), "n_rows": len(rows)}
            self.log("%s → %d 行" % (t, len(rows)))
        self.meta["produced"] = produced
        meta_path = os.path.join(outdir, "adapter_meta.json")
        with open(meta_path, "w", encoding="utf-8") as f:
            json.dump(self.meta, f, ensure_ascii=False, indent=1)
        self.log("meta → %s" % meta_path)
        if self.meta["dropped_rows"]:
            self.meta["warnings"].append("丢弃 %d 行（字段不足）" % self.meta["dropped_rows"])
        return self.meta


CLI_HELP = "适配器把外部来源原始表转换为 pro 标准表（contracts/standard_tables.md）"
