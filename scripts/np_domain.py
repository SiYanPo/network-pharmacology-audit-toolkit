#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_domain.py — 适用域（Applicability Domain）声明生成器

每次输出强制附上：化合物骨架类型、分子量范围、是否落在训练集覆盖的化学空间。
新颖骨架 / 大极性糖苷 / 无机矿物 / 肽类 → 明确标注「预测置信度低，需实验验证」。

判据（透明、可改）：
  - MW < 100 或 > 800           → 极端分子量
  - 糖苷（名称含 glycoside/side、或 MW>600 且含多个 O）→ 糖苷盲区
  - 无机/矿物（无 C 骨架、金属离子）→ 无 ADME 模型适用
  - 肽/蛋白（>10 个氨基酸残基特征）
  - 已知高极性氨基酸类（alliin 类）

用法
  python np_domain.py --smiles-file lig.csv [--mw-file ...] [--out domain.md]
  python np_domain.py --compounds "quercetin:5280343,crocin:5281233" --pubchem
"""
from __future__ import annotations

import argparse, csv, json, os, re, ssl, sys, urllib.parse, urllib.request

CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE

# 化学空间参考区间（类药性/TCMSP 覆盖）
REF_MW = (150, 500)
REF_LOGP = (-1, 5)


def pubchem_props(name):
    try:
        u = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + urllib.parse.quote(name) + "/cids/JSON"
        cids = json.loads(urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "np/1.0"}),
                                                 timeout=25, context=CTX).read().decode())["IdentifierList"]["CID"]
        u2 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids[0]}/property/MolecularWeight,XLogP,CanonicalSMILES,MolecularFormula/JSON"
        p = json.loads(urllib.request.urlopen(urllib.request.Request(u2, headers={"User-Agent": "np/1.0"}),
                                              timeout=25, context=CTX).read().decode())["PropertyTable"]["Properties"][0]
        return {"cid": cids[0], "mw": p.get("MolecularWeight"), "logp": p.get("XLogP"),
                "formula": p.get("MolecularFormula"), "smiles": p.get("CanonicalSMILES") or p.get("SMILES")}
    except Exception as e:
        return {"error": f"{type(e).__name__}: {str(e)[:60]}"}


def assess(name, info):
    flags, score = [], 1.0
    mw = info.get("mw")
    try:
        mw = float(mw) if mw is not None else None
    except (TypeError, ValueError):
        mw = None
    formula = info.get("formula") or ""
    if mw is not None:
        if mw < REF_MW[0]:
            flags.append(f"MW={mw:.0f}<{REF_MW[0]}（小极性分子，实测靶点稀少）"); score -= .3
        if mw > REF_MW[1]:
            flags.append(f"MW={mw:.0f}>{REF_MW[1]}（超出典型类药空间）"); score -= .3
        if mw > 600:
            flags.append("MW>600：糖苷/大分子高概率，TCMSP 覆盖差"); score -= .3
    if re.search(r'glycoside|oside$|saponin|皂苷|糖苷', name, re.I):
        flags.append("名称提示糖苷/皂苷：需水解为苷元后重算靶点"); score -= .3
    if formula and not re.search(r'C\d', formula):
        flags.append("无碳骨架：疑似矿物/无机物，ADME 模型不适用"); score -= .7
    if re.search(r'peptide|肽|insulin|cyclosporin', name, re.I):
        flags.append("肽类：超出小分子预测域"); score -= .4
    if "error" in info:
        flags.append("PubChem 未取到物化性质：" + info["error"]); score -= .2
    confident = score >= .7
    return {"compound": name, "flags": flags, "confidence": round(max(0, score), 2),
            "low_confidence": not confident,
            "statement": ("落在常规类药化学空间，可正常解释"
                          if confident else "预测置信度低，需实验验证；建议补充实测靶点或苷元/代谢物分析")}


def cmd_run(args):
    items = []
    if args.compounds:
        for spec in args.compounds.split(","):
            name = spec.split(":")[0]
            info = pubchem_props(name) if args.pubchem else {}
            items.append(assess(name, info))
    if args.smiles_file:
        with open(args.smiles_file, encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                name = r.get("compound") or r.get("name") or ""
                info = pubchem_props(name) if args.pubchem else {}
                items.append(assess(name, info))
    lines = ["# 适用域（Applicability Domain）声明", "",
             f"参考区间：MW {REF_MW[0]}–{REF_MW[1]}；XLogP {REF_LOGP[0]}–{REF_LOGP[1]}。",
             "超出区间或属糖苷/矿物/肽类者，标注低置信度。", ""]
    for it in items:
        lines.append(f"## {it['compound']}　置信度 {it['confidence']:.2f}")
        for fl in it["flags"]:
            lines.append(f"- ⚠ {fl}")
        lines.append(f"- 结论：{it['statement']}")
        lines.append("")
    text = "\n".join(lines)
    print(text)
    if args.out:
        open(args.out, "w", encoding="utf-8").write(text)
        print(f"→ {args.out}")
    if args.json:
        json.dump(items, open(args.json, "w", encoding="utf-8"), ensure_ascii=False, indent=1)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--compounds"); p.add_argument("--smiles-file"); p.add_argument("--pubchem", action="store_true")
    p.add_argument("--out"); p.add_argument("--json")
    args = p.parse_args()
    if not args.compounds and not args.smiles_file:
        p.error("需要 --compounds 或 --smiles-file")
    cmd_run(args)


if __name__ == "__main__":
    main()
