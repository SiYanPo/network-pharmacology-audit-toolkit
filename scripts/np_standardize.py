#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_standardize.py — 靶点名称标准化与可审计映射（第 1 层质控）

为什么单独成阶段
    文献与实测数据表明，"靶点未标准化/同物异名未合并"是网络药理学最常见的致命错误之一：
    同一个蛋白在 TCMSP 里叫 Prothrombin、在富集里要 F2，若不统一，交集与网络会被系统性稀释。
    本阶段把"映射"变成可审计动作：每条记录都留下 命中层级 / 判定理由 / 是否需人工复核。

四层匹配（与业界做法一致，但补充了"非人类噪音"显式判定）
    L1 exact      与内置蛋白名词典规范化后完全一致
    L2 normalized 去标点/括号/希腊字母/大小写/亚型后缀后一致
    L3 structural 拆解 "推荐名 (别名) [Cleaved into: ...] Includes: ..." 后匹配
    L4 online     --online 时用 UniProt REST 校验并回填（需网络）
    -- 未命中     进入 unmatched.csv，并给出疑似类别（非人类 / mRNA 噪音 / 需人工）

用法
    python np_standardize.py -i targets.csv -o out/
    python np_standardize.py -i genes.txt -o out/            # 纯基因列表
    python np_standardize.py -i targets.csv -o out/ --online  # 联网校验
"""
import argparse
import os
import re
import sys
import urllib.parse
import urllib.request

import npcore as C

HERE = os.path.dirname(os.path.abspath(__file__))
# 词典优先级：UniProt 全量（20k 蛋白、11 万条名称）> 手工精简表（325 行）
DICT_FULL = os.path.join(HERE, "..", "data", "protein_symbols_uniprot.tsv.gz")
DICT_SMALL = os.path.join(HERE, "..", "data", "protein_symbols.tsv")
DICT_PATH = DICT_FULL if os.path.exists(DICT_FULL) else DICT_SMALL

GREEK = {"alpha": "a", "beta": "b", "gamma": "g", "delta": "d", "epsilon": "e",
         "zeta": "z", "eta": "e", "theta": "t", "iota": "i", "kappa": "k",
         "lambda": "l", "mu": "m", "nu": "n", "xi": "x", "omicron": "o",
         "pi": "p", "rho": "r", "sigma": "s", "tau": "t", "upsilon": "u",
         "phi": "f", "chi": "c", "psi": "p", "omega": "w",
         "α": "a", "β": "b", "γ": "g", "δ": "d", "ε": "e", "κ": "k", "μ": "m", "π": "p", "σ": "s", "ω": "w"}

NONHUMAN_PAT = re.compile(
    r"(escherichia|staphylococc|pseudomonas|bacill|mycobacter|plasmodi|trypanosom|leishman|"
    r"candida|saccharomyces|yeast|drosophila|caenorhabditis|danio|zebrafish|xenopus|"
    r"bos taurus|equus|sus scrofa|rattus|mouse|murine|mus musculus|chlorocebus|"
    r"rhinovirus|influenza|hepatitis c virus|hiv|sars-cov|virus|cam\b|bacterial|plant|"
    r"arabidopsis|oryza|nicotiana)", re.I)

NOISE_PAT = re.compile(
    r"(mRNA of|mRNA|nucleotide|expressed sequence|clone|fragment|unnamed|unknown|"
    r"hypothetical|putative uncharacterized|^$)", re.I)

SYMBOL_OK = re.compile(r"^[A-Z][A-Z0-9\-]{1,14}$")


def norm(s):
    s = (s or "").strip().lower()
    for k, v in GREEK.items():
        s = s.replace(k, v)
    s = re.sub(r"\(.*?\)|\[.*?\]", " ", s)
    s = re.sub(r"[\-_/,.、，;:;'\"’‘“”]+", " ", s)
    s = re.sub(r"\s+", " ", s)
    return s.strip()


def strip_isoform(s):
    """只剥离真正的亚型后缀（isoform 2 / variant 3），不动 subunit（它是有意义的名字成分）。"""
    return re.sub(r"\s*\(?\b(isoform|variant)\s*[0-9]{1,3}[A-Za-z]?\)?$", "", (s or "").strip(), flags=re.I)


def load_dict(path=None):
    path = path or DICT_PATH
    if not os.path.exists(path):
        return {}, {}
    _h, rows = C.read_table(path)
    exact, normd = {}, {}
    for r in rows:
        pn = (r.get("protein_name") or "").strip()
        sym = (r.get("gene_symbol") or "").strip()
        if not pn:
            continue
        if sym in ("", "—", "-"):
            exact[pn.lower()] = ("", r.get("note") or "已知非人类/噪音条目")
            normd[norm(pn)] = ("", r.get("note") or "已知非人类/噪音条目")
            continue
        exact[pn.lower()] = (sym, r.get("note") or "")
        normd[norm(pn)] = (sym, r.get("note") or "")
    return exact, normd


def structural_candidates(raw):
    """把 "A (B, C) [Cleaved into: D] Includes: E" 拆成候选名。"""
    s = raw
    cands = []
    m = re.search(r"Includes:\s*(.+)$", s, re.I)
    if m:
        cands.append(m.group(1))
    m = re.search(r"Cleaved into:\s*(.+?)(?:\]|$)", s, re.I)
    if m:
        cands.extend(re.split(r"[;,]", m.group(1)))
    for grp in re.findall(r"\(([^()]*)\)", s):
        cands.extend(re.split(r"[;,]", grp))
    cands.append(re.sub(r"\([^()]*\)|\[[^\]]*\]", " ", s))
    for part in re.split(r"[,;]", re.sub(r"\([^()]*\)|\[[^\]]*\]", " ", s)):
        cands.append(part)
    out = []
    for c in cands:
        c = re.sub(r"\s+", " ", c).strip(" .;,")
        if c and len(c) > 2:
            out.append(c)
    return out


def uniprot_lookup(symbol, timeout=25):
    url = ("https://rest.uniprot.org/uniprotkb/search?query="
           + urllib.parse.quote('gene_exact:%s AND organism_id:9606 AND reviewed:true' % symbol)
           + "&fields=accession,gene_primary&format=tsv&size=1")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "np-pro/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            txt = r.read().decode("utf-8", "replace")
        lines = [l for l in txt.splitlines() if l.strip()]
        if len(lines) >= 2:
            parts = lines[1].split("\t")
            if len(parts) >= 2 and parts[1]:
                return parts[0], parts[1]
    except Exception:
        return None
    return None


def main(argv=None):
    ap = C.Cli(description="靶点名称标准化：蛋白名→HGNC 基因符号，可审计、可对照")
    ap.add_argument("-i", "--input", required=True, help="靶点表（含 target_name/gene 列）或纯基因列表")
    ap.add_argument("-o", "--outdir", default="np_out", help="输出目录")
    ap.add_argument("--col", default=None, help="指定列名（默认自动识别 target_name/gene）")
    ap.add_argument("--online", action="store_true", help="联网用 UniProt 校验并回填符号")
    ap.add_argument("--keep-nonhuman", action="store_true", help="保留非人类/噪音条目（默认剔除并单列）")
    ap.add_argument("--alias", default=None, help="追加人工别名表 CSV：raw,symbol")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    man = C.Manifest("01_standardize", a.outdir, vars(a))
    man.input("raw", a.input)

    # ---- 读入 ----
    # 输入若带 organism/species/taxid 列，必须用起来：
    # 鼠/牛/猪源活性记录若只靠"名字里有没有物种字样"来判定，会 100% 漏掉
    # （ChEMBL 的 pref_name 通常是 "Intestinal-type alkaline phosphatase" 这种不带物种的名字）。
    c_org = None
    raw_list = []
    if a.input.lower().endswith((".txt", ".list")) or not a.col:
        try:
            _h, rows = C.read_table(a.input)
            col = a.col or C.resolve_col(_h, "target_name") or C.resolve_col(_h, "gene")
        except SystemExit:
            rows, col = [], None
        if rows and col:
            c_org = C.resolve_col(_h, "organism", extra=["species", "taxid", "organism_name"])
            raw_list = [(r.get(col, ""), r) for r in rows]
        else:
            raw_list = [(g, {}) for g in C.read_gene_list(a.input)]
    else:
        _h, rows = C.read_table(a.input)
        col = a.col
        c_org = C.resolve_col(_h, "organism", extra=["species", "taxid", "organism_name"])
        raw_list = [(r.get(col, ""), r) for r in rows]
    if c_org:
        man.count("organism_column", c_org)

    exact, normd = load_dict()
    alias = {}
    if a.alias and os.path.exists(a.alias):
        _h, arows = C.read_table(a.alias)
        for r in arows:
            k = list(r.values())[0]
            v = list(r.values())[1]
            if k and v:
                alias[norm(k)] = v.strip().upper()
        man.count("alias_rows", len(alias))

    results, unmatched = [], []
    layer_all = {}          # 全部记录的层级计数（含被剔除的非人类/噪音）
    for raw, _row in raw_list:
        raw = (raw or "").strip()
        if not raw or raw.lower() in ("na", "nan", "none", "-", "—"):
            continue
        base = strip_isoform(raw)
        sym, layer, note = "", "", ""

        # 显式物种列优先：名字里往往没有物种信息，只靠正则必然漏
        forced_nonhuman = False
        if c_org:
            org = str(_row.get(c_org) or "").strip()
            if org and not re.search(r"homo sapiens|human|9606", org, re.I):
                forced_nonhuman = True
                layer = "nonhuman"
                note = "非人类来源（organism=%s，按物种列剔除）" % org[:40]

        if forced_nonhuman:
            pass
        elif SYMBOL_OK.match(base) and not re.search(r"\s", base):
            # 看起来已经是基因符号；仍过一遍词典（有些符号就是蛋白名，如 TH）
            sym, layer, note = base, "as_symbol", "输入本身为符号形式"

        if not sym and not forced_nonhuman:
            key = base.lower()
            if key in exact:
                sym, note = exact[key]
                layer = "L1_exact" if sym else "nonhuman"
            elif norm(base) in exact:
                sym, note = exact[norm(base)]
                layer = "L2_normalized" if sym else "nonhuman"
            elif norm(base) in normd:
                sym, note = normd[norm(base)]
                layer = "L2_normalized" if sym else "nonhuman"
            elif norm(base) in alias:
                sym, layer, note = alias[norm(base)], "L4_alias", "人工别名表"
            else:
                for cand in structural_candidates(base):
                    for probe in (cand.lower(), norm(cand)):
                        if probe in exact and exact[probe][0]:
                            sym = exact[probe][0]
                            layer, note = "L3_structural", "由结构拆解命中：%s" % cand
                            break
                        if probe in normd and normd[probe][0]:
                            sym = normd[probe][0]
                            layer, note = "L3_structural", "由结构拆解命中：%s" % cand
                            break
                    if sym:
                        break

        if not sym and not forced_nonhuman:
            if NONHUMAN_PAT.search(base):
                layer, note = "nonhuman", "疑似非人类来源（菌/病毒/动物/植物）"
            elif NOISE_PAT.search(base):
                layer, note = "noise", "疑似数据库噪音条目（mRNA/clone/unknown）"
            else:
                layer, note = "unmatched", "未在内置词典命中，需人工确认或联网映射"

        sym = (sym or "").strip().upper()
        rec = dict(_row) if _row else {}
        rec.update({"raw": raw, "symbol": sym, "layer": layer, "note": note})
        layer_all[layer] = layer_all.get(layer, 0) + 1
        if sym:
            if a.online:
                hit = uniprot_lookup(sym)
                if hit:
                    rec["uniprot"] = hit[0]
                    rec["uniprot_gene"] = hit[1]
                    if hit[1].upper() != sym:
                        rec["note"] = (note + "；UniProt 主基因名=%s" % hit[1]).strip("；")
                        rec["flag"] = "symbol_mismatch"
                else:
                    rec["flag"] = "uniprot_not_found"
            results.append(rec)
        else:
            unmatched.append(rec)
            if a.keep_nonhuman:
                results.append(rec)

    # ---- 输出 ----
    C.ensure_dir(a.outdir)
    extra = [h for h in (_h if "_h" in dir() else []) if h not in ("raw", "symbol")]
    heads = ["raw"] + extra + ["symbol", "layer", "note", "uniprot", "uniprot_gene", "flag"]
    heads = [h for h in heads if h not in ("raw",)]
    heads = ["raw"] + [h for h in heads if h != "raw"]
    p1 = C.write_table(os.path.join(a.outdir, "targets_standardized.csv"), heads, results)
    seen_u, uniq = set(), []
    for r in unmatched:
        if r["raw"] not in seen_u:
            seen_u.add(r["raw"])
            uniq.append(r)
    p2 = C.write_table(os.path.join(a.outdir, "targets_unmatched.csv"), ["raw", "layer", "note"], uniq)

    genes = sorted({r["symbol"] for r in results if r["symbol"]})
    p3 = C.write_text(os.path.join(a.outdir, "targets_gene_list.txt"), "\n".join(genes) + "\n")

    layers = {}
    for r in results + unmatched:
        layers[r["layer"]] = layers.get(r["layer"], 0) + 1
    total = len(results) + len(unmatched)
    mapped = len(genes)
    rate = mapped / max(1, len({r["raw"] for r in results + unmatched}))

    md = ["# 靶点标准化审计", "",
          "- 输入：`%s`" % a.input,
          "- 记录数：%d，去重后基因：%d，映射率（去重基因/去重原始名）：%.1f%%" % (total, mapped, 100 * rate),
          "", "## 命中层级分布", "", "| 层级 | 条数 | 含义 |", "|---|---|---|"]
    MEAN = {"as_symbol": "输入即符号", "L1_exact": "蛋白名精确命中", "L2_normalized": "规范化后命中",
            "L3_structural": "结构拆解命中（建议抽检）", "L4_alias": "人工别名表",
            "nonhuman": "非人类来源，应剔除", "noise": "数据库噪音", "unmatched": "未命中"}
    for k, v in sorted(layers.items(), key=lambda x: -x[1]):
        md.append("| %s | %d | %s |" % (k, v, MEAN.get(k, "")))
    if unmatched:
        md += ["", "## 未命中条目（前 40）", ""]
        for r in unmatched[:40]:
            md.append("- `%s` → %s（%s）" % (r["raw"], r["layer"], r["note"]))
        md += ["", "> 处理建议：非人类/噪音条目直接剔除；疑似真人类靶点补入 `--alias` 别名表后重跑。",
               "> **不得**为了凑数把非人类蛋白强行映射到人类基因——这会造成系统性假阳性。"]
    if any(r.get("flag") for r in results):
        n = sum(1 for r in results if r.get("flag"))
        md += ["", "## 需人工复核：%d 条" % n,
               "UniProt 主基因名与输入符号不一致，或符号在 UniProt 人类 reviewed 库中查不到。"]
    p4 = C.write_text(os.path.join(a.outdir, "standardize_audit.md"), "\n".join(md) + "\n")

    for k, v in [("std_table", p1), ("unmatched", p2), ("gene_list", p3), ("audit", p4)]:
        man.output(k, v)
    man.count("records", total)
    man.count("genes", mapped)
    man.count("map_rate", round(rate, 4))
    man.count("layers", layer_all or layers)
    if rate < 0.8:
        man.warn("映射率低于 80%%，交集与网络可能被系统性稀释，建议补别名表：%s" % p2)
    if unmatched:
        man.warn("存在 %d 条未命中，已在 %s 列出" % (len(unmatched), p2))
    man.save()

    if not a.quiet:
        print("[标准化] 记录 %d → 基因 %d（映射率 %.1f%%）" % (total, mapped, 100 * rate))
        src = layer_all or layers
        print("[标准化] 层级：" + "，".join("%s=%d" % (k, v) for k, v in sorted(src.items(), key=lambda x: -x[1])))
        if unmatched:
            print("[标准化] 未命中/剔除 %d 条（含非人类 %d）→ %s"
              % (len(unmatched), layer_all.get("nonhuman", 0), p2))
        print("[标准化] 输出 → %s" % p1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
