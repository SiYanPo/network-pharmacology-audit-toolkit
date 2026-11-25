#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_alias.py —— 四层别名/同义词映射与双向索引校验（补齐"生僻药/矿物药/错别字"短板）

与 np_standardize.py 的分工（避免重复造轮子）
    np_standardize.py : 蛋白名 → HGNC 基因符号（L1 exact…L4 online，已有）
    np_alias.py       : ① 中药名/成分名 同义词消歧（龟甲=龟板=败龟板；制首乌=何首乌）
                        ② 靶点 accession↔gene 双向索引校验（复刻市场 G skill 的 validate_aliases）
                        两者输入不同，互补。

四层匹配（herb 模式）
    L1 exact       内置同义词表精确命中（规范名→标准名）
    L2 normalized  去炮制前缀（生/制/酒/醋/盐/姜/蜜/焦/炭）、去规格、繁简、大小写、空白
    L3 structural  拆解 "拉丁名 (别名1, 别名2) [炮制]" 取候选逐个匹配
    L4 external    --alias-file 追加人工 CSV（raw,symbol），未命中追加到 alias_pending.csv

四档判定（validate 模式，靶点双向索引）
    PASS  accession→gene 与 gene→accession 双向命中且不冲突
    WARN  仅单向命中 / 两者皆不在人类 reviewed 表
    FAIL  两者都在表但指向不同（真实错配）
    PENDING 空值
    诚实声明：PASS 只证明"自洽"，不证明"语义正确"（把 Thrombin 错配到真实蛋白仍 PASS）。

用法
    python np_alias.py herbs  --input herbs.txt -o out/ [--alias-file manual.csv]
    python np_alias.py validate --input targets.csv --uniprot data/protein_symbols.tsv -o out/
    python np_alias.py selftest
"""
import argparse
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

# ---- 内置中药同义词表（规范名 → 标准名）—— 覆盖审计中暴露的生僻/矿物/动物/炮制痛点 ----
HERB_SYNONYMS = {
    # 动物药 / 矿物药（TCMSP 常缺，最需联想）
    "龟甲": "龟甲", "龟板": "龟甲", "败龟板": "龟甲", "生龟板": "龟甲", "炙龟板": "龟甲",
    "龟壳": "龟甲", "乌龟壳": "龟甲", "carapax testudinis": "龟甲",
    "龙骨": "龙骨", "煅龙骨": "龙骨", "生龙骨": "龙骨", "fossilis ossis": "龙骨",
    "牡蛎": "牡蛎", "煅牡蛎": "牡蛎", "生牡蛎": "牡蛎", "蛎蛤": "牡蛎",
    "磁石": "磁石", "煅磁石": "磁石", "灵磁石": "磁石",
    "朱砂": "朱砂", "辰砂": "朱砂", "丹砂": "朱砂",
    "石膏": "石膏", "生石膏": "石膏", "煅石膏": "石膏",
    "滑石": "滑石", "滑石块": "滑石",
    # 果实/种子类（楮实子等）
    "楮实子": "楮实子", "楮实": "楮实子", "楮桃": "楮实子", "谷树子": "楮实子",
    "构树子": "楮实子", "fructus broussonetiae": "楮实子",
    "枸杞子": "枸杞子", "枸杞": "枸杞子", "甘枸杞": "枸杞子", "西枸杞": "枸杞子",
    "益智仁": "益智仁", "益智": "益智仁", "益智子": "益智仁",
    "山茱萸": "山茱萸", "山萸肉": "山茱萸", "萸肉": "山茱萸", "枣皮": "山茱萸",
    # 炮制系列（何首乌/地黄等）
    "何首乌": "何首乌", "首乌": "何首乌", "生首乌": "何首乌", "制何首乌": "何首乌",
    "制首乌": "何首乌", "夜交藤": "首乌藤", "首乌藤": "首乌藤",
    "熟地黄": "熟地黄", "熟地": "熟地黄", "制地黄": "熟地黄", "酒熟地": "熟地黄",
    "生地黄": "生地黄", "生地": "生地黄", "干地黄": "生地黄",
    # 常用（错别字/别名）
    "丹参": "丹参", "紫丹参": "丹参", "赤丹参": "丹参", "salvia miltiorrhiza": "丹参",
    "川芎": "川芎", "芎䓖": "川芎", "西川芎": "川芎",
    "石菖蒲": "石菖蒲", "菖蒲": "石菖蒲", "九节菖蒲": "石菖蒲",
    "郁金": "郁金", "广郁金": "郁金", "川郁金": "郁金",
    "赤芍": "赤芍", "赤芍药": "赤芍", "山芍药": "赤芍",
    "黄精": "黄精", "制黄精": "黄精", "酒黄精": "黄精", "鸡头黄精": "黄精",
    "桃仁": "桃仁", "山桃仁": "桃仁", "红花": "红花", "草红花": "红花",
    "益智": "益智仁",
}
# 炮制前缀（L2 归一化剥离）；注意 "生/制" 单独也作药名成分，故仅在"前缀+已知标准名"时剥离
PROCESS_PREFIX = ("煅", "生", "制", "酒", "醋", "盐", "姜", "蜜", "焦", "炭", "炒", "炙", "蒸")
TRAD_SIMPLE = {"龜": "龟", "楮": "楮", "實": "实", "黃": "黄", "精": "精", "茱": "茱",
               "萸": "萸", "熟": "熟", "地": "地", "首": "首", "烏": "乌", "靈": "灵",
               "磁": "磁", "煅": "煅", "龍": "龙", "骨": "骨", "牡": "牡", "蠣": "蛎"}


def norm(s):
    """L2 归一化：繁简、大小写、空白、标点；不剥离炮制前缀（那在 match 中按候选做）。"""
    s = (s or "").strip().lower()
    for k, v in TRAD_SIMPLE.items():
        s = s.replace(k, v)
    s = re.sub(r"[\s\-_/,.、，;；:'\"'“”（）()\[\]]+", "", s)
    return s


def structural_candidates(raw):
    """L3 结构拆解：'X (Y, Z) [W]' → [整体, X, Y, Z, W]。"""
    s = (raw or "").strip()
    cands = [s]
    for grp in re.findall(r"[（(]([^()（）]*)[）)]", s):
        cands.extend([x.strip() for x in re.split(r"[,，;；/]", grp)])
    for grp in re.findall(r"[\[【]([^\[\]【】]*)[】\]]", s):
        cands.append(grp.strip())
    cands.append(re.sub(r"[（(][^()（）]*[）)]|\[[^\]]*\]|【[^】]*】", "", s).strip())
    out, seen = [], set()
    for c in cands:
        c = c.strip()
        if c and c not in seen:
            seen.add(c)
            out.append(c)
    return out


def load_alias_file(path):
    d = {}
    if not path or not os.path.exists(path):
        return d
    with open(path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            vals = list(row.values())
            if len(vals) >= 2 and vals[0] and vals[1]:
                d[norm(vals[0])] = vals[1].strip()
    return d


def match_herb(raw, alias_file=None):
    """返回 (standard_name or '', layer, note)。"""
    raw = (raw or "").strip()
    if not raw:
        return "", "empty", "空输入"
    low = raw.lower()
    # L1 精确
    if low in HERB_SYNONYMS:
        return HERB_SYNONYMS[low], "L1_exact", "内置同义词表精确命中"
    # L2 归一化（含炮制前缀剥离候选）
    n = norm(raw)
    if n in HERB_SYNONYMS:
        return HERB_SYNONYMS[n], "L2_normalized", "归一化命中"
    for p in PROCESS_PREFIX:
        if raw.startswith(p) and len(raw) > len(p) + 1:
            stripped = norm(raw[len(p):])
            if stripped in HERB_SYNONYMS:
                return HERB_SYNONYMS[stripped], "L2_normalized", "剥离炮制前缀「%s」后命中" % p
    # L3 结构拆解
    for c in structural_candidates(raw):
        for probe in (c.lower(), norm(c)):
            if probe in HERB_SYNONYMS:
                return HERB_SYNONYMS[probe], "L3_structural", "结构拆解候选「%s」命中" % c
            for p in PROCESS_PREFIX:
                if c.startswith(p) and norm(c[len(p):]) in HERB_SYNONYMS:
                    return HERB_SYNONYMS[norm(c[len(p):])], "L3_structural", "候选剥离炮制前缀后命中"
    # L4 外部别名表
    ext = alias_file or {}
    if n in ext:
        return ext[n], "L4_external", "人工别名表命中"
    for c in structural_candidates(raw):
        if norm(c) in ext:
            return ext[norm(c)], "L4_external", "人工别名表（结构候选）命中"
    return "", "unmatched", "四层均未命中，需追加人工别名表"


def load_uniprot_index(path):
    """读 protein_symbols.tsv 建 acc↔gene 双向索引。"""
    acc2gene, gene2accs = {}, {}
    if not path or not os.path.exists(path):
        return acc2gene, gene2accs
    with open(path, encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for r in rdr:
            acc = (r.get("accession") or "").strip().upper()
            g = (r.get("gene_symbol") or "").strip().upper()
            if not g or g in ("—", "-"):
                continue
            if acc:
                acc2gene[acc] = g
            gene2accs.setdefault(g, set()).add(acc)
    return acc2gene, gene2accs


def validate_pair(accession, gene, acc2gene, gene2accs):
    """四档判定：PASS / WARN / FAIL / PENDING。"""
    acc = (accession or "").strip().upper()
    gl = (gene or "").strip().upper()
    if not acc and not gl:
        return "PENDING", "两者皆空"
    acc_hit = acc in acc2gene
    gene_hit = gl in gene2accs
    if acc_hit and gene_hit:
        if acc2gene[acc] == gl:
            return "PASS", "双向命中且一致"
        return "FAIL", "错配：%s 应为 %s，实际给 %s" % (acc, acc2gene[acc], gl)
    if acc_hit ^ gene_hit:
        return "WARN", "仅单向命中（%s）" % ("accession" if acc_hit else "gene")
    return "WARN", "两者均不在人类 reviewed 表；可能是新符号/非人类/噪音，需人工确认"


def cmd_herbs(a):
    ext = load_alias_file(a.alias_file)
    raws = [l.strip() for l in open(a.input, encoding="utf-8-sig") if l.strip()] \
        if os.path.isfile(a.input) else [x.strip() for x in a.input.split(",") if x.strip()]
    os.makedirs(a.outdir, exist_ok=True)
    rows, pending = [], []
    for raw in raws:
        std, layer, note = match_herb(raw, ext)
        rows.append({"raw": raw, "standard": std or "（需人工）", "layer": layer, "note": note})
        if not std:
            pending.append({"raw": raw, "standard": "", "note": "请填写规范名"})
    out = os.path.join(a.outdir, "herb_alias_result.csv")
    with open(out, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["raw", "standard", "layer", "note"])
        w.writeheader(); w.writerows(rows)
    pend_path = os.path.join(a.outdir, "alias_pending.csv")
    with open(pend_path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["raw", "standard", "note"])
        w.writeheader(); w.writerows(pending)
    hit = sum(1 for r in rows if r["standard"] != "（需人工）")
    layers = {}
    for r in rows:
        layers[r["layer"]] = layers.get(r["layer"], 0) + 1
    print("[herbs] %d 味，命中 %d（%.1f%%），层级 %s" %
          (len(rows), hit, 100.0 * hit / max(1, len(rows)), layers))
    print("[herbs] → %s ; 待补 %s" % (out, pend_path))
    return 0


def cmd_validate(a):
    acc2gene, gene2accs = load_uniprot_index(a.uniprot)
    if not acc2gene:
        print("[WARN] UniProt 索引为空，校验退化为「全 WARN」；请提供 protein_symbols.tsv", file=sys.stderr)
    rows = []
    with open(a.input, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            acc = r.get("accession") or r.get("uniprot") or ""
            gene = r.get("gene") or r.get("gene_symbol") or r.get("symbol") or ""
            status, note = validate_pair(acc, gene, acc2gene, gene2accs)
            rows.append({"accession": acc, "gene": gene, "status": status, "note": note})
    os.makedirs(a.outdir, exist_ok=True)
    out = os.path.join(a.outdir, "target_alias_validation.csv")
    with open(out, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["accession", "gene", "status", "note"])
        w.writeheader(); w.writerows(rows)
    cnt = {}
    for r in rows:
        cnt[r["status"]] = cnt.get(r["status"], 0) + 1
    print("[validate] %d 行：%s" % (len(rows), cnt))
    print("[validate] → %s" % out)
    if cnt.get("FAIL", 0):
        print("[FAIL] 存在 %d 条真实错配，必须人工修正后方可进入交集。" % cnt["FAIL"], file=sys.stderr)
        return 1
    return 0


def cmd_selftest(a):
    cases = [
        ("龟板", "龟甲", "L2_normalized"),
        ("败龟板", "龟甲", "L2_normalized"),
        ("制首乌", "何首乌", "L2_normalized"),
        ("楮实", "楮实子", "L2_normalized"),
        ("Fructus Broussonetiae", "楮实子", "L1_exact"),
        ("煅龙骨", "龙骨", "L2_normalized"),
        ("X (Y, Z) [W]", "", "unmatched"),
    ]
    ok = 0
    for raw, expect, exp_layer in cases:
        got, layer, _ = match_herb(raw)
        good = (got == expect) if expect else (got == "")
        ok += 1 if good else 0
        print("  %-26s → %-8s [%s] %s" % (raw, got or "-", layer, "OK" if good else "FAIL exp=" + expect))
    # 双向索引自测
    a2g, g2a = {"P12345": "TP53"}, {"TP53": {"P12345"}, "MDM2": {"P00918"}}
    for acc, gene, exp in [("P12345", "TP53", "PASS"),
                           ("P12345", "MDM2", "FAIL"),
                           ("P99999", "MDM2", "WARN"),
                           ("", "", "PENDING")]:
        st, _ = validate_pair(acc, gene, a2g, g2a)
        good = st == exp
        ok += 1 if good else 0
        print("  pair %-8s/%-6s → %-7s %s" % (acc or "-", gene or "-", st, "OK" if good else "FAIL exp=" + exp))
    print("selftest: %d/%d PASS" % (ok, len(cases) + 4))
    return 0 if ok == len(cases) + 4 else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_alias: 四层别名映射 + 双向索引校验")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("herbs", help="中药名四层同义词映射")
    p.add_argument("--input", required=True, help="药名文件或逗号分隔药名")
    p.add_argument("--alias-file", default=None, help="追加人工别名 CSV: raw,symbol")
    p.add_argument("-o", "--outdir", default="np_out")
    p.set_defaults(fn=cmd_herbs)

    p = sub.add_parser("validate", help="靶点 accession↔gene 双向索引校验")
    p.add_argument("--input", required=True, help="靶点表 CSV")
    p.add_argument("--uniprot", default=os.path.join(DATA, "protein_symbols.tsv"))
    p.add_argument("-o", "--outdir", default="np_out")
    p.set_defaults(fn=cmd_validate)

    p = sub.add_parser("selftest", help="内置自测")
    p.set_defaults(fn=cmd_selftest)

    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
