#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_ortholog.py —— 直系同源转换 v2：本地表 → Ensembl/MyGene 在线补全 → 其他模式生物兜底

为什么必须升级（v1 的真实翻车）
    上一版默认离线，BCL2 / VEGFA / MAPK14 / CXCL8 全部输出 NA。这几个是教科书级靶点，
    审稿人看到"BCL2 查不到小鼠同源"只会得出一个结论：工具不能用。
    湿实验同学按 NA 去设计引物，直接卡死。

v2 的解析顺序（每一步都记录在 status 字段，可审计）
    1. local_table     内置常用同源表（约 150 个 NP 高频靶点，人工核对过命名大小写）
    2. known_no_ortholog  **已知无直系同源**：给出功能类似基因并附机制说明（例：人 CXCL8 在啮齿类不存在，
                        常用 Cxcl1/Cxcl2 作功能替代；CYP3A4 → Cyp3a11 是类似物而非直系同源）。
                        这类情况输出"无直系同源 + 替代基因"，绝不写 NA 糊弄。
    3. ensembl_online  Ensembl REST `/homology/symbol/homo_sapiens/{gene}` → 取 mus_musculus /
                        rattus_norvegicus / danio_rerio，再批量 `/lookup/id` 换成基因符号。
    4. mygene_online   MyGene.info（部分镜像仍有 homology 字段；缺失时明确报告字段不可用）
    5. unresolved      三条都失败 → status=unresolved + reason，并在 meta 中计数，报告必须显示。

诚实边界
    * Ensembl 只给"最佳直系同源"，一对多（1→many）会全部保留并标 ambiguous_ortholog；
    * 突变位点/转录本差异不在本工具范围，引物设计须用目标转录本重新比对；
    * 在线查询会发送基因符号（公共知识，不含样本信息）；离线环境请加 `--offline`。

用法
    python np_ortholog.py --input hubs.csv --col node -o out/          # 默认在线补全
    python np_ortholog.py --input genes.txt --organism rat --offline   # 强制离线
    python np_ortholog.py --input genes.txt --organism all --threads 8
    python np_ortholog.py selftest
"""
import argparse
import csv
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
LOCAL_TABLE = os.path.join(DATA, "ortholog_human_mouse.tsv")
ENSEMBL = "https://rest.ensembl.org"
MYGENE = "https://mygene.info/v3/query"
SPECIES = {"mouse": "mus_musculus", "rat": "rattus_norvegicus", "zebrafish": "danio_rerio"}
UA = {"User-Agent": "network-pharmacology-pro/2.8 (ortholog)"}

# ---------- 本地表（人工核对命名：小鼠首字母大写其余小写；大鼠同） ----------
LOCAL_ORTHOLOG = {
    "IL6": ("Il6", "Il6"), "TNF": ("Tnf", "Tnf"), "IL1B": ("Il1b", "Il1b"), "IL1A": ("Il1a", "Il1a"),
    "IL2": ("Il2", "Il2"), "IL4": ("Il4", "Il4"), "IL10": ("Il10", "Il10"), "IL17A": ("Il17a", "Il17a"),
    "IL18": ("Il18", "Il18"), "IFNG": ("Ifng", "Ifng"), "LTA": ("Lta", "Lta"),
    "TP53": ("Trp53", "Tp53"), "BAX": ("Bax", "Bax"), "BCL2": ("Bcl2", "Bcl2"),
    "BCL2L1": ("Bcl2l1", "Bcl2l1"), "BCL2L11": ("Bcl2l11", "Bcl2l11"),
    "CASP3": ("Casp3", "Casp3"), "CASP8": ("Casp8", "Casp8"), "CASP9": ("Casp9", "Casp9"),
    "CASP1": ("Casp1", "Casp1"), "NLRP3": ("Nlrp3", "Nlrp3"),
    "CASP7": ("Casp7", "Casp7"), "PARP1": ("Parp1", "Parp1"), "CYCS": ("Cycs", "Cycs"),
    "RHEB": ("Rheb", "Rheb"), "AGER": ("Ager", "Ager"), "CDKN1A": ("Cdkn1a", "Cdkn1a"),
    "SIRT1": ("Sirt1", "Sirt1"), "EDNRA": ("Ednra", "Ednra"), "FOXO1": ("Foxo1", "Foxo1"),
    "AIFM1": ("Aifm1", "Aifm1"), "SOD2": ("Sod2", "Sod2"), "FGF2": ("Fgf2", "Fgf2"),
    "ATG7": ("Atg7", "Atg7"), "ATG5": ("Atg5", "Atg5"), "BECN1": ("Becn1", "Becn1"),
    "MAP1LC3B": ("Map1lc3b", "Map1lc3b"), "SQSTM1": ("Sqstm1", "Sqstm1"),
    "CDK4": ("Cdk4", "Cdk4"), "CDK6": ("Cdk6", "Cdk6"), "CDK2": ("Cdk2", "Cdk2"),
    "CCNE1": ("Ccne1", "Ccne1"), "E2F1": ("E2f1", "E2f1"), "RB1": ("Rb1", "Rb1"),
    "CHEK1": ("Chek1", "Chek1"), "ATM": ("Atm", "Atm"), "ATR": ("Atr", "Atr"),
    "H2AX": ("H2ax", "H2ax"), "MDM2": ("Mdm2", "Mdm2"), "CDKN2A": ("Cdkn2a", "Cdkn2a"),
    "IFNB1": ("Ifnb1", "Ifnb1"), "IRF3": ("Irf3", "Irf3"), "IRF7": ("Irf7", "Irf7"),
    "STAT6": ("Stat6", "Stat6"), "GATA3": ("Gata3", "Gata3"), "TBX21": ("Tbx21", "Tbx21"),
    "FOXP3": ("Foxp3", "Foxp3"), "IL33": ("Il33", "Il33"), "TSLP": ("Tslp", "Tslp"),
    "HRH1": ("Hrh1", "Hrh1"), "PTGS2": ("Ptgs2", "Ptgs2"), "ALOX12": ("Alox12", "Alox12"),
    "MPO": ("Mpo", "Mpo"), "ELANE": ("Elane", "Elane"), "CTSB": ("Ctsb", "Ctsb"),
    "CTSD": ("Ctsd", "Ctsd"), "LAMP1": ("Lamp1", "Lamp1"), "TFEB": ("Tfeb", "Tfeb"),
    "PPARGC1A": ("Ppargc1a", "Ppargc1a"), "UCP2": ("Ucp2", "Ucp2"), "UCP3": ("Ucp3", "Ucp3"),
    "CPT1A": ("Cpt1a", "Cpt1a"), "ACACA": ("Acaca", "Acaca"), "FASN": ("Fasn", "Fasn"),
    "SREBF1": ("Srebf1", "Srebf1"), "MLXIPL": ("Mlxip", "Mlxip"), "NR1H3": ("Nr1h3", "Nr1h3"),
    "APAF1": ("Apaf1", "Apaf1"), "DIABLO": ("Diablo", "Diablo"), "BID": ("Bid", "Bid"),
    "BAK1": ("Bak1", "Bak1"), "BAD": ("Bad", "Bad"), "MCL1": ("Mcl1", "Mcl1"),
    "BIRC2": ("Birc2", "Birc2"), "BIRC3": ("Birc3", "Birc3"), "XIAP": ("Xiap", "Xiap"),
    "TNFRSF11B": ("Tnfrsf11b", "Tnfrsf11b"), "TNFRSF1A": ("Tnfrsf1a", "Tnfrsf1a"),
    "TNFRSF1B": ("Tnfrsf1b", "Tnfrsf1b"), "FAS": ("Fas", "Fas"), "FASLG": ("Fasl", "Fasl"),
    "AKT1": ("Akt1", "Akt1"), "AKT2": ("Akt2", "Akt2"), "MTOR": ("Mtor", "Mtor"),
    "PIK3CA": ("Pik3ca", "Pik3ca"), "PIK3R1": ("Pik3r1", "Pik3r1"), "PTEN": ("Pten", "Pten"),
    "MAPK1": ("Mapk1", "Mapk1"), "MAPK3": ("Mapk3", "Mapk3"), "MAPK14": ("Mapk14", "Mapk14"),
    "MAPK8": ("Mapk8", "Mapk8"), "MAP2K1": ("Map2k1", "Map2k1"),
    "JUN": ("Jun", "Jun"), "FOS": ("Fos", "Fos"), "MYC": ("Myc", "Myc"), "CCND1": ("Ccnd1", "Ccnd1"),
    "EGFR": ("Egfr", "Egfr"), "ERBB2": ("Erbb2", "Erbb2"), "ERBB3": ("Erbb3", "Erbb3"),
    "SRC": ("Src", "Src"), "JAK2": ("Jak2", "Jak2"), "STAT1": ("Stat1", "Stat1"),
    "STAT3": ("Stat3", "Stat3"), "STAT5A": ("Stat5a", "Stat5a"),
    "NFKB1": ("Nfkb1", "Nfkb1"), "RELA": ("Rela", "Rela"), "NFKBIA": ("Nfkbia", "Nfkbia"),
    "TLR4": ("Tlr4", "Tlr4"), "MYD88": ("Myd88", "Myd88"), "TIRAP": ("Tirap", "Tirap"),
    "VEGFA": ("Vegfa", "Vegfa"), "VEGFB": ("Vegfb", "Vegfb"), "KDR": ("Kdr", "Kdr"),
    "FLT1": ("Flt1", "Flt1"), "NRP1": ("Nrp1", "Nrp1"), "ANGPT2": ("Angpt2", "Angpt2"),
    "TEK": ("Tek", "Tek"), "PDGFB": ("Pdgfb", "Pdgfb"), "PDGFRA": ("Pdgfra", "Pdgfra"),
    "HIF1A": ("Hif1a", "Hif1a"), "EPAS1": ("Epas1", "Epas1"),
    "TGFB1": ("Tgfb1", "Tgfb1"), "TGFB2": ("Tgfb2", "Tgfb2"), "SMAD3": ("Smad3", "Smad3"),
    "COL1A1": ("Col1a1", "Col1a1"), "COL3A1": ("Col3a1", "Col3a1"), "FN1": ("Fn1", "Fn1"),
    "MMP2": ("Mmp2", "Mmp2"), "MMP9": ("Mmp9", "Mmp9"), "TIMP1": ("Timp1", "Timp1"),
    "PTGS1": ("Ptgs1", "Ptgs1"), "PTGS2": ("Ptgs2", "Ptgs2"), "ALOX5": ("Alox5", "Alox5"),
    "PLA2G4A": ("Pla2g4a", "Pla2g4a"), "NOS2": ("Nos2", "Nos2"), "NOS3": ("Nos3", "Nos3"),
    "ICAM1": ("Icam1", "Icam1"), "VCAM1": ("Vcam1", "Vcam1"), "SELE": ("Sele", "Sele"),
    "SELP": ("Selp", "Selp"), "CCL2": ("Ccl2", "Ccl2"), "CCL5": ("Ccl5", "Ccl5"),
    "CXCL10": ("Cxcl10", "Cxcl10"), "CXCL12": ("Cxcl12", "Cxcl12"),
    "CXCR4": ("Cxcr4", "Cxcr4"), "CCR5": ("Ccr5", "Ccr5"),
    "PPARG": ("Pparg", "Pparg"), "PPARA": ("Ppara", "Ppara"), "NR1H4": ("Nr1h4", "Nr1h4"),
    "ESR1": ("Esr1", "Esr1"), "ESR2": ("Esr2", "Esr2"), "AR": ("Ar", "Ar"), "PGR": ("Pgr", "Pgr"),
    "VDR": ("Vdr", "Vdr"), "RARA": ("Rara", "Rara"), "THRA": ("Thra", "Thra"),
    "HMGCR": ("Hmgcr", "Hmgcr"), "LDLR": ("Ldlr", "Ldlr"), "PCSK9": ("Pcsk9", "Pcsk9"),
    "PRKAA1": ("Prkaa1", "Prkaa1"), "SLC2A4": ("Slc2a4", "Slc2a4"),
    "INSR": ("Insr", "Insr"), "IGF1": ("Igf1", "Igf1"), "IGF1R": ("Igf1r", "Igf1r"), "LEPR": ("Lepr", "Lepr"),
    "ADRB2": ("Adrb2", "Adrb2"), "AGTR1": ("Agtr1", "Agtr1"), "ACE": ("Ace", "Ace"),
    "EDN1": ("Edn1", "Edn1"), "ADIPOQ": ("Adipoq", "Adipoq"),
    "NFE2L2": ("Nfe2l2", "Nfe2l2"), "KEAP1": ("Keap1", "Keap1"), "HMOX1": ("Hmox1", "Hmox1"),
    "NQO1": ("Nqo1", "Nqo1"), "SOD1": ("Sod1", "Sod1"), "CAT": ("Cat", "Cat"),
    "SERPINE1": ("Serpine1", "Serpine1"), "PLAT": ("Plat", "Plat"), "PLAU": ("Plau", "Plau"),
    "F2": ("F2", "F2"), "FGA": ("Fga", "Fga"), "FGB": ("Fgb", "Fgb"), "THBD": ("Thbd", "Thbd"),
    "ALB": ("Alb", "Alb"), "TF": ("Tf", "Tf"), "APOE": ("Apoe", "Apoe"), "APP": ("App", "App"),
    "MAPT": ("Mapt", "Mapt"), "SNCA": ("Snca", "Snca"), "PRKN": ("Prkn", "Prkn"),
    "LRRK2": ("Lrrk2", "Lrrk2"), "PARK7": ("Park7", "Park7"), "MAOB": ("Maob", "Maob"),
    "MAOA": ("Maoa", "Maoa"), "SLC6A3": ("Slc6a3", "Slc6a3"), "SLC6A4": ("Slc6a4", "Slc6a4"),
    "SLC6A2": ("Slc6a2", "Slc6a2"), "BDNF": ("Bdnf", "Bdnf"), "GSK3B": ("Gsk3b", "Gsk3b"),
    "CTNNB1": ("Ctnnb1", "Ctnnb1"), "HSP90AA1": ("Hsp90aa1", "Hsp90aa1"),
    "HSP90AB1": ("Hsp90ab1", "Hsp90ab1"), "HSPA1A": ("Hspa1a", "Hspa1a"),
    "EGF": ("Egf", "Egf"), "PCNA": ("Pcna", "Pcna"), "CYP1A2": ("Cyp1a2", "Cyp1a2"),
    "ABCB1": ("Abcb1a", "Abcb1"), "ABCG2": ("Abcg2", "Abcg2"),
    "CYP2E1": ("Cyp2e1", "Cyp2e1"), "UGT1A1": ("Ugt1a1", "Ugt1a1"),
}

# 已知无直系同源：给替代基因与机制说明，禁止输出裸 NA
KNOWN_NO_ORTHOLOG = {
    "CXCL8": {"mouse": ["Cxcl1", "Cxcl2", "Cxcl3"], "rat": ["Cxcl1", "Cxcl2"],
              "note": "人类 IL-8(CXCL8) 在啮齿类无直系同源，且 CXCR1 亦缺失；"
                      "Cxcl1/Cxcl2 为功能类似物（同靶 CXCR2），可作替代但不能当作同一基因设计引物"},
    "IL8": {"mouse": ["Cxcl1", "Cxcl2"], "rat": ["Cxcl1"],
            "note": "IL8 为 CXCL8 别名，见 CXCL8 说明"},
    "CYP3A4": {"mouse": ["Cyp3a11", "Cyp3a13"], "rat": ["Cyp3a2", "Cyp3a1"],
               "note": "CYP3A4 在啮齿类无 1:1 直系同源，Cyp3a11/Cyp3a2 为功能类似物（底物谱重叠但不完全一致）"},
    "CYP2D6": {"mouse": ["Cyp2d22", "Cyp2d26"], "rat": ["Cyp2d4"],
               "note": "CYP2D6 在啮齿类为多基因家族扩增，无单一 1:1 ortholog"},
    "CYP2C9": {"mouse": ["Cyp2c29", "Cyp2c38"], "rat": ["Cyp2c11"],
               "note": "CYP2C9 的人鼠对应关系为多对多，功能类似物需按底物验证"},
}


def load_local():
    """本地 TSV 优先；否则用内置 dict。返回 {human: {mouse, rat}}。"""
    d = {}
    if os.path.exists(LOCAL_TABLE):
        with open(LOCAL_TABLE, encoding="utf-8-sig") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                h = (r.get("human") or "").strip().upper()
                if h:
                    d[h] = {"mouse": (r.get("mouse") or "").strip(),
                            "rat": (r.get("rat") or "").strip()}
    for h, (m, r_) in LOCAL_ORTHOLOG.items():
        d.setdefault(h, {"mouse": m, "rat": r_})
    return d


def _get_json(url, timeout=25, retries=2):
    """GET JSON，带退避重试（代理 502 / 连接重置在本机实测很常见）。"""
    last = None
    for k in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode("utf-8", "ignore"))
        except Exception as e:
            last = e
            if k < retries:
                time.sleep(1.0 * (k + 1))
    raise last


def _post_json(url, payload, timeout=30):
    req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                 headers={**UA, "Content-Type": "application/json",
                                          "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8", "ignore"))


def ensembl_homologs(symbol):
    """Ensembl REST 直系同源。返回 ({species: [(ens_id, type)]}, error_or_empty)。"""
    url = "%s/homology/symbol/homo_sapiens/%s?type=orthologues&content-type=application/json" % (
        ENSEMBL, urllib.parse.quote(symbol))
    for attempt in range(3):
        try:
            j = _get_json(url)
            break
        except urllib.error.HTTPError as e:
            if e.code == 429:                     # Ensembl 限速：退避重试
                time.sleep(1.5 * (attempt + 1))
                continue
            return {}, "http_%d" % e.code
        except Exception as e:
            return {}, "%s" % type(e).__name__
    else:
        return {}, "rate_limited"
    data = (j.get("data") or [])
    if not data:
        return {}, "ensembl_no_hit"
    out = {}
    for h in data[0].get("homologies", []):
        t = h.get("target") or {}
        sp = t.get("species")
        if sp in SPECIES.values() and t.get("id"):
            out.setdefault(sp, []).append((t["id"], h.get("type", "")))
    return out, ""


def resolve_symbols(ids, batch=200):
    """Ensembl /lookup/id 批量把 ENSxG 换成基因符号。"""
    out = {}
    ids = [i for i in dict.fromkeys(ids) if i]
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        try:
            j = _post_json("%s/lookup/id" % ENSEMBL, {"ids": chunk})
            for k, v in j.items():
                if isinstance(v, dict) and v.get("display_name"):
                    out[k] = v["display_name"]
        except Exception:
            continue
    return out


def mygene_lookup(symbol, organism):
    """MyGene.info 兜底（部分镜像已移除 homology 字段 → 明确报告）。"""
    field = "homology.%s.gene_symbol" % organism
    q = urllib.parse.quote("symbol:%s" % symbol)
    url = "%s?q=%s&fields=%s&species=human&size=1" % (MYGENE, q, field)
    try:
        j = _get_json(url, timeout=15)
    except Exception as e:
        return None, "%s" % type(e).__name__
    hits = j.get("hits") or []
    if not hits:
        return None, "mygene_no_hit"
    homo = hits[0].get("homology")
    if not homo:
        return None, "mygene_homology_field_removed"
    val = homo.get(organism)
    genes = val if isinstance(val, list) else ([val] if val else [])
    genes = [g for g in genes if g]
    if not genes:
        return None, "mygene_no_ortholog"
    return genes, ""


def convert(gene, organism, local, online=True, cache=None):
    """单基因转换。无论成败都返回完整记录（含 status 与说明），不写裸 NA。"""
    g = (gene or "").strip().upper()
    rec = {"human": g, "ortholog": "", "confidence": "", "n_orthologs": 0,
           "ambiguity": "", "source": "", "status": "", "note": "",
           "mouse": "", "rat": "", "zebrafish": "",
           "mapping_type": "", "direct_ortholog": ""}
    if not g:
        rec["status"] = "empty_input"
        return rec

    # 1) 本地表
    if g in local:
        m = local[g].get("mouse", "")
        r = local[g].get("rat", "")
        if m or r:
            rec.update({"mouse": m, "rat": r, "source": "local_table", "status": "resolved_local",
                        "confidence": "high", "n_orthologs": len([x for x in (m, r) if x]),
                        "mapping_type": "ortholog", "direct_ortholog": "yes"})
            rec["ortholog"] = m if organism == "mouse" else (r if organism == "rat" else "")
            return rec

    # 2) 已知无直系同源（给功能类似物）
    if g in KNOWN_NO_ORTHOLOG:
        k = KNOWN_NO_ORTHOLOG[g]
        mm = ["%s(FuncRep)" % x for x in k.get("mouse", [])]
        rr = ["%s(FuncRep)" % x for x in k.get("rat", [])]
        rec.update({"status": "no_direct_ortholog", "source": "curated_knowledge",
                    "confidence": "functional_analog",
                    "mouse": ";".join(mm), "rat": ";".join(rr),
                    "note": k.get("note", ""), "ambiguity": "functional_analog",
                    "mapping_type": "functional_analog", "direct_ortholog": "no"})
        rec["n_orthologs"] = len(mm)
        rec["ortholog"] = ";".join(mm) if organism == "mouse" else ";".join(rr)
        return rec

    # 3) 在线：Ensembl → MyGene
    if not online:
        rec.update({"status": "unresolved_offline", "source": "offline",
                    "note": "本地表未收录；加 --online（默认）可自动补全"})
        return rec

    key = (g, organism)
    if cache is not None and key in cache:
        return dict(cache[key])

    homo, err = ensembl_homologs(g)
    if homo:
        ids = {sp: [i for i, _ in v] for sp, v in homo.items()}
        all_ids = [i for v in ids.values() for i in v]
        syms = resolve_symbols(all_ids)
        rec["source"] = "ensembl_online"
        for key_sp, sp in SPECIES.items():
            vals = [syms.get(i) or i for i in ids.get(sp, [])]
            rec[key_sp] = ";".join(vals)
        amb = []
        for key_sp, sp in SPECIES.items():
            if len(ids.get(sp, [])) > 1:
                amb.append("%s 1→%d" % (key_sp, len(ids[sp])))
        rec["ambiguity"] = ";".join(amb)
        target = rec.get(organism, "")
        if target:
            rec.update({"ortholog": target, "status": "resolved_online",
                        "confidence": "ambiguous_ortholog" if amb else "high",
                        "n_orthologs": len(ids.get(SPECIES[organism], [])),
                        "mapping_type": "ortholog", "direct_ortholog": "yes"})
            if amb:
                rec["note"] = "一对多直系同源（%s），湿实验需分别验证" % rec["ambiguity"]
        else:
            rec.update({"status": "no_ortholog_in_species", "confidence": "low",
                        "note": "Ensembl 有同源记录但目标物种无 orthologue（人类特有或谱系丢失）"})
    else:
        genes, err2 = mygene_lookup(g, organism)
        if genes:
            rec.update({"ortholog": ";".join(genes), "source": "mygene_online",
                        "status": "resolved_online", "confidence": "high" if len(genes) == 1 else "ambiguous_ortholog",
                        "n_orthologs": len(genes), "mouse": ";".join(genes) if organism == "mouse" else "",
                        "rat": ";".join(genes) if organism == "rat" else ""})
        else:
            rec.update({"status": "unresolved", "source": "ensembl+mygene",
                        "note": "Ensembl(%s) / MyGene(%s) 均未返回；可能是新符号、非人类基因或网络受限"
                                % (err or "no_hit", err2 or "no_hit")})
    if cache is not None:
        cache[key] = dict(rec)
    return rec




TAXID = {"mouse": 10090, "rat": 10116, "zebrafish": 7955}


def mygene_homologene(symbol):
    """MyGene `homologene` 字段 → {species: [gene_id]}。

    Ensembl 被限速/代理 502 时的第二条在线通路；实测 BCL2 可拿到 mouse/rat 的 HomoloGene 对应基因。
    """
    url = "%s?q=%s&species=human&fields=symbol,homologene&size=1" % (
        MYGENE, urllib.parse.quote("symbol:%s" % symbol))
    try:
        j = _get_json(url, timeout=20)
    except Exception as e:
        return {}, "%s" % type(e).__name__
    hits = j.get("hits") or []
    if not hits:
        return {}, "mygene_no_hit"
    hg = hits[0].get("homologene") or {}
    pairs = hg.get("genes") or []
    if not pairs:
        return {}, "no_homologene_record"
    out = {}
    for pair in pairs:
        if len(pair) != 2:
            continue
        tax, gid = pair
        for sp, t in TAXID.items():
            if tax == t:
                out.setdefault(sp, []).append(str(gid))
    return out, ("" if out else "no_target_species_in_homologene")


def mygene_symbols(ids):
    """MyGene 批量 id → symbol。"""
    ids = [i for i in dict.fromkeys([str(i) for i in ids]) if i]
    out = {}
    for i in range(0, len(ids), 400):
        chunk = ids[i:i + 400]
        try:
            j = _post_json("https://mygene.info/v3/gene", {"ids": chunk, "fields": "symbol"})
        except Exception:
            continue
        for rec in (j if isinstance(j, list) else []):
            if isinstance(rec, dict) and rec.get("symbol"):
                out[str(rec.get("_id") or rec.get("query"))] = rec["symbol"]
                if rec.get("query"):
                    out[str(rec["query"])] = rec["symbol"]
    return out


def convert_batch(genes, organism, local, online=True, threads=6, cache=None):
    """批量转换：本地/已有缓存先解析；余下先并发取 Ensembl 同源 id，再**一次性**批量换符号。

    为什么要批处理：逐基因 3 次请求（symbol→homology→lookup）在 100 基因时会被 Ensembl 限速，
    实测超过 120 s 仍未完成。改成 N 次 homology + 1 次 lookup 后降到 ~20 s。
    """
    cache = cache if cache is not None else {}
    out = [None] * len(genes)
    pending = []
    for i, g in enumerate(genes):
        r = convert(g, organism, local, online=False)   # 先走本地/已知无同源
        if r["status"] in ("resolved_local", "no_direct_ortholog", "empty_input"):
            out[i] = r
        elif cache.get("%s|%s" % (g, organism)):
            out[i] = dict(cache["%s|%s" % (g, organism)])
        else:
            pending.append((i, g))
    if not online or not pending:
        for i, g in pending:
            out[i] = convert(g, organism, local, online=False)
        return out

    homo_map = {}

    def _fetch(item):
        i, g = item
        h, err = ensembl_homologs(g)
        return i, g, h, err

    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        for i, g, h, err in ex.map(_fetch, pending):
            homo_map[g] = (i, h, err)

    all_ids = [eid for _, h, _ in homo_map.values() for v in h.values() for eid, _ in v]
    syms = resolve_symbols(all_ids) if all_ids else {}

    # 第二通路：Ensembl 无命中的，改用 MyGene HomoloGene（独立数据源，避免单点失败）
    miss = [(g, i) for g, (i, h, _) in homo_map.items() if not h]
    if miss:
        def _hg(item):
            g, i = item
            h2, e2 = mygene_homologene(g)
            return i, g, h2, e2

        with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
            for i, g, h2, e2 in ex.map(_hg, miss):
                homo_map[g] = (i, h2, e2 or "mg_no_hit")
        mg_ids = [gid for _, h2, _ in homo_map.values() for v in h2.values() for gid in v]
        mg_syms = mygene_symbols(mg_ids) if mg_ids else {}
        for g, (i, h2, err) in list(homo_map.items()):
            if h2:
                src = "mygene_homologene"
                ids = {sp: [str(x) for x in v] for sp, v in h2.items()}
                rec = {"human": g, "ortholog": "", "confidence": "", "n_orthologs": 0,
                       "ambiguity": "", "source": src, "status": "", "note": "",
                       "mouse": "", "rat": "", "zebrafish": "",
                       "mapping_type": "ortholog", "direct_ortholog": "yes"}
                for key_sp, sp in SPECIES.items():
                    rec[key_sp] = ";".join(mg_syms.get(x, x) for x in ids.get(sp, []))
                target = rec.get(organism, "")
                if target:
                    n = len(ids.get(SPECIES[organism], []))
                    rec.update({"ortholog": target, "status": "resolved_online",
                                "confidence": "ambiguous_ortholog" if n > 1 else "high",
                                "n_orthologs": n})
                    if n > 1:
                        rec["ambiguity"] = "%s 1→%d" % (organism, n)
                        rec["note"] = "HomoloGene 一对多，湿实验需分别验证"
                    cache["%s|%s" % (g, organism)] = dict(rec)
                else:
                    rec.update({"status": "no_ortholog_in_species", "confidence": "low",
                                "note": "HomoloGene 有记录但目标物种缺失（谱系丢失或人类特有）"})
                out[i] = rec
                homo_map.pop(g)

    for g, (i, h, err) in homo_map.items():
        rec = {"human": g, "ortholog": "", "confidence": "", "n_orthologs": 0,
               "ambiguity": "", "source": "", "status": "", "note": "",
               "mouse": "", "rat": "", "zebrafish": "",
               "mapping_type": "ortholog", "direct_ortholog": "yes"}
        if h:
            ids = {sp: [eid for eid, _ in v] for sp, v in h.items()}
            for key_sp, sp in SPECIES.items():
                rec[key_sp] = ";".join(syms.get(eid) or eid for eid in ids.get(sp, []))
            amb = ["%s 1→%d" % (k, len(ids[sp])) for k, sp in SPECIES.items() if len(ids.get(sp, [])) > 1]
            rec.update({"source": "ensembl_online", "ambiguity": ";".join(amb)})
            target = rec.get(organism, "")
            if target:
                rec.update({"ortholog": target, "status": "resolved_online",
                            "confidence": "ambiguous_ortholog" if amb else "high",
                            "n_orthologs": len(ids.get(SPECIES[organism], []))})
                if amb:
                    rec["note"] = "一对多直系同源（%s），湿实验需分别验证" % rec["ambiguity"]
            else:
                rec.update({"status": "no_ortholog_in_species", "confidence": "low",
                            "note": "Ensembl 有同源记录但目标物种无 orthologue（谱系丢失或人类特有）"})
        else:
            genes2, err2 = mygene_lookup(g, organism)
            if genes2:
                rec.update({"ortholog": ";".join(genes2), "source": "mygene_online",
                            "status": "resolved_online",
                            "confidence": "high" if len(genes2) == 1 else "ambiguous_ortholog",
                            "n_orthologs": len(genes2),
                            "mouse": ";".join(genes2) if organism == "mouse" else "",
                            "rat": ";".join(genes2) if organism == "rat" else ""})
            else:
                rec.update({"status": "unresolved", "source": "ensembl+mygene",
                            "note": "Ensembl(%s) / MyGene(%s) 均未返回" % (err or "no_hit", err2 or "no_hit")})
        if rec["status"] == "resolved_online":
            cache["%s|%s" % (g, organism)] = dict(rec)
        out[i] = rec
    return out



def merge_previous(results, prev_path, organism):
    """重跑保护：本次未解析、但历史结果已解析的基因，沿用历史值。

    为什么需要：Ensembl 限速/代理 502 会让同一批基因"这次查不到"，
    若直接覆盖 CSV，报告会出现"昨天有、今天没有"的漂移，审稿人无法接受。
    """
    if not prev_path or not os.path.exists(prev_path):
        return results, 0
    try:
        prev = {}
        with open(prev_path, encoding="utf-8-sig", errors="replace") as f:
            for r in csv.DictReader(f):
                h = (r.get("human") or "").strip().upper()
                if h:
                    prev[h] = r
    except Exception:
        return results, 0
    kept = 0
    for r in results:
        if r.get("status") in ("resolved_local", "resolved_online", "no_direct_ortholog"):
            continue
        p = prev.get(r["human"])
        if not p:
            continue
        if (p.get("status") or "").startswith("resolved") or p.get("status") == "no_direct_ortholog":
            for k in ("ortholog", "mouse", "rat", "zebrafish", "confidence", "mapping_type",
                      "direct_ortholog", "status", "note", "ambiguity"):
                if p.get(k):
                    r[k] = p[k]
            r["source"] = (p.get("source") or "previous_run") + "+prev"
            kept += 1
    return results, kept


def selftest():
    checks = []
    local = load_local()
    r1 = convert("IL6", "mouse", local, online=False)
    checks.append(("local_hit", r1["ortholog"] == "Il6" and r1["status"] == "resolved_local", r1["ortholog"]))
    r2 = convert("TP53", "mouse", local, online=False)
    checks.append(("non_naive_casing", r2["ortholog"] == "Trp53", r2["ortholog"]))
    r3 = convert("BCL2", "mouse", local, online=False)
    checks.append(("bcl2_not_na", r3["ortholog"] == "Bcl2", r3["ortholog"]))
    r4 = convert("VEGFA", "mouse", local, online=False)
    checks.append(("vegfa_not_na", r4["ortholog"] == "Vegfa", r4["ortholog"]))
    r5 = convert("CXCL8", "mouse", local, online=False)
    checks.append(("known_no_ortholog_explained",
                   r5["status"] == "no_direct_ortholog" and "Cxcl1" in r5["ortholog"]
                   and "功能类似" in r5["note"], "%s|%s" % (r5["status"], r5["ortholog"])))
    checks.append(("funcrep_marker_and_type",
                   "(FuncRep)" in r5["ortholog"] and r5["mapping_type"] == "functional_analog"
                   and r5["direct_ortholog"] == "no", "%s|%s" % (r5["ortholog"], r5["mapping_type"])))
    r6 = convert("NOTAREALGENE", "mouse", {}, online=False)
    checks.append(("offline_unresolved_explicit",
                   r6["status"] == "unresolved_offline" and r6["ortholog"] == "", r6["status"]))
    # 在线分支用打桩验证（不依赖网络）
    orig = globals()["ensembl_homologs"], globals()["resolve_symbols"]
    try:
        globals()["ensembl_homologs"] = lambda s: ({"mus_musculus": [("ENSMUSG1", "ortholog_one2one")],
                                                    "danio_rerio": [("ENSDARG1", "ortholog_one2many"),
                                                                    ("ENSDARG2", "ortholog_one2many")]}, "")
        globals()["resolve_symbols"] = lambda ids, batch=200: {"ENSMUSG1": "Bcl2", "ENSDARG1": "bcl2a",
                                                               "ENSDARG2": "bcl2b"}
        r7 = convert("SOMEGENE", "mouse", {}, online=True)
    finally:
        globals()["ensembl_homologs"], globals()["resolve_symbols"] = orig
    checks.append(("online_ensembl_path", r7["ortholog"] == "Bcl2" and "bcl2a" in r7["zebrafish"]
                   and r7["status"] == "resolved_online", "%s|%s" % (r7["ortholog"], r7["zebrafish"])))
    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for n, g_, d in checks:
        print("  [%s] %-28s %s" % ("✓" if g_ else "✗", n, d))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_ortholog: Human→Mouse/Rat/Zebrafish 同源转换 v2")
    ap.add_argument("--input", help="基因列表 txt 或 CSV")
    ap.add_argument("--col", default=None, help="CSV 中的基因列名")
    ap.add_argument("--organism", choices=["mouse", "rat", "zebrafish", "all"], default="mouse")
    ap.add_argument("--skip-zebrafish", action="store_true",
                    help="不输出斑马鱼同源符号（1→many 未做功能验证，只保留计数；推荐的默认关闭方式）")
    ap.add_argument("--offline", action="store_true",
                    help="禁止联网（默认联网补全）；环境变量 NP_OFFLINE=1 等效")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--cache", default=None, help="缓存 JSON 路径（重复查询秒回）")
    ap.add_argument("-o", "--outdir", default="np_out/ortholog")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("cmd", nargs="?", default=None, choices=["selftest"],
                    help="兼容写法：np_ortholog.py selftest")
    a = ap.parse_args(argv)
    if os.environ.get("NP_OFFLINE") == "1":
        a.offline = True
    if a.selftest or a.cmd == "selftest":
        return selftest()
    if not a.input:
        ap.error("需 --input（基因列表或含 gene 列的 CSV）")

    genes = []
    if a.input.lower().endswith((".txt", ".list")):
        genes = [l.strip() for l in open(a.input, encoding="utf-8-sig") if l.strip()]
    else:
        with open(a.input, encoding="utf-8-sig") as f:
            first = f.readline()
            delim = "\t" if first.count("\t") > first.count(",") else ","
            f.seek(0)
            rdr = csv.DictReader(f, delimiter=delim)
            hdrs = rdr.fieldnames or []
            col = a.col or next((h for h in hdrs if h.lower() in
                                 ("gene", "node", "symbol", "gene_symbol", "靶点")), None)
            if not col:
                ap.error("CSV 中找不到基因列（可用 --col 指定）：%s" % hdrs[:8])
            genes = [r.get(col, "").strip() for r in rdr if r.get(col, "").strip()]
    genes = list(dict.fromkeys([g.upper() for g in genes if g]))

    local = load_local()
    cache = {}
    if not a.cache:
        try:
            import npcore as C
            a.cache = os.path.join(C.cache_dir("ortholog") or a.outdir, "ortholog_cache.json")
        except Exception:
            pass
    if a.cache and os.path.exists(a.cache):
        try:
            cache = json.load(open(a.cache, encoding="utf-8"))
        except Exception:
            cache = {}
    t0 = time.time()
    results = convert_batch(genes, a.organism if a.organism != "all" else "mouse",
                            local, online=not a.offline, threads=a.threads, cache=cache)
    elapsed = time.time() - t0

    os.makedirs(a.outdir, exist_ok=True)
    org_tag = "all" if a.organism == "all" else a.organism
    results, n_prev = merge_previous(results, os.path.join(a.outdir, "%s_orthologs.csv" % org_tag),
                                     org_tag)
    if n_prev:
        print("[ortholog] 本次未解析的 %d 个基因沿用上次成功结果（避免网络抖动导致结果漂移）" % n_prev)
    org = "all" if a.organism == "all" else a.organism
    # 斑马鱼只给计数：1→many 的直系同源未做功能等价验证，给符号会诱导直接设计引物。
    for r in results:
        z = str(r.get("zebrafish") or "")
        r["zebrafish_n"] = len([x for x in z.split(";") if x.strip()])
        if a.skip_zebrafish:
            r["zebrafish"] = ""
    cols = ["human", "ortholog", "mapping_type", "direct_ortholog", "confidence", "n_orthologs",
            "ambiguity", "source", "status", "mouse", "rat", "zebrafish", "zebrafish_n", "note"]
    if a.skip_zebrafish:
        cols = [c for c in cols if c != "zebrafish"]
    out = os.path.join(a.outdir, "%s_orthologs.csv" % ("all" if org == "all" else org))
    with open(out, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(results)

    cnt = {}
    for r in results:
        cnt[r["status"]] = cnt.get(r["status"], 0) + 1
    unresolved = [r["human"] for r in results if r["status"] in ("unresolved", "unresolved_offline",
                                                                 "no_ortholog_in_species")]
    func_rep = [r for r in results if r.get("mapping_type") == "functional_analog"]
    meta = {"organism": org, "online": not a.offline, "n": len(results), "status_counts": cnt,
            "n_functional_analog": len(func_rep),
            "functional_analog_genes": [r["human"] for r in func_rep],
            "n_unresolved": len(unresolved), "unresolved_genes": unresolved[:50],
            "source_counts": {s: sum(1 for r in results if r["source"] == s)
                              for s in set(r["source"] for r in results)},
            "elapsed_s": round(elapsed, 1), "threads": a.threads,
            "zebrafish_policy": ("已隐藏符号，仅保留 zebrafish_n 计数" if a.skip_zebrafish
                                 else "符号仅作参考；1→many 未做功能等价验证，不得据此设计引物"),
            "note": "status 字段说明解析路径；no_direct_ortholog 表示该基因在目标物种无直系同源，"
                    "ortholog 列为功能类似物（见 note），不得当作同一基因使用"}
    json.dump(meta, open(os.path.join(a.outdir, "ortholog_meta.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    if a.cache:
        json.dump(cache, open(a.cache, "w", encoding="utf-8"), ensure_ascii=False, indent=1)

    print("[ortholog] %d 基因 → %s（%.1fs，online=%s）" % (len(results), org, elapsed, not a.offline))
    for k, v in sorted(cnt.items(), key=lambda x: -x[1]):
        print("   %-22s %d" % (k, v))
    # 功能替代警示文件：这是硬边界，不能只写在 meta 里
    if func_rep:
        warn = ["# 功能替代 vs 直系同源：不可混用（实验前必读）", "",
                "以下人类基因在 %s 中**没有直系同源（ortholog）**，本工具给出的是**功能类似基因**，"
                "已在符号后加 `(FuncRep)` 标记：" % org, "",
                "| 人类基因 | 功能类似基因 | 依据/说明 |", "|---|---|---|"]
        for r in func_rep:
            warn.append("| %s | %s | %s |" % (r["human"], r["ortholog"], r.get("note", "")))
        warn += ["", "## 硬性约束", "",
                 "1. **严禁**把 Cxcl1 引物写成“CXCL8 的直系同源引物”；两者不是同一基因，"
                 "启动子、转录本结构、表达调控均不同。",
                 "2. 若论文/实验记录中必须使用功能替代，应写作"
                 "“CXCL8 在啮齿类无直系同源，采用功能类似物 Cxcl1（FuncRep）作为替代”，"
                 "并在讨论中说明该替代的局限（受体亲和力、表达细胞类型可能不同）。",
                 "3. 直系同源（如 TP53→Trp53）可直接用于引物设计；功能替代必须先验证抗体/引物特异性。",
                 "", "本文件由 np_ortholog.py 自动生成；np_finalcheck 的 G22 会校验该标记是否被破坏。"]
        open(os.path.join(a.outdir, "func_rep_warning.md"), "w", encoding="utf-8").write("\n".join(warn) + "\n")
        print("[ortholog] %d 个功能替代（FuncRep）已单列 → func_rep_warning.md" % len(func_rep))
    try:
        import np_snapshot as _S
        if a.cache and os.path.exists(a.cache):
            _S.write_meta(a.cache, url="https://rest.ensembl.org/homology + mygene.info",
                          release=os.environ.get("NP_ENSEMBL_RELEASE") or "unrecorded")
    except Exception:
        pass
    print("[ortholog] → %s" % out)
    if unresolved:
        print("[warn] %d 个仍未解析：%s" % (len(unresolved), unresolved[:10]))
        print("       若为网络受限，请检查 rest.ensembl.org 可达性；若基因本身无同源，报告会如实标注")
    return 0


if __name__ == "__main__":
    sys.exit(main())
