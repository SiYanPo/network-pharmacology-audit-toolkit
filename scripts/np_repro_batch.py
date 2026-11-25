#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_repro_batch.py —— v1.7.0 批量复现引擎（500 篇语料 → TCMSP 重算 → 差异归因）

把 v1.3 的"逐篇手工归因"升级为流水线：
  1. 读取 np_corpus extract 的方法表（药味/阈值/报告数量）
  2. 用药材 TCMSP 缓存（cn 名索引）重算：活性成分数（按论文阈值或默认 30%/0.18）
     → 成分并集 → 靶点并集（去重）
  3. 与原文报告值对比：比值 + 判定（一致 ±20% / 偏低 / 偏高 / 无法复现）
  4. 归因标签：药材未收录（覆盖率<1）、阈值差异（披露值≠默认）、靶点源差异（报告靶点
     远大于 TCMSP 实测）、原文未披露
  5. 汇总：可复现率、比值中位数、失败模式分布（按疾病/期刊分层）

命令：
  python np_repro_batch.py run --methods methods.csv --herb-cache HERBDIR --out repro.csv
  python np_repro_batch.py stats --in repro.csv --outdir stats
"""
import argparse
import glob
import json
import os
import re

import ast
import time
import pandas as pd
import requests

HERB_ALIAS = {  # 常见简写/异名 → 缓存拼音键（按需扩展）
    "炙甘草": "gancao", "生甘草": "gancao", "甘草梢": "gancao",
    "麸炒白术": "baizhu", "炒白术": "baizhu", "生白术": "baizhu",
    "酒当归": "danggui", "当归尾": "danggui", "全当归": "danggui",
    "赤芍": "baishao", "白芍": "baishao", "炒白芍": "baishao",
    "生地黄": "shudihuang", "熟地黄": "shudihuang", "地黄": "shudihuang",
    "川牛膝": "niuxi", "怀牛膝": "niuxi", "牛膝": "niuxi",
    "法半夏": "banxia", "姜半夏": "banxia", "清半夏": "banxia",
    "陈皮": "chenpi", "橘皮": "chenpi", "广陈皮": "chenpi",
    "山茱萸": "shanzhuyu", "山萸肉": "shanzhuyu",
    "生黄芪": "huangqi", "炙黄芪": "huangqi", "黄芪": "huangqi",
    "紫苏叶": "zisu", "紫苏梗": "zisu", "紫苏子": "zisuzi",
    "川芎": "chuanxiong", "抚芎": "chuanxiong",
    "丹参": "danshen", "紫丹参": "danshen",
    "茯苓": "fuling", "白茯苓": "fuling",
    "泽泻": "zexie", "建泽泻": "zexie",
    "牡丹皮": "mudanpi", "丹皮": "mudanpi",
    "桂枝": "guizhi", "肉桂": "rougui",
    "麻黄": "mahuang", "炙麻黄": "mahuang",
    "杏仁": "xingren", "苦杏仁": "xingren", "炒杏仁": "xingren",
    "黄芩": "huangqin", "枯黄芩": "huangqin",
    "黄连": "huanglian", "川黄连": "huanglian",
    "黄柏": "huangbai", "川黄柏": "huangbai",
    "大黄": "dahuang", "生大黄": "dahuang", "酒大黄": "dahuang",
    "金银花": "jinyinhua", "银花": "jinyinhua",
    "连翘": "lianqiao", "青连翘": "lianqiao",
    "柴胡": "chaihu", "北柴胡": "chaihu", "南柴胡": "chaihu",
    "栀子": "zhizi", "山栀子": "zhizi",
    "葛根": "gegen", "粉葛": "gegen",
    "薄荷": "bohe", "苏薄荷": "bohe",
    "白术": "baizhu", "苍术": "cangzhu",
    "党参": "dangshen", "潞党参": "dangshen",
    "人参": "renshen", "红参": "renshen", "生晒参": "renshen",
    "五味子": "wuweizi", "北五味子": "wuweizi",
    "麦冬": "maidong", "麦门冬": "maidong",
    "桃仁": "taoren", "燀桃仁": "taoren",
    "红花": "honghua", "草红花": "honghua",
    "益母草": "yimucao", "茺蔚": "yimucao",
    "生姜": "shengjiang", "干姜": "ganjiang", "炮姜": "ganjiang",
    "附子": "fuzi", "制附子": "fuzi", "黑顺片": "fuzi",
    "半夏": "banxia", "厚朴": "houpo", "枳实": "zhishi", "枳壳": "zhiqiao",
    "桔梗": "jiegeng", "苦桔梗": "jiegeng",
    "酸枣仁": "suanzaoren", "炒酸枣仁": "suanzaoren",
    "远志": "yuanzhi", "制远志": "yuanzhi",
    "石菖蒲": "shichangpu", "菖蒲": "shichangpu",
    "枸杞子": "gouqizi", "枸杞": "gouqizi",
    "山药": "shanyao", "淮山药": "shanyao",
    "薏苡仁": "yiyiren", "生薏苡仁": "yiyiren",
    "车前子": "cheqianzi", "盐车前子": "cheqianzi",
    "土茯苓": "tufuling", "绵萆薢": "bixie", "萆薢": "bixie",
    "虎杖": "huzhang", "茵陈": "yinchen", "绵茵陈": "yinchen",
    "金钱草": "jinqiancao", "广金钱草": "jinqiancao",
    "三七": "sanqi", "田七": "sanqi",
    "延胡索": "yanhusuo", "元胡": "yanhusuo",
    "郁金": "yujin", "姜黄": "jianghuang",
    "乳香": "ruxiang", "没药": "moyao",
    "鸡血藤": "jixueteng", "牛膝": "niuxi",
    "香附": "xiangfu", "醋香附": "xiangfu",
    "木香": "muxiang", "砂仁": "sharen",
    "藿香": "huoxiang", "佩兰": "peilan",
    "猪苓": "zhuling", "滑石": "huashi",
    "牡蛎": "muli", "煅牡蛎": "muli", "龙骨": "longgu",
    "杜仲": "duzhong", "续断": "xuduan", "桑寄生": "sangjisheng",
    "淫羊藿": "yinyanghuo", "仙灵脾": "yinyanghuo",
    "补骨脂": "buguzhi", "菟丝子": "tusizi",
    "何首乌": "heshouwu", "制何首乌": "heshouwu",
    "女贞子": "nverzhenzi", "旱莲草": "hanliancao", "墨旱莲": "hanliancao",
    "天麻": "tianma", "钩藤": "gouteng", "石决明": "shijueming",
    "地龙": "dilong", "全蝎": "quanxie", "蜈蚣": "wugong", "僵蚕": "jiangcan",
    "水蛭": "shuizhi", "三棱": "sanleng", "莪术": "ezhu",
    "雷公藤": "leigongteng", "青风藤": "qingfengteng",
    "秦艽": "qinjiao", "威灵仙": "weilingxian",
    "防风": "fangfeng", "荆芥": "jingjie", "羌活": "qianghuo", "独活": "duhuo",
    "白芷": "baizhi", "细辛": "xixin", "辛夷": "xinyi", "苍耳子": "cangerzi",
    "桑叶": "sangye", "菊花": "juhua", "薄荷": "bohe",
    "知母": "zhimu", "石膏": "shigao", "芦根": "lugeng", "白茅根": "baimaogen",
    "小蓟": "xiaoji", "大蓟": "daji", "地榆": "diyu", "槐花": "huaihua",
    "侧柏叶": "cebaiye", "蒲黄": "puhuang", "白及": "baiji",
    "吴茱萸": "wuzhuyu", "小茴香": "xiaohuixiang", "丁香": "dingxiang", "花椒": "huajiao",
    "乌梅": "wumei", "诃子": "hezi", "肉豆蔻": "roudoukou",
    "赤石脂": "chishizhi", "五味子": "wuweizi",
    "百合": "baihe", "川贝母": "chuanbeimu", "浙贝母": "zhebeimu", "瓜蒌": "gualou",
    "天冬": "tiandong", "沙参": "shashen", "玉竹": "yuzhu", "黄精": "huangjing",
}


def load_herb_cache(dirs):
    """药材 JSON 缓存 → {拼音键: 记录}，并用 cn 字段建立中文名索引"""
    cache, cnmap = {}, {}
    for d in dirs:
        for f in glob.glob(os.path.join(d, "**", "*.json"), recursive=True):
            try:
                j = json.load(open(f, encoding="utf-8"))
            except Exception:
                continue
            if not isinstance(j, dict) or not j.get("found"):
                continue
            key = os.path.basename(f).replace("herb_", "").replace("tcmsp_", "").replace(".json", "")
            if not (j.get("active") is not None or j.get("targets") is not None):
                continue
            cache[key] = j
            if j.get("cn"):
                cnmap[j["cn"]] = key
    return cache, cnmap


def resolve(herb, cache, cnmap):
    if herb in cache:
        return herb
    if herb in cnmap:
        return cnmap[herb]
    if herb in HERB_ALIAS and HERB_ALIAS[herb] in cache:
        return HERB_ALIAS[herb]
    # 去炮制前缀再试
    for pre in ("炙", "炒", "生", "制", "盐", "醋", "酒", "煅", "法", "姜", "麸炒", "蜜"):
        if herb.startswith(pre):
            k2 = cnmap.get(herb[len(pre):])
            if k2 in cache:
                return k2
    return None


def repro_one(herbs, cache, cnmap, ob, dl):
    ob = ob or 30; dl = dl or 0.18
    hit, miss = [], []
    act, tg, comp_names = set(), set(), set()
    for h in herbs:
        k = resolve(h, cache, cnmap)
        if not k:
            miss.append(h); continue
        hit.append(h)
        rec = cache[k]
        for a in (rec.get("active") or []):
            try:
                if float(a.get("ob") or 0) >= ob and float(a.get("dl") or 0) >= dl:
                    act.add(a.get("mol_id") or a.get("name")); comp_names.add(a.get("name"))
            except (TypeError, ValueError):
                continue
        for t in (rec.get("targets") or []):
            name = t.get("target_name") or t.get("gene")
            if name:
                tg.add(name.strip().lower())
    return {"herbs_total": len(herbs), "herbs_found": len(hit), "herbs_missing": miss,
            "coverage": round(len(hit) / max(1, len(herbs)), 3),
            "n_active_repro": len(act), "n_targets_repro": len(tg),
            "active_names": sorted(x for x in comp_names if x)[:12]}


def judge(report, repro, tol=0.2):
    if not report or not repro:
        return "无法比对", None
    r = repro / report
    if abs(r - 1) <= tol:
        return "一致(±20%)", r
    return ("偏低" if r < 1 else "偏高"), r



GENE_OVERRIDE = {
    "prostaglandin g/h synthase 1": "PTGS1", "prostaglandin g/h synthase 2": "PTGS2",
    "nitric oxide synthase, inducible": "NOS2", "nitric oxide synthase, endothelial": "NOS3",
    "nitric oxide synthase, brain": "NOS1", "thrombin": "F2", "progesterone receptor": "PGR",
    "androgen receptor": "AR", "estrogen receptor": "ESR1", "estrogen receptor beta": "ESR2",
    "peroxisome proliferator-activated receptor gamma": "PPARG",
    "peroxisome proliferator-activated receptor alpha": "PPARA",
    "peroxisome proliferator-activated receptor delta": "PPARD",
    "transcription factor p65": "RELA", "nf-kappa-b inhibitor alpha": "NFKBIA",
    "mitogen-activated protein kinase 1": "MAPK1", "mitogen-activated protein kinase 3": "MAPK3",
    "mitogen-activated protein kinase 14": "MAPK14", "glycogen synthase kinase-3 beta": "GSK3B",
    "phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit, gamma isoform": "PIK3CG",
    "proto-oncogene tyrosine-protein kinase src": "SRC", "caspase-3": "CASP3",
    "caspase-9": "CASP9", "bcl-2-like protein 1": "BCL2L1", "apoptosis regulator bcl-2": "BCL2",
    "tumor necrosis factor": "TNF", "interleukin-6": "IL6", "interleukin-1 beta": "IL1B",
    "vascular endothelial growth factor a": "VEGFA", "matrix metalloproteinase-9": "MMP9",
    "matrix metalloproteinase-2": "MMP2", "c-reactive protein": "CRP",
    "serine/threonine-protein kinase mtor": "MTOR", "akt1": "AKT1",
}


def build_gene_map(cache, map_path):
    """TCMSP target_name → 基因符号（本地缓存 + mygene 批量 + 人工覆盖）"""
    m = {}
    if os.path.exists(map_path):
        for line in open(map_path, encoding="utf-8-sig"):
            p = line.rstrip("\n").split("\t")
            if len(p) == 2 and p[1]:
                m[p[0].lower()] = p[1]
    names = set()
    for v in cache.values():
        for t in (v.get("targets") or []):
            g = (t.get("gene") or "").strip()
            if g:
                m[(t.get("target_name") or "").strip().lower()] = g.upper()
            else:
                nm = (t.get("target_name") or "").strip()
                if nm and nm.lower() not in m:
                    names.add(nm)
    todo = [n for n in names if n.lower() not in GENE_OVERRIDE]
    print(f"基因映射：已有 {len(m)}，待查 {len(todo)}")
    for i in range(0, len(todo), 180):
        chunk = todo[i:i + 180]
        try:
            r = requests.post("https://mygene.info/v3/query",
                              json={"q": chunk, "scopes": "name", "fields": "symbol",
                                    "species": "human", "size": 1}, timeout=120)
            for rec in r.json():
                q = rec.get("query"); sym = rec.get("symbol")
                if q and sym and not str(sym).startswith("LOC"):
                    m[q.lower()] = sym.upper()
        except Exception as e:
            print("mygene fail", i, str(e)[:60])
        time.sleep(0.3)
    for k, v in GENE_OVERRIDE.items():
        m.setdefault(k, v)
    os.makedirs(os.path.dirname(map_path), exist_ok=True)
    with open(map_path, "w", encoding="utf-8-sig") as f:
        for k, v in sorted(m.items()):
            f.write(f"{k}\t{v}\n")
    print(f"基因映射完成 {len(m)} 条 → {map_path}")
    return m


def cmd_run(args):
    m = pd.read_csv(args.methods, encoding="utf-8-sig")
    cache, cnmap = load_herb_cache(args.herb_cache.split(","))
    print(f"药材缓存: {len(cache)} 个（中文索引 {len(cnmap)}）")
    rows = []
    for _, r in m.iterrows():
        herbs = []
        try:
            v = r.get("herbs")
            herbs = eval(v) if isinstance(v, str) else (v or [])
            herbs = [str(x).strip() for x in herbs if str(x).strip()]
        except Exception:
            herbs = []
        rec = repro_one(herbs, cache, cnmap, r.get("OB"), r.get("DL"))
        rec["pmid"] = r.pmid
        for k in ("iso", "if2y", "year"):
            rec[k] = r.get(k)
        rec["report_compounds"] = r.get("n_compounds")
        rec["report_targets"] = r.get("n_targets")
        rec["report_intersection"] = r.get("n_intersection")
        rec["verdict_compounds"], rec["ratio_compounds"] = judge(
            r.get("n_compounds") if pd.notna(r.get("n_compounds")) else None, rec["n_active_repro"])
        rec["verdict_targets"], rec["ratio_targets"] = judge(
            r.get("n_targets") if pd.notna(r.get("n_targets")) else None, rec["n_targets_repro"])
        # 归因标签
        tags = []
        if herbs and rec["coverage"] < 1: tags.append("herb_missing")
        ob_v, dl_v = r.get("OB"), r.get("DL")
        if (pd.notna(ob_v) and ob_v != 30) or (pd.notna(dl_v) and dl_v != 0.18):
            tags.append("threshold_diff")
        if pd.isna(r.get("n_compounds")) and pd.isna(r.get("n_targets")): tags.append("undisclosed")
        if (rec["verdict_targets"] == "偏低"): tags.append("target_source")
        rec["attribution"] = ";".join(tags)
        rows.append(rec)
    out = pd.DataFrame(rows)
    out["herbs_missing"] = out.herbs_missing.map(lambda x: ";".join(x) if isinstance(x, list) else x)
    out["active_names"] = out.active_names.map(lambda x: ";".join(x) if isinstance(x, list) else x)
    out.to_csv(args.out, index=False, encoding="utf-8-sig")
    print(f"复现 {len(out)} 篇 → {args.out}")
    print(f"  药材覆盖率中位: {out.coverage.median():.0%} | 全药材命中: {(out.coverage>=1).sum()} 篇")
    print(f"  成分数可复现(±20%): {(out.verdict_compounds=='一致(±20%)').sum()} 篇 | 有原文值: {out.verdict_compounds.ne('无法比对').sum()}")
    print(f"  靶点数可复现(±20%): {(out.verdict_targets=='一致(±20%)').sum()} 篇 | 有原文值: {out.verdict_targets.ne('无法比对').sum()}")



def cmd_inter(args):
    """交集复现：TCMSP 靶点基因 ∪ ∩ OpenTargets 疾病靶点 → 与原文报告交集数对比"""
    m = pd.read_csv(args.methods, encoding="utf-8-sig")
    cache, cnmap = load_herb_cache(args.herb_cache.split(","))
    dis = json.load(open(args.disease_json, encoding="utf-8"))
    gmap = build_gene_map(cache, args.gene_map)
    rows = []
    for _, r in m.iterrows():
        pmid = str(r.pmid)
        d = dis.get(pmid) or {}
        if not d.get("targets"):
            continue
        herbs = []
        try:
            v = r.get("herbs")
            herbs = ast.literal_eval(v) if isinstance(v, str) else (v or [])
        except Exception:
            pass
        genes = set()
        for h in herbs:
            k = resolve(str(h).strip(), cache, cnmap)
            if not k:
                continue
            for t in (cache[k].get("targets") or []):
                g = (t.get("gene") or "").upper() or gmap.get((t.get("target_name") or "").strip().lower(), "")
                if g:
                    genes.add(g.upper())
        if not genes:
            continue
        inter = genes & set(x.upper() for x in d["targets"])
        rep = r.get("n_intersection")
        rows.append({"pmid": pmid, "disease": d.get("disease"), "n_herb_targets": len(genes),
                     "n_disease_targets": d.get("n_targets"), "n_intersection_repro": len(inter),
                     "n_intersection_report": rep if pd.notna(rep) else None,
                     "ratio": round(len(inter) / rep, 3) if (pd.notna(rep) and rep) else None,
                     "genes": ";".join(sorted(inter)[:30])})
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False, encoding="utf-8-sig")
    both = df[df.n_intersection_report.notna()]
    print(f"可算交集 {len(df)} 篇 | 有原文值 {len(both)} 篇")
    if len(both):
        print(f"  比值中位 {both.ratio.median():.2f} | ±20% 一致 {((both.ratio-1).abs()<=0.2).sum()} 篇")
        print(f"  中位：本机 {both.n_intersection_repro.median():.0f} vs 原文 {both.n_intersection_report.median():.0f}")


def cmd_stats(args):
    df = pd.read_csv(args.in_, encoding="utf-8-sig")
    os.makedirs(args.outdir, exist_ok=True)
    lines = ["# 500 篇批量复现统计", "", f"篇数: {len(df)}", "",
             "## 药材覆盖", "",
             f"- 覆盖率中位: {df.coverage.median():.0%}；完全命中: {(df.coverage>=1).sum()} 篇",
             f"- 平均每篇药味数: {df.herbs_total.mean():.1f}", ""]
    # 缺失药材 top
    from collections import Counter
    c = Counter()
    for s in df.herbs_missing.dropna():
        for h in str(s).split(";"):
            if h: c[h] += 1
    lines += ["## 未被 TCMSP 缓存覆盖的药味 Top 20", ""]
    lines += [f"- {k}: {v} 篇" for k, v in c.most_common(20)]
    lines += ["", "## 成分/靶点复现判定", ""]
    for col in ("verdict_compounds", "verdict_targets"):
        vc = df[col].value_counts().to_dict()
        lines.append(f"- {col}: {vc}")
    lines += ["", "## 归因标签分布", ""]
    tc = Counter()
    for s in df.attribution.dropna():
        for t in str(s).split(";"):
            if t: tc[t] += 1
    lines += [f"- {k}: {v} ({v/len(df):.0%})" for k, v in tc.most_common()]
    lines += ["", "## 按期刊分层（Top 10）", "", "| 期刊 | 篇数 | 药材覆盖中位 |", "|---|---|---|"]
    g = df.groupby("iso").agg(n=("pmid", "size"), cov=("coverage", "median")).sort_values("n", ascending=False).head(10)
    lines += ["| {} | {} | {:.0%} |".format(i, int(r["n"]), float(r["cov"])) for i, r in g.iterrows()]
    out = os.path.join(args.outdir, "repro_stats.md")
    open(out, "w", encoding="utf-8").write("\n".join(lines))
    print("\n".join(lines[:10])); print("→", out)


def main():
    ap = argparse.ArgumentParser(description="np_repro_batch: 批量复现与归因")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("run")
    s.add_argument("--methods", required=True); s.add_argument("--herb-cache", required=True)
    s.add_argument("--out", required=True)
    s = sub.add_parser("inter")
    s.add_argument("--methods", required=True); s.add_argument("--herb-cache", required=True)
    s.add_argument("--disease-json", required=True); s.add_argument("--out", required=True)
    s.add_argument("--gene-map", default="data/protein_symbols_extra.tsv")
    s = sub.add_parser("stats")
    s.add_argument("--in", dest="in_", required=True); s.add_argument("--outdir", default="stats")
    a = ap.parse_args()
    {"run": cmd_run, "inter": cmd_inter, "stats": cmd_stats}[a.cmd](a)


if __name__ == "__main__":
    main()
