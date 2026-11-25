#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_corpus.py —— v1.7.0 文献语料模块（500 篇级筛查 / 方法挖掘 / 归类统计）

背景：23 篇复现（v1.3）是手工小样本；本模块把"文献 → 参数表"做成可批量流水线：
  screen   元数据（PubMed JSON/CSV）→ IF 代理过滤（OpenAlex 2yr_mean_citedness）→ 选 N 篇
  harvest  PMCID → PMC OA 全文 XML → 纯文本（含表格，成分/靶点表常在 table-wrap）
  extract  全文 → 结构化方法参数（药味/成分数/靶点数/阈值/数据库/对接/MD/湿实验）
  classify 参数表 → 疾病谱/期刊/方法组合的统计与交叉表

设计原则（与 v1.3 归因口径一致）：
  - 全部数值抽取保留原文上下文片段（snippet），可人工复核，不猜；
  - 抽取失败标 None 而非默认值；"未披露"与"未使用"分开记录；
  - IF 用 OpenAlex 2yr_mean_citedness 作可复算代理（非 JCR 官方 IF，报告中如实标注）。

命令：
  python np_corpus.py screen --meta meta.json --journal-if jif.json --if-min 5 --top 500 --out sel.csv
  python np_corpus.py harvest --select sel.csv --outdir fulltext
  python np_corpus.py extract --dir fulltext --select sel.csv --out methods.csv
  python np_corpus.py classify --methods methods.csv --outdir stats
"""
import argparse
import json
import os
import re

import pandas as pd
import requests

# ---------------- 词表 ----------------



DB_TERMS = {
    "TCMSP": r"\bTCMSP\b|Traditional Chinese Medicine Systems Pharmacology",
    "BATMAN-TCM": r"BATMAN[- ]?TCM",
    "ETCM": r"\bETCM\b",
    "HERB": r"\bHERB\b(?!\s*database of traditional)",
    "SymMap": r"SymMap",
    "TCMID": r"\bTCMID\b",
    "SwissTargetPrediction": r"SwissTarget",
    "STITCH": r"\bSTITCH\b",
    "PharmMapper": r"PharmMapper",
    "GeneCards": r"GeneCards",
    "OMIM": r"\bOMIM\b",
    "DisGeNET": r"DisGeNET",
    "TTD": r"Therapeutic Target Database|\bTTD\b",
    "DrugBank": r"DrugBank",
    "STRING": r"\bSTRING\b",
    "Cytoscape": r"Cytoscape",
    "DAVID": r"\bDAVID\b",
    "Metascape": r"Metascape",
    "ClusterProfiler": r"clusterProfiler",
    "AutoDock": r"AutoDock",
    "Vina": r"\bVina\b",
    "Schrodinger": r"Schr[oö]dinger|Glide",
    "GROMACS": r"GROMACS",
    "AMBER": r"\bAMBER\b",
    "Discovery Studio": r"Discovery Studio",
    "PyMOL": r"PyMOL",
}

DISEASE_KB = {
    "肿瘤": r"carcinoma|cancer|tumor|tumour|neoplasm|adenocarcinoma|glioma|leukemia|lymphoma|melanoma|HCC|NSCLC",
    "代谢/内分泌": r"diabet|insulin resistance|obes|hyperlipid|NAFLD|fatty liver|thyroid|PCOS|osteoporosis",
    "心脑血管": r"atherosclero|myocardial|ischemi|hypertens|stroke|heart failure|arrhythm|thromb",
    "神经/精神": r"Alzheimer|Parkinson|depress|anxiety|epilep|neuropathic|insomnia|cognitive|stroke",
    "消化/肝": r"hepat|ulcerative colitis|gastr|intestinal|constipation|diarrhea|pancreat",
    "肾/泌尿": r"nephro|kidney|renal|diabetic kidney|prostat|bladder",
    "呼吸/免疫": r"asthma|COPD|pneumonia|lung injur|fibrosis|psoriasis|arthritis|rhinitis|allerg",
    "妇科/骨": r"endometrio|ovarian|uterine|menopaus|osteopor|bone",
    "感染/炎症": r"sepsis|infection|inflammat|COVID|periodont|wound",
}

HERBS = ("黄芪|当归|丹参|川芎|赤芍|白芍|地黄|熟地黄|山药|山茱萸|茯苓|泽泻|牡丹皮|桂枝|麻黄|杏仁|甘草|"
         "人参|党参|西洋参|太子参|白术|苍术|陈皮|半夏|厚朴|枳实|枳壳|柴胡|黄芩|黄连|黄柏|大黄|芒硝|"
         "金银花|连翘|板蓝根|蒲公英|鱼腥草|白头翁|黄连|栀子|知母|石膏|葛根|升麻|羌活|独活|防风|荆芥|"
         "薄荷|桑叶|菊花|枸杞|女贞子|旱莲草|何首乌|淫羊藿|巴戟天|杜仲|续断|牛膝|桃仁|红花|益母草|"
         "三棱|莪术|水蛭|地龙|全蝎|蜈蚣|僵蚕|天麻|钩藤|石决明|牡蛎|龙骨|酸枣仁|远志|石菖蒲|"
         "五味子|麦冬|天冬|沙参|玉竹|百合|川贝母|浙贝母|瓜蒌|桔梗|紫苏|生姜|干姜|肉桂|附子|吴茱萸|"
         "小茴香|丁香|花椒|茵陈|金钱草|虎杖|垂盆草|五味子|乌梅|诃子|肉豆蔻|赤石脂|"
         "雷公藤|青风藤|秦艽|威灵仙|桑寄生|五加皮|狗脊|骨碎补|补骨脂|菟丝子|覆盆子|"
         "麻黄|细辛|白芷|藁本|辛夷|苍耳子|芦根|白茅根|小蓟|大蓟|地榆|槐花|侧柏叶|三七|蒲黄|"
         "延胡索|郁金|姜黄|乳香|没药|五灵脂|鸡血藤|牛膝|穿山甲|王不留行|路路通|"
         "陈皮|青皮|香附|木香|砂仁|白豆蔻|草果|藿香|佩兰|苍术|厚朴|茯苓|猪苓|薏苡仁|车前子|"
         "滑石|通草|瞿麦|萹蓄|海金沙|石韦|冬葵子|灯心草|萆薢|土茯苓|")
HERB_EN = {
    "Astragalus": "黄芪", "Radix Astragali": "黄芪", "Angelica sinensis": "当归", "Angelica":
    "当归", "Salvia miltiorrhiza": "丹参", "Danshen": "丹参", "Ligusticum chuanxiong": "川芎",
    "Chuanxiong": "川芎", "Paeonia lactiflora": "白芍", "Poria cocos": "茯苓", "Atractylodes": "白术",
    "Glycyrrhiza": "甘草", "Licorice": "甘草", "Panax ginseng": "人参", "Ginseng": "人参",
    "Coptis chinensis": "黄连", "Coptis": "黄连", "Scutellaria baicalensis": "黄芩", "Scutellaria": "黄芩",
    "Rheum": "大黄", "Rhubarb": "大黄", "Lonicera": "金银花", "Honeysuckle": "金银花",
    "Forsythia": "连翘", "Pueraria": "葛根", "Kudzu": "葛根", "Rehmannia": "地黄",
    "Cornus officinalis": "山茱萸", "Dioscorea": "山药", "Alisma": "泽泻", "Moutan": "牡丹皮",
    "Cinnamomum cassia": "肉桂", "Aconitum": "附子", "Zingiber": "干姜", "Pinellia": "半夏",
    "Citrus reticulata": "陈皮", "Magnolia officinalis": "厚朴", "Aurantium": "枳壳",
    "Bupleurum": "柴胡", "Gardenia": "栀子", "Anemarrhena": "知母", "Gypsum": "石膏",
    "Epimedium": "淫羊藿", "Eucommia": "杜仲", "Achyranthes": "牛膝", "Prunus persica": "桃仁",
    "Carthamus": "红花", "Leonurus": "益母草", "Tripterygium": "雷公藤", "Celastrus": "雷公藤",
    "Curcuma": "姜黄", "Corydalis": "延胡索", "Spatholobus": "鸡血藤", "Polygonum": "何首乌",
    "Schisandra": "五味子", "Ophiopogon": "麦冬", "Adenophora": "沙参", "Lilium": "百合",
    "Trichosanthes": "瓜蒌", "Platycodon": "桔梗", "Perilla": "紫苏", "Ziziphus": "酸枣仁",
    "Polygala": "远志", "Acorus": "石菖蒲", "Artemisia": "茵陈", "Lysimachia": "金钱草",
    "Polygonum cuspidatum": "虎杖", "Rhodiola": "红景天", "Ginkgo": "银杏", "Camellia": "茶",
    "Lycium": "枸杞", "Morinda": "巴戟天", "Psoralea": "补骨脂", "Cuscuta": "菟丝子",
    "Ephedra": "麻黄", "Asarum": "细辛", "Angelica dahurica": "白芷", "Notopterygium": "羌活",
    "Saposhnikovia": "防风", "Schizonepeta": "荆芥", "Mentha": "薄荷", "Morus": "桑叶",
    "Chrysanthemum": "菊花", "Panax notoginseng": "三七", "Notoginseng": "三七",
    "Uncaria": "钩藤", "Gastrodia": "天麻", "Alisma orientale": "泽泻", "Coix": "薏苡仁",
    "Plantago": "车前子", "Dioscorea nipponica": "穿山龙", "Smilax": "土茯苓",
}

HERB_EN.update({
 "Rheum palmatum": "大黄", "Rhei Radix": "大黄", "Polygonum multiflorum": "何首乌",
 "Scutellariae Radix": "黄芩", "Coptidis Rhizoma": "黄连", "Phellodendron": "黄柏",
 "Gardeniae Fructus": "栀子", "Bupleuri Radix": "柴胡", "Angelicae Sinensis Radix": "当归",
 "Paeoniae Radix Alba": "白芍", "Atractylodis Macrocephalae Rhizoma": "白术",
 "Poria cocos": "茯苓", "Alismatis Rhizoma": "泽泻", "Moutan Cortex": "牡丹皮",
 "Cinnamomi Ramulus": "桂枝", "Ephedrae Herba": "麻黄", "Armeniacae Semen Amarum": "杏仁",
 "Glycyrrhizae Radix": "甘草", "Ginseng Radix": "人参", "Codonopsis": "党参",
 "Astragali Radix": "黄芪", "Rehmanniae Radix": "地黄", "Dioscoreae Rhizoma": "山药",
 "Corni Fructus": "山茱萸", "Chuanxiong Rhizoma": "川芎", "Salviae Miltiorrhizae Radix": "丹参",
 "Paeoniae Radix Rubra": "赤芍", "Persicae Semen": "桃仁", "Carthami Flos": "红花",
 "Leonuri Herba": "益母草", "Curcumae Longae Rhizoma": "姜黄", "Curcumae Radix": "郁金",
 "Corydalis Rhizoma": "延胡索", "Spatholobi Caulis": "鸡血藤", "Achyranthis Bidentatae Radix": "牛膝",
 "Eucommiae Cortex": "杜仲", "Dipsaci Radix": "续断", "Taxilli Herba": "桑寄生",
 "Epimedii Folium": "淫羊藿", "Psoraleae Fructus": "补骨脂", "Cuscutae Semen": "菟丝子",
 "Polygoni Multiflori Radix": "何首乌", "Ligustri Lucidi Fructus": "女贞子",
 "Ecliptae Herba": "墨旱莲", "Gastrodiae Rhizoma": "天麻", "Uncariae Ramulus cum Uncis": "钩藤",
 "Haliotidis Concha": "石决明", "Pheretima": "地龙", "Scorpio": "全蝎", "Scolopendra": "蜈蚣",
 "Bombyx Batryticatus": "僵蚕", "Hirudo": "水蛭", "Sparganii Rhizoma": "三棱",
 "Curcumae Rhizoma": "莪术", "Tripterygii Wilfordii Radix": "雷公藤", "Sinomenii Caulis": "青风藤",
 "Gentianae Macrophyllae Radix": "秦艽", "Clematidis Radix": "威灵仙", "Saposhnikoviae Radix": "防风",
 "Schizonepetae Herba": "荆芥", "Notopterygii Rhizoma": "羌活", "Angelicae Pubescentis Radix": "独活",
 "Angelicae Dahuricae Radix": "白芷", "Asari Radix et Rhizoma": "细辛", "Magnoliae Flos": "辛夷",
 "Xanthii Fructus": "苍耳子", "Mori Folium": "桑叶", "Chrysanthemi Flos": "菊花",
 "Anemarrhenae Rhizoma": "知母", "Gypsum Fibrosum": "石膏", "Phragmitis Rhizoma": "芦根",
 "Imperatae Rhizoma": "白茅根", "Cirsii Herba": "小蓟", "Cirsii Japonici Herba": "大蓟",
 "Sanguisorbae Radix": "地榆", "Sophorae Flos": "槐花", "Platycladi Cacumen": "侧柏叶",
 "Typhae Pollen": "蒲黄", "Bletillae Rhizoma": "白及", "Evodiae Fructus": "吴茱萸",
 "Foeniculi Fructus": "小茴香", "Caryophylli Flos": "丁香", "Zanthoxyli Pericarpium": "花椒",
 "Mume Fructus": "乌梅", "Chebulae Fructus": "诃子", "Myristicae Semen": "肉豆蔻",
 "Halloysitum Rubrum": "赤石脂", "Lilii Bulbus": "百合", "Fritillariae Cirrhosae Bulbus": "川贝母",
 "Fritillariae Thunbergii Bulbus": "浙贝母", "Trichosanthis Fructus": "瓜蒌",
 "Asparagi Radix": "天冬", "Adenophorae Radix": "沙参", "Polygonati Odorati Rhizoma": "玉竹",
 "Polygonati Rhizoma": "黄精", "Moutan": "牡丹皮", "Poria": "茯苓", "Alisma": "泽泻",
 "Schisandra": "五味子", "Ophiopogon": "麦冬", "Ziziphus": "酸枣仁", "Polygala": "远志",
 "Acorus": "石菖蒲", "Lycium": "枸杞", "Ginkgo": "银杏", "Rhodiola": "红景天",
 "Panax notoginseng": "三七", "Salvia": "丹参", "Astragalus membranaceus": "黄芪",
 "Angelica sinensis": "当归", "Rehmannia glutinosa": "地黄", "Paeonia lactiflora": "白芍",
 "Scutellaria baicalensis": "黄芩", "Coptis chinensis": "黄连", "Rheum officinale": "大黄",
})


HERB_EN.update({  # 拼音名（英文论文常见写法）
 "Huangqi": "黄芪", "Huang Qi": "黄芪", "Danggui": "当归", "Dang Gui": "当归",
 "Gancao": "甘草", "Gan Cao": "甘草", "Dahuang": "大黄", "Da Huang": "大黄",
 "Renshen": "人参", "Ren Shen": "人参", "Dangshen": "党参", "Dang Shen": "党参",
 "Baizhu": "白术", "Bai Zhu": "白术", "Fuling": "茯苓", "Fu Ling": "茯苓",
 "Danzhen": "丹参", "Danshen": "丹参", "Dan Shen": "丹参", "Chuanxiong": "川芎",
 "Baishao": "白芍", "Bai Shao": "白芍", "Chishao": "赤芍", "Huangqin": "黄芩",
 "Huang Qin": "黄芩", "Huanglian": "黄连", "Huang Lian": "黄连", "Huangbai": "黄柏",
 "Chaihu": "柴胡", "Chai Hu": "柴胡", "Zhizi": "栀子", "Zhi Zi": "栀子",
 "Jinyinhua": "金银花", "Jin Yin Hua": "金银花", "Lianqiao": "连翘", "Lian Qiao": "连翘",
 "Gegen": "葛根", "Ge Gen": "葛根", "Chenpi": "陈皮", "Chen Pi": "陈皮",
 "Banxia": "半夏", "Ban Xia": "半夏", "Houpo": "厚朴", "Zhishi": "枳实",
 "Zhiqiao": "枳壳", "Jiegeng": "桔梗", "Jie Geng": "桔梗", "Xingren": "杏仁",
 "Mahuang": "麻黄", "Ma Huang": "麻黄", "Guizhi": "桂枝", "Gui Zhi": "桂枝",
 "Rougui": "肉桂", "Fuzi": "附子", "Ganjiang": "干姜", "Gan Jiang": "干姜",
 "Shengjiang": "生姜", "Zexie": "泽泻", "Zhe Xie": "泽泻", "Mudanpi": "牡丹皮",
 "Shanzhuyu": "山茱萸", "Shanyao": "山药", "Shan Yao": "山药", "Shudihuang": "地黄",
 "Dihuang": "地黄", "Di Huang": "地黄", "Yinchen": "茵陈", "Yin Chen": "茵陈",
 "Jianghuang": "姜黄", "Jiang Huang": "姜黄", "Yujin": "郁金", "Yanhusuo": "延胡索",
 "Taoren": "桃仁", "Honghua": "红花", "Hong Hua": "红花", "Yimucao": "益母草",
 "Niuxi": "牛膝", "Niu Xi": "牛膝", "Duzhong": "杜仲", "Xuduan": "续断",
 "Yinyanghuo": "淫羊藿", "Buguzhi": "补骨脂", "Tusizi": "菟丝子", "Heshouwu": "何首乌",
 "Nvzhenzi": "女贞子", "Hanliancao": "墨旱莲", "Tianma": "天麻", "Tian Ma": "天麻",
 "Gouteng": "钩藤", "Dilong": "地龙", "Quanxie": "全蝎", "Wugong": "蜈蚣",
 "Jiangcan": "僵蚕", "Shuizhi": "水蛭", "Sanleng": "三棱", "Ezhu": "莪术",
 "Leigongteng": "雷公藤", "Lei Gong Teng": "雷公藤", "Qinjiao": "秦艽",
 "Weilingxian": "威灵仙", "Fangfeng": "防风", "Fang Feng": "防风", "Jingjie": "荆芥",
 "Qianghuo": "羌活", "Duhuo": "独活", "Baizhi": "白芷", "Baizhi": "白芷",
 "Xixin": "细辛", "Xinyi": "辛夷", "Cangerzi": "苍耳子", "Sangye": "桑叶",
 "Juhua": "菊花", "Zhimu": "知母", "Shigao": "石膏", "Lugen": "芦根",
 "Baimaogen": "白茅根", "Xiaoji": "小蓟", "Daji": "大蓟", "Diyu": "地榆",
 "Huaihua": "槐花", "Cebaiye": "侧柏叶", "Puhuang": "蒲黄", "Wuzhuyu": "吴茱萸",
 "Xiaohuixiang": "小茴香", "Dingxiang": "丁香", "Huajiao": "花椒", "Wumei": "乌梅",
 "Hezi": "诃子", "Roudoukou": "肉豆蔻", "Chishizhi": "赤石脂", "Baihe": "百合",
 "Chuanbeimu": "川贝母", "Zhebeimu": "浙贝母", "Gualou": "瓜蒌", "Tiandong": "天冬",
 "Shashen": "沙参", "Yuzhu": "玉竹", "Huangjing": "黄精", "Wuweizi": "五味子",
 "Maidong": "麦冬", "Suanzaoren": "酸枣仁", "Yuanzhi": "远志", "Shichangpu": "石菖蒲",
 "Gouqizi": "枸杞子", "Sanqi": "三七", "San Qi": "三七", "Tufuling": "土茯苓",
 "Bixie": "萆薢", "Huzhang": "虎杖", "Jinqiancao": "金钱草", "Xiangfu": "香附",
 "Muxiang": "木香", "Sharen": "砂仁", "Huoxiang": "藿香", "Peilan": "佩兰",
 "Zhuling": "猪苓", "Huashi": "滑石", "Muli": "牡蛎", "Longgu": "龙骨",
 "Sangjisheng": "桑寄生", "Nverzhenzi": "女贞子", "Yiyiren": "薏苡仁",
 "Cheqianzi": "车前子", "Lugeng": "芦根", "Baiji": "白及",
})

# ---------------- 抽取规则 ----------------

def _snip(txt, m, pad=90):
    s = max(0, m.start() - pad); e = min(len(txt), m.end() + pad)
    return re.sub(r"\s+", " ", txt[s:e])


def _num(txt, pats):
    for p in pats:
        m = re.search(p, txt, re.I)
        if m:
            try:
                return int(m.group(1)), _snip(txt, m)
            except (ValueError, IndexError):
                continue
    return None, ""


def extract_one(txt):
    low = txt
    rec = {}
    # 阈值
    m = re.search(r"OB\)?\s*(?:≥|>=|>|of|was set to|threshold of)\s*(\d{2})\s*%", low, re.I)
    rec["OB"] = int(m.group(1)) if m else None
    m = re.search(r"DL\)?\s*(?:≥|>=|>|of|was set to|threshold of)\s*(0?\.\d+)", low, re.I)
    rec["DL"] = float(m.group(1)) if m else None
    # 数量
    rec["n_compounds"], rec["s_compounds"] = _num(low, [
        r"(\d{1,4})\s+(?:active\s+)?(?:compounds?|ingredients?|chemical constituents?)\s+(?:were|was)\s+(?:obtained|identified|screened|retained)",
        r"(\d{1,4})\s+(?:active\s+)?compounds?\s+(?:including|such as|in total)",
        r"a\s+total\s+of\s+(\d{1,4})\s+(?:active\s+)?compounds?",
        r"(\d{1,4})\s+(?:active\s+)?ingredients?\s+(?:were|was)",
    ])
    rec["n_targets"], rec["s_targets"] = _num(low, [
        r"(\d{1,5})\s+(?:potential\s+|predicted\s+|putative\s+)?targets?\s+(?:were|was)\s+(?:obtained|identified|predicted|collected)",
        r"a\s+total\s+of\s+(\d{1,5})\s+(?:potential\s+)?targets?",
    ])
    rec["n_intersection"], rec["s_intersection"] = _num(low, [
        r"(\d{1,5})\s+(?:common\s+|overlapping\s+|intersection\s+|shared\s+)?(?:targets?|genes?)\s+(?:were|was)\s+(?:obtained|identified|found|screened)",
        r"(\d{1,5})\s+(?:common|overlapping|shared)\s+(?:targets?|genes?)",
        r"intersection\s+(?:of\s+)?(\d{1,5})",
        r"(\d{1,5})\s+intersection\s+(?:targets?|genes?)",
    ])
    rec["n_ppi_nodes"], rec["s_ppi_nodes"] = _num(low, [
        r"(\d{1,5})\s+(?:nodes?|proteins?)\s+(?:and|,)\s+\d{1,6}\s+edges",
    ])
    rec["n_ppi_edges"], rec["s_ppi_edges"] = _num(low, [
        r"\d{1,5}\s+(?:nodes?|proteins?)\s+(?:and|,)\s+(\d{1,6})\s+edges",
    ])
    # 对接能量（取最负）
    norm = low.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-")
    raw = re.findall(r"(-?\d{1,2}(?:\.\d{1,2})?)\s*(?:kcal\s*[/\u00b7\s]?\s*mol|kcal/mol|kJ/mol)", norm, re.I)
    es = []
    for x in raw:
        try:
            v = float(x)
        except ValueError:
            continue
        if -15.0 <= v <= -1.0:
            es.append(v)
    es = sorted(es)
    rec["dock_min_kcal"] = min(es) if es else None
    rec["n_dock_pairs"] = len(re.findall(r"binding energy", low, re.I)) or None
    # 数据库
    rec["databases"] = [k for k, p in DB_TERMS.items() if re.search(p, low, re.I)]
    # 药味（中文）
    zh = sorted(set(re.findall(HERBS, low)))
    en = sorted({v for k, v in HERB_EN.items() if re.search(re.escape(k), low, re.I)})
    rec["herbs"] = sorted({x for x in (set(zh) | set(en)) if x and x.strip()})
    rec["n_herbs"] = len(rec["herbs"])
    # 疾病
    title = str(rec.get("title") or "")
    if title:
        rec["diseases"] = [k for k, p in DISEASE_KB.items() if re.search(p, title, re.I)]
        m = DISEASE_TITLE.search(title)
        rec["primary_disease"] = m.group(1).strip(" ,.-") if m else None
    else:
        rec["diseases"] = [k for k, p in DISEASE_KB.items() if re.search(p, low, re.I)]
    # 证型
    rec["has_MD"] = bool(re.search(r"molecular dynamics|GROMACS|AMBER|RMSD", low, re.I))
    rec["has_wetlab"] = bool(re.search(r"western blot|qPCR|RT-qPCR|immunohisto|HE staining|ELISA", low, re.I))
    rec["has_animal"] = bool(re.search(r"\bmice\b|\brats?\b|zebrafish|C57BL", low, re.I))
    rec["has_network_toxicology"] = bool(re.search(r"network toxicolog", low, re.I))
    rec["has_survival"] = bool(re.search(r"TCGA|Kaplan|GEPIA|prognos", low, re.I))
    rec["has_GEO"] = bool(re.search(r"\bGSE\d{4,6}\b", low))
    return rec


def cmd_harvest(args):
    """PMC OA 全文 XML → 纯文本（标题/摘要/正文/表格）"""
    import time as _t
    import xml.etree.ElementTree as ET
    os.makedirs(args.outdir, exist_ok=True)
    sel = pd.read_csv(args.select, encoding="utf-8-sig")
    B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    def to_text(xml_bytes):
        """XML 全文：保留标题/摘要/正文段落/表格，段落与单元格换行，内联标签补空格"""
        import html as _html
        s = xml_bytes.decode("utf-8", "ignore")
        s = re.sub(r"<ref-list.*?</ref-list>", " ", s, flags=re.S)
        s = re.sub(r"<(?:p|title|td|th|caption|label|sec|abstract|article-title)[^>]*>", "\n", s)
        s = re.sub(r"<[^>]+>", " ", s)
        s = _html.unescape(s)
        s = re.sub(r"[ \t]+", " ", s)
        s = re.sub(r"\n\s*\n+", "\n", s)
        return s

    n_ok = 0
    for _, r in sel.iterrows():
        pmid, pmcid = str(r.pmid), str(r.pmcid)
        fn = os.path.join(args.outdir, f"{pmid}.txt")
        if os.path.exists(fn) and os.path.getsize(fn) > 1500:
            n_ok += 1; continue
        for k in range(3):
            try:
                xml = requests.post(B + "efetch.fcgi",
                                    data={"db": "pmc", "id": pmcid, "retmode": "xml"}, timeout=120).content
                open(fn, "w", encoding="utf-8").write(to_text(xml))
                open(fn.replace(".txt", ".xml"), "wb").write(xml)
                n_ok += 1
                break
            except Exception:
                _t.sleep(2 * (k + 1))
        _t.sleep(0.35)
    print(f"全文落地 {n_ok}/{len(sel)} → {args.outdir}")


def cmd_extract(args):
    rows = []
    sel = pd.read_csv(args.select, encoding="utf-8-sig") if args.select else None
    meta = {str(r.pmid): r for _, r in sel.iterrows()} if sel is not None else {}
    files = [f for f in os.listdir(args.dir) if f.endswith(".txt")]
    for fn in files:
        pmid = fn.replace(".txt", "")
        txt = open(os.path.join(args.dir, fn), encoding="utf-8", errors="ignore").read()
        rec = extract_one(txt)
        rec["pmid"] = pmid
        rec["chars"] = len(txt)
        if pmid in meta:
            m = meta[pmid]
            for k in ("title", "journal", "iso", "year", "if2y", "pmcid", "china_affil"):
                rec[k] = m.get(k)
        rows.append(rec)
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, index=False, encoding="utf-8-sig")
    print(f"抽取 {len(df)} 篇 → {args.out}")
    for k in ("OB", "DL", "n_compounds", "n_targets", "n_intersection", "dock_min_kcal"):
        print(f"  {k}: 有值 {df[k].notna().sum()}/{len(df)}")
    print(f"  药味抽取成功: {(df.n_herbs>0).sum()}/{len(df)} | MD {df.has_MD.sum()} | 湿实验 {df.has_wetlab.sum()}")


def cmd_classify(args):
    df = pd.read_csv(args.methods, encoding="utf-8-sig")
    os.makedirs(args.outdir, exist_ok=True)
    md = ["# 500 篇语料归类统计", "", f"篇数: {len(df)}", ""]
    md += ["## 期刊分布（Top 15）", "", "| 期刊 | 篇数 | 中位 IF 代理 |", "|---|---|---|"]
    g = df.groupby("iso").agg(n=("pmid", "size"), iff=("if2y", "median")).sort_values("n", ascending=False).head(15)
    md += [f"| {i} | {int(r.n)} | {r.iff:.1f} |" for i, r in g.iterrows()]
    md += ["", "## 疾病谱（多标签）", ""]
    from collections import Counter
    c = Counter()
    if "primary_disease" in df.columns:
        for d in df.primary_disease.dropna():
            c[str(d).lower()] += 1
    md += [f"- {k}: {v}" for k, v in c.most_common(25)]
    md += ["", "## 方法特征", "",
           f"- 有湿实验: {int(df.has_wetlab.sum())} ({df.has_wetlab.mean():.0%})",
           f"- 有动物实验: {int(df.has_animal.sum())} ({df.has_animal.mean():.0%})",
           f"- 有 MD: {int(df.has_MD.sum())} ({df.has_MD.mean():.0%})",
           f"- 有 TCGA/生存: {int(df.has_survival.sum())} ({df.has_survival.mean():.0%})",
           f"- 有 GEO 验证: {int(df.has_GEO.sum())} ({df.has_GEO.mean():.0%})",
           f"- 网络毒理学: {int(df.has_network_toxicology.sum())}",
           "", "## 阈值披露", "",
           f"- OB 披露: {int(df.OB.notna().sum())} ({df.OB.notna().mean():.0%})；分布 {df.OB.value_counts().head(5).to_dict()}",
           f"- DL 披露: {int(df.DL.notna().sum())} ({df.DL.notna().mean():.0%})；分布 {df.DL.value_counts().head(5).to_dict()}",
           "", "## 数据库使用频次", ""]
    dc = Counter()
    for s in df.databases.dropna():
        try:
            for d in eval(s) if isinstance(s, str) else s:
                dc[d] += 1
        except Exception:
            pass
    md += [f"- {k}: {v} ({v/len(df):.0%})" for k, v in dc.most_common(15)]
    # 关键：数量关系
    md += ["", "## 报告数量 vs 阈值（可复现性检查）", ""]
    sub = df[df.n_intersection.notna()]
    md += [f"- 报告交集基因数: {len(sub)} 篇；中位 {sub.n_intersection.median():.0f}（范围 {sub.n_intersection.min():.0f}–{sub.n_intersection.max():.0f}）"]
    md += [f"- 报告成分数: {df.n_compounds.notna().sum()} 篇；中位 {df.n_compounds.median() if df.n_compounds.notna().any() else '—'}"]
    md += [f"- 报告靶点数: {df.n_targets.notna().sum()} 篇；中位 {df.n_targets.median() if df.n_targets.notna().any() else '—'}"]
    open(os.path.join(args.outdir, "corpus_stats.md"), "w", encoding="utf-8").write("\n".join(md))
    print("\n".join(md[:8]))
    print(f"→ {args.outdir}/corpus_stats.md")



OT_API = "https://api.platform.opentargets.org/api/v4/graphql"


def _ot_search(name):
    q = ('{ search(queryString: "%s", entityNames: ["disease"], page: {index: 0, size: 3}) '
         '{ hits { id name } } }') % name.replace('"', " ")
    r = requests.post(OT_API, json={"query": q}, timeout=60).json()
    return (r.get("data", {}).get("search", {}) or {}).get("hits", [])


def _ot_targets(efo_id, size=1500):
    q = ('{ disease(efoId: "%s") { associatedTargets(page: {index: 0, size: %d}) '
         '{ rows { target { approvedSymbol } score } } } }') % (efo_id, size)
    r = requests.post(OT_API, json={"query": q}, timeout=120).json()
    rows = (((r.get("data") or {}).get("disease") or {}).get("associatedTargets") or {}).get("rows", [])
    return {x["target"]["approvedSymbol"].upper(): float(x.get("score") or 0) for x in rows if x.get("target")}


DISEASE_TITLE = re.compile(
    r"((?:[A-Za-z-]+\s){0,4}(?:cancer|carcinoma|adenocarcinoma|tumor|tumour|sarcoma|glioma|leukemia|lymphoma|"
    r"melanoma|myeloma|diabetes|nephropathy|kidney injury|colitis|hepatitis|liver injury|fibrosis|"
    r"atherosclerosis|myocardial infarction|heart failure|hypertension|stroke|ischemia|"
    r"Alzheimer'?s disease|Parkinson'?s disease|depression|anxiety|epilepsy|insomnia|"
    r"osteoporosis|osteoarthritis|rheumatoid arthritis|psoriasis|asthma|COPD|pneumonia|sepsis|"
    r"endometriosis|polycystic ovary syndrome|obesity|hyperlipidemia|fatty liver|"
    r"ulcerative colitis|irritable bowel syndrome|gastric ulcer|constipation|"
    r"cerebral ischemia|spinal cord injury|disc degeneration|retinopathy|periodontitis))",
    re.I)


def cmd_disease(args):
    import time as _t
    df = pd.read_csv(args.methods, encoding="utf-8-sig")
    out = {}
    cache_path = args.cache if args.cache else None
    if cache_path and os.path.exists(cache_path):
        out = json.load(open(cache_path, encoding="utf-8"))
    n = 0
    for _, r in df.iterrows():
        pmid = str(r.pmid)
        if pmid in out:
            continue
        title = str(r.get("title") or "")
        m = DISEASE_TITLE.search(title)
        if not m:
            out[pmid] = {"disease": None, "efo": None, "n_targets": 0}; continue
        name = m.group(1).strip(" ,.-")
        try:
            hits = _ot_search(name)
            if not hits:
                out[pmid] = {"disease": name, "efo": None, "n_targets": 0}; continue
            efo = hits[0]["id"]
            tg = _ot_targets(efo)
            out[pmid] = {"disease": hits[0]["name"], "efo": efo, "n_targets": len(tg),
                         "targets": sorted(tg)[:1200]}
            n += 1
        except Exception as e:
            out[pmid] = {"disease": name, "efo": None, "n_targets": 0, "error": str(e)[:80]}
        _t.sleep(0.2)
    if cache_path:
        json.dump(out, open(cache_path, "w", encoding="utf-8"), ensure_ascii=False)
    print(f"疾病映射完成 {sum(1 for v in out.values() if v.get('efo'))}/{len(out)} 篇有 EFO 靶点集")


def cmd_screen(args):
    meta = json.load(open(args.meta, encoding="utf-8")) if args.meta.endswith(".json") else None
    jif = json.load(open(args.journal_if, encoding="utf-8"))
    df = pd.DataFrame(meta) if meta else pd.read_csv(args.meta, encoding="utf-8-sig")
    df["if2y"] = df.issn.map(lambda i: (jif.get(i) or {}).get("if2y"))
    df["has_pmc"] = df.pmcid.astype(str).str.startswith("PMC")
    sel = df[(df.if2y >= args.if_min)]
    if args.drop_reviews:
        sel = sel[~sel.title.str.contains(r"\breview\b|overview", case=False, na=False)]
    sel = sel.sort_values(["has_pmc", "if2y"], ascending=False).head(args.top)
    sel.to_csv(args.out, index=False, encoding="utf-8-sig")
    print(f"筛选 {len(sel)} 篇（IF代理≥{args.if_min}，OA {sel.has_pmc.sum()}）→ {args.out}")


def main():
    ap = argparse.ArgumentParser(description="np_corpus: 文献语料筛查/方法挖掘/统计")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("screen")
    s.add_argument("--meta", required=True); s.add_argument("--journal-if", required=True)
    s.add_argument("--if-min", type=float, default=5); s.add_argument("--top", type=int, default=500)
    s.add_argument("--drop-reviews", action="store_true"); s.add_argument("--out", required=True)
    s = sub.add_parser("harvest")
    s.add_argument("--select", required=True); s.add_argument("--outdir", default="fulltext")
    s = sub.add_parser("extract")
    s.add_argument("--dir", required=True); s.add_argument("--select", default=None)
    s.add_argument("--out", required=True)
    s = sub.add_parser("disease")
    s.add_argument("--methods", required=True); s.add_argument("--cache", default="disease_targets.json")
    s = sub.add_parser("classify")
    s.add_argument("--methods", required=True); s.add_argument("--outdir", default="stats")
    a = ap.parse_args()
    {"screen": cmd_screen, "harvest": cmd_harvest, "extract": cmd_extract, "disease": cmd_disease, "classify": cmd_classify}[a.cmd](a)


if __name__ == "__main__":
    main()
