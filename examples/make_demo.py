#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_demo.py — 生成演示数据（**全部为构造数据，仅用于跑通流程，不得用于发表**）

数据设计目标（刻意包含"坑"，用来演示各质控模块的作用）：
  * 成分表里混入低 OB / 低 DL 成分 → 演示筛选与剔除理由
  * 成分表里混入 3 个"泛成分"（槲皮素/山柰酚/β-谷甾醇）→ 演示泛成分标记
  * 靶点表用 TCMSP 风格蛋白全名 → 演示标准化四层匹配
  * 靶点表混入 2 条非人类条目 → 演示非人类剔除
  * 疾病靶点表给出 score 与来源 → 演示证据强度分层
  * PPI 边表含自环与重复边 → 演示清洗
  * 富集 GMT 里含 "Pathways in cancer" 这类伞状通路 → 演示三层分级
  * 对接表含 1 个重对接 RMSD=2.6 Å 的失败案例 → 演示门禁

生成方式固定随机种子，结果可复现。
"""
import csv
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(20261119)

HERBS = ["黄芪", "当归", "丹参"]

# 成分：名称, OB, DL, HL, MW, 含量%, 主要来源药材
COMPOUNDS = [
    ("Quercetin", 46.43, 0.28, 14.40, 302.24, 0.042, "黄芪|当归|丹参"),
    ("Kaempferol", 41.88, 0.24, 15.20, 286.24, 0.031, "黄芪|当归"),
    ("Luteolin", 36.16, 0.25, 15.94, 286.24, 0.018, "丹参"),
    ("beta-sitosterol", 36.91, 0.75, 5.32, 414.71, 0.120, "黄芪|当归|丹参"),
    ("Stigmasterol", 43.83, 0.76, 5.62, 412.69, 0.055, "当归"),
    ("Formononetin", 69.67, 0.22, 16.10, 268.26, 0.026, "黄芪"),
    ("Calycosin", 47.75, 0.24, 15.30, 284.26, 0.038, "黄芪"),
    ("Astragaloside IV", 22.31, 0.15, 3.20, 784.97, 0.085, "黄芪"),
    ("Astragaloside I", 18.20, 0.12, 2.90, 869.04, 0.041, "黄芪"),
    ("Fermononetin", 51.20, 0.21, 15.80, 268.26, 0.009, "黄芪"),
    ("Jaranol", 50.83, 0.29, 16.90, 314.29, 0.004, "黄芪"),
    ("Hederagenin", 36.91, 0.75, 5.10, 472.70, 0.007, "黄芪"),
    ("Isorhamnetin", 49.60, 0.31, 14.20, 316.26, 0.012, "黄芪"),
    ("3,9-di-O-methylnissolin", 53.74, 0.48, 9.60, 314.33, 0.005, "黄芪"),
    ("Ferulic acid", 45.32, 0.06, 2.20, 194.19, 0.210, "当归"),
    ("Z-ligustilide", 51.30, 0.07, 1.80, 190.24, 0.480, "当归"),
    ("Senkyunolide A", 62.10, 0.09, 1.50, 192.25, 0.130, "当归"),
    ("n-butylidenephthalide", 58.40, 0.08, 1.20, 188.22, 0.095, "当归"),
    ("Vanillic acid", 35.48, 0.05, 1.10, 168.15, 0.022, "当归"),
    ("Tanshinone IIA", 49.48, 0.29, 8.60, 294.34, 0.180, "丹参"),
    ("Cryptotanshinone", 52.34, 0.40, 7.90, 296.36, 0.120, "丹参"),
    ("Salvianolic acid B", 3.01, 0.05, 1.40, 718.62, 0.650, "丹参"),
    ("Dihydrotanshinone I", 44.10, 0.33, 6.80, 278.26, 0.045, "丹参"),
    ("Rosmarinic acid", 21.50, 0.04, 1.60, 360.31, 0.088, "丹参"),
    ("Lithospermic acid", 9.20, 0.06, 1.30, 538.46, 0.030, "丹参"),
    ("Ursolic acid", 17.60, 0.68, 6.20, 456.70, 0.011, "当归"),
    ("Gallic acid", 31.69, 0.04, 1.05, 170.12, 0.015, "丹参"),
    ("Caffeic acid", 25.30, 0.05, 1.15, 180.16, 0.019, "丹参"),
]

# 蛋白靶点（TCMSP 风格，刻意用全名 + 少量非人类）
TARGETS = [
    "Prostaglandin G/H synthase 2", "Prostaglandin G/H synthase 1", "Tumor necrosis factor",
    "Interleukin-6", "Interleukin-1 beta", "RAC-alpha serine/threonine-protein kinase",
    "Cellular tumor antigen p53", "Caspase-3", "Apoptosis regulator Bcl-2",
    "Apoptosis regulator BAX", "Transcription factor p65", "Nuclear factor NF-kappa-B p105 subunit",
    "Epidermal growth factor receptor", "Vascular endothelial growth factor A",
    "Mitogen-activated protein kinase 1", "Mitogen-activated protein kinase 14",
    "Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha",
    "Serine/threonine-protein kinase mTOR", "Hypoxia-inducible factor 1-alpha",
    "Heme oxygenase 1", "NAD(P)H dehydrogenase quinone 1", "Superoxide dismutase Mn",
    "Nuclear factor erythroid 2-related factor 2", "Peroxisome proliferator-activated receptor gamma",
    "Peroxisome proliferator-activated receptor alpha", "Insulin receptor",
    "Insulin-like growth factor 1 receptor", "Glycogen synthase kinase-3 beta", "Sirtuin 1",
    "Tumor necrosis factor receptor superfamily member 11B", "NF-kappa-B ligand receptor activator",
    "Matrix metalloproteinase-9", "Matrix metalloproteinase-2", "Tissue inhibitor of metalloproteinases 1",
    "Transforming growth factor beta-1", "C-C motif chemokine 2", "C-X-C motif chemokine 8",
    "Intercellular adhesion molecule 1", "Vascular cell adhesion protein 1", "Selectin E",
    "Angiotensin-converting enzyme", "Angiotensin II receptor type 1", "Endothelin-1 receptor",
    "Nitric-oxide synthase endothelial", "Coagulation factor Xa", "Prothrombin",
    "Estrogen receptor", "Androgen receptor", "Glucocorticoid receptor",
    "Cytochrome P450 3A4", "Cytochrome P450 2D6", "ATP-binding cassette sub-family B member 1",
    "Solute carrier organic anion transporter family member 1B1", "Aldose reductase",
    "Advanced glycosylation end product-specific receptor", "Fibronectin", "Collagen alpha-1(I) chain",
    "Cytochrome P450-cam", "mRNA of Protein-tyrosine phosphatase",
]

GENE_POOL = [
    "AKT1", "TNF", "IL6", "IL1B", "TP53", "CASP3", "BCL2", "BAX", "RELA", "NFKB1",
    "EGFR", "VEGFA", "MAPK1", "MAPK14", "PIK3CA", "MTOR", "HIF1A", "HMOX1", "NQO1", "SOD2",
    "NFE2L2", "PPARG", "PPARA", "INSR", "IGF1R", "GSK3B", "SIRT1", "TNFRSF11B", "TNFSF11",
    "MMP9", "MMP2", "TIMP1", "TGFB1", "CCL2", "CXCL8", "ICAM1", "VCAM1", "SELE", "ACE",
    "AGTR1", "EDNRA", "NOS3", "F10", "F2", "ESR1", "AR", "NR3C1", "CYP3A4", "CYP2D6",
    "ABCB1", "SLCO1B1", "AKR1B1", "AGER", "FN1", "COL1A1", "PTGS2", "PTGS1", "STAT3",
    "JUN", "FOS", "MYD88", "TLR4", "NLRP3", "CASP1", "IL18", "CCND1", "CDKN1A", "MYC",
    "CAT", "GPX1", "KEAP1", "PRKAA1", "PPARGC1A", "FOXO1", "FOXO3", "TSC2", "RHEB",
    "SLC2A4", "GCK", "PCK1", "G6PC", "SREBF1", "HMGCR", "LDLR", "APOB", "APOE",
    "SERPINE1", "PLAT", "PLAU", "THBD", "VWF", "KDR", "FLT1", "PGF", "ANGPT2", "TEK",
]

DISEASE_POOL = GENE_POOL + [
    "NEPH1", "NPHS1", "NPHS2", "PODXL", "WT1", "LAMA5", "ITGB1", "CD2AP", "ACTN4", "TRPC6",
    "CTGF", "PDGFB", "PDGFRB", "FGF2", "FGFR1", "JAK2", "STAT1", "SOCS1", "SOCS3", "PIAS1",
    "CASP8", "CASP9", "BCL2L1", "MCL1", "BIRC5", "XIAP", "DIABLO", "APAF1", "CYCS", "AIFM1",
    "PRKCA", "PRKCB", "PRKCD", "PLCG1", "PLCG2", "SOS1", "GRB2", "SHC1", "PTPN11", "CBL",
    "IRS1", "IRS2", "PDK1", "RPS6KB1", "EIF4EBP1", "ULK1", "BECN1", "MAP1LC3B", "ATG5", "ATG7",
    "SQSTM1", "NBR1", "CALCOCO2", "OPTN", "TBK1", "IRF3", "IFNB1", "ISG15", "MX1", "OAS1",
    "IL10", "IL4", "IL13", "IFNG", "TGFB2", "TGFBR1", "SMAD2", "SMAD3", "SMAD4", "SMAD7",
]

PATHWAYS = {
    "PI3K-Akt signaling pathway": ["PIK3CA", "AKT1", "MTOR", "IGF1R", "INSR", "GSK3B", "FOXO1",
                                   "TSC2", "RHEB", "RPS6KB1", "PDK1", "EGFR", "VEGFA", "KDR"],
    "AGE-RAGE signaling pathway in diabetic complications": ["AGER", "RELA", "NFKB1", "TNF", "IL6",
                                                             "CCL2", "VCAM1", "ICAM1", "SELE", "COL1A1",
                                                             "FN1", "TGFB1", "MAPK1", "MAPK14"],
    "TNF signaling pathway": ["TNF", "TNFRSF11B", "RELA", "NFKB1", "JUN", "FOS", "MAPK1", "MAPK14",
                              "PIK3CA", "AKT1", "CASP3", "CASP8", "BCL2", "MMP9"],
    "IL-17 signaling pathway": ["IL17A", "IL6", "CXCL8", "CCL2", "MMP9", "MMP2", "MAPK1", "MAPK14",
                                "RELA", "NFKB1", "FOS", "JUN", "TNF", "IL1B"],
    "HIF-1 signaling pathway": ["HIF1A", "VEGFA", "KDR", "FLT1", "NOS3", "SLC2A4", "GCK",
                                "SERPINE1", "EDNRA", "ANGPT2", "TEK", "MTOR"],
    "NF-kappa B signaling pathway": ["RELA", "NFKB1", "MYD88", "TLR4", "TNF", "IL1B", "NLRP3",
                                     "ICAM1", "VCAM1", "CCL2", "CXCL8", "BCL2L1"],
    "MAPK signaling pathway": ["MAPK1", "MAPK14", "JUN", "FOS", "MYC", "EGFR", "FGF2", "FGFR1",
                               "PDGFB", "PDGFRB", "PRKCA", "PRKCB", "TGFB1", "CASP3"],
    "FoxO signaling pathway": ["FOXO1", "FOXO3", "AKT1", "SIRT1", "PRKAA1", "PPARGC1A", "CAT",
                               "GPX1", "SOD2", "G6PC", "PCK1", "CDKN1A"],
    "Apoptosis": ["CASP3", "CASP8", "CASP9", "BCL2", "BAX", "BCL2L1", "TP53", "CYCS", "APAF1",
                  "DIABLO", "XIAP", "BIRC5", "MCL1", "AIFM1"],
    "Autophagy - animal": ["ULK1", "BECN1", "MAP1LC3B", "ATG5", "ATG7", "SQSTM1", "MTOR",
                           "PRKAA1", "TSC2", "RHEB", "FOXO3", "HIF1A", "BCL2", "CASP3"],
    "Pathways in cancer": ["AKT1", "TP53", "EGFR", "VEGFA", "MAPK1", "PIK3CA", "MTOR", "MYC",
                           "CCND1", "CDKN1A", "JUN", "FOS", "RELA", "NFKB1", "CASP3", "BCL2",
                           "BAX", "STAT3", "TGFB1", "MMP9", "MMP2", "CXCL8", "IL6", "TNF"],
    "Hepatitis B": ["RELA", "NFKB1", "TNF", "IL6", "CASP3", "CASP8", "TP53", "MYD88", "TLR4",
                    "IFNG", "STAT1", "CCL2", "CXCL8", "BCL2"],
    "Focal adhesion": ["FN1", "COL1A1", "ITGB1", "AKT1", "PIK3CA", "VEGFA", "KDR", "PDGFB",
                       "PDGFRB", "TGFB1", "MAPK1", "ACTN4", "FLT1", "TEK"],
    "Diabetic cardiomyopathy": ["AGER", "TGFB1", "FN1", "COL1A1", "MMP2", "MMP9", "TIMP1",
                                "PPARG", "INSR", "AKT1", "CASP3", "NQO1", "HMOX1", "SOD2"],
}


def w(name, rows, heads):
    path = os.path.join(HERE, name)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(heads)
        for r in rows:
            wr.writerow(r)
    print("wrote", name, len(rows))


def main():
    # 1) 成分表
    rows = []
    for i, (nm, ob, dl, hl, mw, ct, herbs) in enumerate(COMPOUNDS, 1):
        rows.append(["MOL%05d" % (1000 + i), nm, ob, dl, hl, mw, ct, herbs, "TCMSP"])
    w("compounds.csv", rows, ["MOL_ID", "compound", "OB", "DL", "HL", "MW", "content", "herb", "source"])

    # 2) 成分-靶点表（构造：每个成分 4–12 个靶点，泛成分连得更多）
    gene_of = {}
    for t in TARGETS:
        gene_of[t] = None
    # 蛋白名 → 基因（演示数据用的近似映射，真实分析请用 np_standardize + UniProt）
    PMAP = {
        "Prostaglandin G/H synthase 2": "PTGS2", "Prostaglandin G/H synthase 1": "PTGS1",
        "Tumor necrosis factor": "TNF", "Interleukin-6": "IL6", "Interleukin-1 beta": "IL1B",
        "RAC-alpha serine/threonine-protein kinase": "AKT1", "Cellular tumor antigen p53": "TP53",
        "Caspase-3": "CASP3", "Apoptosis regulator Bcl-2": "BCL2", "Apoptosis regulator BAX": "BAX",
        "Transcription factor p65": "RELA", "Nuclear factor NF-kappa-B p105 subunit": "NFKB1",
        "Epidermal growth factor receptor": "EGFR", "Vascular endothelial growth factor A": "VEGFA",
        "Mitogen-activated protein kinase 1": "MAPK1", "Mitogen-activated protein kinase 14": "MAPK14",
        "Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha": "PIK3CA",
        "Serine/threonine-protein kinase mTOR": "MTOR", "Hypoxia-inducible factor 1-alpha": "HIF1A",
        "Heme oxygenase 1": "HMOX1", "NAD(P)H dehydrogenase quinone 1": "NQO1",
        "Superoxide dismutase Mn": "SOD2", "Nuclear factor erythroid 2-related factor 2": "NFE2L2",
        "Peroxisome proliferator-activated receptor gamma": "PPARG",
        "Peroxisome proliferator-activated receptor alpha": "PPARA", "Insulin receptor": "INSR",
        "Insulin-like growth factor 1 receptor": "IGF1R", "Glycogen synthase kinase-3 beta": "GSK3B",
        "Sirtuin 1": "SIRT1", "Tumor necrosis factor receptor superfamily member 11B": "TNFRSF11B",
        "NF-kappa-B ligand receptor activator": "TNFSF11", "Matrix metalloproteinase-9": "MMP9",
        "Matrix metalloproteinase-2": "MMP2", "Tissue inhibitor of metalloproteinases 1": "TIMP1",
        "Transforming growth factor beta-1": "TGFB1", "C-C motif chemokine 2": "CCL2",
        "C-X-C motif chemokine 8": "CXCL8", "Intercellular adhesion molecule 1": "ICAM1",
        "Vascular cell adhesion protein 1": "VCAM1", "Selectin E": "SELE",
        "Angiotensin-converting enzyme": "ACE", "Angiotensin II receptor type 1": "AGTR1",
        "Endothelin-1 receptor": "EDNRA", "Nitric-oxide synthase endothelial": "NOS3",
        "Coagulation factor Xa": "F10", "Prothrombin": "F2", "Estrogen receptor": "ESR1",
        "Androgen receptor": "AR", "Glucocorticoid receptor": "NR3C1",
        "Cytochrome P450 3A4": "CYP3A4", "Cytochrome P450 2D6": "CYP2D6",
        "ATP-binding cassette sub-family B member 1": "ABCB1",
        "Solute carrier organic anion transporter family member 1B1": "SLCO1B1",
        "Aldose reductase": "AKR1B1", "Advanced glycosylation end product-specific receptor": "AGER",
        "Fibronectin": "FN1", "Collagen alpha-1(I) chain": "COL1A1",
        "Cytochrome P450-cam": "—", "mRNA of Protein-tyrosine phosphatase": "—",
    }
    ub = {"Quercetin": 30, "Kaempferol": 26, "Luteolin": 22, "beta-sitosterol": 18,
          "Stigmasterol": 14}
    rows = []
    for i, (nm, ob, dl, hl, mw, ct, herbs) in enumerate(COMPOUNDS, 1):
        k = ub.get(nm, rng.randint(4, 11))
        picks = rng.sample(TARGETS, min(k, len(TARGETS)))
        for t in picks:
            rows.append(["MOL%05d" % (1000 + i), nm, t, herbs.split("|")[0], "TCMSP"])
    # 保证非人类条目出现在表里
    rows.append(["MOL01005", "beta-sitosterol", "Cytochrome P450-cam", "黄芪", "TCMSP"])
    rows.append(["MOL01007", "Calycosin", "mRNA of Protein-tyrosine phosphatase", "黄芪", "TCMSP"])
    w("compound_targets.csv", rows, ["MOL_ID", "compound", "target_name", "herb", "source"])

    # 3) 疾病靶点（含 score 与来源，score 长尾分布）
    rows = []
    for g in DISEASE_POOL:
        if rng.random() < 0.18:
            continue
        score = round(min(99.0, max(0.4, rng.gauss(12, 9))), 2)
        src = rng.choice(["GeneCards", "GeneCards", "OMIM", "DisGeNET", "CTD"])
        rows.append([g, score, src])
    w("disease_targets.tsv", rows, ["gene", "score", "source"])

    # 4) PPI 边表（STRING 风格三列 + 自环 + 重复）
    genes = sorted({PMAP[t] for t in TARGETS if PMAP[t] != "—"} | set(list(PATHWAYS.values())[0]))
    edges = set()
    for pname, members in PATHWAYS.items():
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                if rng.random() < 0.30:
                    edges.add(tuple(sorted((members[i], members[j]))))
    for _ in range(60):
        a, b = rng.sample(genes, 2)
        edges.add(tuple(sorted((a, b))))
    rows = []
    for a, b in sorted(edges):
        rows.append([a, b, round(rng.uniform(0.42, 0.99), 3)])
    rows.append([rows[0][0], rows[0][0], 0.99])          # 自环
    rows.append(list(rows[1]))                            # 重复边
    w("string_edges.tsv", rows, ["node1", "node2", "combined_score"])

    # 5) GMT
    with open(os.path.join(HERE, "pathways.gmt"), "w", encoding="utf-8") as f:
        for p, mem in PATHWAYS.items():
            f.write("\t".join([p, "demo"] + mem) + "\n")
    print("wrote pathways.gmt", len(PATHWAYS))

    # 6) 对接结果
    pairs = [("Quercetin", "AKT1", -8.4, 1.32), ("Kaempferol", "TNF", -7.6, 1.55),
             ("Luteolin", "IL6", -6.9, 1.71), ("beta-sitosterol", "PPARG", -9.1, 1.18),
             ("Tanshinone IIA", "AKT1", -8.9, 1.24), ("Cryptotanshinone", "MAPK1", -7.2, 1.63),
             ("Calycosin", "ESR1", -7.8, 1.44), ("Ferulic acid", "PTGS2", -5.6, 1.90),
             ("Z-ligustilide", "NFKB1", -5.1, 2.60), ("Salvianolic acid B", "MMP9", -8.0, 1.51),
             ("Astragaloside IV", "TGFB1", -6.4, 1.80), ("Formononetin", "AGER", -7.4, 1.47)]
    w("docking_scores.csv", [[c, g, s, r] for c, g, s, r in pairs],
      ["compound", "target", "affinity", "redock_rmsd"])

    # 7) 成分-靶点长表（供网络/富集）
    print("done")


if __name__ == "__main__":
    main()
