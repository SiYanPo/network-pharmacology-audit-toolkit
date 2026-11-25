# 富集分析报告

- 输入基因：23；背景：默认全基因组（N=20596）
- 最小命中数：3；显著性：q<0.05（BH 校正）
- 显著条目：141；折叠后代表条目：23（Jaccard 阈值 0.50）

<div style='border-left:4px solid #b91c1c;background:#fdecec;padding:8pt 10pt;margin:8pt 0'><b style='color:#b91c1c'>P-value may be inflated / P 值可能膨胀</b><br>本次富集未使用疾病相关的表达谱（DEGs）或组织表达基因集作为背景（当前：全基因组（未使用疾病/组织表达背景），N=20596）。全基因组背景会系统性放大富集显著性，本节的 P/q 只能作为<b>探索性</b>证据，不得写入结论段。<br>修正方式：<code>python np_enrich.py -g genes.txt -m x.gmt --background degs.txt</code>（首选，实测 DEGs），或 <code>--background-mode hpa --hpa-tissue liver</code>（次选，HPA 组织表达）。</div>

## 三层分级（决定哪些能写进结论）

| 层级 | 条目数 | 使用规则 |
|---|---|---|
| 机制明确 | 2 | 可作为机制假说的锚点 |
| 伞状泛化 | 3 | 只说明功能谱广，**不得**作为机制结论 |
| 背景性（感染类） | 0 | 与核心通路共享分子，属被动命中 |
| 未分类 | 18 | 需人工判断 |

### 机制明确（2）

| 条目 | 命中 | 条目大小 | q 值 | 折叠 |
|---|---|---|---|---|
| IL-17 signaling pathway | 5 | 94 | 1.66e-06 | +3 |
| HIF-1 signaling pathway | 4 | 109 | 2.73e-05 | +1 |

### 伞状泛化（3）

| 条目 | 命中 | 条目大小 | q 值 | 折叠 |
|---|---|---|---|---|
| Disease R-HSA-1643685 | 9 | 1736 | 0.0001 | +5 |
| Metabolism R-HSA-1430728 | 9 | 2049 | 0.0004 |  |
| Metabolism Of Lipids R-HSA-556833 | 4 | 732 | 0.0092 |  |

### 未分类（18）

| 条目 | 命中 | 条目大小 | q 值 | 折叠 |
|---|---|---|---|---|
| Prostate cancer | 7 | 97 | 1.41e-09 | +41 |
| Reversible Hydration Of Carbon Dioxide R-HSA-1475029 | 3 | 12 | 3.93e-06 | +2 |
| Non-small cell lung cancer | 4 | 72 | 8.55e-06 |  |
| ErbB signaling pathway | 4 | 85 | 1.39e-05 | +12 |
| Alzheimer disease | 6 | 369 | 1.41e-05 | +13 |
| Small cell lung cancer | 4 | 92 | 1.63e-05 | +10 |
| Negative Regulation Of Transferase Activity (GO:0051348) | 3 | 37 | 3.90e-05 | +13 |
| Bladder cancer | 3 | 41 | 4.79e-05 | +1 |
| Insulin Receptor Signaling Pathway (GO:0008286) | 3 | 47 | 6.43e-05 | +11 |
| Ovarian steroidogenesis | 3 | 51 | 7.39e-05 |  |
| Protein Autophosphorylation (GO:0046777) | 4 | 162 | 8.43e-05 | +1 |
| Cellular Response To Reactive Oxygen Species (GO:0034614) | 3 | 73 | 0.0002 |  |
| PD-L1 expression and PD-1 checkpoint pathway in cancer | 3 | 89 | 0.0003 | +2 |
| Cellular Response To Heat Stress R-HSA-3371556 | 3 | 99 | 0.0004 | +1 |
| Cholinergic synapse | 3 | 113 | 0.0005 | +2 |
| Positive Regulation Of Phosphorylation (GO:0042327) | 3 | 231 | 0.0027 |  |
| Positive Regulation Of Cellular Process (GO:0048522) | 4 | 594 | 0.0046 |  |
| Regulation Of Cell Population Proliferation (GO:0042127) | 4 | 766 | 0.0106 |  |

## 写法红线

1. **P 值是显著性，不是效应值**：不能写成「作用强度」。效应量看倍数富集与命中数。
2. **GeneRatio / Rich factor 只是比值**（命中数/条目大小），不是通路活性评分。
3. **P 显著 ≠ 通路被激活**：GO/KEGG 的输入是无方向的基因集，没有激活/抑制信息。
4. **不得只挑「好看」的通路放图**：抗炎研究不能只放 TNF/NF-κB，应给完整表（附录）+ 讨论解释。
5. 冗余条目已折叠：折叠前后条目数必须都报告（141 → 23）。

