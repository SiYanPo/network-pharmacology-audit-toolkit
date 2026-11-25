# 交集与成分富集报告（药物-疾病）

| 集合 | 规模 |
|---|---|
| 药物靶点 | 101 |
| 疾病靶点 | 2244 |
| **交集** | **23** |

## 交集靶点的证据构成

| 证据层级 | 靶点数 | 解读 |
|---|---|---|
|  curated 关联 | 23 | 较强（人工审编数据库） |

> 写法提醒：GeneCards 的 score 是**整合相关性**，不是因果证据；GEO 差异表达是「特定样本 vs 对照」的变化，不等于「疾病靶点」。报告中应按上表分层陈述，不能把三者混为「疾病相关基因」。

## 交集靶点（按共享成分数排序，前 15）

| 基因 | 共享成分数 | 疾病打分 | 来源 | 证据层级 |
|---|---|---|---|---|
| DPP4 | 1 | 0.6 | Open Targets MONDO_0005148 |  curated 关联 |
| FTO | 1 | 0.6 | Open Targets MONDO_0005148 |  curated 关联 |
| MGAM2 | 1 | 0.6 | Open Targets MONDO_0005148 |  curated 关联 |
| IGF1R | 1 | 0.5 | Open Targets MONDO_0005148 |  curated 关联 |
| CAMK2B | 1 | 0.5 | Open Targets MONDO_0005148 |  curated 关联 |
| PTGS2 | 1 | 0.4 | Open Targets MONDO_0005148 |  curated 关联 |
| PNLIP | 1 | 0.4 | Open Targets MONDO_0005148 |  curated 关联 |
| CA2 | 1 | 0.4 | Open Targets MONDO_0005148 |  curated 关联 |
| CA4 | 1 | 0.4 | Open Targets MONDO_0005148 |  curated 关联 |
| PIK3R1 | 1 | 0.4 | Open Targets MONDO_0005148 |  curated 关联 |
| CYP19A1 | 1 | 0.3 | Open Targets MONDO_0005148 |  curated 关联 |
| XDH | 1 | 0.3 | Open Targets MONDO_0005148 |  curated 关联 |
| MMP13 | 1 | 0.3 | Open Targets MONDO_0005148 |  curated 关联 |
| ALK | 1 | 0.3 | Open Targets MONDO_0005148 |  curated 关联 |
| BRAF | 1 | 0.3 | Open Targets MONDO_0005148 |  curated 关联 |

## 成分级富集（谁的靶点最集中在疾病上）

- 检验：精确超几何；背景 N=20596（genome），疾病集 K=2244
- 显著成分（q<0.05）：1 / 1

| 成分 | 靶点数 | 共享 | 期望 | 倍数 | q 值 |
|---|---|---|---|---|---|
| Quercetin | 101 | 23 | 11.0 | 2.09 | 0.0004 |

> 注意：这是「靶点是否富集于疾病基因集」的统计，**不是**「该成分有疗效」。靶点数多的成分（常为泛成分）天然更容易显著，解释时须结合含量与实测证据。

