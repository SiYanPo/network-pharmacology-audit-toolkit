# 外部证据核查记录（2026-11-22）

> 对应评审报告 §14 的四个“需外部求证点”。本文件只记录**能查到的证据**与**不能下的结论**，
> 不把网络搜索当作决定性证明。所有措辞已同步到 `references/33-limitations-and-wording.md`。

## 1. ETCM2 与 TCMSP 的谱系关系

**结论：官方一手论文不支持“ETCM2 是 TCMSP 子集”的强断言；但也不证明两者独立。**

一手来源：Zhang Y et al. *ETCM v2.0: an update with comprehensive resource and rich annotations for traditional Chinese medicine*. Acta Pharm Sin B. 2023;13(6):2559-2571. DOI `10.1016/j.apsb.2023.03.012`; PMCID `PMC10326295`; PMID `37425046`。

可核实事实：

- ETCM v2.0 的靶点识别使用 **D3CARP 二维配体相似性搜索**，参考配体库来自 **BindingDB（updated 2021-11-01）**；相似度为 1.0 记为 confirmed，0.8–1.0 记为 potential。
- 论文 Table 2 把 TCMSP 作为**同类数据库**比较：TCMSP 499 味药材、29,384 成分、3,311 靶点、837 疾病；ETCM2 为 2,079 药材、38,298 成分、1,040 靶点、8,045 疾病。
- 论文来源列表**没有 TCMSP**；没有“ETCM2 使用 TCMSP 作为上游”或“ETCM2 是 TCMSP 子集”的陈述。
- 2025 综述（PMID 见 `PMC12833950`）指出二者都从《中国药典》等权威来源与 PubChem/ChEMBL 标准标识取材，说明**共享权威上游**是可能的，但这不是“子集”证明。

技能中的正确处理：

- 保留 `np_split.py lodo` 的本地边集包含检验（本轮实测 ETCM2 502/502 边落在本项目 TCMSP 边表内）——它是**本项目数据切片的观察**。
- 论文措辞：不能写“ETCM2 is an independent database relative to TCMSP”；应写
  “local edge-set containment was observed; the published ETCM v2.0 paper does not establish database lineage, so LODO is reported as cross-source consistency, not independence.”

## 2. STRING v12 的 interaction partners / 低度节点补边偏差

**结论：官方 STRING v12 论文支持“inherent literature bias 存在且被部分降低”，但不支持“interaction_partners 对低度节点有系统性补边偏差”的强断言。**

一手来源：Szklarczyk D et al. *The STRING database in 2023*. Nucleic Acids Res. 2023;51(D1):D638-D646. PMID `36370105`; PMCID `PMC9825434`。

可核实事实：

- STRING v12 有 7 类证据通道（neighborhood/fusion/co-occurrence/co-expression/experiments/database/textmining），各通道有 sub-score，最终合并为 combined confidence score。
- textmining 从 BioBERT 升级为 biomedical RoBERTa-large；v12 与 v11.5 的物理相互作用计数在阈值 0.4/0.7/0.9 均明显变化，说明**版本间分数统计不可直接混用**。
- 官方原文：RNA-Seq “reduced the inherent literature bias of the STRING network”，并提到对低表达基因敏感性更高；这是“降低偏倚”而非“证明无偏”。
- 官方论文没有讨论 `interaction_partners` API 的补边行为、degree 分布或 low-degree bias；文中 `neighborhood` 是基因组邻域证据通道，不是网络邻居扩展。

技能中的正确处理：

- `references/33-limitations-and-wording.md` 已改成：扩展 subnet 的 degree 含研究偏倚、**不参与 hub 判定**；若需要“低度节点系统性补边偏差”的结论，必须自己做 degree-bin 返回率实验，不能引用 STRING v12 论文替我们背书。
- `np_ppi.py` 的 `--expand-isolated` 保留证据分层（experiments/coexpression）与 hub rank 相关；这是把风险量化，而不是把风险消除。

## 3. HPA rna_tissue_consensus 的版本与命名稳定性

**结论：HPA 确实是持续版本化资源，组织/细胞命名、样本集合、资源结构会跨版本变化；必须记录 release 与 Ensembl 版本。**

一手来源：Human Protein Atlas Release history（访问 2026-11-22）：`https://www.proteinatlas.org/about/releases`。

可核实事实：

- 当前页最新版本为 **v26.0**，发布日期 **2026-11-11**，Ensembl 版本 **109**；v25.1 为 2026-05-25，v25.0 为 2025-11-15，v24.1/v24.0/v23.0/v22.0 均有记录。
- v26.0 重排 brain/immune cell 资源，更新 cell type naming，移除 Tabula Sapiens 命名比较，重分类 epicardial cells；v23.0 修正脑区命名与分组；v25.0 替换全部 tissue/cancer 图像并重染。
- v26.0 因质控移除 2,050 个基因的抗体数据、改变 665 个基因 reliability score；这些变化直接影响力“可复现”的表述。

技能中的正确处理：

- `np_snapshot.write_meta()` 已记录 `release` / `fetched_at` / `sha256` / `etag`。
- 新增要求：HPA 相关 meta 必须额外写 `resource=rna_tissue_consensus` 与 `ensembl_release=109`（或实际使用版本）；报告不得只写 “HPA” 而不写版本。
- `heart` 用于 PD 演示仍属组织-疾病错配，演示数据只作管线验证，不作生物学解读。

## 4. TCM 方剂金标准是否存在可引用 Delphi 研究

**结论：存在疾病特异的中医证候—方剂 Delphi 共识，但没有可直接替代本项目通用方剂金标准的“通用 Delphi 金标准”。**

可核实文献：

- Guo Q et al. *A Delphi Study of Expert Consensus on Chinese Medicine Syndrome Differentiation and Herbal Use for Early Breast Cancer*. Integr Cancer Ther. 2023. PMCID `PMC10557418`。30 位高级中医临床专家，3 轮；证候 ≥50% 专家选择即共识，症状/体征中位数 ≥4 且 ≥70% 专家评 4/5。对应方剂是该证候下最常用方剂，**不是独立阈值方剂金标准**；适用早期乳腺癌五个治疗阶段。
- Chinese medicine practitioners’ consensus ... COVID-19 survivors: A Delphi study. J Integr Med. 2024；13 位中医 practitioners，3 轮，consensus = ≥80% agreement 且 IQR≤1；**仅适用于 COVID-19 康复期**。
- 中医古代经典方证据质量评级量表（2025，WPRIM 记录）用 Delphi 形成 3 维度 12 条目量表，CVI 0.98、ICC 0.90；它是**证据质量评价工具**，不是方剂靶点金标准。

技能中的正确处理：

- `references/34-validation-roadmap.md` V1 保留“外部方剂 + 逐条 PMID/URL 金标准”；新增说明：可优先引用疾病特异 Delphi 共识作为**适应症/证候侧证据**，但仍需 DrugBank/TTD/ChEMBL 机制记录构成靶点侧证据，二者不可互相替代。
- 不得把任一 Delphi 研究写成“通用 TCM 方剂 gold standard”。

## COLOC 与 R `coloc` 对齐

见 `audit/coloc_vs_R_coloc_2026-11-22.md`。该文件实测结论：修复 ABF 公式与 H4 组合后，20 位点 r=1.0000 / MAE=0.0000（verdict PASS），状态 **DONE（模拟数据对齐）**；DoD 第 7 条已完成。本文件不替代数值一致性报告。真实 GWAS 抽验仍列在 `references/34-validation-roadmap.md` V3（属研究计划，不是技能缺陷）。
