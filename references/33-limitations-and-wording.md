# 边界与措辞（评审 §8 / §9）

## 必须写进 Methods/Limitations 的五条

1. **COLOC 已按 R `coloc` 源码对齐（2026-11-22）。** 当前先验 `p1=p2=1e-4, p12=1e-5, W=0.15²`；
   20 位点强信号与部分非饱和集均 PASS（r=1.000, MAE=0.000），证据见
   `audit/coloc_vs_R_coloc_2026-11-22.md`。Methods 可写
   "Wakefield ABF with coloc.abf-equivalent combination, validated against coloc 5.2.3 on 20 simulated loci"；
   真实 GWAS 仍需抽位点复核 LD/MAF 与表型类型（quant sd=0.15 / cc sd=0.2）。
2. **HPA 组织背景必须与疾病组织匹配。** `heart` 用于 PD 演示是组织-疾病错配；
   demo 改用 `brain`，或明确标"仅演示，不作生物学解读"。
3. **zebrafish 1→many 不是功能等价。** Ensembl one2many 只可报计数，不得据此给实验建议。
4. **MR `--drug-direction` 只做方向筛查。** 单靶点方向忽略细胞类型、部分激动、反馈；
   不得写"药物有效"，只能写"遗传证据方向与假设一致/矛盾"。
5. **examples 为合成数据。** 任何示例结果只验证管线，不作生物学结论。

## 扩展子网

`--expand-isolated` 的 degree 含 STRING 研究偏倚，**不参与 hub 判定**；主网络与扩展子网
必须分开落盘、分开标题。报告需给出扩展前后 hub rank Spearman；若 ρ 低，说明原 hub 结论不稳健。

**外部求证校正（2026-11-22）**：STRING v12 官方论文（PMID 36370105）只支持
"STRING 网络有 inherent literature bias，新增表达数据与实验评分降低但并未消除偏倚"；
它**没有**证明 `interaction_partners` 对低度节点存在系统性补边偏差。因此：
- 可以写"扩展 degree 受 STRING 证据构成与文献偏倚影响"；
- 不能写"STRING 官方承认低度节点补边偏差"；
- 若要下后一个结论，必须自己做 degree-bin 返回率/分数分布实验。

## 多重性与校准

- 只有 `q_global<0.05` 可写"显著"；`p<0.05` 只能写"提示性"。
- 校准概率只在拟合基准定义域内有效；换疾病/数据库必须重拟合。
- OOF AUC<0.60 或系数全 0 → `UNCALIBRATED`，不得写 `P(gold)`。

## 金标准独立性

ETCM2 与 TCMSP 的边集包含率若 >0.95，必须先查谱系；在查清前，LODO 的"留一库"
对照只能写"跨库一致性"，不能写"独立验证"。

**一手来源校正（2026-11-22）**：ETCM v2.0 论文（DOI 10.1016/j.apsb.2023.03.012,
PMCID PMC10326295）说明其靶点来自 D3CARP 对 BindingDB 的二维相似性搜索，Table 2
把 TCMSP 列为同类比较库；论文**未声明** ETCM2 使用 TCMSP 或为其子集。故：
- 本地观察到"ETCM2 502/502 边落在本项目 TCMSP 边表内"是可报告的**本地切片现象**；
- 不得据此写"ETCM2 is a subset of TCMSP"或"ETCM2 is independent of TCMSP"；
- LODO 统一写 **cross-source consistency**，不得写 independent validation。

TCM 方剂金标准方面：已有疾病特异 Delphi 共识（如早期乳腺癌 PMCID PMC10557418，
30 位专家、3 轮；证候 ≥50% 阈值；方剂为对应的最常用方剂），但**没有通用方剂靶点金标准**。
Delphi 可用于适应症/证候侧证据，不能替代 DrugBank/TTD/ChEMBL 的机制靶点证据。

HPA 版本方面（官方 release history，2026-11-22 访问）：HPA v26.0 / 2026-11-11 / Ensembl 109；
v22→v26 有细胞命名、组织样本、可靠性评分变更。报告必须写 `HPA vX + Ensembl release`，
不得只写 "HPA"。

外部证据全文见 `references/36-external-evidence-audit.md`。
