# 外部验证与湿实验路线（评审 §11）

> 这不是“已完成的验证”，而是可执行、可预注册的验证协议。
> 任何一条未完成前，论文只能写“internal benchmark / hypothesis-generating”。

## V1 外部方剂验证（2–3 个，2–4 周）

- 选择与开发集**不重叠**的方剂：如黄连解毒汤、血府逐瘀汤、连花清瘟（至少 2 个）。
- 金标准来源：DrugBank/TTD/ChEMBL 机制记录 + 指南/共识，逐条记录 PMID/URL。
  Delphi 共识可作为**证候/适应症侧证据**，不能替代机制靶点证据：例如早期乳腺癌
  Delphi（PMCID PMC10557418，30 位专家、3 轮、证候阈值 ≥50%）只适用于该疾病场景；
  详见 `references/36-external-evidence-audit.md` §4。
- 与开发集做 **LODO**：一个库建边，另一个库做金标准；先跑 `np_split.py lodo`
  检查边集包含率 >0.95 的“伪独立”库。
- 报告 `top20_hits`、AUPRC、ECE/Brier、coverage；AUC 只作次要。

## V2 时间留出（temporal hold-out）

- 以论文/数据库年份切分：≤2022 建边，2023–2025 金标准。
- 禁止用未来文献调参；预注册阈值与 baseline（`config/benchmark_prereg.yaml`）。
- 输出时间漂移曲线：每年 top20 命中数与校准误差。

## V3 COLOC 与 R `coloc` 对齐（1 周）— 模拟数据部分 DONE（2026-11-22）

- ✅ 已完成：抽 20 个模拟位点，跑官方 `coloc.abf` 与当前 ABF 近似；修复 ABF 公式与 H4 组合后
  r=1.0000 / MAE=0.0000（verdict PASS），另用 20 个非饱和位点防 0/1 饱和刷分。
  证据 `audit/coloc_vs_R_coloc_2026-11-22.md`。Methods 现可写 “validated against coloc.abf on 20 simulated loci”。
- ✅ 真实 GWAS 抽验（2026-11-23 完成，工程层）：直接用 DoD 真实 MR 的 `exp.tsv`/`out.tsv` 喂
  `np_coloc_compare.py`，3 靶点（IL1B/BCL2/PTGS2）PP.H4 与 R coloc 最大绝对差 < 3×10⁻⁵。
  证据：`DoD收尾/colic_real/coloc_real_gwas_validation.md`。
- ⏳ 仍未做（需整段基因座数据，非技能缺陷）：PASS 级 verdict 需单个基因座 ±500kb 全 SNP 的
  暴露+结局汇总统计（MR top-hits 跨染色体，不构成单一 LD 块）；以及 Bland–Altman 偏差、方向不一致清单。
- 边界：Methods 仍只写 “ABF-based approximation”，不写 “we performed colocalisation”（除非拿到 PASS 级 verdict）。

## V4 通路层基准（2 周）

- 除基因集外，比较 KEGG 通路重叠：Jaccard、超几何 p、排名前 10 通路命中。
- 回应“只对基因不对机制”的审稿意见。

## V5 不确定性传播（2 周）

- Bootstrap 整个流程（重采样成分/边/金标准），输出每个靶点排名 95% CI。
- 报告“排名在前 20 但 CI 跨 20 名以外”的比例；这是比单点排名更诚实的指标。

## V6 湿实验闭环（3–6 月，1 方剂 × 3 靶点）

- 选择 `wetlab_candidate` 且 P(gold) CI 不跨 0.5 的 3 个靶点。
- 至少 qPCR/WB 两项；阴性结果同样登记（负结果默认值见 `references/NEGATIVE_RESULT_DEFAULTS.yaml`）。
- 预注册假设、细胞系、抗体货号、批号；阴性也要写。

## V7 用户研究（1 月，5–10 位研究生）

- 任务：用 3 个案例跑流程并解释“为什么拒答/为什么不确定”。
- 记录误用类型（把排序当概率、把扩展 degree 当 hub、把 p 当 q）。
- 作为软件论文可用性证据；不替代方法学验证。
