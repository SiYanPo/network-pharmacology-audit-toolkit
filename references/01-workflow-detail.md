# 11 阶段逐步详解（判据、输入输出、失败恢复）

> 每个阶段统一四问：**输入是什么 → 做什么 → 产出什么 → 什么情况必须停**。
> 对应脚本见 SKILL.md §7；所有阶段共享 `<outdir>/`，通过固定文件名通信。

---

## 01 靶点标准化 `np_standardize.py`

**输入**：任意含 `target_name` / `gene` 列的表，或纯基因列表（txt/csv/tsv/xlsx）。

**做四层匹配**

| 层 | 规则 | 典型命中 |
|---|---|---|
| as_symbol | 输入已是符号形式（`AKT1`） | 已整理过的表 |
| L1 exact | 与内置词典规范化后完全一致 | `Prothrombin`→F2、`Tumor necrosis factor`→TNF |
| L2 normalized | 去标点/括号/希腊字母/亚型后缀 | `Tumor necrosis factor (TNF-alpha)` |
| L3 structural | 拆 `(别名)`、`[Cleaved into:]`、`Includes:` | `Alcohol dehydrogenase 1B` |
| L4 alias | 人工别名表 | 词典未覆盖的地方药材靶点 |
| nonhuman / noise | 菌/病毒/动物/mRNA/clone | `Cytochrome P450-cam`、`mRNA of ...` |

**输出**：`targets_standardized.csv`（保留原始列 + symbol/layer/note）、`targets_unmatched.csv`、
`targets_gene_list.txt`、`standardize_audit.md`。

**必须停**：映射率 <80% → 补别名表（`--alias my_alias.csv`）重跑；**不得**把非人类蛋白强行映射。

---

## 02 成分筛选 `np_compounds.py`

**输入**：成分表（至少 `compound`, `OB`, `DL`；可选 `HL`, `content`, `herb`）。

**三道闸**

1. **阈值闸**：OB/DL（+可选 HL、含量），`--logic and|or`。
2. **途径闸**：`--route oral|injection|topical|unknown`。注射/外用 + OB 筛 = 逻辑冲突，
   默认**阻断**；`--force` 可继续但会写警告与方法学声明。
3. **泛成分闸**：槲皮素/山柰酚/木犀草素/β-谷甾醇等标记 `ubiquitous=Y`，并算占比。

**输出**：`compounds_kept.csv` / `compounds_exception.csv`（例外纳入，须有实测或文献依据）/
`compounds_rejected.csv`（含剔除理由）/ `compounds_report.md`。

**必须停**：泛成分占比 >50% → 单列 + 用非泛成分做敏感性分析；全被筛掉 → 阈值不适配，
改含量/入血成分口径。

---

## 03 交集与成分级富集 `np_intersect.py`

**输入**：`--drug`（成分-靶点长表）+ `--disease`（疾病靶点表），或 `--sets A.txt,B.txt`。

**做三件事**
1. N 组交集 + 韦恩图（2–4 组椭圆，>4 组退化为条形）。
2. 每个交集靶点带疾病侧**证据强度**（score / 来源数 / 来源类型→4 个层级）。
3. 每个成分做**精确超几何**：它的靶点撞上疾病基因是否超出随机。同时给"全基因组"与
   "可检出空间"两套背景（见 `04-statistics.md` §1.1）。

**输出**：`intersection.csv`、`intersection_genes.txt`、`venn.svg`、
`compound_enrichment.csv`、`intersect_report.md`。

**必须停**：交集 0 或 <5 → 报告内附**六步恢复序列**（ID 空间→物种→疾病集过窄→成分过少→
疾病名不规范→接受阴性结果）。**不要**靠降阈值凑数。

---

## 04 PPI 与 hub `np_ppi.py`

**输入**：边表（STRING API / 网页导出 / 无表头三列 / CSV 均可）。

**处理链**：清洗（自环/重复/阈值）→ 可选最大连通分量 → 8 项拓扑指标 →
多指标秩均值共识 → Louvain 模块 → 边 bootstrap 入选频率 → 节点丢弃 → 度匹配零模型 → SVG。

**输出**：`ppi_edges_clean.tsv`、`topology.csv`、`hubs.csv`（含 verdict）、`network.svg`、`ppi_report.md`。

**必须停**：无稳健 hub（全部 <0.90）→ 正文不得写"核心枢纽"；度匹配零模型 P>0.05 →
hub 可能只是度大，结论降级。

---

## 05 零模型 `np_null.py`

**输入**：交集基因 + PPI 边表 + 疾病表（+ 可选参照成分库）。

**输出**：`null_tests.csv`、`null_report.md`（含可直接粘贴的方法学句子）。

**必须停**：B 组 P≥0.05 → "功能模块"表述作废；C 组不显著 → 紧密性可由度解释。

**陷阱**：不要构造"保持双边际的边交换 + 聚合统计量"这组合（数学上不变，P 恒为 1），
那是假检验。见 `04-statistics.md` §3.2。

---

## 06 富集 `np_enrich.py`

**输入**：基因列表 +（本地 GMT 或 `--gprofiler`）+ 可选 `--background`。

**处理**：超几何 → BH → 最小命中数 → 冗余折叠（Jaccard）→ 三层分级（机制明确/伞状泛化/
背景性/未分类）→ 气泡图 → 报告。

**输出**：`enrich_all.csv`、`enrich_significant.csv`（含 collapsed_n）、`enrich_bubble.svg`、
`enrich_report.md`。

**必须停**：无背景声明 → 只能写"探索性"；g:Profiler 不可用 → 记录为"数据源不可用"，
**不得**写成"该通路不富集"。

---

## 07 剂量可达性 `np_dose.py`

**输入**：饮片日剂量、成分含量、分子量（+可选 F/Vd/fu/提取率/体外浓度）。

**链条**：`g/日 × 含量% × 提取率% → mg → ÷分服 → ×F → ÷(Vd·体重) → µM → ×fu → 游离 µM`。

**输出**：`dose_result.json`、`dose_report.md`（含三档敏感性区间与"未计算部分"声明）。

**必须停**：体外/体内比值 >10 → 机制结论降级；反推所需剂量 >30 g/日 → 该浓度临床不可达。

---

## 08 对接门禁 `np_dock.py`

**输入**：CSV（compound,target,affinity[,redock_rmsd]）或 Vina 输出目录。

**三道门**：G1 重对接 RMSD ≤2 Å；G2 打分口径（预测非实测）；G3 覆盖说明（只对接部分须给依据）。

**输出**：`docking_graded.csv`（含 target_gate）、`docking_bar.svg`、`docking_report.md`。

**必须停**：某靶点 RMSD >2 Å → 该靶点全部打分不得用于机制论证；未做验证 → 标"探索性"。

---

## 09 反向网络药理学 `np_optimize.py`

**输入**：成分-靶点表（含 herb/compound 列）+ 疾病靶点表。

**指标**：coverage（广度）、precision（命中率）、marginal（leave-one-out 互补性）、fold（富集倍数）；
再做贪心最小覆盖 + 协同增量。

**输出**：`reverse_ranking.csv`、`combination_trace.csv`、`coverage_bar.svg`、`optimize_report.md`。

**必须停**：协同增量 = 0 → 不得称"多靶点协同"；所有 precision <0.3 → 靶点噪声过大，
先回去做标准化与多库交叉。

---

## 10 红线审计 `np_audit.py`

对目录内所有报告与 manifest 做 26 条规则判定，输出 `audit/audit_report.md` + `audit.json`。

**评分**：`(PASS + 0.5×WARN) / 26 × 100`。投稿前建议 `--strict`（WARN 也当失败）。

**退出码**：有 FAIL → 2；`--strict` 且有 WARN → 1；否则 0。

---

## 11 报告 `np_report.py`

把各阶段产物合成单文件 HTML：一页结论 → 证据等级图例 → 参数与溯源快照 → 图表（内嵌 SVG）→
分阶段详情（折叠）→ **自动生成的方法学段落** → 审计摘要 → 流程警告 → 局限。

**注意**：方法学段落只保证"参数与流程一致"，数据库版本与访问日期仍需作者按实际填写。

---

## 阶段依赖图

```
targets.csv ─┐
             ├─01→ targets_standardized.csv ─┐
compounds.csv┘                                ├─03→ intersection_genes.txt ─┬─05 (需 edges+disease)
                                              │                             ├─06 (需 gmt/在线)
disease.tsv ──────────────────────────────────┘                             └─09
edges.tsv ──────────────────────04→ hubs.csv ─────────────────────────────────┘
docking.csv ──08（独立）
dose 参数 ────07（独立，但建议放在 03 之后，用交集靶点/成分）
所有产物 ────10→ audit ────11→ np_report.html
```
