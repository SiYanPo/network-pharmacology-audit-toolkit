# 32 · 证据模型与基准改造（v3.0）

本文件记录"为什么这样改"以及"改完怎么读"。所有数字都来自本仓库实测，不是示意图。

---

## 1. 证据账本（`np_evidence_ledger.py`）

### 1.1 为什么需要账本

旧版有三套互不相干的排序：交集表按支持成分数、hub 表按拓扑共识、MR 按 β。
方法学审稿人的第一问是"你的最终排序函数是什么"。账本的作用是**把三套排序降格为特征**，
再合成一个可校准的映射：

```
evidence ledger (gene × evidence_type) → features → calibrated P(gold) → decision
```

### 1.2 字段与分级

| 字段 | 说明 |
|---|---|
| `gene` | 大写 HGNC 符号 |
| `source` | 产生该证据的模块（`compound_target_network` / `string_ppi` / `mendelian_randomization` …） |
| `evidence_type` | 细类（`compound_support`、`disease_association`、`network_topology`、`pathway_membership`、`causal_estimate`、`colocalisation`、`binding_prediction`、`wetlab_feasibility`、`direction_conflict`） |
| `level` | L0 文本挖掘 / L1 计算预测 / L2 数据库注释 / L3 实验支持 / L4 遗传因果 |
| `direction` | 方向（`support`/`association`/`risk_increasing`/`conflict`…）。方向冲突记 `-1` 效应 |
| `effect_size` / `uncertainty` | 该证据的量（成分数、β、PP.H4、对接分…）与不确定度 |
| `provenance_hash` | 来源文件 SHA256 前 16 位——**账本行可回溯到具体文件** |

> 级别是**先验强度**，不是真值。模型会学到各等级的实际权重；若某级系数为负，
> 说明该级证据在本基准上是噪声，报告必须如实展示（见 §1.4）。

### 1.3 校准与"拒绝给概率"

- 特征 13 个：`log_n_compounds`、`degree_pct`、`betweenness_pct`、`hub_robust`、
  `mr_significant`、`mr_direction`、`coloc_pp4`、`tissue_expressed`、`enrich_hit`、
  `ortholog_ok`、`dock_norm`、`n_evidence_types`、`has_L4`。
- 模型：L2 正则逻辑回归，纯 Python IRLS（无 numpy 依赖），特征标准化后拟合。
- 输出：`P(gold)` + bootstrap 95% CI（默认 300 次重抽样，存前 200 组系数）。
- **无模型时不输出概率**：返回 `UNCALIBRATED` + 级别先验排序，并在报告里写
  "未校准不得写作概率"。这是防止"0–1 评分"被误读成概率的关键。

### 1.4 实测（益智胶囊数据，CDDP 金标准 core 20 个）

| 口径 | Brier | ECE | AUC |
|---|---|---|---|
| 样本内 | 0.0945 | 0.0538 | 0.871 |
| **样本外（5 折 OOF）** | **0.1052** | **0.0543** | 0.834 |

OOF ECE 0.054 < 0.10 → 达到 DoD 第 1 条。同时模型暴露了一个反直觉系数：
`enrich_hit` 权重 **−0.733**（落在显著通路的靶点反而降分）。原因：本基准里通路成员
与"高支持成分数"高度共线，控制后者后前者无增量信息。这类结果必须写进 Limitations，
不能只贴正系数。

### 1.5 决策规则（写死在代码里，报告只引用）

| 规则 | 条件 | 含义 |
|---|---|---|
| `wetlab_candidate` | `P ≥ 0.70` 且 `≥2 类独立证据` | 可进入湿实验候选 |
| `hypothesis` | `P ≥ 0.70` 但单一证据类 | 假设生成，需补第二类 |
| `insufficient` | `P < 0.70` | 证据不足，不得写入结论 |

规则随模型落盘到 `09_evidence/decision_rules.json`；报告与门禁都读它，
避免"模板里另写一套阈值"。

---

## 2. 基准改造（`np_benchmark` + `np_split`）

### 2.1 指标口径

| 旧 | 新 | 理由 |
|---|---|---|
| AUC 为主 | **hits@K + EF@K 为主**，AUC 辅助 | 金标准覆盖候选池 42% 时，AUC 对"前 20 名有几个真靶点"不敏感 |
| Recall@K | 降为辅助 | 基率不可比，跨数据集无法直接比 |
| 单基线 | **预注册 B0/B1** | 看到结果后挑基线属探索性分析 |

实测（真实数据）：AUC 0.599→0.705，但 **hits@20 都是 6、EF@20 都是 0.71×**（低于随机）。
这条必须正面写，不能只报 AUC。

### 2.2 三类切分

| 切分 | 做法 | 解释边界 |
|---|---|---|
`random` | 阳性随机分折 | **乐观**，仅作参照 |
| `complex` | 按基因家族整组留出（无 CORUM 时用前缀代理，输出标 `proxy=gene_family`） | 对**固定排序**而言，测的是"性能是否集中在少数家族"，**不是训练泛化误差** |
| `lodo` | 留一数据库：库 A 建排序、库 B 作金标准 | 先做包含检验，>0.95 判"疑似子集"并报警 |

实测：random 0.629±0.018 vs complex 0.628±0.049（均值几乎相同、方差变大）；
命中集中度 0.39（前 20% 家族贡献 39% 命中）。LODO 检测出 **ETCM2 ⊂ TCMSP（包含率>0.95）**，
因此"跨库验证"在该数据上**不成立**——这是必须写进 Limitations 的硬事实。

### 2.3 外部工具对照

`np_split compare --gold-csv g.csv --pred A=a.csv --pred B=b.csv` 支持任意排序表，
同金标准/同候选池/同指标输出 hits@K 与 EF@K。**接口存在 ≠ 已做外部对照**：
论文里若只有本平台内部基线，不得写"与 SwissTargetPrediction 对比"。

---

## 3. 多重性（`np_fdr.py`）

- 汇总族：`MR_IVW`、`ENRICH`、`BENCHMARK_AT_K`、`NULL_MODEL`（`COLOC_PP4` 不进 BH——
  PP.H4 是后验概率不是 p 值，混入会污染口径）。
- 族内 BH + 全局 BH 双口径；`q_global` 是**唯一**可写入结论的显著性口径。
- 生成中英双语 Methods 段，含"检测了多少假设、用什么校正、校正后剩几个"。
- MR 层 `--fdr` 对同批靶点做 BH；q≥α 自动打 `FDR_NOT_SIGNIFICANT`。

---

## 4. COLOC 与 R `coloc` 对齐（`np_coloc_compare.py`）——本轮最重要的修复

### 4.1 发现的真错误

原实现：

```python
l3 = sum(logsumexp(a_i, b_i))        # H3：两个不同变异
l4 = sum(a_i + b_i)                  # H4：✗ 对所有 SNP 求和
```

H4 的定义是"**同一个**变异同时影响两个性状"，正确写法是
`l4 = logsumexp_i(a_i + b_i)`（对候选变异取 logsumexp，而不是对全部 SNP 求和）。
原写法把几十个噪声 SNP 的负 ABF 也累加，系统性压低 H4。

### 4.2 修前/修后（合成位点，共享因果 vs 两个独立变异各半）

| 版本 | Pearson r（PP.H4） | MAE | 最大偏差 | 判定 |
|---|---|---|---|---|
| 修前 | 0.735 | 0.254 | 0.989 | WARN（11/20 位点把 H4 误判为 H3） |
| **修后** | **1.000** | **0.000** | 0.000 | **PASS**（20 位点 + 40 位点两组） |

### 4.3 措辞规则（写进模块输出）

- `PASS` → Methods 可写 "colocalisation analysis using an ABF implementation matching `coloc.abf`"，
  并给出对照数字。
- `WARN/FAIL/NOT_ASSESSED` → 只能写 "an ABF-based colocalisation approximation"。

R 不可用时状态为 `R_UNAVAILABLE`，脚本 `run_coloc.R` 会写到输出目录，**不伪造结果**。

### 4.4 一个坑（值得记）

R 在中文 Windows 上按本地编码解析脚本，**非 ASCII 路径会被截断成乱码**
（实测 `D:/WorkBuddy/eh:+fd;6/...` → `cannot open the connection`）。
解决方式不是让用户改路径，而是把输入/输出/脚本都放进纯 ASCII 临时目录，跑完再搬回。

---

## 5. 可复现（`np_snapshot.py`）

| 能力 | 命令 | 实测 |
|---|---|---|
| 版本探测 | `probe --save` | STRING **12.5** / Ensembl **116** / MyGene **20261117** |
| 快照打包 | `pack -o snapshot.tar.gz` | 10.3 MB，5 文件，含 manifest（逐文件 sha256） |
| 校验 | `verify snapshot.tar.gz` | OK；篡改一个字节 → `CORRUPT` |
| 挂载 | `mount snapshot.tar.gz --to dir` | 生成 `activate_snapshot.sh/.bat`（设 `NP_CACHE`/`NP_SNAPSHOT`） |
| **离线复现** | 挂载后重跑流水线 | 3 个关键产物 SHA256 **逐字节相同** |
| manifest | 自动 | `versions: {STRING:12.5, Ensembl:116, MyGene:20261117}` |

`verify` 会列出缺失的缓存类型（如 `string`），避免"以为快照完整、实际离线仍降级"。

---

## 6. 怎么把这些写进论文

- **贡献句**：we provide a calibrated, evidence-aggregated target ranking with explicit
  uncertainty and multiplicity control, validated against `coloc.abf` and evaluated under
  family-aware splits.
- **Fig**：证据账本示意 → 校准曲线（OOF）→ random vs complex 切分 → 外部对照 → 案例。
- **Table**：门禁看板（23 门）+ 覆盖率 + 快照哈希。
- **Limitations 必写**：ETCM2⊂TCMSP（跨库验证不成立）、hits@20 与随机无差异、
  `enrich_hit` 负权重、complex 切分是家族代理而非真实复合物、MR 未在真实 GWAS 上跑通。
