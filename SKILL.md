---
name: network-pharmacology-pro
display_name: 网络药理学研究平台
display_name_en: Network Pharmacology Research Platform
version: 3.0.19
category: research-methodology
license: MIT
description_zh: 网络药理学全流程研究平台（可审计、可复现、带阴性对照）。覆盖中药/复方/民族药机制、单成分解析、老药新用、网络毒理学、Q-marker、反向网络药理学与组方优化；内置 12 阶段流水线（含真实对接 08a）、29 条红线审计器、四套**可自校准**的零模型阴性对照、三种可追溯的分子门控（OB/DL·类药性·证据）、hub 稳健性检验（bootstrap + 度匹配零模型）、富集三层分级与冗余折叠、剂量可达性核算、对接门禁、自包含 HTML 报告。零第三方依赖，离线可跑。
description: 网络药理学全流程研究平台——把"成分-靶点-通路-疾病"从数据库搬运升级为可审计、可复现、带阴性对照的研究流程。当用户提到网络药理学、network pharmacology、中药/复方/民族药机制、成分-靶点网络、TCMSP、OB/DL 筛选、靶点交集、韦恩图、STRING、PPI、hub 基因、CytoHubba、GO/KEGG 富集、分子对接、AutoDock Vina、结合能、网络毒理学、Q-marker、质量标志物、老药新用、药物重定位、反向网络药理学、组方优化、审稿意见、方法学质疑、剂量可达性、入血成分、阴性对照、零模型、稳健性、可复现性、报告规范、网络药理学评价方法指南、专家共识，或需要判断"这个网络药理学结论靠不靠得住"时使用。内置 12 阶段流水线（含真实对接 08a）、29 条红线审计（结构/提示分离）、四套可自校准零模型、hub 稳健性检验、富集三层分级、剂量可达性核算与自包含 HTML 报告；零依赖、离线可跑、每步留溯源。
description_en: Full-process, auditable network-pharmacology research platform with negative controls, hub-stability testing, tiered enrichment and dose-reachability audit. 12-stage pipeline (incl. real docking), 29-rule auditor, three traceable molecular gates, zero mandatory third-party dependencies.
triggers:
  - 网络药理学
  - network pharmacology
  - 成分-靶点-通路
  - TCMSP
  - OB/DL 筛选
  - PPI hub
  - 富集分析
  - 分子对接
  - 网络毒理学
  - Q-marker
  - 老药新用
  - 反向网络药理学
  - 审稿意见 网络药理学
---

# 网络药理学研究平台

> **一句话定位**：把网络药理学从"跑数据库、画网络图"升级为**可审计、可复现、能过审**的研究流程。
> 核心判断标准不是"图好不好看"，而是：**换个人、换台机器，能不能得到同样的结论；以及这个结论最多能说到几分。**

## 0. 三条不可让渡的原则

1. **网络药理学是假设生成器，不是结论机器。** 任何"成分-靶点-通路"结果的最强表述是"提示/支持"，除非有 L3（体外）或 L4（体内）证据。
2. **没有零模型的"紧密互联"不算证据。** STRING 网络天然偏向研究热点蛋白（eLife 2024, PMID 39660719）；交集基因当然连得紧。不做阴性对照，"核心模块"随时可被一句"这只是度大"击穿。
3. **可复现是方法学的一部分，不是附录。** 阈值、版本、日期、随机种子、脚本命令必须同时存在；缺一项，结论不可核验。

## 0.4 v1.3.0 的关键变更（23 篇复现审计驱动的工程化升级）

> 数据来自 23 篇真实论文复现日志的聚类归因：成分漏检 52%、靶点源差异 48%、原文参数未披露 35%。

| 变更 | 触发证据 | 新工具 |
|---|---|---|
| **多源数据库适配层**（补 TCMSP 盲区：动物药/矿物药/糖苷） | 鹿角胶、牡蛎、蕨类、crocin(MW977) TCMSP 查无或仅 1 靶点 | `np_fetch.py`（HERB/ETCM2/SymMap/UniProt/ChEMBL/TDC 实测状态+注册网址+手工导入+交叉验证） |
| **复现基准集 + 置信度自评 + 归因审计** | 23 篇里 10 篇置信度<0.6，需定位链路断裂点 | `np_bench.py` + `data/benchmark_23.json` |
| **阈值自适应推断** | 4 篇因原文自定义阈值/去重复算不出 | `np_bench.py infer`（从原文文本抽 OB/DL/STRING/GeneCards/p 值） |
| **数据库版本快照与漂移对比** | 二妙散 8 成分消失、HLJDD JUN 消失需归因版本而非算法 | `np_prov.py snapshot/diff/stamp` |
| **标准中间产物五张表** | 复现不可追溯的核心痛点 | `np_export.py all/check` |
| **异常值自动告警** | 成分-靶点对偏差 78%（celastrol）无自动提示 | `np_guard.py run/demo`（>50% 偏差、PPI<70%、功效不足 BLOCK） |
| **适用域声明** | crocin MW977 大极性糖苷被 TCMSP 系统性漏检 | `np_domain.py`（新颖骨架/糖苷/矿物 → 低置信度标注） |
| **不可复现学术解释生成** | 用户需要可直接写入 Discussion 的客观表述 | `np_bench.py explain` |
| **ML 型 ADMET 预测 + 外部验证** | 只有 OB/DL 过滤、无现代 ADMET 评估；需外部 AUC 证据防御审稿人 | `np_admet.py`（TDC 官方基准 23 端点，scaffold split 外部集 + BBBP 跨库验证） |

> **降级策略透明化**：TCMSP 查无此药时不再静默返回空值，而是写 `herbs.csv matched=no` + `fallback_strategy`、
> 触发 guard 告警、强制适用域声明，并给出 HERB/ETCM2 手工导入路径（详见 `14-databases-and-provenance.md` §2）。

## 0.5 v1.1 的关键变更（相对 v1.0，均为已确认缺陷的修复）

| 变更 | 为什么必须改 |
|---|---|
| **方法学段落不再用默认值兜底** | v1.0 会用默认 OB/DL 覆盖用户实际参数（学术诚信级缺陷）。现在参数必须有来源阶段，缺失写「未记录」 |
| **MCC 改为 Bron–Kerbosch 精确实现** | v1.0 的贪心团在稠密小图上退化为常数，等于把 8 指标共识降为 7 指标 |
| **C 组零模型改为按观察集度分布精确匹配** | v1.0 复现的是总体度分布，系统性偏稀疏 → 反保守（实测：零均值 46→120，效应从 3.1× 缩到 1.19×） |
| **零模型必须自校准** | 新增平均\|Δdeg\|、精确匹配率、KS 距离；不达标直接判 C 组结论无效 |
| **审计器结构/提示分离 + 否定检测** | v1.0 对「本文不将 P 视为作用强度」误报 → 告警疲劳 |
| **结论账本 `claims.json`** | 结论↔数值↔来源文件↔参数阶段，全部可追溯 |
| **环境快照 + 复现脚本 + 复现校验器** | 换机可复现从"声称"变为"可验证"（实测 36 个产物逐字节一致） |
| **修复 PYTHONHASHSEED 导致的非确定性** | 同一数据两次运行 hubs/topology 不同（排序并列时随集合迭代顺序漂移） |

## 0.6 v1.2 的关键变更（真实 48 化合物全流程跑出来的问题）

> 本节全部来自一次真实运行：48 个跨 16 类的化合物 × 2 型糖尿病（ChEMBL + Open Targets +
> STRING + KEGG/Reactome/GO），不是构造数据。发现的问题比单元测试能暴露的更根本。

| 变更 | 触发它的现象 | 为什么要改 |
|---|---|---|
| **新增 `--gate evidence`** | 用类药性（QED/Veber）门控时，**阿托伐他汀、二甲双胍、地高辛、哇巴因被判「不合格」** | 类药性描述的是理化性质空间，不是药效。把已上市药物判为不合格，说明判据与问题错位。证据模式以「有实测活性/入血证据」纳入，类药性只作标记 |
| **类药性剔除已上市药物时强制告警** | 同上 | 静默剔除=选择性偏倚。报告单列「门控与有效性的错位」并给出两条修正路径 |
| **门控按行留痕 `gate_used` / `gate_detail`** | 混合队列（中药成分+已上市药物）无法用同一把尺子 | OB/DL 与 QED/Veber 不可比；逐行记录判据与数值，才能复算与分层比较 |
| **新增 08a `np_dockprep.py`：真实对接** | 旧版只「评估」别人给的分数，无法自证协议 | 结构准备 → Vina → **共晶配体回对接 RMSD** → 打分表，全部可复算；RMSD>2 Å 的靶点自动拒绝出分 |
| **统一裁决：B > C > A** | 真实数据里 A 显著（P=0.002）、B 不显著（P=0.44）、C 显著（P=0.004） | 三套对照不一致时，作者会挑有利的那套。现在输出唯一裁决：**B 不显著即不得称「核心模块」** |
| **度校正改为等频分箱** | 真实 PPI 上 `bc_excess_z` 全为 0 | 「度±1」窄窗里常只有 1–2 个同伴，标准差无定义 → z 恒为 0，掩盖了真正的桥接节点（改后 CLK1 度=2 却 z=4.6） |
| **`bc_excess_z` 区分 NA 与 0** | 同上 | 「算不出来」写成 0，读者会理解为「恰好等于同度均值」 |
| **词典换成 UniProt 全量** | ChEMBL 蛋白名映射率仅 87.5%，未命中含 `[NAD(+)]` 类名称 | 129k 条名称（20k 蛋白），映射率→ 可核查；并修正了推荐名截断 bug |
| **全程 ASCII 临时目录** | 中文路径下 RDKit 抛 `OSError: Bad input file` | RDKit/meeko/Vina 用 C++ 流打开文件，非 ASCII 路径直接失败且报错不提编码。中间产物走 ASCII scratch，结果再拷回 |
| **经验 P 的下限写进报告** | 200 次置换下多个条目并列在 P=0.00498 | 那是 `1/(B+1)` 的地板，不是「同样显著」；报告提示调大 B |
| **λ 的删失校正** | 同一数据 B=200 得 λ=0.16、B=1000 得 λ=0.78 | 触地板的经验 P 被右删失，直接取中位数会把「置换次数不够」误读成「解析 P 偏乐观」 |
| **门控描述贯穿全链路** | 报告的方法学段落、结论账本 C1 仍在写「OB≥30」，而实际用类药性门 | 与 v1.0 同类的「伪造方法」错误，已按 gate 动态生成，并加回归 BUG-15 |

## 0.7 v1.2.1 的关键变更（遗留决策点裁决，仅补证据链、不改契约）

- **背景规模敏感性有了真实数字**：同一批 67 基因、同一三库，`--universe` 从 20596→8000，
  `Pathways in cancer`（命中 20/531 不变）P 从 2.2e-16→5.9e-9，显著条目 481→303。
  结论：默认背景下的极小 P 有一部分是"背景设得大"带来的，无实测背景时必须标探索性（详见 `04-statistics.md` §1.1.1）。
- **黄芪甲苷（CID 13943297）换思路解析**：MW 785、QED 0.162，证实即便名称命中也必被类药性门剔除，
  故"未命中"不改变主分析结论；并记下 PubChem `CanonicalSMILES` 返回空的坑（`09-pitfalls.md` 第 55 条）。
- **方法库增补 M29–M31**（Crossref 已核验）：打分排序能力有限、成对比较优于绝对打分、稳健可复现建议。
- **三个学术分歧的裁决**（对接不排序、B 组边缘显著附敏感性、上市药双轨门控）已写入
  `CHANGELOG.md` v1.2.1 与工作区最终分析报告，`SKILL.md` 的硬规则/门禁不变。

## 0.8 v2.7 的关键变更（软件化：Benchmark / Web UI / CI）

从"一套能跑的脚本"变成"可发表、可被引用的软件"，补的是三类以前完全缺失的东西：

- **基准测试**：`np_benchmark.py` + `data/benchmark_goldstandard.json`（4 个经典方剂：复方丹参滴丸 /
  黄连解毒汤 / 连花清瘟 / 血府逐瘀汤，core 15–20 个强证据靶点）。指标含 AUC（Mann-Whitney + bootstrap CI）、
  AUPRC、Recall@K、Precision@K、EF@K、超几何 p；新增**边级跨库一致性**（敏感性/精确率/Jaccard/过预测倍数），
  用于"金标准是候选池子集"时避免靶点级 AUC 被污染。输出 `benchmark_report.md`（含可直接改写的 Methods 段与 Limitations）。
  实测：益智胶囊数据上，朴素基线 AUC 0.599 → promiscuity 降权后 0.705；对 ETCM2 边级敏感性 1.000、过预测 4.10×。
- **一键 UI**：`app.py`（Streamlit，四页签：一键分析 / 基准测试 / 因果校验 / 环境自检）。进度条由
  `np_pipeline` 的节点事件驱动（非假动画）；报告由 `np_uireport.py` 生成，带 `READY_FOR_SUBMISSION` /
  `DRAFT` / `UNVERIFIED` 平铺水印；无 weasyprint 时只出 HTML 并明确提示，不静默失败。
- **因果严谨化**：`np_mr.py` v2 增加 ① 等位基因统一（`harmonise`）② **Steiger Filtering**（r²x>r²y，
  缺样本量时返回 `STEIGER_SKIPPED` 而非假装跑过）③ **不依赖 PLINK 的 LD 去冗余**（真基因型 r² / 外部 r² 矩阵 /
  距离剪枝三档，`distance_pruning` 会自曝"不是真 clumping"）④ **COLOC 式共定位**（ABF 近似，PP.H0–H4）
  ⑤ leave-one-out 单点驱动检测。
- **方向一致性**：`np_intervention.py` 把 MR 方向与成分作用方向对齐——MR β>0（靶点升高→风险升高）要求**抑制**；
  相反即抛 `Target_Intervention_Mismatch`（Grade A → SEVERE，Grade B → WARN，Grade C 不判定），并生成
  可直接写进 Discussion 的英文段落。
- **协同升级**：`np_synergy.py --mode mechanistic`：同通路上下游（BFS 层级）+ 方向兼容才判 `Mechanistic_Synergy`；
  **共享靶点强制 `Competitive_Antagonism`**；一抑一激判 `Directional_Antagonism`；无层级证据只算 `Topological_Only`。
- **LINCS 本地化**：`np_lincs_local.py` 直读 GEO 的 `.gctx`（HDF5），只读签名基因行、逐行 z-score、
  按 perturbagen 取 top-3 聚合 → 彻底摆脱 CLUE API（下载与裁剪见 `references/28-lincs-local.md`）。
- **CI**：`tests/`（47 用例）+ `.github/workflows/test.yml` 三层（unit 矩阵 / 11 个模块 selftest / examples 全流程集成）。

对应细节：`references/27-software-benchmark-and-ui.md`、`references/28-lincs-local.md`。

## 0.9 v2.8 的关键变更（数据流一致性 / 交互图 / 看板与溯源）

- **G21 数据流一致性门禁**（`np_finalcheck.py`）：强制校验"交集基因 ⊆ PPI 节点"，
  coverage = |交集 ∩ PPI| / |交集|，≥0.90 PASS、0.75–0.90 WARN、<0.75 FAIL；
  输出 `dataflow_check.json`（缺失清单 + 交集基因的 hub 排名 + 排序键说明）。
  直接回应"S4 交集第一名不是 S5 hub 第一名"的质疑：**两节排序键不同**（支持成分数 vs 拓扑共识），
  报告标题写明键、4.1 表并排展示、G21 只校验成员一致性，不要求排名一致。
- **交互式网络**（`np_interactive.py`）：Cytoscape.js 单文件 HTML，缩放/拖拽/搜索/点节点看属性/
  只看核心/导出 PNG；报告用 iframe 内嵌，打印时自动改用静态 SVG。`--inline-js` 可完全离线。
- **Checklist Dashboard**：报告顶部列出未通过的每一门（级别/原因/**可复制修复命令**），
  命令复用 `np_finalcheck.fix_command()`，与命令行诊断同源；未跑门禁则显式提示。
- **Data Provenance**：输入文件 SHA256 + 文件时间 + 阶段 tool_version + 环境包版本；
  明确声明"文件时间 ≠ 数据库检索日期"，缺元数据写"未记录"。
- **PDF 零依赖方案**（`np_pdf.py`）：调用系统 Edge/Chrome `--headless=new --print-to-pdf`，
  跨平台自动发现；失败才回退 weasyprint。实测 235 KB HTML → 839 KB PDF。
- **一键基准**（`run_benchmark.py`）：4 方剂批量评测 → `benchmark_summary.csv` +
  Precision/Recall 分组柱状图 + 多曲线 ROC（300 dpi）+ `benchmark_report.md`；
  缺边表的方剂记 SKIP 并说明原因，绝不插值。
- **同源转换 v2**（`np_ortholog.py`）：本地表 150+ → Ensembl REST → MyGene HomoloGene；
  已知无直系同源者（CXCL8/CYP3A4 等）给**功能类似物 + 机制说明**，禁止裸 NA；
  批量两阶段查询把 107 基因从 >120 s 降到 ~20 s；实测 unresolved = 0。
- **HPA 富集背景**（`np_enrich.py`）：未提供 DEGs 时自动取 HPA 组织表达基因
  （`--hpa-tissue`，nTPM≥1）作背景；拿不到则回退全基因组并**强制**打 `P-value may be inflated`。
- **MR 出版级图与解读**：`np_mr_plots.py`（森林图/散点/LocusZoom）、
  `np_interpret.py`（MR OR、COLOC PP.H4、协同分的中英双语标准措辞）；
  `np_mr.py` 支持 `--bfile/--plink` 真 clumping，无 LD 面板时标 `exploratory_only`。

对应细节：`references/29-dataflow-interactive-provenance.md`。

## 0.10 v2.9 的关键变更（孤立靶点扩展 / FuncRep 隔离 / 打印体验 / 自愈）

- **孤立靶点扩展**（`np_ppi.py --core --expand-isolated`）：对 degree<3 及完全缺失的交集靶点，
  调 STRING `interaction_partners` 补一阶邻居，生成独立的 `extended_ppi_*.{tsv,csv,svg,json}`；
  **主网络不变**（避免用扩展度自我抬排名），离线时标 OFFLINE 不伪造。实测 ESR1 2→12、
  CYP2D6 0→10、NFE2L2 0→10，G21 扩展覆盖率 1.000。
- **FuncRep 隔离**（`np_ortholog.py`）：新增 `mapping_type`（ortholog/functional_analog）与
  `direct_ortholog`（yes/no），功能替代强制加 `(FuncRep)` 后缀并生成 `func_rep_warning.md`；
  **G22 硬门禁**校验"有 FuncRep 必须有标记 + 必须有警示文件"。
  另加 merge_previous：重跑时未解析的基因沿用上次成功结果，避免网络抖动导致结果漂移。
- **打印体验**：`@media print` 隐藏平铺水印与右下角标（屏幕仍显示），顶部改为一行 print-only 说明；
  看板允许跨页、修复命令压成单行短形式（完整命令在 title 悬停）；正文对未通过门禁用黄色框提示。
- **交互图 SIF 导出**：顶栏「导出 SIF」只导当前可见子网（SIF + 节点属性 CSV），可直接进 Cytoscape 桌面版。
- **MR 方向一致性**：`np_mr.py --drug-direction inhibit|activate`，与 IVW β 比对；
  矛盾且 p<0.05、F≥10 → `TARGET_INTERVENTION_MISMATCH` 红灯，报告渲染红色警告框。
- **故障自愈**：`np_finalcheck.py --auto-fix [--tissue liver] [--dry-run]`，按序跑同源/扩展/富集/
  报告/PDF，写 `autofix_log.json`；缺输入 SKIP 并说明，不覆盖主分析产物。

对应细节：`references/30-isolated-nodes-and-funcrep.md`。

## 0.11 v2.10 的关键变更（自审 / 可移植性 / 统计口径）

- **自审工具**（`np_selfaudit.py`）：一键产出"模块清单 / 死代码 / 自检矩阵 / 契约核对 / 依赖矩阵"。
  首轮实测：**72 模块、40 进 DAG（56%）、15 被 evals 覆盖（21%）、33 有 selftest（46%）、24 个执行图不可达**。
  这不是"功能多"，而是"必须被管理的复杂度"——报告已据此给出模块注册表与退役策略。
- **修掉 12 个模块的 `--selftest` 入口缺陷**（此前 argparse 直接 exit 2：benchloop/consensus/defense/
  evidence/goldfetch/hpa/liana/ligand_top/lincs/primer/singlecell/wetlab），现统一等价并全部 PASS。
- **核心统计模块补自检**：`np_enrich`（超几何 P 与手算一致、BH q≥p、背景过滤）、
  `np_intersect`（交集精确、分层字段）、`np_ppi`（自环/重复/阈值/共识排序）——共 15 项断言。
- **统计口径修正**：COLOC 的"中等证据"阈值由 PP.H4≥0.5 上调至 **0.6**（文献常用：>0.8 强、0.6–0.8 中，
  PMID 40399559 / 39985383），同步改单测与文档。
- **可移植性**：新增 `np_preflight.py`（Python/数据/依赖/工具/缓存/路径/CLI 七类体检）、
  `requirements-lock.txt`、`PORTABILITY.md`；`npcore.cache_dir` 统一跨运行缓存（`$NP_CACHE`），
  新增 `NP_OFFLINE=1` 全局离线开关。实测冷启动 153/423 s → 缓存后 **28 s**，产物哈希一致。
- **环境补齐**：安装 **PLINK 1.9**（官方二进制，实测 clump 6→2 真跑通）与 **acpype**（pip）；
  scanpy 在 Python 3.13 实测装不上（生态原因），已在 preflight 与文档标注为"全功能需 3.10–3.12"。
- **使用规范**：新增 `references/31-usage-and-promotion.md`（四类用户入口、SOP、参数改动条件、十种误用、工作坊大纲、推广策略）。

对应细节：`references/31-usage-and-promotion.md`、`PORTABILITY.md`、`audit/skill_audit.md`。

## 0.11 v2.11 的同行评审补强（证据账本 / 全局 FDR / 切分 / 快照）

本轮不是再加功能，而是把 v2.9 评审报告要求的“证据模型”落成可运行、可测试、可注册的层：

- **证据账本**：`np_evidence_ledger.py` 从产物汇总 `gene × evidence_type × L0–L4`，输出校准
  `P(gold)` + bootstrap 95% CI + 决策规则；OOF AUC<0.60 或系数全 0 时**拒绝给概率**
  （真实 24 靶点/10 正例实测 OOF AUC=0.471 → `WEAK_CALIBRATION`，不再包装成概率）。
- **全局 FDR**：`np_fdr.py` 做族内+全局 BH，生成 Methods 可直接引用的“多重性声明”；
  `np_mr.py --fdr` 对同批靶点做 BH，q<0.05 才允许写“支持性因果证据”。
- **切分与快照**：`np_split.py` 提供 random/complex-family/LODO 三种切分，复合物缺失时
  明示 `proxy=gene_family`；`np_snapshot.py` 记录 release/etag/sha256 并打包 `NP_SNAPSHOT`。
- **契约与注册表**：`contracts/standard_tables.md` 与 `np_export.STANDARD` 逐列对齐
  （修复 `standard_*.csv` vs 实际 `herbs.csv` 的契约漂移）；`config/module_registry.yaml`
  为 78 个模块登记 core/audit/manual/deprecated，`np_selfaudit` 不再把已登记手动工具误报死代码。
- **测试**：新增 `tests/test_evidence_ledger.py`、`test_fdr`/`test_contracts`/`test_phase1_modules`，
  pytest 从 70 → 113 用例；`np_evidence_ledger`/`np_fdr`/`np_split`/`np_snapshot` 均有 selftest。

对应细节：`references/32-evidence-and-multiplicity.md`（新增）、`audit/REVIEW_v2.9_TASK_CHECKLIST.md`。
边界与措辞（COLOC 近似、HPA 组织、zebrafish、MR 方向、示例数据）见
`references/33-limitations-and-wording.md`；外部验证/湿实验/用户研究协议见
`references/34-validation-roadmap.md`；论文图表与 7 条 DoD 见
`references/35-paper-and-dod.md`；外部求证（ETCM2/TCMSP、STRING、HPA、Delphi）见
`references/36-external-evidence-audit.md`。

## 0.12 v3.0 的关键变更（证据模型 / 基准改造 / 多重性 / 可复现）

这一版不是"再加功能"，而是把**结论的可信度**补上——回应同行评审的七寸。

- **单一排序函数**（`np_evidence_ledger.py`）：证据账本（gene × 证据类型，L0–L4 + provenance）→
  纯 Python 逻辑回归校准 → `P(gold)` ± 95% CI → 显式决策（P≥0.70 且 ≥2 类证据 = 湿实验候选）。
  **没有校准器时拒绝输出概率**（返回 `UNCALIBRATED`）。实测 5 折 OOF：Brier 0.105、ECE **0.054**、AUC 0.834。
- **基准改造**（`np_split.py` + `np_benchmark` 改口径）：主指标从 AUC 改为 **hits@K + EF@K**；
  三类切分（随机 / 家族留组 / 留一数据库）；预注册 B0/B1 基线；外部工具对照接口。
  真实结果：AUC 0.599→0.705，但 **hits@20 都是 6、EF@20=0.71×（低于随机）**——这条负面结果已写进报告。
- **多重性**（`np_fdr.py`）：族内 + 全局 BH，输出 `q_global`；MR 层 `--fdr` 对同批靶点做 BH，
  q≥α 打 `FDR_NOT_SIGNIFICANT`。报告新增"多重性声明"段（检测数/校正法/校正后剩几个）。
- **COLOC 对齐**（`np_coloc_compare.py`）：与 R `coloc::coloc.abf` 逐位点对照。
  **据此修掉一个真 bug**：H4 应为 `logsumexp(ABF₁+ABF₂)`（同一变异），原实现写成对所有 SNP 求和，
  会把共享因果误判成两个独立变异（修前 r=0.74/MAE=0.25；修后 **r=1.000/MAE=0.000**，60 个模拟位点）。
- **可复现**（`np_snapshot.py`）：缓存 `*.meta.json`（url/release/etag/sha256）+ 快照打包/挂载/校验；
  实测 STRING **12.5** / Ensembl **116** / MyGene **20261117**；挂载快照后重跑，3 个关键产物 **SHA256 逐字节相同**。
- **工程质量**：`config/module_registry.yaml`（80 模块：entry/core/audit/manual/deprecated）、
  `np_contracts.py`（24 条产物契约，覆盖 18 个声明产物）、覆盖率**棘轮**（基线 15.22%，只允许上调）、
  CI 五 job（unit / windows / coverage-ratchet / selftest / integration / docker）。
- **可移植**：`np_tools.py` 统一解析外部二进制（显式→`NP_VINA`…→`tools.local.json`→PATH→旧根），
  `D:/MDtools` 只允许出现在这一处；lint 测试拦截回退。
- **网络语义**（`np_ppi.py`）：扩展子网按实验/共表达**分层计数**、报告 hub 排名 Spearman 漂移（实测 ρ=0.862）、
  图标题硬写"degree 不参与 hub 判定"。实测 CYP2D6 的 10 个新邻居 **0 个有实验证据**——这就是"研究热度"的样子。

对应细节：`references/37-evidence-model-and-benchmark-reform.md`、`references/38-module-registry-and-contracts.md`。

## 0.13 DoD 收尾工具（tools/，v3.0 补充）

`scripts/` 之外有四个**收尾专用**工具，专治"接口就绪但没实测"：

| 工具 | 作用 | 命令 |
|---|---|---|
| `tools/stp_prepare.py` | 成分名 → SMILES（PubChem，带缓存与失败清单） | `--compounds c.csv --top 5 -o out/stp` |
| `tools/stp_convert.py` | SwissTargetPrediction 导出 CSV → `gene,score`（max/mean 聚合，可限候选池） | `--csv-dir out/stp/raw -o out/stp --restrict-to genes.txt` |
| `tools/stp_submit.py` | STP 批量提交（协议实测；STP 恢复后用；含限流退避） | `--input stp_input.txt --names a,b,c -o raw/` |
| `tools/ssim_predict.py` | **外部靶点推断管线**（SwissSimilarity ECFP4 筛 ChEMBL actives → ChEMBL/UniProt 人源靶点；STP 停摆时的对照来源，2026-11-22 实测 5 成分全通） | `--input stp_input.txt --names a,b -o out/ --top-hits 15` |
| `tools/dod_status.py` | DoD 就绪度自检（只认磁盘证据；合成数据不算外部对照） | `--dir np_out [--snapshot s.tar.gz]` |
| `tools/check_coverage.py` | 覆盖率棘轮（只允许上调） | `--coverage coverage.json --baseline coverage_baseline.txt` |

配套手册：`DoD收尾_三步操作手册.md`（真实 MR / STP 对照 / Docker 构建的材料与命令）。

## 1. 先判定任务类型（路由）

| 用户意图 | 走哪条路 | 首要用具 |
|---|---|---|
| "某方/某药治某病的机制" | 全流程 01→11 | `np_run.py` |
| "帮我筛成分/靶点/做富集/算 PPI" | 单阶段精细控制 | 对应脚本 |
| "这个结果靠不靠得住 / 审稿人会怎么问" | 审计 + 零模型 | `np_audit.py`、`np_null.py` |
| "只有疾病基因，想知道用哪些药" | 反向网络药理学 | `np_optimize.py` |
| "某成分体外有效，体内可能吗" | 剂量可达性 | `np_dose.py` |
| "对接结果怎么判 / 能不能说结合" | 对接门禁 | `np_dock.py` |
| "帮我把这些化合物真的对接一遍" | 真实对接（准备+Vina+重对接校验） | `np_dockprep.py` |
| "我的化合物没有 OB/DL 怎么办" | 门控选择 | `--gate evidence` 或 `--gate druglike`，见 §2.1 |
| "课题设计 / 标书技术路线" | 设计模式 + 工作量分级 | `references/03-design-patterns.md` |
| "论文怎么写 / 图表怎么排" | 写作规范 | `references/07-writing-reviewer.md` |
| "审稿意见怎么回" | 逐条应答 | `references/07-writing-reviewer.md` §应答 |
| "科普 / 教学" | 叙事框架 | `references/08-teaching.md` |
| "环境跑不起来" | 自检 | `np_doctor.py` |

> 多意图叠加时：先跑主流程拿到数据，再用写作/应答模块产出文本。

## 2. 标准流水线（11 阶段，阶段间只通过文件通信）

```
01 标准化 → 02 成分筛选 → 03 交集 → 04 PPI → 05 零模型
   → 06 富集 → 07 剂量 → 08 对接 → 09 反向 → 10 审计 → 11 报告
```

| 阶段 | 做什么 | 关键产出 | 硬门禁 |
|---|---|---|---|
| 01 | 蛋白名 → HGNC 符号（四层匹配） | `targets_standardized.csv` | 映射率 <80% 必须先补别名 |
| 02 | 门控筛选（OB/DL·类药性·证据三选一）+ 途径一致性 + 泛成分/PAINS 标记 | `compounds_kept/exception/rejected/unassessed.csv` | 注射/外用却用 OB 筛 → **阻断**；类药性剔除已上市药物 → **告警** |
| 03 | N 组交集 + 成分级超几何 + 证据分层 | `intersection.csv`、`venn.svg` | 交集 <5 → 功效不足，改路线 |
| 04 | 拓扑 + 8 指标共识 hub + 模块 + 稳健性 + 度分箱校正（bc_excess_z，区分 NA/0） | `hubs.csv`、`network.svg` | 无稳健 hub → 禁用「核心枢纽」 |
| 05 | 四套零模型 + **自校准** + **统一裁决（B>C>A）** | `null_tests.csv`、`null_ecdf.svg` | B 不显著 → 禁用「核心模块」 |
| 06 | 富集 + 折叠 + 三层分级 + **经验校准**（可选） | `enrich_significant.csv` | 背景过小 → 阻断 |
| 07 | 饮片→成分→游离 Cmax→比值 | `dose_result.json` | 比值 >10 → 机制结论降级 |
| 08a | **真实对接**：结构准备→Vina→共晶回对接 RMSD | `dock/docking_scores.csv`、`docking_validation.csv` | RMSD >2 Å → 该靶点不出分 |
| 08 | 对接打分分级 + 重对接门禁 | `docking_graded.csv` | RMSD >2 Å → 该靶点结果作废 |
| 09 | 覆盖度/互补性/最小覆盖组合 | `reverse_ranking.csv` | 协同增量 0 → 不得称"协同" |
| 10 | 29 条红线审计（结构 18 + 提示 11，后者不计分） | `audit/audit_report.md` | FAIL 必须修 |
| 11 | 自包含 HTML + **参数溯源表 + 结论账本** | `np_report.html`、`claims.json` | 参数「未记录」不得投稿 |

**跑法**

```bash
python scripts/np_doctor.py                 # 先自检（20 秒）
python scripts/np_run.py --demo             # 内置演示数据，离线跑完 01–11
python scripts/np_run.py --targets t.csv --compounds c.csv --disease d.tsv \
    --edges s.tsv --gmt p.gmt --label "黄芪-当归-糖尿病肾病" -o out/
python scripts/np_run.py --list                  # 看阶段
python scripts/np_run.py ... --stages 04,06      # 只跑指定阶段
python scripts/np_run.py ... --gate evidence     # 已上市药物/单体：不做类药性门槛
python scripts/np_run.py ... --dockplan plan.tsv --vina /path/vina   # 真跑对接（08a）

# 单独跑真实对接（先校验协议，再跨对接）
python scripts/np_dockprep.py --plan plan.tsv --outdir dock/ --validate-only
python scripts/np_dockprep.py --plan plan.tsv --ligands compounds.csv --outdir dock/ --reuse-validation
python scripts/np_run.py ... --from 04           # 从 04 跑到尾
python scripts/np_run.py ... --dry-run           # 只校验配置
```

> 阶段隔离的好处：改阈值重跑 04 不必重抓 TCMSP；换 GMT 重跑 06 不必重算 PPI。

## 3. 证据等级（决定措辞，不可越级）

| 等级 | 内容 | 允许措辞 |
|---|---|---|
| **L0** | 数据库注释 / 文献事实 | "已报道"（须带版本与日期） |
| **L1** | 计算预测（靶点、网络、富集） | "提示 / 支持 / 可能与…相关" |
| **L2** | 结构证据（对接、MD） | "预测可结合 / 结合模式提示" |
| **L3** | 体外 / 离体实验 | "表明" |
| **L4** | 体内 / 临床 | "证明"（须给样本量与效应量） |

**常见越级**：
- ❌ "网络药理学**证明**该方通过 PI3K-Akt 通路发挥作用" → ✅ "**提示**该方可能通过 PI3K-Akt 通路发挥作用（L1）"
- ❌ "分子对接**验证**了成分与靶点的结合" → ✅ "分子对接**预测**二者可结合（L2），尚需 SPR/ITC 验证"
- ❌ "hub 基因 AKT1 是**核心靶点**" → ✅ "AKT1 在 200 次 bootstrap 中入选频率 0.98，为**拓扑稳健的 hub 节点**"

## 4. 六个决策门（触发即停，不许自动推进）

| 门 | 触发条件 | 必须做的事 |
|---|---|---|
| **G1 途径不一致** | 给药途径为注射/外用，却用 OB/DL 筛 | 改用含量/入血成分，或声明"可成药性预筛"并降级 |
| **G2 映射失败** | 标准化映射率 <80% 或出现 not_found | 补别名表重跑；**不得**把非人类蛋白强行映射 |
| **G3 零/低交集** | 交集 =0 或 <5 | 按 `intersect_report.md` 的六步恢复序列排查；不得换阈值硬凑 |
| **G4 泛成分主导** | 达标成分中泛成分 >50% | 单列泛成分 + 用非泛成分做敏感性分析 |
| **G5 零模型不通过** | B 组 P≥0.05 或 C 组不显著 | 机制表述降级为"靶点集合"，不得称"功能模块" |
| **G6 暴露不可达** | 体外/体内游离浓度比 >10 | 机制结论降级，或补低浓度实验 |

## 5. 硬规则（违反即返工）

1. **阈值必须写全**：OB、DL、STRING combined score、富集 q、对接阈值，一次列清并给依据。
2. **每个数据库写版本 + 访问日期**（`TCMSP, 2026-11-19 访问`）。
3. **靶点必须标准化**到 HGNC 符号；非人类/噪音条目单列剔除。
4. **疾病靶点必须分来源强度**：OMIM/GWAS（遗传）> DisGeNET/CTD（审编）> GeneCards（整合）> GEO（特定样本表达）。
5. **富集必须报背景 + 多重校正 + 冗余折叠**；伞状通路（Pathways in cancer）不得作机制结论。
6. **中药必须写炮制与基原**；生品与制品是两个受试物。
7. **复方 ≠ 单药相加**：至少讨论溶出、代谢、配伍对成分的影响。
8. **对接必须报口袋定义、质子化、重对接 RMSD**；只对接部分成分须说明选择依据。
9. **剂量链必须闭合**：饮片克数 → 成分量 → 游离浓度 → 与体外浓度比较。
10. **称 hub 必须给入选频率**；<0.60 只能写"度排序前列节点"。
11. **P 值是显著性，不是效应值**；GeneRatio 也不是活性评分。
12. **不得声称通路"激活/抑制"**，除非有磷酸化/表达方向实验。
13. **不得编造**含量、药典条目、文献、批号；不确定标"需核实"。
14. **结论必须带证据等级**；无 L3+ 证据永远停在"提示/支持"。
15. **报告必须有局限性章节**（数据库偏倚、阈值局限、无实验、暴露不可达）。

## 6. 与"常见水刊写法"的对照（自查用）

| 水刊写法 | 问题 | 本平台的处理 |
|---|---|---|
| "TCMSP 筛出 20 个成分、300 个靶点" | 不说来源/阈值/给药途径 | 02 强制写阈值 + 途径一致性门禁 |
| "Degree 最大者为核心靶点" | 度≠重要性 | 04 多指标共识 + bootstrap + 度匹配零模型 |
| "交集基因 PPI 网络连接紧密" | 无对照 | 05 四套零模型 |
| "KEGG 富集到 Pathways in cancer" | 伞状通路当机制 | 06 三层分级自动标记 |
| "对接结合能 −8.5，结合良好" | 预测当实测 | 08 门禁 + 措辞模板 |
| "围绕 AKT1 一个基因做富集" | 富集对象错 | 06 只接受基因集输入，单基因直接报错 |
| "结果表明…通路被激活" | 无方向信息 | 硬规则 12 |
| 只报好看的 3 条通路 | 选择性报告 | 06 输出完整表 + 折叠前后计数 |

## 7. 脚本清单

| 脚本 | 功能 | 典型命令 |
|---|---|---|
| `np_doctor.py` | 环境/数据源自检 + 能力矩阵 | `python np_doctor.py` |
| `np_standardize.py` | 靶点标准化（四层匹配 + 物种列/名称双路判定非人类） | `-i targets.csv -o out/ [--online]` |
| `np_compounds.py` | 成分筛选（三种门控 + 途径 + 泛成分 + PAINS） | `-i compounds.csv -o out/ --gate auto\|evidence` |
| `np_intersect.py` | N 组交集 + 成分级富集 + 证据分层 | `--drug t.csv --disease d.tsv -o out/` |
| `np_ppi.py` | 拓扑 + hub 共识 + 模块 + 稳健性 | `-e edges.tsv -o out/ --top 15 --bootstrap 200` |
| `np_null.py` | 四套零模型 | `-g genes.txt -e edges.tsv -d d.tsv -o out/ -B 1000` |
| `np_enrich.py` | 富集（GMT/g:Profiler）+ 折叠 + 分级 | `-g genes.txt -m p.gmt -o out/` |
| `np_dose.py` | 剂量可达性 | `--daily 9 --content 0.042 --mw 302.24 --invitro 50` |
| `np_dock.py` | 对接分级 + 重对接门禁（读 RMSD 列） | `-i docking.csv -o out/` |
| `np_dockprep.py` | **真实对接：结构准备 → Vina → 共晶回对接自校验** | `--plan plan.tsv --ligands c.csv --outdir dock/` |
| `np_optimize.py` | 反向网络药理学 + 组合优化 | `--drug t.csv --disease d.tsv -o out/` |
| `np_audit.py` | 29 条红线审计（结构 18 计分 + 提示 11） | `--dir out/ [--strict]` |
| `np_report.py` | 自包含 HTML 报告 | `--dir out/ --label "…"` |
| `np_run.py` | 12 阶段编排（01–11 + 08a） | `--demo` / `--targets … -o out/` |
| `npcore.py` | 核心库（统计/图/SVG/溯源/账本） | 被以上脚本导入 |
| `evals/selftest.py` | 算法与 IO 自检（124 项） | `python evals/selftest.py` |
| `evals/regression.py` | **缺陷回归 + 契约测试（63 项）** | `python evals/regression.py` |
| `evals/repro_check.py` | **跨机复现校验（区分内容差异与时间戳差异）** | `--a out1 --b out2` |
| `np_fetch.py` *(v1.3)* | **多源数据库适配**（HERB/ETCM2/UniProt/ChEMBL/TDC；可达性探测+手工导入+交叉验证） | `probe --live` / `cross --source a=x.tsv --min-sources 2` / `import` |
| `np_prov.py` *(v1.3)* | **数据库版本快照 + 漂移对比 + SHA256 清单** | `snapshot --outdir out/` / `diff a.json b.json` / `stamp out/` |
| `np_export.py` *(v1.3)* | **标准化中间产物**（herbs/compounds/targets_mapped/PPI.sif/enrichment 五张表） | `all --dir out/` / `check --dir out/` |
| `np_guard.py` *(v1.3)* | **异常值告警**（数量偏差>50%、PPI<70%、功效不足自动 BLOCK） | `run --local l.json --paper p.json` / `demo` |
| `np_bench.py` *(v1.3)* | **23 篇基准集**：置信度评分/归因聚类/阈值推断/学术解释生成 | `score` / `audit` / `infer --input paper.txt` / `explain --below .6` |
| `np_domain.py` *(v1.3)* | **适用域声明**（MW/骨架/糖苷/矿物 → 低置信度标注） | `--compounds "quercetin,crocin" --pubchem` |
| `np_admet.py` *(v1.4)* | **ML 型 ADMET 预测 + 外部验证**（TDC 官方基准 23 端点，官方 scaffold split，RF×5seeds，双层外部验证） | `download/extract/train/external/predict/validate`（详见 `15-admet-ml.md`） |
| `np_omics.py` *(v1.5)* | **多组学交叉映射**（GEO series matrix→DEGs→∩NP 靶点；官方 annot.gz 自动探针注释） | `scan/deg/cross --gse GSE14520 --targets x.txt` |
| `np_survival.py` *(v1.5)* | **靶点临床相关性**（cBioPortal TCGA 生存：KM-SVG+log-rank+HR） | `studies --kw lihc` / `survival --study LIHC --genes STAT3,MMP9` |
| `np_sens.py` *(v1.5)* | **参数敏感性分析**（pchembl/OB/DL 网格扫描→稳健基因+Jaccard 告警） | `sweep --compounds x.csv --outgrid "6:0,5:0"` |
| `np_tcmviz.py` *(v1.5)* | **性味归经可视化**（四气雷达/桑基/君臣佐使网络 SVG，内置 40 味药库） | `radar/sankey --herbs 黄芪,当归` / `junchen --map 黄芪:君,...` |
| `np_cyto.py` *(v1.5)* | **Cytoscape 工程文件**（XGMML+SIF+成分-靶点-通路超网络 SVG） | `export --edges x.tsv --paths enrich.csv` |
| `np_lab.py` *(v1.5)* | **湿实验方案生成器**（UniProt 功能+引物/抗体权威检索链接+20 病种细胞/动物模型库） | `plan --genes STAT3,PTGS2 --disease 胃癌` |
| `np_toxpath.py` *(v1.5)* | **网络毒理学**（ADMET 风险→毒理通路映射，Safety Discussion 模板） | `run --profile admet_profile.csv --threshold 0.7` |
| `np_qa.py` *(v1.5)* | **一键 QA 报告**（产物完整性+SHA256+原文对比+审稿薄弱环节） | `report --dir out/ --paper paper.json` |
| `np_md.py` *(v1.8)* | **分子动力学模拟**（OpenMM + **GROMACS@WSL 已解锁** + **MM/PBSA**；prep/run/complex/gmx/mmpbsa/analyze/bench） | `check` / `prep` / `run` / `gmx --ps 20000` / `mmpbsa --tpr --traj --index`（详见 `17-md-simulation.md`） |
| `np_figs.py` *(v2.0)* | **投稿级图片自动生成**（气泡/柱图/热图/PPI/PRISMA；PDF+SVG+600dpi PNG+TIFF RGB；期刊尺寸预设） | `bubble/bar/heatmap/ppi/prisma --outdir figs`（详见 `20-publication-figures.md`） |
| `np_figcheck.py` *(v2.0)* | **图片自检**（11 项/图 · C01–C11：分辨率/字号/重叠/越界/字体嵌入/色彩/溯源/格式齐备 → 可投稿判定；`--list-checks` 可枚举） | `--dir figs --out figure_qa.md` |
| `np_db.py` *(v2.0)* | **数据库登记与体检**（44 个库/引擎：可达性、是否需注册/机构 IP） | `list` / `check --json db_status.json` / `guide` |
| `np_corpus.py` *(v1.7)* | **500 篇级语料流水线**（IF 代理筛选/PMC 全文收割/方法参数挖掘/疾病靶点映射/归类统计） | `screen/harvest/extract/disease/classify`（详见 `18-scale-reproduction.md`） |
| `np_repro_batch.py` *(v1.7)* | **批量复现引擎**（TCMSP 重算→与原文对比→差异归因；交集复现+蛋白名→基因映射） | `run --methods m.csv --herb-cache DIR` / `inter --disease-json` / `stats` |

| `np_benchmark.py` *(v2.7)* | **金标准基准**（内置 4 方剂 core/extended；AUC/AUPRC/Recall@K/EF + 边级一致性 + benchmark_report.md） | `list` / `evaluate --formula CDDP --pred r.csv` / `concordance --pred-edges a.csv --gold-edges b.csv` / `report --dir out/` |
| `np_evidence_ledger.py` *(v3.1)* | **证据账本 → 校准 P(gold) → 决策**（L0–L4；OOF Brier/ECE/校准曲线/DCA；OOF AUC<0.60 或系数全 0 拒绝给概率） | `build --dir out/` / `fit --ledger ... --gold gold.txt -o out/09_evidence` / `score --dir out/` / `score_targets()` |
| `np_fdr.py` *(v3.0)* | **全局多重性控制**（族内+全局 BH；q_global 是唯一可写结论口径；自动生成 Methods 声明） | `run --dir out/` / `mr --dir out/` / `statement --dir out/` |
| `np_split.py` *(v3.0)* | **基准切分**（random / complex 家族 / LODO 留一库；复合物代理声明；外部工具对照） | `run --edges pred.csv --gold-csv gold.csv --split complex -o out/` / `lodo --db-edges all.csv` |
| `np_snapshot.py` *(v3.0)* | **版本化缓存快照**（release/etag/sha256；pack/mount/verify；`NP_SNAPSHOT` 离线复现） | `probe` / `pack -o snapshot.tar.gz` / `mount ... --to DIR` / `verify ...` |
| `np_uireport.py` *(v2.7)* | **一键交付报告**（11 章 + 三态水印 + 产物 SHA256 清单；缺产物写"未运行"） | `--dir out/ --label "…" --verdict READY\|DRAFT\|UNVERIFIED [--pdf]` |
| `np_intervention.py` *(v2.7)* | **靶点-干预方向一致性**（MR β 方向 vs 成分 action → Target_Intervention_Mismatch） | `check --mr-dir mr/ --action action.csv -o out/`（SEVERE 退出码 3） |
| `np_lincs_local.py` *(v2.7)* | **LINCS L1000 gctx 本地流式比对**（只读签名行，不载入 10 GB） | `probe --gctx X.gctx` / `query --gctx X.gctx --signature sig.json` |
| `app.py` *(v2.7)* | **Streamlit 一键平台**（4 页签；真实节点进度条；下载 HTML/PDF/ZIP） | `streamlit run app.py` |
| `np_mr.py` *(v2.7)* | 孟德尔随机化 **v2**：Steiger + LD 去冗余（基因型/矩阵/距离）+ COLOC + LOO | `analyze --exposure e.tsv --outcome o.tsv --n-exposure 500 --n-outcome 450000 --genotypes p.tped` |
| `np_synergy.py` *(v2.7)* | 协同 **v2**：机制定向（通路上下游+方向兼容）／共享靶点判竞争性拮抗 | `--mode mechanistic --ppi e.tsv --pathways p.csv --actions a.csv` |

| `np_interactive.py` *(v2.8)* | **交互式 PPI（Cytoscape.js 单文件）**（缩放/拖拽/搜索/节点属性/导出PNG；可内联 JS 离线） | `--edges ppi.tsv --hubs hubs.csv --core genes.txt -o net.html [--inline-js cytoscape.min.js]` |
| `np_pdf.py` *(v2.8)* | **HTML→PDF（系统浏览器无头模式，零 Python 依赖）** | `--html report.html [--pdf out.pdf]` / `--check` |
| `run_benchmark.py` *(v2.8)* | **4 方剂一键基准 + 论文级柱状图/ROC** | `--map CDDP=a.csv,...` / `--edges-dir DIR` / `--demo` |
| `np_mr_plots.py` *(v2.8)* | **MR/COLOC 出版级图**（森林图/散点/LocusZoom） | `all --mr mr_result.json -o figs` / `locus --exposure e.tsv --outcome o.tsv` |
| `np_interpret.py` *(v2.8)* | **结果解读引擎**（中英双语标准措辞，防"显著/不显著"乱写） | `demo`；或 `from np_interpret import interpret_mr` |

| `np_ppi.py` *(v2.9)* | 新增 `--core/--expand-isolated`：孤立靶点补 STRING 一阶邻居，独立扩展子网 | `-e edges.tsv --core genes.txt --expand-isolated -o out/` |
| `np_ortholog.py` *(v2.9)* | **FuncRep 隔离**：mapping_type/direct_ortholog + `(FuncRep)` 后缀 + 警示文件 + 防漂移合并 | `--input hubs.csv --col node -o 11_ortholog` |
| `np_finalcheck.py` *(v2.9)* | 新增 **G22**（同源类型隔离）与 `--auto-fix`（按序自愈 + autofix_log.json） | `--dir out/ --auto-fix --tissue liver` |

| `np_selfaudit.py` *(v2.10)* | **Skill 自审**（模块清单/死代码/自检矩阵/契约）/ 依赖矩阵 | `run --outdir np_out --run-selftests` |
| `np_preflight.py` *(v2.10)* | **换机体检**（Python/数据/依赖/工具/缓存/路径/CLI 七类） | `--json preflight.json [--strict]` |

| `np_evidence_ledger.py` *(v3.0)* | **证据账本 + 校准打分**（唯一排序函数；无模型时拒绝给概率） | `auto --dir np_out` / `fit --ledger L --gold g.txt -o out/` / `score --dir np_out` |
| `np_fdr.py` *(v3.0)* | **全局多重性控制**（族内+全局 BH；生成 Methods 声明） | `run --dir np_out` / `mr --dir np_out` / `statement --dir np_out` |
| `np_split.py` *(v3.0)* | **基准切分**（随机/家族留组/留一库）+ 外部工具对照 | `run --edges e.csv --gold-csv g.csv --split complex` / `lodo --db-edges all.csv` / `compare --gold-csv g.csv --pred A=a.csv` |
| `np_snapshot.py` *(v3.0)* | **版本化快照**（meta/打包/挂载/校验；离线复现） | `probe --save` / `pack -o snap.tar.gz` / `mount snap.tar.gz --to dir` / `verify` |
| `np_contracts.py` *(v3.0)* | **产物契约**（列名/类型/空表；改列名立刻红） | `check --dir np_out` / `coverage` |
| `np_tools.py` *(v3.0)* | **外部工具解析唯一入口**（五级；给修复提示） | `list` / `resolve vina` |
| `np_coloc_compare.py` *(v3.0)* | **ABF vs R coloc 一致性**（PASS 才可写 colocalisation） | `run --loci loci.tsv --r-lib ~/Rlibs -o out/` |

**所有脚本**：`--help` 可用；中文路径与空格路径可用；BOM/CRLF/GB18030 自动识别；输出 CSV 为 `utf-8-sig`（Excel 不乱码）+ SVG（不依赖中文字体）；异常以中文提示退出，不抛原始 traceback。
**自检**：`python evals/selftest.py`（124 项）+ `python evals/regression.py`（63 项缺陷回归与契约，覆盖 v1.2 的 8 组新缺陷）。
**复现校验**：`python evals/repro_check.py --a 旧产物 --b 新产物`。

## 8. 现实预期（先看数字，再谈结论）

这些数字决定了"能说到几分"，写作时应内化：

| 事实 | 数值 | 含义 |
|---|---|---|
| 自动对接流程 top-1 重对接成功率（全 PDB） | **30.1%**（结构过滤后 61.8%） | 对接≠对得准，逐靶点验证不可省 |
| AlphaFold2 + 对接预测抗菌靶点的 auROC | **0.48**（重打分后 0.63） | 前瞻性对接靶点预测接近随机 |
| 多数"活性模块识别"算法在度匹配随机网络上的表现 | 不下降 | 它们主要学到的是度，不是生物学 |
| 2024 年网络药理学相关论文年产出 | 数千篇量级 | 同质化严重，方法学差异才是区分度 |

> 结论不是"别做"，而是：**把力气花在对照、稳健性与暴露核算上，而不是把图做得更漂亮。**

## 9. 资源索引（按需加载）

**references/**
- `01-workflow-detail.md` —— 11 阶段逐步操作、判据、失败恢复
- `02-databases.md` —— 数据库全表（用途/访问要点/坑/替代），含 TCMSP token、UniProt 匹配陷阱
- `03-design-patterns.md` —— 六种研究设计模式（A–F）+ 四档工作量（Lite/Standard/Advanced/Publication+）
- `04-statistics.md` —— 超几何/背景/多重校正/零模型/稳健性，含公式与选值依据
- `05-docking-md.md` —— 对接与 MD 的判据、几何标准、常见错误
- `06-dose-exposure.md` —— 饮片→成分→浓度完整换算与算例
- `07-writing-reviewer.md` —— 论文结构、图件清单、审稿 20 问逐条应答
- `08-teaching.md` —— 科普/教学叙事与误区清单
- `09-pitfalls.md` —— 37 条实战踩坑（来自社区与文献），含复现命令
- `10-herbal-specifics.md` —— 中药/民族药特殊性：炮制、基原、配伍、Q-marker
- `11-evidence-and-claims.md` —— 证据分级细则与措辞库
- `12-methods-bank.tsv` —— 方法学文献库（每条可核验：PMID/DOI/用途）
- `13-usage-and-promotion.md` —— **高专业度用户使用规范、教学方案与推广策略**
- `14-databases-and-provenance.md` —— **(v1.3) 外部数据库扩展/注册网址/降级策略/版本快照/23 篇基准集**
- `15-admet-ml.md` —— **(v1.4) ML 型 ADMET 预测与外部验证**（TDC 官方基准接入、指标约定、跨库验证语义）
- `16-advanced-modules.md` —— **(v1.5) 高级分析模块**（GEO 交叉/TCGA 生存/敏感性/Cytoscape/湿实验指导/网络毒理学/QA 报告/MD 边界）
- `17-md-simulation.md` —— **(v1.6) MD 模拟**（本机环境长期记忆、双后端选择、20ns 流程参数、生产建议）
- `18-scale-reproduction.md` —— **(v1.7) 500 篇规模化复现**（语料构建/相关性硬筛/方法挖掘/归因统计/CNKI 交接）
- `19-500paper-audit-roadmap.md` —— **(v1.8) 500 篇高级手段审计 + V2.0 路线图**（方法频次/图表风格/三档分级）
- `20-publication-figures.md` —— **(v2.0) 投稿级图片自动化与自检**（尺寸预设/11 项机检 · 图数×项数统计/实测全部 PASS）

- `27-software-benchmark-and-ui.md` —— **(v2.7) 软件化三件套**（基准设计与指标口径/金标准构建规则/UI 架构/测试分层/审稿 6 问）
- `28-lincs-local.md` —— **(v2.7) LINCS L1000 本地化**（GEO 下载与断点续传/裁剪脚本/gctx 结构/打分公式与官方差异/常见报错）

- `29-dataflow-interactive-provenance.md` —— **(v2.8) 数据流一致性/交互图/看板与溯源**（G21 设计、Cytoscape 集成、Dashboard、PDF 无头、基准脚本、同源兜底、HPA 背景、MR 出图与解读）

- `30-isolated-nodes-and-funcrep.md` —— **(v2.9) 孤立靶点扩展与同源类型隔离**（STRING 一阶邻居、FuncRep 门禁、打印水印、SIF 导出、MR 方向、自愈）

- `31-usage-and-promotion.md` —— **(v2.10) 高专业度用户规范与推广**（四类入口/SOP/参数改动条件/十种误用/工作坊/支持流程）

- `37-evidence-model-and-benchmark-reform.md` —— **(v3.0) 证据模型与基准改造**（账本字段/校准口径/切分语义/多重性声明/快照规范/COLOC 对齐）
- `38-module-registry-and-contracts.md` —— **(v3.0) 模块注册表与产物契约**（status 定义/manual 准入/契约写法/覆盖率棘轮）

**assets/**
- `report_outline.md` —— 研究报告/结题报告骨架
- `paper_outline.md` —— 论文提纲与图件清单
- `reviewer_qa.md` —— 审稿质疑速查卡（可打印）

**evals/** —— `selftest.py`（121 项断言）
**examples/** —— 演示数据（**构造数据，仅跑通流程，不得用于发表**）+ `make_demo.py` 生成器
**data/** —— `protein_symbols.tsv` + `protein_symbols_uniprot.tsv.gz` 内置词典；
`benchmark_23.json` **23 篇复现基准集**（原文参数/结果/复现结果/归因标签/置信度）；
`fetch_cache/source_registry.json` 数据源实测注册表

## 10. 数据治理与隐私

- **未发表/受控数据**（临床样本、患者来源测序、企业批号数据）**不得**喂给公共 API
  （STRING / g:Profiler / UniProt / PubChem）。本平台所有在线调用只发送基因符号或化合物名，
  但若基因列表本身即敏感信息（如未发表队列的特征基因），仍应先脱敏或改用离线模式。
- 默认离线：`np_enrich.py -m x.gmt`、`np_ppi.py -e 本地边表`、`np_standardize.py`（不带 `--online`）。
- 所有输出文件仅写入用户指定目录，不上传、不缓存到第三方。
- 报告中的患者/样本标识须去除；病例编号用内部代号。

## 11. 执行风格

- **先结论后展开**：第一段给判断（靠得住/需降级/不成立），再给依据。
- **给数字**：`q=3.1e-11，命中 8/14 个基因`，不给"显著富集"四个字了事。
- **给下一步**：做什么、补什么数据、怎么算做完。
- **标注不确定性**：区分"数据库记载""文献报道""我的推断"。
- **不迎合**：如果用户的方案有 G1–G6 门禁问题，直接指出并给替代路径——这是本平台存在的意义。

## 14. 模块注册表：手动工具与已弃用模块（v3.0）

> 这张表回答"你到底用了哪些模块"。`core` 见 `config/pipeline.yaml` 的 13 个 stage，
> `audit` 见 `audit_tools`（33 个）。下表是**需要人工判断何时用**的工具；
> 每条都有最小命令，不允许"存在但没人知道怎么调"。

### 14.1 手动工具（30 个）

| 模块 | 何时用 | 最小命令 |
|---|---|---|
| `np_admet` | ADMET/类药性过滤（Lipinski、肠吸收、PAINS） | `python scripts/np_admet.py -i compounds.csv -o out/` |
| `np_audit` | 项目红线审计：对目录结构与断言做结构性检查 | `python scripts/np_audit.py --dir np_out` |
| `np_bench` | 文献复现基准集管理与归因审计（≠ np_benchmark 的性能评测） | `python scripts/np_bench.py --help` |
| `np_benchloop` | np_bench 的 golden 回归闭环（capture/verify/accept） | `python scripts/np_benchloop.py capture --dir np_out` |
| `np_corpus` | 文献语料库构建与检索（复现研究用） | `python scripts/np_corpus.py --help` |
| `np_cyto` | 细胞因子/炎性介质分类注释 | `python scripts/np_cyto.py -i genes.txt -o out/` |
| `np_db` | 本地数据库缓存管理（TCMSP/HERB 离线副本） | `python scripts/np_db.py status` |
| `np_dock` | 分子对接批量执行（AutoDock Vina） | `python scripts/np_dock.py --receptor r.pdbqt --ligands lig/ -o dock/` |
| `np_dockprep` | 受体/配体准备（去水、加氢、转 pdbqt） | `python scripts/np_dockprep.py --help` |
| `np_domain` | 蛋白结构域注释（Pfam/InterPro 映射） | `python scripts/np_domain.py -i genes.txt -o out/` |
| `np_dose` | 剂量可达性核算：饮片→摄入→游离暴露 vs 体外浓度 | `python scripts/np_dose.py --dose dose.csv -o out/` |
| `np_export` | 中间产物标准化导出（交付/归档用） | `python scripts/np_export.py --dir np_out -o export/` |
| `np_fetch` | 靶点批量抓取（多库合并，带来源标记） | `python scripts/np_fetch.py --genes genes.txt -o out/` |
| `np_figcheck` | 图件规范检查（分辨率/字体/尺寸/色盲友好） | `python scripts/np_figcheck.py --dir figs/` |
| `np_figs` | 出版级图件批量生成（300 dpi PNG + SVG） | `python scripts/np_figs.py --dir np_out -o figs/` |
| `np_guard` | 结论护栏：拦截「计算预测」被写成「疗效证明」 | `python scripts/np_guard.py --md draft.md` |
| `np_herbcache` | 药材名缓存与别名回填（减少重复在线查询） | `python scripts/np_herbcache.py --help` |
| `np_htmlreport` | 稿件式 HTML 报告（偏论文排版，与 np_uireport 交付版不同） | `python scripts/np_htmlreport.py --dir np_out --label X --out manuscript.html` |
| `np_lab` | 湿实验方案生成（分组、剂量、读数、样本量） | `python scripts/np_lab.py --targets hubs.csv -o lab/` |
| `np_md` | 分子动力学前后处理（GROMACS；需 NP_GMX） | `python scripts/np_md.py --help` |
| `np_omics` | 组学数据整合（转录组/蛋白组/代谢组交集） | `python scripts/np_omics.py --help` |
| `np_optimize` | 参数寻优（阈值/权重扫描，输出敏感性面） | `python scripts/np_optimize.py --help` |
| `np_qa` | 结果问答与解释（对报告做定向追问） | `python scripts/np_qa.py --dir np_out --ask 'why IL6'` |
| `np_repro_batch` | 批量复现引擎（500 篇语料 → 重算 → 差异归因） | `python scripts/np_repro_batch.py --help` |
| `np_sens` | 参数敏感性扫描（OB/DL/阈值） | `python scripts/np_sens.py --dir np_out -o sens/` |
| `np_survival` | 生存分析（TCGA 表达 + 临床，KM/Cox） | `python scripts/np_survival.py --help` |
| `np_tcmviz` | 药材-成分-靶点可视化（中医语义分层图） | `python scripts/np_tcmviz.py --dir np_out -o tcm.html` |
| `np_toxpath` | 毒理通路映射（CTD 通路证据） | `python scripts/np_toxpath.py -i genes.txt -o out/` |
| `np_tools` | 外部工具解析的唯一入口（路径/来源/修复提示）；换机时先跑它 | `python scripts/np_tools.py list` |
| `np_coloc_compare` | ABF 共定位 vs R coloc 一致性核验；PASS 才可在 Methods 写 colocalisation | `python scripts/np_coloc_compare.py run --loci loci.tsv --r-lib ~/Rlibs -o out/` |

### 14.2 已弃用（2 个，保留兼容）

| 模块 | 替代者 / 原因 |
|---|---|
| `np_doctor` | 被 np_preflight（七类体检）+ np_envcheck（工具路径）取代 |
| `np_prov` | 被 np_snapshot 取代（版本探测 + 快照 + 离线挂载）；保留仅为兼容旧脚本 |

> 规则：`deprecated` 不得被 `pipeline.yaml` 引用（`tests/test_contracts.py` 会查），
> 也不得出现在新文档的推荐路径里。
