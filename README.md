# network-pharmacology-pro

网络药理学全流程研究平台：**可审计、可复现、带阴性对照**。

## 这是什么

把网络药理学从"跑数据库、画网络图"升级为能过审的研究流程。核心不是更多分析，
而是**每个结论都能回答三个问题**：数据从哪来、为什么这么算、结果最多能说明什么。

## 与其他网络药理学技能的区别

| 维度 | 常见技能 | 本平台 |
|---|---|---|
| 流程 | 步骤清单 + 脚本集合 | 11 阶段，阶段间文件契约，可单步重跑 |
| 质量控制 | 靠使用者自觉 | 26 条红线审计器（`np_audit.py`），机器判定 PASS/WARN/FAIL |
| 统计 | 报 P 值 | 四套零模型 + 双背景 + 冗余折叠 + 三层分级 |
| hub 判定 | Degree 排序 | 8 指标共识 + bootstrap 入选频率 + 度匹配零模型 |
| 暴露 | 不涉及 | 饮片→成分→游离 Cmax→比值分级 |
| 对接 | 报打分 | 重对接门禁 + 分级 + 覆盖说明 |
| 报告 | Markdown | 自包含 HTML（内嵌 SVG，无字体问题）+ 自动方法学段落 |
| 溯源 | 无 | 每阶段 manifest（输入哈希/参数/种子/输出哈希） |
| 依赖 | pandas/networkx/matplotlib | **零第三方依赖**（numpy/openpyxl 可选；UI 仅需 streamlit） |
| 可用性 | 命令行 20+ 脚本 | **Streamlit 一键平台** + 真实节点进度条 + 三态水印报告 |
| 验证 | 无基准 | **hits@20/EF + 校准（Brier/ECE）为主指标**、AUC 仅作辅助；4 方剂金标准 + 6 类基准切分 + 120 项 pytest + CI 五 job（Linux/Windows/coverage 棘轮/Docker） |

## 怎么读性能数字（重要）

**不要只看 AUC。** 在本仓库的真实基准（益智胶囊：TCMSP 边 vs ETCM2 金标准，基率 42%）上：

| 口径 | B0 支持成分数 | B1 promiscuity 降权 | 说明 |
|---|---|---|---|
| AUC | 0.599 | **0.705** | 看着"提升 0.106" |
| hits@20 | 6 | 6 | **一个都没多** |
| EF@20 | 0.71× | 0.71× | 低于随机 |

结论：AUC 的提升没有转化为"前 20 名里更多真靶点"。所以我们把主指标改成 **hits@K + EF@K + 校准（Brier/ECE）**，
并把这条负面结果写进 `benchmark_report.md` 第 1 节，而不是藏起来。详见"局限"。

## 快速开始

### 方式一：网页界面（推荐，零命令行）

```bash
pip install streamlit            # 仅此一个额外依赖
streamlit run app.py             # 浏览器自动打开 http://localhost:8501
```

四个页签：**一键分析**（方剂+疾病 → 进度条 → 下载带水印报告）、**基准测试**（hits@20/EF/校准曲线）、
**因果校验**（Target_Intervention_Mismatch）、**环境自检**。
不装 streamlit 也不影响命令行与论文流程——UI 只是编排层，不重算任何指标。

### 方式二：命令行

```bash
python scripts/np_doctor.py          # 环境自检（20 秒）
python scripts/np_run.py --demo      # 内置演示数据，离线跑完 11 阶段
python evals/selftest.py             # 断言套件
python -m pytest tests -q            # 113 passed / 1 skipped（含契约/证据/FDR/切分）

# 配置驱动的 DAG（推荐用于真实项目）
python scripts/np_pipeline.py --config config/pipeline.yaml --list
python scripts/np_pipeline.py --inputs inputs.json -o np_out --enable-funnel
python scripts/np_uireport.py --dir np_out --label "复方丹参滴丸-冠心病" --verdict DRAFT

# 基准测试（证明"算得准"，不是只有图好看）
python scripts/np_benchmark.py list
python scripts/np_benchmark.py evaluate --formula CDDP --pred np_out/05_ppi/hubs.csv -o np_out/benchmark
python scripts/np_benchmark.py report --dir np_out/benchmark

# 测试
pytest tests -q
```

演示会在 `np_out/` 生成：标准化审计、成分三张清单、交集与韦恩、PPI 与 hub、
零模型、富集分级、剂量核算、对接门禁、反向优化、审计报告、HTML 主报告。

## 真实使用

```bash
python scripts/np_run.py \
  --targets compound_targets.csv \
  --compounds compounds.csv \
  --disease disease_targets.tsv \
  --edges string_edges.tsv \
  --gmt pathways.gmt \
  --label "黄芪-当归-糖尿病肾病" \
  --route oral \
  -o my_project/
```

## 目录

```
network-pharmacology-pro/
├── SKILL.md                  技能入口（路由 + 硬规则 + 决策门）
├── references/               12 份知识库
│   ├── 01-workflow-detail.md     11 阶段逐步详解
│   ├── 02-databases.md           数据库全表与坑
│   ├── 03-design-patterns.md     6 种设计模式 + 4 档工作量
│   ├── 04-statistics.md          背景/零模型/稳健性
│   ├── 05-docking-md.md          对接与 MD 判据
│   ├── 06-dose-exposure.md       剂量可达性
│   ├── 07-writing-reviewer.md    写作与 20 问应答
│   ├── 08-teaching.md            科普与教学
│   ├── 09-pitfalls.md            37 条实战踩坑
│   ├── 10-herbal-specifics.md    中药/民族药特殊性
│   ├── 11-evidence-and-claims.md 证据分级与措辞库
│   └── 12-methods-bank.tsv       方法学文献库（可核验）
├── scripts/                  13 个脚本（12 CLI + 1 核心库）
├── assets/                   报告/论文/审稿速查模板
├── evals/selftest.py         122 + 49 项断言（含参考实现对照）
├── examples/                 演示数据 + 生成器
└── data/protein_symbols.tsv  蛋白名→基因符号词典（约 320 条）
```

## 设计依据

方法学决策均可追溯到可核验文献，见 `references/12-methods-bank.tsv`。几条关键的：

- 多数"活性模块识别"算法在**度序列相同的随机网络**上表现不降（Brief Bioinform 2021, PMID 33782690）
  → hub 判定必须做度匹配零模型。
- PPI 的幂律度分布可由**研究偏倚**解释（eLife 2024, PMID 39660719）→ "度最高=最重要"不成立。
- 全 PDB 自动对接 top-1 重对接成功率仅 **30.1%**（J Chem Inf Model 2024, PMID 38108627）
  → 逐靶点重对接验证不可省。
- AlphaFold2+对接预测抗菌靶点 auROC ≈ **0.48**（Mol Syst Biol 2022, PMID 36065847）
  → 对接不能作为靶点正确性的证据。
- 多数富集方法在随机数据上 P 值偏斜（Brief Bioinform 2024, PMID 38436561）
  → 必须报背景并做零模型校准。

## 边界

- 所有产出均为**计算预测**（L1–L2），不替代实验。
- 不作临床用药建议。
- 演示数据为构造数据，**不得用于发表**。
- 内置词典约 320 条，覆盖常见 TCMSP 蛋白名；未覆盖项会进入未命中表，需人工/联网补全。

## 许可

MIT。方法学引用请同时引用原始文献（见 `12-methods-bank.tsv`），而非仅引用本技能。
