# 29 · 数据流一致性、交互图与数据溯源

> 对应改动：`np_finalcheck.py`（G21）、`np_interactive.py`、`np_uireport.py`（看板/溯源/解读）、
> `np_pdf.py`、`run_benchmark.py`、`np_ortholog.py` v2、`np_enrich.py`（HPA 背景）、
> `np_mr_plots.py`、`np_interpret.py`。

---

## 1 为什么加 G21：报告"自相矛盾"其实不是 bug，但必须被证明不是 bug

现象（用户第一反应）：第 4 节交集表排第一的是 IL6，第 5 节 PPI hub 前 10 却没有 IL6。

真实数据核对（构造数据，107 节点 / 433 边）：

| 基因 | 交集排名（按支持成分数） | PPI hub 排名（按共识中心性） | degree |
|---|---|---|---|
| IL6 | 1（10 个成分） | 18 | 15 |
| MAPK14 | 3 | 5 | 23 |
| HIF1A | 4 | 37 | 10 |
| ESR1 | 5 | 86 | 2 |
| BCL2 | 中游 | **1** | 30 |

结论：**两节用的是两套排序键**——S4 是"药理学证据密度"，S5 是"网络拓扑位置"。
IL6 被 10 个成分支持，但在这个 STRING 子图里只有 15 条边；BCL2 只有少数成分支持，但 degree=30。
两者不可能一致，也不应该一致。

但"不是 bug"需要**机器证明**，否则读者只能信你一面之词。所以加 G21：

```python
check_dataflow(d):
    inter  = intersection_genes.txt          # 交集基因
    nodes  = ppi_edges_clean.tsv 的两端节点   # PPI 节点集合
    missing = inter - nodes
    coverage = 1 - |missing|/|inter|
    verdict = PASS(≥0.90) / WARN(0.75–0.90) / FAIL(<0.75)
```

同时输出 `dataflow_check.json`：
- `missing_genes`：缺失清单（报告会点名）
- `intersection_rank_in_ppi`：交集前 15 个基因的 hub 排名与 degree（报告 4.1 表直接引用）
- `ranking_note`：写明"S4/S5 排序键不同，名单不重叠属预期；本门禁校验成员一致性，不是排名一致性"

**为什么阈值是 0.90 而不是 1.00**：交集基因若在 STRING 中没有任何边（孤立节点），
它本就不会出现在 PPI 网络里，这不是数据错，而是网络覆盖问题。
1.00 会逼着人删基因或伪造边；0.90 允许少量孤立节点，但必须点名。

报告侧同步改了三处：
1. 标题写明排序键：`4 成分-疾病交集靶点（排序键：支持成分数）`、
   `5 PPI 拓扑与 hub（排序键：拓扑共识中心性，与第 4 节不同）`；
2. 新增 4.1 表：交集前 10 基因的 `n_compounds / hub 排名 / degree / 共识判定`；
3. 缺口超过 10% 时，看板里 G21 标红并给出修复命令（降低 `--min-score` 或换边表）。

---

## 2 交互式网络（Cytoscape.js）

静态 SVG 在 400 节点以上就是 hairball：标签互相压、放大只能截图、点不了。
`np_interactive.py` 输出**单文件 HTML**：

- 数据以 JSON 内嵌；Cytoscape.js 默认走 CDN，`--inline-js cytoscape.min.js` 可完全离线；
- 节点大小 = degree，颜色 = 核心/ Top10 hub / 背景；边宽 = 置信度；
- 顶栏：搜索（高亮 + 缩放到命中）、布局切换（cose/concentric/circle/grid/breadthfirst）、
  只看核心、导出 PNG、重置；
- 点节点：右侧面板显示 degree / hub 排名 / 共识判定 / 中介中心性 / 一阶邻居数，
  同时淡化其他节点（审稿人现场最常做的动作）；
- 报告里用 `<iframe srcdoc="...">` 内嵌（单文件可直接分发），打印时自动隐藏，改出静态 SVG。

体积参考：107 节点 / 433 边 → 51 KB；1000 节点建议 `--max-edges 2000` 控制。

---

## 3 Checklist Dashboard：让 DRAFT 不再是"无头苍蝇"

`np_uireport.dashboard_html()` 读 `finalcheck.json`，在报告最上方渲染：

```
投稿门禁看板  [DRAFT]  硬性 1/6 通过 · 软性 2/15 通过 · 共 21 门
┌────┬──────────────────┬──────┬────────┬──────────────────────────┬─────────────────────────┐
│ G1 │ 时间切分盲测      │ hard │ 未通过 │ 缺文件 temporal_valid... │ python np_defense.py ... │
│ G2 │ Apo/Holo 状态标注 │ hard │ 未通过 │ 缺文件 dock_state_...    │ python np_defense.py ... │
└────┴──────────────────┴──────┴────────┴──────────────────────────┴─────────────────────────┘
已通过 3 门（点击展开）
```

三个设计要点：
1. **只列未通过项**（通过项折叠），否则 21 行表格会淹没重点；
2. **修复命令复用 `np_finalcheck.fix_command()`**，与命令行诊断同一套逻辑，不会出现"报告说 A、命令行说 B"；
3. 没跑门禁时显示灰色提示 + 生成命令，**不假装通过**。

---

## 4 Data Provenance：溯源要精确到版本与哈希

报告第 10 节自动生成：

| 角色 | 文件 | 大小 | 文件时间 | SHA256(16) |
|---|---|---|---|---|
| compounds | compounds.csv | 2.7 KB | 2026-11-19 11:53 | 3f9c… |
| disease | disease_targets.tsv | 3.9 KB | … | … |
| gmt | pathways.gmt | … | … | … |

外加：每个阶段的 `tool_version` / `started`（来自 `run_manifest.json`）、运行环境包版本。

**必须诚实的一点**：文件 mtime **不等于**数据库检索日期。
GeneCards 的检索日期、TCMSP 版本号只能由导出环节写进 `db_status.json`；
没有就写"未记录"，不要用"最新版"含糊过去。报告里这句话是显式的。

---

## 5 PDF：放弃 weasyprint，改用系统浏览器无头模式

`np_pdf.py` 的发现顺序：
Windows（Edge → Chrome → Brave）→ macOS（Chrome.app → Edge.app）→ Linux（PATH 上 6 个名字），
也支持 `--browser` 与环境变量 `NP_BROWSER`。

```bash
msedge.exe --headless=new --disable-gpu --no-pdf-header-footer \
  --virtual-time-budget=15000 --print-to-pdf=report.pdf "file:///D:/out/report.html"
```

实测（Windows 11 + Edge 157）：235 KB HTML → **839 KB PDF**，耗时约 6 s。
校验：文件存在、>5 KB、magic 为 `%PDF`；失败时回退 `--headless`（旧版），再失败才回退 weasyprint，
最后给"用浏览器 Ctrl+P"的提示。**任何一步失败都不假装成功**。

---

## 6 一键基准：`run_benchmark.py`

```bash
python run_benchmark.py --map CDDP=cddp_edges.csv,HLJDT=hljdt_edges.csv
python run_benchmark.py --edges-dir data/formula_edges/   # 按 <ID>*.csv 自动找
python run_benchmark.py --pipeline-dir runs/              # 有 inputs_<ID>.json 时先跑流水线
python run_benchmark.py --demo                            # 构造数据，仅验证管线
```

输出：`benchmark_summary.csv/json`（含 skipped 原因）、`fig_benchmark_bars.png/svg`
（Precision@20 / Recall@20 / Recall@50 + 随机基线，300 dpi）、`fig_benchmark_roc.png/svg`、
`benchmark_report.md`。

**缺数据时的行为**：该方剂记 `SKIP` 并写明"缺成分-靶点边表"，**不猜、不插值、不借用别的方剂**。
`--demo` 会先跑 `np_standardize` 把蛋白全名换成基因符号，否则与金标准无交集（真实踩过的坑）。

---

## 7 同源转换 v2：不允许出现裸 NA

解析顺序与 `status` 字段：

| status | 含义 | ortholog 列 |
|---|---|---|
| `resolved_local` | 内置表命中（150+ 高频靶点） | 目标物种符号 |
| `no_direct_ortholog` | **已知无直系同源**（CXCL8、CYP3A4/2D6/2C9） | 功能类似物（如 Cxcl1;Cxcl2）+ note 说明机制 |
| `resolved_online` | Ensembl REST 或 MyGene HomoloGene 命中 | 符号（一对多则分号并列） |
| `no_ortholog_in_species` | 有同源记录但目标物种缺失 | 空 + 说明 |
| `unresolved` | 两条在线通路都失败 | 空 + 明确原因（不写 NA） |

性能：107 个基因从"逐基因 3 次请求（>120 s 超时）"改为
"并发取 homology（N 次）+ 一次批量 lookup（1 次）"，**约 20–25 s**；失败结果不进缓存，
避免一次网络抖动永久污染。

实测：107 个 hub → 72 本地 + 34 在线 + 1 功能类似物（CXCL8），**unresolved = 0**。

---

## 8 富集背景：HPA 组织表达替代全基因组

`np_enrich.py` 新参数：`--background-mode {auto,deg,hpa,genome}`、`--hpa-tissue`、`--hpa-min-ntpm`。

- `auto`：有 `--background` 用 DEG；否则下载/读取 HPA `rna_tissue_consensus.tsv`，
  取该组织 nTPM ≥ 1 的基因作背景；
- 拿不到 HPA（离线/代理失败）→ 回退全基因组，并**强制**在报告顶部写
  `P-value may be inflated / P 值可能膨胀` 的红框；
- 未指定组织时，HPA 背景退化为"任一组组表达"（接近全基因组），也会打警告。

实测（heart）：N = 13,708（vs 全基因组 20,596），显著通路 10 → 折叠后 7。
`enrich_meta.json` 记录 `background_mode / background_source / N / pvalue_inflation_warning`，
报告第 7 节直接读它渲染徽章。

---

## 9 MR 严谨化与出版级图

- `--bfile` + `--plink`：调用 `plink --clump` 做**真 LD clumping**；
  PLINK 缺失 → `PLINK_MISSING` + 安装指引，不静默降级；
- 无 LD 面板时结果多一个 `exploratory_only: true` 与 `EXPLORATORY_ONLY_NO_LD_PANEL` flag，
  `note` 明确写"仅作探索性参考，不得作为投稿的因果证据"；
- `np_mr_plots.py`：森林图（逐 SNP Wald 比 + IVW/Egger/WM 汇总 + OR 文本列）、
  散点图（βy vs βx + 三方法拟合）、LocusZoom（上下双轨、共享 x 轴）；
- `np_interpret.py`：把 β/OR/p/F/Egger/PP.H4/协同分翻译成中英双语标准句，
  COLOC 0.5–0.8 明确写"中等证据，需排除多效性"，协同标注竞争性拮抗时**禁止**写协同。

---

## 10 反模式清单（第八轮新增）

| 反模式 | 后果 | 正确做法 |
|---|---|---|
| 两节排序键不同却都不写 | 读者怀疑数据造假 | 标题写排序键 + 4.1 对齐表 + G21 机器校验 |
| 同源查不到就写 NA | 湿实验直接卡死 | 在线兜底 + 功能类似物 + status 说明 |
| 全基因组背景不加警告 | 审稿人直接拒 | HPA/DEG 背景 + 红框警告 |
| DRAFT 水印但不列原因 | 用户无从下手 | Checklist Dashboard + 可复制命令 |
| 报告不写数据版本 | 无法复现 | SHA256 + 文件时间 + 阶段版本（缺则写未记录） |
| 用 weasyprint 当唯一 PDF 方案 | Windows 装不上 | 系统浏览器无头打印，weasyprint 仅兜底 |
