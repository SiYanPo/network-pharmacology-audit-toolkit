# 30 · 孤立靶点扩展与同源类型隔离

> 对应改动：`np_ppi.py`（STRING 一阶邻居扩展）、`np_ortholog.py`（FuncRep 隔离 + 防漂移合并）、
> `np_finalcheck.py`（G22、`--auto-fix`）、`np_uireport.py`（打印水印/正文黄框/短命令）、
> `np_interactive.py`（SIF 导出）、`np_mr.py`（`--drug-direction`）。

---

## 1 孤立靶点为什么是"结构性硬伤"

G21 只回答"交集基因是否在 PPI 节点集合里"，0.944 就算通过。但对湿实验的人来说，
真正的问题是：**ESR1 只有 2 条边、CYP2D6 一条都没有，你凭什么说它是核心靶点？**

- 对接：没有共晶/口袋信息时，孤立节点只能盲对接（已被 `np_pocket` 明确禁止）；
- 湿实验：没有互作伙伴，就无法设计 Co-IP/通路验证；
- 审稿人：会直接写"this node is topologically unsupported"。

所以要做的是**补结构基础**，而不是把阈值改宽。

### 1.1 算法（`np_ppi.py::expand_isolated`）

```
seeds = {核心靶点中 degree < min_degree} ∪ {核心靶点中完全不在网络里的}
edges = STRING /api/json/interaction_partners?identifiers=<seeds>&species=9606&required_score=400&limit=10
ext   = 原图 ∪ edges
```

三个关键决定：

1. **用 `interaction_partners` 而不是 `network`**：`network` 只返回"给定基因之间"的边，
   对孤立节点等于什么都没加；`interaction_partners` 才会把一阶邻居拉进来。
2. **不改主网络**：原图仍产出 `topology.csv` / `hubs.csv` / `network.svg`；
   扩展图单独写 `extended_ppi_edges.tsv` / `extended_ppi_nodes.csv` / `extended_ppi.svg` /
   `extended_ppi_manifest.json`。理由：用扩展图的度去重排 hub，等于自己给自己加边抬排名。
3. **离线不装死**：`--offline` 或 API 失败 → `status=OFFLINE`，主网络照常，报告写明
   "扩展子网未生成，孤立节点仍缺结构基础"，**绝不伪造边**。

实测（构造数据，8 个种子：ACE/AKR1B1/AR/ESR1/F2/PPARA + CYP2D6/NFE2L2）：

| 基因 | 原始 degree | 扩展后 degree |
|---|---|---|
| ESR1 | 2 | **12** |
| CYP2D6 | 0（不在网络） | **10** |
| NFE2L2 | 0（不在网络） | **10** |
| IL6（对照） | 15 | 15（未被改动） |

新增 69 节点 / 80 边；G21 的扩展覆盖率从 0.944 → **1.000**（176 节点），
`dataflow_check.json` 同时记录 `coverage` 与 `coverage_extended`。

---

## 2 功能替代 vs 直系同源：必须在数据结构上隔离

CXCL8 → Cxcl1 不是同源，是**功能类似物**（人 IL-8 在啮齿类无直系同源，且 CXCR1 缺失）。
把它写进"直系同源引物表"，PCR 一定失败，而且论文会被抓。

### 2.1 数据结构

`np_ortholog.py` 输出新增两列：

| 字段 | 取值 | 含义 |
|---|---|---|
| `mapping_type` | `ortholog` / `functional_analog` | 映射类型 |
| `direct_ortholog` | `yes` / `no` | 是否可用于引物设计 |
| `ortholog` | `Trp53` / `Cxcl1(FuncRep);Cxcl2(FuncRep)` | **FuncRep 后缀强制存在** |

同时生成 `func_rep_warning.md`：逐基因列出替代关系、机制说明与三条硬性约束
（禁止写成直系同源引物、必须写"功能类似物"、引物/抗体需另行验证）。

### 2.2 门禁 G22

```
若存在 status=no_direct_ortholog 或 mapping_type=functional_analog：
    要求 ortholog 字段含 (FuncRep)            否则 FAIL
    要求 mapping_type 不是 ortholog            否则 FAIL（类型污染）
    要求 func_rep_warning.md 存在              否则 FAIL
无同源结果文件 → 通过（不阻断；一旦跑了就必须标对）
```

实测：修正前（旧版 CSV，Cxcl1 无标记）G22 FAIL；重跑后
"1 个功能替代均带 (FuncRep) 标记，且警示文件齐备" → PASS。

### 2.3 防"重跑漂移"

Ensembl 限速/代理 502 会让同一批基因"这次查不到"。若直接覆盖 CSV，报告会出现
"昨天有、今天没有"，审稿人无法接受。`np_ortholog.py` 增加了 **merge_previous**：
本次未解析、历史已解析的基因沿用历史值（`source` 标 `+prev`），并在控制台明确打印
"本次未解析的 N 个基因沿用上次成功结果"。失败结果**不入缓存**，下次仍会重试。

---

## 3 打印体验：水印只在屏幕

```css
@media print{
  .wm{display:none!important}      /* 平铺水印：屏幕专属 */
  .stamp{display:none!important}   /* 右下角标：打印也去掉，避免每页重复 */
  .print-only{display:block}       /* 顶部一行"本报告为 DRAFT 版（时间）" */
  .dash{page-break-inside:auto}    /* 看板允许跨页，避免首页大片空白 */
}
```

修复命令在表格里会被压成一行短形式（`python np_defense.py temporal ...`，
输出目录统一显示为 `np_out`），完整命令放在 `title` 属性里悬停可见。
实测：PDF 10 页，全文 0 次水印，看板首页完整可读。

正文里未通过的门禁改为**黄色高亮框**（`gate_alert()`），只在对应章节出现一次，
不遮挡阅读：S3↔G3、S4/S5↔G21、S8↔G22。

---

## 4 交互图导出闭环（SIF）

顶栏新增「导出 SIF」：只导出**当前可见**的节点与边（用户先用搜索/只看核心筛选），
一次产出两个文件：

- `network_visible.sif`：`source  pp  target`（Cytoscape 标准 SIF）
- `network_visible_attributes.csv`：`id, degree, hub_rank, verdict, class, betweenness`

拖进 Cytoscape 桌面版即可继续排版（SIF 建网络，属性表按 id 合并）。
**不导出 `.cys`**：那是 Cytoscape 会话快照（含样式/视图状态），跨版本兼容性差，
SIF + 属性表是更稳的交换格式。

---

## 5 MR 方向一致性（`--drug-direction`）

```bash
python np_mr.py analyze ... --drug-direction inhibit
```

规则与 `np_intervention` 完全一致，但作用在单次 MR：

| MR β | 起效所需 | 药物方向 | 判定 |
|---|---|---|---|
| >0 | inhibit | inhibit | OK |
| >0 | inhibit | activate | SEVERE（p<0.05 且 F≥10）否则 WARN |
| <0 | activate | inhibit | **SEVERE**（"药到病除"变"雪上加霜"） |

输出 `intervention` 字段（drug_direction / required_action / consistent / severity / message），
并追加 flag `TARGET_INTERVENTION_MISMATCH`；报告第 9 节渲染红色警告框，
明确写"该靶点不得作为主结论"。

措辞已修正过一次真实 bug：早期版本对 β<0 会输出"该靶点降低使风险降低"（双重错误），
现在统一为"遗传证据提示该靶点表达**升高**与疾病风险**降低**相关（β=-0.30）"。

---

## 6 故障自愈：`np_finalcheck.py --auto-fix`

按依赖顺序执行**可自动化**的步骤，每步记录 rc/耗时/输出末行，写 `autofix_log.json`：

| 步骤 | 前提 | 命令 |
|---|---|---|
| 同源转换（在线补全） | 有 hubs.csv | `np_ortholog --input hubs.csv --col node -o 11_ortholog` |
| 孤立靶点扩展子网 | 有 PPI 边表 + 交集 | `np_ppi -e ... --core ... --expand-isolated -o 05_ppi_extended` |
| 富集（HPA 组织背景） | 有交集基因 | `np_enrich -g ... --background-mode hpa [--hpa-tissue]` |
| 重新生成报告 | 总是 | `np_uireport --dir ... --verdict DRAFT` |
| 导出 PDF | 上一步成功 | `np_pdf --html ui_report.html` |

缺输入 → SKIP 并写"缺什么"；不伪造、不覆盖主分析产物（扩展子网写到独立目录）。
实测：5 步全部 OK，35 s；随后门禁从 2/7 硬性 → 仍 2/7（G1/G2 需要真实实验数据，
本来就不在自动化范围），但 G21/G22 与报告/PDF 都已就绪。

**边界**：自动修复只解决"能算但没算"的门禁；G1（时间切分盲测）、G2（Apo/Holo）、
G8（decoy）、G13（共识打分）依赖用户自己的数据/结构，不在自动化范围——报告里也不会假装它们通过。

---

## 7 反模式（第九轮新增）

| 反模式 | 后果 | 正确做法 |
|---|---|---|
| 为让孤立节点进 hub 榜而直接加边 | 自我抬排名，审稿人一查就穿 | 扩展图独立落盘，主网络不动 |
| 把 FuncRep 写进直系同源引物表 | 实验必失败 | `(FuncRep)` 后缀 + G22 硬门禁 |
| 打印版铺满 DRAFT 水印 | 导师没法看草稿 | `@media print` 隐藏水印与角标 |
| 自动修复覆盖主分析产物 | 结果不可追溯 | 扩展/重跑写独立目录 + autofix_log |
| 网络抖动导致重跑结果变差 | 无法复现 | merge_previous 只升不降 |
| 导出 .cys 当交换格式 | 跨版本打不开 | SIF + 属性表 |
