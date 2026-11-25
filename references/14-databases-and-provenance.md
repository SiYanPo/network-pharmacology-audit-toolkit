# 14. 外部数据库扩展、溯源与复现基准（v1.3.0）

> 这一章解决 23 篇复现审计暴露的最大短板：**成分漏检（52%）与靶点源差异（48%）**。
> 全部网址与状态为 2026-11-20 真实探测结果，不编造 API。

## 1. 数据源注册表（实测状态）

| 源 | 实测状态 | 覆盖 | 注册/下载 | 在本 skill 中的角色 |
|---|---|---|---|---|
| TCMSP | api_session（动态 token） | 植物药为主；矿物/动物/糖苷差 | 免注册 | 成分+靶点主源 |
| HERB | spa_manual（HTTP 200） | 7263 味/49258 成分/12933 靶点，**含动物药+矿物药** | http://herb.ac.cn/#/Download/ | 补非植物药；SPA 后端 `/chedi/api/` 需抓 func_name，优先手工下载→`np_fetch.py import` |
| ETCM2 | spa_manual（HTTP 200） | **动物药 44/矿物药 30**，二维配体相似性靶点预测 | http://www.tcmip.cn/ETCM2/front/#/ | 补动物/矿物药靶点；浏览器交互 |
| TCMID | timeout | 含动物/矿物药，无靶点 | 免注册 | 仅作成分名同义词来源 |
| NPASS | ssl_fail | 天然产物+动物药+活性 | 免注册 | 浏览器访问后手工导入 |
| TCM Database@Taiwan | timeout | 矿物/动物药，典籍汇编 | http://tcm.cmu.edu.tw/ | 手工下载 |
| SymMap | web_ok | 证候-中药-成分-靶点，含动物药 | http://www.symmap.org/ | 手工导出 |
| UniProt REST | api_ok | ID→HGNC 符号 | 免注册 | **强制预处理** |
| ChEMBL | api_flaky（偶发 500） | 实测活性靶点 | 免注册 | 实测靶点核对（需指数退避） |
| STRING | api_ok | PPI | 免注册 | POST 长列表 |
| Open Targets | api_ok | 疾病靶点 | 免注册 | GraphQL page:{index,size} |
| DrugBank | **auth_required（403）** | 上市药/ADMET 金标准 | https://go.drugbank.com/pricing | 外部验证（需机构邮箱注册） |
| TDC | web_ok | ADMET/靶点基准 | https://tdcommons.ai/ | 外部验证基准 |
| BATMAN-TCM2 | api_flaky（500） | 中药靶点预测 | 免注册 | 交叉验证可选源 |

## 2. 降级策略（TCMSP 查无此药时）——算法黑盒透明化

当前 v1.3.0 的确定行为，**按顺序**：

1. **同义词消歧**：`np_fetch.py` 已知别名表（半夏→Arum Ternatum、刺蒺藜→Tribulifructus、黄连→Coptidis Rhizoma）。
   - 注意：别名匹配必须精确，历史上 "Coptis" 匹配不到 "Coptidis"，误命中胡黄连。
2. **多记录选择**：herb_all_name 返回多条时，按拉丁名优先词匹配，不再盲取第一条。
3. **查无此药**（如鹿角胶、牡蛎、红景天、蕨类）：
   - **不是返回空值后静默继续**，而是写 `herbs.csv` 的 `matched=no` + `fallback_strategy`；
   - 自动触发 `np_guard.py` 的「成分数偏差」告警；
   - 报告强制写适用域声明（矿物/动物药不在 TCMSP 覆盖域）。
4. **人工兜底**：`np_fetch.py import` 把 HERB/ETCM2 下载表导入为统一 targets 格式；
   HERB 批量下载：`http://herb.ac.cn/#/Download/`（需浏览器点选六个组件分别下载）。

## 3. 多源交叉验证（--cross-check）

命令：
```bash
python np_fetch.py cross --input compounds.csv \
  --source tcmsp=tcmsp_targets.tsv --source herb=herb_targets.tsv \
  --min-sources 2 --output cross_targets.csv
```
23 篇样本模拟结论（见升级报告）：单源→≥2 库交叉，靶点数平均下降约 35–60%；
对预测源主导的论文（celastrol、alliin、crocin），交叉后与原文吻合度**下降**而非上升——
因为原文的"大靶点集合"本就是单库预测产物。因此交叉验证的用途是**提高可信度**，不是提高吻合度。

## 4. 强制预处理：UniProt ID 映射

```bash
python np_fetch.py uniprot --ids TP53 AKT1 EGFR --fromdb Gene_Name --todb Gene_Synonym \
  --output targets_mapped.csv
```
23 篇中因 ID 转换不一致导致的靶点丢失占比：**0%**（本机 5 篇 × 6,309 行 + 11 篇 × 8,102 行 + 4 篇 × 2,134 行，映射率均 100%）。
原因：UniProt 全量词典（129k 名称/20k 蛋白）在 v1.2 已内置。但**外部回归数据集**（如 DrugBank ID）仍需 REST idmapping。

## 5. 复现基准集（data/benchmark_23.json）

用法：
```bash
python np_bench.py score    # 置信度打分与链路断裂点
python np_bench.py audit    # 失败模式聚类
python np_bench.py infer --input paper.txt   # 阈值自适应推断
python np_bench.py explain --below 0.6       # 学术差异解释
python np_bench.py full     # 导出完整基准报告
```
每次升级后必须跑 `evals/regression.py` + 基准集回归，确保"修好一个、不弄坏另一个"。

## 6. 标准中间产物（np_export.py）

`herbs.csv / compounds.csv / targets_mapped.csv / PPI_network.sif / enrichment_results.csv`
列定义见 `np_export.py` 顶部 STANDARD。每次复现产出并用 `check` 校验。

## 7. 适用域声明（np_domain.py）

MW 150–500、XLogP -1–5 为参考；糖苷/矿物/肽类/极端 MW 自动标"预测置信度低，需实验验证"。
crocin（MW 977）实测判 0.40 低置信度——解释了其 TCMSP 仅 1 靶点的根因。
