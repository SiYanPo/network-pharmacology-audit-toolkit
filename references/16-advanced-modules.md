# 16. 高级分析模块使用指南（v1.5.0，23 篇套路驱动）

## 23 篇文献的高级手段频次（复现审计实测）
| 手段 | 频次 | 对应模块 |
|---|---|---|
| 分子对接 | 23/23 | np_dockprep/np_dock |
| WB/qPCR 湿实验 | 18/23 | np_lab |
| TCGA/生存/预后 | 17/23 | np_survival |
| MD 模拟 | 13/23 | 见下"MD 边界" |
| GEO DEGs | 1/23 | np_omics |
| 网络毒理学 | 1/23 | np_toxpath |
| 性味归经/君臣佐使 | 0/23 | np_tcmviz（补白） |
| WGCNA/scRNA | 0-1/23 | Roadmap v2.0 |

## 多组学交叉（np_omics）
```bash
python scripts/np_omics.py scan  --gse GSE14520 --outdir .      # 看分组
python scripts/np_omics.py deg   --gse GSE14520 --outdir . --logfc 0.5 --fdr 0.05
python scripts/np_omics.py cross --gse GSE14520 --targets intersection_genes.txt --outdir .
```
- 数据源：NCBI GEO FTP 直链（series matrix）+ 官方 GPL annot.gz 探针注释；偶发 404 自动重试。
- RNA-seq 类 GSE 无 series matrix → 需手工下载 suppl count 文件（脚本会明确提示）。
- 输出交集 CSV 直接可喂给 np_enrich（--genes）。

## 生存分析（np_survival）
```bash
python scripts/np_survival.py studies --kw lihc
python scripts/np_survival.py survival --study LIHC --genes STAT3,MMP9,PTGS2 --outdir surv
```
- 数据源：cBioPortal API（pan-can-atlas 2018 队列）；TCGA 缩写自动映射研究 id。
- 方法边界：log-rank 为精确 Mantel-Cox；HR 为 O/E 近似（非 Cox 多因素）；表达分组=中位数。
- 实测：LIHC 366 患者，MMP9 高表达 HR=1.49 p=0.023（与文献一致）。

## 敏感性分析（np_sens）
- 输入含 **pchembl 数值列**（ChEMBL 管线）→ 可真正重筛；纯 TCMSP 预筛串仅做 presence 分析。
- 输出：sensitivity_genes.csv（presence_rate=全阈值出现率）+ Jaccard 相似度矩阵 + <0.5 告警。

## Cytoscape（np_cyto）
- `network.xgmml`：Cytoscape File→Import→Network；节点带 type/degree 属性与圆形布局坐标。
- `.cys` 是 Cytoscape 会话二进制格式，**无法脱离 Cytoscape 程序生成**——这是诚实边界，
  XGMML+SIF 是官方可编程替代；排好版后另存 .cys 即可入投稿包。
- `supernetwork.svg`：成分-靶点-通路三层超网络矢量图，零依赖直接进论文。

## 湿实验指导（np_lab）
- 生成：UniProt 功能摘要（实时）+ qPCR 设计参数 + PrimerBank/CiteAb/HPA 检索链接
  + 病种匹配细胞系/动物模型（20 病种内置库）。
- **不编造**抗体 catalog 号与引物序列（幻觉风险）；给权威检索入口。

## 网络毒理学（np_toxpath）
```bash
python scripts/np_toxpath.py run --profile admet_profile.csv --threshold 0.7 --outdir tox
```
- 输入= np_admet predict 宽表；9 个风险端点→机制→KEGG 映射内置表；输出含
  "Safety Discussion 模板句"。药效/毒性同源通路（如 SR-ARE=Nrf2、PPARγ）标注"药毒两重性"。

## QA 报告（np_qa）
```bash
python scripts/np_qa.py report --dir repro_out --paper paper_params.json --outdir qa
```
- 七类产物完整性、SHA256 指纹、与原文数量对比（>30% 偏差 WARN）、审稿薄弱环节清单。

## MD 模拟边界（评估结论，未集成）
- 23 篇中 13 篇提及 MD（catechin 等为真用：GROMACS 2020.6+AMBER99SB+100ns+RMSD）。
- 本环境无 GPU/GROMACS，自动装 GROMACS+拓扑生成链路（CGenFF/acpype）工程量大且易碎。
- **务实替代**：对接后用 pose 复验 + 短程约束最小化（OpenMM pip 可装，CPU 也能跑 100ps 级
  稳定性初筛）→ 已列入 V2.0 Roadmap P2；全量 100ns 建议服务器跑 GROMACS 官方 workflow。

## 缓存与增量更新
- `np_fetch.py cache`：缓存统计/清空；源注册表（source_registry.json）记录每库访问日期与实测状态。
- GEO/UniProt/cBioPortal/DeepChem S3 下载全部落本地缓存，重复运行秒级命中。
