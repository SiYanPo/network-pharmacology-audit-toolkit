# 19. 500 篇高级手段审计与 V2.0 路线图（v1.8.0）

## 一、500 篇全量扫描（精确计数，可复算）

| 手段 | 篇数/500 | 占比 | skill 现状 |
|---|---|---|---|
| 分子对接 | 500 | 100% | ✅ np_dockprep/np_dock |
| 动物实验（严格） | 298 | 60% | ✅ np_lab 模型推荐 |
| 湿实验（WB/qPCR） | ~270 | 54% | ✅ np_lab 方案骨架 |
| **MD 模拟** | **150** | **30%** | ✅ OpenMM + GROMACS(WSL) |
| 其中 GROMACS | 93 | 19% | ✅ 已实测 51 ns/day |
| 其中 AMBER | 45 | 9% | ⚠️ AmberTools 在 WSL（可跑 MMPBSA.py） |
| 其中 Desmond | 52 | 10% | ❌ 商业软件，不集成 |
| **MM-PBSA/GBSA** | **70** | **14%** | ✅ 新增 `mmpbsa`（gmx_MMPBSA 1.7.1） |
| GEO DEGs | 96 | 19% | ✅ np_omics |
| TCGA/生存 | 77 | 15% | ✅ np_survival（cBioPortal） |
| 机器学习 | 45 | 9% | ⚠️ 部分（ADMET 用 RF；对接打分未做 ML） |
| 代谢组 | 80 | 16% | ⚠️ 仅文本统计，无分析模块 |
| 肠道菌群 | 60 | 12% | ⚠️ 同上 |
| 蛋白组 | 45 | 9% | ⚠️ 同上 |
| WGCNA | 23/650 | 4% | ⚠️ R 生态，Python 不成熟 |
| scRNA-seq | 23/650 | 4% | ⚠️ CELLxGENE 可接，数据重 |
| 网络毒理学 | 5/650 | 1% | ✅ np_toxpath |
| 君臣佐使 | 4/650 | 1% | ✅ np_tcmviz（可视化，非预测） |
| 性味归经 | 2/650 | 0% | ✅ np_tcmviz |

### 图表风格（词频）
heatmap 46% ｜ bubble plot 33% ｜ Sankey 9% ｜ chord 3% ｜ Circos 1%
→ skill 已有：气泡图（np_enrich）、网络图（np_ppi）、热图（可加）、桑基（np_tcmviz）。
**缺口**：热图/弦图/Circos 的自动生成（优先级中）。

## 二、V2.0 路线图（三档，按优先级）

### A 档：现在就能做（本轮已做/可立即做）
| 功能 | 状态 | 工具/API |
|---|---|---|
| GROMACS 全流程（WSL） | ✅ 已实测 51 ns/day | gmx 2026.3 @ /opt/md/envs/gmx |
| MM/PBSA 结合自由能 | ✅ 已装 gmx_MMPBSA 1.7.1 | /opt/md/envs/mmpbsa |
| GEO DEGs 交叉 | ✅ | NCBI FTP + annot.gz |
| TCGA 生存 | ✅ | cBioPortal API |
| Cytoscape XGMML/SIF | ✅ | XGMML 1.1 |
| 湿实验方案 | ✅ | UniProt + PrimerBank/CiteAb 链接 |
| 网络毒理学 | ✅ | 内置映射 + np_enrich |
| QA 报告 | ✅ | np_qa |
| 敏感性分析 | ✅ | np_sens |
| 500 篇语料流水线 | ✅ | np_corpus / np_repro_batch |

### B 档：需您提供资源（服务器/账号/软件）
| 功能 | 需要什么 | 实施方案 | 优先级 |
|---|---|---|---|
| **全量 MD 生产（20–100 ns×多副本）** | 本机可行但慢（20ns≈10h@1万原子；膜蛋白需集群） | np_md gmx --ps 20000，夜里跑；集群用同 mdp | P1 |
| **Desmond（10% 论文）** | Schrödinger 商业 license | 不集成；提供 mae 输入转换说明 | P3 |
| **DrugBank 金标准验证** | 注册 go.drugbank.com/pricing | np_admet 加 external-gold 子命令 | P2 |
| **scRNA（4%）** | CELLxGENE ~GB 数据 + 16GB 内存 | scanpy 读 h5ad → 与 NP 靶点交集 | P3 |
| **代谢组/菌群（16%/12%）** | 原始数据（论文 SI） | 需 SI 自动下载+表格解析 | P2 |
| **CNKI 中文核心** | 机构登录（您 2 分钟） | bsk 浏览器接管 kns.cnki.net | P1 |

### C 档：算法/技术瓶颈（诚实标注）
1. **君臣佐使自动预测**：无标注数据集，现仅做可视化与人工分级 → 不做伪实现。
2. **Desmond 轨迹分析**：商业格式，无开源解析器。
3. **WGCNA**：Python 生态不成熟（PyWGCNA 维护弱），建议 R 脚本交接。
4. **对接打分 ML 化**：需训练集与 GPU，超出轻量基线定位。
5. **补充材料（SI）自动解析**：PDF 表格结构千变万化，OCR+布局解析成本高（48% 未披露数字多藏于此）。

## 三、本轮新增的可执行能力（已落地）
```bash
# GROMACS 全流程（WSL，自动路径转换）
python scripts/np_md.py gmx --pdb prep/system.pdb --ps 20000 --nt 6 --outdir md_gmx
# 轨迹分析（mdtraj）
python scripts/np_md.py analyze --dir md_gmx --backend gmx
# MM/PBSA（需复合物 tpr + 轨迹 + 索引分组）
python scripts/np_md.py mmpbsa --tpr md_gmx/md.tpr --traj md_gmx/md.xtc \
  --index md_gmx/index.ndx --receptor-group Protein --ligand-group LIG --outdir mmpbsa
```

## 四、实测记录（2026-11-20）
- GROMACS 2026.3（conda-forge, WSL Ubuntu-24.04）：`/opt/md/envs/gmx/bin/gmx`
- ubiquitin（18,869 原子，6 线程）：em→nvt(50ps, POSRES)→npt(50ps, POSRES)→md(100ps)
  **51.0 ns/day**，RMSD 末帧 0.41 nm（mean 0.37），Rg 2.94 nm
- gmx_MMPBSA 1.7.1：`/opt/md/envs/mmpbsa/bin/python`（AmberTools 在 /opt/ambertools）
- AmberTools：MMPBSA.py / antechamber / parmchk2 可用（GAFF2 配体参数化备用路线）
