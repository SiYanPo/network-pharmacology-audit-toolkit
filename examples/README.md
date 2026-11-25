# 演示数据（构造数据）

> ⚠️ **本目录全部为构造数据（synthetic），仅用于跑通流程与验证环境，不得用于发表、不得作为任何科研结论的依据。**

生成方式：`python make_demo.py`（固定随机种子 20261119，结果可复现）。

| 文件 | 内容 | 刻意包含的"坑" |
|---|---|---|
| `compounds.csv` | 28 个成分（3 味药） | 含低 OB/DL 成分、3 个泛成分（槲皮素/山柰酚/β-谷甾醇） |
| `compound_targets.csv` | 293 条成分-靶点关系 | 靶点用 TCMSP 风格蛋白全名；混入 2 条非人类/噪音 |
| `disease_targets.tsv` | 135 个疾病基因 | score 长尾分布，来源混杂（OMIM/CTD/GeneCards/DisGeNET） |
| `string_edges.tsv` | 435 条 PPI 边 | 含 1 条自环、1 条重复边 |
| `pathways.gmt` | 14 条通路 | 含 "Pathways in cancer"（伞状）与 "Hepatitis B"（感染类） |
| `docking_scores.csv` | 12 对对接结果 | 含 1 个重对接 RMSD=2.6 Å 的失败案例 |

## 为什么故意放"坑"

演示数据的价值不在于"跑得通"，而在于**让每个质控模块都有事可做**：
缺口径会触发非人类剔除、低 OB 会进剔除清单、泛成分会触发警示、
伞状通路会被降级、重对接失败会被门禁拦下。跑一遍就能看懂每个模块在防什么。

## 跑法

```bash
python scripts/np_run.py --demo
```
