# Skill 自审报告（np_selfaudit）

生成时间：2026-11-22 18:44:51　模块数：80

## 1　覆盖率总览

| 指标 | 数量 | 占比 |
|---|---|---|
| 模块总数 | 80 | 100% |
| 进入 DAG（config/pipeline.yaml） | 47 | 59% |
| 被 evals 覆盖 | 15 | 19% |
| 有 selftest | 44 | 55% |
| `--help` 可用 | 0（80 未检测） | 用 --with-help 可全量探测（并发 4，约 +3–8 分钟） |
| 执行图不可达（未登记死代码） | 0 | — |

| 已登记 manual/deprecated | 32 | — |

## 3　执行图不可达模块

这些模块只在文档里出现，DAG、evals、其他脚本都不调用。**不一定是坏事**（高级分析本来就是手动触发），但必须明确它们是「手动工具」还是「遗留死代码」：

（无）

## 5　第三方依赖

| 包 | 状态 |
|---|---|
| h5py | installed |
| matplotlib | installed |
| networkx | installed |
| numpy | installed |
| pandas | installed |
| rdkit | installed |
| scipy | installed |
| sklearn | installed |
| yaml | installed |

## 6　结论与动作

