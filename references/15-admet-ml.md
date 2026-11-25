# 15. ML 型 ADMET 预测与外部验证（v1.4.0+）

## 何时使用
- 需要对网络药理学选出的活性成分（尤其是新化合物/苷元/代谢物）给出 ADMET 风险谱时；
- 审稿人质疑"只有 OB/DL 过滤，没有现代 ADMET 评估"时；
- 需要**外部验证**证明模型泛化能力（而非内部交叉验证）时。

## 数据通路（全部可自动化，无需注册）
| 数据 | 来源 | 状态 |
|---|---|---|
| TDC ADMET/Tox 官方基准（23+ 端点，官方 scaffold split） | HuggingFace `maomlab/TDC`，经 `hf-mirror.com` | ✅ 自动下载 |
| BBBP（跨库外部验证用） | DeepChem S3 `deepchemdata.s3.amazonaws.com` | ✅ 自动下载 |
| Harvard Dataverse（TDC 原始托管） | dataverse.harvard.edu | ❌ 沙箱代理 502 拦截 |
| Zenodo（第三方冻结快照 20180944） | zenodo.org | ❌ 沙箱代理 502 拦截 |

> 说明：hf-mirror 路径拿到的是官方 TDC 数据的净化镜像（含官方 scaffold split 列），
> 与 `PyTDC` 包下载内容同源； zenodo 上的 `tdc_admet_group_2026-03-24.tar.gz` 是版本锁定
> 快照，若代理放行可作为版本固定选项（SHA256 `8ac217bd8c316d…`）。

## 命令速查
```bash
# 1) 下载基准数据（ADME 4MB + Tox 79MB）
python scripts/np_admet.py download --outdir admet_raw

# 2) 按端点导出 CSV
python scripts/np_admet.py extract --parquet admet_raw/ADME.parquet admet_raw/Tox.parquet \
  --outdir endpoints --tasks "bbb_martins,ames,dili,cyp3a4_veith"

# 3) 训练（默认 23 端点 × 5 seeds；RF-200 + Morgan2048 + 12 描述符）
python scripts/np_admet.py train --datadir endpoints --outdir results --seeds 5

# 4) 跨库外部验证（bbb_martins → BBBP 独立测试集）
python scripts/np_admet.py external --datadir endpoints --modeldir results --outdir results

# 5) 预测化合物的 ADMET 谱
python scripts/np_admet.py predict --modeldir results \
  --input compounds.csv --output admet_profile.csv

# 6) 生成 Markdown 报告
python scripts/np_admet.py validate --results results/admet_metrics.csv
```

## 官方指标约定（勿改）
- 均衡二分类 → **AUROC**；失衡二分类（CYP 抑制四端点）→ **AUPRC**；
- 回归 → **MAE**（Caco2/Lipo/PPBR/LD50）；半衰期/清除率 → **Spearman**（受结构外因素影响）。

## 外部验证的分层语义（写论文时照此表述）
1. **官方 scaffold-split test**：训练期间从未接触的骨架外样本 → 报告 mean±std（5 seeds）；
2. **跨库验证**：不同文献来源的独立数据集（BBBP vs bbb_martins）→ 报告 cross-dataset AUROC；
3. **排行榜对照**：与 TDC 官方 leaderboard（https://tdcommons.ai/benchmark/admet_group/overview/）
   顶级模型对比，注意基线模型（RF+FP）≠ 大模型（ChemBERTa 等），定位为"轻量可复现基线"。

## 与 NP 管线的衔接
- `np_domain.py` 的适用域输出（MW/logP/骨架类型）→ 解释 ADMET 预测置信度；
- Tox21 四端点与机制叙事对应：**SR-ARE**（氧化应激/铁死亡）、**NR-PPAR-gamma**（芒柄花素）、
  **NR-AhR**（芳烃受体）、**NR-ER**（雌激素样作用）——预测结果可直接进入 Discussion；
- 糖苷类（如 crocin，MW>900）属训练集化学空间边缘 → `predict` 会给 uncertain 标记，
  此时建议先走苷元两步管线（v1.3 规划）再评估苷元的 ADMET。

## 已知边界
- RF+Morgan 是**轻量基线**：不追求 SOTA，追求"可复现+零 GPU+外部验证完备"；
- `herg_central_*`（30 万行）默认排除，需要时用 `--tasks` 显式加入；
- 沙箱代理对 dataverse/zenodo 的拦截是环境问题，非数据源问题；离线环境可预下载 parquet。
