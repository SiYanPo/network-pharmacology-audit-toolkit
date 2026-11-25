# 证据账本、校准与全局多重性（v2.11）

## 为什么必须有

网络药理学流程常见失败不是“算不出来”，而是**把所有模块的指标并排展示，却没有唯一排序函数**。
审稿人会问：`evidence → score → decision` 的映射是什么？如果答“多维展示供判断”，那是仪表盘，不是方法。

## 1. 证据账本

`np_evidence_ledger.py build --dir <outdir>` 从真实产物汇总，每行一个 `(gene, evidence_type)`：

`gene, source, evidence_type, level(L0–L4), direction, effect_size, uncertainty, provenance_hash`

- L0 文本共现；L1 计算预测；L2 数据库注释；L3 实验支持；L4 遗传因果（MR/COLOC）。
- 缺模块**不补行**：账本必须反映真实存在的证据，不能用 0 填满。
- 每行带来源文件 SHA256 前 16 位，报告可逐项回溯。

## 2. 校准分数（不是拍脑袋加权）

`fit --ledger ... --gold gold.txt -o 09_evidence` 在基准金标准上拟合逻辑回归，输出：

- `P(gold)` + bootstrap 95% CI；
- 样本内/5 折 OOF 的 Brier、ECE、AUC、可靠性曲线、DCA；
- 系数表（含负权重；负权重不是 bug，必须如实展示）。

**可用性闸门**：OOF AUC < 0.60 或全部系数为 0 → `WEAK_CALIBRATION`；
`score` 返回 `UNCALIBRATED`，只给先验排序，**不得把 p_gold 写成概率**。
2026-11-22 实测：真实 24 靶点、10 正例，OOF AUC=0.471、系数全 0 → 正确拒绝。

决策规则写死在 `DECISION_RULES`：

- `P≥0.70` 且 ≥2 类独立证据 → `wetlab_candidate`；
- `P≥0.70` 但单一证据类 → `hypothesis`；
- 其余 → `insufficient`。

## 3. 全局 FDR

`np_fdr.py` 汇总 MR IVW p、富集超几何 p、基准超几何 p、零模型经验 p，
同时做**族内 BH**与**全局 BH**，输出 `q_family` / `q_global`。

- `q_global` 是唯一可写入结论的显著性口径；仅 `p<0.05` 只能写“提示性关联”。
- `np_mr.py --fdr` 对同批靶点做 BH；`SINGLE_HYPOTHESIS` 时 q=p，但不得据此宣称因果显著。
- 报告自动生成中英双语“多重性声明”，含检测数、族数、校正前后显著数。

## 4. 诚实边界

- 校准器只在拟合基准的定义域内有效；换疾病/数据库必须重拟合。
- 证据级别是强度先验，不是真值；模型可能学到负权重。
- 金标准若与输入同源，AUC/校准仍会乐观；必须配 LODO/complex-aware 切分（`np_split.py`）。
