# 盲法二次标注：操作手册（给未参与本项目的标注者）

**目的**：稿件 §"In this cohort, the largest attribution category was target provenance" 的 0/63
结论完全依赖归因标签。本手册让**第二个标注者在 20 分钟内**独立重标 25 篇，从而把"作者自己的
一致性检查"升级为"评分者间信度（Cohen's κ）"。

**不要做的事**：不要打开 `ground_truth.csv`（含作者标签），不要问作者"我这篇标得对不对"。
标注期间任何与作者关于具体条目的讨论都会破坏盲法。

---

## 你会拿到什么

`blind_sheet.csv`，25 行，列含义：

| 列 | 含义 |
|---|---|
| `pmid` | PubMed ID，用来查原文（如需） |
| `iso` / `year` | 期刊缩写 / 年份 |
| `herbs_requested` | 论文里写的药材数 |
| `herbs_found` | 我们能在数据库里找到的药材数 |
| `herbs_missing` | 找不到的药材名单（空 = 无缺失） |
| `n_active_repro` | 我们能重建的活性成分数 |
| `n_targets_repro` | 我们复现出的靶点数 |
| `report_targets` | 论文自报的靶点数 |
| `report_intersection` | 论文自报的"成分-疾病"交集靶点数 |
| `ratio_targets` | `n_targets_repro / report_targets` |
| `coverage` | `herbs_found / herbs_requested` |

**注意**：本表刻意不含作者的任何结论列（attribution、verdict）。

---

## 判断任务

对每一行，从下列 **5 个类别中选一个**，填入 `rater_label`；并给 1–3 的置信度
（1 = 猜的，2 = 大致有把握，3 = 有把握）。

| 标签 | 定义 | 典型证据 |
|---|---|---|
| `prov` | **靶点来源不一致**：论文使用的靶点预测来源与我们无法对齐，导致靶点集合本身不同 | `ratio_targets` 明显 < 1（常见 0.1–0.4），且 `herbs_missing` 为空 |
| `both` | 上一条 + 药材缺失同时存在 | `herbs_missing` 非空 **且** `ratio_targets` 明显 < 1 |
| `herb` | **仅药材缺失**：靶点数高于自报（我们找到的成分更多），但来源本身不冲突 | `herbs_missing` 非空，`ratio_targets` ≥ 1 |
| `thr` | **阈值/参数差异**：来源与药材都不是主因，差异来自筛选阈值、物种映射或 ID 转换 | 上述都不成立，差异幅度小且方向混合 |
| `none` | 记录不到明确差异 | 无上述任何证据，`ratio_targets` 接近 1 |

**判断顺序**（请严格按此顺序，以保证与作者的定义可比）：

1. `herbs_missing` 非空吗？→ 若非空，进入 2；若为空，进入 3。
2. `ratio_targets` < 0.6 吗？→ 是则 `both`；否则 `herb`。
3. `ratio_targets` < 0.6 吗？→ 是则 `prov`；否则进入 4。
4. `ratio_targets` 在 0.6–1.4 之间且无明显模式 → `none`；否则 `thr`。

**这个顺序是"操作化定义"，不是让你猜作者怎么想。** 如果你认为某条目按此顺序得出的标签与
你对原文的判断相冲突，请在 `rater_confidence_1to3` 填 1，并在旁边（另开一列或备注文件）
写一句理由——**这类冲突条目本身是最有价值的数据**。

---

## 完成后

1. 保存 `blind_sheet.csv`（保持文件名）。
2. 运行：

```bash
python kappa.py blind_sheet.csv ground_truth.csv
```

输出：各类别混淆矩阵、percent agreement、Cohen's κ 及 95% 置信区间、按置信度分层的 κ。

3. 把输出贴回给作者。**不要在运行前打开 `ground_truth.csv`**（`kappa.py` 会读它，你不需要看）。

---

## 怎么解读结果（作者预先声明，避免事后解释）

| κ | 解读 | 稿件怎么改 |
|---|---|---|
| ≥ 0.6 | 实质性一致 | Results 加一句 "attribution labels showed substantial independent agreement (κ = …)"，0/63 结论更稳 |
| 0.4–0.6 | 中等一致 | 保留 0/63，但把 "largest attribution category" 降为 "the category we could assign most consistently"，并在 Limitations 说明 |
| < 0.4 | 一致度低 | **不是失败**：改写为 "provenance attribution is itself hard to standardise"，并据此强化"要求机器可读的逐边来源、而非事后归因"的结论 |

三种结果都已经在稿件里预留了落点，所以**这个实验无论结果如何都不会白做**。
