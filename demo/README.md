# 5 分钟可复现 demo

目的：让审稿人/读者在**不联网、不装任何第三方包**的前提下，亲自验证本 skill 的审计层能跑通，
并看到"结论可被判读"这一差异化能力。

## 一、为什么要单独做 demo

完整流水线依赖 TCMSP / STRING / OpenTargets 等在线数据源，离线必然失败。
但本 skill 的增量不在"再多算一遍通路"，而在**审计层**：

- 产物契约校验（列名/类型/行数是否被下游期待的样子）
- 证据账本（L0–L4 加权，每条边带来源）
- SCM0061-2021 评价评分卡（附录 A 的 16 个要素，0/1/2 分）
- 自包含 HTML 报告（含评分卡章节）

这一层**纯标准库**（`csv/json/glob/os/re/html/math`），没有 numpy/pandas 依赖，
所以能在 `python:3.12-slim` 里直接跑——这正是"可复现"该有的样子。

## 二、两种跑法

### 本机（推荐先跑这个）

```bash
bash demo/run_demo.sh
```

预期输出（约 10 秒）：

```
== 1) 产物契约校验 ==
[contracts] PASS：26 个契约，存在 21，硬失败 0
== 2) 证据账本（L0–L4 分布）==
  条目 1185 | 等级分布 {'L1': 281, 'L2': 904}
== 3) SCM0061-2021 评价评分 ==
SCM0061 评分：26/34 = 76.5%  等级 B  已开展口径 76.5%（26/34）
  可靠性: 9/12
  规范性: 10/10
  合理性: 7/12
  未检出：R5, P2, P6
== 4) 生成 HTML 报告（含评分卡章节）==
[报告] 已生成 /tmp/np_demo_report.html（8 KB）
== DEMO OK ==
```

### 容器（验证第三方环境可复现）

```bash
docker build -f Dockerfile.slim -t np-pro-slim:latest .
docker run --rm np-pro-slim:latest            # 默认就是跑 demo
```

镜像实测 185 MB（基础镜像占大头，未装任何 pip 包），构建数十秒（只需拉 `python:3.12-slim`，不装任何 pip 包）。

## 三、案例是什么

`demo/case/` 是一次真实运行（槲皮素 × 结直肠癌，101 靶点 → 23 交集 → 281 PPI 边）的产物，
共 79 个文件，保留了所有中间表与报告。它**不是**为 demo 专门造的漂亮数据——
包括"未满分项 R5/P2/P6"（无阴性对照、无敏感性、无适用域）也原样保留。

这一点很重要：demo 展示的不是"我们满分"，而是"我们能给出分数并说明为什么不满分"。

## 四、怎么用自己的案例

```bash
python scripts/np_contracts.py check --dir /path/to/your/np_out
python scripts/np_scm0061.py      --dir /path/to/your/np_out
python scripts/np_report.py       --dir /path/to/your/np_out -o report.html
```

只要有 `pipeline_manifest.json`，评分卡就能自动判断哪些阶段跑过、哪些是"跑了但降级"。

## 五、评分口径（避免误读）

| 分 | 含义 |
|---|---|
| 2 | 产物存在且含判定所需字段 |
| 1 | 产物存在但缺关键字段 |
| 0 | 未检出该要素产物（可能"未开展"，也可能"跑了但软失败"——报告会区分） |

- `percent`：全要素口径（分母固定 32）；
- `percent_executed_scope`：只统计本次运行开展过的阶段，避免把"主动不做"算成"做错"。

**靶点数一致 ≠ 机制结论正确。** 评分卡衡量的是"证据留痕完整度"，不是"科学正确性"。
