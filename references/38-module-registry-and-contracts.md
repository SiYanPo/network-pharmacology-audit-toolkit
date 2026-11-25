# 33 · 模块注册表与产物契约（v3.0）

## 1. 问题

评审自审结果：80 个模块里只有 42 个进 DAG，23 个"执行图不可达"，没有任何测试断言产物的列名。
后果是两条：

1. **死代码与活代码混在一起**——审稿人问"你到底用了哪些模块"答不清；
2. **接口静默退化**——改一个列名，两周后出图时才炸。

## 2. 注册表（`config/module_registry.yaml`）

### 2.1 五种 status

| status | 含义 | 准入条件 |
|---|---|---|
| `entry` | 入口（`np_pipeline` / `np_run`） | 必须是用户直接调用的编排层 |
| `core` | DAG 主流程节点 | 必须在 `config/pipeline.yaml` 的 `stages` 里 |
| `audit` | 按需审计工具 | 必须在 `audit_tools` 里，且能独立跑 |
| `manual` | 手动工具 | **必须写 `when`（何时用）+ `min_cmd`（最小命令）**，否则 `test_contracts` 红 |
| `deprecated` | 待移除 | 不得被 `pipeline.yaml` 引用（测试会查） |

### 2.2 当前分布（80 模块）

| status | 数量 | 例子 |
|---|---|---|
| core | 13 | alias / standardize / intersect / ppi / enrich / ortholog / report |
| audit | 33 | benchmark / fdr / split / snapshot / contracts / finalcheck |
| manual | 30 | dock / md / primer / survival / coloc_compare / tools |
| deprecated | 2 | `np_prov`（被 np_snapshot 取代）、`np_doctor`（被 np_preflight 取代） |
| entry | 2 | np_pipeline / np_run |

### 2.3 维护规则

```bash
python scripts/np_selfaudit.py run --run-selftests   # 产出 audit/skill_audit.md，含 registry_status
pytest tests/test_contracts.py -q                    # 未登记模块 / manual 缺命令 / deprecated 被引用 → 红
```

新增模块时三件事：写 selftest、进 registry、若进 DAG 则加契约。

## 3. 契约（`np_contracts.py`）

### 3.1 写法

```python
"05_ppi/hubs.csv": _c(
    ["node", "degree"],                              # 逻辑列名（必须存在）
    {"node": ["node", "gene", "symbol"]},            # 物理别名（任一命中即可）
    {"degree": "number", "betweenness": "number"},   # 类型抽查（前 200 行）
    note="node 与 degree 是 G21 与拓扑章节的输入"),
```

四类产物：

| kind | 校验方式 | 用于 |
|---|---|---|
| `table` | 表头 + 别名 + 类型 + 行数 | CSV/TSV |
| `list` | 非空行数 | `intersection_genes.txt` |
| `json` | 可解析 + 必需键 | manifest / decision_rules |
| `report` | 非空 | md / html |

### 3.2 覆盖

`np_contracts coverage` 断言"`pipeline.yaml` 声明的每个 artifact 都有契约"：
当前 **18/18 覆盖，24 条契约**。新增 stage 若忘了写契约，CI 会红——这正是我们要的。

### 3.3 与门禁的分工

- **契约**：文件"形状"对不对（列在不在、类型对不对、空不空）。改代码就能立刻发现。
- **门禁（G1–G23）**：内容"够不够投稿"（有没有真跑、有没有校准、方向是否矛盾）。

两者互补：契约防"下游读不到"，门禁防"读到了但结论不成立"。

## 4. 覆盖率：棘轮而不是口号

评审建议 `--cov-fail-under=55`，但当前真实覆盖率 15.22%（大量逻辑在 `main()` 里，
没拆纯函数）。直接设 55 会让 CI 永远红，最后一定被人关掉——比不设更糟。

因此用**棘轮**（`tools/check_coverage.py` + `coverage_baseline.txt`）：

- CI 只断言"不低于基线 − 0.5 个百分点"；
- 提高覆盖率后 `--write` 更新基线，且**只允许上调**（下调会被拒绝）。

这样"覆盖率不许退步"是可执行的，"半年内到 55%"仍是待办——写在报告里，不装成已完成。

## 5. 新增模块的清单（照着做）

1. 文件头写"为什么需要 + 诚实边界"（不是功能说明书）。
2. `--selftest` 返回 `0/1`，并打印每项断言。
3. 输出写 manifest（输入哈希/参数/种子/结果哈希）。
4. 缺失依赖 → 状态字段（`TOOL_MISSING`/`UNAVAILABLE`），**不伪造结果**。
5. 进 registry；进 DAG 则加契约；进 registry 的 `manual` 必须给最小命令。
6. 路径不许写死机器目录（用 `np_tools`）。
