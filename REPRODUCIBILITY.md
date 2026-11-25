# 可复现性声明与核验清单

本文件是投稿用的可复现性证据包：**每一条主张都对应一个可执行的核验命令**。
读者不需要信任我们的声明，只需按下面的命令运行。

---

## 一、最小核验（2 分钟，无需联网、无需 pip 包）

```bash
# 1) 构建轻量镜像（仅标准库，实测 185 MB）
docker build -f Dockerfile.slim --build-arg BASE_IMAGE=docker.m.daocloud.io/library/python:3.13-slim -t np-pro-slim .

# 2) 跑 demo：契约校验 → 证据账本 → SCM0061 评分 → HTML 报告
docker run --rm np-pro-slim:latest
```

**预期输出（逐字可比对）**

```
[contracts] PASS：26 个契约，存在 21，硬失败 0
条目 1185 | 等级分布 {'L1': 281, 'L2': 904}
SCM0061 评分：26/34 = 76.5%  等级 B  已开展口径 76.5%（26/34）
  可靠性: 9/12
  规范性: 10/10
  合理性: 7/12
  未检出：R5, P2, P6
== DEMO OK ==
```

若输出与上述不一致，说明审计层被修改过——这正是**回归保障**要防的事。

---

## 二、逐项主张 → 核验方式

| # | 论文中的主张 | 核验命令 | 预期 |
|---|---|---|---|
| 1 | 产物契约可被机器校验 | `python scripts/np_contracts.py check --dir demo/case` | `硬失败 0` |
| 2 | 证据账本按 L0–L4 分级且带溯源哈希 | 见 demo 第 2 步 | 1185 条，L1/L2 分布如上 |
| 3 | SCM0061-2021 附录 A 可被自动打分 | `python scripts/np_scm0061.py --dir demo/case` | 26/34，B |
| 4 | 报告自包含且含评分卡章节 | `python scripts/np_report.py --dir demo/case -o /tmp/r.html` | 生成成功，含 `SCM0061` 锚点 |
| 5 | 回归测试覆盖契约与关键路径 | `python -m pytest tests -q` | `121 passed, 1 skipped` |
| 6 | 湿实验结果可按 n/p 定级 | `python scripts/np_wetlab_result.py --csv demo/wetlab_demo.csv` | L3/L2 逐条给理由 |

---

## 三、关键数字的出处（可追溯）

| 数字 | 出处文件 |
|---|---|
| 500 篇复现记录 | `_corpus500/06_repro500.csv`（19 列，含 coverage / verdict_targets / attribution） |
| 门禁 RR 2.85 [2.01,4.05] | `全流程实测/retro/retro_gate.json` |
| 阈值扫描 25/25 稳健 | `全流程实测/retro/retro_sweep.json` |
| 归因拆解 63/63 全偏低 | 同目录 `retro_gate.md` |
| SCM0061 评分 26/34 | `demo/case/scm0061_scorecard.json` |
| MD 重复性 ICC / 块 SEM | `全流程实测/md_reuse/replicate_stats*/reps_summary.json` |
| DrugBank 召回 16/66 | `全流程实测/drugbank_audit/drugbank_recall.json` |

---

## 四、已知限制（主动声明，避免审稿人自行发现）

1. **demo 只覆盖审计层**：完整流水线依赖 TCMSP/STRING/OpenTargets 在线数据源，离线不可跑。
   我们刻意把"可复现"落在审计层，因为那是差异化的部分。
2. **Docker 基础镜像**：默认 `python:3.12-slim`；国内直连 Docker Hub 可能失败，
   需用 `--build-arg BASE_IMAGE=docker.m.daocloud.io/library/python:3.13-slim`。
3. **回归 121 项**覆盖契约、路由与关键计算，不等于全覆盖。
4. **回顾性分析是二次分析**，非前瞻注册试验；382/500 篇因原文未报告靶点数被排除。
5. **靶点数一致 ≠ 机制结论正确**：本工作不主张生物学正确性，只主张过程可复核。

---

## 五、一键体检（可选）

```bash
python scripts/np_doctor.py      # 环境与数据源能力矩阵
python scripts/np_selfaudit.py   # 模块清单 / help 矩阵 / 契约 / 死代码
```

用于确认"工具装了但流程跑不通"这类假绿——构建期冒烟（`Dockerfile.slim` 的 `RUN`）
已把 `硬失败 0` 作为硬性断言。

---

## Reproducibility-and-reusability checklist (journal pilot)

The seven items below are the ones a reproducibility review asks about. Each is answered with a
file in the archive, not with a statement.

| # | Item | Status | Evidence |
|---|---|---|---|
| 1 | LICENSE file in the repository root | yes | `LICENSE` (MIT) |
| 2 | README describing what the software does and how to run it | yes | `README.md`, `ARCHIVE_README.md` |
| 3 | Install documentation | yes | `INSTALL.md` (audit layer with zero dependencies; container; full pipeline) |
| 4 | Public issue tracker | repository setting | `.github/ISSUE_TEMPLATE/{bug_report,reproduction_failure}.yml` — the second template exists specifically for "a documented claim did not reproduce" |
| 5 | Archived release with a DOI | archive prepared, DOI pending | `network-pharmacology-audit-toolkit-v3.0.19.{zip,tar.gz}` + `ARCHIVE_MANIFEST.tsv` (599 files, per-file SHA-256) |
| 6 | Installable / runnable without the authors' environment | yes | `docker build -f Dockerfile.slim` → `DEMO OK`; measured on a clean extraction, not in the working tree |
| 7 | Passes a linter | yes | `ruff check scripts tools tests` with the pinned rule set in `ruff.toml` → **0 findings** |

### What item 7 changed

Running the linter was not a formality. The pinned error-level rules found five real defects that 121
passing tests had not caught, because none of them are on a tested path:

| File | Defect | Consequence |
|---|---|---|
| `scripts/np_md.py` | `prev` used before assignment in the GROMACS prep chain | `NameError` on the first `grompp`, i.e. the whole prep path was dead |
| `scripts/np_md.py` | same bug in the WSL backend, as `inp` | second dead path |
| `scripts/np_mr.py` | `time.sleep` in the retry branch without `import time` | any HTTP 429/5xx retry would raise instead of waiting |
| `scripts/np_omics.py` | `cmd_cross` defined twice, first definition truncated and shadowed | dead code that reads like a real function |
| `scripts/np_funnel.py` | module-level `math` shadowed by a local import | benign, removed |

All five are fixed, and both dead paths now resolve their input structure explicitly
(`em` from `solv_ions`, then `nvt` from `em`, then `npt` from `nvt`). This is the concrete argument for
shipping a linter configuration with a methods paper: a test suite only covers paths someone thought
to test, and the untested paths are exactly where the audit's own credibility lives.
