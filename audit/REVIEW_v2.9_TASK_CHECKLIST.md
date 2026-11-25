# 评审报告 v2.9 → 现行 skill 执行清单

状态：TODO / DOING / DONE / VERIFIED / BLOCKED
图例：[已有]=报告 §10 已做，仅回归；[新增]=拔高新增；[验证]=需真实运行。

| # | 状态 | 类型 | 执行项 | 证据/文件 |
|---|---|---|---|---|
| 1 | DONE | 已有 | 回归 2.10 已修（12 selftest、3 核心自检、COLOC 0.6、缓存/离线、PLINK/acpype、preflight） | pytest / selftest |
| 2 | DONE | 新增 | 证据账本 schema | scripts/np_evidence_ledger.py |
| 3 | DONE | 新增 | 校准分数 P(gold)+95%CI+决策规则 | 同上 |
| 4 | DONE | 新增 | Brier/ECE/校准曲线/DCA 进入报告 | 同上 / np_report |
| 5 | DONE | 新增 | 接入 audit_tools/pipeline/SKILL/测试 | config/pipeline.yaml |
| 6 | DONE | 新增 | 全局 FDR | scripts/np_fdr.py |
| 7 | DONE | 新增 | MR --fdr | scripts/np_mr.py |
| 8 | DONE | 新增 | 多重性声明段 | np_report/np_uireport |
| 9 | VERIFIED | 验证 | top-20 命中数 + 校准为主指标 | np_benchmark |
| 10 | DONE | 新增 | complex-aware 分组切分 | np_benchmark |
| 11 | DONE | 新增 | LODO + ETCM2/TCMSP 独立性 | np_benchmark / np_goldtest |
| 12 | DONE | 新增 | 外部工具对照入口 | np_benchmark |
| 13 | DONE | 新增 | 基线预注册 | config/benchmark |
| 14 | DONE | 新增 | STRING 扩展隔离 / 证据分层 | np_ppi |
| 15 | DONE | 新增 | 扩展前后 hub rank Spearman | np_ppi |
| 16 | DONE | 新增 | API 版本元数据 | np_prov / np_snapshot.meta（release/etag/sha256/fetched_at） |
| 17 | VERIFIED | 新增 | NP_SNAPSHOT 离线快照 | audit/snapshot_demo_20261122.json（pack→mount 哈希一致） |
| 18 | DONE | 新增 | bootstrap 确定性测试 | tests |
| 19 | DONE | 新增 | 模块注册表 | config/module_registry.yaml |
| 20 | DONE | 新增 | 契约测试 | tests/test_contracts.py |
| 21 | VERIFIED | 新增 | CI 覆盖率/双平台 | .github/workflows/test.yml（unit/windows/coverage/selftest/integration） |
| 22 | DONE | 新增 | Python/requirements 策略 | requirements-lock.txt |
| 23 | DONE | 新增 | 工具注册表/去硬编码 | npcore/np_preflight |
| 24 | VERIFIED | 新增 | Docker 实测 | 2026-11-23 闭合：镜像 np-pro:3.0（sha256:22a9c57d, 9.99GB）已构建；容器内 pytest 119 passed/3 skipped、离线 demo 契约 PASS、完整冒烟 BUILD SMOKE OK 全部实测。见 DoD收尾/dod_status.json D7 |
| 25 | DONE | 新增 | HPA/zebrafish/MR 限制 | docs/report |
| 26 | PLANNED | 新增 | 外部验证/通路/湿实验/用户研究 | references/34-validation-roadmap.md（协议已写，未实际执行） |
| 27 | DONE | 新增 | 论文图表 + DoD | references/35；COLOC DoD 已 DONE |

## 全局反查补项（2026-11-22）

- [x] §14 四个外部求证点：ETCM2/TCMSP、STRING、HPA、Delphi → `references/36-external-evidence-audit.md`
- [x] COLOC 公式错误修复 + R 20 位点 PASS → `audit/coloc_vs_R_coloc_2026-11-22.md`
- [x] 外部证据防回归测试 → `tests/test_external_evidence_docs.py`
- [x] 注册新增模块 np_coloc_compare / np_tools → 112 passed / 1 skipped
