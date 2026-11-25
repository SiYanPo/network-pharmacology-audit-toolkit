# 论文图表与 Definition of Done（评审 §12–§13）

## 最小可投稿图表集

| 图 | 内容 | 数据来源 | 现状 |
|---|---|---|---|
| Fig1 | 架构 + 证据账本 `evidence→score→decision` | `np_evidence_ledger` | 可生成 |
| Fig2 | 校准曲线 + Brier/ECE | `calibration_model.json` | 需真实金标准 |
| Fig3 | random vs complex vs LODO | `np_split` | 可生成 |
| Fig4 | 与外部工具 ROC/PR | `np_split compare` | 需外部预测表 |
| Fig5 | 1 方剂案例 + 文献/湿实验验证 | `np_uireport` | 需真实案例 |
| Table1 | 门禁与覆盖率 | `np_finalcheck` / `np_selfaudit` | 可生成 |
| Supp | benchmark 全表、快照哈希、失败案例 | `np_snapshot` / `audit` | 可生成 |

## DoD（7 条，逐条可检查）

- [ ] `score_targets()` 输出 P(gold)+95%CI，且有校准曲线（ECE<0.10）
      当前：`np_evidence_ledger.score_targets` 已有；真实 24 靶点 OOF AUC=0.471 → 闸门拒绝，需扩金标准。
- [ ] 全局 FDR 落地，报告有“多重性声明”段
      当前：`np_fdr` + `np_uireport.multiplicity_block` 已接；pytest 覆盖。
- [ ] 基准含 3 类切分 × 2 个外部工具
      当前：`np_split` random/complex/LODO 已有；外部工具需用户提供预测表（`compare`）。
- [ ] pytest ≥90 用例，含契约测试；CI 双平台
      当前：108 passed / 1 skipped；`test_contracts`；Windows job 已加。
- [ ] `NP_SNAPSHOT` 离线复现 demo，产物哈希一致
      当前：`np_snapshot selftest` 7/7；2026-11-22 已真实 pack/mount，2 个公共文件哈希一致
      （`audit/snapshot_demo_20261122.json`）。
- [ ] 不可达模块全部标记 manual/deprecated
      当前：`config/module_registry.yaml` 78/78；selfaudit 未登记死代码 0。
- [x] COLOC 与 R `coloc` 20 位点一致性报告
      当前：R 4.6.2 / coloc 5.2.3 实测，20 位点 r=1.000、MAE=0.000，见
      `audit/coloc_vs_R_coloc_2026-11-22.md`；真实数据抽验仍见 roadmap V3。
