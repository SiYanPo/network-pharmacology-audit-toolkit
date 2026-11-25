# 投稿前 Final Check

**判定：DRAFT（不可投稿）**

检查时间：2026-11-23 14:17:58　目录：D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full

| 项 | 名称 | 级别 | 状态 | 说明 |
|---|---|---|---|---|
| G1 | 时间切分盲测 | hard | FAIL | 缺文件 temporal_validation.json |
| G2 | Apo/Holo 状态标注 | hard | FAIL | 缺文件 dock_state_annotated.tsv |
| G3 | 剂量权重 | hard | FAIL | 缺文件 dose_weighted_hubs.csv |
| G4 | 配伍禁忌质检 | soft | PASS | D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full\00_safety\cont |
| G5 | 外部金标准盲测 | soft | FAIL | 缺文件 goldtest_result.json |
| G6 | 图件自检 FAIL=0 | soft | FAIL | 缺文件 figure_qa.md |
| G7 | 证据锚点校验 | soft | FAIL | 缺文件 evidence_report.json |
| G8 | Decoy 阴性对照(无 Red Flag) | hard | FAIL | 缺文件 decoy_summary.json |
| G9 | 单细胞投影 | soft | FAIL | 缺文件 singlecell_manifest.json |
| G10 | 亚细胞定位分层 | soft | FAIL | 缺文件 localization_layers.csv |
| G11 | 引物特异性验证 | soft | FAIL | 缺文件 specificity_plan.json |
| G12 | MD 桥接 | soft | FAIL | 缺文件 md_bridge_meta.json |
| G13 | 共识打分UQ(无Unstable) | hard | FAIL | 缺文件 consensus_summary.json |
| G14 | 生理屏障过滤 | soft | PASS | D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full\09_funnel\funn |
| G15 | 组织特异性校验 | soft | FAIL | 缺文件 hpa_summary.json |
| G16 | 细胞通讯(LIANA) | soft | FAIL | 缺文件 liana_status.json |
| G17 | 配体拓扑就绪 | soft | FAIL | 缺文件 ligand_top_result.json |
| G18 | MR 因果证据 | soft | FAIL | 缺文件 mr_result.json |
| G19 | 协同拮抗分析 | soft | FAIL | 缺文件 synergy_manifest.json |
| G20 | LINCS 重定位 | soft | PASS | D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full\repurposing_ca |
| G21 | 基因数据流一致性(交集⊆PPI) | hard | PASS | 22/23 个交集基因在 PPI 节点中（主网络 coverage 0.9565）；缺失：MGAM2 |
| G23 | 证据聚合与校准(ECE<0.10) | soft | FAIL | 有账本但无校准器：概率不可报告（UNCALIBRATED），只能给先验排序 |
| G22 | 同源类型隔离(FuncRep≠直系同源) | hard | PASS | 3 个功能替代均带 (FuncRep) 标记，且警示文件齐备 |

## 硬性项（必须全部通过才能出 Ready 水印）

- G1 时间切分盲测：✗ 缺失 → 补：`python np_defense.py temporal --gold ... --pred ... --cutoff-year 2020`
- G2 Apo/Holo 状态标注：✗ 缺失 → 补：`python np_defense.py dockstate --plan ...`
- G3 剂量权重：✗ 缺失 → 补：`python np_doseweight.py --edges ... --dose ...`
- G8 Decoy 阴性对照(无 Red Flag)：✗ 缺失 → 补：`python np_decoy.py score --core ... --decoy ...；n_red_flag 必须为 0`
- G13 共识打分UQ(无Unstable)：✗ 缺失 → 补：`python np_consensus.py ...；n_unstable 必须为 0，否则标 Unstable_Binding`
- G21 基因数据流一致性(交集⊆PPI)：✓ 通过
- G22 同源类型隔离(FuncRep≠直系同源)：✓ 通过

## 自愈指令（逐项可复制执行）

> 以下每条命令可直接复制运行；运行后重新执行本脚本直到判定 READY。

### G1 时间切分盲测 — 缺文件 temporal_validation.json

**诊断**：时间切分盲测缺失：需 CTD 金标准 + 预测靶点做 temporal split

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_defense.py" temporal --gold gold/gold_pd_genes.tsv --pred disease_targets.tsv --cutoff-year 2020 -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full
```

### G2 Apo/Holo 状态标注 — 缺文件 dock_state_annotated.tsv

**诊断**：缺少 Apo/Holo 构象状态标注：对接结构状态未声明

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_defense.py" dockstate --plan docking_plan.tsv -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full
```

### G3 剂量权重 — 缺文件 dose_weighted_hubs.csv

**诊断**：缺少剂量权重分析：无法证明成分暴露可达

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_doseweight.py" --edges compound_targets.tsv --dose dose.csv -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full
```

### G5 外部金标准盲测 — 缺文件 goldtest_result.json

**诊断**：外部金标准盲测缺失：结果不可证伪

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_goldtest.py" --gold-ctd gold/gold_pd_genes.tsv --pred disease_targets.tsv --decoy gold/decoy.tsv -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full
```

### G6 图件自检 FAIL=0 — 缺文件 figure_qa.md

**诊断**：图件自检存在 FAIL：需修复图片字号/重叠/越界

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_figcheck.py" --dir figs --out figs/figure_qa.md
```

### G7 证据锚点校验 — 缺文件 evidence_report.json

**诊断**：Discussion 存在无证据锚点的断言句

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_evidence.py" validate --text discussion.txt --datadir D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full --out D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full/evidence_report.json  # 按输出逐句补 [Data:]/[PMID:]
```

### G8 Decoy 阴性对照(无 Red Flag) — 缺文件 decoy_summary.json

**诊断**：Decoy 对照触发 Red Flag：对接假阳性风险高

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_decoy.py" score --core docking_scores.csv --decoy decoy_scores.csv -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full  # 检查 Z<-2 才可作为特异结合
```

### G9 单细胞投影 — 缺文件 singlecell_manifest.json

**诊断**：单细胞投影为 DEMO 或缺失：不可投稿

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_singlecell.py" project --h5ad data.h5ad --targets core_targets.txt --cluster-key celltype -o figs/singlecell
```

### G10 亚细胞定位分层 — 缺文件 localization_layers.csv

**诊断**：亚细胞定位缺失：无法排除核内错配

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_wetlab.py" localization --targets core_with_uniprot.csv -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full
```

### G11 引物特异性验证 — 缺文件 specificity_plan.json

**诊断**：引物特异性未验证：可能非特异扩增

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_wetlab.py" specificity --primers D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full/primer/mouse_primers.csv --genome GRCm39.fa -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full
```

### G12 MD 桥接 — 缺文件 md_bridge_meta.json

**诊断**：MD 桥接缺失：静态对接稳定性未验证

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_wetlab.py" md-bridge --receptor dock/receptor.pdbqt -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full/md
```

### G13 共识打分UQ(无Unstable) — 缺文件 consensus_summary.json

**诊断**：共识打分不稳定（Unstable_Binding）：单一打分函数偏差大

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_consensus.py" offline --scores three_scores.csv -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full  # 或补 GNINA 重打分：np_consensus.py scan --gnina ...
```

### G15 组织特异性校验 — 缺文件 hpa_summary.json

**诊断**：组织特异性校验缺失：无法排除细胞系假阳性

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_hpa.py" fetch -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full/hpa && "C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_hpa.py" check --targets core_targets.txt --tissue <疾病组织> --bulk D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full/hpa/rna_tissue_consensus.tsv -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full/hpa
```

### G16 细胞通讯(LIANA) — 缺文件 liana_status.json

**诊断**：细胞通讯分析缺失：LIANA 脚本未跑

```bash
bash D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full/liana_status.json 中给出的 run_liana.py  # 先 pip install scanpy liana
```

### G17 配体拓扑就绪 — 缺文件 ligand_top_result.json

**诊断**：配体 GAFF 拓扑缺失：MD 无法启动

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_ligand_top.py" gen --ligand dock/best_ligand.sdf -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full/ligand_top
```

### G18 MR 因果证据 — 缺文件 mr_result.json

**诊断**：缺少 MR 因果证据：靶点-疾病仅为相关性

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_mr.py" analyze --exposure exp.tsv --outcome out.tsv -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full/mr
```

### G19 协同拮抗分析 — 缺文件 synergy_manifest.json

**诊断**：缺少协同拮抗分析：多成分组合效应未回答

```bash
"C:\Users\huawei\.workbuddy-2\binaries\python\versions\3.13.12\python.exe" "C:/Users/huawei/.workbuddy-2/skills/network-pharmacology-pro/scripts/../scripts/np_synergy.py" --edges compound_targets.tsv --pathways enrich.csv -o D:\WorkBuddy\分身文件\2026-11-19-10-40-46\全流程实测\np_out\Quercetin_full/synergy
```

### G23 证据聚合与校准(ECE<0.10) — 有账本但无校准器：概率不可报告（UNCALIBRATED），只能给先验排序

**诊断**：未知项，请查阅对应模块 --help

```bash
# 无自动修复命令
```


## 水印声明

> 本水印仅表示 network-pharmacology-pro 的**自查清单**通过，
> 不代表期刊接受或同行评议通过。所有结论仍需实验验证与作者负责任解读。
