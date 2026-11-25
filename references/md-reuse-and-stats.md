# MD 轨迹复用与统计功效协议（v3.0.5）

面向"手头有一条别人跑完的真实轨迹，如何合法地拿它当证据"。

## 1. 先盘资产，再决定是否重跑

| 位置 | 典型内容 |
|---|---|
| `D:\WorkBuddy\pharma-tools\mm_work\md\<体系>\` | OpenMM 生产段 `prod*.dcd`、`prod_meta.txt`（已跑 ns 数）、`prod_log.txt`、`*_checkpoint.chk`、`complex_raw.pdb` |
| `<pipeline>\e7_output*\` | MM-GBSA 逐帧 `e7_mmgbsa_perframe.csv` + `summary.csv` + `json` |
| `C:\Users\huawei\.workbuddy\skills\GeneDrug-AutoMD-Skill` | MD 生产/复现框架 |
| WSL `/root/black8/work/<PDBID>/` | Amber 复合物拓扑（com/rec/lig prmtop + 单帧 nc） |

判定"能否复用"的三件事：**① `prod_meta.txt` 的目标 ns 是否达成；② 段文件是否有重叠；③ 拓扑与 DCD 原子顺序是否一致**（本体系 NPROT=4851、NLIG=58、复合物 4909，溶剂化 DCD 每帧 78949 原子）。

## 2. 必须排除续跑重叠段

OpenMM 续跑是"从最近检查点状态重新起跑"，不是"接着检查点之后的帧继续"。
因此 `prod_seg1.dcd` 与 `prod_seg2.dcd` 会在 2.1/2.2 ns 重叠，且因积分器随机流不随
`setState` 恢复，两段是**同一区间的两次独立实现**，不能叠加统计。

- 规则：从 `prod_seg2.dcd` 起用；`prod_seg1.dcd` 一律 drop。
- 自动化：`md_traj_audit.py --drop prod_seg1.dcd`（默认值）；`e7_mmgbsa_Mv2.py` 内置 `DROP` 集合。

## 3. 报告 mean 时必须同时报"这 n 帧值多少"

```bash
python scripts/np_md.py stats --csv <perframe.csv> --col dG_bind --block-ns 1.0
```

输出：`sem_naive`、`sem_block`、`sem_inflation`、`tau_frames`、`ess`、前后半程 `z`、`ci95_block`、`verdict`。

判读阈值：

| 指标 | 含义 | 可用标准 |
|---|---|---|
| `sem_inflation` | 块 SEM / 朴素 SEM | >1.5 说明帧间强自相关，**只能报块 SEM** |
| `tau_frames` | 积分自相关时间 | 块长取 2–5×τ；否则块 SEM 仍偏小 |
| `ess` | 有效独立样本数 | <10 慎称"收敛"；<5 只能称"趋势" |
| `half_diff_z` | 前后半程差 / 块 SEM | \|z\|<2 视为均值已收敛 |

实测（Regaloside M × PTGDR，20 ns / 203 帧）：
- block=0.5/1/2/5 ns → SEM = 0.52 / 0.65 / 0.86 / 1.05 kcal/mol；τ≈1.9 ns，ESS≈10.5
- 前后半程 −50.73 vs −50.25（z=0.74）→ 均值收敛
- **结论：ΔG_bind = −50.5 kcal/mol（95% CI −51.8 ~ −49.2，块 SEM）；但块间 SD≈2.9，
  残余慢模态未被采样，写论文时应把"统计误差"与"采样系统误差"分开报。**

## 4. 轨迹稳定性证据（审稿人必查）

```bash
python tools/md_traj_audit.py --md-dir <Mv2 目录> --outdir <输出> --dg-csv <perframe.csv>
```

产出 `md_traj_perframe.csv` / `md_traj_summary.json` / `md_traj_panels.png`（骨架 RMSD、配体 RMSD + Rg、界面接触与氢键、ΔG 时序 + 1 ns 块均值）。

判读：
- **界面三项（接触数、最小距离、氢键）稳定** → 结合模式成立，即使骨架 RMSD 大。
- **骨架 RMSD 大且可逆** → 多为膜蛋白在无膜体系中的松弛；必须写明体系局限，
  不能拿它当"结构稳定性良好"的证据。
- 纯 numpy 解 DCD，不依赖 mdtraj；若目标环境的 `numpy.linalg.svd` 崩溃（本机 `envs/md`
  即如此），换系统 Python 3.13 + numpy 2.5 运行。

## 5. 不能做的事

- 不要把"同构象复制帧"当独立样本（会得到 Std=0，方差被结构性清零）。
- 不要用单条轨迹的 ΔG 差异去宣称两个配体的活性差异（需同体系多配体 + 重复轨迹）。
- 不要跨体系比较绝对 ΔG 数值（受体/配体/力场/溶剂模型任一不同即不可比）。
