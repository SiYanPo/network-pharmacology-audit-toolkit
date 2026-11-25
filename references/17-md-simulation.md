# 17. MD 模拟模块（v1.6.0，本机环境长期记忆已固化）

## 本机 MD 环境长期记忆（data/md_environment.json，实测 2026-11-20）

| 组件 | 位置/事实 | 状态 |
|---|---|---|
| **OpenMM 8.6.1** | `D:/MDtools/envs/md/python.exe`（含 pdbfixer 1.12/rdkit/numpy/CUDA 运行时） | ✅ np_md.py 主力后端 |
| OpenMM 平台 | Reference/CPU/OpenCL（无独显 → CPU/OpenCL） | ✅ 可用 |
| **GROMACS** | WSL Ubuntu 24.04（`D:/MDtools/wsl-ubuntu/ext4.vhdx`）；**前会话已跑通 20ns 全流程**（em→nvt→npt→md，产物 `D:/MDwork/00_test/run/`） | ⚠️ wsl.exe 在沙箱被安全策略拦截；用户放行或终端手动跑 |
| conda-forge | gromacs 无 win-64 包 → Windows 原生安装不可行 | 已确认 |
| AutoDock Vina | `D:/MDtools/vina.exe` | ✅ |
| MD 工作报告 | `D:/MDwork/MD环境搭建与计算流程_王百合苷-DP1.md`（DP1-GPCR 膜蛋白体系方案） | 参考 |
| 力场建议 | 蛋白 AMBER ff14SB/CHARMM36m；配体 GAFF2+AM1-BCC（文献可追溯）；水 TIP3P；0.15M NaCl | 落实于 np_md.py |

## 基准（前会话实测，写论文 Materials 可引用）
- **11.27 ns/day @ 7 万原子**（i5-12400 6C12T / 8GB / 无独显）
- 内存：7 万原子纯水峰值 133 MB（8GB 无压力）
- 结论：**本机 ≤1ns 快筛 + 参数调试；20ns≈2 天（WSL gmx）；100ns+ 生产→学校集群/郑州超算**

## np_md.py 用法
```bash
# 环境体检（读 md_environment.json 长期记忆）
python scripts/np_md.py check

# 结构准备（pdbfixer：缺失残基/原子修复、去杂原子、加氢 pH7）
python scripts/np_md.py prep --pdb receptor.pdb --outdir md_prep

# 蛋白 MD（OpenMM 后端）
#   快筛：--ps 100（0.1ns）显式溶剂 PME+0.15M NaCl
#   正式：--ps 20000（20ns，同前会话 gmx 规模）；--implicit 为超快定性筛
python scripts/np_md.py run --pdb md_prep/system.pdb --backend openmm --ps 20000 --outdir md_out

# 蛋白-配体复合物 MD（全自动：OpenFF Sage2.2 配体 + amber14 蛋白 + TIP3P/0.15M NaCl）
python scripts/np_md.py complex --pdb md_prep/system.pdb --smiles "C1=...槲皮素" --ps 20000 --outdir md_cx
python scripts/np_md.py complex --pdb md_prep/system.pdb --ligand ligand.sdf --ps 100 --outdir md_cx

# 分析（RMSD/势能诊断 SVG+summary：漂移判定"稳定/仍在漂移"）
python scripts/np_md.py analyze --dir md_out --backend openmm

# 本机基准（写回 md_environment.json）
python scripts/np_md.py bench --backend openmm --pdb md_prep/system.pdb --outdir .
```

## 配体参数化链路（2026-11-20 全部打通，关键坑已修复）
| 组件 | 安装来源 | 状态 |
|---|---|---|
| openff-toolkit 0.18.0 | PyPI 直链（0.18 被 yank，唯一可用版） | ✅ |
| openff-units/utilities/forcefields/interchange/nagl-models | conda-forge（PyPI 该系列已清理） | ✅ |
| OpenMMForceFields 0.16（SystemGenerator） | **GitHub 官方源** `git+https://github.com/openmm/openmmforcefields.git` | ✅ |
| am1bcc NAGL 模型 bug | toolkit 0.18 传 "am1bcc" 裸名给 get_model → BadFileSuffixError；np_md.py 运行时补丁映射到本地 openff-gnn-am1bcc-1.0.0.pt | ✅ 已修 |
| micromamba pkgs 缓存 | 默认 C: 缓存损坏（"directory is not empty"）→ 用 `CONDA_PKGS_DIRS=D:/MDtools/mamba/pkgs2` | ✅ |

实测：槲皮素（SMILES）→ Sage2.2 am1bcc → SystemGenerator 32 粒子参数化成功。

## MD 环境长期记忆位置
- skill 内：`data/md_environment.json`（np_md.py check 读取并打印）
- 原始工作区：`D:/MDwork/`（20ns 产物+工作报告）、`D:/MDtools/`（micromamba/envs/WSL）

## 与高分论文的对接（13/23 篇 MD 套路）
- catechin 式验证：对接复合物 → 20-100ns → **RMSD 稳定平台期 + 氢键数 + RMSF** →
  结论"结合稳定性"支撑对接分数。np_md.py analyze 输出 RMSD 诊断；氢键/RMSF 分析
  用 gmx 后端（`gmx hbond/rmsf`）或 OpenMM 轨迹后处理（mdtraj）。
- 蛋白-配体复合物 MD：配体拓扑需 GAFF2（acpype/AmberTools，WSL 可装）或 CGenFF；
  np_md.py 不臆造配体参数——给出接入点（prep --ligand 预留 + 文档路径）。
- 膜蛋白（GPCR 如 DP1）：CHARMM-GUI 网页搭膜体系 → 下载的 system 含拓扑 → OpenMM
  可直接跑（forcefield 用 charmm36.xml+charmm36/water.xml）。

## 双后端选择决策树
1. 本沙箱/Windows 原生 → **OpenMM**（唯一可自动化路径）
2. 20ns+ 正式生产且用户放行 wsl.exe → GROMACS@WSL（mdp 参数模板已内置 np_md.py）
3. 100ns+/多副本 → 学校集群（gmx），本机只做 prep+分析
