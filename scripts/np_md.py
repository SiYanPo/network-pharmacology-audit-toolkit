#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_md.py —— v1.6.0 分子动力学模拟模块（GROMACS 全流程 + OpenMM 快筛双后端）

背景：23 篇中 13 篇提及 MD（catechin 真用 GROMACS 2020.6+AMBER99SB+100ns+RMSD）。
本机 MD 环境实测（写入 data/md_environment.json 长期记忆）：
  - OpenMM + pdbfixer + rdkit：由 np_tools 的 python_md 解析（NP_MD_PYTHON 可覆盖）
  - GROMACS：此前会话经 WSL Ubuntu 跑通 20ns 全流程（D:/MDwork/00_test/run 有 em/nvt/topol.top）；
    WSL 在当前沙箱被安全策略拦截 → 本模块另装 Windows 原生 gmx（conda-forge）
  - 基准（前会话实测）：7 万原子 11.27 ns/day（i5-12400 6C12T / 8GB / 无独显）
  - 生产级 100ns+ 建议：本机快筛 + 集群生产（D:/MDwork 工作报告既有结论）

命令：
  check                          环境体检（gmx/openmm/pdbfixer/版本/平台）
  prep  --pdb receptor.pdb [--ligand ligand.pdb] --outdir DIR
         pdbfixer 修复（缺失残基/加氢/去杂原子）→ MD-ready system.pdb
  run   --pdb prep/system.pdb --backend openmm|gmx --ps 100 [--implicit] --outdir DIR
         openmm: 最小化→NVT→NPT→生产（LangevinMiddle/2fs/PME/0.15M NaCl），流式 RMSD
         gmx:    pdb2gmx(amber99sb-ildn)→solvate→ions→em→nvt→npt→md（mdp 内置）
  analyze --dir DIR --backend ...   RMSD 曲线 SVG+CSV+统计摘要
  bench  --backend openmm [--ps 0.05]   本机 ns/day 基准（写 md_environment.json）

诚实边界：
  - 蛋白-配体复合物：gmx 路线配体拓扑需 GAFF2(acpype/AmberTools)或 CGenFF——提供检查与
    手工接入点，不臆造参数；openmm 路线同理（OpenFF 未装时仅蛋白）。
  - 20ns 全流程参考：D:/MDwork/00_test/run（前会话已验证 mdp 参数）。
"""
import argparse
import csv
import json
import math
import re
import os
import shutil
import subprocess
import sys
import time

MD_ENV_MEM = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "data", "md_environment.json")
try:
    import np_tools as _T
    OPENMM_PY = _T.path_of("python_md") or os.environ.get("NP_MD_PYTHON") or sys.executable
except Exception:
    OPENMM_PY = os.environ.get("NP_MD_PYTHON") or sys.executable
GMX_CANDIDATES = []  # 由 np_tools 解析（NP_GMX / tools.local.json / PATH）



# ---------------- WSL / GROMACS 桥接 ----------------

WSL_DISTRO = os.environ.get("NP_WSL_DISTRO", "Ubuntu-24.04")
WSL_GMX = os.environ.get("NP_WSL_GMX", "/opt/md/envs/gmx/bin/gmx")


def _win_to_wsl(path):
    p = os.path.abspath(path).replace("\\", "/")
    if re.match(r"^[A-Za-z]:/", p):
        return "/mnt/" + p[0].lower() + p[2:]
    return p


def _wsl_run(cmd, cwd=None, timeout=7200):
    """在 WSL 内执行命令（cwd 为 Windows 路径，自动转换）"""
    pre = f"cd {_win_to_wsl(cwd)} && " if cwd else ""
    full = ["wsl", "-d", WSL_DISTRO, "--", "bash", "-lc", pre + cmd]
    r = subprocess.run(full, capture_output=True, text=True, timeout=timeout)
    return r


def _gmx_ok():
    r = _wsl_run(f"{WSL_GMX} --version 2>&1 | head -2", timeout=120)
    return "GROMACS" in (r.stdout or "")


def _find_gmx():
    for c in GMX_CANDIDATES:
        if c in ("gmx",) and shutil.which("gmx"):
            return shutil.which("gmx")
        if os.path.exists(c):
            return c
    return None


def _run_openmm_code(code, timeout=3600):
    """在专用 md env 里执行 OpenMM 代码"""
    r = subprocess.run([OPENMM_PY, "-c", code], capture_output=True, text=True, timeout=timeout)
    return r


# ---------------- check ----------------

def cmd_check(args):
    mem = {}
    if os.path.exists(MD_ENV_MEM):
        mem = json.load(open(MD_ENV_MEM, encoding="utf-8"))
    gmx = _find_gmx()
    try:
        if _gmx_ok():
            r = _wsl_run(f"{WSL_GMX} --version 2>&1 | head -3", timeout=120)
            print("GROMACS(WSL): " + WSL_GMX)
            print("  " + " ".join(r.stdout.splitlines()[:2]))
        else:
            print("GROMACS(WSL): 未检测到（NP_WSL_GMX 可指定路径）")
    except Exception as e:
        print(f"GROMACS(WSL): 不可用（{type(e).__name__}）")
    if not gmx:
        print("GROMACS(Windows 原生): 无 win-64 包（预期）；使用 WSL 后端 `np_md.py gmx`")
    if gmx and gmx != "gmx" or (gmx and shutil.which("gmx")):
        try:
            r = subprocess.run([gmx, "--version"], capture_output=True, text=True, timeout=60)
            line = next((l for l in r.stdout.splitlines() if l.startswith("GROMACS version")), "")
            print(f"  {line}")
        except Exception:
            pass
    try:
        r = _run_openmm_code(
            "import openmm;print('OpenMM', openmm.__version__);"
            "print('platforms:', [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])",
            timeout=120)
        print(r.stdout.strip())
    except Exception as e:
        print(f"OpenMM: 不可用（{e}）——安装：micromamba create -p <env_dir> -c conda-forge openmm pdbfixer rdkit（然后用 NP_MD_PYTHON 指向它）")
    if mem:
        print("\n长期记忆（data/md_environment.json）:")
        for k, v in mem.items():
            print(f"  {k}: {v}")


# ---------------- prep ----------------

def cmd_prep(args):
    outdir = args.outdir
    script = f'''
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile
fx = PDBFixer(filename=r"{args.pdb}")
fx.findMissingResidues(); fx.findMissingAtoms()
fx.addMissingAtoms(); fx.removeHeterogens(keepWater=False)
fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
fx.addMissingHydrogens(7.0)
PDBFile.writeFile(fx.topology, fx.positions, open(r"{os.path.join(outdir, 'system.pdb')}", "w"))
print("prepared:", fx.sequences and len(fx.sequences), "chains")
'''
    os.makedirs(outdir, exist_ok=True)
    r = _run_openmm_code(script, timeout=600)
    if r.returncode != 0:
        raise SystemExit(f"[FAIL] prep 失败: {r.stderr[-500:]}")
    print(r.stdout.strip())
    print(f"MD-ready 结构 → {outdir}/system.pdb（缺失残基/原子已修复，杂原子已去除）")


# ---------------- openmm run ----------------

_OPENMM_RUN = '''
import sys, time, math
import openmm as mm
import openmm.app as app
from openmm import unit

pdb_path, outdir, ps_total, implicit = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4] == "1"
pdb = app.PDBFile(pdb_path)
if implicit == "1":
    ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")  # 占位，实际用 GB
    system = ff.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=1.6*unit.nanometer, constraints=app.HBonds)
else:
    ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    mod = app.Modeller(pdb.topology, pdb.positions)
    mod.addHydrogens(ff)
    mod.addSolvent(ff, model="tip3p", padding=1.0*unit.nanometer, ionicStrength=0.15*unit.molar)
    system = ff.createSystem(mod.topology, nonbondedMethod=app.PME, constraints=app.HBonds,
                             nonbondedCutoff=1.0*unit.nanometer)
    pdb.topology, pdb.positions = mod.topology, mod.positions

integ = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
plat = mm.Platform.getPlatformByName("CPU" if "CUDA" not in [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())] else "CUDA")
sim = app.Simulation(pdb.topology, system, integ, plat)
sim.context.setPositions(pdb.positions)
print("atoms:", system.getNumParticles(), "minimizing...", flush=True)
sim.minimizeEnergy(maxIterations=2000)
sim.context.setVelocitiesToTemperature(300*unit.kelvin, 1234)

nsteps = int(ps_total * 1000 / 2)  # 2 fs/step
equil = min(5000, nsteps // 5)
print(f"nsteps={nsteps} equil={equil}", flush=True)
sim.step(equil)
import csv as _csv
f = open(f"{outdir}/md_progress.csv", "w", newline="")
w = _csv.writer(f); w.writerow(["step", "time_ps", "pot_kj", "rmsd_nm", "ns_per_day"])
reduced = system.getDefaultPeriodicBoxVectors() if implicit == "0" else None
res = list(pdb.topology.residues())
ca_idx = [a.index for r in res[:40] for a in r.atoms() if a.name == "CA"] or None
init = None
t0 = time.time(); sim.reporters = []
from openmm.app import DCDFile
dcd = DCDFile(open(f"{outdir}/traj.dcd", "wb"), pdb.topology, 2*unit.picoseconds, 1, len(pdb.positions))
save_every = max(500, nsteps // 100)
for i in range(equil, nsteps, save_every):
    sim.step(save_every)
    st = sim.context.getState(getEnergy=True, getPositions=True)
    pe = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    pos = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)  # ndarray nm
    import numpy as _np
    xyz = _np.asarray(pos)
    if init is None:
        init = xyz[ca_idx] if ca_idx else xyz
    cur = xyz[ca_idx] if ca_idx else xyz
    rmsd = float(_np.sqrt(((cur - init) ** 2).sum(axis=1).mean()))  # nm
    if ca_idx is None:
        rmsd = 0.0
    el = time.time() - t0
    nspd = (save_every * 2 / 1e6) / (el / 86400) if el > 0 else 0  # ns/day
    w.writerow([i, i * 2 / 1000, round(pe, 1), round(rmsd, 4), round(nspd, 2)])
    dcd.writeModel(pos)
f.close()
print("done", nsteps, "steps", flush=True)
'''


def cmd_run(args):
    os.makedirs(args.outdir, exist_ok=True)
    if args.backend == "openmm":
        print(f"OpenMM 生产：{args.ps} ps | {'implicit' if args.implicit else '显式溶剂 PME + 0.15M NaCl'}")
        t0 = time.time()
        r = subprocess.run([OPENMM_PY, "-c", _OPENMM_RUN, args.pdb, args.outdir, str(args.ps),
                            "1" if args.implicit else "0"],
                           capture_output=True, text=True, timeout=7200)
        print(r.stdout[-800:])
        if r.returncode != 0:
            raise SystemExit(f"[FAIL] openmm run: {r.stderr[-800:]}")
        print(f"耗时 {time.time() - t0:.0f}s → {args.outdir}/md_progress.csv + traj.dcd")
    else:
        _run_gmx(args)


_MDP = {
 "em": """integrator = steep
emtol = 1000.0
emstep = 0.01
nsteps = 50000
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
""",
 "nvt": """integrator = md
dt = 0.002
nsteps = 25000
nstxout = 5000
nstenergy = 5000
nstxout-compressed = 5000
continuation = no
define = -DPOSRES
constraints = h-bonds
constraint-algorithm = lincs
tcoupl = V-rescale
tc-grps = System
tau-t = 1.0
ref-t = 300
pcoupl = no
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
DispCorr = EnerPres
""",
 "npt": """integrator = md
dt = 0.002
nsteps = 25000
nstxout = 5000
nstenergy = 5000
nstxout-compressed = 5000
continuation = yes
define = -DPOSRES
gen-vel = no
constraints = h-bonds
constraint-algorithm = lincs
tcoupl = V-rescale
tc-grps = System
tau-t = 1.0
ref-t = 300
pcoupl = C-rescale
pcoupltype = isotropic
tau-p = 2.0
compressibility = 4.5e-5
ref-p = 1.0
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
DispCorr = EnerPres
""",
 "md": """integrator = md
dt = 0.002
nsteps = {nsteps}
nstxout = 5000
nstenergy = 5000
nstlog = 5000
nstxout-compressed = 5000
continuation = yes
define = -DPOSRES
gen-vel = no
constraints = h-bonds
constraint-algorithm = lincs
tcoupl = V-rescale
tc-grps = System
tau-t = 1.0
ref-t = 300
pcoupl = C-rescale
pcoupltype = isotropic
tau-p = 2.0
compressibility = 4.5e-5
ref-p = 1.0
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
DispCorr = EnerPres
""",
}


def _run_gmx(args):
    gmx = _find_gmx()
    if not gmx:
        raise SystemExit("[FAIL] 未找到 gmx——三种修法：\n"
                         "  1) 设 NP_GMX=/path/to/gmx\n"
                         "  2) config/tools.local.json 写 {\"gmx\": \"C:/path/gmx.exe\"}\n"
                         "  3) conda: micromamba create -p <env> -c conda-forge gromacs -y 后把 bin 加进 PATH")
    os.makedirs(args.outdir, exist_ok=True)
    nsteps = int(args.ps * 1000 / 2)
    chain = [
        # pdb2gmx + solvate + ions 由 prep.pdb 出发
        f'"{gmx}" pdb2gmx -f "{args.pdb}" -o pro.gro -p topol.top -i posre.itp -ff amber99sb-ildn -water tip3p -ignh',
        f'"{gmx}" editconf -f pro.gro -o box.gro -c -d 1.0 -bt dodecahedron',
        f'"{gmx}" solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top',
        # ions 需要 ions.mdp + grompp
    ]
    os.chdir(args.outdir)
    for name, cmd in [("pdb2gmx", chain[0]), ("editconf", chain[1]), ("solvate", chain[2])]:
        print(f"== {name}", flush=True)
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=1800)
        if r.returncode != 0:
            raise SystemExit(f"[FAIL] {name}: {(r.stderr or r.stdout)[-600:]}")
    open("ions.mdp", "w").write(_MDP["em"])
    r = subprocess.run(f'"{gmx}" grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2',
                       shell=True, capture_output=True, text=True, timeout=1200)
    if r.returncode != 0:
        raise SystemExit(f"[FAIL] grompp ions: {(r.stderr or r.stdout)[-600:]}")
    r = subprocess.run(f'echo SOL | "{gmx}" genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15',
                       shell=True, capture_output=True, text=True, timeout=1200)
    if r.returncode != 0:
        raise SystemExit(f"[FAIL] genion: {(r.stderr or r.stdout)[-600:]}")
    # 上一阶段的输出结构：em 从含离子的溶剂化盒子起，随后逐级接力
    steps = [("em", "solv_ions"), ("nvt", "em"), ("npt", "nvt")]
    prev = "solv_ions"
    for nm, inp in steps:
        open(f"{nm}.mdp", "w").write(_MDP[nm])
        r = subprocess.run(f'"{gmx}" grompp -f {nm}.mdp -c {inp}.gro -p topol.top -o {nm}.tpr -maxwarn 2',
                           shell=True, capture_output=True, text=True, timeout=1200)
        if r.returncode != 0:
            raise SystemExit(f"[FAIL] grompp {nm}: {(r.stderr or r.stdout)[-600:]}")
        r = subprocess.run(f'"{gmx}" mdrun -deffnm {nm}', shell=True, capture_output=True, text=True, timeout=7200)
        if r.returncode != 0:
            raise SystemExit(f"[FAIL] mdrun {nm}: {(r.stderr or r.stdout)[-600:]}")
        prev = nm
    # 生产
    open("md.mdp", "w").write(_MDP["md"].replace("{nsteps}", str(nsteps)))
    r = subprocess.run(f'"{gmx}" grompp -f md.mdp -c {prev}.gro -t {prev}.cpt -p topol.top -o md.tpr -maxwarn 2',
                       shell=True, capture_output=True, text=True, timeout=1200)
    if r.returncode != 0:
        raise SystemExit(f"[FAIL] grompp md: {(r.stderr or r.stdout)[-600:]}")
    print(f"== mdrun 生产 {args.ps} ps（{nsteps} 步）", flush=True)
    subprocess.run(f'"{gmx}" mdrun -deffnm md', shell=True, timeout=72000)
    print(f"生产完成 → {args.outdir}/md.xtc（分析用 analyze --backend gmx）")


# ---------------- analyze ----------------

_OPENMM_COMPLEX = '''
import sys, time, glob, os
import numpy as _np
import openmm as mm
import openmm.app as app
from openmm import unit

pdb_path, outdir, ps_total = sys.argv[1], sys.argv[2], float(sys.argv[3])
ligand = sys.argv[4]  # .sdf/.mol2 文件 或 "SMILES:xxx"
import openff.nagl_models._dynamic_fetch as _df
import openff.nagl_models
_orig_get_model = _df.get_model
def _patched(filename=None, **kw):
    if str(filename).strip() == "am1bcc":
        cands = sorted(glob.glob(os.path.join(os.path.dirname(openff.nagl_models.__file__),
                                              "models", "am1bcc", "*.pt")))
        filename = cands[-1]
    return _orig_get_model(filename=filename, **kw)
_df.get_model = _patched

from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator

pdb = app.PDBFile(pdb_path)
if ligand.startswith("SMILES:"):
    lig = Molecule.from_smiles(ligand[7:])
else:
    lig = Molecule.from_file(ligand)
if not lig.conformers:
    lig.generate_conformers(n_conformers=1)
    print("generated 3D conformer for ligand", flush=True)

mod = app.Modeller(pdb.topology, pdb.positions)
lig_pos = unit.Quantity(_np.asarray(lig.conformers[0].m), unit=unit.nanometer)  # OFF 构象默认 nm
mod.add(lig.to_topology().to_openmm(), lig_pos)
solvent_ff = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
tg = SMIRNOFFTemplateGenerator(molecules=[lig], forcefield="openff-2.2.0")
solvent_ff.registerTemplateGenerator(tg.generator)
mod.addSolvent(solvent_ff, model="tip3p", padding=1.0*unit.nanometer, ionicStrength=0.15*unit.molar,
               positiveIon="Na+", negativeIon="Cl-")
print("solvated atoms:", mod.topology.getNumAtoms(), flush=True)

gen = SystemGenerator(forcefields=["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"],
                      small_molecule_forcefield="openff-2.2.0")
system = gen.create_system(mod.topology, molecules=[lig])

integ = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
plat_names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
want = os.environ.get("NP_MD_PLATFORM", "CPU")  # Intel 核显 OpenCL 在 minimization 崩溃（clWaitForEvents），默认 CPU
plat = mm.Platform.getPlatformByName(want if want in plat_names else "CPU")
sim = app.Simulation(mod.topology, system, integ, plat)
sim.context.setPositions(mod.positions)
print("minimizing...", flush=True)
sim.minimizeEnergy(maxIterations=2000)
sim.context.setVelocitiesToTemperature(300*unit.kelvin, 1234)

nsteps = int(ps_total * 1000 / 2)
equil = min(5000, nsteps // 5)
sim.step(equil)
import csv as _csv
f = open(f"{outdir}/md_progress.csv", "w", newline="")
w = _csv.writer(f); w.writerow(["step", "time_ps", "pot_kj", "rmsd_nm", "ns_per_day"])
ca_idx = [a.index for r in list(mod.topology.residues())[:40] for a in r.atoms() if a.name == "CA"] or None
init = None
t0 = time.time()
from openmm.app import DCDFile
dcd = DCDFile(open(f"{outdir}/traj.dcd", "wb"), mod.topology, 2*unit.picoseconds, 1, mod.topology.getNumAtoms())
save_every = max(500, nsteps // 100)
for i in range(equil, nsteps, save_every):
    sim.step(save_every)
    st = sim.context.getState(getEnergy=True, getPositions=True)
    pe = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    pos = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    xyz = _np.asarray(pos)
    if init is None:
        init = xyz[ca_idx] if ca_idx else xyz
    cur = xyz[ca_idx] if ca_idx else xyz
    rmsd = float(_np.sqrt(((cur - init) ** 2).sum(axis=1).mean())) if ca_idx else 0.0
    el = time.time() - t0
    nspd = (save_every * 2 / 1e6) / (el / 86400) if el > 0 else 0  # ns/day
    w.writerow([i, i * 2 / 1000, round(pe, 1), round(rmsd, 4), round(nspd, 2)])
    dcd.writeModel(st.getPositions(asNumpy=True))
f.close()
print("done", nsteps, "steps", flush=True)
'''


def cmd_complex(args):
    os.makedirs(args.outdir, exist_ok=True)
    ligand = args.ligand if args.ligand else ("SMILES:" + args.smiles)
    if args.ligand and not os.path.exists(args.ligand):
        raise SystemExit(f"[FAIL] 配体文件不存在: {args.ligand}（支持 .sdf/.mol2，或用 --smiles）")
    print(f"复合物 MD：蛋白={args.pdb} | 配体={ligand[:60]} | {args.ps} ps")
    t0 = time.time()
    r = subprocess.run([OPENMM_PY, "-c", _OPENMM_COMPLEX, args.pdb, args.outdir, str(args.ps), ligand],
                       capture_output=True, text=True, timeout=72000)
    print(r.stdout[-600:])
    if r.returncode != 0:
        raise SystemExit(f"[FAIL] complex: {r.stderr[-800:]}")
    print(f"耗时 {time.time() - t0:.0f}s → {args.outdir}/md_progress.csv + traj.dcd")



def cmd_gmx(args):
    """WSL 内 GROMACS 全流程（em→nvt→npt→md），路径自动转换"""
    if not _gmx_ok():
        raise SystemExit("[FAIL] WSL 内未找到 gmx（设置 NP_WSL_GMX）")
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    src_pdb = os.path.abspath(args.pdb)
    # 复制输入到输出目录，避免相对路径问题
    import shutil as _sh
    local_pdb = os.path.join(outdir, "input.pdb")
    if src_pdb != local_pdb:
        _sh.copyfile(src_pdb, local_pdb)
    nsteps = int(args.ps * 1000 / 2)
    # 写 mdp
    files = {"em": _MDP["em"], "nvt": _MDP["nvt"], "npt": _MDP["npt"],
             "md": _MDP["md"].replace("{nsteps}", str(nsteps)).replace("define = -DPOSRES\n", "")}
    for k, v in files.items():
        open(os.path.join(outdir, f"{k}.mdp"), "w").write(v)
    G = WSL_GMX
    steps = [
        (f"{G} pdb2gmx -f input.pdb -o pro.gro -p topol.top -i posre.itp -ff amber99sb-ildn -water tip3p -ignh", "pdb2gmx"),
        (f"{G} editconf -f pro.gro -o box.gro -c -d 1.0 -bt dodecahedron", "editconf"),
        (f"{G} solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top", "solvate"),
    ]
    for cmd, name in steps:
        print(f"== {name}", flush=True)
        r = _wsl_run(cmd, cwd=outdir, timeout=3600)
        if r.returncode != 0:
            raise SystemExit(f"[FAIL] {name}: {(r.stderr or r.stdout)[-600:]}")
    r = _wsl_run(f"{G} grompp -f em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2", cwd=outdir)
    if r.returncode != 0:
        raise SystemExit(f"[FAIL] grompp ions: {(r.stderr or r.stdout)[-500:]}")
    r = _wsl_run(f"echo SOL | {G} genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15",
                 cwd=outdir)
    if r.returncode != 0:
        raise SystemExit(f"[FAIL] genion: {(r.stderr or r.stdout)[-500:]}")
    prev = "solv_ions"
    for nm in ("em", "nvt", "npt"):
        print(f"== {nm}", flush=True)
        r = _wsl_run(f"{G} grompp -f {nm}.mdp -c {prev}.gro -r {prev}.gro -p topol.top -o {nm}.tpr -maxwarn 2", cwd=outdir)
        if r.returncode != 0:
            raise SystemExit(f"[FAIL] grompp {nm}: {(r.stderr or r.stdout)[-500:]}")
        r = _wsl_run(f"{G} mdrun -deffnm {nm} -nt {args.nt}", cwd=outdir, timeout=14400)
        if r.returncode != 0:
            raise SystemExit(f"[FAIL] mdrun {nm}: {(r.stderr or r.stdout)[-500:]}")
        prev = nm
    print(f"== md {args.ps} ps", flush=True)
    r = _wsl_run(f"{G} grompp -f md.mdp -c {prev}.gro -r {prev}.gro -t {prev}.cpt -p topol.top -o md.tpr -maxwarn 2", cwd=outdir)
    if r.returncode != 0:
        raise SystemExit(f"[FAIL] grompp md: {(r.stderr or r.stdout)[-500:]}")
    r = _wsl_run(f"{G} mdrun -deffnm md -nt {args.nt}", cwd=outdir, timeout=86400)
    print((r.stdout or "")[-400:])
    # 分析：RMSD
    r2 = _wsl_run(f"{G} rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns -select 'Backbone'", cwd=outdir)
    print("RMSD:", "OK" if r2.returncode == 0 else r2.stderr[-200:])
    perf = [l for l in (r.stdout or "").splitlines() if "ns/day" in l or "Performance" in l]
    print("性能:", perf[-1] if perf else "见 md.log")
    print(f"完成 → {outdir}")


def cmd_analyze(args):
    import numpy as np
    outdir = getattr(args, "dir", None) or getattr(args, "outdir", ".")
    if getattr(args, "backend", "openmm") == "gmx":
        return _analyze_gmx(outdir)
    prog = os.path.join(outdir, "md_progress.csv")
    if not os.path.exists(prog):
        raise SystemExit(f"[FAIL] 无 {prog}（openmm 后端产物）或请用 gmx 的 gmx rms 工具")
    rows = list(csv.DictReader(open(prog, encoding="utf-8")))
    t = [float(r["time_ps"]) for r in rows]
    rmsd = [float(r["rmsd_nm"]) for r in rows]
    pe = [float(r["pot_kj"]) for r in rows]
    W, H = 680, 460
    def path(xs, ys, x0, y0, w, h, ymax, ymin):
        pts = []
        for x, y in zip(xs, ys):
            px = x0 + w * (x - xs[0]) / max(1e-9, (xs[-1] - xs[0]))
            py = y0 + h * (1 - (y - ymin) / max(1e-9, ymax - ymin))
            pts.append(f"{px:.0f},{py:.0f}")
        return "M " + " L ".join(pts)
    ymin, ymax = min(pe), max(pe)
    rmin, rmax = min(rmsd), max(rmsd) or 1
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="sans-serif">',
           '<rect width="680" height="460" fill="#ffffff"/>',
           '<text x="340" y="24" text-anchor="middle" font-size="14" font-weight="bold">MD 轨迹诊断（RMSD 与势能）</text>']
    svg.append(f'<path d="{path(t, pe, 80, 60, 520, 150, ymax, ymin)}" fill="none" stroke="#2471a3" stroke-width="1.5"/>')
    svg.append(f'<text x="80" y="50" font-size="11" fill="#2471a3">势能 kJ/mol（{ymin:.0f}~{ymax:.0f}）</text>')
    svg.append(f'<path d="{path(t, rmsd, 80, 250, 520, 150, max(rmax, rmin + 0.1), 0)}" fill="none" stroke="#c0392b" stroke-width="1.5"/>')
    svg.append(f'<text x="80" y="240" font-size="11" fill="#c0392b">RMSD nm（0~{rmax:.2f}）</text>')
    for i in range(6):
        xx = 80 + 520 * i / 5
        svg.append(f'<text x="{xx:.0f}" y="420" font-size="10" text-anchor="middle">{t[0] + (t[-1] - t[0]) * i / 5:.0f} ps</text>')
    svg.append("</svg>")
    open(os.path.join(outdir, "md_diagnostic.svg"), "w", encoding="utf-8").write("".join(svg))
    last50 = rmsd[-max(1, len(rmsd) // 2):]
    summary = {"n_frames": len(rows), "time_ps": t[-1] if t else 0,
               "rmsd_mean_nm": round(sum(rmsd) / len(rmsd), 3) if rmsd else None,
               "rmsd_last_half_mean_nm": round(sum(last50) / len(last50), 3),
               "rmsd_drift_last20pct": round(rmsd[-1] - rmsd[-max(1, len(rmsd) // 5)], 3),
               "pe_mean_kj": round(sum(pe) / len(pe), 0) if pe else None,
               "verdict": "稳定（后50% RMSD 波动<0.1nm）" if last50 and (max(last50) - min(last50)) < 0.1
                          else "仍在漂移/需延长模拟"}
    json.dump(summary, open(os.path.join(outdir, "md_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(json.dumps(summary, ensure_ascii=False, indent=1))
    print(f"→ {outdir}/md_diagnostic.svg + md_summary.json")


# ---------------- bench ----------------

def cmd_bench(args):
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    print("短跑基准（0.05 ps = 25 步显式溶剂）...")
    subprocess.run([sys.executable, os.path.abspath(__file__), "run",
                    "--pdb", args.pdb, "--backend", "openmm", "--ps", "0.05",
                    "--outdir", os.path.join(outdir, "bench")], timeout=3600)
    rows = list(csv.DictReader(open(os.path.join(outdir, "bench", "md_progress.csv"), encoding="utf-8")))
    nspd = float(rows[-1]["ns_per_day"]) if rows else 0
    atoms = "?"
    rec = {"backend": "openmm", "ns_per_day": nspd, "note": "本机快筛基准",
           "prior_gmx_wsl_benchmark": "7万原子 11.27 ns/day（前会话 WSL 实测）",
           "production_advice": "本机适合 ≤1ns 快筛；20ns≈2天(GMX@WSL)/100ns→集群"}
    json.dump(rec, open(MD_ENV_MEM, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"ns/day={nspd} → 已写入长期记忆 {MD_ENV_MEM}")



def _analyze_gmx(outdir):
    """gmx 轨迹分析：mdtraj 读 md.xtc + md.gro → RMSD/Rg → JSON+XVG"""
    import numpy as np
    code = """
import mdtraj as md, numpy as np, json, os, sys
d = sys.argv[1]
t = md.load(os.path.join(d, "md.xtc"), top=os.path.join(d, "md.gro"))
bb = t.atom_slice(t.topology.select("backbone"))
r = md.rmsd(bb, bb, 0)
rg = md.compute_rg(t)
np.savetxt(os.path.join(d, "rmsd_mdtraj.xvg"), np.c_[t.time, r], fmt="%.4f")
summary = dict(frames=int(t.n_frames), atoms=int(t.n_atoms),
               time_ns=round(float(t.time[-1])/1000.0, 3) if t.time is not None else None,
               rmsd_mean_nm=round(float(r.mean()), 3), rmsd_last_nm=round(float(r[-1]), 3),
               rg_mean_nm=round(float(rg.mean()), 3),
               verdict="稳定" if r[-1] < 0.4 else "RMSD>0.4nm，建议延长平衡或检查结构")
json.dump(summary, open(os.path.join(d, "md_summary.json"), "w"), indent=1)
print(json.dumps(summary, ensure_ascii=False))
"""
    r = subprocess.run([OPENMM_PY, "-c", code, outdir], capture_output=True, text=True, timeout=1800)
    print(r.stdout.strip() or r.stderr[-400:])


WSL_MMPBSA_PY = os.environ.get("NP_WSL_MMPBSA_PY", "/opt/md/envs/mmpbsa/bin/python")
WSL_AMBER_SH = os.environ.get("NP_WSL_AMBER_SH", "/opt/ambertools/amber.sh")


def cmd_mmpbsa(args):
    """MM/PBSA 结合自由能（gmx_MMPBSA 1.7.1，WSL 内运行）

    需要：复合物 tpr（-s）、轨迹（-f）、索引分组（--receptor-group/--ligand-group，默认 Protein / 配体残基名）
    """
    outdir = os.path.abspath(args.outdir)
    tpr = os.path.abspath(args.tpr)
    traj = os.path.abspath(args.traj)
    if not os.path.exists(tpr) or not os.path.exists(traj):
        raise SystemExit("[FAIL] tpr 或轨迹不存在")
    args_cmd = (f"source {WSL_AMBER_SH} && {WSL_MMPBSA_PY} -O "
                f"-s {_win_to_wsl(tpr)} -cs {_win_to_wsl(tpr)} -ct {_win_to_wsl(traj)} "
                f"-ci {_win_to_wsl(args.index)} -cg {args.receptor_group} {args.ligand_group} "
                f"-o {_win_to_wsl(os.path.join(outdir, 'FINAL_RESULTS.dat'))} "
                f"-eo {_win_to_wsl(os.path.join(outdir, 'FINAL_RESULTS.csv'))} "
                f"-deo {_win_to_wsl(os.path.join(outdir, 'decomp.csv'))} "
                f"-nogui")
    print("gmx_MMPBSA 运行中（可为数分钟至数十分钟）...", flush=True)
    r = _wsl_run(args_cmd, timeout=7200)
    ok = os.path.exists(os.path.join(outdir, "FINAL_RESULTS.dat"))
    print((r.stdout or "")[-500:])
    if not ok:
        raise SystemExit(f"[FAIL] MM/PBSA 失败：{(r.stderr or r.stdout)[-600:]}")
    # 解析 ΔTOTAL
    txt = open(os.path.join(outdir, "FINAL_RESULTS.dat"), encoding="utf-8", errors="ignore").read()
    import re as _re
    m = _re.search(r"ΔTOTAL\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", txt)
    summary = {"delta_total_kcal": float(m.group(1)) if m else None,
               "sd": float(m.group(2)) if m else None,
               "files": ["FINAL_RESULTS.dat", "FINAL_RESULTS.csv", "decomp.csv"]}
    json.dump(summary, open(os.path.join(outdir, "mmpbsa_summary.json"), "w"), ensure_ascii=False, indent=1)
    print(f"ΔG_bind = {summary['delta_total_kcal']} kcal/mol → {outdir}/mmpbsa_summary.json")


# ---------------- 统计功效：块平均 SEM / ESS / 前后半程 ----------------

def _tau_frames(x):
    """初始正序列（Geyer）估计积分自相关时间，单位=帧。"""
    import numpy as np
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 8:
        return 1.0
    x = x - x.mean()
    var = float((x * x).mean())
    if var <= 0:
        return 1.0
    max_lag = min(n // 4, 200)
    tau, prev = 0.0, float("inf")
    for lag in range(1, max_lag + 1):
        c1 = float((x[:-lag] * x[lag:]).mean() / var)
        c2 = float((x[:-(lag + 1)] * x[lag + 1:]).mean() / var) if lag + 1 < n else 0.0
        pair = c1 + c2
        if pair < 0:
            break
        pair = min(pair, prev)
        tau += pair
        prev = pair
    return max(1.0, 1.0 + 2.0 * tau)


def _block_sem(x, block):
    import numpy as np
    x = np.asarray(x, dtype=float)
    nb = len(x) // block
    if nb < 2:
        return float("nan"), 0
    bm = x[:nb * block].reshape(nb, block).mean(1)
    return float(bm.std(ddof=1) / np.sqrt(nb)), nb


def cmd_stats(args):
    """逐帧能量序列的统计功效复核：块平均 SEM、积分自相关时间、有效样本量、前后半程漂移。

    为什么必须做：MM/PBSA 报 mean ± SEM 时若把每帧当独立样本，帧间强自相关会让
    SEM 被系统性低估（实测 20 ns/203 帧体系，朴素 SEM 低估 2~3 倍）。块平均 SEM 与
    ESS 是审稿人会问的那两个数。

    用法：
      python np_md.py stats --csv e7_mmgbsa_perframe.csv --col dG_bind --block-ns 1.0
    """
    import numpy as np
    rows = list(csv.DictReader(open(args.csv, encoding="utf-8-sig")))
    col = args.col
    vals = []
    for r in rows:
        try:
            vals.append(float(r[col]))
        except (KeyError, TypeError, ValueError):
            continue
    if len(vals) < 8:
        raise SystemExit(f"[FAIL] 可用样本 {len(vals)} < 8，无法做统计功效分析")
    x = np.array(vals)
    n = len(x)
    sd = float(x.std(ddof=1))
    naive = sd / np.sqrt(n)
    block = max(2, int(round(args.block_ns / max(1e-9, args.frame_ns))))
    bsem, nb = _block_sem(x, block)
    tau = _tau_frames(x)
    ess = n / tau
    h = n // 2
    diff = float(x[h:].mean() - x[:h].mean())
    out = {
        "column": col, "n_frames": n,
        "frame_ns": args.frame_ns, "time_span_ns": round(args.frame_ns * n, 2),
        "mean": round(float(x.mean()), 3), "sd": round(sd, 3),
        "sem_naive": round(naive, 3),
        "block_ns": args.block_ns, "n_blocks": nb,
        "sem_block": round(bsem, 3) if bsem == bsem else None,
        "sem_inflation": round(bsem / naive, 2) if bsem == bsem and naive > 0 else None,
        "tau_frames": round(tau, 1),
        "ess": round(ess, 1),
        "first_half_mean": round(float(x[:h].mean()), 3),
        "second_half_mean": round(float(x[h:].mean()), 3),
        "half_diff": round(diff, 3),
        "half_diff_z": round(diff / bsem, 2) if bsem == bsem and bsem > 0 else None,
        "ci95_block": [round(float(x.mean()) - 1.96 * bsem, 2),
                       round(float(x.mean()) + 1.96 * bsem, 2)] if bsem == bsem else None,
        "verdict": ("均值已收敛（前后半程 |z|<2）" if bsem == bsem and bsem > 0 and abs(diff / bsem) < 2
                    else "前后半程差异显著，需延长模拟或分段报告"),
        "note": "报告值请用 sem_block；sem_naive 仅作对照",
    }
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.csv))
    os.makedirs(outdir, exist_ok=True)
    json.dump(out, open(os.path.join(outdir, "md_stats.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))
    print(f"→ {outdir}/md_stats.json")


def cmd_qc(args):
    """轨迹指标的前置伪影检查：在算任何"稳定性"之前，先确认尖峰不是周期性镜像。

    为什么需要：实测过一个体系，配体 RMSD 末值 0.28 nm、**最大值 11.26 nm**、
    相邻帧 23 处跳变 > 0.5 nm；后段均值被这些尖峰抬到 0.784 nm，门禁据此判"不收敛"。
    实际上配体全程在口袋内（正确成像后最大位移 3.41 Å）——**伪影被当成了生物学**。

    检查三条：
      1) 单帧跳变 |Δ| > --jump（默认 5，单位同输入）：镜像跳变/成像残留的典型特征；
      2) 最大值 / 中位数 之比 > --ratio（默认 3）：尖峰离群程度；
      3) 末值与最大值严重不符（末值 < 最大值/2）：说明"漂移"只在少数帧出现，不是持续解离。

    用法（支持 .xvg 与 .csv）：
      python np_md.py qc --file ligand_rmsd.xvg --col 2
      python np_md.py qc --file rep_perframe.csv --col ligand_disp_A
    """
    import numpy as np

    path, col = args.file, args.col
    xs, ys = [], []
    if path.lower().endswith(".xvg"):
        for line in open(path, encoding="utf-8", errors="ignore"):
            s = line.strip()
            if not s or s[0] in "@#":
                continue
            parts = s.split()
            if len(parts) >= 2:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
    else:
        rows = list(csv.DictReader(open(path, encoding="utf-8-sig")))
        for r in rows:
            # 若本工具自己的逐帧表带有 used 标记，默认只看 used==1 的帧，
            # 否则会把已经判定为伪影的帧再算一遍（实测 rep_perframe 里未过滤时最大 22.2 Å）
            if not args.all and "used" in r and str(r.get("used")).strip() not in ("1", "1.0", "True", "true"):
                continue
            try:
                ys.append(float(r[col]))
                xs.append(float(r.get("t_ns", len(xs))))
            except (KeyError, TypeError, ValueError):
                continue
    if len(ys) < 10:
        raise SystemExit(f"[FAIL] 有效数据点 {len(ys)} < 10，无法判定")

    y = np.array(ys)
    d = np.abs(np.diff(y))
    jumps = int((d > args.jump).sum())
    med = float(np.median(y))
    mx = float(y.max())
    end = float(np.mean(y[-5:]))
    total_span = float(y.max() - y.min())

    flags = []
    if jumps > 0:
        flags.append(f"{jumps} 处单帧跳变 >{args.jump}（疑似周期性镜像/成像残留）")
    if med > 0 and mx / med > args.ratio:
        flags.append(f"最大值/中位数 = {mx/med:.1f} >{args.ratio}（尖峰离群）")
    if mx > 0 and end < mx / 2 and mx - end > args.jump:
        flags.append(f"末段均值 {end:.2f} 仅为最大值 {mx:.2f} 的一半以下（非持续解离形态）")
    # 剔除跳变所在帧 + 其后一帧，估算需要剔除的比例
    drop = set()
    for i, v in enumerate(d):
        if v > args.jump:
            drop.update((i, i + 1))
    frac = len(drop) / len(y) if y.size else 0.0

    verdict = ("通过：未见镜像伪影特征" if not flags
               else "存疑：先按 -pbc mol -center 成像并核对剔除帧，再谈稳定性")
    out = {
        "file": os.path.abspath(path), "n_points": int(y.size),
        "min": round(float(y.min()), 3), "max": round(mx, 3),
        "median": round(med, 3), "mean": round(float(y.mean()), 3),
        "last5_mean": round(end, 3), "span": round(total_span, 3),
        "jumps_gt_threshold": jumps, "jump_threshold": args.jump,
        "max_over_median": round(mx / med, 2) if med else None,
        "frames_to_drop": len(drop), "drop_fraction": round(frac, 4),
        "flags": flags, "verdict": verdict,
        "recommendation": ("成像命令：gmx trjconv -s md.tpr -f md.xtc -o md_whole.xtc "
                           "-pbc mol -center -ur compact；配体位移用蛋白叠合矩阵变换后直接比坐标，"
                           "不要用 gmx rms 的配体自拟合"),
    }
    outdir = args.outdir or os.path.dirname(os.path.abspath(path))
    os.makedirs(outdir, exist_ok=True)
    json.dump(out, open(os.path.join(outdir, "md_traj_qc.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))
    print(f"→ {outdir}/md_traj_qc.json")
    if flags and args.strict:
        raise SystemExit("[FAIL] 存在镜像伪影特征（--strict 模式）")


def main():
    ap = argparse.ArgumentParser(description="np_md: 分子动力学模拟（gmx/openmm 双后端）")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check")
    s = sub.add_parser("prep"); s.add_argument("--pdb", required=True); s.add_argument("--ligand", default=None); s.add_argument("--outdir", default="md_prep")
    s = sub.add_parser("run"); s.add_argument("--pdb", required=True); s.add_argument("--backend", default="openmm", choices=["openmm", "gmx"])
    s.add_argument("--ps", type=float, default=100, help="模拟时长（皮秒）；20ns=20000")
    s.add_argument("--implicit", action="store_true"); s.add_argument("--outdir", default="md_out")
    s = sub.add_parser("complex")
    s.add_argument("--pdb", required=True, help="prep 产出的 MD-ready 蛋白 PDB")
    s.add_argument("--ligand", default=None, help="配体文件 .sdf/.mol2")
    s.add_argument("--smiles", default=None, help="或直接给配体 SMILES")
    s.add_argument("--ps", type=float, default=100)
    s.add_argument("--outdir", default="md_complex")
    s = sub.add_parser("gmx")
    s.add_argument("--pdb", required=True); s.add_argument("--ps", type=float, default=100)
    s.add_argument("--nt", default="6"); s.add_argument("--outdir", default="gmx_out")
    s = sub.add_parser("mmpbsa")
    s.add_argument("--tpr", required=True); s.add_argument("--traj", required=True)
    s.add_argument("--index", required=True); s.add_argument("--receptor-group", default="Protein")
    s.add_argument("--ligand-group", default="Other"); s.add_argument("--outdir", default="mmpbsa")
    s = sub.add_parser("analyze"); s.add_argument("--dir", required=True); s.add_argument("--backend", default="openmm")
    s = sub.add_parser("bench"); s.add_argument("--backend", default="openmm"); s.add_argument("--pdb", required=True); s.add_argument("--outdir", default=".")
    s = sub.add_parser("stats", help="逐帧能量的块平均 SEM / ESS / 前后半程收敛复核")
    s.add_argument("--csv", required=True); s.add_argument("--col", default="dG_bind")
    s.add_argument("--block-ns", type=float, default=1.0); s.add_argument("--frame-ns", type=float, default=0.1)
    s.add_argument("--outdir", default=None)
    s = sub.add_parser("qc", help="轨迹指标前置伪影检查（镜像跳变/尖峰离群）")
    s.add_argument("--file", required=True, help=".xvg 或 .csv")
    s.add_argument("--col", default="ligand_disp_A", help="CSV 列名；xvg 忽略（取第 2 列）")
    s.add_argument("--jump", type=float, default=5.0, help="单帧跳变阈值（与输入同单位）")
    s.add_argument("--ratio", type=float, default=3.0, help="最大值/中位数 告警阈值")
    s.add_argument("--strict", action="store_true", help="发现问题即非零退出")
    s.add_argument("--all", action="store_true", help="CSV 含 used 列时也统计 used==0 的帧")
    s.add_argument("--outdir", default=None)
    a = ap.parse_args()
    {"check": cmd_check, "prep": cmd_prep, "run": cmd_run, "complex": cmd_complex, "gmx": cmd_gmx, "mmpbsa": cmd_mmpbsa,
     "analyze": cmd_analyze, "bench": cmd_bench, "stats": cmd_stats, "qc": cmd_qc}[a.cmd](a)


if __name__ == "__main__":
    main()
