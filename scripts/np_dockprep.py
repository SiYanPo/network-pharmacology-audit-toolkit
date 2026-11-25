#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_dockprep.py — 真实分子对接：结构准备 → Vina 运行 → 重对接自校验 → 结果表

为什么需要这个阶段
  np_dock.py 只能"评估"别人给的对接结果。问题是：拿到一张 affinity 表时，
  没人知道那个受体的准备是否合理、盒子是否覆盖结合位点、协议能否复原晶体姿态。
  本阶段把"能对接"变成"能自证"：每个靶点先做共晶配体回对接（redocking），
  RMSD ≤ 2.0 Å 才认为该靶点的协议可信，其打分才允许进入结论。

依赖（均可选，缺失时明确报"不可用"而不是静默降级）
  vina         AutoDock Vina 1.2.x 可执行文件（--vina 指定路径）
  meeko        ligand/receptor → PDBQT
  rdkit        构象生成、键级修复、对称校正 RMSD
  gemmi        PDB 解析与拆分

输出（np_dock.py 可直接消费）
  docking_scores.csv      compound,target,affinity,redock_rmsd,protocol,pose_file
  docking_validation.csv  target,pdb,ligand_code,rmsd,passed,note
  dockprep_report.md      人读报告：哪些靶点通过、哪些失败、为什么

用法
  # 1) 只做重对接校验（先确认协议可用）
  python np_dockprep.py --plan plan.tsv --outdir dock/ --validate-only

  # 2) 校验 + 跨对接（用上游成分表的 SMILES）
  python np_dockprep.py --plan plan.tsv --ligands compounds.csv --outdir dock/

plan.tsv 列：target \t pdb \t ligand_code [\t chains]
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

import npcore as C

RCSB_PDB = "https://files.rcsb.org/download/%s.pdb"
RCSB_SDF = "https://files.rcsb.org/ligands/download/%s_ideal.sdf"

# 共晶中常见的"非目标配体"（辅因子/溶剂/离子），选配体时排除
EXCLUDE_HET = {
    "HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "PGE", "MPD", "TRS",
    "ACT", "ACE", "FMT", "CIT", "MES", "EPE", "IMD", "NH4", "CL", "NA", "K", "MG", "CA",
    "ZN", "MN", "FE", "CU", "NI", "CD", "BR", "IOD", "F", "NO3", "BME", "DTT", "GSH",
    "ADP", "ATP", "GDP", "GTP", "AMP", "GMP", "NAD", "NAP", "NDP", "FAD", "FMN", "SAM",
    "SAH", "HEM", "PLP", "TPP", "COA", "ACO", "ANP", "UDP", "C5P", "PRP", "B12",
}


def ascii_scratch():
    """返回一个**纯 ASCII** 的临时工作目录。

    为什么必须这样：RDKit / meeko / Vina 均用 C++ 流打开文件，
    Windows 上遇到非 ASCII 路径（中文用户名、中文目录）会直接失败，
    且报错只说 "Bad input file" 或一段截断的 traceback，完全不提编码。
    中间产物全部走 ASCII 目录，最终结果再拷回用户目录。
    """
    cand = tempfile.mkdtemp(prefix="npdock_")
    if all(ord(c) < 128 for c in cand):
        return cand
    try:
        shutil.rmtree(cand, ignore_errors=True)
    except Exception:
        pass
    root = os.environ.get("SystemDrive", "C:") + os.sep
    alt = os.path.join(root, "npdock_tmp_%d" % (os.getpid() % 100000))
    os.makedirs(alt, exist_ok=True)
    return alt


def find_tool(name, extra=()):
    """在 PATH、当前解释器的 Scripts 目录与常见安装位置找可执行文件。

    为什么需要看 Scripts：Windows 上 pip 安装的 meeko 命令行只落在
    <python>/Scripts 里，而该目录不一定在 PATH 中（本机即如此）。
    """
    roots = []
    try:
        roots.append(os.path.join(os.path.dirname(sys.executable), "Scripts"))
        import sysconfig
        roots.append(sysconfig.get_path("scripts") or "")
    except Exception:
        pass
    roots += ["/usr/local/bin", "/opt/homebrew/bin", os.path.expanduser("~/.local/bin")]
    cands = list(extra) + [name, name + ".exe"]
    for r in roots:
        if r:
            cands += [os.path.join(r, name), os.path.join(r, name + ".exe")]
    for cand in cands:
        if not cand:
            continue
        p = shutil.which(cand) if os.sep not in cand else None
        if p:
            return p
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return os.path.abspath(cand)
    return None


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def fetch(url, dst, timeout=120):
    if os.path.exists(dst) and os.path.getsize(dst) > 200:
        return dst
    import urllib.request
    req = urllib.request.Request(url, headers={"User-Agent": "np-pro/1.2"})
    with urllib.request.urlopen(req, timeout=timeout) as r, open(dst, "wb") as f:
        f.write(r.read())
    return dst


def split_structure(pdb_path, ligand_code, chains, out_prefix):
    """拆出受体（仅聚合物）与共晶配体；返回 (受体 PDB, 配体原子块, 配体理想 SDF 路径)"""
    import gemmi
    st = gemmi.read_structure(pdb_path)
    st.setup_entities()
    m = st[0]

    # AUTO：自动挑选共晶配体（最大且不属于溶剂/离子/辅因子的杂原子残基）
    if ligand_code.upper() == "AUTO":
        cands = []
        for ch in m:
            if chains and ch.name not in chains:
                continue
            for res in ch:
                if res.het_flag != "H":
                    continue
                nm = res.name.strip().upper()
                if nm in EXCLUDE_HET:
                    continue
                nheavy = sum(1 for at in res if at.element.name != "H")
                if 12 <= nheavy <= 90:
                    cands.append((nheavy, nm))
        if not cands:
            raise RuntimeError("结构中找不到可用的共晶配体（AUTO 模式）")
        cands.sort(reverse=True)
        ligand_code = cands[0][1]
        auto_note = "自动选择（重原子 %d，候选 %d 个）" % (cands[0][0], len(cands))
    else:
        auto_note = "指定"

    lig_lines, keep_chain = [], set()
    for ch in m:
        if chains and ch.name not in chains:
            continue
        for res in ch:
            if res.het_flag == "A":
                keep_chain.add(ch.name)
            elif res.het_flag == "H" and res.name.strip().upper() == ligand_code.upper():
                keep_chain.add(ch.name)
                for at in res:
                    if at.element.name == "H":
                        continue
                    lig_lines.append(
                        "HETATM%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
                        % (len(lig_lines) + 1, at.name, res.name, ch.name,
                           res.seqid.num, at.pos.x, at.pos.y, at.pos.z, at.element.name))
    if not lig_lines:
        raise RuntimeError("结构 %s 中未找到配体 %s" % (os.path.basename(pdb_path), ligand_code))

    rec = gemmi.Structure()
    nm = gemmi.Model("1")
    for ch in m:
        if ch.name not in keep_chain:
            continue
        nc = gemmi.Chain(ch.name)
        for res in ch:
            if res.het_flag == "A":       # 仅保留氨基酸残基
                nc.add_residue(res)
        if len(nc):
            nm.add_chain(nc)
    rec.add_model(nm)
    rec_path = out_prefix + "_receptor.pdb"
    rec.write_pdb(rec_path)

    lig_path = out_prefix + "_ligand_block.pdb"
    with open(lig_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lig_lines) + "\nEND\n")
    nres = sum(len(ch) for ch in nm)      # 注意：len(model) 是链数，不是残基数
    return rec_path, lig_path, nres, ligand_code, auto_note


def mol_from_file(path, remove_hs=False, sanitize=True):
    """读 MOL/SDF —— 不把路径交给 RDKit 的原生文件接口。

    为什么：RDKit 在 Windows 上用 C++ 流打开文件，遇到非 ASCII 路径
    （中文用户名/中文目录）会抛 `OSError: Bad input file ...`。
    这是"换台中文系统的电脑就跑不起来"的典型原因，且报错信息完全不提编码。
    改为 Python 侧读成文本，再用 MolFromMolBlock 解析。
    """
    from rdkit import Chem
    txt = open(path, "r", encoding="utf-8", errors="replace").read()
    return Chem.MolFromMolBlock(txt, removeHs=remove_hs, sanitize=sanitize)


def mol_to_file(mol, path):
    """写 MOL 文件 —— 同样绕开 RDKit 原生写接口。"""
    from rdkit import Chem
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(Chem.MolToMolBlock(mol))


def crystal_reference(ligand_block_pdb, ideal_sdf):
    """用理想 SDF 修复晶体配体的键级，得到带正确化学的"晶体姿态"参照。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    blk = open(ligand_block_pdb, encoding="utf-8").read()
    raw = Chem.MolFromPDBBlock(blk, removeHs=True, proximityBonding=True)
    if raw is None:
        raise RuntimeError("晶体配体无法解析为分子（可能缺键信息）")
    tpl = mol_from_file(ideal_sdf, remove_hs=True, sanitize=True)
    if tpl is None:
        raise RuntimeError("理想 SDF 读取失败：%s" % ideal_sdf)
    if tpl.GetNumAtoms() != raw.GetNumAtoms():
        raise RuntimeError("晶体配体（%d 重原子）与理想构象（%d）原子数不一致"
                           % (raw.GetNumAtoms(), tpl.GetNumAtoms()))
    return AllChem.AssignBondOrdersFromTemplate(tpl, raw)


def prep_ligand_from_smiles(smiles, out_pdbqt, workdir, seed=42, tag="lig"):
    """SMILES → 3D 构象 → PDBQT。返回 (pdbqt, 原子数, 备注)"""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, 0, "SMILES 无法解析"
    m = Chem.AddHs(m)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = seed
    if AllChem.EmbedMolecule(m, ps) != 0:
        ps.useRandomCoords = True
        if AllChem.EmbedMolecule(m, ps) != 0:
            return None, 0, "3D 构象生成失败"
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception:
        pass                    # MMFF 参数缺失不影响后续对接
    sdf = os.path.join(workdir, tag + ".sdf")
    mol_to_file(m, sdf)
    mk = find_tool("mk_prepare_ligand")
    if not mk:
        return None, 0, "缺少 mk_prepare_ligand（pip install meeko）"
    r = run([mk, "-i", sdf, "-o", out_pdbqt])
    if not os.path.exists(out_pdbqt):
        err1 = (r.stderr or r.stdout or "").strip()
        # 某些卤代/含硫分子在 Gasteiger 下出现非有限电荷；Vina 打分不使用电荷，
        # 因此允许坏电荷继续（须在报告中记录，因为其他引擎可能依赖电荷）。
        r2 = run([mk, "-i", sdf, "-o", out_pdbqt, "--bad_charge_ok"])
        if os.path.exists(out_pdbqt):
            return out_pdbqt, m.GetNumAtoms(), "bad_charge_ok（Gasteiger 非有限电荷，已放行）"
        return None, 0, "meeko 未产出 PDBQT：%s" % (err1[-160:] or (r2.stderr or "")[-160:])
    return out_pdbqt, m.GetNumAtoms(), "ok"


def vina_version(vina):
    try:
        r = run([vina, "--version"])
        txt = (r.stdout or "") + (r.stderr or "")
        m = re.search(r"(v?\d+\.\d+\.\d+)", txt)
        return m.group(1) if m else "未知"
    except Exception:
        return "未知"


def dock(vina, receptor_pdbqt, ligand_pdbqt, center, size, out_pdbqt, exhaust=8,
         seed=42, num_modes=5, log_path=None):
    cmd = [vina, "--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt,
           "--center_x", "%.3f" % center[0], "--center_y", "%.3f" % center[1],
           "--center_z", "%.3f" % center[2],
           "--size_x", "%.1f" % size[0], "--size_y", "%.1f" % size[1],
           "--size_z", "%.1f" % size[2],
           "--exhaustiveness", str(exhaust), "--seed", str(seed),
           "--num_modes", str(num_modes), "--out", out_pdbqt]
    r = run(cmd, timeout=1800)
    txt = (r.stdout or "") + (r.stderr or "")
    if log_path:
        C.write_text(log_path, txt)
    if not os.path.exists(out_pdbqt):
        raise RuntimeError("Vina 未产出结果：%s" % txt[-300:])
    aff = re.findall(r"^\s*1\s+(-?\d+\.?\d*)\s+", txt, re.M)
    return float(aff[0]) if aff else None


def pose_rmsd(out_pdbqt, crystal_mol):
    """mode 1（以及所有 mode 的最小值）对晶体姿态的对称校正 RMSD。"""
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    PM = PDBQTMolecule.from_file(out_pdbqt, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(PM)
    if not mols:
        raise RuntimeError("无法把对接姿态转回 RDKit 分子")
    ref = Chem.RemoveHs(crystal_mol)
    pm = mols[0]
    nconf = pm.GetNumConformers()
    vals = []
    for i in range(nconf):
        probe = Chem.Mol(pm, confId=i) if False else Chem.Mol(pm)
        probe.RemoveAllConformers()
        probe.AddConformer(pm.GetConformer(i), assignId=True)
        probe = Chem.RemoveHs(probe)
        try:
            vals.append(rdMolAlign.GetBestRMS(Chem.Mol(probe), Chem.Mol(ref)))
        except Exception:
            continue
    if not vals:
        raise RuntimeError("RMSD 计算失败（原子对应关系无法建立）")
    return vals[0], min(vals), len(vals)


def read_plan(path):
    h, rows = C.read_table(path)
    low = {x.lower().strip(): x for x in h}

    def col(name, extra=(), required=False):
        for c in (name,) + tuple(extra):
            if c.lower() in low:
                return low[c.lower()]
        if required:
            raise SystemExit("plan 表缺少【%s】列。现有列：%s" % (name, "、".join(h)))
        return None

    ct = col("target", ("gene", "term", "protein"), required=True)
    cp = col("pdb", ("pdb_id", "structure"), required=True)
    cl = col("ligand_code", ("ligand", "het", "ccd"), required=True)
    cc = col("chains", ("chain",))
    out = []
    for r in rows:
        chains = (r.get(cc) or "").strip() if cc else ""
        out.append({"target": (r.get(ct) or "").strip().upper(),
                    "pdb": (r.get(cp) or "").strip().lower(),
                    "ligand": (r.get(cl) or "").strip().upper(),
                    "chains": [x.strip() for x in chains.replace(",", " ").split() if x.strip()]})
    return out


def coerce_validation_row(row):
    """把从 CSV 读回的校验行（全为字符串）转成数值类型，供报告格式化与下游使用。
    CSV 无类型：residues/rmsd 等读回是字符串，直接 %d/%.2f 会抛异常（reuse 路径 bug）。"""
    def _num(v, cast=float, default=0):
        try:
            return cast(v)
        except Exception:
            return default
    keys = ("target", "pdb", "ligand_code", "ligand_sel", "residues",
            "redock_affinity", "rmsd_mode1", "rmsd_best", "n_poses",
            "passed", "box", "box_size")
    out = {k: row.get(k, "") for k in keys}
    out["residues"] = _num(out.get("residues"), int, 0)
    out["redock_affinity"] = _num(out.get("redock_affinity"), float, float("nan"))
    out["rmsd_mode1"] = _num(out.get("rmsd_mode1"), float, float("nan"))
    out["rmsd_best"] = _num(out.get("rmsd_best"), float, float("nan"))
    out["n_poses"] = _num(out.get("n_poses"), int, 0)
    return out


def fmt_float(x, digits=2):
    """容错浮点格式化：非数值返回 NA。"""
    try:
        return ("%%.%df" % digits) % float(x)
    except Exception:
        return "NA"


def main(argv=None):
    ap = C.Cli(description="真实对接：结构准备 → Vina → 重对接自校验")
    ap.add_argument("--plan", required=True, help="TSV: target / pdb / ligand_code [/ chains]")
    ap.add_argument("--ligands", default=None, help="成分表（compound,smiles），用于跨对接")
    ap.add_argument("--outdir", default="dock")
    ap.add_argument("--struct-dir", default=None, help="结构缓存目录（默认 outdir/structures）")
    ap.add_argument("--vina", default=None, help="Vina 可执行文件路径")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--num-modes", type=int, default=5)
    ap.add_argument("--box", type=float, default=24.0, help="立方盒子边长 Å（默认 24）")
    ap.add_argument("--rmsd-max", type=float, default=2.0, help="重对接通过阈值 Å")
    ap.add_argument("--validate-only", action="store_true", help="只做重对接校验，不做跨对接")
    ap.add_argument("--max-compounds", type=int, default=0, help="每靶点最多对接几个成分（0=全部）")
    ap.add_argument("--force", action="store_true", help="重对接不通过仍继续跨对接（会写警告）")
    ap.add_argument("--reuse-validation", action="store_true",
                    help="复用已有 docking_validation.csv 中通过校验的受体（避免重复重对接）")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    man = C.Manifest("08_dockprep", a.outdir, vars(a))
    C.ensure_dir(a.outdir)
    struct_dir = a.struct_dir or os.path.join(a.outdir, "structures")
    C.ensure_dir(struct_dir)
    scratch = ascii_scratch()
    work = os.path.join(scratch, "work")
    C.ensure_dir(work)
    poses_dir = os.path.join(a.outdir, "poses")
    rec_dir = os.path.join(a.outdir, "receptors")
    C.ensure_dir(poses_dir)
    C.ensure_dir(rec_dir)

    vina = a.vina or find_tool("vina", [
        os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools", "vina.exe"),
        "/usr/local/bin/vina", "/opt/homebrew/bin/vina",
    ])
    if not vina and not a.validate_only:
        pass
    tools = {"vina": vina, "meeko_lig": find_tool("mk_prepare_ligand"),
             "meeko_rec": find_tool("mk_prepare_receptor")}
    man.count("tools", {k: bool(v) for k, v in tools.items()})
    if not vina:
        print("[阻断] 未找到 AutoDock Vina 可执行文件。")
        print("       获取方式（官方）：https://github.com/ccsb-scripps/AutoDock-Vina/releases")
        print("       Windows 下放 vina.exe 到本技能 tools/ 目录，或加 --vina <路径>。")
        return 3
    if not tools["meeko_lig"] or not tools["meeko_rec"]:
        print("[阻断] 缺少 meeko 命令行（mk_prepare_ligand / mk_prepare_receptor）。")
        print("       安装：python -m pip install meeko")
        return 3
    if not a.quiet:
        print("[对接] Vina %s @ %s" % (vina_version(vina), vina))

    plan = read_plan(a.plan)
    man.count("targets", len(plan))
    val_rows, score_rows, failures = [], [], []

    # ---------------- 1) 重对接自校验 ----------------
    receptors = {}
    reuse = {}
    if a.reuse_validation:
        vp = os.path.join(a.outdir, "docking_validation.csv")
        if os.path.exists(vp):
            _h, vrows = C.read_table(vp)
            for r in vrows:
                if (r.get("passed") or "").upper() != "PASS":
                    continue
                rec = os.path.join(rec_dir, "%s.pdbqt" % r["target"])
                if not os.path.exists(rec):
                    continue
                try:
                    cen = [float(x) for x in r["box"].split(";")]
                    size = [float(x) for x in r["box_size"].split(";")]
                except Exception:
                    continue
                reuse[r["target"]] = {"rec": rec, "center": cen, "size": size,
                                      "pdb": r.get("pdb", ""), "ligand": r.get("ligand_code", ""),
                                      "nres": r.get("residues", "")}
                val_rows.append(coerce_validation_row(r))
            if reuse and not a.quiet:
                print("[对接] 复用 %d 个已通过校验的受体：%s"
                      % (len(reuse), "、".join(sorted(reuse))))

    for entry in plan:
        if entry["target"] in reuse:
            receptors[entry["target"]] = reuse[entry["target"]]
            continue
        tg, pdb, code = entry["target"], entry["pdb"], entry["ligand"]
        pre = os.path.join(work, "%s_%s_%s" % (tg, pdb, code))
        try:
            pdb_path = fetch(RCSB_PDB % pdb.upper(), os.path.join(struct_dir, pdb + ".pdb"))
            # 先拆结构（AUTO 会在此确定真正使用的配体代码），再取该配体的理想构象
            rec_pdb, lig_blk, nres, used_code, auto_note = split_structure(
                pdb_path, code, entry["chains"], pre)
            sdf_path = fetch(RCSB_SDF % used_code, os.path.join(struct_dir, used_code + "_ideal.sdf"))
            cry = crystal_reference(lig_blk, sdf_path)

            # 盒子：晶体配体质心 + 固定边长（可复现）
            conf = cry.GetConformer()
            xs = [conf.GetAtomPosition(i) for i in range(cry.GetNumAtoms())]
            cen = (sum(p.x for p in xs) / len(xs), sum(p.y for p in xs) / len(xs),
                   sum(p.z for p in xs) / len(xs))
            span = (max(p.x for p in xs) - min(p.x for p in xs),
                    max(p.y for p in xs) - min(p.y for p in xs),
                    max(p.z for p in xs) - min(p.z for p in xs))
            size = (max(a.box, span[0] + 10), max(a.box, span[1] + 10), max(a.box, span[2] + 10))
            man.count("box_%s" % tg, {"center": [round(x, 3) for x in cen],
                                      "size": [round(s, 1) for s in size]})

            rec_pdbqt = pre + "_receptor.pdbqt"
            r = run([tools["meeko_rec"], "--read_pdb", rec_pdb, "-o", pre + "_rec", "-p"])
            if not os.path.exists(rec_pdbqt):
                # 第二次：丢弃与模板不匹配的残基（altloc/缺失原子/非标准残基）
                # 必须留痕：删了哪些残基会影响打分，属方法学必须披露的内容。
                r2 = run([tools["meeko_rec"], "--read_pdb", rec_pdb, "-o", pre + "_rec",
                          "-p", "-x"])
                if os.path.exists(rec_pdbqt):
                    man.warn("%s 受体准备需丢弃不匹配残基（-x）才能完成：%s"
                             % (tg, (r.stderr or "").strip().splitlines()[-1][:120]
                                if (r.stderr or "").strip() else "见 work 目录日志"))
            if not os.path.exists(rec_pdbqt):
                # meeko 的 -o 基准名规则可能不同，兜底查找
                cand = [f for f in os.listdir(work) if f.startswith(os.path.basename(pre) + "_rec")
                        and f.endswith(".pdbqt")]
                if cand:
                    rec_pdbqt = os.path.join(work, cand[0])
                else:
                    raise RuntimeError("受体 PDBQT 未生成：%s" % (r.stderr or "")[-200:])

            lig_pdbqt, natoms, note = prep_ligand_from_smiles(
                _mol_to_smiles(cry), pre + "_redock.pdbqt", work, seed=a.seed, tag=tg + "_redock")
            if not lig_pdbqt:
                raise RuntimeError("配体准备失败：%s" % note)

            out_pdbqt = pre + "_redock_out.pdbqt"
            aff = dock(vina, rec_pdbqt, lig_pdbqt, cen, size, out_pdbqt,
                       exhaust=a.exhaustiveness, seed=a.seed, num_modes=a.num_modes,
                       log_path=pre + "_redock.log")
            rms1, rmsmin, npose = pose_rmsd(out_pdbqt, cry)
            passed = rms1 <= a.rmsd_max
            try:
                shutil.copyfile(rec_pdbqt, os.path.join(rec_dir, "%s.pdbqt" % tg))
            except Exception:
                pass
            receptors[tg] = {"rec": rec_pdbqt, "center": cen, "size": size,
                             "pdb": pdb, "ligand": used_code, "nres": nres}
            val_rows.append({"target": tg, "pdb": pdb.upper(), "ligand_code": used_code,
                             "ligand_sel": auto_note,
                             "residues": nres, "redock_affinity": aff,
                             "rmsd_mode1": round(rms1, 2), "rmsd_best": round(rmsmin, 2),
                             "n_poses": npose, "passed": "PASS" if passed else "FAIL",
                             "box": "%.3f;%.3f;%.3f" % cen,
                             "box_size": "%.1f;%.1f;%.1f" % size})
            _ = code
            if not a.quiet:
                print("  [重对接] %-10s %s/%s  RMSD=%.2f Å  %s  (affinity %.2f, %d 姿态)"
                      % (tg, pdb.upper(), code, rms1, "通过" if passed else "未通过", aff, npose))
            if not passed:
                man.warn("%s 重对接 RMSD %.2f Å > %.1f Å：该靶点协议不可信，其打分不得进结论"
                         % (tg, rms1, a.rmsd_max))
        except Exception as e:                                  # noqa: BLE001
            failures.append({"target": tg, "stage": "redock", "error": str(e)[:300]})
            if not a.quiet:
                print("  [重对接] %-10s 失败：%s" % (tg, str(e)[:160]))

    p_val = C.write_table(os.path.join(a.outdir, "docking_validation.csv"),
                          ["target", "pdb", "ligand_code", "ligand_sel", "residues", "redock_affinity",
                           "rmsd_mode1", "rmsd_best", "n_poses", "passed", "box", "box_size"],
                          val_rows)

    # ---------------- 2) 跨对接 ----------------
    if not a.validate_only and a.ligands:
        h, rows = C.read_table(a.ligands)
        cc = C.resolve_col(h, "compound", required=True)
        cs = C.resolve_col(h, "smiles")
        if not cs:
            print("[跳过] 成分表无 smiles 列，无法跨对接。")
        else:
            ligs = [(str(r.get(cc, "")).strip(), str(r.get(cs, "")).strip()) for r in rows]
            ligs = [x for x in ligs if x[0] and x[1]]
            if a.max_compounds:
                ligs = ligs[:a.max_compounds]
            for tg, info in receptors.items():
                ok = any(v["target"] == tg and v["passed"] == "PASS" for v in val_rows)
                if not ok and not a.force:
                    if not a.quiet:
                        print("  [跨对接] %-10s 跳过（重对接未通过）" % tg)
                    continue
                for name, smi in ligs:
                    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", "%s__%s" % (tg, name))[:80]
                    pre = os.path.join(work, safe)
                    try:
                        pdbqt, _n, note = prep_ligand_from_smiles(
                            smi, pre + ".pdbqt", work, seed=a.seed, tag=safe)
                        if not pdbqt:
                            failures.append({"target": tg, "compound": name,
                                             "stage": "ligprep", "error": note})
                            continue
                        out_pdbqt = pre + "_out.pdbqt"
                        aff = dock(vina, info["rec"], pdbqt, info["center"], info["size"],
                                   out_pdbqt, exhaust=a.exhaustiveness, seed=a.seed,
                                   num_modes=a.num_modes, log_path=pre + ".log")
                        pose_dst = os.path.join(poses_dir, safe + "_out.pdbqt")
                        try:
                            shutil.copyfile(out_pdbqt, pose_dst)
                        except Exception:
                            pose_dst = out_pdbqt
                        score_rows.append({
                            "compound": name, "target": tg, "affinity": aff,
                            "redock_rmsd": next((v["rmsd_mode1"] for v in val_rows
                                                 if v["target"] == tg), ""),
                            "pdb": info["pdb"].upper(), "protocol": "vina1.2/meeko0.8/rigid",
                            "pose_file": os.path.relpath(pose_dst, a.outdir),
                        })
                        if not a.quiet:
                            print("    %-28s %-10s %s kcal/mol" % (name[:28], tg, aff))
                    except Exception as e:                      # noqa: BLE001
                        failures.append({"target": tg, "compound": name,
                                         "stage": "dock", "error": str(e)[:200]})

    p_sco = None
    if score_rows:
        p_sco = C.write_table(os.path.join(a.outdir, "docking_scores.csv"),
                              ["compound", "target", "affinity", "redock_rmsd",
                               "pdb", "protocol", "pose_file"], score_rows)

    # ---------------- 3) 报告 ----------------
    n_pass = sum(1 for v in val_rows if v["passed"] == "PASS")
    md = ["# 分子对接准备与自校验报告", "",
          "| 项 | 值 |", "|---|---|",
          "| Vina | %s |" % (vina_version(vina) if vina else "未找到"),
          "| 靶点数 | %d |" % len(plan),
          "| 重对接通过 | %d / %d |" % (n_pass, len(val_rows)),
          "| 跨对接受体 | %d |" % len(receptors),
          "| 对接条目 | %d |" % len(score_rows),
          "| 失败 | %d |" % len(failures),
          "", "## 重对接自校验（RMSD ≤ %.1f Å 才允许其打分进入结论）" % a.rmsd_max, "",
          "| 靶点 | PDB | 共晶配体 | 残基数 | 重对接打分 | RMSD(mode1) | 判定 |",
          "|---|---|---|---|---|---|---|"]
    for v in val_rows:
        md.append("| %s | %s | %s | %s | %s | %s | %s |"
                  % (v["target"], v["pdb"], v["ligand_code"], v["residues"],
                     v["redock_affinity"], fmt_float(v["rmsd_mode1"]), v["passed"]))
    md += ["", "> RMSD 为对称校正后的最优原子对应值（RDKit GetBestRMS）。",
           "> 未通过的靶点**不是**数据缺失，而是「本协议在该靶点上不可靠」，其打分不得用于结论。", ""]
    if failures:
        md += ["## 失败明细", "", "| 靶点 | 成分 | 阶段 | 原因 |", "|---|---|---|---|"]
        for f in failures[:30]:
            md.append("| %s | %s | %s | %s |" % (f.get("target", ""), f.get("compound", "—"),
                                                 f.get("stage", ""), f.get("error", "")))
        md.append("")
    p_md = C.write_text(os.path.join(a.outdir, "dockprep_report.md"), "\n".join(md) + "\n")

    try:
        shutil.rmtree(scratch, ignore_errors=True)
    except Exception:
        pass
    man.count("scratch", os.path.basename(scratch))
    man.output("validation", p_val)
    man.output("report", p_md)
    if p_sco:
        man.output("scores", p_sco)
    man.count("redock_pass", n_pass)
    man.count("scores", len(score_rows))
    man.count("failures", len(failures))
    man.save()

    if not a.quiet:
        print("[对接] 重对接通过 %d/%d，打分 %d 条，失败 %d 条"
              % (n_pass, len(val_rows), len(score_rows), len(failures)))
        print("[对接] 报告 → %s" % p_md)
    return 0


def _mol_to_smiles(mol):
    from rdkit import Chem
    return Chem.MolToSmiles(Chem.RemoveHs(mol))


if __name__ == "__main__":
    sys.exit(main())
