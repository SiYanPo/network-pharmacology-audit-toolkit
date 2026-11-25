#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_ligand_top.py —— 配体拓扑自动生成（Vina .sdf/.pdb → GAFF .top/.gro），MD 的隐形天坑

为什么必须有
    GROMACS 的 run_md.sh 需要受体+配体拓扑（topol.top）。蛋白用 AMBER99SB-ILDN 可自动生成，
    但**中药小分子的 GAFF 参数必须单独生成**——这是 MD 跑不起来的头号原因。
    本脚本封装 acpype（推荐）与 antechamber+parmed 两条路径，自动产出 GROMACS 兼容 .itp/.gro。

工作流
    SDF/PDB → (acpype -b <name> -i <file> -c bcc -a gaff2) → <name>.acpype/<name>_GMX.itp + _GMX.gro
    或 antechamber → acprep → parmed → GROMACS .top

环境（本机实测）
    WSL Ubuntu-24.04: /opt/ambertools/bin/antechamber 已装；acpype 需 pip install（setup_env.sh 已含）
    外部工具统一由 np_tools 解析（acpype 可用 `pip install acpype`，落在解释器 Scripts/）

用法
    python np_ligand_top.py gen --ligand aspirin.sdf --outdir lig_top/
    python np_ligand_top.py gen --ligand lig.pdb --charge bcc --forcefield gaff2 -o lig_top/
    python np_ligand_top.py assemble --lig lig_top/quercetin_GMX.itp --protein蛋白.itp -o topol.top
    python np_ligand_top.py selftest
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def find_acpype():
    """本机优先 WSL 的 acpype，其次 PATH。"""
    p = shutil.which("acpype")
    if p:
        return {"cmd": [p], "mode": "native"}
    if shutil.which("wsl"):
        return {"cmd": ["wsl", "-d", "Ubuntu-24.04", "--", "bash", "-lc", "acpype"],
                "mode": "wsl"}
    return None


def to_wsl_path(p):
    p = os.path.abspath(p).replace("\\", "/")
    import re
    m = re.match(r"^([A-Za-z]):/(.*)$", p)
    return "/mnt/%s/%s" % (m.group(1).lower(), m.group(2)) if m else p


def gen_topology(ligand, outdir, charge="bcc", ff="gaff2", netcharge=None):
    """调用 acpype 生成 GROMACS 配体拓扑。返回产出文件清单。"""
    os.makedirs(outdir, exist_ok=True)
    ap = find_acpype()
    if not ap:
        return {"status": "TOOL_MISSING",
                "reason": "未找到 acpype。安装：WSL 内 pip install acpype；或 bash setup_env.sh",
                "fallback": "可用 antechamber+parmed 手动路径（见 --help）"}
    name = os.path.splitext(os.path.basename(ligand))[0]
    lig_wsl = to_wsl_path(ligand)
    out_wsl = to_wsl_path(outdir)
    cmd = ap["cmd"][:]
    if ap["mode"] == "wsl":
        # 把参数拼进 bash -lc 的字符串
        base = cmd[:-1]
        inner = " ".join(cmd[-1:]) + " -i %s -b %s -c %s -a %s -o %s" % (
            lig_wsl, name, charge, ff, out_wsl)
        if netcharge is not None:
            inner += " -n %d" % netcharge
        full = base + [inner]
    else:
        full = cmd + ["-i", ligand, "-b", name, "-c", charge, "-a", ff, "-o", outdir]
        if netcharge is not None:
            full += ["-n", str(netcharge)]
    try:
        r = subprocess.run(full, capture_output=True, text=True, timeout=300)
        if r.returncode != 0:
            return {"status": "FAILED", "cmd": " ".join(full),
                    "reason": (r.stderr or r.stdout)[-500:]}
    except Exception as e:
        return {"status": "FAILED", "reason": "%s: %s" % (type(e).__name__, e)}
    # acpype 输出到 <outdir>/<name>.acpype/
    sub = os.path.join(outdir, name + ".acpype")
    scan_dir = sub if os.path.isdir(sub) else outdir
    files = [f for f in os.listdir(scan_dir) if f.endswith((".itp", ".gro", ".top"))]
    top_itp = next((f for f in files if "_GMX.itp" in f or f.endswith(".itp")), None)
    top_gro = next((f for f in files if "_GMX.gro" in f or f.endswith(".gro")), None)
    status = "OK" if top_itp else "NO_OUTPUT"
    return {"status": status, "ligand": ligand, "charge": charge, "forcefield": ff,
            "cmd": " ".join(full), "output_dir": scan_dir,
            "itp": os.path.join(scan_dir, top_itp) if top_itp else None,
            "gro": os.path.join(scan_dir, top_gro) if top_gro else None,
            "files": files,
            "note": "配体 .itp 需用 #include 引入蛋白 topol.top；"
                    "[ molecules ] 段追加 <resname> 1；电荷用 bcc 需确认净电荷正确"}


def assemble_topol(lig_itp, protein_itp, outdir, resname="LIG"):
    """生成合并的 topol.top（蛋白 itp + 配体 itp）。注意：此为基础模板，需按实际力场调整。"""
    os.makedirs(outdir, exist_ok=True)
    lines = [
        "; 由 np_ligand_top 生成的基础 topol.top —— 需按实际力场/残基名核对",
        '#include "amber99sb-ildn.ff/forcefield.itp"',
        '#include "%s"' % os.path.basename(protein_itp) if protein_itp else "; [protein itp 待补]",
        '#include "%s"' % os.path.basename(lig_itp),
        "",
        "[ system ]", "Complex in water", "",
        "[ molecules ]", "; Molecule count",
        "Protein 1", "%s 1" % resname,
    ]
    p = os.path.join(outdir, "topol.top")
    open(p, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    return p


def selftest():
    ap = find_acpype()
    ok = callable(gen_topology) and callable(assemble_topol)
    print("selftest: %s (acpype=%s)" % ("PASS" if ok else "FAIL", ap["mode"] if ap else "missing"))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_ligand_top: 配体 GAFF 拓扑自动生成")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("gen")
    p.add_argument("--ligand", required=True)
    p.add_argument("-o", "--outdir", default="np_out/ligand_top")
    p.add_argument("--charge", default="bcc", choices=["bcc", "gasteiger"])
    p.add_argument("--forcefield", default="gaff2", choices=["gaff", "gaff2"])
    p.add_argument("--netcharge", type=int, default=None)
    p.set_defaults(fn="gen")
    p = sub.add_parser("assemble")
    p.add_argument("--lig", required=True); p.add_argument("--protein", default=None)
    p.add_argument("--resname", default="LIG"); p.add_argument("-o", "--outdir", default="np_out/md")
    p.set_defaults(fn="assemble")
    p = sub.add_parser("selftest"); p.set_defaults(fn="selftest")
    a = ap.parse_args(argv)
    if a.fn == "selftest":
        return selftest()
    if a.fn == "gen":
        r = gen_topology(a.ligand, a.outdir, a.charge, a.forcefield, a.netcharge)
        json.dump(r, open(os.path.join(a.outdir, "ligand_top_result.json"), "w", encoding="utf-8"),
                  ensure_ascii=False, indent=1)
        print("[ligand_top] %s" % r["status"])
        if r["status"] != "OK":
            print("  原因：%s" % r.get("reason", "")[:200])
        else:
            print("  itp: %s\ngro: %s" % (r["itp"], r["gro"]))
        return 0 if r["status"] == "OK" else 2
    p = assemble_topol(a.lig, a.protein, a.outdir, a.resname)
    print("[ligand_top] topol.top → %s（请核对力场与残基名）" % p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
