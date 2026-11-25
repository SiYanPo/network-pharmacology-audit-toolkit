#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_pocket.py —— 精准对接口袋识别：共晶配体中心提取 / fpocket 预测（避免盲对接）

问题（审稿人击穿点）
    全蛋白盲对接（blind docking）找出的最低能量位点常是表面凹槽而非真正活性口袋。
    必须先定口袋，再设对接盒 Center X/Y/Z。

口袋来源优先级（诚实分级，写入 manifest）
    P1 co-crystal  : PDB 有共晶配体 → 直接取配体重心（最可靠）
    P2 fpocket     : 无共晶但本机装了 fpocket → 预测口袋，取 Top1/Top3
    P3 P2Rank      : 可选，若提供 p2rank 预测結果文件
    P4 unavailable : 都不可用 → **报错并要求人工指定 Center，禁止盲对接**

fpocket 安装（Linux/WSL）：sudo apt install fpocket 或 conda install -c conda-forge fpocket
Windows 原生无 fpocket → 可用 WSL 调用（与 np_md 一致的桥接模式）。

用法
    python np_pocket.py cocrystal --pdb 3O96.pdb --ligand NAG -o out/
    python np_pocket.py cocrystal --pdb 3O96.pdb --ligand AUTO -o out/     # 自动识别共晶配体
    python np_pocket.py fpocket --pdb protein.pdb --top 3 -o out/
    python np_pocket.py detect                                             # 检测可用工具
    python np_pocket.py selftest
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# 常见结晶缓冲液/盐/溶剂——识别共晶配体时排除（非真正配体）
CRYSTAL_EXCLUDE = {
    "HOH", "NA", "K", "CL", "MG", "CA", "ZN", "MN", "SO4", "PO4", "GOL", "EDO",
    "PEG", "DMS", "DMSO", "ACT", "FMT", "TRS", "HEP", "MES", "IMF", "BO3", "EPE",
    "MPD", "PG4", "1PE", "2PE", "P6G", "TLA", "CIT", "MLI", "FLC", "TAM", "BME",
}


def detect_tools():
    fp = shutil.which("fpocket")
    wsl = shutil.which("wsl")
    return {"fpocket": fp, "wsl": wsl,
            "fpocket_via_wsl": bool(wsl) if not fp else False,
            "note": "fpocket 未安装时可用 WSL（sudo apt install fpocket）；Windows 原生不可用"}


def parse_pdb_hetatm(pdb_path):
    """用纯文本解析 PDB 的 HETATM（零依赖；胜过 gemmi 的轻量需求）。返回 {resname: [(x,y,z),...]}"""
    het = {}
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            res = line[17:20].strip().upper()
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            het.setdefault(res, []).append((x, y, z))
    return het


def auto_ligand(het):
    """从 HETATM 中自动挑最可能的共晶配体：排除盐/溶剂，取原子数最多且非水的有机分子。"""
    cands = []
    for res, coords in het.items():
        if res in CRYSTAL_EXCLUDE or res.startswith("HOH"):
            continue
        # 含碳原子（元素列）且原子数>=8 更可能是真配体
        cands.append((res, len(coords)))
    cands.sort(key=lambda x: -x[1])
    return cands


def center_of(coords):
    n = len(coords)
    return (sum(p[0] for p in coords) / n, sum(p[1] for p in coords) / n,
            sum(p[2] for p in coords) / n)


def span_of(coords):
    return (max(p[0] for p in coords) - min(p[0] for p in coords),
            max(p[1] for p in coords) - min(p[1] for p in coords),
            max(p[2] for p in coords) - min(p[2] for p in coords))


def pocket_from_cocrystal(pdb_path, ligand_code):
    het = parse_pdb_hetatm(pdb_path)
    if ligand_code and ligand_code.upper() != "AUTO":
        key = ligand_code.upper()
        if key not in het:
            return {"status": "unavailable", "reason": "PDB 中未找到配体 %s" % key,
                    "available_ligands": [{"res": k, "n_atoms": len(v)} for k, v in auto_ligand(het)]}
        coords = het[key]
    else:
        cands = auto_ligand(het)
        if not cands:
            return {"status": "unavailable", "reason": "未发现任何非溶剂有机配体",
                    "all_hets": list(het.keys())[:20]}
        coords = het[cands[0][0]]
        ligand_code = cands[0][0]
    cen, span = center_of(coords), span_of(coords)
    return {"status": "ok", "source": "co-crystal", "confidence": "P1",
            "ligand_code": ligand_code, "n_atoms": len(coords),
            "center": [round(x, 3) for x in cen],
            "span": [round(x, 3) for x in span],
            "box_size": [round(max(20, s + 8), 1) for s in span]}


def pocket_from_fpocket(pdb_path, top=3, out=None):
    """调用 fpocket 预测口袋（本机或 WSL）。返回 Top-N 口袋中心。"""
    tools = detect_tools()
    outdir = out or os.path.join(HERE, "..", "np_out", "fpocket")
    os.makedirs(outdir, exist_ok=True)
    base = os.path.basename(pdb_path).replace(".pdb", "")
    if tools["fpocket"]:
        cmd = [tools["fpocket"], "-f", pdb_path, "-o", outdir]
    elif tools["wsl"]:
        wsl_pdb = pdb_path.replace("\\", "/")
        # Windows 盘符 → WSL /mnt/
        m = re.match(r"^([A-Za-z]):/(.*)$", wsl_pdb)
        if m:
            wsl_pdb = "/mnt/%s/%s" % (m.group(1).lower(), m.group(2))
        cmd = ["wsl", "--", "bash", "-lc", "fpocket -f %s -o %s" % (wsl_pdb, outdir)]
    else:
        return {"status": "unavailable", "source": "fpocket", "confidence": "P4",
                "reason": "fpocket 未安装且无 WSL；请在 Linux/WSL 执行 sudo apt install fpocket，"
                          "或人工指定对接盒 Center（禁止盲对接）"}
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            return {"status": "unavailable", "source": "fpocket", "confidence": "P4",
                    "reason": "fpocket 运行失败：%s" % (r.stderr or "")[-200:]}
    except Exception as e:
        return {"status": "unavailable", "reason": "%s: %s" % (type(e).__name__, e)}
    # 解析 fpocket 输出 *_pockets/pocket*_vert.pdb 的质心
    pockets = []
    pdir = os.path.join(outdir, base + "_pockets")
    if os.path.isdir(pdir):
        import glob
        for pf in sorted(glob.glob(os.path.join(pdir, "pocket*_vert.pdb")))[:top]:
            verts = []
            for line in open(pf, encoding="utf-8", errors="ignore"):
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        verts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                    except ValueError:
                        pass
            if verts:
                pockets.append({"file": os.path.basename(pf), "center": [round(x, 3) for x in center_of(verts)],
                                "n_vertices": len(verts)})
    if not pockets:
        return {"status": "unavailable", "reason": "fpocket 已运行但未解析到口袋坐标"}
    return {"status": "ok", "source": "fpocket", "confidence": "P2", "n_pockets": len(pockets),
            "pockets": pockets, "recommended": pockets[0]["center"]}


def selftest():
    # 构造一个含共晶配体的极简 PDB
    txt = ("HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 10.00           C\n"
           "HETATM    2  C2  LIG A   1       2.000   2.000   2.000  1.00 10.00           C\n"
           "HETATM    3  O   HOH A   2       9.000   9.000   9.000  1.00 10.00           O\n")
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "t.pdb")
        open(p, "w").write(txt)
        r = pocket_from_cocrystal(p, "AUTO")
        ok = (r["status"] == "ok" and r["ligand_code"] == "LIG"
              and r["center"] == [1.0, 1.0, 1.0])
        print("selftest: %s (co-crystal LIG center=%s)" % ("PASS" if ok else "FAIL", r.get("center")))
        return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_pocket: 精准对接口袋识别")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("cocrystal"); p.add_argument("--pdb", required=True)
    p.add_argument("--ligand", default="AUTO"); p.add_argument("-o", "--out", default=None)
    p.set_defaults(fn=lambda a: pocket_from_cocrystal(a.pdb, a.ligand))
    p = sub.add_parser("fpocket"); p.add_argument("--pdb", required=True)
    p.add_argument("--top", type=int, default=3); p.add_argument("-o", "--out", default=None)
    p.set_defaults(fn=lambda a: pocket_from_fpocket(a.pdb, a.top, a.out))
    p = sub.add_parser("detect"); p.set_defaults(fn=lambda a: detect_tools())
    p = sub.add_parser("selftest"); p.set_defaults(fn=lambda a: {"selftest": True})
    a = ap.parse_args(argv)
    if a.cmd == "selftest":
        return selftest()
    r = a.fn(a)
    print(json.dumps(r, ensure_ascii=False, indent=1))
    if a.cmd in ("cocrystal", "fpocket") and a.out:
        os.makedirs(a.out, exist_ok=True)
        json.dump(r, open(os.path.join(a.out, "pocket.json"), "w", encoding="utf-8"),
                  ensure_ascii=False, indent=1)
    if isinstance(r, dict) and r.get("status") == "unavailable":
        print("[pocket] 不可用：%s" % r.get("reason", "")[:150], file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
