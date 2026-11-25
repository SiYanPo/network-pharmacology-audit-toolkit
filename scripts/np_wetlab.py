#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_wetlab.py —— 湿实验升维三件套：引物特异性 / MD 桥接 / 亚细胞定位

三个子命令：
  specificity  本地 BLAST+/isPcr 验证引物特异性（工具缺失时生成可执行脚本+诚实降级）
  md-bridge    从 Vina 复合物生成 GROMACS 最小化+NVT/NPT+10ns 生产相 bash 模板
  localization UniProt/HPA 亚细胞定位 → 分层渲染 + 膜穿透差成分配核内靶点标红

诚实边界（实测）
    * 本机 blastn/isPcr/bowtie2/gmx/vina **均未安装** → 本模块生成**可直接运行的脚本模板**，
      检测不到工具时不伪造比对/模拟结果，明确标注 TOOL_MISSING。
    * UniProt 亚细胞定位 API 可用（comments→SUBCELLULAR LOCATION），HPA 需另行下载。

用法
    python np_wetlab.py specificity --primers mouse_primers.csv --genome GRCm39.fa -o out/
    python np_wetlab.py md-bridge --receptor rec.pdbqt --ligand lig.pdbqt -o out/
    python np_wetlab.py localization --targets core.csv --compounds compounds.csv -o out/
    python np_wetlab.py selftest
"""
import argparse
import csv
import json
import os
import re
import shutil
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------- 1. 引物本地特异性 ----------
def cmd_specificity(a):
    os.makedirs(a.outdir, exist_ok=True)
    rows = []
    with open(a.primers, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    with open(a.primers, encoding="utf-8-sig") as f:
        hdrs = (csv.DictReader(f).fieldnames or [])
    fc = next((h for h in hdrs if "forward" in h.lower()), "forward")
    rc = next((h for h in hdrs if "reverse" in h.lower()), "reverse")
    pairs = [(r.get("gene", ""), r.get(fc, ""), r.get(rc, "")) for r in rows if r.get(fc) and r.get(rc)]

    tools = {"blastn": shutil.which("blastn"), "isPcr": shutil.which("isPcr"),
             "bowtie2": shutil.which("bowtie2")}
    # 生成 in-silico PCR 批处理脚本（isPcr：直接从基因组扩增，最贴合 PCR）
    script = ["#!/usr/bin/env bash", "# 由 np_wetlab 生成：引物本地特异性验证（GRCm39 基因组）",
              "set -euo pipefail", 'GENOME="%s"' % (a.genome or "GRCm39.fa"),
              'OUT="%s"' % a.outdir, "mkdir -p \"$OUT\""]
    for g, f, r in pairs:
        script.append('isPcr "$GENOME" <(echo ">%s_F"\\n"%s") <(echo ">%s_R"\\n"%s") '
                      '"$OUT/%s.amplicons.fa" -minPerfect=15 -minGood=25 2>>"$OUT/ispcr.log" || true'
                      % (g, f, g, r, g))
    script.append("# 判定：每个引物对扩增产物应恰好 1 条且长度符合预期；多条=非特异，需用 Primer-BLAST 重设计")
    sp = os.path.join(a.outdir, "check_primer_specificity.sh")
    open(sp, "w", encoding="utf-8").write("\n".join(script) + "\n")

    if not any(tools.values()):
        status = "TOOL_MISSING"
        note = ("未检测到 blastn/isPcr/bowtie2。已生成可执行脚本 %s；"
                "安装：conda install -c bioconda blast ispcr bowtie2；下载 GRCm39："
                "https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/" % sp)
    else:
        status = "TOOL_AVAILABLE_SCRIPT_GENERATED"
        note = "检测到比对工具，运行 %s 完成本地特异性验证" % sp
    with open(os.path.join(a.outdir, "specificity_plan.json"), "w", encoding="utf-8") as f:
        json.dump({"status": status, "tools": tools, "n_primers": len(pairs),
                   "script": os.path.basename(sp), "note": note}, f, ensure_ascii=False, indent=1)
    print("[specificity] %d 对引物，工具=%s → %s" % (len(pairs), status, sp))
    return 0


# ---------- 2. MD 桥接 ----------
MD_SCRIPT = r"""#!/usr/bin/env bash
# 由 np_wetlab md-bridge 生成：GROMACS 短程 MD（最小化 → NVT → NPT → 10ns 生产相）
# 前置：复合物已由 Vina 对接得到；需准备受体/配体拓扑（AMBER99SB-ILDN + GAFF2 或 CGenFF）
# 注意：配体拓扑需 acpype/antechamber 生成，脚本不臆造参数
set -euo pipefail
COMPLEX="{complex}"
OUT="{outdir}"
mkdir -p "$OUT"

echo "== 1. 能量最小化 (steep) =="
gmx grompp -f em.mdp -c "$COMPLEX" -p topol.top -o "$OUT/em.tpr" -maxwarn 2
gmx mdrun -deffnm "$OUT/em"

echo "== 2. NVT 平衡 (100 ps) =="
gmx grompp -f nvt.mdp -c "$OUT/em.gro" -r "$OUT/em.gro" -p topol.top -o "$OUT/nvt.tpr"
gmx mdrun -deffnm "$OUT/nvt"

echo "== 3. NPT 平衡 (100 ps) =="
gmx grompp -f npt.mdp -c "$OUT/nvt.gro" -r "$OUT/nvt.gro" -p topol.top -o "$OUT/npt.tpr"
gmx mdrun -deffnm "$OUT/npt"

echo "== 4. 生产相 MD (10 ns) =="
gmx grompp -f md_10ns.mdp -c "$OUT/npt.gro" -p topol.top -o "$OUT/md.tpr"
gmx mdrun -deffnm "$OUT/md" -ntomp 4

echo "== 5. 分析：RMSD / RMSF / MM-PBSA =="
echo "0\n" | gmx rms -s "$OUT/md.tpr" -f "$OUT/md.xtc" -o "$OUT/rmsd.xvg"
echo "0\n" | gmx rmsf -s "$OUT/md.tpr" -f "$OUT/md.xtc" -o "$OUT/rmsf.xvg"
# gmx_MMPBSA 需单独安装（conda install -c conda-forge gmx_mmpbsa）
# gmx_MMPBSA -O -i mmpbsa.in -cs md.tpr -ci index.ndx -cg 1 -ct md.xtc -o "$OUT/mmpbsa.dat"
echo "完成。RMSD 平稳(<3Å 波动)=复合物稳定；RMSF 峰=柔性区；MM-PBSA 给结合自由能。"
"""


def cmd_mdbridge(a):
    os.makedirs(a.outdir, exist_ok=True)
    gmx = shutil.which("gmx")
    p = os.path.join(a.outdir, "run_md.sh")
    open(p, "w", encoding="utf-8").write(
        MD_SCRIPT.format(complex=a.receptor, outdir=a.outdir))
    os.chmod(p, 0o755)
    mdps = {
        "em.mdp": "integrator = steep\nemtol = 1000\nemstep = 0.01\nnsteps = 50000\ncutoff-scheme = Verlet\n",
        "nvt.mdp": "integrator = md\ndt = 0.002\nnsteps = 50000\ntcoupl = v-rescale\ntc-grps = System\nref-t = 300\nconstraints = h-bonds\n",
        "npt.mdp": "integrator = md\ndt = 0.002\nnsteps = 50000\npcoupl = Parrinello-Rahman\nref-p = 1.0\ncontinuation = yes\n",
        "md_10ns.mdp": "integrator = md\ndt = 0.002\nnsteps = 5000000  ; 10 ns\n"
                       "tcoupl = v-rescale\npcoupl = Parrinello-Rahman\nref-t = 300\nref-p = 1.0\n"
                       "xtc-grps = System\nnstxout-compressed = 5000\n",
    }
    for fn, content in mdps.items():
        open(os.path.join(a.outdir, fn), "w", encoding="utf-8").write(content)
    json.dump({"gmx_available": bool(gmx), "script": os.path.basename(p),
               "status": "READY_TO_RUN" if gmx else "TOOL_MISSING",
               "note": "gmx 未安装：脚本与 mdp 已生成，安装 gromacs 后可直接运行；"
                       "配体拓扑需 acpype/antechamber 生成，不得臆造"},
              open(os.path.join(a.outdir, "md_bridge_meta.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[md-bridge] 脚本+4 个 mdp → %s（gmx=%s）" % (a.outdir, "可用" if gmx else "未安装，已生成模板"))
    return 0


# ---------- 3. 亚细胞定位 ----------
UNIPROT = "https://rest.uniprot.org/uniprotkb/%s.json"


def get_localization(accession):
    try:
        req = urllib.request.Request(UNIPROT % accession, headers={"User-Agent": "np-pro-local/1.0"})
        d = json.loads(urllib.request.urlopen(req, timeout=20).read().decode("utf-8", "ignore"))
        locs = []
        for c in (d.get("comments") or []):
            if c.get("commentType") == "SUBCELLULAR LOCATION":
                for x in (c.get("subcellularLocations") or []):
                    v = (x.get("location") or {}).get("value")
                    if v:
                        locs.append(v)
        return sorted(set(locs))
    except Exception:
        return []


def classify_location(locs):
    s = " ".join(locs).lower()
    layers = []
    if any(k in s for k in ("nucleus", "nuclear", "chromatin", "transcription")):
        layers.append("细胞核")
    if any(k in s for k in ("membrane", "plasma membrane", "receptor")):
        layers.append("细胞膜")
    if "mitochondr" in s:
        layers.append("线粒体")
    if any(k in s for k in ("cytoplasm", "cytosol")):
        layers.append("细胞质")
    if any(k in s for k in ("endoplasmic", "golgi", "vesicle")):
        layers.append("细胞器")
    return layers or ["未注释"]


def cmd_localization(a):
    os.makedirs(a.outdir, exist_ok=True)
    targets = []
    with open(a.targets, encoding="utf-8-sig") as f:
        first = f.readline()
        if "," in first or "\t" in first:
            delim = "\t" if "\t" in first else ","
            f.seek(0)
            rdr = csv.DictReader(f, delimiter=delim)
            # v3.0.2: 兼容 np_ppi hubs.csv 的 node 列（跨模块列名契约）
            gc = next((h for h in (rdr.fieldnames or []) if h.lower() in
                       ("gene", "target", "node", "symbol", "hgnc", "靶点")), None)
            ac = next((h for h in (rdr.fieldnames or []) if "uniprot" in h.lower() or "accession" in h.lower()), None)
            targets = [(r.get(gc, ""), r.get(ac, "")) for r in rdr if r.get(gc)]
        else:
            f.seek(0)
            targets = [(l.strip(), "") for l in f if l.strip()]
    rows = []
    for gene, acc in targets[:50]:
        locs = get_localization(acc) if acc else []
        layers = classify_location(locs)
        rows.append({"gene": gene, "uniprot": acc,
                     "locations": ";".join(locs[:4]), "layers": ";".join(layers),
                     "nuclear": "是" if "细胞核" in layers else "否"})
    # 膜穿透性弱的成分（大分子/多糖/高极性）配核内靶点 → 标红警告
    warnings = []
    low_bbb = set()
    if a.compounds and os.path.exists(a.compounds):
        with open(a.compounds, encoding="utf-8-sig") as f:
            rdr = csv.DictReader(f)
            bc = next((h for h in (rdr.fieldnames or []) if "bbb" in h.lower()), None)
            for r in rdr:
                nm = (r.get("name") or r.get("compound") or "").strip()
                if nm and (str(r.get(bc, "")).lower() in ("no", "0", "low", "低") if bc else False):
                    low_bbb.add(nm)
    for r in rows:
        if r["nuclear"] == "是" and low_bbb:
            warnings.append("核内靶点 %s 搭配 BBB 阴性成分（%s）→ "
                            "Physically Unreachable：成分无法透膜到达核内靶点，机制链在空间上断裂，"
                            "不得作为核心机制解释" % (r["gene"], ",".join(sorted(low_bbb)[:3])))
        elif r["nuclear"] == "是":
            warnings.append("核内靶点 %s：需核查成分核膜通透性（未标注 BBB 数据时默认提示）" % r["gene"])
    with open(os.path.join(a.outdir, "localization_layers.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "uniprot", "locations", "layers", "nuclear"])
        w.writeheader(); w.writerows(rows)
    # 分层渲染指令（供 np_tcmviz/Cytoscape 使用）
    json.dump({"note": "按 layers 字段分层渲染：细胞膜/细胞质/细胞核/线粒体；"
                       "nuclear=是 且成分 BBB=low 时标红警告",
               "n_with_location": sum(1 for r in rows if r["locations"]),
               "red_flags": warnings},
              open(os.path.join(a.outdir, "localization_warnings.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[localization] %d 靶点，%d 个有定位，%d 条核内错配警告 → %s" %
          (len(rows), sum(1 for r in rows if r["locations"]), len(warnings), a.outdir))
    return 0


def selftest():
    ok = ("gmx mdrun" in MD_SCRIPT and "10 ns" in MD_SCRIPT
          and classify_location(["Nucleus", "Cytoplasm"]) == ["细胞核", "细胞质"])
    print("selftest: %s (md template %d chars, layer=%s)" %
          ("PASS" if ok else "FAIL", len(MD_SCRIPT), classify_location(["Nucleus"])))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_wetlab: 引物特异性/MD桥接/亚细胞定位")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("specificity"); p.add_argument("--primers", required=True)
    p.add_argument("--genome", default=None); p.add_argument("-o", "--outdir", default="np_out/wetlab")
    p.set_defaults(fn=cmd_specificity)
    p = sub.add_parser("md-bridge"); p.add_argument("--receptor", required=True)
    p.add_argument("--ligand", default="lig.pdbqt"); p.add_argument("-o", "--outdir", default="np_out/md")
    p.set_defaults(fn=cmd_mdbridge)
    p = sub.add_parser("localization"); p.add_argument("--targets", required=True)
    p.add_argument("--compounds", default=None); p.add_argument("-o", "--outdir", default="np_out/localization")
    p.set_defaults(fn=cmd_localization)
    p = sub.add_parser("selftest"); p.set_defaults(fn=lambda a: selftest())
    a = ap.parse_args(argv)
    if a.cmd == "selftest":
        return selftest()
    return a.fn(a) or 0


if __name__ == "__main__":
    sys.exit(main())
