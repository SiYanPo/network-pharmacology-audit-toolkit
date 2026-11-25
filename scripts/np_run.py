#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_run.py — 全流程编排器

设计要点
  * 阶段之间**只通过文件通信**（`<outdir>/` 下的固定文件名），所以任何一步都能单独重跑：
    改阈值重跑 04 不必重抓 TCMSP，改 GMT 重跑 06 不必重算 PPI。
  * 每个阶段自己写 `run_manifest.json`，最后合并为一份完整溯源清单。
  * `--dry-run` 只做输入校验与阶段清单打印，不产生任何计算，适合先确认参数。
  * `--demo` 用内置演示数据（**构造数据**）离线跑完 01–11，用于验证环境与理解产物。

用法
  python np_run.py --demo
  python np_run.py --targets targets.csv --compounds compounds.csv \
                   --disease disease.tsv --edges string_edges.tsv --gmt pathways.gmt \
                   --label 黄芪-当归-糖尿病肾病 -o out/
  python np_run.py --list
  python np_run.py ... --from 04 --stages 04,06
"""
import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import npcore as C  # noqa: E402

STAGES = [
    ("01", "standardize", "靶点标准化（蛋白名→基因符号）"),
    ("02", "compounds", "成分筛选（OB/DL/途径/泛成分）"),
    ("03", "intersect", "交集 + 成分级富集"),
    ("04", "ppi", "PPI 拓扑 + hub 稳健性"),
    ("05", "null", "四套零模型阴性对照"),
    ("06", "enrich", "富集 + 冗余折叠 + 三层分级"),
    ("07", "dose", "剂量可达性核算"),
    ("08", "dock", "对接门禁（可先跑 08a 真实对接）"),
    ("09", "optimize", "反向网络药理学 / 组合优化"),
    ("10", "audit", "项目红线审计"),
    ("11", "report", "自包含 HTML 报告（纳入审计结果）"),
]


def main(argv=None):
    ap = argparse.ArgumentParser(description="网络药理学全流程编排")
    ap.add_argument("--list", action="store_true", help="列出阶段")
    ap.add_argument("--demo", action="store_true", help="用内置演示数据跑全流程（离线）")
    ap.add_argument("--targets", help="原始靶点表（含 target_name/gene 列）")
    ap.add_argument("--compounds", help="成分表（OB/DL…）")
    ap.add_argument("--disease", help="疾病靶点表")
    ap.add_argument("--edges", help="PPI 边表")
    ap.add_argument("--gmt", help="本地 GMT 通路文件")
    ap.add_argument("--docking", help="对接结果表")
    ap.add_argument("--label", default="网络药理学分析", help="项目标签（用于报告标题）")
    ap.add_argument("--route", default="unknown", choices=["oral", "injection", "topical", "unknown"])
    ap.add_argument("--ob", type=float, default=30.0, help="OB 阈值%%，透传给 02 阶段")
    ap.add_argument("--dl", type=float, default=0.18, help="DL 阈值，透传给 02 阶段")
    ap.add_argument("--hl", type=float, default=None, help="HL 阈值 h，透传给 02 阶段")
    ap.add_argument("--logic", default="and", choices=["and", "or"], help="阈值组合逻辑")
    ap.add_argument("--exceptions", default=None, help="例外纳入名单（成分名逗号分隔或 CSV）")
    ap.add_argument("--background", default=None, help="富集背景基因列表，透传给 06 阶段")
    ap.add_argument("--calibrate", type=int, default=0, help="富集经验校准的置换次数（0=关闭）")
    # --- v1.2: 分子层面门控 ---
    ap.add_argument("--gate", default="auto", choices=["auto", "obdl", "druglike", "evidence", "none"],
                    help="02 阶段门控来源（auto=有 OB/DL 用 OB/DL，否则用类药性）")
    ap.add_argument("--qed-min", type=float, default=0.30, help="类药性门 QED 下限")
    ap.add_argument("--max-veber", type=int, default=1, help="类药性门允许的 Veber 违例数")
    ap.add_argument("--obdl-snapshot", default=None, help="OB/DL 快照表（TCMSP 导出）")
    # --- v1.2: 真实对接 ---
    ap.add_argument("--dockplan", default=None, help="对接计划 TSV（target/pdb/ligand_code），启用 08a 真实对接")
    ap.add_argument("--vina", default=None, help="Vina 可执行文件路径")
    ap.add_argument("--rmsd-max", type=float, default=2.0, help="重对接通过阈值 Å")
    ap.add_argument("--validate-only", action="store_true", help="只做重对接校验")
    ap.add_argument("--max-dock-compounds", type=int, default=0, help="每靶点最多对接几个成分")
    ap.add_argument("--unit", default="compound", choices=["compound", "herb"],
                    help="09 阶段优化单位")
    ap.add_argument("--largest-only", action="store_true", help="04 阶段只用最大连通分量")
    ap.add_argument("--min-score", type=float, default=0.4, help="PPI 边阈值")
    ap.add_argument("--top", type=int, default=15, help="hub 取前 k")
    ap.add_argument("--bootstrap", type=int, default=200)
    ap.add_argument("--iterations", type=int, default=500, help="零模型迭代次数")
    ap.add_argument("--daily", type=float, default=None, help="饮片日剂量 g（启用 07）")
    ap.add_argument("--content", type=float, default=None, help="成分含量%%（启用 07）")
    ap.add_argument("--mw", type=float, default=None, help="分子量（启用 07）")
    ap.add_argument("--invitro", type=float, default=None, help="体外浓度 µM（启用 07）")
    ap.add_argument("--compound-name", default="目标成分")
    ap.add_argument("--herb-name", default="受试药")
    ap.add_argument("-o", "--outdir", default="np_out")
    ap.add_argument("--stages", default=None, help="只跑指定阶段，如 03,04")
    ap.add_argument("--from", dest="from_stage", default=None, help="从某阶段跑到尾")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    if a.list:
        print("阶段清单：")
        for sid, mod, desc in STAGES:
            print("  %s  %-12s %s" % (sid, mod, desc))
        return 0

    outdir = os.path.abspath(a.outdir)
    demo = os.path.join(HERE, "..", "examples")

    if a.demo:
        a.targets = os.path.join(demo, "compound_targets.csv")
        a.compounds = os.path.join(demo, "compounds.csv")
        a.disease = os.path.join(demo, "disease_targets.tsv")
        a.edges = os.path.join(demo, "string_edges.tsv")
        a.gmt = os.path.join(demo, "pathways.gmt")
        a.docking = os.path.join(demo, "docking_scores.csv")
        a.label = "演示：黄芪-当归-丹参 vs 糖尿病肾病（构造数据）"
        a.route = "oral"
        a.daily, a.content, a.mw, a.invitro = 9.0, 0.042, 302.24, 50.0
        a.compound_name, a.herb_name = "Quercetin", "黄芪"
        a.bootstrap = 120
        a.iterations = 300

    # ---- 输入校验 ----
    problems = []
    if not a.targets:
        problems.append("缺 --targets（可用 --demo）")
    for label, p in [("targets", a.targets), ("compounds", a.compounds), ("disease", a.disease),
                     ("edges", a.edges), ("gmt", a.gmt), ("docking", a.docking),
                     ("dockplan", a.dockplan), ("background", a.background),
                     ("obdl_snapshot", a.obdl_snapshot), ("exceptions", a.exceptions)]:
        if not p:
            continue
        # --gmt / --exceptions 支持逗号分隔多文件：逐个校验，不能把整串当路径
        parts = [x.strip() for x in str(p).split(",") if x.strip()] if label in ("gmt", "exceptions") else [p]
        for part in parts:
            if part and not os.path.exists(part):
                problems.append("文件不存在：%s = %s" % (label, part))
    if problems:
        print("[配置错误]")
        for p in problems:
            print("  · " + p)
        return 2

    chosen = {s[0] for s in STAGES}
    if a.stages:
        chosen = {x.strip().zfill(2) for x in a.stages.split(",") if x.strip()}
    if a.from_stage:
        fs = a.from_stage.strip().zfill(2)
        chosen = {sid for sid, _m, _d in STAGES if sid >= fs}

    print("=" * 64)
    print(" 网络药理学研究平台 v%s — %s" % (C.VERSION, a.label))
    print(" 输出目录：%s" % outdir)
    print(" 阶段：" + ", ".join(sid for sid, _m, _d in STAGES if sid in chosen))
    print("=" * 64)
    if a.dry_run:
        print("\n[dry-run] 输入检查通过。将执行：")
        for sid, mod, desc in STAGES:
            if sid in chosen:
                print("  %s %s" % (sid, desc))
        print("\n[dry-run] 未执行任何计算。")
        return 0

    C.ensure_dir(outdir)

    def run(sid, mod_name, args):
        import importlib
        mod = importlib.import_module("np_" + mod_name)
        t0 = time.time()
        try:
            rc = mod.main(args)
        except SystemExit as e:
            rc = e.code if isinstance(e.code, int) else 0
        except Exception as e:
            print("  [%s] 异常：%s: %s" % (sid, type(e).__name__, e))
            return 1
        print("  [%s] 用时 %.1fs%s" % (sid, time.time() - t0, "" if rc in (0, None) else "（返回码 %s）" % rc))
        return rc or 0

    std_csv = os.path.join(outdir, "targets_standardized.csv")
    inter_genes = os.path.join(outdir, "intersection_genes.txt")
    pe = os.path.join(outdir, "ppi_edges_clean.tsv")
    q = ["--quiet"] if a.quiet else []

    if "01" in chosen:
        print("\n[01] 靶点标准化")
        run("01", "standardize", ["-i", a.targets, "-o", outdir] + q)

    if "02" in chosen and a.compounds:
        print("\n[02] 成分筛选")
        c_args = ["-i", a.compounds, "-o", outdir, "--route", a.route,
                  "--ob", str(a.ob), "--dl", str(a.dl), "--logic", a.logic,
                  "--force"] + q
        if a.hl is not None:
            c_args += ["--hl", str(a.hl)]
        if a.exceptions:
            c_args += ["--exceptions", a.exceptions]
        c_args += ["--gate", a.gate, "--qed-min", str(a.qed_min), "--max-veber", str(a.max_veber)]
        if a.obdl_snapshot:
            c_args += ["--obdl-snapshot", a.obdl_snapshot]
        run("02", "compounds", c_args)

    if "03" in chosen:
        print("\n[03] 交集 + 成分级富集")
        src = std_csv if os.path.exists(std_csv) else a.targets
        args = ["--disease", a.disease, "-o", outdir, "--label", a.label] + q
        if os.path.exists(std_csv):
            args = ["--drug", std_csv, "--gene-col", "symbol"] + args
        else:
            args = ["--drug", src] + args
        run("03", "intersect", args)

    if "04" in chosen and a.edges:
        print("\n[04] PPI 拓扑与 hub")
        run("04", "ppi", ["-e", a.edges, "-o", outdir, "--min-score", str(a.min_score),
                          "--top", str(a.top), "--bootstrap", str(a.bootstrap),
                          "--label", a.label] + (["--largest-only"] if a.largest_only else []) + q)

    if "05" in chosen and a.edges and os.path.exists(inter_genes):
        print("\n[05] 零模型")
        run("05", "null", ["-g", inter_genes, "-e", a.edges, "-d", a.disease,
                           "-o", outdir, "-B", str(a.iterations)] + q)
    elif "05" in chosen:
        print("\n[05] 跳过：缺少交集基因表（先跑 03）")

    if "06" in chosen and os.path.exists(inter_genes):
        print("\n[06] 富集")
        args = ["-g", inter_genes, "-o", outdir] + q
        if a.gmt:
            args += ["-m", a.gmt]
        if a.background:
            args += ["--background", a.background]
        if a.calibrate:
            args += ["--calibrate", str(a.calibrate)]
        run("06", "enrich", args)
    elif "06" in chosen:
        print("\n[06] 跳过：缺少交集基因表")

    if "07" in chosen and all(v is not None for v in (a.daily, a.content, a.mw)):
        print("\n[07] 剂量可达性")
        args = ["--daily", str(a.daily), "--content", str(a.content), "--mw", str(a.mw),
                "--herb", a.herb_name, "--compound", a.compound_name, "-o", outdir] + q
        if a.invitro:
            args += ["--invitro", str(a.invitro)]
        run("07", "dose", args)
    elif "07" in chosen:
        print("\n[07] 跳过：需要 --daily/--content/--mw")

    dock_table = a.docking
    if "08" in chosen and a.dockplan:
        print("\n[08a] 真实对接（准备 → Vina → 重对接自校验）")
        d_args = ["--plan", a.dockplan, "--outdir", os.path.join(outdir, "dock"),
                  "--rmsd-max", str(a.rmsd_max)]
        if a.vina:
            d_args += ["--vina", a.vina]
        if a.validate_only:
            d_args += ["--validate-only"]
        if a.max_dock_compounds:
            d_args += ["--max-compounds", str(a.max_dock_compounds)]
        if a.compounds:
            d_args += ["--ligands", a.compounds]
        rc8a = run("08a", "dockprep", d_args + q)
        cand = os.path.join(outdir, "dock", "docking_scores.csv")
        if os.path.exists(cand):
            dock_table = cand
            print("  [08a] 打分表 → %s" % cand)

    if "08" in chosen and dock_table:
        print("\n[08] 对接门禁")
        run("08", "dock", ["-i", dock_table, "-o", outdir] + q)
    elif "08" in chosen:
        print("\n[08] 跳过：无对接结果表（--docking 或 --dockplan）")

    if "09" in chosen and os.path.exists(inter_genes):
        print("\n[09] 反向网络药理学")
        if not os.path.exists(std_csv):
            print("  [09] 跳过：缺少标准化靶点表（反向优化不接受未标准化的蛋白名）")
        else:
            o_args = ["--drug", std_csv, "--gene-col", "symbol", "--disease", a.disease,
                      "-o", outdir, "--max-k", "4", "--unit", a.unit] + q
            kept_names = os.path.join(outdir, "compounds_kept_names.txt")
            if os.path.exists(kept_names):
                o_args += ["--kept", kept_names]
            run("09", "optimize", o_args)

    if "10" in chosen or "11" in chosen:
        # 账本必须在审计之前生成（否则「存在账本」会永远 FAIL —— 顺序依赖）
        try:
            cl = C.build_claims(outdir, a.label)
            if not a.quiet:
                print("\n[账本] claims.json：%d 条结论，%d 个参数未记录"
                      % (len(cl["claims"]), len(cl["missing_params"])))
        except Exception as e:
            print("[账本] 生成失败：%s" % e)

    if "10" in chosen:
        print("\n[10] 红线审计")
        run("10", "audit", ["--dir", outdir] + q)

    if "11" in chosen:
        print("\n[11] 生成 HTML 报告（含审计结果）")
        run("11", "report", ["--dir", outdir, "--label", a.label] + q)

    # ---- 可复现脚本：把本次真实命令行与解释器路径固化下来 ----
    argv = [x for x in sys.argv[1:]]
    py = sys.executable
    if sys.platform.startswith("win"):
        cmd = "@echo off\r\nREM 由 network-pharmacology-pro 自动生成；换机后可直接执行\r\n"
        cmd += "REM 若解释器路径不同，请把下一行的 python 换成你的解释器\r\n"
        cmd += "\"%s\" \"%s\" %s\r\n" % (py, os.path.join(HERE, "np_run.py"),
                                            " ".join('"%s"' % x if " " in x else x for x in argv))
        C.write_text(os.path.join(outdir, "reproduce.cmd"), cmd, encoding="utf-8")
    else:
        sh = "#!/usr/bin/env bash\n# 由 network-pharmacology-pro 自动生成\nset -euo pipefail\n"
        sh += '"%s" "%s" %s\n' % (py, os.path.join(HERE, "np_run.py"),
                                   " ".join('"%s"' % x if " " in x else x for x in argv))
        C.write_text(os.path.join(outdir, "reproduce.sh"), sh, encoding="utf-8")
    C.write_text(os.path.join(outdir, "reproduce.readme.md"),
                 "\n".join([
                     "# 如何复现本次结果", "",
                     "1. 安装 Python 3.8+（本项目核心零依赖；读 xlsx 需 `pip install openpyxl`）。",
                     "2. 确认 `environment.json` 中的平台与包版本与你的环境一致（不一致时先看差异项）。",
                     "3. 执行 `reproduce.cmd`（Windows）或 `reproduce.sh`（macOS/Linux）。",
                     "4. 用 `run_manifest.json` 中每个阶段的 `outputs.*.sha256_16` 与本次产物比对；",
                     "   哈希一致 = 完全复现；不一致时先查 environment.json 的包版本与 locale。", "",
                     "> 注意：涉及随机过程的阶段（零模型、bootstrap、布局）均使用固定种子，",
                     "> 因此**同一环境下应逐字节一致**；跨平台浮点差异可能导致末位不同。",
                 ]) + "\n")

    print("\n" + "=" * 64)
    print(" 完成。产物目录：%s" % outdir)
    idx = os.path.join(outdir, "np_report.html")
    if os.path.exists(idx):
        print(" 主报告：%s" % idx)
    print(" 审计：%s" % os.path.join(outdir, "audit", "audit_report.md"))
    print(" 溯源：%s" % os.path.join(outdir, "run_manifest.json"))
    print("=" * 64)
    return 0


if __name__ == "__main__":
    sys.exit(main())
