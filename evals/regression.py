#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
regression.py — 缺陷回归 + 契约测试（v1.0 → v1.1 修复项的守门测试）

原则：每个曾经出现过的缺陷，都要有一个**会在修复被回退时变红**的断言。
契约测试覆盖所有此前"从未被调用"的分支（--sets / --unit / --ref / 目录模式 / --strict 等）。

用法
  python evals/regression.py            # 全部
  python evals/regression.py --fast     # 跳过最慢的端到端
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
S = os.path.join(ROOT, "scripts")
EX = os.path.join(ROOT, "examples")
sys.path.insert(0, S)

import npcore as C  # noqa: E402

PASS = FAIL = 0
FAILS = []
TMP = tempfile.mkdtemp(prefix="np_reg_")


def check(cond, name, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
    else:
        FAIL += 1
        FAILS.append((name, detail))
        print("  ✗ %s %s" % (name, ("— " + detail) if detail else ""))


def run(script, *args, **kw):
    return subprocess.run([sys.executable, os.path.join(S, script)] + list(args),
                          capture_output=True, text=True, timeout=kw.get("timeout", 300),
                          cwd=ROOT)


def out(*p):
    return os.path.join(TMP, *p)


# ---------------------------------------------------------------- 缺陷回归
def bug1_methods_params():
    """v1.0：方法学段落会用默认值覆盖用户实际参数（学术诚信级）。"""
    print("[BUG-1] 方法学段落不得伪造参数")
    d = out("b1")
    r = run("np_compounds.py", "-i", os.path.join(EX, "compounds.csv"), "-o", d,
            "--ob", "45", "--dl", "0.25", "--route", "oral", "--quiet")
    check(r.returncode == 0, "02 阶段运行", r.stderr[-200:])
    r = run("np_report.py", "--dir", d, "--label", "回归", "--quiet")
    check(r.returncode == 0, "报告生成", r.stderr[-200:])
    h = open(os.path.join(d, "np_report.html"), encoding="utf-8").read()
    check("OB ≥ 45%" in h, "报告写入实际 OB=45（不是默认 30）")
    check("OB ≥ 30%" not in h, "报告中不出现被伪造的默认值")
    check("未记录" in h, "缺失参数被显式标为「未记录」")
    j = json.load(open(os.path.join(d, "claims.json"), encoding="utf-8"))
    check(any("ob" in (c.get("value") or {}) for c in j["claims"]), "账本记录 ob 数值")
    check(len(j["missing_params"]) > 0, "账本列出未记录参数", str(j["missing_params"][:2]))


def bug2_mcc():
    """v1.0：MCC 在稠密小图上退化为常数，抹平排序。"""
    print("[BUG-2] MCC 必须等于精确实现")
    g = C.Graph()
    for u, v in [("A", "B"), ("B", "C"), ("C", "A"), ("A", "D"), ("D", "E"), ("E", "A"), ("B", "D")]:
        g.add_edge(u, v, 1.0)
    m, info = g.mcc()
    ref = {"A": 6, "B": 4, "C": 2, "D": 4, "E": 2}
    check(all(abs(m[k] - ref[k]) < 1e-9 for k in ref), "MCC 精确值匹配", str(m))
    check(len(set(m.values())) > 1, "MCC 不再全部并列（有排序信息）")
    check(info["cliques"] >= 3 and info["max_size"] == 3, "极大团枚举信息完整", str(info))


def bug3_degree_null():
    """v1.0：C 组零模型按总体比例分箱，未匹配观察集度分布（反保守）。"""
    print("[BUG-3] 度匹配零模型必须真正匹配")
    h, rows = C.read_table(os.path.join(EX, "string_edges.tsv"))
    g = C.edges_from_rows(rows, "node1", "node2", "combined_score", 0.4)
    deg = g.degree()
    pool = sorted(g.nodes)
    rng = __import__("random").Random(7)
    obs = sorted(deg.values())[-25:]          # 取高度节点，模拟"交集偏向高度"
    samp, q = C.degree_matched_sample(obs, pool, deg, rng)
    check(q["mean_abs_ddeg"] <= 1.0, "平均|Δdeg| ≤1", "%.2f" % q["mean_abs_ddeg"])
    check(q["ks"] <= 0.10, "KS ≤0.10", "%.3f" % q["ks"])
    check(q["exact_frac"] >= 0.7, "精确匹配比例 ≥0.7", "%.2f" % q["exact_frac"])
    # 与"按总体比例分箱"对比：后者必然低估
    from collections import Counter
    cnt = Counter(deg.values())
    tot = sum(cnt.values())
    exp_mean = sum(d * c for d, c in cnt.items()) / tot
    check(abs(sum(deg[x] for x in samp) / len(samp) - sum(obs) / len(obs)) < 1.0,
          "抽样度均值≈观察集（不是总体均值 %.1f）" % exp_mean)


def bug4_audit_negation():
    """v1.0：审计器对否定句误报，且措辞规则与结构规则混算分。"""
    print("[BUG-4] 审计器：否定句不误报 + 结构/提示分离")
    d = out("b4")
    C.ensure_dir(d)
    open(os.path.join(d, "r.md"), "w", encoding="utf-8").write(
        "本文不将 P 值作为作用强度或重要性排名，仅作为显著性。\n"
        "未声称任何通路被激活。\n")
    r = run("np_audit.py", "--dir", d, "--quiet")
    j = json.load(open(os.path.join(d, "audit", "audit.json"), encoding="utf-8"))
    it = {x["id"]: x for x in j["items"]}
    check(it["R16"]["status"] == "PASS", "否定句不触发 R16", it["R16"]["detail"][:40])
    check(it["R17"]["status"] == "PASS", "否定句不触发 R17")
    check(it["R16"]["kind"] == "advisory", "R16 归类为提示性")
    check(any(x["kind"] == "structural" for x in j["items"]), "存在结构性规则")
    check(j["score"] < 50, "只有散文、没有运行产物 → 结构分应很低", str(j["score"]))


def bug5_inf_csv():
    """v1.1 新发现：写 CSV 时 inf 触发 OverflowError。"""
    print("[BUG-5] inf/NaN 不得让 CSV 写出崩溃")
    p = out("b5.csv")
    C.write_table(p, ["a", "b"], [{"a": float("inf"), "b": float("nan")},
                                  {"a": float("-inf"), "b": 1.0}])
    h, rows = C.read_table(p)
    check(rows[0]["a"] == "inf", "正无穷写成 inf", rows[0]["a"])
    check(rows[0]["b"] == "", "NaN 写成空")
    check(rows[1]["a"] == "-inf", "负无穷写成 -inf")


def bug6_ks():
    print("[BUG-6] KS 距离：相同样本必须为 0")
    check(C.ks_distance([1, 2, 2, 3], [3, 2, 2, 1]) == 0.0, "相同多重集 KS=0")
    check(abs(C.ks_distance([1, 2, 3], [1, 2, 4]) - 1 / 3) < 1e-9, "已知值 KS=1/3")
    check(C.ks_distance([1, 1], [5, 5]) == 1.0, "完全不重叠 KS=1")


def bug7_gates():
    """v1.1：未标准化输入、非法背景必须被阻断，而不是静默产出垃圾。"""
    print("[BUG-7] 输入门禁")
    d = out("b7")
    r = run("np_null.py", "-g", os.path.join(EX, "disease_targets.tsv"),
            "-e", os.path.join(EX, "string_edges.tsv"),
            "-d", os.path.join(EX, "disease_targets.tsv"),
            "--drug", os.path.join(EX, "compound_targets.csv"),
            "--ref", os.path.join(EX, "compound_targets.csv"), "-o", d, "-B", "20")
    check(r.returncode == 3 or "阻断" in (r.stdout + r.stderr), "D 组拒绝未标准化输入",
          str(r.returncode))

    gl = out("b7_genes.txt")
    open(gl, "w", encoding="utf-8").write("AKT1\nTNF\nIL6\nTP53\n")
    r2 = run("np_enrich.py", "-g", gl, "-m", os.path.join(EX, "pathways.gmt"),
             "--background", gl, "-o", out("b7e"))
    check(r2.returncode == 3, "背景集过小被阻断", str(r2.returncode))


def bug8_order():
    """v1.1：审计在报告之前跑，导致「存在账本」永远 FAIL（顺序依赖）。"""
    print("[BUG-8] 账本必须先于审计生成")
    d = out("b8")
    r = run("np_run.py", "--demo", "-o", d, "--quiet", "--iterations", "60", "--bootstrap", "30")
    check(r.returncode in (0, 1, 2), "全流程可完成", str(r.returncode))
    check(os.path.exists(os.path.join(d, "claims.json")), "claims.json 存在")
    j = json.load(open(os.path.join(d, "audit", "audit.json"), encoding="utf-8"))
    it = {x["id"]: x for x in j["items"]}
    check(it["R27"]["status"] == "PASS", "R27 账本项通过", it["R27"]["detail"])
    check(it["R28"]["status"] == "PASS", "R28 参数记录项通过", it["R28"]["detail"])
    check(j["structural"]["fail"] == 0, "结构性规则无 FAIL",
          str([x["id"] for x in j["items"] if x["status"] == "FAIL"]))


# ---------------------------------------------------------------- 契约测试
def contracts():
    print("[契约] 此前从未被调用的分支")
    # --logic or / --hl / --exceptions
    d = out("c1")
    r = run("np_compounds.py", "-i", os.path.join(EX, "compounds.csv"), "-o", d,
            "--logic", "or", "--hl", "4", "--exceptions", "Astragaloside IV", "--quiet")
    check(r.returncode == 0, "--logic or + --hl + --exceptions", r.stderr[-150:])
    check(os.path.getsize(os.path.join(d, "compounds_exception.csv")) > 60, "例外清单非空")

    # --keep-nonhuman
    d = out("c2")
    r = run("np_standardize.py", "-i", os.path.join(EX, "compound_targets.csv"),
            "-o", d, "--keep-nonhuman", "--quiet")
    check(r.returncode == 0, "--keep-nonhuman", r.stderr[-150:])
    h, rows = C.read_table(os.path.join(d, "targets_standardized.csv"))
    check(any(r.get("layer") == "nonhuman" for r in rows), "保留非人类条目（用于审计）")

    # --sets 模式
    a = out("c3a.txt")
    b = out("c3b.txt")
    open(a, "w").write("AKT1\nTNF\nIL6\nTP53\nCASP3\n")
    open(b, "w").write("AKT1\nTNF\nBCL2\nCASP3\nMAPK1\n")
    r = run("np_intersect.py", "--sets", a + "," + b, "-o", out("c3"))
    check(r.returncode == 0, "--sets 多组交集", r.stderr[-150:])
    h, rows = C.read_table(out("c3", "intersection.csv"))
    check(len(rows) == 3, "交集数正确（AKT1/TNF/CASP3）", str(len(rows)))

    # --unit compound
    d = out("c4s")
    run("np_standardize.py", "-i", os.path.join(EX, "compound_targets.csv"), "-o", d, "--quiet")
    r = run("np_optimize.py", "--drug", os.path.join(d, "targets_standardized.csv"),
            "--gene-col", "symbol", "--disease", os.path.join(EX, "disease_targets.tsv"),
            "--unit", "compound", "-o", out("c4"))
    check(r.returncode == 0, "--unit compound", r.stderr[-150:])

    # --largest-only
    r = run("np_ppi.py", "-e", os.path.join(EX, "string_edges.tsv"), "-o", out("c5"),
            "--largest-only", "--top", "8", "--bootstrap", "40")
    check(r.returncode == 0, "--largest-only", r.stderr[-150:])

    # 对接目录模式
    vd = out("c6", "AKT1")
    C.ensure_dir(vd)
    open(os.path.join(vd, "quercetin.log"), "w").write(
        "mode |   affinity\n   1       -8.4\n   2       -7.9\n")
    r = run("np_dock.py", "-i", out("c6"), "-o", out("c6o"))
    check(r.returncode == 0, "对接目录模式", r.stderr[-150:])
    h, rows = C.read_table(out("c6o", "docking_graded.csv"))
    check(len(rows) == 1 and abs(float(rows[0]["affinity"]) + 8.4) < 1e-6, "解析最佳打分")

    # --calibrate
    gl = out("c7.txt")
    open(gl, "w").write("\n".join(["AKT1", "TNF", "IL6", "TP53", "CASP3", "BCL2", "RELA",
                                    "MAPK1", "VEGFA", "TGFB1", "MMP9", "CCL2"]))
    r = run("np_enrich.py", "-g", gl, "-m", os.path.join(EX, "pathways.gmt"),
            "--calibrate", "60", "-o", out("c7"))
    check(r.returncode == 0, "--calibrate", r.stderr[-150:])
    hh, rows = C.read_table(out("c7", "enrich_significant.csv"))
    check("p_empirical" in hh, "输出经验 P 列", str(hh[:6]))

    # xlsx 输入（覆盖 _read_xlsx）
    try:
        import openpyxl
        xp = out("c8.xlsx")
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.append(["compound", "OB", "DL", "herb"])
        ws.append(["X1", 50, 0.3, "黄芪"])
        ws.append(["X2", 10, 0.1, "黄芪"])
        wb.save(xp)
        r = run("np_compounds.py", "-i", xp, "-o", out("c8o"), "--quiet")
        check(r.returncode == 0, "xlsx 输入（openpyxl 路径）", r.stderr[-150:])
    except ImportError:
        print("    （跳过 xlsx：未安装 openpyxl）")

    # 参数错误的友好提示（覆盖 Cli.error）
    r = subprocess.run([sys.executable, os.path.join(S, "np_compounds.py")],
                       capture_output=True, text=True, timeout=60, cwd=ROOT)
    check(r.returncode != 0, "缺参数时非零退出")
    check("参数错误" in (r.stderr + r.stdout) or "required" in (r.stderr + r.stdout),
          "给出可读的参数错误提示", (r.stderr or "")[-120:])

    # --strict 与 --from
    r = run("np_audit.py", "--dir", out("b8"), "--strict", "--quiet")
    check(r.returncode in (0, 1, 2), "--strict 可执行", str(r.returncode))
    r = run("np_run.py", "--demo", "-o", out("c9"), "--from", "04", "--dry-run")
    check(r.returncode == 0 and "04" in r.stdout, "--from + --dry-run", r.stdout[-200:])


def bugs_v12():
    """v1.2 真实数据运行暴露的缺陷回归。"""
    print("[BUG-9] 门控：evidence 不得剔除；auto 剔除已上市药物须告警")
    cmpd = os.path.join(TMP, "cmp_coh.csv")
    with open(cmpd, "w", encoding="utf-8") as f:
        f.write("compound,cohort,smiles,ob,dl\n")
        f.write("Quercetin,tcm,C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O,,\n")
        f.write("Atorvastatin,drug,CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O,,\n")
    d = out("v12a")
    r = run("np_compounds.py", "-i", cmpd, "-o", d, "--gate", "evidence", "--quiet")
    h, rows = C.read_table(os.path.join(d, "compounds_kept.csv"))
    check(len(rows) == 2, "evidence 门控不剔除任何有结构的化合物", str(len(rows)))
    d = out("v12b")
    r = run("np_compounds.py", "-i", cmpd, "-o", d, "--gate", "auto")
    h, kept = C.read_table(os.path.join(d, "compounds_kept.csv"))
    h2, rej = C.read_table(os.path.join(d, "compounds_rejected.csv"))
    check(len(rej) == 1 and rej[0]["compound"] == "Atorvastatin",
          "类药性门控确实会剔除阿托伐他汀（这是要告警的事实）")
    check("已上市药物" in (r.stdout or ""), "剔除已上市药物时必须告警", (r.stdout or "")[-160:])

    print("[BUG-10] 对接门禁必须读 RMSD 列，不能把靶点名当 RMSD")
    dk = os.path.join(TMP, "dock_v12.csv")
    with open(dk, "w", encoding="utf-8") as f:
        f.write("compound,target,affinity,redock_rmsd\n")
        f.write("Quercetin,HMGCR,-7.4,1.30\n")
        f.write("Curcumin,HMGCR,-6.0,1.30\n")
    d = out("v12c")
    r = run("np_dock.py", "-i", dk, "-o", d, "--quiet")
    txt = open(os.path.join(d, "docking_report.md"), encoding="utf-8").read()
    check("重对接验证通过 | 1" in txt.replace("  ", " "),
          "RMSD=1.30 的靶点应判为通过", [l for l in txt.splitlines() if "重对接" in l][:2])
    dk2 = os.path.join(TMP, "dock_v12b.csv")
    with open(dk2, "w", encoding="utf-8") as f:
        f.write("compound,target,affinity,redock_rmsd\n")
        f.write("Quercetin,DPP4,-9.1,3.06\n")
    d = out("v12d")
    run("np_dock.py", "-i", dk2, "-o", d, "--quiet")
    txt = open(os.path.join(d, "docking_report.md"), encoding="utf-8").read()
    check("重对接未通过 | 1" in txt.replace("  ", " "), "RMSD=3.06 应判未通过",
          [l for l in txt.splitlines() if "重对接" in l][:2])

    print("[BUG-11] 非 ASCII 路径下 RDKit 读写不得失败")
    import np_dockprep as D
    zh = os.path.join(TMP, "中文目录")
    os.makedirs(zh, exist_ok=True)
    sdf = os.path.join(zh, "x.sdf")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))
    AllChem.EmbedMolecule(m, randomSeed=42)
    D.mol_to_file(m, sdf)
    m2 = D.mol_from_file(sdf, remove_hs=True)
    check(m2 is not None and m2.GetNumAtoms() == 13, "中文路径 MOL 读写往返",
          str(m2.GetNumAtoms() if m2 else None))

    print("[BUG-12] STRING 0–1 分制不得被 int() 归零")
    g = C.edges_from_rows([{"a": "IL6", "b": "TNF", "s": "0.88"}], "a", "b", "s", 0.4)
    check(g.n == 2 and g.m == 1, "0.88 的边应被保留")
    w = list(g.adj["IL6"].values())[0]
    check(abs(w - 0.88) < 1e-9, "权重应保持 0.88", str(w))

    print("[BUG-13] 度分箱 z：方差为 0 时给 NA，不给 0")
    z, info = C.degree_binned_z({"a": 1.0, "b": 2.0, "c": 3.0},
                                {"a": 1, "b": 2, "c": 3}, n_bins=1, min_bin=2)
    check(any(v is not None for v in z.values()), "有方差时应能算出 z", str(z))
    z2, _ = C.degree_binned_z({"a": 1.0, "b": 1.0, "c": 1.0},
                              {"a": 1, "b": 1, "c": 1}, n_bins=1, min_bin=2)
    check(all(v is None for v in z2.values()), "零方差必须是 NA 而不是 0", str(z2))

    print("[BUG-14] resolve_col 必须认表头原名（organism 等不在别名表里）")
    check(C.resolve_col(["organism", "gene"], "organism") == "organism",
          "resolve_col 支持 key 自身")


def bug15_methods_gate():
    """方法学段落必须与实际门控一致（v1.2 真实运行暴露：用类药性门却写 OB≥30）。"""
    print("[BUG-15] 方法学段落不得把类药性门写成 OB/DL")
    cmpd = os.path.join(TMP, "cmp_gate.csv")
    with open(cmpd, "w", encoding="utf-8") as f:
        f.write("compound,cohort,smiles\n")
        for n, smi in [("Quercetin", "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O"),
                       ("Baicalein", "O=c1cc(-c2ccccc2)oc2cc(O)c(O)c(O)c12")]:
            f.write("%s,tcm,%s\n" % (n, smi))
    d = out("g15")
    r = run("np_compounds.py", "-i", cmpd, "-o", d, "--gate", "druglike", "--quiet")
    check(r.returncode == 0, "druglike 门可运行", r.stderr[-120:])
    r = run("np_report.py", "--dir", d, "--label", "门控测试")
    html = open(os.path.join(d, "np_report.html"), encoding="utf-8").read()
    check("QED ≥" in html or "QED &ge;" in html, "方法学段落应写类药性判据")
    check("OB ≥ 30%" not in html, "方法学段落不得出现未使用的 OB≥30%")


def bug16_reuse_validation_types():
    """--reuse-validation 从 CSV 读回的字符串不得让报告 %d/%.2f 崩溃（复现 Shan2026 暴露）"""
    print("[BUG-16] reuse-validation 的字符串字段须容错")
    import np_dockprep as D
    # CSV 读回的典型形状：全字符串，甚至空串
    row = {"target": "P38", "pdb": "6QYX", "ligand_code": "AUTO",
           "residues": "", "redock_affinity": "-7.40", "rmsd_mode1": "1.0",
           "rmsd_best": "0.8", "n_poses": "5", "passed": "PASS",
           "box": "1;2;3", "box_size": "20;20;20"}
    c = D.coerce_validation_row(row)
    check(isinstance(c["residues"], int) and c["residues"] == 0, "空 residues 应转 int 0", str(c["residues"]))
    check(abs(c["rmsd_mode1"] - 1.0) < 1e-9, "rmsd_mode1 字符串应转 float", str(c["rmsd_mode1"]))
    check(D.fmt_float(c["rmsd_mode1"]) == "1.00", "浮点格式化应为 1.00", D.fmt_float(c["rmsd_mode1"]))
    check(D.fmt_float("NA") == "NA", "非数值应格式化为 NA")
    # 报告格式化整行不得抛异常
    try:
        line = "| %s | %s | %s | %s | %s | %s | %s |" % (
            c["target"], c["pdb"], c["ligand_code"], c["residues"],
            c["redock_affinity"], D.fmt_float(c["rmsd_mode1"]), c["passed"])
        check("P38" in line and "1.00" in line, "复用行格式化不得崩溃", line)
    except Exception as e:
        check(False, "复用行格式化不得崩溃", repr(e))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true")
    a = ap.parse_args()
    print("=" * 60)
    print(" 缺陷回归 + 契约测试")
    print("=" * 60)
    bug1_methods_params()
    bug2_mcc()
    bug3_degree_null()
    bug4_audit_negation()
    bug5_inf_csv()
    bug6_ks()
    bug7_gates()
    bugs_v12()
    bug15_methods_gate()
    bug16_reuse_validation_types()
    if not a.fast:
        bug8_order()
        contracts()
    shutil.rmtree(TMP, ignore_errors=True)
    print("-" * 60)
    print(" 通过 %d · 失败 %d" % (PASS, FAIL))
    if FAILS:
        for n, d in FAILS:
            print("  · %s %s" % (n, d))
        return 1
    print(" 全部通过 ✅")
    return 0


if __name__ == "__main__":
    sys.exit(main())
