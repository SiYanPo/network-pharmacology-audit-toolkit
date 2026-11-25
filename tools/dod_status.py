#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
dod_status.py — DoD 就绪度自检（评审 §13 的机器化版本）

为什么需要：DoD 的每一条都容易"以为做完了"。这个脚本只认**磁盘上的证据**，
不看你的记忆：没文件就是没做，文件里字段不对就是没通过，并直接给出下一条要跑的命令。

用法
----
    python tools/dod_status.py --dir D:/WorkBuddy/分身文件/<proj>/DoD收尾
    python tools/dod_status.py --dir np_out --snapshot snapshot_20261122.tar.gz

判据（每条都可反驳）
--------------------
D1 校准：`09_evidence/calibration_model.json` 的 `cv.ece < 0.10`
D2 多重性：`09_evidence/multiplicity_statement.json` 存在且 `n_hypotheses > 0`
D3 真实 MR：任一 `mr_result.json` 满足 `clump.mode == plink_clump` + `status == OK`
            + `exploratory_only == false` + `F_statistic >= 10`
D4 外部工具：任一 `external_comparison*.json` 里出现**非内部**工具名（不以 M3_/Hub/B0/B1 开头）
D5 离线复现：存在 `versions.json`（或 manifest.versions）且快照 tarball 可 `verify`
D6 契约：`np_contracts.check_dir` 为 PASS（若给了产物目录）
D7 容器：本机存在 docker 可执行文件（只证明"能构建"，不证明"构建过"）
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if os.path.join(ROOT, "scripts") not in sys.path:
    sys.path.insert(0, os.path.join(ROOT, "scripts"))

INTERNAL_TOOLS = ("M3_", "Hub", "B0", "B1", "NP_pro", "support", "penalized")
# 合成/演示数据不算"外部工具对照"——否则跑个假 CSV 就能自称 DoD 达成
NON_REAL = ("synthetic", "demo", "fake", "test", "mock", "stub", "placeholder")


def _j(path):
    try:
        return json.load(open(path, encoding="utf-8"))
    except Exception:
        return {}


def d1_calibration(outdir):
    p = os.path.join(outdir, "09_evidence", "calibration_model.json")
    if not os.path.exists(p):
        return False, "无 calibration_model.json", "python scripts/np_evidence_ledger.py fit --ledger <L> --gold <g.txt> -o <out>/09_evidence"
    ece = (_j(p).get("cv") or {}).get("ece")
    if ece is None:
        return False, "校准器缺 OOF 指标", "重跑 fit（会算 5 折 OOF）"
    return (float(ece) < 0.10), "OOF ECE=%.4f" % float(ece), ""


def d2_multiplicity(outdir):
    p = os.path.join(outdir, "09_evidence", "multiplicity_statement.json")
    if not os.path.exists(p):
        return False, "无 multiplicity_statement.json", "python scripts/np_fdr.py run --dir <out>"
    n = _j(p).get("n_hypotheses") or 0
    return (n > 0), "假设数=%s" % n, ("" if n else "python scripts/np_fdr.py run --dir <out>")


def _fdr_lookup(root):
    """从 fdr_table.csv 读 {gene: q_family}，供 mr_result.json 未带 q 时补。"""
    import csv
    for p in glob.glob(os.path.join(root, "**", "fdr_table.csv"), recursive=True):
        out = {}
        try:
            with open(p, encoding="utf-8-sig", errors="replace") as f:
                for r in csv.DictReader(f):
                    if r.get("family") == "MR_IVW":
                        out[(r.get("unit") or "").upper()] = r.get("q_family")
        except Exception:
            continue
        return out
    return {}


def d3_real_mr(root):
    qtab = _fdr_lookup(root)
    best = None
    for p in glob.glob(os.path.join(root, "**", "mr_result.json"), recursive=True):
        j = _j(p)
        if j.get("status") != "OK":
            continue
        clump = j.get("clump") or {}
        # 真 LD 的两种合法来源：本地 PLINK bfile，或 OpenGWAS /ld/matrix（1000G r²）。
        # 距离剪枝不算（exploratory_only 会为 true）。
        mode = clump.get("mode")
        real_ld = mode in ("plink_clump", "ld_matrix")
        # 注意：纯 Python 的 clump() 不返回 status（只有 plink_clump 会），
        # 所以这里把"缺 status"视为成功——否则 ld_matrix 永远判不通过（已踩过一次）。
        clump_ok = clump.get("status") in ("OK", None)
        ok = (real_ld and clump_ok
              and j.get("exploratory_only") is False
              and (j.get("F_statistic") or 0) >= 10)
        item = {"file": p, "ok": ok, "f": j.get("F_statistic"), "clump": mode,
                "exploratory": j.get("exploratory_only"),
                "ld": (j.get("ld_source") or {}).get("source", ""),
                "q": (j.get("fdr") or {}).get("q_value")
                     or qtab.get(os.path.basename(os.path.dirname(p)).upper())}
        if ok:
            return True, "F=%.1f LD=%s q=%s（%s）" % (
                item["f"] or 0, mode, item["q"], os.path.basename(os.path.dirname(p))), ""
        best = best or item
    if best:
        # 注意：这里必须是三个返回值——之前写成 (msg, cmd) 的嵌套元组，
        # 只有 2 个元素，调用方解包 3 个会 ValueError（已由 test_tools_dod 抓到）。
        return (False,
                "找到 MR 但未达标：clump=%s exploratory=%s F=%s" % (
                    best["clump"], best["exploratory"], best["f"]),
                "加 --bfile <g1000_eur> --n-exposure/--n-outcome 重跑")
    return False, "没有任何 mr_result.json", "python scripts/np_mr.py fetch --exposure <eqtl-id> --outcome <gwas-id> ..."


def d4_external(root):
    # 扫描**全部** external_comparison*.json：glob 顺序不定，先撞上合成对照
    # 就返回会让真实外部对照（stp_yz/compare/…）被永远掩埋。
    files = sorted(glob.glob(os.path.join(root, "**", "external_comparison*.json"),
                             recursive=True))
    best_real, fallback = None, None
    for p in files:
        j = _j(p)
        tools = [t.get("tool", "") for t in (j.get("tools") or [])]
        ext = [t for t in tools
               if not any(t.startswith(x) or t == x for x in INTERNAL_TOOLS)
               and not any(k in t.lower() for k in NON_REAL)]
        if ext and best_real is None:
            best_real = (p, ext, tools)
        if tools and fallback is None:
            fallback = (p, tools)
    if best_real:
        p, ext, tools = best_real
        return True, "外部工具：%s（%s）" % (",".join(ext), os.path.basename(os.path.dirname(p))), ""
    if fallback:
        p, tools = fallback
        return False, "只有内部/合成基线：%s（合成数据不计入）" % ",".join(tools), (
            "tool/stp_convert.py 转换 STP 导出后，用 np_split compare --pred 灌入")
    return False, "没有 external_comparison.json", "python scripts/np_split.py compare --gold-csv <g> --pred <name>=<file>"


def d5_snapshot(root, snapshot=None):
    if snapshot:
        try:
            r = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "np_snapshot.py"),
                                "verify", snapshot], capture_output=True, text=True, timeout=180)
            if r.returncode == 0:
                return True, "快照校验通过：%s" % os.path.basename(snapshot), ""
            return False, "快照校验失败：%s" % (r.stdout or r.stderr)[-80:], "重跑 np_snapshot.py pack"
        except Exception as e:
            return False, "verify 调用失败：%s" % type(e).__name__, ""
    for p in glob.glob(os.path.join(root, "**", "versions.json"), recursive=True) + \
            glob.glob(os.path.join(root, "**", "pipeline_manifest.json"), recursive=True):
        j = _j(p)
        v = j.get("versions") if "versions" in j else j
        if isinstance(v, dict) and (v.get("STRING") or v.get("Ensembl")):
            return True, "版本已记录：STRING=%s Ensembl=%s" % (v.get("STRING"), v.get("Ensembl")), ""
    return False, "未记录数据库版本", "python scripts/np_snapshot.py probe --save"


def d6_contracts(outdir):
    try:
        import np_contracts as CT
        r = CT.check_dir(outdir)
        return (r["verdict"] == "PASS"), "%s（硬失败 %d）" % (r["verdict"], r["n_fail"]), \
            ("" if r["verdict"] == "PASS" else "python scripts/np_contracts.py check --dir <out>")
    except Exception as e:
        return False, "契约模块不可用：%s" % type(e).__name__, ""


def d7_docker():
    # 2026-11-22 实测教训：Docker Desktop 装好了但 CLI 不在当前 shell 的 PATH 里，
    # 且引擎只在 Desktop 应用启动后可用（服务 Running != engine ready）。
    # 所以检查三件事：CLI 位置 → daemon 是否应答 → 服务是否安装。
    w = shutil.which("docker") or shutil.which("docker.exe")
    if not w:
        cand = os.path.join(os.environ.get("ProgramFiles", r"C:\Program Files"),
                            "Docker", "Docker", "resources", "bin", "docker.exe")
        if os.path.exists(cand):
            w = cand
    if w:
        try:
            r = subprocess.run([w, "info", "--format", "{{.ServerVersion}}"],
                               capture_output=True, text=True, timeout=20)
            if r.returncode == 0 and r.stdout.strip():
                # 2026-11-23 增强：引擎可用只证明"能构建"，不等于"构建过"。
                # 若本机已有 np-pro 镜像，一并核实，让 D7 反映真实完成度。
                img_detail = ""
                try:
                    ii = subprocess.run(
                        [w, "image", "inspect", "np-pro:3.0",
                         "--format", "{{.Id}} {{.Size}}"],
                        capture_output=True, text=True, timeout=15)
                    if ii.returncode == 0 and ii.stdout.strip():
                        img_detail = "；镜像 np-pro:3.0 已构建（%s）" % ii.stdout.strip()
                except Exception:
                    pass
                return True, ("docker 引擎可用：版本 %s（%s）%s"
                              % (r.stdout.strip(), w, img_detail)), ""
            return False, ("docker CLI 存在但引擎未启动（Desktop 应用没开）：请先启动 Docker Desktop，"
                           "等托盘鲸鱼图标变绿"), "启动 Docker Desktop 后重跑本面板"
        except Exception as e:
            return False, "docker 调用失败：%s" % str(e)[:60], "检查 Docker Desktop 是否启动"
    svc = os.path.join(os.environ.get("SystemRoot", r"C:\Windows"),
                       "System32", "drivers", "etc")  # 仅占位，不做无谓探测
    return False, "本机无 docker CLI（只能证明\"没验证过\"，不是\"不通过\"）", \
        "安装 Docker Desktop 后重跑"


def main(argv=None):
    ap = argparse.ArgumentParser(description="dod_status: DoD 就绪度自检")
    ap.add_argument("--dir", required=True,
                    help="MR / 外部对照证据目录（含 mr_result.json、stp_ranking.csv 那一层）")
    ap.add_argument("--np-out", default=None,
                    help="流水线产物目录（含 09_evidence/、fdr_table.csv、finalcheck.json）；"
                         "缺省时用 --dir（两种证据在同一目录的老用法）")
    ap.add_argument("--snapshot", default=None, help="快照 tarball（可选）")
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    d = os.path.abspath(a.dir)
    # 证据分家：校准/FDR/契约看流水线产物；MR/外部对照看 --dir。
    # 不分家会导致 FDR/契约去查 MR 目录（那里本来就不该有 multiplicity_statement.json），
    # 误报 9 条契约硬失败——D6 检查的是"流水线该有的产物"，不是"任何目录都该有"。
    np_out = os.path.abspath(a.np_out) if a.np_out else d
    rows = []
    for code, fn in (("D1", lambda: d1_calibration(np_out)), ("D2", lambda: d2_multiplicity(np_out)),
                     ("D3", lambda: d3_real_mr(d)), ("D4", lambda: d4_external(d)),
                     ("D5", lambda: d5_snapshot(np_out, a.snapshot)), ("D6", lambda: d6_contracts(np_out)),
                     ("D7", d7_docker)):
        ok, detail, nxt = fn()
        rows.append({"code": code, "ok": ok, "detail": detail, "next": nxt})
    names = {"D1": "校准 P+CI（ECE<0.10）", "D2": "全局 FDR + 声明", "D3": "真实 MR（真 LD + F≥10）",
             "D4": "外部工具对照", "D5": "版本快照/离线复现", "D6": "产物契约",
             "D7": "容器可构建（本机检测）"}
    print("=" * 72)
    print(" DoD 就绪度：MR/外部对照=%s%s" % (d, ("" if np_out == d else "  流水线产物=%s" % np_out)))
    print("=" * 72)
    for r in rows:
        print("  [%s] %-26s %s" % ("✓" if r["ok"] else "·", names[r["code"]], r["detail"]))
        if not r["ok"] and r["next"]:
            print("        下一步：%s" % r["next"])
    n = sum(r["ok"] for r in rows)
    print("\n  合计 %d/%d；未达标项：%s" % (
        n, len(rows), ", ".join(r["code"] for r in rows if not r["ok"]) or "无"))
    if a.json:
        json.dump({"dir": d, "n_ok": n, "results": rows},
                  open(a.json, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return 0 if n == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
