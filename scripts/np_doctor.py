#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_doctor.py — 环境与数据源自检（含"能力矩阵"）

为什么需要
  网络药理学脚本失败，八成不是算法错，而是：依赖没装、数据库连不上、路径含中文/空格、
  表头是中文、xlsx 缺 openpyxl。先花 20 秒自检，比在流程中间报错省一小时。

输出
  * 依赖与解释器状态
  * 在线端点可达性（不可达 = 该数据源"不可用"，**不等于阴性结果**）
  * 能力矩阵：哪些功能现在就能跑，哪些需要装包/联网
"""
import argparse
import importlib
import os
import socket
import sys
import tempfile
import urllib.request

OK, BAD, WARN = "✅", "❌", "⚠"

ENDPOINTS = [
    ("STRING", "https://string-db.org/api/tsv/version"),
    ("g:Profiler", "https://biit.cs.ut.ee/gprofiler/api/util/organisms_list/"),
    ("UniProt", "https://rest.uniprot.org/uniprotkb/search?query=AKT1&format=tsv&size=1"),
    ("RCSB PDB", "https://files.rcsb.org/download/1TQN.pdb"),
    ("PubChem", "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/quercetin/cids/TXT"),
    ("PubMed", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?retmode=json"),
    ("Europe PMC", "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=test&format=json&pageSize=1"),
    ("TCMSP", "https://tcmsp-e.com/"),
]

DEPS = [
    ("numpy", "可选：加速统计", False),
    ("openpyxl", "可选：读取 xlsx", False),
    ("networkx", "可选：交叉验证图算法", False),
    ("matplotlib", "可选：PNG 出图（SVG 已内置）", False),
    ("rdkit", "可选：结构警示筛查", False),
]

CAPS = [
    ("成分筛选（OB/DL/途径一致性）", "np_compounds.py", [], True),
    ("靶点标准化（内置词典）", "np_standardize.py", [], True),
    ("靶点标准化（UniProt 联网校验）", "np_standardize.py --online", ["UniProt"], False),
    ("交集 + 成分级富集", "np_intersect.py", [], True),
    ("PPI 拓扑 + hub 稳健性", "np_ppi.py", [], True),
    ("四套零模型", "np_null.py", [], True),
    ("富集（本地 GMT）", "np_enrich.py -m x.gmt", [], True),
    ("富集（g:Profiler 在线）", "np_enrich.py --gprofiler", ["g:Profiler"], False),
    ("剂量可达性", "np_dose.py", [], True),
    ("对接结果门禁", "np_dock.py", [], True),
    ("反向网络药理学", "np_optimize.py", [], True),
    ("项目红线审计", "np_audit.py", [], True),
    ("读取 xlsx 输入", "任意脚本", ["openpyxl"], False),
]


def check_endpoint(url, timeout=8):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "np-pro-doctor/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status < 400, r.status
    except Exception as e:
        return False, type(e).__name__


def main(argv=None):
    ap = argparse.ArgumentParser(description="环境与数据源自检")
    ap.add_argument("--online", action="store_true", help="检测在线端点（默认开启，--offline 可关闭）")
    ap.add_argument("--offline", action="store_true", help="只做离线检查")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    print("=" * 62)
    print(" 网络药理学研究平台 · 环境自检  v%s" % "1.0.0")
    print("=" * 62)

    print("\n[解释器]")
    print("  Python %s  (%s)" % (sys.version.split()[0], sys.platform))
    print("  可执行文件：%s" % sys.executable)
    major, minor = sys.version_info[:2]
    print("  版本要求 ≥3.8：%s" % (OK if (major, minor) >= (3, 8) else BAD))

    print("\n[可选依赖]")
    have = {}
    for mod, desc, need in DEPS:
        try:
            m = importlib.import_module(mod)
            have[mod] = True
            ver = getattr(m, "__version__", "")
            print("  %s %-12s %s %s" % (OK, mod, ver, "" if not desc else "— " + desc))
        except Exception:
            have[mod] = False
            print("  %s %-12s 未安装 — %s（缺失会自动降级）" % (WARN, mod, desc))

    print("\n[读写与路径]")
    try:
        d = tempfile.mkdtemp(prefix="np_doctor_")
        p = os.path.join(d, "测试 文件.txt")
        open(p, "w", encoding="utf-8").write("ok")
        ok = open(p, encoding="utf-8").read() == "ok"
        print("  %s 中文与空格路径读写" % (OK if ok else BAD))
        import shutil
        shutil.rmtree(d, ignore_errors=True)
    except Exception as e:
        print("  %s 中文与空格路径：%s" % (BAD, e))
    print("  当前工作目录：%s" % os.getcwd())

    online = {}
    if not a.offline:
        print("\n[在线数据源]（不可达 = 该源不可用，**不是**阴性结果）")
        for name, url in ENDPOINTS:
            ok, info = check_endpoint(url)
            online[name] = ok
            print("  %s %-12s %s" % (OK if ok else BAD, name, "" if ok else "不可达（%s）" % info))
        if not any(online.values()):
            print("  %s 全部不可达：请检查网络/代理；离线功能仍然可用。" % WARN)

    print("\n[能力矩阵]")
    for name, cmd, needs, offline_ok in CAPS:
        if offline_ok:
            print("  %s %-34s %s" % (OK, name, cmd))
        else:
            missing = [n for n in needs if (n in online and not online[n]) or
                       (n in have and not have[n])]
            if missing:
                print("  %s %-34s 需要：%s（当前不可用）" % (WARN, name, "、".join(needs)))
            else:
                print("  %s %-34s %s" % (OK, name, cmd))

    print("\n[建议下一步]")
    if not have.get("openpyxl"):
        print("  · 如需读 xlsx：python -m pip install openpyxl")
    if not a.offline and not online.get("STRING"):
        print("  · STRING 不可达：改用本地导出的边表（STRING 网页 → Export → TSV）")
    if not a.offline and not online.get("g:Profiler"):
        print("  · g:Profiler 不可达：用本地 GMT 富集（MSigDB/KEGG 均可导出 GMT）")
    print("  · 首次使用建议跑：python np_run.py --demo  （全流程 1 分钟，离线）")
    print("=" * 62)
    return 0


if __name__ == "__main__":
    sys.exit(main())
