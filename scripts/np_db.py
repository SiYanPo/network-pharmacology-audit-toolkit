#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_db.py —— v2.0.0 数据库与在线资源登记表 + 连通性体检

为什么需要：500 篇复现中 48% 的关键数字藏在"需要登录/需要注册"的库里（GeneCards、DrugBank、
OMIM、SwissTarget、HERB、ETCM…）。把"哪些库免费/需注册/需机构 IP/已失效"写进 skill，
避免每次重复踩坑。

命令：
  list                      打印登记表（按类别）
  check [--json out.json]   逐个探测可达性（HTTP 状态/超时），标注 login/api/状态
  guide                     打印需要注册的库与申请方式
"""
import argparse
import json
import os
import time

import requests

# 字段：name, url, kind(db/api/web), login(0 免费无需注册 / 1 需注册 / 2 需机构IP), note
DBS = [
    # —— 中药与成分 ——
    ("TCMSP", "https://old.tcmsp-e.com/tcmsp.php", "web", 0, "成分/靶点主库；无 API，脚本抓取(本 skill np_fetch)"),
    ("TCMSP-new", "https://tcmsp-e.com/", "web", 1, "新版需注册，老版可用"),
    ("HERB", "http://herb.ac.cn/", "web", 0, "本草组鉴；批量数据需手工下载后 import"),
    ("ETCM", "http://www.tcmip.cn/ETCM2/front/", "web", 0, "中医百科全书 2.0"),
    ("SymMap", "http://www.symmap.org/", "web", 0, "证-症-方-药映射"),
    ("BATMAN-TCM", "http://bionet.ncpsb.org.cn/batman-tcm/", "web", 1, "需邮箱提交任务"),
    ("TCMID", "http://www.megabionet.org/tcmid/", "web", 0, "长期不稳定"),
    ("PubChem", "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/cids/JSON", "api", 0, "结构/属性；PUG REST 免费"),
    ("ChEMBL", "https://www.ebi.ac.uk/chembl/api/data/molecule/search?q=quercetin&format=json", "api", 0, "活性数据；REST 免费"),
    ("BindingDB", "https://www.bindingdb.org/", "web", 1, "亲和力数据；批量需申请"),
    ("ZINC", "https://zinc.docking.org/", "web", 0, "可购化合物库"),
    ("SwissADME", "http://www.swissadme.ch/", "web", 0, "ADME 在线；仅网页"),
    ("ProTox-3", "https://tox.charite.de/protox3/", "web", 0, "毒性预测；网页"),
    ("admetSAR", "http://lmmd.ecust.edu.cn/admetsar2/", "web", 0, "ADMET 在线"),
    # —— 靶点与疾病 ——
    ("SwissTargetPrediction", "http://www.swisstargetprediction.ch/", "web", 0, "反向靶点预测；无公开 API（复现差异主因之一）"),
    ("STITCH", "http://stitch.embl.de/", "web", 0, "化合物-蛋白互作"),
    ("PharmMapper", "http://www.lilab-ecust.cn/pharmmapper/", "web", 0, "反向药效团"),
    ("UniProt", "https://rest.uniprot.org/uniprotkb/search?query=STAT3&format=json&size=1", "api", 0, "蛋白注释；REST 免费"),
    ("RCSB PDB", "https://files.rcsb.org/download/1UBQ.pdb", "api", 0, "结构下载"),
    ("AlphaFold DB", "https://alphafold.ebi.ac.uk/api/prediction/P04637", "api", 0, "预测结构"),
    ("STRING", "https://string-db.org/api/json/network?identifiers=TP53&species=9606", "api", 0, "PPI；REST 免费"),
    ("GeneCards", "https://www.genecards.org/", "web", 0, "疾病靶点；反爬严格，批量需人工/机构"),
    ("OMIM", "https://www.omim.org/", "web", 1, "遗传病；需注册 API key"),
    ("DisGeNET", "https://www.disgenet.org/", "web", 1, "基因-疾病关联；需注册"),
    ("Open Targets", "https://api.platform.opentargets.org/api/v4/graphql", "api", 0, "疾病靶点(本 skill 用此替代 GeneCards)；POST GraphQL"),
    ("CTD", "https://ctdbase.org/", "web", 1, "化学-基因-疾病；批量需申请"),
    ("KEGG", "https://rest.kegg.jp/list/pathway/hsa", "api", 0, "通路；REST 免费(限速)"),
    ("Reactome", "https://reactome.org/ContentService/data/query/R-HSA-109582", "api", 0, "通路"),
    ("GO", "http://geneontology.org/", "web", 0, "本体；可 API"),
    ("DAVID", "https://david.ncifcrf.gov/", "web", 0, "富集；网页提交"),
    ("Metascape", "https://metascape.org/", "web", 0, "富集；网页提交"),
    # —— 组学与临床 ——
    ("GEO", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSE14520&retmode=json", "api", 0, "表达谱；E-utilities/FTP 免费"),
    ("GDC/TCGA", "https://api.gdc.cancer.gov/projects?size=1", "api", 0, "TCGA 原始/处理数据"),
    ("cBioPortal", "https://www.cbioportal.org/api/studies?pageSize=1", "api", 0, "生存/突变；本 skill np_survival 用"),
    ("GTEx", "https://gtexportal.org/api/v2/dataset/tissueSiteDetail", "api", 0, "正常组织表达"),
    ("Human Protein Atlas", "https://www.proteinatlas.org/api/search_download.php?search=STAT3&format=json&columns=g,gd", "api", 0, "蛋白定位/表达；search_download API"),
    ("CELLxGENE", "https://api.cellxgene.cziscience.com/curation/v1/collections", "api", 0, "单细胞图谱"),
    ("DrugBank", "https://go.drugbank.com/", "web", 1, "药物金标准；学术需注册申请"),
    # —— 计算引擎 ——
    ("AutoDock Vina", "https://github.com/ccsb-scripps/AutoDock-Vina", "tool", 0, "路径见 `python scripts/np_tools.py resolve vina`；仓库可达即算 ok"),
    ("GROMACS", "https://www.gromacs.org/", "tool", 0, "本机 WSL /opt/md/envs/gmx"),
    ("OpenMM", "https://openmm.org/", "tool", 0, "路径见 `python scripts/np_tools.py resolve python_md`"),
    ("gmx_MMPBSA", "https://valdes-tresanco-ms.github.io/gmx_MMPBSA/", "tool", 0, "本机 WSL /opt/md/envs/mmpbsa"),
    ("Cytoscape", "https://cytoscape.org/", "tool", 0, "本 skill 输出 XGMML/SIF 供导入"),
    ("PyMOL", "https://pymol.org/", "tool", 0, "开源版可用"),
]


def cmd_list(a):
    hdr = f"{'名称':<24}{'类型':<6}{'登录':<6}说明"
    print(hdr); print("-" * 100)
    for n, u, k, lg, note in DBS:
        lab = {0: "免费", 1: "注册", 2: "机构"}[lg]
        print(f"{n:<24}{k:<6}{lab:<6}{note}")


def cmd_guide(a):
    print("需要注册/申请的库（按优先级）：\n")
    for n, u, k, lg, note in DBS:
        if lg >= 1:
            print(f"- {n}: {u}\n    说明: {note}")
    print("\n提示：GeneCards 可用 Open Targets 替代（本 skill 已内置，免注册且可复算）。")


def cmd_check(a):
    rows = []
    for n, u, k, lg, note in DBS:
        t0 = time.time()
        status, code = "unknown", None
        try:
            if "platform.opentargets.org" in u:
                r = requests.post(u, json={"query": "{ __typename }"}, timeout=15,
                                  headers={"User-Agent": "np-skill-dbcheck"})
                code = r.status_code
            elif "ebi.ac.uk/chembl" in u:
                r = requests.get(u, timeout=40, headers={"User-Agent": "np-skill-dbcheck"})
                code = r.status_code
            elif "github.com" in u:
                r = requests.get(u, timeout=20, headers={"User-Agent": "Mozilla/5.0"}, allow_redirects=True)
                code = 200 if r.status_code < 400 else r.status_code
            elif u.endswith(".pdb") or "/api/" in u or "rest." in u or "api." in u or u.endswith("JSON"):
                r = requests.get(u, timeout=15, headers={"User-Agent": "np-skill-dbcheck"})
                code = r.status_code
            else:
                r = requests.get(u, timeout=12, headers={"User-Agent": "Mozilla/5.0"})
                code = r.status_code
            status = "ok" if code < 400 else ("blocked" if code in (401, 403, 429) else f"http_{code}")
        except requests.exceptions.Timeout:
            status = "timeout"
        except requests.exceptions.ConnectionError:
            status = "conn_error"
        except Exception as e:
            status = type(e).__name__[:18]
        rows.append({"name": n, "url": u, "kind": k, "login": lg, "status": status,
                     "http": code, "ms": int((time.time() - t0) * 1000), "note": note})
        print(f"{n:<24}{status:<12}{str(code):<6}{rows[-1]['ms']:>6} ms")
    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"\n可达 {ok}/{len(rows)}")
    out = os.path.abspath(a.json or "db_status.json")
    json.dump(rows, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"→ {out}")


def main():
    ap = argparse.ArgumentParser(description="np_db: 数据库登记与连通性体检")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    sub.add_parser("guide")
    s = sub.add_parser("check"); s.add_argument("--json", default=None)
    a = ap.parse_args()
    {"list": cmd_list, "guide": cmd_guide, "check": cmd_check}[a.cmd](a)


if __name__ == "__main__":
    main()
