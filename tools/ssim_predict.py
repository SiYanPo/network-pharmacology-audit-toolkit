#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ssim_predict.py — SwissTargetPrediction 不可用时的外部靶点推断管线

为什么存在
----------
2026-11-22 实测：SwissTargetPrediction 的提交管道停摆——任何 SMILES（含阿司匹林）
都会被服务端模板拒绝：alert("Your YYYjobYYY will not submitted...")，占位符 YYYjobYYY
未被替换 = 后端 job 生成已停用；old.swisstargetprediction.ch 亦无法连接。
SuperPred(prediction.charite.de) 子页面全部 404，TargetNet 已被 ADMETlab 顶替。

本管线是**外部、独立于本仓库代码**的替代对照（DoD D4 语义不变：
用外部工具在同一金标准上产生预测排序，与本平台对比）：

  1) SwissSimilarity (SIB, 2021 版) ECFP4 相似性筛选，库 = ChEMBL actives only（85 万）；
     结果直接给 CHEMBL 分子 ID + Tanimoto 分数（无本地推断）。
  2) 每个成分取 top-N 相似分子，经 ChEMBL API 拿其活性记录（pChEMBL≥6）→ 靶点 ChEMBL ID；
  3) 经 UniProt REST 把靶点 ChEMBL ID → UniProt accession → 基因符号；
  4) 基因分数 = 支持它的命中里**最高**的相似性分数（可写论文的保守聚合）。

诚实声明（写报告时要带上）：
  * 这不是 SwissTargetPrediction——是"相似性筛选 + 已知活性注释"的靶点推断，
    原理与 STP 相同（2D 相似性到已知配体），但覆盖与概率口径不同；
  * 相似分子数（--top-hits）与活性阈值（pChEMBL）是两个自由参数，结论须做敏感性检查；
  * ChEMBL 命中可能恰好覆盖金标准库的靶点——对照结果解读时注意"预测来源"差异。

用法
    python tools/ssim_predict.py --input out/stp/stp_input.txt \
        --names "A,B,C" -o out/stp/ssim --top-hits 15
"""
import argparse
import csv
import json
import os
import re
import ssl
import sys
import time
import http.cookiejar
import urllib.parse
import urllib.request

BASE = "https://www.swisssimilarity.ch"
LIBRARY = "CHEMBL_act"          # ChEMBL (actives only), 849,502 分子
METHOD = "1"                    # 1 = ECFP4（结果页会回显 method 名，可核对）
CLASS_IDX = "2"                 # 2 = Bioactive
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
PCHEMBL_MIN = 6.0

UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126 Safari/537.36")


def _ctx():
    c = ssl.create_default_context()
    c.check_hostname = False
    c.verify_mode = ssl.CERT_NONE
    return c


def _opener():
    cj = http.cookiejar.CookieJar()
    return urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj),
                                       urllib.request.HTTPSHandler(context=_ctx()))


def _get(url, op=None, timeout=90, headers=None):
    h = {"User-Agent": UA}
    if headers:
        h.update(headers)
    req = urllib.request.Request(url, headers=h)
    if op is not None:
        r = op.open(req, timeout=timeout)     # cookie 会话（screen.php 用）
    else:
        r = urllib.request.urlopen(req, timeout=timeout, context=_ctx())  # 裸 GET
    return r.read().decode("utf-8", "ignore")


def screen(op, smiles):
    """提交一个分子到 SwissSimilarity，返回 (job_id, method_used, err)。"""
    body = urllib.parse.urlencode({
        "smiles": smiles, "compoundClasses": CLASS_IDX,
        "library_method": "%s___%s" % (LIBRARY, METHOD),
        "computationTime": ""}).encode()
    req = urllib.request.Request(BASE + "/screen.php", data=body, method="POST",
                                 headers={"Content-Type": "application/x-www-form-urlencoded",
                                          "Origin": BASE, "Referer": BASE + "/"})
    try:
        r = op.open(req, timeout=240)
        txt = r.read().decode("utf-8", "ignore")
    except Exception as e:
        return None, "", "NET:%s %s" % (type(e).__name__, str(e)[:80])
    m = re.search(r'location\.replace\("([^"]*result\.php\?job=(\d+)[^"]*)"\)', txt)
    if not m:
        # 库忙/拒绝时页面会带 alert
        a = re.search(r'alert\("([^"]+)"\)', txt)
        return None, "", (a.group(1)[:100] if a else "NO_JOB（响应里无 result.php 跳转）")
    mu = re.search(r'Screening method:\s*([^<\n]+)', txt)
    return m.group(2), (mu.group(1).strip() if mu else ""), ""


def parse_hits(txt):
    """结果页 → [{id, name, score, smiles}]。

    两种单元格格式（实测）：
      DrugBank 类： <a ...>DB08587</a>, <a ...>Sinapic acid<br></a> Score : 0.659
      CHEMBL 类：   <a ...>CHEMBL501235</a><br> Score : 0.523        （无名称）
    """
    hits = []
    for m in re.finditer(
            r'<a[^>]*>((?:DB|CHEMBL)\d+)</a>'
            r'(?:,\s*<a[^>]*>([^<]+?)<br>\s*</a>)?'   # 名称可选
            r'<br>\s*Score\s*:\s*([0-9.]+)', txt, re.S):
        hits.append({"id": m.group(1), "name": (m.group(2) or "").strip(),
                     "score": float(m.group(3))})
    # 去重（保留最高分）
    best = {}
    for h in hits:
        if h["id"] not in best or h["score"] > best[h["id"]]["score"]:
            best[h["id"]] = h
    return sorted(best.values(), key=lambda x: -x["score"])


def fetch_hits(op, job, timeout=120, poll=10, max_poll=240):
    url = "%s/result.php?job=%s" % (BASE, job)
    t0 = time.time()
    while time.time() - t0 < max_poll:
        try:
            txt = _get(url, op=op, headers={"Referer": BASE + "/screen.php"},
                       timeout=timeout)
        except Exception:
            time.sleep(poll)
            continue
        hits = parse_hits(txt)
        if hits:
            return hits, ""
        time.sleep(poll)
    return None, "TIMEOUT_NO_HITS"


def _json_cache(url, cache, key, timeout=90):
    p = os.path.join(cache, key)
    if os.path.exists(p):
        return json.load(open(p, encoding="utf-8"))
    try:
        txt = _get(url, timeout=timeout, headers={"Accept": "application/json"})
        j = json.loads(txt)
        json.dump(j, open(p, "w", encoding="utf-8"))
        time.sleep(0.35)                      # 遵守 API 礼仪
        return j
    except Exception as e:
        return {"__err": "%s:%s" % (type(e).__name__, str(e)[:80])}


def map_targets(chembl_id, cache):
    """CHEMBL 分子 ID → {uniprot_acc: gene_symbol}（pChEMBL≥6 的活性记录）。"""
    key = "mol_%s.json" % chembl_id
    j = _json_cache("%s/activity.json?molecule_chembl_id=%s&limit=1000" % (
        CHEMBL_API, chembl_id), cache, key)
    if "__err" in j:
        return {}, j["__err"]
    tids = set()
    for a in j.get("activities", []):
        pv = a.get("pchembl_value")
        t = a.get("target_chembl_id")
        if t and pv:
            try:
                if float(pv) >= PCHEMBL_MIN:
                    tids.add(t)
            except (TypeError, ValueError):
                pass
    out = {}
    errs = []
    for tid in sorted(tids):
        tj = _json_cache("%s/target/%s.json" % (CHEMBL_API, tid), cache,
                         "tgt_%s.json" % tid)
        if "__err" in tj:
            errs.append("%s:%s" % (tid, tj["__err"][:40]))
            continue
        # 只保留人源靶点：对照的金标准是 human 基因符号；
        # 混入鼠/大鼠条目会出现 ache/Grin1 这类小写符号，污染排序口径。
        if (tj.get("organism") or "") != "Homo sapiens":
            continue
        for comp in tj.get("target_components", []):
            acc = comp.get("accession") or ""
            if acc:
                out[acc] = comp.get("component_descriptions", [{}])[0].get(
                    "description", "") if comp.get("component_descriptions") else ""
                out["_name_" + acc] = tj.get("pref_name", "")
    return out, "; ".join(errs)


def gene_symbol(accession, cache):
    key = "uni_%s.json" % accession
    j = _json_cache("%s/%s.json" % (UNIPROT_API, accession), cache, key)
    if "__err" in j:
        return ""
    try:
        genes = j.get("genes") or []
        # UniProt REST 字段是 genes[0].geneName.value（不是 genePrimary）
        g = (genes[0].get("geneName") or {}).get("value") if genes else None
        if not g:
            g = (genes[0].get("genePrimary") or {}).get("value") if genes else None
        return g or ""
    except (IndexError, AttributeError, TypeError):
        return ""


def predict_one(op, name, smiles, outdir, cache, top_hits=15, top_targets=25):
    """一个成分 → gene,score CSV + 证据 JSON。"""
    job, method, err = screen(op, smiles)
    if not job:
        return {"name": name, "ok": False, "err": err}
    hits, err2 = fetch_hits(op, job)
    if not hits:
        return {"name": name, "ok": False, "err": err2 or "NO_HITS"}
    top = hits[:top_hits]
    prov = {"name": name, "smiles": smiles, "job": job, "method": method,
            "library": LIBRARY, "n_hits_total": len(hits),
            "pchembl_min": PCHEMBL_MIN, "top_hits": [], "errors": []}
    genes = {}
    for h in top:
        tg, errs = map_targets(h["id"], cache)
        if errs:
            prov["errors"].append(errs)
        row = {"id": h["id"], "name": h["name"], "score": h["score"],
               "n_targets": len([k for k in tg if not k.startswith("_")])}
        prov["top_hits"].append(row)
        for acc in tg:
            if acc.startswith("_"):
                continue
            g = gene_symbol(acc, cache)
            if not g:
                continue
            s = h["score"]
            if g not in genes or s > genes[g]["score"]:
                genes[g] = {"score": s, "via": h["id"], "uniprot": acc}
    ranked = sorted(genes.items(), key=lambda kv: -kv[1]["score"])[:top_targets]
    out_csv = os.path.join(outdir, "%s.csv" % name)
    with open(out_csv, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "score"])
        for g, d in ranked:
            w.writerow([g, round(d["score"], 4)])
    prov["n_genes"] = len(ranked)
    prov["genes"] = [{"gene": g, "score": d["score"], "via": d["via"]}
                     for g, d in ranked]
    json.dump(prov, open(os.path.join(outdir, "%s_provenance.json" % name), "w",
                         encoding="utf-8"), ensure_ascii=False, indent=1)
    return {"name": name, "ok": True, "n_hits": len(hits), "n_genes": len(ranked),
            "csv": out_csv}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", help="SMILES 文件（一行一个）")
    ap.add_argument("--smiles", help="单个 SMILES")
    ap.add_argument("--name", default="compound")
    ap.add_argument("--names", help="逗号分隔名字（与 --input 对应）")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--cache", default=os.path.expanduser("~/.np_cache/chembl"))
    ap.add_argument("--top-hits", type=int, default=15)
    ap.add_argument("--top-targets", type=int, default=25)
    ap.add_argument("--cooldown", type=int, default=20)
    a = ap.parse_args(argv)
    if a.smiles:
        names, smis = [a.name], [a.smiles]
    elif a.input:
        smis = [x.strip() for x in open(a.input, encoding="utf-8") if x.strip()]
        names = ([x.strip() for x in a.names.split(",")] if a.names
                 else ["c%d" % i for i in range(len(smis))])
    else:
        print("需要 --input 或 --smiles")
        return 2
    assert len(names) == len(smis)
    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(a.cache, exist_ok=True)
    op = _opener()
    ok = 0
    for i, (n, s) in enumerate(zip(names, smis)):
        r = predict_one(op, n, s, a.outdir, a.cache, a.top_hits, a.top_targets)
        print("  %s %-18s %s" % ("✓" if r["ok"] else "x", n,
                                 r.get("csv") or r["err"]))
        ok += bool(r["ok"])
        if i < len(names) - 1:
            time.sleep(a.cooldown)
    print("完成 %d/%d" % (ok, len(names)))
    return 0 if ok == len(names) else 1


if __name__ == "__main__":
    sys.exit(main())
