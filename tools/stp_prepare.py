#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
stp_prepare.py — 为 SwissTargetPrediction 准备输入（名称 → SMILES）

背景（DoD 收尾第 2 项）
-----------------------
`np_split compare` 需要一个**外部工具**的靶点排序表才算"外部工具对照"。
SwissTargetPrediction（STP）网页只接受 **SMILES**，而 TCMSP 导出的成分表只有名称。
本脚本把"名称 → SMILES"这一步自动化（走 PubChem PUG REST，无需注册），
并把无法解析的成分**如实列出**，不静默跳过。

用法
----
    python tools/stp_prepare.py --compounds examples/compounds.csv --top 5 -o out/stp
    python tools/stp_prepare.py --names "Quercetin,Tanshinone IIA" -o out/stp
    python tools/stp_prepare.py --compounds c.csv --top 5 -o out/ --offline   # 只做表，不联网

产出
----
    stp_input.txt       每行一个 SMILES（直接粘贴到 STP 输入框）
    stp_compounds.csv   name, smiles, pubchem_cid, source（含失败清单）
    stp_manifest.json   解析成功/失败统计 + 哈希
"""

import argparse
import csv
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

VERSION = "1.0"
# 注意：PUG 的 /property/ 端点用 TXT 输出时**只能给一个属性**；
# 同时要 Canonical+Isomeric 会返回 400（实测 5/5 全失败）。这里先取 Canonical，
# 失败再退回 Isomeric。
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/CanonicalSMILES/TXT"
PUBCHEM_ISO = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/IsomericSMILES/TXT"
UA = {"User-Agent": "network-pharmacology-pro/3.0 (stp_prepare)"}


def _read_table(path):
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        f.seek(0)
        rows = list(csv.DictReader(f, delimiter=delim))
        hdrs = list(rows[0].keys()) if rows else []
    return rows, hdrs


def pick_compounds(path, n=5, col=None, rank_by="content"):
    """从成分表选出 top-N。rank_by: content（含量）/ ob（口服利用度）/ order（原顺序）。"""
    rows, hdrs = _read_table(path)
    if not rows:
        return []
    low = {h.lower().strip(): h for h in hdrs}
    ccol = col or low.get("compound") or low.get("name") or hdrs[0]
    rows = [r for r in rows if (r.get(ccol) or "").strip()]
    if rank_by == "ob":
        key = low.get("ob")
        rows.sort(key=lambda r: -_f(r.get(key)) if key else 0)
    elif rank_by == "content":
        key = low.get("content")
        if key:
            rows.sort(key=lambda r: -_f(r.get(key)))
    out, seen = [], set()
    for r in rows:
        name = (r.get(ccol) or "").strip()
        if name.lower() in seen:
            continue
        seen.add(name.lower())
        out.append({"name": name,
                    "mol_id": (r.get(low.get("mol_id", "")) or "").strip(),
                    "content": (r.get(low.get("content", "")) or "").strip()})
        if len(out) >= n:
            break
    return out


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def _get_txt(url, timeout=20):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
        return [x.strip() for x in r.read().decode("utf-8", "ignore").strip().splitlines() if x.strip()]


def resolve_smiles(name, timeout=20, retries=2):
    """PubChem 名称 → (smiles, cid, err)。先 Canonical，失败再 Isomeric；不抛异常。"""
    for tmpl in (PUBCHEM, PUBCHEM_ISO):
        smi, cid, err = _resolve_one(tmpl % urllib.parse.quote(name), timeout, retries)
        if smi:
            return smi, cid, ""
    return "", "", err


def _resolve_one(url, timeout=20, retries=2):
    last = None
    for k in range(retries + 1):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=timeout) as r:
                txt = r.read().decode("utf-8", "ignore").strip().splitlines()
            smiles = txt[0].strip() if txt else ""
            if smiles and " " not in smiles:
                return smiles, "", ""
            # TXT 两列时第二列是 isomeric；取更"具体"的那个
            if txt:
                cand = [t.strip() for t in txt[0].split("\t") if t.strip()]
                if cand:
                    return max(cand, key=len), "", ""
            return "", "", "empty_response"
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return "", "", "not_found_in_pubchem"
            last = "http_%d" % e.code
        except Exception as e:
            last = "%s" % type(e).__name__
        if k < retries:
            time.sleep(1.0 * (k + 1))
    return "", "", last or "unknown"


def resolve_cid(name, timeout=15):
    u = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/cids/TXT"
         % urllib.parse.quote(name))
    try:
        with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=timeout) as r:
            return r.read().decode("utf-8", "ignore").strip().splitlines()[0]
    except Exception:
        return ""


def main(argv=None):
    ap = argparse.ArgumentParser(description="stp_prepare: 成分名 → SMILES（STP 输入）")
    ap.add_argument("--compounds", default=None, help="成分表 CSV/TSV（需含 compound/name 列）")
    ap.add_argument("--names", default=None, help="直接给名字，逗号分隔")
    ap.add_argument("--col", default=None)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--rank-by", default="content", choices=["content", "ob", "order"])
    ap.add_argument("--offline", action="store_true", help="只生成表单，不联网解析 SMILES")
    ap.add_argument("-o", "--outdir", default="out/stp")
    a = ap.parse_args(argv)

    if a.names:
        items = [{"name": x.strip(), "mol_id": "", "content": ""}
                 for x in a.names.split(",") if x.strip()][:a.top]
    elif a.compounds:
        items = pick_compounds(a.compounds, a.top, a.col, a.rank_by)
    else:
        print("[stp] 需要 --compounds 或 --names")
        return 2
    if not items:
        print("[stp] 未选出任何成分")
        return 2

    os.makedirs(a.outdir, exist_ok=True)
    cache_p = os.path.join(a.outdir, "smiles_cache.json")
    cache = {}
    if os.path.exists(cache_p):
        try:
            cache = json.load(open(cache_p, encoding="utf-8"))
        except Exception:
            cache = {}

    ok, fail = [], []
    for it in items:
        nm = it["name"]
        if nm in cache and cache[nm].get("smiles"):
            it.update(cache[nm])
            ok.append(it)
            continue
        if a.offline:
            fail.append({**it, "error": "offline_mode"})
            continue
        smi, _, err = resolve_smiles(nm)
        if smi:
            it["smiles"] = smi
            it["pubchem_cid"] = resolve_cid(nm)
            cache[nm] = {"smiles": smi, "pubchem_cid": it["pubchem_cid"]}
            ok.append(it)
        else:
            fail.append({**it, "error": err})
        time.sleep(0.25)          # 对 PubChem 友好一点

    json.dump(cache, open(cache_p, "w", encoding="utf-8"), ensure_ascii=False, indent=1)

    inp = os.path.join(a.outdir, "stp_input.txt")
    with open(inp, "w", encoding="utf-8") as f:
        for it in ok:
            f.write(it["smiles"] + "\n")
    csvp = os.path.join(a.outdir, "stp_compounds.csv")
    with open(csvp, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "smiles", "pubchem_cid", "mol_id", "content",
                                          "status", "error"], extrasaction="ignore")
        w.writeheader()
        for it in ok:
            w.writerow({**it, "status": "ok"})
        for it in fail:
            w.writerow({**it, "status": "failed"})
    man = {"n_requested": len(items), "n_ok": len(ok), "n_failed": len(fail),
           "failed": [{"name": x["name"], "error": x.get("error")} for x in fail],
           "input_file": inp, "table": csvp,
           "next": ("打开 https://www.swisstargetprediction.ch/ → 把 stp_input.txt 内容粘进 "
                    "SMILES 框 → 选 Homo sapiens → Predict → 右上 Export to CSV"),
           "created_at": time.strftime("%Y-%m-%d %H:%M:%S")}
    json.dump(man, open(os.path.join(a.outdir, "stp_manifest.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)

    print("[stp] %d/%d 个成分解析成功 → %s" % (len(ok), len(items), inp))
    if fail:
        print("  未解析（需手工补 SMILES 或换成分）：")
        for x in fail:
            print("    %-28s %s" % (x["name"], x.get("error")))
    print("  下一步：%s" % man["next"])
    return 0 if ok else 3


if __name__ == "__main__":
    sys.exit(main())
