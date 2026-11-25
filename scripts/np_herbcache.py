#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_herbcache.py —— 修复 fallback.json L1 死链：动态拉取 ETCM2/HERB 成分-靶点并写入本地缓存

背景（真实审计发现）
    policies/fallback.json 的 L1 partial 声明 "TCMSP 失败→本地 72 味缓存 / HERB / ETCM"，
    但 data/herb_cache/*.json 目录实测不存在 —— L1 降级分支不可达，属于死链。

诚实边界（延续 500 篇复现的核心教训：禁止静默返回空/编造）
    * ETCM2  : 已实测真实 API（www.tcmip.cn:18124，POST /home/ingredient_target/），可程序化拉取。
    * HERB   : 站点可达但查询为 SPA 私有接口、批量下载返回 "unavailable path"，
               实测不可程序化获取 —— 本模块仅做"探针 + 如实回报不可用"，绝不伪造数据。
    任何失败都写 status=unavailable + reason，而非空 JSON；下游据此显式降级到 L2 manual。

用法
    python np_herbcache.py fetch --herb ChuShiZi --source etcm        # 拉单味写入缓存
    python np_herbcache.py fetch --herbs 龟甲,楮实子,何首乌 --source etcm
    python np_herbcache.py probe                                       # 探测各源可用性并落盘状态
    python np_herbcache.py status                                      # 查看已缓存药材

输出
    data/herb_cache/<herb>.json   单味缓存（含 provenance + 时间戳 + sha256）
    data/herb_cache/_index.json   缓存索引（药材/记录数/来源/状态）
"""
import argparse
import hashlib
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(HERE, "..", "data", "herb_cache")

# ETCM2 实测 API（2026-11-21 验证可用）
ETCM_BASE = "http://www.tcmip.cn:18124"
ETCM_ENDPOINT = "/home/ingredient_target/"
HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)",
    "Origin": "http://www.tcmip.cn",
    "Referer": "http://www.tcmip.cn/ETCM2/front/",
    "Content-Type": "application/json;charset=UTF-8",
}
# 常见药材中文名 → ETCM2 id（音译驼峰，已在上一轮实测命中）
HERB_ID = {
    "楮实子": "ChuShiZi", "何首乌": "HeShouWu", "龟甲": "GuiJia",
    "熟地黄": "ShuDiHuang", "山茱萸": "ShanZhuYu", "山萸肉": "ShanZhuYu",
    "桃仁": "TaoRen", "红花": "HongHua", "丹参": "DanShen", "川芎": "ChuanXiong",
    "石菖蒲": "ShiChangPu", "郁金": "YuJin", "赤芍": "ChiShao", "黄精": "HuangJing",
    "枸杞子": "GouQiZi", "益智仁": "YiZhiRen",
}


def _now():
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def ensure_cache():
    os.makedirs(CACHE_DIR, exist_ok=True)


def _http_json(url, payload=None, timeout=45):
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    req = urllib.request.Request(url, data=data, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8", "ignore"))


def fetch_etcm(herb_query):
    """从 ETCM2 拉单味成分-靶点。返回统一结构；失败返回 status=unavailable（不返回空列表）。"""
    herb_id = HERB_ID.get(herb_query, herb_query)
    out = {
        "herb": herb_query, "herb_id": herb_id, "source": "ETCM2",
        "url": ETCM_BASE + ETCM_ENDPOINT, "access_date": _now(),
        "status": "unavailable", "reason": "", "records": [], "n_records": 0,
    }
    try:
        payload = {"type": "herb", "id": herb_id, "pageNo": 1, "pageSize": 1000,
                   "orderBy": "", "sequential": ""}
        j = _http_json(ETCM_BASE + ETCM_ENDPOINT, payload)
        data = j.get("data") or []
        if not data:
            out["reason"] = ("ETCM2 返回空 data（该药材可能为动物药/矿物药，无小分子成分-靶点表；"
                             "或 id 未命中）。按诚信原则不编造。")
            return out
        def scalar(v):
            """ETCM2 字段值有时是 list/dict，取第一个字符串。"""
            if isinstance(v, list):
                v = v[0] if v else ""
            if isinstance(v, dict):
                v = v.get("value") or v.get("name") or ""
            return str(v).strip() if v is not None else ""

        recs = []
        for r in data:
            # 实测字段（2026-11-21）：Ingredient Name / Gene Symbol / Target Name /
            #   UniProt Accession Number / References（extract_etcm_full.py 已验证）
            comp = scalar(r.get("Ingredient Name") or r.get("compound") or r.get("compoundName"))
            gene = scalar(r.get("Gene Symbol") or r.get("gene") or r.get("geneSymbol"))
            pmid = scalar(r.get("References") or r.get("pubmed") or r.get("pmid") or r.get("reference"))
            tgt = scalar(r.get("Target Name") or r.get("Target"))
            uniprot = scalar(r.get("UniProt Accession Number") or r.get("uniprot"))
            if comp and gene:
                recs.append({"compound": comp, "gene": gene.upper(), "target_name": tgt,
                             "uniprot": uniprot, "evidence": pmid[:120]})
        out["records"] = recs
        out["n_records"] = len(recs)
        out["status"] = "ok" if recs else "unavailable"
        if not recs:
            out["reason"] = "ETCM2 返回记录但无法解析出 compound/gene 字段，未臆造。"
    except (urllib.error.URLError, urllib.error.HTTPError, ValueError, KeyError) as e:
        out["reason"] = "ETCM2 请求失败：%s: %s" % (type(e).__name__, str(e)[:120])
    return out


def probe_herb():
    """HERB 探针：实测其 SPA 接口与下载端点是否可程序化获取。诚实回报，不绕过验证。"""
    checks = []
    for url in ["http://herb.ac.cn/", "http://herb.ac.cn/download/file/?file_path=HERB_herb_info.txt"]:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": HEADERS["User-Agent"]})
            with urllib.request.urlopen(req, timeout=15) as r:
                body = r.read(200).decode("utf-8", "ignore")
            ok = "unavailable path" not in body.lower() and "安全验证" not in body
            checks.append({"url": url, "http": r.status, "usable": ok,
                           "detail": "" if ok else "下载端点禁用/需交互验证"})
        except Exception as e:
            checks.append({"url": url, "http": None, "usable": False,
                           "detail": "%s: %s" % (type(e).__name__, str(e)[:80])})
    usable = any(c["usable"] and "download" in c["url"] for c in checks)
    return {"source": "HERB", "access_date": _now(), "status": "ok" if usable else "unavailable",
            "reason": "" if usable else "站点可达但查询为 SPA 私有接口、批量下载被禁用；"
                                       "需人工浏览器导出后手动导入（npm_ import 接口）。",
            "checks": checks}


def write_cache(obj):
    ensure_cache()
    herb = obj["herb"]
    path = os.path.join(CACHE_DIR, "%s.json" % herb)
    raw = json.dumps(obj, ensure_ascii=False, indent=1, sort_keys=True)
    obj["sha256_16"] = hashlib.sha256(raw.encode()).hexdigest()[:16]
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1)
    _update_index(herb, obj, path)
    return path


def _update_index(herb, obj, path):
    idx_path = os.path.join(CACHE_DIR, "_index.json")
    index = {}
    if os.path.exists(idx_path):
        try:
            index = json.load(open(idx_path, encoding="utf-8"))
        except Exception:
            index = {}
    index[herb] = {
        "source": obj["source"], "status": obj["status"], "n_records": obj["n_records"],
        "access_date": obj["access_date"], "sha256_16": obj.get("sha256_16", ""),
        "file": os.path.basename(path), "reason": obj.get("reason", ""),
    }
    with open(idx_path, "w", encoding="utf-8") as f:
        json.dump(index, f, ensure_ascii=False, indent=1, sort_keys=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description="np_herbcache: 动态填充 fallback L1 本地缓存")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("fetch", help="拉取药材成分-靶点写入缓存")
    p.add_argument("--herb", default=None, help="单味药材名（中文或 ETCM id）")
    p.add_argument("--herbs", default=None, help="逗号分隔多味")
    p.add_argument("--source", default="etcm", choices=["etcm", "herb"])

    sub.add_parser("probe", help="探测各数据源可用性")
    sub.add_parser("status", help="显示缓存索引")
    a = ap.parse_args(argv)

    if a.cmd == "probe":
        ensure_cache()
        result = {"ETCM2": fetch_etcm("楮实子"), "HERB": probe_herb()}
        # ETCM probe 只需状态，不写缓存；清掉可能的大 records
        result["ETCM2"] = {k: v for k, v in result["ETCM2"].items() if k != "records"}
        p = os.path.join(CACHE_DIR, "_probe.json")
        with open(p, "w", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=False, indent=1)
        for k, v in result.items():
            print("[probe] %-6s %-11s %s" % (k, v["status"], v.get("reason", "")[:70]))
        print("[probe] 状态已写入", p)
        return 0

    if a.cmd == "status":
        idx = os.path.join(CACHE_DIR, "_index.json")
        if not os.path.exists(idx):
            print("[status] 缓存为空（运行 fetch 子命令填充）")
            return 0
        index = json.load(open(idx, encoding="utf-8"))
        print("[status] 缓存 %d 味：" % len(index))
        for h, v in sorted(index.items()):
            print("  %-8s %-11s %4d 条  %s" % (h, v["status"], v["n_records"], v["source"]))
        return 0

    # fetch
    herbs = []
    if a.herb:
        herbs = [a.herb]
    elif a.herbs:
        herbs = [x.strip() for x in a.herbs.split(",") if x.strip()]
    else:
        print("[FAIL] 需 --herb 或 --herbs", file=sys.stderr)
        return 2
    for h in herbs:
        obj = fetch_etcm(h) if a.source == "etcm" else probe_herb()
        if a.source == "herb":
            obj["herb"] = h
        path = write_cache(obj)
        print("[fetch] %-8s %-11s %4d 条 → %s" % (h, obj["status"], obj["n_records"], path))
        if obj["status"] != "ok":
            print("        原因：%s" % obj["reason"][:100])
    return 0


if __name__ == "__main__":
    sys.exit(main())
