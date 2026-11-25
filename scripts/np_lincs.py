#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_lincs.py —— 反向药物重定位（LINCS/CMap）：模拟"逆转效应"的老药新用候选

科学逻辑
    若中药核心成分引起的表达特征（signature）与某疾病特征**负相关**，
    则该成分（或能模拟该特征的已知药物）可能逆转疾病 → Drug Repurposing 候选。
    CMap 的 connectivity score ∈[-1,1]：负分=逆转（候选），正分=模拟疾病（有害）。

数据现实（实测 2026-11-23）
    clue.io 官网公告：**2026-01-31 起全站退役**（含 Tools→Query，资金限制），
    公开数据迁至 GEO。旧 api.clue.io/api/v2.0/query 端点与 template 已失效（保留仅作历史）。
    → 真实 connectivity 查询改用 **L1000FWD 免费 API**（Maayan Lab，无需 key）：
      POST /sig_search → GET /result/topn/{id}，mode opposite=逆转药物（经典 CMap）。
    核心算法（签名构造/连通性打分/候选排序）仍保留，纯 Python 可测，不伪造结果。

流程
    signature       核心靶点 + 方向表 → up/down 基因集
    score           本地 L1000 子集 → CMap 风格 connectivity score
    query-l1000fwd  ★真实查询：up/down 基因 → L1000FWD connectivity（免费，替代 clue.io）
    template        旧 clue.io API 模板（已退役，仅历史参考）

用法
    python np_lincs.py signature --targets core.txt --direction dir.csv -o out/
    python np_lincs.py score --signature sig.json --l1000 profiles.csv -o out/
    python np_lincs.py query-l1000fwd --up up.txt --down down.txt -o out/    # 真实 CMap 替代
    python np_lincs.py selftest
"""
import argparse
import csv
import json
import math
import os
import sys
import time
import urllib.request

L1000FWD_BASE = "https://maayanlab.cloud/L1000FWD"

HERE = os.path.dirname(os.path.abspath(__file__))


def read_direction(path):
    """靶点方向表：gene,direction(up/down)。"""
    d = {}
    with open(path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            g = (r.get("gene") or r.get("target") or "").strip().upper()
            v = (r.get("direction") or r.get("change") or "").strip().lower()
            if g and v in ("up", "down", "上调", "下调"):
                d[g] = "up" if v in ("up", "上调") else "down"
    return d


def build_signature(targets, direction):
    """核心靶点 + 方向 → up/down 基因集 signature。"""
    up = sorted(g for g in targets if direction.get(g) == "up")
    down = sorted(g for g in targets if direction.get(g) == "down")
    return {"up_genes": up, "down_genes": down,
            "n_annotated": len(up) + len(down),
            "n_unannotated": len(targets) - len(up) - len(down),
            "note": "未标注方向的靶点不纳入签名；可用 tf_targets 推断后再补"}


def connectivity(sig, profile):
    """CMap 风格连通性：比较 signature(up/down) 与药物扰动 profile(up/down) 的一致性。
    score = (|sig_up∩prof_up| + |sig_down∩prof_down| - |sig_up∩prof_down| - |sig_down∩prof_up|) / N
    负值 = 逆转效应（药物把疾病signature反过来）→ 候选。
    """
    su, sd = set(sig["up_genes"]), set(sig["down_genes"])
    pu, pd = set(profile.get("up_genes", [])), set(profile.get("down_genes", []))
    N = max(1, len(su) + len(sd))
    score = (len(su & pu) + len(sd & pd) - len(su & pd) - len(sd & pu)) / N
    return round(score, 4)


def score_profiles(sig, profiles, top=20):
    scored = [{"drug": p.get("drug", p.get("name", "")),
               "connectivity": connectivity(sig, p),
               "reversal": bool(connectivity(sig, p) < -0.2)} for p in profiles]
    scored.sort(key=lambda x: x["connectivity"])  # 最负（最强逆转）在前
    return {"n_profiles": len(scored), "n_reversal": sum(1 for x in scored if x["reversal"]),
            "top_reversal": scored[:top],
            "note": "connectivity<−0.2 记为逆转候选；本地小型库仅作方法验证，"
                    "正式候选需 LINCS L1000 全库 + CMap 官方 connectivity score 确认"}


def gen_template(outdir):
    """生成 CMap/LINCS 官方查询脚本模板（用户填 key 后可跑真实比对）。"""
    os.makedirs(outdir, exist_ok=True)
    script = '''#!/usr/bin/env python
# CMap/LINCS 查询模板 —— 需用户填写 API key，本机直接运行不伪造结果
# CMap: https://clue.io/api  (需登录获取 key)
import json, urllib.request, sys
API_KEY = "<YOUR_CMAP_API_KEY>"
SIGNATURE_UP = ["IL6","TNF","PTGS2"]      # 由 np_lincs signature 产出填入
SIGNATURE_DOWN = ["TP53","CASP3"]
url = "https://api.clue.io/api/v2.0/query"
payload = {"api_key": API_KEY, "q": {"up_genes": SIGNATURE_UP, "dn_genes": SIGNATURE_DOWN},
           "filter": {"pert_type": "trt_cp"}}
req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                             headers={"Content-Type":"application/json"})
try:
    resp = json.loads(urllib.request.urlopen(req, timeout=120).read())
    json.dump(resp, open("cmap_result.json","w"), ensure_ascii=False, indent=1)
    print("CMap 返回", len(resp.get("result", [])), "条；负 connectivity score 为逆转候选")
except Exception as e:
    print("查询失败（请检查 key/网络）：", type(e).__name__, str(e)[:120])
'''
    p = os.path.join(outdir, "query_cmap.py")
    open(p, "w", encoding="utf-8").write(script)
    return p


def _urlopen_json(url, data=None, timeout=60, retries=3):
    """带重试的 JSON GET/POST（L1000FWD 免费 API，无需 key）。"""
    import ssl
    ctx = ssl.create_default_context()
    last = None
    for _ in range(retries):
        try:
            if data is None:
                req = urllib.request.Request(url)
            else:
                req = urllib.request.Request(
                    url, data=data, headers={"Content-Type": "application/json"})
            return json.loads(urllib.request.urlopen(req, timeout=timeout, context=ctx)
                              .read().decode("utf-8", "ignore"))
        except Exception as e:
            last = e
            time.sleep(2)
    raise last


def query_l1000fwd(up_file, down_file, outdir, topn=50, enrich=True):
    """用 L1000FWD 免费 API 做真实 connectivity 查询（clue.io 2026-01-31 退役后的替代）。

    clue.io 于 2026-01-31 因资金问题全站退役（Tools→Query 已下线），
    公开 LINCS L1000 数据迁至 GEO。L1000FWD（Maayan Lab）提供免费无 key 的
    signature 相似性检索 REST API，mode=opposite 即经典 CMap 逆转药物查询。
    API: POST /sig_search {up_genes, down_genes} -> {result_id}
         GET  /result/topn/{result_id} -> {similar:[...], opposite:[...]}
    每条 hit 含 scores(connectivity,负=逆转)/pvals/qvals/combined_scores/sig_id。
    """
    os.makedirs(outdir, exist_ok=True)
    up = [l.strip().upper() for l in open(up_file, encoding="utf-8-sig") if l.strip()]
    down = [l.strip().upper() for l in open(down_file, encoding="utf-8-sig") if l.strip()]
    if not up and not down:
        print("[lincs] 错误：up/down 基因均为空")
        return 2
    print("[lincs] L1000FWD 查询：up=%d down=%d（免费 API，无需 key）" % (len(up), len(down)))
    posted = _urlopen_json(
        L1000FWD_BASE + "/sig_search",
        data=json.dumps({"up_genes": up, "down_genes": down}).encode(), timeout=90)
    rid = posted.get("result_id") or posted.get("id")
    if not rid:
        print("[lincs] 未返回 result_id：%s" % str(posted)[:200])
        return 2
    result = None
    for _ in range(8):
        try:
            result = _urlopen_json(
                "%s/result/topn/%s" % (L1000FWD_BASE, rid), timeout=90, retries=1)
            break
        except Exception:
            time.sleep(8)
    if not result or "opposite" not in result:
        print("[lincs] 结果拉取失败（result_id=%s）" % rid)
        return 2
    json.dump({"result_id": rid, "up": up, "down": down, "raw": result},
              open(os.path.join(outdir, "l1000fwd_result_raw.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)

    def _enrich(hits, mode):
        rows = []
        for i, h in enumerate(hits[:topn], 1):
            sid = h.get("sig_id", "")
            meta = {}
            if enrich:
                try:
                    meta = _urlopen_json(L1000FWD_BASE + "/sig/" + sid, timeout=45, retries=2)
                except Exception:
                    pass
            rows.append({
                "rank": i, "mode": mode,
                "drug": meta.get("pert_desc", ""),
                "pert_id": meta.get("pert_id", ""),
                "cell": meta.get("cell_id", ""),
                "dose": meta.get("pert_dose", ""),
                "time_h": meta.get("pert_time", ""),
                "connectivity_score": h.get("scores"),
                "pval": h.get("pvals"), "qval": h.get("qvals"),
                "combined_score": h.get("combined_scores"),
                "sig_id": sid,
                "n_up": len(meta.get("up_genes", [])),
                "n_down": len(meta.get("down_genes", [])),
            })
            if enrich:
                time.sleep(0.3)
        return rows

    opp = _enrich(result.get("opposite", []), "opposite")
    sim = _enrich(result.get("similar", []), "similar")
    rows = opp + sim
    if rows:
        with open(os.path.join(outdir, "l1000fwd_hits.csv"), "w",
                  encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    json.dump(rows, open(os.path.join(outdir, "l1000fwd_hits.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    n_sig = sum(1 for r in opp if (r["qval"] or 1) < 0.05 and (r["connectivity_score"] or 0) < 0)
    print("[lincs] 逆转候选 %d/%d（opposite, q<0.05）→ l1000fwd_hits.csv" % (n_sig, len(opp)))
    print("[lincs] 模拟候选 %d（similar）；负 connectivity=逆转，正=模拟疾病" % len(sim))
    return 0


def selftest():
    sig = {"up_genes": ["IL6", "TNF"], "down_genes": ["TP53"]}
    # 药物 A：下调 IL6/TNF、上调 TP53 → 完全逆转 → 负分
    rev = {"drug": "A", "up_genes": ["TP53"], "down_genes": ["IL6", "TNF"]}
    # 药物 B：模拟疾病 → 正分
    sim = {"drug": "B", "up_genes": ["IL6", "TNF"], "down_genes": ["TP53"]}
    r = score_profiles(sig, [rev, sim])
    ok = (r["top_reversal"][0]["drug"] == "A"
          and r["top_reversal"][0]["connectivity"] < 0
          and r["n_reversal"] == 1)
    print("selftest: %s (reversal=%s score=%.2f; mimic=%.2f)" %
          ("PASS" if ok else "FAIL", r["top_reversal"][0]["drug"],
           r["top_reversal"][0]["connectivity"], connectivity(sig, sim)))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_lincs: 反向药物重定位框架")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("signature")
    p.add_argument("--targets", required=True); p.add_argument("--direction", required=True)
    p.add_argument("-o", "--outdir", default="np_out/lincs"); p.set_defaults(fn="sig")
    p = sub.add_parser("score")
    p.add_argument("--signature", required=True); p.add_argument("--l1000", required=True)
    p.add_argument("-o", "--outdir", default="np_out/lincs"); p.set_defaults(fn="score")
    p = sub.add_parser("template"); p.add_argument("-o", "--outdir", default="np_out/lincs")
    p.set_defaults(fn="template")
    p = sub.add_parser("query-l1000fwd",
                       help="L1000FWD 免费 API 真实 connectivity 查询（clue.io 退役替代）")
    p.add_argument("--up", required=True, help="上调基因文件（每行一个）")
    p.add_argument("--down", required=True, help="下调基因文件（每行一个）")
    p.add_argument("--topn", type=int, default=50)
    p.add_argument("--no-enrich", action="store_true", help="不拉药物名元数据（更快）")
    p.add_argument("-o", "--outdir", default="np_out/lincs"); p.set_defaults(fn="query")
    p = sub.add_parser("selftest"); p.set_defaults(fn="selftest")
    a = ap.parse_args(argv)
    if a.fn == "selftest":
        return selftest()
    os.makedirs(a.outdir, exist_ok=True)
    if a.fn == "sig":
        targets = {l.strip().upper() for l in open(a.targets, encoding="utf-8-sig") if l.strip()}
        sig = build_signature(targets, read_direction(a.direction))
        json.dump(sig, open(os.path.join(a.outdir, "signature.json"), "w", encoding="utf-8"),
                  ensure_ascii=False, indent=1)
        print("[lincs] 签名 up=%d down=%d 未标注=%d → signature.json" %
              (len(sig["up_genes"]), len(sig["down_genes"]), sig["n_unannotated"]))
        return 0
    if a.fn == "query":
        return query_l1000fwd(a.up, a.down, a.outdir, topn=a.topn, enrich=not a.no_enrich)
    if a.fn == "template":
        p = gen_template(a.outdir)
        print("[lincs] CMap 查询模板 → %s（需填 API key，不伪造结果）" % p)
        return 0
    sig = json.load(open(a.signature, encoding="utf-8"))
    profiles = []
    with open(a.l1000, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            profiles.append({"drug": r.get("drug", ""),
                             "up_genes": [x for x in (r.get("up_genes") or "").split(";") if x],
                             "down_genes": [x for x in (r.get("down_genes") or "").split(";") if x]})
    res = score_profiles(sig, profiles)
    json.dump(res, open(os.path.join(a.outdir, "repurposing_candidates.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[lincs] %d 药物谱，%d 个逆转候选 → repurposing_candidates.json" %
          (res["n_profiles"], res["n_reversal"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
