#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""np_scm0061.py — 把 SCM0061-2021《网络药理学评价方法指南》附录 A 的评价要素做成可自动打分的评分卡。

为什么要有它：SCM0061-2021 定义了"该评什么"（可靠性/规范性/合理性 × 数据收集/网络分析/结果验证），
但它是文档，没有判定。本模块把每条要素落到**产物探针**上：产物在不在、里面有没有关键字段，
给出 0/1/2 分与证据路径。这样"是否符合标准"从主观声明变成可复算的 JSON。

评分口径（写死，避免解释空间）：
  0 = 未检出对应产物 → 该要素未开展/未留证
  1 = 产物存在但缺关键字段，或只覆盖部分环节
  2 = 产物存在且含判定所需字段

用法：
  python scripts/np_scm0061.py --dir np_out/Quercetin_full [--json]
输出：<dir>/scm0061_scorecard.json 与 scm0061_scorecard.md
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

# (编号, 维度, 要素名, 说明)
ITEMS = [
    ("R1", "可靠性", "数据源登记与版本",
     "SCM0061 6.1.1：数据出处、版本、收集时点"),
    ("R2", "可靠性", "快照/溯源",
     "SCM0061 6.3.1：数据可追溯（引用/实验细节）"),
    ("R3", "可靠性", "算法与工具版本记录",
     "SCM0061 6.1.2：软件/算法实现的可靠性与出处"),
    ("R4", "可靠性", "零模型/随机化显著性",
     "SCM0061 6.1.2：分析方法的统计可靠性"),
    ("R5", "可靠性", "阴性对照或基准回归",
     "SCM0061 6.1.2：用标准数据集验证方法"),
    ("N1", "规范性", "产物契约校验",
     "SCM0061 6.2.2：分析路径明确、产物一致"),
    ("N2", "规范性", "中间产物可导出",
     "SCM0061 6.2.1：信息提取与转换规范"),
    ("N3", "规范性", "标准化与未匹配披露",
     "SCM0061 6.2.1：异常信息处理"),
    ("N4", "规范性", "物种/同源映射声明",
     "SCM0061 6.2.1：跨库信息转换"),
    ("N5", "规范性", "全流程 manifest",
     "SCM0061 6.2.2：分析路径可溯源"),
    ("P1", "合理性", "成分筛选尺子明确",
     "SCM0061 6.3.1：筛选原则与条件"),
    ("P2", "合理性", "阈值敏感性",
     "SCM0061 6.3.2：网络分析参数的合理性"),
    ("P3", "合理性", "多重性控制",
     "SCM0061 6.3.2：指标选择与判读"),
    ("P4", "合理性", "富集分层/冗余折叠",
     "SCM0061 6.3.2：生物功能注释的合理性"),
    ("P5", "合理性", "证据等级（L0–L4）",
     "SCM0061 6.1.3：验证方法与模型的可靠性"),
    ("P6", "合理性", "适用域声明",
     "SCM0061 6.3.3：模型适用边界"),
    ("R6", "可靠性", "靶点来源可追溯率",
     "SCM0061 6.1.1/6.3.1：数据来源与可溯源（回顾性分析显示这是复现的首要障碍）"),
]


def _files(root: str, pat: str) -> list[str]:
    return sorted(glob.glob(os.path.join(root, "**", pat), recursive=True))


def _has_token(paths: list[str], token: str, limit: int = 400_000) -> bool:
    rx = re.compile(token, re.I)
    for p in paths:
        try:
            with open(p, encoding="utf-8", errors="ignore") as fh:
                if rx.search(fh.read(limit)):
                    return True
        except OSError:
            continue
    return False


def _nonempty(paths: list[str], min_bytes: int = 2) -> list[str]:
    out = []
    for p in paths:
        try:
            if os.path.getsize(p) >= min_bytes:
                out.append(p)
        except OSError:
            continue
    return out


def score_all(root: str) -> list[dict]:
    rows = []

    def add(iid, score, evidence):
        rows.append({"id": iid, "score": score, "evidence": evidence})

    # R1 数据源登记
    env = _nonempty(_files(root, "environment.json"))
    add("R1", 2 if len(env) >= 3 else (1 if env else 0),
        f"{len(env)} 个 environment.json" if env else "未检出 env 记录")

    # R2 快照/溯源
    prov = _nonempty(_files(root, "pipeline_manifest.json"))
    ok = bool(prov) and _has_token(prov, r'"(sha|hash|prov|source_url)"')
    add("R2", 2 if ok else (1 if prov else 0),
        "pipeline_manifest 含哈希/来源字段" if ok else ("manifest 无溯源字段" if prov else "未检出"))

    # R3 版本记录
    ver = env and _has_token(env, r'"(version|gmx|vina|python|package)"')
    add("R3", 2 if ver else (1 if env else 0),
        "environment.json 含版本字段" if ver else "缺版本字段")

    # R4 零模型
    nul = _nonempty(_files(root, "06_null/*")) or _nonempty(_files(root, "np_null*"))
    add("R4", 2 if nul else 0, f"{len(nul)} 个零模型产物" if nul else "未检出零模型")

    # R5 阴性对照 / 基准
    ref = (_nonempty(_files(root, "decoy*/*")) + _nonempty(_files(root, "bench*/*"))
           + _nonempty(_files(root, "gold*/*")))
    add("R5", 2 if ref else 0, f"{len(ref)} 个对照/基准产物" if ref else "未检出")

    # N1 契约校验
    df = _files(root, "dataflow_check.json")
    ok = bool(df) and _has_token(df, r'"(ok|pass|status)"')
    add("N1", 2 if ok else (1 if df else 0),
        "dataflow_check.json 有状态字段" if ok else "未检出契约校验")

    # N2 中间产物
    mid = _nonempty(_files(root, "targets_standardized.csv")) + _nonempty(_files(root, "ppi_edges_clean.tsv"))
    add("N2", 2 if len(mid) >= 2 else (1 if mid else 0), f"{len(mid)} 个关键中间产物")

    # N3 未匹配披露
    un = _files(root, "targets_unmatched.csv")
    add("N3", 2 if un else 0, f"{len(un)} 个未匹配清单" if un else "未披露未匹配项")

    # N4 同源映射
    orth = _nonempty(_files(root, "11_ortholog/*"))
    add("N4", 2 if orth else 0, f"{len(orth)} 个同源映射产物" if orth else "未检出")

    # N5 manifest
    man = _nonempty(_files(root, "pipeline_manifest.json"))
    ok = bool(man) and _has_token(man, r'"(nodes|stages|commands|calls)"')
    add("N5", 2 if ok else (1 if man else 0), "manifest 含节点/命令" if ok else "未检出")

    # P1 筛选尺子
    comp = _files(root, "compounds_kept.csv")
    ok = bool(comp) and _has_token(comp, r"gate_used|gate_detail")
    add("P1", 2 if ok else (1 if comp else 0),
        "compounds_kept.csv 含 gate_used/gate_detail" if ok else "缺筛选尺子字段")

    # P2 敏感性
    sens = _nonempty(_files(root, "*sens*/*")) + _nonempty(_files(root, "*sensitivity*"))
    add("P2", 2 if sens else 0, f"{len(sens)} 个敏感性产物" if sens else "未检出")

    # P3 多重性
    enr = _nonempty(_files(root, "07_enrich/*"))
    ok = bool(enr) and _has_token(enr, r"q[_ ]?value|fdr|adjusted[_ ]?p")
    add("P3", 2 if ok else (1 if enr else 0),
        "富集产物含 FDR/q 值" if ok else "富集产物无多重性字段")

    # P4 富集分层
    add("P4", 2 if enr else 0, f"{len(enr)} 个富集产物")

    # P5 证据等级
    ev = _nonempty(_files(root, "target_scores.csv"))
    ok = bool(ev) and _has_token(ev, r"\bL[0-4]\b|evidence_level|level")
    add("P5", 2 if ok else (1 if ev else 0),
        "target_scores.csv 含等级列" if ok else "缺证据等级列")

    # P6 适用域
    dom = _nonempty(_files(root, "*domain*/*")) + _nonempty(_files(root, "*applicability*"))
    add("P6", 2 if dom else 0, f"{len(dom)} 个适用域产物" if dom else "未检出")

    # R6 靶点来源可追溯率（回顾性分析发现：靶点来源不一致是复现的首要障碍）
    led = _nonempty(_files(root, "evidence_ledger.csv"))
    if not led:
        add("R6", 0, "未检出证据账本")
    else:
        try:
            import csv as _csv
            rows_ = list(_csv.DictReader(open(led[0], encoding="utf-8-sig")))
            n_all = len(rows_)
            n_src = sum(1 for x in rows_ if (x.get("source") or "").strip())
            n_hash = sum(1 for x in rows_ if (x.get("provenance_hash") or "").strip())
            n_ok = sum(1 for x in rows_
                       if (x.get("source") or "").strip() and (x.get("provenance_hash") or "").strip())
            rate = n_ok / n_all if n_all else 0.0
            add("R6", 2 if rate >= 0.95 else (1 if rate >= 0.5 else 0),
                f"{n_ok}/{n_all} 条同时含 source 与 provenance_hash（{rate:.0%}）；"
                f"distinct_source={len({(x.get('source') or '').strip() for x in rows_})}")
        except Exception as exc:
            add("R6", 1, f"账本读取异常：{type(exc).__name__}")

    # 附上说明
    meta = {i: (d, name, note) for i, d, name, note in ITEMS}
    for r in rows:
        d, name, note = meta[r["id"]]
        r.update({"dimension": d, "item": name, "requirement": note})
    return rows


def _manifest_stages(root: str) -> set[str]:
    """读 pipeline_manifest.json 里真正执行过的 stage id。"""
    p = os.path.join(root, "pipeline_manifest.json")
    if not os.path.exists(p):
        return set()
    try:
        d = json.load(open(p, encoding="utf-8"))
        return {n.get("id") for n in (d.get("nodes") or []) if n.get("id")}
    except Exception:
        return set()


# 每个要素对应的 stage（用于区分"未开展"与"已运行但降级"）；None 表示无对应 stage
ITEM_STAGE = {
    "R1": None, "R2": None, "R3": None,
    "R4": "nullmodel", "R5": None,
    "N1": None, "N2": "standardize", "N3": "standardize", "N4": "ortholog", "N5": None,
    "P1": "compounds", "P2": None, "P3": "enrich", "P4": "enrich", "P5": None, "P6": None,
}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="np_scm0061: SCM0061-2021 评价要素自动评分卡")
    ap.add_argument("--dir", required=True, help="一次 NP 运行的输出目录")
    ap.add_argument("--json", action="store_true", help="只打印 JSON")
    a = ap.parse_args(argv)

    root = a.dir
    if not os.path.isdir(root):
        raise SystemExit(f"[FAIL] 目录不存在：{root}")

    rows = score_all(root)
    ran = _manifest_stages(root)
    for r in rows:
        st = ITEM_STAGE.get(r["id"])
        r["stage"] = st
        r["stage_ran"] = (st in ran) if st else None
        if r["score"] == 0 and st and r["stage_ran"]:
            r["evidence"] = f"{r['evidence']}（stage 已运行但无产物：软失败/降级）"
        elif r["score"] == 0 and st and not r["stage_ran"]:
            r["evidence"] = f"{r['evidence']}（该 stage 未在本案例运行）"

    dims = {}
    for r in rows:
        d = dims.setdefault(r["dimension"], {"got": 0, "max": 0, "items": []})
        d["got"] += r["score"]
        d["max"] += 2
        d["items"].append(r)

    total_got = sum(d["got"] for d in dims.values())
    total_max = sum(d["max"] for d in dims.values())
    pct = 100.0 * total_got / total_max if total_max else 0.0
    # 只对"本案例本应开展"的要素算合规率（stage 未运行的排除，避免把"没做"算成"做错"）
    exe_max = sum(2 for r in rows if r["stage_ran"] in (True, None))
    exe_pct = 100.0 * total_got / exe_max if exe_max else 0.0
    grade = "A" if pct >= 85 else ("B" if pct >= 65 else ("C" if pct >= 40 else "D"))

    out = {"dir": os.path.abspath(root), "standard": "SCM0061-2021",
           "total_score": total_got, "total_max": total_max,
           "percent": round(pct, 1),
           "percent_executed_scope": round(exe_pct, 1),
           "executed_max": exe_max,
           "stages_ran": sorted(ran),
           "grade": grade,
           "dimensions": {k: {"got": v["got"], "max": v["max"]} for k, v in dims.items()},
           "items": rows,
           "note": ("0=未检出该要素产物；1=产物存在但缺关键字段；2=产物与字段齐备。"
                    "percent 为全要素口径；percent_executed_scope 只统计本案例实际开展的 stage，"
                    "用于区分「未开展」与「做而不足」。")}

    jp = os.path.join(root, "scm0061_scorecard.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)

    lines = [f"# SCM0061-2021 评分卡\n",
             f"- 目录：`{out['dir']}`",
             f"- 总分：**{total_got}/{total_max}（{pct:.1f}%，等级 {grade}）**\n",
             "| 维度 | 得分 | 满分 |", "|---|---|---|"]
    for k, v in dims.items():
        lines.append(f"| {k} | {v['got']} | {v['max']} |")
    lines += ["", "| 编号 | 维度 | 要素 | 分 | 证据 |", "|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['id']} | {r['dimension']} | {r['item']} | {r['score']} | {r['evidence']} |")
    lines += ["", "0=未检出该要素产物；1=产物存在但缺关键字段；2=齐备。"]
    open(os.path.join(root, "scm0061_scorecard.md"), "w", encoding="utf-8").write("\n".join(lines))

    if a.json:
        print(json.dumps(out, ensure_ascii=False, indent=1))
    else:
        print(f"SCM0061 评分：{total_got}/{total_max} = {pct:.1f}%  等级 {grade}  已开展口径 {exe_pct:.1f}%（{total_got}/{exe_max}）")
        for k, v in dims.items():
            print(f"  {k}: {v['got']}/{v['max']}")
        miss = [r["id"] for r in rows if r["score"] == 0]
        if miss:
            print(f"  未检出：{', '.join(miss)}")
        print(f"→ {jp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
