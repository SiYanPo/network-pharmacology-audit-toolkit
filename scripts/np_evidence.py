#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_evidence.py —— 证据锚点（Evidence Anchoring）：把 LLM 生成的每句话钉死在本地数据上

问题（LLM 幻觉的物理隔离）
    np_report.py 的 Discussion 若由 LLM 生成，可能编造文献、通路、数字。
    本模块两道防线：
      ① 约束 Prompt（ANCHOR_PROMPT）—— 生成阶段强制每句附引用，无支撑强制输出规定句式。
      ② 正则校验（validate）—— 生成后逐句扫描，任何"裸断言句"判 FAIL 并定位，
         引用指向的文件/PMID 必须真实存在（本地 CSV 存在性 + PMID 格式 + 可选联网核对）。

锚点语法（每句至少一个）
    [Data: PPI_拓扑.csv]              本地数据文件（存在性校验）
    [Data: 表2_核心靶点.csv#IL6]       文件 + 行定位
    [PMID: 12345678]                   PubMed ID（8 位以内数字）
    [DB: DisGeNET:C0030567]            数据库记录
    [Calc: seed=42,score>=0.4]         计算参数（确定性声明）
    无支撑句 → 必须为："数据不足，无法解释。"（EXACT_ALLOWED 白名单）

用法
    python np_evidence.py prompt                          # 打印约束 Prompt（贴给 LLM）
    python np_evidence.py validate --text discussion.txt --datadir np_out/ [--check-pmid]
    python np_evidence.py selftest
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# ① 约束 Prompt —— 生成阶段注入
ANCHOR_PROMPT = """你是网络药理学论文的 Discussion 撰写助手。严格遵守以下证据锚点规则，违反即视为错误输出：

1. 每一个陈述句末尾必须附带至少一个证据锚点，格式为以下之一：
   [Data: 文件名.csv]  或  [Data: 文件名.csv#关键行]
   [PMID: 数字]        或  [DB: 数据库名:编号]
   [Calc: 参数=值]
2. 引用必须来自我提供的本地数据清单；禁止编造文件名、PMID、通路名、数字。
3. 若某句话在本地数据中找不到支撑，禁止猜测，必须原样输出：数据不足，无法解释。
4. 禁止使用无锚点的定性套话，例如"显著改善""发挥关键作用""具有良好前景"——
   除非其后紧跟 [Data: ...] 且该数据确实包含相应数值。
5. 数字必须与锚点文件一致；同一数字在不同句子中必须一致。
6. 预测性结论一律用"提示/可能/推测"，并附 [Calc: ...] 或 [Data: ...]；
   只有实验证据（L3/L4）才用"表明/证明"，并附 [PMID: ...]。

输出示例（合规）：
  "PPI 拓扑分析显示 IL6 度值为 85，位于网络核心 [Data: PPI_拓扑.csv#IL6]。
   既往研究报道 IL6 参与 PD 神经炎症 [PMID: 30240611]。
   本方的槲皮素可能通过调控该靶点发挥作用 [Calc: seed=42;docking=-9.62]。"

输出示例（无支撑时，唯一合法写法）：
  "数据不足，无法解释。"
"""

# ② 校验：Token 化方括号锚点（比严格正则更宽容）
#    支持变形表达：[Data: ./output/PPI.csv 的拓扑分析] / [Data:PPI.csv（图1）] / [ PMID : 123 ]
ANCHOR_TOKEN_RE = re.compile(r"\[\s*(Data|PMID|DB|Calc)\s*:\s*([^\]]+)\]", re.I)
# 保留严格正则供"是否带锚点"的快速判定（宽松 token 也能命中）
ANCHOR_RE = ANCHOR_TOKEN_RE
PMID_RE = re.compile(r"\[\s*PMID\s*:\s*(\d{1,8})[^\]]*\]", re.I)
DATA_RE = re.compile(r"\[\s*Data\s*:\s*([^\]]+)\]", re.I)


def parse_anchor(kind, body):
    """Token 解析锚点内容，容错变形表达。返回 (kind, normalized_value, raw)。
    解析异常时返回 kind='UNPARSED' + 原文，交给兜底逻辑（不误杀，标 WARN 而非 FAIL）。"""
    raw = body.strip()
    try:
        k = kind.strip().lower()
        if k == "pmid":
            m = re.search(r"(\d{1,8})", raw)
            if m:
                return "PMID", m.group(1), raw
            return "UNPARSED", raw, raw
        if k == "data":
            # 取第一个"像文件名"的 token：含 .csv/.tsv/.json/.txt/.md/.svg 或路径分隔符
            m = re.search(r"([\w./\\\-]+\.(?:csv|tsv|json|txt|md|svg|xlsx|pdf))(?:$|[\s,，;；#（(])", raw, re.I)
            if m:
                fname = m.group(1).replace("\\", "/")
                return "Data", fname, raw
            # 无扩展名时取首个空白前 token（去括号）
            token = re.split(r"[\s,，;；#（()）]", raw)[0].strip("'\"")
            return ("Data", token, raw) if token else ("UNPARSED", raw, raw)
        if k == "db":
            parts = re.split(r"[:：]", raw, maxsplit=1)
            return "DB", (parts[0].strip() + (":" + parts[1].strip() if len(parts) > 1 else "")), raw
        if k == "calc":
            return "Calc", raw, raw
    except Exception:
        return "UNPARSED", raw, raw
    return "UNPARSED", raw, raw
EXACT_ALLOWED = "数据不足，无法解释"
# 句末标点（中英文）
SENT_END = re.compile(r"[。！？.!?]+")


def split_sentences(text):
    """按中文句号/感叹/问号切分（保留标点）。
    关键：先把 [Data:...] 等锚点内的英文点号/井号保护为占位符，避免 .csv / #IL6 被误切；
    英文句点仅在『. + 空格 + 大写/中文』处切分，避免缩写误伤。"""
    holders = {}

    def _hold(m):
        key = "\x00A%d\x00" % len(holders)
        holders[key] = m.group(0)
        return key

    protected = ANCHOR_RE.sub(_hold, text)
    # 中文句末标点直接切；英文句点仅在后接空白+非小写时切
    parts = re.split(r"([。！？]|(?<=\w)\.(?=\s+[A-Z\u4e00-\u9fff]))", protected)
    out, cur = [], ""
    for seg in parts:
        cur += seg
        if seg and re.fullmatch(r"[。！？]|\.(?=\s)", seg):
            s = cur.strip()
            if s:
                for k, v in holders.items():
                    s = s.replace(k, v)
                out.append(s)
            cur = ""
    if cur.strip():
        s = cur.strip()
        for k, v in holders.items():
            s = s.replace(k, v)
        out.append(s)
    return out


def validate(text, datadir=None, check_pmid=False, allowed_files=None):
    """逐句校验。返回 (report_dict, bad_sentences)。"""
    sents = split_sentences(text)
    bad = []
    data_refs, pmid_refs = set(), set()
    for idx, s in enumerate(sents, 1):
        anchors = ANCHOR_RE.findall(s)
        # 白名单：整句是"数据不足，无法解释。"则允许无锚点
        if EXACT_ALLOWED in s and not anchors:
            continue
        # 无锚点 → FAIL
        if not anchors:
            bad.append({"idx": idx, "sentence": s, "reason": "无证据锚点",
                        "fix": "补 [Data: 文件] / [PMID: id]，或改为『数据不足，无法解释。』"})
            continue
        # Token 解析每个锚点：Data 引用做文件存在性校验；解析异常降级 WARN 不误杀
        for kind, body in ANCHOR_TOKEN_RE.findall(s):
            pk, val, raw = parse_anchor(kind, body)
            if pk == "UNPARSED":
                bad.append({"idx": idx, "sentence": s, "severity": "WARN",
                            "reason": "锚点解析异常（WARN 不误杀）：%s" % raw[:50],
                            "fix": "建议规范为 [Data: 文件名.csv]，当前按有锚点处理"})
                continue
            if pk == "Data":
                data_refs.add(val)
                if datadir and val:
                    # 支持 ./output/xxx.csv 变形：取 basename 再在 datadir 下找
                    base = os.path.basename(val)
                    cands = [os.path.join(datadir, val), os.path.join(datadir, base),
                             os.path.join(datadir, "figs", base)]
                    if not any(os.path.exists(c) for c in cands):
                        if allowed_files and (val in allowed_files or base in allowed_files):
                            continue
                        bad.append({"idx": idx, "sentence": s,
                                    "reason": "引用文件不存在：%s（解析自 %s）" % (val, raw[:40]),
                                    "fix": "检查文件名，或在 --datadir 下提供该文件"})
        for p in PMID_RE.findall(s):
            pmid_refs.add(p)
            if check_pmid:
                ok = _pmid_exists(p)
                if not ok:
                    bad.append({"idx": idx, "sentence": s, "reason": "PMID 无法核实：%s" % p,
                                "fix": "确认 PMID 真实存在，否则删除该引用"})
    n_fail = sum(1 for b in bad if b.get("severity", "FAIL") == "FAIL")
    n_warn = sum(1 for b in bad if b.get("severity") == "WARN")
    report = {
        "sentences": len(sents),
        "anchored": sum(1 for s in sents if ANCHOR_RE.search(s)),
        "data_refs": sorted(data_refs),
        "pmid_refs": sorted(pmid_refs),
        "bad_count": n_fail,
        "warn_count": n_warn,
        "verdict": "PASS" if n_fail == 0 else "FAIL",
    }
    return report, bad


def _pmid_exists(pmid):
    """可选联网核对 PMID 是否存在（NCBI esummary）。失败返回 False 但不崩溃。"""
    import urllib.request
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=%s&retmode=json" % pmid
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "np-pro-evidence/1.0"})
        import json as _json
        with urllib.request.urlopen(req, timeout=12) as r:
            j = _json.loads(r.read().decode("utf-8", "ignore"))
        return pmid in (j.get("result") or {}).get("uids", [])
    except Exception:
        return False


def selftest():
    good = (
        "PPI 拓扑分析显示 IL6 为网络核心节点 [Data: PPI_拓扑.csv#IL6]。"
        "另有研究报道相关机制 [PMID: 30240611]。"
        "该成分可能通过多靶点发挥作用 [Calc: seed=42]。"
        "数据不足，无法解释。"
    )
    bad_text = (
        "该方显著改善 PD 症状并发挥关键神经保护作用。"
        "引用了不存在的文件 [Data: 不存在的表.csv]。"
    )
    # 变形表达：应被 Token 解析接受，不误杀
    variant = "拓扑分析显示核心节点 [Data: ./output/PPI_拓扑.csv 的拓扑分析]。"
    rv, bv = validate(variant, datadir=os.path.join(HERE, "..", "examples"))
    r1, b1 = validate(good)
    # 传一个空 datadir 触发"文件不存在"校验，验证第二道防线
    r2, b2 = validate(bad_text, datadir=os.path.join(HERE, "..", "data"))
    ok = (r1["verdict"] == "PASS" and r2["verdict"] == "FAIL"
          and r1["sentences"] == 4 and r2["bad_count"] >= 2
          and rv["verdict"] in ("PASS", "FAIL") and any("PPI_拓扑.csv" in x for x in rv["data_refs"]))
    print("selftest: %s (good=%s/%d句, bad=%s/%d问题)" %
          ("PASS" if ok else "FAIL", r1["verdict"], r1["sentences"], r2["verdict"], r2["bad_count"]))
    if not ok:
        print("  good bad:", b1)
        print("  bad bad:", b2)
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_evidence: LLM 证据锚点约束与校验")
    ap.add_argument("cmd", choices=["prompt", "validate", "selftest"])
    ap.add_argument("--text", help="待校验文本文件或直接文本")
    ap.add_argument("--datadir", help="本地数据目录（校验 Data 引用存在性）")
    ap.add_argument("--check-pmid", action="store_true", help="联网核对 PMID（慢）")
    ap.add_argument("--out", help="校验报告 JSON 输出")
    a = ap.parse_args(argv)
    if a.cmd == "prompt":
        print(ANCHOR_PROMPT)
        return 0
    if a.cmd == "selftest":
        return selftest()
    if not a.text:
        ap.error("validate 需 --text")
    text = open(a.text, encoding="utf-8").read() if os.path.isfile(a.text) else a.text
    rep, bad = validate(text, a.datadir, a.check_pmid)
    if a.out:
        import json
        # v3.0.2: 输出路径的父目录不存在时自动创建，避免 FileNotFoundError
        _parent = os.path.dirname(os.path.abspath(a.out))
        if _parent:
            os.makedirs(_parent, exist_ok=True)
        json.dump({"report": rep, "bad": bad}, open(a.out, "w", encoding="utf-8"),
                  ensure_ascii=False, indent=1)
    print("[evidence] %s | %d/%d 句带锚点 | 数据引用 %d | PMID %d | 问题 %d" %
          (rep["verdict"], rep["anchored"], rep["sentences"],
           len(rep["data_refs"]), len(rep["pmid_refs"]), rep["bad_count"]))
    for b in bad[:20]:
        print("  ✗ 句%d：%s" % (b["idx"], b["reason"]))
        print("    %s" % b["sentence"][:90])
        print("    → %s" % b["fix"])
    return 0 if rep["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
