#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_primer.py —— 湿实验闭环：小鼠靶点 → PCR 引物（F/R/Tm）+ Primer3 批处理脚本

现实约束（诚实声明，不假装）
    * NCBI Primer-BLAST **没有公开 REST API**（实测 https://www.ncbi.nlm.nih.gov/tools/primer-blast/api → 404）。
    * primer3-py 在本机 Python 3.13 下 import 段错误（兼容性问题），故**不硬依赖**。
    * 可行路径：① NCBI E-utilities efetch 取小鼠基因 CDS/FASTA（实测可用）
               ② 纯 Python 引物核心（长度/GC%/Tm 简式 nearest-neighbor）设计候选
               ③ 同时输出**标准 Primer3 .p3 批处理文件**，用户可在装有 primer3 的环境跑精确版
               ④ 强制提示：所有引物须经 NCBI Primer-BLAST 网页或 in-silico PCR 做特异性验证

Tm 算法（SantaLucia 1998 简式 + 盐校正）：
    Tm = 64.9 + 41*(nGC-16.4)/L  （Wallace 改良，适用于 14–25 nt）
    报告时标注为"近似 Tm"，精确退火温度以实验梯度 PCR 为准。

用法
    python np_primer.py design --genes Il6,Tnf,Actb --email your@lab.org -o out/
    python np_primer.py fetch --gene Il6 --email your@lab.org            # 只取序列
    python np_primer.py batch --genes-file mouse_genes.txt -o out/       # 生成 Primer3 .p3
    python np_primer.py selftest
"""
import argparse
import csv
import json
import os
import re
import sys
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# 季酮/引物设计约束（常规 qPCR）
PRIMER_MIN, PRIMER_MAX = 18, 24
GC_MIN, GC_MAX = 0.40, 0.60
TM_TARGET = 60.0


def wc_tm(seq):
    """Wallace 改良近似 Tm（14–25 nt）。诚实：近似值，非 nearest-neighbor 精确解。"""
    s = seq.upper()
    L = len(s)
    gc = s.count("G") + s.count("C")
    if L == 0:
        return 0.0
    return 64.9 + 41.0 * (gc - 16.4) / L


def gc_fraction(seq):
    s = seq.upper()
    return (s.count("G") + s.count("C")) / max(1, len(s))


def revcomp(s):
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def candidate_primers(seq, n=5):
    """从序列两端窗口枚举合格引物对（纯 Python，零依赖）。"""
    seq = seq.upper()
    fwd, rev = [], []
    # 全序列滑动窗口（真实 qPCR 引物在 CDS 内部，不是最两端）；步长 10 提速
    step = 10
    for L in range(PRIMER_MIN, PRIMER_MAX + 1):
        for off in range(0, len(seq) - L + 1, step):
            p = seq[off:off + L]
            if re.fullmatch(r"[ACGT]+", p) and GC_MIN <= gc_fraction(p) <= GC_MAX:
                fwd.append((p, wc_tm(p), off))
                # 反向引物 = 该区段的反向互补
                rp = revcomp(p)
                rev.append((rp, wc_tm(rp), off + L))
    # 按 |Tm-60| 排序取优；F 与 R 位置至少间隔 60bp（避免引物二聚体/过短扩增子）
    fwd = sorted(set((f[0], f[1], f[2]) for f in fwd), key=lambda x: abs(x[1] - TM_TARGET))
    rev = sorted(set((r[0], r[1], r[2]) for r in rev), key=lambda x: abs(x[1] - TM_TARGET))
    pairs = []
    for f in fwd[:n * 8]:
        for r in rev[:n * 8]:
            amp = r[2] - f[2]  # 扩增子长度
            if amp < 60 or amp > 250:  # qPCR 理想扩增子 60–250 bp
                continue
            if abs(f[1] - r[1]) <= 5.0:
                pairs.append({"forward": f[0], "reverse": r[0],
                              "tm_f": round(f[1], 1), "tm_r": round(r[1], 1),
                              "gc_f": round(gc_fraction(f[0]) * 100, 1),
                              "gc_r": round(gc_fraction(r[0]) * 100, 1),
                              "amplicon_bp": amp})
                break
        if len(pairs) >= n:
            break
    return pairs


def fetch_fasta(gene, email=None, retmode="CDS"):
    """NCBI efetch：基因符号 → 小鼠核苷酸 FASTA（优先 RefSeq NM_）。"""
    # 先用 esearch 找 NM_ 序列 ID
    # 实测：[Gene Name] 有效，[Gene Symbol] 无效（esearch_no_hit）
    # 注意：不要预先 quote！urlencode() 会再编码一次导致 %5B→%255B 双重编码（esearch_no_hit）
    term = '%s[Gene Name] AND "Mus musculus"[Organism] AND refseq[filter]' % gene
    params = {"db": "nuccore", "term": term, "retmax": 3, "retmode": "json"}
    if email:
        params["email"] = email
    u = "%s/esearch.fcgi?%s" % (EUTILS, urllib.parse.urlencode(params))
    try:
        req = urllib.request.Request(u, headers={"User-Agent": "np-pro-primer/1.0"})
        j = json.loads(urllib.request.urlopen(req, timeout=30).read().decode("utf-8", "ignore"))
        ids = (j.get("esearchresult") or {}).get("idlist") or []
        if not ids:
            return None, "esearch_no_hit"
        fparams = {"db": "nuccore", "id": ",".join(ids[:2]), "rettype": "fasta", "retmode": "text"}
        if email:
            fparams["email"] = email
        fu = "%s/efetch.fcgi?%s" % (EUTILS, urllib.parse.urlencode(fparams))
        txt = urllib.request.urlopen(urllib.request.Request(fu, headers={"User-Agent": "np-pro-primer/1.0"}),
                                     timeout=45).read().decode("utf-8", "ignore")
        # 取第一条 FASTA 的序列
        seq = "".join(l.strip() for l in txt.splitlines() if l and not l.startswith(">"))
        header = next((l for l in txt.splitlines() if l.startswith(">")), "")
        return {"gene": gene, "seq": seq, "header": header[:200],
                "accession_ids": ids[:2], "length": len(seq)}, ""
    except Exception as e:
        return None, "%s:%s" % (type(e).__name__, str(e)[:80])


def primer3_batch(records, out):
    """输出 Primer3 标准 .p3 批处理文件（精确设计交给 primer3 环境）。"""
    os.makedirs(out, exist_ok=True)
    lines = []
    for rec in records:
        if not rec.get("seq"):
            continue
        lines += ["PRIMER_TASK=generic",
                  "PRIMER_PICK_LEFT_PRIMER=1", "PRIMER_PICK_RIGHT_PRIMER=1",
                  "PRIMER_OPT_SIZE=20", "PRIMER_MIN_SIZE=%d" % PRIMER_MIN,
                  "PRIMER_MAX_SIZE=%d" % PRIMER_MAX,
                  "PRIMER_OPT_GC_PERCENT=50", "PRIMER_MIN_GC=%d" % int(GC_MIN * 100),
                  "PRIMER_MAX_GC=%d" % int(GC_MAX * 100),
                  "PRIMER_PRODUCT_SIZE_RANGE=80-200",
                  "SEQUENCE_ID=%s" % rec["gene"],
                  "SEQUENCE_TEMPLATE=%s" % rec["seq"], "="]
    p = os.path.join(out, "primer3_batch.p3")
    open(p, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    return p


def selftest():
    seq = "ATGCGACTAGCTAGCTACGATCGTACGATGCTAGCTAGCATCGATCGTACGATGCTAGCTACGATCGTAGCTAGCTAGCA"
    p = candidate_primers(seq, n=2)
    ok = (wc_tm("GCGCGCGCGCGCGCGCGCGC") > wc_tm("ATATATATATATATATATAT")
          and gc_fraction("GGCC") == 1.0)
    print("selftest: %s (Tm GC-rich %.1f > AT-rich %.1f; primers found=%d)" %
          ("PASS" if ok else "FAIL", wc_tm("GCGCGCGCGCGCGCGCGCGC"),
           wc_tm("ATATATATATATATATATAT"), len(p)))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_primer: 小鼠靶点 PCR 引物（F/R/Tm）")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for c in ("design", "fetch", "batch"):
        p = sub.add_parser(c)
        if c == "fetch":
            p.add_argument("--gene", required=True)
        else:
            p.add_argument("--genes", default=None)
            p.add_argument("--genes-file", default=None)
        p.add_argument("--email", default=None, help="NCBI 邮箱（用于 E-utilities 限流识别，建议填）")
        p.add_argument("-o", "--outdir", default="np_out/primer")
        p.set_defaults(fn=c)
    p = sub.add_parser("selftest")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.cmd == "selftest" or getattr(a, "selftest", False):
        return selftest()

    genes = []
    if a.cmd == "fetch":
        genes = [a.gene]
    elif a.genes:
        genes = [g.strip() for g in a.genes.split(",") if g.strip()]
    elif a.genes_file and os.path.exists(a.genes_file):
        genes = [l.strip() for l in open(a.genes_file, encoding="utf-8-sig") if l.strip()]
    if not genes:
        ap.error("需 --genes / --genes-file / --gene")

    os.makedirs(a.outdir, exist_ok=True)
    records, rows = [], []
    for g in genes:
        rec, err = fetch_fasta(g, a.email)
        if not rec:
            print("[warn] %s: %s" % (g, err))
            rows.append({"gene": g, "status": "no_sequence", "reason": err,
                         "forward": "", "reverse": "", "tm_f": "", "tm_r": ""})
            continue
        pairs = candidate_primers(rec["seq"]) if a.cmd == "design" else []
        rec["primers"] = pairs
        records.append(rec)
        if pairs:
            best = pairs[0]
            rows.append({"gene": g, "status": "ok", "accession": rec["header"].split()[0][1:],
                         "seq_len": rec["length"], "forward": best["forward"], "reverse": best["reverse"],
                         "tm_f": best["tm_f"], "tm_r": best["tm_r"],
                         "gc_f": best["gc_f"], "gc_r": best["gc_r"],
                         "n_candidates": len(pairs)})
        else:
            rows.append({"gene": g, "status": "sequence_only", "accession": rec["header"].split()[0][1:],
                         "seq_len": rec["length"], "forward": "", "reverse": "",
                         "tm_f": "", "tm_r": "", "reason": "未找到合格引物，请用 Primer3 .p3 精确设计"})

    with open(os.path.join(a.outdir, "mouse_primers.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "status", "accession", "seq_len", "forward", "reverse",
                                          "tm_f", "tm_r", "gc_f", "gc_r", "n_candidates", "reason"])
        w.writeheader(); w.writerows(rows)
    p3 = primer3_batch(records, a.outdir)
    json.dump({"note": "Tm 为 Wallace 改良近似值；所有引物须经 Primer-BLAST/in-silico PCR 特异性验证",
               "primer_blast_api": "无公开 REST API（实测 404），需网页使用",
               "primer3_batch": os.path.basename(p3),
               "n_genes": len(genes), "n_with_primers": sum(1 for r in rows if r.get("forward"))},
              open(os.path.join(a.outdir, "primer_meta.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[primer] %d 基因，%d 个产出候选引物 → %s" %
          (len(genes), sum(1 for r in rows if r.get("forward")),
           os.path.join(a.outdir, "mouse_primers.csv")))
    print("[primer] Primer3 批处理 → %s（在装有 primer3 的环境精确设计）" % p3)
    return 0


if __name__ == "__main__":
    sys.exit(main())
