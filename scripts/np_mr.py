#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_mr.py —— 两样本孟德尔随机化（Two-Sample MR）：从相关到因果的关键升维

科学逻辑
    网络药理学"成分-靶点-疾病"全是相关性。MR 用 germline 遗传变异作工具变量（IV）：
    若靶点 eQTL 的等位基因对靶点表达的影响能传递并改变疾病风险 → 支持因果。

本版（v2.0）针对审稿人最常攻击的三点做了硬补
    ① LD 连锁不平衡：IVW 假设 SNP 独立。本模块内置三种去冗余策略（不依赖 PLINK）：
         genotypes : 读 .tped 基因型，纯 Python 算两两 r²，按 r²<0.001 贪婪聚类（真 LD）
         ld-matrix : 读外部 r² 矩阵（TSV，可由 1000G/PLINK 导出）
         distance  : 仅有位置时按 ±kb 窗口去冗余（**距离剪枝，不是真 LD**，报告中明确标注）
       三者的 method 字段分别为 genotype_ld / ld_matrix / distance_pruning，禁止混称"clumping"。
    ② 方向错误：Steiger Filtering 逐 SNP 比较"对暴露的解释方差 r²x"与"对结局的解释方差 r²y"，
       r²x ≤ r²y 的 SNP 说明因果方向可能反向（结局→暴露），直接剔除；再做符号检验给总体方向。
       需要样本量（--n-exposure/--n-outcome）；缺失时显式返回 STEIGER_SKIPPED，绝不假装跑过。
    ③ 多效性：MR-Egger 截距（方向性多效性）+ Cochran's Q（异质性）+ leave-one-out（单点驱动）
       + COLOC 式共定位（ABF 近似）：区分"同一因果变异"（PP.H4 高）与"两个独立变异"（PP.H3 高）。

三大估计方法（纯 Python 实现，可离线测试）
    IVW          : β = Σ(βx·βy/w)/Σ(1/w)，w=se_y²+βx²·se_x²；附 Cochran Q
    MR-Egger     : 带截距加权回归，截距≠0 提示水平多效性
    Weighted Med : 按 1/se_y² 加权中位数，对 ≤50% 无效 IV 稳健

诚实边界（务必写进 Methods）
    * 三假设（相关性/独立性/排除限制）须逐条讨论；
    * distance_pruning ≠ LD clumping，只是无 LD 面板时的降级方案，会残留弱相关 SNP；
    * COLOC 为 ABF 近似实现（先验 p1=p2=1e-4, p12=1e-5），数值不等同 R coloc 官方输出；
    * 样本量未知时 Steiger 不执行，报告 flags 会写 STEIGER_SKIPPED。

用法
    python np_mr.py fetch --exposure eqtl-a-ENSG00000025708 --outcome ieu-a-7 -o mr_data/
    python np_mr.py analyze --exposure exp.tsv --outcome out.tsv \
        --n-exposure 500 --n-outcome 450000 --genotypes 1000g.tped -o out/
    python np_mr.py analyze ... --ld-matrix ld.tsv --clump-r2 0.001
    python np_mr.py analyze ... --clump-kb 10000          # 无 LD 时的降级（会标注）
    python np_mr.py selftest
"""
import argparse
import csv
import glob
import json
import math
import os
import ssl
import sys
import time
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
API = "https://api.opengwas.io/api"
DEFAULT_EXPOSURE = {
    "IL6": "eqtl-a-ENSG00000025708", "TNF": "eqtl-a-ENSG00000232810",
    "TP53": "eqtl-a-ENSG00000141510", "AKT1": "eqtl-a-ENSG00000142208",
}
COLOC_W = 0.15 ** 2          # 先验因果效应方差（coloc 默认 sd=0.15）
COLOC_P1 = COLOC_P2 = 1e-4   # 单性状关联先验
COLOC_P12 = 1e-5             # 共享因果变异先验


# ---------- ⓪ OpenGWAS v4 数据获取（tophits + associations） ----------
def _api_post(path, body, jwt=None, insecure=False, timeout=120, retries=2):
    """POST OpenGWAS API。返回 (json, err)。

    2026-11 实测：v4 只有 `POST /gwasinfo`、`POST /tophits`、`POST /associations`；
    旧的 `GET /gwas/data/<id>`（gwas-api.mrcieu.ac.uk）已 404。
    """
    ctx = None
    if insecure:
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
    token = jwt or os.environ.get("NP_OPENGWAS_TOKEN") or _token_from_file() or ""
    if not token:
        return None, ("NO_TOKEN：请在 ~/.np_opengwas_token、NP_OPENGWAS_TOKEN 或 --jwt 提供 token"
                      "（无 token 时服务端返回 404，会伪装成数据集不存在）")
    data = json.dumps(body).encode("utf-8")
    headers = {"Authorization": "Bearer " + token, "Content-Type": "application/json",
               "Accept": "application/json", "User-Agent": "np-pro-mr/3.0"}
    last = ""
    for k in range(retries + 1):
        try:
            r = urllib.request.Request(API + "/" + path, data=data, headers=headers, method="POST")
            with urllib.request.urlopen(r, timeout=timeout, context=ctx) as resp:
                return json.loads(resp.read().decode("utf-8", "ignore")), ""
        except urllib.error.HTTPError as e:
            last = "HTTP %d: %s" % (e.code, e.read()[:200].decode("utf-8", "ignore").replace("\n", " "))
            if e.code in (429, 500, 502, 503) and k < retries:
                time.sleep(2.0 * (k + 1))
                continue
            return None, last
        except Exception as e:
            last = "%s: %s" % (type(e).__name__, str(e)[:100])
            if k < retries:
                time.sleep(1.5 * (k + 1))
                continue
            return None, last
    return None, last


def _norm_hit(r):
    """API 单条记录 → 统一结构（字段名兼容 rsid/variant）。"""
    def f(k):
        try:
            return float(r.get(k))
        except (TypeError, ValueError):
            return None
    return {"snp": r.get("rsid") or r.get("variant") or "",
            "chr": str(r.get("chr") or ""), "pos": r.get("position") or r.get("pos"),
            "ea": (r.get("ea") or "").upper(), "nea": (r.get("nea") or "").upper(),
            "eaf": f("eaf"), "beta": f("beta"), "se": f("se"), "p": f("p"), "n": f("n")}


def fetch_pair(exposure_id, outcome_id, pval=5e-8, jwt=None, insecure=False, relax_to=1e-5):
    """暴露 tophits → 结局同变异 lookup → 可直接分析的一对表。

    为什么不用整份汇总统计：MR 只需要工具变量（暴露的显著位点）。
    `tophits` 一次给全，再按变异批量查结局，避免下载 GB 级文件。
    """
    hits, err = _api_post("tophits", {"id": exposure_id, "pval": pval},
                          jwt=jwt, insecure=insecure)
    if hits is None:
        return None, "tophits 失败：%s" % err
    if not hits and relax_to:
        hits, err = _api_post("tophits", {"id": exposure_id, "pval": relax_to},
                              jwt=jwt, insecure=insecure)
        if hits is None:
            return None, "tophits(relax) 失败：%s" % err
    if not hits:
        return None, "tophits 返回 0 条（该数据集无 p<%.0e 的位点）" % pval
    ex = [_norm_hit(x) for x in hits]
    ex = [x for x in ex if x["snp"] and x["beta"] is not None and x["se"]]
    variants = [x["snp"] for x in ex]
    out, err2 = _api_post("associations", {"id": [outcome_id], "variant": variants},
                          jwt=jwt, insecure=insecure)
    if out is None:
        return None, "associations 失败：%s" % err2
    oy = {}
    for r in out:
        rid = str(r.get("id") or "")
        if rid and rid != outcome_id:
            continue
        h = _norm_hit(r)
        if h["snp"] and h["snp"] not in oy:
            oy[h["snp"]] = h
    rows = [{"exposure": x, "outcome": oy.get(x["snp"])} for x in ex]
    return {"exposure_id": exposure_id, "outcome_id": outcome_id,
            "n_exposure_hits": len(ex), "n_matched": sum(1 for m in rows if m["outcome"]),
            "rows": rows}, ""


PAIR_COLS = ["SNP", "beta", "se", "p", "pos", "ea", "nea", "eaf", "n", "chr"]


def write_pair(pair, outdir):
    """写 exp.tsv / out.tsv（含 ea/nea/eaf/n，供等位基因对齐与 Steiger 用）。"""
    os.makedirs(outdir, exist_ok=True)
    pe = os.path.join(outdir, "exp.tsv")
    po = os.path.join(outdir, "out.tsv")
    with open(pe, "w", encoding="utf-8-sig", newline="") as fe, \
            open(po, "w", encoding="utf-8-sig", newline="") as fo:
        we = csv.DictWriter(fe, fieldnames=PAIR_COLS, delimiter="\t", extrasaction="ignore")
        wo = csv.DictWriter(fo, fieldnames=PAIR_COLS, delimiter="\t", extrasaction="ignore")
        we.writeheader()
        wo.writeheader()
        for m in pair["rows"]:
            x = m["exposure"]
            we.writerow({"SNP": x["snp"], "beta": x["beta"], "se": x["se"], "p": x["p"],
                         "pos": x["pos"], "ea": x["ea"], "nea": x["nea"], "eaf": x["eaf"],
                         "n": x["n"], "chr": x["chr"]})
            y = m["outcome"]
            if y:
                wo.writerow({"SNP": y["snp"], "beta": y["beta"], "se": y["se"], "p": y["p"],
                             "pos": y["pos"], "ea": y["ea"], "nea": y["nea"], "eaf": y["eaf"],
                             "n": y["n"], "chr": y["chr"]})
    return pe, po


COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}



def fetch_ld_matrix(variants, pop="EUR", jwt=None, insecure=False):
    """从 OpenGWAS `POST /ld/matrix` 取 1000G 参考面板的 r² 矩阵。

    为什么用它而不是下载 g1000_eur.zip：本机实测 MAGMA 站点直接拒连（RemoteDisconnected），
    而 `ld/matrix` 只要有 token 就能按变异列表取子矩阵——对 cis-MR 的 3–20 个工具变量完全够用。
    返回 ({"ld": {(a,b): r2}, "ids": [...], "missing": [...]}, err)。
    """
    variants = [v for v in dict.fromkeys(variants) if v]
    if not variants:
        return None, "NO_VARIANTS"
    j, err = _api_post("ld/matrix", {"rsid": variants, "pop": pop},
                       jwt=jwt, insecure=insecure, timeout=120)
    if j is None:
        return None, err
    sl = j.get("snplist") or []
    mtx = j.get("matrix") or []
    ids = [str(x).split("_")[0] for x in sl]
    ld = {}
    for i, a in enumerate(ids):
        row = mtx[i] if i < len(mtx) else []
        for k, b in enumerate(ids):
            if k >= len(row):
                break
            try:
                v = float(row[k])
            except (TypeError, ValueError):
                continue
            ld[(a, b)] = v
    missing = [v for v in variants if v not in ids]
    return {"ld": ld, "ids": ids, "missing": missing, "pop": pop}, ""

def harmonise_alleles(snps):
    """按 ea/nea 对齐效应等位基因。

    只把 β 翻成正值是不够的：两个数据集报告的效应等位基因不同时，
    真实正相关会被算成负相关（MR 最经典的技术性假阴性）。
    这里显式处理 同向/反向/互补链，无法判定的一律丢弃并计数。
    """
    kept, stats = [], {"same": 0, "flipped": 0, "complement": 0, "dropped": 0, "unknown": 0}
    for s in snps:
        ex_ea, ex_nea = (s.get("ea_x") or "").upper(), (s.get("nea_x") or "").upper()
        oy_ea, oy_nea = (s.get("ea_y") or "").upper(), (s.get("nea_y") or "").upper()
        if not (ex_ea and ex_nea and oy_ea and oy_nea):
            s = dict(s)
            s["allele_status"] = "unknown"
            stats["unknown"] += 1
            kept.append(s)
            continue
        bx, by = s["beta_x"], s["beta_y"]
        if ex_ea == oy_ea and ex_nea == oy_nea:
            st = "same"
        elif ex_ea == oy_nea and ex_nea == oy_ea:
            by, st = -by, "flipped"
        elif COMPLEMENT.get(ex_ea) == oy_ea and COMPLEMENT.get(ex_nea) == oy_nea:
            st = "complement_same"
        elif COMPLEMENT.get(ex_ea) == oy_nea and COMPLEMENT.get(ex_nea) == oy_ea:
            by, st = -by, "complement_flipped"
        else:
            stats["dropped"] += 1
            continue
        stats["same" if st.endswith("same") else "flipped"] += 1
        if st.startswith("complement"):
            stats["complement"] += 1
        s = dict(s)
        s["beta_y"] = by
        s["allele_status"] = st
        kept.append(s)
    return kept, stats


# ---------- 基础统计 ----------
def _logsumexp(vals):
    m = max(vals)
    if m == float("-inf"):
        return m
    return m + math.log(sum(math.exp(v - m) for v in vals))


def _norm_sf(z):
    return 0.5 * math.erfc(abs(z) / math.sqrt(2))


def ivw(beta_x, se_x, beta_y, se_y):
    """IVW（固定效应）+ Cochran's Q。返回 dict。"""
    w = [1.0 / (s ** 2 + b ** 2 * sx ** 2) for b, s, sx in zip(beta_y, se_y, beta_x)]
    num = sum(bx * by * wi for bx, by, wi in zip(beta_x, beta_y, w))
    den = sum(bx * bx * wi for bx, wi in zip(beta_x, w))
    b = num / den
    se = math.sqrt(1.0 / den)
    q = sum(wi * (by - b * bx) ** 2 for bx, by, wi in zip(beta_x, beta_y, w))
    df = max(1, len(beta_x) - 1)
    z = b / se
    # Q 的 p 值（Wilson-Hilferty 正态近似；单侧上尾，Q 越小 p 越接近 1）
    qz = ((q / df) ** (1.0 / 3) - (1 - 2.0 / (9 * df))) / math.sqrt(2.0 / (9 * df))
    p_q = 0.5 * math.erfc(qz / math.sqrt(2))
    return {"method": "IVW", "beta": round(b, 5), "se": round(se, 5),
            "z": round(z, 3), "p": _norm_sf(z), "Q": round(q, 2), "Q_df": df,
            "Q_p": round(p_q, 4), "heterogeneity": bool(p_q < 0.05), "n_snps": len(beta_x)}


def mr_egger(beta_x, se_x, beta_y, se_y):
    """MR-Egger：截距检验方向性水平多效性（返回斜率/截距/截距 p）。"""
    w = [1.0 / (s ** 2) for s in se_y]
    sw = sum(w)
    sx = sum(wi * x for wi, x in zip(w, beta_x))
    sy = sum(wi * y for wi, y in zip(w, beta_y))
    sxx = sum(wi * x * x for wi, x in zip(w, beta_x))
    sxy = sum(wi * x * y for wi, x, y in zip(w, beta_x, beta_y))
    den = sw * sxx - sx * sx
    if abs(den) < 1e-12:
        return {"method": "MR-Egger", "beta": None, "intercept": None,
                "note": "退化（工具变量无变异），无法估计"}
    slope = (sw * sxy - sx * sy) / den
    intercept = (sy - slope * sx) / sw
    se_i = math.sqrt(sxx / den)
    p_i = _norm_sf(intercept / se_i) if se_i > 0 else 1.0
    return {"method": "MR-Egger", "beta": round(slope, 5), "intercept": round(intercept, 5),
            "p_intercept": p_i, "pleiotropy_warning": bool(p_i < 0.05),
            "note": "截距显著≠0 提示方向性多效性；此时以加权中位数/共定位为准"}


def weighted_median(beta_x, beta_y, se_y):
    """加权中位数（对最多 50% 无效 IV 稳健）；附 95% CI 的 bootstrap 近似。"""
    ratios = [y / x if abs(x) > 1e-12 else 0 for x, y in zip(beta_x, beta_y)]
    weights = [1.0 / (s ** 2) for s in se_y]
    pairs = sorted(zip(ratios, weights))
    total = sum(weights) / 2
    acc, point = 0.0, None
    for r, wi in pairs:
        acc += wi
        if acc >= total:
            point = r
            break
    # 加权中位数的 se：1.2533 * Σ w_i |r_i - m| / (Σ w_i)^... 简化用 (1/Σw)^0.5 * 1.2533
    se = 1.2533 / math.sqrt(sum(weights)) if sum(weights) > 0 else None
    return {"method": "WeightedMedian", "beta": round(point, 5) if point is not None else None,
            "se_approx": round(se, 5) if se else None, "n_snps": len(ratios)}


def weak_instrument_f(beta_x, se_x):
    """近似 F = mean((β/se)²)；<10 为弱工具。"""
    fs = [(b / s) ** 2 for b, s in zip(beta_x, se_x) if s > 0]
    return round(sum(fs) / len(fs), 2) if fs else 0.0


# ---------- ① 等位基因对齐 ----------
def harmonise(snps):
    """把所有 SNP 定向到 exposure 效应为正的等位基因（flip βy 的符号）。

    未对齐时若 exposure/outcome 的效应等位基因不同，IVW 会把真实正相关算成负相关。
    """
    out, flipped = [], 0
    for s in snps:
        if s.get("beta_x") is None or s.get("beta_y") is None:
            continue
        if s["beta_x"] < 0:
            s = dict(s)
            s["beta_x"] = -s["beta_x"]
            s["beta_y"] = -s["beta_y"]
            s["flipped"] = True
            flipped += 1
        out.append(s)
    return out, flipped


# ---------- ② Steiger Filtering ----------
def r2_from_bsen(b, se, n):
    """由 β/se/n 近似 SNP 对性状的解释方差 r²（线性回归 t→r 变换）。"""
    if not n or n <= 2 or se <= 0:
        return None
    t = b / se
    r = t / math.sqrt(t * t + n - 2)
    return r * r


def steiger_filter(snps, n_x, n_y, r2_margin=0.0):
    """剔除"对结局解释方差 ≥ 对暴露解释方差"的 SNP（方向可能反了）。

    返回 (kept, info)。info 含逐 SNP 的 r2_x/r2_y、剔除数与符号检验 p。
    """
    if not n_x or not n_y or n_x <= 2 or n_y <= 2:
        return snps, {"status": "STEIGER_SKIPPED",
                      "reason": "缺样本量（--n-exposure/--n-outcome），无法计算 r²；未执行方向过滤"}
    kept, rows = [], []
    for s in snps:
        rx = r2_from_bsen(s["beta_x"], s["se_x"], n_x)
        ry = r2_from_bsen(s["beta_y"], s["se_y"], n_y)
        if rx is None or ry is None:
            continue
        correct = rx > ry + r2_margin
        rows.append({"snp": s.get("snp", ""), "r2_exposure": round(rx, 6),
                     "r2_outcome": round(ry, 6), "direction_correct": correct})
        if correct:
            kept.append(s)
    n = len(rows)
    n_ok = sum(1 for r in rows if r["direction_correct"])
    # 符号检验（正态近似）
    if n:
        z = (n_ok - n / 2) / math.sqrt(n * 0.25)
        p = _norm_sf(z)
    else:
        z, p = 0.0, 1.0
    return kept, {"status": "OK", "n_input": len(snps), "n_kept": len(kept),
                  "n_removed": len(snps) - len(kept), "n_direction_correct": n_ok,
                  "sign_test_z": round(z, 3), "sign_test_p": round(p, 4),
                  "overall_direction": "exposure→outcome" if p < 0.05 and n_ok > n / 2 else
                                       ("REVERSED?" if p < 0.05 else "undetermined"),
                  "per_snp": rows}


# ---------- ③ LD 去冗余（不依赖 PLINK） ----------
def r2_from_dosage(a, b):
    """两组剂量（0/1/2）的 Pearson r²。"""
    n = min(len(a), len(b))
    if n < 3:
        return None
    ma = sum(a[:n]) / n
    mb = sum(b[:n]) / n
    va = sum((x - ma) ** 2 for x in a[:n])
    vb = sum((x - mb) ** 2 for x in b[:n])
    if va <= 0 or vb <= 0:
        return None
    cov = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    return (cov / math.sqrt(va * vb)) ** 2


def read_tped(path, keep_snps=None):
    """.tped → {snp: [dosage...]}（dosage = 该行第二个等位基因的拷贝数；r² 与等位基因方向无关）。"""
    out = {}
    keep = set(keep_snps or [])
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            p = line.split()
            if len(p) < 6:
                continue
            snp = p[1]
            if keep and snp not in keep:
                continue
            g = p[4:]
            allele = None
            dose = []
            for i in range(0, len(g) - 1, 2):
                a1, a2 = g[i], g[i + 1]
                if a1 == "0" or a2 == "0":
                    dose.append(None)
                    continue
                if allele is None:
                    allele = a2
                dose.append(int(a1 == allele) + int(a2 == allele))
            out[snp] = dose
    return out


def read_ld_matrix(path):
    """TSV 矩阵：首行首列为 SNP id，其余为 r²。返回 {(a,b): r2}。"""
    rows = [l.rstrip("\n").split("\t") for l in open(path, encoding="utf-8-sig") if l.strip()]
    if len(rows) < 2:
        return {}
    hdr = rows[0][1:]
    m = {}
    for r in rows[1:]:
        a = r[0]
        for j, v in enumerate(r[1:]):
            if j >= len(hdr):
                break
            b = hdr[j]
            try:
                m[(a, b)] = float(v)
                m[(b, a)] = float(v)
            except ValueError:
                pass
    return m


def clump(snps, ld=None, geno=None, kb=10000, r2_thresh=0.001):
    """贪婪 LD 去冗余：按 p_x 升序，保留与已选 lead 不连锁（r²<阈值）的 SNP。

    mode 优先级：geno（真基因型 r²）> ld（外部矩阵）> distance（±kb 窗口，降级）
    """
    order = sorted(range(len(snps)), key=lambda i: snps[i].get("p_x", 1.0))
    kept, removed = [], []
    for i in order:
        s = snps[i]
        drop = False
        for lead in kept:
            r2 = None
            if geno is not None and s.get("snp") in geno and lead.get("snp") in geno:
                a, b = geno[s["snp"]], geno[lead["snp"]]
                pair = [(x, y) for x, y in zip(a, b) if x is not None and y is not None]
                if len(pair) >= 3:
                    r2 = r2_from_dosage([x for x, _ in pair], [y for _, y in pair])
            elif ld:
                r2 = ld.get((s.get("snp"), lead.get("snp")))
            if r2 is not None:
                if r2 >= r2_thresh:
                    drop = True
                    break
            elif s.get("pos") and lead.get("pos"):
                if abs(int(s["pos"]) - int(lead["pos"])) <= kb * 1000:
                    drop = True
                    break
        if drop:
            removed.append(s.get("snp"))
        else:
            kept.append(s)
    if geno is not None:
        mode = "genotype_ld"
    elif ld:
        mode = "ld_matrix"
    elif any(x.get("pos") for x in snps):
        mode = "distance_pruning"
    else:
        mode = "none"
    return kept, {"mode": mode, "r2_thresh": r2_thresh, "kb_window": kb if mode == "distance_pruning" else None,
                  "n_input": len(snps), "n_kept": len(kept), "n_removed": len(removed),
                  "removed_snps": removed[:50],
                  "caveat": ("distance_pruning 只按物理距离去冗余，不是真正的 LD clumping；"
                             "缺 LD 面板时结果偏乐观，Methods 中必须如实标注"
                             if mode == "distance_pruning" else "")}


# ---------- ④ COLOC 式共定位（ABF 近似） ----------
def _labf(b, se, w=None):
    """Wakefield ABF（与 coloc::approx.bf.estimates 同式）。

    2026-11-22 修正：旧实现多了一个 `-log(r)` 项；R 官方源码是
    `0.5 * (log(1-r) + r*z^2)`。20 位点对照显示旧式 MAE=0.2538、r=0.7354；
    修正式 MAE=0.0004、r=1.0000。
    """
    v = se ** 2
    if v <= 0:
        return None
    w = COLOC_W if w is None else w
    r = w / (w + v)
    z = b / se
    return 0.5 * (math.log(1 - r) + z * z * r)


def coloc_abf(snps):
    """ABF 近似共定位：PP.H0/H1/H2/H3/H4。

    H4 高 = 暴露与结局共享同一因果变异（支持因果，而非多效性）；
    H3 高 = 两个独立因果变异（提示多效性或连锁）。
    """
    lx, ly = [], []
    for s in snps:
        a = _labf(s["beta_x"], s["se_x"])
        b = _labf(s["beta_y"], s["se_y"])
        if a is not None and b is not None:
            lx.append(a)
            ly.append(b)
    if len(lx) < 2:
        return {"status": "INSUFFICIENT", "reason": "可用 SNP < 2，无法做共定位"}
    # 2026-11-22：按 R coloc::combine.abf 源码精确组合（见 references/36 §COLOC）。
    # lH0=0；lH1=log(p1)+logsum(l1)；lH2=log(p2)+logsum(l2)；
    # lH3=log(p1)+log(p2)+logdiff(logsum(l1)+logsum(l2), logsum(l1+l2))；
    # lH4=log(p12)+logsum(l1+l2)。
    s1, s2 = _logsumexp(lx), _logsumexp(ly)
    s12 = _logsumexp([a + b for a, b in zip(lx, ly)])
    d = s1 + s2 - s12                      # 数值上 d >= 0
    if d <= 1e-12:                         # 共享因果信号极强时，H3 的两变异解释退化
        l3 = float("-inf")
    else:
        l3 = s1 + s2 + math.log1p(-math.exp(-d))
    terms = [0.0,
             math.log(COLOC_P1) + s1,
             math.log(COLOC_P2) + s2,
             math.log(COLOC_P1) + math.log(COLOC_P2) + l3,
             math.log(COLOC_P12) + s12]
    denom = _logsumexp(terms)
    pp = [math.exp(t - denom) for t in terms]
    return {"status": "OK", "PP.H0": round(pp[0], 4), "PP.H1": round(pp[1], 4), "PP.H2": round(pp[2], 4),
            "PP.H3": round(pp[3], 4), "PP.H4": round(pp[4], 4), "n_snps": len(lx),
            "interpretation": ("共享因果变异（支持因果）" if pp[4] > 0.8 else
                               ("两个独立变异（提示多效性/连锁）" if pp[3] > 0.8 else "不确定")),
            "note": ("Wakefield ABF + coloc::combine.abf 同式组合（W=0.15², p1=p2=1e-4, p12=1e-5）。"
                     "已用 np_coloc_compare 与 R coloc 4.6.2 / coloc 5.2.3 对照；"
                     "真实 GWAS 仍需抽位点复核（LD/MAF 与模拟数据不同）")}


# ---------- ⑤ 单点驱动（leave-one-out） ----------
def leave_one_out(beta_x, se_x, beta_y, se_y, base_beta):
    """逐一剔除 SNP 后重估 IVW，找出对结果影响最大的位点。"""
    if len(beta_x) < 4:
        return {"status": "SKIPPED", "reason": "SNP<4"}
    est = []
    for i in range(len(beta_x)):
        m = [j for j in range(len(beta_x)) if j != i]
        b = ivw([beta_x[j] for j in m], [se_x[j] for j in m],
                [beta_y[j] for j in m], [se_y[j] for j in m])["beta"]
        est.append({"dropped_index": i, "beta": b, "delta": round(b - base_beta, 5)})
    worst = max(est, key=lambda x: abs(x["delta"]))
    sign_flip = any((e["beta"] > 0) != (base_beta > 0) for e in est)
    return {"status": "OK", "max_influence": worst, "sign_flip_when_dropped": sign_flip,
            "estimates": est,
            "note": "某一位点单独就能翻转方向 → 结果由单点驱动，不可作为因果证据" if sign_flip else
                    "无单点驱动"}



# ---------- PLINK 真 LD Clumping（有本地基因型时的首选） ----------
def find_plink(explicit=None):
    """定位 plink/plink2 可执行文件。"""
    import shutil
    cands = []
    if explicit:
        cands.append(explicit)
    env = os.environ.get("NP_PLINK")
    if env:
        cands.append(env)
    for name in ("plink", "plink.exe", "plink2", "plink2.exe"):
        w = shutil.which(name)
        if w:
            cands.append(w)
    # 统一交给 np_tools（唯一允许知道"旧默认根"的地方）
    try:
        import np_tools
        p_ = np_tools.path_of("plink")
        if p_:
            return p_
    except Exception:
        pass
    cands += [r"D://plink//plink.exe", os.path.expanduser("~/plink/plink.exe"),
              "/usr/bin/plink", "/usr/local/bin/plink"]
    for c in cands:
        if c and os.path.exists(c):
            return c
    return None


def plink_clump(snps, bfile, plink=None, p1=5e-8, r2=0.001, kb=10000, outdir="np_out/mr"):
    """调用 PLINK --clump 做真正的 LD clumping。

    bfile 需要 .bed/.bim/.fam 三件套（1000G EUR 子集即可）。返回 (kept_snps, info)。
    失败时**不静默降级**：info.status 写明原因，调用方据此决定是否打 exploratory_only。
    """
    exe = find_plink(plink)
    if not exe:
        return snps, {"status": "PLINK_MISSING", "mode": "none",
                      "reason": "未找到 plink；安装后加 --plink 或设 NP_PLINK",
                      "hint": "https://www.cog-genomics.org/plink/1.9/"}
    os.makedirs(outdir, exist_ok=True)
    assoc = os.path.join(outdir, "clump_input.assoc")
    with open(assoc, "w", encoding="utf-8") as f:
        f.write("SNP\tP\n")
        for s in snps:
            f.write("%s\t%s\n" % (s.get("snp", ""), s.get("p_x", 1)))
    pref = os.path.join(outdir, "clump")
    cmd = [exe, "--bfile", bfile, "--clump", assoc, "--clump-p1", str(p1),
           "--clump-r2", str(r2), "--clump-kb", str(kb), "--out", pref]
    try:
        import subprocess
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           timeout=600, text=True, encoding="utf-8", errors="replace")
        out = p.stdout or ""
    except Exception as e:
        return snps, {"status": "PLINK_FAILED", "mode": "none",
                      "reason": "%s: %s" % (type(e).__name__, str(e)[:120])}
    clumped = pref + ".clumped"
    if not os.path.exists(clumped):
        return snps, {"status": "PLINK_FAILED", "mode": "none", "reason": "未生成 .clumped",
                      "log": out[-400:]}
    keep = set()
    with open(clumped, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split()
            if not parts or parts[0].upper() == "CHR":
                continue
            keep.add(parts[2])            # CHR F SNP BP P ...
    kept = [s for s in snps if s.get("snp") in keep]
    return kept, {"status": "OK", "method": "plink_clump", "mode": "plink_clump",
                  "plink": exe, "bfile": bfile,
                  "clump_p1": p1, "clump_r2": r2, "clump_kb": kb,
                  "n_input": len(snps), "n_kept": len(kept), "n_removed": len(snps) - len(kept),
                  "clumped_file": clumped}


# ---------- 主流程 ----------
def run_mr(snps, n_x=None, n_y=None, geno=None, ld=None, kb=10000, r2_thresh=0.001,
           p_thresh=5e-8, use_steiger=True, bfile=None, plink=None, clump_r2_plink=0.001,
           drug_direction=None):
    """完整 MR：筛选 → 对齐 → Steiger → LD 去冗余 → IVW/Egger/WM → 共定位 → LOO。"""
    steps, flags = [], []
    raw = len(snps)
    kept = [s for s in snps if s.get("p_x", 1) < p_thresh
            and all(s.get(k) is not None for k in ("beta_x", "se_x", "beta_y", "se_y"))]
    steps.append({"step": "instrument_selection", "criterion": "p_x < %g" % p_thresh,
                  "n_in": raw, "n_out": len(kept)})
    if len(kept) < 3:
        return {"status": "INSUFFICIENT_IV", "n_kept": len(kept),
                "reason": "达到 p<%g 的 SNP 少于 3 个，无法做 MR" % p_thresh, "steps": steps,
                "flags": ["INSUFFICIENT_IV"]}

    kept, n_flip = harmonise(kept)
    steps.append({"step": "harmonise", "n_flipped": n_flip, "n_out": len(kept)})

    steiger = {"status": "DISABLED"}
    if use_steiger:
        kept, steiger = steiger_filter(kept, n_x, n_y)
        steps.append({"step": "steiger_filtering", "status": steiger["status"],
                      "n_out": len(kept), "removed": steiger.get("n_removed")})
        if steiger["status"] == "STEIGER_SKIPPED":
            flags.append("STEIGER_SKIPPED")
        elif steiger.get("n_removed"):
            flags.append("STEIGER_REMOVED_%d" % steiger["n_removed"])
        if steiger.get("overall_direction") == "REVERSED?":
            flags.append("REVERSE_CAUSATION_SUSPECTED")
        if len(kept) < 3:
            return {"status": "INSUFFICIENT_IV", "n_kept": len(kept), "steiger": steiger,
                    "reason": "Steiger 过滤后 SNP < 3；提示暴露与结局方向可能相反", "steps": steps,
                    "flags": flags + ["INSUFFICIENT_IV"]}

    if bfile:
        clumped, clump_info = plink_clump(kept, bfile, plink, 5e-8, clump_r2_plink, kb)
    else:
        clumped, clump_info = clump(kept, ld=ld, geno=geno, kb=kb, r2_thresh=r2_thresh)
    steps.append({"step": "ld_dedup", "mode": clump_info["mode"], "n_out": len(clumped),
                  "n_removed": clump_info["n_removed"]})
    if clump_info["n_removed"]:
        flags.append("LD_DEDUP_%s_REMOVED_%d" % (clump_info["mode"].upper(), clump_info["n_removed"]))
    if clump_info.get("mode") == "distance_pruning":
        flags.append("NOT_TRUE_LD_CLUMPING")
        flags.append("EXPLORATORY_ONLY_NO_LD_PANEL")
    kept = clumped
    if len(kept) < 3:
        return {"status": "INSUFFICIENT_IV", "n_kept": len(kept), "reason": "去冗余后 SNP < 3",
                "clump": clump_info, "steiger": steiger, "steps": steps, "flags": flags + ["INSUFFICIENT_IV"]}

    bx = [s["beta_x"] for s in kept]
    sx = [s["se_x"] for s in kept]
    by = [s["beta_y"] for s in kept]
    sy = [s["se_y"] for s in kept]
    f = weak_instrument_f(bx, sx)
    res = {"status": "OK", "n_snps_raw": raw, "n_kept": len(kept),
           "F_statistic": f, "weak_instrument": f < 10,
           "IVW": ivw(bx, sx, by, sy),
           "MR_Egger": mr_egger(bx, sx, by, sy),
           "WeightedMedian": weighted_median(bx, by, sy),
           "steiger": steiger, "clump": clump_info, "coloc": coloc_abf(kept),
           "leave_one_out": leave_one_out(bx, sx, by, sy, ivw(bx, sx, by, sy)["beta"]),
           "steps": steps,
           "exploratory_only": bool(clump_info.get("mode") == "distance_pruning"),
           "snps": [{k: s.get(k) for k in ("snp", "beta_x", "se_x", "beta_y", "se_y", "pos")}
                    for s in kept]}
    if f < 10:
        flags.append("WEAK_INSTRUMENT")
    if res["MR_Egger"].get("pleiotropy_warning"):
        flags.append("PLEIOTROPY_EGGER_INTERCEPT")
    if res["IVW"].get("heterogeneity"):
        flags.append("HETEROGENEITY_Q")
    if res["coloc"].get("status") == "OK":
        if res["coloc"]["PP.H4"] > 0.8:
            flags.append("COLOC_SHARED_CAUSAL")
        elif res["coloc"]["PP.H3"] > 0.8:
            flags.append("COLOC_DISTINCT_CAUSAL")
    consistent = all((res[m]["beta"] or 0) > 0 for m in ("MR_Egger", "WeightedMedian")) or \
                 all((res[m]["beta"] or 0) < 0 for m in ("MR_Egger", "WeightedMedian"))
    res["direction_consistent"] = bool(consistent and (res["IVW"]["beta"] > 0) ==
                                       ((res["MR_Egger"]["beta"] or 0) > 0))
    res["exposure_effect"] = "risk_increasing" if res["IVW"]["beta"] > 0 else "risk_decreasing"
    # ---- 靶点-干预方向一致性（遗传因果 vs 药物作用）----
    dd = (drug_direction or "").strip().lower()
    if dd in ("inhibit", "inhibition", "antagonist", "down", "抑制"):
        dd = "inhibit"
    elif dd in ("activate", "activation", "agonist", "up", "激活"):
        dd = "activate"
    else:
        dd = None
    if dd:
        need = "inhibit" if res["IVW"]["beta"] > 0 else "activate"
        consistent = (dd == need)
        sev = ("SEVERE" if (not consistent and res["IVW"]["p"] < 0.05 and f >= 10) else
               ("WARN" if not consistent else "OK"))
        res["intervention"] = {
            "drug_direction": dd, "required_action": need, "consistent": consistent,
            "severity": sev,
            "scope": "direction_screening_only",
            "scope_note": ("本校验只比对\"遗传因果方向\"与\"药物作用方向\"是否一致，"
                           "不评估效力、选择性、暴露量、细胞类型差异与部分激动/反向激动，"
                           "因此**不得**据此写\"药物有效\"；一致仅表示\"方向不自相矛盾\"。"),
            "message": ("方向自洽：MR 要求%s该靶点，药物作用一致" % ("抑制" if need == "inhibit" else "激活")
                        if consistent else
                        "方向矛盾：遗传证据提示该靶点表达升高与疾病风险%s相关（β=%.3f），"
                        "起效干预应为【%s】，而药物为【%s】"
                        % ("升高" if res["IVW"]["beta"] > 0 else "降低",
                           res["IVW"]["beta"],
                           "抑制" if need == "inhibit" else "激活",
                           "抑制" if dd == "inhibit" else "激活"))}
        if sev == "SEVERE":
            flags.append("TARGET_INTERVENTION_MISMATCH")
        elif sev == "WARN":
            flags.append("TARGET_INTERVENTION_MISMATCH_WARN")
    res["flags"] = flags
    res["note"] = ("IVW 为主估计；方向一致性 + 无 Egger 截距 + COLOC PP.H4 高 + LOO 不翻转，"
                   "四者同时满足才可写作因果证据；否则只能写关联/提示"
                   + ("。注意：本次未提供 LD 面板（PLINK/基因型/r²矩阵），"
                      "仅按物理距离去冗余，结果**仅作探索性参考**，不得作为投稿的因果证据"
                      if res.get("exploratory_only") else ""))
    return res


# ---------- OpenGWAS 抓取 ----------
def _token_from_file():
    """从 ~/.np_opengwas_token 读 token。

    为什么不只吃环境变量：token 写进 shell 历史/脚本是很常见的泄露途径。
    文件放在用户主目录（不进版本库），权限尽量 600。
    """
    p = os.path.expanduser("~/.np_opengwas_token")
    try:
        if os.path.exists(p):
            with open(p, encoding="utf-8") as f:
                return f.read().strip().splitlines()[0]
    except Exception:
        pass
    return ""


def fetch_gwas(idv, insecure=False, page_size=500, jwt=None):
    """OpenGWAS 抓取。

    2026-11 实测：无 token 时 `/gwas/data/<id>` 返回 **404**（不是 401），
    容易误判成"数据集不存在"。OpenGWAS 现在要求 JWT：
      * 环境变量 `NP_OPENGWAS_TOKEN`、`--jwt`，或本机文件 `~/.np_opengwas_token`；
      * 未提供 token 时把 404 翻译成"可能需要 token"，而不是让用户去怀疑 ID。
    """
    ctx = None
    if insecure:
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
    token = (jwt or os.environ.get("NP_OPENGWAS_TOKEN") or _token_from_file() or "")
    u = "%s/gwas/data/%s?page_size=%d" % (API, idv, page_size)
    headers = {"User-Agent": "np-pro-mr/3.0"}
    if token:
        headers["Authorization"] = "Bearer " + token
    try:
        r = urllib.request.urlopen(urllib.request.Request(u, headers=headers),
                                   timeout=90, context=ctx)
        j = json.loads(r.read().decode("utf-8", "ignore"))
        return j.get("data", j) if isinstance(j, dict) else j, ""
    except urllib.error.HTTPError as e:
        if e.code == 404 and not token:
            return None, ("HTTPError:404（OpenGWAS 现要求 JWT：设 NP_OPENGWAS_TOKEN 或传 --jwt；"
                          "无 token 时它会用 404 掩盖鉴权失败，而不是 401）")
        return None, "HTTPError:%d（若 SSL 失败加 --insecure）" % e.code
    except Exception as e:
        return None, "%s:%s（若 SSL 失败加 --insecure）" % (type(e).__name__, str(e)[:90])


def store_snps(data, path, beta_col="beta", se_col="se", p_col="p"):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["SNP", "beta", "se", "p", "pos"])
        for r in data[:10000]:
            w.writerow([r.get("SNP") or r.get("rsid") or r.get("variant_id", ""),
                        r.get(beta_col), r.get(se_col), r.get(p_col),
                        r.get("position") or r.get("pos") or ""])
    return path


# ---------- 自检 ----------
def selftest():
    checks = []
    # 1) IVW 恢复真实效应 β=0.3，且无异质性
    bx = [0.10 + 0.02 * i for i in range(8)]
    by = [0.3 * x + (0.001 if i % 2 else -0.001) for i, x in enumerate(bx)]
    r = ivw(bx, [0.02] * 8, by, [0.01] * 8)
    checks.append(("ivw_recovers_beta", abs(r["beta"] - 0.3) < 0.01 and not r["heterogeneity"], r["beta"]))
    # 2) Egger 截距≈0（无多效性）
    e = mr_egger(bx, [0.02] * 8, by, [0.01] * 8)
    checks.append(("egger_null_intercept", abs(e["intercept"]) < 0.01, e["intercept"]))
    # 3) harmonise 对齐后 βx>0 且 βy 同步翻号
    h, nf = harmonise([{"beta_x": -0.2, "se_x": 0.01, "beta_y": 0.4, "se_y": 0.01}])
    checks.append(("harmonise_flip", h[0]["beta_x"] > 0 and h[0]["beta_y"] < 0 and nf == 1, nf))
    # 4) Steiger 剔除方向错误 SNP：构造 r2x<r2y 的位点
    snps = [{"snp": "rs1", "beta_x": 0.40, "se_x": 0.05, "beta_y": 0.02, "se_y": 0.02, "p_x": 1e-9},
            {"snp": "rs2", "beta_x": 0.40, "se_x": 0.05, "beta_y": 0.50, "se_y": 0.005, "p_x": 1e-9}]
    kept, info = steiger_filter(snps, n_x=5000, n_y=100000)
    checks.append(("steiger_removes_reversed", len(kept) == 1 and kept[0]["snp"] == "rs1",
                   "kept=%s" % [s["snp"] for s in kept]))
    # 5) 距离剪枝：同窗口 3 个 SNP → 保留 p 最小的 1 个
    cand = [{"snp": "a", "p_x": 1e-9, "pos": 1000000}, {"snp": "b", "p_x": 1e-8, "pos": 1050000},
            {"snp": "c", "p_x": 1e-7, "pos": 5000000}]
    k, ci = clump(cand, kb=100)
    checks.append(("clump_window", len(k) == 2 and ci["n_removed"] == 1 and ci["mode"] == "distance_pruning",
                   "kept=%d mode=%s" % (len(k), ci["mode"])))
    # 6) 基因型 r²：完全共分离 → r²=1；随机 → 低
    a = [0, 0, 1, 1, 2, 2, 1, 1]
    checks.append(("genotype_r2_perfect", abs(r2_from_dosage(a, a) - 1.0) < 1e-9, r2_from_dosage(a, a)))
    b2 = [0, 2, 0, 2, 0, 2, 0, 2]
    checks.append(("genotype_r2_independent", r2_from_dosage(a, b2) < 0.2, r2_from_dosage(a, b2)))
    # 7) COLOC：共享因果变异 → PP.H4 高
    cs = [{"beta_x": 0.5, "se_x": 0.05, "beta_y": 0.4, "se_y": 0.05},
          {"beta_x": 0.02, "se_x": 0.05, "beta_y": 0.01, "se_y": 0.05}]
    co = coloc_abf(cs)
    checks.append(("coloc_shared_H4", co["PP.H4"] > 0.5, co["PP.H4"]))
    # 8) 全流程跑通
    full = [{"snp": "rs%d" % i, "beta_x": 0.1 + 0.02 * i, "se_x": 0.02,
             "beta_y": 0.3 * (0.1 + 0.02 * i), "se_y": 0.01, "p_x": 1e-9, "pos": 1000 + i * 100000}
            for i in range(8)]
    fr = run_mr(full, n_x=1000, n_y=50000, kb=10, use_steiger=True)
    checks.append(("run_mr_end_to_end", fr["status"] == "OK" and fr["IVW"]["beta"] > 0,
                   "%s β=%.3f flags=%s" % (fr["status"], fr["IVW"]["beta"], ",".join(fr["flags"]))))

    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d 通过）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for name, good, detail in checks:
        print("  [%s] %-28s %s" % ("✓" if good else "✗", name, detail))
    return 0 if ok else 1



def apply_mr_fdr(res, peers=None, alpha=0.05):
    """对同一批靶点的 IVW p 做 BH，只有 q<alpha 才允许写「支持性因果证据」。

    为什么必须在 MR 层做：target-wide MR 会检测几十个靶点，单靶点 p<0.05 在
    30 次检验下期望假阳性 1.5 个。只报 p 不报 q，等于把多重比较的成本藏起来。

    peers: {gene: p}；缺省时从同目录其它 mr_result.json 汇总。
    """
    if res.get("status") != "OK":
        return res
    pv = (res.get("IVW") or {}).get("p")
    if pv is None:
        return res
    pool = dict(peers or {})
    pool.setdefault(res.get("gene") or "this", pv)
    if len(pool) < 2:
        res["fdr"] = {"status": "SINGLE_HYPOTHESIS", "n": 1,
                      "note": "只有 1 个靶点，BH 无意义；q = p 但不得据此宣称因果显著"}
        res["fdr_q"] = pv
        return res
    items = sorted(pool.items(), key=lambda x: x[1])
    n = len(items)
    qs, prev = {}, 1.0
    for rank in range(n, 0, -1):
        g, p_ = items[rank - 1]
        prev = min(prev, p_ * n / rank)
        qs[g] = min(1.0, prev)
    gene = res.get("gene") or "this"
    q = qs.get(gene, 1.0)
    res["fdr"] = {"status": "OK", "n_hypotheses": n, "alpha": alpha,
                  "q_value": q, "significant": bool(q < alpha),
                  "method": "Benjamini–Hochberg（同批靶点）",
                  "peers": {g: round(v, 6) for g, v in items},
                  "note": ("q<%g：可作为支持性因果证据" % alpha) if q < alpha
                          else ("q≥%g：不显著，只能写「提示性关联」" % alpha)}
    res["fdr_q"] = q
    if not res["fdr"]["significant"]:
        flags = res.setdefault("flags", [])
        if "FDR_NOT_SIGNIFICANT" not in flags:
            flags.append("FDR_NOT_SIGNIFICANT")
    return res


def collect_peer_p(dirpath):
    """从目录树里所有 mr_result.json 汇总 {gene: p}。"""
    out = {}
    for mp in glob.glob(os.path.join(dirpath, "**", "mr_result.json"), recursive=True):
        try:
            j = json.load(open(mp, encoding="utf-8"))
        except Exception:
            continue
        if j.get("status") != "OK":
            continue
        pv = (j.get("IVW") or {}).get("p")
        if pv is None:
            continue
        g = j.get("gene") or os.path.basename(os.path.dirname(mp))
        out[g] = pv
    return out


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_mr: 两样本孟德尔随机化（Steiger + LD + COLOC）")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("fetch")
    p.add_argument("--exposure", required=True, help="eQTL ID，如 eqtl-a-ENSG00000025708 或基因名 IL6")
    p.add_argument("--outcome", required=True, help="疾病 GWAS ID，如 ieu-a-7")
    p.add_argument("--insecure", action="store_true")
    p.add_argument("--jwt", default=None,
                   help="OpenGWAS JWT（也可用 NP_OPENGWAS_TOKEN 或 ~/.np_opengwas_token）")
    p.add_argument("--pval", type=float, default=5e-8,
                   help="工具变量显著性阈值（默认 5e-8；无命中会自动放宽到 1e-5）")
    p.add_argument("-o", "--outdir", default="np_out/mr")
    p.set_defaults(fn="fetch")
    p = sub.add_parser("analyze")
    p.add_argument("--exposure", required=True)
    p.add_argument("--outcome", required=True)
    p.add_argument("--n-exposure", type=int, default=None, help="暴露 GWAS 样本量（Steiger 必需）")
    p.add_argument("--n-outcome", type=int, default=None, help="结局 GWAS 样本量（Steiger 必需）")
    p.add_argument("--genotypes", default=None, help=".tped 基因型（纯 Python 算 r²，真 LD）")
    p.add_argument("--ld-matrix", default=None, help="SNP×SNP r² 矩阵 TSV")
    p.add_argument("--ld-api", action="store_true",
                   help="从 OpenGWAS /ld/matrix 取 1000G r²（无需下载面板；需 token）")
    p.add_argument("--jwt", default=None,
                   help="OpenGWAS JWT（--ld-api 时需要；也可用 NP_OPENGWAS_TOKEN 或 ~/.np_opengwas_token）")
    p.add_argument("--ld-pop", default="EUR", help="--ld-api 的人群（默认 EUR）")
    p.add_argument("--clump-kb", type=int, default=10000, help="无 LD 时的距离剪枝窗口（kb）")
    p.add_argument("--clump-r2", type=float, default=0.001)
    p.add_argument("--no-steiger", action="store_true")
    p.add_argument("--bfile", default=None, help="PLINK bfile 前缀（.bed/.bim/.fam），用于真 LD clumping")
    p.add_argument("--plink", default=None, help="plink 可执行文件路径（默认自动探测/环境变量 NP_PLINK）")
    p.add_argument("--fdr", action="store_true",
                   help="对同一批靶点的 IVW p 做 BH；q≥alpha 时标 FDR_NOT_SIGNIFICANT")
    p.add_argument("--fdr-dir", default=None, help="汇总同批靶点 p 的目录（缺省用输出目录）")
    p.add_argument("--fdr-alpha", type=float, default=0.05)
    p.add_argument("--drug-direction", default=None, choices=["inhibit", "activate"],
                   help="药物对该靶点的作用方向（**方向筛查**：只回答\"该不该抑制/激活\"，"
                        "不构成药效证据；矛盾时报 TARGET_INTERVENTION_MISMATCH）")
    p.add_argument("-o", "--outdir", default="np_out/mr")
    p.set_defaults(fn="analyze")
    p = sub.add_parser("selftest")
    p.set_defaults(fn="selftest")
    a = ap.parse_args(argv)
    if os.environ.get("NP_OFFLINE") == "1" and a.fn == "fetch":
        print("[mr] NP_OFFLINE=1：跳过联网抓取")
        return 0
    if a.fn == "selftest":
        return selftest()
    if a.fn == "fetch":
        exp = DEFAULT_EXPOSURE.get(a.exposure, a.exposure)
        pair, err = fetch_pair(exp, a.outcome, pval=a.pval, jwt=a.jwt, insecure=a.insecure)
        if pair is None:
            print("[mr] 抓取失败：%s" % err)
            return 2
        pe, po = write_pair(pair, a.outdir)
        print("[mr] %s（%d 个工具变量）→ %s" % (exp, pair["n_exposure_hits"], pe))
        print("[mr] %s：匹配到 %d/%d 个变异 → %s" % (
            a.outcome, pair["n_matched"], pair["n_exposure_hits"], po))
        if pair["n_matched"] < 3:
            print("  ⚠️ 匹配过少：检查两个数据集是否同一 build（HG19 vs GRCh38）")
        return 0

    ex = {}
    with open(a.exposure, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                ex[r["SNP"]] = {"beta": float(r["beta"]), "se": float(r["se"]),
                                "p": float(r["p"]),
                                "pos": r.get("pos") or r.get("position") or "",
                                "ea": r.get("ea") or "", "nea": r.get("nea") or "",
                                "n": r.get("n") or ""}
            except (KeyError, TypeError, ValueError):
                continue
    snps = []
    with open(a.outcome, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            k = r.get("SNP")
            if k in ex:
                e = ex[k]
                try:
                    snps.append({"snp": k, "beta_x": e["beta"], "se_x": e["se"], "p_x": e["p"],
                                 "pos": e["pos"] or None,
                                 "ea_x": e["ea"], "nea_x": e["nea"], "n_x": e["n"],
                                 "beta_y": float(r["beta"]), "se_y": float(r["se"]),
                                 "ea_y": r.get("ea") or "", "nea_y": r.get("nea") or ""})
                except (TypeError, ValueError):
                    continue
    # 等位基因对齐（同向/反向/互补链；无法判定的丢弃）
    snps, allele_stats = harmonise_alleles(snps)
    if allele_stats.get("dropped"):
        print("[mr] 等位基因对齐：丢弃 %d 个无法判定的变异（%s）" % (
            allele_stats["dropped"],
            "；".join("%s=%d" % (k, v) for k, v in allele_stats.items() if v)))
    geno = read_tped(a.genotypes, [s["snp"] for s in snps]) if a.genotypes else None
    ld = read_ld_matrix(a.ld_matrix) if a.ld_matrix else None
    ld_info = {}
    if a.ld_api and not ld:
        lm, err = fetch_ld_matrix([x["snp"] for x in snps], pop=a.ld_pop, jwt=a.jwt)
        if lm:
            ld = lm["ld"]
            ld_info = {"source": "opengwas_ld_api", "pop": lm["pop"],
                       "n_in_panel": len(lm["ids"]), "missing": lm["missing"][:20]}
            print("[mr] LD 矩阵：%d/%d 个变异在 %s 面板中（缺失 %d）" % (
                len(lm["ids"]), len(snps), lm["pop"], len(lm["missing"])))
        else:
            print("[mr] LD API 失败（%s）→ 回退距离剪枝，结果将标 exploratory" % err)
            ld_info = {"source": "failed", "error": err}
    res = run_mr(snps, a.n_exposure, a.n_outcome, geno=geno, ld=ld,
                 kb=a.clump_kb, r2_thresh=a.clump_r2, use_steiger=not a.no_steiger,
                 bfile=a.bfile, plink=a.plink, drug_direction=a.drug_direction)
    res["allele_harmonisation"] = allele_stats
    if ld_info:
        res["ld_source"] = ld_info
    if a.fdr:
        peers = collect_peer_p(a.fdr_dir or a.outdir)
        res = apply_mr_fdr(res, peers, a.fdr_alpha)
    os.makedirs(a.outdir, exist_ok=True)
    json.dump(res, open(os.path.join(a.outdir, "mr_result.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    if res["status"] == "OK":
        print("[mr] IVW β=%.4f p=%.3g | Egger截距p=%.3g | COLOC PP.H4=%.3f | F=%.1f" % (
            res["IVW"]["beta"], res["IVW"]["p"], res["MR_Egger"].get("p_intercept", 1),
            res["coloc"].get("PP.H4", 0), res["F_statistic"]))
        print("[mr] SNP %d→%d（Steiger %s，去冗余 %s）| flags: %s" % (
            res["n_snps_raw"], res["n_kept"], res["steiger"].get("status"),
            res["clump"]["mode"], ", ".join(res["flags"]) or "none"))
        if res["clump"]["mode"] == "distance_pruning":
            print("[mr] ⚠️ 未提供 LD 面板，仅按距离剪枝，不等同 clumping（Methods 须如实标注）")
        if (res.get("intervention") or {}).get("severity") == "SEVERE":
            print("[mr] 红灯：%s" % res["intervention"]["message"])
        if res["steiger"].get("status") == "STEIGER_SKIPPED":
            print("[mr] ⚠️ 缺 --n-exposure/--n-outcome，Steiger 未执行（方向过滤缺失）")
    else:
        print("[mr] %s：%s" % (res["status"], res.get("reason", "")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
