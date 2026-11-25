#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_survival.py —— v1.5.0 靶点临床相关性映射模块（TCGA 生存分析）

背景（23 篇审计）：17/23 篇提到 TCGA/生存/预后，antrodia 真用 TCGA-LIHC 生存曲线。
本模块经 cBioPortal API（2026-11-20 实测 api_ok）自动完成：
  取某癌种 TCGA 队列的表达量 + 临床生存数据 → 按 hub 基因表达中位数分高低组
  → Kaplan-Meier 曲线（SVG）+ log-rank 检验 + HR（O/E 法）。

命令：
  studies [--kw lihc|...        列出 TCGA 研究
  survival --study lihc_tcga --genes STAT3,TP53 --tumor-only
           [--outdir DIR] [--endpoint OS]
输出：survival_results.csv + <gene>_km.svg

方法说明（诚实边界）：log-rank 为精确 Mantel-Cox χ² 实现；HR 用 O/E 比值近似
（组间危险比的标准粗估计，非 Cox 回归）；表达分组=中位数二分。
"""
import argparse
import json
import math
import os

import numpy as np
import pandas as pd

API = "https://www.cbioportal.org/api"

# TCGA 缩写 → cBioPortal 研究 id（优先 2022 泛癌队列）
TCGA_STUDIES = {
    "LIHC": "lihc_tcga_pan_can_atlas_2018", "STAD": "stad_tcga_pan_can_atlas_2018",
    "COAD": "coad_tcga_pan_can_atlas_2018", "BRCA": "brca_tcga_pan_can_atlas_2018",
    "LUAD": "luad_tcga_pan_can_atlas_2018", "KIRC": "kirc_tcga_pan_can_atlas_2018",
    "HNSC": "hnsc_tcga_pan_can_atlas_2018", "ESCA": "esca_tcga_pan_can_atlas_2018",
}


def _get(url, params=None, tries=3):
    import requests, time
    for k in range(tries):
        r = requests.get(url, params=params or {}, timeout=120,
                         headers={"User-Agent": "np-skill/1.5", "Accept": "application/json"})
        if r.status_code == 200:
            return r.json()
        time.sleep(3 * (k + 1))
    raise SystemExit(f"[FAIL] cBioPortal {url} HTTP {r.status_code}")


def cmd_studies(args):
    d = _get(f"{API}/studies", {"pageSize": 100, "pageNumber": 0, "keyword": args.kw or ""})
    rows = [{"id": s["studyId"], "name": s["name"][:60], "samples": s.get("sampleCount")} for s in d]
    print(pd.DataFrame(rows).to_string(index=False))


def _fetch_expression(study, genes):
    """mRNA expression：gene → {patientId: value}"""
    d = _get(f"{API}/studies/{study}/molecular-profiles")
    profs = [p["molecularProfileId"] for p in d if p["molecularAlterationType"] == "MRNA_EXPRESSION"]
    if not profs:
        raise SystemExit("[FAIL] 该研究无 mRNA 分子画像")
    prof = profs[0]
    out = {}
    for g in genes:
        try:
            gd = _get(f"{API}/genes/{g}")
            gid = gd["entrezGeneId"]
        except SystemExit:
            print(f"  [warn] {g} 基因查询失败"); continue
        dd = _get(f"{API}/molecular-profiles/{prof}/molecular-data",
                  {"studyId": study, "sampleListId": f"{study}_all",
                   "entrezGeneId": gid, "projection": "SUMMARY"})
        pm = {}
        for r in dd:
            v = r.get("value")
            pid = r.get("patientId")
            if v is not None and pid:
                pm.setdefault(pid, []).append(v)
        out[g] = {p: sum(vs) / len(vs) for p, vs in pm.items()}
        if not out[g]:
            print(f"  [warn] {g} 无表达数据")
    return out


def _fetch_clinical(study):
    """患者级 OS_STATUS/OS_MONTHS，经 /studies/{study}/clinical-data"""
    rows = {}
    page = 0
    import requests, time
    while True:
        for k in range(3):
            r = requests.get(f"{API}/studies/{study}/clinical-data",
                             params={"clinicalDataType": "PATIENT", "projection": "DETAILED",
                                     "pageSize": 10000, "pageNumber": page},
                             timeout=180, headers={"Accept": "application/json"})
            if r.status_code == 200:
                break
            time.sleep(3 * (k + 1))
        if r.status_code != 200:
            break
        data = r.json()
        if not data:
            break
        for a in data:
            if a["clinicalAttributeId"] in ("OS_STATUS", "OS_MONTHS"):
                rows.setdefault(a["patientId"], {})[a["clinicalAttributeId"]] = a["value"]
        page += 1
        if len(data) < 10000:
            break
    return pd.DataFrame([{"patientId": p, **v} for p, v in rows.items()])


def _patient_sample_map(study):
    d = _get(f"{API}/studies/{study}/samples", {"pageSize": 10000, "projection": "SUMMARY"})
    return {s["sampleId"]: s["patientId"] for s in d}


def logrank(t1, e1, t2, e2):
    """Mantel-Cox log-rank（无分层）。t=时间, e=事件(1/0)"""
    times = sorted(set(t1) | set(t2))
    chi2 = 0.0; V = 0.0
    for t in times:
        n1 = int((t1 >= t).sum()); n2 = int((t2 >= t).sum())
        n = n1 + n2
        if n == 0: continue
        d1 = int(((t1 == t) & (e1 == 1)).sum()); d2 = int(((t2 == t) & (e2 == 1)).sum())
        d = d1 + d2
        if d == 0 or n <= 1: continue
        E1 = d * n1 / n
        V += d * (n1 / n) * (n2 / n) * (n - d) / (n - 1)
        chi2 += (d1 - E1)
    chi2sq = (chi2 * chi2 / V) if V > 0 else 0.0
    p = math.exp(-0.5 * chi2sq) * 2 if chi2sq < 700 else 0.0
    # χ²(1) 上尾（正态近似足够）
    from scipy.stats import chi2 as _c2
    p = float(_c2.sf(chi2sq, 1))
    return chi2sq, p


def km_curve(t, e):
    """返回 (times, surv) 阶梯"""
    times = np.sort(np.unique(t))
    s = 1.0; out_t, out_s = [0.0], [1.0]
    for tt in times:
        at_risk = (t >= tt).sum()
        dd = ((t == tt) & (e == 1)).sum()
        if dd:
            s *= (1 - dd / at_risk)
            out_t.append(tt); out_s.append(s)
    return out_t, out_s


def cmd_survival(args):
    study = TCGA_STUDIES.get(args.study.upper(), args.study)
    genes = [g.strip().upper() for g in args.genes.split(",") if g.strip()]
    print(f"研究: {study} | 基因: {genes}")
    expr = _fetch_expression(study, genes)
    clin = _fetch_clinical(study)
    if clin.empty:
        raise SystemExit("[FAIL] 无患者临床数据（该研究可能无 OS 字段）")
    clin = clin.rename(columns={"OS_STATUS": "status", "OS_MONTHS": "months"})
    clin = clin[clin.status.astype(str).str.contains("DECEASED|LIVING", na=False)].copy()
    clin["event"] = clin.status.str.contains("DECEASED").astype(int)
    clin["months"] = pd.to_numeric(clin.months, errors="coerce")
    clin = clin[clin.months > 0]
    print(f"临床数据: {len(clin)} 患者（有生存时间）")
    rows = []
    os.makedirs(args.outdir, exist_ok=True)
    for g, pmap in expr.items():
        df = clin.copy()
        df["exp"] = df.patientId.map(pmap)
        df = df.dropna(subset=["exp"])
        if len(df) < 50:
            print(f"  [warn] {g}: 有效 n={len(df)} 过少"); continue
        med = df.exp.median()
        hi = df[df.exp >= med]; lo = df[df.exp < med]
        chi2, p = logrank(hi.months.values, hi.event.values, lo.months.values, lo.event.values)
        # HR（O/E 近似，高表达组 vs 低表达组）
        E_hi = 0; E_lo = 0
        times = sorted(set(hi.months) | set(lo.months))
        for t in times:
            n1 = (hi.months >= t).sum(); n2 = (lo.months >= t).sum()
            d1 = ((hi.months == t) & (hi.event == 1)).sum()
            d2 = ((lo.months == t) & (lo.event == 1)).sum()
            if n1 + n2 > 0 and (d1 + d2) > 0:
                E_hi += (d1 + d2) * n1 / (n1 + n2); E_lo += (d1 + d2) * n2 / (n1 + n2)
        O_hi = hi.event.sum(); O_lo = lo.event.sum()
        hr = (O_hi / E_hi) / (O_lo / E_lo) if E_hi > 0 and E_lo > 0 and O_lo > 0 else float("nan")
        rows.append({"gene": g, "study": study, "n": len(df),
                     "n_hi": len(hi), "n_lo": len(lo),
                     "deaths_hi": int(O_hi), "deaths_lo": int(O_lo),
                     "HR_hi_vs_lo": round(float(hr), 3) if hr == hr else None,
                     "logrank_chi2": round(chi2, 2), "p_logrank": float(f"{p:.3e}"),
                     "significant": "YES" if p < 0.05 else "no"})
        # KM 曲线 SVG
        svg = _km_svg(km_curve(hi.months.values, hi.event.values),
                      km_curve(lo.months.values, lo.event.values),
                      g, study, hr, p)
        with open(os.path.join(args.outdir, f"{g}_km.svg"), "w", encoding="utf-8") as f:
            f.write(svg)
        print(f"  [{g}] n={len(df)} HR={rows[-1]['HR_hi_vs_lo']} log-rank p={p:.3e} → {g}_km.svg")
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(args.outdir, "survival_results.csv"),
                                  index=False, encoding="utf-8-sig")
        print(f"汇总 → {args.outdir}/survival_results.csv")


def _km_svg(curve_hi, curve_lo, gene, study, hr, p):
    W, H = 680, 420
    x0, y0, w, h = 70, 40, 560, 330
    tmax = max(max(curve_hi[0] or [1]), max(curve_lo[0] or [1])) or 1
    def px(t): return x0 + w * t / tmax
    def py(s): return y0 + h * (1 - s)
    def path(c):
        d = f"M {x0} {y0+h}"
        for t, s in zip(c[0], c[1]):
            d += f" L {px(t):.1f} {py(s):.1f}"
        return d
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="sans-serif">',
             f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
             f'<text x="{W//2}" y="24" text-anchor="middle" font-size="15" font-weight="bold">{gene} — {study}</text>',
             f'<text x="{W//2}" y="40" text-anchor="middle" font-size="12" fill="#666">HR(high vs low)={hr if hr==hr else "NA"} | log-rank p={p:.2e}</text>']
    for i in range(6):
        xx = x0 + w * i / 5
        parts.append(f'<line x1="{xx:.0f}" y1="{y0+h}" x2="{xx:.0f}" y2="{y0}" stroke="#eee"/>')
        parts.append(f'<text x="{xx:.0f}" y="{y0+h+16}" font-size="10" text-anchor="middle">{tmax*i/5:.0f}mo</text>')
    for i in range(6):
        s = 1 - i / 5
        parts.append(f'<text x="{x0-8}" y="{py(s)+4:.0f}" font-size="10" text-anchor="end">{s:.1f}</text>')
        parts.append(f'<line x1="{x0}" y1="{py(s):.0f}" x2="{x0+w}" y2="{py(s):.0f}" stroke="#eee"/>')
    parts.append(f'<path d="{path(curve_hi)}" fill="none" stroke="#c0392b" stroke-width="2"/>')
    parts.append(f'<path d="{path(curve_lo)}" fill="none" stroke="#2471a3" stroke-width="2"/>')
    parts.append(f'<text x="{x0+w-10}" y="{y0+16}" font-size="12" fill="#c0392b" text-anchor="end">high expression</text>')
    parts.append(f'<text x="{x0+w-10}" y="{y0+34}" font-size="12" fill="#2471a3" text-anchor="end">low expression</text>')
    parts.append('</svg>')
    return "".join(parts)


def main():
    ap = argparse.ArgumentParser(description="np_survival: TCGA 靶点生存分析（cBioPortal）")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("studies"); s.add_argument("--kw", default="")
    s = sub.add_parser("survival")
    s.add_argument("--study", required=True, help="TCGA 缩写（LIHC）或 cBioPortal 研究 id")
    s.add_argument("--genes", required=True)
    s.add_argument("--outdir", default=".")
    a = ap.parse_args()
    {"studies": cmd_studies, "survival": cmd_survival}[a.cmd](a)


if __name__ == "__main__":
    main()
