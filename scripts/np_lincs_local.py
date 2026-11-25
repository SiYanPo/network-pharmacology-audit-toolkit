#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_lincs_local.py —— LINCS L1000 本地化流式读取器（彻底摆脱外部 API）

为什么必须本地化
    上一版只能生成"CMap 查询模板"，审稿人一句"那你到底能不能查？"就废了。
    CLUE 官方 API 需授权且经常 502；GEO 上的 Level 5 矩阵是**免费公开**的，
    格式为 .gctx（本质是 HDF5），完全可以在本机做签名比对。

数据从哪来（都免费，无需注册）
    ① GSE92742（LINCS L1000 Phase I, 2017）
       https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
       关键文件：GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx（约 10 GB）
       配套：GSE92742_Broad_LINCS_gene_info.txt.gz（基因注释）、..._sig_info.txt.gz（扰动注释）
    ② GSE70138（Phase II, 2017）
       GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx
    ③ 只要 978 landmark 基因时，可下 Level 3 的 978 维矩阵，体积小得多（约 1–2 GB）

gctx 里有什么
    /data                     float32，12,328 基因 × N 个扰动实例（chunked，可切片读）
    /meta/row_meta            基因维注释（pr_gene_symbol 等）
    /meta/col_meta            实例维注释（pert_iname / cell_id / pert_idose / pert_time）
    /meta/col_meta_group_size 每个扰动聚合了多少个重复

算法（CMap 风格，但用可复现的显式公式）
    1) 只读签名基因对应的行（不载入整个 10 GB 矩阵）；
    2) 对每个基因在该行内做 z-score（跨实例标准化）→ 消除基因基线差异；
    3) 实例得分 S = mean(z_up) − mean(z_dn)；S 越负 = 越能逆转疾病签名 = 候选；
    4) 同一 perturbagen 的多个实例取"最强 3 个的均值"，并记录实例数；
    5) 输出候选表 + manifest（含数据版本、文件大小、sha256 前 16 位）。

诚实边界
    * 本实现不做 CMap 官方的 tau/plate-normalization（需 plate 级信息），
      数值与 CLUE 的 connectivity score 不可直接等价；
    * Level 5 的 MODZ 已是 moderated z-score，再 z-score 一次属"二次标准化"，
      只影响相对排序，不改变符号，Methods 需说明；
    * 候选只说明"表达层面可逆转"，与结合活性无关。

用法
    python np_lincs_local.py probe --gctx GSE92742_....gctx
    python np_lincs_local.py query --gctx X.gctx --signature sig.json -o out/ [--top 30]
    python np_lincs_local.py query --gctx X.gctx --up up.txt --down down.txt --cell A375
    python np_lincs_local.py selftest          # 自动构造小型 HDF5 验证读取与打分
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys

try:
    import numpy as np
except ImportError:  # 允许无 numpy 时打印说明而不是崩栈
    np = None

HERE = os.path.dirname(os.path.abspath(__file__))


def _h5():
    try:
        import h5py
        return h5py
    except ImportError:
        return None


def _decode(arr):
    """h5py 的 byte string 数组 → python str 列表。"""
    out = []
    for x in arr:
        if isinstance(x, (bytes, bytearray)):
            out.append(x.decode("utf-8", "replace"))
        else:
            out.append(str(x))
    return out


def _meta_col(f, group, name):
    """安全读取 meta 子表的一列，不存在返回 None。"""
    g = f.get("/meta/%s" % group)
    if g is None:
        return None
    if name in g:
        return _decode(g[name][:])
    # 部分版本把 meta 存成 compound dataset 的一级字段
    if isinstance(g, h5py_dataset := g):
        try:
            return _decode(g[name][:])
        except Exception:
            return None
    return None


def probe(path):
    """打印 gctx 结构（基因数/实例数/可用注释），不做任何计算。"""
    h5py = _h5()
    if h5py is None:
        return {"status": "H5PY_MISSING", "hint": "pip install h5py"}
    if not os.path.exists(path):
        return {"status": "FILE_MISSING", "path": path}
    info = {"status": "OK", "path": os.path.abspath(path),
            "size_gb": round(os.path.getsize(path) / 1024 ** 3, 2)}
    with h5py.File(path, "r") as f:
        info["keys"] = list(f.keys())
        d = f.get("/data")
        if d is None:
            info["status"] = "NO_DATA"
            return info
        info["data_shape"] = list(d.shape)
        info["data_dtype"] = str(d.dtype)
        info["chunks"] = list(d.chunks) if d.chunks else None
        for grp in ("row_meta", "col_meta"):
            g = f.get("/meta/%s" % grp)
            if g is not None:
                try:
                    info[grp + "_fields"] = list(g.dtype.names or [])
                except Exception:
                    info[grp + "_fields"] = []
        rm = _meta_col(f, "row_meta", "pr_gene_symbol")
        cm = _meta_col(f, "col_meta", "pert_iname")
        info["n_row_meta"] = len(rm) if rm else None
        info["n_col_meta"] = len(cm) if cm else None
        info["example_genes"] = (rm[:5] if rm else None)
        info["example_perts"] = (cm[:5] if cm else None)
    return info


def _load_signature(sig_path=None, up=None, down=None):
    upg, dng = [], []
    if sig_path:
        j = json.load(open(sig_path, encoding="utf-8"))
        upg = [g.strip().upper() for g in (j.get("up_genes") or [])]
        dng = [g.strip().upper() for g in (j.get("down_genes") or [])]
    if up:
        upg += [l.strip().upper() for l in open(up, encoding="utf-8-sig") if l.strip()]
    if down:
        dng += [l.strip().upper() for l in open(down, encoding="utf-8-sig") if l.strip()]
    seen = set()
    upg = [g for g in upg if g and not (g in seen or seen.add(g))]
    seen = set()
    dng = [g for g in dng if g and not (g in seen or seen.add(g))]
    return upg, dng


def query(path, up_genes, down_genes, gene_field="pr_gene_symbol", pert_field="pert_iname",
          cell=None, top=30, min_overlap=0.3, batch=4000):
    """流式读取 gctx，对签名打分。返回 dict（候选表 + 统计 + 诚实标记）。"""
    h5py = _h5()
    if h5py is None:
        return {"status": "H5PY_MISSING", "hint": "pip install h5py"}
    if np is None:
        return {"status": "NUMPY_MISSING", "hint": "pip install numpy"}
    if not os.path.exists(path):
        return {"status": "FILE_MISSING", "path": path}
    up = [g.upper() for g in up_genes]
    dn = [g.upper() for g in down_genes]
    with h5py.File(path, "r") as f:
        d = f.get("/data")
        if d is None:
            return {"status": "NO_DATA"}
        symbols = _meta_col(f, "row_meta", gene_field) or _meta_col(f, "row_meta", "id")
        perts = _meta_col(f, "col_meta", pert_field) or _meta_col(f, "col_meta", "id")
        cells = _meta_col(f, "col_meta", "cell_id")
        doses = _meta_col(f, "col_meta", "pert_idose")
        times = _meta_col(f, "col_meta", "pert_time")
        if not symbols or not perts:
            return {"status": "NO_META", "reason": "缺 row_meta/col_meta 注释"}
        sym2i = {}
        for i, s in enumerate(symbols):
            sym2i.setdefault(s.upper(), i)
        n_row, n_col = len(symbols), len(perts)
        if list(d.shape) == [n_row, n_col]:
            genes_axis = 0
        elif list(d.shape) == [n_col, n_row]:
            genes_axis = 1
        else:
            return {"status": "SHAPE_MISMATCH",
                    "reason": "data.shape=%s 与 meta(%d×%d) 不匹配" % (list(d.shape), n_row, n_col)}

        hit_up = [sym2i[g] for g in up if g in sym2i]
        hit_dn = [sym2i[g] for g in dn if g in sym2i]
        cov_up = len(hit_up) / max(1, len(up))
        cov_dn = len(hit_dn) / max(1, len(dn))
        if not hit_up and not hit_dn:
            return {"status": "NO_OVERLAP",
                    "reason": "签名基因全部不在该 gctx 的注释里（检查基因符号版本：Entrez/Ensembl/HGNC）"}
        if min(cov_up, cov_dn) < min_overlap:
            note = ("签名覆盖率低（up %.0f%% / down %.0f%%）：结果仅作探索，"
                    "建议用 gene_info 把符号映射到 gctx 的 pr_gene_symbol 版本" %
                    (100 * cov_up, 100 * cov_dn))
        else:
            note = ""

        inst_ids = np.arange(n_col)
        if cell and cells:
            mask = np.array([c == cell for c in cells])
            if mask.sum() == 0:
                return {"status": "NO_CELL", "reason": "cell_id=%s 不存在" % cell}
            inst_ids = inst_ids[mask]
        else:
            mask = None

        n_up = len(hit_up)
        n_dn = len(hit_dn)
        # 逐行读取（每行 n_col 个 float，约 2–10 MB，安全）
        rows = []
        for gi in hit_up + hit_dn:
            if genes_axis == 0:
                v = d[gi, :].astype(np.float32)
            else:
                v = d[:, gi].astype(np.float32)
            if mask is not None:
                v = v[mask]
            mu = float(np.nanmean(v))
            sd = float(np.nanstd(v)) or 1.0
            rows.append((v - mu) / sd)
        if not rows:
            return {"status": "NO_OVERLAP"}
        Z = np.vstack(rows)                       # 基因 × 实例
        score = np.zeros(Z.shape[1], dtype=np.float32)
        if n_up:
            score += Z[:n_up].mean(axis=0)
        if n_dn:
            score -= Z[n_up:].mean(axis=0)
        if n_up and n_dn:
            score /= 2.0

    # 按扰动聚合：同一 perturbagen 取最强（最负）3 个实例的均值
    agg = {}
    for k, ii in enumerate(inst_ids):
        p = perts[ii] if ii < len(perts) else "?"
        s = float(score[k])
        agg.setdefault(p, []).append((s, int(ii)))
    rows_out = []
    for p, vals in agg.items():
        vals.sort(key=lambda x: x[0])
        k3 = vals[:3]
        rows_out.append({"pert_iname": p, "score": round(sum(v[0] for v in k3) / len(k3), 4),
                         "n_instances": len(vals),
                         "cell": cells[vals[0][1]] if cells and vals[0][1] < len(cells) else "",
                         "dose": doses[vals[0][1]] if doses and vals[0][1] < len(doses) else "",
                         "time": times[vals[0][1]] if times and vals[0][1] < len(times) else ""})
    rows_out.sort(key=lambda r: r["score"])
    reversal = [r for r in rows_out if r["score"] < -0.2]
    return {"status": "OK", "n_genes_in_matrix": n_row, "n_instances": n_col,
            "n_instances_used": int(len(inst_ids)),
            "signature": {"n_up": len(up), "n_down": len(dn), "hit_up": len(hit_up),
                          "hit_down": len(hit_dn), "coverage_up": round(cov_up, 3),
                          "coverage_down": round(cov_dn, 3)},
            "n_perturbagens": len(rows_out), "n_reversal": len(reversal),
            "top_reversal": rows_out[:top],
            "mode": "local_gctx_streaming", "score_definition": "mean(z_up)-mean(z_dn)，z 在实例维标准化",
            "note": note or "负分=逆转疾病签名（候选）；正分=模拟疾病签名（应避免）",
            "caveat": "非 CMap 官方 connectivity score（无 tau/plate 归一化），数值不可与 CLUE 直接比较"}


def selftest():
    """构造小型 HDF5（20 基因 × 12 实例）验证：读取方向、z-score、逆转排序。"""
    h5py = _h5()
    if h5py is None or np is None:
        print("selftest: SKIP（需要 h5py + numpy）")
        return 0
    import tempfile
    tmp = os.path.join(tempfile.gettempdir(), "np_lincs_selftest.gctx")
    genes = ["G%02d" % i for i in range(20)]
    perts = ["drugA", "drugA", "drugA", "drugB", "drugB", "drugB",
             "drugC", "drugC", "drugC", "drugD", "drugD", "drugD"]
    # 构造：drugA 上调 up 基因/下调 down 基因 → 模拟疾病（正分）；
    #       drugB 完全逆转 → 负分；drugC 无效应；drugD 弱逆转
    base = np.zeros((20, 12), dtype=np.float32)
    rng = np.random.default_rng(7)
    base += rng.normal(0, 0.1, base.shape)
    up_idx, dn_idx = [0, 1, 2], [3, 4, 5]
    for j in (0, 1, 2):                       # drugA: 模拟
        base[up_idx, j] += 1.5
        base[dn_idx, j] -= 1.5
    for j in (3, 4, 5):                       # drugB: 逆转
        base[up_idx, j] -= 1.5
        base[dn_idx, j] += 1.5
    for j in (9, 10, 11):                     # drugD: 弱逆转
        base[up_idx, j] -= 0.6
        base[dn_idx, j] += 0.6
    with h5py.File(tmp, "w") as f:
        f.create_dataset("/data", data=base, chunks=(5, 12))
        rm = f.create_group("/meta/row_meta")
        rm.create_dataset("pr_gene_symbol", data=np.array([g.encode() for g in genes]))
        rm.create_dataset("id", data=np.array([g.encode() for g in genes]))
        cm = f.create_group("/meta/col_meta")
        cm.create_dataset("pert_iname", data=np.array([p.encode() for p in perts]))
        cm.create_dataset("cell_id", data=np.array([b"A375"] * 12))
        cm.create_dataset("pert_idose", data=np.array([b"10 uM"] * 12))
        cm.create_dataset("pert_time", data=np.array([b"6"] * 12))
    pr = probe(tmp)
    r = query(tmp, ["G00", "G01", "G02"], ["G03", "G04", "G05"], top=5)
    order = [x["pert_iname"] for x in r.get("top_reversal", [])]
    ok = (pr.get("status") == "OK" and r.get("status") == "OK"
          and order[:2] == ["drugB", "drugD"]
          and r["top_reversal"][0]["score"] < -0.5)
    print("selftest: %s（shape=%s，top=%s，score=%.2f）" % (
        "PASS" if ok else "FAIL", pr.get("data_shape"), order[:3],
        r["top_reversal"][0]["score"] if r.get("top_reversal") else 0))
    try:
        os.remove(tmp)
    except OSError:
        pass
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--selftest"]:      # 兼容统一写法：--selftest 与 selftest 等价
        argv = ["selftest"]
    ap = argparse.ArgumentParser(description="np_lincs_local: LINCS L1000 gctx 本地流式比对")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("probe"); p.add_argument("--gctx", required=True)
    p.set_defaults(fn="probe")
    q = sub.add_parser("query")
    q.add_argument("--gctx", required=True)
    q.add_argument("--signature", default=None, help="np_lincs signature 产出的 sig.json")
    q.add_argument("--up", default=None); q.add_argument("--down", default=None)
    q.add_argument("--cell", default=None, help="限定 cell_id（如 A375 / MCF7）")
    q.add_argument("--gene-field", default="pr_gene_symbol")
    q.add_argument("--pert-field", default="pert_iname")
    q.add_argument("--top", type=int, default=30)
    q.add_argument("--min-overlap", type=float, default=0.3)
    q.add_argument("-o", "--outdir", default="np_out/lincs_local")
    q.set_defaults(fn="query")
    s = sub.add_parser("selftest"); s.set_defaults(fn="selftest")
    a = ap.parse_args(argv)
    if a.fn == "selftest":
        return selftest()
    if a.fn == "probe":
        info = probe(a.gctx)
        print(json.dumps(info, ensure_ascii=False, indent=1))
        return 0 if info.get("status") == "OK" else 2

    upg, dng = _load_signature(a.signature, a.up, a.down)
    if not upg and not dng:
        raise SystemExit("[FAIL] 需 --signature 或 --up/--down")
    res = query(a.gctx, upg, dng, a.gene_field, a.pert_field, a.cell, a.top, a.min_overlap)
    os.makedirs(a.outdir, exist_ok=True)
    res["source"] = {"gctx": os.path.abspath(a.gctx)}
    if os.path.exists(a.gctx):
        with open(a.gctx, "rb") as f:
            res["source"]["sha256_16"] = hashlib.sha256(f.read(1024 * 1024)).hexdigest()[:16]
        res["source"]["size_gb"] = round(os.path.getsize(a.gctx) / 1024 ** 3, 2)
    json.dump(res, open(os.path.join(a.outdir, "lincs_local_result.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    if res.get("status") == "OK":
        with open(os.path.join(a.outdir, "lincs_candidates.csv"), "w",
                  encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["pert_iname", "score", "n_instances", "cell", "dose", "time"])
            w.writeheader()
            w.writerows(res["top_reversal"])
        print("[lincs_local] %d 基因 × %d 实例；签名命中 up %d/%d、down %d/%d；"
              "%d 个扰动，%d 个逆转候选 → %s" % (
                  res["n_genes_in_matrix"], res["n_instances_used"],
                  res["signature"]["hit_up"], res["signature"]["n_up"],
                  res["signature"]["hit_down"], res["signature"]["n_down"],
                  res["n_perturbagens"], res["n_reversal"], a.outdir))
        for r in res["top_reversal"][:5]:
            print("  %-22s score=%+.3f  n=%d  %s" % (r["pert_iname"], r["score"], r["n_instances"], r["cell"]))
        if res.get("note"):
            print("  ⚠️ " + res["note"])
    else:
        print("[lincs_local] %s：%s" % (res.get("status"), res.get("reason", res.get("hint", ""))))
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
