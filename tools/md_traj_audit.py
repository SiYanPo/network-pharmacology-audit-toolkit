# -*- coding: utf-8 -*-
"""md_traj_audit.py —— 真实多帧复合物 MD 轨迹的稳定性与统计功效复核。

用途
----
把一条显式溶剂生产轨迹（CHARMM DCD，OpenMM 写出）做成可直接进论文的
"轨迹稳定性 + 结合自由能统计"证据，且**不依赖 mdtraj/MDAnalysis**
（纯 numpy + struct，任何装了 numpy 的环境都能跑）。

为什么需要它
------------
MM-GBSA/MM-PBSA 报 "mean ± SEM" 时，若把每帧当独立样本，SEM 会被系统性低估
（帧间强自相关）。本工具同时给出三种误差估计：
  1. naive SEM = SD/sqrt(n)            —— 偏乐观，仅作对照
  2. block SEM（块平均，默认 1 ns/块）  —— 推荐报告值
  3. ESS / 积分自相关时间 tau           —— 说明"203 帧实际等效几个独立样本"

输出
----
  md_traj_perframe.csv    逐帧指标
  md_traj_summary.json    汇总（含三种误差、收敛诊断）
  md_traj_panels.png      四联图（RMSD / 配体 RMSD+Rg / 接触与氢键 / dG 时序）

用法
----
  python md_traj_audit.py --md-dir <含 prod*.dcd 与 complex_raw.pdb 的目录> \
      --outdir <输出目录> [--top complex_raw.pdb] [--drop prod_seg1.dcd] \
      [--dg-csv e7_mmgbsa_perframe.csv] [--block-ns 1.0] [--stride 1]

注意
----
--drop 用于排除"续跑重叠段"。OpenMM 续跑会从检查点重算上一段尾部，
叠加会把同一时间窗计两遍、人为抬高 n 并低估 SEM。详见 md/Mv2/段文件使用说明.md。
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import struct
import sys

import numpy as np

# ---------------------------------------------------------------------------
# DCD / PDB 读取
# ---------------------------------------------------------------------------
def read_dcd(path: str):
    """读 OpenMM 写出的 CHARMM DCD，返回 (nframes, pos[nframes,natom,3]) 单位 nm。"""
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw[:4] != struct.pack("<i", 84) or raw[4:8] != b"CORD":
        raise ValueError(f"{path}: 不是 OpenMM 写出的 DCD")
    nframes = struct.unpack("<i", raw[8:12])[0]
    natom = struct.unpack("<i", raw[268:272])[0]
    pos = np.empty((nframes, natom, 3), dtype=np.float32)
    off = 276
    for i in range(nframes):
        nbyte = struct.unpack("<i", raw[off:off + 4])[0]
        if nbyte != 48:  # 无周期性时可能缺 cell 块
            off -= 4
        else:
            off += 8 + nbyte
        for d in range(3):
            nbyte = struct.unpack("<i", raw[off:off + 4])[0]
            if nbyte != 4 * natom:
                raise ValueError(f"{path} frame {i} axis {d}: block {nbyte} != {4*natom}")
            pos[i, :, d] = np.frombuffer(raw, dtype="<f4", count=natom, offset=off + 4)
            off += 8 + 4 * natom
    return nframes, pos


def read_pdb_atoms(path: str):
    """返回 (names, resnames, resids, elements) 四个 list，原子顺序与 DCD 一致。"""
    names, resnames, resids, elems = [], [], [], []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            names.append(line[12:16].strip())
            resnames.append(line[17:20].strip())
            resids.append(line[22:26].strip())
            e = line[76:78].strip() if len(line) >= 78 else ""
            if not e:  # 回退：由原子名推元素
                n = line[12:16].strip()
                e = n[0] if not n[0].isdigit() else n[1]
            elems.append(e.upper())
    return names, resnames, resids, elems


# ---------------------------------------------------------------------------
# 几何
# ---------------------------------------------------------------------------
def kabsch_rt(P: np.ndarray, Q: np.ndarray):
    """求把 P 叠到 Q 的旋转 R 与平移 t（P@R + t ≈ Q）。"""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, qc - pc @ R


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


def radius_of_gyration(x: np.ndarray) -> float:
    c = x.mean(0)
    return float(np.sqrt(((x - c) ** 2).sum(1).mean()))


def pair_counts(prot: np.ndarray, lig: np.ndarray, contact_cut=0.40, hb_cut=0.35,
                prot_polar=None, lig_polar=None):
    """接触重原子对数（<contact_cut nm）与氢键重原子对数（N/O…N/O < hb_cut nm）。"""
    d = np.linalg.norm(prot[:, None, :] - lig[None, :, :], axis=2)
    n_contact = int((d < contact_cut).sum())
    dmin = float(d.min())
    n_hb = 0
    if prot_polar is not None and len(prot_polar) and len(lig_polar):
        dd = np.linalg.norm(prot_polar[:, None, :] - lig_polar[None, :, :], axis=2)
        n_hb = int((dd < hb_cut).sum())
    return n_contact, dmin, n_hb


# ---------------------------------------------------------------------------
# 统计
# ---------------------------------------------------------------------------
def integrated_autocorr_time(x: np.ndarray, max_lag: int | None = None) -> float:
    """初始正序列估计量（Geyer），返回积分自相关时间 tau（单位：帧）。"""
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 8:
        return 1.0
    x = x - x.mean()
    var = (x * x).mean()
    if var <= 0:
        return 1.0
    max_lag = max_lag or min(n // 4, 200)
    tau = 0.0
    prev = np.inf
    for lag in range(1, max_lag + 1):
        c = float((x[:-lag] * x[lag:]).mean() / var)
        pair = c + (float((x[:-(lag + 1)] * x[lag + 1:]).mean() / var)
                    if lag + 1 < n else 0.0)
        if pair < 0:
            break
        if pair > prev:
            pair = prev
        tau += pair
        prev = pair
    return max(1.0, 1.0 + 2.0 * tau)


def block_sem(x: np.ndarray, block: int):
    """块平均 SEM：把序列切成 floor(n/block) 个块，用块均值估计标准误。"""
    x = np.asarray(x, dtype=float)
    nb = len(x) // block
    if nb < 2:
        return float("nan"), float("nan"), 0
    trimmed = x[: nb * block].reshape(nb, block).mean(1)
    return float(trimmed.std(ddof=1) / np.sqrt(nb)), float(trimmed.std(ddof=1)), nb


def describe(x: np.ndarray, block: int, label: str):
    x = np.asarray(x, dtype=float)
    n = len(x)
    sd = float(x.std(ddof=1)) if n > 1 else 0.0
    naive = sd / np.sqrt(n) if n else float("nan")
    bsem, bsd, nb = block_sem(x, block)
    tau = integrated_autocorr_time(x)
    ess = n / tau
    return {
        "label": label, "n": n,
        "mean": float(x.mean()), "sd": sd,
        "sem_naive": naive,
        "block_size_frames": block, "n_blocks": nb,
        "sem_block": bsem, "sd_between_blocks": bsd,
        "sem_inflation": (bsem / naive) if naive and np.isfinite(naive) and naive > 0 else None,
        "tau_frames": tau, "ess": ess,
        "min": float(x.min()), "max": float(x.max()),
    }


def half_convergence(x: np.ndarray, block: int):
    """前后半程均值差 / 块 SEM —— |比值| < 2 视为已收敛。"""
    x = np.asarray(x, dtype=float)
    h = len(x) // 2
    a, b = x[:h], x[h:]
    bsem, _, _ = block_sem(x, block)
    diff = float(b.mean() - a.mean())
    return {"first_half_mean": float(a.mean()), "second_half_mean": float(b.mean()),
            "diff": diff, "z": (diff / bsem) if bsem and np.isfinite(bsem) and bsem > 0 else None}


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--md-dir", required=True, help="含 prod*.dcd 与复合物 pdb 的目录")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--top", default="complex_raw.pdb", help="与 DCD 原子顺序一致的拓扑 pdb")
    ap.add_argument("--drop", default="prod_seg1.dcd", help="逗号分隔的重叠段文件名，空串表示不排除")
    ap.add_argument("--dg-csv", default="", help="可选：MM-GBSA 逐帧 csv，做同样的统计功效复核")
    ap.add_argument("--dg-col", default="dG_bind")
    ap.add_argument("--block-ns", type=float, default=1.0, help="块平均的块长（ns）")
    ap.add_argument("--frame-ns", type=float, default=0.1, help="帧间隔（ns）")
    ap.add_argument("--stride", type=int, default=1, help="隔帧取样，仅用于加速")
    ap.add_argument("--nprot", type=int, default=4851)
    ap.add_argument("--nlig", type=int, default=58)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    top_path = os.path.join(args.md_dir, args.top)
    names, resnames, resids, elems = read_pdb_atoms(top_path)
    natom_tot = len(names)
    if natom_tot < args.nprot + args.nlig:
        sys.exit(f"拓扑原子数 {natom_tot} < NPROT+NLIG={args.nprot + args.nlig}")

    idx = np.arange(natom_tot)
    is_lig = idx >= args.nprot
    is_prot = ~is_lig
    heavy = np.array([e != "H" for e in elems])
    ca = np.array([n == "CA" for n in names]) & is_prot
    bb = np.array([n in ("N", "CA", "C", "O") for n in names]) & is_prot
    polar = np.array([e in ("N", "O") for e in elems])
    prot_heavy = is_prot & heavy
    lig_heavy = is_lig & heavy
    lig_polar = is_lig & polar
    prot_polar = is_prot & polar
    # 蛋白/配体各自的本地掩码（长度 = NPROT / NLIG），供切片后索引
    ca_p = ca[:args.nprot]; bb_p = bb[:args.nprot]
    ph_p = prot_heavy[:args.nprot]; pp_p = prot_polar[:args.nprot]
    lh_l = lig_heavy[args.nprot:]; lp_l = lig_polar[args.nprot:]

    drop = {s.strip() for s in args.drop.split(",") if s.strip()}
    dcds = sorted(glob.glob(os.path.join(args.md_dir, "prod*.dcd")))
    kept = [p for p in dcds if os.path.basename(p) not in drop]
    if not kept:
        sys.exit("没有可用轨迹段")
    print(f"[scan] 全部段: {[os.path.basename(p) for p in dcds]}")
    print(f"[scan] 排除重叠段: {sorted(drop) or '无'} -> 使用 {len(kept)} 段")

    block = max(2, int(round(args.block_ns / args.frame_ns)))

    # ---- 第一遍：取全局第一帧作叠合参考 ----
    nf0, pos0 = read_dcd(kept[0])
    if nf0 == 0:
        sys.exit("首段无帧")
    ref = pos0[0].astype(np.float64)[:args.nprot + args.nlig] * 0.1
    ref_ca = ref[:args.nprot][ca_p]
    ref_bb = ref[:args.nprot][bb_p]
    ref_lig_h = ref[args.nprot:][lh_l]

    # ---- 第二遍：逐帧指标 ----
    recs = []
    for path in kept:
        nf, pos = read_dcd(path)
        name = os.path.basename(path)
        print(f"[dcd]  {name}: {nf} 帧")
        for i in range(0, nf, args.stride):
            p_all = pos[i].astype(np.float64)[:args.nprot + args.nlig] * 0.1  # Å -> nm
            prot_all = p_all[:args.nprot]
            lig_all = p_all[args.nprot:]
            Rm, tv = kabsch_rt(prot_all[ca_p], ref_ca)
            prot_fit = prot_all @ Rm + tv
            lig_fit = lig_all @ Rm + tv
            brmsd = rmsd(prot_fit[bb_p], ref_bb)
            lrmsd = rmsd(lig_fit[lh_l], ref_lig_h)
            rg = radius_of_gyration(prot_fit[ca_p])
            n_contact, dmin, n_hb = pair_counts(
                prot_fit[ph_p], lig_fit[lh_l],
                prot_polar=prot_fit[pp_p], lig_polar=lig_fit[lp_l])
            recs.append({"segment": name, "frame": i,
                         "t_ns": round(args.frame_ns * (len(recs) + 1), 2),
                         "rmsd_backbone_A": brmsd * 10.0,
                         "rmsd_ligand_A": lrmsd * 10.0,
                         "rg_protein_A": rg * 10.0,
                         "contacts_4A": n_contact,
                         "min_dist_A": dmin * 10.0,
                         "hbonds_3.5A": n_hb})
    if not recs:
        sys.exit("轨迹无帧")

    keys = ["rmsd_backbone_A", "rmsd_ligand_A", "rg_protein_A",
            "contacts_4A", "min_dist_A", "hbonds_3.5A"]
    summary = {k: describe([r[k] for r in recs], block, k) for k in keys}
    summary["n_frames_total"] = len(recs)
    summary["time_span_ns"] = [recs[0]["t_ns"], recs[-1]["t_ns"]]
    summary["convergence"] = {k: half_convergence([r[k] for r in recs], block)
                              for k in ("rmsd_backbone_A", "contacts_4A", "hbonds_3.5A")}

    # ---- 可选的 MM-GBSA 统计功效复核 ----
    dg_block = None
    if args.dg_csv and os.path.isfile(args.dg_csv):
        vals, ts = [], []
        with open(args.dg_csv, encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                try:
                    vals.append(float(row[args.dg_col]))
                    ts.append(float(row.get("t_ns", len(vals) * args.frame_ns)))
                except (KeyError, ValueError):
                    continue
        if len(vals) > 5:
            dg_block = describe(np.array(vals), block, args.dg_col)
            dg_block["convergence"] = half_convergence(np.array(vals), block)
            dg_block["source"] = os.path.basename(args.dg_csv)
            summary["mmgbsa"] = dg_block
            summary["mmgbsa_series"] = list(zip(ts, vals))

    with open(os.path.join(args.outdir, "md_traj_perframe.csv"), "w", newline="",
              encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0].keys()))
        w.writeheader()
        w.writerows(recs)
    with open(os.path.join(args.outdir, "md_traj_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    # ---- 图 ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        tt = [r["t_ns"] for r in recs]
        fig, ax = plt.subplots(2, 2, figsize=(11, 7.2), dpi=160)
        a = ax[0, 0]
        a.plot(tt, [r["rmsd_backbone_A"] for r in recs], lw=0.9, color="#c0392b")
        a.set_xlabel("Time (ns)"); a.set_ylabel("Backbone RMSD (Å)")
        a.set_title("Protein backbone RMSD")
        a = ax[0, 1]
        a.plot(tt, [r["rmsd_ligand_A"] for r in recs], lw=0.9, color="#8e44ad", label="Ligand RMSD")
        a.plot(tt, [r["rg_protein_A"] for r in recs], lw=0.9, color="#2980b9", ls="--", label="Protein Rg")
        a.set_xlabel("Time (ns)"); a.set_ylabel("Å"); a.legend(frameon=False, fontsize=8)
        a.set_title("Ligand RMSD & protein compactness")
        a = ax[1, 0]
        a.plot(tt, [r["contacts_4A"] for r in recs], lw=0.9, color="#16a085", label="Contacts <4 Å")
        a2 = a.twinx()
        a2.plot(tt, [r["hbonds_3.5A"] for r in recs], lw=0.7, color="#e67e22", alpha=.8, label="H-bonds <3.5 Å")
        a.set_xlabel("Time (ns)"); a.set_ylabel("Contacts"); a2.set_ylabel("H-bonds")
        a.set_title("Interface persistence")
        a = ax[1, 1]
        if dg_block:
            ts = [p[0] for p in summary["mmgbsa_series"]]
            vs = [p[1] for p in summary["mmgbsa_series"]]
            a.plot(ts[:len(vs)], vs, lw=0.7, color="#7f8c8d", alpha=.7, label="per frame")
            nb = len(vs) // block
            bm = np.array(vs[:nb * block]).reshape(nb, block).mean(1)
            a.plot([ts[:nb * block:block][i] for i in range(nb)], bm, "o-", ms=3,
                   lw=1.0, color="#2c3e50", label=f"block mean ({args.block_ns:g} ns)")
            a.axhline(np.mean(vs), color="#c0392b", ls=":", lw=1, label="grand mean")
            a.set_xlabel("Time (ns)"); a.set_ylabel("ΔG_bind (kcal/mol)")
            a.set_title("MM-GBSA binding free energy")
            a.legend(frameon=False, fontsize=7)
        else:
            a.axis("off")
        for a in ax.ravel():
            a.grid(alpha=.25, lw=.5)
        fig.tight_layout()
        png = os.path.join(args.outdir, "md_traj_panels.png")
        fig.savefig(png)
        print(f"[fig]  {png}")
    except Exception as exc:  # 图失败不影响数字结论
        print(f"[warn] 绘图跳过: {exc}")

    # ---- 控制台摘要 ----
    print("\n===== 轨迹稳定性 =====")
    for k in keys:
        s = summary[k]
        print(f"  {k:18s} {s['mean']:9.2f} ± {s['sem_block']:6.2f} (block SEM; "
              f"naive {s['sem_naive']:.3f}, 放大 {s['sem_inflation']:.1f}×, "
              f"tau {s['tau_frames']:.1f} 帧, ESS {s['ess']:.1f})")
    if dg_block:
        print("\n===== MM-GBSA 统计功效 =====")
        print(f"  {dg_block['label']}: {dg_block['mean']:.2f} ± {dg_block['sem_block']:.2f} kcal/mol "
              f"(block SEM; naive {dg_block['sem_naive']:.3f}, 放大 {dg_block['sem_inflation']:.1f}×, "
              f"tau {dg_block['tau_frames']:.1f} 帧, ESS {dg_block['ess']:.1f})")
        c = dg_block["convergence"]
        print(f"  前后半程: {c['first_half_mean']:.2f} vs {c['second_half_mean']:.2f} "
              f"(z={c['z']:.2f})")
    print(f"\n[done] -> {args.outdir}")


if __name__ == "__main__":
    main()
