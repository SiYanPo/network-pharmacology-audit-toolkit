#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_admet.py —— v1.4.0 ML 型 ADMET 预测模块（TDC 官方基准 + 外部验证）

数据源（均可自动化，无需注册）：
  - TDC 官方镜像（HuggingFace maomlab/TDC，经 hf-mirror.com 访问）：
    ADME.parquet + Tox.parquet，含官方 scaffold split（train/valid/test）。
  - DeepChem S3 BBBP.csv：与 TDC bbb_martins 不同文献来源的 BBB 数据集，
    用于真正的跨数据集外部验证。

协议（遵循 TDC ADMET Group 官方约定）：
  - scaffold split：train+valid 训练，test 仅评估（外部集）
  - 分类：AUROC（均衡）/ AUPRC（失衡，如 CYP 抑制）
  - 回归：MAE（主流）/ Spearman（半衰期、清除率等受结构外因素影响者）
  - 5 seed 重复：RF(random_state=0..4)，报告 mean±std

命令：
  download  --outdir DIR                     从 hf-mirror 下载 ADME/Tox parquet
  extract   --parquet F --outdir DIR [--tasks a,b,c]   按端点导出 CSV
  train     --datadir DIR --outdir DIR [--endpoints ...] [--seeds 5]
  external  --datadir DIR --modeldir DIR --outdir DIR   BBBP 跨库外部验证
  predict   --modeldir DIR --input SMILES或CSV --output CSV
  validate  --results CSV                    汇总并生成 Markdown 报告
"""
import argparse
import os
import sys
import time
import pickle
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

MIRROR = "https://hf-mirror.com/datasets/maomlab/TDC/resolve/main/single_instance_datasets"
S3_DEEPCHEM = "https://deepchemdata.s3.amazonaws.com/datasets"

# 官方指标约定（TDC ADMET Group / single_pred 文档）
METRIC = {
    "caco2_wang": ("reg", "MAE"),
    "lipophilicity_astrazeneca": ("reg", "MAE"),
    "ppbr_az": ("reg", "MAE"),
    "ld50_zhu": ("reg", "MAE"),
    "half_life_obach": ("reg", "Spearman"),
    "clearance_microsome_az": ("reg", "Spearman"),
    "clearance_hepatocyte_az": ("reg", "Spearman"),
    "hydrationfreeenergy_freesolv": ("reg", "MAE"),
    "vdss_lombardo": ("reg", "Spearman"),
}
AUPRC_ENDPOINTS = {  # 官方规定失衡二分类用 AUPRC
    "cyp2c9_veith", "cyp2d6_veith", "cyp3a4_veith", "cyp1a2_veith",
    "cyp2c19_veith",
}

# ---------------- 特征工程 ----------------

def featurize(smiles_list, verbose=False):
    """RDKit Morgan FP(2048,r2) + 物化描述符。返回 X, n_invalid"""
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    from rdkit.Chem import Descriptors as RDKitDesc
    DESC = ["MolWt", "MolLogP", "TPSA", "NumHAcceptors", "NumHDonors",
            "NumRotatableBonds", "RingCount", "FractionCSP3", "HeavyAtomCount",
            "NHOHCount", "NOCount", "NumAliphaticRings", "NumAromaticRings"]
    fps, dess, keep = [], [], []
    n_bad = 0
    for i, s in enumerate(smiles_list):
        m = Chem.MolFromSmiles(str(s))
        if m is None:
            fps.append(np.zeros(2048, dtype=np.float32)); dess.append(np.zeros(len(DESC)))
            keep.append(False); n_bad += 1
            continue
        fp = gen.GetFingerprint(m)
        arr = np.zeros(2048, dtype=np.float32); DataStructs.ConvertToNumpyArray(fp, arr)
        d = RDKitDesc.CalcMolDescriptors(m)
        fps.append(arr); dess.append([float(d.get(k, 0) or 0) for k in DESC])
        keep.append(True)
        if verbose and (i + 1) % 5000 == 0:
            print(f"    featurize {i+1}", flush=True)
    X = np.hstack([np.array(fps), np.array(dess)])
    return X, np.array(keep), n_bad


def _feat_cache(datadir, task, split_smiles):
    """缓存特征 npz：相同 smiles 集不重复计算"""
    import hashlib
    key = hashlib.md5("\n".join(map(str, split_smiles)).encode()).hexdigest()[:10]
    p = os.path.join(datadir, f"_feat_{task}_{key}.npz")
    if os.path.exists(p):
        z = np.load(p)
        return z["X"], z["keep"]
    X, keep, nbad = featurize(split_smiles, verbose=True)
    np.savez_compressed(p, X=X, keep=keep)
    if nbad:
        print(f"    [warn] {task}: {nbad} 条 SMILES 解析失败（置零特征）")
    return X, keep

# ---------------- 训练 ----------------

def train_endpoint(datadir, task, seeds, modeldir):
    from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
    from sklearn.metrics import roc_auc_score, average_precision_score, \
        mean_absolute_error, f1_score
    from scipy.stats import spearmanr

    csv = os.path.join(datadir, f"{task}.csv")
    if not os.path.exists(csv):
        return None
    df = pd.read_csv(csv)
    need = {"SMILES", "Y", "split"}
    if not need.issubset(df.columns):
        print(f"  [skip] {task}: 缺列 {need - set(df.columns)}")
        return None
    trv = df[df.split.isin(["train", "valid"])]
    te = df[df.split == "test"]
    if len(trv) < 50 or len(te) < 30:
        print(f"  [skip] {task}: 样本过少 train={len(trv)} test={len(te)}")
        return None

    print(f"  [{task}] train={len(trv)} test={len(te)}", flush=True)
    Xtr, ktr = _feat_cache(datadir, task, trv.SMILES.tolist())
    Xte, kte = _feat_cache(datadir, task, te.SMILES.tolist())
    ytr = trv.Y.values[ktr]; yte = te.Y.values[kte]
    Xtr, Xte = Xtr[ktr], Xte[kte]

    is_reg = METRIC.get(task, ("cls",))[0] == "reg"
    metric_name = METRIC.get(task, ("cls", "AUROC" if task not in AUPRC_ENDPOINTS else "AUPRC"))[1]

    scores = []
    for seed in seeds:
        if is_reg:
            m = RandomForestRegressor(n_estimators=200, n_jobs=-1, random_state=seed)
            m.fit(Xtr, ytr)
            p = m.predict(Xte)
            sc = mean_absolute_error(yte, p) if metric_name == "MAE" \
                else spearmanr(yte, p).correlation
        else:
            m = RandomForestClassifier(n_estimators=200, n_jobs=-1,
                                       class_weight="balanced_subsample",
                                       random_state=seed)
            m.fit(Xtr, ytr)
            p = m.predict_proba(Xte)[:, 1]
            if metric_name == "AUPRC":
                sc = average_precision_score(yte, p)
            else:
                sc = roc_auc_score(yte, p)
        scores.append(float(sc))
        if seed == seeds[0]:
            with open(os.path.join(modeldir, f"{task}.pkl"), "wb") as f:
                pickle.dump({"model": m, "task": task, "is_reg": is_reg,
                             "metric": metric_name, "train_n": int(len(ytr))}, f)
    # 次级指标（分类补 F1@0.5；回归补 Spearman/MAE）
    extra = ""
    if not is_reg:
        f1 = f1_score(yte, (p >= 0.5).astype(int))
        extra = f"F1@0.5={f1:.3f}"
    mean, std = float(np.mean(scores)), float(np.std(scores))
    print(f"    {metric_name}={mean:.3f}±{std:.3f} {extra}", flush=True)
    return dict(task=task, type="regression" if is_reg else "classification",
                metric=metric_name, value=mean, std=std, n_train=int(len(ytr)),
                n_test=int(len(yte)), extra=extra)


DEFAULT_ENDPOINTS = [
    # Absorption
    "caco2_wang", "hia_hou", "pgp_broccatelli", "bioavailability_ma",
    "lipophilicity_astrazeneca",
    # Distribution
    "bbb_martins", "ppbr_az",
    # Metabolism
    "cyp3a4_veith", "cyp2d6_veith", "cyp2c9_veith", "cyp1a2_veith",
    # Excretion
    "clearance_microsome_az", "half_life_obach",
    # Toxicity
    "ld50_zhu", "herg", "ames", "dili", "clintox", "skin_reaction",
    # 与复现论文机制关联的 Tox21 端点
    "tox21_SR-ARE", "tox21_NR-PPAR-gamma", "tox21_NR-AhR", "tox21_NR-ER",
]

def cmd_train(args):
    datadir, modeldir = args.datadir, args.modeldir
    os.makedirs(modeldir, exist_ok=True)
    eps = args.endpoints.split(",") if args.endpoints else DEFAULT_ENDPOINTS
    seeds = list(range(args.seeds))
    rows = []
    for t in eps:
        t = t.strip()
        r = train_endpoint(datadir, t, seeds, modeldir)
        if r:
            rows.append(r)
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(args.outdir, "admet_metrics.csv"),
                                  index=False, encoding="utf-8-sig")
        print(f"\n完成 {len(rows)}/{len(eps)} 端点 → {args.outdir}/admet_metrics.csv")

# ---------------- 外部验证（跨数据集 BBB）----------------

def cmd_external(args):
    """训练 bbb_martins 的模型去预测另一文献来源的 BBBP 测试集"""
    import requests, io
    from sklearn.metrics import roc_auc_score
    datadir, modeldir = args.datadir, args.modeldir
    pkl = os.path.join(modeldir, "bbb_martins.pkl")
    if not os.path.exists(pkl):
        print("[FAIL] 无 bbb_martins 模型，请先 train")
        return
    bbp = os.path.join(datadir, "BBBP.csv")
    if not os.path.exists(bbp):
        r = requests.get(f"{S3_DEEPCHEM}/BBBP.csv", timeout=120)
        open(bbp, "wb").write(r.content)
    df = pd.read_csv(bbp)
    col = "p_np" if "p_np" in df.columns else df.columns[-2]
    X, keep, nbad = featurize(df["smiles"].tolist())
    y = df[col].values[keep]
    X = X[keep]
    art = pickle.load(open(pkl, "rb"))
    prob = art["model"].predict_proba(X)[:, 1]
    auc = roc_auc_score(y, prob)
    out = dict(source="BBBP (DeepChem S3, 不同文献来源)", n=len(y),
               positive_rate=float(y.mean()), cross_auc=round(float(auc), 3))
    print("跨库外部验证（bbb_martins 训练 → BBBP 评估）:", out)
    pd.DataFrame([out]).to_csv(os.path.join(args.outdir, "external_validation.csv"),
                               index=False, encoding="utf-8-sig")

# ---------------- 预测 ----------------

def cmd_predict(args):
    modeldir = args.modeldir
    smis = []
    names = []
    if os.path.exists(args.input):
        df = pd.read_csv(args.input, encoding="utf-8-sig")
        scol = "SMILES" if "SMILES" in df.columns else df.columns[-1]
        ncol = "compound" if "compound" in df.columns else None
        for _, r in df.iterrows():
            smis.append(r[scol])
            names.append(str(r[ncol]) if ncol else str(r[scol]))
    else:
        for s in args.input.replace(";", ",").split(","):
            s = s.strip()
            if s:
                smis.append(s); names.append(s)
    X, keep, nbad = featurize(smis)
    rows = []
    for fn in sorted(os.listdir(modeldir)):
        if not fn.endswith(".pkl"):
            continue
        art = pickle.load(open(os.path.join(modeldir, fn), "rb"))
        m = art["model"]
        if art["is_reg"]:
            pred = m.predict(X)
            val = pred
            flag = ""
        else:
            prob = m.predict_proba(X)[:, 1]
            val = prob
            flag = ["high" if v >= 0.7 or v <= 0.3 else "uncertain" for v in prob]
        for i, nm in enumerate(names):
            rows.append({"compound": nm, "endpoint": art["task"],
                         "value": round(float(val[i]), 3),
                         "type": "regression" if art["is_reg"] else "prob",
                         "confidence": (flag[i] if not art["is_reg"] else "—"),
                         "valid_smiles": bool(keep[i])})
    out = pd.DataFrame(rows).pivot_table(index="compound", columns="endpoint",
                                         values="value", aggfunc="first")
    out.to_csv(args.output, encoding="utf-8-sig")
    detail = pd.DataFrame(rows)
    detail.to_csv(args.output.replace(".csv", "_detail.csv"), index=False,
                  encoding="utf-8-sig")
    print(f"预测 {len(names)} 化合物 × {len(set(detail.endpoint))} 端点 → {args.output}")
    if nbad:
        print(f"[warn] {nbad} 条 SMILES 无效")

# ---------------- 报告 ----------------

def cmd_validate(args):
    df = pd.read_csv(args.results, encoding="utf-8-sig")
    lines = ["# TDC ADMET ML 基线复现报告（scaffold split，官方指标）", "",
             f"端点数: {len(df)} | 训练协议: RF-200×5seeds, train+valid 训练, test 外部评估", "",
             "| 端点 | 类型 | 官方指标 | 本机值 (mean±std) | n_train/n_test |", "|---|---|---|---|---|"]
    for _, r in df.iterrows():
        lines.append(f"| {r['task']} | {r['type']} | {r['metric']} | "
                     f"{float(r['value']):.3f}±{float(r['std']):.3f} | {r['n_train']}/{r['n_test']} |")
    ext = os.path.join(os.path.dirname(args.results), "external_validation.csv")
    if os.path.exists(ext):
        e = pd.read_csv(ext).iloc[0]
        lines += ["", "## 跨库外部验证", "",
                  f"bbb_martins 训练模型 → BBBP 独立测试集：AUROC = **{float(e['cross_auc']):.3f}** (n={int(e['n'])}, 正类率={float(e['positive_rate']):.3f})"]
    md = args.results.replace(".csv", "_report.md")
    open(md, "w", encoding="utf-8").write("\n".join(lines))
    print(f"报告 → {md}")

# ---------------- 下载/抽取 ----------------

def cmd_download(args):
    import requests
    os.makedirs(args.outdir, exist_ok=True)
    for name in ["ADME", "Tox"]:
        p = os.path.join(args.outdir, f"{name}.parquet")
        if os.path.exists(p):
            print(f"[skip] {name} 已存在"); continue
        u = f"{MIRROR}/{name}/train-00000-of-00001.parquet"
        print(f"下载 {name} ...", flush=True)
        r = requests.get(u, timeout=1800)
        r.raise_for_status()
        open(p, "wb").write(r.content)
        print(f"  {name}: {len(r.content)//1024//1024} MB")

def cmd_extract(args):
    os.makedirs(args.outdir, exist_ok=True)
    tasks = set(args.tasks.split(",")) if args.tasks else None
    # 通过文件大小自动判断 parquet 里有哪些类别
    import pyarrow.parquet as pq
    # ADME 端点优先从 ADME 类 parquet 抽；Tox 端点从 Tox 类抽
    adme_tasks, tox_tasks = set(), set()
    for f in (args.parquet or []):
        df = pd.read_parquet(f, columns=["Task"])
        ts = set(df.Task.unique())
        if "caco2_wang" in ts:
            adme_tasks |= ts
        else:
            tox_tasks |= ts
    for t in (tasks or (adme_tasks | tox_tasks)):
        src = args.parquet[0] if t in adme_tasks else args.parquet[-1]
        df = pd.read_parquet(src)
        sub = df[df.Task == t][["Drug_ID", "SMILES", "Y", "split"]]
        sub.to_csv(os.path.join(args.outdir, f"{t}.csv"), index=False)
        print(f"  {t}: {len(sub)} 行")
    print(f"端点 CSV → {args.outdir}")


def main():
    ap = argparse.ArgumentParser(description="np_admet: TDC ML ADMET 模块")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("download"); s.add_argument("--outdir", default=".")
    s = sub.add_parser("extract")
    s.add_argument("--parquet", nargs="+")
    s.add_argument("--outdir", default="admet_endpoints")
    s.add_argument("--tasks", default=None)
    s = sub.add_parser("train")
    s.add_argument("--datadir", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--modeldir", default=None)
    s.add_argument("--endpoints", default=None)
    s.add_argument("--seeds", type=int, default=5)
    s = sub.add_parser("external")
    s.add_argument("--datadir", required=True)
    s.add_argument("--modeldir", required=True)
    s.add_argument("--outdir", required=True)
    s = sub.add_parser("predict")
    s.add_argument("--modeldir", required=True)
    s.add_argument("--input", required=True, help="SMILES 串或 CSV 路径")
    s.add_argument("--output", required=True)
    s = sub.add_parser("validate"); s.add_argument("--results", required=True)
    args = ap.parse_args()
    if args.cmd == "train":
        args.modeldir = args.modeldir or args.outdir
    g = {"download": cmd_download, "extract": cmd_extract, "train": cmd_train,
         "external": cmd_external, "predict": cmd_predict, "validate": cmd_validate}
    g[args.cmd](args)


if __name__ == "__main__":
    main()
