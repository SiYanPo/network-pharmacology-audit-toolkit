#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_figcheck.py —— v2.0.0 投稿图片自检（自动化产出 → 机检 → 修改建议）

检查项注册表（v2.1：消除"45 项"虚标——真实为 11 个检查项/图，逐项机读可枚举）
    见下方 CHECKS。每项有稳定 ID（C01–C11），报告按"检查项 × 图"统计，
    不再把"5 图 × 9 类 = 45 次执行"表述为"45 个独立检查项"。

命令：
  python np_figcheck.py --dir figs [--out report.md]
  python np_figcheck.py --list-checks        # 打印权威检查项清单
输出：figure_qa.md + figure_qa.csv（逐图逐项 PASS/WARN/FAIL）
"""
import argparse
import glob
import json
import os
import re

import pandas as pd
from PIL import Image

# 权威检查项注册表：唯一真源。--list-checks 与报告头部均引用此表，杜绝虚标。
CHECKS = [
    ("C01", "尺寸合规", "宽 80–180 mm 且高 ≤240 mm"),
    ("C02", "最小字号", "≥6 pt PASS；5–6 WARN；<5 FAIL"),
    ("C03", "文本重叠", "重叠计数=0 PASS；≤2 WARN；>2 FAIL"),
    ("C04", "文本越界", "越界计数=0 PASS；否则 FAIL"),
    ("C05", "PNG 分辨率", "≥300 dpi PASS（≥600 优秀）"),
    ("C06", "TIFF 分辨率", "≥300 dpi PASS"),
    ("C07", "TIFF 色彩模式", "RGB/L PASS；RGBA WARN"),
    ("C08", "PDF 矢量", "单页且字体内嵌（FontFile）"),
    ("C09", "SVG 字号", "无 <6 px 文本"),
    ("C10", "数据溯源", "fig.json 含 source_sha256_16"),
    ("C11", "文件齐备", "PDF/SVG/PNG/TIFF 四格式齐备"),
]
CHECK_N = len(CHECKS)


def check_one(meta_path):
    d = os.path.dirname(meta_path)
    name = os.path.basename(meta_path).replace(".json", "")
    meta = json.load(open(meta_path, encoding="utf-8"))
    rows = []

    def add(item, status, detail, fix=""):
        rows.append({"figure": name, "check": item, "status": status, "detail": detail, "fix": fix})

    # 尺寸
    w, h = meta.get("fig_size_mm", [None, None])
    if w:
        if 80 <= w <= 180 and h <= 240:
            add("尺寸", "PASS", f"{w}×{h} mm")
        else:
            add("尺寸", "FAIL", f"{w}×{h} mm 超出 80–180/≤240", "调整 --width 预设")
    # 字号
    mn = meta.get("min_font_pt")
    if mn is not None:
        if mn >= 6:
            add("最小字号", "PASS", f"{mn} pt")
        elif mn >= 5:
            add("最小字号", "WARN", f"{mn} pt（期刊下限）", "增大字号或减少标签")
        else:
            add("最小字号", "FAIL", f"{mn} pt 过小", "提高 rcParams 字号")
    # 文本重叠 / 越界
    ov = meta.get("text_overlaps", 0)
    add("文本重叠", "PASS" if ov == 0 else ("WARN" if ov <= 2 else "FAIL"),
        f"{ov} 处", "缩短标签/加大画幅" if ov else "")
    off = meta.get("text_offcanvas", 0)
    add("文本越界", "PASS" if off == 0 else "FAIL", f"{off} 处", "增加边距")
    # 位图
    for fmt, min_dpi in (("png", 300), ("tiff", 300)):
        p = os.path.join(d, f"{name}.{fmt}")
        if not os.path.exists(p):
            add(f"{fmt.upper()} 位图", "WARN", "缺失", "调整 --formats")
            continue
        im = Image.open(p)
        dpi = float(im.info.get("dpi", (0, 0))[0] or 0)
        wpx, hpx = im.size
        ok = dpi >= min_dpi
        add(f"{fmt.upper()} 分辨率", "PASS" if ok else "FAIL",
            f"{dpi:.0f} dpi, {wpx}×{hpx} px, {im.mode}",
            "提高 --dpi 至 600" if not ok else "")
        if fmt == "tiff" and im.mode not in ("RGB", "L"):
            add("TIFF 色彩模式", "WARN", im.mode, "投稿系统多要求 RGB")
    # 矢量
    pdf = os.path.join(d, f"{name}.pdf")
    if os.path.exists(pdf):
        raw = open(pdf, "rb").read()
        pages = raw.count(b"/Type /Page") - raw.count(b"/Type /Pages")
        embedded = b"FontFile" in raw
        add("PDF", "PASS" if pages <= 1 and embedded else "WARN",
            f"{max(pages,1)} 页, 字体{'已' if embedded else '未'}内嵌",
            "" if embedded else "设置 pdf.fonttype=42")
    svg = os.path.join(d, f"{name}.svg")
    if os.path.exists(svg):
        s = open(svg, encoding="utf-8", errors="ignore").read()
        tiny = len(re.findall(r"font-size:\s*([0-5](?:\.\d+)?)px", s))
        add("SVG", "PASS" if tiny == 0 else "WARN", f"过小字号 {tiny} 处",
            "检查标签字号" if tiny else "")
    # 溯源
    add("数据溯源", "PASS" if meta.get("source_sha256_16") else "WARN",
        meta.get("source_file") or "无",
        "缺少输入文件哈希" if not meta.get("source_sha256_16") else "")
    # C11 四格式齐备
    fmts = {ext: os.path.exists(os.path.join(d, "%s.%s" % (name, ext)))
            for ext in ("pdf", "svg", "png", "tiff")}
    missing = [k for k, v in fmts.items() if not v]
    add("文件齐备", "PASS" if not missing else "WARN",
        "四格式齐备" if not missing else ("缺 " + ",".join(missing)),
        "运行 np_figs --formats pdf,svg,png,tiff" if missing else "")
    return rows


def main():
    ap = argparse.ArgumentParser(description="np_figcheck: 投稿图片自检")
    ap.add_argument("--dir", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--list-checks", action="store_true",
                    help="打印权威检查项清单（C01–C%02d）后退出" % CHECK_N)
    a = ap.parse_args()
    if a.list_checks:
        print("np_figcheck 权威检查项（%d 项/图）：" % CHECK_N)
        for cid, name, rule in CHECKS:
            print("  %s  %-10s %s" % (cid, name, rule))
        print("统计口径：报告 = 检查项 × 图数；总条数 = %d × 图数。" % CHECK_N)
        print("历史文档中的『45 项』为 5 图 × 9 类执行次数，非独立检查项数，已废弃该表述。")
        return
    if not a.dir:
        ap.error("需 --dir 或 --list-checks")
    metas = sorted(glob.glob(os.path.join(a.dir, "*.json")))
    metas = [m for m in metas if not m.endswith("figures_manifest.json")]
    if not metas:
        raise SystemExit(f"[FAIL] {a.dir} 下没有 fig*.json（先运行 np_figs 生成）")
    all_rows = []
    for m in metas:
        all_rows += check_one(m)
    df = pd.DataFrame(all_rows)
    out_csv = a.out.replace(".md", ".csv") if a.out else os.path.join(a.dir, "figure_qa.csv")
    df.to_csv(out_csv, index=False, encoding="utf-8-sig")
    n_fig = df.figure.nunique()
    cnt = df.status.value_counts().to_dict()
    n_checks_actual = len({r["check"] for _, r in df.iterrows()})
    md = ["# 投稿图片自检报告", "",
          f"图片数: {n_fig} | 独立检查项: {n_checks_actual}（注册表 C01–C%02d 共 %d 项）| 执行条目: {len(df)}" %
          (CHECK_N, CHECK_N),
          f"PASS {cnt.get('PASS',0)} / WARN {cnt.get('WARN',0)} / FAIL {cnt.get('FAIL',0)}",
          "",
          "> 统计口径：{0} 图 × ≤{1} 项 = {2} 条执行记录；『执行条目』不等于『独立检查项』。".format(
              n_fig, CHECK_N, len(df)),
          ""]
    verdict = "可投稿" if cnt.get("FAIL", 0) == 0 and cnt.get("WARN", 0) <= 2 else (
        "需小修" if cnt.get("FAIL", 0) == 0 else "不可投稿（有 FAIL）")
    md.append(f"**总判定：{verdict}**")
    md += ["", "| 图 | 项目 | 状态 | 详情 | 建议 |", "|---|---|---|---|---|"]
    md += [f"| {r.figure} | {r.check} | {r.status} | {r.detail} | {r.fix} |" for _, r in df.iterrows()]
    # 汇总问题
    bad = df[df.status.isin(["FAIL", "WARN"])]
    md += ["", "## 需要处理的问题", ""]
    if len(bad) == 0:
        md.append("无。")
    else:
        for _, r in bad.iterrows():
            md.append(f"- **{r.figure} {r.check}**（{r.status}）：{r.detail} → {r.fix}")
    out_md = a.out or os.path.join(a.dir, "figure_qa.md")
    open(out_md, "w", encoding="utf-8").write("\n".join(md))
    print(f"{n_fig} 图 / {len(df)} 项 | PASS {cnt.get('PASS',0)} WARN {cnt.get('WARN',0)} FAIL {cnt.get('FAIL',0)} → {verdict}")
    print(f"→ {out_md}")


if __name__ == "__main__":
    main()
