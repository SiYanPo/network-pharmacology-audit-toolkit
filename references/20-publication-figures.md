# 20. 投稿级图片自动化与自检（v2.0.0）

## 为什么要做
500 篇审计：heatmap 46%、bubble 33%、PPI 网络几乎标配，但作者仍需在 Origin/Cytoscape 手工排版，
常见退修理由集中在——分辨率不足、字号过小、标签重叠、RGB/CMYK、矢量缺失。
本模块把"出图 → 自检 → 修正"做成闭环，一次生成四件套。

## 一键出图
```bash
# 富集气泡图（GO/KEGG）
python scripts/np_figs.py bubble --enrich enrich_significant.csv --top 12 --width double \
  --outdir figs --name fig2_bubble

# 通路柱图
python scripts/np_figs.py bar --enrich enrich_significant.csv --top 10 --outdir figs

# 成分-靶点热图（长表：compound,gene[,value]）
python scripts/np_figs.py heatmap --long compound_targets.csv \
  --top-compounds 14 --top-genes 20 --outdir figs --name fig3_heatmap

# PPI 网络（边表：node1,node2[,weight]）
python scripts/np_figs.py ppi --edges ppi_edges_clean.tsv --top 40 --outdir figs

# 研究筛选流程图（PRISMA 风格）
python scripts/np_figs.py prisma --counts prisma_counts.json --outdir figs
```
每个图输出：`fig.pdf`（矢量，投稿主用）、`fig.svg`（可编辑）、`fig.png`（600 dpi）、
`fig.tiff`（600 dpi, LZW, RGB）、`fig.json`（元数据：尺寸/字号/重叠/数据源 SHA256）。

## 期刊尺寸预设
| 预设 | 宽度 | 适用 |
|---|---|---|
| `single` | 85 mm | 单栏 |
| `1.5col` | 114 mm | 1.5 栏 |
| `double` | 174 mm | 双栏（默认） |
| `full` | 180 mm | 满版 |

字体 Arial 7 pt、标题 8 pt；配色 Okabe-Ito（色盲友好）；PDF 字体 Type42 内嵌（可编辑）。

## 自检（投稿前必跑）
```bash
python scripts/np_figcheck.py --dir figs --out figs/figure_qa.md
```
检查 11 项/图（C01–C11，`np_figcheck.py --list-checks` 可枚举）：分辨率（≥300 PASS/≥600 EXCELLENT）、物理尺寸、最小字号（≥6 pt）、
文本重叠、文本越界（区分真裁切与刻度外溢）、PDF 页数与字体嵌入、SVG 过小字号、
TIFF 色彩模式、数据溯源哈希。输出 `figure_qa.md`（人读）+ `figure_qa.csv`（机读）+ 总判定
（可投稿 / 需小修 / 不可投稿）。

## 实测（2026-11-20）
- 5 图 × 11 项 = **55/55 PASS → 可投稿**（`投稿图_v2/`）
- 自检确实抓到真问题并修掉：TIFF 为 RGBA（→ 转 RGB）、
  PPI 标签越界（→ 力导向坐标归一化到 10–90%）、
  柱图末位刻度溢出（→ 轴宽/余量调整）、气泡图行距不足导致标签重叠（→ 每行 6 mm）。

## 一键全流程（推荐）
```bash
python scripts/np_figs.py batch \
  --enrich enrich_significant.csv \
  --long compound_targets.csv \
  --edges ppi_edges_clean.tsv \
  --counts prisma_counts.json \
  --outdir 投稿图   # 默认自动跑 np_figcheck
```
实测：5 图 × 11 项 = 55/55 PASS，总耗时 ~40 s，产物体积 3.6 MB（PDF+SVG+PNG+TIFF 各 5）。

## PPI 布局说明
按连通分量分块（component-grid）：大簇用 spring 布局、2 节点簇紧凑放置、单节点居中，
不再丢弃小簇（v2.0 实测从"仅最大簇 6 节点"修正为"15 节点 23 边全展示"）。

## 与 Cytoscape 的分工
- 布局与美化：本模块直接出投稿版（无需 Cytoscape）
- 交互调整：`np_cyto.py` 输出 XGMML/SIF，可导入 Cytoscape 继续编辑
- `.cys` 为二进制会话格式，无法程序生成（XGMML 为官方可编程替代）

## 常见期刊要求速查
| 项 | 常见要求 | 本模块默认 |
|---|---|---|
| 位图分辨率 | ≥300 dpi（线条图 600–1200） | 600 dpi |
| 字号 | ≥6 pt（8 pt 更稳） | 6–8 pt |
| 字体 | Arial / Helvetica | Arial（缺失回退 DejaVu Sans） |
| 色彩 | RGB 或 CMYK | RGB（TIFF 明确 RGB） |
| 矢量 | PDF/EPS 优先 | PDF + SVG |
| 宽度 | 85 / 114 / 174 mm | 预设可选 |
