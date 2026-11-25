# 18. 500 篇规模化复现（v1.7.0）

## 为什么做
23 篇手工复现（v1.3）只能给方向；要回答"中文/国际 NP 论文整体可复现性如何"，
必须按系统综述口径做规模化语料 + 逐篇参数挖掘 + 批量重算。

## 流水线（4 个命令，全部可复跑）
```bash
# 1) 语料筛选：PubMed 元数据 + OpenAlex 期刊影响指数（IF 代理）
python scripts/np_corpus.py screen --meta meta.json --journal-if jif.json \
  --if-min 5 --top 500 --drop-reviews --out sel.csv

# 2) PMC OA 全文收割（标题/摘要/正文/表格；同时保存 XML 便于复解析）
python scripts/np_corpus.py harvest --select sel.csv --outdir fulltext

# 3) 方法参数挖掘（药味/OB/DL/成分数/靶点数/交集数/对接能量/数据库/MD/湿实验）
python scripts/np_corpus.py extract --dir fulltext --select sel.csv --out methods.csv

# 4) 疾病靶点映射（OpenTargets EFO → 关联靶点，用于交集复现）
python scripts/np_corpus.py disease --methods methods.csv --cache disease_targets.json

# 5) 归类统计
python scripts/np_corpus.py classify --methods methods.csv --outdir stats

# 6) 批量复现 + 归因
python scripts/np_repro_batch.py run --methods methods.csv --herb-cache "DIR1,DIR2" --out repro.csv
python scripts/np_repro_batch.py inter --methods methods.csv --herb-cache DIR \
  --disease-json disease_targets.json --out inter.csv
python scripts/np_repro_batch.py stats --in repro.csv --outdir repro_stats
```

## 本次实测（2026-11-20，可复算）
| 阶段 | 数字 |
|---|---|
| PubMed 标题级检索（2018-2027） | 5,729 篇 |
| IF 代理 ≥5（OpenAlex 2yr_mean_citedness） | 2,423 篇 |
| 其中有 PMC OA 全文 | 1,665 篇 |
| 下载并解析 | 650 篇（缓冲池）→ 相关性硬筛后 **500 篇** |
| 500 篇构成 | IF 代理中位 8.0；中国单位 90%；期刊 Top：DDDT/Int J Mol Sci/Pharmaceuticals/Heliyon/Chin Med |
| 阈值披露 | OB 171 篇（30% 占 162）、DL 179 篇（0.18 占 168）→ 非默认阈值仅 ~3% |
| 对接能量 | 347 篇可抽，中位 −8.10 kcal/mol（范围 −15.0~−1.2） |
| 方法特征 | 湿实验 66%、动物 70%、MD 48%、TCGA 29%、GEO 19%、网络毒理学 4 篇 |
| 设计类型 | 单体 218、小方(2-4 味) 162、大方(≥5 味) 120 |

## 复现结果（诚实数字）
- 药材覆盖：中位 **71%**，全命中 179 篇（≥80% 覆盖 222 篇）
- 成分数可复现（±20%）：12/57（21%）；靶点数：18/118（15%）
- 交集复现（TCMSP∩OpenTargets vs 原文 GeneCards 口径）：39 篇可比，比值中位 **0.17**（中位 19 vs 122）
- 归因：未披露 48% > 药材漏检 42% > 靶点源差异 14% > 阈值差异 3%
- **牡蛎 101 篇未覆盖**——矿物药盲区在 500 篇尺度上被再次量化（v1.3 的定性结论成立）

## 关键方法学结论
1. **靶点源是最大单点差异**：66% 论文用 SwissTargetPrediction 扩靶，TCMSP 实测靶点只有其 1/5~1/3；
   交集数因此系统性偏低（比值 0.17）。这解释了 v1.3 "celastrol 124 vs 22"。
2. **阈值不是主要问题**（3%），但**披露缺失是**（48%）——多数论文的数字在图表/补充材料，正文不可抽。
3. 药材覆盖是可控变量：补 34 味 TCMSP 药材后覆盖率从 20%→71%，说明"数据库扩充"投入产出最高。
4. 单体论文占 44%（218/500）——这类无法用"药材→成分"复现，必须走化合物靶点预测（SwissTarget 类）。

## CNKI 中文核心（当前阻塞与交接）
- 检索入口：https://kns.cnki.net/kns8s/defaultresult/index
- 阻塞原因：WorkBuddy 的 `bsk` 浏览器 CLI 未安装且 GitHub 安装源不稳定；CNKI 下载需机构登录
- 交接步骤：① 安装 BrowserSkill（Tencent/BrowserSkill 官方）或直接用浏览器登录校园网/图书馆；
  ② 登录后按 `references/cnki-search.md` 检索"网络药理学 + 分子对接"、2018-2027、北大核心/CSCD；
  ③ 命中记录用 `np_corpus extract` 同口径解析（中文文本药味识别已内置 150+ 药名）。
- 已可替代的部分：PubMed 索引的中文核心期刊（Chin Med、Chin Herb Med、Chin J Integr Med 等）
  已包含在本批 500 篇内。

## 复用注意
- 抓取遵守 3 req/s（NCBI 无 key）；PMC 全文仅用于文本挖掘，不随 skill 分发。
- IF 为 OpenAlex 2yr_mean_citedness 代理，非 JCR 官方 IF，论文中引用需注明口径。
