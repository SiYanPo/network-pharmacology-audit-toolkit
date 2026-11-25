# COLOC 近似与 R `coloc` 一致性报告（2026-11-22）

## 环境（实测）

- R `R version 4.6.2 (2026-10-29 ucrt)`，路径 `C:/Program Files/R/R-4.6.2/bin/Rscript.exe`。
- R 包 `coloc 5.2.3`，`requireNamespace("coloc", quietly=TRUE)` → `TRUE`。
- 对照入口：`python scripts/np_coloc_compare.py run --loci <loci.tsv> -o <outdir>`。

## 关键缺陷（已修复）

1. **ABF 公式错误**：旧 `_labf` 多了 `−log(r)`；R `coloc::approx.bf.estimates` 源码为
   `0.5*(log(1−r)+r*z²)`。
2. **H4 组合错误**：旧实现 `sum(a+b)` 把全部 SNP 的联合 ABF 相加；R `combine.abf` 对联合 ABF 取
   `logsum`，H0/H1/H2/H3/H4 再统一归一化。
3. **H3 先验与 logdiff**：按官方源码使用 `p1*p2` 与
   `logdiff(logsum(l1)+logsum(l2), logsum(l1+l2))`；强共享信号时 H3 退化为 −inf。

修复前 20 位点：`r=0.7354, MAE=0.2538, max=0.9886`（verdict WARN）。
修复后 20 位点：`r=1.0000, MAE=0.0000, max=0.0000`（verdict PASS）。

## 证据文件

- 强信号 20 位点：`audit/coloc_cmp_20261122_exact/coloc_comparison.{json,csv}`。
- 部分非饱和 20 位点：`audit/coloc_cmp_soft_20261122/coloc_comparison.{json,csv}`；
  R PP4 范围 `~1e-45 – 1.0`，其中 2/20 位点处于 0.01–0.99，防止只用 0/1 饱和值刷相关性。
- 原始 R 输出：`*/r_coloc_pp.csv`；R 脚本：`*/run_coloc.R`。

## 结论与边界

- 结论：当前 `np_mr.coloc_abf` 与 `coloc::coloc.abf` 在**相同汇总统计与 W=0.15² 先验**下数值一致。
- 边界：真实 GWAS 仍受 LD/MAF、表型类型（R quant sd.prior=0.15²，cc sd.prior=0.2²）、
  `p1/p2/p12` 设定影响；Methods 应写清先验与 “validated against coloc.abf on 20 simulated loci”，
  不应写“未对齐”。
- 论文 DoD 第 7 条：状态 **DONE（模拟数据对齐）**；真实数据抽验仍列入
  `references/34-validation-roadmap.md` V3。
