# 28 · LINCS L1000 本地化：从"只能出模板"到"真的能查"

> 对应脚本：`np_lincs_local.py`（gctx/HDF5 流式读取）、`np_lincs.py`（签名构造与在线模板）。
> 解决的问题：审稿人问"你到底能不能做药物重定位查询"，而 CLUE API 需要授权、经常 502。

---

## 1 数据在哪、多大、要不要注册

| 数据集 | 内容 | 关键文件（GEO 补充材料） | 体积 | 注册 |
|---|---|---|---|---|
| **GSE92742**（Phase I, 2017） | 1.3M 实例；Level 5 为按 perturbagen 聚合的 moderated z-score | `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx` | ≈10 GB | 不需要 |
| **GSE70138**（Phase II, 2017） | 118,050 实例，新增化合物/细胞系 | `GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx` | ≈3 GB | 不需要 |
| 注释 | 基因/扰动/细胞系元数据 | `GSE92742_Broad_LINCS_gene_info.txt.gz`、`..._sig_info.txt.gz`、`..._cell_info.txt.gz` | <50 MB | 不需要 |
| CLUE（CMap 2020） | 官方最新版，含 tau 分数 | 需登录 `clue.io` | 数十 GB | 需要 |

> 文件名会随 GEO 版本更新，**不要硬编码**。先列目录再下载：

```bash
# 列出该 series 的全部补充文件（含真实文件名与大小）
curl -s https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/ | grep -o 'GSE92742[^"<]*'
```

国内网络建议用 `aspera` 或断点续传；10 GB 用 `wget -c` 即可：

```bash
mkdir -p data/lincs && cd data/lincs
BASE=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl
wget -c $BASE/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx
wget -c $BASE/GSE92742_Broad_LINCS_gene_info.txt.gz
wget -c $BASE/GSE92742_Broad_LINCS_sig_info.txt.gz
```

**只要 978 landmark 基因**时，可下 Level 3 的 978 维子集（体积小一个量级），
或用下面的"裁剪"脚本把 12,328 基因裁到 978 + 你的签名基因。

---

## 2 gctx 的内部结构

`.gctx` 就是 HDF5（用 `h5dump`/`h5py` 都能开）：

```
/data                      float32, (12,328 基因 × N 实例) 或转置；chunked
/meta/row_meta             基因维：pr_gene_symbol / pr_gene_id / pr_is_lm / pr_is_bing
/meta/col_meta             实例维：pert_iname / pert_id / cell_id / pert_idose / pert_time / pert_type
/meta/col_meta_group_size  每个 perturbagen 聚合了多少重复
```

先探结构（不会读大矩阵）：

```bash
python scripts/np_lincs_local.py probe --gctx data/lincs/GSE92742_....gctx
```

输出样例：

```json
{"status":"OK","size_gb":9.8,"data_shape":[12328,473647],"data_dtype":"float32",
 "chunks":[12328,1000],"row_meta_fields":["id","pr_gene_id","pr_gene_symbol","pr_is_lm"],
 "col_meta_fields":["id","pert_iname","cell_id","pert_idose","pert_time"],
 "n_row_meta":12328,"n_col_meta":473647}
```

`data_shape` 与 meta 长度不一致时脚本会自动判断基因轴（部分版本是转置存储），不匹配则直接报
`SHAPE_MISMATCH` 而不是给出错误结果。

---

## 3 为什么能"流式"而不爆内存

错误做法：`h5py.File(...)['/data'][:]` → 10 GB 全进内存。

正确做法：**只读签名基因对应的行**。

- 每行 = 1 个基因 × 473,647 个实例 × 4 B ≈ 1.9 MB；
- 200 个签名基因 ≈ 380 MB；读完一行就做 z-score，随后只保留 473,647 维的向量；
- 峰值内存 ≈ 380 MB + 一个 473k 向量，普通笔记本可跑。

代码里对应 `query()` 的逐行循环 + `np.vstack`。若签名基因 >1000，建议改成"分块读列"：

```python
for start in range(0, n_col, 20000):
    block = dset[:, start:start+20000]     # 基因 × 2万实例
    ...
```
（`np_lincs_local.py` 的 `batch` 参数预留了这一路径，默认逐行更省事。）

### 3.1 先裁一个小库（推荐）

反复查询时，把 12,328 × 473,647 裁成"landmark 978 + 你的 300 个签名基因 × 你关心的细胞系"，
体积可降到几百 MB，查询从分钟级降到秒级：

```python
import h5py, numpy as np
src = h5py.File("GSE92742_....gctx", "r")
genes = [g.decode() for g in src["/meta/row_meta"]["pr_gene_symbol"][:]]
keep = [i for i, g in enumerate(genes) if g in set(my_genes) | set(landmark_978)]
with h5py.File("lincs_subset.gctx", "w") as dst:
    dst.create_dataset("/data", data=src["/data"][sorted(keep), :], chunks=(len(keep), 2000))
    for grp in ("row_meta", "col_meta"):
        src.copy(f"/meta/{grp}", dst, name=f"/meta/{grp}")
```
（裁剪后 `row_meta` 必须同步裁行，否则 meta 与 data 长度不一致，脚本会报错——这是故意的。）

---

## 4 打分公式（写进 Methods 的版本）

给定签名 $U$（疾病上调）与 $D$（疾病下调），对每个扰动实例 $j$：

1. 对每个基因 $g$ 在实例维做 z-score：$z_{gj} = (x_{gj} - \mu_g)/\sigma_g$；
2. 实例得分：$S_j = \frac{1}{2}\left(\overline{z_{U,j}} - \overline{z_{D,j}}\right)$；
3. 同一 `pert_iname` 取最强（最负）3 个实例的均值，并记录实例数。

$S < 0$ = 该扰动**逆转**疾病签名（候选）；$S > 0$ = 模拟疾病签名（应避免）。
阈值 `S < -0.2` 记 `reversal`（经验值，报告里注明可调）。

**与官方 connectivity score 的差异（必须写在 Limitations）**：

| 项目 | 本实现 | CLUE 官方 |
|---|---|---|
| 归一化 | 基因维 z-score（跨实例） | plate/well 级 z-score + tau |
| 聚合 | top-3 实例均值 | weighted connectivity（tau 阈值） |
| 数值可比性 | 仅相对排序可比 | 可与 CMap 结果直接比 |

所以：**只有符号与相对排序可解释，绝对数值不要和文献里的 CMap score 对比。**

---

## 5 标准用法

```bash
# 1) 由核心靶点 + 方向表构造签名（复用 np_lincs）
python scripts/np_lincs.py signature --targets core.txt --direction dir.csv -o out/

# 2) 本地全库比对（10 GB 文件，200 基因签名约 2–5 分钟，取决于磁盘）
python scripts/np_lincs_local.py query --gctx data/lincs/GSE92742_....gctx \
    --signature out/signature.json --cell A375 --top 30 -o out/lincs

# 3) 也可直接给上下调基因列表
python scripts/np_lincs_local.py query --gctx X.gctx --up up.txt --down down.txt -o out/lincs
```

输出：
- `lincs_local_result.json`：签名覆盖率（hit_up/n_up、hit_down/n_down）、实例数、候选、caveat；
- `lincs_candidates.csv`：`pert_iname, score, n_instances, cell, dose, time`。

### 5.1 结果怎么用（别过度解读）

1. **候选不是结论**：表达层面可逆转 ≠ 结合该靶点，≠ 临床有效。
2. **细胞系依赖**：`--cell` 不限定会混入所有细胞系，建议按靶点组织选（如肝 → HEPG2）。
3. **覆盖率警告要看**：`coverage_up < 30%` 时结果基本无意义，通常是基因符号版本不匹配
   （gctx 用 HGNC symbol，老数据可能是 Entrez/Ensembl）。
4. **与方向校验联动**：逆转候选若作用于 `Target_Intervention_Mismatch` 的靶点，
   优先级应下调——机制自相矛盾的候选不要写进 Abstract。

---

## 6 常见报错与处理

| 报错 | 原因 | 处理 |
|---|---|---|
| `H5PY_MISSING` | 没装 h5py | `pip install h5py` |
| `NO_META` | gctx 版本字段名不同 | `probe` 看 `row_meta_fields`，用 `--gene-field` 指定（如 `pr_gene_id`） |
| `SHAPE_MISMATCH` | data 与 meta 长度不一致（裁剪时只裁了 data） | 同步裁 `row_meta` |
| `NO_OVERLAP` | 符号体系不匹配 | 用 `gene_info.txt.gz` 建映射：Entrez→symbol |
| `NO_CELL` | 细胞系名写法不同（`A375` vs `A375_1`） | 先 `probe` 看 `example_perts`/`cell_id` 取值 |
| 内存爆 | 签名基因上千 | 用 3.1 的裁剪脚本，或按列分块 |

---

## 7 一句话总结

**"API 挂了"不再是借口**：GEO 的 Level 5 矩阵免费、可断点续传、结构就是 HDF5；
只要签名基因数控制在几百，笔记本就能在分钟级完成全库比对。
把 `probe` 的输出、文件大小、sha256 前 16 位写进 Methods，审稿人就能复现你的查询。
