# 27 · 软件化三件套：Benchmark、Web UI 与 CI

> 适用场景：把网络药理学流水线从"能跑"变成"可发表、可被引用、可被复现"的软件。
> 对应脚本：`np_benchmark.py`、`app.py`、`np_uireport.py`、`tests/`、`.github/workflows/test.yml`。

---

## 1 为什么软件类论文必须先有 Benchmark

审稿人不会因为 PPI 图漂亮就相信预测能力。软件类文章（Briefings in Bioinformatics / NAR / GigaScience）
的审稿意见里，"How do you know it works?" 通常以三种形式出现：

1. **和谁比**：有没有基线？朴素基线（支持成分数）打不过就说明加权毫无价值。
2. **怎么定义对**：金标准从哪来？谁标注的？覆盖多少？
3. **换一批数据还行吗**：跨库/跨方剂是否稳定。

没有这三条，方法部分只能写"我们实现了……"，不能写"我们优于……"。

### 1.1 两种基准，别混用

| 类型 | 金标准形态 | 适用指标 | 什么时候用 |
|---|---|---|---|
| **靶点级**（target-level） | 一个方剂 → 一小组强证据靶点（20 个左右） | AUC、AUPRC、Recall@K、EF@K、超几何 p | 评价"排序能力"；金标准必须**稀疏**（占候选池 <15%） |
| **边级**（edge-level） | 另一数据库的（成分, 靶点）边集 | 敏感性、精确率、Jaccard、过预测倍数 | 评价"覆盖与特异性"；当金标准是候选池子集时**必须**用这个 |

> 血泪教训：若金标准覆盖了候选池的 42%，靶点级 Recall@K 会被基率压死
> （随机抽 20 个也能中 8 个），此时 AUC 仍可用，但 Recall@K 与 EF 不可作主指标。
> 我们把这个现象写进了报告正文，而不是换一个好看的金标准——**换金标准适配指标＝学术不端**。

### 1.2 金标准的构建规则（可直接写进 Methods）

1. **只收强证据**：≥2 个独立来源，或有直接实验支持（L3/L4）→ `core`。
2. **弱证据单列**：单一来源或对接预测（L1）→ `extended`，**从阴性池中剔除**，不当负样本。
   把不确定项当负样本会系统性压低所有方法的 AUC。
3. **标注证据级别**：每个靶点带 `L0/L1/L3/L4`，报告里可切换 `--tier strict-l3` 做敏感性分析。
4. **写清 coverage**：金标准里有多少靶点出现在候选池中。coverage < 60% 时 Recall 上限被压死，
   报告必须显式警告（`coverage_warning`）。

### 1.3 指标定义（避免审稿人抓口径）

- **AUC**：Mann–Whitney U 形式，并列分数取平均秩；1000 次 bootstrap 给 95% CI。
- **AUPRC / AP**：`AP = Σ (P@k × ΔRecall)`，正类稀少时比 AUC 敏感。
- **Recall@K**：前 K 个里命中数 / 金标准总数；随机期望 `K/N`，同时给超几何 p。
- **EF@K**：`(hits/K) / (n_pos/N)`。EF≈1 表示与随机无异，**这一点必须写出来**。
- **过预测倍数**：`|pred| / |gold|`，>2 说明库宽、特异性低，需下游过滤。

### 1.4 内置金标准（`data/benchmark_goldstandard.json`）

4 个已上市/经典方剂：复方丹参滴丸（CDDP）、黄连解毒汤（HLJDT）、连花清瘟（LHQW）、血府逐瘀汤（XFZYT）。
每个方剂含 `core`（15–20 个）、`extended`、`kegg` 通路 ID、标志性成分、适应症。

来源与局限写在 JSON 的 `provenance` 段：DrugBank/TTD/ChEMBL 机制记录 + 标志性成分药理实验文献。
**必须承认**：这是人工整理先验，不是系统评价；存在 publication bias（研究得多的靶点更容易进金标准），
会抬高所有方法的 AUC。因此它只用于**版本间回归**与**方法间横向比较**，不能解释为临床准确率。

### 1.5 标准用法

```bash
# 列出内置方剂
python scripts/np_benchmark.py list

# 靶点级：用流水线输出的排序表
python scripts/np_benchmark.py evaluate --formula CDDP \
    --pred np_out/05_ppi/hubs.csv --tier strict --topk 10,20,50 -o np_out/benchmark

# 靶点级：只有边表时用内置基线（support / penalized 两种，预先注册，不做事后挑选）
python scripts/np_benchmark.py evaluate --formula CDDP \
    --edges np_out/03_compounds/compound_target_edges.csv --weight-mode penalized

# 边级：两个数据库比覆盖
python scripts/np_benchmark.py concordance --pred-edges mine.csv --gold-edges other_db.csv \
    --label "TCMSP vs ETCM2" -o np_out/benchmark

# 汇总出报告（含 Methods 段落与 Limitations 四条）
python scripts/np_benchmark.py report --dir np_out/benchmark
```

输出：`benchmark_<ID>.json`、`benchmark_<ID>_hits.csv`（逐靶点 rank/score/是否命中，便于人工核查）、
`roc_<ID>.svg`、`benchmark_report.md`。

### 1.6 实测样例（益智胶囊公开数据，2026-11-22）

| 指标 | 朴素基线 B0（支持成分数） | 基线 B1（多靶点成分降权） |
|---|---|---|
| AUC [95% CI] | 0.599 [0.55–0.65] | **0.705 [0.65–0.76]** |
| AUPRC | 0.471 | **0.526** |
| Recall@20 | 0.04 | 0.04 |
| 边级敏感性（vs ETCM2） | — | **1.000**（502/502） |
| 边级精确率 / 过预测 | — | 0.244 / 4.10× |

三条可写进论文的结论：

1. **promiscuity penalty 确实有用**：同一批边，AUC 从 0.599 提到 0.705（+0.106），说明"多靶点成分降权"
   不是装饰，而是排序质量的决定因素。
2. **覆盖完整、特异性有限**：金标准边一条不落（敏感性 1.00），但边数是金标准的 4.1 倍。
   这就是为什么下游必须上零模型与诱饵对接——**数据库并集必然带假阳性**。
3. **Recall@20 不可作主指标**：金标准覆盖候选池 42%，随机期望就有 8.5 个命中，我们只中 6 个（EF 0.71）。
   这不是算法崩了，而是基率太高导致该指标失去分辨力——报告里必须同时给 coverage 与 EF。

---

## 2 一键式 Web UI（`app.py`，Streamlit）

### 2.1 目标与边界

目标：输入方剂+疾病 → 点一次 → 进度条跑完 → 下载带水印报告。**不是**把 20 个命令搬到网页上。

边界（写进 README，防止用户误解）：
- UI 只做**编排与呈现**，不重算任何指标；报告数字全部读磁盘产物。
- 进度条由 `np_pipeline` 打印的 `[node_id]` 事件驱动，分母来自同一份 `config/pipeline.yaml`，
  **不是假动画**。
- 没跑门禁 = `UNVERIFIED` 水印；跑了但硬性项没过 = `DRAFT`；全过 = `READY_FOR_SUBMISSION`。
- PDF 依赖 `weasyprint`；缺失时只出 HTML，并明确提示用浏览器打印（不静默失败）。

### 2.2 结构

```
app.py
├─ plan_nodes(enable_flags)         # 与引擎同源：load_config → select_nodes → topo_sort
├─ build_inputs(...)                # 组装 inputs.json（含 herbs.txt 落盘）
├─ run_stream(cmd, on_line)         # Popen 流式读 stdout，逐行回调
├─ page_analyze()                   # 页签1：输入 → 跑 DAG → 出报告
├─ page_benchmark()                 # 页签2：金标准评测 / 跨库一致性
├─ page_causal()                    # 页签3：MR + 动作表 → 方向冲突
├─ page_env()                       # 页签4：环境探测 + 单模块 selftest
└─ _make_report(outdir,label,gate)  # finalcheck → np_uireport.build → HTML/PDF/ZIP
```

关键实现细节（都是踩过的坑）：

1. **`sys.exit()` 不能用**：Streamlit 以 `__main__` 执行脚本，`SystemExit` 会打断 runner，
   前端报 `missing ScriptRunContext` / `session_state has no key $$...` 之类怪错。
2. **路径必须加引号**：`--label "{params.label}"`，否则含空格的中文 label 被 `split_cli` 切成多个参数，
   `np_report` 直接 argparse 报错（真实故障：报告节点 failed）。
3. **子进程而不是 import**：`subprocess.Popen([sys.executable, "-u", np_pipeline.py, ...])`，
   避免模块级全局状态污染，也让 UI 崩了不影响引擎。
4. **`st.session_state` 不存大对象**：只存路径，下载时现读字节。

### 2.3 报告生成（`np_uireport.py`）

与 `np_report.py`（阶段审计）、`np_htmlreport.py`（投稿排版）三者的分工：
`np_uireport` 是**交付版**——不挑产物，缺什么写"未运行"，绝不静默省略；带平铺斜向水印与角标。

11 个章节：运行清单 → 配伍禁忌 → 成分 → 交集 → PPI → 零模型 → 富集 → 同源 → 可选模块
（MR/方向/基准）→ 图件 → 产物 SHA256 清单。

水印实现：`position:fixed` 平铺 `READY FOR SUBMISSION` / `DRAFT` / `UNVERIFIED`，
透明度 0.055，`@media print` 下保留；角标用固定定位，打印首末页都可见。

---

## 3 单元测试与 CI

### 3.1 测什么（优先级从高到低）

1. **会静默算错的地方**：MR 的 IVW/Egger/Steiger/LD/COLOC、decoy 的 Z 阈值、AUC 的并列秩、超几何尾概率。
2. **会静默丢数据的地方**：别名四层映射、含糊靶点是否被误当负样本、`--resume` 的假成功判定。
3. **会把方向搞反的地方**：`harmonise`、`Target_Intervention_Mismatch` 的 inhibit/activate 映射。
4. **UI 与引擎的接缝**：`plan_nodes` 的拓扑序、`build_inputs` 的路径可解析。

`tests/` 现有 47 个用例，覆盖 6 个模块。写法约定：
- 指标类断言**必须手算**（如 `AP = 1×0.5 + (2/3)×0.5`），不引用被测函数自证。
- 边界条件必须有：并列分数、decoy<3、零方差、样本量缺失、金标准与候选池无交集。
- 「诚实降级」也要测：`STEIGER_SKIPPED`、`NO_CONTROL`、`NOT_TRUE_LD_CLUMPING` 必须能被断言到。

### 3.2 CI 三层

```yaml
unit        # pytest，py3.10/3.11/3.12 矩阵
selftest    # 11 个模块的 --selftest（回归护栏）
integration # examples/ 构造数据真跑完整 DAG → 生成报告 → 断言 5 个关键产物非空
```

`integration` 是重点：单测全绿但流程崩是最常见的"提交前一天事故"。
断言用 `test -s`（存在且非空），而不是只检查退出码——因为曾出现过"节点返回 0 但产物只有表头"。

### 3.3 本地跑

```bash
pip install pytest numpy networkx pyyaml h5py
pytest tests -q
for m in np_mr np_synergy np_intervention np_benchmark np_lincs_local np_uireport; do
  python scripts/$m.py --selftest
done
```

---

## 4 审稿人可能追问的 6 个问题（附应答口径）

1. **"金标准是你自己整理的吗？"**
   → 是，规则与来源写在 JSON provenance；同时提供 `--gold-csv` 让审稿人换自己的金标准复算。
2. **"为什么不和 SwissTargetPrediction / BATMAN 比？"**
   → 需要同输入同金标准；我们给了 `--pred` 通用入口，任何工具的排序表都能进来，比较协议可复现。
3. **"AUC 0.7 算高吗？"**
   → 在基率 42% 的候选池里，AUC 0.7 对应的是"中等排序能力"；我们不宣称临床可用，
   只宣称"优于朴素基线且可回归"。
4. **"UI 会不会掩盖参数？"**
   → 每个节点命令都写进 `pipeline_manifest.json`，报告第 1 节逐条列出，可复制重跑。
5. **"报告里的水印是什么意思？"**
   → 仅表示本工具自查清单通过，不代表期刊接受或同行评议通过（`finalcheck_report.md` 末尾有声明）。
6. **"换台机器还能跑吗？"**
   → `integration` job 在干净 runner 上跑同一份 examples；`--resume` 有 sha256 复核；
   产物清单带哈希，可逐文件比对。

---

## 5 反模式清单

| 反模式 | 后果 | 正确做法 |
|---|---|---|
| 用测试集调完参再报测试集指标 | 审稿人一旦要求 hold-out 就崩 | 先注册基线（support/penalized），报告全部，不挑 |
| 金标准覆盖 80% 候选池还报 Recall@K | 指标失去分辨力，等于自欺 | 报 coverage + EF；改用边级指标 |
| 把 extended 靶点当负样本 | 所有方法 AUC 被压低，比较失真 | 从阴性池剔除（含糊集） |
| UI 里显示"预估进度"动画 | 一出错就不可信 | 用真实节点事件驱动 |
| 没有门禁也盖 Ready 水印 | 学术诚信风险 | 默认 UNVERIFIED，门禁决定 READY/DRAFT |
| 只跑单测不做集成 | 提交前一天流程崩 | CI 加 integration job |
