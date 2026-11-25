# 投稿用：Code Availability 与 Data Availability 声明（可直接粘贴）

以下文字按"事实陈述、不越界承诺"的原则起草，可直接用于投稿系统或正文末节。
方括号 `[...]` 内为需作者填写的内容。

---

## 一、Code Availability（推荐版）

> The complete codebase of `network-pharmacology-pro`, including all audit-layer modules
> (evidence ledger, product contracts, SCM0061-2021 scorecard, wet-lab result grading),
> configuration files, the bundled example case, and the regression test suite, is available at
> **[repository URL]** under **[license, e.g. MIT]**.
>
> To make the audit layer independently verifiable without installing any third-party package,
> we provide a slim container image definition (`Dockerfile.slim`, standard library only) and a
> demo script that reproduces the key numbers reported in this work:
>
> ```bash
> docker build -f Dockerfile.slim --build-arg BASE_IMAGE=docker.m.daocloud.io/library/python:3.13-slim -t np-pro-slim .
> docker run --rm np-pro-slim:latest
> ```
>
> Expected output on the bundled case: `[contracts] PASS: 26 contracts, 0 hard failures`;
> `SCM0061-2021 score 26/34 = 76.5% (grade B)`; `DEMO OK`. A step-by-step verification
> checklist mapping each claim in this paper to an executable command is provided in
> `REPRODUCIBILITY.md`. The regression suite reports `121 passed, 1 skipped`.
>
> The full pipeline (which additionally requires online data sources such as TCMSP, STRING and
> Open Targets) is not exercised by the offline demo; the container `Dockerfile` at the
> repository root specifies the complete environment.

**写作要点**：
- 明确区分"审计层可离线核验"与"完整流水线需联网"，不要让人误以为 demo 跑的是全流程；
- 给出**预期输出字符串**，这是最有力的可复现承诺（对方一跑就知道有没有被改过）；
- 不要写"fully reproducible"这类无限定词。

---

## 二、Data Availability（推荐版）

> All intermediate products of the case study (79 files for the *Quercetin* case, including
> standardized target tables, intersection sets, PPI topology, enrichment tables, the L0–L4
> evidence ledger with provenance hashes, and the SCM0061-2021 scorecard) are provided in
> `demo/case/` of the repository and are reproduced verbatim in the container demo.
>
> The retrospective analysis used in Section [x] is a secondary analysis of a reproduction
> dataset of 500 published network pharmacology papers; the derived table
> (`06_repro500.csv`, 19 columns) and the analysis scripts (`retro_gate_analysis.py`,
> `retro_sweep.py`) are provided under `全流程实测/retro/` together with the JSON outputs
> containing all effect sizes and confidence intervals.
>
> Public data sources used: TCMSP, STRING, Open Targets, DrugBank (academic access, accessed
> [date]), Europe PMC / PubMed (via REST API, queried [date]). Third-party data are not
> redistributed; retrieval parameters and database versions are recorded in each stage's
> `environment.json`.

---

## 三、必须同时声明的三点（建议放 Limitations 或 Ethics 段）

1. **靶点数一致 ≠ 机制正确。** 本工作的结局指标是"复现的靶点数量是否落在原文 ±20% 内"，
   不是"机制结论是否正确"。所有计算预测均不构成实验证据。
2. **回顾性分析非前瞻注册试验。** 门禁阈值在看到结局前由脚本固定（`retro_gate_analysis.py`
   头部注释可查证），但研究设计本身未注册；382/500 篇因原文未报告靶点数被排除。
3. **与已发表平台的对比是流程级的，不是工具级的。** 我们尝试运行 NeXus v1.2 作同输入对照，
   但按其论文给出的仓库地址（截至 [date]，经三条独立镜像重试）未获取到源代码；
   因此本文的对照基线为"传统四步流程"，由本代码库的模块串联而成，
   **不代表**对 NeXus 或其他第三方工具的性能评价。

---

## 四、不建议写的表述（会被审稿人抓）

| 不建议 | 原因 | 替代写法 |
|---|---|---|
| "our framework is fully reproducible" | 无限定词 | "the audit layer is independently verifiable offline; the full pipeline requires online data sources" |
| "first quality-control framework for network pharmacology" | SCM0061-2021 在先 | "first executable implementation of the SCM0061-2021 evaluation requirements" |
| "outperforms existing platforms" | 无第三方实跑证据 | "adds audit and uncertainty quantification on top of conventional workflows" |
| "NeXus did not release code" | 属推断 | "we were unable to obtain the source code from the stated repository; we cannot exclude other forms of availability" |
