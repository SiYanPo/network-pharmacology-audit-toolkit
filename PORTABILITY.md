# PORTABILITY —— 换机、交付与离线复现

> 目标：把 skill 从一个"作者本机能跑"的项目，变成"同事拷走能跑、审稿人能复现"的软件。
> 先跑 `python scripts/np_preflight.py`，它会把下面 90% 的问题一次性暴露。

---

## 1 支持矩阵（实测，非声明）

| 项目 | 支持 | 说明 |
|---|---|---|
| Windows 10/11 | ✅ 主力 | 当前开发与实测平台（路径含空格/中文已测） |
| Linux (Ubuntu 22.04+) | ✅ 核心 | CI 上跑通 examples 全流程；无 Windows 专有依赖 |
| macOS | ⚠️ 未实测 | 代码无平台专有调用，但未在 CI 覆盖 |
| Python 3.10–3.12 | ✅ 全功能 | 推荐区间（scanpy/anndata 可用） |
| Python 3.13 | ⚠️ 核心功能 | 实测 2026-11：scanpy 安装失败（numba/llvmlite 生态）→ 单细胞模块降级 |
| Python ≤3.9 | ❌ | 使用了 `math.comb`、`str.removeprefix` 等 3.9+ 特性 |

---

## 2 三分钟部署

```bash
git clone <repo> && cd network-pharmacology-pro
python scripts/np_preflight.py          # 先体检：26/36 通过属正常（工具缺失是 WARN）
pip install -r requirements-lock.txt    # 可选增强（核心零第三方依赖）
pip install streamlit                   # 需要 Web UI 时
streamlit run app.py                    # 或走命令行：python scripts/np_pipeline.py --help
```

**不要**直接 `pip install -r requirements.txt`：那是"可选增强清单"，没有版本锁。
`requirements-lock.txt` 是实测通过的版本组合。

---

## 3 环境变量（全部可选，但推荐固定）

| 变量 | 作用 | 默认 |
|---|---|---|
| `NP_CACHE` | 跨运行缓存根目录（STRING/HPA/ortholog） | `~/.np_cache` |
| `NP_OFFLINE` | `1` = 禁止一切联网（STRING/HPA/Ensembl/MR fetch） | 未设置 |
| `NP_PLINK` | plink 可执行文件路径 | 自动探测（含 `D:/MDtools/plink`） |
| `NP_BROWSER` | HTML→PDF 用的 Chrome/Edge 路径 | 自动探测 |
| `NP_WORKDIR` | Streamlit 输出根目录 | `<skill>/np_ui_runs` |

**缓存为什么必须共享**：把缓存放在输出目录里，每次运行都要重新联网；
实测同一份输入冷启动 **153 s / 423 s**（网络波动 2.8×），统一缓存后 **28 s**，且产物哈希完全一致。

---

## 4 外部工具（按需安装）

核心流程（标准化 → 交集 → PPI → 富集 → 报告）**不需要任何外部二进制**。以下按功能可选：

| 工具 | 用途 | 官方来源 | 本机实测 |
|---|---|---|---|
| PLINK 1.9 | MR 真 LD clumping | https://www.cog-genomics.org/plink/1.9/ | ✅ v1.90b7（`D:/MDtools/plink/plink.exe`），clump 实测 6→2 |
| AutoDock Vina | 分子对接 | https://github.com/ccsb-scripps/AutoDock-Vina | ✅ v1.2.5（`D:/MDtools/vina.exe`） |
| Open Babel | 格式转换/配体准备 | https://openbabel.org/ | ✅（conda 环境内） |
| acpype | 配体 GAFF 拓扑 | `pip install acpype` | ✅ 2026.10.15（当前解释器 `Scripts/acpype.exe`） |
| GROMACS | MD 模拟 | https://manual.gromacs.org/ | ❌ 未装（WSL `apt install gromacs`） |
| AmberTools | antechamber 电荷 | https://ambermd.org/AmberTools.php | ❌ 未装 |
| BLAST+ | 引物特异性 | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ | ❌ 未装 |
| isPcr | 基因组特异性 PCR | http://users.soe.ucsc.edu/~kent/src/ | ❌ 未装 |
| fpocket | 口袋预测 | https://github.com/Discngine/fpocket | ❌ 未装（有 Windows 构建） |
| GNINA | CNN 打分对接 | https://github.com/gnina/gnina | ❌ 未装（需 CUDA） |

安装脚本：`bash scripts/setup_env.sh`（WSL/conda 路径），或按上表手动装。
装完后 `python scripts/np_envcheck.py` 应显示 `[OK ]`。

---

## 5 离线复现（审稿人场景）

```bash
# 1) 联网跑一遍，缓存自动落在 ~/.np_cache
python scripts/np_pipeline.py --inputs inputs.json -o np_out

# 2) 打包快照（缓存 + 输入 + manifest）
tar czf snapshot_$(date +%Y%m%d).tar.gz ~/.np_cache np_out/pipeline_manifest.json inputs.json

# 3) 另一台机器：解开快照，指向它并禁网
export NP_CACHE=$PWD/cache NP_OFFLINE=1
python scripts/np_preflight.py --strict          # 应只剩工具类 WARN
python scripts/np_pipeline.py --inputs inputs.json -o np_out2
python scripts/np_benchloop.py verify --a np_out --b np_out2   # 内容级比对（忽略时间戳）
```

**当前缺口（如实记录）**：缓存文件尚未写 `release/etag` 元数据，所以"数据库版本"只能从 URL 推断。
Phase 1 计划给每个缓存文件配 `*.meta.json`（url/release/fetched_at/sha256/etag），
manifest 引用它——这是把"可复现"从口号变成证据的最小改动。

---

## 6 已知不可移植点（换机必看）

1. `D:/MDtools/...`：本机工具目录，换机后由 `NP_*` 环境变量或 `np_envcheck.KNOWN` 覆盖。
2. WSL 路径（`/opt/ambertools/bin/antechamber`）：仅 Linux/WSL 有效；Windows 原生需另装。
3. Python 3.13 + scanpy：装不上（生态问题，不是代码问题）→ 单细胞模块会 `TOOL_MISSING` 降级。
4. `~/.np_cache`：不随仓库分发；离线复现必须显式打包（第 5 节）。
5. 中文字体：报告用系统字体栈；PDF 由浏览器渲染，Linux 无中文字体时会缺字（需装 `fonts-noto-cjk`）。
6. **h5py/HDF5 在 Windows 上偶发 `HDF5: infinite loop closing library`**（本轮实测，可导致解释器卡死）。
   规避：`np_preflight` 已改为在子进程里探测依赖；业务代码里 h5py 只开一个文件句柄、不要并发打开同一文件；
   若反复出现，升级/重装 `h5py`（`pip install --force-reinstall h5py`）。

---

## 7 换机检查清单（打印版）

- [ ] `python scripts/np_preflight.py` 无 FAIL（WARN 可接受）
- [ ] `python scripts/np_pipeline.py --inputs inputs.json --dry-run` 能列出 11 个节点
- [ ] `python scripts/np_finalcheck.py --selftest` → 22 gates PASS
- [ ] `pytest tests -q` → 全绿
- [ ] `streamlit run app.py` 能打开四个页签
- [ ] 路径含空格/中文的目录能正常输出（已内置探针）
- [ ] 需要 MR 时：`plink --version` 可用，且 `--bfile` 指向 1000G 子集
- [ ] 需要 PDF 时：`python scripts/np_pdf.py --check` 找到 Edge/Chrome

---

## 8 版本与升级

- 语义化版本写在 `CHANGELOG.md`；破坏性变更进 `MAJOR`（如产物列名变更）。
- 升级后先跑 `np_selfaudit.py run --run-selftests`，确认没有模块从 PASS 变 FAIL。
- 论文引用请固定 commit hash + `snapshot_YYYYMMDD.tar.gz` 的 sha256。

## v3.0 追加：工具解析、快照、容器（2026-11-22）

### 1. 外部二进制只走一个入口

```bash
python scripts/np_tools.py list        # 每个工具的 路径 / 来源 / 是否可用
python scripts/np_tools.py resolve vina
```

解析顺序（**任何模块都不得再写死机器路径**，`tests/test_tool_registry.py` 会 grep 拦截）：

1. 显式参数（`--vina`）
2. 环境变量：`NP_VINA` / `NP_OBABEL` / `NP_GMX` / `NP_PLINK` / `NP_FPOCKET` / `NP_GNINA` / `NP_ACPYPE` / `NP_ANTECHAMBER` / `NP_MD_PYTHON` / `NP_PDBFIXER`
3. `config/tools.local.json`：`{"vina": "C:/tools/vina.exe"}`（本机配置，不进版本库）
4. `PATH`
5. 旧默认根：`$NP_TOOL_ROOT` 或 Windows 的 `D:/MDtools`（仅兼容）

缺工具时的报错会直接给三种修法，而不是抛 `FileNotFoundError`。

### 2. 版本与离线快照

```bash
python scripts/np_snapshot.py probe --save          # 记 STRING/Ensembl/MyGene 版本
python scripts/np_snapshot.py pack -o snap.tar.gz   # 打包缓存 + manifest（逐文件 sha256）
python scripts/np_snapshot.py verify snap.tar.gz    # 校验（篡改会报 CORRUPT）
python scripts/np_snapshot.py mount snap.tar.gz --to D:/snap
source D:/snap/activate_snapshot.sh                 # 或 .bat；设 NP_CACHE/NP_SNAPSHOT
```

实测：挂载后重跑流水线，`intersection_genes.txt` / `hubs.csv` / `enrich_all.csv`
三个产物 SHA256 **逐字节相同**。

已知不可自动获取的版本：**HPA**。下载时写 `release: unrecorded`，
需要用户在 Methods 里补 `NP_HPA_RELEASE`（例如 24.0）后重跑 `fetch_bulk`。

### 3. 容器

```bash
docker build -t np-pro:latest .
docker run --rm -v $PWD:/work -w /work np-pro:latest \
  python scripts/np_pipeline.py --inputs examples/inputs.demo.json -o /tmp/o --offline
```

Dockerfile 在**构建期**就跑一遍离线流水线 + 契约 + FDR，避免"工具装了但流程跑不通"。
本机（2026-11-22）无 docker，构建验证由 CI job `docker-demo` 承担——**本地未实测**，如实标注。

### 4. 覆盖率棘轮

```bash
python -m coverage run --source=scripts -m pytest tests -q
python -m coverage json -o coverage.json
python tools/check_coverage.py --coverage coverage.json --baseline coverage_baseline.txt
```

当前基线 15.22%（大量逻辑在 `main()` 里）。**棘轮语义**：只允许上调，不允许回落。
不设 `--cov-fail-under=55` 是因为它会让 CI 永远红，最后一定被人关掉。
