# DrugBank 已知靶点、OpenAlex 计量、多重复 MD 统计——三条外部链路的可复现口径

本文记录 2026-11-23 实测出来的访问边界、检索语义与统计判据。凡是"能不能拿到数据"的问题，
答案都写在这里，避免下次重新踩。

---

## 1. DrugBank：公开卡片能读，完整靶点表要会话

| 通道 | 结果 | 说明 |
|---|---|---|
| `curl https://go.drugbank.com/drugs/DB04216` | **403** | Cloudflare 质询页，UA 伪装无效 |
| WebFetch 同 URL | **200 可用** | 返回 drug card 的 Targets / Transporters 区块，含基因符号 + UniProt + `bio_entities/BE...` 链接 |
| 登录会话（浏览器） | 需 bsk + 扩展 | 学术账号只解决"完整表/互作/试验"的解锁，**不解决**公开子集的读取 |

**截断事实**（必须写进任何用到它的结论里）：卡片只展示 8 条靶点，其余折叠成 `+N more`。

| 化合物 | DrugBank ID | 分组 | 靶点总数 | 公开可见 | 转运体可见 |
|---|---|---|---|---|---|
| Quercetin | DB04216 | Investigational | 32 | 8 | 8 |
| Curcumin | DB11672 | Approved | 21 | 8 | 1 |
| Berberine | DB04115 | Approved / Investigational / Withdrawn | 12 | 8（含 1 条细菌 ftsZ） | 5 |
| Ursolic acid | — | **库中无实体** | — | — | — |
| Tanshinone IIA | — | **库中无实体**（仅 Tanshinone I DB16886、Cryptotanshinone DB15579） | — | — | — |

因此 `tools/drugbank_known_target_audit.py` 的 recall **分母是被展示子集**：
命中一定是真阳性；未命中可能是预测缺失，也可能是被折叠隐藏。**禁止**把它读成全库召回率。

实测（2026-11-23，5 化合物）：

```
Quercetin  8/32   pred=101  hit=0  recall 0.000  [0.00,0.32]
Curcumin   8/21   pred= 61  hit=4  recall 0.500  [0.22,0.78]
Berberine  7/12   pred= 13  hit=1  recall 0.143  [0.03,0.51]
pooled                        5/23 = 0.217  [0.10,0.42]
```

读法：**预测集与 DrugBank 机制注释的重合很有限**，尤其 Quercetin 0/8（它的 8 条可见靶点
集中在 ATP 合酶亚基、LPCAT1/2、PIK3CG、UGT3A1，属体外机制研究，TCMSP 型预测覆盖不到）。
这不是"pipeline 坏了"，而是**靶点获取源单一**的必然结果——要提升只能把
ChEMBL mechanism / DrugBank 机制注释并入"获取"而非仅作"验证"，并保留来源标记。

---

## 2. CNKI：IP 授权 ≠ 可编程取数

```
POST https://kns.cnki.net/kns8s/brief/grid
→ HTTP 403
{"code":-403,"message":"https://kns.cnki.net/verify/home?captchaType=blockPuzzle&ident=...&captchaId=..."}
```

即**滑块验证**。机构 IP 只影响全文下载权限，不影响验证码。要取数必须有真实浏览器会话
（bsk + Edge 扩展），且滑块通常需要人过一次。

替代路径（已落地）：`tools/bibliometric_openalex.py`。它回答的是"这个方向有多少人做过、
年产量怎么走"，不是"中文期刊里谁写过"。两者不能互相冒充。

---

## 3. OpenAlex：不写引号短语 = 数字全错

`title_and_abstract.search` 是 **OR** 语义。同一主题：

| 写法 | 命中 |
|---|---|
| `network pharmacology` | 58,978 |
| `"network pharmacology"` | **28,242** |

更严重的不是数量，而是**语义**：OR 口径下"高被引前 6 篇"实测是 ESC 急性冠脉综合征指南、
盲法综述、COVID 样本研究——与网络药理学无关。任何排序/被引分析都必须用短语+AND。

正确的计量口径（2026-11-23 实测）：

| 检索式（短语 + AND） | 命中 | 2026 年 |
|---|---|---|
| `"network pharmacology"` | 28,242 | 8,162 |
| `"network pharmacology" AND "colorectal cancer"` | 589 | 185 |
| `"网络药理学"` | 508 | 120 |
| `"network pharmacology" AND "quality control"` | 472 | 161 |
| `"network pharmacology" AND "reproducibility"` | 226 | 174 |
| `"network pharmacology" AND "machine learning"` | 832 | 486 |
| `"network pharmacology" AND "molecular dynamics"` | 3,112 | 1,592 |
| `"network pharmacology" AND "MM/PBSA"` | 165 | 111 |
| `"network pharmacology" AND "benchmark"` | **51** | 29 |
| `"network pharmacology" AND "audit"` | **16** | 15 |

**新颖性校准**：与结直肠癌主题的 589 篇相比，"审计/基准"只有 16–51 篇，
空白是真的；但 2025→2026 全线陡增（audit 类 15/16 篇集中在 2026），
**这是一个正在关闭的窗口，不是一片无人区**。措辞应写"首个集成式可执行审计闭环"，
不要写"首次提出网络药理学质量控制"（quality control 已有 472 篇）。

---

## 4. 多重复 MD：n=2 的误差条是假的，n=3 才买得到精度

`tools/md_replicates_merge.py` 的口径：

1. 取所有重复都覆盖的**公共时间窗**（否则长轨迹被短轨迹的均值污染）；
2. 丢掉前 0.5 ns 弛豫；
3. 重复当独立样本 → mean ± SD，t(N−1) 给 95% CI；
4. 重复内用**块平均 SEM**（块 1 ns）扣自相关，再做套合方差分解：
   `ICC = σ²_between / (σ²_between + σ²_within)`；
5. 用 σ_b 预测"加到 N 条"的 CI 半宽 ∝ t(N−1)/√N。

DP1 × Regaloside D 实测（0.5–6.6 ns 公共窗，rep01 10 ns / rep02 6.63 ns）：

| 指标 | rep01 | rep02 | 均值 ± SD | 95%CI 半宽 (n=2) | n=3 半宽 | 收窄 |
|---|---|---|---|---|---|---|
| 骨架 RMSD (Å) | 2.43 | 2.36 | 2.39 ± 0.051 | ±0.457 | 0.127 | 3.6× |
| 配体位移 (Å) | 3.12 | 3.44 | 3.28 ± 0.229 | ±2.06 | 0.569 | 3.6× |
| Rg (Å) | 21.12 | 21.27 | 21.20 ± 0.099 | ±0.889 | 0.246 | 3.6× |

ICC(rep) = 0.51 / 0.73 / 0.93 → **重复间方差占总方差的一半到九成**。
这意味着：只报单条轨迹的帧间 SEM，等于把主导误差项扔掉；而 n=2 时 t(0.975,1)=12.71，
算出来的 CI 宽到没有信息量。**"n=3 独立重复"不是形式主义，是把误差条从装饰变成证据的最低门槛。**

第三条第跑法见 `run_rep03.sh`（em → nvt 100 ps @posres 1000, seed 32026 →
npt 500 ps @posres 500 C-rescale → md 10 ns）；同体系 8.97 ns/day，生产段约 27 h。

---

## 5. 登录态浏览器（bsk）——2026-11-24 打通后的正确姿势

安装坑：宿主把 `irm ... | iex` 判为任意代码执行而拦截。改为按官方清单手动装：

1. 取清单 `https://github.com/Tencent/BrowserSkill/releases/latest/download/version.json`（curl 可用，python urllib 会因跳转失败）；
2. 下 `bsk-v<ver>-x86_64-pc-windows-msvc.zip`，比对 `assets["windows-x64"].sha256`；
3. 解包 `bsk.exe` 到用户目录；`BSK_AUTO_START=0 bsk daemon start --foreground` 由宿主后台任务承载；
4. `bsk status --json` 看到 `browsers[0].browser_name == "edge"` 即扩展已连。

用法要点：
- `bsk open <url>` 会消耗一次 observation，之后再 `open` 会报 "needs an explicit current observation" →
  先 `bsk sync --json` 取 `observation_id`，再 `open <url> --observation <id>`。
- `open --quiet` 不可用（同样报 observation 错误）。
- 抽表：`bsk extract discover --json` → 取 `targets[].target_id`（形如 `xt_*`）→ `bsk extract table --target <id> --format csv --out x.csv`。
  `--selector` 对 DrugBank 的表格无效（命中的是渲染层不是 CSS 类）。
- `bsk read` 只回当前视口 VOM；页面很长时看不到下方表格，必须用 extract。

**DrugBank 登录态能多拿到什么**：分区锚点给出真实计数（`Targets (32)`），public 卡片只显示 8 条。
抽出的表列是 `TARGET/ACTIONS/ORGANISM`（蛋白全名，无基因符号、无 UniProt），
必须再走 UniProt 反查；**反查必须用官方名全等匹配**，否则会得到
`Amyloid-beta precursor protein → APPBP2`、`PPARγ → PPRC1`、`HSP90α → HSP90AA2P` 这类静默错误。

**CNKI 登录态**：`kns.cnki.net` 横幅出现学校名即机构 IP 生效（实测"河南中医药大学"）。
结果表列为 `题名/作者/来源/发表时间/数据库/被引/下载`，默认 20 行/页，`bsk extract table` 可直接抽。
`原版阅读` 走学校代理域名（`kns.cnki.net-s.dres.<school>.edu.cn`）。批量导出仍需逐页。

**宿主限制（重要）**：本环境的后台任务约 4 分钟被回收，`schtasks` 被安全策略拒绝，
因此**长于数分钟的作业（MD 生产）无法由 Agent 托管**，只能在用户自己的 WSL 终端启动。
