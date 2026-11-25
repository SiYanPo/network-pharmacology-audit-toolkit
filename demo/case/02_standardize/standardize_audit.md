# 靶点标准化审计

- 输入：`D:/WorkBuddy/分身文件/2026-11-19-10-40-46\全流程实测\inputs\Quercetin\compound_targets.csv`
- 记录数：109，去重后基因：101，映射率（去重基因/去重原始名）：92.7%

## 命中层级分布

| 层级 | 条数 | 含义 |
|---|---|---|
| L1_exact | 98 | 蛋白名精确命中 |
| unmatched | 8 | 未命中 |
| L2_normalized | 3 | 规范化后命中 |

## 未命中条目（前 40）

- `Casein kinase II` → unmatched（未在内置词典命中，需人工确认或联网映射）
- `Proteasome Macropain subunit MB1` → unmatched（未在内置词典命中，需人工确认或联网映射）
- `Glycogen synthase kinase-3` → unmatched（未在内置词典命中，需人工确认或联网映射）
- `METTL3/METTL14` → unmatched（未在内置词典命中，需人工确认或联网映射）
- `Aldo-keto reductase family 1 member C21` → unmatched（未在内置词典命中，需人工确认或联网映射）
- `D(4) dopamine receptor` → unmatched（未在内置词典命中，需人工确认或联网映射）
- `Alpha amylase` → unmatched（未在内置词典命中，需人工确认或联网映射）
- `Menin/Histone-lysine N-methyltransferase MLL` → unmatched（未在内置词典命中，需人工确认或联网映射）

> 处理建议：非人类/噪音条目直接剔除；疑似真人类靶点补入 `--alias` 别名表后重跑。
> **不得**为了凑数把非人类蛋白强行映射到人类基因——这会造成系统性假阳性。
