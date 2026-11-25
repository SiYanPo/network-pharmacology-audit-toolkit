# 标准表契约（NP 技能）

> 2026-11-22 校正：本文件此前写 `standard_*.csv`，而 `np_export.py` 实际产出
> `herbs.csv / compounds.csv / targets_mapped.csv / PPI_network.sif / enrichment_results.csv`。
> 契约测试 `tests/test_contracts.py` 现在逐列比对，二者不一致即失败。
> 改列名 = 破坏性变更（bump major）。

## herbs.csv
`herb, source_db, matched, fallback_strategy`

## compounds.csv
`compound, source_db, filter_rule, ob, dl, kept`

## targets_mapped.csv
`raw, uniprot_id, symbol, organism, layer`

## PPI_network.sif
`node1, interaction, node2, score`

## enrichment_results.csv
`term, source, pvalue, qvalue, n_hit, term_size, term_id`

规则：UTF-8-SIG；缺失留空并写 meta.missing；改列名 = 破坏性变更（bump major）。
