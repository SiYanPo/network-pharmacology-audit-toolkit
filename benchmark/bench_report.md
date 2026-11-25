# 回归基准报告

- 时间：2026-11-20 17:30
- 用例：5，通过 5
- runner：`python evals/integrity.py --input {input} --out {out}`

| 用例 | 类别 | 结果 | 说明 |
|---|---|---|---|
| empty-input | empty | PASS | 空药名表 |
| malformed-row | malformed | PASS | 多余逗号 |
| nonstandard-name | nonstandard | PASS | 矿物药/拉丁名混写 |
| typical | typical | PASS | 典型药对 |
| idempotent-rerun | regression | PASS | 两次一致 |