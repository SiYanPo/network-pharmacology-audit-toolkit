#!/usr/bin/env bash
# 非 NP 最小验证：同一套契约定义语言，不加修改，校验一张分子描述符表。
# 目的：把 "transferable" 从声明变成可检查的断言（Additional file 5）。
set -uo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}; command -v "$PY" >/dev/null 2>&1 || PY=python
SPEC=contract_descriptors.json
for f in descriptors.csv renamed_bad.csv type_bad.csv; do
  echo "== $f"
  "$PY" ../../scripts/np_contracts.py spec --file "$f" --spec "$SPEC"
  echo "   exit=$?"
done
