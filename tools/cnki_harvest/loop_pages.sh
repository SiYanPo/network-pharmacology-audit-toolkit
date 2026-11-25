#!/usr/bin/env bash
# ============================================================================
# loop_pages.sh <start> <end> — 逐页抽 CNKI 结果表并翻页
# 依赖：bsk daemon 已在运行（BSK_AUTO_START=0），当前标签停在该检索的结果第 1 页
# 每页：extract discover 取 xt_* → extract table 落 CSV → click a#PageNext
# 断点续跑：已存在的 pN.csv 且行数>1 则跳过抽取，但仍需翻页（用 --resume-skip）
# ============================================================================
set -uo pipefail
BSK=/c/WorkBuddy-Home2/.local/bin/bsk.exe
OUT=/c/WorkBuddy-Home2/npwork/cnki
S=${1:-1}; E=${2:-13}
mkdir -p "$OUT"
export BSK_AUTO_START=0

for i in $(seq "$S" "$E"); do
  # 若目标页已存在，直接翻页跳过
  if [ -s "$OUT/p$i.csv" ] && [ "$(wc -l < "$OUT/p$i.csv")" -gt 1 ]; then
    echo "page $i exists, skip extract"
  else
    timeout 45 "$BSK" extract discover --json > "$OUT/disc_$i.json" 2>&1
    TID=$(python /c/WorkBuddy-Home2/npwork/parse_disc.py "$OUT/disc_$i.json" table 2>/dev/null | head -1 | cut -f2)
    if [ -z "$TID" ]; then
      echo "page $i: NO_TABLE (可能被验证拦截或未加载)"; break
    fi
    timeout 60 "$BSK" extract table --target "$TID" --format csv --out "$OUT/p$i.csv" --overwrite --max-rows 60 >/dev/null 2>&1
    ROWS=$(python /c/WorkBuddy-Home2/npwork/count_rows.py "$OUT/p$i.csv")
    echo "page $i ok rows=$ROWS"
    if [ "$ROWS" -le 0 ]; then echo "page $i: EMPTY"; break; fi
  fi
  # 翻页
  timeout 30 "$BSK" click --selector "a#PageNext" >/dev/null 2>&1 || { echo "page $i: click next failed"; break; }
  sleep 2
done
echo "done $S-$E; files:"; ls "$OUT" | grep -c '^p[0-9]*\.csv$'
