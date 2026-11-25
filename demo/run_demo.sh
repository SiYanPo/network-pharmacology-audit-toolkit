#!/usr/bin/env bash
# ============================================================================
# run_demo.sh — 5 分钟离线 demo：对一个真实案例跑"审计层"。
#
# 为什么 demo 只跑审计层而不跑全流程：
#   全流程需要 TCMSP/STRING/OpenTargets 等在线数据源，无网络环境必然失败；
#   而本 skill 的差异化恰恰在**审计层**（契约 + 证据账本 + SCM0061 评分 + 报告），
#   这一层纯标准库、零第三方依赖、离线可跑——正好用来证明"结论可被判读"。
#
# 为什么先复制到临时目录再跑：
#   np_scm0061.py 与 np_report.py 会把 scorecard / claims.json 写回 --dir。
#   若直接对 demo/case 跑，第二次运行时"冻结"的案例产物已被自己改过，
#   归档清单的逐文件 sha256 就对不上了（实测踩到过 3 个文件失配）。
#   因此 demo 一律在临时副本上跑，demo/case 保持只读语义。
#
# 预计耗时：< 30 秒（不含 python 启动）
# 用法： bash demo/run_demo.sh            # 在 skill 根目录执行
# ============================================================================
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
CASE=${1:-demo/case}
PY=${PY:-python3}
command -v "$PY" >/dev/null 2>&1 || PY=python

echo "== 0) 环境 =="
"$PY" -c "import sys;print('python',sys.version.split()[0])"

# 在临时副本上跑，避免把审计产物写回冻结的案例目录
WORK="${TMPDIR:-/tmp}/np_demo_case.$$"
rm -rf "$WORK"
mkdir -p "$WORK"
cp -r "$CASE"/. "$WORK"/
trap 'rm -rf "$WORK"' EXIT
echo "   工作副本：$WORK"

echo "== 1) 产物契约校验 =="
"$PY" scripts/np_contracts.py check --dir "$WORK" 2>&1 | tail -5
echo
echo "== 2) 证据账本（L0–L4 分布）=="
if [ -s "$WORK/09_evidence/evidence_ledger.csv" ]; then
  "$PY" - "$WORK" <<'PYEOF'
import csv,sys
from collections import Counter
rows=list(csv.DictReader(open(sys.argv[1]+"/09_evidence/evidence_ledger.csv",encoding="utf-8-sig")))
c=Counter((r.get("level") or r.get("evidence_level") or "?").strip() for r in rows)
print("  条目",len(rows),"| 等级分布",dict(sorted(c.items())))
PYEOF
else
  echo "  (该案例无证据账本，跳过)"
fi
echo
echo "== 3) SCM0061-2021 评价评分 =="
"$PY" scripts/np_scm0061.py --dir "$WORK" 2>&1 | head -6
echo
echo "== 4) 生成 HTML 报告（含评分卡章节）=="
OUT=${OUT:-/tmp/np_demo_report.html}
"$PY" scripts/np_report.py --dir "$WORK" -o "$OUT" 2>&1 | tail -2
echo "  报告：$OUT  ($(wc -c < "$OUT" 2>/dev/null || echo 0) 字节)"
echo
echo "== DEMO OK =="
echo "下一步（可选）：docker build -f Dockerfile.slim -t np-pro-slim:latest . && \\"
echo "                docker run --rm np-pro-slim:latest bash demo/run_demo.sh"
