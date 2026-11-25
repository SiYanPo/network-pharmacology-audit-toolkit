#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""parse_disc.py — 从 `bsk extract discover --json` 的输出里取出表格名与 target_id。

bsk 会在 JSON 之后追加日志行，所以不能直接 json.load；用 raw_decode 取第一个完整对象。
用法： python parse_disc.py <disc.json> [表格名关键字]
"""
import json
import sys

p = sys.argv[1]
kw = sys.argv[2] if len(sys.argv) > 2 else ""
raw = open(p, encoding="utf-8", errors="ignore").read()
i = raw.find("{")
d = json.JSONDecoder().raw_decode(raw[i:])[0]
for t in d.get("targets", []):
    name = t.get("name", "")
    if kw and kw not in name:
        continue
    print(f"{name}\t{t['target_id']}\t{','.join(t.get('columns', []))[:70]}")
