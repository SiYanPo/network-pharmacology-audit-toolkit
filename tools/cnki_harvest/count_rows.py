#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""count_rows.py <csv> — 统计 CNKI 结果 CSV 的数据行数（跳过表头，容忍单元格内换行）。"""
import csv
import sys

try:
    rows = list(csv.reader(open(sys.argv[1], encoding="utf-8-sig")))
    print(max(0, len(rows) - 1))
except Exception:
    print(0)
