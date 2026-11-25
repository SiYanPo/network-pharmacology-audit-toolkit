#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_archive.py — confirm that an extracted tree matches the archive manifest.

Usage:
    python verify_archive.py <ARCHIVE_MANIFEST.tsv> <extracted-root>

The manifest is shipped next to the archive (not inside it) and holds one row per file:
`sha256 <TAB> bytes <TAB> <top-level-dir>/<relative-path>`.

Exit code 0 = every file matches; 1 = at least one mismatch; 2 = manifest missing.
"""
import csv
import hashlib
import os
import sys


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    man, root = sys.argv[1], sys.argv[2]
    if not os.path.exists(man):
        print(f"manifest not found: {man}")
        return 2
    bad = missing = n = 0
    for r in csv.DictReader(open(man, encoding="utf-8"), delimiter="\t"):
        rel = r["path"].split("/", 1)[1]
        p = os.path.join(root, rel)
        if not os.path.exists(p):
            missing += 1
            print("MISSING   ", rel)
            continue
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        n += 1
        if h != r["sha256"]:
            bad += 1
            print("MISMATCH  ", rel)
    print(f"checked {n} files: {bad} mismatched, {missing} missing")
    if bad or missing:
        return 1
    print("verified OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
