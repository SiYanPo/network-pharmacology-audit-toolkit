#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_prov.py — 数据库版本快照与可复现性日志（provenance）

每次复现自动记录：每个外部数据源的 访问日期 / 版本号 / 请求参数 / 端点 / 重试次数 / 结果行数。
当结果与原文不一致时，可用本报告精确归因于「数据库版本漂移」而非算法错误。

子命令
  snapshot   从 np_fetch 注册表 + 用户提供的版本信息生成 provenance.json
  diff       对比两次快照，输出漂移的数据源和字段
  stamp      给结果目录打时间戳+SHA256 清单（FAIR 审计）
"""
from __future__ import annotations

import argparse, datetime as dt, hashlib, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REG = os.path.join(ROOT, "data", "fetch_cache", "source_registry.json")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cmd_snapshot(args):
    reg = json.load(open(REG, encoding="utf-8")) if os.path.exists(REG) else {}
    versions = json.load(open(args.versions, encoding="utf-8")) if args.versions else {}
    snap = {
        "schema": "np-provenance/1.0",
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "skill_version": open(os.path.join(ROOT, "SKILL.md"), encoding="utf-8").read().split("version: ")[1].split("\n")[0].strip(),
        "sources": {}
    }
    for name, meta in reg.items():
        v = versions.get(name, {})
        snap["sources"][name] = {
            "url": meta.get("url"),
            "access_mode": meta.get("status"),
            "probe_status": meta.get("status_at_probe"),
            "version": v.get("version", "未记录"),
            "access_date": v.get("date", snap["generated"][:10]),
            "request_params": v.get("params", {}),
            "n_records": v.get("n_records", "未记录"),
            "register_url": meta.get("register"),
            "note": meta.get("note")
        }
    out = args.output or os.path.join(args.outdir or ".", "provenance.json")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    json.dump(snap, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"溯源快照 → {out}")
    print(f"记录数据源 {len(snap['sources'])} 个；skill {snap['skill_version']}")
    miss = [n for n, s in snap["sources"].items() if s["version"] == "未记录"]
    if miss:
        print(f"⚠ 未记录版本：{', '.join(miss)}（建议用 --versions 补，否则无法归因版本漂移）")


def cmd_diff(args):
    a = json.load(open(args.left, encoding="utf-8"))["sources"]
    b = json.load(open(args.right, encoding="utf-8"))["sources"]
    print("=== 数据库漂移对比 ===")
    drift = 0
    for name in sorted(set(a) | set(b)):
        x, y = a.get(name, {}), b.get(name, {})
        for field in ("version", "access_date", "probe_status", "n_records"):
            if x.get(field) != y.get(field):
                print(f"  {name}.{field}: {x.get(field)!r} → {y.get(field)!r}")
                drift += 1
    print(f"共 {drift} 处漂移" if drift else "无漂移")


def cmd_stamp(args):
    items = []
    for root, dirs, files in os.walk(args.dir):
        dirs[:] = [d for d in dirs if d not in {".git", "__pycache__"}]
        for fn in files:
            p = os.path.join(root, fn)
            rel = os.path.relpath(p, args.dir)
            items.append({"path": rel, "size": os.path.getsize(p), "sha256": sha256(p)})
    out = os.path.join(args.dir, "MANIFEST.sha256.json")
    json.dump({"generated": dt.datetime.now().isoformat(timespec="seconds"),
               "files": items}, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"已签名 {len(items)} 个文件 → {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("snapshot"); s.add_argument("--versions"); s.add_argument("--outdir"); s.add_argument("--output"); s.set_defaults(fn=cmd_snapshot)
    d = sub.add_parser("diff"); d.add_argument("left"); d.add_argument("right"); d.set_defaults(fn=cmd_diff)
    t = sub.add_parser("stamp"); t.add_argument("dir"); t.set_defaults(fn=cmd_stamp)
    args = p.parse_args(); args.fn(args)


if __name__ == "__main__":
    main()
