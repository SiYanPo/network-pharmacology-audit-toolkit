#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_snapshot.py — 版本化快照与离线复现（v3.0）

为什么需要（评审报告 §5，七寸之五）
------------------------------------
审稿人问"结果可追溯到哪个数据库版本"时，只记 URL 是不够的：
STRING 有 v12.0、Ensembl 有 release 112、HPA 有 release 24.x。没有版本号的 Methods
只能写"latest"，而"latest"在投稿周期内会变，等于不可复现。

本模块做三件事
--------------
1. **probe**：联网探测 STRING / Ensembl 的当前 release，并把每个缓存文件的
   `{url, release, fetched_at, sha256, etag}` 写成同名 `*.meta.json`。
2. **pack**：把整棵缓存目录打成 `snapshot_YYYYMMDD.tar.gz`，内含 manifest（版本 + 哈希）。
   随论文/仓库发布，几十 MB，可审稿人自证。
3. **mount**：把快照解到本地并把 `NP_CACHE` 指过去；所有模块经 `npcore.cache_dir()`
   自动读它，从而实现"无网复现"。

诚实的边界
----------
* HPA 的 release 若下载时未记录，只能写 `unrecorded`——**不猜测**。
* 快照只覆盖"缓存过的外部响应"，不覆盖用户输入数据（那些本来就要随论文提供）。
* 若某模块在无缓存时会回退全基因组/离线降级，快照里必须有对应文件，否则仍会降级；
  `verify` 会列出缺失项。

用法
----
    python np_snapshot.py probe                 # 探测在线版本 + 补齐缓存 meta
    python np_snapshot.py pack -o snapshot_20261122.tar.gz
    python np_snapshot.py mount snapshot_20261122.tar.gz --to D:/np_snap
    python np_snapshot.py verify snapshot_20261122.tar.gz
    python np_snapshot.py selftest
"""

import argparse
import glob
import hashlib
import json
import os
import sys
import tarfile
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

VERSION = "3.0"
UA = {"User-Agent": "network-pharmacology-pro/3.0 (snapshot probe)"}


def cache_root():
    env = os.environ.get("NP_SNAPSHOT") or os.environ.get("NP_CACHE")
    if env and os.path.isdir(env):
        return env
    return os.path.join(os.path.expanduser("~"), ".np_cache")


def _get_json(url, timeout=20, retries=1):
    last = None
    for k in range(retries + 1):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                etag = r.headers.get("ETag", "")
                body = r.read().decode("utf-8", "ignore")
                return json.loads(body), etag
        except Exception as e:
            last = e
            if k < retries:
                time.sleep(1.0)
    return None, str(last)


def probe_versions(online=True):
    """探测外部数据库版本。返回 {name: {release, source, fetched_at, status}}。"""
    out = {}
    if not online:
        return {"status": "OFFLINE"}
    j, err = _get_json("https://string-db.org/api/json/version")
    if isinstance(j, list) and j:
        out["STRING"] = {"release": j[0].get("string_version") or j[0].get("stringVersion"),
                         "source": "https://string-db.org/api/json/version",
                         "status": "OK"}
    elif isinstance(j, dict):
        out["STRING"] = {"release": j.get("string_version") or j.get("version"),
                         "source": "https://string-db.org/api/json/version", "status": "OK"}
    else:
        out["STRING"] = {"release": None, "status": "FAIL", "error": str(err)[:100]}

    j2, err2 = _get_json("https://rest.ensembl.org/info/data?content-type=application/json")
    if isinstance(j2, dict) and j2.get("releases"):
        out["Ensembl"] = {"release": max(j2["releases"]),
                          "source": "https://rest.ensembl.org/info/data", "status": "OK"}
    else:
        out["Ensembl"] = {"release": None, "status": "FAIL", "error": str(err2)[:100]}

    # MyGene（无版本号，记录 API 版本字段）
    j3, err3 = _get_json("https://mygene.info/v3/metadata")
    out["MyGene"] = {"release": (j3 or {}).get("build_version") if isinstance(j3, dict) else None,
                     "source": "https://mygene.info/v3/metadata",
                     "status": "OK" if isinstance(j3, dict) else "FAIL"}

    for k in out:
        out[k]["fetched_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    return out


def write_meta(path, url="", release=None, extra=None):
    """给缓存文件写同名 .meta.json（sha256 + 大小 + 时间 + 版本）。"""
    if not path or not os.path.exists(path):
        return None
    try:
        h = hashlib.sha256(open(path, "rb").read()).hexdigest()
    except Exception:
        return None
    meta = {"file": os.path.basename(path), "size": os.path.getsize(path),
            "sha256": h, "sha256_16": h[:16],
            "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "mtime": time.strftime("%Y-%m-%d %H:%M:%S",
                                   time.localtime(os.path.getmtime(path))),
            "url": url, "release": release or "unrecorded"}
    if extra:
        meta.update(extra)
    mp = path + ".meta.json"
    try:
        json.dump(meta, open(mp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        return mp
    except Exception:
        return None


def _walk_meta(root):
    rows = []
    for p in glob.glob(os.path.join(root, "**", "*"), recursive=True):
        if not os.path.isfile(p) or p.endswith(".meta.json"):
            continue
        mp = p + ".meta.json"
        if os.path.exists(mp):
            try:
                j = json.load(open(mp, encoding="utf-8"))
            except Exception:
                j = {}
        else:
            j = {}
        if not j.get("sha256"):
            try:
                h = hashlib.sha256(open(p, "rb").read()).hexdigest()
            except Exception:
                continue
            j = {"sha256": h, "sha256_16": h[:16], "size": os.path.getsize(p),
                 "release": "unrecorded", "url": ""}
        rows.append({"path": os.path.relpath(p, root).replace("\\", "/"),
                     "size": j.get("size"), "sha256": j.get("sha256"),
                     "release": j.get("release", "unrecorded"),
                     "url": j.get("url", ""),
                     "mtime": j.get("mtime", "")})
    return rows


def pack(out_path, root=None, versions=None, note=""):
    root = root or cache_root()
    if not os.path.isdir(root):
        return {"status": "NO_CACHE", "root": root}
    files = _walk_meta(root)
    if not files:
        return {"status": "EMPTY_CACHE", "root": root}
    man = {"tool": "network-pharmacology-pro", "snapshot_version": VERSION,
           "packed_at": time.strftime("%Y-%m-%d %H:%M:%S"),
           "cache_root": "<snapshot>", "n_files": len(files),
           "total_bytes": sum(f.get("size") or 0 for f in files),
           "files": files,
           "versions": versions if versions is not None else probe_versions(False),
           "note": note or ("离线快照：解压后设 NP_CACHE/NP_SNAPSHOT 指向该目录即可无网复现。")}
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    man_tmp = out_path + ".manifest.json"
    json.dump(man, open(man_tmp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    with tarfile.open(out_path, "w:gz") as tf:
        tf.add(root, arcname="cache")
        tf.add(man_tmp, arcname="snapshot_manifest.json")
    os.remove(man_tmp)
    h = hashlib.sha256(open(out_path, "rb").read()).hexdigest()
    man["tarball"] = os.path.basename(out_path)
    man["tarball_sha256"] = h
    man["tarball_sha256_16"] = h[:16]
    json.dump(man, open(out_path + ".meta.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    return {"status": "OK", "path": out_path, "n_files": len(files),
            "bytes": os.path.getsize(out_path), "sha256_16": h[:16]}


def mount(tar_path, to_dir):
    """解压快照并写出环境变量脚本（NP_CACHE / NP_SNAPSHOT）。"""
    if not os.path.exists(tar_path):
        return {"status": "NOT_FOUND"}
    os.makedirs(to_dir, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as tf:
        try:
            tf.extractall(to_dir, filter="data")     # 3.12+：拒绝越界/绝对路径成员
        except TypeError:
            tf.extractall(to_dir)
    cache = os.path.join(to_dir, "cache")
    man = os.path.join(to_dir, "snapshot_manifest.json")
    info = {}
    if os.path.exists(man):
        try:
            info = json.load(open(man, encoding="utf-8"))
        except Exception:
            info = {}
    sh = os.path.join(to_dir, "activate_snapshot.sh")
    open(sh, "w", encoding="utf-8").write(
        "# 由 np_snapshot.py 生成：source 本文件即可离线复现\n"
        "export NP_CACHE=\"%s\"\nexport NP_SNAPSHOT=\"%s\"\nexport NP_OFFLINE=0\n"
        % (cache.replace("\\", "/"), cache.replace("\\", "/")))
    bat = os.path.join(to_dir, "activate_snapshot.bat")
    open(bat, "w", encoding="utf-8").write(
        "@echo off\r\nREM 由 np_snapshot.py 生成：call 本文件即可离线复现\r\n"
        "set NP_CACHE=%s\r\nset NP_SNAPSHOT=%s\r\n" % (cache, cache))
    return {"status": "OK", "cache": cache, "n_files": info.get("n_files"),
            "versions": info.get("versions"), "activate_sh": sh, "activate_bat": bat}


def verify(tar_path):
    """校验快照内的文件哈希，并列出"缺失但可能需要"的缓存类型。"""
    if not os.path.exists(tar_path):
        return {"status": "NOT_FOUND"}
    with tarfile.open(tar_path, "r:gz") as tf:
        names = tf.getnames()
        man = None
        if "snapshot_manifest.json" in names:
            man = json.loads(tf.extractfile("snapshot_manifest.json").read().decode("utf-8"))
    if not man:
        return {"status": "NO_MANIFEST"}
    bad = []
    tmp = tar_path + ".verify"
    os.makedirs(tmp, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as tf:
        for f in man.get("files", []):
            try:
                data = tf.extractfile("cache/" + f["path"]).read()
            except Exception:
                bad.append({**f, "reason": "missing_in_tar"})
                continue
            if hashlib.sha256(data).hexdigest() != f["sha256"]:
                bad.append({**f, "reason": "sha256_mismatch"})
    kinds = {f["path"].split("/")[0] for f in man.get("files", [])}
    expected = {"string", "hpa", "ortholog"}
    return {"status": "OK" if not bad else "CORRUPT", "n_files": man.get("n_files"),
            "cache_kinds": sorted(kinds), "missing_kinds": sorted(expected - kinds),
            "bad": bad[:10], "versions": man.get("versions")}



def versions_used(root=None):
    """汇总"本次运行实际用到的数据库版本"：优先 versions.json，其次缓存 meta。"""
    root = root or cache_root()
    out = {"cache_root": root, "source": "unrecorded"}
    vp = os.path.join(root, "versions.json")
    if os.path.exists(vp):
        try:
            v = json.load(open(vp, encoding="utf-8"))
            out.update({k: (val.get("release") if isinstance(val, dict) else val)
                        for k, val in v.items() if isinstance(val, dict)})
            out["source"] = "versions.json"
            out["probed_at"] = max((val.get("fetched_at", "") for val in v.values()
                                    if isinstance(val, dict)), default="")
            return out
        except Exception:
            pass
    rel = {}
    for f in _walk_meta(root):
        kind = f["path"].split("/")[0]
        if f.get("release") and f["release"] != "unrecorded":
            rel.setdefault(kind, set()).add(f["release"])
    if rel:
        out.update({k: ",".join(sorted(v)) for k, v in rel.items()})
        out["source"] = "cache_meta"
    else:
        out["note"] = ("未记录任何数据库版本。请运行 `python np_snapshot.py probe --save`，"
                       "否则 Methods 只能写「latest」，投稿周期内不可复现。")
    return out


def selftest():
    checks = []
    import tempfile
    d = tempfile.mkdtemp(prefix="np_snap_")
    cache = os.path.join(d, "cache")
    os.makedirs(os.path.join(cache, "string"), exist_ok=True)
    p = os.path.join(cache, "string", "resp.json")
    open(p, "w", encoding="utf-8").write('{"ok":1}')
    mp = write_meta(p, url="https://example.org/string", release="12.0")
    checks.append(("meta_written", mp is not None and os.path.exists(mp), str(mp)))

    tar = os.path.join(d, "snap.tar.gz")
    r = pack(tar, root=cache, versions={"STRING": {"release": "12.0"}})
    checks.append(("pack_ok", r.get("status") == "OK" and os.path.getsize(tar) > 100,
                   str(r.get("status"))))

    m = mount(tar, os.path.join(d, "mnt"))
    checks.append(("mount_ok", m.get("status") == "OK"
                   and os.path.exists(os.path.join(m["cache"], "string", "resp.json")), ""))
    checks.append(("activate_files", os.path.exists(m.get("activate_sh", ""))
                   and os.path.exists(m.get("activate_bat", "")), ""))

    v = verify(tar)
    checks.append(("verify_ok", v.get("status") == "OK", str(v.get("status"))))
    checks.append(("verify_missing_kinds", "hpa" in (v.get("missing_kinds") or []),
                   str(v.get("missing_kinds"))))

    # 篡改检测
    bad_tar = os.path.join(d, "bad.tar.gz")
    import shutil
    shutil.copy(tar, bad_tar)
    with tarfile.open(bad_tar, "w:gz") as tf:
        tmpf = os.path.join(d, "resp_tampered.json")
        open(tmpf, "w", encoding="utf-8").write('{"ok":2}')
        tf.add(tmpf, arcname="cache/string/resp.json")
        with tarfile.open(tar, "r:gz") as src:
            tf.addfile(src.getmember("snapshot_manifest.json"),
                       src.extractfile("snapshot_manifest.json"))
    v2 = verify(bad_tar)
    checks.append(("tamper_detected", v2.get("status") == "CORRUPT", str(v2.get("status"))))

    ok = all(c[1] for c in checks)
    print("selftest: %s（%d/%d）" % ("PASS" if ok else "FAIL", sum(c[1] for c in checks), len(checks)))
    for n, g, d_ in checks:
        print("  [%s] %-24s %s" % ("✓" if g else "✗", n, d_))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["selftest"], ["--selftest"]):
        return selftest()
    ap = argparse.ArgumentParser(description="np_snapshot: 版本化快照与离线复现")
    sub = ap.add_subparsers(dest="fn")
    p = sub.add_parser("probe")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--write-meta", action="store_true", help="给现有缓存文件补 .meta.json")
    p.add_argument("--save", action="store_true", help="把版本写入 <cache>/versions.json 供 manifest 引用")
    p = sub.add_parser("pack")
    p.add_argument("-o", "--out", default=None)
    p.add_argument("--cache", default=None)
    p.add_argument("--note", default="")
    p = sub.add_parser("mount")
    p.add_argument("tar")
    p.add_argument("--to", required=True)
    p = sub.add_parser("verify")
    p.add_argument("tar")
    a = ap.parse_args(argv)

    if a.fn == "probe":
        v = probe_versions(not a.offline)
        print(json.dumps(v, ensure_ascii=False, indent=1))
        if a.write_meta:
            root = cache_root()
            n = 0
            for f in _walk_meta(root):
                pp = os.path.join(root, f["path"].replace("/", os.sep))
                if not os.path.exists(pp + ".meta.json"):
                    rel = (v.get("STRING") or {}).get("release") if "string" in f["path"] else None
                    if write_meta(pp, release=rel):
                        n += 1
            print("[snapshot] 补写 meta：%d 个文件" % n)
        if a.save:
            root = cache_root()
            os.makedirs(root, exist_ok=True)
            json.dump(v, open(os.path.join(root, "versions.json"), "w", encoding="utf-8"),
                      ensure_ascii=False, indent=1)
            print("[snapshot] 版本 → %s" % os.path.join(root, "versions.json"))
        return 0
    if a.fn == "pack":
        out = a.out or ("snapshot_%s.tar.gz" % time.strftime("%Y%m%d"))
        v = probe_versions(True)
        r = pack(out, root=a.cache, versions=v, note=a.note)
        print(json.dumps(r, ensure_ascii=False))
        return 0 if r.get("status") == "OK" else 2
    if a.fn == "mount":
        r = mount(a.tar, a.to)
        print(json.dumps(r, ensure_ascii=False, indent=1))
        return 0 if r.get("status") == "OK" else 2
    if a.fn == "verify":
        r = verify(a.tar)
        print(json.dumps(r, ensure_ascii=False, indent=1))
        return 0 if r.get("status") == "OK" else 2
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
