#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
stp_submit.py — SwissTargetPrediction 提交器（非官方；遵守其服务器策略）

为什么需要它：
  STP 一次只接受一个分子、无 API、且对 IP 限流（"submitting too many calculations"）。
  手工粘贴容易出错（用户实测把非法字符串粘进去：弹窗原样回显 "YYYYjbYYYY"）。
  本脚本按 ToolUniverse 开源实现确认的协议提交，带限流退避，逐化合物落盘 CSV。

协议（2026-11 实测确认）：
  1) POST /predict.php  body: organism=Homo_sapiens&smiles=<urlenc>&Example=&ioi=2
     （ioi=2 必须保留，否则拒绝）
  2) 从响应 url/text 提取 job=ID
  3) GET /result.php?job=X&organism=Y，轮询直到出现 id="resultTable"
  4) 解析表格 → <compound>.csv

用法
    python tools/stp_submit.py --input out/stp/stp_input.txt --names a,b,c -o out/stp/raw
    python tools/stp_submit.py --smiles "CC(=O)Oc1ccccc1C(=O)O" --name Aspirin -o raw/
"""
import argparse
import csv
import os
import re
import sys
import time
import http.cookiejar
import urllib.parse
import urllib.request
import ssl

BASE = "https://www.swisstargetprediction.ch"
ORG = "Homo_sapiens"
UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126 Safari/537.36")
COLS = ["Target", "Common name", "Uniprot ID", "ChEMBL ID", "Target Class",
        "Probability", "Known actives (3D)", "Known actives (2D)"]
RATE_MARKS = ("too many calculations", "will not be submitted")


def _ctx():
    c = ssl.create_default_context()
    c.check_hostname = False
    c.verify_mode = ssl.CERT_NONE
    return c


def _opener():
    cj = http.cookiejar.CookieJar()
    return urllib.request.build_opener(
        urllib.request.HTTPCookieProcessor(cj),
        urllib.request.HTTPSHandler(context=_ctx()))


def _hdr(extra=None):
    h = {"User-Agent": UA, "Accept-Language": "en-US,en;q=0.9",
         "Upgrade-Insecure-Requests": "1"}
    if extra:
        h.update(extra)
    return h


def submit(op, smiles, timeout=90):
    """提交一个分子。返回 (job_id, err)。限流/拒绝时 err 以 RATE_LIMIT/REJECTED 开头。"""
    body = urllib.parse.urlencode({"organism": ORG, "smiles": smiles,
                                   "Example": "", "ioi": "2"}).encode()
    req = urllib.request.Request(BASE + "/predict.php", data=body, headers=_hdr(
        {"Content-Type": "application/x-www-form-urlencoded",
         "Origin": BASE, "Referer": BASE + "/"}), method="POST")
    try:
        r = op.open(req, timeout=timeout)
        txt = r.read().decode("utf-8", "ignore")
        final = r.geturl()
    except Exception as e:
        return None, "NET:%s %s" % (type(e).__name__, str(e)[:80])
    if "too many calculations" in txt:
        return None, "RATE_LIMIT"
    if "will not be submitted" in txt:
        m = re.search(r'Your\s+(.+?)\s+will not be submitted', txt)
        return None, "REJECTED:%s" % (m.group(1)[:60] if m else "?")
    m = re.search(r"job=(\d+)", final) or re.search(r"job=(\d+)", txt)
    if not m:
        return None, "NO_JOB"
    return m.group(1), ""


def fetch_table(op, job, timeout=40, poll=15, max_poll=150):
    """轮询结果页直到出现 resultTable。返回 (rows, err)。"""
    url = "%s/result.php?job=%s&organism=%s" % (BASE, job, ORG)
    t0 = time.time()
    while time.time() - t0 < max_poll:
        try:
            r = op.open(urllib.request.Request(url, headers=_hdr(
                {"Referer": BASE + "/predict.php"})), timeout=timeout)
            txt = r.read().decode("utf-8", "ignore")
        except Exception as e:
            time.sleep(poll)
            continue
        if "too many calculations" in txt:
            return None, "RATE_LIMIT"
        m = re.search(r'<table[^>]*id="resultTable"[^>]*>(.*?)</table>', txt, re.S | re.I)
        if m:
            return _parse_table(m.group(1)), ""
        time.sleep(poll)
    return None, "TIMEOUT"


def _clean(s):
    return re.sub(r"<[^>]+>", "", s).replace("&amp;", "&").strip()


def _parse_table(frag):
    rows = []
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", frag, re.S | re.I):
        tds = re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S | re.I)
        if len(tds) < 6:
            continue
        prob = _clean(tds[5]).replace("*", "")
        k3 = _clean(tds[6]) if len(tds) > 6 else ""
        k2 = _clean(tds[7]) if len(tds) > 7 else ""
        rows.append({"Target": _clean(tds[0]), "Common name": _clean(tds[1]),
                     "Uniprot ID": _clean(tds[2]), "ChEMBL ID": _clean(tds[3]),
                     "Target Class": _clean(tds[4]),
                     "Probability": prob.split()[0] if prob else "",
                     "Known actives (3D)": k3.split("/")[0].strip(),
                     "Known actives (2D)": k2.split("/")[0].strip()})
    return rows


def save_csv(rows, path):
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def run(names, smiles_list, outdir, cooldown=90, rate_wait=180, max_rate_waits=15):
    """逐个提交。限流时等 rate_wait 秒再试（最多 max_rate_waits 次）。"""
    os.makedirs(outdir, exist_ok=True)
    op = _opener()
    done, failed = [], []
    for i, (name, smi) in enumerate(zip(names, smiles_list)):
        out_csv = os.path.join(outdir, "%s.csv" % name)
        if os.path.exists(out_csv):
            print("  = %-18s 已存在，跳过" % name)
            done.append(name)
            continue
        waits = 0
        while True:
            job, err = submit(op, smi)
            if job:
                break
            if err.startswith("RATE_LIMIT") and waits < max_rate_waits:
                waits += 1
                print("  ! %-18s 限流，等 %ds 再试（%d/%d）" % (name, rate_wait, waits, max_rate_waits))
                time.sleep(rate_wait)
                continue
            print("  x %-18s 提交失败：%s" % (name, err))
            failed.append((name, err))
            break
        if not job:
            continue
        rows, err2 = fetch_table(op, job)
        if rows is None:
            print("  x %-18s 取结果失败：%s（job=%s）" % (name, err2, job))
            failed.append((name, err2))
            continue
        save_csv(rows, out_csv)
        done.append(name)
        print("  ✓ %-18s job=%s  %d 个靶点 → %s" % (name, job, len(rows), out_csv))
        if i < len(names) - 1:
            time.sleep(cooldown)          # 遵守服务器策略：分子间冷却
    print("\n完成 %d/%d" % (len(done), len(names)))
    if failed:
        print("失败：", "; ".join("%s(%s)" % f for f in failed))
    return 0 if len(done) == len(names) else 1


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", help="SMILES 文件（一行一个，顺序与 --names 对应）")
    ap.add_argument("--smiles", help="单个 SMILES（与 --name 搭配）")
    ap.add_argument("--name", default="compound")
    ap.add_argument("--names", help="逗号分隔的化合物名（与 --input 对应）")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--cooldown", type=int, default=90, help="分子间隔秒数")
    ap.add_argument("--rate-wait", type=int, default=180, help="限流后退避秒数")
    ap.add_argument("--max-rate-waits", type=int, default=15)
    a = ap.parse_args(argv)
    if a.smiles:
        names, smis = [a.name], [a.smiles]
    elif a.input:
        smis = [x.strip() for x in open(a.input, encoding="utf-8") if x.strip()]
        names = ([x.strip() for x in a.names.split(",")] if a.names
                 else ["c%d" % i for i in range(len(smis))])
    else:
        print("需要 --input 或 --smiles")
        return 2
    assert len(names) == len(smis), "名字数与 SMILES 数不一致"
    return run(names, smis, a.outdir, a.cooldown, a.rate_wait, a.max_rate_waits)


if __name__ == "__main__":
    sys.exit(main())
