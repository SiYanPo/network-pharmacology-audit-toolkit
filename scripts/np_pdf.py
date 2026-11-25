#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_pdf.py —— 不依赖任何 Python 第三方库的 HTML→PDF（调用系统 Chrome/Edge 无头模式）

为什么放弃 weasyprint
    weasyprint 在 Windows 上要 GTK3/Pango/Cairo 一整套 C 库，装不上的概率远高于装得上；
    而"报告只能出 HTML"对合作者是不可接受的。系统里几乎一定有 Edge（Windows）或 Chrome，
    它们自带 `--headless --print-to-pdf`，效果比 weasyprint 更接近浏览器所见（CSS 变量、
    flex、@media print 全部支持），且**零 Python 依赖**。

跨平台发现顺序
    Windows : Edge(Program Files / Program Files (x86) / LocalAppData) → Chrome(同上) → Brave
    macOS   : Google Chrome.app → Microsoft Edge.app → Chromium.app
    Linux   : google-chrome / google-chrome-stable / chromium / chromium-browser / microsoft-edge
    另可用环境变量 NP_BROWSER 或 --browser 显式指定。

用法
    python np_pdf.py --html report.html [--pdf report.pdf] [--browser "...msedge.exe"]
    python np_pdf.py --check          # 只探测可用浏览器
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time


def _win_candidates():
    roots = [os.environ.get("ProgramFiles", r"C:\Program Files"),
             os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)"),
             os.environ.get("LOCALAPPDATA", "")]
    rel = [
        r"Microsoft\Edge\Application\msedge.exe",
        r"Google\Chrome\Application\chrome.exe",
        r"BraveSoftware\Brave-Browser\Application\brave.exe",
        r"Chromium\Application\chrome.exe",
    ]
    return [os.path.join(r, x) for r in roots if r for x in rel]


def _mac_candidates():
    return ["/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
            "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
            "/Applications/Chromium.app/Contents/MacOS/Chromium",
            "/Applications/Brave Browser.app/Contents/MacOS/Brave Browser"]


def _linux_candidates():
    names = ["google-chrome", "google-chrome-stable", "chromium", "chromium-browser",
             "microsoft-edge", "microsoft-edge-stable", "brave-browser"]
    return [shutil.which(n) for n in names if shutil.which(n)]


def find_browser(explicit=None):
    """返回 (path, source) 或 (None, reason)。"""
    if explicit:
        return (explicit, "explicit") if os.path.exists(explicit) else (None, "explicit_not_found:%s" % explicit)
    env = os.environ.get("NP_BROWSER")
    if env and os.path.exists(env):
        return env, "env:NP_BROWSER"
    if sys.platform.startswith("win"):
        cands = _win_candidates()
    elif sys.platform == "darwin":
        cands = _mac_candidates()
    else:
        return (cands[0], "PATH") if (cands := _linux_candidates()) else (None, "no_chromium_on_path")
    for p in cands:
        if p and os.path.exists(p):
            return p, "auto"
    return None, "not_found（可安装 Edge/Chrome，或设置环境变量 NP_BROWSER）"


def html_to_pdf(html_path, pdf_path=None, browser=None, timeout=120):
    """用无头浏览器打印 PDF。返回 (pdf_path|None, info_dict)。"""
    html_path = os.path.abspath(html_path)
    if not os.path.exists(html_path):
        return None, {"status": "HTML_MISSING", "path": html_path}
    pdf_path = os.path.abspath(pdf_path or os.path.splitext(html_path)[0] + ".pdf")
    exe, src = find_browser(browser)
    if not exe:
        return None, {"status": "NO_BROWSER", "reason": src,
                      "hint": "安装 Edge/Chrome，或 --browser 指定路径 / 设 NP_BROWSER"}
    os.makedirs(os.path.dirname(pdf_path) or ".", exist_ok=True)
    if os.path.exists(pdf_path):
        try:
            os.remove(pdf_path)
        except OSError:
            pass
    uri = "file:///" + html_path.replace("\\", "/").lstrip("/") if sys.platform.startswith("win") \
        else "file://" + html_path
    base = [exe, "--headless=new", "--disable-gpu", "--no-first-run", "--no-default-browser-check",
            "--disable-extensions", "--disable-background-networking", "--virtual-time-budget=15000",
            "--no-pdf-header-footer", "--print-to-pdf-no-header", "--print-to-pdf=%s" % pdf_path]
    attempts = []
    for flags in (base[:], [f for f in base if f != "--headless=new"] + ["--headless"]):
        cmd = flags + [uri]
        t0 = time.time()
        try:
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               timeout=timeout, text=True, encoding="utf-8", errors="replace")
            rc, out = p.returncode, (p.stdout or "")
        except subprocess.TimeoutExpired:
            rc, out = -1, "timeout"
        except Exception as e:
            rc, out = -2, "%s:%s" % (type(e).__name__, str(e)[:120])
        attempts.append({"rc": rc, "flags": os.path.basename(flags[1]) if len(flags) > 1 else "",
                         "sec": round(time.time() - t0, 1), "out": out[-300:]})
        if os.path.exists(pdf_path) and os.path.getsize(pdf_path) > 5120:
            with open(pdf_path, "rb") as f:
                magic = f.read(4)
            if magic == b"%PDF":
                return pdf_path, {"status": "OK", "browser": exe, "browser_source": src,
                                  "size_kb": round(os.path.getsize(pdf_path) / 1024, 1),
                                  "attempts": attempts}
    return None, {"status": "FAILED", "browser": exe, "attempts": attempts,
                  "hint": "浏览器返回非 PDF。常见原因：文件路径含特殊字符、浏览器被策略禁用无头模式"}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(description="np_pdf: 用系统浏览器无头模式把 HTML 打成 PDF")
    ap.add_argument("--html", help="输入 HTML")
    ap.add_argument("--pdf", default=None)
    ap.add_argument("--browser", default=None)
    ap.add_argument("--check", action="store_true", help="只探测浏览器")
    a = ap.parse_args(argv)
    if a.check:
        exe, src = find_browser(a.browser)
        if exe:
            print("[pdf] 可用浏览器：%s（%s）" % (exe, src))
            try:
                v = subprocess.run([exe, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   timeout=15, text=True)
                print("[pdf] 版本：%s" % (v.stdout or "").strip())
            except Exception:
                pass
            return 0
        print("[pdf] 未找到 Chrome/Edge：%s" % src)
        return 2
    if not a.html:
        ap.error("需 --html（或用 --check）")
    p, info = html_to_pdf(a.html, a.pdf, a.browser)
    if p:
        print("[pdf] → %s（%.0f KB，%s）" % (p, info["size_kb"], os.path.basename(info["browser"])))
        return 0
    print("[pdf] 失败：%s" % info.get("status"))
    for k in ("reason", "hint"):
        if info.get(k):
            print("      %s" % info[k])
    for at in info.get("attempts", []):
        print("      rc=%s %.1fs %s" % (at["rc"], at["sec"], (at["out"] or "").splitlines()[-1][:120]
                                        if at["out"] else ""))
    return 2


if __name__ == "__main__":
    sys.exit(main())
