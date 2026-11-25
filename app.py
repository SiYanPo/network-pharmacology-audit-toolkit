#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
app.py —— network-pharmacology-pro 一站式 Web 平台（Streamlit）

    streamlit run app.py
    # 或指定端口：streamlit run app.py --server.port 8501

四个页签
    1 一键分析   ：输入方剂/药味/疾病 → 带进度条调度 np_pipeline.py → 生成带水印报告并下载
    2 基准测试   ：内置 4 个经典方剂金标准，上传排序表/边表 → AUC / Recall@K / EF + benchmark_report.md
    3 因果校验   ：上传 MR 结果 + 成分动作表 → Target_Intervention_Mismatch 判定
    4 环境自检   ：探测 vina/obabel/rdkit/gmx 等真实可用性（不猜 PATH）

设计约束（为什么这么写）
    * 进度条来自**真实节点执行**：解析 np_pipeline 打印的 [node_id] 推进，不做假动画；
    * 报告数字全部读磁盘产物（np_uireport 只读不重算），不显示任何未落盘的数字；
    * 水印由 np_finalcheck 结论决定；没跑门禁 = UNVERIFIED，绝不默认 Ready；
    * 任何外部依赖缺失都走"诚实降级"（PDF 无 weasyprint 就只给 HTML 并说明）。
"""
from __future__ import annotations

import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(ROOT, "scripts")
sys.path.insert(0, SCRIPTS)

try:
    import streamlit as st
except ImportError:  # 允许在没有 streamlit 的环境里 import 本文件做单测
    st = None

import np_pipeline  # noqa: E402  复用配置加载/拓扑排序，避免 UI 与引擎两套逻辑
import np_uireport  # noqa: E402

CONFIG = os.path.join(ROOT, "config", "pipeline.yaml")
EXAMPLES = os.path.join(ROOT, "examples")
WORKDIR = os.environ.get("NP_WORKDIR", os.path.join(ROOT, "np_ui_runs"))
PLUGIN_FLAGS = [("enable-funnel", "先筛选后对接漏斗"), ("enable-doseweight", "剂量加权拓扑"),
                ("enable-louvain", "Louvain 子网络")]


# ============================ 基础工具 ============================
def plan_nodes(enable_flags):
    """与引擎同源的执行计划（含拓扑序）。UI 的进度条分母来自这里。"""
    cfg, _ = np_pipeline.load_config(CONFIG)
    enable = {"--" + f for f in enable_flags}
    sel = np_pipeline.select_nodes(cfg, enable, set(), set())
    by_id = {s["id"]: s for s in (cfg.get("stages") or [])}
    order = np_pipeline.topo_sort([by_id[n["id"]] for n in sel])
    return [n["id"] for n in order if n["id"] in {x["id"] for x in sel}]


def build_inputs(outdir, herbs, disease_label, gmt=None, files=None, tissue=""):
    """组装 inputs.json。files 为空时用 examples/ 的构造数据（开箱即用）。"""
    os.makedirs(outdir, exist_ok=True)
    herbs_file = os.path.join(outdir, "herbs.txt")
    with open(herbs_file, "w", encoding="utf-8") as f:
        f.write("\n".join([h.strip() for h in re.split(r"[,，;；\s]+", herbs) if h.strip()]) + "\n")
    if files:
        paths = files
    else:
        paths = {"compounds": os.path.join(EXAMPLES, "compounds.csv"),
                 "compound_edges": os.path.join(EXAMPLES, "compound_targets.csv"),
                 "disease": os.path.join(EXAMPLES, "disease_targets.tsv"),
                 "edges": os.path.join(EXAMPLES, "string_edges.tsv")}
        gmt = gmt or os.path.join(EXAMPLES, "pathways.gmt")
    cfg = {"inputs": {"herbs": ",".join([h.strip() for h in re.split(r"[,，;；\s]+", herbs) if h.strip()]),
                      "herbs_file": herbs_file,
                      "compounds": os.path.abspath(paths["compounds"]),
                      "compound_edges": os.path.abspath(paths["compound_edges"]),
                      "targets": os.path.abspath(paths["compound_edges"]),
                      "disease": os.path.abspath(paths["disease"]),
                      "edges": os.path.abspath(paths["edges"]),
                      "dose": "", "gmt": os.path.abspath(gmt) if gmt else ""},
           "params": {"label": "%s — %s" % (herbs or "复方", disease_label or "疾病"),
                      "organism": "mouse", "tissue": tissue or ""}}
    p = os.path.join(outdir, "inputs.json")
    json.dump(cfg, open(p, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return p


def run_stream(cmd, on_line=None, cwd=None, env=None):
    """逐行流式执行子进程；返回 (rc, 输出全文)。"""
    e = os.environ.copy()
    e["PYTHONUNBUFFERED"] = "1"
    e["PYTHONIOENCODING"] = "utf-8"
    if env:
        e.update(env)
    p = subprocess.Popen(cmd, cwd=cwd or ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, encoding="utf-8", errors="replace", bufsize=1, env=e)
    lines = []
    for line in p.stdout:
        line = line.rstrip("\n")
        lines.append(line)
        if on_line:
            on_line(line)
    p.wait()
    return p.returncode, "\n".join(lines)


# ============================ 页签 1：一键分析 ============================
def page_analyze():
    st.subheader("一键分析：方剂 → 靶点 → 通路 → 报告")
    c1, c2 = st.columns([1, 1])
    with c1:
        formula = st.text_input("方剂名称", value="复方丹参滴丸")
        herbs = st.text_input("药味（逗号分隔）", value="丹参,三七,冰片")
        disease = st.text_input("疾病 / 证候", value="冠心病稳定型心绞痛")
        tissue = st.text_input(
            "病变组织（HPA 富集背景）", value="",
            placeholder="heart / brain / liver / kidney …（留空则用全基因组，报告会强制打膨胀警告）",
            help="背景基因集必须与该疾病组织匹配。填错组织比留空更危险："
                 "留空至少会触发 P-value may be inflated 警告，填错则连警告都没有。")
    with c2:
        mode = st.radio("数据来源", ["内置演示数据（开箱即用）", "上传自己的数据", "只重出报告（已有输出目录）"],
                        help="演示数据为构造数据，仅用于跑通流程；正式分析请上传自己的表")
        run_label = re.sub(r"[\\/:*?\"<>|\s]+", "_", "%s_%s" % (formula, time.strftime("%Y%m%d_%H%M")))[:60]
        outdir = os.path.join(WORKDIR, run_label)
        plugins = st.multiselect("可选插件", [p[0] for p in PLUGIN_FLAGS], default=[])
        do_gate = st.checkbox("跑投稿门禁（np_finalcheck，决定水印为 READY 还是 DRAFT）", value=False)

    uploads = {}
    if mode == "上传自己的数据":
        st.caption("需要 4 个文件：成分表（compounds）、成分-靶点边表、疾病靶点表、PPI 边表；"
                   "列名要求见 examples/README.md。")
        u1 = st.file_uploader("成分表 CSV（MOL_ID,compound,OB,DL,herb）", type=["csv", "tsv"])
        u2 = st.file_uploader("成分-靶点边表（compound,target_name/gene）", type=["csv", "tsv"])
        u3 = st.file_uploader("疾病靶点表（gene,score,source）", type=["csv", "tsv"])
        u4 = st.file_uploader("PPI 边表（node1,node2,combined_score）", type=["csv", "tsv"])
        u5 = st.file_uploader("通路 GMT（可选，用于离线富集）", type=["gmt"])
        if all([u1, u2, u3, u4]):
            os.makedirs(os.path.join(outdir, "_in"), exist_ok=True)
            saved = {}
            for key, u in (("compounds", u1), ("compound_edges", u2), ("disease", u3), ("edges", u4)):
                p = os.path.join(outdir, "_in", key + os.path.splitext(u.name)[1])
                open(p, "wb").write(u.getbuffer())
                saved[key] = p
            gmt = None
            if u5:
                gmt = os.path.join(outdir, "_in", "pathways.gmt")
                open(gmt, "wb").write(u5.getbuffer())
            uploads = {"files": saved, "gmt": gmt}

    if mode == "只重出报告（已有输出目录）":
        existing = st.text_input("已有输出目录（含 pipeline_manifest.json）", value=outdir)
        if st.button("重新生成报告", type="primary"):
            _make_report(existing, "%s — %s" % (formula, disease), do_gate)
        return

    go = st.button("一键运行", type="primary", use_container_width=True)
    if not go:
        st.info("点击后按真实节点顺序执行：安全前置 → 别名 → 标准化 → 成分 → 交集 → PPI → 零模型 → 富集 → 同源 → 报告。"
                "进度条由节点事件推进，不是假动画。")
        return
    if mode == "上传自己的数据" and not uploads:
        st.error("请先上传 4 个必需文件（成分表、成分-靶点边、疾病靶点、PPI 边）。")
        return

    nodes = plan_nodes(plugins)
    pbar = st.progress(0.0)
    status = st.empty()
    logbox = st.empty()
    logs = []
    idx = {"i": 0}
    t0 = time.time()

    inp = build_inputs(outdir, herbs, disease,
                       gmt=uploads.get("gmt"), files=uploads.get("files"), tissue=tissue)

    def on_line(line):
        logs.append(line)
        m = re.match(r"^\[([a-z0-9_]+)\]\s", line.strip())
        if m and m.group(1) in nodes:
            i = nodes.index(m.group(1))
            idx["i"] = max(idx["i"], i)
            pbar.progress(min(1.0, (i) / max(1, len(nodes))))
            status.markdown("**正在运行：`%s`**（%d/%d）" % (m.group(1), i + 1, len(nodes)))
        if len(logs) % 6 == 0:
            logbox.code("\n".join(logs[-14:]), language="text")

    cmd = [sys.executable, "-u", os.path.join(SCRIPTS, "np_pipeline.py"),
           "--config", CONFIG, "--inputs", inp, "-o", outdir]
    for f in plugins:
        cmd.append("--" + f)
    status.markdown("**启动流水线…**")
    rc, out = run_stream(cmd, on_line)
    pbar.progress(1.0 if rc == 0 else min(0.98, (idx["i"] + 1) / max(1, len(nodes))))
    logbox.code("\n".join(logs[-40:]), language="text")
    status.markdown("**%s**（%.0fs）" % ("全部节点通过" if rc == 0 else "有节点失败，见日志", time.time() - t0))

    if rc != 0:
        st.error("流水线未全部通过。可修复后点击『一键运行』重跑，或勾选 --resume 思路："
                 "命令行 `python scripts/np_pipeline.py --inputs %s -o %s --resume`" % (inp, outdir))
    _make_report(outdir, "%s — %s" % (formula, disease), do_gate, show_spinner=False)


def _make_report(outdir, label, do_gate, show_spinner=True):
    """产出 HTML（+ 可选 PDF）并提供下载。"""
    if not os.path.isdir(outdir):
        st.error("目录不存在：%s" % outdir)
        return
    verdict = "UNVERIFIED"
    gate_note = "未运行投稿门禁"
    if do_gate:
        with st.spinner("运行投稿门禁（20 项自查）…"):
            rc, out = run_stream([sys.executable, os.path.join(SCRIPTS, "np_finalcheck.py"),
                                  "--dir", outdir, "--json", os.path.join(outdir, "finalcheck.json")])
        vp = os.path.join(outdir, "finalcheck.json")
        if os.path.exists(vp):
            verdict = json.load(open(vp, encoding="utf-8")).get("verdict", "DRAFT")
            gate_note = "门禁结论：%s" % verdict
        else:
            gate_note = "门禁未产出结果（%s）" % out.strip().splitlines()[-1][:120] if out.strip() else "门禁失败"

    with st.spinner("生成带水印报告…"):
        html = np_uireport.build(outdir, label, verdict)
        rp = os.path.join(outdir, "ui_report.html")
        open(rp, "w", encoding="utf-8").write(html)
        pdf = np_uireport.to_pdf(rp)

    st.success("报告已生成：%s（%s）" % (os.path.basename(rp), gate_note))
    c = np_uireport.counts(np_uireport.collect(outdir))
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("通过成分", c.get("compounds", "—"))
    m2.metric("交集靶点", c.get("targets", "—"))
    m3.metric("稳健 hub", c.get("hubs", "—"))
    m4.metric("显著通路", c.get("pathways", "—"))

    # 门禁摘要（让用户一眼知道哪几门没过）
    fj = os.path.join(outdir, "finalcheck.json")
    if os.path.exists(fj):
        try:
            jf = json.load(open(fj, encoding="utf-8"))
            fails = [r for r in (jf.get("results") or []) if not r.get("pass")]
            hard = [r["id"] for r in fails if r.get("level") == "hard"]
            st.markdown("**门禁：%s** ｜ 未通过 %d 门%s" % (
                jf.get("verdict"), len(fails), "（硬性：%s）" % "、".join(hard) if hard else ""))
            if fails:
                with st.expander("查看未通过项与修复命令"):
                    for r in fails[:12]:
                        st.markdown("- **%s %s**：%s" % (r.get("id"), r.get("name"), r.get("reason", "")))
                        try:
                            import np_finalcheck as FC
                            st.code(FC.fix_command(r.get("id"), outdir)[1], language="bash")
                        except Exception:
                            pass
        except Exception:
            pass

    d1, d2, d3, d4 = st.columns(4)
    d1.download_button("下载 HTML 报告", open(rp, "rb").read(), file_name="np_report.html",
                       mime="text/html", use_container_width=True)
    if pdf:
        d2.download_button("下载 PDF 报告", open(pdf, "rb").read(), file_name="np_report.pdf",
                           mime="application/pdf", use_container_width=True)
    else:
        d2.button("PDF 不可用（未找到 Chrome/Edge）", disabled=True, use_container_width=True)
    z = _zip_key_outputs(outdir)
    if z:
        d3.download_button("下载结果包 ZIP", open(z, "rb").read(), file_name="np_results.zip",
                           mime="application/zip", use_container_width=True)
    fc = os.path.join(outdir, "finalcheck_report.md")
    if os.path.exists(fc):
        d4.download_button("下载门禁清单", open(fc, "rb").read(), file_name="finalcheck_report.md",
                           use_container_width=True)
    st.caption("PDF 由系统 Edge/Chrome 的无头模式直接打印（零 Python 依赖，CSS/分页与浏览器所见一致）；"
               "若机器上没有 Chromium 内核浏览器，用 HTML + Ctrl+P 亦可。")
    with st.expander("报告预览（内嵌）"):
        import streamlit.components.v1 as components
        components.html(html, height=620, scrolling=True)


def _zip_key_outputs(outdir):
    import zipfile
    keep = []
    for pat in ("pipeline_manifest.json", "finalcheck*.json", "finalcheck_report.md",
                "*/*.csv", "*/*.md", "*/*.svg", "ui_report.html"):
        keep.extend(glob.glob(os.path.join(outdir, pat)))
    keep = [p for p in keep if os.path.isfile(p)]
    if not keep:
        return None
    z = os.path.join(outdir, "np_results.zip")
    with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in keep:
            zf.write(p, os.path.relpath(p, outdir))
    return z


# ============================ 页签 2：基准测试 ============================
def page_benchmark():
    st.subheader("基准测试：AUC / Recall@K / 富集倍数")
    gold = json.load(open(os.path.join(ROOT, "data", "benchmark_goldstandard.json"), encoding="utf-8"))
    opts = ["%s（%s，core %d）" % (f["name_zh"], f["id"], len(f["core"])) for f in gold["formulas"]]
    pick = st.selectbox("内置金标准方剂", opts)
    f = gold["formulas"][opts.index(pick)]
    st.caption("适应症：%s｜药味：%s｜金标准来源：DrugBank/TTD/ChEMBL + 标志性成分实验文献（见 provenance）" % (
        f["indication_zh"], "、".join(f["herbs"])))

    c1, c2 = st.columns([1, 1])
    with c1:
        tier = st.radio("口径", ["strict（仅 core）", "lenient（core+extended）"], horizontal=True)
        tier_key = "strict" if tier.startswith("strict") else "lenient"
        src = st.radio("排序来源", ["上传靶点排序表", "上传成分-靶点边表"], horizontal=True)
    with c2:
        topk = st.text_input("K 值（逗号分隔）", value="10,20,50")
        boot = st.number_input("Bootstrap 次数", 0, 5000, 500, step=100)

    pred = edges = None
    if src == "上传靶点排序表":
        u = st.file_uploader("靶点排序表（gene,score；单列基因列表亦可）", type=["csv", "tsv", "txt"])
        if u:
            os.makedirs(WORKDIR, exist_ok=True)
            pred = os.path.join(WORKDIR, "ranking_" + u.name)
            open(pred, "wb").write(u.getbuffer())
    else:
        u = st.file_uploader("成分-靶点边表（compound,gene）", type=["csv", "tsv"])
        if u:
            os.makedirs(WORKDIR, exist_ok=True)
            edges = os.path.join(WORKDIR, "edges_" + u.name)
            open(edges, "wb").write(u.getbuffer())
        wm = st.radio("排序基线", ["penalized（多靶点成分降权，默认）", "support（支持成分数）"], horizontal=True)
        st.session_state["_weight_mode"] = "penalized" if wm.startswith("penalized") else "support"

    if st.button("运行基准评测", type="primary", disabled=not (pred or edges)):
        outb = os.path.join(WORKDIR, "benchmark")
        cmd = [sys.executable, os.path.join(SCRIPTS, "np_benchmark.py"), "evaluate",
               "--formula", f["id"], "--tier", tier_key, "--topk", topk,
               "--bootstrap", str(boot), "-o", outb]
        cmd += ["--pred", pred] if pred else ["--edges", edges,
                                              "--weight-mode", st.session_state.get("_weight_mode", "penalized")]
        rc, outtxt = run_stream(cmd)
        st.code(outtxt[-1500:], language="text")
        jp = os.path.join(outb, "benchmark_%s.json" % f["id"])
        if os.path.exists(jp):
            r = json.load(open(jp, encoding="utf-8"))
            if r.get("status") == "OK":
                a, b, c, d = st.columns(4)
                a.metric("AUC", "%.3f" % r["AUC"], help="95% CI %s" % (r.get("AUC_95CI") or "NA"))
                b.metric("AUPRC", "%.3f" % (r["AUPRC"] or 0))
                c.metric("Recall@20", "%.2f" % r["at_k"]["20"]["recall"],
                         help="随机期望 %.2f" % r["at_k"]["20"]["expected_hits"])
                d.metric("EF@20", "%.2fx" % r["at_k"]["20"]["enrichment_factor"])
                if r.get("coverage_warning"):
                    st.warning(r["coverage_warning"])
                st.dataframe([{"K": k, **v} for k, v in r["at_k"].items()], use_container_width=True)
                svg = os.path.join(outb, "roc_%s.svg" % f["id"])
                if os.path.exists(svg):
                    st.image(svg, caption="ROC（AUC=%.3f）" % r["AUC"])
                run_stream([sys.executable, os.path.join(SCRIPTS, "np_benchmark.py"),
                            "report", "--dir", outb])
                rep = os.path.join(outb, "benchmark_report.md")
                if os.path.exists(rep):
                    st.download_button("下载 benchmark_report.md", open(rep, "rb").read(),
                                       file_name="benchmark_report.md", use_container_width=True)
            else:
                st.error(r.get("reason", "无交集"))

    st.divider()
    st.markdown("**跨库边级一致性（推荐用于比较两个数据库）**")
    e1 = st.file_uploader("预测边表", type=["csv", "tsv"], key="ce1")
    e2 = st.file_uploader("独立金标准边表", type=["csv", "tsv"], key="ce2")
    if e1 and e2 and st.button("计算敏感性/精确率"):
        os.makedirs(WORKDIR, exist_ok=True)
        p1 = os.path.join(WORKDIR, "pred_" + e1.name)
        p2 = os.path.join(WORKDIR, "gold_" + e2.name)
        open(p1, "wb").write(e1.getbuffer())
        open(p2, "wb").write(e2.getbuffer())
        outb = os.path.join(WORKDIR, "benchmark")
        rc, outtxt = run_stream([sys.executable, os.path.join(SCRIPTS, "np_benchmark.py"),
                                 "concordance", "--pred-edges", p1, "--gold-edges", p2,
                                 "--label", "%s vs %s" % (e1.name, e2.name), "-o", outb])
        st.code(outtxt, language="text")


# ============================ 页签 3：因果校验 ============================
def page_causal():
    st.subheader("靶点-干预方向一致性（Target_Intervention_Mismatch）")
    st.caption("规则：MR β>0（靶点升高→风险升高）→ 起效干预必须是【抑制】；MR β<0 → 必须是【激活】。"
               "方向相反即 mismatch：SEVERE 仅在 Grade A 证据下抛出，Grade B 降为 WARN。")
    c1, c2 = st.columns(2)
    with c1:
        mrs = st.file_uploader("MR 结果 JSON（可多选）", type=["json"], accept_multiple_files=True)
    with c2:
        act = st.file_uploader("成分动作表（gene,action 或 gene,compound,action）", type=["csv", "tsv"])
    if st.button("运行校验", type="primary", disabled=not (mrs and act)):
        outd = os.path.join(WORKDIR, "intervention")
        mrd = os.path.join(outd, "mr")
        os.makedirs(mrd, exist_ok=True)
        for u in mrs:
            open(os.path.join(mrd, u.name), "wb").write(u.getbuffer())
        ap = os.path.join(outd, "action" + os.path.splitext(act.name)[1])
        open(ap, "wb").write(act.getbuffer())
        rc, outtxt = run_stream([sys.executable, os.path.join(SCRIPTS, "np_intervention.py"),
                                 "check", "--mr-dir", mrd, "--action", ap, "-o", outd])
        st.code(outtxt, language="text")
        csvp = os.path.join(outd, "intervention_check.csv")
        if os.path.exists(csvp):
            import csv as _csv
            rows = list(_csv.DictReader(open(csvp, encoding="utf-8-sig")))
            st.dataframe(rows, use_container_width=True)
            st.download_button("下载 intervention_check.csv", open(csvp, "rb").read(),
                               file_name="intervention_check.csv")
        if rc == 3:
            st.error("检出 SEVERE 方向冲突：这些靶点不能写成因果卖点，详见 Discussion 段落。")


# ============================ 页签 4：环境自检 ============================
def page_env():
    st.subheader("环境自检（纠正 TOOL_MISSING 误判）")
    st.caption("很多「工具缺失」其实是装在非 PATH 目录（如 D:/MDtools）或 WSL 里。本页直接探测真实路径。")
    if st.button("运行探测"):
        rc, out = run_stream([sys.executable, os.path.join(SCRIPTS, "np_envcheck.py")])
        st.code(out[-4000:], language="text")
    st.divider()
    st.markdown("**单模块自检**：下面每个按钮都会真实执行对应模块的 selftest 并回显结果。")
    mods = ["np_mr", "np_synergy", "np_benchmark", "np_intervention",
            "np_contraindication", "np_decoy", "np_alias", "np_pocket", "np_finalcheck"]
    cols = st.columns(3)
    for i, m in enumerate(mods):
        if cols[i % 3].button(m, use_container_width=True):
            rc, out = run_stream([sys.executable, os.path.join(SCRIPTS, m + ".py"), "selftest"])
            st.code("%s（exit %d）\n%s" % (m, rc, out[-900:]), language="text")


# ============================ 主入口 ============================
def main():
    if st is None:
        print("未安装 streamlit。请先：pip install streamlit 然后：streamlit run app.py")
        return 1
    st.set_page_config(page_title="network-pharmacology-pro", layout="wide")
    st.title("network-pharmacology-pro")
    st.caption("证据感知 · 因果导向的中药网络药理学计算平台　|　"
               "安全前置 → 别名/标准化 → 成分 → 交集 → PPI → 零模型 → 富集 → 同源 → 门禁报告")
    with st.sidebar:
        st.markdown("### 运行设置")
        st.write("输出根目录：`%s`" % WORKDIR)
        st.write("配置文件：`config/pipeline.yaml`")
        st.write("引擎版本：%s" % (np_pipeline.load_config(CONFIG)[0].get("version", "?")))
        st.divider()
        st.markdown("**诚实边界**")
        st.caption("本平台输出的是假说，不是疗效结论。所有图表均可追溯到磁盘产物；"
                   "未跑门禁的报告一律盖 UNVERIFIED 水印。")
        if st.button("清理本次会话输出目录", use_container_width=True):
            if os.path.isdir(WORKDIR):
                shutil.rmtree(WORKDIR, ignore_errors=True)
            st.success("已清理 %s" % WORKDIR)

    t1, t2, t3, t4 = st.tabs(["一键分析", "基准测试", "因果校验", "环境自检"])
    with t1:
        page_analyze()
    with t2:
        page_benchmark()
    with t3:
        page_causal()
    with t4:
        page_env()
    return 0


if __name__ == "__main__":
    # 注意：Streamlit 会以 __main__ 执行本文件，这里**不能** sys.exit()，
    # 否则脚本被 SystemExit 中断，前端会报 "missing ScriptRunContext" 之类的怪错。
    main()
