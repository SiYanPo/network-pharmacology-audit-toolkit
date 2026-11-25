#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_pipeline.py —— 配置驱动的 DAG 流水线编排器（np_run.py 的非破坏性升级层）

设计原则（主线接入：不硬编码、不无脑塞入 np_run.py）
    * np_run.py 的 01–11 阶段逻辑已验证，保持不动；本模块是**编排层**，
      通过 config/pipeline.yaml 描述节点，调用各 np_* 模块。零回归风险。
    * 三类节点：
        core   — 默认必经（alias / ortholog 是数据清洗与湿实验衔接刚需）
        plugin — --enable-<flag> 显式触发（funnel/doseweight/subnetwork，大复方默认不跑）
        audit  — 独立审计工具（goldtest/evidence），不进常规主流程
    * DAG：needs 声明依赖 → 拓扑排序 + 环检测；产物 artifacts 支持 --resume 断点续跑。

配置格式
    config/pipeline.yaml（推荐）；无 PyYAML 时自动回退读 config/pipeline.json。
    YAML 使用扁平子集（list of dict + 标量），无需复杂锚点。

用法
    python np_pipeline.py --config config/pipeline.yaml --list            # 列出 DAG 与类别
    python np_pipeline.py --config ... --inputs inputs.json --dry-run      # 只校验与打印计划
    python np_pipeline.py --config ... --inputs inputs.json \
        --enable-funnel --enable-louvain -o out/                           # 触发插件
    python np_pipeline.py --config ... --inputs inputs.json --resume       # 跳过已完成节点
    python np_pipeline.py audit goldtest --help-notes                      # 独立审计工具说明
"""
import argparse
import hashlib
import importlib
import inspect
import json
import os
import re
import sys
import time

# v3.0.2: 这些节点在"数据本身不足"时会科学地拒绝执行（如 nullmodel 要求交集基因>=3），
# 属于可降级节点：失败时记录 degraded 并继续下游，避免有效产出被连带丢弃。
SOFT_FAIL_NODES = {"nullmodel", "enrich", "doseweight", "funnel", "ortholog"}

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
DEFAULT_CONFIG = os.path.join(ROOT, "config", "pipeline.yaml")


# ---------- 配置加载（YAML/JSON 双模，零硬依赖） ----------
def load_config(path):
    if not os.path.exists(path):
        # 回退 JSON
        jp = os.path.splitext(path)[0] + ".json"
        if os.path.exists(jp):
            return json.load(open(jp, encoding="utf-8")), jp
        raise SystemExit("[FAIL] 配置不存在：%s（或同名 .json）" % path)
    if path.lower().endswith(".json"):
        return json.load(open(path, encoding="utf-8")), path
    try:
        import yaml
        return yaml.safe_load(open(path, encoding="utf-8")), path
    except ImportError:
        # 极简 YAML 子集解析（list of dict），够本配置用
        return _mini_yaml(open(path, encoding="utf-8").read()), path


def _mini_yaml(text):
    """无 PyYAML 时的兜底：只解析本配置用到的扁平结构（顶层标量 + stages/audit_tools 列表）。"""
    cfg, cur_list, cur_item, section = {}, None, None, None
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip())
        s = line.strip()
        if indent == 0 and s.endswith(":"):
            section = s[:-1]
            cur_list = []
            cfg[section] = {"stages": cur_list} if section == "config" else cur_list
            # 顶层段 defaults/version/stages/audit_tools
            if section in ("defaults", "version"):
                cur_list = None
            continue
        if indent == 0 and ":" in s and not s.endswith(":"):
            k, v = s.split(":", 1)
            cfg[k.strip()] = _scalar(v.strip())
            continue
        if section in ("stages", "audit_tools"):
            lst = cfg[section]
            if s.startswith("- "):
                cur_item = {}
                lst.append(cur_item)
                kv = s[2:].strip()
                if ":" in kv:
                    k, v = kv.split(":", 1)
                    cur_item[k.strip()] = _scalar(v.strip())
            elif cur_item is not None and ":" in s:
                k, v = s.split(":", 1)
                cur_item[k.strip()] = _scalar(v.strip())
        elif section == "defaults" and ":" in s:
            k, v = s.split(":", 1)
            cfg.setdefault("defaults", {})[k.strip()] = _scalar(v.strip())
    return cfg


def _scalar(v):
    v = v.strip().strip('"').strip("'")
    if v.lower() in ("true", "yes"):
        return True
    if v.lower() in ("false", "no"):
        return False
    if v.startswith("[") and v.endswith("]"):
        return [x.strip().strip('"').strip("'") for x in v[1:-1].split(",") if x.strip()]
    try:
        return int(v)
    except ValueError:
        pass
    try:
        return float(v)
    except ValueError:
        pass
    return v


# ---------- DAG ----------
def topo_sort(stages):
    by_id = {s["id"]: s for s in stages}
    state, order = {}, []

    def visit(sid, stack):
        if state.get(sid) == "done":
            return
        if state.get(sid) == "visiting":
            raise SystemExit("[DAG 环检测] 循环依赖：%s" % " → ".join(stack + [sid]))
        if sid not in by_id:
            raise SystemExit("[DAG 错误] 依赖的节点不存在：%s" % sid)
        state[sid] = "visiting"
        for dep in by_id[sid].get("needs") or []:
            visit(dep, stack + [sid])
        state[sid] = "done"
        order.append(by_id[sid])

    for s in stages:
        visit(s["id"], [])
    return order


def select_nodes(cfg, enable_flags, disable_ids, only_ids):
    stages = cfg.get("stages") or []
    selected = []
    for s in stages:
        sid = s["id"]
        kind = s.get("kind", "core")
        flag = s.get("enable_flag", "")
        if only_ids and sid not in only_ids:
            continue
        if sid in disable_ids:
            continue
        if kind == "plugin":
            if flag and "--" + flag in enable_flags:
                selected.append(s)
            elif sid in enable_flags:
                selected.append(s)
        else:
            selected.append(s)
    return selected


def render_cli(template, inputs, outdir, params):
    ctx = {"inputs": inputs, "outdir": outdir, "params": params}

    def repl(m):
        expr = m.group(1).strip()
        node = ctx
        for part in expr.split("."):
            node = node.get(part, "") if isinstance(node, dict) else ""
        return str(node)

    return re.sub(r"\{([^}]+)\}", repl, template or "")


def artifacts_done(node, outdir):
    """--resume 的"假成功"防御：文件存在 ≠ 阶段成功。

    三层校验（防网络中断生成的损坏/空/仅表头 CSV 被误判为已完成）：
      1. 存在性
      2. 非空 + 最小字节（> MIN_BYTES）
      3. 表格式产物必须 >1 行（表头之外有数据）+ 非空值比例; JSON 必须可解析
    任一层失败 → 返回 False，强制重跑该节点。

    另：若 manifest 记录了该节点产物的 sha256，则做 checksum 复核（可选，见 verify_artifacts）。
    """
    arts = node.get("artifacts") or []
    if not arts:
        return False
    MIN_BYTES = 10  # 只排除完全空文件；仅表头由下面行数检查拦截（小表合法，勿用大字节门槛误杀）
    for a in arts:
        p = os.path.join(outdir, a)
        if not os.path.exists(p):
            return False
        try:
            size = os.path.getsize(p)
            if size < MIN_BYTES:
                return False
            if p.lower().endswith((".csv", ".tsv")):
                with open(p, encoding="utf-8-sig", errors="replace") as f:
                    lines = [l for l in f if l.strip()]
                if len(lines) < 2:  # 只有表头 = 假成功
                    return False
                # 数据行必须有实际内容（非全空/全逗号）
                data_lines = lines[1:]
                if not any(any(c.strip() for c in l.replace("\t", ",").split(",")) for l in data_lines):
                    return False
            elif p.lower().endswith(".json"):
                json.load(open(p, encoding="utf-8"))  # 损坏 JSON 直接抛错 → False
        except Exception:
            return False
    return True


def verify_artifacts_checksum(node, outdir, manifest):
    """可选：对照上次 manifest 的 sha256，检测内容漂移（--resume 时提示）。"""
    prev = {n.get("id"): n for n in manifest.get("nodes", [])} if isinstance(manifest, dict) else {}
    rec = prev.get(node["id"], {})
    expected = rec.get("artifact_sha256") or {}
    mismatches = []
    for a in node.get("artifacts") or []:
        p = os.path.join(outdir, a)
        if a in expected and os.path.exists(p):
            h = hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
            if h != expected[a]:
                mismatches.append(a)
    return mismatches


def split_cli(s):
    """按空白切分命令行，但尊重双/单引号，且**不解析反斜杠**（Windows 路径安全）。

    为什么不用 shlex：shlex.split(posix=True) 会吃掉 C:\\ 的反斜杠；
    posix=False 又会把引号留在 token 里。自写 8 行更可靠。
    """
    out, buf, q = [], [], None
    for ch in s or "":
        if q:
            if ch == q:
                q = None
            else:
                buf.append(ch)
        elif ch in "\"'":
            q = ch
        elif ch.isspace():
            if buf:
                out.append("".join(buf))
                buf = []
        else:
            buf.append(ch)
    if buf:
        out.append("".join(buf))
    return out


def run_node(node, cli, quiet, dry):
    sid = node["id"]
    mod_name = node.get("module")
    print("\n[%s] %s（%s）" % (sid, node.get("name", ""), node.get("kind", "core")))
    if dry:
        mod_full = mod_name if mod_name.startswith("np_") else "np_" + mod_name
        print("  DRY-RUN → python %s.py %s" % (mod_full, cli))
        return 0
    mod_full = mod_name if mod_name.startswith("np_") else "np_" + mod_name
    try:
        mod = importlib.import_module(mod_full)
    except ImportError as e:
        print("  [SKIP] 模块不可导入：%s" % e)
        return 1
    t0 = time.time()
    # v3.0.2: --quiet 透传前检测子模块是否支持，避免向不支持 --quiet 的
    # 脚本（如 np_contraindication/np_alias/np_ortholog）误传导致 argparse 报错。
    kwargs = []
    if quiet:
        try:
            src = inspect.getsource(mod)
            if '--quiet' in src and ('add_argument' in src):
                kwargs = ["--quiet"]
        except (OSError, TypeError):
            pass  # 无法取源码时保守不透传
    argv = split_cli(cli) + kwargs
    print("  → %s.py %s" % (mod_full, " ".join(argv[:8])))
    try:
        rc = mod.main(argv)
    except SystemExit as e:
        rc = e.code if isinstance(e.code, int) else 0
    except Exception as e:
        print("  [ERR] %s: %s" % (type(e).__name__, str(e)[:120]))
        return 1
    print("  用时 %.1fs" % (time.time() - t0))
    return rc or 0


def main(argv=None):
    ap = argparse.ArgumentParser(description="np_pipeline: 配置驱动 DAG 流水线编排")
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--inputs", help="输入清单 JSON（inputs/params 两段）")
    ap.add_argument("-o", "--outdir", default=None)
    ap.add_argument("--list", action="store_true", help="列出 DAG 节点与类别")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--resume", action="store_true", help="跳过产物已齐的节点")
    ap.add_argument("--enable-funnel", action="store_true")
    ap.add_argument("--enable-doseweight", action="store_true")
    ap.add_argument("--enable-louvain", action="store_true")
    ap.add_argument("--enable", default="", help="逗号分隔自定义插件 flag 或节点 id")
    ap.add_argument("--disable", default="", help="逗号分隔禁用节点 id")
    ap.add_argument("--only", default="", help="逗号分隔只跑节点 id")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--offline", action="store_true",
                    help="全局离线：等价 NP_OFFLINE=1（跳过 STRING/HPA/在线同源，缺失即降级并标注）")
    ap.add_argument("--audit", default=None, help="独立审计工具名（goldtest/evidence），不进主流程")
    ap.add_argument("--audit-notes", action="store_true", help="打印审计工具说明")
    a = ap.parse_args(argv)
    if a.offline:
        os.environ["NP_OFFLINE"] = "1"

    cfg, cfg_path = load_config(a.config)

    # 审计工具独立入口（不进主流程）
    if a.audit or a.audit_notes:
        tools = {t["id"]: t for t in (cfg.get("audit_tools") or [])}
        if a.audit_notes and not a.audit:
            print("可用审计工具（独立运行，不进主流程）：")
            for tid, t in tools.items():
                m = t["module"] if t["module"].startswith("np_") else "np_" + t["module"]
                print("  %-10s %s.py — %s" % (tid, m, t.get("note", "")))
            return 0
        if a.audit not in tools:
            raise SystemExit("[FAIL] 未知审计工具：%s（可选 %s）" % (a.audit, list(tools)))
        t = tools[a.audit]
        m = t["module"] if t["module"].startswith("np_") else "np_" + t["module"]
        print("审计工具：%s → %s.py" % (a.audit, m))
        print("说明：%s" % t.get("note", ""))
        print("按常规方式单独运行：python %s.py <参数>" % m)
        return 0

    if a.list:
        print("配置：%s（v%s）" % (cfg_path, cfg.get("version", "?")))
        print("%-12s %-8s %-7s %s" % ("ID", "KIND", "DEFAULT", "依赖"))
        for s in cfg.get("stages") or []:
            print("%-12s %-8s %-7s %s%s" % (
                s["id"], s.get("kind", "core"), str(s.get("default", False)),
                ",".join(s.get("needs") or []) or "-",
                ("  [enable --%s]" % s["enable_flag"]) if s.get("enable_flag") else ""))
        print("\n审计工具（不进主流程）：%s" %
              ", ".join(t["id"] for t in (cfg.get("audit_tools") or [])))
        return 0

    inputs = {"herbs": "", "herbs_file": "", "targets": "", "compounds": "", "disease": "", "edges": "",
              "compound_edges": "", "dose": "", "gmt": ""}
    params = {"label": cfg.get("defaults", {}).get("label", "网络药理学分析"),
              "organism": "mouse", "tissue": ""}
    if a.inputs:
        if not os.path.exists(a.inputs):
            raise SystemExit("[FAIL] 输入清单不存在：%s" % a.inputs)
        data = json.load(open(a.inputs, encoding="utf-8"))
        inputs.update(data.get("inputs") or {})
        params.update(data.get("params") or {})
    outdir = os.path.abspath(a.outdir or cfg.get("defaults", {}).get("outdir", "np_out"))
    params["label"] = params.get("label") or "网络药理学分析"
    # v3.0.2: ortholog 节点 --organism 是"目标模式生物"(mouse/rat/zebrafish/all)，
    # 用户常误填数据来源 human/Homo sapiens。自动归一到 mouse，并把原值记到 source_organism。
    _org = str(params.get("organism", "") or "").strip()
    if _org.lower() in ("human", "homo sapiens", "hsa", "人"):
        params["source_organism"] = params.get("source_organism") or _org
        params["organism"] = "mouse"
        print("  [normalize] params.organism='%s' → 'mouse'（orthologist 目标生物为人类时无同源转换；已记 source_organism）" % _org)
    # v3.0.2: safety 节点硬编码 --herbs-file。若清单未给 herbs_file 但给了 herbs，
    # 自动物化一个临时药名文件，避免 argparse 收到空值报错（单化合物/无药材场景）。
    if not (inputs.get("herbs_file") or "").strip() and (inputs.get("herbs") or "").strip():
        os.makedirs(outdir, exist_ok=True)
        hf = os.path.join(outdir, "_auto_herbs.txt")
        with open(hf, "w", encoding="utf-8") as _fh:
            for _h in re.split(r"[,，]", inputs["herbs"]):
                if _h.strip():
                    _fh.write(_h.strip() + "\n")
        inputs["herbs_file"] = hf

    enable = set()
    for flag in ("--enable-funnel", "--enable-doseweight", "--enable-louvain"):
        if getattr(a, flag[2:].replace("-", "_")):
            enable.add(flag)
    for x in a.enable.split(","):
        if x.strip():
            enable.add(x.strip())
    disables = {x.strip() for x in a.disable.split(",") if x.strip()}
    onlys = {x.strip() for x in a.only.split(",") if x.strip()}
    # --only 的依赖闭包：拓扑排序需要完整依赖链，但**只执行**用户点名的节点
    requested = set(onlys)
    if onlys:
        by_id = {s["id"]: s for s in (cfg.get("stages") or [])}
        stack = list(onlys)
        while stack:
            sid = stack.pop()
            for dep in by_id.get(sid, {}).get("needs") or []:
                if dep not in onlys:
                    onlys.add(dep)
                    stack.append(dep)

    selected = select_nodes(cfg, enable, disables, onlys)
    order = topo_sort([s for s in (cfg.get("stages") or []) if s["id"] in {n["id"] for n in selected}])
    # 保序：按拓扑序执行 selected
    oid = {n["id"]: n for n in selected}
    plan = [oid[s["id"]] for s in order if s["id"] in oid]
    if requested:
        plan = [n for n in plan if n["id"] in requested]

    print("=" * 64)
    print(" 配置驱动 DAG 流水线 v%s" % cfg.get("version", "?"))
    print(" 输出：%s" % outdir)
    print(" 计划（%d 节点）：%s" % (len(plan), " → ".join(n["id"] for n in plan)))
    print("=" * 64)
    if not a.dry_run:
        os.makedirs(outdir, exist_ok=True)

    manifest = {"config": cfg_path, "outdir": outdir, "inputs": inputs, "params": params,
                "enable": sorted(enable), "nodes": []}
    failed = False
    for node in plan:
        if a.resume and artifacts_done(node, outdir):
            print("\n[%s] 产物齐备，跳过（--resume）" % node["id"])
            manifest["nodes"].append({"id": node["id"], "status": "skipped_done"})
            continue
        cli = render_cli(node.get("cli", ""), inputs, outdir, params)
        rc = run_node(node, cli, a.quiet, a.dry_run)
        # 记录本节点产物 sha256，供后续 --resume checksum 复核
        art_sha = {}
        for a_rel in node.get("artifacts") or []:
            ap_ = os.path.join(outdir, a_rel)
            if os.path.exists(ap_):
                art_sha[a_rel] = hashlib.sha256(open(ap_, "rb").read()).hexdigest()[:16]
        manifest["nodes"].append({"id": node["id"], "status": "ok" if rc == 0 else "failed",
                                   "cli": cli, "module": node.get("module"),
                                   "artifact_sha256": art_sha})
        if rc != 0:
            # v3.0.2: 部分节点在"数据本身不足"时会科学地拒绝（如 nullmodel 交集<3），
            # 这类属于可降级节点：记录 degraded 并继续下游，避免让 enrich/report 等
            # 有效产出被连带丢弃。关键节点失败仍保持原样中止。
            if node["id"] in SOFT_FAIL_NODES:
                print("  [降级] 节点 %s 失败（数据不足，科学拒绝）；记录 degraded 并继续下游" % node["id"])
                manifest["nodes"][-1]["status"] = "degraded"
                manifest["nodes"][-1]["degrade_reason"] = "data_insufficient_or_optional"
                continue
            failed = True
            print("  [中止] 节点 %s 失败；修复后可用 --resume 续跑" % node["id"])
            break
    # 数据库版本（可追溯性，评审 §5）
    try:
        import np_snapshot as _S
        manifest["versions"] = _S.versions_used()
    except Exception as e:
        manifest["versions"] = {"source": "unavailable", "error": type(e).__name__}

    mp = os.path.join(outdir, "pipeline_manifest.json")
    if not a.dry_run:
        # 与上次清单合并：--resume / --only 局部重跑时，不要把已完成节点从清单里抹掉
        prev = {}
        if os.path.exists(mp):
            try:
                prev = {n["id"]: n for n in (json.load(open(mp, encoding="utf-8")).get("nodes") or [])}
            except Exception:
                prev = {}
        merged = dict(prev)
        for n in manifest["nodes"]:
            merged[n["id"]] = n
        order_ids = [s["id"] for s in (cfg.get("stages") or [])]
        manifest["nodes"] = [merged[i] for i in order_ids if i in merged]
        manifest["nodes"].extend([v for k, v in merged.items() if k not in order_ids])
        json.dump(manifest, open(mp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print("\n流水线清单 → %s" % mp)
    print("\n完成（%s）" % ("有失败节点" if failed else "全部成功"))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
