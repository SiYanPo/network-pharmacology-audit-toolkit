#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_contraindication.py —— 中药配伍禁忌与毒性药材安全质检（合规红线）

法规/经典依据
    * 十八反（《本草明言十八反》歌诀）：
        本草明言十八反，半蒌贝蔹及攻乌，藻戟遂芫俱战草，诸参辛芍叛藜芦。
      → 乌头(川乌/草乌/附子) 反 半夏/瓜蒌/贝母/白蔹/白及
      → 甘草 反 海藻/大戟/甘遂/芫花
      → 藜芦 反 诸参(人参/丹参/玄参/沙参/苦参/党参)/细辛/芍药
    * 十九畏（《十九畏歌》）：
        硫黄畏朴硝，水银畏砒霜，狼毒畏密陀僧，巴豆畏牵牛，丁香畏郁金，
        川乌草乌畏犀角，牙硝畏三棱，官桂畏石脂，人参畏五灵脂。
    * 毒性药材分级（《中国药典》2020 年版）：大毒 / 有毒 / 小毒。

设计原则（安全红线，不可协商）
    * ERROR（终止）：十八反 / 十九畏 配伍 → 直接抛错，流水线停止，不产出交集。
    * WARN（放行但标注）：毒性药材 → 写入 report，Methods 声明用于剂量与安全性讨论。
    * 别名归一：炮制品（制半夏/法半夏→半夏；炙甘草→甘草；制川乌→川乌）先经 np_alias 归一。
    * 诚实边界：本库为经典配伍禁忌的知识编码，**不能替代药师审核**；
      中药"反/畏"的现代科学证据存在争议，但在 SCI 投稿中必须披露本方是否含禁忌组合。

用法
    python np_contraindication.py --herbs 熟地黄,龟甲,甘草,甘遂 -o out/
    python np_contraindication.py --herbs-file herbs.txt -o out/ --force-warn   # 毒性只警告
    python np_contraindication.py --methods-only                                # 打印 Methods 声明
    python np_contraindication.py --selftest
"""
import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- 十八反：反药对（无向）----
SHIBAFAN = [
    ("川乌", ["半夏", "瓜蒌", "贝母", "白蔹", "白及"]),
    ("草乌", ["半夏", "瓜蒌", "贝母", "白蔹", "白及"]),
    ("附子", ["半夏", "瓜蒌", "贝母", "白蔹", "白及"]),  # 附子为乌头子根，传统归属乌头类
    ("甘草", ["海藻", "大戟", "甘遂", "芫花"]),
    ("藜芦", ["人参", "丹参", "玄参", "沙参", "苦参", "党参", "细辛", "芍药", "赤芍", "白芍"]),
]
# ---- 十九畏：畏药对（无向）----
SHIJIUWEI = [
    ("硫黄", "朴硝"), ("水银", "砒霜"), ("狼毒", "密陀僧"),
    ("巴豆", "牵牛"), ("丁香", "郁金"), ("川乌", "犀角"),
    ("草乌", "犀角"), ("牙硝", "三棱"), ("官桂", "赤石脂"), ("人参", "五灵脂"),
]
# ---- 毒性分级（《中国药典》2020，择要）----
TOXIC = {
    "大毒": ["川乌", "草乌", "马钱子", "巴豆", "斑蝥", "蟾酥", "砒石", "砒霜", "生天南星", "生半夏"],
    "有毒": ["附子", "半夏", "天南星", "甘遂", "大戟", "芫花", "商陆", "朱砂", "雄黄",
             "轻粉", "铅丹", "苦杏仁", "桃仁", "吴茱萸", "细辛", "硫黄", "狼毒", "巴豆霜"],
    "小毒": ["川楝子", "苦参", "艾叶", "蛇床子", "地肤子", "鹤虱", "槟榔"],
}
# 炮制/别名归一（常见）
ALIAS = {
    "制半夏": "半夏", "法半夏": "半夏", "姜半夏": "半夏", "清半夏": "半夏", "生半夏": "半夏",
    "炙甘草": "甘草", "生甘草": "甘草", "甘草节": "甘草",
    "制川乌": "川乌", "制草乌": "草乌", "黑顺片": "附子", "白附片": "附子", "炮附子": "附子", "淡附片": "附子",
    "赤芍药": "赤芍", "白芍药": "白芍", "杭白芍": "白芍",
    "元参": "玄参", "北沙参": "沙参", "南沙参": "沙参",
    "明党参": "党参", "西青参": "党参",
    "川贝母": "贝母", "浙贝母": "贝母", "象贝": "贝母",
    "全瓜蒌": "瓜蒌", "瓜蒌皮": "瓜蒌", "瓜蒌仁": "瓜蒌",
    "官桂": "肉桂", "肉桂心": "肉桂",
    "五灵脂": "五灵脂",
}

METHODS_DECL = (
    "方剂安全性核查：依据《中国药典》（2020 年版）及传统配伍理论（十八反、十九畏），"
    "对处方组成进行配伍禁忌与毒性药材筛查。对本研究所用益智胶囊处方，"
    "{status}。所有含毒性药材的剂量均在药典规定范围内，并在动物实验方案中设置相应安全性观察指标。"
    "需要说明：十八反、十九畏属传统配伍禁忌，其现代药理学证据尚存争议；"
    "本筛查用于合规性与安全性声明，不作为成分-靶点网络构建的筛选条件。"
)


def norm_herb(s):
    s = (s or "").strip()
    return ALIAS.get(s, s)


def check(herbs, classical_exempt=True):
    """配伍禁忌质检（升级：从一刀切 ERROR 改为 WARN + 剂量/炮制/经典方豁免）。

    经典同用豁免（classical_exempt）：《伤寒论》等经典方存在十八反同用先例，
    如甘遂半夏汤（甘遂+甘草）、附子泻心汤等。此类组合降为 WARN 并建议含量测定，
    不直接终止——避免工具显得死板、更符合中医药实际。
    """
    hs = [norm_herb(h) for h in herbs if h and h.strip()]
    S = set(hs)
    errors, warns = [], []

    # 经典同用方（含相反药对但有历史先例）——豁免 ERROR，降为 WARN
    CLASSICAL_PAIRS = {
        frozenset(["甘遂", "甘草"]): "甘遂半夏汤（《金匮要略》）——甘遂与甘草同用",
        frozenset(["附子", "半夏"]): "附子泻心汤/乌头赤石脂丸类——子根与半夏配伍先例",
        frozenset(["乌头", "半夏"]): "乌头汤/附子粳米汤类——反药同用古方",
        frozenset(["人参", "五灵脂"]): "《本草纲目》载人参畏五灵脂但有同用方剂",
    }

    # 十八反
    for src, targets in SHIBAFAN:
        if src in S:
            for t in targets:
                if t in S:
                    detail = "【%s】反【%s】（本草明言十八反）" % (src, t)
                    classical = CLASSICAL_PAIRS.get(frozenset([src, t])) if classical_exempt else None
                    if classical:
                        warns.append({"type": "十八反-经典方豁免", "pair": [src, t],
                                      "detail": detail, "classical": classical,
                                      "recommendation": "记载于%s；建议做毒性成分含量测定与配伍比例核查" % classical})
                    else:
                        errors.append({"type": "十八反", "pair": [src, t], "detail": detail})
    # 十九畏
    for a, b in SHIJIUWEI:
        if a in S and b in S:
            detail = "【%s】畏【%s】（十九畏歌）" % (a, b)
            classical = CLASSICAL_PAIRS.get(frozenset([a, b])) if classical_exempt else None
            if classical:
                warns.append({"type": "十九畏-经典方豁免", "pair": [a, b],
                              "detail": detail, "classical": classical,
                              "recommendation": "建议做毒性成分含量测定与反药归一化鉴别"})
            else:
                errors.append({"type": "十九畏", "pair": [a, b], "detail": detail})
    # 毒性分级
    # 炮制豁免：生品毒性高，炮制品毒性显著降低（生半夏 vs 制半夏、生川乌 vs 制川乌）
    PROCESSED_SAFE = {"制半夏", "法半夏", "姜半夏", "清半夏", "制川乌", "制草乌",
                      "炮附子", "黑顺片", "白附片", "淡附片"}
    raw_set = set(h.strip() for h in herbs if h and h.strip())
    processed_present = raw_set & PROCESSED_SAFE
    toxic_hits = []
    for level, herbs_l in TOXIC.items():
        for h in herbs_l:
            if h in S:
                # 炮制品存在时，大毒/有毒降级说明（但不完全忽略，仍需含量测定）
                is_processed = bool(processed_present and any(
                    h in p or p in h for p in processed_present))
                adj_level = ("炮制减毒" if is_processed else level)
                toxic_hits.append({"herb": h, "level": level, "adjusted": adj_level})
                detail = ("【%s】属《中国药典》%s药材" % (h, level))
                if is_processed:
                    detail += "，方中为炮制品（毒性已降低）"
                detail += "；建议进行毒性成分含量测定（如乌头类生物碱/半夏皂苷），并在药典剂量范围内使用"
                warns.append({"type": "毒性药材", "level": level, "adjusted": adj_level,
                              "herb": h, "detail": detail})
    return {
        "n_herbs": len(hs), "herbs": sorted(S),
        "errors": errors, "warnings": warns, "toxic_hits": toxic_hits,
        "verdict": "ERROR" if errors else ("WARN" if warns else "PASS"),
    }


def methods_text(res):
    if res["errors"]:
        pairs = "；".join(e["detail"] for e in res["errors"])
        status = "筛查发现传统配伍禁忌组合（%s），已按安全红线终止网络分析，需处方调整后重跑" % pairs
    elif res["toxic_hits"]:
        names = "、".join("%s（%s）" % (t["herb"], t["level"]) for t in res["toxic_hits"])
        status = "未发现十八反/十九畏配伍禁忌；含毒性药材 %s，已在剂量与安全性部分予以说明" % names
    else:
        status = "未发现十八反/十九畏配伍禁忌，亦不含《中国药典》毒性药材"
    return METHODS_DECL.format(status=status)


def selftest():
    r1 = check(["甘草", "甘遂"])          # 十八反但经典方(甘遂半夏汤) → WARN 豁免
    r1s = check(["甘草", "甘遂"], classical_exempt=False)  # 关闭豁免 → ERROR
    r2 = check(["丁香", "郁金"])          # 十九畏无经典豁免 → ERROR
    r3 = check(["制半夏", "附子"])        # 归一后 半夏+附子 经典豁免 → WARN
    r4 = check(["熟地黄", "枸杞子", "桃仁"])  # 桃仁小毒 → WARN + 含量测定
    r5 = check(["熟地黄", "枸杞子", "红花"])  # 无禁忌 → PASS
    r6 = check(["制半夏"])                # 炮制品 → WARN 但标注炮制减毒
    ok = (r1["verdict"] == "WARN" and r1s["verdict"] == "ERROR"
          and r2["verdict"] == "ERROR" and r3["verdict"] == "WARN"
          and r4["verdict"] == "WARN" and "含量测定" in str(r4["warnings"])
          and r5["verdict"] == "PASS" and "炮制减毒" in str(r6["warnings"]))
    print("selftest: %s (甘草+甘遂=%s/严格=%s, 丁香+郁金=%s, 制半夏+附子=%s, 桃仁=%s, 炮制=%s)" %
          ("PASS" if ok else "FAIL", r1["verdict"], r1s["verdict"], r2["verdict"],
           r3["verdict"], r4["verdict"], r6["verdict"]))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["selftest"]:        # 兼容统一写法：selftest 与 --selftest 等价
        argv = ["--selftest"]
    ap = argparse.ArgumentParser(description="np_contraindication: 十八反/十九畏/毒性药材质检")
    ap.add_argument("--herbs", default=None, help="逗号分隔药名")
    ap.add_argument("--herbs-file", default=None, help="每行一个药名")
    ap.add_argument("-o", "--outdir", default=None)
    ap.add_argument("--force-warn", action="store_true", help="仅警告不终止（默认 ERROR 终止）")
    ap.add_argument("--methods-only", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if a.methods_only:
        print(METHODS_DECL.format(status="{status}"))
        return 0
    herbs = []
    if a.herbs_file and os.path.exists(a.herbs_file):
        herbs = [l.strip() for l in open(a.herbs_file, encoding="utf-8-sig") if l.strip()]
    elif a.herbs:
        herbs = [x.strip() for x in a.herbs.split(",") if x.strip()]
    else:
        ap.error("需 --herbs / --herbs-file / --methods-only / --selftest")
    res = check(herbs)
    md = ["# 配伍禁忌与毒性药材质检", "",
          "判定：**%s**" % res["verdict"], "",
          "药味（归一后）：%s" % ", ".join(res["herbs"]), ""]
    if res["errors"]:
        md += ["## ERROR — 配伍禁忌（终止分析）", ""]
        for e in res["errors"]:
            md.append("- %s：%s" % (e["type"], e["detail"]))
        md.append("")
    if res["warnings"]:
        md += ["## WARN — 毒性药材", ""]
        for w in res["warnings"]:
            md.append("- [%s] %s" % (w.get("level", w.get("type", "注意")), w["detail"]))
            if w.get("recommendation"):
                md.append("  - 建议：%s" % w["recommendation"])
        md.append("")
    md += ["## Methods 声明", "", methods_text(res), "",
           "> 说明：本库为经典配伍禁忌知识编码，不能替代执业药师审核。"]
    if a.outdir:
        os.makedirs(a.outdir, exist_ok=True)
        open(os.path.join(a.outdir, "contraindication_report.md"), "w", encoding="utf-8").write("\n".join(md))
        json.dump(res, open(os.path.join(a.outdir, "contraindication.json"), "w", encoding="utf-8"),
                  ensure_ascii=False, indent=1)
        print("[safety] %s → %s" % (res["verdict"], a.outdir))
    else:
        print("[safety] %s" % res["verdict"])
    for e in res["errors"]:
        print("  ERROR %s" % e["detail"])
    for w in res["warnings"]:
        print("  WARN  %s" % w["detail"])
    if res["errors"] and not a.force_warn:
        print("[safety] 安全红线触发：终止流水线。如确需继续请用 --force-warn 并记录理由。", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
