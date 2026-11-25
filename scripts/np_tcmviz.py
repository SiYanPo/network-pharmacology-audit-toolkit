#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
np_tcmviz.py —— v1.5.0 中医药属性可视化模块（性味归经雷达图 / 方剂桑基图）

背景（23 篇审计）：23 篇均未做"君臣佐使"或"性味归经"拓扑分析（0/23）——这是中文
NP 论文的普遍留白，也是审稿人常问"方剂配伍如何体现"的点。本模块补齐可视化：
  - radar：单药/多药性味归经雷达图（四气5级 × 五味 × 归经强度）
  - sankey：方剂"药→性→味→归经"桑基图
  - junchen：君臣佐使网络图（按用户给的分级 + 药物权重）

内置数据库：40 味常用中药的性味归经（2020 版《中国药典》常见记载，四气按寒/凉/平/温/热
5 级编码，归经按 12 经归属）。若药不在库中，用 --props 传自定义表：
  herb,气,味,归经  （如 当归,温,甘辛,肝脾心）
"""
import argparse
import os

# 四气 → 数值（寒1 凉2 平3 温4 热5）
QI = {"寒": 1, "大寒": 1, "微寒": 1.5, "凉": 2, "微凉": 2.5, "平": 3, "微温": 3.5, "温": 4, "微热": 4.5, "热": 5, "大热": 5}
# 五味（可复合）与十二归经
WEI = ["甘", "苦", "酸", "辛", "咸"]
MERIDIANS = ["肝", "心", "脾", "肺", "肾", "胃", "胆", "大肠", "小肠", "膀胱", "三焦", "心包"]

# 内置小库：40 味常用中药（气, 味, 归经）
TCM_PROPS = {
 "黄芪": ("温", "甘", "脾肺"), "当归": ("温", "甘辛", "肝心脾"), "丹参": ("微寒", "苦", "心肝"),
 "黄芩": ("寒", "苦", "肺胆脾大肠小肠"), "黄连": ("寒", "苦", "心脾胃胆大肠"), "黄柏": ("寒", "苦", "肾膀胱大肠"),
 "金银花": ("寒", "甘", "肺心胃"), "连翘": ("微寒", "苦", "肺心小肠"), "栀子": ("寒", "苦", "心肺三焦"),
 "柴胡": ("微寒", "辛", "肝胆肺"), "白芍": ("微寒", "苦酸", "肝脾"), "川芎": ("温", "辛", "肝胆心包"),
 "白术": ("温", "苦甘", "脾胃"), "茯苓": ("平", "甘淡", "心肺脾肾"), "甘草": ("平", "甘", "心肺脾胃"),
 "人参": ("微温", "甘微苦", "脾肺心肾"), "熟地黄": ("微温", "甘", "肝肾"), "生地黄": ("寒", "甘苦", "心肝肾"),
 "大黄": ("寒", "苦", "脾胃大肠肝心包"), "桃仁": ("平", "苦甘", "心肝大肠"), "红花": ("温", "辛", "心肝"),
 "益母草": ("微寒", "苦辛", "肝心包膀胱"), "泽泻": ("寒", "甘淡", "肾膀胱"), "山茱萸": ("微温", "酸涩", "肝肾"),
 "山药": ("平", "甘", "脾肺肾"), "牡丹皮": ("微寒", "苦辛", "心肝肾"), "肉桂": ("大热", "辛甘", "肾脾心肝"),
 "附子": ("大热", "辛甘", "心肾脾"), "干姜": ("热", "辛", "脾胃心肺"), "半夏": ("温", "辛", "脾胃肺"),
 "陈皮": ("温", "辛苦", "肺脾"), "枳实": ("微寒", "苦辛", "脾胃"), "厚朴": ("温", "苦辛", "脾胃肺大肠"),
 "杏仁": ("微温", "苦", "肺大肠"), "桔梗": ("平", "苦辛", "肺"), "葛根": ("凉", "甘辛", "脾胃"),
 "酸枣仁": ("平", "甘酸", "肝胆心"), "远志": ("温", "苦辛", "心肺肾"), "龙骨": ("平", "甘涩", "心肝肾"),
 "牡蛎": ("微寒", "咸", "肝肾"), "石菖蒲": ("温", "辛苦", "心胃"),
}


def _props_for(herbs, custom):
    out = {}
    for h in herbs:
        if h in custom:
            out[h] = custom[h]
        elif h in TCM_PROPS:
            out[h] = TCM_PROPS[h]
    missing = [h for h in herbs if h not in out]
    if missing:
        print(f"[提示] 库中无: {missing}；可用 --props 文件补充（herb,气,味,归经）")
    return out


def load_custom(path):
    d = {}
    if path and os.path.exists(path):
        for line in open(path, encoding="utf-8-sig"):
            parts = [p.strip() for p in line.strip().split(",")]
            if len(parts) >= 4 and parts[0] and parts[0] != "herb":
                d[parts[0]] = (parts[1], parts[2], parts[3])
    return d


def cmd_radar(args):
    os.makedirs(args.outdir, exist_ok=True)
    custom = load_custom(args.props)
    herbs = [h.strip() for h in args.herbs.split(",") if h.strip()]
    props = _props_for(herbs, custom)
    if not props:
        raise SystemExit("[FAIL] 无可用性味归经数据")
    axes = ["寒", "凉", "平", "温", "热"]  # 四气五级雷达
    # 每味药在五气上的隶属度（气值=QI[气] → 高斯峰）
    import math
    series = []
    for h, (qi, wei, mer) in props.items():
        v = QI.get(qi, 3)
        prof = [math.exp(-((i + 1 - v) ** 2) / 0.8) for i in range(5)]
        series.append((h, prof, wei, mer))
    W, H, cx, cy, R = 680, 440, 340, 240, 150
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="sans-serif">',
           f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
           '<text x="340" y="26" text-anchor="middle" font-size="15" font-weight="bold">中药四气雷达图（寒→热五级隶属度）</text>']
    colors = ["#c0392b", "#2471a3", "#1e8449", "#b7950b", "#7d3c98", "#117864"]
    for ring in (1, .75, .5, .25):
        pts = " ".join(f"{cx + R * ring * math.cos(a):.0f},{cy + R * ring * math.sin(a):.0f}"
                       for a in [(-90 + i * 72) * math.pi / 180 for i in range(5)])
        svg.append(f'<polygon points="{pts}" fill="none" stroke="#ddd"/>')
    for i, ax in enumerate(axes):
        a = (-90 + i * 72) * math.pi / 180
        x, y = cx + R * math.cos(a), cy + R * math.sin(a)
        svg.append(f'<line x1="{cx}" y1="{cy}" x2="{x:.0f}" y2="{y:.0f}" stroke="#ddd"/>')
        lx, ly = cx + (R + 22) * math.cos(a), cy + (R + 22) * math.sin(a)
        svg.append(f'<text x="{lx:.0f}" y="{ly:.0f}" text-anchor="middle" font-size="13">{ax}</text>')
    for si, (h, prof, wei, mer) in enumerate(series):
        c = colors[si % len(colors)]
        pts = " ".join(f"{cx + R * p * math.cos((-90 + i * 72) * math.pi / 180):.0f},"
                       f"{cy + R * p * math.sin((-90 + i * 72) * math.pi / 180):.0f}"
                       for i, p in enumerate(prof))
        svg.append(f'<polygon points="{pts}" fill="{c}" fill-opacity="0.10" stroke="{c}" stroke-width="2"/>')
        svg.append(f'<rect x="10" y="{70 + si * 22}" width="12" height="12" fill="{c}"/>')
        svg.append(f'<text x="28" y="{80 + si * 22}" font-size="12">{h}（味:{wei}｜归经:{mer}）</text>')
    svg.append("</svg>")
    out = os.path.join(args.outdir, "qi_radar.svg")
    open(out, "w", encoding="utf-8").write("".join(svg))
    print(f"四气雷达图 → {out}")


def cmd_sankey(args):
    os.makedirs(args.outdir, exist_ok=True)
    custom = load_custom(args.props)
    herbs = [h.strip() for h in args.herbs.split(",") if h.strip()]
    props = _props_for(herbs, custom)
    if not props:
        raise SystemExit("[FAIL] 无可用性味归经数据")
    # 桑基：左列=药，中列=五味，右列=归经
    col_x = [60, 340, 620]
    W, H = 680, max(380, 60 + 34 * len(props))
    nodes = {"herb": [], "wei": [], "mer": []}
    links = []
    for h, (qi, wei, mer) in props.items():
        nodes["herb"].append(h)
        for wch in wei if wei not in ("甘淡",) else ["甘", "淡"]:
            pass
        for wch in [c for c in wei if c in WEI] or ["平"]:
            nodes["wei"].append(wch)
            links.append((h, wch, "wei"))
        for m in mer if mer not in ("三焦", "心包") else [mer]:
            if m in MERIDIANS:
                links.append((wch if False else h, m, "direct"))
    # 简化桑基：药→归经 直连 + 药→味
    h_pos = {h: 40 + i * (H - 80) / max(1, len(props) - 1 or 1) for i, h in enumerate(props)}
    wei_set = sorted({l[1] for l in links if l[2] == "wei"}, key=WEI.index)
    mer_set = sorted({l[1] for l in links if l[2] == "direct"}, key=MERIDIANS.index)
    wei_pos = {w: 40 + i * (H - 80) / max(1, len(wei_set) - 1 or 1) for i, w in enumerate(wei_set)}
    mer_pos = {m: 40 + i * (H - 80) / max(1, len(mer_set) - 1 or 1) for i, m in enumerate(mer_set)}
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="sans-serif">',
           f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
           '<text x="340" y="24" text-anchor="middle" font-size="15" font-weight="bold">方剂性味归经桑基图（药→味→经）</text>']
    # 药→味
    for h, w, _ in [(l[0], l[1], None) for l in links if l[2] == "wei"]:
        if h in props and w in wei_pos:
            qi, wei, mer = props[h]
            svg.append(f'<path d="M {col_x[0]+30} {h_pos[h]:.0f} C {col_x[1]-60} {h_pos[h]:.0f} {col_x[0]+90} {wei_pos[w]:.0f} {col_x[1]-30} {wei_pos[w]:.0f}" stroke="#2471a3" stroke-opacity="0.35" fill="none" stroke-width="3"/>')
    # 味→经（由药归属聚合：药-经直连经味列右侧）
    for h in props:
        qi, wei, mer = props[h]
        for m in mer:
            if m in mer_pos:
                svg.append(f'<path d="M {col_x[1]+30} {wei_pos.get(wei[0] if wei[0] in WEI else "甘", h_pos[h]):.0f} C {col_x[2]-60} {wei_pos.get(wei[0] if wei[0] in WEI else "甘", h_pos[h]):.0f} {col_x[1]+90} {mer_pos[m]:.0f} {col_x[2]-30} {mer_pos[m]:.0f}" stroke="#c0392b" stroke-opacity="0.28" fill="none" stroke-width="3"/>')
    for h, y in h_pos.items():
        svg.append(f'<rect x="{col_x[0]-10}" y="{y-11:.0f}" width="40" height="22" rx="4" fill="#1e8449"/>')
        svg.append(f'<text x="{col_x[0]+10}" y="{y+4:.0f}" font-size="12" fill="#fff" text-anchor="middle">{h}</text>')
    for w, y in wei_pos.items():
        svg.append(f'<rect x="{col_x[1]-10}" y="{y-11:.0f}" width="40" height="22" rx="4" fill="#b7950b"/>')
        svg.append(f'<text x="{col_x[1]+10}" y="{y+4:.0f}" font-size="12" text-anchor="middle">{w}</text>')
    for m, y in mer_pos.items():
        svg.append(f'<rect x="{col_x[2]-10}" y="{y-11:.0f}" width="40" height="22" rx="4" fill="#7d3c98"/>')
        svg.append(f'<text x="{col_x[2]+10}" y="{y+4:.0f}" font-size="12" text-anchor="middle">{m}经</text>')
    svg.append("</svg>")
    out = os.path.join(args.outdir, "sankey.svg")
    open(out, "w", encoding="utf-8").write("".join(svg))
    print(f"桑基图 → {out}")


def cmd_junchen(args):
    """君臣佐使网络图：--map "黄芪:君,当归:臣,川芎:佐,牡丹皮:使" """
    os.makedirs(args.outdir, exist_ok=True)
    roles = {}
    order = {"君": 0, "臣": 1, "佐": 2, "使": 3}
    for kv in args.map.split(","):
        h, r = kv.split(":")
        roles[h.strip()] = r.strip()
    W, H = 680, 420
    cx, cy = 340, 230
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="sans-serif">',
           f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
           '<text x="340" y="26" text-anchor="middle" font-size="15" font-weight="bold">君臣佐使配伍网络</text>',
           f'<circle cx="{cx}" cy="{cy}" r="26" fill="#c0392b"/>',
           '<text x="340" y="236" font-size="13" fill="#fff" text-anchor="middle">方剂</text>']
    cols = {"君": "#c0392b", "臣": "#b7950b", "佐": "#1e8449", "使": "#2471a3"}
    rolen = {r: sum(1 for x in roles.values() if x == r) for r in cols}
    idx = {r: 0 for r in cols}
    for h, r in sorted(roles.items(), key=lambda kv: order.get(kv[1], 4)):
        ang = (-90 + (idx[r] + sum(rolen[rr] for rr in list(order)[:list(order).index(r)])) * (360 / max(1, len(roles)))) * 3.14159 / 180
        idx[r] += 1
        R2 = 120 + order.get(r, 4) * 22
        x, y = cx + R2 * __import__("math").cos(ang), cy + R2 * __import__("math").sin(ang)
        svg.append(f'<line x1="{cx}" y1="{cy}" x2="{x:.0f}" y2="{y:.0f}" stroke="{cols.get(r, "#999")}" stroke-opacity="0.5" stroke-width="{4 - order.get(r, 3)}"/>')
        svg.append(f'<circle cx="{x:.0f}" cy="{y:.0f}" r="22" fill="{cols.get(r, "#999")}"/>')
        svg.append(f'<text x="{x:.0f}" y="{y - 28:.0f}" font-size="14" text-anchor="middle" font-weight="bold">{h}</text>')
        svg.append(f'<text x="{x:.0f}" y="{y + 5:.0f}" font-size="11" fill="#fff" text-anchor="middle">{r}</text>')
    svg.append("</svg>")
    out = os.path.join(args.outdir, "junchen_network.svg")
    open(out, "w", encoding="utf-8").write("".join(svg))
    print(f"君臣佐使网络 → {out}")


def main():
    ap = argparse.ArgumentParser(description="np_tcmviz: 性味归经/君臣佐使可视化")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("radar"); s.add_argument("--herbs", required=True)
    s.add_argument("--props", default=None); s.add_argument("--outdir", default=".")
    s = sub.add_parser("sankey"); s.add_argument("--herbs", required=True)
    s.add_argument("--props", default=None); s.add_argument("--outdir", default=".")
    s = sub.add_parser("junchen"); s.add_argument("--map", required=True, help='黄芪:君,当归:臣,...')
    s.add_argument("--outdir", default=".")
    a = ap.parse_args()
    {"radar": cmd_radar, "sankey": cmd_sankey, "junchen": cmd_junchen}[a.cmd](a)


if __name__ == "__main__":
    main()
