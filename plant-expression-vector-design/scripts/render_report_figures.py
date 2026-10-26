"""Render bounded conceptual/report diagrams, never sequence maps or biological evidence."""
from __future__ import annotations

import argparse
import json
import re
from html import escape
from pathlib import Path

PALETTE = {
    "input": ("#cce5ff", "#3366cc"),
    "core": ("#ffd8a8", "#b85c00"),
    "output": ("#d3f9d8", "#2b8a3e"),
    "control": ("#eeeeee", "#666666"),
    "risk": ("#ffe1e1", "#b42318"),
}
STATUS = {
    "direct": "直接证据",
    "indirect": "间接证据",
    "computed": "计算结果",
    "hypothesis": "待验证",
    "missing": "信息缺失",
    "workflow": "流程顺序",
    "concept": "概念顺序",
}
SOLID = {"direct", "computed", "workflow"}


def short_text(value, name, limit):
    if not isinstance(value, str) or not value.strip() or len(value) > limit:
        raise ValueError(f"{name}: expected nonempty text of at most {limit} characters")
    if any(ord(c) < 32 for c in value):
        raise ValueError(f"{name}: control characters are not allowed")


def validate(data):
    if not isinstance(data, dict) or data.get("purpose") not in {"template", "project"}:
        raise ValueError("purpose must be template or project")
    figures = data.get("figures")
    if not isinstance(figures, list) or not 1 <= len(figures) <= 30:
        raise ValueError("figures must contain 1–30 items")
    ids = set()
    for f in figures:
        if not isinstance(f, dict):
            raise ValueError("figure must be an object")
        fid = f.get("id", "")
        if not isinstance(fid, str) or not re.fullmatch(r"[a-z0-9][a-z0-9-]{0,63}", fid):
            raise ValueError("invalid figure id")
        if fid in ids:
            raise ValueError("duplicate figure id")
        ids.add(fid)
        short_text(f.get("title"), "title", 32)
        nodes, links, layout = f.get("nodes"), f.get("links"), f.get("layout")
        if layout not in {"flow", "decision"} or not isinstance(nodes, list):
            raise ValueError("layout must be flow or decision with nodes")
        if not 2 <= len(nodes) <= 4 or (layout == "decision" and len(nodes) != 4):
            raise ValueError("use 2–4 nodes for flow, exactly 4 for decision; split dense figures")
        if not isinstance(links, list) or len(links) != len(nodes) - 1:
            raise ValueError("one link per adjacent flow pair or decision branch is required")
        for node in nodes:
            if not isinstance(node, dict) or node.get("role") not in PALETTE:
                raise ValueError("invalid node role")
            short_text(node.get("label"), "label", 28)
        for link in links:
            if not isinstance(link, dict) or link.get("status") not in STATUS:
                raise ValueError("invalid link evidence status")
            short_text(link.get("basis"), "basis", 2000)
            if link["status"] in {"direct", "indirect", "computed"} and not link.get("source"):
                raise ValueError("evidence-bearing links need source URL or local record identifier")
            if "source" in link:
                short_text(link["source"], "source", 2000)
    return data


def text(x, y, value, size=18, color="#263238", anchor="middle", weight=400):
    return (f'<text x="{x}" y="{y}" fill="{color}" font-size="{size}" '
            f'font-weight="{weight}" text-anchor="{anchor}">{escape(value)}</text>')


def wrap(value, units=22):
    """Conservative width estimate: CJK=2 units, other glyphs=1."""
    lines, line, used = [], "", 0
    for char in value:
        width = 2 if ord(char) > 255 else 1
        if used + width > units:
            lines.append(line)
            line, used = "", 0
        line += char
        used += width
    return lines + ([line] if line else [])


def render_figure(f, purpose):
    decision = f["layout"] == "decision"
    flow_width = 1100 / len(f["nodes"]) - 120
    box_height = max(110, max(len(wrap(n["label"], int((flow_width - 32) / 11)))
                              for n in f["nodes"]) * 26 + 36)
    height = 500 if decision else 124 + box_height + 96
    boxes = ([(455, 92, 290, 84), (55, 278, 290, 112),
              (455, 278, 290, 112), (855, 278, 290, 112)] if decision else
             [(110 + i * (1100 / len(f["nodes"])), 124,
               flow_width, box_height) for i in range(len(f["nodes"]))])
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 {height}" '
        f'width="1200" height="{height}" role="img">',
        f'<title>{escape(f["title"])}</title>',
        '<desc>概念图，非比例、坐标未定。线型与文字共同标识证据；颜色只表示角色。</desc>',
        '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" '
        'markerWidth="7" markerHeight="7" orient="auto">'
        '<path d="M1 1L9 5L1 9" fill="none" stroke="#64748b" stroke-width="1.4"/>'
        '</marker></defs>',
        '<g font-family="Noto Sans CJK SC, PingFang SC, Microsoft YaHei, sans-serif">',
        f'<rect width="1200" height="{height}" fill="white"/>',
        text(36, 39, f["title"], 25, anchor="start", weight=600),
        text(36, 71, "格式示例｜非真实设计" if purpose == "template" else
             "项目概念示意｜非比例、坐标未定", 16, "#64748b", "start"),
    ]
    for i, link in enumerate(f["links"]):
        status = link["status"]
        dash = "" if status in SOLID else ' stroke-dasharray="7 5"'
        x, y, w, h = boxes[0 if decision else i]
        nx, ny, nw, nh = boxes[i + 1]
        if decision:
            sx, sy, ex, ey = x + w / 2, y + h, nx + nw / 2, ny
            path = f"M{sx},{sy}V206H{ex}V{ey-3}"
            tx, ty = ex, 252
        else:
            sx, sy, ex, ey = x + w, y + h / 2, nx, ny + nh / 2
            path = f"M{sx},{sy}H{ex-3}"
            tx, ty = (sx + ex) / 2, sy - 16
        parts.append(f'<path d="{path}" fill="none" stroke="#64748b" stroke-width="1.6"'
                     f'{dash} marker-end="url(#arrow)"/>')
        parts.append(text(tx, ty, STATUS[status], 15, "#475569"))
    for node, (x, y, w, h) in zip(f["nodes"], boxes):
        fill, stroke = PALETTE[node["role"]]
        parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="5" '
                     f'fill="{fill}" stroke="{stroke}" stroke-width="1"/>')
        lines = wrap(node["label"], int((w - 32) / 11))
        for j, line in enumerate(lines):
            parts.append(text(x + w / 2, y + h / 2 + 7 + (j-(len(lines)-1)/2)*26, line, 21))
    parts.extend([
        text(36, height-47, "角色：蓝＝输入/感知　橙＝执行　绿＝输出/检测　灰＝对照　红＝风险/阻断",
             16, "#475569", "start"),
        text(36, height-19, "线型：实线＝直接证据/计算结果/流程顺序；虚线＝间接/待验证/缺失/概念。以箭头文字及图注为准。",
             15, "#64748b", "start"),
        "</g></svg>",
    ])
    return "".join(parts)


def generate(data, out, png=False):
    validate(data)  # Validate the entire manifest before creating output.
    converter = None
    if png:
        try:
            from cairosvg import svg2png
            converter = svg2png
        except (ImportError, OSError) as exc:
            raise RuntimeError("PNG needs CairoSVG and system Cairo; SVG-only export remains available") from exc
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    index = {"purpose": data["purpose"], "figures": []}
    for f in data["figures"]:
        svg = render_figure(f, data["purpose"])
        # Text placeholders and facts are input data, not interpreted as code.
        svg_path = out / (f["id"] + ".svg")
        svg_path.write_text(svg, encoding="utf-8")
        files = [svg_path.name]
        if converter:
            png_path = out / (f["id"] + ".png")
            converter(bytestring=svg.encode("utf-8"), write_to=str(png_path), scale=2)
            files.append(png_path.name)
        index["figures"].append({**f, "files": files})
    (out / "figure-index.json").write_text(
        json.dumps(index, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return index


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--png", action="store_true", help="also export PNG using optional CairoSVG")
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    index = generate(data, args.out, args.png)
    print(json.dumps({"purpose": index["purpose"], "figure_count": len(index["figures"]),
                      "formats": ["svg", "png"] if args.png else ["svg"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
