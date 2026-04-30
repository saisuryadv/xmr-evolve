#!/usr/bin/env python3
"""
Render docs/mr3_tree_gallery_fortran.html from
experiments/mr3_tree_fortran_log.json.

Each matrix's tree appears TWICE side-by-side:
  - left  SVG:  nodes colored by NCD_eigs   = (eig_hi - eig_lo) / max(|eig|)
  - right SVG:  nodes colored by NCD_pivots = (|D|max - |D|min) / |D|max

All numbers come from the actual Fortran kernel — `dlaxrf_traced.f` writes
one log line per *successful child representation* it builds. The "root"
representation built by `dlaxre` is NOT logged here; the tree shown is the
forest of child reps, ordered by visit order, with depth as written.

Tooltips (SVG <title>) show every node's full Fortran-side state.
"""
import os
import sys
import json
import math
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN_JSON = os.path.join(ROOT, "experiments", "mr3_tree_fortran_log.json")
OUT_HTML = os.path.join(ROOT, "docs", "mr3_tree_gallery_fortran.html")


CSS = """
body { font-family: -apple-system, "Segoe UI", Helvetica, sans-serif;
       max-width: 1320px; margin: 1em auto; color: #222; line-height: 1.5;
       padding: 0 1em 4em; }
h1 { color: #06335c; border-left: 4px solid #06c; padding-left: 0.5em; }
h2 { color: #06335c; margin-top: 1.6em; border-bottom: 1px solid #ccc;
     padding-bottom: 0.2em; font-size: 1.2em; }
.controls { background: #f6f9fc; padding: 0.6em 0.9em; border-radius: 4px;
            font-size: 0.92em; margin-bottom: 1em; }
.controls .formula { background: #eef; padding: 0.4em 0.6em; border-radius: 3px;
                     font-family: SF Mono, Menlo, Consolas, monospace;
                     margin-top: 0.4em; font-size: 0.88em; }
.tree-card { background: #fff; border: 1px solid #ddd; border-radius: 5px;
             padding: 0.5em 0.7em; margin-bottom: 1em; }
.tree-card h3 { margin: 0 0 0.3em; font-family: SF Mono, Menlo, Consolas, monospace;
                font-size: 0.93em; color: #024; }
.tree-row { display: grid; grid-template-columns: 1fr 1fr; gap: 0.5em; }
.tree-svg { width: 100%; height: auto; background: #fbfbfd; }
text.tree-title { font-family: SF Mono, Menlo, Consolas, monospace;
                  font-size: 10px; fill: #555; }
text.node-label { font-family: SF Mono, Menlo, Consolas, monospace;
                  font-size: 10px; fill: #fff; text-anchor: middle;
                  pointer-events: none; font-weight: bold;
                  paint-order: stroke; stroke: #00000060; stroke-width: 2px; }
text.ncd-label  { font-family: SF Mono, Menlo, Consolas, monospace;
                  font-size: 9px; fill: #222; text-anchor: middle;
                  pointer-events: none;
                  paint-order: stroke; stroke: #ffffffaa; stroke-width: 2px; }
text.edge-label { font-family: SF Mono, Menlo, Consolas, monospace;
                  font-size: 8.5px; fill: #777; text-anchor: middle;
                  pointer-events: none; }
.edge { fill: none; stroke: #888; stroke-width: 1.2; }
.detail-box { font-size: 0.74em; max-height: 130px; overflow-y: auto;
              background: #f6f8fa; padding: 0.4em 0.6em; border-radius: 3px;
              font-family: SF Mono, Menlo, Consolas, monospace;
              white-space: pre; margin-top: 0.4em; }
.banner { background: #fef9e7; border-left: 4px solid #c80; padding: 0.6em 0.9em;
          border-radius: 3px; margin-top: 1em; font-size: 0.92em; }
.legend { display: flex; gap: 1em; align-items: center; font-size: 0.86em;
          margin-bottom: 0.8em; flex-wrap: wrap; }
.legend .swatch { display: inline-block; width: 14px; height: 14px;
                  border-radius: 50%; vertical-align: middle; margin-right: 4px;
                  border: 1px solid #444; }
"""


def color_for(v, log_scale, vmin, vmax):
    """Return an rgb(...) string. v in [vmin..vmax]."""
    if log_scale:
        v = max(v, 1e-300); vmin = max(vmin, 1e-300); vmax = max(vmax, 1e-300)
        t = (math.log10(v) - math.log10(vmin)) / max(
            math.log10(vmax) - math.log10(vmin), 1e-9)
    else:
        t = (v - vmin) / max(vmax - vmin, 1e-9)
    t = max(0.0, min(1.0, t))
    # blue (0) -> green (0.5) -> red (1)
    r = int(40 + (250 - 40) * t)
    g = int(40 + (220 - 40) * (1 - abs(t - 0.5) * 2))
    b = int(220 - (220 - 60) * t)
    return f"rgb({r},{g},{b})"


def _xml_escape(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace("\"", "&quot;").replace("'", "&apos;"))


def render_svg(label, nodes, color_field, max_width_px=620,
               node_min_r=8, node_max_r=34):
    """Render the Fortran tree.

    The Fortran trace gives us a flat list of child nodes ordered by
    visit; we infer parent links by depth (each node's parent is the
    most-recently-pushed node at depth-1 with a not-yet-resolved icbeg/icend
    that contains the current node's range). For our visualization purposes,
    we use a simpler heuristic: each node at depth d > 0 is a child of the
    most-recent node at depth d-1 in scan order.
    """
    if not nodes:
        return f'<div class="empty">no descent (single root, no children).</div>'

    # Assign IDs in scan order
    for i, n in enumerate(nodes):
        n["id"] = i

    # Infer parent: most recent node at (depth-1) before this one.
    # Depth-0 nodes have no parent.
    last_at_depth = {}
    for n in nodes:
        d = n["depth"]
        if d == 0:
            n["parent"] = -1
        else:
            n["parent"] = last_at_depth.get(d - 1, -1)
        last_at_depth[d] = n["id"]

    # Build child lookup
    children_of = defaultdict(list)
    for n in nodes:
        if n["parent"] >= 0:
            children_of[n["parent"]].append(n["id"])

    # Color range
    vals = [max(n[color_field], 1e-300) for n in nodes]
    vmin, vmax = min(vals), max(vals)

    # Sizes by cluster size
    sizes = [n["n_in_node"] for n in nodes]
    max_n = max(sizes)
    def radius(n_in):
        if max_n <= 1:
            return node_min_r
        return node_min_r + (node_max_r - node_min_r) * math.sqrt(n_in / max_n)

    # X positions: depth-0 nodes are spread across; deeper nodes nest under
    # their parent.
    leaf_count = {}
    def count_leaves(nid):
        ch = children_of.get(nid, [])
        if not ch:
            leaf_count[nid] = 1; return 1
        total = sum(count_leaves(c) for c in ch)
        leaf_count[nid] = total
        return total

    # Synthetic forest: depth-0 nodes are roots
    roots = [n["id"] for n in nodes if n["depth"] == 0]
    if not roots:
        # all depth >=1 → single virtual root
        roots = []
    for r in roots:
        count_leaves(r)

    pos = {}; counter = [0]
    def assign_x(nid):
        ch = children_of.get(nid, [])
        if not ch:
            x = counter[0]; counter[0] += 1; pos[nid] = x; return x
        xs = [assign_x(c) for c in ch]
        x = (xs[0] + xs[-1]) / 2.0
        pos[nid] = x; return x
    for r in roots:
        assign_x(r)

    total_leaves = max(counter[0], 1)
    max_depth = max(n["depth"] for n in nodes)

    margin_x, margin_y = 36, 32
    inner_w = max_width_px - 2 * margin_x
    level_h = 110  # extra room for NCD label under each node
    height = margin_y * 2 + (max_depth + 1) * level_h
    width = max_width_px

    def tx(x):
        if total_leaves <= 1:
            return margin_x + inner_w / 2
        return margin_x + (x / max(total_leaves - 1, 1)) * inner_w
    def ty(d):
        return margin_y + d * level_h

    # ncd_recon spans many orders of magnitude → log color scale
    log_scale = (color_field == "ncd_recon")
    svg = [f'<svg class="tree-svg" viewBox="0 0 {width} {height}" '
           f'xmlns="http://www.w3.org/2000/svg">',
           f'<text x="{margin_x}" y="18" class="tree-title">'
           f'{label} ({len(nodes)} nodes, depth 0..{max_depth}, '
           f'color={color_field})</text>']

    # edges
    for n in nodes:
        for cid in children_of[n["id"]]:
            x0, y0 = tx(pos[n["id"]]), ty(n["depth"])
            x1, y1 = tx(pos[cid]), ty(n["depth"] + 1)
            cy = y0 + (y1 - y0) * 0.5
            svg.append(
                f'<path d="M{x0:.1f},{y0:.1f} C{x0:.1f},{cy:.1f} '
                f'{x1:.1f},{cy:.1f} {x1:.1f},{y1:.1f}" class="edge" />')
            child = next((nn for nn in nodes if nn["id"] == cid), None)
            if child is not None:
                svg.append(
                    f'<text x="{(x0+x1)/2:.1f}" y="{cy:.1f}" '
                    f'class="edge-label">τ={child["tau"]:.2g}</text>')

    # nodes
    for n in nodes:
        x, y = tx(pos[n["id"]]), ty(n["depth"])
        r = radius(n["n_in_node"])
        col = color_for(n[color_field], log_scale, vmin, vmax)
        title = (
            f'id={n["id"]}  depth={n["depth"]}  '
            f'icbeg..icend={n["icbeg"]}..{n["icend"]}  '
            f'(n_in_node={n["n_in_node"]} of n_block={n["n"]})\n'
            f'tau={n["tau"]:.4g}\n'
            f'reconstructed diag(LDL^T) range:\n'
            f'  recon_dmin = {n["recon_dmin"]:+.4g}\n'
            f'  recon_dmax = {n["recon_dmax"]:+.4g}\n'
            f'NCD_recon = (recon_dmax - recon_dmin) / |recon_dmin|\n'
            f'          = {n["ncd_recon"]:.4g}'
        )
        # Format the NCD compactly (e.g. 1.0e+11)
        ncd_v = n["ncd_recon"]
        if ncd_v == 0.0:
            ncd_str = "0"
        elif abs(ncd_v) >= 1e3 or abs(ncd_v) < 1e-2:
            ncd_str = f"{ncd_v:.1e}".replace("e+0", "e").replace("e-0", "e-")
        else:
            ncd_str = f"{ncd_v:.2g}"
        svg.append(
            f'<g class="node">'
            f'<title>{_xml_escape(title)}</title>'
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" '
            f'fill="{col}" stroke="#333" stroke-width="1" />'
            f'<text x="{x:.1f}" y="{y+3.5:.1f}" class="node-label">'
            f'{n["n_in_node"]}</text>'
            f'<text x="{x:.1f}" y="{y+r+11:.1f}" class="ncd-label">'
            f'{ncd_str}</text>'
            f'</g>'
        )
    svg.append('</svg>')
    return "\n".join(svg)


def card_html(label, info):
    nodes = info["nodes"]
    if not nodes:
        return (f'<div class="tree-card"><h3>{label}</h3>'
                f'<div style="color:#888; font-size:0.86em;">'
                f'No descent occurred — every eigenvalue resolved as a '
                f'singleton or zero-width batch at the root level. '
                f'(Tree depth 0; <code>dlaxrf</code> not invoked or never '
                f'produced a successful child.)</div></div>')
    md = info["max_depth"]
    ncd_vals = sorted(n["ncd_recon"] for n in nodes)
    ncd_min = ncd_vals[0]
    ncd_max = ncd_vals[-1]
    ncd_med = ncd_vals[len(ncd_vals) // 2]
    head = (f'{label} — {len(nodes)} child-nodes, max depth {md} '
            f'(matrix n={info["n_in_matrix"]}) &nbsp;|&nbsp; '
            f'NCD_recon min={ncd_min:.3e} '
            f'median={ncd_med:.3e} '
            f'max={ncd_max:.3e}')
    detail = []
    for n in nodes[:30]:
        detail.append(
            f'd={n["depth"]}  ic=[{n["icbeg"]:>3},{n["icend"]:>3}]  '
            f'n={n["n_in_node"]:>4}  tau={n["tau"]:+.3e}  '
            f'NCD_recon={n["ncd_recon"]:.3e}  '
            f'recon=[{n["recon_dmin"]:+.2e},{n["recon_dmax"]:+.2e}]'
        )
    if len(nodes) > 30:
        detail.append(f'… ({len(nodes) - 30} more)')
    detail_block = "<pre class='detail-box'>" + "\n".join(detail) + "</pre>"
    svg_recon = render_svg(label, nodes, "ncd_recon")
    return (f'<div class="tree-card"><h3>{head}</h3>'
            f'<div class="tree-row" style="grid-template-columns:1fr;">{svg_recon}</div>'
            f'{detail_block}</div>')


def main():
    if not os.path.exists(IN_JSON):
        print(f"ERROR: {IN_JSON} not found. Run run_fortran_traced.py first.")
        sys.exit(1)
    with open(IN_JSON) as f:
        data = json.load(f)

    # Sort matrices into sections for readability
    easy = ["diagonal_only@50", "constant@50", "arithmetic_progression@50"]
    medium = ["step_function@200", "two_clusters@400", "three_clusters@400",
              "demmel_S1pe_k4@400", "demmel_S1pe_k8@400", "saw_tooth@400"]
    hard = ["one_big_cluster@400", "random_clustered_5@400",
            "constant_d_graded_e@400", "gl_clustered_at_eps@400", "pd_T0@100"]
    wilk = ["wilkinson_exact@200", "wilkinson_like@200", "glued_wilkinson@200",
            "glued_wilkinson_tight@200", "gl_wilkw@400", "gl_wilkm@400"]
    stcoll = [k for k in data.keys() if k.startswith("ST_")]

    sections = [
        ("Easy: well-separated eigenvalues", easy),
        ("Medium: clustered subspaces", medium),
        ("Hard: deep clusters / wide dynamic range", hard),
        ("Wilkinson family", wilk),
        ("STCollection (real-world)", stcoll),
    ]

    body = ['<h1>MR3 Representation Tree — captured from the Fortran kernel</h1>',
        '<div class="controls">'
        '<p>Each tree is the <b>actual</b> output of '
        '<code>mr3gk_fortran/mr3gk_run_traced</code>, which links a copy of '
        '<code>dlaxrf.f</code> instrumented to write one log line per '
        'successful child representation. <b>No Python re-implementation '
        'is involved</b> in producing these trees — every value below is '
        'computed inside the Fortran kernel (reconstruction in '
        'quad precision via libquadmath) and emitted via '
        '<code>WRITE(99, ...)</code>. The Python harness only parses the '
        'resulting log file.</p>'
        '<div class="formula"><b>recon_d<sub>i</sub></b> = D[i] + e[i−1]<sup>2</sup> / D[i−1]'
        '&nbsp;&nbsp;(reconstructs the diagonal of L·D·L<sup>T</sup>)<br>'
        '<b>NCD_recon</b> = (recon_d<sub>max</sub> − recon_d<sub>min</sub>) / |recon_d<sub>min</sub>|<br>'
        '&nbsp;&nbsp;denominator gets +1e-16 only when |recon_d<sub>min</sub>| ≤ 1e-16</div>'
        '<p style="margin-top:0.5em;">Hover over any node circle to see '
        'recon range and NCD_recon. Color reflects log<sub>10</sub>(NCD_recon).</p>'
        '</div>',
        '<div class="legend">'
        '<span><span class="swatch" style="background:rgb(40,40,220)"></span>low NCD</span>'
        '<span><span class="swatch" style="background:rgb(150,220,150)"></span>medium</span>'
        '<span><span class="swatch" style="background:rgb(250,40,60)"></span>high NCD (≈1)</span>'
        '</div>']

    for sec_title, labels in sections:
        body.append(f'<h2>{sec_title}</h2>')
        for label in labels:
            if label not in data:
                continue
            body.append(card_html(label, data[label]))

    body.append('<div class="banner">'
        '<b>Note on root-level representation.</b> <code>dlaxrf</code> '
        'is invoked only when descending into a sub-cluster — the root '
        'representation is built by <code>dlaxre_gk</code> (not '
        '<code>dlaxrf</code>), so trees with <i>0 child-nodes</i> '
        '(diagonal_only, constant, arithmetic_progression) had every '
        'eigenvalue resolved as a singleton without descent and therefore '
        'show no <code>dlaxrf</code> traces. NCD measurement of the root '
        'rep itself would require instrumenting <code>dlaxre_gk.f</code> '
        'similarly — out of scope for this gallery.</div>')

    html = (f"<!doctype html><html><head><meta charset='utf-8'>"
            f"<title>MR3 Tree Gallery (Fortran-traced)</title>"
            f"<style>{CSS}</style></head>"
            f"<body>{chr(10).join(body)}</body></html>")
    os.makedirs(os.path.dirname(OUT_HTML), exist_ok=True)
    with open(OUT_HTML, "w") as f:
        f.write(html)
    sz = os.path.getsize(OUT_HTML) // 1024
    print(f"wrote {OUT_HTML}  ({sz} KB, {len(data)} trees)")


if __name__ == "__main__":
    main()
