#!/usr/bin/env python3
"""
MR3 representation-tree gallery.

Runs the instrumented mr3_block on a curated set of easy + medium + hard
test matrices, captures each tree, and renders side-by-side SVG trees in
docs/mr3_tree_gallery.html (self-contained HTML; no JS, no CDN).

Each tree node:
  - circle, area ∝ cluster size, color = log10(NCD_min)
  - label: id, n, NCD_max, NCD_min, outcome
  - hover tooltip (SVG <title>) with full per-node detail
  - edges to children labeled with the shift τ

Usage:
  python3 experiments/mr3_tree_gallery.py
  python3 experiments/mr3_tree_gallery.py --add demmel_S1pe_k8 100
"""
import os
import sys
import io
import json
import math
import argparse
from collections import defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

# Suppress full_eval's import-time test loop
_saved = sys.stdout; sys.stdout = io.StringIO()
try:
    import mr3_gk
    from mr3_gk import (
        bisect_evals, classify, is_zero_width,
        select_shift_geometric, refine_eval, compute_singleton,
        GAPTOL, MAX_DEPTH, SAFMIN,
    )
    from full_eval import make
finally:
    sys.stdout = _saved

OUT_HTML = os.path.join(ROOT, "docs", "mr3_tree_gallery.html")
OUT_JSON = os.path.join(ROOT, "experiments", "mr3_tree_gallery.json")


# ---------- tracing ----------

def ncd_max(lo, hi):
    sp = hi - lo; d = max(abs(lo), abs(hi))
    return sp / d if d > 0 else 0.0


def ncd_min(lo, hi):
    sp = hi - lo; d = max(min(abs(lo), abs(hi)), SAFMIN)
    return sp / d


def trace_mr3_block(c, e_off, nn, pos_evals, spdiam):
    """Same instrumented driver as cluster_tree_trace.py.
    Returns list of node dicts ordered by visit (pre-order, since it's a stack)."""
    nev = len(pos_evals)
    result_evals = np.zeros(nev)
    result_evecs = np.zeros((nn, nev))
    next_id = [0]

    def new_id():
        i = next_id[0]; next_id[0] = i + 1; return i

    root_id = new_id()
    stack = [(c[:nn].copy(), e_off[:max(nn-1, 0)].copy(), nn,
              pos_evals.copy(), list(range(nev)), 0.0, 0, root_id, -1)]
    nodes = []
    while stack:
        bc, be, bnn, local_ev, out_idx, shift, depth, my_id, parent_id = stack.pop()
        if len(local_ev) == 0:
            continue
        lo = float(np.min(local_ev)); hi = float(np.max(local_ev))
        node = {
            "id": my_id, "parent": parent_id, "depth": depth,
            "n_in_node": len(local_ev), "shift": float(shift),
            "lo": lo, "hi": hi,
            "ncd_max": ncd_max(lo, hi),
            "ncd_min": ncd_min(lo, hi),
            "outcome": None,
            "n_singletons": 0, "n_zerowidth_batches": 0, "n_descend_groups": 0,
            "children": [],
            "tau_to_children": [],
            # eigenvalues sample (cap 50 for plot)
            "eigvals_sample": (
                local_ev[:50].tolist()
                if len(local_ev) <= 50
                else local_ev[::max(1, len(local_ev)//50)][:50].tolist()
            ),
        }

        if depth >= MAX_DEPTH:
            node["outcome"] = "max_depth_singletons"
            mult = len(local_ev)
            for j in range(mult):
                lam_ref, z = compute_singleton(
                    bc, be, bnn, local_ev[j], all_evals=local_ev,
                    eval_idx=j, depth=depth, twist_rank=j, multiplicity=mult)
                result_evals[out_idx[j]] = lam_ref + shift
                result_evecs[:bnn, out_idx[j]] = z
            nodes.append(node)
            continue

        groups = classify(local_ev, GAPTOL)
        for gs, ge in groups:
            gev = local_ev[gs:ge+1]; gidx = out_idx[gs:ge+1]
            if gs == ge:
                node["n_singletons"] += 1
                lam_ref, z = compute_singleton(
                    bc, be, bnn, gev[0], all_evals=local_ev,
                    eval_idx=gs, depth=depth)
                result_evals[gidx[0]] = lam_ref + shift
                result_evecs[:bnn, gidx[0]] = z
                continue
            if is_zero_width(gev):
                node["n_zerowidth_batches"] += 1
                mult = len(gev)
                for j in range(mult):
                    lam_ref, z = compute_singleton(
                        bc, be, bnn, gev[j], all_evals=local_ev,
                        eval_idx=gs+j, depth=depth,
                        twist_rank=j, multiplicity=mult)
                    result_evals[gidx[j]] = lam_ref + shift
                    result_evecs[:bnn, gidx[j]] = z
                continue
            node["n_descend_groups"] += 1
            tau, cc, ee = select_shift_geometric(bc, be, bnn, gev, spdiam, local_ev)
            shifted_ev = gev - tau
            refined = np.array([refine_eval(cc, ee, bnn, s) for s in shifted_ev])
            si = np.argsort(refined); refined = refined[si]
            sgidx = [gidx[i] for i in si]
            child_id = new_id()
            node["children"].append(child_id)
            node["tau_to_children"].append(float(tau))
            stack.append((cc, ee, bnn, refined, sgidx,
                          shift + tau, depth + 1, child_id, my_id))

        node["outcome"] = "internal" if node["n_descend_groups"] > 0 else "leaf"
        nodes.append(node)
    return nodes


def build_tgk_offdiag(d, e):
    n = len(d); m = 2 * n
    c = np.zeros(m, dtype=np.float64)
    e_off = np.zeros(m - 1, dtype=np.float64)
    for i in range(n):
        e_off[2*i] = abs(float(d[i]))
    for i in range(n - 1):
        e_off[2*i+1] = abs(float(e[i]))
    spdiam = float(np.max(np.abs(e_off))) if m > 1 else 0.0
    return c, e_off, n, spdiam


def trace_matrix(name, n):
    d, e = make(name, n)
    c, e_off, n_nonneg, spdiam = build_tgk_offdiag(d, e)
    pos_evals = bisect_evals(c, e_off, len(c), n + 1, 2 * n)
    if len(pos_evals) != n:
        return None
    return trace_mr3_block(c, e_off, len(c), pos_evals, spdiam)


# ---------- SVG rendering ----------

def render_tree_svg(label, nodes, max_width_px=560, node_min_r=10, node_max_r=42):
    """Render a single tree to an SVG string."""
    if not nodes:
        return f'<div class="empty">no nodes for {label}</div>'

    # Build child lookup by id
    by_id = {n["id"]: n for n in nodes}
    children_of = defaultdict(list)
    for n in nodes:
        if n["parent"] >= 0:
            children_of[n["parent"]].append(n["id"])

    # Compute node sizes by cluster size
    sizes = [n["n_in_node"] for n in nodes]
    max_n = max(sizes)
    def radius(n_in):
        if max_n <= 1:
            return node_min_r
        return node_min_r + (node_max_r - node_min_r) * math.sqrt(n_in / max_n)

    # Color: log10(NCD_min) → viridis-like 0..1
    ncd_logs = [math.log10(max(n["ncd_min"], 1e-300)) for n in nodes]
    cmin = min(ncd_logs); cmax = max(ncd_logs)
    span = max(cmax - cmin, 1e-9)
    def color(ncd_min):
        v = (math.log10(max(ncd_min, 1e-300)) - cmin) / span
        # interpolate dark-blue (0) → cyan (0.5) → yellow (1.0)
        if v < 0.5:
            t = v * 2.0
            r = int(40 + (60 - 40) * t)
            g = int(40 + (180 - 40) * t)
            b = int(140 + (220 - 140) * t)
        else:
            t = (v - 0.5) * 2.0
            r = int(60 + (250 - 60) * t)
            g = int(180 + (220 - 180) * t)
            b = int(220 + (60 - 220) * t)
        return f"rgb({r},{g},{b})"

    # Layout: walk tree in DFS, assign X positions based on subtree leaf counts.
    # We compute leaf-counts (descendant leaves), then X is based on leaf order.
    leaf_counts = {}
    def count_leaves(nid):
        ch = children_of.get(nid, [])
        if not ch:
            leaf_counts[nid] = 1
            return 1
        total = 0
        for c in ch:
            total += count_leaves(c)
        leaf_counts[nid] = total
        return total

    # Find all root ids (parent == -1)
    roots = [n["id"] for n in nodes if n["parent"] == -1]
    for r in roots:
        count_leaves(r)

    # Assign X by inorder leaf positions
    pos = {}
    counter = [0]
    def assign_x(nid):
        ch = children_of.get(nid, [])
        if not ch:
            x = counter[0]; counter[0] += 1
            pos[nid] = x
            return x
        xs = [assign_x(c) for c in ch]
        x = (xs[0] + xs[-1]) / 2.0
        pos[nid] = x
        return x
    for r in roots:
        assign_x(r)

    total_leaves = max(counter[0], 1)
    max_depth = max(n["depth"] for n in nodes)

    # Coordinate system
    margin_x = 40
    margin_y = 40
    inner_w = max_width_px - 2 * margin_x
    level_h = 110
    height = margin_y * 2 + (max_depth + 1) * level_h
    width = max_width_px

    def tx(x):
        if total_leaves <= 1:
            return margin_x + inner_w / 2
        return margin_x + (x / max(total_leaves - 1, 1)) * inner_w
    def ty(d):
        return margin_y + d * level_h

    # SVG body
    svg = [f'<svg class="tree-svg" viewBox="0 0 {width} {height}" '
           f'xmlns="http://www.w3.org/2000/svg">']
    # title bar
    svg.append(f'<text x="{margin_x}" y="20" class="tree-title">'
               f'{label}: {len(nodes)} nodes, max depth {max_depth}</text>')

    # edges
    for n in nodes:
        for ci, cid in enumerate(n["children"]):
            x0, y0 = tx(pos[n["id"]]), ty(n["depth"])
            x1, y1 = tx(pos[cid]), ty(n["depth"] + 1)
            # gentle curve
            cy = y0 + (y1 - y0) * 0.5
            svg.append(
                f'<path d="M{x0:.1f},{y0:.1f} C{x0:.1f},{cy:.1f} {x1:.1f},{cy:.1f} {x1:.1f},{y1:.1f}" '
                f'class="edge" />')
            tau = n["tau_to_children"][ci] if ci < len(n["tau_to_children"]) else None
            if tau is not None:
                svg.append(
                    f'<text x="{(x0+x1)/2:.1f}" y="{cy:.1f}" class="edge-label">'
                    f'τ={tau:.2g}</text>')

    # nodes
    for n in nodes:
        x, y = tx(pos[n["id"]]), ty(n["depth"])
        r = radius(n["n_in_node"])
        col = color(n["ncd_min"])
        title = (
            f'id={n["id"]}  parent={n["parent"]}\n'
            f'depth={n["depth"]}  n_in_node={n["n_in_node"]}\n'
            f'lo={n["lo"]:.4g}  hi={n["hi"]:.4g}\n'
            f'NCD_max={n["ncd_max"]:.3g}\n'
            f'NCD_min={n["ncd_min"]:.3g}\n'
            f'shift={n["shift"]:.3g}\n'
            f'outcome={n["outcome"]}\n'
            f'singletons/zerowidth/descend = '
            f'{n["n_singletons"]}/{n["n_zerowidth_batches"]}/{n["n_descend_groups"]}'
        )
        svg.append(
            f'<g class="node">'
            f'<title>{_xml_escape(title)}</title>'
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" '
            f'fill="{col}" stroke="#333" stroke-width="1" />'
            f'<text x="{x:.1f}" y="{y+4:.1f}" class="node-label">'
            f'{n["n_in_node"]}</text>'
            f'</g>'
        )
    svg.append('</svg>')
    return "\n".join(svg)


def _xml_escape(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace("\"", "&quot;").replace("'", "&apos;"))


# ---------- HTML composition ----------

CSS = """
body { font-family: -apple-system, "Segoe UI", Helvetica, sans-serif;
       max-width: 1280px; margin: 1em auto; color: #222; line-height: 1.5;
       padding: 0 1em 4em; }
h1 { color: #06335c; border-left: 4px solid #06c; padding-left: 0.5em; }
h2 { color: #06335c; margin-top: 2em; border-bottom: 1px solid #ccc; padding-bottom: 0.2em; }
h3 { color: #244; margin-top: 1.5em; }
.controls { background: #f6f9fc; padding: 0.6em 0.9em; border-radius: 4px;
            font-size: 0.9em; margin-bottom: 1em; }
.tree-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 1.2em; }
@media (max-width: 880px) { .tree-grid { grid-template-columns: 1fr; } }
.tree-card { background: #fff; border: 1px solid #ddd; border-radius: 5px;
             padding: 0.7em 0.9em; }
.tree-card h3 { margin: 0 0 0.4em 0; color: #024; font-size: 1em;
                font-family: SF Mono, Menlo, Consolas, monospace; }
.tree-card .stats { font-size: 0.85em; color: #444; margin-bottom: 0.4em;
                    font-family: SF Mono, Menlo, Consolas, monospace; }
.tree-svg { width: 100%; height: auto; background: #fbfbfd; }
text.tree-title { font-family: SF Mono, Menlo, Consolas, monospace; font-size: 11px; fill: #555; }
text.node-label { font-family: SF Mono, Menlo, Consolas, monospace; font-size: 11px;
                  fill: #fff; text-anchor: middle; pointer-events: none;
                  font-weight: bold;
                  paint-order: stroke; stroke: #00000060; stroke-width: 2px; }
text.edge-label { font-family: SF Mono, Menlo, Consolas, monospace; font-size: 9px;
                  fill: #777; text-anchor: middle; pointer-events: none; }
.edge { fill: none; stroke: #888; stroke-width: 1.2; }
.legend { display: flex; gap: 1em; align-items: center; font-size: 0.88em;
          margin-bottom: 1em; flex-wrap: wrap; }
.legend .swatch { display: inline-block; width: 14px; height: 14px;
                  border-radius: 50%; vertical-align: middle; margin-right: 4px;
                  border: 1px solid #444; }
.detail-box { font-size: 0.78em; max-height: 130px; overflow-y: auto;
              background: #f6f8fa; padding: 0.4em 0.6em; border-radius: 3px;
              font-family: SF Mono, Menlo, Consolas, monospace; }
.banner { background: #fef9e7; border-left: 4px solid #c80; padding: 0.6em 0.9em;
          border-radius: 3px; margin-top: 1em; }
"""


def card_html(label, nodes):
    if not nodes:
        return f'<div class="tree-card"><h3>{label}</h3><div class="empty">no tree</div></div>'
    max_depth = max(n["depth"] for n in nodes)
    n_nodes = len(nodes)
    n_leaves = sum(1 for n in nodes if not n["children"])
    root = nodes[0]
    stats = (f'nodes={n_nodes} &nbsp; depth={max_depth} &nbsp; leaves={n_leaves} &nbsp;'
             f' root NCD_max={root["ncd_max"]:.3g} &nbsp; '
             f'root NCD_min={root["ncd_min"]:.3g}')
    detail = []
    for n in nodes[:25]:
        detail.append(
            f'id={n["id"]:>3}  d={n["depth"]}  n={n["n_in_node"]:>4}  '
            f'NCD_max={n["ncd_max"]:.2e}  NCD_min={n["ncd_min"]:.2e}  '
            f'sing={n["n_singletons"]:<3}  zw={n["n_zerowidth_batches"]:<2}  '
            f'desc={n["n_descend_groups"]:<2}  out={n["outcome"]}'
        )
    if len(nodes) > 25:
        detail.append(f'... ({len(nodes) - 25} more)')
    detail_html = "<pre class='detail-box'>" + "\n".join(detail) + "</pre>"
    svg = render_tree_svg(label, nodes)
    return (f'<div class="tree-card"><h3>{label}</h3>'
            f'<div class="stats">{stats}</div>{svg}{detail_html}</div>')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--add", nargs=2, action="append", default=[],
                    metavar=("NAME", "N"),
                    help="add an extra (name, n) case (can repeat)")
    args = ap.parse_args()

    cases_easy = [
        ("diagonal_only",          50),
        ("constant",               50),
        ("arithmetic_progression", 50),
        ("random_uniform",         50),
        ("exponential_graded",     50),
        ("alternating_sign",       50),
    ]
    cases_medium = [
        ("step_function",          200),
        ("two_clusters",           400),
        ("three_clusters",         400),
        ("demmel_S1pe_k8",         400),
        ("demmel_S1pe_k4",         400),
        ("saw_tooth",              400),
    ]
    cases_hard = [
        ("one_big_cluster",        400),
        ("one_big_cluster",        200),
        ("random_clustered_5",     400),
        ("constant_d_graded_e",    400),
        ("gl_clustered_at_eps",    400),
        ("gl_clustered_at_pmeps",  400),
    ]
    cases_wilk = [
        ("wilkinson_exact",        200),    # canonical Wilkinson W21+
        ("wilkinson_like",         200),    # synthetic Wilkinson-like
        ("glued_wilkinson",        200),    # glued copies (paper testset)
        ("glued_wilkinson_tight",  200),    # glued with sqrt(eps) gluing
        ("gl_wilkm",               400),    # Großer-Lang Wilkinson minus
        ("gl_wilkw",               400),    # Großer-Lang Wilkinson wide
        ("gl_wilk2w",              200),    # double-wide; root-only (good contrast)
    ]
    cases_stcoll = [
        ("ST", "B_Kimura_429"),       # tag, basename
        ("ST", "B_gg_30_1D-5"),
        ("ST", "B_glued_09b"),
        ("ST", "B_20_graded"),
        ("ST", "B_40_graded"),
        ("ST", "B_bug316_gesdd"),
    ]
    extra = [(name, int(n)) for (name, n) in args.add]

    sections = [
        ("Easy: well-separated eigenvalues — single-root trees", cases_easy),
        ("Medium: clustered subspaces — depth-1 / 2 trees", cases_medium),
        ("Hard: deep clusters — depth-2 / 3 trees", cases_hard),
        ("Wilkinson family", cases_wilk),
        ("STCollection (real-world)", cases_stcoll),
    ]
    if extra:
        sections.append(("User-added cases", extra))

    all_traces = {}
    html_sections = []

    # Lazy STCollection loader
    _stcoll_de_cache = {}
    def _trace_stcoll(basename):
        import glob
        if basename in _stcoll_de_cache:
            d, e = _stcoll_de_cache[basename]
        else:
            paths = glob.glob(os.path.join(ROOT, "stcollection",
                                           f"{basename}.dat"))
            if not paths:
                return None, None
            saved = sys.stdout; sys.stdout = io.StringIO()
            try:
                from evaluate import load_stcoll
            finally:
                sys.stdout = saved
            _, d, e = load_stcoll(paths[0])
            _stcoll_de_cache[basename] = (d, e)
        n = len(d)
        c, e_off, _, sp = build_tgk_offdiag(d, e)
        pos = bisect_evals(c, e_off, len(c), n + 1, 2 * n)
        if len(pos) != n:
            return None, n
        return trace_mr3_block(c, e_off, len(c), pos, sp), n

    for sec_title, cases in sections:
        cards = []
        for entry in cases:
            if entry[0] == "ST":
                # STCollection lookup
                _, basename = entry
                try:
                    t, n = _trace_stcoll(basename)
                except Exception as ex:
                    print(f"  ST_{basename}: ERROR {ex}")
                    continue
                if t is None:
                    print(f"  ST_{basename}: skipped (n={n})")
                    continue
                label = f"ST_{basename}@{n}"
            else:
                name, n = entry
                label = f"{name}@{n}"
                try:
                    t = trace_matrix(name, n)
                except Exception as ex:
                    print(f"  {label}: ERROR {ex}")
                    continue
                if t is None:
                    print(f"  {label}: bisect mismatch — skipped")
                    continue
            print(f"  {label}: {len(t)} nodes, depth {max(nd['depth'] for nd in t)}")
            all_traces[label] = t
            cards.append(card_html(label, t))
        html_sections.append((sec_title, cards))

    body_parts = []
    body_parts.append('<h1>MR3 Representation-Tree Gallery</h1>')
    body_parts.append('<div class="controls">'
        'Trees produced by the instrumented Python <code>mr3_block</code> '
        '(<code>mr3_gk.py:418-460</code>), which is algorithmically equivalent '
        'to the Fortran kernel (verified by bit-identity on 21,375 specs). '
        'Each circle is a representation-tree node. <b>Area</b> ∝ cluster '
        'size; <b>color</b> = log<sub>10</sub>(NCD_min). Hover over a node '
        'to see full per-node detail (lo, hi, NCD_max, NCD_min, shift, outcome).</div>')
    body_parts.append('<div class="legend">'
        '<span><span class="swatch" style="background:rgb(40,40,140)"></span>'
        'low NCD_min (tight cluster)</span>'
        '<span><span class="swatch" style="background:rgb(60,180,220)"></span>'
        'medium</span>'
        '<span><span class="swatch" style="background:rgb(250,220,60)"></span>'
        'high NCD_min (wide dynamic range)</span>'
        '</div>')

    for sec_title, cards in html_sections:
        body_parts.append(f'<h2>{sec_title}</h2>')
        body_parts.append('<div class="tree-grid">')
        body_parts.extend(cards)
        body_parts.append('</div>')

    body_parts.append('<div class="banner"><b>Note on tree depth.</b> Most '
        'matrices in our test corpus descend at most one level: at the root, '
        '<code>classify</code> finds enough relative gap (NCD ≥ GAPTOL = 1e-3) '
        'between adjacent eigenvalues to dispatch most evals as singletons, '
        'and the remaining sub-clusters resolve via <code>compute_singleton</code> '
        'or <code>is_zero_width</code> at depth 1. Deeper-tree behaviour shows '
        'on the harder STCollection matrices (<code>T_Godunov_*</code>, '
        '<code>ev3_ec1_*</code>) — out of this gallery&#39;s scope but '
        'auditable via <code>experiments/cluster_tree_trace.py</code> on the '
        'STCollection bidiagonals.</div>')

    body = "\n".join(body_parts)
    html = (f"<!doctype html><html><head><meta charset='utf-8'>"
            f"<title>MR3 Tree Gallery</title>"
            f"<style>{CSS}</style></head>"
            f"<body>{body}</body></html>")

    os.makedirs(os.path.dirname(OUT_HTML), exist_ok=True)
    with open(OUT_HTML, "w") as f:
        f.write(html)
    with open(OUT_JSON, "w") as f:
        json.dump(all_traces, f, indent=1, default=float)
    print(f"\nwrote {OUT_HTML}  ({os.path.getsize(OUT_HTML)//1024} KB, "
          f"{len(all_traces)} trees)")
    print(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
