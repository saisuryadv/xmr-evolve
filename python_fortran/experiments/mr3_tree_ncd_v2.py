#!/usr/bin/env python3
"""
MR3 representation-tree gallery v2 — every node logs TWO NCD metrics:

  NCD_eigs := (max_eig - min_eig) / max(|max_eig|, |min_eig|)
              ("(max-min)/max" of eigenvalue range at the node — what XMR uses)

  NCD_pivots := (max|D_i| - min|D_i|) / max(|D_i|)
                from the LDL^T factorization at the node, where
                T_node = (T_GK - shift*I) at the node's depth, and
                T_node = L · D · L^T (top-down LDL).
                This is the "(max-min)/max" of the diagonal of D after the
                LDL^T factorization, per the user's spec.

Renders side-by-side SVG trees in
  docs/mr3_tree_gallery_ncd.html
with hover tooltips showing both NCDs and the pivot vector summary.
"""
import os
import sys
import io
import json
import math
from collections import defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

_saved = sys.stdout; sys.stdout = io.StringIO()
try:
    import mr3_gk
    from mr3_gk import (
        bisect_evals, classify, is_zero_width,
        select_shift_geometric, refine_eval, compute_singleton,
        ldl_top_down,
        GAPTOL, MAX_DEPTH, SAFMIN,
    )
    from full_eval import make
    from evaluate import load_stcoll
finally:
    sys.stdout = _saved

OUT_HTML = os.path.join(ROOT, "docs", "mr3_tree_gallery_ncd.html")
OUT_JSON = os.path.join(ROOT, "experiments", "mr3_tree_gallery_ncd.json")


# ---------- NCD metrics ----------

def ncd_max_eigs(lo, hi):
    """(hi-lo) / max(|hi|, |lo|) — eigenvalue-range NCD, MR3 paper-style."""
    sp = hi - lo
    d = max(abs(lo), abs(hi))
    return sp / d if d > 0 else 0.0


def ncd_pivots(c, e_off, n):
    """Run top-down LDL^T factorization on (c, e_off) and return
       (max|D_i| - min|D_i|) / max(|D_i|).
       This is the user-requested NCD: "(diag_max - diag_min) / diag_max
       of the matrix after multiplying LDL^T of T_GK"."""
    if n <= 0:
        return 0.0, [], 0.0, 0.0
    d_pivots, _ = ldl_top_down(c, e_off, n)
    abs_d = np.abs(d_pivots[:n])
    if len(abs_d) == 0:
        return 0.0, [], 0.0, 0.0
    dmax = float(np.max(abs_d))
    dmin = float(np.min(abs_d))
    if dmax <= 0:
        return 0.0, d_pivots[:n].tolist(), dmin, dmax
    ncd = (dmax - dmin) / dmax
    return ncd, d_pivots[:n].tolist(), dmin, dmax


# ---------- tracing ----------

def trace_mr3_block(c, e_off, nn, pos_evals, spdiam):
    """Same instrumented driver, now also logging NCD_pivots at every node."""
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

        # Pivot NCD: T_node = (T_GK - shift*I) → LDL^T at this node.
        # bc, be already represent the shifted block (the recursion shifts
        # them by tau before passing the child).
        ncd_piv, d_pivots, dmin, dmax = ncd_pivots(bc, be, bnn)

        node = {
            "id": my_id, "parent": parent_id, "depth": depth,
            "n_in_node": len(local_ev), "shift": float(shift),
            "lo": lo, "hi": hi,
            "ncd_eigs": ncd_max_eigs(lo, hi),
            "ncd_pivots": ncd_piv,
            "pivot_dmin": dmin, "pivot_dmax": dmax,
            "n_pivots": bnn,
            "outcome": None,
            "n_singletons": 0, "n_zerowidth_batches": 0, "n_descend_groups": 0,
            "children": [],
            "tau_to_children": [],
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


def trace_stcoll(basename):
    import glob
    paths = glob.glob(os.path.join(ROOT, "stcollection", f"{basename}.dat"))
    if not paths:
        return None, None
    saved = sys.stdout; sys.stdout = io.StringIO()
    try:
        _, d, e = load_stcoll(paths[0])
    finally:
        sys.stdout = saved
    n = len(d)
    c, e_off, _, sp = build_tgk_offdiag(d, e)
    pos = bisect_evals(c, e_off, len(c), n + 1, 2 * n)
    if len(pos) != n:
        return None, n
    return trace_mr3_block(c, e_off, len(c), pos, sp), n


# ---------- SVG rendering ----------

def render_tree_svg(label, nodes, color_by, max_width_px=560,
                    node_min_r=10, node_max_r=42):
    """color_by ∈ {'ncd_eigs', 'ncd_pivots'}."""
    if not nodes:
        return f'<div class="empty">no nodes for {label}</div>'

    by_id = {n["id"]: n for n in nodes}
    children_of = defaultdict(list)
    for n in nodes:
        if n["parent"] >= 0:
            children_of[n["parent"]].append(n["id"])

    sizes = [n["n_in_node"] for n in nodes]
    max_n = max(sizes)
    def radius(n_in):
        if max_n <= 1:
            return node_min_r
        return node_min_r + (node_max_r - node_min_r) * math.sqrt(n_in / max_n)

    vals = [max(n[color_by], 1e-300) for n in nodes]
    if color_by == "ncd_pivots":
        # NCD_pivots ∈ [0, 1] linear; values often near 1 for nodes that have
        # been shifted, near 0 for unshifted nodes
        def color(v):
            v = max(0.0, min(1.0, v))
            r = int(40 + (250 - 40) * v)
            g = int(40 + (220 - 40) * (1 - abs(v - 0.5) * 2))
            b = int(220 - (220 - 60) * v)
            return f"rgb({r},{g},{b})"
    else:
        # NCD_eigs spans many orders of magnitude → log scale
        log_vals = [math.log10(v) for v in vals]
        cmin = min(log_vals); cmax = max(log_vals)
        span = max(cmax - cmin, 1e-9)
        def color(v):
            t = (math.log10(max(v, 1e-300)) - cmin) / span
            t = max(0.0, min(1.0, t))
            r = int(40 + (250 - 40) * t)
            g = int(40 + (220 - 40) * (1 - abs(t - 0.5) * 2))
            b = int(220 - (220 - 60) * t)
            return f"rgb({r},{g},{b})"

    leaf_counts = {}
    def count_leaves(nid):
        ch = children_of.get(nid, [])
        if not ch:
            leaf_counts[nid] = 1; return 1
        total = sum(count_leaves(c) for c in ch)
        leaf_counts[nid] = total
        return total

    roots = [n["id"] for n in nodes if n["parent"] == -1]
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

    margin_x, margin_y = 40, 40
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

    svg = [f'<svg class="tree-svg" viewBox="0 0 {width} {height}" '
           f'xmlns="http://www.w3.org/2000/svg">',
           f'<text x="{margin_x}" y="20" class="tree-title">'
           f'{label}: {len(nodes)} nodes, max depth {max_depth} (color={color_by})</text>']

    for n in nodes:
        for ci, cid in enumerate(n["children"]):
            x0, y0 = tx(pos[n["id"]]), ty(n["depth"])
            x1, y1 = tx(pos[cid]), ty(n["depth"] + 1)
            cy = y0 + (y1 - y0) * 0.5
            svg.append(
                f'<path d="M{x0:.1f},{y0:.1f} C{x0:.1f},{cy:.1f} '
                f'{x1:.1f},{cy:.1f} {x1:.1f},{y1:.1f}" class="edge" />')
            tau = n["tau_to_children"][ci] if ci < len(n["tau_to_children"]) else None
            if tau is not None:
                svg.append(
                    f'<text x="{(x0+x1)/2:.1f}" y="{cy:.1f}" '
                    f'class="edge-label">τ={tau:.2g}</text>')

    for n in nodes:
        x, y = tx(pos[n["id"]]), ty(n["depth"])
        r = radius(n["n_in_node"])
        col = color(n[color_by])
        title = (
            f'id={n["id"]}  parent={n["parent"]}\n'
            f'depth={n["depth"]}  n_in_node={n["n_in_node"]}\n'
            f'eig range: [{n["lo"]:.4g} .. {n["hi"]:.4g}]\n'
            f'NCD_eigs   = (eig_max - eig_min) / max(|eig|) = {n["ncd_eigs"]:.3g}\n'
            f'NCD_pivots = (|D|max - |D|min) / |D|max       = {n["ncd_pivots"]:.3g}\n'
            f'  pivot |D|min={n["pivot_dmin"]:.3g}  |D|max={n["pivot_dmax"]:.3g}  '
            f'(over {n["n_pivots"]} pivots)\n'
            f'shift τ_total={n["shift"]:.3g}\n'
            f'outcome={n["outcome"]}\n'
            f'sing/zw/desc = {n["n_singletons"]}/{n["n_zerowidth_batches"]}/{n["n_descend_groups"]}'
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
       max-width: 1300px; margin: 1em auto; color: #222; line-height: 1.5;
       padding: 0 1em 4em; }
h1 { color: #06335c; border-left: 4px solid #06c; padding-left: 0.5em; }
h2 { color: #06335c; margin-top: 1.6em; border-bottom: 1px solid #ccc; padding-bottom: 0.2em; }
h3 { color: #244; margin-top: 1.2em; }
.controls { background: #f6f9fc; padding: 0.6em 0.9em; border-radius: 4px;
            font-size: 0.92em; margin-bottom: 1em; }
.tree-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 1em; }
@media (max-width: 880px) { .tree-grid { grid-template-columns: 1fr; } }
.tree-card { background: #fff; border: 1px solid #ddd; border-radius: 5px;
             padding: 0.5em 0.7em; }
.tree-card h3 { margin: 0 0 0.3em 0; font-family: SF Mono, Menlo, Consolas, monospace;
                font-size: 0.95em; color: #024; }
.tree-row { display: grid; grid-template-columns: 1fr 1fr; gap: 0.5em; }
.tree-svg { width: 100%; height: auto; background: #fbfbfd; }
text.tree-title { font-family: SF Mono, Menlo, Consolas, monospace; font-size: 10px; fill: #555; }
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
.detail-box { font-size: 0.74em; max-height: 130px; overflow-y: auto;
              background: #f6f8fa; padding: 0.4em 0.6em; border-radius: 3px;
              font-family: SF Mono, Menlo, Consolas, monospace; }
.banner { background: #fef9e7; border-left: 4px solid #c80; padding: 0.6em 0.9em;
          border-radius: 3px; margin-top: 1em; }
.formula { background: #eef; padding: 0.4em 0.6em; border-radius: 3px;
           font-family: SF Mono, Menlo, Consolas, monospace; }
"""


def card_html(label, nodes):
    if not nodes:
        return f'<div class="tree-card"><h3>{label}</h3><div class="empty">no tree</div></div>'
    md = max(n["depth"] for n in nodes)
    n_leaves = sum(1 for n in nodes if not n["children"])
    root = nodes[0]
    stats = (f'<b>{label}</b>: nodes={len(nodes)} depth={md} leaves={n_leaves} &nbsp; '
             f'root NCD_eigs={root["ncd_eigs"]:.3g} &nbsp; '
             f'root NCD_pivots={root["ncd_pivots"]:.3g}')
    detail = []
    for n in nodes[:25]:
        detail.append(
            f'id={n["id"]:>3}  d={n["depth"]}  n={n["n_in_node"]:>4}  '
            f'NCD_eigs={n["ncd_eigs"]:.2e}  '
            f'NCD_pivots={n["ncd_pivots"]:.2e}  '
            f'|D|min={n["pivot_dmin"]:.2e}  |D|max={n["pivot_dmax"]:.2e}  '
            f'sing={n["n_singletons"]:<3}  zw={n["n_zerowidth_batches"]:<2}  '
            f'desc={n["n_descend_groups"]:<2}  out={n["outcome"]}'
        )
    if len(nodes) > 25:
        detail.append(f'... ({len(nodes) - 25} more)')
    detail_html = "<pre class='detail-box'>" + "\n".join(detail) + "</pre>"
    svg_eigs = render_tree_svg(label + " — color: NCD_eigs",
                                nodes, "ncd_eigs")
    svg_pivs = render_tree_svg(label + " — color: NCD_pivots",
                                nodes, "ncd_pivots")
    return (f'<div class="tree-card"><h3>{stats}</h3>'
            f'<div class="tree-row">{svg_eigs}{svg_pivs}</div>'
            f'{detail_html}</div>')


def main():
    cases_easy = [
        ("diagonal_only",          50),
        ("constant",               50),
        ("arithmetic_progression", 50),
    ]
    cases_medium = [
        ("step_function",          200),
        ("two_clusters",           400),
        ("three_clusters",         400),
        ("demmel_S1pe_k4",         400),
        ("demmel_S1pe_k8",         400),
        ("saw_tooth",              400),
    ]
    cases_hard = [
        ("one_big_cluster",        400),
        ("one_big_cluster",        200),
        ("random_clustered_5",     400),
        ("constant_d_graded_e",    400),
        ("gl_clustered_at_eps",    400),
        ("pd_T0",                  100),
    ]
    cases_wilk = [
        ("wilkinson_exact",        200),
        ("wilkinson_like",         200),
        ("glued_wilkinson",        200),
        ("glued_wilkinson_tight",  200),
        ("gl_wilkw",               400),
        ("gl_wilkm",               400),
    ]
    cases_stcoll = [
        ("ST", "B_Kimura_429"),
        ("ST", "B_gg_30_1D-5"),
        ("ST", "B_20_graded"),
        ("ST", "B_40_graded"),
    ]

    sections = [
        ("Easy: well-separated — single-root trees", cases_easy),
        ("Medium: clustered subspaces", cases_medium),
        ("Hard: deep clusters", cases_hard),
        ("Wilkinson family", cases_wilk),
        ("STCollection (real-world)", cases_stcoll),
    ]

    all_traces = {}
    html_sections = []
    for sec_title, cases in sections:
        cards = []
        for entry in cases:
            if entry[0] == "ST":
                _, basename = entry
                t, n = trace_stcoll(basename)
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
                    continue
            md = max(nd["depth"] for nd in t)
            print(f"  {label}: {len(t)} nodes, depth {md}, "
                  f"root NCD_eigs={t[0]['ncd_eigs']:.3g}, "
                  f"root NCD_pivots={t[0]['ncd_pivots']:.3g}")
            all_traces[label] = t
            cards.append(card_html(label, t))
        html_sections.append((sec_title, cards))

    body_parts = ['<h1>MR3 Representation-Tree Gallery — dual NCD</h1>',
        '<div class="controls">'
        '<p>Each tree node is shown twice — left tree colored by '
        '<b>NCD_eigs</b>, right tree colored by <b>NCD_pivots</b>. Hover any '
        'circle to see both NCDs plus the LDL^T pivot summary.</p>'
        '<div class="formula"><b>NCD_eigs</b> &nbsp;=&nbsp; '
        '(λ<sub>max</sub> − λ<sub>min</sub>) / max(|λ|) &nbsp;'
        '(eigenvalue range; what XMR uses for split decisions)</div>'
        '<div class="formula"><b>NCD_pivots</b> &nbsp;=&nbsp; '
        '(|D|<sub>max</sub> − |D|<sub>min</sub>) / |D|<sub>max</sub>'
        '&nbsp; where T<sub>node</sub> = LDL<sup>T</sup> '
        '(pivot dynamic range after the LDL<sup>T</sup> factorization '
        'of the shifted T_GK at the node)</div></div>',
        '<div class="legend">'
        '<span><span class="swatch" style="background:rgb(40,40,220)"></span>low NCD</span>'
        '<span><span class="swatch" style="background:rgb(150,220,150)"></span>medium</span>'
        '<span><span class="swatch" style="background:rgb(250,40,60)"></span>high NCD (≈1)</span>'
        '</div>']
    for sec_title, cards in html_sections:
        body_parts.append(f'<h2>{sec_title}</h2>')
        for c in cards:
            body_parts.append(c)
    body_parts.append('<div class="banner">'
        '<b>Why two metrics?</b> NCD_eigs (eigenvalue range) is what XMR '
        'actually consumes via <code>dlaxrb_clssfy</code> when deciding '
        'which sub-cluster to split. NCD_pivots is the user-requested '
        '"diag-range of LDL<sup>T</sup>" view — it measures the same '
        'clustered structure from the other side: tightly-clustered '
        'eigenvalues produce widely-varying pivots in the LDL<sup>T</sup> '
        'factorization (large pivot range → high NCD_pivots ≈ 1).</div>')

    body = "\n".join(body_parts)
    html = (f"<!doctype html><html><head><meta charset='utf-8'>"
            f"<title>MR3 Tree Gallery — dual NCD</title>"
            f"<style>{CSS}</style></head>"
            f"<body>{body}</body></html>")

    os.makedirs(os.path.dirname(OUT_HTML), exist_ok=True)
    with open(OUT_HTML, "w") as f:
        f.write(html)
    with open(OUT_JSON, "w") as f:
        json.dump(all_traces, f, indent=1, default=float)
    sz = os.path.getsize(OUT_HTML) // 1024
    print(f"\nwrote {OUT_HTML}  ({sz} KB, {len(all_traces)} trees, dual-NCD)")


if __name__ == "__main__":
    main()
