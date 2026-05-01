#!/usr/bin/env python3
"""
Cluster-tree NCD trace on a highly-clustered bidiagonal matrix.

Drives the Python `mr3_block` MR3 fallback on a chosen test matrix and logs:
  * EVERY tree node (depth, parent, n_in_node, lo, hi, ncd_max, ncd_min,
    classification outcome).
  * EVERY adjacent pair within every node (the actual NCD values that
    `classify` consumes).

Outputs in experiments/:
  cluster_tree_log.json         full per-node + per-pair log
  cluster_tree_log.csv          flat CSV (rows = nodes; columns include the
                                two NCD aggregates and adjacent-pair stats)
  cluster_tree.png              tree visualization (nodes coloured by
                                log10(NCD_min); size ∝ n_in_node)
  ncd_per_adjacent_pair.png     scatter of all adjacent-pair NCDs
                                (max-form vs min-form, coloured by depth)
  ncd_progression.png           per-depth aggregate scatter (legacy panel)

Usage:
  python3 experiments/cluster_tree_trace.py                            # default cases
  python3 experiments/cluster_tree_trace.py gl_clustered_at_eps 200    # one case
"""
import os
import sys
import io
import csv
import json
import math

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

# Suppress full_eval / mr3_gk import-time prints
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

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def ncd_max(lo, hi):
    spread = hi - lo
    denom = max(abs(lo), abs(hi))
    return spread / denom if denom > 0 else 0.0


def ncd_min(lo, hi):
    spread = hi - lo
    denom = max(min(abs(lo), abs(hi)), SAFMIN)
    return spread / denom


def trace_mr3_block(c, e_off, nn, pos_evals, spdiam):
    """Instrumented copy of mr3_gk.mr3_block.

    Returns a list of node dicts with full per-pair NCD logs and
    classification outcomes.
    """
    nev = len(pos_evals)
    result_evals = np.zeros(nev)
    result_evecs = np.zeros((nn, nev))
    next_id = [0]

    def new_id():
        i = next_id[0]; next_id[0] = i + 1
        return i

    root_id = new_id()
    stack = [(c[:nn].copy(), e_off[:max(nn - 1, 0)].copy(), nn,
              pos_evals.copy(), list(range(nev)), 0.0, 0, root_id, -1)]
    nodes = []
    while stack:
        bc, be, bnn, local_ev, out_idx, shift, depth, my_id, parent_id = stack.pop()
        if len(local_ev) == 0:
            continue
        lo = float(np.min(local_ev)); hi = float(np.max(local_ev))
        # per-adjacent-pair NCD
        pairs = []
        for k in range(len(local_ev) - 1):
            a = float(local_ev[k]); b = float(local_ev[k + 1])
            pairs.append({
                "i": k,
                "lo": a, "hi": b,
                "ncd_max": ncd_max(a, b),
                "ncd_min": ncd_min(a, b),
                "would_split": ncd_max(a, b) >= GAPTOL,  # mirrors mr3_gk.classify
            })

        node = {
            "id": my_id,
            "parent": parent_id,
            "depth": depth,
            "n_in_node": len(local_ev),
            "shift": shift,
            "lo": lo, "hi": hi,
            "ncd_max": ncd_max(lo, hi),
            "ncd_min": ncd_min(lo, hi),
            "outcome": None,        # filled in below
            "n_singletons": 0,
            "n_zerowidth_batches": 0,
            "n_descend_groups": 0,
            "children": [],
            "adjacent_pairs": pairs,
        }

        if depth >= MAX_DEPTH:
            node["outcome"] = "max_depth_singletons"
            mult = len(local_ev)
            for j in range(mult):
                lam_ref, z = compute_singleton(
                    bc, be, bnn, local_ev[j],
                    all_evals=local_ev, eval_idx=j, depth=depth,
                    twist_rank=j, multiplicity=mult)
                result_evals[out_idx[j]] = lam_ref + shift
                result_evecs[:bnn, out_idx[j]] = z
            nodes.append(node)
            continue

        groups = classify(local_ev, GAPTOL)
        node["n_groups_after_classify"] = len(groups)
        for gs, ge in groups:
            gev = local_ev[gs:ge + 1]
            gidx = out_idx[gs:ge + 1]
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
                        eval_idx=gs + j, depth=depth,
                        twist_rank=j, multiplicity=mult)
                    result_evals[gidx[j]] = lam_ref + shift
                    result_evecs[:bnn, gidx[j]] = z
                continue
            node["n_descend_groups"] += 1
            tau, cc, ee = select_shift_geometric(
                bc, be, bnn, gev, spdiam, local_ev)
            shifted_ev = gev - tau
            refined = np.array([refine_eval(cc, ee, bnn, s)
                                for s in shifted_ev])
            si = np.argsort(refined)
            refined = refined[si]
            sgidx = [gidx[i] for i in si]
            child_id = new_id()
            node["children"].append(child_id)
            stack.append((cc, ee, bnn, refined, sgidx,
                          shift + tau, depth + 1, child_id, my_id))

        node["outcome"] = (
            "leaf"
            if (node["n_descend_groups"] == 0)
            else "internal"
        )
        nodes.append(node)
    return nodes


def build_tgk_offdiag(d, e):
    n = len(d); m = 2 * n
    c = np.zeros(m, dtype=np.float64)
    e_off = np.zeros(m - 1, dtype=np.float64)
    for i in range(n):
        e_off[2 * i] = abs(float(d[i]))
    for i in range(n - 1):
        e_off[2 * i + 1] = abs(float(e[i]))
    spdiam = float(np.max(np.abs(e_off))) if m > 1 else 0.0
    return c, e_off, n, spdiam


def trace_matrix(name, n):
    d, e = make(name, n)
    c, e_off, n_nonneg, spdiam = build_tgk_offdiag(d, e)
    pos_evals = bisect_evals(c, e_off, len(c), n + 1, 2 * n)
    if len(pos_evals) != n:
        print(f"  WARNING ({name}@{n}): bisect returned {len(pos_evals)} evs, expected {n}")
        return None
    return trace_mr3_block(c, e_off, len(c), pos_evals, spdiam)


# ---------- writers / plots ----------

def write_json(traces, path):
    safe = {}
    for label, nodes in traces.items():
        safe[label] = nodes  # already JSON-friendly
    with open(path, "w") as f:
        json.dump(safe, f, indent=1, default=float)


def write_csv(traces, path):
    rows = []
    for label, nodes in traces.items():
        for n in nodes:
            rows.append({
                "matrix": label,
                "id": n["id"],
                "parent": n["parent"],
                "depth": n["depth"],
                "n_in_node": n["n_in_node"],
                "shift": n["shift"],
                "lo": n["lo"], "hi": n["hi"],
                "ncd_max": n["ncd_max"],
                "ncd_min": n["ncd_min"],
                "outcome": n["outcome"],
                "n_singletons": n["n_singletons"],
                "n_zerowidth_batches": n["n_zerowidth_batches"],
                "n_descend_groups": n["n_descend_groups"],
                "n_adjacent_pairs": len(n["adjacent_pairs"]),
                "max_pair_ncd_max": max((p["ncd_max"] for p in n["adjacent_pairs"]), default=0.0),
                "min_pair_ncd_max": min((p["ncd_max"] for p in n["adjacent_pairs"]), default=0.0),
                "max_pair_ncd_min": max((p["ncd_min"] for p in n["adjacent_pairs"]), default=0.0),
                "min_pair_ncd_min": min((p["ncd_min"] for p in n["adjacent_pairs"]), default=0.0),
            })
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else
                           ["matrix", "id", "depth"])
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _layout_tree(nodes):
    """Return dict id -> (x, y) by simple depth-leveled layout."""
    by_depth = {}
    for n in nodes:
        by_depth.setdefault(n["depth"], []).append(n)
    pos = {}
    for d, ns in by_depth.items():
        ns.sort(key=lambda n: n["id"])
        m = len(ns)
        for i, n in enumerate(ns):
            x = (i + 1) / (m + 1)
            y = -d
            pos[n["id"]] = (x, y)
    return pos


def plot_tree(traces, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(f"  matplotlib missing; skipping {path}")
        return
    items = list(traces.items())
    if not items:
        return
    fig, axes = plt.subplots(1, len(items), figsize=(7 * len(items), 6),
                             squeeze=False)
    for ax, (label, nodes) in zip(axes[0], items):
        pos = _layout_tree(nodes)
        # edges
        for n in nodes:
            for cid in n["children"]:
                if cid in pos:
                    x0, y0 = pos[n["id"]]; x1, y1 = pos[cid]
                    ax.plot([x0, x1], [y0, y1], "k-", alpha=0.3, lw=0.7)
        # nodes
        ncd_min_log = []
        sizes = []
        xs, ys = [], []
        for n in nodes:
            x, y = pos[n["id"]]
            xs.append(x); ys.append(y)
            v = n["ncd_min"] if n["ncd_min"] > 0 else 1e-300
            ncd_min_log.append(math.log10(v))
            sizes.append(20 + 6 * n["n_in_node"])
        sc = ax.scatter(xs, ys, c=ncd_min_log, s=sizes, cmap="viridis",
                        edgecolor="k", linewidth=0.4, alpha=0.85)
        for n in nodes:
            x, y = pos[n["id"]]
            ax.annotate(f"id{n['id']}\nn={n['n_in_node']}",
                        (x, y), fontsize=6, ha="center", va="center")
        ax.set_title(f"{label}\n{len(nodes)} nodes, max depth {max(n['depth'] for n in nodes)}")
        ax.set_xlabel("x (within depth)")
        ax.set_ylabel("-depth")
        ax.set_xlim(-0.05, 1.05)
        plt.colorbar(sc, ax=ax, label="log10 NCD_min")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    print(f"  wrote {path}")


def plot_pair_ncd(traces, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for label, nodes in traces.items():
        depth = []
        nm = []; nx = []
        for n in nodes:
            for p in n["adjacent_pairs"]:
                depth.append(n["depth"])
                nx.append(max(p["ncd_max"], 1e-300))
                nm.append(max(p["ncd_min"], 1e-300))
        if not depth:
            continue
        axes[0].scatter(depth, nx, alpha=0.4, s=12, label=label)
        axes[1].scatter(depth, nm, alpha=0.4, s=12, label=label)
    axes[0].axhline(GAPTOL, color="r", linestyle="--", alpha=0.5,
                    label=f"GAPTOL={GAPTOL:g}")
    axes[1].axhline(GAPTOL, color="r", linestyle="--", alpha=0.5,
                    label=f"GAPTOL={GAPTOL:g}")
    for ax, t in zip(axes,
                     ["adjacent-pair NCD_max = (hi-lo)/max(|lo|,|hi|)",
                      "adjacent-pair NCD_min = (hi-lo)/min(|lo|,|hi|)"]):
        ax.set_xlabel("depth")
        ax.set_ylabel(t)
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
        ax.set_title(t)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    print(f"  wrote {path}")


def plot_progression(traces, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for label, nodes in traces.items():
        depths = [n["depth"] for n in nodes if n["n_in_node"] > 1]
        nx = [n["ncd_max"] for n in nodes if n["n_in_node"] > 1]
        nm = [n["ncd_min"] for n in nodes if n["n_in_node"] > 1]
        axes[0].scatter(depths, nx, alpha=0.6, s=18, label=label)
        axes[1].scatter(depths, nm, alpha=0.6, s=18, label=label)
    for ax, t in zip(axes,
                     ["per-node NCD_max = (max-min)/max",
                      "per-node NCD_min = (max-min)/min"]):
        ax.set_xlabel("depth")
        ax.set_ylabel(t)
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=9)
        ax.set_title(t)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    print(f"  wrote {path}")


def main():
    if len(sys.argv) == 3:
        cases = [(sys.argv[1], int(sys.argv[2]))]
    else:
        cases = [
            ("gl_clustered_at_eps", 200),
            ("step_function", 200),
            ("demmel_S1pe_k8", 100),
            ("pd_T0", 100),
        ]
    traces = {}
    for name, n in cases:
        try:
            t = trace_matrix(name, n)
        except Exception as ex:
            print(f"  {name}@{n}: ERROR {ex}")
            continue
        if t is not None:
            label = f"{name}@{n}"
            traces[label] = t
            print(f"  {label}: {len(t)} nodes, depth 0..{max(n['depth'] for n in t)}")

    if not traces:
        print("no traces produced")
        return

    write_json(traces, os.path.join(OUT_DIR, "cluster_tree_log.json"))
    write_csv(traces, os.path.join(OUT_DIR, "cluster_tree_log.csv"))
    plot_tree(traces, os.path.join(OUT_DIR, "cluster_tree.png"))
    plot_pair_ncd(traces, os.path.join(OUT_DIR, "ncd_per_adjacent_pair.png"))
    plot_progression(traces, os.path.join(OUT_DIR, "ncd_progression.png"))


if __name__ == "__main__":
    main()
