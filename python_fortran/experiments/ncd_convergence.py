#!/usr/bin/env python3
"""
NCD convergence analysis on clustered bidiagonal matrices.

For each matrix we drive the Python `mr3_block` MR3 fallback (algorithmically
equivalent to the Fortran kernel) and record:
  - parent → child NCD shrink (both (max-min)/max and (max-min)/min)
  - per-leaf getvec ratio = |residual / gap|  (predicts ortU breakdown)
  - actual ortU on the SVD output

Outputs in experiments/:
  ncd_convergence_log.json   per-matrix per-edge log
  ncd_convergence.png        4-panel figure (shrink, getvec, ortU, summary)
"""
import os
import sys
import io
import json
import math

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

# Suppress import-time chatter
_saved = sys.stdout; sys.stdout = io.StringIO()
try:
    import mr3_gk
    from mr3_gk import (
        bisect_evals, classify, is_zero_width,
        select_shift_geometric, refine_eval, compute_singleton, bidiag_svd,
        GAPTOL, MAX_DEPTH, SAFMIN, EPS,
    )
    from full_eval import make
finally:
    sys.stdout = _saved

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def ncd_max(lo, hi):
    sp = hi - lo; d = max(abs(lo), abs(hi))
    return sp / d if d > 0 else 0.0


def ncd_min(lo, hi):
    sp = hi - lo; d = max(min(abs(lo), abs(hi)), SAFMIN)
    return sp / d


def trace_with_metrics(c, e_off, nn, pos_evals, spdiam):
    """Variant of mr3_block that records:
      - one record per node visit (depth, lo, hi, ncd_*)
      - one record per parent->child edge (parent_ncd, child_ncd, shift)
      - one record per leaf singleton (residual, gap, ratio)
    """
    nev = len(pos_evals)
    result_evals = np.zeros(nev)
    result_evecs = np.zeros((nn, nev))
    next_id = [0]

    def new_id():
        i = next_id[0]; next_id[0] = i + 1; return i

    nodes = []
    edges = []
    leaves = []

    root_id = new_id()
    stack = [(c[:nn].copy(), e_off[:max(nn-1, 0)].copy(), nn,
              pos_evals.copy(), list(range(nev)), 0.0, 0, root_id, -1,
              None)]  # parent NCD pair (lo, hi)
    while stack:
        bc, be, bnn, local_ev, out_idx, shift, depth, my_id, parent_id, parent_lohi = stack.pop()
        if len(local_ev) == 0:
            continue
        lo = float(np.min(local_ev)); hi = float(np.max(local_ev))
        node_rec = {
            "id": my_id, "parent": parent_id, "depth": depth,
            "n": len(local_ev),
            "lo": lo, "hi": hi,
            "ncd_max": ncd_max(lo, hi),
            "ncd_min": ncd_min(lo, hi),
            "shift": shift,
        }
        nodes.append(node_rec)

        if depth >= MAX_DEPTH:
            for j in range(len(local_ev)):
                lam_ref, z = compute_singleton(
                    bc, be, bnn, local_ev[j],
                    all_evals=local_ev, eval_idx=j, depth=depth,
                    twist_rank=j, multiplicity=len(local_ev))
                result_evals[out_idx[j]] = lam_ref + shift
                result_evecs[:bnn, out_idx[j]] = z
            continue

        groups = classify(local_ev, GAPTOL)
        for gs, ge in groups:
            gev = local_ev[gs:ge+1]; gidx = out_idx[gs:ge+1]
            if gs == ge:
                lam_ref, z = compute_singleton(
                    bc, be, bnn, gev[0], all_evals=local_ev,
                    eval_idx=gs, depth=depth)
                result_evals[gidx[0]] = lam_ref + shift
                result_evecs[:bnn, gidx[0]] = z
                # Compute residual and gap for getvec ratio
                tz = bc[:bnn] * z
                if bnn > 1:
                    tz[:bnn-1] += be[:bnn-1] * z[1:]
                    tz[1:]     += be[:bnn-1] * z[:bnn-1]
                resid = float(np.linalg.norm(tz - lam_ref * z))
                # local gap: smallest |lam_ref - other|
                others = np.delete(local_ev, gs)
                if len(others) > 0:
                    gap = float(np.min(np.abs(others - lam_ref)))
                else:
                    gap = float("inf")
                leaves.append({
                    "node_id": my_id, "depth": depth,
                    "lambda": float(lam_ref + shift),
                    "residual": resid, "gap": gap,
                    "getvec_ratio": resid / max(gap, SAFMIN),
                })
                continue
            if is_zero_width(gev):
                for j in range(len(gev)):
                    lam_ref, z = compute_singleton(
                        bc, be, bnn, gev[j], all_evals=local_ev,
                        eval_idx=gs+j, depth=depth,
                        twist_rank=j, multiplicity=len(gev))
                    result_evals[gidx[j]] = lam_ref + shift
                    result_evecs[:bnn, gidx[j]] = z
                continue
            # Non-singleton, not zero-width → descend
            tau, cc, ee = select_shift_geometric(bc, be, bnn, gev, spdiam, local_ev)
            shifted_ev = gev - tau
            refined = np.array([refine_eval(cc, ee, bnn, s) for s in shifted_ev])
            si = np.argsort(refined)
            refined = refined[si]
            sgidx = [gidx[i] for i in si]
            child_id = new_id()
            child_lo = float(refined[0]); child_hi = float(refined[-1])
            edges.append({
                "parent": my_id, "child": child_id, "depth": depth,
                "tau": float(tau),
                "parent_ncd_max": ncd_max(float(gev[0]), float(gev[-1])),
                "parent_ncd_min": ncd_min(float(gev[0]), float(gev[-1])),
                "child_ncd_max": ncd_max(child_lo, child_hi),
                "child_ncd_min": ncd_min(child_lo, child_hi),
                "n": len(refined),
            })
            stack.append((cc, ee, bnn, refined, sgidx,
                          shift + tau, depth + 1, child_id, my_id,
                          (float(gev[0]), float(gev[-1]))))
    return nodes, edges, leaves


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


def measure_ortu_ortv(d, e):
    """Run the full bidiag_svd pipeline; return (res, ortU, ortV) in n*EPS units."""
    n = len(d)
    sigma, U, V, info = bidiag_svd(d.copy(), e.copy())
    B = np.diag(d) + np.diag(e, 1)
    bnorm = abs(d[-1]) if n > 0 else 0.0
    for i in range(n - 1):
        bnorm = max(bnorm, abs(d[i]) + abs(e[i]))
    if bnorm == 0.0:
        bnorm = 1.0
    res = float(np.max(np.abs(B - U @ np.diag(sigma) @ V.T)) / (bnorm * n * EPS))
    ou  = float(np.max(np.abs(U.T @ U - np.eye(n))) / (n * EPS))
    ov  = float(np.max(np.abs(V.T @ V - np.eye(n))) / (n * EPS))
    return res, ou, ov


def analyze_one(name, n):
    d, e = make(name, n)
    c, e_off, n_nonneg, spdiam = build_tgk_offdiag(d, e)
    pos_evals = bisect_evals(c, e_off, len(c), n + 1, 2 * n)
    if len(pos_evals) != n:
        return None
    nodes, edges, leaves = trace_with_metrics(c, e_off, len(c), pos_evals, spdiam)
    res, ou, ov = measure_ortu_ortv(d, e)
    return {
        "name": name, "n": n,
        "max_depth": max((nd["depth"] for nd in nodes), default=0),
        "n_nodes": len(nodes), "n_edges": len(edges), "n_leaves": len(leaves),
        "median_ncd_max_shrink": float(np.median(
            [ed["parent_ncd_max"] / max(ed["child_ncd_max"], SAFMIN) for ed in edges]
        )) if edges else float("nan"),
        "median_ncd_min_shrink": float(np.median(
            [ed["parent_ncd_min"] / max(ed["child_ncd_min"], SAFMIN) for ed in edges]
        )) if edges else float("nan"),
        "max_getvec_ratio": float(max((l["getvec_ratio"] for l in leaves), default=0.0)),
        "median_getvec_ratio": float(np.median(
            [l["getvec_ratio"] for l in leaves]
        )) if leaves else float("nan"),
        "res_neps": res, "ortU_neps": ou, "ortV_neps": ov,
        "nodes": nodes, "edges": edges, "leaves": leaves,
    }


def make_plot(records, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(f"  matplotlib unavailable; skipping {path}")
        return
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # Panel A: parent NCD vs child NCD (NCD_min) — log-log scatter, color = depth
    axA = axes[0, 0]
    for rec in records:
        for ed in rec["edges"]:
            axA.scatter(max(ed["parent_ncd_min"], 1e-300),
                        max(ed["child_ncd_min"], 1e-300),
                        c=ed["depth"], cmap="viridis",
                        s=22, alpha=0.7, vmin=0, vmax=3)
    lo, hi = 1e-3, 1e18
    axA.plot([lo, hi], [lo, hi], "k--", alpha=0.4, lw=0.7,
             label="y=x (no shrink)")
    axA.set_xscale("log"); axA.set_yscale("log")
    axA.set_xlim(lo, hi); axA.set_ylim(lo, hi)
    axA.set_xlabel("parent NCD_min"); axA.set_ylabel("child NCD_min")
    axA.set_title("Parent→Child NCD_min (color = parent depth)")
    axA.grid(True, which="both", alpha=0.3); axA.legend(fontsize=8)

    # Panel B: getvec ratio vs cluster size, log
    axB = axes[0, 1]
    for rec in records:
        depths = [l["depth"] for l in rec["leaves"]]
        ratios = [max(l["getvec_ratio"], 1e-300) for l in rec["leaves"]]
        if depths:
            axB.scatter(depths, ratios, alpha=0.5, s=14,
                        label=f'{rec["name"]}@{rec["n"]}')
    axB.axhline(100 * EPS, color="r", linestyle="--", alpha=0.6,
                label="100·EPS (orth target)")
    axB.set_yscale("log")
    axB.set_xlabel("leaf depth")
    axB.set_ylabel("getvec ratio = |residual| / gap")
    axB.set_title("Per-leaf getvec ratio  (predicts ortU)")
    axB.grid(True, which="both", alpha=0.3); axB.legend(fontsize=7)

    # Panel C: ortU achieved vs max getvec ratio per matrix
    axC = axes[1, 0]
    xs = [max(rec["max_getvec_ratio"], 1e-300) for rec in records]
    ys = [max(rec["ortU_neps"], 1e-3) for rec in records]
    labels = [f'{rec["name"]}@{rec["n"]}' for rec in records]
    axC.scatter(xs, ys, s=40, alpha=0.7, edgecolor="k")
    for x, y, lab in zip(xs, ys, labels):
        axC.annotate(lab, (x, y), fontsize=7, alpha=0.8,
                     xytext=(3, 3), textcoords="offset points")
    axC.set_xscale("log"); axC.set_yscale("log")
    axC.axhline(5.0, color="g", linestyle="--", alpha=0.6,
                label="evaluate.py ortU threshold = 5 nε")
    axC.set_xlabel("max getvec ratio across all leaves")
    axC.set_ylabel("ortU (n·EPS units)")
    axC.set_title("Predicted vs achieved orthogonality")
    axC.grid(True, which="both", alpha=0.3); axC.legend(fontsize=8)

    # Panel D: summary table as text
    axD = axes[1, 1]; axD.axis("off")
    rows = [["matrix", "depth", "edges", "leaves",
             "med shrink (max)", "med shrink (min)",
             "max getvec",
             "ortU"]]
    for rec in records:
        rows.append([f'{rec["name"]}@{rec["n"]}',
                     str(rec["max_depth"]),
                     str(rec["n_edges"]),
                     str(rec["n_leaves"]),
                     f'{rec["median_ncd_max_shrink"]:.2g}',
                     f'{rec["median_ncd_min_shrink"]:.2g}',
                     f'{rec["max_getvec_ratio"]:.2e}',
                     f'{rec["ortU_neps"]:.2g}'])
    tbl = axD.table(cellText=rows[1:], colLabels=rows[0],
                    loc="center", cellLoc="center",
                    colColours=["#ddd"] * len(rows[0]))
    tbl.auto_set_font_size(False); tbl.set_fontsize(7)
    tbl.scale(1.0, 1.4)
    axD.set_title("Per-matrix summary", y=0.96)

    fig.tight_layout()
    fig.savefig(path, dpi=120)
    print(f"  wrote {path}")


def main():
    cases = [
        ("gl_clustered_at_eps", 200),
        ("step_function", 200),
        ("demmel_S1pe_k8", 100),
        ("pd_T0", 100),
        ("three_clusters", 200),
        ("demmel_S1pe_k4", 200),
        ("gl_abcon3", 100),
        ("gl_wilkw", 200),
    ]
    records = []
    for name, n in cases:
        try:
            r = analyze_one(name, n)
        except Exception as ex:
            print(f"  {name}@{n}: ERROR {ex}")
            continue
        if r is None:
            print(f"  {name}@{n}: bisect mismatch — skipped")
            continue
        records.append(r)
        print(f"  {r['name']}@{r['n']}: depth={r['max_depth']}  edges={r['n_edges']}  "
              f"leaves={r['n_leaves']}  med_shrink_max={r['median_ncd_max_shrink']:.2g}  "
              f"max_getvec={r['max_getvec_ratio']:.2e}  ortU={r['ortU_neps']:.2g}")

    out_log = os.path.join(OUT_DIR, "ncd_convergence_log.json")
    with open(out_log, "w") as f:
        json.dump(records, f, indent=1, default=float)
    print(f"\nwrote {out_log}")
    make_plot(records, os.path.join(OUT_DIR, "ncd_convergence.png"))


if __name__ == "__main__":
    main()
