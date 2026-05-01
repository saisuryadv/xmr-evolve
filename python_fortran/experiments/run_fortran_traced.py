#!/usr/bin/env python3
"""
Drive the instrumented Fortran SVD `mr3gk_run_traced` across a curated
list of test bidiagonals and parse the per-tree-node NCD logs.

Output: experiments/mr3_tree_fortran_log.json — keyed by "{name}@{n}",
each entry is a list of node dicts captured straight from the Fortran
kernel (one per successful child representation built by dlaxrf_traced).

Per node fields (all from Fortran-side computation):
  depth, n, icbeg, icend, tau,
  pivot_dmin, pivot_dmax,
  eig_lo, eig_hi,
  ncd_eigs   = (eig_hi - eig_lo) / max(|eig_lo|, |eig_hi|)
  ncd_pivots = (pivot_dmax - pivot_dmin) / pivot_dmax
"""
import os
import re
import sys
import io
import json
import struct
import subprocess
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

VARIANTS = {
    "blocks":   os.path.join(ROOT, "mr3gk_fortran", "mr3gk_run_traced"),
    "noblocks": os.path.join(ROOT, "mr3gk_fortran", "mr3gk_run_traced_noblk"),
}
MR3GK_RUN_TRACED = VARIANTS["blocks"]   # default for back-compat
OUT_JSON = os.path.join(ROOT, "experiments", "mr3_tree_fortran_log.json")


def silent_imports():
    saved = sys.stdout; sys.stdout = io.StringIO()
    try:
        from full_eval import make
        from evaluate import load_stcoll
    finally:
        sys.stdout = saved
    return make, load_stcoll


_NODE_RE = re.compile(
    r"^NODE\s+depth=(?P<depth>-?\d+)"
    r"\s+n=(?P<n>\d+)"
    r"\s+icbeg=(?P<icbeg>\d+)"
    r"\s+icend=(?P<icend>\d+)"
    r"\s+tau=(?P<tau>[+\-\d.E ]+?)"
    r"\s+recon_dmin=(?P<recon_dmin>[+\-\d.E ]+?)"
    r"\s+recon_dmax=(?P<recon_dmax>[+\-\d.E ]+?)"
    r"\s+ncd_recon=(?P<ncd_recon>[+\-\d.E ]+?)"
    r"\s+nblocks=(?P<nblocks>\d+)"
    r"\s+twist=(?P<twist>\d+)\s*$"
)


def parse_trace(path):
    """Parse the trace log written by mr3gk_run_traced."""
    nodes = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            m = _NODE_RE.match(ln)
            if not m:
                continue
            d = m.groupdict()
            nodes.append({
                "depth": int(d["depth"]),
                "n":     int(d["n"]),
                "icbeg": int(d["icbeg"]),
                "icend": int(d["icend"]),
                "tau":           float(d["tau"]),
                "recon_dmin":    float(d["recon_dmin"]),
                "recon_dmax":    float(d["recon_dmax"]),
                "ncd_recon":     float(d["ncd_recon"]),
                "nblocks":       int(d["nblocks"]),
                "twist":         int(d["twist"]),
                "n_in_node":     int(d["icend"]) - int(d["icbeg"]) + 1,
            })
    return nodes


def write_bin(path, d, e):
    n = len(d)
    with open(path, "wb") as f:
        f.write(struct.pack("i", n))
        f.write(np.asarray(d, dtype=np.float64).tobytes())
        if n > 1:
            f.write(np.asarray(e, dtype=np.float64).tobytes())


def run_traced(d, e, label, variant="blocks", timeout=300):
    """Run mr3gk_run_traced[_noblk] on (d, e). Return parsed nodes."""
    binary = VARIANTS[variant]
    in_path = "/tmp/mr3gk_traced_in.bin"
    out_path = "/tmp/mr3gk_traced_out.bin"
    safe = label.replace('/','_').replace('@','_at_')
    trace_path = f"/tmp/mr3gk_traced_{variant}_{safe}.log"
    write_bin(in_path, d, e)
    try:
        r = subprocess.run([binary, in_path, out_path, trace_path],
                           capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return [], "timeout"
    if r.returncode != 0:
        return [], f"crash rc={r.returncode}"
    if not os.path.exists(trace_path):
        return [], "no_trace"
    nodes = parse_trace(trace_path)
    return nodes, None


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", choices=["blocks", "noblocks", "both"],
                    default="both")
    args = ap.parse_args()
    variants_to_run = ["blocks", "noblocks"] if args.variant == "both" else [args.variant]

    for v in variants_to_run:
        if not os.path.exists(VARIANTS[v]):
            print(f"ERROR: {VARIANTS[v]} not built.")
            sys.exit(1)

    make, load_stcoll = silent_imports()

    cases = [
        ("diagonal_only",          50),
        ("constant",               50),
        ("arithmetic_progression", 50),
        ("step_function",          200),
        ("two_clusters",           400),
        ("three_clusters",         400),
        ("demmel_S1pe_k4",         400),
        ("demmel_S1pe_k8",         400),
        ("saw_tooth",              400),
        ("one_big_cluster",        400),
        ("random_clustered_5",     400),
        ("constant_d_graded_e",    400),
        ("gl_clustered_at_eps",    400),
        ("pd_T0",                  100),
        ("wilkinson_exact",        200),
        ("wilkinson_like",         200),
        ("glued_wilkinson",        200),
        ("glued_wilkinson_tight",  200),
        ("gl_wilkw",               400),
        ("gl_wilkm",               400),
    ]
    stcoll = ["B_Kimura_429", "B_gg_30_1D-5", "B_20_graded", "B_40_graded"]

    all_traces = {}  # {label: {variant: {nodes, ...}}}

    def add(label, n_mat, d_arr, e_arr):
        per_variant = {}
        for v in variants_to_run:
            nodes, err = run_traced(d_arr, e_arr, label, variant=v)
            if err:
                print(f"  {label} [{v}]: {err}")
                continue
            md = max((nd["depth"] for nd in nodes), default=0)
            per_variant[v] = {
                "nodes": nodes,
                "n_in_matrix": int(n_mat),
                "n_nodes": len(nodes),
                "max_depth": md,
            }
            print(f"  {label} [{v}]: {len(nodes)} child-nodes, "
                  f"depth 0..{md}")
        if per_variant:
            all_traces[label] = per_variant

    for name, n in cases:
        try:
            d, e = make(name, n)
        except Exception as ex:
            print(f"  {name}@{n}: gen failed: {ex}")
            continue
        add(f"{name}@{n}", n, d, e)

    for basename in stcoll:
        path = os.path.join(ROOT, "stcollection", f"{basename}.dat")
        if not os.path.exists(path):
            print(f"  ST_{basename}: missing")
            continue
        sname, d, e = load_stcoll(path)
        n = len(d)
        add(f"ST_{basename}@{n}", n, d, e)

    with open(OUT_JSON, "w") as f:
        json.dump(all_traces, f, indent=1, default=float)
    print(f"\nwrote {OUT_JSON}  ({len(all_traces)} matrices x {len(variants_to_run)} variants)")


if __name__ == "__main__":
    main()
