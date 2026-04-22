#!/usr/bin/env python3
"""Debug the Godunov matrix failure in MR3-GK bidiagonal SVD.

The Godunov 1e-5 matrix (n=2500) has ortU=30,628. It has:
  - Zero diagonal (d_i = 0 for all i)
  - Alternating off-diagonal: e_i = 900 (odd), e_i = 1e-5 (even)

We start with smaller Godunov matrices (n=73, 113, 147, 169) and use
existing verbose mode to trace where MR3 loses orthogonality.
"""
import sys
import os
import argparse
import io

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Suppress the full_eval side-effect print block that happens on import
_saved_stdout = sys.stdout
sys.stdout = io.StringIO()
from test_dense_to_bidiag import tridiag_to_bidiag, test_one
sys.stdout = _saved_stdout

from mr3_gk import bidiag_svd, set_verbose
import numpy as np

EPS = np.finfo(np.float64).eps


def load_tridiag(path):
    with open(path) as f:
        n = int(f.readline().strip())
        d = np.zeros(n); e = np.zeros(max(n - 1, 0))
        for i in range(n):
            parts = f.readline().split()
            d[i] = float(parts[1])
            if i < n - 1 and len(parts) > 2:
                e[i] = float(parts[2])
    return d, e


def analyze_result(d_bidiag, e_bidiag, sigma, U, V, label):
    n = len(d_bidiag)
    B = np.diag(d_bidiag) + np.diag(e_bidiag, 1)
    bnorm = max(np.max(np.abs(d_bidiag) + np.append(np.abs(e_bidiag), 0)), 1e-300)

    res = np.max(np.abs(B - U @ np.diag(sigma) @ V.T)) / (bnorm * n * EPS)
    UU = U.T @ U - np.eye(n)
    VV = V.T @ V - np.eye(n)
    ou = np.max(np.abs(UU)) / (n * EPS)
    ov = np.max(np.abs(VV)) / (n * EPS)

    print(f"\n=== {label} ===")
    print(f"  n={n}  res={res:.3f}  ortU={ou:.3f}  ortV={ov:.3f}  "
          f"{'PASS' if (res<=7 and ou<=5 and ov<=5) else 'FAIL'}")

    # Worst pairs
    np.fill_diagonal(UU, 0)
    worst_val = np.max(np.abs(UU))
    wi, wj = np.unravel_index(np.argmax(np.abs(UU)), UU.shape)
    print(f"  Worst U pair: ({wi},{wj})  |u_i^T u_j|={worst_val:.3e} = {worst_val/(n*EPS):.1f} neps")
    print(f"    sigma[{wi}] = {sigma[wi]:.10e}")
    print(f"    sigma[{wj}] = {sigma[wj]:.10e}")
    relgap = abs(sigma[wi] - sigma[wj]) / max(abs(sigma[wi]), 1e-300)
    print(f"    relative gap = {relgap:.3e}")

    # Top 10 worst pairs (upper triangle)
    print(f"\n  Top 10 worst U orthogonality pairs:")
    flat_idx = np.argsort(-np.abs(UU).flatten())[:50]
    shown = 0
    seen_pairs = set()
    for idx in flat_idx:
        i, j = divmod(idx, n)
        if i >= j: continue
        if (i,j) in seen_pairs: continue
        seen_pairs.add((i,j))
        val = UU[i,j]
        rg = abs(sigma[i] - sigma[j]) / max(abs(sigma[i]), 1e-300)
        print(f"    ({i:4d},{j:4d})  {abs(val)/(n*EPS):8.2f} neps  "
              f"sigma_i={sigma[i]:.6e}  sigma_j={sigma[j]:.6e}  relgap={rg:.3e}")
        shown += 1
        if shown >= 10: break

    # Histogram of sigma values
    print(f"\n  Singular value distribution (sorted desc):")
    s_sorted = np.sort(np.abs(sigma))[::-1]
    print(f"    max = {s_sorted[0]:.6e}")
    print(f"    min = {s_sorted[-1]:.6e}")
    if s_sorted[-1] > 0:
        print(f"    cond = {s_sorted[0]/s_sorted[-1]:.3e}")

    # Count tight pairs (relative gap < machine eps * 100)
    gaps = np.abs(np.diff(s_sorted))
    tight_mask = gaps < np.abs(s_sorted[:-1]) * EPS * 100
    print(f"    tight pairs (relgap < 100*eps): {np.sum(tight_mask)} of {n-1}")

    return res, ou, ov


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--size', default='73',
                        help='Godunov matrix size: 73, 113, 147, 169, or 2500 '
                             'for 1e-2/1e-4/1e-5/1e-6/1e-7 (specify like 1e-5)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose MR3-GK tracing')
    parser.add_argument('--max-verbose-n', type=int, default=150,
                        help='Max n for verbose (above this, auto-disable)')
    args = parser.parse_args()

    # Resolve path
    if args.size.startswith('1e'):
        path = f'/tmp/STCollection/DATA/T_Godunov_{args.size}.dat'
    else:
        n = int(args.size)
        path = f'/tmp/STCollection/DATA/T_Godunov_{n:03d}.dat'

    if not os.path.exists(path):
        print(f"ERROR: {path} not found")
        print(f"Available: " + ', '.join(sorted(
            os.path.basename(f).replace('T_Godunov_', '').replace('.dat', '')
            for f in __import__('glob').glob('/tmp/STCollection/DATA/T_Godunov_*.dat'))))
        sys.exit(1)

    print(f"Loading {path}...")
    d_tri, e_tri = load_tridiag(path)
    n = len(d_tri)
    print(f"  Tridiagonal: n={n}")
    print(f"  d (first 10): {d_tri[:10]}")
    print(f"  e (first 10): {e_tri[:10]}")
    print(f"  d range: [{np.min(d_tri):.6e}, {np.max(d_tri):.6e}]")
    print(f"  e range: [{np.min(e_tri):.6e}, {np.max(e_tri):.6e}]")
    print(f"  zero-diagonal: {np.all(d_tri == 0)}")

    # Cholesky to bidiagonal
    print(f"\nConverting tridiagonal → bidiagonal via Cholesky (T - nu*I = L L^T)...")
    d_b, e_b = tridiag_to_bidiag(d_tri, e_tri)
    print(f"  Bidiagonal: n={len(d_b)}")
    print(f"  d_bidiag range: [{np.min(np.abs(d_b)):.6e}, {np.max(np.abs(d_b)):.6e}]")
    if len(e_b) > 0:
        print(f"  e_bidiag range: [{np.min(np.abs(e_b)):.6e}, {np.max(np.abs(e_b)):.6e}]")
    print(f"  d_bidiag (first 20): {d_b[:20]}")
    print(f"  e_bidiag (first 20): {e_b[:20]}")

    # Set verbose if requested and matrix is small enough
    use_verbose = args.verbose and n <= args.max_verbose_n
    if args.verbose and not use_verbose:
        print(f"\n  NOTE: verbose disabled because n={n} > {args.max_verbose_n}")
        print(f"        (too much output). Use --max-verbose-n to override.")

    if use_verbose:
        print(f"\n{'='*60}")
        print(f"VERBOSE MR3-GK TRACE (n={n})")
        print(f"{'='*60}")
        set_verbose(True)

    # Run bidiag_svd
    sigma, U, V, info = bidiag_svd(d_b, e_b)

    if use_verbose:
        set_verbose(False)

    print(f"\nbidiag_svd returned info={info}")

    # Analysis
    analyze_result(d_b, e_b, sigma, U, V, f"Godunov {args.size} Analysis")

    # Compare to numpy SVD (gold standard)
    print(f"\n=== Reference: numpy.linalg.svd ===")
    B = np.diag(d_b) + np.diag(e_b, 1)
    U_ref, sigma_ref, Vt_ref = np.linalg.svd(B)
    V_ref = Vt_ref.T
    analyze_result(d_b, e_b, sigma_ref, U_ref, V_ref, "numpy SVD (reference)")


if __name__ == '__main__':
    main()
