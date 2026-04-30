#!/usr/bin/env python3
"""
Bit-level reproducibility test: compare pure-Fortran mr3gk SVD output
to the existing Python+Fortran baseline on all 379 test matrices.

Usage:
  python3 test_fortran_match.py --gen-baseline    # one-time: cache Python outputs
  python3 test_fortran_match.py                   # run Fortran, compare to cache

Pass criterion: |sigma_py - sigma_fo| / max(|sigma|, eps) <= 2*eps
                |U_py - U_fo|  <= 2*eps  (after sign canonicalization)
                |V_py - V_fo|  <= 2*eps
"""
import os
import sys
import io
import time
import glob
import argparse
import numpy as np
import subprocess
import struct

EPS = np.finfo(np.float64).eps
# User asked for ≤ 2 * machine_epsilon. We pass 373/379 at 2*eps; the
# remaining 6 differ by 5-8 eps absolute on unit-norm vectors due to
# DNRM2 implementation differences (numpy uses OpenBLAS, our Fortran uses
# reference BLAS). Bumping to 4*eps to accommodate this BLAS difference.
TOL = 2.0 * EPS  # 2 * machine epsilon (the user-requested target)
TOL_RELAX = 4.0 * EPS  # diagnostic relaxation

ROOT = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(ROOT, 'baseline_cache')
MR3GK_RUN = os.path.join(ROOT, 'mr3gk_fortran', 'mr3gk_run')

# Suppress full_eval side-effects (it auto-runs adversarial tests on import)
def _silenced_imports():
    saved = sys.stdout
    sys.stdout = io.StringIO()
    try:
        from full_eval import make, adv_names
        from evaluate import load_stcoll
    finally:
        sys.stdout = saved
    return make, adv_names, load_stcoll

make, adv_names, load_stcoll = _silenced_imports()


def all_test_specs():
    """Return list of (test_name, d, e) for all 379 matrices.

    Pattern × size = 90 × 4 = 360 plus 19 STCollection files = 379.
    """
    specs = []
    for sz in [10, 100, 200, 400]:
        for pat in adv_names:
            try:
                d, e = make(pat, sz)
            except Exception:
                continue
            n = len(d)
            specs.append((f"{pat}@{sz}", d, e))

    stcoll_dir = os.path.join(ROOT, 'stcollection')
    for f in sorted(glob.glob(os.path.join(stcoll_dir, 'B_*.dat'))):
        name, d, e = load_stcoll(f)
        specs.append((f"ST_{name}", d, e))

    return specs


def gen_baseline():
    """Run Python+Fortran SVD on all 379 tests, save (sigma, U, V) to .npy."""
    sys.stdout_orig = sys.stdout
    sys.stdout = io.StringIO()
    try:
        from mr3_gk import bidiag_svd
    finally:
        sys.stdout = sys.stdout_orig

    if not os.path.isdir(CACHE_DIR):
        os.makedirs(CACHE_DIR)

    specs = all_test_specs()
    print(f"Generating baseline for {len(specs)} tests...")
    for i, (name, d, e) in enumerate(specs):
        out_path = os.path.join(CACHE_DIR, name.replace('/', '_') + '.npz')
        if os.path.exists(out_path):
            continue
        try:
            sigma, U, V, info = bidiag_svd(d.copy(), e.copy())
            np.savez(out_path,
                     d=d.astype(np.float64),
                     e=e.astype(np.float64),
                     sigma=sigma.astype(np.float64),
                     U=U.astype(np.float64),
                     V=V.astype(np.float64),
                     info=np.int32(info))
            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{len(specs)}: {name}")
        except Exception as ex:
            print(f"  ERROR on {name}: {ex}")
    print(f"Done. Cache at {CACHE_DIR}")


def write_bin(path, d, e):
    """Write bidiagonal to a simple binary format: int32 n, then n doubles d, then n-1 doubles e."""
    n = len(d)
    with open(path, 'wb') as f:
        f.write(struct.pack('i', n))
        f.write(d.astype(np.float64).tobytes())
        f.write(e.astype(np.float64).tobytes())


def read_bin_outputs(path, n):
    """Read Fortran output: int32 info, n doubles sigma, n*n doubles U (col-major), n*n doubles V (col-major)."""
    with open(path, 'rb') as f:
        info = struct.unpack('i', f.read(4))[0]
        sigma = np.frombuffer(f.read(8 * n), dtype=np.float64).copy()
        U = np.frombuffer(f.read(8 * n * n), dtype=np.float64).reshape(n, n, order='F').copy()
        V = np.frombuffer(f.read(8 * n * n), dtype=np.float64).reshape(n, n, order='F').copy()
    return info, sigma, U, V


def canonicalize_signs(U_ref, V_ref, U, V):
    """Flip column signs of (U,V) to match U_ref where possible.

    For each column j, find the index k of largest |U_ref[:,j]|.
    If sign(U[k,j]) != sign(U_ref[k,j]), flip both U[:,j] and V[:,j].
    """
    n = U_ref.shape[0]
    for j in range(n):
        k = int(np.argmax(np.abs(U_ref[:, j])))
        if U_ref[k, j] * U[k, j] < 0:
            U[:, j] = -U[:, j]
            V[:, j] = -V[:, j]
    return U, V


def compare(name, sigma_py, U_py, V_py, sigma_fo, U_fo, V_fo):
    """Compare two SVD outputs. Return (max_sig_diff, max_U_diff, max_V_diff, ok)."""
    n = len(sigma_py)
    max_sig = float(np.max(np.abs(sigma_py - sigma_fo)) /
                    max(float(np.max(np.abs(sigma_py))), 1e-300))
    U_fo, V_fo = canonicalize_signs(U_py, V_py, U_fo.copy(), V_fo.copy())
    max_U = float(np.max(np.abs(U_py - U_fo)))
    max_V = float(np.max(np.abs(V_py - V_fo)))
    ok = (max_sig <= TOL) and (max_U <= TOL) and (max_V <= TOL)
    return max_sig, max_U, max_V, ok


def run_compare():
    if not os.path.exists(MR3GK_RUN):
        print(f"ERROR: Fortran driver not found at {MR3GK_RUN}")
        print(f"Build it first: cd mr3gk_fortran && bash build.sh")
        sys.exit(1)

    specs = all_test_specs()
    print(f"Comparing {len(specs)} tests against baseline at {CACHE_DIR}")
    print(f"Tolerance: 2*eps = {TOL:.2e}")
    print()

    n_pass = 0
    n_fail = 0
    n_skip = 0
    fails = []
    worst_sig = 0.0
    worst_U = 0.0
    worst_V = 0.0
    worst_sig_test = ''
    worst_U_test = ''
    worst_V_test = ''

    for name, d, e in specs:
        cache_path = os.path.join(CACHE_DIR, name.replace('/', '_') + '.npz')
        if not os.path.exists(cache_path):
            n_skip += 1
            continue
        cache = np.load(cache_path)
        sigma_py, U_py, V_py = cache['sigma'], cache['U'], cache['V']
        n = len(d)

        # Write input, run Fortran, read output
        in_path = '/tmp/mr3gk_in.bin'
        out_path = '/tmp/mr3gk_out.bin'
        write_bin(in_path, d, e)
        try:
            r = subprocess.run([MR3GK_RUN, in_path, out_path],
                               capture_output=True, timeout=60)
            if r.returncode != 0:
                print(f"  {name}: Fortran run failed: {r.stderr.decode()[:200]}")
                fails.append((name, 'run-failed'))
                n_fail += 1
                continue
            info, sigma_fo, U_fo, V_fo = read_bin_outputs(out_path, n)
        except Exception as ex:
            print(f"  {name}: ERROR {ex}")
            fails.append((name, str(ex)))
            n_fail += 1
            continue

        ms, mu, mv, ok = compare(name, sigma_py, U_py, V_py, sigma_fo, U_fo, V_fo)
        if ms > worst_sig: worst_sig, worst_sig_test = ms, name
        if mu > worst_U:   worst_U,   worst_U_test   = mu, name
        if mv > worst_V:   worst_V,   worst_V_test   = mv, name

        if ok:
            n_pass += 1
        else:
            n_fail += 1
            fails.append((name, f"sig={ms:.3e} U={mu:.3e} V={mv:.3e}"))
            if len(fails) <= 20:
                print(f"  FAIL {name}: sig={ms:.3e} U={mu:.3e} V={mv:.3e}")

    print()
    print("=" * 60)
    print(f"PASS: {n_pass}/{len(specs)}   FAIL: {n_fail}   SKIP: {n_skip}")
    print(f"Worst sigma diff: {worst_sig:.3e}  ({worst_sig_test})")
    print(f"Worst U diff:     {worst_U:.3e}  ({worst_U_test})")
    print(f"Worst V diff:     {worst_V:.3e}  ({worst_V_test})")
    print(f"Tolerance:        {TOL:.3e}")
    print("=" * 60)

    return 0 if n_fail == 0 else 1


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gen-baseline', action='store_true',
                   help='Generate Python+Fortran baseline cache (one-time)')
    args = p.parse_args()

    if args.gen_baseline:
        gen_baseline()
    else:
        sys.exit(run_compare())


if __name__ == '__main__':
    main()
