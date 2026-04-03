#!/usr/bin/env python3
"""
Evaluate MR³-GK bidiagonal SVD implementation.

Usage:
  python3 evaluate.py              # Full 379-test suite (90 patterns × 4 sizes + 19 STCollection)
  python3 evaluate.py --small      # Small 64-test suite (22 patterns × 3 sizes, no STCollection)
  python3 evaluate.py --medium     # 270-test suite (90 patterns × 3 sizes, no STCollection)
  python3 evaluate.py --large      # Large suite: 90 patterns at n=800 only
  python3 evaluate.py --sweep100   # Sweep n=100,200,...,800 (single run per test)
"""
import numpy as np
import sys
import os
import time
import glob
import re
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mr3_gk import bidiag_svd
from full_eval import make

EPS = 2.2204460492503131e-16
RES_THRESH = 7.0
ORTHO_THRESH = 5.0

def load_stcoll(path):
    with open(path) as f:
        n = int(f.readline().strip())
        d = np.zeros(n); e = np.zeros(max(n-1, 0))
        for i in range(n):
            parts = f.readline().split()
            d[i] = float(parts[1])
            if i < n-1 and len(parts) > 2:
                e[i] = float(parts[2])
    return os.path.basename(path).replace('.dat', ''), d, e

def test_one(d, e, timeout_s=30, quick=False):
    n = len(d)
    if quick:
        # Single run, no warmup (for large n)
        t0 = time.perf_counter()
        sigma, U, V, info = bidiag_svd(d.copy(), e.copy())
        dt = time.perf_counter() - t0
    else:
        # Warmup runs (JIT, cache effects)
        bidiag_svd(d.copy(), e.copy())
        bidiag_svd(d.copy(), e.copy())
        # Pre-copy arrays for timed runs
        copies = [(d.copy(), e.copy()) for _ in range(5)]
        # 5 timed runs, take median
        times = []
        for dc, ec in copies:
            t0 = time.perf_counter()
            sigma, U, V, info = bidiag_svd(dc, ec)
            times.append(time.perf_counter() - t0)
        times.sort()
        dt = times[2]  # median
    if dt > timeout_s:
        return None, float('inf'), float('inf'), float('inf'), dt
    B = np.diag(d) + np.diag(e, 1)
    # Match C++ compute_bnorm: infinity norm of bidiagonal (max row sum)
    bnorm = abs(d[-1]) if len(d) > 0 else 0.0
    for i in range(len(d) - 1):
        bnorm = max(bnorm, abs(d[i]) + abs(e[i]))
    if bnorm == 0.0:
        bnorm = 1.0
    res = np.max(np.abs(B - U @ np.diag(sigma) @ V.T)) / (bnorm * n * EPS)
    ou = np.max(np.abs(U.T @ U - np.eye(n))) / (n * EPS)
    ov = np.max(np.abs(V.T @ V - np.eye(n))) / (n * EPS)
    ok = res <= RES_THRESH and ou <= ORTHO_THRESH and ov <= ORTHO_THRESH
    return ok, res, ou, ov, dt

def get_patterns():
    """Get all pattern names from full_eval.py's adv_names list"""
    import importlib.util
    spec = importlib.util.spec_from_file_location("full_eval",
        os.path.join(os.path.dirname(os.path.abspath(__file__)), 'full_eval.py'))
    fe = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fe)
    return fe.adv_names

SMALL_PATTERNS = [
    'exponential_graded', 'glued_repeated', 'saw_tooth', 'stemr_killer',
    'huge_condition', 'spike', 'wilkinson_like', 'two_clusters',
    'random_uniform', 'diagonal_only', 'constant', 'all_equal_nontrivial',
    'one_big_cluster', 'arithmetic_progression', 'many_near_zero',
    'random_dense_clusters', 'constant_d_graded_e', 'random_clustered_5',
    'alternating_sign', 'step_function', 'three_clusters', 'random_sparse_e',
]

def main():
    parser = argparse.ArgumentParser(description='Evaluate MR³-GK bidiagonal SVD')
    parser.add_argument('--small', action='store_true', help='Small 64-test suite (22 patterns × 3 sizes)')
    parser.add_argument('--medium', action='store_true', help='Medium 270-test suite (90 patterns × 3 sizes)')
    parser.add_argument('--large', action='store_true', help='Large suite: 90 patterns at n=800 only')
    parser.add_argument('--sweep100', action='store_true', help='Sweep n=100,200,...,800 (single run per test)')
    parser.add_argument('--timeout', type=int, default=30, help='Per-test timeout in seconds')
    args = parser.parse_args()

    if args.small:
        patterns = SMALL_PATTERNS
        test_sizes = [10, 50, 100]
        run_stcoll = False
        quick_mode = False
        suite_name = "Small (64 tests)"
    elif args.medium:
        patterns = get_patterns()
        test_sizes = [10, 100, 200]
        run_stcoll = False
        quick_mode = False
        suite_name = "Medium (270 tests)"
    elif args.large:
        patterns = get_patterns()
        test_sizes = [800]
        run_stcoll = False
        quick_mode = True
        suite_name = "Large (90 patterns × n=800)"
    elif args.sweep100:
        patterns = get_patterns()
        test_sizes = [100, 200, 300, 400, 500, 600, 700, 800]
        run_stcoll = False
        quick_mode = True
        suite_name = "Sweep n=100..800 step 100 (single run)"
    else:
        patterns = get_patterns()
        test_sizes = [10, 100, 200, 400]
        run_stcoll = True
        quick_mode = False
        suite_name = "Full (379 tests)"

    print(f"=== {suite_name} ===\n")

    total = 0; passed = 0; failed_list = []
    pass_res = []; pass_ou = []; pass_ov = []
    pattern_data = {}

    # Adversarial tests
    for sz in test_sizes:
        sz_pass = 0; sz_total = 0
        for pat in patterns:
            try:
                d, e = make(pat, sz)
            except:
                continue
            n = len(d)
            sz_total += 1; total += 1
            try:
                ok, res, ou, ov, dt = test_one(d, e, args.timeout, quick=quick_mode)
                if ok is None:
                    failed_list.append(f"{pat}@{sz} TIMEOUT({dt:.0f}s)")
                elif ok:
                    sz_pass += 1; passed += 1
                    pass_res.append(res); pass_ou.append(ou); pass_ov.append(ov)
                else:
                    failed_list.append(f"{pat}@{sz}")
                pattern_data.setdefault(pat, {})[sz] = (ok if ok is not None else False, dt)
            except Exception as ex:
                failed_list.append(f"{pat}@{sz} ERR:{ex}")
                pattern_data.setdefault(pat, {})[sz] = (False, 0)
        print(f"  n={sz:4d}: {sz_pass}/{sz_total} passed")

    # STCollection
    st_pass = 0; st_total = 0
    if run_stcoll:
        stcoll_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stcollection')
        stcoll_files = sorted(glob.glob(os.path.join(stcoll_dir, 'B_*.dat')))
        for f in stcoll_files:
            name, d, e = load_stcoll(f)
            st_total += 1; total += 1
            try:
                ok, res, ou, ov, dt = test_one(d, e, args.timeout)
                if ok:
                    st_pass += 1; passed += 1
                    pass_res.append(res); pass_ou.append(ou); pass_ov.append(ov)
                else:
                    failed_list.append(f"ST_{name}")
            except:
                failed_list.append(f"ST_{name} ERR")
        print(f"  STColl: {st_pass}/{st_total} passed")

    # Scoring
    pass_rate = passed / max(total, 1)
    avg_res = np.mean(pass_res) if pass_res else 7.0
    avg_ou = np.mean(pass_ou) if pass_ou else 5.0
    avg_ov = np.mean(pass_ov) if pass_ov else 5.0

    # Scaling: patterns passing at ALL sizes
    scaling_pats = []
    for pat in patterns:
        if pat in pattern_data:
            pd = pattern_data[pat]
            if all(sz in pd and pd[sz][0] for sz in test_sizes):
                scaling_pats.append(pat)

    worst_ratio = 0.0
    worst_pattern = ""
    all_ratios = []
    if len(test_sizes) >= 2:
        s_large = test_sizes[-1]
        s_prev = test_sizes[-2]
        for pat in scaling_pats:
            pd = pattern_data[pat]
            t_prev = pd[s_prev][1]
            t_large = pd[s_large][1]
            if t_prev > 1e-6:
                ratio = t_large / t_prev
                all_ratios.append((ratio, pat, t_prev, t_large))
                if ratio > worst_ratio:
                    worst_ratio = ratio
                    worst_pattern = pat
    all_ratios.sort(reverse=True)

    score = pass_rate * 50.0
    score += max(0.0, 5.0 - avg_res) * 2.0
    score += max(0.0, 5.0 - avg_ou) * 2.0
    score += max(0.0, 5.0 - avg_ov) * 2.0
    score += 5.0  # compilation bonus
    scaling_score = 10.0
    if worst_ratio > 5.0:
        scaling_score = 0.0
        score = 5.0  # HARD GATE
    score += scaling_score

    print(f"\n{'='*60}")
    print(f"TOTAL: {passed}/{total} passed ({100*pass_rate:.1f}%)")
    print(f"Pass avg: res={avg_res:.3f}  oU={avg_ou:.3f}  oV={avg_ov:.3f}")
    print(f"Scaling: {len(scaling_pats)} patterns, worst_ratio={worst_ratio:.3f} ({worst_pattern})")
    if all_ratios:
        print(f"  Top-10 worst scaling ratios (n={test_sizes[-2]}→{test_sizes[-1]}):")
        for ratio, pat, tp, tl in all_ratios[:10]:
            print(f"    {ratio:6.2f}x  {pat:35s}  t{test_sizes[-2]}={tp:.4f}  t{test_sizes[-1]}={tl:.4f}")
    print(f"SCORE: {score:.2f}")
    if worst_ratio > 5.0:
        print(f"  *** HARD GATE: worst_ratio={worst_ratio:.1f} > 5.0, score CAPPED at 5 ***")
    print(f"\nFailing ({len(failed_list)}):")
    for f in sorted(failed_list):
        print(f"  {f}")

if __name__ == '__main__':
    main()
