#!/usr/bin/env python3
"""
Evaluate MR³-GK bidiagonal SVD implementation.

Usage:
  python3 evaluate.py              # Full 379-test suite (90 patterns × 4 sizes + 19 STCollection)
  python3 evaluate.py --small      # Small 64-test suite (22 patterns × 3 sizes, no STCollection)
  python3 evaluate.py --medium     # 270-test suite (90 patterns × 3 sizes, no STCollection)
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

def test_one(d, e, timeout_s=30):
    n = len(d)
    t0 = time.time()
    sigma, U, V, info = bidiag_svd(d, e)
    dt = time.time() - t0
    if dt > timeout_s:
        return None, float('inf'), float('inf'), float('inf'), dt
    B = np.diag(d) + np.diag(e, 1)
    bn = max(np.max(np.abs(B)), 1e-300)
    res = np.max(np.abs(B - U @ np.diag(sigma) @ V.T)) / (bn * n * EPS)
    ou = np.max(np.abs(U.T @ U - np.eye(n))) / (n * EPS)
    ov = np.max(np.abs(V.T @ V - np.eye(n))) / (n * EPS)
    ok = res <= RES_THRESH and ou <= ORTHO_THRESH and ov <= ORTHO_THRESH
    return ok, res, ou, ov, dt

def get_patterns():
    """Extract all pattern names from full_eval.py"""
    src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'full_eval.py')).read()
    return list(dict.fromkeys(re.findall(r"(?:if|elif)\s+name\s*==\s*'(\w+)'", src)))

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
    parser.add_argument('--timeout', type=int, default=30, help='Per-test timeout in seconds')
    args = parser.parse_args()

    if args.small:
        patterns = SMALL_PATTERNS
        test_sizes = [10, 50, 100]
        run_stcoll = False
        suite_name = "Small (64 tests)"
    elif args.medium:
        patterns = get_patterns()
        test_sizes = [10, 100, 200]
        run_stcoll = False
        suite_name = "Medium (270 tests)"
    else:
        patterns = get_patterns()
        test_sizes = [10, 100, 200, 400]
        run_stcoll = True
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
                ok, res, ou, ov, dt = test_one(d, e, args.timeout)
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
    if len(test_sizes) >= 2:
        s_large = test_sizes[-1]
        s_prev = test_sizes[-2]
        for pat in scaling_pats:
            pd = pattern_data[pat]
            t_prev = pd[s_prev][1]
            t_large = pd[s_large][1]
            if t_prev > 1e-6:
                worst_ratio = max(worst_ratio, t_large / t_prev)

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
    print(f"Scaling: {len(scaling_pats)} patterns, worst_ratio={worst_ratio:.3f}")
    print(f"SCORE: {score:.2f}")
    if worst_ratio > 5.0:
        print(f"  *** HARD GATE: worst_ratio={worst_ratio:.1f} > 5.0, score CAPPED at 5 ***")
    print(f"\nFailing ({len(failed_list)}):")
    for f in sorted(failed_list):
        print(f"  {f}")

if __name__ == '__main__':
    main()
