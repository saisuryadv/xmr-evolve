#!/usr/bin/env python3
"""
Evaluate MR³-GK bidiagonal SVD implementation.

Usage:
  python3 evaluate.py                        # Full 379-test suite
  python3 evaluate.py --small                # Small 64-test suite
  python3 evaluate.py --medium               # 270-test suite
  python3 evaluate.py --large                # 90 patterns at n=800
  python3 evaluate.py --sweep100             # n=100,200,...,800
  python3 evaluate.py --ablation blockfactor # Ablation: block factorization on vs off
"""
import numpy as np
import sys
import os
import time
import glob
import re
import argparse
import subprocess

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
        t0 = time.perf_counter()
        sigma, U, V, info = bidiag_svd(d.copy(), e.copy())
        dt = time.perf_counter() - t0
    else:
        bidiag_svd(d.copy(), e.copy())
        bidiag_svd(d.copy(), e.copy())
        copies = [(d.copy(), e.copy()) for _ in range(5)]
        times = []
        for dc, ec in copies:
            t0 = time.perf_counter()
            sigma, U, V, info = bidiag_svd(dc, ec)
            times.append(time.perf_counter() - t0)
        times.sort()
        dt = times[2]
    if dt > timeout_s:
        return None, float('inf'), float('inf'), float('inf'), dt
    B = np.diag(d) + np.diag(e, 1)
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
    fe = __import__('full_eval')
    return fe.adv_names

SMALL_PATTERNS = [
    'exponential_graded', 'glued_repeated', 'saw_tooth', 'stemr_killer',
    'huge_condition', 'spike', 'wilkinson_like', 'two_clusters',
    'random_uniform', 'diagonal_only', 'constant', 'all_equal_nontrivial',
    'one_big_cluster', 'arithmetic_progression', 'many_near_zero',
    'random_dense_clusters', 'constant_d_graded_e', 'random_clustered_5',
    'alternating_sign', 'step_function', 'three_clusters', 'random_sparse_e',
]

# ======================================================================
#  Run a test suite, return per-test results
# ======================================================================
def run_suite(patterns, test_sizes, run_stcoll, quick_mode, timeout_s=30, silent=False):
    results = {}
    for sz in test_sizes:
        sz_pass = 0; sz_total = 0
        for pat in patterns:
            try:
                d, e = make(pat, sz)
            except:
                continue
            sz_total += 1
            key = f"{pat}@{sz}"
            try:
                ok, res, ou, ov, dt = test_one(d, e, timeout_s, quick=quick_mode)
                if ok is None:
                    results[key] = (False, float('inf'), float('inf'), float('inf'), dt)
                else:
                    results[key] = (ok, res, ou, ov, dt)
                    if ok: sz_pass += 1
            except:
                results[key] = (False, float('inf'), float('inf'), float('inf'), 0)
        if not silent:
            print(f"  n={sz:4d}: {sz_pass}/{sz_total} passed")

    if run_stcoll:
        stcoll_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stcollection')
        stcoll_files = sorted(glob.glob(os.path.join(stcoll_dir, 'B_*.dat')))
        st_pass = 0; st_total = 0
        for f in stcoll_files:
            name, d, e = load_stcoll(f)
            st_total += 1
            key = f"ST_{name}"
            try:
                ok, res, ou, ov, dt = test_one(d, e, timeout_s, quick=quick_mode)
                results[key] = (ok if ok is not None else False, res, ou, ov, dt)
                if ok: st_pass += 1
            except:
                results[key] = (False, float('inf'), float('inf'), float('inf'), 0)
        if not silent:
            print(f"  STColl: {st_pass}/{st_total} passed")

    return results

def print_scoring(results, patterns, test_sizes):
    total = len(results)
    passed = sum(1 for v in results.values() if v[0])
    pass_rate = passed / max(total, 1)
    pass_res = [v[1] for v in results.values() if v[0]]
    pass_ou  = [v[2] for v in results.values() if v[0]]
    pass_ov  = [v[3] for v in results.values() if v[0]]
    avg_res = np.mean(pass_res) if pass_res else 7.0
    avg_ou = np.mean(pass_ou) if pass_ou else 5.0
    avg_ov = np.mean(pass_ov) if pass_ov else 5.0

    pattern_data = {}
    for k, v in results.items():
        if '@' in k:
            pat, sz_s = k.rsplit('@', 1)
            sz = int(sz_s)
            pattern_data.setdefault(pat, {})[sz] = (v[0], v[4])

    scaling_pats = [p for p in patterns
                    if p in pattern_data and
                    all(sz in pattern_data[p] and pattern_data[p][sz][0] for sz in test_sizes)]

    worst_ratio = 0.0; worst_pattern = ""; all_ratios = []
    if len(test_sizes) >= 2:
        s_large = test_sizes[-1]; s_prev = test_sizes[-2]
        for pat in scaling_pats:
            pd = pattern_data[pat]
            t_prev = pd[s_prev][1]; t_large = pd[s_large][1]
            if t_prev > 1e-6:
                ratio = t_large / t_prev
                all_ratios.append((ratio, pat, t_prev, t_large))
                if ratio > worst_ratio:
                    worst_ratio = ratio; worst_pattern = pat
    all_ratios.sort(reverse=True)

    score = pass_rate * 50.0
    score += max(0.0, 5.0 - avg_res) * 2.0
    score += max(0.0, 5.0 - avg_ou) * 2.0
    score += max(0.0, 5.0 - avg_ov) * 2.0
    score += 5.0
    scaling_score = 10.0
    if worst_ratio > 5.0:
        scaling_score = 0.0; score = 5.0
    score += scaling_score

    failed_list = [k for k, v in sorted(results.items()) if not v[0]]

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
    for f in failed_list:
        print(f"  {f}")

# ======================================================================
#  Ablation: rebuild library with a feature toggled
# ======================================================================
def rebuild_lib(use_blocks=True):
    mydir = os.path.dirname(os.path.abspath(__file__))
    # Use a unique .so name to force ctypes to load the new library
    # (Linux won't reload a .so that's already mmap'd)
    tag = 'blk' if use_blocks else 'noblk'
    so_name = f'libxmr_{tag}.so'
    so_path = os.path.join(mydir, so_name)

    if not use_blocks:
        src = os.path.join(mydir, 'xmr_src', 'dlaxrs_stat.f')
        with open(src) as f:
            code = f.read()
        code = code.replace('USEBLOCKS = .TRUE.', 'USEBLOCKS = .FALSE.')
        tmpf = os.path.join(mydir, '_tmp_noblk.f')
        with open(tmpf, 'w') as f:
            f.write(code)
        subprocess.run([
            'gfortran', '-fPIC', '-O2', '-std=legacy', '-w',
            '-fno-second-underscore', '-c', tmpf,
            '-o', os.path.join(mydir, 'fortran_objects', 'dlaxrs_stat.o')
        ], check=True, capture_output=True)
        os.remove(tmpf)
    else:
        # Recompile with blocks from source
        subprocess.run([
            'gfortran', '-fPIC', '-O2', '-std=legacy', '-w',
            '-fno-second-underscore', '-c',
            os.path.join(mydir, 'xmr_src', 'dlaxrs_stat.f'),
            '-o', os.path.join(mydir, 'fortran_objects', 'dlaxrs_stat.o')
        ], check=True, capture_output=True)

    # Link to unique .so name
    objdir = os.path.join(mydir, 'fortran_objects')
    objs = [
        'dlaxrr', 'dlaxrs_stat', 'dlaxrs_prog', 'dlaxrs',
        'dlaxrt_stat', 'dlaxrt_prog', 'dlaxrt',
        'dlaxrn_stat', 'dlaxrn', 'dlaxrn0',
        'dlaxrc', 'dlaxrg', 'dlaxrx',
        'dlaxrb_clssfy', 'dlaxrb_refcls', 'dlaxrb_refsng',
        'dlaxrf_selshf', 'dlaxrf_seltw_part', 'dlaxrf_seltw',
        'dlaxrf_cob', 'dlaxrf_iib', 'dlaxrf_env', 'dlaxrf_grpenv',
        'dlaxrl_refine', 'dlaxrl_update', 'dlaxrl_reset', 'dlaxrf',
        'dlaxrm_stat2', 'dlaxrm_stat4', 'dlaxrm_stat8',
        'dlaxrm_stat16', 'dlaxrm_stat32', 'dlaxrm_stat64', 'dlaxrm',
        'dlaxre_initewldqds', 'dlaxre_gk', 'dlaxrv', 'xmr_wrapper',
    ]
    subprocess.run([
        'gcc', '-shared', '-o', so_path,
        *[os.path.join(objdir, f'{o}.o') for o in objs],
        '-lgfortran', '-llapack', '-lblas', '-lm'
    ], check=True, capture_output=True)

    # Point xmr_ctypes at the new .so via env var and reload
    os.environ['XMR_LIB_PATH'] = so_path
    import importlib, xmr_ctypes, mr3_gk
    importlib.reload(xmr_ctypes)
    importlib.reload(mr3_gk)
    global bidiag_svd
    from mr3_gk import bidiag_svd

def run_ablation(name, timeout_s=30):
    import subprocess, re

    mydir = os.path.dirname(os.path.abspath(__file__))

    if name != 'blockfactor':
        print(f"Unknown ablation: {name}. Available: blockfactor")
        sys.exit(1)

    def parse_results(output):
        results = {}
        for line in output.split('\n'):
            m = re.match(r'\s+(\S+)\s+n=\s*(\d+)\s+res=\s*(\S+)\s+ortU=\s*(\S+)\s+ortV=\s*(\S+)\s+t=\S+\s+(PASS|FAIL)', line)
            if m:
                results[f"{m.group(1)}@{m.group(2)}"] = {
                    'pass': m.group(6) == 'PASS',
                    'res': float(m.group(3)), 'ortU': float(m.group(4)),
                }
        return results

    print("=" * 70)
    print("ABLATION: Block Factorization (USEBLOCKS=TRUE vs FALSE)")
    print("=" * 70)

    # --- WITH blocks ---
    print("\n[1/3] Running full suite WITH block factorization...")
    subprocess.run(['bash', 'build.sh'], cwd=mydir, capture_output=True)
    r = subprocess.run([sys.executable, 'evaluate.py'],
                       cwd=mydir, capture_output=True, text=True, timeout=600)
    res_with = parse_results(r.stdout)
    n_with = sum(1 for v in res_with.values() if v['pass'])
    print(f"  {n_with}/{len(res_with)} passed")

    # --- WITHOUT blocks ---
    print("\n[2/3] Rebuilding WITHOUT blocks and running...")
    src = os.path.join(mydir, 'xmr_src', 'dlaxrs_stat.f')
    with open(src) as f:
        code = f.read().replace('USEBLOCKS = .TRUE.', 'USEBLOCKS = .FALSE.')
    with open('/tmp/_noblk.f', 'w') as f:
        f.write(code)
    subprocess.run(['gfortran', '-fPIC', '-O2', '-std=legacy', '-w',
                    '-fno-second-underscore', '-c', '/tmp/_noblk.f',
                    '-o', os.path.join(mydir, 'fortran_objects', 'dlaxrs_stat.o')],
                   check=True, capture_output=True)
    subprocess.run(['bash', '-c',
        'gcc -shared -o libxmr.so fortran_objects/*.o -lgfortran -llapack -lblas -lm'],
        cwd=mydir, capture_output=True)
    r = subprocess.run([sys.executable, 'evaluate.py'],
                       cwd=mydir, capture_output=True, text=True, timeout=600)
    res_without = parse_results(r.stdout)
    n_without = sum(1 for v in res_without.values() if v['pass'])
    print(f"  {n_without}/{len(res_without)} passed")

    # --- Restore ---
    print("\n[3/3] Restoring...")
    subprocess.run(['bash', 'build.sh'], cwd=mydir, capture_output=True)

    # --- Compare ---
    print(f"\n{'=' * 70}")
    print(f"  WITH blocks:    {n_with}/{len(res_with)}")
    print(f"  WITHOUT blocks: {n_without}/{len(res_without)}")

    helped = []
    for key in sorted(set(list(res_with.keys()) + list(res_without.keys()))):
        w = res_with.get(key, {"pass": False, "ortU": 999})
        wo = res_without.get(key, {"pass": False, "ortU": 999})
        if w['pass'] and not wo['pass']:
            helped.append((key, w, wo))

    if helped:
        print(f"\n  Block factorization HELPS ({len(helped)} tests):")
        print(f"    {'test':35s}  {'ortU(blk)':>10s}  {'ortU(noblk)':>12s}")
        for key, w, wo in helped:
            print(f"    {key:35s}  {w['ortU']:10.3f}  {wo['ortU']:12.3f}")
    else:
        print("\n  No difference.")

    print(f"\n  Summary: block factorization saves {len(helped)} tests")
    print(f"{'=' * 70}")

# ======================================================================
#  Main
# ======================================================================
def main():
    parser = argparse.ArgumentParser(description='Evaluate MR³-GK bidiagonal SVD')
    parser.add_argument('--small', action='store_true', help='Small 64-test suite (22 patterns × 3 sizes)')
    parser.add_argument('--medium', action='store_true', help='Medium 270-test suite (90 patterns × 3 sizes)')
    parser.add_argument('--large', action='store_true', help='Large suite: 90 patterns at n=800 only')
    parser.add_argument('--sweep100', action='store_true', help='Sweep n=100,200,...,800 (single run per test)')
    parser.add_argument('--ablation', type=str, metavar='NAME', help='Ablation study (e.g. blockfactor)')
    parser.add_argument('--timeout', type=int, default=30, help='Per-test timeout in seconds')
    args = parser.parse_args()

    if args.ablation:
        run_ablation(args.ablation, args.timeout)
        return

    if args.small:
        patterns = SMALL_PATTERNS
        test_sizes = [10, 50, 100]
        run_stcoll = False; quick_mode = False
        suite_name = "Small (64 tests)"
    elif args.medium:
        patterns = get_patterns()
        test_sizes = [10, 100, 200]
        run_stcoll = False; quick_mode = False
        suite_name = "Medium (270 tests)"
    elif args.large:
        patterns = get_patterns()
        test_sizes = [800]
        run_stcoll = False; quick_mode = True
        suite_name = "Large (90 patterns × n=800)"
    elif args.sweep100:
        patterns = get_patterns()
        test_sizes = [100, 200, 300, 400, 500, 600, 700, 800]
        run_stcoll = False; quick_mode = True
        suite_name = "Sweep n=100..800 step 100 (single run)"
    else:
        patterns = get_patterns()
        test_sizes = [10, 100, 200, 400]
        run_stcoll = True; quick_mode = False
        suite_name = "Full (379 tests)"

    print(f"=== {suite_name} ===\n")
    results = run_suite(patterns, test_sizes, run_stcoll, quick_mode, args.timeout)
    print_scoring(results, patterns, test_sizes)

if __name__ == '__main__':
    main()
