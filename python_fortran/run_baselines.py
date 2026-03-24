#!/usr/bin/env python3
"""
Baseline comparison: DBDSQR, TGK+DSTEMR vs MR³-GK.
All evaluated with same methodology: 2 warmup + median of 5, perf_counter, pre-copied arrays.

Usage:
  python3 run_baselines.py
"""
import numpy as np, sys, os, time, re, glob, ctypes, ctypes.util

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from full_eval import make

EPS = 2.2204460492503131e-16
RES_THRESH = 7.0
ORTHO_THRESH = 5.0

# ============================================================
# LAPACK baselines via ctypes
# ============================================================
_lapack = ctypes.CDLL(ctypes.util.find_library('lapack'))

def _ptr(arr):
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
def _iptr(arr):
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
def _int(v):
    return ctypes.byref(ctypes.c_int(v))
def _char(c):
    return ctypes.c_char_p(c.encode())

def bidiag_svd_dbdsqr(d, e):
    """DBDSQR: QR iteration (O(n^3) worst case, gold standard accuracy)."""
    d = np.asarray(d, dtype=np.float64).copy()
    e = np.asarray(e, dtype=np.float64).copy()
    n = len(d)
    if n == 0: return np.array([]), np.zeros((0,0)), np.zeros((0,0)), 0
    if n == 1:
        return (np.array([abs(d[0])]),
                np.array([[1.0 if d[0]>=0 else -1.0]]),
                np.array([[1.0]]), 0)
    sigma = d.copy()
    ee = np.zeros(n); ee[:n-1] = e.copy()
    VT = np.eye(n, dtype=np.float64, order='F')
    U = np.eye(n, dtype=np.float64, order='F')
    C = np.zeros((1,1), dtype=np.float64, order='F')
    work = np.zeros(4*n, dtype=np.float64)
    info = ctypes.c_int(0)
    _lapack.dbdsqr_(_char('U'), _int(n), _int(n), _int(n), _int(0),
                     _ptr(sigma), _ptr(ee), _ptr(VT), _int(n),
                     _ptr(U), _int(n), _ptr(C), _int(1),
                     _ptr(work), ctypes.byref(info))
    return np.abs(sigma), U, VT.T, info.value

# ============================================================
# Evaluation framework
# ============================================================
def get_patterns():
    src = open(os.path.join(os.path.dirname(__file__), 'full_eval.py')).read()
    return list(dict.fromkeys(re.findall(r"(?:if|elif)\s+name\s*==\s*'(\w+)'", src)))

def load_stcoll(path):
    with open(path) as f:
        n = int(f.readline().strip())
        d = np.zeros(n); e = np.zeros(max(n-1, 0))
        for i in range(n):
            parts = f.readline().split()
            d[i] = float(parts[1])
            if i < n-1 and len(parts) > 2: e[i] = float(parts[2])
    return os.path.basename(path).replace('.dat', ''), d, e

def compute_bnorm(d, e):
    """Infinity norm of bidiagonal matrix, matching C++ evaluate.cpp."""
    bnorm = abs(d[-1]) if len(d) > 0 else 0.0
    for i in range(len(d) - 1):
        bnorm = max(bnorm, abs(d[i]) + abs(e[i]))
    return bnorm if bnorm != 0.0 else 1.0

def evaluate_algo(name, bidiag_svd_fn):
    patterns = get_patterns()
    test_sizes = [10, 100, 200, 400]

    total = 0; passed = 0
    pass_res = []; pass_ou = []; pass_ov = []
    pattern_data = {}

    for sz in test_sizes:
        for pat in patterns:
            try: d, e = make(pat, sz)
            except: continue
            n = len(d); total += 1
            try:
                # 2 warmup runs
                bidiag_svd_fn(d.copy(), e.copy())
                bidiag_svd_fn(d.copy(), e.copy())
                # 5 timed runs with pre-copied arrays, take median
                copies = [(d.copy(), e.copy()) for _ in range(5)]
                times = []
                for dc, ec in copies:
                    t0 = time.perf_counter()
                    sigma, U, V, info = bidiag_svd_fn(dc, ec)
                    times.append(time.perf_counter() - t0)
                times.sort()
                dt = times[2]

                B = np.diag(d) + np.diag(e, 1)
                bnorm = compute_bnorm(d, e)
                res = np.max(np.abs(B - U@np.diag(sigma)@V.T)) / (bnorm*n*EPS)
                ou = np.max(np.abs(U.T@U - np.eye(n))) / (n*EPS)
                ov = np.max(np.abs(V.T@V - np.eye(n))) / (n*EPS)
                ok = res<=RES_THRESH and ou<=ORTHO_THRESH and ov<=ORTHO_THRESH
                if ok:
                    passed += 1
                    pass_res.append(res); pass_ou.append(ou); pass_ov.append(ov)
                pattern_data.setdefault(pat, {})[sz] = (ok, dt)
            except:
                pattern_data.setdefault(pat, {})[sz] = (False, 0)

    # STCollection
    st_pass = 0; st_total = 0
    stdir = os.path.join(os.path.dirname(__file__), 'stcollection')
    if os.path.isdir(stdir):
        for f in sorted(glob.glob(os.path.join(stdir, '*.dat'))):
            _, d, e = load_stcoll(f); n = len(d); st_total += 1; total += 1
            try:
                sigma, U, V, info = bidiag_svd_fn(d.copy(), e.copy())
                B = np.diag(d) + np.diag(e, 1)
                bnorm = compute_bnorm(d, e)
                res = np.max(np.abs(B - U@np.diag(sigma)@V.T)) / (bnorm*n*EPS)
                ou = np.max(np.abs(U.T@U - np.eye(n))) / (n*EPS)
                ov = np.max(np.abs(V.T@V - np.eye(n))) / (n*EPS)
                if res<=RES_THRESH and ou<=ORTHO_THRESH and ov<=ORTHO_THRESH:
                    st_pass += 1; passed += 1
                    pass_res.append(res); pass_ou.append(ou); pass_ov.append(ov)
            except: pass

    # Scaling
    scaling_pats = [pat for pat in patterns
                    if pat in pattern_data and
                    all(sz in pattern_data[pat] and pattern_data[pat][sz][0] for sz in test_sizes)]
    worst_ratio = 0.0; worst_pat = ""
    ratios = []
    for pat in scaling_pats:
        t200 = pattern_data[pat][200][1]
        t400 = pattern_data[pat][400][1]
        if t200 > 1e-6:
            r = t400/t200
            ratios.append((r, pat))
            if r > worst_ratio:
                worst_ratio = r; worst_pat = pat
    ratios.sort(reverse=True)

    avg_res = np.mean(pass_res) if pass_res else 999
    avg_ou = np.mean(pass_ou) if pass_ou else 999
    avg_ov = np.mean(pass_ov) if pass_ov else 999

    score = (passed/max(total,1))*50.0
    score += max(0, 5-avg_res)*2 + max(0, 5-avg_ou)*2 + max(0, 5-avg_ov)*2
    score += 5.0
    if worst_ratio <= 5.0:
        score += 10.0
    else:
        score = 5.0

    print(f"\n{'='*60}")
    print(f"{name}: {passed}/{total} passed ({100*passed/max(total,1):.1f}%)")
    print(f"  STColl: {st_pass}/{st_total}")
    print(f"  Pass avg: res={avg_res:.3f}  oU={avg_ou:.3f}  oV={avg_ov:.3f}")
    print(f"  Scaling: {len(scaling_pats)} patterns, worst={worst_ratio:.2f}x ({worst_pat})")
    if ratios:
        print(f"  Top-5:")
        for r, p in ratios[:5]:
            print(f"    {r:.2f}x  {p}")
    print(f"  SCORE: {score:.2f}")

if __name__ == '__main__':
    from mr3_gk import bidiag_svd as bidiag_svd_mr3gk

    for algo_name, fn in [
        ("DBDSQR (LAPACK, O(n^3))", bidiag_svd_dbdsqr),
        ("MR3-GK (ours, O(n^2))", bidiag_svd_mr3gk),
    ]:
        print(f"\n>>> Running {algo_name}...")
        evaluate_algo(algo_name, fn)
