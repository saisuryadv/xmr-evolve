#!/usr/bin/env python3
"""Evaluate BidiagonalSVD_TGK (DBDSVDMR3, advisor's code, from origin/main
commit 15a946d) on every bidiagonal-SVD test suite in this repo.

Suites:
  pract            — self-contained-fortran-bidiagsvd/test_synth_pract.generate_pract
  synth            — self-contained-fortran-bidiagsvd/test_synth_pract.generate_synth
  379              — lapack-dbdsvr/full_eval.{adv_names,make} × {10,100,200,400}
                     + lapack-dbdsvr/stcollection/B_*.dat
  dense_to_bidiag  — self-contained-fortran-bidiagsvd/test_dense_to_bidiag
                     (DENSE_TESTS + PAPER_TESTS) × {10,100,200,400}

For each matrix (name, d, e): write a temp STCollection bidiag .dat,
subprocess BidiagonalSVD_TGK/test_stcoll_alloc, parse res/(n.eps.||B||),
orthU/(n.eps), orthV/(n.eps). Threshold: res<=7, ortU/ortV<=5 (same as
test_dense_to_bidiag.test_one).
"""
import argparse
import glob
import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, 'BidiagonalSVD_TGK', 'test_stcoll_alloc')
SC = os.path.join(HERE, 'self-contained-fortran-bidiagsvd')
LDBD = os.path.join(HERE, 'lapack-dbdsvr')

RES_THRESH = 7.0
ORTHO_THRESH = 5.0

LINE2_RE = re.compile(
    r'(?<!paper\.)resid/\(n\.eps\.\|\|B\|\|\)=\s*([-\d.E+]+)\s+'
    r'orthU/\(n\.eps\)=\s*([-\d.E+]+)\s+'
    r'orthV/\(n\.eps\)=\s*([-\d.E+]+)'
)
ABS_RE = re.compile(
    r'(?<!paper\.)rel\.resid=\s*([-\d.E+]+)\s+'
    r'orthU=\s*([-\d.E+]+)\s+'
    r'orthV=\s*([-\d.E+]+)'
)
PAPER_NEPS_RE = re.compile(
    r'paper\.resid/\(n\.eps\.\|\|B\|\|\)=\s*([-\d.E+]+)\s+'
    r'paper\.orthU/\(n\.eps\)=\s*([-\d.E+]+)\s+'
    r'paper\.orthV/\(n\.eps\)=\s*([-\d.E+]+)'
)
PAPER_ABS_RE = re.compile(
    r'paper\.rel\.resid=\s*([-\d.E+]+)\s+'
    r'paper\.orthU=\s*([-\d.E+]+)\s+'
    r'paper\.orthV=\s*([-\d.E+]+)'
)
INFO_RE = re.compile(r'\bINFO=\s*(-?\d+)\b')
EPS = 2.2204460492503131e-16


# ----- subprocess core -----

def write_dat(path, d, e):
    n = len(d)
    with open(path, 'w') as fh:
        fh.write(f'{n}\n')
        for i in range(n):
            ei = e[i] if i < n - 1 else 0.0
            fh.write(f'{i+1} {d[i]:.17e} {ei:.17e}\n')


def _parse_advisor_output(out, n, paper_norms):
    """Return (res, ortU, ortV) in n.eps units, or None on parse-fail."""
    neps_re, abs_re = (PAPER_NEPS_RE, PAPER_ABS_RE) if paper_norms \
                      else (LINE2_RE, ABS_RE)
    m = neps_re.search(out)
    if m and '*' not in m.group(0):
        return float(m.group(1)), float(m.group(2)), float(m.group(3))
    am = abs_re.search(out)
    if not am:
        return None
    scale = max(n * EPS, 1e-300)
    return (float(am.group(1)) / scale,
            float(am.group(2)) / scale,
            float(am.group(3)) / scale)


def run_one_advisor(name, d, e, timeout, paper_norms):
    n = len(d)
    if n < 2:
        return True, 0.0, 0.0, 0.0, 0.0, 'n<2'
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat',
                                     prefix='eval_', delete=False) as tf:
        tmp = tf.name
    try:
        write_dat(tmp, d, e)
        t0 = time.perf_counter()
        try:
            r = subprocess.run([BIN, tmp], capture_output=True, text=True,
                               timeout=timeout)
        except subprocess.TimeoutExpired:
            return False, float('inf'), float('inf'), float('inf'), timeout, 'TIMEOUT'
        dt = time.perf_counter() - t0
        out = r.stdout + r.stderr
        if r.returncode != 0:
            return False, float('inf'), float('inf'), float('inf'), dt, f'rc={r.returncode}'
        info_m = INFO_RE.search(out)
        info = int(info_m.group(1)) if info_m else 0
        parsed = _parse_advisor_output(out, n, paper_norms)
        if parsed is None:
            return False, float('inf'), float('inf'), float('inf'), dt, 'parse-fail'
        res, ou, ov = parsed
        ok = (info == 0 and res <= RES_THRESH
              and ou <= ORTHO_THRESH and ov <= ORTHO_THRESH)
        note = '' if info == 0 else f'INFO={info}'
        return ok, res, ou, ov, dt, note
    finally:
        try:
            os.unlink(tmp)
        except OSError:
            pass


_SC_bidiag_svd = None


def _get_selfcontained_solver():
    """Lazy import of self-contained-fortran-bidiagsvd.mr3_gk.bidiag_svd."""
    global _SC_bidiag_svd
    if _SC_bidiag_svd is None:
        sys.path.insert(0, SC)
        from mr3_gk import bidiag_svd  # noqa: E402
        _SC_bidiag_svd = bidiag_svd
    return _SC_bidiag_svd


def _compute_metrics_py(d, e, sigma, U, V, paper_norms):
    """Compute (res, ortU, ortV) in n.eps units in NumPy. Matches the
    Fortran driver: paper-norms uses per-triplet 2-norm residual; the
    orthogonality metric is the elementwise max in both modes (this is
    what Willems-Lang 2012 Table 5.1 actually reports)."""
    n = len(d)
    if sigma.size == 0:
        return 0.0, 0.0, 0.0
    # Normalize: use only the first n-1 super-diagonal entries; pad to length n
    # so [:-1] slicing always gives length n-1 regardless of input.
    e_full = np.zeros(n, dtype=np.float64)
    e_full[:min(len(e), n - 1)] = e[:min(len(e), n - 1)]
    e = e_full
    Bnorm = max(float(sigma.max()), 1e-300)
    scale = max(n * EPS, 1e-300)
    # orthogonality (elementwise max — paper convention)
    GU = U.T @ U
    GV = V.T @ V
    np.fill_diagonal(GU, GU.diagonal() - 1.0)
    np.fill_diagonal(GV, GV.diagonal() - 1.0)
    ortU = float(np.max(np.abs(GU))) / scale
    ortV = float(np.max(np.abs(GV))) / scale
    # residual
    if paper_norms:
        # Per-triplet: max_i max(||B v_i - u_i s_i||_2, ||B^T u_i - v_i s_i||_2)
        # B is upper bidiag: (B V)_i = d_i V_i + e_i V_{i+1}
        BV = d[:, None] * V
        if n > 1:
            BV[:-1, :] += e[:-1, None] * V[1:, :]
        BtU = d[:, None] * U
        if n > 1:
            BtU[1:, :] += e[:-1, None] * U[:-1, :]
        R1 = BV - U * sigma
        R2 = BtU - V * sigma
        rn1 = np.linalg.norm(R1, axis=0)
        rn2 = np.linalg.norm(R2, axis=0)
        res = max(float(rn1.max()), float(rn2.max())) / Bnorm
    else:
        # Reconstruction max-norm
        recon = (U * sigma) @ V.T
        B = np.diag(d).astype(np.float64)
        if n > 1:
            idx = np.arange(n - 1)
            B[idx, idx + 1] = e[:-1]
        res = float(np.max(np.abs(recon - B))) / Bnorm
    return res / scale, ortU, ortV


def run_one_selfcontained(name, d, e, timeout, paper_norms):
    n = len(d)
    if n < 2:
        return True, 0.0, 0.0, 0.0, 0.0, 'n<2'
    solver = _get_selfcontained_solver()
    t0 = time.perf_counter()
    try:
        sigma, U, V, info = solver(d, e)
    except subprocess.TimeoutExpired:
        return False, float('inf'), float('inf'), float('inf'), \
               time.perf_counter() - t0, 'TIMEOUT'
    except Exception as ex:
        return False, float('inf'), float('inf'), float('inf'), \
               time.perf_counter() - t0, f'err:{type(ex).__name__}'
    dt = time.perf_counter() - t0
    if info != 0:
        return False, float('inf'), float('inf'), float('inf'), dt, f'INFO={info}'
    # Detect non-finite output before metric loop (paper-norm would NaN otherwise).
    if not (np.all(np.isfinite(sigma)) and np.all(np.isfinite(U))
            and np.all(np.isfinite(V))):
        return False, float('inf'), float('inf'), float('inf'), dt, 'non-finite'
    res, ou, ov = _compute_metrics_py(np.asarray(d, dtype=np.float64),
                                       np.asarray(e, dtype=np.float64),
                                       sigma, U, V, paper_norms)
    ok = res <= RES_THRESH and ou <= ORTHO_THRESH and ov <= ORTHO_THRESH
    return ok, res, ou, ov, dt, ''


def run_one(name, d, e, timeout, solver='advisor', paper_norms=False):
    if solver == 'advisor':
        return run_one_advisor(name, d, e, timeout, paper_norms)
    if solver == 'selfcontained':
        return run_one_selfcontained(name, d, e, timeout, paper_norms)
    raise ValueError(solver)


# ----- generator adapters -----

def gen_pract():
    sys.path.insert(0, SC)
    from test_synth_pract import generate_pract
    for name, d, e in generate_pract():
        yield name, np.asarray(d), np.asarray(e)


def gen_synth(mode='paper'):
    sys.path.insert(0, SC)
    from test_synth_pract import generate_synth
    if mode == 'quick':
        max_dim, step, econd = 50, 10, [1, 4]
    elif mode == 'medium':
        max_dim, step, econd = 100, 5, [1, 2, 3, 4, 5, 6]
    elif mode == 'paper':
        max_dim, step, econd = 100, 1, [1, 4]
    elif mode == 'full':
        max_dim, step, econd = 100, 1, [1, 2, 3, 4, 5, 6]
    else:
        raise ValueError(mode)
    for name, d, e in generate_synth(max_dim, step, econd):
        yield name, np.asarray(d), np.asarray(e)


def gen_dense_to_bidiag():
    sys.path.insert(0, SC)
    from test_dense_to_bidiag import (
        make_dense, make_paper, DENSE_TESTS, PAPER_TESTS, FIXED_SIZE_TESTS)
    sizes = [10, 100, 200, 400]
    for name in DENSE_TESTS:
        for sz in sizes:
            try:
                d, e = make_dense(name, sz)
            except Exception as ex:
                yield f'dense:{name}@{sz}', None, str(ex)
                continue
            yield f'dense:{name}@{sz}', np.asarray(d), np.asarray(e)
    for name in PAPER_TESTS:
        if name in FIXED_SIZE_TESTS:
            try:
                d, e = make_paper(name, 0)
            except Exception as ex:
                yield f'paper:{name}', None, str(ex)
                continue
            yield f'paper:{name}', np.asarray(d), np.asarray(e)
        else:
            for sz in sizes:
                try:
                    d, e = make_paper(name, sz)
                except Exception as ex:
                    yield f'paper:{name}@{sz}', None, str(ex)
                    continue
                yield f'paper:{name}@{sz}', np.asarray(d), np.asarray(e)


def _load_full_eval_make():
    """Exec only the prefix of lapack-dbdsvr/full_eval.py up to (but not
    including) the module-level adversarial test loop, so we get `make` and
    `adv_names` without triggering 270 calls into libxmr.so at import."""
    src_path = os.path.join(LDBD, 'full_eval.py')
    with open(src_path) as fh:
        src = fh.read()
    cut = src.find('\ntest_sizes = sorted(')
    if cut > 0:
        src = src[:cut]
    src = src.replace('from mr3_gk import bidiag_svd', '# stripped')
    ns = {'__name__': '__not_main__'}
    exec(compile(src, src_path, 'exec'), ns)
    return ns['make'], ns['adv_names']


def _load_stcoll(path):
    """Mirror of lapack-dbdsvr/evaluate.py:load_stcoll."""
    with open(path) as fh:
        n = int(fh.readline().strip())
        d = np.zeros(n); e = np.zeros(max(n - 1, 0))
        for i in range(n):
            parts = fh.readline().split()
            d[i] = float(parts[1])
            if i < n - 1 and len(parts) > 2:
                e[i] = float(parts[2])
    return os.path.basename(path).replace('.dat', ''), d, e


def gen_379():
    make, adv_names = _load_full_eval_make()
    sizes = [10, 100, 200, 400]
    for sz in sizes:
        for name in adv_names:
            try:
                d, e = make(name, sz)
            except Exception as ex:
                yield f'adv:{name}@{sz}', None, str(ex)
                continue
            yield f'adv:{name}@{sz}', np.asarray(d), np.asarray(e)
    for path in sorted(glob.glob(os.path.join(LDBD, 'stcollection', 'B_*.dat'))):
        try:
            name, d, e = _load_stcoll(path)
        except Exception as ex:
            yield f'stcoll:{os.path.basename(path)}', None, str(ex)
            continue
        yield f'stcoll:{name}', np.asarray(d), np.asarray(e)


# ----- driver -----

SUITE_DISPATCH = {
    'pract':           ('generate_pract',
                        'self-contained-fortran-bidiagsvd/test_synth_pract.py:358'),
    'synth':           ('generate_synth',
                        'self-contained-fortran-bidiagsvd/test_synth_pract.py:268'),
    '379':             ('full_eval.{make,adv_names} + stcollection/B_*.dat',
                        'lapack-dbdsvr/full_eval.py:23, evaluate.py:30'),
    'dense_to_bidiag': ('DENSE_TESTS+PAPER_TESTS',
                        'self-contained-fortran-bidiagsvd/test_dense_to_bidiag.py:800'),
}


def get_generator(suite, synth_mode):
    if suite == 'pract':
        return gen_pract()
    if suite == 'synth':
        return gen_synth(synth_mode)
    if suite == '379':
        return gen_379()
    if suite == 'dense_to_bidiag':
        return gen_dense_to_bidiag()
    raise ValueError(suite)


def run_suite(suite, synth_mode, max_n, timeout, limit, smoke, only=None,
              solver='advisor', paper_norms=False):
    gen_name, gen_src = SUITE_DISPATCH[suite]
    if solver == 'advisor':
        solver_label = 'BidiagonalSVD_TGK / DBDSVDMR3 (advisor)'
        source = 'python_fortran/BidiagonalSVD_TGK/  (advisor, origin/main 15a946d)'
    else:
        solver_label = 'self-contained-fortran-bidiagsvd / mr3gk_run (lab)'
        source = 'python_fortran/self-contained-fortran-bidiagsvd/'
    norm_label = ('Willems-Lang 2012 paper-norm (res=per-triplet 2-norm)'
                  if paper_norms else 'reconstruction max-norm (default)')
    print('=' * 110)
    print(f"BidiagonalSVD evaluation -- suite={suite}  solver={solver}  "
          f"norms={'paper' if paper_norms else 'maxnorm'}"
          + (f"  synth-mode={synth_mode}" if suite == 'synth' else ''))
    print('=' * 110)
    print(f"Solver:     {solver_label}")
    print(f"Source:     {source}")
    print(f"Metric:     {norm_label}")
    print(f"Generator:  {gen_name}  ({gen_src})")
    print(f"Thresholds: res<=7 n.eps,  ortU/ortV<=5 n.eps")
    print(f"Timeout:    {timeout}s/matrix")
    if smoke:
        print(f"SMOKE MODE: first 5 matrices only")
    if max_n:
        print(f"max-n:      skip matrices with n > {max_n}")
    if limit:
        print(f"limit:      first {limit} matrices")
    print('=' * 110)
    print()
    print(f"{'Matrix':<54} {'n':>6} {'res':>10} {'ortU':>10} {'ortV':>10} {'dt(s)':>8}  status  note")
    print('-' * 130)

    n_pass = n_total = 0
    fails = []
    t_start = time.time()
    cap = 5 if smoke else (limit if limit else None)

    for idx, item in enumerate(get_generator(suite, synth_mode)):
        if cap and idx >= cap:
            break
        name, d, e = item
        if only and name not in only:
            continue
        if d is None:
            # generator error
            print(f"{name:<54}     -            ERROR: {e}")
            continue
        n = len(d)
        if n < 2:
            continue
        if max_n and n > max_n:
            print(f"{name:<54} {n:>6}  (skipped: n > {max_n})")
            continue
        n_total += 1
        ok, res, ou, ov, dt, note = run_one(name, d, e, timeout,
                                            solver=solver,
                                            paper_norms=paper_norms)
        status = 'PASS' if ok else 'FAIL'
        if ok:
            n_pass += 1
        else:
            fails.append((name, n, res, ou, ov, note))
        if (n_total % 100 == 0) or smoke or (suite != 'synth'):
            print(f"{name:<54} {n:>6} {res:>10.3f} {ou:>10.3f} {ov:>10.3f} {dt:>8.3f}  {status}  {note}")
            sys.stdout.flush()
        elif suite == 'synth' and n_total % 100 == 0:
            elapsed = time.time() - t_start
            print(f"  [progress] {n_pass}/{n_total} passed  ({elapsed:.1f}s elapsed)")
            sys.stdout.flush()

    print('-' * 130)
    elapsed = time.time() - t_start
    print(f"PASS: {n_pass}/{n_total}    ({elapsed:.1f}s wall)")
    if fails:
        print()
        print(f"Failing ({len(fails)}), sorted by max(ortU,ortV) desc:")
        fails_sorted = sorted(fails, key=lambda x: -max(x[3], x[4]))
        for name, n, res, ou, ov, note in fails_sorted:
            print(f"  {name:<54} n={n:>6}  res={res:>10.3f}  ortU={ou:>10.3f}  ortV={ov:>10.3f}  {note}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--suite', required=True,
                    choices=['pract', 'synth', '379', 'dense_to_bidiag', 'all'])
    ap.add_argument('--synth-mode', default='paper',
                    choices=['quick', 'medium', 'paper', 'full'])
    ap.add_argument('--max-n', type=int, default=0,
                    help='Skip matrices with n>MAX_N (0 = no limit).')
    ap.add_argument('--timeout', type=float, default=600.0,
                    help='Per-matrix timeout in seconds.')
    ap.add_argument('--limit', type=int, default=0,
                    help='Run only first LIMIT matrices.')
    ap.add_argument('--smoke', action='store_true',
                    help='Run only first 5 matrices (smoke test).')
    ap.add_argument('--only', default='',
                    help='Comma-separated subset of matrix names to run.')
    ap.add_argument('--solver', default='advisor',
                    choices=['advisor', 'selfcontained', 'both'],
                    help='advisor: BidiagonalSVD_TGK; selfcontained: '
                         'self-contained-fortran-bidiagsvd; both: run both.')
    ap.add_argument('--paper-norms', action='store_true',
                    help='Use Willems-Lang 2012 Table 5.1 residual '
                         '(per-triplet 2-norm) instead of reconstruction max-norm.')
    args = ap.parse_args()
    only = set(s.strip() for s in args.only.split(',') if s.strip()) or None

    if args.solver in ('advisor', 'both') and not os.path.exists(BIN):
        print(f"advisor binary not found: {BIN}")
        print(f"build with: bash {os.path.dirname(BIN)}/build.sh")
        sys.exit(1)
    sc_bin = os.path.join(SC, 'mr3gk_fortran', 'mr3gk_run')
    if args.solver in ('selfcontained', 'both') and not os.path.exists(sc_bin):
        print(f"self-contained binary not found: {sc_bin}")
        print(f"build with: bash {os.path.dirname(sc_bin)}/build.sh")
        sys.exit(1)

    solvers = ['advisor', 'selfcontained'] if args.solver == 'both' \
              else [args.solver]
    suites = ('pract', 'synth', '379', 'dense_to_bidiag') \
             if args.suite == 'all' else (args.suite,)

    for solver in solvers:
        for s in suites:
            run_suite(s, args.synth_mode, args.max_n, args.timeout,
                      args.limit, args.smoke, only,
                      solver=solver, paper_norms=args.paper_norms)
            print()


if __name__ == '__main__':
    main()
