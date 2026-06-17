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
    r'resid/\(n\.eps\.\|\|B\|\|\)=\s*([-\d.E+]+)\s+'
    r'orthU/\(n\.eps\)=\s*([-\d.E+]+)\s+'
    r'orthV/\(n\.eps\)=\s*([-\d.E+]+)'
)
ABS_RE = re.compile(
    r'rel\.resid=\s*([-\d.E+]+)\s+'
    r'orthU=\s*([-\d.E+]+)\s+'
    r'orthV=\s*([-\d.E+]+)'
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


def run_one(name, d, e, timeout):
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
        m = LINE2_RE.search(out)
        if m and '*' not in m.group(0):
            res = float(m.group(1)); ou = float(m.group(2)); ov = float(m.group(3))
        else:
            # F11.1 overflow on line 2 (catastrophic orth >~ 1e9 n.eps) — fall
            # back to the absolute-value line (E10.3 format) and rescale.
            am = ABS_RE.search(out)
            if not am:
                return False, float('inf'), float('inf'), float('inf'), dt, 'parse-fail'
            scale = max(n * EPS, 1e-300)
            res = float(am.group(1)) / scale
            ou  = float(am.group(2)) / scale
            ov  = float(am.group(3)) / scale
        ok = (info == 0 and res <= RES_THRESH
              and ou <= ORTHO_THRESH and ov <= ORTHO_THRESH)
        note = '' if info == 0 else f'INFO={info}'
        return ok, res, ou, ov, dt, note
    finally:
        try:
            os.unlink(tmp)
        except OSError:
            pass


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


def run_suite(suite, synth_mode, max_n, timeout, limit, smoke, only=None):
    gen_name, gen_src = SUITE_DISPATCH[suite]
    print('=' * 110)
    print(f"BidiagonalSVD_TGK (DBDSVDMR3) -- suite={suite}"
          + (f"  synth-mode={synth_mode}" if suite == 'synth' else ''))
    print('=' * 110)
    print(f"Binary:     {BIN}")
    print(f"Source:     python_fortran/BidiagonalSVD_TGK/  (advisor's code, origin/main 15a946d)")
    print(f"Driver:     DBDSVDMR3  (dbdsvdmr3.f)")
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
        ok, res, ou, ov, dt, note = run_one(name, d, e, timeout)
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
    args = ap.parse_args()
    only = set(s.strip() for s in args.only.split(',') if s.strip()) or None

    if not os.path.exists(BIN):
        print(f"binary not found: {BIN}")
        print(f"build with: bash {os.path.dirname(BIN)}/build.sh")
        sys.exit(1)

    if args.suite == 'all':
        for s in ('pract', 'synth', '379', 'dense_to_bidiag'):
            run_suite(s, args.synth_mode, args.max_n, args.timeout, args.limit, args.smoke, only)
            print()
    else:
        run_suite(args.suite, args.synth_mode, args.max_n, args.timeout,
                  args.limit, args.smoke, only)


if __name__ == '__main__':
    main()
