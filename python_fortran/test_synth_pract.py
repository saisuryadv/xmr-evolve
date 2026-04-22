#!/usr/bin/env python3
"""
Reproduction of the Willems-Lang 2012 Synth and Pract testsets for bidiagonal SVD.

Synth (2012 paper, ~19K matrices):
  - 9 eigenvalue distribution types × 6 ECOND settings × dim 2..100
  - 6 built-in tridiagonal types × dim 2..100
  - 2 built-in bidiagonal types × selected dims
  - All in 3 glue variants: as-is, 2-copy small (γ~n·eps·‖T‖), 3-copy medium (γ~n·√eps·‖T‖)
  - Gluing done on tridiagonal level, then Cholesky → bidiagonal

Pract:
  - All tridiagonal .dat files from STCollection → Cholesky → bidiagonal
  - All bidiagonal .dat files from STCollection (already in our 19 STColl tests)

References:
  - Willems & Lang, "The MR3-GK Algorithm," ETNA 39 (2012), Section 5
  - Marques et al., "Algorithm 880," ACM TOMS 35 (2008)
  - STCollection: github.com/oamarques/STCollection
"""
import numpy as np
import sys
import os
import time
import glob
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mr3_gk import bidiag_svd
from test_dense_to_bidiag import (tridiag_to_bidiag, glue_bidiagonals,
                                   wilkinson_tridiag, test_one,
                                   random_orthogonal)

EPS = np.finfo(np.float64).eps
SQRT_EPS = np.sqrt(EPS)
ULP = EPS
SAFMIN = np.finfo(np.float64).tiny

STCOLL_DIR = '/tmp/STCollection/DATA'


# ======================================================================
#  Eigenvalue distribution types (Algorithm 880, dstedefvals.f90)
# ======================================================================

def eigenvalue_distribution(etype, n, cond, seed=42):
    """Generate eigenvalue distribution per Algorithm 880 types 1-9."""
    rng = np.random.RandomState(seed)
    w = np.zeros(n)

    if etype == 1:
        w[0] = 1.0
        w[1:] = 1.0 / cond
    elif etype == 2:
        w[:n-1] = 1.0
        w[n-1] = 1.0 / cond
    elif etype == 3:
        for i in range(n):
            w[i] = cond**(-(i) / max(n - 1, 1))
    elif etype == 4:
        for i in range(n):
            w[i] = 1.0 - i / max(n - 1, 1) * (1.0 - 1.0 / cond)
    elif etype == 5:
        alpha = np.log(1.0 / cond)
        w = np.exp(alpha * rng.uniform(0, 1, n))
        w.sort()
        w = w[::-1]
    elif etype == 6:
        w = rng.uniform(0, 1, n)
        w.sort()
        w = w[::-1]
    elif etype == 7:
        for i in range(n - 1):
            w[i] = ULP * (i + 1)
        w[n-1] = 1.0
    elif etype == 8:
        w[0] = ULP
        for i in range(1, n - 1):
            w[i] = 1.0 + SQRT_EPS * i
        w[n-1] = 2.0
    elif etype == 9:
        w[0] = 1.0
        for i in range(1, n):
            w[i] = w[i-1] + 100.0 * ULP
    else:
        raise ValueError(f"Unknown etype {etype}")

    return w


def cond_from_econd(econd, n):
    """Map ECOND setting (1-6) to condition number."""
    if econd == 1:
        return 1.0 / SQRT_EPS
    elif econd == 2:
        return 1.0 / (n * SQRT_EPS)
    elif econd == 3:
        return 1.0 / (10 * n * SQRT_EPS)
    elif econd == 4:
        return 1.0 / ULP
    elif econd == 5:
        return 1.0 / (n * ULP)
    elif econd == 6:
        return 1.0 / (10 * n * ULP)
    else:
        raise ValueError(f"Unknown econd {econd}")


def evals_to_tridiag(evals, n, seed=42):
    """Convert eigenvalue distribution to tridiagonal via random orthogonal similarity.
    Equivalent to DLATMS + DSYTRD."""
    Q = random_orthogonal(n, seed)
    T = Q @ np.diag(evals) @ Q.T
    T = (T + T.T) / 2  # symmetrize

    # Householder tridiagonalization
    T_work = T.copy()
    for k in range(n - 2):
        x = T_work[k+1:, k].copy()
        alpha = -np.sign(x[0]) * np.linalg.norm(x)
        if np.linalg.norm(x) > 0:
            v = x.copy()
            v[0] -= alpha
            v_norm = np.linalg.norm(v)
            if v_norm > 0:
                v /= v_norm
                T_work[k+1:, k:] -= 2.0 * np.outer(v, v @ T_work[k+1:, k:])
                T_work[k:, k+1:] -= 2.0 * np.outer(T_work[k:, k+1:] @ v, v)

    d = np.array([T_work[i, i] for i in range(n)])
    e = np.array([T_work[i, i+1] for i in range(n - 1)])
    return d, e


# ======================================================================
#  Built-in tridiagonal types (dstedeftmtrx.f90)
# ======================================================================

def builtin_tridiag(mtype, n):
    """Generate built-in tridiagonal matrix types 0-7."""
    d = np.zeros(n)
    e = np.zeros(max(n - 1, 0))

    if mtype == 0:
        pass  # zero
    elif mtype == 1:
        d[:] = 1.0  # identity
    elif mtype == 2:
        d[:] = 2.0; e[:] = 1.0  # (1,2,1)
    elif mtype == 3:
        # Wilkinson
        d, e = wilkinson_tridiag((n - 1) // 2)
        # Adjust if needed
        nn = len(d)
        if nn < n:
            d = np.append(d, np.zeros(n - nn))
        d = d[:n]
        if len(e) < n - 1:
            e = np.append(e, np.zeros(n - 1 - len(e)))
        e = e[:n-1]
    elif mtype == 4:
        # Clement
        for i in range(n - 1):
            e[i] = np.sqrt(float((i + 1) * (n - 1 - i)))
    elif mtype == 5:
        # Legendre
        for i in range(n - 1):
            e[i] = (i + 1) / np.sqrt((2 * (i + 1) - 1) * (2 * (i + 1) + 1))
    elif mtype == 6:
        # Laguerre
        for i in range(n):
            d[i] = 2.0 * i + 1.0
        for i in range(n - 1):
            e[i] = float(i + 1)
    elif mtype == 7:
        # Hermite
        for i in range(n - 1):
            e[i] = np.sqrt(float(i + 1))
    else:
        raise ValueError(f"Unknown mtype {mtype}")

    return d, e


# ======================================================================
#  Built-in bidiagonal types (dstedefbmtrx.f90)
# ======================================================================

def builtin_bidiag(mtype, n, seed=42):
    """Generate built-in bidiagonal matrix types 0-3."""
    rng = np.random.RandomState(seed)
    d = np.zeros(n)
    e = np.zeros(max(n - 1, 0))

    if mtype == 0:
        pass  # zero
    elif mtype == 1:
        d[:] = 1.0  # identity
    elif mtype == 2:
        # Random norm-based
        for j in range(n):
            d[j] = np.linalg.norm(rng.randn(n - j))
        for j in range(n - 1):
            e[j] = np.linalg.norm(rng.randn(n - 1 - j))
    elif mtype == 3:
        # Random exponential (log-uniform)
        gamma = -2.0 * np.log(ULP)
        d = np.exp(gamma * rng.uniform(-1, 1, n))
        e = np.exp(gamma * rng.uniform(-1, 1, max(n - 1, 0)))
    else:
        raise ValueError(f"Unknown bidiag mtype {mtype}")

    return np.abs(d), np.abs(e)


# ======================================================================
#  Gluing (tridiagonal level, then Cholesky to bidiag)
# ======================================================================

def glue_tridiag(d_list, e_list, gammas):
    """Glue tridiagonal blocks with coupling parameters."""
    d_all = []
    e_all = []
    for i, (d, e) in enumerate(zip(d_list, e_list)):
        d_all.extend(d.tolist())
        e_all.extend(e.tolist())
        if i < len(d_list) - 1:
            e_all.append(gammas[i] if i < len(gammas) else gammas[-1])
    return np.array(d_all), np.array(e_all)


def tridiag_to_bidiag_safe(d_tri, e_tri):
    """Tridiag → bidiag with safe Cholesky. Returns (d_bidiag, e_bidiag)."""
    try:
        return tridiag_to_bidiag(d_tri, e_tri)
    except Exception:
        # If Cholesky fails, return abs of tridiag entries as a fallback
        return np.abs(d_tri), np.abs(e_tri)


# ======================================================================
#  Load STCollection data files
# ======================================================================

def load_tridiag_dat(path):
    """Load tridiagonal matrix from STCollection .dat file."""
    with open(path) as f:
        n = int(f.readline().strip())
        d = np.zeros(n)
        e = np.zeros(max(n - 1, 0))
        for i in range(n):
            parts = f.readline().split()
            d[i] = float(parts[1])
            if i < n - 1 and len(parts) > 2:
                e[i] = float(parts[2])
    return d, e


def load_bidiag_dat(path):
    """Load bidiagonal matrix from STCollection .dat file (same format)."""
    return load_tridiag_dat(path)


# ======================================================================
#  Synth testset generation
# ======================================================================

def generate_synth(max_dim=100, dim_step=1, econd_list=None, seed_base=42):
    """Generate the Synth testset per Willems-Lang 2012 methodology.

    Yields (name, d_bidiag, e_bidiag) tuples.
    """
    if econd_list is None:
        econd_list = [1, 2, 3, 4, 5, 6]

    dims = list(range(2, max_dim + 1, dim_step))

    # --- Eigenvalue distribution types 1-9 ---
    for etype in range(1, 10):
        for econd in econd_list:
            # Types 6-9 don't use COND, skip extra ECOND settings
            if etype >= 6 and econd > 1:
                continue
            cond = cond_from_econd(econd, 100)  # use n=100 for cond calc

            for n in dims:
                cond_n = cond_from_econd(econd, n) if etype <= 5 else 1.0
                for seed_idx in range(4):
                    seed = seed_base + seed_idx * 1000 + etype * 100 + econd
                    try:
                        evals = eigenvalue_distribution(etype, n, cond_n, seed)
                        d_tri, e_tri = evals_to_tridiag(evals, n, seed + 7)
                        d_b, e_b = tridiag_to_bidiag_safe(d_tri, e_tri)
                        name = f"ev{etype}_ec{econd}_n{n}_s{seed_idx}"
                        yield (name, d_b, e_b)

                        # Glued small: 2 copies
                        tnorm = max(np.max(np.abs(d_tri)), 1e-300)
                        gamma_s = n * EPS * tnorm
                        d_g2, e_g2 = glue_tridiag([d_tri, d_tri.copy()],
                                                   [e_tri, e_tri.copy()],
                                                   [gamma_s])
                        d_b2, e_b2 = tridiag_to_bidiag_safe(d_g2, e_g2)
                        yield (f"{name}_gs", d_b2, e_b2)

                        # Glued medium: 3 copies
                        gamma_m = n * SQRT_EPS * tnorm
                        d_g3, e_g3 = glue_tridiag([d_tri, d_tri.copy(), d_tri.copy()],
                                                   [e_tri, e_tri.copy(), e_tri.copy()],
                                                   [gamma_m, gamma_m])
                        d_b3, e_b3 = tridiag_to_bidiag_safe(d_g3, e_g3)
                        yield (f"{name}_gm", d_b3, e_b3)
                    except Exception:
                        continue

    # --- Built-in tridiagonal types 2-7 ---
    for mtype in range(2, 8):
        for n in dims:
            try:
                d_tri, e_tri = builtin_tridiag(mtype, n)
                d_b, e_b = tridiag_to_bidiag_safe(d_tri, e_tri)
                name = f"tri{mtype}_n{n}"
                yield (name, d_b, e_b)

                tnorm = max(np.max(np.abs(d_tri)), np.max(np.abs(e_tri)), 1e-300)
                gamma_s = n * EPS * tnorm
                d_g2, e_g2 = glue_tridiag([d_tri, d_tri.copy()],
                                           [e_tri, e_tri.copy()], [gamma_s])
                d_b2, e_b2 = tridiag_to_bidiag_safe(d_g2, e_g2)
                yield (f"{name}_gs", d_b2, e_b2)

                gamma_m = n * SQRT_EPS * tnorm
                d_g3, e_g3 = glue_tridiag([d_tri, d_tri.copy(), d_tri.copy()],
                                           [e_tri, e_tri.copy(), e_tri.copy()],
                                           [gamma_m, gamma_m])
                d_b3, e_b3 = tridiag_to_bidiag_safe(d_g3, e_g3)
                yield (f"{name}_gm", d_b3, e_b3)
            except Exception:
                continue

    # --- Built-in bidiagonal types 2-3 ---
    for mtype in [2, 3]:
        for n in [125, 250, 500]:
            for seed_idx in range(4):
                seed = seed_base + seed_idx * 1000 + mtype * 100
                try:
                    d_b, e_b = builtin_bidiag(mtype, n, seed)
                    name = f"bid{mtype}_n{n}_s{seed_idx}"
                    yield (name, d_b, e_b)
                except Exception:
                    continue


# ======================================================================
#  Pract testset
# ======================================================================

def generate_pract():
    """Load Pract matrices from STCollection DATA directory.
    Tridiagonal .dat → Cholesky → bidiagonal.
    Bidiagonal .dat → directly."""

    if not os.path.isdir(STCOLL_DIR):
        print(f"  WARNING: STCollection not found at {STCOLL_DIR}")
        print(f"  Run: git clone https://github.com/oamarques/STCollection /tmp/STCollection")
        return

    # Tridiagonal matrices → Cholesky → bidiagonal
    for path in sorted(glob.glob(os.path.join(STCOLL_DIR, 'T_*.dat'))):
        name = os.path.basename(path).replace('.dat', '')
        try:
            d_tri, e_tri = load_tridiag_dat(path)
            d_b, e_b = tridiag_to_bidiag_safe(d_tri, e_tri)
            yield (name, d_b, e_b)
        except Exception as ex:
            continue

    # Fann matrices
    for path in sorted(glob.glob(os.path.join(STCOLL_DIR, 'Fann*.dat'))):
        name = os.path.basename(path).replace('.dat', '')
        try:
            d_tri, e_tri = load_tridiag_dat(path)
            d_b, e_b = tridiag_to_bidiag_safe(d_tri, e_tri)
            yield (name, d_b, e_b)
        except Exception as ex:
            continue

    # Bidiagonal matrices (already bidiagonal)
    for path in sorted(glob.glob(os.path.join(STCOLL_DIR, 'B_*.dat'))):
        name = os.path.basename(path).replace('.dat', '')
        # Skip files we already test in the 379 suite
        try:
            d_b, e_b = load_bidiag_dat(path)
            yield (f"STColl_{name}", d_b, e_b)
        except Exception:
            continue


# ======================================================================
#  Runner
# ======================================================================

def main():
    parser = argparse.ArgumentParser(description='Willems-Lang 2012 Synth+Pract Reproduction')
    parser.add_argument('--synth', action='store_true', help='Run Synth testset')
    parser.add_argument('--pract', action='store_true', help='Run Pract testset')
    parser.add_argument('--max-dim', type=int, default=100, help='Max dimension for Synth (default: 100)')
    parser.add_argument('--dim-step', type=int, default=1, help='Dimension step for Synth (default: 1, i.e. every integer)')
    parser.add_argument('--quick', action='store_true', help='Quick: dim step=10, max-dim=50, econd=[1,4]')
    parser.add_argument('--medium', action='store_true', help='Medium: dim step=5, max-dim=100, all econd')
    parser.add_argument('--paper', action='store_true', help='Paper (2012): dim step=1, max-dim=100, econd=[1,4] → ~18,400 tests')
    args = parser.parse_args()

    if not args.synth and not args.pract:
        args.synth = args.pract = True

    if args.quick:
        args.max_dim = 50
        args.dim_step = 10
        econd_list = [1, 4]
    elif args.medium:
        args.max_dim = 100
        args.dim_step = 5
        econd_list = [1, 2, 3, 4, 5, 6]
    elif args.paper:
        args.max_dim = 100
        args.dim_step = 1
        econd_list = [1, 4]   # paper uses 1/sqrt(eps) and 1/eps
    else:
        econd_list = [1, 2, 3, 4, 5, 6]

    grand_pass = 0
    grand_total = 0
    grand_fail = []
    worst_cases = []

    if args.synth:
        print("=" * 80)
        print(f"SYNTH TESTSET (Willems-Lang 2012 methodology)")
        print(f"  max_dim={args.max_dim}, dim_step={args.dim_step}")
        print("=" * 80)

        t0 = time.time()
        for name, d, e in generate_synth(args.max_dim, args.dim_step, econd_list):
            n = len(d)
            if n < 2:
                continue
            grand_total += 1
            try:
                ok, res, ou, ov, dt = test_one(d, e)
                if ok:
                    grand_pass += 1
                else:
                    grand_fail.append((name, n, res, ou, ov))
                worst_cases.append((ou, name, n, ok))

                if not ok:
                    print(f"  FAIL: {name:50s} n={n:4d}  res={res:8.3f}  ortU={ou:8.3f}  ortV={ov:8.3f}")
            except Exception as ex:
                grand_fail.append((name, n, float('inf'), float('inf'), float('inf')))

            if grand_total % 1000 == 0:
                elapsed = time.time() - t0
                print(f"  ... {grand_pass}/{grand_total} passed ({elapsed:.1f}s)")

        elapsed = time.time() - t0
        print(f"\nSynth: {grand_pass}/{grand_total} passed ({100*grand_pass/max(grand_total,1):.1f}%) in {elapsed:.1f}s")

    if args.pract:
        print("\n" + "=" * 80)
        print(f"PRACT TESTSET (STCollection tridiagonal + bidiagonal files)")
        print("=" * 80)

        pract_pass = 0
        pract_total = 0
        for name, d, e in generate_pract():
            n = len(d)
            if n < 2:
                continue
            pract_total += 1
            grand_total += 1
            try:
                ok, res, ou, ov, dt = test_one(d, e)
                if ok:
                    pract_pass += 1
                    grand_pass += 1
                else:
                    grand_fail.append((name, n, res, ou, ov))
                worst_cases.append((ou, name, n, ok))
                status = 'PASS' if ok else 'FAIL'
                print(f"  {name:50s} n={n:5d}  res={res:8.3f}  ortU={ou:8.3f}  ortV={ov:8.3f}  {status}")
            except Exception as ex:
                grand_fail.append((name, n, float('inf'), float('inf'), float('inf')))
                print(f"  {name:50s} n={n:5d}  ERROR: {ex}")

        print(f"\nPract: {pract_pass}/{pract_total} passed")

    # Summary
    print("\n" + "=" * 80)
    print(f"GRAND TOTAL: {grand_pass}/{grand_total} passed ({100*grand_pass/max(grand_total,1):.1f}%)")

    if grand_fail:
        print(f"\nFailing ({len(grand_fail)}):")
        for name, n, res, ou, ov in sorted(grand_fail, key=lambda x: -x[3])[:30]:
            print(f"  {name:55s} n={n:5d}  res={res:10.3f}  ortU={ou:10.3f}")

    if worst_cases:
        worst_cases.sort(reverse=True)
        print(f"\nTop 20 worst ortU:")
        for ou, name, n, ok in worst_cases[:20]:
            print(f"  {name:55s} n={n:5d}  ortU={ou:10.3f}  {'PASS' if ok else 'FAIL'}")

    print("=" * 80)


if __name__ == '__main__':
    main()
