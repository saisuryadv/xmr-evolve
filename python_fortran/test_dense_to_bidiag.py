#!/usr/bin/env python3
"""
New test suite for MR3-GK bidiagonal SVD.

PART 1: Dense-to-Bidiagonal tests — construct dense matrices with prescribed
        singular value distributions, reduce to bidiagonal via Householder,
        test bidiag_svd on the result.

PART 2: Missing paper test cases — bidiagonal matrices from the literature
        not covered in the existing 379-test suite.

Sources: LAPACK dchkbd (Types 1-16), Willems-Lang 2012/2013, Demmel 2008,
         Marques 2020, Grosser-Lang 2001 (DMATGEN.f), Dhillon-Parlett-Vomel 2005.
"""
import numpy as np
import sys
import os
import time
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mr3_gk import bidiag_svd
from full_eval import make, DetRNG

EPS = np.finfo(np.float64).eps
SQRT_EPS = np.sqrt(EPS)
SAFMIN = np.finfo(np.float64).tiny
RES_THRESH = 7.0
ORTHO_THRESH = 5.0


# ======================================================================
#  Utilities
# ======================================================================

def householder_bidiag(A):
    """Reduce dense m x n matrix A (m >= n) to upper bidiagonal form.
    Returns (d, e) where d is diagonal (length n) and e is superdiagonal (length n-1).
    Uses Householder reflections (equivalent to LAPACK DGEBRD)."""
    m, n = A.shape
    assert m >= n, "Requires m >= n"
    A = A.copy().astype(np.float64)
    for k in range(n):
        # Left Householder to zero out A[k+1:, k]
        x = A[k:, k].copy()
        alpha = -np.sign(x[0]) * np.linalg.norm(x)
        if alpha == 0:
            A[k, k] = 0.0
        else:
            v = x.copy()
            v[0] -= alpha
            v_norm = np.linalg.norm(v)
            if v_norm > 0:
                v /= v_norm
                A[k:, k:] -= 2.0 * np.outer(v, v @ A[k:, k:])
        # Right Householder to zero out A[k, k+2:]
        if k < n - 1:
            x = A[k, k+1:].copy()
            alpha = -np.sign(x[0]) * np.linalg.norm(x)
            if alpha == 0:
                A[k, k+1] = 0.0
            else:
                v = x.copy()
                v[0] -= alpha
                v_norm = np.linalg.norm(v)
                if v_norm > 0:
                    v /= v_norm
                    A[k:, k+1:] -= 2.0 * np.outer(A[k:, k+1:] @ v, v)
    d = np.array([A[i, i] for i in range(n)])
    e = np.array([A[i, i+1] for i in range(n - 1)])
    return d, e


def random_orthogonal(n, seed=42):
    """Generate a random orthogonal matrix via QR of random Gaussian."""
    rng = np.random.RandomState(seed)
    H = rng.randn(n, n)
    Q, R = np.linalg.qr(H)
    # Fix sign to make decomposition unique
    Q = Q @ np.diag(np.sign(np.diag(R)))
    return Q


def dense_from_svs(sigma, seed_u=42, seed_v=137):
    """Construct dense A = U @ diag(sigma) @ V^T with random orthogonal U, V."""
    n = len(sigma)
    U = random_orthogonal(n, seed_u)
    V = random_orthogonal(n, seed_v)
    return U @ np.diag(sigma) @ V.T


def dense_to_bidiag(A):
    """Reduce dense matrix to bidiagonal via Householder and return (d, e)."""
    m, n = A.shape
    if m < n:
        # Work with A^T for lower bidiagonal, then transpose result
        d, e = householder_bidiag(A.T)
    else:
        d, e = householder_bidiag(A)
    # Make entries non-negative (our bidiag_svd handles signs)
    return np.abs(d), np.abs(e)


def tridiag_to_bidiag(diag, offdiag):
    """Convert symmetric tridiagonal to bidiagonal via Cholesky.
    T - nu*I = L*L^T, B = L^T (upper bidiagonal).
    nu is chosen to make T - nu*I positive definite."""
    n = len(diag)
    # Gershgorin lower bound
    gl = diag[0] - abs(offdiag[0]) if n > 1 else diag[0]
    for i in range(1, n - 1):
        gl = min(gl, diag[i] - abs(offdiag[i-1]) - abs(offdiag[i]))
    if n > 1:
        gl = min(gl, diag[n-1] - abs(offdiag[n-2]))
    # Shift to make positive definite
    nu = gl - 0.001 * max(abs(gl), 1.0)
    # Cholesky of T - nu*I
    d_shifted = np.array(diag, dtype=np.float64) - nu
    e_off = np.array(offdiag, dtype=np.float64)
    # L*D*L^T factorization
    d_chol = np.zeros(n)
    l_chol = np.zeros(n - 1)
    d_chol[0] = d_shifted[0]
    for i in range(n - 1):
        if d_chol[i] <= 0:
            # Not positive definite, shift more
            nu -= 1.0
            d_shifted = np.array(diag, dtype=np.float64) - nu
            d_chol[0] = d_shifted[0]
            for j in range(i):
                l_chol[j] = e_off[j] / d_chol[j]
                d_chol[j+1] = d_shifted[j+1] - l_chol[j]**2 * d_chol[j]
        l_chol[i] = e_off[i] / d_chol[i]
        d_chol[i+1] = d_shifted[i+1] - l_chol[i]**2 * d_chol[i]
    # B = L^T: d_bidiag[i] = sqrt(d_chol[i]), e_bidiag[i] = l_chol[i] * sqrt(d_chol[i])
    d_bidiag = np.sqrt(np.abs(d_chol))
    e_bidiag = np.array([l_chol[i] * d_bidiag[i] for i in range(n - 1)])
    return np.abs(d_bidiag), np.abs(e_bidiag)


def wilkinson_tridiag(m):
    """Construct W_{2m+1}^+: diag = |m, m-1, ..., 0, ..., m-1, m|, offdiag = 1."""
    n = 2 * m + 1
    d = np.array([abs(m - i) for i in range(n)], dtype=np.float64)
    e = np.ones(n - 1, dtype=np.float64)
    return d, e


def glue_bidiagonals(d_list, e_list, gamma):
    """Glue multiple bidiagonal matrices together with coupling gamma.
    d_list: list of diagonal arrays, e_list: list of superdiagonal arrays.
    gamma: coupling parameter placed between blocks."""
    d_all = []
    e_all = []
    for i, (d, e) in enumerate(zip(d_list, e_list)):
        d_all.extend(d.tolist())
        e_all.extend(e.tolist())
        if i < len(d_list) - 1:
            e_all.append(gamma)
    return np.array(d_all), np.array(e_all)


def test_one(d, e):
    """Run bidiag_svd and compute metrics. Returns (ok, res, ortU, ortV, dt)."""
    n = len(d)
    if n == 0:
        return True, 0.0, 0.0, 0.0, 0.0
    t0 = time.perf_counter()
    try:
        sigma, U, V, info = bidiag_svd(d.copy(), e.copy())
    except Exception as ex:
        return False, float('inf'), float('inf'), float('inf'), time.perf_counter() - t0
    dt = time.perf_counter() - t0

    B = np.diag(d) + np.diag(e, 1)
    bnorm = abs(d[-1]) if n > 0 else 0.0
    for i in range(n - 1):
        bnorm = max(bnorm, abs(d[i]) + abs(e[i]))
    if bnorm == 0.0:
        bnorm = 1.0

    res = np.max(np.abs(B - U @ np.diag(sigma) @ V.T)) / (bnorm * n * EPS)
    ou = np.max(np.abs(U.T @ U - np.eye(n))) / (n * EPS)
    ov = np.max(np.abs(V.T @ V - np.eye(n))) / (n * EPS)
    ok = res <= RES_THRESH and ou <= ORTHO_THRESH and ov <= ORTHO_THRESH
    return ok, res, ou, ov, dt


# ======================================================================
#  PART 1: Dense-to-Bidiagonal Test Generators
# ======================================================================

def make_dense(name, n):
    """Generate a dense matrix and reduce to bidiagonal. Returns (d, e)."""

    if name == 'dlatms_mode1_eps':
        # LAPACK Mode 1: sigma(1)=1, sigma(2:n)=eps
        sigma = np.full(n, EPS)
        sigma[0] = 1.0
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dlatms_mode2_eps':
        # LAPACK Mode 2: sigma(1:n-1)=1, sigma(n)=eps
        sigma = np.ones(n)
        sigma[-1] = EPS
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dlatms_mode3_eps':
        # LAPACK Mode 3: geometric from 1 to eps
        kappa = 1.0 / EPS
        sigma = np.array([kappa**(-(i) / (n - 1)) for i in range(n)]) if n > 1 else np.array([1.0])
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dlatms_mode4_eps':
        # LAPACK Mode 4: arithmetic from 1 to eps
        sigma = np.array([1.0 - i / (n - 1) * (1.0 - EPS) for i in range(n)]) if n > 1 else np.array([1.0])
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dlatms_mode5_random':
        # LAPACK Mode 5: log-uniform random in (eps, 1)
        rng = np.random.RandomState(42)
        sigma = np.exp(rng.uniform(np.log(EPS), 0.0, n))
        sigma.sort()
        sigma = sigma[::-1]
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dense_three_clusters':
        # Three well-separated clusters
        k = n // 3
        sigma = np.concatenate([
            np.full(k, 100.0),
            np.full(k, 1.0),
            np.full(n - 2 * k, EPS)
        ])
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dense_one_tight_cluster':
        # All SVs within n*eps of 1
        sigma = np.array([1.0 + i * EPS for i in range(n)])
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dense_half_rank':
        # Rank n/2
        half = n // 2
        sigma = np.zeros(n)
        for i in range(half):
            sigma[i] = 10.0**(-(i) * 8.0 / max(half - 1, 1))
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dense_repeated_sv':
        # Exact multiplicity: half at 1, half at eps
        half = n // 2
        sigma = np.concatenate([np.ones(half), np.full(n - half, EPS)])
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dense_hilbert':
        # Hilbert matrix
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                A[i, j] = 1.0 / (i + j + 1)
        return dense_to_bidiag(A)

    elif name == 'dense_kahan':
        # Kahan matrix
        s = 0.3  # sin(theta)
        c = np.sqrt(1 - s * s)
        S = np.diag([c**i for i in range(n)])
        T = np.eye(n)
        for i in range(n):
            for j in range(i + 1, n):
                T[i, j] = -s
        A = S @ T
        return dense_to_bidiag(A)

    elif name == 'dense_condition_1e20':
        # Geometric from 1 to 1e-20
        sigma = np.array([10.0**(-(i) * 20.0 / (n - 1)) for i in range(n)]) if n > 1 else np.array([1.0])
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dense_wilkinson_sv':
        # SVs = eigenvalues of W_{2m+1}^+ (tight pairs)
        m = (n - 1) // 2
        d_w, e_w = wilkinson_tridiag(m)
        nn = len(d_w)
        T = np.diag(d_w) + np.diag(e_w, 1) + np.diag(e_w, -1)
        evals = np.sort(np.linalg.eigvalsh(T))[::-1]
        sigma = np.abs(evals[:n])
        sigma = np.sort(sigma)[::-1]
        sigma = np.maximum(sigma, SAFMIN)
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dense_wl_example41':
        # Willems-Lang 2012 Example 4.1: 20 prescribed singular values
        # sigma_13=0.9, sigma_14=1-1e-7, sigma_15=1+1e-7, sigma_16=1.1
        # sigma_i = sigma_{i+4}/100 for i=12..1, sigma_i = 100*sigma_{i-4} for i=17..20
        if n < 20:
            # For small n, use a scaled version
            sigma = np.array([10.0**(-(i) * 8.0 / (n - 1)) for i in range(n)])
        else:
            sv = np.zeros(20)
            sv[12] = 0.9
            sv[13] = 1.0 - 1e-7
            sv[14] = 1.0 + 1e-7
            sv[15] = 1.1
            for i in range(11, -1, -1):
                sv[i] = sv[i + 4] / 100.0
            for i in range(16, 20):
                sv[i] = 100.0 * sv[i - 4]
            # Extend to n if needed
            if n > 20:
                sigma = np.zeros(n)
                sigma[:20] = sv
                for i in range(20, n):
                    sigma[i] = sv[0] * 10.0**(-(i - 19) * 2)
            else:
                sigma = sv[:n]
        sigma = np.sort(np.abs(sigma))[::-1]
        sigma = np.maximum(sigma, SAFMIN)
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dense_near_overflow':
        # Near overflow: geometric from 1e150 to 1e140
        sigma = np.array([10.0**(150.0 - i * 10.0 / (n - 1)) for i in range(n)]) if n > 1 else np.array([1e150])
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dense_near_underflow':
        # Near underflow: geometric from 1e-150 to 1e-160
        sigma = np.array([10.0**(-150.0 - i * 10.0 / (n - 1)) for i in range(n)]) if n > 1 else np.array([1e-150])
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'dense_random':
        # Random dense matrix
        rng = np.random.RandomState(42)
        A = rng.uniform(-1.0, 1.0, (n, n))
        return dense_to_bidiag(A)

    elif name == 'dense_vandermonde':
        # Vandermonde on Chebyshev nodes
        nodes = np.cos(np.pi * np.arange(n) / (n - 1)) if n > 1 else np.array([1.0])
        A = np.vander(nodes, n)
        return dense_to_bidiag(A)

    elif name == 'dense_toeplitz':
        # Symmetric Toeplitz
        from scipy.linalg import toeplitz
        first_row = np.arange(n, 0, -1, dtype=np.float64)
        A = toeplitz(first_row)
        return dense_to_bidiag(A)

    elif name == 'dense_frank':
        # Frank matrix: F(i,j) = min(i,j) + 1
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                A[i, j] = float(min(i, j) + 1)
        return dense_to_bidiag(A)

    elif name == 'dense_companion_exp':
        # Companion matrix of truncated exp(z) Taylor series
        # Triggers LAPACK bug #316 (DGESDD/DBDSDC) for n >= 26
        # c = [1, 1, 1/2!, 1/3!, ..., 1/n!]
        c = np.ones(n + 1)
        for i in range(2, n + 1):
            c[i] = c[i-1] / i
        # Companion matrix (Frobenius form)
        A = np.zeros((n, n))
        for i in range(1, n):
            A[i, i-1] = 1.0
        # First row = -reversed(c[1:])/c[0]
        for j in range(n):
            A[0, j] = -c[n - j] / c[0]
        return dense_to_bidiag(A)

    elif name == 'dense_moler_like':
        # Moler-like convergence stress test for bidiagonal SVD
        # Bidiagonal with d_i near +/-1 with random signs, e_i decaying
        # Inspired by Cleve Moler's matrix that exposed DLASQ1 convergence bug
        rng = np.random.RandomState(42)
        d = np.array([(-1)**i * (1.0 - rng.uniform(0, 0.05)) for i in range(n)])
        e = np.array([0.5 * (1.0 - i / n)**2 for i in range(n - 1)])
        return np.abs(d), np.abs(e)

    else:
        raise ValueError(f"Unknown dense test: {name}")


# ======================================================================
#  PART 2: Missing Paper Test Generators
# ======================================================================

def make_paper(name, n):
    """Generate bidiagonal matrices from papers. Returns (d, e)."""

    # --- 2A: Proper Glued Wilkinson ---
    if name == 'glued_wilk_3x21_sqrteps':
        m = 10  # W_21^+
        d_w, e_w = wilkinson_tridiag(m)
        d_b, e_b = tridiag_to_bidiag(d_w, e_w)
        d_list = [d_b.copy() for _ in range(3)]
        e_list = [e_b.copy() for _ in range(3)]
        return glue_bidiagonals(d_list, e_list, SQRT_EPS)

    elif name == 'glued_wilk_5x21_sqrteps':
        m = 10
        d_w, e_w = wilkinson_tridiag(m)
        d_b, e_b = tridiag_to_bidiag(d_w, e_w)
        d_list = [d_b.copy() for _ in range(5)]
        e_list = [e_b.copy() for _ in range(5)]
        return glue_bidiagonals(d_list, e_list, SQRT_EPS)

    elif name == 'glued_wilk_3x21_1e14':
        m = 10
        d_w, e_w = wilkinson_tridiag(m)
        d_b, e_b = tridiag_to_bidiag(d_w, e_w)
        d_list = [d_b.copy() for _ in range(3)]
        e_list = [e_b.copy() for _ in range(3)]
        return glue_bidiagonals(d_list, e_list, 1e-14)

    elif name == 'glued_wilk_10x11_sqrteps':
        m = 5  # W_11^+
        d_w, e_w = wilkinson_tridiag(m)
        d_b, e_b = tridiag_to_bidiag(d_w, e_w)
        d_list = [d_b.copy() for _ in range(10)]
        e_list = [e_b.copy() for _ in range(10)]
        return glue_bidiagonals(d_list, e_list, SQRT_EPS)

    # --- 2B: Glued versions of existing patterns ---
    elif name.endswith('_glued_small') or name.endswith('_glued_medium'):
        is_small = name.endswith('_glued_small')
        base_name = name.replace('_glued_small', '').replace('_glued_medium', '')
        copies = 3 if is_small else 2
        block_n = max(n // copies, 3)
        d_base, e_base = make(base_name, block_n)
        bnorm = max(np.max(np.abs(d_base)), 1e-300)
        if is_small:
            gamma = block_n * EPS * bnorm
        else:
            gamma = block_n * SQRT_EPS * bnorm
        d_list = [d_base.copy() for _ in range(copies)]
        e_list = [e_base.copy() for _ in range(copies)]
        return glue_bidiagonals(d_list, e_list, gamma)

    # --- 2C: CHKBD LAPACK general ---
    elif name == 'chkbd_lapack_general':
        rng = np.random.RandomState(42)
        gamma = -2.0 * np.log(EPS)  # ~ 73.6
        d = np.exp(gamma * rng.uniform(-1, 1, n))
        e = np.exp(gamma * rng.uniform(-1, 1, max(n - 1, 0)))
        return np.abs(d), np.abs(e)

    # --- 2D: Corrected Grosser-Lang ---
    elif name == 'gl_gro0_proper':
        # Paper ID 240: d[0]=1, rest=alpha~eps, e=alpha
        alpha = EPS
        d = np.full(n, alpha)
        d[0] = 1.0
        e = np.full(max(n - 1, 0), alpha)
        return d, e

    elif name == 'gl_gro1_proper':
        # Paper ID 241: d[0:2]=1, rest=alpha
        alpha = EPS
        d = np.full(n, alpha)
        d[0] = 1.0
        if n > 1:
            d[1] = 1.0
        e = np.full(max(n - 1, 0), alpha)
        return d, e

    elif name == 'gl_gro2_proper':
        # Paper ID 242: d[0:4]=1, rest=alpha
        alpha = EPS
        d = np.full(n, alpha)
        for i in range(min(4, n)):
            d[i] = 1.0
        e = np.full(max(n - 1, 0), alpha)
        return d, e

    elif name == 'gl_gro3_proper':
        # Paper ID 243: d[0:4]=1, e[0:4]=1, rest=alpha
        alpha = EPS
        d = np.full(n, alpha)
        for i in range(min(4, n)):
            d[i] = 1.0
        e = np.full(max(n - 1, 0), alpha)
        for i in range(min(4, n - 1)):
            e[i] = 1.0
        return d, e

    elif name == 'gl_abcon1_proper':
        # Paper ID 201: d[:]=alpha=2, d[0]=alpha-beta=1, e[:]=beta=1
        d = np.full(n, 2.0)
        d[0] = 1.0
        e = np.ones(max(n - 1, 0))
        return d, e

    elif name == 'gl_abcon2_proper':
        # Paper ID 202: d[:]=2, d[0]=1, d[-1]=3, e[:]=1
        d = np.full(n, 2.0)
        d[0] = 1.0
        if n > 1:
            d[-1] = 3.0
        e = np.ones(max(n - 1, 0))
        return d, e

    elif name == 'gl_abcon3_proper':
        # Paper ID 203: d[:]=2, d[0]=d[-1]=3, e[:]=1
        d = np.full(n, 2.0)
        d[0] = 3.0
        if n > 1:
            d[-1] = 3.0
        e = np.ones(max(n - 1, 0))
        return d, e

    # --- 2E: Marques STEXR failure ---
    elif name.startswith('marques_stexr_'):
        nn = int(name.split('_')[-1])
        # Tridiag with eigenvalues c^((i-1)/(nn-1)), c = 1/sqrt(eps)
        c = 1.0 / SQRT_EPS
        evals = np.array([c**((i) / (nn - 1)) for i in range(nn)]) if nn > 1 else np.array([1.0])
        # Construct tridiagonal via DLATMS-like: diag matrix with these eigenvalues,
        # similarity transform to get tridiagonal
        D = np.diag(evals)
        Q = random_orthogonal(nn, seed=77)
        T = Q @ D @ Q.T
        T = (T + T.T) / 2  # Ensure symmetry
        # Tridiagonalize
        from scipy.linalg import hessenberg
        # Use eigvals to get proper tridiag
        evals_check, evecs = np.linalg.eigh(T)
        # Rebuild as tridiagonal
        from scipy.linalg import eigvalsh_tridiagonal
        # Just directly build tridiag with these eigenvalues using DSYTRD equivalent
        # Simplest: use scipy's schur for symmetric -> tridiag
        d_tri = np.diag(T).copy()
        e_tri = np.array([T[i, i+1] for i in range(nn - 1)])
        # Actually use the Householder tridiag reduction
        T_work = T.copy()
        for k in range(nn - 2):
            x = T_work[k+1:, k].copy()
            alpha_h = -np.sign(x[0]) * np.linalg.norm(x)
            if np.linalg.norm(x) > 0:
                v = x.copy()
                v[0] -= alpha_h
                v_norm = np.linalg.norm(v)
                if v_norm > 0:
                    v /= v_norm
                    T_work[k+1:, k:] -= 2.0 * np.outer(v, v @ T_work[k+1:, k:])
                    T_work[k:, k+1:] -= 2.0 * np.outer(T_work[k:, k+1:] @ v, v)
        d_tri = np.array([T_work[i, i] for i in range(nn)])
        e_tri = np.array([T_work[i, i+1] for i in range(nn - 1)])
        return tridiag_to_bidiag(d_tri, e_tri)

    # --- 2F: Tridiagonal-to-Bidiagonal ---
    elif name == 'wilkinson_cholesky_m10':
        d_w, e_w = wilkinson_tridiag(10)  # W_21^+
        return tridiag_to_bidiag(d_w, e_w)

    elif name == 'wilkinson_cholesky_m50':
        d_w, e_w = wilkinson_tridiag(50)  # W_101^+
        return tridiag_to_bidiag(d_w, e_w)

    elif name == 'legendre_bidiag':
        # Jacobi matrix for Legendre polynomials: d=0, e[i]=(i+1)/sqrt((2(i+1)-1)*(2(i+1)+1))
        d = np.zeros(n)
        e = np.array([(i + 1) / np.sqrt((2 * (i + 1) - 1) * (2 * (i + 1) + 1)) for i in range(n - 1)])
        return tridiag_to_bidiag(d, e)

    elif name == 'laguerre_bidiag':
        # Jacobi matrix for Laguerre: d[i]=2*i+1, e[i]=i+1
        d = np.array([2.0 * i + 1.0 for i in range(n)])
        e = np.array([float(i + 1) for i in range(n - 1)])
        return tridiag_to_bidiag(d, e)

    elif name == 'hermite_bidiag':
        # Jacobi matrix for Hermite: d=0, e[i]=sqrt(i+1)
        d = np.zeros(n)
        e = np.array([np.sqrt(float(i + 1)) for i in range(n - 1)])
        return tridiag_to_bidiag(d, e)

    # --- 2G: Missing variants ---
    elif name == 'demmel_S1pe_k16':
        d = np.full(n, 1e-16)
        d[0] = 1.0
        e = np.array([np.sqrt(abs(d[i] * d[i+1])) * 0.1 for i in range(n - 1)])
        return d, e

    elif name == 'wl_element_growth':
        # Willems-Lang 2012 Example 4.8: element growth vulnerability
        alpha = EPS
        d = np.array([alpha**i for i in range(n)])
        d[0] = 1.0
        e = np.array([alpha**i for i in range(n - 1)])
        e[0] = 1.0
        return d, e

    elif name == 'random_zero_e_10pct':
        rng = DetRNG(42)
        d = np.array([0.5 + rng.uniform() * 1.5 for _ in range(n)])
        e = np.array([0.1 + rng.uniform() * 0.9 for _ in range(n - 1)])
        # Zero out 10% of e
        for i in range(0, n - 1, 10):
            e[i] = 0.0
        return d, e

    elif name == 'random_zero_e_50pct':
        rng = DetRNG(42)
        d = np.array([0.5 + rng.uniform() * 1.5 for _ in range(n)])
        e = np.array([0.1 + rng.uniform() * 0.9 for _ in range(n - 1)])
        # Zero out 50% of e (every other)
        for i in range(0, n - 1, 2):
            e[i] = 0.0
        return d, e

    # --- 2H: From Dhillon-Parlett-Vomel 2005 / LAWN 163/166 / Dhillon Thesis ---
    elif name == 'bidiag_coupling_failure':
        # Grosser-Lang 2005 eq (3.2) / LAWN 166: one large SV, rest clustered at ~10^-27
        # B = diag([1, alpha, ..., alpha]) + diag([alpha, ..., alpha], +1)
        # alpha = 200 * eps. Orthogonality > 10^10 without couplings.
        alpha = 200.0 * EPS
        d = np.full(n, alpha)
        d[0] = 1.0
        e = np.full(max(n - 1, 0), alpha)
        return d, e

    elif name == 'dhillon_uniform_sqrteps_apart':
        # Dhillon Thesis Type 2: eigenvalues spaced sqrt(eps) apart
        # lambda_1=eps, lambda_i=1+(i-1)*sqrt(eps), lambda_n=2
        # Boundary case for MRRR clustering threshold (gap = sqrt(eps))
        evals = np.zeros(n)
        evals[0] = EPS
        for i in range(1, n - 1):
            evals[i] = 1.0 + i * SQRT_EPS
        evals[-1] = 2.0
        # Construct as dense matrix with prescribed eigenvalues, reduce to bidiag
        D = np.diag(evals)
        Q = random_orthogonal(n, seed=99)
        T = Q @ D @ Q.T
        T = (T + T.T) / 2
        # Tridiagonalize via Householder
        T_work = T.copy()
        for k in range(n - 2):
            x = T_work[k+1:, k].copy()
            alpha_h = -np.sign(x[0]) * np.linalg.norm(x)
            if np.linalg.norm(x) > 0:
                v = x.copy()
                v[0] -= alpha_h
                v_norm = np.linalg.norm(v)
                if v_norm > 0:
                    v /= v_norm
                    T_work[k+1:, k:] -= 2.0 * np.outer(v, v @ T_work[k+1:, k:])
                    T_work[k:, k+1:] -= 2.0 * np.outer(T_work[k:, k+1:] @ v, v)
        d_tri = np.array([T_work[i, i] for i in range(n)])
        e_tri = np.array([T_work[i, i+1] for i in range(n - 1)])
        return tridiag_to_bidiag(d_tri, e_tri)

    elif name == 'dhillon_121_toeplitz':
        # Dhillon Thesis Type 12: tridiag(1,2,1), analytically known eigenvalues
        # eigenvalues = 4*sin^2(k*pi/(2*(n+1))), k=1..n
        d_tri = np.full(n, 2.0)
        e_tri = np.ones(n - 1)
        return tridiag_to_bidiag(d_tri, e_tri)

    elif name == 'glued_wilk201_5x_sqrteps':
        # Dhillon-Parlett-Vomel 2005: THE actual failure matrix
        # 5 copies of W_{201}^+ (m=100) glued with gamma=sqrt(eps), total n=1005
        # Triggers MRRR failure: FP underflow, rep tree depth 21+
        m = 100  # W_{201}^+
        d_w, e_w = wilkinson_tridiag(m)
        d_b, e_b = tridiag_to_bidiag(d_w, e_w)
        d_list = [d_b.copy() for _ in range(5)]
        e_list = [e_b.copy() for _ in range(5)]
        return glue_bidiagonals(d_list, e_list, SQRT_EPS)

    # --- 2I: STCollection eigenvalue distribution types 7/8/9 ---
    elif name == 'stcoll_evdist_type7':
        # Type 7: W(i)=ULP*i for i=1..N-1, W(N)=1 — linearly growing from ULP
        evals = np.array([EPS * (i + 1) for i in range(n - 1)] + [1.0])
        sigma = np.sort(np.abs(evals))[::-1]
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'stcoll_evdist_type8':
        # Type 8: W(1)=ULP, W(i)=1+sqrt(ULP)*i, W(N)=2 — one tiny + tight cluster at 1
        evals = np.zeros(n)
        evals[0] = EPS
        for i in range(1, n - 1):
            evals[i] = 1.0 + SQRT_EPS * i
        evals[-1] = 2.0
        sigma = np.sort(np.abs(evals))[::-1]
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    elif name == 'stcoll_evdist_type9':
        # Type 9: W(1)=1, W(i)=W(i-1)+100*ULP — extremely tight arithmetic cluster
        # All eigenvalues within 100*n*eps of each other
        evals = np.array([1.0 + 100.0 * EPS * i for i in range(n)])
        sigma = np.sort(np.abs(evals))[::-1]
        A = dense_from_svs(sigma)
        return dense_to_bidiag(A)

    # --- 2J: Skew-Wilkinson with grading ---
    elif name == 'skew_wilk21_g1':
        # Skew-Wilkinson W21 with grading factor 1 (standard)
        m = 10
        d_w, e_w = wilkinson_tridiag(m)
        return tridiag_to_bidiag(d_w, e_w)

    elif name == 'skew_wilk21_g1e3':
        # Skew-Wilkinson W21 with grading factor 10^3
        # Off-diagonals scaled: e[i] = 10^3 for i < m, e[i] = 10^{-3} for i >= m
        m = 10
        d_w, e_w = wilkinson_tridiag(m)
        nn = len(d_w)
        for i in range(nn - 1):
            if i < m:
                e_w[i] *= 1e3
            else:
                e_w[i] *= 1e-3
        return tridiag_to_bidiag(d_w, e_w)

    elif name == 'skew_wilk21_g1e6':
        # Skew-Wilkinson W21 with grading factor 10^6
        m = 10
        d_w, e_w = wilkinson_tridiag(m)
        nn = len(d_w)
        for i in range(nn - 1):
            if i < m:
                e_w[i] *= 1e6
            else:
                e_w[i] *= 1e-6
        return tridiag_to_bidiag(d_w, e_w)

    # --- 2K: W21 with varying off-diagonal gamma ---
    elif name == 'wilk21_gamma_1e4':
        # W21 with off-diag scaled by 1e-4
        m = 10
        d_w, e_w = wilkinson_tridiag(m)
        e_w *= 1e-4
        return tridiag_to_bidiag(d_w, e_w)

    elif name == 'wilk21_gamma_1e8':
        # W21 with off-diag scaled by 1e-8
        m = 10
        d_w, e_w = wilkinson_tridiag(m)
        e_w *= 1e-8
        return tridiag_to_bidiag(d_w, e_w)

    elif name == 'wilk21_gamma_1e14':
        # W21 with off-diag scaled by 1e-14
        m = 10
        d_w, e_w = wilkinson_tridiag(m)
        e_w *= 1e-14
        return tridiag_to_bidiag(d_w, e_w)

    # --- 2L: Willems 2x2 block LDL* edge case ---
    elif name == 'willems_block_ldl_edge':
        # T = GK(B) - alpha*I where alpha = O(eps)
        # 4x4 tridiagonal that forces 2x2 blocks in LDL* factorization
        # Extend to larger n by repeating the pattern
        alpha = EPS
        block_size = 4
        num_blocks = max(n // block_size, 1)
        d_all = []
        e_all = []
        for b in range(num_blocks):
            # Each block: T = [-alpha, 1, 0, 0; 1, -alpha, 1, 0; 0, 1, -alpha, alpha; 0, 0, alpha, -alpha]
            d_all.extend([alpha, alpha, alpha, alpha])
            if b < num_blocks - 1:
                e_all.extend([1.0, 1.0, alpha, SQRT_EPS])  # last e is glue
            else:
                e_all.extend([1.0, 1.0, alpha])
        d = np.array(d_all[:n])
        e = np.array(e_all[:max(n - 1, 0)])
        return np.abs(d), np.abs(e)

    else:
        raise ValueError(f"Unknown paper test: {name}")


# ======================================================================
#  Test Lists
# ======================================================================

DENSE_TESTS = [
    'dlatms_mode1_eps', 'dlatms_mode2_eps', 'dlatms_mode3_eps',
    'dlatms_mode4_eps', 'dlatms_mode5_random',
    'dense_three_clusters', 'dense_one_tight_cluster',
    'dense_half_rank', 'dense_repeated_sv',
    'dense_hilbert', 'dense_kahan', 'dense_condition_1e20',
    'dense_wilkinson_sv', 'dense_wl_example41',
    'dense_near_overflow', 'dense_near_underflow',
    'dense_random', 'dense_vandermonde', 'dense_toeplitz', 'dense_frank',
    'dense_companion_exp', 'dense_moler_like',
]

PAPER_TESTS = [
    # 2A: Proper Glued Wilkinson (fixed size)
    'glued_wilk_3x21_sqrteps', 'glued_wilk_5x21_sqrteps',
    'glued_wilk_3x21_1e14', 'glued_wilk_10x11_sqrteps',
    # 2B: Glued versions of existing patterns
    'exponential_graded_glued_small', 'exponential_graded_glued_medium',
    'stemr_killer_glued_small', 'constant_glued_small',
    'all_equal_nontrivial_glued_small', 'demmel_S1pe_glued_small',
    'pd_T0_glued_medium', 'chkbd_glued_small',
    # 2C: CHKBD LAPACK general
    'chkbd_lapack_general',
    # 2D: Corrected Grosser-Lang
    'gl_gro0_proper', 'gl_gro1_proper', 'gl_gro2_proper', 'gl_gro3_proper',
    'gl_abcon1_proper', 'gl_abcon2_proper', 'gl_abcon3_proper',
    # 2E: Marques STEXR failure
    'marques_stexr_10', 'marques_stexr_20', 'marques_stexr_50',
    # 2F: Tridiagonal-to-Bidiagonal
    'wilkinson_cholesky_m10', 'wilkinson_cholesky_m50',
    'legendre_bidiag', 'laguerre_bidiag', 'hermite_bidiag',
    # 2G: Missing variants
    'demmel_S1pe_k16', 'wl_element_growth',
    'random_zero_e_10pct', 'random_zero_e_50pct',
    # 2H: From papers (Dhillon-Parlett-Vomel 2005 / LAWN 163/166 / Dhillon Thesis)
    'bidiag_coupling_failure', 'dhillon_uniform_sqrteps_apart',
    'dhillon_121_toeplitz', 'glued_wilk201_5x_sqrteps',
    # 2I: STCollection eigenvalue distribution types 7/8/9
    'stcoll_evdist_type7', 'stcoll_evdist_type8', 'stcoll_evdist_type9',
    # 2J: Skew-Wilkinson with grading
    'skew_wilk21_g1', 'skew_wilk21_g1e3', 'skew_wilk21_g1e6',
    # 2K: W21 with varying off-diagonal gamma
    'wilk21_gamma_1e4', 'wilk21_gamma_1e8', 'wilk21_gamma_1e14',
    # 2L: Willems 2x2 block LDL* edge case
    'willems_block_ldl_edge',
]

# Fixed-size tests (glued Wilkinson has intrinsic size)
FIXED_SIZE_TESTS = {
    'glued_wilk_3x21_sqrteps', 'glued_wilk_5x21_sqrteps',
    'glued_wilk_3x21_1e14', 'glued_wilk_10x11_sqrteps',
    'wilkinson_cholesky_m10', 'wilkinson_cholesky_m50',
    'marques_stexr_10', 'marques_stexr_20', 'marques_stexr_50',
    'glued_wilk201_5x_sqrteps',
    'skew_wilk21_g1', 'skew_wilk21_g1e3', 'skew_wilk21_g1e6',
    'wilk21_gamma_1e4', 'wilk21_gamma_1e8', 'wilk21_gamma_1e14',
}


# ======================================================================
#  Runner
# ======================================================================

def run_tests(test_names, generator, sizes, fixed_size_set, label, timeout_s=60):
    """Run a set of tests and return results dict."""
    results = {}
    total_pass = 0
    total = 0

    for name in test_names:
        if name in fixed_size_set:
            # Fixed size test — run once
            try:
                d, e = generator(name, 0)  # n is ignored for fixed-size
            except Exception as ex:
                key = name
                results[key] = (False, float('inf'), float('inf'), float('inf'), 0)
                total += 1
                print(f"  {name:45s} n={len(d) if 'd' in dir() else '?':>4s}  ERROR: {ex}")
                continue
            n = len(d)
            key = f"{name}@{n}"
            total += 1
            try:
                ok, res, ou, ov, dt = test_one(d, e)
                results[key] = (ok, res, ou, ov, dt)
                if ok:
                    total_pass += 1
                status = 'PASS' if ok else 'FAIL'
                print(f"  {name:45s} n={n:4d}  res={res:8.3f}  ortU={ou:8.3f}  ortV={ov:8.3f}  t={dt:.4f}s  {status}")
            except Exception as ex:
                results[key] = (False, float('inf'), float('inf'), float('inf'), 0)
                print(f"  {name:45s} n={n:4d}  ERROR: {ex}")
        else:
            for sz in sizes:
                try:
                    d, e = generator(name, sz)
                except Exception as ex:
                    key = f"{name}@{sz}"
                    results[key] = (False, float('inf'), float('inf'), float('inf'), 0)
                    total += 1
                    print(f"  {name:45s} n={sz:4d}  ERROR: {ex}")
                    continue
                n = len(d)
                key = f"{name}@{n}"
                total += 1
                try:
                    ok, res, ou, ov, dt = test_one(d, e)
                    results[key] = (ok, res, ou, ov, dt)
                    if ok:
                        total_pass += 1
                    status = 'PASS' if ok else 'FAIL'
                    print(f"  {name:45s} n={n:4d}  res={res:8.3f}  ortU={ou:8.3f}  ortV={ov:8.3f}  t={dt:.4f}s  {status}")
                except Exception as ex:
                    results[key] = (False, float('inf'), float('inf'), float('inf'), 0)
                    print(f"  {name:45s} n={n:4d}  ERROR: {ex}")

    return results, total_pass, total


def main():
    parser = argparse.ArgumentParser(description='Dense-to-Bidiag + Paper Test Suite for MR3-GK')
    parser.add_argument('--quick', action='store_true', help='Quick mode: n=10,50 only')
    parser.add_argument('--part1', action='store_true', help='Run only PART 1 (dense-to-bidiag)')
    parser.add_argument('--part2', action='store_true', help='Run only PART 2 (paper tests)')
    parser.add_argument('--timeout', type=int, default=60, help='Per-test timeout in seconds')
    args = parser.parse_args()

    if args.quick:
        sizes = [10, 100]
    else:
        sizes = [10, 100, 200, 400]

    run_part1 = not args.part2 or args.part1
    run_part2 = not args.part1 or args.part2
    if not args.part1 and not args.part2:
        run_part1 = run_part2 = True

    all_results = {}
    grand_pass = 0
    grand_total = 0

    if run_part1:
        print("=" * 80)
        print("PART 1: Dense-to-Bidiagonal Tests")
        print(f"  {len(DENSE_TESTS)} patterns x {len(sizes)} sizes = {len(DENSE_TESTS) * len(sizes)} tests")
        print("=" * 80)
        r, p, t = run_tests(DENSE_TESTS, make_dense, sizes, set(), "Dense", args.timeout)
        all_results.update(r)
        grand_pass += p
        grand_total += t
        print(f"\n  PART 1 subtotal: {p}/{t} passed")

    if run_part2:
        print("\n" + "=" * 80)
        print("PART 2: Missing Paper Test Cases")
        n_fixed = len([t for t in PAPER_TESTS if t in FIXED_SIZE_TESTS])
        n_sized = len(PAPER_TESTS) - n_fixed
        total_expected = n_fixed + n_sized * len(sizes)
        print(f"  {len(PAPER_TESTS)} patterns ({n_fixed} fixed + {n_sized} x {len(sizes)} sizes) = {total_expected} tests")
        print("=" * 80)
        r, p, t = run_tests(PAPER_TESTS, make_paper, sizes, FIXED_SIZE_TESTS, "Paper", args.timeout)
        all_results.update(r)
        grand_pass += p
        grand_total += t
        print(f"\n  PART 2 subtotal: {p}/{t} passed")

    # Summary
    print("\n" + "=" * 80)
    print(f"GRAND TOTAL: {grand_pass}/{grand_total} passed ({100 * grand_pass / max(grand_total, 1):.1f}%)")
    failed = [(k, v) for k, v in sorted(all_results.items()) if not v[0]]
    if failed:
        print(f"\nFailing ({len(failed)}):")
        for key, (ok, res, ou, ov, dt) in failed:
            print(f"  {key:55s}  res={res:10.3f}  ortU={ou:10.3f}  ortV={ov:10.3f}")
    else:
        print("\nAll tests passed!")
    print("=" * 80)


if __name__ == '__main__':
    main()
