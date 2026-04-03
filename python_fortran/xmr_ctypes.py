"""
Python ctypes interface to the XMR block factorization library (libxmr.so).

Provides:
  xmr_build_repr  — Build block representation from G, OMEGA, E
  xmr_shift       — Block-aware shift transformation (dlaxrs)
  xmr_twist       — Block-aware twisted factorization (dlaxrt)
  xmr_eigvec      — Eigenvector extraction from twisted factorization
  xmr_singleton   — Full singleton computation (twist + eigvec + RQI)
"""

import ctypes
import numpy as np
import os

# Load the shared library
# Support override via _LIB_PATH_OVERRIDE for ablation studies
_libpath = os.environ.get('XMR_LIB_PATH',
    os.path.join(os.path.dirname(os.path.abspath(__file__)), 'libxmr.so'))
_lib = ctypes.CDLL(_libpath)

# C function signatures
_lib.xmr_init.restype = None
_lib.xmr_init.argtypes = []

_lib.xmr_build_repr.restype = None
_lib.xmr_build_repr.argtypes = [
    ctypes.c_int,                          # n
    ctypes.c_int,                          # k
    ctypes.POINTER(ctypes.c_double),       # g (1-based)
    ctypes.POINTER(ctypes.c_int),          # omega (0-based)
    ctypes.POINTER(ctypes.c_double),       # e (1-based)
    ctypes.c_double,                       # pivbase
    ctypes.POINTER(ctypes.c_double),       # repr (output)
    ctypes.POINTER(ctypes.c_int),          # repi (output)
]

_lib.xmr_shift.restype = None
_lib.xmr_shift.argtypes = [
    ctypes.c_int,                          # n
    ctypes.c_int,                          # ia
    ctypes.c_int,                          # ie
    ctypes.POINTER(ctypes.c_double),       # repr
    ctypes.POINTER(ctypes.c_int),          # repi
    ctypes.c_double,                       # tau
    ctypes.POINTER(ctypes.c_double),       # shfprt
    ctypes.POINTER(ctypes.c_double),       # dplus (output)
    ctypes.POINTER(ctypes.c_int),          # omgadp (output)
    ctypes.POINTER(ctypes.c_double),       # rplus (output)
    ctypes.POINTER(ctypes.c_int),          # omgarp (output)
    ctypes.POINTER(ctypes.c_double),       # gammap (output)
    ctypes.POINTER(ctypes.c_int),          # twistok (output)
    ctypes.POINTER(ctypes.c_double),       # rwork
]

_lib.xmr_twist.restype = None
_lib.xmr_twist.argtypes = [
    ctypes.c_int,                          # n
    ctypes.c_int,                          # j1
    ctypes.c_int,                          # j2
    ctypes.POINTER(ctypes.c_double),       # repr
    ctypes.POINTER(ctypes.c_int),          # repi
    ctypes.c_double,                       # tau
    ctypes.POINTER(ctypes.c_double),       # dplus (output)
    ctypes.POINTER(ctypes.c_double),       # rplus (output)
    ctypes.POINTER(ctypes.c_double),       # gammap (output)
    ctypes.POINTER(ctypes.c_int),          # twistok (output)
    ctypes.POINTER(ctypes.c_double),       # rwork
]

# Initialize the COMMON block once
_lib.xmr_init()

EPS = np.finfo(np.float64).eps
SAFMIN = np.finfo(np.float64).tiny
BIGG = np.finfo(np.float64).max


def _ptr(arr):
    """Get ctypes pointer to numpy array data."""
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def _iptr(arr):
    """Get ctypes pointer to numpy int32 array data."""
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int))


class XMRRepr:
    """XMR representation: wraps the Fortran REPR/REPI arrays."""

    def __init__(self, n, k=None):
        self.n = n
        if k is None:
            k = n
        self.k = k
        self.repr = np.zeros(4 * n + 3, dtype=np.float64)
        self.repi = np.zeros(6 + n + n // 2, dtype=np.int32)
        self.e = np.zeros(n, dtype=np.float64)  # 1-based off-diags
        self.shfprt = np.ones(n, dtype=np.float64)

    @property
    def nb(self):
        return self.repi[2]  # Fortran REPI(3) = C index 2


def build_repr_from_tridiag(c, e_off, n, k=None):
    """Build XMR block representation from tridiagonal (c, e_off).

    c: 0-based diagonal [0:n-1]
    e_off: 0-based off-diagonal [0:n-2]
    k: twist index (1-based). Default: n.

    Returns XMRRepr.
    """
    if k is None:
        k = n

    rep = XMRRepr(n, k)

    # Compute LDL^T pivots with block creation at small pivots
    g = np.zeros(n + 1, dtype=np.float64)  # 1-based
    omega = np.zeros(n + 2, dtype=np.int32)  # 0-based

    # e in 1-based format
    e_1b = np.zeros(n, dtype=np.float64)
    for i in range(n - 1):
        e_1b[i + 1] = e_off[i]
    rep.e[:] = e_1b

    # Top-down LDL^T with block creation (for positions 1..k)
    KBLOCK = 1.0 / 8.0
    d = np.zeros(n + 1, dtype=np.float64)
    d[1] = c[0]
    i = 1
    while i < k:
        if i + 1 > n:
            break
        esq = e_1b[i] ** 2 if i < n else 0.0
        if abs(d[i]) < SAFMIN:
            omega[i + 1] = 1
            d[i + 1] = c[i]
            bdet = d[i] * c[i] - esq
            i += 2
            if i <= k and i <= n:
                esq_next = e_1b[i - 1] ** 2 if i - 1 < n else 0.0
                d[i] = c[i - 1] - esq_next * d[i - 2] / bdet if abs(bdet) > SAFMIN else c[i - 1]
        elif esq > 0 and abs(d[i] * (c[i] - esq / d[i])) < KBLOCK * esq:
            omega[i + 1] = 1
            d[i + 1] = c[i]
            bdet = d[i] * c[i] - esq
            i += 2
            if i <= k and i <= n:
                esq_next = e_1b[i - 1] ** 2 if i - 1 < n else 0.0
                d[i] = c[i - 1] - esq_next * d[i - 2] / bdet if abs(bdet) > SAFMIN else c[i - 1]
        else:
            d[i + 1] = c[i] - esq / d[i] if abs(d[i]) > SAFMIN else c[i]
            i += 1

    # Bottom-up for positions n..k+1
    r = np.zeros(n + 1, dtype=np.float64)
    r[n] = c[n - 1]
    i = n
    while i > k:
        if i - 1 < 1:
            break
        esq = e_1b[i - 1] ** 2 if i - 1 >= 1 else 0.0
        if abs(r[i]) < SAFMIN:
            omega[i - 1] = 1
            r[i - 1] = c[i - 2]
            bdet = r[i] * c[i - 2] - esq
            i -= 2
            if i >= k and i >= 1:
                esq_prev = e_1b[i] ** 2 if i >= 1 else 0.0
                r[i] = c[i - 1] - esq_prev * r[i + 2] / bdet if abs(bdet) > SAFMIN else c[i - 1]
        elif esq > 0 and abs(r[i] * (c[i - 2] - esq / r[i])) < KBLOCK * esq:
            omega[i - 1] = 1
            r[i - 1] = c[i - 2]
            bdet = r[i] * c[i - 2] - esq
            i -= 2
            if i >= k and i >= 1:
                esq_prev = e_1b[i] ** 2 if i >= 1 else 0.0
                r[i] = c[i - 1] - esq_prev * r[i + 2] / bdet if abs(bdet) > SAFMIN else c[i - 1]
        else:
            r[i - 1] = c[i - 2] - esq / r[i] if abs(r[i]) > SAFMIN else c[i - 2]
            i -= 1

    # Merge
    for i in range(1, k + 1):
        g[i] = d[i]
    for i in range(k + 1, n + 1):
        g[i] = r[i]

    # Call Fortran dlaxrr to finalize REPR/REPI
    _lib.xmr_build_repr(n, k, _ptr(g), _iptr(omega), _ptr(e_1b),
                         0.0, _ptr(rep.repr), _iptr(rep.repi))
    return rep


def xmr_twist_and_eigvec(rep, tau):
    """Compute eigenvector at eigenvalue tau using block twisted factorization.

    rep: XMRRepr from build_repr_from_tridiag
    tau: eigenvalue approximation

    Returns (gamma_best, twist_idx, z) where z is 0-based [0:n-1].
    """
    n = rep.n

    dplus = np.zeros(n, dtype=np.float64)
    rplus = np.zeros(n, dtype=np.float64)
    gammap = np.zeros(n, dtype=np.float64)
    twistok = np.zeros(n, dtype=np.int32)
    rwork = np.zeros(2 * n, dtype=np.float64)

    _lib.xmr_twist(n, 1, n, _ptr(rep.repr), _iptr(rep.repi),
                    tau, _ptr(dplus), _ptr(rplus),
                    _ptr(gammap), _iptr(twistok), _ptr(rwork))

    # Find best twist (1-based indices stored in 0-based arrays)
    # gammap[i-1] = gamma at position i (Fortran 1-based → C 0-based)
    best_r = -1
    best_g = BIGG
    for i in range(n):
        if twistok[i] == 1 and abs(gammap[i]) < best_g:
            best_g = abs(gammap[i])
            best_r = i + 1  # convert to 1-based

    if best_r < 0:
        best_r = n  # fallback to twist position

    # Extract eigenvector using DPLUS/RPLUS and e_off
    # z[r] = 1
    # z[i] = -(e[i] / DPLUS[i]) * z[i+1]   for i = r-1, ..., 1  (1-based)
    # z[i] = -(e[i-1] / RPLUS[i]) * z[i-1]  for i = r+1, ..., N  (1-based)
    #
    # In 0-based arrays:
    #   DPLUS[i-1] corresponds to Fortran DPLUS(i) for i=1..N-1
    #   RPLUS[i-1] corresponds to Fortran RPLUS(i+1) for i=0..N-2
    #   e_1b[i] is 1-based: e_1b[i] = e_off[i-1]

    z = np.zeros(n + 1, dtype=np.float64)  # 1-based
    z[best_r] = 1.0

    e_1b = rep.e

    # Top-down: r-1 down to 1
    for i in range(best_r - 1, 0, -1):
        # DPLUS(i) in Fortran = dplus[i-1] in C (0-based)
        dp = dplus[i - 1]
        if abs(dp) > SAFMIN:
            z[i] = -(e_1b[i] / dp) * z[i + 1]
        else:
            z[i] = 0.0

    # Bottom-up: r+1 up to N
    for i in range(best_r + 1, n + 1):
        # RPLUS(i) in Fortran = rplus[i-2] in C (since RPLUS declared as (2:N))
        rp = rplus[i - 2]
        if abs(rp) > SAFMIN:
            z[i] = -(e_1b[i - 1] / rp) * z[i - 1]
        else:
            z[i] = 0.0

    z_out = z[1:n + 1].copy()
    nrm = np.linalg.norm(z_out)
    if nrm > SAFMIN:
        z_out /= nrm

    return best_g, best_r, z_out


def xmr_singleton(c, e_off, n, lam, rep=None, max_rqi=5):
    """Compute eigenpair at eigenvalue lam using XMR block factorizations.

    c: 0-based diagonal [0:n-1] of the current node's tridiagonal
    e_off: 0-based off-diagonal [0:n-2]
    n: matrix size
    lam: eigenvalue approximation
    rep: XMRRepr (if None, builds from c, e_off)
    max_rqi: max Rayleigh quotient iterations

    Returns (lam_refined, z) where z is 0-based [0:n-1].
    """
    if rep is None:
        rep = build_repr_from_tridiag(c, e_off, n)

    for rqi_iter in range(max_rqi):
        gamma, twist, z = xmr_twist_and_eigvec(rep, lam)

        # Rayleigh quotient
        tz = c * z
        if n > 1:
            tz[:n - 1] += e_off * z[1:]
            tz[1:] += e_off * z[:n - 1]
        rq = np.dot(z, tz)
        res_norm = np.linalg.norm(tz - rq * z)

        if res_norm < 10 * n * EPS * (abs(rq) + SAFMIN):
            lam = rq
            break
        lam = rq

    return lam, z


def xmr_shift_repr(rep, tau):
    """Compute shifted representation using block-aware dlaxrs.

    Returns (DPLUS, RPLUS, OMGADP, OMGARP, GAMMAP, TWISTOK).
    All arrays are 0-based (Fortran 1-based mapped to C 0-based).
    """
    n = rep.n
    dplus = np.zeros(n, dtype=np.float64)
    rplus = np.zeros(n, dtype=np.float64)
    omgadp = np.zeros(n, dtype=np.int32)
    omgarp = np.zeros(n, dtype=np.int32)
    gammap = np.zeros(n, dtype=np.float64)
    twistok = np.zeros(n, dtype=np.int32)
    rwork = np.zeros(2 * n, dtype=np.float64)

    _lib.xmr_shift(n, 1, n, _ptr(rep.repr), _iptr(rep.repi),
                    tau, _ptr(rep.shfprt),
                    _ptr(dplus), _iptr(omgadp), _ptr(rplus), _iptr(omgarp),
                    _ptr(gammap), _iptr(twistok), _ptr(rwork))

    return dplus, rplus, omgadp, omgarp, gammap, twistok


# Register xmr_negcount
_lib.xmr_negcount.restype = ctypes.c_int
_lib.xmr_negcount.argtypes = [
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_double,
]


def xmr_negcount(rep, tau):
    """Compute inertia (negative pivot count) using block-aware representation.

    Returns xi = 2*negcount + issing.
    xi/2 gives the number of eigenvalues strictly below tau.
    """
    return _lib.xmr_negcount(rep.n, _ptr(rep.repr), _iptr(rep.repi), tau)


def xmr_refine_eigenvalues(rep, approxs, rel_tol=None, max_iter=100):
    """Refine eigenvalue approximations to relative precision using dlaxrn.

    Uses bisection with the Repr-based inertia counts and a relative
    stopping criterion: width <= rel_tol * |midpoint|.

    rep: XMRRepr
    approxs: array of eigenvalue approximations (0-based)
    rel_tol: relative tolerance (default 4*n*eps)

    Returns refined eigenvalue array.
    """
    n = rep.n
    if rel_tol is None:
        rel_tol = 4 * n * EPS

    refined = np.zeros(len(approxs), dtype=np.float64)

    for j, approx in enumerate(approxs):
        # Get target inertia: xi just above approx
        xi_at = xmr_negcount(rep, approx)
        # Target eigenvalue index: the one approx is closest to
        # xi = 2*negc + issing. If xi is odd, we're ON an eigenvalue.
        # If xi is even, we're between ew xi/2 and ew xi/2+1.
        # The target eigenvalue is the one we're closest to.
        target_ew = (xi_at + 1) // 2  # 1-based eigenvalue index

        # Find bracket
        lo = approx
        hi = approx
        step = max(abs(approx) * 1e-3, 100 * EPS)

        for _ in range(60):
            lo = approx - step
            xi_lo = xmr_negcount(rep, lo)
            if xi_lo < 2 * target_ew - 1:
                break
            step *= 2

        step = max(abs(approx) * 1e-3, 100 * EPS)
        for _ in range(60):
            hi = approx + step
            xi_hi = xmr_negcount(rep, hi)
            if xi_hi >= 2 * target_ew - 1:
                break
            step *= 2

        # Bisect with relative stopping criterion
        for it in range(max_iter):
            mid = (lo + hi) * 0.5
            width = hi - lo
            if width <= rel_tol * max(abs(mid), EPS) or mid == lo or mid == hi:
                break
            xi_mid = xmr_negcount(rep, mid)
            if xi_mid < 2 * target_ew - 1:
                lo = mid
            else:
                hi = mid

        refined[j] = (lo + hi) * 0.5

    return refined


# Register xmr_find_child_repr
_lib.xmr_find_child_repr.restype = ctypes.c_int
_lib.xmr_find_child_repr.argtypes = [
    ctypes.c_int,                          # n
    ctypes.POINTER(ctypes.c_double),       # e (1-based)
    ctypes.POINTER(ctypes.c_double),       # father_repr
    ctypes.POINTER(ctypes.c_int),          # father_repi
    ctypes.c_int,                          # icbeg
    ctypes.c_int,                          # icend
    ctypes.POINTER(ctypes.c_double),       # evals
    ctypes.c_double,                       # lgap
    ctypes.c_double,                       # ugap
    ctypes.c_double,                       # gaptol
    ctypes.c_double,                       # spdiam
    ctypes.c_double,                       # taubar
    ctypes.POINTER(ctypes.c_double),       # son_repr (output)
    ctypes.POINTER(ctypes.c_int),          # son_repi (output)
    ctypes.POINTER(ctypes.c_double),       # tau_out (output)
]


def xmr_find_child(father_rep, e_off, cluster_evals, icbeg, icend,
                    lgap, ugap, gaptol=1e-3, spdiam=None, taubar=0.0):
    """Find child representation for a cluster using XMR's dlaxrf.

    father_rep: XMRRepr of the parent
    e_off: 0-based off-diagonal array [0:n-2]
    cluster_evals: eigenvalue approximations for the cluster
    icbeg, icend: 1-based eigenvalue index range in the parent
    lgap, ugap: gaps to left/right of the cluster
    gaptol: gap tolerance
    spdiam: spectral diameter (default: 2*max(|evals|))
    taubar: accumulated shift

    Returns (child_repr, tau, info) where child_repr is an XMRRepr.
    """
    n = father_rep.n
    if spdiam is None:
        spdiam = 2.0 * max(abs(cluster_evals[0]), abs(cluster_evals[-1]), SAFMIN)

    # e in 1-based format
    e_1b = np.zeros(n, dtype=np.float64)
    for i in range(n - 1):
        e_1b[i + 1] = e_off[i]

    evals_arr = np.array(cluster_evals, dtype=np.float64)
    son_repr = np.zeros(4 * n + 3, dtype=np.float64)
    son_repi = np.zeros(6 + n + n // 2, dtype=np.int32)
    tau_out = np.zeros(1, dtype=np.float64)

    info = _lib.xmr_find_child_repr(
        n, _ptr(e_1b), _ptr(father_rep.repr), _iptr(father_rep.repi),
        icbeg, icend, _ptr(evals_arr),
        lgap, ugap, gaptol, spdiam, taubar,
        _ptr(son_repr), _iptr(son_repi), _ptr(tau_out))

    child = XMRRepr(n)
    child.repr[:] = son_repr
    child.repi[:] = son_repi
    child.e[:] = e_1b

    return child, tau_out[0], info


# Register xmr_eigenvectors
_lib.xmr_eigenvectors.restype = ctypes.c_int
_lib.xmr_eigenvectors.argtypes = [
    ctypes.c_int,                          # n
    ctypes.POINTER(ctypes.c_double),       # e_1based
    ctypes.POINTER(ctypes.c_double),       # rootr
    ctypes.POINTER(ctypes.c_int),          # rooti
    ctypes.POINTER(ctypes.c_double),       # evals_in (all n eigenvalues)
    ctypes.c_int,                          # wil
    ctypes.c_int,                          # wiu
    ctypes.c_double,                       # spdiam
    ctypes.c_double,                       # gaptol
    ctypes.POINTER(ctypes.c_double),       # w_out
    ctypes.POINTER(ctypes.c_double),       # z_out
]


def xmr_eigenvectors(root_repr, e_off, all_evals, wil, wiu,
                      spdiam=None, gaptol=1e-3):
    """Compute eigenvectors wil..wiu using full XMR pipeline (dlaxrv).

    root_repr: XMRRepr of the root representation
    e_off: 0-based off-diagonals [0..n-2]
    all_evals: ALL n eigenvalue approximations (0-based, ascending)
    wil, wiu: wanted eigenvalue range (1-based)
    spdiam: spectral diameter
    gaptol: gap tolerance

    Returns (evals, evecs, info) where:
      evals: eigenvalues [0..wiu-wil] (0-based)
      evecs: eigenvectors [n, wiu-wil+1] column-major (0-based)
    """
    n = root_repr.n
    nwant = wiu - wil + 1
    if spdiam is None:
        spdiam = all_evals[-1] - all_evals[0]

    e_1b = np.zeros(n, dtype=np.float64)
    for i in range(n - 1):
        e_1b[i + 1] = e_off[i]

    evals_arr = np.array(all_evals, dtype=np.float64)
    w_out = np.zeros(nwant, dtype=np.float64)
    z_out = np.zeros(n * nwant, dtype=np.float64)

    info = _lib.xmr_eigenvectors(
        n, _ptr(e_1b), _ptr(root_repr.repr), _iptr(root_repr.repi),
        _ptr(evals_arr), wil, wiu, spdiam, gaptol,
        _ptr(w_out), _ptr(z_out))

    evecs = z_out.reshape((n, nwant), order='F')
    return w_out, evecs, info


def xmr_eigenvectors_v(root_repr, e_off, all_evals, wil, wiu,
                        spdiam=None, gaptol=1e-3):
    """Like xmr_eigenvectors but also returns the dlaxre shift tau_re."""
    n = root_repr.n
    nwant = wiu - wil + 1
    if spdiam is None:
        spdiam = all_evals[-1] - all_evals[0]

    # Register the C function if not already done
    if not hasattr(_lib, '_xmr_eigvec_v_registered'):
        _lib.xmr_eigenvectors_v.restype = ctypes.c_int
        _lib.xmr_eigenvectors_v.argtypes = [
            ctypes.c_int,                          # n
            ctypes.POINTER(ctypes.c_double),       # e_1based
            ctypes.POINTER(ctypes.c_double),       # rootr
            ctypes.POINTER(ctypes.c_int),          # rooti
            ctypes.POINTER(ctypes.c_double),       # evals_in
            ctypes.c_int,                          # wil
            ctypes.c_int,                          # wiu
            ctypes.c_double,                       # spdiam
            ctypes.c_double,                       # gaptol
            ctypes.POINTER(ctypes.c_double),       # w_out
            ctypes.POINTER(ctypes.c_double),       # z_out
            ctypes.POINTER(ctypes.c_double),       # tau_re_out
        ]
        _lib._xmr_eigvec_v_registered = True

    e_1b = np.zeros(n, dtype=np.float64)
    for i in range(n - 1):
        e_1b[i + 1] = e_off[i]

    evals_arr = np.array(all_evals, dtype=np.float64)
    w_out = np.zeros(nwant, dtype=np.float64)
    z_out = np.zeros(n * nwant, dtype=np.float64)
    tau_re_out = np.zeros(1, dtype=np.float64)

    info = _lib.xmr_eigenvectors_v(
        n, _ptr(e_1b), _ptr(root_repr.repr), _iptr(root_repr.repi),
        _ptr(evals_arr), wil, wiu, spdiam, gaptol,
        _ptr(w_out), _ptr(z_out), _ptr(tau_re_out))

    evecs = z_out.reshape((n, nwant), order='F')
    return w_out, evecs, info, tau_re_out[0]
