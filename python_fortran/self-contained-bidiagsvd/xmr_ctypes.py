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
