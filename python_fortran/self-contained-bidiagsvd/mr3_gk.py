"""
MR3-GK: Algorithm 4.1 (Willems & Lang, ETNA 2012)

Implements from paper and XMR:
- Geometric shift spacing, inside shifts, last resort midpoint
- EVAL = DELG + DRCOND for shift ranking
- Inertia check, COB, fudging
- Bidirectional one-sided recovery
- Block LDL^T factorizations (ETNA 2011)
- Full 5-requirement metrics
- No reorthogonalization, no fallbacks
"""

import numpy as np
from scipy.linalg.lapack import _flapack
import ctypes as _ct
import ctypes.util as _ctutil
# LAPACK DLARNV used for reproducible random vectors in GS completion.
# Using the same call (and same fixed seed) on the Fortran side guarantees
# bit-identical results.
_lapack_lib = _ct.CDLL(_ctutil.find_library('lapack'))

def _dlarnv_normal(n, seed):
    """Generate n N(0,1) samples via LAPACK DLARNV (IDIST=3).
    seed: tuple of 4 ints (must be in 0..4095, last must be odd)."""
    iseed = (_ct.c_int * 4)(*seed)
    n_c = _ct.c_int(n)
    idist = _ct.c_int(3)
    x = np.zeros(n, dtype=np.float64)
    _lapack_lib.dlarnv_(_ct.byref(idist), iseed, _ct.byref(n_c),
                        x.ctypes.data_as(_ct.POINTER(_ct.c_double)))
    return x

# Unified DNRM2 via system BLAS (bit-matches Fortran's dnrm2 call).
_blas_lib = _ct.CDLL(_ctutil.find_library('blas') or _ctutil.find_library('lapack'))
_blas_lib.dnrm2_.restype = _ct.c_double
_blas_lib.dnrm2_.argtypes = [_ct.POINTER(_ct.c_int),
                              _ct.POINTER(_ct.c_double),
                              _ct.POINTER(_ct.c_int)]

def _system_dnrm2(x):
    """Vector 2-norm via system BLAS DNRM2 — matches the Fortran side bit-for-bit."""
    x_arr = np.ascontiguousarray(x, dtype=np.float64)
    n_c = _ct.c_int(x_arr.size)
    inc = _ct.c_int(1)
    return _blas_lib.dnrm2_(_ct.byref(n_c),
                            x_arr.ctypes.data_as(_ct.POINTER(_ct.c_double)),
                            _ct.byref(inc))

def _system_dnrm2_axis0(M):
    """Column-wise DNRM2 of a 2D array, matching Fortran's per-column dnrm2 calls."""
    M_arr = np.ascontiguousarray(M, dtype=np.float64)
    rows, cols = M_arr.shape
    out = np.zeros(cols, dtype=np.float64)
    for j in range(cols):
        # Column j: contiguous memory in column-major (Fortran-order),
        # but our M may be row-major; use vector slice.
        col = M_arr[:, j]
        col_c = np.ascontiguousarray(col)
        n_c = _ct.c_int(rows)
        inc = _ct.c_int(1)
        out[j] = _blas_lib.dnrm2_(_ct.byref(n_c),
                                   col_c.ctypes.data_as(_ct.POINTER(_ct.c_double)),
                                   _ct.byref(inc))
    return out

EPS = np.finfo(np.float64).eps
SAFMIN = np.finfo(np.float64).tiny
GAPTOL = 1e-3
MAX_DEPTH = 50
ELG_THRESH = 8
MAXRELCOND = 10

# Verbose mode
_verbose = False
def set_verbose(v=True):
    global _verbose
    _verbose = v

def zero_shift_qr_sweep(d, e):
    """One downward sweep of implicit zero-shift QR on upper bidiagonal.
    Modifies d, e in place. Returns (right_rots, left_rots).
    Zero-shift QR = inverse iteration on B^TB (Demmel-Kahan 1990, Sec 3).
    When sigma_min = 0, converges in one sweep."""
    from scipy.linalg.lapack import dlartg
    n = len(d)
    if n <= 1:
        return [], []
    right_rots = []; left_rots = []
    cs = 1.0; oldcs = 1.0; oldsn = 0.0
    for i in range(n - 1):
        d1 = d[i] * cs
        cs, sn, r = dlartg(d1, e[i])
        right_rots.append((cs, sn, i))
        if i > 0:
            e[i-1] = oldsn * r
        d1 = oldcs * r
        d2 = d[i+1] * sn
        oldcs, oldsn, d[i] = dlartg(d1, d2)
        left_rots.append((oldcs, oldsn, i))
    h = d[n-1] * cs
    d[n-1] = h * oldcs
    e[n-2] = h * oldsn
    return right_rots, left_rots

def _apply_givens_to_rows(M, rots, row_offset=0):
    """Apply G_0 @ G_1 @ ... @ G_{m-1} to M from the left (reversed order).
    G_i = [cs, -sn; sn, cs] at rows (row_offset+i, row_offset+i+1)."""
    for cs, sn, i in reversed(rots):
        ri = row_offset + i
        r0 = M[ri, :].copy()
        r1 = M[ri+1, :].copy()
        M[ri, :]   = cs * r0 - sn * r1
        M[ri+1, :] = sn * r0 + cs * r1

# ============================================================
# Bisection
# ============================================================
def bisect_evals(c, e_off, n, il, iu, tol=0.0):
    if n <= 0 or il > iu: return np.array([])
    if n == 1: return np.array([c[0]])
    try:
        m, w, _, _, info = _flapack.dstebz(
            c[:n].astype(np.float64), e_off[:n-1].astype(np.float64),
            2, 0.0, 0.0, int(il), int(iu), tol, b'E')
        if info == 0 and m > 0: return w[:m].copy()
    except Exception: pass
    return np.array([])

# ============================================================
# Splitting
# ============================================================
def split_bidiag(d, e, n):
    if n <= 1: return [(0, n-1)]

    # Phase 1: Relative split (preserves relative accuracy)
    rel_blocks = []; beg = 0
    for i in range(n-1):
        rel_thresh = EPS*(abs(d[i])+abs(d[i+1]))
        if abs(e[i]) <= max(rel_thresh, SAFMIN):
            rel_blocks.append((beg, i)); beg = i+1
    rel_blocks.append((beg, n-1))

    # Phase 2: Absolute split within each relative block per condition (3.5)
    blocks = []
    for rb_beg, rb_end in rel_blocks:
        rb_k = rb_end - rb_beg + 1
        if rb_k <= 1:
            blocks.append((rb_beg, rb_end))
            continue
        # Compute ||B_sub||_∞ for this relative block
        bnorm = 0.0
        for i in range(rb_beg, rb_end + 1):
            s = abs(d[i])
            if i < rb_end: s += abs(e[i])
            bnorm = max(bnorm, s)
        abs_thresh = rb_k * EPS * bnorm  # condition (3.5): nε||B_sub||

        sub_beg = rb_beg
        for i in range(rb_beg, rb_end):
            if abs(e[i]) <= max(abs_thresh, SAFMIN):
                blocks.append((sub_beg, i)); sub_beg = i+1
        blocks.append((sub_beg, rb_end))

    return blocks

def _solve_tgk_block(bc, bte, m, bk_hint=None):
    """Solve a single T_GK eigenproblem of size m via XMR.
    Returns (w, Z) — non-negative eigenvalues and the corresponding
    eigenvectors of T_GK = perfect_shuffle of (bc, bte).
    """
    from xmr_ctypes import xmr_eigenvectors, xmr_eigenvectors_v, build_repr_from_tridiag
    n_nonneg = (m + 1) // 2
    if bk_hint is not None:
        n_nonneg = bk_hint
    wil = m - n_nonneg + 1
    wiu = m

    all_evals = np.sort(bisect_evals(bc, bte, m, 1, m))
    spdiam = all_evals[-1] - all_evals[0] if len(all_evals) > 0 else 0.0

    if _verbose:
        print(f'  [_solve_tgk_block] m={m}, n_nonneg={n_nonneg}, wil={wil}, wiu={wiu}')
        print(f'  [_solve_tgk_block] spdiam={spdiam:.10e}')
        print(f'  [_solve_tgk_block] evals range=[{all_evals[0]:.10e}, {all_evals[-1]:.10e}]')

    root = build_repr_from_tridiag(bc, bte, m, k=m)
    block_gaptol = max(GAPTOL, 0.02 / n_nonneg)

    if _verbose:
        w, Z, info, tau_re = xmr_eigenvectors_v(root, bte, all_evals,
                                                  wil=wil, wiu=wiu,
                                                  spdiam=spdiam, gaptol=block_gaptol)
        print(f'  [_solve_tgk_block] dlaxre shift tau_re={tau_re:.10e}')
        print(f'  [_solve_tgk_block] XMR info={info}')
        if info == 0:
            print(f'  [_solve_tgk_block] XMR output: {len(w)} evals, range=[{w[0]:.10e}, {w[-1]:.10e}]')
    else:
        w, Z, info = xmr_eigenvectors(root, bte, all_evals,
                                        wil=wil, wiu=wiu,
                                        spdiam=spdiam, gaptol=block_gaptol)

    if info != 0:
        raise RuntimeError(f"XMR returned info={info} on T_GK block of size {m}")

    return w[:n_nonneg], Z[:m, :n_nonneg]


def mr3_tgk_multiblock(d_bidiag, e_bidiag, n, bbeg, bend):
    """Run MR³ on a single multi-element block. Returns (evals, evecs) of size (k, 2k×k)."""
    from xmr_ctypes import xmr_eigenvectors, build_repr_from_tridiag
    bk = bend - bbeg + 1
    m2k = 2 * bk
    bd_abs = np.abs(d_bidiag[bbeg:bend+1])
    be = np.abs(e_bidiag[bbeg:bend])
    bc = np.zeros(m2k)
    bte = np.zeros(m2k - 1)
    bte[0::2] = bd_abs
    if bk > 1:
        bte[1::2] = be[:bk-1]

    if _verbose:
        print(f'\n[mr3_tgk_multiblock] block [{bbeg},{bend}] bk={bk} m2k={m2k}')
        if bk <= 20:
            print(f'  input bidiag: d[{bbeg}:{bend+1}] = {np.array2string(d_bidiag[bbeg:bend+1], precision=10, separator=", ")}')
            print(f'  input bidiag: e[{bbeg}:{bend}]   = {np.array2string(e_bidiag[bbeg:bend], precision=10, separator=", ")}')
        else:
            print(f'  input bidiag: d range=[{np.min(np.abs(d_bidiag[bbeg:bend+1])):.3e}, {np.max(np.abs(d_bidiag[bbeg:bend+1])):.3e}]')
            print(f'  input bidiag: e range=[{np.min(np.abs(e_bidiag[bbeg:bend])):.3e}, {np.max(np.abs(e_bidiag[bbeg:bend])):.3e}]')
        print(f'  T_GK off-diagonal bte ({len(bte)} entries), range [{np.min(bte):.3e}, {np.max(bte):.3e}]')

    w, Z = _solve_tgk_block(bc, bte, m2k, bk_hint=bk)

    if _verbose:
        print(f'  [mr3_tgk_multiblock] result: {bk} eigenvalues, range [{w[0]:.10e}, {w[min(bk-1,len(w)-1)]:.10e}]')
        even_norms = np.array([np.linalg.norm(Z[0::2, j]) for j in range(bk)])
        odd_norms = np.array([np.linalg.norm(Z[1::2, j]) for j in range(bk)])
        max_norm_diff = np.max(np.abs(even_norms - odd_norms))
        print(f'    GK structure: max|even_norm - odd_norm| = {max_norm_diff:.3e}')

    return w[:bk], Z[:m2k, :bk]

    # =========================================================
    # T_GK has zeros in off-diagonal — split into sub-blocks
    # =========================================================
    # bte[zp]=0 means T_GK[zp, zp+1] = 0, so rows [start..zp] are decoupled
    # from rows [zp+1..end].
    boundaries = []
    start = 0
    for zp in zero_positions:
        boundaries.append((start, zp + 1))  # [start, zp+1) = rows start..zp
        start = zp + 1
    if start < m2k:
        boundaries.append((start, m2k))

    print(f'  TGK_SPLIT: block [{bbeg},{bend}] bk={bk} m2k={m2k}')
    print(f'  TGK_SPLIT: bte zeros at positions {zero_positions.tolist()}')
    print(f'  TGK_SPLIT: {len(boundaries)} sub-blocks:')
    for i, (sb_s, sb_e) in enumerate(boundaries):
        sb_m = sb_e - sb_s
        sb_nonneg = (sb_m + 1) // 2
        sb_bte = bte[sb_s:sb_e-1] if sb_m > 1 else np.array([])
        print(f'    sub[{i}]: T_GK rows [{sb_s},{sb_e}) size={sb_m} ceil(m/2)={sb_nonneg}')
        print(f'      bte = {np.round(sb_bte, 6)}')
        # Show which bidiag components this maps to
        # T_GK row 2i -> v_{i+1}, row 2i+1 -> u_{i+1}
        row_labels = []
        for r in range(sb_s, sb_e):
            if r % 2 == 0:
                row_labels.append(f'v{r//2}')
            else:
                row_labels.append(f'u{r//2}')
        print(f'      maps to: {row_labels}')

    # Solve each sub-block independently
    all_w = []
    all_Z_entries = []  # list of (sub_start, sub_end, w_sub, Z_sub)

    for i, (sb_start, sb_end) in enumerate(boundaries):
        sb_m = sb_end - sb_start
        if sb_m == 0:
            continue

        sb_bc = np.zeros(sb_m)
        sb_bte = bte[sb_start:sb_end-1].copy() if sb_m > 1 else np.array([])

        if sb_m == 1:
            # Size-1 sub-block: single eigenvalue = 0
            w_sub = np.array([0.0])
            Z_sub = np.ones((1, 1))
            print(f'    sub[{i}]: size=1, eigenvalue=0 (trivial)')
        else:
            w_sub, Z_sub = _solve_tgk_block(sb_bc, sb_bte, sb_m)
            n_ev = len(w_sub)
            print(f'    sub[{i}]: size={sb_m}, solved -> {n_ev} eigenpairs')
            print(f'      eigenvalues: {np.round(w_sub, 8)}')
            # Check GK structure for each eigenvector
            for j in range(n_ev):
                z = np.zeros(m2k)
                z[sb_start:sb_end] = Z_sub[:, j]
                even_norm = np.linalg.norm(z[0::2])
                odd_norm = np.linalg.norm(z[1::2])
                print(f'      evec[{j}]: eval={w_sub[j]:.8e}  even_norm={even_norm:.6f}  odd_norm={odd_norm:.6f}')

        all_Z_entries.append((sb_start, sb_end, w_sub, Z_sub))

    # Collect all eigenvalues and their associated eigenvectors
    # Total non-negative eigenvalues across sub-blocks may exceed bk
    # (each odd sub-block contributes an extra zero eigenvalue)
    # We want the bk largest eigenvalues
    total_evals = []
    for sb_start, sb_end, w_sub, Z_sub in all_Z_entries:
        for j in range(len(w_sub)):
            total_evals.append((w_sub[j], sb_start, sb_end, Z_sub[:, j]))

    print(f'  TGK_SPLIT: total eigenpairs collected = {len(total_evals)}, need bk = {bk}')

    # Sort by eigenvalue descending, take top bk
    total_evals.sort(key=lambda x: -x[0])
    selected = total_evals[:bk]
    # Sort selected back by eigenvalue ascending (standard ordering)
    selected.sort(key=lambda x: x[0])

    w_out = np.zeros(bk)
    Z_out = np.zeros((m2k, bk))
    for j, (ev, sb_start, sb_end, z_sub) in enumerate(selected):
        w_out[j] = ev
        Z_out[sb_start:sb_end, j] = z_sub

    print(f'  TGK_SPLIT: selected {bk} eigenvalues: {np.round(w_out, 8)}')
    # Verify: show how many zeros we kept vs discarded
    n_zero_kept = np.sum(np.abs(w_out) < SAFMIN)
    n_zero_total = sum(1 for ev, _, _, _ in total_evals if abs(ev) < SAFMIN)
    print(f'  TGK_SPLIT: zeros: {n_zero_kept} kept / {n_zero_total} total')

    return w_out, Z_out

# ============================================================
# Main entry
# ============================================================
def _bidiag_matvec_batch(d, e, M):
    """Compute B @ M where B = diag(d) + superdiag(e), batch over columns."""
    BM = d[:, None] * M
    if len(e) > 0:
        BM[:-1] += e[:, None] * M[1:]
    return BM

def _bidiagT_matvec_batch(d, e, M):
    """Compute B^T @ M where B = diag(d) + superdiag(e), batch over columns."""
    BTM = d[:, None] * M
    if len(e) > 0:
        BTM[1:] += e[:, None] * M[:-1]
    return BTM

def bidiag_svd(d, e):
    d = np.asarray(d, dtype=np.float64).copy()
    e = np.asarray(e, dtype=np.float64).copy()
    n = len(d)
    if n <= 0: return np.array([]), np.zeros((0,0)), np.zeros((0,0)), 0
    if n == 1:
        return (np.array([abs(d[0])]),
                np.array([[1.0 if d[0]>=0 else -1.0]]),
                np.array([[1.0]]), 0)

    if _verbose:
        print(f'\n{"="*60}')
        print(f'[bidiag_svd] n={n}')
        if n <= 20:
            print(f'  d = {np.array2string(d, precision=10, separator=", ")}')
            print(f'  e = {np.array2string(e, precision=10, separator=", ")}')
        else:
            print(f'  d range=[{np.min(d):.3e}, {np.max(d):.3e}], e range=[{np.min(e):.3e}, {np.max(e):.3e}]')

    # Pre-scaling: scale so max entry ~ 1. Doesn't affect U,V since
    # B = U·Σ·V^T  =>  (αB) = U·(αΣ)·V^T
    bnorm_scale = max(np.max(np.abs(d)), np.max(np.abs(e)) if n > 1 else 0.0)
    if bnorm_scale > 0.0:
        scale_factor = 1.0 / bnorm_scale
        d *= scale_factor
        e *= scale_factor
    else:
        scale_factor = 1.0

    if _verbose:
        print(f'  bnorm_scale={bnorm_scale:.10e}, scale_factor={scale_factor:.10e}')
        if scale_factor != 1.0:
            print(f'  d (scaled) = {np.array2string(d, precision=10, separator=", ")}')
            print(f'  e (scaled) = {np.array2string(e, precision=10, separator=", ")}')

    # Split bidiagonal and zero split entries
    blocks = split_bidiag(d, e, n)
    for bbeg, bend in blocks:
        if bend < n-1:
            e[bend] = 0.0

    if _verbose:
        print(f'  bidiag splits: {len(blocks)} blocks: {blocks}')

    # D1, D2 sign matrices (O(n) scalar loop — sequential dependency)
    d1 = np.ones(n); d2 = np.ones(n)
    d1[0] = 1.0; d2[0] = np.sign(d[0]) if d[0] != 0 else 1.0
    for i in range(n-1):
        s_e = np.sign(e[i]) if e[i] != 0 else 1.0
        d2[i+1] = s_e / d1[i] if abs(d1[i]) > 0 else s_e
        s_d = np.sign(d[i+1]) if d[i+1] != 0 else 1.0
        d1[i+1] = s_d / d2[i+1] if abs(d2[i+1]) > 0 else s_d

    if _verbose:
        print(f'  D1 = {d1}')
        print(f'  D2 = {d2}')

    sigma = np.zeros(n)
    U = np.zeros((n, n))
    V = np.zeros((n, n))

    # Classify blocks into singletons and multi-blocks
    # For multi-blocks: apply one zero-shift QR sweep to detect zero SVs
    # (Demmel-Kahan 1990: zero-shift QR = inverse iteration, converges in 1 step for σ=0)
    sing_rows = []; sing_cols = []; multi_blocks = []
    qr_deflated_cols = []
    zero_sv_thresh = n * EPS  # after prescaling, ||B||≈1
    col = 0
    for bbeg, bend in blocks:
        bk = bend - bbeg + 1
        if bk <= 0: continue
        if bk == 1:
            sing_rows.append(bbeg)
            sing_cols.append(col)
            col += 1
        else:
            # Sweep-then-check: one QR sweep to detect zero SV
            # (Demmel-Kahan 1990: zero-shift QR = inverse iteration,
            #  drives both d[-1] AND e[-1] to zero when σ_min = 0)
            d_sw = d[bbeg:bend+1].copy()
            e_sw = e[bbeg:bend].copy()
            right_rots, left_rots = zero_shift_qr_sweep(d_sw, e_sw)
            # Clean split requires BOTH d[-1] and e[-1] near zero
            clean_split = abs(d_sw[-1]) < zero_sv_thresh and abs(e_sw[-1]) < zero_sv_thresh
            if clean_split:
                # Zero SV detected — deflate
                if _verbose:
                    print(f'  [QR deflation] block [{bbeg},{bend}] bk={bk}: d[-1]={d_sw[-1]:.3e} < thresh={zero_sv_thresh:.3e}')
                # Solve (bk-1) sub-problem (unscale for recursive call)
                sigma_sub, U_sub, V_sub, _ = bidiag_svd(
                    d_sw[:-1] / scale_factor, e_sw[:-1] / scale_factor)
                # Embed into bk×bk: top-left = sub-solution, bottom-right = 1
                U_blk = np.eye(bk); U_blk[:-1, :-1] = U_sub
                V_blk = np.eye(bk); V_blk[:-1, :-1] = V_sub
                # Apply QR rotations: B_orig = U_qr @ B_swept @ V_qr^T
                _apply_givens_to_rows(V_blk, right_rots)
                _apply_givens_to_rows(U_blk, left_rots)
                # Store results (sigma in prescaled form for final /= scale_factor)
                c_slice = slice(col, col + bk)
                sigma[c_slice] = np.concatenate([sigma_sub * scale_factor, [abs(d_sw[-1])]])
                U[bbeg:bend+1, c_slice] = U_blk
                V[bbeg:bend+1, c_slice] = V_blk
                qr_deflated_cols.extend(range(col, col + bk))
                col += bk
            else:
                # No zero SV — use MR³ on original (d, e unchanged)
                multi_blocks.append((bbeg, bend, bk, col))
                col += bk

    if _verbose:
        print(f'  singletons: {len(sing_rows)}')
        print(f'  multi-blocks: {len(multi_blocks)} -> {[(b,e,k) for b,e,k,_ in multi_blocks][:10]}{"..." if len(multi_blocks)>10 else ""}')

    # --- Singletons: O(n) total, no intermediate matrix ---
    if sing_rows:
        sr = np.array(sing_rows)
        sc = np.array(sing_cols)
        sigma[sc] = np.abs(d[sr])
        V[sr, sc] = d2[sr]   # ±1, already unit vectors
        U[sr, sc] = d1[sr]   # ±1, already unit vectors

    # --- Multi-element blocks: XMR eigensolver ---
    multi_col_list = []
    for bbeg, bend, bk, col_offset in multi_blocks:
        w, Z = mr3_tgk_multiblock(d, e, n, bbeg, bend)
        nc = min(bk, n - col_offset)
        c_slice = slice(col_offset, col_offset + nc)
        r_slice = slice(bbeg, bend + 1)

        # Extract V, U from even/odd rows of T_GK eigenvectors, apply D1/D2
        V[r_slice, c_slice] = Z[0::2, :nc] * d2[bbeg:bend+1, None]
        U[r_slice, c_slice] = Z[1::2, :nc] * d1[bbeg:bend+1, None]
        sigma[c_slice] = np.abs(w[:nc])
        multi_col_list.extend(range(col_offset, col_offset + nc))

        if _verbose:
            print(f'  [extraction] block [{bbeg},{bend}] -> cols [{col_offset},{col_offset+nc}], sigma range [{np.min(sigma[c_slice]):.3e}, {np.max(sigma[c_slice]):.3e}]')

    # --- Post-processing: only on multi-block columns ---
    if multi_col_list:
        mc = np.array(multi_col_list)

        if _verbose:
            print(f'\n[post-processing] {len(mc)} multi-block columns')

        # Normalize (use system DNRM2 for bit-match with Fortran)
        nv_mc = _system_dnrm2_axis0(V[:, mc])
        nu_mc = _system_dnrm2_axis0(U[:, mc])

        if _verbose:
            print(f'  [normalize] V norms: {np.array2string(nv_mc, precision=6, separator=", ")}')
            print(f'  [normalize] U norms: {np.array2string(nu_mc, precision=6, separator=", ")}')

        max_norms = np.maximum(nv_mc, nu_mc)
        norm_thresh = max_norms * 1e-4
        v_good = (nv_mc > norm_thresh) & (nv_mc > SAFMIN)
        u_good = (nu_mc > norm_thresh) & (nu_mc > SAFMIN)
        scale_v = np.where(v_good, 1.0 / np.maximum(nv_mc, SAFMIN), 0.0)
        V[:, mc] *= scale_v[None, :]
        scale_u = np.where(u_good, 1.0 / np.maximum(nu_mc, SAFMIN), 0.0)
        U[:, mc] *= scale_u[None, :]

        if _verbose:
            bad_v = mc[~v_good]
            bad_u = mc[~u_good]
            if len(bad_v): print(f'  [normalize] V bad (zeroed): cols {bad_v.tolist()}')
            if len(bad_u): print(f'  [normalize] U bad (zeroed): cols {bad_u.tolist()}')

        # Bv recovery
        sigma_max = np.max(sigma)
        sig_thresh = n * EPS * sigma_max
        sig_mc = sigma[mc]
        need_recovery = sig_mc > sig_thresh
        recover_u = need_recovery & (~u_good) & v_good
        recover_v = need_recovery & (~v_good) & u_good
        both_good = need_recovery & v_good & u_good

        if np.any(recover_u):
            idx = mc[recover_u]
            if _verbose: print(f'  [Bv recovery] recovering U for cols {idx.tolist()}')
            BV = _bidiag_matvec_batch(d, e, V[:, idx])
            U_new = BV / sigma[idx][None, :]
            norms = np.maximum(_system_dnrm2_axis0(U_new), SAFMIN)
            U[:, idx] = U_new / norms[None, :]

        if np.any(recover_v):
            idx = mc[recover_v]
            if _verbose: print(f'  [Bv recovery] recovering V for cols {idx.tolist()}')
            BTU = _bidiagT_matvec_batch(d, e, U[:, idx])
            V_new = BTU / sigma[idx][None, :]
            norms = np.maximum(_system_dnrm2_axis0(V_new), SAFMIN)
            V[:, idx] = V_new / norms[None, :]

        if np.any(both_good):
            idx = mc[both_good]
            BV = _bidiag_matvec_batch(d, e, V[:, idx])
            BTU = _bidiagT_matvec_batch(d, e, U[:, idx])
            sig_idx = sigma[idx]
            res_v = _system_dnrm2_axis0(BV - sig_idx[None, :] * U[:, idx])
            res_u = _system_dnrm2_axis0(BTU - sig_idx[None, :] * V[:, idx])
            bnorm_val = max(np.max(np.abs(d)), np.max(np.abs(e[:n-1])) if n > 1 else 0.0)
            thresh_arr = 10 * n * EPS * np.maximum(sig_idx, bnorm_val)
            fix_v = (res_v > thresh_arr) & (res_v > res_u)
            fix_u = (res_u > thresh_arr) & (res_u > res_v) & (~fix_v)
            if _verbose:
                n_fix_v = np.sum(fix_v)
                n_fix_u = np.sum(fix_u)
                if n_fix_v or n_fix_u:
                    print(f'  [Bv check] {len(idx)} cols checked, {n_fix_v} fix_v, {n_fix_u} fix_u, max_res_v={np.max(res_v):.3e}, max_res_u={np.max(res_u):.3e}')
            if np.any(fix_v):
                fv_idx = idx[fix_v]
                if _verbose: print(f'  [Bv recovery] fixing V for cols {fv_idx.tolist()}')
                V_new = BTU[:, fix_v] / sigma[fv_idx][None, :]
                norms = np.maximum(_system_dnrm2_axis0(V_new), SAFMIN)
                V[:, fv_idx] = V_new / norms[None, :]
            if np.any(fix_u):
                fu_idx = idx[fix_u]
                if _verbose: print(f'  [Bv recovery] fixing U for cols {fu_idx.tolist()}')
                U_new = BV[:, fix_u] / sigma[fu_idx][None, :]
                norms = np.maximum(_system_dnrm2_axis0(U_new), SAFMIN)
                U[:, fu_idx] = U_new / norms[None, :]

        # Sign fix — only multi-block columns
        nonzero_sig = sig_mc >= SAFMIN
        if np.any(nonzero_sig):
            idx = mc[nonzero_sig]
            BV = _bidiag_matvec_batch(d, e, V[:, idx])
            dots = np.sum(BV * U[:, idx], axis=0)
            sign_flip = np.where(dots < 0, -1.0, 1.0)
            if _verbose:
                flipped = idx[sign_flip < 0]
                if len(flipped): print(f'  [sign fix] flipped U sign for cols {flipped.tolist()}')
            U[:, idx] *= sign_flip[None, :]

        # GS completion — at most 1 zero singular value per bidiag split,
        # so at most 1 zero-V col and 1 zero-U col. The missing vector
        # is the unique orthogonal complement of the other n-1 columns.
        all_v_norms = _system_dnrm2_axis0(V)
        all_u_norms = _system_dnrm2_axis0(U)
        zero_v_mc = mc[all_v_norms[mc] < 0.5]
        zero_u_mc = mc[all_u_norms[mc] < 0.5]
        if len(zero_v_mc) > 0 or len(zero_u_mc) > 0:
            if _verbose:
                print(f'  [GS completion] {len(zero_v_mc)} zero-V cols {zero_v_mc.tolist()}, {len(zero_u_mc)} zero-U cols {zero_u_mc.tolist()}')
        for j in zero_v_mc:
            good_mask = np.arange(n) != j
            good_mask &= all_v_norms > 0.5
            V_good = V[:, good_mask]
            # Two-pass MGS with random start. Deterministic seed = (1, 3, 5, 7+2j+1)
            # so we get reproducibility (and bit-match between Python+Fortran).
            v_cand = _dlarnv_normal(n, (1, 3, 5, (2*int(j) + 1) & 4095 | 1))
            v_cand -= V_good @ (V_good.T @ v_cand)  # pass 1
            v_cand -= V_good @ (V_good.T @ v_cand)  # pass 2
            nrm = _system_dnrm2(v_cand)
            if nrm > 1e-14:
                V[:,j] = v_cand / nrm
                if _verbose: print(f'    filled V[:,{j}] via orthogonal complement')
        for j in zero_u_mc:
            good_mask = np.arange(n) != j
            good_mask &= all_u_norms > 0.5
            U_good = U[:, good_mask]
            # Two-pass MGS with random start. Different fixed seed than V loop.
            u_cand = _dlarnv_normal(n, (2, 4, 6, (2*int(j) + 3) & 4095 | 1))
            u_cand -= U_good @ (U_good.T @ u_cand)  # pass 1
            u_cand -= U_good @ (U_good.T @ u_cand)  # pass 2
            nrm = _system_dnrm2(u_cand)
            if nrm > 1e-14:
                U[:,j] = u_cand / nrm
                if _verbose: print(f'    filled U[:,{j}] via orthogonal complement')

    # Scale sigma back if we pre-scaled
    if scale_factor != 1.0:
        sigma /= scale_factor

    if _verbose:
        print(f'\n[bidiag_svd] FINAL OUTPUT:')
        print(f'  sigma range=[{np.min(sigma):.10e}, {np.max(sigma):.10e}]')
        ort_U = np.max(np.abs(U.T @ U - np.eye(n)))
        ort_V = np.max(np.abs(V.T @ V - np.eye(n)))
        print(f'  max|U^TU-I|={ort_U:.3e} ({ort_U/(n*EPS):.1f} neps)')
        print(f'  max|V^TV-I|={ort_V:.3e} ({ort_V/(n*EPS):.1f} neps)')
        print(f'{"="*60}\n')

    return sigma, U, V, 0
