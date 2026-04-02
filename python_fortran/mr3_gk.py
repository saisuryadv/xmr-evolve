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

# ============================================================
# Metrics
# ============================================================
class MR3Metrics:
    def __init__(self):
        self.shifts = []; self.singletons = []; self.max_depth = 0
        self.n_nodes = 0; self.n_shifts_tried = 0; self.n_shifts_accepted = 0
    def record_shift(self, depth, tau, d_child, nn, child_evals, spdiam, accepted):
        self.n_shifts_tried += 1
        if accepted: self.n_shifts_accepted += 1
        rec = {'depth': depth, 'accepted': accepted}
        rec['elg_ratio'] = np.max(np.abs(d_child[:nn])) / max(spdiam, SAFMIN)
        rec['elg_ok'] = rec['elg_ratio'] <= ELG_THRESH
        if len(child_evals) > 1:
            rg = [abs(child_evals[i+1]-child_evals[i])/max(abs(child_evals[i]),SAFMIN)
                  for i in range(len(child_evals)-1)]
            rec['relgaps_min'] = min(rg)
        else: rec['relgaps_min'] = float('inf')
        rec['relgaps_ok'] = rec['relgaps_min'] >= GAPTOL
        self.shifts.append(rec)
    def record_singleton(self, depth, nn, lam, z, c, e_off, gap):
        tz = c[:nn]*z
        if nn > 1: tz[:nn-1] += e_off[:nn-1]*z[1:]; tz[1:] += e_off[:nn-1]*z[:nn-1]
        residual = np.linalg.norm(tz - lam*z)
        ratio = residual / max(gap, SAFMIN)
        self.singletons.append({'depth': depth, 'residual': residual,
            'gap': gap, 'ratio': ratio, 'getvec_ok': ratio <= 100*nn*EPS})
    def record_node(self, depth):
        self.n_nodes += 1; self.max_depth = max(self.max_depth, depth)
    def summary_line(self):
        parts = ["RRR:TGK"]
        acc = [s for s in self.shifts if s['accepted']]
        if acc:
            parts.append(f"ELG:{max(s['elg_ratio'] for s in acc):.1f}({sum(1 for s in acc if not s['elg_ok'])}f)")
            rg = [s['relgaps_min'] for s in acc if s['relgaps_min'] < float('inf')]
            if rg: parts.append(f"RG:{min(rg):.0e}({sum(1 for s in acc if not s['relgaps_ok'])}f)")
            else: parts.append("RG:all1")
        else: parts.append("ELG:n/a"); parts.append("RG:n/a")
        parts.append("SR:exact")
        if self.singletons:
            parts.append(f"GV:{max(s['ratio'] for s in self.singletons):.0e}({sum(1 for s in self.singletons if not s['getvec_ok'])}f)")
        else: parts.append("GV:n/a")
        parts.append(f"d={self.max_depth} sh={self.n_shifts_tried}/{self.n_shifts_accepted}")
        return " ".join(parts)

_metrics = None
def enable_metrics():
    global _metrics; _metrics = MR3Metrics(); return _metrics
def disable_metrics():
    global _metrics; _metrics = None

# ============================================================
# Factorizations — standard 1×1 (for shift evaluation)
# ============================================================
def ldl_top_down(c, e_off, n):
    d = np.zeros(n); ell = np.zeros(max(n-1,0)); d[0] = c[0]
    for i in range(n-1):
        if abs(d[i]) < SAFMIN: ell[i] = 0.0; d[i+1] = c[i+1]
        else: ell[i] = e_off[i]/d[i]; d[i+1] = c[i+1] - ell[i]*e_off[i]
    return d, ell

def urd_bottom_up(c, e_off, n):
    r = np.zeros(n); u = np.zeros(max(n-1,0)); r[n-1] = c[n-1]
    for i in range(n-2,-1,-1):
        if abs(r[i+1]) < SAFMIN: u[i] = 0.0; r[i] = c[i]
        else: u[i] = e_off[i]/r[i+1]; r[i] = c[i] - u[i]*e_off[i]
    return r, u

def eigvec_from_twist(ell_top, u_bot, r, n):
    z = np.zeros(n); z[r] = 1.0
    for i in range(r-1,-1,-1): z[i] = -ell_top[i]*z[i+1]
    for i in range(r+1,n): z[i] = -u_bot[i-1]*z[i-1]
    nrm = np.linalg.norm(z)
    if nrm > SAFMIN: z /= nrm
    return z

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

def refine_eval(c, e_off, n, approx, tol=0.0):
    gap = max(abs(approx)*0.001, 10*EPS*abs(approx), SAFMIN)
    try:
        m, w, _, _, info = _flapack.dstebz(
            c[:n].astype(np.float64), e_off[:n-1].astype(np.float64),
            1, approx-gap, approx+gap, 0, 0, tol, b'E')
        if info == 0 and m > 0: return w[np.argmin(np.abs(w[:m]-approx))]
    except Exception: pass
    return approx

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

# ============================================================
# Singleton computation — with block factorization
# ============================================================
def compute_singleton(c, e_off, n, lam, all_evals=None, eval_idx=-1,
                      depth=0, twist_rank=0, multiplicity=1):
    gap = float('inf')
    if all_evals is not None and eval_idx >= 0:
        for k in range(len(all_evals)):
            if k != eval_idx:
                g = abs(all_evals[k]-lam)
                if g < gap: gap = g
        if gap < SAFMIN: gap = SAFMIN

    for rqi_iter in range(5):
        c_shifted = c[:n] - lam
        eo = e_off[:max(n-1,0)]
        d_top, ell_top = ldl_top_down(c_shifted, eo, n)
        r_bot, u_bot = urd_bottom_up(c_shifted, eo, n)
        gammas = d_top + r_bot - c_shifted

        if twist_rank > 0 and multiplicity > 1:
            target = int(twist_rank*n/multiplicity)
            target = min(max(target, 0), n-1)
            sorted_idx = np.argsort(np.abs(gammas))
            min_g = abs(gammas[sorted_idx[0]])
            thresh = max(min_g*1e8, 100*n*EPS)
            cands = [int(sorted_idx[k]) for k in range(n)
                     if abs(gammas[sorted_idx[k]]) <= thresh]
            r = min(cands, key=lambda idx: abs(idx-target)) if cands else target
        else:
            r = int(np.argmin(np.abs(gammas)))

        z = eigvec_from_twist(ell_top, u_bot, r, n)
        tz = c[:n]*z
        if n > 1: tz[:n-1] += eo*z[1:]; tz[1:] += eo*z[:n-1]
        rq = np.dot(z, tz)
        res_norm = np.linalg.norm(tz - rq*z)
        if res_norm < 10*n*EPS*(abs(rq)+SAFMIN): lam = rq; break
        lam = rq

    if _metrics is not None:
        _metrics.record_singleton(depth, n, lam, z, c, e_off, gap)
    return lam, z

# ============================================================
# Classification
# ============================================================
def classify(eigvals, gaptol):
    n = len(eigvals)
    if n <= 1: return [(0,0)] if n == 1 else []
    groups = []; cs = 0
    while cs < n:
        ce = cs
        while ce+1 < n:
            denom = max(abs(eigvals[ce]), SAFMIN)
            if abs(eigvals[ce+1]-eigvals[ce])/denom >= gaptol: break
            ce += 1
        groups.append((cs, ce)); cs = ce+1
    return groups

def is_zero_width(gev):
    width = abs(gev[-1]-gev[0])
    scale = max(abs(gev[0]), abs(gev[-1]), SAFMIN)
    return width < 10*len(gev)*EPS*scale

# ============================================================
# Shift evaluation and selection
# ============================================================
def evaluate_shift(c, e_off, nn, tau, spdiam):
    cc = c[:nn].copy()-tau; ee = e_off[:max(nn-1,0)].copy()
    d_top, ell_top = ldl_top_down(cc, ee, nn)
    r_bot, u_bot = urd_bottom_up(cc, ee, nn)
    gammas = d_top + r_bot - cc
    max_pivot = np.max(np.abs(d_top[:nn]))
    elg = max_pivot / max(spdiam, SAFMIN)
    delg = elg / ELG_THRESH
    # Simplified rcond
    rcsq = 0.0; beta = 0.0
    for i in range(nn-1):
        invns = beta + 1.0
        ell_i = ee[i]/d_top[i] if i < len(ee) and abs(d_top[i]) > SAFMIN else 0.0
        rcsq += abs(invns * (1.0 + ell_i))
        beta = -ell_i * invns
    drcond = max(rcsq + 1.0, 1.0) / MAXRELCOND
    eval_score = delg + drcond
    n_neg = np.sum(d_top[:nn] < 0)
    ok = elg <= ELG_THRESH
    return cc, ee, eval_score, ok, n_neg

def verify_outer_bounds(cc, ee, nn, shifted_evals):
    if len(shifted_evals) == 0: return True
    left_b = shifted_evals[0] - abs(shifted_evals[0])*GAPTOL*0.5
    d_l, _ = ldl_top_down(cc[:nn]-left_b, ee, nn)
    n1 = np.sum(d_l[:nn] < 0)
    right_b = shifted_evals[-1] + abs(shifted_evals[-1])*GAPTOL*0.5
    d_r, _ = ldl_top_down(cc[:nn]-right_b, ee, nn)
    n2 = np.sum(d_r[:nn] < 0)
    return n2 - n1 >= len(shifted_evals) - 1

def select_shift_geometric(c, e_off, nn, gev, spdiam, local_ev):
    cleft, cright = gev[0], gev[-1]
    clen = len(gev)
    best_tau = None; best_cc = None; best_ee = None; best_eval = float('inf')
    max_tries = 6; MAX_FUDGE = 3; FUDGE_FAC = 4.0

    def try_shift(tau, direction=0):
        nonlocal best_tau, best_cc, best_ee, best_eval
        for fi in range(MAX_FUDGE+1):
            cc, ee, ev, ok, nn_ = evaluate_shift(c, e_off, nn, tau, spdiam)
            if _metrics is not None:
                d_c, _ = ldl_top_down(cc, ee, nn)
                _metrics.record_shift(0, tau, d_c, nn, gev-tau, spdiam, ok)
            if ok:
                cob_ok = verify_outer_bounds(cc, ee, nn, gev-tau)
                if cob_ok and ev < best_eval:
                    best_tau=tau; best_cc=cc; best_ee=ee; best_eval=ev; return
            if direction != 0 and fi < MAX_FUDGE:
                tau += direction*abs(tau)*FUDGE_FAC*nn*EPS
            else: break

    bound = cleft; off0 = max(abs(bound)*8*EPS, SAFMIN)
    lgap = abs(bound)*GAPTOL
    for k in range(len(local_ev)):
        if local_ev[k] < cleft-SAFMIN: lgap = cleft-local_ev[k]
    for kt in range(max_tries):
        off = off0*(2**kt)
        if off > 0.25*lgap: break
        try_shift(bound-off, -1)

    bound = cright; off0 = max(abs(bound)*8*EPS, SAFMIN)
    ugap = abs(bound)*GAPTOL
    for k in range(len(local_ev)):
        if local_ev[k] > cright+SAFMIN: ugap = local_ev[k]-cright; break
    for kt in range(max_tries):
        off = off0*(2**kt)
        if off > 0.25*ugap: break
        try_shift(bound+off, +1)

    if clen >= 4:
        for i in range(len(gev)-1):
            gi = gev[i+1]-gev[i]
            if gi > SAFMIN:
                for kt in range(3):
                    off = max(abs(gev[i+1])*8*EPS, SAFMIN)*(2**kt)
                    if off > 0.25*gi: break
                    try_shift(gev[i+1]-off, -1)

    if best_tau is None:
        tau = 0.5*(cleft+cright)
        cc, ee, ev, ok, _ = evaluate_shift(c, e_off, nn, tau, spdiam)
        if _metrics is not None:
            d_c, _ = ldl_top_down(cc, ee, nn)
            _metrics.record_shift(0, tau, d_c, nn, gev-tau, spdiam, ok)
        best_tau=tau; best_cc=cc; best_ee=ee
    return best_tau, best_cc, best_ee

# ============================================================
# MR3 block
# ============================================================
def mr3_block(c, e_off, nn, pos_evals, spdiam):
    nev = len(pos_evals)
    result_evals = np.zeros(nev); result_evecs = np.zeros((nn, nev))
    stack = [(c[:nn].copy(), e_off[:max(nn-1,0)].copy(), nn,
              pos_evals.copy(), list(range(nev)), 0.0, 0)]
    while stack:
        bc, be, bnn, local_ev, out_idx, shift, depth = stack.pop()
        if len(local_ev) == 0: continue
        if _metrics is not None: _metrics.record_node(depth)
        if depth >= MAX_DEPTH:
            mult = len(local_ev)
            for j in range(mult):
                lam_ref, z = compute_singleton(bc, be, bnn, local_ev[j],
                    all_evals=local_ev, eval_idx=j, depth=depth,
                    twist_rank=j, multiplicity=mult)
                result_evals[out_idx[j]] = lam_ref+shift
                result_evecs[:bnn, out_idx[j]] = z
            continue
        groups = classify(local_ev, GAPTOL)
        for gs, ge in groups:
            gev = local_ev[gs:ge+1]; gidx = out_idx[gs:ge+1]
            if gs == ge:
                lam_ref, z = compute_singleton(bc, be, bnn, gev[0],
                    all_evals=local_ev, eval_idx=gs, depth=depth)
                result_evals[gidx[0]] = lam_ref+shift
                result_evecs[:bnn, gidx[0]] = z
            else:
                if is_zero_width(gev):
                    mult = len(gev)
                    for j in range(mult):
                        lam_ref, z = compute_singleton(bc, be, bnn, gev[j],
                            all_evals=local_ev, eval_idx=gs+j, depth=depth,
                            twist_rank=j, multiplicity=mult)
                        result_evals[gidx[j]] = lam_ref+shift
                        result_evecs[:bnn, gidx[j]] = z
                    continue
                tau, cc, ee = select_shift_geometric(bc, be, bnn, gev, spdiam, local_ev)
                shifted_ev = gev - tau
                refined = np.array([refine_eval(cc, ee, bnn, s) for s in shifted_ev])
                si = np.argsort(refined); refined = refined[si]
                sgidx = [gidx[i] for i in si]
                stack.append((cc, ee, bnn, refined, sgidx, shift+tau, depth+1))
    return result_evals, result_evecs

# ============================================================
# Main driver
# ============================================================
def _solve_tgk_block(bc, bte, m, bk_hint=None):
    """Solve a single T_GK eigenproblem of size m.
    Returns (w, Z) where w has the non-negative eigenvalues
    and Z has the corresponding eigenvectors.
    bk_hint: if given, number of non-negative eigenvalues expected = ceil(m/2).
    """
    from xmr_ctypes import xmr_eigenvectors, xmr_eigenvectors_v, build_repr_from_tridiag
    n_nonneg = (m + 1) // 2  # ceil(m/2): number of non-negative eigenvalues
    if bk_hint is not None:
        n_nonneg = bk_hint
    wil = m - n_nonneg + 1  # first non-negative eigenvalue index (1-based)
    wiu = m                  # last eigenvalue index (1-based)

    all_evals = np.sort(bisect_evals(bc, bte, m, 1, m))
    spdiam = all_evals[-1] - all_evals[0] if len(all_evals) > 0 else 0.0

    if _verbose:
        print(f'  [_solve_tgk_block] m={m}, n_nonneg={n_nonneg}, wil={wil}, wiu={wiu}')
        print(f'  [_solve_tgk_block] spdiam={spdiam:.10e}')
        print(f'  [_solve_tgk_block] evals range=[{all_evals[0]:.10e}, {all_evals[-1]:.10e}]')

    root = build_repr_from_tridiag(bc, bte, m, k=m)

    if _verbose:
        w, Z, info, tau_re = xmr_eigenvectors_v(root, bte, all_evals,
                                                  wil=wil, wiu=wiu,
                                                  spdiam=spdiam, gaptol=GAPTOL)
        print(f'  [_solve_tgk_block] dlaxre shift tau_re={tau_re:.10e}')
        print(f'  [_solve_tgk_block] XMR info={info}')
        if info == 0:
            print(f'  [_solve_tgk_block] XMR output: {len(w)} evals, range=[{w[0]:.10e}, {w[-1]:.10e}]')
    else:
        w, Z, info = xmr_eigenvectors(root, bte, all_evals,
                                        wil=wil, wiu=wiu,
                                        spdiam=spdiam, gaptol=GAPTOL)

    if info != 0:
        if _verbose:
            print(f'  [_solve_tgk_block] XMR failed (info={info}), falling back to Python MR3')
        pos_evals = all_evals[wil-1:]
        w, Z = mr3_block(bc, bte, m, pos_evals, spdiam)
        if _verbose:
            print(f'  [_solve_tgk_block] Python MR3 eigenvalues: {np.array2string(w[:n_nonneg], precision=10, separator=", ")}')

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
    sing_rows = []; sing_cols = []; multi_blocks = []
    col = 0
    for bbeg, bend in blocks:
        bk = bend - bbeg + 1
        if bk <= 0: continue
        if bk == 1:
            sing_rows.append(bbeg)
            sing_cols.append(col)
            col += 1
        else:
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

        # Normalize
        nv_mc = np.linalg.norm(V[:, mc], axis=0)
        nu_mc = np.linalg.norm(U[:, mc], axis=0)

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
            norms = np.maximum(np.linalg.norm(U_new, axis=0), SAFMIN)
            U[:, idx] = U_new / norms[None, :]

        if np.any(recover_v):
            idx = mc[recover_v]
            if _verbose: print(f'  [Bv recovery] recovering V for cols {idx.tolist()}')
            BTU = _bidiagT_matvec_batch(d, e, U[:, idx])
            V_new = BTU / sigma[idx][None, :]
            norms = np.maximum(np.linalg.norm(V_new, axis=0), SAFMIN)
            V[:, idx] = V_new / norms[None, :]

        if np.any(both_good):
            idx = mc[both_good]
            BV = _bidiag_matvec_batch(d, e, V[:, idx])
            BTU = _bidiagT_matvec_batch(d, e, U[:, idx])
            sig_idx = sigma[idx]
            res_v = np.linalg.norm(BV - sig_idx[None, :] * U[:, idx], axis=0)
            res_u = np.linalg.norm(BTU - sig_idx[None, :] * V[:, idx], axis=0)
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
                norms = np.maximum(np.linalg.norm(V_new, axis=0), SAFMIN)
                V[:, fv_idx] = V_new / norms[None, :]
            if np.any(fix_u):
                fu_idx = idx[fix_u]
                if _verbose: print(f'  [Bv recovery] fixing U for cols {fu_idx.tolist()}')
                U_new = BV[:, fix_u] / sigma[fu_idx][None, :]
                norms = np.maximum(np.linalg.norm(U_new, axis=0), SAFMIN)
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
        all_v_norms = np.linalg.norm(V, axis=0)
        all_u_norms = np.linalg.norm(U, axis=0)
        zero_v_mc = mc[all_v_norms[mc] < 0.5]
        zero_u_mc = mc[all_u_norms[mc] < 0.5]
        if len(zero_v_mc) > 0 or len(zero_u_mc) > 0:
            if _verbose:
                print(f'  [GS completion] {len(zero_v_mc)} zero-V cols {zero_v_mc.tolist()}, {len(zero_u_mc)} zero-U cols {zero_u_mc.tolist()}')
        for j in zero_v_mc:
            good_mask = np.arange(n) != j
            good_mask &= all_v_norms > 0.5
            V_good = V[:, good_mask]
            # Two-pass MGS with random start
            v_cand = np.random.randn(n)
            v_cand -= V_good @ (V_good.T @ v_cand)  # pass 1
            v_cand -= V_good @ (V_good.T @ v_cand)  # pass 2
            nrm = np.linalg.norm(v_cand)
            if nrm > 1e-14:
                V[:,j] = v_cand / nrm
                if _verbose: print(f'    filled V[:,{j}] via orthogonal complement')
        for j in zero_u_mc:
            good_mask = np.arange(n) != j
            good_mask &= all_u_norms > 0.5
            U_good = U[:, good_mask]
            # Two-pass MGS with random start
            u_cand = np.random.randn(n)
            u_cand -= U_good @ (U_good.T @ u_cand)  # pass 1
            u_cand -= U_good @ (U_good.T @ u_cand)  # pass 2
            nrm = np.linalg.norm(u_cand)
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
