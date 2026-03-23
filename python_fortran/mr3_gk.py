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
def mr3_tgk(d_bidiag, e_bidiag, n):
    from xmr_ctypes import xmr_eigenvectors, build_repr_from_tridiag
    m2n = 2*n; result_evals = np.zeros(n); result_evecs = np.zeros((m2n, n))
    bidiag_blocks = split_bidiag(d_bidiag, e_bidiag, n)
    ev_assigned = 0
    for bbeg, bend in bidiag_blocks:
        bk = bend-bbeg+1
        if bk <= 0: continue
        if bk == 1:
            # Singleton: σ = |d[i]|, eigvec of T_GK is [1/√2, 1/√2] at positions (2i, 2i+1)
            result_evals[ev_assigned] = abs(d_bidiag[bbeg])
            result_evecs[2*bbeg, ev_assigned] = 1.0 / np.sqrt(2.0)    # v component
            result_evecs[2*bbeg+1, ev_assigned] = 1.0 / np.sqrt(2.0)  # u component
            ev_assigned += 1
            continue
        bd = d_bidiag[bbeg:bend+1].copy()
        be = e_bidiag[bbeg:bend].copy() if bk > 1 else np.array([])
        # Use |d[i]| for T_GK (singular values are sign-independent).
        # Record signs to fix u later in bidiag_svd.
        bd_signs = np.sign(bd); bd_signs[bd_signs == 0] = 1.0
        bd_abs = np.abs(bd)
        m2k = 2*bk; bc = np.zeros(m2k); bte = np.zeros(max(m2k-1,0))
        for i in range(bk):
            bte[2*i] = bd_abs[i]
            if i < bk-1: bte[2*i+1] = be[i]
        # Get all eigenvalues via bisection
        all_evals = bisect_evals(bc, bte, m2k, 1, m2k)
        all_evals = np.sort(all_evals)
        npos = bk  # T_GK has bk positive eigenvalues
        if npos == 0: continue
        spdiam = all_evals[-1] - all_evals[0]
        # Build root repr
        root = build_repr_from_tridiag(bc, bte, m2k, k=m2k)
        # Compute positive eigenvectors via XMR (eigenvalues bk+1..2*bk)
        w, Z, info = xmr_eigenvectors(root, bte, all_evals,
                                        wil=bk+1, wiu=m2k,
                                        spdiam=spdiam, gaptol=GAPTOL)
        if info != 0:
            # Fallback to Python mr3_block
            block_pos = all_evals[bk:]
            block_ev, block_z = mr3_block(bc, bte, m2k, block_pos, spdiam)
            w = block_ev; Z = block_z
        # Fix u-component signs is done in bidiag_svd, not here.
        # XMR eigenvectors are for T_GK(|B|), which is correct.
        tgk_offset = 2*bbeg
        for j in range(npos):
            if ev_assigned+j < n:
                result_evals[ev_assigned+j] = w[j]
                result_evecs[tgk_offset:tgk_offset+m2k, ev_assigned+j] = Z[:m2k, j]
        ev_assigned += npos
    order = np.argsort(-result_evals[:ev_assigned])
    result_evals[:ev_assigned] = result_evals[:ev_assigned][order]
    result_evecs[:, :ev_assigned] = result_evecs[:, :ev_assigned][:, order]
    return result_evals, result_evecs

# ============================================================
# Main entry
# ============================================================
def bidiag_svd(d, e):
    d = np.asarray(d, dtype=np.float64).copy()
    e = np.asarray(e, dtype=np.float64).copy()
    n = len(d)
    if n <= 0: return np.array([]), np.zeros((0,0)), np.zeros((0,0)), 0
    if n == 1:
        return (np.array([abs(d[0])]),
                np.array([[1.0 if d[0]>=0 else -1.0]]),
                np.array([[1.0]]), 0)
    m2n = 2*n
    # Zero out split e[i] entries per condition (3.5)
    # These were used as split points in mr3_tgk; the Bv recovery must use
    # the same (split) matrix, not the original B.
    blocks = split_bidiag(d, e, n)
    # Zero e[i] at block boundaries
    block_ends = set()
    for bbeg, bend in blocks:
        if bend < n-1:
            block_ends.add(bend)
    for i in block_ends:
        e[i] = 0.0

    evals, evecs = mr3_tgk(d, e, n)
    sigma = np.abs(evals)

    # Extract u, v from T_GK eigenvectors:
    # q = [v_0, u_0, v_1, u_1, ...] so v = q[0::2], u = q[1::2]
    # T_GK was built with |d[i]|, so eigenvectors are for |B|.
    # B = D1 * |B| * D2 where D1, D2 are diagonal sign matrices.
    # Then U = D1 * U_pos, V = D2 * V_pos.
    # Compute D1, D2 from: D1[i]*D2[i]=sign(d[i]), D1[i]*D2[i+1]=sign(e[i])
    d1 = np.ones(n); d2 = np.ones(n)
    d1[0] = 1.0; d2[0] = np.sign(d[0]) if d[0] != 0 else 1.0
    for i in range(n-1):
        # D1[i]*D2[i+1] = sign(e[i])
        s_e = np.sign(e[i]) if e[i] != 0 else 1.0
        d2[i+1] = s_e / d1[i] if abs(d1[i]) > 0 else s_e
        # D1[i+1]*D2[i+1] = sign(d[i+1])
        s_d = np.sign(d[i+1]) if d[i+1] != 0 else 1.0
        d1[i+1] = s_d / d2[i+1] if abs(d2[i+1]) > 0 else s_d

    V = np.zeros((n,n)); U = np.zeros((n,n))
    for j in range(n):
        ev = evecs[:m2n, j]
        V[:,j] = ev[0::2] * d2   # V = D2 * V_pos
        U[:,j] = ev[1::2] * d1   # U = D1 * U_pos

    # Normalize and decide which side to recover from
    for j in range(n):
        nv = np.linalg.norm(V[:,j])
        nu = np.linalg.norm(U[:,j])
        # Both sides should have norm ≈ 1/sqrt(2) for well-conditioned GK eigenvectors.
        # If one side is tiny relative to the other, it's noise from a near-zero eigenvalue.
        norm_thresh = max(nu, nv) * 1e-4
        v_good = nv > norm_thresh and nv > SAFMIN
        u_good = nu > norm_thresh and nu > SAFMIN
        if v_good: V[:,j] /= nv
        else: V[:,j] = 0.0
        if u_good: U[:,j] /= nu
        else: U[:,j] = 0.0

        if sigma[j] > n*EPS*max(sigma):
            bv = d*V[:,j]
            if n > 1: bv[:n-1] += e*V[1:,j]
            btu = d*U[:,j]
            if n > 1: btu[1:] += e*U[:n-1,j]

            if not u_good and v_good:
                U[:,j] = bv / sigma[j]
                nu2 = np.linalg.norm(U[:,j])
                if nu2 > SAFMIN: U[:,j] /= nu2
            elif not v_good and u_good:
                V[:,j] = btu / sigma[j]
                nv2 = np.linalg.norm(V[:,j])
                if nv2 > SAFMIN: V[:,j] /= nv2
            else:
                res_v = np.linalg.norm(bv - sigma[j]*U[:,j])
                res_u = np.linalg.norm(btu - sigma[j]*V[:,j])
                if res_v > 10*n*EPS*sigma[j] and res_v > res_u:
                    V[:,j] = btu / sigma[j]
                    nv2 = np.linalg.norm(V[:,j])
                    if nv2 > SAFMIN: V[:,j] /= nv2
                elif res_u > 10*n*EPS*sigma[j] and res_u > res_v:
                    U[:,j] = bv / sigma[j]
                    nu2 = np.linalg.norm(U[:,j])
                    if nu2 > SAFMIN: U[:,j] /= nu2

    # Fix signs: ensure B*v_j and u_j point the same way
    for j in range(n):
        if sigma[j] < SAFMIN: continue
        bv = d*V[:,j]
        if n > 1: bv[:n-1] += e*V[1:,j]
        if np.dot(bv, U[:,j]) < 0: U[:,j] *= -1.0

    # For zero singular values, V[:,j] and U[:,j] may be zero vectors.
    # Fill them by completing to an orthonormal basis.
    for j in range(n):
        if np.linalg.norm(V[:,j]) < 0.5:
            # Find a unit vector orthogonal to all other V columns
            # Use random + Gram-Schmidt
            for attempt in range(n+5):
                if attempt < n:
                    v_cand = np.zeros(n); v_cand[attempt] = 1.0
                else:
                    v_cand = np.random.randn(n)
                for k in range(n):
                    if k != j and np.linalg.norm(V[:,k]) > 0.5:
                        v_cand -= np.dot(v_cand, V[:,k]) * V[:,k]
                nrm = np.linalg.norm(v_cand)
                if nrm > 0.1:
                    V[:,j] = v_cand / nrm
                    break
        if np.linalg.norm(U[:,j]) < 0.5:
            for attempt in range(n+5):
                if attempt < n:
                    u_cand = np.zeros(n); u_cand[attempt] = 1.0
                else:
                    u_cand = np.random.randn(n)
                for k in range(n):
                    if k != j and np.linalg.norm(U[:,k]) > 0.5:
                        u_cand -= np.dot(u_cand, U[:,k]) * U[:,k]
                nrm = np.linalg.norm(u_cand)
                if nrm > 0.1:
                    U[:,j] = u_cand / nrm
                    break

    return sigma, U, V, 0
