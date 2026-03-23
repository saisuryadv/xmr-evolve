# MR³-GK Changes Summary: Full History (0 → 300/379)

## Phase 1: Architecture (Sessions 1-4, ~22/64)
Built the entire MR³-GK pipeline from scratch:
- Studied Willems-Lang 2012 paper and XMR Fortran source (~15K lines)
- Implemented Python `mr3_gk.py`: T_GK construction, LAPACK bisection for eigenvalues, bidiag splitting, sign recovery (D1/D2 diagonal matrices), Bv recovery, orthogonal completion
- Compiled all XMR Fortran objects (~30 .o files) into `libxmr.so`
- Built C wrapper (`xmr_wrapper.c`) calling `dlaxre_` + `dlaxrv_` via ctypes
- Python calls: `dstebz` (bisection) → `xmr_eigenvectors()` (C wrapper) → Fortran XMR pipeline

## Phase 2: Negative Off-Diagonal Fix (Session 5, 22→23/64)

**The bug**: XMR's `dlaxre` requires positive off-diagonal entries. When the bidiagonal has negative d[i], T_GK gets negative off-diagonals at positions 2i. XMR builds a wrong representation → eigenvectors are garbage (residuals 0.86 instead of 1e-15).

**Fix**: Diagonal similarity transform before calling XMR.
```python
# Make all off-diags positive: D*T*D where D = diag(±1)
signs = np.ones(m2n)
for i in range(m2n-1):
    if tgk_e[i] < 0: signs[i+1] = -signs[i]
    else: signs[i+1] = signs[i]
tgk_e_pos = np.abs(tgk_e)
# Call XMR on positive-offdiag matrix, then Z[:,j] *= signs
```

## Phase 3: GK Detection in dlaxre (Session 6, 23→142/379)

**The bug**: XMR's `dlaxre.f` had GK structure detection **disabled**: `IF(.FALSE.) THEN` guarded the zero-diagonal check. Without GK detection, dlaxre used standard symmetric tridiagonal root representations instead of GK-specific ones, losing the GK structure that ensures ||u||=||v||=1/√2.

**Fix**: Modified `dlaxre_gk.f` to detect zero diagonals (T_GK always has c[i]=0) and enable the GK-specific root representation path. This includes the NEGL guard that checks whether the GK root actually satisfies quality criteria, falling back to PD (positive definite) shift when it doesn't.

**Impact**: 23/64 → 142/379 (massive jump, as GK structure preservation is the core of the algorithm)

## Phase 4: Sign Recovery in bidiag_svd (Session 6, included in 142)

**The fix**: T_GK uses |d[i]| (absolute values), so its eigenvectors are for |B|, not B. Need diagonal sign matrices D1, D2 where B = D1·|B|·D2 to recover correct U, V:
```python
d1[0] = 1.0; d2[0] = sign(d[0])
for i in range(n-1):
    d2[i+1] = sign(e[i]) / d1[i]
    d1[i+1] = sign(d[i+1]) / d2[i+1]
V[:,j] = ev[0::2] * d2;  U[:,j] = ev[1::2] * d1
```

## Phase 5: Bv Recovery and Orthogonal Completion (Sessions 6-7, up to 290)

**Bv recovery**: For each singular triplet, check residuals B·v - σ·u and Bᵀ·u - σ·v. If one side has larger residual, recompute it from the other. This compensates for cases where XMR's eigenvectors have uneven quality on the two GK halves.

**Orthogonal completion via Gram-Schmidt**: For zero singular values, one side of the eigenvector is zero. GS fills in the missing side to maintain U'U = I. (Later optimized with singleton handling to avoid O(n³) pathology.)

**Result**: 142/379 → 290/379 (76.5%)

---

## Phase 6: 290 → 300 (Sessions 8-9)

### Starting Point
- **290/379 tests passing** (76.5%)
- Key failures: chkbd (oU=4.5×10¹³), gl_gradm/gl_gradp@200 (oU=2.25×10¹³), many_near_zero, gl_gro2, demmel_S1pe, demmel_G1, stemr_killer@10 (oU=8.4×10⁶)

### Files Changed
- `mr3_gk.py` — 3 changes (Python only)
- `dlaxre_gk.f` — unchanged (NEGL guard kept)
- `libxmr.so` — unchanged

### Change 4: Two-Phase Bidiagonal Splitting

**Location**: `split_bidiag()` function

**Before**: Single-pass relative splitting only.

**After**: Two-phase splitting per the Willems-Lang 2012 paper.
- **Phase 1 (relative)**: Same as before — preserves relative accuracy.
- **Phase 2 (absolute)**: Within each relative-split block, compute `||B_sub||_∞` and split any `e[i]` below `n_sub × ε × ||B_sub||`. This enforces condition (3.5).

**Why per-block norm matters**: Global `||B||` over-splits matrices with large dynamic range (e.g., B_bug316_gesdd has d[1]=-2.3×10²⁵).

### Change 5: Consistency Between Splitting and Bv Recovery

**Location**: `bidiag_svd()`, before calling `mr3_tgk()`

**The bug**: Split points determined block boundaries, but the Bv recovery used the **original unsplit** `e` array. Cross-block contamination through unsplit e[i] destroyed orthogonality.

**Fix**: Zero out `e[i]` at every block boundary before passing to `mr3_tgk` and Bv recovery.

### Change 6: Direct Singleton Block Handling

**Location**: `mr3_tgk()`, at the start of the block loop

**The bug**: For singleton blocks after splitting, XMR was called on trivial 2×2 T_GK matrices. XMR returned eigenvectors with ||v||=1, ||u||=0 (only populating even positions). Then `bidiag_svd`'s orthogonal completion saw hundreds of columns with ||u||<0.5 and ran Gram-Schmidt on each: O(n³) → timeout.

**Fix**: Handle singletons directly — σ=|d[i]|, both eigenvector sides set to 1/√2.

**Impact**: chkbd@400 from timeout (>60s) to 0.025s.

### Result
- **300/379 tests passing** (79.2%), up from 290/379 (76.5%)
- +10 tests gained, 0 regressions

## What We Investigated But Did NOT Change

### NEGL Guard Removal
Removing the guard caused regressions on stemr_killer where PD shift gives better results. The guard correctly routes different matrix structures to different algorithms.

### Gram-Schmidt Replacement with Null-Space Recurrence
O(n) null-space computation works for genuine zero-σ from GK blocks, but fails for PD-fallback matrices where σ≈0 from numerical noise. GS remains necessary for those cases, but only triggers once per unsplit block (O(n²) total).
