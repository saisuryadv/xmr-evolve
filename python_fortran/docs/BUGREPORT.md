# Bug Report — Modifications to Upstream XMR for MR3-GK Bidiagonal SVD

**Repo:** `xmr-evolve` (this directory: `python_fortran/`)
**Upstream:** Willems' XMR Fortran source (`xmr_src/`, ~14 K LOC, untouched)
**As of:** 2026-04-29

This is the consolidated bug report for advisor review. It lists every modification we made on top of upstream XMR, what each one fixes, and **empirical evidence** from running our 379-spec evaluation suite under each combination of patches.

For full diffs and the ablation details: `docs/02_xmr_modifications.md`. For test counts per variant: `experiments/ablation_results.json`. For the verified 21,375-spec Python ↔ Fortran bit-match: `docs/01_fortran_test_report.md`.

---

## 1. Summary of changes

| # | Change | Type | File | LOC | Required for |
|---|---|---|---|---|---|
| 1 | Activate GK-path in root representation builder | XMR Fortran fix | `dlaxre_gk.f` (replaces `xmr_src/dlaxre.f` at link time) | 7 | Every bidiagonal-SVD spec |
| 2 | `AVGTHRESH > 0` guard in cluster classifier | XMR Fortran fix | `dlaxrb_clssfy_fix.f` (replaces `xmr_src/dlaxrb_clssfy.f`) | ~25 | Tight-cluster specs (49 in our suite) |
| 3 | Reproducible random vectors via LAPACK DLARNV | Python helper | `mr3_gk.py:_dlarnv_normal` (lines 18-32) | 15 | Py ↔ Fo bit-identity on GS-completion path |
| 4 | Unified system BLAS DNRM2 via ctypes | Python helper | `mr3_gk.py:_system_dnrm2*` (lines 34-65) | 30 | Py ↔ Fo bit-identity on every normalization |
| 5 | Pure-Fortran orchestration layer | New Fortran | `mr3gk_fortran/*.f90` (8 files) | ~1,400 | Replaces Python wrapper; no algorithmic change |

**Total upstream XMR Fortran changed: 32 lines across 2 files.** Everything else is new code that orchestrates the unmodified XMR kernel.

---

## 2. Bug #1 — `dlaxre_gk.f`: GK-path was disabled

### Defect (in upstream)

`xmr_src/dlaxre.f:269` contains:
```fortran
      IF( .FALSE. )THEN
C        Support for GK-type matrices deactivated.
C        Reason: One should call a specific SVD routine in this case,
*         see below. We ran into problems for GRO-type synthetics,
*         mainly because full blocks in the root restrict possible
*         choices for blocks in the children too much, due to
*         overlap sequences.
```

The branch that handles **constant-diagonal tridiagonals** (T_GK matrices, i.e. the perfect-shuffle of a bidiagonal `B` whose diagonal is identically 0) is dead code. T_GK matrices fall through to the post-GK Gershgorin path, which produces an unusable root representation.

### Fix

Replace `IF( .FALSE. )` with the standard zero-diagonal test:
```fortran
      IF( ABS(DMAX)+ABS(DMIN) .LT. EPS*EMAX )THEN
C        GK-type matrix: constant (zero) diagonal.
C        Use T_GK directly as root per Algorithm 4.1 (Willems-Lang 2012).
C        No outside shift needed — preserves GK structure.
```

### Reference

Willems & Lang, ETNA 2012, Algorithm 4.1; also LAWN 163 §5.4.

### Empirical evidence

Running the 379-spec evaluation with this patch reverted (variant **revert-dlaxre**):

- **Pass: 233 / 379  (146 failures)**
- Failure mode: orthogonality blow-up to ~10¹¹ × n·eps. The singular values are still computed correctly via DSTEBZ; only the eigenvectors are wrong.
- 5 worst regressions (res / ortU / ortV — patched-baseline values in parentheses):
  - `exponential_graded@10` — 0.41 / **1.32 × 10¹²** / **1.02 × 10¹²**  (0.22 / 0.15)
  - `gl_gradp@100` — 0.64 / **5.89 × 10¹¹** / **6.19 × 10¹¹**  (0.09 / 0.03)
  - `huge_condition@100` — 0.90 / **4.57 × 10¹¹** / **4.31 × 10¹¹**  (0.36 / 0.28)
  - `exponential_graded@100` — 0.79 / **4.05 × 10¹¹** / **3.92 × 10¹¹**  (0.17 / 0.30)
  - `gl_gradp@200` — 0.64 / **2.22 × 10¹¹** / **2.20 × 10¹¹**  (0.05 / 0.01)

**Conclusion:** required for any MR3-GK bidiagonal SVD. Without it, only ~60 % of the suite passes.

---

## 3. Bug #2 — `dlaxrb_clssfy_fix.f`: cluster classifier's `MIN(0, x)` collapse at deep tree nodes

### Defect (in upstream)

`xmr_src/dlaxrb_clssfy.f:340-361` (and two near-identical sites at 414, 445):
```fortran
      AVGTHRESH = TWO * AVGTOL
      GAPTHRESH = TWO * GAPTOL
      IF( .NOT. DOAVG )  AVGTHRESH = SPDIAM
      ...
      IF( SPREAD .GE. MIN( AVGTHRESH, ABSMAX*GAPTHRESH ) )THEN
         RGINFO(J) = GI_FULLGAP
      ELSE
         RGINFO(J) = GI_NOFULL
      ENDIF
```

`AVGTHRESH = 2·AVGTOL`, where AVGTOL is scaled by `2^(-DEPTH)` and **underflows to 0** at deep tree nodes. Then `MIN(0, ABSMAX·GAPTHRESH) = 0`, the test becomes `SPREAD ≥ 0` — **trivially true** — and every adjacent eigenvalue pair is classified as a full gap. Clusters never form, and the MR3 representation tree degenerates.

### Fix

Explicit case-split when AVGTHRESH has underflowed:
```fortran
      LOGICAL  ::  GAPOK
      ...
      IF( AVGTHRESH .GT. ZERO )THEN
         GAPOK = SPREAD .GE. MIN( AVGTHRESH, ABSMAX*GAPTHRESH )
      ELSE
         GAPOK = SPREAD .GE. ABSMAX*GAPTHRESH
      ENDIF
      IF( GAPOK )THEN
         RGINFO(J) = GI_FULLGAP
      ELSE
         RGINFO(J) = GI_NOFULL
      ENDIF
```
Applied at lines 358-366, 414-417, 445-453.

### Empirical evidence

Variant **revert-clssfy** (only this patch reverted):

- **Pass: 330 / 379  (49 failures)**
- Failure mode: orthogonality 10-35 × n·eps on tight-cluster specs.
- 5 worst regressions:
  - `ST_B_20_graded` — res **592.4** / ortU **2782.5** / ortV **2782.5**  (patched: 0.24 / 0.26)
  - `demmel_S1pe_k4@100` — ortU **33.23** / ortV **33.22**  (patched: 3.62 / 3.64)
  - `pd_T0@100` — res **5.62** / ortU **32.97** / ortV **32.97**  (patched: 0.34 / 0.72 / 0.72)
  - `gl_abcon3@100` — ortU **24.28** / ortV **24.28**  (patched: 1.11 / 1.11)
  - `three_clusters@100` — ortU **21.62** / ortV **21.63**  (patched: 0.23 / 0.23)

**Conclusion:** required for tight-cluster matrices (Demmel S1/S2, three_clusters, pd_T0, gl_abcon, gl_wilkw, etc.). The README's "fundamental MR³ tight-cluster limit" section was misattributing these failures — they are fixable by this 25-line classification correction, not a fundamental algorithm limit.

---

## 4. Bug #3 — `mr3_gk.py:_dlarnv_normal`: non-reproducible RNG broke Py ↔ Fortran bit-match

### Defect

The Python orchestration's Gram-Schmidt completion path used `np.random.randn(n)` with no seeding contract that the Fortran orchestration could match. The two paths produced different random vectors → different completion bases on degenerate σ = 0 specs → 1-2 ULP differences in U/V columns.

### Fix

Bind LAPACK `DLARNV` (IDIST=3, normal distribution) via ctypes and use deterministic 4-tuple seeds derived from the column index:
```python
_lapack_lib = _ct.CDLL(_ctutil.find_library('lapack'))
def _dlarnv_normal(n, seed):
    iseed = (_ct.c_int * 4)(*seed)
    n_c = _ct.c_int(n); idist = _ct.c_int(3)
    x = np.zeros(n, dtype=np.float64)
    _lapack_lib.dlarnv_(_ct.byref(idist), iseed, _ct.byref(n_c),
                        x.ctypes.data_as(_ct.POINTER(_ct.c_double)))
    return x
```
Used at `mr3_gk.py:969` (V completion, seed `(1,3,5, (2j+1)&4095|1)`) and `mr3_gk.py:981` (U completion, seed `(2,4,6, (2j+3)&4095|1)`).

The Fortran side calls `DLARNV(3, iseed, n, …)` with the **identical** seed sequence (`mr3gk_postproc.f90:371-377, 439-443`).

### Empirical evidence

Before this fix: ~6 / 379 specs differed by 1-2 ULPs between Python and Fortran orchestrations.
After: bit-identical or ≤ 1 ULP. **No impact on the standalone evaluation suite** (Python alone always passes; this is purely cross-implementation harmonization).

---

## 5. Bug #4 — `mr3_gk.py:_system_dnrm2*`: numpy DNRM2 ≠ system DNRM2 broke Py ↔ Fortran bit-match

### Defect

`np.linalg.norm` calls numpy's bundled OpenBLAS. The Fortran orchestration calls system reference BLAS via `-lblas`. The two implementations differ at the ~1 ULP level — enough to break the 2·eps bit-match test on a handful of `gl_wilkw` / `gl_wilk2w` adversarial specs.

### Fix

Bind the system `dnrm2_` symbol via ctypes (same `.so` the Fortran links against) and route every Python normalization through it:
```python
_blas_lib = _ct.CDLL(_ctutil.find_library('blas') or _ctutil.find_library('lapack'))
_blas_lib.dnrm2_.restype = _ct.c_double
_blas_lib.dnrm2_.argtypes = [_ct.POINTER(_ct.c_int),
                              _ct.POINTER(_ct.c_double),
                              _ct.POINTER(_ct.c_int)]
def _system_dnrm2(x): ...
def _system_dnrm2_axis0(M): ...    # column-wise
```
Used at 10 sites in `bidiag_svd`'s post-processing (lines 865-984).

### Empirical evidence

Combined with Bug #3 fix: 379 / 379 of the bit-match test now pass at 2·eps; worst U diff = 1.110e-16, worst V diff = 2.272e-19. **No impact on the standalone evaluation suite**; purely cross-implementation harmonization.

---

## 6. New code — `mr3gk_fortran/`: pure-Fortran orchestration layer

This is **not a bug fix** — it is a re-implementation of the Python wrapper (`mr3_gk.py:bidiag_svd`, ~640 lines) in Fortran so the orchestration layer doesn't depend on Python at all.

| File | Purpose |
|---|---|
| `mr3gk.f90` | Public entry `dmr3gk_svd`: pre-scaling, splitting, sign matrices D1/D2, block dispatch, σ rescale |
| `mr3gk_split.f90` | Two-phase bidiagonal splitting (Phase 1 relative, Phase 2 absolute) |
| `mr3gk_qrsweep.f90` | Zero-shift QR sweep (Demmel-Kahan 1990) for zero-σ deflation |
| `mr3gk_tgk.f90` | T_GK construction + per-block solve via `dlaxre_gk_` + `dlaxrv_` |
| `mr3gk_postproc.f90` | Normalize, Bv recovery, sign fix, GS completion |
| `mr3gk_utils.f90` | DSTEBZ wrapper, vec_norm |
| `mr3gk_consts.f90` | EPS, SAFMIN, GAPTOL, MAX_DEPTH, etc. |
| `mr3gk_run.f90` | CLI driver (binary I/O, used by `test_fortran_match.py`) |

Compile flags `-O0 -fno-fast-math` to keep FP operation order identical to the Python wrapper. See `docs/03_build_flags.md`.

### Empirical evidence

21,375 / 21,375 specs (379 + 1,130 + 19,866 across multiple suites) match the Python+Fortran baseline at 2·eps; σ bit-identical, worst U/V diffs ≤ 1 ULP. See `docs/01_fortran_test_report.md`.

---

## 7. Combined empirical impact

Same 379-spec evaluation, four library variants:

| Variant | `dlaxre` | `dlaxrb_clssfy` | Pass | Fail |
|---|---|---|---|---|
| **patched** (current) | our patch | our patch | **379 / 379** | 0 |
| revert-dlaxre | upstream | our patch | 233 / 379 | 146 |
| revert-clssfy | our patch | upstream | 330 / 379 | 49 |
| revert-both | upstream | upstream | 189 / 379 | 190 |

**The two XMR patches are synergistic:** 13 specs pass under either single revert but fail under double revert. Applying only one is not safe.

---

## 8. Reproducibility

Every claim in this report is reproducible:

```bash
cd python_fortran

# Verify pass-rate of patched build:
bash build.sh
python3 evaluate.py
# → 379/379 PASS

# Verify Py ↔ Fortran bit-match:
cd mr3gk_fortran && bash build.sh && cd ..
python3 test_fortran_match.py
# → PASS: 379/379  Worst sigma diff: 0  Worst U diff: 1.110e-16  Worst V diff: 2.272e-19

# Verify ablation table (Section 7) — ~12 min:
python3 experiments/ablation_xmr_patches.py
# writes experiments/ablation_results.json
```

---

## 9. Changelog

| Date | Change |
|---|---|
| 2026-04 | Bug #1 (dlaxre_gk activate GK branch) |
| 2026-04 | Bug #2 (dlaxrb_clssfy AVGTHRESH guard) |
| 2026-04 | Pure-Fortran orchestration layer (`mr3gk_fortran/`) |
| 2026-04-29 | Bug #3 (DLARNV) and Bug #4 (system DNRM2) for Py ↔ Fo bit-match |
| 2026-04-29 | Empirical ablation harness + this bug report |

---

## 10. Files of interest (with line refs)

**Upstream (read-only, unchanged):**
- `xmr_src/dlaxre.f:269` — disabled GK branch (Bug #1)
- `xmr_src/dlaxrb_clssfy.f:340-361, 414-417, 445-453` — `MIN(0, x)` collapse (Bug #2)

**Our patches:**
- `dlaxre_gk.f:269-275` — Bug #1 fix
- `dlaxrb_clssfy_fix.f:192, 358-366, 414-417, 445-453` — Bug #2 fix
- `mr3_gk.py:18-32` — Bug #3 fix (`_dlarnv_normal`)
- `mr3_gk.py:34-65` — Bug #4 fix (`_system_dnrm2`, `_system_dnrm2_axis0`)
- `mr3gk_fortran/*.f90` — pure-Fortran orchestration

**Build:**
- `build.sh` — links `libxmr.so` from `xmr_src/` + our 2 patched `.f` files
- `mr3gk_fortran/build.sh` — links `mr3gk_run` from `mr3gk_fortran/*.f90` + `libxmr.so`

**Validation:**
- `evaluate.py` — 379-spec MR3-GK SVD evaluation suite (res, ortU, ortV thresholds)
- `test_fortran_match.py` — Py ↔ Fortran bit-identity (2·eps tolerance)
- `experiments/ablation_xmr_patches.py` — 4-variant ablation
- `experiments/cluster_tree_trace.py` — representation-tree NCD trace

**Docs:**
- `docs/01_fortran_test_report.md` — combined 21,375-spec test report
- `docs/02_xmr_modifications.md` — extended version of this report (full diffs + per-failure detail)
- `docs/03_build_flags.md` — compile-flag audit
- `docs/04_ncd_audit.md` — NCD `(max−min)/min` vs `(max−min)/max` audit
- `docs/05_tgk_path.md` — T_GK call graph: XMR vs ours
- `docs/06_cluster_tree_trace.md` — cluster-tree NCD trace + visualization
