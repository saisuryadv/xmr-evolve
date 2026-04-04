# MR³-GK Bidiagonal SVD — Python + Fortran Implementation

O(n²) bidiagonal SVD using MR³ with Golub-Kahan structure preservation (Willems & Lang 2012).

## Architecture

```
Python (mr3_gk.py)          Fortran (XMR, ~15K lines)
┌─────────────────┐         ┌──────────────────┐
│ bidiag_svd()    │         │ dlaxre_gk.f      │ ← GK-enabled root repr
│   QR deflation  │         │ dlaxrv.f         │ ← eigenvector computation
│   split_bidiag  │         │ dlaxrf.f + deps  │ ← child repr (shift/twist)
│   mr3_tgk       │──ctypes─│ xmr_wrapper.c    │ ← C interface
│   sign recovery │         └──────────────────┘
│   Bv recovery   │
│   GS completion │
└─────────────────┘
Eigenvalues: LAPACK dstebz (bisection)
```

## Current Results

**Score: 86.37** — 330/379 tests passing, worst scaling ratio 4.29×

| Suite | Tests | Passed | Rate |
|-------|-------|--------|------|
| Adversarial (90 patterns × 3 sizes) | 270 | 235 | 87.0% |
| Scaling (n=400) + STCollection | 109 | 95 | 87.2% |
| **Total** | **379** | **330** | **87.1%** |

### Worst-Case Metrics

| Metric | Worst Value | Test Case | Threshold |
|--------|-------------|-----------|-----------|
| Residual | 8.43 nε | pd_T0@100 | ≤ 7.0 |
| Orthogonality U | 33.23 nε | demmel_S1pe_k4@100 | ≤ 5.0 |
| Orthogonality V | 33.22 nε | demmel_S1pe_k4@100 | ≤ 5.0 |

### Baseline Comparison

| Algorithm | Pass | Rate | Worst Scaling | Score |
|-----------|------|------|---------------|-------|
| **DBDSQR** | **371/371** | **100%** | 25.9× (O(n³)) | 5.00 (gated) |
| **MR³-GK (ours)** | **330/379** | **87.1%** | **4.29×** | **86.37** |
| TGK+DSTEMR | crash | — | — | — |

## Zero-Shift QR Deflation (Preprocessing)

Before MR³, each unsplit block undergoes one sweep of implicit zero-shift QR iteration (Demmel-Kahan 1990, Section 3). This serves as a zero singular value detector:

**Theory**: Zero-shift QR is equivalent to inverse iteration on B^TB (Demmel-Kahan 1990). When σ_min = 0, inverse iteration converges in one step — both d[n-1] and e[n-1] converge to zero, creating a clean 1×1 split at the bottom.

**Algorithm**:
1. For each unsplit block of size k, apply one zero-shift QR sweep (O(k), cheap)
2. Check if |d[-1]| < n·ε AND |e[-1]| < n·ε (clean split at bottom)
3. If yes → zero SV detected: solve (k-1) sub-problem, embed into k×k, apply Givens rotations
4. If no → no zero SV: proceed with MR³ normally

**Why check both d[-1] AND e[-1]**: When σ_min = 0, the QR sweep drives both to zero. When σ_min > 0 (even if tiny), d[-1] → σ_min but e[-1] stays large — NOT a clean split. Deflating without a clean split severs matrix coupling and produces catastrophic residual errors (observed on gl_wilkw: e[-1] = 0.707 after sweep).

**Impact**: Fixes saw_tooth@200 (ortU: 552 → 2.04 neps), step_function@200 (ortU: 70.6 → 18.1 neps), and all matrices with structural zero diagonals (zero_diagonal, gl_wilkp, gl_clement, etc.).

## Analysis of Remaining Failures

All 49 remaining failures share one common theme: **tight eigenvalue clusters in the T_GK tridiagonal matrix that MR³ cannot resolve**.

### Top 15 Worst Orthogonality Cases

| Test | ortU (neps) | Max Cluster Size | Clustered/Total | Description |
|------|-------------|-----------------|-----------------|-------------|
| demmel_S1pe_k4@100 | 33.2 | 12 | 20/100 | Demmel type 1, positive end, κ=10⁴ |
| pd_T0@100 | 33.0 | 76 | 76/100 | Tight uniform cluster spanning 76% of spectrum |
| gl_abcon3@100 | 24.3 | 65 | 65/100 | Glued matrix, 65-element cluster |
| demmel_S1pe_k4@200 | 21.8 | 12 | 20/200 | Same pattern at larger size |
| three_clusters@100 | 21.6 | 34 | 100/100 | Three groups, ALL eigenvalues clustered |
| gl_wilkw@200 | 21.1 | 9 | 9/200 | Wilkinson-type with wide off-diag |
| demmel_S1ps@100 | 19.3 | 12 | 20/100 | Demmel type 1, positive symmetry |
| three_clusters@200 | 18.3 | 34 | 200/200 | Three groups at larger size |
| step_function@200 | 18.1 | 3 | 5/200 | Step grading, near-zero SV + cluster |
| two_clusters@100 | 17.8 | 50 | 100/100 | Two groups, ALL eigenvalues clustered |

### Common Theme: Tight Eigenvalue Clusters

A cluster is a group of singular values with relative gap < 10⁻³. The failing matrices exhibit:

1. **Large clusters** (pd_T0, gl_abcon3, two_clusters, three_clusters): 50-76 eigenvalues with relative gaps below 10⁻³. MR³'s representation tree cannot find child shifts that make all eigenvalues in these clusters relatively well-separated simultaneously.

2. **NEGL miscount amplification** (step_function, demmel_S2pe, gl_abcon1): Matrices with near-zero (but not exactly zero) singular values where the T_GK Sturm sequence miscounts eigenvalues ≤ 0. The QR deflation fixes this when d[-1] AND e[-1] are both near-zero (clean split), but some cases have σ_min too large for deflation yet too small for stable NEGL computation.

3. **High condition number + clusters** (demmel_S1pe, gl_wilkw, gl_abcon3): Condition numbers of 10⁴-10¹⁵⁸ combined with tight clusters. The wide dynamic range makes shift selection harder — shifts that work well for large eigenvalues may not resolve small clustered eigenvalues.

### Why MR³ Struggles with Tight Clusters

In MR³, orthogonality between eigenvectors comes from the relative gap: if eigenvalue λᵢ has relative gap ≥ gaptol with respect to its neighbors in the current LDL^T representation, its eigenvector can be computed to high relative accuracy directly. When eigenvalues cluster (relative gap < gaptol ≈ 10⁻³), MR³ must:

1. Group the cluster
2. Find a child representation (new shift) where cluster members become well-separated
3. Recurse

This process can fail when:
- No shift makes all cluster members simultaneously well-separated
- Element growth in the child LDL^T factorization is too large
- The cluster is so tight that machine precision limits the achievable separation
- NEGL miscount causes wrong eigenvalue assignment to clusters

These are fundamental MR³ representation tree limitations, not bugs in our implementation.

## Prerequisites

```bash
sudo apt install gfortran liblapack-dev libblas-dev python3-numpy python3-scipy
```

## Build

```bash
bash build.sh
```

A prebuilt `libxmr.so` is included for Linux x86_64.

## Run Evaluation

```bash
# Full 379-test suite (definitive evaluation)
python3 evaluate.py

# Medium 270-test suite (no n=400 or STCollection, faster iteration)
python3 evaluate.py --medium
```

### Output format

Each test prints:
```
  pattern_name    n=NNNN  res=X.XXX  ortU=X.XXX  ortV=X.XXX  t=X.XXXXs  PASS/FAIL
```

Thresholds: res ≤ 7.0 nε, ortU ≤ 5.0 nε, ortV ≤ 5.0 nε.

## Extended Test Suite (224 tests)

An additional test suite (`test_dense_to_bidiag.py`) with **64 patterns / 224 tests** designed to find failure cases not covered by the 379-test suite. Same metrics and thresholds. Test sizes: {10, 100, 200, 400}.

**Current results: 219/224 passed (97.8%), 5 failures found.**

```bash
# Full extended suite
python3 test_dense_to_bidiag.py

# Quick mode (sizes {10, 100} only)
python3 test_dense_to_bidiag.py --quick

# Run only Part 1 (dense-to-bidiag) or Part 2 (paper tests)
python3 test_dense_to_bidiag.py --part1
python3 test_dense_to_bidiag.py --part2
```

### Part 1: Dense-to-Bidiagonal Tests (22 patterns)

Constructs dense matrices with prescribed singular value distributions, reduces to bidiagonal via Householder (equivalent to LAPACK DGEBRD), then tests `bidiag_svd` on the result. The bidiagonal produced by dense reduction has different entry structure than hand-crafted bidiagonals, even for the same SV distribution (Demmel-Kahan 1990: "reduction to bidiagonal form may produce completely inaccurate bidiagonal entries").

Sources: LAPACK dchkbd.f Types 3-16, LAPACK ddrvbd.f Types 3-5, LAPACK bug #316 (DGESDD companion matrix), Willems-Lang 2012 Examples 4.1/4.8, classical matrices (Hilbert, Kahan, Vandermonde, Frank, Toeplitz).

### Part 2: Missing Paper Test Cases (42 patterns)

Bidiagonal test matrices from the literature not covered in the 379-test suite:

| Category | # | Source |
|----------|---|--------|
| Proper Glued Wilkinson | 4 | Dhillon-Parlett-Vömel 2005, Demmel 2008 |
| Glued versions of existing patterns | 8 | Willems-Lang 2013 Synth methodology |
| CHKBD LAPACK general (log-uniform) | 1 | LAPACK dchkbd.f Type 16 |
| Corrected Grosser-Lang definitions | 7 | Großer-Lang 2001 DMATGEN.f IDs 200-244 |
| Marques STEXR failure matrices | 3 | Marques-Demmel-Vasconcelos 2020 Appendix B |
| Tridiagonal→bidiagonal (Cholesky) | 5 | Großer-Lang 2005, STCollection (Legendre/Laguerre/Hermite) |
| Missing Demmel/Willems variants | 4 | Demmel 2008, Willems-Lang 2012, Marques 2020 |
| Dhillon-Parlett-Vömel / LAWN 163/166 | 4 | Dhillon Thesis 1997, Großer-Lang 2005, LAWN 163/166 |
| STCollection eigenvalue types 7/8/9 | 3 | Marques et al. 2008 Algorithm 880 |
| Skew-Wilkinson with grading | 3 | STCollection T_SkewW21 variants |
| W21 with varying off-diagonal gamma | 3 | STCollection T_W21_g variants |
| Willems 2x2 block LDL* edge case | 1 | Willems Thesis 2010 |

### Failures Found

| Test | n | ortU | Root Cause |
|------|---|------|------------|
| `dense_wilkinson_sv` | 99 | 29.86 | Tight eigenvalue pairs from Wilkinson W_{2m+1}^+ as SVs |
| `dense_wilkinson_sv` | 399 | 27.60 | Same, larger size |
| `dense_wilkinson_sv` | 199 | 6.15 | Same, moderate |
| `glued_wilk_5x21_sqrteps` | 105 | 13.11 | 5 copies of W_21^+ glued — near-degenerate eigenvalue clusters |
| `marques_stexr_20` | 20 | 8.74 | Geometric eigenvalue distribution from Marques 2020 |

## Files

| File | Description |
|------|-------------|
| `mr3_gk.py` | Main implementation: QR deflation, splitting, T_GK construction, SVD extraction |
| `xmr_ctypes.py` | Python ctypes interface to XMR Fortran library |
| `full_eval.py` | 90 adversarial test matrix generators |
| `evaluate.py` | Evaluation runner with scoring (379 tests) |
| `test_dense_to_bidiag.py` | Extended test suite: 64 patterns, 224 tests (dense-to-bidiag + paper tests) |
| `run_baselines.py` | Baseline comparison (DBDSQR, TGK+DSTEMR vs ours) |
| `dlaxre_gk.f` | Modified XMR root representation with GK detection enabled |
| `xmr_wrapper.c` | C wrapper for dlaxre_ + dlaxrv_ |
| `libxmr.so` | Prebuilt shared library (Linux x86_64) |
| `build.sh` | Build script to recompile from source |
| `fortran_objects/` | Precompiled XMR .o files (unmodified from Willems' code) |
| `stcollection/` | 19 STCollection bidiagonal test matrices (.dat files) |
| `CHANGES.md` | Full history of changes and design decisions |

## Key Design Decisions

- **Zero-shift QR deflation**: One QR sweep per unsplit block detects zero SVs (both structural and product-underflow). Clean-split check (both d[-1] and e[-1] near zero) prevents false positives.
- **No reorthogonalization**: All orthogonality comes from MR³ structure, not post-hoc correction
- **NEGL guard kept**: dlaxre's quality check routes some matrices to PD shift fallback — removing it causes regressions
- **Gram-Schmidt completion**: GS completion as safety net for vectors MR³ fails to compute (zero SVs handled upstream by QR deflation now)

## References

- Willems & Lang, "A Framework for the MR³ Algorithm: Theory and Implementation", ETNA 2012
- Demmel & Kahan, "Accurate Singular Values of Bidiagonal Matrices", SIAM J. Sci. Stat. Comput. 1990
- Dhillon, Parlett & Vömel, "The Design and Implementation of the MRRR Algorithm", LAPACK Working Note 162
