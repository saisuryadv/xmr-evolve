# Baseline Algorithm Comparison

All baselines compiled and run as **pure C/C++** — no Fortran, no gfortran, no Apple Accelerate.

## Test Configuration

- **90 adversarial patterns** × **4 sizes** (n=10, 100, 200, 400) + **19 STCollection** = **379 tests** (stage 3)
- **Pass thresholds** (in units of n·eps): residual ≤ **7.0**, orthoU ≤ **5.0**, orthoV ≤ **5.0**
- **Scaling** = worst time doubling ratio (largest size transition only: n=200→400) across patterns that **pass at all sizes**
- **Avg** = over all tests (per-value capped at 1e6)

## Evaluator Stages

| Stage | adv_size | Sizes | STColl | Total Tests | Timeout |
|-------|----------|-------|--------|-------------|---------|
| 1 | 50 | {10, 50, 100} | No | 270 | 30s |
| 2 | 100 | {10, 100, 200} | Yes | 289 | 120s |
| 3 | 200 | {10, 100, 200, 400} | Yes | 379 | 600s |

Stage 3 only runs if stage 2 pass_rate ≥ 0.6.

## Stage 3 Results (Full Run — 379 tests)

Scaling measured on patterns that **pass at all sizes**. Worst ratio used for scoring.

| Algorithm | Pass | Avg Res | Avg OrtU | Avg OrtV | Worst Pass Scaling | #Pass Patterns | Score |
|-----------|------|---------|----------|----------|--------------------|----------------|-------|
| **DBDSQR** | **379/379** | 0.43 | 0.37 | 0.35 | 8.53x (O(n³)) | 90/90 | **85.64** |
| **DBDSVDX** | 371/379 | 13193 | 15852 | 13213 | 13.21x (O(n³)) | 87/90 | 53.94 |
| **HGBSVD** | 174/379 | 533121 | 532982 | 533112 | 4.35x (O(n²)) | 36/90 | 37.96 |
| **TGK+DSTEMR** | 93/379 | 52935 | 492760 | 451378 | 4.82x (O(n²)) | 18/90 | 27.27 |
| **TGK+DSTEXR** | 61/379 | 163642 | 567171 | 520277 | 5.73x | 11/90 | 21.58 |

## Top-5 Worst Pass Scaling (Stage 3, n=200→400 only)

### DBDSQR (O(n³) confirmed)
| Ratio | Pattern |
|-------|---------|
| 8.53x | diagonal_only |
| 7.94x | gl_uniform_sqrteps |
| 7.77x | demmel_W3 |
| 7.52x | demmel_G1_k12 |
| 7.39x | gl_abcon0 |

### DBDSVDX (O(n³) confirmed)
| Ratio | Pattern |
|-------|---------|
| 13.21x | demmel_W3 |
| 12.86x | demmel_S2ps |
| 12.83x | gl_clustered_at_pm1 |
| 12.74x | demmel_S2pe |
| 12.66x | gl_clustered_at_1 |

### HGBSVD (O(n²) confirmed)
| Ratio | Pattern |
|-------|---------|
| 4.35x | stemr_killer_k5 |
| 4.31x | stemr_killer_k10 |
| 4.25x | huge_condition_k5 |
| 4.22x | gl_abcon0 |
| 4.19x | exponential_graded_k12 |

### TGK+DSTEMR (O(n²) confirmed)
| Ratio | Pattern |
|-------|---------|
| 4.82x | constant |
| 4.76x | random_uniform |
| 4.73x | mixed_signs |
| 4.57x | gl_abcon0 |
| 4.55x | gl_ones |

### TGK+DSTEXR
| Ratio | Pattern |
|-------|---------|
| 5.73x | constant |
| 4.99x | glued_repeated |
| 4.95x | diagonal_only |
| 4.89x | gl_ones |
| 4.58x | gl_random_spectrum |

## Algorithm Details

### DBDSQR (reference baseline, 379/379)
- **Method**: LAPACK QR iteration for bidiagonal SVD. Implicit zero-shift QR with convergence tests.
- **Complexity**: O(n²) typical (few iterations), O(n³) worst case (n iterations for convergence). Empirically O(n³) on adversarial patterns — worst ratio 8.53x.
- **Accuracy**: Perfect on all 379 matrices. Gold standard: avg res=0.43, ortU=0.37, ortV=0.35.
- **Why it's reliable**: Direct bidiagonal algorithm — no reduction to tridiagonal, no extraction step. U and V updated simultaneously per QR step, maintaining coupling.
- Source: `src/bidiag_dbdsqr.h`

### DBDSVDX (371/379)
- **Method**: LAPACK bidiagonal SVD via DSTEVX (bisection + inverse iteration) on the 2n×2n TGK tridiagonal. After DSTEIN computes TGK eigenvectors, DBDSVDX extracts u/v parts.
- **Complexity**: O(n²) for eigenvalues (bisection), O(n³) worst case for eigenvectors (DSTEIN reortho). Empirically O(n³) — worst ratio 13.21x.
- **Failures**: 8 tests fail — all from chkbd and checkerboard families. Root cause: DBDSVDX's post-extraction reorthogonalization trigger is norm-based (`ABS(NRMU-ORTOL)*SQRT2.GT.ONE`) and never fires for these matrices even though orthogonality is O(10^13). This is INDEX.md Bug #7, unfixed in LAPACK 3.12.1.
- **Not the same as DSTEIN's bug**: DSTEIN itself uses eigenvalue-separation-based reortho and works correctly. The bug is in DBDSVDX's own post-processing.
- Source: `src/bidiag_dbdsvdx.h`

### HGBSVD (174/379)
- **Method**: Großer-Lang (2001) coupling-based bidiagonal SVD via DBDSGR. Works on B^TB tridiagonal using MR³ (dstemr_ internally) with coupling transforms to simultaneously compute U and V. Avoids TGK entirely.
- **Complexity**: O(n²) on passing matrices (worst 4.35x).
- **Accuracy on passing tests**: Excellent — pass max res=2.18, ortU=4.0.
- **Failures**: 205/379 fail with INFO!=0. DBDSGR returns INFO=3/4/5 for deep recursion, tight clusters, or MGS-requiring clusters. Only positive definite couplings implemented. Splittings not supported (INFO=33). DLARRI (eigenvalue refinement) missing from XMR distribution.
- **f2c calling convention**: Requires `ftnlen` args `1, 1` at end for CHARACTER*1 parameters.
- Source: `src/bidiag_hgbsvd.h`

### TGK+DSTEMR (93/379)
- **Method**: LAPACK MR³ (DSTEMR) on the 2n×2n TGK tridiagonal, then extract U from odd rows (1,3,5,...) and V from even rows (0,2,4,...) of eigenvectors.
- **Complexity**: O(n²) confirmed (worst 4.82x on passing matrices).
- **Failures**: 75% fail. Eigenvector extraction destroys orthogonality because DSTEMR's representation tree doesn't know about GK structure — shifts can cluster +σ and -σ eigenvalues, producing eigenvectors whose u/v parts lack the equal-norm property.
- **Post-processing**: Normalization, sign consistency, one-sided recovery (U=BV/σ for small σ), chunked MGS reortho. Helps some cases but can't fix fundamental extraction failures.
- Source: `src/bidiag_tgk_stemr.h`

### TGK+DSTEXR (61/379)
- **Method**: Willems XMR improved MR³ (DSTEXR) on TGK. Same extraction as DSTEMR.
- **Complexity**: Mostly O(n²), one outlier at 5.73x (`constant` pattern).
- **Failures**: Worse pass rate than DSTEMR despite theoretically better eigensolver. GK-form support deactivated in dlaxre.f (INDEX.md Bug #10). AVGAPFAC bug fixed (INDEX.md Bug #1) but graded-spectrum failures persist (Bug #2).
- Source: `src/bidiag_tgk_stexr.h`

### DSTEGR = DSTEMR
DSTEGR is just a thin wrapper that calls DSTEMR. They are identical.

## The Fundamental Tradeoff

| Property | DBDSQR/DBDSVDX | HGBSVD | TGK+STEMR/STEXR |
|----------|----------------|--------|------------------|
| Accuracy | Excellent | Good (when it works) | Poor |
| Pass Rate | 97-100% | 46% | 16-25% |
| Scaling | O(n³) | O(n²) | ~O(n²) |

**The goal**: achieve DBDSQR-level accuracy (379/379) with O(n²) scaling (≤5x ratio). This is an open problem.

## Build & Link (Pure C/C++)

```bash
g++ -std=c++17 -O2 -Isrc/clapack -o evaluate src/evaluate.cpp -Llib -lxmr_c -lm
```

No Fortran compiler needed. All LAPACK/BLAS routines are pre-converted C from CLAPACK 3.2.1.
XMR (DSTEXR) and hgbsvd (DBDSGR) converted from Fortran via f2c.

## Reference Implementation Files

All in `src/`:
- `bidiag_dbdsqr.h` — DBDSQR wrapper (O(n³) baseline)
- `bidiag_dbdsvdx.h` — DBDSVDX wrapper (f2c from LAPACK 3.12.1)
- `bidiag_hgbsvd.h` — HGBSVD/DBDSGR wrapper (Großer-Lang, O(n²))
- `bidiag_tgk_stemr.h` — TGK + DSTEMR with post-processing
- `bidiag_tgk_stexr.h` — TGK + DSTEXR with post-processing
- `bidiag_tgk_common.h` — Shared: TGK construction, U/V extraction, normalization, sign fix, one-sided recovery, chunked MGS reorthogonalization
- `bidiag_svd.h` — **THE FILE BEING EVOLVED** (hybrid HGBSVD + TGK+DSTEMR fallback)
- `fortran_interface.h` — C++ extern "C" declarations for all routines
