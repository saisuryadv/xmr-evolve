# Evaluation Framework Documentation

## Overview

Two components:
- **`src/evaluate.cpp`** — C++ driver that compiles with `src/bidiag_svd.h`, runs all tests, computes metrics, prints results
- **`evaluator.py`** — Python wrapper for OpenEvolve: copies evolved `bidiag_svd.h` to temp dir, compiles, runs evaluate binary, parses metrics, implements cascade stages

## Metrics (all normalized by n·eps)

- **residual**: `max_{i,j} |B(i,j) - (UΣV')(i,j)| / (||B||_max · n · eps)` — where `||B||_max = max_i(|d[i]| + |e[i]|)` (evaluate.cpp:44-76)
- **ortho_u**: `max_{i≤j} |δ_{ij} - (UU')_{ij}| / (n · eps)` — lower triangle of I - UU' (evaluate.cpp:79-96)
- **ortho_v**: same as ortho_u but for V

## Pass Thresholds

```cpp
const double RES_THRESH = 7.0;   // evaluate.cpp:897
const double ORTHO_THRESH = 5.0; // evaluate.cpp:898
```

A test passes iff: `info == 0 && residual ≤ 7.0 && ortho_u ≤ 5.0 && ortho_v ≤ 5.0`.

## Data Layout Convention

`bidiag_svd()` returns `BidiagSVDResult` with column-major matrices:
- `sigma[k]` = k-th singular value (descending order)
- `U[k*n + i]` = i-th component of k-th left singular vector
- `V[k*n + j]` = j-th component of k-th right singular vector
- Relationship: `B(i,j) = Σ_k U[k*n+i] · σ_k · V[k*n+j]`

**Important**: DBDSQR returns VT (transpose of V), must be transposed. DBDSVDX returns U in rows 0..n-1 and V in rows n..2n-1 of Z, in ascending singular value order.

## Test Sizes

Built from `{10, 100, adv_size, adv_size*2}`, deduplicated and sorted (evaluate.cpp:929-933).

| adv_size | Resulting sizes | Adversarial tests | + STColl | Total |
|----------|-----------------|-------------------|----------|-------|
| 50 | {10, 50, 100} | 90×3 = 270 | 0 | 270 |
| 100 | {10, 100, 200} | 90×3 = 270 | 19 | 289 |
| 200 | {10, 100, 200, 400} | 90×4 = 360 | 19 | 379 |

## Complete Pattern List (90 patterns)

### Group 1: Original Adversarial (22 patterns)
| # | Pattern | Description |
|---|---------|-------------|
| 1 | exponential_graded | d[i]=10^(-16i/n), e=d/2. κ~10^16 |
| 2 | glued_repeated | blocks of 10 with d≈1, e≈0.5, glue=1e-15 |
| 3 | saw_tooth | d alternates 1.0/1e-8, e=0.5 |
| 4 | stemr_killer | d[i]=10^(-20i/n), e=√(d[i]·d[i+1]) |
| 5 | huge_condition | d[i]=10^(-15i/(n-1)), e=√(d·d)·0.9 |
| 6 | spike | d=0.01 except d[n/2]=100, e=0.01 |
| 7 | wilkinson_like | d[i]=\|i-n/2\|+1, e=1 |
| 8 | two_clusters | d[0..n/2]=1, d[n/2..n]=1e-8, e=1e-10 |
| 9 | random_uniform | random d,e in (0,1) |
| 10 | diagonal_only | d random, e=0 |
| 11 | constant | d=1, e=0.5 |
| 12 | all_equal_nontrivial | d=1, e=1 |
| 13 | one_big_cluster | d≈1+O(1e-12), e≈O(1e-12) |
| 14 | arithmetic_progression | d[i]=1+i/(n-1) |
| 15 | many_near_zero | d≈1e-15, e≈1e-15 except a few |
| 16 | random_dense_clusters | random cluster centers, tight spacing |
| 17 | constant_d_graded_e | d=1, e[i]=10^(-16i/(n-1)) |
| 18 | random_clustered_5 | 5 cluster centers, d random within eps |
| 19 | alternating_sign | d alternates +/- |
| 20 | step_function | d steps from 1 to 1e-8 at midpoint |
| 21 | three_clusters | 3 separated clusters |
| 22 | random_sparse_e | random d, most e=0 |

### Group 2: Marques 2020 (2 patterns)
| # | Pattern | Description |
|---|---------|-------------|
| 23 | chkbd | d[i]=10^(-(2i+1)), e[i]=10^(-2(i+1)). DBDSVDX reortho bug trigger |
| 24 | marques_graded | d[i]=√eps^(i/(n-1)), e=√(d·d)·0.3 |

### Group 3: Grosser-Lang Entry-Based (17 patterns)
| # | Pattern | Description |
|---|---------|-------------|
| 25-28 | gl_abcon0–gl_abcon3 | ABCON heat equation tridiagonal (α=2, β=1) |
| 29 | gl_random | random bidiagonal entries |
| 30-31 | gl_gradp, gl_gradm | graded diagonal (increasing/decreasing) |
| 32-35 | gl_wilkp, gl_wilkm, gl_wilkw, gl_wilk2w | Wilkinson-type variants |
| 36 | gl_clement | Clement matrix |
| 37-40 | gl_gro0–gl_gro3 | near-constant with outliers |
| 41 | gl_bwilkp | blocked Wilkinson+ |

### Group 4: Grosser-Lang Spectrum-Based (10 patterns)
| # | Pattern | Description |
|---|---------|-------------|
| 42 | gl_ones | all-ones spectrum |
| 43 | gl_uniform_eps | uniform, eps-apart |
| 44 | gl_uniform_sqrteps | uniform, √eps-apart |
| 45 | gl_uniform_eps_to_1 | uniform from eps to 1 |
| 46 | gl_geometric_eps_to_1 | geometric from eps to 1 |
| 47 | gl_random_spectrum | random spectrum |
| 48-49 | gl_clustered_at_1, gl_clustered_at_pm1 | clustered at 1 / ±1 |
| 50-51 | gl_clustered_at_eps, gl_clustered_at_pmeps | clustered at eps / ±eps |

### Group 5: Demmel 2008 Strongly Clustered (4 patterns)
| # | Pattern |
|---|---------|
| 52-53 | demmel_S1pe, demmel_S1ps |
| 54-55 | demmel_S2pe, demmel_S2ps |

### Group 6: Demmel 2008 Weakly Clustered (3 patterns)
| # | Pattern |
|---|---------|
| 56-58 | demmel_W1, demmel_W2, demmel_W3 |

### Group 7: Demmel 2008 Geometric/Uniform (4 patterns)
| # | Pattern |
|---|---------|
| 59-60 | demmel_G1, demmel_G1s |
| 61-62 | demmel_U1, demmel_U1s |

### Group 8: Exact Wilkinson / Glued Wilkinson (3 patterns)
| # | Pattern |
|---|---------|
| 63 | wilkinson_exact |
| 64-65 | glued_wilkinson, glued_wilkinson_tight |

### Group 9: Willems-Lang (1 pattern)
| # | Pattern |
|---|---------|
| 66 | wl_example48 |

### Group 10: Parlett-Dhillon Diagnostics (2 patterns)
| # | Pattern |
|---|---------|
| 67-68 | pd_T0, pd_T1 |

### Group 11: Stress Tests (6 patterns)
| # | Pattern |
|---|---------|
| 69 | zero_diagonal |
| 70 | single_element |
| 71-72 | near_overflow, near_underflow |
| 73 | mixed_signs |
| 74 | checkerboard |

### Group 12: Condition Number Variants (16 patterns)
Same structural patterns at reduced κ to match Willems-Lang Synth methodology:
| # | Pattern |
|---|---------|
| 75-77 | exponential_graded_k4, _k8, _k12 |
| 78-79 | stemr_killer_k5, _k10 |
| 80-81 | huge_condition_k5, _k10 |
| 82-84 | demmel_G1_k4, _k8, _k12 |
| 85-86 | demmel_S1pe_k4, _k8 |
| 87-88 | chkbd_4, chkbd_16 |
| 89-90 | marques_graded_k4, _k8 |

## Scoring Formula

### evaluate.cpp scoring (evaluate.cpp:1119-1128)

```
score = pass_rate × 50                              // up to 50 pts
      + max(0, 5 - pass_avg_res) × 2.0              // up to 10 pts
      + max(0, 5 - pass_avg_ortU) × 2.0             // up to 10 pts
      + max(0, 5 - pass_avg_ortV) × 2.0             // up to 10 pts
      + 5.0                                          // compilation bonus
      + scaling_score                                // up to 10 pts
```

Where `scaling_score`:
- `worst_ratio ≤ 5.0` → 10.0
- `5.0 < worst_ratio ≤ 10.0` → `10.0 × (10.0 - worst_ratio) / 5.0`
- `worst_ratio > 10.0` → 0.0

Max from evaluate.cpp alone: **95 pts**.

### evaluator.py additional scoring (evaluator.py:189-199)

If stage 2 `pass_rate ≥ 0.6`, evaluator.py runs stage 3 and adds:
- `large_pass_rate × 10.0` — up to 10 pts
- Large scaling bonus — up to 5 pts (same formula as above, scaled to 5)

Max from evaluator.py: **110 pts**.

## Evaluator Cascade (evaluator.py)

| Stage | adv_size | STColl | Total Tests | Timeout | Gate |
|-------|----------|--------|-------------|---------|------|
| 1 | 50 | No | 270 | 30s | Always runs |
| 2 | 100 | Yes | 289 | 120s | Always runs (in `evaluate()`) |
| 3 | 200 | Yes | 379 | 600s | Stage 2 pass_rate ≥ 0.6 |

## Scaling Test

All patterns that pass at **every** test size are included. The worst time doubling ratio from the **largest size transition only** (e.g., n=200→400 for stage 3) is used. This is NOT limited to a fixed set of 5 patterns — all passing patterns contribute (evaluate.cpp:1025-1043).

Target: worst ratio ≤ 5.0 for O(n²). O(n²) → ratio ≈ 4.0, O(n³) → ratio ≈ 8.0.

## Available Routines (via extern "C")

All from `lib/libxmr_c.a` — **pure C** (CLAPACK 3.2.1 + f2c-converted XMR/hgbsvd). No Fortran compiler, no Apple Accelerate.

| Routine | Source | Purpose |
|---------|--------|---------|
| dbdsqr_ | CLAPACK | QR iteration bidiag SVD (O(n³)) |
| dstemr_ | CLAPACK | MR³ tridiag eigensolver (O(n²)) |
| dstexr_ | XMR (f2c) | Willems XMR improved MR³ (O(n²)) |
| dbdsgr_ | hgbsvd (f2c) | Großer-Lang coupling SVD (O(n²)) |
| dlamch_ | CLAPACK | Machine parameters (eps, safmin) |
| BLAS | CLAPACK | dnrm2, ddot, dscal, dcopy, daxpy, dgemv, etc. |

**No DBDSVDX** — it requires full LAPACK (DSTEVX, DSTEIN), not in libxmr_c.a.

## Fortran Calling Convention

`dbdsgr_` (f2c-generated from Fortran with CHARACTER*1 args) requires `ftnlen` arguments at the end of each call. Pass `1, 1` for JOBZ and RANGE character lengths:

```cpp
dbdsgr_("V", "A", &n, d, e, &vl, &vu, &il, &iu, &abstol, &m, w,
        z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info,
        /*ftnlen*/ 1, 1);
```

`dstemr_` and `dstexr_` do **NOT** need ftnlen args (their f2c wrappers handle this internally).

## Build Command

```bash
g++ -std=c++17 -O2 -Isrc/clapack -o evaluate src/evaluate.cpp -Llib -lxmr_c -lm
```

No Fortran compiler needed. All routines are pre-converted C.
