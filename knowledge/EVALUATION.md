# Evaluation Framework Documentation

## Overview

`src/evaluate.cpp` compiles with `src/bidiag_svd.h` and runs 270 adversarial tests + scaling tests.
`evaluator.py` wraps this for OpenEvolve (compile → run → parse metrics).

## Metrics (all normalized by n·eps)

- **residual**: ||B - U·Σ·V'||_max / (||B|| · n · eps)
- **ortho_u**: ||I - U·U'||_max / (n · eps)
- **ortho_v**: ||I - V·V'||_max / (n · eps)

Pass threshold: residual ≤ 7.0, orthoU ≤ 5.0, orthoV ≤ 5.0.

## Data Layout Convention

`bidiag_svd()` returns `BidiagSVDResult` with column-major matrices:
- `sigma[k]` = k-th singular value (descending order)
- `U[k*n + i]` = i-th component of k-th left singular vector
- `V[k*n + j]` = j-th component of k-th right singular vector
- Relationship: B(i,j) = Σ_k U[k*n+i] · σ_k · V[k*n+j]

**Important**: DBDSQR returns VT (transpose of V), must be transposed. DBDSVDX returns U in rows 0..n-1 and V in rows n..2n-1 of Z, in ascending singular value order.

## Scoring (0-100 composite)

- 5 pts: compiles successfully
- 50 pts: pass rate × 50 (fraction of 270 tests passing)
- 15 pts: accuracy bonus (lower avg residual + ortho = more points)
- 10 pts: O(n²) scaling (worst ratio ≤ 5.0 → full 10 pts)
- 10 pts: large adversarial pass rate (separate from main 270)
- 10 pts: bonus for very good accuracy (all metrics < 2)

## Test Matrix Groups (270 total = 90 patterns × 3 sizes)

Each pattern tested at n=10, n=100, n=200 (3 tests per pattern).

### Group 1: Original Adversarial (22 patterns)
| Pattern | Description | Why Hard |
|---------|-------------|----------|
| exponential_graded | d[i]=10^(-16i/n), e=d/2 | κ~10^16, graded spectrum |
| glued_repeated | blocks of 10 with d≈1, e≈0.5, glue=1e-15 | repeated eigenvalue clusters |
| saw_tooth | d alternates 1.0 / 1e-8, e=0.5 | B^TB well-conditioned but d graded |
| stemr_killer | d[i]=10^(-20i/n), e=√(d[i]·d[i+1]) | κ~10^20, geometric mean coupling |
| huge_condition | d[i]=10^(-15i/(n-1)), e=√(d·d)·0.9 | κ~10^15, near-geometric coupling |
| spike | d=0.01 except d[n/2]=100, e=0.01 | one huge outlier singular value |
| wilkinson_like | d[i]=|i-n/2|+1, e=1 | Wilkinson eigenvalue pairs |
| two_clusters | d[0..n/2]=1, d[n/2..n]=1e-8, e=1e-10 | two scale-separated clusters |
| all_equal_nontrivial | d=1, e=1 | all singular values equal, max clustering |
| one_big_cluster | d≈1+O(1e-12), e≈O(1e-12) | one huge cluster of nearly equal σ |
| constant_d_graded_e | d=1, e[i]=10^(-16i/(n-1)) | graded off-diagonal only |
| random_clustered_5 | 5 cluster centers, d random within eps | mixed tight clusters |
| checkerboard | d alternates 1/eps, e=√(d·d) | alternating scale with coupling |
| ... | (plus diagonal_only, constant, random_uniform, etc.) | |

### Group 2: Marques 2020 (2 patterns)
| Pattern | Description | Why Hard |
|---------|-------------|----------|
| chkbd | d[i]=10^(-(2i+1)), e[i]=10^(-2(i+1)) | DBDSVDX reortho bug trigger |
| marques_graded | d[i]=√eps^(i/(n-1)), e=√(d·d)·0.3 | Graded from 1 to √eps |

### Group 3: Grosser-Lang Entry-Based (16 patterns)
gl_abcon0-3, gl_random, gl_gradp/m, gl_wilkp/m/w/2w, gl_clement, gl_gro0-3, gl_bwilkp

### Group 4: Grosser-Lang Spectrum-Based (10 patterns)
gl_ones, gl_uniform_eps/sqrteps/eps_to_1, gl_geometric_eps_to_1, gl_random_spectrum, gl_clustered_at_1/pm1/eps/pmeps

### Group 5-7: Demmel 2008 (10 patterns)
demmel_S1pe/S1ps, demmel_S2pe/S2ps, demmel_W1/W2/W3, demmel_G1/G1s, demmel_U1/U1s

### Group 8-11: Wilkinson, Diagnostic, Stress (12 patterns)
wilkinson_exact, glued_wilkinson/tight, wl_example48, pd_T0/T1, zero_diagonal, single_element, near_overflow/underflow, mixed_signs, checkerboard

### Group 12: Condition Number Variants (18 patterns)
Same structural patterns at κ=10^4, 10^8, 10^12: exponential_graded_k4/k8/k12, stemr_killer_k5/k10, huge_condition_k5/k10, demmel_G1_k4/k8/k12, demmel_S1pe_k4/k8, chkbd_4/chkbd_16, marques_graded_k4/k8

## Scaling Test

5 patterns tested at n=100, 200, 400. Doubling ratio = time(2n)/time(n).
- Patterns: exponential_graded, all_equal_nontrivial, stemr_killer, glued_wilkinson, chkbd
- Target: worst ratio ≤ 5.0 for O(n²)
- O(n²) → ratio ≈ 4.0, O(n³) → ratio ≈ 8.0

## Available Fortran Routines (via extern "C")

All available through `lib/libxmr.a` + Accelerate framework:
- **DBDSQR**: QR iteration bidiagonal SVD (O(n³))
- **DSTEMR**: MR³ tridiagonal eigensolver (O(n²)) — the key building block
- **DSTEXR**: Willems XMR improved MR³ (O(n²)) — in libxmr.a
- **DLAMCH**: Machine parameters (eps, safmin, etc.)
- **BLAS**: dnrm2, ddot, dscal, dcopy, daxpy, dgemv, etc.

DBDSVDX requires separate gfortran-compiled LAPACK (`lapack/build/lib/liblapack.a`).

## Fortran String Length Convention

gfortran on arm64 uses `size_t` (8 bytes) for hidden CHARACTER string lengths, not `int` (4 bytes). This matters when declaring `extern "C"` wrappers for Fortran routines with CHARACTER arguments. Accelerate framework uses `int`. Mixing them causes heap corruption.
