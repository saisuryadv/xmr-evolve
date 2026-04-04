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

**HARD GATE on scaling**: if `worst_ratio > 5.0`, **score is capped at 5** (compilation only).
Otherwise, full 10 pts scaling bonus is awarded.

Max from evaluate.cpp alone: **95 pts** (if worst_ratio ≤ 5.0).

### evaluator.py additional scoring (evaluator.py:188-225)

**HARD GATE** applies in evaluator.py too: if `worst_ratio > 5.0`, score is capped at 5.

If scaling passes AND stage 2 `pass_rate ≥ 0.6`, evaluator.py runs stage 3 and adds:
- `large_pass_rate × 10.0` — up to 10 pts
- Large scaling bonus — up to 5 pts (if large worst ratio ≤ 5.0)

Max from evaluator.py: **110 pts** (if scaling passes everywhere).

## Evaluator Cascade (evaluator.py)

| Stage | adv_size | STColl | Total Tests | Timeout | Gate |
|-------|----------|--------|-------------|---------|------|
| 1 | 50 | No | 270 | 30s | Always runs |
| 2 | 100 | Yes | 289 | 120s | Always runs (in `evaluate()`) |
| 3 | 200 | Yes | 379 | 600s | Stage 2 pass_rate ≥ 0.6 |

## Scaling Test

All patterns that pass at **every** test size are included. The worst time doubling ratio from the **largest size transition only** (e.g., n=200→400 for stage 3) is used. This is NOT limited to a fixed set of 5 patterns — all passing patterns contribute (evaluate.cpp:1025-1043).

Target: worst ratio ≤ 5.0 for O(n²). O(n²) → ratio ≈ 4.0, O(n³) → ratio ≈ 8.0.

## Extended Test Suite: Dense-to-Bidiag + Paper Tests (224 tests)

File: `python_fortran/test_dense_to_bidiag.py`

A supplementary test suite of **64 patterns / 224 total tests** designed to find failure cases in the MR3-GK bidiagonal SVD. Uses the same metrics and thresholds as the 379-test suite (res ≤ 7.0, ortU ≤ 5.0, ortV ≤ 5.0 in units of n·eps). Test sizes: {10, 100, 200, 400} (matching the 379 suite). Quick mode: {10, 100}.

**Current results: 219/224 passed (97.8%), 5 failures.**

### PART 1: Dense-to-Bidiagonal Tests (22 patterns × 4 sizes = 88 tests)

Construct dense matrices with prescribed singular value distributions, reduce to bidiagonal via Householder bidiagonalization (equivalent to LAPACK's DGEBRD), test `bidiag_svd` on the resulting bidiagonal. The key insight from Demmel-Kahan 1990: "reduction to bidiagonal form may produce completely inaccurate bidiagonal entries" — so the bidiagonal from dense reduction can have very different structure than hand-crafted bidiagonals, even for the same singular value distribution.

| # | Pattern | Reference | Characteristic Feature |
|---|---------|-----------|----------------------|
| 1 | `dlatms_mode1_eps` | LAPACK TESTING/EIG/dchkbd.f Type 5; Demmel et al. 2008 "S1pe" type | **n-1 clustered SVs at eps, 1 outlier at 1.** Tests cluster resolution when nearly all SVs are identical and tiny. After dense reduction, the bidiagonal has graded entries reflecting the eps-level cluster. |
| 2 | `dlatms_mode2_eps` | LAPACK dchkbd.f Type 5 variant; Demmel et al. 2008 "S2pe" type | **n-1 clustered SVs at 1, 1 outlier at eps.** Opposite of Mode 1 — cluster at top. Tests deflation of one tiny SV from a large cluster. |
| 3 | `dlatms_mode3_eps` | LAPACK dchkbd.f Type 4 (DLATMS MODE=3, COND=1/ULP); Demmel et al. 2008 "G1" type | **Geometric SV spacing from 1 to eps.** κ~10^16. Tests relative accuracy across 16 decades of dynamic range. The bidiagonal inherits grading from dense reduction. |
| 4 | `dlatms_mode4_eps` | LAPACK dchkbd.f Type 3 (DLATMS MODE=4, COND=1/ULP); Demmel et al. 2008 "U1" type | **Arithmetic SV spacing from 1 to eps.** Linear decay — SVs are evenly spaced, no large gaps. Tests accuracy when no natural cluster structure exists. |
| 5 | `dlatms_mode5_random` | LAPACK MATGEN/dlatms.f MODE=5 | **Log-uniform random SVs in (eps, 1).** Tests general robustness on random graded SVs without specific structure. Not used in LAPACK dchkbd (Modes 5/6 are absent from the 16 types). |
| 6 | `dense_three_clusters` | Inspired by Demmel et al. 2008 cluster classification | **Three well-separated SV clusters: {100}×k, {1}×k, {eps}×k.** Tests MR3's ability to independently resolve three clusters at different scales after dense reduction. |
| 7 | `dense_one_tight_cluster` | Inspired by Willems-Lang 2012 Theorem 4.5 orthogonality bound | **All SVs within n·eps of 1: σ(i) = 1 + i·eps.** Maximum clustering — every pair has relative gap ~eps. Tests whether representation tree can resolve O(n) eigenvalues that agree to ~15 digits. |
| 8 | `dense_half_rank` | Standard numerical linear algebra rank-deficiency test | **Rank n/2: σ(1:n/2) geometric, σ(n/2+1:n)=0.** Tests handling of exact zero SVs produced by dense reduction. Zero SVs should be detected by QR deflation sweep. |
| 9 | `dense_repeated_sv` | Inspired by LAPACK dchkbd Type 2 (identity = all SVs equal) | **Exact multiplicity: σ = [1]×n/2 + [eps]×n/2.** Two clusters with exact repeated values within each. Tests deflation and cluster resolution with exact degeneracy. |
| 10 | `dense_hilbert` | Classical; Hilbert 1894; standard SVD accuracy benchmark | **Hilbert matrix H(i,j)=1/(i+j+1).** SVs decay exponentially (κ grows as e^{3.5n}). The poster child for ill-conditioning — small SVs are determined to high relative accuracy by the matrix entries. |
| 11 | `dense_kahan` | Kahan 1966 "Numerical linear algebra"; Higham 2002 Ch. 8 | **Kahan matrix: R=diag(1,s,s²,...)·upper_tri(c), s=sin(θ), c=cos(θ).** Known to fool rank-revealing QR. After bidiag reduction, produces a graded bidiagonal where small SVs may lose relative accuracy. |
| 12 | `dense_condition_1e20` | Extreme variant of DLATMS MODE=3 | **Geometric SVs from 1 to 10^{-20}.** κ=10^20, spanning 20 decades. Tests whether dense reduction + bidiag SVD maintains relative accuracy across extreme dynamic range. |
| 13 | `dense_wilkinson_sv` | Willems-Lang 2012 Example 4.1 concept; Wilkinson 1965 | **SVs = eigenvalues of Wilkinson W_{2m+1}^+.** Eigenvalues come in tight pairs that differ by O(1/m!). **FAILS** at n=99 (ortU=29.9), n=199 (ortU=6.2), n=399 (ortU=27.6). The tight pairs from Wilkinson structure persist through dense reduction and stress MR3's cluster resolution. |
| 14 | `dense_wl_example41` | Willems-Lang 2012 Section 4, Example 4.1 (ETNA Vol. 39, pp. 1-21) | **Exact 20-SV construction: σ₁₃=0.9, σ₁₄=1-10⁻⁷, σ₁₅=1+10⁻⁷, σ₁₆=1.1, geometric extensions.** This is the matrix that demonstrates naive MR3 on T_GK fails — the initial outside shift clusters +σ and -σ for small σ. Tests whether dense reduction preserves this pathological structure. |
| 15 | `dense_near_overflow` | LAPACK dchkbd.f Types 6,11 (KMAGN=2) | **Geometric SVs from 10^{150} to 10^{140}.** Tests overflow handling in Householder reduction and subsequent bidiag SVD. Intermediate products in QR sweeps may overflow. |
| 16 | `dense_near_underflow` | LAPACK dchkbd.f Types 7,12 (KMAGN=3) | **Geometric SVs from 10^{-150} to 10^{-160}.** Tests gradual underflow in bidiag entries. Small entries may underflow to zero during QR deflation, creating false splits. |
| 17 | `dense_random` | LAPACK dchkbd.f Type 13 (KTYPE=9); standard stress test | **Random Uniform(-1,1) dense matrix.** No prescribed structure — tests general robustness. The SV distribution is unknown a priori (typically well-conditioned, κ~n). |
| 18 | `dense_vandermonde` | Classical; Gautschi 1975 "Norm estimates for inverses of Vandermonde matrices" | **Vandermonde on Chebyshev nodes.** Extremely ill-conditioned (κ grows exponentially). SVs cluster near zero. Tests accuracy on a matrix where bidiag reduction introduces large errors. |
| 19 | `dense_toeplitz` | Classical; Grenander-Szego 1958 "Toeplitz Forms and Applications" | **Symmetric Toeplitz, first row=[n,n-1,...,1].** Eigenvalues cluster in specific patterns determined by the generating function. Tests behavior on structured dense matrices with moderate clustering. |
| 20 | `dense_frank` | Frank 1958; Higham 2005 "The Test Matrix Toolbox" | **Frank matrix: F(i,j)=min(i,j)+1 (upper Hessenberg).** Has eigenvalues that come in nearly-reciprocal pairs. Tests bidiag SVD on structured matrices with paired singular values. |
| 21 | `dense_companion_exp` | LAPACK GitHub Issue #316 (2019); EigTool companion_demo.m | **Companion matrix of truncated Taylor series of exp(z): coeffs=[1,1,1/2!,...,1/n!].** Triggers DGESDD/DBDSDC bug for n≥26 (D&C produces wrong SVs with vectors). The bidiagonal from DGEBRD has entries spanning ~284 decades (matching B_bug316_gesdd.dat in STCollection). Tests extreme grading from a real-world matrix structure. |
| 22 | `dense_moler_like` | Inspired by Cleve Moler's 200×200 matrix (LAPACK Known Issues, DLASQ1 convergence); Marques improved dqds 2008 | **Bidiagonal with d≈±1 random signs, e quadratically decaying.** The original Moler matrix needed 8,073 QR iterations without flip (vs 970 with). Tests convergence robustness on oscillatory diagonal patterns. |

### PART 2: Missing Paper Test Cases (42 patterns, 136 tests)

#### 2A. Proper Glued Wilkinson (4 fixed-size tests)

**Reference:** Dhillon, Parlett & Vömel, "Glued Matrices and the MRRR Algorithm," SIAM J. Sci. Comput. 27:1 (2005), and LAPACK Working Note 163. Also Demmel et al. 2008: "Most difficult class for MR3."

**Characteristic:** p copies of Wilkinson W_{2m+1}^+ connected by off-diagonal γ. Each Wilkinson eigenvalue splits into p copies spaced O(γ²/gap) apart. The tightest pair has p copies spaced O(γ²·m!) — for m=10, these are equal to working accuracy, forcing deep representation trees.

| Pattern | Construction | n | Feature |
|---------|-------------|---|---------|
| `glued_wilk_3x21_sqrteps` | 3×W_21^+, γ=√eps | 63 | Moderate: 3 copies, eigenvalue triplication |
| `glued_wilk_5x21_sqrteps` | 5×W_21^+, γ=√eps | 105 | **FAILS** (ortU=13.1): 5 copies create O(5) near-degenerate eigenvalues per pair |
| `glued_wilk_3x21_1e14` | 3×W_21^+, γ=10⁻¹⁴ | 63 | Very weak coupling — eigenvalue splitting negligible |
| `glued_wilk_10x11_sqrteps` | 10×W_11^+, γ=√eps | 110 | Many small blocks — broad representation tree |

#### 2B. Glued Versions of Existing Patterns (8 patterns × 4 sizes)

**Reference:** Willems & Lang, "A Framework for the MR3 Algorithm: Theory and Implementation," SIAM J. Sci. Comput. 35:2 (2013), Section 8.1. The 116,874 Synth testset uses each base type in 3 variants: as-is, glued-small (~n·eps·‖T‖), glued-medium (~n·√eps·‖T‖). Our 379-test suite only has the "as-is" variants.

**Characteristic:** Gluing replicates the base spectrum 2-3 times, creating near-degenerate eigenvalue clusters. Small glue (γ~n·eps·‖B‖) makes copies indistinguishable to working precision. Medium glue (γ~n·√eps·‖B‖) creates marginally resolvable copies.

| Pattern | Base | Glue | Feature |
|---------|------|------|---------|
| `exponential_graded_glued_small` | 3 copies | n·eps·‖B‖ | Graded spectrum triplicated at eps-level coupling |
| `exponential_graded_glued_medium` | 2 copies | n·√eps·‖B‖ | Graded spectrum duplicated at √eps-level coupling |
| `stemr_killer_glued_small` | 3 copies | n·eps·‖B‖ | Already-hard geometric grading + triplication |
| `constant_glued_small` | 3 copies | n·eps·‖B‖ | Identical constant spectra barely separated |
| `all_equal_nontrivial_glued_small` | 3 copies | n·eps·‖B‖ | d=1, e=1 triplicated — massive cluster |
| `demmel_S1pe_glued_small` | 2 copies | n·eps·‖B‖ | Tight cluster at eps duplicated |
| `pd_T0_glued_medium` | 2 copies | n·√eps·‖B‖ | Uniform cluster duplicated at √eps coupling |
| `chkbd_glued_small` | 2 copies | n·eps·‖B‖ | Extreme grading (CHKBD) + duplication |

#### 2C. CHKBD LAPACK General (1 pattern × 4 sizes)

**Reference:** LAPACK TESTING/EIG/dchkbd.f Type 16; also STCollection bidiagonal Type 3. Identical to Demmel-Kahan 1990 numerical experiment format.

**Characteristic:** d,e = exp(γ·r), γ=-2·log(ULP), r~Uniform(-1,1). Entries span [ULP², ULP⁻²] ≈ 64 orders of magnitude. The widest dynamic range in any standard bidiagonal test. Different from existing `chkbd` (which uses d[i]=10^{-(2i+1)}, deterministic powers of 10).

| Pattern | Construction |
|---------|-------------|
| `chkbd_lapack_general` | Random exponential bidiagonal spanning ~64 decades |

#### 2D. Corrected Grosser-Lang (7 patterns × 4 sizes)

**Reference:** Großer & Lang, "An O(n²) Algorithm for the Bidiagonal SVD," Lin. Alg. Appl. 358 (2003), pp. 45-70. Test matrix generator DMATGEN.f, IDs 200-244. Verified by code comparison that our existing `gl_abcon1-3` and `gl_gro0-3` do NOT match the paper's definitions.

**Characteristic:** GRO0-3 test element growth in coupling: outlier diagonal entries at 1 with bulk at eps create extreme condition. ABCON variants are 1D heat equation tridiagonals with boundary modifications.

| Pattern | Paper ID | Construction | Feature |
|---------|----------|-------------|---------|
| `gl_gro0_proper` | 240 | d[0]=1, d[1:]=eps, e[:]=eps | **1 outlier SV near 1, rest at ~10⁻²⁷.** Tests element growth when shifting near the large SV |
| `gl_gro1_proper` | 241 | d[0:2]=1, rest=eps, e[:]=eps | **2 outliers.** Two large SVs create a small cluster at top + massive cluster at bottom |
| `gl_gro2_proper` | 242 | d[0:4]=1, rest=eps, e[:]=eps | **4 outliers.** Larger cluster at top |
| `gl_gro3_proper` | 243 | d[0:4]=1, e[0:4]=1, rest=eps | **Outlier block with strong coupling.** First 4 entries are fully coupled (e=1), rest at eps |
| `gl_abcon1_proper` | 201 | d[:]=2, d[0]=1, e[:]=1 | Heat equation with modified left boundary |
| `gl_abcon2_proper` | 202 | d[:]=2, d[0]=1, d[-1]=3, e[:]=1 | Modified both boundaries |
| `gl_abcon3_proper` | 203 | d[:]=2, d[0]=d[-1]=3, e[:]=1 | Symmetric boundary modification |

#### 2E. Marques STEXR Failure (3 fixed-size tests)

**Reference:** Marques, Demmel & Vasconcelos, "Bidiagonal SVD Computation via an Associated Tridiagonal Eigenproblem," ACM TOMS 46:2 (2020), Appendix B. Table B.1 gives exact entries.

**Characteristic:** Tridiagonal with geometrically spaced eigenvalues λᵢ=c^((i-1)/(n-1)), c=1/√eps. Creates eigenvalues spanning ~8 decades (1.5×10⁻⁸ to 1.0). STEXR (Willems-Lang MR3) produces orthogonality 3.95×10⁴·n·eps — first four eigenvectors linearly dependent. STEMR and STEVX pass the same matrix.

| Pattern | n | Feature |
|---------|---|---------|
| `marques_stexr_10` | 10 | Original failure size from the paper |
| `marques_stexr_20` | 20 | **FAILS** (ortU=8.7): scaled version exposes our implementation |
| `marques_stexr_50` | 50 | Larger version, passes (clustering structure changes) |

#### 2F. Tridiagonal-to-Bidiagonal via Cholesky (5 patterns, mixed sizes)

**Characteristic:** Construct symmetric tridiagonal T, Cholesky-factor T-ν·I = LL^T (ν from Gershgorin lower bound), set B = L^T. The bidiagonal inherits the eigenvalue structure of T as squared singular values: σᵢ² = λᵢ(T) - ν.

| Pattern | Reference | Construction | Feature |
|---------|-----------|-------------|---------|
| `wilkinson_cholesky_m10` | Großer & Lang 2005 "On Symmetric Eigenproblems..." | W_21^+ → Cholesky → bidiag | **Tight SV pairs from Wilkinson eigenvalue pairs.** σ pairs differ by O(1/10!) |
| `wilkinson_cholesky_m50` | Dhillon, Parlett & Vömel 2005 (LAWN 163) | W_101^+ → Cholesky → bidiag | **FP underflow threshold.** For m≥50, eigenvectors of W decay so rapidly that FP computation underflows. MRRR returns bisectors instead of true eigenvectors |
| `legendre_bidiag` | STCollection Type 5; Abramowitz & Stegun 1964 Ch. 22 | Jacobi(Legendre): d=0, e[i]=i/√(4i²-1) → Cholesky | **Eigenvalues = zeros of Legendre polynomials.** Cluster near ±1 for large n. Well-separated interior |
| `laguerre_bidiag` | STCollection Type 6; Abramowitz & Stegun 1964 Ch. 22 | Jacobi(Laguerre): d[i]=2i+1, e[i]=i+1 → Cholesky | **Eigenvalues spread on (0,∞).** Graded structure with increasing spacing |
| `hermite_bidiag` | STCollection Type 7; Abramowitz & Stegun 1964 Ch. 22 | Jacobi(Hermite): d=0, e[i]=√(i+1) → Cholesky | **Eigenvalues = zeros of Hermite polynomials.** Symmetric about 0, growing spacing |

#### 2G. Missing Variants (4 patterns × 4 sizes)

| Pattern | Reference | Construction | Feature |
|---------|-----------|-------------|---------|
| `demmel_S1pe_k16` | Demmel et al. 2008, S1pe variant | d[0]=1, d[1:]=10⁻¹⁶, e=√(d·d)·0.1 | **Higher condition (κ=10¹⁶) than existing k4/k8.** Pushes cluster at 10⁻¹⁶ further from outlier at 1 |
| `wl_element_growth` | Willems-Lang 2012 Example 4.8 (ETNA Vol. 39, p. 17) | d=[1,α,α²,...], e=[1,α,...], α=eps | **Element growth vulnerability.** Shifting T_GK by -α produces D(2)~1/eps~10¹⁶. Tests whether NCD check prevents catastrophic element growth at deeper tree levels |
| `random_zero_e_10pct` | Marques et al. 2020 Section 9 (zero-entry injection methodology) | Random d, 10% of e[i]=0 | **Splitting test: ~10 independent sub-problems.** Tests splitting detection and reassembly of partial SVDs |
| `random_zero_e_50pct` | Marques et al. 2020 Section 9 | Random d, 50% of e[i]=0 | **Heavy splitting: many 1×1 and 2×2 sub-problems.** Most SVs are just |d[i]|, tests degenerate case handling |

#### 2H. From Dhillon-Parlett-Vömel 2005 / LAWN 163/166 / Dhillon Thesis (4 patterns)

| Pattern | Reference | Construction | Feature |
|---------|-----------|-------------|---------|
| `bidiag_coupling_failure` | Großer & Lang 2005 eq (3.2); LAPACK Working Note 166 (Willems, Lang & Vömel 2005) | d[0]=1, d[1:]=200·eps, e[:]=200·eps | **Coupling failure stress test.** One SV near 1, rest clustered at ~10⁻²⁷. Without couplings, TGK eigenvector extraction gives orthogonality > 10¹⁰. Tests whether our MR3-GK Algorithm 4.1 avoids this via NCD-aware shifting |
| `dhillon_uniform_sqrteps_apart` | Dhillon 1997 PhD Thesis (UC Berkeley), Section 6.4.1 Type 2 | λ_1=eps, λ_i=1+i·√eps, λ_n=2 → dense → bidiag | **Eigenvalue gaps = √eps — right at MRRR clustering boundary.** gaptol~10⁻³ means these gaps are above threshold, but barely. Tests the transition between "singleton" and "cluster" classification |
| `dhillon_121_toeplitz` | Dhillon 1997 Thesis Type 12; classical (1,2,1) Toeplitz tridiagonal | tridiag(1,2,1) → Cholesky → bidiag | **Analytically known eigenvalues: 4·sin²(kπ/(2(n+1))).** Provides ground-truth verification. Eigenvalues cluster near 0 and 4 for large n |
| `glued_wilk201_5x_sqrteps` | Dhillon, Parlett & Vömel 2005 (SIAM J. Sci. Comput. 27:1); LAWN 163 | 5×W_201^+ glued γ=√eps, n=1005 | **THE actual MRRR failure matrix from the literature.** 5 copies of W with m=100. Tightest pairs have 5 eigenvalues equal to working accuracy. Representation tree grows to depth 21+. FP vector computation underflows for m=100. Passes our implementation (ortU=0.14), suggesting our NCD-aware approach handles it |

#### 2I. STCollection Eigenvalue Distribution Types 7/8/9 (3 patterns × 4 sizes)

**Reference:** Marques, Demmel, Parlett & Vömel, "A Testing Infrastructure for Symmetric Tridiagonal Eigensolvers," ACM TOMS 35:1 (2008), Algorithm 880, Table II.

| Pattern | Construction | Feature |
|---------|-------------|---------|
| `stcoll_evdist_type7` | W(i)=ULP·i, W(N)=1 → dense → bidiag | **Linearly growing from ULP.** n-1 tiny eigenvalues scaling as i·eps, one outlier at 1. Tests accuracy of tiny SVs that grow linearly |
| `stcoll_evdist_type8` | W(1)=ULP, W(i)=1+√ULP·i, W(N)=2 → dense → bidiag | **One tiny outlier + tight cluster near 1.** n-2 eigenvalues spaced √eps apart around 1, one at eps, one at 2. Tests simultaneous handling of isolated tiny SV + tight cluster |
| `stcoll_evdist_type9` | W(1)=1, W(i)=W(i-1)+100·ULP → dense → bidiag | **Extremely tight arithmetic cluster.** All n eigenvalues within 100·n·eps of each other. Gaps = 100·eps — just above machine precision. Maximum stress on cluster resolution |

#### 2J. Skew-Wilkinson with Grading (3 fixed-size tests)

**Reference:** STCollection (Marques et al.), T_SkewW21 variants. Wilkinson W_21^+ with asymmetric off-diagonal scaling.

**Characteristic:** Standard Wilkinson has e[:]=1. Skew variants scale e[i] by different factors above/below the center, breaking the diagonal symmetry of the tridiagonal. This changes which eigenvalue pairs are tight and where eigenvector weight concentrates.

| Pattern | Construction | Feature |
|---------|-------------|---------|
| `skew_wilk21_g1` | W_21^+ standard → Cholesky | Baseline (no grading). n=21 fixed |
| `skew_wilk21_g1e3` | W_21^+ with e[i<m]×10³, e[i≥m]×10⁻³ | **Asymmetric 10³ grading.** Upper half strongly coupled, lower half weakly coupled |
| `skew_wilk21_g1e6` | W_21^+ with e[i<m]×10⁶, e[i≥m]×10⁻⁶ | **Extreme asymmetric 10⁶ grading.** 12-decade contrast between upper/lower coupling |

#### 2K. W21 with Varying Off-diagonal Gamma (3 fixed-size tests)

**Reference:** STCollection, T_W21_g_* variants. Single W_21^+ with all off-diagonals uniformly scaled.

**Characteristic:** Scaling off-diagonals by γ changes the eigenvalue gaps: for small γ, eigenvalue pairs become nearly degenerate (gap ~ γ²). For γ=10⁻¹⁴, pairs agree to ~28 digits — far beyond machine precision.

| Pattern | Construction | Feature |
|---------|-------------|---------|
| `wilk21_gamma_1e4` | W_21^+ with e×10⁻⁴ | Moderate decoupling. Pairs differ by ~10⁻⁸ |
| `wilk21_gamma_1e8` | W_21^+ with e×10⁻⁸ | Strong decoupling. Pairs differ by ~10⁻¹⁶ (≈eps) |
| `wilk21_gamma_1e14` | W_21^+ with e×10⁻¹⁴ | **Extreme decoupling.** Pairs agree to ~28 digits. Nearly diagonal matrix. Tests splitting detection |

#### 2L. Willems 2x2 Block LDL* Edge Case (1 pattern × 4 sizes)

**Reference:** Willems 2010 PhD Thesis (Bergische Universität Wuppertal); Parlett's slides on Willems Thesis, pp. 186-228.

**Characteristic:** T = GK(B) - α·I with α=O(eps) produces LDL* with entries requiring 2×2 blocks in the diagonal (Ω set nonempty in Willems' block factorization framework). Standard scalar LDL* factorization fails to be an RRR for these representations. Tests whether our implementation's block factorization handling works correctly on the specific edge case from the thesis.

| Pattern | Construction | Feature |
|---------|-------------|---------|
| `willems_block_ldl_edge` | Repeated 4×4 blocks: d=α, e=[1,1,α], glued by √eps | Forces 2×2 blocks in LDL* factorization. The GK structure with near-zero diagonal creates large element growth (D(2)~1/eps) after one shift |

### Failures Found (5/224)

| Test | n | res | ortU | ortV | Root Cause |
|------|---|-----|------|------|------------|
| `dense_wilkinson_sv` | 99 | 4.01 | **29.86** | 29.86 | Tight eigenvalue pairs from W_{2m+1}^+ as SVs. Pairs differ by O(1/m!) — MR3 cannot resolve these clusters |
| `dense_wilkinson_sv` | 399 | 1.23 | **27.60** | 27.60 | Same mechanism at larger n |
| `dense_wilkinson_sv` | 199 | 0.77 | **6.15** | 6.15 | Same, moderate — close to threshold |
| `glued_wilk_5x21_sqrteps` | 105 | 1.21 | **13.11** | 13.09 | 5 copies of W_21^+ eigenvalues replicated — O(5) near-degenerate values per pair that MR3 cannot separate |
| `marques_stexr_20` | 20 | 1.49 | **8.74** | 8.71 | Geometric eigenvalue distribution λᵢ=c^((i-1)/19) at n=20 creates clustering that our MR3-GK doesn't fully resolve |

### Running the Extended Suite

```bash
cd python_fortran
python3 test_dense_to_bidiag.py          # Full suite: 64 patterns, 224 tests
python3 test_dense_to_bidiag.py --quick  # Quick mode: sizes {10, 100} only
python3 test_dense_to_bidiag.py --part1  # Dense-to-bidiag tests only
python3 test_dense_to_bidiag.py --part2  # Paper tests only
```

### Source Provenance

| Source | Paper/Code Reference | # Patterns |
|--------|---------------------|-----------|
| LAPACK dchkbd.f / ddrvbd.f | LAPACK TESTING/EIG/dchkbd.f Types 3-16; TESTING/EIG/ddrvbd.f Types 3-5 | 8 |
| LAPACK bug #316 | GitHub Reference-LAPACK/lapack Issue #316 (2019); EigTool companion_demo.m | 1 |
| LAPACK DBDSQR convergence | LAPACK Known Issues (netlib.org); Cleve Moler's 200×200 matrix | 1 |
| Willems-Lang 2012 | "The MR3-GK Algorithm for the Bidiagonal SVD," ETNA Vol. 39, pp. 1-21 (2012). Examples 4.1, 4.8 | 2 |
| Willems-Lang 2013 | "A Framework for the MR3 Algorithm," SIAM J. Sci. Comput. 35:2 (2013). Section 8.1 Synth methodology | 8 |
| Willems Thesis 2010 | PhD Thesis, Bergische Universität Wuppertal. 2×2 block LDL* edge case | 1 |
| Demmel et al. 2008 | "Performance and Accuracy of LAPACK's Symmetric Tridiagonal Eigensolvers," SIAM J. Sci. Comput. 30:3 (2008) | 2 |
| Dhillon-Parlett-Vömel 2005 | "Glued Matrices and the MRRR Algorithm," SIAM J. Sci. Comput. 27:1 (2005); LAPACK Working Notes 163, 166 | 4 |
| Dhillon Thesis 1997 | PhD Thesis, UC Berkeley. Section 6.4.1 Types 2, 12 | 2 |
| Großer-Lang 2001 | "An O(n²) Algorithm for the Bidiagonal SVD," Lin. Alg. Appl. 358 (2003). DMATGEN.f IDs 200-244 | 7 |
| Marques et al. 2020 | "Bidiagonal SVD via Associated Tridiagonal Eigenproblem," ACM TOMS 46:2 (2020). Appendix B, Section 9 | 5 |
| STCollection (Marques et al.) | "A Testing Infrastructure for Symmetric Tridiagonal Eigensolvers," ACM TOMS 35:1 (2008). Algorithm 880 | 14 |
| Classical matrices | Hilbert (1894), Kahan (1966), Vandermonde/Gautschi (1975), Frank (1958), Toeplitz/Grenander-Szego (1958) | 5 |
| Custom | Inspired by cluster/rank/gap analysis from the above papers | 4 |

---

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
