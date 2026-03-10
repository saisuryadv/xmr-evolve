# Demmel, Marques, Parlett, Vomel 2008 -- Complete Extraction

**Paper**: "Performance and Accuracy of LAPACK's Symmetric Tridiagonal Eigensolvers"
SIAM J. Sci. Comput., Vol. 30, No. 3, pp. 1508-1526 (2008). DOI: 10.1137/070688778

**Authors**: James W. Demmel, Osni A. Marques, Beresford N. Parlett, Christof Vomel

**Algorithms compared**:
1. QR iteration -- STEV (LAPACK driver)
2. Bisection and inverse iteration -- STEVX (BI)
3. Divide and conquer -- STEVD (DC)
4. Multiple relatively robust representations -- STEVR (MR)

---

## 1. Synthetic Test Matrices -- Complete Specification

All synthetic matrices are generated using LAPACK's `DLATMS` routine, which creates
a symmetric tridiagonal matrix with a specified eigenvalue distribution. The eigenvalue
distribution is controlled by `MODE` and `COND` parameters.

### DLATMS MODE Definitions (eigenvalue distribution)

Given condition number kappa and dimension n:

| MODE | Eigenvalue distribution D(i), i=1..n |
|------|---------------------------------------|
| 0    | D is provided as input |
| 1    | D(1) = 1, D(2:n) = 1/kappa -- one large, rest small |
| 2    | D(1:n-1) = 1, D(n) = 1/kappa -- rest large, one small |
| 3    | D(i) = kappa^(-(i-1)/(n-1)) -- geometric spacing |
| 4    | D(i) = 1 - (i-1)/(n-1) * (1 - 1/kappa) -- uniform (arithmetic) spacing |
| 5    | D = random in (1/kappa, 1), logarithms uniformly distributed |
| 6    | D = random from same distribution as rest of matrix |
| -k   | Same as MODE=k but eigenvalues in reversed order |

For symmetric matrices (`SYM='S'`), eigenvalues are randomly negated (sign multiplied
by +1 or -1) when MODE is not 0, 6, or -6. The suffix "p" means positive (no negation),
"n" means negated (all flipped).

### kappa (Condition Number) Choices

- Default: kappa = 1/eps (approximately 4.5e15 in double precision)
- For strongly clustered (S): kappa = 1/eps ("e" suffix) or kappa = 1/sqrt(eps) ("s" suffix)

### S1, S2: Strongly Clustered Eigenvalues

**S1**: MODE = 1, kappa as specified
- D(1) = 1, D(2) = D(3) = ... = D(n) = 1/kappa
- Spectrum: n-1 eigenvalues clustered at 1/kappa, 1 eigenvalue at 1
- Fraction deflated in DC: 0.88 to 1.00 (nearly all deflated)

**S2**: MODE = 2, kappa as specified
- D(1) = D(2) = ... = D(n-1) = 1, D(n) = 1/kappa
- Spectrum: n-1 eigenvalues clustered at 1, 1 eigenvalue at 1/kappa
- Fraction deflated in DC: 0.88 to 0.98

Suffixes:
- **S1pe** / **S1ne**: MODE = 1 / MODE = -1, kappa = 1/eps, positive/negated eigenvalues
- **S1ps** / **S1ns**: MODE = 1 / MODE = -1, kappa = 1/sqrt(eps)
- **S2pe** / **S2ne**: MODE = 2 / MODE = -2, kappa = 1/eps
- **S2ps** / **S2ns**: MODE = 2 / MODE = -2, kappa = 1/sqrt(eps)

### W1, W2, W3: Weakly Clustered Eigenvalues

All use kappa = 1/eps by default.

**W1**: Eigenvalues at 1 and 1/kappa * [1 : n-1]
- D = {1/kappa, 2/kappa, ..., (n-1)/kappa, 1}
- Spectrum: arithmetically spaced from 1/kappa to (n-1)/kappa, plus one eigenvalue at 1
- Fraction deflated in DC: 0.34, 0.81

**W2**: Eigenvalues at 1/kappa, 2, and 1 + 1/sqrt(kappa) * [1 : n-2]
- D = {1/kappa, 1 + 1/sqrt(kappa), 1 + 2/sqrt(kappa), ..., 1 + (n-2)/sqrt(kappa), 2}
- Spectrum: one outlier at 1/kappa, (n-2) eigenvalues weakly clustered around 1, plus 2
- Fraction deflated in DC: 0.01, 0.05

**W3**: Eigenvalues at 1 + 100/kappa * [1 : n]
- D = {1 + 100/kappa, 1 + 200/kappa, ..., 1 + 100n/kappa}
- Spectrum: very tightly clustered around 1 (spacing ~100*eps)
- Fraction deflated in DC: 0.06, 0.09

### G1, G2: Geometric Distributions

Both use kappa = 1/eps.

**G1**: MODE = 3, exactly geometric spacing
- D(i) = kappa^(-(i-1)/(n-1))
- Spectrum: logarithmically spaced from 1 to 1/kappa
- Fraction deflated in DC: 0.16, 0.36

**G2**: MODE = -3, randomly geometric spacing
- Same geometric spacing but in reversed order (plus random negation for symmetric)
- Fraction deflated in DC: 0.16, 0.19

### U1, U2: Uniform (Arithmetic) Distributions

Both use kappa = 1/eps.

**U1**: MODE = 4, exactly uniform spacing
- D(i) = 1 - (i-1)/(n-1) * (1 - 1/kappa)
- Spectrum: arithmetically spaced from 1 to 1/kappa
- Fraction deflated in DC: 0 (no deflation)

**U2**: MODE = -4, randomly uniform spacing
- Same arithmetic spacing but in reversed order (plus random negation for symmetric)
- Fraction deflated in DC: 0 (no deflation)

### Wi: Wilkinson W_{m+1}^+ Matrices

The classic Wilkinson matrix of order n = 2m+1:

```
W_{2m+1}^+ = tridiag(1, [m, m-1, ..., 1, 0, 1, ..., m-1, m], 1)
```

- Diagonal: d(i) = |m - (i-1)| for i=1,...,2m+1 (i.e., m, m-1, ..., 1, 0, 1, ..., m-1, m)
- All off-diagonal elements = 1
- Eigenvalues come in close pairs near the top of the spectrum
- Eigenvalues of the leading/trailing m x m submatrix are Ritz values
- Fraction deflated in DC: 0.35, 0.84
- Notoriously difficult for MR3 due to:
  - Tight eigenvalue pairs at top of spectrum
  - Eigenvectors with "essential" nonzero part at top and bottom, decaying toward middle
  - Underflow in FP vector computation for large m
  - MRRR may return bisectors instead of true eigenvectors

### GW: Glued Wilkinson Matrices

Multiple copies of W_{m+1}^+ "glued" together with a coupling parameter gamma:

```
T(gamma) = block_diag(T_1, ..., T_p) + gamma * sum_{i=1}^{p-1} v_i * v_i^T
```

where:
- T_1, ..., T_p are copies of W_{2m+1}^+ (all identical)
- v_i is a vector with two nonzero entries of size 1: at the last row of T_i and first row of T_{i+1}
- In tridiagonal form: this means the off-diagonal element connecting block i to block i+1 is gamma
- Typical parameter: gamma = sqrt(eps)

Construction in practice (for p copies of W_{2m+1}^+, glued by gamma):
```
n = p * (2m+1)
diagonal = [m, m-1, ..., 0, ..., m-1, m, m, m-1, ..., 0, ..., m-1, m, ...]  (repeat p times)
off-diag = [1, 1, ..., 1, gamma, 1, 1, ..., 1, gamma, ...]  (1s within blocks, gamma between blocks)
```

- Fraction deflated in DC: 0.59, 0.78
- **Most difficult class for MR3**: eigenvalue clusters come in groups of small size but extreme tightness
- The representation tree becomes very broad
- The glue gamma is small enough that eigenvalue clusters from different copies are barely split
- The effect of glue on eigenvalue separation is tied to top/bottom entries of eigenvectors
- When gamma's effect is negligible, MR3 faces effectively multiplicity-p eigenvalues

---

## 2. Practical Matrix Sources

### Source 1: NWChem Computational Chemistry (Fann)

- **Software**: NWChem package (refs [3, 27] in paper)
- **Application**: Computational chemistry, analysis of molecules
- **Property**: Matrices have clustered eigenvalues requiring a large number of reorthogonalizations in BI
- **Historical significance**: These matrices motivated the development of MR3 (Dhillon, Fann and Parlett, refs [12, 14])
- **Reduction method**: For smaller matrices, tridiagonal form obtained via LAPACK's `sytrd` routine. For larger matrices, a simple Lanczos algorithm without reorthogonalization was used, which produces copies of eigenvalues as clusters in finite precision.

### Source 2: Sparse Matrix Collections (BCSSTRUC1 set)

- **Collections**: BCSSTRUC1 set (refs [20, 21, 22]), plus matrices from:
  - Alemdar
  - National Aeronautics and Space Administration (NASA)
  - Cannizzo sets (ref [9], Tim Davis's University of Florida sparse matrix collection)
- **Applications**: Power system networks, shallow wave equations, finite-element problems
- **Spectral property**: Typically a "continuously" varying part of the spectrum plus several isolated large clusters of eigenvalues of varying tightness

### Reduction to Tridiagonal Form

- **Small matrices**: LAPACK's `sytrd` (Householder tridiagonal reduction)
- **Large matrices**: Simple Lanczos without reorthogonalization (produces artificial eigenvalue clusters)

---

## 3. Bugs and Failure Modes by Algorithm

### QR Iteration (STEV)
- **No accuracy failures documented**
- Most accurate algorithm (with DC)
- Observed accuracy: O(sqrt(n) * eps)
- **Performance weakness**: O(n^3) complexity; BLAS 1-like operations on eigenvector matrix cause MFlop rate drop when eigenvector matrix exceeds L3 cache

### Bisection and Inverse Iteration (STEVX / BI)
- **Known theoretical failure**: BI may completely fail to guarantee orthogonal eigenvectors (ref [13], Dhillon 1998: "Current inverse iteration software can fail")
- **Observed**: Did not occur catastrophically on the test matrices in this study, but known to be rare
- **Accuracy**: O(n*eps) -- less accurate than QR and DC
- **Orthogonality issues for subset computations**: A few subset computations showed issues (Table 6.1 shows worst-case orthogonality loss up to 440 * n*eps on SUN platform)
- **Performance issue**: When eigenvalues are tightly clustered, Gram-Schmidt reorthogonalization increases cost to O(n*k^2), potentially O(n^3) in worst case
- **Specific failures**: Worst-case residuals of 100-331 * n*eps across platforms; worst-case orthogonality loss of 53-440 * n*eps

### Divide and Conquer (STEVD / DC)
- **No accuracy failures documented**
- Most accurate algorithm (with QR)
- Observed accuracy: O(sqrt(n) * eps)
- **Performance issue**: O(n^3) worst case; strongly matrix-dependent. When deflation fraction is low (U1, U2, W2, W3), DC does close to O(n^3) work. When deflation is high (S1, S2, Wi, GW), DC is much faster than O(n^2)
- **Memory**: Uses the most memory: 2n^2 + O(n) reals, plus 3 + 5n integers

### Multiple Relatively Robust Representations (STEVR / MR)
- **LAPACK 3.0 version had major failures**:
  - Returned very large errors (residuals and orthogonality up to 10^14 * n*eps)
  - **Completely failed** to return any answer on 22 test matrices, including 9 practical matrices
  - These failures are documented in Figure 6.2 of the paper
- **LAPACK 3.1 version improvements**:
  - Fixed a bug where no off-diagonal entries satisfy the splitting criterion but eigenvalues can still be numerically indistinguishable down to the underflow threshold (ref [17])
  - Tightened the internal accuracy threshold on relative gaps from 1/n to a fixed 10^{-3} in double precision, making MR more accurate on larger matrices
- **Remaining accuracy limitations in 3.1**:
  - Accuracy is O(n*eps), not O(sqrt(n)*eps) like QR/DC
  - Worst-case residuals: 14-80 * n*eps across platforms
  - Worst-case orthogonality loss: 70-190 * n*eps across platforms
- **Specific difficult matrices for MR**:
  - **Wilkinson matrices**: Eigenvalues in tight pairs; FP vector computation corrupted by underflow; MRRR returns bisectors instead of true eigenvectors for large m
  - **Glued Wilkinson matrices (GW)**: Most difficult class. Eigenvalue clusters of extreme tightness; representation tree becomes very broad; glue parameter makes eigenvalue separation negligible for some eigenvalues. MR overhead for tree generation is considerable.
  - **Strongly clustered with kappa = 1/eps (S1ne, S2ne, etc.)**: Eigenvalues agree to ~15 digits, requiring deep representation trees
- **Performance characteristic**: MR always does fewest flops but at lowest MFlop rate (due to high fraction of divides vs. multiply/add)
- **MR3 Failure Mechanisms** (from WN 163, Dhillon-Parlett-Vomel 2005):
  1. **Assumption 4.1 failure** (finding large relative gaps): By suitable shifting, at least one eigenvalue in a tight cluster should get a relative gap exceeding any threshold. Fails for glued matrices because glue has negligible effect on eigenvalue separation.
  2. **Assumption 4.2 failure** (finding an RRR): It should always be possible to find an RRR close enough to a cluster endpoint to yield a singleton. Fails because even-indexed eigenvalues of W_{2m+1}^+ are Ritz values, and LDL^T factorizations break down near Ritz values.
  3. **Assumption 4.4 failure** (computing the FP vector): The eigenvalue approximation should have high relative accuracy. Fails due to underflow: eigenvectors of Wilkinson matrices decay rapidly from their "essential" nonzero parts at top/bottom toward the middle. For large m, the FP vector computation (equation 2.12) underflows, producing bisectors instead of true eigenvectors.
  4. **Random perturbation remedy**: Small componentwise relative random perturbation of the root representation. Preserves O(n^2) complexity and introduces artificial roundoff effects that break the special structure causing failures.

---

## 4. Accuracy Results

### Table 6.1: Worst-Case Errors (multiples of n*eps) -- All Matrices, All Platforms

| Architecture | Residual QR | Residual BI | Residual DC | Residual MR | Orth QR | Orth BI | Orth DC | Orth MR |
|---|---|---|---|---|---|---|---|---|
| SP3 (Power 3) | 0.23 | 100 | 0.13 | 22 | 0.34 | 140 | 0.25 | 70 |
| SP5 (Power 5) | 0.30 | 59 | 0.13 | 18 | 0.39 | 70 | 0.25 | 163 |
| SUN | 0.20 | 331 | 0.11 | 14 | 0.45 | 440 | 0.19 | 92 |
| SGI | 0.30 | 210 | 0.13 | 14 | 0.46 | 280 | 0.19 | 160 |
| ITN2 (Itanium 2) | 0.30 | 240 | 0.13 | 29 | 0.45 | 320 | 0.19 | 190 |
| P4 (Pentium 4 Xeon) | 0.30 | 39 | 0.13 | 33 | 0.38 | 53 | 0.19 | 140 |
| X1 (Cray X1) | 0.30 | 34 | 0.11 | 80 | 0.46 | 45 | 0.19 | 160 |
| OPT (Opteron) | 0.30 | 100 | 0.11 | 14 | 0.46 | 130 | 0.19 | 160 |

### Key Accuracy Observations

- **QR and DC**: Most accurate. Errors decrease as n increases (trend line slopes negative in log-log). Observed accuracy O(sqrt(n) * eps).
- **BI**: Least accurate. Worst-case residuals up to 331 * n*eps, worst orthogonality loss up to 440 * n*eps. Known to theoretically fail completely (ref [13]), though rare.
- **MR (LAPACK 3.1)**: O(n*eps) accuracy. Never exceeded 190 * n*eps orthogonality loss or 80 * n*eps residual on any platform.
- **MR (LAPACK 3.0)**: Had errors up to 10^14 * n*eps and complete failures on 22 matrices (Figure 6.2).

### Accuracy Trend Line Slopes (Opteron, Figure 6.1)

For practical matrices (residuals vs n on log-log):
- QR: slope -0.3
- BI: slope -0.5
- DC: slope -0.5
- MR: slope 0.4

For practical matrices (orthogonality loss vs n on log-log):
- QR: slope -0.3
- BI: slope -0.7
- DC: slope -0.4
- MR: slope 0.5

For testing matrices:
- QR: residual slope -0.4, orth slope -0.3
- BI: residual slope -0.5, orth slope -0.7
- DC: residual slope -0.7, orth slope -0.7
- MR: residual slope -0.3, orth slope -0.2

**Interpretation**: Negative slopes mean errors DECREASE relative to n*eps as n grows (i.e., actual error grows slower than n*eps). QR and DC have true accuracy of approximately O(sqrt(n)*eps). Positive slopes for MR mean error approaches the n*eps bound more tightly as n grows.

---

## 5. Performance Scaling Results

### Table 4.1: Slope of Time Trends (t = c * n^alpha) -- Practical Matrices, All Architectures

| Architecture | QR alpha | BI alpha | DC alpha | MR alpha |
|---|---|---|---|---|
| SP3 | 3.3 | 2.1 | 2.8 | 2.5 |
| SP5 | 3.0 | 2.6 | 2.5 | 2.3 |
| SUN | 3.8 | 2.0 | 2.6 | 2.4 |
| SGI | 3.5 | 3.2 | 2.7 | 2.3 |
| ITN2 | 3.0 | 2.4 | 2.5 | 2.3 |
| P4 | 3.0 | 2.6 | 2.5 | 2.4 |
| X1 | 2.4 | 2.0 | 1.9 | 2.2 |
| OPT | 2.9 | 2.9 | 2.5 | 2.2 |

### Opteron Detailed Slopes (Figure 4.1, 4.2)

By time: QR = O(n^{2.9}), BI = O(n^{2.9}), DC = O(n^{2.5}), MR = O(n^{2.2})
By flops: QR = O(n^{3.0}), BI = O(n^{2.9}), DC = O(n^{2.8}), MR = O(n^{2.2}) (flop count)

### Per-Matrix-Class Slopes (Time/n^2 trend on Opteron)

**DC slopes** (Figure 4.5, time/n^2):
- S1pe, S1ps, S1ne, S1ns, S2pe, S2ps, S2ne, S2ns: nearly O(n) => total ~O(n) [massive deflation]
- W1p, W1n: ~O(n^1) in time/n^2 (moderate deflation, 34-81%)
- GW, Wi: intermediate (59-84% deflation)
- G1p, G1n, G2p, G2n: moderate
- U1p, U2p, W2p, W2n, W3p, Prac: nearly O(n^3) [minimal deflation]

**MR slopes** (Figure 4.7, time/n^2):
- S1ne, S1pe, S1ps: decreasing => ~O(n^2) or less [singletons, trivial for MR]
- S1ns: flat line => O(n^2) [moderate clustering]
- Wi, GW: highest time/n^2 [most difficult for MR]
- Most matrix classes: roughly flat or slowly increasing

### Time Ratios: Algorithm vs MR (practical matrices, n >= 500)

| Metric | Min | Median | Max |
|---|---|---|---|
| QR/MR | 0.38 | varied | 710 |
| BI/MR | 0.74 | varied | 300 |
| DC/MR | 0.018 | varied | 12 |

### Flop Count Ratios vs MR (Opteron)
- MR always does fewest flops (up to 390x fewer than QR, 380x fewer than BI)
- MR does up to 66x fewer flops than DC, never more than ~2x more than DC
- Median: MR does 8.3x fewer flops than DC for large practical matrices
- **But MR has lowest MFlop rate** due to high fraction of divide operations

### GFlop Rates (Opteron, for large practical matrices)

| Algorithm | GFlop Rate (range) |
|---|---|
| QR | 1.0 - 2.0 |
| BI | 0.6 - 1.3 |
| DC | 1.5 - 4.1 |
| MR | 0.5 - 0.8 |

DC achieves highest GFlop rate because it uses BLAS 3 (GEMM). MR is lowest because of many scalar divides.

---

## 6. Condition Numbers Tested

- **Default**: kappa = 1/eps ~ 4.5e15 (double precision)
- **Strongly clustered "s" variant**: kappa = 1/sqrt(eps) ~ 6.7e7
- **Strongly clustered "e" variant**: kappa = 1/eps ~ 4.5e15

The paper uses kappa = 1/eps as the default for all matrix classes unless otherwise noted. The "e" and "s" suffixes apply only to the S1/S2 class to distinguish between two levels of clustering tightness.

For Wilkinson and Glued Wilkinson matrices, the condition number is not a parameter -- it is determined by the matrix structure (W_{2m+1}^+ has condition number that grows with m).

---

## 7. The Stetester Framework

The testing infrastructure is described in technical report LBNL-61831 (ref [11]):
> J. W. Demmel, O. A. Marques, B. N. Parlett, and C. Vomel, "A testing infrastructure for symmetric tridiagonal eigensolvers," Technical report LBNL-61831, Lawrence Berkeley National Laboratory, Berkeley, CA, 2006.

### Framework Components

1. **Matrix Generation**:
   - Practical matrices: reduced from sparse matrices via `sytrd` (small) or Lanczos (large)
   - Synthetic matrices: generated via LAPACK's `DLATMS` with specified eigenvalue distributions
   - Wilkinson and Glued Wilkinson: constructed directly

2. **Test Execution**:
   - Tests run on 8 architectures (Power 3, Power 5, Sun UltraSparc 2i, MIPS R12000, Itanium 2, Pentium 4 Xeon, Cray X1, Opteron)
   - Each architecture has specific compiler, timer (PAPI, ETIME, CPU_TIME), and BLAS library

3. **Metrics** (Section 3.3):
   - **Orthogonality loss**: O(Z) = max_{i != j} |z_i^T z_j| / (n * eps)
   - **Residual norm**: R(Lambda, Z) = max_i ||T z_i - lambda_i z_i|| / (||T|| * n * eps)
   - Both normalized by n*eps; a "satisfactorily accurate" algorithm should have these bounded by a modest constant

4. **Performance Measurement**:
   - Wall-clock time via platform-specific timers
   - Flop counts via PAPI hardware performance counters (on platforms that support it)
   - Time trends fitted as t = c * n^alpha on log-log plots for n >= 363 (Opteron L3 cache boundary) or n >= 500 (cross-platform)
   - Deflation fraction fr_defl tracked for DC

5. **Platforms** (Table 3.1):

| Architecture | Symbol | MHz | OS | Compiler | Timer | BLAS |
|---|---|---|---|---|---|---|
| Power 3 | SP3 | 375 | AIX | IBM xlf90 -O3 | PAPI | ESSL |
| Power 5 | SP5 | 1900 | AIX | IBM xlf90 -O3 | PAPI | ESSL |
| Sun UltraSparc 2i | SUN | 650 | Solaris | SUN f90 forte 7.0 -O4 | CPU_TIME | SUNPERF |
| MIPS R12000 | SGI | 600 | IRIX | MIPS pro 7.3.1.3m -O2 | ETIME | SCS |
| Itanium 2 | ITN2 | 1400 | Linux | Intel ifort 9.0 -O2 | ETIME | MKL |
| Pentium 4 Xeon | P4 | 4000 | Linux | Intel ifort 9.0 -O3 | ETIME | MKL |
| Cray X1 | X1 | 800 | UNICOS/mp | Cray ftn 5.4.0.4 -O2 | CPU_TIME | LIBSCI |
| Opteron | OPT | 2200 | Linux | Pathscale pathf90 2.1 -O3 | CPU_TIME | ACML |

---

## 8. Matrices That Broke Specific Algorithms

### Matrices That Break MR (LAPACK 3.0)
- 22 test matrices caused complete failure (no answer returned), including 9 practical matrices
- Large residuals (up to 10^14 * n*eps) and orthogonality losses (up to 10^14 * n*eps) on others
- Root causes: (a) no off-diag satisfies splitting but eigenvalues still numerically equal down to underflow threshold; (b) internal relative gap threshold too loose

### Matrices Most Difficult for MR (LAPACK 3.1)
1. **Glued Wilkinson (GW)**: Worst performance and most overhead. Eigenvalue clusters from p copies of W_{2m+1}^+ glued by sqrt(eps). Representation tree becomes very broad. Time/n^2 is highest among all matrix classes.
2. **Wilkinson (Wi)**: Large Wi matrices cause underflow in FP vector computation. MRRR returns bisectors instead of true eigenvectors. Tight eigenvalue pairs with extreme relative gaps.
3. **Strongly clustered S2ne, S2ns**: All eigenvalues near 1 except one at 1/kappa. Deep representation trees needed.
4. **Strongly clustered S1ns**: MODE=-1, kappa=1/sqrt(eps), with flop count/n^2 slowly increasing (approaching O(n^2) * log factor)

### Matrices Where DC Is Slowest
- **U1, U2** (uniform): Zero or near-zero deflation (fr_defl ~ 0 to 0.03). DC must do full O(n^3) work via BLAS 3 GEMM.
- **W2, W3** (weakly clustered): Very low deflation (0.01-0.09). DC approaches O(n^3).
- These are exactly the matrices where MR excels (eigenvalues not strongly clustered).

### Matrices Where DC Is Fastest
- **S1** (strongly clustered, MODE=1): fr_defl = 0.88 to 1.00. Nearly all eigenvalues deflated. DC faster than O(n^2).
- **GW** (glued Wilkinson): fr_defl = 0.59 to 0.78. High deflation + fast BLAS 3.
- **Wi** (Wilkinson): fr_defl = 0.35 to 0.84.

### Matrices Where BI Struggles
- **NWChem practical matrices**: Clustered eigenvalues require massive reorthogonalization, making BI very slow (up to 300x slower than MR).
- **Highest orthogonality loss**: SUN platform, BI shows orthogonality loss of 440 * n*eps (worst across all algorithms/platforms).

### The DC vs MR Tradeoff (Summary)
- DC wins when deflation is high (S1, GW, Wi): deflation reduces work to O(n) or less per eigenvalue
- MR wins when eigenvalues are well-separated (U1, U2, G1, G2, Prac): O(n^2) with low constant
- Neither dominates: performance depends on eigenvalue distribution

---

## 9. Summary Conclusions from Paper

1. DC and MR are generally much faster than QR and BI on large matrices.
2. MR almost always does the fewest flops but at lowest MFlop rate (scalar divides dominate).
3. Performance of MR and DC strongly depends on the matrix. DC benefits from deflation; MR benefits from well-separated eigenvalues.
4. QR and DC are most accurate: O(sqrt(n)*eps). BI and MR: O(n*eps).
5. MR is preferable to BI for subset computations.
6. LAPACK 3.1 MR is significantly more reliable than 3.0 (fixed splitting criterion bug and tightened accuracy threshold from 1/n to 10^{-3}).

---

## 10. Complete List of Synthetic Test Matrix Symbols

| Symbol | Description | DLATMS MODE | kappa | Eigenvalue Signs |
|---|---|---|---|---|
| S1pe | Strongly clustered, n-1 at 1/kappa, 1 at 1 | 1 | 1/eps | positive |
| S1ne | Same but negated | -1 | 1/eps | negated |
| S1ps | Same but sqrt condition | 1 | 1/sqrt(eps) | positive |
| S1ns | Same but sqrt + negated | -1 | 1/sqrt(eps) | negated |
| S2pe | Strongly clustered, n-1 at 1, 1 at 1/kappa | 2 | 1/eps | positive |
| S2ne | Same but negated | -2 | 1/eps | negated |
| S2ps | Same but sqrt condition | 2 | 1/sqrt(eps) | positive |
| S2ns | Same but sqrt + negated | -2 | 1/sqrt(eps) | negated |
| W1p | Weakly clustered type 1 | custom | 1/eps | positive |
| W1n | Same but negated | custom | 1/eps | negated |
| W2p | Weakly clustered type 2 | custom | 1/eps | positive |
| W2n | Same but negated | custom | 1/eps | negated |
| W3p | Weakly clustered type 3 | custom | 1/eps | positive |
| W3n | Same but negated | custom | 1/eps | negated |
| G1p | Geometric, exact | 3 | 1/eps | positive |
| G1n | Same but negated | -3 | 1/eps | negated |
| G2p | Geometric, random | 3 | 1/eps | random/positive |
| G2n | Same but negated | -3 | 1/eps | negated |
| U1p | Uniform, exact | 4 | 1/eps | positive |
| U1n | Same but negated | -4 | 1/eps | negated |
| U2p | Uniform, random | 4 | 1/eps | random/positive |
| U2n | Same but negated | -4 | 1/eps | negated |
| Wi | Wilkinson W_{m+1}^+ | direct construction | intrinsic | intrinsic |
| GW | Glued Wilkinson (p copies, gamma=sqrt(eps)) | direct construction | intrinsic | intrinsic |

Note: "positive" suffix "p" means all eigenvalues are positive (no random sign flipping).
"negated" suffix "n" means eigenvalues are negated (all signs flipped).
For G2 and U2, the "randomly" refers to random ordering of eigenvalues (MODE reversal plus random signs for symmetric).

---

## 11. Subset Performance (Table 5.1)

Performance of BI relative to MR for subset computations (Time_BI / Time_MR):

| Architecture | By Index (min/med/max) | By Interval (min/med/max) |
|---|---|---|
| SP3 | 0.24 / 1.3 / 8.5 | 0.22 / 1.3 / 9.0 |
| SP5 | 0.22 / 1.2 / 6.5 | 0.20 / 1.2 / 6.8 |
| SUN | 0.38 / 1.7 / 25.2 | 0.32 / 1.7 / 27.4 |
| SGI | 0.30 / 1.3 / 15.2 | 0.29 / 1.4 / 14.2 |
| ITN2 | 0.24 / 1.1 / 4.7 | 0.22 / 1.1 / 4.9 |
| P4 | 0.19 / 1.1 / 7.2 | 0.16 / 1.1 / 7.8 |
| X1 | 0.39 / 1.4 / 4.3 | 0.31 / 1.5 / 4.5 |
| OPT | 0.30 / 1.3 / 21.1 | 0.28 / 1.3 / 16.7 |

MR is faster than BI on average (median ratio > 1 means BI takes longer). BI can be up to 6x faster for "easy" subsets but up to 27x slower for NWChem-type matrices.

---

## 12. Deflation Fractions in DC (Table 4.3)

| Matrix Class | Symbol | fr_defl min | fr_defl max |
|---|---|---|---|
| Strongly clustered S1 | S1 | 0.98 | 1.00 |
| Strongly clustered S2 | S2 | 0.88 | 0.98 |
| Weakly clustered W1 | W1 | 0.34 | 0.81 |
| Weakly clustered W2 | W2 | 0.01 | 0.05 |
| Weakly clustered W3 | W3 | 0.06 | 0.09 |
| Geometric G1 | G1 | 0.16 | 0.36 |
| Geometric G2 | G2 | 0.16 | 0.19 |
| Uniform U1 | U1 | 0 | 0.03 |
| Uniform U2 | U2 | 0 | 0.03 |
| Wilkinson | Wi | 0.35 | 0.84 |
| Glued Wilkinson | GW | 0.59 | 0.78 |
| Practical matrices | Prac | -- | max 0.502 (median 0.125) |
