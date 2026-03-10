# Grosser-Lang 2001: "An O(n^2) Algorithm for the Bidiagonal SVD" -- Complete Analysis

**Paper**: Linear Algebra and its Applications 358 (2003) 45--70
**Authors**: Benedikt Grosser (Wuppertal), Bruno Lang (RWTH Aachen)
**Code**: `hgbsvd/v2/` -- Fortran 77 implementation (circa 1999--2001)

---

## 1. Complete Algorithm from the 2001 Paper

### 1.1 Problem Statement

Given an unreduced n x n upper bidiagonal matrix B = U Sigma V^T, compute singular values sigma_1, ..., sigma_n and singular vector pairs (u_j, v_j).

### 1.2 Key Insight: Why Black-Box RRR Fails for bSVD

Using MR3 (RRR) as a black box on B^T B, BB^T, or T_GK separately **cannot** fulfill the bSVD accuracy requirements:

- **Normal equations (B^T B, BB^T)**: Applying RRR to B^T B = V Sigma^2 V^T and BB^T = U Sigma^2 U^T separately produces U and V that are individually orthogonal, but the resulting pairs may have large residuals ||Bv - sigma*u|| because the two eigensolves use different shift paths and build different subspace bases for tight clusters.

- **Golub-Kahan (T_GK)**: The 2n x 2n matrix T_GK has eigenvalues {-sigma_j, +sigma_j}. Factorizing T_GK - mu*I near small singular values produces large element growth in LDL^T pivots (every 1-minor of T_GK is singular). Also, inner loops run over 2n instead of n, doubling work.

### 1.3 The Coupling Solution

The algorithm works **explicitly on B^T B** for the RRR tree but keeps **BB^T factorizations implicitly** via coupling transformations. This ensures:

1. path(B^T B, j) = path(BB^T, j) = path(T_GK, j) -- same shift sequence for both
2. Coupling uses only multiplications and divisions (no additions/subtractions), making it backward stable
3. Inner loops run over n (not 2n)

### 1.4 Algorithm 4.1: The RRR Algorithm with Embedded Couplings (Step by Step)

```
Input:  B^T B = L_hat D_hat L_hat^T and eigenvalues lambda_hat_1 <= ... <= lambda_hat_n
        (lambda_hat_j = sigma_j^2)
Output: (u, v) -- singular vector pair approximation for each j

For each eigenvalue lambda_hat_j:
  path(B^T B, j) := (nu^(1), ..., nu^(rec_j))

  1. rec = 1
  2. WHILE lambda_hat_j lies in a cluster lambda_hat_f, ..., lambda_hat_l DO:
     3. Choose shift parameter nu^(rec)
     4. Factorize L_hat D_hat L_hat^T - nu^(rec) I = L_hat^+ D_hat^+ (L_hat^+)^T  [EXPLICITLY]
     5. IF rec == 1 THEN
     6.    Use Corollary 2.4 to determine [D_check^+, L_check^+] IMPLICITLY
     7. ELSE
     8.    Use Algorithm 3.2 to determine [D_hat^+, L_hat^+, D_check^+, L_check^+, S_check] IMPLICITLY
     9. END IF
    10. Determine eigenvalue approximations lambda_hat_f^+ <= ... <= lambda_hat_l^+
        for L_hat^+ D_hat^+ (L_hat^+)^T
    11. Determine eigenvalue approximations lambda_check_f^+ <= ... <= lambda_check_l^+
        for L_check^+ D_check^+ (L_check^+)^T
    12. Assess coupling quality by comparing lambda_hat_j^+ and lambda_check_j^+
    13. IF coupling quality is not sufficient THEN
    14.    Go back to step 3 (try different shift)
    15. END IF
    16. Update D_hat <- D_hat^+, L_hat <- L_hat^+, lambda_hat_j <- lambda_hat_j^+
    17. Update D_check <- D_check^+, L_check <- L_check^+, D_hat_tilde <- D_hat_tilde^+, L_hat_tilde <- L_hat_tilde^+
    18. rec = rec + 1
  19. END WHILE
  20. nu^(rec) = lambda_hat_j  (final shift = eigenvalue approximation)

  {Set up twisted factorizations and solve for vectors}
  21. L_hat D_hat L_hat^T - nu^(rec) I = N_hat_k^+ G_hat_k^+ (N_hat_k^+)^T
  22. Solve (N_hat_k^+ G_hat_k^+ (N_hat_k^+)^T) * v = gamma_hat_k^+ * e_k  and normalize v
  23. L_check D_check L_check^T - nu^(rec) I = N_check_k^+ G_check_k^+ (N_check_k^+)^T
  24. Solve (N_check_k^+ G_check_k^+ (N_check_k^+)^T) * u = gamma_check_k^+ * e_k  and normalize u
```

### 1.5 Coupling Formulas

**For positive definite initial matrices (rec == 1, Section 2):**

Using Algorithm 2.1 (differential qd), compute factorizations of B^T B - mu^2 I and BB^T - mu^2 I simultaneously from B and mu. The auxiliary variables s_hat_i (stationary) and p_hat_i (progressive) from differential qd enable coupling with **no additions or subtractions**:

Corollary 2.4 gives direct coupling Z_hat <-> Z_check:
```
d_check_i = s_hat_{i+1} * d_hat_i / s_hat_i
p_check_i = -mu^2 * d_hat_i / s_hat_i
r_check_i = p_hat_i * r_hat_{i+1} / p_hat_{i+1}
s_check_i = -mu^2 * r_hat_{i+1} / p_hat_{i+1}
gamma_check_i * s_check_i = gamma_hat_i * s_hat_i  (twist element relation)
```

**For indefinite initial matrices (rec > 1, Section 3):**

Given LDL^T factorizations of translates, coupling via Lemma 3.1:
```
d_hat_i^+ = -d_tilde_{2i-1}^+ * d_tilde_{2i}^+    (diagonal pivots)
l_hat_i^+ = -l_tilde_{2i-1}^+ * l_tilde_{2i}^+    (off-diagonal elements)
```
Corollary 3.2 gives twisted factorization couplings:
```
r_hat_{i+1}^+ = -r_tilde_{2i}^+ * r_tilde_{2i+1}^+
u_hat_i^+ = -u_tilde_{2i-1}^+ * u_tilde_{2i}^+
gamma_check_i^+ = (mu_tilde + mu) * gamma_tilde_{2i-1}^+  (twist elements)
```

**The X_k matrix (Remark 3.3):**
Instead of u = (1/sigma) B v, use u = X_k v where:
```
X_k = diag( l_tilde_{2k-1}^+/l_tilde_1^+, ..., l_tilde_{2k-1}^+/l_tilde_{2k-3}^+,
             1,
             u_tilde_{2k+1}^+/u_tilde_{2k-1}^+, ..., u_tilde_{2n-1}^+/u_tilde_{2k-1}^+ )
```
This diagonal scaling provides componentwise relative accuracy for u given v, unlike the naive u = Bv/sigma which loses accuracy for small sigma.

### 1.6 Why Couplings Produce Good Results

1. **Same shift paths**: path(B^T B, j) = path(BB^T, j) guarantees aligned subspaces
2. **Backward stable couplings at rec=1**: Corollary 2.4 uses only mults/divs
3. **Quality assessment at deeper levels**: Comparing eigenvalue approximations from explicit vs. coupled factorizations detects coupling failures
4. **X_k conversion**: Gives u with small relative error component-by-component

---

## 2. How DLARRI Works

**File**: `dlarri.f`
**Purpose**: Replaces DLARRE from LAPACK's DSTEGR. Computes eigenvalues of B^T B (or BB^T) and sets up Gershgorin intervals.

### What DLARRI Computes

1. **Splitting points**: Scans for |e_i| <= TOL to split B into blocks (but returns INFO=33 if any splitting or negative entry is found -- splittings NOT supported)

2. **Auxiliary variables** (precomputed for later coupling use):
   - `D2(i) = d_i^2` -- squared diagonal entries
   - `E2(i) = e_i^2` -- squared off-diagonal entries
   - `DE00(i) = d_i * e_i` -- off-diagonal of B^T B
   - `DE10(i) = d_{i+1} * e_i` -- off-diagonal of BB^T

3. **Gershgorin intervals**: For WANTR or WANTA, forms T = B^T B explicitly and computes intervals. For WANTL or WANTA, forms T = BB^T. For WANTA, unifies both sets of intervals.

4. **LDL^T factorization**:
   - If WANTL: computes BB^T = LDL^T using the differential stationary qd:
     ```
     tmp = d_1^2
     for i = 1 to n-1:
       work(i) = tmp + e_i^2     (diagonal of D)
       work(n+i) = d_{i+1}*e_i / work(i)  (subdiagonal of L)
       tmp = d_{i+1}^2 * tmp / work(i)
     work(n) = tmp
     ```
   - If WANTR or WANTA: computes B^T B = LDL^T:
     ```
     work(i) = d_i^2             (diagonal of D)
     work(n+i) = e_i / d_i       (subdiagonal of L)
     ```

5. **Eigenvalue computation**: Calls DLASQ2 (dqds algorithm) on the LDL^T representation to find all eigenvalues to high relative accuracy. Returns them sorted in ascending order.

### Key Difference from DLARRE

- DLARRI does NOT compute initial shift (the matrices are already positive definite)
- DLARRI computes and stores auxiliary variables D2, E2, DE00, DE10 needed for coupling
- DLARRI checks that all d_i > 0 and e_i > 0 (no negative entries allowed) -- returns INFO=33 otherwise
- No Sturm count or bisection; eigenvalues come purely from dqds

---

## 3. How DBDSGR.f Implements the Algorithm

### Call Flow

```
DBDSGR (main driver)
  |
  +-- Sign handling: make all d_i, e_i positive (store signs for later)
  +-- Save originals: INDAOR = a_orig, INDBOR = b_orig
  +-- Precompute: INDD2 = d^2, INDE2 = e^2, INDDE0 = d*e, INDDE1 = d_{i+1}*e
  |
  +-- DLARRI (eigenvalue computation + LDL^T + Gershgorin)
  |     |
  |     +-- DLASQ2 (dqds for eigenvalues)
  |
  +-- DLARRV (eigenvector computation with embedded couplings)
  |     |
  |     +-- DLARRB (bisection refinement for deeper levels)
  |     +-- DLARRC (child representation WITH coupling -- replaces DLARRF for COUP=true)
  |     |     |
  |     |     +-- Computes DCOUP, LCOUP from explicit [D+, L+] using coupling formulas
  |     |
  |     +-- DLARRF (child representation WITHOUT coupling -- standard RRR)
  |     +-- DLAR1V (twisted factorization + eigenvector via inverse iteration)
  |     |     |
  |     |     +-- DLATS1 (solve twisted system N G N^T z = gamma e_R)
  |     |
  |     +-- DLATS2 (eigenvector refinement: solve N G N^T x = z, then G z = x, then N^T x = z)
  |     +-- DLACSV (compute coupled singular vector u from v)
  |           |
  |           +-- IF NDEPTH==0: DLAS2P -> DLAP2S -> DLAG2G -> DLATS1 -> DLATS2
  |           +-- IF NDEPTH>0:  DLAR1V -> DLATS2
  |           +-- Sign arrangement (make sign(u_R) match sign((Bv)_R))
  |
  +-- Take sqrt of eigenvalues to get singular values
  +-- Rescale if matrix was scaled
  +-- Apply stored signs to U and V
  +-- Sort singular values/vectors in descending order
```

### Storage Layout

- `Z(1:N, 1:M)` -- right singular vectors V (rows 1..N)
- `Z(LDZ/2+1:LDZ/2+N, 1:M)` -- left singular vectors U (rows LDZ/2+1..LDZ/2+N)
- Workspace: 26*N doubles, 10*N integers
- Workspace pointers:
  - `INDGRS  =  1` (2N: Gershgorin intervals)
  - `INDSGL  =  2N+1` (N: signs of d_i)
  - `INDSGR  =  3N+1` (N: signs of e_i)
  - `INDAOR  =  4N+1` (N: original |d_i|)
  - `INDBOR  =  5N+1` (N-1: original |e_i|)
  - `INDD2   =  6N+1` (N: d_i^2)
  - `INDE2   =  7N+1` (N-1: e_i^2)
  - `INDDE0  =  8N+1` (N-1: d_i * e_i)
  - `INDDE1  =  9N+1` (N-1: d_{i+1} * e_i)
  - `INDDC   = 10N+1` (N: coupled D diagonal)
  - `INDLC   = 11N+1` (N-1: coupled L subdiagonal)
  - `INDWRK  = 12N+1` (14N: DLARRV workspace)

### Sign Handling

DBDSGR makes B nonnegative before processing:
```
For i = 1 to n-1:
  sgl_i = sign(d_i);  d_i = |d_i|
  e_i = sgl_i * e_i
  sgr_i = sign(e_i);  e_i = |e_i|
  d_{i+1} = sgr_i * d_{i+1}
```
After computing vectors, signs are restored by negating rows of U or V as needed.

---

## 4. Test Matrices in DMATGEN.f

### Spectrum-Based Matrices (IDs 100-121)

Construct B^T B with prescribed eigenvalues via DLATMS (Householder reduction of diagonal to bidiagonal). Epsilon = machine precision.

| ID  | Name | Formula |
|-----|------|---------|
| 110 | Ones | lambda_j = 1 for all j |
| 111 | Uniform (eps apart) | lambda_j = j * eps, j=1..n-1; lambda_n = 1 |
| 112 | Uniform (sqrt(eps) apart) | lambda_1 = eps; lambda_j = 1 + j*eps^(1/4), j=2..n-1; lambda_n = 2 |
| 113 | Uniform (eps to 1) | lambda_j = eps + (j-1)*(1-eps)/(n-1) |
| 114 | Uniform (eps to 1, random signs) | Same as 113 but with random sign flips |
| 115 | Geometric (eps to 1) | lambda_j = eps^((n-j)/(n-1)) |
| 116 | Geometric (eps to 1, random signs) | Same as 115 but with random sign flips |
| 117 | Random | lambda_j = DLARAN(seed), uniform in (0,1) |
| 118 | Clustered at 1 | lambda_1 = eps^2; lambda_j = (10n + rand)/(10n) for j >= 2 |
| 119 | Clustered at 1, -1 | lambda_1 = eps; lambda_j = +/- 1 randomly for j >= 2 |
| 120 | Clustered at eps | lambda_j = eps^2 * (10n + rand)/(10n) for j=1..n-1; lambda_n = 1 |
| 121 | Clustered at eps, -eps | lambda_j = +/- eps randomly for j=1..n-1; lambda_n = 1 |

DLATMS modes 1-6 follow LAPACK LAWN 9:
| Mode | Description |
|------|-------------|
| 1 | D(1)=1, D(2:N)=1/COND |
| 2 | D(1:N-1)=1, D(N)=1/COND |
| 3 | D(I) = COND^(-(I-1)/(N-1)) (geometric) |
| 4 | D(I) = 1 - (I-1)/(N-1)*(1-1/COND) (arithmetic) |
| 5 | Random in (1/COND, 1), log-uniform |
| 6 | Random from same distribution as matrix |

### Entry-Based Matrices (IDs 200-244)

These directly specify bidiagonal entries:

| ID  | Name | Diagonal a_i | Off-diagonal b_i | Notes |
|-----|------|-------------|------------------|-------|
| 200 | ABCON0 | alpha (=2) | beta (=1) | 1D heat equation tridiagonal |
| 201 | ABCON1 | alpha, a_1=alpha-beta | beta | |
| 202 | ABCON2 | alpha, a_1=alpha-beta, a_n=alpha+beta | beta | |
| 203 | ABCON3 | alpha, a_1=a_n=alpha+beta | beta | |
| 210 | RANDOM | dlaran(seed) | dlaran(seed) | Direct bidiagonal entries |
| 220 | GRADP | a_n=1, a_i=beta*a_{i+1} | b_i=a_i | Graded, biggest at top |
| 221 | GRADM | a_1=1, a_{i+1}=beta*a_i | b_i=a_{i+1} | Graded, biggest at bottom |
| 222 | WILKP | Wilkinson+ | b_i=1 | Symmetric W-shape |
| 223 | WILKM | Wilkinson- | b_i=1 | Antisymmetric W-shape |
| 224 | WILKW | |1-n/2+i| | b_i=1 | Fernando-Parlett Wilkinson variant |
| 225 | WILK2W | Double WILKW | b_i=1 | Two copies of WILKW spectrum |
| 230 | CLEM | a_i=0 | sqrt(i*(n-i)) | Clement matrix |
| 240 | GRO0 | a_1=1, a_i=alpha (i>1) | b_i=alpha | Near-constant with one outlier |
| 241 | GRO1 | a_1=a_2=1, rest=alpha | b_i=alpha | Two outliers |
| 242 | GRO2 | a_1..a_4=1, rest=alpha | b_i=alpha | Four outliers |
| 243 | GRO3 | a_1..a_4=1, rest=alpha | b_1..b_4=1, rest=alpha | Outlier block |
| 244 | bWILK+ | Wilkinson+ + 1 | b_i=1 | Shifted Wilkinson bidiagonal |

**Tridiagonal-to-bidiagonal conversion** (IDs 200-205, 222, 223, 230):
For matrices given as tridiagonal T, compute Cholesky T - val*I = B^T B where val = 1.001 * gl (if T is indefinite) or val = 0 (if T is positive definite).

**Gluing**: Base matrices can be repeated with `gluef` copies, joined by a small element `small` at b_{k*basedim}.

---

## 5. Evaluation Criteria in DTEST.f

### Metrics Computed by DCHECK_RESULTS

1. **Orthogonality of U** (orthu):
   ```
   orthu = || U^T U - I ||_1    (1-norm)
   maxu  = max |[U^T U - I]_{ij}|
   ```
   Reported as: `orthu / (m * eps)` -- should be O(1)

2. **Orthogonality of V** (orthv):
   ```
   orthv = || V^T V - I ||_1    (1-norm)
   maxv  = max |[V^T V - I]_{ij}|
   ```
   Reported as: `orthv / (n * eps)` -- should be O(1)

3. **Residual** (res):
   ```
   For each i = 1..n:
     r_i = B * v_i - sigma_i * u_i
     res = max_i || r_i ||_2
   ```
   Reported as: `res / (min(m,n) * eps * sigma_max)` -- should be O(1)

4. **Componentwise checks**:
   - `maxr = max_{i,j} |[U^T B V - Sigma]_{ij}|`
   - `maxrc = maxr / (n * eps * sigma_max)`

### Paper's Results Table (Table 3)

Residual norms `max_j {||B v_j - sigma_j u_j|| / (n * eps * |lambda_n|)}`:

| Matrix | n | DSTEGR | RRR/coup | DBDSDC | DBDSQR |
|--------|---|--------|----------|--------|--------|
| Geometric | 500 | 0.019 | 0.077 | 0.018 | 0.073 |
| Geometric | 2000 | 0.017 | 0.060 | 0.007 | 0.055 |
| Arithmetic | 500 | 0.043 | 0.018 | 0.021 | 0.082 |
| Internal cluster | 500 | 0.013 | 0.278 | 0.032 | 0.131 |
| Internal cluster | 2000 | 0.010 | 0.113 | 0.019 | 0.090 |
| 1-2-1 | 2000 | 0.051 | 0.101 | 0.008 | 0.056 |
| Wilkinson | 501 | NaN | 0.062 | 0.040 | 0.105 |
| Wilkinson | 2001 | NaN | 0.108 | 0.007 | 0.106 |

Key observations:
- DSTEGR (black-box on T_GK) produces NaN for Wilkinson matrices
- RRR/coup residuals are moderate (< 0.3) but higher than DBDSDC for clustered cases
- Internal cluster (type 3) is the hardest case -- residuals reach 0.278

---

## 6. Known Bugs and Limitations

### In the Paper

1. **Only positive definite couplings implemented**: The paper derives coupling formulas for both positive definite (Section 2, rec=1) and indefinite (Section 3, rec>1) initial matrices, but the code comment says: *"Only couplings for so-called positive definite initial matrices are implemented. Thus the routine may return no results for very tight clusters of singular values."* This is flagged by INFO=3, 4, or 5.

2. **Coupling quality at deeper recursion levels**: For indefinite initial matrices (deeper recursion), Algorithm 3.2 uses a mixed explicit/coupled approach where odd-numbered elements of D_check^+ require one dstqds step. This can lose information: even if relative error in even-numbered elements is small, odd-numbered elements can have large deviations (Remark 3.4).

3. **Internal cluster residuals**: The RRR/coup method shows higher residuals (0.278) for internal cluster matrices. The paper acknowledges this is the hardest case.

4. **DSTEGR NaN issues**: The LAPACK DSTEGR (at the time) produced NaN for Wilkinson-type matrices.

### In the Code

1. **INFO=33 -- Splittings not supported**: DLARRI returns INFO=33 if any d_i <= 0 or e_i <= 0 or any element is below tolerance. This means the code CANNOT handle matrices with near-zero diagonal or off-diagonal entries.

2. **INFO=3 -- No backward stable coupling at deeper levels**: DLARRC returns INFO=3 when NDEPTH > 0, meaning no backward-stable coupling transformation is available for indefinite initial matrices. The code falls back to... nothing (returns error).

3. **INFO=4, 5 -- No couplings for DSTEIN or MGS branches**: When DLARRF fails to find a new RRR and the minimum relative gap is too small, or when MGS is needed for a small cluster, the code has no coupling implementation and returns INFO=4 or INFO=5.

4. **RANGE='A' only**: Only computing all singular values is supported. RANGE='V' and RANGE='I' are not implemented.

5. **Hardcoded sign negation in dtest.f**: The test driver contains `A_ORIG(3) = -A_ORIG(3)` and `B_ORIG(6) = -B_ORIG(6)` -- this modifies the test matrix before solving, which is a testing artifact/bug left in the code.

6. **Maximum matrix size**: MMAX=2100 in dtest.f limits testing.

7. **Debug write statement left in dbdsgr.f**: Line 430 `write(*,*) TIMNG( 7 )` prints timing info unconditionally.

---

## 7. How This Differs from Willems-Lang 2013

### Grosser-Lang 2001 (this paper/code)

- **Approach**: Work on B^T B explicitly with RRR, couple to BB^T implicitly
- **Eigenvalue computation**: dqds (DLASQ2) on B^T B or BB^T
- **Eigenvector computation**: Standard DLARRV with embedded coupling calls
- **Coupling mechanism**: Analytical formulas from differential qd auxiliary variables (s_i, p_i)
- **Positive definite case**: Fully backward stable (Corollary 2.4 -- no additions/subtractions)
- **Indefinite case**: Partially backward stable (some dstqds steps needed for odd elements)
- **Limitation**: Fails for very tight clusters where coupling quality degrades

### Willems-Lang 2013 (Algorithm 4.1 in that paper)

- **Approach**: NCD (Narrow Cluster Deflation) -- explicitly form Golub-Kahan T_GK at deeper recursion levels when normal equations fail
- **Key innovation**: When B^T B has a tight cluster of eigenvalues, switch to T_GK representation which naturally separates +sigma from -sigma eigenvalues
- **Dual eigensolver**: Can solve on either B^T B or T_GK depending on cluster properties
- **Wider applicability**: Handles the tight-cluster failure mode that Grosser-Lang cannot
- **Implementation status**: The XMR code has DSTEXR (tridiagonal eigensolver) but DBDSGR was left incomplete -- the DLARRI subroutine was missing from the XMR distribution

### Critical Difference

Grosser-Lang 2001 relies on coupling transformations to get both U and V from a single eigensolve on B^T B. When clusters are very tight (relative gap < machine precision), the coupling quality degrades and the method fails gracefully (returns INFO=3/4/5).

Willems-Lang 2013 proposes NCD as the solution: for tight clusters, form T_GK - mu*I = LDL^T explicitly and work directly on the 2n x 2n system, where the +-sigma pairs are naturally separated. This avoids the coupling instability entirely but requires different machinery.

---

## 8. Complete Call Graph

```
DBDSGR(JOBZ, RANGE, N, D, E, ...)
  |
  +-- LSAME (character comparison)
  +-- DLAMCH ('Precision')
  +-- DLANST ('M', max-norm of B)
  +-- DSCAL (scaling)
  +-- DCOPY (save originals, signs)
  +-- DSECND (timing)
  |
  +-- DLARRI(N, WANTL, WANTR, WANTA, D, E, TOL, ...)
  |     |
  |     +-- DCOPY (save/restore d, e)
  |     +-- DLASQ2(IN, WORK, INFO)  -- dqds eigenvalues
  |
  +-- DLARRV(N, D, L, ISPLIT, M, W, IBLOCK, GERSCH, TOL, COUP, Z, ...)
  |     |
  |     +-- DLAMCH, DLASET, DCOPY
  |     |
  |     +-- DLARRB(IN, D, L, LD, LLD, ...)  -- bisection refinement
  |     |
  |     +-- [COUP=true] DLARRC(IN, D, L, LD, ...)  -- coupled child rep
  |     |     |
  |     |     +-- DCOPY
  |     |
  |     +-- [COUP=false] DLARRF(IN, D, L, LD, LLD, ...)  -- standard child rep
  |     |     |
  |     |     +-- DLAMCH, DCOPY
  |     |
  |     +-- DLAR1V(B1, BN, SIGMA, D, L, LD, LLD, ...)  -- twisted factorization + eigenvector
  |     |     |
  |     |     +-- DLAMCH, DSECND
  |     |     +-- DLATS1(Z, ZTZ, L, U, LD, B1, R, BN, ...)  -- solve twisted system
  |     |
  |     +-- DLATS2(Z, ZTZ, MINGMA, L, U, LD, ...)  -- eigenvector refinement
  |     |
  |     +-- [COUP=true] DLACSV(N, B1, R, BN, NDEPTH, ...)  -- coupled singular vector
  |     |     |
  |     |     +-- [NDEPTH=0] DLAS2P(...)  -- coupling L+ -> Lcoup (stationary)
  |     |     +-- [NDEPTH=0] DLAP2S(...)  -- coupling U- -> Ucoup (progressive)
  |     |     +-- [NDEPTH=0] DLAG2G(...)  -- find coupled twist position
  |     |     +-- [NDEPTH=0] DLATS1(...)  -- solve coupled twisted system
  |     |     +-- [NDEPTH=0] DLATS2(...)  -- refine coupled eigenvector
  |     |     +-- [NDEPTH>0] DLAR1V(...)  -- direct from coupled representation
  |     |     +-- [NDEPTH>0] DLATS2(...)  -- refine
  |     |     +-- Sign arrangement
  |     |
  |     +-- [COUP=true, MGS cluster] -> INFO=5 (not implemented)
  |     +-- [COUP=false, tight cluster] DSTEIN  -- fallback inverse iteration
  |     +-- DDOT, DAXPY, DNRM2, DSCAL  -- Gram-Schmidt for small clusters
  |
  +-- SQRT(W) -- eigenvalues -> singular values
  +-- DSCAL (rescale)
  +-- DSCAL (sign arrangement for U rows and V rows)
  +-- Selection sort (descending order)
  +-- DSWAP (reorder vectors)
```

### Routine Purposes Summary

| Routine | Purpose | New/Modified |
|---------|---------|-------------|
| **DBDSGR** | Main bSVD driver (replaces DSTEGR for bSVD) | New |
| **DLARRI** | Eigenvalues + LDL^T + Gershgorin (replaces DLARRE) | New |
| **DLARRV** | Eigenvectors with coupling support | Modified |
| **DLARRC** | Child representation with coupled BB^T (replaces DLARRF for COUP) | New |
| **DLARRF** | Standard child representation (unchanged from LAPACK) | Unchanged |
| **DLAR1V** | Twisted factorization + eigenvector | Modified (FACTRL, split LPL/UMN/SPL/PMN) |
| **DLACSV** | Master coupling routine: v -> u | New |
| **DLAS2P** | Coupling: stationary LDL^T -> coupled LDL^T | New |
| **DLAP2S** | Coupling: progressive URU^T -> coupled URU^T | New |
| **DLAG2G** | Find coupled twist position and element | New |
| **DLATS1** | Solve twisted system N G N^T z = gamma e_R | New (extracted from DLAR1V) |
| **DLATS2** | Eigenvector refinement (one step of iterative refinement) | New |
| **DLARRB** | Bisection refinement | Unchanged |

---

## 9. Detailed Subroutine Analysis

### DLAS2P: Stationary Coupling (B^T B -> BB^T, top-to-bottom)

Given the stationary qd factorization of B^T B - mu I = L^+ D^+ (L^+)^T with auxiliary variables S (from dstqds), computes the coupled representation for BB^T - mu I.

Normal code path (no NaN):
```
For i = B1 to R2-1:
  dfact = AB00(i) / L(i)       -- a_i * b_i / l_i^+
  t = dfact / (S(i) - mu)      -- scaling factor
  dcoup = (S(i+1) - mu) * t    -- coupled diagonal
  pcoup(i) = -mu * t            -- coupled auxiliary
  lcoup(i) = AB10(i) / dcoup   -- coupled off-diagonal = a_{i+1}*b_i / dcoup
```

NaN-proof code handles zero pivots by resetting to b_i^2 and restarting.

### DLAP2S: Progressive Coupling (B^T B -> BB^T, bottom-to-top)

Given the progressive qd factorization with auxiliary variables P (from dqds), computes coupled representation from bottom to top.

```
rcoup = P(BN)
scoup(BN) = -mu
t = 1/P(BN)
For i = BN-1 down to R1:
  ucoup(i) = AB10(i) / rcoup
  rfact = AB00(i) / U(i)
  rcoup = rfact * P(i) * t
  t = 1/P(i)
  scoup(i) = -mu * rcoup * t
```

### DLAG2G: Find Coupled Twist Position

After DLAS2P and DLAP2S produce the coupled stationary (PCOUP via SFACT) and progressive (SCOUP via PFACT) auxiliary variables, DLAG2G finds an acceptable twist position RC for the coupled representation.

The coupled twist element at position i is:
```
gcoup(i) = gfact(i) * scoup(i) / (sfact(i) - mu)
```
where `gfact(i) = sfact(i) + pfact(i)` is the explicit twist element.

Acceptance criterion: `|gcoup(RC)| <= 2 * |gfact(RF)|` where RF is the explicit twist position.

### DLATS2: Eigenvector Refinement

Performs one step of iterative refinement on the twisted factorization N G N^T:
1. Solve N * x = z (forward/backward substitution using L and U)
2. Solve G * z = x (diagonal scaling by 1/gamma_i)
3. Solve N^T * x = z (transpose solve)
4. Copy x -> z, compute z^T z

This improves orthogonality significantly (noted in DLARRV comments).

---

## 10. Complexity Analysis

From Table 2 in the paper, execution times scale approximately as:

| Matrix | n=500 | n=1000 | n=2000 | Ratio 1000/500 | Ratio 2000/1000 |
|--------|-------|--------|--------|----------------|-----------------|
| Geometric | 0.31 | 1.34 | 5.51 | 4.3x | 4.1x |
| Arithmetic | 0.35 | 1.54 | 6.44 | 4.4x | 4.2x |
| Internal cluster | 2.03 | 8.82 | 36.62 | 4.3x | 4.2x |
| 1-2-1 | 0.51 | 2.22 | 9.01 | 4.4x | 4.1x |
| Wilkinson | 0.47 | 1.96 | 7.94 | 4.2x | 4.1x |

Doubling ratios are consistently ~4x, confirming O(n^2) scaling. The internal cluster case is ~6x more expensive than geometric (but still O(n^2)).

For comparison, DBDSQR (O(n^3)) shows ratios of 13-14x, and DBDSDC shows ratios of 5-7x (O(n^2) with larger constant).

---

## 11. Key Insights for Implementation

### What Makes This Work

1. **Differential qd provides the coupling "glue"**: The auxiliary variables s_i (stationary) and p_i (progressive) from dstqds/dqds are the exact quantities needed for backward-stable coupling at the first recursion level.

2. **Same shift path = same subspace**: By using equivalent shifts for B^T B and BB^T, the computed U and V naturally pair up.

3. **Twisted factorization for both sides**: The twist position k for B^T B gives v; the coupled twist position (possibly different) gives u.

4. **Sign arrangement**: After computing u from the coupled eigenproblem, the sign is adjusted so that sign(u_R) matches sign((Bv)_R).

### What Limits This

1. **Deeper recursion coupling is fragile**: At rec > 1, the coupling uses mixed explicit/coupled steps that can amplify errors.

2. **Cannot handle zero/negative entries**: The whole approach assumes B is unreduced with all positive entries (after sign extraction). Splittings are not handled.

3. **Incomplete implementation**: INFO=3/4/5 return paths indicate the algorithm gives up when clusters are too tight or when fallback paths (DSTEIN, MGS) would be needed.

4. **No NCD**: The key innovation of Willems-Lang 2013 (switching to T_GK for tight clusters) was not available in 2001.

### What DLARRI Contains That Was Missing from XMR

The XMR (Willems) distribution's DBDSGR referenced DLARRI but didn't include it. The hgbsvd DLARRI provides:
- Splitting detection (with INFO=33 for unsupported cases)
- LDL^T factorization of B^T B or BB^T using differential qd
- Eigenvalue computation via DLASQ2
- Gershgorin interval computation for both B^T B and BB^T
- Precomputation of auxiliary arrays D2, E2, DE00, DE10

This is the "eigenvalue refinement" step that was identified as missing from the XMR DBDSGR implementation.
