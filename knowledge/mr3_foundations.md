# MR3 Foundations: Bugs, Test Matrices, and Fixes

Synthesized from four foundational papers:
1. Dhillon & Parlett 2004 — "Orthogonal Eigenvectors and Relative Gaps"
2. Dhillon, Parlett & Vomel 2005 — "Glued Matrices and the MRRR Algorithm"
3. Parlett & Dhillon 2000 — "Relatively Robust Representations of Symmetric Tridiagonals"
4. Grosser & Lang 2005 — "On Symmetric Eigenproblems Induced by the Bidiagonal SVD"

---

## 1. Core Theory: Relatively Robust Representations (RRR)

**Source: Parlett & Dhillon 2000**

### Definition
An LDL^t factorization of T - sigma*I is a **Relatively Robust Representation (RRR)** if small relative perturbations to L and D cause only small relative perturbations to eigenvalues. Formally, if L_hat = L*(I + dL), D_hat = D*(I + dD), with |dL|, |dD| <= eps, then each eigenvalue nu_j of T - sigma*I satisfies:

    |nu_j_hat - nu_j| / |nu_j| <= c(n) * eps

where c(n) is a modest function of n.

### Relative Condition Number (relcond)
For eigenvalue lambda_j with shift tau:

    relcond(lambda_j - tau) = || |D_+|^{1/2} * L_+^t * s_j ||^2 / |lambda_j - tau|

where L_+ D_+ L_+^t = T - tau*I and s_j is the corresponding eigenvector.

**Key property**: relcond is the same for ALL n twisted factorizations of T - tau*I (Theorem 1). This is remarkable --- twist index doesn't affect conditioning.

### When is LDL^t an RRR?
**Theorem 3**: LDL^t is an RRR iff every eigenvector s of L D L^t satisfies:

    || |D|^{1/2} L^t s || = O(1)

Equivalently, relcond(nu_j) = O(1) for all j.

### Element Growth and RRR Quality
**Theorem 5** (Cluster relcond bound): For a cluster {lambda_j, ..., lambda_k} with gap_left and gap_right:

    relcond(lambda_j - tau) <= (sum |d_i| + sum |d_i * l_i^2|) / min(gap_left, gap_right)

This means: if the LDL^t factors have large element growth (large |d_i| or |l_i|), then relcond is large, and the representation is NOT robust.

**Theorem 6** (via rho values): For twisted factorization at index r:

    relcond(nu) = nu * sum_{k=1}^{n} |rho_k|^{-1}

where rho_k are the "gamma" values from twisted factorizations at each index k.

### Differential qd Transforms
Three key transforms for computing LDL^t factorizations:

1. **dstqds** (stationary quotient-difference with shift): Computes L_+ D_+ L_+^t = LDL^t - s*I
   - d_plus[1] = d[1] - s
   - For i = 1 to n-1: l_plus[i] = d[i]*l[i] / d_plus[i]; d_plus[i+1] = d[i+1] - s - l[i]*d_plus[i]*l_plus[i] (actually uses l_plus)
   - Correction: d_plus[i+1] = d[i+1]*l[i]^2 * (d[i]/d_plus[i]) - s
   - **Mixed stability**: relative perturbation bound of ~3 ulps per element

2. **dqds** (progressive qd): More stable variant
   - d_plus[1] = d[1] - s, t[1] = d[1]/d_plus[1]
   - l_plus[i] = l[i]*t[i], t[i+1] = d[i+1]/d_plus[i]*(something)
   - Also ~3 ulps per element

3. **dtwqds** (twisted qd): Computes twisted factorization for eigenvector
   - Runs dstqds from left (getting L_+, D_+) and from right (getting L_-, D_-)
   - Combines at twist index r where |gamma_r| is minimized
   - gamma_r = d_plus[r] + d_minus[r] - (T[r,r] - lambda)

### Fernando-Parlett (FP) Vector
Given twisted factorization N_r D_r N_r^t = T - lambda*I at twist index r:
- Solve N_r^t * z = e_r (i.e., z(r) = 1, propagate up/down)
- ||z||_2 = 1/|gamma_r|^{1/2} approximately
- Residual: ||(T - lambda*I)*z|| = |gamma_r| * |z(r)| = |gamma_r|

**Guarantee** (Theorem 15): The angle between computed z_hat and true eigenvector v satisfies:

    |sin angle(z_hat, v)| <= 5n*eps + |lambda_bar - lambda_hat| / (|v_bar(r)| * gap) + eps_rel * relcond(v) / (1 - eps_rel)

where eps_rel ~ 3n*eps (from qd roundoff). The bound has three terms:
1. Direct roundoff in FP vector computation
2. Eigenvalue error amplified by 1/(|v(r)| * gap)
3. Representation quality (relcond)

---

## 2. The Representation Tree (MR3 Algorithm Structure)

**Source: Dhillon & Parlett 2004**

### Algorithm Getvec (simplified)
```
Input: T (tridiagonal), eigenvalues lambda_1, ..., lambda_n
1. Compute root representation: L_0 D_0 L_0^t = T - sigma_0*I
   (sigma_0 chosen so all eigenvalues are well-represented)
2. Build representation tree:
   - Root node: all eigenvalues, represented by L_0 D_0 L_0^t
   - For each cluster of eigenvalues (relative gap < threshold):
     a. Choose shift tau (left or right end of cluster)
     b. Compute child: L_c D_c L_c^t = L_0 D_0 L_0^t - (tau - sigma_0)*I via dstqds
     c. Recurse: if sub-clusters exist, shift again
   - Leaves: singletons (isolated eigenvalues in their representation)
3. For each leaf eigenvalue lambda_j:
   a. Refine lambda_j in its local representation using bisection/Rayleigh
   b. Find twist index r where |gamma_r| is minimized
   c. Compute FP vector z_j
   d. Normalize z_j
   e. Transform back to original basis using the chain of L factors
```

### Relative Gap Definition
    relgap(lambda_j) = min(|lambda_j - lambda_{j-1}|, |lambda_j - lambda_{j+1}|) / |lambda_j|

Threshold: relgap >= 10^{-3} (typical) means "well separated" --- can compute eigenvector directly.

### Shift Selection
- Shift to LEFT end of cluster: tau = lambda_left - relgap_left * |lambda_left|
- Shift to RIGHT end of cluster: tau = lambda_right + relgap_right * |lambda_right|
- Goal: after shifting, eigenvalues within cluster are separated in RELATIVE sense

---

## 3. Documented Bugs and Failure Modes

### Bug 1: Singleton Assumption Failure (Assumption 4.1)
**Source: Dhillon, Parlett & Vomel 2005**

**The assumption**: After O(1) levels of the representation tree, every cluster can be broken into singletons (isolated eigenvalues with relative gap >= threshold).

**How it fails**: For "glued" matrices (defined below), eigenvalues can be so tightly clustered that NO amount of shifting produces relative gaps. The eigenvalues are spaced at O(epsilon) apart, and each shift only reveals O(1) of them as singletons, requiring O(n) levels.

**Consequence**: Algorithm may loop indefinitely or require O(n) tree depth, violating O(n^2) complexity.

### Bug 2: RRR Existence Failure (Assumption 4.2)
**Source: Dhillon, Parlett & Vomel 2005**

**The assumption**: For each cluster, there exists a shift tau such that L D L^t - tau*I is an RRR for that cluster (i.e., small element growth, small relcond).

**How it fails**: When an eigenvalue lambda is also a Ritz value of the leading or trailing principal submatrix with multiplicity > 1, then for shift tau = lambda:
- The factorization T - tau*I has a zero pivot (d_plus[k] = 0 for some k)
- Any nearby shift produces enormous element growth
- The representation is NOT robust

**Example**: In 5 copies of W_201^+ glued by sqrt(eps), eigenvalue 0 of each block is also a Ritz value. Shifting to eigenvalue 0 causes breakdown.

### Bug 3: FP Vector Underflow (Assumption 4.4)
**Source: Dhillon, Parlett & Vomel 2005**

**The assumption**: The FP vector z computed from the twisted factorization accurately represents the eigenvector.

**How it fails**: For large Wilkinson matrices W_{2m+1}^+, eigenvector components decay exponentially away from the "peak." For the pair of eigenvalues near 0:
- Both have their peak at the center of the matrix
- Components decay as ~(1/m)^k away from center
- For m >= 50 (n >= 101), components underflow to zero in double precision
- The FP vector becomes a "bisector" (nonzero only near center), not the true eigenvector
- **Both vectors of the pair get the SAME twist index**, so they end up nearly parallel

**Root cause**: The twisted factorization at the optimal twist index r produces z(r) = 1, but z propagation via z(i) = -l(i)*z(i+1) causes underflow when |l(i)| < 1 for many consecutive i.

### Bug 4: Separate B^TB/BB^T Factorizations for Bidiag SVD
**Source: Grosser & Lang 2005**

**The problem**: The naive approach to bidiag SVD is:
1. Form T1 = B^T B (or T2 = B B^T)
2. Run MR3 on T1 to get eigenvalues (= sigma^2) and eigenvectors V
3. Recover U = B V / sigma (one-sided recovery)

**What goes wrong**: When B^T B has a cluster of eigenvalues, MR3 finds an RRR L D L^t for the cluster. But this LDL^t representation of B^T B does NOT correspond to any valid bidiagonal matrix. The singular vectors V are computed accurately (orthogonal to each other), but the recovered U = B V / sigma may NOT be orthogonal.

**Quantitative bound (Theorem 4.4)**: For separate factorizations, the absolute eigenvalue deviation is:

    |lambda_hat - lambda_check| = O(mu^2 * epsilon)

where mu = sigma_max / sigma_j. For ill-conditioned B (mu >> 1), this is much larger than the relative accuracy O(epsilon * sigma_j^2) that coupling provides.

**Geometric interpretation (stretched ellipsoid)**: Consider B = diag(sigma_1, sigma_2) with sigma_1 >> sigma_2. The set of all orthogonal V-pairs maps through B to a stretched ellipse. Only the V-pair where BV is also orthogonal gives valid singular vectors. Separate B^TB eigenvectors are "somewhere on the ellipse" but generally NOT at the orthogonal point.

### Bug 5: TGK Matrix Factorization Instability
**Source: Grosser & Lang 2005**

**The problem**: The 2n x 2n TGK (Golub-Kahan) matrix:
```
        [ 0   a_1                    ]
        [ a_1  0   b_1               ]
T_GK =  [      b_1  0   a_2          ]
        [           a_2  0   b_2     ]
        [                ...         ]
```
has eigenvalues +/- sigma_j. Its LDL^t factorization has **strongly alternating pivots**:

    d_1 = 0 (or near-zero), d_2 = -(a_1^2/d_1), d_3 = -(b_1^2/d_2), ...

For ill-conditioned B, the pivots alternate between very large and very small values, with ratios up to kappa(B)^2. This makes the LDL^t factorization unlikely to be an RRR.

**Consequence**: Running standard MR3 on T_GK directly may fail for ill-conditioned bidiagonal matrices.

---

## 4. Test Matrices That Expose These Bugs

### Wilkinson Matrix W_{2m+1}^+
**Definition**: Tridiagonal with:
- Diagonal: (m, m-1, m-2, ..., 1, 0, 1, ..., m-1, m)
- Off-diagonal: all ones

**Key properties**:
- Eigenvalues come in pairs: lambda_j and lambda_{n+1-j} differ by O(1/m!)
- The pair nearest 0 has gap ~ 2/(m! * 4^m)
- Eigenvectors peak at center, decay exponentially toward edges

**What it exposes**: Bug 3 (FP vector underflow). For m >= 50, the eigenvector components underflow, causing MR3 to return bisectors instead of eigenvectors.

**Specific sizes tested in papers**:
- W_21^+ (m=10): manageable, eigenvectors don't underflow
- W_201^+ (m=100): severe underflow, FP vectors become bisectors

### Glued Wilkinson Matrix
**Definition**: Block diagonal matrix with p copies of W_{2m+1}^+ connected by rank-1 modifications:

    T(gamma) = blkdiag(T_1, ..., T_p) + gamma * sum_{i=1}^{p-1} v_i v_i^t

where v_i = e_{last_of_block_i} + e_{first_of_block_{i+1}} (gluing adjacent blocks).

**Key parameters**:
- gamma = sqrt(machine_epsilon) ~ 1.5e-8 in double precision
- p = 5 copies of W_201^+

**What it exposes**: Bugs 1 and 2 simultaneously. The gluing creates O(p) copies of the tight Wilkinson eigenvalue clusters, connected by O(sqrt(eps)) coupling. The clusters cannot be resolved by shifting, and the RRR may not exist.

**Eigenvalue structure**: Each Wilkinson eigenvalue lambda_j splits into p copies spaced O(gamma^2/gap) apart, where gap is the original Wilkinson gap. For the tightest pair, the p copies are spaced at O(gamma^2 * m!) ~ O(eps * m!), which for m=100 means spacing >> eps * lambda but splitting requires impractical tree depth.

### T0, T1, T2 (Small diagnostic matrices)
**Source: Parlett & Dhillon 2000**

**T0** (3x3, benign):
```
T0 = tridiag([1, 2, 1], [0.5, 0.5])
```
- LDL^t: d = (1, 1.75, 6/7), l = (0.5, 2/7)
- relcond for all eigenvalues: O(1)
- Conclusion: LDL^t IS an RRR

**T1** (3x3, hard):
```
T1 = tridiag([1, 1+eps, 1], [1, 1]) where eps ~ 1e-15
```
- Eigenvalues: 1+eps+sqrt(2), 1+eps-sqrt(2), and one near 0
- LDL^t has element growth: d_3 = O(eps), l_2 = O(1/eps)
- relcond = O(1/eps) --- NOT robust
- Conclusion: LDL^t is NOT an RRR (need to shift first)

**T2** (4x4): Used to illustrate that shifted factorization T - tau*I can be an RRR even when T's own factorization is not.

### Bidiagonal Test: Cholesky of W_21^+
**Source: Grosser & Lang 2005**

Take W_21^+ and compute B such that B^T B = W_21^+ (i.e., B is the Cholesky-like bidiagonal factor). This gives a bidiagonal matrix where:
- sigma_j^2 are the eigenvalues of W_21^+
- The tight eigenvalue pairs of W_21^+ create tight singular value pairs in B
- Separate B^TB eigenvectors give poor singular vectors (Bug 4)

---

## 5. Fixes and Remedies

### Fix 1: Random Perturbation of Root Representation
**Source: Dhillon, Parlett & Vomel 2005, Section 5**

**What it fixes**: Bugs 1 and 2 (glued matrix failures)

**Method**: Before starting the representation tree, randomly perturb the initial LDL^t factorization:

    d_hat[i] = d[i] * (1 + mu_i * k * epsilon)
    l_hat[j] = l[j] * (1 + nu_j * k * epsilon)

where:
- mu_i, nu_j are random in [-1, 1]
- k ~ 4 (a small integer constant)
- epsilon = machine epsilon

**Why it works**: The perturbation breaks the exact algebraic relationships that cause Ritz value coincidences (Bug 2). It shifts eigenvalues by O(k * eps * |lambda|), which is within the unavoidable roundoff anyway. But it generically destroys the multiplicity structure that prevents RRR existence.

**Preservation**: The perturbed matrix has the same eigenvalues to machine precision (perturbation is O(eps)), so accuracy is not affected.

**Complexity**: O(n) work --- does not affect O(n^2) bound.

### Fix 2: Selective Inverse Iteration with Modified Gram-Schmidt
**Source: Dhillon, Parlett & Vomel 2005, Section 5.2**

**What it fixes**: Bug 3 (FP vector underflow) and any remaining cluster that cannot be resolved

**Method**: When the FP vector fails (detected by checking ||z||_2 or orthogonality), fall back to:
1. Solve (T - lambda*I) * x = random_rhs using bidiagonal solve
2. MGS-orthogonalize x against previously computed eigenvectors in the cluster
3. Normalize

**Complexity**: O(n * k^2) for a cluster of size k. Total O(n^2) if sum of k^2 over all clusters is O(n). This is NOT guaranteed for arbitrary matrices --- glued matrices can have a single cluster of size O(n), giving O(n^3).

**Safeguard**: Limit cluster size. If a cluster has k > sqrt(n) eigenvalues, further subdivision is attempted before falling back.

### Fix 3: Submatrix Approach
**Source: Dhillon, Parlett & Vomel 2005, Section 5.3**

**What it fixes**: Isolated failures where a single eigenvalue's representation breaks

**Method**: Identify the submatrix of T containing the problematic eigenvalue. Solve the smaller eigenproblem directly. Pad the eigenvector with zeros.

**Complexity**: O(n * k) where k is submatrix size. Only works when the eigenvalue is "isolated" to a submatrix.

### Fix 4: Rayleigh-Ritz Projection
**Source: Dhillon, Parlett & Vomel 2005, Section 5.4**

**What it fixes**: Multiple FP vector failures in the same cluster

**Method**: For a cluster of k eigenvalues where multiple FP vectors fail:
1. Compute whatever approximate vectors you can (even bisectors)
2. Form k x k projected matrix H = Z^t T Z
3. Solve the k x k eigenproblem for H
4. Rotate: V = Z * eigvecs(H)

**Complexity**: O(n * k^2). Same O(n^2) concern as Fix 2.

### Fix 5: Coupling Transformations for Bidiag SVD
**Source: Grosser & Lang 2005, Sections 4-5**

**What it fixes**: Bug 4 (separate factorizations give inaccurate singular vectors) and Bug 5 (TGK instability)

**Core idea**: Instead of working with B^TB, BB^T, or T_GK separately, use **coupled** representations that maintain the relationship between left and right singular vectors.

#### Coupling Relations (Lemma 5.1, 5.3)
Given B with diagonal a_i and superdiagonal b_i, and T_GK's LDL^t factorization with pivots d_bar_i and multipliers l_bar_i:

**From T_GK to B^TB**: The LDL^t pivots of B^T B - sigma^2 * I can be computed from T_GK's factorization via:

    d_hat_i = -d_bar_{2i-1} * d_bar_{2i}       (coupling of adjacent TGK pivots)
    l_hat_i = l_bar_{2i} * l_bar_{2i+1}          (coupling of adjacent TGK multipliers)

**From B^TB to BB^T**: Given L_1 D_1 L_1^t = B^T B - sigma^2 I, the factorization of B B^T - sigma^2 I is:

    d2_i = d1_i * (b_i / a_i)^2 * (something involving a, b, d1)
    l2_i = l1_i * (a_{i+1} / b_i)

(Exact formulas in Lemma 5.1.)

**Why coupling helps (Theorem 5.2)**: Coupled factorizations preserve RELATIVE accuracy:

    |lambda_hat - lambda_check| = O(epsilon * lambda)

versus O(mu^2 * epsilon) for separate factorizations. The improvement factor is mu^2 = (sigma_max/sigma_j)^2, which can be enormous for ill-conditioned B.

#### Algorithm 4.1 (Grosser-Lang Coupled Bidiag SVD)
1. Form T_GK from B
2. Compute initial LDL^t of T_GK
3. For each cluster of singular values:
   a. Choose shift tau
   b. Compute COUPLED factorizations of B^TB - tau^2*I AND BB^T - tau^2*I simultaneously, using the coupling relations
   c. Compute right singular vector v from B^TB factorization (twisted FP vector)
   d. Compute left singular vector u from BB^T factorization (twisted FP vector)
   e. Both u and v are computed in their natural representations, avoiding the ill-conditioned one-sided recovery u = Bv/sigma

**Advantage over one-sided**: When sigma is tiny, u = Bv/sigma amplifies errors by 1/sigma. The coupled approach computes u directly from BB^T, avoiding this amplification.

**Advantage over T_GK directly**: Avoids the alternating-pivot instability (Bug 5) by working with B^TB/BB^T representations derived from T_GK via coupling.

---

## 6. Orthogonality Guarantees and Their Limits

### What MR3 Guarantees (Theorem 15, Dhillon & Parlett 2004)
For well-separated eigenvalues (relgap >= threshold) with good representations (relcond = O(1)):

    |sin angle(z_hat, v)| <= 5n * eps + O(eps)

This gives eigenvectors accurate to O(n * eps), which implies:

    |z_i^t z_j| <= O(n * eps) for i != j

i.e., orthogonality to near machine precision.

### What MR3 Does NOT Guarantee
1. **Clustered eigenvalues**: If relgap < threshold, the representation tree must shift to create relative gaps. Success depends on the RRR existing (which can fail --- Bugs 1, 2).

2. **Large n with degenerate spectrum**: For Wilkinson-type matrices with very tight pairs, FP vectors underflow (Bug 3) and orthogonality is completely lost.

3. **Bidiag SVD**: Even if V (right singular vectors) are orthogonal, U = BV/sigma may not be orthogonal. Coupling (Fix 5) is needed.

4. **Guaranteed O(n^2)**: Orthogonality is guaranteed only when the representation tree has O(1) depth. Pathological matrices (glued Wilkinson) can force O(n) depth, breaking both the complexity and accuracy guarantees.

### The Fundamental Tradeoff
MR3 achieves O(n^2) AND orthogonality for "most" matrices, but the guarantees degrade for:
- Matrices with eigenvalue clusters that cannot be resolved by O(1) levels of shifting
- Matrices where LDL^t factorizations have inherent element growth
- Bidiagonal SVD when sigma_min/sigma_max is tiny (requires coupling)

This is why the O(n^2) bidiagonal SVD via MR3 with full accuracy remains an open problem.

---

## 7. Implications for Bidiag SVD Implementation

### What Works (Current Approach K)
- B^TB eigensolver via scipy's stemr (compiled MR3) --- fast, accurate eigenvalues/vectors
- One-sided recovery U = BV/sigma --- works when sigma is not too small relative to sigma_max
- TGK fallback when B^TB is ill-conditioned
- Chunked MGS for O(n^2) reorthogonalization

### What's Missing for a Complete Solution
1. **Coupled factorizations (Algorithm 4.1)**: Would fix the one-sided recovery problem for tiny singular values. Currently handled by zeroing unreliable sigmas, which is a workaround.

2. **NCD shifts inside eigensolver**: Would guarantee that representation tree depth stays O(1) even for pathological spectra. Not accessible from outside scipy.

3. **Random perturbation of root representation**: Could help with glued-matrix-type failures. Currently not applied (scipy's stemr may or may not do this internally).

4. **Relative-accuracy eigenvalues for B^TB**: dqds-based refinement gives eigenvalues to relative precision, which is crucial for coupling. Bisection (used by mr3_tridiag.py) only gives absolute accuracy.

### Priority Order for Implementation
1. **Coupling transformations** (highest impact): Directly addresses the #1 accuracy failure (Bug 4). Requires:
   - LDL^t factorization of T_GK (straightforward)
   - Coupling formulas to derive B^TB and BB^T factorizations (Lemma 5.1, 5.3)
   - Twisted factorization + FP vector for BOTH u and v
   - Shift selection that works for coupled representations

2. **Random perturbation** (medium impact): Simple to implement (O(n) code), helps with edge cases.

3. **NCD shifts** (highest complexity, open research): Requires modifying the core eigensolver. The Willems-Lang approach is the most promising framework but was never fully implemented.

---

## 8. Key Formulas Reference

### Sturm Count
Number of eigenvalues of T less than lambda equals the number of negative pivots in L D L^t = T - lambda*I.

### dstqds (Stationary QD with Shift)
```
d_plus[1] = d[1] - s
for i = 1 to n-1:
    l_plus[i] = (d[i] * l[i]) / d_plus[i]
    d_plus[i+1] = d[i+1] * (l[i]^2) * (d[i] / d_plus[i]) - s
```
Computes L_+ D_+ L_+^t = L D L^t - s*I.

### Twisted Factorization at Index r
```
Left sweep (i = 1 to r):   L_+ D_+ via dstqds with shift lambda
Right sweep (i = n to r+1): L_- D_- via dstqds from right with shift lambda
gamma_r = d_plus[r] + d_minus[r] - (T[r,r] - lambda)
```
Choose r to minimize |gamma_r|.

### FP Vector from Twisted Factorization
```
z[r] = 1
for i = r-1 downto 1: z[i] = -l_plus[i] * z[i+1]   (left part)
for i = r+1 to n:     z[i] = -l_minus[i-1] * z[i-1]  (right part)
normalize z
```

### Coupling: T_GK pivots to B^TB pivots
```
Given T_GK's LDL^t with pivots d_bar[1..2n] and multipliers l_bar[1..2n-1]:
d_hat[i] = -d_bar[2i-1] * d_bar[2i]        for i = 1..n
l_hat[i] = l_bar[2i] * l_bar[2i+1]          for i = 1..n-1
```
Then L_hat D_hat L_hat^t = B^T B - sigma^2 I (with appropriate shift).

### Coupling: B^TB to BB^T
```
Given L1 D1 L1^t = B^TB - sigma^2 I:
d2[i] = d1[i] * (b[i]/a[i])^2 * ... (see Lemma 5.1)
l2[i] = l1[i] * (a[i+1]/b[i])
```
Then L2 D2 L2^t = BB^T - sigma^2 I.
