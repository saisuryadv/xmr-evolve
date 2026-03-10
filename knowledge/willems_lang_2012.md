# Willems & Lang 2012: "The MR3-GK Algorithm for the Bidiagonal SVD"

**Full citation:** Paul R. Willems and Bruno Lang, "The MR3-GK Algorithm for the Bidiagonal SVD," *Electronic Transactions on Numerical Analysis (ETNA)*, Volume 39, pp. 1-21, 2012. Published online March 5, 2012.

**Core claim:** Running MR3 on the Golub-Kahan matrix T_GK(B) directly -- with just a minor modification (NCD check on shifts) -- is a valid O(n^2) solution strategy for the bidiagonal SVD. This **refutes** the long-standing opinion that MR3 on T_GK is fundamentally flawed.

---

## Table of Contents

1. [Problem Setup and Notation](#1-problem-setup-and-notation)
2. [Three Associated Tridiagonal Problems](#2-three-associated-tridiagonal-problems)
3. [The Golub-Kahan Matrix](#3-the-golub-kahan-matrix)
4. [Why Naive Black-Box MR3 on T_GK Fails](#4-why-naive-black-box-mr3-on-tgk-fails)
5. [Algorithm 4.1 -- Full Detail](#5-algorithm-41----full-detail)
6. [Definition 4.3: GK Structure](#6-definition-43-gk-structure)
7. [Lemma 4.4: Equal Norm Property](#7-lemma-44-equal-norm-property)
8. [Theorem 4.5: Correctness Proof](#8-theorem-45-correctness-proof)
9. [Definition 4.6: NCD (Nearly Constant Diagonal)](#9-definition-46-ncd-nearly-constant-diagonal)
10. [The NCD Check for Shifts](#10-the-ncd-check-for-shifts)
11. [Example 4.8: Element Growth Vulnerability](#11-example-48-element-growth-vulnerability)
12. [Block Factorizations and Element Growth](#12-block-factorizations-and-element-growth)
13. [Representations: e-representations and Z-representations](#13-representations-e-representations-and-z-representations)
14. [How MR3 Handles Clusters](#14-how-mr3-handles-clusters)
15. [Orthogonality and Residual Bounds](#15-orthogonality-and-residual-bounds)
16. [XMR-TGK vs XMR-CPL](#16-xmr-tgk-vs-xmr-cpl)
17. [Section 5 Numerical Results -- Table 5.1](#17-section-5-numerical-results----table-51)
18. [Bugs, Limitations, and Open Problems](#18-bugs-limitations-and-open-problems)
19. [Key Constants and Parameters](#19-key-constants-and-parameters)
20. [Five Requirements for MR3 Correctness](#20-five-requirements-for-mr3-correctness)
21. [Preprocessing: Splitting Strategy](#21-preprocessing-splitting-strategy)
22. [Implementation-Critical Details](#22-implementation-critical-details)

---

## 1. Problem Setup and Notation

**(Page 8, Section 3.1)**

**Input:** Upper bidiagonal matrix B in R^{n x n} with diagonal entries a_1, ..., a_n and superdiagonal entries b_1, ..., b_{n-1}:

```
B = diag(a_1, ..., a_n) + diag_{+1}(b_1, ..., b_{n-1})
```

**Goal:** Compute the full SVD:

```
B = U * Sigma * V^*
```

where U^*U = V^*V = I, Sigma = diag(sigma_1, ..., sigma_n), sigma_1 <= ... <= sigma_n.

**Quality requirements (Equations 3.2 and 3.3, Page 8):**

- Orthogonality: `max{|U^*U - I|, |V^*V - I|} = O(n * eps)`
- Residuals: `max_i{||B*v_i - u_i*sigma_i||, ||B^*u_i - v_i*sigma_i||} = O(||B|| * n * eps)`

where eps = machine epsilon = 2^{-53} ~ 1.1e-16 for IEEE double precision (page 5).

---

## 2. Three Associated Tridiagonal Problems

**(Pages 8-9, Section 3.3)**

The BSVD can be cast as TSEPs (tridiagonal symmetric eigenproblems) via the **normal equations**:

```
BB^* = U * Sigma^2 * U^*     (left normal equation)
B^*B = V * Sigma^2 * V^*     (right normal equation)
```

with explicit entries:

```
BB^* = diag(a_1^2 + b_1^2, ..., a_{n-1}^2 + b_{n-1}^2, a_n^2) + diag_{+-1}(a_2*b_1, ..., a_n*b_{n-1})
B^*B = diag(a_1^2, a_2^2 + b_1^2, ..., b_{n-1}^2) + diag_{+-1}(a_1*b_1, ..., a_{n-1}*b_{n-1})
```

**Why the black-box normal equations approach fails (Page 9):**
Computing U from BB^* and V from B^*B independently yields orthogonal U and V separately, but the residuals `||B*v_i - u_i*sigma_i||` may be O(sigma_i) for clustered singular values -- **unacceptable for large sigma_i**. The u_i and v_i don't "fit together" because they were computed independently.

The **coupling approach** (references [19,20,21,39]) runs MR3 on one normal equation and "simulates" the other using coupling relations. This works but is more complex.

---

## 3. The Golub-Kahan Matrix

**(Pages 9-10, Section 3.3.2)**

The **Golub-Kahan (GK) matrix** or **Golub-Kahan form** is defined as:

```
T_GK(B) := P_ps * [0  B; B^* 0] * P_ps^*
```

where P_ps is the **perfect shuffle permutation** on R^{2n} defined by:

```
P_ps * x = [x(n+1), x(1), x(n+2), x(2), ..., x(2n), x(n)]^*
```

This yields a symmetric tridiagonal matrix of size 2n x 2n with **zero diagonal** and the entries of B interleaved on the off-diagonals:

```
T_GK(B) = diag_{+-1}(a_1, b_1, a_2, b_2, ..., a_{n-1}, b_{n-1}, a_n)
```

**Eigenpair-to-SVD relationship (Page 10):**

```
(sigma, u, v) is a singular triplet of B with ||u|| = ||v|| = 1
    iff
(+-sigma, q) are eigenpairs of T_GK(B), where ||q|| = 1,
    q = (1/sqrt(2)) * P_ps * [u; +-v]
```

The eigenvector q interleaves u and v components:

```
q = (1/sqrt(2)) * [v(1), u(1), v(2), u(2), ..., v(n), u(n)]^*
```

**Key property:** v makes up the odd-numbered entries in q (indices 1,3,5,...) and u the even-numbered ones (indices 2,4,6,...). Extracting singular vectors:

```
[u_i; v_i] := sqrt(2) * P_ps^* * q_i
```

---

## 4. Why Naive Black-Box MR3 on T_GK Fails

**(Pages 12-14, Section 4, Example 4.1 and Example 4.2)**

### Example 4.1 (Page 12)

A bidiagonal matrix was constructed using LAPACK's DLATMS with 20 singular values:
- sigma_13 = 0.9
- sigma_14 = 1 - 10^{-7}
- sigma_15 = 1 + 10^{-7}
- sigma_16 = 1.1
- sigma_i = sigma_{i+4}/100 for i = 12, 11, ..., 1
- sigma_i = 100 * sigma_{i-4} for i = 17, ..., 20

Singular values range from 0.9 * 10^{-8} to 110.

DSTEMR was called on T_GK(B) in R^{40x40} for the upper 20 eigenpairs. **Result:**
- DSTEMR solves the TSEP correctly (left plot of Figure 4.1 shows good eigenvector orthogonality and residuals)
- But the **extracted singular vectors** are far from orthogonal (right plot)
- Small singular values are the problem: the u and v components have lost their equal-norm property
- Normalizing norms doesn't help because orthogonality is already ruined

### Root cause -- the clustering problem (Pages 12-13, Figure 4.2)

MR3 starts by shifting T_GK(B) by tau ~ -sigma_n to get a positive definite root representation L_0*D_0*L_0^*. This shift maps the eigenvalue pairs (+sigma_i, -sigma_i) to:

```
lambda_{+i} = sigma_i - tau    (positive side)
lambda_{-i} = -sigma_i - tau   (negative side)
```

For small sigma_i, the **relative distance** between lambda_{+i} and lambda_{-i} is:

```
|lambda_{+i} - lambda_{-i}| / max{|lambda_{+i}|, |lambda_{-i}|} = 2*sigma_i / (sigma_i - tau) < gaptol
```

when sigma_i is small compared to sigma_n. This means **both +sigma_i and -sigma_i get clustered together** in MR3. Further shifts cannot reliably separate them, potentially producing eigenvectors q with identical u or v components.

### Example 4.2 (Page 14)

Exact GK eigenvectors:

```
P_ps^* q_i = (1/sqrt(2)) * [u_i; v_i] = (1/2)*[1; 1; 1; -1]
P_ps^* q_j = (1/sqrt(2)) * [u_j; v_j] = (1/2)*[1; -1; 1; 1]
```

These form a cluster. A rotation G_rot in the (i,j) plane with c^2 + s^2 = 1 produces:

```
sqrt(2)*u_i = [1; c+s],  sqrt(2)*u_j = [1; s-c]
sqrt(2)*v_i = [c-s; -1], sqrt(2)*v_j = [c+s; 1]
```

Orthogonality levels: |u_i^* u_j| = |v_i^* v_j| = s^2. If s^2 is large, orthogonality is destroyed even though the GK eigenvectors q_i, q_j are perfectly orthogonal.

### The solution (Page 13-14)

The initial outside shift is unnecessary! Take M_0 := T_GK(B) directly as root. For child representations M^+ = T_GK(B) - mu, use a **special shifting routine** that exploits the zero diagonal. This handles small singular values in one step because the positive shift never smears +sigma_i and -sigma_i together.

---

## 5. Algorithm 4.1 -- Full Detail

**(Page 14)**

```
Algorithm 4.1: MR3 on the Golub-Kahan Matrix
============================================

Input:  Upper bidiagonal B in R^{n x n}, index set I_0 subset {1,...,n}
Output: Singular triplets (sigma_i, u_i, v_i), i in I_0

Step 1: Execute the MR3 algorithm for TSEP (Algorithm 2.1), but take
        M_0 := T_GK(B) as root representation in step 1, using the
        entries of B directly.
        This gives eigenpairs (sigma_i, q_i), i in I_0.

Step 2: Extract the singular vectors via
        [u_i; v_i] := sqrt(2) * P_ps^* * q_i
```

**Critical implementation details for Step 1:**

1. **Root representation:** T_GK(B) is represented by its 2n-1 entries directly (not via LDL^* factorization). This is a **heterogeneous** representation tree -- the root uses entries while children use factorizations.

2. **First-level shifts:** For computing M^+ = T_GK(B) - mu, a **special routine exploiting the zero diagonal** must be used. This is much easier than standard dtwqds because the zero diagonal allows a positive shift to be handled in one step. References: [13, 25] and remarks in [38, Sect. 8.3].

3. **Why this works for small singular values:** With T_GK(B) as root, the first shift can be a small positive value near the cluster of small +sigma_i's. The eigenvalues -sigma_i are far away on the negative side (relative distance >> gaptol), so they don't get lumped together.

4. **No initial outside shift needed:** Unlike standard MR3 which shifts outside the spectrum to get a positive definite root, Algorithm 4.1 uses T_GK(B) itself as root. The shifts to level one handle positive definiteness.

5. **General MR3 framework:** The general setup of MR3 and its proof in [35, 37] can handle this heterogeneous representation tree. All levels in the representation tree (Figure 2.1) can be handled uniformly.

---

## 6. Definition 4.3: GK Structure

**(Page 15)**

> **Definition 4.3.** A subspace S of R^{2n x 2n} with orthonormal basis (q_i)_{i in I} is said to have **GK structure** if the systems (u_i)_{i in I} and (v_i)_{i in I} of vectors extracted according to
>
> `[u_i; v_i] := sqrt(2) * P_ps^* * q_i,  i in I`
>
> are **orthonormal each**.

**Key properties:**

1. **Shift-invariant:** Any invariant subspace belonging to (at most) the first or second half of the spectrum of a GK matrix has GK structure. This is because eigenvectors are shift-invariant, and T_GK(B) - mu for suitable B is just a symmetric tridiagonal matrix of even dimension with a constant diagonal.

2. **Automatically satisfied at level one:** By Remark 4.7 (page 17), all representations on level one of the representation tree will automatically be ncd and hence have subspaces near GK structure.

3. **The crucial insight:** GK structure means that extracting u_i and v_i from q_i yields orthonormal sets. If the invariant subspaces at each node "nearly" have GK structure, then the extracted singular vectors will be nearly orthonormal.

---

## 7. Lemma 4.4: Equal Norm Property

**(Page 15)**

> **Lemma 4.4.** Let the subspace S in R^{2n x 2n} have GK structure. Then for each s in S,
>
> `sqrt(2) * s = P_ps * [s_u; s_v]  with  ||s_u|| = ||s_v||`

**Proof:** Since S has GK structure, it has an orthonormal basis (q_1, ..., q_m) with orthonormal u_i and v_i. Any s in S can be written as s = alpha_1*q_1 + ... + alpha_m*q_m, giving s_u = alpha_1*u_1 + ... + alpha_m*u_m and s_v = alpha_1*v_1 + ... + alpha_m*v_m. Since the u_i and v_j are orthonormal respectively, ||s_u||^2 = sum(alpha_i^2) = ||s_v||^2.

**Practical implication (Lemma 3.1, Page 10):** For approximate GK eigenvectors, if the angle phi from the exact eigenvector is small, then ||u'|| and ||v'|| deviate from 1/sqrt(2) by at most sin(phi) + O(sin^2(phi)). This bounds the norm deviation of extracted singular vectors.

---

## 8. Theorem 4.5: Correctness Proof

**(Pages 15-17)**

### Statement

> **Theorem 4.5 (Proof of correctness for Algorithm 4.1).** Let Algorithm 4.1 be executed such that the representation tree built by MR3 satisfies all five requirements listed in Section 2.5. Furthermore, let each node (M, I) have the property that a suitable perturbation M_GK = erp(M, xi_GK) can be found such that the subspace Q_I[M_GK] has GK structure. Finally, let resid_GK and orth_GK denote the right-hand side bounds from Theorem 2.4 and from the second inequality in Theorem 2.5, respectively. Then the computed singular triplets will satisfy:
>
> **Orthogonality:**
> ```
> max{cos(angle(u_i, u_j)), cos(angle(v_i, v_j))} <= 2*sqrt(2)*A,  i != j
> ```
>
> **Norm deviation:**
> ```
> max{| ||u_i|| - 1 |, | ||v_i|| - 1 |} <= sqrt(2)*A + O(A^2)
> ```
>
> **Residuals:**
> ```
> max{||B*v_i - u_i*sigma_i||, ||B^*u_i - v_i*sigma_i||} <= sqrt(2) * resid_GK
> ```
>
> where `A := orth_GK + C_vecs * n * xi_GK / gaptol`.

### Key proof steps (Page 16)

**Step 1 -- Residual bound:**
Using the definition of T_GK(B):

```
T_GK(B)*q_i - q_i*sigma_i = (1/sqrt(2)) * P_ps * [B*v_i - u_i*sigma_i; B^*u_i - v_i*sigma_i]
```

Therefore:

```
||B*v_i - u_i*sigma_i||^2 + ||B^*u_i - v_i*sigma_i||^2 = 2*||T_GK(B)*q_i - q_i*sigma_i||^2 <= 2*resid_GK^2
```

**Step 2 -- Orthogonality bound:**
For indices i and j, let (M, N) be the last common ancestor node in the representation tree for i in I and j in J (I intersect J = empty). The angle bound from Theorem 2.5 gives:

```
sin(angle(q_i, Q_I[M])) <= orth_GK
```

Link q_i to the nearby perturbed matrix M_GK via triangle inequality:

```
sin(angle(q_i, Q_I[M_GK])) <= orth_GK + C_vecs * n * xi_GK / gaptol = A
```

Find unit vector q in Q_I[M_GK] with sin(angle(q_i, q)) <= A. Since Q_I[M_GK] has GK structure, by Lemma 4.4:

```
sqrt(2)*q = P_ps * [u; v]  with  ||u|| = ||v|| = 1
```

Apply Lemma 3.1 to get sin(angle(u_i, u)) <= sqrt(2)*A. Similarly for j. Since Q_N[M_GK] has GK structure and I intersect J = empty, the u-spaces U_I and U_J are orthogonal (x perp y). Using Cauchy-Schwarz:

```
|u_i^* u_j| = |x^*(y+s) + r^*u_j| <= ||x||*||s|| + ||r||*||u_j||
```

This gives:

```
cos(angle(u_i, u_j)) <= ||s||/||u_j|| + ||r||/||u_i|| <= 2*sqrt(2)*A
```

**Step 3 -- Norm bounds:** From Lemma 3.1, | ||u_i|| - 1 | <= sqrt(2)*A + O(A^2).

### What makes it work

The key additional requirement beyond standard MR3 is that **representations M can be perturbed to have local invariant subspaces with GK structure**. This is captured by the NCD condition (Definition 4.6).

---

## 9. Definition 4.6: NCD (Nearly Constant Diagonal)

**(Page 17)**

> **Definition 4.6.** If for a given representation of symmetric tridiagonal M there exists an elementwise relative perturbation
>
> `M_tilde = erp(M, xi)  such that  M_tilde(i,i) = c (constant)`
>
> then we say that M has a **nearly constant diagonal (ncd)**, in short M is ncd, or more precisely M in ncd(c, xi).

**Key points:**

1. Being ncd does **not** necessarily mean all diagonal entries are about equal. For example, LDL^* can be ncd even if |d_i| >> |(LDL^*)(i,i)| for some index i (see Example 4.8).

2. **NCD implies GK structure:** Any even-dimensioned symmetric tridiagonal matrix with a constant diagonal is just a shifted Golub-Kahan matrix T_GK(B) - mu for suitable B. So trivially each subspace (within one half) has GK structure. If the diagonal is only *nearly* constant, the subspaces are *nearly* GK.

3. **Practical verification for LDL^* factorization (Page 17):** Check the condition:

```
|(LDL^*)(i,i) - const| = O(eps) * max{|d_i|, |ell_{i-1}^2 * d_{i-1}|}  for all i > 1
```

4. **Vulnerability to successive shifts:** The successively shifted descendants of a GK matrix can only violate the ncd property if there was large local element growth at some diagonal entries.

---

## 10. The NCD Check for Shifts

**(Pages 17, 19)**

### The requirement

A shift candidate tau for a cluster must produce a child representation M^+ = M - tau that is ncd. From Remark 4.7 (page 17):

> Since Theorem 4.5 needs SHIFTREL anyway, the shifts T_GK(B) - mu = M^+ to get to level one must be executed with mixed relative stability. Therefore, **all representations on level one will automatically be ncd** and as such always fulfill the extra requirement of having subspaces near GK structure, independent of element growth or relative condition numbers.

### The practical check in XMR (Page 19)

In the XMR implementation, a shift candidate must satisfy:

```
M in ncd(-mu_bar, 32*n*eps)
```

where mu_bar is the mean of the cluster, in order to be considered acceptable a priori. This is the **NCD acceptance criterion** for shifts.

### What happens if NCD is violated

If no shift candidate satisfies the NCD condition, the algorithm may need to:
- Try different shift strategies
- Use block factorizations (see Section 12) to handle element growth
- In XMR-CPL, use consistency checks with Sturm counts as a fallback

### Connection to element growth

The NCD condition is closely related to element growth. If D(2) in the LDL^* factorization has huge element growth (as in Example 4.8), the representation may still be ncd. But if shifted again, the ncd property would probably be lost. Block factorizations solve this by preventing element growth from propagating.

---

## 11. Example 4.8: Element Growth Vulnerability

**(Page 17)**

### Setup

Let alpha << 1 (e.g., alpha ~ eps ~ 1.1e-16). Consider the 2x2 bidiagonal matrix:

```
B = [1  1]
    [   alpha]
```

with singular values sigma_1 ~ alpha, sigma_2 ~ 2.

### Shifting T_GK(B) by -alpha

The GK matrix is 4x4. After shifting by -alpha:

```
T_GK(B) - (-alpha)*I = [-alpha    1              ]
                        [1        -alpha   1      ]
                        [         1       -alpha  alpha]
                        [                  alpha  -alpha]
```

### LDL^* factorization

```
= LDL^*
```

with:

```
D = diag(-alpha,  (1-alpha^2)/alpha,  -alpha*(2-alpha^2)/(1-alpha^2),  -alpha*(1/(2-alpha^2)))
```

**The critical entry: D(2) = (1-alpha^2)/alpha ~ 1/alpha ~ 1/eps ~ 10^{16}**

This is **huge local element growth** -- D(2) is of order 1/eps while the matrix entries are O(1).

### NCD property

Despite the huge D(2), the LDL^* is **still ncd** because the diagonal entries of the reconstructed tridiagonal (LDL^*)(i,i) = d_i + ell_{i-1}^2 * d_{i-1} can still be close to a constant.

### The danger

If we had to shift this LDL^* representation again (to resolve a sub-cluster), the ncd property would **probably be lost completely**. The huge D(2) would amplify perturbations and destroy the nearly-constant-diagonal structure.

### Exact numbers for reference

For alpha = eps (~ 1.1e-16):
- D(1) = -1.1e-16
- D(2) ~ 9.1e+15  (element growth factor ~ 1/eps)
- D(3) ~ -1.1e-16 * 2
- D(4) ~ -1.1e-16 * 0.5

---

## 12. Block Factorizations and Element Growth

**(Page 18, with references [35, 36])**

### The problem

Golub-Kahan matrices are **highly vulnerable to element growth** when confronted with a tiny shift. Even one level of shifting can produce D entries of order 1/eps. A second shift would lose all structure.

### The solution: block factorizations

References [35, 36] investigate a **generalization to twisted factorizations called block factorizations**:

> "In [35, 36] a generalization to twisted factorizations called block factorizations is investigated. The latter are especially suited for shifting Golub-Kahan matrices and essentially render the above concerns obsolete."

Block factorizations:
- Work with small (e.g., 2x2) matrix blocks instead of scalar entries
- Prevent element growth from propagating across blocks
- Are "especially suited for shifting Golub-Kahan matrices"
- Can be used to increase accuracy in the combined XMR-TGK + coupling approach

### Current status

Block factorizations are described in Willems' thesis [35] and the block factorization paper [36] (Willems & Lang, ETNA 38, 2011, pp. 363-400). They are mentioned as promising but the XMR-TGK implementation described in the 2012 paper uses standard (scalar) factorizations with the NCD check as the primary safeguard.

---

## 13. Representations: e-representations and Z-representations

**(Pages 18-19)**

### Background: what is a representation?

**(Definition 2.1, Page 4)** A representation M of T in R^{n x n} is a set of m <= 2n-1 scalars (the **primary data**) together with a mapping f: R^m -> R^{2n-1} that generates the entries of T.

### Twisted factorizations

**(Pages 3-4, Equation 2.1)** A symmetric tridiagonal T = N_k * G_k * N_k^* where:

```
N_k = lower-triangular with 1's on diagonal, ell_1..ell_{k-1} below diagonal,
      u_{k+1}..u_n above diagonal (from row k+1)
G_k = diag(d_1, ..., d_{k-1}, gamma_k, r_{k+1}, ..., r_n)
```

The twist element gamma_k = d_k + r_k - T(k,k).

### Standard DSTEMR representation

DSTEMR (LAPACK) uses twisted factorizations T = N_k * G_k * N_k^* represented by the nontrivial entries d_1, ..., d_{k-1}, gamma_k, r_{k+1}, ..., r_n from G_k and offdiagonal entries ell_1, ..., ell_{k-1}, u_{k+1}, ..., u_n from N_k.

### e-representation (XMR)

**(Page 18)** XMR uses the same entries from G_k, together with the n-1 **offdiagonal entries of the tridiagonal matrix T**:

```
ell_1*d_1, ..., ell_{k-1}*d_{k-1}, u_{k+1}*r_{k+1}, ..., u_n*r_n
```

These are the products T(i,i+1) = ell_i * d_i (for the L part) and T(i,i+1) = u_{i+1} * r_{i+1} (for the U part).

**Advantage:** The e-representation "provides somewhat smaller error bounds at comparable cost" (page 18). It captures the actual off-diagonal matrix entries rather than the factorization entries, which provides better relative sensitivity.

### Z-representation (XMR-TGK)

**(Pages 18-19)** For **child nodes** in XMR-TGK, Z-representations are used. These use the entries of G_k together with the n-1 quantities:

```
ell_1^2 * d_1, ..., ell_{k-1}^2 * d_{k-1}, u_{k+1}^2 * r_{k+1}, ..., u_n^2 * r_n
```

**Advantage:** Z-representations "provide even sharper error bounds, albeit at higher cost" (page 18-19). They are quadratic in the factorization off-diagonals, which means relative perturbations in Z-data cause even smaller absolute perturbations in the matrix entries.

**Motivation for XMR-TGK:** Z-representations were adopted to "cushion the effect of moderate element growth on the diagonal" (page 18). Since T_GK(B) has a zero diagonal and is vulnerable to element growth after shifting, Z-representations help maintain the ncd condition at child nodes.

### Summary table

| Representation | Primary Data | Error Bounds | Used In |
|---|---|---|---|
| Standard (DSTEMR) | G_k entries + N_k entries | Baseline | LAPACK DSTEMR |
| e-representation | G_k entries + T off-diagonals | Tighter | XMR (root/central) |
| Z-representation | G_k entries + squared off-diag products | Tightest | XMR-TGK (children) |

---

## 14. How MR3 Handles Clusters

**(Pages 2-4, Section 2.1, Algorithm 2.1)**

### The representation tree

MR3 builds a **representation tree** (Figure 2.1, page 3) where:
- **Root:** Original matrix T (or T_GK(B) in Algorithm 4.1) with representation M_0
- **Nodes:** Shifted matrices representing sub-clusters
- **Leaves:** Individual eigenpairs computed by RQI

### Cluster resolution process (Algorithm 2.1, Page 4)

1. **Initialize:** S = {(M_0, I_0, tau_bar = 0)} where I_0 is the index set of desired eigenvalues

2. **Main loop:** While S is non-empty, remove a node (M, I, tau_bar):
   - **Approximate eigenvalues** lambda_i^{loc} of M for i in I
   - **Classify** into singletons and clusters based on relative gaps >= gaptol
   - For each singleton {i}: compute eigenvector via RQI, set lambda_i = lambda_i^{loc} + tau_bar
   - For each cluster I_r:
     - Optionally refine boundary eigenvalue approximations
     - **Choose shift tau** near the cluster
     - Compute child representation M^+ = M - tau
     - Add (M^+, I_r, tau_bar + tau) to S

3. **Key parameters:**
   - **gaptol** ~ 0.001 to 0.01: threshold for classifying singletons vs clusters
   - Eigenvalues agree to ~3 leading digits when gaptol ~ 10^{-3}

### For the GK-specific case

The cluster handling is identical to standard MR3 except:
- The root is T_GK(B) represented by its entries (not LDL^*)
- First-level shifts use the special zero-diagonal-aware routine
- The NCD check (Definition 4.6) must be satisfied at each shift
- Z-representations are used for children to mitigate element growth

### Depth of the tree

The tree depth d_max affects the orthogonality bound (see Theorem 2.5). For well-separated singular values, d_max = 1 (all singletons at level 1). For highly clustered problems, d_max can be larger. The bounds grow with d_max but gaptol ensures relative gaps open up at each level.

---

## 15. Orthogonality and Residual Bounds

### Theorem 2.4: Residual norms for MR3 (Page 7)

For given index j in I_0, let d = depth(j) be the depth of the node where q = q_j was computed. Let M_0, M_1, ..., M_d be representations along the path, with shifts tau_i. Then:

```
||(M_0 - lambda^*) * q|| <= (||r^{leaf}|| + gamma * spdiam[M_0]) * (1 + beta_dagger) / (1 - beta_dagger) =: resid_{M_0}
```

where:
- lambda^* = tau_0 + ... + tau_{d-1} + lambda^{leaf}
- gamma = C_{elg} * n * (d*(alpha_down + alpha_up) + 2*(d+1)*beta_dagger^2)

### Theorem 2.5: Orthogonality bound (Page 7)

For each node (M, I) in the tree with child index set J subset I, the computed vectors q_j obey:

```
sin(angle(q_j, Q_J[M])) <= C_{vecs} * (alpha_dagger + (depth(j) - depth(M)) * (alpha_down + alpha_up)) * n / gaptol + kappa
```

where kappa = R_{gv} * n * eps + beta_dagger. For any two computed vectors q_i, q_j (i != j):

```
(1/2) * |q_i^* q_j| <= C_{vecs} * (alpha_dagger + d_max * (alpha_down + alpha_up)) * n / gaptol + kappa =: orth_{M_0}
```

### Theorem 4.5 bounds for BSVD (Page 15-16)

Combining with GK structure, the **BSVD-specific bounds** are:

```
Orthogonality: max{cos(angle(u_i, u_j)), cos(angle(v_i, v_j))} <= 2*sqrt(2)*A
Norm deviation: max{| ||u_i||-1 |, | ||v_i||-1 |} <= sqrt(2)*A + O(A^2)
Residuals:      max{||B*v_i - u_i*sigma_i||, ||B^*u_i - v_i*sigma_i||} <= sqrt(2)*resid_GK
```

where `A = orth_GK + C_{vecs} * n * xi_GK / gaptol`.

### Practical expectations (Page 7)

- C_{vecs} and C_{elg} should be of moderate size (~10)
- alpha_down, alpha_up, alpha_dagger should be O(eps)
- beta_dagger = O(n*eps)
- R_{gv} may become as large as O(1/gaptol)
- **Overall:** resid_{M_0} = O(n*eps*||M_0||/gaptol) and orth_{M_0} = O(n*eps/gaptol)
- **For BSVD:** resid scales as O(||B||*n*eps) and orth scales as O(n*eps)

---

## 16. XMR-TGK vs XMR-CPL

**(Pages 18-19)**

### XMR-TGK

- **What:** Prototype MR3 TSEP solver (XMR) adapted to BSVD by using T_GK(B) as root representation
- **Implements:** Algorithm 4.1
- **Representation strategy:**
  - Root: T_GK(B) entries directly
  - Children: Z-representations (to cushion element growth)
- **NCD check:** Shift candidates must satisfy ncd(-mu_bar, 32*n*eps)
- **Matrix size:** All internal matrices are size 2n (the GK matrix dimension)
- **Advantages:**
  - Simpler to implement and analyze than coupling
  - All levels in the representation tree handled uniformly
  - Benefits directly from improvements to the underlying MR3 algorithm
  - Better orthogonality than XMR-CPL

### XMR-CPL

- **What:** MR3 on the GK matrix using the **coupling approach** (Section 3.3.1)
- **Implements:** The normal equations coupling strategy from [19, 20, 21, 39]
- **How it works:**
  - Runs XMR on the GK matrix ("central layer")
  - Uses **coupling relations** to implicitly run MR3 simultaneously on BB^* and B^*B ("outer layers")
  - The representations in the central layer also use Z-representations with the NCD condition
  - The outer representations only need to be checked for relative condition (not the central ones)
- **Eigenvalue refinement:** Done on the side (BB^* or B^*B) that gives the better a priori bound for relative condition
- **SHIFTREL replacement:** Cannot prove SHIFTREL always holds for coupled approach; uses Sturm count consistency checks as fallback
- **Matrix size:** Outer layer matrices are size n (half the GK size)
- **Advantages:**
  - ~20-30% faster than XMR-TGK (outer matrices are half the size)
  - Uses more bisections and RQI steps but on smaller matrices
- **Disadvantages:**
  - More complex implementation
  - 24 of 19240 Synth cases failed (up to 2.04% of triplets not computed) due to no shift candidate satisfying consistency bounds
  - Cannot prove SHIFTREL in all cases

### Head-to-head comparison

| Feature | XMR-TGK | XMR-CPL |
|---|---|---|
| Matrix size | 2n | n (outer), 2n (central) |
| Speed | Baseline | 20-30% faster |
| Orthogonality (Pract AVG) | 5.35 | 10.71 |
| Orthogonality (Pract MAX) | 48.40 | 154 |
| Residual (Pract MAX) | 4.19 | 453 |
| Failures | None | 24/19240 Synth cases |
| Correctness proof | Complete | Has subtle error (ref [21]) |
| Implementation complexity | Simpler | More complex |

### Proposed hybrid (Page 20)

A third method combining both: use MR3 on T_GK(B) like Algorithm 4.1, but employ coupling relations to outsource expensive eigenvalue refinements to smaller matrices of half the size. This retains XMR-TGK's accuracy advantage at reduced cost. **Status:** "new and sounds promising" but proof in [21] contains a subtle error.

---

## 17. Section 5 Numerical Results -- Table 5.1

**(Page 20)**

### Test sets

1. **Pract** (75 cases): Bidiagonal matrices from tridiagonal problems from various applications. Dimensions up to 6245.

2. **Synth** (19240 cases): Artificially generated from tridiagonal problems. Standard types (Wilkinson, DLATMS distributions), up to dimension 100, with 2 or 3 copies glued together.

### Table 5.1 -- Exact numbers

**Orthogonality levels: max{|U^*U - I|, |V^*V - I|} / (n*eps)**

| Metric | XMR-TGK (Pract) | XMR-CPL (Pract) | XMR-TGK (Synth) | XMR-CPL (Synth) |
|--------|-----------------|-----------------|-----------------|-----------------|
| AVG    | 5.35            | 10.71           | 5.34            | 6.33            |
| MED    | 2.71            | 2.44            | 1.38            | 1.01            |
| MAX    | 48.40           | 154             | 3095            | 27729           |
| 0...10 | 81.33%          | 82.67%          | 92.59%          | 91.04%          |
| 10...100 | 14.67%        | 14.67%          | 7.04%           | 8.61%           |
| 100...200 | 2.67%        |                 | 0.12%           | 0.21%           |
| 200...500 |              |                 | 0.11%           | 0.10%           |
| 500...1000 |             |                 | 0.07%           | 0.02%           |
| 1000...1e6 |             |                 | 0.06%           | 0.03%           |

**Residual norms: max_i{||B*v_i - u_i*sigma_i||, ||B^*u_i - v_i*sigma_i||} / (||B||*n*eps)**

| Metric | XMR-TGK (Pract) | XMR-CPL (Pract) | XMR-TGK (Synth) | XMR-CPL (Synth) |
|--------|-----------------|-----------------|-----------------|-----------------|
| AVG    | 0.35            | 15.78           | 0.45            | 3.14            |
| MED    | 0.07            | 1.37            | 0.13            | 0.72            |
| MAX    | 4.19            | 453             | 118             | 6873            |
| 0...1  | 92.00%          | 34.67%          | 84.96%          | 57.45%          |
| 1...10 | 8.00%           | 50.67%          | 15.03%          | 35.50%          |
| 10...100 |               | 8.00%           |                 | 7.00%           |
| >100   |                | 6.67%           | 0.01%           | 0.06%           |

### Key observations from the paper (Page 19-20)

1. **XMR-TGK works "amazingly well"** -- better orthogonality than LAPACK DSTEMR on B^*B alone, and residuals not much worse than XMR.

2. **XMR-CPL works well on Pract** but has "undeniable problems" on Synth: 24 of 19240 cases had failures (up to 2.04% of singular triplets not computed) due to no shift candidate satisfying consistency bounds for SHIFTREL.

3. **Example 4.1 matrix:** Both XMR-TGK and XMR-CPL solve it with orthogonality of 1.15*n*eps and BSVD residual norms of 0.68*||B||*n*eps. The computed vectors are identical for both methods. By contrast, black-box MR3 on T_GK gives orthogonality ~80*n*eps.

4. **Speed:** XMR-CPL is expected to be about 20-30% faster than XMR-TGK because the outer layer matrices in the coupling approach are size n vs 2n.

---

## 18. Bugs, Limitations, and Open Problems

### Acknowledged limitations

1. **Element growth vulnerability (Page 17-18):** GK matrices are "highly vulnerable to element growth when confronted with a tiny shift." Example 4.8 shows D(2) ~ 1/eps after just one shift. A second shift would destroy the ncd property.

2. **Block factorizations not yet integrated (Page 18):** Block factorizations [35, 36] solve the element growth problem but are only mentioned as a promising direction, not implemented in the tested code.

3. **XMR-CPL failures (Page 19):** 24 of 19240 Synth cases had clusters where no shift candidate satisfied the eigenvalue consistency bounds needed to replace SHIFTREL. These are flagged as "not computed" rather than errors.

4. **Subtle proof error in coupling (Page 20):** The combined approach (XMR-TGK + coupling for efficiency) relies on Theorem 5.2 in [21], "but its proof contains a subtle error."

5. **XMR-TGK worst case orthogonality:** MAX of 3095 on Synth (page 20, Table 5.1). This is ~3095*n*eps, which for large n could be concerning.

6. **XMR-TGK worst case residual:** MAX of 118 on Synth (page 20). This is 118*||B||*n*eps, exceeding the ideal O(1) level.

7. **No way to test for GK structure directly (Page 17):** "At the moment we do not see a way to specifically test for this property." The NCD condition is used as a sufficient proxy.

8. **No guarantee that NCD is preserved at deeper levels:** If element growth is large at level 1, level 2 shifts may violate NCD even if level 1 was ncd.

### What the paper does NOT solve

1. **Production-quality implementation:** XMR-TGK is a "prototype implementation" (page 18), not a production solver.

2. **The speed question:** XMR-TGK operates on 2n-size matrices, making it inherently slower per operation than approaches working with n-size normal equations. The 20-30% speed advantage of XMR-CPL motivates the hybrid approach.

3. **Relative accuracy of singular values:** The paper notes (page 8) that recomputing singular values via dqds after computing eigenvectors achieves relative accuracy. This step is not part of Algorithm 4.1 itself.

4. **Subset computation efficiency:** While MR3 can compute subsets, the paper notes that getting "a consistent mapping of triplet indices between the blocks" after splitting "is not entirely trivial" (page 12).

---

## 19. Key Constants and Parameters

| Symbol | Meaning | Typical Value | Reference |
|--------|---------|---------------|-----------|
| eps (epsilon_o) | Machine epsilon (IEEE double) | 2^{-53} ~ 1.1e-16 | Page 5 |
| gaptol | Relative gap tolerance for singleton/cluster classification | 0.001 to 0.01 | Page 6 |
| C_vecs | Constant for RRR requirement (eigenvector sensitivity) | ~10 | Page 7 |
| C_elg | Constant for element growth requirement | ~10 | Page 7 |
| alpha_down, alpha_up | Shift relation perturbation bounds | O(eps) | Page 7 |
| beta_dagger | Vector perturbation bound | O(n*eps) | Page 7 |
| R_gv | Rayleigh quotient iteration constant | up to O(1/gaptol) | Page 7 |
| xi_GK | NCD perturbation bound | ~32*n*eps (in XMR) | Page 19 |
| d_max | Maximum depth of representation tree | Problem-dependent | Page 7 |

---

## 20. Five Requirements for MR3 Correctness

**(Page 6-7, Section 2.5)**

These five requirements, from [37], together guarantee correctness of Algorithm 2.1:

### Requirement RRR (Relatively Robust Representations)

There exists constant C_vecs such that for any perturbation M_tilde = erp(M, alpha) at node (M, I), the effect on eigenvectors is bounded:

```
sin(angle(Q_J[M], Q_J[M_tilde])) <= C_vecs * n * alpha / relgap_M(J)
```

for all J in {I, I_1, ..., I_r} with |J| < n.

This implies singleton eigenvalues and cluster boundary eigenvalues are relatively robust.

### Requirement ELG (Conditional Element Growth)

There exists constant C_elg such that element growth is bounded:

```
||M_tilde - M|| <= spdiam[M_0]
||(M_tilde - M)*q_i|| <= C_elg * n * alpha * spdiam[M_0]  for each i in I
```

Large element growth is permissible overall, but only where the local eigenvectors have tiny entries.

### Requirement RELGAPS (Relative Gaps)

The classification into child index sets ensures:

```
relgap_M(I_r) >= gaptol  (if |I_r| < n)
```

for r = 1, ..., m. Cannot be fulfilled if relgap_M(I) < gaptol at the node's creation -- must be checked during shift evaluation.

### Requirement SHIFTREL (Shift Relation)

There exist constants alpha_down, alpha_up such that for every child node with matrix H computed from parent M using shift tau:

```
M_tilde = erp(M, alpha_down)  and  H_tilde = erp(H, alpha_up)
```

exist satisfying the exact shift relation M_tilde - tau = H_tilde.

This connects nodes in the tree and ensures shifts are computed in a mixed relatively stable way. Fulfilled when using twisted factorizations with qd-type transformations.

### Requirement GETVEC (Computation of Eigenvectors)

There exist constants alpha_dagger, beta_dagger, and R_gv such that the computed eigenvector q at node (M, I) has small residual norm:

```
||r^{leaf}|| := ||(M_tilde - lambda^{leaf})*q|| / ||q|| <= R_gv * n * eps * gap_{M_tilde}({i}; lambda^{leaf})
```

This captures that RQI must produce small residuals. The keys are qd-type transformations to compute twisted factorizations M - lambda = N_k * G_k * N_k^* and then solving N_k * G_k * N_k^* * q = gamma_k * e_k.

---

## 21. Preprocessing: Splitting Strategy

**(Pages 11-12, Section 3.4)**

### Scaling and sign normalization

1. If input is lower bidiagonal, work with B^* and swap U/V roles
2. Multiply by diagonal signature matrices to make all entries nonneg
3. Scale to get the largest elements into proper range

### Splitting criterion (Equation 3.5)

```
n * eps * ||B|| < min{a_i, b_i}
```

**Issue:** Setting tiny entries to zero to achieve (3.5) impedes computing singular values to high relative accuracy.

### Two-phase splitting strategy

The paper proposes a 2-phase split:

**Phase 1 -- Relative split:** Split B as much as possible **without spoiling relative accuracy**. Uses splitting criteria from the zero-shift QR algorithm [4, p. 18] and Li [28, 32] which retain relative accuracy. This gives blocks B_rs^{(1)}, ..., B_rs^{(N)}.

**Phase 2 -- Absolute split:** Split each B_rs^{(i)} further aggressively into blocks B_as^{(i,1)}, ..., B_as^{(i,n_i)} using the criterion (3.5).

Then:
3. Solve BSVD for each block in the absolute split independently
4. Use bisection to refine computed singular values of each B_as^{(i,j)} to high relative accuracy with respect to parent block B_rs^{(i)}

**If dqds is used to precompute singular values:** Steps 1 and 4 can be skipped entirely since dqds already provides relative accuracy and the absolute-split singular values are discarded anyway.

### Handling zero diagonal entries

If a diagonal element a_i = 0, B is singular. Deflation via one sweep of implicit zero-shift QR method yields B' with b'_{i-1} = b'_{n-1} = a'_n = 0. This reveals one zero singular value and splits B into three parts.

---

## 22. Implementation-Critical Details

### The perfect shuffle permutation P_ps

**(Page 9)**

Maps x in R^{2n} to:
```
P_ps * x = [x(n+1), x(1), x(n+2), x(2), ..., x(2n), x(n)]^*
```

Inverse:
```
P_ps^* * x = [x(2), x(4), ..., x(2n), x(1), x(3), ..., x(2n-1)]^*
```

### Extracting singular vectors from GK eigenvectors (Equation 3.4)

```
q = (1/sqrt(2)) * [v(1), u(1), v(2), u(2), ..., v(n), u(n)]^*
```

So:
- v(j) = sqrt(2) * q(2j-1)  for j = 1, ..., n  (odd indices)
- u(j) = sqrt(2) * q(2j)    for j = 1, ..., n  (even indices)

### Special first-level shift routine

For M^+ = T_GK(B) - mu where T_GK(B) has zero diagonal:
- Much easier than standard dtwqds
- See references [13, 25] and [38, Sect. 8.3]
- Handles positive shifts to isolate positive singular values from negative ones
- Small singular values can be handled by a (positive) shift without danger of "spoiling them by unwanted contributions from the negative counterparts" (page 14)

### XMR differences from DSTEMR (Page 18)

1. Uses e-representations instead of standard twisted factorization data
2. Even if RRR and moderate ELG cannot be guaranteed before performing a shift, sufficient a priori criteria are available (improved in XMR)
3. Several modifications for robustness: interplay of RQI and bisection, bisection strategy itself

### NCD acceptance criterion in XMR (Page 19)

```
shift candidate must satisfy: ncd(-mu_bar, 32*n*eps)
```

where mu_bar is the mean/center of the cluster being resolved.

### Normalization (Page 17)

Theorem 4.5 shows it does not matter whether singular vectors are extracted by multiplying q subvectors by sqrt(2) or by normalizing explicitly. Both approaches yield the same bounds.

### Lemma 3.1: Angle-to-component bound (Page 10)

For non-orthogonal unit vectors q, q' with conforming partition into u,v components:

```
max{||u||*sin(phi_u), ||v||*sin(phi_v)} <= sin(phi)
max{| ||u'||-||u|| |, | ||v'||-||v|| |} <= (sin(phi) + (1-cos(phi))) / cos(phi)
```

where phi_u = angle(u,u'), phi_v = angle(v,v'), phi = angle(q,q'). For small phi, this means component norms deviate from 1/sqrt(2) by basically sin(phi) + O(sin^2(phi)).

---

## Appendix: Key References from the Paper

- **[4] Demmel & Kahan 1990:** Accurate singular values of bidiagonal matrices. Foundation for relative accuracy and dqds.
- **[6] Dhillon 1997:** The O(n^2) MR3 algorithm (PhD thesis).
- **[7] Dhillon & Parlett 2004a:** Multiple representations for orthogonal eigenvectors.
- **[8] Dhillon & Parlett 2004b:** Orthogonal eigenvectors and relative gaps.
- **[13] Fernando 1998:** Accurately counting singular values via skew-symmetric tridiagonals. (Key for the special first-level shift routine.)
- **[19] Grosser 2001:** The original O(n^2) BSVD via coupling (PhD thesis).
- **[20] Grosser & Lang 2003:** The coupling approach paper.
- **[25] Kahan 1966:** Accurate eigenvalues of a symmetric tridiagonal matrix. (Used for first-level shifts.)
- **[35] Willems 2010:** PhD thesis with the comprehensive MR3 framework, XMR implementation, block factorizations.
- **[36] Willems & Lang 2011a:** Block factorizations and qd-type transformations for MR3.
- **[37] Willems & Lang 2011b:** Framework for MR3: theory and implementation.
- **[38] Willems & Lang 2011c:** Twisted factorizations and qd-type transformations -- new representations and analysis.
- **[39] Willems, Lang & Vomel 2007:** Computing the bidiagonal SVD using multiple relatively robust representations. (Earlier version of coupling approach.)
