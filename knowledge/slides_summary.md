# Slides Summary: MR^3 for Bidiagonal SVD

## Source Files
- `BeBOP Talk: Revisiting MR^3 for the Bidiagonal SVD.pdf` (88 slides) -- Ryan Schneider, UC Berkeley (ryan.schneider@berkeley.edu, Evans 895)
- `Slides on Willems Thesis.pdf` (16 slides) -- Beresford Parlett, UC Berkeley, September 2010, LAPACK seminar
- `Holy Grail.png` -- joke image used on final slide of BeBOP talk

---

## 1. DBDSVDX Bug (from BeBOP Talk)

**There is no explicit DBDSVDX bug described in these slides.** The slides discuss the general failure mode of using MR^3 as a black box for the bidiagonal SVD, and present a specific failure example from [Marques, Demmel, & Vasconcelos 2020].

### The Failure Example (Slide 87, "Failure Example")

**Matrix:** T is a 10 x 10 tridiagonal matrix with eigenvalues:

```
lambda_j = (1/sqrt(eps))^(-(j-1)/9)
```

where eps = machine precision. This creates a geometrically graded spectrum spanning the full range from 1 to 1/sqrt(eps).

**Result:** Computed eigenvectors V satisfy only:

```
||I - V^T V||_2 / (n * eps) = 3.95 x 10^4
```

This is a catastrophic orthogonality failure (should be O(1), got ~40,000 times worse).

**Attribution:** [Marques, Demmel, & Vasconcelos 2020] -- this is described as a failure of the existing MR^3-based bidiagonal SVD approach (xBDSCR/xBDSVDX), not a specific code bug.

### How to Construct a Test Matrix Triggering This Failure

The matrix is constructed with eigenvalues lambda_j = (1/sqrt(eps))^(-(j-1)/9) for j=1..10. For double precision (eps ~ 2.22e-16):
- lambda_1 = 1 (j=1: exponent 0)
- lambda_10 = (1/sqrt(eps))^(-1) = sqrt(eps) ~ 1.49e-8
- Intermediate values geometrically graded between 1 and sqrt(eps)

This creates a tridiagonal T with condition number ~1/sqrt(eps) ~ 6.7e7, where eigenvalues span roughly 8 orders of magnitude. This wide relative spread causes the MR^3 representation tree to require many shifts, and the accumulated errors from squaring (B^TB or BB^T approach) or from the TGK zero-structure (Golub-Kahan approach) lead to loss of orthogonality.

---

## 2. The Core Bug/Issue: Pitfalls of Using MR^3 as a Black Box (Slides 77-80)

This is the central technical issue documented in the talk. There are two options for computing the bidiagonal SVD via MR^3, and **neither gives both residual accuracy AND orthogonality simultaneously**:

### Option 1: Apply MR^3 to B^TB (for V) and BB^T (for U)
- **Orthogonality: PASS** -- |v_i^T v_j| = |u_i^T u_j| = O(eps*n) for i != j
- **Residual: FAIL** -- ||Bv_i - sigma_i * u_i||_2 = O(eps*n*||B||_2)... but actually FAILS because U and V are computed independently from different matrices, so they are not properly coupled.

### Option 2: Form Golub-Kahan matrix T_GK, apply MR^3 to T_GK, extract U and V
- **Residual: PASS** -- ||Bv_i - sigma_i * u_i||_2 = O(eps*n*||B||_2)
- **Orthogonality: FAIL** -- |v_i^T v_j| = |u_i^T u_j| = O(eps*n) FAILS because the eigenvalue structure of T_GK (with +/- sigma pairs) creates challenges for MR^3's representation tree and shift strategy.

**This is the fundamental dilemma**: No black-box approach to using MR^3 for the bidiagonal SVD achieves both properties simultaneously. This is what makes the problem hard and still open.

---

## 3. The "Holy Grail" (PNG and Slide 88)

The Holy Grail PNG is a humorous parody of the "Monty Python and the Holy Grail" DVD cover, photoshopped to read **"Beresford Parlett and the Holy Grail"** with the tagline **"Pure Algorithmic Genius"** at the bottom.

**What it represents:** The "Holy Grail" is achieving an O(n^2) bidiagonal SVD algorithm via MR^3 that simultaneously delivers:
1. Small residual: ||Bv_i - sigma_i * u_i||_2 = O(eps*n*||B||_2)
2. Orthogonal left singular vectors: |u_i^T u_j| = O(eps*n)
3. Orthogonal right singular vectors: |v_i^T v_j| = O(eps*n)

This has been Beresford Parlett's long-standing research quest -- making MR^3 work reliably for the bidiagonal SVD, which remains an open problem despite decades of effort by Parlett, Dhillon, Willems, Lang, Grosser, Marques, Demmel, Vasconcelos, and others.

---

## 4. Algorithms Compared in the BeBOP Talk

The talk compares these approaches for solving the bidiagonal SVD:

### Tridiagonal SEP solvers (for the eigenvalue step):
1. **QR** (implicit QR iteration)
2. **Divide-and-conquer**
3. **Bisection + Inverse Iteration**
   - Pros: simple, can compute any subset, O(n^2) eigenvalues
   - Cons: explicit re-orthogonalization may be necessary (O(n^3) worst case)
4. **MR^3 (MRRR)** -- the main subject
   - O(kn) for k eigenpairs, no explicit orthogonalization
   - History: Dhillon 1997 thesis -> Parlett & Dhillon 2000 -> Dhillon & Parlett 2004 -> LAPACK xSTEGR 2005

### Bidiagonal SVD approaches compared:
1. **Option 1**: Apply MR^3 separately to B^TB and BB^T
   - Gets orthogonality, loses residual coupling
2. **Option 2**: Apply MR^3 to the Golub-Kahan (TGK) matrix
   - Gets residual, loses orthogonality
3. **xBDSCR** (Grosser & Lang 2003, LAPACK working note 166, 2005)
   - NCD-aware approach, but the [Marques, Demmel, Vasconcelos 2020] failure example shows it can fail catastrophically

### Three reductions from bidiagonal SVD to tridiagonal SEP (Slide 76):
1. B^TB = V Sigma^2 V^T (gives V, tridiagonal)
2. BB^T = U Sigma^2 U^T (gives U, tridiagonal)
3. TGK = [0 B; B^T 0] with perfect shuffle -> tridiagonal with +/- sigma eigenvalues (gives both U and V)

---

## 5. Test Results Shown

### Failure Example (Slide 87)
- **Matrix**: 10x10 tridiagonal with geometrically graded eigenvalues lambda_j = (1/sqrt(eps))^(-(j-1)/9)
- **Orthogonality failure**: ||I - V^T V||_2 / (n*eps) = 3.95 x 10^4
- **Source**: [Marques, Demmel, & Vasconcelos 2020]
- This demonstrates that existing implementations (xBDSCR and presumably xBDSVDX which is based on similar ideas) can fail on matrices with wide eigenvalue spread.

### History Timeline (Slides 81-85)
Key timeline entries in red (SVD-specific):
- 2003: Grosser & Lang, "An O(n^2) Algorithm for the Bidiagonal SVD" (LAA)
- 2005: xBDSCR proposed for LAPACK (working note 166)
- 2012-13: Theory revisited by Willems & Lang (+ new implementation)

---

## 6. Implementation Insights from Willems Thesis Slides

### Representations of T (Slide 4, "Representations of T (without twists)")
Four representations listed:
1. **(T)** Matrix entries {c_i, e_i}
2. **(N)** {d_i, l_i}, T = LDL*
   - Pro: defines tiny eigenvalues to high relative accuracy
   - Con: does not always exist. Element growth.
3. **(e)** {d_i, e_i}, e_i = l_i * d_i
4. **(Z)** {d_i, lld_i}, lld_i = d_i * l_i^2 = Schur complements
   - c_{i+1} = d_{i+1} + lld_i
   - Can convert between representations with **no adds or subtracts**
   - Square roots needed for (Z)

### Comparison of Representations (Slide 5)
- **Accuracy**: (Z) is best. Max error is 3 ulps for (Z) versus 4 for (e).
- **Speed**: (e) is fastest when properly optimized. (N) almost as good.
- **Conclusion**: Always use (Z) but switch to (N) or (e) when a node contains a singleton.

### Error Analysis Notation (Slides 2-3)
- Willems introduces epsilon^[k](n) := n*eps / (1 - k*n*eps)
- This is cleaner than Higham's gamma_n = n*eps/(1-n*eps) because:
  - (1 + eps^[k](n))^s = 1 + eps^[k/s](sn) for 0 < s < 1
  - (1 + eps^[k](n))^{-1} = 1 + eps^[k+1](n)

### Shift of Origin (Slide 7)
Key recurrence for shifted factorization:
```
L_+ D_+ (L_+)^T = LDL^T - tau*I
s_i = d^+_i - d_i  (where s_i - tau)
s_{i+1} = e_i^2 (s_i - tau) / (d_i * d^+_i)
```
This is "an important recurrence" -- the basis of the dstqds/dtwqds transforms.

### Error Analysis of dstqds (Slide 8)
Algorithm (unblocked, e-rep):
1. d_i^+ = d_i + (s_i - tau)
2. e_i^+ = e_i
3. s_{i+1} = e_i^2 (s_i - tau) / (d_i * d_i^+)

"Keep the computed s_i sacred."

Error bound: e_tilde_i = e_i [1 + eps^[4](3)], using Willems notation.

### Mixed Relative Error Analysis for dtwqds (Slide 9)
- Diagram showing computed vs exact factorizations
- N_k G_k N_k* --dtwqds/computed--> N_t^+ G_t^+ (N_t^+)*
- Exact relationship: -tau X where X = I + delta_k e_k e_k* + delta_t e_t e_t*, with delta_k, delta_t = O(eps*1)

### Perturbing the Shift (Slide 10)
- Outer perturbations: NGN* -> DNGN*D, D approx I
- Inner perturbations: NGN* -> NDGDN*, D approx I
- Ostrowski's Theorem: |lambda(FAF*) - lambda(A)| <= |lambda(A)| * ||FF* - I||
- **Outer perturbations make tiny relative changes** -- this is key for MR^3

### Preserving Tridiagonal Form (Slide 11)
- Allow 2x2 blocks in D of LDL^T
- Example: 4x4 case with L having entries k_1, l_2, l_3
- Key issue: "Cannot perturb k_1 and l_2 independently and preserve tridiagonal form"
- This is a complication when the matrix has 2x2 blocks (i.e., near-degenerate eigenvalue pairs)

### Blocks in LDL* (Slide 12)
- Example: T = GK - alpha*I, alpha = O(eps)
- Shows concrete 4x4 case with L, D having 2x2 block structure
- Willems representation includes: D, e, and Omega = {i where a 2x2 block ends}
- Secondary data: k_i, l_j, Delta_i = d_i*c_{i+1} - e_i^2, inv_D(i)
- **The s_{i+1} computation becomes "much more complicated (9 cases)"**

### Table of 9 Cases for s_{i+1} (Slide 13)
Cases S1-S9 depending on whether i and i+1 are in Omega/Omega^+ (block boundaries):
- S1: standard case, s_{i+1} = e_i^2/d_i - e_i^2/d_i^+
- S2: s_{i+1} = 0
- S3-S9: various formulas involving Delta_{i-1}, d_{i-1}^+, etc.

### Error Bounds Table (Slide 14)
For concrete parameters R=3, K_square=1/8:
- d_i -> d_tilde_i: 1 ulp (outside block), 2 ulps (outside), 5 ulps (inside block)
- c_i -> c_tilde_i: 4 ulps
- e_i -> e_tilde_i: 3 ulps (outside), 10 ulps (inside)
- tau -> tau_tilde: 0 ulps (outside), 12 ulps (inside)

"Only first-order bounds shown, i.e., entry p stands for p*eps_0 + O(eps_0^2)"

### Recap on MRRR (Slide 16)
Key points summarized:
- Users want orthogonality. Constraint: No Gram-Schmidt (distributive computing).
- Need ||Tz - z*lambda|| = O(n*eps)|lambda|. Need relgap(lambda) >= 10^{-3}.
- Make relgaps large by shifting origin (to clusters). Hence Multiple Representations.
- Organize computation in a Representation Tree.
- Twisted factors permit residual norms proportional to |lambda|.
  - Solve N_k G_k N_k^T z = e_k gamma_k. N_k^T z = e_k.
  - "It is e_k that yields z by products only."
  - Can check |gamma_k| = O(n*eps)|lambda|.
- **Differential qd algorithms allow for roundoff between representations.**
- **Eigenvectors invariant under exact shifts.**

---

## 7. Key Open Questions (Slide 87, "Open questions")

1. Can we finally get MR^3 working for the bidiagonal SVD?
   - Debug existing implementations
   - Revisit the theory
   - Randomize?
2. Other potential applications of MR^3:
   - Preconditioned Jacobi's method for the bidiagonal SVD?

---

## 8. Key Historical References

| Year | Work | Slide |
|------|------|-------|
| 1997 | Dhillon PhD thesis: O(n^2) algorithm for symmetric tridiagonal | 24-28 |
| 2000 | Parlett & Dhillon: Relatively Robust Representations (LAA) | 25-28 |
| 2003 | Grosser & Lang: O(n^2) Algorithm for Bidiagonal SVD (LAA) | 81 |
| 2004 | Dhillon & Parlett: Orthogonal Eigenvectors and Relative Gaps (SIMAX) | 26-28 |
| 2005 | MR^3 added to LAPACK as xSTEGR (working note 162) | 27-28 |
| 2005 | xBDSCR proposed for LAPACK (working note 166) | 82 |
| 2012-13 | Willems & Lang: theory revisited + new implementation | 84 |
| 2020 | Marques, Demmel, & Vasconcelos: failure example | 85-86 |

---

## 9. Summary of Critical Findings for Our Project

1. **The fundamental dilemma**: Option 1 (B^TB/BB^T) gets orthogonality but loses residual coupling. Option 2 (TGK) gets residual but loses orthogonality. No black-box MR^3 approach achieves both. This is exactly what we observe in our Approach K.

2. **The failure matrix**: lambda_j = (1/sqrt(eps))^(-(j-1)/9) for a 10x10 matrix gives orthogonality error of ~40,000x machine precision. This is a wide-condition-number graded matrix. We should test this pattern.

3. **Willems' 9-case s_{i+1} formula** is the key to NCD-aware MR^3 for the bidiagonal SVD. The standard dstqds/dtwqds must be extended to handle 2x2 blocks in the LDL^T factorization, which arise from the +/- eigenvalue structure of the Golub-Kahan matrix.

4. **The problem remains open** as of the date of these slides. The BeBOP talk ends with "Can we finally get MR^3 working in this setting?" as an open question.

5. **Representation choice matters**: Use (Z) representation (Schur complements) for best accuracy (3 ulps vs 4), but switch to (N) or (e) for singletons.

6. **The "Holy Grail"** is explicitly identified as Beresford Parlett's quest for a reliable O(n^2) bidiagonal SVD via MR^3. The joke image confirms this is viewed as a legendary, possibly unattainable goal in the numerical linear algebra community.
