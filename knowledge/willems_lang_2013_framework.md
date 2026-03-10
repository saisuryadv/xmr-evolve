# Willems & Lang 2013: A Framework for the MR3 Algorithm -- Theory and Implementation

**Citation:** Paul R. Willems and Bruno Lang, "A Framework for the MR3 Algorithm: Theory and Implementation," SIAM J. Sci. Comput., Vol. 35, No. 2, pp. A740-A766, 2013. DOI: 10.1137/110834020

**Source PDF:** `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Papers/2013 - A Framework for the MR3 Algorithm Theory and Implementation - Willems and Lang.pdf`

---

## 1. How This Paper Extends the 2012 Paper

The 2012 Willems-Lang paper [39] focused specifically on **twisted factorizations and qd-type transformations** -- the low-level building blocks for computing representations and shifting. This 2013 paper operates at a higher level of abstraction:

1. **Modular, representation-agnostic proof of correctness.** The original MR3 proofs by Dhillon and Parlett [8,9] were tightly intertwined with LDL* bidiagonal factorizations and twisted factorizations. This paper **disentangles** the core recursive algorithm from the specific decomposition type. The proof works for *any* representation that satisfies five abstract requirements.

2. **Supports block factorizations.** The 2012 paper introduced twisted block factorizations; this paper shows they fit cleanly into the framework. Block LDL* factorizations (with 1x1 or 2x2 diagonal blocks) were made viable for MR3 by the mixed relatively stable shifting algorithm in [38,39].

3. **Streamlined proof structure.** Rather than the two-paper proof spread across [8] and [9], this paper gives a complete, self-contained correctness proof in one place, organized around five clearly stated requirements.

4. **New implementation XMR.** Describes a complete new MR3 implementation that uses insights from the framework, particularly twisted block factorizations, and compares it against LAPACK's DSTEMR.

5. **Practical robustness improvements.** Documents specific bugs/weaknesses in DSTEMR discovered through exhaustive testing, and shows how the framework's requirements guide fixes.

---

## 2. Implementation Details (XMR vs DSTEMR)

### 2.1 DSTEMR (LAPACK 3.2.2) Design

- Uses plain LDL* factorizations only (no blocks).
- Shift selection (line 12 of Algorithm MR3): refines eigenvalues at cluster borders to full precision, then tries to shift close to one border, backing off if initial tries fail.
- Evaluates shifts by checking element growth only: condition (5.3) with C_elg = 8. If no candidate passes, takes the one with minimal element growth.
- Envelope information: uses a rough approximation -- the "rightmost" eigenvector of the cluster serves as envelope on one side, vice versa on the other. Only used when plain element growth check passes.
- Bisection: one interval per eigenvalue, refined independently. Wasteful for close eigenvalues.
- Child eigenvalue bounds initialized as: [lambda_i^+] := [lambda_i] - tau (shifted parent bounds).

### 2.2 XMR Design Principles

The guiding principle was **reliability** -- accumulate as much information as possible about shift candidates before committing:

**Key differences from DSTEMR:**

1. **Twisted block factorizations** at the nodes. Block LDL* with 1x1 or 2x2 diagonal blocks D_i. The twist pivot gamma_k is always 1x1. Twist index chosen on case-by-case basis to minimize both element growth (5.3) and relative robustness (5.4) simultaneously. Eigenvectors still computed using plain (unblocked) twisted factorizations.

2. **A priori relative robustness checking** via condition (5.4): the envelope-based bound ||Upsilon N* s|| ||Upsilon^{-1} N^{-1} s|| <= maxrc (parameter of moderate size, e.g., 10) is checked in addition to element growth. This catches cases where element growth alone is misleading.

3. **Restructured bisection** incorporating shift selection. Once a shift is deemed worthy (passes a priori tests or fails less badly than alternatives), the shifted representation is computed and initial child bounds are set as:
   ```
   [lambda_i^+] := ([lambda_i](1 +/- alpha) - tau)(1 +/- alpha)
   ```
   where alpha = 100*n*eps. Inflated bounds are verified with two Sturm counts each. If bounds are inconsistent, the candidate is immediately discarded.

4. **Cluster substructure exploitation.** XMR refines eigenvalues within clusters to detect interior gaps exceeding sqrt(eps). Envelopes are computed per subgroup using techniques from [35], giving useful envelope information in many cases where DSTEMR cannot.

5. **Inside shifts.** XMR can shift to interior gaps within a cluster (near a boundary, not at the midpoint), keeping representations nearly singular. Deactivated in released version as all test cases solvable without it.

6. **Smarter RQI/bisection tolerances.** Singletons switch from bisection to RQI earlier. Acceptance thresholds in RQI relaxed compared to DSTEMR (which was too strict, causing unnecessary fallback bisection).

7. **Underflow handling.** XMR is completely underflow-safe and IEEE 754 compliant -- no NaN or exception dependence. A specific DSTEMR bug was identified: underflow in qd-transformations causes DSTEMR to return a zero vector without error, leading to a wrong twist index and a vector with huge norm.

---

## 3. Representation Types

### 3.1 Definition of a Representation

A **representation** is a set of m <= 2n-1 scalars (the "data") together with a mapping f: R^m -> R^{2n-1} that defines the entries of T. A **relatively robust representation (RRR)** is one where small relative changes in the data (bounded by xi) cause only proportionally small relative changes in eigenvalues and eigenvectors.

### 3.2 Entrywise Representation

The obvious one: store c_1,...,c_n, e_1,...,e_{n-1} directly. The representational mapping is the identity. This will rarely be an RRR, at least not for the small eigenvalues.

### 3.3 LDL* (Classical Bidiagonal) Factorizations

```
T = LDL*    where L = unit lower bidiagonal, D = diagonal
T = URU*    where U = unit upper bidiagonal, R = diagonal
```
- L = diag(1,...,1) + diag_{-1}(l_1,...,l_{n-1}), D = diag(d_1,...,d_n)
- U = diag(1,...,1) + diag_{+1}(u_1,...,u_{n-1}), R = diag(r_1,...,r_n)

These are the representations used in DSTEMR. The seminal paper [5] showed bidiagonal factorizations of *definite* matrices always give RRRs. Evidence from [7] suggests LDL* usually gives at least a partial RRR for the smallest eigenvalues, which is why they were chosen for the original MR3. However, once definiteness is lost, the robustness guarantee is lost.

### 3.4 Twisted Factorizations

Combine the upper part of an LDL* and the lower part of a URU*:

```
T = N_k G_k N_k*
```

where N_k is a unit bidiagonal-like matrix:

```
N_k = | 1                    |
      | l_1  1               |
      |      ...  ...        |
      |       l_{k-1} 1 u_{k+1} |
      |              ...  ... |
      |                  1 u_n|
      |                     1 |
```

and G_k = diag(d_1,...,d_{k-1}, gamma_k, r_{k+1},...,r_n).

The twist index k is chosen so that |gamma_k| = O(n|lambda_bar - lambda|), exposing the near-singularity of T - lambda_bar. The system N_k G_k N_k* q_bar = gamma_k e_k is solved in O(n) time to get the eigenvector.

Twisted factorizations perform better than plain LDL* because of the extra flexibility in choosing the twist index k. The 2012 paper [39] provides the detailed qd-type transformations for computing these.

### 3.5 Block LDL* Factorizations

```
T = LDL*    where D is block diagonal with 1x1 or 2x2 blocks D_i
             L is block lower bidiagonal with conformant blocks L_i = (1) or (a 1; 0 1)
```

Blocked twisted factorizations: combine upper blocked LDL* with lower blocked URU*, with a single element (1x1) as the twist pivot gamma_k.

**Key advantage:** Element growth can be countered dynamically by choosing when to use 2x2 blocks. The 2012 paper [38,39] showed how to shift these with mixed relative stability.

**For eigenvector computation:** Even when M has blocks, eigenvectors are still computed using plain (unblocked) twisted factorizations -- G_k is diagonal.

---

## 4. Shift Computation Strategies

### 4.1 The Fundamental Tension

Shifting is the most critical and subtle step. Two conflicting goals:
- **Shift as close as possible** to (or even inside) the cluster to maximize relative gap amplification: relgap_{T-tau}(lambda - tau) = relgap_T(lambda) * |lambda|/|lambda - tau| >> relgap_T(lambda)
- **Ensure the shifted representation satisfies all five requirements** (RRR, ELG, RELGAPS, SHIFTREL, GETVEC) -- and this cannot always be verified a priori.

### 4.2 DSTEMR's Strategy

1. Refine eigenvalues at cluster borders to full precision.
2. Try to shift close to one border (left or right).
3. Evaluate candidate by element growth check only (condition 5.3, C_elg = 8).
4. If no candidate passes, take the one with minimal element growth.
5. Back off (shift less aggressively) if first attempts fail.

### 4.3 XMR's Strategy

1. Refine eigenvalue intervals at cluster borders.
2. Also refine interior of cluster to detect substructure (interior gaps > sqrt(eps)).
3. Compute envelopes per subgroup using [35] techniques.
4. For each shift candidate, check **both**:
   - Element growth: ||Gs|| + ||Us|| + ||Vs|| <= C_elg * spdiam[M_0] (condition 5.3)
   - Relative robustness a priori: ||Upsilon N*s|| ||Upsilon^{-1} N^{-1} s|| <= maxrc (condition 5.4)
5. The candidate passing both tests best (or failing least badly) is selected.
6. Immediately compute shifted representation and set child bounds with inflation: [lambda_i^+] := ([lambda_i](1 +/- alpha) - tau)(1 +/- alpha), alpha = 100*n*eps.
7. Verify child bounds with two Sturm counts each -- reject candidate immediately if inconsistent. This is the **first sign** that something is wrong.
8. If bounds pass, verify outer relative gaps with two more Sturm counts.

### 4.4 Inside Shifts (Experimental)

When cluster substructure reveals interior gaps, XMR can shift to the interior gap boundary (not the cluster midpoint -- that causes fatal element growth). Place the shift near a boundary, close to being singular. Implemented but deactivated in released code.

---

## 5. Twisted Factorization Details

### 5.1 Computing the Eigenvector

Given the twisted factorization T - lambda_bar = N_k G_k N_k*, the eigenvector is found by solving:

```
N_k G_k N_k* q_bar = gamma_k e_k
```

This factors into two O(n) solves (forward and backward substitution on bidiagonal parts). The key property: the canonical vector e_k as the right-hand side is near the position where the true eigenvector has a large entry [12,29]. One step of inverse iteration with this starting vector delivers a residual norm satisfying (3.5), at most a factor sqrt(n) from optimal [19].

### 5.2 Choosing the Twist Index k

For a given eigenvalue approximation lambda_bar, one can always find k such that:
- |gamma_k| = O(n|lambda_bar - lambda|) -- exposing near-singularity
- gamma_k is the minimum pivot in magnitude

The twist index is computed in O(n) using the qd-type transformations from [9,29,39]. In XMR with block factorizations, the twist index is chosen to minimize both conditions (5.3) and (5.4) simultaneously.

### 5.3 Residual and Convergence

After computing the eigenvector, the residual r^leaf = ||(M_tilde - lambda_bar^leaf) q_bar|| / ||q_bar|| must satisfy:

```
r^leaf <= R_gv * n * eps * gap_{M_tilde}({i}; lambda_bar^leaf)
```

If not small enough, Rayleigh Quotient Iteration (RQI) is applied. This updates both the eigenvalue approximation and the eigenvector, with cubic convergence. Fallback to bisection if RQI fails.

### 5.4 Element Growth in Twisted Factorizations

For T = NGN* with twist index k, write:

```
(N - I)G(N - I)* = U + V
```

where U and V are diagonal matrices with nonzero entries:
- U(i,i) = l_{i-1}^2 * d_{i-1} for 1 < i <= k
- V(i,i) = u_{i+1}^2 * r_{i+1} for k <= i < n

The ELG condition (5.3) becomes:
```
||Gs|| + ||Us|| + ||Vs|| <= C_elg * spdiam[M_0]
```
where s is the eigenvector envelope. This can be evaluated in O(n).

---

## 6. Bugs and Numerical Difficulties

### 6.1 DSTEMR Bugs

1. **Underflow in qd-transformations.** The only safeguard was deflecting small pivots from zero using pivmin threshold. A new failure mode was found: underflow during a *sequence* of multiplications and divisions (not in the pivots themselves). This causes DSTEMR to return a **zero vector without error**, which then leads to a completely wrong twist index and a vector with huge norm.

2. **Overly strict RQI acceptance thresholds.** DSTEMR's thresholds for accepting an eigenvector from RQI are too strict, causing unnecessary continuation of iteration or fallback to bisection.

3. **Element growth as sole quality metric.** For some rare synthetic test cases, two shift candidates had comparable element growth, but only one was relatively robust. Element growth alone is insufficient.

4. **Wrong results (WRONG) without warning.** On the Synth testset, DSTEMR produced eigenpairs with orthogonality levels exceeding 10^9 in 0.636% of cases -- effectively wrong answers with no error flag.

5. **FAIL cases.** DSTEMR reported explicit failures (error codes) in 0.062% of Synth cases.

### 6.2 General Numerical Difficulties

1. **Accuracy vs. robustness trade-off.** MR3 produces somewhat higher average residuals than DSTEDC (divide-and-conquer), but much better orthogonality. This is inherent: MR3 computes eigenvectors independently.

2. **Hidden constants in O-terms.** The theoretical bounds ||Tq - q*lambda||/||T||n*eps are guaranteed to be O(1), but for some test cases the hidden constant can be in the hundreds [observation (c) from [6]].

3. **Glued matrices.** Matrices "glued" together (coupling two tridiagonal matrices with a small off-diagonal) create extremely tight clusters that stress-test MR3. These are included in the Synth testset.

4. **The gaptol parameter.** With gaptol = 0.001, orthogonality levels of O(1000*n*eps) are possible even when all five requirements are perfectly satisfied. This is inherent to the worst-case bounds.

---

## 7. Child Representation Computation

### 7.1 The Process (Line 12 of Algorithm MR3)

Given a parent node (M, I, tau_bar) with a cluster I_r = {c,...,d}:

1. **Select a shift tau** near the cluster (see Section 4 above).
2. **Compute M^+ = M - tau** as a new representation using qd-type transformations. This must satisfy SHIFTREL: there exist perturbations M_tilde = erp(M, alpha_down) and H_hat = erp(H, alpha_up) with the exact shift relation M_tilde - tau = H_hat.
3. **Set the child's accumulated shift:** tau_bar_c = tau_bar + tau.
4. **Initialize child eigenvalue bounds:** [lambda_i^+] := [lambda_i] - tau (DSTEMR) or with inflation (XMR).
5. **Add new node (M^+, I_r, tau_bar_c) to set S.**

### 7.2 The Shift Relation (SHIFTREL)

The critical requirement: the shifting must be done in a componentwise mixed relatively stable way. For any node with matrix H computed as child of M using shift tau:

```
M_tilde = erp(M, alpha_down)    and    H_hat = erp(H, alpha_up)
```

such that the exact shift relation M_tilde - tau = H_hat is attained. The constants alpha_down, alpha_up are O(eps) for both twisted factorizations [9] and block factorizations [38,39].

### 7.3 Representation at Child vs Parent

The child representation may be of a **different type** than the parent. For example, the parent might use a block LDL* with 2x2 blocks, while the child uses different block sizes depending on where element growth arises. The framework allows this flexibility because the requirements are stated abstractly.

### 7.4 The Cascade of Perturbations (Figure 4.1)

From root M_0 to the leaf where an eigenpair is computed, the path passes through:

```
M_0 -> M_1 -> M_2 -> ... -> M_d  (exact matrices)
```

via shifts tau_0, tau_1, ..., tau_{d-1}. But in floating point, what actually happens is:

```
M_0 --erp(alpha_down)--> M_0_tilde --(-tau_0)--> M_1_hat --erp(alpha_up)--> M_1
M_1 --erp(alpha_down)--> M_1_tilde --(-tau_1)--> M_2_hat --erp(alpha_up)--> M_2
...
```

The total eigenvalue approximation is lambda_bar_i = tau_0 + tau_1 + ... + tau_{d-1} + lambda_bar^leaf. The telescoping sum in the residual proof (Theorem 4.1) links these perturbations across all levels.

---

## 8. Test Matrices and Results

### 8.1 Test Sets

**Pract (56 matrices):** Real-world matrices from various applications, dimensions up to 6245. Obtained from Osni Marques and Christof Voemel. Same set used to evaluate LAPACK eigensolvers in [6].

**Synth (116,874 matrices):** Synthetic stress-test with very tight clusters. Generated as:
1. Base collection of synthetic types (Wilkinson, Clement, matrices from given eigenvalue distributions, all artificial types from LAPACK's DLATMS test generator [23]), dimensions 2,...,200.
2. Random types: 4 versions per type with different seeds.
3. Each matrix in 3 variants: as-is, glued with small glue (~n*eps*||T||), glued twice with medium glue (~n*sqrt(eps)*||T||).

**Large (72 matrices):** For timing only. Two groups: (a) for each Pract matrix with n > 1000, a DLATMS-generated matrix of same size and eigenvalue distribution; (b) for each synthetic type, a matrix of dimension ~5000.

### 8.2 Results Summary (Table 6.1)

#### Orthogonality |Q*Q - I|/(n*eps)

| Metric | DSTEDC | DSTEMR | XMR-nb | XMR |
|--------|--------|--------|--------|-----|
| **Pract (56 cases)** | | | | |
| AVG | 0.09 | 13.74 | -- | 3.48 |
| MED | 0.04 | 3.84 | -- | 1.71 |
| MAX | 0.50 | 313.36 | -- | 33.95 |
| **Synth (116874 cases)** | | | | |
| AVG | 1617 | 49030 | 234.19 | 1.44 |
| MED | 0.06 | 0.80 | 0.49 | 0.48 |
| MAX | 16.7e6 | 99.5e6 | 9.68e6 | 9473 |
| WRONG (>=10^9) | -- | 0.636% | -- | -- |
| FAIL | -- | 0.012% | 0.062% | 0.062% |

#### Residual ||Tq - q*lambda||/(||T||*n*eps)

| Metric | DSTEDC | DSTEMR | XMR-nb | XMR |
|--------|--------|--------|--------|-----|
| **Pract** | | | | |
| AVG | 0.07 | 0.44 | -- | 0.41 |
| MED | 0.01 | 0.03 | -- | 0.07 |
| MAX | 0.60 | 3.73 | -- | 4.67 |
| **Synth** | | | | |
| AVG | 0.07 | 1.71 | 0.63 | 0.63 |
| MED | 0.05 | 0.07 | 0.13 | 0.13 |
| MAX | 1.34 | 82526 | 13.48 | 13.48 |
| WRONG (>=10^5) | -- | 0.006% | -- | -- |

### 8.3 Key Findings

1. **XMR never broke on Synth.** Zero WRONG results (DSTEMR: 0.636% WRONG). Zero FAIL for eigenpair quality (0.062% FAIL are cases where the routine returned an error code, same rate as XMR-nb).

2. **XMR worst-case orthogonality: 9473** vs DSTEMR's 99.5 million. A reduction of 4 orders of magnitude in worst case.

3. **XMR average orthogonality on Synth: 1.44** vs DSTEMR's 49030. Block factorizations provide significant improvement over XMR-nb (234.19).

4. **Residuals slightly higher** for XMR than DSTEMR on average (Pract: 0.41 vs 0.44), but worst case dramatically better (Synth MAX: 13.48 vs 82526).

5. **DSTEDC is most accurate** for residuals (as expected -- it computes with higher intermediate precision via matrix-matrix multiplications). But DSTEDC has worst-case orthogonality of 16.7 million on Synth.

### 8.4 Timing Results (Large testset)

- XMR average: 3.99 seconds
- DSTEMR average: 4.09 seconds
- DSTEDC average: 8.07 seconds
- XMR was faster in 32 out of 72 cases
- MR3-based methods trend faster for larger matrices
- For n=27069 (SCF application): DSTEDC 1561s, DSTEMR 231s, XMR 245s (sequential). With PARALLEL=12: DSTEDC 220s, XMR parallelizable into 12 subsets with max 28s each.

---

## 9. Key Theorems About Element Growth and Block Factorizations

### 9.1 Theorem 3.1 (Gap Theorem)

For any symmetric A in R^{n x n}, unit vector x, scalar mu, index set I with gap_A(I; mu) != 0:

```
sin angle(x, Q_I[A]) <= ||Ax - x*mu|| / gap_A(I; mu)
```

This is the fundamental result connecting residual norms to eigenvector accuracy. It shows that if the residual ||Ax - x*mu|| is small relative to the gap, the computed vector is close to the true invariant subspace.

### 9.2 Theorem 4.1 (Residual Norms for MR3)

If the representation tree satisfies ELG, SHIFTREL, and GETVEC, then for any eigenpair computed at depth d:

```
||(M_0 - lambda*) q_bar|| <= (||r^leaf|| + gamma * spdiam[M_0]) * (1 + beta_dagger)/(1 - beta_dagger)
```

where gamma = C_elg * n * (d(alpha_down + alpha_up) + alpha_dagger) + 2(d+1)*beta_dagger. This shows the residual grows linearly with depth d, not exponentially.

### 9.3 Theorem 4.5 (Orthogonality of MR3)

If all five requirements hold, with R as defined in (4.2) and d_max the maximum depth:

```
(1/2)|q_bar_i* q_bar_j| <= R*n*eps + C_vecs * d_max * (alpha_down + alpha_up) * n / gaptol
```

This is the core result: orthogonality scales with n*eps times a factor that depends on depth and gaptol. With typical d_max ~ log(n) or O(1), this gives near-machine-precision orthogonality.

### 9.4 Element Growth (ELG) for Twisted Factorizations

For T = NGN* with envelope s:

```
||Gs|| + ||Us|| + ||Vs|| <= C_elg * spdiam[M_0]
```

where U(i,i) = l_{i-1}^2 * d_{i-1} and V(i,i) = u_{i+1}^2 * r_{i+1}. This is evaluable in O(n) given the envelope.

### 9.5 Relative Condition Number and RRR Testing

The relative condition number of eigenvalue lambda_j with respect to twisted factorization T = NGN* is:

```
kappa_rel(lambda_j) := (q_j* N|G|N* q_j) / |q_j* NGN* q_j|
```

Since q_j* NGN* q_j = lambda_j, this simplifies to kappa_rel = |q_j* N Omega N^{-1} q_j| where Omega = diag(+/-1) is the sign matrix of G. This leads to the a priori checkable bound (5.4):

```
||Upsilon N* s|| * ||Upsilon^{-1} N^{-1} s|| <= maxrc
```

for any nonsingular diagonal Upsilon, where s is the envelope. This can be evaluated in O(n).

### 9.6 Block Factorizations: Why They Help

Block factorizations generalize twisted factorizations. The key advantage: **element growth can be countered dynamically**. When a 1x1 block would cause large element growth (large l_i or u_i), using a 2x2 block absorbs the growth. The mixed relatively stable shifting algorithm [38,39] makes this practical.

Empirical evidence (Table 6.1): XMR with blocks achieves average orthogonality 1.44 vs XMR-nb's 234.19 on Synth -- block factorizations reduce worst-case orthogonality by 2 orders of magnitude.

---

## 10. How the Representation Tree Works

### 10.1 Structure

The computation is organized as a **tree traversal**:

- **Root node:** (M_0, I_0, tau_bar = 0) where M_0 is a suitable (preferably definite) representation of T, and I_0 is the index set of desired eigenpairs.
- **Each node** has three features:
  1. A representation of a tridiagonal matrix M
  2. An index set I of "local" eigenpairs to compute
  3. The accumulated shift tau_bar from the root matrix

### 10.2 Processing a Node (Algorithm MR3, Figure 3.1)

```
Input:  T in R^{n x n}, index set I_0 subset {1,...,n}
Output: Eigenpairs (lambda_bar_i, q_bar_i), i in I_0
Params: gaptol (gap tolerance)

1. Find suitable representation M_0 for T (preferably definite, possibly by shifting T)
2. S := {(M_0, I_0, tau_bar = 0)}
3. while S != empty do
4.   Remove one node (M, I, tau_bar) from S
5.   Approximate eigenvalues [lambda_i^loc], classify into singletons and clusters:
     I = I_1 union ... union I_m
6.   for r = 1 to m do
7.     if I_r = {i} then          // singleton
8.       Refine eigenvalue, compute eigenvector q_bar_i via twisted factorization.
         If necessary, iterate (RQI) until residual is small enough.
9.       lambda_bar_i := lambda_i^loc + tau_bar
10.    else                        // cluster
11.      Refine eigenvalue approx at borders (and/or inside) the cluster
12.      Choose shift tau near cluster, compute child M^+ = M - tau
13.      Add (M^+, I_r, tau_bar + tau) to S
14.    endif
15.  endfor
16. endwhile
```

### 10.3 Classification (Line 5)

Eigenvalues are classified using the **gap tolerance** gaptol (e.g., 0.001 for double precision):
- **Singleton/isolated:** relgap_M(lambda_i) > gaptol -- eigenvector can be computed directly with sufficient accuracy
- **Cluster:** eigenvalues with relative gaps below gaptol are grouped together

Care is needed: relgap_T(lambda_i) > gaptol does NOT necessarily imply relgap_T(lambda_{i+1}) > gaptol, so clusters must be formed carefully.

### 10.4 Tree Depth and Complexity

- Each non-root node covers at least 2 eigenvalues (clusters have |I_r| >= 2).
- **Leaves** are singleton eigenpairs -- conceptually nodes but no representation is stored.
- The depth d(j) of eigenpair j equals the number of shifts applied.
- In practice, d_max is typically small (1-3 for most problems).
- **O(kn) total complexity** for computing k eigenpairs, since each node does O(n) work (eigenvalue refinement, shift computation, eigenvector computation).

### 10.5 Tree Traversal Strategies

**Breadth-first** (original MR3 design [11]):
- Saves workspace: representation data (<=2n-1 values) stored in eigenvector matrix, since each node covers >=2 eigenvalues.
- No backtracking possible.

**Depth-first** (proposed in [17]):
- Enables backtracking if a subtree encounters insurmountable problems (e.g., no suitable shift).
- Avoids copying representation data to/from eigenvector matrix.
- Requires capping max depth (e.g., 10) so all path representations fit in reserved workspace.

### 10.6 Example (Figure 3.2)

For 8 eigenvalues with clusters {3:6} and {7:8}:

```
(M_0, [1:8], tau=0)
 |-- singleton: compute q_1, q_2 directly
 |-- cluster {3:6}: shift by tau_{3:6} -> (M_1, [3:6], tau_{3:6})
 |    |-- singleton: compute q_3
 |    |-- cluster {4:6}: shift by tau_{4:6} -> (M_3, [4:6], tau_{3:6} + tau_{4:6})
 |         |-- singletons: compute q_4, q_5, q_6
 |-- cluster {7:8}: shift by tau_{7:8} -> (M_2, [7:8], tau_{7:8})
      |-- singletons: compute q_7, q_8
```

---

## Appendix: The Five Requirements (Summary)

| # | Name | Statement | Purpose |
|---|------|-----------|---------|
| 1 | **RRR** | sin angle(Q_J[M], Q_J[M_tilde]) <= C_vecs * n * alpha / relgap_M(J) for all child index sets J | Small relative perturbations in data cause small changes in invariant subspaces |
| 2 | **ELG** | ||M_tilde - M|| <= spdiam[M_0] AND ||(M_tilde - M) q_bar_i|| <= C_elg * n * alpha * spdiam[M_0] | Element growth from relative perturbations is bounded, especially in directions of local eigenvectors |
| 3 | **RELGAPS** | relgap_M(I_r) >= gaptol for each child index set I_r with |I_r| < n | Each child index set has sufficient relative gap to the rest of the spectrum |
| 4 | **SHIFTREL** | There exist M_tilde = erp(M, alpha_down), H_hat = erp(H, alpha_up) with M_tilde - tau = H_hat | Shifting is done in a componentwise mixed relatively stable way |
| 5 | **GETVEC** | ||r^leaf|| <= R_gv * n * eps * gap_{M_tilde}({i}; lambda_bar^leaf) with bounded perturbations to M and q_bar | Eigenvectors computed at leaves have residuals small relative to the gap |

### Constants (Table 4.1)

| Constant | Source | Meaning | Expected Size |
|----------|--------|---------|---------------|
| C_vecs | RRR | Controls eigenvector sensitivity | ~10 |
| C_elg | ELG | Element growth factor | ~10 |
| alpha_down | SHIFTREL | erp at parent to link with child | O(eps) |
| alpha_up | SHIFTREL | erp at child to link with parent | O(eps) |
| alpha_dagger | GETVEC | erp for vector computations | O(eps) |
| beta_dagger | GETVEC | relative change to vector entries | O(n*eps) |
| R_gv | GETVEC | residual bound factor | up to O(1/gaptol) |

---

## Relevance to Bidiagonal SVD via MR3

This paper addresses the **symmetric tridiagonal eigenproblem** only, not the bidiagonal SVD directly. However, the framework is directly relevant because:

1. The TGK (tridiagonal Golub-Kahan) matrix converts bidiagonal SVD to a symmetric tridiagonal eigenproblem of dimension 2n.
2. The B^T B and BB^T approaches convert to symmetric tridiagonal eigenproblems of dimension n.
3. The five requirements and the representation tree structure apply to any of these formulations.
4. The block factorization techniques could be especially valuable for the TGK matrix, where the zero diagonal creates definiteness issues that plague standard LDL* representations.
5. Willems' thesis [37] explicitly covers "MR3-type Algorithms for the Bidiagonal SVD" but the NCD-aware bidiagonal SVD was never fully implemented (as noted in the project memory).
