# Paper-by-Paper Progression: MR3 for Bidiagonal SVD

## 1990 -- Demmel & Kahan: Accurate Singular Values of Bidiagonal Matrices
- **Contributes:** Zero-shift QR algorithm for bidiagonal SVD with guaranteed relative accuracy. Foundation for dqds. Splitting criteria for negligible off-diagonals.
- **Bugs found in previous work:** N/A (foundational paper).
- **Bugs/limitations in this paper:** O(n^2) per sweep, O(n) sweeps typical but O(n^2) worst case = O(n^3) total. Full SVD only, no subset computation.
- **Fixed by:** dqds refinements (Fernando-Parlett 1994), MR3 approach replaces iterative sweeps with direct eigenvector computation.

## 1997 -- Dhillon: A New O(n^2) Algorithm for the Symmetric Tridiagonal Eigenvalue/Eigenvector Problem (PhD thesis)
- **Contributes:** The MR3 algorithm: representation tree with shifted LDL^T factorizations, twisted factorizations for O(n) eigenvectors, relative gap-based classification. First O(n^2) symmetric tridiagonal eigensolver.
- **Bugs found in previous work:** Bisection+inverse iteration is O(n^3) worst case due to reorthogonalization.
- **Bugs/limitations in this paper:** Three assumptions (4.1: singletons always achievable, 4.2: RRR always exists, 4.4: FP vector always accurate) fail for Wilkinson and glued matrices. Proof tightly coupled to LDL^T representation.
- **Fixed by:** Parlett-Dhillon 2000 (RRR theory), Dhillon-Parlett 2004 (orthogonality proof), Dhillon-Parlett-Vomel 2005 (failure analysis + fixes), Willems-Lang 2013 (framework + block factorizations).

## 2000 -- Parlett & Dhillon: Relatively Robust Representations of Symmetric Tridiagonals
- **Contributes:** Formal definition of RRR. relcond theory: relcond is invariant across all twisted factorizations (Theorem 1). Conditions for when LDL^T is/isn't an RRR (Theorem 3). Differential qd transforms (dstqds, dqds, dtwqds) with mixed relative stability.
- **Bugs found in previous work:** Clarified that LDL^T of indefinite matrices is generally NOT an RRR (element growth destroys robustness).
- **Bugs/limitations in this paper:** No analysis of what happens when RRR doesn't exist. No bidiagonal SVD treatment.
- **Fixed by:** Dhillon-Parlett-Vomel 2005 (random perturbation fix), Grosser-Lang 2001/2005 (SVD coupling), Willems-Lang 2013 (block factorizations).

## 2001/2003 -- Grosser & Lang: An O(n^2) Algorithm for the Bidiagonal SVD
- **Contributes:** First O(n^2) bidiagonal SVD via MR3. Coupling transformations: run RRR on B^TB, derive BB^T factorizations implicitly via differential qd auxiliary variables (s_i, p_i). Backward-stable coupling at rec=1 (Corollary 2.4, no additions/subtractions). X_k diagonal scaling for componentwise-accurate u from v. O(n^2) confirmed with ~4x doubling ratios.
- **Bugs found in previous work:** Showed that separate B^TB/BB^T eigensolves fail residual coupling. Showed naive MR3 on T_GK has alternating-pivot instability.
- **Bugs/limitations in this paper:** Only positive-definite couplings implemented (INFO=3/4/5 for deep recursion). Splittings not supported. Coupling quality degrades for very tight clusters. Code incomplete (RANGE='A' only).
- **Fixed by:** Willems-Lang 2012 (NCD on T_GK as alternative to coupling for tight clusters), but complete integration never achieved.

## 2004 -- Dhillon & Parlett: Orthogonal Eigenvectors and Relative Gaps
- **Contributes:** Complete proof that MR3 produces orthogonal eigenvectors when five requirements hold. Fernando-Parlett (FP) vector: solve N_k G_k N_k^T z = gamma_k e_k in O(n). Angle bound: sin(angle) <= 5n*eps + eigenvalue_error/(|v(r)|*gap) + eps_rel*relcond.
- **Bugs found in previous work:** Dhillon 1997 proof was incomplete for clusters requiring multiple recursion levels.
- **Bugs/limitations in this paper:** Three assumptions from 1997 still not proven to hold universally. Bound has hidden O(1/gaptol) constant.
- **Fixed by:** Dhillon-Parlett-Vomel 2005 (explicit failure analysis), Willems-Lang 2013 (modular five-requirement framework).

## 2005 -- Dhillon, Parlett & Vomel: Glued Matrices and the MRRR Algorithm
- **Contributes:** First systematic analysis of MR3 failure modes. Three explicit assumption failures: (4.1) glued matrices prevent singleton creation, (4.2) Ritz value coincidence prevents RRR existence, (4.4) FP vector underflow for Wilkinson m>=50. Four fixes: random perturbation of root rep (O(n), breaks algebraic structure), selective inverse iteration + MGS for FP failures, submatrix approach for isolated failures, Rayleigh-Ritz projection for multiple FP failures.
- **Bugs found in previous work:** LAPACK 3.0 DSTEMR: 22 complete failures + errors up to 10^14*n*eps. Bug: splitting criterion wrong when no off-diagonal satisfies it but eigenvalues still numerically equal down to underflow. Internal gaptol = 1/n was too loose.
- **Bugs/limitations in this paper:** MGS fallback is O(n*k^2), potentially O(n^3) for large clusters. Random perturbation doesn't fix FP underflow. Glued Wilkinson with p copies still requires O(p) tree levels.
- **Fixed by:** LAPACK 3.1 (splitting fix, gaptol tightened to 10^{-3}). Willems-Lang 2013 (block factorizations reduce element growth, envelope localization helps FP vectors).

## 2005 -- Grosser & Lang: On Symmetric Eigenproblems Induced by the Bidiagonal SVD (LAPACK WN 166)
- **Contributes:** Formal coupling theory: T_GK pivots to B^TB pivots (d_hat_i = -d_bar_{2i-1}*d_bar_{2i}). B^TB to BB^T coupling (Lemma 5.1). Proved coupled factorizations preserve relative accuracy: |lambda_hat - lambda_check| = O(eps*lambda) vs O(mu^2*eps) for separate. Geometric "stretched ellipsoid" interpretation. Proposed xBDSCR for LAPACK.
- **Bugs found in previous work:** Separate B^TB/BB^T eigensolves have absolute (not relative) eigenvalue agreement, amplified by mu^2 = (sigma_max/sigma_j)^2. T_GK LDL^T has strongly alternating pivots for ill-conditioned B.
- **Bugs/limitations in this paper:** xBDSCR never successfully integrated into LAPACK. Coupling formulas complex for indefinite case.
- **Fixed by:** Willems 2010 (NCD approach avoids coupling entirely for some cases), Marques 2020 (BDSVDX uses STEVX instead of MRRR).

## 2008 -- Demmel, Marques, Parlett, Vomel: Performance and Accuracy of LAPACK's Symmetric Tridiagonal Eigensolvers
- **Contributes:** Comprehensive benchmark of QR, BI, DC, MR on 8 architectures. MR scaling: O(n^2.2) on Opteron. QR/DC accuracy O(sqrt(n)*eps), MR accuracy O(n*eps). Test matrix taxonomy: S1/S2/W1-W3/G1-G2/U1-U2/Wi/GW with exact specifications. DC deflation fractions by matrix class. MR flop rate lowest (0.5-0.8 GFlop/s) due to scalar divides.
- **Bugs found in previous work:** LAPACK 3.0 MR3: 22 failures, errors 10^14. LAPACK 3.1: fixed but MR worst-case orthogonality 190*n*eps vs QR's 0.46*n*eps. BI worst-case 440*n*eps.
- **Bugs/limitations in this paper:** Only evaluates tridiagonal eigensolvers, not bidiagonal SVD. MR positive slope in error vs n means accuracy approaches n*eps bound as n grows.
- **Fixed by:** Willems-Lang 2013 (XMR reduces worst-case from 99.5M to 9473). BDSVDX (Marques 2020) addresses the SVD problem.

## 2010 -- Willems: PhD thesis (Bergische Universitaet Wuppertal)
- **Contributes:** Comprehensive MR3 framework. Block factorizations (2x2 blocks in LDL^T). NCD concept for bidiagonal SVD via T_GK. Z-representations for tighter error bounds. XMR implementation. Theory for XMR-TGK and XMR-CPL approaches.
- **Bugs found in previous work:** DSTEMR element-growth-only shift evaluation is insufficient. Coupling approach (BDSCR) has twisted factorization mismatch (p.128, 149).
- **Bugs/limitations in this paper:** DBDSGR implementation incomplete (DLARRI missing). GK-form support in DLAXRE deactivated (caused problems). Block factorization proof for coupling contains subtle error (Theorem 5.2).
- **Fixed by:** Willems-Lang 2012 (published NCD theory), Willems-Lang 2013 (published XMR framework).

## 2012 -- Willems & Lang: The MR3-GK Algorithm for the Bidiagonal SVD (ETNA 39)
- **Contributes:** Algorithm 4.1: MR3 on T_GK with NCD shifts. Theorem 4.5: correctness proof with GK structure. Definition 4.3 (GK structure), Definition 4.6 (NCD), Lemma 4.4 (equal-norm). XMR-TGK vs XMR-CPL comparison. XMR-TGK: zero failures, orth avg 5.35 (Pract), resid avg 0.35. Refuted the "MR3 on T_GK is fundamentally flawed" opinion.
- **Bugs found in previous work:** Black-box MR3 on T_GK fails because initial outside shift clusters +/-sigma pairs (Example 4.1). Coupling approach proof error in Theorem 5.2 of thesis.
- **Bugs/limitations in this paper:** XMR-TGK worst-case orth 3095 (Synth), resid 118. XMR-CPL: 24/19240 failures. Element growth after one T_GK shift can be 1/eps. Block factorizations mentioned but not integrated. Prototype only, not production.
- **Fixed by:** Partially by Willems-Lang 2013 (improved XMR framework). Complete solution remains open.

## 2013 -- Willems & Lang: A Framework for the MR3 Algorithm (SIAM J. Sci. Comput. 35)
- **Contributes:** Modular, representation-agnostic proof of MR3 correctness via five abstract requirements. Block LDL^T factorizations with dynamic 2x2 blocks. XMR implementation: orth avg 1.44 (Synth), worst 9473 (vs DSTEMR's 99.5M). A priori relative robustness checking (condition 5.4) beyond just element growth. Structured bisection with shift selection. Zero WRONG results.
- **Bugs found in previous work:** DSTEMR underflow in qd-transformations causes zero vector return. DSTEMR RQI thresholds too strict. Element growth alone insufficient for shift quality (rare Synth failures).
- **Bugs/limitations in this paper:** AVGAPFAC bug (0.1 instead of 0.3, missing DEPTH scaling) discovered later by Marques. STEXR fails on graded-spectrum matrix. Framework is for tridiagonal eigenproblem only; SVD via T_GK requires additional NCD machinery from 2012 paper. XMR not integrated into LAPACK.
- **Fixed by:** Willems provided AVGAPFAC fix (email 2014, 4.5 years after discovery). STEXR abandoned; no maintained successor.

## 2020 -- Marques, Demmel, Vasconcelos: Bidiagonal SVD via Associated Tridiagonal Eigenproblem (ACM TOMS 46)
- **Contributes:** BDSVDX: production bidiagonal SVD using T_GK + STEVX (bisection+inverse iteration). 250-matrix test suite with conditions 1 to 10^200. CHKBD matrix analysis exposing reorthogonalization insufficiency. Documented STEXR failure on n=10 graded matrix. Plans to replace STEVX with MRRR "when available."
- **Bugs found in previous work:** STEXR orthogonality 3.95e4 on 10x10 graded matrix. BDSCR abandoned. STEMR NaN safeguard issues. DBDSVDX's own post-extraction reorthogonalization trigger (`ORTOL=SQRT(ULP); ABS(NRMU-ORTOL)*SQRT2.GT.ONE` in dbdsvdx.f lines 344, 645-646) fails on CHKBD — norm test passes but extracted u/v vectors not orthogonal. Note: this is NOT a DSTEIN bug (DSTEIN already uses eigenvalue-separation-based reortho); it is DBDSVDX's norm-based trigger that is insufficient.
- **Bugs/limitations in this paper:** BDSVDX outliers: CHKBD matrices have orthU/orthV ~ 10^13. Matrix 36 (dim 1260, Lanczos) and matrix 225 (dim 396) show results >10. STEVX is O(n^2) per eigenvector worst case, not O(n) like MRRR. No O(n^2) guarantee for full SVD. Paper proposes trigger should use eigenvalue separation and magnitude, not just norms, but does not implement the fix. DBDSVDX source states plan to replace STEVX with MRRR "when available" (citing Willems-Lang 2013).
- **Fixed by:** Open problem. No subsequent publication has achieved O(n^2) bidiagonal SVD with full accuracy guarantees. DBDSVDX reortho trigger unchanged in LAPACK 3.12.1 (Aug 2025). MRRR replacement never implemented.
