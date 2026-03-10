# MR3-Bidiagonal SVD: Condensed Knowledge Index

**See also:**
- [BASELINES.md](BASELINES.md) — Comparison of 4 baseline algorithms (DBDSQR, DBDSVDX, TGK+STEMR, TGK+STEXR) with empirical results
- [EVALUATION.md](EVALUATION.md) — How evaluation works: metrics, scoring, all 90 test patterns, data layout
- [PRIOR_APPROACHES.md](PRIOR_APPROACHES.md) — 12 prior approaches (A-L), key lessons, why each failed
- [PROGRESSION.md](PROGRESSION.md) — Paper-by-paper timeline (1990-2020)

## The Algorithm: Willems-Lang Algorithm 4.1 (MR3 on T_GK)

**Input:** Upper bidiagonal B (n x n), index set I_0. **Output:** Singular triplets (sigma_i, u_i, v_i).

1. Form T_GK(B): 2n x 2n symmetric tridiagonal with zero diagonal and off-diagonals interleaving a_1, b_1, a_2, b_2, ..., a_n. Eigenvalues are +/-sigma_i.
2. Take M_0 := T_GK(B) as root representation directly (NOT via LDL^T factorization -- heterogeneous tree).
3. For first-level shifts M^+ = T_GK(B) - mu, use a special routine exploiting the zero diagonal. Small positive shift near +sigma cluster keeps -sigma far away, avoiding the clustering disaster.
4. Build representation tree (depth-first): classify eigenvalues by relative gap >= gaptol (~0.001). Singletons get eigenvectors via RQI on twisted factorizations. Clusters get shifted child representations.
5. At each shift, check NCD condition: ncd(-mu_bar, 32*n*eps). Children use Z-representations for tighter error bounds.
6. Extract singular vectors: v(j) = sqrt(2) * q(2j-1), u(j) = sqrt(2) * q(2j) from eigenvector q of T_GK.

**Why no initial outside shift:** Standard MR3 shifts outside the spectrum to get positive definiteness. Algorithm 4.1 uses T_GK directly. The first shift handles definiteness, and +sigma/-sigma never get lumped together.

## Theorem 4.5: What It Guarantees

If the representation tree satisfies five requirements (RRR, ELG, RELGAPS, SHIFTREL, GETVEC) AND each node's invariant subspace nearly has GK structure (via NCD), then:
- **Orthogonality:** max{cos(angle(u_i, u_j)), cos(angle(v_i, v_j))} <= 2*sqrt(2)*A, where A = orth_GK + C_vecs * n * xi_GK / gaptol
- **Norm deviation:** max{| ||u_i||-1 |, | ||v_i||-1 |} <= sqrt(2)*A + O(A^2)
- **Residuals:** max{||B*v_i - u_i*sigma_i||} <= sqrt(2) * resid_GK

Practical expectation: orth = O(n*eps), resid = O(||B||*n*eps). The critical extra requirement beyond standard MR3 is that representations can be perturbed to have local invariant subspaces with GK structure, captured by NCD.

## NCD (Nearly Constant Diagonal)

A representation M has NCD if an elementwise relative perturbation M_tilde = erp(M, xi) exists such that all diagonal entries of M_tilde equal a constant c. Even-dimensioned symmetric tridiagonals with constant diagonal are shifted GK matrices T_GK(B) - mu, so NCD implies near-GK structure for invariant subspaces. Practical check for LDL^T: |(LDL^T)(i,i) - const| = O(eps) * max{|d_i|, |l_{i-1}^2 * d_{i-1}|}. All level-one representations are automatically NCD (from SHIFTREL). The danger is at deeper levels: if element growth was large at level 1, level 2 may lose NCD. XMR acceptance criterion: ncd(-mu_bar, 32*n*eps).

## GK Structure

A subspace S of R^{2n} has GK structure if extracting u_i and v_i from its orthonormal basis via the perfect shuffle yields orthonormal sets {u_i} and {v_i} separately. Key properties: (1) shift-invariant -- any invariant subspace of T_GK belonging to one half of the spectrum has GK structure; (2) equal-norm -- any vector in a GK-structured subspace has ||u|| = ||v||; (3) no direct test exists -- NCD is used as a sufficient proxy.

## Element Growth Problem

GK matrices are highly vulnerable to element growth from shifting. Example: B = [1, 1; 0, alpha] with alpha ~ eps. Shifting T_GK(B) by -alpha produces LDL^T with D(2) ~ 1/alpha ~ 1/eps ~ 10^16, while matrix entries are O(1). Despite this, the LDL^T may still be NCD. But a second shift on this LDL^T would almost certainly destroy NCD. Block factorizations (2x2 blocks in D) counter element growth dynamically by absorbing growth into blocks. Willems-Lang 2011 describes the theory; XMR-TGK uses scalar factorizations with NCD check as safeguard. Block factorizations reduce worst-case orthogonality from ~234 to ~1.44 (average, Synth testset).

## ALL Known Bugs

**1. DSTEXR (XMR) AVGAPFAC bug** -- `dlaxrb_clssfy.f`. Average gap tolerance too aggressive (0.1 instead of 0.3) without depth scaling. Causes premature cluster splitting, loss of orthogonality on STCollection matrices. **Fixed**: AVGAPFAC=0.3, AVGTOL scales with DEPTH.

**2. DSTEXR orthogonality failure on graded spectrum** -- Appendix B of Marques 2020. 10x10 tridiagonal with lambda_j = c^((j-1)/9), c=1/sqrt(eps). Orthogonality: 3.95e4 * n*eps. First four eigenvectors linearly dependent. STEVX and STEMR pass the same matrix. Root cause: MRRR implementation failing on very close eigenvalues. STEXR is no longer maintained.

**3. DSTEMR underflow in qd-transformations** -- Willems-Lang 2013. Underflow during multiplication/division sequences (not in pivots) causes DSTEMR to return a zero vector without error flag. Wrong twist index chosen, producing huge-norm vector. Additional: overly strict RQI thresholds cause unnecessary fallback to bisection.

**4. DSTEMR NaN safeguard issues** -- LAPACK GitHub (2018). Safeguards against NaN may cause early termination even for benign matrices.

**5. DSTEMR wrong results without warning** -- 0.636% of Synth cases produced orthogonality > 10^9. LAPACK 3.0 had 22 complete failures; 3.1 fixed splitting criterion and tightened gaptol from 1/n to 10^{-3}.

**6. BDSCR (coupled B^TB/BB^T MRRR) abandoned** -- Willems et al. 2006. "Potential mismatch in twisted factorizations" caused large SVD residuals. Abandoned entirely.

**7. DBDSVDX reorthogonalization insufficient** -- Marques 2020. The bug is in DBDSVDX's own post-extraction reorthogonalization (NOT in DSTEIN — DSTEIN already uses eigenvalue-separation-based reortho). After DSTEIN returns TGK eigenvectors, DBDSVDX extracts u/v parts and checks whether to reorthogonalize using a norm-based trigger: `ORTOL = SQRT(ULP); IF (NRMU.NE.ONE .AND. ABS(NRMU-ORTOL)*SQRT2.GT.ONE)` (dbdsvdx.f lines 344, 645-646). For CHKBD matrix B_8x8 (a_i=10^(-(2i-1)), b_i=10^(-(2i-2))), the extracted u/v norms look fine (~1/sqrt(2)), so the trigger never fires, yet orthU/orthV ~ O(10^5) to O(10^13). Fix proposed: trigger must account for eigenvalue separation and magnitude, not just norms. **Not fixed** as of LAPACK 3.12.1 (Aug 2025) — code unchanged since 2015. DBDSVDX source explicitly states DSTEVX is "to be replaced with a version of the Multiple Relative Robust Representation algorithm" (referencing Willems-Lang 2013), but this has not happened.

**8. DBDSGR (Grosser-Lang) incomplete** -- Only positive definite couplings implemented. INFO=3/4/5 returned for deep recursion, tight clusters, or MGS-requiring clusters. DLARRI missing from XMR distribution. Splittings not supported (INFO=33).

**9. XMR dlaxrt.f method 2** -- Commented out with note "something is wrong with method 2." Fallback (method 1) in use.

**10. XMR dlaxre.f GK-form** -- GK-form support deactivated; caused problems with GRO synthetics.

## ALL Test Matrix Construction Formulas

**DLATMS Modes (eigenvalue distribution):**
MODE 1: D(1)=1, D(2:n)=1/kappa. MODE 2: D(1:n-1)=1, D(n)=1/kappa. MODE 3: D(i)=kappa^(-(i-1)/(n-1)) (geometric). MODE 4: D(i)=1-(i-1)/(n-1)*(1-1/kappa) (arithmetic). MODE 5: random in (1/kappa,1), log-uniform. MODE 6: random from matrix distribution.

**Wilkinson W_{2m+1}^+:** diag = |m-(i-1)|, off-diag = 1. Eigenvalue pairs differ by O(1/m!). FP vectors underflow for m>=50.

**Glued Wilkinson:** p copies of W_{2m+1}^+ connected by off-diagonal gamma (typically sqrt(eps)). Creates O(p) copies of each tight pair.

**CHKBD (LAPACK SVD tester):** a_i = 10^(-(2i-1)), b_i = 10^(-(2i-2)). Condition ~ 10^22 for n=8. Reorthogonalization trigger fails.

**Grosser-Lang test matrices (IDs 110-244):**
110: all-ones spectrum. 111: uniform eps-apart. 112: uniform sqrt(eps)-apart. 113: uniform eps-to-1. 115: geometric eps-to-1. 117: random. 118: clustered at 1. 119: clustered at +/-1. 120: clustered at eps. 121: clustered at +/-eps. 200-203: ABCON (heat equation tridiagonal with alpha=2, beta=1). 210: random bidiagonal. 220/221: graded (GRADP/GRADM). 222-225: Wilkinson variants. 230: Clement. 240-243: GRO (near-constant with outliers). 244: bWILK+.

**Demmel 2008 synthetic types:** S1pe/S1ne/S1ps/S1ns (strongly clustered, n-1 at 1/kappa), S2pe/S2ne/S2ps/S2ns (n-1 at 1), W1-W3 (weakly clustered), G1/G2 (geometric), U1/U2 (uniform), Wi (Wilkinson), GW (glued Wilkinson). Default kappa=1/eps.

**Marques 2020 failure matrix:** lambda_j = (1/sqrt(eps))^(-(j-1)/9), n=10. Orthogonality failure: 3.95e4.

## The DBDSVDX Bug (from slides)

There is no single "DBDSVDX bug." The fundamental dilemma is that **no black-box MR3 approach achieves both orthogonality AND residual accuracy simultaneously**:
- Option 1 (B^TB/BB^T separately): Orthogonality PASS, Residual FAIL (u,v computed independently, not coupled).
- Option 2 (T_GK directly): Residual PASS, Orthogonality FAIL (+/- sigma pairs create representation tree challenges).
BDSVDX uses STEVX (bisection + inverse iteration) on T_GK as a workaround, not MRRR. The CHKBD matrix (a_i=10^(-(2i-1))) exposes DBDSVDX's post-extraction reorthogonalization insufficiency: after DSTEIN computes TGK eigenvectors, DBDSVDX extracts u/v parts and applies a norm-based trigger (`ABS(NRMU-ORTOL)*SQRT2.GT.ONE`) to decide whether to reorthogonalize. For CHKBD, the norms look correct but vectors are not orthogonal — the trigger never fires. Note: DSTEIN itself uses eigenvalue-separation-based reortho and is not the source of this bug. The "Holy Grail" is achieving both properties in O(n^2). DBDSVDX source code states the plan to replace STEVX with MRRR "when available" (citing Willems-Lang 2013), but this remains unimplemented as of LAPACK 3.12.1.

## XMR Code Structure (Condensed Call Graph)

```
DSTEXR (driver)
  DLAXRA    - split into irreducible blocks, sign/scale
  DLAXRI    - map global indices to per-block (uses DLAXRK for bisection)
  DLAXRE    - root rep: LDL^T + DQDS eigenvalues (DLASQ1), perturbation (8 ULPs)
  DLAXRV    - depth-first rep tree traversal (MAXDEPTH=10)
    DLAXRB_CLSSFY  - classify singletons/clusters (BUGFIX: AVGAPFAC=0.3)
    DLAXRB_REFSNG  - refine singletons via vectorized negcount (DLAXRM)
    DLAXRB_REFCLS  - refine cluster bounds
    DLAXRF_ENV     - eigenvector envelope (DLAXRF_GRPENV -> DLAXRT -> DLAXRG)
    DLAXRF         - find shift for child rep
      DLAXRF_SELSHF  - generate shift candidates
      DLAXRS         - blocked shift factorization (DLAXRS_STAT, DLAXRS_PROG)
      DLAXRF_SELTW   - select twist (MAXRELCOND=10, MAXGROWTH=8)
      DLAXRF_COB/IIB - check/init child bounds
    DLAXRX         - RQI + bisection for singletons (DLAXRT -> DLAXRG)
  DLAXRO    - sort eigenpairs
```

**45 files, ~12K lines Fortran + C wrapper. Key innovations vs DSTEMR:** block-aware 2x2 LDL^T, vectorized multi-shift negcount (2/4/8/16/32/64 unrolled), envelope localization, Z-representations for children.

**Representation data:** REPR (4N+3 reals: G, BDET, NGN, GNSQ, PIVBASE), REPI (6+N+N/2 ints: TYPE, TWIST, NB, OMEGA, LBBEGK). Eigenvalue list: EWL_LU (bracketing intervals), EWL_AE (shared-bound indices).

## What Works, What Doesn't, and Why

**Works:**
- XMR-TGK on practical matrices: orthogonality avg 5.35, max 48.40 (in units of n*eps). Residuals: avg 0.35, max 4.19. Zero failures.
- XMR (generic tridiag): orthogonality avg 1.44 on Synth, worst 9473 (vs DSTEMR's 49030 avg, 99.5M worst). Zero WRONG results.
- Block factorizations: reduce worst-case orthogonality by 2 orders of magnitude over non-blocked XMR.
- Coupling approach (Grosser-Lang): O(n^2) with ~4x doubling ratios. Works for well-separated and moderately clustered singular values.

**Doesn't work:**
- STEXR on graded-spectrum matrices (orthogonality 3.95e4). No longer maintained.
- BDSCR coupled MRRR: abandoned due to twisted factorization mismatch.
- XMR-CPL: 24/19240 Synth cases fail (no shift satisfying consistency bounds). Proof contains subtle error.
- Naive MR3 on T_GK: initial shift clusters +/- sigma pairs for small sigma. Eigenvectors lose equal-norm property.
- Normal equations independently: orthogonal u,v separately but ||Bv-sigma*u|| can be O(sigma) for clusters.
- BDSVDX on CHKBD: orthogonality O(10^13) because DBDSVDX's post-extraction reorthogonalization trigger is norm-based only (not DSTEIN's fault — DSTEIN uses eigenvalue separation). Unfixed as of LAPACK 3.12.1.
- Deep recursion coupling (Grosser-Lang): at rec>1, coupling uses mixed explicit/coupled steps that amplify errors.

**Why it remains open:** The fundamental tension is that MR3's representation tree builds shift sequences optimized for one formulation (B^TB, BB^T, or T_GK), but the SVD requires coupling between left and right singular vectors. NCD-aware T_GK solves the coupling problem but introduces element growth and Z-representation overhead. Block factorizations help but are not yet integrated with T_GK for production use. The "Holy Grail" requires all of: NCD shifts inside the eigensolver, block factorizations for element growth, and either coupling or GK-structure preservation -- no single implementation has achieved all three.

## Key Numerical Thresholds and Parameters

| Parameter | Value | Purpose |
|-----------|-------|---------|
| gaptol | 0.001 | Singleton/cluster classification threshold |
| xi_GK (NCD tol) | 32*n*eps | NCD acceptance for shifts in XMR-TGK |
| AVGAPFAC | 0.3 | Average gap factor (bugfix from 0.1) |
| MAXRELCOND | 10 | Max relative condition for twist acceptance |
| MAXGROWTH | 8 | Max element growth for twist acceptance |
| C_elg | ~8-10 | Element growth constant |
| EVALOK | 5.0 | Child rep quality threshold |
| MAXDEPTH | 10 | Max representation tree depth |
| NULPROOTPERT | 8 | Root rep perturbation in ULPs |
| KBLOCK | 1/8 | Block creation threshold in shift factorization |
| sqrt(eps) reorth | ~1.5e-8 | BDSVDX reorthogonalization trigger (insufficient for CHKBD) |
| ENVGAPFAC | 100 | Cluster width < mingap/100 for envelope computation |
| RESDECFAC | 0.7 | Required residual decrease per RQI step |
