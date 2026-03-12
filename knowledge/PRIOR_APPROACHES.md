# Prior Approaches and What We Learned

## Strategies Tried

| Strategy | Accuracy | Scaling | Fatal Flaw |
|----------|----------|---------|------------|
| TGK + DSTEMR (extract U/V from eigvecs) | ~24% pass | O(n²) | Extraction destroys orthogonality |
| TGK + DSTEMR + post-processing | ~24% pass | O(n²) | Post-processing can't fix fundamental extraction failures |
| HGBSVD (Großer-Lang coupling via dbdsgr_) | ~50% pass | O(n²) | Returns INFO!=0 on many matrices |
| TGK + DSTEXR (XMR eigensolver) | ~16% pass | O(n²) | GK-form deactivated, same extraction problem |
| TGK + full MGS reorthogonalization | ~100% pass | **O(n³)** | Global MGS is O(n³) — HARD GATE kills score |
| B^TB eigsolve + one-sided U=BV/σ | ~50% pass | O(n²) | Fails when B^TB is ill-conditioned |
| Hybrid HGBSVD/TGK dispatch | ~60% pass | O(n²) | Current bidiag_svd.h approach |

## Critical Lessons

### 1. TGK eigenvector extraction is THE core problem
Extracting U from odd rows and V from even rows of TGK eigenvectors produces vectors that are individually unit-norm but NOT mutually orthogonal. This is because:
- MR³ (DSTEMR) guarantees eigenvector orthogonality for the FULL 2n-vector
- The perfect shuffle (split into even/odd) does NOT preserve orthogonality
- Only if the eigenvectors have **GK structure** does the split preserve orthogonality
- GK structure requires **NCD-aware shifts** inside the eigensolver (Willems-Lang Theorem 4.5)
- No existing LAPACK routine (DSTEMR, DSTEXR) implements NCD-aware shifts for TGK

### 2. Post-processing helps but isn't sufficient
Adding normalize + sign fix + one-sided recovery + chunked MGS to TGK+DSTEMR improves pass rate but can't fix catastrophic orthogonality failures where eigenvectors are fundamentally wrong. The orthogonality error from extraction can be 10^6 × n·eps.

### 3. HGBSVD (Großer-Lang coupling) works well but has INFO!=0 failures
- `dbdsgr_` returns INFO!=0 on ~50% of test matrices (convergence failure at deep recursion)
- When it works, accuracy is excellent (residual < 1.0, orthoU/V < 1.0)
- The coupling approach itself is mathematically sound — the issue is implementation robustness
- Requires `ftnlen` args: pass `1, 1` at end of each `dbdsgr_` call

### 4. One-sided recovery (U = BV/σ) is powerful but limited
- For well-conditioned singular values: U_i = B·V_i / σ_i gives excellent residual
- Breaks down for tiny σ_i (amplifies errors in V_i by 1/σ_i)
- Threshold: σ_i < sqrt(n·eps)·σ_max → unreliable, better to zero out
- Does NOT fix orthogonality — if V columns aren't orthogonal, U columns won't be either

### 5. Global MGS is O(n³) and forbidden
Full Gram-Schmidt on n vectors of length n is O(n³). Even "selective" MGS on clusters of size k is O(k²·n) — if k=O(n), it's O(n³). Chunked MGS with fixed chunk size (e.g., 32) ensures O(n²) but may miss cross-chunk orthogonality for clusters larger than the chunk.

### 6. Hybrid dispatch is the current practical solution
The current `bidiag_svd.h` tries HGBSVD first, falls back to TGK+DSTEMR:
- HGBSVD handles ~50% of matrices with high accuracy
- TGK+DSTEMR with post-processing catches some of the rest
- Neither path alone achieves >60% pass rate

### 7. What remains unsolved (the open problem)
- No O(n²) algorithm passes all 270 tests
- The root cause: DSTEMR doesn't know about GK structure, so extraction always risks orthogonality loss
- **Willems-Lang Algorithm 4.1** (NCD-aware MR³ on TGK) is the theoretical solution but was never fully implemented — not in LAPACK, not in Willems' XMR code, not anywhere
- Implementing Algorithm 4.1 from scratch requires: NCD detection, NCD-aware shift computation, modified representation tree, GK-structure-preserving eigenvector computation
- Block factorizations for element growth control (Willems-Lang 2012 §5) could help HGBSVD convergence but aren't integrated

## What NOT to try (dead ends)
- **DBDSQR as fallback**: O(n³) → hard gate caps score at 5, even if all tests pass
- **Full MGS reortho**: O(n³) → same hard gate problem
- **Increasing chunk size in chunked MGS**: larger chunks → closer to O(n³), doesn't fix root cause
- **BB^T eigsolve**: same problems as B^TB, no additional benefit
- **Pure inverse iteration on TGK**: O(n²) per vector only for well-separated eigenvalues; for clusters it's O(n³) total
