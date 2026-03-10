# Prior Approaches and What We Learned

12 Python approaches (A–L) were tested in `../bidiag-algo/`. Key lessons:

## Approach Summary

| Approach | Strategy | STCollection Pass | Scaling | Fatal Flaw |
|----------|----------|-------------------|---------|------------|
| A | TGK + eigh_tridiagonal(stemr) | 15/19 | O(n²) | Orthogonality fails on graded/clustered |
| B | B^TB eigsolve only | 12/19 | O(n²) | No coupling → residuals fail |
| C | BB^T eigsolve only | 12/19 | O(n²) | Same as B, transposed |
| D | B^TB + BB^T separately | 17/19 | O(n²) | Coupling still broken |
| E | TGK + full MGS reortho | 19/19 | **O(n³)** | Global MGS is O(n³) |
| F | TGK + selective reortho | 18/19 | O(n²) | Missed one cluster |
| G | B^TB + one-sided recovery | 18/19 | O(n²) | Ill-conditioned B^TB fails |
| H | Custom MR³ solver | 14/19 | O(n²) | Bisection-only, no relative accuracy |
| I | Hybrid B^TB/TGK dispatch | 18/19 | O(n²) | Dispatch heuristic wrong for some matrices |
| J | Approach I + better dispatch | 19/19 | ~O(n²) | Some adversarial patterns > 5x |
| **K** | **Hybrid B^TB/TGK + fixes** | **19/19** | **O(n²)** | **Best result (all 22 adversarial ≤ 5.5x)** |
| L | B^TB only with aggressive reortho | 17/19 | O(n²) | Ill-conditioned B^TB fails |

## Critical Lessons

### 1. TGK eigenvector extraction is the core problem
Extracting U from odd rows and V from even rows of TGK eigenvectors produces vectors that are individually unit-norm but NOT mutually orthogonal. This is because:
- MR³ guarantees eigenvector orthogonality for the FULL 2n-vector
- The perfect shuffle (split into even/odd) does NOT preserve orthogonality
- Only if the eigenvectors have **GK structure** does the split preserve orthogonality
- GK structure requires **NCD-aware shifts** inside the eigensolver (Willems-Lang Theorem 4.5)

### 2. Post-processing helps but isn't sufficient
The TGK+STEMR approach with post-processing (normalize, sign fix, one-sided recovery, chunked MGS) improves from ~20/270 to 65/270 pass. But it can't fix catastrophic orthogonality failures where eigenvectors are fundamentally wrong.

### 3. B^TB path works for well-conditioned matrices
When B^TB's condition number < 1/eps (~4.5×10^15), solving B^TB for eigenvalues → V, then one-sided recovery U = BV/σ, works well. Fails when B^TB condition exceeds this.

### 4. Hybrid dispatch is the practical solution
Approach K uses:
- B^TB path when B^TB diagonal condition < 1/eps
- TGK path as fallback for ill-conditioned matrices
- One-sided recovery (U=BV/σ) for all non-null singular values
- Chunked MGS (MAX_CHUNK=32) for O(n²) reorthogonalization

### 5. Global MGS is O(n³) and forbidden
Full Gram-Schmidt on n vectors of length n is O(n³). Even "selective" MGS on clusters of size k is O(k²·n) — if k=O(n), it's O(n³). Chunked MGS with fixed chunk size (32) ensures O(n²) but may miss cross-chunk orthogonality.

### 6. What remains unsolved
- TGK+STEMR fails on ~76% of adversarial matrices (orthogonality)
- No O(n²) algorithm passes all 270 tests
- The Willems-Lang Algorithm 4.1 (NCD-aware MR³ on TGK) was never fully implemented
- Block factorizations for element growth control exist in theory but not integrated with TGK

## Approach K Architecture (current best)

```
bidiag_svd(n, d, e):
  1. Split at near-zero e[i] into blocks
  2. For each block:
     a. Estimate B^TB condition: btb_max/btb_min
     b. If condition < 1/eps:
        - Solve B^TB eigenvalues → V via eigh_tridiagonal(stemr)
        - One-sided recovery: U = BV/σ
        - Normalize, sign fix
     c. Else (ill-conditioned):
        - Build TGK, solve with eigh_tridiagonal(stemr)
        - Extract U/V from even/odd rows
        - One-sided recovery for small σ
  3. Chunked MGS reortho (MAX_CHUNK=32) within clusters
  4. Fill null columns for zero singular values
```

Source: `../bidiag-algo/approach_k_ncd.py` (Python, 36 KB, 1163 lines)
