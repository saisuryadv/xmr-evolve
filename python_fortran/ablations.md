# Ablation Studies

## 1. Block Factorization (USEBLOCKS)

### What it is

Block factorizations are a generalization of standard LDL* twisted factorizations
introduced in Willems-Lang 2011 ("Block factorizations and qd-type transformations
for the MR3 algorithm", ETNA 38, pp. 363-400). Instead of D being diagonal (1×1
pivots), D may contain 2×2 blocks. L is block lower bidiagonal with conforming
structure.

During the shift T_GK - τI = L·D·L*, whenever a pivot d_plus(i) would become
dangerously small (signaling element growth), a 2×2 block is created dynamically
to absorb the growth. The criterion is:

```
|d_plus(i) * c_plus| < (1/8) * gnsq(i)
```

where gnsq(i) = e(i)^2 is the squared off-diagonal.

### Why it matters

The 2012 paper (Willems-Lang, "Computing the bidiagonal SVD using MR3",
Example 4.8) shows that shifting a GK matrix by α ≪ 1 produces element growth
~1/α² in the standard LDL*. For B = [[1,1],[0,α]], shifting by -α gives
D(2) = (1-α²)/(α·(2-α²)) ≈ 1/α. This destroys the Relatively Robust
Representation (RRR) property needed for MR3 orthogonality guarantees.

The 2013 framework paper (Table 6.1) compares XMR (with blocks) vs XMR-nb
(without): "the blocked factorizations provide a noticeable further improvement"
in worst-case orthogonality on the Synth testset.

### How to reproduce

The switch is a compile-time parameter in `xmr_src/dlaxrs_stat.f`, line ~68:

```fortran
LOGICAL, PARAMETER :: USEBLOCKS = .TRUE.   ! default: blocks enabled
```

To disable:

```bash
cd python_fortran

# Create modified source
sed 's/USEBLOCKS = .TRUE./USEBLOCKS = .FALSE./' \
    xmr_src/dlaxrs_stat.f > /tmp/dlaxrs_stat_noblk.f

# Compile just the modified file
gfortran -fPIC -O2 -std=legacy -w -fno-second-underscore \
    -c /tmp/dlaxrs_stat_noblk.f -o fortran_objects/dlaxrs_stat.o

# Relink (don't run build.sh — it would recompile from xmr_src/)
gcc -shared -o libxmr.so fortran_objects/*.o \
    -lgfortran -llapack -lblas -lm

# Run evaluation
python3 evaluate.py
```

To restore:

```bash
bash build.sh   # recompiles everything from xmr_src/ with USEBLOCKS=TRUE
```

### Results (full 379-test suite)

| Configuration | Pass rate | Score |
|---|---|---|
| **USEBLOCKS=TRUE** (default) | **379/379 (100%)** | **92.58** |
| USEBLOCKS=FALSE | 373/379 (98.4%) | 91.89 |

### Failing tests without block factorization

| Test | ortU (blocks) | ortU (no blocks) | Degradation |
|---|---|---|---|
| spike@200 | 0.451 | **15.197** | 34× |
| demmel_S1ps@100 | 3.100 | **7.375** | 2.4× |
| demmel_S1ps@400 | 1.282 | **26.815** | 21× |
| demmel_S1pe_k4@400 | 1.476 | **22.367** | 15× |
| demmel_S1pe_k8@100 | 3.944 | **7.777** | 2.0× |
| demmel_S1pe_k8@400 | 1.624 | **29.635** | 18× |

All 6 failures share the same pattern: matrices with **clusters of tiny singular
values** where MR3 must shift T_GK by a tiny τ. Without 2×2 blocks to absorb the
element growth (~1/τ²), the shifted LDL* representation loses its RRR property.

The degradation grows with n: `demmel_S1ps` goes from ortU=7.4 at n=100 to 26.8
at n=400, while staying under 3.1 with blocks at both sizes.

### Interpretation

Block factorizations are essential for robustness on matrices with tiny clustered
singular values. They prevent catastrophic element growth during the shift step of
MR3's representation tree, preserving the relative robustness property that
guarantees orthogonal eigenvectors. The 2×2 blocks act as numerical "shock
absorbers" — they store the product d_plus(i)·d_plus(i+1) as a single 2×2
determinant instead of two separate pivots, avoiding the cancellation that destroys
accuracy.
