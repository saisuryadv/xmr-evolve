# Self-Contained MR3-GK Bidiagonal SVD

O(n²) bidiagonal SVD using MR³ with Golub-Kahan structure preservation
(Willems & Lang, ETNA 2012). This folder is a **self-contained copy** of
the complete implementation — it builds and runs independently of the rest
of the repository.

**Result: 379 / 379 tests passing** on the standard evaluation suite.

## What's here

This implementation has **two equivalent paths** to the same SVD:

### Python + Fortran hybrid
Python orchestration (`mr3_gk.py`) calls the Fortran XMR kernel
(`libxmr.so`) via ctypes for the eigenvector computation.

### Pure-Fortran driver
Fortran orchestration (`mr3gk_fortran/*.f90`) calls the same Fortran
XMR kernel directly. Produces **bit-identical** singular values and
U/V agreeing to ≤ 1 ULP vs the Python hybrid.

Both paths use the same `libxmr.so` (the XMR kernel + our 2 patches).

## File inventory (68 files, 676 KB)

```
self-contained-bidiagsvd/
├── build.sh                  # builds libxmr.so from xmr_src/ + patches
├── xmr_src/                  # 29 critical-path Willems XMR Fortran files
├── dlaxre_gk.f               # patch: activate GK branch (Bug #1 fix)
├── dlaxrb_clssfy_fix.f       # patch: AVGTHRESH=0 guard (Bug #2 fix)
├── xmr_wrapper.c             # C glue (dlaxre_ + dlaxrv_ interface)
├── xmr_ctypes.py             # Python ctypes binding to libxmr.so
├── mr3_gk.py                 # Python SVD orchestration (638 lines)
├── mr3gk_fortran/             # Pure-Fortran SVD orchestration
│   ├── build.sh              #   builds mr3gk_run
│   ├── mr3gk.f90             #   public entry: dmr3gk_svd
│   ├── mr3gk_split.f90       #   two-phase bidiagonal splitting
│   ├── mr3gk_qrsweep.f90     #   zero-shift QR deflation
│   ├── mr3gk_tgk.f90         #   T_GK construction + per-block solve
│   ├── mr3gk_postproc.f90    #   normalize, Bv recovery, sign fix, GS
│   ├── mr3gk_utils.f90       #   DSTEBZ wrapper
│   ├── mr3gk_consts.f90      #   EPS, GAPTOL, etc.
│   └── mr3gk_run.f90         #   CLI driver (binary I/O)
├── evaluate.py               # 379-test scoring harness
├── full_eval.py              # 90 adversarial pattern generators
├── run_evaluate_fortran.py   # routes evaluate.py through Fortran driver
├── test_fortran_match.py     # Py↔Fortran bit-identity test
└── stcollection/             # 19 STCollection real-world bidiagonals
```

## Prerequisites

```bash
sudo apt install gfortran gcc liblapack-dev libblas-dev python3-numpy python3-scipy
```

## Build

```bash
# Step 1: build the XMR kernel (libxmr.so)
bash build.sh

# Step 2 (optional): build the pure-Fortran driver
cd mr3gk_fortran && bash build.sh && cd ..
```

## Run — Python + Fortran

```python
from mr3_gk import bidiag_svd
import numpy as np

d = np.array([4.0, 3.0, 2.0, 1.0])   # bidiagonal diagonal
e = np.array([0.5, 0.5, 0.5])         # bidiagonal super-diagonal
sigma, U, V, info = bidiag_svd(d, e)
print("singular values:", sigma)

# verify: B ≈ U · diag(sigma) · V^T
B = np.diag(d) + np.diag(e, 1)
print("residual:", np.max(np.abs(B - U @ np.diag(sigma) @ V.T)))
```

## Run — pure Fortran

```bash
./mr3gk_fortran/mr3gk_run input.bin output.bin
```

**Input format:** `int32 n`, then `n × float64` diagonal, then `(n−1) × float64` super-diagonal.
**Output format:** `int32 info`, then `n × float64` sigma, then `n×n × float64` U (column-major), then `n×n × float64` V (column-major).

## Run — Fortran via evaluate.py harness

```bash
python3 run_evaluate_fortran.py --suite evaluate
# expect: TOTAL: 379/379 passed
```

## Evaluate

```bash
# Python+Fortran on all 379 tests:
python3 evaluate.py
# expect: TOTAL: 379/379

# Fortran-only on all 379 tests:
python3 run_evaluate_fortran.py --suite evaluate
# expect: TOTAL: 379/379
```

## Test thresholds

A test passes iff:
- residual `‖B − UΣVᵀ‖ / (‖B‖·n·ε) ≤ 7.0`
- orthogonality `‖UᵀU − I‖ / (n·ε) ≤ 5.0`
- orthogonality `‖VᵀV − I‖ / (n·ε) ≤ 5.0`

## References

- Willems & Lang, "A Framework for the MR³ Algorithm," ETNA 39, 2012
- Willems & Lang, "Block factorizations and qd-type transformations," ETNA 38, 2011
- Demmel & Kahan, "Accurate Singular Values of Bidiagonal Matrices," SIAM J. Sci. Stat. Comput. 1990
