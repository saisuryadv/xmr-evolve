# Self-Contained Fortran-Only MR3-GK Bidiagonal SVD

O(n²) bidiagonal SVD using MR³ with Golub-Kahan structure preservation
(Willems & Lang, ETNA 2012). **Pure Fortran** — no Python needed to
compute an SVD. Python is included only for running the test suite.

**Result: 379 / 379 tests passing.**

## What's here

The SVD is computed entirely in Fortran:

- **`mr3gk_fortran/*.f90`** (~1,400 lines) — Fortran orchestration:
  T_GK construction, two-phase splitting, sign matrices, zero-shift QR
  deflation, Bv recovery, sign fix, Gram-Schmidt completion, σ rescaling.

- **`xmr_src/`** (29 files) — Paul Willems' XMR eigenvector kernel
  (critical-path subset only).

- **`dlaxre_gk.f`**, **`dlaxrb_clssfy_fix.f`** — our two patches to the
  XMR kernel (Bug #1: activate GK branch; Bug #2: AVGTHRESH=0 guard).

- **`xmr_wrapper.c`** — thin C glue exposing `dlaxre_` + `dlaxrv_` for
  the Fortran orchestration to link against.

- **`mr3_gk.py`** — a minimal stub (50 lines) that calls the Fortran
  executable via subprocess. Exists only so `evaluate.py` can import
  `bidiag_svd`. **Not used to compute anything.**

## Prerequisites

```bash
sudo apt install gfortran gcc liblapack-dev libblas-dev python3-numpy python3-scipy
```

## Build

```bash
bash build.sh                              # → libxmr.so (94 KB)
cd mr3gk_fortran && bash build.sh && cd .. # → mr3gk_run
```

## Run

```bash
./mr3gk_fortran/mr3gk_run input.bin output.bin
```

**Input:** `int32 n` + `n × float64` diagonal + `(n−1) × float64` super-diagonal.
**Output:** `int32 info` + `n × float64` sigma + `n×n × float64` U (col-major) + `n×n × float64` V (col-major).

To call from Fortran code directly, link against `mr3gk_fortran/*.o` +
`libxmr.so` and call `dmr3gk_svd(n, d, e, sigma, U, V, ldu, info)`.

## Evaluate

```bash
python3 evaluate.py
# expect: TOTAL: 379/379 passed
```

## References

- Willems & Lang, ETNA 39 (2012) — Algorithm 4.1
- Willems & Lang, ETNA 38 (2011) — block factorizations
- Demmel & Kahan, SIAM J. Sci. Stat. Comput. (1990) — zero-shift QR
