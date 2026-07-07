# Offline PR: DBDSVR alongside DBDSVDX

This directory is a self-contained bundle intended for offline review
before opening a Reference-LAPACK pull request.  It layouts new files
under the same tree as `Reference-LAPACK`:

```
SRC/dbdsvr.f                          <- new wrapper (LAPACK style)
SRC/{dbdsvdmr3,dbdtgk,dlar1v_tgk,dlarrf_tgk,dlarrv_tgk}.f
SRC/stegr_ID/dlarrb.f                 <- advisor's TGK-MR^3 engine (unchanged)
TESTING/EIG/dchkbsvr.f                <- new residual-sweep driver
TESTING/EIG/dbsvrmg.f                 <- new matrix generator
TESTING/EIG/derrbsvr.f                <- new argument error checker
TESTING/dbsvr_test.f                  <- new main program
TESTING/svrtest.in                    <- input file (default: N=2,10,30,60,100)
Makefile / build.sh                   <- standalone build (not for merge)
```

`SRC/dbdsvr.f` mirrors `DBDSVDX`'s signature; under the hood it stages
private copies of `D`, `E`, invokes the advisor's `DBDSVDMR3`, and
packs the requested slice into `S`, `Z` following DBDSVDX's `Z = [U; V]`
convention.

## Test coverage

`DCHKBSVR` runs the **DBDSVDX residual battery from `dchkbd.f`
tests 25–34** verbatim (residual formulas, sortedness check, and
σ-comparison denominator all lifted from
`/tmp/lapack-ref/TESTING/EIG/dchkbd.f` lines 1330–1490), and extends
it with an analogous RANGE='A' battery numbered 20–24 for a full-
spectrum tie-breaker.  Per matrix, fifteen residuals are computed:

Tests 20–24 (`RANGE='A'`):
  * (20) `‖SA - Uᵀ B Vᵀᵀ‖ / (n·‖B‖·ulp)` — `DBDT04`
  * (21) `‖I - Uᵀ U‖ / (n·ulp)` — `DORT01('Columns')`
  * (22) `‖I - Vᵀ Vᵀᵀ‖ / (n·ulp)` — `DORT01('Rows')`
  * (23) `1/ulp` if `SA` not non-increasing / non-negative
  * (24) `maxⱼ |SA(j)-S₂(j)|` σ-comparison against a companion `JOBZ='N'` call

Tests 25–29 (`RANGE='I'`) and 30–34 (`RANGE='V'`): same shape.

`VL, VU` for the RANGE='V' block are picked from the RANGE='A'
spectrum `SA` per `dchkbd.f:1374-1388` byte-identical.

Threshold `THRESH = 20.0`, mirroring `dchkbd`'s `bd.in`.

`TSTREF=T` in the input file makes the driver *also* run `DBDSVDX` on
the same generated matrices with the same metrics — an in-place
side-by-side.

Matrix types (see `TESTING/EIG/dbsvrmg.f` for details):

| id | source                       | recipe                                           |
|---:|------------------------------|--------------------------------------------------|
|  1 | dchkbd type 1                | Zero bidiagonal                                  |
|  2 | dchkbd type 2                | Identity                                         |
|  3 | dchkbd type 3  (DLATMS mode 4) | Arithmetic-spectrum bidiag                     |
|  4 | dchkbd type 4  (DLATMS mode 3) | Clustered-spectrum bidiag                      |
|  5 | dchkbd type 16               | Log-distributed bidiag on [ulp, 1/ulp]           |
|  6 | dchkst type 8  (mode 4)      | Symm-tridiag arith → Cholesky → bidiag           |
|  7 | dchkst type 9  (mode 3)      | Symm-tridiag clustered → Cholesky → bidiag       |
|  8 | dchkst type 10 (mode 1)      | Symm-tridiag random-log → Cholesky → bidiag      |
|  9 | Wilkinson W_{2k+1}^+         | Hand-built T → shift → Cholesky → bidiag         |
| 10 | Glued Wilkinson              | 2×W blocks + tiny glue → Cholesky → bidiag       |
| 11 | dchkbd type 5                | Arith-spectrum × √overflow                       |
| 12 | dchkbd type 6                | Arith-spectrum × √underflow                      |

Cholesky recipe for types 6–10 is
`T_shift = T - (Gershgorin_lower - 0.001·max(|gl|,1))·I`,
`T_shift = LDL^T` via `DPTTRF`, then
`B = L·√D_chol` so that `B^T B = T_shift`.  This mirrors
`python_fortran/self-contained-fortran-bidiagsvd/test_dense_to_bidiag.py`.

## Result on the default sweep (N ∈ {2,10,30,60,100})

With `TSTREF=T`, `TSTERR=T`, `THRESH=20.0`:

```
BSR argument error checks passed.
DCHKBSVR summary:  1800 residuals tested,     5 failures at THRESH.
```

Per solver:

| solver   | tested | failed | max residual |
|----------|-------:|-------:|-------------:|
| DBDSVR   |    900 |      2 |         22.0 |
| DBDSVDX  |    900 |      3 |         4.50e13 |

DBDSVR: both failures are at N=60 type 6 (arithmetic-tridiag → Cholesky
→ bidiag), tests 21 and 22 — U/Vᵀ orthogonality drift of 21.8 and 22.0
(just above THRESH=20).  Reconstruction (test 20) and every RANGE='I'
and RANGE='V' residual pass.

DBDSVDX: three failures, all at N=100.  Type 4 (clustered
`DLATMS`) sees `‖I − UᵀU‖ ≈ 4.19e12` in the RANGE='A' block only —
reconstruction and Vᵀ orthogonality are fine (0.027, 0.062).  Type 8
(Cholesky-lifted random-log tridiagonal) shows U/Vᵀ orthogonality
≈4.5e13 in the RANGE='I' block only — the RANGE='A' and RANGE='V'
blocks pass.  These are both interior-eigenvector orthogonality issues
inside `DSTEVX` on ill-conditioned inputs, and give an independent
motivation for the MR³-based `DBDSVR`.

## Build and run

```sh
cd python_fortran/lapack_pr_dbdsvr
bash build.sh
```

`build.sh` builds `/tmp/lapack-ref` if absent (needs `git` and `cmake`),
compiles this bundle, then runs
`./build/dbsvr_test < TESTING/svrtest.in`.

`gfortran >= 8` is required (`ieee_arithmetic` intrinsic module is used
inside the advisor's engine).  On RHEL/CentOS 7, `source
/opt/rh/devtoolset-11/enable` first.

## What is NOT in this bundle

- CMake integration inside the LAPACK tree (`SRC/CMakeLists.txt`,
  `TESTING/EIG/CMakeLists.txt` updates).  Should be added when opening
  the real PR.
- Doc updates in `DOCS/`.
- Complex counterpart `ZBDSVR`.

## Attribution

Everything under `SRC/{dbdsvdmr3,dbdtgk,dlar1v_tgk,dlarrf_tgk,dlarrv_tgk}.f`
and `SRC/stegr_ID/` is copied byte-identical from the advisor's
`BidiagonalSVD_TGK.zip` (Willems–Lang TGK-rooted MR^3 driver).  The
wrapper `SRC/dbdsvr.f` and every file under `TESTING/` are new for this
PR.
