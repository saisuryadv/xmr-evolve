# Offline PR: DBDSVR alongside DBDSVDX

This directory is a self-contained bundle intended for offline review
before opening a Reference-LAPACK pull request.  It layouts new files
under the same tree as `Reference-LAPACK`:

```
SRC/dbdsvr.f                          <- new wrapper (LAPACK style)
SRC/{dbdsvdmr3,dbdtgk,dlar1v_tgk,dlarrf_tgk,dlarrv_tgk}.f
SRC/stegr_ID/dlarrb.f                 <- advisor bisection + progress guard
TESTING/EIG/dchkbsvr.f                <- new residual-sweep driver
TESTING/EIG/dbsvrmg.f                 <- new matrix generator
TESTING/EIG/derrbsvr.f                <- new argument error checker
TESTING/dbsvr_test.f                  <- new main program
TESTING/dbdsvr_fallback_test.f90      <- INFO/nonfinite fallback regression
TESTING/svrtest-dbdsvr-only.in        <- canonical DBDSVR-only input
TESTING/svrtest.in                    <- DBDSVR/DBDSVDX comparison input
Makefile / build.sh                   <- standalone build (not for merge)
```

`SRC/dbdsvr.f` mirrors `DBDSVDX`'s signature; under the hood it stages
private copies of `D`, `E`, invokes the advisor's `DBDSVDMR3`, and
packs the requested slice into `S`, `Z` following DBDSVDX's `Z = [U; V]`
convention.

DBDSVR uses `DBDSDC` as a result-based safety path. Every order tries
DBDSVDMR3 first. If MR3 returns a positive `INFO`, an incomplete spectrum,
or a NaN/infinity in the requested S/U/V output, DBDSVR restores the inputs
and retries with DBDSDC. A solver-only probe showed that all 11 large cases
historically labelled `TIMEOUT` return from MR3; their old 600-second limit
covered the entire executable, including cubic-cost validation. There is
therefore no size-only backend cutoff. On successful nontrivial calls,
`IWORK(2)` records the accepted path as documented in `SRC/dbdsvr.f`.
The bundled `DLARRB` additionally rejects a nonfinite, wrong-direction, or
non-progressing bracket expansion and caps both bracket growth and refinement.
This converts the observed negative-`WERR` infinite loop into positive `INFO`,
which reaches the same DBDSDC fallback.

## Test coverage

`DCHKBSVR` runs the **DBDSVDX residual battery from `dchkbd.f`
tests 25–34** verbatim (residual formulas, sortedness check, and
σ-comparison denominator all lifted from
Reference LAPACK's `TESTING/EIG/dchkbd.f` lines 1330–1490), and extends
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

| id | source                         | recipe                                         |
|---:|--------------------------------|------------------------------------------------|
|  1 | dchkbd type 1                  | Zero bidiagonal                                |
|  2 | dchkbd type 2                  | Identity bidiagonal                            |
|  3 | dchkbd type 3 (mode 4)         | Arithmetic-spectrum bidiag                     |
|  4 | dchkbd type 4 (mode 3)         | Geometric-spectrum bidiag                      |
|  5 | dchkbd type 16                 | Log-distributed bidiag on [ulp, 1/ulp]         |
|  6 | dchkst type 8 (mode 4)         | Symm-tridiag arithmetic → Cholesky → bidiag    |
|  7 | dchkst type 9 (mode 3)         | Symm-tridiag geometric → Cholesky → bidiag     |
|  8 | dchkst type 10 (mode 1)        | Symm-tridiag clustered → Cholesky → bidiag     |
|  9 | dchkbd type 5                  | Arithmetic bidiag × √overflow                  |
| 10 | dchkbd type 6                  | Arithmetic bidiag × √underflow                 |
| 11 | dchkst type 3 (mode 4)         | Diagonal arithmetic                            |
| 12 | dchkst type 4 (mode 3)         | Diagonal geometric                             |
| 13 | dchkst type 5 (mode 1)         | Diagonal clustered                             |
| 14 | dchkst type 6                  | Diagonal arithmetic × √overflow                |
| 15 | dchkst type 7                  | Diagonal arithmetic × √underflow               |
| 16 | dchkst type 21 (mode 3)        | SPD tridiag geometric → Cholesky → bidiag      |

Cholesky recipe for types 6–10 is
`T_shift = T - (Gershgorin_lower - 0.001·max(|gl|,1))·I`,
`T_shift = LDL^T` via `DPTTRF`, then
`B = L·√D_chol` so that `B^T B = T_shift`.  This mirrors
`python_fortran/self-contained-fortran-bidiagsvd/test_dense_to_bidiag.py`.

## Result on the canonical sweep (N ∈ {2,10,30,60,100})

With the DBDSVR-only input, GNU Fortran 12.2.0, and Reference LAPACK
3.11.0 on the c2 Linux reproducer:

```text
DCHKBSVR summary:  1200 residuals tested,    10 failures at THRESH.
```

All ten failures occur at N=30: eight metrics for type 7 (maximum 25.8)
and two metrics for type 16 (maximum 30.2). Every solver call returns
normally. See `../docs/dbdsvr_investigation_2026-08-04.md` for the exact
output, environment, four-solver benchmark, and cross-suite analysis.

## Build and run

```sh
cd python_fortran/lapack_pr_dbdsvr
bash build.sh
```

`build.sh` checks out Reference LAPACK v3.11.0 under
`/tmp/lapack-ref-3.11.0` if absent (needs `git` and `cmake`), builds the
required LAPACK, BLAS, and TMG libraries, compiles this bundle from clean
objects, verifies the positive-`INFO` and nonfinite-output fallback paths,
then runs
`./build/dbsvr_test < TESTING/svrtest-dbdsvr-only.in`. The numerical failure
count is compiler/platform-sensitive; the exact emailed count uses the c2
environment documented in the investigation report.

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

The TGK-rooted MR^3 engine originated in the advisor's
`BidiagonalSVD_TGK.zip` (Willems–Lang). This bundle adds DBDSVR result-based
fallback, propagates `DLARRB` failures through `DLARRV_TGK`, and adds bounded
progress checks to the bundled `DLARRB`. The wrapper `SRC/dbdsvr.f` and every
file under `TESTING/` are new for this PR.
