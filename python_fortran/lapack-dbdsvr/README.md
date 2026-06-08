# DBDSVR — bidiagonal SVD driver for LAPACK

`DBDSVR` is a drop-in `DBDSVDX`-API bidiagonal SVD driver:

- **Singular values** via `DLASQ1` (dqds): high relative accuracy on every
  singular value, regardless of cluster structure.
- **Singular vectors** via Willems' MR³-GK algorithm: a multiple
  relative-robust-representation eigensolver specialised to the
  Golub–Kahan tridiagonal of the bidiagonal matrix.

The motivation: `DBDSVDX` tridiagonalises B^T·B implicitly via
Golub–Kahan and then runs MRRR on the resulting tridiagonal. Sigma
accuracy there is limited by tridiagonal-eigenvalue accuracy. `DBDSVR`
decouples the two — DLASQ1 gets the sigmas to ulp accuracy directly,
and MR³-GK is used only for the vector pairs. The two outputs are
reconciled by **positional pairing**: sort MR³-GK sigmas descending,
pair the k-th MR³ column with `S(k)`. This keeps clustered triples
intact where greedy nearest-match would mis-pair.

## API

```fortran
SUBROUTINE DBDSVR( UPLO, JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
                   NS, S, Z, LDZ, WORK, LWORK, IWORK, INFO )
```

Identical to `DBDSVDX` apart from the `LWORK` argument added to make
the workspace query well-defined. Only `JOBZ='V'`, `RANGE='A'` are
supported in this initial version. `LWORK = -1` returns required
workspace (`2*N*N + 7*N`) in `WORK(1)` and `N` in `IWORK(1)`.

`INFO = N+1` signals a soft warning: a sigma reported by DLASQ1 and
the corresponding MR³-GK sigma disagree by more than
`16*N*EPS*||B||`. Output is still produced.

## Accuracy

Measured on the project's adversarial test suite (270 cases across
`n ∈ {10, 100, 200, 400}` × pathological matrix families: tight
clusters, glued Wilkinson, huge condition number, Demmel S1ps, …) plus
109 cases from the `stcollection` reference matrices:

```
TOTAL: 270/270 passed (100.0%)
TOTAL: 379/379 passed (100.0%)
```

A passing case has residual `||B − U·diag(S)·V^T||/||B|| ≤ N·EPS`,
column-wise orthogonality `max|U^T·U − I|, max|V^T·V − I| ≤ 5.0` on
the suite's normalised metric, and `INFO == 0`.

## Thread safety

Run `./dbdsvr_thread_test` (built by the `build.sh` chain together with
`libxmr.so` and `mr3gk_run`):

```
OMP_NUM_THREADS=4 LD_LIBRARY_PATH=. ./dbdsvr_thread_test
  case 1 (n=32): OK
  case 2 (n=48): OK
  case 3 (n=64): OK
  case 4 (n=80): OK
PASS: 4 parallel cases bit-identical to serial
```

This is bit-identity, not float-tolerance. To get there, the legacy
`COMMON /XMRSTATS/` block — 24 profiling counters with `SAVE` shared
across all calls — was removed from every `xmr_src/*.f` file plus
`dlaxre_gk.f` and `dlaxrb_clssfy_fix.f`. The counters were write-only
and never read by the SVD path. No other `COMMON` or `SAVE` state
remains anywhere in the call graph.

## Building

```
bash build.sh                          # libxmr.so + parent patches
cd mr3gk_fortran && bash build.sh      # mr3gk modules, dbdsvr.f90,
                                       # mr3gk_run CLI
python3 evaluate.py                    # 270 adversarial + 109 stcoll
```

`mr3gk_run` is a binary-protocol CLI driver around `DBDSVR` used by
`evaluate.py` (Python harness builds B, calls the CLI, scores residual
/ orthogonality).

## Known follow-ups for the LAPACK PR

1. **F90 → F77.** Reference-LAPACK's `SRC/` is still primarily F77
   fixed-form. `dbdsvr.f90` uses F90 features (`use … only`,
   allocatable locals). A straight port replacing allocatables with
   slices of `WORK(*)` is needed before PR submission.
2. **`JOBZ='N'`, `RANGE='V'`, `RANGE='I'`** are not yet implemented.
   Trapped via `XERBLA` for now.
3. **Auxiliary-routine doc headers.** `dmr3gk_svd`,
   `mr3gk_tgk_multiblock`, and the `xmr_src/dlax*` files still use the
   original Willems comment style. LAPACK style (`*> \brief`,
   `*> \par Purpose`, `\param`) on each is a deferred polish step.

## Layout

```
dbdsvr.f90                  # the new LAPACK-style driver
dbdsvr_thread_test.f90      # OpenMP bit-identity smoke test
dlaxre_gk.f                 # GK-aware root representation (Willems patch)
dlaxrb_clssfy_fix.f         # depth-0 classification bug fix
xmr_src/                    # Willems XMR Fortran source (vendored)
mr3gk_fortran/              # F90 MR³-GK pipeline modules
   mr3gk*.f90 + mr3gk_run.f90
build.sh                    # libxmr.so build
mr3gk_fortran/build.sh      # mr3gk + dbdsvr + mr3gk_run build
evaluate.py / full_eval.py  # Python test harness
stcollection/               # reference test matrices
```
