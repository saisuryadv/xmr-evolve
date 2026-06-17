# BidiagonalSVD_TGK (DBDSVDMR3) — cross-suite evaluation (2026-06-12)

Driver: `DBDSVDMR3` from advisor's `BidiagonalSVD_TGK/` (origin/main 15a946d).
Binary: `python_fortran/BidiagonalSVD_TGK/test_stcoll_alloc` (allocatable
STCollection driver; only addition to advisor tree).
Thresholds (same as the rest of the repo): res ≤ 7 n·ε, ortU/ortV ≤ 5 n·ε.

## PASS table

| Suite | Generator | PASS / Total | Pass rate | Wall |
|---|---|---|---|---|
| pract | `test_synth_pract.generate_pract` (STCollection bidiag) | 73 / 127 | 57.5% | — |
| synth (mode=paper) | `test_synth_pract.generate_synth` (Willems-Lang 2012) | 13960 / 18438 | 75.7% | 788.9 s |
| 379 | `lapack-dbdsvr/full_eval` 90 adv × 4 sizes + 19 STColl | 335 / 375 | 89.3% | 56.1 s |
| dense_to_bidiag | `test_dense_to_bidiag` DENSE + PAPER patterns × 4 sizes | 192 / 224 | 85.7% | 67.6 s |

## Worst-case observations

- **pract**: catastrophic on `T_W21_g_1e+12` (res ≈ 1.98e6 n·ε, ortU ≈ 8.3e6 n·ε)
  — strong-glue Wilkinson, the documented accuracy floor in `BENCHMARKS.txt`.
- **synth**: 7 entries overflow even the absolute-line formatter (parse-fail =
  catastrophic). The worst finite cases are the `tri*_gm` glued-tridiagonal
  family (e.g. `tri4_n3_gm` n=9 → ortU ≈ 2.6e6 n·ε).
- **379**: worst orth on `adv:gl_wilkm@100` ≈ 1.46e13 n·ε (only visible after
  fallback-parsing the absolute-value line; the n·ε column overflowed Fortran's
  F11.1 format).
- **dense_to_bidiag**: worst orth on `glued_wilk_5x21_sqrteps` ≈ 240 n·ε; all
  other failures are mild relative to the synth/pract glued cases.

## Log files

- `docs/pract_dbdsvdmr3_2026-06-12.log`
- `docs/synth_dbdsvdmr3_2026-06-12.log`
- `docs/lapack_dbdsvr_379_dbdsvdmr3_2026-06-12.log`
- `docs/dense_to_bidiag_dbdsvdmr3_2026-06-12.log`

## Orchestrator

`python_fortran/eval_dbdsvdmr3.py --suite {pract|synth|379|dense_to_bidiag|all}`

Subprocess pattern (write `.dat` → run `test_stcoll_alloc` → parse the
`resid/(n.eps.||B||)= … orthU/(n.eps)= … orthV/(n.eps)= …` line; fall back to
the absolute-values line and divide by n·ε whenever the F11.1 column overflows
to `***`).
