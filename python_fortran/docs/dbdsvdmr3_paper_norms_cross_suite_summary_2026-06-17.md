# DBDSVDMR3 paper-norm sweep — cross-suite, two solvers (2026-06-17)

Same 4 suites as the 2026-06-12 max-norm run, but now using the
**Willems-Lang 2012 Table 5.1** metrics:
- orthogonality = `max(||UᵀU − I||₂, ||VᵀV − I||₂) / (n·ε)` (spectral norm via DSYEV / numpy)
- residual = `max_i max(||B·vᵢ − uᵢ·σᵢ||₂, ||Bᵀ·uᵢ − vᵢ·σᵢ||₂) / (||B||·n·ε)`
- thresholds unchanged: `res ≤ 7 n·ε`, `ortU/ortV ≤ 5 n·ε`

Two solvers, identical inputs:
- **advisor** — `DBDSVDMR3` from `BidiagonalSVD_TGK/` (origin/main 15a946d). Metrics added inside `dev/test_stcoll_alloc.f` (DSYEV + per-triplet DNRM2).
- **selfcontained** — `mr3_gk.bidiag_svd` from `self-contained-fortran-bidiagsvd/` (Fortran via subprocess + binary I/O). Metrics computed in Python.

## PASS table (paper norms)

| Suite | advisor | selfcontained |
|---|---|---|
| pract           | 64 / 127  (50.4%) | 115 / 127  (90.6%) |
| synth (paper)   | 13400 / 18438 (72.7%) | 17644 / 18438 (95.7%) |
| 379             | 317 / 375 (84.5%) | 368 / 375 (98.1%) |
| dense_to_bidiag | 178 / 224 (79.5%) | 221 / 224 (98.7%) |
| **all**         | **13959 / 19164 (72.8%)** | **18348 / 19164 (95.7%)** |

## Wall-clock

| Suite | advisor | selfcontained |
|---|---|---|
| pract           | 15492.1 s (4.3 h, many TIMEOUTs) | 753.9 s |
| synth (paper)   | 1210.7 s | 600.4 s |
| 379             | 1275.1 s | 33.5 s |
| dense_to_bidiag | 84.4 s | 49.6 s |

The advisor pract run is dominated by `T_*` matrices with `n ∈ {3258..9941}`
that exceed the 600 s per-matrix timeout. Self-contained handles the same
matrices in seconds because it uses the in-process MR³-GK binary that scales
better and shares a process for the metric step.

## Worst cases (top 1 per suite × solver, paper-norm units of n·ε)

| Suite | solver | worst matrix | ortU | ortV | res |
|---|---|---|---|---|---|
| pract           | advisor       | T_Alemdar_1 (n=6245)            | TIMEOUT | TIMEOUT | TIMEOUT |
| pract           | selfcontained | T_TSC_OPF_300 (n=9774)          | TIMEOUT | TIMEOUT | TIMEOUT |
| synth           | advisor       | ev2_ec1_n35_s0_gm (n=105)       | 4.29e13 | 4.29e13 | 0.048   |
| synth           | selfcontained | ev8_ec1_n72_s3_gm (n=216)       | inf (ValueError) | inf | inf |
| 379             | advisor       | adv:gl_wilkp@400                | TIMEOUT | TIMEOUT | TIMEOUT |
| 379             | selfcontained | adv:gl_wilkw@200 (n=200)        | 2.25e13 | 1.690 | 0.073 |
| dense_to_bidiag | advisor       | dense:dense_three_clusters@10   | 0.716 | 4.50e14 | 4.50e14 |
| dense_to_bidiag | selfcontained | dense:dense_wilkinson_sv@100    | 31.16 | 31.15 | 0.054 |

Self-contained failures are concentrated on **Wilkinson-glue** patterns
(`gl_wilkw`, `gl_wilk2w`, `dense_wilkinson_sv`), exactly the family the paper
flags as the hardest case. Advisor failures are broader and include
catastrophic blow-ups on benign-looking cases (`dense_three_clusters@10`).

## Log files (all in `python_fortran/docs/`)

### Advisor
- `pract_paper_advisor_2026-06-17.log`
- `synth_paper_advisor_2026-06-17.log`
- `lapack_dbdsvr_379_paper_advisor_2026-06-17.log`
- `dense_to_bidiag_paper_advisor_2026-06-17.log`

### Self-contained
- `pract_paper_selfcontained_2026-06-17.log`
- `synth_paper_selfcontained_2026-06-17.log`
- `lapack_dbdsvr_379_paper_selfcontained_2026-06-17.log`
- `dense_to_bidiag_paper_selfcontained_2026-06-17.log`

## Reproducer

```sh
cd python_fortran
bash BidiagonalSVD_TGK/build.sh        # builds advisor's test_stcoll_alloc
for solver in advisor selfcontained; do
  for suite in 379 dense_to_bidiag pract synth; do
    python3 eval_dbdsvdmr3.py --suite $suite --solver $solver --paper-norms \
        2>&1 | tee docs/${suite}_paper_${solver}_2026-06-17.log
  done
done
```

`--paper-norms` switches both the Fortran driver (advisor) and the Python
adapter (selfcontained) to the Willems-Lang Table 5.1 metric definitions.
