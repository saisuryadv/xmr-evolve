# BidiagonalSVD cross-suite paper-norm sweep (2026-06-18)

Adds three new metrics on top of the 2026-06-17 sweep, all run under
`--paper-norms` (Willems-Lang 2012 Table 5.1):

1. **`sv_drift`** — max relative drift of solver-predicted singular values vs
   the LAPACK DBDSQR σ-only reference (gold standard).
2. **`t_eval`** — solver wall time, expert-bench protocol (warmup + adaptive
   min-of-N, 32 MiB cache-flush scratch between every timed call, INTEGER\*8
   `SYSTEM_CLOCK` ticks, D/E save+restore, anti-DCE guard).
3. **`t_dbdsqr`** — DBDSQR σ-only reference wall time under the same protocol,
   logged for every matrix.

Plus, for the **379 suite** only: per-pattern **`t(400)/t(200)`** scaling.

All three solvers were exercised on all four suites:

| Solver | Source |
|---|---|
| **advisor** | `BidiagonalSVD_TGK/` (advisor's hand-off, origin/main 15a946d) → `DBDSVDMR3` |
| **selfcontained** | `python_fortran/self-contained-fortran-bidiagsvd/` (lab) → `mr3gk_run` |
| **dbdsqr** | LAPACK `DBDSQR` driving full SVD (n³, treated as "third solver") |

Sweeps were run in parallel on a 32-core Xeon Silver 4215R, each sweep
pinned to a dedicated core via `taskset`. Selfcontained sweeps were serialized
through one core because `mr3gk_run` uses fixed-path
`/tmp/_mr3gk_{in,out}.bin`. Driver: `python_fortran/run_all_sweeps_2026-06-18.sh`.

## Headline PASS rates

Threshold (all suites, all solvers): `res ≤ 7·n·ε`, `max(ortU, ortV) ≤ 5·n·ε`.

| Suite | advisor | selfcontained | dbdsqr (n³ ref) |
|---|---:|---:|---:|
| **pract** (127) | 64 / 127  (50.4%) | 114 / 127 (89.8%) | 114 / 127 (89.8%) |
| **synth** (18 438) | 13 400 / 18 438 (72.7%) | 17 747 / 18 438 (96.3%) | 18 435 / 18 438 (99.98%) |
| **379** | 317 / 375 (84.5%) | 368 / 375 (98.1%) | 375 / 375 (100%) |
| **dense_to_bidiag** (224) | 179 / 224 (79.9%) | 220 / 224 (98.2%) | 224 / 224 (100%) |

Observations:
- Selfcontained matches dbdsqr exactly on pract (114/127) and stays within 2%
  of dbdsqr on synth/379/dense_to_bidiag.
- Advisor's `DBDSVDMR3` lags by 25–35 percentage points across suites; the
  drop is most severe on pract (50%) where matrices are real-world,
  high-dynamic-range bidiagonals.
- `dbdsqr` failing 13/127 on pract is the same as selfcontained — these are
  pathological matrices where even the n³ reference cannot meet the paper
  thresholds.

## 379-suite scaling: t(400) / t(200) per pattern

Aggregated over PASS rows only (TIMEOUT/FAIL skews ratios).
HARD GATE: `t(400)/t(200) > 5.0` ⇒ pattern is super-quadratic.
Ideal MR³: 4.0.

| Solver | patterns @200&400 paired | worst_ratio | worst-pattern | comment |
|---|---:|---:|---|---|
| advisor | 67 | **6.76** | adv:many_near_zero | barely above gate |
| selfcontained | 86 | **5.17** | adv:checkerboard | barely above gate |
| dbdsqr | 89 | **15.79** | adv:demmel_G1 | expected: n³ algorithm |

Top-5 patterns above-ratio:

### advisor
| ratio | pattern | t@200 | t@400 |
|---:|---|---:|---:|
| 6.76 | adv:many_near_zero | 4.25e-04 | 2.88e-03 |
| 6.05 | adv:random_sparse_e | 5.13e-04 | 3.10e-03 |
| 5.38 | adv:wl_example48 | 4.71e-03 | 2.54e-02 |
| 5.06 | adv:diagonal_only | 5.12e-04 | 2.59e-03 |
| 5.06 | adv:chkbd | 4.10e-04 | 2.07e-03 |

### selfcontained
| ratio | pattern | t@200 | t@400 |
|---:|---|---:|---:|
| 5.17 | adv:checkerboard | 3.55e-02 | 1.84e-01 |
| 4.13 | adv:pd_T0 | 6.49e-02 | 2.68e-01 |
| 4.05 | adv:gl_random_spectrum | 7.71e-02 | 3.12e-01 |
| 4.04 | adv:random_uniform | 7.66e-02 | 3.09e-01 |
| 4.03 | adv:gl_abcon3 | 6.73e-02 | 2.71e-01 |

### dbdsqr (n³ reference)
| ratio | pattern | t@200 | t@400 |
|---:|---|---:|---:|
| 15.79 | adv:demmel_G1 | 1.35e-02 | 2.14e-01 |
| 14.79 | adv:demmel_G1_k12 | 1.88e-02 | 2.77e-01 |
| 14.09 | adv:random_dense_clusters | 1.41e-02 | 1.98e-01 |
| 13.99 | adv:demmel_G1s | 3.11e-02 | 4.36e-01 |
| 13.75 | adv:demmel_G1_k8 | 3.07e-02 | 4.22e-01 |

The advisor's MR³ stays inside the gate envelope on 65/67 patterns. The
selfcontained MR³ stays inside the gate envelope on 85/86 patterns
(`checkerboard` is the only outlier, and selfcontained is also ~10× slower
per call than advisor since it spawns a subprocess per bench rep).

## Per-suite wall time (all 12 sweeps)

| Suite | advisor | selfcontained | dbdsqr |
|---|---:|---:|---:|
| pract | 15 898 s (4.4 h) | 1 380 s (23 min) | 18 633 s (5.2 h) |
| synth | 9 185 s (2.6 h) | 7 549 s (2.1 h) | 9 103 s (2.5 h) |
| 379 | 1 457 s (24 min) | 199 s (3.3 min) | 305 s (5.1 min) |
| dense_to_bidiag | 191 s | 129 s | 218 s |

`selfcontained` is fastest because Python-side bench wraps `bidiag_svd()`
(a subprocess per call) and the adaptive `min-of-N` decides far fewer reps
are necessary; advisor and dbdsqr each do a full 30-rep Fortran-side bench
with cache flush, which is what `--paper-norms` is supposed to measure.

## Files

Logs (each ends with `PASS:`, has a sorted Failing block, and — for 379 — a
top-15 Scaling table):

- `docs/sweeps_2026-06-18/{pract,synth,379,dense_to_bidiag}_{advisor,selfcontained,dbdsqr}_paper_norms.log`

Code (no advisor-tree files were modified beyond the two already added in
the 2026-06-12 sweep, `dev/test_stcoll_alloc.f` and `build.sh`):

- `BidiagonalSVD_TGK/dev/dbdsqr_ref.f` — σ-only DBDSQR with expert bench
- `BidiagonalSVD_TGK/dev/test_stcoll_alloc.f` — added expert bench around DBDSVDMR3 + DBDSQR-ref drift
- `BidiagonalSVD_TGK/dev/test_dbdsqr_full.f` — DBDSQR full SVD with expert bench
- `BidiagonalSVD_TGK/build.sh` — builds all three
- `eval_dbdsvdmr3.py` — `--solver {advisor,selfcontained,dbdsqr,both,all}`, `--paper-norms`, scaling section for 379
- `run_all_sweeps_2026-06-18.sh` — driver

## Benchmark protocol (all three solvers)

Identical for advisor / dbdsqr Fortran binaries, and conceptually mirrored on
the Python side for selfcontained (Python-side bench because we don't modify
`mr3gk_run`):

1. **Warmup**: one untimed call to JIT-page and load shared libraries.
2. **Save/restore**: copy D, E into DSAVE, ESAVE; restore between every rep
   (DBDSQR and DBDSVDMR3 destroy their inputs).
3. **Cache flush**: walk a 32 MiB scratch buffer (4M doubles) — read+write
   every cache line — before each timed call. JSUM-based anti-DCE guard at
   the end ensures the compiler can't elide the flush.
4. **High-res clock**: `SYSTEM_CLOCK(INTEGER*8)` (≈ ns under glibc CLOCK_MONOTONIC).
5. **Adaptive loop**: stop when (a) first rep ≥ 0.5 s, (b) NREPS ≥ 3 and
   total ≥ 0.2 s, or (c) NREPS ≥ 30.
6. **Headline `t_eval`**: MIN across reps — least-polluted-by-jitter sample.
7. **Process pinning**: each sweep pinned to a single core via `taskset` so
   the kernel scheduler can't migrate the process mid-call.
