# Claude debugging handoff: why current DBDSVR loses to DBDSDC

Copy everything below this line into Claude Code while checked out on branch
`agent/dbdsvr-paper29-debug-handoff`.

---

You are debugging a LAPACK-style bidiagonal SVD implementation. Work as a
senior numerical-linear-algebra and single-core performance engineer. Do not
give a generic complexity explanation: inspect the Fortran, reproduce or use
the supplied measurements, profile the actual phases, and identify concrete,
line-level causes.

## Objective

Determine why the current `DBDSVR`/`DBDSVDMR3` is slower than Netlib `DBDSDC`
on most cases in our reconstructed LAW 166-style benchmark, even at
N = 3000, while being 6×–65× faster on a smaller set of matrix families.

Produce:

1. an evidence-ranked root-cause analysis;
2. a phase/operation profile comparing representative winners and losers;
3. the highest-impact actionable optimizations that preserve all observable
   solver behavior and bitwise output;
4. if a large safe optimization exists, implement it and rerun the exact c2
   benchmark before and after;
5. a clear separation between algorithmic disadvantages, implementation
   redundancies, wrapper overhead, and fallback overhead.

## Non-negotiable constraints

- Performance measurements must run on c2, not on a Mac.
- Use exactly one CPU: CPU request = limit = 1, pin the process to one allowed
  CPU, and set OMP/OpenBLAS/MKL/BLIS/vecLib thread counts to 1.
- Keep the repository's `-O2` optimization setting. Do not explain gains by
  switching to `-O3`, `-Ofast`, architecture-specific flags, a different BLAS,
  or multithreading.
- Focus on `RANGE = 'A'`: all singular values and both full singular-vector
  sets.
- Reuse `TESTING/dbdsv_benchmark.f90` and
  `TESTING/EIG/dbsvrmg_paper29.f`. Do not replace the evaluation with a new
  benchmark whose timing boundary or matrices differ.
- Any optimized solver must produce bitwise-identical `S` and `Z`, the same
  `INFO`, `MFOUND`, and effective fallback result as the current solver for the
  tested inputs. Validate with bytewise hashes, not only residual tolerances.
- Preserve the result-based safety policy. Do not obtain speed by silently
  disabling the nonfinite-output audit or DBDSDC fallback.
- Diagnostic instrumentation may be conditional and may print extra data, but
  it must be disabled for the final timing run.
- Do not count matrix generation, reference-spectrum setup, or cubic
  residual/orthogonality validation as solver time.
- Do not call the reconstructed suite the exact Figure 6.1 suite. LAW 166 did
  not publish its individual matrix manifest.

## High-level result already obtained

Define:

`speedup = DBDSDC seconds / DBDSVR seconds`

A value above 1 favors DBDSVR; below 1 favors DBDSDC.

| N | DBDSVR wins | DBDSDC wins | Median speedup | Minimum–maximum | Routed DBDSDC fallbacks |
|---:|---:|---:|---:|---:|---:|
| 500 | 12/29 | 17/29 | 0.428× | 0.070×–11.68× | 2 |
| 1000 | 12/29 | 17/29 | 0.338× | 0.058×–13.34× | 2 |
| 2000 | 13/29 | 16/29 | 0.441× | 0.070×–39.24× | 2 |
| 3000 | 13/29 | 16/29 | 0.537× | 0.073×–65.02× | 2 |
| All | 50/116 | 66/116 | 0.430× | 0.058×–65.02× | 8 |

At N = 3000, the median case makes routed DBDSVR about 1.86× slower than
DBDSDC. Excluding the two fallback families, native MR3 wins 13/27 cases and
has median speedup 0.798×. Therefore fallback hurts but does not explain the
whole result.

All 116 routed DBDSVR calls ended with `INFO = 0` and `MFOUND = N`.

## Strongly bimodal behavior

DBDSVR wins at all four sizes on these 12 families:

- types 3–5, DMATGEN 112–114: three spread prescribed spectra;
- types 13–16, IDs 200–203: ABCON0–ABCON3;
- type 24, ID 230: Clement;
- types 25–28, IDs 240–243: GRO0–GRO3.

Geometric type 6, ID 115, changes from a DBDSDC win to a DBDSVR win at
N = 2000 and N = 3000.

DBDSDC wins all four sizes on the other 16 families, including the explicit
clusters 118–121, random, graded, and Wilkinson families.

Representative N = 3000 timings:

| Type | Family | DBDSDC s | DBDSVR s | Speedup | Interpretation |
|---:|---|---:|---:|---:|---|
| 24 | Clement | 22.7893 | 0.3505 | 65.02× | huge DBDSVR win |
| 16 | ABCON3 | 23.4434 | 0.4625 | 50.68× | huge DBDSVR win |
| 13 | ABCON0 | 23.5422 | 0.5045 | 46.67× | huge DBDSVR win |
| 4 | uniform ε → 1 | 16.2032 | 0.3781 | 42.85× | huge DBDSVR win |
| 25 | GRO0 | 3.0015 | 0.0926 | 32.42× | huge DBDSVR win |
| 2 | uniform ε apart | 0.8180 | 1.0255 | 0.798× | near crossover |
| 8 | random spectrum | 0.2503 | 1.2012 | 0.208× | DBDSDC 4.8× faster |
| 9 | clustered at 1 | 0.1946 | 1.1670 | 0.167× | DBDSDC 6.0× faster |
| 11 | clustered at ε | 0.00913 | 0.06536 | 0.140× | DBDSDC 7.2× faster |
| 20 | Wilkinson+ | 0.1600 | 2.1946 | 0.0729× | DBDSDC 13.7× faster |
| 22 | Wilkinson W | 0.1614 | 0.4255 | 0.379× | MR3 output audit fallback |
| 23 | double Wilkinson W | 0.1343 | 0.8144 | 0.165× | MR3 output audit fallback |

This split is a central clue. Determine what DBDSDC deflation statistics and
what MR3 cluster/tree/refinement statistics distinguish these groups.

## The eight fallbacks

`DBDSVR` attempts `DBDSVDMR3` first for every matrix order. It falls back to
`DBDSDC` only when MR3 returns positive `INFO`, returns `MFOUND ≠ N`, or leaves
a NaN/infinity in requested `S`, `U`, or `VT` output. `IWORK(2)` reports:

- 1: native DBDSVDMR3 accepted;
- 2: positive-INFO fallback to DBDSDC;
- 3: output-audit fallback to DBDSDC.

The audit found path 3 for exactly these eight calls:

- Wilkinson W, ID 224, at N = 500, 1000, 2000, and 3000;
- double Wilkinson W, ID 225, at the same four sizes.

No path-2 fallback occurred. Find the first nonfinite value, the output array
and index in which it appears, and the internal MR3 phase that creates it.
Determine whether these are all manifestations of one defect. Do not merely
state that the output audit caught the problem.

The routed timing for a fallback includes the failed MR3 attempt, the O(N²)
output scan, input restoration, and DBDSDC retry. If there is a provably safe
way to detect the same inevitable failure earlier and still return the exact
same DBDSDC bits and documented path, quantify it. Do not hard-code benchmark
type numbers.

## Exact benchmark environment

The completed run used:

- c2 Kubernetes context;
- AMD EPYC 9454 48-Core Processor;
- pod CPU request = limit = 1;
- cgroup `cpu.max = 100000 100000`;
- `taskset -c 0` inside the pod;
- 8 GiB memory request and limit;
- GNU Fortran 12.3.0;
- unchanged benchmark flags `-O2 -g`;
- Netlib LAPACK 3.11.0 at commit
  `7866626840f5d5e7e27f027a55182da8b3303872`;
- Netlib reference serial BLAS;
- all library thread controls set to 1;
- adaptive timing target ≥0.5 seconds per case, capped at 64 repetitions;
- one untimed solver call before the timed repetition batch.

The benchmark was compiled separately for each N so `BENCH_MAXN = N`; this
avoids introducing a larger leading dimension for smaller matrices.

The paper's environment was very different: Pentium 4 at 2.8 GHz with 512 KiB
cache, Intel Fortran 8.1 with `-O3 -tpp7 -mp`, and LAPACK 3.0. It does not state
the OS, BLAS implementation, thread count, Figure 6.1 repetition policy, or the
29 matrix formulas. See:

https://www.netlib.org/lapack/lawnspdf/lawn166.pdf

Treat compiler/platform differences as possible constant-factor contributors,
not as an adequate root-cause conclusion.

## What the paper actually says about its matrices

Figure 6.1 has 116 plotted cases in four blocks N = 500, 1000, 2000, 3000, so
29 cases per size is inferred. LAW 166 says the matrices were robustness tests,
many having very tight and sometimes large clusters of singular values. It says
tight clusters can favor DBDSDC through heavy deflation and can force bidiagonal
MRRR down multiple representation-tree levels. It does not enumerate the
individual cases and defers details to an unpublished 2005 reference.

The supplied proxy is the best recoverable direct-bSVD reconstruction: the
29-case Groesser–Lang DMATGEN family, IDs 110–121 and 200–244. It contains:

- 12 prescribed-spectrum cases: ones, near-ε uniform, ε¹ᐟ⁴-spaced, arithmetic,
  signed arithmetic, geometric, signed geometric, random, and four clustered
  distributions near ±1 or ±ε;
- 17 entry-defined cases: ABCON0–3, random bidiagonal, GRADP/GRADM, Wilkinson
  variants, Clement, GRO0–3, and bidiagonal Wilkinson+.

Spectrum-defined and symmetric-tridiagonal inputs are converted to upper
bidiagonal form by a Gershgorin-shifted Cholesky lift. Direct bidiagonals are
passed through unchanged. The exact manifest is supplied in the results.

## Timing boundary and fairness details

Read `TESTING/dbdsv_benchmark.f90` before drawing conclusions.

- Matrix generation occurs outside the timing interval.
- The driver calls a DBDSDC reference-spectrum helper before the measured
  method, but this is outside the reported time.
- `solve_case` starts and stops `system_clock` around the requested solver call.
- `TIMING_ONLY` reports the sum of those native solver intervals divided by the
  repetition count, not total executable time.
- The initial untimed solver call is excluded from the reported average.
- For DBDSVR, input copies, MR3, result audit, RANGE packing, and any fallback
  are intentionally inside the solver call because they are part of DBDSVR.
- Residual and orthogonality calculations are disabled in `TIMING_ONLY`; they
  are cubic at these sizes and caused misleading historical executable
  timeouts.
- DBDSDC and DBDSVR run as separate processes, but both begin with the same
  fixed LAPACK seed and traverse the same generator path, producing paired
  inputs.

Check for cache-state or repetition bias, but preserve this primary benchmark
when reporting before/after results. Any alternate microbenchmark must be
clearly secondary.

## Files to inspect first

Implementation:

- `python_fortran/lapack_pr_dbdsvr/SRC/dbdsvr.f`
  - LAPACK-style wrapper, copies, O(N²) output audit, fallback, ordering and
    RANGE packing.
- `python_fortran/lapack_pr_dbdsvr/SRC/dbdsvdmr3.f`
  - MR3 bidiagonal SVD engine entry point.
- `python_fortran/lapack_pr_dbdsvr/SRC/dbdtgk.f`
  - Golub–Kahan representation construction.
- `python_fortran/lapack_pr_dbdsvr/SRC/dlarrv_tgk.f`
  - representation tree, cluster handling, vector construction.
- `python_fortran/lapack_pr_dbdsvr/SRC/dlar1v_tgk.f`
  - twisted-factorization recurrence and vector work.
- `python_fortran/lapack_pr_dbdsvr/SRC/dlarrf_tgk.f`
  - child representations.
- `python_fortran/lapack_pr_dbdsvr/SRC/stegr_ID/dlarrb.f`
  - eigenvalue refinement and bounded progress checks.

Benchmark and generator:

- `python_fortran/lapack_pr_dbdsvr/TESTING/dbdsv_benchmark.f90`
- `python_fortran/lapack_pr_dbdsvr/TESTING/EIG/dbsvrmg_paper29.f`
- `python_fortran/lapack_pr_dbdsvr/Makefile`
- `python_fortran/scripts/run_dbdsvr_paper29_c2.sh`
- `python_fortran/scripts/plot_dbdsvr_figure61.py`

Evidence:

- `python_fortran/results/dbdsvr_paper29_c2_1cpu_2026-08-07/README.md`
- `python_fortran/results/dbdsvr_paper29_c2_1cpu_2026-08-07/figure61_points.csv`
- `python_fortran/results/dbdsvr_paper29_c2_1cpu_2026-08-07/summary.json`
- `python_fortran/results/dbdsvr_paper29_c2_1cpu_2026-08-07/matrix_manifest.csv`
- `python_fortran/results/dbdsvr_paper29_c2_1cpu_2026-08-07/raw/`
- `python_fortran/results/dbdsvr_paper29_c2_1cpu_2026-08-07/backend_audit_n*.csv`
- `knowledge/grosser_lang_2001_hgbsvd.md`

## Reproduction commands

Build Reference LAPACK 3.11.0 and the bundle as documented in
`python_fortran/lapack_pr_dbdsvr/build.sh`. For a specific benchmark size,
clean first because Make does not track a changed preprocessor macro:

```sh
cd python_fortran/lapack_pr_dbdsvr
make clean
make -j1 FC=gfortran-12 BENCH_MAXN=3000 paper29-benchmark
```

Then run one method on the pinned CPU:

```sh
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

taskset -c 0 ./build/dbdsv_paper29_benchmark \
  0 3000 A DBDSVR 0.5 64 "" TIMING_ONLY

taskset -c 0 ./build/dbdsv_paper29_benchmark \
  0 3000 A DBDSDC 0.5 64 "" TIMING_ONLY
```

Arguments are `type_filter`, `N`, `RANGE`, `method`, timing target, repetition
cap, optional dump path, and `TIMING_ONLY`. Type filter 0 runs all 29 cases.

`backend_path` is emitted by the current benchmark driver. A cheap one-call
route audit uses timing target `0.000001` and repetition cap `1`.

## Required investigation

### 1. Verify the baseline

- Reproduce at least representative N = 3000 winners and losers on c2.
- Confirm exact compiler, flags, CPU quota/affinity, BLAS, and source hashes.
- Confirm paired matrices are byte-identical across methods.
- Confirm the timing boundary with source inspection and a small external
  sanity check.

### 2. Add phase counters and timers

For each representative matrix, measure at minimum:

- wrapper input copies and output packing;
- DBDSVDMR3 setup and scaling;
- singular-value computation/refinement;
- cluster discovery and representation-tree construction;
- number and sizes of clusters by level;
- maximum and total representation-tree depth;
- calls and iterations in `DLARRB`;
- calls/retries in `DLARRF_TGK`;
- calls and work in `DLAR1V_TGK`;
- left/right singular-vector recurrence work;
- output audit scan;
- fallback restoration and DBDSDC retry where applicable.

For DBDSDC, capture split/deflation behavior and the size of surviving secular
subproblems if practicable. We need evidence that explains why the same N can
take DBDSDC 0.009 seconds on one case and 23 seconds on another.

Use Linux `perf` counters if permitted, but do not depend exclusively on them.
Conditional Fortran counters are acceptable. Preserve `-O2` for comparative
timings.

### 3. Find true redundancies

Look for large repeated O(N²) work, repeated refinements, recomputation of the
same recurrence data, avoidable full-vector sweeps, poor leading-dimension
access, or transformations that can be reused without changing floating-point
operation order. Quantify each candidate in seconds and percent of runtime on
both a fast DBDSVR case and a slow one.

Do not recommend tiny changes as the main answer. A proposed optimization must
either remove a dominant measured phase or explain why no large
bitwise-preserving optimization exists because the loss is algorithmic.

The O(N²) output audit and wrapper copies are obvious suspects, but measure
them. Removing the audit is forbidden. A bitwise-safe fusion or reuse is
acceptable only if it preserves detection and documented outputs.

### 4. Diagnose the nonfinite Wilkinson outputs

Instrument the native MR3 result before fallback. Record:

- first nonfinite array and index;
- singular value and cluster/representation identifiers;
- tree depth;
- pivots, shifts, gaps, and recurrence quantities immediately before the
  nonfinite value appears;
- whether the failure is deterministic across compiler runs;
- whether type 224 and 225 share the same mechanism.

Explain the numerical root cause, not just the final NaN/Inf symptom.

### 5. Compare with LAW 166 carefully

The paper timed `xBDSCR`, not this exact executable, against LAPACK 3.0. Inspect
whether the current code contains extra work or lacks an implementation detail
described by LAW 166. State what is demonstrated from source and what remains
speculation. Do not infer the unpublished 29-case manifest from the plot alone.

### 6. Validate any change

For any proposed patch:

- rerun the 116 c2 timing points under the identical protocol;
- show per-N wins and median/min/max speedup before and after;
- compare `S`, packed `Z`, `INFO`, `MFOUND`, and backend outcome byte for byte;
- run the existing fallback test and canonical email-derived test suite;
- verify Wilkinson W and double-W still return the same successful DBDSDC
  fallback results unless you prove an earlier equivalent route;
- report any case that changes bits as a failed constraint, even if residuals
  improve.

## Expected final answer

Lead with a high-level conclusion in plain language. Then provide:

1. a ranked root-cause table with measured evidence;
2. a winner-versus-loser phase comparison;
3. the precise Wilkinson nonfinite-output cause;
4. actionable optimization candidates, expected speedup, risk, and whether
   bitwise identity is possible;
5. before/after c2 results if code was changed;
6. exact commands and artifact paths so another engineer can reproduce it.

If the evidence says DBDSDC's advantage is fundamentally deflation-driven and
cannot be beaten without selecting a different algorithm or changing output
bits, say that explicitly. Do not manufacture an implementation optimization
to satisfy the premise.
