# DBDSVR versus DBDSDC on the reconstructed 29-case suite

## Result

On this controlled c2 run, DBDSDC wins more reconstructed cases at every
size.  DBDSVR nevertheless has very large wins on 12 matrix families, reaching
65.0× at N = 3000.  The result is therefore strongly bimodal rather than a
uniform advantage for either routine.

Speedup below means:

`speedup = DBDSDC seconds / DBDSVR seconds`

Thus, values above 1 favor DBDSVR and values below 1 favor DBDSDC.

| N | DBDSVR wins | DBDSDC wins | Median speedup | Range | DBDSDC fallbacks inside DBDSVR |
|---:|---:|---:|---:|---:|---:|
| 500 | 12/29 | 17/29 | 0.428× | 0.070×–11.68× | 2 |
| 1000 | 12/29 | 17/29 | 0.338× | 0.058×–13.34× | 2 |
| 2000 | 13/29 | 16/29 | 0.441× | 0.070×–39.24× | 2 |
| 3000 | 13/29 | 16/29 | 0.537× | 0.073×–65.02× | 2 |
| **All** | **50/116** | **66/116** | **0.430×** | **0.058×–65.02×** | **8** |

At N = 3000, the median ratio 0.537× means that DBDSVR takes about 1.86×
as long as DBDSDC on the median case.  This is not caused only by fallback:
after excluding the two fallback families, DBDSVR wins 13/27 native-MR3 cases
and the median ratio is 0.798×.

## Which cases favor each algorithm

DBDSVR wins at all four sizes on these 12 families:

- 112–114: three spread prescribed-spectrum cases
- 200–203: ABCON0–ABCON3
- 230: Clement
- 240–243: GRO0–GRO3

At N = 3000, its largest wins are Clement (65.0×), ABCON0–ABCON3
(44.1×–50.7×), the GRO family (31.3×–32.4×), and prescribed-spectrum cases
112–114 (6.2×–42.8×).  Geometric case 115 begins favoring DBDSVR only at
N = 2000 and N = 3000.

DBDSDC wins at all four sizes on the remaining 16 families.  In particular,
the explicitly clustered cases 118–121 favor DBDSDC by about 2.4×–7.2× at
N = 3000.  This agrees with LAW 166's explanation: heavy deflation helps
divide-and-conquer, while tight clusters can force MRRR through more levels of
its representation tree.

The DBDSVR path audit found eight successful result-based fallbacks:

- type 224, Wilkinson W, at N = 500, 1000, 2000, and 3000
- type 225, double Wilkinson W, at N = 500, 1000, 2000, and 3000

For each, DBDSVDMR3 returned nonfinite requested output, the output audit
selected path 3, and DBDSDC completed successfully.  These timings correctly
measure the current routed DBDSVR call, including the failed MR3 attempt,
audit, and DBDSDC retry.  Every one of the 116 measured DBDSVR calls ended with
INFO = 0 and MFOUND = N.

## Why LAW 166 looks more favorable

The comparison cannot be made as an exact reproduction:

1. LAW 166 does not enumerate the Figure 6.1 input matrices.  It shows 116
   points in four size blocks, implying 29 cases per size, and describes them
   only as robustness matrices, many with tight and sometimes large clusters.
   The authors say the environment was adapted from Osni Marques's xSTEGR test
   suite and defer details to an unpublished 2005 reference.
2. This run therefore uses the best recoverable direct-bSVD reconstruction:
   the 29-case Groesser–Lang DMATGEN family, IDs 110–121 and 200–244.  It is a
   source-grounded proxy, not proof of the unpublished Figure 6.1 manifest.
3. LAW 166 timed its xBDSCR implementation against LAPACK 3.0 DBDSDC.  This run
   times the current DBDSVR/DBDSVDMR3 code against LAPACK 3.11.0 DBDSDC.  These
   are related MRRR implementations, but they are not the same executable.
4. The hardware and compiler are radically different.  The old paper used a
   Pentium 4 and Intel Fortran 8.1 at `-O3 -tpp7 -mp`; this run uses an AMD EPYC
   9454 and GNU Fortran 12.3 at the repository's unchanged `-O2` setting.

The strongest explanation is therefore suite and implementation sensitivity,
with platform/compiler constants as an additional source of movement.  The new
data actually reproduce the paper's qualitative observation: MRRR can win by a
large and growing factor on some distributions, while deflation-friendly or
tightly clustered matrices still strongly favor DBDSDC.  What does not
reproduce is the paper's claim that MRRR wins *most* cases by N = 3000.

## Reconstructed varieties

The full formulas and construction paths are in
[`matrix_manifest.csv`](matrix_manifest.csv).  The 29 cases are:

- prescribed spectra 110–121: ones; near-ε uniform; ε¹ᐟ⁴-spaced; arithmetic;
  signed arithmetic; geometric; signed geometric; random; and four clustered
  spectra near ±1 or ±ε
- entry-defined 200–244: ABCON0–3; random bidiagonal; GRADP/GRADM;
  Wilkinson+, Wilkinson−, W, and double-W; Clement; GRO0–3; and bidiagonal
  Wilkinson+

Spectrum-defined and symmetric-tridiagonal cases are reduced and then lifted
to an upper bidiagonal B using a Gershgorin-shifted Cholesky factorization.
Direct bidiagonal cases are passed through unchanged.  Every run uses the same
fixed LAPACK random seed.

## Controlled environment and timing method

- Kubernetes context: c2
- node CPU: AMD EPYC 9454 48-Core Processor
- pod resources: CPU request = limit = 1; cgroup `cpu.max = 100000 100000`
- process affinity: CPU 0 only
- memory request = limit = 8 GiB
- compiler: GNU Fortran 12.3.0
- benchmark flags: `-O2 -g`; no optimization-level change
- baseline: Netlib LAPACK 3.11.0, commit
  `7866626840f5d5e7e27f027a55182da8b3303872`, with reference serial BLAS
- thread controls: OMP, OpenBLAS, MKL, BLIS, and vecLib thread counts all 1
- operation: RANGE = A, all singular values and both full singular-vector sets
- sizes: N = 500, 1000, 2000, 3000
- timing interval: native solver call only; generation and metrics excluded
- adaptive repetitions: at least 0.5 seconds per case, capped at 64 repeats

LAW 166 reports a Pentium 4 at 2.8 GHz with 512 KiB cache, Intel Fortran 8.1,
and flags `-O3 -tpp7 -mp`.  It does not report the operating system, BLAS
implementation, thread count, or Figure 6.1 repetition policy.

## Artifacts

- [`dbdsvr_vs_dbdsdc_figure61_paper29.png`](dbdsvr_vs_dbdsdc_figure61_paper29.png): Figure 6.1-style plot; orange X marks a DBDSDC fallback
- [`figure61_points.csv`](figure61_points.csv): all 116 paired ratios and backend paths
- [`summary.json`](summary.json): machine-readable aggregate results
- [`raw/`](raw/): eight raw timing CSVs, one per size and method
- [`backend_audit_n500.csv`](backend_audit_n500.csv),
  [`backend_audit_n1000.csv`](backend_audit_n1000.csv),
  [`backend_audit_n2000.csv`](backend_audit_n2000.csv), and
  [`backend_audit_n3000.csv`](backend_audit_n3000.csv): one-call route audit
- [`source_sha256.txt`](source_sha256.txt): exact timed-source hashes
- [`backend_audit_source_sha256.txt`](backend_audit_source_sha256.txt): exact route-audit source hashes
- [`build_logs/`](build_logs/): one build log per matrix size

## Sources

- [LAPACK Working Note 166](https://www.netlib.org/lapack/lawnspdf/lawn166.pdf)
- [LAPACK Working Note 183](https://www.netlib.org/lapack/lawnspdf/lawn183.pdf)
