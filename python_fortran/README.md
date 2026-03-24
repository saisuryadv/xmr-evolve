# MR³-GK Bidiagonal SVD — Python + Fortran Implementation

O(n²) bidiagonal SVD using MR³ with Golub-Kahan structure preservation (Willems & Lang 2012).

## Architecture

```
Python (mr3_gk.py)          Fortran (XMR, ~15K lines)
┌─────────────────┐         ┌──────────────────┐
│ bidiag_svd()    │         │ dlaxre_gk.f      │ ← GK-enabled root repr
│   split_bidiag  │         │ dlaxrv.f         │ ← eigenvector computation
│   mr3_tgk       │──ctypes─│ dlaxrf.f + deps  │ ← child repr (shift/twist)
│   sign recovery │         │ xmr_wrapper.c    │ ← C interface
│   Bv recovery   │         └──────────────────┘
│   GS completion │
└─────────────────┘
Eigenvalues: LAPACK dstebz (bisection)
```

## Current Results

| Suite | Tests | Passed | Rate |
|-------|-------|--------|------|
| Full (90 patterns × 4 sizes + 19 STCollection) | 379 | 308 | 81.3% |
| Adversarial only (90 patterns × 4 sizes) | 360 | 299 | 83.1% |
| STCollection | 19 | 9 | 47.4% |

Pass averages (over passing tests): res=0.303 nε, oU=0.509 nε, oV=0.529 nε

### Baseline Comparison

All algorithms evaluated with the same Python evaluator (`run_baselines.py`): 2 warmup runs + median of 5 timed runs, `perf_counter`, pre-copied arrays, infinity-norm bnorm matching C++ evaluator.

| Algorithm | Pass | Rate | Avg Res | Avg oU | Avg oV | Worst Scaling | Score |
|-----------|------|------|---------|--------|--------|---------------|-------|
| **DBDSQR** | **371/371** | **100%** | 0.194 | 0.190 | 0.199 | 25.9× (O(n³)) | 5.00 (gated) |
| **MR³-GK (ours)** | **300/371** | **80.9%** | 0.303 | 0.509 | 0.529 | **4.90×** | **82.75** |
| TGK+DSTEMR | crash | — | — | — | — | — | — |

DBDSQR is the gold standard for accuracy (100% pass) but is O(n³), so the scaling gate caps it at 5.0. Our MR³-GK is the only algorithm with an ungated score. TGK+DSTEMR crashes due to memory corruption in LAPACK's dstemr (via f2c/ctypes). From the C++ evaluator, HGBSVD passes 172/379 and TGK+DSTEXR passes 147/379 — both well below our 300.

To reproduce:

```bash
python3 run_baselines.py
```

### Scoring

```
score = pass_rate×50 + max(0, 5-avg_res)×2 + max(0, 5-avg_oU)×2 + max(0, 5-avg_oV)×2 + 5 + scaling_score
```

- HARD GATE: if worst_ratio (t400/t200 for patterns passing all sizes) > 5.0, score capped at 5.0
- Current worst_ratio ≈ 4.9× (under 5.0), score: **82.75**

### Evaluators

There are two evaluators that can score this code:

1. **Python evaluator** (`evaluate.py`) — runs `bidiag_svd()` from `mr3_gk.py` directly. Uses 2 warmup runs + median of 5 timed runs for stable scaling measurement. This is what we use day-to-day.
2. **C++ evaluator** (`../src/evaluate.cpp`) — compiled against a `bidiag_svd.h` header. Designed for C++ implementations but can call our Python code via a pipe bridge (see below). Uses slightly different matrix generation (xorshift RNG vs numpy), which causes ~23 test differences. Single-run timing (no warmup), so scaling ratios are noisier.

Both use identical scoring formulas and thresholds (res≤7nε, oU≤5nε, oV≤5nε).

Results are saved to `../results/`:
- `mr3_gk_python.txt` — output from the Python evaluator
- `dbdsqr.txt`, `tgk_stemr.txt`, etc. — baselines from C++ evaluator with other algorithms

## Prerequisites

```bash
sudo apt install gfortran liblapack-dev libblas-dev python3-numpy
```

## Build

```bash
# Build libxmr.so from Fortran objects + modified dlaxre_gk.f
bash build.sh
```

A prebuilt `libxmr.so` is included for Linux x86_64.

## Run Evaluation

```bash
# Full 379-test suite (90 patterns × {10,100,200,400} + 19 STCollection)
# This is the definitive evaluation. Output includes per-test metrics and final score.
python3 evaluate.py

# Save results to the results directory:
python3 evaluate.py | tee ../results/mr3_gk_python.txt

# Small ~66-test suite (22 patterns × {10,50,100}) — fast iteration
python3 evaluate.py --small

# Medium 270-test suite (90 patterns × {10,100,200}) — no n=400 or STCollection
python3 evaluate.py --medium
```

### Running with C++ evaluator (pipe bridge)

To verify against the C++ evaluator's matrix generation and scoring:

```bash
# 1. Create the Python server script (bsvd_server.py) that reads d,e from stdin
#    and writes sigma,U,V to stdout in binary format.

# 2. Build evaluate.cpp with bidiag_svd_python.h (pipe bridge header)
cd ../src
sed 's|#include "bidiag_svd.h"|#include "bidiag_svd_python.h"|' evaluate.cpp > evaluate_py.cpp
g++ -O2 -std=c++17 -o evaluate_cpp evaluate_py.cpp -llapack -lblas -lgfortran -lm

# 3. Run (increase timeouts in evaluate_py.cpp since Python is slower)
./evaluate_cpp ../python_fortran/stcollection 200
```

### Output format

Each test prints:
```
  pattern_name    n=NNNN  res=X.XXX  ortU=X.XXX  ortV=X.XXX  t=X.XXXXs  PASS/FAIL
```

Thresholds: res ≤ 7.0 nε, ortU ≤ 5.0 nε, ortV ≤ 5.0 nε.

Final summary includes pass counts, averages, worst scaling ratio, and composite score.

## Files

| File | Description |
|------|-------------|
| `mr3_gk.py` | Main implementation: splitting, T_GK construction, SVD extraction |
| `xmr_ctypes.py` | Python ctypes interface to XMR Fortran library |
| `full_eval.py` | 90 adversarial test matrix generators |
| `evaluate.py` | Evaluation runner with scoring |
| `run_baselines.py` | Baseline comparison (DBDSQR, TGK+DSTEMR vs ours) |
| `dlaxre_gk.f` | Modified XMR root representation with GK detection enabled |
| `xmr_wrapper.c` | C wrapper for dlaxre_ + dlaxrv_ |
| `libxmr.so` | Prebuilt shared library (Linux x86_64) |
| `build.sh` | Build script to recompile from source |
| `fortran_objects/` | Precompiled XMR .o files (unmodified from Willems' code) |
| `stcollection/` | 19 STCollection bidiagonal test matrices (.dat files) |
| `CHANGES.md` | Full history of changes and design decisions |

## Key Design Decisions

- **No reorthogonalization**: All orthogonality comes from MR³ structure, not post-hoc correction
- **NEGL guard kept**: dlaxre's quality check routes some matrices to PD shift fallback — removing it causes regressions
- **Gram-Schmidt for zero σ only**: GS completion triggers at most once per unsplit block for zero singular values
- **Two-phase splitting**: Relative split for accuracy, absolute split per-block for condition (3.5)
