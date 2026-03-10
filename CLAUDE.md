# XMR-Evolve: Bidiagonal SVD via MR³

## What This Is
An OpenEvolve project evolving a C++ bidiagonal SVD algorithm. The target file
is `src/bidiag_svd.h`. The evaluator compiles it with `src/evaluate.cpp` and
runs 270 adversarial tests + scaling tests.

## The Goal
Find a **strictly O(n²)** bidiagonal SVD that passes **all 270 tests** (matching
DBDSQR's 270/270 accuracy but with O(n²) scaling instead of O(n³)).
**This is an open problem** — no known implementation achieves this fully.

## Project Structure
- `src/bidiag_svd.h` — **THE FILE BEING EVOLVED** (self-contained hybrid HGBSVD+DSTEMR)
- `src/evaluate.cpp` — Evaluation driver (DO NOT MODIFY)
- `src/bidiag_dbdsqr.h` — DBDSQR reference (270/270, O(n³))
- `src/bidiag_hgbsvd.h` — HGBSVD reference (174/379, O(n²))
- `src/bidiag_tgk_stemr.h` — TGK+STEMR reference (93/379, O(n²))
- `src/bidiag_tgk_stexr.h` — TGK+DSTEXR reference (61/379, O(n²))
- `src/bidiag_tgk_common.h` — Shared TGK utilities
- `src/fortran_interface.h` — C++ declarations for Fortran routines
- `evaluator.py` — OpenEvolve Python evaluator (saves variants to `evolved_variants/`)
- `claude_agent_llm.py` — Claude Agent SDK adapter (single-phase, directory-scoped writes)
- `run_xmr.py` — OpenEvolve entry point
- `lib/libxmr_c.a` — Pure C library (CLAPACK + XMR + hgbsvd, no Fortran)
- `config.yaml` — OpenEvolve configuration
- `scratch/` — Agent workspaces (scratch/agent_NNN/ per call, gitignored)

## Evolved File Requirements
- `bidiag_svd.h` MUST be **self-contained** — no `#include` of project headers
  (OpenEvolve copies only this file to a temp dir for evaluation)
- Only standard library + extern "C" LAPACK/BLAS declarations allowed
- `dbdsgr_` (f2c) requires `ftnlen` args: pass `1, 1` at end of each call

## Knowledge Base (READ THESE FIRST)
- `knowledge/INDEX.md` — **Master reference**: Algorithm 4.1, NCD, GK structure, ALL known bugs (#1-#10), test matrix formulas, numerical thresholds
- `knowledge/BASELINES.md` — **Comparison table** of 5 baseline algorithms with pass rates, accuracy, scaling
- `knowledge/EVALUATION.md` — **How evaluation works**: metrics, scoring, all 90 test patterns, data layout, available routines
- `knowledge/PRIOR_APPROACHES.md` — **What's been tried**: 12 approaches (A-L), key lessons, why each failed
- `knowledge/PROGRESSION.md` — Paper-by-paper timeline (1990-2020)
- `knowledge/RESOURCES.md` — **Index** of all MRRR Resources: 15 papers, code trees, STCollection, presentations
- `knowledge/willems_lang_2012.md` — Algorithm 4.1 (MR³ on TGK with NCD shifts) — the target algorithm
- `knowledge/willems_lang_2013_framework.md` — Five-requirement MR³ framework
- `knowledge/xmr_code_documentation.md` — Willems XMR 45-file Fortran codebase documentation
- `knowledge/grosser_lang_2001_hgbsvd.md` — Coupling-based O(n²) SVD approach
- `knowledge/mr3_foundations.md` — Dhillon 1997, Parlett-Dhillon 2000 (RRR theory)
- `knowledge/marques_2020_bugs_matrices.md` — DBDSVDX bugs, CHKBD analysis
- `knowledge/demmel_2008_test_matrices.md` — Demmel 2008 test matrix taxonomy

## Build & Test
```bash
# Pure C/C++ build — no Fortran, no gfortran, no Apple Accelerate
g++ -std=c++17 -O2 -Isrc/clapack -o evaluate src/evaluate.cpp -Llib -lxmr_c -lm
./evaluate       # Run all 270 tests (no STCollection dir = adversarial only)
```

## Scoring (0-100)
- 5 pts: compiles
- 50 pts: test pass rate (270 tests, threshold: res ≤ 7, orthoU ≤ 5, orthoV ≤ 5)
- 15 pts: accuracy (residual + orthogonality)
- 10 pts: O(n²) scaling (worst doubling ratio ≤ 5.0)
- 10 pts: large adversarial pass rate
- 10 pts: bonus for very good accuracy (all metrics < 2)

## Current Baselines (Pure C/C++ — No Fortran)
90 patterns × 4 sizes (n=10,100,200,400) + 19 STCollection = 379 tests (stage 3).
Pass: res ≤ 7.0, orthoU ≤ 5.0, orthoV ≤ 5.0 (n·eps). Avg over all tests (capped at 1e6).
Scaling = worst time doubling ratio (n=200→400) across patterns that pass at ALL sizes.

| Algorithm | Pass | Avg Res | Avg OrtU | Avg OrtV | Worst Pass Scaling | Score |
|-----------|------|---------|----------|----------|--------------------|-------|
| DBDSQR | 379/379 | 0.43 | 0.37 | 0.35 | 8.53x (O(n³)) | 85.64 |
| DBDSVDX | 371/379 | 13193 | 15852 | 13213 | 13.21x (O(n³)) | 53.94 |
| HGBSVD | 174/379 | 533121 | 532982 | 533112 | 4.35x (O(n²)) | 37.96 |
| TGK+STEMR | 93/379 | 52935 | 492760 | 451378 | 4.82x (O(n²)) | 27.27 |
| TGK+STEXR | 61/379 | 163642 | 567171 | 520277 | 5.73x | 21.58 |

## Agent Architecture
Single-phase agent with full tools and isolated per-call workspace:
- **Workspace**: `scratch/agent_NNN/` — isolated directory per agent call
  - `program/bidiag_svd.h` — the output file (pre-copied, agent edits it, read back as result)
  - `evaluate.cpp` — pre-copied for compilation
- **Write scope**: only `scratch/agent_NNN/` (enforced via `can_use_tool` callback)
- **Read scope**: unrestricted (all project files, knowledge base, references)
- **Tools**: Read, Write, Edit, Bash, Grep, Glob, WebSearch, WebFetch, Agent
- **Subagents**: `numerical-analyst` (stability analysis), `bug-hunter` (known bugs)
- **Result capture**: file-based (read `program/bidiag_svd.h` from disk, not text extraction)
- **Budget**: $3.50/call default (configurable via `--budget`)

## Agent Autonomy Rules
- **NEVER stop to ask whether an approach or direction is correct.** Just try it, evaluate it, and report results.
- **NEVER ask for permission or confirmation** before implementing, compiling, running tests, or changing `bidiag_svd.h`. Just do it.
- If an approach fails, **diagnose why it failed** (which matrices, which metrics, what numerical quantities are off), then try a fix or a different approach. Do not stop and ask the user what to do next.
- If you're unsure between two approaches, **try both** and compare results. Don't ask which one to try.
- The only reason to stop and ask is if you need information that isn't in the codebase or knowledge files (e.g., a missing dependency, a broken build environment).
- **Be relentless**: compile, test, analyze failures, iterate. The evaluation framework gives you everything you need to judge whether something works.

## Critical Rules
- **O(n²) worst-case required** — no global MGS, no dense n×n multiplies
- The fundamental tradeoff: accuracy vs scaling. DBDSQR has accuracy but O(n³). TGK+STEMR has O(n²) but fails 76% of tests.
- TGK eigenvector extraction (U from odd rows, V from even rows) destroys orthogonality unless eigenvectors have GK structure (requires NCD-aware shifts inside the eigensolver)
- One-sided recovery (U = BV/σ) helps for small σ but doesn't fix extraction failures
- Chunked MGS (chunk size 32) is O(n²) but may miss cross-chunk orthogonality
- Read `knowledge/PRIOR_APPROACHES.md` before proposing new ideas — 12 approaches already tested

## Key Reference Code
- XMR Fortran: `../bidiag-algo/MRRR Resources/Code/tdsolver/xmr/SRC/O/`
- hgbsvd Fortran: `../bidiag-algo/MRRR Resources/Code/hgbsvd/hgbsvd/v2/`
- LAPACK source: `lapack/SRC/` (dbdsqr.f, dbdsvdx.f, dstemr.f, dlarrv.f, dstein.f)
- STCollection: `../bidiag-algo/MRRR Resources/STCollection/DATA/`
- Papers: `../bidiag-algo/MRRR Resources/Papers/`
- Prior Python approaches: `../bidiag-algo/approach_*.py` (esp. approach_k_ncd.py — best result)
