# How to reproduce the python_fortran MR3-GK SVD with an LLM agent

This document gives you the **prompt sequence** that, when fed to a Claude (or
similar) coding agent in order, will evolve a starting state into the final
working 379/379 implementation. It also documents what initial code we built
on top of (Paul Willems' XMR Fortran library + the seed Python wrapper).

There are two checkpoints you can start from:

- **Checkpoint A** — full reproduction from upstream Willems XMR + the test
  scaffolding (no Python implementation at all). 9 prompts.
- **Checkpoint B** — start from the seed commit `c6d73b2` (300/379 already
  passing) and replay only the 300 → 379 transition. 6 prompts. Faster and
  cheaper to validate; primary recipe.

## What's in the initial state

### Checkpoint A — upstream Willems XMR (~early March 2026)

You start from these files (none of which we wrote):

```
python_fortran/
├── xmr_src/                  # Paul Willems' XMR Fortran source
│   ├── dlaxre.f              #   root representation builder
│   ├── dlaxrf.f              #   representation-tree forest
│   ├── dlaxrv.f              #   eigenvector solver
│   ├── dlaxrb_clssfy.f       #   cluster classifier
│   ├── dlaxr{r,s,t,n,m}*.f   #   qd-style transformation kernels
│   └── … (44 files total, ~14 K LOC; from tdsolver/xmr/SRC/O)
├── stcollection/             # 19 STCollection real-world bidiagonals
│   └── B_*.dat
├── full_eval.py              # 90 adversarial pattern generators
└── evaluate.py               # 379-test scoring harness
                              # (res ≤ 7 nε, ortU/V ≤ 5 nε)
```

XMR is **eigenvalue-only**: no SVD entry point, no T_GK construction, no
preprocessing or post-processing for bidiagonals. You will build all of that.

XMR has two latent bugs that block bidiagonal SVD entirely:

1. **`dlaxre.f` line 269**: GK-detection branch is `IF(.FALSE.)` — disabled.
2. **`dlaxrb_clssfy.f` lines ~358, 414, 445**: classifier uses
   `MIN(AVGTHRESH, ABSMAX·GAPTHRESH)`, but at deep tree nodes
   `AVGTHRESH = 2·AVGTOL · 2^(−DEPTH)` underflows to 0, making the test
   trivially true and over-fragmenting clusters.

You will fix both.

### Checkpoint B — seed at commit `c6d73b2` (300/371)

You start from the python_fortran/ first commit, which already contains:

- `libxmr.so` built from upstream `xmr_src/` (with our XMR Bug #1 fix
  patched in — `dlaxre.f` → `dlaxre_gk.f`).
- `xmr_wrapper.c` C glue + `xmr_ctypes.py` ctypes binding.
- `mr3_gk.py` Python orchestration (~640 lines): T_GK construction, splitting,
  sign matrices D1/D2, Bv recovery, basic Gram-Schmidt completion.
- `evaluate.py`, `full_eval.py`, the 19 STCollection files.

At this checkpoint, `python3 evaluate.py` reports **300/371** passing.
The "/371" (rather than "/379") is a **regex bug** in the seed's
`get_patterns()` — it parses `full_eval.py` source with
`r"if\s+name\s*==\s*'(\w+)'"`, which catches 88 of the 90 patterns
(misses `chkbd_4`, `chkbd_16` that come from a parameterized loop).
4 sizes × 88 patterns + 19 STCollection = 371. The fix in prompt 1
switches to `importlib` to read `adv_names` directly, recovering all
4 × 90 + 19 = 379.

`evaluate.py` also reports a **SCORE** (pass-rate × O(n²) timing-ratio
gate). At the seed it's hard-gated at `SCORE = 5.00` by timing
outliers. Subsequent prompts both improve the pass count *and* lift
the SCORE gate.

The remaining 71 failures fall into three groups:

- ~~Singleton-block timeout~~ (~70 patterns time out in O(n³) GS)
- Tight clusters (`demmel_S1pe_k4`, `pd_T0`, `gl_abcon3`, …)
- Zero singular values from product underflow (`saw_tooth`, `step_function`)

The 6 prompts below take you from 300/379 to **379/379**.

---

## Prompt sequence — Checkpoint B (300 → 379)

Each prompt is a single message to the agent. Use **a fresh agent session per
prompt** so each is self-contained — the agent has no memory of prior prompts,
only the current state of the working tree.

The agent should have read/write tooling on `python_fortran/` and shell access
to run `python3 evaluate.py`, `bash build.sh`, and Fortran compilation.

### Prompt 1 — get_patterns regex + singleton-block bypass + timing robustness (300/371 → ~310/379)

```
Run `python3 evaluate.py` and observe the current score
(you should see TOTAL: 300/371 — the suite is documented as 379 but
only 371 tests run because of a regex bug; we'll fix that here too).

Three changes in this single prompt:

(a) get_patterns() in evaluate.py uses the regex
    r"if\s+name\s*==\s*'(\w+)'" to extract pattern names from
    full_eval.py source. This misses chkbd_4 and chkbd_16 (generated
    from a parameterized loop, not literal `==` comparisons). Replace
    the regex parsing with importlib: load full_eval.py as a module
    and read its `adv_names` list directly. After this fix the suite
    grows from 88 × 4 + 19 = 371 to 90 × 4 + 19 = 379.

(b) Singleton-block bypass in mr3_gk.py:mr3_tgk. When the bidiagonal
    splitter produces a singleton block (k=1), the code still calls
    XMR on a trivial 2x2 T_GK matrix. XMR returns degenerate
    eigenvectors with ||v||=1, ||u||=0 (only even rows populated).
    The post-processing then runs O(n^3) Gram-Schmidt on hundreds
    of such columns, causing chkbd_* and saw_tooth patterns at
    n >= 200 to time out.

    Add a fast path in mr3_tgk that handles singleton blocks
    directly, bypassing the XMR call:
       - sigma = |d[i]|
       - both u and v components are 1/sqrt(2) at the corresponding row
       - apply sign matrices D1, D2 (constructed earlier in
         bidiag_svd) to recover correct signs

(c) Also add timing robustness to test_one(): two warmup runs (JIT,
    cache effects), then 5 timed runs and take the median. This
    stabilizes the SCORE gate (which is sensitive to single-run
    timing noise).

Do NOT modify dlaxre_gk.f or any Fortran file. Re-run evaluate.py.
You should see TOTAL: ~310/379 (a few specs that newly appeared from
the get_patterns fix may fail at first) and SCORE much higher than
5.00 (the timing gate stops capping; expect ~80s).
```

### Prompt 2 — Bv recovery threshold + GS vectorization (~330 → ~340)

```
Some matrices show one-sided eigenvector quality issues — only u or
only v has the correct norm. Look at gl_gradm@200 (ortU about 33 in
the current state) and chkbd@200.

Two fixes in mr3_gk.py:
(a) The Bv recovery threshold is hardcoded; make it scale with
    sigma_max so it adapts to matrix dynamic range.
(b) The Gram-Schmidt completion runs one column at a time in a Python
    loop. Vectorize it: build a matrix of all candidate vectors and
    do the projection in a single numpy operation so each MGS pass
    is one matmul instead of n matvecs.

Verify on gl_gradm and chkbd — ortU should drop to ~1 nε. Don't
worry about other failures yet.
```

### Prompt 3 — T_GK sub-splitting + negative-e fix (~340 → ~340)

```
Two refinements in mr3_gk.py:bidiag_svd, near where T_GK is built:

(a) Zero diagonal entries: when the bidiagonal has d[i]=0, the
    corresponding rows of T_GK have zero off-diagonals on both sides
    — there's no eigenvalue coupling across that point. Sub-split
    T_GK at such positions so MR^3 doesn't waste effort on
    decoupled blocks.

(b) Negative e[i]: the perfect-shuffle T_GK construction gives
    negative off-diagonals when e[i] < 0, but XMR requires positive
    off-diagonals. Apply a diagonal-similarity sign normalization
    before forming T_GK (track the sign flips and re-apply them when
    extracting U, V).

These don't change the headline score by much but eliminate edge-case
failures. Verify the score doesn't regress.
```

### Prompt 4 — zero-σ deflation via zero-shift QR (~340 → ~370)

```
The Demmel–Kahan 1990 paper (Sec. 3) shows that one zero-shift QR
sweep on a bidiagonal is equivalent to one step of inverse iteration
on B^T B. When sigma_min = 0, BOTH d[-1] AND e[-1] are driven to
zero in a single sweep — clean detection of zero singular values.

Implement this as preprocessing in mr3_gk.py:bidiag_svd, applied to
each unsplit block of size >= 2 BEFORE the MR^3 dispatch:

1. Implement zero_shift_qr_sweep(d, e) using LAPACK dlartg:
   - Walk i = 0..n-2, generate a Givens rotation per step,
     update d, e in place,
     remember the rotations so we can undo them later.
   - Returns (right_rots, left_rots) — Givens-rotation lists.

2. After the sweep, check the clean-split criterion:
     |d[-1]| < n * eps  AND  |e[-1]| < n * eps
   If TRUE: zero singular value detected. Solve the (k-1) sub-problem
     recursively (on the leading (k-1)×(k-1) part), embed the
     solution into k×k with sigma=0 in the bottom-right, and apply
     the saved Givens rotations to recover the original basis.
   If FALSE: discard the swept (d, e), use the original block as
     input to MR^3 normally.

CRITICAL: the clean-split check must require BOTH d[-1] AND e[-1]
to be small. Don't deflate when only d[-1] is small — that severs
matrix coupling and produces catastrophic residuals (we measured
e[-1] = 0.707 after sweep on gl_wilkw, which would mean wrong
results if we deflated).

This fixes saw_tooth, step_function, gl_wilkp, gl_clement, and any
matrix with a structural zero diagonal. Score should jump to ~360.
```

### Prompt 5 — dlaxrb_clssfy AVGTHRESH=0 fix (~370 → ~378)

```
A subtle bug in upstream xmr_src/dlaxrb_clssfy.f. At deep tree
nodes, AVGTHRESH = 2·AVGTOL collapses to 0 because AVGTOL is scaled
by 2^(-DEPTH) and underflows. The classifier test is

   IF( SPREAD .GE. MIN(AVGTHRESH, ABSMAX*GAPTHRESH) ) gap is full

— with AVGTHRESH=0 this becomes SPREAD >= 0, trivially true, so
EVERY adjacent eigenvalue pair is labeled a full gap and the rep
tree degenerates into singletons. Tight-cluster matrices fail
catastrophically.

Steps:

1. Create python_fortran/dlaxrb_clssfy_fix.f as a copy of
   xmr_src/dlaxrb_clssfy.f.
2. Add a logical local variable GAPOK to the declarations.
3. At all THREE test sites in the file (around lines 358, 414, 445
   — search for "MIN( AVGTHRESH, ABSMAX*GAPTHRESH )"), replace
   each with:

       IF( AVGTHRESH .GT. ZERO )THEN
          GAPOK = SPREAD .GE. MIN(AVGTHRESH, ABSMAX*GAPTHRESH)
       ELSE
          GAPOK = SPREAD .GE. ABSMAX*GAPTHRESH
       ENDIF
       IF( GAPOK )THEN
          ! existing 'gap is full' code path
       ELSE
          ! existing 'no full gap' code path
       ENDIF

4. Modify build.sh to compile dlaxrb_clssfy_fix.f and link it in
   place of the upstream dlaxrb_clssfy.o (rename or override).
   Rebuild libxmr.so via `bash build.sh`.

5. Re-run evaluate.py. demmel_S1pe_k4, three_clusters, pd_T0,
   gl_wilkw, ST_B_20_graded should all flip from FAIL to PASS.
   Score should be ~378.
```

### Prompt 6 — adaptive GAPTOL for tiny blocks (378 → 379)

```
One last failure remains: two_clusters@10 (n=10). The cluster's
relative gap is right at the GAPTOL boundary, but at small n the
hard-coded GAPTOL = 1e-3 is too coarse — small absolute differences
look 'large enough' to satisfy the threshold even when they are not.

In mr3_gk.py:classify (the per-node gap classifier), replace every
use of the literal `GAPTOL` with `max(GAPTOL, 0.02 / n)` where n is
the current block size. This is adaptive: large n keeps the standard
1e-3 threshold; tiny n tightens it.

Verify: python3 evaluate.py should now report TOTAL: 379/379.
```

**Final state:** 379/379 tests passing on the standard `evaluate.py` suite.

---

## Optional: Bug #3 and Bug #4 — Py↔Fortran bit-identity (only needed for the pure-Fortran port)

Bugs #3 (DLARNV) and #4 (system BLAS DNRM2) in BUGREPORT.md are
**harmonization fixes** for matching a future pure-Fortran orchestration to
the Python orchestration bit-for-bit. They do not affect the standalone
`evaluate.py` suite. You add them only if you also do the Python→Fortran
port (round 3 work). Two further prompts:

```
Prompt 7 — DLARNV-seeded random vectors:
The Gram-Schmidt completion uses np.random.randn(n) for zero-σ
fill-in. To make Python and Fortran orchestrations produce
bit-identical completion vectors, replace np.random.randn with a
ctypes binding to LAPACK's DLARNV (IDIST=3, normal). Use a 4-tuple
seed derived deterministically from the column index, e.g.
(1, 3, 5, (2*j+1)&4095|1) for V-completion and
(2, 4, 6, (2*j+3)&4095|1) for U-completion.
```

```
Prompt 8 — system BLAS DNRM2:
np.linalg.norm calls numpy's bundled OpenBLAS, which differs from
system reference BLAS at the ~1 ULP level. Bind dnrm2_ via ctypes
from the same libblas the Fortran links against, and route every
post-processing normalization through it. Verify with
python3 test_fortran_match.py that all 379 specs match the Fortran
driver to within 2*eps.
```

---

## Checkpoint A (full reproduction from upstream Willems) — sketch

If you want to reproduce from scratch (no seed code), here are the 9 prompts
that, applied in order, evolve from upstream XMR + test scaffolding to
379/379:

1. **Read the paper.** "Willems & Lang 2012 ETNA Algorithm 4.1 describes an
   O(n²) bidiagonal SVD via T_GK + MR³. Read the paper, summarize the steps,
   identify what the existing xmr_src/ provides and what's missing."
2. **Build libxmr.so + ctypes glue.** "Compile every .f in xmr_src/ with
   `-fPIC -O2 -std=legacy`. Write `xmr_wrapper.c` that calls `dlaxre_` then
   `dlaxrv_`. Write `xmr_ctypes.py` that loads `libxmr.so` and exposes
   `xmr_eigenvectors(d, e, …)`."
3. **First-pass `mr3_gk.py:bidiag_svd`.** "T_GK perfect-shuffle, LAPACK
   dstebz for eigenvalues, ctypes call to XMR for eigenvectors. Return
   (sigma, U, V). Goal: ~22/64 small-test pass."
4. **Sign matrices D1, D2.** "T_GK uses |d|,|e|. Build D1, D2 such that
   B = D1·|B|·D2, apply them when extracting U, V from the T_GK eigenvectors."
5. **Activate dlaxre GK branch.** Bug #1 from BUGREPORT.md (replace
   `IF(.FALSE.)` with `IF(ABS(DMAX)+ABS(DMIN) < EPS·EMAX)`). Score ~142/379.
6. **Bv recovery + Gram-Schmidt completion.** Score ~290/379.
7. **Two-phase splitting (relative + absolute condition 3.5).** "Phase 1
   relative: split where `|e[i]| ≤ EPS·(|d[i]|+|d[i+1]|)`. Phase 2 absolute:
   within each Phase-1 block, split where `|e[i]| ≤ k·EPS·||B_block||`.
   Zero out e at split points." Score ~300/379.
8. **Apply Checkpoint B prompts 1–6 in order** to reach 379/379.

---

## Validating the recipe with `reproduce_with_agent.py`

The harness `python_fortran/experiments/reproduce_with_agent.py`
automates the whole thing:

```
cd python_fortran
python3 experiments/reproduce_with_agent.py --checkpoint B
```

What it does:

1. Creates a git worktree at the seed commit (`c6d73b2`) in `/tmp/repro_…/`.
2. For each prompt 1..6:
   - Runs `python3 evaluate.py` in the worktree to get the **before** score.
   - Spawns a Claude agent with `claude --print --max-turns 50` and the
     prompt as input, with the worktree's `python_fortran/` as the
     working directory.
   - After the agent finishes, runs `python3 evaluate.py` again to get
     the **after** score.
   - Records the result.
3. Writes `experiments/reproduce_results.json` with per-prompt
   (score_before, score_after, wall_time, files_changed).

Expected progression:

| Prompt | Pass before | Pass after (expected) | What also moves |
|---|---|---|---|
| 1 | 300/371 | ~310/379 | Singletons no longer time out; get_patterns regex fixed → suite grows to 379; SCORE gate lifts (5.00 → ~83) |
| 2 | ~310/379 | ~325/379 | gl_gradm/chkbd ortU drops |
| 3 | ~325/379 | ~325/379 | edge-case fixes (zero-d, neg-e) |
| 4 | ~325/379 | ~360/379 | saw_tooth/step_function/gl_wilkp pass via zero-shift QR deflation |
| 5 | ~360/379 | ~378/379 | Tight clusters pass (demmel_S1pe_*, three_clusters, pd_T0) via dlaxrb_clssfy fix |
| 6 | ~378/379 | **379/379** | two_clusters@10 passes via adaptive GAPTOL |

The exact numbers may shift by a few specs depending on agent execution
nuances; the validation harness records the actual progression.

Cost estimate: ~30–60 min wall, depending on agent latency. Each prompt
is independent so you can run them sequentially or interrupt and resume.

---

## Bug-report mapping (for cross-reference)

This table maps each prompt to the corresponding entry in `BUGREPORT.md`:

| Prompt | What it adds | BUGREPORT entry |
|---|---|---|
| Seed `c6d73b2` | dlaxre_gk activation | **Bug #1** |
| 1 | Singleton-block fast path | (Python orchestration) |
| 2 | Bv recovery + vectorized GS | (Python orchestration) |
| 3 | T_GK sub-splitting, neg-e | (Python orchestration) |
| 4 | Zero-shift QR + clean split | (Python preprocessing) |
| 5 | dlaxrb_clssfy AVGTHRESH guard | **Bug #2** |
| 6 | Adaptive GAPTOL | (Python orchestration tuning) |
| 7 (optional) | DLARNV seeded RNG | **Bug #3** |
| 8 (optional) | System BLAS DNRM2 | **Bug #4** |

The two XMR Fortran patches (Bugs #1 and #2) are the **only changes to
upstream Willems source**. Everything else is new Python orchestration
or new Fortran orchestration code that lives in our own files.
