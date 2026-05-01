# How to reproduce the python_fortran MR3-GK SVD with an LLM agent

This document gives you the **prompt sequence** that, when fed to a Claude (or
similar) coding agent in order, will evolve a naive starting state into the
final working 379/379 implementation. It also documents what initial code we
built on top of (Paul Willems' XMR Fortran library) and the historical
evolution that produced our seed commit.

## Historical evolution (22/64 → 300/379)

The `python_fortran/` directory's first commit `c6d73b2` already passed
300/379 — but it represents the *end-state* of about two weeks of evolution
on top of upstream Willems XMR. Each step is documented in the seed commit's
`CHANGES.md`. The full arc:

| Step | Pass rate | Key change |
|---|---|---|
| 0 | 22/64 | Naive Python wrapper: T_GK perfect-shuffle + LAPACK `dstebz` for eigenvalues + ctypes call to XMR for eigenvectors. No sign handling, no recovery, no splitting. |
| 1 | 23/64 | **Negative off-diag fix**: T_GK off-diagonals can be negative when bidiagonal d_i < 0; XMR rejects them. Apply diagonal-similarity sign normalization before XMR. |
| 2 | 142/379 | **Activate dlaxre GK branch** (Bug #1): `IF(.FALSE.)` guard on the GK-detection branch in `dlaxre.f` was disabled. Replace with `IF(ABS(DMAX)+ABS(DMIN) < EPS·EMAX)`. Also add the **D1/D2 sign matrices** for U/V recovery from |B|'s eigenvectors. |
| 3 | 290/379 | **Bv recovery + Gram-Schmidt completion**: when one half of a T_GK eigenvector is corrupt, recompute via `u = Bv/σ` (or `v = Bᵀu/σ`); for σ=0 columns, fill in the orthogonal complement via Gram-Schmidt. |
| 4 | **300/379** | **Two-phase splitting + singleton handling**: phase-1 relative split + phase-2 absolute split (Willems-Lang condition 3.5). Singletons (k=1 blocks) directly produce ±1 columns. ← **Seed commit `c6d73b2`** |

The 6 prompts in this recipe take the 300/379 seed to 379/379. The 10 prompts
in the optional "Checkpoint A" appendix walk through the entire 22/64 →
379/379 arc.

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

### Checkpoint B — seed at commit `c6d73b2` + modern evaluator (300/379)

You start from the python_fortran/ first commit (which is the day the
`python_fortran/` directory was first created), but with **the latest
`evaluate.py` and `full_eval.py` overlaid from project HEAD**. The
algorithm code being evolved (`mr3_gk.py`, `dlaxre_gk.f`, `libxmr.so`,
`xmr_wrapper.c`, etc.) is exactly the seed-commit version; only the
test scoreboard is upgraded.

Why the overlay: the seed-commit `evaluate.py` has a regex parsing
bug in `get_patterns()` that misses 2 of 90 adversarial patterns
(`chkbd_4`, `chkbd_16`, generated from a parameterized loop), so its
suite reports 371 tests instead of 379. The modern evaluator (which
fixed that regex with `importlib`) reports the full 379 against the
same algorithm code. Using the modern scoreboard from the start means
the recipe is measured against the standard 379-test suite all the
way through, which makes the progression numbers cleaner and matches
the headline 379/379 you're trying to reproduce.

After the overlay, the worktree contains:

- **From `c6d73b2` (the algorithm and substrate):**
  - `libxmr.so` built from upstream `xmr_src/` (with our XMR Bug #1 fix
    patched in — `dlaxre.f` → `dlaxre_gk.f`).
  - `xmr_wrapper.c` C glue + `xmr_ctypes.py` ctypes binding.
  - `mr3_gk.py` Python orchestration (~640 lines): T_GK construction,
    splitting, sign matrices D1/D2, Bv recovery, basic GS completion.
  - 19 STCollection `.dat` files.

- **From HEAD (the test scoreboard only):**
  - `evaluate.py` (modern; reports the 379-test full suite).
  - `full_eval.py` (90 adversarial pattern generators).

`python3 evaluate.py` on this configuration reports **300/379**
passing. `evaluate.py` also reports a **SCORE** (pass-rate × O(n²)
timing-ratio gate). At this state it's hard-gated at `SCORE = 5.00`
by single-run timing noise; the modern evaluator's median-of-5 timing
helps stabilize it once the prompt-1 singleton bypass also drops the
timing outliers from the slowest tests.

The remaining 79 failures fall into three groups:

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

### Prompt 1 — singleton-block bypass (300/379 → ~310/379)

```
Run `python3 evaluate.py` and observe the current score (you will see
TOTAL: 300/379 — the chkbd_*, saw_tooth, and step_function patterns
at n >= 200 are timing out).

The root cause is in mr3_gk.py:mr3_tgk: when the bidiagonal splitter
produces a singleton block (k=1), the code still calls XMR on a
trivial 2x2 T_GK matrix. XMR returns degenerate eigenvectors with
||v||=1, ||u||=0 (only even rows populated). The post-processing
then runs O(n^3) Gram-Schmidt on hundreds of such columns.

Add a fast path in mr3_tgk that handles singleton blocks directly,
bypassing the XMR call:
   - sigma = |d[i]|
   - both u and v components are 1/sqrt(2) at the corresponding row
   - apply sign matrices D1, D2 (constructed earlier in bidiag_svd)
     to recover correct signs

Do NOT modify dlaxre_gk.f or any Fortran file. Re-run evaluate.py.
The chkbd patterns must no longer time out. Score should jump to
roughly 305-320/379 and SCORE should rise above 5.00 once the
timing outliers are gone.
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

## Checkpoint A (full reproduction from upstream Willems — 10 prompts)

If you want to reproduce the *entire* 22/64 → 379/379 arc starting from
upstream Willems XMR + the test scaffolding only (no seed Python code),
here are the 10 prompts in order. Each is a single message to a fresh
agent session.

The starting state for Checkpoint A is:

```
python_fortran_a/
├── xmr_src/                  # upstream Willems Fortran (44 .f files)
├── stcollection/             # 19 STCollection .dat files
├── full_eval.py              # 90 adversarial pattern generators
└── evaluate.py               # 379-test scoring harness
```

(No `mr3_gk.py`, no `dlaxre_gk.f`, no `libxmr.so`, no C wrapper, no
ctypes glue — the agent writes all of those.)

### Prompt A0 — read the paper and survey

```
Inside python_fortran_a/ you have:
  - xmr_src/ : Paul Willems' XMR Fortran source (44 files, ~14K LOC)
  - stcollection/ : 19 .dat files of bidiagonal test matrices
  - full_eval.py : 90 adversarial pattern generators (function `make`,
    list `adv_names`)
  - evaluate.py : 379-test scoring harness

Read the Willems-Lang 2012 ETNA preprint (Algorithm 4.1, T_GK
construction, the MR^3 representation tree). Then survey xmr_src/
and identify:
  - the public entry points (DSTEXR, the dlaxr* family),
  - what XMR provides (tridiagonal eigenpairs only — no SVD),
  - what we'll need to build on top of XMR for bidiagonal SVD.

Output: a short SURVEY.md summarizing the paper's algorithm and what
the user-side wrapper has to provide.
```

### Prompt A1 — naive Python wrapper (target: ~22/64 on small subset)

```
Build the basic Python+Fortran pipeline:

(1) Compile xmr_src/*.f into libxmr.so:
    gfortran -fPIC -O2 -std=legacy -w -fno-second-underscore -c
    For each .f, then link with `gcc -shared -o libxmr.so *.o
    -lgfortran -llapack -lblas -lm`. Put a build.sh script in place.

(2) Write xmr_wrapper.c — a minimal C function `xmr_eigenvectors`
    that calls `dlaxre_` (root representation) followed by
    `dlaxrv_` (eigenvectors) and returns (eigenvalues, eigenvectors)
    by pointer.

(3) Write xmr_ctypes.py — a thin ctypes wrapper exposing
    `xmr_eigenvectors(d, e, n_eigvecs)` from libxmr.so to Python.

(4) Write mr3_gk.py:bidiag_svd(d, e):
    - Form T_GK = perfect_shuffle([0]*2n with off-diag [|d_0|, |e_0|,
      |d_1|, |e_1|, ...]).
    - Call scipy.linalg.lapack.dstebz to get the n positive
      eigenvalues of T_GK (these are the singular values).
    - Call xmr_eigenvectors(0..0, off_diag, n) to get eigenvectors
      of T_GK.
    - De-interleave: V[i,j] = z[2i, j], U[i,j] = z[2i+1, j].
    - Return (sigma, U, V, info=0).

No sign handling, no Bv recovery, no GS completion, no splitting.

Run `python3 evaluate.py --small` (or the smallest subset). Expect
about 22/64 passing — the negative-bidiag matrices and tight clusters
will fail.
```

### Prompt A2 — negative off-diagonal sign normalization (target: 23/64)

```
Several test bidiagonals have d[i] < 0 or e[i] < 0. The T_GK
perfect-shuffle propagates those signs into off-diagonals, which
XMR rejects. Normalize before calling XMR:

  build the T_GK off-diagonal vector with absolute values: |d_i|, |e_i|
  track sign flips
  remember to undo them when extracting U, V.

Score should creep to 23/64.
```

### Prompt A3 — activate dlaxre GK branch + sign matrices D1, D2 (target: 142/379)

```
Two changes for the big jump:

(a) Activate the GK branch in dlaxre.f.
    Copy xmr_src/dlaxre.f to python_fortran_a/dlaxre_gk.f.
    Around line 269 you'll find:
        IF( .FALSE. )THEN
        C    Support for GK-type matrices deactivated.
    Replace with:
        IF( ABS(DMAX)+ABS(DMIN) .LT. EPS*EMAX )THEN
        C    GK-type matrix: constant (zero) diagonal.
    This activates the GK-aware root representation per Willems-Lang
    Algorithm 4.1.

    Modify build.sh to compile dlaxre_gk.f and link its .o in place
    of xmr_src/dlaxre.o. Rebuild libxmr.so.

(b) Sign matrices D1, D2 in mr3_gk.py.
    T_GK uses |d|, |e|; its eigenvectors are for |B|, not B itself.
    To recover U, V for the original B, build:
        d1[0] = 1.0
        d2[0] = sign(d[0]) if d[0] != 0 else 1.0
        for i = 0..n-2:
            s_e = sign(e[i]) if e[i] != 0 else 1.0
            d2[i+1] = s_e / d1[i] if abs(d1[i]) > 0 else s_e
            s_d = sign(d[i+1]) if d[i+1] != 0 else 1.0
            d1[i+1] = s_d / d2[i+1] if abs(d2[i+1]) > 0 else s_d

    Apply when extracting:
        V[:,j] = z[0::2,j] * d2
        U[:,j] = z[1::2,j] * d1

Run `python3 evaluate.py`. Score should jump to ~142/379. This is
the largest single improvement in the entire arc (GK structure
preservation is the core of the algorithm).
```

### Prompt A4 — Bv recovery + Gram-Schmidt completion (target: ~290/379)

```
After step A3, many singular triplets pass but some have one corrupt
half (e.g., ||u_j|| ≈ 0 while ||v_j|| ≈ 1, or vice versa).

(a) Bv recovery. After extracting U, V, sigma:
    - For each column j with sigma[j] > sigma_max * eps:
      compute B*v[:,j] (one bidiagonal matvec) and B^T*u[:,j].
      If ||B*v - sigma*u|| > threshold, recompute u = B*v/sigma.
      If ||B^T*u - sigma*v|| > threshold, recompute v = B^T*u/sigma.
      Both branches: re-normalize.

(b) Gram-Schmidt completion. For columns j with sigma[j] near zero:
    one half (u_j or v_j) may be all-zero; fill in via Gram-Schmidt
    against all the other already-good columns to maintain
    orthogonality. Two passes of MGS give numerical stability.

Run evaluate.py. Score should reach ~290/379. The remaining failures
are clustered/graded matrices and tight zero-σ cases.
```

### Prompt A5 — two-phase splitting + singleton fast-path (target: 300/379)

```
Two algorithmic refinements that bring us to the seed-commit state:

(a) Two-phase bidiagonal splitting per Willems-Lang condition (3.5).
    Phase 1 (relative): split where |e[i]| <= eps * (|d[i]| + |d[i+1]|).
    Phase 2 (absolute): within each phase-1 block of size k_sub,
      compute ||B_sub||_inf, split where |e[i]| <= k_sub * eps * ||B_sub||.
    Zero out e[i] at every split point so cross-block contamination
    can't happen during Bv recovery.

(b) Singleton-block fast path: when splitting produces a 1x1 block,
    don't call XMR. Direct write:
      sigma = |d[i]|
      U[i, col] = d1[i] (already +/-1)
      V[i, col] = d2[i] (already +/-1)

Re-run evaluate.py. Should now report 300/379. This matches the
state of seed commit c6d73b2 in our git history.
```

### Prompts A6 – A10 — same as Checkpoint B prompts 1 – 6

From here, follow Checkpoint B's prompts 1 through 6 verbatim.
They take 300/379 → 379/379.

| Prompt | Title | Target |
|---|---|---|
| A6 = B1 | (skip — already done in A5; or apply the BV-threshold tightening if you want a cleaner singleton fast-path) | – |
| A7 = B2 | Bv recovery threshold + GS vectorization | ~325/379 |
| A8 = B3 | T_GK sub-splitting + negative-e fix | ~325/379 |
| A9 = B4 | Zero-shift QR deflation | ~360/379 |
| A10 = B5 | dlaxrb_clssfy AVGTHRESH=0 fix | ~378/379 |
| A11 = B6 | Adaptive GAPTOL | **379/379** |

Total: ~10 distinct prompts to evolve from "naive Python+Fortran
wrapper" to 379/379. The prompts A0-A5 are not validated by the
agent harness in this repo — they're documented for reproducibility
clarity but constructing the 22/64 starting state is itself a
research project. The harness validates A6-A11 (= Checkpoint B)
which is the empirically tested portion.

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

Expected progression (against the 379-test scoreboard from the start):

| Prompt | Pass before | Pass after (expected) | What also moves |
|---|---|---|---|
| 1 | 300/379 | ~310/379 | chkbd_*/saw_tooth no longer time out; SCORE gate lifts (5.00 → ~83) |
| 2 | ~310/379 | ~325/379 | gl_gradm/chkbd ortU drops |
| 3 | ~325/379 | ~325/379 | edge-case fixes (zero-d, neg-e) |
| 4 | ~325/379 | ~360/379 | saw_tooth/step_function/gl_wilkp pass via zero-shift QR deflation |
| 5 | ~360/379 | ~378/379 | Tight clusters pass (demmel_S1pe_*, three_clusters, pd_T0) via dlaxrb_clssfy fix |
| 6 | ~378/379 | **379/379** | two_clusters@10 passes via adaptive GAPTOL |

The exact numbers may shift by a few specs depending on agent
execution nuances; the validation harness records the actual
progression.

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
