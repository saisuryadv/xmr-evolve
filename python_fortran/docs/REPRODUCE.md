# How to reproduce the python_fortran MR3-GK SVD with an LLM agent

This document gives you the **prompt sequence** that, when fed to a Claude (or
similar) coding agent in order, will evolve a naive starting state — Paul
Willems' upstream XMR Fortran kernel + the test scaffolding only — into the
final working **379/379** MR3-GK bidiagonal SVD implementation.

## Historical evolution (22/64 → 379/379)

The full arc the agent will retrace:

| Step | Pass rate | Key change |
|---|---|---|
| 0 | 22/64 | Naive Python wrapper: T_GK perfect-shuffle + LAPACK `dstebz` for eigenvalues + ctypes call to XMR for eigenvectors. No sign handling, no recovery, no splitting. |
| 1 | 23/64 | Negative off-diag fix: T_GK off-diagonals can be negative when `d_i < 0`; XMR rejects them. Apply diagonal-similarity sign normalization before XMR. |
| 2 | 142/379 | Activate `dlaxre` GK branch (**Bug #1**): `IF(.FALSE.)` guard on the GK-detection branch in `dlaxre.f` was disabled. Replace with `IF(ABS(DMAX)+ABS(DMIN) < EPS·EMAX)`. Also add **D1/D2 sign matrices** for U/V recovery from |B|'s eigenvectors. |
| 3 | 290/379 | Bv recovery + Gram-Schmidt completion: when one half of a T_GK eigenvector is corrupt, recompute via `u = Bv/σ` (or `v = Bᵀu/σ`); for σ=0 columns, fill the orthogonal complement via Gram-Schmidt. |
| 4 | 300/379 | Two-phase splitting + singleton handling: Phase 1 relative split + Phase 2 absolute split (Willems-Lang condition 3.5). Singletons (k=1 blocks) directly produce ±1 columns. |
| 5 | ~310/379 | Bv recovery threshold + vectorized GS. |
| 6 | ~325/379 | T_GK sub-splitting for zero diagonals + neg-e edge case. |
| 7 | ~360/379 | Zero-shift QR deflation (Demmel-Kahan 1990) for σ_min = 0. |
| 8 | ~378/379 | `dlaxrb_clssfy` `AVGTHRESH=0` fix (**Bug #2**). |
| 9 | **379/379** | Adaptive GAPTOL `max(GAPTOL, 0.02/n)` for tiny blocks. |

## Initial state

You start from Paul Willems' upstream code + the test scaffolding only — no
Python implementation, no SVD entry point, no orchestration layer.

The starting state is the folder **`python_fortran/initial_code/`** in this
repo. It contains exactly:

```
initial_code/
├── README.md                 # explains contents + verification steps
├── xmr_src/                  # Paul Willems' Fortran XMR source — UNTOUCHED
│                             #   44 .f files; both bugs from BUGREPORT.md
│                             #   are present in their original buggy form:
│                             #   * dlaxre.f:268      — IF(.FALSE.) (Bug #1)
│                             #   * dlaxrb_clssfy.f:361,412,438
│                             #                       — unguarded MIN (Bug #2)
├── stcollection/             # 19 .dat files (Marques et al., untouched)
├── full_eval.py              # 90 adversarial pattern generators
└── evaluate.py               # 379-test scoring harness
```

To run the recipe, copy this folder to a fresh workspace:

```
cp -r python_fortran/initial_code /tmp/repro_workspace
cd /tmp/repro_workspace
# feed the prompts below in order to a fresh `claude --print` agent for each
```

XMR is **eigenvalue-only** — there is no SVD entry point, no T_GK
construction, no preprocessing or post-processing for bidiagonals. The agent
will build all of that. It also has to re-create the two Fortran patches
(`dlaxre_gk.f`, `dlaxrb_clssfy_fix.f`) on top of the unmodified upstream
kernel.

The agent has to produce:

- `mr3_gk.py` — Python orchestration entry point (`bidiag_svd(d, e)`)
- `xmr_wrapper.c` — C glue calling the Fortran kernel
- `xmr_ctypes.py` — Python ctypes binding
- `build.sh` — compiles `xmr_src/` + wrappers + patches into `libxmr.so`
- `dlaxre_gk.f` — patched copy of `xmr_src/dlaxre.f`, Bug #1 fix
- `dlaxrb_clssfy_fix.f` — patched copy of `xmr_src/dlaxrb_clssfy.f`, Bug #2 fix

XMR has two latent bugs that block bidiagonal SVD entirely (also in
`docs/BUGREPORT.md`):

1. **`dlaxre.f` line 268**: GK-detection branch is `IF(.FALSE.)` — disabled.
2. **`dlaxrb_clssfy.f` lines 361, 412, 445**: classifier uses
   `MIN(AVGTHRESH, ABSMAX·GAPTHRESH)`, but at deep tree nodes
   `AVGTHRESH = 2·AVGTOL · 2^(−DEPTH)` underflows to 0, making the test
   trivially true and over-fragmenting clusters.

The prompt sequence below fixes both.

---

## Prompt sequence (10 prompts)

Each prompt is a single message to the agent. Use a fresh agent session per
prompt so each is self-contained — the agent has no memory of prior prompts,
only the current state of the working tree. The agent should have read/write
tooling and shell access (run `python3 evaluate.py`, `bash build.sh`, and
Fortran compilation).

### Prompt A0 — read the paper and survey

```
Inside the workspace you have:
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
    Copy xmr_src/dlaxre.f to dlaxre_gk.f.
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
Two algorithmic refinements:

(a) Two-phase bidiagonal splitting per Willems-Lang condition (3.5).
    Phase 1 (relative): split where |e[i]| <= eps * (|d[i]| + |d[i+1]|).
    Phase 2 (absolute): within each Phase-1 block of size k_sub,
      compute ||B_sub||_inf, split where |e[i]| <= k_sub * eps * ||B_sub||.
    Zero out e[i] at every split point so cross-block contamination
    can't happen during Bv recovery.

(b) Singleton-block fast path: when splitting produces a 1x1 block,
    don't call XMR. Direct write:
      sigma = |d[i]|
      U[i, col] = d1[i]   # already +/-1
      V[i, col] = d2[i]   # already +/-1

Re-run evaluate.py. Should now report ~300/379.
```

### Prompt A6 — Bv recovery threshold + GS vectorization (~310/379)

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

Verify on gl_gradm and chkbd — ortU should drop to ~1 nε.
```

### Prompt A7 — T_GK sub-splitting + negative-e fix (~325/379)

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

### Prompt A8 — zero-σ deflation via zero-shift QR (~360/379)

```
The Demmel-Kahan 1990 paper (Sec. 3) shows that one zero-shift QR
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
     recursively (on the leading (k-1)x(k-1) part), embed the
     solution into kxk with sigma=0 in the bottom-right, and apply
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

### Prompt A9 — dlaxrb_clssfy AVGTHRESH=0 fix (~378/379)

```
A subtle bug in upstream xmr_src/dlaxrb_clssfy.f. At deep tree
nodes, AVGTHRESH = 2*AVGTOL collapses to 0 because AVGTOL is scaled
by 2^(-DEPTH) and underflows. The classifier test is

   IF( SPREAD .GE. MIN(AVGTHRESH, ABSMAX*GAPTHRESH) ) gap is full

— with AVGTHRESH=0 this becomes SPREAD >= 0, trivially true, so
EVERY adjacent eigenvalue pair is labeled a full gap and the rep
tree degenerates into singletons. Tight-cluster matrices fail
catastrophically.

Steps:

1. Create dlaxrb_clssfy_fix.f as a copy of
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

### Prompt A10 — adaptive GAPTOL for tiny blocks (379/379)

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

## Optional: Bug #3 and Bug #4 — Py↔Fortran bit-identity (only for the pure-Fortran port)

Bugs #3 (DLARNV) and #4 (system BLAS DNRM2) in `BUGREPORT.md` are
**harmonization fixes** for matching a future pure-Fortran orchestration to
the Python orchestration bit-for-bit. They do not affect the standalone
`evaluate.py` suite. You add them only if you also do the Python→Fortran
port. Two further prompts:

### Prompt A11 — DLARNV-seeded random vectors (Bug #3)

```
The Gram-Schmidt completion uses np.random.randn(n) for zero-σ
fill-in. To make Python and Fortran orchestrations produce
bit-identical completion vectors, replace np.random.randn with a
ctypes binding to LAPACK's DLARNV (IDIST=3, normal). Use a 4-tuple
seed derived deterministically from the column index, e.g.
(1, 3, 5, (2*j+1)&4095|1) for V-completion and
(2, 4, 6, (2*j+3)&4095|1) for U-completion.
```

### Prompt A12 — system BLAS DNRM2 (Bug #4)

```
np.linalg.norm calls numpy's bundled OpenBLAS, which differs from
system reference BLAS at the ~1 ULP level. Bind dnrm2_ via ctypes
from the same libblas the Fortran links against, and route every
post-processing normalization through it. Verify with
python3 test_fortran_match.py that all 379 specs match the Fortran
driver to within 2*eps.
```

---

## Bug-report mapping (for cross-reference)

This table maps each prompt to the corresponding entry in `BUGREPORT.md`:

| Prompt | What it adds | BUGREPORT entry |
|---|---|---|
| A0 | Read paper + survey | – |
| A1 | Naive Python wrapper | (Python orchestration scaffold) |
| A2 | Negative off-diag sign normalization | (Python orchestration) |
| A3 | Activate dlaxre GK branch + D1/D2 sign matrices | **Bug #1** |
| A4 | Bv recovery + Gram-Schmidt completion | (Python orchestration) |
| A5 | Two-phase splitting + singleton fast-path | (Python orchestration) |
| A6 | Bv recovery threshold + vectorized GS | (Python orchestration) |
| A7 | T_GK sub-splitting + negative-e fix | (Python orchestration) |
| A8 | Zero-shift QR + clean-split deflation | (Python preprocessing) |
| A9 | dlaxrb_clssfy AVGTHRESH guard | **Bug #2** |
| A10 | Adaptive GAPTOL | (Python orchestration tuning) |
| A11 (optional) | DLARNV seeded RNG | **Bug #3** |
| A12 (optional) | System BLAS DNRM2 | **Bug #4** |

The two XMR Fortran patches (Bugs #1 and #2) are the **only changes to
upstream Willems source**. Everything else is new Python orchestration
or new Fortran orchestration code that lives in our own files.
