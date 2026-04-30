# 05 — T_GK Path: Upstream XMR vs Our Code

This document traces the bidiagonal-SVD code path through both upstream XMR (`xmr_src/`) and our wrapper (`mr3gk_fortran/` / `mr3_gk.py`). Deliverable for advisor task **(f)**.

---

## TL;DR

| | Upstream XMR (`xmr_src/`) | Our code (`mr3gk_fortran/` + `mr3_gk.py`) |
|---|---|---|
| Public bidiagonal-SVD entry point | **none** | `dmr3gk_svd` (Fortran) / `bidiag_svd` (Python) |
| Public eigenpair entry point | `dstexr` (selected eigenpairs of symmetric tridiagonal) | (re-uses XMR via `dlaxre_gk_` + `dlaxrv_`) |
| T_GK construction | not implemented | `mr3gk_tgk.f90:139-144` (Fortran), `mr3_gk.py:580-584` (Python) |
| GK-aware root representation | gated off (`IF( .FALSE. )`) | activated (`dlaxre_gk.f:269`) |
| Two-phase splitting | not implemented | `mr3gk_split.f90` |
| Zero-σ deflation pre-pass | not implemented | `mr3gk_qrsweep.f90` (Demmel-Kahan zero-shift QR) |
| Sign matrices D1, D2 | not implemented | inline in `mr3gk.f90` |
| Post-processing (de-interleave, recover, sign-fix, GS-completion) | not implemented | `mr3gk_postproc.f90` |
| Eigenvector solve | `dlaxrv` | (same `dlaxrv` — we link the same `libxmr.so`) |

The upstream XMR codebase ships **no driver** for bidiagonal SVD. It exposes only the symmetric tridiagonal eigenpair solver `dstexr` (and helper layers `dlaxra/dlaxre/dlaxrf/dlaxrv`). Everything between "I have a bidiagonal `(d, e)`" and "I have `(σ, U, V)`" is in our code.

---

## What is on the T_GK path in XMR (i.e. nothing about SVD specifically)

`xmr_src/` contains 44 Fortran files. The eigenpair stack from outermost in:

```
dstexr                              public entry (selected eigenpairs of T)
  ├─ dlaxra                         block split / sign normalize / scale
  ├─ dlaxri                         locate block containing WIL:WIU
  └─ for each block:
       └─ dlaxre                    build root representation (RRR) from (D,E)
                                    >>> upstream branch IF( .FALSE. ) at line 269
                                        means the GK-detection path is dead.
            ├─ dlaxre_initewldqds   initial eigenvalue list from DQDS bisection
            └─ dlaxrl_update        update gap info on the EWL
       └─ dlaxrv                    compute eigenpairs WIL:WIU
            ├─ dlaxrf               representation tree forest
            │    ├─ dlaxrb_clssfy   classify clusters / detect full gaps  >>> our patch lives here
            │    ├─ dlaxrf_selshf   select shift candidate
            │    ├─ dlaxrf_seltw    select twist position
            │    ├─ dlaxrf_iib      inner-interval bisection
            │    ├─ dlaxrf_cob      change-of-basis verification
            │    └─ dlaxrr / dlaxrs / dlaxrt   build / shift representations
            ├─ dlaxrn / dlaxrn0     representation norms
            └─ dlaxrm               twisted-LDL-T eigenvector assembly
```

There is **no SVD entry point**. To reach a bidiagonal SVD using only upstream XMR, one would have to:
1. Build T_GK = perfect-shuffle of `[0,…,0; β_1,β_2,…]` (a 2n × 2n symmetric tridiagonal) by hand.
2. Call `dstexr` with `WIL = n+1, WIU = 2n` to get the n positive eigenvalues = singular values, plus their eigenvectors of length 2n.
3. De-interleave each eigenvector `z` into `[u; v]` halves.
4. Renormalize each half to unit length, with bookkeeping for the overall sign.
5. Detect zero singular values (which produce degenerate ±0 eigenvalue pairs in T_GK), and Gram-Schmidt-complete the missing left/right vectors.
6. Re-undo the sign matrices D1, D2 (which our code applies before forming T_GK to make the off-diagonals nonneg).
7. Undo any pre-scaling.

None of these steps are in upstream XMR.

---

## What is on the T_GK path in our code

### Fortran orchestration (`mr3gk_fortran/`)

```
dmr3gk_svd                              (mr3gk.f90)            public entry
  ├─ pre-scaling                        (mr3gk.f90:68-107)     factor = max(|d|,|e|); divide
  ├─ dmr3gk_split                       (mr3gk_split.f90)      two-phase split into blocks
  └─ for each block [bbeg..bend]:
       ├─ build sign matrices D1, D2   (mr3gk.f90:119-158)   make all e ≥ 0
       ├─ dmr3gk_qrsweep                (mr3gk_qrsweep.f90)    zero-shift QR pass
       │     └─ DLARTG (LAPACK)         per-Givens rotation
       │     └─ apply_givens_rows       reverse-order accumulation into U
       ├─ if clean split detected:      (mr3gk.f90:201-215)
       │     └─ recursive dmr3gk_svd    on the remaining sub-block
       └─ otherwise:
            └─ mr3_tgk_multiblock       (mr3gk_tgk.f90:123)
                 ├─ build T_GK         (mr3gk_tgk.f90:135-144)   2k×2k tridiag, [0,...,0; |d|,|e| interleaved]
                 └─ solve_tgk_block    (mr3gk_tgk.f90:22-118)
                      ├─ bisect_evals  (mr3gk_utils.f90)   call DSTEBZ for the n positive evals
                      ├─ DLAXRE_GK     (dlaxre_gk.f, OUR PATCH)   build RRR; activates GK-detection branch
                      ├─ wsreq_xrv     (xmr_src)                  workspace query
                      └─ DLAXRV        (xmr_src)                  compute eigenpairs
                           └─ DLAXRF → DLAXRB_CLSSFY (OUR PATCH)
       └─ post-processing               (mr3gk_postproc.f90)
            ├─ de-interleave even/odd halves of each z
            ├─ DNRM2 column-by-column normalization
            ├─ Bv recovery for one-sided-good vectors
            ├─ residual check, sign fix
            └─ GS completion via DLARNV (zero-σ slots)
  └─ rescale σ                          (mr3gk.f90:307-311)
```

### Python orchestration (`mr3_gk.py:bidiag_svd`)

Identical structure (it's the algorithm we ported from). Calls into the same `libxmr.so` via ctypes. Lines:
- pre-scaling: 723–746
- splitting: `split_bidiag` 225–257
- sign matrices: 764–770
- QR sweep + clean-split: 802–809
- T_GK multiblock: `mr3_tgk_multiblock` 523–655 → ctypes call to `xmr_eigenvectors` → `dlaxre_gk_` → `dlaxrv_`
- post-processing: 857–987

Both orchestrators produce identical bytes (verified to 2·eps on 21,375 specs).

---

## Where our two XMR patches sit on this path

| Patch | Path location | What it does |
|---|---|---|
| `dlaxre_gk.f:269` (activate GK branch) | inside `dlaxre`, called once per block at the top of `solve_tgk_block` | When all diagonal entries are zero (T_GK), use the matrix directly as a root representation instead of falling through to the post-GK Gershgorin path. Without this, the eigenvector quality degrades catastrophically for any bidiagonal-SVD problem. |
| `dlaxrb_clssfy_fix.f:358-453` (AVGTHRESH guard) | inside `dlaxrb_clssfy`, called recursively from `dlaxrf` (representation-tree traversal) | When the depth-scaled `AVGTHRESH` underflows to 0, fall back to the relative threshold only. Without this, deep tree nodes mis-classify every adjacent eigenvalue pair as a "full gap," fragmenting clusters. |

Neither patch changes the kernel arithmetic; they only fix dead-code branches and a degenerate `MIN(0, x)` test.

---

## Diagrammatic summary

```
   bidiagonal (d, e)                  symmetric tridiagonal (Dt, Et)
     │                                       │
     │   <-- our code lives here -->         │
     │                                       │
     ▼                                       ▼
   ┌──────────────────────┐            ┌────────────────┐
   │ dmr3gk_svd / bidiag  │            │ dstexr (XMR    │
   │ (mr3gk.f90,          │            │  public entry) │
   │  mr3_gk.py)          │            └─────┬──────────┘
   └─────────┬────────────┘                  │
             │  scale, split, signs,         │
             │  zero-shift QR,               │
             │  build T_GK ────────────────► │
             │                               │
             │                               ▼
             │                       ┌──────────────┐
             │                       │ dlaxre_gk    │  ← OUR PATCH (activate GK branch)
             │                       └──────┬───────┘
             │                              │
             │                              ▼
             │                       ┌──────────────┐
             │                       │ dlaxrv (XMR) │
             │                       │  └─ dlaxrf   │
             │                       │     └─ dlaxrb_clssfy_fix │  ← OUR PATCH (AVGTHRESH guard)
             │                       └──────┬───────┘
             │                              │
             ▼                              ▼
   ┌──────────────────────┐         eigenpairs of T_GK
   │ post-processing      │ ◄────── (z of length 2n)
   │ (mr3gk_postproc.f90, │
   │  mr3_gk.py)          │
   │ - de-interleave      │
   │ - normalize          │
   │ - Bv recover         │
   │ - sign fix           │
   │ - GS completion      │
   └──────────┬───────────┘
              │
              ▼
        (σ, U, V)
```

The dashed boundary "our code lives here" is the entire bidiagonal-SVD layer. It is ~1,400 lines of new Fortran (`mr3gk_fortran/`) plus the same logic in Python (`mr3_gk.py`). The two XMR patches are tiny (~32 lines combined) and live inside the kernel layer, sitting on the path that our orchestrator drives.
