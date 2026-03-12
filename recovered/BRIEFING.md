# Research Briefing: Prior Agent Work on XMR/DSTEXR Algorithm 4.1

**READ THIS ENTIRE FILE BEFORE DOING ANYTHING ELSE.**

This document summarizes what a previous agent accomplished, what worked, what didn't,
and exactly what remains to be done. Do not repeat work that has already been done.

## Previous Agent Session

- **JSONL transcript**: `/Users/saisurya/.claude/projects/-Users-saisurya-MRRR-xmr-evolve/b05f34bb-76b2-42f5-8184-65ab5d3937ea.jsonl`
- **Recovered variants**: `/Users/saisurya/MRRR/xmr-evolve/recovered/SUMMARY.md` — full details of all 19 versions
- **Recovered code**: `/Users/saisurya/MRRR/xmr-evolve/recovered/version_NNN_passXXX.h` — 19 snapshots

## What the Previous Agent Did

### 1. Modified Fortran: Re-enabled GK-form in dlaxre.f (Bug #10 fix)

The most important change. In `dlaxre.f` (root representation builder), GK-form
support was deactivated with `IF(.FALSE.)THEN`. The agent changed this to:
```fortran
IF( (DMAX-DMIN) .LE. 2*ABSERR )THEN
```
This re-enables Algorithm 4.1's GK-form at the root level — when the TGK matrix
has a nearly constant diagonal (which it always does since TGK has zero diagonal),
the representation is flagged as GK-form so the eigensolver can preserve GK structure.

**Impact**: Improved pass rate from 61/379 (baseline STEXR) to 114-138/379.

The modified Fortran is at: `/Users/saisurya/MRRR/xmr-evolve/recovered/other_files/xmr_src_dlaxre.f`
The corresponding C translation is still on disk at: `/Users/saisurya/MRRR/xmr-evolve/src/xmr_c/dlaxre.c`

### 2. Added C runtime toggle for GK-form

Added `int xmr_disable_gkform_ = 0;` global flag in `dlaxre.c` so the C++ code
can enable/disable GK-form at runtime via `extern int xmr_disable_gkform_;`.

### 3. Iterated on C++ post-processing (19 versions)

The agent tried many post-processing strategies in `bidiag_svd.h`:
- Normalization
- Bidirectional one-sided recovery (U=Bv/σ and V=B^Tu/σ with residual comparison)
- Sign consistency (flip u if dot(Bv, u) < 0)
- Chunked MGS reorthogonalization (MAX_CHUNK=32, O(n²))
- Basis completion for zero singular values
- Quick orthogonality check for quality assessment

## Key Results — AND WHY PASS RATE ALONE IS MISLEADING

**IMPORTANT: The SCORE is what matters, not the pass rate.** You can have a high pass rate
and still get a terrible score. The previous agent's best pass rate was 207/289 but its
score was only **5.0** (the minimum) because scaling was 19.5x, triggering the hard gate.
Meanwhile v013 with only 138/289 pass rate scored **64.3** — 13x better — because it had
good O(n²) scaling. A solution with 100/289 pass rate and 4.0x scaling would outscore
207/289 with 6.0x scaling.

| Version | Pass Rate | Score | Scaling | Key Feature |
|---------|-----------|-------|---------|-------------|
| v013 (best balanced) | 138/289 | **64.3** | **4.29x** | GK-form + chunked reortho, O(n²) compliant |
| v017 (best accuracy) | 207/289 | **5.0** | 19.50x | Added basis completion — O(n³), SCORE KILLED |

**The previous agent wasted significant effort on post-processing and O(n³) techniques
that boosted pass rate but destroyed the score.** Do not repeat this mistake.

### Why 138→207 required O(n³)

The 69 extra tests that v017 passes over v013 require **basis completion** — when
DSTEXR returns eigenvectors with zero V or U components (for zero/near-zero singular
values), v017 fills them by orthogonalizing canonical basis vectors against ALL existing
columns. This is O(n²) per bad column × n columns worst case = O(n³).

v013 avoids this by limiting one-sided recovery to n≤80 and using only chunked MGS.

### The fundamental gap

138 tests pass with pure GK-form + O(n²) post-processing.
207 tests pass when you add O(n³) basis completion.
The 69-test gap represents cases where **DSTEXR's eigenvectors fundamentally lack
GK structure** — the extracted U/V are garbage and no amount of O(n²) post-processing
can fix them. Only O(n³) reorthogonalization can rescue these.

## What Has NOT Been Done (The Actual Work Remaining)

The previous agent only modified `dlaxre.f` at the **root level**. But Algorithm 4.1
requires NCD enforcement at **every level** of the representation tree. Here is exactly
what's missing:

### IMPORTANT: NCD was PLANNED but NEVER IMPLEMENTED in original XMR

A search of the entire XMR codebase reveals:
- `dlaxrf_seltw.f` has NCD documented as an output parameter: "returned in ELG, RCOND, NCD"
- RWORK(5) and RWORK(6) are allocated as "reserved for ncd"
- **BUT the actual NCD computation code was never written** — these values always remain at -1
- GK-form block infrastructure EXISTS in `dlaxrm_stat*.f` files (LBBEGK parameter for GK-aware
  block traversal: `IF( I .EQ. LBBEGK(IBB) )`)
- So the scaffolding is there, but the core NCD logic is missing

This means implementing NCD is not just "fixing a bug" — it's **completing unfinished work**
that Willems planned but never finished. The infrastructure (RWORK slots, block traversal)
is already there waiting for the NCD computation to be filled in.

### 1. NCD-aware shift selection in dlaxrf.f / dlaxrf_selshf.f

**Status**: NOT IMPLEMENTED (infrastructure exists but NCD computation missing).

When the representation tree splits a cluster and creates child representations via
shifting, the shift must satisfy the NCD condition (Definition 4.6 in Willems-Lang 2012):
after shifting by μ, the child L⁺D⁺(L⁺)^T must have a nearly constant diagonal:
  max(D⁺) - min(D⁺) ≤ ξ_GK    where ξ_GK = 32·n·eps

DSTEXR picks shifts that minimize element growth (which is good for generic tridiagonal
problems) but doesn't check whether the shift preserves GK structure. When NCD fails,
the child eigenvectors lose GK structure, and U/V extraction produces garbage.

**What to implement**: After computing a candidate shift μ and the resulting L⁺D⁺(L⁺)^T,
check the NCD condition. If it fails:
- Try alternative shift candidates (different Ritz values, bisection midpoints)
- Use blocked (2×2) factorizations (Section 4.3) which control element growth AND
  can maintain NCD through careful pivot selection

**Where to add it**: `dlaxrf_seltw.f` already has RWORK(5)/RWORK(6) reserved for NCD output.
Fill in the NCD computation: after computing shifted LDL^T, set NCD = max(D⁺) - min(D⁺).
Then in `dlaxrf_selshf.f`, use NCD as a selection criterion alongside element growth.

**Key files**: `dlaxrf_selshf.f` (shift candidate generation), `dlaxrf.f` (child rep
construction), `dlaxrf_seltw.f` (twist selection — has NCD placeholder), `dlaxrs.f`
(blocked shift factorizations)

### 2. NCD checking in representation tree traversal (dlaxrv.f)

**Status**: NOT IMPLEMENTED.

The tree traversal in `dlaxrv.f` builds child representations and computes eigenvectors.
At each node, it should verify NCD before accepting the child representation. If NCD
fails, it should request a different shift from `dlaxrf.f`.

**What to implement**: After each shift produces a child LDL^T, compute
max(D)-min(D) and compare against ξ_GK. If NCD fails, try alternative shifts or
fall back to blocked factorizations.

### 3. GK-structure-preserving eigenvector computation (dlaxrg.f)

**Status**: PARTIALLY IMPLEMENTED (root level only via dlaxre.f fix).

When GK-form is active AND NCD holds through the tree, eigenvectors naturally have
GK structure: the q vector has the pattern [v₁, u₁, v₂, u₂, ...] where u and v
components are cleanly separated. Currently this only works at the root because
NCD isn't enforced at deeper levels.

The GK-form block infrastructure already exists in `dlaxrm_stat*.f` files (LBBEGK
parameter, GK-aware block traversal). This infrastructure just needs the NCD checks
upstream to actually produce GK-form child representations.

### 4. The 9 shift cases from Willems thesis (Section 4.2)

**Status**: NOT IMPLEMENTED.

Algorithm 4.1 describes 9 distinct cases for the first-level shift from the TGK root,
depending on where the eigenvalue cluster falls relative to the spectrum. Each case
has a specific shift formula that maintains GK structure. DSTEXR uses a generic shift
strategy that doesn't distinguish these cases.

### 5. Alternative: C++ reimplementation of MR³ with NCD

If modifying the Fortran proves too complex (45 interdependent files), an alternative
is to reimplement the core MR³ loop in C++ directly in `bidiag_svd.h` with NCD checks
built in from the start. The key routines needed:
- LDL^T factorization with shift (dqds-based)
- Sturm count / bisection for eigenvalue refinement
- Twisted factorization for eigenvector computation
- NCD condition check after each shift
- The 9-case shift selection for TGK

This is a significant effort (~500-1000 lines of C++) but avoids the complexity of
modifying deeply interconnected Fortran code. The prior agent's `mr3_tridiag.py`
(pure Python MR³, ~1050 lines) could serve as a reference, though it uses absolute
accuracy (bisection) rather than relative accuracy (dqds).

## Which Fortran Files Need Changes

| File | What | Priority |
|------|------|----------|
| `dlaxrf_seltw.f` | Fill in NCD computation (RWORK(5/6) placeholder exists) | **HIGH** |
| `dlaxrf_selshf.f` | Use NCD as shift selection criterion | **HIGH** |
| `dlaxrf.f` | Try alternative shifts when NCD fails | **HIGH** |
| `dlaxrv.f` | Add NCD verification in tree traversal | **HIGH** |
| `dlaxrs.f` | Extend blocked factorizations for NCD | MEDIUM |
| `dlaxrb_clssfy.f` | AVGAPFAC=0.1→0.3 fix (Bug #1, may already be applied) | LOW |
| `dlaxre.f` | Already fixed (GK-form re-enabled at root) | DONE |

## DO NOT: Post-processing and O(n³) Approaches

**Post-processing is a dead end. O(n³) algorithms are a dead end. Do not go there.**

The previous agent spent most of its 19 iterations tweaking post-processing:
normalization, one-sided recovery, sign fixes, basis completion, chunked MGS thresholds.
The result: pass rate went up but the score stayed at 5.0 because scaling broke.

Specifically, these are PROVEN FAILURES:
- **Full MGS reorthogonalization** → O(n³), score = 5
- **Basis completion for zero columns** → O(n³) for large n, score = 5
- **One-sided recovery U=Bv/σ for all vectors** → doesn't fix extraction failures, helps marginally
- **Larger chunk sizes in chunked MGS** → breaks O(n²) for large clusters
- **Aggressive cluster detection + reortho** → either O(n³) or insufficient
- **Any DBDSQR/HGBSVD fallback** → O(n³) or hybrid, score = 5

The only path to a high score is fixing the eigensolver so it produces eigenvectors
with GK structure in the first place. Then extraction works cleanly and O(n²) post-processing
suffices.

## Starting Point Recommendation

Use `version_013_pass138.h` as the starting `bidiag_svd.h` — it's the best O(n²)-compliant
version (138/289, score=64.3, 4.29x scaling).

The key insight: **improving beyond 138/289 requires fixing the Fortran eigensolver,
not the C++ post-processing.** The post-processing has been exhaustively optimized.
The remaining failures are caused by DSTEXR producing eigenvectors without GK structure,
which happens because NCD is not enforced in the representation tree below the root.

## Diff Between v013 and v017

v017 adds two things over v013:
1. **Basis completion** (lines 109-178 in v017): O(n³) orthogonal completion for zero-component columns
2. Removes the `n≤80` gate on one-sided recovery

These additions fix 69 more tests but destroy O(n²) scaling. The goal is to make
DSTEXR produce eigenvectors good enough that these O(n³) band-aids aren't needed.

## Test Infrastructure Notes

- Build: `g++ -std=c++17 -O2 -Isrc/clapack -I<program_dir> -o evaluate evaluate.cpp -Llib -lxmr_c -lm`
- Run: `./evaluate "<STCollection_path>" 200` for 379 tests
- Pass thresholds: residual ≤ 7.0, orthoU ≤ 5.0, orthoV ≤ 5.0 (in units of n·eps)
- Scaling hard gate: worst doubling ratio > 5.0x → score capped at 5
- The `xmr_disable_gkform_` toggle is already in the compiled `lib/libxmr_c.a`
- To test Fortran changes: copy XMR sources to workspace, edit, rebuild with gfortran, relink
