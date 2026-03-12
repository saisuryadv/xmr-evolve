## THE PROBLEM

Given an upper bidiagonal matrix B (n×n) with diagonal entries d[0..n-1] and
superdiagonal entries e[0..n-2], compute the **singular value decomposition**:
  B = U · Σ · V^T
where Σ = diag(σ_1 ≥ σ_2 ≥ ... ≥ σ_n ≥ 0), U and V are orthogonal.

**The challenge**: Do this in **O(n²) worst-case time** while achieving
LAPACK-quality accuracy (residual ≤ 7·n·ε·‖B‖, orthogonality ≤ 5·n·ε).
This is an **open problem** — DBDSQR (QR iteration) achieves the accuracy but
is O(n³); MR³-based methods achieve O(n²) but fail accuracy on adversarial inputs.

**Why it's hard**: The MR³ algorithm computes eigenvectors of the 2n×2n
Golub-Kahan (TGK) tridiagonal matrix, then extracts u/v from the eigenvectors.
This extraction requires the eigenvector subspace to have "GK structure" —
a property that shifting (core to MR³) can destroy. The NCD (Nearly Constant
Diagonal) condition ensures GK structure is preserved, but verifying and
maintaining NCD through the representation tree is the unsolved challenge.

You are debugging the XMR (DSTEXR) implementation — the Willems code that
implements Algorithm 4.1 from Willems-Lang 2012.

---

## BEFORE ANYTHING ELSE: ENTER PLAN MODE

**Your very first action must be to enter plan mode** (use the EnterPlanMode tool).
Do ALL of your reading, analysis, and planning inside plan mode. Only exit plan mode
when you have a complete implementation plan and are ready to write code.

## STEP 0: READ PRIOR AGENT RESEARCH (MANDATORY — DO THIS BEFORE EVERYTHING)

A previous agent already spent an entire session on this exact task. It made real progress
but also proved several approaches are dead ends. You MUST read these files FIRST:

1. **`recovered/BRIEFING.md`** — **MUST READ. DO NOT SKIP.** Contains:
   - What was done (GK-form root fix, 19 post-processing iterations)
   - What failed and why (post-processing, O(n³) approaches — all dead ends)
   - **THE REMAINING WORK** — exactly which Fortran files need NCD implementation,
     where the placeholder code already exists (RWORK(5/6) in dlaxrf_seltw.f),
     and what Algorithm 4.1 requires at each tree level
   - Why SCORE matters more than pass rate (207/289 pass = score 5.0 vs 138/289 = score 64.3)
2. **`recovered/SUMMARY.md`** — All 19 code versions with pass rates and scores
3. **`recovered/version_013_pass138.h`** — Best version: 138/289, score=64.3, 4.29x scaling

**Your plan should be built around the remaining work described in BRIEFING.md.**
The prior agent identified that NCD was PLANNED but NEVER IMPLEMENTED in the original
XMR code — the infrastructure (RWORK slots, LBBEGK block traversal) exists but the
actual NCD computation was never written. Your job is to complete this unfinished work.

**CRITICAL LESSONS FROM PRIOR AGENT:**
- Post-processing and O(n³) algorithms are DEAD ENDS — do not go there
- **SCORE is what matters, not pass rate.** High pass rate + bad scaling = score of 5
- The GK-form fix in `dlaxre.f` (root level) is ALREADY DONE — don't redo it
- The remaining work is NCD-aware shifts in deeper tree levels (dlaxrf_selshf.f, dlaxrv.f)

The prior agent's full conversation transcript (for deep context if needed):
`/Users/saisurya/.claude/projects/-Users-saisurya-MRRR-xmr-evolve/b05f34bb-76b2-42f5-8184-65ab5d3937ea.jsonl`

## STEP 0.5: READ THE PAPERS

Read these PDFs yourself (using the Read tool directly, NOT via subagents):

1. knowledge/paper_willems_lang_2012_mr3gk.pdf — THE KEY PAPER: Algorithm 4.1, Theorem 4.5, NCD
2. knowledge/paper_willems_lang_2013_framework.pdf — Five requirements, XMR implementation
3. knowledge/paper_willems_thesis_2010.pdf — Full XMR details (Ch. 4-6)

Then read these failure analysis resources:
4. knowledge/slides_bebop_mr3_bsvd.pdf — 88 slides with failure examples
5. knowledge/communication_willems_xmr_failure.pdf — Known XMR failure cases
6. knowledge/paper_marques_demmel_2020.pdf — Modern bSVD failures

Then read these knowledge files:
- knowledge/PRIOR_APPROACHES.md — 12 approaches already tried
- knowledge/INDEX.md — All known bugs, numerical thresholds
- knowledge/xmr_code_documentation.md — XMR Fortran codebase documentation

Only after you have read and understood these, proceed to the mission below.

---

## YOUR MISSION: DEBUG AND FIX THE XMR (DSTEXR) IMPLEMENTATION

According to **Theorem 4.5** (Willems-Lang 2012), if the representation tree
satisfies five requirements (RRR, ELG, RELGAPS, SHIFTREL, GETVEC) AND each
node's invariant subspace nearly has GK structure (via NCD), then:
- Orthogonality: O(n·eps)
- Norm deviation: O(n·eps)
- Residuals: O(||B||·n·eps)

**The theory says Algorithm 4.1 should work. But the existing implementation (DSTEXR)
is only a generic eigensolver — it's missing the NCD-aware bidiagonal SVD layer that
Algorithm 4.1 requires. Your job is to implement those missing pieces.**

This is an implementation task guided by existing theory. Algorithm 4.1 and Theorem 4.5
describe exactly what's needed — the challenge is making it work in finite-precision
arithmetic.

**CRITICAL CONSTRAINTS:**
- **SCORE matters, not pass rate.** The prior agent got 207/289 pass rate but scored 5.0
  (minimum!) because O(n³) scaling. 138/289 with O(n²) scaling scored 64.3. Optimize for SCORE.
- Do NOT add post-processing or O(n³) algorithms. The prior agent exhaustively tried these
  across 19 iterations — normalization, one-sided recovery, basis completion, chunked MGS,
  sign fixes — and proved they are dead ends. Post-processing cannot fix eigenvectors that
  fundamentally lack GK structure. The problem is INSIDE the eigensolver.
- Do NOT use hybrid/fallback approaches. No HGBSVD + DSTEXR, no DBDSQR fallback,
  no "try X, if it fails try Y". Fix the core algorithm.
- The GK-form fix at the root level (`dlaxre.f`) is ALREADY DONE. What's missing is
  NCD-aware shifts at DEEPER tree levels (`dlaxrf_selshf.f`, `dlaxrv.f`).
  This is the actual remaining work.

---

## STEP 1: READ THE THEORY (MANDATORY — before any code)

Read these using the Read tool (PDFs, NOT via subagents):

**Algorithm & Theorem:**
```
knowledge/paper_willems_lang_2012_mr3gk.pdf    — Algorithm 4.1, Theorem 4.5, NCD
knowledge/paper_willems_lang_2013_framework.pdf — Five requirements, XMR implementation
knowledge/paper_willems_thesis_2010.pdf         — Full XMR details (Ch. 4-6)
```

**Failure analysis:**
```
knowledge/slides_bebop_mr3_bsvd.pdf             — 88 slides with failure examples
knowledge/communication_willems_xmr_failure.pdf — Known XMR failure cases
knowledge/paper_marques_demmel_2020.pdf         — Modern bSVD failures
```

**MR³ foundations (read as needed):**
```
knowledge/paper_dhillon_thesis_1997.pdf         — Original MR³
knowledge/paper_parlett_dhillon_2000.pdf         — RRR theory
knowledge/paper_dhillon_parlett_2004.pdf         — Orthogonal eigenvectors
knowledge/paper_dhillon_parlett_vomel_2005_glued.pdf — Glued matrices, failure modes
```

Write `/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/THEORY_NOTES.md` summarizing:
- The five requirements and what each guarantees
- The NCD condition: what it checks, how it's enforced, what threshold
- GK structure: definition, why extraction needs it, when it breaks
- Algorithm 4.1 step by step: root rep → first shift → rep tree → extraction
- Block factorizations: what they are, how they control element growth

## STEP 1.5: WRITE YOUR PLAN (MANDATORY — before any code or diagnostics)

Now that you've read the papers and theory, write a detailed implementation plan to
`/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/PLAN.md` covering:

1. **What you learned**: Key findings from papers — which Algorithm 4.1 requirements
   are missing from DSTEXR, what NCD is and why it matters
2. **Specific changes**: Which files you'll modify (`dlaxre.f`, `dlaxrv.f`, etc.),
   what exact changes you'll make, and which Algorithm 4.1 requirement each change addresses
3. **Diagnosis plan**: What diagnostics you'll run to identify which bugs trigger on failing matrices
4. **Testing strategy**: How you'll verify each change (rebuild library, run tests, compare)
5. **Expected impact**: Which failing matrices you expect each change to fix, and why
6. **Fallback plan**: What you'll try if your primary approach doesn't improve results

**DO NOT skip this step.** Write the plan file BEFORE running any diagnostics or writing any C++ code.

## STEP 2: READ THE REFERENCE CODE

**Your starting code** (the file you'll edit): `/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/program/bidiag_svd.h`
This is the self-contained C++ TGK+DSTEXR implementation. Read it first.

**C++ baselines** (for comparison):
- `src/bidiag_tgk_stexr.h` + `src/bidiag_tgk_common.h` — original DSTEXR wrapper
- `src/bidiag_tgk_stemr.h` — same TGK, different eigensolver
- `src/bidiag_hgbsvd.h` — Großer-Lang coupling approach
- `src/bidiag_dbdsqr.h` — QR iteration (gold standard accuracy)

**XMR Fortran** (for understanding dstexr_ internals):
`../bidiag-algo/MRRR Resources/Code/tdsolver/xmr/SRC/O/`
- `dstexr.f` — Main driver (generic improved MR³ eigensolver)
- `dlaxrv.f` — Representation tree traversal
- `dlaxrf.f` — Child shift computation
- `dlaxrt.f` — Twisted factorization (method 2 commented out — bug)
- See `knowledge/xmr_code_documentation.md` for full file listing

**IMPORTANT — What DSTEXR is and isn't:**

DSTEXR (Willems XMR) is an **improved generic MR³ tridiagonal eigensolver**. It does NOT
implement Algorithm 4.1 for bidiagonal SVD. It has no awareness that it's being called
on a TGK matrix.

**What DSTEXR implements (improvements over DSTEMR):**
- Blocked shift factorizations (`dlaxrs.f`) — controls element growth
- Better shift candidate selection (`dlaxrf_selshf.f`)
- Better twist index selection (`dlaxrf_seltw.f`)
- Improved singleton/cluster classification (`dlaxrb_clssfy.f`)

**What was coded but deactivated (NOW FIXED by prior agent):**
- GK-form support in `dlaxre.f` — was turned off, prior agent re-enabled it at root level

**What Algorithm 4.1 requires but is NOT in DSTEXR:**
- NCD condition checking during the representation tree traversal
- GK-structure-preserving shift selection (shifts that maintain the nearly-constant-diagonal property)
- The bidiagonal SVD-specific layer that ensures U/V extraction from TGK eigenvectors preserves orthogonality

Your job: implement the missing Algorithm 4.1 pieces in C++ in `bidiag_svd.h` — either
as pre/post-processing around `dstexr_()`, or by reimplementing the full MR³ algorithm
with NCD-aware shifts directly in C++.

**CRITICAL: Do NOT use hybrid/fallback approaches.** Do not combine HGBSVD + DSTEXR,
or DBDSQR fallback, or any "try X, if it fails try Y" pattern. The solution must be
a single DSTEXR-based pipeline. Fix the core extraction problem, don't route around it.

## STEP 3: RUN DIAGNOSTICS

Run TGK+STEXR baseline on all tests:
```bash
g++ -std=c++17 -O2 -Isrc/clapack -I/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/stexr -o /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/stexr/evaluate /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/evaluate.cpp -Llib -lxmr_c -lm
/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/stexr/evaluate "/Users/saisurya/MRRR/xmr-evolve/../bidiag-algo/MRRR Resources/STCollection/DATA" 200
```

Compare against TGK+STEMR (same TGK, different eigensolver):
```bash
g++ -std=c++17 -O2 -Isrc/clapack -I/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/stemr -o /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/stemr/evaluate /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/evaluate.cpp -Llib -lxmr_c -lm
/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/stemr/evaluate "/Users/saisurya/MRRR/xmr-evolve/../bidiag-algo/MRRR Resources/STCollection/DATA" 200
```

For each failing matrix, determine:
1. Does DSTEXR return a nonzero INFO? What does the info code mean?
2. Are eigenvalues correct (compare σ² from DSTEXR vs DSTEMR)?
3. Are full 2n-eigenvectors orthogonal before extraction?
4. Does extraction (U from odd, V from even rows) destroy orthogonality?
5. Which of the five requirements is being violated?

Write diagnostic C++ programs in your workspace to instrument specific failures.

## STEP 4: IDENTIFY AND DOCUMENT BUGS

Write `/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/BUG_REPORT.md` documenting each bug:

**Known bugs to verify (from knowledge/INDEX.md):**
1. `dlaxrb_clssfy.f`: AVGAPFAC=0.1→0.3, AVGTOL depth scaling
2. `dlaxrt.f`: method 2 commented out ("something is wrong")
3. `dlaxre.f`: GK-form at root level — **ALREADY FIXED by prior agent** (see BRIEFING.md)
4. Underflow in qd-transformations → zero vectors without error
5. NaN safeguard → early termination on benign matrices

For each bug: does it actually trigger on our failing matrices?
What happens if you fix it?

## STEP 5: IMPLEMENT ALGORITHM 4.1

**DO NOT just tweak post-processing.** Adding better extraction pairing, cluster tolerances,
one-sided recovery, or reorthogonalization CANNOT fix the core problem: DSTEXR's eigenvectors
lack GK structure because it doesn't check NCD conditions or use GK-structure-preserving shifts.
Post-processing is polishing garbage — 12 prior approaches tried this and all failed
(see knowledge/PRIOR_APPROACHES.md).

**You must implement the actual Algorithm 4.1 from Willems-Lang 2012.** This means:

1. **NCD condition checking** (Definition 4.6): After each shift μ, verify that the shifted
   representation L⁺D⁺L⁺ᵀ has a nearly constant diagonal. If NCD fails, the shift destroys
   GK structure and eigenvectors will be wrong. Use a different shift or a blocked factorization.

2. **GK-structure-preserving shifts** (Section 4): The first shift from the root T_GK must
   be chosen so that the child representation maintains GK structure. This is NOT what DSTEXR
   does — it picks shifts for generic tridiagonal matrices without awareness of the zero-diagonal
   GK structure.

3. **Representation tree with NCD enforcement** (Algorithm 4.1, steps 3-5): At each node,
   check NCD before accepting a child representation. If NCD fails, try alternative shifts
   or use blocked (2×2) factorizations to control element growth.

4. **Extend GK-form to deeper tree levels**: The prior agent already re-enabled GK-form
   at the root level in `dlaxre.f` (DONE). But deeper tree levels still use generic shifts
   that break GK structure. Extend NCD enforcement through `dlaxrf.f` and `dlaxrv.f`.

**Implementation approach: Fix the XMR Fortran code directly.**

You are NOT writing anything from scratch. The XMR codebase (45 Fortran files) already
implements MR³. Your job is to edit the existing Fortran to add the missing Algorithm 4.1
pieces (NCD checks, GK-form re-enablement), then rebuild and link.

```bash
# 1. Copy XMR Fortran to your workspace
mkdir -p /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/xmr_src
cp "../bidiag-algo/MRRR Resources/Code/tdsolver/xmr/SRC/O/"*.f /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/xmr_src/

# 2. Edit the Fortran files — add NCD checks to dlaxrf.f, dlaxrv.f, etc.
#    NOTE: dlaxre.f GK-form fix is ALREADY in lib/libxmr_c.a (prior agent did this)
# Use the Edit tool on /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/xmr_src/dlaxrv.f, dlaxrf.f, dlaxrf_selshf.f, etc.

# 3. Rebuild the library (MUST be named libxmr_fixed.a in xmr_src/)
cd /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/xmr_src && gfortran -O2 -c *.f && ar rcs libxmr_fixed.a *.o

# 4. Test manually with the fixed library
g++ -std=c++17 -O2 -Isrc/clapack -I/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/program -o /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/evaluate /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/evaluate.cpp -L/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/xmr_src -lxmr_fixed -lm

# NOTE: When your session ends, libxmr_fixed.a is automatically copied to lib/libxmr_c.a
# so the evaluator uses your fixed version for scoring.
```

The key files to modify (see `knowledge/xmr_code_documentation.md` for full docs):
- `dlaxre.f` — GK-form at root level **ALREADY FIXED** by prior agent. No changes needed here.
- `dlaxrv.f` — Add NCD condition checking in representation tree traversal
- `dlaxrf.f` — Modify shift selection to be GK-structure-preserving
- `dlaxrs.f` — Blocked shift factorizations (may need NCD-aware extension)
- `dlaxrb_clssfy.f` — Cluster classification (AVGAPFAC fix already applied)

For each change: cite which requirement from Theorem 4.5 it addresses (RRR, ELG, RELGAPS,
SHIFTREL, GETVEC) and explain how NCD is maintained.

**DO NOT reimplement MR³ from scratch in C++.** The Fortran code is already there — edit it.

---

## YOUR WORKSPACE

Your isolated workspace is: `/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002`

```
/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/
  program/bidiag_svd.h  ← YOUR OUTPUT FILE — edit this, final version is read from here
  evaluate.cpp           ← pre-copied, do not modify
```

The workspace has been pre-populated with the current bidiag_svd.h and evaluate.cpp.
You can ONLY write/edit files inside `/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/`. All other files are read-only.

## BUILD & TEST (use these exact commands)

```bash
g++ -std=c++17 -O2 -Isrc/clapack -I/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/program -o /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/evaluate /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/evaluate.cpp -Llib -lxmr_c -lm
/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/evaluate "/Users/saisurya/MRRR/xmr-evolve/../bidiag-algo/MRRR Resources/STCollection/DATA" 200
```

**IMPORTANT**: Always pass the STCollection directory and adv_size=200 to the evaluate binary.
This runs 379 tests (90 patterns × 4 sizes + 19 STCollection). Without these args you only get ~270 tests and miss critical matrices.

## PRIOR RESEARCH

No prior research available yet.

## YOUR CAPABILITIES

You have full Claude Code capabilities: read any file, run any command, search code,
explore the codebase, launch subagents for parallel investigation, plan complex tasks.
The ONLY restriction: Write/Edit are scoped to your workspace directory.

Use subagents liberally — e.g., spawn one to trace a failing matrix through the Fortran
code while you read the relevant paper. Spawn one to compile and test while you investigate.

## BASELINE IMPLEMENTATIONS (compile and compare against these)

Baseline headers are pre-copied into `/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/`. Just compile and run:
```bash
# DBDSQR — gold standard, 379/379 pass, O(n³), score=5 (hard gate)
g++ -std=c++17 -O2 -Isrc/clapack -I/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/dbdsqr -o /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/dbdsqr/evaluate /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/evaluate.cpp -Llib -lxmr_c -lm && /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/dbdsqr/evaluate "/Users/saisurya/MRRR/xmr-evolve/../bidiag-algo/MRRR Resources/STCollection/DATA" 200

# HGBSVD — coupling-based, 174/379 pass, O(n²), score=5 (hard gate)
g++ -std=c++17 -O2 -Isrc/clapack -I/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/hgbsvd -o /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/hgbsvd/evaluate /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/evaluate.cpp -Llib -lxmr_c -lm && /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/hgbsvd/evaluate "/Users/saisurya/MRRR/xmr-evolve/../bidiag-algo/MRRR Resources/STCollection/DATA" 200

# TGK+STEMR — MR³ on TGK, 93/379 pass, O(n²), score=5 (hard gate)
g++ -std=c++17 -O2 -Isrc/clapack -I/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/stemr -o /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/stemr/evaluate /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/evaluate.cpp -Llib -lxmr_c -lm && /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/stemr/evaluate "/Users/saisurya/MRRR/xmr-evolve/../bidiag-algo/MRRR Resources/STCollection/DATA" 200

# TGK+STEXR — XMR on TGK (vanilla=61/379, with GK-form fix=~138) — THE BASELINE YOU'RE IMPROVING
g++ -std=c++17 -O2 -Isrc/clapack -I/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/stexr -o /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/stexr/evaluate /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/evaluate.cpp -Llib -lxmr_c -lm && /Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/stexr/evaluate "/Users/saisurya/MRRR/xmr-evolve/../bidiag-algo/MRRR Resources/STCollection/DATA" 200
```
Your working file at `/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/program/bidiag_svd.h` is never touched by baseline tests.

**Top evolved variants** (if any exist from prior iterations) are also pre-copied into
`/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/baselines/top_variants/`. Each has its own subdir with `bidiag_svd.h` and
`RESEARCH.md`. Compile and run them the same way.

Baseline source code (read to understand their techniques):
- `src/bidiag_dbdsqr.h` — DBDSQR (QR iteration, how it achieves perfect accuracy)
- `src/bidiag_hgbsvd.h` — HGBSVD (Großer-Lang coupling, where/why it returns INFO!=0)
- `src/bidiag_tgk_stemr.h` — TGK+STEMR (extraction logic, where orthogonality breaks)
- `src/bidiag_tgk_stexr.h` — TGK+STEXR (XMR-based, NCD handling)
- `src/bidiag_tgk_common.h` — Shared TGK utilities (sign fixes, one-sided recovery)

## CURRENT BASELINES (379 tests, adv_size=200)
| Algorithm | Pass | Worst Scaling | Score |
|-----------|------|---------------|-------|
| DBDSQR | 379/379 | 8.53x (O(n³)) | **5** (hard gate) |
| HGBSVD | 174/379 | 5.37x | **5** (hard gate) |
| TGK+STEMR | 93/379 | 5.96x | **5** (hard gate) |
| TGK+STEXR (vanilla) | 61/379 | 5.73x | **5** (hard gate) |
| **TGK+STEXR (GK-form fix)** | **~138/289** | **4.29x** | **64.3** (prior agent's best) |

The prior agent already fixed STEXR from 61→138 by re-enabling GK-form at the root level.
Your goal: improve beyond 138 while keeping scaling ≤ 5.0x. Fix DSTEXR's deeper tree levels.

## EVOLVED FILE REQUIREMENTS
- `bidiag_svd.h` MUST be **self-contained** — no `#include` of project headers
  (OpenEvolve copies only this file to a temp dir for evaluation)
- Only standard library + extern "C" LAPACK/BLAS declarations allowed
- `dbdsgr_` (f2c) requires `ftnlen` args: pass `1, 1` at end of each call

## SCORING (0-95)
Exact formula from `src/evaluate.cpp`:
- **50 pts**: `pass_rate * 50` — fraction of tests passing (res <= 7, orthoU <= 5, orthoV <= 5)
- **10 pts**: `max(0, 5.0 - pass_avg_residual) * 2` — accuracy bonus for low residual
- **10 pts**: `max(0, 5.0 - pass_avg_orthoU) * 2` — accuracy bonus for low orthoU
- **10 pts**: `max(0, 5.0 - pass_avg_orthoV) * 2` — accuracy bonus for low orthoV
- **5 pts**: compilation bonus (always awarded if code compiles)
- **10 pts**: O(n^2) scaling bonus (awarded if worst doubling ratio <= 5.0)
- **HARD GATE**: if worst scaling ratio > 5.0x, score is capped at 5 (compilation only)

## EVALUATION TIMEOUTS
- **30-second cumulative timeout**: If total evaluation exceeds 30 seconds, score = 0 and exit.
  DBDSQR (the slowest baseline) finishes all tests in ~5 seconds. Your O(n²) code should be faster.
- **2-second per-test timeout**: Any single test taking >2s is auto-FAIL with metrics = 1e10.
- **Early abort**: If a pattern catastrophically fails at small n (metric > 1000), larger sizes skipped.
- If your algorithm triggers these timeouts, it has a bug (infinite loop, O(n⁴), etc.). Fix it.

## RULES
- NEVER ask questions or stop for confirmation. Just investigate, diagnose, fix, test.
- **O(n²) IS NON-NEGOTIABLE.** If your algorithm's worst scaling ratio > 5x when n doubles,
  **your score is capped at 5** (compilation-only credit). This is a hard gate in the evaluator.
- NEVER add O(n³) operations: no global MGS over all n vectors, no dense n×n matrix multiply,
  no full SVD/eigendecomposition of n×n matrices, no DBDSQR. Chunked operations (chunk ≤ 32) are OK.
- **NO FALLBACKS OR HYBRID ALGORITHMS.** Do NOT combine HGBSVD+DSTEMR, HGBSVD+DSTEXR,
  DBDSQR fallback, or any other "try X, if it fails try Y" hybrid approach. Your algorithm
  must be a single coherent DSTEXR-based pipeline. The goal is to fix and improve the
  DSTEXR path itself — not to paper over its failures with a different algorithm.
  Do NOT use `dbdsgr_` or `dbdsqr_` anywhere in your code.
- You can spawn subagents (Agent tool) for parallel investigation or focused tasks.
- You can create scratch files in `/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/` for notes, intermediate results, etc.
- You can ONLY write/edit files inside `/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/`. All other files are read-only.

## AVAILABLE ROUTINES (pure C, lib/libxmr_c.a)
- dstexr_: Willems XMR improved MR³ (O(n²)) — THE ROUTINE YOU'RE IMPROVING
- dstemr_: LAPACK MR³ tridiag eigensolver (O(n²)) — for comparison only
- dlamch_: machine parameters
- BLAS: dnrm2_, ddot_, dscal_, dcopy_, daxpy_, dgemv_
- dbdsqr_: **DO NOT USE** — O(n³), score capped at 5
- dbdsgr_: **DO NOT USE** — no hybrid/fallback approaches allowed

## KNOWLEDGE BASE
- `knowledge/INDEX.md` — Master reference: algorithms, bugs, test matrices, thresholds
- `knowledge/RESOURCES.md` — Index of all papers, code, test data with paths
- `knowledge/PRIOR_APPROACHES.md` — 12 approaches tried, what worked and failed
- `knowledge/BASELINES.md` — Comparison of 5 baseline algorithms
- `knowledge/EVALUATION.md` — Scoring formula, test patterns, metrics
- `knowledge/willems_lang_2012.md` — Algorithm 4.1 (MR³ on TGK with NCD shifts)
- `knowledge/grosser_lang_2001_hgbsvd.md` — Coupling-based O(n²) SVD
- `knowledge/mr3_foundations.md` — Dhillon 1997, Parlett-Dhillon 2000 (RRR theory)
- `knowledge/marques_2020_bugs_matrices.md` — DBDSVDX bugs, CHKBD analysis
- `knowledge/xmr_code_documentation.md` — Willems XMR 45-file Fortran codebase documentation

## REFERENCE CODE
- Baseline C++ headers: `src/bidiag_dbdsqr.h`, `src/bidiag_hgbsvd.h`, `src/bidiag_tgk_stemr.h`, `src/bidiag_tgk_stexr.h`
- XMR Fortran: `../bidiag-algo/MRRR Resources/Code/tdsolver/xmr/SRC/O/` (45 files, `dstexr.f` is master)
- HGBSVD Fortran: `../bidiag-algo/MRRR Resources/Code/hgbsvd/hgbsvd/v2/` (24 files, `dbdsgr.f` is master)
- LAPACK source: `lapack/SRC/` (dbdsqr.f, dbdsvdx.f, dstemr.f, dlarrv.f)
- STCollection: `../bidiag-algo/MRRR Resources/STCollection/DATA/B_*.dat` (19 matrices)

## DO NOT READ (irrelevant to this project)
- **DO NOT** read `approach_*.py` files — these are Python prototypes from a different experiment
- **DO NOT** read files in `../bidiag-algo/` except Papers, Code, and STCollection directories
- **DO NOT** read `evaluate_alphaevolve.py`, `mr3_tridiag.py`, `problem_statement_alphaevolve.md`
- **DO NOT** read `accuracy_tables.txt`, `results_approach_*.txt`
- **DO NOT** read `.claude/` memory files or auto-memory — these are for a different agent
- Focus on **C++ implementation**, **papers**, **Fortran reference code**, and **knowledge base**

## Evolution History (2 variants so far)

Previous variants saved at: `/Users/saisurya/MRRR/xmr-evolve/evolved_variants/exp_20260311_174919/`
- Code: `/Users/saisurya/MRRR/xmr-evolve/evolved_variants/exp_20260311_174919/v{N}_{score}.h` — Read any variant to see its full code
- Research: `/Users/saisurya/MRRR/xmr-evolve/evolved_variants/exp_20260311_174919/v{N}_RESEARCH.md` — Structured research report

# Evolution Changelog

Each entry = one evaluated variant. Variants saved in this directory.
Full research reports: see v{N}_RESEARCH.md files.

| # | Score | Pass | Avg Res | Avg OrtU | Scaling | Approach | Failing Tests | Time |
|---|-------|------|---------|----------|---------|----------|---------------|------|
| 1 | 5.0 | 178/379 | 0.64 | 0.82 | 7.18x | (no report) | B_bug316_gesdd, B_gg_30_1D-5, B_16_smallsv, B_05_d5eq0, B_11_splits_a, ... | 2026-03-11 17:49 |
| 2 | 5.0 | 178/379 | 0.64 | 0.82 | 8.86x | (no report) | B_bug316_gesdd, B_gg_30_1D-5, B_16_smallsv, B_05_d5eq0, B_11_splits_a, ... | 2026-03-11 18:07 |



---

## NOW: GET INTO PLAN MODE, THEN DEBUG SYSTEMATICALLY

Do NOT jump straight to implementation. You must work in three phases:

**PHASE 1 — PLAN MODE: THEORY & UNDERSTANDING (mandatory, ~25% of your time)**
Enter plan mode. Read the papers, reference code, and baselines.
Write `/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/THEORY_NOTES.md` summarizing:
- Algorithm 4.1 step by step
- The five requirements and what each guarantees
- NCD condition and GK structure
- Where the implementation could diverge from the theory

Do NOT write any C++ code during Phase 1. Only read, analyze, and write notes.

**PHASE 2 — DIAGNOSIS (~35%)**
Run failing tests. Trace through code. Create synthetic test matrices.
Write `/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/BUG_REPORT.md` documenting each divergence from theory.

**PHASE 3 — FIX & VERIFY (~40%)**
Fix bugs. Compile. Test. Compare against DSTEXR baseline.
Write `/Users/saisurya/MRRR/xmr-evolve/scratch/exp_20260311_174919/agent_5671_0002/RESEARCH.md` with final results.

### RESEARCH.md Template

```markdown
# XMR Debug Report: [One-line summary of bugs found and fixed]

## Abstract
[2-3 sentences: what bugs you found, what you fixed, key result (pass rate, score)]

## Theory Summary
[Brief summary of Algorithm 4.1, Theorem 4.5, and the five requirements.
What does the theory guarantee and under what conditions?]

## Bug Analysis
[For each bug found:
- What symptom was observed (which tests fail, what metrics blow up)
- What requirement/condition is violated
- Root cause in the implementation (which Fortran routine, which code path)
- How it diverges from what the theory expects
- Your fix and why it restores the theoretical guarantee]

## Synthetic Test Cases
[Small matrices you constructed to isolate specific numerical issues.
Include the matrix values and what they demonstrate.]

## Experiments
[Test results: pass rate, score, avg residual, avg orthoU, avg orthoV.
Comparison table: DSTEXR baseline → your fix. Which tests improved? Regressed?]

## Remaining Issues
[What tests still fail? Why? What numerical phenomena does your fix not address?]

## Future Work
[What should be investigated next? What deeper issues remain?]
```

The implementation SHOULD work according to Theorem 4.5. Find why it doesn't and fix it.

## CRITICAL MINDSET RULES

**DO NOT resort to hacks, workarounds, or fallback strategies.**
- Do NOT add HGBSVD, DBDSQR, B^TB eigsolve, or any other algorithm as a fallback.
  The goal is to make DSTEXR/Algorithm 4.1 itself work correctly.
- Do NOT add heuristic thresholds, ad-hoc reorthogonalization, or band-aids.
  If orthogonality is bad, find the ROOT CAUSE in the algorithm — which requirement
  is violated? Which shift destroyed NCD? Which factorization has element growth?
**DO find the actual numerical issues and fix them at the source.**
- Read the theory (Theorem 4.5, five requirements, NCD). Understand the guarantees.
  Study the PROOFS — understand not just what the theorem states but WHY it holds.
  The proof reveals which intermediate quantities must be bounded, which cancellations
  are relied upon, and where finite-precision arithmetic can violate the assumptions.
- Read the XMR Fortran code. Trace through failing cases step by step.
- When something fails, ask: which requirement broke? Why? What floating-point
  issue caused it? Is it underflow? Element growth? Wrong shift selection?
  Premature cluster splitting? Missing NCD check?
- Create your own small synthetic test matrices that isolate the specific numerical
  issue you're investigating. For example: if you suspect element growth from shifting
  near-zero eigenvalues, construct a 6x6 bidiagonal with that exact structure and
  trace every step of the algorithm on it.
- Understand the gap between theory and implementation deeply. The proofs assume
  exact arithmetic in certain steps, or assume bounds on intermediate quantities
  that may not hold in practice. Identify exactly WHERE the proof's assumptions
  break in floating-point. For example:
  - Does a dqds transformation accumulate relative error beyond what the proof allows?
  - Does a shift selection rely on a gap estimate that is inaccurate for this spectrum?
  - Does element growth in LDL^T exceed the bound assumed by the NCD check?
  - Does the twisted factorization lose accuracy when the twist index is suboptimal?
- Once you understand the numerical issue, fix it in a way that restores the
  theoretical guarantee. Search the literature and online resources for known
  techniques to handle the specific numerical problem you've identified — there
  may be better shift strategies, more stable factorizations, or tighter error
  bounds published in other papers or numerical computing resources.

**Our goal is to bring the theory to fruition through a strong implementation.**
The theory (Algorithm 4.1, Theorem 4.5) describes what's needed. DSTEXR provides a
strong foundation (improved MR³) but is missing the NCD-aware bidiagonal SVD layer.
Implement it.
