"""
Claude Agent SDK adapter for OpenEvolve

Single-phase agentic pipeline: one agent with full tool access and
directory-scoped write permissions via the SDK's can_use_tool callback.

Each agent call gets an isolated workspace:
  scratch/agent_NNN/
    program/bidiag_svd.h  ← the output file (agent edits this, we read it back)
    evaluate.cpp           ← pre-copied for compilation
    evaluate               ← compiled binary

Write permissions restricted to scratch/agent_NNN/ only.
Read permissions: unrestricted (read anything from anywhere).
Subagent spawning: enabled (Agent tool).

Result capture: we read scratch/agent_NNN/program/bidiag_svd.h after the agent
finishes, instead of parsing ```cpp blocks from conversational output.
"""

import json
import logging
import os
import re
import shutil
from typing import Any, Dict, List, Optional

from openevolve.llm.base import LLMInterface

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Lazy import
# ---------------------------------------------------------------------------
_agent_sdk = None


def _ensure_sdk():
    global _agent_sdk
    if _agent_sdk is None:
        try:
            import claude_agent_sdk as sdk
            _agent_sdk = sdk
        except ImportError:
            raise ImportError(
                "claude-agent-sdk is required. Install with:\n"
                "  pip install claude-agent-sdk\n"
                "Also ensure the Claude Code CLI is available:\n"
                "  npm install -g @anthropic-ai/claude-code"
            )
    return _agent_sdk


# ---------------------------------------------------------------------------
# Agent settings persistence (survives multiprocessing serialization)
# ---------------------------------------------------------------------------
# OpenEvolve serializes LLMModelConfig via dataclasses.asdict() which drops
# custom attributes. We persist agent settings to a JSON file that the
# factory function reads in the worker process.

_AGENT_SETTINGS_FILE = None


def save_agent_settings(scratch_dir: str, settings: dict) -> str:
    """Save agent settings to a JSON file. Call this before starting OpenEvolve.

    Returns the path to the settings file.
    """
    global _AGENT_SETTINGS_FILE
    os.makedirs(scratch_dir, exist_ok=True)
    path = os.path.join(scratch_dir, ".agent_settings.json")
    # Filter out non-serializable values
    serializable = {}
    for k, v in settings.items():
        if isinstance(v, (str, int, float, bool, list, dict, type(None))):
            serializable[k] = v
    with open(path, "w") as f:
        json.dump(serializable, f)
    _AGENT_SETTINGS_FILE = path
    return path


def _load_agent_settings(scratch_dir: str = None) -> dict:
    """Load agent settings from the JSON file."""
    global _AGENT_SETTINGS_FILE

    # Try the cached path first
    if _AGENT_SETTINGS_FILE and os.path.exists(_AGENT_SETTINGS_FILE):
        with open(_AGENT_SETTINGS_FILE) as f:
            return json.load(f)

    # Try the scratch_dir
    if scratch_dir:
        path = os.path.join(scratch_dir, ".agent_settings.json")
        if os.path.exists(path):
            _AGENT_SETTINGS_FILE = path
            with open(path) as f:
                return json.load(f)

    # Auto-discover: look relative to this module's directory
    module_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(module_dir, "scratch", ".agent_settings.json")
    if os.path.exists(path):
        _AGENT_SETTINGS_FILE = path
        with open(path) as f:
            return json.load(f)

    logger.warning("Agent settings file not found — using defaults")
    return {}


# ---------------------------------------------------------------------------
# Unified prompt template
# ---------------------------------------------------------------------------

# Tag: exploration_prompt — general-purpose "explore and evolve" prompt
EXPLORATION_PROMPT_TEMPLATE = """\
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

You are evolving a C++ implementation in `bidiag_svd.h`.

---

## BEFORE ANYTHING ELSE: READ THE PAPERS

You MUST read these PDFs yourself (using the Read tool directly, NOT via subagents) before
writing ANY code or running ANY tests. Do not delegate this to a subagent. Do not skip this.

Read these 4 must-read papers NOW:
1. knowledge/paper_willems_lang_2012_mr3gk.pdf — THE KEY PAPER: MR³-GK Algorithm 4.1
2. knowledge/paper_grosser_lang_2001_on2.pdf — O(n²) coupling approach
3. knowledge/paper_willems_lang_2013_framework.pdf — MR³ framework, 5 requirements
4. knowledge/paper_marques_demmel_2020.pdf — Modern bSVD, DBDSVDX bugs

Then read these knowledge files:
- knowledge/PRIOR_APPROACHES.md — 12 approaches already tried
- knowledge/INDEX.md — All known bugs, numerical thresholds

Additional papers available in knowledge/ (read as needed):
- paper_demmel_kahan_1990.pdf, paper_dhillon_thesis_1997.pdf, paper_parlett_dhillon_2000.pdf
- paper_dhillon_parlett_2004.pdf, paper_willems_thesis_2010.pdf, slides_bebop_mr3_bsvd.pdf

---

## YOUR WORKSPACE

Your isolated workspace is: `{agent_dir}`

```
{agent_dir}/
  program/bidiag_svd.h  ← YOUR OUTPUT FILE — edit this, final version is read from here
  evaluate.cpp           ← pre-copied, do not modify
```

The workspace has been pre-populated with the current bidiag_svd.h and evaluate.cpp.
You can ONLY write/edit files inside `{agent_dir}/`. All other files are read-only.

## BUILD & TEST (use these exact commands)

```bash
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/program -o {agent_dir}/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm
{agent_dir}/evaluate "{stcoll_dir}" 200
```

**IMPORTANT**: Always pass the STCollection directory and adv_size=200 to the evaluate binary.
This runs 379 tests (90 patterns × 4 sizes + 19 STCollection). Without these args you only get ~270 tests and miss critical matrices.

## PRIOR RESEARCH

{prior_research}

## YOUR CAPABILITIES

You have full Claude Code capabilities: read any file, run any command, search code,
explore the codebase, launch subagents for parallel investigation, plan complex tasks.
The ONLY restriction: Write/Edit are scoped to your workspace directory.

Use subagents liberally — e.g., spawn one to analyze a failing matrix while you read
the relevant paper. Spawn one to compile and test a variant while you investigate another.

## BASELINE IMPLEMENTATIONS (compile and compare against these)

You can compile and run ANY baseline to compare results on specific matrices:

Baseline headers are pre-copied into `{agent_dir}/baselines/`. Just compile and run:
```bash
# DBDSQR — gold standard, 379/379 pass, O(n³), score=5 (hard gate)
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/dbdsqr -o {agent_dir}/baselines/dbdsqr/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm && {agent_dir}/baselines/dbdsqr/evaluate "{stcoll_dir}" 200

# HGBSVD — coupling-based, 174/379 pass, O(n²), score=5 (hard gate)
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/hgbsvd -o {agent_dir}/baselines/hgbsvd/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm && {agent_dir}/baselines/hgbsvd/evaluate "{stcoll_dir}" 200

# TGK+STEMR — MR³ on TGK, 93/379 pass, O(n²), score=5 (hard gate)
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/stemr -o {agent_dir}/baselines/stemr/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm && {agent_dir}/baselines/stemr/evaluate "{stcoll_dir}" 200

# TGK+STEXR — XMR on TGK (vanilla=61/379, with GK-form fix=~138)
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/stexr -o {agent_dir}/baselines/stexr/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm && {agent_dir}/baselines/stexr/evaluate "{stcoll_dir}" 200
```
Your working file at `{agent_dir}/program/bidiag_svd.h` is never touched by baseline tests.

**Top evolved variants** (if any exist from prior iterations) are also pre-copied into
`{agent_dir}/baselines/top_variants/`. Each has its own subdir with `bidiag_svd.h` and
`RESEARCH.md`. Compile and run them the same way:
```bash
# Example: test top variant v3_75
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/top_variants/v3_75 -o {agent_dir}/baselines/top_variants/v3_75/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm
{agent_dir}/baselines/top_variants/v3_75/evaluate
```

Baseline source code (read to understand their techniques):
- `src/bidiag_dbdsqr.h` — DBDSQR (QR iteration, how it achieves perfect accuracy)
- `src/bidiag_hgbsvd.h` — HGBSVD (Großer-Lang coupling, where/why it returns INFO!=0)
- `src/bidiag_tgk_stemr.h` — TGK+STEMR (extraction logic, where orthogonality breaks)
- `src/bidiag_tgk_stexr.h` — TGK+STEXR (XMR-based, NCD handling)
- `src/bidiag_tgk_common.h` — Shared TGK utilities (sign fixes, one-sided recovery)

## YOUR WORKFLOW

You have a single unified session with no budget limit. Research deeply, implement,
test, iterate, and document — in whatever order makes sense.

### STEP 0a: READ PRIOR AGENT RESEARCH (MANDATORY — DO THIS FIRST)

A previous agent already made significant progress. Read these files BEFORE reading papers:

1. **`recovered/BRIEFING.md`** — What was done, what worked, what's left to do
2. **`recovered/SUMMARY.md`** — All 19 code versions with pass rates and scores
3. **`recovered/version_013_pass138.h`** — Best version: 138/289, score=64.3, 4.29x scaling (O(n²) compliant)

Do NOT repeat the prior agent's post-processing work. It has been exhaustively optimized.
Post-processing and O(n³) approaches are dead ends (prior agent proved this — 207/289 pass
rate but score=5.0 because scaling was 19.5x). Focus on the Fortran eigensolver internals.

### STEP 0b: READ THE PAPERS (MANDATORY)

This is not optional. You MUST read these papers before writing any code, running any
tests, or even looking at the source code. The papers contain the algorithms, theorems,
and mathematical foundations you need. Without reading them, you will waste time
reinventing approaches that have already been tried and failed.

Read these papers using the Read tool (they are PDFs in the knowledge/ directory).

**Must-read (read ALL of these before writing any code):**
```
knowledge/paper_willems_lang_2012_mr3gk.pdf    — THE KEY PAPER: MR³-GK Algorithm 4.1, NCD, GK structure
knowledge/paper_grosser_lang_2001_on2.pdf       — O(n²) coupling approach for bSVD
knowledge/paper_willems_lang_2013_framework.pdf — MR³ framework: 5 requirements, XMR implementation
knowledge/paper_marques_demmel_2020.pdf         — Modern bSVD: DBDSVDX bugs, CHKBD matrices
```

**Foundation papers (read if you need deeper understanding):**
```
knowledge/paper_demmel_kahan_1990.pdf              — Accurate singular values (foundational)
knowledge/paper_dhillon_thesis_1997.pdf            — Original MR³ algorithm
knowledge/paper_parlett_dhillon_2000.pdf           — Relatively robust representations
knowledge/paper_dhillon_parlett_2004.pdf           — Orthogonal eigenvectors and relative gaps
knowledge/paper_dhillon_parlett_vomel_2005_glued.pdf — Glued matrices and MRRR
knowledge/paper_willems_thesis_2010.pdf            — Willems thesis: XMR full details
```

**Additional references:**
```
knowledge/paper_barlow_2002.pdf                    — More accurate bidiagonal reduction
knowledge/paper_lapack_wn163_2005.pdf              — LAPACK Working Note 163
knowledge/paper_lapack_wn166_2005.pdf              — LAPACK Working Note 166
knowledge/paper_symmetric_eigenproblems_bsvd_2005.pdf — Eigenproblems from bidiagonal SVD
knowledge/paper_demmel_marques_2008.pdf            — LAPACK eigensolver performance
```

**Slides & communications:**
```
knowledge/slides_bebop_mr3_bsvd.pdf                — BeBOP talk: MR³ for bSVD (88 slides)
knowledge/slides_willems_thesis.pdf                — Willems thesis slides (16 slides)
knowledge/communication_willems_xmr_failure.pdf    — XMR failure cases discussion
```

After reading the papers, also read:
- `knowledge/PRIOR_APPROACHES.md` — 12 approaches already tried, what worked and failed
- `knowledge/INDEX.md` — All known bugs, test matrix formulas, numerical thresholds

Take notes on: What is Algorithm 4.1? What is the NCD condition? What is GK structure?
What are the five MR³ requirements? Why does coupling fail at deep recursion? What are
block factorizations and how do they reduce element growth?

Only after you have read and understood the papers, proceed to Step 1.

### STEP 1: IDENTIFY FAILING TESTS

Compile and run the current code. Group failures by mode (INFO!=0, residual blow-up,
orthoU/V blow-up). Note WHICH specific matrices fail and at WHICH sizes.

### STEP 2: COMPARE AGAINST BASELINES AND TOP VARIANTS

For the failing matrices, run the baselines AND top evolved variants on the same tests:
   - Which baselines/variants pass where yours fails? What do they do differently?
   - Which baselines/variants fail on the same matrices? This reveals fundamental hardness.
   - Read the source code (and RESEARCH.md for variants) to understand their techniques.
   - If a top variant already solves a failure mode, DON'T reinvent — build on its approach.

### STEP 3: CONNECT FAILURES TO THE PAPERS

For the most impactful failure mode, go back to the papers you read in Step 0:
   - **What algebraic property breaks?** Orthogonality, norm deviation, residual
     coupling, representation quality (NCD/RRR)?
   - **What structural property of the failing matrix triggers it?** Clustered
     singular values, extreme grading, near-zero entries, specific condition pattern?
   - **Which paper addresses this?** Read the relevant theorem/algorithm. Understand
     the assumptions and where they break for the failing case.
   - **Read the reference Fortran code.** See how XMR (dstexr.f) or HGBSVD (dbdsgr.f)
     handle the same case. What technique do they use? Can it be adapted?

### STEP 4: FORMULATE A PRINCIPLED FIX

Design a solution that addresses the root cause, not the symptom. Know which
theorem/property guarantees correctness. Reference the specific paper and theorem.

### STEP 5: VERIFY NOVELTY BEFORE IMPLEMENTING

This is critical. Before writing code,
   check whether your proposed approach is ALREADY implemented by an existing
   baseline or prior variant:
   - **Spawn a subagent** to search all RESEARCH.md files, baseline source code,
     knowledge/PRIOR_APPROACHES.md, and the relevant papers for the technique you're
     proposing. Is it the same as what HGBSVD does? Is it what a prior variant tried?
   - **Example of what NOT to do**: Proposing "use B^TB eigensolve with coupling" when
     HGBSVD (Großer-Lang) already does exactly this. Read `src/bidiag_hgbsvd.h` and
     `knowledge/grosser_lang_2001_hgbsvd.md` first.
   - **Your approach must be genuinely different** from all existing baselines and prior
     variants. Explain specifically what is new and why prior attempts at the same idea
     didn't work.

### STEP 6: PREDICT FAILURES BEFORE IMPLEMENTING

Before writing code:
   - **Which matrices will your fix NOT help?** Identify cases where your approach's
     assumptions break (e.g., "this assumes well-separated singular values, so
     `two_clusters` at n=400 will still fail because...").
   - **Could your fix REGRESS anything?** If you change the fallback logic, do
     previously-passing HGBSVD cases still pass? If you change the TGK extraction,
     do easy matrices regress?
   - **What's the complexity?** Count the operations. Will it stay O(n²)?

### STEP 7: IMPLEMENT

Edit `{agent_dir}/program/bidiag_svd.h`.

### STEP 8: COMPILE & TEST

Build and run. Compare against the parent variant AND against baselines:
   - Did you improve on the parent? (pass rate, score, specific matrices)
   - Did you overcome any baseline bugs/failures? (e.g., HGBSVD INFO!=0 cases)
   - Did anything regress? If so, diagnose why and fix.
   - Keep iterating until pass rate is maximized.

### STEP 9: VERIFY IMPROVEMENT

Run your final version side-by-side with:
   - The **baselines** (DBDSQR, HGBSVD, TGK+STEMR) on failing matrices
   - The **top evolved variants** (in `baselines/top_variants/`) on the same matrices
   Confirm you actually beat them, not just match. If a top variant already handles
   a case you're working on, study HOW it does it (read its code and RESEARCH.md)
   and build on that approach rather than reinventing. Report a comparison table.

### STEP 10: DOCUMENT
Write `{agent_dir}/RESEARCH.md`. Future agents read this to
    build on your work. Follow this structure:

```markdown
# Variant N: [One-line title of your approach]

## Abstract
[2-3 sentences: what problem you addressed, what you did, key result (pass rate, score)]

## Introduction
[What specific failures motivated this work? Which test patterns fail and why?
What is the parent variant's score and pass rate?]

## Related Work
[What prior variants (from the changelog) attempted similar fixes? What approaches
from knowledge/PRIOR_APPROACHES.md are relevant? Why did they succeed or fail?]

## Proposed Method
[Detailed description of your code changes. Which functions modified? What conditions
added/changed? Include the key formulas/thresholds. If this is a modification of an
existing technique, clearly state what's different.]

### Novelty
[What is new about your approach compared to all prior variants? Be specific.]

## Experiments
[Test results: pass rate, score, avg residual, avg orthoU, avg orthoV.
Comparison table vs parent variant. Which tests improved? Which regressed?
Include actual numbers from ./evaluate output.]

## Limitations
[What tests still fail? Why? What numerical phenomena does your approach not handle?]

## Future Work
[Specific suggestions for the next agent. What should be tried next based on
your findings? What promising directions did you not have budget to explore?]
```

## EVALUATION TIMEOUTS
- **30-second cumulative timeout**: If total evaluation exceeds 30 seconds, score = 0 and exit.
  DBDSQR (the slowest baseline) finishes all 270 tests in ~5 seconds. Your O(n²) code should be faster.
- **2-second per-test timeout**: Any single test taking >2s is auto-FAIL with metrics = 1e10.
- **Early abort**: If a pattern catastrophically fails at small n (metric > 1000), larger sizes skipped.
- If your algorithm triggers these timeouts, it has a bug (infinite loop, O(n⁴), etc.). Fix it.

## RULES
- NEVER ask questions or stop for confirmation. Just implement, compile, test, iterate.
- **O(n²) IS NON-NEGOTIABLE.** If your algorithm's worst scaling ratio > 5x when n doubles,
  **your score is capped at 5** (compilation-only credit). This is a hard gate in the evaluator.
  Do NOT use DBDSQR or any O(n³) routine as a fallback — it will zero your score even if
  every test passes. The entire point is to achieve accuracy WITHOUT O(n³) operations.
- NEVER add O(n³) operations: no global MGS over all n vectors, no dense n×n matrix multiply,
  no full SVD/eigendecomposition of n×n matrices, no DBDSQR. Chunked operations (chunk ≤ 32) are OK.
- You can spawn subagents (Agent tool) for parallel investigation or focused tasks.
- You can create scratch files in `{agent_dir}/` for notes, intermediate results, etc.
- You can ONLY write/edit files inside `{agent_dir}/`. All other files are read-only.

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

## CURRENT BASELINES (379 tests: 90 patterns × 4 sizes + 19 STCollection, adv_size=200)
| Algorithm | Pass | Worst Scaling | Score |
|-----------|------|---------------|-------|
| DBDSQR | 379/379 | 8.53x (O(n^3)) | **5** (hard gate) |
| HGBSVD | 174/379 | 5.37x | **5** (hard gate) |
| TGK+STEMR | 93/379 | 5.96x | **5** (hard gate) |
| TGK+STEXR (vanilla) | 61/379 | 5.73x | **5** (hard gate) |
| TGK+STEXR (GK-form fix) | ~138/289 | 4.29x | **64.3** (prior agent's best) |
| Current hybrid | 213/379 | ~4.7x | 71.4 |

## CRITICAL ALGORITHMIC RULES
- **O(n^2) worst-case required** — no global MGS, no dense n*n multiplies
- The fundamental tradeoff: accuracy vs scaling. DBDSQR has accuracy but O(n^3). TGK+STEMR has O(n^2) but fails 76% of tests.
- TGK eigenvector extraction (U from odd rows, V from even rows) destroys orthogonality unless eigenvectors have GK structure (requires NCD-aware shifts inside the eigensolver)
- One-sided recovery (U = BV/sigma) helps for small sigma but does not fix extraction failures
- Chunked MGS (chunk size 32) is O(n^2) but may miss cross-chunk orthogonality

## AVAILABLE ROUTINES (pure C, lib/libxmr_c.a)
- dbdsqr_: QR iteration bidiag SVD — **O(n³), DO NOT USE** (score capped at 5 if scaling > 5x)
- dstemr_: LAPACK MR³ tridiag eigensolver (O(n²))
- dstexr_: Willems XMR improved MR³ (O(n²))
- dbdsgr_: Großer-Lang coupling SVD (O(n²), needs ftnlen 1,1)
- dlamch_: machine parameters
- BLAS: dnrm2_, ddot_, dscal_, dcopy_, daxpy_, dgemv_

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

{changelog}

---

## NOW: GET INTO PLAN MODE, EXPLORE, THEN EXECUTE

Do NOT jump straight to implementation. You must work in three phases:

**PHASE 1 — PLAN MODE (mandatory, ~40% of your time)**
Enter plan mode. Read the papers, baselines, prior approaches, and failing tests.
Write a detailed analysis document to `{agent_dir}/ANALYSIS.md` covering:
1. Which tests fail on the current code and WHY (group by failure mode)
2. What each paper says about the relevant failure modes (cite specific theorems/sections)
3. What prior approaches tried and why they failed (from PRIOR_APPROACHES.md)
4. What the baselines do differently on these failing cases
5. Your proposed approach — with specific hypotheses about what will fix each failure category
6. Predicted failures — which matrices your approach will NOT help and why
7. Complexity analysis — why your approach is O(n²)

Do NOT write any C++ code during Phase 1. Only read, analyze, and write ANALYSIS.md.

**PHASE 2 — IMPLEMENT (after ANALYSIS.md is complete)**
Now implement the approach described in your analysis. Edit `{agent_dir}/program/bidiag_svd.h`.
Compile and test. Compare against baselines on failing matrices.

**PHASE 3 — ITERATE & DOCUMENT**
If tests fail, diagnose WHY (connect back to your analysis), fix, retest.
Write `{agent_dir}/RESEARCH.md` with final results and comparison tables.

The quality of your ANALYSIS.md determines the quality of your implementation.
Rushing to code without deep analysis produces crappy solutions. Think first, code second.
"""


# Tag: xmr_debug_prompt — focused on debugging the XMR/DSTEXR implementation
# to make Algorithm 4.1 (Theorem 4.5) work correctly.
XMR_DEBUG_PROMPT_TEMPLATE = """\
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

Write `{agent_dir}/THEORY_NOTES.md` summarizing:
- The five requirements and what each guarantees
- The NCD condition: what it checks, how it's enforced, what threshold
- GK structure: definition, why extraction needs it, when it breaks
- Algorithm 4.1 step by step: root rep → first shift → rep tree → extraction
- Block factorizations: what they are, how they control element growth

## STEP 1.5: WRITE YOUR PLAN (MANDATORY — before any code or diagnostics)

Now that you've read the papers and theory, write a detailed implementation plan to
`{agent_dir}/PLAN.md` covering:

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

**Your starting code** (the file you'll edit): `{agent_dir}/program/bidiag_svd.h`
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
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/stexr -o {agent_dir}/baselines/stexr/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm
{agent_dir}/baselines/stexr/evaluate "{stcoll_dir}" 200
```

Compare against TGK+STEMR (same TGK, different eigensolver):
```bash
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/stemr -o {agent_dir}/baselines/stemr/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm
{agent_dir}/baselines/stemr/evaluate "{stcoll_dir}" 200
```

For each failing matrix, determine:
1. Does DSTEXR return a nonzero INFO? What does the info code mean?
2. Are eigenvalues correct (compare σ² from DSTEXR vs DSTEMR)?
3. Are full 2n-eigenvectors orthogonal before extraction?
4. Does extraction (U from odd, V from even rows) destroy orthogonality?
5. Which of the five requirements is being violated?

Write diagnostic C++ programs in your workspace to instrument specific failures.

## STEP 4: IDENTIFY AND DOCUMENT BUGS

Write `{agent_dir}/BUG_REPORT.md` documenting each bug:

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
mkdir -p {agent_dir}/xmr_src
cp "../bidiag-algo/MRRR Resources/Code/tdsolver/xmr/SRC/O/"*.f {agent_dir}/xmr_src/

# 2. Edit the Fortran files — add NCD checks to dlaxrf.f, dlaxrv.f, etc.
#    NOTE: dlaxre.f GK-form fix is ALREADY in lib/libxmr_c.a (prior agent did this)
# Use the Edit tool on {agent_dir}/xmr_src/dlaxrv.f, dlaxrf.f, dlaxrf_selshf.f, etc.

# 3. Rebuild the library (MUST be named libxmr_fixed.a in xmr_src/)
cd {agent_dir}/xmr_src && gfortran -O2 -c *.f && ar rcs libxmr_fixed.a *.o

# 4. Test manually with the fixed library
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/program -o {agent_dir}/evaluate {agent_dir}/evaluate.cpp -L{agent_dir}/xmr_src -lxmr_fixed -lm

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

Your isolated workspace is: `{agent_dir}`

```
{agent_dir}/
  program/bidiag_svd.h  ← YOUR OUTPUT FILE — edit this, final version is read from here
  evaluate.cpp           ← pre-copied, do not modify
```

The workspace has been pre-populated with the current bidiag_svd.h and evaluate.cpp.
You can ONLY write/edit files inside `{agent_dir}/`. All other files are read-only.

## BUILD & TEST (use these exact commands)

```bash
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/program -o {agent_dir}/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm
{agent_dir}/evaluate "{stcoll_dir}" 200
```

**IMPORTANT**: Always pass the STCollection directory and adv_size=200 to the evaluate binary.
This runs 379 tests (90 patterns × 4 sizes + 19 STCollection). Without these args you only get ~270 tests and miss critical matrices.

## PRIOR RESEARCH

{prior_research}

## YOUR CAPABILITIES

You have full Claude Code capabilities: read any file, run any command, search code,
explore the codebase, launch subagents for parallel investigation, plan complex tasks.
The ONLY restriction: Write/Edit are scoped to your workspace directory.

Use subagents liberally — e.g., spawn one to trace a failing matrix through the Fortran
code while you read the relevant paper. Spawn one to compile and test while you investigate.

## BASELINE IMPLEMENTATIONS (compile and compare against these)

Baseline headers are pre-copied into `{agent_dir}/baselines/`. Just compile and run:
```bash
# DBDSQR — gold standard, 379/379 pass, O(n³), score=5 (hard gate)
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/dbdsqr -o {agent_dir}/baselines/dbdsqr/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm && {agent_dir}/baselines/dbdsqr/evaluate "{stcoll_dir}" 200

# HGBSVD — coupling-based, 174/379 pass, O(n²), score=5 (hard gate)
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/hgbsvd -o {agent_dir}/baselines/hgbsvd/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm && {agent_dir}/baselines/hgbsvd/evaluate "{stcoll_dir}" 200

# TGK+STEMR — MR³ on TGK, 93/379 pass, O(n²), score=5 (hard gate)
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/stemr -o {agent_dir}/baselines/stemr/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm && {agent_dir}/baselines/stemr/evaluate "{stcoll_dir}" 200

# TGK+STEXR — XMR on TGK (vanilla=61/379, with GK-form fix=~138) — THE BASELINE YOU'RE IMPROVING
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/stexr -o {agent_dir}/baselines/stexr/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm && {agent_dir}/baselines/stexr/evaluate "{stcoll_dir}" 200
```
Your working file at `{agent_dir}/program/bidiag_svd.h` is never touched by baseline tests.

**Top evolved variants** (if any exist from prior iterations) are also pre-copied into
`{agent_dir}/baselines/top_variants/`. Each has its own subdir with `bidiag_svd.h` and
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
- You can create scratch files in `{agent_dir}/` for notes, intermediate results, etc.
- You can ONLY write/edit files inside `{agent_dir}/`. All other files are read-only.

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

{changelog}

---

## NOW: GET INTO PLAN MODE, THEN DEBUG SYSTEMATICALLY

Do NOT jump straight to implementation. You must work in three phases:

**PHASE 1 — PLAN MODE: THEORY & UNDERSTANDING (mandatory, ~25% of your time)**
Enter plan mode. Read the papers, reference code, and baselines.
Write `{agent_dir}/THEORY_NOTES.md` summarizing:
- Algorithm 4.1 step by step
- The five requirements and what each guarantees
- NCD condition and GK structure
- Where the implementation could diverge from the theory

Do NOT write any C++ code during Phase 1. Only read, analyze, and write notes.

**PHASE 2 — DIAGNOSIS (~35%)**
Run failing tests. Trace through code. Create synthetic test matrices.
Write `{agent_dir}/BUG_REPORT.md` documenting each divergence from theory.

**PHASE 3 — FIX & VERIFY (~40%)**
Fix bugs. Compile. Test. Compare against DSTEXR baseline.
Write `{agent_dir}/RESEARCH.md` with final results.

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
"""


# ---------------------------------------------------------------------------
# Markdown section extraction
# ---------------------------------------------------------------------------

def _extract_section(markdown: str, section_name: str) -> str:
    """Extract the content of a ## section from markdown text."""
    lines = markdown.split("\n")
    capture = False
    result = []
    for line in lines:
        if line.startswith("## ") and section_name.lower() in line.lower():
            capture = True
            continue
        elif line.startswith("## ") and capture:
            break
        elif line.startswith("# ") and capture:
            break
        elif capture:
            result.append(line)
    return "\n".join(result).strip()


# ---------------------------------------------------------------------------
# Code extraction (fallback only — primary capture is file-based)
# ---------------------------------------------------------------------------

def _extract_code_from_response(response: str, language: str = "cpp") -> str:
    """
    Extract evolved code from agent's conversational response.
    Used as fallback when the agent didn't write the output file.

    Strategy: find last ```cpp block with #pragma once, else largest block.
    """
    if not response or not response.strip():
        return response

    # ```cpp blocks
    pattern = r"```" + re.escape(language) + r"\s*\n(.*?)```"
    matches = re.findall(pattern, response, re.DOTALL)
    if matches:
        full_files = [m for m in matches if "#pragma once" in m and len(m.strip()) > 500]
        if full_files:
            best = full_files[-1]
        else:
            best = max(matches, key=len)
        if len(best.strip()) > 100:
            return f"```{language}\n{best.strip()}\n```"

    # Generic ``` blocks
    pattern = r"```\s*\n(.*?)```"
    matches = re.findall(pattern, response, re.DOTALL)
    if matches:
        full_files = [m for m in matches if "#pragma once" in m and len(m.strip()) > 500]
        if full_files:
            best = full_files[-1]
        else:
            best = max(matches, key=len)
        if len(best.strip()) > 100:
            return f"```{language}\n{best.strip()}\n```"

    logger.warning(
        "No substantial code block found in agent response "
        f"({len(response)} chars). Returning raw response."
    )
    return response


# ---------------------------------------------------------------------------
# AsyncIterable prompt wrapper (required by can_use_tool callback)
# ---------------------------------------------------------------------------

async def _make_prompt_iterable(prompt: str):
    """Wrap a string prompt into an AsyncIterable for can_use_tool support."""
    yield {
        "type": "user",
        "session_id": "",
        "message": {"role": "user", "content": prompt},
        "parent_tool_use_id": None,
    }


# ---------------------------------------------------------------------------
# Agent runner
# ---------------------------------------------------------------------------

def _log_message(log_file, message):
    """Write a human-readable log entry for an SDK message."""
    import datetime
    ts = datetime.datetime.now().strftime("%H:%M:%S")

    try:
        # AssistantMessage — agent's thoughts and tool calls
        if hasattr(message, "content") and hasattr(message, "model"):
            for block in message.content:
                if hasattr(block, "text") and block.text.strip():
                    text = block.text.strip()
                    if len(text) > 5000:
                        text = text[:5000] + f"\n... [{len(block.text)} chars total]"
                    log_file.write(f"[{ts}] 🤖 ASSISTANT: {text}\n")
                elif hasattr(block, "name"):
                    # Tool use
                    inp = str(block.input)
                    if len(inp) > 2000:
                        inp = inp[:2000] + "..."
                    log_file.write(f"[{ts}] 🔧 TOOL CALL: {block.name}({inp})\n")
                elif hasattr(block, "thinking") and block.thinking.strip():
                    text = block.thinking.strip()
                    if len(text) > 5000:
                        text = text[:5000] + f"\n... [{len(block.thinking)} chars total]"
                    log_file.write(f"[{ts}] 💭 THINKING: {text}\n")

        # ToolResultBlock in UserMessage
        elif hasattr(message, "tool_use_result") and message.tool_use_result:
            content = str(message.tool_use_result)
            if len(content) > 3000:
                content = content[:3000] + f"\n... [{len(str(message.tool_use_result))} chars total]"
            log_file.write(f"[{ts}] 📋 TOOL RESULT: {content}\n")

        # ResultMessage — session end
        elif hasattr(message, "result") and hasattr(message, "num_turns"):
            log_file.write(
                f"[{ts}] ✅ SESSION END: turns={message.num_turns}, "
                f"cost=${message.total_cost_usd or 0:.4f}, "
                f"stop={message.stop_reason}\n"
            )

        # SystemMessage — task progress
        elif hasattr(message, "subtype"):
            if hasattr(message, "last_tool_name"):
                log_file.write(
                    f"[{ts}] 📊 PROGRESS: tool={message.last_tool_name}, "
                    f"tokens={message.usage.get('total_tokens', '?')}\n"
                )

        log_file.flush()
    except Exception:
        pass  # Never crash on logging


async def _run_agent_session(
    prompt: str,
    options,
    label: str = "agent",
    log_path: str = None,
) -> str:
    """Run a single agent session and collect the text result."""
    sdk = _ensure_sdk()
    result_text = ""
    cost_usd = 0.0

    # can_use_tool requires AsyncIterable prompt
    if options.can_use_tool:
        actual_prompt = _make_prompt_iterable(prompt)
    else:
        actual_prompt = prompt

    # Open session log file if path provided
    log_file = None
    if log_path:
        try:
            log_file = open(log_path, "w")
            log_file.write(f"=== {label} session started ===\n")
            log_file.flush()
        except Exception:
            log_file = None

    try:
        import asyncio
        # 2-hour timeout to prevent infinite hangs after max_turns
        async def _consume():
            nonlocal result_text, cost_usd
            async for message in sdk.query(prompt=actual_prompt, options=options):
                if log_file:
                    _log_message(log_file, message)
                if hasattr(message, "result") and message.result:
                    result_text = message.result
                if hasattr(message, "total_cost_usd") and message.total_cost_usd:
                    cost_usd = message.total_cost_usd
                if hasattr(message, "content") and not result_text:
                    for block in getattr(message, "content", []):
                        if hasattr(block, "text"):
                            result_text += block.text

        await asyncio.wait_for(_consume(), timeout=7200)

    except asyncio.TimeoutError:
        logger.warning(f"{label} session timed out after 2 hours")
        if log_file:
            log_file.write(f"⏰ SESSION TIMEOUT: 2 hour limit reached\n")
    except Exception as e:
        # Don't re-raise — the agent may have written useful output to disk
        # before crashing (e.g., max_turns reached, SDK exit code -9, etc.).
        # The caller will read program/bidiag_svd.h from disk regardless.
        logger.error(f"{label} session failed: {e}")
        if log_file:
            log_file.write(f"❌ SESSION ERROR: {e}\n")
        logger.info(f"{label} session error is non-fatal — will try to read output from disk")
    finally:
        if log_file:
            log_file.close()

    if cost_usd > 0:
        logger.info(f"{label} session cost: ${cost_usd:.4f}")

    return result_text


# ---------------------------------------------------------------------------
# Main adapter class
# ---------------------------------------------------------------------------

class ClaudeAgentLLM(LLMInterface):
    """
    Single-phase agentic LLM for OpenEvolve with per-agent isolated workspaces.

    Each generate call creates scratch/agent_NNN/ with:
      program/bidiag_svd.h  — the file to evolve (pre-copied, agent edits it)
      evaluate.cpp          — pre-copied for compilation

    Write/Edit restricted to the agent's directory via can_use_tool callback.
    Result captured by reading program/bidiag_svd.h from disk (not from text output).
    """

    def __init__(
        self,
        model: str = "claude-opus-4-6",
        cwd: Optional[str] = None,
        scratch_dir: Optional[str] = None,
        variants_dir: Optional[str] = None,
        max_budget_usd: Optional[float] = None,
        effort: Optional[str] = "high",
        enable_tools: bool = True,
        enable_checkpointing: bool = False,
        load_project_settings: bool = True,
        permission_mode: str = "bypassPermissions",
        custom_agents: Optional[dict] = None,
        language: str = "cpp",
        prompt_template: str = "xmr_debug",
        initial_program: Optional[str] = None,
    ):
        # Unset CLAUDECODE env var to allow nested Claude Code sessions
        # (we may be running inside a Claude Code session)
        os.environ.pop("CLAUDECODE", None)

        _ensure_sdk()
        self.model = model
        self.prompt_template = prompt_template
        self.cwd = cwd
        self.max_budget_usd = max_budget_usd
        self.effort = effort
        self.enable_tools = enable_tools
        self.enable_checkpointing = enable_checkpointing
        self.load_project_settings = load_project_settings
        self.permission_mode = permission_mode
        self.custom_agents = custom_agents or {}
        self.language = language

        # Initial program to seed the agent workspace with
        # If set, this file is copied as program/bidiag_svd.h instead of src/bidiag_svd.h
        self._initial_program = initial_program

        # Agent counter (PID-prefixed to avoid collisions across worker processes)
        self._agent_counter = 0
        self._pid = os.getpid()

        # Scratch root — use explicit scratch_dir if given, else default
        if scratch_dir:
            self._scratch_dir = os.path.abspath(scratch_dir)
        elif self.cwd:
            self._scratch_dir = os.path.abspath(os.path.join(self.cwd, "scratch"))
        else:
            self._scratch_dir = None
        if self._scratch_dir:
            os.makedirs(self._scratch_dir, exist_ok=True)

        # Variants directory (research reports, changelog, saved variants)
        if variants_dir:
            self._variants_dir = os.path.abspath(variants_dir)
        elif self.cwd:
            self._variants_dir = os.path.join(self.cwd, "evolved_variants")
        else:
            self._variants_dir = None

        # Current agent's directory (set per-call for the permission callback)
        self._current_agent_dir: Optional[str] = None

        budget_str = f"${max_budget_usd}" if max_budget_usd else "unlimited"
        logger.info(
            f"ClaudeAgentLLM initialized: model={model}, "
            f"effort={effort}, tools={enable_tools}, budget={budget_str}, "
            f"scratch={self._scratch_dir}"
        )

    # -------------------------------------------------------------------
    # Agent workspace setup
    # -------------------------------------------------------------------

    def _create_agent_workspace(self) -> str:
        """
        Create an isolated workspace for this agent call.

        Returns the absolute path to scratch/agent_NNN/.
        Pre-populates with:
          program/bidiag_svd.h  — copied from src/
          evaluate.cpp          — copied from src/
          baselines/            — pre-built baseline dirs for comparison
        """
        self._agent_counter += 1
        agent_name = f"agent_{self._pid}_{self._agent_counter:04d}"
        agent_dir = os.path.join(self._scratch_dir, agent_name)
        program_dir = os.path.join(agent_dir, "program")
        os.makedirs(program_dir, exist_ok=True)

        # Copy initial program into workspace as bidiag_svd.h
        src_dir = os.path.join(self.cwd, "src")
        if self._initial_program and os.path.exists(self._initial_program):
            shutil.copy2(self._initial_program, os.path.join(program_dir, "bidiag_svd.h"))
        else:
            shutil.copy2(
                os.path.join(src_dir, "bidiag_svd.h"),
                os.path.join(program_dir, "bidiag_svd.h"),
            )
        shutil.copy2(
            os.path.join(src_dir, "evaluate.cpp"),
            os.path.join(agent_dir, "evaluate.cpp"),
        )

        # Pre-populate baseline dirs for comparison testing
        baselines = {
            "dbdsqr": ["bidiag_dbdsqr.h"],
            "hgbsvd": ["bidiag_hgbsvd.h"],
            "stemr": ["bidiag_tgk_stemr.h", "bidiag_tgk_common.h"],
            "stexr": ["bidiag_tgk_stexr.h", "bidiag_tgk_common.h"],
        }
        for name, headers in baselines.items():
            bdir = os.path.join(agent_dir, "baselines", name)
            os.makedirs(bdir, exist_ok=True)
            for hdr in headers:
                src_file = os.path.join(src_dir, hdr)
                if os.path.exists(src_file):
                    dst_name = "bidiag_svd.h" if hdr.startswith("bidiag_") and hdr != "bidiag_tgk_common.h" else hdr
                    shutil.copy2(src_file, os.path.join(bdir, dst_name))

        # Copy top variants from THIS experiment's variants directory
        self._copy_top_variants(agent_dir)

        logger.info(f"Created agent workspace: {agent_dir}")
        return agent_dir

    def _copy_top_variants(self, agent_dir: str, max_variants: int = 5) -> list:
        """
        Copy the top N evolved variants into the agent workspace for comparison.

        Parses variant filenames (v{N}_{score}.h), sorts by score descending,
        copies the top max_variants into baselines/top_variants/.
        Returns list of (name, score) tuples copied.
        """
        if not self._variants_dir or not os.path.exists(self._variants_dir):
            return []

        import glob as globmod

        # Find all variant files
        variant_files = globmod.glob(
            os.path.join(self._variants_dir, "v*_*.h")
        )
        if not variant_files:
            return []

        # Parse scores from filenames: v{N}_{score}.h
        scored = []
        for path in variant_files:
            base = os.path.basename(path)  # e.g. "v3_75.h"
            try:
                parts = base[:-2].split("_", 1)  # ["v3", "75"]
                score = float(parts[1])
                scored.append((score, base, path))
            except (ValueError, IndexError):
                pass

        if not scored:
            return []

        # Sort by score descending, take top N
        scored.sort(key=lambda x: x[0], reverse=True)
        top = scored[:max_variants]

        # Copy into workspace
        top_dir = os.path.join(agent_dir, "baselines", "top_variants")
        os.makedirs(top_dir, exist_ok=True)
        result = []
        for score, name, path in top:
            # Each variant gets its own subdir so it can be compiled independently
            vdir = os.path.join(top_dir, name[:-2])  # strip .h
            os.makedirs(vdir, exist_ok=True)
            shutil.copy2(path, os.path.join(vdir, "bidiag_svd.h"))
            result.append((name, score))

            # Also copy its RESEARCH.md if available
            variant_num = name.split("_")[0]  # "v3"
            research_pattern = os.path.join(
                self._variants_dir, f"{variant_num}_RESEARCH.md"
            )
            research_files = globmod.glob(research_pattern)
            if research_files:
                shutil.copy2(
                    research_files[0],
                    os.path.join(vdir, "RESEARCH.md"),
                )

        if result:
            logger.info(
                f"Copied {len(result)} top variants to workspace: "
                + ", ".join(f"{n}({s:.0f})" for n, s in result)
            )

        return result

    def _read_agent_output(self, agent_dir: str) -> Optional[str]:
        """
        Read the agent's output file from disk.

        Returns the file content wrapped in ```cpp fences for OpenEvolve's
        parse_full_rewrite, or None if the file wasn't written.
        """
        output_file = os.path.join(agent_dir, "program", "bidiag_svd.h")
        if not os.path.exists(output_file):
            return None

        try:
            with open(output_file) as f:
                content = f.read()
        except Exception as e:
            logger.warning(f"Failed to read agent output: {e}")
            return None

        if not content.strip() or len(content.strip()) < 100:
            return None

        # Verify it looks like a valid evolved file
        if "#pragma once" not in content:
            logger.warning("Agent output missing #pragma once — may be incomplete")

        # Check for fixed XMR library — if the agent rebuilt the Fortran,
        # replace lib/libxmr_c.a so the evaluator uses the fixed version
        for lib_name in ["libxmr_fixed.a", "libxmr_debug.a"]:
            fixed_lib = os.path.join(agent_dir, "xmr_src", lib_name)
            if os.path.exists(fixed_lib):
                lib_dest = os.path.join(self.cwd, "lib", "libxmr_c.a")
                shutil.copy2(fixed_lib, lib_dest)
                logger.info(
                    f"Replaced lib/libxmr_c.a with agent's {lib_name} "
                    f"({os.path.getsize(fixed_lib)} bytes)"
                )
                break

        return f"```{self.language}\n{content}\n```"

    # -------------------------------------------------------------------
    # Research report staging
    # -------------------------------------------------------------------

    def _stage_research(self, agent_dir: str) -> None:
        """
        Read the agent's RESEARCH.md and stage it for the evaluator.

        The evaluator's save_variant_and_log() picks up the staged file,
        saves it as v{N}_RESEARCH.md, and extracts the abstract for the
        changelog one-liner.
        """
        research_file = os.path.join(agent_dir, "RESEARCH.md")
        if not os.path.exists(research_file):
            logger.debug("No RESEARCH.md found in agent workspace")
            return

        try:
            with open(research_file) as f:
                research = f.read().strip()
        except Exception:
            return

        if not research:
            return

        # Stage for evaluator to pick up
        os.makedirs(self._variants_dir, exist_ok=True)
        pending = os.path.join(self._variants_dir, ".pending_research.md")
        try:
            with open(pending, "w") as f:
                f.write(research)
            logger.info(f"Staged RESEARCH.md ({len(research)} chars)")
        except Exception as e:
            logger.warning(f"Failed to stage RESEARCH.md: {e}")

    # -------------------------------------------------------------------
    # Prior research summaries (for inclusion in prompt)
    # -------------------------------------------------------------------

    def _build_prior_research_summary(self) -> str:
        """
        Build a summary of prior agents' research for inclusion in the prompt.

        Reads the abstract from the last N RESEARCH.md files and lists all
        available research files so the agent can read any in full.
        """
        if not self._variants_dir:
            return "No prior research available."

        if not os.path.exists(self._variants_dir):
            return "No prior research available (first iteration)."

        # Find all RESEARCH.md files, sorted by variant number
        import glob as globmod
        research_files = globmod.glob(os.path.join(self._variants_dir, "v*_RESEARCH.md"))
        if not research_files:
            return "No prior research available yet."

        # Parse variant numbers and sort
        numbered = []
        for path in research_files:
            base = os.path.basename(path)
            try:
                num = int(base.split("_")[0][1:])
                numbered.append((num, path))
            except (ValueError, IndexError):
                pass

        numbered.sort(key=lambda x: x[0])

        # Extract abstracts from the last 5 research papers
        max_recent = 5
        recent = numbered[-max_recent:]
        parts = []

        parts.append(
            f"{len(numbered)} prior research reports available at `{self._variants_dir}/`.\n"
            f"You can Read any of them in full for detailed analysis.\n"
        )

        # List all available research files
        if len(numbered) > max_recent:
            parts.append("**All reports:** " + ", ".join(
                f"v{num}" for num, _ in numbered
            ) + "\n")

        # Show abstracts of recent ones
        parts.append("**Recent abstracts:**\n")
        for num, path in recent:
            try:
                with open(path) as f:
                    content = f.read()
            except Exception:
                continue

            # Extract title (first # line)
            title = f"v{num}"
            for line in content.split("\n"):
                if line.startswith("# "):
                    title = line[2:].strip()
                    break

            # Extract abstract section
            abstract = _extract_section(content, "Abstract")
            if not abstract:
                abstract = "(no abstract)"
            # Truncate long abstracts
            if len(abstract) > 300:
                abstract = abstract[:300] + "..."

            parts.append(f"- **{title}**: {abstract}")

        return "\n".join(parts)

    # -------------------------------------------------------------------
    # Write permission callback (scoped to current agent directory)
    # -------------------------------------------------------------------

    async def _write_permission_callback(
        self,
        tool_name: str,
        tool_input: dict,
        context,
    ):
        """
        Directory-scoped write permission callback.

        - Write/Edit: only allowed inside the current agent's directory
        - ExitPlanMode: auto-approved (workaround for bug #30463 where
          ExitPlanMode always prompts even in bypassPermissions mode)
        - All other tools: allowed unconditionally
        - AskUserQuestion: blocked via disallowed_tools
        """
        sdk = _ensure_sdk()

        # Auto-approve ExitPlanMode — workaround for SDK bug #30463
        # where it always prompts for approval even in bypassPermissions mode
        if tool_name == "ExitPlanMode":
            return sdk.PermissionResultAllow()

        # Only gate Write and Edit tools
        if tool_name not in ("Write", "Edit"):
            return sdk.PermissionResultAllow()

        # Get the file path from tool input
        file_path = tool_input.get("file_path", "")
        if not file_path:
            return sdk.PermissionResultDeny(
                message="Write/Edit requires a file_path"
            )

        # Resolve to absolute path
        abs_path = os.path.abspath(file_path)

        # Check: is it inside the current agent's directory?
        if self._current_agent_dir:
            agent_prefix = self._current_agent_dir + os.sep
            if abs_path.startswith(agent_prefix) or abs_path == self._current_agent_dir:
                return sdk.PermissionResultAllow()

        return sdk.PermissionResultDeny(
            message=(
                f"Write/Edit denied for: {file_path}\n"
                f"You can only write to your workspace: {self._current_agent_dir}/\n"
                f"Edit {self._current_agent_dir}/program/bidiag_svd.h for code changes.\n"
                f"Use {self._current_agent_dir}/ for notes and scratch files."
            )
        )

    # -------------------------------------------------------------------
    # Build agent options
    # -------------------------------------------------------------------

    def _build_options(self):
        """Build ClaudeAgentOptions for an agent session."""
        sdk = _ensure_sdk()

        # Full tool access — no whitelist, no blacklist. All Claude Code tools
        # available (Read, Write, Edit, Bash, Grep, Glob, WebSearch, WebFetch,
        # Agent, NotebookEdit, planning, exploration, etc.).
        # Write restrictions enforced solely via can_use_tool callback.

        # Disable project settings (CLAUDE.md) — everything the agent needs
        # is in the user prompt template.
        setting_sources = []

        # Use Claude Code's built-in system prompt (tool instructions, safety,
        # coding guidelines). Without this preset, the SDK uses a minimal
        # system prompt that lacks Claude Code's full behavior.
        system_prompt_config = {"type": "preset", "preset": "claude_code"}

        # Convert custom agents to AgentDefinition if needed
        agents = None
        if self.custom_agents:
            agents = {}
            for name, agent_def in self.custom_agents.items():
                if isinstance(agent_def, dict):
                    agents[name] = sdk.AgentDefinition(**agent_def)
                else:
                    agents[name] = agent_def

        # Build options — all tools available except AskUserQuestion
        # (agent must be fully autonomous, no human-in-the-loop).
        # Write restrictions enforced via can_use_tool callback.
        # PermissionRequest hook auto-approves ExitPlanMode (bypasses CLI's
        # interactive prompt which fails headlessly with "Stream closed").
        async def _permission_request_hook(input_data, tool_use_id, context):
            tool_name = input_data.get("tool_name", "")
            if tool_name == "ExitPlanMode":
                return {
                    "continue_": True,
                    "hookSpecificOutput": {
                        "hookEventName": "PermissionRequest",
                        "decision": {"allow": True},
                    },
                }
            return {"continue_": True}

        opts_kwargs = dict(
            model=self.model,
            permission_mode=self.permission_mode,
            enable_file_checkpointing=self.enable_checkpointing,
            can_use_tool=self._write_permission_callback,
            system_prompt=system_prompt_config,
            disallowed_tools=["AskUserQuestion"],
            max_turns=200,  # prevent infinite loops
            hooks={
                "PermissionRequest": [
                    sdk.HookMatcher(
                        matcher="ExitPlanMode",
                        hooks=[_permission_request_hook],
                    ),
                ],
            },
        )

        if not self.enable_tools:
            opts_kwargs["allowed_tools"] = []

        if self.max_budget_usd is not None:
            opts_kwargs["max_budget_usd"] = self.max_budget_usd

        if self.effort:
            opts_kwargs["effort"] = self.effort

        if self.cwd:
            opts_kwargs["cwd"] = self.cwd

        if setting_sources:
            opts_kwargs["setting_sources"] = setting_sources

        if agents:
            opts_kwargs["agents"] = agents

        return sdk.ClaudeAgentOptions(**opts_kwargs)

    # -------------------------------------------------------------------
    # Evolution changelog (read from evaluator's persistent log)
    # -------------------------------------------------------------------

    def _read_changelog(self) -> str:
        """
        Read the evolution changelog maintained by evaluator.py.
        Returns a concise summary for inclusion in prompts, or empty string.
        """
        if not self._variants_dir:
            return ""

        changelog_path = os.path.join(self._variants_dir, "CHANGELOG.md")
        if not os.path.exists(changelog_path):
            return ""

        try:
            with open(changelog_path) as f:
                content = f.read()
        except Exception:
            return ""

        if not content.strip():
            return ""

        # Count variants
        lines = [
            l for l in content.split("\n")
            if l.startswith("| ") and not l.startswith("| #") and not l.startswith("|---")
        ]
        num_variants = len(lines)
        if num_variants == 0:
            return ""

        # Truncate if too many entries — show header + last 20
        all_lines = content.split("\n")
        if num_variants > 20:
            header = all_lines[:4]
            data_lines = [l for l in all_lines[4:] if l.strip()]
            truncated = (
                header
                + [f"| ... | ({num_variants - 20} earlier variants omitted) | ... | ... | ... | ... | ... | ... | ... |"]
                + data_lines[-20:]
            )
            content = "\n".join(truncated)

        return (
            f"## Evolution History ({num_variants} variants so far)\n\n"
            f"Previous variants saved at: `{self._variants_dir}/`\n"
            f"- Code: `{self._variants_dir}/v{{N}}_{{score}}.h` — Read any variant to see its full code\n"
            f"- Research: `{self._variants_dir}/v{{N}}_RESEARCH.md` — Structured research report\n\n"
            f"{content}\n"
        )

    # -------------------------------------------------------------------
    # Generation (single-phase)
    # -------------------------------------------------------------------

    async def generate(self, prompt: str, **kwargs) -> str:
        """Generate text from a prompt."""
        return await self.generate_with_context(
            system_message=kwargs.get("system_message", ""),
            messages=[{"role": "user", "content": prompt}],
            **kwargs,
        )

    async def generate_with_context(
        self, system_message: str, messages: List[Dict[str, str]], **kwargs
    ) -> str:
        """
        Single-phase generation with isolated agent workspace.

        1. Create scratch/agent_NNN/ with pre-copied files
        2. Run agent with write permissions scoped to that directory
        3. Read program/bidiag_svd.h from disk as the result
        4. Fallback: extract code from agent's text response
        """
        # Build the combined prompt from OpenEvolve messages
        prompt_parts = []
        for msg in messages:
            role = msg.get("role", "user")
            content = msg.get("content", "")
            if role == "user":
                prompt_parts.append(content)
            elif role == "assistant":
                prompt_parts.append(f"[Previous assistant response]\n{content}")

        # OpenEvolve passes the current code + "rewrite to improve" but we
        # don't use it — our prompt templates provide all context. The agent
        # reads the code from its workspace (program/bidiag_svd.h).
        # Inter-iteration info comes from changelog + prior_research + top_variants.

        # Read evolution changelog and prior research
        changelog = self._read_changelog()
        prior_research = self._build_prior_research_summary()

        # Create isolated workspace for this agent
        agent_dir = self._create_agent_workspace()
        self._current_agent_dir = agent_dir

        try:
            # Build unified prompt with workspace paths
            # STCollection path for evaluate binary
            stcoll_dir = os.path.join(
                self.cwd, "..", "bidiag-algo",
                "MRRR Resources", "STCollection", "DATA"
            )
            # Select prompt template
            if self.prompt_template == "exploration":
                template = EXPLORATION_PROMPT_TEMPLATE
            else:
                template = XMR_DEBUG_PROMPT_TEMPLATE
            prompt = template.format(
                agent_dir=agent_dir,
                stcoll_dir=stcoll_dir,
                prior_research=prior_research,
                changelog=changelog,
            )

            options = self._build_options()

            log_path = os.path.join(agent_dir, "session.log")
            logger.info(f"Running evolution agent in {agent_dir}")
            logger.info(f"Live session log: {log_path}")
            print(f"\n{'='*60}")
            print(f"  Agent workspace: {agent_dir}")
            print(f"  Session log:     {log_path}")
            print(f"  Monitor with:    tail -f {log_path}")
            print(f"{'='*60}\n", flush=True)

            result_text = await _run_agent_session(
                prompt=prompt,
                options=options,
                label="Evolution",
                log_path=log_path,
            )

            # Stage the agent's RESEARCH.md for the evaluator
            self._stage_research(agent_dir)

            # Primary: read the output file from disk
            file_result = self._read_agent_output(agent_dir)
            if file_result:
                logger.info(
                    f"Read agent output from {agent_dir}/program/bidiag_svd.h "
                    f"({len(file_result)} chars)"
                )
                return file_result

            # Fallback: extract code from conversational response
            logger.warning(
                "Agent output file not found or empty — "
                "falling back to text extraction"
            )
            if not result_text:
                logger.warning("Agent returned empty result")
                return ""

            processed = _extract_code_from_response(result_text, self.language)
            logger.debug(
                f"Text extraction: {len(result_text)} chars raw -> "
                f"{len(processed)} chars processed"
            )
            return processed

        finally:
            self._current_agent_dir = None


# ---------------------------------------------------------------------------
# Factory function for OpenEvolve's init_client hook
# ---------------------------------------------------------------------------

def create_claude_agent_client(model_cfg) -> ClaudeAgentLLM:
    """
    Factory function compatible with OpenEvolve's LLMModelConfig.init_client.

    Custom attributes (cwd, max_budget_usd, etc.) are NOT preserved by
    OpenEvolve's config serialization (dataclasses.asdict drops non-field
    attributes). We read them from a persisted JSON file instead.
    """
    # Try model_cfg first (works in parent process), fall back to settings file
    settings = _load_agent_settings()

    def _get(key, default=None):
        """Get from model_cfg attr, then settings file, then default."""
        val = getattr(model_cfg, key, None)
        if val is not None:
            return val
        return settings.get(key, default)

    return ClaudeAgentLLM(
        model=getattr(model_cfg, "name", None) or settings.get("model", "claude-opus-4-6"),
        cwd=_get("cwd"),
        scratch_dir=_get("scratch_dir"),
        variants_dir=_get("variants_dir"),
        max_budget_usd=_get("max_budget_usd"),
        effort=_get("effort", "high"),
        enable_tools=_get("enable_tools", True),
        enable_checkpointing=_get("enable_checkpointing", False),
        load_project_settings=_get("load_project_settings", True),
        permission_mode=_get("permission_mode", "bypassPermissions"),
        custom_agents=_get("custom_agents"),
        language=_get("language", "cpp"),
        prompt_template=_get("prompt_template", "xmr_debug"),
        initial_program=_get("initial_program"),
    )
