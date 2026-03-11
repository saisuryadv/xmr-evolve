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

UNIFIED_PROMPT_TEMPLATE = """\
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
{agent_dir}/evaluate
```

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
# DBDSQR — gold standard, 270/270 pass, O(n³), score=5 (hard gate)
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/dbdsqr -o {agent_dir}/baselines/dbdsqr/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm && {agent_dir}/baselines/dbdsqr/evaluate

# HGBSVD — coupling-based, 135/270 pass, O(n²), score=68.6
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/hgbsvd -o {agent_dir}/baselines/hgbsvd/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm && {agent_dir}/baselines/hgbsvd/evaluate

# TGK+STEMR — MR³ on TGK, 65/270 pass, O(n²), score=51.2
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/stemr -o {agent_dir}/baselines/stemr/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm && {agent_dir}/baselines/stemr/evaluate

# TGK+STEXR — XMR on TGK, 43/270 pass, O(n²), score=41.0
g++ -std=c++17 -O2 -Isrc/clapack -I{agent_dir}/baselines/stexr -o {agent_dir}/baselines/stexr/evaluate {agent_dir}/evaluate.cpp -Llib -lxmr_c -lm && {agent_dir}/baselines/stexr/evaluate
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

### STEP 0: READ THE PAPERS (MANDATORY — DO THIS BEFORE ANYTHING ELSE)

This is not optional. You MUST read these papers before writing any code, running any
tests, or even looking at the source code. The papers contain the algorithms, theorems,
and mathematical foundations you need. Without reading them, you will waste time
reinventing approaches that have already been tried and failed.

Read these 3 papers now using the Read tool (they are PDFs in the knowledge/ directory):

```
knowledge/paper_willems_lang_2012_mr3gk.pdf
knowledge/paper_grosser_lang_2001_on2.pdf
knowledge/paper_willems_lang_2013_framework.pdf
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

## WHAT TO READ (in priority order)

### 1. PAPERS (PDFs) — READ THESE FIRST, they contain the actual math
The knowledge/ files are only summaries. The real algorithms, theorems, and proofs are in the papers.
Use the Read tool on these PDFs (you can read PDFs directly):

**Must-read papers** (read before implementing anything):
- `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Papers/2012 - The MR3-GK Algorithm for the Bidiagonal SVD - Willems and Lang.pdf`
  → Algorithm 4.1 (the target algorithm), NCD condition, GK structure, Theorem 4.5
- `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Papers/2001 - An O(n^2) Algorithm for the Bidiagonal SVD - Großer and Lang.pdf`
  → Coupling B^TB/BB^T, how HGBSVD works, why it fails at deep recursion
- `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Papers/2013 - A Framework for the MR3 Algorithm Theory and Implementation - Willems and Lang.pdf`
  → Five requirements for MR³, block factorizations, DSTEMR underflow bugs

**Read when investigating specific failures:**
- `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Papers/2020 - Bidiagonal SVD Computation via an Associated Tridiagonal Eigenproblem - Marques Demmel and Vasconcelos.pdf`
  → DBDSVDX bugs, CHKBD matrix analysis, modern failure catalog
- `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Papers/1997 - Dhillon Thesis.pdf`
  → Original MR³: representation tree, twisted factorizations, assumptions that fail
- `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Papers/2005 - Glued Matrices and the MRRR Algorithm - Dhillon Parlett and Vömel.pdf`
  → Wilkinson/glued matrix failures, random perturbation fixes
- `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Papers/2000 - Relatively Robust Representations of Symmetric Tridiagonals - Parlett and Dhillon.pdf`
  → RRR theory, when LDL^T is/isn't an RRR
- `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Papers/1990 - Accurate Singular Values of Bidiagonal Matrices - Demmel and Kahan.pdf`
  → Foundational accuracy theory, splitting criteria

### 2. Knowledge base (summaries of the above papers + project-specific info)
- `knowledge/INDEX.md` — Master reference: all known bugs, test matrix formulas, thresholds
- `knowledge/PRIOR_APPROACHES.md` — 12 approaches tried, what worked and failed
- `knowledge/BASELINES.md` — 5 baseline algorithms compared

### 3. Reference Fortran code (production implementations of the algorithms in the papers)
- XMR: `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Code/tdsolver/xmr/SRC/O/` (45 files, `dstexr.f` master)
- HGBSVD: `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Code/hgbsvd/hgbsvd/v2/` (24 files, `dbdsgr.f` master)

### 4. Baseline C++ headers
- `src/bidiag_dbdsqr.h`, `src/bidiag_hgbsvd.h`, `src/bidiag_tgk_stemr.h`, `src/bidiag_tgk_stexr.h`

### 5. Top evolved variants
- `{agent_dir}/baselines/top_variants/` (read their RESEARCH.md files)

## DO NOT READ (irrelevant to this project)
- **DO NOT** read `approach_*.py` files — these are Python prototypes from a different experiment
- **DO NOT** read files in `../bidiag-algo/` except Papers, Code, and STCollection directories
- **DO NOT** read `evaluate_alphaevolve.py`, `mr3_tridiag.py`, `problem_statement_alphaevolve.md`
- **DO NOT** read `accuracy_tables.txt`, `results_approach_*.txt`
- **DO NOT** read `.claude/` memory files or auto-memory — these are for a different agent
- Focus on **C++ implementation**, **papers**, **Fortran reference code**, and **knowledge base**

{changelog}

---

{openevolve_prompt}
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
                    # Truncate very long text blocks
                    text = block.text.strip()
                    if len(text) > 500:
                        text = text[:500] + f"... [{len(block.text)} chars total]"
                    log_file.write(f"[{ts}] 🤖 ASSISTANT: {text}\n")
                elif hasattr(block, "name"):
                    # Tool use
                    inp = str(block.input)
                    if len(inp) > 300:
                        inp = inp[:300] + "..."
                    log_file.write(f"[{ts}] 🔧 TOOL CALL: {block.name}({inp})\n")
                elif hasattr(block, "thinking") and block.thinking.strip():
                    text = block.thinking.strip()
                    if len(text) > 300:
                        text = text[:300] + "..."
                    log_file.write(f"[{ts}] 💭 THINKING: {text}\n")

        # ToolResultBlock in UserMessage
        elif hasattr(message, "tool_use_result") and message.tool_use_result:
            content = str(message.tool_use_result)
            if len(content) > 200:
                content = content[:200] + "..."
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
        async for message in sdk.query(prompt=actual_prompt, options=options):
            # Log every message for live monitoring
            if log_file:
                _log_message(log_file, message)

            if hasattr(message, "result") and message.result:
                result_text = message.result
            if hasattr(message, "total_cost_usd") and message.total_cost_usd:
                cost_usd = message.total_cost_usd
            # Fallback: capture text blocks from assistant messages
            if hasattr(message, "content") and not result_text:
                for block in getattr(message, "content", []):
                    if hasattr(block, "text"):
                        result_text += block.text

    except Exception as e:
        logger.error(f"{label} session failed: {e}")
        if log_file:
            log_file.write(f"❌ SESSION ERROR: {e}\n")
        raise
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
        extra_system_prompt: Optional[str] = None,
        custom_agents: Optional[dict] = None,
        language: str = "cpp",
    ):
        # Unset CLAUDECODE env var to allow nested Claude Code sessions
        # (we may be running inside a Claude Code session)
        os.environ.pop("CLAUDECODE", None)

        _ensure_sdk()
        self.model = model
        self.cwd = cwd
        self.max_budget_usd = max_budget_usd
        self.effort = effort
        self.enable_tools = enable_tools
        self.enable_checkpointing = enable_checkpointing
        self.load_project_settings = load_project_settings
        self.permission_mode = permission_mode
        self.extra_system_prompt = extra_system_prompt
        self.custom_agents = custom_agents or {}
        self.language = language

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

        # Copy source files into workspace
        src_dir = os.path.join(self.cwd, "src")
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

        # Pre-populate top evolved variants for comparison testing
        # Find the top 5 variants by score from the variants directory
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

        Allows Write/Edit only for files inside the current agent's directory
        (scratch/agent_NNN/). All other tools allowed unconditionally.
        """
        sdk = _ensure_sdk()

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

    def _build_options(
        self,
        system_prompt: Optional[str] = None,
    ):
        """Build ClaudeAgentOptions for an agent session."""
        sdk = _ensure_sdk()

        # Full tool access — no whitelist, no blacklist. All Claude Code tools
        # available (Read, Write, Edit, Bash, Grep, Glob, WebSearch, WebFetch,
        # Agent, NotebookEdit, planning, exploration, etc.).
        # Write restrictions enforced solely via can_use_tool callback.

        # Settings sources
        setting_sources = ["user", "project"] if self.load_project_settings else []

        # System prompt
        full_system_prompt = system_prompt or ""
        if self.extra_system_prompt:
            full_system_prompt += "\n\n" + self.extra_system_prompt

        # Convert custom agents to AgentDefinition if needed
        agents = None
        if self.custom_agents:
            agents = {}
            for name, agent_def in self.custom_agents.items():
                if isinstance(agent_def, dict):
                    agents[name] = sdk.AgentDefinition(**agent_def)
                else:
                    agents[name] = agent_def

        # Build options — no allowed_tools whitelist so all tools are available
        opts_kwargs = dict(
            model=self.model,
            permission_mode=self.permission_mode,
            enable_file_checkpointing=self.enable_checkpointing,
            can_use_tool=self._write_permission_callback,
        )

        if not self.enable_tools:
            opts_kwargs["allowed_tools"] = []

        if self.max_budget_usd is not None:
            opts_kwargs["max_budget_usd"] = self.max_budget_usd

        if self.effort:
            opts_kwargs["effort"] = self.effort

        if full_system_prompt.strip():
            opts_kwargs["system_prompt"] = full_system_prompt

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

        openevolve_prompt = "\n\n".join(prompt_parts)

        # Read evolution changelog and prior research
        changelog = self._read_changelog()
        prior_research = self._build_prior_research_summary()

        # Create isolated workspace for this agent
        agent_dir = self._create_agent_workspace()
        self._current_agent_dir = agent_dir

        try:
            # Build unified prompt with workspace paths
            prompt = UNIFIED_PROMPT_TEMPLATE.format(
                agent_dir=agent_dir,
                prior_research=prior_research,
                changelog=changelog,
                openevolve_prompt=openevolve_prompt,
            )

            options = self._build_options(system_prompt=system_message)

            log_path = os.path.join(agent_dir, "session.log")
            logger.info(f"Running evolution agent in {agent_dir}")
            logger.info(f"Live session log: {log_path}")

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
        extra_system_prompt=_get("extra_system_prompt"),
        custom_agents=_get("custom_agents"),
        language=_get("language", "cpp"),
    )
