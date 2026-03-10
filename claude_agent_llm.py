"""
Claude Agent SDK adapter for OpenEvolve

Drop-in replacement for OpenEvolve's OpenAI LLM backend that uses the
Claude Agent SDK instead.  Each "generation" is a TWO-PHASE agent pipeline:

Phase 1 — RESEARCH AGENT (read-only, cheap model):
  Analyzes the parent program's test failures by reading the knowledge base,
  reference code, and documented bugs.  Produces a focused "research brief"
  with specific clues about what technique to try and why.

Phase 2 — MUTATION AGENT (full tools, main model):
  Receives the OpenEvolve prompt PLUS the research brief.  Has all tools
  (read, write, bash, grep, subagents) to implement the suggested changes.
  Returns the complete evolved file in ```cpp fences.

This two-phase design means the mutation agent starts each session with
domain-specific guidance instead of having to rediscover failure modes
from scratch.
"""

import asyncio
import logging
import os
import re
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
# Phase 1: Research agent prompt
# ---------------------------------------------------------------------------

RESEARCH_SYSTEM_PROMPT = """\
You are a numerical linear algebra research analyst for bidiagonal SVD algorithms.
Your job is to analyze test failures and produce actionable research briefs that
guide a code mutation agent.

RULES:
- You MUST NOT modify any files. You are read-only.
- NEVER ask questions or request clarification. Produce your analysis and stop.
- Be CONCRETE: name specific functions, line numbers, formulas — not vague advice.
- Read AT MOST 3 files from knowledge/ — pick the most relevant ones, don't read all 13.
  Start with knowledge/INDEX.md (master reference), then pick 1-2 based on the failure pattern:
  * Orthogonality failures → knowledge/PRIOR_APPROACHES.md
  * HGBSVD INFO!=0 → knowledge/grosser_lang_2001_hgbsvd.md
  * Specific matrix failures → knowledge/EVALUATION.md (pattern descriptions)
  * Scaling problems → knowledge/BASELINES.md
"""

RESEARCH_PROMPT_TEMPLATE = """\
## Task: Failure Analysis & Research Brief

Below is the current evolved program with its evaluation results.
Produce a SHORT, ACTIONABLE research brief (max 500 words).

{openevolve_prompt}

---

## Required Analysis (be concrete, not vague):

1. **FAILING TESTS** — List the specific pattern names that fail (e.g., "stemr_killer at n=200:
   orthoU=847.3"). Group by failure mode: INFO!=0, residual blow-up, orthoU/V blow-up.

2. **ROOT CAUSE** — What specific code path triggers the failure? Name the function and the
   problematic logic. Examples of good root causes:
   - "The TGK extraction in extract_uv() doesn't normalize when eigenvector has near-zero
     odd-row norm, causing orthoU > 1000"
   - "dbdsgr_ returns INFO=5 for stemr_killer because B^TB has condition > 10^20"
   Bad root causes (too vague, DO NOT write these):
   - "orthogonality is lost due to numerical issues"
   - "the algorithm doesn't handle clustered eigenvalues well"

3. **SPECIFIC FIX** — One concrete code change. Must include:
   - Which function to modify
   - What condition to add/change (with formula)
   - Why this fixes the identified root cause
   - O(n²) verification: state whether the fix is O(n), O(n²), or O(n³)
   Example: "In post_process_tgk(), after line ~450, add: if (vnorm < sqrt(n*eps)*sigma_max)
   recompute u[k] = B*v[k]/sigma[k] using one-sided recovery. This is O(n) per vector."

4. **DEAD ENDS** — Read knowledge/PRIOR_APPROACHES.md. List 1-3 techniques already tried
   that failed for THIS specific failure pattern. One line each.

## Output: Plain text, start with "RESEARCH BRIEF:". No code blocks. No JSON.
"""


# ---------------------------------------------------------------------------
# Phase 2: Mutation output format instruction
# ---------------------------------------------------------------------------

OUTPUT_INSTRUCTION = """

## MANDATORY WORKFLOW (follow this exact sequence)

1. **READ** the current bidiag_svd.h (it's in the prompt above) and the research brief
2. **IMPLEMENT** your changes by editing src/bidiag_svd.h directly using Edit/Write tools
3. **COMPILE** — run: g++ -std=c++17 -O2 -Isrc/clapack -o evaluate src/evaluate.cpp -Llib -lxmr_c -lm
4. **TEST** — run: ./evaluate (no args = adversarial only, fast). Check PASS/FAIL counts.
5. **ITERATE** if compile fails or tests regress — fix and re-test. Do NOT give up after one try.
6. **OUTPUT** the final working version as described below.

IMPORTANT RULES:
- NEVER ask questions or stop for confirmation. Just implement, compile, test, iterate.
- Do NOT read knowledge files or reference implementations to "understand the problem."
  The research agent already did that — its diagnosis is in the RESEARCH BRIEF above.
  Your first tool call should be editing src/bidiag_svd.h, not reading other files.
- NEVER add O(n³) operations: no global MGS over all n vectors, no dense n×n matrix multiply,
  no full SVD/eigendecomposition of n×n matrices. Chunked operations (chunk ≤ 32) are OK.

## CRITICAL OUTPUT FORMAT

After testing, output the COMPLETE final bidiag_svd.h in a SINGLE ```cpp block.
This is how the evolution framework captures your work.

- The ```cpp block MUST contain the ENTIRE file — every line, every function
- Do NOT use "// ... rest unchanged" or "// [omitted]" — output EVERYTHING
- Do NOT output multiple ```cpp blocks — only ONE, and it must be the LAST code block
- The file must start with #pragma once and contain the complete bidiag_svd() function

```cpp
#pragma once
// YOUR COMPLETE EVOLVED FILE HERE — EVERY LINE
```
"""


# ---------------------------------------------------------------------------
# Code extraction
# ---------------------------------------------------------------------------

def _extract_code_from_response(response: str, language: str = "cpp") -> str:
    """
    Extract the evolved code from an agent's conversational response.

    Strategy (in order):
    1. Find all ```cpp blocks. Among those with #pragma once (full files),
       take the LAST one (matching our "LAST code block" instruction).
       If none have #pragma once, take the largest.
    2. Find all generic ``` blocks, take the largest
    3. Return raw response (OpenEvolve's parse_full_rewrite will handle it)
    """
    if not response or not response.strip():
        return response

    # Strategy 1: ```cpp blocks
    pattern = r"```" + re.escape(language) + r"\s*\n(.*?)```"
    matches = re.findall(pattern, response, re.DOTALL)
    if matches:
        # Prefer blocks that look like complete files (#pragma once)
        full_files = [m for m in matches if "#pragma once" in m and len(m.strip()) > 500]
        if full_files:
            # Take the LAST complete file (our instruction says "LAST code block")
            best = full_files[-1]
        else:
            # Fallback: take the largest block
            best = max(matches, key=len)

        if len(best.strip()) > 100:
            return f"```{language}\n{best.strip()}\n```"

    # Strategy 2: generic ``` blocks
    pattern = r"```\s*\n(.*?)```"
    matches = re.findall(pattern, response, re.DOTALL)
    if matches:
        # Same logic: prefer last complete file, then largest
        full_files = [m for m in matches if "#pragma once" in m and len(m.strip()) > 500]
        if full_files:
            best = full_files[-1]
        else:
            best = max(matches, key=len)

        if len(best.strip()) > 100:
            return f"```{language}\n{best.strip()}\n```"

    # Strategy 3: raw fallback
    logger.warning(
        "No substantial code block found in agent response "
        f"({len(response)} chars). Returning raw response."
    )
    return response


# ---------------------------------------------------------------------------
# Agent runner (shared by both phases)
# ---------------------------------------------------------------------------

async def _run_agent_session(
    prompt: str,
    options,
    label: str = "agent",
) -> str:
    """Run a single agent session and collect the text result."""
    sdk = _ensure_sdk()
    result_text = ""
    cost_usd = 0.0

    try:
        async for message in sdk.query(prompt=prompt, options=options):
            if hasattr(message, "result") and message.result:
                result_text = message.result
            if hasattr(message, "cost_usd") and message.cost_usd:
                cost_usd = message.cost_usd
            # Fallback: capture text blocks from assistant messages
            if hasattr(message, "content") and not result_text:
                for block in getattr(message, "content", []):
                    if hasattr(block, "text"):
                        result_text += block.text

    except Exception as e:
        logger.error(f"{label} session failed: {e}")
        raise

    if cost_usd > 0:
        logger.info(f"{label} session cost: ${cost_usd:.4f}")

    return result_text


# ---------------------------------------------------------------------------
# Main adapter class
# ---------------------------------------------------------------------------

class ClaudeAgentLLM(LLMInterface):
    """
    Two-phase agentic LLM for OpenEvolve:
      Phase 1: Research agent analyzes failures (read-only, cheap)
      Phase 2: Mutation agent implements changes (full tools)
    """

    def __init__(
        self,
        model: str = "claude-sonnet-4-20250514",
        cwd: Optional[str] = None,
        max_budget_usd: float = 1.0,
        enable_tools: bool = True,
        enable_explore: bool = True,
        enable_checkpointing: bool = False,
        load_project_settings: bool = True,
        permission_mode: str = "acceptEdits",
        extra_system_prompt: Optional[str] = None,
        custom_agents: Optional[dict] = None,
        language: str = "cpp",
        # Research agent settings
        research_model: Optional[str] = None,
        research_budget_usd: float = 0.50,
        enable_research: bool = True,
    ):
        _ensure_sdk()
        self.model = model
        self.cwd = cwd
        self.max_budget_usd = max_budget_usd
        self.enable_tools = enable_tools
        self.enable_explore = enable_explore
        self.enable_checkpointing = enable_checkpointing
        self.load_project_settings = load_project_settings
        self.permission_mode = permission_mode
        self.extra_system_prompt = extra_system_prompt
        self.custom_agents = custom_agents or {}
        self.language = language

        # Research agent config
        self.enable_research = enable_research
        self.research_model = research_model or model
        self.research_budget_usd = research_budget_usd

        logger.info(
            f"ClaudeAgentLLM initialized: model={model}, "
            f"tools={enable_tools}, explore={enable_explore}, "
            f"budget=${max_budget_usd}, "
            f"research={'ON' if enable_research else 'OFF'} "
            f"(model={self.research_model}, budget=${research_budget_usd})"
        )

    def _build_options(
        self,
        system_prompt: Optional[str] = None,
        *,
        read_only: bool = False,
        model_override: Optional[str] = None,
        budget_override: Optional[float] = None,
        tools_override: Optional[list] = None,
    ):
        """Build ClaudeAgentOptions for an agent session."""
        sdk = _ensure_sdk()

        # Tools
        if tools_override is not None:
            allowed_tools = tools_override
        elif read_only:
            allowed_tools = ["Read", "Grep", "Glob", "Bash"]
        elif self.enable_tools:
            allowed_tools = [
                "Read", "Write", "Edit", "MultiEdit",
                "Bash", "Grep", "Glob",
                "WebSearch", "WebFetch",
            ]
            if self.enable_explore:
                allowed_tools.append("Task")
        else:
            allowed_tools = []

        # Settings sources
        setting_sources = ["user", "project"] if self.load_project_settings else []

        # System prompt
        full_system_prompt = system_prompt or ""
        if self.extra_system_prompt and not read_only:
            full_system_prompt += "\n\n" + self.extra_system_prompt

        # Build options
        opts_kwargs = dict(
            model=model_override or self.model,
            allowed_tools=allowed_tools if allowed_tools else None,
            permission_mode="plan" if read_only else self.permission_mode,
            max_budget_usd=budget_override or self.max_budget_usd,
            enable_file_checkpointing=self.enable_checkpointing,
        )

        if full_system_prompt.strip():
            opts_kwargs["system_prompt"] = full_system_prompt

        if self.cwd:
            opts_kwargs["cwd"] = self.cwd

        if setting_sources:
            opts_kwargs["setting_sources"] = setting_sources

        if self.custom_agents and not read_only:
            opts_kwargs["agents"] = self.custom_agents

        return sdk.ClaudeAgentOptions(**opts_kwargs)

    # -------------------------------------------------------------------
    # Phase 1: Research
    # -------------------------------------------------------------------

    async def _run_research(self, openevolve_prompt: str) -> str:
        """
        Run the research agent to analyze failures and produce a brief.

        Args:
            openevolve_prompt: The full prompt OpenEvolve built (code + metrics + artifacts)

        Returns:
            Research brief text, or empty string if research fails/is disabled
        """
        if not self.enable_research:
            return ""

        research_prompt = RESEARCH_PROMPT_TEMPLATE.format(
            openevolve_prompt=openevolve_prompt
        )

        options = self._build_options(
            system_prompt=RESEARCH_SYSTEM_PROMPT,
            read_only=True,
            model_override=self.research_model,
            budget_override=self.research_budget_usd,
        )

        logger.info("Phase 1: Running research agent...")

        try:
            result = await _run_agent_session(
                prompt=research_prompt,
                options=options,
                label="Research",
            )
        except Exception as e:
            logger.warning(f"Research agent failed (non-fatal): {e}")
            return ""

        if not result or not result.strip():
            logger.warning("Research agent returned empty result")
            return ""

        # Truncate if too long (don't overwhelm the mutation prompt)
        if len(result) > 3000:
            result = result[:3000] + "\n... [research brief truncated]"

        logger.info(f"Research brief: {len(result)} chars")
        return result

    # -------------------------------------------------------------------
    # Phase 2: Mutation
    # -------------------------------------------------------------------

    async def generate(self, prompt: str, **kwargs) -> str:
        """Generate text from a prompt using the two-phase pipeline."""
        return await self.generate_with_context(
            system_message=kwargs.get("system_message", ""),
            messages=[{"role": "user", "content": prompt}],
            **kwargs,
        )

    async def generate_with_context(
        self, system_message: str, messages: List[Dict[str, str]], **kwargs
    ) -> str:
        """
        Two-phase generation:
          1) Research agent analyzes failures (read-only)
          2) Mutation agent implements changes (full tools)
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

        # --- Phase 1: Research ---
        research_brief = await self._run_research(openevolve_prompt)

        # --- Phase 2: Mutation ---
        # Structure: research brief FIRST (so agent sees guidance before the code),
        # then the OpenEvolve prompt (code + metrics), then output instructions
        mutation_parts = []

        if research_brief:
            mutation_parts.append(
                "## RESEARCH BRIEF — READ THIS FIRST\n"
                "A research agent analyzed the current failures, knowledge base, and prior\n"
                "approaches. Follow its specific recommendations. Do NOT re-investigate\n"
                "what it already analyzed — go straight to implementing the suggested fix.\n\n"
                + research_brief
            )

        mutation_parts.append(openevolve_prompt)
        mutation_parts.append(OUTPUT_INSTRUCTION)

        mutation_prompt = "\n\n".join(mutation_parts)

        options = self._build_options(system_prompt=system_message)

        logger.info("Phase 2: Running mutation agent...")

        result_text = await _run_agent_session(
            prompt=mutation_prompt,
            options=options,
            label="Mutation",
        )

        if not result_text:
            logger.warning("Mutation agent returned empty result")
            return ""

        # Post-process: extract code from conversational response
        processed = _extract_code_from_response(result_text, self.language)

        logger.debug(
            f"Mutation response: {len(result_text)} chars raw -> "
            f"{len(processed)} chars processed"
        )

        return processed


# ---------------------------------------------------------------------------
# Factory function for OpenEvolve's init_client hook
# ---------------------------------------------------------------------------

def create_claude_agent_client(model_cfg) -> ClaudeAgentLLM:
    """
    Factory function compatible with OpenEvolve's LLMModelConfig.init_client.
    """
    return ClaudeAgentLLM(
        model=getattr(model_cfg, "name", "claude-sonnet-4-20250514"),
        cwd=getattr(model_cfg, "cwd", None),
        max_budget_usd=getattr(model_cfg, "max_budget_usd", 1.0),
        enable_tools=getattr(model_cfg, "enable_tools", True),
        enable_explore=getattr(model_cfg, "enable_explore", True),
        enable_checkpointing=getattr(model_cfg, "enable_checkpointing", False),
        load_project_settings=getattr(model_cfg, "load_project_settings", True),
        permission_mode=getattr(model_cfg, "permission_mode", "acceptEdits"),
        extra_system_prompt=getattr(model_cfg, "extra_system_prompt", None),
        custom_agents=getattr(model_cfg, "custom_agents", None),
        language=getattr(model_cfg, "language", "cpp"),
        # Research agent settings
        research_model=getattr(model_cfg, "research_model", None),
        research_budget_usd=getattr(model_cfg, "research_budget_usd", 0.50),
        enable_research=getattr(model_cfg, "enable_research", True),
    )
