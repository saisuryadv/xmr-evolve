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

## YOUR WORKFLOW

### Phase 1: Analyze (~20% of budget)
1. **READ PRIOR RESEARCH** — Check the prior research summaries above. If any are
   listed, read the most relevant RESEARCH.md files to understand what was tried,
   what worked, and what failed. Don't repeat failed approaches.

2. **IDENTIFY FAILING TESTS** — From the evaluation output below, list which
   patterns fail and group by failure mode (INFO!=0, residual blow-up, orthoU/V blow-up).

3. **DIAGNOSE ROOT CAUSE** — For the top 1-3 failure modes, identify the specific
   code path in bidiag_svd.h that triggers the failure. Be concrete: name the function,
   the condition, the formula. Read up to 2 knowledge files if needed (start with
   knowledge/INDEX.md). Use WebSearch for LAPACK/MR³ techniques if local knowledge
   is insufficient.

### Phase 2: Implement (~60% of budget)
4. **EDIT** `{agent_dir}/program/bidiag_svd.h` with your fix.
5. **COMPILE** and **TEST** using the commands above.
6. **ITERATE** if compile fails or tests regress — fix and re-test.

### Phase 3: Document & Finalize (~20% of budget)
7. Verify your final version compiles and passes more tests than the parent.
8. Your final `{agent_dir}/program/bidiag_svd.h` IS your output — no need to
   paste it in chat. Just make sure it compiles and passes tests.
9. **Write `{agent_dir}/RESEARCH.md`** — a structured research report. This is
   critical: future agents will read your RESEARCH.md to build on your work.
   Follow this structure:

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

## RULES
- NEVER ask questions or stop for confirmation. Just implement, compile, test, iterate.
- NEVER add O(n³) operations: no global MGS over all n vectors, no dense n×n matrix multiply,
  no full SVD/eigendecomposition of n×n matrices. Chunked operations (chunk ≤ 32) are OK.
- You can spawn subagents (Agent tool) for parallel investigation or focused tasks.
- You can create scratch files in `{agent_dir}/` for notes, intermediate results, etc.
- You can ONLY write/edit files inside `{agent_dir}/`. All other files are read-only.

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

async def _run_agent_session(
    prompt: str,
    options,
    label: str = "agent",
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

    try:
        async for message in sdk.query(prompt=actual_prompt, options=options):
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
        raise

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
        model: str = "claude-sonnet-4-20250514",
        cwd: Optional[str] = None,
        max_budget_usd: float = 3.50,
        enable_tools: bool = True,
        enable_checkpointing: bool = False,
        load_project_settings: bool = True,
        permission_mode: str = "bypassPermissions",
        extra_system_prompt: Optional[str] = None,
        custom_agents: Optional[dict] = None,
        language: str = "cpp",
    ):
        _ensure_sdk()
        self.model = model
        self.cwd = cwd
        self.max_budget_usd = max_budget_usd
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

        # Scratch root
        if self.cwd:
            self._scratch_dir = os.path.abspath(os.path.join(self.cwd, "scratch"))
            os.makedirs(self._scratch_dir, exist_ok=True)
        else:
            self._scratch_dir = None

        # Current agent's directory (set per-call for the permission callback)
        self._current_agent_dir: Optional[str] = None

        logger.info(
            f"ClaudeAgentLLM initialized: model={model}, "
            f"tools={enable_tools}, budget=${max_budget_usd}, "
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

        logger.info(f"Created agent workspace: {agent_dir}")
        return agent_dir

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
        variants_dir = os.path.join(self.cwd, "evolved_variants")
        os.makedirs(variants_dir, exist_ok=True)
        pending = os.path.join(variants_dir, ".pending_research.md")
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
        if not self.cwd:
            return "No prior research available."

        variants_dir = os.path.join(self.cwd, "evolved_variants")
        if not os.path.exists(variants_dir):
            return "No prior research available (first iteration)."

        # Find all RESEARCH.md files, sorted by variant number
        import glob as globmod
        research_files = globmod.glob(os.path.join(variants_dir, "v*_RESEARCH.md"))
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
            f"{len(numbered)} prior research reports available at `{variants_dir}/`.\n"
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

        # Full tool access — write restrictions enforced via can_use_tool
        allowed_tools = [
            "Read", "Write", "Edit", "Bash", "Grep", "Glob",
            "WebSearch", "WebFetch", "Agent",
        ]

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

        # Build options
        opts_kwargs = dict(
            model=self.model,
            allowed_tools=allowed_tools if self.enable_tools else [],
            permission_mode=self.permission_mode,
            max_budget_usd=self.max_budget_usd,
            enable_file_checkpointing=self.enable_checkpointing,
            can_use_tool=self._write_permission_callback,
        )

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
        if not self.cwd:
            return ""

        changelog_path = os.path.join(self.cwd, "evolved_variants", "CHANGELOG.md")
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

        variants_dir = os.path.join(self.cwd, "evolved_variants")

        return (
            f"## Evolution History ({num_variants} variants so far)\n\n"
            f"Previous variants saved at: `{variants_dir}/`\n"
            f"- Code: `{variants_dir}/v{{N}}_{{score}}.h` — Read any variant to see its full code\n"
            f"- Research: `{variants_dir}/v{{N}}_RESEARCH.md` — Structured research report\n\n"
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

            logger.info(f"Running evolution agent in {agent_dir}")

            result_text = await _run_agent_session(
                prompt=prompt,
                options=options,
                label="Evolution",
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
    """
    return ClaudeAgentLLM(
        model=getattr(model_cfg, "name", "claude-sonnet-4-20250514"),
        cwd=getattr(model_cfg, "cwd", None),
        max_budget_usd=getattr(model_cfg, "max_budget_usd", 3.50),
        enable_tools=getattr(model_cfg, "enable_tools", True),
        enable_checkpointing=getattr(model_cfg, "enable_checkpointing", False),
        load_project_settings=getattr(model_cfg, "load_project_settings", True),
        permission_mode=getattr(model_cfg, "permission_mode", "bypassPermissions"),
        extra_system_prompt=getattr(model_cfg, "extra_system_prompt", None),
        custom_agents=getattr(model_cfg, "custom_agents", None),
        language=getattr(model_cfg, "language", "cpp"),
    )
