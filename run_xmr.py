"""
run_xmr.py — Run OpenEvolve with Claude Agent SDK for bidiagonal SVD evolution

Single-phase agentic pipeline per generation:
  One agent with full tools (read, write, edit, bash, grep, glob, search, fetch,
  subagents) and directory-scoped write permissions (scratch/ + src/bidiag_svd.h).
  The agent handles both research/analysis and code mutation in a single session.

Usage:
    python run_xmr.py --iterations 50 --budget 3.50
    python run_xmr.py --iterations 10 --budget 5.0 --model claude-opus-4-6
"""

import argparse
import asyncio
import os
import sys
from pathlib import Path

from openevolve.config import Config, LLMModelConfig, load_config
from openevolve.controller import OpenEvolve

# Import the Claude Agent SDK adapter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from claude_agent_llm import create_claude_agent_client


async def main():
    parser = argparse.ArgumentParser(
        description="OpenEvolve + Claude Agent SDK for Bidiagonal SVD"
    )
    parser.add_argument("--config", default="config.yaml", help="YAML config file")
    parser.add_argument("--iterations", type=int, default=50)
    parser.add_argument("--checkpoint", default=None, help="Resume from checkpoint")

    # Claude Agent SDK settings
    parser.add_argument(
        "--model", default="claude-sonnet-4-20250514",
        help="Claude model for evolution agent"
    )
    parser.add_argument(
        "--budget", type=float, default=3.50,
        help="Max USD per evolution agent call"
    )
    parser.add_argument(
        "--permission-mode", default="bypassPermissions",
        choices=["plan", "acceptEdits", "bypassPermissions"],
        help="Agent permission mode"
    )

    args = parser.parse_args()

    # Load base config
    config_path = os.path.join(os.path.dirname(__file__), args.config)
    if os.path.exists(config_path):
        config = load_config(config_path)
    else:
        config = Config()

    # Disable OpenEvolve's built-in LLM feedback — the agent does its own
    # analysis as part of the unified workflow
    config.evaluator.use_llm_feedback = False

    # Build system prompt with knowledge base
    knowledge_dir = os.path.join(os.path.dirname(__file__), "knowledge")
    system_prompt = build_system_prompt(knowledge_dir)

    # Configure the evolution agent (Claude Agent SDK)
    agent_model = LLMModelConfig(
        name=args.model,
        weight=1.0,
        init_client=create_claude_agent_client,
        api_base="unused",
        api_key="unused",
        temperature=0.7,
        max_tokens=16000,
        timeout=600,
        retries=2,
        retry_delay=10,
    )

    # Agent-specific settings
    agent_model.cwd = os.path.dirname(os.path.abspath(__file__))
    agent_model.max_budget_usd = args.budget
    agent_model.enable_tools = True
    agent_model.enable_checkpointing = False
    agent_model.load_project_settings = True
    agent_model.permission_mode = args.permission_mode
    agent_model.extra_system_prompt = system_prompt
    agent_model.language = "cpp"

    agent_model.custom_agents = {
        "numerical-analyst": {
            "description": "Analyze numerical algorithm correctness and stability",
            "prompt": (
                "You are a numerical linear algebra expert. Analyze the given "
                "code for numerical stability issues: element growth, loss of "
                "orthogonality, catastrophic cancellation. Reference Willems-Lang "
                "2012 and the MR³ literature. Report specific failure modes."
            ),
            "tools": ["Read", "Bash", "Grep", "WebSearch"],
            "model": "claude-sonnet-4-20250514",
        },
        "bug-hunter": {
            "description": "Search for known bugs and edge cases in bidiagonal SVD",
            "prompt": (
                "You are a bug hunter for bidiagonal SVD algorithms. Read the "
                "knowledge base files in the knowledge/ directory to find all "
                "documented bugs. Check if the current implementation has any "
                "of these bugs. Propose specific test matrices that would expose "
                "failures."
            ),
            "tools": ["Read", "Bash", "Grep", "Glob"],
            "model": "claude-haiku-4-5-20251001",
        },
    }

    # Set models in config
    config.llm.models = [agent_model]

    # Paths
    initial_program = os.path.join(os.path.dirname(__file__), "src", "bidiag_svd.h")
    evaluator_file = os.path.join(os.path.dirname(__file__), "evaluator.py")

    # Run
    print(f"{'='*60}")
    print(f"  OpenEvolve + Claude Agent SDK — Bidiagonal SVD")
    print(f"  Model:       {args.model} (${args.budget:.2f}/call)")
    print(f"  Permission:  {args.permission_mode}")
    print(f"  Iterations:  {args.iterations}")
    print(f"  Write scope: scratch/agent_NNN/ (isolated per call)")
    print(f"  Subagents:   numerical-analyst, bug-hunter")
    print(f"{'='*60}")

    evolve = OpenEvolve(
        initial_program_path=initial_program,
        evaluation_file=evaluator_file,
        config=config,
    )

    if args.checkpoint:
        evolve._load_checkpoint(args.checkpoint)

    best = await evolve.run(iterations=args.iterations)

    if best:
        print(f"\nBest program metrics:")
        for name, value in best.metrics.items():
            if isinstance(value, float):
                print(f"  {name}: {value:.4f}")
            else:
                print(f"  {name}: {value}")
    else:
        print("\nNo valid programs found.")


def build_system_prompt(knowledge_dir: str) -> str:
    """Build the system prompt from knowledge base files."""
    prompt_parts = [
        "# Bidiagonal SVD via MR³ — Evolution Context\n",
        "You are evolving a C++ bidiagonal SVD implementation in bidiag_svd.h.\n",
        "The algorithm uses a hybrid HGBSVD + TGK+DSTEMR (MR³) approach.\n",
        "\n## AUTONOMY RULES (MANDATORY)\n",
        "- NEVER ask questions or request clarification. NEVER stop for confirmation.\n",
        "- NEVER say 'should I proceed?' or 'would you like me to...'. Just DO IT.\n",
        "- If unsure between two approaches, try the more promising one. If it fails, try the other.\n",
        "- Your ONLY job: analyze failures → implement changes → compile → test → output final file.\n",
        "\n## Your Goal\n",
        "Improve the `bidiag_svd()` function to:\n",
        "1. Pass ALL test matrices (STCollection + adversarial) with res ≤ 7, orthoU ≤ 5, orthoV ≤ 5\n",
        "2. Maintain O(n²) worst-case complexity (scaling ratio ≤ 5x when n doubles)\n",
        "3. Handle extreme condition numbers (up to 10^15)\n",
        "\n## Write Permissions\n",
        "You have an isolated workspace directory (provided in the prompt).\n",
        "You can ONLY write/edit files inside your workspace. All other files are READ-ONLY.\n",
        "Your output file is `<workspace>/program/bidiag_svd.h` — edit that, not src/.\n",
        "\n## Key Constraints\n",
        "- The file is bidiag_svd.h — a single self-contained C++ header (no #include of project headers)\n",
        "- You can call Fortran routines via extern \"C\" (DSTEMR, DBDSGR, DBDSQR, etc.)\n",
        "- dbdsgr_ (f2c) requires ftnlen args: pass 1, 1 at end of each call\n",
        "- dstemr_ does NOT need ftnlen args\n",
        "\n## O(n²) Checklist — FORBIDDEN Operations\n",
        "- Global MGS/reorthogonalization over all n vectors → O(n³)\n",
        "- Dense n×n matrix multiply (DGEMM on full U or V) → O(n³)\n",
        "- Full eigendecomposition of an n×n matrix → O(n³)\n",
        "- Nested loops: for i in n { for j in n { for k in n } } → O(n³)\n",
        "- OK: chunked MGS with chunk_size ≤ 32 → O(n² · 32) = O(n²)\n",
        "- OK: O(n) one-sided recovery per vector → O(n²) total\n",
        "- OK: DSTEMR on 2n×2n TGK → O(n²) (MR³ is O(n) per eigenvector)\n",
        "\n## Available Routines (pure C, lib/libxmr_c.a)\n",
        "- dbdsqr_: QR iteration bidiag SVD (O(n³) — use only as last resort)\n",
        "- dstemr_: LAPACK MR³ tridiag eigensolver (O(n²))\n",
        "- dstexr_: Willems XMR improved MR³ (O(n²))\n",
        "- dbdsgr_: Großer-Lang coupling SVD (O(n²), needs ftnlen 1,1)\n",
        "- dlamch_: machine parameters. BLAS: dnrm2, ddot, dscal, dcopy, daxpy, dgemv\n",
        "\n## Build & Test\n",
        "```bash\n",
        "g++ -std=c++17 -O2 -Isrc/clapack -o evaluate src/evaluate.cpp -Llib -lxmr_c -lm\n",
        "./evaluate\n",
        "```\n",
        "\n## Current Architecture (hybrid HGBSVD + TGK+DSTEMR)\n",
        "1. Try HGBSVD (dbdsgr_) first — O(n²), good accuracy on passing tests\n",
        "2. If HGBSVD returns INFO!=0, fall back to TGK+DSTEMR with post-processing\n",
        "3. Post-processing: normalize, sign-fix, one-sided recovery (U=BV/σ), chunked MGS\n",
        "\n## Known Failure Modes\n",
        "- stemr_killer, wilkinson_like: TGK extraction loses GK structure\n",
        "- two_clusters, one_big_cluster: orthogonality fails within clusters\n",
        "- exponential_graded: near-zero σ extraction diverges\n",
        "- spike, saw_tooth: degenerate spectrum trips HGBSVD (INFO!=0)\n",
        "- chkbd: HGBSVD returns info!=0, TGK fallback has reortho issues\n",
        "\n",
    ]

    # Include condensed knowledge if available
    index_file = os.path.join(knowledge_dir, "INDEX.md")
    if os.path.exists(index_file):
        with open(index_file) as f:
            content = f.read()
        if len(content) > 15000:
            content = content[:15000] + "\n... [truncated]\n"
        prompt_parts.append("## Knowledge Base Summary\n")
        prompt_parts.append(content)
        prompt_parts.append("\n")

    progression_file = os.path.join(knowledge_dir, "PROGRESSION.md")
    if os.path.exists(progression_file):
        with open(progression_file) as f:
            content = f.read()
        if len(content) > 5000:
            content = content[:5000] + "\n... [truncated]\n"
        prompt_parts.append("## Paper Progression\n")
        prompt_parts.append(content)

    return "".join(prompt_parts)


if __name__ == "__main__":
    asyncio.run(main())
