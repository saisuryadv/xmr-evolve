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
from claude_agent_llm import create_claude_agent_client, save_agent_settings


async def main():
    parser = argparse.ArgumentParser(
        description="OpenEvolve + Claude Agent SDK for Bidiagonal SVD"
    )
    parser.add_argument("--config", default="config.yaml", help="YAML config file")
    parser.add_argument("--iterations", type=int, default=50)
    parser.add_argument("--checkpoint", default=None, help="Resume from checkpoint")

    # Claude Agent SDK settings
    parser.add_argument(
        "--model", default="claude-opus-4-6",
        help="Claude model for evolution agent"
    )
    parser.add_argument(
        "--budget", type=float, default=None,
        help="Max USD per evolution agent call (default: unlimited)"
    )
    parser.add_argument(
        "--effort", default="high",
        choices=["low", "medium", "high", "max"],
        help="Thinking effort level"
    )
    parser.add_argument(
        "--permission-mode", default="bypassPermissions",
        choices=["plan", "acceptEdits", "bypassPermissions"],
        help="Agent permission mode"
    )
    parser.add_argument(
        "--experiment", default=None,
        help="Experiment name (organizes scratch/ and output under this name)"
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
    agent_model.effort = args.effort
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

    # Persist agent settings to disk so worker processes can read them
    # (OpenEvolve serializes config via asdict() which drops custom attributes)
    project_dir = os.path.dirname(os.path.abspath(__file__))

    # Organize scratch, variants, and output by experiment name
    if args.experiment:
        scratch_dir = os.path.join(project_dir, "scratch", args.experiment)
        variants_dir = os.path.join(project_dir, "evolved_variants", args.experiment)
        config.log_dir = os.path.join("src", "openevolve_output", args.experiment, "logs")
        output_base = os.path.join(project_dir, "src", "openevolve_output", args.experiment)
        os.makedirs(output_base, exist_ok=True)
    else:
        scratch_dir = os.path.join(project_dir, "scratch")
        variants_dir = os.path.join(project_dir, "evolved_variants")

    # Set variants dir for evaluator (reads env var at module level)
    os.environ["XMR_VARIANTS_DIR"] = variants_dir

    agent_model.scratch_dir = scratch_dir
    agent_model.variants_dir = variants_dir
    agent_settings = {
        "model": args.model,
        "cwd": project_dir,
        "scratch_dir": scratch_dir,
        "variants_dir": variants_dir,
        "max_budget_usd": args.budget,
        "effort": args.effort,
        "enable_tools": True,
        "enable_checkpointing": False,
        "load_project_settings": True,
        "permission_mode": args.permission_mode,
        "extra_system_prompt": system_prompt,
        "custom_agents": agent_model.custom_agents,
        "language": "cpp",
    }
    # Save to experiment scratch dir
    save_agent_settings(scratch_dir, agent_settings)
    # Also save to canonical location (scratch/.agent_settings.json) so worker
    # processes can auto-discover it — OpenEvolve's ProcessPoolExecutor drops
    # custom attributes during serialization, and the worker's auto-discover
    # only looks at scratch/.agent_settings.json
    canonical_scratch = os.path.join(project_dir, "scratch")
    if canonical_scratch != scratch_dir:
        save_agent_settings(canonical_scratch, agent_settings)

    # Set models in config
    config.llm.models = [agent_model]

    # Paths
    initial_program = os.path.join(os.path.dirname(__file__), "src", "bidiag_svd.h")
    evaluator_file = os.path.join(os.path.dirname(__file__), "evaluator.py")

    # Run
    print(f"{'='*60}")
    print(f"  OpenEvolve + Claude Agent SDK — Bidiagonal SVD")
    budget_str = f"${args.budget:.2f}/call" if args.budget else "unlimited"
    print(f"  Model:       {args.model} (effort={args.effort}, {budget_str})")
    print(f"  Permission:  {args.permission_mode}")
    print(f"  Iterations:  {args.iterations}")
    exp_label = args.experiment or "(default)"
    print(f"  Experiment:  {exp_label}")
    print(f"  Scratch:     {scratch_dir}")
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
        "\n## THE PROBLEM\n",
        "Given an upper bidiagonal matrix B (n×n) with diagonal entries d[0..n-1] and\n",
        "superdiagonal entries e[0..n-2], compute the **singular value decomposition**:\n",
        "  B = U · Σ · V^T\n",
        "where Σ = diag(σ_1 ≥ σ_2 ≥ ... ≥ σ_n ≥ 0), U and V are orthogonal.\n",
        "\n",
        "**The challenge**: Do this in **O(n²) worst-case time** while achieving\n",
        "LAPACK-quality accuracy (residual ≤ 7·n·ε·‖B‖, orthogonality ≤ 5·n·ε).\n",
        "This is an **open problem** — DBDSQR (QR iteration) achieves the accuracy but\n",
        "is O(n³); MR³-based methods achieve O(n²) but fail accuracy on adversarial inputs.\n",
        "\n",
        "**Why it's hard**: The MR³ algorithm computes eigenvectors of the 2n×2n\n",
        "Golub-Kahan (TGK) tridiagonal matrix, then extracts u/v from the eigenvectors.\n",
        "This extraction requires the eigenvector subspace to have \"GK structure\" —\n",
        "a property that shifting (core to MR³) can destroy. The NCD (Nearly Constant\n",
        "Diagonal) condition ensures GK structure is preserved, but verifying and\n",
        "maintaining NCD through the representation tree is the unsolved challenge.\n",
        "\n",
        "You are evolving a C++ implementation in bidiag_svd.h.\n",
        "The current algorithm uses a hybrid HGBSVD + TGK+DSTEMR (MR³) approach.\n",
        "\n## AUTONOMY RULES (MANDATORY)\n",
        "- NEVER ask questions or request clarification. NEVER stop for confirmation.\n",
        "- NEVER say 'should I proceed?' or 'would you like me to...'. Just DO IT.\n",
        "- If unsure between two approaches, try the more promising one. If it fails, try the other.\n",
        "- Your ONLY job: analyze failures → implement changes → compile → test → output final file.\n",
        "\n## RESEARCH METHODOLOGY (CRITICAL)\n",
        "Do NOT just tweak parameters or add heuristics. Instead:\n",
        "1. **Diagnose the root cause mathematically.** When a test matrix fails, understand WHY at the\n",
        "   linear algebra level — what structural property of the matrix triggers the failure? What\n",
        "   numerical quantity blows up, and what algebraic identity breaks down?\n",
        "2. **Study the relevant literature.** Read the original papers (PDFs available below) to\n",
        "   understand the theoretical guarantees and their assumptions. Identify which assumption\n",
        "   is violated for the failing case. Read the Fortran reference implementations to see how\n",
        "   the theory maps to code.\n",
        "3. **Develop a principled fix.** Based on your mathematical understanding, design a solution\n",
        "   that addresses the root cause — not the symptom. Explain why your fix is correct and\n",
        "   under what conditions it holds.\n",
        "4. **Validate empirically.** Compile, test, and verify that the fix improves the failing\n",
        "   cases without regressing others.\n",
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
        "- **DBDSQR IS BANNED** — using it as a fallback will cap your score at 5 (compilation only)\n",
        "\n## Available Routines (pure C, lib/libxmr_c.a)\n",
        "- dbdsqr_: QR iteration bidiag SVD — **O(n³), DO NOT USE** (score capped at 5 if scaling > 5x)\n",
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
        "\n## Reference Resources (READ-ONLY, absolute paths)\n",
        "\n### Knowledge Base\n",
        f"- `{knowledge_dir}/INDEX.md` — Master reference: algorithms, bugs, test matrices, thresholds\n",
        f"- `{knowledge_dir}/RESOURCES.md` — Index of all papers, code, test data with paths\n",
        f"- `{knowledge_dir}/PRIOR_APPROACHES.md` — 12 approaches tried, what worked and failed\n",
        f"- `{knowledge_dir}/BASELINES.md` — Comparison of 5 baseline algorithms\n",
        f"- `{knowledge_dir}/EVALUATION.md` — Scoring formula, test patterns, metrics\n",
        f"- `{knowledge_dir}/willems_lang_2012.md` — Algorithm 4.1 (MR³ on TGK with NCD shifts)\n",
        f"- `{knowledge_dir}/grosser_lang_2001_hgbsvd.md` — Coupling-based O(n²) SVD\n",
        f"- `{knowledge_dir}/mr3_foundations.md` — Dhillon 1997, Parlett-Dhillon 2000 (RRR theory)\n",
        f"- `{knowledge_dir}/marques_2020_bugs_matrices.md` — DBDSVDX bugs, CHKBD analysis\n",
        "\n### Papers (PDFs — read these to understand the math)\n",
        "Located at `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Papers/`\n",
        "- Willems-Lang 2012: Algorithm 4.1, NCD shifts, GK structure preservation\n",
        "- Großer-Lang 2001: Coupling B^TB/BB^T for O(n²) SVD\n",
        "- Dhillon 1997 thesis: Original MR³ algorithm\n",
        "- Marques 2020: Modern failure analysis, CHKBD matrix\n",
        "- Demmel-Kahan 1990: Foundational bidiag SVD accuracy\n",
        "\n### Reference Fortran Code\n",
        "- XMR (Willems): `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Code/tdsolver/xmr/SRC/O/` (45 files, `dstexr.f` is master)\n",
        "- HGBSVD (Großer-Lang): `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Code/hgbsvd/hgbsvd/v2/` (24 files, `dbdsgr.f` is master)\n",
        "- LAPACK: `/Users/saisurya/MRRR/xmr-evolve/lapack/SRC/` (dbdsqr.f, dstemr.f, dlarrv.f, etc.)\n",
        "\n### Reference C++ Implementations\n",
        "- DBDSQR baseline: `src/bidiag_dbdsqr.h` (270/270 pass, O(n³), score=5 hard gate)\n",
        "- HGBSVD baseline: `src/bidiag_hgbsvd.h` (135/270 pass, O(n²), score=68.6)\n",
        "- TGK+STEMR baseline: `src/bidiag_tgk_stemr.h` (65/270 pass, O(n²), score=51.2)\n",
        "\n### Test Matrices\n",
        "- STCollection: `/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/STCollection/DATA/B_*.dat` (19 matrices)\n",
        "\n### DO NOT READ (irrelevant)\n",
        "- Do NOT read `approach_*.py` files — Python prototypes from a different experiment\n",
        "- Do NOT read `.claude/` memory files, `evaluate_alphaevolve.py`, `mr3_tridiag.py`\n",
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
