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
    parser.add_argument(
        "--prompt", default="xmr_debug",
        choices=["exploration", "xmr_debug"],
        help="Prompt template to use (default: xmr_debug)"
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
    # Always create an experiment directory (auto-generate timestamp if not provided)
    if args.experiment:
        exp_name = args.experiment
    else:
        from datetime import datetime
        exp_name = datetime.now().strftime("exp_%Y%m%d_%H%M%S")

    scratch_dir = os.path.join(project_dir, "scratch", exp_name)
    variants_dir = os.path.join(project_dir, "evolved_variants", exp_name)
    config.log_dir = os.path.join("src", "openevolve_output", exp_name, "logs")
    output_base = os.path.join(project_dir, "src", "openevolve_output", exp_name)
    os.makedirs(output_base, exist_ok=True)

    # Set variants dir for evaluator (reads env var at module level)
    os.environ["XMR_VARIANTS_DIR"] = variants_dir

    agent_model.scratch_dir = scratch_dir
    agent_model.variants_dir = variants_dir

    # Determine initial program path
    if args.prompt == "xmr_debug":
        initial_prog = os.path.join(project_dir, "recovered", "version_013_pass138.h")
    else:
        initial_prog = os.path.join(project_dir, "src", "bidiag_svd.h")
    agent_model.initial_program = initial_prog

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
        "custom_agents": agent_model.custom_agents,
        "language": "cpp",
        "prompt_template": args.prompt,
        "initial_program": initial_prog,
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
    # XMR debug starts from v013 (best O(n²)-compliant version from prior agent, 138/289, score=64.3);
    # exploration starts from the current hybrid bidiag_svd.h
    if args.prompt == "xmr_debug":
        initial_program = os.path.join(os.path.dirname(__file__), "recovered", "version_013_pass138.h")
    else:
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
    print(f"  Prompt:      {args.prompt}")
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



if __name__ == "__main__":
    asyncio.run(main())
