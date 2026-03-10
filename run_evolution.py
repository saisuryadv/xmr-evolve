"""
run_evolution.py — Run OpenEvolve with Claude Agent SDK backend

This script shows how to wire the Claude Agent SDK into OpenEvolve so that
each "generation" step is a full agentic session (with file reading, bash,
explore subagents, etc.) instead of a single API call.

Setup
-----
    pip install openevolve claude-agent-sdk
    npm install -g @anthropic-ai/claude-code

    # For Vertex AI:
    export CLAUDE_CODE_USE_VERTEX=1
    export CLOUD_ML_REGION=global
    export ANTHROPIC_VERTEX_PROJECT_ID=your-project-id
    gcloud auth application-default login

    # Or for direct Anthropic API:
    export ANTHROPIC_API_KEY=sk-ant-...

Usage
-----
    python run_evolution.py \\
        examples/function_minimization/initial_program.py \\
        examples/function_minimization/evaluator.py \\
        --iterations 50
"""

import argparse
import asyncio
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# 1. Import OpenEvolve internals
# ---------------------------------------------------------------------------
from openevolve.config import Config, LLMModelConfig, load_config
from openevolve.controller import OpenEvolve

# 2. Import our adapter
from claude_agent_llm import create_claude_agent_client


async def main():
    parser = argparse.ArgumentParser(
        description="OpenEvolve + Claude Agent SDK"
    )
    parser.add_argument("initial_program", help="Path to initial program")
    parser.add_argument("evaluator", help="Path to evaluator script")
    parser.add_argument("--config", default=None, help="YAML config file")
    parser.add_argument("--iterations", type=int, default=50)
    parser.add_argument("--checkpoint", default=None, help="Resume from checkpoint")

    # Claude Agent SDK specific
    parser.add_argument(
        "--model", default="claude-sonnet-4-20250514",
        help="Claude model (default: claude-sonnet-4-20250514)"
    )
    parser.add_argument(
        "--budget", type=float, default=1.0,
        help="Max USD per generation call (default: 1.0)"
    )
    parser.add_argument(
        "--no-tools", action="store_true",
        help="Disable agentic tools (bare text generation only)"
    )
    parser.add_argument(
        "--no-explore", action="store_true",
        help="Disable Explore/Plan subagents"
    )
    parser.add_argument(
        "--permission-mode", default="acceptEdits",
        choices=["plan", "acceptEdits", "bypassPermissions"],
        help="Agent permission mode (default: acceptEdits)"
    )

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # 3. Load base config (from YAML if provided, else defaults)
    # ------------------------------------------------------------------
    if args.config:
        config = load_config(args.config)
    else:
        config = Config()

    # ------------------------------------------------------------------
    # 4. Override the LLM config to use Claude Agent SDK
    # ------------------------------------------------------------------
    # Build a model config that uses our adapter factory
    agent_model_cfg = LLMModelConfig(
        name=args.model,
        weight=1.0,
        init_client=create_claude_agent_client,  # <-- THE KEY LINE
        # Standard OpenEvolve fields (not used by our adapter, but required)
        api_base="unused",
        api_key="unused",
        temperature=0.7,
        max_tokens=16000,
        timeout=300,
        retries=2,
        retry_delay=5,
    )

    # Attach Agent SDK specific settings as extra attributes
    agent_model_cfg.cwd = str(Path(args.initial_program).parent.resolve())
    agent_model_cfg.max_budget_usd = args.budget
    agent_model_cfg.enable_tools = not args.no_tools
    agent_model_cfg.enable_explore = not args.no_explore
    agent_model_cfg.enable_checkpointing = False
    agent_model_cfg.load_project_settings = True
    agent_model_cfg.permission_mode = args.permission_mode
    agent_model_cfg.extra_system_prompt = None
    agent_model_cfg.custom_agents = None

    # Replace the model list in config
    config.llm.models = [agent_model_cfg]

    # Also set evaluator models (can use same agent or a cheaper one)
    if config.llm.evaluator_models:
        # Keep evaluator models as-is if you want separate eval LLM
        pass
    else:
        config.llm.evaluator_models = [agent_model_cfg]

    # ------------------------------------------------------------------
    # 5. Run OpenEvolve with the agent-powered LLM
    # ------------------------------------------------------------------
    print(f"{'='*60}")
    print(f"  OpenEvolve + Claude Agent SDK")
    print(f"  Model:       {args.model}")
    print(f"  Budget/call: ${args.budget:.2f}")
    print(f"  Tools:       {'enabled' if not args.no_tools else 'disabled'}")
    print(f"  Explore:     {'enabled' if not args.no_explore else 'disabled'}")
    print(f"  Permission:  {args.permission_mode}")
    print(f"  Iterations:  {args.iterations}")
    print(f"{'='*60}")

    evolve = OpenEvolve(
        initial_program_path=args.initial_program,
        evaluation_file=args.evaluator,
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
