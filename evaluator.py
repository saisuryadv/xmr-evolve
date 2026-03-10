"""
OpenEvolve evaluator for the bidiagonal SVD via TGK + MR³.
Compiles the C++ initial program, runs the evaluation, parses metrics.

Cascade evaluation:
  Stage 1: Compile check + small matrices (n≤50) — fast filter
  Stage 2: All STCollection (19 matrices) + small adversarial (n=100)
  Stage 3: Full adversarial (n=200) + scaling test
"""

import subprocess
import os
import re
import tempfile
import shutil
import time
import glob as globmod
from datetime import datetime

from openevolve.evaluation_result import EvaluationResult

# Paths
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.join(PROJECT_DIR, "lib")
STCOLL_DIR = os.path.join(PROJECT_DIR, "..", "bidiag-algo", "MRRR Resources", "STCollection", "DATA")
EVAL_SRC = os.path.join(PROJECT_DIR, "src", "evaluate.cpp")
VARIANTS_DIR = os.environ.get("XMR_VARIANTS_DIR", os.path.join(PROJECT_DIR, "evolved_variants"))
CHANGELOG_FILE = os.path.join(VARIANTS_DIR, "CHANGELOG.md")

# Compiler settings — pure C/C++ (no Fortran, no Accelerate)
CXX = "g++"
CLAPACK_DIR = os.path.join(PROJECT_DIR, "src", "clapack")
CXXFLAGS = ["-std=c++17", "-O2", f"-I{CLAPACK_DIR}"]
LDFLAGS = [f"-L{LIB_DIR}", "-lxmr_c", "-lm"]


# Default metrics for failed evaluations — must include all MAP-Elites
# feature dimensions so the database doesn't crash
FAIL_METRICS = {
    "compile_success": 0,
    "pass_rate": 0.0,
    "pass_worst_scaling": 99.0,
    "pass_avg_residual": 999999.0,
    "pass_avg_ortho_u": 999999.0,
    "pass_avg_ortho_v": 999999.0,
}


def compile_program(program_path: str, output_binary: str) -> dict:
    """Compile the evolved bidiag_svd.h and evaluate.cpp into a binary."""
    build_dir = os.path.dirname(output_binary)

    # Copy both evaluate.cpp and the evolved program into build_dir.
    # evaluate.cpp does #include "bidiag_svd.h" — GCC searches the directory
    # of the including file first, so both must be in the same directory.
    # OpenEvolve may name the temp file something other than bidiag_svd.h,
    # so we always copy it as bidiag_svd.h.
    eval_copy = os.path.join(build_dir, "evaluate.cpp")
    shutil.copy2(EVAL_SRC, eval_copy)

    svd_copy = os.path.join(build_dir, "bidiag_svd.h")
    shutil.copy2(program_path, svd_copy)

    cmd = [CXX] + CXXFLAGS + [
        "-o", output_binary,
        eval_copy,
    ] + LDFLAGS

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    return {
        "success": result.returncode == 0,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "returncode": result.returncode,
    }


def run_evaluation(binary: str, stcoll_dir: str, adv_size: int, timeout: int = 120) -> dict:
    """Run the compiled evaluation binary and parse metrics."""
    cmd = [binary, stcoll_dir, str(adv_size)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "timeout", "metrics": {}}

    if result.returncode != 0:
        return {"success": False, "error": "crash", "stderr": result.stderr, "metrics": {}}

    # Parse metrics from output
    metrics = {}
    output = result.stdout

    # Parse individual matrix results
    pass_count = 0
    fail_count = 0
    for line in output.split("\n"):
        if "PASS" in line:
            pass_count += 1
        elif "FAIL" in line and "res=" in line:
            fail_count += 1

    # Parse OpenEvolve metrics section
    for line in output.split("\n"):
        for key in ["pass_rate", "avg_residual", "avg_ortho_u", "avg_ortho_v",
                    "max_residual", "max_ortho_u", "max_ortho_v",
                    "pass_avg_residual", "pass_avg_ortho_u", "pass_avg_ortho_v",
                    "pass_max_residual", "pass_max_ortho_u", "pass_max_ortho_v",
                    "pass_worst_scaling",
                    "composite_score"]:
            if line.strip().startswith(f"{key}="):
                try:
                    metrics[key] = float(line.strip().split("=")[1])
                except (ValueError, IndexError):
                    pass

    # Parse scaling metrics
    if "pass_worst_scaling" in metrics:
        metrics["max_scaling_ratio"] = metrics["pass_worst_scaling"]

    metrics["pass_count"] = pass_count
    metrics["fail_count"] = fail_count
    metrics["total_count"] = pass_count + fail_count

    return {
        "success": True,
        "metrics": metrics,
        "stdout": output,
    }


def _make_result(score: float, metrics: dict, artifacts: dict = None) -> EvaluationResult:
    """Build an EvaluationResult with score + flat float-only metrics + artifacts.

    OpenEvolve expects metrics to be Dict[str, float]. Non-float values are
    stripped. The 'score' key is always included in metrics.
    """
    flat = {"score": float(score)}
    for k, v in metrics.items():
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            flat[k] = float(v)
    return EvaluationResult(metrics=flat, artifacts=artifacts or {})


def evaluate(program_path: str) -> EvaluationResult:
    """
    Main evaluation function called by OpenEvolve.

    Args:
        program_path: Path to the evolved bidiag_svd.h file

    Returns:
        EvaluationResult with flat metrics dict and artifacts
    """
    artifacts = {}
    metrics = {}

    # Create temp directory for build
    build_dir = tempfile.mkdtemp(prefix="xmr_eval_")
    binary = os.path.join(build_dir, "evaluate")

    try:
        # === Stage 1: Compile check ===
        compile_result = compile_program(program_path, binary)
        if not compile_result["success"]:
            artifacts["compile_error"] = compile_result["stderr"][:2000]
            return _make_result(0.0, FAIL_METRICS, artifacts)
        metrics["compile_success"] = 1

        # === Stage 2: Full evaluation (adv_size=200 + STCollection) ===
        # With 30s cumulative timeout + 2s per-test timeout + early abort in
        # evaluate.cpp, this completes in <5s for correct O(n²) algorithms.
        eval_result = run_evaluation(binary, STCOLL_DIR, 200, timeout=60)
        if not eval_result["success"]:
            artifacts["eval_error"] = eval_result.get("error", "unknown")
            return _make_result(5.0, {**FAIL_METRICS, "compile_success": 1}, artifacts)

        metrics.update(eval_result["metrics"])
        artifacts["eval_output"] = eval_result["stdout"][:10000]

        # === Compute composite score ===
        # Single-stage scoring — matches evaluate.cpp formula exactly
        score = 0.0
        pass_rate = metrics.get("pass_rate", 0)

        # Compile success: 5 points
        score += 5.0

        # Pass rate: up to 50 points
        score += pass_rate * 50.0

        # Accuracy bonus: up to 30 points (use PASS averages, not all-test averages)
        avg_res = metrics.get("pass_avg_residual", 100)
        avg_ortU = metrics.get("pass_avg_ortho_u", 100)
        avg_ortV = metrics.get("pass_avg_ortho_v", 100)
        score += max(0, 5 - min(avg_res, 5)) * 2.0
        score += max(0, 5 - min(avg_ortU, 5)) * 2.0
        score += max(0, 5 - min(avg_ortV, 5)) * 2.0

        # O(n²) HARD GATE: if worst scaling ratio > 5.0x, score capped at 5
        # This prevents gaming by using O(n³) fallbacks (e.g., DBDSQR)
        worst_ratio = metrics.get("pass_worst_scaling", 99)
        if worst_ratio > 5.0:
            score = 5.0  # compilation credit only
        else:
            score += 10.0  # full scaling bonus

        metrics["combined_score"] = score

        # Save variant and log to changelog
        try:
            save_variant_and_log(program_path, score, metrics, artifacts)
        except Exception:
            pass  # Non-fatal — don't fail evaluation because of logging

        return _make_result(score, metrics, artifacts)

    finally:
        # Cleanup
        shutil.rmtree(build_dir, ignore_errors=True)


# ============================================================
# Cascade evaluation stages (for OpenEvolve cascade mode)
# ============================================================

def evaluate_stage1(program_path: str) -> EvaluationResult:
    """Stage 1: Compile + small test (fast filter, ~5 seconds)."""
    build_dir = tempfile.mkdtemp(prefix="xmr_s1_")
    binary = os.path.join(build_dir, "evaluate")
    try:
        compile_result = compile_program(program_path, binary)
        if not compile_result["success"]:
            return _make_result(0.0, FAIL_METRICS,
                                {"compile_error": compile_result["stderr"][:2000]})
        # Run only small adversarial (n=50, no STCollection)
        eval_result = run_evaluation(binary, "", 50, timeout=30)
        if not eval_result["success"]:
            return _make_result(5.0, {**FAIL_METRICS, "compile_success": 1})
        score = 5.0 + eval_result["metrics"].get("pass_rate", 0) * 30.0
        return _make_result(score, {**eval_result["metrics"], "compile_success": 1})
    finally:
        shutil.rmtree(build_dir, ignore_errors=True)


def evaluate_stage2(program_path: str) -> EvaluationResult:
    """Stage 2: STCollection + adversarial n=100 (~30 seconds)."""
    build_dir = tempfile.mkdtemp(prefix="xmr_s2_")
    binary = os.path.join(build_dir, "evaluate")
    try:
        compile_result = compile_program(program_path, binary)
        if not compile_result["success"]:
            return _make_result(0.0, FAIL_METRICS)
        eval_result = run_evaluation(binary, STCOLL_DIR, 100, timeout=120)
        if not eval_result["success"]:
            return _make_result(5.0, {**FAIL_METRICS, "compile_success": 1})
        metrics = eval_result["metrics"]
        score = 5.0 + metrics.get("pass_rate", 0) * 50.0
        avg_res = metrics.get("pass_avg_residual", 100)
        avg_ortU = metrics.get("pass_avg_ortho_u", 100)
        avg_ortV = metrics.get("pass_avg_ortho_v", 100)
        score += max(0, 5 - min(avg_res, 5)) * 2.0
        score += max(0, 5 - min(avg_ortU, 5)) * 2.0
        score += max(0, 5 - min(avg_ortV, 5)) * 2.0
        return _make_result(score, {**metrics, "compile_success": 1},
                            {"eval_output": eval_result["stdout"][:10000]})
    finally:
        shutil.rmtree(build_dir, ignore_errors=True)


def evaluate_stage3(program_path: str) -> EvaluationResult:
    """Stage 3: Full evaluation with large matrices + scaling."""
    return evaluate(program_path)


# ============================================================
# Variant persistence + changelog
# ============================================================

def _next_variant_number() -> int:
    """Get the next variant number by scanning evolved_variants/."""
    os.makedirs(VARIANTS_DIR, exist_ok=True)
    existing = globmod.glob(os.path.join(VARIANTS_DIR, "v*.h"))
    if not existing:
        return 1
    nums = []
    for f in existing:
        base = os.path.basename(f)
        # Format: v{N}_{score}.h
        try:
            nums.append(int(base.split("_")[0][1:]))
        except (ValueError, IndexError):
            pass
    return max(nums, default=0) + 1


def _extract_failing_tests(eval_output: str, max_tests: int = 5) -> list:
    """Extract names of failing tests from evaluation output."""
    fails = []
    for line in eval_output.split("\n"):
        if "FAIL" in line and "res=" in line:
            # Lines look like: "  stemr_killer          n= 200  res= 847.3  ...  FAIL"
            parts = line.strip().split()
            if parts:
                fails.append(parts[0])
    # Deduplicate (same pattern fails at multiple sizes)
    seen = set()
    unique = []
    for f in fails:
        if f not in seen:
            seen.add(f)
            unique.append(f)
    return unique[:max_tests]


def _read_pending_research() -> str:
    """Read and consume the agent's staged RESEARCH.md."""
    pending = os.path.join(VARIANTS_DIR, ".pending_research.md")
    if not os.path.exists(pending):
        return ""
    try:
        with open(pending) as f:
            research = f.read().strip()
        os.remove(pending)
        return research
    except Exception:
        return ""


def _extract_abstract(research: str) -> str:
    """Extract the Abstract section from a RESEARCH.md."""
    lines = research.split("\n")
    capture = False
    result = []
    for line in lines:
        if line.startswith("## ") and "abstract" in line.lower():
            capture = True
            continue
        elif line.startswith("## ") and capture:
            break
        elif line.startswith("# ") and capture:
            break
        elif capture:
            result.append(line)
    return " ".join(result).strip()


def save_variant_and_log(program_path: str, score: float, metrics: dict,
                         artifacts: dict) -> str:
    """Save evolved variant to disk and append to CHANGELOG.md."""
    os.makedirs(VARIANTS_DIR, exist_ok=True)
    num = _next_variant_number()
    variant_name = f"v{num}_{score:.0f}.h"
    variant_path = os.path.join(VARIANTS_DIR, variant_name)

    # Save the variant file
    try:
        shutil.copy2(program_path, variant_path)
    except Exception:
        return ""

    # Read agent's RESEARCH.md (staged by claude_agent_llm.py)
    research = _read_pending_research()

    # Save full research report alongside variant file
    if research:
        research_path = os.path.join(VARIANTS_DIR, f"v{num}_RESEARCH.md")
        try:
            with open(research_path, "w") as f:
                f.write(research)
        except Exception:
            pass

    # Build changelog entry
    pass_rate = metrics.get("pass_rate", 0)
    pass_count = metrics.get("pass_count", 0)
    total_count = metrics.get("total_count", 0)
    avg_res = metrics.get("pass_avg_residual", -1)
    avg_ortU = metrics.get("pass_avg_ortho_u", -1)
    avg_ortV = metrics.get("pass_avg_ortho_v", -1)
    scaling = metrics.get("pass_worst_scaling", -1)
    large_pass = metrics.get("large_pass_rate", -1)

    # Extract failing test names
    eval_output = artifacts.get("eval_output", "") + artifacts.get("large_eval_output", "")
    failing = _extract_failing_tests(eval_output)
    fail_str = ", ".join(failing) if failing else "(none)"
    if len(failing) == 5:
        fail_str += ", ..."

    # Extract abstract one-liner for table (full report in v{N}_RESEARCH.md)
    approach_summary = ""
    if research:
        abstract = _extract_abstract(research)
        if abstract:
            approach_summary = abstract[:80]
            if len(abstract) > 80:
                approach_summary += "..."

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Initialize changelog if it doesn't exist
    if not os.path.exists(CHANGELOG_FILE):
        with open(CHANGELOG_FILE, "w") as f:
            f.write("# Evolution Changelog\n\n")
            f.write("Each entry = one evaluated variant. Variants saved in this directory.\n")
            f.write("Full research reports: see v{N}_RESEARCH.md files.\n\n")
            f.write("| # | Score | Pass | Avg Res | Avg OrtU | Scaling | Approach | Failing Tests | Time |\n")
            f.write("|---|-------|------|---------|----------|---------|----------|---------------|------|\n")

    # Append entry
    with open(CHANGELOG_FILE, "a") as f:
        f.write(
            f"| {num} | {score:.1f} | {pass_count}/{total_count} "
            f"| {avg_res:.2f} | {avg_ortU:.2f} | {scaling:.2f}x "
            f"| {approach_summary or '(no report)'} | {fail_str} | {timestamp} |\n"
        )

    return variant_path


# ============================================================
# Variant persistence + changelog
# ============================================================

def _next_variant_number() -> int:
    """Get the next variant number by scanning evolved_variants/."""
    os.makedirs(VARIANTS_DIR, exist_ok=True)
    existing = globmod.glob(os.path.join(VARIANTS_DIR, "v*.h"))
    if not existing:
        return 1
    nums = []
    for f in existing:
        base = os.path.basename(f)
        # Format: v{N}_{score}.h
        try:
            nums.append(int(base.split("_")[0][1:]))
        except (ValueError, IndexError):
            pass
    return max(nums, default=0) + 1


def _extract_failing_tests(eval_output: str, max_tests: int = 5) -> list:
    """Extract names of failing tests from evaluation output."""
    fails = []
    for line in eval_output.split("\n"):
        if "FAIL" in line and "res=" in line:
            # Lines look like: "  stemr_killer          n= 200  res= 847.3  ...  FAIL"
            parts = line.strip().split()
            if parts:
                fails.append(parts[0])
    # Deduplicate (same pattern fails at multiple sizes)
    seen = set()
    unique = []
    for f in fails:
        if f not in seen:
            seen.add(f)
            unique.append(f)
    return unique[:max_tests]


def _read_pending_research() -> str:
    """Read and consume the agent's staged RESEARCH.md."""
    pending = os.path.join(VARIANTS_DIR, ".pending_research.md")
    if not os.path.exists(pending):
        return ""
    try:
        with open(pending) as f:
            research = f.read().strip()
        os.remove(pending)
        return research
    except Exception:
        return ""


def _extract_abstract(research: str) -> str:
    """Extract the Abstract section from a RESEARCH.md."""
    lines = research.split("\n")
    capture = False
    result = []
    for line in lines:
        if line.startswith("## ") and "abstract" in line.lower():
            capture = True
            continue
        elif line.startswith("## ") and capture:
            break
        elif line.startswith("# ") and capture:
            break
        elif capture:
            result.append(line)
    return " ".join(result).strip()


def save_variant_and_log(program_path: str, score: float, metrics: dict,
                         artifacts: dict) -> str:
    """Save evolved variant to disk and append to CHANGELOG.md."""
    os.makedirs(VARIANTS_DIR, exist_ok=True)
    num = _next_variant_number()
    variant_name = f"v{num}_{score:.0f}.h"
    variant_path = os.path.join(VARIANTS_DIR, variant_name)

    # Save the variant file
    try:
        shutil.copy2(program_path, variant_path)
    except Exception:
        return ""

    # Read agent's RESEARCH.md (staged by claude_agent_llm.py)
    research = _read_pending_research()

    # Save full research report alongside variant file
    if research:
        research_path = os.path.join(VARIANTS_DIR, f"v{num}_RESEARCH.md")
        try:
            with open(research_path, "w") as f:
                f.write(research)
        except Exception:
            pass

    # Build changelog entry
    pass_rate = metrics.get("pass_rate", 0)
    pass_count = metrics.get("pass_count", 0)
    total_count = metrics.get("total_count", 0)
    avg_res = metrics.get("pass_avg_residual", -1)
    avg_ortU = metrics.get("pass_avg_ortho_u", -1)
    avg_ortV = metrics.get("pass_avg_ortho_v", -1)
    scaling = metrics.get("pass_worst_scaling", -1)
    large_pass = metrics.get("large_pass_rate", -1)

    # Extract failing test names
    eval_output = artifacts.get("eval_output", "") + artifacts.get("large_eval_output", "")
    failing = _extract_failing_tests(eval_output)
    fail_str = ", ".join(failing) if failing else "(none)"
    if len(failing) == 5:
        fail_str += ", ..."

    # Extract abstract one-liner for table (full report in v{N}_RESEARCH.md)
    approach_summary = ""
    if research:
        abstract = _extract_abstract(research)
        if abstract:
            approach_summary = abstract[:80]
            if len(abstract) > 80:
                approach_summary += "..."

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Initialize changelog if it doesn't exist
    if not os.path.exists(CHANGELOG_FILE):
        with open(CHANGELOG_FILE, "w") as f:
            f.write("# Evolution Changelog\n\n")
            f.write("Each entry = one evaluated variant. Variants saved in this directory.\n")
            f.write("Full research reports: see v{N}_RESEARCH.md files.\n\n")
            f.write("| # | Score | Pass | Avg Res | Avg OrtU | Scaling | Approach | Failing Tests | Time |\n")
            f.write("|---|-------|------|---------|----------|---------|----------|---------------|------|\n")

    # Append entry
    with open(CHANGELOG_FILE, "a") as f:
        f.write(
            f"| {num} | {score:.1f} | {pass_count}/{total_count} "
            f"| {avg_res:.2f} | {avg_ortU:.2f} | {scaling:.2f}x "
            f"| {approach_summary or '(no report)'} | {fail_str} | {timestamp} |\n"
        )

    return variant_path


# For standalone testing
if __name__ == "__main__":
    import sys
    program = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PROJECT_DIR, "src", "bidiag_svd.h")
    result = evaluate(program)
    print(f"\nScore: {result['score']:.2f}")
    print("Metrics:")
    for k, v in sorted(result["metrics"].items()):
        print(f"  {k}: {v}")
