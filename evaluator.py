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

# Paths
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.join(PROJECT_DIR, "lib")
STCOLL_DIR = os.path.join(PROJECT_DIR, "..", "bidiag-algo", "MRRR Resources", "STCollection", "DATA")
EVAL_SRC = os.path.join(PROJECT_DIR, "src", "evaluate.cpp")

# Compiler settings — pure C/C++ (no Fortran, no Accelerate)
CXX = "g++"
CLAPACK_DIR = os.path.join(PROJECT_DIR, "src", "clapack")
CXXFLAGS = ["-std=c++17", "-O2", f"-I{CLAPACK_DIR}"]
LDFLAGS = [f"-L{LIB_DIR}", "-lxmr_c", "-lm"]


def compile_program(program_path: str, output_binary: str) -> dict:
    """Compile the evolved bidiag_svd.h and evaluate.cpp into a binary."""
    build_dir = os.path.dirname(output_binary)

    # Copy evaluate.cpp to build dir so #include "bidiag_svd.h" finds the
    # evolved version (GCC searches the directory of the including file first,
    # so compiling the original src/evaluate.cpp always finds src/bidiag_svd.h)
    eval_copy = os.path.join(build_dir, "evaluate.cpp")
    shutil.copy2(EVAL_SRC, eval_copy)

    cmd = [CXX] + CXXFLAGS + [
        f"-I{os.path.dirname(program_path)}",  # Evolved file's dir
        "-o", output_binary,
        eval_copy,  # Compile the copy, not the original
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


def evaluate(program_path: str) -> dict:
    """
    Main evaluation function called by OpenEvolve.

    Args:
        program_path: Path to the evolved bidiag_svd.h file

    Returns:
        dict with:
          - "score": float (0-100, higher is better)
          - "metrics": dict of named metrics
          - "artifacts": dict of string artifacts (evaluation output, etc.)
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
            return {
                "score": 0.0,
                "metrics": {"compile_success": 0},
                "artifacts": artifacts,
            }
        metrics["compile_success"] = 1

        # === Stage 2: STCollection + small adversarial ===
        eval_result = run_evaluation(binary, STCOLL_DIR, 100, timeout=120)
        if not eval_result["success"]:
            artifacts["eval_error"] = eval_result.get("error", "unknown")
            return {
                "score": 5.0,  # Gets some credit for compiling
                "metrics": metrics,
                "artifacts": artifacts,
            }

        metrics.update(eval_result["metrics"])
        artifacts["eval_output"] = eval_result["stdout"][:10000]

        # === Stage 3: Larger adversarial + scaling (if stage 2 passed well) ===
        pass_rate = metrics.get("pass_rate", 0)
        if pass_rate >= 0.6:
            eval_large = run_evaluation(binary, STCOLL_DIR, 200, timeout=600)
            if eval_large["success"]:
                large_metrics = eval_large["metrics"]
                metrics["large_pass_rate"] = large_metrics.get("pass_rate", 0)
                metrics["large_max_ratio"] = large_metrics.get("max_scaling_ratio", 99)
                artifacts["large_eval_output"] = eval_large["stdout"][:5000]

        # === Compute composite score ===
        # Match evaluate.cpp's formula but use pass_avg metrics (not all-test avg)
        score = 0.0

        # Compile success: 5 points
        score += 5.0

        # Pass rate: up to 50 points
        score += pass_rate * 50.0

        # Accuracy bonus: up to 15 points (use PASS averages, not all-test averages)
        avg_res = metrics.get("pass_avg_residual", 100)
        avg_ortU = metrics.get("pass_avg_ortho_u", 100)
        avg_ortV = metrics.get("pass_avg_ortho_v", 100)
        score += max(0, 5 - min(avg_res, 5)) * 2.0
        score += max(0, 5 - min(avg_ortU, 5)) * 2.0
        score += max(0, 5 - min(avg_ortV, 5)) * 2.0

        # Scaling bonus: up to 10 points (from evaluate.cpp formula)
        worst_ratio = metrics.get("pass_worst_scaling", 99)
        if worst_ratio <= 5.0:
            score += 10.0
        elif worst_ratio <= 10.0:
            score += 10.0 * (10.0 - worst_ratio) / 5.0

        # Large adversarial pass rate bonus: up to 10 points
        if "large_pass_rate" in metrics:
            score += metrics["large_pass_rate"] * 10.0

        # Large adversarial scaling bonus: up to 5 points
        if "large_max_ratio" in metrics:
            ratio = metrics["large_max_ratio"]
            if ratio <= 5.0:
                score += 5.0
            elif ratio <= 10.0:
                score += 5.0 * (10.0 - ratio) / 5.0

        metrics["combined_score"] = score

        return {
            "score": score,
            "metrics": metrics,
            "artifacts": artifacts,
        }

    finally:
        # Cleanup
        shutil.rmtree(build_dir, ignore_errors=True)


# ============================================================
# Cascade evaluation stages (for OpenEvolve cascade mode)
# ============================================================

def evaluate_stage1(program_path: str) -> dict:
    """Stage 1: Compile + small test (fast filter, ~5 seconds)."""
    build_dir = tempfile.mkdtemp(prefix="xmr_s1_")
    binary = os.path.join(build_dir, "evaluate")
    try:
        compile_result = compile_program(program_path, binary)
        if not compile_result["success"]:
            return {
                "score": 0.0,
                "metrics": {"compile_success": 0},
                "artifacts": {"compile_error": compile_result["stderr"][:2000]},
            }
        # Run only small adversarial (n=50, no STCollection)
        eval_result = run_evaluation(binary, "", 50, timeout=30)
        if not eval_result["success"]:
            return {"score": 5.0, "metrics": {"compile_success": 1}}
        score = 5.0 + eval_result["metrics"].get("pass_rate", 0) * 30.0
        return {
            "score": score,
            "metrics": {**eval_result["metrics"], "compile_success": 1},
        }
    finally:
        shutil.rmtree(build_dir, ignore_errors=True)


def evaluate_stage2(program_path: str) -> dict:
    """Stage 2: STCollection + adversarial n=100 (~30 seconds)."""
    build_dir = tempfile.mkdtemp(prefix="xmr_s2_")
    binary = os.path.join(build_dir, "evaluate")
    try:
        compile_result = compile_program(program_path, binary)
        if not compile_result["success"]:
            return {"score": 0.0, "metrics": {"compile_success": 0}}
        eval_result = run_evaluation(binary, STCOLL_DIR, 100, timeout=120)
        if not eval_result["success"]:
            return {"score": 5.0, "metrics": {"compile_success": 1}}
        metrics = eval_result["metrics"]
        score = 5.0 + metrics.get("pass_rate", 0) * 50.0
        # Use pass averages for accuracy bonus
        avg_res = metrics.get("pass_avg_residual", 100)
        avg_ortU = metrics.get("pass_avg_ortho_u", 100)
        avg_ortV = metrics.get("pass_avg_ortho_v", 100)
        score += max(0, 5 - min(avg_res, 5)) * 2.0
        score += max(0, 5 - min(avg_ortU, 5)) * 2.0
        score += max(0, 5 - min(avg_ortV, 5)) * 2.0
        return {
            "score": score,
            "metrics": {**metrics, "compile_success": 1},
            "artifacts": {"eval_output": eval_result["stdout"][:10000]},
        }
    finally:
        shutil.rmtree(build_dir, ignore_errors=True)


def evaluate_stage3(program_path: str) -> dict:
    """Stage 3: Full evaluation with large matrices + scaling."""
    return evaluate(program_path)


# For standalone testing
if __name__ == "__main__":
    import sys
    program = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PROJECT_DIR, "src", "bidiag_svd.h")
    result = evaluate(program)
    print(f"\nScore: {result['score']:.2f}")
    print("Metrics:")
    for k, v in sorted(result["metrics"].items()):
        print(f"  {k}: {v}")
