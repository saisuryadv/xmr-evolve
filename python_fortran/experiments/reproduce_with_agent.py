#!/usr/bin/env python3
"""
Validate the prompt sequence in docs/REPRODUCE.md by running Claude Code
agents against a clean checkpoint and feeding the prompts in order.

Usage:
    # Smoke test on prompt 1 only, against a worktree of the seed commit
    python3 experiments/reproduce_with_agent.py --checkpoint B --only 1

    # Full Checkpoint B run (all 6 prompts)
    python3 experiments/reproduce_with_agent.py --checkpoint B

    # Use an existing worktree (skip the checkout step)
    python3 experiments/reproduce_with_agent.py --worktree /tmp/repro_xyz --only 1

For each prompt, the harness:
  1. Records the current evaluate.py score (before).
  2. Spawns `claude --print --dangerously-skip-permissions --max-turns 60`
     with the prompt as input, cwd = worktree/python_fortran.
  3. Records the score after (re-running evaluate.py).
  4. Logs to experiments/reproduce_results.json.

Each prompt runs a *fresh* Claude session (no context carry-over between
prompts) so each prompt is self-contained.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent  # python_fortran/
REPO_ROOT = ROOT.parent  # xmr-evolve/
SEED_COMMIT = "c6d73b2"  # Checkpoint B starting point
DEFAULT_OUT = ROOT / "experiments" / "reproduce_results.json"

CHECKPOINT_B_PROMPTS = [
    {
        "id": 1,
        "title": "singleton-block bypass",
        "expected_after": (305, 320),
        "prompt": """Run `python3 evaluate.py` and observe the current score (you should see TOTAL: ~300/379). Several adversarial patterns time out — particularly chkbd_*, saw_tooth, and step_function at n >= 200. The root cause is in mr3_gk.py:mr3_tgk: when the bidiagonal splitter produces a singleton block (k=1), the code still calls XMR on a trivial 2x2 T_GK matrix. XMR returns degenerate eigenvectors with ||v||=1, ||u||=0 (only even rows populated). The post-processing then runs O(n^3) Gram-Schmidt on hundreds of such columns, causing timeout.

Add a fast path in mr3_tgk that handles singleton blocks directly, bypassing the XMR call:
  - sigma = |d[i]|
  - both u and v components are 1/sqrt(2) at the corresponding row
  - apply sign matrices D1, D2 (constructed earlier in bidiag_svd) to recover correct signs

Do NOT modify dlaxre_gk.f or any Fortran file. Re-run evaluate.py. The chkbd patterns must no longer time out. Score should jump to roughly 305-320/379 and SCORE much higher than 5.00.""",
    },
    {
        "id": 2,
        "title": "Bv recovery threshold + GS vectorization",
        "expected_after": (330, 360),
        "prompt": """Some matrices show one-sided eigenvector quality issues — only u or only v has the correct norm. Look at gl_gradm@200 (ortU about 33 currently) and chkbd@200.

Two fixes in mr3_gk.py:
(a) The Bv recovery threshold is hardcoded — make it scale with sigma_max so it adapts to matrix dynamic range.
(b) The Gram-Schmidt completion runs one column at a time in a Python loop. Vectorize it: build a matrix of all candidate vectors and do the projection in a single numpy operation so each MGS pass is one matmul instead of n matvecs.

Verify on gl_gradm and chkbd — ortU should drop below threshold. Don't worry about other failures yet. Re-run evaluate.py.""",
    },
    {
        "id": 3,
        "title": "T_GK sub-splitting + negative-e fix",
        "expected_after": (335, 360),
        "prompt": """Two refinements in mr3_gk.py:bidiag_svd, near where T_GK is built:

(a) Zero diagonal entries: when the bidiagonal has d[i]=0, the corresponding rows of T_GK have zero off-diagonals on both sides — there's no eigenvalue coupling across that point. Sub-split T_GK at such positions so MR^3 doesn't waste effort on decoupled blocks.

(b) Negative e[i]: the perfect-shuffle T_GK construction gives negative off-diagonals when e[i] < 0, but XMR requires positive off-diagonals. Apply a diagonal-similarity sign normalization before forming T_GK (track the sign flips and re-apply them when extracting U, V).

These don't change the headline score by much but eliminate edge-case failures. Verify the score doesn't regress. Re-run evaluate.py.""",
    },
    {
        "id": 4,
        "title": "zero-shift QR deflation",
        "expected_after": (350, 378),
        "prompt": """Demmel-Kahan 1990 (Sec. 3) shows that one zero-shift QR sweep on a bidiagonal is equivalent to one step of inverse iteration on B^T B. When sigma_min = 0, BOTH d[-1] AND e[-1] are driven to zero in a single sweep — clean detection of zero singular values.

Implement this as preprocessing in mr3_gk.py:bidiag_svd, applied to each unsplit block of size >= 2 BEFORE the MR^3 dispatch:

1. Implement zero_shift_qr_sweep(d, e) using LAPACK dlartg (from scipy.linalg.lapack):
   - Walk i = 0..n-2, generate a Givens rotation per step, update d, e in place.
   - Return (right_rots, left_rots) — Givens-rotation lists.

2. After the sweep, check the clean-split criterion:
     |d[-1]| < n * eps  AND  |e[-1]| < n * eps
   If TRUE: zero singular value detected. Solve the (k-1) sub-problem recursively, embed the solution into k×k with sigma=0 in the bottom-right, and apply the saved Givens rotations to recover the original basis.
   If FALSE: discard the swept (d, e), use the original block for MR^3 normally.

CRITICAL: the clean-split check must require BOTH d[-1] AND e[-1] to be small. Don't deflate when only d[-1] is small — that severs matrix coupling and produces catastrophic residuals.

This fixes saw_tooth, step_function, gl_wilkp, gl_clement, and any matrix with a structural zero diagonal. Re-run evaluate.py; score should jump significantly.""",
    },
    {
        "id": 5,
        "title": "dlaxrb_clssfy AVGTHRESH=0 fix",
        "expected_after": (370, 379),
        "prompt": """A subtle bug in upstream xmr_src/dlaxrb_clssfy.f. At deep tree nodes, AVGTHRESH = 2*AVGTOL collapses to 0 because AVGTOL is scaled by 2^(-DEPTH) and underflows. The classifier test is

   IF( SPREAD .GE. MIN(AVGTHRESH, ABSMAX*GAPTHRESH) ) gap is full

— with AVGTHRESH=0 this becomes SPREAD >= 0, trivially true, so EVERY adjacent eigenvalue pair is labeled a full gap and the rep tree degenerates into singletons. Tight-cluster matrices fail catastrophically.

Steps:
1. Create python_fortran/dlaxrb_clssfy_fix.f as a copy of xmr_src/dlaxrb_clssfy.f.
2. Add a logical local variable GAPOK to the declarations.
3. At all THREE test sites in the file (search for "MIN( AVGTHRESH, ABSMAX*GAPTHRESH )"), wrap each with a case-split:

       IF( AVGTHRESH .GT. ZERO )THEN
          GAPOK = SPREAD .GE. MIN(AVGTHRESH, ABSMAX*GAPTHRESH)
       ELSE
          GAPOK = SPREAD .GE. ABSMAX*GAPTHRESH
       ENDIF
       IF( GAPOK )THEN
          ! existing 'gap is full' code path
       ELSE
          ! existing 'no full gap' code path
       ENDIF

4. Modify build.sh to compile dlaxrb_clssfy_fix.f and link its .o in place of the upstream dlaxrb_clssfy.o (rename or override). Rebuild libxmr.so via `bash build.sh`.

5. Re-run evaluate.py. demmel_S1pe_k4, three_clusters, pd_T0, gl_wilkw, ST_B_20_graded should all flip from FAIL to PASS.""",
    },
    {
        "id": 6,
        "title": "adaptive GAPTOL for tiny blocks",
        "expected_after": (379, 379),
        "prompt": """One last failure remains: two_clusters@10 (n=10). The cluster's relative gap is right at the GAPTOL boundary, but at small n the hard-coded GAPTOL = 1e-3 is too coarse — small absolute differences look 'large enough' even when they shouldn't.

In mr3_gk.py:classify (the per-node gap classifier), replace every use of the literal `GAPTOL` with `max(GAPTOL, 0.02 / n)` where n is the current block size. This is adaptive: large n keeps the standard 1e-3; tiny n tightens it.

Verify: python3 evaluate.py should now report TOTAL: 379/379.""",
    },
]


def run_evaluate(workdir: Path, timeout: int = 600):
    """Run python3 evaluate.py in workdir; return (passed, total, score, err).
    Captures the LAST `TOTAL: X/Y` line (the full-suite total, not any
    intermediate per-size subtotal) and the SCORE line if present.
    """
    try:
        r = subprocess.run(
            ["python3", "evaluate.py"],
            cwd=str(workdir),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return None, None, None, "timeout"
    out = r.stdout + "\n" + r.stderr
    # Last "TOTAL: X/Y passed" wins (final summary)
    matches = re.findall(r"TOTAL:\s*(\d+)\s*/\s*(\d+)", out)
    score = None
    sm = re.search(r"SCORE:\s*([\-+\d.]+)", out)
    if sm:
        try:
            score = float(sm.group(1))
        except ValueError:
            pass
    if matches:
        passed, total = matches[-1]
        return int(passed), int(total), score, None
    return None, None, score, "no_total_line"


def make_worktree(commit: str, dest: Path, overlay_eval_from_head: bool = True):
    """Create a fresh git worktree at the given commit.

    If overlay_eval_from_head=True, copy evaluate.py and full_eval.py from the
    current project tree (HEAD) into the worktree's python_fortran/ directory.
    The seed-commit evaluator has a regex bug that misses 2 of 90 patterns
    (suite reports 371 instead of 379); the modern evaluator reports the full
    379 against the same seed-commit code. We want the recipe measured against
    the 379-test scoreboard from the very first run.
    """
    if dest.exists():
        shutil.rmtree(str(dest))
    dest.parent.mkdir(parents=True, exist_ok=True)
    r = subprocess.run(
        ["git", "worktree", "add", "--detach", str(dest), commit],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        raise RuntimeError(f"git worktree add failed: {r.stderr}")

    if overlay_eval_from_head:
        wf = dest / "python_fortran"
        head_pf = REPO_ROOT / "python_fortran"
        for fname in ("evaluate.py", "full_eval.py"):
            src = head_pf / fname
            dst = wf / fname
            if src.exists() and wf.exists():
                shutil.copy2(str(src), str(dst))
                print(f"  overlaid {fname} from HEAD into worktree")
    return dest


def remove_worktree(dest: Path):
    """Best-effort worktree cleanup."""
    subprocess.run(
        ["git", "worktree", "remove", "--force", str(dest)],
        cwd=str(REPO_ROOT),
        capture_output=True,
    )


def run_agent(prompt: str, workdir: Path, max_turns: int = 60,
              timeout: int = 1800):
    """Spawn `claude --print` with the prompt; return (stdout_text, returncode, wall_s)."""
    cmd = [
        "claude",
        "--print",
        "--dangerously-skip-permissions",
        "--max-turns", str(max_turns),
        "--add-dir", str(workdir),
    ]
    t0 = time.perf_counter()
    try:
        r = subprocess.run(
            cmd,
            input=prompt,
            cwd=str(workdir),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as ex:
        return ex.stdout or "", 124, time.perf_counter() - t0
    return r.stdout, r.returncode, time.perf_counter() - t0


def list_changed_files(workdir: Path):
    r = subprocess.run(
        ["git", "diff", "--name-only", "HEAD"],
        cwd=str(workdir),
        capture_output=True,
        text=True,
    )
    return [ln for ln in r.stdout.splitlines() if ln.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", choices=["B"], default="B",
                    help="B = start from c6d73b2 (300/379)")
    ap.add_argument("--only", type=int, default=0,
                    help="Run only prompt N (1-6), 0 = all")
    ap.add_argument("--worktree", type=str, default=None,
                    help="Use existing worktree path; skip checkout step")
    ap.add_argument("--max-turns", type=int, default=60)
    ap.add_argument("--timeout", type=int, default=1800,
                    help="Per-prompt agent timeout (seconds)")
    ap.add_argument("--out", type=str, default=str(DEFAULT_OUT))
    ap.add_argument("--keep-worktree", action="store_true",
                    help="Don't delete the worktree at the end")
    args = ap.parse_args()

    prompts = CHECKPOINT_B_PROMPTS
    if args.only:
        prompts = [p for p in prompts if p["id"] == args.only]
        if not prompts:
            sys.exit(f"no such prompt id: {args.only}")

    # Set up worktree
    if args.worktree:
        worktree = Path(args.worktree)
        if not worktree.exists():
            sys.exit(f"worktree {worktree} does not exist")
        print(f"Using existing worktree: {worktree}")
    else:
        ts = int(time.time())
        worktree = Path(f"/tmp/repro_{ts}")
        print(f"Creating worktree at {worktree} from commit {SEED_COMMIT}")
        make_worktree(SEED_COMMIT, worktree)

    workdir = worktree / "python_fortran"
    if not workdir.exists():
        sys.exit(f"{workdir} not found in worktree")

    # Build libxmr.so in the worktree (so evaluate.py can run)
    print(f"\n=== Building libxmr.so in {workdir} ===")
    r = subprocess.run(["bash", "build.sh"], cwd=str(workdir),
                       capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        print(f"build.sh failed: {r.stderr}")
        if not args.keep_worktree and not args.worktree:
            remove_worktree(worktree)
        sys.exit(2)

    # Initial baseline
    print("\n=== Initial evaluate.py ===")
    p0, t0_total, sc0, err0 = run_evaluate(workdir)
    print(f"Score: {p0}/{t0_total}  evaluate.py SCORE={sc0}  (err={err0})")

    results = {
        "checkpoint": args.checkpoint,
        "seed_commit": SEED_COMMIT,
        "worktree": str(worktree),
        "initial_score": {"passed": p0, "total": t0_total,
                           "evaluate_score": sc0, "err": err0},
        "steps": [],
    }

    for spec in prompts:
        pid = spec["id"]
        title = spec["title"]
        prompt = spec["prompt"]
        print(f"\n=== Prompt {pid}: {title} ===")
        score_before = run_evaluate(workdir)
        print(f"Score before: {score_before[0]}/{score_before[1]}  "
              f"SCORE={score_before[2]}")

        agent_out, rc, wall = run_agent(prompt, workdir,
                                         max_turns=args.max_turns,
                                         timeout=args.timeout)
        print(f"Agent finished: rc={rc}  wall={wall:.1f}s "
              f"out_chars={len(agent_out)}")

        score_after = run_evaluate(workdir)
        changed = list_changed_files(worktree)
        print(f"Score after:  {score_after[0]}/{score_after[1]}  "
              f"SCORE={score_after[2]}")
        print(f"Files changed: {len(changed)}")

        step = {
            "id": pid,
            "title": title,
            "prompt": prompt,
            "score_before": {"passed": score_before[0], "total": score_before[1],
                              "evaluate_score": score_before[2]},
            "score_after":  {"passed": score_after[0],  "total": score_after[1],
                              "evaluate_score": score_after[2]},
            "agent_returncode": rc,
            "agent_wall_s": wall,
            "agent_stdout_tail": agent_out[-2000:] if agent_out else "",
            "files_changed": changed,
            "expected_after_range": list(spec["expected_after"]),
        }
        results["steps"].append(step)

        # Persist after each prompt so we don't lose progress on a crash
        with open(args.out, "w") as f:
            json.dump(results, f, indent=1, default=str)

    if not args.keep_worktree and not args.worktree:
        print(f"\n=== Cleaning up worktree {worktree} ===")
        remove_worktree(worktree)

    print(f"\n=== Summary ===")
    print(f"{'#':>2}  {'before':>10}  {'after':>10}  {'sc_before':>9}  {'sc_after':>8}  {'wall':>7}  title")
    for s in results["steps"]:
        b = f"{s['score_before']['passed']}/{s['score_before']['total']}"
        a = f"{s['score_after']['passed']}/{s['score_after']['total']}"
        scb = f"{s['score_before']['evaluate_score']}"
        sca = f"{s['score_after']['evaluate_score']}"
        print(f"{s['id']:>2}  {b:>10}  {a:>10}  {scb:>9}  {sca:>8}  "
              f"{s['agent_wall_s']:>6.0f}s  {s['title']}")
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
