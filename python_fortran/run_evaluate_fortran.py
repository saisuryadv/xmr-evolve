#!/usr/bin/env python3
"""
Run the Fortran-only MR3-GK driver through the same evaluation harnesses the
README documents for the Python+Fortran code:

  python3 run_evaluate_fortran.py                     # all four README suites
  python3 run_evaluate_fortran.py --suite evaluate    # just evaluate.py
  python3 run_evaluate_fortran.py --suite dense       # test_dense_to_bidiag.py
  python3 run_evaluate_fortran.py --suite glued       # test_glued_synth.py
  python3 run_evaluate_fortran.py --suite synth       # test_synth_pract.py

For each suite this script:
  1. Imports the suite's Python module.
  2. Replaces its `bidiag_svd` symbol with a thin wrapper that calls
     `mr3gk_fortran/mr3gk_run` via subprocess + binary I/O.
  3. Calls the suite's existing `main()` so the printed PASS/FAIL output is
     IDENTICAL in format to the README invocation — same thresholds, same
     test ordering, same per-test format string.

This means the four sub-runs each produce the same console output the README
references for Python+Fortran, but the SVD work is done by the pure-Fortran
driver. Their final tallies are then summarised at the end.
"""
import argparse
import importlib
import io
import os
import struct
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
MR3GK_RUN = os.path.join(ROOT, "mr3gk_fortran", "mr3gk_run")


# ---------- Fortran-driver shim ----------

def _write_bin(path, d, e):
    n = len(d)
    with open(path, "wb") as f:
        f.write(struct.pack("i", n))
        f.write(d.astype(np.float64).tobytes())
        f.write(e.astype(np.float64).tobytes())


def _read_bin(path, n):
    with open(path, "rb") as f:
        info = struct.unpack("i", f.read(4))[0]
        sigma = np.frombuffer(f.read(8 * n), dtype=np.float64).copy()
        U = np.frombuffer(f.read(8 * n * n), dtype=np.float64).reshape(n, n, order="F").copy()
        V = np.frombuffer(f.read(8 * n * n), dtype=np.float64).reshape(n, n, order="F").copy()
    return info, sigma, U, V


def fortran_bidiag_svd(d, e):
    """Drop-in replacement for mr3_gk.bidiag_svd that calls the Fortran driver.
    Returns (sigma, U, V, info) — same signature."""
    n = len(d)
    if n == 0:
        return (np.zeros(0), np.zeros((0, 0)), np.zeros((0, 0)), 0)
    in_path = "/tmp/mr3gk_eval_in.bin"
    out_path = "/tmp/mr3gk_eval_out.bin"
    _write_bin(in_path, np.asarray(d, dtype=np.float64),
               np.asarray(e, dtype=np.float64))
    r = subprocess.run([MR3GK_RUN, in_path, out_path],
                       capture_output=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError(
            f"mr3gk_run failed (rc={r.returncode}): "
            f"{r.stderr.decode(errors='ignore')[:200]}")
    info, sigma, U, V = _read_bin(out_path, n)
    return sigma, U, V, int(info)


# ---------- per-suite invocation ----------

SUITES = ["evaluate", "dense", "glued", "synth"]


def patch_bidiag_svd_in(module_name):
    """Import `module_name`; replace its bidiag_svd attribute with the Fortran shim.

    Suppresses the import-time test loop output of full_eval/test_glued_synth
    (which run their own quick suite on import) so we don't print twice.
    """
    saved = sys.stdout; sys.stdout = io.StringIO()
    try:
        mod = importlib.import_module(module_name)
    finally:
        sys.stdout = saved
    # Replace the symbol the suite uses for the SVD work.
    if hasattr(mod, "bidiag_svd"):
        mod.bidiag_svd = fortran_bidiag_svd
    return mod


def run_suite(suite, extra_argv):
    print("\n" + "=" * 78)
    print(f"=== Suite: {suite}  (Fortran driver via mr3gk_run) ===")
    print("=" * 78)
    if not os.path.exists(MR3GK_RUN):
        print(f"ERROR: {MR3GK_RUN} not built. cd mr3gk_fortran && bash build.sh")
        return None

    # Stash original argv so each suite sees only its own flags.
    saved_argv = sys.argv

    if suite == "evaluate":
        # README: `python3 evaluate.py`  /  `python3 evaluate.py --medium`
        # `evaluate.py` imports `from mr3_gk import bidiag_svd`; we patch
        # that name so the Fortran shim is what actually runs.
        mod_mr3 = patch_bidiag_svd_in("mr3_gk")  # ensure module loaded
        mod = patch_bidiag_svd_in("evaluate")
        # Re-bind: evaluate.py imports `bidiag_svd` *into its own namespace*,
        # so patching mr3_gk.bidiag_svd alone is not enough. Patch the alias too.
        mod.bidiag_svd = fortran_bidiag_svd
        sys.argv = ["evaluate.py"] + list(extra_argv)
        try:
            t0 = time.perf_counter()
            mod.main()
            wall = time.perf_counter() - t0
            print(f"\n[suite=evaluate]  wall: {wall:.1f}s")
        finally:
            sys.argv = saved_argv
        return wall

    if suite == "dense":
        # README: `python3 test_dense_to_bidiag.py [--quick|--part1|--part2]`
        mod = patch_bidiag_svd_in("test_dense_to_bidiag")
        mod.bidiag_svd = fortran_bidiag_svd
        sys.argv = ["test_dense_to_bidiag.py"] + list(extra_argv)
        try:
            t0 = time.perf_counter()
            mod.main()
            wall = time.perf_counter() - t0
            print(f"\n[suite=dense]  wall: {wall:.1f}s")
        finally:
            sys.argv = saved_argv
        return wall

    if suite == "glued":
        # README: `python3 test_glued_synth.py [--quick|--no-large]`
        mod = patch_bidiag_svd_in("test_glued_synth")
        mod.bidiag_svd = fortran_bidiag_svd
        sys.argv = ["test_glued_synth.py"] + list(extra_argv)
        try:
            t0 = time.perf_counter()
            mod.main()
            wall = time.perf_counter() - t0
            print(f"\n[suite=glued]  wall: {wall:.1f}s")
        finally:
            sys.argv = saved_argv
        return wall

    if suite == "synth":
        # README: `python3 test_synth_pract.py [--quick|--medium|--paper]`
        mod = patch_bidiag_svd_in("test_synth_pract")
        mod.bidiag_svd = fortran_bidiag_svd
        sys.argv = ["test_synth_pract.py"] + list(extra_argv)
        try:
            t0 = time.perf_counter()
            mod.main()
            wall = time.perf_counter() - t0
            print(f"\n[suite=synth]  wall: {wall:.1f}s")
        finally:
            sys.argv = saved_argv
        return wall

    raise ValueError(f"unknown suite: {suite}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", choices=SUITES + ["all"], default="all")
    ap.add_argument("rest", nargs=argparse.REMAINDER,
                    help="extra args forwarded to the inner suite (after `--`)")
    args = ap.parse_args()
    extra = args.rest
    if extra and extra[0] == "--":
        extra = extra[1:]
    suites = SUITES if args.suite == "all" else [args.suite]
    walls = {}
    for s in suites:
        try:
            walls[s] = run_suite(s, extra)
        except Exception as ex:
            print(f"\n[suite={s}]  ERROR {type(ex).__name__}: {ex}")
            walls[s] = None

    print("\n" + "=" * 78)
    print("=== overall wall times (Fortran driver) ===")
    for s, w in walls.items():
        print(f"  {s:<10}  {w if w is None else f'{w:.1f}s'}")


if __name__ == "__main__":
    main()
