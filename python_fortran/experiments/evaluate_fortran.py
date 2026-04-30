#!/usr/bin/env python3
"""
Run the pure-Fortran MR3-GK driver (mr3gk_fortran/mr3gk_run) on every
test bidiagonal in our suite and apply the same evaluate.py thresholds:

  res  ≤ 7.0 nε
  ortU ≤ 5.0 nε
  ortV ≤ 5.0 nε

Then compare per-test pass/fail to the Python+Fortran baseline.

Suites covered:
  * evaluate.py 379-spec corpus (90 adv patterns × 4 sizes + 19 STCollection)
  * test_dense_to_bidiag.py   (~88 specs)
  * test_glued_synth.py       (~344 specs)
  * test_synth_pract.py quick (subset; full sweep is 18K)

Usage:
  python3 experiments/evaluate_fortran.py                  # all suites
  python3 experiments/evaluate_fortran.py --suite evaluate # 379 only
  python3 experiments/evaluate_fortran.py --quick          # n ≤ 100 only
"""
import os
import sys
import io
import json
import glob
import struct
import argparse
import subprocess
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

EPS = 2.2204460492503131e-16
RES_THRESH = 7.0
ORTHO_THRESH = 5.0

MR3GK_RUN = os.path.join(ROOT, "mr3gk_fortran", "mr3gk_run")


def silent_imports():
    saved = sys.stdout; sys.stdout = io.StringIO()
    try:
        from full_eval import make, adv_names
        from evaluate import load_stcoll
    finally:
        sys.stdout = saved
    return make, adv_names, load_stcoll


def silent_extended_imports():
    """Pull in spec generators from test_dense_to_bidiag, test_glued_synth, test_synth_pract."""
    saved = sys.stdout; sys.stdout = io.StringIO()
    out = {}
    try:
        # Each test_*.py runs its main suite at import time; we suppress and
        # then look at module-level helpers if exposed.
        try:
            import test_dense_to_bidiag as _t1
            out["dense"] = _t1
        except Exception:
            out["dense"] = None
        try:
            import test_glued_synth as _t2
            out["glued"] = _t2
        except Exception:
            out["glued"] = None
        try:
            import test_synth_pract as _t3
            out["synth"] = _t3
        except Exception:
            out["synth"] = None
    finally:
        sys.stdout = saved
    return out


def write_bin(path, d, e):
    n = len(d)
    with open(path, "wb") as f:
        f.write(struct.pack("i", n))
        f.write(d.astype(np.float64).tobytes())
        f.write(e.astype(np.float64).tobytes())


def read_bin(path, n):
    with open(path, "rb") as f:
        info = struct.unpack("i", f.read(4))[0]
        sigma = np.frombuffer(f.read(8 * n), dtype=np.float64).copy()
        U = np.frombuffer(f.read(8 * n * n), dtype=np.float64).reshape(n, n, order="F").copy()
        V = np.frombuffer(f.read(8 * n * n), dtype=np.float64).reshape(n, n, order="F").copy()
    return info, sigma, U, V


def metrics(d, e, sigma, U, V):
    n = len(d)
    B = np.diag(d) + np.diag(e, 1)
    bnorm = abs(d[-1]) if n > 0 else 0.0
    for i in range(n - 1):
        bnorm = max(bnorm, abs(d[i]) + abs(e[i]))
    if bnorm == 0.0:
        bnorm = 1.0
    res = float(np.max(np.abs(B - U @ np.diag(sigma) @ V.T)) / (bnorm * n * EPS))
    ou = float(np.max(np.abs(U.T @ U - np.eye(n))) / (n * EPS))
    ov = float(np.max(np.abs(V.T @ V - np.eye(n))) / (n * EPS))
    ok = (res <= RES_THRESH) and (ou <= ORTHO_THRESH) and (ov <= ORTHO_THRESH)
    return ok, res, ou, ov


def run_fortran_one(d, e, timeout=120):
    n = len(d)
    in_path = "/tmp/mr3gk_eval_in.bin"
    out_path = "/tmp/mr3gk_eval_out.bin"
    write_bin(in_path, d, e)
    t0 = time.perf_counter()
    try:
        r = subprocess.run([MR3GK_RUN, in_path, out_path],
                           capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "dt": float("nan")}
    dt = time.perf_counter() - t0
    if r.returncode != 0:
        return {"status": "crash", "rc": r.returncode, "dt": dt,
                "stderr_tail": r.stderr.decode(errors="ignore").strip().splitlines()[-3:]}
    try:
        info, sigma, U, V = read_bin(out_path, n)
    except Exception as ex:
        return {"status": "noresult", "err": str(ex), "dt": dt}
    ok, res, ou, ov = metrics(d, e, sigma, U, V)
    return {"status": "ok", "ok": ok, "info": int(info),
            "res": res, "ortU": ou, "ortV": ov, "dt": dt}


# -------- spec enumeration --------

def specs_evaluate(quick):
    make, adv_names, load_stcoll = silent_imports()
    sizes = [10, 100] if quick else [10, 100, 200, 400]
    out = []
    for sz in sizes:
        for pat in adv_names:
            try:
                d, e = make(pat, sz)
            except Exception:
                continue
            out.append(("evaluate", f"{pat}@{sz}", d, e))
    if not quick:
        for f in sorted(glob.glob(os.path.join(ROOT, "stcollection", "B_*.dat"))):
            name, d, e = load_stcoll(f)
            out.append(("evaluate", f"ST_{name}", d, e))
    return out


def specs_dense(quick):
    """test_dense_to_bidiag suites:
       DENSE_TESTS (22 patterns) via t.make(name, n) — dense → Householder → bidiag
       PAPER_TESTS (46 patterns) via t.make(name, n) — bidiagonal directly
    """
    sizes = [100] if quick else [100, 200, 400]
    out = []
    saved = sys.stdout; sys.stdout = io.StringIO()
    try:
        import test_dense_to_bidiag as t
        for pat in t.DENSE_TESTS:
            for sz in sizes:
                try:
                    d, e = t.make(pat, sz)
                except Exception:
                    continue
                out.append(("dense", f"{pat}@{sz}", d, e))
        for pat in t.PAPER_TESTS:
            for sz in sizes:
                try:
                    d, e = t.make(pat, sz)
                except Exception:
                    continue
                out.append(("dense", f"paper_{pat}@{sz}", d, e))
    except Exception as ex:
        sys.stdout = saved
        print(f"  dense suite import failed: {ex}")
    finally:
        sys.stdout = saved
    return out


def specs_glued(quick):
    """test_glued_synth.make_glued(base_name, n_block, glue_type) — glue 5 copies."""
    block_sizes = [50] if quick else [50, 100]
    glue_types = ["sqrteps", "1e-10"]
    out = []
    saved = sys.stdout; sys.stdout = io.StringIO()
    try:
        import test_glued_synth as t
        skip = getattr(t, "SKIP_PATTERNS", set())
        for pat in t.adv_names:
            if pat in skip:
                continue
            for nb in block_sizes:
                for gt in glue_types:
                    try:
                        d, e = t.make_glued(pat, nb, gt)
                    except Exception:
                        continue
                    n = len(d)
                    out.append(("glued",
                                f"glued_{pat}_{nb}x{gt}@{n}", d, e))
    except Exception as ex:
        sys.stdout = saved
        print(f"  glued suite import failed: {ex}")
    finally:
        sys.stdout = saved
    return out


def specs_synth_quick():
    """test_synth_pract: small synth + pract subset.
       generate_synth(max_dim=100, dim_step=1, econd_list=None, seed_base=42)
       generate_pract() yields STCollection-derived bidiagonals.
    """
    out = []
    saved = sys.stdout; sys.stdout = io.StringIO()
    try:
        import test_synth_pract as t
        # synth: cap dimensions to keep the run snappy
        n_synth = 0
        for spec in t.generate_synth(max_dim=100, dim_step=20):
            name, d, e = spec
            if len(d) > 200:
                continue
            out.append(("synth", f"S_{name}", d, e))
            n_synth += 1
            if n_synth >= 250:
                break
        # pract: take all
        for spec in t.generate_pract():
            name, d, e = spec
            if len(d) > 600:
                continue
            out.append(("synth", f"P_{name}", d, e))
    except Exception as ex:
        sys.stdout = saved
        print(f"  synth_pract import failed: {ex}")
    finally:
        sys.stdout = saved
    return out


def load_python_baseline_pass(suite_name):
    """Return {test_name: pass_bool} for the patched Python+Fortran baseline,
    if available. For the evaluate suite we pull from
    experiments/ablation_results.json (variant=patched). For other suites we
    don't have a cached pass map so we return {} (the per-test diff column
    will then read 0 — meaning we have no Py+Fo reference, not zero diffs)."""
    if suite_name != "evaluate":
        return {}
    path = os.path.join(ROOT, "experiments", "ablation_results.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        d = json.load(f)
    pat = next((v for v in d["variants"] if v["variant"] == "patched"), None)
    if not pat:
        return {}
    return {r["name"]: (r["status"] == "ok" and r["ok"]) for r in pat["results"]
            if r.get("status") == "ok"}


def run_python_one(d, e, timeout=60):
    """Subprocess-isolated Python+Fortran run on a single (d, e), returning
    the same dict shape as run_fortran_one (status / res / ortU / ortV / ok)."""
    import tempfile
    fd, npz_path = tempfile.mkstemp(suffix=".npz")
    os.close(fd)
    try:
        np.savez(npz_path,
                 d=np.asarray(d, dtype=np.float64),
                 e=np.asarray(e, dtype=np.float64))
        src = (
            "import os, sys, json, time, io, numpy as np\n"
            f"sys.path.insert(0, {ROOT!r})\n"
            "saved=sys.stdout; sys.stdout=io.StringIO()\n"
            "from mr3_gk import bidiag_svd\n"
            "sys.stdout=saved\n"
            f"z = np.load({npz_path!r})\n"
            "d = z['d'].astype(np.float64).copy()\n"
            "e = z['e'].astype(np.float64).copy()\n"
            "n = len(d)\n"
            "t0 = time.perf_counter()\n"
            "sigma, U, V, info = bidiag_svd(d.copy(), e.copy())\n"
            "dt = time.perf_counter() - t0\n"
            "EPS = 2.2204460492503131e-16\n"
            "B = np.diag(d) + np.diag(e, 1)\n"
            "bnorm = abs(d[-1]) if n>0 else 0.0\n"
            "for i in range(n-1):\n"
            "    bnorm = max(bnorm, abs(d[i]) + abs(e[i]))\n"
            "if bnorm == 0.0: bnorm = 1.0\n"
            "try:\n"
            "    res = float(np.max(np.abs(B - U @ np.diag(sigma) @ V.T)) / (bnorm*n*EPS))\n"
            "except Exception:\n"
            "    res = float('nan')\n"
            "try:\n"
            "    ou = float(np.max(np.abs(U.T @ U - np.eye(n))) / (n*EPS))\n"
            "    ov = float(np.max(np.abs(V.T @ V - np.eye(n))) / (n*EPS))\n"
            "except Exception:\n"
            "    ou = ov = float('nan')\n"
            "ok = bool(np.isfinite(res) and np.isfinite(ou) and np.isfinite(ov)\n"
            "          and res<=7.0 and ou<=5.0 and ov<=5.0)\n"
            "print('__OK__:' + json.dumps({'res':res,'ortU':ou,'ortV':ov,'ok':ok,'dt':dt,'info':int(info)}))\n"
        )
        try:
            r = subprocess.run([sys.executable, "-c", src],
                               capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return {"status": "timeout"}
        if r.returncode != 0:
            return {"status": "crash", "rc": r.returncode}
        line = next((ln for ln in r.stdout.splitlines() if ln.startswith("__OK__:")), None)
        if line is None:
            return {"status": "noresult"}
        payload = json.loads(line[len("__OK__:"):])
        payload["status"] = "ok"
        return payload
    finally:
        try: os.unlink(npz_path)
        except OSError: pass


def cross_check_failures(suite_summary, max_specs=20):
    """For Fortran-failures in any suite (where we don't have a Py+Fo cached
    pass map), re-run those specs through Py+Fortran and report whether they
    fail there too (i.e. confirm 'similar errors')."""
    print("\n=== cross-checking Fortran-failures vs Python+Fortran ===")
    # we need to access the original (d, e) too, which we didn't store in the
    # results json; re-run the spec generators to recover them.
    name_to_de = {}
    # re-enumerate every spec in the same order used originally
    for s in suite_summary:
        # we passed `specs_*` lists in run_suite; but those (d, e) are not in
        # `s['results']`. So pull the failing names and re-derive (d,e) by the
        # spec generators.
        pass
    # Instead, we'll handle this in main(): recompute failing specs there.
    return None


# -------- main loop --------

def run_suite_with_de(specs, suite_label, log_each=False):
    """Variant of run_suite that keeps (d, e) attached to each result so we
    can cross-check failing specs through Python+Fortran."""
    print(f"\n=== suite: {suite_label}  ({len(specs)} specs) ===")
    if not specs:
        return {"suite": suite_label, "n": 0, "results": []}, []
    py_pass = load_python_baseline_pass(suite_label)
    out = []
    de_by_name = {}
    n_pass = n_fail = n_crash = n_timeout = 0
    n_diff_with_py = 0
    t_start = time.perf_counter()
    for i, (suite, name, d, e) in enumerate(specs):
        de_by_name[name] = (d, e)
        r = run_fortran_one(np.asarray(d, dtype=np.float64),
                            np.asarray(e, dtype=np.float64))
        r["name"] = name; r["suite"] = suite
        out.append(r)
        if r["status"] == "ok":
            if r["ok"]:
                n_pass += 1
            else:
                n_fail += 1
            if name in py_pass and py_pass[name] != r["ok"]:
                n_diff_with_py += 1
        elif r["status"] == "timeout":
            n_timeout += 1
        else:
            n_crash += 1
        if log_each or (i + 1) % 50 == 0:
            print(f"  [{i+1}/{len(specs)}] pass={n_pass} fail={n_fail} "
                  f"crash={n_crash} timeout={n_timeout} diff_with_py={n_diff_with_py}")
    wall = time.perf_counter() - t_start
    print(f"  total: pass={n_pass} fail={n_fail} crash={n_crash} "
          f"timeout={n_timeout}  ({wall:.1f}s)")
    if py_pass:
        py_p = sum(1 for v in py_pass.values() if v)
        py_f = len(py_pass) - py_p
        print(f"  Python+Fortran (cached) on same suite: pass={py_p} fail={py_f}")
        print(f"  Per-test PASS/FAIL diffs vs Py+Fo: {n_diff_with_py}")
    summary = {"suite": suite_label, "n": len(specs), "n_pass": n_pass,
               "n_fail": n_fail, "n_crash": n_crash, "n_timeout": n_timeout,
               "n_diff_with_py": n_diff_with_py, "wall_s": wall, "results": out}
    return summary, de_by_name


def run_suite(specs, suite_label, log_each=False):
    print(f"\n=== suite: {suite_label}  ({len(specs)} specs) ===")
    if not specs:
        return {"suite": suite_label, "n": 0, "results": []}
    py_pass = load_python_baseline_pass(suite_label)
    out = []
    n_pass = n_fail = n_crash = n_timeout = 0
    n_diff_with_py = 0
    t_start = time.perf_counter()
    for i, (suite, name, d, e) in enumerate(specs):
        r = run_fortran_one(np.asarray(d, dtype=np.float64),
                            np.asarray(e, dtype=np.float64))
        r["name"] = name
        r["suite"] = suite
        out.append(r)
        if r["status"] == "ok":
            if r["ok"]:
                n_pass += 1
            else:
                n_fail += 1
            if name in py_pass and py_pass[name] != r["ok"]:
                n_diff_with_py += 1
        elif r["status"] == "timeout":
            n_timeout += 1
        else:
            n_crash += 1
        if log_each or (i + 1) % 50 == 0:
            print(f"  [{i+1}/{len(specs)}] pass={n_pass} fail={n_fail} "
                  f"crash={n_crash} timeout={n_timeout} diff_with_py={n_diff_with_py}")
    wall = time.perf_counter() - t_start
    print(f"  total: pass={n_pass} fail={n_fail} crash={n_crash} "
          f"timeout={n_timeout}  ({wall:.1f}s)")
    if py_pass:
        py_p = sum(1 for v in py_pass.values() if v)
        py_f = len(py_pass) - py_p
        print(f"  Python+Fortran (same suite, evaluate.py thresholds): "
              f"pass={py_p} fail={py_f}")
        print(f"  Per-test PASS/FAIL diffs vs Py+Fo: {n_diff_with_py}")
    return {"suite": suite_label, "n": len(specs), "n_pass": n_pass,
            "n_fail": n_fail, "n_crash": n_crash, "n_timeout": n_timeout,
            "n_diff_with_py": n_diff_with_py, "wall_s": wall, "results": out}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", choices=["evaluate", "dense", "glued", "synth", "all"],
                    default="all")
    ap.add_argument("--quick", action="store_true",
                    help="restrict each suite to small sizes")
    ap.add_argument("--out",
                    default=os.path.join(ROOT, "experiments",
                                         "evaluate_fortran_results.json"))
    ap.add_argument("--log-each", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(MR3GK_RUN):
        print(f"ERROR: {MR3GK_RUN} not built. cd mr3gk_fortran && bash build.sh")
        sys.exit(1)

    summary = []
    de_maps = {}
    if args.suite in ("evaluate", "all"):
        s, de = run_suite_with_de(specs_evaluate(args.quick), "evaluate",
                                  args.log_each)
        summary.append(s); de_maps["evaluate"] = de
    if args.suite in ("dense", "all"):
        s, de = run_suite_with_de(specs_dense(args.quick), "dense",
                                  args.log_each)
        summary.append(s); de_maps["dense"] = de
    if args.suite in ("glued", "all"):
        s, de = run_suite_with_de(specs_glued(args.quick), "glued",
                                  args.log_each)
        summary.append(s); de_maps["glued"] = de
    if args.suite in ("synth", "all"):
        s, de = run_suite_with_de(specs_synth_quick(), "synth",
                                  args.log_each)
        summary.append(s); de_maps["synth"] = de

    # Cross-check Fortran failures vs Python+Fortran
    cross_rows = []
    for s in summary:
        de = de_maps.get(s["suite"], {})
        fails = [r for r in s["results"]
                 if r["status"] == "ok" and not r.get("ok", False)]
        if not fails:
            continue
        fails.sort(key=lambda r: -max(r["res"], r["ortU"], r["ortV"]))
        print(f"\n=== cross-check Py+Fo on {len(fails)} Fortran failures "
              f"(suite={s['suite']}) ===")
        for r in fails[:30]:
            d, e = de.get(r["name"], (None, None))
            if d is None:
                continue
            pr = run_python_one(d, e)
            same = (pr.get("status") == "ok" and pr.get("ok") == r["ok"])
            status_match = "MATCH" if same else "DIFF"
            row = {"suite": s["suite"], "name": r["name"],
                   "fortran": {"res": r["res"], "ortU": r["ortU"], "ortV": r["ortV"], "ok": r["ok"]},
                   "python":  pr, "status_match": status_match}
            cross_rows.append(row)
            py_str = (f"py: res={pr.get('res', float('nan')):8.2f} "
                      f"ortU={pr.get('ortU', float('nan')):8.2f} "
                      f"ortV={pr.get('ortV', float('nan')):8.2f}"
                      if pr.get("status") == "ok" else f"py: {pr.get('status')}")
            print(f"  [{status_match}] {r['name']:42s}  "
                  f"fo: res={r['res']:7.2f} ortU={r['ortU']:7.2f} ortV={r['ortV']:7.2f}  |  {py_str}")

    with open(args.out, "w") as f:
        json.dump({"suites": summary,
                   "thresholds": {"res": RES_THRESH, "ortU": ORTHO_THRESH,
                                  "ortV": ORTHO_THRESH},
                   "cross_check": cross_rows},
                  f, indent=1, default=float)
    print(f"\nwrote {args.out}")

    print("\n=== summary ===")
    print(f"{'suite':<10}  {'n':>6}  {'pass':>5}  {'fail':>5}  {'crash':>5}  "
          f"{'timeout':>7}  {'diff_py':>7}  {'wall_s':>7}")
    for s in summary:
        print(f"{s['suite']:<10}  {s['n']:>6}  "
              f"{s.get('n_pass', 0):>5}  {s.get('n_fail', 0):>5}  "
              f"{s.get('n_crash', 0):>5}  {s.get('n_timeout', 0):>7}  "
              f"{s.get('n_diff_with_py', 0):>7}  {s.get('wall_s', 0):>7.1f}")


if __name__ == "__main__":
    main()
