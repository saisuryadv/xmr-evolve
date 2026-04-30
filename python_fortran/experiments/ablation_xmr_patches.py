#!/usr/bin/env python3
"""
Ablation: run the MR3-GK SVD evaluation suite under four library variants
to quantify the empirical impact of each XMR patch.

Variants (linked from python_fortran/):
  patched      libxmr_both_patch.so   our dlaxre_gk.f          + our dlaxrb_clssfy_fix.f
  origre       libxmr_origre.so       upstream xmr_src/dlaxre.f + our dlaxrb_clssfy_fix.f
  origcls      libxmr_origcls.so      our dlaxre_gk.f           + upstream xmr_src/dlaxrb_clssfy.f
  origboth     libxmr_origboth.so     upstream dlaxre.f         + upstream dlaxrb_clssfy.f

For each variant, we run the same 379-spec corpus that evaluate.py uses
(90 adversarial patterns × 4 sizes [10,100,200,400] + 19 STCollection files)
and record per-spec (res, ortU, ortV, status). A subprocess wrapper survives
hangs/crashes (e.g. the upstream branch of dlaxre will hit the dead `IF(.FALSE.)`
on T_GK matrices and likely produce NaN / orthogonality blowup).

Workflow:
  1. (Once) Generate all (d, e) for the 379 corpus by importing full_eval and
     dumping each spec to /tmp/ablation_specs/<safe_name>.npz.
  2. For each variant, set XMR_LIB_PATH and run a lightweight worker subprocess
     per spec; each worker only loads its .npz, calls mr3_gk.bidiag_svd, and
     prints one JSON line.

Usage:
  python3 experiments/ablation_xmr_patches.py            # all variants, full 379
  python3 experiments/ablation_xmr_patches.py --quick    # n=100 only
  python3 experiments/ablation_xmr_patches.py --variant origboth
  python3 experiments/ablation_xmr_patches.py --regen    # force regenerate spec cache
"""
import os
import sys
import json
import argparse
import subprocess
import glob
import io
import time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

VARIANTS = {
    "patched":  "libxmr_both_patch.so",
    "origre":   "libxmr_origre.so",
    "origcls":  "libxmr_origcls.so",
    "origboth": "libxmr_origboth.so",
}

CACHE_DIR = "/tmp/ablation_specs"


def safe(name):
    return name.replace("/", "_").replace("@", "_at_")


def regenerate_cache(quick):
    """Generate all (d, e) for the 379 corpus once, dumping to /tmp/ablation_specs/.
    Suppresses full_eval's import-time test loop output."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    saved = sys.stdout
    sys.stdout = io.StringIO()
    try:
        from full_eval import make, adv_names
        from evaluate import load_stcoll
    finally:
        sys.stdout = saved
    sizes = [100] if quick else [10, 100, 200, 400]
    specs = []
    for sz in sizes:
        for pat in adv_names:
            try:
                d, e = make(pat, sz)
            except Exception:
                continue
            name = f"{pat}@{sz}"
            path = os.path.join(CACHE_DIR, safe(name) + ".npz")
            np.savez(path,
                     d=np.asarray(d, dtype=np.float64),
                     e=np.asarray(e, dtype=np.float64))
            specs.append((name, path))
    if not quick:
        for f in sorted(glob.glob(os.path.join(ROOT, "stcollection", "B_*.dat"))):
            name, d, e = load_stcoll(f)
            full = f"ST_{name}"
            path = os.path.join(CACHE_DIR, safe(full) + ".npz")
            np.savez(path,
                     d=np.asarray(d, dtype=np.float64),
                     e=np.asarray(e, dtype=np.float64))
            specs.append((full, path))
    return specs


def load_specs_from_cache(quick):
    sizes = [100] if quick else [10, 100, 200, 400]
    specs = []
    saved = sys.stdout; sys.stdout = io.StringIO()
    try:
        from full_eval import adv_names
    finally:
        sys.stdout = saved
    for sz in sizes:
        for pat in adv_names:
            name = f"{pat}@{sz}"
            path = os.path.join(CACHE_DIR, safe(name) + ".npz")
            if os.path.exists(path):
                specs.append((name, path))
    if not quick:
        for f in sorted(glob.glob(os.path.join(ROOT, "stcollection", "B_*.dat"))):
            base = os.path.basename(f).replace(".dat", "")
            full = f"ST_{base}"
            path = os.path.join(CACHE_DIR, safe(full) + ".npz")
            if os.path.exists(path):
                specs.append((full, path))
    return specs


WORKER_SRC = r'''
import os, sys, json, time, io
import numpy as np
sys.path.insert(0, "%(root)s")
saved = sys.stdout; sys.stdout = io.StringIO()
from mr3_gk import bidiag_svd
sys.stdout = saved

z = np.load("%(npz)s")
d = z["d"].astype(np.float64).copy()
e = z["e"].astype(np.float64).copy()
n = len(d)

t0 = time.perf_counter()
sigma, U, V, info = bidiag_svd(d.copy(), e.copy())
dt = time.perf_counter() - t0

EPS = 2.2204460492503131e-16
B = np.diag(d) + np.diag(e, 1)
bnorm = abs(d[-1]) if n > 0 else 0.0
for i in range(n - 1):
    bnorm = max(bnorm, abs(d[i]) + abs(e[i]))
if bnorm == 0.0: bnorm = 1.0

# Avoid NaN: catch overflow
try:
    res = float(np.max(np.abs(B - U @ np.diag(sigma) @ V.T)) / (bnorm * n * EPS))
except Exception:
    res = float("nan")
try:
    ou = float(np.max(np.abs(U.T @ U - np.eye(n))) / (n * EPS))
    ov = float(np.max(np.abs(V.T @ V - np.eye(n))) / (n * EPS))
except Exception:
    ou = ov = float("nan")

ok = bool(np.isfinite(res) and np.isfinite(ou) and np.isfinite(ov)
          and res <= 7.0 and ou <= 5.0 and ov <= 5.0)
print("__OK__:" + json.dumps({"res": res, "ortU": ou, "ortV": ov,
                              "ok": ok, "dt": dt, "info": int(info)}))
'''


def run_one(name, npz, lib_path, timeout):
    src = WORKER_SRC % {"root": ROOT, "npz": npz}
    env = dict(os.environ)
    env["XMR_LIB_PATH"] = lib_path
    env.pop("PYTHONSTARTUP", None)
    try:
        r = subprocess.run([sys.executable, "-c", src],
                           capture_output=True, text=True,
                           timeout=timeout, env=env)
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "name": name}
    if r.returncode != 0:
        return {"status": "crash", "name": name,
                "rc": r.returncode,
                "stderr_tail": r.stderr.strip().splitlines()[-3:] if r.stderr else []}
    line = next((ln for ln in r.stdout.splitlines() if ln.startswith("__OK__:")), None)
    if line is None:
        return {"status": "noresult", "name": name,
                "stdout_tail": r.stdout.strip().splitlines()[-3:]}
    payload = json.loads(line[len("__OK__:"):])
    payload["status"] = "ok"
    payload["name"] = name
    return payload


def run_variant(variant_name, specs, timeout, log_each=False):
    lib = os.path.join(ROOT, VARIANTS[variant_name])
    if not os.path.exists(lib):
        raise SystemExit(f"missing library: {lib}")
    print(f"\n=== variant: {variant_name}  ({lib}) ===  specs={len(specs)}")
    t0 = time.perf_counter()
    out = []
    n_pass = n_fail = n_crash = n_timeout = 0
    for i, (name, npz) in enumerate(specs):
        r = run_one(name, npz, lib, timeout)
        out.append(r)
        if r["status"] == "ok":
            if r["ok"]:
                n_pass += 1
            else:
                n_fail += 1
        elif r["status"] == "timeout":
            n_timeout += 1
        else:
            n_crash += 1
        if log_each or (i + 1) % 25 == 0:
            print(f"  [{i+1}/{len(specs)}] pass={n_pass} fail={n_fail} "
                  f"crash={n_crash} timeout={n_timeout}  ({name} : {r['status']})")
    wall = time.perf_counter() - t0
    print(f"  total: pass={n_pass}  fail={n_fail}  crash={n_crash}  "
          f"timeout={n_timeout}  ({wall:.1f}s)")
    return {"variant": variant_name, "wall_s": wall, "results": out,
            "n_pass": n_pass, "n_fail": n_fail,
            "n_crash": n_crash, "n_timeout": n_timeout}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--quick", action="store_true",
                   help="n=100 only (fewer specs)")
    p.add_argument("--variant", choices=list(VARIANTS),
                   help="run only one variant")
    p.add_argument("--timeout", type=int, default=60,
                   help="per-spec timeout seconds")
    p.add_argument("--regen", action="store_true",
                   help="regenerate /tmp/ablation_specs cache")
    p.add_argument("--out", default=os.path.join(os.path.dirname(__file__),
                                                  "ablation_results.json"))
    p.add_argument("--log-each", action="store_true")
    args = p.parse_args()

    if args.regen or not os.path.exists(CACHE_DIR) or len(os.listdir(CACHE_DIR)) == 0:
        print("Generating /tmp/ablation_specs cache...")
        specs = regenerate_cache(args.quick)
    else:
        specs = load_specs_from_cache(args.quick)
        if len(specs) == 0:
            specs = regenerate_cache(args.quick)
    print(f"Specs: {len(specs)}  (quick={args.quick})  cache={CACHE_DIR}")

    variants = [args.variant] if args.variant else list(VARIANTS)

    summary = []
    for v in variants:
        summary.append(run_variant(v, specs, args.timeout, args.log_each))

    with open(args.out, "w") as f:
        json.dump({"specs_count": len(specs), "quick": args.quick,
                   "variants": summary}, f, indent=1)
    print(f"\nwrote {args.out}")

    print("\n=== summary ===")
    print(f"{'variant':<10}  {'pass':>5}  {'fail':>5}  {'crash':>5}  {'timeout':>7}  {'wall_s':>7}")
    for s in summary:
        print(f"{s['variant']:<10}  {s['n_pass']:>5}  {s['n_fail']:>5}  "
              f"{s['n_crash']:>5}  {s['n_timeout']:>7}  {s['wall_s']:>7.1f}")


if __name__ == "__main__":
    main()
