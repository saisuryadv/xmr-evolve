# initial_code/ — the starting point for `docs/REPRODUCE.md`

This directory is the **canonical Checkpoint A starting state** referenced
by `python_fortran/docs/REPRODUCE.md`. To replay the prompt sequence that
evolves a naive Python+Fortran wrapper into the working 379/379 MR3-GK
bidiagonal SVD, copy this folder to a fresh location and feed the prompts
to a coding agent against that copy.

**Important:** `xmr_src/` here is **Paul Willems' code, untouched**. Both
XMR bugs we describe in `docs/BUGREPORT.md` are present in this source as
they originally shipped:

- `xmr_src/dlaxre.f:268` still has the disabled GK branch
  (`IF( .FALSE. )THEN`) — Bug #1.
- `xmr_src/dlaxrb_clssfy.f:361,412,438` still has the unguarded
  `MIN( AVGTHRESH, ABSMAX*GAPTHRESH )` test — Bug #2.

Neither of our patched files (`dlaxre_gk.f`, `dlaxrb_clssfy_fix.f`) is
present here. The agent re-creates them in prompts A3 and A10
respectively, faithfully reproducing the bug-fixing arc.

## What's in here

| Path | What it is | Provenance |
|---|---|---|
| `xmr_src/` | Paul Willems' XMR Fortran source (44 `.f` files, ~14K LOC) | **Upstream, untouched.** Source: Willems' `tdsolver/xmr/SRC/O` directory (the "Optimized" build flavor). Verified byte-identical to the in-tree `python_fortran/xmr_src/` and to Willems' published archive — every one of the 44 `.f` files matches exactly. (Willems' archive also ships `Dcl/` 50-file and `Oc/` 47-file variants; we don't use those.) |
| `stcollection/` | 19 STCollection real-world bidiagonal matrices (`.dat`) | From github.com/oamarques/STCollection (Marques et al.). Untouched. |
| `evaluate.py` | 379-test scoring harness | The current `evaluate.py` from project HEAD. Reports `TOTAL: X/379` against the `mr3_gk.bidiag_svd` symbol that the agent will create. |
| `full_eval.py` | 90 adversarial bidiagonal-matrix generators | The current `full_eval.py` from HEAD. Defines `make(name, n)` and `adv_names`. |

## What's NOT in here (the agent must produce these)

| What | Where it goes | Which prompt creates it |
|---|---|---|
| `mr3_gk.py` | `bidiag_svd(d, e)` Python entry point | A1 |
| `xmr_wrapper.c` | C wrapper exposing `xmr_eigenvectors(...)` | A1 |
| `xmr_ctypes.py` | Python ctypes binding to `libxmr.so` | A1 |
| `build.sh` | Compiles `xmr_src/` + wrapper + patches into `libxmr.so` | A1 |
| `dlaxre_gk.f` | Patched copy of `xmr_src/dlaxre.f` activating GK branch | A3 |
| `dlaxrb_clssfy_fix.f` | Patched copy of `xmr_src/dlaxrb_clssfy.f` with AVGTHRESH guard | A10 |

## Verifying the upstream source is pristine

```
diff -qr xmr_src ../xmr_src
# (empty output = byte-identical to in-tree xmr_src/)
```

To confirm the upstream source compiles standalone:

```
cd initial_code
mkdir -p _objs
for f in xmr_src/*.f; do
   gfortran -fPIC -O2 -std=legacy -w -fno-second-underscore \
      -c "$f" -o "_objs/$(basename $f .f).o"
done
ls _objs/*.o | wc -l   # → 44
```

Expect zero errors and zero warnings (the `-w` flag suppresses them).

## Running the recipe

See `python_fortran/docs/REPRODUCE.md` for the full 10-prompt sequence
(prompts A0–A11). The agent is launched against a copy of this directory:

```
cp -r initial_code /tmp/repro_workspace
cd /tmp/repro_workspace
# feed prompts in order to a fresh `claude --print` agent for each
```

After all prompts the workspace should contain everything needed to run
`python3 evaluate.py` and report **TOTAL: 379/379**.

## Why is `xmr_src/` already here, but `dlaxre_gk.f` isn't?

The upstream Willems code is the *unmodified algorithmic kernel* — XMR is
eigenvalue-only, so even untouched it compiles into a working tridiagonal
eigenpair solver. The two files we ship as patches (`dlaxre_gk.f` and
`dlaxrb_clssfy_fix.f`) are *replacements* for two specific upstream files,
applied at link time. They represent the two XMR bugs documented in
`docs/BUGREPORT.md` and the agent re-creates them in prompts A3 and A10
respectively.

So the cleanest "what came from Willems vs what we wrote" distinction is:

- **Willems gave us:** `xmr_src/` (the kernel).
- **Marques gave us:** `stcollection/` (the real-world test matrices).
- **We wrote:** evaluator + adversarial generators (already in this folder
  as `evaluate.py` / `full_eval.py`), plus everything in the "NOT in here"
  table above (the agent re-creates those).

The agent's job is to produce the orchestration layer + the two Fortran
patches, on top of the unmodified kernel.
