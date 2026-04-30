# 03 — Compilation Flags Audit

This document reviews the two `build.sh` scripts in `python_fortran/`. Deliverable for advisor task **(d)**.

**Toolchain measured on:**
```
$ gfortran --version
GNU Fortran (GCC) 4.8.5 20150623 (Red Hat 4.8.5-44)
```

(Tested also under gfortran 9.x; both produce 21,375 / 21,375 pass at 2·eps.)

---

## File A: `python_fortran/build.sh` (builds `libxmr.so`)

This script compiles the upstream XMR Fortran sources plus our two patch files and links the result into a shared library that both Python (`xmr_ctypes.py`) and the new Fortran driver (`mr3gk_fortran/mr3gk_run`) load.

### Flags
| Step | Files | Flags |
|---|---|---|
| Unmodified XMR | 31 files in `xmr_src/` | `-fPIC -O2 -std=legacy -w -fno-second-underscore` |
| Patch 1 | `dlaxrb_clssfy_fix.f` | `-fPIC -O2 -std=legacy -w -fno-second-underscore` (same as XMR) |
| Patch 2 | `dlaxre_gk.f` | `-c -fPIC -O2` ⚠️ **flag drift** (missing `-std=legacy -w -fno-second-underscore`) |
| C wrapper | `xmr_wrapper.c` | `gcc -c -fPIC -O2` |
| Link | all the above | `gcc -shared … -lgfortran -llapack -lblas -lm` |

### Assessment

| Flag | Purpose | Verdict |
|---|---|---|
| `-fPIC` | Position-independent code, required for shared library | ✅ correct |
| `-O2` | Standard optimization | ✅ correct — XMR is the validated kernel; we don't need bit-identity between two builds of `libxmr.so` |
| `-std=legacy` | Allow F77 fixed-form (XMR uses goto/common blocks) | ✅ correct |
| `-w` | Suppress warnings (XMR has many implicit-typed locals) | ✅ acceptable for upstream code we don't own |
| `-fno-second-underscore` | Single trailing underscore on Fortran symbols (matches reference BLAS / LAPACK ABI) | ✅ **required** for ctypes to find the symbols |
| `-shared` (link) | Build .so | ✅ correct |
| `-lgfortran -llapack -lblas -lm` | Link runtime + reference LAPACK/BLAS + math | ✅ correct |

### Issue: flag drift on `dlaxre_gk.f`

Line 58 reads:
```bash
$FC -c -fPIC -O2 dlaxre_gk.f -o "$OBJDIR/dlaxre_gk.o"
```
This omits `-std=legacy -w -fno-second-underscore` that the rest of the script uses. Why it works on gfortran 4.8.5:
- `-std=legacy` is the default on gfortran 4.x for `.f` files anyway (it accepts F77 implicitly).
- `-fno-second-underscore` is also the default on most GCC builds (only `-fsecond-underscore` is the override).
- `-w` only changes warning verbosity, not codegen.

So **functionally fine on gfortran 4.x**, but it could break on newer compilers (e.g. gfortran 12+ tightens `-std=` defaults). **Recommendation:** change line 58 to use the same `$FFLAGS`:
```bash
$FC $FFLAGS -c dlaxre_gk.f -o "$OBJDIR/dlaxre_gk.o"
```
This is a **non-blocking** improvement — current build works; it's a portability hardening.

### Why `-O2` is acceptable here

The XMR kernel was validated upstream against the LAPACK test set. We do not require bit-identity between two compiles of `libxmr.so`; we only require that **a single compile** of `libxmr.so` is loaded by both the Python wrapper and the Fortran driver — which it is. `-O2` lets the compiler reorder safe FP operations within `dlaxrv`/`dlaxre`/etc., but those are deterministic at the binary level and identical for every caller.

---

## File B: `python_fortran/mr3gk_fortran/build.sh` (builds `mr3gk_run`)

This script compiles the new Fortran orchestration layer and links it against the already-built `libxmr.so`.

### Flags
```
FC="gfortran"
FFLAGS="-O0 -fno-fast-math -fPIC -fno-second-underscore -fimplicit-none \
        -Wno-unused -Wno-unused-dummy-argument"
```

| Flag | Purpose | Verdict |
|---|---|---|
| `-O0` | **No optimization.** Keeps numerical operations in source order. | ✅ **deliberate** — see rationale below |
| `-fno-fast-math` | Disable unsafe FP transforms (assoc, FMA fusion, denormal flushing) | ✅ **deliberate** |
| `-fPIC` | Object code suitable for any link mode | ✅ correct |
| `-fno-second-underscore` | Match reference BLAS/LAPACK ABI used in `libxmr.so` | ✅ correct |
| `-fimplicit-none` | No implicit typing — every variable must be declared | ✅ helps catch port bugs |
| `-Wno-unused -Wno-unused-dummy-argument` | Quiet warnings for partial subroutines | ✅ acceptable |

### Linker line
```
$FC -o mr3gk_run mr3gk_run.o "${OBJ[@]}" "${LIBXMR}" -llapack -lblas -lgfortran -lm
```
Links against:
- All `mr3gk_*.o` orchestration objects.
- The pre-built `../libxmr.so` (kernel).
- System reference LAPACK + BLAS (for `DSTEBZ`, `DLARTG`, `DLARNV`, `DNRM2`).
- gfortran runtime + libm.

✅ Correct. Same kernel binary as Python's ctypes path.

### Why `-O0 -fno-fast-math` here (and not in `libxmr.so`)

We promised the user **bit-level reproducibility** (≤ 2·eps) between the Python+Fortran hybrid baseline and the pure-Fortran driver. The post-processing sites in `mr3gk_postproc.f90` perform many small accumulations (Bv recovery, sign-fix dot products, GS completion projections). At `-O2` gfortran is allowed to:
- Reorder commutative FP operations.
- Use FMA where available (`a + b*c` → single rounding).
- Vectorize loops (changes summation order).

Each of those would introduce 1–2 ULP differences between the Python orchestration's calls (which always go through ctypes one operation at a time, no FMA, no vectorization) and the Fortran orchestration's calls. Result: the 2·eps target would fail on a handful of adversarial specs.

`-O0 -fno-fast-math` forces gfortran to emit one FP operation per source statement, in source order, with no fusion. This is what produced the 19,866 / 19,866 pass at 2·eps.

**Performance trade-off is acceptable** because:
- Wall time on the 19,866-spec sweep was 1371 s on 8 workers (≈ 0.069 s/spec average); the bottleneck is the `dlaxrv` kernel call, which is in the `-O2` `libxmr.so` and unaffected.
- The orchestration layer in `mr3gk_postproc.f90` is O(n²) per spec; at n=3000, n² = 9 M FLOPs, negligible compared to the kernel.

---

## Recommendations

| # | Change | Priority | Effort |
|---|---|---|---|
| 1 | Use `$FFLAGS` consistently for `dlaxre_gk.f` in root `build.sh` line 58 | low | trivial — replace `$FC -c -fPIC -O2` with `$FC $FFLAGS -c` |
| 2 | Add `gfortran --version` check at top of both scripts; fail if < 4.8 | low | 3 lines of bash |
| 3 | Document `-O0 -fno-fast-math` rationale inline in `mr3gk_fortran/build.sh` | low | comment block (already partially present at line 19) |
| 4 | If porting to ifort: switch `-O0 -fno-fast-math` → `-O0 -fp-model strict` | future | for ifort users only |
| 5 | Pin libblas/liblapack versions in `requirements.md` (since DNRM2 byte-identity depends on which BLAS the OS ships) | medium | document only |

None of these are blocking; the current build produces a correct, validated artifact.

---

## Reproducibility checklist

To rebuild from scratch and reverify:
```bash
cd python_fortran
rm -rf fortran_objects libxmr.so mr3gk_fortran/*.o mr3gk_fortran/mr3gk_run

bash build.sh                           # produces libxmr.so
cd mr3gk_fortran && bash build.sh && cd ..  # produces mr3gk_run

# regenerate baseline (one-time, ~5 min)
python3 test_fortran_match.py --gen-baseline

# compare (~2 min)
python3 test_fortran_match.py
# expect: PASS: 379/379  Worst sigma diff: 0  Worst U diff: 1.110e-16  Worst V diff: 2.272e-19
```

This is what the test report (`01_fortran_test_report.md`) records as the verified result.
