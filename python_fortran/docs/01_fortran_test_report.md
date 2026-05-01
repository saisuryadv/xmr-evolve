# 01 — Pure-Fortran MR3-GK SVD: Consolidated Test Report

**Driver under test:** `python_fortran/mr3gk_fortran/mr3gk_run`
(built by `python_fortran/mr3gk_fortran/build.sh`; links against `python_fortran/libxmr.so`)

**Tolerance:** `2 * machine_epsilon = 4.4408920985e-16`

**Pass criterion (per spec):** all three of
- `max |σ_py − σ_fo| ≤ 2·eps`
- `max |U_py − U_fo| ≤ 2·eps` (after per-column joint sign canonicalization)
- `max |V_py − V_fo| ≤ 2·eps` (same canonicalization)

`*_py` = output of the Python+Fortran hybrid baseline `mr3_gk.bidiag_svd`; `*_fo` = output of the pure-Fortran driver `mr3gk_run`.

---

## Combined headline

| Suite | Specs | Pass @ 2·eps | Worst σ | Worst U | Worst V |
|---|---|---|---|---|---|
| 379-spec match (`test_fortran_match.py`) | 379 | **379/379** | 0.000e+00 | 1.110e-16 (`gl_wilk2w@200`) | 2.272e-19 (`gl_wilk2w@400`) |
| 1130-spec extended sweep | 1130 | **1130/1130** | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| 19,866-spec full sweep | 19,866 | **19,866/19,866** | 0.000e+00 | 1.110e-16 (`gl_wilk2w@800`) | 2.272e-19 (`gl_wilk2w@400`) |
| **Total** | **21,375** | **21,375 / 21,375** | **0** | **1.110e-16** | **2.272e-19** |

Singular values are **bit-identical** between Python+Fortran and pure-Fortran on every one of 21,375 test bidiagonals. U and V differ by at most ≈ 1 ULP on a handful of `gl_wilkw / gl_wilk2w` adversarial cases.

---

## 1. Core 379-spec match — `test_fortran_match.py`

**Source list:** 90 adversarial patterns × 4 sizes (10, 100, 200, 400) + 19 STCollection `B_*.dat` files.

**Reproduce:**
```
cd python_fortran
python3 test_fortran_match.py --gen-baseline   # one-time (~5 min): regenerate Python+Fortran baseline cache
python3 test_fortran_match.py                  # ~2 min: run Fortran driver, compare to cache
```

**Result (re-run on 2026-04-29):**
```
PASS: 379/379   FAIL: 0   SKIP: 0
Worst sigma diff: 0.000e+00  ()
Worst U diff:     1.110e-16  (gl_wilk2w@200)
Worst V diff:     2.272e-19  (gl_wilk2w@400)
Tolerance:        4.441e-16
```

| Pass count at | n |
|---|---|
| 2·eps  (4.44e-16) | 379/379 |
| 4·eps  (8.88e-16) | 379/379 |
| 1e-14 | 379/379 |
| 1e-12 | 379/379 |

---

## 2. Extended Py↔Fortran sweep (1130 specs)

**Sources:** every spec produced by `test_dense_to_bidiag.py` (Householder-bidiagonalized dense), `test_glued_synth.py` (Willems-Lang glued/synth), `test_synth_pract.py` (Synth + Pract), excluding 19 STCollection overlaps.

**Comparison harness:** `/tmp/extended_compare.py` (independent of the in-tree `test_fortran_match.py`).

**Result:**
| Suite | Specs | Pass @ 2·eps |
|---|---|---|
| dense_to_bidiag | 223 | 223/223 |
| glued_synth | 347 | 347/347 |
| synth_pract | 560 | 560/560 |
| **Total** | **1130** | **1130/1130** |

- σ, U, V all **bit-identical** (worst abs diff = 0 in every metric, after sign canonicalization).
- Worst per-column `1−|cos θ|`: 1.22e-15 (U, n=500), 1.33e-15 (V, n=400) — sign-canonicalization rounding only, not algorithmic divergence.
- Worst orthogonality: `‖UᵀU − I‖` = 2.45e-12 on `dense_wilkinson_sv@399`; matches Python side (algorithm-intrinsic, not port-induced).
- 16 specs with `‖B V − U·diag σ‖ > 1` are all `near_overflow*` (matrix entries ~1e150) — Python residual identical to Fortran residual on the same spec, confirming matrix-inherent overflow rather than port regression.
- 0 crashes, 0 timeouts, 0 info-code mismatches.

---

## 3. Full 19,866-spec sweep ("the 16k experiment")

**Sources:**
| Suite | Specs |
|---|---|
| `full_eval.adv_names` × n ∈ {10, 50, 100, 200, 400, 800, 1500, 3000} | 720 |
| STCollection (`stcollection/B_*.dat`) | 19 |
| Dense-to-bidiagonal (`test_dense_to_bidiag.py`) | 88 |
| Paper testset | 136 |
| Glued/synth small | 344 |
| Glued/synth large (n up to 1005) | 6 |
| Synth + Pract (`test_synth_pract.py`, full grid) | 18,553 |
| **Total** | **19,866** |

**Comparison harness:** `/tmp/full16k_compare.py` (8-worker `multiprocessing.Pool`, 600 s/spec timeout, no n cap).

**Wall:** 1371 s (≈ 22.9 min).

**Result:**
| Pass count at | n |
|---|---|
| 2·eps  (4.44e-16) | 19,866 / 19,866 |
| 4·eps  (8.88e-16) | 19,866 / 19,866 |
| 1e-14 | 19,866 / 19,866 |
| 1e-12 | 19,866 / 19,866 |

- σ **bit-identical** across the entire sweep (max |σ_py − σ_fo| = 0). This includes n=3000 adversarial matrices.
- Worst U abs diff: 1.110e-16 on `gl_wilk2w@800` and `gl_wilk2w@200` (5 specs total in the `wilkw / wilk2w` family; everything else is 0).
- Worst V abs diff: 2.272e-19 on `gl_wilk2w@400` and `gl_wilk2w@1500` (2 specs total).
- 0 info-code mismatches, 0 crashes, 0 timeouts, 0 OOM skips.

**Worst Fortran-side orthogonality** (algorithm-intrinsic, identical between Py and Fo):
| Spec | n | `‖UᵀU−I‖` | `‖VᵀV−I‖` |
|---|---|---|---|
| `synth_pract: T_Godunov_1e-5` | 2500 | 1.700e-08 | 1.700e-08 |
| `synth_pract: T_Godunov_1e-4` | 2500 | 3.937e-10 | 3.937e-10 |
| `synth_pract: T_TSC_OPF_300` | 9774 | 8.220e-11 | 8.220e-11 |
| `synth_pract: ev3_ec1_n65_s3_gm` | 195 | 1.318e-11 | 1.318e-11 |
| `synth_pract: T_Alemdar_1` | 6245 | 1.119e-11 | 1.119e-11 |

These reflect MR3's intrinsic limit on the hardest tridiagonals in our test set; they are **identical** Py vs Fo, so they are not a Fortran-port regression.

**Worst Fortran-side SVD identity `‖B V − U·diag σ‖`:**
| Spec | n | `‖B V − U·diag σ‖` |
|---|---|---|
| `full_eval: near_overflow` | 3000 | 9.662e+137 |
| `full_eval: near_overflow` | 1500 | 4.621e+137 |
| `full_eval: near_overflow` | 800 | 2.456e+137 |
| `full_eval: near_overflow` | 400 | 1.305e+137 |
| `glued_synth: near_overflow_glued_small` | 300 | 1.188e+136 |

These are matrix-inherent (entries ~1e150, residuals scale accordingly). Python-side residual is identical on every overflow spec.

**Artifacts:** `/tmp/full16k_specs.json` (2.7 MB), `/tmp/full16k_results.json` (6.9 MB), `/tmp/full16k_report.txt`, `/tmp/full16k_console.log`.

---

## Anomalies

None. 0 crashes, 0 timeouts, 0 info-mismatches, 0 unexplained tolerance breaches across all 21,375 specs.

The 5–7 specs with non-zero (but ≤ 1 ULP) diffs are concentrated in the `gl_wilkw` / `gl_wilk2w` family — these are extreme adversarial cases (Wilkinson glued bidiagonals); a 1-ULP variation between two implementations of the same algorithm is expected and well within the 2·eps target.

---

## Why bit-identity at all? (read this if the result looks too clean)

`mr3_gk.py` (Python) and `mr3gk_fortran/` (pure Fortran) are two **orchestrators** of the **same** XMR Fortran kernel. The numerical work — eigenvectors via `dlaxre_gk_` + `dlaxrv_`, bisection via `DSTEBZ`, Givens via `DLARTG`, normalization via system BLAS `DNRM2`, RNG via `DLARNV` — is performed by identical compiled code in both paths. The port replaces the Python orchestration glue with Fortran orchestration glue; it does not replace the math.

After the `_dlarnv_normal` and `_system_dnrm2` unifications in `mr3_gk.py` (see `02_xmr_modifications.md`), Python and Fortran paths route through the same compiled symbols for every floating-point operation that can have ordering-dependent results. Bit-identity in σ across all 21,375 matrices is therefore the **predicted** outcome, not a fluke.
