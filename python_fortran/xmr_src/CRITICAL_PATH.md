# XMR Critical Path — which Fortran files are actually needed

Of the 44 `.f` files in Paul Willems' upstream `xmr_src/` (from
`tdsolver/xmr/SRC/O`), only **29 are on the critical path** of our
`libxmr.so`. We also add 2 patched files (`dlaxre_gk.f`,
`dlaxrb_clssfy_fix.f`) and a C wrapper (`xmr_wrapper.c`), giving
**32 compiled sources** total.

## What's compiled and why

### Root representation (2 files)
| File | Source | Purpose |
|---|---|---|
| `dlaxre_gk.f` | **Our patch** | GK-aware root rep builder (Bug #1 fix: activates `IF(.FALSE.)` branch) |
| `dlaxre_initewldqds.f` | upstream | Init eigenvalue list from dqds |

### Representation-tree forest (8 files)
| File | Purpose |
|---|---|
| `dlaxrf.f` | Main forest driver: find shift → build child rep |
| `dlaxrf_selshf.f` | Select candidate shifts |
| `dlaxrf_seltw.f` | Select valid twists |
| `dlaxrf_seltw_part.f` | Partial twist selection helper |
| `dlaxrf_cob.f` | Consecutive outer-bounds check |
| `dlaxrf_iib.f` | Interior inner-bounds check |
| `dlaxrf_env.f` | Envelope computation |
| `dlaxrf_grpenv.f` | Group-envelope computation |

### Cluster classification (4 files)
| File | Source | Purpose |
|---|---|---|
| `dlaxrb_clssfy_fix.f` | **Our patch** | Cluster classifier (Bug #2 fix: AVGTHRESH=0 guard) |
| `dlaxrb_refcls.f` | upstream | Refine classified eigenvalues |
| `dlaxrb_refsng.f` | upstream | Refine singleton eigenvalues |
| `dlaxrc.f` | upstream | Bisection for classified eigenvalues |

### Eigenvalue-list management (3 files)
| File | Purpose |
|---|---|
| `dlaxrl_refine.f` | Refine eigenvalue bounds |
| `dlaxrl_reset.f` | Reset eigenvalue list |
| `dlaxrl_update.f` | Update list after inertia call |

### qd-style operations (9 files)
| File | Purpose |
|---|---|
| `dlaxrr.f` | Initialize representation data structure from (G, Omega, E) |
| `dlaxrs.f` | Shift transformation dispatcher (calls stat + prog) |
| `dlaxrs_stat.f` | Stationary qd shift |
| `dlaxrs_prog.f` | Progressive qd shift |
| `dlaxrt.f` | Twisted factorization dispatcher |
| `dlaxrt_stat.f` | Stationary twisted factorization |
| `dlaxrt_prog.f` | Progressive twisted factorization |
| `dlaxrn.f` | Sturm count (inertia via twisted factorization) |
| `dlaxrn_stat.f` | Static inertia helper |

### Sturm-count batching (2 files)
| File | Purpose |
|---|---|
| `dlaxrm.f` | Multiple-shift Sturm count dispatcher |
| `dlaxrm_stat2.f` | Vectorized 2-shift Sturm count (the only variant used at runtime; `MAXPARNEG = 2`) |

### Eigenvectors (2 files)
| File | Purpose |
|---|---|
| `dlaxrv.f` | Main MRRR eigenvector computation (depth-first tree traversal) |
| `dlaxrx.f` | Singleton eigenvector via RQI + bisection |

### Misc (1 file)
| File | Purpose |
|---|---|
| `dlaxrg.f` | RQI eigenvector computation + Gershgorin discs |

### C glue (1 file)
| File | Purpose |
|---|---|
| `xmr_wrapper.c` | Exposes `xmr_eigenvectors`, `xmr_build_repr`, etc. for Python ctypes |

## What's NOT compiled (15 files)

### Never compiled by `build.sh` (8 upstream files)
These are part of Willems' public-API pipeline that we bypass:

| File | Purpose | Why we skip it |
|---|---|---|
| `dstexr.f` | Public eigenpair driver | We call `dlaxre_`/`dlaxrv_` directly |
| `dlaxra.f` | Block split / sign normalize | We do our own in `mr3_gk.py` / `mr3gk_fortran/` |
| `dlaxre.f` | Root rep builder (unpatched) | Replaced by `dlaxre_gk.f` |
| `dlaxri.f` | Locate block containing WIL:WIU | We handle indexing in our orchestration |
| `dlaxrk.f` | Pivoting / pinning for RRR | Not in call path |
| `dlaxro.f` | Reorder eigenvector output | We sort in our orchestration |
| `reptools.f` | Representation utilities | Not needed |
| `stats.f` | Statistics tracking | Not linked |

### Removed from `build.sh` (6 files, dead code)
These were previously compiled but never reached at runtime:

| File | Why dead |
|---|---|
| `dlaxrn0.f` | Never called by any reachable subroutine |
| `dlaxrm_stat4.f` | `dlaxrm.f` has `MAXPARNEG = 2` (compile-time parameter at line 151); only `dlaxrm_stat2` is entered |
| `dlaxrm_stat8.f` | Same — `IF(MAXPARNEG.GE.8)` at line 347 is dead |
| `dlaxrm_stat16.f` | Same — `IF(MAXPARNEG.GE.16)` at line 315 |
| `dlaxrm_stat32.f` | Same — `IF(MAXPARNEG.GE.32)` at line 283 |
| `dlaxrm_stat64.f` | Same — `IF(MAXPARNEG.GE.64)` at line 251 |

### `README.md` (1 file)
Provenance documentation — not Fortran source.

## Verification

After pruning `build.sh` to compile only the 29 critical upstream files
(+ 2 patches + 1 C wrapper = 32 total), `libxmr.so` shrinks from 184 KB
to 94 KB and all tests pass unchanged:

```
python3 evaluate.py              → TOTAL: 379/379 passed
python3 test_dense_to_bidiag.py  → GRAND TOTAL: 219/224 passed
python3 run_evaluate_fortran.py  → TOTAL: 379/379 passed (Fortran-only)
```

## Call-graph summary

```
xmr_wrapper.c (C entry points)
  ├─ dlaxre_gk (root rep)
  │    ├─ dlaxre_initewldqds
  │    ├─ dlaxrl_update
  │    └─ dlaxrr
  ├─ dlaxrv (eigenvectors)
  │    ├─ dlaxrb_clssfy_fix → dlaxrc → dlaxrl_refine, dlaxrm
  │    ├─ dlaxrb_refcls     → dlaxrc
  │    ├─ dlaxrb_refsng     → dlaxrl_refine, dlaxrm
  │    ├─ dlaxrf (forest)
  │    │    ├─ dlaxrf_selshf
  │    │    ├─ dlaxrf_seltw → dlaxrf_seltw_part
  │    │    ├─ dlaxrf_cob   → dlaxrl_reset, dlaxrl_update
  │    │    ├─ dlaxrf_iib   → dlaxrl_refine
  │    │    ├─ dlaxrr
  │    │    └─ dlaxrs → dlaxrs_stat, dlaxrs_prog
  │    ├─ dlaxrf_env → dlaxrf_grpenv → dlaxrg, dlaxrt
  │    ├─ dlaxrx → dlaxrg, dlaxrt
  │    └─ dlaxrr
  ├─ dlaxrn → dlaxrn_stat
  ├─ dlaxrm → dlaxrm_stat2
  ├─ dlaxrs → dlaxrs_stat, dlaxrs_prog
  ├─ dlaxrt → dlaxrt_stat, dlaxrt_prog
  └─ dlaxrr
```
