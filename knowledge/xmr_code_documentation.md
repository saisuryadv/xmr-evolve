# XMR Code Documentation: Willems' MR3 Implementation

Complete analysis of all 45 source files in:
`/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Code/tdsolver/xmr/SRC/O/`

---

## Table of Contents

1. [Overview](#overview)
2. [The Willems Bug and Fix](#the-willems-bug-and-fix)
3. [Files Grouped by Function Category](#files-grouped-by-function-category)
4. [Complete Call Graph](#complete-call-graph)
5. [Per-File Documentation](#per-file-documentation)
6. [Data Structures](#data-structures)
7. [Key Parameters and Tolerances](#key-parameters-and-tolerances)
8. [Safe vs Buggy Classification](#safe-vs-buggy-classification)

---

## Overview

XMR is Paul Willems' research implementation of the MR3 (Multiple Relatively Robust Representations) algorithm for computing selected eigenpairs of a real symmetric tridiagonal matrix. It was developed as part of Willems' PhD work under Bruno Lang at the University of Wuppertal.

The codebase consists of 44 Fortran 77/90 files and 1 C interface file (45 total). The top-level driver is `DSTEXR`, which computes eigenpairs WIL:WIU for a symmetric tridiagonal matrix T = diag(D) + diag(E,+-1).

**Key design features compared to LAPACK's DSTEMR:**
- Block-aware LDL^T factorizations (2x2 blocks for numerical stability)
- Compressed representation data structure (REPR/REPI arrays)
- Vectorized negcount (Sturm count) with unrolled variants for 2,4,8,16,32,64 simultaneous shifts
- Envelope localization (Parlett/Voemel strategy) for bounding eigenvector support
- Consistency guarantees for "naive" parallelization (non-overlapping index sets produce orthogonal vectors)
- DQDS integration for high relative accuracy eigenvalue initialization

**Total subroutine count:** ~55 subroutines/functions across the 45 files. Several files contain multiple routines (e.g., `dlaxrg.f` has DLAXRG0 and DLAXRG; `reptools.f` has XMR_INITREP, XMR_REPSIZE_REAL, XMR_REPSIZE_INT; `stats.f` has ~15 accessor functions).

---

## The Willems Bug and Fix

**Source:** Email communication between Osni Marques (LBNL), Bruno Lang, Paul Willems, and Jim Demmel, June 2014.

**Problem:** DSTEXR was producing non-orthogonal eigenvectors on certain STCollection matrices. LAPACK's DSTEGR worked correctly on the same inputs.

**Root cause:** In `dlaxrb_clssfy.f`, the "average gap recognition" feature used a tolerance that was too aggressive:
- Original: `AVGAPFAC = 0.1D0`, with `AVGTOL = MAX(SPDIAM, REPELG) * (AVGAPFAC / (N-1))`
- This caused eigenvalue clusters to be split too aggressively at shallow tree depths, leading to insufficient refinement and loss of orthogonality.

**Fix (two changes in `dlaxrb_clssfy.f`):**
1. Changed `AVGAPFAC = 0.1D0` to `AVGAPFAC = 0.3D0`
2. Changed AVGTOL formula to scale with DEPTH: `AVGTOL = MAX(SPDIAM, REPELG) * ((DEPTH * AVGAPFAC) / (N-1))`

**Note from Willems:** The fix was 4.5 years old at the time of the email but had been forgotten in the version distributed via Bruno Lang's link.

**Status in this codebase:** The fix IS included in the version analyzed here (AVGAPFAC=0.3D0 and DEPTH scaling are present in the code).

---

## Files Grouped by Function Category

### 1. Driver (1 file, 1 C wrapper)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dstexr.f` | DSTEXR | Top-level driver: splitting, root rep, rep tree, sorting |
| `dstexr-i.c` | dstexr_i() | C interface wrapper with workspace management |

### 2. Matrix Splitting & Scaling (1 file)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxra.f` | DLAXRA | Sign matrix, scaling, splitting into irreducible blocks |

### 3. Index Range Determination (2 files)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxri.f` | DLAXRI | Maps global WIL:WIU to per-block local indices |
| `dlaxrk.f` | DLAXRK | Bisection step for split matrix (Sturm count across blocks) |

### 4. Root Representation & Eigenvalue Initialization (2 files)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxre.f` | DLAXRE | Builds root representation and initializes eigenvalue list |
| `dlaxre_initewldqds.f` | DLAXRE_INITEWLDQDS | Initializes EW list from DQDS high-accuracy values |

### 5. Representation Tree Traversal (1 file)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxrv.f` | WSREQ_XRV, DLAXRV | Main MR3 depth-first loop: classify, refine, recurse |

### 6. Eigenvalue Classification & Refinement (4 files)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxrb_clssfy.f` | DLAXRB_CLSSFY | Classify eigenvalues as singletons/clusters **(BUGFIX FILE)** |
| `dlaxrb_refcls.f` | DLAXRB_REFCLS | Refine cluster bounds, reveal internal gaps |
| `dlaxrb_refsng.f` | DLAXRB_REFSNG | Refine singleton bounds via vectorized negcounts |
| `dlaxrc.f` | DLAXRC | Core bisection engine for eigenvalue refinement |

### 7. Child Representation (Shift Computation) (5 files)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxrf.f` | WSREQ_XRF, DLAXRF | Main child rep routine: find shift TAU for RRR |
| `dlaxrf_selshf.f` | DLAXRF_SELSHF | Generate and prioritize shift candidates |
| `dlaxrf_seltw.f` | DLAXRF_SELTW | Select twist index with optimal growth/condition |
| `dlaxrf_seltw_part.f` | DLAXRF_SELTW_PART | Partial twist evaluation (one direction) |
| `dlaxrf_cob.f` | DLAXRF_COB | Check outer bounds of child representation |
| `dlaxrf_iib.f` | DLAXRF_IIB | Initialize inner bounds from father's bounds |

### 8. Envelope Computation (2 files)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxrf_env.f` | DLAXRF_ENV | Determine eigenvector envelope for cluster |
| `dlaxrf_grpenv.f` | DLAXRF_GRPENV | Group envelope via Parlett/Voemel strategy |

### 9. Eigenvector Computation (2 files)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxrg.f` | DLAXRG0, DLAXRG | Compute FP-vector (zero/nonzero eigenvalue) |
| `dlaxrx.f` | WSREQ_XRX, DLAXRX | RQI with bisection backup for singletons |

### 10. Representation Data Structure (2 files)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxrr.f` | DLAXRR | Build derived quantities from primary LDL^T data |
| `reptools.f` | XMR_INITREP, XMR_REPSIZE_REAL, XMR_REPSIZE_INT | Convenience tools for rep data structure |

### 11. Negcount / Sturm Count (9 files)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxrn.f` | DLAXRN | Single negcount (returns xi = 2*negc + issing) |
| `dlaxrn0.f` | DLAXRN0 | Zero inertia from base D/R/OMEGA data |
| `dlaxrn_stat.f` | DLAXRN_STAT | Stationary negcount from twist to one end |
| `dlaxrm.f` | DLAXRM | Vectorized multi-negcount dispatcher |
| `dlaxrm_stat2.f` | DLAXRM_STAT2 | Unrolled negcount for 2 shifts |
| `dlaxrm_stat4.f` | DLAXRM_STAT4 | Unrolled negcount for 4 shifts |
| `dlaxrm_stat8.f` | DLAXRM_STAT8 | Unrolled negcount for 8 shifts |
| `dlaxrm_stat16.f` | DLAXRM_STAT16 | Unrolled negcount for 16 shifts |
| `dlaxrm_stat32.f` | DLAXRM_STAT32 | Unrolled negcount for 32 shifts |
| `dlaxrm_stat64.f` | DLAXRM_STAT64 | Unrolled negcount for 64 shifts |

### 12. Shift Factorization (3 files)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxrs.f` | DLAXRS | Blocked shift factorization (N G N' - tau) |
| `dlaxrs_prog.f` | DLAXRS_PROG | Progressive shift factorization |
| `dlaxrs_stat.f` | DLAXRS_STAT | Stationary shift factorization with full block handling |

### 13. Twisted Factorization (3 files)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxrt.f` | DLAXRT | Blocked-to-nonblocked twisted factorization |
| `dlaxrt_prog.f` | DLAXRT_PROG | Progressive twisted factorization |
| `dlaxrt_stat.f` | DLAXRT_STAT | Stationary twisted factorization |

### 14. Eigenvalue List (EWL) Operations (3 files)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxrl_refine.f` | DLAXRL_REFINE | Bisect interval using inertia |
| `dlaxrl_reset.f` | DLAXRL_RESET | Reset EWL to single interval |
| `dlaxrl_update.f` | DLAXRL_UPDATE | Incorporate negcount sample into EWL |

### 15. Sorting (1 file)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `dlaxro.f` | DLAXRO | Selection sort eigenpairs into ascending order |

### 16. Statistics Infrastructure (1 file)
| File | Subroutine(s) | Purpose |
|------|--------------|---------|
| `stats.f` | XMRS_INIT, XMRS_CHECK, XMRS_NNODES, XMRS_MAXDEPTH, XMRS_NUMFN, XMRS_NUMFT, XMRS_NUMGV, XMRS_NUMGV0, XMRS_NUMFS_2, XMRS_NUMFS_K, XMRS_NBIS_INIT, XMRS_NBIS_COB, XMRS_NBIS_IIB, XMRS_NBIS_CLASS, XMRS_NBIS_SNG, XMRS_NBIS_CLB | Statistics counters via COMMON /XMRSTATS/ |

---

## Complete Call Graph

```
DSTEXR (dstexr.f) — Top-level driver
|
|-- DLAXRA (dlaxra.f) — Split & scale
|   |-- DSCAL (LAPACK)
|   |-- DLAMCH (LAPACK)
|
|-- DLAXRI (dlaxri.f) — Index ranges per block
|   |-- DLAXRK (dlaxrk.f) — Bisection across blocks
|
|-- DLAXRE (dlaxre.f) — Root representation & EW init [per block]
|   |-- DLAXRN (dlaxrn.f) — Single negcount
|   |   |-- DLAXRN_STAT (dlaxrn_stat.f) — Stationary negcount half
|   |-- DLAXRR (dlaxrr.f) — Build representation data
|   |-- DLAXRE_INITEWLDQDS (dlaxre_initewldqds.f) — Init EWL from DQDS
|   |-- DLAXRL_UPDATE (dlaxrl_update.f) — Update EWL with sample
|   |   |-- DLAXRL_REFINE (dlaxrl_refine.f) — Bisect interval
|   |-- DLARNV (LAPACK) — Random perturbation
|   |-- DLASQ1 (LAPACK) — DQDS eigenvalue computation
|
|-- DLAXRV (dlaxrv.f) — Representation tree traversal [per block]
|   |
|   |-- DLAXRB_CLSSFY (dlaxrb_clssfy.f) — Classify singletons/clusters
|   |   |-- DLAXRC (dlaxrc.f) — Bisection refinement engine
|   |       |-- DLAXRL_REFINE (dlaxrl_refine.f)
|   |       |-- DLAXRM (dlaxrm.f) — Vectorized negcount
|   |           |-- DLAXRM_STAT2..64 (dlaxrm_stat{2,4,8,16,32,64}.f)
|   |           |-- DLAXRN (dlaxrn.f) — Fallback single negcount
|   |
|   |-- DLAXRB_REFSNG (dlaxrb_refsng.f) — Refine singletons
|   |   |-- DLAXRL_REFINE (dlaxrl_refine.f)
|   |   |-- DLAXRM (dlaxrm.f)
|   |
|   |-- DLAXRB_REFCLS (dlaxrb_refcls.f) — Refine clusters
|   |   |-- DLAXRC (dlaxrc.f)
|   |
|   |-- DLAXRF_ENV (dlaxrf_env.f) — Envelope computation
|   |   |-- DLAXRF_GRPENV (dlaxrf_grpenv.f) — Group envelope
|   |       |-- DLAXRT (dlaxrt.f) — Twisted factorization
|   |       |   |-- DLAXRT_STAT (dlaxrt_stat.f)
|   |       |   |-- DLAXRT_PROG (dlaxrt_prog.f)
|   |       |-- DLAXRG0 (dlaxrg.f) — FP-vector (zero eigenvalue)
|   |       |-- DLAXRG (dlaxrg.f) — FP-vector (nonzero eigenvalue)
|   |
|   |-- DLAXRF (dlaxrf.f) — Child representation (find shift)
|   |   |-- DLAXRF_SELSHF (dlaxrf_selshf.f) — Select shift candidates
|   |   |-- DLAXRF_SELTW (dlaxrf_seltw.f) — Select twist index
|   |   |   |-- DLAXRF_SELTW_PART (dlaxrf_seltw_part.f)
|   |   |-- DLAXRF_COB (dlaxrf_cob.f) — Check outer bounds
|   |   |   |-- DLAXRL_RESET (dlaxrl_reset.f)
|   |   |   |-- DLAXRL_UPDATE (dlaxrl_update.f)
|   |   |   |-- DLAXRN (dlaxrn.f)
|   |   |-- DLAXRF_IIB (dlaxrf_iib.f) — Init inner bounds
|   |   |   |-- DLAXRL_REFINE (dlaxrl_refine.f)
|   |   |   |-- DLAXRN (dlaxrn.f)
|   |   |-- DLAXRS (dlaxrs.f) — Shift factorization
|   |   |   |-- DLAXRS_STAT (dlaxrs_stat.f)
|   |   |   |-- DLAXRS_PROG (dlaxrs_prog.f)
|   |   |-- DLAXRR (dlaxrr.f)
|   |   |-- DLAXRN (dlaxrn.f)
|   |   |-- DLAXRL_REFINE (dlaxrl_refine.f)
|   |
|   |-- DLAXRX (dlaxrx.f) — RQI for singletons
|   |   |-- DLAXRT (dlaxrt.f) — Twisted factorization
|   |   |-- DLAXRG (dlaxrg.f) — FP-vector
|   |   |-- DLAXRG0 (dlaxrg.f) — FP-vector (zero)
|   |
|   |-- DLAXRR (dlaxrr.f) — Build child representation
|   |-- DLARNV (LAPACK) — Random perturbation for gap retry
|
|-- DLAXRO (dlaxro.f) — Sort eigenpairs
|   |-- DSWAP (LAPACK)

Utility/Infrastructure (not in main call chain):
|-- DLAXRN0 (dlaxrn0.f) — Zero inertia from base data
|-- reptools.f — XMR_INITREP, XMR_REPSIZE_REAL, XMR_REPSIZE_INT
|-- stats.f — XMRS_INIT, XMRS_CHECK, XMRS_NNODES, ... (15+ accessor functions)
|-- dstexr-i.c — C interface wrapper
```

---

## Per-File Documentation

### dstexr.f
**Subroutine:** `DSTEXR`
**Purpose:** Top-level driver for computing selected eigenpairs WIL:WIU of a symmetric tridiagonal matrix. Guarantees consistency for parallel (non-overlapping index sets produce orthogonal vectors).

**Parameters:**
- `N` (IN): Matrix size
- `D(N)` (INOUT): Diagonal elements (destroyed on output)
- `E(N-1)` (INOUT): Off-diagonal elements (destroyed on output)
- `WIL, WIU` (IN): Indices of wanted eigenvalues
- `M` (OUT): Number of eigenvalues found
- `W(N)` (OUT): Eigenvalues
- `Z(LDZ,N)` (OUT): Eigenvectors
- `ISUPPZ(2*N)` (OUT): Support of each eigenvector
- `WORK, IWORK` (workspace)
- `INFO` (OUT): Error code

**Key algorithmic steps:**
1. Call DLAXRA to split T into irreducible blocks
2. Call DLAXRI to determine per-block eigenvalue index ranges
3. For each block: call DLAXRE (root rep + EW init), then DLAXRV (rep tree)
4. Call DLAXRO to sort eigenpairs

**Key parameters:** QRDIM=0, GAPTOL=0.001, USEDQD=.TRUE.

**Calls:** DLAXRA, DLAXRI, DLAXRE, DLAXRV, DLAXRO, DLAMCH, DSCAL, XMRS_INIT

**TODO/HACK:** None.

---

### dstexr-i.c
**Function:** `dstexr_i()`
**Purpose:** C interface wrapper for Fortran DSTEXR. Handles workspace query, memory allocation, and Fortran calling convention.

**Parameters:** Same as DSTEXR but with C types and pointer passing.

**Calls:** DSTEXR (Fortran), malloc, free.

---

### dlaxra.f
**Subroutine:** `DLAXRA`
**Purpose:** Matrix splitting and scaling. Three steps: (1) sign matrix S to make E positive, (2) scale into numerical range, (3) split into irreducible blocks.

**Parameters:**
- `N` (IN): Matrix size
- `D(N)` (INOUT): Diagonal (modified by sign matrix and scaling)
- `E(N-1)` (INOUT): Off-diagonal (modified)
- `NBLCKS` (OUT): Number of irreducible blocks
- `ISPLIT(N)` (OUT): Block boundaries
- `SGNMTX(N)` (OUT): Sign matrix entries
- `SCLF` (OUT): Scaling factor

**Key algorithmic steps:**
1. Build sign matrix S such that S*T*S has positive off-diagonals
2. Scale matrix so that ||T||_inf is in [smlnum, bignum]
3. Split at positions where |e_i| < TOLFAC * N * EPS * max(|GL|, |GU|, GU-GL)

**Key parameters:** TOLFAC=3

**Calls:** DSCAL, DLAMCH

**TODO/HACK:** None.

---

### dlaxri.f
**Subroutine:** `DLAXRI`
**Purpose:** Maps global eigenpair indices WIL:WIU to per-block local indices. Computes Gershgorin bounds and value separators per block.

**Parameters:**
- `N, NBLCKS` (IN): Matrix size, number of blocks
- `D(N), E(N-1)` (IN): Matrix data
- `ISPLIT(N)` (IN): Block boundaries
- `WIL, WIU` (IN): Global eigenvalue indices
- `BWIL(NBLCKS), BWIU(NBLCKS)` (OUT): Per-block local indices
- `BGL(NBLCKS), BGU(NBLCKS)` (OUT): Per-block Gershgorin bounds
- `SVALS(WIL:WIU)` (OUT): Value separators

**Key algorithmic steps:**
1. Compute Gershgorin bounds for each block
2. If partial spectrum, use bisection (DLAXRK) to find value separators
3. Map global indices to per-block ranges

**Calls:** DLAXRK

**TODO/HACK:** None.

---

### dlaxrk.f
**Subroutine:** `DLAXRK`
**Purpose:** One bisection step using Sturm count across all blocks of a split matrix.

**Parameters:**
- `N, NBLCKS` (IN): Matrix size, block count
- `D(N), E(N-1), E2(N-1)` (IN): Matrix data
- `ISPLIT(N)` (IN): Block boundaries
- `GL, GU` (INOUT): Current bounds
- `BGL, BGU` (IN): Per-block bounds
- `WINDEX` (IN): Target eigenvalue index
- `PIVMIN` (IN): Pivot minimum

**Key algorithmic steps:**
1. Bisect at midpoint MID = (GL+GU)/2
2. Count eigenvalues <= MID across all blocks (skip blocks where MID is outside Gershgorin bounds)
3. Update GL or GU based on count vs WINDEX

**Calls:** None (uses inline Sturm count with PIVMIN).

**TODO/HACK:** None.

---

### dlaxre.f
**Subroutine:** `DLAXRE`
**Purpose:** Builds root representation (LDL^T factorization) and initializes the eigenvalue list for one irreducible block.

**Parameters:**
- `N` (IN): Block size
- `D(N), E(N-1)` (IN): Block diagonal/off-diagonal
- `IL, IU` (IN): Eigenvalue index range
- `WIL, WIU` (IN): Wanted eigenvalue range
- `GL, GU` (IN): Gershgorin bounds
- `REPR(4*N+3)` (OUT): Representation real data
- `REPI(6+N+N/2)` (OUT): Representation integer data
- `EWL_AE, EWL_LU` (OUT): Eigenvalue list
- `INFO` (OUT): Error code

**Key algorithmic steps:**
1. Try shifts at GL/GU quarter points; prefer side with more eigenvalues
2. Up to NTRYMAX=10 attempts to find semidefinite LDL^T factorization
3. Perturb root rep primary data by NULPROOTPERT=8 ULPs
4. If EMODE='d': use DLASQ1 (DQDS) for O(n^2) high-accuracy eigenvalues
5. If EMODE='o': use Gershgorin disc estimates (O(n))
6. Initialize eigenvalue list (EWL)

**Key parameters:** NTRYMAX=10, NULPROOTPERT=8, EMODE='d' (default)

**Calls:** DLAXRN, DLAXRR, DLAXRE_INITEWLDQDS, DLAXRL_UPDATE, DLARNV, DLASQ1

**TODO/HACK:** GK-form support is deactivated (commented out; caused problems with GRO synthetics).

---

### dlaxre_initewldqds.f
**Subroutine:** `DLAXRE_INITEWLDQDS`
**Purpose:** Initialize eigenvalue list from high-accuracy DQDS approximations.

**Parameters:**
- `N` (IN): Matrix size
- `REPR, REPI` (IN): Representation data
- `IL, IU` (IN): Eigenvalue index range
- `WIL, WIU` (IN): Wanted range
- `QDVALS(IL:IU)` (IN): DQDS eigenvalue approximations (ascending order)
- `EWL_AE(2*IL-1:2*IU)` (OUT): Eigenvalue list indices
- `EWL_LU(2*IL-1:2*IU)` (OUT): Eigenvalue list bounds

**Key algorithmic steps:**
1. Group eigenvalues by relative gap: gap > |max| * sqrt(PREC) separates groups
2. For each group, set per-eigenvalue bounds equal to the DQDS value
3. Relax outer bounds of groups by RELOFF = log(N) * PREC

**Note:** Currently only used with WIL=1, WIU=N (all eigenvalues wanted per block).

**Calls:** DLAMCH

**TODO/HACK:** Comment notes future versions without consistency guarantees might use partial feature.

---

### dlaxrv.f
**Subroutines:** `WSREQ_XRV` (workspace query), `DLAXRV` (main)
**Purpose:** Main MR3 loop: depth-first traversal of the representation tree.

**Parameters:**
- `N` (IN): Block size
- `REPR, REPI` (IN): Root representation
- `IL, IU, WIL, WIU` (IN): Eigenvalue ranges
- `EWL_AE, EWL_LU` (INOUT): Eigenvalue list
- `GAPTOL` (IN): Gap tolerance for classification
- `W(WIL:WIU)` (OUT): Eigenvalues
- `Z(LDZ,WIL:WIU)` (OUT): Eigenvectors
- `ISUPPZ(2*WIL-1:2*WIU)` (OUT): Eigenvector support
- `WORK, IWORK` (workspace)
- `INFO` (OUT): Error code

**Key algorithmic steps:**
1. At each node: classify eigenvalues (DLAXRB_CLSSFY)
2. Refine singletons (DLAXRB_REFSNG) then compute eigenvectors (DLAXRX)
3. For clusters: refine bounds (DLAXRB_REFCLS), compute envelope (DLAXRF_ENV), find child representation (DLAXRF)
4. Push child clusters onto stack, continue depth-first
5. MAXDEPTH=10 controls maximum tree depth (and workspace)

**Key parameters:** MAXDEPTH=10, NULPSHIFTPERT=0 (deactivated), NULPGRTPERT=16

**Calls:** DLAXRB_CLSSFY, DLAXRB_REFSNG, DLAXRB_REFCLS, DLAXRF_ENV, DLAXRF, DLAXRX, DLAXRR, DLARNV

**TODO/HACK:** None.

---

### dlaxrb_clssfy.f
**Subroutine:** `DLAXRB_CLSSFY`
**Purpose:** Refine eigenvalue bounds and classify them into singletons (isolated) and clusters (close together). **This is the file containing the Willems bug fix.**

**Parameters:**
- `N` (IN): Block size
- `REPR, REPI` (IN): Representation
- `IL, IU` (IN): Index range
- `GAPTOL` (IN): Gap tolerance (default 0.001)
- `DEPTH` (IN): Current tree depth
- `EWL_AE, EWL_LU, LGAP, UGAP` (INOUT): Eigenvalue list and gaps
- `NSINGLES, NCLUSTERS` (OUT): Classification counts
- `SINGLES, CLUSTERS` (OUT): Classification arrays
- `INFO` (OUT): Error code

**Key algorithmic steps:**
1. Refine eigenvalue bounds via bisection (DLAXRC)
2. Compute gaps between consecutive eigenvalue bounds
3. Apply gap tolerance: gap > GAPTOL * max(|lambda_i|, |lambda_{i+1}|) means "well-separated"
4. **Average gap check (BUGFIX):** AVGTOL = MAX(SPDIAM, REPELG) * (DEPTH * AVGAPFAC) / (N-1)
5. Classify: singletons have full gaps on both sides; otherwise cluster
6. Handle "fringes" - clusters intersecting but not contained in wanted range

**Key parameters:** AVGAPFAC=0.3D0 (was 0.1D0 -- BUG FIX), MAXAVGAPDEPTH=10

**Gap info constants:** GI_NOFULL=-2, GI_UNKNOWN=-1, GI_NOGAP=0, GI_INTGAP=1, GI_FULLGAP=2

**Calls:** DLAXRC

**TODO/HACK:** AVGAPFAC is the key bugfix parameter. The DEPTH scaling was the second part of the fix.

---

### dlaxrb_refcls.f
**Subroutine:** `DLAXRB_REFCLS`
**Purpose:** Refine eigenvalue bounds within a cluster to reveal internal gaps (subclusters).

**Parameters:**
- `N` (IN): Block size
- `REPR, REPI` (IN): Representation
- `IL, IU, CI, CJ` (IN): Index range and cluster extent
- `EWL_AE, EWL_LU, LGAP, UGAP` (INOUT): Eigenvalue list

**Key algorithmic steps:**
1. Use SCLTOL = 2^LOGSCLTOLFAC * sqrt(PREC) for subcluster gap tolerance
2. Refine via DLAXRC to reveal internal gaps
3. Extra in-cluster refinement is commented out (efficiency concerns for subset case)

**Calls:** DLAXRC

**TODO/HACK:** Extra in-cluster refinement commented out.

---

### dlaxrb_refsng.f
**Subroutine:** `DLAXRB_REFSNG`
**Purpose:** Refine singleton eigenvalue bounds after classification using vectorized negcounts.

**Parameters:**
- `N` (IN): Block size
- `REPR, REPI` (IN): Representation
- `IL, IU` (IN): Index range
- `NSINGLES` (IN): Number of singletons
- `SINGLES` (IN): Singleton indices
- `EWL_AE, EWL_LU, LGAP, UGAP` (INOUT): Eigenvalue list

**Key algorithmic steps:**
1. Set accuracy targets: GAPACCSNG = 0.01/N, RELACCSNG = 0.01/N
2. Collect all singleton midpoints into a batch
3. Compute vectorized negcounts via DLAXRM
4. Refine bounds via DLAXRL_REFINE

**Calls:** DLAXRL_REFINE, DLAXRM

**TODO/HACK:** None.

---

### dlaxrc.f
**Subroutine:** `DLAXRC`
**Purpose:** Core bisection refinement engine. Refines eigenvalues by index to relative/absolute tolerance.

**Parameters:**
- `N` (IN): Block size
- `REPR, REPI` (IN): Representation
- `IL, IU, CI, CJ` (IN): Ranges
- `RELACCFAC, ABSACCFAC` (IN): Accuracy factors
- `MAXITS` (IN): Max bisection iterations
- `EWL_AE, EWL_LU, LGAP, UGAP` (INOUT): Eigenvalue list

**Key algorithmic steps:**
1. Loop over intervals needing refinement
2. For each: bisect at midpoint, compute negcount (via DLAXRM)
3. Update bounds via DLAXRL_REFINE
4. Handle duplicate indices and already-processed intervals efficiently

**Calls:** DLAXRL_REFINE, DLAXRM

**TODO/HACK:** None.

---

### dlaxrf.f
**Subroutines:** `WSREQ_XRF` (workspace query), `DLAXRF` (main, ~1000 lines)
**Purpose:** Central routine for finding shift TAU such that the child representation T - TAU = N_p G_p N_p^T is a Relatively Robust Representation (RRR).

**Parameters:**
- `N` (IN): Block size
- `REPR, REPI` (IN): Father representation
- `IL, IU, CI, CJ` (IN): Ranges and cluster extent
- `GAPTOL` (IN): Gap tolerance
- `EWL_AE, EWL_LU, LGAP, UGAP` (INOUT): Father's eigenvalue list
- `SONREPR, SONREPI` (OUT): Child representation
- `SONEWL_AE, SONEWL_LU, SONLGAP, SONUGAP` (OUT): Child's eigenvalue list
- `SONIL, SONIU, SONCJ, SONCI` (OUT): Child's ranges
- `WORK, IWORK` (workspace)
- `INFO` (OUT): Error code

**Key algorithmic steps:**
1. Generate shift candidates (DLAXRF_SELSHF) at outside/inside locations
2. For each candidate: compute shifted factorization (DLAXRS)
3. Select optimal twist index (DLAXRF_SELTW) with minimal element growth
4. Check outer bounds (DLAXRF_COB) and inner bounds (DLAXRF_IIB)
5. Evaluate quality (EVALOK=5.0 threshold); accept or try next candidate
6. FORCEGAP=.TRUE.: require child eigenvalues to be well-separated enough
7. MAXCPL_OUTSIDE=3, MAXCPL_INSIDE=2: max candidates per location

**Key parameters:** EVALOK=5.0, FORCEGAP=.TRUE., FUDGEFAC=2, MAXFUDGE=1, MAXCPL_OUTSIDE=3, MAXCPL_INSIDE=2

**Calls:** DLAXRS, DLAXRF_SELSHF, DLAXRF_SELTW, DLAXRF_COB, DLAXRF_IIB, DLAXRR, DLAXRN, DLAXRL_REFINE

**TODO/HACK:** Line 279: "TODO, see dlaxrt" regarding block ending at source twist.

---

### dlaxrf_selshf.f
**Subroutine:** `DLAXRF_SELSHF`
**Purpose:** Generate and prioritize shift candidates around a cluster.

**Parameters:**
- `IL, IU, CI, CJ` (IN): Index range and cluster extent
- `EWL_AE, EWL_LU, LGAP, UGAP` (IN): Eigenvalue list
- `MODE` (IN): 1=both outside/inside, 2=outside only, 3=inside only
- `MAXCPL` (IN): Max candidates per location
- `NCAND` (OUT): Number of candidates generated
- `CANDVAL, CANDLOC, CANDDIR, CANDFLG` (OUT): Candidate arrays

**Key algorithmic steps:**
1. Classify locations as primary (flag=1), secondary (flag=2), ternary (flag=3)
2. Inside shifts require min partition MININPRIMPART=0.4 of cluster
3. Last resort: throw shift in middle of tight cluster (DOLASTRESORTTHROWIN=.TRUE.)
4. Organize into batches and compress

**Key parameters:** MININPRIMPART=0.4, DOLASTRESORTTHROWIN=.TRUE.

**Calls:** None (pure computation).

**TODO/HACK:** Lines 179-181: "prefer side closer to min/max ew (Parlett)" (not yet implemented).

---

### dlaxrf_seltw.f
**Subroutine:** `DLAXRF_SELTW`
**Purpose:** Select twist index with optimal element growth and relative condition number.

**Parameters:**
- `N` (IN): Block size
- Various factorization arrays (IN)
- `TWISTOK(N)` (IN): Per-twist feasibility flags from DLAXRS
- `BESTTWIST` (OUT): Optimal twist index
- `BESTRC` (OUT): Relative condition at best twist
- `BESTGROWTH` (OUT): Element growth at best twist

**Key algorithmic steps:**
1. Compute per-twist inertias from two sweeps (top-to-bottom, bottom-to-top)
2. Set tiny gammas (< EPS*|TAU|) to zero for consistency
3. Accept twist a priori if MAXRELCOND=10 and MAXGROWTH=8 thresholds met
4. Special case: definite shifts always pass
5. If no a priori acceptance, evaluate all twists via DLAXRF_SELTW_PART

**Key parameters:** MAXRELCOND=10, MAXGROWTH=8

**Calls:** DLAXRF_SELTW_PART

**TODO/HACK:** None.

---

### dlaxrf_seltw_part.f
**Subroutine:** `DLAXRF_SELTW_PART`
**Purpose:** Compute relative condition (RC) and element growth estimates for one direction of a partial twist evaluation.

**Parameters:**
- Various factorization arrays (IN)
- Direction, twist index (IN)
- `ARCSQ, ABETA, AGMAX, AGSMAX, ABRCFAC` (OUT): Growth/condition arrays

**Key algorithmic steps:**
1. Sweep from one end towards twist
2. Accumulate ARCSQ (sum of squares for RC), AGMAX (max element growth)
3. Handle block structure boundaries

**Calls:** None (pure computation).

**TODO/HACK:** None.

---

### dlaxrf_cob.f
**Subroutine:** `DLAXRF_COB`
**Purpose:** Check outer bounds of child representation: verify no eigenvalue underflow, check left/right bounds and outer gaps.

**Parameters:**
- `N` (IN): Block size
- `SONREPR, SONREPI` (IN): Child representation
- `CI, CJ` (IN): Cluster extent
- `EWL_LU, LGAP, UGAP` (IN): Father's bounds
- `SONEWL_AE, SONEWL_LU, SONLGAP, SONUGAP` (OUT): Child's EWL
- `STATUS` (OUT): Result code

**Status codes:**
- 0: All checks passed
- +/-1: Eigenvalue underflow
- +/-2: Two very small eigenvalues
- +/-3: Cannot verify bound
- +/-5: Cannot verify gap

**Key parameters:** LOGFBNDRELAXFAC=7 (bounds relaxed by 2^7 * N * EPS)

**Calls:** DLAXRL_RESET, DLAXRL_UPDATE, DLAXRN

**TODO/HACK:** None.

---

### dlaxrf_iib.f
**Subroutine:** `DLAXRF_IIB`
**Purpose:** Initialize inner bounds of child's eigenvalue list from father's bounds.

**Parameters:**
- Father and child representation data
- Father's EWL bounds
- Child's EWL (OUT)

**Key algorithmic steps:**
1. Map father's inner bounds to child's coordinate system (subtract shift)
2. Verify via negcount (DLAXRN)
3. Three relaxation levels: 2*EPS, BISACCFAC*N*EPS, 2^LOGFBNDRELAXFAC*N*EPS
4. Refine via DLAXRL_REFINE

**Calls:** DLAXRL_REFINE, DLAXRN

**TODO/HACK:** None.

---

### dlaxrf_env.f
**Subroutine:** `DLAXRF_ENV`
**Purpose:** Determine eigenvector envelope (support bounds) for a cluster using Parlett/Voemel ideas.

**Parameters:**
- `N` (IN): Block size
- `REPR, REPI` (IN): Representation
- `CI, CJ` (IN): Cluster extent
- `EWL_LU, LGAP, UGAP` (IN): Eigenvalue list
- `ENV(N)` (OUT): Envelope (per-component bound)
- `MODE` (IN): full, outer only, or improve outer with inner

**Key algorithmic steps:**
1. Check prerequisites: ENVGAPFAC=100 (cluster width < mingap/100), cluster < N/3 eigenvalues
2. Call DLAXRF_GRPENV for subclusters
3. Three modes: full computation, outer bounds only, improve outer with inner data

**Key parameters:** ENVGAPFAC=100

**Calls:** DLAXRF_GRPENV

**TODO/HACK:** None.

---

### dlaxrf_grpenv.f
**Subroutine:** `DLAXRF_GRPENV`
**Purpose:** Compute envelope for a subcluster using two twisted factorizations at cluster endpoints.

**Parameters:**
- `N` (IN): Block size
- `REPR, REPI` (IN): Representation
- Subcluster bounds, eigenvalue data (IN)
- `ENV(N)` (INOUT): Envelope (updated)

**Key algorithmic steps:**
1. For singletons: compute eigenvector approximation + sin bound
2. For subgroups: Parlett/Voemel strategy with TAUFAC=2 backoff from cluster edges
3. Compute twisted factorizations at both endpoints
4. Combine to get envelope bounds per component
5. Force minimal envelope entry: MAXRELCOND * N * EPS / relgap

**Key parameters:** TAUFAC=2, MAXRELCOND (from dlaxrf_seltw)

**Calls:** DLAXRT, DLAXRG0, DLAXRG

**TODO/HACK:** None.

---

### dlaxrg.f
**Subroutines:** `DLAXRG0`, `DLAXRG`
**Purpose:**
- DLAXRG0: Compute FP-vector for zero eigenvalue using block-aware back-substitution
- DLAXRG: Compute FP-vector for nonzero eigenvalue via twisted factorization back-substitution

**Parameters (DLAXRG):**
- `N` (IN): Block size
- Factorization data (DPLUS, RPLUS, etc.) (IN)
- `TWIST` (IN): Twist index
- `Z(N)` (OUT): Normalized eigenvector approximation
- `NORMZ` (OUT): Norm of Z before normalization
- `RESID` (OUT): Residual = |gamma|/norm
- `RQCORR` (OUT): Rayleigh quotient correction = gamma/normsq

**Key algorithmic steps:**
1. Start from twist index with z(twist) = 1
2. Solve upper triangular system upward, lower triangular downward
3. Handle creeping underflow and support cuts
4. Normalize and compute residual/RQ correction

**Calls:** None (pure computation).

**TODO/HACK:** None.

---

### dlaxrx.f
**Subroutines:** `WSREQ_XRX` (workspace query), `DLAXRX`
**Purpose:** Rayleigh Quotient Iteration (RQI) with bisection backup for computing singleton eigenvectors.

**Parameters:**
- `N` (IN): Block size
- `REPR, REPI` (IN): Representation
- `LAMBDA` (INOUT): Eigenvalue approximation (refined on output)
- `Z(N)` (OUT): Eigenvector
- `RESID` (OUT): Residual norm
- `ISUPPZ(2)` (OUT): Support
- `WORK, IWORK` (workspace)
- `INFO` (OUT): Error code

**Key algorithmic steps:**
1. Compute twisted factorization at lambda (DLAXRT)
2. Compute FP-vector (DLAXRG/DLAXRG0)
3. If residual decreasing by RESDECFAC=0.7: apply RQ correction and repeat
4. If residual stalls: switch to bisection burst (PBBLEN=4 steps)
5. Track best vector via double-buffered storage
6. Special case: lambda=0 uses DLAXRG0 directly

**Key parameters:** ACCTOLFAC=2, RESTOLFAC=2, CUTTOLFAC=0.5, PBBLEN=4, RESDECFAC=0.7

**Calls:** DLAXRT, DLAXRG, DLAXRG0

**TODO/HACK:** None.

---

### dlaxrr.f
**Subroutine:** `DLAXRR`
**Purpose:** Build derived quantities from primary LDL^T data (G and Omega arrays).

**Parameters:**
- `N` (IN): Block size
- `K` (IN): Twist index
- `TYPE` (IN): Representation type
- `E(N-1)` (IN): Off-diagonal elements
- `PIVBASE` (IN): Pivot base for numerical stability
- `REPR(4*N+3)` (INOUT): Real data (G at offset 0, derived data filled)
- `REPI(6+N+N/2)` (INOUT): Integer data (block boundaries, etc.)

**Data layout:**
- `REPR`: IXG=0 (G values), IXBDET=N (block determinants), IXNGN=2N+1 (e^2/g), IXGNSQ=3N+2 (e^2), PIVBASE at 4N+3
- `REPI`: TYPE at 1, TWIST at 2, NB at 3, OMEGA at 4, LBBEGK at N+6

**Computed derived quantities:**
- GNSQ(i) = e_i^2
- NGN(i) = e_i^2 / g_i
- BDET(i) = block determinants (for 2x2 blocks)
- LBBEGK(j) = block boundary positions

**Calls:** None (pure computation).

**TODO/HACK:** None.

---

### reptools.f
**Subroutines:** `XMR_INITREP`, `XMR_REPSIZE_REAL`, `XMR_REPSIZE_INT`
**Purpose:** Convenience tools for representation data structure, primarily for C++ integration.

- `XMR_INITREP(N, K, TYPE, G, OMEGA, E, PIVBASE, REPR, REPI)`: Set up representation from explicitly given primary data. Copies G and OMEGA, then calls DLAXRR.
- `XMR_REPSIZE_REAL(N)`: Returns 4*N+3 (real array size)
- `XMR_REPSIZE_INT(N)`: Returns 6+N+N/2 (integer array size)

**Calls:** DLAXRR

**TODO/HACK:** None.

---

### dlaxrn.f
**Function:** `DLAXRN` (returns INTEGER)
**Purpose:** Single negcount (Sturm count). Returns xi = 2*negc + issing, encoding both the negative count and whether TAU is an eigenvalue.

**Parameters:**
- `N` (IN): Block size
- `REPR(4*N+3), REPI(6+N+N/2)` (IN): Representation
- `TAU` (IN): Shift value
- Returns: xi (integer encoding inertia)

**Key algorithmic steps:**
1. Special handling for TAU=0 with block-aware counting
2. General case: call DLAXRN_STAT for top half (1..K) and bottom half (K..N)
3. Combine: gamma = G(K) + (aux_top + aux_bottom) - TAU
4. xi = 2*(negc_top + negc_bottom) + {0 if gamma>0, 1 if gamma==0, 2 if gamma<0}

**Calls:** DLAXRN_STAT

**TODO/HACK:** None.

---

### dlaxrn0.f
**Function:** `DLAXRN0`
**Purpose:** Specialized zero inertia computation from raw D, R, OMEGA arrays (before representation is fully built).

**Parameters:**
- `N` (IN): Matrix size
- `D(N), R(N)` (IN): Diagonal and off-diagonal data
- `OMEGA(N)` (IN): Sign array

**Calls:** None.

**TODO/HACK:** None.

---

### dlaxrn_stat.f
**Subroutine:** `DLAXRN_STAT`
**Purpose:** Compute stationary part of negcount from twist to one end of the matrix.

**Parameters:**
- `N, K` (IN): Block size, twist index
- `DIR` (IN): Direction (+1 or -1)
- `G, GNSQ, NGN, BDET` (IN): Representation data
- `IBB` (INOUT): Block boundary index
- `LBBEGK` (IN): Block boundary positions
- `PIVMIN, TAU` (IN): Pivot minimum, shift
- `NEGC` (OUT): Negative count
- `AUX` (OUT): Auxiliary accumulator

**Key algorithmic steps:**
1. Sweep from end toward twist K
2. At each step: accumulate aux = aux - tau, compute g + aux, count negatives
3. Handle block boundaries with R_brk computation (PRBRK=8)

**Key parameters:** PRBRK=8

**Calls:** None (pure computation).

**TODO/HACK:** None.

---

### dlaxrm.f
**Subroutine:** `DLAXRM`
**Purpose:** Vectorized multi-shift negcount dispatcher. Computes Sturm counts for M shifts simultaneously against the same representation.

**Parameters:**
- `N` (IN): Block size
- `REPR, REPI` (IN): Representation
- `M` (IN): Number of shifts
- `ATAU(M)` (IN): Array of shift values
- `AXI(M)` (OUT): Array of inertias (xi values)

**Key algorithmic steps:**
1. Process shifts in batches of 64, 32, 16, 8, 4, 2 using unrolled routines
2. Each batch: call DLAXRM_STATXX twice (top-to-bottom, bottom-to-top)
3. Combine results at twist: gamma = G(K) + (aux_top + aux_bottom) - tau
4. Remaining shifts (< 2): use single DLAXRN
5. MAXPARNEG=2 controls maximum batch size used (conservative default)

**Key parameters:** MAXPARNEG=2 (only uses DLAXRM_STAT2 by default)

**Calls:** DLAXRM_STAT{2,4,8,16,32,64}, DLAXRN

**TODO/HACK:** None.

---

### dlaxrm_stat2.f through dlaxrm_stat64.f (6 files)
**Subroutines:** `DLAXRM_STAT2`, `DLAXRM_STAT4`, `DLAXRM_STAT8`, `DLAXRM_STAT16`, `DLAXRM_STAT32`, `DLAXRM_STAT64`
**Purpose:** Manually unrolled stationary negcount for 2/4/8/16/32/64 shifts simultaneously. These are performance-critical inner loops.

**Parameters (all identical pattern):**
- `N, K, DIR` (IN): Block size, twist, direction
- `G, GNSQ, NGN, BDET` (IN): Representation data arrays
- `IBB` (INOUT): Block boundary index
- `LBBEGK` (IN): Block boundary positions
- `PIVMIN` (IN): Pivot minimum
- `TAU(XX)` (IN): Array of XX shift values
- `ANEGC(XX)` (OUT): Negative counts
- `AAUX(XX)` (OUT): Auxiliary accumulators

**Key algorithmic steps (same for all, just XX-way unrolled):**
1. Sweep from end toward twist
2. For each matrix element, perform XX independent negcount updates in parallel (register-level)
3. Handle block boundaries with PRBRK=8 branching
4. Two branches at block boundaries: "broken" branch (AUX small or sign change) vs. "continuation" branch

**Key parameters:** PRBRK=8 (all variants)

**Calls:** None (pure computation, innermost loops).

**TODO/HACK:** None. These are auto-generated/templated code.

---

### dlaxrs.f
**Subroutine:** `DLAXRS`
**Purpose:** Blocked shift factorization: compute N_p G_p N_p^T = N G N^T - TAU with block structure tracking.

**Parameters:**
- `N` (IN): Block size
- `REPR, REPI` (IN): Source representation
- `TAU` (IN): Shift value
- `DPLUS(N)` (OUT): Top-to-bottom factorization diagonal
- `RPLUS(N)` (OUT): Bottom-to-top factorization diagonal
- `GAMMAP(N)` (OUT): Gamma values at each twist
- `TWISTOK(N)` (OUT): Feasibility flags per twist
- `OMGADP(N), OMGARP(N)` (OUT): Omega arrays for shifted factorization

**Key algorithmic steps:**
1. Compute progressive factorization top-to-bottom (DPLUS, via DLAXRS_PROG)
2. Compute progressive factorization bottom-to-top (RPLUS, via DLAXRS_PROG)
3. Handle block structure changes during shift
4. At each position: determine if twist is feasible (TWISTOK)

**Calls:** DLAXRS_STAT, DLAXRS_PROG

**TODO/HACK:** Line 279: "TODO, see dlaxrt" regarding block ending at source twist.

---

### dlaxrs_prog.f
**Subroutine:** `DLAXRS_PROG`
**Purpose:** Progressive shift factorization from J1 to J2 with shift perturbation SHFPRT. Non-blocked implementation only.

**Parameters:**
- `N, J1, J2` (IN): Block size, sweep range
- Factorization data (IN/OUT)
- `SHFPRT` (IN): Shift perturbation

**Calls:** None (pure computation).

**TODO/HACK:** None.

---

### dlaxrs_stat.f
**Subroutine:** `DLAXRS_STAT`
**Purpose:** Stationary shift factorization with full block structure handling. Most complex of the shift factorization routines.

**Parameters:**
- `N, K, DIR` (IN): Block size, twist, direction
- Source representation data (IN)
- `TAU` (IN): Shift
- Various output arrays for factorization data

**Key algorithmic steps:**
1. Sweep from end toward twist K
2. At each step: create/break/overlap blocks as needed
3. Block creation: when accumulator is small relative to diagonal
4. Block breaking: when block determinant condition fails
5. Overlap control: PMAXOSEQLEN=1 (max overlap sequence length)

**Key parameters:** KBLOCK=1/8, PK1=KBLOCK*0.999, PK2=PK1/3.01, PMAXOSEQLEN=1, PRBRK=5, PRCRE=5, PROSQ=0.25

**Calls:** None (pure computation).

**TODO/HACK:** None.

---

### dlaxrt.f
**Subroutine:** `DLAXRT`
**Purpose:** Twisted factorization from blocked to non-blocked form for all twists in J1:J2.

**Parameters:**
- `N, J1, J2` (IN): Block size, twist range
- `DPLUS, RPLUS` (IN): From DLAXRS
- Source representation data (IN)
- `GAMMA(J1:J2)` (OUT): Gamma values
- `D(N), L(N)` (OUT): Non-blocked LDL^T factors

**Key algorithmic steps:**
1. For each twist in J1:J2: combine DPLUS (top) and RPLUS (bottom)
2. Special handling for twists where source has block-end
3. Method 1 (primary) and fallback for block-end twists
4. Method 2 is commented out ("something is wrong with method 2")

**Calls:** DLAXRT_STAT, DLAXRT_PROG

**TODO/HACK:** "something is wrong with method 2" (commented out alternative).

---

### dlaxrt_prog.f
**Subroutine:** `DLAXRT_PROG`
**Purpose:** Progressive twisted factorization from J1 to J2.

**Parameters:**
- `N, J1, J2` (IN): Block size, sweep range
- Various factorization arrays (IN/OUT)

**Key parameters:** PRBRK=8

**Calls:** None (pure computation).

**TODO/HACK:** None.

---

### dlaxrt_stat.f
**Subroutine:** `DLAXRT_STAT`
**Purpose:** Stationary twisted factorization from end to twist J. Simpler than dlaxrs_stat (no block creation, only breaking).

**Parameters:**
- `N, J, DIR` (IN): Block size, twist, direction
- Various factorization arrays (IN/OUT)

**Key parameters:** PRBRK=8

**Calls:** None (pure computation).

**TODO/HACK:** None.

---

### dlaxrl_refine.f
**Subroutine:** `DLAXRL_REFINE`
**Purpose:** Bisect an interval [I,J] in the eigenvalue list at lambda using inertia xi.

**Parameters:**
- `IL, IU` (IN): Global index range
- `I, J` (IN): Current interval indices
- `LAMBDA` (IN): Bisection point
- `XI` (IN): Inertia at lambda
- `EWL_AE, EWL_LU, LGAP, UGAP` (INOUT): Eigenvalue list

**Key algorithmic steps:**
1. Determine which side of the interval lambda falls on (using xi)
2. Update bounds accordingly
3. Handle non-monotonicity in Sturm counts
4. Update LGAP/UGAP if boundary eigenvalues affected

**Calls:** None (pure computation).

**TODO/HACK:** None.

---

### dlaxrl_reset.f
**Subroutine:** `DLAXRL_RESET`
**Purpose:** Reset eigenvalue list to a single interval [L, U] with given inertias.

**Parameters:**
- `IL, IU` (IN): Index range
- `L, U` (IN): Lower and upper bounds
- `XIL, XIU` (IN): Inertias at L and U
- `EWL_AE, EWL_LU, LGAP, UGAP` (INOUT): Eigenvalue list

**Key algorithmic steps:**
1. Set all eigenvalue bounds to [L, U]
2. Set all interval indices to [IL, IU]
3. Handle odd inertias at boundaries (singleton splitting)

**Calls:** None (pure computation).

**TODO/HACK:** None.

---

### dlaxrl_update.f
**Subroutine:** `DLAXRL_UPDATE`
**Purpose:** Incorporate a negcount sample (LAMBDA, XI) into the eigenvalue list.

**Parameters:**
- `IL, IU` (IN): Index range
- `LAMBDA` (IN): Sample point
- `XI` (IN): Inertia at sample point
- `EWL_AE, EWL_LU, LGAP, UGAP` (INOUT): Eigenvalue list

**Key algorithmic steps:**
1. Look for intervals containing LAMBDA (using XI to speed up search)
2. Refine matching intervals via DLAXRL_REFINE
3. Update outer gaps (LGAP/UGAP) if sample falls outside all intervals

**Calls:** DLAXRL_REFINE

**TODO/HACK:** None.

---

### dlaxro.f
**Subroutine:** `DLAXRO`
**Purpose:** Selection sort eigenpairs (W, Z, ISUPPZ) into ascending eigenvalue order.

**Parameters:**
- `N, M, LDZ` (IN): Matrix size, number of eigenpairs, leading dimension
- `W(M)` (INOUT): Eigenvalues (sorted on output)
- `Z(LDZ,M)` (INOUT): Eigenvectors (reordered on output)
- `ISUPPZ(2*M)` (INOUT): Support (reordered on output)
- `REVORD(M)` (OUT): Reverse permutation

**Key algorithmic steps:**
1. Standard selection sort: find minimum in W(I:M), swap with position I
2. Swap corresponding columns of Z and entries of ISUPPZ
3. Track reverse permutation in REVORD

**Calls:** DSWAP (LAPACK)

**TODO/HACK:** None.

---

### stats.f
**Subroutines:** `XMRS_INIT`, `XMRS_CHECK`, and 14 accessor functions (`XMRS_NNODES`, `XMRS_MAXDEPTH`, `XMRS_NUMFN`, `XMRS_NUMFT`, `XMRS_NUMGV`, `XMRS_NUMGV0`, `XMRS_NUMFS_2`, `XMRS_NUMFS_K`, `XMRS_NBIS_INIT`, `XMRS_NBIS_COB`, `XMRS_NBIS_IIB`, `XMRS_NBIS_CLASS`, `XMRS_NBIS_SNG`, `XMRS_NBIS_CLB`)
**Purpose:** Statistics infrastructure using COMMON /XMRSTATS/ block. Tracks operation counts, tree depth, bisection steps, RQI iterations, etc.

**COMMON block contents:**
- `XTIME1, XTIME2, XTIME3`: Times for stages 1-3
- `XDDDDD`: Data integrity marker (double, set to 666.66)
- `XNBLCKS`: Number of blocks
- `XNNODES, XMAXDEPTH`: Tree node count, max depth
- `XNUMFN, XNUMFT, XNUMGV, XNUMGV0`: Call counts (negcount, twisted fac, gv, gv0)
- `XNUMFS_2, XNUMFS_K`: Shift factorization counts
- `XNBIS_INIT, XNBIS_COB, XNBIS_IIB, XNBIS_CLASS, XNBIS_SNG, XNBIS_CLB`: Bisection step counts
- `XNBIS_WASTED`: Wasted bisection steps
- `XNRQI, XNRQIBIS, XMAXNRQI, XMAXNRQIBIS`: RQI counts
- `XNENVGV, XNENVTF`: Envelope computation counts
- `XIIIII`: Data integrity marker (integer, set to -1234567)
- `XSTEALTHMODE`: When true, suppresses counting in internal routines

**Note:** Comment warns "SYNCHRONIZE ANY CHANGES HERE WITH xmr.h"

**Calls:** None.

**TODO/HACK:** None.

---

## Data Structures

### Representation (REPR / REPI)

The compressed representation stores an LDL^T factorization with 2x2 block structure.

**REPR (real array, size 4*N+3):**
```
Offset 0:      G(1:N)        -- Diagonal of D in N G N^T (primary data)
Offset N:      BDET(1:N)     -- Block determinants (derived)
Offset 2N+1:   NGN(1:N)      -- e_i^2 / g_i (derived)
Offset 3N+2:   GNSQ(1:N)     -- e_i^2 (derived)
Position 4N+3: PIVBASE       -- Pivot base for numerical stability
```

**REPI (integer array, size 6+N+N/2):**
```
Position 1:    TYPE          -- Representation type
Position 2:    TWIST (K)     -- Twist index
Position 3:    NB            -- Number of blocks
Position 4:    OMEGA(1:N)    -- Sign array (+1/-1)
Offset N+6:    LBBEGK(1:NB+1) -- Block boundary positions
```

### Eigenvalue List (EWL)

The eigenvalue list stores bracketing intervals for eigenvalues IL through IU.

**EWL_LU (real array, indices 2*IL-1 : 2*IU):**
```
EWL_LU(2*i-1) = lower bound for eigenvalue i
EWL_LU(2*i)   = upper bound for eigenvalue i
```

**EWL_AE (integer array, indices 2*IL-1 : 2*IU):**
```
EWL_AE(2*i-1) = leftmost eigenvalue index sharing this lower bound
EWL_AE(2*i)   = rightmost eigenvalue index sharing this upper bound
```
When EWL_AE(2*i-1) == EWL_AE(2*i) == i, eigenvalue i has its own isolated interval (singleton).

**LGAP, UGAP (scalars):**
- LGAP: gap to the left of eigenvalue IL (distance to next eigenvalue below)
- UGAP: gap to the right of eigenvalue IU (distance to next eigenvalue above)

---

## Key Parameters and Tolerances

| Parameter | Value | File | Purpose |
|-----------|-------|------|---------|
| GAPTOL | 0.001 | dstexr.f | Gap tolerance for singleton classification |
| AVGAPFAC | 0.3 | dlaxrb_clssfy.f | Average gap factor **(BUG FIX: was 0.1)** |
| MAXAVGAPDEPTH | 10 | dlaxrb_clssfy.f | Max depth for AVGTOL scaling |
| TOLFAC | 3 | dlaxra.f | Splitting tolerance factor |
| NTRYMAX | 10 | dlaxre.f | Max attempts for root representation |
| NULPROOTPERT | 8 | dlaxre.f | ULPs for root rep perturbation |
| MAXDEPTH | 10 | dlaxrv.f | Max representation tree depth |
| NULPGRTPERT | 16 | dlaxrv.f | Gap retry perturbation ULPs |
| EVALOK | 5.0 | dlaxrf.f | Child rep quality threshold |
| MAXRELCOND | 10 | dlaxrf_seltw.f | Max relative condition for twist acceptance |
| MAXGROWTH | 8 | dlaxrf_seltw.f | Max element growth for twist acceptance |
| ENVGAPFAC | 100 | dlaxrf_env.f | Envelope gap requirement |
| LOGFBNDRELAXFAC | 7 | dlaxrf_cob.f | Bound relaxation (2^7 * N * EPS) |
| MAXCPL_OUTSIDE | 3 | dlaxrf.f | Max shift candidates per outside location |
| MAXCPL_INSIDE | 2 | dlaxrf.f | Max shift candidates per inside location |
| MININPRIMPART | 0.4 | dlaxrf_selshf.f | Min partition for inside primary shifts |
| ACCTOLFAC | 2 | dlaxrx.f | RQI eigenvalue accuracy tolerance |
| RESTOLFAC | 2 | dlaxrx.f | RQI residual tolerance |
| CUTTOLFAC | 0.5 | dlaxrx.f | RQI support cutting tolerance |
| PBBLEN | 4 | dlaxrx.f | Bisection burst length in RQI |
| RESDECFAC | 0.7 | dlaxrx.f | Required residual decrease per RQI step |
| MAXPARNEG | 2 | dlaxrm.f | Max parallel negcounts (conservative) |
| PRBRK | 8 | dlaxrn_stat.f, dlaxrm_stat*.f | Block breaking branch parameter |
| KBLOCK | 1/8 | dlaxrs_stat.f | Block creation threshold |

---

## Safe vs Buggy Classification

### Known Buggy (Fixed)
| File | Issue | Status |
|------|-------|--------|
| `dlaxrb_clssfy.f` | AVGAPFAC too aggressive (0.1 instead of 0.3), missing DEPTH scaling in AVGTOL | **FIXED** in this version |

### Contains TODO / Known Issues
| File | Issue | Risk |
|------|-------|------|
| `dlaxrf.f` | "TODO, see dlaxrt" on line 279 (block ending at source twist) | LOW - edge case |
| `dlaxrf_selshf.f` | "prefer side closer to min/max ew (Parlett)" not implemented (lines 179-181) | LOW - optimization only |
| `dlaxrs.f` | "TODO, see dlaxrt" on line 279 (same issue as dlaxrf.f) | LOW - edge case |
| `dlaxrt.f` | "something is wrong with method 2" (commented out) | NONE - fallback in use |
| `dlaxre.f` | GK-form support deactivated (caused problems) | NONE - deactivated |
| `dlaxrb_refcls.f` | Extra in-cluster refinement commented out | LOW - efficiency tradeoff |

### Safe (No Known Issues)
All other 38 files are clean with no TODO/BUG/FIXME/HACK comments and no known issues:

`dstexr.f`, `dstexr-i.c`, `dlaxra.f`, `dlaxri.f`, `dlaxrk.f`, `dlaxre_initewldqds.f`, `dlaxrv.f`, `dlaxrb_refsng.f`, `dlaxrc.f`, `dlaxrf_seltw.f`, `dlaxrf_seltw_part.f`, `dlaxrf_cob.f`, `dlaxrf_iib.f`, `dlaxrf_env.f`, `dlaxrf_grpenv.f`, `dlaxrg.f`, `dlaxrx.f`, `dlaxrr.f`, `reptools.f`, `dlaxrn.f`, `dlaxrn0.f`, `dlaxrn_stat.f`, `dlaxrm.f`, `dlaxrm_stat2.f`, `dlaxrm_stat4.f`, `dlaxrm_stat8.f`, `dlaxrm_stat16.f`, `dlaxrm_stat32.f`, `dlaxrm_stat64.f`, `dlaxrs_prog.f`, `dlaxrs_stat.f`, `dlaxrt_prog.f`, `dlaxrt_stat.f`, `dlaxrl_refine.f`, `dlaxrl_reset.f`, `dlaxrl_update.f`, `dlaxro.f`, `stats.f`

---

## Summary Statistics

- **Total files:** 45 (44 Fortran + 1 C)
- **Total subroutines/functions:** ~55
- **Lines of code (approx):** ~12,000 Fortran + ~100 C
- **Largest file:** `dlaxrf.f` (~1000 lines), `dlaxrm_stat16.f` (~590 lines)
- **Known bugs:** 1 (fixed: AVGAPFAC in dlaxrb_clssfy.f)
- **TODOs:** 4 (all low risk)
- **LAPACK dependencies:** DLAMCH, DSCAL, DSWAP, DLARNV, DLASQ1
- **Key innovation vs LAPACK DSTEMR:** Block-aware LDL^T factorizations with 2x2 blocks, vectorized multi-shift negcount, envelope localization, consistency guarantees for parallelization
