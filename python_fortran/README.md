# MR³-GK Bidiagonal SVD — Python + Fortran Implementation

O(n²) bidiagonal SVD using MR³ with Golub-Kahan structure preservation (Willems & Lang 2012).

## Architecture

```
Python (mr3_gk.py)          Fortran (XMR, ~15K lines)
┌─────────────────┐         ┌──────────────────┐
│ bidiag_svd()    │         │ dlaxre_gk.f      │ ← GK-enabled root repr
│   split_bidiag  │         │ dlaxrv.f         │ ← eigenvector computation
│   mr3_tgk       │──ctypes─│ dlaxrf.f + deps  │ ← child repr (shift/twist)
│   sign recovery │         │ xmr_wrapper.c    │ ← C interface
│   Bv recovery   │         └──────────────────┘
│   GS completion │
└─────────────────┘
Eigenvalues: LAPACK dstebz (bisection)
```

## Current Results

| Suite | Tests | Passed | Rate |
|-------|-------|--------|------|
| Small (22 patterns × 3 sizes) | ~66 | ~53 | 80% |
| Full (90 patterns × 4 sizes + 19 STCollection) | 379 | 300 | 79% |

## Prerequisites

```bash
sudo apt install gfortran liblapack-dev libblas-dev python3-numpy
```

## Build

```bash
# Build libxmr.so from Fortran objects + modified dlaxre_gk.f
bash build.sh
```

A prebuilt `libxmr.so` is included for Linux x86_64.

## Run Evaluation

```bash
# Full 379-test suite (90 patterns × {10,100,200,400} + 19 STCollection)
python3 evaluate.py

# Small ~66-test suite (22 patterns × {10,50,100})
python3 evaluate.py --small

# Medium 270-test suite (90 patterns × {10,100,200})
python3 evaluate.py --medium
```

## Files

| File | Description |
|------|-------------|
| `mr3_gk.py` | Main implementation: splitting, T_GK construction, SVD extraction |
| `xmr_ctypes.py` | Python ctypes interface to XMR Fortran library |
| `full_eval.py` | 90 adversarial test matrix generators |
| `evaluate.py` | Evaluation runner with scoring |
| `dlaxre_gk.f` | Modified XMR root representation with GK detection enabled |
| `xmr_wrapper.c` | C wrapper for dlaxre_ + dlaxrv_ |
| `libxmr.so` | Prebuilt shared library (Linux x86_64) |
| `build.sh` | Build script to recompile from source |
| `fortran_objects/` | Precompiled XMR .o files (unmodified from Willems' code) |
| `stcollection/` | 19 STCollection bidiagonal test matrices (.dat files) |
| `CHANGES.md` | Full history of changes and design decisions |

## Key Design Decisions

- **No reorthogonalization**: All orthogonality comes from MR³ structure, not post-hoc correction
- **NEGL guard kept**: dlaxre's quality check routes some matrices to PD shift fallback — removing it causes regressions
- **Gram-Schmidt for zero σ only**: GS completion triggers at most once per unsplit block for zero singular values
- **Two-phase splitting**: Relative split for accuracy, absolute split per-block for condition (3.5)
