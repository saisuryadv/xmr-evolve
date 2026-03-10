# MRRR Resources Index

All resources are in `../bidiag-algo/MRRR Resources/`.

## Papers (15 PDFs)

Located in `Papers/`:

| Year | Authors | Title | Key Contribution | Knowledge File |
|------|---------|-------|------------------|----------------|
| 1990 | Demmel, Kahan | Accurate Singular Values of Bidiagonal Matrices | Foundational bidiag SVD accuracy theory, zero-shift QR | — |
| 1997 | Dhillon | PhD Thesis | MR³ algorithm for symmetric tridiagonal eigenproblem | [mr3_foundations.md](mr3_foundations.md) |
| 2000 | Parlett, Dhillon | Relatively Robust Representations of Symmetric Tridiagonals | RRR theory, representation tree, relative gaps | [mr3_foundations.md](mr3_foundations.md) |
| 2001 | Großer, Lang | An O(n²) Algorithm for the Bidiagonal SVD | Coupling-based SVD via B^TB with coupled U/V recovery | [grosser_lang_2001_hgbsvd.md](grosser_lang_2001_hgbsvd.md) |
| 2002 | Barlow | More Accurate Bidiagonal Reduction for Computing the SVD | Bidiag reduction accuracy analysis | — |
| 2004 | Dhillon, Parlett | Orthogonal Eigenvectors and Relative Gaps | Relative gap theory, singleton/cluster classification | [mr3_foundations.md](mr3_foundations.md) |
| 2005 | — | LAPACK Working Note 163 | xSTEGR implementation details (MR³ in LAPACK) | — |
| 2005 | — | LAPACK Working Note 166 | xBDSCR proposal for LAPACK (bidiag SVD via coupling) | — |
| 2005 | — | On Symmetric Eigenproblems Induced by the Bidiagonal SVD | B^TB vs BB^T vs TGK reductions analysis | — |
| 2005 | Dhillon, Parlett, Vömel | Glued Matrices and the MRRR Algorithm | Glued Wilkinson matrices, MR³ stress tests | — |
| 2008 | Demmel, Marques, Parlett, Vömel | Performance and Accuracy of LAPACK's Symmetric Tridiagonal Eigensolvers | Test matrix taxonomy (S1/S2/W/G/U types), benchmark methodology | [demmel_2008_test_matrices.md](demmel_2008_test_matrices.md) |
| 2010 | Willems | PhD Thesis | XMR algorithm: block factorizations, Z-representations, envelope localization | [xmr_code_documentation.md](xmr_code_documentation.md) |
| 2012 | Willems, Lang | The MR³-GK Algorithm for the Bidiagonal SVD | Algorithm 4.1: NCD-aware MR³ on TGK. The target algorithm | [willems_lang_2012.md](willems_lang_2012.md) |
| 2013 | Willems, Lang | A Framework for the MR³ Algorithm: Theory and Implementation | Five-requirement MR³ framework, DSTEMR underflow bugs | [willems_lang_2013_framework.md](willems_lang_2013_framework.md) |
| 2020 | Marques, Demmel, Vasconcelos | Bidiagonal SVD Computation via an Associated Tridiagonal Eigenproblem | DBDSVDX bugs, CHKBD analysis, failure examples | [marques_2020_bugs_matrices.md](marques_2020_bugs_matrices.md) |

## Code

### Großer-Lang HGBSVD
`Code/hgbsvd/hgbsvd/v2/` — 24 Fortran files

Master routine: `dbdsgr.f` (coupling-based bidiag SVD). Supporting: `dlarrv.f` (eigenvector computation), `dlarrb.f` (bisection refinement), `dlarrf.f` (shift selection), `dlar1v.f` (twisted factorization), `dlasq1-6.f` (dqds eigenvalue computation), `dlacsv.f` (coupling transforms), `dlag2g.f` (GK construction), `dlap2s.f` (positive-definiteness), `dlas2p.f`, `dlats1-2.f` (test matrix generation).

Missing: `dlarri.f` (eigenvalue refinement by index) — referenced but not included in distribution. This causes INFO!=0 failures on some inputs.

### Willems XMR
`Code/tdsolver/xmr/SRC/O/` — 45 Fortran files (production variant)

Master routine: `dstexr.f` (improved MR³ eigensolver). Call graph:
```
dstexr → dlaxra (split/scale) → dlaxri (index map) → dlaxre (root rep)
       → dlaxrv (rep tree traversal)
           → dlaxrb_clssfy (classify singletons/clusters)
           → dlaxrb_refsng/refcls (refine)
           → dlaxrf_env (eigenvector envelope)
           → dlaxrf (child shift) → dlaxrs (blocked shift factorization)
           → dlaxrx (RQI + bisection for singletons)
       → dlaxro (sort eigenpairs)
```

Key innovations vs DSTEMR: block-aware 2×2 LDL^T, vectorized multi-shift negcount (2/4/8/16/32/64 unrolled via `dlaxrm.f`), envelope localization, Z-representations for children.

`Code/tdsolver/xmr/SRC/Dcl/` — 52 Fortran files (debug variant with logging)

### Other Code
- `Code/First SVD Attempt - dbdsgr.f` — Early version of Großer-Lang
- `Code/stegr_ID/` — DSTEGR identity/wrapper code
- `Code/Communication with Willems about XMR Failure.pdf` — Discussion of XMR failure cases
- `Code/Repository of test bi_tridiagonals.docx` — Test matrix documentation

## Test Data (STCollection)

`STCollection/DATA/` — Test matrices for bidiag SVD and symmetric tridiag eigensolvers:

### Bidiagonal Matrices (19 files)
`B_*.dat` files, sizes n=3 to n=429:
- `B_03.dat` through `B_Kimura_429.dat`
- Include splits, graded spectra, near-zero d entries, LAPACK bug triggers (bug316_gesdd, bug414)
- Used by evaluate.cpp as STCollection test set

### Symmetric Tridiagonal Matrices
`T_*.eig` and `T_*.dat` files — 20 tridiagonal test matrices:
- Godunov matrices, Wilkinson W21 with various glue values, bus/nasa/sts4098 from real applications
- Not used by bidiag SVD evaluator (different problem), but useful for testing tridiag eigensolvers

### Other Data
- `dist*.eig`, `Fann*.eig`, `sinc*.eig` — Additional eigenvalue distribution test cases
- Various `.jpeg` files — Log-scale eigenvalue distribution plots

### Evaluation Harness
- `STCollection/stetester/` — Fortran evaluation harness
- `STCollection/test_stetester.f90` — Test driver using DBDCHKRSLT metrics

## Presentations

Located in `Slides/`:

| File | Description |
|------|-------------|
| `BeBOP Talk: Revisiting MR^3 for the Bidiagonal SVD.pdf` | 88 slides by Ryan Schneider (UC Berkeley). Covers Option 1 vs Option 2 dilemma, failure examples (Marques 2020 graded matrix), MR³ history 1997-2020, algorithms compared. Ends with "Can we finally get MR³ working for the bidiag SVD?" |
| `Slides on Willems Thesis.pdf` | 16 slides by Beresford Parlett (UC Berkeley, 2010). XMR representation details, 9-case blocked shift formula, error analysis notation, Z/N/e representation comparison. |
| `Holy Grail.png` | Parody image — "Beresford Parlett and the Holy Grail." Represents the quest for reliable O(n²) bidiag SVD via MR³. |

## Prior Python Approaches

Located in `../bidiag-algo/`:

| File | Description |
|------|-------------|
| `approach_a_gk_split.py` | TGK with splitting |
| `approach_b_coupling.py` | B^TB/BB^T coupling |
| `approach_c_dc_hybrid.py` | Divide-and-conquer hybrid |
| `approach_d_split_coupling.py` | Split + coupling |
| `approach_e_hybrid.py` | Hybrid approach |
| `approach_f_tgk_stegr.py` | TGK + stemr |
| `approach_g_optimized.py` | Optimized variant |
| `approach_h_on2.py` | O(n²) attempt |
| `approach_i_willems.py` | Willems-inspired |
| `approach_j_on2.py` | O(n²) attempt #2 |
| `approach_k_ncd.py` | NCD-aware hybrid (best result: 19/19 STColl, all adversarial pass) |
| `approach_l_btb.py` | B^TB direct |
| `evaluate_alphaevolve.py` | Python evaluation framework (19 embedded matrices, metrics) |
| `mr3_tridiag.py` | Pure Python MR³ solver (~1050 lines, tinkerable reference) |

See [PRIOR_APPROACHES.md](PRIOR_APPROACHES.md) for detailed analysis of each approach's results and failure modes.
