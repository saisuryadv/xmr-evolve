# XMR-Evolve: Master Execution Plan

## Track A: Fortran → C++ Conversion
- Standard LAPACK (F77): Use f2c (/tmp/src/f2c) → generates .c files → compile as C++
- XMR files (F90 with INTENT): Strip F90 features → f2c, OR manual structured conversion
- hgbsvd files (dbdsgr, dlarri): Same as standard LAPACK (F77 compatible)
- DBDSVDX from LAPACK: Already in lapack-double-src/, convert via f2c
- Need: f2c.h header, libf2c library

## Track B: Evaluation Framework
- 19 STCollection matrices (already have .dat files)
- 22 adversarial patterns (from test_scaling_parallel.py)
- Paper-sourced bugs and matrices (agents researching)
- DBDSVDX bug from slides
- Online LAPACK bug reports

## Track C: Initial C++ Program (STEXR-based bidiagonal SVD)
- Based on Willems-Lang 2012 Algorithm 4.1
- Core: dstexr.f → C++ as tridiagonal eigensolver
- Wrapper: TGK matrix construction + U/V extraction
- hgbsvd code: dlarri (eigenvalue refinement for Großer-Lang)
- Test against Fortran for correctness

## Track D: Knowledge Organization for OpenEvolve
- Per-paper summaries with: contribution, bugs found, bugs in this paper, fixes
- Progression narrative linking papers chronologically
- Code documentation and call graphs
- Slide content extraction

## Track E: OpenEvolve + Claude Agent SDK Integration
- Install openevolve
- Integrate claude_agent_llm.py adapter
- config.yaml with cascading evaluation
- Prompt templates for code mutation, feedback, brainstorming
- LLM-driven test case discovery
