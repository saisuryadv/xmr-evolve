# XMR Source Files

These 44 Fortran source files are from Christof Willems' XMR (eXtended MR³) implementation:

```
tdsolver/xmr/SRC/O/
```

Extracted from `tdsolver.zip` (Willems-Lang 2012 reference implementation).

## Compilation

```bash
gfortran -fPIC -O2 -std=legacy -w -fno-second-underscore -c <file>.f
```

## Modified files (not compiled from here)

- `dlaxrb_clssfy.f` → replaced by `../dlaxrb_clssfy_fix.f` (depth-0 classification fix)
- `dlaxre.f` → replaced by `../dlaxre_gk.f` (GK-aware root representation)

See `../build.sh` for the full build process.
