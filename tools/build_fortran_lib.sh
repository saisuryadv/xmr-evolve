#!/bin/bash
# Build all Fortran code (LAPACK + XMR + hgbsvd) into a static library
# that C++ can link against via extern "C"

set -e

LAPACK_SRC="/Users/saisurya/MRRR/bidiag-algo/lapack-double-src"
XMR_SRC="/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Code/tdsolver/xmr/SRC/O"
HGBSVD_SRC="/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Code/hgbsvd/hgbsvd/v2"
DBDSGR_SRC="/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Code/First SVD Attempt - dbdsgr.f"
OUT_DIR="/Users/saisurya/MRRR/xmr-evolve/lib"

mkdir -p "$OUT_DIR/obj"

FC="gfortran"
FFLAGS="-O2 -fPIC -std=legacy -w -fno-second-underscore"

echo "=== Compiling LAPACK routines ==="

# Core MR³ chain
LAPACK_FILES="dstemr dlarre dlarrv dlar1v dlarrf dlarrb dlarrc dlarra dlarrj
  dlasq2 dlasq3 dlasq4 dlasq5 dlasq6
  dbdsqr dbdsvdx dbdsdc
  dlamch dlae2 dlaev2 dlaneg dlanst dlapy2 dlarnv dlascl
  dlas2 dlartg dlassq dlasr dlaset dlasrt dlasda dlasdq dlasdt
  dlasd0 dlasd1 dlasd2 dlasd3 dlasd4 dlasd5 dlasd6 dlasd7 dlasd8
  dlacpy dlansp dlaruv dlaebz dlagts dlagtf
  dlacn2 dlaneg
  dswap dcopy dscal drot dnrm2 ddot daxpy dgemv dger dgemm
  dlabad dlacn2 dlagtf dlagts
  lsame xerbla ilaenv ieeeck iparmq disnan dlaisnan"

COMPILED=0
for name in $LAPACK_FILES; do
    src="$LAPACK_SRC/${name}.f"
    if [ -f "$src" ]; then
        $FC $FFLAGS -c "$src" -o "$OUT_DIR/obj/${name}.o" 2>/dev/null && COMPILED=$((COMPILED + 1)) || echo "  WARN: $name failed"
    fi
done
echo "  LAPACK: $COMPILED files compiled"

echo "=== Compiling XMR routines ==="
COMPILED=0
for src in "$XMR_SRC"/*.f; do
    name=$(basename "$src" .f)
    $FC $FFLAGS -c "$src" -o "$OUT_DIR/obj/xmr_${name}.o" 2>/dev/null && COMPILED=$((COMPILED + 1)) || echo "  WARN: $name failed"
done
echo "  XMR: $COMPILED files compiled"

echo "=== Compiling hgbsvd routines ==="
COMPILED=0
for src in "$HGBSVD_SRC"/*.f; do
    name=$(basename "$src" .f)
    $FC $FFLAGS -c "$src" -o "$OUT_DIR/obj/hgbsvd_${name}.o" 2>/dev/null && COMPILED=$((COMPILED + 1)) || echo "  WARN: $name failed"
done
echo "  hgbsvd: $COMPILED files compiled"

echo "=== Compiling dbdsgr (standalone) ==="
if [ -f "$DBDSGR_SRC" ]; then
    $FC $FFLAGS -c "$DBDSGR_SRC" -o "$OUT_DIR/obj/dbdsgr_standalone.o" 2>/dev/null && echo "  dbdsgr: compiled" || echo "  dbdsgr: FAILED"
fi

echo "=== Building static library ==="
cd "$OUT_DIR"
ar rcs libxmrlapack.a obj/*.o
ranlib libxmrlapack.a
echo "  Built: $OUT_DIR/libxmrlapack.a ($(ls -lh libxmrlapack.a | awk '{print $5}'))"
echo "  Objects: $(ls obj/*.o | wc -l | tr -d ' ') files"
