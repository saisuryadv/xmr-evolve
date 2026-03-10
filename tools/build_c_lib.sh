#!/bin/bash
# Build pure C library from CLAPACK + f2c-converted XMR + hgbsvd
# NO Fortran compiler needed. NO gfortran. NO Accelerate.
set -e

ROOT="/Users/saisurya/MRRR/xmr-evolve"
OUT_DIR="$ROOT/lib"
OBJ_DIR="$OUT_DIR/obj_c"

rm -rf "$OBJ_DIR"
mkdir -p "$OBJ_DIR"

CC="cc"
CFLAGS="-O2 -fPIC -w -I$ROOT/src/clapack"

echo "=== Compiling CLAPACK (LAPACK + BLAS) ==="
COMPILED=0
FAILED=0
for src in "$ROOT/src/clapack"/*.c; do
    name=$(basename "$src" .c)
    if $CC $CFLAGS -c "$src" -o "$OBJ_DIR/clapack_${name}.o" 2>/dev/null; then
        COMPILED=$((COMPILED + 1))
    else
        echo "  FAIL: $name"
        $CC $CFLAGS -c "$src" -o "$OBJ_DIR/clapack_${name}.o" 2>&1 | head -3
        FAILED=$((FAILED + 1))
    fi
done
echo "  CLAPACK: $COMPILED compiled, $FAILED failed"

echo "=== Compiling XMR (f2c-converted) ==="
COMPILED=0
FAILED=0
for src in "$ROOT/src/xmr_c"/*.c; do
    name=$(basename "$src" .c)
    if $CC $CFLAGS -I"$ROOT/src/xmr_c" -c "$src" -o "$OBJ_DIR/xmr_${name}.o" 2>/dev/null; then
        COMPILED=$((COMPILED + 1))
    else
        echo "  FAIL: $name"
        $CC $CFLAGS -I"$ROOT/src/xmr_c" -c "$src" -o "$OBJ_DIR/xmr_${name}.o" 2>&1 | head -3
        FAILED=$((FAILED + 1))
    fi
done
echo "  XMR: $COMPILED compiled, $FAILED failed"

echo "=== Compiling hgbsvd (f2c-converted) ==="
COMPILED=0
FAILED=0
for src in "$ROOT/src/hgbsvd_c"/*.c; do
    name=$(basename "$src" .c)
    if $CC $CFLAGS -c "$src" -o "$OBJ_DIR/hgbsvd_${name}.o" 2>/dev/null; then
        COMPILED=$((COMPILED + 1))
    else
        echo "  FAIL: $name"
        $CC $CFLAGS -c "$src" -o "$OBJ_DIR/hgbsvd_${name}.o" 2>&1 | head -3
        FAILED=$((FAILED + 1))
    fi
done
echo "  hgbsvd: $COMPILED compiled, $FAILED failed"

echo "=== Building static library ==="
cd "$OUT_DIR"
ar rcs libxmr_c.a obj_c/*.o
ranlib libxmr_c.a
echo "  Built: $OUT_DIR/libxmr_c.a ($(ls -lh libxmr_c.a | awk '{print $5}'))"
echo "  Objects: $(ls obj_c/*.o | wc -l | tr -d ' ') files"
echo ""
echo "Link with: -lxmr_c -lm (NO gfortran, NO Accelerate needed)"
