#!/bin/bash
# Build ONLY XMR + hgbsvd code into a static library.
# Standard LAPACK/BLAS comes from Apple Accelerate framework.

set -e

XMR_SRC="/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Code/tdsolver/xmr/SRC/O"
HGBSVD_SRC="/Users/saisurya/MRRR/bidiag-algo/MRRR Resources/Code/hgbsvd/hgbsvd/v2"
OUT_DIR="/Users/saisurya/MRRR/xmr-evolve/lib"

rm -rf "$OUT_DIR/obj"
mkdir -p "$OUT_DIR/obj"

FC="gfortran"
FFLAGS="-O2 -fPIC -std=legacy -w"
# NOTE: NO -fno-second-underscore so names match Accelerate's convention

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
    # Skip dmatgen.f and dtest.f (test infrastructure, not needed)
    if [ "$name" = "dmatgen" ] || [ "$name" = "dtest" ]; then
        continue
    fi
    $FC $FFLAGS -c "$src" -o "$OUT_DIR/obj/hgbsvd_${name}.o" 2>/dev/null && COMPILED=$((COMPILED + 1)) || echo "  WARN: $name failed"
done
echo "  hgbsvd: $COMPILED files compiled"

echo "=== Adding dsecnd stub ==="
cat > /tmp/dsecnd_stub.f << 'EOF'
      DOUBLE PRECISION FUNCTION DSECND()
      DSECND = 0.0D0
      RETURN
      END
EOF
$FC $FFLAGS -c /tmp/dsecnd_stub.f -o "$OUT_DIR/obj/dsecnd_stub.o"

echo "=== Building static library ==="
cd "$OUT_DIR"
ar rcs libxmr.a obj/*.o
ranlib libxmr.a
echo "  Built: $OUT_DIR/libxmr.a ($(ls -lh libxmr.a | awk '{print $5}'))"
echo "  Objects: $(ls obj/*.o | wc -l | tr -d ' ') files"
echo ""
echo "Link with: -lxmr -framework Accelerate -lgfortran"
