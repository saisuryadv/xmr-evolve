#!/bin/bash
# Batch convert LAPACK Fortran 77 files to C using f2c
# Usage: ./convert_lapack.sh

F2C="/tmp/src/f2c"
LAPACK_SRC="/Users/saisurya/MRRR/bidiag-algo/lapack-double-src"
OUT_DIR="/Users/saisurya/MRRR/xmr-evolve/converted_fortran/lapack"
mkdir -p "$OUT_DIR"

# Copy f2c.h
cp /tmp/libf2c/f2c.h "$OUT_DIR/"

# Core MR³ eigensolver chain
MR3_FILES="dstemr dlarre dlarrv dlar1v dlarrf dlarrb dlarrc dlarra dlarrj dlasq2 dlasq3 dlasq4 dlasq5 dlasq6"

# Bidiag SVD drivers
BSVD_FILES="dbdsqr dbdsvdx dbdsdc"

# Utility routines needed by the above
UTIL_FILES="dlamch dlamc3 dlae2 dlaev2 dlaneg dlanst dlapy2 dlarnv dlascl dlasd0 dlasd1 dlasd2 dlasd3 dlasd4 dlasd5 dlasd6 dlasd7 dlasd8 dlasda dlasdq dlasdt dlasrt dlasrt2 dlas2 dlartg dlassq dlaneg ilaenv disnan dlaisnan dlaset dlasr dswap dcopy dscal drot dnrm2 ddot daxpy lsame xerbla ieeeck iparmq"

CONVERTED=0
FAILED=0

for group in "$MR3_FILES" "$BSVD_FILES" "$UTIL_FILES"; do
    for name in $group; do
        src="$LAPACK_SRC/${name}.f"
        if [ -f "$src" ]; then
            cd "$OUT_DIR"
            $F2C "$src" 2>/dev/null
            if [ -f "${name}.c" ]; then
                CONVERTED=$((CONVERTED + 1))
            else
                echo "FAILED: $name"
                FAILED=$((FAILED + 1))
            fi
        else
            # Try uppercase
            src="$LAPACK_SRC/${name^^}.f"
            if [ -f "$src" ]; then
                cd "$OUT_DIR"
                $F2C "$src" 2>/dev/null
                if [ -f "${name^^}.c" ]; then
                    mv "${name^^}.c" "${name}.c"
                    CONVERTED=$((CONVERTED + 1))
                fi
            fi
        fi
    done
done

echo ""
echo "Converted: $CONVERTED files"
echo "Failed: $FAILED files"
echo "Output in: $OUT_DIR"
ls "$OUT_DIR"/*.c 2>/dev/null | wc -l | xargs echo "Total .c files:"
