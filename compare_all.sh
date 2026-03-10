#!/bin/bash
# Compare all bidiagonal SVD variants
# Each variant is in src/bidiag_<name>.h, gets copied to src/bidiag_svd.h for compilation

set -e

CXX="g++"
CXXFLAGS="-std=c++17 -O2"
LDFLAGS="-Llib -lxmr -framework Accelerate -L/opt/homebrew/lib/gcc/15 -lgfortran -lm"
# DBDSVDX needs gfortran-compiled LAPACK (not Accelerate) to avoid string length convention mismatch
LDFLAGS_LAPACK="-Llib -lxmr -Llapack/build/lib -llapack -lblas -L/opt/homebrew/lib/gcc/15 -lgfortran -lm"

VARIANTS=("dbdsqr" "dbdsvdx" "tgk_stemr" "tgk_stexr")
LABELS=("DBDSQR" "DBDSVDX" "TGK+STEMR" "TGK+STEXR")

RESULTS_DIR="results"
mkdir -p "$RESULTS_DIR"

# Save original
cp src/bidiag_svd.h src/bidiag_svd.h.bak

for i in "${!VARIANTS[@]}"; do
    variant="${VARIANTS[$i]}"
    label="${LABELS[$i]}"
    header="src/bidiag_${variant}.h"
    outfile="$RESULTS_DIR/${variant}.txt"

    echo "=== Building $label ($header) ==="

    if [ ! -f "$header" ]; then
        echo "  SKIP: $header not found"
        continue
    fi

    cp "$header" src/bidiag_svd.h
    # DBDSVDX needs gfortran-compiled LAPACK + ASan to avoid heap corruption
    if [ "$variant" = "dbdsvdx" ]; then
        BUILD_LDFLAGS="$LDFLAGS_LAPACK"
        BUILD_CXXFLAGS="$CXXFLAGS -fsanitize=address"
    else
        BUILD_LDFLAGS="$LDFLAGS"
        BUILD_CXXFLAGS="$CXXFLAGS"
    fi
    if $CXX $BUILD_CXXFLAGS -o "evaluate_${variant}" src/evaluate.cpp $BUILD_LDFLAGS 2>"$RESULTS_DIR/${variant}_build.log"; then
        echo "  Running..."
        "./evaluate_${variant}" > "$outfile" 2>&1 || true
        echo "  Done -> $outfile"
    else
        echo "  BUILD FAILED (see $RESULTS_DIR/${variant}_build.log)"
        cat "$RESULTS_DIR/${variant}_build.log"
    fi
done

# Restore original
cp src/bidiag_svd.h.bak src/bidiag_svd.h
rm -f src/bidiag_svd.h.bak

echo ""
echo "============================================"
echo "         COMPARISON TABLE"
echo "============================================"
printf "%-12s  %7s  %8s  %8s  %8s  %6s  %6s\n" "Algorithm" "Pass" "AvgRes" "AvgOrtU" "AvgOrtV" "WScale" "Score"
echo "------------- -------  --------  --------  --------  ------  ------"

for i in "${!VARIANTS[@]}"; do
    variant="${VARIANTS[$i]}"
    label="${LABELS[$i]}"
    outfile="$RESULTS_DIR/${variant}.txt"

    if [ ! -f "$outfile" ]; then
        printf "%-12s  %7s\n" "$label" "N/A"
        continue
    fi

    pass=$(grep "^Pass:" "$outfile" 2>/dev/null | sed 's/Pass: //')
    avg_res=$(grep "^avg_residual=" "$outfile" 2>/dev/null | sed 's/avg_residual=//')
    avg_ortu=$(grep "^avg_ortho_u=" "$outfile" 2>/dev/null | sed 's/avg_ortho_u=//')
    avg_ortv=$(grep "^avg_ortho_v=" "$outfile" 2>/dev/null | sed 's/avg_ortho_v=//')
    worst_sc=$(grep "^worst_scaling_ratio=" "$outfile" 2>/dev/null | sed 's/worst_scaling_ratio=//')
    score=$(grep "^composite_score=" "$outfile" 2>/dev/null | sed 's/composite_score=//')

    printf "%-12s  %7s  %8s  %8s  %8s  %6s  %6s\n" \
        "$label" "$pass" "$avg_res" "$avg_ortu" "$avg_ortv" "$worst_sc" "$score"
done

# Cleanup binaries
rm -f evaluate_dbdsqr evaluate_dbdsvdx evaluate_tgk_stemr evaluate_tgk_stexr

echo ""
echo "Done."
