#!/usr/bin/env bash
# Reproducer: build reference LAPACK (if needed), then this offline PR.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LAPACK_REF="${LAPACK_REF:-/tmp/lapack-ref}"

# Enable devtoolset-11 if present (gfortran 4.8.x lacks ieee_arithmetic).
if [ -f /opt/rh/devtoolset-11/enable ]; then
    # shellcheck disable=SC1091
    source /opt/rh/devtoolset-11/enable
fi

if [ ! -f "${LAPACK_REF}/build/lib/liblapack.a" ] || \
   [ ! -f "${LAPACK_REF}/build/lib/libblas.a" ]; then
    echo "=== Building reference LAPACK at ${LAPACK_REF} ==="
    if [ ! -d "${LAPACK_REF}" ]; then
        git clone --depth 1 https://github.com/Reference-LAPACK/lapack.git \
            "${LAPACK_REF}"
    fi
    cmake -S "${LAPACK_REF}" -B "${LAPACK_REF}/build" \
          -DBUILD_INDEX_64_EXT_API=OFF \
          -DBUILD_INDEX_64=OFF \
          -DCBLAS=OFF -DLAPACKE=OFF \
          -DBUILD_TESTING=OFF \
          -DCMAKE_BUILD_TYPE=Release
    cmake --build "${LAPACK_REF}/build" --target lapack blas -j 4
fi

echo "=== Building offline DBDSVR PR ==="
cd "${HERE}"
make LAPACK_LIB="${LAPACK_REF}/build/lib/liblapack.a" \
     BLAS_LIB="${LAPACK_REF}/build/lib/libblas.a"

echo "=== Running canonical residual sweep ==="
./build/dbsvr_test < TESTING/svrtest.in
