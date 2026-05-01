#!/bin/bash
# Build the pure-Fortran mr3gk SVD library and CLI driver.
# Reuses the existing libxmr.so for the inner XMR Fortran routines.
#
# Output: mr3gk_run (CLI executable)

set -e
cd "$(dirname "$0")"

PARENT="$(cd .. && pwd)"
LIBXMR="${PARENT}/libxmr.so"

if [ ! -f "${LIBXMR}" ]; then
   echo "ERROR: ${LIBXMR} not found. Run ../build.sh first."
   exit 1
fi

FC="gfortran"
# -O0 to keep numerical operations in source order (helps bit-match Python)
FFLAGS="-O0 -fno-fast-math -fPIC -fno-second-underscore -fimplicit-none -Wno-unused -Wno-unused-dummy-argument"

OBJ=()
for f in mr3gk_consts mr3gk_split mr3gk_qrsweep mr3gk_utils mr3gk_tgk mr3gk_postproc mr3gk; do
   echo "Compiling ${f}.f90..."
   $FC $FFLAGS -c ${f}.f90 -o ${f}.o
   OBJ+=("${f}.o")
done

echo "Compiling mr3gk_run.f90..."
$FC $FFLAGS -c mr3gk_run.f90 -o mr3gk_run.o

echo "Linking mr3gk_run..."
$FC -o mr3gk_run mr3gk_run.o "${OBJ[@]}" "${LIBXMR}" \
    -llapack -lblas -lgfortran -lm

echo "Done: $(pwd)/mr3gk_run"
