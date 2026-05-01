#!/bin/bash
# Build the instrumented `mr3gk_run_traced` executable, which emits a
# per-tree-node NCD trace via dlaxrf_traced.f.
#
# Strategy:
#   - Compile our dlaxrf_traced.f with the same flags as the original
#     XMR kernel objects (../fortran_objects/*.o).
#   - Link a traced shared library libxmr_traced.so by reusing every
#     object in ../fortran_objects/ EXCEPT dlaxrf.o, plus our
#     dlaxrf_traced.o.
#   - Build mr3gk_run_traced linking against libxmr_traced.so plus the
#     same orchestration objects (mr3gk_*.o) as production.
#
# The production libxmr.so and mr3gk_run remain untouched.

set -e
cd "$(dirname "$0")"

PARENT="$(cd .. && pwd)"
OBJDIR="${PARENT}/fortran_objects"

if [ ! -d "${OBJDIR}" ]; then
   echo "ERROR: ${OBJDIR} not found. Run ../build.sh first."
   exit 1
fi

FC="gfortran"

# Match the kernel build flags from ../build.sh
KERNEL_FFLAGS="-fPIC -O2 -std=legacy -w -fno-second-underscore"

# Match the orchestration build flags from ./build.sh
ORCH_FFLAGS="-O0 -fno-fast-math -fPIC -fno-second-underscore -fimplicit-none -Wno-unused -Wno-unused-dummy-argument"

echo "[1/3] Compiling dlaxrf_traced.f..."
$FC $KERNEL_FFLAGS -c dlaxrf_traced.f -o dlaxrf_traced.o

echo "[2/3] Compiling orchestration objects (if needed)..."
# Reuse the existing build.sh outputs if present; otherwise compile.
for f in mr3gk_consts mr3gk_split mr3gk_qrsweep mr3gk_utils mr3gk_tgk mr3gk_postproc mr3gk; do
   if [ ! -f "${f}.o" ]; then
      $FC $ORCH_FFLAGS -c ${f}.f90 -o ${f}.o
   fi
done

$FC $ORCH_FFLAGS -c mr3gk_run_traced.f90 -o mr3gk_run_traced.o

echo "[3/3] Linking mr3gk_run_traced..."
# Trick: link our dlaxrf_traced.o BEFORE libxmr.so on the command line.
# The static linker resolves dlaxrf_ from our object first; the .so version
# is shadowed but its other 30+ symbols still resolve normally.
KERNEL_OBJS=$(ls "${OBJDIR}"/*.o | grep -v -E "/dlaxrf\.o$" | tr '\n' ' ')
$FC -o mr3gk_run_traced mr3gk_run_traced.o \
    mr3gk.o mr3gk_postproc.o mr3gk_tgk.o mr3gk_utils.o \
    mr3gk_qrsweep.o mr3gk_split.o mr3gk_consts.o \
    dlaxrf_traced.o \
    ${KERNEL_OBJS} \
    -llapack -lblas -lgfortran -lquadmath -lm

echo
echo "Done: $(pwd)/mr3gk_run_traced"
echo "      ${PARENT}/libxmr_traced.so"
