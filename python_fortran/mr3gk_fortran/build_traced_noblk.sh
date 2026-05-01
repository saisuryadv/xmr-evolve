#!/bin/bash
# Build the instrumented `mr3gk_run_traced_noblk` executable, which is
# identical to mr3gk_run_traced EXCEPT it uses USEBLOCKS = .FALSE.
# in dlaxrs_stat (every position becomes a scalar pivot, no 2x2 blocks).
#
# Same trace output (NCD_recon via N*G*N^T matmul in quad precision).

set -e
cd "$(dirname "$0")"

PARENT="$(cd .. && pwd)"
OBJDIR="${PARENT}/fortran_objects"

if [ ! -d "${OBJDIR}" ]; then
   echo "ERROR: ${OBJDIR} not found. Run ../build.sh first."
   exit 1
fi

FC="gfortran"
KERNEL_FFLAGS="-fPIC -O2 -std=legacy -w -fno-second-underscore"
ORCH_FFLAGS="-O0 -fno-fast-math -fPIC -fno-second-underscore -fimplicit-none -Wno-unused -Wno-unused-dummy-argument"

echo "[1/4] Compiling dlaxrf_traced.f..."
$FC $KERNEL_FFLAGS -c dlaxrf_traced.f -o dlaxrf_traced.o

echo "[2/4] Compiling dlaxrs_stat_noblk.f (USEBLOCKS=.FALSE.)..."
$FC $KERNEL_FFLAGS -c dlaxrs_stat_noblk.f -o dlaxrs_stat_noblk.o

echo "[3/4] Compiling orchestration objects (if needed)..."
for f in mr3gk_consts mr3gk_split mr3gk_qrsweep mr3gk_utils mr3gk_tgk mr3gk_postproc mr3gk; do
   if [ ! -f "${f}.o" ]; then
      $FC $ORCH_FFLAGS -c ${f}.f90 -o ${f}.o
   fi
done

# Re-use mr3gk_run_traced.f90 as-is (same CLI, same trace logging path).
$FC $ORCH_FFLAGS -c mr3gk_run_traced.f90 -o mr3gk_run_traced_noblk.o

echo "[4/4] Linking mr3gk_run_traced_noblk..."
# Override BOTH dlaxrf.o (for trace logging) AND dlaxrs_stat.o (for noblocks).
KERNEL_OBJS=$(ls "${OBJDIR}"/*.o | grep -v -E "/dlaxrf\.o$|/dlaxrs_stat\.o$" | tr '\n' ' ')
$FC -o mr3gk_run_traced_noblk mr3gk_run_traced_noblk.o \
    mr3gk.o mr3gk_postproc.o mr3gk_tgk.o mr3gk_utils.o \
    mr3gk_qrsweep.o mr3gk_split.o mr3gk_consts.o \
    dlaxrf_traced.o dlaxrs_stat_noblk.o \
    ${KERNEL_OBJS} \
    -llapack -lblas -lgfortran -lquadmath -lm

echo
echo "Done: $(pwd)/mr3gk_run_traced_noblk"
