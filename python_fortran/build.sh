#!/bin/bash
# Build libxmr.so from Willems' XMR source (tdsolver/xmr/SRC/O) + our modifications.
#
# Prerequisites: gfortran, gcc, liblapack-dev, libblas-dev
#
# Source: xmr_src/        — original Willems XMR Fortran source
# Fixes:  dlaxrb_clssfy_fix.f — depth-0 classification bug fix
#         dlaxre_gk.f          — GK-aware root representation
#         xmr_wrapper.c        — C wrapper with debug callbacks
#
# Compilation flags must match original: -fPIC -O2 -std=legacy -w -fno-second-underscore
set -e

cd "$(dirname "$0")"

FC="gfortran"
FFLAGS="-fPIC -O2 -std=legacy -w -fno-second-underscore"
OBJDIR="fortran_objects"

mkdir -p "$OBJDIR"

# XMR source files compiled unmodified from xmr_src/
XMR_UNMODIFIED="
  dlaxrr dlaxrs_stat dlaxrs_prog dlaxrs
  dlaxrt_stat dlaxrt_prog dlaxrt
  dlaxrn_stat dlaxrn dlaxrn0
  dlaxrc dlaxrg dlaxrx
  dlaxrb_refcls dlaxrb_refsng
  dlaxrf_selshf dlaxrf_seltw_part dlaxrf_seltw
  dlaxrf_cob dlaxrf_iib
  dlaxrf_env dlaxrf_grpenv
  dlaxrl_refine dlaxrl_update dlaxrl_reset
  dlaxrf
  dlaxrm_stat2 dlaxrm_stat4 dlaxrm_stat8
  dlaxrm_stat16 dlaxrm_stat32 dlaxrm_stat64
  dlaxrm
  dlaxre_initewldqds
  dlaxrv
"

echo "Compiling XMR source files from xmr_src/..."
COMPILED=0
for name in $XMR_UNMODIFIED; do
    src="xmr_src/${name}.f"
    if [ -f "$src" ]; then
        $FC $FFLAGS -c "$src" -o "$OBJDIR/${name}.o"
        COMPILED=$((COMPILED + 1))
    else
        echo "  WARNING: $src not found"
    fi
done
echo "  Compiled $COMPILED unmodified XMR files"

echo "Compiling dlaxrb_clssfy_fix.f (depth-0 classification fix)..."
$FC $FFLAGS -c dlaxrb_clssfy_fix.f -o "$OBJDIR/dlaxrb_clssfy.o"

echo "Compiling dlaxre_gk.f (GK-aware root representation)..."
$FC -c -fPIC -O2 dlaxre_gk.f -o "$OBJDIR/dlaxre_gk.o"

echo "Compiling xmr_wrapper.c..."
gcc -c -fPIC -O2 xmr_wrapper.c -o "$OBJDIR/xmr_wrapper.o"

echo "Linking libxmr.so..."
gcc -shared -o libxmr.so \
  $OBJDIR/dlaxrr.o \
  $OBJDIR/dlaxrs_stat.o $OBJDIR/dlaxrs_prog.o $OBJDIR/dlaxrs.o \
  $OBJDIR/dlaxrt_stat.o $OBJDIR/dlaxrt_prog.o $OBJDIR/dlaxrt.o \
  $OBJDIR/dlaxrn_stat.o $OBJDIR/dlaxrn.o $OBJDIR/dlaxrn0.o \
  $OBJDIR/dlaxrc.o $OBJDIR/dlaxrg.o $OBJDIR/dlaxrx.o \
  $OBJDIR/dlaxrb_clssfy.o $OBJDIR/dlaxrb_refcls.o $OBJDIR/dlaxrb_refsng.o \
  $OBJDIR/dlaxrf_selshf.o $OBJDIR/dlaxrf_seltw_part.o $OBJDIR/dlaxrf_seltw.o \
  $OBJDIR/dlaxrf_cob.o $OBJDIR/dlaxrf_iib.o \
  $OBJDIR/dlaxrf_env.o $OBJDIR/dlaxrf_grpenv.o \
  $OBJDIR/dlaxrl_refine.o $OBJDIR/dlaxrl_update.o $OBJDIR/dlaxrl_reset.o \
  $OBJDIR/dlaxrf.o \
  $OBJDIR/dlaxrm_stat2.o $OBJDIR/dlaxrm_stat4.o $OBJDIR/dlaxrm_stat8.o \
  $OBJDIR/dlaxrm_stat16.o $OBJDIR/dlaxrm_stat32.o $OBJDIR/dlaxrm_stat64.o \
  $OBJDIR/dlaxrm.o \
  $OBJDIR/dlaxre_initewldqds.o $OBJDIR/dlaxre_gk.o \
  $OBJDIR/dlaxrv.o \
  $OBJDIR/xmr_wrapper.o \
  -lgfortran -llapack -lblas -lm

echo "Built libxmr.so ($(wc -c < libxmr.so) bytes)"
