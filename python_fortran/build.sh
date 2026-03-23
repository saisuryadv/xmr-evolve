#!/bin/bash
# Build libxmr.so from Fortran objects + C wrapper
# Prerequisites: gfortran, gcc, liblapack-dev, libblas-dev
set -e

cd "$(dirname "$0")"

echo "Compiling dlaxre_gk.f (GK-enabled root representation)..."
gfortran -c -fPIC -O2 dlaxre_gk.f -o fortran_objects/dlaxre_gk.o

echo "Compiling xmr_wrapper.c..."
gcc -c -fPIC -O2 xmr_wrapper.c -o fortran_objects/xmr_wrapper.o

echo "Linking libxmr.so..."
gcc -shared -o libxmr.so \
  fortran_objects/dlaxrr.o \
  fortran_objects/dlaxrs_stat.o fortran_objects/dlaxrs_prog.o fortran_objects/dlaxrs.o \
  fortran_objects/dlaxrt_stat.o fortran_objects/dlaxrt_prog.o fortran_objects/dlaxrt.o \
  fortran_objects/dlaxrn_stat.o fortran_objects/dlaxrn.o fortran_objects/dlaxrn0.o \
  fortran_objects/dlaxrc.o fortran_objects/dlaxrg.o fortran_objects/dlaxrx.o \
  fortran_objects/dlaxrb_clssfy.o fortran_objects/dlaxrb_refcls.o fortran_objects/dlaxrb_refsng.o \
  fortran_objects/dlaxrf_selshf.o fortran_objects/dlaxrf_seltw_part.o fortran_objects/dlaxrf_seltw.o \
  fortran_objects/dlaxrf_cob.o fortran_objects/dlaxrf_iib.o \
  fortran_objects/dlaxrf_env.o fortran_objects/dlaxrf_grpenv.o \
  fortran_objects/dlaxrl_refine.o fortran_objects/dlaxrl_update.o fortran_objects/dlaxrl_reset.o \
  fortran_objects/dlaxrf.o \
  fortran_objects/dlaxrm_stat2.o fortran_objects/dlaxrm_stat4.o fortran_objects/dlaxrm_stat8.o \
  fortran_objects/dlaxrm_stat16.o fortran_objects/dlaxrm_stat32.o fortran_objects/dlaxrm_stat64.o \
  fortran_objects/dlaxrm.o \
  fortran_objects/dlaxre_initewldqds.o fortran_objects/dlaxre_gk.o \
  fortran_objects/dlaxrv.o \
  fortran_objects/xmr_wrapper.o \
  -lgfortran -llapack -lblas -lm

echo "Built libxmr.so ($(wc -c < libxmr.so) bytes)"
