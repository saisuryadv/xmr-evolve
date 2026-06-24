#!/bin/sh
# Build BidiagonalSVD_TGK with the allocatable STCollection driver.
# Mirrors run_tests.sh but only produces test_stcoll_alloc (the pract runner).
set -e
cd "$(dirname "$0")"
LAPACK="${LAPACK:--llapack -lblas}"
echo "Compiling core objects (LAPACK flags: $LAPACK) ..."
gfortran -c -O2 -std=legacy -fno-automatic \
    dbdtgk.f dlar1v_tgk.f dlarrf_tgk.f dlarrv_tgk.f dbdsvdmr3.f \
    stegr_ID/dlarrb.f
gfortran -O2 -std=legacy -fno-automatic \
    dev/test_stcoll_alloc.f *.o $LAPACK -o test_stcoll_alloc
gfortran -O2 -std=legacy -fno-automatic \
    dev/dbdsqr_ref.f $LAPACK -o dbdsqr_ref
gfortran -O2 -std=legacy -fno-automatic \
    dev/test_dbdsqr_full.f $LAPACK -o test_dbdsqr_full
rm -f *.o
echo "Built: $(pwd)/test_stcoll_alloc"
echo "Built: $(pwd)/dbdsqr_ref"
echo "Built: $(pwd)/test_dbdsqr_full"
