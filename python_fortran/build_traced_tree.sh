#!/bin/bash
# Build traced and no-additive libraries for Godunov tree visualization.
#
# Produces:
#   libxmr_traced.so   — additive fix + dlaxrv/dlaxrx debug tracing to stderr
#   libxmr_noadd.so    — original multiplicative perturbation (no additive fix)
#
# The production libxmr.so is NOT touched.

set -e
cd "$(dirname "$0")"

OBJDIR="fortran_objects"
FC="gfortran"
KFLAGS="-fPIC -O2 -std=legacy -w -fno-second-underscore"

if [ ! -d "$OBJDIR" ]; then
    echo "ERROR: $OBJDIR not found. Run bash build.sh first."
    exit 1
fi

echo "=== Building libxmr_traced.so ==="
# Compile traced versions (with dbg_ calls)
$FC $KFLAGS -c xmr_src/dlaxrv_traced.f -o /tmp/dlaxrv_traced.o
$FC $KFLAGS -c xmr_src/dlaxrx_traced.f -o /tmp/dlaxrx_traced.o
# Rebuild xmr_wrapper.c (has dbg_singleton_rqi_)
gcc -fPIC -O2 -c xmr_wrapper.c -o /tmp/xmr_wrapper_traced.o

# Link: all objects except dlaxrv.o, dlaxrx.o, xmr_wrapper.o + traced versions
OBJS=$(ls $OBJDIR/*.o | grep -v -E "dlaxrv\.o$|dlaxrx\.o$|xmr_wrapper\.o$" | tr '\n' ' ')
$FC -shared -o libxmr_traced.so \
    /tmp/dlaxrv_traced.o /tmp/dlaxrx_traced.o /tmp/xmr_wrapper_traced.o \
    $OBJS -llapack -lblas
echo "  Built libxmr_traced.so ($(stat -c%s libxmr_traced.so) bytes)"

echo "=== Building libxmr_noadd.so ==="
# Compile original dlaxre (pure multiplicative, no additive fix)
$FC $KFLAGS -c initial_code/xmr_src/dlaxre.f -o /tmp/dlaxre_noadd.o

# Link: all objects except dlaxre_gk.o + original dlaxre
OBJS_NOADD=$(ls $OBJDIR/*.o | grep -v "dlaxre_gk\.o$" | tr '\n' ' ')
$FC -shared -o libxmr_noadd.so /tmp/dlaxre_noadd.o $OBJS_NOADD -llapack -lblas
echo "  Built libxmr_noadd.so ($(stat -c%s libxmr_noadd.so) bytes)"

echo "=== Done ==="
