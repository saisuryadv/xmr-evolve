CXX = g++
CXXFLAGS = -std=c++17 -O2 -Isrc/clapack

# Pure C library — no gfortran, no Accelerate
LDFLAGS_C = -Llib -lxmr_c -lm

# Legacy Fortran library (kept for comparison)
LDFLAGS_FORTRAN = -Llib -lxmr -framework Accelerate -L/opt/homebrew/lib/gcc/15 -lgfortran -lm

# Default to pure C
LDFLAGS = $(LDFLAGS_C)

STCOLL = /Users/saisurya/MRRR/bidiag-algo/MRRR\ Resources/STCollection/DATA

.PHONY: all clean test lib lib_c lib_fortran evaluate evaluate_fortran

all: evaluate

lib: lib_c

lib_c:
	bash tools/build_c_lib.sh

lib_fortran:
	bash tools/build_xmr_lib.sh

evaluate: src/evaluate.cpp src/bidiag_svd.h lib/libxmr_c.a
	$(CXX) $(CXXFLAGS) -o $@ src/evaluate.cpp $(LDFLAGS_C)

evaluate_fortran: src/evaluate.cpp src/bidiag_svd.h
	$(CXX) $(CXXFLAGS) -o $@ src/evaluate.cpp $(LDFLAGS_FORTRAN)

test: evaluate
	./evaluate "$(STCOLL)" 100

clean:
	rm -f evaluate evaluate_fortran test_fortran_lib
