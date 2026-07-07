# Standalone Makefile for the offline DBDSVR PR bundle.
# Requires a prebuilt LAPACK + BLAS pair.  Defaults point at the
# /tmp/lapack-ref build produced by cmake.

FC       ?= gfortran
FFLAGS   ?= -O2 -g -std=legacy -ffree-line-length-none
LAPACK_LIB ?= /tmp/lapack-ref/build/lib/liblapack.a
BLAS_LIB   ?= /tmp/lapack-ref/build/lib/libblas.a
TMG_LIB    ?= /tmp/lapack-ref/build/lib/libtmglib.a
LAPACK_ROOT ?= /tmp/lapack-ref
LIBS       ?= $(TMG_LIB) $(LAPACK_LIB) $(BLAS_LIB)

# LAPACK test residual/orthogonality helpers not shipped in the static
# libraries: pull them straight from the reference sources.
EIG_HELPERS := \
    $(LAPACK_ROOT)/TESTING/EIG/dbdt04.f \
    $(LAPACK_ROOT)/TESTING/EIG/dort01.f

BUILD_DIR := build
SRC_DIR   := SRC
EIG_DIR   := TESTING/EIG
TEST_DIR  := TESTING

SRC_FILES := \
    $(SRC_DIR)/dbdsvr.f \
    $(SRC_DIR)/dbdsvdmr3.f \
    $(SRC_DIR)/dbdtgk.f \
    $(SRC_DIR)/dlar1v_tgk.f \
    $(SRC_DIR)/dlarrf_tgk.f \
    $(SRC_DIR)/dlarrv_tgk.f \
    $(SRC_DIR)/stegr_ID/dlarrb.f

TEST_FILES := \
    $(EIG_DIR)/dbsvrmg.f \
    $(EIG_DIR)/dchkbsvr.f \
    $(EIG_DIR)/derrbsvr.f \
    $(TEST_DIR)/dbsvr_test.f

EIG_HELPER_OBJS := $(patsubst $(LAPACK_ROOT)/%.f,$(BUILD_DIR)/lref/%.o,$(EIG_HELPERS))
ALL_OBJS := $(patsubst %.f,$(BUILD_DIR)/%.o,$(SRC_FILES) $(TEST_FILES)) \
            $(EIG_HELPER_OBJS)

$(BUILD_DIR)/lref/%.o: $(LAPACK_ROOT)/%.f
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY: all clean run

all: $(BUILD_DIR)/dbsvr_test

$(BUILD_DIR)/dbsvr_test: $(ALL_OBJS)
	@mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -o $@ $(ALL_OBJS) $(LIBS)

$(BUILD_DIR)/%.o: %.f
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c $< -o $@

run: $(BUILD_DIR)/dbsvr_test
	./$(BUILD_DIR)/dbsvr_test < $(TEST_DIR)/svrtest.in

clean:
	rm -rf $(BUILD_DIR)
