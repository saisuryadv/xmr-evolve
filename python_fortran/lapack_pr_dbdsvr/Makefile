# Standalone Makefile for the offline DBDSVR PR bundle.
# Requires a prebuilt LAPACK + BLAS pair.  Defaults point at the
# /tmp/lapack-ref-3.11.0 build produced by build.sh.

ifeq ($(origin FC), default)
FC       = gfortran
endif
FFLAGS   ?= -O2 -g -std=legacy -ffree-line-length-none
FREEFLAGS ?= -O2 -g -std=f2008
BENCH_MAXN ?= 400
LAPACK_LIB ?= /tmp/lapack-ref-3.11.0/build/lib/liblapack.a
BLAS_LIB   ?= /tmp/lapack-ref-3.11.0/build/lib/libblas.a
TMG_LIB    ?= /tmp/lapack-ref-3.11.0/build/lib/libtmglib.a
LAPACK_ROOT ?= /tmp/lapack-ref-3.11.0
LIBS       ?= $(TMG_LIB) $(LAPACK_LIB) $(BLAS_LIB)

# Keep DBDSVR itself on the configured reference BLAS while allowing the
# post-solve email-metric GEMMs to use CBLAS.  On macOS Accelerate supplies
# cblas_dgemm without replacing the already-linked Fortran BLAS symbols.
ifeq ($(shell uname -s),Darwin)
CASE_METRIC_LIBS ?= $(LIBS) -framework Accelerate
else
CASE_METRIC_LIBS ?= $(LIBS) -lopenblas
endif

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

.PHONY: all benchmark paper29-benchmark clean fallback-test run probe

all: $(BUILD_DIR)/dbsvr_test $(BUILD_DIR)/dbdsv_solver_probe \
	$(BUILD_DIR)/dbdsv_benchmark $(BUILD_DIR)/dbdsvr_case_metrics \
	$(BUILD_DIR)/dbdsvr_fallback_test

$(BUILD_DIR)/dbsvr_test: $(ALL_OBJS)
	@mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -o $@ $(ALL_OBJS) $(LIBS)

$(BUILD_DIR)/%.o: %.f
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/TESTING/dbdsv_solver_probe.o: TESTING/dbdsv_solver_probe.f90
	@mkdir -p $(dir $@)
	$(FC) $(FREEFLAGS) -c $< -o $@

$(BUILD_DIR)/TESTING/dbdsv_benchmark.o: TESTING/dbdsv_benchmark.f90
	@mkdir -p $(dir $@)
	$(FC) $(FREEFLAGS) -cpp -DBENCH_MAXN=$(BENCH_MAXN) -c $< -o $@

$(BUILD_DIR)/TESTING/dbdsv_paper29_benchmark.o: TESTING/dbdsv_benchmark.f90
	@mkdir -p $(dir $@)
	$(FC) $(FREEFLAGS) -cpp -DBENCH_MAXN=$(BENCH_MAXN) \
		-DBENCH_NTYPES=29 -DPAPER29 -c $< -o $@

$(BUILD_DIR)/TESTING/dbdsvr_case_metrics.o: TESTING/dbdsvr_case_metrics.f90
	@mkdir -p $(dir $@)
	$(FC) $(FREEFLAGS) -c $< -o $@

$(BUILD_DIR)/TESTING/dbdsvr_fallback_test.o: TESTING/dbdsvr_fallback_test.f90
	@mkdir -p $(dir $@)
	$(FC) $(FREEFLAGS) -c $< -o $@

$(BUILD_DIR)/dbdsv_solver_probe: \
		$(BUILD_DIR)/TESTING/dbdsv_solver_probe.o \
		$(patsubst %.f,$(BUILD_DIR)/%.o,$(filter $(SRC_DIR)/%,$(SRC_FILES)))
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

probe: $(BUILD_DIR)/dbdsv_solver_probe

$(BUILD_DIR)/dbdsv_benchmark: \
		$(BUILD_DIR)/TESTING/dbdsv_benchmark.o \
		$(BUILD_DIR)/TESTING/EIG/dbsvrmg.o \
		$(EIG_HELPER_OBJS) \
		$(patsubst %.f,$(BUILD_DIR)/%.o,$(filter $(SRC_DIR)/%,$(SRC_FILES)))
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

benchmark: $(BUILD_DIR)/dbdsv_benchmark

$(BUILD_DIR)/dbdsv_paper29_benchmark: \
		$(BUILD_DIR)/TESTING/dbdsv_paper29_benchmark.o \
		$(BUILD_DIR)/TESTING/EIG/dbsvrmg_paper29.o \
		$(EIG_HELPER_OBJS) \
		$(patsubst %.f,$(BUILD_DIR)/%.o,$(filter $(SRC_DIR)/%,$(SRC_FILES)))
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

paper29-benchmark: $(BUILD_DIR)/dbdsv_paper29_benchmark

$(BUILD_DIR)/dbdsvr_case_metrics: \
		$(BUILD_DIR)/TESTING/dbdsvr_case_metrics.o \
		$(patsubst %.f,$(BUILD_DIR)/%.o,$(filter $(SRC_DIR)/%,$(SRC_FILES)))
	$(FC) $(FFLAGS) -o $@ $^ $(CASE_METRIC_LIBS)

$(BUILD_DIR)/dbdsvr_fallback_test: \
		$(BUILD_DIR)/TESTING/dbdsvr_fallback_test.o \
		$(BUILD_DIR)/SRC/dbdsvr.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

fallback-test: $(BUILD_DIR)/dbdsvr_fallback_test
	./$(BUILD_DIR)/dbdsvr_fallback_test

run: $(BUILD_DIR)/dbsvr_test
	./$(BUILD_DIR)/dbsvr_test < $(TEST_DIR)/svrtest-dbdsvr-only.in

clean:
	rm -rf $(BUILD_DIR)
