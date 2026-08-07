#!/usr/bin/env bash
set -euo pipefail

bundle=/work/xmr-evolve/python_fortran/lapack_pr_dbdsvr
results=/work/results/dbdsvr_paper29_c2_1cpu_2026-08-07

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

mkdir -p "${results}/raw" "${results}/build_logs"

for n in 500 1000 2000 3000; do
  echo "BUILD n=${n}"
  make -C "${bundle}" clean >"${results}/build_logs/n${n}.log" 2>&1
  make -C "${bundle}" -j1 FC=gfortran-12 BENCH_MAXN="${n}" \
    paper29-benchmark >>"${results}/build_logs/n${n}.log" 2>&1

  for method in DBDSDC DBDSVR; do
    method_lower=${method,,}
    echo "RUN n=${n} method=${method}"
    taskset -c 0 "${bundle}/build/dbdsv_paper29_benchmark" \
      0 "${n}" A "${method}" 0.5 64 "" TIMING_ONLY \
      >"${results}/raw/n${n}_${method_lower}.csv"
    rows=$(wc -l <"${results}/raw/n${n}_${method_lower}.csv")
    echo "DONE n=${n} method=${method} rows=${rows}"
  done
done

gfortran-12 --version >"${results}/compiler.txt"
cat /sys/fs/cgroup/cpu.max >"${results}/cpu.max"
grep -m1 'model name' /proc/cpuinfo >"${results}/cpu_model.txt"
grep Cpus_allowed_list /proc/self/status >"${results}/cpus_allowed.txt"
sha256sum \
  "${bundle}/Makefile" \
  "${bundle}/TESTING/dbdsv_benchmark.f90" \
  "${bundle}/TESTING/EIG/dbsvrmg_paper29.f" \
  "${bundle}/SRC/dbdsvr.f" \
  "${bundle}/SRC/dbdsvdmr3.f" \
  >"${results}/source_sha256.txt"

echo "ALL DONE"
