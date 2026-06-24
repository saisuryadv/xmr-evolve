#!/bin/bash
# Launch all 12 paper-norm sweeps (3 solvers x 4 suites) in the background.
# Advisor + dbdsqr sweeps run concurrently on dedicated cores (taskset) so
# they don't contend with each other; selfcontained must be serialized
# because mr3gk_run uses fixed paths /tmp/_mr3gk_in.bin and _mr3gk_out.bin.

set -u
cd "$(dirname "$0")"
OUT=docs/sweeps_2026-06-18
mkdir -p "$OUT"

# Suite list (kept in array form so we can loop without subshells).
SUITES=(pract synth 379 dense_to_bidiag)

# advisor: cores 0..3
# dbdsqr:  cores 4..7
# selfcontained serial driver: core 8

launch_one () {
    local solver="$1" suite="$2" core="$3"
    local log="$OUT/${suite}_${solver}_paper_norms.log"
    taskset -c "$core" python3 -u eval_dbdsvdmr3.py \
        --suite "$suite" --solver "$solver" --paper-norms \
        > "$log" 2>&1 &
    echo "[launch] solver=$solver suite=$suite core=$core pid=$! log=$log"
}

for i in 0 1 2 3; do
    launch_one advisor "${SUITES[$i]}" "$i"
done
for i in 0 1 2 3; do
    launch_one dbdsqr "${SUITES[$i]}" "$((4 + i))"
done

# Selfcontained serial driver (one suite at a time, pinned to core 8).
(
    for s in "${SUITES[@]}"; do
        log="$OUT/${s}_selfcontained_paper_norms.log"
        echo "[serial] selfcontained suite=$s start $(date -Iseconds)" >> "$OUT/selfcontained_serial.driver.log"
        taskset -c 8 python3 -u eval_dbdsvdmr3.py \
            --suite "$s" --solver selfcontained --paper-norms \
            > "$log" 2>&1
        echo "[serial] selfcontained suite=$s end   $(date -Iseconds)" >> "$OUT/selfcontained_serial.driver.log"
    done
) &
echo "[launch] selfcontained serial driver pid=$!"

wait
echo "[done] all 12 sweeps complete $(date -Iseconds)"
