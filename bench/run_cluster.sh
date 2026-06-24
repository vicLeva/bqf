#!/usr/bin/env bash
#
# [BENCH] Build + run the resize-strategy comparison on a cluster node, sized to
# the available RAM, then render the plot.
#
# Usage:
#   bench/run_cluster.sh [RAM_GB] [REPEATS] [Q_MIN] [OUTDIR]
#
# Defaults: RAM_GB=256  REPEATS=3  Q_MIN=14  OUTDIR=bench_out
#
# The benchmark self-limits q_max: it stops before the estimated peak memory
# exceeds (RAM_GB - HEADROOM_GB). On 256 GB the largest point reached is q=31
# (~2.0e9 elements, ~157 GB peak). The whole sweep can take several hours,
# dominated by the largest q (fill + enumerate). Run it under nohup / a job
# scheduler. Results are flushed per-row, so partial CSVs are usable if killed.
#
# Run from the repository root.
set -euo pipefail

RAM_GB="${1:-256}"
REPEATS="${2:-3}"
Q_MIN="${3:-14}"
OUTDIR="${4:-bench_out}"

HEADROOM_GB=48                       # OS + allocator fragmentation safety margin
BUDGET=$(awk "BEGIN{print $RAM_GB - $HEADROOM_GB}")
Q_MAX=40                             # upper bound; the bench stops on the budget

BUILD_DIR=build_bench

echo "[run_cluster] RAM_GB=$RAM_GB  budget=${BUDGET}GB  repeats=$REPEATS  q_min=$Q_MIN  outdir=$OUTDIR"

# --- build (Release) ---
cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" --target resize_bench -j"$(nproc)"

mkdir -p "$OUTDIR"
CSV="$OUTDIR/resize_results.csv"
LOG="$OUTDIR/resize_results.log"
PLOT_PREFIX="$OUTDIR/resize"        # -> <prefix>_time.png, <prefix>_memory.png

# --- run sweep (self-limits q_max via the RAM budget) ---
echo "[run_cluster] running sweep -> $CSV (log: $LOG)"
"./$BUILD_DIR/bin/resize_bench" "$Q_MIN" "$Q_MAX" "$REPEATS" 42 "$BUDGET" \
    > "$CSV" 2> "$LOG" || { echo "[run_cluster] bench failed; see $LOG"; exit 1; }

cat "$LOG"

# --- plot (skip gracefully if matplotlib is unavailable on the node) ---
if python3 -c "import matplotlib" 2>/dev/null; then
    python3 bench/plot_resize.py "$CSV" "$PLOT_PREFIX"
    echo "[run_cluster] plots -> ${PLOT_PREFIX}_time.png, ${PLOT_PREFIX}_memory.png"
else
    echo "[run_cluster] matplotlib not found; copy $CSV locally and run:"
    echo "             python3 bench/plot_resize.py $CSV $PLOT_PREFIX"
fi

echo "[run_cluster] done."
