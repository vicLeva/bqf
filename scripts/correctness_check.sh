#!/usr/bin/env bash
# correctness_check.sh — build dev + main, run same queries, diff outputs.
#
# Usage:
#   ./scripts/correctness_check.sh               # fast: uses pre-counted ecoli data
#   ./scripts/correctness_check.sh --large        # counts fresh 31-mers on 5 ecoli genomes
#   ./scripts/correctness_check.sh --skip-build   # skip cmake/make
#   ./scripts/correctness_check.sh --skip-count   # reuse previously counted file (--large only)
set -euo pipefail

TUNA="$HOME/documents/giulio_colab/softs/tuna/build/tuna"
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$PROJECT/scripts/.correctness_tmp"
JOBS=$(nproc)

# ── Flags ──────────────────────────────────────────────────────────────────
SKIP_BUILD=0; SKIP_COUNT=0; LARGE=0
for arg in "$@"; do
    [[ "$arg" == "--skip-build" ]] && SKIP_BUILD=1
    [[ "$arg" == "--skip-count" ]] && SKIP_COUNT=1
    [[ "$arg" == "--large"      ]] && LARGE=1
done

mkdir -p "$TMP"

# ── 1. Build dev ────────────────────────────────────────────────────────────
DEV_BIN="$PROJECT/build/bin/bqf"
if [[ $SKIP_BUILD -eq 0 ]]; then
    echo "[1/4] Building dev..."
    cmake -S "$PROJECT" -B "$PROJECT/build" -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=OFF -Wno-dev > /dev/null 2>&1
    make -C "$PROJECT/build" -j"$JOBS" bqf 2>&1 | grep -E "error:|Built target bqf|Linking" || true
fi
[[ -x "$DEV_BIN" ]] || { echo "ERROR: dev binary not found at $DEV_BIN"; exit 1; }

# ── 2. Build main in a worktree ─────────────────────────────────────────────
MAIN_WT="$TMP/main_worktree"
MAIN_BIN="$MAIN_WT/build/bin/bqf"
if [[ $SKIP_BUILD -eq 0 ]]; then
    echo "[2/4] Building main..."
    if [[ ! -d "$MAIN_WT" ]]; then
        git -C "$PROJECT" worktree add "$MAIN_WT" main 2>/dev/null
    fi
    cmake -S "$MAIN_WT" -B "$MAIN_WT/build" -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=OFF -Wno-dev > /dev/null 2>&1
    make -C "$MAIN_WT/build" -j"$JOBS" bqf 2>&1 | grep -E "error:|Built target bqf|Linking" || true
fi
[[ -x "$MAIN_BIN" ]] || { echo "ERROR: main binary not found at $MAIN_BIN"; exit 1; }

# ── 3. Prepare counted k-mers ───────────────────────────────────────────────
if [[ $LARGE -eq 1 ]]; then
    # tuna is compiled for k=31; use z=1 so K=32, S=31
    K=32; Z=1; S=31; Q=20; C=5
    INPUT_LIST="$HOME/tmp/data/ecoli_fof_5.list"
    QUERIES="$PROJECT/examples/data/queries.fasta"
    COUNTED="$TMP/counted_31mers.tsv"
    if [[ $SKIP_COUNT -eq 0 || ! -f "$COUNTED" ]]; then
        echo "[3/4] Counting ${S}-mers with tuna (5 ecoli genomes)..."
        "$TUNA" -k "$S" -hp @"$INPUT_LIST" "$COUNTED"
        echo "      $(wc -l < "$COUNTED") distinct ${S}-mers counted"
    else
        echo "[3/4] Reusing existing count file ($COUNTED)"
    fi
else
    # Fast path: use pre-counted ecoli 28-mers from examples/
    K=32; Z=4; S=28; Q=18; C=5
    COUNTED="$PROJECT/examples/data/ecoli_28_counted.txt"
    QUERIES="$PROJECT/examples/data/queries.fasta"
    echo "[3/4] Using pre-counted ecoli 28-mers (examples/data/ecoli_28_counted.txt)"
fi

# ── 4. Build + query + diff ─────────────────────────────────────────────────
echo "[4/4] Build → query → diff..."

for branch in dev main; do
    BIN="$DEV_BIN"; [[ "$branch" == "main" ]] && BIN="$MAIN_BIN"
    "$BIN" build -q "$Q" -c "$C" -k "$K" -z "$Z" -i "$COUNTED" -o "$TMP/index_$branch" 2>/dev/null
    "$BIN" query -b "$TMP/index_$branch" -i "$QUERIES" -o "$TMP/out_$branch.txt" 2>/dev/null
done

echo ""
if diff -q "$TMP/out_dev.txt" "$TMP/out_main.txt" > /dev/null; then
    echo "PASS — dev and main outputs are identical"
    cat "$TMP/out_dev.txt"
else
    echo "INFO — outputs differ (showing diff):"
    diff "$TMP/out_dev.txt" "$TMP/out_main.txt" | head -40 || true
    echo ""
    # Check structural fields only (min, max, presence ratio) — ignore average
    # because dev intentionally fixed integer division in average computation.
    sed 's/average:[0-9.]*/average:X/g' "$TMP/out_dev.txt"  > "$TMP/out_dev_norm.txt"
    sed 's/average:[0-9.]*/average:X/g' "$TMP/out_main.txt" > "$TMP/out_main_norm.txt"
    if diff -q "$TMP/out_dev_norm.txt" "$TMP/out_main_norm.txt" > /dev/null; then
        echo "PASS — min/max/presence-ratio match; average differs only due to fixed integer division bug"
    else
        echo "FAIL — structural difference in min/max/presence-ratio:"
        diff "$TMP/out_dev_norm.txt" "$TMP/out_main_norm.txt"
        exit 1
    fi
fi
