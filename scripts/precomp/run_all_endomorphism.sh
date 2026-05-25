#!/usr/bin/env bash
# Run endomorphism precomputation for all 3 security levels.
# Designed to run overnight on a GCP VM with SageMath installed.
#
# Prerequisites:
#   - SageMath >= 10.0 (ideally 10.5+)
#   - 16+ GB RAM (for PARI allocation)
#   - The sqisign-rs repo checked out at the same state
#   - If using micromamba: activate the sage env first
#
# Usage:
#   cd /path/to/sqisign-rs
#   bash scripts/precomp/run_all_endomorphism.sh
#
# Output:
#   crates/precomp/src/level{1,3,5}/endomorphism_action.rs
#
# Each level takes 1-3+ hours depending on CPU speed.
# The script logs progress and saves results incrementally.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

LOG_DIR="$REPO_ROOT/scripts/precomp/logs"
mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)

for level in 1 3 5; do
    OUT_DIR="$REPO_ROOT/crates/precomp/src/level${level}"
    LOG_FILE="$LOG_DIR/endomorphism_level${level}_${TIMESTAMP}.log"

    echo "=== Level $level === ($(date))"
    echo "  Output: $OUT_DIR/endomorphism_action.rs"
    echo "  Log:    $LOG_FILE"

    # Run with unbuffered Python output so log shows real-time progress
    PYTHONUNBUFFERED=1 sage "$SCRIPT_DIR/run_endomorphism_precomp.sage" "$level" "$OUT_DIR" 2>&1 | tee "$LOG_FILE"

    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "[OK] Level $level completed successfully at $(date)"
    else
        echo "[FAIL] Level $level failed at $(date). See $LOG_FILE"
        echo "  Continuing to next level..."
    fi
    echo
done

echo "=== All levels done at $(date) ==="
echo "Check generated files:"
for level in 1 3 5; do
    f="$REPO_ROOT/crates/precomp/src/level${level}/endomorphism_action.rs"
    if [ -f "$f" ]; then
        echo "  [OK]   $f ($(wc -l < "$f") lines)"
    else
        echo "  [MISS] $f"
    fi
done
