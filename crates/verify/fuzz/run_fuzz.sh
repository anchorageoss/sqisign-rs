#!/usr/bin/env bash
set -euo pipefail

FUZZ_DIR="$(cd "$(dirname "$0")" && pwd)"
DURATION=60        # seconds per target
LEVELS="l1"        # default: l1 only
TARGETS="all"      # default: all target types

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

Run SQIsign verification fuzz targets.

Options:
  -d, --duration SECS   Seconds per target (default: 60)
  -l, --levels LEVELS   Comma-separated levels: l1,l3,l5 or "all" (default: l1)
  -t, --targets TARGETS Comma-separated target types or "all" (default: all)
                        Types: verify, signature, pubkey, compressed, expanded, any
  -h, --help            Show this help

Examples:
  # Fuzz all L1 targets for 60s each (~6 min total)
  $0

  # Fuzz all levels, all targets, 120s each (~36 min total)
  $0 -l all -d 120

  # Fuzz only verify targets across all levels, 300s each (~15 min total)
  $0 -l all -t verify -d 300

  # Fuzz L1 and L3 verify+signature targets, 180s each (~12 min total)
  $0 -l l1,l3 -t verify,signature -d 180
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--duration) DURATION="$2"; shift 2 ;;
        -l|--levels)   LEVELS="$2"; shift 2 ;;
        -t|--targets)  TARGETS="$2"; shift 2 ;;
        -h|--help)     usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

if [[ "$LEVELS" == "all" ]]; then
    LEVELS="l1,l3,l5"
fi

if [[ "$TARGETS" == "all" ]]; then
    TARGETS="verify,signature,pubkey,compressed,expanded,any"
fi

IFS=',' read -ra LEVEL_ARR <<< "$LEVELS"
IFS=',' read -ra TARGET_ARR <<< "$TARGETS"

TARGET_MAP_verify="fuzz_verify"
TARGET_MAP_signature="fuzz_signature_from_bytes"
TARGET_MAP_pubkey="fuzz_pubkey_from_bytes"
TARGET_MAP_compressed="fuzz_compressed_from_bytes"
TARGET_MAP_expanded="fuzz_expanded_from_bytes"
TARGET_MAP_any="fuzz_any_signature"

get_fuzz_name() {
    local type="$1" level="$2"
    local varname="TARGET_MAP_${type}"
    local base="${!varname}"
    if [[ "$level" == "l1" ]]; then
        echo "$base"
    else
        echo "${base}_${level}"
    fi
}

FUZZ_TARGETS=()
for level in "${LEVEL_ARR[@]}"; do
    for type in "${TARGET_ARR[@]}"; do
        FUZZ_TARGETS+=("$(get_fuzz_name "$type" "$level")")
    done
done

TOTAL=$((${#FUZZ_TARGETS[@]} * DURATION))
echo "Fuzzing ${#FUZZ_TARGETS[@]} targets for ${DURATION}s each (${TOTAL}s total)"
echo "Levels: ${LEVELS} | Targets: ${TARGETS}"
echo ""

cd "$FUZZ_DIR"

FAILED=()
for target in "${FUZZ_TARGETS[@]}"; do
    SEED_DIR="seeds/$target"
    SEED_FLAG=""
    if [[ -d "$SEED_DIR" ]] && [[ -n "$(ls -A "$SEED_DIR" 2>/dev/null)" ]]; then
        SEED_FLAG="$SEED_DIR"
    fi

    echo "=== $target (${DURATION}s) ==="
    if cargo +nightly fuzz run "$target" $SEED_FLAG -- \
        -max_total_time="$DURATION" \
        -print_final_stats=1 2>&1 | tail -8; then
        echo ""
    else
        echo "FAILED: $target"
        FAILED+=("$target")
        echo ""
    fi
done

echo "==============================="
if [[ ${#FAILED[@]} -eq 0 ]]; then
    echo "All ${#FUZZ_TARGETS[@]} targets completed successfully."
else
    echo "FAILURES: ${FAILED[*]}"
    exit 1
fi
