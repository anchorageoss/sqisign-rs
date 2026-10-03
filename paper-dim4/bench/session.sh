#!/usr/bin/env bash
# The paper's one session: the Rust harness (both crate versions) and the C
# reference's own benchmark apps, round 3 (`nist-v3`) and round 2
# (`nist-v2`), both built `broadwell`, on this machine. Prints Markdown for
# paper-dim4/BENCH.md's session section. Usage:
#   SQISIGN_V3_BW=/path/the-sqisign-v3/build-bw SQISIGN_V2_BW=/path/the-sqisign-v2/build-bw ./session.sh
set -euo pipefail
cd "$(dirname "$0")"
echo "Machine: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //'); $(rustc --version); $(uname -sr)"
echo
cargo build --release --quiet
./target/release/dim4-session
c_rows() { # $1 = label, $2 = round, $3 = build dir, $4... = "prime level app"
  local label=$1 round=$2 dir=$3; shift 3
  for spec in "$@"; do
    set -- $spec; local level=$1 app=$2 iters=$3
    "$dir/apps/$app" --iterations="$iters" 2>/dev/null | awk -v L="$label" -v R="$round" -v LV="$level" -v IT="$iters" '
      /^ +(keypair|sign|verify) +\|/ {
        op=$1; if (op=="keypair") op="keygen"; if (op=="verify") op="verify from bytes";
        for (i=1;i<=NF;i++){ if($i=="average") mean=$(i+1); if($i=="median") med=$(i+1); if($i=="min") mn=$(i+1); if($i=="max") mx=$(i+1) }
        printf "| %s | %s | 2 | %s | %s | %s | %s | %s | %s | %s | |\n", L, R, LV, op, IT, med, mean, mn, mx }'
  done
}
echo
echo "#### The C reference (\`broadwell\` builds; its own benchmark, medians over iterations)"
echo
echo "| implementation | round | dimension | level | operation | runs | median (Mcycles) | mean | min | max | median (ms) |"
echo "|---|---|---|---|---|---|---|---|---|---|---|"
[ -n "${SQISIGN_V3_BW:-}" ] && c_rows "SQIsign (C reference, nist-v3, broadwell)" 3 "$SQISIGN_V3_BW" "I benchmark_p324_3 300" "III benchmark_p500_27 100" "V benchmark_p664_17 100"
[ -n "${SQISIGN_V2_BW:-}" ] && c_rows "SQIsign (C reference, nist-v2, broadwell)" 2 "$SQISIGN_V2_BW" "I benchmark_lvl1 300" "III benchmark_lvl3 100" "V benchmark_lvl5 100"
true
