#!/usr/bin/env bash
# The paper's benchmark session (Section "Setup" of main_standalone.tex):
# one session, one build, one pinned core.
#
#   paper_session.sh <out-dir> [core]
#
# Records the machine, builds the two binaries from clean (release, fat
# LTO, no target-cpu flag), then runs passes: the Rust binary (`paper`),
# the C reference's own benchmark (the-sqisign, --iterations=300 at level
# I, 100 at III and V) and, once, the heap binary. Two passes; a third
# when a median of the first two differs by more than 3 percent
# (paper_report.py decides and reports which). Raw files and the exact
# commands go under <out-dir>.
set -u
OUT="${1:?out dir}"
CORE="${2:-11}"
REF="${SQISIGN_REF:-$HOME/the-sqisign}"
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUT"
LOG="$OUT/session.log"
say() { echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) $*" | tee -a "$LOG"; }
run() { say "+ $*"; "$@"; }

{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host kernel: $(uname -srm)"
  echo "rustc: $(rustc --version)"
  echo "cargo: $(cargo --version)"
  echo "RUSTFLAGS: '${RUSTFLAGS:-}' CARGO_BUILD_RUSTFLAGS: '${CARGO_BUILD_RUSTFLAGS:-}'"
  echo "machine type: $(curl -s -m 3 -H 'Metadata-Flavor: Google' http://metadata.google.internal/computeMetadata/v1/instance/machine-type 2>/dev/null | sed 's#.*/##')"
  echo "pinned core: $CORE (siblings: $(cat /sys/devices/system/cpu/cpu$CORE/topology/thread_siblings_list 2>/dev/null))"
  echo "governor: $(cat /sys/devices/system/cpu/cpu$CORE/cpufreq/scaling_governor 2>/dev/null || echo 'none exposed')"
  echo "load before: $(cat /proc/loadavg)"
  echo "--- lscpu"; lscpu
  echo "--- /proc/meminfo (head)"; head -5 /proc/meminfo
  echo "--- repository"
  git -C "$HERE" rev-parse HEAD
  git -C "$HERE" status --short | head -20
  echo "--- C reference"
  git -C "$REF" describe --tags 2>/dev/null; git -C "$REF" rev-parse HEAD
  grep -E 'CMAKE_BUILD_TYPE|SQISIGN_BUILD_TYPE|CMAKE_C_COMPILER:' "$REF/build/CMakeCache.txt"
  echo "--- busiest processes"; ps -eo pcpu,psr,comm --sort=-pcpu | head -8
} > "$OUT/machine.txt" 2>&1

say "build from clean: release, lto fat, codegen-units 1 (Cargo.toml), no target-cpu"
( cd "$HERE" && run cargo clean && run cargo build --release --bin paper --bin paper_heap ) >> "$LOG" 2>&1 || { say "build failed"; exit 1; }
sha256sum "$HERE/target/release/paper" "$HERE/target/release/paper_heap" "$REF"/build/apps/benchmark_p* | tee -a "$LOG" > "$OUT/binaries.sha256"

pass() {
  local p="$OUT/pass$1"
  mkdir -p "$p"
  say "pass $1: Rust"
  run taskset -c "$CORE" "$HERE/target/release/paper" --out "$p" 2>> "$LOG"
  say "pass $1: C reference"
  run taskset -c "$CORE" "$REF/build/apps/benchmark_p324_3" --iterations=300 > "$p/c-ref-I.txt" 2>&1
  run taskset -c "$CORE" "$REF/build/apps/benchmark_p500_27" --iterations=100 > "$p/c-ref-III.txt" 2>&1
  run taskset -c "$CORE" "$REF/build/apps/benchmark_p664_17" --iterations=100 > "$p/c-ref-V.txt" 2>&1
  echo "load after pass $1: $(cat /proc/loadavg)" >> "$OUT/machine.txt"
}

say "session start"
pass 1
pass 2
if ! python3 "$HERE/paper_report.py" "$OUT" --check-noise >> "$LOG" 2>&1; then
  say "a median moved by more than 3 percent between the two passes: third pass"
  pass 3
fi
say "heap"
run taskset -c "$CORE" "$HERE/target/release/paper_heap" --out "$OUT"
say "session end"
python3 "$HERE/paper_report.py" "$OUT" | tee -a "$LOG"
