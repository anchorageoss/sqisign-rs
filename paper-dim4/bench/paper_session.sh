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
# when a median of the first two differs by more than PAPER_NOISE (default
# 0.03; paper_report.py decides and reports which). Raw files and the exact
# commands go under <out-dir>.
#
# With QUIET=1 the script waits before the build is used and before every
# pass until the 1-minute load average has stayed under QUIET_LOAD (default
# 0.10) for five minutes, sampled every ten seconds, and stops with exit 3
# when that does not happen within thirty minutes. OTHER_CPUS (a taskset
# list) is where the build and the script's own helpers run, away from the
# measured core and its sibling. Steal time is logged by vmstat throughout.
set -u
OUT="${1:?out dir}"
CORE="${2:-11}"
REF="${SQISIGN_REF:-$HOME/the-sqisign}"
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUT"
export PAPER_NOISE="${PAPER_NOISE:-0.03}"
OTHER_CPUS="${OTHER_CPUS:-}"
other() { if [ -n "$OTHER_CPUS" ]; then taskset -c "$OTHER_CPUS" "$@"; else "$@"; fi; }
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
( cd "$HERE" && run other cargo clean && run other cargo build --release --bin paper --bin paper_heap ) >> "$LOG" 2>&1 || { say "build failed"; exit 1; }
sha256sum "$HERE/target/release/paper" "$HERE/target/release/paper_heap" "$REF"/build/apps/benchmark_p* | tee -a "$LOG" > "$OUT/binaries.sha256"

# Wait until the 1-minute load average has been under the bound for five
# minutes straight; exit 3 after thirty minutes without that.
wait_quiet() {
  [ "${QUIET:-0}" = 1 ] || return 0
  local bound="${QUIET_LOAD:-0.10}" ok=0 waited=0
  say "waiting for quiet: 1-minute load under $bound for 300 s"
  while [ $ok -lt 30 ]; do
    if awk -v b="$bound" '{exit !($1 < b)}' /proc/loadavg; then ok=$((ok + 1)); else ok=0; fi
    waited=$((waited + 10))
    if [ $waited -gt 1800 ] && [ $ok -lt 30 ]; then
      say "not quiet after 30 minutes: load $(cat /proc/loadavg)"
      { echo "--- not quiet at $(date -u +%FT%TZ): load $(cat /proc/loadavg)"; ps -eo pcpu,psr,pid,etime,args --sort=-pcpu | head -25; } >> "$OUT/machine.txt"
      exit 3
    fi
    sleep 10
  done
  say "quiet: load $(cat /proc/loadavg) after ${waited} s"
}

pass() {
  local p="$OUT/pass$1"
  mkdir -p "$p"
  wait_quiet
  echo "load before pass $1: $(cat /proc/loadavg)" >> "$OUT/machine.txt"
  say "pass $1: Rust"
  run taskset -c "$CORE" "$HERE/target/release/paper" --out "$p" 2>> "$LOG"
  say "pass $1: C reference"
  run taskset -c "$CORE" "$REF/build/apps/benchmark_p324_3" --iterations=300 > "$p/c-ref-I.txt" 2>&1
  run taskset -c "$CORE" "$REF/build/apps/benchmark_p500_27" --iterations=100 > "$p/c-ref-III.txt" 2>&1
  run taskset -c "$CORE" "$REF/build/apps/benchmark_p664_17" --iterations=100 > "$p/c-ref-V.txt" 2>&1
  echo "load after pass $1: $(cat /proc/loadavg)" >> "$OUT/machine.txt"
}

other vmstat -t 5 > "$OUT/vmstat.log" 2>&1 &
VMSTAT=$!
trap 'kill $VMSTAT 2>/dev/null' EXIT
say "session start"
pass 1
pass 2
if ! python3 "$HERE/paper_report.py" "$OUT" --check-noise >> "$LOG" 2>&1; then
  say "a median moved by more than $PAPER_NOISE between the two passes: third pass"
  pass 3
fi
wait_quiet
say "heap"
run taskset -c "$CORE" "$HERE/target/release/paper_heap" --out "$OUT"
say "session end"
python3 "$HERE/paper_report.py" "$OUT" | tee -a "$LOG"
