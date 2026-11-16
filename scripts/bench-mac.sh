#!/usr/bin/env bash
# The dimension-4 paper's benchmark on Apple Silicon (macOS, aarch64), one
# command:
#
#   scripts/bench-mac.sh [--with-c-reference]
#
# Builds the harness in release mode (fat LTO, no target-cpu flag), runs two
# full passes with the wall-clock timer and the portable field kernels, then
# the heap measurement, and writes results/<UTC>-m4/ with the raw logs,
# MACHINE.md and TABLES.md (milliseconds, the row layout of the paper's
# Table 3, the two-point over four-point ratio and the "times dimension 2"
# column). It prints the path of one tar.gz of that folder at the end.
#
# --with-c-reference also clones the-sqisign at nist-v3 into the results
# folder, builds its portable ("ref") variant with CMake and runs its own
# benchmark (--iterations=300 at level I, 100 at III and V). Without the
# flag the C rows read "not run".
#
# Needs: rustup with a stable toolchain, the Xcode command line tools, and
# python3 (shipped with the command line tools); with --with-c-reference
# also git and cmake. It installs nothing and needs no sudo. It also runs on
# Linux, for testing the script itself.
set -euo pipefail

WITH_C=0
for arg in "$@"; do
  case "$arg" in
    --with-c-reference) WITH_C=1 ;;
    -h | --help)
      sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *)
      echo "bench-mac.sh: unknown argument: $arg" >&2
      exit 2
      ;;
  esac
done

die() {
  echo "bench-mac.sh: $*" >&2
  exit 1
}

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BENCH="$ROOT/paper-dim4/bench"
OS="$(uname -s)"

command -v rustup > /dev/null 2>&1 || die "rustup not found. Install it from https://rustup.rs and run this again."
command -v cargo > /dev/null 2>&1 || die "cargo not found. Run: rustup toolchain install stable"
command -v rustc > /dev/null 2>&1 || die "rustc not found. Run: rustup toolchain install stable"
command -v python3 > /dev/null 2>&1 || die "python3 not found. On macOS it comes with the Xcode command line tools: xcode-select --install"
if [ "$OS" = "Darwin" ]; then
  xcode-select -p > /dev/null 2>&1 || die "the Xcode command line tools are missing (the linker). Run: xcode-select --install"
fi
command -v cc > /dev/null 2>&1 || die "no C compiler (cc) found; the Rust build needs a linker."
if [ "$WITH_C" = 1 ]; then
  command -v git > /dev/null 2>&1 || die "--with-c-reference needs git."
  command -v cmake > /dev/null 2>&1 || die "--with-c-reference needs cmake (for example: brew install cmake). Nothing was installed."
fi
if [ -n "${RUSTFLAGS:-}" ]; then
  echo "bench-mac.sh: note: RUSTFLAGS is set to '$RUSTFLAGS'; the paper's build uses none." >&2
fi

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
OUT="$ROOT/results/$STAMP-m4"
mkdir -p "$OUT"
LOG="$OUT/session.log"
say() { echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) $*" | tee -a "$LOG"; }

{
  echo "# MACHINE"
  echo
  echo "- date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  if [ "$OS" = "Darwin" ]; then
    echo "- cpu: $(sysctl -n machdep.cpu.brand_string)"
    echo "- hw.ncpu: $(sysctl -n hw.ncpu)"
    echo "- hw.memsize: $(sysctl -n hw.memsize)"
    echo "- sw_vers: $(sw_vers -productName) $(sw_vers -productVersion) ($(sw_vers -buildVersion))"
    echo "- power: $(pmset -g batt 2> /dev/null | head -1 || echo unknown)"
  else
    echo "- cpu: $(grep -m1 'model name' /proc/cpuinfo 2> /dev/null | cut -d: -f2- || echo unknown)"
    echo "- cpus: $(getconf _NPROCESSORS_ONLN)"
    echo "- memory: $(grep MemTotal /proc/meminfo 2> /dev/null || echo unknown)"
  fi
  echo "- uname: $(uname -srm)"
  echo "- rustc: $(rustc --version)"
  echo "- cargo: $(cargo --version)"
  echo "- git rev: $(git -C "$ROOT" rev-parse HEAD 2> /dev/null || echo unknown)"
  echo "- git status: $(git -C "$ROOT" status --short 2> /dev/null | wc -l | tr -d ' ') changed files"
  echo "- build: release, lto fat, codegen-units 1 (paper-dim4/bench/Cargo.toml), no target-cpu flag"
  echo "- C reference: $([ "$WITH_C" = 1 ] && echo "the-sqisign nist-v3, portable (ref) build, built without a cycle counter so that its benchmark reports nanoseconds" || echo "not run")"
} > "$OUT/MACHINE.md"

echo
echo "Run this with nothing else open and on mains power. Two passes follow; each takes several minutes."
echo

say "build: cargo build --release --bin paper --bin paper_heap"
(cd "$BENCH" && cargo build --release --bin paper --bin paper_heap) >> "$LOG" 2>&1 || die "the build failed; see $LOG"
PAPER="$BENCH/target/release/paper"
HEAP="$BENCH/target/release/paper_heap"

C_BUILD=""
if [ "$WITH_C" = 1 ]; then
  say "C reference: clone the-sqisign at nist-v3 and build its portable variant"
  git clone -q --depth 1 --branch nist-v3 https://github.com/SQISign/the-sqisign "$OUT/the-sqisign" >> "$LOG" 2>&1 || die "could not clone the-sqisign; see $LOG"
  C_BUILD="$OUT/the-sqisign/build"
  # -DNO_CYCLE_COUNTER: the reference's benchmark then reads the wall clock
  # in nanoseconds, so the figures it labels megacycles are milliseconds.
  cmake -S "$OUT/the-sqisign" -B "$C_BUILD" -DCMAKE_BUILD_TYPE=Release -DSQISIGN_BUILD_TYPE=ref -DCMAKE_C_FLAGS=-DNO_CYCLE_COUNTER >> "$LOG" 2>&1 || die "cmake failed for the C reference; see $LOG"
  cmake --build "$C_BUILD" -j 4 >> "$LOG" 2>&1 || die "the C reference did not build; see $LOG"
  export C_REF_UNIT=ms
fi

pass() {
  local p="$OUT/pass$1"
  mkdir -p "$p"
  say "pass $1: Rust"
  "$PAPER" --out "$p" 2>&1 | tee -a "$LOG"
  if [ -n "$C_BUILD" ]; then
    say "pass $1: C reference"
    "$C_BUILD/apps/benchmark_p324_3" --iterations=300 > "$p/c-ref-I.txt" 2>&1
    "$C_BUILD/apps/benchmark_p500_27" --iterations=100 > "$p/c-ref-III.txt" 2>&1
    "$C_BUILD/apps/benchmark_p664_17" --iterations=100 > "$p/c-ref-V.txt" 2>&1
  fi
}

say "session start"
pass 1
pass 2
say "heap"
"$HEAP" --out "$OUT" >> "$LOG" 2>&1
say "session end"

python3 "$BENCH/paper_report.py" "$OUT" > "$OUT/report.txt"
{
  echo "# TABLES: $(head -3 "$OUT/MACHINE.md" | tail -1 | sed 's/^- //')"
  echo
  echo "Field kernels: $(cat "$OUT/pass1/backend.txt"). Timer: $(head -1 "$OUT/pass1/timer.txt")."
  echo "Rows that moved by more than 3 percent between the two passes are listed in numbers.json."
  echo
  python3 "$BENCH/paper_tables_ms.py" "$OUT"
} > "$OUT/TABLES.md"

rm -rf "$OUT/the-sqisign"
TAR="$ROOT/results/$STAMP-m4.tar.gz"
tar -czf "$TAR" -C "$ROOT/results" "$STAMP-m4"
echo
cat "$OUT/TABLES.md"
echo
echo "Results: $OUT"
echo "Send back this file: $TAR"
