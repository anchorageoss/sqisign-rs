#!/usr/bin/env bash
# The measurement session of the dimension-4 two-point work (SECENG-2423,
# SECENG-2424), by the method of the ECLIPSE session of 2026-10-27
# (SECENG-2322): one script, one pinned core, one build of one commit,
# every binary built before the start and not rebuilt during the session,
# raw per-iteration cycles on disk, timestamps around every block, the
# package log and the process list checked before and after, the C
# reference's own benchmark in the same session on the same core.
#
#   paper-dim4/bench/session.sh <out dir> [core] [rounds]
#
# Blocks run round-robin: in each round every Rust block once, then one C
# reference level (I, III, V in turn, once each over the first three rounds);
# the fp2-ops block (the per-operation costs) runs once per round after the
# profile blocks.
# Run after `cargo build --release` in paper-dim4/bench; the script refuses
# to run if a cargo, LaTeX or Sage job is running. The legacy all-in-one
# run of 2026-09-26 is `dim4-session legacy` (feature round2).
set -u
OUT="${1:?out dir}"; CORE="${2:-3}"; ROUNDS="${3:-4}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
BIN="$HERE/target/release/dim4-session"
SQ="${SQISIGN_V3_BW:-$HOME/the-sqisign/build/apps}"
mkdir -p "$OUT/raw"
LOG="$OUT/session.log"
stamp() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { printf '%s %s\n' "$(stamp)" "$*" | tee -a "$LOG"; }
pin() { taskset -c "$CORE" "$@"; }
block() { # $1 name, $2... command
  local name=$1; shift
  log "block $name start"
  pin "$@" >> "$OUT/blocks.stdout" 2>&1
  log "block $name end rc=$?"
}

# Per-round iteration counts: the session totals are ROUNDS times these.
R2_KEYGEN=13; R2_SIGN=13; R2_VERIFY=50          # dimension 2, per level
R4_KEYGEN=13; R4_SIGN=25; R4_VERIFY=50; R4_PROFILE=25; SIGS=5
C_ITERS="p324_3:300 p500_27:100 p664_17:100"

log "session start; out $OUT; core $CORE; rounds $ROUNDS; per round: dim2 keygen/sign/verify $R2_KEYGEN/$R2_SIGN/$R2_VERIFY, dim4 keygen/sign/verify/profile $R4_KEYGEN/$R4_SIGN/$R4_VERIFY/$R4_PROFILE over $SIGS signatures, fp2-ops once; C iterations $C_ITERS"
log "commit $(git -C "$ROOT" rev-parse HEAD) ($(git -C "$ROOT" status --short | wc -l) uncommitted paths)"
log "binaries: $(ls -l --time-style=+%Y-%m-%dT%H:%M:%SZ "$BIN" | awk '{print $6, $7}') sha256 $(sha256sum "$BIN" | cut -c1-16); C reference: the-sqisign $(git -C "$HOME/the-sqisign" rev-parse --short HEAD) $SQ ($(ls -l --time-style=+%Y-%m-%dT%H:%M:%SZ "$SQ/benchmark_p324_3" | awk '{print $6}'))"
OTHER="$(pgrep -af 'cargo|rustc|pdflatex|bibtex|sage-venv|/usr/bin/sage' | grep -v "$$" | grep -v pgrep || true)"
if [ -n "$OTHER" ]; then log "REFUSED: other jobs running: $OTHER"; exit 2; fi
log "process check before: clean"
DPKG_START="$(stamp)"
log "dpkg.log lines before the window: $(wc -l < /var/log/dpkg.log)"
log "apt timers: $(systemctl is-active apt-daily.timer apt-daily-upgrade.timer | tr '\n' ' ')"
log "cpu $(lscpu | awk -F: '/Model name/{print $2}' | xargs); tsc $(grep -o 'constant_tsc\|nonstop_tsc\|tsc_known_freq' /proc/cpuinfo | sort -u | tr '\n' ' '); kernel $(uname -r); tsc from dmesg: $(sudo -n dmesg 2>/dev/null | grep -o 'Detected [0-9.]* MHz processor' | head -1)"

C_LEVELS=($C_ITERS)
for ((r = 1; r <= ROUNDS; r++)); do
  log "round $r start"
  for lvl in I III V; do
    block "dim2-$lvl r$r" "$BIN" "dim2-$lvl" --out "$OUT" --round "$r" --runs-keygen "$R2_KEYGEN" --runs-sign "$R2_SIGN" --runs-verify "$R2_VERIFY"
  done
  block "dim4-keygen r$r" "$BIN" dim4-keygen --out "$OUT" --round "$r" --runs-keygen "$R4_KEYGEN"
  block "dim4-sign r$r" "$BIN" dim4-sign --out "$OUT" --round "$r" --runs-sign "$R4_SIGN"
  for mode in four two-check two-nocheck; do
    block "dim4-verify-$mode r$r" "$BIN" "dim4-verify-$mode" --out "$OUT" --round "$r" --runs-verify "$R4_VERIFY" --sigs "$SIGS"
  done
  for mode in two-check four; do
    block "profile-$mode r$r" "$BIN" "profile-$mode" --out "$OUT" --round "$r" --runs "$R4_PROFILE" --sigs "$SIGS"
  done
  block "fp2-ops r$r" "$BIN" fp2-ops --out "$OUT" --round "$r"
  idx=$((r - 1))
  if [ "$idx" -lt "${#C_LEVELS[@]}" ]; then
    spec="${C_LEVELS[$idx]}"; prime="${spec%%:*}"; it="${spec##*:}"
    log "block sqisign-c $prime start ($it iterations)"
    pin "$SQ/benchmark_$prime" --iterations="$it" --seed=1 > "$OUT/sqisign-c-$prime.txt" 2>&1
    log "block sqisign-c $prime end rc=$?"
  fi
  log "round $r end"
done

log "process check after: $(pgrep -af 'cargo|rustc|pdflatex|bibtex|sage' | grep -v pgrep | grep -v "$$" || echo clean)"
log "dpkg.log lines after the window: $(wc -l < /var/log/dpkg.log); entries in the window:"
awk -v s="${DPKG_START:0:10} ${DPKG_START:11:8}" '$1" "$2 >= s' /var/log/dpkg.log | tee -a "$LOG" | head -5
log "session end"
