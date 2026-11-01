#!/usr/bin/env bash
# The machine record of a session (SECENG-2322 method): written before the
# session into <out dir>/machine.txt. Usage: machine.sh <out dir> [core]
set -u
OUT="${1:?out dir}"; CORE="${2:-3}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
SQ="${SQISIGN_V3_BW:-$HOME/the-sqisign/build/apps}"
mkdir -p "$OUT"
{
  echo "machine record, $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "== lscpu"; lscpu
  echo "== tsc flags: $(grep -o 'constant_tsc\|nonstop_tsc\|tsc_known_freq\|rdtscp' /proc/cpuinfo | sort -u | tr '\n' ' ')"
  echo "== dmesg tsc lines"; sudo -n dmesg 2>/dev/null | grep -i 'tsc' || echo "(dmesg not readable)"
  echo "== kernel: $(uname -a)"
  echo "== pinned core $CORE; the session runs every block as: taskset -c $CORE <binary>"
  echo "== cpufreq governor of core $CORE: $(cat /sys/devices/system/cpu/cpu$CORE/cpufreq/scaling_governor 2>/dev/null || echo 'no cpufreq interface (virtual machine)')"
  echo "== toolchain: $(rustc --version); $(cargo --version); $(gcc --version | head -1); $(cmake --version | head -1)"
  echo "== this crate: $(git -C "$ROOT" rev-parse HEAD) on $(git -C "$ROOT" branch --show-current) ($(git -C "$ROOT" status --short | wc -l) uncommitted paths)"
  echo "== session binary: $(ls -l --time-style=+%Y-%m-%dT%H:%M:%SZ "$HERE/target/release/dim4-session" | awk '{print $6, $7}'); profile: release, lto fat, codegen-units 1, no target-cpu flag"
  echo "== the-sqisign: $(git -C "$HOME/the-sqisign" rev-parse HEAD) ($(git -C "$HOME/the-sqisign" describe --tags 2>/dev/null)); build dir $SQ; $(grep -E '^CMAKE_BUILD_TYPE|^SQISIGN_BUILD_TYPE|^CMAKE_C_FLAGS_RELEASE' "$SQ/../CMakeCache.txt" | tr '\n' ' ')"
  echo "== the-sqisign binaries: $(ls -l --time-style=+%Y-%m-%dT%H:%M:%SZ "$SQ"/benchmark_p324_3 "$SQ"/benchmark_p500_27 "$SQ"/benchmark_p664_17 | awk '{print $6, $7}' | tr '\n' ';')"
  echo "== dpkg.log lines: $(wc -l < /var/log/dpkg.log); apt timers: $(systemctl is-active apt-daily.timer apt-daily-upgrade.timer | tr '\n' ' ')"
  echo "== memory: $(grep MemTotal /proc/meminfo)"
} > "$OUT/machine.txt"
echo "wrote $OUT/machine.txt"
