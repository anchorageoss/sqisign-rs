#!/bin/bash
# Runs the paper_timing test under each build configuration, sequentially,
# waiting for the host to be quiet (1-minute load average below 1.5, at most
# 5 minutes) before each run, and recording load before and after.
S=/tmp/claude-1000/-home-user-sqisign-rs/f14775d5-0441-4505-aadf-f508f743ea2b/scratchpad
cd "$S/r2b" || exit 1

wait_quiet() {
  for _ in $(seq 1 30); do
    l=$(cut -d' ' -f1 /proc/loadavg)
    if awk "BEGIN{exit !($l < 1.5)}"; then return 0; fi
    sleep 10
  done
  echo "WARNING: proceeding although load is $(cat /proc/loadavg)"
}

run_cfg() {
  name=$1; log=$2; shift 2
  wait_quiet
  {
    echo "config: $name"
    echo "start $(date -u +%FT%TZ) load: $(cat /proc/loadavg)"
    env "$@" cargo test -p sqisign-kat --release --features kat-compat --test paper_timing -- --nocapture --test-threads=1 2>&1
    echo "exit=${PIPESTATUS[0]}"
    echo "end $(date -u +%FT%TZ) load: $(cat /proc/loadavg)"
  } > "$log" 2>&1
}

run_cfg "primary (run 2): fat LTO + 1 CGU via env, .cargo/config.toml target-cpu=native in effect" \
  "$S/run_timing_primary2.log" CARGO_PROFILE_RELEASE_LTO=fat CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1
run_cfg "tag default profile (run 2): no profile overrides, .cargo/config.toml target-cpu=native in effect" \
  "$S/run_timing_tagdefault2.log" CARGO_TARGET_DIR="$S/target_tagdefault"
run_cfg "portable: fat LTO + 1 CGU via env, RUSTFLAGS empty (no target-cpu=native)" \
  "$S/run_timing_portable.log" RUSTFLAGS= CARGO_PROFILE_RELEASE_LTO=fat CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 CARGO_TARGET_DIR="$S/target_portable"
run_cfg "primary (run 3): fat LTO + 1 CGU via env, .cargo/config.toml target-cpu=native in effect" \
  "$S/run_timing_primary3.log" CARGO_PROFILE_RELEASE_LTO=fat CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1
echo "all done $(date -u +%FT%TZ)"
