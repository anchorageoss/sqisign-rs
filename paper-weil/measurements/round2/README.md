# Round-2 measurements behind paper-weil/main.tex

Every `% source:` comment in the paper that names `paper_pairing_audit.rs`,
`paper_stats.rs` or `paper_timing.rs` points here. The round-3 numbers come
from `crates/kat/tests/kat_compressed.rs` in this repository.

The round-2 code is sqisign-rs at tag `v0.4.31` (the last round-2 release),
which this tree no longer contains. The vectors are the official NIST
round-2 KATs of the-sqisign at tag `nist-v2` (commit `91e9e46`), 100 per
level. Setup used on 2026-10-22:

```sh
git worktree add /path/r2 v0.4.31
cd /path/r2 && rm -rf reference
git clone --branch nist-v2 --depth 1 https://github.com/SQISign/the-sqisign.git reference
cp tests/paper_pairing_audit.rs crates/verify/tests/
cp tests/paper_stats.rs tests/paper_timing.rs crates/kat/tests/
git apply tests/verify_visibility_for_timing.patch   # 13 pub(crate) -> pub, timing only
```

Commands (release only; the KAT tests at that tag need `kat-compat`):

```sh
# formula audit, checks A and B, Tate sweep (variants a, b, c), timings
cargo test -p sqisign-verify --release --features kat-compat --test paper_pairing_audit -- --nocapture --test-threads=1 task1_formula_audit
cargo test -p sqisign-verify --release --features kat-compat --test paper_pairing_audit -- --nocapture --test-threads=1 task2_tate_sweep
PAPER_TATE_M=188 PAPER_TATE_MC=187 cargo test -p sqisign-verify --release --features kat-compat --test paper_pairing_audit -- --nocapture --test-threads=1 timing_l1_kat0
# KAT statistics (a) to (g) and the cost breakdown; the tag has no release
# profile, so fat LTO and one codegen unit are supplied by the environment
CARGO_PROFILE_RELEASE_LTO=fat CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 cargo test -p sqisign-kat --release --features kat-compat --test paper_stats -- --nocapture --test-threads=1
CARGO_PROFILE_RELEASE_LTO=fat CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 cargo test -p sqisign-kat --release --features kat-compat --test paper_timing -- --nocapture --test-threads=1
```

`logs/` holds the printed output of each run as the paper quotes it
(`task2_tate_sweep_full.txt.gz` is the complete per-m sweep,
`task2_tate_sweep_summary.txt` its run-length view). Machine:
`logs/lscpu.txt` and `logs/environment.txt`.

## The fixed round-2 format (branch `fix/compressed-nbt-round2`)

`tests/fix_compressed_nbt_round2.patch` is the diff of that branch against
`v0.4.31` (the `[r:4 | h:1 | n_bt:3]` metadata byte, the one-bit hint, the
rejection of `n_bt >= 8` and `r >= 16`, and the tests);
`tests/kat_compressed_round2.rs` is its KAT test and
`logs/kat_compressed_fix.log` the run the paper quotes (300/300 round trips,
level III entry 97 with `n_bt = 4` included, the hint histograms and the
contingency counts):

```sh
git apply tests/fix_compressed_nbt_round2.patch        # on the v0.4.31 worktree
cargo test -p sqisign-kat --release --features kat-compat --test kat_compressed -- --nocapture --test-threads=1
```
