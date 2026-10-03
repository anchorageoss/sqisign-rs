# The paper's session harness

`session.sh` runs the one session `paper-dim4/main.tex` reads: the Rust
harness (this crate for round 3, the `v0.4.31` crate for round 2, one
process, `rdtsc` and wall-clock medians) and the C reference's own
benchmark apps for both rounds, and prints Markdown in the format of
paper-dim4/BENCH.md's session section.

```sh
SQISIGN_V3_BW=/path/to/the-sqisign-v3/build-bw \
SQISIGN_V2_BW=/path/to/the-sqisign-v2/build-bw ./session.sh
```

Each C tree is built with `cmake -DSQISIGN_BUILD_TYPE=broadwell
-DCMAKE_BUILD_TYPE=Release` (`nist-v2` needs GMP headers). Without the
environment variables only the Rust tables are printed.

## The M4 Pro column of Table 4

On the Mac, `cargo run --release` in this directory prints the same
tables; `rdtsc` is not available there, so the cycle columns read 0 and
the `median (ms)` column is the one to keep. Copy those thirteen
`median (ms)` values into the "M4 Pro" table of paper-dim4/BENCH.md's session
section (same rows, in place of `pending`) and run `make figures` in
`paper-dim4`: Table 4 and `\MFourAvailable` pick them up.
