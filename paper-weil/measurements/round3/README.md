# Round-3 measurements behind paper-weil/main.tex

Every `% source:` comment in the paper that names `paper_pairing_audit.rs`,
`paper_timing.rs` or `kat_compressed.rs` points here. The code is sqisign-rs
0.7 on branch `fix/compressed-nbt` of this repository; the vectors are the
official NIST round-3 KATs of the-sqisign at tag `nist-v3` (commit
`6d01770`), 100 per level, in `crates/kat/kat`. Runs of 2026-10-26.

The two harnesses live here and are copied into `crates/kat/tests` to run
(the `sqisign-kat` crate has the vectors, `rayon` and `hex`); they use only
the public API of `sqisign-verify` and copy the two private routines they
time (`check_canonical_basis_change_matrix`, `compressed::reconstruct`):

```sh
cp tests/paper_pairing_audit.rs tests/paper_timing.rs crates/kat/tests/
# theorem checks (Section "What was checked"), the Tate sweep and the
# pairing-stage timing (Section "Pairing selection")
cargo test -p sqisign-kat --release --test paper_pairing_audit -- --nocapture --test-threads=1 task1
cargo test -p sqisign-kat --release --test paper_pairing_audit -- --nocapture --test-threads=1 task2
cargo test -p sqisign-kat --release --test paper_pairing_audit -- --nocapture --test-threads=1 timing_l1_kat0
# the decompressor on one Tate ladder per curve at the bound and one below, end to end
cargo test -p sqisign-kat --release --test paper_pairing_audit -- --nocapture --test-threads=1 task3
# the cost table (Section "Costs"), level I over the 100 vectors: run three times
cargo test -p sqisign-kat --release --test paper_timing -- --nocapture --test-threads=1
# the compressed format on the vectors (the sizes and the round trips)
cargo test -p sqisign-kat --release --test kat_compressed -- --nocapture --test-threads=1
```

The release profile of `Cargo.toml` (fat LTO, one codegen unit, no
`target-cpu` flag; the x86-64 assembly field kernels are selected at run
time) is the one measured.

`logs/` holds the printed output of each run as the paper quotes it:
`task1_formula_audit.txt`, `task2_tate_sweep_full.txt.gz` (the complete
per-`m` sweep) with `task2_tate_sweep_summary.txt` (its threshold lines),
`task3_end_to_end_tate.txt` (the end-to-end Tate check, not quoted in
the paper), `timing_l1_kat0.txt` (the pairing stage on entry 0, Weil
against Tate),
`run_timing_1.log` to `run_timing_3.log` (the cost table, medians over
the 100 level-I vectors; the paper quotes run 1), `kat_compressed.log`.
Machine: `logs/lscpu.txt` and `logs/environment.txt`.
