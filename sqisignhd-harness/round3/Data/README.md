# Compact format vectors, level I

Self-generated known-answer vectors of the compact (dimension-4) format at
`p = 3·2^324 − 1`, `e = 174`, `r = 89`, `λ = 128`, in the SQIsignHD
library's text format (three lines per public key, eight per signature).
Provenance: `crates/sqisign-rs/tests/compact_vectors.rs` with the seeds
`"compact level I vector <i>"` / domain `kat`, messages `"message <i>"`,
five keys and one signature each; `COMPACT_VECTORS_OUT=<dir> cargo test
-p sqisign-rs --features compact --release --test compact_vectors`
rewrites them. The parameters are ours (docs/COMPACT_R3.md); no reference
implementation exists at this prime. `Signatures_r3lvl1_bad.txt` holds
one tampering per vector (a scalar, the degree, the challenge, in turn)
that every verifier must reject. The two `*_TABLE_r3lvl1.txt` files are
the non-residue tables of the library's hint convention for this prime.
