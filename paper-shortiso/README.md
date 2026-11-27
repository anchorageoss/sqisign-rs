# paper-shortiso

Rust checks of the results of "Shortest isogenies between supersingular
elliptic curves with complex multiplication" (Barreto, Ray). The paper is
not in this repository. Findings: `REPORT.md`.

## Run everything

    sh paper-shortiso/run.sh

It builds five test binaries of `sqisign-rs` in release mode, runs them,
writes one log per item to `results/` and exits non-zero if any assertion
fails. Measured on a 32-core x86-64 machine on 2026-11-27: 56 s for the
cold release build, then 15 to 50 s for the run (about 2 minutes of CPU
time). A debug build passes too and takes about 2.5 minutes of wall time.

The same tests run under plain `cargo test -p sqisign-rs --release`. They
compare what they compute with the records checked in under `classpoly/`
and `endpoints/`; `SHORTISO_WRITE=1 sh paper-shortiso/run.sh` regenerates
those records instead.

## Where things are

| Item | Code (`crates/sqisign-rs/tests/`) | Log |
|---|---|---|
| 1, 2: the three SQIsign primes on the crate's field, curve, basis, stored action, quaternion layer, ideal-to-isogeny | `shortiso_sqisign.rs` | `results/item1.log`, `results/item2.log` |
| 3: other primes over a generic field | `shortiso_generic.rs`, `shortiso_common/` | `results/item3.log` |
| 4: the norm equations, integers only | `shortiso_equations.rs` | `results/item4.log` |
| 5: what stays with Sage | none | `results/item5.log` |
| 6: the (D, m, h(-D)) table, class polynomials | `shortiso_table.rs` | `results/item6.log`, `SOURCES.md` |
| 7: minimality by exact enumeration in the ideal classes | `shortiso_minimum.rs`, `shortiso_common/quat.rs` | `results/item7.log` |

`shortiso_common/` is test-only: `field.rs` (F_p and F_{p^2} over any
prime p = 3 mod 4, on num-bigint), `curve.rs` (Jacobian points, iota and
Frobenius on points, Velu isogenies of degree 2 and of odd prime degree,
discrete logarithms in smooth cyclic groups), `classpoly.rs` (reduced
forms, class numbers, class polynomials by fixed-point q-expansions),
`quat.rs` (the quaternion algebra in exact integers: lattices in Hermite
normal form, the norm form of a left ideal of O0, the exact minimum of a
quaternary form), `mod.rs` (the kernel of Corollary A, the chain walk, the
records).

## Records

- `results/primes.tsv`: every prime checked, with D, m, the number of
  steps, the endpoint (a j-invariant or the class polynomial it is a root
  of) and pass or fail.
- `classpoly/H_-D.txt`: the class polynomials used, D = 68, 164, 308, 468,
  1268, 2372, 2532, computed by `shortiso_common/classpoly.rs` and
  recomputed at every run.
- `endpoints/level{1,3,5}.txt`: a basis of the right order of
  I = O0 theta + O0 m at each SQIsign level, over (1, i, j, k). It is the
  endomorphism ring of the endpoint of the chain.
- `SOURCES.md`: file and section references for the claims about SQIsign
  and MIKE.
- `notes/iso_isogeny_2chain_from-doc.patch`: a one-line fix to a doc
  comment of the library, for a separate pull request; not applied here.

Nothing here uses SageMath or Python. The paper's Sage scripts were not
available when this was written; `REPORT.md` says which claims that leaves
uncovered.
