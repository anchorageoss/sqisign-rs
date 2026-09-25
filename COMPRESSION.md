# The compressed signature format

sqisign-rs carries a round-3 signature in two formats: the specification's
(200 / 306 / 406 bytes at levels I / III / V) and a compressed one of
176 / 269 / 353 bytes, 12 % smaller, which drops one of the four entries
of the basis-change matrix and lets the verifier recover it from two Weil
pairings it can compute on data it already has. **The compressed format is
not in the SQIsign specification.** It is this repository's construction
(the method of the ECLIPSE encoding in prism-rs, carried to round 3), and
a compressed signature verifies only with this crate's verifier.

## API

```rust
use sqisign_rs::{generate, CompressedSignature, Level1, Verifier};

let (pk, sk) = generate::<Level1>(&mut rng);
let sig = sk.sign(msg, &mut rng)?;
let c = sig.compress();                       // CompressedSignature<Level1>, 176 bytes
let bytes = c.to_bytes();
let c = CompressedSignature::<Level1>::from_bytes(&bytes)?;
pk.verify_compressed(msg, &c)?;               // or Verifier::verify(&pk, msg, &c)
let sig_again = c.decompress(&pk)?;           // the standard signature, needs the key
```

`sqisign_verify::compressed` has the same on slices for `no_std` users:
`compress`, `compressed_to_bytes`, `compressed_from_bytes`,
`verify_compressed`, `verify_compressed_prepared` (a `PreparedPublicKey`
does the per-key work once), `decompress`, and `compressed_bytes(params)`.
Compression needs only the signature; decompression and verification need
the public key.

## Layout

```text
A_aux | m00 | m01 | m_var | challenge | bits | hint_aux | hint_chall
```

| field | bytes (I / III / V) | content |
|---|---|---|
| `A_aux` | 82 / 128 / 168 | the auxiliary curve's Montgomery coefficient, as in the standard format |
| `m00`, `m01`, `m_var` | 3 x (25 / 38 / 50) | the low `RESPONSE_BITS` bits (196 / 302 / 400) of three matrix entries; `m_var` is `m10` when `m00` is odd, `m11` otherwise |
| `challenge` | 16 / 24 / 32 | as in the standard format |
| `bits` | 1 | bit `RESPONSE_BITS` of `m00, m01, m10, m11` in bits 0 to 3 |
| `hint_aux`, `hint_chall` | 2 | as in the standard format |

`2 fp + 3 ceil(RESPONSE_BITS / 8) + challenge + 3` = 176 / 269 / 353. The
standard format spends `4 ceil((RESPONSE_BITS + 2) / 8) + 2` on the
matrix and hints: 102 / 154 / 206 against 78 / 116 / 153 here.

The decoder rejects a wrong length, a field element at or above `p`, an
entry at or above `2^RESPONSE_BITS`, a `bits` byte with more than its four
bits, and a matrix whose first row is even (no odd pivot, so the
determinant would be even). Everything else is decided by the verifier.

## Recovery

Write `t = RESPONSE_BITS + 2`. The verifier computes the canonical bases
`(P, Q)` of `E_chall[2^t]` and `(P_aux, Q_aux)` of `E_aux[2^t]` from the
two hints, then `P' = [m00] P + [m10] Q`, `Q' = [m01] P + [m11] Q`, and
runs the `(2^RESPONSE_BITS, 2^RESPONSE_BITS)`-isogeny with kernel
`<[4](P', P_aux), [4](Q', Q_aux)>`. That kernel is maximal isotropic for
the Weil pairing on the product, so

```text
e((P', P_aux), (Q', Q_aux))^4 = 1,   i.e.   g^(4 det M) = h^-4
```

with `g = e(P, Q)` and `h = e(P_aux, Q_aux)`, both of exact order `2^t`
(the bases generate the full torsion). Hence `det M = log_g(h^-1) (mod
2^RESPONSE_BITS)`: two Weil pairings at `2^t` and one discrete logarithm in
the `2^t`-subgroup of `F_p^2^*` (`fp2_dlog_2e`, the verifier's own routine)
give the determinant modulo `2^RESPONSE_BITS`, and nothing above it, since
the relation only holds to the fourth power.

The pivot is odd because `det M` is odd. If `m00` is odd, `m11 = (det M +
m01 m10) m00^-1`; else `m01` is odd and `m10 = (m00 m11 - det M) m01^-1`,
both modulo `2^RESPONSE_BITS` (inverses by Newton iteration). The four
`bits` are set on the result, bit `RESPONSE_BITS + 1` of every entry is
left zero, and the standard verification runs on the matrix so obtained,
including the specification's range and sign check.

## Why four bits are transmitted

The isogeny depends on the matrix entries modulo `2^RESPONSE_BITS` only
(`[4] P'` does), but the chain that computes it also consumes the points
`P'`, `Q'` themselves, so the verifier's acceptance depends on the bits
above. Measured on generated signatures against this crate's verifier and
the reference implementation (the-sqisign at `nist-v3`, `sqisign_verify`
of the `broadwell` build), on all three levels:

- bit `RESPONSE_BITS + 1` of every entry is free: clearing or setting it
  changes nothing (where the specification's range check on the first
  entry allows the value);
- of the sixteen values of bit `RESPONSE_BITS` across the four entries,
  exactly two are accepted: the signer's, and the one with the first
  column `(m00, m10)` multiplied by `1 + 2^RESPONSE_BITS`; the analogous
  change of the second column is rejected;
- the determinant modulo `2^RESPONSE_BITS` from the pairings agreed with
  the signature's on every signature.

So the pairings recover the dropped entry below bit `RESPONSE_BITS`, and
the format carries bit `RESPONSE_BITS` of all four entries and drops the
bit above. Two bits would suffice in principle (one is fixed by the
pairing relation at `2^(RESPONSE_BITS + 1)`, one by the column-scaling
freedom), but they would not change the byte count, and transmitting the
four keeps the reconstruction a copy rather than a search.

The compressed verifier accepts a compressed signature exactly when the
standard verifier accepts the signature it reconstructs. The freedoms of
the standard format in those bits carry over unchanged (the column scaling
is one of them); the format adds no rule for uniqueness of the encoding,
which is not a goal of this crate (SQIsign is not strongly unforgeable).

## Measurement

`crates/sqisign-rs/tests/compressed.rs` runs the round trip at each level,
the decoder's rejections, single-bit tampering, and the acceptance survey:
`SQISIGN_SURVEY_SIGS=1000 cargo test -p sqisign-rs --release --test
compressed survey -- --nocapture` signs 1000 messages per level and checks
that every signature compresses, verifies compressed, and decompresses to
a standard signature the standard verifier accepts. Result of the
2026-09-25 run:

| level | signatures | compressed and verified | decompressed, accepted by the standard verifier | the four bits |
|---|---|---|---|---|
| I | 1000 | 1000 | 1000 | all sixteen values, 53 to 76 each |
| III | 1000 | 1000 | 1000 | all sixteen values, 53 to 72 each |
| V | 1000 | 1000 | 1000 | all sixteen values, 54 to 77 each |

The four bits are uniform: the signer's matrix leaves them random, and
they cannot be dropped without a search.

## Cost

Recovery adds two Weil pairings at `2^t` and one variable-base discrete
logarithm to a verification; the standard verification then runs
unchanged. Session of 2026-09-25 ([BENCH.md](BENCH.md)), prepared key,
millions of cycles:

| level | standard | compressed | recovery | |
|---|---|---|---|---|
| I | 11.9 | 14.1 | 2.2 | 18 % |
| III | 29.9 | 35.1 | 5.2 | 17 % |
| V | 86.2 | 97.2 | 11.0 | 13 % |

Against the C reference's `broadwell` cold verification the compressed
verification is 1.24 / 1.22 / 1.22x (the standard one 1.09 / 1.07 /
1.12x).

## History

sqisign-rs 0.4 had a compressed format for SQIsign round 2 (129 / 196 /
257 bytes, with the round-2 signature's backtracking and response-length
fields packed into a metadata byte and the hints recomputed). Round 2 is
not supported from 0.6 on, and that format has no reader here; see
[CHANGELOG.md](CHANGELOG.md).
