# Signature Compression

This document describes how SQIsign compressed signatures work: the decompression formula, precision analysis, adaptive pivot, wire formats, and approaches that were investigated and ruled out.

## The decompression formula

The dropped matrix entry is recovered using the determinant of the basis-change matrix, which is computed from Weil pairings on the challenge and auxiliary curves:

```
y = fp2_dlog_2e(omega_aux^{-1}, omega_f^{-1}, N)
det(M) = y + 2^(N-2) * u  (mod 2^N),  u in {1, 3}
det(M) mod 2^(N-1) = (y + 2^(N-2)) mod 2^(N-1)
```

Where:
- `M` is the 2x2 basis-change matrix from the signature
- `N = E_RSP + HD_EXTRA_TORSION - backtracking` is the entry width: the matrix acts on the `2^N`-torsion of the challenge curve, so every entry is an `N`-bit integer
- `pow_dim2 = E_RSP - two_resp_length - backtracking`
- `omega_f = weil(N, P_chall, Q_chall, P_chall+Q_chall, E_chall)` on the `2^N`-torsion of the challenge curve
- `omega_aux = weil(g, P_aux, Q_aux, P_aux+Q_aux, E_aux)` on the `2^g`-torsion of the auxiliary curve, where `g = pow_dim2 + HD_EXTRA_TORSION`
- `fp2_dlog_2e` solves the discrete logarithm in the `2^N`-subgroup of GF(p^2)*, so `y` is known at the full order `2^N`

**This formula is not in the SQIsign specification.** It was discovered empirically during this implementation and the relation between `y` and `det(M)` is proved in the accompanying paper. The specification defines the standard and expanded formats but does not describe a compressed format with pairing-based determinant recovery.

## Precision

For an honest signature `det(M) = y + 2^(N-2) * u (mod 2^N)` with `u` in `{1, 3}`. Since `u` is odd, reducing modulo `2^(N-1)` removes it: `det(M) mod 2^(N-1) = (y + 2^(N-2)) mod 2^(N-1)` with no hint at all. The two candidates for the dropped entry `x` (`u = 1` and `u = 3`), both solutions of `pivot * x = det(M) + (product of the other two entries) mod 2^N`, differ exactly in bit `N-1`. One hint bit `h`, bit `N-1` of the dropped entry, therefore recovers the entry completely. `h` carries the same single bit of information as `u` (`u = 1` as 0, `u = 3` as 1 would be the direct encoding), but encoding `u` itself would require the two pairings at compression time, which the entry bit avoids.

This is a stronger result than the naive analysis would suggest. The naive bound gives `pow_dim2` bits of precision from the auxiliary curve's `2^g`-torsion alone, leaving `two_resp_length + 2` unknown bits. The full `2^N`-torsion pairing on the challenge curve provides the additional `two_resp_length` bits for free, and the structure of `u` provides one more.

## Compression journey

| Size (L1) | Technique |
|---|---|
| 148 | Standard format (baseline) |
| 133 | Drop one matrix entry, Weil pairing det recovery (1-byte hint for up to `trl+2` unknown bits) |
| 132 | Precision analysis: always exactly 2 unknown bits, not `2+trl` (hint shrinks to 2 bits) |
| 130 | Recompute canonical basis hints from curves instead of storing them (drop `hint_aux` + `hint_chall`) |
| 129 | Pack n_bt, the hint and trl (4 bits) into 1 metadata byte (drop separate bt + trl bytes) |
| 129 | Determinant structure `u` in `{1, 3}`: the hint shrinks to 1 bit, which frees a third bit for n_bt (`[r:4 \| h:1 \| n_bt:3]`); n_bt up to 7 is now representable without any truncation |

## Adaptive pivot

Recovery requires dividing by a pivot value mod `2^(N-1)`, which requires the pivot to be odd (invertible mod 2). The compressor chooses which second-row entry to drop based on `M[0][0]` parity:

**M[0][0] odd** (common case, ~72% of signatures):
- Drop `M[1][1]`, store `M[1][0]` as `mat_var`
- Recover: `M[1][1] = (det + M[0][1] * M[1][0]) * M[0][0]^{-1} mod 2^(N-1)`, then set bit `N-1` to `h`

**M[0][0] even** (~28% of signatures):
- Drop `M[1][0]`, store `M[1][1]` as `mat_var`
- Recover: `M[1][0] = (M[0][0] * M[1][1] - det) * M[0][1]^{-1} mod 2^(N-1)`, then set bit `N-1` to `h`

The three transmitted entries are taken modulo `2^(N-1)` in these formulas; on the wire they stay `N` bits wide. If both `M[0][0]` and `M[0][1]` are even, no pivot exists and the signature is rejected. This cannot happen for honestly-generated signatures since the matrix has odd determinant.

## Wire formats

### All formats, all levels

| Format | L1 | L3 | L5 |
|---|---|---|---|
| Standard | 148 B | 224 B | 292 B |
| Expanded | 212 B | 316 B | 420 B |
| Compressed | 129 B | 196 B | 257 B |
| Public key | 65 B | 97 B | 129 B |

Format detection is purely length-based (each format has a unique byte count per level).

### Standard layout

```
| e_aux_a (Fp2) | bt | trl | M[0][0] | M[0][1] | M[1][0] | M[1][1] | challenge | h_aux | h_chl |
```

### Expanded layout

```
| e_aux_a (Fp2) | bt+flags | trl | challenge | P_chl_x (Fp2) | Q_chl_x (Fp2) | h_aux | h_chl |
```

The `bt+flags` byte packs: bit 7 = `kernel_is_q`, bit 6 = `pmq_sign_hint`, bits 0-5 = backtracking.

### Compressed layout

```
| e_aux_a (Fp2) | packed_meta | M[0][0] | M[0][1] | M[var] | challenge |
```

Metadata byte, most significant bit first: `[r:4 | h:1 | n_bt:3]`. Bits 7..4 hold `two_resp_length` (0..15), bit 3 holds the hint bit `h` (bit `N-1` of the dropped matrix entry), bits 2..0 hold `backtracking` (0..7). `compress()` returns an error for `backtracking >= 8` or `two_resp_length >= 16` instead of truncating. No basis hint bytes are stored.

### Field sizes

| Field | L1 | L3 | L5 |
|---|---|---|---|
| Fp2 | 64 B | 96 B | 128 B |
| Matrix entry | 16 B | 25 B | 32 B |
| Challenge | 16 B | 24 B | 32 B |

Matrix entry size is `floor((E_RSP + 9) / 8)` bytes. Challenge size is `LAMBDA / 8` bytes.

## Approaches investigated and ruled out

**Tate pairing instead of Weil for speed.** The biextension Tate pairing produces cyclic subgroup generators that are incompatible with the discrete log computation. The Tate pairing's output lives in a quotient group GF(p^2)* / (GF(p^2)*)^{2^e}, and after exponentiation to the reduced Tate pairing, the resulting root of unity does not align with the generator produced by the Weil pairing. The dlog fails. Weil is required.

**Full N-bit precision (eliminate the hint entirely).** The pairing discrete logarithm at order `2^N` determines `det(M)` modulo `2^(N-1)` and nothing more: `u` in `{1, 3}` is exactly the bit it cannot see. No amount of additional torsion data eliminates that bit without increasing the isogeny degree by a factor of 4.

**Signer retry to force a specific hint value.** The hint bit `h` is close to uniformly distributed over the KAT signatures (the per-level histograms are printed by `crates/kat/tests/kat_compressed.rs`). Eliminating it by retrying signing until `h = 0` would cost ~2x expected signing attempts. Not viable given that signing already takes ~185 ms.

**Signer-side matrix normalization for 128 bytes.** Folding the metadata (n_bt, h, r) into the matrix entries by canonicalizing the matrix form would eliminate the metadata byte entirely, reaching 128 bytes. Analysis of 100 KAT signatures: only 35/100 have `M[0][0]` with a leading 1 at the expected position, and 0/100 allow unique recovery of `(n_bt, two_resp_length)` from the matrix structure alone. The current signer does not produce canonical form. Filed upstream as SECENG-928.

**Infer n_bt and two_resp_length from matrix structure.** Attempted to recover `(n_bt, r')` by scanning for the leading 1 in `M[0][0]` and testing candidates. 0/100 KAT signatures had a unique recovery. The matrix entries do not encode these values in a recoverable way without signer cooperation.

## Path to 128 bytes

128 bytes = `Fp2 + 3 * matrix_entry + challenge` with zero metadata overhead. This requires signer-side matrix normalization: the signer must produce a canonical matrix form from which `(backtracking, two_resp_length, h)` can be uniquely inferred by the verifier. This is a signer protocol change, not a verifier-only optimization. Filed upstream as SECENG-928.
