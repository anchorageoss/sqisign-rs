# DPA countermeasures for the signing path (SECENG-1003)

Audit of physical side-channel (Differential/Simple Power Analysis, EM)
countermeasures for the elliptic-curve and theta arithmetic in the signing
path: projective-coordinate randomization and scalar blinding. Verification
operates on public data and is out of scope.

> **Status:** **implemented** behind the `dpa-protect` feature (default off).
> A significant countermeasure was **already present** unconditionally
> (theta-splitting randomization, §3). The residual gap (§4) — projective
> randomization of the `ec_biscalar_mul` basis ladders and the secondary u/v
> theta chain — is now closed under `dpa-protect`; see §6 for the realized
> design and §9 for validation results. The randomization helpers live entirely
> in the signing crate (`id2iso/dpa.rs`), so the RNG-free `verify` crate is
> untouched. With the feature off the signing path is byte-identical and the
> 300 NIST KAT vectors still reproduce exactly.

## 1. Threat model

An attacker measuring power draw or EM emanation of the signing device while it
runs the secret-dependent EC/theta arithmetic. Two classic attacks:

- **SPA / single-trace:** read the square-and-multiply (or ladder swap) pattern
  off one trace to recover scalar bits. Defeated by *coordinate randomization*
  (each run's intermediate representatives differ).
- **DPA / multi-trace:** correlate many traces (same secret, fresh randomness)
  to average out noise and recover the scalar. Defeated by *scalar blinding*
  (the scalar's bit pattern differs each run) on top of randomization.

Secrets: the secret key (endomorphism-ring generators), the per-signature
commitment randomness, and everything derived (commitment/secret ideals,
response quaternion, the secret-curve torsion bases the ladders act on).

## 2. Inventory of EC / theta scalar multiplications in signing

| Site | Operation | Point(s) | Scalar | Randomized today? |
|------|-----------|----------|--------|-------------------|
| `theta/chain.rs` `theta_chain_compute_and_eval_randomized` (called `sign_side.rs:1023`, the **response** (2,2)-isogeny) | (2,2)-isogeny theta chain | secret kernel | degree (structural) | **YES** — random normalization (§3) |
| `theta/chain.rs` `theta_chain_compute_and_eval` (called `sign_side.rs:781`, fixed-degree u/v isogeny) | (2,2)-isogeny theta chain | secret-derived | structural | **NO** |
| `ec/point.rs` `ec_biscalar_mul` via `matrix_application_even_basis` (`sign_side.rs:127/131/137`) | 2-D Montgomery ladder | secret-curve `E[2ᶠ]` basis | published `M` coeffs | NO |
| `ec/point.rs` `ec_biscalar_mul` via `verify_side.rs:47/50/70` (prover recomputes the response basis) | 2-D Montgomery ladder | secret-curve basis | published `M` coeffs | NO |
| `ec` `ec_dbl_iter_basis` (`sign_side.rs:749`) | repeated doubling (fixed count `f`) | secret-derived basis | none (pure doubling) | NO (but fixed pattern) |
| `ec/pairing.rs` `ec_dlog_2_tate` (`sign_side.rs:240/260`, `verify_side.rs:122/144`) | 2-adic discrete log | secret-curve points | — | output = published `M` (see note) |

Notes:
- The **theta chains are the dominant secret-dependent EC work** in SQIsign2D
  signing; the response chain is already randomized.
- The `ec_biscalar_mul` ladders act on **secret-curve** torsion bases but with
  the **published** basis-change matrix `M` as the scalar. So the high-value
  countermeasure there is *coordinate randomization* (protects the secret
  point/curve representation); scalar blinding is lower-value (the scalar is
  public once the signature ships) but still good practice pre-publication.
- `ec_dlog_2_tate` produces `M`/the canonical-basis hints, which are published
  in the signature; its leakage is bounded by published data (see
  CACHE_TIMING_AUDIT.md §2), so DPA on it reveals little new.

## 3. Already implemented: theta-splitting coordinate randomization

`theta/splitting.rs` (`randomize == true`, signing only) multiplies the splitting
step's base-change matrix by a randomly chosen normalization transform:

```text
secret_index = sample_random_index(rng)          // consumes the signing RNG
m_random     = select one of 6 normalization_transforms   // CONSTANT-TIME select
out_m        = m_random · out_m                  // randomizes the theta coords
```

Two good properties:

- It is the theta-model analogue of projective-coordinate randomization: the
  null point and downstream theta coordinates are in a randomized representative
  each signing, so single-trace SPA on the final splitting sees fresh values.
- The choice among the 6 transforms uses a **constant-time masked select**
  (`select_base_change_matrix` with `ct_eq`), so it is not itself a
  secret-indexed table lookup (consistent with the cache-timing audit).

It consumes the signing DRBG, so it is part of the **byte-exact, KAT-matching**
algorithm (the C reference randomizes here too) — it is not optional hardening
that can be toggled off without diverging from the KATs.

## 4. The residual gap

1. **`ec_biscalar_mul` ladders** (`matrix_application_even_basis`,
   `verify_side` basis change) run on secret-curve torsion bases with **no
   projective randomization**: the basis points enter the ladder in a fixed
   representative, so the first ladder steps have a reproducible power
   signature across signatures. Add `(X:Z) → (rX:rZ)` per input point before
   the ladder.
2. **No scalar blinding anywhere.** The `ec_biscalar_mul` scalars are the
   published `M`, so blinding is low-value there; but if blinding is added it is
   `s → s + r·2ᶠ` (the bases have order `2ᶠ`, so the result is unchanged), with
   `kbits` extended by the ~64 blinding bits (the ladder's fixed buffers,
   `k_t:[u64;16]` and `r:[u64;1024]`, accommodate this — ≤16 words, <512 bits).
3. **The `sign_side.rs:781` non-randomized theta chain** could use the
   `_randomized` variant.

## 5. Architectural blockers (why this is not a one-line change)

- **`Fp2` has no random constructor.** The `verify` crate is `no_std` and
  RNG-free by design; `Fp2::random_nonzero` would be new code (sample bytes →
  field element via the existing `Fp2::decode`, reject zero). `rand_core` is
  already a dependency of `verify`'s `theta` module, so this is feasible, but it
  extends the verifier's surface.
- **The ladder call sites don't thread an RNG.** `matrix_application_even_basis`
  and its callers take no `rng`; adding randomization means threading
  `&mut impl CryptoRngCore` through several `id2iso` signatures (or it must be
  done one frame up, where `rng` is in scope — e.g. the `dim2id2iso` functions
  at `sign_side.rs:706/803/837/1078` already carry `rng`).
- **2-D ladder correctness.** `xdblmul` consumes `P, Q, P−Q`; the three must be
  randomized to *consistent* projective representatives of the same affine
  points. `(rX:rZ) ~ (X:Z)` is valid per point, but the optimized differential
  addition must tolerate a projective difference point — to be confirmed, with a
  **sign→verify round-trip test (randomization ON)** as the correctness gate.

## 6. Realized design (implemented, gated)

Implemented in `crates/sqisign-rs/src/id2iso/dpa.rs` plus RNG plumbing:

1. **No `verify` changes.** `Fp2` and `EcPoint` already expose public fields
   (`re`/`im`, `x`/`z`) and `Fp::decode_reduce`, so the random scaling factor is
   built **in the signing crate** — the RNG-free `verify` crate is untouched
   (resolving blocker §5's "new `Fp2` random in `verify`" concern). The helper
   `fp2_random_nonzero` reduces fresh random bytes mod `p` per component
   (no byte-boundary rejection bias) and rejects the zero element.
2. `randomize_point_projective(&mut EcPoint, rng)` does `x*=r; z*=r`;
   `maybe_randomize_basis(&mut EcBasis, rng)` randomizes `P`, `Q`, and `P-Q`
   with **independent** factors (valid: the x-only ladder is homogeneous in each
   argument — confirmed against `xadd`/`xdbladd`, both fully projective in the
   difference point).
3. Feature `dpa-protect` on `sqisign-rs`, **default off** so the KATs stay
   DRBG-deterministic (`cargo test` runs with it off); production enables it.
   The existing theta randomization (§3) is unconditional regardless, so the
   dominant path is protected even with `dpa-protect` off.
4. `maybe_randomize_basis` is a **no-op** when the feature is off (it still
   takes the `&mut` borrows, so call-site `mut` bindings stay justified and no
   randomness is drawn). `rng` is threaded through `matrix_application_even_basis`,
   `endomorphism_application_even_basis` (and its three callers), and
   `compute_and_set_basis_change_matrix` ← `protocols_sign`. `sign_side.rs:781`
   switches to `theta_chain_compute_and_eval_randomized` under the flag.
5. **Scalar blinding was not added.** The `ec_biscalar_mul` scalars are the
   published basis-change matrix `M`, so blinding is low-value (§4 note 2); the
   high-value countermeasure there is coordinate randomization, which is done.
6. **`verify_side.rs` ladders excluded.** Its `matrix_application_even_basis` /
   `ec_biscalar_mul_ibz` have **no callers** in `protocols_sign` (only
   `mp_invert_matrix` from that file is used, by a cross-validation test), so
   they are not on the active signing path and were left unchanged.

## 9. Validation results

| Check | `dpa-protect` off (default) | `dpa-protect` on |
|-------|-----------------------------|------------------|
| `cargo build -p sqisign-rs --release` | pass | pass |
| `cargo clippy -p sqisign-rs --release -- -D warnings` | pass | pass |
| `no_std` (`--no-default-features`) + `x86_64-unknown-none` | pass | pass |
| NIST KAT vectors (keygen/sign/verify × L1/L3/L5) | **300/300 byte-exact** | n/a (randomness diverges by design) |
| `sign_roundtrip` (sign → verify) | 7/7 | **7/7** (math preserved) |

The feature-off KAT pass proves the change is byte-transparent; the feature-on
round-trip pass proves the randomized representatives still produce valid
signatures (the affine result is unchanged).

## 7. Cost

Per protected ladder: one `Fp2::random_nonzero` + 2 `Fp2` muls (randomization),
optionally +~64 ladder steps (blinding). The existing theta randomization adds
one masked select over 6 matrices + one matrix multiply per signature. Total EC
+ theta work is ~8 ms (CT profiling); these additions are well under 1%.

## 8. Conclusion

The dominant secret-dependent isogeny computation (the response (2,2)-chain) is
**already coordinate-randomized**, with a constant-time transform selection. The
residual gap is the `ec_biscalar_mul` basis-construction ladders (and one
secondary theta chain), where the secret-sensitive value is the curve/point
representation (randomization) more than the scalar (which is the published
`M`). Closing it is the §6 plan: feasible and cheap, but it touches the RNG-free
`verify` crate and several `id2iso` signatures and must be round-trip-validated,
so it is specified here for a deliberate follow-up rather than applied blind.
