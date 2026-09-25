# Security Policy

## Security properties

### Verification (`sqisign-verify`)

- `no_std`, heap-free, checked on `thumbv7em-none-eabihf`.
- Designed to be constant-time in the public-key and signature data it
  processes: the field kernels are straight-line (inline assembly with
  `mulx`/`adcx`/`adox` on x86-64 with BMI2 and ADX, portable kernels of
  the same limb layout otherwise, the generated radix backends
  elsewhere); conditional operations use `subtle`; the basis search is
  bounded and fails closed on crafted curves.
- The x86-64 kernels were measured with dudect in the prism-rs
  repository, which shares them (its `CT_STATUS.md`, Sections 6 and 7,
  2026-09-25): no distinguishable difference on the assembly path; the
  portable path showed a value-dependent difference on `Fp2::mul` with
  two fixed inputs and is the correctness fallback for CPUs without BMI2
  and ADX, not a path for timing claims.
- Every decoder checks the length before arithmetic and rejects a field
  element at or above `p`; the verifier rejects a non-canonical
  basis-change matrix before any curve operation (the specification's
  check).
- A formal constant-time audit has not been completed.

### Signing (variable-time by design)

- SQIsign's signing algorithm is variable-time (lattice reduction,
  integer representation and the search for the auxiliary path have
  data-dependent iteration counts), as in the C reference. Timing side
  channels exist; sign in an environment where they are acceptable.
- Signing uses fixed-precision integers (`Ibz<N>`) with a constant-time
  layer where the reference has one; no `num-bigint`.

### Zeroization

- `SigningKey<L>`, the protocol-level `SecretKey` and the compact
  `CompactSecretKey` zeroize on drop; their `Debug` output is redacted.
  `SigningKey::to_bytes` returns the secret-key encoding in a
  `Zeroizing<Vec<u8>>`, allocated once at its final size.
- Inside key generation and signing (standard and compact), the PRNG
  seed, the commitment and response ideals, the connecting quaternions,
  the challenge pulled back through the secret basis change, the pushed
  bases and the inverses of the secret matrix and norm are held in
  `Zeroizing` wrappers and cleared when they go out of scope. The
  SHAKE256 state behind the protocol's random streams is zeroized on
  drop (`sha3/zeroize`).
- Not covered: the stack temporaries inside the lattice reduction and
  the isogeny routines, the SHAKE reader's last output block, and the
  copies the compiler makes of the `Copy` integer types. There are no
  heap big-integer temporaries.

### Compact format (feature `compact`)

- Experimental, level I, not for production. Our instantiation of
  SQIsignHD on the round-3 prime; the parameters and the security posture
  are in `docs/COMPACT_R3.md`. It is not verified against an independent
  implementation beyond the SQIsignHD library's Sage verifier accepting
  our vectors at these parameters.
- The dimension-4 verifier checks the middle-codomain match of the two
  half-chains and the paper's final image condition, as the 0.4.x verifier
  did; the layer allocates (`alloc`), and its integer arithmetic on the
  response degree is variable-time on public data.

### Formats

- The specification's format and this crate's compressed format
  ([COMPRESSION.md](COMPRESSION.md)) are both accepted; a compressed
  signature verifies only with this crate. SQIsign is not strongly
  unforgeable: a valid signature has other byte strings the verifier
  accepts (the specification's own format has them), and this crate does
  not claim or add uniqueness of encodings.

## Reporting

Please do not open public issues for security vulnerabilities. Use
GitHub's private vulnerability reporting:
<https://github.com/anchorageoss/sqisign-rs/security/advisories/new>.

## Supported versions

Only the latest release is supported with security updates.
