# Security Policy

## Security Properties

### Verification

- The verification path is designed to be constant-time, using the `subtle` crate for conditional operations (select, cswap)
- No BigInt operations in the verification path
- The underlying Fp arithmetic was translated from modarith-generated C code which is designed for constant-time execution
- A formal constant-time audit has not yet been completed; data-dependent branches may exist in the Rust translation
- Verified `no_std` on `thumbv7em-none-eabihf`

### Signing (variable-time by design)

- SQIsign's signing algorithm is inherently variable-time
- LLL, represent_integer, and find_uv have data-dependent iteration counts
- This is shared with the C reference implementation
- Timing side channels exist; use appropriate deployment mitigations

### Zeroization

- `SecretKey` implements `ZeroizeOnDrop` (zeroed on drop)
- Intermediate signing values (`QuatLeftIdeal`, `QuatAlgElem`, torsion bases) are explicitly zeroized after use
- Known limitation: `num-bigint`'s `BigInt` does not expose its backing `Vec<u64>`, so `ibz_zeroize()` replaces the value with zero but cannot scrub the freed heap allocation
- Optional: `sqisign-alloc::ZeroizingAllocator` zeros ALL heap memory on deallocation (measured overhead: <1%)
- The C reference (GMP) has the same heap residue limitation

## Supported Versions

Only the latest release is supported with security updates.
