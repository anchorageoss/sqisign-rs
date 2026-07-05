//! Canonical (non-malleable) basis-change matrix encoding.
//!
//! By default the signer emits the canonical representative of the challenge
//! basis-change matrix and the verifier rejects the non-canonical (negated)
//! twin, closing the signature malleability described in ePrint 2026/1305.
//! These tests exercise that default behaviour and are compiled out under the
//! `kat-compat` feature, which restores the legacy (malleable) encoding for
//! byte-exact C reference KAT compatibility.

#![cfg(not(feature = "kat-compat"))]

use hybrid_array::typenum::Unsigned;
use sqisign_rs::keygen::keypair;
use sqisign_rs::params::Level1;
use sqisign_rs::sign::sign;
use sqisign_rs::verify::is_canonical_basis_change_matrix;
use sqisign_rs::{SecurityLevel, Signature, Verifier};

type L1 = Level1;

/// `(matrix_start_offset, bytes_per_entry)` in the standard wire encoding:
/// the matrix follows the `E_aux_A` field and the two metadata bytes
/// (`backtracking`, `two_resp_length`).
fn matrix_layout() -> (usize, usize) {
    let fp2 = <L1 as SecurityLevel>::Fp2EncodedBytes::USIZE;
    let entry_bytes = (<L1 as SecurityLevel>::E_RSP as usize + 9) / 8;
    (fp2 + 2, entry_bytes)
}

/// Negate a little-endian integer stored in `buf` modulo `2^n_bits`, in place.
/// Computes the two's complement over the full byte window (which equals
/// `2^(8·len) − value`) and then masks off the bits at or above `n_bits`,
/// yielding `(2^n_bits − value) mod 2^n_bits`. A zero value maps to zero.
fn negate_le_mod_2n(buf: &mut [u8], n_bits: u32) {
    let mut carry = 1u16;
    for b in buf.iter_mut() {
        let v = (!*b) as u16 + carry;
        *b = v as u8;
        carry = v >> 8;
    }

    let len = buf.len();
    let n = n_bits as usize;
    let full = n / 8;
    let rem = n % 8;
    if rem != 0 {
        if full < len {
            buf[full] &= (1u8 << rem) - 1;
        }
        for b in buf[core::cmp::min(full + 1, len)..].iter_mut() {
            *b = 0;
        }
    } else {
        for b in buf[core::cmp::min(full, len)..].iter_mut() {
            *b = 0;
        }
    }
}

/// Negate every entry of the basis-change matrix encoded in a standard
/// signature's wire bytes, modulo `2^N` with `N = E_RSP + HD_EXTRA_TORSION −
/// backtracking`. This produces the malleable twin an attacker could forge.
fn negate_matrix_in_bytes(bytes: &mut [u8]) {
    let fp2 = <L1 as SecurityLevel>::Fp2EncodedBytes::USIZE;
    let backtracking = bytes[fp2] as u32;
    let n_bits = <L1 as SecurityLevel>::E_RSP + sqisign_rs::theta::HD_EXTRA_TORSION - backtracking;

    let (mat_start, entry_bytes) = matrix_layout();
    for k in 0..4 {
        let off = mat_start + k * entry_bytes;
        negate_le_mod_2n(&mut bytes[off..off + entry_bytes], n_bits);
    }
}

#[test]
fn test_canonical_rejects_negated() {
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<L1>(&mut rng);
    let msg = b"canonical matrix malleability test";

    let sig = sign::<L1>(&sk, &pk, msg, &mut rng).expect("signing must succeed");

    // The honestly-produced (canonical) signature verifies.
    assert!(
        pk.verify(msg, &sig).is_ok(),
        "canonical signature must verify"
    );

    // Its negated twin is still cryptographically valid but non-canonical, so
    // the default verifier must reject it.
    let mut bytes = sig.to_bytes().as_slice().to_vec();
    negate_matrix_in_bytes(&mut bytes);
    let negated =
        Signature::<L1>::from_bytes(&bytes).expect("negated matrix must still deserialize");

    assert!(
        pk.verify(msg, &negated).is_err(),
        "non-canonical (negated) signature must be rejected"
    );
}

#[test]
fn test_canonical_signatures_are_canonical() {
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<L1>(&mut rng);
    let msg = b"canonical matrix normalization test";

    let sig = sign::<L1>(&sk, &pk, msg, &mut rng).expect("signing must succeed");

    // Re-parse to inspect the matrix through the public predicate.
    let parsed =
        Signature::<L1>::from_bytes(sig.to_bytes().as_slice()).expect("signature must round-trip");

    assert!(
        is_canonical_basis_change_matrix(&parsed),
        "signer must emit a canonical basis-change matrix"
    );
}
