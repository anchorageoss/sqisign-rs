//!
//! The prime is `p = 65 * 2^376 - 1`. Field elements are 7-limb arrays
//! `[u64; 7]` storing values in unsaturated radix `2^55`: each limb
//! holds a roughly-55-bit quantity, leaving 9 carry bits at the top.
//! Arithmetic operates in an internal Montgomery form; canonical
//! little-endian byte encoding (via `fp_encode` / `fp_decode`)
//! returns the standard integer representation in `[0, p)`.
//!
//! The free functions below implement the limb-level primitives;
//! the [`FpBackend`] trait impl for [`Level3`] at the bottom of the
//! file wires them into the generic [`super::Fp`] API.

use super::FpBackend;
use crate::params::Level3;
use hybrid_array::Array;
use subtle::Choice;

/// Number of 64-bit limbs in an `Fp` element.
pub const NLIMBS: usize = 7;

/// Bit width of each unsaturated limb.
pub const RADIX: u32 = 55;

/// `2^RADIX - 1`: keeps the low 55 bits of a limb.
pub const MASK: u64 = (1u64 << RADIX) - 1;

/// The contribution of `p + 1 = 65 * 2^376` to limb 6 in radix-`2^55`:
/// `65 * 2^376 / 2^(55*6) = 65 * 2^46`. Used as the Montgomery folding
/// constant inside `modmul` and `modsqr`.
pub const P6: u64 = 0x10400000000000u64;

/// `2 * P6`, used inside `modadd`, `modsub`, and `modneg` to
/// keep partial sums in `[-p, p)` before the final reduction step.
pub const TWO_P6: u64 = 0x20800000000000u64;

/// `0` as 7 radix-`2^55` limbs.
pub const ZERO: [u64; NLIMBS] = [0; NLIMBS];

/// Internal Montgomery form of `1 mod p`.
pub const ONE: [u64; NLIMBS] = [
    0x0000000000000007,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x000e400000000000,
];

/// Internal Montgomery form of 2⁻¹ mod p.
pub const TWO_INV: [u64; NLIMBS] = [
    0x0000000000000003,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x000f400000000000,
];

/// Internal Montgomery form of 3⁻¹ mod p.
pub const THREE_INV: [u64; NLIMBS] = [
    0x0055555555555557,
    0x002aaaaaaaaaaaaa,
    0x0055555555555555,
    0x002aaaaaaaaaaaaa,
    0x0055555555555555,
    0x002aaaaaaaaaaaaa,
    0x000f955555555555,
];

/// Internal Montgomery form of `2^384 mod p`. Used by
/// `fp_decode_reduce` when folding 48-byte input blocks.
pub const R2: [u64; NLIMBS] = [
    0x0007e07e07e07e26,
    0x007c0fc0fc0fc0fc,
    0x0001f81f81f81f81,
    0x003f03f03f03f03f,
    0x00607e07e07e07e0,
    0x000fc0fc0fc0fc0f,
    0x000e9f81f81f81f8,
];

/// Conversion constant for `nres`: multiplying by `R2_NRES` and
/// dividing by Montgomery `R` takes a canonical integer into internal
/// Montgomery form (`R2_NRES = R^2 mod p` in the internal Montgomery
/// scheme).
pub const R2_NRES: [u64; NLIMBS] = [
    0xfc0fc0fc0fc4d,
    0x781f81f81f81f8,
    0x3f03f03f03f03,
    0x7e07e07e07e07e,
    0x40fc0fc0fc0fc0,
    0x1f81f81f81f81f,
    0xcff03f03f03f0,
];

/// Propagate carries between limbs, leaving each of limbs 0..5 in
/// `[0, 2^55)` and limb 6 holding the accumulated high bits. Returns
/// `0xFFFF_FFFF_FFFF_FFFF` if the value (interpreted as two's-complement
/// via the top bits of limb 6) is negative, else `0`.
#[inline]
pub(crate) fn prop(n: &mut [u64; NLIMBS]) -> u64 {
    let mask = MASK;
    let mut carry: i64 = n[0] as i64;
    carry >>= RADIX;
    n[0] &= mask;
    for limb in n.iter_mut().take(6).skip(1) {
        carry = carry.wrapping_add(*limb as i64);
        *limb = (carry as u64) & mask;
        carry >>= RADIX;
    }
    n[6] = n[6].wrapping_add(carry as u64);
    // Sign mask: 0 if non-negative, all-ones if negative.
    // Limb 6 is not masked to radix width; the sign occupies the high bits.
    let sign = (n[6] >> 1) >> 62;
    0u64.wrapping_sub(sign)
}

/// [`prop`] followed by a conditional add of `p` if the propagation
/// detected a negative intermediate, then a second [`prop`] pass.
/// Returns `1` if the value was originally negative (and was just
/// fixed up), `0` otherwise.
#[inline]
pub(crate) fn flatten(n: &mut [u64; NLIMBS]) -> u32 {
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(1 & carry);
    n[6] = n[6].wrapping_add(0x10400000000000u64 & carry);
    let _ = prop(n);
    (carry & 1) as u32
}

/// Final subtract of `p`: tries `n -= p` and adds `p` back if the
/// result went negative. Returns `1` if the original value was `< p`
/// (and was preserved), `0` if it was `>= p` (and was reduced).
#[inline]
pub(crate) fn modfsb(n: &mut [u64; NLIMBS]) -> u32 {
    n[0] = n[0].wrapping_add(1);
    n[6] = n[6].wrapping_sub(0x10400000000000u64);
    flatten(n)
}

/// `n <- a + b mod 2p` (lazily reduced).
#[inline]
pub(crate) fn modadd(n: &mut [u64; NLIMBS], a: &[u64; NLIMBS], b: &[u64; NLIMBS]) {
    n[0] = a[0].wrapping_add(b[0]);
    n[1] = a[1].wrapping_add(b[1]);
    n[2] = a[2].wrapping_add(b[2]);
    n[3] = a[3].wrapping_add(b[3]);
    n[4] = a[4].wrapping_add(b[4]);
    n[5] = a[5].wrapping_add(b[5]);
    n[6] = a[6].wrapping_add(b[6]);
    n[0] = n[0].wrapping_add(2);
    n[6] = n[6].wrapping_sub(TWO_P6);
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2 & carry);
    n[6] = n[6].wrapping_add(TWO_P6 & carry);
    let _ = prop(n);
}

/// `n <- a - b mod 2p` (lazily reduced).
#[inline]
pub(crate) fn modsub(n: &mut [u64; NLIMBS], a: &[u64; NLIMBS], b: &[u64; NLIMBS]) {
    n[0] = a[0].wrapping_sub(b[0]);
    n[1] = a[1].wrapping_sub(b[1]);
    n[2] = a[2].wrapping_sub(b[2]);
    n[3] = a[3].wrapping_sub(b[3]);
    n[4] = a[4].wrapping_sub(b[4]);
    n[5] = a[5].wrapping_sub(b[5]);
    n[6] = a[6].wrapping_sub(b[6]);
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2 & carry);
    n[6] = n[6].wrapping_add(TWO_P6 & carry);
    let _ = prop(n);
}

/// `n <- -b mod 2p` (lazily reduced).
#[inline]
pub(crate) fn modneg(n: &mut [u64; NLIMBS], b: &[u64; NLIMBS]) {
    n[0] = 0u64.wrapping_sub(b[0]);
    n[1] = 0u64.wrapping_sub(b[1]);
    n[2] = 0u64.wrapping_sub(b[2]);
    n[3] = 0u64.wrapping_sub(b[3]);
    n[4] = 0u64.wrapping_sub(b[4]);
    n[5] = 0u64.wrapping_sub(b[5]);
    n[6] = 0u64.wrapping_sub(b[6]);
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2 & carry);
    n[6] = n[6].wrapping_add(TWO_P6 & carry);
    let _ = prop(n);
}

/// `c <- a * b mod 2p`. Schoolbook 7x7 multiplication in radix `2^55`
/// with interleaved Montgomery folding: each newly-emitted low limb
/// `v_i` is added back at limb position `i + 6` as `v_i * P6`, which
/// represents `v_i * 65 * 2^376 = v_i * (p + 1) = v_i (mod p)` and
/// effectively divides the final result by Montgomery `R = 2^385`.
#[inline]
pub(crate) fn modmul(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS], b: &[u64; NLIMBS]) {
    let mask = MASK;
    let p6: u128 = P6 as u128;

    // t accumulates the partial sum at the current limb position. After
    // emitting each limb (`t & mask`) we shift right by RADIX.
    let mut t: u128;

    t = (a[0] as u128) * (b[0] as u128);
    let v0 = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[0] as u128) * (b[1] as u128))
        .wrapping_add((a[1] as u128) * (b[0] as u128));
    let v1 = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[0] as u128) * (b[2] as u128))
        .wrapping_add((a[1] as u128) * (b[1] as u128))
        .wrapping_add((a[2] as u128) * (b[0] as u128));
    let v2 = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[0] as u128) * (b[3] as u128))
        .wrapping_add((a[1] as u128) * (b[2] as u128))
        .wrapping_add((a[2] as u128) * (b[1] as u128))
        .wrapping_add((a[3] as u128) * (b[0] as u128));
    let v3 = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[0] as u128) * (b[4] as u128))
        .wrapping_add((a[1] as u128) * (b[3] as u128))
        .wrapping_add((a[2] as u128) * (b[2] as u128))
        .wrapping_add((a[3] as u128) * (b[1] as u128))
        .wrapping_add((a[4] as u128) * (b[0] as u128));
    let v4 = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[0] as u128) * (b[5] as u128))
        .wrapping_add((a[1] as u128) * (b[4] as u128))
        .wrapping_add((a[2] as u128) * (b[3] as u128))
        .wrapping_add((a[3] as u128) * (b[2] as u128))
        .wrapping_add((a[4] as u128) * (b[1] as u128))
        .wrapping_add((a[5] as u128) * (b[0] as u128));
    let v5 = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[0] as u128) * (b[6] as u128))
        .wrapping_add((a[1] as u128) * (b[5] as u128))
        .wrapping_add((a[2] as u128) * (b[4] as u128))
        .wrapping_add((a[3] as u128) * (b[3] as u128))
        .wrapping_add((a[4] as u128) * (b[2] as u128))
        .wrapping_add((a[5] as u128) * (b[1] as u128))
        .wrapping_add((a[6] as u128) * (b[0] as u128))
        .wrapping_add((v0 as u128) * p6);
    let v6 = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[1] as u128) * (b[6] as u128))
        .wrapping_add((a[2] as u128) * (b[5] as u128))
        .wrapping_add((a[3] as u128) * (b[4] as u128))
        .wrapping_add((a[4] as u128) * (b[3] as u128))
        .wrapping_add((a[5] as u128) * (b[2] as u128))
        .wrapping_add((a[6] as u128) * (b[1] as u128))
        .wrapping_add((v1 as u128) * p6);
    c[0] = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[2] as u128) * (b[6] as u128))
        .wrapping_add((a[3] as u128) * (b[5] as u128))
        .wrapping_add((a[4] as u128) * (b[4] as u128))
        .wrapping_add((a[5] as u128) * (b[3] as u128))
        .wrapping_add((a[6] as u128) * (b[2] as u128))
        .wrapping_add((v2 as u128) * p6);
    c[1] = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[3] as u128) * (b[6] as u128))
        .wrapping_add((a[4] as u128) * (b[5] as u128))
        .wrapping_add((a[5] as u128) * (b[4] as u128))
        .wrapping_add((a[6] as u128) * (b[3] as u128))
        .wrapping_add((v3 as u128) * p6);
    c[2] = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[4] as u128) * (b[6] as u128))
        .wrapping_add((a[5] as u128) * (b[5] as u128))
        .wrapping_add((a[6] as u128) * (b[4] as u128))
        .wrapping_add((v4 as u128) * p6);
    c[3] = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[5] as u128) * (b[6] as u128))
        .wrapping_add((a[6] as u128) * (b[5] as u128))
        .wrapping_add((v5 as u128) * p6);
    c[4] = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[6] as u128) * (b[6] as u128))
        .wrapping_add((v6 as u128) * p6);
    c[5] = (t as u64) & mask;
    t >>= RADIX;

    c[6] = t as u64;
}

/// `c <- a * a mod 2p` (specialized squaring). Each cross-term
/// `a_i * a_j` (`i != j`) is computed once and doubled, saving roughly
/// half the partial-product work compared to a general multiply.
#[inline]
pub(crate) fn modsqr(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS]) {
    let mask = MASK;
    let p6: u128 = P6 as u128;

    let mut t: u128;
    let mut tot: u128;

    tot = (a[0] as u128) * (a[0] as u128);
    t = tot;
    let v0 = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[0] as u128) * (a[1] as u128);
    tot = tot.wrapping_mul(2);
    t = t.wrapping_add(tot);
    let v1 = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[0] as u128) * (a[2] as u128);
    tot = tot.wrapping_mul(2);
    tot = tot.wrapping_add((a[1] as u128) * (a[1] as u128));
    t = t.wrapping_add(tot);
    let v2 = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[0] as u128) * (a[3] as u128);
    tot = tot.wrapping_add((a[1] as u128) * (a[2] as u128));
    tot = tot.wrapping_mul(2);
    t = t.wrapping_add(tot);
    let v3 = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[0] as u128) * (a[4] as u128);
    tot = tot.wrapping_add((a[1] as u128) * (a[3] as u128));
    tot = tot.wrapping_mul(2);
    tot = tot.wrapping_add((a[2] as u128) * (a[2] as u128));
    t = t.wrapping_add(tot);
    let v4 = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[0] as u128) * (a[5] as u128);
    tot = tot.wrapping_add((a[1] as u128) * (a[4] as u128));
    tot = tot.wrapping_add((a[2] as u128) * (a[3] as u128));
    tot = tot.wrapping_mul(2);
    t = t.wrapping_add(tot);
    let v5 = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[0] as u128) * (a[6] as u128);
    tot = tot.wrapping_add((a[1] as u128) * (a[5] as u128));
    tot = tot.wrapping_add((a[2] as u128) * (a[4] as u128));
    tot = tot.wrapping_mul(2);
    tot = tot.wrapping_add((a[3] as u128) * (a[3] as u128));
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v0 as u128) * p6);
    let v6 = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[1] as u128) * (a[6] as u128);
    tot = tot.wrapping_add((a[2] as u128) * (a[5] as u128));
    tot = tot.wrapping_add((a[3] as u128) * (a[4] as u128));
    tot = tot.wrapping_mul(2);
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v1 as u128) * p6);
    c[0] = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[2] as u128) * (a[6] as u128);
    tot = tot.wrapping_add((a[3] as u128) * (a[5] as u128));
    tot = tot.wrapping_mul(2);
    tot = tot.wrapping_add((a[4] as u128) * (a[4] as u128));
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v2 as u128) * p6);
    c[1] = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[3] as u128) * (a[6] as u128);
    tot = tot.wrapping_add((a[4] as u128) * (a[5] as u128));
    tot = tot.wrapping_mul(2);
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v3 as u128) * p6);
    c[2] = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[4] as u128) * (a[6] as u128);
    tot = tot.wrapping_mul(2);
    tot = tot.wrapping_add((a[5] as u128) * (a[5] as u128));
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v4 as u128) * p6);
    c[3] = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[5] as u128) * (a[6] as u128);
    tot = tot.wrapping_mul(2);
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v5 as u128) * p6);
    c[4] = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[6] as u128) * (a[6] as u128);
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v6 as u128) * p6);
    c[5] = (t as u64) & mask;
    t >>= RADIX;

    c[6] = t as u64;
}

/// `c <- a`.
#[inline]
pub(crate) fn modcpy(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS]) {
    *c = *a;
}

/// Square `a` in-place `n` times.
#[inline]
pub(crate) fn modnsqr(a: &mut [u64; NLIMBS], n: u32) {
    for _ in 0..n {
        let mut tmp = [0u64; NLIMBS];
        modsqr(&mut tmp, a);
        *a = tmp;
    }
}

/// Square root progenitor: `z <- w^((p-3)/4) mod p` via a fixed
/// addition chain whose structure encodes the binary representation
/// of `(p-3)/4` for this prime.
#[inline]
pub(crate) fn modpro(z: &mut [u64; NLIMBS], w: &[u64; NLIMBS]) {
    let mut x = [0u64; NLIMBS];
    let mut t0 = [0u64; NLIMBS];
    let mut t1: [u64; NLIMBS];
    let mut t2 = [0u64; NLIMBS];
    let mut t3 = [0u64; NLIMBS];
    let mut t4 = [0u64; NLIMBS];
    let mut t5 = [0u64; NLIMBS];

    modcpy(&mut x, w);
    modsqr(z, &x);
    {
        let z_copy = *z;
        modsqr(&mut t0, &z_copy);
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &x, &t0);
        t1 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, z, &t1);
        *z = tmp;
    }
    {
        let z_copy = *z;
        modsqr(&mut t0, &z_copy);
    }
    {
        modsqr(&mut t3, &t0);
    }
    {
        modsqr(&mut t4, &t3);
    }
    {
        modsqr(&mut t2, &t4);
    }
    modcpy(&mut t5, &t2);
    modnsqr(&mut t5, 3);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t2, &t5);
        t2 = tmp;
    }
    modcpy(&mut t5, &t2);
    modnsqr(&mut t5, 6);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t2, &t5);
        t2 = tmp;
    }
    modcpy(&mut t5, &t2);
    modnsqr(&mut t5, 2);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t4, &t5);
        t5 = tmp;
    }
    modnsqr(&mut t5, 13);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t2, &t5);
        t2 = tmp;
    }
    modcpy(&mut t5, &t2);
    modnsqr(&mut t5, 2);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t4, &t5);
        t4 = tmp;
    }
    modnsqr(&mut t4, 28);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t2, &t4);
        t2 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modsqr(&mut tmp, &t2);
        t4 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t3, &t4);
        t3 = tmp;
    }
    modnsqr(&mut t3, 59);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t2, &t3);
        t2 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t1, &t2);
        t1 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, z, &t1);
        *z = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t0, z);
        t0 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t1, &t0);
        t1 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modsqr(&mut tmp, &t1);
        t2 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t1, &t2);
        t2 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modsqr(&mut tmp, &t2);
        t2 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t1, &t2);
        t2 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t0, &t2);
        t0 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, z, &t0);
        *z = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modsqr(&mut tmp, z);
        t2 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, z, &t2);
        t2 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t0, &t2);
        t0 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t1, &t0);
        t1 = tmp;
    }
    modcpy(&mut t2, &t1);
    modnsqr(&mut t2, 128);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t1, &t2);
        t1 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t0, &t1);
        t0 = tmp;
    }
    modnsqr(&mut t0, 125);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, z, &t0);
        *z = tmp;
    }
}

/// `modinv`: `z <- 1 / x mod p`. If `h` is `Some(h)` the precomputed
/// progenitor `x^((p-3)/4)` is used; otherwise computed via `modpro`.
#[inline]
pub(crate) fn modinv(z: &mut [u64; NLIMBS], x: &[u64; NLIMBS], h: Option<&[u64; NLIMBS]>) {
    let mut s = [0u64; NLIMBS];
    let mut t = [0u64; NLIMBS];
    match h {
        None => modpro(&mut t, x),
        Some(h) => modcpy(&mut t, h),
    }
    modcpy(&mut s, x);
    modnsqr(&mut t, 2);
    modmul(z, &s, &t);
}

/// Convert from canonical integer form to internal Montgomery form
/// (multiplies by `R2_NRES` and divides by Montgomery `R`).
#[inline]
pub(crate) fn nres(n: &mut [u64; NLIMBS], m: &[u64; NLIMBS]) {
    modmul(n, m, &R2_NRES);
}

/// Convert from internal Montgomery form back to canonical integer
/// form, fully reducing modulo `p`.
#[inline]
pub(crate) fn redc(m: &mut [u64; NLIMBS], n: &[u64; NLIMBS]) {
    let mut c = [0u64; NLIMBS];
    c[0] = 1;
    modmul(m, n, &c);
    let _ = modfsb(m);
}

/// Returns `1` if `a == 0 mod p`, `0` otherwise.
#[inline]
pub(crate) fn modis0(a: &[u64; NLIMBS]) -> u32 {
    let mut c = [0u64; NLIMBS];
    redc(&mut c, a);
    let mut d: u64 = 0;
    for limb in c.iter() {
        d |= *limb;
    }
    ((1u64) & ((d.wrapping_sub(1)) >> RADIX)) as u32
}

/// Returns `1` if `a == 1 mod p`, `0` otherwise.
#[inline]
pub(crate) fn modis1(a: &[u64; NLIMBS]) -> u32 {
    let mut c = [0u64; NLIMBS];
    redc(&mut c, a);
    let mut d: u64 = 0;
    for limb in c.iter().skip(1) {
        d |= *limb;
    }
    let c0 = c[0];
    ((1u64) & ((d.wrapping_sub(1)) >> RADIX) & (((c0 ^ 1).wrapping_sub(1)) >> RADIX)) as u32
}

/// Returns `1` if `a == b mod p`, `0` otherwise. Constant-time: folds
/// the per-limb differences into a single bit using masked XOR.
#[inline]
pub(crate) fn modcmp(a: &[u64; NLIMBS], b: &[u64; NLIMBS]) -> u32 {
    let mut c = [0u64; NLIMBS];
    let mut d = [0u64; NLIMBS];
    redc(&mut c, a);
    redc(&mut d, b);
    let mut eq: u64 = 1;
    for i in 0..NLIMBS {
        eq &= ((c[i] ^ d[i]).wrapping_sub(1)) >> RADIX;
    }
    eq &= 1;
    eq as u32
}

/// Returns `1` if `x` is a quadratic residue (or zero), `0`
/// otherwise. If `h` is `Some`, it is interpreted as the precomputed
/// progenitor `x^((p-3)/4)`; otherwise the progenitor is computed.
#[inline]
pub(crate) fn modqr(x: &[u64; NLIMBS], h: Option<&[u64; NLIMBS]>) -> u32 {
    let mut r = [0u64; NLIMBS];
    match h {
        None => {
            modpro(&mut r, x);
            let mut r2 = [0u64; NLIMBS];
            modsqr(&mut r2, &r);
            r = r2;
        }
        Some(h) => {
            modsqr(&mut r, h);
        }
    }
    let mut r2 = [0u64; NLIMBS];
    modmul(&mut r2, &r, x);
    let r = r2;
    modis1(&r) | modis0(x)
}

/// `r <- sqrt(x) mod p`. Since `p = 3 mod 4`, the square root is
/// `x^((p+1)/4) = x * x^((p-3)/4)`, so the progenitor (cached as `h`
/// when available) is multiplied by `x` to get the answer. Output is
/// well-defined modulo a sign.
#[inline]
pub(crate) fn modsqrt(r: &mut [u64; NLIMBS], x: &[u64; NLIMBS], h: Option<&[u64; NLIMBS]>) {
    let mut y = [0u64; NLIMBS];
    let mut s = [0u64; NLIMBS];
    match h {
        None => modpro(&mut y, x),
        Some(h) => modcpy(&mut y, h),
    }
    modmul(&mut s, &y, x);
    modcpy(r, &s);
}

/// Set `a` to the small integer `x` in internal Montgomery form.
#[inline]
pub(crate) fn modint(a: &mut [u64; NLIMBS], x: u64) {
    a[0] = x;
    a[1] = 0;
    a[2] = 0;
    a[3] = 0;
    a[4] = 0;
    a[5] = 0;
    a[6] = 0;
    let a_copy = *a;
    nres(a, &a_copy);
}

/// `c <- a * x mod 2p` for a small integer `x`.
#[inline]
pub(crate) fn modmli(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS], x: u64) {
    let mut t = [0u64; NLIMBS];
    modint(&mut t, x);
    modmul(c, a, &t);
}

/// Shift `a` left by `n < 55` bits (per-limb, with carry into the
/// next limb).
#[inline]
pub(crate) fn modshl(a: &mut [u64; NLIMBS], n: u32) {
    a[6] = (a[6] << n) | (a[5] >> (RADIX - n));
    for i in (1..=5).rev() {
        a[i] = ((a[i] << n) & MASK) | (a[i - 1] >> (RADIX - n));
    }
    a[0] = (a[0] << n) & MASK;
}

/// Shift `a` right by `n < 55` bits, returning the low `n`
/// shifted-out bits.
#[inline]
pub(crate) fn modshr(a: &mut [u64; NLIMBS], n: u32) -> u64 {
    let r = a[0] & ((1u64 << n) - 1);
    for i in 0..6 {
        a[i] = (a[i] >> n) | ((a[i + 1] << (RADIX - n)) & MASK);
    }
    a[6] >>= n;
    r
}

/// Constant-time conditional swap of `f` and `g` if `b == 1`, no-op
/// if `b == 0`.
#[inline]
pub(crate) fn modcsw(b: u64, g: &mut [u64; NLIMBS], f: &mut [u64; NLIMBS]) {
    let r: u64 = 0x3cc3_c33c_5aa5_a55a;
    let c0 = (1u64.wrapping_sub(b)).wrapping_add(r);
    let c1 = b.wrapping_add(r);
    for i in 0..NLIMBS {
        let s = g[i];
        let t = f[i];
        let w = r.wrapping_mul(t.wrapping_add(s));
        let new_f = c0
            .wrapping_mul(t)
            .wrapping_add(c1.wrapping_mul(s))
            .wrapping_sub(w);
        let new_g = c0
            .wrapping_mul(s)
            .wrapping_add(c1.wrapping_mul(t))
            .wrapping_sub(w);
        f[i] = new_f;
        g[i] = new_g;
    }
}

/// Bytes per encoded Fp element for Level 3.
const ENCODED_BYTES: usize = 48;

/// Serialize `a` (in internal Montgomery form) to its canonical
/// 48-byte little-endian representation: convert back to canonical
/// integer form via [`redc`], then peel off one byte at a time.
#[inline]
pub(crate) fn fp_encode(out: &mut [u8], a: &[u64; NLIMBS]) {
    debug_assert!(out.len() >= ENCODED_BYTES);
    let mut c = [0u64; NLIMBS];
    redc(&mut c, a);
    for byte in out.iter_mut().take(ENCODED_BYTES) {
        *byte = (c[0] & 0xff) as u8;
        let _ = modshr(&mut c, 8);
    }
}

/// Parse 48 canonical little-endian bytes into an Fp element.
/// Returns `0xFFFF_FFFF` on in-range input, `0` otherwise; on failure
/// `out` is zeroed.
#[inline]
pub(crate) fn fp_decode(out: &mut [u64; NLIMBS], bytes: &[u8]) -> u32 {
    if bytes.len() < ENCODED_BYTES {
        *out = [0; NLIMBS];
        return 0;
    }
    *out = [0; NLIMBS];
    // Build the integer by shifting in one byte at a time, MSB first.
    for i in (0..ENCODED_BYTES).rev() {
        modshl(out, 8);
        out[0] = out[0].wrapping_add(bytes[i] as u64);
    }
    // res is all-ones if the value was in `[0, p)`, all-zeros otherwise.
    let res_u64 = 0u64.wrapping_sub(modfsb(out) as u64);
    let res_u32 = res_u64 as u32;
    let out_copy = *out;
    nres(out, &out_copy);
    for limb in out.iter_mut() {
        *limb &= res_u64;
    }
    res_u32
}

/// Partial reduction of a 6-limb 64-bit-saturated accumulator: split
/// off the top byte of limb 5 and fold it back into the low limbs
/// using the identity `65 * 2^376 = 1 (mod p)`. Used by
/// [`fp_decode_reduce`] to absorb each 48-byte chunk of a long input.
fn partial_reduce_6(out: &mut [u64; 6], src: &[u64; 6]) {
    let h = src[5] >> 56;
    let l = src[5] & 0x00FF_FFFF_FFFF_FFFF;
    // Add floor(h/65) + (h mod 65) * 2^376 to the low part.
    let quo = (h.wrapping_mul(0xFC1)) >> 18;
    let rem = h.wrapping_sub(65u64.wrapping_mul(quo));
    let (r0, c0) = src[0].overflowing_add(quo);
    let (r1, c1) = src[1].overflowing_add(c0 as u64);
    let (r2, c2) = src[2].overflowing_add(c1 as u64);
    let (r3, c3) = src[3].overflowing_add(c2 as u64);
    let (r4, c4) = src[4].overflowing_add(c3 as u64);
    let r5 = l.wrapping_add(rem << 56).wrapping_add(c4 as u64);
    out[0] = r0;
    out[1] = r1;
    out[2] = r2;
    out[3] = r3;
    out[4] = r4;
    out[5] = r5;
}

/// Parse a little-endian byte string of arbitrary length, reducing it
/// modulo `p`. Always succeeds. Used to map hash output uniformly into
/// Fp.
#[inline]
pub(crate) fn fp_decode_reduce(out: &mut [u64; NLIMBS], bytes: &[u8]) {
    *out = [0; NLIMBS];
    if bytes.is_empty() {
        return;
    }

    let mut len = bytes.len();
    let rem = len % ENCODED_BYTES;
    if rem != 0 {
        // Decode a partial trailing block (zero-padded to 48 bytes;
        // the value is already < p so no reduction is needed).
        let k = len - rem;
        let mut tmp = [0u8; ENCODED_BYTES];
        tmp[..(len - k)].copy_from_slice(&bytes[k..]);
        let _ = fp_decode(out, &tmp);
        len = k;
    }

    while len > 0 {
        // Shift the accumulator left by 2^384 in the Montgomery sense:
        // multiplying by `R2` (which represents 2^384) makes room for
        // the next 48-byte chunk in the high bits.
        let out_copy = *out;
        modmul(out, &out_copy, &R2);
        len -= ENCODED_BYTES;
        let mut t = [0u64; 6];
        for (j, t_j) in t.iter_mut().enumerate() {
            let off = len + j * 8;
            *t_j = u64::from_le_bytes(
                bytes[off..off + 8]
                    .try_into()
                    .expect("invariant: slice length is exactly 8"),
            );
        }
        let mut reduced = [0u64; 6];
        partial_reduce_6(&mut reduced, &t);
        let mut tmp_bytes = [0u8; ENCODED_BYTES];
        for (j, r) in reduced.iter().enumerate() {
            tmp_bytes[j * 8..(j + 1) * 8].copy_from_slice(&r.to_le_bytes());
        }
        let mut a = [0u64; NLIMBS];
        let _ = fp_decode(&mut a, &tmp_bytes);
        let out_copy = *out;
        modadd(out, &out_copy, &a);
    }
}

/// Borrow a `&Array<u64, U7>` as `&[u64; 7]` for the limb-level helpers.
#[inline]
fn as_arr(a: &Array<u64, <Level3 as crate::params::SecurityLevel>::FpLimbs>) -> &[u64; NLIMBS] {
    // SAFETY: Array<u64, U7> has exactly NLIMBS=7 elements, matching [u64; 7]
    <&[u64; NLIMBS]>::try_from(&a[..]).expect("invariant: Level3 FpLimbs == U7")
}

/// Mutable borrow of `&mut Array<u64, U7>` as `&mut [u64; 7]`.
#[inline]
fn as_arr_mut(
    a: &mut Array<u64, <Level3 as crate::params::SecurityLevel>::FpLimbs>,
) -> &mut [u64; NLIMBS] {
    // SAFETY: Array<u64, U7> has exactly NLIMBS=7 elements, matching [u64; 7]
    <&mut [u64; NLIMBS]>::try_from(&mut a[..]).expect("invariant: Level3 FpLimbs == U7")
}

impl FpBackend for Level3 {
    #[inline]
    fn set_zero(out: &mut Array<u64, Self::FpLimbs>) {
        *as_arr_mut(out) = ZERO;
    }

    #[inline]
    fn set_one(out: &mut Array<u64, Self::FpLimbs>) {
        *as_arr_mut(out) = ONE;
    }

    #[inline]
    fn set_small(out: &mut Array<u64, Self::FpLimbs>, val: u64) {
        modint(as_arr_mut(out), val);
    }

    #[inline]
    fn is_equal(a: &Array<u64, Self::FpLimbs>, b: &Array<u64, Self::FpLimbs>) -> Choice {
        Choice::from(modcmp(as_arr(a), as_arr(b)) as u8)
    }

    #[inline]
    fn is_zero(a: &Array<u64, Self::FpLimbs>) -> Choice {
        Choice::from(modis0(as_arr(a)) as u8)
    }

    #[inline]
    fn copy(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>) {
        modcpy(as_arr_mut(out), as_arr(a));
    }

    #[inline]
    fn add(
        out: &mut Array<u64, Self::FpLimbs>,
        a: &Array<u64, Self::FpLimbs>,
        b: &Array<u64, Self::FpLimbs>,
    ) {
        modadd(as_arr_mut(out), as_arr(a), as_arr(b));
    }

    #[inline]
    fn sub(
        out: &mut Array<u64, Self::FpLimbs>,
        a: &Array<u64, Self::FpLimbs>,
        b: &Array<u64, Self::FpLimbs>,
    ) {
        modsub(as_arr_mut(out), as_arr(a), as_arr(b));
    }

    #[inline]
    fn neg(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>) {
        modneg(as_arr_mut(out), as_arr(a));
    }

    #[inline]
    fn mul(
        out: &mut Array<u64, Self::FpLimbs>,
        a: &Array<u64, Self::FpLimbs>,
        b: &Array<u64, Self::FpLimbs>,
    ) {
        modmul(as_arr_mut(out), as_arr(a), as_arr(b));
    }

    #[inline]
    fn sqr(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>) {
        modsqr(as_arr_mut(out), as_arr(a));
    }

    #[inline]
    fn inv(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>) {
        let a_copy = *as_arr(a);
        modinv(as_arr_mut(out), &a_copy, None);
    }

    #[inline]
    fn sqrt(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>) {
        let a_copy = *as_arr(a);
        modsqrt(as_arr_mut(out), &a_copy, None);
    }

    #[inline]
    fn is_square(a: &Array<u64, Self::FpLimbs>) -> Choice {
        Choice::from(modqr(as_arr(a), None) as u8)
    }

    #[inline]
    fn half(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>) {
        modmul(as_arr_mut(out), &TWO_INV, as_arr(a));
    }

    #[inline]
    fn div3(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>) {
        modmul(as_arr_mut(out), &THREE_INV, as_arr(a));
    }

    #[inline]
    fn exp3div4(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>) {
        modpro(as_arr_mut(out), as_arr(a));
    }

    #[inline]
    fn mul_small(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>, val: u32) {
        modmli(as_arr_mut(out), as_arr(a), val as u64);
    }

    #[inline]
    fn encode(out: &mut [u8], a: &Array<u64, Self::FpLimbs>) {
        fp_encode(out, as_arr(a));
    }

    #[inline]
    fn decode(out: &mut Array<u64, Self::FpLimbs>, bytes: &[u8]) -> Choice {
        Choice::from((fp_decode(as_arr_mut(out), bytes) & 1) as u8)
    }

    #[inline]
    fn decode_reduce(out: &mut Array<u64, Self::FpLimbs>, bytes: &[u8]) {
        fp_decode_reduce(as_arr_mut(out), bytes);
    }

    #[inline]
    fn cswap(a: &mut Array<u64, Self::FpLimbs>, b: &mut Array<u64, Self::FpLimbs>, ctl: Choice) {
        modcsw(ctl.unwrap_u8() as u64, as_arr_mut(a), as_arr_mut(b));
    }

    #[inline]
    fn select(
        out: &mut Array<u64, Self::FpLimbs>,
        a0: &Array<u64, Self::FpLimbs>,
        a1: &Array<u64, Self::FpLimbs>,
        ctl: Choice,
    ) {
        let cw = 0u64.wrapping_sub(ctl.unwrap_u8() as u64);
        let a0r = as_arr(a0);
        let a1r = as_arr(a1);
        let out_r = as_arr_mut(out);
        for i in 0..NLIMBS {
            out_r[i] = a0r[i] ^ (cw & (a0r[i] ^ a1r[i]));
        }
    }
}
