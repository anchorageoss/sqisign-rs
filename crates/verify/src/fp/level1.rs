//!
//! The prime is `p = 5 * 2^248 - 1`. Field elements are 5-limb arrays
//! `[u64; 5]` storing values in unsaturated radix `2^51`: each limb
//! holds a roughly-51-bit quantity, leaving 13 carry bits at the top.
//! Arithmetic operates in an internal Montgomery form; canonical
//! little-endian byte encoding (via `fp_encode` / `fp_decode`)
//! returns the standard integer representation in `[0, p)`.
//!
//! The free functions below implement the limb-level primitives;
//! the [`FpBackend`] trait impl for [`Level1`] at the bottom of the
//! file wires them into the generic [`super::Fp`] API.

use super::FpBackend;
use crate::params::Level1;
use hybrid_array::Array;
use subtle::Choice;

/// Number of 64-bit limbs in an `Fp` element.
pub const NLIMBS: usize = 5;

/// Bit width of each unsaturated limb.
pub const RADIX: u32 = 51;

/// `2^RADIX - 1`: keeps the low 51 bits of a limb.
pub const MASK: u64 = (1u64 << RADIX) - 1;

/// The contribution of `p + 1 = 5 * 2^248` to limb 4 in radix-`2^51`:
/// `5 * 2^248 / 2^(51*4) = 5 * 2^44`. Used as the Montgomery folding
/// constant inside `modmul` and `modsqr`.
pub const P4: u64 = 0x5000_0000_0000;

/// `2 * P4`, used inside `modadd`, `modsub`, and `modneg` to
/// keep partial sums in `[-p, p)` before the final reduction step.
pub const TWO_P4: u64 = 0xa000_0000_0000;

/// `0` as 5 radix-`2^51` limbs.
pub const ZERO: [u64; NLIMBS] = [0; NLIMBS];

/// Internal Montgomery form of `1 mod p`.
pub const ONE: [u64; NLIMBS] = [
    0x0000_0000_0000_0019,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_3000_0000_0000,
];

/// Internal Montgomery form of 2⁻¹ mod p.
pub const TWO_INV: [u64; NLIMBS] = [
    0x0000_0000_0000_000c,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_4000_0000_0000,
];

/// Internal Montgomery form of 3⁻¹ mod p.
pub const THREE_INV: [u64; NLIMBS] = [
    0x0005_5555_5555_555d,
    0x0002_aaaa_aaaa_aaaa,
    0x0005_5555_5555_5555,
    0x0002_aaaa_aaaa_aaaa,
    0x0000_4555_5555_5555,
];

/// Internal Montgomery form of `2^256 mod p`. Used by
/// `fp_decode_reduce` when folding 32-byte input blocks.
pub const R2: [u64; NLIMBS] = [
    0x0001_9999_9999_9eb8,
    0x0003_3333_3333_3333,
    0x0006_6666_6666_6666,
    0x0004_cccc_cccc_cccc,
    0x0000_1999_9999_9999,
];

/// Conversion constant for `nres`: multiplying by `R2_NRES` and
/// dividing by Montgomery `R` takes a canonical integer into internal
/// Montgomery form (`R2_NRES = R^2 mod p` in the internal Montgomery
/// scheme).
pub const R2_NRES: [u64; NLIMBS] = [
    0x0004_cccc_cccc_cf5c,
    0x0001_9999_9999_9999,
    0x0003_3333_3333_3333,
    0x0006_6666_6666_6666,
    0x0000_0ccc_cccc_cccc,
];

/// Propagate carries between limbs, leaving each of limbs 0..3 in
/// `[0, 2^51)` and limb 4 holding the accumulated high bits. Returns
/// `0xFFFF_FFFF_FFFF_FFFF` if the value (interpreted as two's-complement
/// via the top bit of limb 4) is negative, else `0`.
#[inline]
pub(crate) fn prop(n: &mut [u64; NLIMBS]) -> u64 {
    let mask = MASK;
    let mut carry: i64 = n[0] as i64;
    carry >>= RADIX;
    n[0] &= mask;
    for limb in n.iter_mut().take(4).skip(1) {
        carry = carry.wrapping_add(*limb as i64);
        *limb = (carry as u64) & mask;
        carry >>= RADIX;
    }
    n[4] = n[4].wrapping_add(carry as u64);
    // Sign mask: 0 if non-negative, all-ones if negative.
    let sign = n[4] >> 63;
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
    n[4] = n[4].wrapping_add(0x5000_0000_0000 & carry);
    let _ = prop(n);
    (carry & 1) as u32
}

/// Final subtract of `p`: tries `n -= p` and adds `p` back if the
/// result went negative. Returns `1` if the original value was `< p`
/// (and was preserved), `0` if it was `>= p` (and was reduced).
#[inline]
pub(crate) fn modfsb(n: &mut [u64; NLIMBS]) -> u32 {
    n[0] = n[0].wrapping_add(1);
    n[4] = n[4].wrapping_sub(0x5000_0000_0000);
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
    n[0] = n[0].wrapping_add(2);
    n[4] = n[4].wrapping_sub(TWO_P4);
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2 & carry);
    n[4] = n[4].wrapping_add(TWO_P4 & carry);
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
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2 & carry);
    n[4] = n[4].wrapping_add(TWO_P4 & carry);
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
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2 & carry);
    n[4] = n[4].wrapping_add(TWO_P4 & carry);
    let _ = prop(n);
}

/// `c <- a * b mod 2p`. Schoolbook 5x5 multiplication in radix `2^51`
/// with interleaved Montgomery folding: each newly-emitted low limb
/// `v_i` is added back at limb position `i + 4` as `v_i * P4`, which
/// represents `v_i * 5 * 2^248 = v_i * (p + 1) = v_i (mod p)` and
/// effectively divides the final result by Montgomery `R = 2^255`.
#[inline]
pub(crate) fn modmul(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS], b: &[u64; NLIMBS]) {
    let mask = MASK;
    let p4: u128 = P4 as u128;

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
        .wrapping_add((a[4] as u128) * (b[0] as u128))
        .wrapping_add((v0 as u128) * p4);
    let v4 = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[1] as u128) * (b[4] as u128))
        .wrapping_add((a[2] as u128) * (b[3] as u128))
        .wrapping_add((a[3] as u128) * (b[2] as u128))
        .wrapping_add((a[4] as u128) * (b[1] as u128))
        .wrapping_add((v1 as u128) * p4);
    c[0] = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[2] as u128) * (b[4] as u128))
        .wrapping_add((a[3] as u128) * (b[3] as u128))
        .wrapping_add((a[4] as u128) * (b[2] as u128))
        .wrapping_add((v2 as u128) * p4);
    c[1] = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[3] as u128) * (b[4] as u128))
        .wrapping_add((a[4] as u128) * (b[3] as u128))
        .wrapping_add((v3 as u128) * p4);
    c[2] = (t as u64) & mask;
    t >>= RADIX;

    t = t
        .wrapping_add((a[4] as u128) * (b[4] as u128))
        .wrapping_add((v4 as u128) * p4);
    c[3] = (t as u64) & mask;
    t >>= RADIX;

    c[4] = t as u64;
}

/// `c <- a * a mod 2p` (specialized squaring). Each cross-term
/// `a_i * a_j` (`i != j`) is computed once and doubled, saving roughly
/// half the partial-product work compared to a general multiply.
#[inline]
pub(crate) fn modsqr(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS]) {
    let mask = MASK;
    let p4: u128 = P4 as u128;

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
    t = t.wrapping_add((v0 as u128) * p4);
    let v4 = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[1] as u128) * (a[4] as u128);
    tot = tot.wrapping_add((a[2] as u128) * (a[3] as u128));
    tot = tot.wrapping_mul(2);
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v1 as u128) * p4);
    c[0] = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[2] as u128) * (a[4] as u128);
    tot = tot.wrapping_mul(2);
    tot = tot.wrapping_add((a[3] as u128) * (a[3] as u128));
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v2 as u128) * p4);
    c[1] = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[3] as u128) * (a[4] as u128);
    tot = tot.wrapping_mul(2);
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v3 as u128) * p4);
    c[2] = (t as u64) & mask;
    t >>= RADIX;

    tot = (a[4] as u128) * (a[4] as u128);
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v4 as u128) * p4);
    c[3] = (t as u64) & mask;
    t >>= RADIX;

    c[4] = t as u64;
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
    let mut t1 = [0u64; NLIMBS];
    let mut t2 = [0u64; NLIMBS];
    let mut t3 = [0u64; NLIMBS];
    let mut t4 = [0u64; NLIMBS];

    modcpy(&mut x, w);
    modsqr(z, &x);
    modmul(&mut t0, &x, z);
    modsqr(z, &t0);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &x, z);
        *z = tmp;
    }
    modsqr(&mut t1, z);
    modsqr(&mut t3, &t1);
    modsqr(&mut t2, &t3);
    modcpy(&mut t4, &t2);
    modnsqr(&mut t4, 3);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t2, &t4);
        t2 = tmp;
    }
    modcpy(&mut t4, &t2);
    modnsqr(&mut t4, 6);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t2, &t4);
        t2 = tmp;
    }
    modcpy(&mut t4, &t2);
    modnsqr(&mut t4, 2);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t3, &t4);
        t3 = tmp;
    }
    modnsqr(&mut t3, 13);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t2, &t3);
        t2 = tmp;
    }
    modcpy(&mut t3, &t2);
    modnsqr(&mut t3, 27);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t2, &t3);
        t2 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, z, &t2);
        *z = tmp;
    }
    modcpy(&mut t2, z);
    modnsqr(&mut t2, 4);
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
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t1, &t0);
        t1 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t0, &t1);
        t0 = tmp;
    }
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t1, &t0);
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
    modnsqr(&mut t1, 63);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t0, &t1);
        t1 = tmp;
    }
    modnsqr(&mut t1, 64);
    {
        let mut tmp = [0u64; NLIMBS];
        modmul(&mut tmp, &t0, &t1);
        t0 = tmp;
    }
    modnsqr(&mut t0, 57);
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

/// Shift `a` left by `n < 51` bits (per-limb, with carry into the
/// next limb).
#[inline]
pub(crate) fn modshl(a: &mut [u64; NLIMBS], n: u32) {
    a[4] = (a[4] << n) | (a[3] >> (RADIX - n));
    for i in (1..=3).rev() {
        a[i] = ((a[i] << n) & MASK) | (a[i - 1] >> (RADIX - n));
    }
    a[0] = (a[0] << n) & MASK;
}

/// Shift `a` right by `n < 51` bits, returning the low `n`
/// shifted-out bits.
#[inline]
pub(crate) fn modshr(a: &mut [u64; NLIMBS], n: u32) -> u64 {
    let r = a[0] & ((1u64 << n) - 1);
    for i in 0..4 {
        a[i] = (a[i] >> n) | ((a[i + 1] << (RADIX - n)) & MASK);
    }
    a[4] >>= n;
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

/// Bytes per encoded Fp element for Level 1.
const ENCODED_BYTES: usize = 32;

/// Serialize `a` (in internal Montgomery form) to its canonical
/// 32-byte little-endian representation: convert back to canonical
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

/// Parse 32 canonical little-endian bytes into an Fp element.
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

/// Partial reduction of a 4-limb 64-bit-saturated accumulator: split
/// off the top byte of limb 3 and fold it back into the low limbs
/// using the identity `5 * 2^248 = 1 (mod p)`. Used by
/// [`fp_decode_reduce`] to absorb each 32-byte chunk of a long input.
fn partial_reduce_4(out: &mut [u64; 4], src: &[u64; 4]) {
    let h = src[3] >> 56;
    let l = src[3] & 0x00FF_FFFF_FFFF_FFFF;
    // Add floor(h/5) + (h mod 5) * 2^248 to the low part.
    let quo = (h.wrapping_mul(0xCD)) >> 10;
    let rem = h.wrapping_sub(5u64.wrapping_mul(quo));
    let (r0, c0) = src[0].overflowing_add(quo);
    let (r1, c1) = src[1].overflowing_add(c0 as u64);
    let (r2, c2) = src[2].overflowing_add(c1 as u64);
    let r3 = l.wrapping_add(rem << 56).wrapping_add(c2 as u64);
    out[0] = r0;
    out[1] = r1;
    out[2] = r2;
    out[3] = r3;
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
        // Decode a partial trailing block (zero-padded to 32 bytes;
        // the value is already < p so no reduction is needed).
        let k = len - rem;
        let mut tmp = [0u8; ENCODED_BYTES];
        tmp[..(len - k)].copy_from_slice(&bytes[k..]);
        let _ = fp_decode(out, &tmp);
        len = k;
    }

    while len > 0 {
        // Shift the accumulator left by 2^256 in the Montgomery sense:
        // multiplying by `R2` (which represents 2^256) makes room for
        // the next 32-byte chunk in the high bits.
        let out_copy = *out;
        modmul(out, &out_copy, &R2);
        len -= ENCODED_BYTES;
        let mut t = [0u64; 4];
        for (j, t_j) in t.iter_mut().enumerate() {
            let off = len + j * 8;
            // SAFETY: slice length is exactly 8, matching the [u8; 8] target
            *t_j = u64::from_le_bytes(
                bytes[off..off + 8]
                    .try_into()
                    .expect("invariant: slice length is exactly 8"),
            );
        }
        let mut reduced = [0u64; 4];
        partial_reduce_4(&mut reduced, &t);
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

/// Borrow a `&Array<u64, U5>` as `&[u64; 5]` for the limb-level helpers.
#[inline]
fn as_arr(a: &Array<u64, <Level1 as crate::params::SecurityLevel>::FpLimbs>) -> &[u64; NLIMBS] {
    <&[u64; NLIMBS]>::try_from(&a[..]).expect("invariant: Level1 FpLimbs == U5")
}

/// Mutable borrow of `&mut Array<u64, U5>` as `&mut [u64; 5]`.
#[inline]
fn as_arr_mut(
    a: &mut Array<u64, <Level1 as crate::params::SecurityLevel>::FpLimbs>,
) -> &mut [u64; NLIMBS] {
    <&mut [u64; NLIMBS]>::try_from(&mut a[..]).expect("invariant: Level1 FpLimbs == U5")
}

impl FpBackend for Level1 {
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
