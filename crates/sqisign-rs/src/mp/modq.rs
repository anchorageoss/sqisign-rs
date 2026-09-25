//! Montgomery arithmetic modulo an odd integer in an unsaturated radix-61
//! representation (the reference's `modqx`), and the modular functions
//! built on it: exponentiation, square roots of `-1`, square roots modulo
//! a prime.

use super::rand::{DefaultDomain, Rng};
use super::{Ibz, IBZ_MAX_LIMBS};

/// Bits per unsaturated limb.
pub(crate) const RADIX: u32 = 61;
/// Limb mask.
pub(crate) const LIMBMASK: u64 = (1u64 << RADIX) - 1;
/// Largest number of radix-61 limbs any container needs.
pub(crate) const MAX_MODQ_LIMBS: usize = (IBZ_MAX_LIMBS * 64 - 1 + 60) / 61;

const WINDOW_SIZE: u32 = 3;
const TABLE_SIZE: usize = 1 << WINDOW_SIZE;

/// Storage for one residue.
pub(crate) type Spint = [u64; MAX_MODQ_LIMBS];

#[inline]
pub(crate) fn spint_zero() -> Spint {
    [0u64; MAX_MODQ_LIMBS]
}

/// Propagate carries; returns all-ones when the result is negative.
fn prop(n: &mut [u64], numwords: usize) -> u64 {
    let mut carry = (n[0] as i64) >> RADIX;
    if numwords == 1 {
        return 0u64.wrapping_sub((n[0] >> 1) >> 62);
    }
    n[0] &= LIMBMASK;
    for i in 1..numwords - 1 {
        carry += n[i] as i64;
        n[i] = (carry as u64) & LIMBMASK;
        carry >>= RADIX;
    }
    n[numwords - 1] = n[numwords - 1].wrapping_add(carry as u64);
    0u64.wrapping_sub((n[numwords - 1] >> 1) >> 62)
}

/// Propagate carries and add `p` if negative.
fn flatten(n: &mut [u64], p: &[u64], numwords: usize) -> u64 {
    let carry = prop(n, numwords);
    for i in 0..numwords {
        n[i] = n[i].wrapping_add(p[i] & carry);
    }
    let _ = prop(n, numwords);
    carry & 1
}

/// Montgomery final subtraction.
fn modfsb(n: &mut [u64], p: &[u64], numwords: usize) -> u64 {
    for i in 0..numwords {
        n[i] = n[i].wrapping_sub(p[i]);
    }
    flatten(n, p, numwords)
}

/// `n = a + b mod 2p`.
pub(crate) fn modadd(a: &[u64], b: &[u64], n: &mut [u64], two_p: &[u64], numwords: usize) {
    for i in 0..numwords {
        n[i] = a[i].wrapping_add(b[i]).wrapping_sub(two_p[i]);
    }
    let carry = prop(n, numwords);
    for i in 0..numwords {
        n[i] = n[i].wrapping_add(two_p[i] & carry);
    }
    let _ = prop(n, numwords);
}

/// `n = a - b mod 2p`.
pub(crate) fn modsub(a: &[u64], b: &[u64], n: &mut [u64], two_p: &[u64], numwords: usize) {
    for i in 0..numwords {
        n[i] = a[i].wrapping_sub(b[i]);
    }
    let carry = prop(n, numwords);
    for i in 0..numwords {
        n[i] = n[i].wrapping_add(two_p[i] & carry);
    }
    let _ = prop(n, numwords);
}

/// Montgomery product `c = a * b`, output below `2p`.
pub(crate) fn modmul(a: &[u64], b: &[u64], c: &mut [u64], p: &[u64], ndash: u64, numwords: usize) {
    let mut t: u128 = 0;
    let mut v = spint_zero();
    for i in 0..numwords {
        for j in 0..i {
            t += (a[j] as u128) * (b[i - j] as u128) + (v[j] as u128) * (p[i - j] as u128);
        }
        t += (a[i] as u128) * (b[0] as u128);
        v[i] = ((t as u64).wrapping_mul(ndash)) & LIMBMASK;
        t += (v[i] as u128) * (p[0] as u128);
        t >>= RADIX;
    }
    for i in numwords..2 * numwords - 1 {
        let jlo = i - numwords + 1;
        for j in jlo..numwords {
            t += (a[j] as u128) * (b[i - j] as u128) + (v[j] as u128) * (p[i - j] as u128);
        }
        c[i - numwords] = (t as u64) & LIMBMASK;
        t >>= RADIX;
    }
    c[numwords - 1] = t as u64;
}

/// Montgomery square.
pub(crate) fn modsqr(a: &[u64], c: &mut [u64], p: &[u64], ndash: u64, numwords: usize) {
    let mut t: u128 = 0;
    let mut v = spint_zero();
    for i in 0..numwords {
        let mut tot: u128 = 0;
        let mut j = 0;
        while 2 * j < i {
            tot += (a[j] as u128) * (a[i - j] as u128);
            j += 1;
        }
        tot *= 2;
        if i & 1 == 0 {
            let mid = i / 2;
            tot += (a[mid] as u128) * (a[mid] as u128);
        }
        t += tot;
        for j in 0..i {
            t += (v[j] as u128) * (p[i - j] as u128);
        }
        v[i] = ((t as u64).wrapping_mul(ndash)) & LIMBMASK;
        t += (v[i] as u128) * (p[0] as u128);
        t >>= RADIX;
    }
    for i in numwords..2 * numwords - 1 {
        let jlo = i - numwords + 1;
        let mut tot: u128 = 0;
        let mut j = jlo;
        while 2 * j < i {
            tot += (a[j] as u128) * (a[i - j] as u128);
            j += 1;
        }
        tot *= 2;
        if i & 1 == 0 {
            let mid = i / 2;
            tot += (a[mid] as u128) * (a[mid] as u128);
        }
        t += tot;
        for j in jlo..numwords {
            t += (v[j] as u128) * (p[i - j] as u128);
        }
        c[i - numwords] = (t as u64) & LIMBMASK;
        t >>= RADIX;
    }
    c[numwords - 1] = t as u64;
}

/// Into Montgomery form knowing only `2p`.
pub(crate) fn nresx(m: &[u64], n: &mut [u64], two_p: &[u64], numwords: usize) {
    let mut tmp = spint_zero();
    modadd(m, m, &mut tmp, two_p, numwords);
    n[..numwords].copy_from_slice(&tmp[..numwords]);
    for _ in 0..RADIX as usize * numwords - 1 {
        let cur = *n_as_spint(n);
        modadd(&cur, &cur, n, two_p, numwords);
    }
}

#[inline]
fn n_as_spint(n: &[u64]) -> &Spint {
    // slices passed here always come from a Spint
    n.try_into().expect("residue buffer")
}

/// Out of Montgomery form, fully reduced.
pub(crate) fn redc(n: &[u64], m: &mut [u64], p: &[u64], ndash: u64, numwords: usize) {
    let mut c = spint_zero();
    c[0] = 1;
    modmul(n, &c, m, p, ndash, numwords);
    let _ = modfsb(m, p, numwords);
}

/// Whether the residue is zero.
pub(crate) fn modis0(a: &[u64], p: &[u64], ndash: u64, numwords: usize) -> bool {
    let mut c = spint_zero();
    redc(a, &mut c, p, ndash, numwords);
    let mut d = 0u64;
    for &x in &c[..numwords] {
        d |= x;
    }
    (1 & (d.wrapping_sub(1) >> RADIX)) == 1
}

/// `R mod 2p` by Knuth division, for `two_p` with `active` significant
/// limbs.
fn modone_div(a: &mut [u64], two_p: &[u64], numwords: usize, active: usize) {
    for x in a[..numwords].iter_mut() {
        *x = 0;
    }
    let mask = LIMBMASK;
    let base: u128 = 1u128 << RADIX;
    let mut v = [0u64; MAX_MODQ_LIMBS + 1];
    let mut u = [0u64; MAX_MODQ_LIMBS + 2];
    let top = two_p[active - 1] & mask;
    debug_assert!(top != 0);
    let shift = top.leading_zeros() - (64 - RADIX);
    let mut carry = 0u64;
    if shift == 0 {
        for i in 0..active {
            v[i] = two_p[i] & mask;
        }
    } else {
        for i in 0..active {
            let z = (((two_p[i] & mask) as u128) << shift) | carry as u128;
            v[i] = (z as u64) & mask;
            carry = (z >> RADIX) as u64;
        }
        debug_assert_eq!(carry, 0);
    }
    u[numwords] = 1u64 << shift;
    let jmax = (numwords + 1) - active;
    for j in (0..=jmax).rev() {
        let numerator = ((u[j + active] as u128) << RADIX) | u[j + active - 1] as u128;
        let vtop = v[active - 1] as u128;
        let mut qhat = numerator / vtop;
        let mut rhat = numerator - qhat * vtop;
        let cap = (qhat >= base) as u128;
        qhat -= cap;
        rhat += cap * vtop;
        for _ in 0..2 {
            let lhs = qhat * v[active - 2] as u128;
            let rhs = (rhat << RADIX) + u[j + active - 2] as u128;
            let decrement = ((rhat < base) as u128) & ((lhs > rhs) as u128);
            qhat -= decrement;
            rhat += decrement * vtop;
        }
        let mut borrow: u128 = 0;
        for i in 0..active {
            let product = qhat * v[i] as u128 + borrow;
            let lo = (product as u64) & mask;
            let hi = (product >> RADIX) as u64;
            let ui = u[j + i] & mask;
            let under = (ui < lo) as u64;
            u[j + i] = ui.wrapping_sub(lo) & mask;
            borrow = hi as u128 + under as u128;
        }
        let top_u = u[j + active] & mask;
        let borrow_digit = borrow as u64;
        let negative = (top_u < borrow_digit) as u64;
        u[j + active] = top_u.wrapping_sub(borrow_digit) & mask;
        let add_mask = 0u64.wrapping_sub(negative);
        let mut addcarry: u128 = 0;
        for i in 0..active {
            let z = (u[j + i] & mask) as u128 + (v[i] & add_mask) as u128 + addcarry;
            u[j + i] = (z as u64) & mask;
            addcarry = z >> RADIX;
        }
        u[j + active] = (u[j + active].wrapping_add(addcarry as u64)) & mask;
    }
    if shift == 0 {
        for i in 0..active {
            a[i] = u[i] & mask;
        }
    } else {
        for i in 0..active {
            let lo = (u[i] & mask) >> shift;
            let hi = if i + 1 < active {
                ((u[i + 1] & mask) << (RADIX - shift)) & mask
            } else {
                0
            };
            a[i] = (lo | hi) & mask;
        }
    }
}

/// Montgomery representation of 1.
pub(crate) fn modone(a: &mut [u64], two_p: &[u64], numwords: usize) {
    if numwords > 1 && (two_p[numwords - 1] & LIMBMASK) != 0 {
        modone_div(a, two_p, numwords, numwords);
        return;
    }
    if numwords > 2 && (two_p[numwords - 2] & LIMBMASK) != 0 {
        modone_div(a, two_p, numwords, numwords - 1);
        return;
    }
    a[0] = 1;
    for x in a[1..numwords].iter_mut() {
        *x = 0;
    }
    for _ in 0..numwords * RADIX as usize {
        let cur = *n_as_spint(a);
        modadd(&cur, &cur, a, two_p, numwords);
    }
}

/// Conditional move `f = g` when `b == 1`, with the reference's masking
/// trick.
fn modcmv(b: u64, g: &[u64], f: &mut [u64], numwords: usize) {
    let w = super::ct::barrier(0x3cc3_c33c_5aa5_a55a);
    let c0 = (!b) & w.wrapping_add(1);
    let c1 = b.wrapping_add(w);
    for i in 0..numwords {
        let s = g[i];
        let t = f[i];
        let aux = c0.wrapping_mul(t).wrapping_add(c1.wrapping_mul(s));
        f[i] = aux.wrapping_sub(w.wrapping_mul(t.wrapping_add(s)));
    }
}

/// Conditional swap when `b == 1`.
fn modcsw(b: u64, g: &mut [u64], f: &mut [u64], numwords: usize) {
    let w = super::ct::barrier(0x3cc3_c33c_5aa5_a55a);
    let c0 = (!b) & w.wrapping_add(1);
    let c1 = b.wrapping_add(w);
    for i in 0..numwords {
        let s = g[i];
        let t = f[i];
        let v = w.wrapping_mul(t.wrapping_add(s));
        let aux = c0.wrapping_mul(t).wrapping_add(c1.wrapping_mul(s));
        f[i] = aux.wrapping_sub(v);
        let aux = c0.wrapping_mul(s).wrapping_add(c1.wrapping_mul(t));
        g[i] = aux.wrapping_sub(v);
    }
}

/// Shift left by less than a limb.
pub(crate) fn modshl(n: u32, a: &mut [u64], numwords: usize) {
    if numwords == 1 {
        a[0] <<= n;
        return;
    }
    a[numwords - 1] = (a[numwords - 1] << n).wrapping_add(a[numwords - 2] >> (RADIX - n));
    for i in (1..numwords - 1).rev() {
        a[i] = ((a[i] << n) & LIMBMASK).wrapping_add(a[i - 1] >> (RADIX - n));
    }
    a[0] = (a[0] << n) & LIMBMASK;
}

/// Shift right by less than a limb; returns the bits shifted out.
#[allow(dead_code)]
pub(crate) fn modshr(n: u32, a: &mut [u64], numwords: usize) -> u64 {
    let r = a[0] & ((1u64 << n) - 1);
    for i in 0..numwords - 1 {
        a[i] = (a[i] >> n).wrapping_add((a[i + 1] << (RADIX - n)) & LIMBMASK);
    }
    a[numwords - 1] >>= n;
    r
}

/// Equality of residues.
pub(crate) fn modcmp(a: &[u64], b: &[u64], p: &[u64], ndash: u64, numwords: usize) -> bool {
    let mut c = spint_zero();
    let mut d = spint_zero();
    redc(a, &mut c, p, ndash, numwords);
    redc(b, &mut d, p, ndash, numwords);
    let mut eq = 1u64;
    for i in 0..numwords {
        eq &= (((c[i] ^ d[i]).wrapping_sub(1)) >> RADIX) & 1;
    }
    eq == 1
}

/// `(-p)^-1 mod 2^61`.
pub(crate) fn getndash(p0: u64) -> u64 {
    let q = 1u64 << RADIX;
    let mask = q - 1;
    let a = q.wrapping_sub(p0 & mask);
    let mut x = 1u64;
    for _ in 0..6 {
        let t = 2u64.wrapping_sub(a.wrapping_mul(x));
        x = x.wrapping_mul(t);
    }
    x & mask
}

/// Constant-time `x^e` (Montgomery ladder over the whole width).
#[allow(dead_code)]
pub(crate) fn modxpowe(
    p: &[u64],
    one: &[u64],
    ndash: u64,
    x: &[u64],
    e: &[u64],
    s: &mut [u64],
    numwords: usize,
) {
    let mut r0 = spint_zero();
    let mut r1 = spint_zero();
    let mut ee = spint_zero();
    r0[..numwords].copy_from_slice(&one[..numwords]);
    r1[..numwords].copy_from_slice(&x[..numwords]);
    ee[..numwords].copy_from_slice(&e[..numwords]);
    for _ in 0..numwords * RADIX as usize {
        let b = (ee[numwords - 1] >> (RADIX - 1)) & 1;
        modcsw(b, &mut r0, &mut r1, numwords);
        let r0c = r0;
        let r1c = r1;
        modmul(&r0c, &r1c, &mut r1, p, ndash, numwords);
        modsqr(&r0c, &mut r0, p, ndash, numwords);
        modcsw(b, &mut r0, &mut r1, numwords);
        modshl(1, &mut ee, numwords);
    }
    s[..numwords].copy_from_slice(&r0[..numwords]);
}

/// Constant-time table lookup.
fn modtablelookup(table: &[Spint; TABLE_SIZE], index: usize, a: &mut [u64], numwords: usize) {
    for (i, entry) in table.iter().enumerate() {
        let mut b = ((index ^ i) & (TABLE_SIZE - 1)) as u64;
        for k in 0..WINDOW_SIZE {
            b |= b >> k;
        }
        let b = 1 - (b & 1);
        modcmv(b, entry, a, numwords);
    }
}

/// Constant-time windowed `x^e`; `table[0]` holds the Montgomery one.
pub(crate) fn modxpowe_windowed(
    p: &[u64],
    ndash: u64,
    x: &[u64],
    e: &[u64],
    table: &mut [Spint; TABLE_SIZE],
    s: &mut [u64],
    numwords: usize,
) {
    let mut ee = spint_zero();
    let mut tmp = spint_zero();
    let total_bits = RADIX as usize * numwords;
    let num_windows = (total_bits + WINDOW_SIZE as usize - 1) / WINDOW_SIZE as usize;
    let shiftr =
        (num_windows * WINDOW_SIZE as usize) as u32 - RADIX * (numwords as u32 - 1) - WINDOW_SIZE;
    for i in 1..TABLE_SIZE {
        let prev = table[i - 1];
        modmul(&prev, x, &mut table[i], p, ndash, numwords);
    }
    s[..numwords].copy_from_slice(&table[0][..numwords]);
    ee[..numwords].copy_from_slice(&e[..numwords]);
    for _ in 0..num_windows {
        for _ in 0..WINDOW_SIZE {
            let sc = *n_as_spint(s);
            modsqr(&sc, s, p, ndash, numwords);
        }
        let index = ((ee[numwords - 1] >> shiftr) & ((1u64 << WINDOW_SIZE) - 1)) as usize;
        modtablelookup(table, index, &mut tmp, numwords);
        let sc = *n_as_spint(s);
        modmul(&sc, &tmp, s, p, ndash, numwords);
        modshl(WINDOW_SIZE, &mut ee, numwords);
    }
}

/// Variable-time sliding-window `s = x^e` with `x` in Montgomery form.
pub(crate) fn modxpowe_non_ct<const N: usize>(
    s: &mut [u64],
    x: &[u64],
    e: &Ibz<N>,
    n: &[u64],
    two_n: &[u64],
    ndash: u64,
    numwords: usize,
    window: u32,
) {
    debug_assert!((1..=6).contains(&window));
    let table_size = 1usize << (window - 1);
    let mut table = [spint_zero(); 32];
    table[0][..numwords].copy_from_slice(&x[..numwords]);
    if window > 1 {
        let mut x2 = spint_zero();
        modsqr(x, &mut x2, n, ndash, numwords);
        for i in 1..table_size {
            let prev = table[i - 1];
            modmul(&prev, &x2, &mut table[i], n, ndash, numwords);
        }
    }
    let nb = e.bitsize();
    if nb == 0 {
        modone(s, two_n, numwords);
        return;
    }
    let mut i = nb - 1;
    let mut wlen = if i + 1 < window as i32 {
        i + 1
    } else {
        window as i32
    };
    while wlen > 1 && e.bit_at(i - wlen + 1) == 0 {
        wlen -= 1;
    }
    let mut wval = 0usize;
    for k in 0..wlen {
        wval = (wval << 1) | e.bit_at(i - k) as usize;
    }
    s[..numwords].copy_from_slice(&table[(wval - 1) / 2][..numwords]);
    i -= wlen;
    while i >= 0 {
        if e.bit_at(i) == 0 {
            let sc = *n_as_spint(s);
            modsqr(&sc, s, n, ndash, numwords);
            i -= 1;
            continue;
        }
        wlen = if i + 1 < window as i32 {
            i + 1
        } else {
            window as i32
        };
        while wlen > 1 && e.bit_at(i - wlen + 1) == 0 {
            wlen -= 1;
        }
        wval = 0;
        for k in 0..wlen {
            wval = (wval << 1) | e.bit_at(i - k) as usize;
        }
        for _ in 0..wlen {
            let sc = *n_as_spint(s);
            modsqr(&sc, s, n, ndash, numwords);
        }
        let sc = *n_as_spint(s);
        modmul(&sc, &table[(wval - 1) / 2], s, n, ndash, numwords);
        i -= wlen;
    }
}

/// Montgomery context for an odd modulus.
pub(crate) struct ModCtx {
    pub numwords: usize,
    pub n: Spint,
    pub two_n: Spint,
    pub one: Spint,
    pub ndash: u64,
}

impl ModCtx {
    /// Words needed for a modulus of the given bound, with the two extra
    /// bits of headroom the `< 2p` output bound of the product needs.
    pub(crate) fn numwords_for(bitlen: i32) -> usize {
        ((bitlen as u32 + RADIX + 1) / RADIX) as usize
    }

    pub(crate) fn new<const N: usize>(p: &Ibz<N>) -> Self {
        let numwords = Self::numwords_for(p.bitlen);
        debug_assert!(numwords < MAX_MODQ_LIMBS);
        let n = p.to_modq();
        let mut two_n = n;
        modshl(1, &mut two_n, numwords);
        let ndash = getndash(n[0]);
        let mut one = spint_zero();
        modone(&mut one, &two_n, numwords);
        Self {
            numwords,
            n,
            two_n,
            one,
            ndash,
        }
    }
}

impl<const N: usize> Ibz<N> {
    /// Non-negative value into radix-61 limbs.
    #[allow(clippy::wrong_self_convention)]
    pub(crate) fn to_modq(&self) -> Spint {
        debug_assert!(self.is_positive());
        let nlimbs = self.n();
        let mut a = spint_zero();
        for (i, ai) in a.iter_mut().enumerate() {
            let bitpos = (i as u32) * RADIX;
            if bitpos as i32 >= self.bitlen {
                *ai = 0;
                continue;
            }
            let limb = (bitpos / 64) as usize;
            let offset = bitpos % 64;
            let mut v = (self.limbs[limb] as u128) >> offset;
            if offset != 0 && offset + RADIX > 64 && limb + 1 < nlimbs {
                v |= (self.limbs[limb + 1] as u128) << (64 - offset);
            }
            *ai = (v as u64) & LIMBMASK;
        }
        a
    }

    /// Radix-61 limbs into an integer with bound `numwords * 61 + 1`.
    pub(crate) fn from_modq(a: &[u64], numwords: usize) -> Self {
        debug_assert!(numwords < MAX_MODQ_LIMBS);
        debug_assert!((numwords as i32) * (RADIX as i32) < Self::MAX_BITS - 1);
        let mut x = Self::zero();
        for i in (0..numwords).rev() {
            x = x.mul_2exp(RADIX);
            x.limbs[0] |= a[i] & LIMBMASK;
        }
        x.bitlen = (numwords as i32) * RADIX as i32 + 1;
        x
    }

    /// `x^e mod p` for `0 <= x, e < p`, `p` odd; constant time in `x` and
    /// `e` (windowed ladder over the full width).
    pub fn pow_mod(&self, e: &Self, p: &Self) -> Self {
        debug_assert!(p.is_positive() && self.is_positive() && e.is_positive());
        debug_assert!(*self < *p && *e < *p);
        let ctx = ModCtx::new(p);
        let mut xx = self.to_modq();
        let xc = xx;
        nresx(&xc, &mut xx, &ctx.two_n, ctx.numwords);
        let ee = e.to_modq();
        let mut table = [spint_zero(); TABLE_SIZE];
        table[0] = ctx.one;
        let mut out = spint_zero();
        modxpowe_windowed(
            &ctx.n,
            ctx.ndash,
            &xx,
            &ee,
            &mut table,
            &mut out,
            ctx.numwords,
        );
        let mut res = spint_zero();
        redc(&out, &mut res, &ctx.n, ctx.ndash, ctx.numwords);
        let mut r = Self::from_modq(&res, ctx.numwords);
        r.bitlen = p.bitlen;
        r
    }

    /// A square root of `-1` modulo a prime `p = 1 mod 4`, by random
    /// exponentiation. Variable time.
    pub fn sqrt_m1_mod(p: &Self, rng: &mut impl Rng) -> Option<Self> {
        debug_assert!(p.is_positive() && (p.get() & 3) == 1);
        debug_assert!(p.bitlen + 64 < Self::MAX_BITS);
        let window = if p.bitsize() <= 160 { 4 } else { 6 };
        let ctx = ModCtx::new(p);
        let e = p.div_2exp(2);
        let mut xx;
        loop {
            let x = Self::rand_interval(&Self::zero(), p, &mut DefaultDomain(&mut *rng))?;
            xx = x.to_modq();
            let xc = xx;
            modxpowe_non_ct(
                &mut xx,
                &xc,
                &e,
                &ctx.n,
                &ctx.two_n,
                ctx.ndash,
                ctx.numwords,
                window,
            );
            let mut yy = spint_zero();
            modsqr(&xx, &mut yy, &ctx.n, ctx.ndash, ctx.numwords);
            let yc = yy;
            modadd(&yc, &ctx.one, &mut yy, &ctx.two_n, ctx.numwords);
            if modis0(&yy, &ctx.n, ctx.ndash, ctx.numwords) {
                break;
            }
        }
        let mut res = spint_zero();
        redc(&xx, &mut res, &ctx.n, ctx.ndash, ctx.numwords);
        let mut r = Self::from_modq(&res, ctx.numwords);
        r.bitlen = p.bitlen;
        Some(r)
    }

    /// Tonelli-Shanks modulo `p = 1 mod 4`, variable time. `None` for a
    /// non-square.
    fn sqrt_mod_p_1mod4(&self, p: &Self, rng: &mut impl Rng) -> Option<Self> {
        let window = if p.bitsize() <= 160 { 4 } else { 6 };
        let pm1 = p.sub(&Self::one());
        let s = pm1.two_adic();
        let q = pm1.div_2exp(s as u32);
        let qp1o2 = q.add(&Self::one()).div_2exp(1);
        let mut z;
        loop {
            z = Self::rand_interval(&Self::one(), &pm1, &mut DefaultDomain(&mut *rng))?;
            if z.legendre(p) == -1 {
                break;
            }
        }
        let ctx = ModCtx::new(p);
        let nw = ctx.numwords;
        let mut aa = self.to_modq();
        let ac = aa;
        nresx(&ac, &mut aa, &ctx.two_n, nw);
        let mut zz = z.to_modq();
        let zc = zz;
        nresx(&zc, &mut zz, &ctx.two_n, nw);
        let mut cc = spint_zero();
        let mut tt = spint_zero();
        let mut rr = spint_zero();
        modxpowe_non_ct(&mut cc, &zz, &q, &ctx.n, &ctx.two_n, ctx.ndash, nw, window);
        modxpowe_non_ct(&mut tt, &aa, &q, &ctx.n, &ctx.two_n, ctx.ndash, nw, window);
        modxpowe_non_ct(
            &mut rr, &aa, &qp1o2, &ctx.n, &ctx.two_n, ctx.ndash, nw, window,
        );
        let mut m = s;
        while !modcmp(&tt, &ctx.one, &ctx.n, ctx.ndash, nw) {
            let mut tmp = tt;
            let mut i = 0;
            while !modcmp(&tmp, &ctx.one, &ctx.n, ctx.ndash, nw) {
                let tc = tmp;
                modsqr(&tc, &mut tmp, &ctx.n, ctx.ndash, nw);
                i += 1;
                if i >= m {
                    return None;
                }
            }
            let mut bb = cc;
            for _ in 0..m - i - 1 {
                let bc = bb;
                modsqr(&bc, &mut bb, &ctx.n, ctx.ndash, nw);
            }
            let rc = rr;
            modmul(&rc, &bb, &mut rr, &ctx.n, ctx.ndash, nw);
            let mut bsq = spint_zero();
            modsqr(&bb, &mut bsq, &ctx.n, ctx.ndash, nw);
            let tc = tt;
            modmul(&tc, &bsq, &mut tt, &ctx.n, ctx.ndash, nw);
            cc = bsq;
            m = i;
        }
        let mut res = spint_zero();
        redc(&rr, &mut res, &ctx.n, ctx.ndash, nw);
        let mut r = Self::from_modq(&res, nw);
        r.bitlen = p.bitlen;
        Some(r)
    }

    /// Square root modulo `p = 3 mod 4` by exponentiation, variable time.
    fn sqrt_mod_p_3mod4(&self, p: &Self) -> Option<Self> {
        let window = if p.bitsize() <= 160 { 4 } else { 6 };
        let e = p.add(&Self::one()).div_2exp(2);
        let ctx = ModCtx::new(p);
        let nw = ctx.numwords;
        let mut aa = self.to_modq();
        let ac = aa;
        nresx(&ac, &mut aa, &ctx.two_n, nw);
        let mut ss = spint_zero();
        modxpowe_non_ct(&mut ss, &aa, &e, &ctx.n, &ctx.two_n, ctx.ndash, nw, window);
        let mut sq = spint_zero();
        modsqr(&ss, &mut sq, &ctx.n, ctx.ndash, nw);
        let ok = modcmp(&sq, &aa, &ctx.n, ctx.ndash, nw);
        let mut res = spint_zero();
        redc(&ss, &mut res, &ctx.n, ctx.ndash, nw);
        let mut r = Self::from_modq(&res, nw);
        r.bitlen = p.bitlen;
        if ok {
            Some(r)
        } else {
            None
        }
    }

    /// Square root modulo an odd prime, `None` for a non-square (or an even
    /// modulus). Variable time. The random non-residue search for `p = 1
    /// mod 4` draws from `rng`.
    pub fn sqrt_mod_p(&self, p: &Self, rng: &mut impl Rng) -> Option<Self> {
        debug_assert!(p.is_positive() && self.is_positive() && *self < *p);
        debug_assert!(p.bitlen + 64 < Self::MAX_BITS);
        match p.get() & 3 {
            3 => self.sqrt_mod_p_3mod4(p),
            1 => self.sqrt_mod_p_1mod4(p, rng),
            _ => None,
        }
    }

    /// A verified square root of `-1` modulo `n = 1 mod 4` from a small
    /// Jacobi base and one exponentiation; `None` if none is found (the
    /// input is then composite or the bases ran out). Variable time.
    pub fn sqrt_m1_mod_verified(n: &Self) -> Option<Self> {
        debug_assert!(n.is_positive() && (n.limbs[0] & 3) == 1);
        let window = if n.bitsize() <= 160 { 4 } else { 6 };
        let a = super::prime::find_jacobi_minus_one(n)?;
        let ctx = ModCtx::new(n);
        let nw = ctx.numwords;
        // a in Montgomery form by double-and-add on the Montgomery one
        let mut x = ctx.one;
        let mut acc = spint_zero();
        let mut aa = a;
        let mut tmp = spint_zero();
        while aa != 0 {
            if aa & 1 == 1 {
                modadd(&acc, &x, &mut tmp, &ctx.two_n, nw);
                acc = tmp;
            }
            aa >>= 1;
            if aa != 0 {
                let xc = x;
                modadd(&xc, &xc, &mut tmp, &ctx.two_n, nw);
                x = tmp;
            }
        }
        x = acc;
        let e = n.div_2exp(2);
        let xc = x;
        modxpowe_non_ct(&mut x, &xc, &e, &ctx.n, &ctx.two_n, ctx.ndash, nw, window);
        let mut y = spint_zero();
        modsqr(&x, &mut y, &ctx.n, ctx.ndash, nw);
        let yc = y;
        modadd(&yc, &ctx.one, &mut y, &ctx.two_n, nw);
        if !modis0(&y, &ctx.n, ctx.ndash, nw) {
            return None;
        }
        let mut res = spint_zero();
        redc(&x, &mut res, &ctx.n, ctx.ndash, nw);
        let mut r = Self::from_modq(&res, nw);
        r.bitlen = n.bitlen;
        Some(r)
    }
}
