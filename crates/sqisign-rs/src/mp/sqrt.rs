//! Variable-time floor square root: Newton's reciprocal square root at
//! doubling widths on the normalised input, then an exact correction.

use super::arith::{scratch, shiftl, shiftr};
use super::div::{add_raw, compare, mul_raw, sub_raw};
use super::{nlimbs, Ibz};

/// Top `keep` limbs of `a * b`, one-sided (never above the true value).
fn mul_high(r: &mut [u64], a: &[u64], la: usize, b: &[u64], lb: usize, keep: usize) {
    let c0 = (la + lb).saturating_sub(keep + 2);
    let next = la + lb - c0;
    let mut ext = scratch();
    for i in 0..la {
        let j0 = c0.saturating_sub(i);
        if j0 >= lb {
            continue;
        }
        let mut carry = 0u128;
        for j in j0..lb {
            let sum = (a[i] as u128) * (b[j] as u128) + ext[i + j - c0] as u128 + carry;
            ext[i + j - c0] = sum as u64;
            carry = sum >> 64;
        }
        ext[i + lb - c0] = carry as u64;
    }
    r[..keep].copy_from_slice(&ext[next - keep..next]);
}

/// Reciprocal square root seed of one normalised limb, in Q2.62.
fn digit_rsqrt_seed(a1: u64) -> u64 {
    debug_assert!(a1 >= 1u64 << 62);
    let mut y = 1u64 << 62;
    for _ in 0..8 {
        let p = (((y as u128) * (y as u128)) >> 64) as u64;
        let q = (((p as u128) * (a1 as u128)) >> 64) as u64;
        let h = (3u64 << 60).wrapping_sub(q);
        y = (((y as u128) * (h as u128)) >> 61) as u64;
    }
    y - 16
}

fn bitlen_non_ct(x: &[u64], nwords: usize) -> i32 {
    for i in (0..nwords).rev() {
        if x[i] != 0 {
            return (i as i32) * 64 + (64 - x[i].leading_zeros() as i32);
        }
    }
    0
}

fn shl1_bit(x: &mut [u64], nwords: usize, bit_in: u64) {
    let mut carry = bit_in;
    for v in x[..nwords].iter_mut() {
        let old = *v;
        *v = (old << 1) | carry;
        carry = old >> 63;
    }
}

impl<const N: usize> Ibz<N> {
    /// `floor(sqrt(a))` for non-negative `a`. Variable time.
    pub fn sqrt_floor(&self) -> Self {
        debug_assert!(self.is_positive());
        let b_eff = bitlen_non_ct(&self.limbs, self.n());
        if b_eff == 0 {
            return Self::zero();
        }
        let n_eff = nlimbs(b_eff);
        let m = (n_eff + 1) / 2;

        // Step 1: normalise into A (2m + 2 limbs, two guard limbs below the value)
        let z = 2 * 64 * m as i32 - b_eff;
        let t = z & !1;
        let th = z >> 1;
        let off = 2 * m - n_eff;
        let tb = t - 64 * off as i32;
        let mut a = scratch();
        a[2 + off..2 + off + n_eff].copy_from_slice(&self.limbs[..n_eff]);
        if tb != 0 {
            shiftl(&mut a[2..2 + 2 * m], tb as u32);
        }
        debug_assert!(a[2 * m + 1] >= 1u64 << 62);

        // Step 2: seed at L = 1 embedded into L = 2
        let mut y = scratch();
        y[0] = 0;
        y[1] = digit_rsqrt_seed(a[2 * m + 1]);
        let mut l = 2usize;

        // Step 3: Newton stages
        let mut s = scratch();
        let mut tt = scratch();
        let mut h = scratch();
        let mut u = scratch();
        let margin = [64u64];
        let mut lout = 2usize;
        loop {
            debug_assert!(lout < 2 * l && lout <= m + 1);
            mul_raw(&mut s, &y, l, &y, l);
            let a_t = &a[2 * m + 2 - 2 * l..2 * m + 2];
            mul_high(&mut tt, a_t, 2 * l, &s, 2 * l, 2 * l + 2);
            for v in h[..2 * l + 2].iter_mut() {
                *v = 0;
            }
            h[2 * l + 1] = 3u64 << 60;
            debug_assert!(compare(&h, &tt, 2 * l + 2) != core::cmp::Ordering::Less);
            let hc = h;
            sub_raw(&mut h, 2 * l + 2, &hc, 2 * l + 2, &tt, 2 * l + 2);
            mul_raw(&mut u, &y, l, &h[2..], 2 * l);
            let win_start = 3 * l - lout - 1;
            shiftr(&mut u[win_start..win_start + lout + 1], 61);
            debug_assert_eq!(u[win_start + lout], 0);
            y[..lout].copy_from_slice(&u[win_start..win_start + lout]);
            let yc = y;
            sub_raw(&mut y, lout, &yc, lout, &margin, 1);
            l = lout;
            lout = if 2 * l - 1 < m + 1 { 2 * l - 1 } else { m + 1 };
            if l > m {
                break;
            }
        }

        // Step 4: multiply through, then correct
        let mut p = scratch();
        mul_raw(&mut p, &a[m + 1..2 * m + 2], m + 1, &y, m + 1);
        shiftr(&mut p[m + 1..2 * m + 2], 62);
        debug_assert_eq!(p[2 * m + 1], 0);
        let mut sh = scratch();
        sh[..m].copy_from_slice(&p[m + 1..2 * m + 1]);

        let mut av = scratch();
        av[..2 * m].copy_from_slice(&a[2..2 + 2 * m]);
        let mut ssq = scratch();
        let mut d = scratch();
        mul_raw(&mut ssq, &sh, m, &sh, m);
        debug_assert!(compare(&av, &ssq, 2 * m) != core::cmp::Ordering::Less);
        let avc = av;
        sub_raw(&mut av, 2 * m, &avc, 2 * m, &ssq, 2 * m);
        d[..m].copy_from_slice(&sh[..m]);
        shl1_bit(&mut d, 2 * m, 1);
        let one = [1u64];
        let two = [2u64];
        let mut corrections = 0;
        while compare(&av, &d, 2 * m) != core::cmp::Ordering::Less {
            let avc = av;
            sub_raw(&mut av, 2 * m, &avc, 2 * m, &d, 2 * m);
            let shc = sh;
            add_raw(&mut sh, m, &shc, m, &one, 1);
            let dc = d;
            add_raw(&mut d, 2 * m, &dc, 2 * m, &two, 1);
            corrections += 1;
        }
        debug_assert!(corrections <= 2);

        // Step 5: denormalise
        if th != 0 {
            shiftr(&mut sh[..m], th as u32);
        }
        let out = Self::from_bits(&sh[..m], (b_eff + 1) / 2);
        #[cfg(debug_assertions)]
        {
            let chk = out.mul(&out);
            let r = self.sub(&chk);
            debug_assert!(r.is_positive());
            let r2 = r.sub(&out).sub(&out);
            debug_assert!(!r2.is_positive() || r2.is_zero());
        }
        out
    }
}
