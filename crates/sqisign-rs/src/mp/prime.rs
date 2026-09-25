//! Probable-prime testing: the reference's small-prime sieve, a strong
//! Miller-Rabin round in base 2, a strong Lucas-Selfridge test, and extra
//! random Miller-Rabin rounds beyond the first 24 (variable time).

use super::modq::{
    modadd, modmul, modsqr, modsub, modxpowe_non_ct, nresx, spint_zero, ModCtx, Spint, LIMBMASK,
    RADIX,
};
use super::rand::{DefaultDomain, Rng};
use super::Ibz;

/// Small odd primes up to 1021, the sieving table.
pub(crate) const AGGRESSIVE_PRIMES: [u32; 171] = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
    431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
    661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
    809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929,
    937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021,
];

const JACOBI_BASES: [u32; 34] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139,
];

struct SieveBatch {
    d: u64,
    di: u64,
    shift: u32,
    first: usize,
    last: usize,
}

const SPLIT_A: u64 = 0xf600_59af; // 3*5*7*11*13*17*19*23*37
const SPLIT_B: u64 = 0xeb52_e3f3; // 29*31*41*43*47*53

const SIEVE_COMMON: [SieveBatch; 10] = [
    SieveBatch {
        d: 0x6329899ea9f2714b,
        di: 0x4a72c477c0963cdb,
        shift: 1,
        first: 15,
        last: 25,
    },
    SieveBatch {
        d: 0x58edcb4c9ed39c8b,
        di: 0x707965cc3eac798e,
        shift: 1,
        first: 25,
        last: 34,
    },
    SieveBatch {
        d: 0x09966ff94fd516fb,
        di: 0xab37686a5d6e1835,
        shift: 4,
        first: 34,
        last: 42,
    },
    SieveBatch {
        d: 0x3bd7632c1f36eb51,
        di: 0x11ca639fffb1f62a,
        shift: 2,
        first: 42,
        last: 50,
    },
    SieveBatch {
        d: 0x00fd14b3c90d88a9,
        di: 0x02f3ead63d941e0e,
        shift: 8,
        first: 50,
        last: 57,
    },
    SieveBatch {
        d: 0x02ad3dbe0cca85ff,
        di: 0x7e8ea9abc0c5a2d6,
        shift: 6,
        first: 57,
        last: 64,
    },
    SieveBatch {
        d: 0x0787f9a02c3388a7,
        di: 0x0fefe739af03fd10,
        shift: 5,
        first: 64,
        last: 71,
    },
    SieveBatch {
        d: 0x1113c5cc6d101657,
        di: 0xdfb3f0f9b8c8b7b1,
        shift: 3,
        first: 71,
        last: 78,
    },
    SieveBatch {
        d: 0x2456c94f936bdb15,
        di: 0xc2dd8758870f9711,
        shift: 2,
        first: 78,
        last: 85,
    },
    SieveBatch {
        d: 0x4236a30b85ffe139,
        di: 0xeee270a9d031bebc,
        shift: 1,
        first: 85,
        last: 92,
    },
];
const SIEVE_509: [SieveBatch; 1] = [SieveBatch {
    d: 0x0000000e9aef58cb,
    di: 0x1872a19b47be8fb4,
    shift: 28,
    first: 92,
    last: 96,
}];
const SIEVE_1021: [SieveBatch; 13] = [
    SieveBatch {
        d: 0x805437b38eada69d,
        di: 0xfeaffe452003a47e,
        shift: 0,
        first: 92,
        last: 99,
    },
    SieveBatch {
        d: 0x00723e97bddcd2af,
        di: 0x1ed2cc34bc0c181d,
        shift: 9,
        first: 99,
        last: 105,
    },
    SieveBatch {
        d: 0x00a5a792ee239667,
        di: 0x8b9e4e9bcd5b9458,
        shift: 8,
        first: 105,
        last: 111,
    },
    SieveBatch {
        d: 0x00e451352ebca269,
        di: 0x1f0a0b7149843601,
        shift: 8,
        first: 111,
        last: 117,
    },
    SieveBatch {
        d: 0x013a7955f14b7805,
        di: 0xa0cc308fcfbd3e0d,
        shift: 7,
        first: 117,
        last: 123,
    },
    SieveBatch {
        d: 0x01d37cbd653b06ff,
        di: 0x1860242e63e661d6,
        shift: 7,
        first: 123,
        last: 129,
    },
    SieveBatch {
        d: 0x0288fe4eca4d7cdf,
        di: 0x93ec8aa98a84de0b,
        shift: 6,
        first: 129,
        last: 135,
    },
    SieveBatch {
        d: 0x039fddb60d3af63d,
        di: 0x1a8606753dac09b2,
        shift: 6,
        first: 135,
        last: 141,
    },
    SieveBatch {
        d: 0x04cd73f19080fb03,
        di: 0xaa70a958c486716c,
        shift: 5,
        first: 141,
        last: 147,
    },
    SieveBatch {
        d: 0x0639c390b9313f05,
        di: 0x48f64f21c7ab90db,
        shift: 5,
        first: 147,
        last: 153,
    },
    SieveBatch {
        d: 0x08a1c420d25d388f,
        di: 0xda84dcd209b64d3f,
        shift: 4,
        first: 153,
        last: 159,
    },
    SieveBatch {
        d: 0x0b4b5322977db499,
        di: 0x6aa9aee2cab2fcf9,
        shift: 4,
        first: 159,
        last: 165,
    },
    SieveBatch {
        d: 0x0e94c170a802ee29,
        di: 0x18e97b7f459899fe,
        shift: 4,
        first: 165,
        last: 171,
    },
];

/// `n mod m` for a single-limb `m`.
pub(crate) fn mod_digit<const N: usize>(n: &Ibz<N>, m: u64) -> u64 {
    debug_assert!(m != 0 && n.is_positive());
    let mut r = 0u64;
    for i in (0..n.n()).rev() {
        let x = ((r as u128) << 64) + n.limbs[i] as u128;
        r = (x % m as u128) as u64;
    }
    r
}

/// A small prime `a` with Jacobi symbol `(a / n) = -1`, or `None`.
pub(crate) fn find_jacobi_minus_one<const N: usize>(n: &Ibz<N>) -> Option<u32> {
    for &a in JACOBI_BASES.iter() {
        if a == 2 {
            if (n.limbs[0] & 7) == 5 {
                return Some(2);
            }
            continue;
        }
        let r = mod_digit(n, a as u64) as u32;
        if r == 0 {
            continue;
        }
        let mut e = (a - 1) >> 1;
        let mut base = r as u64;
        let mut acc = 1u64;
        while e != 0 {
            if e & 1 == 1 {
                acc = (acc * base) % a as u64;
            }
            base = (base * base) % a as u64;
            e >>= 1;
        }
        if acc == (a - 1) as u64 {
            return Some(a);
        }
    }
    None
}

fn mod_split<const N: usize>(n: &Ibz<N>, m: u64) -> u32 {
    let mut r = 0u64;
    for i in (0..n.n()).rev() {
        let w = n.limbs[i];
        r = ((r << 32) | (w >> 32)) % m;
        r = ((r << 32) | (w & 0xffff_ffff)) % m;
    }
    r as u32
}

fn rem_2by1_preinv(nh: u64, nl: u64, d: u64, di: u64) -> u64 {
    let p = (nh as u128) * (di as u128);
    let mut qh = (p >> 64) as u64;
    let mut ql = p as u64;
    let suml = ql.wrapping_add(nl);
    let carry = (suml < ql) as u64;
    ql = suml;
    qh = qh.wrapping_add(nh).wrapping_add(1).wrapping_add(carry);
    let mut r = nl.wrapping_sub(qh.wrapping_mul(d));
    if r > ql {
        r = r.wrapping_add(d);
    }
    if r >= d {
        r -= d;
    }
    r
}

fn mod_u64_preinv<const N: usize>(n: &Ibz<N>, b: &SieveBatch) -> u64 {
    let mut r = 0u64;
    let s = b.shift;
    let dn = b.d << s;
    for i in (0..n.n()).rev() {
        let w = n.limbs[i];
        let (nh, nl) = if s == 0 {
            (r, w)
        } else {
            ((r << s) | (w >> (64 - s)), w << s)
        };
        let rn = rem_2by1_preinv(nh, nl, dn, b.di);
        r = if s == 0 { rn } else { rn >> s };
    }
    r
}

fn sieve_fast<const N: usize>(n: &Ibz<N>, sieve_bound: u32) -> bool {
    let ra = mod_split(n, SPLIT_A);
    if [3u32, 5, 7, 11, 13, 17, 19, 23, 37]
        .iter()
        .any(|&q| ra % q == 0)
    {
        return false;
    }
    let rb = mod_split(n, SPLIT_B);
    if [29u32, 31, 41, 43, 47, 53].iter().any(|&q| rb % q == 0) {
        return false;
    }
    let tail: &[SieveBatch] = if sieve_bound == 509 {
        &SIEVE_509
    } else {
        &SIEVE_1021
    };
    for b in SIEVE_COMMON.iter().chain(tail.iter()) {
        let r = mod_u64_preinv(n, b);
        for k in b.first..b.last {
            if r % AGGRESSIVE_PRIMES[k] as u64 == 0 {
                return false;
            }
        }
    }
    true
}

fn sieve<const N: usize>(n: &Ibz<N>, sieve_bound: u32) -> bool {
    if n.bitsize() > 32 && (sieve_bound == 509 || sieve_bound == 1021) {
        return sieve_fast(n, sieve_bound);
    }
    let count = AGGRESSIVE_PRIMES.len();
    let mut i = 0;
    while i < count && AGGRESSIVE_PRIMES[i] <= sieve_bound {
        let mut product = 1u64;
        let mut j = i;
        while j < count && AGGRESSIVE_PRIMES[j] <= sieve_bound {
            let q = AGGRESSIVE_PRIMES[j] as u64;
            if product > u64::MAX / q {
                break;
            }
            product *= q;
            j += 1;
        }
        let r = mod_digit(n, product);
        for k in i..j {
            let q = AGGRESSIVE_PRIMES[k] as u64;
            if r % q == 0 {
                if n.bitsize() <= 32 && n.limbs[0] == q {
                    continue;
                }
                return false;
            }
        }
        i = j;
    }
    true
}

struct PrimeCtx<const N: usize> {
    m: ModCtx,
    minus_one: Spint,
    d: Ibz<N>,
    s: i32,
}

impl<const N: usize> PrimeCtx<N> {
    fn new(n: &Ibz<N>) -> Option<Self> {
        let numwords = ModCtx::numwords_for(n.bitlen);
        if numwords == 0 || numwords >= super::modq::MAX_MODQ_LIMBS {
            return None;
        }
        let m = ModCtx::new(n);
        let zero = spint_zero();
        let mut minus_one = spint_zero();
        modsub(&zero, &m.one, &mut minus_one, &m.two_n, m.numwords);
        let nm1 = n.sub(&Ibz::one());
        let s = nm1.two_adic();
        let d = nm1.div_2exp(s as u32);
        Some(Self { m, minus_one, d, s })
    }

    /// Fully reduce a residue below `n`.
    fn reduce_once(&self, v: &mut [u64]) {
        let nw = self.m.numwords;
        let mut ge = 0i32;
        for i in (0..nw).rev() {
            let vi = v[i] & LIMBMASK;
            let ni = self.m.n[i] & LIMBMASK;
            if vi > ni {
                ge = 1;
                break;
            }
            if vi < ni {
                ge = -1;
                break;
            }
        }
        if ge >= 0 {
            let base: u128 = 1u128 << RADIX;
            let mut borrow = 0u128;
            for i in 0..nw {
                let vi = (v[i] & LIMBMASK) as u128;
                let ni = (self.m.n[i] & LIMBMASK) as u128 + borrow;
                if vi >= ni {
                    v[i] = ((vi - ni) as u64) & LIMBMASK;
                    borrow = 0;
                } else {
                    v[i] = ((base + vi - ni) as u64) & LIMBMASK;
                    borrow = 1;
                }
            }
        }
    }

    fn mont_equal(&self, a: &[u64], b: &[u64]) -> bool {
        let nw = self.m.numwords;
        let mut aa = spint_zero();
        let mut bb = spint_zero();
        aa[..nw].copy_from_slice(&a[..nw]);
        bb[..nw].copy_from_slice(&b[..nw]);
        self.reduce_once(&mut aa);
        self.reduce_once(&mut bb);
        (0..nw).all(|i| (aa[i] & LIMBMASK) == (bb[i] & LIMBMASK))
    }

    fn mr_round(&self, x: &mut Spint, window: u32) -> bool {
        let nw = self.m.numwords;
        let xc = *x;
        modxpowe_non_ct(
            x,
            &xc,
            &self.d,
            &self.m.n,
            &self.m.two_n,
            self.m.ndash,
            nw,
            window,
        );
        if self.mont_equal(x, &self.m.one) || self.mont_equal(x, &self.minus_one) {
            return true;
        }
        for _ in 1..self.s {
            let xc = *x;
            modsqr(&xc, x, &self.m.n, self.m.ndash, nw);
            if self.mont_equal(x, &self.minus_one) {
                return true;
            }
            if self.mont_equal(x, &self.m.one) {
                return false;
            }
        }
        false
    }

    fn half_mod(&self, a: &mut [u64]) {
        let nw = self.m.numwords;
        self.reduce_once(a);
        if a[0] & 1 == 1 {
            let mut carry = 0u64;
            for i in 0..nw {
                let z =
                    (a[i] & LIMBMASK) as u128 + (self.m.n[i] & LIMBMASK) as u128 + carry as u128;
                a[i] = (z as u64) & LIMBMASK;
                carry = (z >> RADIX) as u64;
            }
        }
        let mut carry = 0u64;
        for i in (0..nw).rev() {
            let next = a[i] & 1;
            a[i] = ((a[i] >> 1) | (carry << (RADIX - 1))) & LIMBMASK;
            carry = next;
        }
    }

    fn scale_small(&self, a: &[u64], k: i64, out: &mut [u64]) {
        let nw = self.m.numwords;
        let mut acc = spint_zero();
        let mut base = spint_zero();
        let mut tmp = spint_zero();
        base[..nw].copy_from_slice(&a[..nw]);
        let mut mag = k.unsigned_abs();
        while mag != 0 {
            if mag & 1 == 1 {
                modadd(&acc, &base, &mut tmp, &self.m.two_n, nw);
                acc = tmp;
            }
            mag >>= 1;
            if mag != 0 {
                let bc = base;
                modadd(&bc, &bc, &mut tmp, &self.m.two_n, nw);
                base = tmp;
            }
        }
        if k < 0 {
            let zero = spint_zero();
            modsub(&zero, &acc, out, &self.m.two_n, nw);
        } else {
            out[..nw].copy_from_slice(&acc[..nw]);
        }
    }

    fn lucas_double(&self, v: &mut Spint, qk: &mut Spint) {
        let nw = self.m.numwords;
        let mut v2 = spint_zero();
        let mut twoq = spint_zero();
        let mut q2 = spint_zero();
        modsqr(v, &mut v2, &self.m.n, self.m.ndash, nw);
        modadd(qk, qk, &mut twoq, &self.m.two_n, nw);
        modsub(&v2, &twoq, v, &self.m.two_n, nw);
        modsqr(qk, &mut q2, &self.m.n, self.m.ndash, nw);
        *qk = q2;
    }

    fn strong_lucas(&self, n: &Ibz<N>) -> bool {
        if perfect_square(n) {
            return false;
        }
        let mut dabs = 5u64;
        let q: i64;
        loop {
            let r = mod_digit(n, dabs);
            if r == 0 {
                return false;
            }
            let j = jacobi_u64(r, dabs);
            if j == 0 {
                return false;
            }
            if j == -1 {
                q = if dabs & 2 != 0 {
                    (dabs >> 2) as i64 + 1
                } else {
                    -((dabs >> 2) as i64)
                };
                break;
            }
            dabs += 2;
            if dabs > 1_000_001 {
                return false;
            }
        }
        let np1 = n.add(&Ibz::one());
        let s = np1.two_adic();
        let d = np1.div_2exp(s as u32);
        let nw = self.m.numwords;
        let mut u = self.m.one;
        let mut v = self.m.one;
        let mut qk = spint_zero();
        self.scale_small(&self.m.one, q, &mut qk);
        let qmont = qk;
        let zero = spint_zero();
        let nb = d.bitsize();
        for pos in (0..nb - 1).rev() {
            let mut u2 = spint_zero();
            modmul(&u, &v, &mut u2, &self.m.n, self.m.ndash, nw);
            self.lucas_double(&mut v, &mut qk);
            u = u2;
            if d.bit_at(pos) == 1 {
                let old_u = u;
                let mut sum = spint_zero();
                modadd(&u, &v, &mut sum, &self.m.two_n, nw);
                u = sum;
                self.half_mod(&mut u);
                let mut term = spint_zero();
                self.scale_small(&old_u, -2 * q, &mut term);
                let mut vnew = spint_zero();
                modadd(&u, &term, &mut vnew, &self.m.two_n, nw);
                v = vnew;
                let mut qnew = spint_zero();
                modmul(&qk, &qmont, &mut qnew, &self.m.n, self.m.ndash, nw);
                qk = qnew;
            }
        }
        if self.mont_equal(&u, &zero) {
            return true;
        }
        for r in 0..s {
            if self.mont_equal(&v, &zero) {
                return true;
            }
            if r + 1 < s {
                self.lucas_double(&mut v, &mut qk);
            }
        }
        false
    }
}

fn jacobi_u64(mut a: u64, mut n: u64) -> i32 {
    if n == 0 || n & 1 == 0 {
        return 0;
    }
    a %= n;
    let mut sign = 1;
    while a != 0 {
        while a & 1 == 0 {
            a >>= 1;
            let r = n & 7;
            if r == 3 || r == 5 {
                sign = -sign;
            }
        }
        core::mem::swap(&mut a, &mut n);
        if (a & 3) == 3 && (n & 3) == 3 {
            sign = -sign;
        }
        a %= n;
    }
    if n == 1 {
        sign
    } else {
        0
    }
}

fn perfect_square<const N: usize>(n: &Ibz<N>) -> bool {
    for &m in &[8u32, 3, 5, 7, 11, 13] {
        let r = mod_digit(n, m as u64) as u32;
        if !(0..m).any(|x| (x * x) % m == r) {
            return false;
        }
    }
    let root = n.sqrt_floor();
    root.mul(&root) == *n
}

impl<const N: usize> Ibz<N> {
    /// Probable-prime test with `reps` rounds: sieve, Miller-Rabin base 2,
    /// strong Lucas-Selfridge, then `max(reps - 24, 0)` random Miller-Rabin
    /// rounds drawn from `rng`. `n` must be at least 66 bits below the
    /// container. Variable time.
    pub fn probab_prime(&self, reps: u32, rng: &mut impl Rng) -> bool {
        debug_assert!(reps > 0);
        debug_assert!(self.bitlen + 65 < Self::MAX_BITS);
        let mut arg = *self;
        if self.is_positive() {
            let actual_bits = self.bitsize();
            if actual_bits + 1 < self.bitlen {
                arg = Self::from_bits(&self.limbs, actual_bits);
            }
        }
        bpsw(&arg, reps, rng)
    }
}

fn bpsw<const N: usize>(n: &Ibz<N>, reps: u32, rng: &mut impl Rng) -> bool {
    if *n < Ibz::two() {
        return false;
    }
    if *n == Ibz::two() || *n == Ibz::three() {
        return true;
    }
    if n.is_even() {
        return false;
    }
    let modqbits = Ibz::<N>::MAX_BITS - 1;
    let (sieve_bound, window) = if modqbits < 1200 {
        (509u32, 4u32)
    } else if modqbits < 1800 {
        (1021, 6)
    } else {
        (509, 6)
    };
    if !sieve(n, sieve_bound) {
        return false;
    }
    if n.bitsize() <= 32 && n.limbs[0] <= sieve_bound as u64 {
        return true;
    }
    let Some(ctx) = PrimeCtx::new(n) else {
        return false;
    };
    let nw = ctx.m.numwords;
    let mut x = spint_zero();
    modadd(&ctx.m.one, &ctx.m.one, &mut x, &ctx.m.two_n, nw);
    if !ctx.mr_round(&mut x, window) {
        return false;
    }
    if !ctx.strong_lucas(n) {
        return false;
    }
    let extra = reps.saturating_sub(24);
    if extra == 0 {
        return true;
    }
    let mut r2 = spint_zero();
    nresx(&ctx.m.one, &mut r2, &ctx.m.two_n, nw);
    let n_m2 = n.sub(&Ibz::two());
    for _ in 0..extra {
        let Some(a) = Ibz::rand_interval(&Ibz::two(), &n_m2, &mut DefaultDomain(&mut *rng)) else {
            return false;
        };
        let normal = a.to_modq();
        let mut witness = spint_zero();
        modmul(&normal, &r2, &mut witness, &ctx.m.n, ctx.m.ndash, nw);
        if !ctx.mr_round(&mut witness, window) {
            return false;
        }
    }
    true
}
