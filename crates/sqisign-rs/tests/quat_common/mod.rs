//! Helpers for the quaternion-layer tests: parsing the harness's records
//! and the seeded PRNG domains.

#![allow(dead_code)]

use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use sqisign_rs::mp::{Ibz, Rng, ShakeRng};
use sqisign_rs::quat::integers::cornacchia_prime;
use sqisign_rs::quat::{Mat4x4, QuatAlg, QuatAlgElem, QuatIdeal, QuatLattice, Vec4};

/// The harness's `seed_default_domain`.
pub fn default_domain(label: &str) -> ShakeRng {
    let mut h = Shake256::default();
    h.update(b"prism-quat-prng-seed/");
    h.update(label.as_bytes());
    let mut seed = [0u8; 48];
    h.finalize_xof().read(&mut seed);
    ShakeRng::new(&seed, b"dom")
}

pub fn parse_ibz<const N: usize>(tok: &str) -> Ibz<N> {
    let (h, b) = tok.split_once(':').expect("hex:bitlen");
    let mut x =
        Ibz::<N>::from_str_radix(h, 16).unwrap_or_else(|| panic!("bad integer token `{tok}`"));
    x.set_bound(b.parse().expect("bitlen"));
    x
}

pub fn fmt_ibz<const N: usize>(x: &Ibz<N>) -> String {
    let mut buf = [0u8; 2 * 64 * 16 + 2];
    format!("{}:{}", x.to_hex(&mut buf), x.get_bound())
}

/// A cursor over the whitespace-separated tokens of a record.
pub struct Toks<'a> {
    toks: Vec<&'a str>,
    pos: usize,
}

impl<'a> Toks<'a> {
    pub fn new(s: &'a str) -> Self {
        Self {
            toks: s.split_whitespace().collect(),
            pos: 0,
        }
    }

    pub fn str(&mut self) -> &'a str {
        let t = self.toks[self.pos];
        self.pos += 1;
        t
    }

    pub fn int(&mut self) -> i64 {
        self.str().parse().expect("int")
    }

    pub fn ibz<const N: usize>(&mut self) -> Ibz<N> {
        parse_ibz(self.str())
    }

    pub fn elem<const N: usize>(&mut self) -> QuatAlgElem<N> {
        let denom = self.ibz();
        let c = [self.ibz(), self.ibz(), self.ibz(), self.ibz()];
        QuatAlgElem {
            denom,
            coord: Vec4(c),
        }
    }

    pub fn ideal<const N: usize>(&mut self) -> QuatIdeal<N> {
        QuatIdeal {
            x: self.ibz(),
            y: self.ibz(),
            norm: self.ibz(),
        }
    }

    pub fn lattice<const N: usize>(&mut self) -> QuatLattice<N> {
        let denom = self.ibz();
        QuatLattice {
            denom,
            basis: self.mat4(),
        }
    }

    pub fn mat4<const N: usize>(&mut self) -> Mat4x4<N> {
        let mut m = Mat4x4::zero();
        for i in 0..4 {
            for j in 0..4 {
                m.0[i][j] = self.ibz();
            }
        }
        m
    }

    pub fn done(&self) -> bool {
        self.pos >= self.toks.len()
    }
}

pub fn fmt_elem<const N: usize>(e: &QuatAlgElem<N>) -> String {
    let mut s = fmt_ibz(&e.denom);
    for c in &e.coord.0 {
        s.push(' ');
        s.push_str(&fmt_ibz(c));
    }
    s
}

pub fn fmt_ideal<const N: usize>(i: &QuatIdeal<N>) -> String {
    format!("{} {} {}", fmt_ibz(&i.x), fmt_ibz(&i.y), fmt_ibz(&i.norm))
}

pub fn fmt_mat4<const N: usize>(m: &Mat4x4<N>) -> String {
    let mut s = String::new();
    for i in 0..4 {
        for j in 0..4 {
            if !s.is_empty() {
                s.push(' ');
            }
            s.push_str(&fmt_ibz(&m.0[i][j]));
        }
    }
    s
}

pub fn fmt_lattice<const N: usize>(l: &QuatLattice<N>) -> String {
    format!("{} {}", fmt_ibz(&l.denom), fmt_mat4(&l.basis))
}

/// Port of the reference test helper `quat_represent_integer_even`
/// (`src/quaternion/ref/lvlx/test/helpers.c`): RepresentInteger for an even
/// target norm, drawing from `rng` exactly like the reference draws from its
/// default domain. Used to build challenge ideals of norm `2^f` the way the
/// reference protocol test does.
pub fn represent_integer_even<const N: usize>(
    n_gamma: &Ibz<N>,
    alg: &QuatAlg<N>,
    rng: &mut impl Rng,
) -> Option<QuatAlgElem<N>> {
    let p = &alg.p;
    let pbits = p.bitsize();
    let adjusted = n_gamma.mul_2exp(2);
    let (mut sq_bound, _) = adjusted.div(p);
    sq_bound = sq_bound.sub(&Ibz::one());
    let mut bound = sq_bound.sqrt_floor();
    bound.set_bound((adjusted.get_bound() - pbits) / 2 + 3);
    let temp = p.mul(p).sqrt_floor();
    let (mut counter, _) = adjusted.div(&temp);
    counter.set_bound(pbits * 3 / 2);
    let mut coeffs = Vec4::<N>::zero();
    let mut gamma = QuatAlgElem::<N>::zero();
    let mut found = false;
    while !found && !counter.is_zero() {
        counter = counter.sub(&Ibz::one());
        counter.set_bound(pbits * 3 / 2);
        coeffs.0[2] = Ibz::rand_interval(&Ibz::one(), &bound, rng)?;
        let mut target = coeffs.0[2].mul(&coeffs.0[2]);
        let mut temp = adjusted.sub(&target.mul(p));
        temp = temp.div(p).0;
        temp.set_bound(adjusted.get_bound() - pbits + 2);
        temp = temp.sqrt_floor();
        if temp <= Ibz::zero() {
            continue;
        }
        coeffs.0[3] = Ibz::rand_interval(&Ibz::one(), &temp, rng)?;
        temp = coeffs.0[3].mul(&coeffs.0[3]);
        target = adjusted.sub(&target.add(&temp).mul(p));
        target.set_bound(n_gamma.get_bound() + 2);
        assert!(target > Ibz::zero());
        found = false;
        if target.probab_prime(alg.primality_num_iter, rng) {
            if let Some((x, y)) = cornacchia_prime(&target) {
                coeffs.0[0] = x;
                coeffs.0[1] = y;
                found = true;
            }
        }
        if found {
            for c in coeffs.0[..2].iter_mut() {
                c.set_bound(target.get_bound() / 2 + 2);
            }
            gamma.coord = coeffs;
            gamma.denom = Ibz::set(1, 2);
            let (prim, content) = gamma.make_primitive(&alg.o0);
            coeffs = prim;
            found = content == Ibz::two();
        }
    }
    if !found {
        return None;
    }
    gamma.coord = alg.o0.basis.eval(&coeffs);
    gamma.denom = alg.o0.denom;
    Some(gamma)
}
