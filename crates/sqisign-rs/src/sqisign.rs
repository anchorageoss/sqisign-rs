//! Round-3 key generation and signing (`keygen.c`, `sign.c`,
//! `encode_secret.c`), statement for statement; verification and the
//! public encodings are [`sqisign_verify::sqisign`].

use crate::id2iso::{
    arbitrary_isogeny_evaluation, change_of_basis_matrix_tate, change_of_basis_matrix_tate_invert,
    endomorphism_application_even_basis, ideal_to_isogeny_qlapoty, kernel_dlogs_to_ideal_even,
    matrix_application_even_basis, E0Actions, Id2IsoStats,
};
use crate::mp::{Ibz, Rng, ShakeRng, SplitRng};
use crate::quat::protocol::response_element;
use crate::quat::{Mat2x2, QuatAlg, QuatIdeal, Vec2};
use alloc::vec::Vec;
use sqisign_verify::ec::basis::ec_curve_to_basis_2f_to_hint;
use sqisign_verify::ec::point::{ec_dbl_iter_basis, ec_mul};
use sqisign_verify::ec::{EcBasis, EcCurve, MAX_ORDER_WORDS};
use sqisign_verify::fp::FpBackend;
use sqisign_verify::params::Prime;
use sqisign_verify::precomp::{PrimePrecomp, EXTRA_TORSION};
use sqisign_verify::sqisign::{
    hash_to_challenge, public_key_from_bytes, public_key_to_bytes, PublicKey, Scalar, Signature,
    VerifyParams, HD_EXTRA_TORSION, MAX_PUBLICKEY_BYTES,
};
use sqisign_verify::theta::chain::theta_chain_compute_and_eval;
use sqisign_verify::theta::{
    ChainMode, ThetaCoupleCurve, ThetaCouplePoint, ThetaKernelCouplePoints,
};
use zeroize::{Zeroize, ZeroizeOnDrop, Zeroizing};

/// Everything a level's key generation and signing need: the verifier's
/// constants plus the quaternion algebra and the endomorphism data.
pub struct Params<const N: usize> {
    /// The quaternion algebra and its lattice-reduction configuration.
    pub alg: QuatAlg<N>,
    /// The action of `O0` on the precomputed basis of `E0[2^f]`.
    pub act: E0Actions<N>,
    /// `SEC_DEGREE = COM_DEGREE = D_mix`, with the reference's bound.
    pub sec_degree: Ibz<N>,
    /// The verifier's constants and the encoded sizes.
    pub verify: VerifyParams,
}

macro_rules! level {
    ($(#[$m:meta])* $name:ident, $modname:ident, $qlevel:ident, $vlevel:ident) => {
        $(#[$m])*
        pub fn $name() -> Params<{ crate::precomp::$modname::IBZ_NLIMBS }> {
            use sqisign_verify::params::sqisign_v3::$modname as l;
            // the reference declares SEC_DEGREE with bound DEGREE_NORM_BITS + 1
            let mut sec_degree = Ibz::from_le_bytes(&l::DEG_MIX_LE);
            sec_degree.set_bound(l::DEGREE_NORM_BITS as i32 + 1);
            Params {
                alg: crate::quat::params::$qlevel(),
                act: crate::id2iso::params::$qlevel(),
                sec_degree,
                verify: sqisign_verify::sqisign::$vlevel(),
            }
        }
    };
}

level!(
    /// NIST level I, `p324_3`.
    level1, p324_3, level1, level1
);
level!(
    /// NIST level III, `p500_27`.
    level3, p500_27, level3, level3
);
level!(
    /// NIST level V, `p664_17`.
    level5, p664_17, level5, level5
);

/// `secret_key_t` (the canonical basis is not kept: signing recomputes
/// nothing from it). Zeroized on drop; `Debug` is redacted.
#[derive(Clone)]
pub struct SecretKey<L: Prime, const N: usize> {
    /// `E_pk`.
    pub curve: EcCurve<L>,
    /// The secret ideal in inert form.
    pub secret_ideal: QuatIdeal<N>,
    /// Change of basis from the canonical basis of `E_pk[2^(CHALLENGE_BITS
    /// + 2)]` to the image of `E0`'s basis, reduced modulo
    /// `2^CHALLENGE_BITS` and normalised.
    pub mat_bacan_to_ba0_two: Mat2x2<N>,
}

impl<L: Prime, const N: usize> Zeroize for SecretKey<L, N> {
    fn zeroize(&mut self) {
        // the curve is the public key; the ideal and the matrix are secret
        self.secret_ideal.zeroize();
        self.mat_bacan_to_ba0_two.zeroize();
    }
}

impl<L: Prime, const N: usize> Drop for SecretKey<L, N> {
    fn drop(&mut self) {
        self.zeroize();
    }
}

impl<L: Prime, const N: usize> ZeroizeOnDrop for SecretKey<L, N> {}

impl<L: Prime, const N: usize> core::fmt::Debug for SecretKey<L, N> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("SecretKey([REDACTED])")
    }
}

/// `prng_seed`: the 48-byte root seed, from which the domain streams are
/// derived (`prng_domain_seed`).
fn prng_seed(entropy: &mut impl Rng) -> Option<Zeroizing<[u8; 48]>> {
    let mut seed = Zeroizing::new([0u8; 48]);
    entropy.fill(&mut *seed).then_some(seed)
}

fn domain(seed: &[u8; 48], dom: &[u8; 3]) -> ShakeRng {
    ShakeRng::new(seed, dom)
}

fn digits<const N: usize>(x: &Ibz<N>) -> Scalar {
    let mut d = [0u64; MAX_ORDER_WORDS];
    x.to_digits(&mut d);
    d
}

/// `protocols_keygen`. `entropy` stands for `randombytes`: it supplies the
/// 48-byte root seed and nothing else. `None` only on an entropy failure
/// or an impossible discrete logarithm.
pub fn keygen<L: FpBackend + PrimePrecomp, const N: usize>(
    params: &Params<N>,
    entropy: &mut impl Rng,
) -> Option<(PublicKey<L>, SecretKey<L, N>)> {
    let alg = &params.alg;
    let act = &params.act;
    let vp = &params.verify;
    let seed = prng_seed(entropy)?;
    let mut def = domain(&seed, b"def");
    let mut key = SplitRng {
        domain: domain(&seed, b"key"),
        default: &mut def,
    };
    debug_assert_eq!(vp.torsion_even_power, L::TWO_ADIC_EXPONENT);

    // iterating until a solution has been found
    let (mut secret_ideal, mut curve, b_0_two) = loop {
        let Some(ideal) = QuatIdeal::random_given_prime_norm(&params.sec_degree, alg, &mut key)
            .map(Zeroizing::new)
        else {
            continue;
        };
        // replacing the secret key ideal by a shorter equivalent one
        let Some((beta, ideal)) = ideal.small_equivalent_coprime(Some(&Ibz::zero()), alg, &mut key)
        else {
            continue;
        };
        drop(Zeroizing::new(beta));
        let ideal = Zeroizing::new(ideal);
        let Some((curve, basis)) = arbitrary_isogeny_evaluation::<L, N>(&ideal, alg, act, &mut key)
        else {
            continue;
        };
        break (ideal, curve, Zeroizing::new(basis));
    };

    // a deterministic basis with a hint, and the change of basis to the
    // image of E0's basis
    let f_chall = vp.challenge_bits + EXTRA_TORSION;
    let (canonical_basis, hint_pk) = ec_curve_to_basis_2f_to_hint(&mut curve, f_chall, 0)?;
    let mut mat: Zeroizing<Mat2x2<N>> = Zeroizing::new(change_of_basis_matrix_tate(
        &canonical_basis,
        &b_0_two,
        &mut curve,
        f_chall,
    )?);
    for row in mat.0.iter_mut() {
        for x in row.iter_mut() {
            *x = x.mod2exp(vp.challenge_bits);
        }
    }
    mat.normalize(vp.challenge_bits);

    // the public key is the curve without precomputation
    let mut pk_curve = curve.clone();
    pk_curve.is_a24_computed_and_normalized = false;
    debug_assert!(bool::from(pk_curve.c.ct_is_one()));

    // the known bounds
    let fp_bits = 8 * vp.fp_encoded_bytes as i32;
    secret_ideal.norm.set_bound(fp_bits + 1);
    secret_ideal.x.set_bound(fp_bits + 1);
    secret_ideal.y.set_bound(fp_bits + 1);
    for row in mat.0.iter_mut() {
        for x in row.iter_mut() {
            x.set_bound(8 * vp.challenge_bytes as i32 + 1);
        }
    }

    Some((
        PublicKey {
            curve: pk_curve,
            hint_pk,
        },
        SecretKey {
            curve,
            secret_ideal: *secret_ideal,
            mat_bacan_to_ba0_two: *mat,
        },
    ))
}

/// `protocols_sign`. `None` where the reference returns 0 (the negligible
/// non-coprimality events, a failed chain) or on an entropy failure.
pub fn sign<L: FpBackend + PrimePrecomp, const N: usize>(
    params: &Params<N>,
    pk: &PublicKey<L>,
    sk: &SecretKey<L, N>,
    msg: &[u8],
    entropy: &mut impl Rng,
) -> Option<Signature<L>> {
    let alg = &params.alg;
    let act = &params.act;
    let vp = &params.verify;
    let f = L::TWO_ADIC_EXPONENT;
    let reduced_order = vp.response_bits + HD_EXTRA_TORSION;
    let seed = prng_seed(entropy)?;
    let mut def = domain(&seed, b"def");

    // the commitment: a prime-norm ideal, its small equivalent coprime to
    // two, and its isogeny with the images of E0's basis
    let (mut e_com, b_com, ideal_commit) = {
        let mut com = SplitRng {
            domain: domain(&seed, b"com"),
            default: &mut def,
        };
        let ideal = Zeroizing::new(QuatIdeal::random_given_prime_norm(
            &params.sec_degree,
            alg,
            &mut com,
        )?);
        let (beta, ideal) = ideal.small_equivalent_coprime(Some(&Ibz::two()), alg, &mut com)?;
        drop(Zeroizing::new(beta));
        let ideal = Zeroizing::new(ideal);
        let (e_com, b_com) = arbitrary_isogeny_evaluation::<L, N>(&ideal, alg, act, &mut com)?;
        (e_com, Zeroizing::new(b_com), ideal)
    };

    // the challenge: a scalar, the kernel P + [s] Q in the canonical basis
    // of E_pk, pulled back to E0 through the secret key's basis change
    let chall_coeff = hash_to_challenge(vp, pk, &e_com, msg);
    let vec = Vec2([
        Ibz::set(1, 2),
        Ibz::from_digits(&chall_coeff[..L::ORDER_WORDS]),
    ]);
    let vec = Zeroizing::new(sk.mat_bacan_to_ba0_two.eval(&vec));
    let (ideal_chall_two, chall_split) =
        kernel_dlogs_to_ideal_even(&vec, vp.challenge_bits, alg, act)?;
    let (ideal_chall_two, chall_split) =
        (Zeroizing::new(ideal_chall_two), Zeroizing::new(chall_split));

    // the response and the auxiliary ideal
    let (ideal_skchall, resp_quat, deg_odd, aux_ideal, mut aux_split) = {
        let mut res = SplitRng {
            domain: domain(&seed, b"res"),
            default: &mut def,
        };
        let (ideal_skchall, resp_quat, deg_odd, sk_chall_quat) = response_element(
            &sk.secret_ideal,
            &ideal_chall_two,
            Some(&*chall_split),
            &ideal_commit,
            vp.response_bits - 1,
            alg,
            &mut res,
        )?;
        let (ideal_skchall, resp_quat, deg_odd, sk_chall_quat) = (
            Zeroizing::new(ideal_skchall),
            Zeroizing::new(resp_quat),
            Zeroizing::new(deg_odd),
            Zeroizing::new(sk_chall_quat),
        );
        if !ideal_skchall.norm.gcd(&ideal_commit.norm).is_one() {
            // "Non-coprime resp norms. This should never happen."
            return None;
        }
        let remain = Ibz::<N>::one().mul_2exp(vp.response_bits);
        let mut random_aux_norm = Zeroizing::new(remain.sub(&deg_odd));
        random_aux_norm.set_bound(vp.response_bits as i32 + 1);
        let (aux_ideal, aux_split) = QuatIdeal::random_given_arbitrary_odd_norm(
            &random_aux_norm,
            &sk_chall_quat,
            alg,
            &mut res,
        )?;
        (
            ideal_skchall,
            resp_quat,
            deg_odd,
            Zeroizing::new(aux_ideal),
            Zeroizing::new(aux_split),
        )
    };
    if !aux_ideal.norm.gcd(&ideal_skchall.norm).is_one()
        || !aux_ideal.norm.gcd(&sk.secret_ideal.norm).is_one()
    {
        // "Non-coprime aux norms. This should never happen."
        return None;
    }
    let ideal_skchall_aux = Zeroizing::new(QuatIdeal::intersect_o0(
        &aux_ideal,
        &ideal_skchall,
        Some(&*aux_split),
    ));

    // 2^(RESPONSE_BITS + HD_EXTRA_TORSION): the torsion above the kernel
    let remain = Ibz::<N>::one().mul_2exp(vp.response_bits + HD_EXTRA_TORSION);

    // the isogeny of the intersection, and the response applied to the
    // images of E0's basis
    let out = Zeroizing::new(ideal_to_isogeny_qlapoty::<L, N>(
        &ideal_skchall_aux,
        alg,
        act,
        &mut def,
        &mut Id2IsoStats::default(),
    )?);
    let mut e_aux = out.codomain.clone();
    let mut b_aux = Zeroizing::new(out.basis.clone());
    let degree_resp_inv = Zeroizing::new(ideal_skchall.norm.invmod(&remain)?);
    for i in 0..2 {
        aux_split.coord.0[i] = aux_split.coord.0[i].mul(&degree_resp_inv).modulo(&remain);
    }
    for i in 2..4 {
        aux_split.coord.0[i].set_bound(remain.get_bound());
    }
    let resp_aux_quat = Zeroizing::new(aux_split.mul(&resp_quat, alg));
    endomorphism_application_even_basis(&mut b_aux, &e_aux, &resp_aux_quat, f, false, alg, act)?;

    // reduce both bases to the relevant order
    let b_com = Zeroizing::new(ec_dbl_iter_basis(
        &b_com,
        (f - reduced_order) as usize,
        &mut e_com,
    ));
    let b_aux = Zeroizing::new(ec_dbl_iter_basis(
        &b_aux,
        (f - reduced_order) as usize,
        &mut e_aux,
    ));

    // the (2^RESPONSE_BITS, 2^RESPONSE_BITS)-isogeny E_com x E_aux ->
    // E_aux2 x E_chall2 with kernel <(B_com.P, [1/deg] B_aux.P), (B_com.Q,
    // [1/deg] B_aux.Q)>, pushing B_com
    let degree_resp_inv = Zeroizing::new(deg_odd.invmod(&remain)?);
    let mut e12 = ThetaCoupleCurve {
        e1: e_com,
        e2: e_aux,
    };
    let mut ker = ThetaKernelCouplePoints::from_bases(&b_com, &b_aux);
    let scalar = Zeroizing::new(digits(&degree_resp_inv));
    let scalar = &scalar[..L::ORDER_WORDS];
    ker.t1.p2 = ec_mul(&ker.t1.p2, scalar, reduced_order as usize, &mut e12.e2);
    ker.t2.p2 = ec_mul(&ker.t2.p2, scalar, reduced_order as usize, &mut e12.e2);
    let push = [
        ThetaCouplePoint::on_e1(b_com.p.clone()),
        ThetaCouplePoint::on_e1(b_com.q.clone()),
        ThetaCouplePoint::on_e1(b_com.pmq.clone()),
    ];
    let out = theta_chain_compute_and_eval(
        vp.response_bits as u16,
        &mut e12,
        &ker,
        &push,
        ChainMode::Both,
    )?;
    let img = |i: usize| out.images[i].as_ref().expect("pushed point");
    // the auxiliary curve is the second factor, the challenge curve the first
    let mut e_aux_2 = out.codomain.e2.clone();
    let mut e_chall = out.codomain.e1.clone();
    let b_aux_2 = EcBasis::new(img(0).p2.clone(), img(1).p2.clone(), img(2).p2.clone());
    let mut b_chall_2 = EcBasis::new(img(0).p1.clone(), img(1).p1.clone(), img(2).p1.clone());

    // set_aux_curve_signature
    e_aux_2.normalize();
    let e_aux_a = e_aux_2.a.clone();

    // compute_and_set_basis_change_matrix
    let (b_can_chall, hint_chall) = ec_curve_to_basis_2f_to_hint(&mut e_chall, f, reduced_order)?;
    let (b_aux_2_can, hint_aux) = ec_curve_to_basis_2f_to_hint(&mut e_aux_2, f, reduced_order)?;
    let mut mat_baux2_to_baux2_can: Mat2x2<N> =
        change_of_basis_matrix_tate_invert(&b_aux_2_can, &b_aux_2, &mut e_aux_2, reduced_order)?;
    matrix_application_even_basis(
        &mut b_chall_2,
        &e_chall,
        &mut mat_baux2_to_baux2_can,
        reduced_order,
        true,
    )?;
    let mut mat_bchall_can_to_bchall: Mat2x2<N> =
        change_of_basis_matrix_tate(&b_chall_2, &b_can_chall, &mut e_chall, reduced_order)?;
    mat_bchall_can_to_bchall.normalize(reduced_order);
    let mut mat = [[[0u64; MAX_ORDER_WORDS]; 2]; 2];
    for (row, src) in mat.iter_mut().zip(mat_bchall_can_to_bchall.0.iter()) {
        for (x, s) in row.iter_mut().zip(src.iter()) {
            debug_assert!(s.bitsize() as u32 <= reduced_order);
            *x = digits(s);
        }
    }

    Some(Signature {
        e_aux_a,
        mat_bchall_can_to_bchall: mat,
        chall_coeff,
        hint_aux,
        hint_chall,
    })
}

/// `sig_ibz_from_bytes`: `ceil(nbytes / 8)` digits copied with
/// `ibz_copy_digits`, so the bound is `64 * ndigits + 1`.
fn ibz_from_bytes<const N: usize>(bytes: &[u8]) -> Ibz<N> {
    let d = sqisign_verify::sqisign::digits_from_bytes(bytes);
    Ibz::from_digits(&d[..bytes.len().div_ceil(8)])
}

fn ibz_to_bytes<const N: usize>(x: &Ibz<N>, out: &mut [u8]) {
    assert!(x.is_positive() && x.bitsize() <= 8 * out.len() as i32);
    x.to_le_bytes(out);
}

/// `secret_key_to_bytes`: the public key, the inert form of the secret
/// ideal, the basis-change matrix.
pub fn secret_key_to_bytes<L: FpBackend, const N: usize>(
    params: &Params<N>,
    sk: &SecretKey<L, N>,
    pk: &PublicKey<L>,
) -> Zeroizing<Vec<u8>> {
    let vp = &params.verify;
    let mut pkb = [0u8; MAX_PUBLICKEY_BYTES];
    let n = public_key_to_bytes(pk, &mut pkb);
    // one allocation of the final size: a reallocation would leave an
    // unscrubbed copy behind
    let mut out = Zeroizing::new(Vec::with_capacity(vp.secretkey_bytes));
    out.extend_from_slice(&pkb[..n]);
    let fp = vp.fp_encoded_bytes;
    let cb = vp.challenge_bytes;
    let mut push = |x: &Ibz<N>, n: usize| {
        let start = out.len();
        out.resize(start + n, 0);
        ibz_to_bytes(x, &mut out[start..]);
    };
    push(&sk.secret_ideal.norm, fp);
    push(&sk.secret_ideal.x, fp);
    push(&sk.secret_ideal.y, fp);
    for i in 0..2 {
        for j in 0..2 {
            push(&sk.mat_bacan_to_ba0_two.0[i][j], cb);
        }
    }
    assert_eq!(out.len(), vp.secretkey_bytes);
    out
}

/// `secret_key_from_bytes`, also returning the embedded public key.
pub fn secret_key_from_bytes<L: FpBackend, const N: usize>(
    params: &Params<N>,
    bytes: &[u8],
) -> Option<(SecretKey<L, N>, PublicKey<L>)> {
    let vp = &params.verify;
    if bytes.len() != vp.secretkey_bytes {
        return None;
    }
    let pk = public_key_from_bytes::<L>(vp, &bytes[..vp.publickey_bytes])?;
    let fp = vp.fp_encoded_bytes;
    let cb = vp.challenge_bytes;
    let mut pos = vp.publickey_bytes;
    let mut next = |n: usize| {
        let x = ibz_from_bytes::<N>(&bytes[pos..pos + n]);
        pos += n;
        x
    };
    let norm = next(fp);
    let x = next(fp);
    let y = next(fp);
    let mut mat = Mat2x2::<N>::zero();
    for i in 0..2 {
        for j in 0..2 {
            mat.0[i][j] = next(cb);
        }
    }
    Some((
        SecretKey {
            curve: pk.curve.clone(),
            secret_ideal: QuatIdeal { x, y, norm },
            mat_bacan_to_ba0_two: mat,
        },
        pk,
    ))
}
