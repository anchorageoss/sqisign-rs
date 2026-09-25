//! The compact (dimension-4) format's key generation and signing, level I,
//! experimental (`docs/COMPACT_R3.md`; the verifier is
//! `sqisign_verify::hd`).
//!
//! Round 3 up to the response: the secret ideal and the commitment are the
//! round-3 signer's, the challenge is the round-3 hash shape over the
//! compact public key (`SHAKE256("SQI" || pk || A_com || msg)`, `λ` bits)
//! with kernel `P + [c] Q` on the library's hinted basis of `E_pk`. Then
//! SQIsignHD's response: a quaternion `γ` of the connecting lattice with
//! degree `q < 2^e` such that `2^e − q` is a prime `≡ 1 (mod 4)`
//! (`q ≡ 3 mod 4`, sampled and rejected), and its action on the
//! `2^r`-torsion of `E_com` as three scalars modulo `2^r` plus `q`.
//!
//! # The action matrix
//!
//! Write `J` for the commitment ideal (`ψ: E_0 → E_com`), `K = I_sk ∩
//! I_chall` (`φ_chall ∘ τ: E_0 → E_chall`), and let `response_element` return
//! `K' ∼ K` of small odd norm with connecting element `β̄` (`β = φ̂_K' ∘ φ_K`)
//! and `γ' ∈ conj(K') J` of degree `q` (`γ' = φ̂_K' ∘ σ ∘ ψ`). Then `ε = β̄ γ'
//! / n(K') = φ̂_K ∘ σ ∘ ψ = τ̂ ∘ φ̂_chall ∘ σ ∘ ψ`, an element of `O_0`, so for
//! a point `X` of `E_0[2^f]`:
//!
//! ```text
//! coords_{B_pk_can}( φ̂_chall(σ(ψ(X))) ) = n(sk)^-1 · M_sk^-1 · A_ε · coords_{B_0}(X)
//! ```
//!
//! with `A_ε` the action matrix of `ε` on `B_0`, `M_sk` the change of basis
//! from the hinted basis of `E_pk` to `τ(B_0)`. With `X` running over the
//! hinted basis of `E_com` written in `ψ(B_0)` (`M_com`), `M = M_sk^-1 A_ε
//! M_com n(sk)^-1 (mod 2^f)` holds `φ̂_chall(σ(P_com))` and
//! `φ̂_chall(σ(Q_com))` in the public key's hinted basis. Applying `φ_chall`
//! and passing to the `2^r`-torsion, with `P_resc = φ(P) + [c] φ(Q)` and
//! `Q_resc = [2^λ] φ(Q)` as the verifier builds them, gives the signature
//! scalars `a = M_00`, `b = (M_10 − c M_00) / 2^λ`, `c_or_d = M_01` when `a`
//! is odd, else `(M_11 − c M_01) / 2^λ`, all modulo `2^r`. This is the
//! round-2 compact signer's formula (the SQIsignHD library's) with round 3's
//! reduced connecting ideal accounted for by `β̄` and `n(K')`.

use sqisign_verify::ec::jacobian::jac_add;
use sqisign_verify::ec::{EcBasis, EcCurve};
use sqisign_verify::fp::Fp2;
use sqisign_verify::hd::{
    canonical_hints, encode_public_key, encode_signature, hd_challenge, hd_challenge_len,
    hd_torsion_basis, pk_wire_bytes, sig_wire_bytes, HdLevel, MAX_CHAL_BYTES, MAX_PK_WIRE_BYTES,
    U4,
};

use crate::id2iso::{
    arbitrary_isogeny_evaluation, change_of_basis_matrix_tate, endomorphism_action_matrix,
    kernel_dlogs_to_ideal_even,
};
use crate::mp::{DefaultDomain, Ibz, Rng, ShakeRng};
use crate::quat::{Mat2x2, QuatIdeal, Vec2};
use crate::sqisign::Params;

/// A compact public key: the curve and the library's two basis hints.
#[derive(Clone, Debug)]
pub struct CompactPublicKey<L: HdLevel> {
    /// Montgomery coefficient of `E_pk`.
    pub a_pk: Fp2<L>,
    /// Hint of `P`.
    pub hint_p: u32,
    /// Hint of `Q`.
    pub hint_q: u32,
}

/// A compact secret key.
#[derive(Clone)]
pub struct CompactSecretKey<L: HdLevel, const N: usize> {
    /// `E_pk`, for the size-checked encoding of the key (unused by signing).
    #[allow(dead_code)]
    curve: EcCurve<L>,
    secret_ideal: QuatIdeal<N>,
    /// `M_sk`: coordinates in `τ(B_0)` of a point given in the hinted basis
    /// of `E_pk[2^f]`, modulo `2^f`.
    mat_bpkcan_to_bpk0: Mat2x2<N>,
}

impl<L: HdLevel, const N: usize> zeroize::Zeroize for CompactSecretKey<L, N> {
    fn zeroize(&mut self) {
        self.secret_ideal.x.zeroize();
        self.secret_ideal.y.zeroize();
        self.secret_ideal.norm.zeroize();
        for row in self.mat_bpkcan_to_bpk0.0.iter_mut() {
            for x in row.iter_mut() {
                x.zeroize();
            }
        }
    }
}

impl<L: HdLevel, const N: usize> Drop for CompactSecretKey<L, N> {
    fn drop(&mut self) {
        zeroize::Zeroize::zeroize(self);
    }
}

/// Why signing failed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CompactSignError {
    /// The entropy source failed.
    Entropy,
    /// No hinted basis for the commitment curve (a crafted or degenerate curve).
    CommitmentBasis,
    /// The challenge kernel gave no ideal.
    Challenge,
    /// No response of good degree within the sampling budget.
    NoGoodResponse,
    /// A matrix or norm was not invertible modulo `2^f`.
    NotInvertible,
    /// The challenge division of the first column was not exact.
    NotDivisibleB,
    /// The challenge division of the second column was not exact.
    NotDivisibleD,
    /// The signature did not encode.
    Encode,
}

/// How many response samples to try before giving up on a signature.
const MAX_RESPONSE_TRIES: usize = 4096;
/// Miller–Rabin rounds on `2^e − q`.
const PRIME_ROUNDS: u32 = 40;

fn seeded(entropy: &mut impl Rng, dom: &[u8; 3]) -> Option<ShakeRng> {
    let mut seed = [0u8; 48];
    entropy.fill(&mut seed).then_some(())?;
    Some(ShakeRng::new(&seed, dom))
}

/// The hinted basis of `E[2^f]` as an x-only basis, with its hints.
fn hinted_basis<L: HdLevel>(curve: &EcCurve<L>) -> Option<(EcBasis<L>, u32, u32)> {
    let (hp, hq) = canonical_hints::<L>(&curve.a)?;
    let (p, q) = hd_torsion_basis::<L>(&curve.a, hp, hq)?;
    let pmq = jac_add(&p, &q.neg(), curve);
    Some((EcBasis::new(p.to_xz(), q.to_xz(), pmq.to_xz()), hp, hq))
}

/// Key generation: the round-3 secret ideal and curve, the public key's
/// hinted basis and its change of basis to the `E_0` images at the full
/// torsion. `None` on an entropy failure.
pub fn compact_keygen<L: HdLevel, const N: usize>(
    params: &Params<N>,
    entropy: &mut impl Rng,
) -> Option<(CompactPublicKey<L>, CompactSecretKey<L, N>)> {
    let alg = &params.alg;
    let act = &params.act;
    let f = L::TWO_ADIC_EXPONENT;
    let mut rng = seeded(entropy, b"ckg")?;
    loop {
        let Some(ideal) = QuatIdeal::random_given_prime_norm(&params.sec_degree, alg, &mut rng)
        else {
            continue;
        };
        let Some((_, ideal)) = ideal.small_equivalent_coprime(Some(&Ibz::zero()), alg, &mut rng)
        else {
            continue;
        };
        let Some((mut curve, b_0_two)) =
            arbitrary_isogeny_evaluation::<L, N>(&ideal, alg, act, &mut rng)
        else {
            continue;
        };
        curve.normalize();
        let Some((basis, hp, hq)) = hinted_basis::<L>(&curve) else {
            continue;
        };
        let Some(mat) = change_of_basis_matrix_tate::<L, N>(&basis, &b_0_two, &mut curve, f) else {
            continue;
        };
        let mut pk_curve = curve.clone();
        pk_curve.is_a24_computed_and_normalized = false;
        return Some((
            CompactPublicKey {
                a_pk: pk_curve.a.clone(),
                hint_p: hp,
                hint_q: hq,
            },
            CompactSecretKey {
                curve,
                secret_ideal: ideal,
                mat_bpkcan_to_bpk0: mat,
            },
        ));
    }
}

/// Encode a compact public key; returns the length.
pub fn compact_public_key_to_bytes<L: HdLevel>(
    pk: &CompactPublicKey<L>,
    out: &mut [u8],
) -> Option<usize> {
    encode_public_key::<L>(&pk.a_pk, pk.hint_p, pk.hint_q, out)
}

fn to_u128<const N: usize>(x: &Ibz<N>) -> u128 {
    let mut d = [0u64; 2];
    let mut all = [0u64; 32];
    x.to_digits(&mut all[..N.min(32)]);
    d.copy_from_slice(&all[..2]);
    (d[0] as u128) | ((d[1] as u128) << 64)
}

fn to_u4<const N: usize>(x: &Ibz<N>) -> U4 {
    let mut all = [0u64; 32];
    x.to_digits(&mut all[..N.min(32)]);
    U4([all[0], all[1], all[2], all[3]])
}

/// `x / 2^k` when `2^k` divides the non-negative `x`, else `None`.
fn exact_div_2exp<const N: usize>(x: &Ibz<N>, k: u32) -> Option<Ibz<N>> {
    if !x.mod2exp(k).is_zero() {
        return None;
    }
    Some(x.div_2exp(k))
}

/// Sign `msg`; writes the signature into `out` (at least
/// [`sig_wire_bytes`] long) and returns its length. `None` on an entropy
/// failure or when no good response was found within the sampling budget.
pub fn compact_sign<L: HdLevel, const N: usize>(
    params: &Params<N>,
    pk: &CompactPublicKey<L>,
    sk: &CompactSecretKey<L, N>,
    msg: &[u8],
    entropy: &mut impl Rng,
    out: &mut [u8],
) -> Result<usize, CompactSignError> {
    use CompactSignError as E;
    let alg = &params.alg;
    let act = &params.act;
    let f = L::TWO_ADIC_EXPONENT;
    let two_f = Ibz::<N>::one().mul_2exp(f);
    let two_r = Ibz::<N>::one().mul_2exp(L::R);
    let two_e = Ibz::<N>::one().mul_2exp(L::E_EMBED);
    let four = Ibz::<N>::set(4, 4);
    let mut rng = seeded(entropy, b"csg").ok_or(E::Entropy)?;

    // commitment (round 3)
    let (mut e_com, b_com, ideal_commit) = loop {
        let Some(ideal) = QuatIdeal::random_given_prime_norm(&params.sec_degree, alg, &mut rng)
        else {
            continue;
        };
        let Some((_, ideal)) = ideal.small_equivalent_coprime(Some(&Ibz::two()), alg, &mut rng)
        else {
            continue;
        };
        if let Some((e, b)) = arbitrary_isogeny_evaluation::<L, N>(&ideal, alg, act, &mut rng) {
            break (e, b, ideal);
        }
    };
    e_com.normalize();
    let a_com = e_com.a.clone();
    let (b_com_can, hp_com, hq_com) = hinted_basis::<L>(&e_com).ok_or(E::CommitmentBasis)?;
    let m_com: Mat2x2<N> = change_of_basis_matrix_tate::<L, N>(&b_com_can, &b_com, &mut e_com, f)
        .ok_or(E::CommitmentBasis)?;

    // challenge: SHAKE256("SQI" || pk || A_com || msg), λ bits
    let mut pkb = [0u8; MAX_PK_WIRE_BYTES];
    let n = compact_public_key_to_bytes(pk, &mut pkb).ok_or(E::Encode)?;
    let mut chal = [0u8; MAX_CHAL_BYTES];
    let clen = hd_challenge_len::<L>();
    hd_challenge::<L>(&pkb[..n], &a_com, msg, &mut chal[..clen]);
    let mut k = Ibz::<N>::from_le_bytes(&chal[..clen]);
    k.set_bound(L::LAMBDA as i32 + 1);

    // the challenge ideal: kernel P + [k] Q on the hinted basis, in τ(B_0) coordinates
    let mut vec = sk.mat_bpkcan_to_bpk0.eval(&Vec2([Ibz::one(), k]));
    for x in vec.0.iter_mut() {
        *x = x.mod2exp(L::LAMBDA);
    }
    let (ideal_chall_two, chall_split) =
        kernel_dlogs_to_ideal_even(&vec, L::LAMBDA, alg, act).ok_or(E::Challenge)?;

    // the response: response_element's construction once (the intersection
    // with the secret ideal, its small equivalent K' with the connecting
    // element β̄, the product lattice conj(K') J), then samples from the
    // ball until the degree is good: q odd, ≡ 3 mod 4, 2^e − q a probable
    // prime. (`response_element` itself would redo the lattice reduction on
    // every sample.)
    let inter = QuatIdeal::intersect_o0(&ideal_chall_two, &sk.secret_ideal, Some(&chall_split));
    let (beta_bar, k_ideal) = {
        let (beta, k_ideal) = inter
            .small_equivalent_coprime(Some(&Ibz::zero()), alg, &mut DefaultDomain(&mut rng))
            .ok_or(E::NoGoodResponse)?;
        (beta.conj(), k_ideal)
    };
    if !k_ideal.norm.gcd(&ideal_commit.norm).is_one() {
        return Err(E::NoGoodResponse);
    }
    let prod = QuatIdeal::mul_o0(&k_ideal, &ideal_commit);
    let mut bound = Ibz::<N>::one().mul_2exp(L::E_EMBED).sub(&Ibz::one());
    bound.set_bound(L::E_EMBED as i32 + 1);
    let mut found = None;
    for _ in 0..MAX_RESPONSE_TRIES {
        let Some((gamma, mut q)) = prod.sample_from_ball(&bound, alg, &mut rng) else {
            continue;
        };
        q.set_bound(L::E_EMBED as i32 + 1);
        if q.is_zero() || !q.is_odd() || q.modulo(&four).cmp_i32(3) != core::cmp::Ordering::Equal {
            continue;
        }
        if q >= two_e {
            continue;
        }
        let mut nn = two_e.sub(&q);
        nn.set_bound(L::E_EMBED as i32 + 1);
        if !nn.probab_prime(PRIME_ROUNDS, &mut rng) {
            continue;
        }
        found = Some((gamma, q));
        break;
    }
    let (gamma, q) = found.ok_or(E::NoGoodResponse)?;

    // ε = β̄ γ' (the n(K') factor sits in the content) and its action matrix
    // The challenge ideal is `I_c · s` with `s` the split element (`1`, or
    // `1 − i` when the kernel is above `(0, 0)` on `E_0`), so its isogeny is
    // `φ_{I_c} ∘ s` and the dual picks up `conj(s)` after `ε`.
    let eps = chall_split.conj().mul(&beta_bar.mul(&gamma, alg), alg);
    let mut a_eps = endomorphism_action_matrix(&eps, f, alg, act);
    let nk_inv = k_ideal.norm.invmod(&two_f).ok_or(E::NotInvertible)?;
    for row in a_eps.0.iter_mut() {
        for x in row.iter_mut() {
            *x = x.mul(&nk_inv).modulo(&two_f);
        }
    }
    let (inv_sk, ok) = sk.mat_bpkcan_to_bpk0.inv_mod(&two_f);
    if !ok {
        return Err(E::NotInvertible);
    }
    let nsk_inv = sk
        .secret_ideal
        .norm
        .invmod(&two_f)
        .ok_or(E::NotInvertible)?;
    let mut m = inv_sk.mul_mod(&a_eps, &two_f).mul_mod(&m_com, &two_f);
    for row in m.0.iter_mut() {
        for x in row.iter_mut() {
            *x = x.mul(&nsk_inv).modulo(&two_f);
        }
    }

    // the scalars: through the challenge, onto the 2^r-torsion
    let a = m.0[0][0].modulo(&two_r);
    let v_b = m.0[1][0].sub(&m.0[0][0].mul(&k)).modulo(&two_f);
    let b = exact_div_2exp(&v_b, L::LAMBDA)
        .ok_or(E::NotDivisibleB)?
        .modulo(&two_r);
    let c_or_d = if a.is_odd() {
        m.0[0][1].modulo(&two_r)
    } else {
        let v_d = m.0[1][1].sub(&m.0[0][1].mul(&k)).modulo(&two_f);
        exact_div_2exp(&v_d, L::LAMBDA)
            .ok_or(E::NotDivisibleD)?
            .modulo(&two_r)
    };

    encode_signature::<L>(
        &a_com,
        to_u128(&a) as i128,
        to_u128(&b) as i128,
        to_u128(&c_or_d) as i128,
        &to_u4(&q),
        hp_com,
        hq_com,
        out,
    )
    .ok_or(E::Encode)
}

/// The wire size of a compact signature at this level.
pub const fn compact_signature_bytes<L: HdLevel>() -> usize {
    sig_wire_bytes::<L>()
}

/// The wire size of a compact public key at this level.
pub const fn compact_public_key_bytes<L: HdLevel>() -> usize {
    pk_wire_bytes::<L>()
}

// ---- the typed API ----------------------------------------------------------

use sqisign_verify::compact::{
    CompactLevel, CompactPublicKey as TypedCompactPublicKey, CompactSignature,
};

/// A compact signing key (level I), the typed API over [`compact_sign`].
pub struct CompactSigningKey<L: CompactLevel, const N: usize> {
    sk: CompactSecretKey<L, N>,
    pk: CompactPublicKey<L>,
    params: Params<N>,
}

impl<L: CompactLevel, const N: usize> CompactSigningKey<L, N> {
    /// Sign `msg` with randomness from `rng`.
    pub fn sign(
        &self,
        msg: &[u8],
        rng: &mut impl Rng,
    ) -> Result<CompactSignature<L>, CompactSignError> {
        let mut out = [0u8; sqisign_verify::hd::MAX_SIG_WIRE_BYTES];
        let n = compact_sign(&self.params, &self.pk, &self.sk, msg, rng, &mut out)?;
        CompactSignature::from_bytes(&out[..n]).map_err(|_| CompactSignError::Encode)
    }

    /// The public key.
    pub fn public_key(&self) -> Result<TypedCompactPublicKey<L>, CompactSignError> {
        let mut pkb = [0u8; MAX_PK_WIRE_BYTES];
        let n = compact_public_key_to_bytes(&self.pk, &mut pkb).ok_or(CompactSignError::Encode)?;
        TypedCompactPublicKey::from_bytes(&pkb[..n]).map_err(|_| CompactSignError::Encode)
    }
}

impl<L: CompactLevel, const N: usize> core::fmt::Debug for CompactSigningKey<L, N> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("CompactSigningKey(..)")
    }
}

/// Generate a compact key pair at level I.
pub fn generate_compact<L: CompactLevel, const N: usize>(
    params: Params<N>,
    rng: &mut impl Rng,
) -> Option<(TypedCompactPublicKey<L>, CompactSigningKey<L, N>)> {
    let (pk, sk) = compact_keygen::<L, N>(&params, rng)?;
    let key = CompactSigningKey { sk, pk, params };
    let typed = key.public_key().ok()?;
    Some((typed, key))
}
