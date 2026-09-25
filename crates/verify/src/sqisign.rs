//! Round-3 verification and the public encodings (`verify.c`, `common.c`,
//! `encode_public.c`), statement for statement, on the round-3 layers of
//! this module. The signing side (`keygen.c`, `sign.c`, `encode_secret.c`)
//! is `sqisign_rs::sqisign`.

use crate::ec::basis::ec_curve_to_basis_2f_from_hint;
use crate::ec::isogeny::iso_isogeny_2chain;
use crate::ec::point::{ec_biscalar_mul_verif, ec_ladder3pt};
use crate::ec::{EcBasis, EcCurve, EcIsogEven, MAX_ORDER_WORDS};
use crate::fp::{Fp2, FpBackend};
use crate::params::Prime;
use crate::precomp::{PrimePrecomp, EXTRA_TORSION};
use crate::theta::chain::theta_chain_compute_and_eval;
use crate::theta::{ChainMode, ThetaCoupleCurve, ThetaKernelCouplePoints};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

/// `scalar_t`: `NWORDS_ORDER` little-endian digits (zero-padded here).
pub type Scalar = [u64; MAX_ORDER_WORDS];

/// `HD_EXTRA_TORSION`: the chains of `sign` and `verify` use two more bits
/// of torsion than their length.
pub const HD_EXTRA_TORSION: u32 = 2;

/// The largest encoded public key of the three levels (`p664_17`: 169).
pub const MAX_PUBLICKEY_BYTES: usize = 169;
/// The largest encoded signature of the three levels (`p664_17`: 406).
pub const MAX_SIGNATURE_BYTES: usize = 406;
/// The largest challenge scalar (`p664_17`: 32 bytes).
pub const MAX_CHALLENGE_BYTES: usize = 32;

/// The per-level constants the verifier and the encodings need: the
/// reference's `SQISIGN_*` sizes, no quaternion data.
#[derive(Clone, Debug)]
pub struct VerifyParams {
    /// `CHALLENGE_BITS = lambda`.
    pub challenge_bits: u32,
    /// `RESPONSE_BITS`.
    pub response_bits: u32,
    /// Bytes of a basis-change matrix entry in a signature.
    pub response_bytes: usize,
    /// Bytes of the challenge scalar.
    pub challenge_bytes: usize,
    /// `f`, the 2-adic exponent of `p + 1`.
    pub torsion_even_power: u32,
    /// Bytes of a canonical `Fp` element.
    pub fp_encoded_bytes: usize,
    /// Encoded sizes.
    pub publickey_bytes: usize,
    /// Encoded sizes.
    pub secretkey_bytes: usize,
    /// Encoded sizes.
    pub signature_bytes: usize,
}

macro_rules! level {
    ($(#[$m:meta])* $name:ident, $modname:ident) => {
        $(#[$m])*
        pub const fn $name() -> VerifyParams {
            use crate::params::sqisign_v3::$modname as l;
            VerifyParams {
                challenge_bits: l::LAMBDA,
                response_bits: l::RESPONSE_BITS,
                response_bytes: l::RESPONSE_BYTES,
                challenge_bytes: l::CHALLENGE_BYTES,
                torsion_even_power: l::TORSION_EVEN_POWER,
                fp_encoded_bytes: (l::PUBLICKEY_BYTES - 1) / 2,
                publickey_bytes: l::PUBLICKEY_BYTES,
                secretkey_bytes: l::SECRETKEY_BYTES,
                signature_bytes: l::SIGNATURE_BYTES,
            }
        }
    };
}

level!(
    /// NIST level I, `p324_3` ([`crate::P324_3`]).
    level1, p324_3
);
level!(
    /// NIST level III, `p500_27` ([`crate::P500_27`]).
    level3, p500_27
);
level!(
    /// NIST level V, `p664_17` ([`crate::P664_17`]).
    level5, p664_17
);

/// `public_key_t`.
#[derive(Clone, Debug)]
pub struct PublicKey<L: Prime> {
    /// `E_pk`, with `C = 1` once encoded.
    pub curve: EcCurve<L>,
    /// Hint for the canonical basis of `E_pk[2^(CHALLENGE_BITS + 2)]`.
    pub hint_pk: u8,
}

/// `signature_t`.
#[derive(Clone, Debug)]
pub struct Signature<L: Prime> {
    /// Montgomery coefficient of the auxiliary curve.
    pub e_aux_a: Fp2<L>,
    /// Change of basis on the challenge curve, entries in
    /// `[0, 2^(RESPONSE_BITS + 2))`, normalised.
    pub mat_bchall_can_to_bchall: [[Scalar; 2]; 2],
    /// The challenge scalar.
    pub chall_coeff: Scalar,
    /// Basis hint of the auxiliary curve.
    pub hint_aux: u8,
    /// Basis hint of the challenge curve.
    pub hint_chall: u8,
}

/// Little-endian bytes of a scalar's low `out.len()` bytes.
pub fn digits_to_bytes(x: &Scalar, out: &mut [u8]) {
    for (i, b) in out.iter_mut().enumerate() {
        *b = (x[i / 8] >> (8 * (i % 8))) as u8;
    }
}

/// A scalar from little-endian bytes (at most `8 * MAX_ORDER_WORDS`).
pub fn digits_from_bytes(bytes: &[u8]) -> Scalar {
    let mut x = [0u64; MAX_ORDER_WORDS];
    for (i, b) in bytes.iter().enumerate() {
        x[i / 8] |= (*b as u64) << (8 * (i % 8));
    }
    x
}

/// `public_key_to_bytes`: `A / C` then the hint, into `out` (at least
/// `2 * fp_encoded_bytes + 1` long); returns the length written.
pub fn public_key_to_bytes<L: FpBackend>(pk: &PublicKey<L>, out: &mut [u8]) -> usize {
    let a = pk.curve.a.mul(&pk.curve.c.inv());
    let enc = a.encode();
    let n = enc.len();
    out[..n].copy_from_slice(&enc[..]);
    out[n] = pk.hint_pk;
    n + 1
}

/// `public_key_from_bytes`. `None` if the length or a field element is
/// not canonical.
pub fn public_key_from_bytes<L: FpBackend>(
    params: &VerifyParams,
    bytes: &[u8],
) -> Option<PublicKey<L>> {
    if bytes.len() != params.publickey_bytes {
        return None;
    }
    let n = 2 * params.fp_encoded_bytes;
    let a = Fp2::<L>::decode(&bytes[..n])?;
    Some(PublicKey {
        curve: EcCurve {
            a,
            ..EcCurve::default()
        },
        hint_pk: bytes[n],
    })
}

/// `signature_to_bytes`: `A_aux`, the four matrix entries, the challenge,
/// the two hints, into `out` (at least `signature_bytes` long); returns
/// the length written, `None` if `out` is too short.
pub fn signature_to_bytes<L: FpBackend>(
    params: &VerifyParams,
    sig: &Signature<L>,
    out: &mut [u8],
) -> Option<usize> {
    if out.len() < params.signature_bytes {
        return None;
    }
    let enc = sig.e_aux_a.encode();
    let mut pos = enc.len();
    out[..pos].copy_from_slice(&enc[..]);
    let rb = params.response_bytes;
    for row in sig.mat_bchall_can_to_bchall.iter() {
        for x in row.iter() {
            digits_to_bytes(x, &mut out[pos..pos + rb]);
            pos += rb;
        }
    }
    digits_to_bytes(
        &sig.chall_coeff,
        &mut out[pos..pos + params.challenge_bytes],
    );
    pos += params.challenge_bytes;
    out[pos] = sig.hint_aux;
    out[pos + 1] = sig.hint_chall;
    pos += 2;
    debug_assert_eq!(pos, params.signature_bytes);
    Some(pos)
}

/// `signature_from_bytes`. `None` if the length or `A_aux` is not canonical.
pub fn signature_from_bytes<L: FpBackend>(
    params: &VerifyParams,
    bytes: &[u8],
) -> Option<Signature<L>> {
    if bytes.len() != params.signature_bytes {
        return None;
    }
    let mut pos = 2 * params.fp_encoded_bytes;
    let e_aux_a = Fp2::<L>::decode(&bytes[..pos])?;
    let rb = params.response_bytes;
    let mut next = |n: usize| {
        let x = digits_from_bytes(&bytes[pos..pos + n]);
        pos += n;
        x
    };
    let mut mat = [[[0u64; MAX_ORDER_WORDS]; 2]; 2];
    for row in mat.iter_mut() {
        for x in row.iter_mut() {
            *x = next(rb);
        }
    }
    let chall_coeff = next(params.challenge_bytes);
    let hint_aux = bytes[pos];
    let hint_chall = bytes[pos + 1];
    Some(Signature {
        e_aux_a,
        mat_bchall_can_to_bchall: mat,
        chall_coeff,
        hint_aux,
        hint_chall,
    })
}

/// `hash_to_challenge`: `SHAKE256("SQI" || pk || A_com || msg)`, the first
/// `CHALLENGE_BYTES` as a little-endian scalar. `A_com` is taken as
/// stored, without normalisation, as in the reference.
pub fn hash_to_challenge<L: FpBackend>(
    params: &VerifyParams,
    pk: &PublicKey<L>,
    com_curve: &EcCurve<L>,
    msg: &[u8],
) -> Scalar {
    let mut h = Shake256::default();
    h.update(b"SQI");
    let mut pkb = [0u8; MAX_PUBLICKEY_BYTES];
    let n = public_key_to_bytes(pk, &mut pkb);
    h.update(&pkb[..n]);
    h.update(&com_curve.a.encode()[..]);
    h.update(msg);
    let mut bytes = [0u8; MAX_CHALLENGE_BYTES];
    h.finalize_xof().read(&mut bytes[..params.challenge_bytes]);
    digits_from_bytes(&bytes[..params.challenge_bytes])
}

/// Bits of a scalar.
pub(crate) fn bit_length(x: &Scalar) -> u32 {
    for i in (0..MAX_ORDER_WORDS).rev() {
        if x[i] != 0 {
            return 64 * i as u32 + 64 - x[i].leading_zeros();
        }
    }
    0
}

/// Whether `x == 2^k`.
pub(crate) fn is_power_of_two_at(x: &Scalar, k: u32) -> bool {
    bit_length(x) == k + 1 && x.iter().map(|w| w.count_ones()).sum::<u32>() == 1
}

/// `check_canonical_basis_change_matrix`: every entry below
/// `2^(RESPONSE_BITS + 2)` and the matrix normalised. Scanning the entries
/// in order with `mid = 2^(RESPONSE_BITS + 1)`, every entry before the
/// first one that is neither `0` nor `mid` must be at most `mid`, and that
/// first entry must be below `mid`.
pub(crate) fn check_canonical_basis_change_matrix<L: FpBackend>(
    params: &VerifyParams,
    sig: &Signature<L>,
) -> bool {
    let mid_bits = params.response_bits + HD_EXTRA_TORSION - 1; // mid = 2^mid_bits
    let mut normal = false;
    for row in sig.mat_bchall_can_to_bchall.iter() {
        for x in row.iter() {
            let bits = bit_length(x);
            // x >= bound = 2^(mid_bits + 1)
            if bits > mid_bits + 1 {
                return false;
            }
            if normal {
                continue;
            }
            // x > mid: the same bit length as mid and not mid itself
            if bits == mid_bits + 1 && !is_power_of_two_at(x, mid_bits) {
                return false;
            }
            // normal = x < mid && x != 0
            normal = bits <= mid_bits && bits != 0;
        }
    }
    true
}

/// `protocols_verify`.
pub fn verify<L: FpBackend + PrimePrecomp>(
    params: &VerifyParams,
    pk: &PublicKey<L>,
    sig: &Signature<L>,
    msg: &[u8],
) -> bool {
    let Some(key) = PreparedPublicKey::new(params, pk) else {
        return false;
    };
    verify_inner(params, &key, sig, msg).is_some()
}

/// A public key with the work every verification under it repeats done
/// once: `E_pk` normalised and its canonical basis of
/// `E_pk[2^(CHALLENGE_BITS + 2)]` rebuilt from the hint (about 3 % of a
/// verification).
#[derive(Clone, Debug)]
pub struct PreparedPublicKey<L: FpBackend> {
    pub(crate) pk: PublicKey<L>,
    pub(crate) curve: EcCurve<L>,
    pub(crate) basis: EcBasis<L>,
}

impl<L: FpBackend + PrimePrecomp> PreparedPublicKey<L> {
    /// Prepare `pk`; `None` on an invalid coefficient or hint.
    pub fn new(params: &VerifyParams, pk: &PublicKey<L>) -> Option<Self> {
        if !EcCurve::<L>::verify_a(&pk.curve.a) {
            return None;
        }
        debug_assert!(
            bool::from(pk.curve.c.ct_is_one()) && !pk.curve.is_a24_computed_and_normalized
        );
        let f_chall = params.challenge_bits + EXTRA_TORSION;
        let mut curve = pk.curve.clone();
        let basis = ec_curve_to_basis_2f_from_hint(&mut curve, f_chall, pk.hint_pk)?;
        Some(Self {
            pk: pk.clone(),
            curve,
            basis,
        })
    }

    /// The key.
    pub fn public_key(&self) -> &PublicKey<L> {
        &self.pk
    }
}

/// `protocols_verify` under a prepared key.
pub fn verify_prepared<L: FpBackend + PrimePrecomp>(
    params: &VerifyParams,
    key: &PreparedPublicKey<L>,
    sig: &Signature<L>,
    msg: &[u8],
) -> bool {
    verify_inner(params, key, sig, msg).is_some()
}

/// The verifier's state after the challenge: `E_chall` with its canonical
/// basis of `E_chall[2^(RESPONSE_BITS + 2)]`, and `E_aux` with its.
pub(crate) struct ChallengeState<L: Prime> {
    pub(crate) e_chall: EcCurve<L>,
    pub(crate) b_chall_can: EcBasis<L>,
    pub(crate) e_aux: EcCurve<L>,
    pub(crate) b_aux_can: EcBasis<L>,
}

/// `compute_challenge_verify` and the canonical bases of
/// `challenge_and_aux_basis_verify`: everything of a verification that
/// does not depend on the basis-change matrix.
pub(crate) fn challenge_state<L: FpBackend + PrimePrecomp>(
    params: &VerifyParams,
    key: &PreparedPublicKey<L>,
    e_aux_a: &Fp2<L>,
    chall_coeff: &Scalar,
    hint_chall: u8,
    hint_aux: u8,
) -> Option<ChallengeState<L>> {
    let mut e_aux = EcCurve::<L>::from_a(e_aux_a)?;

    // compute_challenge_verify: the 2^(CHALLENGE_BITS + 2)-isogeny with
    // kernel P + [chall_coeff] Q on the canonical basis of E_pk
    let f_chall = params.challenge_bits + EXTRA_TORSION;
    let e_pk = key.curve.clone();
    let bas_ea = &key.basis;
    let kernel = ec_ladder3pt(
        &chall_coeff[..L::ORDER_WORDS],
        &bas_ea.p,
        &bas_ea.q,
        &bas_ea.pmq,
        &e_pk,
    )?;
    let mut e_chall = iso_isogeny_2chain(&EcIsogEven {
        curve: e_pk,
        kernel,
        length: f_chall,
    })?;

    // challenge_and_aux_basis_verify, the canonical bases
    let reduced_order = params.response_bits + HD_EXTRA_TORSION;
    let b_chall_can = ec_curve_to_basis_2f_from_hint(&mut e_chall, reduced_order, hint_chall)?;
    let b_aux_can = ec_curve_to_basis_2f_from_hint(&mut e_aux, reduced_order, hint_aux)?;
    Some(ChallengeState {
        e_chall,
        b_chall_can,
        e_aux,
        b_aux_can,
    })
}

/// The rest of a verification: the basis change on `E_chall`, the
/// `(2^RESPONSE_BITS, 2^RESPONSE_BITS)`-chain to the commitment curve, the
/// challenge recomputed and compared.
pub(crate) fn finish<L: FpBackend + PrimePrecomp>(
    params: &VerifyParams,
    key: &PreparedPublicKey<L>,
    st: ChallengeState<L>,
    m: &[[Scalar; 2]; 2],
    chall_coeff: &Scalar,
    msg: &[u8],
) -> Option<()> {
    let ChallengeState {
        mut e_chall,
        mut b_chall_can,
        e_aux,
        b_aux_can,
    } = st;
    let reduced_order = params.response_bits + HD_EXTRA_TORSION;
    let tmp = b_chall_can.clone();
    let w = L::ORDER_WORDS;
    b_chall_can.p = ec_biscalar_mul_verif(
        &m[0][0][..w],
        &m[1][0][..w],
        reduced_order as usize,
        &tmp,
        &mut e_chall,
    )?;
    b_chall_can.q = ec_biscalar_mul_verif(
        &m[0][1][..w],
        &m[1][1][..w],
        reduced_order as usize,
        &tmp,
        &mut e_chall,
    )?;

    // compute_commitment_curve_verify: E_chall x E_aux -> E_com x E_aux'
    let mut e12 = ThetaCoupleCurve {
        e1: e_chall,
        e2: e_aux,
    };
    let ker = ThetaKernelCouplePoints::from_bases(&b_chall_can, &b_aux_can);
    let out = theta_chain_compute_and_eval(
        params.response_bits as u16,
        &mut e12,
        &ker,
        &[],
        ChainMode::Verify,
    )?;
    let e_com = out.codomain.e1;

    // recompute the challenge and compare
    let chk = hash_to_challenge(params, &key.pk, &e_com, msg);
    (chk == *chall_coeff).then_some(())
}

fn verify_inner<L: FpBackend + PrimePrecomp>(
    params: &VerifyParams,
    key: &PreparedPublicKey<L>,
    sig: &Signature<L>,
    msg: &[u8],
) -> Option<()> {
    if !check_canonical_basis_change_matrix(params, sig) {
        return None;
    }
    let st = challenge_state(
        params,
        key,
        &sig.e_aux_a,
        &sig.chall_coeff,
        sig.hint_chall,
        sig.hint_aux,
    )?;
    finish(
        params,
        key,
        st,
        &sig.mat_bchall_can_to_bchall,
        &sig.chall_coeff,
        msg,
    )
}

/// One item of [`verify_batch`]: a prepared key, a message and an encoded
/// signature.
pub struct BatchItem<'a, L: FpBackend> {
    /// The key, prepared once.
    pub key: &'a PreparedPublicKey<L>,
    /// The message.
    pub msg: &'a [u8],
    /// The encoded signature.
    pub signature: &'a [u8],
}

/// Verify a batch of encoded signatures under prepared keys: `out[i]` is
/// the verdict of `items[i]` (decoding failures are rejections); returns
/// the number accepted. A plain loop with no allocation, the shape a
/// caller spreads across threads.
pub fn verify_batch<L: FpBackend + PrimePrecomp>(
    params: &VerifyParams,
    items: &[BatchItem<'_, L>],
    out: &mut [bool],
) -> usize {
    let mut accepted = 0;
    for (item, verdict) in items.iter().zip(out.iter_mut()) {
        *verdict = match signature_from_bytes::<L>(params, item.signature) {
            Some(sig) => verify_prepared(params, item.key, &sig, item.msg),
            None => false,
        };
        accepted += *verdict as usize;
    }
    accepted
}
