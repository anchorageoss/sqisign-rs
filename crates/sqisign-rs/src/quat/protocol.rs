//! QuaternionResponseComputation: the response ideal of SQIsign, kept for
//! the substrate check (P07).

use super::{QuatAlg, QuatAlgElem, QuatIdeal};
use crate::mp::{DefaultDomain, Ibz, Rng};

/// The response of SQIsign: `(sk_chall_ideal, resp_quat, norm,
/// sk_chall_quat)`. `sk_chall_ideal` is a small prime-norm ideal equivalent
/// to `sk_ideal ∩ ideal_chall_two`; `resp_quat` lies in
/// `conj(sk_chall_ideal) ideal_commit` with norm `norm * n(commit) *
/// n(sk_chall)` and `norm < 2^e`; `sk_chall_quat` is the (conjugated)
/// connecting element. `None` on randomness failure.
#[allow(clippy::too_many_arguments)]
pub fn response_element<const N: usize>(
    sk_ideal: &QuatIdeal<N>,
    ideal_chall_two: &QuatIdeal<N>,
    chall_split: Option<&QuatAlgElem<N>>,
    ideal_commit: &QuatIdeal<N>,
    e: u32,
    alg: &QuatAlg<N>,
    rng: &mut impl Rng,
) -> Option<(QuatIdeal<N>, QuatAlgElem<N>, Ibz<N>, QuatAlgElem<N>)> {
    let inter = QuatIdeal::intersect_o0(ideal_chall_two, sk_ideal, chall_split);
    let (mut sk_chall_quat, sk_chall_ideal) =
        inter.small_equivalent_coprime(Some(&Ibz::zero()), alg, &mut DefaultDomain(&mut *rng))?;
    sk_chall_quat = sk_chall_quat.conj();
    let prod = QuatIdeal::mul_o0(&sk_chall_ideal, ideal_commit);
    let mut bound = Ibz::<N>::one().mul_2exp(e).sub(&Ibz::one());
    bound.set_bound(e as i32 + 1);
    let (resp_quat, mut norm) = prod.sample_from_ball(&bound, alg, rng)?;
    norm.set_bound(e as i32 + 1);
    Some((sk_chall_ideal, resp_quat, norm, sk_chall_quat))
}
