# The compact (dimension-4) format on the round-3 primes: derivation (P29, step 1)

This is the document the plan asks for before any code: the parameters of
SQIsign round 3 with the response embedded in dimension 4, derived at the
round-3 primes, with every choice marked **theirs** (SQIsignHD, ePrint
2023/436, or the round-3 specification) or **ours**. The format is the
round-2 compact format of sqisign-rs 0.4 (SQIsignHD's FastVerify on a
SQIsign prime, 108 bytes at `p = 5·2^248 − 1`, history at `987d1e2`) with
`p` changed and the exponents derived again. The SQIsignHD authors have
published parameters only for their own primes (`13·2^126·3^78 − 1` and the
level III and V analogues); everything at the round-3 primes is our
instantiation, in the posture of ECLIPSE.

**Verdict of step 1: the predicted level I signature is 142 bytes
(with the challenge recomputed from the transmitted commitment curve, as
SQIsignHD and the round-2 format do) or 158 bytes with the
challenge on the wire; both clear the 200-byte criterion, as do
214 / 280 (or 238 / 312) against 306 / 406 at levels III and V.
Steps 2 to 5 may proceed.**

## 0. What the format is

Round-3 SQIsign up to the response (**theirs**, the specification: key
generation, the commitment `E_com`, the challenge `φ: E_pk → E_chall` of
degree `2^λ` with kernel `P + [c] Q` on the canonical basis, the hint
convention, the challenge hash over `(E_pk, E_com, msg)`), then SQIsignHD's
response in place of the specification's two-dimensional one (**theirs**,
ePrint 2023/436 Section 4): the signer samples an ideal `I` equivalent to
the response class with norm `q < 2^e` such that `2^e − q` is a prime
`≡ 1 (mod 4)` ("`2^e`-good"), and transmits the action of the isogeny
`σ: E_com → E_chall` of degree `q` on the `2^r`-torsion as three of the four
scalars of a `2x2` matrix over `Z/2^r`, plus `q`. The verifier writes
`2^e − q = a_1^2 + a_2^2` (Cornacchia), builds the kernel of the
`(2^e, 2^e, 2^e, 2^e)`-isogeny `F` on `E_com^2 × E_chall^2` from `a_1, a_2`
and the transmitted action (Kani, ePrint 2023/436 Lemma 4 and Section 4.4),
and computes `F` as two half-chains of `⌈e/2⌉` steps from both ends that must
meet in the middle (Section 4.3). No auxiliary isogeny and no auxiliary
curve: the commitment curve is transmitted instead.

## 1. The exponent `e` and the torsion `r`

**Theirs.** The response degree is bounded by `q < 2^e` (Section 4.2).
Lemma 12: every left ideal class contains an ideal of norm at most
`2√(2p)/π ≈ 0.90 √p`, so `2^e = Ω(√p)` is necessary, and the smallest
values of `q` are close to `√p`. Sufficiency is heuristic: "assuming that
`q_J(α)` behaves like a random integer, we should expect to find a suitable
`α ∈ J` with probability `O(1/log p)`. Hence, taking `ℓ^e` a few bits over
`√p` might be sufficient" (Section 4.2); a proof would need
`ℓ^e = ω(p^2)` (their Appendix, RigorousSQIsignHD). Their Example 22 at
NIST-I takes `e = 142` for `p ≈ 2^253.3`, that is 15.3 bits over
`log_2 √p`; the RUGDIO of Definition 20 (the oracle the HVZK simulator of
Theorem 21 uses) returns isogenies of `2^e`-good degree, so `e` is where
the simulation argument's "good degree" lives, and nothing in the argument
fixes its margin over `√p` beyond existence.

**Theirs, the SQIsignHD library's rule on SQIsign primes.** The
reference library (SQISignHD-lib, `Signature/scripts/parameters.py` and
`Verification/Verify.py`, master at `0287566`) is instantiated on the
round-2 SQIsign primes `5·2^248 − 1`, `65·2^376 − 1`, `27·2^500 − 1`, and
this is what the round-2 compact format ported. Its exponents are
`e = λ + ⌈log_2(2λ)⌉ = 136 / 201 / 265` and `r = ⌈e/2⌉ + 2 = 70 / 103 / 135`,
with a `2^λ` challenge. Written in `λ`, the rule leans on `p ≈ 2^(2λ)` at
those primes, where it is `log_2 √p + 10.8 / 10.0 / 12.6` bits. It does not
transfer to the round-3 primes as a formula: `3·2^324 − 1` has 70 bits
more than `2^(2λ)` at level I (the two-dimensional verification's torsion),
and `2^136` lies far below the smallest norm any class contains
(`0.90 √p ≈ 2^162.6`, Lemma 12), so no response would exist. What
transfers is its content, `√p` plus about ten bits, which is what
Section 4.2 of the paper says and what we take below.

**Theirs, the round-3 specification's two-dimensional rule.**
`e_rsp = 1 + ⌈(λ + 2 log_2 p)/4⌉ = 196 / 302 / 400`. That response is a
`(2^e_rsp, 2^e_rsp)`-isogeny whose kernel is built from the response *and* an
auxiliary isogeny of degree `2^e_rsp − q`, and the auxiliary curve is sent;
`e_rsp ≈ (log_2 p)/2 + λ/4` because the auxiliary ideal must be found in a
class of norm about `√p` with `λ/4` bits of freedom. In dimension 4 the
auxiliary is `a_1^2 + a_2^2 = 2^e − q` on the *same* curves (Zarhin), so `e`
needs only to exceed `q`, not to leave room for a second ideal: `e` sits a
few bits over `(log_2 p)/2` and no `λ/4` term appears.

**Ours: `e = ⌈(log_2 p)/2⌉ + 11`.** The margin of 11 bits is the one the
library's rule amounts to at the round-2 primes (10.8 / 10.0 / 12.6 bits;
its signer found a good `q` on every signature there) and is the smaller
of the two values in evidence (the paper's own prime has 15.3). Two checks
pin it from below:

- the expected number of ideal-lattice vectors of norm below `2^e` grows as
  `2^(2·margin)` (a rank-4 form of determinant `Θ(p^2)`); with 11 bits that is
  about `2^22` candidates against one good norm in a few hundred (Section 2
  below), so a good `q` exists with overwhelming probability;
- Nakagawa–Onuki (ePrint 2025/1602) recommend a degree of at least `4λ/3`
  bits for a revealed isogeny of large prime-like degree
  (prism-rs `docs/PRISM_4D_FEASIBILITY.md` Section 3): `4λ/3 = 170.7 / 256.0 / 341.3`
  against `e = 174 / 264 / 346`; the margin clears it at every
  level, and 8 bits would not at level I.

Taking SQIsignHD's own 15 bits instead (`e = 178 / 268 / 350`) costs one byte
in `q` and none to one in each scalar; the verdict below does not change.

| level | `p` | `log_2 p` | `log_2 √p` | Lemma 12 minimum (bits) | `e` (ours) | `r = ⌈e/2⌉ + 2` | half-chain steps `⌈e/2⌉` | full chain would need `e + 2` | available `2`-torsion |
|---|---|---|---|---|---|---|---|---|---|
| I | `3·2^324 − 1` | 325.58 | 162.79 | 162.6 | **174** | **89** | 87 | 176 | `2^324` |
| III | `27·2^500 − 1` | 504.75 | 252.38 | 252.2 | **264** | **134** | 132 | 266 | `2^500` |
| V | `17·2^664 − 1` | 668.09 | 334.04 | 333.9 | **346** | **175** | 173 | 348 | `2^664` |

**Theirs: the torsion the verifier needs.** Section 4.3 and Remark 4.2:
`F` of degree `2^e` is computed as `F = F_2 ∘ F_1` from both ends, each half
from points of order `2^(⌈e/2⌉ + 2)` (the two bits above the kernel are the
theta chain's extra torsion; their Section 6.1 writes `f_1 = ⌈e/2⌉ + 2` in
Example 22, `73` for `e = 142`). So the response is transmitted modulo
`2^r`, `r = ⌈e/2⌉ + 2`, and the commitment and challenge bases are rescaled
to order `2^r` (the round-2 format: `r = 70` for `e = 136`). The round-3
primes have `2^324 / 2^500 / 2^664`: a single chain from one end
(`e + 2 = 176 / 266 / 348`) would also fit, unlike at SQIsignHD's
primes, but it would double the scalars' width; **ours** is to keep the
half-chains and `r` as above, for the size.

## 2. Response sampling

**Theirs.** FastRespond (Algorithm 2): sample `I ∼ J` uniformly among
ideals of norm `q < 2^e` (RandomEquivalentIdeal, a uniform sample of a
rank-4 quadratic form below `2^e`), reject unless `q` is `2^e`-good, that
is `2^e − q` prime and `≡ 1 (mod 4)` (Definition 6), and prime to the
challenge degree. The reference signer's `is_good_norm` (as the round-2 port
had it): `q ≡ 3 (mod 4)`, `q < 2^e`, `2^e − q` probably prime (40 rounds).
The verifier repeats the test (FastVerify line 2) and decomposes `2^e − q`
by Cornacchia; a prime `≡ 1 (mod 4)` is a sum of two squares in one way,
so no factoring is needed on either side.

**Ours.** The same test, with the round-3 challenge: `Dφ = 2^λ`, so
"prime to `Dφ`" is `q` odd, which `q ≡ 3 (mod 4)` already gives (the
paper's `q ∧ ℓ' = 1` for its `3`-power challenge has no analogue). We do
not take the two-squares relaxation (accept any `2^e − q = a_1^2 + a_2^2`):
it raises the density by the Landau–Ramanujan factor to
`K/√(ln N) ≈ 0.070 / 0.056 / 0.049`, but the verifier would have to factor `2^e − q`
or receive `(a_1, a_2)`; the P18 note reached the same conclusion for
PRISM-4D.

**Rejection rate.** For `N = 2^e − q` behaving as a random integer of `e`
bits, `P(N ≡ 1 mod 4) · P(N prime | N ≡ 1 mod 4) = (1/4)(2/ln N) = 1/(2 ln N)`:

| level | `e` | `1/(2 e ln 2)` | expected candidates per signature | primality tests (only `q ≡ 3 mod 4` reach one) |
|---|---|---|---|---|
| I | 174 | 1/241 | 241 | 60 |
| III | 264 | 1/366 | 366 | 91 |
| V | 346 | 1/480 | 480 | 120 |

A candidate costs one norm (a few `Ibz` products) and, for one in four, a
Miller–Rabin test on an `e`-bit number whose composites fail at the first
round: 60 / 91 / 120 modular exponentiations of
174 / 264 / 346 bits, under a million cycles at level I against a
round-3 signing of 119 Mcycles (BENCH.md), so the loop is not visible in
signing time. This is where the format is cheap where PRISM-4D was not: the
signer chooses `q` among `2^22` lattice vectors and needs one primality
event, whereas PRISM-4D's hash had to hit a prime `q` with `2^a − q` prime at
once, a second-order event measured at 94.7 / 240.8 / 447.6 Mcycles per
signature (P18, prism-rs `docs/PRISM_4D_FEASIBILITY.md` Section 2).

## 3. Sizes

Fields (**theirs** unless marked): the commitment curve `A_com` (`2·FP`),
the degree `q` (`⌈e/8⌉` bytes), three of the four response scalars modulo
`2^r` (`⌈r/8⌉` bytes each; the fourth from `a·d − b·c ≡ k·q (mod 2^r)`,
Section 6.1 equation (2), with `k = log e(P_com, Q_com) / e(P_chall, Q_chall)`
at `2^r`, the selector implicit in `a`'s parity as in the round-2 format),
and two basis-hint bytes for `E_com` (**theirs**: the library's
basis-from-hint convention, a table index per basis point, kept so that the
unmodified Sage verifier is the oracle; the round-2 format packed them into
spare bits of `A_com`, which the 326-bit prime leaves only two of per
component). `E_chall`'s basis is derived from the challenge, not hinted. The
challenge is not transmitted: the verifier
recomputes it from `(E_pk, E_com, msg)` exactly as the round-3 challenge
hash does, and `E_com` is on the wire because `F` starts from
`E_com^2 × E_chall^2`. The plan's field list carries the challenge as well;
the table gives both.

| level | `A_com` | `q` | 3 scalars | hints | **total** | with the challenge (`λ/8`) | if all four scalars were needed | round-3 signature | criterion |
|---|---|---|---|---|---|---|---|---|---|
| I | 82 | 22 | 3 × 12 = 36 | 2 | **142** | 158 | 154 | 200 | < 200: clears |
| III | 128 | 33 | 3 × 17 = 51 | 2 | **214** | 238 | 231 | 306 | < 306: clears |
| V | 168 | 44 | 3 × 22 = 66 | 2 | **280** | 312 | 302 | 406 | < 406: clears |

The public key carries `A_pk` and its two basis hints in the library's
convention (**theirs**), 84 / 130 / 170 bytes: one byte more than the
round-3 key, whose single hint follows the specification's convention. The
two conventions give different bases, so a compact key is not a round-3 key
(as at round 2); the `CompactPublicKey` of the plan.

**The fourth scalar.** The relation `e_{2^r}(σ(P), σ(Q)) = e_{2^r}(P, Q)^q`
holds exactly for points of order `2^r` and an isogeny of degree `q`, so
`a·d − b·c ≡ k·q (mod 2^r)` determines the dropped scalar completely:
zero residual bits, as in ECLIPSE and as the round-2 format used (three
scalars of exactly `r` bits, no hint). The two-dimensional round-3 format's
four residual bits came from a relation that holds only to the fourth power
(`COMPRESSION.md`); nothing of that kind is present here. Step 4 measures it
on 1000 signatures rather than assumes it; the "all four" column is the
fallback and still clears.

## 4. Security statement

In the register of SECURITY.md.

**Theirs, proven in ePrint 2023/436.** Special soundness of the
identification protocol (Proposition 17) given `q ∧ Dφ = 1`, hence a proof
of knowledge of a non-scalar endomorphism of `E_pk` with knowledge error
`1/#C`; computational HVZK in the RUGDIO model assuming the commitment is
computationally indistinguishable from a uniform supersingular curve
(Theorem 21); EUF-CMA of the Fiat–Shamir signature in the random-oracle and
RUGDIO model under the endomorphism-ring problem (Section 6, opening).
Completeness of the dimension-4 representation (Lemma 4, `q ∧ 2 = 1`).

**Ours, assumed or changed.**

1. `p` is the round-3 prime. The arguments are stated for
   `p = c·2^f·3^f' − 1` with a `3`-power challenge; they use of the prime
   only the accessible `2`-torsion and `p = Θ(2^2λ)`. Both hold here
   (`2^324 / 2^500 / 2^664` against `r = 89 / 134 / 175`).
2. The challenge is round 3's `2^λ`-isogeny, so `#C = 2^λ` (soundness
   `λ` bits as in the specification) and `q ∧ Dφ = 1` is `q` odd, given by
   the goodness rule. The paper's requirement that `Dφ` be prime to `2` is
   for its signer's torsion evaluation along an alternate path; our signer
   evaluates the response ideal with the round-3 ideal-to-isogeny
   machinery instead, as the round-2 format did.
3. The commitment is round 3's, whose own zero-knowledge argument is what
   makes `E_com` close to uniform; Theorem 21's assumption is inherited
   from the specification rather than from SQIsignHD's commitment.
4. `e` and the margin are ours (Section 1). The existence of a good `q` in
   every class is heuristic in the paper and remains so here; we add the
   Nakagawa–Onuki check that `e ≥ 4λ/3` at every level.
5. The signer's distribution: FastRespond outputs uniform good ideals given
   termination; our sampler will be the round-2 port's (random lattice
   vectors of the response ideal, rejection on goodness), whose uniformity
   over ideals of norm below `2^e` is what the RUGDIO simulator models. If
   step 3 changes the sampler, this paragraph changes.
6. The verifier's acceptance is the half-chain middle-codomain match, the
   Cornacchia decomposition and the recomputed challenge, plus the
   goodness test on `q`. The round-2 port's own note records that the
   paper's final image condition (its "stage 6", the evaluation of `F` on
   a point to confirm the embedded isogeny is `σ` and not another isogeny
   with the same kernel data) was out of scope there; step 3 states what
   the round-3 verifier checks and closes that item or records why the
   middle match suffices.
7. The format is not strongly unforgeable and does not try to be
   (SECURITY.md, "Formats"): the scalars carry the
   sign freedom of Kummer points and the hints are not verified as
   canonical, as in the round-2 format.

## 5. Memory

A dimension-4 theta point is 16 `F_p^2` coordinates. The half-chain
verifier of the round-2 port stores, per step, the kernel basis (4 points),
the codomain null point (1) and the image precomputation (16 `F_p^2`, one
point's worth), for both halves, and walks an optimal strategy whose stack
of kernel bases is bounded by the chain length.

| level | `F_p^2` (limbs × 2 × 8 B) | theta point | per step (6 points) | half chain of `⌈e/2⌉` steps | strategy stack, worst case (`⌈e/2⌉` × 4 points) |
|---|---|---|---|---|---|
| I | 96 B | 2 KiB | 9 KiB | 783 KiB | 522 KiB |
| III | 128 B | 2 KiB | 12 KiB | 1.55 MiB | 1.03 MiB |
| V | 176 B | 3 KiB | 16 KiB | 2.79 MiB | 1.86 MiB |

Peak with both half-chains kept as the round-2 code keeps them: about
1.53 MiB / 3.09 MiB / 5.58 MiB; the second half needs only its last
codomain for the match, so a verifier that discards it as it goes peaks at
one half chain plus the stack. None of this is a constraint on a host; the
dimension-4 layer allocates (it is the one `alloc` user of the verify
crate), so `thumbv7em` stays a build check, not a memory claim.

## 6. Predicted verification cost (to be measured in step 4)

The round-2 dimension-4 verifier took 97.1 Mcycles at `p = 5·2^248 − 1`
(`e = 136`, two half-chains of 68 steps, 4-limb field) against 12.6 for the
dimension-2 verifier on the same box, 7.7x (P18). At round-3 level I the
half-chains have 87 steps (×1.28) on a 6-limb field
(×2.25 by limb count, less with the P26 assembly kernels the old layer did
not have), so a first estimate is 150 to 280 Mcycles, 12 to 23x the
round-3 dimension-2 verification of 12.3 Mcycles (BENCH.md). The number to
report is the measured one.

## 7. What step 2 can be

The library's Sage is a verifier only (`Verification/Verify.py` with the
`Theta_dim4` package; the signer is C, on its own field backends for the
round-2 primes). Re-parameterising the Sage verifier to `3·2^324 − 1` is a
level entry (`p`, `c`, the `2`-adic exponent), `e = 174`, `r = 89`, and the
two non-residue tables (`Data/NQR_TABLE`, `Data/Z_NQR_TABLE`, twenty
elements each) regenerated for the new prime. It cannot sign, so the plan's
"keygen, sign, verify in Sage" becomes what the round-2 format did:
signatures come from our signer (the round-3 flow with the response
sampler of Section 2 and the `2^r` action matrix), the Sage verifier is the
oracle that accepts them and rejects the negatives, and the toy-prime
exhaustive checks of the plan (embedding, sampling, splitting) run in Sage
where the arithmetic is small.

### Step 2, status (2026-09-25)

Done: the library cloned with its `Theta_dim4` submodule and self-tested;
the two non-residue tables regenerated for `3·2^324 − 1` with the library's
rule; the oracle driver `sqisignhd-harness/round3/verify_r3.py` (the
unmodified verifier under our parameters, reading the library's vector
format); and the real-prime check of the dimension-4 half-chain at `e =
174`, `r = 89`: a 173-bit response degree (`3^109`) on a generic domain,
half-chains of 87 + 87 steps in 1.5 s of Sage, middle-codomain match and
image check both true on two seeds (`sqisignhd-harness/round3/README.md`).
Two facts recorded there: the goodness rule forces `q ≡ 3 (mod 4)`, so test
isogenies must have such a degree (odd powers of 3 do; `E_0`'s cheap
endomorphisms do not), and chains through `E_0` hit products while walks
from a generic curve do not. Open: signatures for the oracle (they need the
signer of step 3), the negatives, and the toy-prime exhaustive runs.

### Step 3 to 5, status (2026-09-25): level I built

Built behind the `compact` feature, level I only; levels III and V stay
derived. What the build settled beyond the derivation:

- **The split element.** With the library's hinted basis the challenge
  kernel `P + [c] Q` lies above `(0, 0)` on `E_0` for even `c`, a class the
  round-3 basis never produces; there the challenge ideal is `I_c · s`
  with `s = 1 − i`, its isogeny `φ_{I_c} ∘ s`, and the response's connecting
  element is `δ = conj(s) β̄ γ' / n(K')` (the signer's docs). Without `conj(s)`
  every even challenge failed the challenge division.
- **Square roots.** The library's basis lift, difference point and `α` use
  the even-normalised root (`sqrt_Fp2_det`, the crate's
  `sqrt_canonical_even`); its challenge images are lifted with the FESTA
  root (`sqrt_Fp2`). With the round-3 crate's own root the oracle accepted
  two of five vectors; with the library's, five of five.
- **The image check.** The 0.4.x verifier already ran the paper's final
  image condition (`hd_image_check`); its module header saying otherwise
  was stale. The check stays.
- **Sizes.** 142-byte signatures, 84-byte keys, as predicted.
- **Cost** (own run, level I): keygen 42.8, sign 538 (mean; a rejection
  loop of about 240 response samples), verify 122.2 Mcycles, 9.9x the
  two-dimensional verifier (the prediction of Section 6 was 150 to 280).
- **Oracle.** Five self-generated vectors in the library's text format
  are accepted by its unmodified Sage verifier at `e = 174`, `r = 89`,
  five tampered ones rejected (`sqisignhd-harness/round3`, workflow
  `compact-oracle`).

## 8. Verdict

Predicted 142 bytes at level I (158 with the challenge on the wire)
against the 200-byte criterion, 214 / 280 against 306 / 406: the
format clears at every level with the parameters above. Level I is built
and measured at 142 bytes; III and V are derived only.
