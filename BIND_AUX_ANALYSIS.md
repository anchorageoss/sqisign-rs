# Can SQIsign v2.0 bind the auxiliary curve `E_aux` into the Fiat–Shamir hash?

**Feasibility analysis. No code changed.** Sources: the SQIsign v2.0 specification
(*Version 2.0.1*, July 7 2025, `sqisign.org/spec/sqisign-20250707.pdf`) — signing
is **Algorithm 4.2**, §4.4; verification §4.5; security §10.1–10.2 — and the PRISM
v2 reference (`isogeny_swrk/sqisign_prism_v2/prism.py`).

Motivation: SQIsign v2.0 is not SUF-CMA. A third party rerouting the terminal
`l`-isogeny step of the auxiliary isogeny produces a second valid signature on a
new auxiliary curve (see [`LUCA_SINGLE_FACTOR.md`](LUCA_SINGLE_FACTOR.md)). The
proposed mitigation — bind `E_aux` into the challenge, `c = H(pk, j(E_com),
E_aux, msg)` — is examined here against the actual signing dependency structure.

**Notation map** (question → spec): `q = q_rsp` (odd part of the response degree),
`a = e'_rsp` (2-D isogeny exponent), auxiliary degree `2^a − q = 2^{e'_rsp} −
q_rsp`, `A_aux = E_aux`, `sigma = φ_rsp`, `E_com = E_com`.

---

## Q1. Does `E_aux` depend on the challenge `c`? — **Yes, causally and unavoidably.**

Algorithm 4.2 (`SQIsign.Sign`, §4.4, p.40) fixes the order of operations. Tracing
`E_aux` back to the challenge, line by line:

```
 7:  E_com, P_com, Q_com ← IdealToIsogeny(I_com)          // commitment curve
10:  chl ← HASH(pk || j(E_com) || msg)                    // CHALLENGE (hashes E_com, NOT E_aux)
11:  (c1,c2) ← M_sk · (1, chl)                            // depends on chl
12:  I'_chl ← KernelDecomposedToIdeal2f(c1,c2)            // depends on chl
13:  I_chl ← [I_sk]_* I'_chl                              // depends on chl
14:  α_rsp ← RandomEquivalentQuaternion(I_com ∩ I_sk·I_chl)   // depends on I_chl → chl
15:  α_rsp, nbt ← ComputeBacktrackingAndNormalize(α_rsp)
16:  d_rsp ← nrd(α_rsp)/D_mix² · 2^(f−nbt)                // response degree
17:  r_rsp ← DyadicValuation(d_rsp)
18:  q_rsp ← d_rsp / 2^{r_rsp}                            // = "q"
20:  e'_rsp ← e_rsp − r_rsp − nbt                         // = "a"
23:  I_aux ← RandomIdealGivenNorm(2^{e'_rsp} − q_rsp, false)   // AUX ideal, norm depends on q_rsp, e'_rsp
24:  E'_aux, P'_aux, Q'_aux ← IdealToIsogeny(I_com,rsp ∩ I_aux)
27:  E_aux, … ← SplitAuxiliaryIsogeny(E_com, E'_aux, …, q_rsp, e'_rsp, r_rsp)   // AUX CURVE produced here
38:  σ ← (E_aux, nbt, r_rsp, M_chl, chl, hint_aux, hint_chl)
```

The dependency chain is a strict causal cascade:

```
E_com ──(line 10)──▶ chl ──▶ I_chl ──▶ α_rsp ──▶ d_rsp ──▶ (q_rsp, r_rsp, e'_rsp) ──▶ I_aux ──▶ E'_aux ──▶ E_aux
     commitment      challenge  response ideal   response  degree decomposition   aux ideal        aux curve
```

- **`E_aux` depends on the challenge `c`:** yes — four ideal-theoretic steps
  downstream of line 10.
- **`E_aux` depends on `E_com`:** yes — both directly (`SplitAuxiliaryIsogeny`
  input, line 27, and `I_com,rsp = O_0 α_rsp + O_0(q_rsp D_mix)`, line 19) and
  transitively through `chl`.
- **`E_aux` depends on the response `sigma = φ_rsp`:** yes — it is a **byproduct**
  of embedding `φ_rsp` into the 2-D (Kani) isogeny. §4.4.3 (p.41–42): the auxiliary
  ideal `I_aux` is sampled with norm `2^{e'_rsp} − q_rsp` *precisely to complete the
  response degree* to a power of two, then `SplitAuxiliaryIsogeny` (Algorithm 4.5,
  p.43) evaluates the 2-D isogeny `Φ` whose kernel is fixed by `φ_rsp`'s torsion.
  `E_aux` is not chosen freely; its very *degree* `2^{e'_rsp} − q_rsp` is dictated by
  the response (`q_rsp` is the odd part of `nrd(α_rsp)`, and `e'_rsp` its 2-adic
  bookkeeping).

**The hash hashes `E_com`, never `E_aux`.** Confirmed three times: signing line 10
(`HASH(pk || j(E_com) || msg)`); verification (§4.5, p.45, line 29:
`chl' ← HASH(pk || j(F1) || msg)` where `F1 ≅ E_com` is reconstructed by the 2-D
isogeny); and §10.2.5 (p.~68, "inputs of the form `pk || j(E_com) || msg`").

**This is exactly the SUF-CMA gap.** In verification, `E_aux` is an *input* used to
reconstruct `F1 ≅ E_com`, but `E_aux` itself is never hashed. Any alternative
`E_aux''` that reconstructs the same `F1` yields the same `chl` and verifies — which
is what the terminal-`l`-step reroute produces.

---

## Q2. Could `E_aux` be moved to the commitment phase? — **No** (for v2.0 as specified).

To include `E_aux` in the challenge the signer would have to know `E_aux` **at line
10**, before `chl` exists. But `E_aux` is produced at line 27, and its computation
consumes `chl` (line 10) → `I_chl` (13) → `α_rsp` (14) → `q_rsp, e'_rsp` (18, 20) →
`I_aux` (23). This is a genuine **circular dependency**:

```
        need E_aux to compute chl  ◀── (proposed) hash includes E_aux
        need chl to compute E_aux  ◀── (actual) lines 10→27
```

A fixed-point iteration (`chl_0 = H(…no aux…)` → response → `E_aux_0` → `chl_1 =
H(…E_aux_0…)` → …) does not help: each new `chl` induces a fresh lattice
`I_com ∩ I_sk·I_chl`, a fresh short vector `α_rsp`, a fresh degree `d_rsp`, and hence
an unrelated `E_aux`. There is no reason for the map to have a fixed point, and no
efficient way to find one.

Crucially, **even the *degree* of `E_aux` is unknown until after the response is
sampled.** `I_aux` has norm `2^{e'_rsp} − q_rsp` (line 23), where `q_rsp` is the odd
part of `nrd(α_rsp)` and `e'_rsp = e_rsp − r_rsp − nbt` depends on the 2-adic
valuation and backtracking of that same random `α_rsp` (lines 15–20). The signer
cannot pre-sample an ideal — let alone fix a curve — of a norm it does not yet know.

---

## Q3 / Q4. What restructuring would be needed, and the fundamental obstacle.

Two restructuring routes, matching the plan's option (a)/(b)/(c):

**Route A — pre-pick `I_aux` at commitment, then force the response to match.**
Sample `I_aux` (hence `E_aux`) during the commitment phase, hash it, then given
`chl` compute `α_rsp` such that the induced `(q_rsp, e'_rsp)` satisfy
`2^{e'_rsp} − q_rsp = nrd(I_aux)`.

*Obstacle (algebraic, fundamental).* `α_rsp` is a short vector of the
challenge-dependent lattice `I_com ∩ I_sk·I_chl` inside a norm ball
(`RandomEquivalentQuaternion`, Algorithm 4.3, `nrd(α) < D_rsp·D_mix²·2^f`). Its
reduced norm — and therefore `q_rsp` (odd part) and `e'_rsp` (2-adic part) — is
whatever the lattice sampler returns; it is **not a free parameter the signer can
dictate**. Forcing `q_rsp` and `e'_rsp` to values that hit a *pre-committed*
`nrd(I_aux)` is a constrained norm-equation problem: find `α ∈ I_com ∩ I_sk·I_chl`,
`nrd(α) ≤ bound`, with `nrd(α)/D_mix²` having a *prescribed* odd part and 2-adic
valuation. This is strictly harder than the free sampling the scheme relies on and
is not known to be efficiently solvable — it would demand a new `RepresentInteger`/
`IdealToIsogeny` variant that accepts a target auxiliary norm as a hard constraint.
The spec offers no such subroutine, and the difficulty of sampling isogenies of a
*prescribed non-smooth degree* without `End` is exactly the assumption invoked in
§10.2.3 (Signature forgery). **This is the fundamental obstacle.**

**Route B — change *what* is bound (bind the response, not the commitment).** Give
up the commitment-curve challenge and instead make the hash fix the *response
degree*, then verify the response against `pk` directly. This is a different scheme
— and is essentially what PRISM does (Q6). It does not "move `E_aux` to the
commitment"; it removes the free auxiliary curve altogether.

**Summary answers.**
1. Does `E_aux` depend on `c`? **Yes** — chain in Q1.
2. Move `E_aux` to commitment? **No** — circular dependency; its degree is
   response-derived.
3. Changes needed for a "pre-pick" fix: a constrained-norm ideal-to-isogeny that
   targets a fixed auxiliary norm — not available and not known to be efficient.
4. Fundamental obstacle: the auxiliary degree `2^{e'_rsp} − q_rsp` is a function of
   the *reduced norm of the challenge-dependent response quaternion*; it cannot be
   fixed ahead of the challenge without solving a hard prescribed-degree isogeny
   problem.

---

## Q5. What does the spec itself say?

- **Only EUF-CMA is claimed or proved.** p.2 ("proven EUF-CMA secure"); §10.1
  (p.~65–67): Theorem reduces **EUF-CMA** to `q-unif-hint-OneEnd`. There is **no
  SUF-CMA / strong-unforgeability statement anywhere**, and no discussion of
  auxiliary-curve malleability or a second valid signature on the same message.
  The gap is out of scope of the proof, not closed by it.
- **§4.4.3 / Figure 4 (Auxiliary isogeny diagram, p.42)** make the auxiliary an
  explicit *byproduct* of the Kani embedding of the response: `I_aux` exists only to
  pad the response degree `q_rsp` up to `2^{e'_rsp}`, and `SplitAuxiliaryIsogeny`
  (Algorithm 4.5) derives `E_aux` from the response's 2-D isogeny. The auxiliary is
  **constrained, not free** — its degree is `2^{e'_rsp} − q_rsp`.
- **§10.2.5 (Attacking the Fiat–Shamir transform, p.~68)** states the hash input is
  `pk || j(E_com) || msg` and reasons only about collision/second-preimage and a
  "deniability attack." It explicitly assumes *`E_com` is the curve that must be
  fixed* ("the curve `E_com` would a priori need to be fixed before starting the
  collision search"). The designers' FS threat model binds `E_com` and treats the
  response (hence `E_aux`) as canonical — precisely the assumption the reroute
  violates.

---

## Q6. How PRISM handles this differently

PRISM (`prism.py`) uses a structurally different challenge and has **no separately
transmitted auxiliary curve**, which removes the reroute handle at its root.

**Hash to a *prime* response degree, with a salt — no commitment curve.**
`hash_to_prime(msg, E_pk, r)` (lines 24–54) hashes `enc(E_pk) || msg || r` and
reduces mod `2^a`, **grinding the salt `r` until the output `q` is prime**
(`is_pseudoprime(q)`, lines 42/53). The challenge is the *prime degree* `q`, derived
from `(pk, msg, salt)` — there is no commitment curve in the hash at all.

**The signature is `(E_rsp, (P_rsp, Q_rsp), r)` (lines 138–140).** Signing (lines
101–145): `q, r ← hash_to_prime(msg, E_pk)`; response ideal norm
`n_rsp = q·(2^a − q)` (line 122); one isogeny `E_rsp ← IdealToIsogeny(I_cra)` (line
131). There is **one response curve `E_rsp`**, not a `(E_com, E_aux)` pair — the
degree-`(2^a − q)` complement lives only inside the 2-D embedding and is never a
free, transmitted curve.

**Verification is a degree/pairing check against `pk`, not a commitment
re-hash.** `PRISM_verify` (lines 148–192): recompute `q` from `(msg, E_pk, r)`;
build the 2-D isogeny `Φ` from kernel `((P_pk, P_rsp), (Q_pk, Q_rsp))` of degree
`2^a` (line 171); check `e(Φ(P), Φ(Q)) ∈ {pair_pk^q, pair_pk^{−q}}` (lines 182–186)
— i.e. verify the response has degree **exactly the prime `q`** relative to `E_pk`.

**Why this addresses the SUF issue at the level relevant here.**
1. **The dependency is inverted.** In SQIsign v2.0 the constraint (`E_aux`, its
   degree) is *downstream* of the challenge; in PRISM the challenge **is** the degree
   `q`, fixed by the hash *before* the response is computed. PRISM already solves
   the "determine the constraint before hashing" problem — for `q`.
2. **The transmitted response has *prime* degree.** The SQIsign v2.0 reroute needs a
   small factor `l | (2^a − q)` of the transmitted **auxiliary** isogeny. PRISM
   transmits the degree-`q` response with `q` **prime** and pins that degree by
   pairing — a prime-degree isogeny has no small terminal step to reroute, and the
   degree-`(2^a − q)` complement is not a transmitted curve one can swap. The exact
   handle the v2.0 forgery uses is absent by construction.
3. **The salt `r` is in the signature and re-hashed**, giving deterministic,
   verifier-recomputable challenge derivation without a commitment curve.

*Caveat (not over-claiming).* PRISM does not hash `E_rsp` either, and a full
SUF-CMA proof for PRISM is beyond this spec-reading analysis. The point of the
comparison is mechanistic: PRISM's salted **hash-to-prime-degree** design fixes the
response-degree constraint ahead of the response and eliminates the free
smooth-degree auxiliary curve — the two structural features that make "bind `E_aux`
into the hash" both necessary and impossible inside SQIsign v2.0's
commitment→challenge→response→auxiliary cascade.

---

## Bottom line

Binding `E_aux` into the SQIsign v2.0 Fiat–Shamir hash is **not possible without
redesigning the response-generation subroutine**, because `E_aux` — down to its
degree `2^{e'_rsp} − q_rsp` — is a deterministic byproduct of the
challenge-dependent response quaternion `α_rsp` (Algorithm 4.2, lines 10→27). The
signer cannot know `E_aux` at hashing time (line 10), and cannot force it to a
pre-committed value without solving a prescribed-non-smooth-degree isogeny problem
(the §10.2.3 hardness assumption, working *against* the honest signer). The spec
only ever targets EUF-CMA and its FS analysis (§10.2.5) assumes `E_com` is the sole
bound curve. A real fix changes *what* is bound rather than adding `E_aux` to the
existing hash — PRISM's salted hash-to-prime-degree response is the concrete
alternative, fixing the degree before the response and removing the transmitted
smooth auxiliary curve entirely.
