# Experiment A — verification of equation (5) at toy parameters

**Claim (ePrint 2026/1305, eq. 5).** For the Kani diamond with component
`a: E_A -> E_com` (degree `q`) and auxiliary side `c': E_aux -> E_com`
(degree `2^a - q`), the anti-isometry `psi': E_A[2^a] -> E_aux[2^a]` whose graph
is `ker(Phi)` (with `Phi: E_A x E_aux -> E_com x E'` the `(2^a,2^a)`-gluing) is

```
    psi'(P) = - hat(c') * [(2^a - q)]^{-1} * a(P) ,      P in E_A[2^a].     (5)
```

This is the load-bearing algebraic identity of the SUF-CMA argument. We verify
it at toy parameters two independent ways.

## Verdict

**Equation (5) is VALIDATED.** It matched the geometric anti-isometry in
**every** tested instance — 90/90 in the pointwise sweep (540 point checks) and
3/3 in the theta-chain corroboration — with **no disagreements**.

One precision detail surfaced and must be respected (see *Precision note*): the
inverse `[(2^a - q)]^{-1}` is taken **mod `2^a`**, and (5) is an identity on the
`2^a`-torsion only — not on the `2^(a+2)`-torsion carried through the gluing.

## Setup

Tooling: SageMath 10.9 + the SQIsign2D-West theta `(2,2)`-isogeny library
(`ThetaIsogenies/two-isogenies`, Dartois–Maino–Pope–Robert), which computes the
full `(2^a,2^a)`-chain between elliptic products and returns the split codomain
factors.

Ground truth ("actual kernel"): a genuine diamond from a **common source**
`E_com`, given by two isogenies `alpha: E_com -> E_A` (degree `q`) and
`beta: E_com -> E_aux` (degree `2^a - q`). Then `a = alpha.dual()` and
`hat(c') = beta`, and the kernel of the gluing on `E_A x E_aux` is
`graph{(alpha(R), beta(R)) : R in E_com[2^a]}`, so the geometric anti-isometry is

```
    psi'_true(alpha(R)) = beta(R)          ("actual kernel").
```

## Check 1 — pointwise identity (model-independent), E_A != E_com

For a spread of torsion points we compare `psi'_formula(P)` from (5) against
`psi'_true(P)`, on a basis of `E_A[2^a]` (six points per instance). This is pure
isogeny/point arithmetic and does not touch the theta model.

Swept: primes `p in {191, 1279}` (with `2^a | p+1`), powers `a in {4,5,6,8}`,
degrees `q in {3,5,7,9,15,27}` (covering both `q < 2^a/2` and `q > 2^a/2`),
three distinct commitment curves `E_com`, two auxiliary isogeny choices — all
with `E_A != E_com` (genuine isogeny component, not an endomorphism).

```
SUMMARY: eq.(5) pointwise identity holds on 90/90 instances (6 torsion points each)
VERDICT: VALIDATED
```

## Check 2 — geometric corroboration via the theta chain, E_A = E_com

In the Montgomery-model diamond the library accepts (`E0`, `j = 1728`; component
`a = X + Y*iota` an endomorphism of degree `C = X^2 + Y^2 = q`; auxiliary
`phi_B` of degree `B = 2^a - q`), we build the kernel **from eq. (5)** and run
the real `(2^a,2^a)`-chain.

| param | prime | gluing | (i) `a.dual()∘a = [q]` | (ii) (5) reproduces honest images on `E_A[2^a]` | (iii) (5)-kernel splits, `E_com` factor |
|------:|------:|--------|:---:|:---:|:---:|
| 0 | 19-bit  | `(2^9,2^9)`     | ✓ | ✓ | ✓ |
| 1 | 139-bit | `(2^72,2^72)`   | ✓ | ✓ | ✓ |
| 2 | 254-bit | `(2^126,2^126)` | ✓ | ✓ | ✓ |

So the kernel produced by eq. (5) is exactly the honest anti-isometry on the
`2^a`-torsion, and its theta-chain codomain splits with `E_com` as a factor.

## Precision note (matters for the write-up)

`psi'` is a map on `E_A[2^a]`, so `[(2^a - q)]^{-1}` in (5) is the inverse of
`2^a - q` **modulo `2^a`**. Concretely `-[(2^a-q)]^{-1}·q ≡ 1 (mod 2^a)`, which
is why `psi'_formula(alpha(R)) = beta(R)`.

The gluing is computed with points of order `2^(a+2)` (the theta model needs
"torsion above the kernel"). On those extended points the same expression does
**not** reduce to the identity — it differs in the top 2 bits, because
`-[(2^a-q)]^{-1}·q ≡ 1 (mod 2^a)` but `≢ 1 (mod 2^(a+2))`. This is harmless (the
kernel proper is the `2^a`-part, recovered after the `×4` scaling the chain
applies), but the paper should state that the inverse and the identity are taken
**mod `2^a`**. (This is the same "2 extra torsion bits" phenomenon documented for
the compression decompressor.)

## Reproduction

```
# SageMath >= 10 (used 10.9 via conda-forge)
git clone https://github.com/ThetaIsogenies/two-isogenies.git
cd two-isogenies/Theta-SageMath
# Check 1 (pointwise, broad): copy scripts/experiments/verify_eq5_broad.py here
sage verify_eq5_broad.py
# Check 2 (theta chain): copy scripts/experiments/verify_eq5_theta.py here
sage verify_eq5_theta.py
```
