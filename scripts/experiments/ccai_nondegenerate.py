"""
Experiment B: non-degenerate CCAI splitting verification (E_A != E_com).

The theta-chain CCAI experiment (Experiment 1) had E_A = E_com = E0 because its
prime only carried 3-power odd torsion, forcing the component to be an
endomorphism. Here we pick primes whose odd torsion has TWO coprime smooth
factors summing (with the 2-power) appropriately, so BOTH the component and the
auxiliary are genuine smooth-degree isogenies from a common source E_com:

    alpha : E_com -> E_A     (degree q,       FIXED)   => E_A != E_com
    beta  : E_com -> E_aux   (degree 2^a - q, VARIED)  => many E_aux

The gluing on E_A x E_aux has kernel graph{(alpha(R), beta(R))}. By Kani the
source E_com is a codomain factor. We confirm this on the real theta (2,2)-chain
for many auxiliaries, with E_A != E_com throughout, plus a non-diamond control.

Primes:
  a=3: p=479   (p+1 = 2^5*3*5),    q=3,  aux=5
  a=5: p=51839 (p+1 = 2^7*3^3*5),  q=5,  aux=27  (~ up to 36 distinct E_aux)
"""

import sys
sys.path.insert(0, ".")
from sage.all import *
from richelot_isogenies.richelot_isogenies import compute_richelot_chain
from montgomery_isogenies.isogenies_x_only import (
    isogeny_from_scalar_x_only, evaluate_isogeny_x_only,
)
from utilities.supersingular import torsion_basis, fix_torsion_basis_renes
from utilities.strategy import optimised_strategy


def run(p, a, q, aux, n_aux_samples):
    p, a, q, aux = ZZ(p), ZZ(a), ZZ(q), ZZ(aux)
    assert q + aux == 2**a and gcd(q, aux) == 1
    F = GF(p**2, name="i", modulus=[1, 0, 1])
    E_com = EllipticCurve(F, [1, 0])          # commitment curve (j=1728)
    assert E_com.is_supersingular()
    M = 2**(a + 2)

    # 2^(a+2) torsion basis on E_com (used to carry the kernel through the chain)
    T1, T2 = torsion_basis(E_com, M)
    T1, T2 = fix_torsion_basis_renes(T1, T2, a + 2)

    # FIXED component alpha: E_com -> E_A of degree q  (genuine isogeny)
    Pq, Qq = torsion_basis(E_com, q)
    alpha, E_A = isogeny_from_scalar_x_only(E_com, q, ZZ(1), basis=(Pq, Qq))
    aT1, aT2 = evaluate_isogeny_x_only(alpha, T1, T2, M, q)   # on E_A
    jA, jcom = E_A.j_invariant(), E_com.j_invariant()
    assert jA != jcom, "component must be non-degenerate (E_A != E_com)"

    N_constant = F(aux + q) / F(aux - q)      # N1=deg beta side, N2=deg alpha side
    strategy = optimised_strategy(a - 1)

    # VARY auxiliary beta: E_com -> E_aux of degree aux
    Pa, Qa = torsion_basis(E_com, aux)
    seen = {}
    results = []
    scalars = list(range(aux)) + [None]       # <P+sQ> for s, plus <Q>
    for s in scalars:
        if len([r for r in results]) >= n_aux_samples:
            break
        try:
            if s is None:
                beta, E_aux = isogeny_from_scalar_x_only(E_com, aux, ZZ(1), basis=(Qa, Pa))
            else:
                beta, E_aux = isogeny_from_scalar_x_only(E_com, aux, ZZ(s), basis=(Pa, Qa))
        except Exception:
            continue
        jaux = E_aux.j_invariant()
        if jaux in seen:
            continue
        seen[jaux] = True
        bT1, bT2 = evaluate_isogeny_x_only(beta, T1, T2, M, aux)   # on E_aux
        ker = (4 * aT1, 4 * aT2, 4 * bT1, 4 * bT2)
        try:
            Phi, (F1, F2) = compute_richelot_chain(ker, a, N_constant, strategy)
            if F1 is None or F2 is None:
                split, ecom = False, None
            else:
                split = True
                ecom = bool(F1.is_isomorphic(E_com) or F2.is_isomorphic(E_com))
        except Exception as e:
            split, ecom = False, f"err:{type(e).__name__}"
        results.append((str(jaux), split, ecom, str(jA != jaux)))

    # ---- non-diamond control: random anti-isometry-shaped kernel ----
    # Replace the component-side images by a random 2-torsion combination that is
    # NOT alpha's image, so the kernel is not a genuine diamond.
    beta0, E_aux0 = isogeny_from_scalar_x_only(E_com, aux, ZZ(2), basis=(Pa, Qa))
    b0T1, b0T2 = evaluate_isogeny_x_only(beta0, T1, T2, M, aux)
    rT1 = 3 * aT1 + 5 * aT2   # scrambled component side (not a genuine image)
    rT2 = 2 * aT1 + 7 * aT2
    ker_ctrl = (4 * rT1, 4 * rT2, 4 * b0T1, 4 * b0T2)
    try:
        Phi, (F1, F2) = compute_richelot_chain(ker_ctrl, a, N_constant, strategy)
        ctrl_split = not (F1 is None or F2 is None)
    except Exception:
        ctrl_split = False

    n = len(results)
    all_split_ecom = all(r[1] is True and r[2] is True for r in results)
    print(f"\n===== p={p} (p+1={factor(p+1)}), a={a}, gluing (2^{a},2^{a}); "
          f"component deg q={q}, auxiliary deg={aux} =====")
    print(f"  j(E_A)   = {jA}")
    print(f"  j(E_com) = {jcom}   [E_A != E_com: {jA != jcom}]")
    for (jaux, split, ecom, distinct) in results:
        print(f"    E_aux j={jaux[:24]:<24} split={split}  E_com_factor={ecom}  "
              f"E_aux!=E_A:{distinct}")
    print(f"  distinct E_aux tested: {n};  ALL split with E_com factor: {all_split_ecom}")
    print(f"  CONTROL (non-diamond kernel) splits: {ctrl_split}  (expected: False)")
    return dict(p=int(p), a=int(a), q=int(q), aux=int(aux), jA=str(jA), jcom=str(jcom),
                EA_ne_Ecom=bool(jA != jcom), n_aux=n,
                all_split_with_Ecom=bool(all_split_ecom), control_split=bool(ctrl_split))


import json
out = []
out.append(run(479, 3, 3, 5, n_aux_samples=6))
out.append(run(51839, 5, 5, 27, n_aux_samples=20))
print("\n=== JSON ===")
print(json.dumps(out, indent=2))
