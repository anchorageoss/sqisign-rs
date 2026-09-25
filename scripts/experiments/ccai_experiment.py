"""
CCAI exhaustive-ish search at toy parameters, using the SQIsign2D-West theta
(2,2)-isogeny machinery (ThetaIsogenies/two-isogenies, Dartois-Maino-Pope-Robert).

Setup (ePrint 2026/1305 Assumption 5.7, mapped onto a computable Kani diamond):
The (2^a,2^a) gluing Phi acts on codomain(component) x codomain(auxiliary),
where both the component a and the auxiliary tau emanate from a common source S.
By Kani (1997, Thm 2.3) the codomain of such a gluing splits as  S x E3 , i.e.
the SOURCE reappears as a codomain factor. In the SQIsign2D construction the
source is the commitment curve E_com.

CCAI asks: fixing the component a (hence E_com), can a DIFFERENT auxiliary
tau' (-> E_aux' != E_aux) still yield a split codomain with E_com as a factor?

Experiment: fix the component (aux_endomorphism X + Y*iota, degree C, an
endomorphism of E0 so here E_A = E_com = E0), and vary the auxiliary phi_B over
many secrets -> many distinct E_aux. For each, build the Kani kernel exactly as
the reference does, run the real (2,2)-chain, and check whether E_com (= E0) is
a codomain factor. Also run a CONTROL with a non-diamond (random) anti-isometry
kernel to confirm that splitting is NOT automatic for arbitrary kernels.
"""

import sys
sys.path.insert(0, ".")
from sage.all import *
from isogeny_diamond import DIAMONDS
from richelot_isogenies.richelot_isogenies import compute_richelot_chain
from montgomery_isogenies.isogenies_x_only import (
    isogeny_from_scalar_x_only,
    evaluate_isogeny_x_only,
)
from utilities.supersingular import torsion_basis, fix_torsion_basis_renes
from utilities.strategy import optimised_strategy


def build_diamond(param_index):
    f, ea, eb, X, Y = DIAMONDS[param_index]
    A = ZZ(2**ea); B = ZZ(3**eb); C = A - B
    assert C == X**2 + Y**2
    p = f * 4 * A * B - 1
    F = GF(p**2, name="i", modulus=[1, 0, 1])
    E0 = EllipticCurve(F, [1, 0])
    assert E0.is_supersingular()
    iota = E0.automorphisms()[2]
    P2, Q2 = torsion_basis(E0, 4 * A)
    P2, Q2 = fix_torsion_basis_renes(P2, Q2, ea + 2)
    P3, Q3 = torsion_basis(E0, B)
    N_constant = F(B + C) / F(B - C)
    strategy = optimised_strategy(ea - 1)
    return dict(f=f, ea=ea, eb=eb, X=X, Y=Y, A=A, B=B, C=C, p=p, F=F, E0=E0,
                iota=iota, P2=P2, Q2=Q2, P3=P3, Q3=Q3,
                N_constant=N_constant, strategy=strategy)


def aux_endo(D, P):
    return D["X"] * P + D["Y"] * D["iota"](P)


def run_auxiliary(D, bob_secret):
    """Honest-shape Kani kernel with component = aux_endo (fixed, deg C) and
    auxiliary = phi_B for the given secret (deg B, -> E_aux = EB). Returns
    (E_aux_j, split?, factor_js) from the real (2,2)-chain."""
    E0, A, B = D["E0"], D["A"], D["B"]
    phiB, EB = isogeny_from_scalar_x_only(E0, B, bob_secret, basis=(D["P3"], D["Q3"]))
    phi_P0, phi_Q0 = evaluate_isogeny_x_only(phiB, D["P2"], D["Q2"], 4 * A, B)
    P1 = aux_endo(D, D["P2"]); Q1 = aux_endo(D, D["Q2"])
    ker = (4 * P1, 4 * Q1, 4 * phi_P0, 4 * phi_Q0)
    try:
        Phi, (F1, F2) = compute_richelot_chain(ker, D["ea"], D["N_constant"], D["strategy"])
    except Exception as e:
        return (EB.j_invariant(), False, None, f"chain-error:{type(e).__name__}")
    if F1 is None or F2 is None:
        return (EB.j_invariant(), False, None, "no-split")
    js = (F1.j_invariant(), F2.j_invariant())
    ecom_is_factor = bool(F1.is_isomorphic(D["E0"]) or F2.is_isomorphic(D["E0"]))
    return (EB.j_invariant(), True, js, ecom_is_factor)


def run_control_random(D, seed_pt_scalar):
    """Control: a NON-diamond anti-isometry kernel. We keep the auxiliary side
    (phi_B images) but replace the component side by a random 2^(ea+2)-torsion
    image (scaled multiple of the basis) that is not aux_endo. Expect: does not
    split, or does not recover E_com."""
    E0, A, B = D["E0"], D["A"], D["B"]
    phiB, EB = isogeny_from_scalar_x_only(E0, B, 12345, basis=(D["P3"], D["Q3"]))
    phi_P0, phi_Q0 = evaluate_isogeny_x_only(phiB, D["P2"], D["Q2"], 4 * A, B)
    # Random component-side images: a "random" endomorphism-like combination that
    # does NOT satisfy the norm equation, so the kernel is not a genuine diamond.
    s = seed_pt_scalar
    P1 = s * D["P2"] + (s + 1) * D["iota"](D["Q2"])
    Q1 = (s + 2) * D["Q2"] + s * D["iota"](D["P2"])
    ker = (4 * P1, 4 * Q1, 4 * phi_P0, 4 * phi_Q0)
    try:
        Phi, (F1, F2) = compute_richelot_chain(ker, D["ea"], D["N_constant"], D["strategy"])
    except Exception as e:
        return (False, None, f"chain-error:{type(e).__name__}")
    if F1 is None or F2 is None:
        return (False, None, "no-split")
    js = (F1.j_invariant(), F2.j_invariant())
    ecom_is_factor = bool(F1.is_isomorphic(D["E0"]) or F2.is_isomorphic(D["E0"]))
    return (True, js, ecom_is_factor)


def main():
    import json
    n_samples = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    indices = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0, 1, 2]
    report = []
    for pi in indices:
        D = build_diamond(pi)
        print(f"\n===== param_index={pi}: p = {D['f']}*4*2^{D['ea']}*3^{D['eb']}-1 "
              f"({D['p'].nbits()} bits), component deg C={D['C']}, aux deg B={D['B']}, "
              f"gluing (2^{D['ea']},2^{D['ea']}) =====")
        print(f"E_A = E_com = E0 (j={D['E0'].j_invariant()}); "
              f"component a = X+Y*iota is an endomorphism, deg {D['C']}")
        distinct_aux = {}
        compat = 0
        for k in range(n_samples):
            secret = 1 + k * 7 + 3  # spread over secrets, deterministic
            secret = secret % max(2, D["B"])
            eaux_j, split, js, ecom = run_auxiliary(D, secret)
            key = str(eaux_j)
            if key in distinct_aux:
                continue
            distinct_aux[key] = (split, ecom)
            if split and ecom is True:
                compat += 1
            print(f"  aux secret={secret:>6}: E_aux j={str(eaux_j)[:22]:<22} "
                  f"split={split} E_com_is_factor={ecom}")
        # control
        ctrl_split, ctrl_js, ctrl_ecom = run_control_random(D, 3)
        print(f"  CONTROL (non-diamond kernel): split={ctrl_split} "
              f"E_com_is_factor={ctrl_ecom}")
        n_distinct = len(distinct_aux)
        all_compat = all(v[0] and v[1] is True for v in distinct_aux.values())
        print(f"  => distinct E_aux tested: {n_distinct}; "
              f"all split with E_com as factor: {all_compat}; "
              f"compatible tau count = {compat}")
        report.append(dict(param_index=pi, p_bits=int(D['p'].nbits()),
                           component_deg=int(D['C']), aux_deg=int(D['B']),
                           gluing_2power=int(D['ea']),
                           distinct_Eaux=n_distinct, compatible=compat,
                           all_compatible=bool(all_compat),
                           control_split=bool(ctrl_split),
                           control_ecom_factor=(None if ctrl_ecom is None else bool(ctrl_ecom))))
    print("\n=== JSON ===")
    print(json.dumps(report, indent=2))


main()
