"""
Experiment A (Check 2): geometric corroboration of eq. (5) via the theta chain.

Build the (2^a,2^a) kernel *from eq. (5)* and confirm the real theta
(2,2)-chain (compute_richelot_chain) splits with E_com as a codomain factor.

We work in the Montgomery-model diamond the theta library accepts:
E0 (=E_com=E_A, j=1728), component a = aux_endo = X + Y*iota (deg C = q, an
endomorphism), auxiliary phi_B: E0 -> EB (deg B = 2^a - q).

- hat(c') = phi_B ,  a = aux_endo ,  a.dual() = X - Y*iota (conjugate, deg C).
- eq. (5):  psi'_formula(P) = - phi_B( [B^{-1}] * a.dual()(P) ),  P in E_A[2^a].
We evaluate phi_B at arbitrary points via BiDLP in the (P2,Q2) basis.

Checks: (i) a.dual() o a == [C];  (ii) psi'_formula(aux_endo(P2)) equals the
honest image phi_B(P2) as an ACTUAL point (so the eq.(5) kernel = honest
kernel); (iii) the theta chain on the eq.(5) kernel splits with E_com = E0.
"""

import sys
sys.path.insert(0, ".")
from sage.all import *
from isogeny_diamond import DIAMONDS
from richelot_isogenies.richelot_isogenies import compute_richelot_chain
from montgomery_isogenies.isogenies_x_only import (
    isogeny_from_scalar_x_only, evaluate_isogeny_x_only,
)
from utilities.supersingular import torsion_basis, fix_torsion_basis_renes
from utilities.discrete_log import BiDLP
from utilities.strategy import optimised_strategy


def run(param_index):
    f, ea, eb, X, Y = DIAMONDS[param_index]
    A = ZZ(2**ea); B = ZZ(3**eb); C = A - B
    assert C == X**2 + Y**2
    p = f * 4 * A * B - 1
    F = GF(p**2, name="i", modulus=[1, 0, 1])
    E0 = EllipticCurve(F, [1, 0])
    iota = E0.automorphisms()[2]
    P2, Q2 = torsion_basis(E0, 4 * A)
    P2, Q2 = fix_torsion_basis_renes(P2, Q2, ea + 2)
    P3, Q3 = torsion_basis(E0, B)
    M = 4 * A                       # point order 2^(ea+2)
    assert iota(iota(P2)) == -P2    # iota^2 = -1

    def aux_endo(P):        # component a, degree C
        return X * P + Y * iota(P)

    def aux_endo_dual(P):   # a.dual() = conjugate, degree C
        return X * P - Y * iota(P)

    # honest auxiliary phi_B and its action on the (P2,Q2) basis
    phiB, EB = isogeny_from_scalar_x_only(E0, B, 7, basis=(P3, Q3))
    phi_P0, phi_Q0 = evaluate_isogeny_x_only(phiB, P2, Q2, M, B)

    def phiB_eval(W):       # phi_B at an arbitrary M-torsion point via BiDLP
        u, v = BiDLP(W, P2, Q2, M)
        return u * phi_P0 + v * phi_Q0

    # (i) a.dual() o a == [C]
    chk_dual = all(aux_endo_dual(aux_endo(T)) == C * T for T in (P2, Q2))

    Binv = int(inverse_mod(B % M, M))

    def psi_formula(P):     # eq. (5)
        return -phiB_eval((Binv * aux_endo_dual(P)) )

    # (ii) eq.(5) is the anti-isometry on E_A[2^a]; test on the 2^a-torsion.
    N = 2**ea
    BinvN = int(inverse_mod(B % N, N))
    P2c, Q2c = aux_endo(P2), aux_endo(Q2)   # basis of E_A[M] (E_A = E0)
    def psi_formula_N(P):                   # eq.(5) restricted to 2^a-torsion
        return -phiB_eval(BinvN * aux_endo_dual(P))
    match_P = (psi_formula_N(4 * P2c) == 4 * phi_P0)
    match_Q = (psi_formula_N(4 * Q2c) == 4 * phi_Q0)

    # (iii) theta chain on the eq.(5)-built kernel splits with E_com = E0
    KP, KQ = psi_formula(P2c), psi_formula(Q2c)
    ker = (4 * P2c, 4 * Q2c, 4 * KP, 4 * KQ)   # E1=E_A parts, then E2=EB parts
    N_constant = F(B + C) / F(B - C)
    strategy = optimised_strategy(ea - 1)
    try:
        Phi, (F1, F2) = compute_richelot_chain(ker, ea, N_constant, strategy)
        split = F1 is not None and F2 is not None
        ecom = bool(F1.is_isomorphic(E0) or F2.is_isomorphic(E0)) if split else None
    except Exception as e:
        split, ecom = False, f"err:{type(e).__name__}"

    print(f"param_index={param_index} (p={p.nbits()} bits, q=C={C}, 2^a-q=B={B}, gluing 2^{ea}):")
    print(f"  (i)   a.dual() o a == [C]:                 {chk_dual}")
    print(f"  (ii)  eq.(5) psi' == honest phi_B images:  P:{match_P}  Q:{match_Q}")
    print(f"  (iii) theta chain on eq.(5) kernel splits: {split}  E_com factor: {ecom}")
    return chk_dual and match_P and match_Q and split and ecom is True


ok = True
for pi in [0, 1, 2]:
    ok = run(pi) and ok
print("\nCheck 2 VERDICT:", "VALIDATED" if ok else "PROBLEM")
