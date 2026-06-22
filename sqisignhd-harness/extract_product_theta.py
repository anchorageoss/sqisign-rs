#!/usr/bin/env sage
r"""
extract_product_theta.py  --  additive oracle for Phase 5b.5 (dim-1 and product
theta structures forming the dim-4 gluing's domain).

Each half-chain's gluing (`KaniGluingIsogenyChainDim4Half`) builds its domain as

    Theta1 = ThetaStructureDim1(E1, T1, T2)          # E1 = E_com
    Theta2 = ThetaStructureDim1(E2, U1, U2)          # E2 = E_chal
    Theta12 = ProductThetaStructureDim2(Theta1, Theta2)        # dim-1 x dim-1 -> dim-2
    ... m dim-2 (2,2)-isogeny steps -> codomain Phi_codomain (dim-2) ...
    domain_product = ProductThetaStructureDim2To4(Phi_codomain, Phi_codomain)  # dim-2 x dim-2 -> dim-4
    domain_base_change = domain_product.base_change_struct(N_dim4)             # gluing domain

This dumps, per signature and per half-chain (F1, F2_dual), the theta-null
points of every structure plus the inputs needed to rebuild them with the
product wrappers, and a generic-point conversion to exercise
`montgomery_point_to_theta_point`:

  * T1_xz / U1_xz: (X:Z) of the canonical 4-torsion P used by each dim-1 null;
  * dim1_null_1 / dim1_null_2: (X+Z, X-Z) (`torsion_to_theta_null_point`);
  * theta12_null: the dim-2 product null (4 coords);
  * dim2_codomain_null: Phi_codomain null (4 coords) -- input to the dim2->dim4 product;
  * domain_product_null: the dim-4 product null (16 coords) = ProductThetaStructureDim2To4 null;
  * conv_pt_1/2: a generic E1/E2 point (X:Z) and its dim-1 theta image, and
    conv_prod: the dim-2 product theta image of the pair.

Captured by wrapping `KaniGluingIsogenyChainDim4Half.__init__` at runtime; no
reference file and no Phase 0 deliverable is modified.

Usage:  sage extract_product_theta.py [--out product_theta_l1.json] [--n N] [--lib PATH]
Env  :  SQISIGNHD_LIB
"""
import sys
import os
import json
import argparse


def find_lib(cli_path):
    cands = []
    if cli_path:
        cands.append(cli_path)
    if os.environ.get("SQISIGNHD_LIB"):
        cands.append(os.environ["SQISIGNHD_LIB"])
    cands += [os.path.join(os.getcwd(), "SQISignHD-lib"),
              os.path.expanduser("~/SQISignHD-lib"), "/home/user/SQISignHD-lib"]
    for c in cands:
        if c and os.path.isdir(os.path.join(c, "Verification", "Theta_dim4",
                                            "Theta_dim4_sage")):
            return c
    raise SystemExit("Could not locate SQISignHD-lib; set --lib or $SQISIGNHD_LIB")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--lib", default=None)
    args = ap.parse_args()
    out_path = os.path.abspath(args.out or "product_theta_l1.json")

    LIB = find_lib(args.lib)
    VERIF = os.path.join(LIB, "Verification")
    os.chdir(VERIF)
    sys.path.insert(0, VERIF)

    import Verify
    from Verify import SQIsignHD, SQIsignHD_verif
    import Theta_dim4.Theta_dim4_sage.pkg.isogenies.Kani_gluing_isogeny_chain_dim4 as kgmod

    pp = SQIsignHD(1)
    p = int(pp.p)
    Fp2 = pp.Fp2

    def fp2_pair(x):
        c = list(Fp2(x).list())
        a = int(c[0]) % p if len(c) > 0 else 0
        b = int(c[1]) % p if len(c) > 1 else 0
        return ["0x%x" % a, "0x%x" % b]

    def coords(c):
        return [fp2_pair(x) for x in c]

    def pt_xz(P):
        return [fp2_pair(P[0]), fp2_pair(P[2])]

    captured = []  # one record per KaniGluingIsogenyChainDim4Half built

    orig_init = kgmod.KaniGluingIsogenyChainDim4Half.__init__

    def patched_init(self, points_m, a1, a2, q, m, Theta12, M_product_dim2,
                     M_start_dim4, M_gluing_dim4, e4, dual=False, strategy_dim2=None):
        orig_init(self, points_m, a1, a2, q, m, Theta12, M_product_dim2,
                  M_start_dim4, M_gluing_dim4, e4, dual, strategy_dim2)
        T1s, T2s = self.Theta12._theta_structures
        # generic points to exercise montgomery_point_to_theta_point
        R1 = T1s.P + T1s.Q
        R2 = T2s.P + T2s.Q
        rec = {
            "dual": bool(dual),
            "T1_xz": pt_xz(T1s.P),
            "U1_xz": pt_xz(T2s.P),
            "dim1_null_1": coords(T1s.null_point().coords()),
            "dim1_null_2": coords(T2s.null_point().coords()),
            "theta12_null": coords(self.Theta12.zero().coords()),
            "dim2_codomain_null": coords(self._isogenies_dim2._codomain.null_point().coords()),
            "domain_product_null": coords(self.domain_product.null_point().coords()),
            "conv_pt_1_xz": pt_xz(R1),
            "conv_theta_1": coords(T1s(R1).coords()),
            "conv_pt_2_xz": pt_xz(R2),
            "conv_theta_2": coords(T2s(R2).coords()),
            "conv_prod": coords(self.Theta12.product_theta_point(
                [T1s(R1), T2s(R2)]).coords()),
        }
        captured.append(rec)

    kgmod.KaniGluingIsogenyChainDim4Half.__init__ = patched_init

    out = {"level": 1, "prime_decimal": str(p),
           "field": "F_{p^2}=F_p[i]/(i^2+1); a+b*i serialized as [hex(a),hex(b)]",
           "note": "dim-1 theta null = (X+Z, X-Z); product dim2 = t[k%2]*u[k//2]; "
                   "product dim2->4 = s1[k%4]*s2[k//4]. F1 then F2_dual per vector.",
           "vectors": []}

    for idx in range(args.n):
        captured.clear()
        v = SQIsignHD_verif(pp, idx)
        v.recover_pk_and_com()
        v.recover_chal()
        v.image_response()
        v.compute_HD()
        assert len(captured) == 2, "expected 2 half-chains, got %d" % len(captured)
        captured[0]["chain"] = "F1"
        captured[1]["chain"] = "F2_dual"
        out["vectors"].append({"index": idx, "half_chains": list(captured)})
        print("vector %d: captured %d half-chains" % (idx, len(captured)))

    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote %s (%d vectors)" % (out_path, len(out["vectors"])))


if __name__ == "__main__":
    main()
