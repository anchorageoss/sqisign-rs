#!/usr/bin/env sage
r"""
extract_theta_arith.py  --  supplementary oracle for Phase 1 (theta4 crate).

The Phase 0 vectors (test_vectors_l1.json) capture isogeny *codomain* theta null
points, which validate projective equality, normalisation and the Hadamard
transform.  They do NOT contain (point, 2*point) pairs, so generic-point theta
*doubling* cannot be checked against them.  This script produces that missing
ground truth, without touching any Phase 0 deliverable.

It builds a genuine dimension-4 theta structure as the product
ProductThetaStructureDim1To4(E, E, E, E) of a real level-1 supersingular curve
E = E_pk (from public key 0), takes deterministic points on E^4, maps them into
the theta model, and records their doubles / quadruples computed by the
reference `ThetaPointDim4.double`.  It also records the 2-torsion action
`act_point` on the null point.

Usage:  sage extract_theta_arith.py [--out theta_arith_vectors.json] [--lib PATH]
Env  :  SQISIGNHD_LIB  (path to the SQISignHD-lib checkout)
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
    cands += [
        os.path.join(os.getcwd(), "SQISignHD-lib"),
        os.path.expanduser("~/SQISignHD-lib"),
        "/home/user/SQISignHD-lib",
    ]
    for c in cands:
        if c and os.path.isdir(os.path.join(c, "Verification", "Theta_dim4",
                                            "Theta_dim4_sage")):
            return c
    raise SystemExit("Could not locate SQISignHD-lib; set --lib or $SQISIGNHD_LIB")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--lib", default=None)
    args = ap.parse_args()
    out_path = os.path.abspath(args.out or "theta_arith_vectors.json")

    LIB = find_lib(args.lib)
    VERIF = os.path.join(LIB, "Verification")
    os.chdir(VERIF)
    sys.path.insert(0, VERIF)

    from sage.all import ZZ
    import Verify
    from Verify import SQIsignHD, SQIsignHD_verif
    from Theta_dim4.Theta_dim4_sage.pkg.theta_structures.Theta_dim1 import ThetaStructureDim1
    from Theta_dim4.Theta_dim4_sage.pkg.theta_structures.Theta_dim4 import ProductThetaStructureDim1To4
    from Theta_dim4.Theta_dim4_sage.pkg.theta_structures.Tuple_point import TuplePoint

    pp = SQIsignHD(1)
    p = int(pp.p)
    Fp2 = pp.Fp2

    def fp2_pair(x):
        c = list(Fp2(x).list())
        a = int(c[0]) % p if len(c) > 0 else 0
        b = int(c[1]) % p if len(c) > 1 else 0
        return ["0x%x" % a, "0x%x" % b]

    def theta(tp):
        co = list(tp.coords())
        assert len(co) == 16
        return [fp2_pair(c) for c in co]

    # A real level-1 supersingular curve and a 2^f-torsion basis.
    v = SQIsignHD_verif(pp, 0)
    v.recover_pk_and_com()
    E = v.E_pk
    P, Q = v.P_pk, v.Q_pk          # order 2^f

    # Dimension-1 theta structure (auto canonical 4-torsion basis), then the
    # product over E^4.
    Th1 = ThetaStructureDim1(E)
    Th4 = ProductThetaStructureDim1To4(Th1, Th1, Th1, Th1)

    suitable = bool(Th4.has_suitable_doubling())
    print("has_suitable_doubling:", suitable)

    # Deterministic points on E (scalar combinations of the basis), grouped into
    # tuple points on E^4.
    scalars = [
        (1, 0), (0, 1), (1, 1), (3, 5), (7, 2), (123, 456),
        (1000003, 7), (5, 1000003), (2, 3), (9, 8), (17, 31), (40, 41),
        (12345, 0), (0, 67890), (111, 222), (98765, 43210),
    ]
    pts = [ZZ(a) * P + ZZ(b) * Q for (a, b) in scalars]

    cases = []
    for idx in range(0, len(pts) - 3, 4):
        tp = TuplePoint(pts[idx], pts[idx + 1], pts[idx + 2], pts[idx + 3])
        thP = Th4(tp)
        d1 = thP.double()
        d2 = thP.double_iter(2)
        cases.append({
            "scalars": [scalars[idx + r] for r in range(4)],
            "P": theta(thP),
            "double": theta(d1),
            "quad": theta(d2),
        })
        print("  case %d: P -> 2P -> 4P captured" % len(cases))

    # 2-torsion action on the null point (act_null), a few (i, j) pairs.
    null = Th4.null_point()
    act_cases = []
    for (i, j) in [(0, 0), (1, 0), (0, 1), (3, 5), (10, 6), (15, 15)]:
        a = null.act_point(tuple((i >> b) & 1 for b in range(4)),
                           tuple((j >> b) & 1 for b in range(4)))
        act_cases.append({"i": i, "j": j, "act": [fp2_pair(c) for c in a]})

    import sage.version
    doc = {
        "purpose": "Phase 1 theta-arithmetic oracle: generic-point doubling and "
                   "act_point on a product dim-4 theta structure E_pk^4.",
        "level": 1,
        "prime_decimal": str(p),
        "field": "F_{p^2}=F_p[i]/(i^2+1); a+b*i serialized as [hex(a),hex(b)]",
        "sage_version": str(sage.version.version),
        "structure": "ProductThetaStructureDim1To4(E_pk, E_pk, E_pk, E_pk)",
        "suitable_doubling": suitable,
        "null_point": theta(null),
        "cases": cases,
        "act_null_cases": act_cases,
    }
    with open(out_path, "w") as fh:
        json.dump(doc, fh, indent=1)
    print("wrote %d doubling cases + %d act cases to %s"
          % (len(cases), len(act_cases), out_path))


if __name__ == "__main__":
    main()
