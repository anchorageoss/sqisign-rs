#!/usr/bin/env sage
r"""
extract_gluing.py  --  Phase 4 oracle for the dim-4 gluing isogeny.

Each half-chain begins with one dim-4 gluing isogeny (`GluingIsogenyDim4`),
whose codomain is the first codomain of the chain (`chain_codomains[0]` in
test_vectors_l1.json) and whose dual theta-null point has zero coordinates. This
script dumps, for both gluing objects of all 5 Level-1 signatures:

  * the (base-changed) product domain theta-null point,
  * the 5 kernel 8-torsion theta points `L_K_8` and their multi-index
    directions `L_K_8_ind` (= [1,2,4,8,3]),
  * the codomain theta-null point and its dual (with zeros),
  * a handful of `special_image` (P, translates, output) triples.

It also counts how often the reference's doubling base-change fallback
(`ThetaStructureDim4._arithmetic_base_change`) triggers, to confirm whether the
"change of theta coordinates" path is exercised at Level 1.

Captured by wrapping reference methods at runtime; no reference file and no
Phase 0 deliverable modified.

Usage:  sage extract_gluing.py [--out gluing_vectors.json] [--n NVEC] [--lib PATH]
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
    ap.add_argument("--max_image_pairs", type=int, default=6)
    args = ap.parse_args()
    out_path = os.path.abspath(args.out or "gluing_vectors.json")

    LIB = find_lib(args.lib)
    VERIF = os.path.join(LIB, "Verification")
    os.chdir(VERIF)
    sys.path.insert(0, VERIF)

    import Verify
    from Verify import SQIsignHD, SQIsignHD_verif
    import Theta_dim4.Theta_dim4_sage.pkg.isogenies.gluing_isogeny_dim4 as glumod
    import Theta_dim4.Theta_dim4_sage.pkg.theta_structures.Theta_dim4 as thmod
    from Theta_dim4.Theta_dim4_sage.pkg.theta_structures.theta_helpers_dim4 import multindex_to_index

    pp = SQIsignHD(1)
    p = int(pp.p)
    Fp2 = pp.Fp2

    def fp2_pair(x):
        c = list(Fp2(x).list())
        a = int(c[0]) % p if len(c) > 0 else 0
        b = int(c[1]) % p if len(c) > 1 else 0
        return ["0x%x" % a, "0x%x" % b]

    def theta(coords):
        co = list(coords)
        assert len(co) == 16
        return [fp2_pair(c) for c in co]

    gluings = []        # captured GluingIsogenyDim4 records (per construction)
    gluing_objs = []    # the corresponding GluingIsogenyDim4 instances (for identity)
    image_pairs = []    # captured special_image triples
    counters = {"base_change": 0, "special_image_calls": 0}

    # ---- wrap GluingIsogenyDim4.__init__ ----
    orig_init = glumod.GluingIsogenyDim4.__init__
    def patched_init(self, domain, L_K_8, L_K_8_ind, coerce=None):
        orig_init(self, domain, L_K_8, L_K_8_ind, coerce)
        gluing_objs.append(self)
        dual = list(self._codomain.null_point_dual())
        n_zeros = sum(1 for x in dual if x == 0)
        gluings.append({
            "domain_null": theta(domain.null_point().coords()),
            "L_K_8": [theta(k.coords()) for k in L_K_8],
            "L_K_8_ind": [int(multindex_to_index(ind)) for ind in L_K_8_ind],
            "codomain_null": theta(self._codomain.null_point().coords()),
            "codomain_dual_null": theta(dual),
            "dual_zero_count": int(n_zeros),
        })
    glumod.GluingIsogenyDim4.__init__ = patched_init

    # ---- wrap special_image (capture a few triples) ----
    orig_special = glumod.GluingIsogenyDim4.special_image
    def patched_special(self, P, L_trans, L_trans_ind):
        out = orig_special(self, P, L_trans, L_trans_ind)
        counters["special_image_calls"] += 1
        # Identify which gluing object this call belongs to (both gluings are
        # built before any special_image is called, so index-by-identity).
        gi = next(i for i, o in enumerate(gluing_objs) if o is self)
        if sum(1 for ip in image_pairs if ip["gluing_index"] == gi) < args.max_image_pairs:
            image_pairs.append({
                "gluing_index": gi,
                "P": theta(P.coords()),
                "L_trans": [theta(Q.coords()) for Q in L_trans],
                "L_trans_ind": [int(t) for t in L_trans_ind],
                "out": theta(out.coords()),
            })
        return out
    glumod.GluingIsogenyDim4.special_image = patched_special

    # ---- wrap the doubling base-change fallback to count triggers ----
    orig_bc = thmod.ThetaStructureDim4._arithmetic_base_change
    def patched_bc(self, max_iter=50):
        counters["base_change"] += 1
        return orig_bc(self, max_iter)
    thmod.ThetaStructureDim4._arithmetic_base_change = patched_bc

    all_vectors = []
    for idx in range(args.n):
        gluings.clear()
        gluing_objs.clear()
        image_pairs.clear()
        counters["base_change"] = 0
        counters["special_image_calls"] = 0

        v = SQIsignHD_verif(pp, idx)
        v.recover_pk_and_com()
        v.recover_chal()
        v.image_response()
        v.compute_HD()

        # Two gluing objects per signature (F1 then F2_dual half-chains).
        assert len(gluings) == 2, "expected 2 gluing objects, got %d" % len(gluings)
        for g, chain in zip(gluings, ("F1", "F2_dual")):
            g["chain"] = chain
        all_vectors.append({
            "index": idx,
            "gluings": list(gluings),
            "special_image_pairs": list(image_pairs),
            "ref_base_change_triggers": counters["base_change"],
            "ref_special_image_calls": counters["special_image_calls"],
        })
        print("vector %d: %d gluing objs, dual_zero_counts=%s, base_change_triggers=%d, "
              "special_image_calls=%d (captured %d pairs)"
              % (idx, len(gluings), [g["dual_zero_count"] for g in gluings],
                 counters["base_change"], counters["special_image_calls"], len(image_pairs)))

    import sage.version
    doc = {
        "purpose": "Phase 4 oracle: dim-4 gluing isogeny internals (domain, "
                   "5-point kernel + directions, codomain with zero dual null) "
                   "and special_image triples.",
        "level": 1,
        "prime_decimal": str(p),
        "field": "F_{p^2}=F_p[i]/(i^2+1); a+b*i serialized as [hex(a),hex(b)]",
        "sage_version": str(sage.version.version),
        "note": "L_K_8_ind directions are [1,2,4,8,3]; the 5th (3 = bits 0+1) is "
                "the diagonal edge that lets the spanning tree cover all 16 "
                "indices despite zero-denominator single-bit edges. "
                "ref_base_change_triggers counts ThetaStructureDim4 doubling "
                "base-change fallbacks (the 'change of theta coordinates' path).",
        "n_vectors": len(all_vectors),
        "vectors": all_vectors,
    }
    with open(out_path, "w") as fh:
        json.dump(doc, fh, indent=1)
    sz = os.path.getsize(out_path)
    total_bc = sum(v["ref_base_change_triggers"] for v in all_vectors)
    print("wrote %d vectors to %s (%.2f MB)" % (len(all_vectors), out_path, sz / 1e6))
    print("TOTAL doubling base-change triggers across all vectors: %d" % total_bc)


if __name__ == "__main__":
    main()
