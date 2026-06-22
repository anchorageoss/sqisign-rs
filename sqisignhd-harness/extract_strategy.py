#!/usr/bin/env sage
r"""
extract_strategy.py  --  additive oracle for Phase 5b.6 (the dim-4
optimal-strategy chain loop).

The reference dim-4 half-chain `IsogenyChainDim4(B_K, gluing, e, m)` runs one
gluing step (k=0) then a sequence of plain `(2,2,2,2)`-steps (k=1..n-1, where
n=e-m), generating each plain step's order-8 kernel by an optimal-strategy walk
of doublings and pushforwards starting from `B_K`. This script captures, per
signature and per half-chain (F1, F2_dual), the inputs the Rust strategy loop
needs to *re-derive* those per-step kernels and reproduce the chain:

  * post_glue_basis: the four kernel generators on the gluing codomain
    (`[gluing(T) for T in B_K]`, dim-4 theta, full order 2^(n-1));
  * glue_codomain_null: the gluing codomain theta null point (the structure the
    basis lives on, == chain_codomains[0]);
  * e, m, n, n_plain (= n-1), and the reference strategy + its doubling /
    image-eval counts (the completion-independent structure to cross-check).

The per-step plain kernels and codomains the loop must reproduce are already in
`chain_vectors.json` (F1_kernels / F2_dual_kernels) and `test_vectors_l1.json`
(chain_codomains); this only adds the loop's starting point.

Captured by wrapping `IsogenyChainDim4.isogeny_chain` at runtime; no reference
file and no Phase 0 deliverable is modified.

Usage:  sage extract_strategy.py [--out strategy_vectors.json] [--n N] [--lib PATH]
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
    out_path = os.path.abspath(args.out or "strategy_vectors.json")

    LIB = find_lib(args.lib)
    VERIF = os.path.join(LIB, "Verification")
    os.chdir(VERIF)
    sys.path.insert(0, VERIF)

    import Verify
    from Verify import SQIsignHD, SQIsignHD_verif
    import Theta_dim4.Theta_dim4_sage.pkg.isogenies.isogeny_chain_dim4 as chmod
    import Theta_dim4.Theta_dim4_sage.pkg.isogenies.isogeny_dim4 as isogmod
    import Theta_dim4.Theta_dim4_sage.pkg.theta_structures.Theta_dim4 as thmod

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

    captured = []
    counters = {"doublings": 0, "image_calls": 0}

    # Count the reference's dim-4 doublings and plain image-evals.
    orig_di = thmod.ThetaPointDim4.double_iter
    def patched_di(self, n):
        counters["doublings"] += int(n)
        return orig_di(self, n)
    thmod.ThetaPointDim4.double_iter = patched_di

    orig_img = isogmod.IsogenyDim4.image
    def patched_img(self, P):
        counters["image_calls"] += 1
        return orig_img(self, P)
    isogmod.IsogenyDim4.image = patched_img

    orig_chain = chmod.IsogenyChainDim4.isogeny_chain
    def patched_chain(self, B_K, first_isogenies):
        # full-order kernel generators pushed onto the gluing codomain
        post = [first_isogenies(T) for T in B_K]
        glue_codomain = first_isogenies._codomain
        rec = {
            "e": int(self.e), "m": int(self.m), "n": int(self.e - self.m),
            "n_plain": int(self.e - self.m - 1),
            "strategy": [int(s) for s in self.strategy],
            "post_glue_basis": [theta(P.coords()) for P in post],
            "glue_codomain_null": theta(glue_codomain.null_point().coords()),
        }
        captured.append(rec)
        return orig_chain(self, B_K, first_isogenies)
    chmod.IsogenyChainDim4.isogeny_chain = patched_chain

    out = {"level": 1, "prime_decimal": str(p),
           "field": "F_{p^2}=F_p[i]/(i^2+1); a+b*i serialized as [hex(a),hex(b)]",
           "note": "post_glue_basis: 4 dim-4 theta generators (full order 2^(n-1)) "
                   "on the gluing codomain. Walking the plain optimal strategy from "
                   "it re-derives F1/F2_dual plain kernels (chain_vectors.json) and "
                   "codomains (test_vectors_l1.json chain_codomains[1:]).",
           "vectors": []}

    for idx in range(args.n):
        captured.clear()
        counters["doublings"] = 0
        counters["image_calls"] = 0
        v = SQIsignHD_verif(pp, idx)
        v.recover_pk_and_com()
        v.recover_chal()
        v.image_response()
        v.compute_HD()
        assert len(captured) == 2, "expected 2 half-chains, got %d" % len(captured)
        captured[0]["chain"] = "F1"
        captured[1]["chain"] = "F2_dual"
        out["vectors"].append({
            "index": idx,
            "ref_dim4_doublings": counters["doublings"],
            "ref_plain_image_calls": counters["image_calls"],
            "half_chains": list(captured),
        })
        print("vector %d: n=%s strat_len=%s doublings=%d images=%d"
              % (idx, [hc["n"] for hc in captured],
                 [len(hc["strategy"]) for hc in captured],
                 counters["doublings"], counters["image_calls"]))

    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote %s (%d vectors)" % (out_path, len(out["vectors"])))


if __name__ == "__main__":
    main()
