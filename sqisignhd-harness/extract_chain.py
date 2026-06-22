#!/usr/bin/env sage
r"""
extract_chain.py  --  per-chain oracle for Phase 3 (dim-4 isogeny chain).

For each Level-1 signature this dumps, in chain order, the kernel `K_8` of every
*plain* `IsogenyDim4` step of both half-chains F1 and F2_dual. Driving
`IsogenyDim4::from_kernel` over these reproduces the full sequence of per-step
codomains (validated against `test_vectors_l1.json`'s `chain_codomains`) and,
via the last codomains, the middle-codomain match.

It also counts the work the *reference* chain performs to derive those kernels -
the dimension-4 point doublings (`double_iter`) and plain isogeny image
evaluations (`IsogenyDim4.image`) - so the Rust side can report an honest
full-chain timing estimate (those steps are the optimal-strategy bookkeeping,
not implemented in Phase 3).

The gluing step (index 0 of each half-chain) overrides `__init__`/uses
`special_image`, so it is neither captured nor counted here (Phase 4).

Captured at runtime by wrapping reference methods; no reference file and no
Phase 0/2 deliverable is modified.

Usage:  sage extract_chain.py [--out chain_vectors.json] [--n NVEC] [--lib PATH]
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
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--lib", default=None)
    args = ap.parse_args()
    out_path = os.path.abspath(args.out or "chain_vectors.json")

    LIB = find_lib(args.lib)
    VERIF = os.path.join(LIB, "Verification")
    os.chdir(VERIF)
    sys.path.insert(0, VERIF)

    import Verify
    from Verify import SQIsignHD, SQIsignHD_verif
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

    # ---- instrumentation ----
    kernels = []        # ordered K_8 of each plain step (F1 then F2_dual)
    counters = {"doublings": 0, "image_calls": 0}

    orig_init = isogmod.IsogenyDim4.__init__
    def patched_init(self, domain, K_8, codomain=None, precomputation=None):
        orig_init(self, domain, K_8, codomain, precomputation)
        if K_8 is not None:
            kernels.append([theta(k.coords()) for k in K_8])
    isogmod.IsogenyDim4.__init__ = patched_init

    orig_image = isogmod.IsogenyDim4.image
    def patched_image(self, P):
        counters["image_calls"] += 1
        return orig_image(self, P)
    isogmod.IsogenyDim4.image = patched_image

    orig_double_iter = thmod.ThetaPointDim4.double_iter
    def patched_double_iter(self, n):
        counters["doublings"] += int(n)
        return orig_double_iter(self, n)
    thmod.ThetaPointDim4.double_iter = patched_double_iter

    all_vectors = []
    for idx in range(args.n):
        kernels.clear()
        counters["doublings"] = 0
        counters["image_calls"] = 0

        v = SQIsignHD_verif(pp, idx)
        v.recover_pk_and_com()
        v.recover_chal()
        v.image_response()
        v.compute_HD()

        n1 = len(v.F.F1._isogenies)
        n2 = len(v.F.F2_dual._isogenies)
        f1_plain = n1 - 1          # plain steps in F1 (index 0 is gluing)
        f2_plain = n2 - 1
        assert len(kernels) == f1_plain + f2_plain, \
            "captured %d, expected %d" % (len(kernels), f1_plain + f2_plain)

        all_vectors.append({
            "index": idx,
            "n1": n1, "n2": n2,
            "f1_plain_steps": f1_plain,
            "f2_plain_steps": f2_plain,
            "ref_dim4_doublings": counters["doublings"],
            "ref_plain_image_calls": counters["image_calls"],
            "F1_kernels": kernels[:f1_plain],
            "F2_dual_kernels": kernels[f1_plain:],
        })
        print("vector %d: F1 %d plain + F2_dual %d plain steps; "
              "ref dim4 doublings=%d, plain image-evals=%d"
              % (idx, f1_plain, f2_plain, counters["doublings"],
                 counters["image_calls"]))

    import sage.version
    doc = {
        "purpose": "Phase 3 oracle: ordered per-step kernels (K_8) for the plain "
                   "IsogenyDim4 steps of both half-chains, plus reference "
                   "operation counts.",
        "level": 1,
        "prime_decimal": str(p),
        "field": "F_{p^2}=F_p[i]/(i^2+1); a+b*i serialized as [hex(a),hex(b)]",
        "sage_version": str(sage.version.version),
        "note": "Plain steps only (gluing is index 0 of each half-chain, Phase 4). "
                "ref_dim4_doublings / ref_plain_image_calls count the dimension-4 "
                "double_iter steps and IsogenyDim4.image calls the reference makes "
                "to build both half-chains (optimal-strategy bookkeeping), for the "
                "Rust full-chain timing estimate.",
        "n_vectors": len(all_vectors),
        "vectors": all_vectors,
    }
    with open(out_path, "w") as fh:
        json.dump(doc, fh, indent=1)
    sz = os.path.getsize(out_path)
    print("wrote %d vectors to %s (%.1f MB)" % (len(all_vectors), out_path, sz / 1e6))


if __name__ == "__main__":
    main()
