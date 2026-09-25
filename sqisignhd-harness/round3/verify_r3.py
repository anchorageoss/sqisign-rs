#!/usr/bin/env sage -python
r"""
The SQIsignHD reference verifier (Pierrick-Dartois/SQISignHD-lib,
Verification/Verify.py with its Theta_dim4 submodule), unmodified, run at the
round-3 level I prime with the parameters of docs/COMPACT_R3.md:

    p = 3 * 2**324 - 1,  e (embedding exponent, the library's `f`) = 174,
    r = ceil(e/2) + 2 = 89,  lambda = 128,

and the two non-residue tables regenerated for this prime with the library's
own rule (gen_tables_r3.sage, seed 0). Nothing in the library is edited: a
subclass supplies the parameters and the table paths, and the vector files
are the library's text format (Public_keys / Signatures, one entry per 3 /
8 lines), written by our signer.

    sage -python verify_r3.py --lib /path/to/SQISignHD-lib \
        --pk Public_keys_r3lvl1.txt --sig Signatures_r3lvl1.txt [-i N | -n N]

Exit status 0 when every checked signature verifies (middle-codomain match
and the image check, as the library's `verify()` decides).
"""
import argparse, os, sys
from sage.all import ZZ, GF, ceil, log, BinaryQF, proof
proof.all(False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lib", default=os.environ.get("SQISIGNHD_LIB", "SQISignHD-lib"))
    ap.add_argument("--pk", required=True)
    ap.add_argument("--sig", required=True)
    ap.add_argument("-i", "--index", type=int)
    ap.add_argument("-n", "--n_samples", type=int, default=1)
    ap.add_argument("--expect-reject", action="store_true",
                    help="the vectors are negatives: succeed when every one is rejected")
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    verif_dir = os.path.join(os.path.abspath(args.lib), "Verification")
    sys.path.insert(0, verif_dir)
    os.chdir(verif_dir)   # the library reads its Data/ tables by relative path
    import Verify

    class SQIsignHD_R3(Verify.SQIsignHD):
        """The library's parameter object at the round-3 level I prime."""
        def __init__(self):
            self.p = 3 * 2**324 - 1
            self.c = 3
            self.e = 324                      # the 2-adic exponent of p + 1
            self.Fp2 = GF(self.p**2, 'i', modulus=[1, 0, 1], proof=False)
            self.i = self.Fp2.gen()
            self.NQR_TABLE = self.read_gf_table(os.path.join(here, "Data", "NQR_TABLE_r3lvl1.txt"))
            self.Z_NQR_TABLE = self.read_gf_table(os.path.join(here, "Data", "Z_NQR_TABLE_r3lvl1.txt"))
            for x in self.NQR_TABLE:
                assert not x.is_square()
            for x in self.Z_NQR_TABLE:
                assert x.is_square() and not (x - 1).is_square()
            self.lvl = 1
            self.lamb = 128
            self.n_bytes = 2 * self.lamb // 8
            self.f = 174                      # docs/COMPACT_R3.md Section 1 (ours)
            self.r = ceil(self.f / 2) + 2     # 89
            self.Q0 = BinaryQF([1, 0, 1])

    pp = SQIsignHD_R3()
    pk_file = os.path.abspath(os.path.join(here, args.pk)) if not os.path.isabs(args.pk) else args.pk
    sig_file = os.path.abspath(os.path.join(here, args.sig)) if not os.path.isabs(args.sig) else args.sig
    indices = [args.index] if args.index is not None else list(range(args.n_samples))
    ok = True
    for n in indices:
        v = Verify.SQIsignHD_verif.__new__(Verify.SQIsignHD_verif)
        v.params = pp
        v.A_pk, v.h_pk_P, v.h_pk_Q = pp.read_public_key(pk_file, n)
        v.A_com, v.a, v.b, v.c_or_d, v.q, v.h_com_P, v.h_com_Q, v.chal = pp.read_signature(sig_file, n)
        try:
            accepted = bool(v.verify(verbose=False))
        except Exception as exc:  # a malformed negative may raise inside the library
            accepted = False
            import traceback
            tb = traceback.extract_tb(exc.__traceback__)[-2:]
            where = "; ".join(f"{os.path.basename(f.filename)}:{f.lineno} {f.line}" for f in tb)
            print(f"vector {n}: rejected by exception: {type(exc).__name__}: {exc} [{where}]")
        print(f"vector {n}: {'accepted' if accepted else 'rejected'}")
        ok &= (not accepted) if args.expect_reject else accepted
    print("ALL " + ("REJECTED" if args.expect_reject else "ACCEPTED") if ok else "MISMATCH")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
