proof.all(False)
import sys, os

level = int(sys.argv[1])
output_dir = sys.argv[2]

os.chdir(f"/home/user/sqisign-rs/reference/src/precomp/ref/lvl{level}")
sys.path.insert(0, "/home/user/sqisign-rs/reference/scripts/precomp")

from parameters import p, f

# The C reference uses non-standard radix per limb (not 64-bit full limbs).
# These are the "ref" radix values from cformat.py:
# Level 1 (p5248): radix=51, 5 limbs -> R = 2^(51*5) = 2^255
# Level 3 (p65376): radix=55, 7 limbs -> R = 2^(55*7) = 2^385
# Level 5 (p27500): radix=57, 9 limbs -> R = 2^(57*9) = 2^513
radix_map = {1: 51, 3: 55, 5: 57}
nwords_map = {1: 5, 3: 7, 5: 9}

radix = radix_map[level]
nwords = nwords_map[level]
n_limbs = 1 + floor(log(p, 2^radix))
R_mont = 2^(radix * ceil(log(p, 2^radix)))

print(f"Level {level}: radix={radix}, nwords={nwords}, n_limbs={n_limbs}")
print(f"R = 2^({radix} * {ceil(log(p, 2^radix))}) = 2^{radix * ceil(log(p, 2^radix))}")

Fp2.<i> = GF((p,2), modulus=[1,0,1])
E0 = EllipticCurve(Fp2, [1, 0])

from torsion_basis import even_torsion_basis_E0
P, Q = even_torsion_basis_E0(E0, f)

def fp_to_mont_limbs(val):
    """Convert an Fp element to Montgomery-form limbs using the reference's radix."""
    v = (ZZ(val) * R_mont) % p
    limbs = []
    for idx in range(n_limbs):
        limbs.append(int(v) & ((1 << radix) - 1))
        v >>= radix
    return limbs

def format_limbs(limbs):
    return ', '.join(f'0x{l:x}' for l in limbs)

def get_fp2_coeffs(el):
    coeffs = Fp2(el).polynomial().list()
    re = ZZ(coeffs[0]) if len(coeffs) > 0 else ZZ(0)
    im = ZZ(coeffs[1]) if len(coeffs) > 1 else ZZ(0)
    return re, im

Px_re, Px_im = get_fp2_coeffs(P.x())
Qx_re, Qx_im = get_fp2_coeffs(Q.x())

out_path = os.path.join(output_dir, 'e0_basis.rs')
with open(out_path, 'w') as out:
    out.write('//! E0 torsion basis point coordinates.\n')
    out.write('//!\n')
    out.write('//! Fp elements are stored as Montgomery-form limb arrays.\n')
    out.write('//! Generated from SageMath precompute scripts. DO NOT EDIT.\n\n')
    out.write(f'pub const NWORDS_FIELD: usize = {nwords};\n\n')

    for name, val in [('BASIS_E0_P_X_RE', Px_re), ('BASIS_E0_P_X_IM', Px_im),
                       ('BASIS_E0_Q_X_RE', Qx_re), ('BASIS_E0_Q_X_IM', Qx_im)]:
        limbs = fp_to_mont_limbs(val)
        out.write(f'pub const {name}: [u64; {nwords}] = [{format_limbs(limbs)}];\n')

print(f"  wrote {out_path}")

# Cross-check: verify our limbs match the C reference
print(f"\n  Sage computed limbs (verify against C non-broadwell path):")
for name, val in [('P_X_RE', Px_re), ('P_X_IM', Px_im), ('Q_X_RE', Qx_re), ('Q_X_IM', Qx_im)]:
    limbs = fp_to_mont_limbs(val)
    print(f"    {name}: {{{format_limbs(limbs)}}}")
