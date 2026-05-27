"""
Wrapper to run the C reference's precompute_endomorphism_action.sage
with deuring2d dependencies on the path, then output Rust source
instead of C source.
"""
proof.all(False)
import sys, os

level = int(sys.argv[1])
output_dir = sys.argv[2]

# Add deuring2d and its dependencies to path
deps_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "deps")
sys.path.insert(0, os.path.join(deps_dir, "deuring-2D"))
sys.path.insert(0, os.path.join(deps_dir, "two-isogenies", "Theta-SageMath"))
sys.path.insert(0, os.path.join(deps_dir, "deuring-2D", "qlapoti"))

# Change to the reference precomp level dir for parameters.py
ref_precomp_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'reference', 'src', 'precomp', 'ref', f'lvl{level}')
ref_precomp_dir = os.path.realpath(ref_precomp_dir)
os.chdir(ref_precomp_dir)

# Add the reference scripts dir to path
ref_scripts_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'reference', 'scripts', 'precomp')
ref_scripts_dir = os.path.realpath(ref_scripts_dir)
sys.path.insert(0, ref_scripts_dir)

from sage.misc.banner import require_version
if not require_version(10, 0, print_message=True):
    exit('')

from parameters import p, f
from torsion_basis import even_torsion_basis_E0
from maxorders import orders

print(f"Level {level}: p has {int(p).bit_length()} bits, f={f}")
print(f"Orders: q = {[q for q,_,_,_,_,_ in orders]}")

from sage.groups.generic import order_from_multiple
pari.allocatemem(1 << 34)

if p % 4 != 3:
    raise NotImplementedError('requires p = 3 (mod 4)')
assert (1 << f).divides(p + 1)

Fp2.<i> = GF((p,2), modulus=[1,0,1])
sqrtm1 = min(Fp2(-1).sqrt(all=True))

# Montgomery form constants (from cformat.py)
radix_map = {1: 51, 3: 55, 5: 57}
nwords_map = {1: 5, 3: 7, 5: 9}
radix = radix_map[level]
nwords = nwords_map[level]
n_limbs = 1 + floor(log(p, 2^radix))
R_mont = 2^(radix * ceil(log(p, 2^radix)))

def fp_to_mont_limbs(val):
    v = (ZZ(val) * R_mont) % p
    limbs = []
    for _ in range(n_limbs):
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

# The compute function from the reference script
def compute(q, mat, idl, iso1q):
    print(f'Computing q = {q}...')
    E0 = EllipticCurve(Fp2, [1,0])
    E0.set_order((p+1)^2)

    if q == 1:
        E1 = E0
        P1, Q1 = even_torsion_basis_E0(E1, f)
    else:
        Quat.<ii,jj,kk> = QuaternionAlgebra(-1, -p)
        I = Quat.ideal(map(Quat, idl))
        O0 = Quat.quaternion_order(list(map(Quat, orders[0][2])))
        O1 = I.right_order()
        assert I.left_order() == O0
        assert O0.is_maximal() and O1.is_maximal()
        assert I.norm() % 2

        from deuring2d import Deuring2D
        ctx = Deuring2D(p)
        assert ctx.O0.order == O0
        assert ctx.E0 == E0
        ctx.sqrtm1 = sqrtm1

        P0, Q0 = data[0][1]

        for deg in range(1,10):
            print(f'  trying deg = {deg}...')
            ctx.e = E0.cardinality(extension_degree=2^deg).sqrt().valuation(2) - 1

            Fbig.<U> = Fp2.extension(2^deg)
            E0_big = E0.change_ring(Fbig)
            ctx.E0 = E0_big
            ctx.F = Fbig
            ctx.P = P0.change_ring(Fbig)
            ctx.Q = Q0.change_ring(Fbig)
            # Re-derive iota and pi over the extended field so all
            # endomorphisms act on the same base ring as ctx.E0.
            ctx.iota = E0_big.automorphisms()[-1]
            ctx.pi = E0_big.frobenius_isogeny()
            assert ctx.e == ctx.E0.order().sqrt().valuation(2) - 1
            for _ in range(ctx.e - f):
                ctx.P = ctx.P.division_points(2)[0]
                ctx.Q = ctx.Q.division_points(2)[0]
            ctx.P.set_order(multiple=2^ctx.e)
            ctx.Q.set_order(multiple=2^ctx.e)

            try:
                E1, P1, Q1 = ctx.IdealToIsogeny(I)
                break
            except Deuring2D.Failure:
                continue
        else:
            raise NotImplementedError('Deuring2D failed')

        E1 = E1.change_ring(Fp2)
        j = GF(p)(E1.j_invariant())
        X = polygen(GF(p))
        for A,_ in sorted((256*(X^2-3)^3 - (X^2-4)*j).roots()):
            E1_ = EllipticCurve(Fp2, [0,A,0,1,0])
            try:
                iso = min(E1.isomorphisms(E1_))
                break
            except ValueError:
                pass
        E1 = iso.codomain()
        P1 = iso._eval(P1)
        Q1 = iso._eval(Q1)

        P1 *= ctx.P.order() // P0.order()
        Q1 *= ctx.Q.order() // Q0.order()
        P1 = P1.change_ring(Fp2)
        Q1 = Q1.change_ring(Fp2)
        P1.set_order(P0.order())
        Q1.set_order(Q0.order())
        assert P0.order() == Q0.order() == P1.order() == Q1.order() == 2^f
        assert P1.weil_pairing(Q1,2^f) == P0.weil_pairing(Q0,2^f)^I.norm()

    if q == 1:
        endo_i, = (a for a in E1.automorphisms() if a.scaling_factor() == sqrtm1)
    else:
        iso = E1.isomorphism(min(Fp2(-q).sqrt(all=True)), is_codomain=True)
        try:
            endo_i = iso * E1.isogeny(None, codomain=iso.domain(), degree=q)
        except ValueError:
            assert False

    endo_1 = E1.scalar_multiplication(1)
    endo_j = E1.frobenius_isogeny()
    endo_k = endo_i * endo_j

    denom = mat.denominator()
    coprime = denom.prime_to_m_part(lcm(P1.order(), Q1.order()))
    P1d, Q1d = (inverse_mod(coprime, T.order()) * T for T in (P1, Q1))
    denom //= coprime

    extdeg = next(d for d in range(1,denom+1) if ((denom<<f)^2).divides(E1.order(extension_degree=d)))
    if extdeg == 1:
        Fbig = Fp2
    else:
        Fbig.<U> = Fp2.extension(extdeg)

    P1d, Q1d = (T.change_ring(Fbig) for T in (P1d, Q1d))
    P1d.set_order(multiple=denom<<f)
    for l,m in denom.factor():
        for ii in range(m):
            assert l.divides(P1d.order())
            P1d = P1d.division_points(l)[0]
            P1d.set_order(multiple=denom<<f)
            for Q1d_ in Q1d.division_points(l):
                o = order_from_multiple(P1d.weil_pairing(Q1d_, P1d.order()), denom<<f, operation='*')
                if o == P1d.order():
                    Q1d = Q1d_
                    break
            else:
                assert False
    assert hasattr(P1d, '_order')
    Q1d.set_order(multiple=denom<<f)

    denom *= coprime

    PQ1d = P1d, Q1d
    mati = matrix(Zmod(1<<f), [endo_i._eval(T).log(PQ1d) for T in PQ1d])
    matj = matrix(Zmod(1<<f), [endo_j._eval(T).log(PQ1d) for T in PQ1d])
    matk = matj * mati

    gens = []
    for row in denom * mat:
        endo = sum(ZZ(c)*e for c,e in zip(row, (endo_1,endo_i,endo_j,endo_k)))
        gens.append(endo)
    gen1, gen2, gen3, gen4 = gens

    assert mat[0] == vector((1,0,0,0))
    mat2 = matrix(ZZ, [gen2._eval(T).log(PQ1d) for T in PQ1d]) / denom
    mat3 = matrix(ZZ, [gen3._eval(T).log(PQ1d) for T in PQ1d]) / denom
    mat4 = matrix(ZZ, [gen4._eval(T).log(PQ1d) for T in PQ1d]) / denom
    mat2, mat3, mat4 = (M.change_ring(Zmod(1<<f)) for M in (mat2,mat3,mat4))

    A = E1.a2()
    assert E1.a_invariants() == (0,A,0,1,0)

    return (A, (A+2)/4), (P1, Q1), (mati,matj,matk), (mat2,mat3,mat4)

# Run computation
todo = [(q, mat*iso1q, idl, iso1q) for q,iso1q,mat,_,idl,_ in orders]
data = [None] * len(todo)

assert todo[0][0] == 1
data[0] = compute(*todo[0])
print(f'[+] finished precomputation for q = {todo[0][0]}.')

# Compute the rest sequentially (parallel fails due to sys.path not being inherited by workers)
for idx, inp in enumerate(todo[1:], 1):
    data[idx] = compute(*inp)
    print(f'[+] finished precomputation for q = {inp[0]}.')

# Output Rust source
out_path = os.path.join(output_dir, 'endomorphism_action.rs')
with open(out_path, 'w') as out:
    out.write('//! Endomorphism action precomputed data: curves with endomorphism rings.\n')
    out.write('//!\n')
    out.write('//! Contains curves, torsion bases, and endomorphism action matrices\n')
    out.write('//! for the standard and alternate starting curves.\n')
    out.write('//!\n')
    out.write('//! Fp elements are stored as Montgomery-form limb arrays.\n')
    out.write('//! Generated from SageMath precompute scripts. DO NOT EDIT.\n\n')
    out.write('use num_bigint::BigInt;\n')
    out.write('use once_cell::sync::Lazy;\n\n')
    out.write('type Ibz = BigInt;\n\n')
    out.write('fn ibz(s: &str) -> Ibz {\n')
    out.write('    s.parse::<Ibz>().unwrap()\n')
    out.write('}\n\n')
    out.write(f'pub const NWORDS_FIELD: usize = {nwords};\n')
    out.write(f'pub const NUM_CURVES_WITH_ENDOMORPHISMS: usize = {len(data)};\n\n')

    for entry_idx, ((A, A24), (P1, Q1), (mati, matj, matk), (mat2, mat3, mat4)) in enumerate(data):
        out.write(f'// ==== Curve {entry_idx} ====\n\n')

        # Curve: A.re, A.im, C.re, C.im (C is always 1),
        #        A24.x.re, A24.x.im, A24.z.re, A24.z.im (A24 = (A+2)/4 : 1)
        A_re, A_im = get_fp2_coeffs(A)
        A24_re, A24_im = get_fp2_coeffs(A24)
        one_re, one_im = get_fp2_coeffs(1)

        fp_entries = [
            ('CURVE_A_RE', A_re), ('CURVE_A_IM', A_im),
            ('CURVE_C_RE', one_re), ('CURVE_C_IM', one_im),
            ('CURVE_A24_X_RE', A24_re), ('CURVE_A24_X_IM', A24_im),
            ('CURVE_A24_Z_RE', one_re), ('CURVE_A24_Z_IM', one_im),
        ]
        for name, val in fp_entries:
            limbs = fp_to_mont_limbs(val)
            out.write(f'pub const ENDOMORPHISM_{entry_idx}_{name}: [u64; {nwords}] = [{format_limbs(limbs)}];\n')
        out.write('\n')

        # Basis: P, Q, P-Q in projective (X:Z) form
        # Sage points use T[0]=x, T[2]=z in projective coordinates
        PmQ = P1 - Q1
        basis_points = [(P1, 'P'), (Q1, 'Q'), (PmQ, 'PMQ')]
        basis_entries = []
        for pt, prefix in basis_points:
            x_re, x_im = get_fp2_coeffs(pt[0])  # X coordinate
            z_re, z_im = get_fp2_coeffs(pt[2])   # Z coordinate
            basis_entries.extend([
                (f'BASIS_{prefix}_X_RE', x_re),
                (f'BASIS_{prefix}_X_IM', x_im),
                (f'BASIS_{prefix}_Z_RE', z_re),
                (f'BASIS_{prefix}_Z_IM', z_im),
            ])
        for name, val in basis_entries:
            limbs = fp_to_mont_limbs(val)
            out.write(f'pub const ENDOMORPHISM_{entry_idx}_{name}: [u64; {nwords}] = [{format_limbs(limbs)}];\n')
        out.write('\n')

        # Action matrices: i, j, k, gen2, gen3, gen4 (transposed/column-major)
        matrices = [
            ('ACTION_I', mati), ('ACTION_J', matj), ('ACTION_K', matk),
            ('ACTION_GEN2', mat2), ('ACTION_GEN3', mat3), ('ACTION_GEN4', mat4),
        ]
        for mat_name, mat_val in matrices:
            # Column-major: [0][0], [1][0], [0][1], [1][1]
            entries = [int(mat_val[0,0]), int(mat_val[1,0]), int(mat_val[0,1]), int(mat_val[1,1])]
            out.write(f'pub static ENDOMORPHISM_{entry_idx}_{mat_name}: [Lazy<Ibz>; 4] = [\n')
            for dec in entries:
                out.write(f'    Lazy::new(|| ibz("{dec}")),\n')
            out.write(f'];\n')
        out.write('\n')

print(f'\nwrote {out_path} ({len(data)} curves)')
