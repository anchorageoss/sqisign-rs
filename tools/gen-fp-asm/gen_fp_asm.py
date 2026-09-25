#!/usr/bin/env python3
"""Emit the x86-64 field backends `crates/prism-verify/src/fp/<prime>_asm.rs`
(inline assembly with BMI2 and ADX, or portable kernels of the same
schedule, chosen at run time by `fp::dispatch`) for the SQIsign round-3 primes
`p = c 2^e - 1`: saturated 64-bit limbs, Montgomery form with `R = 2^(64 n)`,
the multiplication schedule of the reference's `broadwell` backend (one
`mulx` row per word of `b` accumulated with two carry chains, one fold per
word because `p + 1 = K 2^(64 f)` has a single non-zero limb and `-p^-1 = 1
mod 2^64`), fused `Fp2` products, and the rest in Rust.

Register plans, as the reference's: 6 limbs fit r8-r14 + r15/rax with the
three pointers in rsi/rcx/rdi; 8 limbs borrow rbx and rbp (pushed); 11 limbs
park both operands and the output pointer on the stack and use rsi, rcx and
rdi as accumulators too.

Also emits `<dir>/../mont_asm.rs`, the variable-modulus kernels of the
Miller-Rabin test (6 to 11 words).

Usage: gen_fp_asm.py [--check] [--crate-path <prefix>] [--without-mont] <dir>
(prism-rs: default dir crates/prism-verify/src/fp, prefix `crate`; sqisign-rs:
`--crate-path crate::v3 --without-mont crates/verify/src/v3/fp`, no
Miller-Rabin kernels since that crate has no hash-to-prime)
"""
import os
import sys

PRIMES = [
    # name, c, e, limbs, encoded bytes, origin
    ("p324_3", 3, 324, 6, 41, "SQIsign round 3, NIST level I"),
    ("p500_27", 27, 500, 8, 64, "SQIsign round 3, NIST level III"),
    ("p664_17", 17, 664, 11, 84, "SQIsign round 3, NIST level V"),
]

GP = ["r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15"]

# The path of the field module's parent crate root in `use` lines of the
# generated code: `crate` in prism-rs, `crate::v3` in sqisign-rs.
CRATE = "crate"


def limbs(x, n):
    return [(x >> (64 * i)) & (2**64 - 1) for i in range(n)]


def rs_array(name, x, n, doc):
    body = ", ".join(f"0x{l:016x}" for l in limbs(x, n))
    return f"/// {doc}\nconst {name}: [u64; {n}] = [{body}];\n"


def d32(r):
    """The 32-bit alias, for zeroing with `xor`."""
    m = {"rax": "eax", "rbx": "ebx", "rcx": "ecx", "rdx": "edx", "rsi": "esi", "rdi": "edi", "rbp": "ebp"}
    return m.get(r, r + "d")


class Plan:
    """Registers and operand addressing for `n` limbs."""

    def __init__(self, n, fold_idx):
        self.n = n
        self.fold_idx = fold_idx
        if n <= 7:
            self.acc = GP[: n + 1]
            self.t0, self.t1 = GP[n + 1] if n + 1 < 8 else "rax", "rax"
            if n + 1 < 8:
                self.t0, self.t1 = GP[n + 1], "rax"
            self.saved = []
            self.stack = False
        elif n == 8:
            self.acc = GP[:8] + ["rax"]
            self.t0, self.t1 = "rbx", "rbp"
            self.saved = ["rbx", "rbp"]
            self.stack = False
        else:
            self.acc = GP[:8] + ["rax", "rbx", "rbp", "rsi"]
            assert len(self.acc) == n + 1
            self.t0, self.t1 = "rcx", "rdi"
            self.saved = ["rbx", "rbp"]
            self.stack = True

    # operand addressing --------------------------------------------------
    def a(self, i, fp2=False, part=0):
        off = 8 * (i + part * self.n)
        return f"qword ptr [rsp + {off}]" if self.stack else f"qword ptr [rsi + {off}]"

    def b(self, i, fp2=False, part=0):
        if self.stack:
            base = 8 * (2 * self.n if fp2 else self.n)
            return f"qword ptr [rsp + {base + 8 * (i + part * self.n)}]"
        return f"qword ptr [rcx + {8 * (i + part * self.n)}]"

    def scratch(self, i):
        """The `2p - b_im` buffer of the fused real part: in the frame when
        the operands are parked, else the caller's scratch array in rdi."""
        if self.stack:
            return f"qword ptr [rsp + {8 * (4 * self.n + i)}]"
        return f"qword ptr [rdi + {8 * i}]"

    def frame(self, fp2):
        return 8 * (5 * self.n) if fp2 else 8 * (2 * self.n)


def emit_asm(plan, K_sym, fp2_kind, p2_needed):
    """Lines of one asm! body. fp2_kind: None (Fp mul), 're' or 'im'."""
    n = plan.n
    Z = plan.acc[:]
    T0, T1 = plan.t0, plan.t1
    out = []

    def e(s):
        out.append(f'        "{s}",')

    # prologue
    for r in plan.saved:
        e(f"push {r}")
    if plan.stack:
        fp2 = fp2_kind is not None
        e(f"sub rsp, {plan.frame(fp2)}")
        # copy the operands (rax is free here)
        width = 2 * n if fp2 else n
        for i in range(width):
            e(f"mov rax, qword ptr [rsi + {8 * i}]")
            e(f"mov qword ptr [rsp + {8 * i}], rax")
        for i in range(width):
            e(f"mov rax, qword ptr [rcx + {8 * i}]")
            e(f"mov qword ptr [rsp + {8 * (width + i)}], rax")
    # operand selectors for the rows
    if fp2_kind == "re":
        # z = a_re b_re[i] + a_im (2p - b_im)[i]
        first = (lambda i: plan.b(i, True, 0), lambda k: plan.a(k, True, 0))
        second = (lambda i: plan.scratch(i), lambda k: plan.a(k, True, 1))
        # 2p - b_im into the scratch, with the first n+1 accumulators as temps
        e(f"mov {Z[0]}, qword ptr [rip + {{p2}}]")
        e(f"mov {Z[1]}, qword ptr [rip + {{p2}} + 8]")
        for k in range(2, n - 1):
            e(f"mov {Z[k]}, {Z[1]}")
        e(f"mov {Z[n - 1]}, qword ptr [rip + {{p2}} + {8 * (n - 1)}]")
        e(f"sub {Z[0]}, {plan.b(0, True, 1)}")
        for k in range(1, n):
            e(f"sbb {Z[k]}, {plan.b(k, True, 1)}")
        for k in range(n):
            e(f"mov {plan.scratch(k)}, {Z[k]}")
    elif fp2_kind == "im":
        first = (lambda i: plan.b(i, True, 1), lambda k: plan.a(k, True, 0))
        second = (lambda i: plan.b(i, True, 0), lambda k: plan.a(k, True, 1))
    else:
        first = (lambda i: plan.b(i), lambda k: plan.a(k))
        second = None

    def row0(bsel, asel):
        e(f"mov rdx, {bsel(0)}")
        e(f"mulx {Z[1]}, {Z[0]}, {asel(0)}")
        e(f"xor {d32(Z[n])}, {d32(Z[n])}")
        for k in range(1, n):
            e(f"mulx {Z[k + 1]}, {T1}, {asel(k)}")
            e(f"adcx {Z[k]}, {T1}")
        e(f"adc {Z[n]}, 0")

    def muladd(z, i, bsel, asel, clear):
        e(f"mov rdx, {bsel(i)}")
        e(f"xor {d32(clear)}, {d32(clear)}")
        e(f"mulx {T0}, {T1}, {asel(0)}")
        e(f"adox {z[0]}, {T1}")
        e(f"adox {z[1]}, {T0}")
        for k in range(1, n):
            e(f"mulx {T0}, {T1}, {asel(k)}")
            e(f"adcx {z[k]}, {T1}")
            e(f"adox {z[k + 1]}, {T0}")
        e(f"adc {z[n]}, 0")

    def fold(z):
        e(f"mov rdx, {z[0]}")
        e(f"xor {d32(T0)}, {d32(T0)}")
        e(f"mulx {T0}, {T1}, qword ptr [rip + {{pp1}}]")
        e(f"adox {z[plan.fold_idx]}, {T1}")
        e(f"adox {z[plan.fold_idx + 1]}, {T0}")

    z = Z[:]
    row0(*first)
    if second:
        muladd(z, 0, second[0], second[1], T0)
    fold(z)
    for i in range(1, n):
        l = z[1:] + [z[0]]
        muladd(l, i, first[0], first[1], l[n])
        if second:
            muladd(l, i, second[0], second[1], T0)
        fold(l)
        z = l
    res = z[1:]
    # epilogue: the result stays in registers (the compiler stores it where
    # it is needed and keeps no zeroed temporary); registers Rust reserves
    # (rbx, rbp) are moved into free operand-capable ones first
    if plan.stack:
        e(f"add rsp, {plan.frame(fp2_kind is not None)}")
    free = [r for r in ["rdx", T0, T1, z[0]] if r not in ("rbx", "rbp") and r not in res]
    outs = []
    for r in res:
        if r in ("rbx", "rbp"):
            dst = free.pop(0)
            e(f"mov {dst}, {r}")
            outs.append(dst)
        else:
            outs.append(r)
    for r in reversed(plan.saved):
        e(f"pop {r}")
    return "\n".join(out), outs


# ---------------------------------------------------------------------------
# Addition, subtraction and the multiplication's lazy operands, with the
# result in registers (P26 item 1, second attempt: with the result behind a
# pointer the compiler's copies around the block cost what the block saved).
# Baseline x86-64 only, no dispatch. `p`'s limbs below the top are all ones,
# so a masked `p` is the mask itself there.
# ---------------------------------------------------------------------------


class AddPlan:
    def __init__(self, n):
        pool = GP + ["rdi", "rdx", "rax"]
        self.n = n
        self.z = pool[:n]
        temps = [r for r in ["rax", "rdx"] if r not in self.z]
        self.saved = []
        for r in ["rbx", "rbp"]:
            if len(temps) < 2:
                temps.append(r)
                self.saved.append(r)
        self.t0, self.t1 = temps[0], temps[1]


def emit_add_asm(plan, kind, top_shift, p_top, p2):
    """`kind`: 'add' (two conditional `-p`), 'sub' (`+p` on borrow, `+p`
    when negative), 'add_lazy' (plain sum), 'sub_lazy_2p' (`a - b + 2p`)."""
    n, z, T0, T1 = plan.n, plan.z, plan.t0, plan.t1
    top = z[n - 1]
    out = []

    def e(s):
        out.append(f'        "{s}",')

    for r in plan.saved:
        e(f"push {r}")
    for k in range(n):
        e(f"mov {z[k]}, qword ptr [rsi + {8 * k}]")
    if kind in ("add", "add_lazy"):
        e(f"add {z[0]}, qword ptr [rcx]")
        for k in range(1, n):
            e(f"adc {z[k]}, qword ptr [rcx + {8 * k}]")
        if kind == "add":
            for _ in range(2):
                e(f"mov {T0}, {top}")
                e(f"shr {T0}, {top_shift}")
                e(f"neg {T0}")
                e(f"mov {T1}, 0x{p_top:x}")
                e(f"and {T1}, {T0}")
                e(f"sub {z[0]}, {T0}")
                for k in range(1, n - 1):
                    e(f"sbb {z[k]}, {T0}")
                e(f"sbb {top}, {T1}")
    else:
        if kind == "sub":
            e(f"xor {d32(T0)}, {d32(T0)}")
        e(f"sub {z[0]}, qword ptr [rcx]")
        for k in range(1, n):
            e(f"sbb {z[k]}, qword ptr [rcx + {8 * k}]")
        if kind == "sub":
            e(f"sbb {T0}, 0")
            for round_ in range(2):
                if round_ == 1:
                    e(f"mov {T0}, {top}")
                    e(f"sar {T0}, {top_shift}")
                e(f"mov {T1}, 0x{p_top:x}")
                e(f"and {T1}, {T0}")
                e(f"add {z[0]}, {T0}")
                for k in range(1, n - 1):
                    e(f"adc {z[k]}, {T0}")
                e(f"adc {top}, {T1}")
        else:
            # + 2p: limb 0 is -2, the middle limbs -1, the top 2K - 1
            e(f"add {z[0]}, -2")
            for k in range(1, n - 1):
                e(f"adc {z[k]}, -1")
            e(f"mov {T0}, 0x{p2[n - 1]:x}")
            e(f"adc {top}, {T0}")
    for r in reversed(plan.saved):
        e(f"pop {r}")
    return "\n".join(out)


def add_operands(plan):
    lines = ['        in("rsi") a,', '        in("rcx") b,']
    for i, r in enumerate(plan.z):
        lines.append(f'        lateout("{r}") c{i},')
    for r in (plan.t0, plan.t1):
        if r not in ("rbx", "rbp"):
            lines.append(f'        out("{r}") _,')
    lines.append("        options()," if plan.saved else "        options(nostack),")
    return "\n".join(lines)


def operands(plan, with_p2, outs, fp2_kind):
    """The operand list of an asm! block for this plan: the pointer inputs
    (a, b, and the scratch buffer of the fused real part when the operands
    are not parked), the result registers as outputs, the rest clobbered."""
    lines = ["        pp1 = sym P_PLUS_1_TOP,"]
    if with_p2:
        lines.append("        p2 = sym P2,")
    ptrs = {"rsi": "a", "rcx": "b"}
    if fp2_kind == "re" and not plan.stack:
        ptrs["rdi"] = "scratch"
    for i, r in enumerate(outs):
        if r in ptrs:
            lines.append(f'        inout("{r}") {ptrs[r]} => c{i},')
        else:
            lines.append(f'        lateout("{r}") c{i},')
    for r in ["rsi", "rcx", "rdi", "rax", "rdx"] + GP:
        if r in outs:
            continue
        if r in ptrs:
            lines.append(f'        inout("{r}") {ptrs[r]} => _,' if plan.stack else f'        in("{r}") {ptrs[r]},')
        elif r == "rdi" and fp2_kind != "re" and not plan.stack:
            continue  # untouched (it is a temporary in the parked plans)
        else:
            lines.append(f'        out("{r}") _,')
    if plan.saved or plan.stack:
        lines.append("        options(),")
    else:
        lines.append("        options(nostack),")
    return "\n".join(lines)


def result_decls(n):
    return "\n".join(f"    let c{i}: u64;" for i in range(n))


def result_array(n):
    return "[" + ", ".join(f"c{i}" for i in range(n)) + "]"



def generate(name, c, e, n, enc, origin):
    p = c * 2**e - 1
    R = 2 ** (64 * n)
    bits = p.bit_length()
    fold_idx = e // 64
    K = c << (e - 64 * fold_idx)
    assert p + 1 == K * 2 ** (64 * fold_idx) and K < 2**64
    assert pow(-p, -1, 2**64) == 1
    assert fold_idx + 1 == n, "the fold must end at the top accumulator"
    top_shift = bits - 64 * (n - 1)
    chunk = enc - 1
    plan = Plan(n, fold_idx)
    consts = (
        rs_array("P", p, n, "`p`.")
        + rs_array("ONE", R % p, n, "`R mod p`, the Montgomery form of 1.")
        + rs_array("R2", R * R % p, n, "`R^2 mod p`, the constant that takes an integer to Montgomery form.")
        + rs_array("THREE_INV", pow(3, -1, p) * R % p, n, "`3^-1` in Montgomery form.")
        + rs_array("CHUNK_R", (2 ** (8 * chunk)) * R % p, n, f"`2^{8 * chunk}` in Montgomery form: one {chunk}-byte chunk of `decode_reduce`.")
    )
    p2 = ", ".join(f"0x{l:016x}" for l in limbs(2 * p, n))
    upper = name.upper()
    cm1 = c - 1
    mul_body, mul_outs = emit_asm(plan, "pp1", None, False)
    re_body, re_outs = emit_asm(plan, "pp1", "re", True)
    im_body, im_outs = emit_asm(plan, "pp1", "im", False)
    mul_ops = operands(plan, False, mul_outs, None)
    aplan = AddPlan(n)
    p_top = limbs(p, n)[n - 1]
    assert all(l == 2**64 - 1 for l in limbs(p, n)[: n - 1]), "the limbs of p below the top must be all ones"
    p2_limbs = limbs(2 * p, n)
    assert p2_limbs[0] == 2**64 - 2 and all(l == 2**64 - 1 for l in p2_limbs[1 : n - 1])
    add_kernels = ""
    for kind, doc in [
        ("add", "`a + b`, below `2^BITS` for inputs below `2^BITS`: the sum, then twice `-p` when the result is at or above `2^BITS`"),
        ("sub", "`a - b`, non-negative and below `2^BITS` for inputs below `2^BITS`: `+p` on borrow, `+p` again when still negative"),
        ("add_lazy", "unreduced `a + b`, an input of the multiplication only"),
        ("sub_lazy_2p", "unreduced `a - b + 2p` (positive), an input of the multiplication only"),
    ]:
        add_kernels += f'''/// {doc}. The result leaves in registers.
///
/// # Safety
/// Live `NLIMBS`-limb arrays.
#[inline(always)]
unsafe fn {kind}_asm(a: *const u64, b: *const u64) -> [u64; NLIMBS] {{
{result_decls(n)}
    core::arch::asm!(
{emit_add_asm(aplan, kind, top_shift, p_top, p2_limbs)}
{add_operands(aplan)}
    );
    {result_array(n)}
}}

'''
    re_ops = operands(plan, True, re_outs, "re")
    im_ops = operands(plan, False, im_outs, "im")
    decls = result_decls(n)
    result = result_array(n)
    saved_note = (
        "borrows `rbx` and `rbp` (pushed and popped inside the block)" if n == 8
        else "parks both operands and the output pointer on the stack and accumulates in `rsi`, `rcx`, `rdi` as well" if n > 8
        else "keeps everything in `r8`-`r15`, `rax` and the three pointer registers"
    )
    return f'''//! `{name}` (`p = {c} * 2^{e} - 1`, {origin}): the x86-64 backend, saturated
//! 64-bit limbs with the Montgomery multiplication in inline assembly
//! (`mulx`, `adcx`, `adox`) when the CPU has BMI2 and ADX and in portable
//! Rust otherwise, chosen once at run time (`fp::dispatch`, P26 item 0; the
//! assembly was P25's opt-in `asm` feature).
//!
//! GENERATED by `tools/gen-fp-asm/gen_fp_asm.py`. Do not edit; change the
//! generator.
//!
//! Representation: {n} little-endian 64-bit limbs, Montgomery form with
//! `R = 2^{64 * n}`, every stored value below `2^{bits}` ({64 * n - bits} spare bits).
//! Addition subtracts `p` conditionally twice and returns below `2^{bits}`;
//! subtraction adds `p` conditionally twice and returns non-negative and
//! below `2^{bits}`; the product of two values below `2^{bits + 1}` is below
//! `p + 2^{2 * bits + 2 - 64 * n}`. Comparisons and the encoding reduce to `[0, p)` first.
//!
//! The multiplication is the schedule of the SQIsign round-3 reference's
//! `broadwell` backend (`src/gf/broadwell/{name}/fp_asm.S`, Apache-2.0):
//! one row of `mulx` per word of `b` accumulated with the two carry chains
//! (`adcx`/`adox`), and after each row one Montgomery fold, a single
//! multiplication because `p + 1 = 0x{K:x} * 2^{64 * fold_idx}` has one non-zero
//! limb and `-p^-1 = 1 mod 2^64`. This size {saved_note}. The fused `Fp2`
//! product computes `a_re b_re + a_im (2p - b_im)` and `a_re b_im + a_im
//! b_re` with one reduction each; squaring uses two products of unreduced
//! sums. The exponentiations use `(p - 3) / 4 = {cm1} 2^{e - 2} + (2^{e - 2} - 1)`:
//! the all-ones chain for `2^{e - 2} - 1` and a short tail.
//!
//! The assembly needs BMI2 and ADX (Intel Broadwell or AMD Zen and later);
//! `fp::dispatch` reads CPUID once (in `no_std`, through `core::arch`) and
//! every kernel below branches on the cached answer, a relaxed atomic load.
//! The portable kernels of this module share the representation, so the
//! two paths agree bit for bit; `tests/fp_backends_differential.rs`
//! compares both against big-integer arithmetic, and the known-answer tests
//! run under both (`PRISM_FORCE_PORTABLE=1`). On other architectures the
//! generated radix backend (`{name}.rs`) is the implementation.

use super::dispatch;
use super::{{Fp2, FpBackend}};
use {CRATE}::params::{upper};
use hybrid_array::Array;
use subtle::Choice;

const NLIMBS: usize = {n};
type Limbs = Array<u64, <{upper} as {CRATE}::params::Prime>::FpLimbs>;

{consts}/// `2p`, for the fused multiplication's `2p - b_im` and the lazy subtraction.
static P2: [u64; {n}] = [{p2}];
/// The one non-zero limb of `p + 1`, at index {fold_idx}.
static P_PLUS_1_TOP: u64 = 0x{K:x};
/// Its index: the fold adds `m P_PLUS_1_TOP` at this limb.
const FOLD: usize = {fold_idx};
/// Canonical encoding length.
const ENCODED_BYTES: usize = {enc};
const CHUNK_BYTES: usize = ENCODED_BYTES - 1;
/// `e - 2`: the all-ones part of `(p - 3) / 4` has this many bits.
const ONES_BITS: u32 = {e - 2};
/// `c - 1`, the small exponent of the tail of `(p - 3) / 4`.
const TAIL: u64 = {cm1};

#[inline]
fn as_arr(a: &Limbs) -> &[u64; NLIMBS] {{
    <&[u64; NLIMBS]>::try_from(&a[..]).expect("invariant: FpLimbs == NLIMBS")
}}

#[inline]
fn as_arr_mut(a: &mut Limbs) -> &mut [u64; NLIMBS] {{
    <&mut [u64; NLIMBS]>::try_from(&mut a[..]).expect("invariant: FpLimbs == NLIMBS")
}}

/// `a b / R mod p`, the assembly kernel. The result leaves the block in registers and
/// is stored into `c` by the compiler and drops the
/// zeroed temporary an output pointer would have forced it to keep (a
/// third of the multiplication's cost, measured, P26 item 2).
///
/// # Safety
/// Live `NLIMBS`-limb arrays; the CPU has BMI2 and ADX
/// (`dispatch::asm_available`).
#[inline(always)]
unsafe fn mont_mul_asm(a: *const u64, b: *const u64) -> [u64; NLIMBS] {{
{decls}
    core::arch::asm!(
{mul_body}
{mul_ops}
    );
    {result}
}}

/// `(a_re b_re - a_im b_im) / R`: the real part of an `Fp2` product, inputs
/// as `2 NLIMBS` contiguous limbs (`re` then `im`), one reduction.
/// `scratch` is an `NLIMBS`-limb buffer for `2p - b_im` (unused when the
/// kernel parks its operands in its own frame).
///
/// # Safety
/// As [`mont_mul_asm`], with `2 NLIMBS`-limb inputs and a live scratch.
#[inline(always)]
unsafe fn fp2_mul_re_asm(a: *const u64, b: *const u64, scratch: *mut u64) -> [u64; NLIMBS] {{
    let _ = &scratch;
{decls}
    core::arch::asm!(
{re_body}
{re_ops}
    );
    {result}
}}

/// `(a_re b_im + a_im b_re) / R`: the imaginary part, one reduction.
///
/// # Safety
/// As [`mont_mul_asm`], with `2 NLIMBS`-limb inputs.
#[inline(always)]
unsafe fn fp2_mul_im_asm(a: *const u64, b: *const u64) -> [u64; NLIMBS] {{
{decls}
    core::arch::asm!(
{im_body}
{im_ops}
    );
    {result}
}}

/// `c <- a b / R mod p` in portable Rust: operand scanning with `u128`
/// products and the same one-limb fold per word as the assembly
/// (`-p^-1 = 1 mod 2^64`, so the Montgomery multiplier is the low limb
/// itself). Inputs below `2^(BITS + 1)`, output below `p + 2^(2 BITS + 2 - 64 NLIMBS)`.
#[inline]
fn mont_mul_portable(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS], b: &[u64; NLIMBS]) {{
    let mut t = [0u64; NLIMBS + 2];
    for bi in b.iter() {{
        let mut carry = 0u128;
        for j in 0..NLIMBS {{
            let p = (t[j] as u128) + (a[j] as u128) * (*bi as u128) + carry;
            t[j] = p as u64;
            carry = p >> 64;
        }}
        let s = (t[NLIMBS] as u128) + carry;
        t[NLIMBS] = s as u64;
        t[NLIMBS + 1] = (s >> 64) as u64;
        fold_and_shift(&mut t);
    }}
    debug_assert_eq!(t[NLIMBS], 0);
    c.copy_from_slice(&t[..NLIMBS]);
}}

/// `c <- (a0 b0 + a1 b1) / R mod p` in portable Rust, the schedule of the
/// fused assembly kernels: per word `i`, the two rows `a0 b0[i]` and
/// `a1 b1[i]`, then one fold. Bit for bit the assembly's result.
#[inline]
fn mont_mul2_portable(
    c: &mut [u64; NLIMBS],
    a0: &[u64; NLIMBS],
    b0: &[u64; NLIMBS],
    a1: &[u64; NLIMBS],
    b1: &[u64; NLIMBS],
) {{
    let mut t = [0u64; NLIMBS + 2];
    for i in 0..NLIMBS {{
        for (a, bi) in [(a0, b0[i]), (a1, b1[i])] {{
            let mut carry = 0u128;
            for j in 0..NLIMBS {{
                let p = (t[j] as u128) + (a[j] as u128) * (bi as u128) + carry;
                t[j] = p as u64;
                carry = p >> 64;
            }}
            let s = (t[NLIMBS] as u128) + carry;
            t[NLIMBS] = s as u64;
            t[NLIMBS + 1] = t[NLIMBS + 1].wrapping_add((s >> 64) as u64);
        }}
        fold_and_shift(&mut t);
    }}
    debug_assert_eq!(t[NLIMBS], 0);
    c.copy_from_slice(&t[..NLIMBS]);
}}

/// One Montgomery fold of the accumulator and the division by `2^64`:
/// `m = t[0]`, `t <- (t + m p) / 2^64`, where `m p = m K 2^(64 FOLD) - m`
/// and `t[0] - m = 0`.
#[inline(always)]
fn fold_and_shift(t: &mut [u64; NLIMBS + 2]) {{
    let m = t[0];
    let prod = (m as u128) * (P_PLUS_1_TOP as u128);
    let s = (t[FOLD] as u128) + (prod as u64 as u128);
    t[FOLD] = s as u64;
    let s = (t[FOLD + 1] as u128) + ((prod >> 64) as u64 as u128) + (s >> 64);
    t[FOLD + 1] = s as u64;
    let mut carry = s >> 64;
    for w in t.iter_mut().skip(FOLD + 2) {{
        let s = (*w as u128) + carry;
        *w = s as u64;
        carry = s >> 64;
    }}
    for j in 0..NLIMBS + 1 {{
        t[j] = t[j + 1];
    }}
    t[NLIMBS + 1] = 0;
}}

/// The portable real part: `(a_re b_re + a_im (2p - b_im)) / R`.
#[inline]
fn fp2_mul_re_portable(
    c: &mut [u64; NLIMBS],
    a_re: &[u64; NLIMBS],
    a_im: &[u64; NLIMBS],
    b_re: &[u64; NLIMBS],
    b_im: &[u64; NLIMBS],
) {{
    let mut s = [0u64; NLIMBS];
    let mut borrow = 0u64;
    for i in 0..NLIMBS {{
        let (r, bw) = sbb(P2[i], b_im[i], borrow);
        s[i] = r;
        borrow = bw;
    }}
    mont_mul2_portable(c, a_re, b_re, a_im, &s);
}}

#[inline(always)]
fn mont_mul(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS], b: &[u64; NLIMBS]) {{
    if dispatch::asm_available() {{
        // SAFETY: live arrays; BMI2 and ADX detected.
        *c = unsafe {{ mont_mul_asm(a.as_ptr(), b.as_ptr()) }};
    }} else {{
        mont_mul_portable(c, a, b)
    }}
}}

/// The fused product on two `Fp2` values in memory (`re` then `im`, the
/// `repr(C)` layout of [`{CRATE}::fp::Fp2`]): no copies of the operands.
#[inline(always)]
fn fp2_mul_contiguous(out_re: &mut [u64; NLIMBS], out_im: &mut [u64; NLIMBS], a: *const u64, b: *const u64) {{
    let mut scratch = [0u64; NLIMBS];
    // SAFETY: the caller passes `2 NLIMBS` live limbs per operand; BMI2 and
    // ADX detected by the caller.
    unsafe {{
        *out_re = fp2_mul_re_asm(a, b, scratch.as_mut_ptr());
        *out_im = fp2_mul_im_asm(a, b);
    }}
}}

#[inline]
fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {{
    let t = (a as u128) + (b as u128) + (carry as u128);
    (t as u64, (t >> 64) as u64)
}}

#[inline]
fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {{
    let t = (a as u128).wrapping_sub(b as u128).wrapping_sub(borrow as u128);
    (t as u64, (t >> 127) as u64)
}}

#[inline]
fn sub_p_masked(a: &mut [u64; NLIMBS], mask: u64) {{
    let mut borrow = 0u64;
    for i in 0..NLIMBS {{
        let (r, b) = sbb(a[i], P[i] & mask, borrow);
        a[i] = r;
        borrow = b;
    }}
}}

{add_kernels}/// `c <- a + b`, below `2^BITS` for inputs below `2^BITS`.
#[inline(always)]
fn fp_add(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS], b: &[u64; NLIMBS]) {{
    // SAFETY: live arrays of NLIMBS limbs.
    *c = unsafe {{ add_asm(a.as_ptr(), b.as_ptr()) }};
}}

/// `c <- a - b`, non-negative and below `2^BITS` for inputs below `2^BITS`.
#[inline(always)]
fn fp_sub(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS], b: &[u64; NLIMBS]) {{
    // SAFETY: live arrays of NLIMBS limbs.
    *c = unsafe {{ sub_asm(a.as_ptr(), b.as_ptr()) }};
}}

/// Unreduced `a + b`, an input of the multiplication only.
#[inline(always)]
fn add_lazy(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS], b: &[u64; NLIMBS]) {{
    // SAFETY: live arrays of NLIMBS limbs.
    *c = unsafe {{ add_lazy_asm(a.as_ptr(), b.as_ptr()) }};
}}

/// Unreduced `a - b + 2p` (positive), an input of the multiplication only.
#[inline(always)]
fn sub_lazy_2p(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS], b: &[u64; NLIMBS]) {{
    // SAFETY: live arrays of NLIMBS limbs.
    *c = unsafe {{ sub_lazy_2p_asm(a.as_ptr(), b.as_ptr()) }};
}}

/// All ones when `a >= b` as `NLIMBS`-limb integers.
#[inline]
fn ge_mask(a: &[u64; NLIMBS], b: &[u64; NLIMBS]) -> u64 {{
    let mut borrow = 0u64;
    for i in 0..NLIMBS {{
        let (_, bw) = sbb(a[i], b[i], borrow);
        borrow = bw;
    }}
    borrow.wrapping_sub(1)
}}

/// Fully reduce a stored value (below `2^BITS < 2p`) to `[0, p)`.
#[inline]
fn canonical(a: &[u64; NLIMBS]) -> [u64; NLIMBS] {{
    let mut c = *a;
    let m = ge_mask(&c, &P);
    sub_p_masked(&mut c, m);
    let m = ge_mask(&c, &P);
    sub_p_masked(&mut c, m);
    c
}}

#[inline]
fn is_zero_limbs(a: &[u64; NLIMBS]) -> Choice {{
    let mut acc = 0u64;
    for &w in a {{
        acc |= w;
    }}
    Choice::from((((acc | acc.wrapping_neg()) >> 63) ^ 1) as u8)
}}

#[inline]
fn eq_limbs(a: &[u64; NLIMBS], b: &[u64; NLIMBS]) -> Choice {{
    let mut acc = 0u64;
    for i in 0..NLIMBS {{
        acc |= a[i] ^ b[i];
    }}
    Choice::from((((acc | acc.wrapping_neg()) >> 63) ^ 1) as u8)
}}

#[inline]
fn sqr_n(x: &mut [u64; NLIMBS], n: u32) {{
    for _ in 0..n {{
        let t = *x;
        mont_mul(x, &t, &t);
    }}
}}

/// `w^(2^k - 1)` by the binary all-ones chain: about `k` squarings and
/// `2 log2 k` multiplications.
fn ones(w: &[u64; NLIMBS], k: u32) -> [u64; NLIMBS] {{
    let mut acc = *w;
    let mut cur = 1u32;
    for bit in (0..(31 - k.leading_zeros())).rev() {{
        // acc = w^(2^cur - 1) -> w^(2^(2 cur) - 1)
        let prev = acc;
        sqr_n(&mut acc, cur);
        let t = acc;
        mont_mul(&mut acc, &t, &prev);
        cur *= 2;
        if (k >> bit) & 1 == 1 {{
            let t = acc;
            let mut sq = [0u64; NLIMBS];
            mont_mul(&mut sq, &t, &t);
            mont_mul(&mut acc, &sq, w);
            cur += 1;
        }}
    }}
    debug_assert_eq!(cur, k);
    acc
}}

/// `z <- w^((p-3)/4)` with `(p-3)/4 = TAIL 2^ONES_BITS + (2^ONES_BITS - 1)`.
fn exp3div4(z: &mut [u64; NLIMBS], w: &[u64; NLIMBS]) {{
    let t = ones(w, ONES_BITS);
    // u = w^(2^ONES_BITS) = t w; v = u^TAIL; z = v t
    let mut u = [0u64; NLIMBS];
    mont_mul(&mut u, &t, w);
    let mut v = ONE;
    for bit in (0..(64 - TAIL.leading_zeros())).rev() {{
        let s = v;
        mont_mul(&mut v, &s, &s);
        if (TAIL >> bit) & 1 == 1 {{
            let s = v;
            mont_mul(&mut v, &s, &u);
        }}
    }}
    mont_mul(z, &v, &t);
}}

/// Little-endian bytes to limbs, `Choice(1)` iff the value is below `p`;
/// the value itself (not yet in Montgomery form), zero on failure.
#[inline]
fn bytes_to_int(out: &mut [u64; NLIMBS], bytes: &[u8]) -> Choice {{
    *out = [0; NLIMBS];
    let mut buf = [0u8; 8 * NLIMBS];
    buf[..ENCODED_BYTES].copy_from_slice(&bytes[..ENCODED_BYTES]);
    for i in 0..NLIMBS {{
        out[i] = u64::from_le_bytes(buf[8 * i..8 * i + 8].try_into().expect("8 bytes"));
    }}
    let lt = Choice::from(((!ge_mask(out, &P)) & 1) as u8);
    let keep = 0u64.wrapping_sub(lt.unwrap_u8() as u64);
    for w in out.iter_mut() {{
        *w &= keep;
    }}
    lt
}}

impl FpBackend for {upper} {{
    #[inline]
    fn set_zero(out: &mut Limbs) {{
        *as_arr_mut(out) = [0; NLIMBS];
    }}

    #[inline]
    fn set_one(out: &mut Limbs) {{
        *as_arr_mut(out) = ONE;
    }}

    #[inline]
    fn set_small(out: &mut Limbs, val: u64) {{
        let mut v = [0u64; NLIMBS];
        v[0] = val;
        mont_mul(as_arr_mut(out), &v, &R2);
    }}

    #[inline]
    fn is_equal(a: &Limbs, b: &Limbs) -> Choice {{
        eq_limbs(&canonical(as_arr(a)), &canonical(as_arr(b)))
    }}

    #[inline]
    fn is_zero(a: &Limbs) -> Choice {{
        is_zero_limbs(&canonical(as_arr(a)))
    }}

    #[inline]
    fn copy(out: &mut Limbs, a: &Limbs) {{
        *as_arr_mut(out) = *as_arr(a);
    }}

    #[inline]
    fn add(out: &mut Limbs, a: &Limbs, b: &Limbs) {{
        fp_add(as_arr_mut(out), as_arr(a), as_arr(b));
    }}

    #[inline]
    fn sub(out: &mut Limbs, a: &Limbs, b: &Limbs) {{
        fp_sub(as_arr_mut(out), as_arr(a), as_arr(b));
    }}

    #[inline]
    fn neg(out: &mut Limbs, a: &Limbs) {{
        fp_sub(as_arr_mut(out), &[0; NLIMBS], as_arr(a));
    }}

    #[inline]
    fn mul(out: &mut Limbs, a: &Limbs, b: &Limbs) {{
        mont_mul(as_arr_mut(out), as_arr(a), as_arr(b));
    }}

    #[inline]
    fn sqr(out: &mut Limbs, a: &Limbs) {{
        let a = as_arr(a);
        mont_mul(as_arr_mut(out), a, a);
    }}

    #[inline]
    fn fp2_mul(
        out_re: &mut Limbs,
        out_im: &mut Limbs,
        a_re: &Limbs,
        a_im: &Limbs,
        b_re: &Limbs,
        b_im: &Limbs,
    ) {{
        if dispatch::asm_available() {{
            let mut a = [0u64; 2 * NLIMBS];
            let mut b = [0u64; 2 * NLIMBS];
            a[..NLIMBS].copy_from_slice(&a_re[..]);
            a[NLIMBS..].copy_from_slice(&a_im[..]);
            b[..NLIMBS].copy_from_slice(&b_re[..]);
            b[NLIMBS..].copy_from_slice(&b_im[..]);
            fp2_mul_contiguous(as_arr_mut(out_re), as_arr_mut(out_im), a.as_ptr(), b.as_ptr());
        }} else {{
            let (a_re, a_im, b_re, b_im) = (as_arr(a_re), as_arr(a_im), as_arr(b_re), as_arr(b_im));
            fp2_mul_re_portable(as_arr_mut(out_re), a_re, a_im, b_re, b_im);
            mont_mul2_portable(as_arr_mut(out_im), a_re, b_im, a_im, b_re);
        }}
    }}

    /// `a <- a b` with the product stored straight into `a` (P28): the
    /// kernel reads `a` through a raw pointer while it runs and the store
    /// follows it, so no temporary is copied.
    #[inline]
    fn mul_assign(a: &mut Limbs, b: &Limbs) {{
        if dispatch::asm_available() {{
            let dst = as_arr_mut(a);
            let pa = dst.as_mut_ptr() as *const u64;
            // SAFETY: live arrays; `pa` is read by the kernel before `dst`
            // is written; BMI2 and ADX detected.
            let r = unsafe {{ mont_mul_asm(pa, b.as_ptr()) }};
            *dst = r;
        }} else {{
            let mut t = [0u64; NLIMBS];
            mont_mul_portable(&mut t, as_arr(a), as_arr(b));
            *as_arr_mut(a) = t;
        }}
    }}

    #[inline]
    fn sqr_assign(a: &mut Limbs) {{
        if dispatch::asm_available() {{
            let dst = as_arr_mut(a);
            let pa = dst.as_mut_ptr() as *const u64;
            // SAFETY: as `mul_assign`.
            let r = unsafe {{ mont_mul_asm(pa, pa) }};
            *dst = r;
        }} else {{
            let mut t = [0u64; NLIMBS];
            let x = *as_arr(a);
            mont_mul_portable(&mut t, &x, &x);
            *as_arr_mut(a) = t;
        }}
    }}

    #[inline]
    fn add_assign(a: &mut Limbs, b: &Limbs) {{
        let dst = as_arr_mut(a);
        let pa = dst.as_mut_ptr() as *const u64;
        // SAFETY: live arrays; the kernel reads `pa` before `dst` is written.
        let r = unsafe {{ add_asm(pa, b.as_ptr()) }};
        *dst = r;
    }}

    #[inline]
    fn sub_assign(a: &mut Limbs, b: &Limbs) {{
        let dst = as_arr_mut(a);
        let pa = dst.as_mut_ptr() as *const u64;
        // SAFETY: as `add_assign`.
        let r = unsafe {{ sub_asm(pa, b.as_ptr()) }};
        *dst = r;
    }}

    #[inline]
    fn fp2_mul_assign(a: &mut Fp2<Self>, b: &Fp2<Self>) {{
        if dispatch::asm_available() {{
            let pa = a as *mut Fp2<Self> as *mut u64;
            let pb = b as *const Fp2<Self> as *const u64;
            let mut scratch = [0u64; NLIMBS];
            // SAFETY: `2 NLIMBS` live limbs behind each pointer (the layout
            // assertion of `fp2_mul_pair`); both kernels read `a` before
            // either half is stored; BMI2 and ADX detected.
            unsafe {{
                let re = fp2_mul_re_asm(pa as *const u64, pb, scratch.as_mut_ptr());
                let im = fp2_mul_im_asm(pa as *const u64, pb);
                *as_arr_mut(&mut a.re.limbs) = re;
                *as_arr_mut(&mut a.im.limbs) = im;
            }}
        }} else {{
            let (a_re, a_im) = (*as_arr(&a.re.limbs), *as_arr(&a.im.limbs));
            let (b_re, b_im) = (as_arr(&b.re.limbs), as_arr(&b.im.limbs));
            fp2_mul_re_portable(as_arr_mut(&mut a.re.limbs), &a_re, &a_im, b_re, b_im);
            mont_mul2_portable(as_arr_mut(&mut a.im.limbs), &a_re, b_im, &a_im, b_re);
        }}
    }}

    #[inline]
    fn fp2_sqr_assign(a: &mut Fp2<Self>) {{
        let (a_re, a_im) = (*as_arr(&a.re.limbs), *as_arr(&a.im.limbs));
        let mut sum = [0u64; NLIMBS];
        let mut diff = [0u64; NLIMBS];
        add_lazy(&mut sum, &a_re, &a_im);
        sub_lazy_2p(&mut diff, &a_re, &a_im);
        mont_mul(as_arr_mut(&mut a.re.limbs), &sum, &diff);
        let mut twice = [0u64; NLIMBS];
        add_lazy(&mut twice, &a_re, &a_re);
        mont_mul(as_arr_mut(&mut a.im.limbs), &twice, &a_im);
    }}

    #[inline]
    fn fp2_mul_pair(out: &mut Fp2<Self>, a: &Fp2<Self>, b: &Fp2<Self>) {{
        const _: () = assert!(core::mem::size_of::<Fp2<{upper}>>() == 16 * NLIMBS);
        if dispatch::asm_available() {{
            // `Fp2` is `repr(C)` of two `repr(transparent)` limb arrays: `2
            // NLIMBS` contiguous limbs, `re` first (the assertion above)
            let pa = a as *const Fp2<Self> as *const u64;
            let pb = b as *const Fp2<Self> as *const u64;
            let (re, im) = (&mut out.re.limbs, &mut out.im.limbs);
            fp2_mul_contiguous(as_arr_mut(re), as_arr_mut(im), pa, pb);
        }} else {{
            let (a_re, a_im) = (as_arr(&a.re.limbs), as_arr(&a.im.limbs));
            let (b_re, b_im) = (as_arr(&b.re.limbs), as_arr(&b.im.limbs));
            fp2_mul_re_portable(as_arr_mut(&mut out.re.limbs), a_re, a_im, b_re, b_im);
            mont_mul2_portable(as_arr_mut(&mut out.im.limbs), a_re, b_im, a_im, b_re);
        }}
    }}

    #[inline]
    fn fp2_sqr(out_re: &mut Limbs, out_im: &mut Limbs, a_re: &Limbs, a_im: &Limbs) {{
        let (a_re, a_im) = (as_arr(a_re), as_arr(a_im));
        let mut sum = [0u64; NLIMBS];
        let mut diff = [0u64; NLIMBS];
        add_lazy(&mut sum, a_re, a_im);
        sub_lazy_2p(&mut diff, a_re, a_im);
        mont_mul(as_arr_mut(out_re), &sum, &diff);
        let mut twice = [0u64; NLIMBS];
        add_lazy(&mut twice, a_re, a_re);
        mont_mul(as_arr_mut(out_im), &twice, a_im);
    }}

    #[inline]
    fn inv(out: &mut Limbs, a: &Limbs) {{
        let a = *as_arr(a);
        let mut z = [0u64; NLIMBS];
        exp3div4(&mut z, &a);
        let mut z2 = [0u64; NLIMBS];
        mont_mul(&mut z2, &z, &z);
        let mut z4 = [0u64; NLIMBS];
        mont_mul(&mut z4, &z2, &z2);
        mont_mul(as_arr_mut(out), &a, &z4);
    }}

    #[inline]
    fn sqrt(out: &mut Limbs, a: &Limbs) {{
        let a = *as_arr(a);
        let mut z = [0u64; NLIMBS];
        exp3div4(&mut z, &a);
        mont_mul(as_arr_mut(out), &a, &z);
    }}

    #[inline]
    fn is_square(a: &Limbs) -> Choice {{
        let a = *as_arr(a);
        let mut z = [0u64; NLIMBS];
        exp3div4(&mut z, &a);
        let mut z2 = [0u64; NLIMBS];
        mont_mul(&mut z2, &z, &z);
        let mut s = [0u64; NLIMBS];
        mont_mul(&mut s, &a, &z2);
        eq_limbs(&canonical(&s), &ONE) | is_zero_limbs(&canonical(&a))
    }}

    #[inline]
    fn half(out: &mut Limbs, a: &Limbs) {{
        let mut c = *as_arr(a);
        let odd = 0u64.wrapping_sub(c[0] & 1);
        let mut carry = 0u64;
        for i in 0..NLIMBS {{
            let (r, cy) = adc(c[i], P[i] & odd, carry);
            c[i] = r;
            carry = cy;
        }}
        let out = as_arr_mut(out);
        for i in 0..NLIMBS - 1 {{
            out[i] = (c[i] >> 1) | (c[i + 1] << 63);
        }}
        out[NLIMBS - 1] = (c[NLIMBS - 1] >> 1) | (carry << 63);
    }}

    #[inline]
    fn div3(out: &mut Limbs, a: &Limbs) {{
        mont_mul(as_arr_mut(out), as_arr(a), &THREE_INV);
    }}

    #[inline]
    fn exp3div4(out: &mut Limbs, a: &Limbs) {{
        let a = *as_arr(a);
        exp3div4(as_arr_mut(out), &a);
    }}

    #[inline]
    fn mul_small(out: &mut Limbs, a: &Limbs, val: u32) {{
        let mut v = [0u64; NLIMBS];
        let mut vi = [0u64; NLIMBS];
        vi[0] = val as u64;
        mont_mul(&mut v, &vi, &R2);
        mont_mul(as_arr_mut(out), as_arr(a), &v);
    }}

    #[inline]
    fn encode(out: &mut [u8], a: &Limbs) {{
        let mut one = [0u64; NLIMBS];
        one[0] = 1;
        let mut x = [0u64; NLIMBS];
        mont_mul(&mut x, as_arr(a), &one);
        let x = canonical(&x);
        let mut buf = [0u8; 8 * NLIMBS];
        for i in 0..NLIMBS {{
            buf[8 * i..8 * i + 8].copy_from_slice(&x[i].to_le_bytes());
        }}
        out[..ENCODED_BYTES].copy_from_slice(&buf[..ENCODED_BYTES]);
    }}

    #[inline]
    fn decode(out: &mut Limbs, bytes: &[u8]) -> Choice {{
        let mut x = [0u64; NLIMBS];
        let ok = bytes_to_int(&mut x, bytes);
        mont_mul(as_arr_mut(out), &x, &R2);
        ok
    }}

    #[inline]
    fn decode_reduce(out: &mut Limbs, bytes: &[u8]) {{
        let out = as_arr_mut(out);
        *out = [0; NLIMBS];
        let mut len = bytes.len();
        let rem = len % CHUNK_BYTES;
        if rem != 0 {{
            let mut tmp = [0u8; ENCODED_BYTES];
            tmp[..rem].copy_from_slice(&bytes[len - rem..]);
            let mut x = [0u64; NLIMBS];
            let _ = bytes_to_int(&mut x, &tmp);
            mont_mul(out, &x, &R2);
            len -= rem;
        }}
        while len > 0 {{
            len -= CHUNK_BYTES;
            let acc = *out;
            mont_mul(out, &acc, &CHUNK_R);
            let mut tmp = [0u8; ENCODED_BYTES];
            tmp[..CHUNK_BYTES].copy_from_slice(&bytes[len..len + CHUNK_BYTES]);
            let mut x = [0u64; NLIMBS];
            let _ = bytes_to_int(&mut x, &tmp);
            let mut chunk = [0u64; NLIMBS];
            mont_mul(&mut chunk, &x, &R2);
            let acc = *out;
            fp_add(out, &acc, &chunk);
        }}
    }}

    #[inline]
    fn cswap(a: &mut Limbs, b: &mut Limbs, ctl: Choice) {{
        let m = 0u64.wrapping_sub(ctl.unwrap_u8() as u64);
        let (a, b) = (as_arr_mut(a), as_arr_mut(b));
        for i in 0..NLIMBS {{
            let t = m & (a[i] ^ b[i]);
            a[i] ^= t;
            b[i] ^= t;
        }}
    }}

    #[inline]
    fn select(out: &mut Limbs, a0: &Limbs, a1: &Limbs, ctl: Choice) {{
        let cw = 0u64.wrapping_sub(ctl.unwrap_u8() as u64);
        let (a0, a1, out) = (as_arr(a0), as_arr(a1), as_arr_mut(out));
        for i in 0..NLIMBS {{
            out[i] = a0[i] ^ (cw & (a0[i] ^ a1[i]));
        }}
    }}
}}
'''


# ---------------------------------------------------------------------------
# Variable-modulus Montgomery multiplication for the Miller-Rabin test of
# the hash (P26 item 7): CIOS with `mulx`/`adcx`/`adox`, one kernel per word
# count 6..11 (the degrees `d` of the two profiles have 6 to 11 words). The
# modulus is passed as `NW` words followed by `-n^-1 mod 2^64`; the Rust side
# guarantees `bitlen(n) <= 64 NW - 2` and `a, b < n`, so every accumulator
# fits `NW + 1` words and the top word never overflows (`z < 2n <
# 2^(64 NW - 1)`, `a b_i < 2^(64 NW + 62)`, `m n < 2^(64 NW + 62)`).
# ---------------------------------------------------------------------------

MONT_WORDS = list(range(6, 12))


class MontPlan:
    """Registers and addressing for the `NW`-word variable-modulus kernel.

    Operands arrive in rsi (a), rcx (b), rdi (c) and rdx (n, moved to rbx
    behind a push because rbx cannot be an asm operand); accumulators are
    r8-r15 first, then rdi, rsi, rcx, rbx as the word count grows, each freed
    by parking the pointer it held (c) or copying the operand it addressed
    (a, b, n) to the stack.
    """

    def __init__(self, nw):
        self.nw = nw
        acc_extra = ["rdi", "rsi", "rcx", "rbx"][: max(0, nw + 1 - 8)]
        self.acc = GP[: min(nw + 1, 8)] + acc_extra
        assert len(self.acc) == nw + 1
        if nw == 6:
            self.t0, self.t1 = "r15", "rax"
            self.saved = ["rbx"]
        else:
            self.t0, self.t1 = "rax", "rbp"
            self.saved = ["rbx", "rbp"]
        self.park_c = "rdi" in acc_extra
        self.copy_a = "rsi" in acc_extra
        self.copy_b = "rcx" in acc_extra
        self.copy_n = "rbx" in acc_extra
        # frame: [c] [a...] [b...] [n..., ninv]
        off = 0
        self.off_c = off
        off += 8 if self.park_c else 0
        self.off_a = off
        off += 8 * nw if self.copy_a else 0
        self.off_b = off
        off += 8 * nw if self.copy_b else 0
        self.off_n = off
        off += 8 * (nw + 1) if self.copy_n else 0
        self.frame = off
        # keep the stack 16-byte aligned modulo the pushes (not required by
        # the kernel, harmless)

    def a(self, k):
        return f"qword ptr [rsp + {self.off_a + 8 * k}]" if self.copy_a else f"qword ptr [rsi + {8 * k}]"

    def b(self, i):
        return f"qword ptr [rsp + {self.off_b + 8 * i}]" if self.copy_b else f"qword ptr [rcx + {8 * i}]"

    def n(self, k):
        return f"qword ptr [rsp + {self.off_n + 8 * k}]" if self.copy_n else f"qword ptr [rbx + {8 * k}]"

    def ninv(self):
        return self.n(self.nw)


def emit_mont_asm(plan):
    """Lines of the asm! body of the `NW`-word variable-modulus kernel."""
    nw = plan.nw
    Z = plan.acc[:]
    T0, T1 = plan.t0, plan.t1
    out = []

    def e(s):
        out.append(f'        "{s}",')

    for r in plan.saved:
        e(f"push {r}")
    e("mov rbx, rdx")
    if plan.frame:
        e(f"sub rsp, {plan.frame}")
        if plan.park_c:
            e(f"mov qword ptr [rsp + {plan.off_c}], rdi")
        for flag, src, off, count in [
            (plan.copy_a, "rsi", plan.off_a, nw),
            (plan.copy_b, "rcx", plan.off_b, nw),
            (plan.copy_n, "rbx", plan.off_n, nw + 1),
        ]:
            if flag:
                for k in range(count):
                    e(f"mov rax, qword ptr [{src} + {8 * k}]")
                    e(f"mov qword ptr [rsp + {off + 8 * k}], rax")

    def row0(z):
        e(f"mov rdx, {plan.b(0)}")
        e(f"mulx {z[1]}, {z[0]}, {plan.a(0)}")
        e(f"xor {d32(z[nw])}, {d32(z[nw])}")
        for k in range(1, nw):
            e(f"mulx {z[k + 1]}, {T1}, {plan.a(k)}")
            e(f"adcx {z[k]}, {T1}")
        e(f"adc {z[nw]}, 0")

    def accumulate(z, sel, clear):
        """z += rdx * operand(sel), two carry chains; `clear` is zeroed."""
        e(f"xor {d32(clear)}, {d32(clear)}")
        e(f"mulx {T0}, {T1}, {sel(0)}")
        e(f"adox {z[0]}, {T1}")
        e(f"adox {z[1]}, {T0}")
        for k in range(1, nw):
            e(f"mulx {T0}, {T1}, {sel(k)}")
            e(f"adcx {z[k]}, {T1}")
            e(f"adox {z[k + 1]}, {T0}")
        e(f"adc {z[nw]}, 0")

    def fold(z):
        """m = z0 * ninv; z += m n (z0 becomes 0)."""
        e(f"mov rdx, {z[0]}")
        e(f"imul rdx, {plan.ninv()}")
        accumulate(z, plan.n, T0)

    z = Z[:]
    row0(z)
    fold(z)
    for i in range(1, nw):
        l = z[1:] + [z[0]]
        e(f"mov rdx, {plan.b(i)}")
        accumulate(l, plan.a, l[nw])
        fold(l)
        z = l
    res = z[1:]
    if plan.park_c:
        e(f"mov rdx, qword ptr [rsp + {plan.off_c}]")
        for k, r in enumerate(res):
            e(f"mov qword ptr [rdx + {8 * k}], {r}")
    else:
        for k, r in enumerate(res):
            e(f"mov qword ptr [rdi + {8 * k}], {r}")
    if plan.frame:
        e(f"add rsp, {plan.frame}")
    for r in reversed(plan.saved):
        e(f"pop {r}")
    return "\n".join(out)


def mont_operands():
    lines = [
        '        inout("rsi") a => _,',
        '        inout("rcx") b => _,',
        '        inout("rdi") c => _,',
        '        inout("rdx") n => _,',
        '        out("rax") _,',
    ]
    for r in GP:
        lines.append(f'        out("{r}") _,')
    lines.append("        options(),")
    return "\n".join(lines)


def generate_mont():
    kernels = []
    arms = []
    for nw in MONT_WORDS:
        plan = MontPlan(nw)
        kernels.append(
            f"""/// `{nw}` words.
///
/// # Safety
/// As [`mont_mul_var`]'s contract, with `nw == {nw}`.
#[inline]
unsafe fn mont_mul_var{nw}_asm(c: *mut u64, a: *const u64, b: *const u64, n: *const u64) {{
    core::arch::asm!(
{emit_mont_asm(plan)}
{mont_operands()}
    );
}}
"""
        )
        arms.append(f"        {nw} => mont_mul_var{nw}_asm(c.as_mut_ptr(), a.as_ptr(), b.as_ptr(), n.as_ptr()),")
    lo, hi = MONT_WORDS[0], MONT_WORDS[-1]
    return f'''//! Montgomery multiplication modulo a run-time odd modulus of {lo} to {hi}
//! 64-bit words in x86-64 inline assembly (`mulx`, `adcx`, `adox`), for the
//! Miller-Rabin rounds of the hash's membership test (`hash::Mont`, P26
//! item 7): CIOS, one row of `mulx` per word of `b` accumulated with the two
//! carry chains, then the fold `m = z_0 (-n^-1)`, `z += m n`, a second row.
//!
//! GENERATED by `tools/gen-fp-asm/gen_fp_asm.py`. Do not edit; change the
//! generator.
//!
//! The kernels keep the accumulator in `NW + 1` registers and never carry
//! out of the top one; that is sound when `n < 2^(64 NW - 2)` and `a, b <
//! n` (then `z < 2n < 2^(64 NW - 1)`, and each row adds less than
//! `2^(64 NW + 62)`), which the caller guarantees. Registers: `r8`-`r15`
//! first, then `rdi`, `rsi`, `rcx`, `rbx` as the word count grows, each
//! freed by parking the pointer it held or copying the operand it
//! addressed onto the stack; `rbx` and `rbp` are pushed inside the block
//! since Rust reserves them. Run only when `fp::dispatch::asm_available()`.

/// `c <- a b / 2^(64 nw) mod n`, reduced to below `2n` (the caller
/// subtracts `n` once conditionally). `n` is `nw` words followed by
/// `-n^-1 mod 2^64`. Returns `false`, computing nothing, when `nw` has no
/// kernel.
///
/// Contract: `bitlen(n) <= 64 nw - 2`, `a, b < n` (their words above `nw`
/// zero), and the CPU has BMI2 and ADX.
#[inline]
pub(crate) fn mont_mul_var(
    nw: usize,
    c: &mut [u64; crate::hash::MAX_DEGREE_WORDS],
    a: &[u64; crate::hash::MAX_DEGREE_WORDS],
    b: &[u64; crate::hash::MAX_DEGREE_WORDS],
    n: &[u64; crate::hash::MAX_DEGREE_WORDS + 1],
) -> bool {{
    // SAFETY: live arrays of at least `nw + 1` words each (`nw <= {hi} <=
    // MAX_DEGREE_WORDS`), the output aliasing no input; the arithmetic
    // contract is the caller's (`hash::Mont::new` checks it) and BMI2/ADX
    // are detected before this is chosen.
    unsafe {{
        match nw {{
{chr(10).join(arms)}
            _ => return false,
        }}
    }}
    true
}}

{"".join(kernels)}'''


def rustfmt(text):
    """Format as `cargo fmt` would, so the checked-in files are fmt-stable."""
    import subprocess
    r = subprocess.run(["rustfmt", "--edition", "2021", "--emit", "stdout"], input=text, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("rustfmt failed:\n" + r.stderr)
    return r.stdout


def main():
    args = sys.argv[1:]
    check = "--check" in args
    without_mont = "--without-mont" in args
    global CRATE
    if "--crate-path" in args:
        i = args.index("--crate-path")
        CRATE = args[i + 1]
        del args[i : i + 2]
    args = [a for a in args if a not in ("--check", "--without-mont")]
    d = args[0] if args else "crates/prism-verify/src/fp"
    bad = 0
    outputs = [(os.path.join(d, f"{name}_asm.rs"), lambda t=(name, c, e, n, enc, origin): generate(*t)) for name, c, e, n, enc, origin in PRIMES]
    if not without_mont:
        outputs.append((os.path.normpath(os.path.join(d, "..", "mont_asm.rs")), generate_mont))
    for path, gen in outputs:
        text = rustfmt(gen())
        if check:
            if not os.path.exists(path) or open(path).read() != text:
                print(f"{path}: out of date", file=sys.stderr)
                bad += 1
        else:
            open(path, "w").write(text)
    if check:
        print(f"all {len(outputs)} generated x86-64 modules up to date" if not bad else f"{bad} out of date")
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
