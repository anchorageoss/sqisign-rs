//! Computes the constants for one prime or one level and renders them as
//! Rust source.

use crate::e0::{actions, check_actions, e0_basis, x_of_difference, Actions, E0Basis, Mat};
use crate::field::prime;
use crate::numtheory::{ei_box, isqrt, next_prime, ri_cofactor};
use hybrid_array::typenum::Unsigned;
use num_bigint::BigUint;
use num_traits::One;
use sqisign_verify::fp::FpBackend;
use sqisign_verify::params::{Prime, P324_3, P500_27, P664_17};
use std::fmt::Write as _;

/// Pipe generated source through `rustfmt`.
pub fn rustfmt(src: String) -> String {
    use std::io::Write as _;
    use std::process::{Command, Stdio};
    let mut child = Command::new("rustfmt")
        .args(["--edition", "2021", "--emit", "stdout"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("rustfmt must be installed (rustup component add rustfmt)");
    child
        .stdin
        .take()
        .expect("stdin")
        .write_all(src.as_bytes())
        .expect("feed rustfmt");
    let out = child.wait_with_output().expect("rustfmt exit");
    assert!(out.status.success(), "rustfmt rejected generated code");
    String::from_utf8(out.stdout).expect("utf8")
}

fn le_bytes(x: &BigUint, len: usize) -> Vec<u8> {
    let mut v = x.to_bytes_le();
    assert!(v.len() <= len, "value too large for {len} bytes");
    v.resize(len, 0);
    v
}

fn byte_array(out: &mut String, doc: &str, name: &str, bytes: &[u8]) {
    writeln!(out, "/// {doc}").unwrap();
    writeln!(out, "pub const {name}: [u8; {}] = [", bytes.len()).unwrap();
    for chunk in bytes.chunks(16) {
        let row: Vec<String> = chunk.iter().map(|b| format!("0x{b:02x}")).collect();
        writeln!(out, "    {},", row.join(", ")).unwrap();
    }
    writeln!(out, "];\n").unwrap();
}

fn matrix_const(out: &mut String, doc: &str, name: &str, m: &Mat, len: usize) {
    writeln!(out, "/// {doc}").unwrap();
    writeln!(out, "pub const {name}: Mat2x2<ENTRY_BYTES> = [").unwrap();
    for row in m {
        writeln!(out, "    [").unwrap();
        for entry in row {
            writeln!(out, "        [").unwrap();
            for chunk in le_bytes(entry, len).chunks(16) {
                let r: Vec<String> = chunk.iter().map(|b| format!("0x{b:02x}")).collect();
                writeln!(out, "            {},", r.join(", ")).unwrap();
            }
            writeln!(out, "        ],").unwrap();
        }
        writeln!(out, "    ],").unwrap();
    }
    writeln!(out, "];\n").unwrap();
}

/// Everything derived from one prime.
pub struct PrimeData<L: FpBackend> {
    pub name: &'static str,
    pub basis: E0Basis<L>,
    pub actions: Actions,
}

impl<L: FpBackend> PrimeData<L> {
    pub fn compute(name: &'static str) -> Self {
        eprintln!("{name}: E0 basis ...");
        let basis = e0_basis::<L>();
        eprintln!("{name}: endomorphism actions ...");
        let acts = actions::<L>(&basis);
        check_actions::<L>(basis.f, &acts);
        Self {
            name,
            basis,
            actions: acts,
        }
    }

    fn f(&self) -> u64 {
        self.basis.f
    }

    fn torsion_bytes(&self) -> usize {
        (self.f() as usize).div_ceil(8)
    }

    fn log2p(&self) -> u64 {
        prime::<L>().bits()
    }

    /// The largest `lambda` this prime is used at (round 3, and round 2 for
    /// the primes kept for the PRISM_v2 comparison).
    fn max_lambda(&self) -> u64 {
        match self.name {
            "p324_3" | "p248_5" => 128,
            "p376_65" => 192,
            "p500_27" | "p664_17" => 256,
            _ => 256,
        }
    }

    /// Limbs of the fixed-precision integers: the reference's `5 * ceil(log2
    /// p / 64)`, or more if a reduced ideal needs it. The constant-time
    /// reduction of an ideal of norm bound `nb` needs `2 (nb + 3) + log2 p +
    /// 3` bits; salt-PRISM reduces the intersection of the secret ideal
    /// (norm bound `EQUIV_NORM_BITS + 1`) with the challenge ideal of norm
    /// `q (2^a - q) < 2^(2a)`, `a = lambda + 64`, and key generation reduces
    /// the ideal of norm `D_mix` (`log2 p + 2 lambda` bits). Only the
    /// round-2 primes, whose `lambda` is large relative to `log2 p`, exceed
    /// the reference's count (`p248_5`: 21, `p376_65`: 31).
    fn ibz_nlimbs(&self) -> u64 {
        let reference = 5 * self.log2p().div_ceil(64);
        let lambda = self.max_lambda();
        let bits_for = |nb: u64| (2 * (nb + 3) + self.log2p() + 3).div_ceil(64);
        // salt-PRISM's intersection I_sk ∩ I_chall
        let a = lambda + 64;
        let equiv_norm_bits = self.log2p() / 2 + 15;
        let intersection = bits_for(equiv_norm_bits + 1 + 2 * a + 2);
        // key generation's reduction of the ideal of norm D_mix = next_prime(p 2^(2 lambda))
        let keygen = bits_for(self.log2p() + 2 * lambda + 3);
        reference.max(intersection).max(keygen)
    }

    fn header(&self, what: &str) -> String {
        format!(
            "//! {what} for `p = {c} * 2^{e} - 1` ({bits} bits).\n\
             //!\n\
             //! GENERATED by `tools/gen-precomp`. Do not edit; change the generator.\n\
             //! Derivations are documented in PRECOMP.md.\n",
            c = L::COFACTOR,
            e = self.f(),
            bits = self.log2p(),
        )
    }

    /// Constants for `sqisign-verify`: the `E0[2^f]` basis and small sizes.
    pub fn verify_source(&self) -> String {
        let mut o = self.header("Precomputed curve constants");
        writeln!(o).unwrap();
        writeln!(
            o,
            "/// Two-adic exponent `f`: the rational torsion is `E0[2^f]`."
        )
        .unwrap();
        writeln!(o, "pub const TWO_ADIC_EXPONENT: u32 = {};\n", self.f()).unwrap();
        writeln!(
            o,
            "/// Byte length of `2^f` itself, `ceil((f + 1) / 8)`, as the reference's\n\
             /// `TORSION_2POWER_BYTES` computes it. The spec's Table 10 prints `ceil(f / 8)`,\n\
             /// which differs at level V (83 vs 84). Integers reduced modulo `2^f` need one\n\
             /// byte less when `8 | f`; see `ENTRY_BYTES` in the signing crate."
        )
        .unwrap();
        writeln!(
            o,
            "pub const TORSION_2POWER_BYTES: usize = {};\n",
            (self.f() as usize + 1).div_ceil(8)
        )
        .unwrap();
        writeln!(o, "/// Bit length of the odd cofactor `(p + 1) / 2^f`.").unwrap();
        writeln!(
            o,
            "pub const COFACTOR_BITLENGTH: usize = {};\n",
            64 - L::COFACTOR.leading_zeros()
        )
        .unwrap();
        writeln!(
            o,
            "/// 64-bit words needed to hold an integer of `log2 p` bits (`NWORDS_ORDER`)."
        )
        .unwrap();
        writeln!(
            o,
            "pub const NWORDS_ORDER: usize = {};\n",
            self.log2p().div_ceil(64)
        )
        .unwrap();
        writeln!(o, "/// `ceil(log2 log2 p)`.").unwrap();
        writeln!(
            o,
            "pub const LOG2P: u32 = {};\n",
            64 - (self.log2p() - 1).leading_zeros()
        )
        .unwrap();
        let n2 = 2 * L::FpEncodedBytes::USIZE;
        let px = self.basis.p.as_ref().expect("finite").0.encode();
        let qx = self.basis.q.as_ref().expect("finite").0.encode();
        let pmqx = x_of_difference(&self.basis).encode();
        byte_array(
            &mut o,
            &format!(
                "Canonical `Fp2` encoding ({n2} bytes) of `x(P0)`, the first generator of the\n\
                 /// deterministic basis of `E0[2^f]` (spec Algorithm B.1)."
            ),
            "E0_BASIS_PX",
            &px,
        );
        byte_array(
            &mut o,
            "Canonical `Fp2` encoding of `x(Q0)`, the second generator; `[2^(f-1)] Q0 = (0, 0)`.",
            "E0_BASIS_QX",
            &qx,
        );
        byte_array(
            &mut o,
            "Canonical `Fp2` encoding of `x(P0 - Q0)`.",
            "E0_BASIS_PMQX",
            &pmqx,
        );
        o
    }

    /// Constants for `sqisign-rs` (signing side): action matrices and quaternion sizes.
    pub fn sign_source(&self) -> String {
        let mut o = self.header("Precomputed endomorphism-ring constants");
        let tb = self.torsion_bytes();
        let p = prime::<L>();
        let sqrt_p = isqrt(&p);
        writeln!(o).unwrap();
        writeln!(o, "use super::Mat2x2;\n").unwrap();
        writeln!(
            o,
            "/// Each matrix entry is an integer modulo `2^f`, little-endian, this many bytes."
        )
        .unwrap();
        writeln!(o, "pub const ENTRY_BYTES: usize = {tb};\n").unwrap();
        let doc = |what: &str| {
            format!(
                "Action of `{what}` on the basis `(P0, Q0)` of `E0[2^f]`: entry `[0][0]` is `a(P0)`,\n\
                 /// `[0][1]` is `a(Q0)`, `[1][0]` is `b(P0)`, `[1][1]` is `b(Q0)`, where\n\
                 /// `alpha(R) = a(R) P0 + b(R) Q0`."
            )
        };
        matrix_const(
            &mut o,
            &doc("i: (x, y) -> (-x, i y)"),
            "ACTION_I",
            &self.actions.i,
            tb,
        );
        matrix_const(
            &mut o,
            &doc("j: (x, y) -> (x^p, y^p)"),
            "ACTION_J",
            &self.actions.j,
            tb,
        );
        matrix_const(&mut o, &doc("k = i j"), "ACTION_K", &self.actions.k, tb);
        matrix_const(
            &mut o,
            &doc("the O0 generator i"),
            "ACTION_GEN2",
            &self.actions.gen2,
            tb,
        );
        matrix_const(
            &mut o,
            &doc("the O0 generator (i + j)/2"),
            "ACTION_GEN3",
            &self.actions.gen3,
            tb,
        );
        matrix_const(
            &mut o,
            &doc("the O0 generator (1 + k)/2"),
            "ACTION_GEN4",
            &self.actions.gen4,
            tb,
        );
        writeln!(o, "/// `log2 p`.").unwrap();
        writeln!(o, "pub const P_BITS: u32 = {};\n", self.log2p()).unwrap();
        writeln!(
            o,
            "/// Limbs of the fixed-precision integers of this parameter set: `5 * ceil(log2 p / 64)`."
        )
        .unwrap();
        writeln!(o, "pub const IBZ_NLIMBS: usize = {};\n", self.ibz_nlimbs()).unwrap();
        writeln!(o, "/// Bit length of `isqrt(p)`.").unwrap();
        writeln!(o, "pub const SQRT_P_BITS: u32 = {};\n", sqrt_p.bits()).unwrap();
        byte_array(
            &mut o,
            "`isqrt(p)`, little-endian.",
            "SQRT_P_LE",
            &le_bytes(&sqrt_p, L::FpEncodedBytes::USIZE),
        );
        writeln!(
            o,
            "/// Power of two the ideal-to-isogeny translation may use: `f - 2`."
        )
        .unwrap();
        writeln!(
            o,
            "pub const QLAPOTY_USED_POWER_OF_TWO: u32 = {};\n",
            self.f() - 2
        )
        .unwrap();
        writeln!(o, "/// Bit size of ideal norms after equivalent-ideal reduction: `floor(log2 p / 2) + 15`.").unwrap();
        writeln!(
            o,
            "pub const EQUIV_NORM_BITS: u32 = {};\n",
            self.log2p() / 2 + 15
        )
        .unwrap();
        writeln!(
            o,
            "/// Bit size of the two-ideal product norm: `2 * EQUIV_NORM_BITS`."
        )
        .unwrap();
        writeln!(
            o,
            "pub const PROD_NORM_BITS: u32 = {};",
            2 * (self.log2p() / 2 + 15)
        )
        .unwrap();
        o
    }
}

/// Human-readable dump of one prime's values, for debugging.
pub fn prime_dump<L: FpBackend>() -> String {
    let data = PrimeData::<L>::compute("dump");
    let mut o = String::new();
    let hex = |m: &Mat| {
        format!(
            "[[{:x}, {:x}], [{:x}, {:x}]]",
            m[0][0], m[0][1], m[1][0], m[1][1]
        )
    };
    writeln!(o, "f = {}", data.f()).unwrap();
    writeln!(
        o,
        "x(P0) = {}",
        crate::field::fp2_hex(&data.basis.p.as_ref().unwrap().0)
    )
    .unwrap();
    writeln!(
        o,
        "x(Q0) = {}",
        crate::field::fp2_hex(&data.basis.q.as_ref().unwrap().0)
    )
    .unwrap();
    writeln!(o, "i    = {}", hex(&data.actions.i)).unwrap();
    writeln!(o, "j    = {}", hex(&data.actions.j)).unwrap();
    writeln!(o, "k    = {}", hex(&data.actions.k)).unwrap();
    writeln!(o, "gen3 = {}", hex(&data.actions.gen3)).unwrap();
    writeln!(o, "gen4 = {}", hex(&data.actions.gen4)).unwrap();
    o
}

/// A SQIsign round-3 parameter set.
struct Level {
    module: &'static str,
    lambda: u64,
    prime: BigUint,
    fp_bytes: usize,
    f: u64,
}

fn level<L: Prime>(module: &'static str, lambda: u64) -> Level {
    Level {
        module,
        lambda,
        prime: prime::<L>(),
        fp_bytes: L::FpEncodedBytes::USIZE,
        f: L::TWO_ADIC_EXPONENT as u64,
    }
}

/// The round-3 level constants for all three parameter sets, as one module.
pub fn sqisign_v3_levels() -> String {
    let levels = [
        level::<P324_3>("p324_3", 128),
        level::<P500_27>("p500_27", 192),
        level::<P664_17>("p664_17", 256),
    ];
    levels_module(
        &levels,
        "//! SQIsign round-3 level constants (spec Tables 8, 9, 11 and Appendix B),\n\
         //! derived from `(p, lambda)`. Needed by the substrate check that runs\n\
         //! SQIsign v3 on the ported layers; PRISM defines its own level constants.\n\
         //!\n\
         //! GENERATED by `tools/gen-precomp`. Do not edit; change the generator.\n\
         //! Derivations are documented in PRECOMP.md.\n\n",
    )
}

fn levels_module(levels: &[Level], header: &str) -> String {
    let mut o = String::from(header);
    for lv in levels {
        let log2p = lv.prime.bits();
        // e_rsp = 1 + ceil(log2 p / 2 + lambda / 4)
        let e_rsp = 1 + (2 * log2p + lv.lambda).div_ceil(4);
        let response_bytes = (e_rsp + 2).div_ceil(8);
        let challenge_bytes = lv.lambda / 8;
        let fp = lv.fp_bytes as u64;
        let pk_bytes = 2 * fp + 1;
        let sk_bytes = pk_bytes + 3 * fp + 4 * challenge_bytes;
        let sig_bytes = 2 * fp + 4 * response_bytes + challenge_bytes + 2;
        let deg_mix = next_prime(&(&lv.prime << (2 * lv.lambda)));
        let ri = ri_cofactor(&lv.prime, log2p, e_rsp, lv.lambda);
        let eibox = ei_box(log2p, lv.lambda);
        let primality_iter = lv.lambda.div_ceil(2);
        let deg_mix_bytes = le_bytes(&deg_mix, (deg_mix.bits() as usize).div_ceil(8));
        let ri_bytes = le_bytes(&ri, (ri.bits() as usize).div_ceil(8));
        writeln!(
            o,
            "/// NIST level with `lambda = {}` on `{}`.",
            lv.lambda, lv.module
        )
        .unwrap();
        writeln!(o, "pub mod {} {{", lv.module).unwrap();
        writeln!(
            o,
            "    /// Security parameter `lambda`; also the challenge length in bits."
        )
        .unwrap();
        writeln!(o, "    pub const LAMBDA: u32 = {};", lv.lambda).unwrap();
        writeln!(
            o,
            "    /// Response exponent `e_rsp = 1 + ceil((lambda + 2 log2 p) / 4)`."
        )
        .unwrap();
        writeln!(o, "    pub const RESPONSE_BITS: u32 = {e_rsp};").unwrap();
        writeln!(o, "    /// `ceil((e_rsp + 2) / 8)`.").unwrap();
        writeln!(o, "    pub const RESPONSE_BYTES: usize = {response_bytes};").unwrap();
        writeln!(o, "    /// `lambda / 8`.").unwrap();
        writeln!(
            o,
            "    pub const CHALLENGE_BYTES: usize = {challenge_bytes};"
        )
        .unwrap();
        writeln!(o, "    /// Hash-to-challenge iterations.").unwrap();
        writeln!(o, "    pub const HASH_ITERATIONS: u32 = 1;").unwrap();
        writeln!(o, "    /// Two-adic exponent of the prime.").unwrap();
        writeln!(o, "    pub const TORSION_EVEN_POWER: u32 = {};", lv.f).unwrap();
        writeln!(o, "    /// `2 FP_ENCODED_BYTES + 1`.").unwrap();
        writeln!(o, "    pub const PUBLICKEY_BYTES: usize = {pk_bytes};").unwrap();
        writeln!(
            o,
            "    /// `PUBLICKEY_BYTES + 3 FP_ENCODED_BYTES + 4 CHALLENGE_BYTES`."
        )
        .unwrap();
        writeln!(o, "    pub const SECRETKEY_BYTES: usize = {sk_bytes};").unwrap();
        writeln!(
            o,
            "    /// `2 FP_ENCODED_BYTES + 4 RESPONSE_BYTES + CHALLENGE_BYTES + 2`."
        )
        .unwrap();
        writeln!(o, "    pub const SIGNATURE_BYTES: usize = {sig_bytes};").unwrap();
        writeln!(
            o,
            "    /// Miller-Rabin rounds for the quaternion layer: `ceil(lambda / 2)`."
        )
        .unwrap();
        writeln!(
            o,
            "    pub const PRIMALITY_NUM_ITER: u32 = {primality_iter};"
        )
        .unwrap();
        writeln!(o, "    /// EIBox (spec Algorithm B.2).").unwrap();
        writeln!(o, "    pub const EQUIV_BOUND_COEFF: u32 = {eibox};").unwrap();
        writeln!(o, "    /// Bit length of `D_mix`.").unwrap();
        writeln!(
            o,
            "    pub const DEGREE_NORM_BITS: u32 = {};",
            deg_mix.bits()
        )
        .unwrap();
        writeln!(
            o,
            "    /// `D_mix = next_prime(p * 2^(2 lambda))`, the secret-key and commitment ideal\n\
             \x20   /// norm, little-endian."
        )
        .unwrap();
        writeln!(
            o,
            "    pub const DEG_MIX_LE: [u8; {}] = [",
            deg_mix_bytes.len()
        )
        .unwrap();
        for chunk in deg_mix_bytes.chunks(16) {
            let r: Vec<String> = chunk.iter().map(|b| format!("0x{b:02x}")).collect();
            writeln!(o, "        {},", r.join(", ")).unwrap();
        }
        writeln!(o, "    ];").unwrap();
        writeln!(
            o,
            "    /// RICofactor (spec Algorithm B.3): the prime cofactor used by\n\
             \x20   /// RandomIdealGivenNorm, little-endian."
        )
        .unwrap();
        writeln!(
            o,
            "    pub const RI_COFACTOR_LE: [u8; {}] = [",
            ri_bytes.len()
        )
        .unwrap();
        for chunk in ri_bytes.chunks(16) {
            let r: Vec<String> = chunk.iter().map(|b| format!("0x{b:02x}")).collect();
            writeln!(o, "        {},", r.join(", ")).unwrap();
        }
        writeln!(o, "    ];").unwrap();
        writeln!(o, "}}\n").unwrap();
    }
    let _ = BigUint::one();
    o
}
