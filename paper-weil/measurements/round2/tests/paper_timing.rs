//! Verification cost breakdown for the compression paper (Task 2).
//!
//! Instant-based timing on KAT entry 0 of each level. Level 1 gets the full
//! per-component table; Levels 3 and 5 get full standard-versus-compressed
//! verification plus the decompress and protocols_verify checkpoints.
//!
//! Run with:
//!   CARGO_PROFILE_RELEASE_LTO=fat CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 \
//!   cargo test -p sqisign-kat --release --features kat-compat \
//!     --test paper_timing -- --nocapture --test-threads=1

#![cfg(feature = "kat-compat")]

use std::fmt::Write as _;
use std::hint::black_box;
use std::time::Instant;

use sqisign_kat::kat_parser;
use sqisign_verify::ec::pairing::{fp2_dlog_2e_pub, weil};
use sqisign_verify::ec::point::{ec_dbl_iter_basis, xadd};
use sqisign_verify::ec::EcCurve;
use sqisign_verify::formats::{basis_to_hint, hensel_inv_mod_2e, mp_mod_2exp, mp_mul_mod, set_high_bits};
use sqisign_verify::hash::hash_to_challenge;
use sqisign_verify::theta::HD_EXTRA_TORSION;
use sqisign_verify::verify::{
    basis_from_hint, check_canonical_basis_change_matrix, compute_challenge_curve,
    compute_commitment_curve_verify, matrix_scalar_application_even_basis, mp_compare, mp_is_even,
    protocols_verify, two_response_isogeny_verify_inner, verify_canonical_hint,
};
use sqisign_verify::{
    CompressedSignature, FpBackend, Level1, Level3, Level5, LevelPrecomp, PublicKey, Scalar,
    Signature, Verifier,
};

const KAT_L1: &str = include_str!("../../../reference/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp");
const KAT_L3: &str = include_str!("../../../reference/KAT/PQCsignKAT_529_SQIsign_lvl3.rsp");
const KAT_L5: &str = include_str!("../../../reference/KAT/PQCsignKAT_701_SQIsign_lvl5.rsp");

const ITERS: usize = 200;

#[derive(Clone, Copy)]
struct Stat {
    median_us: f64,
    min_us: f64,
}

/// Time `f` over `ITERS` iterations. `setup` runs outside the timed region
/// and produces a fresh input for each iteration (curves are mutated by
/// several of the routines, so each iteration starts from the captured state).
fn measure<I, T>(mut setup: impl FnMut() -> I, mut f: impl FnMut(I) -> T) -> Stat {
    for _ in 0..5 {
        black_box(f(setup()));
    }
    let mut samples = Vec::with_capacity(ITERS);
    for _ in 0..ITERS {
        let input = setup();
        let t0 = Instant::now();
        let out = f(input);
        let dt = t0.elapsed();
        black_box(out);
        samples.push(dt.as_secs_f64() * 1e6);
    }
    samples.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_us = if ITERS % 2 == 0 {
        (samples[ITERS / 2 - 1] + samples[ITERS / 2]) / 2.0
    } else {
        samples[ITERS / 2]
    };
    Stat {
        median_us,
        min_us: samples[0],
    }
}

/// The dropped-entry recovery exactly as `CompressedSignature::decompress`
/// performs it after the discrete logarithm: Hensel inverse of the pivot,
/// one product, one add or sub, one product, then the two hint bits.
fn recover_dropped_entry<L: FpBackend>(
    mat: &[[Scalar<L>; 2]; 2],
    det: &[u64],
    det_precision: usize,
    det_hint: u8,
) -> Scalar<L> {
    let nw = det.len();
    let m00_odd = mat[0][0].digits()[0] & 1 != 0;
    let (pivot, other, mat_var) = if m00_odd {
        (&mat[0][0], &mat[0][1], &mat[1][0])
    } else {
        (&mat[0][1], &mat[0][0], &mat[1][1])
    };
    let mut inv_pivot = Scalar::<L>::default();
    hensel_inv_mod_2e(
        Signature::<L>::scalar_digits_mut(&mut inv_pivot),
        pivot.digits(),
        det_precision,
    );
    let mut product = Scalar::<L>::default();
    mp_mul_mod(
        Signature::<L>::scalar_digits_mut(&mut product),
        other.digits(),
        mat_var.digits(),
        det_precision,
    );
    let mut numerator = Scalar::<L>::default();
    {
        let nd = Signature::<L>::scalar_digits_mut(&mut numerator);
        let pd = product.digits();
        if m00_odd {
            let mut carry: u64 = 0;
            for i in 0..nw {
                let sum = (det[i] as u128) + (pd[i] as u128) + (carry as u128);
                nd[i] = sum as u64;
                carry = (sum >> 64) as u64;
            }
        } else {
            let mut borrow: u64 = 0;
            for i in 0..nw {
                let (diff, b1) = pd[i].overflowing_sub(det[i]);
                let (diff2, b2) = diff.overflowing_sub(borrow);
                nd[i] = diff2;
                borrow = (b1 as u64) + (b2 as u64);
            }
        }
        mp_mod_2exp(nd, det_precision);
    }
    let mut recovered = Scalar::<L>::default();
    mp_mul_mod(
        Signature::<L>::scalar_digits_mut(&mut recovered),
        numerator.digits(),
        inv_pivot.digits(),
        det_precision,
    );
    set_high_bits(
        Signature::<L>::scalar_digits_mut(&mut recovered),
        det_hint,
        det_precision,
        HD_EXTRA_TORSION as usize,
    );
    recovered
}

struct Row {
    name: &'static str,
    std_count: usize,
    comp_count: usize,
    stat: Stat,
}

fn row(name: &'static str, std_count: usize, comp_count: usize, stat: Stat) -> Row {
    Row {
        name,
        std_count,
        comp_count,
        stat,
    }
}

fn breakdown<L: FpBackend + LevelPrecomp>(name: &str, rsp: &str) -> String {
    let entries = kat_parser::parse_rsp(rsp);
    let e = &entries[0];
    assert_eq!(e.count, 0);
    let pk = PublicKey::<L>::from_bytes(&e.pk).unwrap();
    let sig_len = e.sm.len() - e.msg.len();
    let sig_bytes = &e.sm[..sig_len];
    let msg: &[u8] = &e.msg;
    let sig = Signature::<L>::from_bytes(sig_bytes).unwrap();
    let comp = sig.compress();
    let cbytes = comp.to_bytes();

    let bt = sig.backtracking() as usize;
    let trl = sig.two_resp_length() as usize;
    let hd = HD_EXTRA_TORSION as usize;
    let pow_dim2 = L::E_RSP as usize - trl - bt;
    let pow_i32 = pow_dim2 as i32;
    let f = pow_dim2 + hd + trl;
    let g = pow_dim2 + hd;
    let det_precision = pow_dim2 + trl;
    let f_chr = L::F_CHR as usize;
    let dbl_chall = f_chr - f;
    let dbl_aux = f_chr - g;
    let nw = sig.chall_coeff().digits().len();
    let fp2_len = sig.e_aux_a().encode().len();
    let hint_aux = sig.hint_aux();
    let hint_chall = sig.hint_chall();

    let mut out = String::new();
    let w = &mut out;
    writeln!(w, "==== {name}: KAT entry 0, bt={bt}, trl={trl}, pow_dim2={pow_dim2}, N=f={f}, g={g}, det_precision={det_precision}, F_CHR={f_chr}, msg {} bytes, sig {} bytes, compressed {} bytes, {ITERS} iterations ====", msg.len(), sig_len, cbytes.len()).unwrap();

    // ---- Full measurements from bytes (public key decoded once, outside the loop).
    let full_std = measure(|| (), |_| {
        let s = Signature::<L>::from_bytes(sig_bytes).unwrap();
        pk.verify(msg, &s).unwrap()
    });
    let full_comp = measure(|| (), |_| {
        let c = CompressedSignature::<L>::from_bytes(cbytes.as_slice()).unwrap();
        pk.verify(msg, &c).unwrap()
    });
    let pk_decode = measure(|| (), |_| PublicKey::<L>::from_bytes(&e.pk).unwrap());
    let pv_only = measure(|| sig.clone(), |s| protocols_verify(&pk, msg, &s).unwrap());
    let decompress_only = measure(|| comp.clone(), |c| c.decompress(&pk).unwrap());

    // ---- Standard-path components, captured state by state.
    let c_decode_std = measure(|| (), |_| Signature::<L>::from_bytes(sig_bytes).unwrap());
    let c_decode_comp = measure(|| (), |_| CompressedSignature::<L>::from_bytes(cbytes.as_slice()).unwrap());
    let c_precheck = measure(|| (), |_| {
        check_canonical_basis_change_matrix(&sig).unwrap();
        EcCurve::<L>::verify_a(&pk.curve().a)
    });

    let e_aux0 = EcCurve::<L>::from_a(sig.e_aux_a()).unwrap();
    let c_from_a = measure(|| (), |_| EcCurve::<L>::from_a(sig.e_aux_a()).unwrap());

    let mut e_aux1 = e_aux0.clone();
    assert!(verify_canonical_hint::<L>(&mut e_aux1, hint_aux));
    let c_tohint_aux = measure(|| e_aux0.clone(), |mut c| verify_canonical_hint::<L>(&mut c, hint_aux));

    let e_chall0 = compute_challenge_curve(sig.chall_coeff(), sig.backtracking(), pk.curve(), pk.hint_pk()).unwrap();
    let c_chall = measure(|| (), |_| {
        compute_challenge_curve(sig.chall_coeff(), sig.backtracking(), pk.curve(), pk.hint_pk()).unwrap()
    });

    let mut e_chall1 = e_chall0.clone();
    assert!(verify_canonical_hint::<L>(&mut e_chall1, hint_chall));
    let c_tohint_chall = measure(|| e_chall0.clone(), |mut c| verify_canonical_hint::<L>(&mut c, hint_chall));

    let mut e_chall2 = e_chall1.clone();
    let b_chall_raw = basis_from_hint(&mut e_chall2, L::F_CHR, hint_chall).unwrap();
    let c_fromhint_chall = measure(|| e_chall1.clone(), |mut c| basis_from_hint(&mut c, L::F_CHR, hint_chall).unwrap());

    let mut e_chall3 = e_chall2.clone();
    let b_chall_dbl = ec_dbl_iter_basis(&b_chall_raw, dbl_chall, &mut e_chall3);
    let c_dbl_chall = measure(|| e_chall2.clone(), |mut c| ec_dbl_iter_basis(&b_chall_raw, dbl_chall, &mut c));

    let mut e_aux2 = e_aux1.clone();
    let b_aux_raw = basis_from_hint(&mut e_aux2, L::F_CHR, hint_aux).unwrap();
    let c_fromhint_aux = measure(|| e_aux1.clone(), |mut c| basis_from_hint(&mut c, L::F_CHR, hint_aux).unwrap());

    let mut e_aux3 = e_aux2.clone();
    let b_aux_dbl = ec_dbl_iter_basis(&b_aux_raw, dbl_aux, &mut e_aux3);
    let c_dbl_aux = measure(|| e_aux2.clone(), |mut c| ec_dbl_iter_basis(&b_aux_raw, dbl_aux, &mut c));

    let mut b_chall_mat = b_chall_dbl.clone();
    matrix_scalar_application_even_basis(&mut b_chall_mat, &e_chall3, sig.mat(), f).unwrap();
    let c_matrix = measure(|| b_chall_dbl.clone(), |mut b| {
        matrix_scalar_application_even_basis(&mut b, &e_chall3, sig.mat(), f).unwrap();
        b
    });

    let kernel_is_q = mp_is_even::<L>(&sig.mat()[0][0]) && mp_is_even::<L>(&sig.mat()[1][0]);
    let (e_chall4, b_chall_fin, c_small) = if trl > 0 {
        let mut ec = e_chall3.clone();
        let mut b = b_chall_mat.clone();
        two_response_isogeny_verify_inner(&mut ec, &mut b, kernel_is_q, trl as u8, pow_i32).unwrap();
        let st = measure(
            || (e_chall3.clone(), b_chall_mat.clone()),
            |(mut ec, mut b)| {
                two_response_isogeny_verify_inner(&mut ec, &mut b, kernel_is_q, trl as u8, pow_i32).unwrap();
                (ec, b)
            },
        );
        (ec, b, st)
    } else {
        (e_chall3.clone(), b_chall_mat.clone(), Stat { median_us: 0.0, min_us: 0.0 })
    };

    let e_com = compute_commitment_curve_verify(&b_chall_fin, &b_aux_dbl, &e_chall4, &e_aux3, pow_i32).unwrap();
    let c_theta = measure(|| (), |_| {
        compute_commitment_curve_verify(&b_chall_fin, &b_aux_dbl, &e_chall4, &e_aux3, pow_i32).unwrap()
    });

    let chk = hash_to_challenge(&pk, &e_com, msg).unwrap();
    assert_eq!(mp_compare::<L>(sig.chall_coeff(), &chk), 0);
    let c_hash = measure(|| (), |_| {
        let c = hash_to_challenge(&pk, &e_com, msg).unwrap();
        mp_compare::<L>(sig.chall_coeff(), &c)
    });

    // ---- Decompression-only components (CompressedSignature::decompress).
    let mut e_chall_d = e_chall0.clone();
    let (b_chall_th, h_chall) = basis_to_hint::<L>(&mut e_chall_d, L::F_CHR).unwrap();
    assert_eq!(h_chall, hint_chall);
    let c_basis_to_hint_chall = measure(|| e_chall0.clone(), |mut c| basis_to_hint::<L>(&mut c, L::F_CHR).unwrap());

    let mut e_chall_d2 = e_chall_d.clone();
    let b_chall_f = ec_dbl_iter_basis(&b_chall_th, f_chr - f, &mut e_chall_d2);
    let c_dbl_chall_dec = measure(|| e_chall_d.clone(), |mut c| ec_dbl_iter_basis(&b_chall_th, f_chr - f, &mut c));

    let ppq_chall = xadd(&b_chall_f.p, &b_chall_f.q, &b_chall_f.pmq);
    let c_xadd = measure(|| (), |_| xadd(&b_chall_f.p, &b_chall_f.q, &b_chall_f.pmq));

    let mut e_chall_d3 = e_chall_d2.clone();
    let omega_f = weil::<L>(f as u32, &b_chall_f.p, &b_chall_f.q, &ppq_chall, &mut e_chall_d3);
    let c_weil_chall = measure(|| e_chall_d2.clone(), |mut c| {
        weil::<L>(f as u32, &b_chall_f.p, &b_chall_f.q, &ppq_chall, &mut c)
    });

    let mut e_aux_d = e_aux0.clone();
    let (b_aux_th, h_aux) = basis_to_hint::<L>(&mut e_aux_d, L::F_CHR).unwrap();
    assert_eq!(h_aux, hint_aux);
    let c_basis_to_hint_aux = measure(|| e_aux0.clone(), |mut c| basis_to_hint::<L>(&mut c, L::F_CHR).unwrap());

    let mut e_aux_d2 = e_aux_d.clone();
    let b_aux_g = ec_dbl_iter_basis(&b_aux_th, f_chr - g, &mut e_aux_d2);
    let c_dbl_aux_dec = measure(|| e_aux_d.clone(), |mut c| ec_dbl_iter_basis(&b_aux_th, f_chr - g, &mut c));

    let ppq_aux = xadd(&b_aux_g.p, &b_aux_g.q, &b_aux_g.pmq);
    let mut e_aux_d3 = e_aux_d2.clone();
    let omega_aux = weil::<L>(g as u32, &b_aux_g.p, &b_aux_g.q, &ppq_aux, &mut e_aux_d3);
    let c_weil_aux = measure(|| e_aux_d2.clone(), |mut c| {
        weil::<L>(g as u32, &b_aux_g.p, &b_aux_g.q, &ppq_aux, &mut c)
    });

    let omega_f_inv = omega_f.inv();
    let omega_aux_inv = omega_aux.inv();
    let c_inv = measure(|| (), |_| (omega_f.inv(), omega_aux.inv()));

    let mut det = [0u64; 8];
    fp2_dlog_2e_pub::<L>(&mut det[..nw], &omega_aux_inv, &omega_f_inv, f as u32).unwrap();
    mp_mod_2exp(&mut det[..nw], det_precision);
    let c_dlog = measure(|| [0u64; 8], |mut d| {
        fp2_dlog_2e_pub::<L>(&mut d[..nw], &omega_aux_inv, &omega_f_inv, f as u32).unwrap();
        mp_mod_2exp(&mut d[..nw], det_precision);
        d
    });

    let det_hint = (cbytes[fp2_len] >> 2) & 0x03;
    let recovered = recover_dropped_entry::<L>(sig.mat(), &det[..nw], det_precision, det_hint);
    let m00_odd = sig.mat()[0][0].digits()[0] & 1 != 0;
    let dropped = if m00_odd { &sig.mat()[1][1] } else { &sig.mat()[1][0] };
    assert_eq!(recovered.digits(), dropped.digits(), "recovery must reproduce the dropped entry");
    let c_recover = measure(|| (), |_| recover_dropped_entry::<L>(sig.mat(), &det[..nw], det_precision, det_hint));

    // ---- Table.
    let rows = vec![
        row("decode standard signature (Signature::from_bytes)", 1, 0, c_decode_std),
        row("decode compressed signature (CompressedSignature::from_bytes)", 0, 1, c_decode_comp),
        row("matrix range check + verify_a(A_pk)", 1, 1, c_precheck),
        row("E_aux from A_aux (EcCurve::from_a)", 1, 2, c_from_a),
        row("canonical basis to-hint on E_aux, hint check (verify_canonical_hint)", 1, 1, c_tohint_aux),
        row("canonical basis to-hint on E_aux, decompress (basis_to_hint)", 0, 1, c_basis_to_hint_aux),
        row("challenge derivation E_pk -> E_chall (compute_challenge_curve)", 1, 2, c_chall),
        row("canonical basis to-hint on E_chall, hint check (verify_canonical_hint)", 1, 1, c_tohint_chall),
        row("canonical basis to-hint on E_chall, decompress (basis_to_hint)", 0, 1, c_basis_to_hint_chall),
        row("canonical basis from-hint on E_chall (basis_from_hint)", 1, 1, c_fromhint_chall),
        row("basis doubling on E_chall, verify (ec_dbl_iter_basis, F_CHR-N)", 1, 1, c_dbl_chall),
        row("basis doubling on E_chall, decompress (ec_dbl_iter_basis, F_CHR-f)", 0, 1, c_dbl_chall_dec),
        row("canonical basis from-hint on E_aux (basis_from_hint)", 1, 1, c_fromhint_aux),
        row("basis doubling on E_aux, verify (ec_dbl_iter_basis, F_CHR-g)", 1, 1, c_dbl_aux),
        row("basis doubling on E_aux, decompress (ec_dbl_iter_basis, F_CHR-g)", 0, 1, c_dbl_aux_dec),
        row("P+Q (xadd), per curve", 0, 2, c_xadd),
        row("Weil pairing on E_chall, 2^f", 0, 1, c_weil_chall),
        row("Weil pairing on E_aux, 2^g", 0, 1, c_weil_aux),
        row("two Fp2 inversions of the pairing values", 0, 1, c_inv),
        row("discrete logarithm (fp2_dlog_2e) + mask", 0, 1, c_dlog),
        row("modular-inverse entry recovery (Hensel + 2 products + hint bits)", 0, 1, c_recover),
        row("basis change, apply M (matrix_scalar_application_even_basis)", 1, 1, c_matrix),
        row("2^trl small isogeny chain (two_response_isogeny_verify_inner)", 1, 1, c_small),
        row("(2,2)-isogeny theta chain (compute_commitment_curve_verify)", 1, 1, c_theta),
        row("final hash + compare (hash_to_challenge, mp_compare)", 1, 1, c_hash),
    ];

    writeln!(w, "{:<76} {:>4} {:>4} {:>11} {:>11} {:>11} {:>11}", "component", "xS", "xC", "median us", "min us", "S us", "C us").unwrap();
    let mut sum_s_med = 0.0;
    let mut sum_c_med = 0.0;
    let mut sum_s_min = 0.0;
    let mut sum_c_min = 0.0;
    for r in &rows {
        let s_us = r.std_count as f64 * r.stat.median_us;
        let c_us = r.comp_count as f64 * r.stat.median_us;
        sum_s_med += s_us;
        sum_c_med += c_us;
        sum_s_min += r.std_count as f64 * r.stat.min_us;
        sum_c_min += r.comp_count as f64 * r.stat.min_us;
        writeln!(w, "{:<76} {:>4} {:>4} {:>11.1} {:>11.1} {:>11.1} {:>11.1}", r.name, r.std_count, r.comp_count, r.stat.median_us, r.stat.min_us, s_us, c_us).unwrap();
    }
    writeln!(w, "{:<76} {:>4} {:>4} {:>11} {:>11} {:>11.1} {:>11.1}", "sum of component medians", "", "", "", "", sum_s_med, sum_c_med).unwrap();
    writeln!(w, "{:<76} {:>4} {:>4} {:>11} {:>11} {:>11.1} {:>11.1}", "sum of component minima", "", "", "", "", sum_s_min, sum_c_min).unwrap();
    writeln!(w, "{:<76} {:>4} {:>4} {:>11} {:>11} {:>11.1} {:>11.1}", "FULL verification from bytes, median", "", "", "", "", full_std.median_us, full_comp.median_us).unwrap();
    writeln!(w, "{:<76} {:>4} {:>4} {:>11} {:>11} {:>11.1} {:>11.1}", "FULL verification from bytes, min", "", "", "", "", full_std.min_us, full_comp.min_us).unwrap();
    let res_s = full_std.median_us - sum_s_med;
    let res_c = full_comp.median_us - sum_c_med;
    writeln!(w, "residual (full median - component sum): standard {res_s:.1} us ({:.2}%), compressed {res_c:.1} us ({:.2}%)", 100.0 * res_s / full_std.median_us, 100.0 * res_c / full_comp.median_us).unwrap();
    writeln!(w, "checkpoints (median/min us): protocols_verify on decoded Signature {:.1}/{:.1}; decompress() alone {:.1}/{:.1}; PublicKey::from_bytes (not part of verify) {:.1}/{:.1}", pv_only.median_us, pv_only.min_us, decompress_only.median_us, decompress_only.min_us, pk_decode.median_us, pk_decode.min_us).unwrap();
    writeln!(w, "compressed - standard full medians: {:.1} us; decompress() alone: {:.1} us; sum of C-only component medians: {:.1} us", full_comp.median_us - full_std.median_us, decompress_only.median_us, sum_c_med - sum_s_med + c_decode_std.median_us - c_decode_comp.median_us).unwrap();
    out
}

fn full_only<L: FpBackend + LevelPrecomp>(name: &str, rsp: &str) -> String {
    let entries = kat_parser::parse_rsp(rsp);
    let e = &entries[0];
    assert_eq!(e.count, 0);
    let pk = PublicKey::<L>::from_bytes(&e.pk).unwrap();
    let sig_len = e.sm.len() - e.msg.len();
    let sig_bytes = &e.sm[..sig_len];
    let msg: &[u8] = &e.msg;
    let sig = Signature::<L>::from_bytes(sig_bytes).unwrap();
    let comp = sig.compress();
    let cbytes = comp.to_bytes();

    let full_std = measure(|| (), |_| {
        let s = Signature::<L>::from_bytes(sig_bytes).unwrap();
        pk.verify(msg, &s).unwrap()
    });
    let full_comp = measure(|| (), |_| {
        let c = CompressedSignature::<L>::from_bytes(cbytes.as_slice()).unwrap();
        pk.verify(msg, &c).unwrap()
    });
    let pv_only = measure(|| sig.clone(), |s| protocols_verify(&pk, msg, &s).unwrap());
    let decompress_only = measure(|| comp.clone(), |c| c.decompress(&pk).unwrap());
    let c_decode_std = measure(|| (), |_| Signature::<L>::from_bytes(sig_bytes).unwrap());
    let c_decode_comp = measure(|| (), |_| CompressedSignature::<L>::from_bytes(cbytes.as_slice()).unwrap());

    let mut out = String::new();
    let w = &mut out;
    writeln!(w, "==== {name}: KAT entry 0, bt={}, trl={}, sig {} bytes, compressed {} bytes, {ITERS} iterations ====", sig.backtracking(), sig.two_resp_length(), sig_len, cbytes.len()).unwrap();
    writeln!(w, "{:<50} {:>11} {:>11}", "measurement", "median us", "min us").unwrap();
    for (n, s) in [
        ("FULL standard verification from bytes", full_std),
        ("FULL compressed verification from bytes", full_comp),
        ("protocols_verify on decoded Signature", pv_only),
        ("decompress() alone", decompress_only),
        ("Signature::from_bytes", c_decode_std),
        ("CompressedSignature::from_bytes", c_decode_comp),
    ] {
        writeln!(w, "{:<50} {:>11.1} {:>11.1}", n, s.median_us, s.min_us).unwrap();
    }
    writeln!(w, "compressed - standard full medians: {:.1} us; decompress() alone: {:.1} us", full_comp.median_us - full_std.median_us, decompress_only.median_us).unwrap();
    out
}

#[test]
fn paper_timing_level1_breakdown() {
    print!("{}", breakdown::<Level1>("Level 1", KAT_L1));
}

#[test]
fn paper_timing_level3_full() {
    print!("{}", full_only::<Level3>("Level 3", KAT_L3));
}

#[test]
fn paper_timing_level5_full() {
    print!("{}", full_only::<Level5>("Level 5", KAT_L5));
}
