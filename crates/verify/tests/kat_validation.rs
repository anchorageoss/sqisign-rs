//! KAT (Known Answer Test) validation for SQIsign verification.
//!
//! Parses the NIST PQC KAT response files and verifies that the Rust
//! implementation produces identical results across all three security
//! levels and all signature formats (standard, expanded, compressed).
//!
//! These vectors come from the C reference, which does not canonicalize the
//! basis-change matrix, so ~half of them are non-canonical and are rejected by
//! the default (canonical) verifier. The suite therefore only compiles under
//! the `kat-compat` feature, which restores the legacy encoding. Run it with:
//!   cargo test -p sqisign-verify --release --features kat-compat --test kat_validation

#![cfg(feature = "kat-compat")]

use sqisign_verify::fp::FpBackend;
use sqisign_verify::params::{Level1, Level3, Level5, SecurityLevel};
use sqisign_verify::precomp::LevelPrecomp;
use sqisign_verify::{PublicKey, Signature, Verifier};

type L1 = Level1;

const L1_SIG_BYTES: usize = 148;
const L3_SIG_BYTES: usize = 224;
const L5_SIG_BYTES: usize = 292;

// KAT response files, embedded at compile time. Paths are anchored to the crate
// manifest dir (CARGO_MANIFEST_DIR) and factored into single constants so a
// path change only has to be made once.
const KAT_LVL1: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp"
));
const KAT_LVL3: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/KAT/PQCsignKAT_529_SQIsign_lvl3.rsp"
));
const KAT_LVL5: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/KAT/PQCsignKAT_701_SQIsign_lvl5.rsp"
));

fn hex_to_bytes(hex: &str) -> Vec<u8> {
    (0..hex.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).unwrap())
        .collect()
}

struct KatEntry {
    pk: Vec<u8>,
    sm: Vec<u8>,
    mlen: usize,
}

fn parse_kat_entries(content: &str, max_entries: usize) -> Vec<KatEntry> {
    let mut entries = Vec::new();
    let mut pk = None;
    let mut sm = None;
    let mut mlen = None;

    for line in content.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("mlen = ") {
            mlen = Some(rest.parse::<usize>().unwrap());
        } else if let Some(rest) = line.strip_prefix("pk = ") {
            pk = Some(hex_to_bytes(rest));
        } else if let Some(rest) = line.strip_prefix("sm = ") {
            sm = Some(hex_to_bytes(rest));
        }

        if pk.is_some() && sm.is_some() && mlen.is_some() {
            entries.push(KatEntry {
                pk: pk.take().unwrap(),
                sm: sm.take().unwrap(),
                mlen: mlen.take().unwrap(),
            });
            if entries.len() >= max_entries {
                break;
            }
        }
    }
    entries
}

fn verify_kat_entry(entry: &KatEntry, index: usize) -> bool {
    let pk = PublicKey::<L1>::from_bytes(&entry.pk)
        .unwrap_or_else(|_| panic!("KAT {}: failed to decode public key", index));

    assert!(
        entry.sm.len() >= L1_SIG_BYTES,
        "KAT {}: sm too short ({} < {})",
        index,
        entry.sm.len(),
        L1_SIG_BYTES
    );

    let sig_bytes = &entry.sm[..L1_SIG_BYTES];
    let msg = &entry.sm[L1_SIG_BYTES..];

    assert_eq!(
        msg.len(),
        entry.mlen,
        "KAT {}: message length mismatch",
        index
    );

    let sig = Signature::<L1>::from_bytes(sig_bytes)
        .unwrap_or_else(|_| panic!("KAT {}: failed to decode signature", index));

    pk.verify(msg, &sig).is_ok()
}

#[test]
fn test_kat_entry_0() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    assert!(!entries.is_empty(), "no KAT entries found");
    assert!(verify_kat_entry(&entries[0], 0), "KAT entry 0 failed");
}

#[test]
fn test_kat_entry_1() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 2);
    assert!(entries.len() >= 2, "not enough KAT entries");
    assert!(verify_kat_entry(&entries[1], 1), "KAT entry 1 failed");
}

#[test]
fn test_kat_entry_2() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 3);
    assert!(entries.len() >= 3, "not enough KAT entries");
    assert!(verify_kat_entry(&entries[2], 2), "KAT entry 2 failed");
}

#[test]
fn test_kat_entry_3() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 4);
    assert!(entries.len() >= 4, "not enough KAT entries");
    assert!(verify_kat_entry(&entries[3], 3), "KAT entry 3 failed");
}

#[test]
fn test_kat_entry_4() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 5);
    assert!(entries.len() >= 5, "not enough KAT entries");
    assert!(verify_kat_entry(&entries[4], 4), "KAT entry 4 failed");
}

// --- Expanded format tests ---

#[test]
fn test_expand_and_verify_kat_entry_0() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    let entry = &entries[0];

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();
    let msg = &entry.sm[L1_SIG_BYTES..];

    let expanded = sig
        .expand(&pk)
        .expect("failed to expand KAT entry 0 signature");

    pk.verify(msg, &expanded)
        .expect("expanded verification failed on KAT entry 0");
}

#[test]
fn test_expand_and_verify_all_kat_entries() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 5);

    for (i, entry) in entries.iter().enumerate() {
        let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
        let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();
        let msg = &entry.sm[L1_SIG_BYTES..];

        let expanded = sig
            .expand(&pk)
            .unwrap_or_else(|e| panic!("KAT {}: expand failed: {:?}", i, e));

        pk.verify(msg, &expanded)
            .unwrap_or_else(|e| panic!("KAT {}: expanded verify failed: {:?}", i, e));
    }
}

#[test]
fn test_expanded_serialization_roundtrip_kat() {
    use sqisign_verify::formats::ExpandedSignature;

    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    let entry = &entries[0];

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();
    let msg = &entry.sm[L1_SIG_BYTES..];

    let expanded = sig.expand(&pk).expect("expand failed");
    let wire = expanded.to_bytes();
    let decoded = ExpandedSignature::<L1>::from_bytes(&wire[..ExpandedSignature::<L1>::WIRE_BYTES])
        .expect("expanded from_bytes failed");

    pk.verify(msg, &decoded)
        .expect("verify after roundtrip failed");
}

#[test]
fn test_expanded_rejects_wrong_message() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    let entry = &entries[0];

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();

    let expanded = sig.expand(&pk).expect("expand failed");
    assert!(
        pk.verify(b"wrong message", &expanded).is_err(),
        "expanded verify should reject wrong message"
    );
}

#[test]
fn test_expanded_any_signature_dispatch() {
    use sqisign_verify::formats::{AnySignature, SignatureFormat};

    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    let entry = &entries[0];

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();
    let msg = &entry.sm[L1_SIG_BYTES..];

    let expanded = sig.expand(&pk).expect("expand failed");
    let any = AnySignature::Expanded(expanded);
    assert_eq!(any.format(), SignatureFormat::Expanded);
    pk.verify(msg, &any)
        .expect("AnySignature::verify failed for expanded");
}

// --- Negative tests ---

#[test]
fn test_corrupted_signature() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    let entry = &entries[0];

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let mut sm_corrupt = entry.sm.clone();
    sm_corrupt[10] ^= 0x01;

    let sig = Signature::<L1>::from_bytes(&sm_corrupt[..L1_SIG_BYTES]).unwrap();
    let msg = &sm_corrupt[L1_SIG_BYTES..];

    assert!(
        pk.verify(msg, &sig).is_err(),
        "corrupted signature should be rejected"
    );
}

#[test]
fn test_wrong_message() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    let entry = &entries[0];

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();
    let wrong_msg = b"this is not the right message";

    assert!(
        pk.verify(wrong_msg, &sig).is_err(),
        "wrong message should be rejected"
    );
}

#[test]
fn test_corrupted_public_key() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    let entry = &entries[0];

    let mut pk_corrupt = entry.pk.clone();
    pk_corrupt[5] ^= 0x01;

    // Corrupting the curve coefficient invalidates the hint_pk byte,
    // so from_bytes rejects the non-canonical encoding.
    assert!(
        PublicKey::<L1>::from_bytes(&pk_corrupt).is_err(),
        "corrupted public key should be rejected at decode"
    );
}

// --- Bit-flip mutation tests ---

#[test]
fn test_signature_bitflip_mutation() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    let entry = &entries[0];

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let sig_bytes = &entry.sm[..L1_SIG_BYTES];
    let msg = &entry.sm[L1_SIG_BYTES..];

    let orig_sig = Signature::<L1>::from_bytes(sig_bytes).unwrap();
    assert!(pk.verify(msg, &orig_sig).is_ok(), "original must verify");

    let mut accepted_bits = Vec::new();
    let total_bits = L1_SIG_BYTES * 8;
    for bit in 0..total_bits {
        let mut mutated = sig_bytes.to_vec();
        mutated[bit / 8] ^= 1 << (bit % 8);

        let accepted = Signature::<L1>::from_bytes(&mutated)
            .ok()
            .and_then(|sig| pk.verify(msg, &sig).ok())
            .is_some();

        if accepted {
            accepted_bits.push(bit);
        }
    }

    // Matrix entries are allocated (E_RSP + 9) / 8 bytes = 128 bits for L1.
    // The signer uses up to 128 bits but verification only checks the low
    // E_RSP bits, so the highest 2 bits per entry are unused. These are the
    // only bits that may survive a flip without detection.
    let fp2_bytes = 64usize;
    let mat_entry_bytes = 16usize;
    let mat_start = fp2_bytes + 2; // after Fp2 + backtracking + two_resp_length
    for &bit in &accepted_bits {
        let byte = bit / 8;
        let in_matrix = byte >= mat_start && byte < mat_start + 4 * mat_entry_bytes;
        assert!(
            in_matrix,
            "non-matrix bit flip at bit {} (byte {}) was accepted",
            bit, byte,
        );
        let offset_in_entry = (byte - mat_start) % mat_entry_bytes;
        assert!(
            offset_in_entry == mat_entry_bytes - 1,
            "non-high-byte bit flip at bit {} accepted",
            bit,
        );
    }
    assert!(
        accepted_bits.len() <= 8,
        "too many accepted bit flips: {:?}",
        accepted_bits,
    );
}

#[test]
fn test_public_key_bitflip_mutation() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    let entry = &entries[0];

    let pk_bytes = &entry.pk;
    let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();
    let msg = &entry.sm[L1_SIG_BYTES..];

    let orig_pk = PublicKey::<L1>::from_bytes(pk_bytes).unwrap();
    assert!(orig_pk.verify(msg, &sig).is_ok(), "original must verify");

    let mut rejected = 0u32;
    let total_bits = pk_bytes.len() * 8;
    for bit in 0..total_bits {
        let mut mutated = pk_bytes.to_vec();
        mutated[bit / 8] ^= 1 << (bit % 8);

        let accepted = PublicKey::<L1>::from_bytes(&mutated)
            .ok()
            .and_then(|pk| pk.verify(msg, &sig).ok())
            .is_some();

        if !accepted {
            rejected += 1;
        }
    }

    assert_eq!(
        rejected,
        total_bits as u32,
        "all {} single-bit flips must be rejected, but {} were accepted",
        total_bits,
        total_bits as u32 - rejected,
    );
}

// --- Full 100-entry KAT validation across all formats ---

fn verify_all_formats_generic<L: FpBackend + LevelPrecomp>(
    content: &str,
    sig_bytes: usize,
    label: &str,
) {
    let entries = parse_kat_entries(content, 100);
    assert_eq!(entries.len(), 100, "{}: expected 100 KAT entries", label);

    let mut pass_standard = 0u32;
    let mut pass_expanded = 0u32;
    let mut pass_compressed = 0u32;

    for (idx, entry) in entries.iter().enumerate() {
        let pk = PublicKey::<L>::from_bytes(&entry.pk)
            .unwrap_or_else(|_| panic!("{} KAT {}: pk decode failed", label, idx));
        let sig = Signature::<L>::from_bytes(&entry.sm[..sig_bytes])
            .unwrap_or_else(|_| panic!("{} KAT {}: sig decode failed", label, idx));
        let msg = &entry.sm[sig_bytes..];

        // Standard
        {
            assert!(
                pk.verify(msg, &sig).is_ok(),
                "{} KAT {} standard failed",
                label,
                idx
            );
            pass_standard += 1;
        }

        // Expanded
        {
            let expanded = sig
                .expand(&pk)
                .unwrap_or_else(|e| panic!("{} KAT {}: expand failed: {:?}", label, idx, e));
            pk.verify(msg, &expanded).unwrap_or_else(|e| {
                panic!("{} KAT {}: expanded verify failed: {:?}", label, idx, e)
            });
            pass_expanded += 1;
        }

        // Compressed verification
        {
            let compressed = sig.compress();

            pk.verify(msg, &compressed).unwrap_or_else(|e| {
                let pow_dim2 =
                    L::E_RSP as i32 - sig.two_resp_length() as i32 - sig.backtracking() as i32;
                panic!(
                    "{} KAT {}: compressed verify failed: {:?} (bt={}, trl={}, pow_dim2={})",
                    label,
                    idx,
                    e,
                    sig.backtracking(),
                    sig.two_resp_length(),
                    pow_dim2,
                )
            });
            pass_compressed += 1;
        }
    }

    assert_eq!(pass_standard, 100, "{}: standard", label);
    assert_eq!(pass_expanded, 100, "{}: expanded", label);
    assert_eq!(pass_compressed, 100, "{}: compressed", label);
}

#[test]
fn test_all_100_kats_all_formats_l1() {
    let content = KAT_LVL1;
    verify_all_formats_generic::<Level1>(content, L1_SIG_BYTES, "L1");
}

#[test]
fn test_all_100_kats_all_formats_l3() {
    let content = KAT_LVL3;
    verify_all_formats_generic::<Level3>(content, L3_SIG_BYTES, "L3");
}

#[test]
fn test_all_100_kats_all_formats_l5() {
    let content = KAT_LVL5;
    verify_all_formats_generic::<Level5>(content, L5_SIG_BYTES, "L5");
}

// --- Determinant bit-precision diagnostic ---

#[test]
fn diagnostic_det_bit_precision() {
    use hybrid_array::typenum::Unsigned;
    use sqisign_verify::ec::pairing::{fp2_dlog_2e_pub, weil};
    use sqisign_verify::ec::point::{ec_dbl_iter_basis, xadd};
    use sqisign_verify::ec::EcCurve;
    use sqisign_verify::theta::HD_EXTRA_TORSION;
    use sqisign_verify::verify::{basis_from_hint, compute_challenge_curve};

    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 100);
    assert_eq!(entries.len(), 100);

    let nw = <L1 as SecurityLevel>::MpLimbs::USIZE; // 4 limbs = 256 bits
    let wide = 2 * nw; // 8 limbs = 512 bits for full-precision products

    let mut match_pow_dim2 = 0u32;
    let mut match_g = 0u32;
    let mut match_f = 0u32;
    let mut mismatch_histogram = [0u32; 256];

    for (idx, entry) in entries.iter().enumerate() {
        let pk = PublicKey::<L1>::from_bytes(&entry.pk)
            .unwrap_or_else(|_| panic!("KAT {}: pk decode failed", idx));
        let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES])
            .unwrap_or_else(|_| panic!("KAT {}: sig decode failed", idx));

        let pow_dim2 =
            L1::E_RSP as usize - sig.two_resp_length() as usize - sig.backtracking() as usize;
        let g = pow_dim2 + HD_EXTRA_TORSION as usize;
        let f = pow_dim2 + HD_EXTRA_TORSION as usize + sig.two_resp_length() as usize;

        // --- Ground truth det(M) at full precision (512-bit intermediates) ---
        // det_true = M[0][0]*M[1][1] - M[0][1]*M[1][0]
        let mut prod1 = [0u64; 8]; // M00 * M11
        for i in 0..nw {
            let mut carry: u64 = 0;
            for j in 0..nw {
                if i + j >= wide {
                    break;
                }
                let p = (sig.mat()[0][0].digits()[i] as u128)
                    * (sig.mat()[1][1].digits()[j] as u128)
                    + (prod1[i + j] as u128)
                    + (carry as u128);
                prod1[i + j] = p as u64;
                carry = (p >> 64) as u64;
            }
            if i + nw < wide {
                prod1[i + nw] = carry;
            }
        }

        let mut prod2 = [0u64; 8]; // M01 * M10
        for i in 0..nw {
            let mut carry: u64 = 0;
            for j in 0..nw {
                if i + j >= wide {
                    break;
                }
                let p = (sig.mat()[0][1].digits()[i] as u128)
                    * (sig.mat()[1][0].digits()[j] as u128)
                    + (prod2[i + j] as u128)
                    + (carry as u128);
                prod2[i + j] = p as u64;
                carry = (p >> 64) as u64;
            }
            if i + nw < wide {
                prod2[i + nw] = carry;
            }
        }

        // det_true = prod1 - prod2 (full width, then mask to f bits)
        let mut det_true = [0u64; 8];
        {
            let mut borrow: u64 = 0;
            for i in 0..wide {
                let (d1, b1) = prod1[i].overflowing_sub(prod2[i]);
                let (d2, b2) = d1.overflowing_sub(borrow);
                det_true[i] = d2;
                borrow = (b1 as u64) + (b2 as u64);
            }
        }
        mask_bits(&mut det_true, f);

        // --- Pairing-based det at order 2^f ---
        let mut e_chall = compute_challenge_curve::<L1>(
            sig.chall_coeff(),
            sig.backtracking(),
            pk.curve(),
            pk.hint_pk(),
        )
        .unwrap_or_else(|| panic!("KAT {}: challenge curve failed", idx));

        let mut b_chall = basis_from_hint::<L1>(&mut e_chall, L1::F_CHR, sig.hint_chall())
            .unwrap_or_else(|| panic!("KAT {}: basis chall failed", idx));
        b_chall = ec_dbl_iter_basis(&b_chall, L1::F_CHR as usize - f, &mut e_chall);
        let ppq_chall = xadd(&b_chall.p, &b_chall.q, &b_chall.pmq);
        let omega_f = weil::<L1>(f as u32, &b_chall.p, &b_chall.q, &ppq_chall, &mut e_chall);

        let mut e_aux =
            EcCurve::<L1>::from_a(sig.e_aux_a()).unwrap_or_else(|| panic!("KAT {}: e_aux", idx));
        let mut b_aux = basis_from_hint::<L1>(&mut e_aux, L1::F_CHR, sig.hint_aux())
            .unwrap_or_else(|| panic!("KAT {}: basis aux failed", idx));
        b_aux = ec_dbl_iter_basis(&b_aux, L1::F_CHR as usize - g, &mut e_aux);
        let ppq_aux = xadd(&b_aux.p, &b_aux.q, &b_aux.pmq);
        let omega_aux = weil::<L1>(g as u32, &b_aux.p, &b_aux.q, &ppq_aux, &mut e_aux);

        let omega_f_inv = omega_f.inv();
        let omega_aux_inv = omega_aux.inv();
        let mut det_pair = [0u64; 8];
        fp2_dlog_2e_pub::<L1>(&mut det_pair[..nw], &omega_aux_inv, &omega_f_inv, f as u32)
            .unwrap_or_else(|| panic!("KAT {}: dlog failed", idx));
        mask_bits(&mut det_pair, f);

        // --- XOR to find exact mismatch position ---
        let mut diff = [0u64; 8];
        for i in 0..wide {
            diff[i] = det_true[i] ^ det_pair[i];
        }
        mask_bits(&mut diff, f);

        let first_mismatch = lowest_set_bit(&diff, wide);

        // Check match at each level
        let m_pdim2 = first_mismatch.map_or(true, |b| b >= pow_dim2);
        let m_g = first_mismatch.map_or(true, |b| b >= g);
        let m_f = first_mismatch.is_none();

        if m_pdim2 {
            match_pow_dim2 += 1;
        }
        if m_g {
            match_g += 1;
        }
        if m_f {
            match_f += 1;
        }

        let mm_str = match first_mismatch {
            Some(b) => {
                if b < 256 {
                    mismatch_histogram[b] += 1;
                }
                format!("{}", b)
            }
            None => "NONE (exact match)".to_string(),
        };

        eprintln!(
            "KAT {:3}: bt={} trl={} pow_dim2={:3} g={:3} f={:3} first_mismatch={:<20} match_pdim2={} match_g={} match_f={}",
            idx, sig.backtracking(), sig.two_resp_length(),
            pow_dim2, g, f, mm_str, m_pdim2, m_g, m_f
        );
    }

    eprintln!();
    eprintln!("=== SUMMARY ===");
    eprintln!("match at pow_dim2: {}/100", match_pow_dim2);
    eprintln!("match at g:        {}/100", match_g);
    eprintln!("match at f:        {}/100", match_f);
    eprintln!();

    eprintln!("first_mismatch_bit histogram:");
    for (bit, &count) in mismatch_histogram.iter().enumerate() {
        if count > 0 {
            eprintln!("  bit {:3}: {} entries", bit, count);
        }
    }
    let no_mismatch = 100 - mismatch_histogram.iter().sum::<u32>();
    if no_mismatch > 0 {
        eprintln!("  exact:   {} entries", no_mismatch);
    }

    eprintln!();
    if match_g == 100 && match_f < 100 {
        eprintln!("CONCLUSION: pairing gives det(M) mod 2^g (not just pow_dim2).");
        eprintln!("  Hint shrinks from (2+trl) bits to (trl) bits.");
        eprintln!("  For trl=0 entries ({}%), hint is unnecessary.", match_f);
    } else if match_f == 100 {
        eprintln!("CONCLUSION: pairing gives det(M) mod 2^f. No hint needed. 132 bytes.");
    } else if match_pow_dim2 == 100 {
        eprintln!("CONCLUSION: pairing gives det(M) mod 2^pow_dim2 only. Full hint needed.");
    } else {
        eprintln!("CONCLUSION: unexpected, some entries fail even at pow_dim2.");
    }

    assert_eq!(
        match_pow_dim2, 100,
        "pairing formula broken at pow_dim2 level"
    );
}

// ---------------------------------------------------------------------------
// Pairing / point order diagnostic.
//
// Prints, for the first 10 Level-1 KAT vectors:
//   Order(P_chall), Order(Q_chall)         (challenge basis, reduced to 2^f)
//   Order(omega_f), Order(omega_g)         (Weil pairings at exponents f, g)
//   order of the full-torsion Tate value   before and after ^(p-1)
//
// Number theory for Level 1 (p = 5*2^248 - 1):
//   p + 1   = 5 * 2^248                      (2- and 5-smooth)
//   p - 1   = 2 * (5*2^247 - 1)              (odd part large / unfactored)
//   p^2 - 1 = 2^249 * 5 * M,   M = 5*2^247-1 (coprime to 10)
// So an element's order after ^(p-1) is 2^a * 5^b (fully computable), while the
// raw value before ^(p-1) can also carry a factor of the unfactored M.
// ---------------------------------------------------------------------------

/// Smallest `k` with `x^(2^k) == 1`, or `None` if not reached within `cap`.
fn fp2_ord_2adic<L: FpBackend>(x: &sqisign_verify::fp::Fp2<L>, cap: u32) -> Option<u32> {
    if bool::from(x.ct_is_one()) {
        return Some(0);
    }
    let mut t = x.clone();
    for k in 1..=cap {
        t = t.sqr();
        if bool::from(t.ct_is_one()) {
            return Some(k);
        }
    }
    None
}

/// Multiplicative order of a 2-power torsion point, as the exponent `k` in
/// `2^k` (smallest `k` with `[2^k]P = O`), or `None` if not reached.
fn point_ord_2adic<L: FpBackend>(
    p: &sqisign_verify::ec::EcPoint<L>,
    curve: &sqisign_verify::ec::EcCurve<L>,
    cap: u32,
) -> Option<u32> {
    use sqisign_verify::ec::point::ec_dbl;
    if bool::from(p.is_zero()) {
        return Some(0);
    }
    let mut t = p.clone();
    for k in 1..=cap {
        t = ec_dbl(&t, curve);
        if bool::from(t.is_zero()) {
            return Some(k);
        }
    }
    None
}

/// Order of `x` given `ord(x) | (p+1) = 5 * 2^248`. Returns `(a, b)` for
/// `ord = 2^a * 5^b`.
fn order_in_p_plus_1<L: FpBackend>(x: &sqisign_verify::fp::Fp2<L>) -> (u32, u32) {
    const E_2F: &[u64] = &[0, 0, 0, 72057594037927936]; // 2^248
    // 5-part: raise to 2^248 to kill the 2-part; remainder has order | 5.
    let after2 = x.pow_vartime(E_2F);
    let b = if bool::from(after2.ct_is_one()) { 0 } else { 1 };
    // 2-part: raise to 5 to kill the 5-part, then measure the 2-adic order.
    let after5 = x.pow_vartime(&[5]);
    let a = fp2_ord_2adic(&after5, 250).expect("2-part exceeds 2^250");
    (a, b)
}

/// Order of `x` given `ord(x) | (p^2-1) = 2^249 * 5 * M`, `M = 5*2^247-1`
/// coprime to 10. Returns `(a, b, m_nontrivial)` for
/// `ord = 2^a * 5^b * m` with `m | M`, `m > 1` iff `m_nontrivial`.
fn order_in_p2_minus_1<L: FpBackend>(x: &sqisign_verify::fp::Fp2<L>) -> (u32, u32, bool) {
    // (p^2-1)/M = 2^249 * 5
    const E_M: &[u64] = &[0, 0, 0, 720575940379279360];
    // (p^2-1)/2^249 = 5*M  (odd)
    const E_2: &[u64] = &[
        18446744073709551611,
        18446744073709551615,
        18446744073709551615,
        900719925474099199,
    ];
    // (p^2-1)/5 = 2^249 * M
    const E_5: &[u64] = &[
        0,
        0,
        0,
        18302628885633695744,
        18446744073709551615,
        18446744073709551615,
        18446744073709551615,
        1407374883553279,
    ];
    // M-part trivial iff x^((p^2-1)/M) == 1.
    let xm = x.pow_vartime(E_M);
    let m_nontrivial = !bool::from(xm.ct_is_one());
    // 2-part: raise to (p^2-1)/2^249, remainder is the 2-Sylow component.
    let u = x.pow_vartime(E_2);
    let a = fp2_ord_2adic(&u, 250).expect("2-part exceeds 2^250");
    // 5-part.
    let v = x.pow_vartime(E_5);
    let b = if bool::from(v.ct_is_one()) { 0 } else { 1 };
    (a, b, m_nontrivial)
}

#[test]
fn diagnostic_pairing_orders() {
    use sqisign_verify::ec::pairing::{reduced_tate_stages, weil};
    use sqisign_verify::ec::point::{ec_dbl_iter_basis, xadd};
    use sqisign_verify::ec::EcCurve;
    use sqisign_verify::theta::HD_EXTRA_TORSION;
    use sqisign_verify::verify::{basis_from_hint, compute_challenge_curve};

    let f_chr = L1::F_CHR; // 248
    const COFAC: &[u64] = &[5]; // (p+1)/2^248

    let entries = parse_kat_entries(KAT_LVL1, 10);
    assert_eq!(entries.len(), 10);

    eprintln!("Level 1: p = 5*2^248 - 1 ; p+1 = 5*2^248 ; p-1 = 2*(5*2^247-1)");
    eprintln!("F_CHR = {} (full 2-power torsion). Orders shown as 2^a[*5^b].\n", f_chr);

    let fmt = |o: Option<u32>| match o {
        Some(k) => format!("2^{}", k),
        None => ">2^cap (unexpected)".to_string(),
    };

    for (idx, entry) in entries.iter().enumerate() {
        let pk = PublicKey::<L1>::from_bytes(&entry.pk)
            .unwrap_or_else(|_| panic!("KAT {}: pk decode failed", idx));
        let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES])
            .unwrap_or_else(|_| panic!("KAT {}: sig decode failed", idx));

        let pow_dim2 =
            L1::E_RSP as usize - sig.two_resp_length() as usize - sig.backtracking() as usize;
        let g = pow_dim2 + HD_EXTRA_TORSION as usize;
        let f = pow_dim2 + HD_EXTRA_TORSION as usize + sig.two_resp_length() as usize;

        // --- challenge curve + full 2^F_CHR torsion basis ---
        let mut e_chall = compute_challenge_curve::<L1>(
            sig.chall_coeff(),
            sig.backtracking(),
            pk.curve(),
            pk.hint_pk(),
        )
        .unwrap_or_else(|| panic!("KAT {}: challenge curve failed", idx));
        let b_chall_full = basis_from_hint::<L1>(&mut e_chall, f_chr, sig.hint_chall())
            .unwrap_or_else(|| panic!("KAT {}: basis chall failed", idx));

        // Challenge basis reduced to order 2^f (the points that feed omega_f).
        let b_chall_f = ec_dbl_iter_basis(&b_chall_full, f_chr as usize - f, &mut e_chall);
        let ppq_f = xadd(&b_chall_f.p, &b_chall_f.q, &b_chall_f.pmq);
        let omega_f = weil::<L1>(f as u32, &b_chall_f.p, &b_chall_f.q, &ppq_f, &mut e_chall);

        // --- aux curve + basis; omega_g (=omega_aux) at exponent g ---
        let mut e_aux =
            EcCurve::<L1>::from_a(sig.e_aux_a()).unwrap_or_else(|| panic!("KAT {}: e_aux", idx));
        let b_aux_full = basis_from_hint::<L1>(&mut e_aux, f_chr, sig.hint_aux())
            .unwrap_or_else(|| panic!("KAT {}: basis aux failed", idx));
        let b_aux_g = ec_dbl_iter_basis(&b_aux_full, f_chr as usize - g, &mut e_aux);
        let ppq_g = xadd(&b_aux_g.p, &b_aux_g.q, &b_aux_g.pmq);
        let omega_g = weil::<L1>(g as u32, &b_aux_g.p, &b_aux_g.q, &ppq_g, &mut e_aux);

        // --- full-torsion (e = F_CHR) reduced Tate on the challenge basis,
        //     exposing the raw value and the value after ^(p-1) ---
        let ppq_full = xadd(&b_chall_full.p, &b_chall_full.q, &b_chall_full.pmq);
        let (raw, after_pm1, reduced) = reduced_tate_stages::<L1>(
            f_chr,
            &b_chall_full.p,
            &b_chall_full.q,
            &ppq_full,
            &mut e_chall,
            f_chr,
            COFAC,
        );

        // --- orders (e_chall / e_aux are a24-normalized by the pairing calls) ---
        let ord_p = point_ord_2adic(&b_chall_f.p, &e_chall, 260);
        let ord_q = point_ord_2adic(&b_chall_f.q, &e_chall, 260);
        let ord_p_full = point_ord_2adic(&b_chall_full.p, &e_chall, 260);
        let ord_omega_f = fp2_ord_2adic(&omega_f, 260);
        let ord_omega_g = fp2_ord_2adic(&omega_g, 260);
        let ord_reduced = fp2_ord_2adic(&reduced, 260);

        let (a_raw, b_raw, m_raw) = order_in_p2_minus_1(&raw);
        let (a_aft, b_aft) = order_in_p_plus_1(&after_pm1);

        let smooth = |a: u32, b: u32| {
            if b > 0 {
                format!("2^{} * 5^{}", a, b)
            } else {
                format!("2^{}", a)
            }
        };
        let raw_str = {
            let mut s = smooth(a_raw, b_raw);
            if m_raw {
                s.push_str(" * m   (m>1, m | (5*2^247-1), not factored)");
            }
            s
        };

        eprintln!("KAT {:2}:  f={}  g={}  pow_dim2={}", idx, f, g, pow_dim2);
        eprintln!("   Order(P_chall)            = {}   (2^f basis)", fmt(ord_p));
        eprintln!("   Order(Q_chall)            = {}   (2^f basis)", fmt(ord_q));
        eprintln!("   Order(omega_f)  [Weil,e=f]= {}", fmt(ord_omega_f));
        eprintln!("   Order(omega_g)  [Weil,e=g]= {}", fmt(ord_omega_g));
        eprintln!("   full-torsion basis pt ord = {}   (sanity, expect 2^{})", fmt(ord_p_full), f_chr);
        eprintln!("   Tate value BEFORE ^(p-1)  = {}", raw_str);
        eprintln!("   Tate value AFTER  ^(p-1)  = {}", smooth(a_aft, b_aft));
        eprintln!("   Tate fully reduced        = {}   (expect 2^{})", fmt(ord_reduced), f_chr);
        eprintln!();
    }
}

fn mask_bits(a: &mut [u64], nbits: usize) {
    let q = nbits / 64;
    let r = nbits % 64;
    if q < a.len() {
        if r != 0 {
            a[q] &= (1u64 << r) - 1;
        } else {
            a[q] = 0;
        }
        for limb in a[q + 1..].iter_mut() {
            *limb = 0;
        }
    }
}

fn lowest_set_bit(a: &[u64], nlimbs: usize) -> Option<usize> {
    for (i, &limb) in a.iter().enumerate().take(nlimbs) {
        if limb != 0 {
            return Some(i * 64 + limb.trailing_zeros() as usize);
        }
    }
    None
}

// ===========================================================================
// Paulo's short-ladder Tate experiment.
//
// Idea: replace the two Weil ladders (2 x ~f biextension steps) with a single
// order-2^f Tate ladder (~f steps) plus a *partial* final exponentiation
// ^((p-1)*c), where c = (p+1)/2^F_CHR is the small odd cofactor (5 @ L1,
// 65 @ L3, 27 @ L5). This deliberately omits the (F_CHR - f) squarings that
// the standard reduction ^((p^2-1)/2^f) would apply.
// ===========================================================================

/// Ground-truth det(M) = M00*M11 - M01*M10 at full precision, into `out`
/// (needs `out.len() >= 2*nw`). Not reduced mod 2^f (caller masks).
fn ground_truth_det<L: FpBackend>(sig: &Signature<L>, nw: usize, out: &mut [u64]) {
    let wide = 2 * nw;
    let mut prod1 = [0u64; 16]; // M00 * M11
    for i in 0..nw {
        let mut carry: u64 = 0;
        for j in 0..nw {
            if i + j >= wide {
                break;
            }
            let p = (sig.mat()[0][0].digits()[i] as u128) * (sig.mat()[1][1].digits()[j] as u128)
                + (prod1[i + j] as u128)
                + (carry as u128);
            prod1[i + j] = p as u64;
            carry = (p >> 64) as u64;
        }
        if i + nw < wide {
            prod1[i + nw] = carry;
        }
    }
    let mut prod2 = [0u64; 16]; // M01 * M10
    for i in 0..nw {
        let mut carry: u64 = 0;
        for j in 0..nw {
            if i + j >= wide {
                break;
            }
            let p = (sig.mat()[0][1].digits()[i] as u128) * (sig.mat()[1][0].digits()[j] as u128)
                + (prod2[i + j] as u128)
                + (carry as u128);
            prod2[i + j] = p as u64;
            carry = (p >> 64) as u64;
        }
        if i + nw < wide {
            prod2[i + nw] = carry;
        }
    }
    let mut borrow: u64 = 0;
    for i in 0..wide {
        let (d1, b1) = prod1[i].overflowing_sub(prod2[i]);
        let (d2, b2) = d1.overflowing_sub(borrow);
        out[i] = d2;
        borrow = (b1 as u64) + (b2 as u64);
    }
}

/// First differing bit (below `nbits`) of two little-endian integers, or None
/// if they agree on all `nbits` low bits.
fn first_mismatch(a: &[u64], b: &[u64], nbits: usize) -> Option<usize> {
    let mut d = [0u64; 16];
    let n = a.len().min(b.len()).min(16);
    for i in 0..n {
        d[i] = a[i] ^ b[i];
    }
    mask_bits(&mut d, nbits);
    lowest_set_bit(&d, n)
}

/// One order-2^e Tate ladder + partial exponentiation ^((p-1)*c). Returns
/// `(reduced, paulo)`: `reduced` is the standard order-2^e reduced Tate value
/// (= `paulo` squared down (F_CHR - e) times); `paulo = raw^((p-1)*c)` has
/// order 2^F_CHR.
fn short_tate<L: FpBackend + LevelPrecomp>(
    e: u32,
    p: &sqisign_verify::ec::EcPoint<L>,
    q: &sqisign_verify::ec::EcPoint<L>,
    ppq: &sqisign_verify::ec::EcPoint<L>,
    curve: &mut sqisign_verify::ec::EcCurve<L>,
) -> (sqisign_verify::fp::Fp2<L>, sqisign_verify::fp::Fp2<L>) {
    use sqisign_verify::ec::pairing::{clear_cofac, reduced_tate_stages};
    let f_chr = <L as SecurityLevel>::F_CHR;
    let cofac = L::p_cofactor_for_2f();
    // reduced_tate_stages returns (raw, raw^(p-1), fully_reduced).
    // Return (raw, paulo) where paulo = raw^((p-1)*c).
    let (raw, pm1, _) = reduced_tate_stages::<L>(e, p, q, ppq, curve, f_chr, cofac);
    let paulo = clear_cofac::<L>(&pm1, cofac);
    (raw, paulo)
}

/// Generic det-recovery experiment for one security level.
fn paulo_det_recovery<L: FpBackend + LevelPrecomp>(kat: &str, sig_bytes: usize, n: usize, label: &str) {
    use hybrid_array::typenum::Unsigned;
    use sqisign_verify::ec::pairing::{clear_cofac, fp2_dlog_2e_pub, weil};
    use sqisign_verify::ec::point::{ec_dbl_iter_basis, xadd};
    use sqisign_verify::ec::EcCurve;
    use sqisign_verify::theta::HD_EXTRA_TORSION;
    use sqisign_verify::verify::{basis_from_hint, compute_challenge_curve};

    let nw = <L as SecurityLevel>::MpLimbs::USIZE;
    let f_chr = <L as SecurityLevel>::F_CHR;
    let cofac = L::p_cofactor_for_2f();

    let entries = parse_kat_entries(kat, n);
    let total = entries.len();

    // Reference (Weil) match vs ground truth, at the three bit precisions.
    let (mut weil_pd, mut weil_g, mut weil_f) = (0u32, 0u32, 0u32);
    // Single short-Tate ladder (Paulo's proposal), depth f, both arg orders,
    // counting exact agreement with the Weil-recovered det (mod 2^f).
    let (mut single_pq_eq, mut single_qp_eq) = (0u32, 0u32);
    // pow_dim2-level match vs ground truth for the single-ladder (best order).
    let mut single_pd = 0u32;
    // Ratio of two short-Tate ladders  t(Q,P)/t(P,Q)  (antisymmetric = Weil),
    // depth f: match vs ground truth and exact agreement with Weil.
    let (mut ratio_pd, mut ratio_g, mut ratio_f, mut ratio_eq) = (0u32, 0u32, 0u32, 0u32);
    // Relationships: is the single value == omega_f^{±c}?  is the ratio == omega_f^{±c}?
    let (mut rel_single_pos, mut rel_single_neg) = (0u32, 0u32);
    let (mut rel_ratio_pos, mut rel_ratio_neg) = (0u32, 0u32);
    // RAW-monodromy ratio  raw(Q,P)/raw(P,Q)  (no (p-1) reduction) = Weil?
    let (mut rawratio_pd, mut rawratio_g, mut rawratio_f, mut rawratio_eq) =
        (0u32, 0u32, 0u32, 0u32);
    let (mut rel_rawratio_pos, mut rel_rawratio_neg) = (0u32, 0u32);

    for entry in entries.iter() {
        let pk = PublicKey::<L>::from_bytes(&entry.pk).unwrap();
        let sig = Signature::<L>::from_bytes(&entry.sm[..sig_bytes]).unwrap();

        let pow_dim2 =
            L::E_RSP as usize - sig.two_resp_length() as usize - sig.backtracking() as usize;
        let g = pow_dim2 + HD_EXTRA_TORSION as usize;
        let f = pow_dim2 + HD_EXTRA_TORSION as usize + sig.two_resp_length() as usize;

        let mut det_true = [0u64; 16];
        ground_truth_det::<L>(&sig, nw, &mut det_true);
        mask_bits(&mut det_true, f);

        // --- challenge basis (reduced to 2^f) and omega_f ---
        let mut e_chall = compute_challenge_curve::<L>(
            sig.chall_coeff(),
            sig.backtracking(),
            pk.curve(),
            pk.hint_pk(),
        )
        .unwrap();
        let b_chall_full = basis_from_hint::<L>(&mut e_chall, f_chr, sig.hint_chall()).unwrap();
        let b_chall_f = ec_dbl_iter_basis(&b_chall_full, f_chr as usize - f, &mut e_chall);
        let ppq_f = xadd(&b_chall_f.p, &b_chall_f.q, &b_chall_f.pmq);
        let omega_f = weil::<L>(f as u32, &b_chall_f.p, &b_chall_f.q, &ppq_f, &mut e_chall);

        // --- aux basis (reduced to 2^g) and omega_g ---
        let mut e_aux = EcCurve::<L>::from_a(sig.e_aux_a()).unwrap();
        let b_aux_full = basis_from_hint::<L>(&mut e_aux, f_chr, sig.hint_aux()).unwrap();
        let b_aux_g = ec_dbl_iter_basis(&b_aux_full, f_chr as usize - g, &mut e_aux);
        let ppq_g = xadd(&b_aux_g.p, &b_aux_g.q, &b_aux_g.pmq);
        let omega_g = weil::<L>(g as u32, &b_aux_g.p, &b_aux_g.q, &ppq_g, &mut e_aux);

        // --- reference Weil det ---
        let mut det_weil = [0u64; 16];
        let _ = fp2_dlog_2e_pub::<L>(&mut det_weil[..nw], &omega_g.inv(), &omega_f.inv(), f as u32);
        mask_bits(&mut det_weil, f);
        match first_mismatch(&det_true, &det_weil, f) {
            None => {
                weil_pd += 1;
                weil_g += 1;
                weil_f += 1;
            }
            Some(b) => {
                if b >= pow_dim2 {
                    weil_pd += 1;
                }
                if b >= g {
                    weil_g += 1;
                }
            }
        }

        // --- short Tate values, BOTH argument orders. Each is one order-2^f
        // ladder + ^((p-1)*c); empirically the result already has order 2^f
        // (no squaring down). ---
        let (raw_c_pq, tate_c_pq) =
            short_tate::<L>(f as u32, &b_chall_f.p, &b_chall_f.q, &ppq_f, &mut e_chall);
        let (raw_c_qp, tate_c_qp) =
            short_tate::<L>(f as u32, &b_chall_f.q, &b_chall_f.p, &ppq_f, &mut e_chall);
        let (raw_a_pq, tate_a_pq) =
            short_tate::<L>(g as u32, &b_aux_g.p, &b_aux_g.q, &ppq_g, &mut e_aux);
        let (raw_a_qp, tate_a_qp) =
            short_tate::<L>(g as u32, &b_aux_g.q, &b_aux_g.p, &ppq_g, &mut e_aux);

        // (1) Paulo's proposal: SINGLE ladder, dlog at depth f, both orders.
        let mut det_pq = [0u64; 16];
        let _ = fp2_dlog_2e_pub::<L>(&mut det_pq[..nw], &tate_a_pq.inv(), &tate_c_pq.inv(), f as u32);
        mask_bits(&mut det_pq, f);
        let mut det_qp = [0u64; 16];
        let _ = fp2_dlog_2e_pub::<L>(&mut det_qp[..nw], &tate_a_qp.inv(), &tate_c_qp.inv(), f as u32);
        mask_bits(&mut det_qp, f);
        if first_mismatch(&det_weil, &det_pq, f).is_none() {
            single_pq_eq += 1;
        }
        if first_mismatch(&det_weil, &det_qp, f).is_none() {
            single_qp_eq += 1;
        }
        let single_best = first_mismatch(&det_true, &det_pq, f)
            .map_or(f, |b| b)
            .max(first_mismatch(&det_true, &det_qp, f).map_or(f, |b| b));
        if single_best >= pow_dim2 {
            single_pd += 1;
        }

        // (2) Antisymmetric RATIO of two ladders  t(Q,P)/t(P,Q)  (= Weil), depth f.
        let ratio_c = tate_c_qp.mul(&tate_c_pq.inv());
        let ratio_a = tate_a_qp.mul(&tate_a_pq.inv());
        let mut det_r = [0u64; 16];
        let _ = fp2_dlog_2e_pub::<L>(&mut det_r[..nw], &ratio_a.inv(), &ratio_c.inv(), f as u32);
        mask_bits(&mut det_r, f);
        match first_mismatch(&det_true, &det_r, f) {
            None => {
                ratio_pd += 1;
                ratio_g += 1;
                ratio_f += 1;
            }
            Some(bit) => {
                if bit >= pow_dim2 {
                    ratio_pd += 1;
                }
                if bit >= g {
                    ratio_g += 1;
                }
            }
        }
        if first_mismatch(&det_weil, &det_r, f).is_none() {
            ratio_eq += 1;
        }

        // (3) RAW-monodromy ratio  raw(Q,P)/raw(P,Q)  with NO (p-1) reduction.
        // This is what the Weil pairing actually is (ratio of the two cubical
        // ladders). Should reproduce Weil and recover det.
        let rawratio_c = raw_c_qp.mul(&raw_c_pq.inv());
        let rawratio_a = raw_a_qp.mul(&raw_a_pq.inv());
        let mut det_rr = [0u64; 16];
        let _ =
            fp2_dlog_2e_pub::<L>(&mut det_rr[..nw], &rawratio_a.inv(), &rawratio_c.inv(), f as u32);
        mask_bits(&mut det_rr, f);
        match first_mismatch(&det_true, &det_rr, f) {
            None => {
                rawratio_pd += 1;
                rawratio_g += 1;
                rawratio_f += 1;
            }
            Some(bit) => {
                if bit >= pow_dim2 {
                    rawratio_pd += 1;
                }
                if bit >= g {
                    rawratio_g += 1;
                }
            }
        }
        if first_mismatch(&det_weil, &det_rr, f).is_none() {
            rawratio_eq += 1;
        }

        // Relationships to the Weil value omega_f (order 2^f).
        let omega_f_c = clear_cofac::<L>(&omega_f, cofac);
        let omega_f_neg_c = clear_cofac::<L>(&omega_f.inv(), cofac);
        if bool::from(tate_c_pq.ct_equal(&omega_f_c)) {
            rel_single_pos += 1;
        }
        if bool::from(tate_c_pq.ct_equal(&omega_f_neg_c)) {
            rel_single_neg += 1;
        }
        if bool::from(ratio_c.ct_equal(&omega_f_c)) {
            rel_ratio_pos += 1;
        }
        if bool::from(ratio_c.ct_equal(&omega_f_neg_c)) {
            rel_ratio_neg += 1;
        }
        if bool::from(rawratio_c.ct_equal(&omega_f)) {
            rel_rawratio_pos += 1;
        }
        if bool::from(rawratio_c.ct_equal(&omega_f.inv())) {
            rel_rawratio_neg += 1;
        }
    }

    eprintln!("================ {} ({} vectors) ================", label, total);
    eprintln!(
        "Weil (reference, 2 ladders)  vs ground truth: pow_dim2 {}/{}  g {}/{}  f {}/{}",
        weil_pd, total, weil_g, total, weil_f, total
    );
    eprintln!(
        "Paulo SINGLE Tate ladder     exact-match w/ Weil det: (P,Q) {}/{}  (Q,P) {}/{}  | pow_dim2 vs truth {}/{}",
        single_pq_eq, total, single_qp_eq, total, single_pd, total
    );
    eprintln!(
        "REDUCED-Tate ratio t(Q,P)/t(P,Q) [with ^(p-1)c]  vs truth: pow_dim2 {}/{}  g {}/{}  f {}/{}  | ==Weil {}/{}",
        ratio_pd, total, ratio_g, total, ratio_f, total, ratio_eq, total
    );
    eprintln!(
        "RAW-monodromy ratio raw(Q,P)/raw(P,Q) [no reduction]  vs truth: pow_dim2 {}/{}  g {}/{}  f {}/{}  | ==Weil {}/{}",
        rawratio_pd, total, rawratio_g, total, rawratio_f, total, rawratio_eq, total
    );
    eprintln!(
        "Relationships: single==omega_f^±c {}/{},{}  reduced-ratio==omega_f^±c {}/{},{}  raw-ratio==omega_f^±1 {}/{},{}",
        rel_single_pos, rel_single_neg, total,
        rel_ratio_pos, rel_ratio_neg, total,
        rel_rawratio_pos, rel_rawratio_neg, total
    );
    eprintln!();
}

#[test]
fn paulo_short_tate_l1() {
    paulo_det_recovery::<Level1>(KAT_LVL1, L1_SIG_BYTES, 100, "Level 1  (p=5*2^248-1, c=5)");
}

#[test]
fn paulo_short_tate_l3() {
    paulo_det_recovery::<Level3>(KAT_LVL3, L3_SIG_BYTES, 100, "Level 3  (p=65*2^376-1, c=65)");
}

#[test]
fn paulo_short_tate_l5() {
    paulo_det_recovery::<Level5>(KAT_LVL5, L5_SIG_BYTES, 100, "Level 5  (p=27*2^500-1, c=27)");
}

/// Step 1: order of the order-2^f raw Tate value and of `raw^((p-1)*c)`, for
/// the first 10 Level-1 vectors. (Exact-order helpers are L1-specific.)
#[test]
fn paulo_short_tate_orders_l1() {
    use sqisign_verify::ec::pairing::{clear_cofac, reduced_tate_stages};
    use sqisign_verify::ec::point::{ec_dbl_iter_basis, xadd};
    use sqisign_verify::theta::HD_EXTRA_TORSION;
    use sqisign_verify::verify::{basis_from_hint, compute_challenge_curve};

    let f_chr = L1::F_CHR;
    let cofac: &[u64] = &[5];
    let entries = parse_kat_entries(KAT_LVL1, 10);

    eprintln!("Order at 2^f: raw Tate ladder output and raw^((p-1)*c). F_CHR={}\n", f_chr);
    for (idx, entry) in entries.iter().enumerate() {
        let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
        let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();
        let pow_dim2 =
            L1::E_RSP as usize - sig.two_resp_length() as usize - sig.backtracking() as usize;
        let f = pow_dim2 + HD_EXTRA_TORSION as usize + sig.two_resp_length() as usize;

        let mut e_chall = compute_challenge_curve::<L1>(
            sig.chall_coeff(),
            sig.backtracking(),
            pk.curve(),
            pk.hint_pk(),
        )
        .unwrap();
        let b_full = basis_from_hint::<L1>(&mut e_chall, f_chr, sig.hint_chall()).unwrap();
        let b_f = ec_dbl_iter_basis(&b_full, f_chr as usize - f, &mut e_chall);
        let ppq = xadd(&b_f.p, &b_f.q, &b_f.pmq);

        let (raw, pm1, reduced) =
            reduced_tate_stages::<L1>(f as u32, &b_f.q, &b_f.p, &ppq, &mut e_chall, f_chr, cofac);
        let paulo = clear_cofac::<L1>(&pm1, cofac);

        // cross-checks: Weil value at 2^f, and the fully-reduced (.2) Tate value
        let omega = {
            use sqisign_verify::ec::pairing::weil;
            weil::<L1>(f as u32, &b_f.p, &b_f.q, &ppq, &mut e_chall)
        };
        let (a_red, b_red) = order_in_p_plus_1(&reduced);
        let (a_om, b_om) = order_in_p_plus_1(&omega);

        let (a_raw, b_raw, m_raw) = order_in_p2_minus_1(&raw);
        let (a_p, b_p) = order_in_p_plus_1(&paulo);
        let _ = (b_red, b_om);
        let raw_s = {
            let mut s = if b_raw > 0 {
                format!("2^{} * 5^{}", a_raw, b_raw)
            } else {
                format!("2^{}", a_raw)
            };
            if m_raw {
                s.push_str(" * m (m>1, m|(5*2^247-1))");
            }
            s
        };
        let paulo_s = if b_p > 0 {
            format!("2^{} * 5^{}", a_p, b_p)
        } else {
            format!("2^{}", a_p)
        };
        eprintln!(
            "KAT {:2}: f={}  raw={}  paulo=raw^((p-1)c)={}  fully_reduced(.2)=2^{}  weil(2^f)=2^{}",
            idx, f, raw_s, paulo_s, a_red, a_om
        );
    }
}

// --- Matrix bit structure diagnostic ---

#[test]
fn diagnostic_matrix_bit_structure() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 100);

    use hybrid_array::typenum::Unsigned;

    let e_rsp = L1::E_RSP as usize;
    let nw = <L1 as SecurityLevel>::MpLimbs::USIZE;

    let mut count_leading_one = 0u32;
    let mut count_unique = 0u32;
    let mut count_ambiguous = 0u32;
    let mut bt_histogram = [0u32; 4];
    let mut rp_histogram = [0u32; 16];
    let mut high_bits_width_histogram = [0u32; 16];

    for (idx, entry) in entries.iter().enumerate() {
        let _pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
        let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();

        let n_bt = sig.backtracking() as usize;
        let r_prime = sig.two_resp_length() as usize;
        let pow_dim2 = e_rsp - r_prime - n_bt;
        let high_bits_count = r_prime + n_bt;

        bt_histogram[n_bt] += 1;
        if r_prime < rp_histogram.len() {
            rp_histogram[r_prime] += 1;
        }
        if high_bits_count < high_bits_width_histogram.len() {
            high_bits_width_histogram[high_bits_count] += 1;
        }

        // Extract M[0][0] and M[0][1] as digit arrays
        let m00 = sig.mat()[0][0].digits();
        let m01 = sig.mat()[0][1].digits();
        let m10 = sig.mat()[1][0].digits();
        let m11 = sig.mat()[1][1].digits();

        // Compute actual bit widths
        let bit_width = |digits: &[u64]| -> usize {
            for i in (0..nw).rev() {
                if digits[i] != 0 {
                    return i * 64 + (64 - digits[i].leading_zeros() as usize);
                }
            }
            0
        };

        let m00_bw = bit_width(m00);
        let m01_bw = bit_width(m01);
        let m10_bw = bit_width(m10);
        let m11_bw = bit_width(m11);

        // Extract high bits of M[0][0] and M[0][1] (shift right by pow_dim2)
        let extract_high = |digits: &[u64], shift: usize| -> u64 {
            let word = shift / 64;
            let bit = shift % 64;
            if word >= nw {
                return 0;
            }
            let lo = digits[word] >> bit;
            if bit > 0 && word + 1 < nw {
                lo | (digits[word + 1] << (64 - bit))
            } else {
                lo
            }
        };

        let m00_high = extract_high(m00, pow_dim2);
        let m01_high = extract_high(m01, pow_dim2);
        let m10_high = extract_high(m10, pow_dim2);
        let m11_high = extract_high(m11, pow_dim2);

        // Mask to just high_bits_count bits
        let mask = if high_bits_count >= 64 {
            u64::MAX
        } else if high_bits_count == 0 {
            0
        } else {
            (1u64 << high_bits_count) - 1
        };
        let m00_high = m00_high & mask;
        let m01_high = m01_high & mask;
        let m10_high = m10_high & mask;
        let m11_high = m11_high & mask;

        // Check if M[0][0]'s high part has a leading 1
        let m00_has_leading_one = if high_bits_count > 0 {
            (m00_high >> (high_bits_count - 1)) & 1 == 1
        } else {
            false
        };

        if m00_has_leading_one {
            count_leading_one += 1;
        }

        eprintln!(
            "KAT {:3}: n_bt={} r'={} pow_dim2={:3} high_bits={}",
            idx, n_bt, r_prime, pow_dim2, high_bits_count,
        );
        eprintln!(
            "  M[0][0] bw={:3}/{} high({} bits)=0b{:0>width$b} = 0x{:x} leading_1={}",
            m00_bw,
            e_rsp,
            high_bits_count,
            m00_high,
            m00_high,
            m00_has_leading_one,
            width = high_bits_count.max(1),
        );
        eprintln!(
            "  M[0][1] bw={:3}/{} high({} bits)=0b{:0>width$b} = 0x{:x}",
            m01_bw,
            e_rsp,
            high_bits_count,
            m01_high,
            m01_high,
            width = high_bits_count.max(1),
        );
        eprintln!(
            "  M[1][0] bw={:3}/{} high({} bits)=0b{:0>width$b} = 0x{:x}",
            m10_bw,
            e_rsp,
            high_bits_count,
            m10_high,
            m10_high,
            width = high_bits_count.max(1),
        );
        eprintln!(
            "  M[1][1] bw={:3}/{} high({} bits)=0b{:0>width$b} = 0x{:x}",
            m11_bw,
            e_rsp,
            high_bits_count,
            m11_high,
            m11_high,
            width = high_bits_count.max(1),
        );

        // If n_bt > 0, show the top n_bt bits above r'
        if n_bt > 0 {
            let top_nbt_m00 = m00_high >> r_prime;
            let top_nbt_m01 = m01_high >> r_prime;
            let nbt_mask = (1u64 << n_bt) - 1;
            eprintln!(
                "  Top {} bits above r': M00=0b{:0>width$b} M01=0b{:0>width$b}",
                n_bt,
                top_nbt_m00 & nbt_mask,
                top_nbt_m01 & nbt_mask,
                width = n_bt,
            );
        }

        // Check if r' bits (middle portion of high bits) have structure
        if r_prime > 0 {
            let rp_mask = (1u64 << r_prime) - 1;
            let m00_rp_bits = m00_high & rp_mask;
            let m01_rp_bits = m01_high & rp_mask;
            eprintln!(
                "  r' bits (low {} of high): M00=0b{:0>width$b} M01=0b{:0>width$b}",
                r_prime,
                m00_rp_bits,
                m01_rp_bits,
                width = r_prime,
            );
        }

        // Try all candidate (n_bt, r') pairs to check uniqueness
        let mut candidates = Vec::new();
        for candidate_nbt in 0u8..=3 {
            for candidate_rp in 0u8..=8 {
                let cn = candidate_nbt as usize;
                let cr = candidate_rp as usize;
                let candidate_hbc = cn + cr;
                if candidate_hbc == 0 {
                    continue;
                }
                let candidate_pd = e_rsp as i32 - cn as i32 - cr as i32;
                if candidate_pd < 2 {
                    continue;
                }

                let candidate_m00_high = extract_high(m00, candidate_pd as usize);
                let cmask = if candidate_hbc >= 64 {
                    u64::MAX
                } else {
                    (1u64 << candidate_hbc) - 1
                };
                let candidate_m00_high = candidate_m00_high & cmask;
                let has_leading_one = (candidate_m00_high >> (candidate_hbc - 1)) & 1 == 1;
                let is_correct = cn == n_bt && cr == r_prime;

                if has_leading_one {
                    candidates.push((candidate_nbt, candidate_rp, is_correct));
                    if is_correct || candidates.len() <= 5 {
                        eprintln!(
                            "  Candidate n_bt={} r'={}: M00_high=0b{:0>width$b} leading_1=true {}",
                            candidate_nbt,
                            candidate_rp,
                            candidate_m00_high,
                            if is_correct { "<<< CORRECT" } else { "" },
                            width = candidate_hbc,
                        );
                    }
                }
            }
        }

        let correct_found = candidates.iter().any(|(_, _, c)| *c);
        let is_unique = candidates.len() == 1 && correct_found;

        if is_unique {
            count_unique += 1;
        }
        if candidates.len() > 1 {
            count_ambiguous += 1;
        }
        if !correct_found {
            eprintln!("  WARNING: correct (n_bt, r') NOT among leading-1 candidates!");
        }
        eprintln!(
            "  => {} candidates, unique={}, correct_found={}",
            candidates.len(),
            is_unique,
            correct_found,
        );
        eprintln!();
    }

    eprintln!("\n=== SUMMARY ===");
    eprintln!(
        "Entries where M[0][0] high has leading 1: {}/100",
        count_leading_one
    );
    eprintln!(
        "Entries where (n_bt, r') is uniquely recoverable from leading-1 scan: {}/100",
        count_unique
    );
    eprintln!(
        "Entries where multiple candidates match: {}/100",
        count_ambiguous
    );

    eprintln!("\nbacktracking (n_bt) histogram:");
    for (bt, &count) in bt_histogram.iter().enumerate() {
        if count > 0 {
            eprintln!("  n_bt={}: {} entries", bt, count);
        }
    }

    eprintln!("\ntwo_resp_length (r') histogram:");
    for (rp, &count) in rp_histogram.iter().enumerate() {
        if count > 0 {
            eprintln!("  r'={}: {} entries", rp, count);
        }
    }

    eprintln!("\nhigh_bits_count (n_bt + r') histogram:");
    for (hbc, &count) in high_bits_width_histogram.iter().enumerate() {
        if count > 0 {
            eprintln!("  high_bits={}: {} entries", hbc, count);
        }
    }

    if count_unique == 100 {
        eprintln!("\nCONCLUSION: (n_bt, r') uniquely recoverable from M[0][0] leading-1 scan.");
        eprintln!("No separate bytes needed. 128 bytes achievable.");
    } else if count_leading_one == 100 {
        eprintln!(
            "\nCONCLUSION: M[0][0] always has leading 1, but {} entries are ambiguous.",
            count_ambiguous
        );
        eprintln!("Additional disambiguation needed (check M[0][1], det, etc.).");
    } else {
        eprintln!(
            "\nCONCLUSION: M[0][0] does NOT always have leading 1 ({}/100).",
            count_leading_one
        );
        eprintln!("Canonical form not present in current signer output.");
        eprintln!("128 bytes requires signer changes to normalize matrix.");
    }
}

#[test]
fn diagnostic_det_hint_distribution() {
    use hybrid_array::typenum::Unsigned;

    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 100);

    let e_rsp = L1::E_RSP as usize;
    let mut hint_hist = [0u32; 4];
    let mut hint_by_bt = [[0u32; 4]; 4]; // [bt][hint]
    let mut hint_by_parity = [[0u32; 4]; 2]; // [m00_odd][hint]
    let mut hint_by_trl = [[0u32; 4]; 16]; // [trl][hint]

    for (idx, entry) in entries.iter().enumerate() {
        let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();
        let compressed = sig.compress();

        let bt = sig.backtracking() as usize;
        let trl = sig.two_resp_length() as usize;
        let pow_dim2 = e_rsp - trl - bt;
        let det_precision = pow_dim2 + trl;
        let m00_odd = sig.mat()[0][0].digits()[0] & 1 != 0;

        // Read hint from the wire format: packed meta byte is right after Fp2
        let wire = compressed.to_bytes();
        let fp2_len = <L1 as SecurityLevel>::Fp2EncodedBytes::USIZE;
        let packed = wire[fp2_len];
        let hint = (packed >> 2) & 0x03;
        hint_hist[hint as usize] += 1;
        hint_by_bt[bt][hint as usize] += 1;
        hint_by_parity[m00_odd as usize][hint as usize] += 1;
        if trl < 16 {
            hint_by_trl[trl][hint as usize] += 1;
        }

        eprintln!(
            "KAT {:3}: bt={} trl={} m00_odd={} det_prec={:3} hint={}",
            idx, bt, trl, m00_odd as u8, det_precision, hint,
        );
    }

    eprintln!("\n=== DET_HINT DISTRIBUTION ===");
    for (h, &count) in hint_hist.iter().enumerate() {
        eprintln!("  hint={}: {:3}/100 ({:2}%)", h, count, count);
    }

    eprintln!("\n=== BY BACKTRACKING ===");
    for (bt, row) in hint_by_bt.iter().enumerate() {
        let total: u32 = row.iter().sum();
        if total > 0 {
            eprintln!(
                "  bt={}: hint[0]={} hint[1]={} hint[2]={} hint[3]={} (total={})",
                bt, row[0], row[1], row[2], row[3], total,
            );
        }
    }

    eprintln!("\n=== BY M[0][0] PARITY ===");
    for (p, row) in hint_by_parity.iter().enumerate() {
        let total: u32 = row.iter().sum();
        if total > 0 {
            eprintln!(
                "  m00_{}:  hint[0]={} hint[1]={} hint[2]={} hint[3]={} (total={})",
                if p == 0 { "even" } else { "odd " },
                row[0],
                row[1],
                row[2],
                row[3],
                total,
            );
        }
    }

    eprintln!("\n=== BY TWO_RESP_LENGTH ===");
    for (trl, row) in hint_by_trl.iter().enumerate() {
        let total: u32 = row.iter().sum();
        if total > 0 {
            eprintln!(
                "  trl={}: hint[0]={} hint[1]={} hint[2]={} hint[3]={} (total={})",
                trl, row[0], row[1], row[2], row[3], total,
            );
        }
    }

    let uniform_expected = 25.0f64;
    let chi2: f64 = hint_hist
        .iter()
        .map(|&c| {
            let diff = c as f64 - uniform_expected;
            diff * diff / uniform_expected
        })
        .sum();
    eprintln!(
        "\nChi-squared vs uniform: {:.2} (critical value at p=0.05, df=3: 7.81)",
        chi2,
    );
    if chi2 < 7.81 {
        eprintln!("Distribution is consistent with uniform random (cannot reject H0).");
    } else {
        eprintln!("Distribution is NOT uniform, there may be exploitable structure.");
    }

    if hint_hist.contains(&100) {
        eprintln!("\nALL entries have the same hint value! The hint is unnecessary.");
        eprintln!("128 bytes achievable without any signer changes.");
    } else if hint_hist.contains(&0) {
        eprintln!("\nSome hint values never appear, partial structure exists.");
    }
}

// --- Tate vs Weil pairing diagnostic ---

#[test]
fn diagnostic_tate_vs_weil_for_compression() {
    use hybrid_array::typenum::Unsigned;
    use sqisign_verify::ec::pairing::{fp2_dlog_2e_pub, reduced_tate, weil};
    use sqisign_verify::ec::point::{ec_dbl_iter_basis, xadd};
    use sqisign_verify::ec::EcCurve;
    use sqisign_verify::precomp::LevelPrecomp;
    use sqisign_verify::theta::HD_EXTRA_TORSION;

    use sqisign_verify::verify::{basis_from_hint, compute_challenge_curve};

    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 10);
    let nw = <L1 as SecurityLevel>::MpLimbs::USIZE;
    let e_rsp = L1::E_RSP as usize;
    let wide = 2 * nw;

    let cofactor = Level1::p_cofactor_for_2f();
    let torsion_even_power = Level1::torsion_even_power();

    let mut tate_direct_match = 0u32;
    let mut tate_swapped_match = 0u32;

    for (idx, entry) in entries.iter().enumerate() {
        let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
        let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();

        let n_bt = sig.backtracking() as usize;
        let r_prime = sig.two_resp_length() as usize;
        let pow_dim2 = e_rsp - r_prime - n_bt;
        let trl = r_prime;
        let det_precision = pow_dim2 + trl;
        let f = det_precision + HD_EXTRA_TORSION as usize;
        let g = pow_dim2 + HD_EXTRA_TORSION as usize;

        // --- Ground truth det(M) ---
        let mut prod1 = [0u64; 8];
        for i in 0..nw {
            let mut carry: u64 = 0;
            for j in 0..nw {
                if i + j >= wide {
                    break;
                }
                let p = (sig.mat()[0][0].digits()[i] as u128)
                    * (sig.mat()[1][1].digits()[j] as u128)
                    + (prod1[i + j] as u128)
                    + (carry as u128);
                prod1[i + j] = p as u64;
                carry = (p >> 64) as u64;
            }
            if i + nw < wide {
                prod1[i + nw] = carry;
            }
        }
        let mut prod2 = [0u64; 8];
        for i in 0..nw {
            let mut carry: u64 = 0;
            for j in 0..nw {
                if i + j >= wide {
                    break;
                }
                let p = (sig.mat()[0][1].digits()[i] as u128)
                    * (sig.mat()[1][0].digits()[j] as u128)
                    + (prod2[i + j] as u128)
                    + (carry as u128);
                prod2[i + j] = p as u64;
                carry = (p >> 64) as u64;
            }
            if i + nw < wide {
                prod2[i + nw] = carry;
            }
        }
        let mut det_true = [0u64; 8];
        {
            let mut borrow: u64 = 0;
            for i in 0..wide {
                let (d1, b1) = prod1[i].overflowing_sub(prod2[i]);
                let (d2, b2) = d1.overflowing_sub(borrow);
                det_true[i] = d2;
                borrow = (b1 as u64) + (b2 as u64);
            }
        }
        mask_bits(&mut det_true, det_precision);

        // --- Compute bases ---
        let mut e_chall = compute_challenge_curve::<L1>(
            sig.chall_coeff(),
            sig.backtracking(),
            pk.curve(),
            pk.hint_pk(),
        )
        .unwrap();
        let mut b_chall = basis_from_hint::<L1>(&mut e_chall, L1::F_CHR, sig.hint_chall()).unwrap();
        b_chall = ec_dbl_iter_basis(&b_chall, L1::F_CHR as usize - f, &mut e_chall);
        let ppq_chall = xadd(&b_chall.p, &b_chall.q, &b_chall.pmq);

        let mut e_aux = EcCurve::<L1>::from_a(sig.e_aux_a()).unwrap();
        let mut b_aux = basis_from_hint::<L1>(&mut e_aux, L1::F_CHR, sig.hint_aux()).unwrap();
        b_aux = ec_dbl_iter_basis(&b_aux, L1::F_CHR as usize - g, &mut e_aux);
        let ppq_aux = xadd(&b_aux.p, &b_aux.q, &b_aux.pmq);

        // --- Weil pairings (known correct) ---
        let omega_f_weil = weil::<L1>(f as u32, &b_chall.p, &b_chall.q, &ppq_chall, &mut e_chall);
        let omega_aux_weil = weil::<L1>(g as u32, &b_aux.p, &b_aux.q, &ppq_aux, &mut e_aux);

        let mut det_weil = [0u64; 8];
        let ok_weil = fp2_dlog_2e_pub::<L1>(
            &mut det_weil[..nw],
            &omega_aux_weil.inv(),
            &omega_f_weil.inv(),
            f as u32,
        )
        .is_some();
        mask_bits(&mut det_weil, det_precision);

        // --- Tate pairings ---
        // Variant 1: tate(P, Q) for both
        let tate_f = reduced_tate::<L1>(
            f as u32,
            &b_chall.p,
            &b_chall.q,
            &ppq_chall,
            &mut e_chall,
            torsion_even_power,
            cofactor,
        );
        let tate_aux = reduced_tate::<L1>(
            g as u32,
            &b_aux.p,
            &b_aux.q,
            &ppq_aux,
            &mut e_aux,
            torsion_even_power,
            cofactor,
        );

        // Variant 2: tate(Q, P) swapped args
        let tate_f_swap = reduced_tate::<L1>(
            f as u32,
            &b_chall.q,
            &b_chall.p,
            &ppq_chall,
            &mut e_chall,
            torsion_even_power,
            cofactor,
        );
        let tate_aux_swap = reduced_tate::<L1>(
            g as u32,
            &b_aux.q,
            &b_aux.p,
            &ppq_aux,
            &mut e_aux,
            torsion_even_power,
            cofactor,
        );

        // Check tate is a proper root of unity
        let mut check = tate_f.clone();
        for _ in 0..f {
            check = check.sqr();
        }
        let one = sqisign_verify::Fp2::<L1>::one();
        let tate_f_is_root = bool::from(check.ct_equal(&one));

        let mut check = tate_aux.clone();
        for _ in 0..g {
            check = check.sqr();
        }
        let tate_aux_is_root = bool::from(check.ct_equal(&one));

        eprintln!(
            "KAT {:2}: bt={} trl={} f={} g={} det_prec={}",
            idx, n_bt, r_prime, f, g, det_precision
        );
        eprintln!(
            "  Weil dlog ok={}, det matches ground truth={}",
            ok_weil,
            det_weil[..nw] == det_true[..nw]
        );
        eprintln!(
            "  Tate(P,Q)_f root_of_unity={}, Tate(P,Q)_aux root_of_unity={}",
            tate_f_is_root, tate_aux_is_root
        );

        // Try direct: dlog(tate_aux^-1, tate_f^-1, f)
        let mut det_tate_direct = [0u64; 8];
        let ok_direct = fp2_dlog_2e_pub::<L1>(
            &mut det_tate_direct[..nw],
            &tate_aux.inv(),
            &tate_f.inv(),
            f as u32,
        );
        mask_bits(&mut det_tate_direct, det_precision);
        let direct_match = ok_direct.is_some() && det_tate_direct[..nw] == det_true[..nw];
        if direct_match {
            tate_direct_match += 1;
        }
        eprintln!(
            "  Direct tate(P,Q): dlog_ok={:?} det_match={}",
            ok_direct, direct_match
        );

        // Try swapped: dlog(tate_aux_swap^-1, tate_f_swap^-1, f)
        let mut det_tate_swap = [0u64; 8];
        let ok_swap = fp2_dlog_2e_pub::<L1>(
            &mut det_tate_swap[..nw],
            &tate_aux_swap.inv(),
            &tate_f_swap.inv(),
            f as u32,
        );
        mask_bits(&mut det_tate_swap, det_precision);
        let swap_match = ok_swap.is_some() && det_tate_swap[..nw] == det_true[..nw];
        if swap_match {
            tate_swapped_match += 1;
        }
        eprintln!(
            "  Swapped tate(Q,P): dlog_ok={:?} det_match={}",
            ok_swap, swap_match
        );

        // Try: dlog(tate_aux^-1, tate_f_swap^-1, f)  (mixed)
        let mut det_mix1 = [0u64; 8];
        let ok_mix1 = fp2_dlog_2e_pub::<L1>(
            &mut det_mix1[..nw],
            &tate_aux.inv(),
            &tate_f_swap.inv(),
            f as u32,
        );
        mask_bits(&mut det_mix1, det_precision);
        let mix1_match = ok_mix1.is_some() && det_mix1[..nw] == det_true[..nw];
        eprintln!(
            "  Mixed tate_aux(P,Q) / tate_f(Q,P): dlog_ok={:?} det_match={}",
            ok_mix1, mix1_match
        );

        // Try: dlog(tate_aux_swap^-1, tate_f^-1, f)
        let mut det_mix2 = [0u64; 8];
        let ok_mix2 = fp2_dlog_2e_pub::<L1>(
            &mut det_mix2[..nw],
            &tate_aux_swap.inv(),
            &tate_f.inv(),
            f as u32,
        );
        mask_bits(&mut det_mix2, det_precision);
        let mix2_match = ok_mix2.is_some() && det_mix2[..nw] == det_true[..nw];
        eprintln!(
            "  Mixed tate_aux(Q,P) / tate_f(P,Q): dlog_ok={:?} det_match={}",
            ok_mix2, mix2_match
        );

        // Try without inversions
        let mut det_noinv1 = [0u64; 8];
        let ok_noinv1 = fp2_dlog_2e_pub::<L1>(&mut det_noinv1[..nw], &tate_aux, &tate_f, f as u32);
        mask_bits(&mut det_noinv1, det_precision);
        let noinv1_match = ok_noinv1.is_some() && det_noinv1[..nw] == det_true[..nw];
        eprintln!(
            "  No-inv tate(P,Q): dlog_ok={:?} det_match={}",
            ok_noinv1, noinv1_match
        );

        let mut det_noinv2 = [0u64; 8];
        let ok_noinv2 = fp2_dlog_2e_pub::<L1>(
            &mut det_noinv2[..nw],
            &tate_aux_swap,
            &tate_f_swap,
            f as u32,
        );
        mask_bits(&mut det_noinv2, det_precision);
        let noinv2_match = ok_noinv2.is_some() && det_noinv2[..nw] == det_true[..nw];
        eprintln!(
            "  No-inv tate(Q,P): dlog_ok={:?} det_match={}",
            ok_noinv2, noinv2_match
        );

        // Also try: tate_f * tate_f_swap should equal weil_f (if weil = tate(P,Q)/tate(Q,P))
        let product_f = tate_f.mul(&tate_f_swap.inv());
        let product_matches_weil = bool::from(product_f.re.ct_equal(&omega_f_weil.re))
            && bool::from(product_f.im.ct_equal(&omega_f_weil.im));
        let product_f_neg = tate_f_swap.mul(&tate_f.inv());
        let neg_matches_weil = bool::from(product_f_neg.re.ct_equal(&omega_f_weil.re))
            && bool::from(product_f_neg.im.ct_equal(&omega_f_weil.im));
        eprintln!(
            "  tate(P,Q)/tate(Q,P) == weil? {}  tate(Q,P)/tate(P,Q) == weil? {}",
            product_matches_weil, neg_matches_weil
        );

        eprintln!();
    }

    eprintln!("=== SUMMARY (10 entries) ===");
    eprintln!("Tate direct (same args as Weil): {}/10", tate_direct_match);
    eprintln!(
        "Tate swapped (Q,P instead of P,Q): {}/10",
        tate_swapped_match
    );
}

// --- Compressed format tests ---

#[test]
fn test_compress_and_verify_kat_entry_0() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    let entry = &entries[0];

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();
    let msg = &entry.sm[L1_SIG_BYTES..];

    let compressed = sig.compress();
    pk.verify(msg, &compressed)
        .expect("compressed verification failed on KAT entry 0");
}

#[test]
fn test_compress_and_verify_all_kat_entries() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 5);

    for (i, entry) in entries.iter().enumerate() {
        let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
        let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();
        let msg = &entry.sm[L1_SIG_BYTES..];

        let compressed = sig.compress();
        pk.verify(msg, &compressed)
            .unwrap_or_else(|e| panic!("KAT {}: compressed verify failed: {:?}", i, e));
    }
}

#[test]
fn test_compressed_serialization_roundtrip_kat() {
    use sqisign_verify::formats::CompressedSignature;

    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    let entry = &entries[0];

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();
    let msg = &entry.sm[L1_SIG_BYTES..];

    let compressed = sig.compress();
    let wire = compressed.to_bytes();
    let decoded =
        CompressedSignature::<L1>::from_bytes(&wire[..CompressedSignature::<L1>::WIRE_BYTES])
            .expect("compressed from_bytes failed");

    pk.verify(msg, &decoded)
        .expect("verify after roundtrip failed");
}

#[test]
fn test_compressed_rejects_wrong_message() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    let entry = &entries[0];

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();

    let compressed = sig.compress();
    assert!(
        pk.verify(b"wrong message", &compressed).is_err(),
        "compressed verify should reject wrong message"
    );
}

#[test]
fn test_compressed_any_signature_dispatch() {
    use sqisign_verify::formats::{AnySignature, SignatureFormat};

    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 1);
    let entry = &entries[0];

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();
    let msg = &entry.sm[L1_SIG_BYTES..];

    let compressed = sig.compress();
    let any = AnySignature::Compressed(compressed);
    assert_eq!(any.format(), SignatureFormat::Compressed);
    pk.verify(msg, &any)
        .expect("AnySignature::verify failed for compressed");
}

#[test]
fn test_compress_decompress_byte_identical_all_100() {
    let content = KAT_LVL1;
    let entries = parse_kat_entries(content, 100);
    assert_eq!(entries.len(), 100);

    for (i, entry) in entries.iter().enumerate() {
        let pk = PublicKey::<L1>::from_bytes(&entry.pk)
            .unwrap_or_else(|_| panic!("KAT {}: failed to decode public key", i));
        let sig = Signature::<L1>::from_bytes(&entry.sm[..L1_SIG_BYTES]).unwrap();
        let compressed = sig.compress();
        let recovered = compressed
            .decompress(&pk)
            .unwrap_or_else(|e| panic!("KAT {}: decompress failed: {:?}", i, e));

        let orig_bytes = sig.to_bytes();
        let rec_bytes = recovered.to_bytes();
        assert_eq!(
            orig_bytes.as_slice(),
            rec_bytes.as_slice(),
            "KAT {}: compress/decompress roundtrip not byte-identical (bt={}, trl={})",
            i,
            sig.backtracking(),
            sig.two_resp_length(),
        );
    }
}
