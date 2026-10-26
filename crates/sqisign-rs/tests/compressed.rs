//! The compressed signature format (three matrix entries, the fourth
//! recovered from two Weil pairings, the fields bit-packed): the round trip
//! at each level, the decoder's rejections, tampering, and the acceptance measurement on
//! generated signatures (`SQISIGN_SURVEY_SIGS` per level, default 3: every
//! honest signature must compress, verify compressed, and decompress to a
//! signature the standard verifier accepts).

mod common;

use common::{all_levels, rng};
use signature::{SignatureEncoding, Verifier};
use sqisign_rs::{generate, CompressedSignature, Error, Level1, Level3, Level5, SigningLevel};
use sqisign_verify::compressed::{compressed_bits, compressed_bytes, entry_bits};
use sqisign_verify::sqisign::verify;
use sqisign_verify::Level;

fn roundtrip<L: SigningLevel>(expect_len: usize) {
    let mut r = rng(L::NAME, b"cmp");
    let (pk, sk) = generate::<L>(&mut r);
    let msg = b"compressed";
    let sig = sk.sign(msg, &mut r).expect("sign");
    let c = sig.compress();
    let cb = c.to_bytes();
    assert_eq!(cb.len(), expect_len, "{}", L::NAME);
    assert_eq!(compressed_bytes(&L::PARAMS), expect_len);
    let c2 = CompressedSignature::<L>::from_bytes(&cb).expect("decode");
    assert_eq!(c2.to_bytes(), cb);
    assert!(pk.verify_compressed(msg, &c2).is_ok(), "{}", L::NAME);
    assert!(Verifier::verify(&pk, msg, &c2).is_ok());
    assert_eq!(
        pk.verify_compressed(b"other", &c2),
        Err(Error::InvalidSignature)
    );
    // the recovered standard signature verifies and compresses back to the
    // same bytes (it equals the signer's up to the free top bits)
    let back = c2.decompress(&pk).expect("decompress");
    assert!(pk.verify(msg, &back).is_ok());
    assert_eq!(back.compress().to_bytes(), cb);
    // another key does not recover it
    let (pk_other, _) = generate::<L>(&mut r);
    assert!(pk_other.verify_compressed(msg, &c2).is_err());
    let repr: <CompressedSignature<L> as SignatureEncoding>::Repr = c2.clone().to_bytes();
    assert_eq!(repr.len(), expect_len);
    assert!(format!("{c2:?}").starts_with("CompressedSignature<"));
}

#[test]
fn roundtrip_level1() {
    roundtrip::<Level1>(174);
}

#[test]
fn roundtrip_levels_3_and_5() {
    if !all_levels() {
        eprintln!("skipped in a debug build (SQISIGN_ALL_LEVELS=1 to force)");
        return;
    }
    roundtrip::<Level3>(266);
    roundtrip::<Level5>(352);
}

#[test]
fn decoder_and_tampering() {
    let mut r = rng("strict", b"cmp");
    let (pk, sk) = generate::<Level1>(&mut r);
    let msg = b"strict";
    let sig = sk.sign(msg, &mut r).expect("sign");
    let cb = sig.compress().to_bytes();
    let n = cb.len();
    let params = &Level1::PARAMS;
    // the bit layout: A_aux (two coordinates), m00, m01, m_var, chl, hints
    let fb = params.fp_bits as usize;
    let eb = entry_bits(params);
    let rb = params.response_bits as usize;
    let (off_m00, off_m01, off_var) = (2 * fb, 2 * fb + eb, 2 * fb + 2 * eb);
    let total = compressed_bits(params);
    assert!(total < 8 * n, "level I has padding bits");
    let flip = |b: &mut [u8], i: usize| b[i / 8] ^= 1 << (i % 8);
    let clear = |b: &mut [u8], i: usize| b[i / 8] &= !(1 << (i % 8));
    let set = |b: &mut [u8], i: usize| b[i / 8] |= 1 << (i % 8);
    let decode = |b: &[u8]| CompressedSignature::<Level1>::from_bytes(b).map(|_| ());
    // length
    assert_eq!(decode(&cb[..n - 1]), Err(Error::InvalidLength));
    let mut long = cb.to_vec();
    long.push(0);
    assert_eq!(decode(&long), Err(Error::InvalidLength));
    // a padding bit set
    {
        let mut t = cb;
        set(&mut t, total);
        assert_eq!(decode(&t), Err(Error::MalformedInput));
    }
    // both pivots even
    {
        let mut t = cb;
        clear(&mut t, off_m00);
        clear(&mut t, off_m01);
        assert_eq!(decode(&t), Err(Error::MalformedInput));
    }
    // a non-canonical field element: the real part all ones, above p
    {
        let mut t = cb;
        for i in 0..fb {
            set(&mut t, i);
        }
        assert_eq!(decode(&t), Err(Error::MalformedInput));
    }
    // every single-bit change outside bit RESPONSE_BITS of the three
    // entries is rejected by the decoder or the verifier
    let free = [off_m00 + rb, off_m01 + rb, off_var + rb];
    for i in (0..8 * n).filter(|i| !free.contains(i)) {
        let mut t = cb;
        flip(&mut t, i);
        let rejected = match CompressedSignature::<Level1>::from_bytes(&t) {
            Err(_) => true,
            Ok(c) => pk.verify_compressed(msg, &c).is_err(),
        };
        assert!(rejected, "bit {i} flipped and accepted");
    }
    // bit RESPONSE_BITS of the entries: the compressed verifier accepts a
    // value exactly when the standard verifier accepts the signature it
    // stands for (the standard format's own freedom in those bits, nothing
    // more)
    let mut accepted = 0;
    for v in 0..8u8 {
        let mut t = cb;
        for (k, &i) in free.iter().enumerate() {
            clear(&mut t, i);
            if (v >> k) & 1 == 1 {
                set(&mut t, i);
            }
        }
        let c = CompressedSignature::<Level1>::from_bytes(&t).unwrap();
        let lift = c.decompress(&pk).unwrap();
        let standard = pk.verify(msg, &lift).is_ok();
        assert_eq!(
            pk.verify_compressed(msg, &c).is_ok(),
            standard,
            "bits {v:#x}"
        );
        accepted += standard as usize;
    }
    assert!(
        (1..=2).contains(&accepted),
        "bits values accepted: {accepted}"
    );
}

/// Acceptance on generated signatures.
fn survey<L: SigningLevel>(count: usize) {
    let params = &L::PARAMS;
    let mut r = rng(&format!("survey {}", L::NAME), b"cmp");
    let (pk, sk) = generate::<L>(&mut r);
    let mut bits_hist = [0usize; 8];
    for i in 0..count {
        let msg = (i as u64).to_le_bytes();
        let sig = sk.sign(&msg, &mut r).expect("sign");
        let cs = sig.compress();
        let r = params.response_bits;
        let top = |x: &[u64]| ((x[(r / 64) as usize] >> (r % 64)) & 1) as usize;
        let c = cs.inner();
        bits_hist[top(&c.m00) | top(&c.m01) << 1 | top(&c.m_var) << 2] += 1;
        assert!(
            pk.verify_compressed(&msg, &cs).is_ok(),
            "{}: signature {i}",
            L::NAME
        );
        let lift = cs.decompress(&pk).unwrap();
        assert!(
            verify(params, pk.inner(), lift.inner(), &msg),
            "{}: signature {i}: the recovered signature fails the standard verifier",
            L::NAME
        );
    }
    eprintln!(
        "{}: {count} signatures compressed, verified compressed, decompressed to standard \
         signatures the standard verifier accepts; histogram of bit RESPONSE_BITS of (m00, m01, m_var) {bits_hist:?}",
        L::NAME
    );
}

fn survey_count() -> usize {
    std::env::var("SQISIGN_SURVEY_SIGS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(3)
}

#[test]
fn survey_level1() {
    survey::<Level1>(survey_count());
}

#[test]
fn survey_levels_3_and_5() {
    if !all_levels() {
        eprintln!("skipped in a debug build (SQISIGN_ALL_LEVELS=1 to force)");
        return;
    }
    survey::<Level3>(survey_count());
    survey::<Level5>(survey_count());
}
