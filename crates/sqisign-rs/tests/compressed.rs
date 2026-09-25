//! The compressed signature format (three matrix entries, the fourth
//! recovered from two Weil pairings): the round trip at each level, the
//! decoder's rejections, tampering, and the acceptance measurement on
//! generated signatures (`SQISIGN_SURVEY_SIGS` per level, default 3: every
//! honest signature must compress, verify compressed, and decompress to a
//! signature the standard verifier accepts).

mod common;

use common::{all_levels, rng};
use signature::{SignatureEncoding, Verifier};
use sqisign_rs::{generate, CompressedSignature, Error, Level1, Level3, Level5, SigningLevel};
use sqisign_verify::compressed::compressed_bytes;
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
    roundtrip::<Level1>(176);
}

#[test]
fn roundtrip_levels_3_and_5() {
    if !all_levels() {
        eprintln!("skipped in a debug build (SQISIGN_ALL_LEVELS=1 to force)");
        return;
    }
    roundtrip::<Level3>(269);
    roundtrip::<Level5>(353);
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
    let fp = params.fp_encoded_bytes;
    let eb = (params.response_bits as usize).div_ceil(8);
    let decode = |b: &[u8]| CompressedSignature::<Level1>::from_bytes(b).map(|_| ());
    // length
    assert_eq!(decode(&cb[..n - 1]), Err(Error::InvalidLength));
    let mut long = cb.to_vec();
    long.push(0);
    assert_eq!(decode(&long), Err(Error::InvalidLength));
    // an entry at or above 2^RESPONSE_BITS (196 bits in 25 bytes: 4 slack bits)
    for e in 0..3 {
        let mut t = cb;
        t[2 * fp + (e + 1) * eb - 1] |= 0x80;
        assert_eq!(decode(&t), Err(Error::MalformedInput), "entry {e}");
    }
    // both pivots even
    {
        let mut t = cb;
        t[2 * fp] &= !1;
        t[2 * fp + eb] &= !1;
        assert_eq!(decode(&t), Err(Error::MalformedInput));
    }
    // the bits byte outside its four bits
    {
        let mut t = cb;
        t[n - 3] |= 0x10;
        assert_eq!(decode(&t), Err(Error::MalformedInput));
    }
    // a non-canonical field element
    {
        let mut t = cb;
        for b in t[..fp].iter_mut() {
            *b = 0xff;
        }
        assert_eq!(decode(&t), Err(Error::MalformedInput));
    }
    // every single-bit change outside the bits byte is rejected by the
    // decoder or the verifier
    for i in (0..n).filter(|i| *i != n - 3) {
        for bit in 0..8 {
            let mut t = cb;
            t[i] ^= 1 << bit;
            let rejected = match CompressedSignature::<Level1>::from_bytes(&t) {
                Err(_) => true,
                Ok(c) => pk.verify_compressed(msg, &c).is_err(),
            };
            assert!(rejected, "byte {i} bit {bit} flipped and accepted");
        }
    }
    // the bits byte: the compressed verifier accepts a value exactly when
    // the standard verifier accepts the signature it stands for (the
    // standard format's own freedom in those bits, nothing more)
    let mut accepted = 0;
    for v in 0..16u8 {
        let mut t = cb;
        t[n - 3] = v;
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
    let mut bits_hist = [0usize; 16];
    for i in 0..count {
        let msg = (i as u64).to_le_bytes();
        let sig = sk.sign(&msg, &mut r).expect("sign");
        let cs = sig.compress();
        bits_hist[cs.inner().bits as usize] += 1;
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
         signatures the standard verifier accepts; the four transmitted bits' histogram {bits_hist:?}",
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
