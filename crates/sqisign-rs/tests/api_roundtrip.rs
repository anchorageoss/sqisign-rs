//! The typed API end to end at level I: generate, sign, verify; the byte
//! codecs of key, signature and signing key; the RustCrypto traits; and
//! the rejections a verifier must make (wrong message, tampered bytes, the
//! negated basis-change matrix of ePrint 2026/1305).

use signature::{RandomizedSigner, SignatureEncoding};
use sqisign_rs::mp::ShakeRng;
use sqisign_rs::sqisign::Params;
use sqisign_rs::{generate, Error, Level1, PublicKey, Signature, SigningKey, Verifier};
use sqisign_verify::sqisign::{level1, signature_from_bytes, signature_to_bytes, verify};
use sqisign_verify::{Level, MAX_SIGNATURE_BYTES, P324_3};

/// A deterministic `rand_core` generator over the crate's SHAKE stream.
struct Det(ShakeRng);
impl rand_core::RngCore for Det {
    fn next_u32(&mut self) -> u32 {
        let mut b = [0u8; 4];
        self.fill_bytes(&mut b);
        u32::from_le_bytes(b)
    }
    fn next_u64(&mut self) -> u64 {
        let mut b = [0u8; 8];
        self.fill_bytes(&mut b);
        u64::from_le_bytes(b)
    }
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        use sqisign_verify::rng::Rng;
        assert!(self.0.fill(dest));
    }
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.fill_bytes(dest);
        Ok(())
    }
}
impl rand_core::CryptoRng for Det {}

fn rng(label: &str) -> Det {
    Det(ShakeRng::new(label.as_bytes(), b"tst"))
}

#[test]
fn generate_sign_verify_and_codecs() {
    let mut r = rng("api");
    let (pk, sk): (PublicKey<Level1>, SigningKey<Level1>) = generate(&mut r);
    let msg = b"round three";
    let sig = sk.sign(msg, &mut r).expect("sign");
    assert!(pk.verify(msg, &sig).is_ok());
    assert_eq!(pk.verify(b"round two", &sig), Err(Error::InvalidSignature));

    // public key and signature codecs
    let pkb = pk.to_bytes();
    assert_eq!(pkb.len(), 83);
    let pk2 = PublicKey::<Level1>::from_bytes(&pkb).expect("pk decode");
    assert_eq!(pk2.to_bytes(), pkb);
    let sigb = sig.to_bytes();
    assert_eq!(sigb.len(), 200);
    let sig2 = Signature::<Level1>::from_bytes(&sigb).expect("sig decode");
    assert_eq!(sig2.to_bytes(), sigb);
    assert!(pk2.verify(msg, &sig2).is_ok());
    assert_eq!(
        PublicKey::<Level1>::from_bytes(&pkb[..82]).map(|_| ()),
        Err(Error::InvalidLength)
    );
    assert_eq!(
        Signature::<Level1>::from_bytes(&sigb[..199]).map(|_| ()),
        Err(Error::InvalidLength)
    );

    // the RustCrypto traits
    let sig3 = sk.try_sign_with_rng(&mut r, msg).expect("RandomizedSigner");
    assert!(Verifier::verify(&pk, msg, &sig3).is_ok());
    let repr: <Signature<Level1> as SignatureEncoding>::Repr = sig3.clone().to_bytes();
    assert_eq!(repr.len(), 200);

    // the signing key's codec
    let skb = sk.to_bytes();
    assert_eq!(skb.len(), Level1::SK_BYTES);
    let sk2 = SigningKey::<Level1>::from_bytes(&skb).expect("sk decode");
    assert_eq!(sk2.public_key().to_bytes(), pkb);
    let sig4 = sk2.sign(msg, &mut r).expect("sign with decoded key");
    assert!(pk.verify(msg, &sig4).is_ok());
    assert_eq!(
        SigningKey::<Level1>::from_bytes(&skb[..skb.len() - 1]).map(|_| ()),
        Err(Error::InvalidLength)
    );

    // a prepared key verifies the same
    let prepared = pk.prepare().expect("prepare");
    assert!(sqisign_verify::verify_bytes_prepared(&prepared, msg, &sigb).is_ok());
    assert!(sqisign_verify::verify_bytes_prepared(&prepared, b"other", &sigb).is_err());

    // every single-byte flip of the signature is rejected or decodes to a
    // signature that fails (a sample of positions, the tail included)
    for i in [0usize, 41, 82, 90, 120, 150, 182, 198, 199] {
        let mut t = sigb.to_vec();
        t[i] ^= 0x01;
        match Signature::<Level1>::from_bytes(&t) {
            Ok(s) => assert!(pk.verify(msg, &s).is_err(), "flip at {i} accepted"),
            Err(e) => assert_eq!(e, Error::MalformedInput),
        }
    }
}

/// The negated basis-change matrix (the second representative of the same
/// isogeny, ePrint 2026/1305) is not in canonical form and the verifier
/// rejects it, as the reference does.
#[test]
fn negated_matrix_is_rejected() {
    let params: Params<{ sqisign_rs::precomp::p324_3::IBZ_NLIMBS }> = sqisign_rs::sqisign::level1();
    let vp = &params.verify;
    let mut r = rng("canonical");
    let (pk, sk): (PublicKey<Level1>, SigningKey<Level1>) = generate(&mut r);
    let msg = b"canonical form";
    let sig = sk.sign(msg, &mut r).expect("sign");
    let mut inner = sig.inner().clone();
    // negate every entry modulo 2^(RESPONSE_BITS + 2): the same isogeny
    let bits = vp.response_bits + 2;
    for row in inner.mat_bchall_can_to_bchall.iter_mut() {
        for x in row.iter_mut() {
            let mut borrow = 0u64;
            let mut neg = [0u64; sqisign_verify::ec::MAX_ORDER_WORDS];
            for (k, (n, xi)) in neg.iter_mut().zip(x.iter()).enumerate() {
                let (d, b1) = 0u64.overflowing_sub(*xi);
                let (d, b2) = d.overflowing_sub(borrow);
                *n = d;
                borrow = (b1 | b2) as u64;
                let _ = k;
            }
            // reduce modulo 2^bits
            let top = (bits / 64) as usize;
            for (k, n) in neg.iter_mut().enumerate() {
                if k > top {
                    *n = 0;
                } else if k == top && bits % 64 != 0 {
                    *n &= (1u64 << (bits % 64)) - 1;
                }
            }
            *x = neg;
        }
    }
    let mut buf = [0u8; MAX_SIGNATURE_BYTES];
    let n = signature_to_bytes(vp, &inner, &mut buf).expect("encode");
    let negated = signature_from_bytes::<P324_3>(vp, &buf[..n]).expect("decode");
    assert!(
        !verify(vp, pk.inner(), &negated, msg),
        "the negated matrix verified"
    );
    assert!(verify(vp, pk.inner(), sig.inner(), msg));
    let _ = level1();
}
