use criterion::{criterion_group, criterion_main, Criterion};
use sqisign_rs::id2iso::sign_precomp::SigningPrecomp;
use sqisign_rs::keygen::keygen::protocols_keygen;
use sqisign_rs::sign::sign::protocols_sign;
use sqisign_rs::Level1;
use sqisign_rs::Verifier;

type L = Level1;

fn make_precomp() -> SigningPrecomp<L> {
    SigningPrecomp::<L>::level1()
}

fn bench_keygen(c: &mut Criterion) {
    let precomp = make_precomp();
    let mut rng = rand::thread_rng();

    c.bench_function("keygen (Level 1)", |b| {
        b.iter(|| {
            let (pk, sk) = protocols_keygen::<L>(&mut rng, &precomp);
            std::hint::black_box((&pk, &sk));
        })
    });
}

fn bench_sign(c: &mut Criterion) {
    let precomp = make_precomp();
    let mut rng = rand::thread_rng();

    let (pk, sk) = protocols_keygen::<L>(&mut rng, &precomp);
    let msg = b"benchmark message for SQIsign signing";

    c.bench_function("sign (Level 1)", |b| {
        b.iter(|| {
            let sig = protocols_sign(&pk, &sk, msg, &mut rng).expect("signing must succeed");
            std::hint::black_box(&sig);
        })
    });
}

fn bench_keygen_sign_verify(c: &mut Criterion) {
    let precomp = make_precomp();
    let mut rng = rand::thread_rng();
    let msg = b"end-to-end benchmark";

    c.bench_function("keygen+sign+verify (Level 1)", |b| {
        b.iter(|| {
            let (pk, sk) = protocols_keygen::<L>(&mut rng, &precomp);
            let sig = protocols_sign(&pk, &sk, msg, &mut rng).expect("signing must succeed");
            assert!(pk.verify(msg, &sig).is_ok());
        })
    });
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = bench_keygen, bench_sign, bench_keygen_sign_verify,
}
criterion_main!(benches);
