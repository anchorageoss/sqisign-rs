#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use sqisign_rs::id2iso::sign_precomp::SigningPrecomp;
use sqisign_rs::Level1;
use sqisign_rs::keygen::keygen::protocols_keygen;
use sqisign_rs::sign::sign::protocols_sign;
use std::time::Instant;

fn main() {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(997)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build()
        .expect("failed to build profiler");

    eprintln!("Initializing precomputed constants...");
    let t0 = Instant::now();
    let precomp = SigningPrecomp::<Level1>::level1();
    eprintln!("  precomp init: {:?}", t0.elapsed());

    let mut rng = rand::thread_rng();

    eprintln!("Generating keypair...");
    let t0 = Instant::now();
    let (pk, sk) = protocols_keygen::<Level1>(&mut rng, &precomp);
    eprintln!("  keygen: {:?}", t0.elapsed());

    let sk_bytes = sk.to_bytes().expect("sk encoding must succeed");
    let mut sk = sqisign_rs::SecretKey::<Level1>::from_bytes(&sk_bytes)
        .expect("sk round-trip must succeed");
    sk.populate_from_pk(&pk);

    let num_signs = 5;
    eprintln!("Signing {} messages...", num_signs);
    let t_total = Instant::now();

    for i in 0..num_signs {
        let msg = format!("profile message {}", i);
        let t0 = Instant::now();
        let sig = protocols_sign::<Level1>(&pk, &sk, msg.as_bytes(), &mut rng)
            .expect("signing must succeed");
        let sign_time = t0.elapsed();

        let t0 = Instant::now();
        let result = sig.verify(&pk, msg.as_bytes());
        let verify_time = t0.elapsed();

        let valid = result.is_ok();
        eprintln!(
            "  msg {}: sign {:?}, verify {:?}, valid={}",
            i, sign_time, verify_time, valid
        );
        assert!(valid);
    }
    eprintln!("Total signing+verify: {:?}", t_total.elapsed());

    match guard.report().build() {
        Ok(report) => {
            let file = std::fs::File::create("flamegraph.svg").unwrap();
            report.flamegraph(file).unwrap();
            eprintln!("Flamegraph written to flamegraph.svg");
        }
        Err(e) => {
            eprintln!("Failed to generate flamegraph: {:?}", e);
        }
    }
}
