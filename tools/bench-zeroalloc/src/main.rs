use std::alloc::{GlobalAlloc, Layout, System};

struct ZeroizingAllocator;

unsafe impl GlobalAlloc for ZeroizingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        System.alloc(layout)
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        System.alloc_zeroed(layout)
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        core::ptr::write_bytes(ptr, 0, layout.size());
        System.dealloc(ptr, layout)
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        if new_size < layout.size() {
            core::ptr::write_bytes(ptr.add(new_size), 0, layout.size() - new_size);
        }
        System.realloc(ptr, layout, new_size)
    }
}

#[cfg(feature = "zero-alloc")]
#[global_allocator]
static ALLOC: ZeroizingAllocator = ZeroizingAllocator;

use sqisign_rs::id2iso::sign_precomp::SigningPrecomp;
use sqisign_rs::{Level1, Verifier};
use sqisign_rs::keygen::keygen::protocols_keygen;
use sqisign_rs::sign::sign::protocols_sign;
use std::time::Instant;

fn main() {
    #[cfg(feature = "zero-alloc")]
    eprintln!("Mode: ZeroizingAllocator ENABLED");
    #[cfg(not(feature = "zero-alloc"))]
    eprintln!("Mode: System allocator (baseline)");

    let precomp = SigningPrecomp::<Level1>::level1();
    let mut rng = rand::thread_rng();

    let (pk, sk) = protocols_keygen::<Level1>(&mut rng, &precomp);
    let sk_bytes = sk.to_bytes().expect("sk encoding must succeed");
    let mut sk = sqisign_rs::SecretKey::<Level1>::from_bytes(&sk_bytes)
        .expect("sk round-trip must succeed");
    sk.populate_from_pk(&pk);

    let num_signs = 10;
    let mut sign_times = Vec::with_capacity(num_signs);
    let mut verify_times = Vec::with_capacity(num_signs);

    for i in 0..num_signs {
        let msg = format!("bench message {}", i);

        let t0 = Instant::now();
        let sig = protocols_sign::<Level1>(&pk, &sk, msg.as_bytes(), &mut rng)
            .expect("signing must succeed");
        sign_times.push(t0.elapsed());

        let t0 = Instant::now();
        let result = pk.verify(msg.as_bytes(), &sig);
        verify_times.push(t0.elapsed());

        assert!(result.is_ok());
    }

    let sign_total: std::time::Duration = sign_times.iter().sum();
    let verify_total: std::time::Duration = verify_times.iter().sum();
    let sign_avg = sign_total / num_signs as u32;
    let verify_avg = verify_total / num_signs as u32;

    eprintln!("Results ({} iterations):", num_signs);
    eprintln!("  Sign avg:   {:?}", sign_avg);
    eprintln!("  Verify avg: {:?}", verify_avg);
    eprintln!("  Sign total: {:?}", sign_total);

    println!("{}", sign_avg.as_millis());
}
