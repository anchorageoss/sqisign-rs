//! The heap figures of the paper (Section "Results", last paragraph): the
//! peak of the heap during one verification from bytes, measured with a
//! counting allocator around the call, for the compact verifier with four
//! points and with two at the three levels, and for the dimension-2
//! verifier.
//!
//!     paper-heap --out <dir>
//!
//! Writes `<dir>/heap.csv`: `row,signature,peak_bytes,allocations`. The peak
//! is the high-water mark of live bytes above the level at the start of the
//! call.

use std::alloc::{GlobalAlloc, Layout, System};
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

use sqisign_rs::compact::{
    compact_keygen, compact_public_key_bytes, compact_public_key_to_bytes, compact_sign,
    compact_signature_bytes,
};
use sqisign_rs::mp::ShakeRng;
use sqisign_rs::sqisign::{level1, level3, level5, Params};
use sqisign_verify::hd::{hd_verify_bytes_with, HdLevel, HdVerifyOptions, KernelPoints};
use sqisign_verify::{P324_3, P500_27, P664_17};

struct Counting;

static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);
static COUNT: AtomicUsize = AtomicUsize::new(0);

// SAFETY: every call is forwarded to the system allocator unchanged; the
// counters are updated around it.
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let p = System.alloc(l);
        if !p.is_null() {
            let live = LIVE.fetch_add(l.size(), Relaxed) + l.size();
            PEAK.fetch_max(live, Relaxed);
            COUNT.fetch_add(1, Relaxed);
        }
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        System.dealloc(p, l);
        LIVE.fetch_sub(l.size(), Relaxed);
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, new: usize) -> *mut u8 {
        let q = System.realloc(p, l, new);
        if !q.is_null() {
            if new >= l.size() {
                let live = LIVE.fetch_add(new - l.size(), Relaxed) + (new - l.size());
                PEAK.fetch_max(live, Relaxed);
            } else {
                LIVE.fetch_sub(l.size() - new, Relaxed);
            }
            COUNT.fetch_add(1, Relaxed);
        }
        q
    }
}

#[global_allocator]
static A: Counting = Counting;

/// The peak of live heap bytes above the starting level during `f`, and the
/// number of allocations.
fn peak<T>(f: impl FnOnce() -> T) -> (T, usize, usize) {
    let base = LIVE.load(Relaxed);
    PEAK.store(base, Relaxed);
    let c0 = COUNT.load(Relaxed);
    let r = f();
    (
        r,
        PEAK.load(Relaxed) - base,
        COUNT.load(Relaxed) - c0,
    )
}

fn compact<L: HdLevel, const N: usize>(level: &str, params: &Params<N>, out: &mut fs::File) {
    let mut r = ShakeRng::new(format!("paper heap {level}").as_bytes(), b"hea");
    let (pk, sk) = compact_keygen::<L, N>(params, &mut r).unwrap();
    let mut pkb = vec![0u8; compact_public_key_bytes::<L>()];
    compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
    for j in 0..5 {
        let msg = format!("heap {j}").into_bytes();
        let mut sig = vec![0u8; compact_signature_bytes::<L>()];
        compact_sign::<L, N>(params, &pk, &sk, &msg, &mut r, &mut sig).unwrap();
        for (name, mode) in [
            ("four", KernelPoints::Four),
            ("two", KernelPoints::TwoWithImageCheck),
        ] {
            let opts = HdVerifyOptions {
                kernel_points: mode,
            };
            let ((verdict, _), bytes, n) =
                peak(|| hd_verify_bytes_with::<L>(&sig, &pkb, &msg, &opts));
            assert_eq!(verdict, Ok(()));
            writeln!(out, "dim4-{level}-{name},{j},{bytes},{n}").unwrap();
        }
    }
}

fn dim2<L: sqisign_rs::SigningLevel>(level: &str, out: &mut fs::File) {
    use rand::SeedableRng;
    use sqisign_rs::{generate, PublicKey, Signature};
    let mut rng = rand::rngs::StdRng::seed_from_u64(7);
    let (pk, sk) = generate::<L>(&mut rng);
    let sig = sk.sign(b"heap", &mut rng).unwrap();
    let (pkb, sigb) = (pk.to_bytes(), sig.to_bytes());
    let (ok, bytes, n) = peak(|| {
        let pk = PublicKey::<L>::from_bytes(&pkb).unwrap();
        let sig = Signature::<L>::from_bytes(&sigb).unwrap();
        pk.verify(b"heap", &sig).is_ok()
    });
    assert!(ok);
    writeln!(out, "dim2-{level}-verify,0,{bytes},{n}").unwrap();
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let out = match argv.as_slice() {
        [flag, dir] if flag == "--out" => PathBuf::from(dir),
        _ => panic!("usage: paper-heap --out <dir>"),
    };
    fs::create_dir_all(&out).unwrap();
    let mut f = fs::File::create(out.join("heap.csv")).unwrap();
    writeln!(f, "# peak of live heap bytes above the level at the start of one verification from bytes, counting allocator around the call").unwrap();
    writeln!(f, "row,signature,peak_bytes,allocations").unwrap();
    const N1: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;
    const N3: usize = sqisign_rs::precomp::p500_27::IBZ_NLIMBS;
    const N5: usize = sqisign_rs::precomp::p664_17::IBZ_NLIMBS;
    compact::<P324_3, N1>("I", &level1(), &mut f);
    compact::<P500_27, N3>("III", &level3(), &mut f);
    compact::<P664_17, N5>("V", &level5(), &mut f);
    dim2::<sqisign_rs::Level1>("I", &mut f);
    dim2::<sqisign_rs::Level3>("III", &mut f);
    dim2::<sqisign_rs::Level5>("V", &mut f);
}
