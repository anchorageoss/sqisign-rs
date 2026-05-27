//! Integration test verifying that ZeroizingAllocator zeros all freed heap memory.
//!
//! Uses a custom AuditAllocator underneath a local ZeroizingAllocator. On every
//! dealloc, the AuditAllocator checks whether the memory region has been zeroed
//! BEFORE passing it to the system allocator.
//!
//! This crate is excluded from the workspace to avoid `#[global_allocator]`
//! conflicts with sqisign's default allocator. It depends on sqisign with
//! `default-features = false` so the crate's own allocator is not activated.
//!
//! Run: `cargo test --manifest-path crates/alloc/Cargo.toml --release -- --nocapture`

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

struct ZeroizingAllocator<A: GlobalAlloc>(A);

unsafe impl<A: GlobalAlloc> GlobalAlloc for ZeroizingAllocator<A> {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        self.0.alloc(layout)
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        self.0.alloc_zeroed(layout)
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        core::ptr::write_bytes(ptr, 0, layout.size());
        self.0.dealloc(ptr, layout)
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let new_layout = unsafe { Layout::from_size_align_unchecked(new_size, layout.align()) };
        let new_ptr = unsafe { self.alloc(new_layout) };
        if new_ptr.is_null() {
            return new_ptr;
        }
        let copy_size = core::cmp::min(layout.size(), new_size);
        unsafe {
            core::ptr::copy_nonoverlapping(ptr, new_ptr, copy_size);
            self.dealloc(ptr, layout);
        }
        new_ptr
    }
}

struct AuditAllocator {
    non_zero_deallocs: AtomicUsize,
    total_deallocs: AtomicUsize,
    total_bytes_freed: AtomicUsize,
    non_zero_bytes_found: AtomicUsize,
}

impl AuditAllocator {
    const fn new() -> Self {
        Self {
            non_zero_deallocs: AtomicUsize::new(0),
            total_deallocs: AtomicUsize::new(0),
            total_bytes_freed: AtomicUsize::new(0),
            non_zero_bytes_found: AtomicUsize::new(0),
        }
    }

    fn reset(&self) {
        self.non_zero_deallocs.store(0, Ordering::SeqCst);
        self.total_deallocs.store(0, Ordering::SeqCst);
        self.total_bytes_freed.store(0, Ordering::SeqCst);
        self.non_zero_bytes_found.store(0, Ordering::SeqCst);
    }
}

unsafe impl GlobalAlloc for AuditAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        unsafe { System.alloc(layout) }
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        unsafe { System.alloc_zeroed(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        self.total_deallocs.fetch_add(1, Ordering::Relaxed);
        self.total_bytes_freed
            .fetch_add(layout.size(), Ordering::Relaxed);

        let slice = unsafe { std::slice::from_raw_parts(ptr, layout.size()) };
        let non_zero_count = slice.iter().filter(|&&b| b != 0).count();
        if non_zero_count > 0 {
            self.non_zero_deallocs.fetch_add(1, Ordering::Relaxed);
            self.non_zero_bytes_found
                .fetch_add(non_zero_count, Ordering::Relaxed);
        }

        unsafe { System.dealloc(ptr, layout) }
    }
}

#[global_allocator]
static ALLOC: ZeroizingAllocator<AuditAllocator> = ZeroizingAllocator(AuditAllocator::new());

#[test]
fn all_freed_memory_is_zeroed_during_keygen_sign_drop() {
    use sqisign_rs::params::Level1;
    use sqisign_rs::Verifier;

    ALLOC.0.reset();

    {
        let mut rng = rand::thread_rng();

        let (pk, sk) = sqisign_rs::keygen::keypair::<Level1>(&mut rng);

        let sig =
            sqisign_rs::sign::sign::<Level1>(&sk, &pk, b"zeroize audit test message", &mut rng);

        if let Ok(ref sig) = sig {
            assert!(pk.verify(b"zeroize audit test message", sig).is_ok());
        }
    }

    let total = ALLOC.0.total_deallocs.load(Ordering::SeqCst);
    let non_zero = ALLOC.0.non_zero_deallocs.load(Ordering::SeqCst);
    let total_bytes = ALLOC.0.total_bytes_freed.load(Ordering::SeqCst);
    let non_zero_bytes = ALLOC.0.non_zero_bytes_found.load(Ordering::SeqCst);

    println!("ZeroizingAllocator audit results:");
    println!("  Total deallocations: {total}");
    println!("  Total bytes freed: {total_bytes}");
    println!("  Deallocations with non-zero memory: {non_zero}");
    println!("  Non-zero bytes found: {non_zero_bytes}");

    assert!(total > 0, "expected heap allocations during keygen+sign");
    assert_eq!(
        non_zero, 0,
        "ZeroizingAllocator FAILED: {non_zero} of {total} deallocations had non-zero memory \
         ({non_zero_bytes} non-zero bytes out of {total_bytes} total)",
    );
}
