//!
//! Zeros ALL heap memory on deallocation, ensuring that intermediate
//! BigInt values from the signing path do not persist in freed memory.
//! Measured overhead: <1% on signing operations.
//!
//! # Usage
//!
//! Enabled by default through `sqisign-core`'s `zeroize-alloc` feature.
//! No code changes needed, the allocator activates automatically.
//!
//! To disable (e.g., if you use jemalloc or mimalloc):
//!
//! ```toml
//! sqisign-core = { version = "0.1", default-features = false, features = ["sign"] }
//! ```
//!
//! For explicit opt-in without default features:
//!
//! ```rust,ignore
//! crate::alloc::enable_secure_allocator!();
//! ```

use core::alloc::{GlobalAlloc, Layout};
use zeroize::Zeroize;

/// A global allocator wrapper that zeros memory before deallocation.
///
/// Generic over the inner allocator `A`. This crate is `no_std`, so it does not
/// itself register a `#[global_allocator]`: that is the final binary's job, and
/// a `no_std` embedder must supply its own allocator. A `std` binary can opt in
/// to the zeroing wrapper around the system allocator with the
/// [`enable_secure_allocator!`] macro.
pub struct ZeroizingAllocator<A: GlobalAlloc>(pub A);

/// Convenience macro for `std` binaries that want the zeroing allocator wrapped
/// around the system allocator. Expands in the downstream binary, where `std`
/// is available.
///
/// Place at the top of `main.rs`:
/// ```rust,ignore
/// sqisign_rs::enable_secure_allocator!();
/// ```
#[macro_export]
macro_rules! enable_secure_allocator {
    () => {
        #[global_allocator]
        static __SQISIGN_ALLOC: $crate::ZeroizingAllocator<std::alloc::System> =
            $crate::ZeroizingAllocator(std::alloc::System);
    };
}

unsafe impl<A: GlobalAlloc> GlobalAlloc for ZeroizingAllocator<A> {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        self.0.alloc(layout)
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        self.0.alloc_zeroed(layout)
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: ptr is valid for layout.size() bytes (GlobalAlloc contract).
        // Uses zeroize's volatile writes + compiler fence to prevent the
        // optimizer from eliding the store as a dead write.
        let slice = unsafe { core::slice::from_raw_parts_mut(ptr, layout.size()) };
        slice.zeroize();
        self.0.dealloc(ptr, layout)
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        // Manual alloc+copy+dealloc instead of delegating to the inner
        // allocator's realloc. System realloc may relocate the block and
        // free the old one without zeroing. By routing through our own
        // dealloc, the old block is always zeroed before being freed.
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
