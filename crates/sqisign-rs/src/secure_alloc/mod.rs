//! A global-allocator wrapper that zeros all heap memory on deallocation, so
//! intermediate `BigInt` values from the signing path do not persist in freed
//! pages. Measured overhead: <1% on signing.
//!
//! With the default `std` feature the crate registers this as the process
//! [`#[global_allocator]`](https://doc.rust-lang.org/std/alloc/trait.GlobalAlloc.html)
//! automatically (secure by default, no caller action). Build with
//! `default-features = false` for `no_std`: the crate then registers no global
//! allocator (the embedder supplies one), and a `std` binary that disabled the
//! feature can still opt in with the
//! [`enable_secure_allocator!`](crate::enable_secure_allocator) macro.
//!
//! Because a library setting `#[global_allocator]` is a process-wide singleton,
//! an application that registers its own allocator (jemalloc, mimalloc, ...) or
//! depends on another crate that does will get a "multiple global allocators"
//! error. Disable the `std` feature in that case and manage the allocator
//! yourself (optionally wrapping it in [`ZeroizingAllocator`]).

use core::alloc::{GlobalAlloc, Layout};
use zeroize::Zeroize;

/// A global allocator wrapper that zeros memory before deallocation.
///
/// Generic over the inner allocator `A`. With the default `std` feature the
/// crate registers `ZeroizingAllocator<System>` as the global allocator
/// automatically (see the module docs); a `no_std` embedder wraps its own
/// allocator instead, either directly or via
/// [`enable_secure_allocator!`](crate::enable_secure_allocator).
pub struct ZeroizingAllocator<A: GlobalAlloc>(pub A);

/// Secure-by-default: under the `std` feature, scrub all freed heap memory by
/// making the zeroing wrapper around the system allocator the process global
/// allocator. No effect on `no_std` builds (the feature is off there).
#[cfg(feature = "std")]
#[global_allocator]
static SQISIGN_GLOBAL_ALLOC: ZeroizingAllocator<std::alloc::System> =
    ZeroizingAllocator(std::alloc::System);

/// Register the zeroing allocator around the system allocator in a downstream
/// `std` binary.
///
/// Only needed when the default `std` feature is **disabled** — with it enabled
/// the crate already registers the allocator, and invoking this as well is a
/// second `#[global_allocator]` (a compile error). Place at the top of
/// `main.rs`:
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
