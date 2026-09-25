//! Run-time choice of the x86-64 field kernels (P26, item 0): the inline
//! assembly of the `<prime>_asm` modules needs BMI2 (`mulx`) and ADX
//! (`adcx`/`adox`); CPUID leaf 7 says whether the CPU has them, read once
//! through `core::arch` (this crate is `no_std`) and cached in an atomic.
//! Every kernel branches on the cached answer, one relaxed load. The
//! portable kernels share the representation, so nothing else changes.
//!
//! Tests and CI force the portable path with `force_portable` (the test
//! harnesses read `PRISM_FORCE_PORTABLE=1`); `force_auto` re-detects.

use core::sync::atomic::{AtomicU8, Ordering};

const AUTO: u8 = 0;
const ASM: u8 = 1;
const PORTABLE: u8 = 2;

static MODE: AtomicU8 = AtomicU8::new(AUTO);

/// Whether the assembly kernels may run: BMI2 and ADX present and not
/// forced off.
#[inline(always)]
pub(crate) fn asm_available() -> bool {
    match MODE.load(Ordering::Relaxed) {
        ASM => true,
        PORTABLE => false,
        _ => detect(),
    }
}

// `__cpuid` / `__cpuid_count` are `unsafe fn` on the MSRV (1.80) and safe
// on newer toolchains, so the block below is required on one and flagged
// as unnecessary on the other; the allow keeps both warning-free.
#[allow(unused_unsafe)]
#[cold]
fn detect() -> bool {
    // SAFETY: `cpuid` is available on every x86-64 CPU; leaf 7 is queried
    // only when the maximum basic leaf allows it.
    let ok = unsafe {
        let max = core::arch::x86_64::__cpuid(0).eax;
        if max >= 7 {
            let r = core::arch::x86_64::__cpuid_count(7, 0);
            (r.ebx >> 8) & 1 == 1 && (r.ebx >> 19) & 1 == 1
        } else {
            false
        }
    };
    MODE.store(if ok { ASM } else { PORTABLE }, Ordering::Relaxed);
    ok
}

/// Use the portable kernels from now on, whatever the CPU has.
pub fn force_portable() {
    MODE.store(PORTABLE, Ordering::Relaxed);
}

/// Detect again on the next operation (undoes `force_portable`).
pub fn force_auto() {
    MODE.store(AUTO, Ordering::Relaxed);
}

/// The kernels in use, for benchmark headers and logs.
pub fn backend() -> &'static str {
    if asm_available() {
        "x86-64 assembly (mulx, adcx, adox)"
    } else {
        "portable (saturated limbs, u128)"
    }
}
