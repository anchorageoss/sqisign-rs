//! The harness's timer, behind one trait. On x86-64 a measurement is the
//! time-stamp counter (`rdtsc`) and the monotonic wall clock around the
//! call, as in every session so far. On every other architecture (aarch64,
//! Apple Silicon) it is the monotonic wall clock alone, in nanoseconds:
//! there is no cycle column, and none is made up. The wall-clock timer can
//! be chosen on x86-64 too (`--timer wall`), which is how the path a Mac
//! takes is exercised on the x86-64 machine.
//!
//! Reading Apple's kperf cycle counters needs a private framework and root;
//! it is not implemented here.

use std::time::Instant;

/// One measurement: nanoseconds always, cycles where a cycle counter exists.
#[derive(Clone, Copy, Debug)]
pub struct Sample {
    pub cycles: Option<u64>,
    pub ns: u128,
}

impl Sample {
    /// The CSV cell of the cycle column: the count, or `NA`.
    pub fn cycles_cell(&self) -> String {
        self.cycles.map_or_else(|| "NA".to_string(), |c| c.to_string())
    }
}

/// A timer: a started measurement that can be stopped.
pub trait Timer: Sized {
    /// What the timer reads, for file headers.
    const NAME: &'static str;
    fn start() -> Self;
    fn stop(self) -> Sample;
    /// A monotonic counter in the timer's own unit (cycles, or nanoseconds
    /// since the first call), for the verifier's phase profiler.
    fn counter() -> u64;
}

/// `rdtsc` and the wall clock (x86-64).
#[cfg(target_arch = "x86_64")]
pub struct TscTimer {
    t0: Instant,
    c0: u64,
}

#[cfg(target_arch = "x86_64")]
fn rdtsc() -> u64 {
    // SAFETY: rdtsc is available on every x86-64 CPU.
    unsafe { core::arch::x86_64::_rdtsc() }
}

#[cfg(target_arch = "x86_64")]
impl Timer for TscTimer {
    const NAME: &'static str = "rdtsc and the monotonic wall clock";
    #[inline(always)]
    fn start() -> Self {
        let t0 = Instant::now();
        TscTimer { t0, c0: rdtsc() }
    }
    #[inline(always)]
    fn stop(self) -> Sample {
        let c1 = rdtsc();
        Sample {
            cycles: Some(c1 - self.c0),
            ns: self.t0.elapsed().as_nanos(),
        }
    }
    #[inline(always)]
    fn counter() -> u64 {
        rdtsc()
    }
}

/// The monotonic wall clock alone (every architecture).
pub struct WallTimer {
    t0: Instant,
}

impl Timer for WallTimer {
    const NAME: &'static str = "the monotonic wall clock (std::time::Instant), nanoseconds";
    #[inline(always)]
    fn start() -> Self {
        WallTimer { t0: Instant::now() }
    }
    #[inline(always)]
    fn stop(self) -> Sample {
        Sample {
            cycles: None,
            ns: self.t0.elapsed().as_nanos(),
        }
    }
    fn counter() -> u64 {
        use std::sync::OnceLock;
        static EPOCH: OnceLock<Instant> = OnceLock::new();
        EPOCH.get_or_init(Instant::now).elapsed().as_nanos() as u64
    }
}

/// Keep the calling thread on a performance core where the system offers a
/// way that needs no privilege. macOS: the user-interactive quality-of-
/// service class (there is no affinity call). Elsewhere: nothing; the
/// session script pins with `taskset`. Returns a line for the log.
pub fn prefer_performance_core() -> &'static str {
    #[cfg(target_os = "macos")]
    {
        extern "C" {
            fn pthread_set_qos_class_self_np(qos_class: u32, relative_priority: i32) -> i32;
        }
        const QOS_CLASS_USER_INTERACTIVE: u32 = 0x21;
        // SAFETY: a libSystem call on the current thread with a valid class.
        let rc = unsafe { pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0) };
        if rc == 0 {
            "thread set to the user-interactive QoS class (performance cores). Run with nothing else open and on mains power."
        } else {
            "could not set the user-interactive QoS class. Run with nothing else open and on mains power."
        }
    }
    #[cfg(not(target_os = "macos"))]
    {
        "no QoS class on this system; pin the process with taskset (paper_session.sh does)."
    }
}
