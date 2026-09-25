//! A criterion `Measurement` in `rdtsc` cycles (constant-rate TSC), so the
//! statistics criterion computes are in the unit the C references report.

use criterion::measurement::{Measurement, ValueFormatter};
use criterion::Throughput;

pub struct Cycles;

fn rdtsc() -> u64 {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        core::arch::x86_64::_rdtsc()
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        0
    }
}

impl Measurement for Cycles {
    type Intermediate = u64;
    type Value = u64;
    fn start(&self) -> u64 {
        rdtsc()
    }
    fn end(&self, i: u64) -> u64 {
        rdtsc() - i
    }
    fn add(&self, v1: &u64, v2: &u64) -> u64 {
        v1 + v2
    }
    fn zero(&self) -> u64 {
        0
    }
    fn to_f64(&self, v: &u64) -> f64 {
        *v as f64
    }
    fn formatter(&self) -> &dyn ValueFormatter {
        &CyclesFormatter
    }
}

pub struct CyclesFormatter;

impl ValueFormatter for CyclesFormatter {
    fn scale_values(&self, typical: f64, values: &mut [f64]) -> &'static str {
        let (factor, unit) = if typical >= 1e9 {
            (1e-9, "Gcycles")
        } else if typical >= 1e6 {
            (1e-6, "Mcycles")
        } else if typical >= 1e3 {
            (1e-3, "kcycles")
        } else {
            (1.0, "cycles")
        };
        for v in values.iter_mut() {
            *v *= factor;
        }
        unit
    }
    fn scale_throughputs(
        &self,
        _typical: f64,
        _t: &Throughput,
        _values: &mut [f64],
    ) -> &'static str {
        "elem/cycle"
    }
    fn scale_for_machines(&self, _values: &mut [f64]) -> &'static str {
        "cycles"
    }
}
