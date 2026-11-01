//! Phase timing of the compact verifier for the benchmark session
//! (SECENG-2424). The verifier's stages take a [`Profiler`]; the shipped
//! entry points pass [`NoProfiler`], whose methods are empty and compile
//! away, and [`hd_verify_profiled`](crate::hd::self_contained::hd_verify_profiled)
//! passes a [`CycleProfiler`] fed by the caller's cycle counter. The phase
//! boundaries are the functions named on each slot.

/// Slots of a profile: cycles accumulated per phase, both half-chains added.
pub const SLOT_BASIS: usize = 0;
pub const SLOT_SETUP: usize = 1;
pub const SLOT_DIM2: usize = 2;
pub const SLOT_GLUING: usize = 3;
pub const SLOT_PLAIN_PUSHES: usize = 4;
pub const SLOT_PLAIN_CODOMAIN: usize = 5;
pub const SLOT_COMPARE: usize = 6;
pub const SLOT_IMAGE: usize = 7;
/// The number of slots.
pub const SLOTS: usize = 8;

/// The names of the slots, in order: what each one encloses.
pub const SLOT_NAMES: [&str; SLOTS] = [
    "basis generation: recover_challenge, recover_response, hd_torsion_basis",
    "Cornacchia and kernel setup: norm_equation_2f_minus_q, make_canonical, product null, starting_two_symplectic_matrices, gluing base changes, b_kpp",
    "dimension-2 steps: IsogenyChainDim2::new, both halves",
    "gluing: KaniGluingChainHalf::new after the dimension-2 chain, and the gluing images of the generators, both halves",
    "plain dimension-4 steps, point pushes: double_iter and image of the saved bases, both halves",
    "plain dimension-4 steps, codomain: IsogenyDim4::from_kernel and the next structure's precompute, both halves",
    "final comparisons: middle_codomain_matches",
    "image check: hd_image_check",
];

/// The finer slots inside the two setup phases (SECENG-2427, SECENG-2428):
/// the challenge recovery, the response recovery and the kernel setup,
/// in the order of the computation. Working evidence; the record's slots
/// above are unchanged by them.
pub const SUB_SLOTS: usize = 16;

/// The names of the sub-slots, in order: what each one encloses.
pub const SUB_NAMES: [&str; SUB_SLOTS] = [
    "chal: hd_torsion_basis on E_pk",
    "chal: x-only doublings of the basis and ec_ladder3pt",
    "chal: the challenge isogeny chain with its image points",
    "chal: the lifts of the image points",
    "chal: the Weil pairing(s), the lambda squarings and the sign of phi Q",
    "chal: P_resc and Q_resc",
    "resp: hd_torsion_basis on E_com",
    "resp: the (f - r) doublings to R_com and S_com",
    "resp: the Weil pairing at 2^r, the inverse and fp2_dlog_2e",
    "resp: recover_response_cd and the images phi R, phi S",
    "setup: the second hd_torsion_basis on E_com (P_com)",
    "setup: norm_equation_2f_minus_q (Cornacchia)",
    "setup: the 4-torsion doublings, make_canonical x2, weil4, the nulls, points_m",
    "matrices: starting_two_symplectic_matrices and the two gluing base changes",
    "b_kpp, both halves: lambda, lambda phi S and point_matrix_product_k",
    "ThetaStructureDim4::new of the start, both halves",
];

/// A sink for phase timings.
pub trait Profiler {
    /// The current count of the caller's cycle counter; `0` when unused.
    #[inline(always)]
    fn now(&self) -> u64 {
        0
    }
    /// Add `cycles` to `slot`.
    #[inline(always)]
    fn add(&mut self, _slot: usize, _cycles: u64) {}
    /// Add `cycles` to the sub-slot `sub`.
    #[inline(always)]
    fn sub_add(&mut self, _sub: usize, _cycles: u64) {}
}

/// The profiler of the shipped entry points: does nothing.
#[derive(Clone, Copy, Debug, Default)]
pub struct NoProfiler;

impl Profiler for NoProfiler {}

/// A profiler that reads the caller's cycle counter and accumulates per slot.
pub struct CycleProfiler<'a> {
    clock: &'a dyn Fn() -> u64,
    /// Cycles per slot.
    pub slots: [u64; SLOTS],
    /// Cycles per sub-slot.
    pub subs: [u64; SUB_SLOTS],
}

impl<'a> CycleProfiler<'a> {
    /// A profiler over `clock`, every slot at zero.
    pub fn new(clock: &'a dyn Fn() -> u64) -> Self {
        Self {
            clock,
            slots: [0; SLOTS],
            subs: [0; SUB_SLOTS],
        }
    }
}

impl Profiler for CycleProfiler<'_> {
    #[inline(always)]
    fn now(&self) -> u64 {
        (self.clock)()
    }
    #[inline(always)]
    fn add(&mut self, slot: usize, cycles: u64) {
        self.slots[slot] += cycles;
    }
    #[inline(always)]
    fn sub_add(&mut self, sub: usize, cycles: u64) {
        self.subs[sub] += cycles;
    }
}
