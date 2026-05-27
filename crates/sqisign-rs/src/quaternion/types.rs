//!
//! Provides vectors, matrices, quaternion elements, lattices, and ideals.
//! All types default to zero-initialized values.

use super::intbig::Ibz;
use num_traits::{One, Zero};
use std::ops::{Index, IndexMut};
use zeroize::Zeroize;

/// 2-element vector of big integers.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IbzVec2(pub [Ibz; 2]);

impl Default for IbzVec2 {
    fn default() -> Self {
        Self([Ibz::zero(), Ibz::zero()])
    }
}

impl Index<usize> for IbzVec2 {
    type Output = Ibz;
    fn index(&self, i: usize) -> &Ibz {
        &self.0[i]
    }
}

impl IndexMut<usize> for IbzVec2 {
    fn index_mut(&mut self, i: usize) -> &mut Ibz {
        &mut self.0[i]
    }
}

/// 4-element vector of big integers.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IbzVec4(pub [Ibz; 4]);

impl Default for IbzVec4 {
    fn default() -> Self {
        Self([Ibz::zero(), Ibz::zero(), Ibz::zero(), Ibz::zero()])
    }
}

impl Index<usize> for IbzVec4 {
    type Output = Ibz;
    fn index(&self, i: usize) -> &Ibz {
        &self.0[i]
    }
}

impl IndexMut<usize> for IbzVec4 {
    fn index_mut(&mut self, i: usize) -> &mut Ibz {
        &mut self.0[i]
    }
}

/// 2×2 matrix of big integers. Indexed as `mat[row][col]`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IbzMat2x2(pub [[Ibz; 2]; 2]);

impl Default for IbzMat2x2 {
    fn default() -> Self {
        Self([[Ibz::zero(), Ibz::zero()], [Ibz::zero(), Ibz::zero()]])
    }
}

impl Index<usize> for IbzMat2x2 {
    type Output = [Ibz; 2];
    fn index(&self, i: usize) -> &[Ibz; 2] {
        &self.0[i]
    }
}

impl IndexMut<usize> for IbzMat2x2 {
    fn index_mut(&mut self, i: usize) -> &mut [Ibz; 2] {
        &mut self.0[i]
    }
}

/// 4×4 matrix of big integers. Indexed as `mat[row][col]`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IbzMat4x4(pub [[Ibz; 4]; 4]);

impl Default for IbzMat4x4 {
    fn default() -> Self {
        Self([
            [Ibz::zero(), Ibz::zero(), Ibz::zero(), Ibz::zero()],
            [Ibz::zero(), Ibz::zero(), Ibz::zero(), Ibz::zero()],
            [Ibz::zero(), Ibz::zero(), Ibz::zero(), Ibz::zero()],
            [Ibz::zero(), Ibz::zero(), Ibz::zero(), Ibz::zero()],
        ])
    }
}

impl Index<usize> for IbzMat4x4 {
    type Output = [Ibz; 4];
    fn index(&self, i: usize) -> &[Ibz; 4] {
        &self.0[i]
    }
}

impl IndexMut<usize> for IbzMat4x4 {
    fn index_mut(&mut self, i: usize) -> &mut [Ibz; 4] {
        &mut self.0[i]
    }
}

/// Quaternion algebra `B_{p,∞}` defined by the prime `p`.
///
/// The algebra has basis `{1, i, j, ij}` where `i^2 = -1`, `j^2 = -p`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuatAlg {
    pub p: Ibz,
}

impl QuatAlg {
    pub fn new(p: &Ibz) -> Self {
        Self { p: p.clone() }
    }
}

impl Default for QuatAlg {
    fn default() -> Self {
        Self { p: Ibz::zero() }
    }
}

/// Element of a quaternion algebra, represented as `coord / denom`
/// where `coord` is a 4-vector of numerators in basis `{1, i, j, ij}`
/// and `denom` is a common denominator.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuatAlgElem {
    pub coord: IbzVec4,
    pub denom: Ibz,
}

impl Default for QuatAlgElem {
    fn default() -> Self {
        Self {
            coord: IbzVec4::default(),
            denom: Ibz::one(),
        }
    }
}

/// Lattice in the quaternion algebra, represented by an integer basis
/// matrix divided by a common denominator.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuatLattice {
    pub basis: IbzMat4x4,
    pub denom: Ibz,
}

impl Default for QuatLattice {
    fn default() -> Self {
        Self {
            basis: IbzMat4x4::default(),
            denom: Ibz::one(),
        }
    }
}

/// Left ideal of a maximal order, represented as a lattice with its
/// norm and a copy of its parent order lattice.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuatLeftIdeal {
    pub lattice: QuatLattice,
    pub norm: Ibz,
    pub parent_order: QuatLattice,
}

impl Default for QuatLeftIdeal {
    fn default() -> Self {
        Self {
            lattice: QuatLattice::default(),
            norm: Ibz::zero(),
            parent_order: QuatLattice::default(),
        }
    }
}

/// p-extremal maximal order, storing the order lattice and two
/// distinguished elements `z` and `t`, plus the level `q` of the order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuatPExtremalMaximalOrder {
    pub order: QuatLattice,
    pub z: QuatAlgElem,
    pub t: QuatAlgElem,
    pub q: i64,
}

impl Default for QuatPExtremalMaximalOrder {
    fn default() -> Self {
        Self {
            order: QuatLattice::default(),
            z: QuatAlgElem::default(),
            t: QuatAlgElem::default(),
            q: 1,
        }
    }
}

/// Parameters for the `represent_integer` algorithm.
pub struct QuatRepresentIntegerParams<'a> {
    pub algebra: &'a QuatAlg,
    pub order: &'a QuatPExtremalMaximalOrder,
    pub primality_test_iterations: u32,
}

// Zeroize impls
//
// SECURITY: BigInt does not expose its backing Vec<u64>, so ibz_zeroize
// replaces the value with zero but cannot scrub the freed heap allocation.
// Use a zeroing allocator for comprehensive heap hygiene.

impl Zeroize for IbzVec2 {
    fn zeroize(&mut self) {
        for v in self.0.iter_mut() {
            super::intbig::ibz_zeroize(v);
        }
    }
}

impl Zeroize for IbzVec4 {
    fn zeroize(&mut self) {
        for v in self.0.iter_mut() {
            super::intbig::ibz_zeroize(v);
        }
    }
}

impl Zeroize for IbzMat2x2 {
    fn zeroize(&mut self) {
        for row in self.0.iter_mut() {
            for v in row.iter_mut() {
                super::intbig::ibz_zeroize(v);
            }
        }
    }
}

impl Zeroize for IbzMat4x4 {
    fn zeroize(&mut self) {
        for row in self.0.iter_mut() {
            for v in row.iter_mut() {
                super::intbig::ibz_zeroize(v);
            }
        }
    }
}

impl Zeroize for QuatAlgElem {
    fn zeroize(&mut self) {
        self.coord.zeroize();
        super::intbig::ibz_zeroize(&mut self.denom);
    }
}

impl Zeroize for QuatLattice {
    fn zeroize(&mut self) {
        self.basis.zeroize();
        super::intbig::ibz_zeroize(&mut self.denom);
    }
}

impl Zeroize for QuatLeftIdeal {
    fn zeroize(&mut self) {
        self.lattice.zeroize();
        super::intbig::ibz_zeroize(&mut self.norm);
        self.parent_order.zeroize();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    #[test]
    fn test_quat_alg() {
        let p = BigInt::from(7);
        let alg = QuatAlg::new(&p);
        assert_eq!(alg.p, p);
    }

    #[test]
    fn test_quat_alg_elem() {
        let mut elem = QuatAlgElem::default();
        for i in 0..4 {
            elem.coord[i] = BigInt::from(i as i32);
        }
        elem.denom = BigInt::from(1);
        assert!(elem.denom.is_one());
        for i in 0..4 {
            assert_eq!(elem.coord[i], BigInt::from(i as i32));
        }
    }

    #[test]
    fn test_ibz_vec_2() {
        let mut vec = IbzVec2::default();
        for i in 0..2 {
            vec[i] = BigInt::from(i as i32);
        }
        for i in 0..2 {
            assert_eq!(vec[i], BigInt::from(i as i32));
        }
    }

    #[test]
    fn test_ibz_vec_4() {
        let mut vec = IbzVec4::default();
        for i in 0..4 {
            vec[i] = BigInt::from(i as i32);
        }
        for i in 0..4 {
            assert_eq!(vec[i], BigInt::from(i as i32));
        }
    }

    #[test]
    fn test_ibz_mat_2x2() {
        let mut mat = IbzMat2x2::default();
        for i in 0..2 {
            for j in 0..2 {
                mat[i][j] = BigInt::from((i + j) as i32);
            }
        }
        for i in 0..2 {
            for j in 0..2 {
                assert_eq!(mat[i][j], BigInt::from((i + j) as i32));
            }
        }
    }

    #[test]
    fn test_ibz_mat_4x4() {
        let mut mat = IbzMat4x4::default();
        for i in 0..4 {
            for j in 0..4 {
                mat[i][j] = BigInt::from((i + j) as i32);
            }
        }
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(mat[i][j], BigInt::from((i + j) as i32));
            }
        }
    }

    #[test]
    fn test_quat_lattice() {
        let mut lat = QuatLattice::default();
        for i in 0..4 {
            for j in 0..4 {
                lat.basis[i][j] = BigInt::from((i + j) as i32);
            }
        }
        lat.denom = BigInt::from(1);
        assert!(lat.denom.is_one());
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(lat.basis[i][j], BigInt::from((i + j) as i32));
            }
        }
    }

    #[test]
    fn test_quat_left_ideal() {
        let mut lideal = QuatLeftIdeal::default();
        for i in 0..4 {
            for j in 0..4 {
                lideal.lattice.basis[i][j] = BigInt::from((i + j) as i32);
            }
        }
        lideal.lattice.denom = BigInt::from(1);
        lideal.norm = BigInt::from(5);
        assert!(lideal.parent_order.denom.is_one());
        assert!(lideal.lattice.denom.is_one());
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(lideal.lattice.basis[i][j], BigInt::from((i + j) as i32));
            }
        }
        assert_eq!(lideal.norm, BigInt::from(5));
    }
}
