//! The antipode of a graded connected Hopf algebra.
//!
//! The antipode $\chi$ is the map that makes the dual of a left module into a left module again:
//! since $A$ is not commutative, transposing the action of $A$ on $M$ gives a module over
//! $A^{\mathrm{op}}$, and $\chi \colon A \to A^{\mathrm{op}}$ is what identifies the two. It is
//! therefore what [`crate::module::FDModule::dual`] is built on.

use std::sync::Arc;

use fp::vector::{FpSlice, FpSliceMut, FpVector};
use once::OnceVec;

use crate::algebra::{Algebra, Bialgebra};

/// The antipode of a graded connected Hopf algebra, computed on demand and cached by degree.
///
/// $\chi$ is determined by the Hopf algebra axioms: it is the identity in degree 0, and for $x$ of
/// positive degree
///
/// $$ \sum \chi(x') x'' = \varepsilon(x) = 0, $$
///
/// the sum being over the coproduct $\Delta x = \sum x' \otimes x''$. Splitting off the single term
/// whose right factor is $1$ leaves an expression for $\chi(x)$ in terms of $\chi$ in strictly lower
/// degrees, which is the recursion used here. Since [`Bialgebra::coproduct`] is only defined on the
/// elements [`Bialgebra::decompose`] produces — the individual $Sq^i$ in the Adem basis — the
/// recursion is applied to those, and $\chi$ is extended to the rest of the basis by
/// anti-multiplicativity:
///
/// $$ \chi(xy) = (-1)^{|x||y|} \chi(y)\chi(x). $$
pub struct Antipode<A: Algebra + Bialgebra> {
    algebra: Arc<A>,
    /// `chi[d][i]` is $\chi$ of the `i`th basis element of degree `d`, as an element of degree `d`.
    chi: OnceVec<Vec<FpVector>>,
}

impl<A: Algebra + Bialgebra> Antipode<A> {
    pub fn new(algebra: Arc<A>) -> Self {
        Self {
            algebra,
            chi: OnceVec::new(),
        }
    }

    pub fn algebra(&self) -> Arc<A> {
        Arc::clone(&self.algebra)
    }

    /// Compute $\chi$ on every basis element up to and including `degree`.
    ///
    /// Calling this is optional: the accessors below do it themselves. It is public because
    /// computing a range up front is cheaper than one degree at a time.
    pub fn compute_through_degree(&self, degree: i32) {
        if degree < 0 {
            return;
        }
        self.algebra.compute_basis(degree);
        // Each degree needs only strictly lower degrees, so a single increasing pass suffices.
        self.chi
            .extend(degree as usize, |d| self.compute_degree(d as i32));
    }

    /// $\chi$ of a basis element.
    pub fn on_basis_element(&self, degree: i32, idx: usize) -> &FpVector {
        assert!(degree >= 0, "the antipode is only defined in degrees >= 0");
        self.compute_through_degree(degree);
        &self.chi[degree as usize][idx]
    }

    /// Add `coeff` times $\chi$ of `element` to `result`.
    pub fn apply(&self, mut result: FpSliceMut, coeff: u32, degree: i32, element: FpSlice) {
        let p = self.algebra.prime();
        for (idx, value) in element.iter_nonzero() {
            result.add(
                self.on_basis_element(degree, idx).as_slice(),
                (coeff * value) % p,
            );
        }
    }

    fn compute_degree(&self, degree: i32) -> Vec<FpVector> {
        (0..self.algebra.dimension(degree))
            .map(|idx| self.compute_basis_element(degree, idx))
            .collect()
    }

    fn compute_basis_element(&self, degree: i32, idx: usize) -> FpVector {
        let p = self.algebra.prime();
        let mut result = FpVector::new(p, self.algebra.dimension(degree));

        if degree == 0 {
            // The algebra is connected, so degree 0 is spanned by the unit and `chi` fixes it.
            result.set_entry(idx, 1);
            return result;
        }

        let factors = self.algebra.decompose(degree, idx);
        if factors.len() < 2 {
            // An element that does not factor further is one `coproduct` accepts, so the Hopf
            // recursion applies to it directly. Every term but the one with a degree-0 right factor
            // involves a strictly lower degree on the left, which is already computed.
            for (left_degree, left_idx, right_degree, right_idx) in
                self.algebra.coproduct(degree, idx)
            {
                if right_degree == 0 {
                    // This is the `chi(x) * 1` term, i.e. the one being solved for.
                    continue;
                }
                debug_assert!(left_degree < degree);
                self.algebra.multiply_element_by_basis_element(
                    result.as_slice_mut(),
                    p - 1,
                    left_degree,
                    self.chi[left_degree as usize][left_idx].as_slice(),
                    right_degree,
                    right_idx,
                );
            }
            return result;
        }

        // `decompose` lists the factors in the order they act on a module element, so the product is
        // the reverse: `x = a_k ... a_1`. Anti-multiplicativity then gives
        // `chi(x) = (-1)^s chi(a_1) ... chi(a_k)` with `s` the sum of `|a_i||a_j|` over pairs `i < j`,
        // which is symmetric in the factors.
        let mut sign_exponent = 0;
        for (i, (first, _)) in factors.iter().enumerate() {
            for (second, _) in &factors[i + 1..] {
                sign_exponent += first * second;
            }
        }

        let mut product = FpVector::new(p, self.algebra.dimension(0));
        product.set_entry(0, 1);
        let mut product_degree = 0;
        for &(factor_degree, factor_idx) in &factors {
            debug_assert!(factor_degree < degree);
            let next_degree = product_degree + factor_degree;
            let mut next = FpVector::new(p, self.algebra.dimension(next_degree));
            self.algebra.multiply_element_by_element(
                next.as_slice_mut(),
                1,
                product_degree,
                product.as_slice(),
                factor_degree,
                self.chi[factor_degree as usize][factor_idx].as_slice(),
            );
            product = next;
            product_degree = next_degree;
        }
        debug_assert_eq!(product_degree, degree);

        result.add(&product, if sign_exponent % 2 == 0 { 1 } else { p - 1 });
        result
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::algebra::{AdemAlgebra, MilnorAlgebra};

    /// The full coproduct of a basis element, assembled from `decompose` and `coproduct` the same way
    /// `TensorModule` does it: the coproduct of a product is the product of the coproducts.
    ///
    /// Returns a list of `(coefficient, left degree, left index, right degree, right index)`.
    fn full_coproduct<A: Algebra + Bialgebra>(
        algebra: &A,
        degree: i32,
        idx: usize,
    ) -> Vec<(u32, i32, usize, i32, usize)> {
        let p = algebra.prime();
        // Start from `1 (x) 1` and multiply in each factor's coproduct.
        let mut terms = vec![(1u32, 0, 0, 0, 0)];
        // `decompose` lists factors in the order they act, so the product is the reverse.
        for &(factor_degree, factor_idx) in algebra.decompose(degree, idx).iter().rev() {
            let mut next: Vec<(u32, i32, usize, i32, usize)> = Vec::new();
            for &(coeff, ld, li, rd, ri) in &terms {
                for (fld, fli, frd, fri) in algebra.coproduct(factor_degree, factor_idx) {
                    // Multiply `(left (x) right)` by `(fl (x) fr)` in `A (x) A`, which carries the
                    // Koszul sign `(a (x) b)(c (x) d) = (-1)^{|b||c|} ac (x) bd`. That sign is not
                    // optional bookkeeping here: `b P1 b` at an odd prime has two Bockstein factors
                    // and gets a genuine sign from it.
                    let sign = if (rd * fld) % 2 == 0 { 1 } else { p - 1 };
                    let mut left = FpVector::new(p, algebra.dimension(ld + fld));
                    algebra.multiply_basis_elements(left.as_slice_mut(), 1, ld, li, fld, fli);
                    let mut right = FpVector::new(p, algebra.dimension(rd + frd));
                    algebra.multiply_basis_elements(right.as_slice_mut(), 1, rd, ri, frd, fri);
                    for (lidx, lv) in left.iter_nonzero() {
                        for (ridx, rv) in right.iter_nonzero() {
                            next.push((
                                (coeff * lv * rv * sign) % p,
                                ld + fld,
                                lidx,
                                rd + frd,
                                ridx,
                            ));
                        }
                    }
                }
            }
            terms = next;
        }
        terms
    }

    /// The defining identity of the antipode: `sum chi(x') x'' = 0` for every basis element of
    /// positive degree.
    #[rstest]
    #[case(2, 20)]
    #[case(3, 30)]
    #[case(5, 30)]
    fn defining_identity_adem(#[case] p: u32, #[case] max_degree: i32) {
        let p = fp::prime::ValidPrime::new(p);
        let algebra = Arc::new(AdemAlgebra::new(p, false));
        algebra.compute_basis(max_degree);
        let antipode = Antipode::new(Arc::clone(&algebra));
        antipode.compute_through_degree(max_degree);

        for degree in 1..=max_degree {
            for idx in 0..algebra.dimension(degree) {
                let mut total = FpVector::new(p, algebra.dimension(degree));
                for (coeff, ld, li, rd, ri) in full_coproduct(&*algebra, degree, idx) {
                    antipode.algebra.multiply_element_by_basis_element(
                        total.as_slice_mut(),
                        coeff,
                        ld,
                        antipode.on_basis_element(ld, li).as_slice(),
                        rd,
                        ri,
                    );
                }
                assert!(
                    total.is_zero(),
                    "sum chi(x') x'' != 0 at p = {p} for {} (degree {degree}, index {idx}): got \
                     {total}",
                    algebra.basis_element_to_string(degree, idx),
                );
            }
        }
    }

    /// `chi` is an involution, which is a strong check that is independent of the recursion used to
    /// compute it.
    #[rstest]
    #[case(2, 24)]
    #[case(3, 36)]
    #[case(5, 40)]
    fn involution_adem(#[case] p: u32, #[case] max_degree: i32) {
        let p = fp::prime::ValidPrime::new(p);
        let algebra = Arc::new(AdemAlgebra::new(p, false));
        algebra.compute_basis(max_degree);
        let antipode = Antipode::new(Arc::clone(&algebra));

        for degree in 0..=max_degree {
            for idx in 0..algebra.dimension(degree) {
                let mut twice = FpVector::new(p, algebra.dimension(degree));
                antipode.apply(
                    twice.as_slice_mut(),
                    1,
                    degree,
                    antipode.on_basis_element(degree, idx).as_slice(),
                );
                let mut expected = FpVector::new(p, algebra.dimension(degree));
                expected.set_entry(idx, 1);
                assert_eq!(
                    twice,
                    expected,
                    "chi^2 != id at p = {p} on {} (degree {degree}, index {idx})",
                    algebra.basis_element_to_string(degree, idx),
                );
            }
        }
    }

    /// The Milnor basis only at the prime 2: `MilnorAlgebra`'s [`Bialgebra`] implementation asserts
    /// `p == 2` in `coproduct`, so the antipode is unavailable there at odd primes. That is a
    /// pre-existing limitation of the algebra, not of this recursion, and it does not affect the
    /// module builder, which works in the Adem basis.
    #[rstest]
    #[case(2, 20)]
    fn involution_milnor(#[case] p: u32, #[case] max_degree: i32) {
        let p = fp::prime::ValidPrime::new(p);
        let algebra = Arc::new(MilnorAlgebra::new(p, false));
        algebra.compute_basis(max_degree);
        let antipode = Antipode::new(Arc::clone(&algebra));

        for degree in 0..=max_degree {
            for idx in 0..algebra.dimension(degree) {
                let mut twice = FpVector::new(p, algebra.dimension(degree));
                antipode.apply(
                    twice.as_slice_mut(),
                    1,
                    degree,
                    antipode.on_basis_element(degree, idx).as_slice(),
                );
                let mut expected = FpVector::new(p, algebra.dimension(degree));
                expected.set_entry(idx, 1);
                assert_eq!(
                    twice, expected,
                    "chi^2 != id at p = {p} in the Milnor basis (degree {degree}, index {idx})"
                );
            }
        }
    }

    /// The classical values at the prime 2: `chi(Sq1) = Sq1`, `chi(Sq2) = Sq2`,
    /// `chi(Sq3) = Sq2 Sq1`, `chi(Sq4) = Sq4 + Sq3 Sq1`.
    #[test]
    fn known_values_at_two() {
        let p = fp::prime::ValidPrime::new(2);
        let algebra = Arc::new(AdemAlgebra::new(p, false));
        algebra.compute_basis(10);
        let antipode = Antipode::new(Arc::clone(&algebra));

        let chi_of = |name: &str| {
            let (degree, idx) = algebra.basis_element_from_string(name).unwrap();
            algebra.element_to_string(degree, antipode.on_basis_element(degree, idx).as_slice())
        };

        assert_eq!(chi_of("Sq1"), "Sq1");
        assert_eq!(chi_of("Sq2"), "Sq2");
        assert_eq!(chi_of("Sq3"), "Sq2 Sq1");
        assert_eq!(chi_of("Sq4"), "Sq4 + Sq3 Sq1");
    }

    /// At an odd prime the Bockstein is the only odd-degree generator, and `chi(b) = -b`.
    #[test]
    fn bockstein_at_three() {
        let p = fp::prime::ValidPrime::new(3);
        let algebra = Arc::new(AdemAlgebra::new(p, false));
        algebra.compute_basis(10);
        let antipode = Antipode::new(Arc::clone(&algebra));

        let (degree, idx) = algebra.basis_element_from_string("b").unwrap();
        assert_eq!(degree, 1);
        assert_eq!(
            algebra.element_to_string(degree, antipode.on_basis_element(degree, idx).as_slice()),
            "2 * b"
        );
    }
}
