//! Extensions of one module by another.
//!
//! Given modules $N$ and $Q$, this enumerates the extensions
//!
//! $$ 0 \to N \to M \to Q \to 0 $$
//!
//! exactly, by linear algebra alone — no resolution and no dependency on `ext`.
//!
//! # Why it is linear
//!
//! As a graded vector space $M_d = N_d \oplus Q_d$. Since $N$ is a submodule the action of an algebra
//! generator $g$ on the $N$ summand stays in $N$, and since $Q$ is the quotient the action on the $Q$
//! summand is the $Q$-action plus an unknown component
//!
//! $$ \theta_g \colon Q_d \to N_{d + |g|}. $$
//!
//! Those $\theta_g$ are the only free data. A composite of generator actions applied to an element of
//! the $Q$ summand either stays in $Q$ the whole way or crosses into $N$ at exactly one step and stays
//! there, because nothing maps $N$ back to $Q$. So every composite is affine in $\theta$, with at most
//! one factor of $\theta$ per term: a product of two $\theta$s would need a map $N \to Q$ and is zero.
//!
//! Hence the value of each Adem relation is *linear* in $\theta$, and vanishes at $\theta = 0$, where
//! $M$ is the split extension $N \oplus Q$ and is certainly a module. Writing $F$ for the map sending
//! $\theta$ to the tuple of relation values:
//!
//! - the valid $\theta$ are $\ker F$, the cocycles;
//! - changing the splitting by a degree-0 map $s \colon Q \to N$ replaces $\theta_g$ by
//!   $\theta_g + g_N \circ s - s \circ g_Q$, so the coboundaries are the image of $\delta$;
//! - $\Ext^1_A(Q, N) = \ker F / \operatorname{im} \delta$.
//!
//! $F$ and $\delta$ are computed one basis vector at a time: each column is obtained by building the
//! module for a single elementary $\theta$ and reading off what the relations come to.

use std::sync::Arc;

use algebra::{
    AdemAlgebra, Algebra, GeneratedAlgebra,
    module::{FDModule, Module},
};
use anyhow::{Context, ensure};
use bivec::BiVec;
use fp::{
    matrix::{AugmentedMatrix, Matrix, Subquotient, Subspace},
    prime::{Prime, ValidPrime},
    vector::FpVector,
};

/// Where each coordinate of $\theta$ lives.
///
/// A coordinate is a single matrix entry of one $\theta_g$: the map sending the `source_idx`th basis
/// element of $Q_{source\_degree}$ to the `target_idx`th basis element of
/// $N_{source\_degree + op\_degree}$.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Coordinate {
    op_degree: i32,
    op_idx: usize,
    source_degree: i32,
    source_idx: usize,
    target_idx: usize,
}

/// The extensions of `quotient` by `sub`, together with everything needed to build one.
pub struct Extensions {
    algebra: Arc<AdemAlgebra>,
    sub: FDModule<AdemAlgebra>,
    quotient: FDModule<AdemAlgebra>,
    /// The coordinates of $\theta$, in the order the vectors below are indexed by.
    coordinates: Vec<Coordinate>,
    /// A basis of $\Ext^1(Q, N)$, each element a representative cocycle $\theta$.
    classes: Vec<FpVector>,
    min_degree: i32,
    top_degree: i32,
}

impl Extensions {
    /// The dimension of $\Ext^1_A(Q, N)$, i.e. the number of independent extension classes.
    pub fn dimension(&self) -> usize {
        self.classes.len()
    }

    pub fn prime(&self) -> ValidPrime {
        self.algebra.prime()
    }

    /// The extension corresponding to the given combination of the basis classes.
    ///
    /// All-zero coefficients give the split extension $N \oplus Q$.
    pub fn realise(&self, coefficients: &[u32]) -> anyhow::Result<FDModule<AdemAlgebra>> {
        ensure!(
            coefficients.len() == self.classes.len(),
            "Expected {} coefficients, one per class, but got {}",
            self.classes.len(),
            coefficients.len()
        );
        let mut theta = FpVector::new(self.prime(), self.coordinates.len());
        for (coefficient, class) in coefficients.iter().zip(&self.classes) {
            theta.add(class, *coefficient % self.prime().as_u32());
        }
        Ok(self.build(&theta))
    }

    /// The dimension of `sub` in a degree, which is the size of the $N$ block there.
    fn sub_dimension(&self, degree: i32) -> usize {
        self.sub.dimension(degree)
    }

    /// Assemble $M$ for a given $\theta$.
    ///
    /// The basis of $M_d$ is the basis of $N_d$ followed by the basis of $Q_d$, so the $N$ block is an
    /// initial segment and the quotient map is "forget the first `sub_dimension(d)` coordinates".
    fn build(&self, theta: &FpVector) -> FDModule<AdemAlgebra> {
        let mut graded_dimension = BiVec::with_capacity(self.min_degree, self.top_degree + 1);
        for degree in self.min_degree..=self.top_degree {
            graded_dimension.push(self.sub.dimension(degree) + self.quotient.dimension(degree));
        }
        let mut module = FDModule::new(Arc::clone(&self.algebra), String::new(), graded_dimension);

        for degree in self.min_degree..=self.top_degree {
            let offset = self.sub_dimension(degree);
            for idx in 0..offset {
                module.set_basis_element_name(
                    degree,
                    idx,
                    self.sub.basis_element_to_string(degree, idx),
                );
            }
            for idx in 0..self.quotient.dimension(degree) {
                module.set_basis_element_name(
                    degree,
                    offset + idx,
                    self.quotient.basis_element_to_string(degree, idx),
                );
            }
        }

        // The two diagonal blocks: `N` acts as it does in `N`, `Q` as it does in `Q`.
        for source_degree in self.min_degree..=self.top_degree {
            for op_degree in 1..=(self.top_degree - source_degree) {
                let Some(op_idx) = generator_index(&self.algebra, op_degree) else {
                    continue;
                };
                let target_degree = source_degree + op_degree;
                if module.dimension(target_degree) == 0 {
                    continue;
                }
                let target_offset = self.sub_dimension(target_degree);

                for idx in 0..self.sub.dimension(source_degree) {
                    if self.sub.dimension(target_degree) > 0 {
                        let image = self
                            .sub
                            .action(op_degree, op_idx, source_degree, idx)
                            .clone();
                        let row = module.action_mut(op_degree, op_idx, source_degree, idx);
                        for (target_idx, value) in image.iter_nonzero() {
                            row.add_basis_element(target_idx, value);
                        }
                    }
                }

                let source_offset = self.sub_dimension(source_degree);
                for idx in 0..self.quotient.dimension(source_degree) {
                    if self.quotient.dimension(target_degree) > 0 {
                        let image = self
                            .quotient
                            .action(op_degree, op_idx, source_degree, idx)
                            .clone();
                        let row = module.action_mut(
                            op_degree,
                            op_idx,
                            source_degree,
                            source_offset + idx,
                        );
                        for (target_idx, value) in image.iter_nonzero() {
                            row.add_basis_element(target_offset + target_idx, value);
                        }
                    }
                }
            }
        }

        // The off-diagonal part, which is what `theta` parameterises.
        for (index, coordinate) in self.coordinates.iter().enumerate() {
            let value = theta.entry(index);
            if value == 0 {
                continue;
            }
            let source_offset = self.sub_dimension(coordinate.source_degree);
            module
                .action_mut(
                    coordinate.op_degree,
                    coordinate.op_idx,
                    coordinate.source_degree,
                    source_offset + coordinate.source_idx,
                )
                .add_basis_element(coordinate.target_idx, value);
        }

        for source_degree in (self.min_degree..=self.top_degree).rev() {
            for target_degree in (source_degree + 1)..=self.top_degree {
                module.extend_actions(source_degree, target_degree);
            }
        }
        module
    }

    /// The value of every generating relation on every basis element, concatenated.
    ///
    /// This is `FDModule::check_validity` computing vectors instead of strings: the module satisfies
    /// the relations exactly when this is zero. The layout depends only on the graded dimensions and
    /// the algebra, so it is the same for every `theta` and the results can be compared entry by
    /// entry.
    fn relation_values(&self, module: &FDModule<AdemAlgebra>) -> FpVector {
        let p = self.prime();
        let mut values = Vec::new();
        let mut output = FpVector::new(p, 0);
        let mut intermediate = FpVector::new(p, 0);

        for input_degree in self.min_degree..=self.top_degree {
            for output_degree in (input_degree + 1)..=self.top_degree {
                let op_degree = output_degree - input_degree;
                let output_dimension = module.dimension(output_degree);
                output.set_scratch_vector_size(output_dimension);
                for input_idx in 0..module.dimension(input_degree) {
                    for relation in self.algebra.generating_relations(op_degree) {
                        output.set_to_zero();
                        for (coeff, (deg_1, idx_1), (deg_2, idx_2)) in relation {
                            intermediate
                                .set_scratch_vector_size(module.dimension(input_degree + deg_2));
                            module.act_on_basis(
                                intermediate.as_slice_mut(),
                                1,
                                deg_2,
                                idx_2,
                                input_degree,
                                input_idx,
                            );
                            module.act(
                                output.as_slice_mut(),
                                coeff,
                                deg_1,
                                idx_1,
                                input_degree + deg_2,
                                intermediate.as_slice(),
                            );
                        }
                        values.extend(output.iter());
                    }
                }
            }
        }
        FpVector::from_slice(p, &values)
    }
}

/// The single algebra generator in `op_degree`, if there is one.
fn generator_index(algebra: &Arc<AdemAlgebra>, op_degree: i32) -> Option<usize> {
    if op_degree <= 0 {
        return None;
    }
    algebra.compute_basis(op_degree);
    algebra.generators(op_degree).first().copied()
}

/// Enumerate the extensions of `quotient` by `sub`.
pub fn compute(
    algebra: Arc<AdemAlgebra>,
    sub: FDModule<AdemAlgebra>,
    quotient: FDModule<AdemAlgebra>,
) -> anyhow::Result<Extensions> {
    let p = algebra.prime();
    let sub_min = sub.min_degree();
    let quotient_min = quotient.min_degree();
    let sub_top = sub.max_degree().context("The submodule must be bounded")?;
    let quotient_top = quotient
        .max_degree()
        .context("The quotient must be bounded")?;

    // A zero module has `max_degree() == min_degree() - 1`, so guard before taking the span.
    let min_degree = if sub.total_dimension() == 0 {
        quotient_min
    } else if quotient.total_dimension() == 0 {
        sub_min
    } else {
        sub_min.min(quotient_min)
    };
    let top_degree = sub_top.max(quotient_top);
    ensure!(
        top_degree >= min_degree,
        "Both modules are zero, so there is nothing to extend"
    );
    algebra.compute_basis(top_degree - min_degree);

    // The coordinates of theta.
    let mut coordinates = Vec::new();
    for source_degree in min_degree..=top_degree {
        for op_degree in 1..=(top_degree - source_degree) {
            let Some(op_idx) = generator_index(&algebra, op_degree) else {
                continue;
            };
            let target_degree = source_degree + op_degree;
            for source_idx in 0..quotient.dimension(source_degree) {
                for target_idx in 0..sub.dimension(target_degree) {
                    coordinates.push(Coordinate {
                        op_degree,
                        op_idx,
                        source_degree,
                        source_idx,
                        target_idx,
                    });
                }
            }
        }
    }

    let mut extensions = Extensions {
        algebra: Arc::clone(&algebra),
        sub,
        quotient,
        coordinates,
        classes: Vec::new(),
        min_degree,
        top_degree,
    };

    let width = extensions.coordinates.len();
    // The split extension is a module, so the relation values vanish there and `F` is linear rather
    // than merely affine. Its length also fixes the number of columns.
    let zero = FpVector::new(p, width);
    let baseline = extensions.relation_values(&extensions.build(&zero));
    debug_assert!(baseline.is_zero(), "the split extension must be a module");
    let height = baseline.len();

    if width == 0 {
        // Nothing to vary: the split extension is the only one.
        extensions.classes = Vec::new();
        return Ok(extensions);
    }

    // `F`, one row per coordinate, augmented with the identity so that `compute_kernel` can express
    // the kernel back in terms of theta.
    let mut matrix = AugmentedMatrix::<2>::new(p, width, [height, width]);
    for index in 0..width {
        let mut elementary = FpVector::new(p, width);
        elementary.set_entry(index, 1);
        let values = extensions.relation_values(&extensions.build(&elementary));
        let mut row = matrix.row_segment_mut(index, 0, 0);
        for (column, value) in values.iter_nonzero() {
            row.set_entry(column, value);
        }
    }
    matrix.segment(1, 1).add_identity();
    matrix.row_reduce();
    let cocycles = matrix.compute_kernel();

    // The coboundaries: `delta(s)_g = g_N . s - s . g_Q` for a degree-0 map `s: Q -> N`.
    let coboundaries = extensions.coboundaries();

    let ext = Subquotient::from_parts(cocycles, coboundaries);
    extensions.classes = ext.gens().map(|class| class.to_owned()).collect();
    Ok(extensions)
}

impl Extensions {
    /// The image of the coboundary map, as a subspace of the space of thetas.
    ///
    /// A degree-0 map $s \colon Q \to N$ is a choice of splitting, and two thetas differing by
    /// $\delta s$ give isomorphic extensions.
    fn coboundaries(&self) -> Subspace {
        let p = self.prime();
        let width = self.coordinates.len();

        // One generator of the source for each matrix entry of a degree-preserving `s`.
        let mut splittings = Vec::new();
        for degree in self.min_degree..=self.top_degree {
            for source_idx in 0..self.quotient.dimension(degree) {
                for target_idx in 0..self.sub.dimension(degree) {
                    splittings.push((degree, source_idx, target_idx));
                }
            }
        }
        if splittings.is_empty() || width == 0 {
            return Subspace::new(p, width);
        }

        let mut matrix = Matrix::new(p, splittings.len(), width);
        for (row_idx, &(degree, source_idx, target_idx)) in splittings.iter().enumerate() {
            let mut row = matrix.row_mut(row_idx);
            for (column, coordinate) in self.coordinates.iter().enumerate() {
                let mut value = 0;
                // `g_N . s`: `s` sends this source element into `N`, then `g` acts inside `N`.
                if coordinate.source_degree == degree && coordinate.source_idx == source_idx {
                    let image = self.sub.action(
                        coordinate.op_degree,
                        coordinate.op_idx,
                        degree,
                        target_idx,
                    );
                    value += image.entry(coordinate.target_idx);
                }
                // `s . g_Q`: `g` acts inside `Q` first, then `s` sends the result into `N`. This lands
                // on `target_idx` only when `s` was chosen there, i.e. in the target degree.
                if coordinate.source_degree + coordinate.op_degree == degree
                    && coordinate.target_idx == target_idx
                {
                    let image = self.quotient.action(
                        coordinate.op_degree,
                        coordinate.op_idx,
                        coordinate.source_degree,
                        coordinate.source_idx,
                    );
                    value += (p.as_u32() - 1) * image.entry(source_idx);
                }
                if value % p.as_u32() != 0 {
                    row.set_entry(column, value % p.as_u32());
                }
            }
        }
        Subspace::from_matrix(matrix)
    }
}

#[cfg(test)]
mod tests {
    use serde_json::{Value, json};

    use super::*;
    use crate::builder::Builder;

    fn builder(json: &Value) -> Builder {
        Builder::from_json(json).unwrap()
    }

    /// A single cell in the given degree, i.e. a shifted copy of the sphere.
    fn cell(p: u32, degree: i32) -> Builder {
        let mut builder = Builder::new(ValidPrime::new(p));
        builder.add_generator(degree, Some("y")).unwrap();
        builder
    }

    /// `Ext^{1,t}(F_p, F_p)` is spanned by the `h_i`, which sit exactly in the degrees of the algebra
    /// generators. Enumerating the extensions of a cell in degree 0 by a cell in degree `t` must see
    /// the same thing: one class when `t` is a generator degree, none otherwise.
    #[test]
    fn extensions_of_a_sphere_see_the_h_i() {
        for (degree, expected) in [(1, 1), (2, 1), (3, 0), (4, 1), (5, 0), (6, 0), (8, 1)] {
            let quotient = cell(2, 0);
            let extensions = quotient.extensions(&cell(2, degree)).unwrap();
            assert_eq!(
                extensions.dimension(),
                expected,
                "Ext^1 into degree {degree} should have dimension {expected}"
            );
        }
    }

    /// The non-split extension of the sphere by its own suspension is `C(2)`.
    #[test]
    fn the_nonsplit_extension_is_c2() {
        let quotient = cell(2, 0);
        let extensions = quotient.extensions(&cell(2, 1)).unwrap();
        assert_eq!(extensions.dimension(), 1);

        let module = extensions.realise(&[1]).unwrap();
        let expected = builder(&json!({
            "p": 2,
            "type": "finite dimensional module",
            "gens": { "x0": 0, "x1": 1 },
            "actions": ["Sq1 x0 = x1"],
        }));
        let (expected, _) = expected.build();
        module.test_equal(&expected).unwrap();

        // The zero class is the direct sum, which has no action at all.
        let split = extensions.realise(&[0]).unwrap();
        let mut json = json!({});
        split.to_json(&mut json);
        assert_eq!(json["actions"], json!([]));
    }

    /// Extending by a cell two degrees up gives `C(eta)`.
    #[test]
    fn extending_by_two_degrees_gives_ceta() {
        let quotient = cell(2, 0);
        let extensions = quotient.extensions(&cell(2, 2)).unwrap();
        let module = extensions.realise(&[1]).unwrap();

        let (expected, _) = builder(&json!({
            "p": 2,
            "type": "finite dimensional module",
            "gens": { "x0": 0, "x2": 2 },
            "actions": ["Sq2 x0 = x2"],
        }))
        .build();
        module.test_equal(&expected).unwrap();
    }

    /// Coboundaries have to be quotiented out, not just counted.
    ///
    /// Extending the sphere by `C(2)` has a one dimensional space of cocycles — the choice of
    /// `Sq1 x0`, landing in the top cell of `C(2)` — but it is entirely a coboundary, since changing
    /// the splitting by the degree-0 map into the bottom cell of `C(2)` moves it. So every such
    /// extension splits and `Ext^1` is zero.
    #[test]
    fn coboundaries_are_quotiented_out() {
        let quotient = cell(2, 0);
        let sub = builder(&json!({
            "p": 2,
            "type": "finite dimensional module",
            "gens": { "x0": 0, "x1": 1 },
            "actions": ["Sq1 x0 = x1"],
        }));
        assert_eq!(quotient.extensions(&sub).unwrap().dimension(), 0);
    }

    /// Every extension enumerated must actually be a module.
    #[test]
    fn every_extension_satisfies_the_relations() {
        let quotient = builder(&json!({
            "p": 2,
            "type": "finite dimensional module",
            "gens": { "x0": 0, "x1": 1 },
            "actions": ["Sq1 x0 = x1"],
        }));
        let sub = cell(2, 3);
        let extensions = quotient.extensions(&sub).unwrap();

        // At the prime 2 the classes are indexed by subsets, so check them all.
        for mask in 0..(1u32 << extensions.dimension()) {
            let coefficients: Vec<u32> = (0..extensions.dimension())
                .map(|i| (mask >> i) & 1)
                .collect();
            let module = extensions.realise(&coefficients).unwrap();
            let min = module.min_degree();
            let top = module.max_degree().unwrap();
            for input in min..=top {
                for output in (input + 1)..=top {
                    let failures = module.check_validity_all(input, output);
                    assert!(
                        failures.is_empty(),
                        "extension {coefficients:?} fails a relation: {failures:?}"
                    );
                }
            }
        }
    }

    /// At an odd prime the Bockstein plays the role `Sq1` does at 2, so extending a cell by the next
    /// degree up gives the mod `p` Moore space.
    #[test]
    fn odd_prime_moore_space() {
        for p in [3, 5] {
            let quotient = cell(p, 0);
            let extensions = quotient.extensions(&cell(p, 1)).unwrap();
            assert_eq!(extensions.dimension(), 1, "at p = {p}");

            let module = extensions.realise(&[1]).unwrap();
            let mut json = json!({});
            module.to_json(&mut json);
            assert_eq!(json["actions"], json!(["b y = y"]), "at p = {p}");

            // Degree 2 carries no generator at an odd prime, so there is nothing there.
            assert_eq!(cell(p, 0).extensions(&cell(p, 2)).unwrap().dimension(), 0);
        }
    }

    /// A module with a `profile`, or a mismatched prime, must be refused rather than answered wrongly.
    #[test]
    fn refuses_what_it_cannot_answer() {
        let quotient = cell(2, 0);
        assert!(quotient.extensions(&cell(3, 1)).is_err());

        let profiled = builder(&json!({
            "p": 2,
            "algebra": ["milnor"],
            "profile": { "truncated": true, "p_part": [3, 2, 1] },
            "type": "finite dimensional module",
            "gens": { "x0": 0, "x2": 2 },
            "actions": ["Sq2 x0 = x2"],
        }));
        assert!(quotient.extensions(&profiled).is_err());
        assert!(profiled.extensions(&quotient).is_err());
    }

    /// `apply_extension` replaces the module in place and names the result after its parts.
    #[test]
    fn apply_extension_replaces_the_module() {
        let mut quotient = cell(2, 0);
        quotient.set_name("S");
        let mut sub = cell(2, 1);
        sub.set_name("S1");

        quotient.apply_extension(&sub, &[1]).unwrap();
        let json = quotient.to_json();
        // Both cells were called `y`; names must be unique across the module, so the one met second
        // — the submodule's, in the higher degree — is renamed.
        assert_eq!(json["actions"], json!(["Sq1 y = x1"]));
        assert_eq!(json["name"], json!("S1 . S"));
    }
}
