//! The module builder proper.
//!
//! This module holds all of the logic and deliberately does not mention `wasm_bindgen`, so that it
//! can be unit tested on the host.
//!
//! # The source of truth
//!
//! The state we keep is the list of basis elements and the list of *generator* actions the user has
//! entered, not an [`FDModule`]. Every query rebuilds the [`FDModule`] from scratch, in the same
//! order [`FDModule::from_json`] does it. This is because [`FDModule::extend_actions`] writes into
//! the same action table it reads from, so mutating a module in place after an edit is hard to
//! reason about, whereas a rebuild is exactly reproducible. At the sizes involved — tens of basis
//! elements — a rebuild costs microseconds.
//!
//! Actions are keyed by generator *name* rather than by index, so that adding, removing and
//! renaming basis elements cannot invalidate them.

use std::{collections::BTreeMap, sync::Arc};

use algebra::{
    AdemAlgebra, Algebra, GeneratedAlgebra,
    module::{FDModule, Module},
    steenrod_evaluator::SteenrodEvaluator,
};
use anyhow::{Context, anyhow, bail, ensure};
use bivec::BiVec;
use fp::{
    prime::{Prime, ValidPrime},
    vector::FpVector,
};
use serde_json::{Map, Value, json};

/// The json fields the builder owns. Every other field of a loaded file is preserved verbatim so
/// that `profile`, `products`, `self_maps`, `shift` and `cofiber` survive a load and save cycle.
const OWNED_FIELDS: [&str; 5] = ["p", "type", "gens", "actions", "name"];

/// The value of the `type` field we read and write.
const MODULE_TYPE: &str = "finite dimensional module";

/// A relation the module fails to satisfy, decorated with the names the interface needs in order to
/// point at the cells responsible.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Failure {
    /// The degree of the basis element the relation was applied to.
    pub input_degree: i32,
    /// The index of that basis element within its degree.
    pub input_idx: usize,
    /// The name of that basis element.
    pub input_name: String,
    /// The total degree of the relation.
    pub op_degree: i32,
    /// The relation that ought to act as zero.
    pub relation: String,
    /// The value it takes instead.
    pub value: String,
}

impl Failure {
    fn to_json(&self) -> Value {
        json!({
            "input_degree": self.input_degree,
            "input_idx": self.input_idx,
            "input_name": self.input_name,
            "op_degree": self.op_degree,
            "relation": self.relation,
            "value": self.value,
        })
    }
}

/// An interactively editable finite dimensional module over the Steenrod algebra.
pub struct Builder {
    p: ValidPrime,
    /// The Adem algebra. We edit in the Adem basis because its generators are the minimal generating
    /// set — the $Sq^{2^k}$, or $\beta$ and $P^{p^k}$ at odd primes — and because its
    /// `generating_relations` are the Adem relations. The generator names it produces (`Sq1`, `P1`,
    /// `b`) are also accepted by `MilnorAlgebra::basis_element_from_string`, so the files we write
    /// stay loadable in either basis.
    algebra: Arc<AdemAlgebra>,
    /// Used to evaluate arbitrary Steenrod expressions typed by the user. Its `AdemAlgebra` is a
    /// separate instance from `algebra`, but the Adem basis is generated deterministically from the
    /// prime, so the two agree on basis indices.
    evaluator: SteenrodEvaluator,
    name: String,
    /// The names of the basis elements, indexed by internal degree. `gens[d][i]` names the `i`th
    /// basis element of degree `d`. Empty degrees at either end are trimmed, so `min_degree()` and
    /// `max_degree()` are tight.
    gens: BiVec<Vec<String>>,
    /// The user-entered action of the algebra generators, keyed by `(op_degree, source name)`, with
    /// the value a list of `(coefficient, target name)` pairs.
    actions: BTreeMap<(i32, String), Vec<(u32, String)>>,
    /// Fields of a loaded json file that we do not interpret and must not lose.
    extra: Map<String, Value>,
}

impl Builder {
    pub fn new(p: ValidPrime) -> Self {
        Self {
            p,
            algebra: Arc::new(AdemAlgebra::new(p, false)),
            evaluator: SteenrodEvaluator::new(p),
            name: String::new(),
            gens: BiVec::new(0),
            actions: BTreeMap::new(),
            extra: Map::new(),
        }
    }

    pub fn prime(&self) -> ValidPrime {
        self.p
    }

    /// Replace the prime.
    ///
    /// The basis elements are kept, but the actions are discarded: which operations are algebra
    /// generators depends on the prime, so an action entered at one prime is generally meaningless
    /// at another.
    pub fn set_prime(&mut self, p: ValidPrime) {
        self.p = p;
        self.algebra = Arc::new(AdemAlgebra::new(p, false));
        self.evaluator = SteenrodEvaluator::new(p);
        self.actions.clear();
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn set_name(&mut self, name: &str) {
        self.name = name.to_string();
    }

    /// The lowest degree carrying a basis element, or 0 if there are none.
    pub fn min_degree(&self) -> i32 {
        self.gens.min_degree()
    }

    /// The highest degree carrying a basis element. This is `min_degree() - 1` when the module is
    /// zero, matching `BiVec::max_degree`.
    pub fn top_degree(&self) -> i32 {
        self.gens.max_degree()
    }

    pub fn is_zero(&self) -> bool {
        self.gens.is_empty()
    }

    pub fn dimension(&self, degree: i32) -> usize {
        self.gens.get(degree).map_or(0, Vec::len)
    }

    /// The largest operation degree that can act non-trivially, i.e. the width of the module.
    fn max_op_degree(&self) -> i32 {
        if self.is_zero() {
            0
        } else {
            self.top_degree() - self.min_degree()
        }
    }

    fn name_exists(&self, name: &str) -> bool {
        self.gens.iter().any(|v| v.iter().any(|n| n == name))
    }

    /// The degree and index of the basis element with the given name.
    fn lookup(&self, name: &str) -> Option<(i32, usize)> {
        self.gens
            .iter_enum()
            .find_map(|(d, v)| v.iter().position(|n| n == name).map(|i| (d, i)))
    }

    fn require(&self, degree: i32, idx: usize) -> anyhow::Result<&str> {
        self.gens
            .get(degree)
            .and_then(|v| v.get(idx))
            .map(String::as_str)
            .ok_or_else(|| anyhow!("No basis element with index {idx} in degree {degree}"))
    }

    /// The index of the single algebra generator in `op_degree`, if there is one.
    ///
    /// `AdemAlgebra` has at most one generator per degree: $Sq^{2^k}$ at the prime 2, and $\beta$ in
    /// degree 1 together with $P^{p^k}$ at odd primes.
    fn generator_index(&self, op_degree: i32) -> Option<usize> {
        if op_degree <= 0 {
            return None;
        }
        self.algebra.compute_basis(op_degree);
        self.algebra.generators(op_degree).first().copied()
    }

    /// The operation degrees that carry an algebra generator and can act non-trivially on this
    /// module, together with the generator's name.
    pub fn generator_degrees(&self) -> Vec<(i32, String)> {
        (1..=self.max_op_degree())
            .filter_map(|d| {
                let idx = self.generator_index(d)?;
                Some((d, self.algebra.generator_to_string(d, idx)))
            })
            .collect()
    }

    // ---------------------------------------------------------------- editing

    /// Add a basis element in `degree`.
    ///
    /// If `name` is `None` a fresh name is generated. Returns the name used.
    pub fn add_generator(&mut self, degree: i32, name: Option<&str>) -> anyhow::Result<String> {
        let name = match name {
            Some(name) => {
                let name = name.trim();
                validate_name(name)?;
                ensure!(
                    !self.name_exists(name),
                    "There is already a basis element called {name}"
                );
                name.to_string()
            }
            None => self.default_name(degree),
        };

        if self.is_zero() {
            self.gens = BiVec::new(degree);
            self.gens.push(Vec::new());
        } else {
            self.gens.extend_negative(degree, Vec::new());
            self.gens.extend_with(degree, |_| Vec::new());
        }
        self.gens[degree].push(name.clone());
        Ok(name)
    }

    /// Remove a basis element, along with every action mentioning it.
    pub fn remove_generator(&mut self, degree: i32, idx: usize) -> anyhow::Result<()> {
        let name = self.require(degree, idx)?.to_string();
        self.gens[degree].remove(idx);

        self.actions.retain(|(_, source), _| source != &name);
        for targets in self.actions.values_mut() {
            targets.retain(|(_, target)| target != &name);
        }
        self.actions.retain(|_, targets| !targets.is_empty());

        self.normalize();
        Ok(())
    }

    /// Rename a basis element, rewriting the actions that refer to it.
    pub fn rename_generator(
        &mut self,
        degree: i32,
        idx: usize,
        new_name: &str,
    ) -> anyhow::Result<()> {
        let new_name = new_name.trim();
        validate_name(new_name)?;
        let old_name = self.require(degree, idx)?.to_string();
        if old_name == new_name {
            return Ok(());
        }
        ensure!(
            !self.name_exists(new_name),
            "There is already a basis element called {new_name}"
        );

        self.gens[degree][idx] = new_name.to_string();
        self.actions = std::mem::take(&mut self.actions)
            .into_iter()
            .map(|((op_degree, source), targets)| {
                let source = if source == old_name {
                    new_name.to_string()
                } else {
                    source
                };
                let targets = targets
                    .into_iter()
                    .map(|(coeff, target)| {
                        if target == old_name {
                            (coeff, new_name.to_string())
                        } else {
                            (coeff, target)
                        }
                    })
                    .collect();
                ((op_degree, source), targets)
            })
            .collect();
        Ok(())
    }

    /// Set the action of the generator in `op_degree` on a basis element.
    ///
    /// `coeffs` is a coefficient vector for the basis of `source_degree + op_degree`.
    pub fn set_action(
        &mut self,
        op_degree: i32,
        source_degree: i32,
        source_idx: usize,
        coeffs: &[u32],
    ) -> anyhow::Result<()> {
        ensure!(
            self.generator_index(op_degree).is_some(),
            "Degree {op_degree} carries no algebra generator, so there is no action to set. At \
             the prime {} the generators are {}.",
            self.p,
            self.generator_description(),
        );
        let source = self.require(source_degree, source_idx)?.to_string();

        let output_degree = source_degree + op_degree;
        let output_dim = self.dimension(output_degree);
        ensure!(
            coeffs.len() == output_dim,
            "Expected {output_dim} coefficients for degree {output_degree}, got {}",
            coeffs.len()
        );

        let targets: Vec<(u32, String)> = coeffs
            .iter()
            .enumerate()
            .filter(|&(_, &c)| c % self.p.as_u32() != 0)
            .map(|(i, &c)| (c % self.p.as_u32(), self.gens[output_degree][i].clone()))
            .collect();

        if targets.is_empty() {
            self.actions.remove(&(op_degree, source));
        } else {
            self.actions.insert((op_degree, source), targets);
        }
        Ok(())
    }

    /// Add `coeff` to one entry of an action, reducing mod `p`.
    ///
    /// This is what dragging an arc from one cell to another does: at the prime 2 it toggles the
    /// entry, and at odd primes repeated use cycles through the coefficients.
    pub fn add_to_action(
        &mut self,
        op_degree: i32,
        source_degree: i32,
        source_idx: usize,
        target_idx: usize,
        coeff: u32,
    ) -> anyhow::Result<()> {
        let output_degree = source_degree + op_degree;
        let output_dim = self.dimension(output_degree);
        ensure!(
            target_idx < output_dim,
            "No basis element with index {target_idx} in degree {output_degree}"
        );

        let mut coeffs = self.action_coefficients(op_degree, source_degree, source_idx)?;
        coeffs[target_idx] = (coeffs[target_idx] + coeff) % self.p.as_u32();
        self.set_action(op_degree, source_degree, source_idx, &coeffs)
    }

    /// The currently stored coefficient vector of an action, as a dense vector.
    pub fn action_coefficients(
        &self,
        op_degree: i32,
        source_degree: i32,
        source_idx: usize,
    ) -> anyhow::Result<Vec<u32>> {
        let source = self.require(source_degree, source_idx)?;
        let output_degree = source_degree + op_degree;
        let mut coeffs = vec![0; self.dimension(output_degree)];
        if let Some(targets) = self.actions.get(&(op_degree, source.to_string())) {
            for (coeff, target) in targets {
                // A target that no longer exists in this degree has been removed from `gens`; the
                // stored action is stale and simply contributes nothing.
                if let Some(idx) = self.gens[output_degree].iter().position(|n| n == target) {
                    coeffs[idx] = (coeffs[idx] + coeff) % self.p.as_u32();
                }
            }
        }
        Ok(coeffs)
    }

    /// Replace every action with the ones parsed from `text`, one action per line.
    ///
    /// Lines are in the same syntax as the `actions` field of a module file, e.g. `Sq1 x0 = x1` or
    /// `P1 x0 = 2 x4`. Blank lines and lines starting with `#` are ignored. Parsing is done by
    /// [`FDModule::parse_action`], so the syntax accepted here is exactly the syntax accepted in a
    /// file.
    pub fn set_actions_text(&mut self, text: &str) -> anyhow::Result<()> {
        let lines: Vec<&str> = text
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .collect();
        self.set_actions(&lines)
    }

    /// Replace every action with the ones in `actions`, each an entry of a file's `actions` array.
    pub fn set_actions(&mut self, actions: &[impl AsRef<str>]) -> anyhow::Result<()> {
        let mut module = self.skeleton();
        let lookup = self.lookup_table();
        for action in actions {
            let action = action.as_ref();
            module
                .parse_action(&lookup, action, false)
                .with_context(|| format!("Failed to parse action: {action}"))?;
        }
        self.harvest(&module);
        Ok(())
    }

    /// Trim degrees that no longer carry a basis element off either end.
    fn normalize(&mut self) {
        let occupied: Vec<i32> = self
            .gens
            .iter_enum()
            .filter(|(_, names)| !names.is_empty())
            .map(|(d, _)| d)
            .collect();
        let (Some(&min), Some(&max)) = (occupied.first(), occupied.last()) else {
            self.gens = BiVec::new(0);
            return;
        };
        if min == self.min_degree() && max == self.top_degree() {
            return;
        }
        let mut trimmed = BiVec::with_capacity(min, max + 1);
        for d in min..=max {
            trimmed.push(self.gens.get(d).cloned().unwrap_or_default());
        }
        self.gens = trimmed;
    }

    fn default_name(&self, degree: i32) -> String {
        // A leading `-` would not be a valid identifier, and the action syntax splits on
        // whitespace, so negative degrees get an underscore instead.
        let suffix = if degree < 0 {
            format!("_{}", -degree)
        } else {
            degree.to_string()
        };
        let base = format!("x{suffix}");
        if !self.name_exists(&base) {
            return base;
        }
        (1..)
            .map(|i| format!("{base}_{i}"))
            .find(|candidate| !self.name_exists(candidate))
            .expect("the range 1.. is infinite, so some candidate name is free")
    }

    /// Why this module cannot be checked against the Adem relations, if it cannot.
    ///
    /// A `profile` field means the module is defined over a sub-Hopf-algebra of the Steenrod algebra
    /// rather than over the whole thing: `ext/steenrod_modules/tmf2_sm_DA1.json` is an $A(2)$-module,
    /// on which $Sq^8$ simply does not act. Checking the full Adem relations against such a module
    /// reports failures that are not mistakes, so we decline to check it and say why instead.
    ///
    /// Editing over a proper subalgebra needs the profile-aware Milnor generating set — which
    /// includes the $Q_k$ and $P^s_t$, not just the $Sq^{2^k}$ — so it is not supported yet. Such a
    /// module can still be viewed, edited and saved without losing its `profile`.
    pub fn restriction(&self) -> Option<String> {
        if self.extra.get("profile").is_some_and(|v| !v.is_null()) {
            Some(
                "This module is defined over a sub-Hopf-algebra of the Steenrod algebra (it has a \
                 `profile`), so the Adem relations for the full algebra do not apply to it and \
                 are not checked. It can be edited and saved, but the relation checker is \
                 unavailable."
                    .to_string(),
            )
        } else {
            None
        }
    }

    /// A human-readable description of the algebra generators, for error messages.
    fn generator_description(&self) -> String {
        if self.p == 2 {
            "the Sq^{2^k}".to_string()
        } else {
            format!("b, and P^({}^k)", self.p)
        }
    }

    // --------------------------------------------------------------- building

    /// The module with the right graded dimensions and names, but no actions set.
    fn skeleton(&self) -> FDModule<AdemAlgebra> {
        let mut graded_dimension = BiVec::with_capacity(self.min_degree(), self.gens.len());
        for names in self.gens.iter() {
            graded_dimension.push(names.len());
        }
        let mut module = FDModule::new(
            Arc::clone(&self.algebra),
            self.name.clone(),
            graded_dimension,
        );
        for (degree, names) in self.gens.iter_enum() {
            for (idx, name) in names.iter().enumerate() {
                module.set_basis_element_name(degree, idx, name.clone());
            }
        }
        module
    }

    /// A `name -> (degree, index)` lookup of the shape [`FDModule::parse_action`] wants.
    ///
    /// The table is owned rather than borrowing `self`, so that the caller can go on to mutate the
    /// builder while the closure is alive.
    fn lookup_table(&self) -> impl for<'a> Fn(&'a str) -> anyhow::Result<(i32, usize)> + use<> {
        let table: BTreeMap<String, (i32, usize)> = self
            .gens
            .iter_enum()
            .flat_map(|(degree, names)| {
                names
                    .iter()
                    .enumerate()
                    .map(move |(idx, name)| (name.clone(), (degree, idx)))
            })
            .collect();
        move |name| {
            table
                .get(name)
                .copied()
                .ok_or_else(|| anyhow!("Invalid generator: {name}"))
        }
    }

    /// Read the generator actions of `module` back into our name-keyed storage.
    fn harvest(&mut self, module: &FDModule<AdemAlgebra>) {
        self.actions.clear();
        for (source_degree, names) in self.gens.iter_enum() {
            for op_degree in 1..=(self.top_degree() - source_degree) {
                let Some(op_idx) = self.generator_index(op_degree) else {
                    continue;
                };
                let output_degree = source_degree + op_degree;
                if self.dimension(output_degree) == 0 {
                    continue;
                }
                for (source_idx, source) in names.iter().enumerate() {
                    let vector = module.action(op_degree, op_idx, source_degree, source_idx);
                    let targets: Vec<(u32, String)> = vector
                        .iter_nonzero()
                        .map(|(idx, coeff)| (coeff, self.gens[output_degree][idx].clone()))
                        .collect();
                    if !targets.is_empty() {
                        self.actions.insert((op_degree, source.clone()), targets);
                    }
                }
            }
        }
    }

    /// Build the module, together with every relation it fails.
    ///
    /// The `extend_actions` and `check_validity_all` loop mirrors [`FDModule::from_json`]: the input
    /// degree descends, because extending at `(input, output)` consumes the already-extended actions
    /// of *larger* input degrees.
    pub fn build(&self) -> (FDModule<AdemAlgebra>, Vec<Failure>) {
        let mut module = self.skeleton();

        for ((op_degree, source), targets) in &self.actions {
            let Some((source_degree, source_idx)) = self.lookup(source) else {
                continue;
            };
            let Some(op_idx) = self.generator_index(*op_degree) else {
                continue;
            };
            let output_degree = source_degree + op_degree;
            if self.dimension(output_degree) == 0 {
                continue;
            }
            let mut coeffs = vec![0; self.dimension(output_degree)];
            for (coeff, target) in targets {
                if let Some(idx) = self.gens[output_degree].iter().position(|n| n == target) {
                    coeffs[idx] = (coeffs[idx] + coeff) % self.p.as_u32();
                }
            }
            module.set_action(*op_degree, op_idx, source_degree, source_idx, &coeffs);
        }

        let mut failures = Vec::new();
        // Deriving decomposable actions and checking relations both use the full Steenrod algebra,
        // which is the wrong algebra for a module carrying a profile. See `Self::restriction`.
        if !self.is_zero() && self.restriction().is_none() {
            let (min, top) = (self.min_degree(), self.top_degree());
            for input_degree in (min..=top).rev() {
                for output_degree in (input_degree + 1)..=top {
                    module.extend_actions(input_degree, output_degree);
                    for failure in module.check_validity_all(input_degree, output_degree) {
                        failures.push(Failure {
                            input_degree,
                            input_idx: failure.input_idx,
                            input_name: self.gens[input_degree][failure.input_idx].clone(),
                            op_degree: output_degree - input_degree,
                            relation: failure.relation,
                            value: failure.value,
                        });
                    }
                }
            }
        }

        (module, failures)
    }

    // -------------------------------------------------------- serialised state

    /// The module as json, in exactly the format used by `ext/steenrod_modules`.
    pub fn to_json(&self) -> Value {
        let (module, _) = self.build();
        let mut json = json!({ "p": self.p.as_u32() });
        module.to_json(&mut json);
        for (key, value) in &self.extra {
            json[key] = value.clone();
        }
        json
    }

    /// Read a module file.
    ///
    /// Unlike [`FDModule::from_json`] this does *not* reject a module that fails the Adem relations:
    /// a file with a mistake in it should be loadable so that the mistake can be fixed. Failures are
    /// reported through [`Self::state`] instead. Syntax errors are still errors.
    pub fn from_json(json: &Value) -> anyhow::Result<Self> {
        let object = json
            .as_object()
            .context("A module file must be a json object")?;

        let p = json["p"]
            .as_u64()
            .context("Module file is missing an integer `p` field")?;
        let p = u32::try_from(p)
            .ok()
            .and_then(|p| ValidPrime::try_from(p).ok())
            .with_context(|| format!("{p} is not a supported prime"))?;

        match json["type"].as_str() {
            Some(MODULE_TYPE) => {}
            Some(other) => bail!(
                "The module builder can only load modules of type \"{MODULE_TYPE}\", but this \
                 file has type \"{other}\""
            ),
            None => bail!("Module file is missing a `type` field"),
        }

        let mut builder = Self::new(p);
        builder.name = json["name"].as_str().unwrap_or_default().to_string();

        let gens = json["gens"]
            .as_object()
            .context("Module file is missing a `gens` object")?;
        for (name, degree) in gens {
            let degree = degree
                .as_i64()
                .with_context(|| format!("The degree of generator {name} is not an integer"))?;
            let degree = i32::try_from(degree)
                .with_context(|| format!("The degree of generator {name} is out of range"))?;
            builder.add_generator(degree, Some(name))?;
        }

        let actions: Vec<String> = match &json["actions"] {
            Value::Null => Vec::new(),
            Value::Array(actions) => actions
                .iter()
                .map(|action| {
                    action
                        .as_str()
                        .map(str::to_string)
                        .context("Every entry of `actions` must be a string")
                })
                .collect::<anyhow::Result<_>>()?,
            _ => bail!("The `actions` field must be an array of strings"),
        };
        builder.set_actions(&actions)?;

        for (key, value) in object {
            if !OWNED_FIELDS.contains(&key.as_str()) {
                builder.extra.insert(key.clone(), value.clone());
            }
        }

        Ok(builder)
    }

    /// Everything the interface needs to draw the module.
    pub fn state(&self) -> Value {
        let (module, failures) = self.build();

        let basis: Vec<Value> = self
            .gens
            .iter_enum()
            .map(|(degree, names)| json!({ "degree": degree, "names": names }))
            .collect();

        let generators: Vec<Value> = self
            .generator_degrees()
            .into_iter()
            .map(|(degree, name)| json!({ "degree": degree, "name": name }))
            .collect();

        // Every non-zero action, of every basis element of the algebra rather than only the
        // generators, so the interface can show what the generator actions imply.
        let mut arcs = Vec::new();
        for (source_degree, names) in self.gens.iter_enum() {
            for op_degree in 1..=(self.top_degree() - source_degree) {
                let output_degree = source_degree + op_degree;
                if self.dimension(output_degree) == 0 {
                    continue;
                }
                self.algebra.compute_basis(op_degree);
                let generator_idx = self.generator_index(op_degree);
                for op_idx in 0..self.algebra.dimension(op_degree) {
                    let is_generator = generator_idx == Some(op_idx);
                    for source_idx in 0..names.len() {
                        let vector = module.action(op_degree, op_idx, source_degree, source_idx);
                        if vector.is_zero() {
                            continue;
                        }
                        let targets: Vec<Value> = vector
                            .iter_nonzero()
                            .map(|(idx, coeff)| json!({ "idx": idx, "coeff": coeff }))
                            .collect();
                        arcs.push(json!({
                            "op_degree": op_degree,
                            "op_name": if is_generator {
                                self.algebra.generator_to_string(op_degree, op_idx)
                            } else {
                                self.algebra.basis_element_to_string(op_degree, op_idx)
                            },
                            "is_generator": is_generator,
                            "source_degree": source_degree,
                            "source_idx": source_idx,
                            "targets": targets,
                        }));
                    }
                }
            }
        }

        let json = self.to_json();
        let actions_text = json["actions"]
            .as_array()
            .map(|actions| {
                actions
                    .iter()
                    .filter_map(Value::as_str)
                    .collect::<Vec<_>>()
                    .join("\n")
            })
            .unwrap_or_default();

        let restriction = self.restriction();
        json!({
            "p": self.p.as_u32(),
            "name": self.name,
            "min_degree": self.min_degree(),
            "top_degree": self.top_degree(),
            "is_zero": self.is_zero(),
            "basis": basis,
            "generators": generators,
            "arcs": arcs,
            "failures": failures.iter().map(Failure::to_json).collect::<Vec<_>>(),
            // `null` rather than `true` when the relations were not checked at all, so the interface
            // can distinguish "no failures" from "not checked".
            "valid": if restriction.is_some() { Value::Null } else { Value::Bool(failures.is_empty()) },
            "restriction": restriction,
            "actions_text": actions_text,
            "json": json,
        })
    }

    /// Evaluate an arbitrary Steenrod expression applied to module elements, e.g. `Sq2*Sq1 x0` or
    /// `Sq2 x0 + Sq1 x1`.
    ///
    /// This is the assertion box: it lets the user ask what an operation that is *not* an algebra
    /// generator does, given the generator actions they have entered.
    pub fn evaluate(&self, expr: &str) -> anyhow::Result<String> {
        let expr = expr.trim();
        if expr.is_empty() {
            return Ok("0".to_string());
        }
        let (module, _) = self.build();

        // (op_degree, op_vector, source_degree, source_idx), with `None` meaning the identity.
        let mut terms = Vec::new();
        let mut output_degree = None;
        for (op_expr, name) in split_terms(expr)? {
            let (source_degree, source_idx) = self
                .lookup(name)
                .with_context(|| format!("Unknown basis element: {name}"))?;
            let operation = if op_expr.is_empty() {
                None
            } else {
                Some(self.evaluator.evaluate_algebra_adem(op_expr)?)
            };
            let op_degree = operation.as_ref().map_or(0, |(degree, _)| *degree);
            let degree = source_degree + op_degree;
            match output_degree {
                None => output_degree = Some(degree),
                Some(previous) => ensure!(
                    previous == degree,
                    "The terms of this expression have different degrees, {previous} and {degree}"
                ),
            }
            terms.push((operation, source_degree, source_idx));
        }
        let output_degree = output_degree.expect("split_terms returns at least one term");

        let mut result = FpVector::new(self.p, module.dimension(output_degree));
        for (operation, source_degree, source_idx) in terms {
            match operation {
                None => result.add_basis_element(source_idx, 1),
                Some((op_degree, op_vector)) => {
                    self.algebra.compute_basis(op_degree);
                    for (op_idx, coeff) in op_vector.iter_nonzero() {
                        module.act_on_basis(
                            result.as_slice_mut(),
                            coeff,
                            op_degree,
                            op_idx,
                            source_degree,
                            source_idx,
                        );
                    }
                }
            }
        }
        Ok(module.element_to_string(output_degree, result.as_slice()))
    }
}

/// Split an expression such as `Sq2*Sq1 x0 + Sq1 x1` into its `(operation, basis element)` terms.
///
/// The Steenrod algebra's own module parser requires an explicit `*` between the operation and the
/// element (`Sq2 * x0`), whereas the `actions` syntax writes them separated by a space
/// (`Sq2 x0 = ...`). Users should not have to remember which box wants which, so we accept either,
/// as well as a bare element name for the identity operation.
fn split_terms(expr: &str) -> anyhow::Result<Vec<(&str, &str)>> {
    expr.split('+')
        .map(|term| {
            let term = term.trim();
            ensure!(!term.is_empty(), "Empty term in expression `{expr}`");
            let split = term
                .rsplit_once(char::is_whitespace)
                .or_else(|| term.rsplit_once('*'));
            match split {
                // A bare name is the element itself, i.e. the identity operation.
                None => Ok(("", term)),
                Some((operation, name)) => {
                    let operation = operation.trim().trim_end_matches('*').trim();
                    let name = name.trim();
                    ensure!(
                        !name.is_empty(),
                        "Expected a basis element after the operation in `{term}`"
                    );
                    Ok((operation, name))
                }
            }
        })
        .collect()
}

/// Reject names that the `actions` syntax could not represent.
///
/// The grammar [`FDModule::parse_action`] implements splits an action on `" = "`, separates the
/// operation from the element at the last space, and splits the right-hand side on `" + "` and then
/// on a space to peel off a coefficient. So a name may not contain whitespace, `+` or `=`, and must
/// not begin with a digit or it would be read as a coefficient.
///
/// Anything else is allowed. In particular `*` is permitted, because the library uses
/// tensor-product names such as `x0*x0` in `ext/steenrod_modules/C2_sm_Ceta.json`, and a name we
/// reject is a file we could not load.
fn validate_name(name: &str) -> anyhow::Result<()> {
    let first = name
        .chars()
        .next()
        .context("A basis element name cannot be empty")?;
    ensure!(
        first.is_alphabetic(),
        "A basis element name must start with a letter, but {name} does not"
    );
    for c in name.chars() {
        ensure!(
            !c.is_whitespace(),
            "A basis element name cannot contain whitespace: {name}"
        );
        ensure!(
            c != '+' && c != '=',
            "A basis element name cannot contain {c:?}, which separates the parts of an action: \
             {name}"
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p2() -> ValidPrime {
        ValidPrime::new(2)
    }

    /// The Joker, which is the standard test module: it has a non-trivial `Sq2` and needs
    /// `extend_actions` to fill in `Sq3 = Sq1 Sq2`.
    fn joker() -> Value {
        json!({
            "p": 2,
            "type": "finite dimensional module",
            "gens": { "x0": 0, "x1": 1, "x2": 2, "x3": 3, "x4": 4 },
            "actions": [
                "Sq1 x0 = x1",
                "Sq2 x1 = x3",
                "Sq1 x3 = x4",
                "Sq2 x0 = x2",
                "Sq2 x2 = x4",
            ],
        })
    }

    #[test]
    fn empty_builder_is_zero() {
        let builder = Builder::new(p2());
        assert!(builder.is_zero());
        let state = builder.state();
        assert_eq!(state["valid"], json!(true));
        assert_eq!(state["json"]["gens"], json!({}));
        assert_eq!(state["json"]["actions"], json!([]));
    }

    #[test]
    fn build_c2_by_hand() {
        let mut builder = Builder::new(p2());
        assert_eq!(builder.add_generator(0, None).unwrap(), "x0");
        assert_eq!(builder.add_generator(1, None).unwrap(), "x1");
        builder.set_action(1, 0, 0, &[1]).unwrap();

        let state = builder.state();
        assert_eq!(state["valid"], json!(true));
        assert_eq!(state["json"]["actions"], json!(["Sq1 x0 = x1"]));
        assert_eq!(state["generators"][0]["name"], json!("Sq1"));
    }

    /// Round-tripping the Joker through json must be stable, and the resulting module must equal the
    /// one the library builds from the same file.
    #[test]
    fn joker_round_trip() {
        let builder = Builder::from_json(&joker()).unwrap();
        let state = builder.state();
        assert_eq!(state["valid"], json!(true));

        let again = Builder::from_json(&builder.to_json()).unwrap();
        assert_eq!(builder.to_json(), again.to_json());

        let algebra = Arc::new(AdemAlgebra::new(p2(), false));
        algebra.compute_basis(10);
        let expected = FDModule::from_json(algebra, &joker()).unwrap();
        let (actual, failures) = builder.build();
        assert!(failures.is_empty());
        actual.test_equal(&expected).unwrap();
    }

    /// `Sq3` is not a generator, so it cannot be set directly, but it is derived: on the Joker
    /// `Sq3 x1 = Sq1 Sq2 x1 = Sq1 x3 = x4`.
    #[test]
    fn derived_actions_are_reported() {
        let mut builder = Builder::from_json(&joker()).unwrap();
        assert!(builder.set_action(3, 0, 0, &[1]).is_err());

        // Both the `actions` spelling and the calculator's `*` spelling are accepted.
        assert_eq!(builder.evaluate("Sq3 x1").unwrap(), "x4");
        assert_eq!(builder.evaluate("Sq3 * x1").unwrap(), "x4");
        assert_eq!(builder.evaluate("Sq3*x1").unwrap(), "x4");
        assert_eq!(builder.evaluate("Sq1*Sq2 x1").unwrap(), "x4");
        assert_eq!(builder.evaluate("Sq3 x0").unwrap(), "0");
        // A bare name is the element itself, and sums of terms work.
        assert_eq!(builder.evaluate("x1").unwrap(), "x1");
        assert_eq!(builder.evaluate("Sq2 x2 + Sq1 x3").unwrap(), "0");
        assert!(builder.evaluate("Sq1 x0 + Sq1 x1").is_err());
        assert!(builder.evaluate("Sq1 nonsense").is_err());

        let arcs = builder.state()["arcs"].as_array().unwrap().clone();
        let sq3 = arcs
            .iter()
            .find(|arc| arc["op_name"] == json!("Sq3") && arc["source_degree"] == json!(1))
            .expect("Sq3 acting on x1 should be reported");
        assert_eq!(sq3["is_generator"], json!(false));
    }

    /// `Sq1 Sq1 = 0` is the simplest Adem relation, and violating it must be reported against the
    /// basis element responsible rather than just failing.
    #[test]
    fn sq1_squared_is_reported() {
        let mut builder = Builder::new(p2());
        builder.add_generator(0, None).unwrap();
        builder.add_generator(1, None).unwrap();
        builder.add_generator(2, None).unwrap();
        builder.set_action(1, 0, 0, &[1]).unwrap();
        builder.set_action(1, 1, 0, &[1]).unwrap();

        let (_, failures) = builder.build();
        assert_eq!(failures.len(), 1);
        assert_eq!(failures[0].relation, "1 * Sq1 * Sq1");
        assert_eq!(failures[0].input_name, "x0");
        assert_eq!(failures[0].op_degree, 2);
        assert_eq!(failures[0].value, "x2");

        assert_eq!(builder.state()["valid"], json!(false));
    }

    /// A file that fails the Adem relations must still load, so that it can be repaired. This is the
    /// one place we deliberately differ from `FDModule::from_json`.
    #[test]
    fn invalid_module_still_loads() {
        let json = json!({
            "p": 2,
            "type": "finite dimensional module",
            "gens": { "x0": 0, "x1": 1, "x2": 2 },
            "actions": ["Sq1 x0 = x1", "Sq1 x1 = x2"],
        });
        assert!(FDModule::from_json(Arc::new(AdemAlgebra::new(p2(), false)), &json).is_err());

        let builder = Builder::from_json(&json).unwrap();
        assert_eq!(builder.state()["valid"], json!(false));
    }

    #[test]
    fn removing_a_generator_drops_its_actions() {
        let mut builder = Builder::from_json(&joker()).unwrap();
        // Remove x3, which is both a target (of Sq2 x1) and a source (of Sq1 x3).
        let (degree, idx) = builder.lookup("x3").unwrap();
        builder.remove_generator(degree, idx).unwrap();

        let json = builder.to_json();
        let actions = json["actions"].as_array().unwrap();
        assert!(actions.iter().all(|a| !a.as_str().unwrap().contains("x3")));
        assert_eq!(json["gens"].as_object().unwrap().len(), 4);
    }

    #[test]
    fn renaming_rewrites_actions() {
        let mut builder = Builder::from_json(&joker()).unwrap();
        let (degree, idx) = builder.lookup("x1").unwrap();
        builder.rename_generator(degree, idx, "y").unwrap();

        let json = builder.to_json();
        assert!(json["gens"].get("y").is_some());
        let actions = json["actions"].as_array().unwrap();
        assert!(actions.iter().any(|a| a.as_str().unwrap() == "Sq1 x0 = y"));
        assert!(actions.iter().any(|a| a.as_str().unwrap() == "Sq2 y = x3"));
    }

    #[test]
    fn rejects_duplicate_and_malformed_names() {
        let mut builder = Builder::new(p2());
        builder.add_generator(0, Some("x0")).unwrap();
        assert!(builder.add_generator(1, Some("x0")).is_err());
        assert!(builder.add_generator(1, Some("1x")).is_err());
        assert!(builder.add_generator(1, Some("x 0")).is_err());
        assert!(builder.add_generator(1, Some("x+y")).is_err());
        assert!(builder.add_generator(1, Some("")).is_err());
        // Tensor-product names are used by the library and must be accepted.
        assert!(builder.add_generator(1, Some("x0*y0")).is_ok());
    }

    /// Removing the only element of the top degree must shrink the module rather than leave an empty
    /// degree hanging off the end.
    #[test]
    fn removal_trims_empty_degrees() {
        let mut builder = Builder::new(p2());
        builder.add_generator(0, None).unwrap();
        builder.add_generator(3, None).unwrap();
        assert_eq!(builder.top_degree(), 3);

        let (degree, idx) = builder.lookup("x3").unwrap();
        builder.remove_generator(degree, idx).unwrap();
        assert_eq!(builder.min_degree(), 0);
        assert_eq!(builder.top_degree(), 0);
    }

    #[test]
    fn negative_degrees_work() {
        let mut builder = Builder::new(p2());
        assert_eq!(builder.add_generator(-2, None).unwrap(), "x_2");
        builder.add_generator(-1, None).unwrap();
        builder.set_action(1, -2, 0, &[1]).unwrap();

        let state = builder.state();
        assert_eq!(state["valid"], json!(true));
        assert_eq!(state["min_degree"], json!(-2));
        assert_eq!(state["json"]["actions"], json!(["Sq1 x_2 = x_1"]));
    }

    #[test]
    fn drag_toggles_an_arc() {
        let mut builder = Builder::new(p2());
        builder.add_generator(0, None).unwrap();
        builder.add_generator(1, None).unwrap();

        builder.add_to_action(1, 0, 0, 0, 1).unwrap();
        assert_eq!(builder.action_coefficients(1, 0, 0).unwrap(), vec![1]);
        // At the prime 2 dragging the same arc again removes it.
        builder.add_to_action(1, 0, 0, 0, 1).unwrap();
        assert_eq!(builder.action_coefficients(1, 0, 0).unwrap(), vec![0]);
        assert_eq!(builder.to_json()["actions"], json!([]));
    }

    #[test]
    fn actions_text_round_trips() {
        let mut builder = Builder::from_json(&joker()).unwrap();
        let text = builder.state()["actions_text"]
            .as_str()
            .unwrap()
            .to_string();
        builder.set_actions_text(&text).unwrap();
        assert_eq!(
            builder.to_json(),
            Builder::from_json(&joker()).unwrap().to_json()
        );

        assert!(builder.set_actions_text("Sq1 x0 = nonsense").is_err());
        assert!(builder.set_actions_text("Sq1 x0").is_err());
        // A degree mismatch is caught by the library parser.
        assert!(builder.set_actions_text("Sq1 x0 = x3").is_err());
    }

    #[test]
    fn preserves_unknown_fields() {
        let json = json!({
            "p": 2,
            "type": "finite dimensional module",
            "gens": { "x0": 0, "x1": 1 },
            "actions": ["Sq1 x0 = x1"],
            "cofiber": { "idx": 0, "s": 4, "t": 12 },
            "self_maps": [
                { "hom_deg": 4, "int_deg": 12, "map_data": [[1]], "name": "v_1^4" }
            ],
        });
        let builder = Builder::from_json(&json).unwrap();
        let out = builder.to_json();
        assert_eq!(out["cofiber"], json["cofiber"]);
        assert_eq!(out["self_maps"], json["self_maps"]);
    }

    /// A module over a sub-Hopf-algebra must load and save without loss, but must not be reported as
    /// failing relations that do not apply to it. This is `ext/steenrod_modules/tmf2_sm_DA1.json`,
    /// which is an $A(2)$-module.
    #[test]
    fn profiled_module_is_not_checked() {
        let json = json!({
            "p": 2,
            "algebra": ["milnor"],
            "profile": { "truncated": true, "p_part": [3, 2, 1] },
            "type": "finite dimensional module",
            "gens": { "x0": 0, "x2": 2, "x4": 4, "x8": 8, "x10": 10, "x12": 12 },
            "actions": [
                "Sq2 x0 = x2",
                "Sq4 x0 = x4",
                "Sq4 x4 = x8",
                "Sq4 x8 = x12",
                "Sq2 x10 = x12",
            ],
        });
        let builder = Builder::from_json(&json).unwrap();
        assert!(builder.restriction().is_some());

        let state = builder.state();
        // Not checked, rather than valid or invalid.
        assert_eq!(state["valid"], Value::Null);
        assert!(state["restriction"].is_string());
        assert_eq!(state["failures"], json!([]));

        // The profile survives, and the actions still round-trip.
        let out = builder.to_json();
        assert_eq!(out["profile"], json["profile"]);
        assert_eq!(out["actions"], json["actions"]);
    }

    #[test]
    fn rejects_other_module_types() {
        let json = json!({ "p": 2, "type": "real projective space", "min": 1 });
        assert!(Builder::from_json(&json).is_err());
    }

    // ---------------------------------------------------------- odd primes

    /// `C(5, v_1)` at the prime 5, which is `ext/steenrod_modules/C5v1.json`. The generators at odd
    /// primes are `b` and the `P^{p^k}`.
    #[test]
    fn odd_prime_module() {
        let json = json!({
            "p": 5,
            "type": "finite dimensional module",
            "gens": { "x0": 0, "x1": 1, "x9": 9, "x10": 10 },
            "actions": ["b x0 = x1", "P1 x1 = x9", "b x9 = x10"],
        });
        let builder = Builder::from_json(&json).unwrap();
        let state = builder.state();
        assert_eq!(state["valid"], json!(true));
        assert_eq!(
            state["json"]["actions"],
            json!(["b x0 = x1", "P1 x1 = x9", "b x9 = x10"])
        );

        // Degree 1 is the Bockstein and degree 8 = q is P^1; degree 2 carries no generator.
        let generators = state["generators"].as_array().unwrap();
        assert_eq!(generators[0], json!({ "degree": 1, "name": "b" }));
        assert!(
            generators
                .iter()
                .any(|g| g == &json!({ "degree": 8, "name": "P1" }))
        );
        assert!(!generators.iter().any(|g| g["degree"] == json!(2)));
    }

    /// `b^2 = 0` is the odd-prime edge case in `generating_relations`.
    #[test]
    fn odd_prime_bockstein_squared_is_reported() {
        let mut builder = Builder::new(ValidPrime::new(3));
        builder.add_generator(0, None).unwrap();
        builder.add_generator(1, None).unwrap();
        builder.add_generator(2, None).unwrap();
        builder.set_action(1, 0, 0, &[1]).unwrap();
        builder.set_action(1, 1, 0, &[1]).unwrap();

        let (_, failures) = builder.build();
        assert!(!failures.is_empty());
        assert_eq!(failures[0].input_name, "x0");
        assert!(failures[0].relation.contains('b'));
    }

    /// Coefficients other than 1 have to survive a round trip at odd primes.
    #[test]
    fn odd_prime_coefficients() {
        let mut builder = Builder::new(ValidPrime::new(3));
        builder.add_generator(0, None).unwrap();
        builder.add_generator(1, None).unwrap();
        builder.set_action(1, 0, 0, &[2]).unwrap();

        assert_eq!(builder.to_json()["actions"], json!(["b x0 = 2 x1"]));
        let again = Builder::from_json(&builder.to_json()).unwrap();
        assert_eq!(again.action_coefficients(1, 0, 0).unwrap(), vec![2]);

        // Coefficients are reduced mod p, and p times anything is no action at all.
        builder.set_action(1, 0, 0, &[3]).unwrap();
        assert_eq!(builder.to_json()["actions"], json!([]));
    }

    #[test]
    fn changing_prime_clears_actions() {
        let mut builder = Builder::from_json(&joker()).unwrap();
        builder.set_prime(ValidPrime::new(3));
        assert_eq!(builder.prime(), ValidPrime::new(3));
        assert_eq!(builder.to_json()["actions"], json!([]));
        assert_eq!(builder.to_json()["p"], json!(3));
        // The basis elements survive.
        assert_eq!(builder.to_json()["gens"].as_object().unwrap().len(), 5);
    }
}
