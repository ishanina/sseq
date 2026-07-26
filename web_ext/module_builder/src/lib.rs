//! A web interface for building modules over the Steenrod algebra.
//!
//! The interesting code is in [`builder`]; this module is the `wasm_bindgen` boundary. Every method
//! takes and returns json encoded as a string, which keeps the boundary trivial and lets the logic
//! be tested without a browser.

use algebra::module::FDModule;
use fp::prime::{Prime, ValidPrime};
use wasm_bindgen::prelude::*;

pub mod builder;

use builder::Builder;

/// Convert an [`anyhow::Error`] into something JS can throw, including the whole context chain
/// rather than only the outermost message.
fn to_js(error: anyhow::Error) -> JsValue {
    let mut message = error.to_string();
    for cause in error.chain().skip(1) {
        message.push_str(&format!("\n  caused by: {cause}"));
    }
    JsValue::from(message)
}

fn parse(json: &str) -> Result<serde_json::Value, JsValue> {
    serde_json::from_str(json).map_err(|e| JsValue::from(format!("Invalid json: {e}")))
}

#[wasm_bindgen]
pub struct ModuleBuilder(Builder);

#[wasm_bindgen]
impl ModuleBuilder {
    /// Create an empty module over the Steenrod algebra at the prime `p`.
    ///
    /// Returns `None` if `p` is not a prime the library supports.
    pub fn new(p: u32) -> Option<Self> {
        Some(Self(Builder::new(ValidPrime::try_from(p).ok()?)))
    }

    /// Load a module file. See [`Builder::from_json`] for what is accepted.
    pub fn load(json: &str) -> Result<ModuleBuilder, JsValue> {
        Builder::from_json(&parse(json)?).map(Self).map_err(to_js)
    }

    /// Everything the interface needs to draw the module, as json.
    pub fn state(&self) -> String {
        self.0.state().to_string()
    }

    /// The module in the format used by `ext/steenrod_modules`, pretty-printed for saving to a file.
    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(&self.0.to_json()).expect("a module always serialises")
    }

    /// The module as compact json, for putting in a URL.
    pub fn to_json_compact(&self) -> String {
        self.0.to_json().to_string()
    }

    pub fn prime(&self) -> u32 {
        self.0.prime().as_u32()
    }

    /// Change the prime. This keeps the basis elements but discards the actions, since which
    /// operations are algebra generators depends on the prime.
    pub fn set_prime(&mut self, p: u32) -> Result<String, JsValue> {
        let p = ValidPrime::try_from(p).map_err(|e| JsValue::from(e.to_string()))?;
        self.0.set_prime(p);
        Ok(self.state())
    }

    pub fn set_name(&mut self, name: &str) -> String {
        self.0.set_name(name);
        self.state()
    }

    /// Add a basis element in `degree`. An empty `name` means "choose one for me".
    pub fn add_generator(&mut self, degree: i32, name: &str) -> Result<String, JsValue> {
        let name = if name.trim().is_empty() {
            None
        } else {
            Some(name)
        };
        self.0.add_generator(degree, name).map_err(to_js)?;
        Ok(self.state())
    }

    pub fn remove_generator(&mut self, degree: i32, idx: usize) -> Result<String, JsValue> {
        self.0.remove_generator(degree, idx).map_err(to_js)?;
        Ok(self.state())
    }

    pub fn rename_generator(
        &mut self,
        degree: i32,
        idx: usize,
        name: &str,
    ) -> Result<String, JsValue> {
        self.0.rename_generator(degree, idx, name).map_err(to_js)?;
        Ok(self.state())
    }

    /// Set the action of the algebra generator in `op_degree` on a basis element, as a coefficient
    /// vector for the basis of `source_degree + op_degree`.
    pub fn set_action(
        &mut self,
        op_degree: i32,
        source_degree: i32,
        source_idx: usize,
        coeffs: Vec<u32>,
    ) -> Result<String, JsValue> {
        self.0
            .set_action(op_degree, source_degree, source_idx, &coeffs)
            .map_err(to_js)?;
        Ok(self.state())
    }

    /// Add `coeff` to a single entry of an action, reducing mod `p`. This is what dragging an arc
    /// between two cells does.
    pub fn add_to_action(
        &mut self,
        op_degree: i32,
        source_degree: i32,
        source_idx: usize,
        target_idx: usize,
        coeff: u32,
    ) -> Result<String, JsValue> {
        self.0
            .add_to_action(op_degree, source_degree, source_idx, target_idx, coeff)
            .map_err(to_js)?;
        Ok(self.state())
    }

    /// Replace every action with the ones parsed from `text`, one per line.
    pub fn set_actions_text(&mut self, text: &str) -> Result<String, JsValue> {
        self.0.set_actions_text(text).map_err(to_js)?;
        Ok(self.state())
    }

    /// Evaluate a Steenrod expression applied to module elements, e.g. `Sq2*Sq1 x0`.
    pub fn evaluate(&self, expr: &str) -> Result<String, JsValue> {
        self.0.evaluate(expr).map_err(to_js)
    }
}

/// Check a module file without loading it into the editor, for the library listing.
///
/// Returns the empty string if the module is valid, and a description of the first failure
/// otherwise. Unlike [`ModuleBuilder::load`] this uses [`FDModule::from_json`], so it reports exactly
/// what `ext` would say about the file.
#[wasm_bindgen]
pub fn check_module(json: &str) -> Result<String, JsValue> {
    use std::sync::Arc;

    use algebra::AdemAlgebra;

    let json = parse(json)?;
    let p = json["p"]
        .as_u64()
        .and_then(|p| u32::try_from(p).ok())
        .and_then(|p| ValidPrime::try_from(p).ok())
        .ok_or_else(|| JsValue::from("Module file is missing a valid `p` field"))?;
    let algebra = Arc::new(AdemAlgebra::new(p, false));
    match FDModule::from_json(algebra, &json) {
        Ok(_) => Ok(String::new()),
        Err(e) => Ok(e.to_string()),
    }
}
