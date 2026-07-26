//! Check the builder against every finite dimensional module in `ext/steenrod_modules`.
//!
//! These files are the format the builder has to be compatible with, so they are the natural test
//! corpus: whatever is in the library must load, round-trip, and produce the same module the library
//! itself produces.

use std::{path::PathBuf, sync::Arc};

use algebra::{AdemAlgebra, Algebra, module::FDModule};
use fp::prime::ValidPrime;
use module_builder::builder::Builder;
use serde_json::Value;

fn module_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../ext/steenrod_modules")
}

/// Every `finite dimensional module` in the library, as `(file name, json)`.
fn library() -> Vec<(String, Value)> {
    let mut modules = Vec::new();
    for entry in std::fs::read_dir(module_dir()).expect("the module library should be readable") {
        let path = entry.expect("directory entries should be readable").path();
        if path.extension().is_none_or(|ext| ext != "json") {
            continue;
        }
        let contents = std::fs::read_to_string(&path).expect("module files should be readable");
        let json: Value = serde_json::from_str(&contents)
            .unwrap_or_else(|e| panic!("{} is not valid json: {e}", path.display()));
        if json["type"].as_str() != Some("finite dimensional module") {
            continue;
        }
        let name = path
            .file_name()
            .expect("a file we just read has a name")
            .to_string_lossy()
            .into_owned();
        modules.push((name, json));
    }
    assert!(
        modules.len() > 10,
        "expected to find the module library, found {} files",
        modules.len()
    );
    modules.sort_by(|a, b| a.0.cmp(&b.0));
    modules
}

/// Loading and saving must be idempotent: the second save has to equal the first.
///
/// The first save need not equal the file on disk, since the file may write its actions in a
/// different but equivalent order from `FDModule::to_json`'s canonical one.
#[test]
fn round_trip() {
    for (name, json) in library() {
        let builder = Builder::from_json(&json).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let once = builder.to_json();
        let again = Builder::from_json(&once)
            .unwrap_or_else(|e| panic!("{name} does not reload after saving: {e:?}"));
        assert_eq!(once, again.to_json(), "{name} does not round-trip");
    }
}

/// The module the builder produces must be the module `ext` produces from the same file.
///
/// Modules carrying a `profile` are excluded: they are modules over a sub-Hopf-algebra, so `ext`
/// itself would not load them with the full Adem algebra either. They are covered by
/// [`profiled_modules_are_not_checked`] instead.
#[test]
fn agrees_with_fdmodule() {
    let mut compared = 0;
    for (name, json) in library() {
        let builder = Builder::from_json(&json).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        if builder.restriction().is_some() || builder.is_zero() {
            continue;
        }

        let (built, failures) = builder.build();
        assert!(
            failures.is_empty(),
            "{name} is in the library but fails the Adem relations: {failures:?}"
        );

        let p = ValidPrime::try_from(json["p"].as_u64().unwrap() as u32).unwrap();
        let algebra = Arc::new(AdemAlgebra::new(p, false));
        algebra.compute_basis(builder.top_degree() - builder.min_degree());
        let expected = FDModule::from_json(algebra, &json)
            .unwrap_or_else(|e| panic!("{name} does not load in ext: {e:?}"));
        built
            .test_equal(&expected)
            .unwrap_or_else(|e| panic!("{name} disagrees with ext's module: {e}"));
        compared += 1;
    }
    assert!(compared > 10, "only compared {compared} modules");
}

/// The library's sub-Hopf-algebra modules — `tmf2`, `ko`, the `y(n)` — must load and be reported as
/// unchecked rather than as failing the Adem relations, which do not apply to them.
#[test]
fn profiled_modules_are_not_checked() {
    let mut seen = 0;
    for (name, json) in library() {
        if json["profile"].is_null() {
            continue;
        }
        let builder = Builder::from_json(&json).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert!(
            builder.restriction().is_some(),
            "{name} has a profile but is not flagged as restricted"
        );
        let state = builder.state();
        assert!(state["valid"].is_null(), "{name} should not be checked");
        assert_eq!(state["failures"], serde_json::json!([]), "{name}");
        seen += 1;
    }
    assert!(seen > 0, "no library module carries a profile");
}

/// Saving must preserve the fields the builder does not interpret, since several library modules
/// carry `cofiber`, `products`, `self_maps` or `profile` and would otherwise be silently downgraded.
#[test]
fn preserves_extra_fields() {
    const EXTRA: [&str; 5] = ["cofiber", "products", "self_maps", "profile", "algebra"];
    let mut seen = 0;
    for (name, json) in library() {
        let saved = Builder::from_json(&json)
            .unwrap_or_else(|e| panic!("{name}: {e:?}"))
            .to_json();
        for field in EXTRA {
            if !json[field].is_null() {
                seen += 1;
                assert_eq!(saved[field], json[field], "{name} lost its {field} field");
            }
        }
    }
    assert!(seen > 0, "no library module exercised the preserved fields");
}
