# `module_builder` — design

An interactive web interface for building modules over the Steenrod algebra, and
more generally bounded chain complexes of such ("derived modules"), producing
files that are drop-in compatible with `ext/steenrod_modules/`.

This document is the plan. Phases 1 and 2 are implemented, and phase 3a — the
$\Ext^1$ extension enumeration — with it; the rest of phase 3 and phase 4 are
not. See the [Phasing](#phasing) section at the end for what that covers, and
`README.md` for how to build and run what exists.

## Goals

1. Add and remove cells (basis elements) in arbitrary internal degrees.
2. Draw the action of the Steenrod algebra generators between cells — the
   $Sq^{2^k}$ at $p = 2$, and $\beta$ and $P^{p^k}$ at odd primes.
3. Continuously check that the Adem relations hold, and say precisely which
   relation fails and where when they do not.
4. Save and load, in exactly the JSON format the rest of the repo uses.
5. Operations: shift, dual, tensor, direct sum, submodule, quotient,
   truncation, and extensions.
6. Handle bounded chain complexes of finite-dimensional modules, and take the
   cofibre or fibre of a class in $\Ext^i(X, Y)$ between two such.
7. Work at every prime the library supports, not only at 2.
8. Hand a finished module straight to `sseq_gui` to compute its Adams
   spectral sequence.

## What the repository already provides

Most of the mathematics exists; the work is a UI, a serialisation format for
complexes, and three or four genuinely missing pieces of algebra.

| Capability | Location |
| --- | --- |
| The module type | `FDModule` = `algebra::module::FiniteDimensionalModule`, `ext/crates/algebra/src/module/finite_dimensional_module.rs` |
| Adem relation checking | `FDModule::check_validity`, which walks `AdemAlgebra::generating_relations` (`adem_algebra.rs:384`) — literally the inadmissible-pair Adem relations |
| Filling in decomposable operations | `FDModule::extend_actions`, via `GeneratedAlgebra::decompose_basis_element` |
| Load/save in the repo's format | `FDModule::from_json` / `to_json` |
| Text action syntax `Sq1 x0 = x1` | `FDModule::parse_action` |
| Module equality (for tests) | `FDModule::test_equal`, plus `PartialEq` |
| Shift | `algebra::module::SuspensionModule` |
| Tensor product | `algebra::module::TensorModule`, plus `FDModule::from(&_)`; see `ext/examples/tensor.rs` |
| Quotients | `algebra::module::QuotientModule` |
| Steenrod algebra parsing/evaluation, Adem <-> Milnor | `algebra::steenrod_evaluator` (what `steenrod_calculator` already wraps) |
| Bounded complexes of modules | `ext::chain_complex::FiniteChainComplex`, `FiniteAugmentedChainComplex` |
| Explicit maps of modules from matrices | `FullModuleHomomorphism::from_matrices` |
| Realising an Ext class as an iterated extension | `ext::yoneda::{yoneda_representative, yoneda_representative_element, try_yoneda_representative_element}` — returns a `FiniteAugmentedChainComplex` of `FDModule`s |
| Computing `Ext(M, N)` for two finite modules | `HomModule` + `HomPullback` + `Subquotient`, assembled in `ext/examples/ext_m_n.rs` |
| An existing interactive builder, CLI only | `ext/examples/define_module.rs` |
| The deployment pattern to copy | `web_ext/steenrod_calculator` (crate + `files/` + `justfile` + a CI job) |

Two convenient hooks:

- `sseq_gui` accepts `?module_json=<urlencoded json>`
  (`web_ext/sseq_gui/interface/index.js:286`), so we can hand off a module with
  no file round-trip.
- `web_ext/sseq_gui/justfile` already copies `ext/steenrod_modules` into `dist/`,
  so serving the module library from the site is a one-line recipe.

## Derived modules

The interesting question is whether the tool can handle objects of the derived
category rather than just modules, so that arbitrary Ext classes are
expressible. It can, and to a large extent the repository already does this.

`ext/src/utils.rs:219-282` reads a `cofiber` field from a module's JSON:

```json
{
    "p": 2,
    "type": "finite dimensional module",
    "gens": { "x0": 0, "x1": 1 },
    "actions": ["Sq1 x0 = x1"],
    "cofiber": { "idx": 0, "s": 4, "t": 12 }
}
```

That is `ext/steenrod_modules/C2v14.json`. Reading it resolves `C(2)`, picks out
the generator of $\Ext^{4,12}$ (which is $v_1^4$), builds the corresponding
`ChainMap`, calls `yoneda::yoneda_representative` to get a
`FiniteAugmentedChainComplex` of `FDModule`s, pops the augmentation, and resolves
*that* complex. `C4`, `C9`, `Ceta2`, `C3v1b1` are defined the same way.

So the object the repository actually resolves is `CCC = FiniteChainComplex<SteenrodModule>`
(`ext/src/lib.rs:181`), a plain module being the special case
`FiniteChainComplex::ccdz(module)`. That is the right native object for the
workspace:

```
DerivedModule = FiniteChainComplex<FDModule<SteenrodAlgebra>,
                                   FullModuleHomomorphism<FDModule<SteenrodAlgebra>>>
```

Consequences:

- A module is a complex concentrated in homological degree 0, and saves as
  exactly today's `{gens, actions}` format. No compatibility break.
- Shift means shift in either grading (internal degree, or homological degree).
- Extensions are cones. For $\alpha \in \Ext^1_A(Q, N) = \operatorname{Hom}_{D(A)}(Q, N[1])$,
  the extension $0 \to N \to M \to Q \to 0$ is the fibre of $\alpha$.
- $\Ext^{s}$ classes give $s$-fold extensions, which `yoneda_representative_element`
  already constructs out of `FDModule`s.

Choosing this object model from the start avoids reworking the data model when
the Ext features land in phase 3. Note that the `cofiber` field is currently
undocumented in the module-specification section of `ext/src/lib.rs`; part of
this work is to document it.

## File format

Compatibility is the hard requirement: anything saved must be loadable by
`resolve`, `sseq_gui`, and the rest of `ext`.

**Single modules.** Byte-for-byte the current format — `FDModule::to_json`
output, with the optional `name`, `algebra`, `profile`, `products`, `self_maps`,
`shift`, `cofiber` fields preserved across a load/save cycle. `products` and
`self_maps` drive `sseq_gui`'s display, so they get a metadata editor rather
than being silently dropped.

**Cofibres.** The builder can emit the existing `cofiber` field, which keeps the
file loadable by today's `ext` unchanged. It can also expand a `cofiber` spec
into an explicit complex for display and editing (by running the same
`yoneda_representative` call `construct_standard` does).

**General complexes.** These need a new format. Proposal:

```json
{
    "p": 2,
    "type": "chain complex of finite dimensional modules",
    "modules": [ { "gens": {...}, "actions": [...] }, ... ],
    "differentials": [
        { "degree_shift": 0, "matrices": { "4": [[1, 0], [0, 1]] } }
    ]
}
```

`modules[s]` is the module in homological degree `s`; `differentials[s]` is the
map `modules[s+1] -> modules[s]`. Following
`FullModuleHomomorphism::apply_to_basis_element`, `matrices` is keyed by
*target* internal degree `t`, with rows indexed by the source basis in degree
`t + degree_shift` and columns by the target basis in degree `t`.

Two rules keep this honest:

- Writing a complex concentrated in homological degree 0 must produce exactly
  the single-module format, and reading a single-module file must produce
  `ccdz` of it. Enforced by a round-trip test.
- A file using the new `type` is *not* loadable by today's `ext`. Loading it
  needs a corresponding arm in `steenrod_module::from_json` /
  `construct_standard`; that is a small additive change and is listed below. The
  UI will warn when an object can only be saved in the new format, and offer the
  `cofiber` encoding instead when the object happens to be a cofibre.

## Backend: the wasm crate

Layout mirrors `web_ext/steenrod_calculator` exactly: `Cargo.toml`
(`crate-type = ["cdylib", "rlib"]`, deps on `algebra`, `bivec`, `fp`, and — from
phase 3 — `ext`), `src/lib.rs` holding the `#[wasm_bindgen]` surface, static
assets in `files/`, and a `justfile` cloned from the calculator's
(`setup-wasm`, `wasm-lib`, `wasm-bindgen-step`, `all`).

The algebra is `SteenrodAlgebra` rather than bare `AdemAlgebra`, so that
`TensorModule` (which needs `Bialgebra`) works and so that Milnor-basis display
is available. This is what `ext/examples/tensor.rs` does.

### The source of truth is the edit list, not the `FDModule`

State is

```rust
struct ModuleSpec {
    p: ValidPrime,
    gens: Vec<(String, i32)>,
    actions: Vec<(i32, usize, i32, usize, Vec<u32>)>, // op_deg, op_idx, in_deg, in_idx, coeffs
}
```

and every edit rebuilds an `FDModule` from scratch: `FDModule::new` ->
`set_basis_element_name` -> `set_action` -> `extend_actions` -> `check_validity`,
in the same order as `FDModule::from_json`.

Reasons:

- `extend_actions` writes into the table it reads from, so incremental mutation
  after an edit is hard to reason about; a rebuild is exactly reproducible.
- Undo/redo and "revert this arc" become trivial.
- It avoids `FDModule::add_generator`, which cannot extend below the current
  `min_degree`: it calls `graded_dimension.extend_with` but never
  `extend_negative` (`finite_dimensional_module.rs:266`). Duals put generators
  in negative degrees and would hit this immediately. Either fix that method or
  keep avoiding it; the rebuild path means we are not blocked on the decision.
- At the sizes involved (tens of cells) a rebuild is microseconds.

### wasm surface

`rebuild() -> Json` returning graded dimensions, names, the full derived action
table (every Adem basis element, not just generators), and the list of relation
failures; `add_generator`, `remove_generator`, `rename_generator`, `set_action`,
`parse_actions_text`, `to_json`, `from_json`; `eval_action(expr, element)` for
querying the action of an arbitrary Steenrod expression via
`SteenrodEvaluator`; and one entry point per operation.

The wasm runs in a worker, using the same `no-modules` + `importScripts` pattern
as `files/steenrod_calculator_worker.js`, so a panic or a long Ext computation
cannot freeze the page.

## Editing generator actions

The editable data is the action of the algebra *generators*: at $p = 2$ the
$Sq^{2^k}$ (`AdemAlgebra::generators` returns empty unless
`degree.count_ones() == 1`, `adem_algebra.rs:357`), and at odd primes $\beta$ and
$P^{p^k}$. Everything else is forced by `extend_actions`. This is also exactly
what the JSON `actions` array stores, and what
`ext/examples/define_module.rs` prompts for.

So arcs are draggable only across generator degree gaps. Dragging across any
other gap explains why rather than silently failing: "Sq3 is not a generator; it
is determined by Sq2 Sq1."

Two read-only conveniences on top:

1. A derived view of the action of every $Sq^i$ and every Adem basis element, so
   one can see what a choice of generator actions implies.
2. An assertion box where any operation may be typed — `Sq3 x0 = x3`, or
   `Sq2*Sq1 x0 = x3` — and is checked against the derived value, reporting
   agreement or the actual forced value. This reuses `SteenrodEvaluator`, the
   same engine as `steenrod_calculator`.

## Odd primes

Odd primes are a first-class requirement, not an afterthought, and the library
supports them throughout. Concretely:

- **Feature flags.** `Cargo.toml` must use `default-features = false, features =
  ["odd-primes"]` on `fp` and `algebra`, exactly as
  `web_ext/steenrod_calculator/Cargo.toml` does. The prime selector offers
  2, 3, 5, 7, as the calculator's does.
- **Generators and syntax.** At odd primes the generators are $\beta$ in degree 1
  and $P^{p^k}$ in degree $q p^k$, written `b` and `P{n}` by
  `AdemAlgebra::generator_to_string`. So the action strings look like
  `["b x0 = x1", "P1 x1 = x9", "b x9 = x10"]` — that is
  `ext/steenrod_modules/C5v1.json` verbatim. Coefficients in $\F_p$ appear on
  both the arcs and the right-hand sides, e.g. `P1 x0 = 2 x4`.
- **Relations.** `AdemAlgebra::generating_relations` returns $\beta^2 = 0$ as an
  explicit edge case in degree 2, then $P^i P^j$ for $i < pj$ and $P^i \beta P^j$
  for $i < pj + 1$ (`adem_algebra.rs:382-410`). `check_validity` therefore checks
  exactly the right thing at odd primes with no extra work.
- **Signs.** `TensorModule` already carries the Koszul signs
  (`fp::prime::minus_one_to_the_n`), and the antipode implementation will carry
  them too rather than assuming $\beta$ is the only odd-degree factor.
- **Resolutions.** Nassau's algorithm is $p = 2$ only and `construct_nassau`
  rejects odd primes outright (`ext/src/utils.rs:182`), so the phase 3 Ext
  features must go through `construct_standard`, which is the default path.
- **Coverage.** Odd-prime resolutions are already benchmarked in CI
  (`resolve-C3`, `resolve-S_3`), and `C9.json` / `C3v1b1.json` use `cofiber` at
  $p = 3$. But every `yoneda-*` benchmark is at $p = 2$, so the Yoneda path at
  odd primes is supported in principle and untested in practice. Adding an
  odd-prime Yoneda test is part of this work, not an optional extra.

### Choice of basis

The generating set, and hence which arcs are editable, depends on which basis the
algebra is instantiated with. `AdemAlgebra::generators` returns only the
$Sq^{2^k}$ (resp. $\beta$, $P^{p^k}$), whereas `MilnorAlgebra::generators`
additionally returns the $Q_k$ and the $P^s_t$, and is profile-aware
(`milnor_algebra.rs:696`).

The builder therefore instantiates the **Adem** basis for editing: it is the
minimal generating set, it matches `ext/examples/define_module.rs`, and its
`generating_relations` are the Adem relations the user wants checked. Crucially,
the generator names it emits — `Sq{n}`, `P{n}`, `b` — are also accepted by
`MilnorAlgebra::basis_element_from_string` (`milnor_algebra.rs:569`), which is why
a file like `Calpha.json` containing `"P1 x0 = x4"` loads correctly even though
the default basis for resolving is Milnor. Sticking to those forms is a hard
constraint on what we write.

One consequence: a `profile` field is honoured only by the Milnor branch of
`SteenrodAlgebra::from_json` and silently ignored by the Adem branch. Modules
carrying a profile (e.g. `y1_3.json`) will be loaded, displayed, and re-saved
with the field preserved, but *editing* over a proper sub-Hopf-algebra requires
the Milnor generating set and is deferred; the UI will say so rather than
pretending the profile is in effect.

## Operations

Each produces a new workspace entry; inputs are never mutated.

| Operation | Route |
| --- | --- |
| Shift (internal) | `SuspensionModule` -> `FDModule::from` |
| Shift (homological) | Re-index the complex |
| Direct sum | Direct; this is also the split extension |
| Tensor | `TensorModule` -> `FDModule::from`, as in `ext/examples/tensor.rs`; for complexes, totalise the double complex |
| Dual | New — see below |
| Submodule generated by selected cells | Close the span under the action using `fp::matrix::Subspace` |
| Quotient by that submodule | Complementary basis, emitted as a fresh `FDModule` with readable names |
| Truncate above/below a degree | `QuotientModule`, or directly |
| Cone / fibre of a map | `FiniteChainComplex` bookkeeping over `FullModuleHomomorphism` |
| Cohomology of a complex | Kernel/image via `ModuleHomomorphism::{kernel, image}` |
| Extensions | See below |

### Dual

Nothing in the repository dualises a module (`grep -ri dual` finds only
comments). The construction is $(M^\vee)_{-n} = (M_n)^*$ with

$$\langle a f, x \rangle = \langle f, \chi(a) x \rangle,$$

where $\chi$ is the antipode. The naive transpose is **wrong**: transposing
gives a module over $A^{\mathrm{op}}$, and $A$ is not commutative; $\chi$ is
what identifies the two.

The antipode is not in the crate either, and is worth adding on its own merits.
It is cheap to compute:

- For a single square, $\chi(Sq^n) = -\sum_{j<n} \chi(Sq^j) Sq^{n-j}$, which
  needs only `Algebra::multiply_basis_elements`.
- For a general basis element, use `Bialgebra::decompose` to write it as a
  product of squares (Adem basis) and multiply the $\chi$'s in reverse order,
  with the Koszul sign $(-1)^{|x||y|}$. At odd primes the only odd-degree factor
  is $\beta$ and $\beta^2 = 0$, but the signs will be implemented in general
  rather than assumed away.

Tests: the defining identity $\sum \chi(x')x'' = \varepsilon$ and $\chi^2 =
\mathrm{id}$, checked exhaustively in low degrees at $p = 2, 3, 5$ (the full
coproduct of a basis element being assembled from `decompose` + `coproduct` the
same way `TensorModule::act_helper` does it); the classical values
$\chi(Sq^3) = Sq^2 Sq^1$ and $\chi(Sq^4) = Sq^4 + Sq^3 Sq^1$; and
$M^{\vee\vee} = M$.

Two corrections to what this section originally said, both found while
implementing it:

- **The Joker is not self-dual here.** It is self-dual as an $A(1)$-module, which
  is the familiar statement, but this is the full Steenrod algebra. Over $A$ the
  library's Joker has $Sq^4 x_0 = 0$, while its dual has $Sq^4 x_4^* = x_0^*$,
  precisely because $\chi(Sq^4) = Sq^4 + Sq^3 Sq^1$ and $Sq^3 Sq^1 x_0 = x_4$.
  Both degrees are one dimensional, so this is basis independent: the dual really
  is not isomorphic to the Joker shifted. It makes a sharper test than
  self-duality would have — a naive transpose gives $Sq^4 = 0$ and fails it.
- **The Koszul signs are not only in the antipode.** The coproduct of a product
  needs them too, since multiplication in $A \otimes A$ satisfies
  $(a \otimes b)(c \otimes d) = (-1)^{|b||c|} ac \otimes bd$. The test that checks
  the defining identity has to carry that sign, and `b P1 b` at $p = 3$ and
  $p = 5$ is where it shows.

One limitation, which is the algebra's rather than the antipode's:
`MilnorAlgebra`'s `Bialgebra::coproduct` asserts $p = 2$, so $\chi$ is
unavailable in the Milnor basis at odd primes. This does not affect the builder,
which edits in the Adem basis.

Once the dual's action table is filled in for *all* basis elements, rather than
just generators, `check_validity` should pass automatically — a useful
self-check, and `to_json` will still emit only the generator actions.

### Extensions and arbitrary Ext

Two routes, deliberately in this order.

**Phase 3a — $\Ext^1$ by linear algebra, no resolutions.** Given a submodule $N$
and quotient $Q$, an extension is $N \oplus Q$ as a graded vector space with the
known actions of $N$ and of $Q$ installed, plus an unknown off-diagonal
$\theta \colon Q \to N$ for each generator. Because $\theta$ lands in $N$ and $N$
maps to $0$ in $Q$, any product of two $\theta$'s vanishes, so **the Adem
relation constraints on $\theta$ are linear**. Hence:

- valid $\theta$ = kernel of an explicit $\F_p$-linear map,
- coboundaries = the $\theta$'s obtained by changing the splitting,
- $\Ext^1_A(Q, N)$ = the quotient, computed with `fp::matrix::{Matrix, Subspace, Subquotient}`.

This enumerates all extensions exactly, produces explicit `FDModule`s, needs no
dependency on `ext`, and is fast. It also gives a nice UI: a list of extension
classes, each previewable as a cell diagram. As a fallback it doubles as a manual
mode, where the user fills in the off-diagonal arcs by hand with live checking —
the check being precisely the cocycle condition.

**Phase 3b — arbitrary $\Ext^{i}(X, Y)$ between derived modules, and the
cofibre/fibre of a class.** This is the general goal: given two bounded complexes
$X$ and $Y$ of finite-dimensional modules and a class $\alpha \in
\operatorname{Hom}_{D(A)}(X, Y[i])$, produce the cofibre and the fibre as explicit
complexes.

The existing API reaches further towards this than it first appears:

- **The source may already be a complex.** `MuResolution::new_with_save` takes a
  `CCC = FiniteChainComplex<SteenrodModule>`, and `construct_standard` resolves a
  multi-module complex — the popped Yoneda complex — as a matter of course
  (`ext/src/utils.rs:281-284`). So resolving $X$ needs nothing new.
- **The target may already be a complex.** `yoneda_representative` takes
  `map: ChainMap<FreeModuleHomomorphism<_>>`, and `ChainMap::chain_maps` is a
  `Vec` indexed by `s - s_shift` (`yoneda.rs:248, 258, 375, 472`). The
  single-module case that `construct_standard` uses is `chain_maps: vec![map]`;
  a chain map into a bounded complex $Y$ is the same call with a longer vector.
- **The construction of the cofibre is already spelled out.**
  `yoneda_representative` returns a finite quasi-isomorphic quotient of the
  resolution through which the map factors; `construct_standard` then does
  `FiniteChainComplex::from(yoneda)` followed by `.pop()`, and that is the
  cofibre. The same recipe generalises verbatim.
- The fibre is the cofibre shifted, so one operation plus a homological shift
  covers both.

What must be written:

- **$\Ext^{i}(X, Y)$ for complexes.** The cochain complex in
  `ext/examples/ext_m_n.rs` takes a `FreeChainComplex` source and a single
  *module* target. For a complex target it becomes the Hom double complex,
  totalised. That machinery should be lifted out of the example's private
  `mod hom_cochain_complex` into `ext` proper and generalised there, so the
  example and this crate share it.
- **Assembling the `ChainMap` from a chosen class.** Given a class in the
  computed $\Ext^i$, build the corresponding `FreeModuleHomomorphism`s degree by
  degree — the same `add_generators_from_matrix_rows` / `extend_by_zero` pattern
  as `ext/src/utils.rs:253-271`, one per homological degree of $Y$.
- A UI for picking $\alpha$: a clickable $\Ext^{i}(X, Y)$ chart, which in the
  case $Y = \F_p$ is the ordinary Adams $E_2$ chart.
- Saving: an explicit complex in the new format, or — when the object is a
  cofibre of a class in $\Ext^{s,t}(M, M)$ — the existing `cofiber` field, which
  keeps the file loadable by today's `ext`.

Risks, stated plainly: `ext::yoneda` is one of the heavier parts of the
codebase, uses `std::any::Any` downcasting for its operation-rating heuristic,
and can be slow — hence the `try_` variant and the existing
`ext/tests/try_yoneda.rs`. Its output is not small either: the $g$-class
representative at $(20, 4)$ has modules of dimension $7, 15, 12, 4, 1$ in the
Adem basis (`ext/examples/benchmarks/yoneda-g-S_2-adem`). The UI must run it in
the worker behind an explicit "compute" button with a bounded range, surface
errors rather than panicking, and allow cancellation. Depending on `ext` also
grows the wasm binary substantially; `sseq_gui` already does exactly this and is
deployed, so it is proven to work, but the size should be tracked in CI the way
the other two sites' are.

## Additions to the library crates

All additive; nothing changes an existing signature.

1. `FDModule::check_validity_all(in_deg, out_deg) -> Vec<RelationFailure>` next
   to `check_validity`, where `RelationFailure` records the input index as well
   as the relation and value. The current method returns only the first failure
   and does not say *which* basis element failed, which is what the UI needs in
   order to highlight the offending arcs.
2. `antipode` for the Steenrod algebra, in `algebra::algebra`, with Koszul signs
   (see above). Useful well beyond this tool.
3. A dual for `FDModule`.
4. Lift `hom_cochain_complex` out of `ext/examples/ext_m_n.rs` into `ext`, and
   generalise it from a single module target to a bounded complex target (the Hom
   double complex, totalised) so that $\Ext^i(X, Y)$ between derived modules can
   be computed.
5. A `"chain complex of finite dimensional modules"` arm in
   `steenrod_module::from_json` and `construct_standard`, so files the builder
   writes are loadable by `resolve`.
6. Document the `cofiber` field in the module-specification section of
   `ext/src/lib.rs`.
7. An odd-prime Yoneda test, since CI currently exercises that path only at
   $p = 2$.
8. Optionally, fix `FDModule::add_generator` for degrees below `min_degree`.

## Frontend

Static HTML/CSS/JS in `files/`, no build step, KaTeX and Bootstrap from CDN as
in the calculator, formatted to the repository's `.prettierrc` (four spaces,
single quotes).

- **Centre — cell diagram, in SVG.** Internal degree on the vertical axis, one
  dot per basis element, curved arcs for the generator actions: the standard
  picture, as for the Joker. Click an empty degree to add a cell; drag dot to
  dot to toggle the target in the corresponding generator action; coefficient
  labels at odd primes; arcs styled by which generator they represent.
- **Right — relation status.** A live indicator, and on failure a list of the
  violated Adem relations, each naming the offending basis element and the
  nonzero value it produces, with the implicated arcs highlighted in the
  diagram.
- **Text mode**, two-way synced with the diagram, showing the raw `actions`
  array. Parsed with `FDModule::parse_action`, so what is displayed is exactly
  what will be saved.
- **Left — workspace and library.** Several objects open at once, which is
  needed for tensor and extension operations; the ~45 files from
  `ext/steenrod_modules/` copied into `dist/` as a library; file load, save,
  drag-and-drop, and paste; `localStorage` autosave; and a `?module_json=`
  permalink.
- **Complex view.** For objects with more than one homological degree, a column
  per degree with the differentials shown between them, and the ability to drill
  into any single module.
- **"Compute Ext"** links to `../?module_json=...`, opening the object in
  `sseq_gui`, which is where the deployed site puts it. That is not where it is
  when the builder is served on its own, so the page checks that a viewer is
  really there and takes `?viewer=` for one served elsewhere.

## Build, CI, deployment

`just all` builds the wasm and copies `files/` plus `ext/steenrod_modules` into
`dist/`. A `module_builder` job in `.github/workflows/ext.yaml`, cloned from the
existing `calculator` job (lint, `setup-wasm`, `just all`, wasm size benchmark,
artifact upload), added to `deploy`'s `needs` list and downloaded into a
`module_builder/` path, so it publishes to
`https://spectralsequences.github.io/sseq/module_builder/`.

Links from the root `README.md`, from the calculator, and from the "Custom"
section of `sseq_gui`'s `index.html`.

## Testing

- Round-trip every file in `ext/steenrod_modules/`: `from_json -> to_json ->
  from_json`, compared with `FDModule::test_equal`.
- Round-trip a one-term complex through the new format and check it comes back
  as the plain single-module format.
- Operations: dual (including double dual and the Joker's self-duality), tensor
  cross-checked against `ext/examples/tensor.rs`, shift, submodule and quotient
  of the Joker and of `A-mod-Sq1-Sq2`.
- Negative cases: `Sq1 x0 = x1, Sq1 x1 = x2` must be rejected, naming the
  `Sq1 * Sq1` relation and the element `x0`.
- Extensions: the $\Ext^1$ enumeration must reproduce, for instance, $C(4)$ and
  $C(\eta^2)$, which the library also defines via `cofiber` — a genuine
  cross-check of the two routes against each other.
- Odd primes throughout, not only at $p = 2$: round-trip and relation checks at
  $p = 3, 5, 7$; $\beta^2 = 0$ as a rejection case; the dual and the antipode at
  odd primes; and a Yoneda/cofibre test at $p = 3$, which CI does not currently
  cover at all. `C9.json` and `C3v1b1.json` are the natural fixtures, since both
  are defined by a `cofiber` at $p = 3$.
- `node --test` for the pure-JS parts (state reducer, arc geometry), in the style
  of `web_ext/sseq_gui/wasm/worker_panic.test.mjs`.
- Optionally a selenium test in the style of `web_ext/sseq_gui/tests`.

## Phasing

1. **Core — done.** Crate, wasm API, `check_validity_all`, cell-diagram editor,
   text mode, live Adem checking, save/load, module library, `sseq_gui` hand-off,
   CI job and deployment.
2. **Module operations — done.** Shift, dual, truncate, direct sum, tensor,
   submodule, quotient, plus the antipode they rest on. Multi-cell selection and
   an undo stack came with them, since the operations are destructive.
3. **Derived modules.** In three parts:
   - **a — done.** $\Ext^1$ enumeration by linear algebra, and building the
     extension it picks out.
   - **b.** Complexes as the native object, cone and fibre, cohomology, and the
     new file format.
   - **c.** Resolution-backed $\Ext^i(X, Y)$ charts and the cofibre or fibre of a
     chosen class, via Yoneda.
4. **Polish.** Remaining tests, further documentation.

Phases 1 and 2 are self-contained and useful on their own, and neither depends
on `ext`. Phase 3 is where the wasm binary grows and where the schedule risk
lives.

### What phase 1 changed, and what it taught us

The CI job and deployment were pulled forward out of phase 4, since a site that
is not deployed is not testable by anyone else.

Two things came out of building it that the plan above did not anticipate:

- **Sub-Hopf-algebra modules are not modules over the Steenrod algebra.** Thirteen
  of the thirty-six finite dimensional modules in `ext/steenrod_modules` carry a
  `profile`; `tmf2_sm_DA1.json` is an $A(2)$-module on which $Sq^8$ does not act,
  so checking it against the full Adem relations reports three failures that are
  not mistakes. The builder now detects a `profile`, declines to check the module,
  and says why. The plan had already deferred *editing* over a subalgebra, but not
  noticed that *checking* had to be handled too.
- **Generator names may contain more than `[A-Za-z0-9_]`.** The library uses
  tensor-product names such as `x0*x0` in `C2_sm_Ceta.json`. Validation now
  rejects only what the `actions` grammar genuinely cannot represent — whitespace,
  `+`, `=`, and a leading digit — because a name we reject is a file we cannot
  load.

One incidental fix outside the crate: `sseq_gui`'s home page rewrote the `href`
of every anchor in a module section into a `?module=` link, so it could not
contain an ordinary link. It now skips anchors with no `data` attribute.

### What phase 2 changed

The antipode corrections are recorded under [Dual](#dual) above. Three further
decisions worth writing down:

- **Operations drop metadata.** `cofiber`, `products` and `self_maps` describe the
  module they were written for. The cofibre spec of $M$ says nothing about
  $M^\vee$ or $M \otimes N$, so carrying those fields across an operation would
  attach a false claim to the result. Every operation clears them. `shift` and
  `truncate` are the exception in spirit — they too clear, since a `self_maps`
  entry's internal degree no longer matches after a suspension.
- **Operations refuse a profiled module.** Everything here works over the full
  Steenrod algebra, so dualising or tensoring an $A(2)$-module would silently give
  the wrong answer. Those operations now fail with an explanation instead. `shift`
  and `truncate` are allowed, since neither consults the algebra.
- **Suspension needs no rebuild.** Because actions are keyed by generator name
  rather than by degree and index, shifting is just a re-grading of the basis; the
  action list is untouched. This fell out of the phase 1 data model rather than
  being designed for.

The destructive nature of the operations also forced two interface changes that
were not in the original plan: selection became a *set* of cells, since submodule
and quotient act on a set, and an undo stack was added, since a mis-clicked
`dualise` would otherwise lose work.

### What phase 3a changed

The $\Ext^1$ computation went in as `src/extension.rs`, working on `FDModule`s
rather than on `Builder`s so that it needs none of the builder's internals.

The linearity argument in the plan held up, and the implementation leans on it
twice: $F$ is evaluated one elementary $\theta$ at a time to get its matrix, and
the split extension is used as the check that $F(0) = 0$, i.e. that $F$ really is
linear rather than merely affine. Relation *values* rather than relation
*failures* were needed, so there is a small vector-valued twin of
`check_validity` in that file; the string-producing one in the algebra crate
cannot serve, since the linear algebra needs coordinates.

One interface gap only showed up in the browser: a class needs the submodule to
sit above the quotient, and the operand picker offered library modules
unshifted, so every pair a user could reach had $\Ext^1 = 0$. The operand now
takes a shift. Suspending a file is just adding to each degree in `gens`, since
actions name their generators, so this is done in the page rather than in the
wasm — though it does have to drop `cofiber`, `self_maps` and `products`, which
are stated in the unshifted degrees.

The result reproduces the $h_i$: extending the sphere by its own suspension
$\Sigma^t S$ gives $\Ext^1$ of dimension 1 exactly when $t$ is a generator
degree — 1, 2, 4, 8 and not 3, 5 or 6 — and the non-split extensions in degrees
1 and 2 are $C(2)$ and $C(\eta)$ on the nose.
