# Steenrod Module Builder

An interactive web interface for building finite dimensional modules over the
Steenrod algebra. Add cells, drag the action of the algebra generators between
them, and see the Adem relations checked as you go. Modules save in the same
JSON format as `ext/steenrod_modules`, so the result can be resolved by `ext` or
opened directly in the Adams spectral sequence viewer.

Works at every prime the library supports: at the prime 2 the generators you can
draw are the $Sq^{2^k}$, and at an odd prime they are the Bockstein $\beta$ and
the $P^{p^k}$.

`DESIGN.md` records the plan this is being built to, including the phases that
are not implemented yet — the module operations (shift, dual, tensor, submodule,
quotient) and the derived-module and Ext features.

## Build and run

Requires [`just`](https://github.com/casey/just). To install the matching
`wasm-bindgen-cli` and the wasm target:

```console
 $ just setup-wasm
```

Then build the site into `dist/` and serve it:

```console
 $ just serve
```

`just all` builds without serving. `wasm-opt` is used to shrink the binary if
`binaryen` is installed, and skipped otherwise.

## Tests

```console
 $ just test
```

The unit tests in `src/builder.rs` cover editing, the relation checker and odd
primes. `tests/steenrod_modules.rs` runs the whole of `ext/steenrod_modules`
through the builder and checks that every finite dimensional module there loads,
round-trips, and produces exactly the module `FDModule::from_json` produces.

## Notes on what is editable

The free parameters of a Steenrod module are the actions of the algebra
*generators*; every other operation is determined by them, and non-linearly, so
there is nothing to solve for. `Sq3` cannot be drawn because it is `Sq1 Sq2`.
The interface therefore only lets you drag arcs across generator degrees, and
offers the derived actions read-only: turn on *show derived actions*, or ask the
evaluator for a specific one.

A module whose file carries a `profile` field, such as `tmf2` or `ko`, is a
module over a sub-Hopf-algebra rather than over the whole Steenrod algebra. Such
a file loads, displays and saves with its profile intact, but the relation
checker for the full algebra does not apply to it and is switched off. Editing
over a proper subalgebra needs the profile-aware Milnor generating set and is
not implemented.
