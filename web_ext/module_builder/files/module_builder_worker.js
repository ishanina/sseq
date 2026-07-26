'use strict';
importScripts('./module_builder_wasm.js');

const { ModuleBuilder, check_module } = wasm_bindgen;

const ready = wasm_bindgen('./module_builder_wasm_bg.wasm').catch(e => {
    self.postMessage({
        id: null,
        ok: false,
        error: `Failed to load wasm: ${e}`,
    });
    throw e;
});

/// The single module being edited. Every command that changes it returns the new state, so the main
/// thread never has to ask for it separately.
let builder = null;

function requireBuilder() {
    if (builder === null) {
        throw new Error('No module is open');
    }
    return builder;
}

/// Commands return either `{ state }` or `{ value }`. A command that mutates returns the new state,
/// which is what `ModuleBuilder`'s methods hand back.
const commands = {
    create: ({ p }) => {
        const created = ModuleBuilder.new(p);
        if (created === undefined || created === null) {
            throw new Error(`${p} is not a supported prime`);
        }
        builder = created;
        return { state: builder.state() };
    },
    load: ({ json }) => {
        // Only replace the open module once the new one has parsed, so a failed load leaves the
        // editor as it was.
        builder = ModuleBuilder.load(json);
        return { state: builder.state() };
    },
    state: () => ({ state: requireBuilder().state() }),
    setName: ({ name }) => ({ state: requireBuilder().set_name(name) }),
    setPrime: ({ p }) => ({ state: requireBuilder().set_prime(p) }),
    addGenerator: ({ degree, name }) => ({
        state: requireBuilder().add_generator(degree, name ?? ''),
    }),
    removeGenerator: ({ degree, idx }) => ({
        state: requireBuilder().remove_generator(degree, idx),
    }),
    renameGenerator: ({ degree, idx, name }) => ({
        state: requireBuilder().rename_generator(degree, idx, name),
    }),
    setAction: ({ opDegree, sourceDegree, sourceIdx, coeffs }) => ({
        state: requireBuilder().set_action(
            opDegree,
            sourceDegree,
            sourceIdx,
            new Uint32Array(coeffs),
        ),
    }),
    addToAction: ({ opDegree, sourceDegree, sourceIdx, targetIdx, coeff }) => ({
        state: requireBuilder().add_to_action(
            opDegree,
            sourceDegree,
            sourceIdx,
            targetIdx,
            coeff,
        ),
    }),
    setActionsText: ({ text }) => ({
        state: requireBuilder().set_actions_text(text),
    }),
    shift: ({ by }) => ({ state: requireBuilder().shift(by) }),
    dual: () => ({ state: requireBuilder().dual() }),
    truncate: ({ min, max }) => ({
        // `undefined` rather than `null` is what wasm-bindgen maps to `Option::None`.
        state: requireBuilder().truncate(min ?? undefined, max ?? undefined),
    }),
    tensor: ({ other }) => ({ state: requireBuilder().tensor(other) }),
    directSum: ({ other }) => ({ state: requireBuilder().direct_sum(other) }),
    submodule: ({ cells }) => ({
        state: requireBuilder().submodule(new Int32Array(cells)),
    }),
    quotient: ({ cells }) => ({
        state: requireBuilder().quotient(new Int32Array(cells)),
    }),
    evaluate: ({ expr }) => ({ value: requireBuilder().evaluate(expr) }),
    toJson: () => ({ value: requireBuilder().to_json() }),
    toJsonCompact: () => ({ value: requireBuilder().to_json_compact() }),
    checkModule: ({ json }) => ({ value: check_module(json) }),
};

self.onmessage = async ev => {
    await ready;
    const { id, cmd, args } = ev.data;
    const handler = commands[cmd];
    if (handler === undefined) {
        self.postMessage({ id, ok: false, error: `Unknown command: ${cmd}` });
        return;
    }
    try {
        const result = handler(args ?? {});
        self.postMessage({ id, ok: true, ...result });
    } catch (e) {
        // Errors thrown from wasm are plain strings; errors thrown here are `Error`s.
        self.postMessage({
            id,
            ok: false,
            error: e instanceof Error ? e.message : `${e}`,
        });
    }
};
