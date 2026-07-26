'use strict';

// ---------------------------------------------------------------- worker RPC

const worker = new Worker('./module_builder_worker.js');
const pending = new Map();
let nextId = 0;

worker.addEventListener('message', ev => {
    const { id, ok, error, state, value } = ev.data;
    // A message with no id is a failure to start up at all, which nobody is waiting on.
    if (id === null || id === undefined) {
        showToast(error, 'error');
        setStatus('loading', 'failed');
        return;
    }
    const request = pending.get(id);
    if (request === undefined) {
        return;
    }
    pending.delete(id);
    if (ok) {
        request.resolve(state === undefined ? value : JSON.parse(state));
    } else {
        request.reject(new Error(error));
    }
});

function call(cmd, args) {
    const id = nextId++;
    worker.postMessage({ id, cmd, args });
    return new Promise((resolve, reject) =>
        pending.set(id, { resolve, reject }),
    );
}

/// Run a command that returns a new module state, and render it.
///
/// Errors are shown to the user rather than thrown: every one of them is a mistake in what was
/// asked for (a bad name, an action in a degree that carries no generator), not a bug.
///
/// The previous module is pushed onto the undo stack first, and only if the command succeeds, so a
/// rejected operation does not leave a useless entry behind.
async function run(cmd, args) {
    const previous = state === null ? null : JSON.stringify(state.json);
    try {
        const next = await call(cmd, args);
        if (previous !== null) {
            undoStack.push(previous);
            if (undoStack.length > UNDO_LIMIT) {
                undoStack.shift();
            }
        }
        setState(next);
        return true;
    } catch (e) {
        showToast(e.message, 'error');
        return false;
    }
}

async function undo() {
    const previous = undoStack.pop();
    if (previous === undefined) {
        return;
    }
    try {
        setState(await call('load', { json: previous }));
    } catch (e) {
        showToast(`Could not undo: ${e.message}`, 'error');
    }
}

// ------------------------------------------------------------- module state

/// The last state received from the worker. The worker owns the module; this is a snapshot for
/// rendering.
let state = null;

/// The range of degrees the diagram shows. It always covers the module, but the user can extend it
/// in order to add cells outside the current range.
let view = { min: 0, max: 4 };

/// The selected cells, as a set of `"degree,idx"` keys. Several may be selected at once, since the
/// submodule and quotient operations act on a set of cells.
let selected = new Set();

/// Previous module states, most recent last, so that a destructive operation can be undone.
const undoStack = [];
const UNDO_LIMIT = 64;

function cellKey({ degree, idx }) {
    return `${degree},${idx}`;
}

function selectedCells() {
    return [...selected].map(key => {
        const [degree, idx] = key.split(',');
        return {
            degree: Number.parseInt(degree, 10),
            idx: Number.parseInt(idx, 10),
        };
    });
}

function setState(next) {
    state = next;
    if (!state.is_zero) {
        view.min = Math.min(view.min, state.min_degree);
        view.max = Math.max(view.max, state.top_degree);
    }
    // Drop selections that no longer point at a cell, which happens after any operation that
    // reshapes the module.
    for (const cell of selectedCells()) {
        if (dimension(cell.degree) <= cell.idx) {
            selected.delete(cellKey(cell));
        }
    }
    save();
    render();
}

function dimension(degree) {
    const row = state.basis.find(d => d.degree === degree);
    return row === undefined ? 0 : row.names.length;
}

function names(degree) {
    const row = state.basis.find(d => d.degree === degree);
    return row === undefined ? [] : row.names;
}

/// The algebra generator whose degree is `gap`, or undefined if there is none.
function generatorAt(gap) {
    return state.generators.find(g => g.degree === gap);
}

// ------------------------------------------------------------------ rendering

const SVG_NS = 'http://www.w3.org/2000/svg';
const ROW_HEIGHT = 52;
/// Room on the left for the degree numbers.
const GUTTER = 46;
const CELL_X0 = GUTTER + 40;
/// Enough for a name like `x0*x0` to sit above its cell without touching its neighbour's.
const CELL_DX = 84;
const CELL_R = 5.5;
/// Cell names go above their cell, leaving the whole space to the right of the column free for arcs,
/// which is the only way a single-file module like the Joker stays legible.
const LABEL_DY = -13;

function el(name, attrs = {}, text = null) {
    const node = document.createElementNS(SVG_NS, name);
    for (const [key, value] of Object.entries(attrs)) {
        node.setAttribute(key, value);
    }
    if (text !== null) {
        node.textContent = text;
    }
    return node;
}

function y(degree) {
    return (view.max - degree) * ROW_HEIGHT + ROW_HEIGHT / 2;
}

function x(idx) {
    return CELL_X0 + idx * CELL_DX;
}

/// Split an operation name into text and superscript parts, so `Sq1` renders as Sq¹, `b` as β and
/// `Q1` as Q₁. Products in the Adem basis are written with spaces, e.g. `Sq2 Sq1`.
function opLabelParts(name) {
    const parts = [];
    const re = /(Sq|P|Q)(\d+)|(b)/g;
    let last = 0;
    let match;
    while ((match = re.exec(name)) !== null) {
        if (match.index > last) {
            parts.push({ text: name.slice(last, match.index) });
        }
        if (match[3] !== undefined) {
            parts.push({ text: 'β' });
        } else {
            parts.push({ text: match[1] });
            parts.push({ text: match[2], shift: match[1] === 'Q' ? 3 : -4 });
        }
        last = re.lastIndex;
    }
    if (last < name.length) {
        parts.push({ text: name.slice(last) });
    }
    return parts;
}

function opLabel(name, attrs) {
    const text = el('text', attrs);
    // `dy` on a tspan is cumulative, so each shift has to be undone by the next part.
    let shift = 0;
    for (const part of opLabelParts(name)) {
        const span = el('tspan', {}, part.text);
        const wanted = part.shift ?? 0;
        if (wanted !== shift) {
            span.setAttribute('dy', wanted - shift);
            shift = wanted;
        }
        if (wanted !== 0) {
            span.setAttribute('font-size', '78%');
        }
        text.appendChild(span);
    }
    return text;
}

function render() {
    const svg = document.getElementById('diagram');
    svg.textContent = '';

    const rows = view.max - view.min + 1;
    const widest = Math.max(1, ...state.basis.map(d => d.names.length));
    // Leave room for the widest arc, which is the one spanning the whole module.
    const span = state.is_zero ? 1 : state.top_degree - state.min_degree;
    const width = Math.max(
        x(widest) + arcBulge(span, widest - 1, false) + 70,
        380,
    );
    const height = rows * ROW_HEIGHT;
    svg.setAttribute('width', width);
    svg.setAttribute('height', height);
    svg.setAttribute('viewBox', `0 0 ${width} ${height}`);

    const implicated = new Set(
        state.failures.map(f => `${f.input_degree},${f.input_idx}`),
    );
    const showDerived = document.getElementById('show-derived').checked;

    // Degree rows: the label, a gridline, and a click target for adding a cell.
    for (let degree = view.min; degree <= view.max; degree++) {
        const top = (view.max - degree) * ROW_HEIGHT;
        const row = el('rect', {
            class: 'degree-row',
            x: GUTTER,
            y: top,
            width: width - GUTTER,
            height: ROW_HEIGHT,
        });
        row.dataset.degree = degree;
        svg.appendChild(row);
        svg.appendChild(
            el('line', {
                class: 'gridline',
                x1: GUTTER,
                y1: y(degree),
                x2: width,
                y2: y(degree),
            }),
        );
        svg.appendChild(
            el(
                'text',
                { class: 'degree-label', x: 8, y: y(degree) },
                `${degree}`,
            ),
        );
    }

    for (const arc of state.arcs) {
        if (!arc.is_generator && !showDerived) {
            continue;
        }
        for (const target of arc.targets) {
            drawArc(svg, arc, target);
        }
    }

    for (const { degree, names: row } of state.basis) {
        row.forEach((name, idx) => {
            const classes = ['cell'];
            if (selected.has(cellKey({ degree, idx }))) {
                classes.push('selected');
            }
            if (implicated.has(`${degree},${idx}`)) {
                classes.push('implicated');
            }
            const cell = el('circle', {
                class: classes.join(' '),
                cx: x(idx),
                cy: y(degree),
                r: CELL_R,
            });
            cell.dataset.degree = degree;
            cell.dataset.idx = idx;
            svg.appendChild(cell);
            svg.appendChild(
                el(
                    'text',
                    { class: 'cell-label', x: x(idx), y: y(degree) + LABEL_DY },
                    name,
                ),
            );
        });
    }

    renderPanel();
}

/// How far to the right the control point of an arc sits.
///
/// Without this every arc of a single-file module like the Joker would lie along the same vertical
/// line. Arcs are separated by how many generators are below the operation rather than by the degree
/// itself: the generators are Sq1, Sq2, Sq4, Sq8, … so that ordinal stays small even for a module
/// twenty degrees wide, whereas the degree does not. A quadratic curve reaches only half way to its
/// control point, so these numbers are twice the visible offset.
function arcBulge(opDegree, indexGap, isGenerator) {
    const rank = state.generators.filter(g => g.degree < opDegree).length;
    // Nudge derived arcs outwards so that, say, Sq3 and Sq4 do not coincide.
    return 56 + 40 * (rank + (isGenerator ? 0 : 0.5)) + 26 * indexGap;
}

function drawArc(svg, arc, target) {
    const sourceDegree = arc.source_degree;
    const targetDegree = sourceDegree + arc.op_degree;
    const x1 = x(arc.source_idx);
    const y1 = y(sourceDegree);
    const x2 = x(target.idx);
    const y2 = y(targetDegree);
    const bulge = arcBulge(
        arc.op_degree,
        Math.abs(target.idx - arc.source_idx),
        arc.is_generator,
    );
    const cx = Math.max(x1, x2) + bulge;
    const cy = (y1 + y2) / 2;
    const kind = arc.is_generator ? 'arc' : 'arc derived';
    svg.appendChild(
        el('path', {
            class: kind,
            d: `M ${x1} ${y1} Q ${cx} ${cy} ${x2} ${y2}`,
        }),
    );

    // The midpoint of a quadratic Bézier, where the label goes. It sits just outside the apex, so a
    // label never lands on the curve it names or on the column of cells.
    const mx = 0.25 * x1 + 0.5 * cx + 0.25 * x2;
    const my = 0.25 * y1 + 0.5 * cy + 0.25 * y2;
    const label =
        target.coeff === 1 ? arc.op_name : `${target.coeff}·${arc.op_name}`;
    svg.appendChild(
        opLabel(label, {
            class: arc.is_generator ? 'arc-label' : 'arc-label derived',
            x: mx + 5,
            y: my,
        }),
    );
}

function renderPanel() {
    if (state.valid === null) {
        setStatus('unchecked', 'not checked');
    } else if (state.valid) {
        setStatus('valid', 'relations hold');
    } else {
        const n = state.failures.length;
        setStatus('invalid', `${n} relation${n === 1 ? '' : 's'} fail`);
    }

    const relations = document.getElementById('relations');
    relations.textContent = '';
    if (state.restriction !== null) {
        const notice = document.createElement('p');
        notice.className = 'notice';
        notice.textContent = state.restriction;
        relations.appendChild(notice);
    } else if (state.failures.length === 0) {
        const ok = document.createElement('p');
        ok.className = 'hint';
        ok.textContent = state.is_zero
            ? 'The zero module satisfies every relation. Add a cell to begin.'
            : 'Every Adem relation is satisfied.';
        relations.appendChild(ok);
    } else {
        for (const failure of state.failures) {
            const div = document.createElement('div');
            div.className = 'failure';
            const relation = document.createElement('code');
            relation.textContent = tidyRelation(failure.relation);
            const value = document.createElement('code');
            value.textContent = failure.value;
            div.append(
                relation,
                document.createTextNode(' applied to '),
                strong(failure.input_name),
                document.createTextNode(' should be 0, but is '),
                value,
                document.createTextNode('.'),
            );
            relations.appendChild(div);
        }
    }

    const count = selected.size;
    document.getElementById('selection-hint').textContent =
        count === 0
            ? 'Nothing selected. Click cells in the diagram to choose the generators of a submodule.'
            : `${count} cell${count === 1 ? '' : 's'} selected: ` +
              selectedCells()
                  .map(({ degree, idx }) => names(degree)[idx])
                  .join(', ');
    document.getElementById('undo').disabled = undoStack.length === 0;

    document.getElementById('module-name').value = state.name;
    document.getElementById('prime').value = `${state.p}`;
    const actions = document.getElementById('actions');
    // Do not fight the user for the cursor while they are typing in the box.
    if (document.activeElement !== actions) {
        actions.value = state.actions_text;
    }
}

function strong(text) {
    const node = document.createElement('strong');
    node.textContent = text;
    return node;
}

/// `check_validity` writes relations as `1 * Sq1 * Sq1` and leaves an empty factor where one side is
/// the identity, e.g. `1 * Sq10 * `. Tidy that up for display.
function tidyRelation(relation) {
    return relation
        .split('+')
        .map(
            term =>
                term
                    .trim()
                    .split('*')
                    .map(factor => factor.trim())
                    .filter(factor => factor !== '' && factor !== '1')
                    .join(' ') || '1',
        )
        .join(' + ');
}

function setStatus(kind, text) {
    const status = document.getElementById('status');
    status.dataset.state = kind;
    status.textContent = text;
}

let toastTimer = null;
function showToast(message, kind = 'info') {
    const toast = document.getElementById('toast');
    toast.textContent = message;
    toast.dataset.kind = kind;
    toast.hidden = false;
    clearTimeout(toastTimer);
    toastTimer = setTimeout(
        () => {
            toast.hidden = true;
        },
        kind === 'error' ? 7000 : 3000,
    );
}

// ---------------------------------------------------------------- interaction

let drag = null;

/// The cell under the pointer, wherever the event was delivered.
///
/// This cannot use `event.target`: the drag captures the pointer on the `<svg>`, so from then on every
/// event is delivered to the `<svg>` itself rather than to the circle under the cursor. Arcs, labels
/// and the rubber band are all `pointer-events: none`, so what is left under the cursor is either a
/// cell or a degree row.
function cellUnder(event) {
    return cellFrom(document.elementFromPoint(event.clientX, event.clientY));
}

function cellFrom(node) {
    if (node instanceof Element && node.classList.contains('cell')) {
        return {
            degree: Number.parseInt(node.dataset.degree, 10),
            idx: Number.parseInt(node.dataset.idx, 10),
        };
    }
    return null;
}

function degreeFrom(node) {
    if (node instanceof Element && node.dataset.degree !== undefined) {
        return Number.parseInt(node.dataset.degree, 10);
    }
    return null;
}

function setUpDiagram() {
    const svg = document.getElementById('diagram');

    svg.addEventListener('pointerdown', event => {
        // A click on empty space in a row adds a cell there, but only if it turns out not to be the
        // start of a drag, so the decision is deferred to `pointerup`.
        drag = {
            from: cellFrom(event.target),
            fromDegree: degreeFrom(event.target),
            moved: false,
            rubber: null,
        };
        if (drag.from !== null) {
            event.preventDefault();
            svg.setPointerCapture(event.pointerId);
        }
    });

    svg.addEventListener('pointermove', event => {
        if (drag === null || drag.from === null) {
            return;
        }
        drag.moved = true;
        const point = svgPoint(svg, event);
        if (drag.rubber === null) {
            drag.rubber = el('line', { class: 'rubber' });
            svg.appendChild(drag.rubber);
        }
        drag.rubber.setAttribute('x1', x(drag.from.idx));
        drag.rubber.setAttribute('y1', y(drag.from.degree));
        drag.rubber.setAttribute('x2', point.x);
        drag.rubber.setAttribute('y2', point.y);

        // Highlight the rows a generator could reach from the source.
        for (const row of svg.querySelectorAll('.degree-row')) {
            const degree = Number.parseInt(row.dataset.degree, 10);
            const reachable =
                generatorAt(degree - drag.from.degree) !== undefined;
            row.classList.toggle('drop-target', reachable);
        }
    });

    svg.addEventListener('pointerup', async event => {
        if (drag === null) {
            return;
        }
        const { from, fromDegree, moved, rubber } = drag;
        drag = null;
        rubber?.remove();
        for (const row of svg.querySelectorAll('.degree-row')) {
            row.classList.remove('drop-target');
        }

        if (from === null) {
            if (!moved && fromDegree !== null) {
                await run('addGenerator', { degree: fromDegree });
            }
            return;
        }

        const target = cellUnder(event);
        if (
            target === null ||
            (target.degree === from.degree && target.idx === from.idx)
        ) {
            // A click on a cell, or a drag that ended nowhere in particular: toggle its
            // selection, so a set of cells can be built up for the submodule and quotient
            // operations.
            const key = cellKey(from);
            if (selected.has(key)) {
                selected.delete(key);
            } else {
                selected.add(key);
            }
            render();
            return;
        }
        await toggleAction(from, target);
    });

    // A drag released outside the window never produces `pointerup`.
    svg.addEventListener('pointercancel', () => {
        drag?.rubber?.remove();
        drag = null;
        for (const row of svg.querySelectorAll('.degree-row')) {
            row.classList.remove('drop-target');
        }
    });

    svg.addEventListener('dblclick', event => {
        const cell = cellFrom(event.target);
        if (cell !== null) {
            renameCell(cell);
        }
    });
}

/// Convert a pointer event to the SVG's own coordinates, which differ from client coordinates as
/// soon as the diagram is scrolled.
function svgPoint(svg, event) {
    const rect = svg.getBoundingClientRect();
    return { x: event.clientX - rect.left, y: event.clientY - rect.top };
}

async function toggleAction(from, to) {
    const gap = to.degree - from.degree;
    if (gap <= 0) {
        showToast(
            'An operation raises degree, so drag from the lower cell to the higher one.',
        );
        return;
    }
    const generator = generatorAt(gap);
    if (generator === undefined) {
        const available = state.generators.map(
            g => `${g.name} (degree ${g.degree})`,
        );
        showToast(
            `No algebra generator has degree ${gap}, so there is no action to set across that ` +
                `gap. Operations of other degrees are determined by the generators — Sq3 is ` +
                `Sq1 Sq2, for instance — so they cannot be drawn.` +
                (available.length > 0
                    ? `\nGenerators that fit this module: ${available.join(
                          ', ',
                      )}.`
                    : ''),
        );
        return;
    }
    await run('addToAction', {
        opDegree: gap,
        sourceDegree: from.degree,
        sourceIdx: from.idx,
        targetIdx: to.idx,
        coeff: 1,
    });
}

/// Remove every selected cell.
///
/// Removal is by name rather than index: deleting one cell renumbers the ones after it in the same
/// degree, so a list of indices collected beforehand would go stale mid-loop.
async function removeSelected() {
    const doomed = selectedCells().map(({ degree, idx }) => names(degree)[idx]);
    selected.clear();
    for (const name of doomed) {
        const found = findCell(name);
        if (found !== null) {
            if (!(await run('removeGenerator', found))) {
                return;
            }
        }
    }
}

function findCell(name) {
    for (const { degree, names: row } of state.basis) {
        const idx = row.indexOf(name);
        if (idx !== -1) {
            return { degree, idx };
        }
    }
    return null;
}

function renameCell({ degree, idx }) {
    const current = names(degree)[idx];
    const name = window.prompt(`Rename ${current}`, current);
    if (name !== null && name !== current) {
        run('renameGenerator', { degree, idx, name });
    }
}

function setUpKeyboard() {
    document.addEventListener('keydown', event => {
        // Never steal keys from a text field.
        const tag = document.activeElement?.tagName;
        if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') {
            return;
        }
        if (event.key === 'z' && (event.ctrlKey || event.metaKey)) {
            event.preventDefault();
            undo();
            return;
        }
        if (selected.size === 0) {
            return;
        }
        if (event.key === 'Delete' || event.key === 'Backspace') {
            event.preventDefault();
            removeSelected();
        } else if (event.key === 'F2') {
            event.preventDefault();
            const cells = selectedCells();
            if (cells.length === 1) {
                renameCell(cells[0]);
            } else {
                showToast('Select a single cell to rename it.');
            }
        } else if (event.key === 'Escape') {
            selected.clear();
            render();
        }
    });
}

// ------------------------------------------------------------------- controls

function setUpControls() {
    document.getElementById('help-button').addEventListener('click', () => {
        const help = document.getElementById('help');
        help.hidden = !help.hidden;
    });

    document.getElementById('show-derived').addEventListener('change', render);

    document.getElementById('module-name').addEventListener('change', event => {
        run('setName', { name: event.target.value });
    });

    document.getElementById('prime').addEventListener('change', event => {
        const p = Number.parseInt(event.target.value, 10);
        if (p === state.p) {
            return;
        }
        const hasActions = state.actions_text.trim() !== '';
        if (
            hasActions &&
            !window.confirm(
                'Changing the prime discards the actions, since which operations are algebra ' +
                    'generators depends on the prime. The cells are kept. Continue?',
            )
        ) {
            event.target.value = `${state.p}`;
            return;
        }
        run('setPrime', { p });
    });

    document.getElementById('add-cell').addEventListener('click', () => {
        const degree = Number.parseInt(
            document.getElementById('new-degree').value,
            10,
        );
        if (Number.isNaN(degree)) {
            showToast('Enter a degree first.', 'error');
            return;
        }
        run('addGenerator', { degree });
    });

    document.getElementById('clear').addEventListener('click', () => {
        if (!state.is_zero && !window.confirm('Discard the current module?')) {
            return;
        }
        selected.clear();
        view = { min: 0, max: 4 };
        run('create', { p: state.p });
    });

    document.getElementById('extend-up').addEventListener('click', () => {
        view.max += 1;
        render();
    });

    document.getElementById('extend-down').addEventListener('click', () => {
        view.min -= 1;
        render();
    });

    setUpOperations();

    document.getElementById('evaluate').addEventListener('click', evaluate);
    document
        .getElementById('evaluate-input')
        .addEventListener('keydown', event => {
            if (event.key === 'Enter') {
                evaluate();
            }
        });

    document
        .getElementById('apply-actions')
        .addEventListener('click', async () => {
            const error = document.getElementById('actions-error');
            error.textContent = '';
            error.classList.remove('error');
            try {
                setState(
                    await call('setActionsText', {
                        text: document.getElementById('actions').value,
                    }),
                );
            } catch (e) {
                error.textContent = e.message;
                error.classList.add('error');
            }
        });

    document
        .getElementById('upload')
        .addEventListener('change', async event => {
            const file = event.target.files[0];
            if (file === undefined) {
                return;
            }
            // Allow the same file to be picked twice in a row.
            event.target.value = '';
            await load(await file.text(), file.name.replace(/\.json$/, ''));
        });

    document.getElementById('download').addEventListener('click', async () => {
        const json = await call('toJson');
        const name = state.name.trim() === '' ? 'module' : state.name.trim();
        const url = URL.createObjectURL(
            new Blob([`${json}\n`], { type: 'application/json' }),
        );
        const link = document.createElement('a');
        link.href = url;
        link.download = `${name}.json`;
        link.click();
        URL.revokeObjectURL(url);
    });

    document.getElementById('copy').addEventListener('click', async () => {
        await copy(await call('toJson'), 'JSON copied to the clipboard.');
    });

    document.getElementById('permalink').addEventListener('click', async () => {
        const json = await call('toJsonCompact');
        const url = new URL(window.location.href);
        url.search = `?module_json=${encodeURIComponent(json)}`;
        await copy(url.href, 'Link copied to the clipboard.');
    });

    document
        .getElementById('compute-ext')
        .addEventListener('click', async () => {
            if (state.is_zero) {
                showToast('There is nothing to resolve yet.', 'error');
                return;
            }
            if (state.valid === false) {
                showToast(
                    'This is not a module yet: some Adem relations fail, so it cannot be resolved.',
                    'error',
                );
                return;
            }
            const json = await call('toJsonCompact');
            // The spectral sequence viewer is deployed one level up from this page.
            window.open(
                `../?module_json=${encodeURIComponent(json)}`,
                '_blank',
            );
        });
}

async function copy(text, message) {
    try {
        await navigator.clipboard.writeText(text);
        showToast(message);
    } catch (e) {
        showToast(`Could not copy: ${e.message}`, 'error');
    }
}

async function evaluate() {
    const output = document.getElementById('evaluate-output');
    const expr = document.getElementById('evaluate-input').value;
    output.classList.remove('error');
    try {
        const value = await call('evaluate', { expr });
        output.textContent =
            expr.trim() === '' ? '' : `${expr.trim()} = ${value}`;
    } catch (e) {
        output.textContent = e.message;
        output.classList.add('error');
    }
}

// -------------------------------------------------------------- operations

function setUpOperations() {
    document.getElementById('undo').addEventListener('click', undo);

    document.getElementById('op-dual').addEventListener('click', () => {
        run('dual');
    });

    document.getElementById('op-shift').addEventListener('click', () => {
        const by = Number.parseInt(
            document.getElementById('shift-by').value,
            10,
        );
        if (Number.isNaN(by)) {
            showToast('Enter how far to shift.', 'error');
            return;
        }
        run('shift', { by });
    });

    document.getElementById('op-truncate').addEventListener('click', () => {
        // A blank bound means "no bound", which is why these are nullable rather than defaulted.
        const bound = id => {
            const raw = document.getElementById(id).value.trim();
            return raw === '' ? null : Number.parseInt(raw, 10);
        };
        const [min, max] = [bound('truncate-min'), bound('truncate-max')];
        if (Number.isNaN(min) || Number.isNaN(max)) {
            showToast('The bounds must be whole numbers, or blank.', 'error');
            return;
        }
        if (min === null && max === null) {
            showToast('Give at least one bound.', 'error');
            return;
        }
        run('truncate', { min, max });
    });

    document.getElementById('op-submodule').addEventListener('click', () => {
        withSelection(cells => run('submodule', { cells }));
    });

    document.getElementById('op-quotient').addEventListener('click', () => {
        withSelection(cells => run('quotient', { cells }));
    });

    document.getElementById('op-tensor').addEventListener('click', () => {
        withOperand(json => run('tensor', { other: json }));
    });

    document.getElementById('op-sum').addEventListener('click', () => {
        withOperand(json => run('directSum', { other: json }));
    });
}

/// Call `f` with the selection as the flat `[degree, idx, ...]` list the worker expects.
function withSelection(f) {
    if (selected.size === 0) {
        showToast(
            'Select the cells to generate the submodule first, by clicking them in the diagram.',
            'error',
        );
        return;
    }
    f(selectedCells().flatMap(({ degree, idx }) => [degree, idx]));
}

async function withOperand(f) {
    const name = document.getElementById('operand').value;
    if (name === '') {
        showToast('Pick a module to operate with.', 'error');
        return;
    }
    try {
        const response = await fetch(`./steenrod_modules/${name}.json`);
        if (!response.ok) {
            throw new Error(`${response.status}`);
        }
        const spec = JSON.parse(await response.text());
        // Most library files carry no `name`, and the result is named after its factors, so fall back
        // to the file name the way opening a module does.
        if (spec.name === undefined) {
            spec.name = name;
        }
        await f(JSON.stringify(spec));
    } catch (e) {
        showToast(`Could not load ${name}: ${e.message}`, 'error');
    }
}

// ------------------------------------------------------------- library & files

async function setUpLibrary() {
    const select = document.getElementById('library');
    const hint = document.getElementById('library-hint');
    let modules;
    try {
        const response = await fetch('./modules.json');
        if (!response.ok) {
            throw new Error(`${response.status}`);
        }
        modules = await response.json();
    } catch (e) {
        hint.textContent = `Could not load the module library (${e.message}).`;
        return;
    }

    const loadable = modules.filter(
        m => m.type === 'finite dimensional module',
    );
    // The same list serves the library and the operand of a tensor or direct sum.
    const operand = document.getElementById('operand');
    for (const module of loadable) {
        const option = document.createElement('option');
        option.value = module.name;
        option.textContent = `${module.name}  (p = ${module.p})`;
        operand.appendChild(option.cloneNode(true));
        select.appendChild(option);
    }
    const skipped = modules.length - loadable.length;
    hint.textContent =
        skipped === 0
            ? ''
            : `${skipped} module${skipped === 1 ? '' : 's'} in the library ` +
              'are not finite dimensional, and are not listed.';

    const open = async () => {
        const name = select.value;
        if (name === '') {
            return;
        }
        const response = await fetch(`./steenrod_modules/${name}.json`);
        await load(await response.text(), name);
    };
    document.getElementById('load-library').addEventListener('click', open);
    select.addEventListener('dblclick', open);
}

/// Replace the open module with the one described by `json`.
async function load(json, name) {
    try {
        const next = await call('load', { json });
        selected.clear();
        view = {
            min: Math.min(0, next.min_degree),
            max: Math.max(4, next.top_degree),
        };
        if (next.name === '' && name !== undefined) {
            setState(await call('setName', { name }));
        } else {
            setState(next);
        }
        showToast(`Opened ${next.name === '' ? name ?? 'module' : next.name}.`);
    } catch (e) {
        showToast(e.message, 'error');
    }
}

const STORAGE_KEY = 'module_builder.module';

function save() {
    // Best effort: private browsing modes and full quotas both throw here, and neither is worth
    // interrupting the user over.
    try {
        window.localStorage.setItem(STORAGE_KEY, JSON.stringify(state.json));
    } catch (e) {
        /* ignore */
    }
}

/// Restore the module from the URL, then from the last session, then fall back to empty.
async function restore() {
    const params = new URLSearchParams(window.location.search);

    const fromUrl = params.get('module_json');
    if (fromUrl !== null) {
        await load(fromUrl);
        return;
    }

    const fromLibrary = params.get('module');
    if (fromLibrary !== null) {
        try {
            const response = await fetch(
                `./steenrod_modules/${fromLibrary}.json`,
            );
            if (response.ok) {
                await load(await response.text(), fromLibrary);
                return;
            }
        } catch (e) {
            showToast(`Could not open ${fromLibrary}: ${e.message}`, 'error');
        }
    }

    const stored = window.localStorage.getItem(STORAGE_KEY);
    if (stored !== null) {
        try {
            setState(await call('load', { json: stored }));
            return;
        } catch (e) {
            // A stored module we can no longer read is not worth a complaint; start fresh.
        }
    }

    setState(await call('create', { p: 2 }));
}

// ------------------------------------------------------------------- start up

(async () => {
    setUpDiagram();
    setUpKeyboard();
    setUpControls();
    const library = setUpLibrary();
    try {
        await restore();
    } catch (e) {
        showToast(`Could not start: ${e.message}`, 'error');
        setStatus('loading', 'failed');
        return;
    }
    await library;
})();
