# Review: shared input validation (`check_inputs`, `check_symbol_clash`)

Range `5a48f91..735cfba` (70f6db4 "Reject inputs and params that share a name with a
different symbol", 735cfba "Share input validation between to_standard_form and
linearize"). Tests: `tests/test_inputs_adversarial.py` (93 pass, 11 strict xfail). Full
suite with the new file: 3101 passed, 13 xfailed (11 new + 2 existing), ruff clean.

**Verdict: GO WITH CHANGES.** The core fix holds: a constant, another class or other flags
under an input's or param's name is rejected by all three functions, both input functions
give the same type and message, every system builder of the test suite still passes, and
the cost is negligible. The findings are small (one false positive on an unusual size
spelling, two message hints, two gaps that let a confusable input through). None changes a
derived equation. F4 and F5 are cheap to fix and worth doing before release; F1 to F3 can go
to `specs/follow-ups.md`.

## Findings (ranked)

### F4 (medium-low): an input printed like a node of the equations passes `check_inputs`

`check_symbol_clash` collects only leaves (`geomech/dynamics/eom.py:187`), so a plain leaf
printed like a non-leaf node of the equations is not compared.

```python
x, F, c = Vector("x"), Vector("F"), Scalar("c", attr=["Constant"])
eom = compute_eom(½ m ẋ·ẋ, Dot(x.delta(), F - c*x.t_diff()), V)   # holds d/dt(x)
to_standard_form(eom, V, [F, Vector(str(x.t_diff()))])   # accepted
linearize(eom, V, [F, Vector(str(x.t_diff()))])          # accepted, zero B column
```

Observed: both accept; `linearize` returns a zero column `\delta{\frac{d}{dt}(x)}`.
`StateLayout` (`geomech/codegen/layout.py:189`) then rejects the same input ("the input
'\frac{d}{dt}(x)' prints like the velocity"). Expected: `check_inputs` rejects it, so the
three entry points agree. The states check (`_input_state_symbols`) does not catch it because
`δ(ẋ)` is rewritten to `d/dt(δx)` and the input's variation prints `\delta{\frac{d}{dt}(x)}`.
Fix: in `check_inputs`, also reject an input whose name equals `str()` of a velocity or
acceleration of a configuration variable (or collect printed names of every node in the
walk). Test: `test_f4_input_printed_like_the_velocity_is_rejected`.

### F5 (medium-low): a set or dict of inputs is accepted; input order then depends on the hash seed

`check_inputs` accepts any iterable (`eom.py:254-260`). For inputs, order matters (columns
of `B`, `lin.inputs`, the codegen `u` vector).

```python
s = vbl_systems.quadrotor(); eom = compute_eom(s.L, s.dW, s.variables)
linearize(eom, s.variables, set(s.inputs)).inputs
# PYTHONHASHSEED=1: [δf, δM]   PYTHONHASHSEED=2: [δM, δf]
```

A dict is iterated by keys (`{F: 1}` acts as `[F]`). Expected: `TypeError` for `set`/`dict`
(as for one leaf), or document that the order is the iteration order. Docs
(`docs/conventions.md`) list "a list (or tuple, generator)". Pre-existing in `linearize`,
now shared. Test: `test_f5_unordered_inputs_rejected`.

### F1 (low): a size spelled as a list or int is "a different symbol"

`_symbol_key` (`eom.py:150`) compares `size` as stored. `Vector("F", size=[3])` and
`Vector("F", size=3)` store `[3]`/`3`, the default leaf `(3,)`.

```python
eom = compute_eom(½ m ẋ·ẋ, Dot(x.delta(), Vector("F")), V)
to_standard_form(eom, V, [Vector("F", size=[3])])
# ValueError: input 'F' (Vector) has the name of the Vector 'F' in the equations, a different symbol; ...
```

Same for a `Matrix` param with `size=[3, 3]` in `compute_eom`. Expected: the same symbol
(normalise with `tuple(np.atleast_1d(size))`, or normalise `size` in the leaf
constructors as `_zero_leaf` already does with `tuple(leaf.size)`). No call in the repo
spells a size as a list, so impact is low. Tests: `test_f1_*`.

### F2 (low, message): a size-only clash does not say what differs

`Vector("F", size=(6,))` against `Vector("F")`: "input 'F' (Vector) has the name of the
Vector 'F' in the equations, a different symbol; … build the input with the class and attr
it has in the equations" (`eom.py:204-209`). Both sides are described identically and the
hint omits size. Expected: `_describe_leaf` adds the size when it is not the default (or
when the sizes differ), and the hint says "class, attr and size". Test:
`test_f2_size_clash_message_names_the_sizes`.

### F3 (low, message): wrong hint against a named numeric leaf

`L` holds `Scalar("k", value=2.0)` (a named numeric symbol, decision D1); input
`Scalar("k")`: "… a symbol is either a constant parameter or an input (rename one, or drop
attr=['Constant'] if 'k' is an input)" (`eom.py:198-203`). `k` has no `Constant` attribute;
its `value=` makes it constant. Same for a `compute_eom` param. Expected: the hint names
`value=` when the leaf is numeric. Test: `test_f3_numeric_leaf_hint_names_value`.

## What held

1. **False positives.** Every builder of `tests/codegen_systems.ALL_SYSTEMS` (oracle systems,
   the nine `vbl_systems`, the phase-1 adversarial systems) runs compute_eom →
   to_standard_form → linearize with its own inputs without an input error. Also accepted: a
   `UnitNorm` input used as the same symbol in `δW`; a rebuilt input with `attr` as a tuple;
   flags given in another order (flags compare as a dataclass, order-free); an input that is
   also a param of `L` and a force in `δW` (`G[u] = −e₁`, matching `B`); a param used only
   inside `d/dt(l)`, also passed as a generator of rebuilt leaves; a rebuilt `Vector` input
   with `value=` (the key ignores `value`; `Vector` sets no flag for it).
2. **False negatives rejected.** `Transpose(F)`, `Matrix`, a string, a list of strings, a
   dict of strings (`TypeError`); `Scalar("F")` against `Vector("F")`; another `attr`;
   a literal of the equations (`Scalar("(-2)")` against the `(-2)` of `−k x·x`, "constant");
   `Vector("0v")` against `ZeroVector`; a Vector param against a Scalar `l` that occurs only
   under `d/dt` (the walk enters calculus nodes). An input used nonlinearly gives
   `NonLinearError` naming `u` in `to_standard_form` and the correct `2u e₁` column in
   `linearize`. An input named like a tangent vector of an `S2` that is not in the system
   cannot clash in practice: such a `TS2` cannot reach the equations through `compute_eom`.
3. **Key comparison.** Equal printed name and class with different `value=` on two Scalars
   is not distinguished, but both are numeric and an input or param cannot be numeric, so
   this is unreachable through the public functions.
4. **Message parity.** Same type and message after the prefix in both functions, including
   the new duplicate-of-another-class, other-size and other-class cases. The four new
   troubleshooting rows (`docs/troubleshooting.md`) are exercised by
   `tests/test_troubleshooting_doc.py` and match. Not a regression but noted: an input that
   is also a declared param of `compute_eom` still fails in the engine with an
   `AlgebraicError` about `d/dt(l)` that does not name the cause (as before the range); and
   `to_standard_form` on a non-dict `eom` still fails with `AttributeError` (now from
   `eom_dict.values()` instead of `.items()`), where `linearize` raises `TypeError`.
5. **Performance.** Best of N, this machine, base `5a48f91` vs `735cfba`:

   | system | compute_eom | to_standard_form | linearize | `check_inputs` | clash on `L`,`δW` |
   |---|---|---|---|---|---|
   | quadrotor with load | 12.7 → 12.9 ms | 1.0 → 1.1 ms | 7.3 → 7.8 ms | 0.06 ms | 0.02 ms |
   | double rigid pendulum | 19.7 → 20.4 ms | 1.9 → 2.1 ms | 21.8 → 22.4 ms | 0.10 ms | 0.02 ms |

   The differences are within noise; the check itself is under 5 % of `to_standard_form`.
   The walk is not memoised, so a heavily shared DAG is walked once per path, but so is
   `str()` of the same tree, which the engine needs anyway.

## Fix log

Fixed on `feature/codegen-simulation` on top of a4a849c (not committed). All 11 strict
xfails of `tests/test_inputs_adversarial.py` flipped by removing the markers only; no
assertion changed. One docstring of a passing test was updated
(`test_input_that_is_also_a_param_fails_in_the_engine`: the AlgebraicError now comes from
`check_inputs`; its assertions are unchanged).

| Finding | Fix | Message |
|---|---|---|
| F1 | `geomech/core/base/expressions.py:66` `_leaf_size`, called by `Vector.__init__` (:536, :539) and `Matrix.__init__` (:654): `size` is stored as a tuple of ints (`3`, `[3]` → `(3,)`). `_symbol_key` (`eom.py:150`) is unchanged and now compares equal sizes. | none (accepted) |
| F2 | `eom.py:250` `_clash_against`: when the classes match and the sizes differ, `_describe_leaf(..., with_size=True)` appends the sizes and the hint says "class, attr and size". The class-only hint ("class and attr") is unchanged. | `input 'F' (Vector of size (6,)) has the name of the Vector of size (3,) 'F' in the equations, a different symbol; … build the input with the class, attr and size it has in the equations)` |
| F3 | `eom.py:240-249`: the hint depends on the leaf in the equations: `attr=['Constant']` → drop it; a number literal → rename the input; a named numeric leaf → drop `value=`. Same path for `compute_eom` params. | `input 'k' has the name of the constant 'k' in the equations; a symbol is either a constant parameter or an input (rename one, or drop value= from the leaf 'k' if it is an input)` |
| F4 | `eom.py:187` `_printed_nodes`: one walk, a dict by `str` of every node (leaves and operations); `_clash_against` (`eom.py:225`) rejects a symbol printed like an operation node. Operations whose printed length is bounded below (sum of their operands' lengths) by more than the longest looked-up name are not printed, so the walk stays linear. `check_inputs` also rejects an input printed like the velocity or acceleration of a configuration variable (`_kinematic_names`, `eom.py:310`, check at :400) before the walk, so `StateLayout`, which has no equations, agrees: `StateLayout.__init__` now calls `check_inputs(inputs, variables, (), "StateLayout")` (`geomech/codegen/layout.py:165`) instead of its own type check. | `input '\frac{d}{dt}(x)' prints like the velocity '\frac{d}{dt}(x)' of the configuration variable 'x'; leaves are matched by printed name, so the two would be confused (rename the input)`; for another node: `input '\frac{d}{dt}(l)' has the printed name of the TimeDerivative node '\frac{d}{dt}(l)' in the equations, not of a leaf; … (rename the input)` |
| F5 | `eom.py:366`: after the per-entry type check, a `collections.abc.Set` or `Mapping` raises `TypeError`; generators are still accepted. The per-entry type check runs first so that a dict of strings keeps its "got str" message (an existing passing test). | `inputs must be an ordered sequence (list or tuple), got set; a set has no order, so the columns of B would depend on the hash seed` / `… got dict; a mapping is iterated by its keys (pass a list of the inputs)` |
| input that is a param | `eom.py:415-424`: an input whose `d/dt(u)` is a node of the equations raises `AlgebraicError` in `check_inputs` (both callers). `to_standard_form` cannot know the params (they are not in the `eom` dict), so the check uses the equations: a param whose rate does not occur may still be an input (`test_input_in_the_lagrangian_and_the_virtual_work`). The type stays `AlgebraicError` with "calculus operator", as the existing test asserts. | `input 'l' occurs in the equations inside the calculus operator '\frac{d}{dt}(l)'; an input enters the equations as itself, not through its rate (is 'l' a prescribed parameter of compute_eom(..., params=[...])? It can be an input only if its rate does not occur)` |
| non-dict `eom` | `geomech/dynamics/standard_form.py:95` | `to_standard_form: eom must be the dict returned by compute_eom, got <class 'list'>` |

Also changed: the duplicate-input message adds "(an entry that prints like 'f' is the
same symbol)", so `StateLayout`'s duplicate case keeps matching its existing test
("'f'.*prints like").

Not done: rejecting a `set` for `compute_eom(..., params=...)`. The params are order-free
(only their printed names are used) and `tests/test_phase2_adversarial.py`
`test_sequence_forms_for_params_are_equivalent` asserts that `params={l}` gives the same
equations; it stays accepted and documented as such.

Docs: `docs/troubleshooting.md` gains seven rows (set of inputs, velocity look-alike,
operation-node look-alike, param with a rate as input, size-only clash, named numeric
leaf, non-dict `eom` in `to_standard_form`) and the duplicate row's new text; each has a
trigger in `tests/test_phase4_adversarial.py` `TRIGGERS`. `docs/conventions.md` (inputs),
the docstrings of `to_standard_form`, `linearize`, `check_inputs`, `check_symbol_clash`,
`StateLayout` and `CHANGELOG.md` are updated.

Checks: `ruff format --check` and `ruff check` on `geomech/ tests/ scripts/` clean; `mypy
geomech --ignore-missing-imports` clean; full suite with `-W error`: 3192 passed, 16
xfailed (was 3174 passed, 27 xfailed at a4a849c: 11 flipped, 7 new trigger cases);
`mkdocs build --strict` passes; `scripts/regenerate_generated.py --check` clean (18
modules). Pinned snapshot against a4a849c: 5983 outputs unchanged, 0 strings changed, 0
value changes, `check_eom` 243 compared and none lower; the 6 "failures" are the six
former-xfail calls of F4/F5 that now raise as intended. Cost (best of 50, this machine):
`check_inputs` 0.11 ms (quadrotor with load) and 0.17 ms (double rigid pendulum), against
0.105/0.093 ms before and `to_standard_form` 1.1/2.1 ms.
