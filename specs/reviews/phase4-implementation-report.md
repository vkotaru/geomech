# Phase 4 implementation report: debuggability

Branch `fix/phase4-debuggability`, created from `57232da` (`fix/phase3-operators-printing`).
Nothing is committed. Plan: `specs/fix-plan.md`, phase 4, items 4.1 to 4.7, decision D5
(trace is a keyword argument), plus the phase-4 follow-ups of `phase1-final.md`,
`phase2-final.md` (F-1, F-2, F-6) and `phase3-final.md` (F1, F2, F7).

The goal: "when I debug in 6 months or 1 year from now, it should be easy to debug it."

## Summary

| Check | Result |
|---|---|
| Suite, `-W error` | **2011 passed, 5 xfailed** (base 1931 passed, 5 xfailed). 76 tests in five new files, plus 4 from the two new oracle systems (`test_eom_oracle` and the `test_phase1_adversarial` loop over its `SYSTEMS`). No existing test function was edited. |
| Suite, `GEOMECH_DEBUG=1 -W error` | **2011 passed, 5 xfailed** (every pipeline call also runs `check_tree` on its input and output, and logs every stage). |
| `check_tree` in the suite | `tests/conftest.py` (session, autouse) wraps the pipelines in every `tests.*` module: one full run checks the trees of **435 `compute_eom`, 187 `to_standard_form` and 112 `linearize` calls**; none violates an invariant. |
| `str` byte-identity | **Identical.** All 65 argument-free builders of `tests/vbl_systems.py`, `tests/phase1_systems.py` and `tests/test_eom_oracle.py` (the 61 of earlier phases plus the 4 with default arguments; the 2 new systems excluded) through `str(L)`, `str(dW)`, δS, `apply_manifold_rules`, `full_simplify`, `integrate_by_parts`, `expand`, `compute_eom`, `to_standard_form` (every field) and `linearize` (every attribute), in a worktree at `57232da` and in this tree: 379 556 bytes at base. The only 2 differing lines are the two expected `AlgebraicError` texts (`so3_/s2_outer_product_on_velocity`), which now start with `compute_eom: while integrating by parts (stage ibp, …):` followed by the unchanged engine message; with the error lines removed the dumps are byte-identical. The worktree was removed. |
| Runtime, happy path | `compute_eom` on the quadrotor with a cable-suspended load (median of 40, three alternating rounds): **base 10.9 / 10.9 / 10.9 ms, now 10.9 / 10.9 / 10.9 ms**. Suite wall time: base 6.83–6.90 s (1931 tests); now 7.48–7.67 s (2011 tests), and **6.89–7.12 s for the same 1933 pre-existing tests** (the difference is the conftest `check_tree` pass, test time only). |
| `ruff format --check` / `ruff check` on `geomech/ tests/` | clean (92 files) |
| `ruff check examples/` | 111 errors, all pre-existing (star imports; identical count at base); the contributor notes exempt the example scripts. No example was changed. |
| `mkdocs build --strict` | exit 0; the only INFO line left is the existing `reference/` link (the `manifold-implementation-notes.md` "not in nav" line is gone). |
| `engine_fuzz.py --quick -W error` | exit 0, 36 distinct findings, identical to base except the `collect.py` line numbers in "where:". 0 wrong-value. |
| Examples (6), `-W error` | all exit 0, output byte-identical to base. |
| Notebooks (4), `nbconvert --execute` | all execute; every cell output identical to base (18, 28, 12, 9 outputs). |

New test files: `tests/test_trace.py` (13), `tests/test_invariants.py` (24),
`tests/test_rule_names.py` (19), `tests/test_logging.py` (7), `tests/test_error_context.py`
(13); new `tests/conftest.py`; two systems appended to `tests/test_eom_oracle.py`
(`row_kinetic_energy`, `row_kinetic_energy_so3`, added to `SYSTEMS`).

## 4.1 Trace objects

- `geomech/utils/trace.py:171` `_Trace`: ordered `stages` list of `(name, expr)`;
  `__getitem__` accepts a stage name or an attribute name; `format_summary()` (`:201`, stage
  name, node count, first 200 characters of `repr`), `summary()` (`:217`, prints it),
  `show(stage)` (`:221`, `print_tree(style="indent")`, dict stages printed per key),
  `str(trace)` is the summary.
- `trace.py:238` `EOMTrace`: `dL`, `dW`, `dS`, `after_manifold_rules`, `after_simplify`,
  `variation_vectors`, `ibp_targets`, `after_ibp_prepare`, `ibp` (per target:
  `collect`/`integrate`/`simplify`), `after_ibp`, `after_expand`, `raw_equations`,
  `equations`.
- `trace.py:309` `LinearizationTrace`: per row `dE`, `buckets` (each bucket's terms summed),
  `raw_blocks` (before `full_simplify` and before the mass negation), `blocks`.
- `geomech/dynamics/eom.py:339` `compute_eom(..., *, trace=None)` (D5: keyword argument);
  stages recorded at `eom.py:428-466`; per-target IBP stages through the new
  `geomech/core/math/ibp.py:196` `_integrate_by_parts(expr, targets, record)` (the public
  `integrate_by_parts` delegates to it with `record=None`); extraction stages in
  `eom.py:518` `_separate_variations`.
- `geomech/dynamics/vbl.py:183` `linearize(..., params=(), *, trace=None)`; records at
  `vbl.py:294-319`.
- The recorder `trace.py:140` `stage_recorder(trace, where)` returns `None` when there is
  no trace and the logger is not at `DEBUG`; every recording site is
  `if record is not None:` — nothing is built on the happy path.
- Exports: `geomech/__init__.py:70` `apply_manifold_rules`, `:116` `EOMTrace`,
  `LinearizationTrace`, `set_trace`; `apply_manifold_rules` also from
  `geomech.core.transformations`.
- A trace is filled as the pipeline runs, so after an exception it holds every stage
  before the failure (`tests/test_trace.py::test_trace_keeps_the_stages_before_a_failure`).
  `test_each_stage_reproduces_from_the_previous_one` rebuilds every stage of the rigid
  pendulum from the previous one with public functions.

## 4.2 Errors with pipeline context

All re-raises keep the exception type and chain the original (`raise type(err)(…) from err`).

- `eom.py:454` integration by parts: `AlgebraicError`/`NonLinearError` from `collect` is
  prefixed `compute_eom: while integrating by parts (stage ibp, targets …):` (the phase-1
  follow-up for `collect.py:55-65`).
- `eom.py:540` `_extraction_error`: when a time derivative of a listed variation vector is
  left in δS, the message names the variation vector being extracted, the derivative,
  the term of δS that holds it (preferring a term holding the vector being extracted), says
  that `collect` could not expose it as the left operand of a `Dot`, points to
  `print_tree` and the `ibp.collect[...]` trace stage, and keeps the engine text after
  `(engine: …)`. If the term holds a row inner product it adds "write it as
  Dot(v, A * w)". Any other extraction failure is prefixed
  `compute_eom: while extracting the equation for '…' from δS (stage expand):`.
  `separate_variations` (public) gets the same context.
- `geomech/dynamics/standard_form.py:166` `_has_in_term`: names the equation and the
  term; for an acceleration it says "a third or higher time derivative of the
  configuration (a jerk): the equation is not second order", for an input that it occurs
  only under a derivative. `:189` `_extract` wraps the two extractions; the
  bilinear-acceleration `NonLinearError` names the equation. The jerk and the
  surviving-derivative messages now start differently (scenario (a) below).
- `vbl.py:589` `_check_coefficient_free` keeps the inner `AlgebraicError` text and chains
  it.
- Phase-2 F-1 (`vbl.py:497` `_no_perturbation_hint`, `:470` `_param_names`): the "no
  perturbation symbol" message names the stray variation (`It carries '\delta{l}' (of
  'l')`) and lists every remedy, including "linearize does not support prescribed
  parameters yet". With `linearize(..., params=[l])` (new optional argument, validated
  like `compute_eom`'s) the message says only that prescribed parameters are not
  supported, without the input/`Constant` suggestions. The prefix "contains no
  perturbation symbol" that existing tests match is unchanged.
- Phase-3 F1, `geomech/core/math/collect.py:63`: `case VVMul() if expr.has(vec)` rewrites
  an inner product holding the target to `Dot(column_of(row), column)` (rule
  `collect.inner_product_to_dot`) and continues with the `Dot` rules. Only inner products
  that contain the target are touched; before, any such product made extraction fail, so
  **no pinned string changes** (dump identical). `½ ẋᵀAẋ − ½ xᵀKx` gives
  `(A\frac{d}{dt}(\frac{d}{dt}(x))(-0.5)+(A)'\frac{d}{dt}(\frac{d}{dt}(x))(-0.5)+Kx(-1.0)+f)`
  (`assert_eom` 1.7e-9), and `½ ΩᵀJΩ` prints exactly the free-rigid-body equation of
  `½ Ω·JΩ` (Lee et al. Eq. 6.16; 5.1e-9). Both are in the oracle `SYSTEMS` list and pinned
  in `tests/test_error_context.py`. The stale "write Dot(v, K * v) in a kinetic energy"
  comments in `README.md` and `docs/getting-started.md`, and the CHANGELOG line saying rows
  in a kinetic energy are unsupported, were updated.

## 4.3 Logging

- `trace.py:34` `logger = logging.getLogger("geomech")`; one `DEBUG` line per stage
  (`compute_eom: <stage> [<n> nodes] <repr ≤ 200 chars>`) from `compute_eom`,
  `integrate_by_parts` (per target), `to_standard_form` (`M[...]`, `f[...]`, `G[...]` per
  equation) and `linearize` (`dE`, `bucket`, `raw`, `block` per row).
- `trace.py:126` `configure_from_env()`, called at the end of `geomech/__init__.py:121`:
  `GEOMECH_DEBUG=1` sets the level to `DEBUG` and attaches a `stderr` handler if neither the
  logger nor `geomech` has one; `trace.py:49` `debug_enabled()` is read at pipeline entry
  and turns on `check_tree` (inputs and equations of `compute_eom`, `eom.py:424-427` and
  after extraction; equations and blocks of `to_standard_form`; equations and blocks of
  `linearize`). Values `""`, `0`, `false`, `no` mean off.
- No `print` in library code: `render_eom`'s two status prints became
  `logger.warning`/`logger.info` (`geomech/utils/printing.py:890`). The remaining `print`s
  are the display functions themselves (`print_eom`, `print_tree`, `_Trace.summary`,
  `_Trace.show`); `tests/test_logging.py::test_no_print_in_library_code` enforces this by AST.

## 4.4 Invariants

- `geomech/utils/invariants.py:77` `check_tree(expr)` (iterative walk, shared subtrees
  once) and `geomech/utils/errors.py:98` `InvariantError(Exception)` (the package has no
  common base error class; exported from `geomech` and `geomech.utils`). Checks: every
  child is an expression; unary nodes have one operand; every node has a type; sum
  children have the sum's type and are not the same sum class (nested `Add` via a list
  argument); binary operand types per class (`MVMul` (M, V) or (row, M)); `VVMul` exactly
  one row; a column times a matrix and a matrix times a row; `Dot`/`Cross`/`Hat` take
  columns (phase-3 F2); `VAdd` does not mix rows and columns, zero addends exempt (phase-3
  F2); no `TimeDerivative`/`Variation` directly over an `S2`/`SO3` leaf, no δ of a δ, no δ of
  a variation vector; `TS2`/`TSO3` have their point and are named like its tangent or
  variation vector; `Scalar.value is not None` iff `is_numeric`, a real value, and a
  parenthesised numeric name (`(-1)`) carries that value. Messages name the invariant,
  the path (`root.1.R.L`) and the subtree.
- Phase-3 F2 at construction: I added the row checks to the `Dot`/`Cross`/`Hat`
  constructors and the orientation check to `VAdd` first. Two blockers:
  `tests/test_geometry.py::TestTransposeOperators::test_add_vector`/`test_sub_vector`
  build `Transpose(x) + y` and assert a `VAdd`, and `tests/test_latex_precedence.py:50`
  builds `Cross(a, Transpose(b))` at collection. `Dot`/`Hat` had no blocker, but the
  check costs ~0.4 ms (4 %) of `compute_eom` on the quadrotor with load (`is_row_vector`
  on every `Dot` built), which breaks the zero-cost constraint. **All four are therefore
  in `check_tree` only**, consistently; the constructors are unchanged.
- Suite: `tests/conftest.py` (autouse, session) replaces `compute_eom`,
  `to_standard_form` and `linearize` in every loaded `tests.*` module, on `geomech`,
  `geomech.dynamics` and in the defining modules, with wrappers that run `check_tree` on
  the inputs and every output tree after a *successful* call (tests that expect a raise
  are unaffected); `CHECKED_CALLS` counts them.
- Phase-2 F-2, `eom.py:87-92`: a `SystemVariables` entry with children (`x + y`,
  `m * z`, `Cross(a, b)`) raises `TypeError: … entry '(x+y)' is a VAdd expression, not a
  single symbol; a configuration variable is one Scalar, Vector, S2 or SO3 leaf`. It runs
  after the existing `TimeDerivative`/`Variation` check, so those keep their message.

## 4.5 Named rules

- `trace.py:89` `rule(name, before, after)`: returns `after`; when tracing is on it logs
  `name: before ==> after` on `geomech.rules` (repr, each side cut at 300 characters),
  skipping calls whose output prints like the input (the `*_fold` rules are called on
  every visit of their node class). Off, it is one global read.
- `trace.py:63` `set_trace(on)` and `GEOMECH_TRACE=1` (read at import).
- 76 names in 7 modules. Each module docstring has a "Rules" block (name, identity,
  hypothesis); `tests/test_rule_names.py` parses the docstring with `ast` and the
  literal first argument of every `rule(...)` call (and the name forwarded through
  `expand._distribute`) and asserts the sets are equal per module; a second test asserts
  every `rule()` call has a literal name.
- Firing counts for one `compute_eom` of the quadrotor with load: 877 `rule` calls;
  overhead not measurable (table above).

### Rule-name table

| Module | Rules | Hypotheses the code checks |
|---|---|---|
| `vector_rules.py` (12) | `vector.dot_zero`, `vector.dot_cross_orthogonal`, `vector.dot_unit_norm`, `vector.dot_tangent_orthogonal`, `vector.cross_self`, `vector.cross_zero`, `vector.cross_scaled_self`, `vector.s2_tangent_recovery`, `vector.s2_tangent_negated`, `vector.s2_double_cross`, `vector.transpose_transpose`, `vector.hat_zero` | zero flag; operand equal to a Cross operand; both unit norm **and** equal; `TS2` whose S² point equals `q` by name; tangency = `TS2` of `q` or `Cross(·, q)` (never `ω̇`); the docstring states that the pass recurses only through `Add`/`Mul`/`VAdd`/`SVMul`/`Cross` (audit C1) |
| `simplify.py` `_eliminate` (34) | `simplify.add_fold`, `vadd_fold`, `madd_fold`, `mul_fold`, `svmul_zero`, `svmul_one`, `smmul_zero`, `smmul_one`, `mvmul_zero`, `mvmul_identity`, `mvmul_hat_to_cross`, `mvmul_reassociate`, `mvmul_combine`, `mmmul_zero`, `mmmul_identity`, `orthogonal_transpose`, `vvmul_zero`, `dot_zero`, `cross_zero`, `hat_zero`, `hat_vee_skew`, `vee_zero`, `vee_hat`, `transpose_transpose`, `transpose_symmetric`, `transpose_hat_skew`, `transpose_scaled`, `transpose_scaled_row`, `transpose_product`, `variation_zero`, `time_derivative_zero`, `time_derivative_integral`, `time_integral_zero`, `time_integral_derivative` (all `simplify.`) | zero/one/identity flags; orthogonal flag and equal by name; skew-symmetric (`Hat`, `SkewSymmetricMatrix` leaf, sum/scaling/transpose of such); symmetric flag; float cancellation `1e-12` |
| `pull.py` (7) | `pull.dot_scaled`, `pull.svmul_nested`, `pull.cross_scaled`, `pull.mvmul_scaled`, `pull.smmul_nested`, `pull.mmmul_scaled`, `pull.vvmul_scaled_row` | operand classes only; scaled row `SVMul(row, s)` on the left |
| `expand.py` (13) | `expand.mul_over_add`, `dot_over_vadd`, `cross_over_vadd`, `svmul_over_sum`, `row_svmul_over_sum`, `mvmul_over_sum`, `row_mvmul_over_sum`, `smmul_over_sum`, `mmmul_over_madd`, `vvmul_outer_over_vadd`, `vvmul_inner_over_vadd`, `hat_over_vadd`, `vee_over_vadd` | none (bilinearity); `vee_over_vadd` is documented as unreachable (the operand of `Vee` is a matrix; Vee over `MAdd` is phase 5.2) |
| `manifold_rules.py` (2) | `manifold.s2_variation_omega`, `manifold.so3_variation_Omega` | the tangent vector belongs to a listed point (by printed name) |
| `collect.py` (5) | `collect.inner_product_to_dot`, `collect.flip`, `collect.triple_product`, `collect.peel_matrix`, `collect.strip_rotation` | target inside the inner product; target in the vector operand, not the matrix; both rotations `SO3` and equal by name |
| `ibp.py` (3) | `ibp.dot`, `ibp.dot_prefactor`, `ibp.scalar_factor` | target is the left `Dot` operand; prefactor constant or not; target a direct `Mul` factor |

## 4.6 Un-shadowed submodules

- `geomech/core/transformations/__init__.py:13` imports the submodules themselves (`from
  geomech.core.transformations import expand, manifold_rules, pull, simplify,
  vector_rules`) and no longer rebinds the names to functions; `geomech/core/math/__init__.py:10`
  does the same for `collect`.
- `geomech/utils/callable_module.py:15` `CallableModule`: the four transformation modules
  and `collect` get this `ModuleType` subclass, whose `__call__` forwards to the function
  of the module's own name. So `import geomech.core.transformations.simplify as m` gives
  the module, and every existing import path keeps working: `from
  geomech.core.transformations import expand, full_simplify, simplify` (used by
  `tests/test_operators.py:40`, `tests/test_simplification.py:3`, `tests/test_operations.py:14`,
  `tests/test_simplify.py:726`) returns callables. `geomech.expand` etc. are imported
  from the submodules in `geomech/__init__.py:66-72` and are the plain functions. Renaming
  the modules (`_expand.py`) was rejected: it would make `geomech.core.transformations.simplify`
  the function, which is exactly what the item forbids. Grepped: repo, examples,
  notebooks, docs, a downstream application, the accompanying paper; no other import of
  these names.
- `manifold_rules.py` now imports `SystemVariables` under `TYPE_CHECKING` only (the package
  `__init__` importing it created an import cycle through `geomech.dynamics`).

## 4.7 Docs

- `docs/pipeline.md` (new): the code's stage order, one section per stage (validation,
  variation, manifold rules, simplify with the three passes, IBP with collect and the
  IBP pass per target, expand, extraction, the stray check, `to_standard_form` and
  `linearize`), how to see each stage (`EOMTrace`, `LinearizationTrace`, the logger,
  `GEOMECH_DEBUG`), real `summary()` output for the rigid pendulum, the rule trace, which
  rules fire for the rigid pendulum and a prefix → module → stage table.
- `docs/troubleshooting.md` (new): symptom → cause → fix table with every silent case of
  phases 1–3 marked *fixed* (attr typo, unlisted variable, `−ḟ` from `δẋ` in δW, two `S2`
  objects, half-strength gyroscopic terms, `collect` false identity and missing
  cross-inertia, IBP prefactor, `Hat(Vee(X))`, silent extraction zeros, duplicate entries,
  tangent vector as variable, composite entries, `Number(0)`, float residue,
  `_rebuild`, `full_simplify` cap, bool, aliased `attr`, reserved printed names) and the
  ones that remain (`Matrix("J1")` twice, unreduced terms under `Dot`, tangent-plane
  solve, negative `M`); the two `AlgebraicError` causes with their full text; the exact
  text of every typed error a user can hit (building expressions, `compute_eom`,
  `to_standard_form`, `linearize`, engine functions, `check_tree`), produced by a script
  that triggers each; the tools (trace, `GEOMECH_DEBUG`, `GEOMECH_TRACE`, `print_tree`,
  `check_tree`, `check_eom`, the importable submodules).
- `docs/conventions.md` (new): `E = −(Euler–Lagrange) + forcing`, `M = −(inertia)`, VBL
  sign, S² tangent plane, row forms (`(v)'s` vs `(vs)'`, inner/outer, the row kinetic
  energy, rows rejected by `check_tree`), zeros (`0`, `0v`, `(0)`, `(0v)'`, and the constant
  vector/matrix variations printing `0`/`O`), `0 + x` vs `x + 0` (`x + 0` is an error for
  a vector), numeric folding, the `1e-12` cancellation tolerance, a named `(-1)`, `str` is
  the identity, glossary (configuration variable, tangent vector, variation vector,
  perturbation, params, inputs, constant).
- `mkdocs.yml` nav: Pipeline, Conventions, Troubleshooting, and
  `manifold-implementation-notes.md` as "Manifold Rules: Derivation". I kept that page
  (the contributor notes tell readers to read it before changing the S²/SO(3) rules, and I
  may not edit them); its stale rule lists (one listed `Vee(a + b) → Vee(a) + Vee(b)`, which
  never fires) were replaced by the named rules that depend on the manifold structure.
- Same pipeline order in `README.md` ("The pipeline" step 3 and "How it works", with the
  `_eliminate` snippet refreshed: it showed the pre-phase-1 `Hat(Vee(M)) = M`),
  `docs/getting-started.md` ("How it works", six numbered stages), `docs/index.md`, and the
  `eom.py` module and `compute_eom` docstrings; all link to `docs/pipeline.md`.
- Stale text the reviews listed: `variables.py:12` ("not yet used in EOM pipeline") and the
  `extract.py` header (`(query, expr)` signature, "Assumes query appears linearly") were
  already fixed in phase 1/2 (`865699f`, `feeb476`); nothing left to change.
- **`paper_context.md`**: `.gitignore:7` ignores it while the contributor notes cite it as a
  design document next to `docs/manifold-implementation-notes.md`. A fresh clone has no such
  file. Not changed (per instructions); the author should either track it or drop the
  reference.

## 4.x `full_simplify` cap warning (phase-2 F-6)

`geomech/core/transformations/simplify.py:131`: "full_simplify: no confirmed fixed point
after N passes (the last pass still changed the expression, so no fixed point after N
passes was confirmed by an unchanged pass); result may be unsimplified: <repr>", and for
`max_iter=0` "(max_iter=0: no simplify pass ran)". The parenthetical is longer than I would
write from scratch: the existing `tests/test_validation.py:604` matches
`"no fixed point after 1 passes"`, which the requested wording alone does not contain, and
existing tests are not edited; the explanation keeps that substring.

## Scenario walkthroughs again (debuggability audit (i))

### (a) The user gets an exception

The audit's A-table rows A3–A33 were addressed in phases 2–3. The phase-4 cases, run on
`57232da` (before) and on this tree (after):

| Case | Before (`57232da`) | After |
|---|---|---|
| C3: `L = w·Vee(Hat(Ω) A)` | `AlgebraicError: extract_linear_coeff: '\eta_{R}' appears only inside the calculus operator '\frac{d}{dt}(\eta_{R})' (a higher order derivative, a variation or an integral of the query); such a term has no linear coefficient` | `AlgebraicError: compute_eom: while extracting the equation for '\eta_{R}': the time derivative '\frac{d}{dt}(\eta_{R})' of the variation vector '\eta_{R}' survived integration by parts in the term 'Dot(w,Vee((Hat(\frac{d}{dt}(\eta_{R}))A+Hat(Cross(\Omega_{R},\eta_{R}))A)))' of δS; collect could not expose it as the left operand of a Dot (Dot(\frac{d}{dt}(\eta_{R}), f)), so the term was not integrated by parts and has no coefficient. Inspect the term with print_tree(term, style='indent'), or run compute_eom(..., trace=EOMTrace()) and look at the stage 'ibp.collect[\frac{d}{dt}(\eta_{R})]'. (engine: extract_linear_coeff: …)` |
| A26: jerk in `to_standard_form` | `AlgebraicError: extract_linear_coeff: '\frac{d}{dt}(\frac{d}{dt}(xs))' appears only inside the calculus operator '\frac{d}{dt}(\frac{d}{dt}(\frac{d}{dt}(xs)))' …` (same shape as C3) | `AlgebraicError: to_standard_form: in equation '\delta{xs}', term '\frac{d}{dt}(\frac{d}{dt}(\frac{d}{dt}(xs)))': the acceleration '\frac{d}{dt}(\frac{d}{dt}(xs))' occurs only inside '…', a third or higher time derivative of the configuration (a jerk): the equation is not second order, so it has no standard form M a + f + G u = 0 (engine: …)` |
| `SystemVariables(vectors=[x + y])` | `ValueError: compute_eom: δS still contains the variation symbols '\delta{x}', '\delta{y}' of the symbols 'x', 'y', which are neither entries of SystemVariables nor declared in params; …` | `TypeError: compute_eom: SystemVariables.vectors entry '(x+y)' is a VAdd expression, not a single symbol; a configuration variable is one Scalar, Vector, S2 or SO3 leaf (list the symbols it is built from instead)` |
| `linearize` on a params equation | `… the linearization is incomplete. Is an input missing from \`inputs\`, or a symbol missing the 'Constant' attribute?` (both wrong for a param) | without `params`: `… It carries '\delta{l}' (of 'l'), which is neither a configuration, velocity nor input variation. If the symbol is an input, list it in \`inputs\`; if it is a constant, give it attr=['Constant']; if it is a prescribed parameter (compute_eom(..., params=[...])), linearize does not support prescribed parameters yet (pass params=[...] to linearize to have it say so).`; with `params=[l]`: `… It carries the variation of the prescribed parameter(s) 'l' declared in params: linearize does not support prescribed parameters yet (…)` |
| `L = ½ ẋᵀAẋ − ½ xᵀKx` | `AlgebraicError: extract_linear_coeff: '\delta{x}' appears only inside the calculus operator '\frac{d}{dt}(\delta{x})' …` | derives: `(A\frac{d}{dt}(\frac{d}{dt}(x))(-0.5)+(A)'\frac{d}{dt}(\frac{d}{dt}(x))(-0.5)+Kx(-1.0)+f)` |
| `½ Ω·JΩ + c₁·((Ω c₂ᵀ) c₃)` | `AlgebraicError: collect: the target … occurs inside the matrix factor …` (no pipeline stage) | `AlgebraicError: compute_eom: while integrating by parts (stage ibp, targets \frac{d}{dt}(\eta_{R})): collect: the target …` |

What the debugger sees now: every pipeline error says which function, which stage, which
variable and which term, and the two "appears only inside the calculus operator" causes
are distinguishable from the first words. Each error chains the engine exception
(`__cause__`), so the traceback still ends at the engine frame.

### (b) The equation is wrong: can I look inside?

Before (audit): eight hand-written steps, one function not exported
(`apply_manifold_rules`), one private (`_apply_ibp`) to see the IBP pass without the
surrounding `expand`/`full_simplify`/`collect`, and the knowledge that
`integrate_by_parts` re-expands internally. No `trace=`, logger or stages object.

After:

```python
tr = EOMTrace()
eom = compute_eom(L, dW, vars_, trace=tr)
tr.summary()          # 13 lines for the rigid pendulum, see docs/pipeline.md
tr.show("ibp")        # print_tree(style="indent")
tr.ibp["\\frac{d}{dt}(\\eta_{R})"]["collect"]   # collect alone, IBP pass alone
```

`apply_manifold_rules` is exported; every stage is reproducible from the previous one
with public names (`tests/test_trace.py`). On failure the trace stops where the pipeline
stopped. For C3 the summary shows at a glance that `ibp.collect`, `ibp.integrate` and
`ibp.simplify` all equal `ibp.prepare` (15 nodes each, same repr): collect and the IBP
pass did nothing to the `Vee` term:

```
simplify                                  15 nodes  <v:w*, vee((hat(d/dt(v:\eta_{R}))·M:A* + hat(cross(v:\Omega_{R}, v:\eta_{R}))·M:A*))>
ibp.prepare                               15 nodes  <v:w*, vee((hat(d/dt(v:\eta_{R}))·M:A* + …))>
ibp.collect[\frac{d}{dt}(\eta_{R})]       15 nodes  <v:w*, vee((hat(d/dt(v:\eta_{R}))·M:A* + …))>
ibp.integrate[\frac{d}{dt}(\eta_{R})]     15 nodes  <v:w*, vee((hat(d/dt(v:\eta_{R}))·M:A* + …))>
ibp.simplify[\frac{d}{dt}(\eta_{R})]      15 nodes  …
ibp                                       15 nodes  …
expand                                    15 nodes  …
```

Without changing code: `GEOMECH_DEBUG=1 python script.py` prints the same stages on
stderr (`geomech: compute_eom: dS [33 nodes] …`).

### (c) A rule misfires

Before: no rule names, no trace, and `import geomech.core.transformations.simplify as S`
bound the *function*.

After, with `GEOMECH_TRACE=1` on the audit's C1/C6 probes:

```
geomech.rules: vector.s2_tangent_recovery: cross(S2:q, cross(v:\omega_{q}, S2:q)) ==> v:\omega_{q}
top level: \omega_{q}
under Dot: Dot(a,Cross(q,Cross(\omega_{q},q)))          ← no line logged: no rule fired
geomech.rules: vector.dot_tangent_orthogonal: <v:\omega_{q}, S2:q> ==> 0
second S2 object: 0
```

The silence for `Dot(a, q × (ω × q))` is now the diagnosis, and the `vector_rules`
docstring says why next to the rule: the pass recurses only through `Add`, `Mul`, `VAdd`,
`SVMul` and `Cross` (C1 itself is phase 5.1). C6 (a second `S2("q")`) was fixed in phase 2;
the trace shows the rule firing for it. `import geomech.core.transformations.vector_rules
as V` gives the module (`CallableModule`), so a breakpoint or a patch of a helper works
without `sys.modules`.

## Follow-ups and notes

- Phase 5.1 (C1) is unchanged: `vector_rules` and `pull` do not descend into `Dot`,
  `MVMul`, `Hat`, `Transpose` or calculus nodes. The trace makes it visible, it does not
  fix it.
- The rows in `Dot`/`Cross`/`Hat` and mixed `VAdd` are rejected by `check_tree` only (see
  4.4 for the two blocking tests and the 4 % constructor cost). If the author wants them
  rejected at construction, the two `test_geometry.py` tests and the
  `test_latex_precedence.py:50` parameter would have to change, and `Dot`'s check needs a
  cheaper orientation test.
- `Vector('c', Constant).delta()` prints `0` (a zero vector leaf), not `ZeroVector`'s `0v`,
  and a constant matrix's `δ` prints `O`; documented in Conventions, not changed (it would
  change printed output; phase 5.5's `zero_like`).
- `expand.vee_over_vadd` is dead code (Vee's operand is a matrix); named and documented as
  unreachable so the rule-name test covers it; phase 5.2 replaces it with Vee over `MAdd`.
- `_extract_single` in `standard_form.py` gained a `key` argument (private).
- `linearize` has a new optional `params` argument that only improves the error message;
  supporting params in VBL remains open (phase-2 F-1, "phase 4/5").
