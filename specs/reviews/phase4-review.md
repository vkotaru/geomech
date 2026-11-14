# Phase 4 adversarial review: debuggability

Branch `fix/phase4-debuggability`, commit `bffba73`, range `57232da..bffba73`. Reviewed
against `specs/fix-plan.md` phase 4 (4.1–4.7), `specs/reviews/debuggability-audit.md` and
`specs/reviews/phase4-implementation-report.md`. Added: `tests/test_phase4_adversarial.py`
(94 passing tests, 9 strict xfails F1–F9) and this file. Nothing else was changed and
nothing is committed.

The question asked: will debugging this code in 6–12 months be easy?

**Verdict: GO WITH CHANGES.** The stage trace, the named rules and the docs did their job
in a real test: three blind debuggers each found an injected engine bug in under a
minute. The changes below are small. Two are about trust, because the tooling must never
point the wrong way: the misdiagnosed second-order error (F1) and fold rules taking the
blame for a child's rewrite (F3). The other three are a trace that silently mixes two runs
(F2), lines printed twice under `logging.basicConfig` (F4), and the callable-module shim
(F5).

## 1. Bug hunts (main deliverable)

Each bug was injected into a plain copy of `bffba73` (no git history) inside the
temporary directory. Each copy got a `scenario.py` whose `check_eom` fails. A *fresh* agent, which
knew nothing about the bug or the code, then debugged it. It was told to use only the
phase-4 tooling (`EOMTrace`, `GEOMECH_DEBUG`, `GEOMECH_TRACE`/`set_trace`, `check_tree`,
`print_tree`, `pipeline.md`, `troubleshooting.md`, the rule docstrings), and to read code
only after the tooling had pointed at a suspect. Running the test suite, reading tests and
diffing against other copies were forbidden. The copies and worktrees have been removed.

### (i) Sign of one S² substitution term

- **Injected:** `manifold_rules._NEG2 = Scalar("2", value=2)`, so `δω → ξ̇ + 2ω×ξ`.
- **System:** spherical pendulum with a gyroscopic term `c·ω`. `check_eom` = 2.8 (correct
  code: 5.9e-9). Like the other two bugs, it also fails 15+ pinned tests and is caught
  by the oracle.
- **Path:**
  1. Read and ran the scenario. A hand derivation showed the `ω×c` term with the wrong
     sign.
  2. Read `pipeline.md`. Stage 2 states `δω → ξ̇ − 2ω×ξ`.
  3. Narrowed L to `c·ω` and ran `tr.summary(chars=300)`. The `manifold_rules` line reads
     `(d/dt(ξ) + cross(ω, ξ)·2)`, which contradicts the doc.
  4. Ran `set_trace(True)` with `apply_manifold_rules`. The log line
     `manifold.s2_variation_omega: δ(ω) ==> (d/dt(ξ) + cross(ω, ξ)·2)` names the rule.
  5. Read lines 52–110. Found `_NEG2` with `value=2`.
- **Cost:** 7 tool steps, about 1 minute.
- **What helped:** `summary()` showed the first stage that disagrees with the documented
  formula. The rule name plus the formula in the docstring turned "stage 2 is wrong" into
  "this rule disagrees with its spec" without reading code.
- **What failed:**
  - The troubleshooting row "half-strength gyroscopic terms (*fixed*)" is a near miss
    that could reassure: the output does contain a 2.
  - The docs do not suggest bisecting the Lagrangian term by term. The agent thought of
    it alone.
  - Nothing checks a rule against its documented identity.

### (ii) A wrong `collect` rotation

- **Injected:** `collect._collect_dot`, case `(_, Cross())`, `r.left == vec` returned
  `Dot(vec, Cross(l, r.right))` (a·(v×c) → v·(a×c), a sign flip). The comment above it
  still states the right identity.
- **System:** rigid pendulum about the centre of mass. `check_eom` = 2.43 (correct: 2.9e-8).
- **Path:**
  1. Ran the scenario. A hand derivation found one wrong term, `−m(Ω̇×ρ)×ρ`.
  2. Printed `tr.ibp[t]["collect"]`. The term `Dot(R(η̇×ρ), R(Ω×ρ))m` became
     `Dot(η̇, (Ω×ρ)×ρ)m`, which has the wrong sign, at stage `ibp.collect`.
  3. Ran `collect` on that term alone with `set_trace(True)`. The log showed
     `collect.strip_rotation → collect.flip → collect.triple_product: <Ω×ρ, η̇×ρ> ==> <η̇, (Ω×ρ)×ρ>`.
  4. Grepped the rule name: 8 sites.
  5. Read the code and found the line where code and comment disagree.
- **Cost:** 6 steps, about 35 s.
- **What helped:** the per-target `ibp.collect`/`integrate`/`simplify` sub-stages, which
  were decisive, and calling `collect` by hand.
- **What failed:**
  - `collect.triple_product` names 8 call sites, so the name localises to a family, not a
    line.
  - Under `logging.basicConfig()` every rule line printed twice (F4).
  - `pipeline.md`'s sample summary is for a different L. That is harmless.

### (iii) A vector rule firing without its hypothesis

- **Injected:** `_is_tangent_orthogonal` lost its parent check (`a.S2 == b`), so
  `Dot(ω_{q1}, q2) → 0`.
- **System:** two spheres with `½k(ω₁·q₂)²`. The whole coupling vanishes. `check_eom` = 1.08
  (correct: 9.7e-8). No oracle system catches this; only three unit tests do.
- **Path:**
  1. Ran the scenario.
  2. Read `pipeline.md`.
  3. Ran `summary()` and grepped each stage for `k`. It is present through
     `manifold_rules` and gone at `simplify` (114 → 33 nodes).
  4. Reran `full_simplify(tr["manifold_rules"])` under `set_trace`. 99 lines, including
     `vector.dot_tangent_orthogonal: <ω_{q1}, S2:q2> ==> 0` (false on sight).
  5. Read the docstring hypothesis ("point equals q by name") and the helper. The check
     is missing.
- **Cost:** 9 steps, about 45 s.
- **What helped:** the node-count drop in `summary()`, and `before ==> after` with the
  hypothesis next to the name.
- **What failed:**
  - The rule log is noisy: 36 `simplify.mul_fold` lines buried 20 `vector.*` lines. Part
    of this is F3, part is that one rewrite is logged on every pass and for every copy
    of the subterm.
  - The trace switch is global, so the stage had to be rerun by hand to scope the log.
  - `troubleshooting.md` has no row for "a term vanished". Its "compared by name
    (*fixed*)" row points away from the bug.

### Summary of the hunts

| Bug | Steps | Time | Localised by | `check_tree` / `GEOMECH_DEBUG` |
|---|---|---|---|---|
| (i) S² sign | 7 | ~1 min | `summary()` stage vs doc formula, then rule trace | silent (semantic bug) |
| (ii) collect rotation | 6 | ~35 s | `tr.ibp[t]["collect"]`, then rule trace | silent |
| (iii) vector rule hypothesis | 9 | ~45 s | node-count drop in `summary()`, then rule trace | silent |

All three hunts needed a hand derivation to know which stage was *wrong*. The trace shows
*where* each change happens, not which change is false. The missing piece is already in
the repo but not in the debugging docs. `specs/reviews/scripts/engine_fuzz.py --quick`
on each bugged copy reported `[wrong-value]` findings with shrunk counterexamples within
seconds:

- `manifold_rules[full delta]` for (i);
- `collect[target occurs once]` for (ii);
- `vector_rules: z_v Dot(q, ω_p) → z_v 0` for (iii).

No human mathematics was needed. All three agents asked for exactly this ("a numeric
spot-check per rule firing"). See F7.

## 2. Findings, ranked

The order is: tooling that misleads, then tooling that fails to localise, then runtime
cost, then docs mismatch, then cosmetic. Each finding has a strict xfail in
`tests/test_phase4_adversarial.py`, except F6 and F7, which are recommendations.

### Tooling that misleads

**F1. A Lagrangian in accelerations is reported as a `collect` failure, with a stage name
that does not exist.**
- **Where:** `geomech/dynamics/eom.py:540-569` (`_extraction_error`).
- **Reproducer:** `L = ½ ẍ·ẍ`, `SystemVariables(vectors=[x])`.
- **What it says:** "the time derivative `d/dt(d/dt(δx))` … survived integration by parts
  in the term `Dot(d/dt(d/dt(δx)), ẍ)` …; collect could not expose it as the left operand
  of a Dot". It then tells the user to look at stage `ibp.collect[\frac{d}{dt}(\frac{d}{dt}(\delta{x}))]`.
- **Why it misleads:**
  - The derivative *is* already the left operand.
  - The real cause: IBP runs once per first derivative, so a second derivative of a
    variation (L depends on accelerations, a higher-order Lagrangian) is never removed.
  - The named stage does not exist. The stages are keyed by the first derivative
    `ibp.collect[\frac{d}{dt}(\delta{x})]`.
  - The scalar variant says "the left operand of a Dot" for a `Mul`.
- **Fix:**
  - When the surviving derivative is of order ≥ 2, say "L depends on the acceleration of
    x (higher-order Lagrangian); not supported".
  - Build the stage name from the first-order target.
- **Test:** `test_second_order_lagrangian_error_names_the_real_cause`.

**F2. A reused `EOMTrace`/`LinearizationTrace` silently accumulates.**
- **Where:** `trace.py:171-185` (`_Trace._add` appends; nothing resets).
- **Reproducer:**
  ```python
  tr = EOMTrace()
  compute_eom(*point_mass, trace=tr)
  compute_eom(*rigid_pendulum, trace=tr)
  ```
- **Result:** 26 stages in `summary()`, two runs interleaved. `tr.equations` holds `δx`
  from the first run and `η_R` from the second. `tr["dL"]` silently returns the second
  run's value. `linearize` behaves the same: 10 → 20 stages.
- **Why it matters:** a notebook user who re-executes a cell gets a trace that is wrong
  without any sign of it.
- **Fix:** reset the trace at pipeline entry, or raise "trace already used".
- **Tests:** `test_reused_eom_trace_is_reset_or_rejected`,
  `test_reused_linearization_trace_is_reset_or_rejected`.

**F3. Fold rules are logged for a change made by a child.**
- **Where:** `simplify.py:149-160` (`rule("simplify.add_fold", expr, …)`, `mul_fold`,
  `vadd_fold`, `madd_fold`). *before* is the node before its children were eliminated.
- **Reproducer:** `_eliminate(k + m·Dot(a, Hat(w)·v))` logs three lines:
  - `simplify.mvmul_hat_to_cross`;
  - `simplify.mul_fold: m·<a, hat(w)·v> ==> m·<a, cross(w, v)>`;
  - `simplify.add_fold: (k + …hat…) ==> (k + …cross…)`.
- **Why it matters:**
  - Every real firing is echoed once per enclosing `Add`/`Mul`. The two fold lines claim
    rewrites that did not happen.
  - The firing counts in `pipeline.md` ("`simplify.mul_fold`", 7 lines for the rigid
    pendulum) are inflated by this.
  - In hunt (iii) this noise was 36 of 99 lines.
- **Fix:** pass the node rebuilt from the eliminated children as *before*
  (`Add(*elim_children)`), so a fold logs only its own work.
- **Test:** `test_fold_rules_are_not_logged_for_a_child_change`.

### Tooling that fails to localise

**F4. Every `geomech` line is printed twice when the application uses
`logging.basicConfig()`.**
- **Where:** `trace.py:53-59` (`_ensure_handler`) attaches a `stderr` handler to
  `geomech`, and leaves `propagate=True`.
- **When it happens:**
  - With `set_trace(True)` after `basicConfig()`. This is exactly the recipe that
    `pipeline.md` and `troubleshooting.md` give for notebooks.
  - With `GEOMECH_DEBUG=1` plus a root handler.
- **Effect:** each line appears once from `geomech` and once more through root.
  Debugger (ii) read this as "the rule fires twice".
- **Fix:** do not add a handler when the root logger has one, or set `propagate=False`
  on the handler-owning logger (the library-logging convention is a `NullHandler` only).
- **Test:** `test_set_trace_does_not_duplicate_lines_under_basic_config`.

**F5. Callable modules (4.6) break pickling, `inspect` and `functools.wraps`.**
- **Where:** `geomech/utils/callable_module.py`, `core/transformations/__init__.py:13-21`,
  `core/math/__init__.py`.
- **Effect of `from geomech.core.transformations import simplify`** (it used to be a
  function, now it is a `CallableModule`):
  - `pickle.dumps(simplify)` / `copy.deepcopy` → `TypeError: cannot pickle
    'CallableModule'`;
  - `ProcessPoolExecutor().map(simplify, …)` fails;
  - `inspect.signature` → `(*args, **kwargs)`;
  - `isinstance(simplify, FunctionType)` and `inspect.isfunction` are `False`;
    `inspect.ismodule` is `True`;
  - `help(simplify)` and `__doc__` show the module docstring;
  - `functools.wraps(simplify)` names the wrapper `geomech.core.transformations.simplify`.
- **What still works:** mkdocs (`--strict` builds; griffe is static), pickling
  `geomech.simplify`, and monkeypatching `module.simplify`.
- **Why it matters:** two names that look alike (`geomech.simplify` and
  `geomech.core.transformations.simplify`) now have different kinds. That is exactly the
  surprise item 4.6 meant to remove.
- **Recommendation:** rename the four modules (and `collect`) to names that differ from
  their functions, and keep the package exporting the *functions* under the old names:
  - `expand` → `expansion`, `simplify` → `simplification`, `pull` → `pulling`,
    `vector_rules` → `vector_identities`, `collect` → `collecting`;
  - or prefix the modules with `_`.
  - Then nothing is shadowed and nothing needs a `__class__` swap. The four tests that
    import from the package keep working unchanged. A breakpoint goes in
    `geomech.core.transformations.simplification`.
  - The report says "renaming would make `geomech.core.transformations.simplify` the
    function, which is exactly what the item forbids". That is a misreading: the item
    forbids a function *shadowing its module*, and with distinct module names there is
    no shadowing.
  - Delete `callable_module.py`.
- **Tests:** `test_package_simplify_is_picklable`,
  `test_package_simplify_has_the_function_signature`.

**F6. The rule-name sync test cannot see a rewrite that bypasses `rule()`.**
- **Where:** `tests/test_rule_names.py:65-90`.
- **What I checked:** adding `rule("vector.transpose_zero", …)` without a docstring
  entry fails the test, which is good (`test_sync_test_catches_a_named_rule_missing_from_the_docstring`).
  Adding `case Transpose() if expr.expr.is_zero: return expr.expr` with no `rule()` passes
  it (`test_sync_test_is_blind_to_a_rewrite_that_bypasses_rule`).
- **Why it matters:** the docstring is the index of the rules only as long as everyone
  remembers to call `rule()`.
- **Cheap guard:** the fuzz harness already rewrites thousands of trees. Run it with
  `set_trace(True)` and assert that every output differing from its input has at least
  one logged rule line, for each transformation.
- **Status:** recommendation only, no xfail.

**F7. The numeric localiser is not part of the debugging story.**
- **What:** `specs/reviews/scripts/engine_fuzz.py` found all three injected bugs in
  seconds with minimal counterexamples. `pipeline.md` and `troubleshooting.md` never
  mention it, and `GEOMECH_DEBUG` checks structure only, which was silent for all three
  bugs.
- **Recommendation:**
  - Add a "the equation is wrong and no row matches" recipe to `troubleshooting.md`:
    1. bisect with `tr.summary()`;
    2. run `engine_fuzz.py --quick`;
    3. rerun the suspect stage under `set_trace`.
  - Consider moving the fuzzer's jet evaluator into `geomech.numeric` as
    `check_rewrite(before, after)`.
- **Status:** recommendation only, no xfail.

### Runtime cost

None found.
- **`compute_eom`, median of 40, three alternating rounds, same machine, idle:**

  | System | `57232da` | `bffba73` |
  |---|---|---|
  | Quadrotor with load | 11.07 / 11.10 / 11.04 ms | 10.87 / 10.91 / 10.83 ms |
  | Double rigid pendulum | 17.94 / 18.26 / 18.05 ms | 17.87 / 17.38 / 17.28 ms |

- **`rule()` with tracing off:** one Python call that reads one global, about 40 ns. It
  is called 888 times per quadrotor `compute_eom`. `_eliminate` calls it on every
  `Add`/`Mul` visit, not only on rewrites. That totals ≈ 35 µs, ≈ 0.3 %, below noise.
- **Lazy formatting:**
  - `stage_recorder` returns `None` unless a trace is given or `geomech` is at `DEBUG`,
    and every recording site is guarded.
  - Logger calls use %-formatting; an AST test asserts there are no f-strings.
  - `rule()` builds no `repr` when off (test).
- **Cosmetic:** with tracing on but the `geomech.rules` logger raised above `DEBUG`,
  `rule()` still builds both `repr`s (`trace.py:99`). It could guard on `isEnabledFor`.

### Docs mismatch

**F8. The `SystemVariables()` row cannot be reached with a real Lagrangian.**
- **Where:** `eom.py:500`; the check runs last, after the stray-symbol check.
- **What the user gets:** `troubleshooting.md` says the empty set gives "no configuration
  variables listed". Any Lagrangian that holds a symbol actually gives the stray message
  (`δS still contains '\delta{x}'…`). The documented text appears only for a constant L.
- **Fix:** move the check up, or fix the row.
- **Test:** `test_empty_system_variables_message_for_a_real_lagrangian`.

**F9. `integrate_by_parts` does not log.** `trace.py:6` (module docstring) and the report
say it logs "one DEBUG line per stage of … `integrate_by_parts` (per target)". The public
function passes `record=None` and logs nothing; only the `compute_eom` path logs. Test:
`test_public_integrate_by_parts_logs_its_stages`.

### Cosmetic

- **Troubleshooting coverage.** The table covers every silent case of phase 1 (S² rule,
  `collect` identity, IBP prefactor, `Hat(Vee)`, extraction zero) and phase 2 (2.1–2.12).
  Phase 3's silent cases are missing:
  - 3.1/3.2: `Transpose(v) * s` built a double transpose, and the row encoding was lost
    under `expand`;
  - 3.5: ambiguous LaTeX such as `a × b × c` and `(b·c) a`, which made a correct
    equation print as a wrong one.

  All three bug hunts also found no row for "wrong sign" or "a term vanished".
- **The `(stage expand)` label.** The nonlinear-extraction message says "(stage expand)"
  for a failure in extraction (`eom.py:568`).
- **`tr["ibp"]` vs `tr.ibp`.** `tr["ibp"]` is `after_ibp`, while the attribute `tr.ibp` is
  the per-target dict. `show("ibp")` therefore cannot show the per-target steps.
- **No failure marker in a trace.** After an exception the summary just stops. A final
  `failed: <ExcType>` line would make "stopped here" explicit. Similarly, a trace from a
  run that failed in the stray check looks like a success.
- **`GEOMECH_DEBUG` checks only the pipeline's inputs and outputs.** Running `check_tree`
  inside `record()` in debug mode would check every stage at no extra code.
- **`expand.vee_over_vadd` is unreachable.** `Vee(a + b)` is rejected by `Vee`'s
  constructor (test). It should be **deleted**, not documented: a rule index that lists
  dead rules teaches the reader to distrust the index. Phase 5.2 adds `vee_over_madd`.
- **Rule names cover families, not lines.** `collect.triple_product` and `collect.flip`
  name 8 and 6 call sites. Adding the caller's line to the trace line (`stacklevel`) or
  splitting the names (`triple_product.left`, `.right`) would finish the localisation
  that hunt (ii) had to do by grep.

## 3. What held

- **Trace fidelity.** For the point mass, rigid pendulum, spherical pendulum and
  quadrotor, `after_simplify`, `after_ibp`, `after_expand` and every `equations[k]` are
  the *same objects* (`is`) that the next pipeline function received or that
  `compute_eom` returned. They are not recomputations; I checked this by spying on the
  calls.
- **Trace after a failure.**
  - After an extraction failure, the trace ends at `expand`, and `raw_equations` and
    `equations` are `{}`.
  - After a validation failure, the trace is empty.
  - Stages are recorded before the next stage runs.
- **Error chaining.**
  - Every re-raise keeps the exception type and sets `__cause__` (`AlgebraicError`,
    `NonLinearError`).
  - The engine text appears exactly once, and `compute_eom:` appears exactly once.
  - The pipelines never nest (`linearize` does not call `compute_eom`), so no message is
    duplicated.
  - A term with 60 addends is cut: the message is under 2000 characters.
  - The nonlinear-`dW` message names the real cause ("appears in both sides of Mul").
- **Named rules.**
  - The sync test catches a named rule missing from the docstring.
  - I sampled ten rules and checked each against the hypothesis next to it, firing and
    not firing: `dot_tangent_orthogonal`, `s2_tangent_recovery`, `dot_unit_norm`,
    `transpose_symmetric`, `orthogonal_transpose`, `hat_vee_skew`, `cross_self`,
    `cross_scaled_self`, `mvmul_hat_to_cross`, `s2_double_cross`. All ten match the
    code.
  - `pipeline.md`'s list of rules that fire for the rigid pendulum matches the code.
- **`check_tree`.**
  - Each invariant I violated raises and names the path: nested sum, wrong addend type,
    row in `Cross`/`Hat`, mixed `VAdd`, `TS2` without a point, `(-1)` without a value,
    `d/dt` over an `S2` at `root.R.R.R`.
  - No false positive on 18 000 fuzz trees: 3000 `engine_fuzz.Gen` trees, seed 4, each
    with its `expand`, `pull`, `vector_rules`, `simplify` and `full_simplify` output.
  - The autouse fixture really checks: a Lagrangian with a row in a `Dot` derives
    silently (`(ẍ(-1)+(a)')`), but the patched `compute_eom` raises `InvariantError`
    in the suite.
- **Docs.**
  - **61 of 63** typed-error rows of `troubleshooting.md` reproduce **exactly**,
    compared as type plus message. The other two are the `full_simplify` warning, whose
    prefix I checked, and the `linearize(params=…)` row.
  - Both long `AlgebraicError` blocks (IBP survivor, jerk) match.
  - `pipeline.md`: the `summary(chars=70)` output for the rigid pendulum reproduces
    byte for byte. Its snippets run, and its stage names, rule names and module table
    are accurate.
  - `conventions.md` has no code blocks.
- **Regression.**

  | Check | Result |
  |---|---|
  | Suite, `-W error` | 2105 passed, 14 xfailed (2011 + 94 new; 5 + 9 xfails), 8.3 s |
  | Suite, `GEOMECH_DEBUG=1` | same, 9.1 s |
  | Suite, `GEOMECH_TRACE=1` | same, 10.6 s. With `-s` the rule log is 5.2 MB (31 600 lines); captured, so usable, but a failing test under `GEOMECH_TRACE` prints its share |
  | `engine_fuzz.py --cases 2000 --seed 4` | exit 0, 0 `wrong-value` findings |
  | Examples (6), `-W error` | byte-identical to `57232da` |
  | Notebooks (4), executed | every code-cell output identical to the committed ones |
  | `str` dump of the 65 systems vs `57232da` | `compute_eom`, `to_standard_form` (every field) and `linearize` (M, A, B) identical; only the two expected `AlgebraicError` prefixes differ |
  | `ruff format --check`, `ruff check` (`geomech/`, `tests/`) | clean |
  | `mkdocs build --strict` | passes |

## 4. Changes before merge

1. **F1:** order-aware message and the correct stage name in `_extraction_error`.
2. **F3:** fold rules log against the rebuilt node.
3. **F2:** reset or reject a used trace.
4. **F4:** no handler when root has one, or `propagate=False`.
5. **F5:** replace `CallableModule` with distinct module names.
6. **F8, F9:** fix the docs. Delete `expand.vee_over_vadd`. Add a "no row matches"
   recipe that names `engine_fuzz.py` (F7).

F6 and the cosmetic items can go into phase 5.

## Fix log

Fixes on `fix/phase4-debuggability` on top of `6239da7`, applying the decisions taken on
this review. Nothing is committed.

### Findings

| Finding | Fix | xfail |
|---|---|---|
| F1 | `_extraction_error` (`geomech/dynamics/eom.py`) diagnoses before it explains. A surviving derivative of order ≥ 2 on a variation, or an acceleration of a configuration variable in `L` (`ẍ`, `ω̇`, `Ω̇`), gives "the Lagrangian depends on accelerations (…); higher-order Lagrangians are not supported", pointing to the recorded stage `dL`, with no mention of `collect`. Only a first-order derivative that is not the left operand of its `Dot` (for a scalar variation: not a direct factor of its product) gets "collect could not expose it", pointing to `ibp.collect[<first-order target>]`; one already in that shape points to `ibp.collect[t]`/`ibp.integrate[t]`. Stage names are built from `str(v.t_diff())`, the key the IBP stages use. `compute_eom` keeps calling `_separate_variations(dS, vecs, record)` with three arguments (the trace-fidelity spy pins that call) and retells the error with `L` when `L` holds an acceleration. The `(stage expand)` label of other extraction errors is now `(extraction, after the stage 'expand')`. | flipped |
| F2 | `_Trace._reset()` clears every field at the entry of `compute_eom` and `linearize`; new fields `failed_at: str \| None` and `error: str \| None`, set when the pipeline raises. `failed_at` is the stage that raised (`validation`, `dL`, `manifold_rules`, `simplify`, `ibp.prepare`/`ibp.collect[t]`/…, `expand`, `raw[ξ]`/`equation[ξ]`, `stray_check`; for `linearize` `dE[row]`, `bucket[row]`, `raw[row][sym]`, `block[row][sym]`). `summary()` starts with `FAILED at <stage>: <exception> (the stages below ran before it)`. Documented in `trace.py`, the `compute_eom`/`linearize` docstrings, `pipeline.md`, `troubleshooting.md`. | both flipped |
| F3 | The four fold rules log against `cls(*eliminated children)` (`_folded_from`), built only while tracing is on (`if _trace._TRACE_RULES:` before the `rule(...)` call, so the happy path is one module attribute read instead of a `rule()` call). Rigid pendulum: 42 → 38 lines (`mul_fold` 7 → 4, `vadd_fold` 2 → 1); quadrotor with load: 261 → 248. `pipeline.md` now lists the per-rule counts of the rigid pendulum; its `summary(chars=70)` block is unchanged (re-run, byte-identical). | flipped |
| F4 | `_ensure_handler` adds the `stderr` handler only if `log.hasHandlers()` is false (the logger, `geomech` and the root logger). With `basicConfig(level=DEBUG)` first, each line prints once (subprocess test; also checked by hand with `GEOMECH_DEBUG=1`). Docs say to configure logging before `set_trace(True)`. | flipped |
| F5 | Modules renamed: `expand.py` → `expansion.py`, `pull.py` → `pulling.py`, `simplify.py` → `simplification.py`, `vector_rules.py` → `vector_identities.py`, `core/math/collect.py` → `collecting.py`. The packages export the functions under the old names; `callable_module.py` is deleted. `pickle.dumps(simplify)` works, `inspect.signature` is `(expr)`, `import geomech.core.transformations.simplification as m` gives the module, mkdocs generates `reference/geomech/core/transformations/{expansion,pulling,simplification,vector_identities}` and `core/math/collecting`. | both flipped |
| F6 | `tests/test_rule_coverage.py`: 1500 fuzz trees (seed 4) through `expand`, `pull`, `vector_rules`, `simplify`, `full_simplify`, and 1500 (seed 5) through `collect` on their first vector leaf, with tracing on; every output that prints differently from its input must have logged a rule line. No unnamed rewrite today. Checked that it catches one: replacing the `rule(...)` of `simplify.mvmul_hat_to_cross` by a bare return fails it. | — |
| F7 | `troubleshooting.md`: new section "The equation is wrong and no row matches" (confirm with `check_eom`; bisect with `summary()`; `engine_fuzz.py --quick`, look for `[wrong-value]`; rerun the suspect stage under `set_trace`) and a row pointing to it; `pipeline.md`: new section "When the equation is wrong" naming the fuzzer. The rows of fixed bugs moved to "History: silent bugs fixed in earlier versions", which says they are not suspects today. | — |
| F8 | Not flipped (see below). Docs fixed: the `SystemVariables()` row says the message appears only for an `L` without variables, and that a real `L` gets the stray-symbol message; the `compute_eom` docstring says the check runs last. | **kept** |
| F9 | The public `integrate_by_parts` passes `stage_recorder(None, "integrate_by_parts")`: one `DEBUG` line per stage when the `geomech` logger is at `DEBUG`, nothing otherwise. Docstring updated. | flipped |

### Cosmetic

- `expand.vee_over_vadd` deleted (code and docstring entry); the docstring says `Vee` is
  not distributed.
- The `(stage expand)` label is fixed (F1 row).
- `tr["ibp"]` vs `tr.ibp`: documented in the `EOMTrace` docstring and `pipeline.md`
  (stage `"ibp"` = `after_ibp`; attribute `ibp` = per-target dict; print one step with
  `show("ibp.collect[t]")`). Not renamed: `tests/test_trace.py` uses `tr.ibp[...]`.
  The examples in the docs now use `show("simplify")`.
- `collect.triple_product` split into eight names: `collect.triple_product_left`,
  `_right`, `_left_nested`, `_right_nested` and the `collect.triple_product_mvmul_*`
  variants of `Dot(MVMul, Cross)`; docstring and `pipeline.md` updated.
- History section: phase-3 silent cases 3.1/3.2 (double transpose of `Transpose(v) * s`;
  row encoding lost under `expand`) and 3.5 (ambiguous LaTeX) added.
- No failure marker: covered by `failed_at` (F2).
- `tests/test_troubleshooting_doc.py`: the `linearize(params=[l])` row, the new third
  `AlgebraicError` block (Lagrangian in accelerations), a completeness check (every
  typed-error row has a trigger in `TRIGGERS` or here, except the `full_simplify` warning
  checked by prefix elsewhere), the failed trace's stage names, and that every `(#anchor)`
  link resolves.

### Not done, needs a decision

1. **F8 cannot be flipped.** `tests/test_phase2_adversarial.py::
   test_empty_variables_with_a_real_lagrangian_raises_stray` pins the stray-symbol
   message for `compute_eom(point_mass L, dW, SystemVariables())`, and the F8 test wants
   the "no configuration variables" message for the same kind of call, identical to the
   doc row that the constant-`L` trigger also matches. Both cannot hold, so the check
   order is unchanged and the docs were fixed instead (as the decision allowed). The F8
   xfail stays; delete it, or retire the phase-2 test, to choose.
2. **`tests/test_rule_names.py::test_submodule_is_the_module` (5 cases) fails.** It pins
   the removed 4.6 mechanism: `getattr(package, "simplify") is module` and
   `callable(module)`. Under decision F5 the package attribute is the function, so these
   two assertions are wrong by design. Only its parametrize data (module paths) was
   updated; the assertions were not edited, per the instruction. The fix is to replace the
   two assertions with `getattr(package, func_name) is getattr(module, func_name)` (the
   rest of the test holds).

### Files changed

- Code: `geomech/__init__.py`, `geomech/core/transformations/__init__.py`,
  `geomech/core/math/__init__.py`, the five renamed modules, `geomech/core/math/ibp.py`,
  `geomech/dynamics/eom.py`, `geomech/dynamics/vbl.py`,
  `geomech/dynamics/standard_form.py` (imports), `geomech/utils/trace.py`;
  `geomech/utils/callable_module.py` deleted.
- Tests, import lines and module-path strings only: `test_collect.py` (also its module
  docstring), `test_collect_fix.py`, `test_collect_termination.py`, `test_expand.py`,
  `test_expand_vvmul.py`, `test_phase1_adversarial.py`, `test_pull.py`,
  `test_simplify.py`, `test_vector_rules.py`, `test_rule_names.py` (`RULE_MODULES` and
  the parametrize list), `test_phase4_adversarial.py` (imports, the source path and
  module name of the two sync-test checks, the import in the subprocess snippet, and the
  eight xfail markers removed). No assertion was edited. New: `test_rule_coverage.py`,
  `test_troubleshooting_doc.py`.
- Scripts: `specs/reviews/scripts/double_rigid_pendulum_numeric_check.py` (import).
  `engine_fuzz.py` needed no change.
- Docs: `troubleshooting.md`, `pipeline.md`, `conventions.md`,
  `manifold-implementation-notes.md`, `CHANGELOG.md`. Examples and notebooks import only
  from `geomech`; nothing to change. The historical review files under `specs/reviews/`
  keep the old module names on purpose.

### Verification

| Check | Result |
|---|---|
| Suite, `-W error` | 2115 passed, 6 xfailed (5 pre-existing + F8), **5 failed** (`test_submodule_is_the_module`, item 2 above) |
| Suite, `GEOMECH_DEBUG=1` / `GEOMECH_TRACE=1` | same counts, 9.8 s / 11.1 s |
| `ruff format --check`, `ruff check` (`geomech/`, `tests/`) | clean |
| `str` dump of the 65 systems | byte-identical to `6239da7`; against `57232da` only the two known `AlgebraicError` prefixes differ |
| `mkdocs build --strict` | exit 0 |
| `engine_fuzz.py --quick -W error` | exit 0, 0 `wrong-value`, the same 36 findings as `6239da7` (apart from `where:` file names) |
| Examples (6), `-W error` | exit 0, output byte-identical to `6239da7` |
| Notebooks (4), executed | every code-cell output identical to the committed ones |
| `compute_eom`, quadrotor with load, median of 40, alternating | `6239da7` 10.8–11.1 ms, now 10.9–11.1 ms (noise) |

### Arbitration (author decisions, 2026-11-13)

- `tests/test_rule_names.py::test_submodule_is_the_module`: the author chose to rewrite
  the assertion to the renamed-module design; it now asserts that the package exports the
  module's function under the old name.
- F8: the author kept the phase 2 behaviour (the stray-symbol message, which names the
  symbol to add, comes before the empty-SystemVariables check); the F8 xfail test was
  removed and the docs describe this order.
Suite: 2120 passed, 5 xfailed (pre-existing), also under GEOMECH_DEBUG=1.
