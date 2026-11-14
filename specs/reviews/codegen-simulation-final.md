# Final review: input validation + code generation + simulation

Branch `feature/codegen-simulation`, range `5a48f91..6f7ffc2` (12 commits), base `geomech`
at 5a48f91. Repo venv: Python 3.12.13, numpy 2.5.4, scipy 1.18.1, casadi 3.8.2. Read-only on
`geomech/`, `tests/`, `examples/`; nothing committed.

## Verdict: GO WITH CHANGES

Each feature holds on its own, and the three fit together on the normal path. The user
story below runs end to end from the docs. The reduced `A`, `B` that `lqr` builds from the
generated linearization agree with a finite difference of `accelerations`, the S² projected
solve: 6.5e-10 abs on `|A| ≈ 12` and 1.3e-12 abs on `|B| ≈ 50`, for the quadrotor with a
load at hover. Every gate passes.

The joins between the three features have one silent failure (R1) and two traps in
`collect_params`/`StateLayout` (R2, R3), and a student following the docs hits both. Fix
R1 to R3 before merging. R4 to R6 are documentation and bookkeeping.

## Required changes

**R1 (silent wrong controller): `generate_module` of a `LinearizedSystem` ignores the
layout's input order and variable order.**

- What goes wrong:
  - `geomech/codegen/module.py:202-220` (`_linearized`) lays out the `B_lin` columns in
    `lin.inputs` order and the rows/states in `lin.rows`/`lin.states` order.
    `_linearization_data` (`:223-235`) records that order.
  - The standard-form path reorders its rows by the layout (`_eom_rows`, `:122-128`), and
    `input_matrix` puts its columns in the layout's input order.
  - So `StateLayout(V, inputs=[M, f])` together with `linearize(eom, V, inputs=[f, M])` is
    accepted. `B_lin` column 0 is then `δf`, while `u[0:3]` is `M`.
- Reproducer (quadrotor with a load, hover, 10° off):
  - `lqr(..., module=plant)` accepts, because the equilibrium check runs on the
    standard-form module, which is consistent with the layout.
  - It reports closed-loop eigenvalues with max real part −0.60.
  - The closed loop diverges: `|x_L| = 35 m` and `|δx| ≈ 7e3` after 3 s.
- A different variable order (`SystemVariables(vectors=[q, x_L])` for `linearize` against a
  layout built from `[x_L, q]`) is also accepted. `lqr` then fails with the misleading
  `LinAlgError: lqr: the reduced M_lin is singular at x_eq (condition number inf)`.
- Fix: in `_linearized`, require

  ```python
  [str(s) for s in lin.inputs] == [str(u.delta()) for u in layout.inputs]
  [str(r) for r in lin.rows]   == ["kin:" + str(v) for v in layout.config_variables]
                                  + [variation keys in layout order]
  ```

  and raise `ValueError` otherwise, for example:
  `generate_module: the linearization's inputs ['\delta{f}', '\delta{M}'] are not in the
  layout's input order ['M', 'f']; build the StateLayout with the inputs and
  SystemVariables passed to linearize`. Reordering the columns instead would also work.
- Defence in depth: `_reduction` (`geomech/simulation/control.py:273-298`) maps
  `states[k]` to `model.variables[k % nv]` by position only. It should compare
  `LINEARIZATION["states"]`/`["inputs"]` symbols with the model's variables and inputs.
- Tests: swapped inputs and swapped vector variables both raise at `generate_module`.

**R2: a variation vector is accepted as a parameter.**

- `collect_params([L, dW], V, inputs)`, the natural call after reading `docs/codegen.md:65`,
  returns `['J', '\eta_{R}', '\xi_{q}', 'e_3', 'g', 'l', 'm_L', 'm_Q']`. The ξ and η come
  from `δW`.
- `StateLayout._build_params` (`geomech/codegen/layout.py:246-262`) accepts them. Every
  generated function then requires `p["\eta_{R}"]`.
- The first symptom appears only inside `lqr`:
  `KeyError: "missing parameter '\\eta_{R}' (expected [...])"`.
- Fix:
  - `collect_params` (`layout.py:411-428`) leaves out the tangent and variation vectors of
    the layout's variables.
  - `_build_params` raises `ValueError` for them, with the wording of `check_inputs`
    ("is a tangent or variation vector of a configuration variable").

**R3: `collect_params` silently returns `[]` for non-`Expr` entries.**

- `collect_params(sf.values(), V, inputs)` gives `[]`, because `StandardFormEquation` has
  no `nodes` (`_symbols`, `layout.py:398-409`).
- Fix: raise `TypeError` for an entry that is not an `Expr`, or accept the three system
  types directly (`compute_eom` dict, `to_standard_form` dict, `LinearizedSystem`).
- `docs/codegen.md:65` should say what to pass, for example `collect_params([L, dW], V,
  inputs)` once R2 is in.

**R4 (docs): say that both modules share one `StateLayout`.**

- In `docs/codegen.md` (end of "The generated module") and `docs/simulation.md` (LQR
  section, before "The equilibrium is checked…"), add: "Generate the standard-form and the
  linearization module from the same `StateLayout`, built with the same `SystemVariables`
  and the same `inputs` list (same order) as `to_standard_form` and `linearize`." Add the
  R1 error to `docs/troubleshooting.md`.
- `docs/codegen.md:224`: "plus `q qᵀ` on its own `ω̇` column" → "plus a multiple `s q qᵀ`
  on its own `ω̇` column (`s` = the tangent inertia, see Simulation)". The code uses `s q qᵀ`
  since the F1 fix.

**R5 (CI): the LQR example is not linted or run in CI.**

- `.github/workflows/tests.yml:24,27` lint `geomech/ tests/ scripts/ examples/generated/`
  only, and no test runs `examples/inverted_rigid_pendulum_lqr.py`.
- Add the file to both ruff lines, and add a test step `python
  examples/inverted_rigid_pendulum_lqr.py` (1.7 s, asserts its own result). If that is
  deferred, record it in `specs/follow-ups.md`.

**R6 (bookkeeping): record the deferred items in `specs/follow-ups.md`.** See "Open items".

Nits (optional):

- `geomech/simulation/model.py:225`: "Examples: a zero mass, …" sits under `ValueError`. It
  belongs to `LinAlgError`.
- `docs/vbl.md:186`: "It does not synthesise controllers, discretise, simulate or generate
  code" is still true of `linearize`, but should point to [Code generation](codegen.md) and
  [Simulation](simulation.md).
- CHANGELOG `lqr` entry: add `eq_tol`.

## User story (public API and docs only)

The script lives in the reviewer's scratch directory. The equivalent steps:

1. **Model.** `x_L: Vector`, `q: S2`, `R: SO3`; `x_Q = x_L − l q`;
   `L = ½ m_Q|ẋ_Q|² + ½ m_L|ẋ_L|² + ½ Ω·JΩ − m_Q g e₃·x_Q − m_L g e₃·x_L`;
   `δW = Dot(x_Q.delta(), R*e3*f) + Dot(η, M)`; `V = SystemVariables(vectors=[x_L, q],
   matrices=[R])`. From `docs/examples.md`. `.delta()` of a composite expression worked
   without being shown anywhere.
2. **Derivation.** `compute_eom`, `to_standard_form(eom, V, inputs=[f, M])`,
   `linearize(eom, V, inputs=[f, M])`: OK.
3. **Parameters.** `collect_params([L, dW], V, [f, M])` gives ξ and η as parameters
   (**R2**). I had to read `collect_params` to find what to pass. `sf.values()` silently
   gives `[]` (**R3**). `collect_params([L], …)` works.
4. **Modules.** `StateLayout(V, inputs=[f, M], params=…)`, `generate_module(sf, layout,
   path, lagrangian=L)`, `generate_module(lin, layout, path)`, `load_module`: OK. The docs
   never say that both modules need the same layout and input order (**R4**), and nothing
   checks it (**R1**).
5. **State.** Packed by hand as `[x_L, q, R row-major, ẋ_L, ω, Ω]` from the layout table.
   This is documented, but there is no helper that packs named slot values; `state_from_env`
   takes only `Evaluator` envs (friction, follow-up).
6. **LQR at hover.** `q_eq = −e₃`, `R = I`, `u_eq = [(m_Q+m_L) g, 0, 0, 0]`;
   `lqr(linmod, x_eq, u_eq, p, Q=np.eye(16), R=np.eye(4), module=plant)`.
   - The reduced size 16 (3+2+3+3+2+3) has to be worked out before `ctrl.labels` exists. A
     scalar `Q` avoids that.
   - Result: K is 4×16, closed-loop max Re λ = −0.60.
7. **Closed loop from 10° off.** `simulate(plant, x0, 10.0, p, u=ctrl)`: angle to `q_eq`
   2e-4°, `|x_L|` 8.4e-5, `|δx|` 9.5e-5 at 10 s.
8. **Unforced energy.**
   - Relative drift over 10 s: 9.8e-11 (rkmk4, dt = 1e-3) and 1.4e-12 (solve_ivp, rtol =
     1e-10).
   - `constraint_drift` keys: `|q|-1`, `q.omega_q`, `R^T R-I`.
   - The whole script runs in 10 s under `-W error`.

Source was read only for: what `collect_params` expects (step 3). Guesses: `Scalar("f")`
for a single scalar (`getScalars` is only shown for several names), and whether one layout
serves both modules (step 4).

## Consistency (item 2)

| Join | Result |
|---|---|
| `StateLayout` order vs `LINEARIZATION` states/rows vs `lqr` error order | Consistent when both come from the same `SystemVariables`: `[δx_L, ξ_q, η_R, δẋ_L, δω_q, δΩ_R]` against layout `[x_L, q, R]`. Not checked when they differ (**R1**) |
| Input order `check_inputs` → `StateLayout` → `G`/`B` columns → `u` | `input_matrix` follows the layout; `B_lin` follows `lin.inputs`. They agree only if the user passes the same list (**R1**) |
| S² tangent basis in `lqr` vs S² projected solve in `accelerations` | Agree: reduced `A`/`B` vs central differences of `accelerations` + `kinematics` in the coordinates `q = exp((Bξ)^) q_eq`, `R = exp(η̂)`, `ω = Bν`. Differences 6.5e-10 and 1.3e-12 abs |
| `simulate` state vs `lqr` error | Same layout slices (`Model`); restarting from `traj.state[-1]` works (checked on the spherical pendulum) |

## Gates (item 3)

| Gate | Result |
|---|---|
| `python -m pytest tests/ -q -W error` | 3271 passed, 2 xfailed (the two older `test_scalar_adversarial.py` ones), 61.3 s; slowest test 5.3 s (double rigid pendulum energy, rkmk4) |
| `ruff format --check` / `ruff check` on `geomech/ tests/ scripts/ examples/generated/ examples/inverted_rigid_pendulum_lqr.py` | 149 files formatted; all checks passed |
| `mypy geomech --ignore-missing-imports` | no issues in 46 files |
| `mkdocs build --strict` | exit 0 (only the upstream MkDocs-2 notice) |
| `scripts/regenerate_generated.py --check` | 18 modules up to date; also exits 0 with scipy and casadi blocked (numpy-only, as in the lint job) |
| `examples/inverted_rigid_pendulum_lqr.py` under `-W error` | 1.7 s, `|η|` = 4.846e-08 rad at 5 s, matches `docs/simulation.md` |
| CI dry read | Lint job: ruff, mypy, numpy; the regeneration script puts the repo root on `sys.path` itself, so it works without an install. Test job: `.[dev]` includes scipy (simulation tests), 3.12 adds `.[casadi]`. Docs job: `.[docs]`. The gap is R5. No 3.11+ syntax found (`tomllib`, `Self`, `except*`, …); 3.10/3.11 were not run locally |

`import geomech, geomech.codegen, geomech.simulation` in a fresh interpreter from `/tmp`
loads neither scipy nor casadi (item 4).

## Docs and CHANGELOG (item 5)

- The Unreleased section covers all three features:
  - inputs under Changed (two bullets, plus F1–F5) and the size fix under Fixed;
  - codegen under Added, with F1–F6 and the CSE advisory under Fixed;
  - simulation under Added, with F1 scaling, F2 `module=`/warning, F3 restart, F4–F6
    input errors, `output_every`.
- Every behaviour change in the three fix logs is listed. The one omission is `eq_tol`
  (nit).
- Snippets run:
  - `docs/codegen.md` worked example and "Used directly": they run, and the generated
    `mass`/`bias` are byte-identical to `examples/generated/spherical_pendulum_s2.py`.
  - `docs/simulation.md` lead snippet, on the spherical pendulum: `accelerations` (q·ω̇ =
    0), `simulate`, `traj["q"]`, and the restart from `traj.state[-1]`.
  - The example output block matches.
- Consistency between pages:
  - `docs/conventions.md` (inputs) and `docs/troubleshooting.md` (seven new rows) match
    `check_inputs`.
  - `docs/codegen.md:224` is out of date (R4).
  - No troubleshooting rows exist for codegen/simulation errors. Acceptable, but the R1
    error should get one.

## Open items (item 6)

Deferred by the reviews and implementation reports and **not** in `specs/follow-ups.md`:

1. Codegen `energy` imports private `linearize` helpers (`_perturbation_variables`,
   `_partition_terms`). Make them public.
2. `energy_function` does not support prescribed parameters (`compute_eom(...,
   params=[l])`).
3. `TimeDerivative(Transpose(x))` is rejected as "TimeDerivative of a non-leaf"
   (conforming; no engine output holds it).
4. Simulation F7 remainder: stepper overhead of 35–55 µs per stage (preallocation, fused
   `mass`/`bias`/`input_matrix` call); `Trajectory.energy` is not vectorised.
5. Anisotropic `Q` on an S² block depends on `tangent_basis`, which jumps 90° under tiny
   changes of `q_eq`. Neither documented nor recorded; accept a user basis or document it.
6. The LQR example is not in CI (R5).
7. Suite time grew from about 17 s to 61 s (rkmk4 10 s rows). Possible `slow` marker.
8. A helper to pack a state from named slot values (user story step 5).
9. R1 to R3 here, if not fixed before merge.

Recorded or closed:

- Inputs F1–F5 are fixed.
- `set` for `compute_eom` params stays accepted on purpose (documented in
  `docs/conventions.md`).
- The `lqr` ω_eq = 0 assumption is now enforced by `module=`.

## PII scan (item 7)

`git log 5a48f91..HEAD -p` (18 403 lines) and the 12 commit messages:

- **Author/committer:** every commit is `vkotaru <17078723+vkotaru@users.noreply.github.com>`
  (the noreply address). No other email address appears.
- **None found:** home paths, `/tmp` or scratch paths, host names, IPs, tokens, keys, or
  employer/affiliation.
- **Hits:**
  1. `specs/codegen.md:8`: named two sibling projects (a chat front-end and the paper
     repo) as consumers. Suggested text: "that simulation, control design and other tools
     can call".
  2. `specs/reviews/simulation-implementation-report.md:185`: referred to the
     contributors' local notes file outside the repo by path.
     Suggested text: "which CI does not lint (examples are scripts, exempt from the
     unused-import rule)".
  3. `specs/follow-ups.md:28`: "The accompanying paper's S² variation paragraph". This was
     already in the base (`5a48f91`, line 30); the range only renumbers it. Generalise it
     if the paper should stay unnamed.
  4. All 12 commit messages end with a `Claude-Session: https://claude.ai/code/session_…`
     trailer. This is the configured attribution, left in on purpose; it is a session
     link, not PII. The user decides.

## Fix log

Fixed on `feature/codegen-simulation` on top of 6f7ffc2; nothing committed. Line numbers
are those of the working tree after the fix.

**R1 reproduced first** (a script from the public API, the quadrotor with a load, hover,
`Q = I₁₆`, `R = I₄`, 10° off):

- Before: inputs `[M, f]` for `linearize` against `StateLayout(V, [f, M])`: `lqr` accepted
  (max Re λ = −0.596) and the closed loop reached `|x_L| = 51 m` after 3 s.
  `SystemVariables(vectors=[q, x_L])` for `linearize`: `LinAlgError: lqr: the reduced M_lin
  is singular at x_eq (condition number 2.38e+17)`. The dynamic rows stay in the
  `compute_eom` order while the kinematic rows and states follow the other order.
- After: both raise at `generate_module`, before any module exists:
  - `ValueError: generate_module: the linearization's inputs ['\delta{M}', '\delta{f}'] are
    not the variations of the layout's inputs ['f', 'M'] in that order (expected
    ['\delta{f}', '\delta{M}']); build the StateLayout with the SystemVariables and the
    inputs (same order) passed to linearize and to_standard_form, and generate both modules
    from that one StateLayout`
  - `ValueError: generate_module: the linearization's configuration variables ['q', 'x_L',
    'R'] (rows ['kin:q', 'kin:x_L', 'kin:R', '\delta{x_L}', '\xi_{q}', '\eta_{R}']) are not
    the layout's variables ['x_L', 'q', 'R'] in that order (rows [...]); build the
    StateLayout …`

| Item | Change | Where |
|---|---|---|
| R1 | `_check_linearization_order`: `lin.rows` must be `kin:<var>` per layout variable then the variation keys in layout order, `lin.states` the `perturbation_variables` of the layout, `lin.inputs` `[u.delta() for u in layout.inputs]`; `ValueError` naming both orders | `geomech/codegen/module.py:213-238`, called at `:242`; docstring `Raises` |
| R1 | `_reduction` compares `LINEARIZATION` kinematic row labels with the `LAYOUT` variables and the input symbols with `\delta{u}` of the `LAYOUT` inputs (`ValueError` with both orders) | `geomech/simulation/control.py:309-331` |
| R1 | `lqr(module=...)`: `_check_same_layout` compares `LAYOUT` sizes, configuration, velocity, acceleration and input entries (parameters may differ) before anything is evaluated; the message names both variable and input orders: `lqr: module does not have the layout of lin_module: lin_module has the variables ['x'] and inputs ['f', 'T'], module has the variables ['x'] and inputs ['T', 'f']; build one StateLayout …` | `geomech/simulation/control.py:208`, `:253-280` |
| R2 | `collect_params` leaves out `TS2`/`TSO3` leaves (`ξ`, `η`; `ω`, `Ω` were already slots); `StateLayout._build_params` raises `ValueError: StateLayout: parameter '\xi_{q}' is a tangent or variation vector of a configuration variable; …`. Only the leaves: a `TimeDerivative` of one (a prescribed rate, as the codegen fuzz test generates) is unchanged | `geomech/codegen/layout.py:269-275`, `:545` |
| R2 | The code generator names a stray tangent/variation vector as such instead of advising `params` (which now rejects it) | `geomech/codegen/emitter.py:267-274` |
| R3 | `collect_params` walks an `Expr`, `StandardFormEquation`, `LinearizedSystem`, a dict (its values) or any iterable of these, and raises `TypeError: collect_params: cannot find symbols in str 'dW'; pass the expressions of the system ([L, dW], or [L] plus the equations), the dict of compute_eom or of to_standard_form, a LinearizedSystem, or an iterable of these` | `geomech/codegen/layout.py:495-547` |
| Friction (step 5) | `StateLayout.pack(values=None, /, **slots)` (slot identifiers `x_L`, `q`, `R`, `x_L_dot`, `omega_q`, `Omega_R`, or printed symbols as dict keys; configuration required, velocities default to zero; `R` 3×3 or 9) and `unpack(state)` (also `(N, n_state)`) | `geomech/codegen/layout.py:371-430` |
| Open item 1 | `perturbation_variables`/`partition_terms` public in `geomech/dynamics/vbl.py:470`, `:556` (old names aliased at `:653-654`); `geomech/codegen/energy.py:34` imports the public names | |
| R4 | "One layout for everything" (`docs/codegen.md:118`), `collect_params` and `pack`/`unpack` (`:65-87`), `s q qᵀ` (`:253`); "One layout for both modules" with a runnable quadrotor snippet using `pack` (`docs/simulation.md:210`), `pack` in "The state" (`:46-52`), anisotropic `Q` on S² (`:274`); `docs/troubleshooting.md:258` "Code generation and simulation", five rows, each with a trigger in `TRIGGERS` (`tests/test_phase4_adversarial.py:948-961`) | |
| R5 | `.github/workflows/tests.yml:24,27` lint the example; the example now uses `collect_params([L, dW], …)` and `layout.pack` (`examples/inverted_rigid_pendulum_lqr.py:56-59`, `:73`, `:86`); the test runs it in a subprocess under `-W error` (< 10 s, 1.7 s here) and checks `|eta|` at 5 s < 1e-6 | |
| R6 | `specs/follow-ups.md` items 24-29: energy and prescribed parameters, `TimeDerivative(Transpose(x))`, stepper overhead and `Trajectory.energy`, anisotropic `Q` and `tangent_basis` (now documented), suite time and a `slow` marker, a `pack` for inputs | `specs/follow-ups.md:55` |
| Nits | `model.py:221-224` "Examples …" moved under `LinAlgError`; `docs/vbl.md:186-189` points to Code generation and Simulation; CHANGELOG `lqr` entry has `eq_tol`; CHANGELOG Fixed entry for R1-R6 and Added entry for `pack`/`unpack` | |
| PII | `specs/codegen.md:7-8` → "that simulation, control design and other tools can call"; `specs/reviews/simulation-implementation-report.md:185` → "the example scripts are exempt from unused-import lint by project convention"; the PII section of this file no longer quotes the two hits verbatim | |

New tests, `tests/test_codegen_simulation_integration.py` (23 tests, 4.5 s): both R1
reproducers raise at `generate_module` with both orders in the message; `lqr` rejects a
plant with another input order and a tampered `LINEARIZATION` (inputs swapped, kinematic
rows relabelled); the user story from the public API (`collect_params([L, dW])` gives
exactly `J, e3, g, l, m_L, m_Q`; `lqr` K 4×16, max Re λ = −0.596; closed loop from 10°
off for 4 s: error and angle below 20 %, `|x_L| < 0.1`, drift < 1e-10; unforced energy
over 1 s to 1e-9 and the drift keys); R2/R3 (`collect_params` of `eom`, `eom.values()`,
`sf`, `sf.values()`, `lin`, `L`, a mixed list; `TypeError` for `3`, `"L"`, `None`,
`object()`, `[1.0]`, also inside a list; `ξ`/`η` rejected as params; a stray `ξ` named by
the generator); `pack`/`unpack` (slices, printed keys, R as 9 entries, defaults, `(N, n)`,
the four errors, scalar slots); the example run.

Grep of the whole tree (outside `venv/`, `.git/`, `site/`) for references to the local
development setup, home or temporary paths: one hit,
`specs/reviews/phase3-review.md:372` ("No stored output contains `/home/`, `/tmp/` or a
user name"), a statement of absence, left. `specs/reviews/README.md` has none. Not
changed: `specs/follow-ups.md:28` ("The accompanying paper's S² variation paragraph",
hit 3 above, from the base) and the `Claude-Session` trailers (hit 4).

Gates after the fix:

| Gate | Result |
|---|---|
| `python -m pytest tests/ -q -W error` | 3299 passed, 2 xfailed, 63.9 s (was 3271 + 2 xfailed; +23 integration tests, +5 troubleshooting rows) |
| `ruff format --check` / `ruff check` on `geomech/ tests/ scripts/ examples/generated/ examples/inverted_rigid_pendulum_lqr.py` | 150 files formatted; all checks passed |
| `mypy geomech --ignore-missing-imports` | no issues in 46 files |
| `mkdocs build --strict` | exit 0 (only the upstream MkDocs-2 notice); the new anchors `#code-generation-and-simulation`, `#the-state-layout`, `#the-generated-module` exist |
| `scripts/regenerate_generated.py --check` | 18 modules up to date (no generated module changed) |
| `python -W error examples/inverted_rigid_pendulum_lqr.py` | 1.7 s, `|η|` = 4.846e-08 rad at 5 s, unchanged |
