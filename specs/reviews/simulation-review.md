# Simulation (`geomech.simulation`): adversarial review

Range `cc1f58f..8be087f` on `feature/codegen-simulation` ("Add simulation and validation").
Contract: `specs/simulation.md`. Tests: `tests/test_simulation_adversarial.py` (20 pass,
8 strict xfail, 3.6 s). Full suite with the new file: 3121 passed, 21 xfailed in 76 s.
Ruff is clean. Nothing is committed.

Every reference is plain NumPy with no geomech call:

- S² accelerations: a constrained Lagrangian in R³ with one multiplier per sphere,
  `M q̈ = F + Cᵀλ` with `q·q̈ = −|q̇|²`, mapped by `ω̇ = q × q̈`. The generalized force of
  `δW = ξ·f` is `f × q`.
- Rigid bodies: Euler's equations with `Ṙ = R Ω̂` in R¹².
- LQR gains: central differences of those dynamics in the coordinates `R = R_eq exp(η̂)`,
  `q = exp((Bξ)^) q_eq`, `ω = Bν`, then `scipy.linalg.solve_continuous_are`.

## Verdict: GO WITH CHANGES

The dynamics are right. The S² solve, the RKMK4 stepper and the LQR gains agree with the
independent references to round-off or to the finite-difference error. No finding produces
wrong equations of motion at physical parameter values.

The findings concern robustness, as ranked below:

- F1: the S² stabiliser does not scale with the inertia.
- F2: `lqr` does not check that `x_eq` is an equilibrium.
- F3: a trajectory's final state cannot be passed back as `x0`.
- F4–F6: input-handling gaps.
- F7: performance.

Fix F1–F4 before merging. They are small changes. F5–F7 can follow.

## What held (passing tests)

| Check | Result |
|---|---|
| **S² solve vs multiplier solve.** Double spherical pendulum (both S² rows couple to the other ω̇), 10 random states and inputs | rel. err ≤ 1e-14 |
| **Same, quadrotor with load.** The S² row couples to ẍ_L, and its mass block is rank 2 | ≤ 1e-14 |
| **Same, gyroscopic spherical pendulum.** `L = ½ m l²‖ω‖² + b e₃·ω − m g l e₃·q`; the reference has the force `2b e₃ × q̇` | ≤ 2e-15 |
| **Mass ratio 1e6.** Double spherical pendulum with m₂/m₁ = 1e-6 | no false alarm, err 1.8e-11 |
| **RKMK4 order 4, SO(3) with a time-varying torque.** Free body, Ω₀ = (2, −3, 1.5) | ratios 15.86, 15.94 |
| **Same, S² with a gyroscopic term and a time-varying torque** | 16.79, 16.41 |
| **Same, S² × S² double spherical pendulum, time-varying inputs** | 14.46, 15.28 (dt 1e-2 → 2.5e-3) |
| **Same, R³ × S² × SO(3) quadrotor with load, time-varying thrust and torque** | 16.33, 16.18 |
| **Sign of the right action.** The odd `dexp⁻¹` term is flipped for SO(3) in a mutation test | the order drops to 2 (ratio 4.00). The shipped sign `+½ r×g` is correct |
| **Sign on S²** | the order is 4 with either sign and about 3 with no bracket. At leading order `[V, G] ∝ ω × ω̇ ∥ q`, which is the isotropy of q, so the sign cannot be observed. The shipped `−½` is the correct left-action `dexp⁻¹` |
| **S² update.** With g = 0, q(t) = exp(tω̂) q₀ and ω is constant | exact to 1e-13 (left rotation about ω, consistent with q̇ = ω × q) |
| **SO(3) update.** With J = I and no torque, R(t) = R₀ exp(tΩ̂) | exact to 1e-13 (R·exp(hΩ), consistent with Ṙ = R Ω̂) |
| **LQR, inverted rigid pendulum.** A, B, K vs FD + CARE | ≤ 1.2e-12 |
| **LQR, spherical pendulum at q = +e₃ and −e₃** | K ≤ 3e-12; the basis is (e₁, ±e₂), with no degeneracy at the poles |
| **LQR, quadrotor with load at hover.** 16 coordinates, coupled S² row | A 6.5e-12, K 1.6e-9 |
| **Basis switch with isotropic Q.** `tangent_basis` jumps 90° between two equilibria 1e-9 apart; Q is isotropic per S² block | u(state) agrees to 1e-6 |
| **solve_ivp drift.** `constraint_drift` is computed on the raw states | the outputs have R projected to 1e-13 |
| **Input sizes.** Wrong-length constant or callable `u`; `u` given to a system with no inputs | `ValueError` with the expected message |
| **Feedback at stages.** `u(t, state)` is evaluated at RKMK4 stage states | `u = −kω` gives monotone energy decay |
| **x0 tolerance.** `|q|` off by 5e-10 / 2e-9 | accepted / `ValueError` naming `q` |

## Findings, ranked

### F1 — S² solve: the stabiliser `q qᵀ` has unit scale (accuracy loss and false alarms) — *dynamics, robustness*

`geomech/simulation/model.py:133` adds `np.outer(q, q)` with coefficient 1 to a block whose
tangent part scales like `m l²`. The condition check `MAX_COND = 1e12` (`model.py:24`,
`:142`) runs on that unscaled matrix. The result depends on the units, not only on the
physics.

- **Silent accuracy loss.** Take the double spherical pendulum with l₁ = l₂ = 1 and the
  state in `test_s2_solve_light_link_accurate`. A reference with each sphere's rows
  divided by its own inertia has condition number 2.6. Relative error of `accelerations`
  against it:

  | m₂/m₁ | Error | `accelerations` |
  |---|---|---|
  | 1e-6 | 1.8e-11 | accepted |
  | 1e-9 | 3.1e-8 | accepted (silent) |
  | 1e-12 | 4.5e-6 | raises |
  | 1e-14 | 8.4e-4 | raises |

  The unscaled KKT solve stays at 2e-16 throughout, so the problem is well posed.
- **False alarm.** A 1 t mass on a 100 km tether (`m = 1e3`, `l = 1e5`, SI units;
  tethered-satellite scale) gives cond 9.97e12. `simulate` raises `LinAlgError ... is
  singular`, but the solve is exact: ω̇ = 8.25e-5 = (g/l) sin θ. A pendulum with 1 g and
  10 µm also raises.
- **Fix, checked in scratch.** Add `s·q qᵀ` with `s = max|A[rows, acc]|`, then
  row-equilibrate A and b (divide each row by its max) before the solve and before
  `cond`. With that change:
  - the light-link error is 1e-16 at every ratio down to 1e-14;
  - the scaled cond is 2.9;
  - the tether's cond is 1.1;
  - `adv_r3_s2_pendulum_with_gyro` still reports 5.9e16, and the quadrotor with load at
    m_Q/m_L = 1e12 (ill-posed in the (x_L, q) coordinates) still reports 4.2e12. So
    genuine singularities are still caught.
- **Tests:** `test_s2_solve_light_link_accurate` and `test_s2_solve_tether_not_rejected`
  are xfail.

### F2 — `lqr` accepts a non-equilibrium `x_eq` — *wrong controller, silent*

`geomech/simulation/control.py:170` checks only that `x_eq` lies on the manifold, and
`:172` sets `a_eq = 0`. The report lists this under "not done". Its effect is worse than
the report suggests.

- **Reproducer.** Rigid pendulum about the pivot, `x_eq = R_up exp(30° ê₁)` at rest,
  `u_eq = 0`, Q = diag(100 I, 10 I), R = I.
  - At that point `accelerations` gives Ω̇ = (2.46, −0.12, 0.008), so it is not an
    equilibrium.
  - `lqr` returns a gain whose `closed_loop_eigenvalues` all have real part ≤ −1.37.
  - The closed loop from `x_eq` settles at ‖δx‖ = 0.22 rad away from `x_eq`, with a
    constant `u ≈ (−3.3, 0, 0)`.
  - Nothing warns. The reported eigenvalues suggest convergence to `x_eq`.
- **Fix.** The linearization module has no residual function, so `lqr` cannot check by
  itself. Two options:
  - an optional `module=` argument (the standard-form module) and a check that
    `‖accelerations(module, x_eq, u_eq, p) − a_eq‖` is small, else `ValueError`;
  - emit a residual function in the linearization module.

  At least, state in the docstring and docs that `x_eq` is assumed to be an equilibrium.
- **Test:** `test_lqr_rejects_a_non_equilibrium` is xfail.

### F3 — a trajectory's final state is not a valid `x0`; `project_state` ignores ω — *silent drift, then crash*

Two causes:

- `project_state` (`geomech/simulation/integrators.py:125-138`) normalises q and takes the
  polar factor of R. It leaves ω with its normal component, so `traj.state` from
  `solve_ivp` is not on TS². This contradicts the docstring ("projected onto the
  manifold").
- RKMK4 translates ω (`integrators.py:271`) while it rotates q. The drift in q·ω therefore
  grows with the truncation error.

Reproducers:

| Case | Drift in q·ω | Restart from `traj.state[-1]` |
|---|---|---|
| spherical pendulum, `solve_ivp`, rtol = 1e-6, 5 s | 1.2e-8 after projection | `ValueError: velocity slot ... must be orthogonal` |
| spherical pendulum, `rkmk4`, dt = 0.02, 20 s | 5.4e-9 | the same `ValueError` |
| quadrotor with load, `solve_ivp`, rtol = 1e-6 | 4.2e-6 | the same `ValueError` |

Restarting from the last state is the natural way to chain segments, for example to
switch controllers.

Fix:

- `project_state` should also set `ω ← ω − (q·ω) q`;
- RKMK4 can apply the same O(error) projection to ω after each step, which does not
  change the order;
- alternatively, `check_state` could use a relative tolerance on q·ω.

Tests: `test_solve_ivp_final_state_is_a_valid_x0` and `test_rkmk4_final_state_is_a_valid_x0`
are xfail.

### F4 — a callable `u` that returns `None` gives zero input silently — *wrong control, silent*

`_input_function` (`integrators.py:153`) passes the return value to `Model.as_input`.
`as_input` maps `None` to zeros (`model.py:82`), because `None` is also the "no inputs"
argument. A controller missing its `return` statement therefore runs open loop: the
trajectory equals the zero-input one exactly (max difference 0.0). Fix: in `ufun`, reject
a `None` return with a `TypeError`. Test: `test_callable_input_returning_none_is_rejected`
is xfail.

### F5 — a NaN input is blamed on the mass matrix — *cosmetic*

`u = [20, 0, 0, nan]` raises `LinAlgError: the mass system gives non-finite accelerations
(condition number 22.2)`. `as_input` should check that the inputs are finite and name `u`.
For a callable, the check can run once at `t0`, since `simulate` already evaluates
`ufun(t0, x0)`. Test: `test_nan_input_is_reported_as_an_input_error` is xfail.

### F6 — `t_span = inf` raises `OverflowError` — *cosmetic*

`int(np.ceil(inf))` at `integrators.py:205`. `_span` should require a finite `tf`. Test:
`test_infinite_t_span_is_a_value_error` is xfail.

### F7 — performance: the cost is in the generated `np.cross`, then in stepper overhead

Profile of a 2 s double-rigid-pendulum run (rkmk4, dt = 1e-3), and the per-stage cost of
each system:

| System | µs/stage | generated `mass`+`bias`+`G` | rest of `solve_accelerations` | stepper overhead |
|---|---|---|---|---|
| double rigid pendulum | 344 | 286 (`bias` alone 286) | ~5 | small |
| quadrotor with load | 172 | 100 | 18 | 55 |
| free rigid body | 54 | 14 | 6 | 34 |

Findings:

- **The generated code dominates.** The generated `bias` of the double pendulum calls
  `np.cross` 26 times on 3-vectors (it is already CSE'd). Each call costs about 10 µs:
  `np.cross` goes through `moveaxis` and `normalize_axis_tuple`. In the profile, `cross`
  and its helpers take 4.9 s of 6.2 s.
- **One-line fix.** `NumpyTarget.cross` (`geomech/codegen/numpy_target.py:115-116`) could
  emit a module-level `_cross(a, b)` helper (three multiply-subtracts, like
  `integrators._cross`). Patching the module's `np.cross` this way gives:
  - `bias` 286 → 72 µs;
  - 1 s of simulation 1.44 → 0.57 s (2.5×);
  - an identical trajectory (max diff 0.0).
- **Stepper overhead** is 35–55 µs per stage, from:
  - tuple and list construction in `_lin` and `dexpinv`;
  - `np.concatenate` and `np.array([...])` in `generator`;
  - the per-call `_check_params`/`_arr`/`_mat` preludes of three generated functions.

  Preallocated work arrays and a single fused `mass`/`bias`/`input_matrix` call with
  bound parameters would bring it near 10 µs.
- **`Trajectory.energy`** calls the generated `energy` in a Python loop over 10 001
  states. It takes about 1.7 s of the 15.9 s double-pendulum test, also mostly
  `np.cross`.

The `_cross` helper alone gives about 2.5× on the slow rows. With the stepper changes, the
total is about 4–5× (double pendulum about 344 → 80 µs per stage), with no change to
tolerances, dt or horizons. Using `solve_ivp` instead is already 10–20× faster at the
spec's tolerances, but the table requires `rkmk4` at dt = 1e-3.

## Notes (no test)

- **`traj.energy()` for `solve_ivp` uses the projected states** (q normalised, R polar).
  The `Trajectory` docstring says the states are projected. `docs/simulation.md` does not
  say which states the energy table uses. At rtol = 1e-10 the difference is about 1e-10,
  below every tolerance. At looser rtol the projected energy hides part of the drift.
  State this in the docs.
- **The callable `u` is evaluated at rejected `solve_ivp` trial states** and at
  non-monotone times. This is harmless for `LQRController`. A stateful controller (an
  integral term) would be corrupted. Document it.
- **A 2-D constant `u`** (e.g. shape (2, 2) for 4 inputs) is silently flattened.
- **Anisotropic Q on an S² block depends on the basis.** `tangent_basis` switches its
  reference axis discontinuously in `q_eq`, so the meaning of `xi_q[0]` jumps by 90° under
  a 1e-9 change of `q_eq`. Isotropic Q per S² block is invariant (tested). Consider
  accepting a user basis, or document the caveat.
- **`simulate` stores every step**, so 1000 s at dt = 1e-3 for the quadrotor with load is
  about 190 MB. An output stride option would help long runs.
- **The S² reduction in `lqr` assumes `ω_eq = 0`** (`q·δω = −(ξ×q)·ω_eq` otherwise). This
  is consistent with F2: only equilibria are supported.

## Fix log

On `feature/codegen-simulation` after 2db6623; nothing committed. Seven of the eight strict
xfails flipped with their markers removed and their assertions unchanged. The F2 xfail
stays in place, for the reason given under F2.

**Gates.**

- Full suite, `-W error`: 3270 passed, 3 xfailed in 59.8 s (57 s at 2db6623). The xfails
  are the F2 test here and two unrelated `test_scalar_adversarial.py` ones.
- Every energy/momentum row and both cross-checks pass at their tolerances.
- `ruff format --check` and `ruff check` are clean on geomech/, tests/, scripts/ and the
  example.
- `mypy geomech --ignore-missing-imports` (as in CI): clean. Plain `mypy` reports only the
  missing scipy stubs, as at 2db6623.
- `mkdocs build --strict`: OK. `scripts/regenerate_generated.py --check`: 18 modules up
  to date.
- `examples/inverted_rigid_pendulum_lqr.py` under `-W error`: output unchanged
  (|η| = 4.846e-08 rad at 5 s).

**F1 (stabiliser scale).**

- `geomech/simulation/model.py:158`: the stabiliser is `s q qᵀ` with
  `s = ‖P_q D‖_F/√2`, where D is the row's own diagonal block. Both derived forms
  (`−m l² I` and `Hat(q)Hat(q) m l²`) project to `−m l² P_q`, and `‖P_q‖_F = √2`, so
  `s = m l²` exactly. For an anisotropic tangent inertia, s is the RMS of the two tangent
  singular values.
- `:161`: every row of A and b is then divided by its max |entry| before `solve` and
  `cond`. A zero row is left as is, so it stays singular.
- Double spherical pendulum, relative error against the scaled reference:

  | m₂/m₁ | Error |
  |---|---|
  | 1e-6 | 6e-17 |
  | 1e-9 | 2.7e-16 |
  | 1e-12 | 2.7e-16 |
  | 1e-14 | 2.2e-16 |

  The condition number is 2.87 at 1e-14.
- Tether (1 t on 100 km): relative error 1.4e-20, condition number 1.1. A 1 g mass on
  10 µm is accepted.
- Genuine singularities still raise:
  - `adv_r3_s2_pendulum_with_gyro`: condition 4.0e16–1.1e17;
  - quadrotor with load at m_Q/m_L = 1e12: 4.2e12;
  - the zero-mass and singular-J tests still pass.
- Documented in the `accelerations` docstring and in docs/simulation.md under "Scaling".

**F2 (equilibrium).**

- `geomech/simulation/control.py:140`: new keywords `module=None` and `eq_tol=1e-8`.
- `_check_equilibrium` (`:235`):
  - checks that the layout matches `lin_module`;
  - requires the velocities to be zero (to `tol`) (`:250`);
  - requires the residual r (`mass·a_eq + bias + G u_eq`, S² rows projected) to satisfy
    `‖r‖∞ ≤ eq_tol · max(‖bias‖∞, ‖G u_eq‖∞, ‖A_lin T‖∞)`. The `A_lin T` term is needed
    because the bias is exactly zero at the upright pendulum;
  - otherwise raises `ValueError` with the residual and `accelerations(module, …)`
    (`:266`).
- Without `module=`, `lqr` emits a `UserWarning` (`:221`) after the argument checks.
- Measured:
  - 30° pendulum: residual 2.45, relative 0.58, a = (2.465, −0.124, 0.008);
  - quadrotor with load at 1.01× hover thrust: relative 9.9e-3;
  - upright pendulum, hover and q = ±e₃: accepted.
- **Stopped, not flipped.** `test_lqr_rejects_a_non_equilibrium` calls `lqr` without
  `module=` and expects `ValueError`. The requested design warns in that case: the
  linearization module carries no residual, so `lqr` cannot check by itself. Flipping
  the test needs one of two decisions: add `module=` to that call (the new test
  `tests/test_simulation_lqr.py:165` asserts exactly that, plus the warning, a moving
  `x_eq` and a layout mismatch), or make the missing check an error. The test is left
  as a strict xfail; under `-W error` it fails with the UserWarning.
- Call sites updated to pass `module=`: three in tests/test_simulation_lqr.py, three in
  tests/test_simulation_adversarial.py and the example. No assertion changed.
- `test_lqr_feedback_independent_of_tangent_basis_for_isotropic_weights` uses a
  non-equilibrium on purpose. It got `@pytest.mark.filterwarnings` for that warning
  (tests/test_simulation_adversarial.py:549).

**F3 (restart).**

- `project_state` (`integrators.py:137`) now also sets `ω ← ω − (q·ω) q`.
- RKMK4 does the same after each step (`:405`).
- Measured:
  - RKMK4, spherical pendulum, 20 s at dt = 0.02: q·ω max 2.8e-17 (was 5.4e-9);
  - solve_ivp at rtol = 1e-6, spherical pendulum: raw drift 2.5e-8, projected 3.5e-18;
  - quadrotor with load: raw 2.5e-4, projected 8.9e-16. The restart works.
- The order ratios stay at about 16:

  | Case | Ratios |
  |---|---|
  | SO(3) | 15.86, 15.94 |
  | S² gyroscopic | 15.67, 15.84 |
  | S²×S² | 13.50, 14.86 |
  | quadrotor with load | 16.31, 16.17 |

- New test `tests/test_simulation.py:338` covers both methods.

**F4 (input).**

- `integrators.py:163`: a callable returning `None` raises `TypeError` naming t.
- `Model.as_input` (`model.py:92`) rejects `ndim > 1` with the shape, and still accepts
  a scalar for a single input.

**F5 and F6.**

- F5: `model.py:101` reports non-finite inputs by index. For a callable, the time is
  appended (`integrators.py:170`).
- F6: `_span` requires a finite `t0, tf` (`integrators.py:261`).
- New test `tests/test_simulation.py:312`.

**Notes.**

- docs/simulation.md now states that:
  - `traj.state`, `energy()` and the slots of `solve_ivp` are projected;
  - a callable `u` is evaluated at RKMK4 stages and at rejected or non-monotone
    `solve_ivp` trial steps, so stateful controllers must not mutate state in `u`.
- New `simulate(..., output_every=k)` (`integrators.py:231–235`) stores every k-th step
  plus the last, for both methods. Test at `tests/test_simulation.py:355`.
- `Trajectory.energy` is still a loop over the generated `energy`, now via `np.fromiter`.
  Measured: 0.233 s for 10 001 states of the double rigid pendulum, against 4.6 s for the
  run. The generated function takes one state at a time, so it is not vectorised.
- Energy table in the docs re-measured; the changes are at round-off:

  | System | Method | Before | After |
  |---|---|---|---|
  | free rigid body | rkmk4 | 4.6e-15 | 3.0e-15 |
  | rigid pendulum about pivot | solve_ivp | 4.9e-10 | 4.8e-10 |
  | quadrotor with load | rkmk4 | 1.8e-11 | 2.6e-11 |
  | quadrotor with load | solve_ivp | 6.5e-11 | 5.2e-11 |

**F7 (re-measured).** Double rigid pendulum, 10 s, rkmk4, dt = 1e-3 (40 000 stages), best
of 3, same machine and session:

- 2db6623: 4.44 s (111 µs/stage);
- now: 4.64 s (+4%).

The cost is the row equilibration, about 3 µs/stage, plus the per-step ω projection and
storage bookkeeping. The reviewer's pre-F7 figure was 344 µs/stage (about 13.8 s).

### Arbitration (orchestrator)

`test_lqr_rejects_a_non_equilibrium` called `lqr` without `module=`; the equilibrium
check needs the standard-form module (the linearization module has no residual), and
without it `lqr` warns rather than raises so that callers without a standard-form module
keep working. The test now passes `module=`; its assertion (ValueError mentioning the
equilibrium) is unchanged, and its xfail marker is removed.
