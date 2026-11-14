# Simulation (`specs/simulation.md`): implementation report

Branch `feature/codegen-simulation`, on top of cc1f58f. Nothing is committed.

## What was built

`geomech/simulation/` has four files. Simulation calls only the functions of the generated modules (`LAYOUT`, `mass`, `bias`, `input_matrix`, `kinematics`, `energy`; `M_lin`, `A_lin`, `B_lin`, `LINEARIZATION`). It never evaluates a tree.

| File | Contents |
|---|---|
| `model.py` | `Model`: the module's `LAYOUT`, parsed once per module and cached in a `WeakKeyDictionary`. Also `accelerations`, `solve_accelerations` (the inner-loop form), `check_state`, `STATE_TOL = 1e-9`, `MAX_COND = 1e12` |
| `integrators.py` | `simulate`, `Trajectory`, the RKMK4 stepper (`_LieRK`, `_rkmk4`), the `solve_ivp` path, `constraint_drift`, `project_state` |
| `control.py` | `lqr`, `LQRController`, `tangent_basis` |
| `__init__.py` | the public names |

Other changes:

- **`geomech.numeric.log_so3`**: the principal logarithm, returned as a `(3,)` vector. The angle comes from `atan2(|skew|, (tr R − 1)/2)`. Near 0 it uses a series for θ/sin θ. Near π it takes the axis from the symmetric part, with the sign taken from the skew part.
- **`tests/codegen_systems.tangent_accelerations` was promoted, not duplicated.** The helper is deleted, and `tests/test_codegen_energy.py` now calls `geomech.simulation.accelerations`.
- **Spec correction about rank 2.** The texts that claimed the S² mass block is always rank 2 are fixed: `specs/simulation.md`, `specs/codegen.md` (marked as a correction) and the `generate_module` docstring. `docs/codegen.md` already stated the rank correctly; its last sentence now links to `docs/simulation.md`.
- **Example:** `examples/inverted_rigid_pendulum_lqr.py`. It derives the pendulum from its Lagrangian, generates both modules in memory, runs `lqr` and a closed-loop simulation, and asserts the result.
- **Docs:** `docs/simulation.md` with an mkdocs nav entry ("Simulation"). It covers the layout, the S² solve with the rank table, both integrators with the RKMK4 derivation, energy checks with the achieved numbers, LQR with the S² reduction, and the example.
- **README:** pipeline step 6, the `[sim]` install line and the package tree.
- **CHANGELOG:** an entry under Unreleased/Added.
- **New test helper:** `tests/simulation_systems.py` holds physical parameters and on-manifold initial states for the nine VBL systems.

API:

```python
from geomech.simulation import accelerations, check_state, simulate, Trajectory, lqr
a = accelerations(module, state, u=None, p)              # flat a, LAYOUT order
traj = simulate(module, x0, t_span, p, u=None, method="rkmk4", dt=1e-3,
                rtol=1e-10, atol=1e-12, ivp_method="DOP853", tol=1e-9)
traj.t, traj.state, traj.constraint_drift, traj.energy(fn=None), traj["q"], traj.max_drift()
ctrl = lqr(lin_module, x_eq, u_eq, p, Q, R, a_eq=None)   # LQRController
ctrl.K, ctrl.P, ctrl.A, ctrl.B, ctrl.T, ctrl.labels, ctrl.closed_loop_eigenvalues
ctrl.error(state); ctrl(t, state)                         # u_eq − K δx
```

Decisions the spec left open:

- **`accelerations` takes the module first:** `accelerations(module, state, u, p)`.
- **Condition-number check.** A singular system is detected by its condition number exceeding 1e12. `np.linalg.solve` alone does not raise on a numerically singular matrix: the degenerate `adv_r3_s2_pendulum_with_gyro`, with `L = ½ m |ẋ + l q̇|²`, solved without error at cond ≈ 6e16. The SVD costs about 10 µs, so:
  - the public `accelerations` checks it on every call;
  - `simulate` checks it once at `x0`;
  - the inner loop raises only on a failed or non-finite solve.
- **`t_span`** is `(t0, tf)` or `tf`. `dt` is shortened so that it divides the span. For `solve_ivp`, `dt` sets the output spacing (`t_eval`).
- **`constraint_drift`** is a dict of `(N,)` arrays keyed `"|q|-1"`, `"q.omega_q"` and `"R^T R-I"` (Frobenius norm), using the generated slot names. The velocity constraint `q·ω` is reported as well, because it is not preserved exactly (see below).
- **Q and R for `lqr`** are given in the reduced coordinates, as a scalar, a diagonal or a matrix. `ctrl.labels` names each coordinate (`xi_q[0]`, `eta_R[2]`, `delta_x_L_dot[1]`).
- **scipy is imported only** inside `_solve_ivp` and `lqr`, with an `ImportError` that gives the `geomech[sim]` hint. A subprocess test checks that `import geomech, geomech.simulation` loads neither scipy nor casadi.

## Integrator derivation (RKMK4)

The state is acted on by the product group:

- each S² point is rotated from the left, `q ← expm(u) q`, with generator ω (`q̇ = ω × q`);
- each SO(3) point is rotated from the right, `R ← R expm(u)`, with generator Ω (`Ṙ = R Ω̂`);
- every other slot (the R and R³ configurations and all velocities) is translated.

Write `y = exp(U)·yₙ`. Then `U` solves `U̇ = dexp⁻¹_U(A(exp(U)·yₙ))` with `U(tₙ) = 0`. One step applies classical RK4 to that ODE:

```
K₁ = h F(0),  K₂ = h F(K₁/2),  K₃ = h F(K₂/2),  K₄ = h F(K₃)
U  = (K₁ + 2K₂ + 2K₃ + K₄)/6,      yₙ₊₁ = exp(U)·yₙ
F(V) = dexp⁻¹_V(A(exp(V)·yₙ))  at  t, t+h/2, t+h/2, t+h
```

The `dexp⁻¹` series is truncated after the ad² term:

- left action: `dexp⁻¹_V(G) ≈ G − ½ V×G + (1/12) V×(V×G)`;
- right action on SO(3): `dexp⁻¹_{−V}`, the same with the sign of the odd term flipped;
- the vector part commutes, so there the step is exactly classical RK4.

The order argument: V = O(h) and B₃ = 0, so the truncation changes F by O(h⁴) and each Kᵢ by O(h⁵). The local error stays O(h⁵), which gives order 4. The docstring cites Munthe-Kaas, Appl. Numer. Math. 29 (1999), and Hairer, Lubich & Wanner, *GNI* §IV.8.

The truncated-series form was chosen over the "minimal commutator" variant because its order argument is a two-line proof. It costs only a few cross products per stage.

Order verified numerically. The table gives the error at t = 2 s against a DOP853 reference at rtol = 1e-13, for dt = 0.04, 0.02, 0.01 and 0.005:

| System | Error ratios per halving |
|---|---|
| free rigid body, Ω₀ = (2, −3, 1.5) | 15.88, 15.96, 15.99 |
| spherical pendulum | 15.27, 15.66, 15.84 |
| rigid pendulum about the pivot | 17.02, 16.53, 16.28 |
| quadrotor with load | 16.00, 16.11, then 4.8 (the reference's accuracy floor, 3e-12) |
| free rigid body, brackets removed | 4.00, 4.00, 4.00 (order 2) |

The last row is evidence that the `dexp⁻¹` correction is needed on SO(3). On the spherical pendulum the brackets made no difference: the method stays order 4 with or without them. The test asserts a ratio between 13 and 19 on the free rigid body and on the spherical pendulum.

## S² solve

Implemented as written in the spec, on both kinds of block:

- the spherical pendulum's mass block is `−m l² I` (rank 3);
- the quadrotor with load's S² block is `Hat(q)Hat(q) m_Q l²` (rank 2), and its full `mass` is singular.

The test runs on every S² system of the oracle gate (35 builders of `ALL_SYSTEMS`), at 5 sample states with physical parameters and random inputs. It checks `|q·ω̇| < 1e-12·max(1, |ω̇|)`, the tangent-plane residual `‖P_q(M a + f + G u)‖ < 1e-10·scale`, and the residual of the other rows below 1e-10.

The exception is `adv_r3_s2_pendulum_with_gyro`. Its kinetic energy `½ m |ẋ + l q̇|²` is degenerate for every parameter value, so the test asserts that it raises `LinAlgError` with the condition number.

## Achieved numbers vs the spec table

All runs are 10 s with zero inputs, from `tests/simulation_systems.SETUPS`. The energy is the generated `energy`. Drift is `max |E − E₀| / |E₀|`, and every |E₀| is above 0.3.

| System | Spec tol | rkmk4, dt = 1e-3 | solve_ivp DOP853, rtol = 1e-10 | rkmk4 time | solve_ivp time (nfev) |
|---|---|---|---|---|---|
| spherical pendulum | 1e-8 | 4.0e-13 | 6.5e-10 | 2.6 s | 0.3 s (2438) |
| free rigid body, E | 1e-9 | 4.6e-15 | 3.7e-11 | 2.3 s | 0.1 s (1061) |
| free rigid body, R J Ω (abs, \|π\| = 3.1) | 1e-9 | 1.1e-13 | 9.7e-11 | | |
| rigid pendulum about the pivot | 1e-8 | 2.8e-14 | 4.8e-10 | 2.8 s | 0.2 s (1727) |
| double rigid pendulum | 1e-7 | 1.9e-14 | 9.4e-11 | 14.3 s | 0.8 s (1856) |
| elastic spherical pendulum | 1e-7 | 3.4e-12 | 8.5e-10 | 2.7 s | 0.2 s (3278) |
| quadrotor with load, u = 0 | 1e-7 | 1.8e-11 | 6.5e-11 | 7.1 s | 1.8 s (12125) |

Every tolerance is met with at least 15× margin for solve_ivp and at least 5000× for rkmk4, so none was loosened.

With zero inputs the quadrotor with load falls freely. A test checks that its centre of mass follows the parabola under −g e₃ to 1e-6.

### Constraint drift

Max over 10 s:

| Method | \|q\| − 1 | ‖RᵀR − I‖ | \|q·ω\| |
|---|---|---|---|
| rkmk4 | ≤ 4.8e-14 | ≤ 3.5e-14 | ≤ 2.7e-11 |
| solve_ivp, before projection | up to 2.7e-10 | up to 1.1e-9 | up to 1.5e-10 |

For rkmk4, the q·ω value of 2.7e-11 is for the quadrotor with load. q·ω is a velocity constraint: the S² solve preserves it exactly in continuous time, and the integrator keeps it to its truncation error. No velocity projection was added; the spec asks for round-off only for |q| and RᵀR. The test asserts:

- below 1e-12 for |q| − 1 and ‖RᵀR − I‖;
- below 1e-9 for q·ω.

The solve_ivp outputs are projected onto the manifold.

### Closed-form cross-checks over 2 s

The references are written in NumPy and integrated with DOP853 at rtol = 1e-13 and atol = 1e-14. Max deviation over all outputs and state entries:

| Reference | rkmk4 | solve_ivp, rtol = 1e-12 |
|---|---|---|
| Spherical pendulum: Lee et al. Eq. 5.14 in R³ × R³ (`m l² ω̇ = m g l e₃ × q` for this Lagrangian's potential, `q̇ = ω × q`) | 8.8e-12 | 6.9e-12 |
| Free rigid body: Euler's equations `J Ω̇ = J Ω × Ω`, `Ṙ = R Ω̂`, in R¹² | 4.0e-13 | 3.2e-12 |

The gate is 1e-8. For the cross-checks the geomech solve_ivp path uses rtol = 1e-12. At the table's rtol = 1e-10 the error would be about 1e-9, which still passes.

## LQR results

- **Rigid pendulum about the pivot, inverted** (`R_eq = diag(1, −1, −1)`, ρ = −0.5 e₃, J as in the tests, Q = diag(100·I₃, 10·I₃), R = I):
  - open-loop eigenvalues ±2.23, ±1.57 and 0, 0 (the rotation about ρ feels no gravity);
  - closed loop: −3.26, −2.16, −1.39 (pairs).
  - From 20° about random axes (seeds 0, 1, 2), ‖η(5 s)‖ = 1.5e-4, 4.9e-5 and 8.6e-5. The gate is 1e-3.
  - Without control the tilt of Rρ from e₃ more than doubles in 1 s.
- **The example** (J = diag(0.3, 0.3, 0.02)):
  - ‖η‖ = 0.349 → 1.5e-2 (1 s) → 6.5e-4 (2 s) → 4.8e-8 (5 s);
  - ‖RᵀR − I‖ ≤ 1.9e-14;
  - runtime 2.2 s.
- **Inverted spherical pendulum** (q₀ 15° from +e₃, ω₀ = 0, Q = diag(10, 10, 1, 1), R = I):
  - K is 3×4, and its row for the input along q is zero;
  - closed-loop eigenvalues −3.25 and −3.17 (pairs);
  - after 5 s, |q − e₃| = 4.9e-7 and |ω| = 1.5e-6 (gates 1e-4 and 1e-3);
  - uncontrolled, |q − e₃| > 0.2 after 2 s.
- **S² reduction tests:**
  - The 3-D `Ā, B̄` of the spherical pendulum fail the PBH rank test at λ = 0 (uncontrollable normal mode); the reduced 4-D pair passes.
  - T is orthonormal and tangent.
  - `tangent_basis` is orthonormal, with b₁ × b₂ = q, at four points including ±e₃.
- **Extra test, quadrotor with load about hover:**
  - the full `M_lin` has rank 17 of 18, from the `Hat(q)Hat(q)` block; the reduced one is invertible;
  - with Q = 100 and R = I, K is 4×16 and the slowest closed-loop mode is about −1/s;
  - after 5 s the error has dropped more than 10×.
- **State error:** η is exact for R = R_eq exp(η̂), and ξ = Bᵀ(q_eq × q).

## Tests and checks

New tests: 127.

| File | Tests | What it checks |
|---|---|---|
| `test_simulation.py` | 87 | The S² solve on every S² system; the rank-2 and rank-3 blocks; closed-form accelerations for a point mass and a rigid body; the singular error from both `accelerations` and `simulate`; size and module-type errors; eight x0 violations naming the slot (S² norm, ω ⊥ q, SO(3) orthogonality, det −1, second SO(3) of the double pendulum, quadrotor-with-load ω, NaN) plus size and tolerance; argument errors; constant and callable `u` against exact polynomials (both methods); `dt` shortening; `Trajectory` slot access and energy; 1 s on all nine VBL systems with both methods (18 tests); solve_ivp projection and drift; the scipy-missing hint; no scipy at import (subprocess); `log_so3` round trips (9 angles up to π − 1e-8) and at exactly π |
| `test_simulation_conservation.py` | 27 | The energy table (6 systems × 2 methods); R J Ω (2); free fall; rkmk4 drift (6); the two cross-checks × 2 methods; RKMK4 order (2) |
| `test_simulation_lqr.py` | 13 | The two spec LQR tests (pendulum with 3 seeds; spherical pendulum), the open-loop instability, the PBH reduction test, `tangent_basis` (4), the quadrotor-with-load reduction, the error coordinates, argument errors |

Commands run:

- Full suite `python -m pytest tests/ -q -W error`: **3008 passed, 2 xfailed in 77 s**. Before this feature it was 2881 passed in about 17 s.
- `ruff format --check` and `ruff check` on `geomech/ tests/ scripts/ examples/generated/ examples/inverted_rigid_pendulum_lqr.py`: clean. `ruff check examples/` still reports the existing star-import errors of the older example scripts (`sanity_check.py` and others), which CI does not lint and which `../the contributor notes exempts. No existing file was reformatted.
- `mypy geomech --ignore-missing-imports`: no issues (46 files).
- `mkdocs build --strict`: passes. griffe asked for annotations on the public signatures, so `ArrayLike` annotations were added.
- `python scripts/regenerate_generated.py --check`: 18 modules, no drift. The generated modules are unchanged; the only edit to `module.py` is a docstring.
- `python examples/inverted_rigid_pendulum_lqr.py`: runs in 2.2 s and passes its own assertion.

## Timing and cost

- **Suite time grows by about 60 s.** The new tests take about 60 s, most of it in the 10 s rkmk4 table runs:
  - double rigid pendulum 16.6 s;
  - quadrotor with load 7.8 s;
  - each of the others about 2.5–3 s.
- **Why rkmk4 is slow:** at dt = 1e-3 it makes 40 000 stage evaluations. Each costs one `mass` + `bias` + `input_matrix` + a solve, plus about 35 µs of Python overhead in the stepper. The double pendulum's generated blocks alone take about 320 µs per call (codegen report).
- **solve_ivp is 10–20× faster** at these tolerances.
- **Optimisations made:** a hand-written `_cross` and `_expm` in the stepper, since `np.cross` dominated the first profile.
- **Possible follow-up:** if the suite time matters, mark the double-pendulum and quadrotor-with-load rkmk4 rows as slow, or run them at T = 5 s.

## Not done or open

- **CI does not lint `examples/inverted_rigid_pendulum_lqr.py`.** The file is ruff-clean, but the lint job covers only `scripts/` and `examples/generated/`. It is not run in CI either.
- **No velocity projection for q·ω in rkmk4** (see the drift section).
- **The `x_eq` check in `lqr` does not test whether the point is an equilibrium**, only that it lies on the manifold. A non-equilibrium point gives the LQR of the linearization at that point.
- **The S² reduction assumes ω_eq = 0**, which holds at any equilibrium on S². The basis is `tangent_basis(q_eq)`, fixed at the equilibrium.
- **The `energy` helper of codegen still uses private `linearize` helpers.** This is unchanged from the codegen report.
- **Non-goals of the spec** were not attempted: contact, extra constraints, variable-step Lie-group integrators, MPC and real-time performance.
