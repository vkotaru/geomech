# Spec: simulation and validation (`geomech.simulation`)

Roadmap feature 6. Branch `feature/codegen-simulation`, after codegen.

## Goal

Integrate the equations produced by the pipeline forward in time on the manifold, and use
the result as an independent physical check: conserved quantities must be conserved, and a
controller designed from `linearize` must stabilise the nonlinear system. Simulation
consumes only `geomech.codegen` output (no tree evaluation in the inner loop).

Optional dependency: `scipy` (extra `geomech[sim]`). Importing `geomech` must not require it.

## Accelerations from standard form

The generated module gives `mass(state, p)`, `bias(state, p)`, `input_matrix(state, p)`
with `mass @ a + bias + input_matrix @ u = 0` (signs as in `docs/conventions.md`).

- R, R³, SO(3) rows: the blocks are square and invertible for a physical system.
- S² rows: the equation is meaningful only on the tangent plane (`docs/conventions.md`).
  The derived mass block is not always rank 2: it is `Hat(q)Hat(q) m l² = −m l² (I − qqᵀ)`
  (rank 2) for the quadrotor with a suspended load, but `−m l² I` (full rank) for the
  spherical pendulum, where the engine already used q·ω̇ = 0. Either way the unprojected
  3×3 solve is not the equation of motion. The
  acceleration ω̇ also satisfies the tangency constraint q·ω̇ = 0, from
  d/dt(q·ω) = (ω×q)·ω + q·ω̇ = 0. So each S² row block is replaced by

      (P_q · row) + q qᵀ on its own ω̇ column,    with right-hand side −P_q (bias + G u),

  where P_q = I − q qᵀ. The result is square and full rank. The other columns of an S²
  row (couplings to other variables' accelerations) are projected by P_q as well.

`accelerations(state, u, p)` builds this system and solves it with `np.linalg.solve`. A
singular system raises `np.linalg.LinAlgError` with the condition number and the state.
Test: on every S² system, the solved ω̇ is tangent (|q·ω̇| < 1e-12) and satisfies the
tangent-plane equation to 1e-10.

## Integrators

State: the flat layout of `StateLayout` (configuration slots, then velocity slots).

1. `method="rkmk4"` (default): fixed-step Runge–Kutta–Munthe-Kaas of order 4.
   - SO(3): R ← R · expm_so3(h Ω̄), with the RK stages combined in the Lie algebra.
   - S²: q ← expm_so3(h ω̄) · q (rotation of q about ω), which preserves |q| = 1 exactly.
   - R, R³ and all velocity slots: classical RK4.
   Constraint drift (|q| − 1, ‖RᵀR − I‖) stays at round-off without projection.
2. `method="solve_ivp"`: `scipy.integrate.solve_ivp` (RK45/DOP853) on the embedded state,
   with `rtol`/`atol` passed through, and projection of the outputs onto the manifold
   (q/|q|; R by polar decomposition); the drift before projection is reported.

API:

```python
from geomech.simulation import simulate, Trajectory
traj = simulate(module, x0, t_span, p, u=None, method="rkmk4", dt=1e-3)
traj.t, traj.state, traj.constraint_drift, traj.energy(L_fn)   # energy if given
```

`u` is `None`, a constant array, or a callable `u(t, state)`. `x0` is validated: |q| = 1,
ω ⊥ q, R orthogonal with det +1 (tolerance 1e-9), else `ValueError` naming the slot.

## Energy and momentum

`geomech.codegen` gains `energy_function(L, layout)`: the total energy
E = Σ p·v − L evaluated from the Lagrangian, generated through the same path. Simulation
tests assert:

| System | Conserved quantity | Tolerance over 10 s, rkmk4, dt = 1e-3 |
|---|---|---|
| spherical pendulum | E | 1e-8 relative |
| free rigid body | E, inertial angular momentum R J Ω | 1e-9 |
| rigid pendulum about pivot | E | 1e-8 |
| double rigid pendulum | E | 1e-7 |
| elastic spherical pendulum | E | 1e-7 |
| quadrotor with suspended load, zero inputs | E | 1e-7 |

and the same with `method="solve_ivp"`, `rtol=1e-10`.

Independent cross-checks against closed forms written in NumPy in the test (not through
geomech): spherical pendulum vs Lee et al. eq. 5.14 integrated in R³; free rigid body vs
Euler's equations J Ω̇ = JΩ × Ω with Ṙ = R Ω̂. Trajectories agree to 1e-8 over 2 s.

## Control demo (VBL → LQR)

`geomech.simulation.lqr(lin_module, x_eq, u_eq, p, Q, R)`: evaluates the generated A, B,
M blocks at an equilibrium, forms Ā = M⁻¹A, B̄ = M⁻¹B, solves the continuous algebraic
Riccati equation (`scipy.linalg.solve_continuous_are`) and returns K. State error on
SO(3) is η = vee(log(R_eqᵀ R)), with `log_so3` added to `geomech.numeric`. On S² the
error is projected to the tangent plane, and the gain is computed in a 2-D tangent basis
(the 3-D S² perturbation is redundant); this is documented and tested.

Demo and test: a rigid pendulum about its pivot (SO(3)), inverted equilibrium, starting
20° off with a random axis, is stabilised by u = −K δx; ‖η‖ < 1e-3 after 5 s. A second
test: an inverted spherical pendulum (q near +e₃) with a tangent torque input converges
to the upright point.

## Definition of done

- All energy/momentum rows and both cross-checks pass for both integrators.
- The two LQR tests pass; the demo is `examples/inverted_rigid_pendulum_lqr.py` and a
  docs page `docs/simulation.md` (layout, integrators, S² solve, energy checks, LQR).
- `simulate` on every system in `tests/vbl_systems.py` runs 1 s without error.
- No scipy import at `import geomech`; ruff, mypy clean; CHANGELOG entry.

## Non-goals

Contact, constraints beyond the manifold, variable-step Lie-group integrators, MPC,
real-time performance.
