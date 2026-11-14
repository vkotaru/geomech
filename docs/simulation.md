# Simulation

`geomech.simulation` integrates the modules written by
[`geomech.codegen`](codegen.md) forward in time on the manifold, and designs LQR gains
from a generated linearization. It calls only the generated functions (`mass`, `bias`,
`input_matrix`, `kinematics`, `energy`; `M_lin`, `A_lin`, `B_lin`), never the expression
trees. The integrators are also an independent check of the derived equations:
conserved quantities stay conserved, and the trajectories agree with closed forms written
directly in NumPy.

```python
from geomech.simulation import accelerations, simulate, lqr

a = accelerations(module, state, u, p)              # mass @ a + bias + G @ u = 0, S² on T_qS²
traj = simulate(module, x0, (0.0, 10.0), p, u=None, method="rkmk4", dt=1e-3)
traj.t, traj.state, traj.constraint_drift, traj.energy()
ctrl = lqr(lin_module, x_eq, u_eq, p, Q, R)         # LQRController: K, error(state), u(t, state)
traj = simulate(module, x0, 5.0, p, u=ctrl)
```

scipy is an optional dependency (`pip install "geomech[sim]"`), imported only by
`method="solve_ivp"` and by `lqr`; `import geomech.simulation` does not import it, and
`method="rkmk4"` needs only NumPy.

## The state

The state is the flat array of the module's `LAYOUT` (see
[the state layout](codegen.md#the-state-layout)): every configuration slot in variable
order, then every velocity slot.

| Variable | Configuration | Velocity | Constraint |
|---|---|---|---|
| `s: Scalar` | `s` (1) | `ṡ` (1) | |
| `x: Vector` | `x` (3) | `ẋ` (3) | |
| `q: S2` | `q` (3) | `ω` (3) | `‖q‖ = 1`, `q · ω = 0` (`q̇ = ω × q`) |
| `R: SO3` | `R` (9, row-major) | `Ω` (3, body frame) | `RᵀR = I`, `det R = 1` (`Ṙ = R Ω̂`) |

`simulate` checks `x0` against these constraints (tolerance 1e-9) and raises
`ValueError` naming the slot, for example
`x0: the S2 slot 'q' must have unit norm, got |q| = 1.01` or
`x0: the velocity slot '\omega_{q}' must be orthogonal to 'q'`. `check_state(module,
state, tol)` runs the same check. A trajectory gives the values of a slot by its printed
symbol or its generated name: `traj["q"]` is `(N, 3)`, `traj["R"]` is `(N, 3, 3)`,
`traj["omega_q"]` is `(N, 3)`.

## Accelerations and the S² solve

The generated equations read `mass @ a + bias + input_matrix @ u = 0` with
`M = −(physical inertia)` ([Conventions](conventions.md#signs)). For R, R³ and SO(3)
rows the blocks are square and invertible for a physical system.

An S² equation is meaningful only on the tangent plane of its point
([Conventions](conventions.md#s2-equations-live-on-the-tangent-plane)), and the derived
mass block may be singular or not:

| System | S² mass block | Rank |
|---|---|---|
| quadrotor with a suspended load | `Hat(q) Hat(q) m_Q l² = −m_Q l² (I − qqᵀ)` | 2 |
| spherical pendulum | `−m l² I` (the engine already used `q · ω̇ = 0`) | 3 |

In both cases the unprojected 3×3 solve is not the equation of motion: it either fails or
gives `q · ω̇ ≠ 0`, and the trajectory leaves the sphere. The missing equation is the
tangency of the acceleration: differentiating `q · ω = 0` gives
`(ω × q) · ω + q · ω̇ = 0`, and the first term vanishes, so `q · ω̇ = 0`. `accelerations`
replaces each S² row block (all its columns, including couplings to other variables'
accelerations) by

    P_q · row  +  q qᵀ on its own ω̇ column,        right-hand side  −P_q (bias + G u),

with `P_q = I − q qᵀ`, and solves the square system with `np.linalg.solve`. Multiplying
the new rows by `qᵀ` gives `q · ω̇ = 0`; multiplying by `P_q` gives the tangent-plane
equation. The test suite checks both (to 1e-12 and 1e-10) on every S² system of the
oracle gate, including the two of the table.

A singular system (condition number above 1e12: a zero mass, a singular inertia, a
kinetic energy that ignores a velocity such as `½ m |ẋ + l q̇|²`) raises
`np.linalg.LinAlgError` with the condition number and the state. `accelerations` checks
the condition number on every call; `simulate` checks it at `x0` and afterwards only
catches a failed or non-finite solve.

## Integrators

### `method="rkmk4"` (default): Runge–Kutta–Munthe-Kaas of order 4

A fixed-step Lie-group integrator. The state is moved by a group element
`exp(U)` of the product group that acts on it:

- each S² point by a rotation from the left, `q ← expm(u) q` (so `q̇ = ω × q` has generator `ω`);
- each SO(3) point by a rotation from the right, `R ← R expm(u)` (generator `Ω`);
- every other slot (R, R³ configuration and all velocities) by a translation.

One step of size `h` applies classical RK4 to the ODE of `U` in the Lie algebra,
`U̇ = dexp⁻¹_U(A(exp(U)·yₙ))`, `U(tₙ) = 0` (Munthe-Kaas, *High order Runge–Kutta methods
on manifolds*, Appl. Numer. Math. 29, 1999; Hairer, Lubich, Wanner, *Geometric Numerical
Integration*, §IV.8):

    K₁ = h F(0),   K₂ = h F(K₁/2),   K₃ = h F(K₂/2),   K₄ = h F(K₃),
    U = (K₁ + 2K₂ + 2K₃ + K₄) / 6,   yₙ₊₁ = exp(U) · yₙ,
    F(V) = dexp⁻¹_V(A(exp(V) · yₙ))    at t, t + h/2, t + h/2, t + h.

`dexp⁻¹_V(G) = G − ½[V, G] + 1/12 [V, [V, G]] + …` is truncated after the `ad²` term.
On so(3) ≅ R³ the bracket is the cross product; for the right action on SO(3) the sign of
the odd term flips (`dexp⁻¹_{−V}`); on the vector slots every bracket vanishes and the step
is exactly classical RK4. Because `V = O(h)` and the `ad³` coefficient is zero, the
truncation changes each `Kᵢ` by `O(h⁵)` and the method keeps order 4. Without the
brackets the method drops to order 2 on a tumbling rigid body (error ratio 4.0 per
halving instead of 16).

`q` and `R` are only ever multiplied by rotation matrices, so `|q| − 1` and `‖RᵀR − I‖`
stay at round-off (≈ 5e-14 over 10 s) without projection. `q · ω = 0` is a velocity
constraint that the dynamics preserve (`d/dt (q·ω) = 0` from the S² solve); the
integrator keeps it to its truncation error (≤ 3e-11 over 10 s at `dt = 1e-3`).

`dt` is shortened, if needed, so that it divides `tf − t0`; the trajectory holds every
step.

### `method="solve_ivp"`

`scipy.integrate.solve_ivp` (`ivp_method="DOP853"` by default, `rtol`/`atol` passed
through) on the embedded state with `ẏ = (kinematics(y), accelerations(y))`. The solver
does not know the manifold, so the outputs (at spacing `dt`) are projected onto it: `q/|q|`
and the polar factor of `R` (`U Vᵀ` of its SVD). `traj.constraint_drift` reports the
drift *before* projection: `|q|−1` reaches about 3e-10 over 10 s at `rtol = 1e-10`.

### Inputs

`u` is `None` (zero inputs), a constant array of `LAYOUT["n_inputs"]` entries, or a
callable `u(t, state)` evaluated at every stage (an `LQRController` is one).

## Energy and momentum checks

`traj.energy()` evaluates the module's `energy(state, p)` (generated with
`generate_module(..., lagrangian=L)`, see [Total energy](codegen.md#total-energy)) at every
output state; `traj.energy(fn)` uses another function `fn(state, p)`, such as
`energy_function(L, layout)`.

`tests/test_simulation_conservation.py` runs the systems of `tests/vbl_systems.py` for
10 s with zero inputs from the states of `tests/simulation_systems.py`. Relative energy
drift, `max |E(t) − E(0)| / |E(0)|`:

| System | Tolerance | rkmk4, dt = 1e-3 | solve_ivp, rtol = 1e-10 |
|---|---|---|---|
| spherical pendulum | 1e-8 | 4.0e-13 | 6.5e-10 |
| free rigid body | 1e-9 | 4.6e-15 | 3.7e-11 |
| rigid pendulum about pivot | 1e-8 | 2.8e-14 | 4.9e-10 |
| double rigid pendulum | 1e-7 | 1.9e-14 | 9.4e-11 |
| elastic spherical pendulum | 1e-7 | 3.4e-12 | 8.5e-10 |
| quadrotor with suspended load | 1e-7 | 1.8e-11 | 6.5e-11 |

The inertial angular momentum `R J Ω` of the free rigid body is conserved to 1e-9 by both
methods. Two closed forms, written in NumPy without geomech and integrated with DOP853 at
`rtol = 1e-13`, agree with both methods to 1e-8 over 2 s: the spherical pendulum (Lee et
al. Eq. 5.14, `m l² ω̇ = m g l e₃ × q` for the potential `m g l e₃·q`, `q̇ = ω × q`, in
R³ × R³) and the free rigid body (Euler's equations `J Ω̇ = J Ω × Ω`, `Ṙ = R Ω̂`).

## LQR from the linearization

`lqr(lin_module, x_eq, u_eq, p, Q, R, a_eq=None)` evaluates the generated `M_lin`,
`A_lin`, `B_lin` ([Variation-Based Linearization](vbl.md)) at an equilibrium, forms
`Ā = M⁻¹A`, `B̄ = M⁻¹B` in reduced coordinates, solves the continuous algebraic Riccati
equation (`scipy.linalg.solve_continuous_are`) and returns an `LQRController` with the
gain `K`, the Riccati solution `P`, the reduced `A`, `B`, the closed-loop eigenvalues and
the state error. Calling the controller, `ctrl(t, state) = u_eq − K δx(state)`, is the
input signature of `simulate`.

The state error uses the perturbation coordinates of `linearize`:

| Variable | Configuration error | Velocity error |
|---|---|---|
| R, R³ | `x − x_eq` | `ẋ − ẋ_eq` |
| SO(3) | `η = vee(log(R_eqᵀ R))` (`geomech.numeric.log_so3`) | `Ω − Ω_eq` |
| S² | `Bᵀ (q_eq × q)` (2 entries) | `Bᵀ (ω − ω_eq)` (2 entries) |

The SO(3) error follows from `δR = R η̂`: `R = R_eq exp(η̂)`. On S², `δq = ξ × q` with
`ξ ⊥ q`, so `ξ = q_eq × q` to first order.

**The S² tangent basis.** The linearization of an S² point has 3-D perturbations `ξ` and
`δω`, but at an equilibrium (where `ω = 0`) both are orthogonal to `q`: the normal
direction is a zero mode that no input reaches, so the Riccati equation has no
stabilising solution in 3-D. `lqr` therefore works in the 2-D basis
`B = tangent_basis(q_eq)` (orthonormal, `b₁ × b₂ = q_eq`): with `T` block diagonal (`B`
for each S² perturbation, `I` elsewhere), it uses `Tᵀ M T`, `Tᵀ A T`, `Tᵀ B`, so the S²
rows are projected onto the tangent plane like the S² solve above. A singular S² block of
`M_lin` (`Hat(q) Hat(q) m l²` for the quadrotor with a load, where the full `M_lin` has
rank 17 of 18) becomes invertible. `Q` is given in the reduced coordinates
(`ctrl.labels` names them: `xi_q[0]`, `xi_q[1]`, `eta_R[0]`, …), as a scalar, a diagonal
or a matrix.

### Example: the inverted rigid pendulum

`examples/inverted_rigid_pendulum_lqr.py` derives a rigid pendulum about its pivot on
SO(3) (Lee et al. Eq. 6.8), generates its standard form and linearization, designs an LQR
gain about the inverted equilibrium `R_eq = diag(1, −1, −1)` and simulates the closed
loop from 20° off about a random axis:

```
closed loop   : [-158.082  -11.554  -11.554   -3.213   -3.213   -3.163]

open loop: |eta| goes from 20.0 deg to 137.0 deg in 1 s
closed loop:
  t = 0.0 s   |eta| = 3.491e-01 rad
  t = 1.0 s   |eta| = 1.545e-02 rad
  t = 2.0 s   |eta| = 6.498e-04 rad
  t = 5.0 s   |eta| = 4.846e-08 rad
max |RᵀR − I| over the run: 1.9e-14
```

The tests also stabilise an inverted spherical pendulum (`q` 15° from `+e₃`, the maximum
of `m g l e₃·q`, with its tangent input) to `|q − e₃| < 1e-4` in 5 s, and the quadrotor
with a suspended load about hover.

## Not covered

Contact, constraints beyond the manifold, variable-step Lie-group integrators, MPC and
real-time performance. Each RKMK4 stage costs one call of `mass`, `bias` and
`input_matrix` plus a small dense solve: 10 s at `dt = 1e-3` (40 000 stages) takes about
2.5 s for the single-body systems and 15 s for the double rigid pendulum, whose generated
blocks are the slowest.
