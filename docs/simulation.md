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
traj.u                                              # inputs at the stored states, (N, n_inputs)
ctrl = lqr(lin_module, x_eq, u_eq, p, Q, R, module=module)   # checks (x_eq, u_eq)
traj = simulate(module, x0, 5.0, p, u=ctrl)
more = simulate(module, traj.state[-1], (5.0, 10.0), p)      # continue from the last state
```

`traj.u` holds the inputs at the stored states (zeros for `u=None`), recorded during
the run: with RKMK4 they are the values the first stage of each stored step used, so a
controller is never called twice at one state (with `solve_ivp`, whose stored states
are interpolated, it is called once more at each of them after the run). An exception
raised by a controller propagates as soon as it is raised, unchanged, with a note giving
the time. `geomech.viz.Report` reads `traj.u` to label a forced run's energy as "energy
change" rather than drift; a hand-built `Trajectory` without `u` gets the neutral
"energy change … inputs not recorded".

scipy is an optional dependency (`pip install "geomech[sim]"`), imported only by
`method="solve_ivp"` and by `lqr`; `import geomech.simulation` does not import it, and
`method="rkmk4"` needs only NumPy.

## The state

The state is the flat array of the module's `LAYOUT` (see
[the state layout](codegen.md#the-state-layout)): every configuration slot in variable
order, then every velocity slot.

| Variable | Configuration | Velocity | Constraint |
|---|---|---|---|
| `s: Scalar` | $s$ (1) | $\dot s$ (1) | |
| `x: Vector` | $x$ (3) | $\dot x$ (3) | |
| `q: S1` | $q$ (3) | $\omega_q$ (1) | $\lVert q \rVert = 1$, $q \cdot e = 0$ ($\dot q = \omega_q\, e \times q$, axis $e$ in `p`) |
| `q: S2` | $q$ (3) | $\omega_q$ (3) | $\lVert q \rVert = 1$, $q \cdot \omega_q = 0$ ($\dot q = \omega_q \times q$) |
| `R: SO3` | $R$ (9, row-major) | $\Omega_R$ (3, body frame) | $R^T R = I$, $\det R = 1$ ($\dot R = R \widehat{\Omega_R}$) |

`simulate` checks `x0` against these constraints (tolerance 1e-9) and raises
`ValueError` naming the slot, for example
`x0: the S2 slot 'q' must have unit norm, got |q| = 1.01` or
`x0: the velocity slot '\omega_{q}' must be orthogonal to 'q'`. `check_state(module,
state, tol, p=p)` runs the same check ($q \cdot e$ of an $S^1$ point only when `p` holds its
axis). A trajectory gives the values of a slot by its printed
symbol or its generated name: `traj["q"]` is `(N, 3)`, `traj["R"]` is `(N, 3, 3)`,
`traj["omega_q"]` is `(N, 3)`. The `StateLayout` the module was generated from builds a
state from named slots and splits one
([`pack`/`unpack`](codegen.md#the-state-layout)):

```python
x0 = layout.pack(q=q0, omega_q=w0)          # velocities left out are zero
layout.unpack(traj.state[-1])["omega_q"]    # (3,)
```

Every stored state is on the manifold (see the two integrators below), so
`traj.state[-1]` passes the check and can start the next segment, for example to switch
controllers.

## Accelerations and the S² solve

The generated equations read `mass @ a + bias + input_matrix @ u = 0` with
$M = -(\text{physical inertia})$ ([Conventions](conventions.md#signs)). For $\mathbb{R}$, $\mathbb{R}^3$ and $SO(3)$
rows the blocks are square and invertible for a physical system.

An $S^2$ equation is meaningful only on the tangent plane of its point
([Conventions](conventions.md#s2-equations-live-on-the-tangent-plane)), and the derived
mass block may be singular or not:

| System | $S^2$ mass block | Rank |
|---|---|---|
| quadrotor with a suspended load | $\widehat{q}\,\widehat{q}\, m_Q l^2 = -m_Q l^2 (I - q q^T)$ | 2 |
| spherical pendulum | $-m l^2 I$ (the engine already used $q \cdot \dot{\omega_q} = 0$) | 3 |

In both cases the unprojected $3 \times 3$ solve is not the equation of motion: it either fails or
gives $q \cdot \dot{\omega_q} \neq 0$, and the trajectory leaves the sphere. The missing equation is the
tangency of the acceleration: differentiating $q \cdot \omega_q = 0$ gives
$(\omega_q \times q) \cdot \omega_q + q \cdot \dot{\omega_q} = 0$, and the first term vanishes, so $q \cdot \dot{\omega_q} = 0$. `accelerations`
replaces each $S^2$ row block (all its columns, including couplings to other variables'
accelerations) by

$$
P_q \cdot \text{row} + s\, q q^T \ \text{on its own } \dot{\omega_q} \text{ column}, \qquad
\text{right-hand side } -P_q (\text{bias} + G u),
$$

with $P_q = I - q q^T$. Multiplying the new rows by $q^T$ gives $q \cdot \dot{\omega_q} = 0$; multiplying
by $P_q$ gives the tangent-plane equation. The test suite checks both (to 1e-12 and
1e-10) on every $S^2$ system of the numerical check gate, including the two of the table.

**Scaling.** The result must not depend on the units. The stabiliser weight is
$s = \lVert P_q D \rVert_F / \sqrt{2}$, with $D$ the row's own diagonal mass block. Both forms in the table
project to $P_q D = -m l^2 P_q$, and $\lVert P_q \rVert_F = \sqrt{2}$, so $s = m l^2$: the normal direction
gets the magnitude of the two tangent directions (for an anisotropic tangent inertia,
their root mean square). Then every row of the system and its right-hand side is divided
by the row's largest entry (row equilibration) and the system is solved with
`np.linalg.solve`. The solution is the same; the rounding and the condition number no
longer depend on how the masses and lengths compare with 1. With a unit stabiliser and no
equilibration, a double spherical pendulum with $m_2/m_1 = 10^{-9}$ lost 8 digits silently,
and a 1 t mass on a 100 km tether (SI units, $m l^2 = 10^{13}$) was rejected as singular. Both
are now solved to round-off (condition numbers 2.9 and 1.1).

A singular system (condition number above 1e12 after equilibration: a zero mass, a
singular inertia, a kinetic energy that ignores a velocity such as $\tfrac12 m \lvert \dot x + l \dot q \rvert^2$,
condition number ≈ 4e16) raises `np.linalg.LinAlgError` with the condition number and the
state. `accelerations` checks the condition number on every call; `simulate` checks it at
`x0` and afterwards only catches a failed or non-finite solve.

## Integrators

### `method="rkmk4"` (default): Runge–Kutta–Munthe-Kaas of order 4

A fixed-step Lie-group integrator. The state is moved by a group element
$\exp(U)$ of the product group that acts on it:

- each $S^2$ point by a rotation from the left, $q \leftarrow \exp(u)\, q$ (so $\dot q = \omega_q \times q$ has generator $\omega_q$);
- each $S^1$ point by a rotation about its axis `e`, $q \leftarrow \exp(u)\, q$ with $u$ a multiple of
  $e$ (generator $\omega_q e$; the axis is read from `p`). These rotations commute, so the step
  turns $q$ by the RK4 average of $\omega_q h$ exactly; $q$ is then re-projected onto the circle
  ($q \leftarrow q - (q \cdot e)\, e$, normalised), a round-off correction;
- each $SO(3)$ point by a rotation from the right, $R \leftarrow R \exp(u)$ (generator $\Omega_R$);
- every other slot ($\mathbb{R}$, $\mathbb{R}^3$ configuration and all velocities) by a translation.

One step of size $h$ applies classical RK4 to the ODE of $U$ in the Lie algebra,
$\dot U = \operatorname{dexp}^{-1}_U\big(A(\exp(U) \cdot y_n)\big)$, $U(t_n) = 0$ (Munthe-Kaas, *High order Runge–Kutta methods
on manifolds*, Appl. Numer. Math. 29, 1999; Hairer, Lubich, Wanner, *Geometric Numerical
Integration*, §IV.8):

$$
\begin{aligned}
&K_1 = h F(0), \quad K_2 = h F(K_1/2), \quad K_3 = h F(K_2/2), \quad K_4 = h F(K_3), \\
&U = (K_1 + 2K_2 + 2K_3 + K_4) / 6, \quad y_{n+1} = \exp(U) \cdot y_n, \\
&F(V) = \operatorname{dexp}^{-1}_V\big(A(\exp(V) \cdot y_n)\big) \quad \text{at } t,\ t + h/2,\ t + h/2,\ t + h.
\end{aligned}
$$

$\operatorname{dexp}^{-1}_V(G) = G - \tfrac12 [V, G] + \tfrac1{12} [V, [V, G]] + \dots$ is truncated after the $\operatorname{ad}^2$ term.
On $\mathfrak{so}(3) \cong \mathbb{R}^3$ the bracket is the cross product; for the right action on $SO(3)$ the sign of
the odd term flips ($\operatorname{dexp}^{-1}_{-V}$); on the vector slots every bracket vanishes and the step
is exactly classical RK4. Because $V = O(h)$ and the $\operatorname{ad}^3$ coefficient is zero, the
truncation changes each $K_i$ by $O(h^5)$ and the method keeps order 4. Without the
brackets the method drops to order 2 on a tumbling rigid body (error ratio 4.0 per
halving instead of 16).

$q$ and $R$ are only ever multiplied by rotation matrices, so $\lvert q \rvert - 1$ and $\lVert R^T R - I \rVert$
stay at round-off (≈ 5e-14 over 10 s) without projection. $q \cdot \omega_q = 0$ is a velocity
constraint that the dynamics preserve ($\tfrac{d}{dt} (q \cdot \omega_q) = 0$ from the $S^2$ solve), but $\omega_q$ is a
vector slot that is translated while $q$ is rotated, so the integrator alone would keep
it only to its truncation error. After each step $\omega_q \leftarrow \omega_q - (q \cdot \omega_q)\, q$, a change of the size
of the local error that keeps order 4, so $q \cdot \omega_q$ stays at round-off too.

`dt` is shortened, if needed, so that it divides $t_f - t_0$. The trajectory holds every
step, or every `output_every`-th step (and the last) with `simulate(...,
output_every=k)`; the integration is the same, only the storage is thinned (1000 s of the
quadrotor with a load at `dt = 1e-3` is about 190 MB in full).

### `method="solve_ivp"`

`scipy.integrate.solve_ivp` (`ivp_method="DOP853"` by default, `rtol`/`atol` passed
through) on the embedded state with $\dot y = \big(\mathtt{kinematics}(y), \mathtt{accelerations}(y)\big)$. The solver
does not know the manifold, so the outputs (at spacing `dt`, thinned by `output_every`)
are projected onto it: $q / \lvert q \rvert$, $\omega_q$ made tangent ($\omega_q - (q \cdot \omega_q)\, q$), an $S^1$ point moved into
the plane of its axis before normalising, and the polar factor
of $R$ ($U V^T$ of its SVD). `traj.constraint_drift` reports the drift *before*
projection: $\lvert q \rvert - 1$ reaches about 3e-10 over 10 s at `rtol = 1e-10`.

`traj.state`, and so `traj.energy()`, `traj["q"]` and the controller errors computed from
it, are the projected states. At `rtol = 1e-10` projection changes the energy by about
1e-10, below every tolerance of the table below; at a loose `rtol` it hides part of the
drift, so check `traj.constraint_drift` as well.

### Inputs

`u` is `None` (zero inputs), a constant 1-D array of `LAYOUT["n_inputs"]` entries (a
scalar for one input), or a callable `u(t, state)` evaluated at every stage (an
`LQRController` is one). A value of another shape, such as a `(2, 2)` array for four
inputs, raises `ValueError` rather than being flattened; so do non-finite entries
(`u has non-finite entries at [3]`, with the time for a callable). A callable that
returns `None`, typically a controller missing its `return`, raises `TypeError` naming
the time instead of running open loop.

A callable `u` must be a pure function of `(t, state)`. RKMK4 evaluates it at the stage
states (`t + h/2` twice, at two different states), and `solve_ivp` at trial steps that
it may reject and at times that are not monotone. A controller with internal state (an
integral term, a filter, a log) must not update that state inside `u`; integrate it as
part of the model instead, or update it between `simulate` calls of one step each.

## Energy and momentum checks

`traj.energy()` evaluates the module's `energy(state, p)` (generated with
`generate_module(..., lagrangian=L)`, see [Total energy](codegen.md#total-energy)) at every
output state; `traj.energy(fn)` uses another function `fn(state, p)`, such as
`energy_function(L, layout)`.

`tests/test_simulation_conservation.py` runs the systems of `tests/vbl_systems.py` for
10 s with zero inputs from the states of `tests/simulation_systems.py`. Relative energy
drift, $\max \lvert E(t) - E(0) \rvert / \lvert E(0) \rvert$:

| System | Tolerance | rkmk4, dt = 1e-3 | solve_ivp, rtol = 1e-10 |
|---|---|---|---|
| spherical pendulum | 1e-8 | 4.0e-13 | 6.5e-10 |
| free rigid body | 1e-9 | 3.0e-15 | 3.7e-11 |
| rigid pendulum about pivot | 1e-8 | 2.8e-14 | 4.8e-10 |
| double rigid pendulum | 1e-7 | 1.9e-14 | 9.4e-11 |
| elastic spherical pendulum | 1e-7 | 3.4e-12 | 8.5e-10 |
| quadrotor with suspended load | 1e-7 | 2.6e-11 | 5.2e-11 |

The inertial angular momentum $R J \Omega_R$ of the free rigid body is conserved to 1e-9 by both
methods. Two closed forms, written in NumPy without geomech and integrated with DOP853 at
`rtol = 1e-13`, agree with both methods to 1e-8 over 2 s: the spherical pendulum (Lee et
al. eq. 5.14, $m l^2 \dot{\omega_q} = m g l\, e_3 \times q$ for the potential $m g l\, e_3 \cdot q$, $\dot q = \omega_q \times q$, in
$\mathbb{R}^3 \times \mathbb{R}^3$) and the free rigid body (Euler's equations $J \dot{\Omega_R} = J \Omega_R \times \Omega_R$, $\dot R = R \widehat{\Omega_R}$).

## LQR from the linearization

`lqr(lin_module, x_eq, u_eq, p, Q, R, a_eq=None, module=None)` evaluates the generated `M_lin`,
`A_lin`, `B_lin` ([Linearization](vbl.md)) at an equilibrium, forms
$\bar A = M^{-1} A$, $\bar B = M^{-1} B$ in reduced coordinates, solves the continuous algebraic Riccati
equation (`scipy.linalg.solve_continuous_are`) and returns an `LQRController` with the
gain `K`, the Riccati solution `P`, the reduced `A`, `B`, the closed-loop eigenvalues and
the state error. Calling the controller, `ctrl(t, state)` $= u_\mathrm{eq} - K\, \delta x(\text{state})$, is the
input signature of `simulate`.

**One layout for both modules.** Generate the standard-form and the linearization module
from the same `StateLayout`, built with the same `SystemVariables` and the same `inputs`
list (same order) as `to_standard_form` and `linearize`; build `x_eq` and `x0` with
`layout.pack`:

```python
V = SystemVariables(vectors=[x_L, q], matrices=[R])
eom = compute_eom(L, dW, V)
layout = StateLayout(V, inputs=[f, M], params=collect_params([L, dW], V, [f, M]))
plant = load_module(generate_module(to_standard_form(eom, V, inputs=[f, M]), layout,
                                    lagrangian=L))
lin_module = load_module(generate_module(linearize(eom, V, inputs=[f, M]), layout))
x_eq = layout.pack(x_L=np.zeros(3), q=[0, 0, -1], R=np.eye(3))   # hover, at rest
u_eq = np.r_[(p["m_Q"] + p["m_L"]) * p["g"], 0, 0, 0]              # [f, M]
ctrl = lqr(lin_module, x_eq, u_eq, p, Q=1.0, R=1.0, module=plant)
```

`generate_module` raises `ValueError` for a linearization whose variables or inputs are
in another order than the layout's (the columns of `B_lin` would otherwise pair with the
wrong entries of `u`, and the loop diverges under a gain whose eigenvalues look stable),
and `lqr` checks that `LINEARIZATION` matches its `LAYOUT` and, with `module=`, that the
two modules have the same variables and inputs in the same order.

**The equilibrium is checked when `module=` is given.** The linearization module has no
residual, so by itself `lqr` cannot tell an equilibrium from any other point, and at a
non-equilibrium it still returns a gain whose closed-loop eigenvalues look stable: for the
rigid pendulum 30° from upright ($\dot{\Omega_R} = (2.46, -0.12, 0.008)$ there) the gain's
eigenvalues have real parts $\le -1.37$, and the nonlinear loop settles 0.22 rad away from
`x_eq` with a constant input. Pass the standard-form module of the same system
(`module=plant`) and `lqr` requires

- every velocity of `x_eq` to be zero (to `tol`; the $S^2$ reduction assumes $\omega_{q,\mathrm{eq}} = 0$),
- the equation residual $r = \mathtt{mass}\, a_\mathrm{eq} + \mathtt{bias} + G\, u_\mathrm{eq}$ ($S^2$ rows projected) to satisfy
  $\lVert r \rVert_\infty \le \mathtt{eq\_tol} \cdot \mathtt{scale}$ (`eq_tol = 1e-8`), where `scale` is the largest of $\lVert \mathtt{bias} \rVert_\infty$,
  $\lVert G\, u_\mathrm{eq} \rVert_\infty$ and the linearization's $\lVert A_\mathrm{lin} T \rVert_\infty$. The last term matters at an
  equilibrium held by a potential, where the bias itself is zero.

Otherwise it raises `ValueError` with the residual and `accelerations(module, x_eq,
u_eq, p)`. Without `module=`, `lqr` emits a `UserWarning` that the equilibrium was not
checked (an error under `-W error`).

The state error uses the perturbation coordinates of `linearize`:

| Variable | Configuration error | Velocity error |
|---|---|---|
| $\mathbb{R}$, $\mathbb{R}^3$ | $x - x_\mathrm{eq}$ | $\dot x - \dot x_\mathrm{eq}$ |
| $SO(3)$ | $\eta_R = \log(R_\mathrm{eq}^T R)^\vee$ (`geomech.numeric.log_so3`) | $\Omega_R - \Omega_{R,\mathrm{eq}}$ |
| $S^1$ | $\operatorname{atan2}\big(e \cdot (q_\mathrm{eq} \times q),\, q_\mathrm{eq} \cdot q\big)$, the angle about $e$ (1 entry) | $\omega_q - \omega_{q,\mathrm{eq}}$ (1 entry) |
| $S^2$ | $B^T (q_\mathrm{eq} \times q)$ (2 entries) | $B^T (\omega_q - \omega_{q,\mathrm{eq}})$ (2 entries) |

The $SO(3)$ error follows from $\delta R = R \widehat{\eta_R}$: $R = R_\mathrm{eq} \exp(\widehat{\eta_R})$. On $S^2$, $\delta q = \xi_q \times q$ with
$\xi_q \perp q$, so $\xi_q = q_\mathrm{eq} \times q$ to first order.

**The $S^2$ tangent basis.** The linearization of an $S^2$ point has 3-D perturbations $\xi_q$ and
$\delta\omega_q$, but at an equilibrium (where $\omega_q = 0$) both are orthogonal to $q$: the normal
direction is a zero mode that no input reaches, so the Riccati equation has no
stabilising solution in 3-D. `lqr` therefore works in the 2-D basis
`B = tangent_basis(q_eq)` (orthonormal, $b_1 \times b_2 = q_\mathrm{eq}$): with $T$ block diagonal ($B$
for each $S^2$ perturbation, $I$ elsewhere), it uses $T^T M T$, $T^T A T$, $T^T B$, so the $S^2$
rows are projected onto the tangent plane like the $S^2$ solve above. A singular $S^2$ block of
`M_lin` ($\widehat{q}\,\widehat{q}\, m l^2$ for the quadrotor with a load, where the full `M_lin` has
rank 17 of 18) becomes invertible. `Q` is given in the reduced coordinates
(`ctrl.labels` names them: `xi_q[0]`, `xi_q[1]`, `eta_R[0]`, …), as a scalar, a diagonal
or a matrix.

**Anisotropic weights on an $S^2$ block depend on the basis.** `xi_q[0]` and `xi_q[1]` are
the components along $b_1$ and $b_2$ of `tangent_basis(q_eq)`, and $b_1$ is the coordinate
axis least aligned with `q_eq`, projected. A tiny change of `q_eq` can switch that axis
and turn the basis by 90° (near $q_\mathrm{eq} = -e_3$, where the two smallest components tie:
`(2e-9, 1e-9, -1)` gives $b_1 = e_2$ instead of $e_1$), so different
weights on the two entries of an $S^2$ block weight different physical directions. Give each
$S^2$ block (both `xi_q` entries, and both `delta_omega_q` entries) equal weights, or check
`ctrl.bases[str(q)]` for the basis that was used. A scalar `Q` is always safe.

### The reduced linearization without a gain

`reduced_linearization(lin_module, x, u, p, *, a=None, basis=tangent_basis)` returns the
matrices `lqr` designs with, without solving a Riccati equation: a
`ReducedLinearization` with `A`, `B` ($\delta\dot x = A\,\delta x + B\,\delta u$ in the reduced coordinates),
`T`, `labels`, `bases` and the full `M_lin`, `A_lin`, `B_lin` at the point. `lqr` calls
it, so the two always agree (`ctrl.A`, `ctrl.B`, `ctrl.T`, `ctrl.labels`). Use it to check
controllability before designing (an uncontrollable mode on the imaginary axis or in the
right half-plane makes `lqr`'s Riccati solve fail), or to design another controller in the
same coordinates as `LQRController.error`. *basis* replaces `tangent_basis` for the $S^2$
blocks (any function returning an orthonormal $3 \times 2$ basis of the tangent plane).

The reduction assumes $\delta\omega_q \perp q$, a constant basis and a zero $S^2$ equation residual, which
hold at an equilibrium. Along a moving reference three terms are missing: a time-varying
basis adds $-T^T M \dot T$ to the reduced $A$; the normal component of the velocity
perturbation, $q \cdot \delta\omega_q = -(q \times \omega_q) \cdot \xi_q$, couples the velocity rows to $\xi_q$; and the $S^2$ equation
$B^T E = 0$, satisfied by the reference with a normal residual $E_r = q\, (q \cdot E_r)$, gains
$(q \cdot E_r)\, B^T \widehat{q} B$ in its $\xi_q$ columns when it is linearized (from $\delta B^T q = -B^T \delta q$).
`examples/control/tracking.py` builds that time-varying version and checks it against
finite differences of the nonlinear error dynamics. `model_of(module)`
(also exported) parses a generated module's `LAYOUT` into its variables and slices.

### Example: the inverted rigid pendulum

`examples/inverted_rigid_pendulum_lqr.py` derives a rigid pendulum about its pivot on
$SO(3)$ (Lee et al. eq. 6.8), generates its standard form and linearization, designs an LQR
gain about the inverted equilibrium $R_\mathrm{eq} = \operatorname{diag}(1, -1, -1)$ and simulates the closed
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

The tests also stabilise an inverted spherical pendulum ($q$ 15° from $+e_3$, the maximum
of $m g l\, e_3 \cdot q$, with its tangent input) to $\lvert q - e_3 \rvert < 10^{-4}$ in 5 s, and the quadrotor
with a suspended load about hover.

## Not covered

Contact, constraints beyond the manifold, variable-step Lie-group integrators, MPC and
real-time performance. Each RKMK4 stage costs one call of `mass`, `bias` and
`input_matrix` plus a small dense solve: 10 s at `dt = 1e-3` (40 000 stages) takes about
4.6 s for the double rigid pendulum, whose generated blocks are the slowest, and less for
the single-body systems. `traj.energy()` calls the generated `energy` once per stored
state (about 0.23 s for those 10 001 states, 5 % of the run), so it is not vectorised.
