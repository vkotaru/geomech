# Closed-loop control examples

LQR stabilisation of equilibria of the systems in `geomech.systems`, designed from the
variation-based linearization (`linearize`) and simulated with the nonlinear generated
model on the manifold. The systems are used as they are in the library.

```bash
python examples/control/run.py --list
python examples/control/run.py quadrotor_suspended_load_hover
python examples/control/run.py --all        # about 3 min; rewrites results/ and table.md
python examples/control/tracking.py --all   # time-varying LQR tracking (stretch), about 1 min
```

Needs scipy (`pip install "geomech[sim]"`).

## The chain (`control_helper.py`)

1. `compute_eom` → `to_standard_form` and `linearize`; NumPy code for both from one
   `StateLayout`.
2. The equilibrium: a configuration chosen per example (velocities zero) and the input
   `u_eq` that solves `bias + G u = 0` in the least-squares sense (S² rows on their
   tangent planes). The residual is at round-off (≤ 1e-15) for every example.
3. The reduced linearization `δẋ = A δx + B δu` from
   `geomech.simulation.reduced_linearization`, in the coordinates of `lqr`: S²
   perturbations in a 2-D tangent basis, S¹ errors as the angle about the axis, SO(3)
   errors as `vee(log(R_eqᵀR))`.
4. Controllability: an orthonormal basis of the controllable subspace (block Krylov
   iteration) and the eigenvalues of `A` on its complement (the uncontrollable modes),
   cross-checked with the PBH test. A fully controllable system is designed with the
   library's `lqr` (which also checks the equilibrium against the nonlinear model);
   otherwise LQR is designed on the controllable subspace.
5. `u = u_eq − K δx(state)`, no input saturation, simulated with RKMK4 at dt = 1 ms.
   Closed-loop eigenvalues are kept inside the RKMK4 stability region
   (`|λ|max · dt ≤ 0.91` for every example; the stiffest are the links with a small
   axial inertia).
6. Metrics in `results/<name>.json`, the tables below in `results/table.md` (regenerated
   from the JSON files by `run.py --all`; a test checks it); time histories (100 Hz) in
   `results/<name>.npz`, which `examples/control/.gitignore` keeps out of git (about
   5 MB; `--all` regenerates them).

**Metrics.** `n`: reduced state dimension; `rank`: dimension of the controllable
subspace; `max Re λ`: the slowest controlled closed-loop eigenvalue; `min ζ`: the smallest
damping ratio `−Re λ/|λ|` of the controlled modes; `e₀ pos.`: the norm of the initial
position error (R³ and length coordinates, m); `e₀ att.`: the norm of the initial
attitude error (S¹, S² and SO(3) coordinates, rad). Convergence uses the normalized error
`ν(t) = √(V(z(t)) / V(z(0)))`, with `V = zᵀPz` the LQR cost-to-go (`P` the Riccati
solution on the controllable subspace): a weighted norm of the whole reduced error that
the linear closed loop decreases monotonically. `t_s` is the last exit time of ν from the
2 % band; converged means `ν ≤ 2 %` over the whole last 10 % of the horizon. The peak
input is the largest deviation from `u_eq` of any component, per unit (N for forces and
thrusts, N·m for moments and torques). Drift is the largest manifold constraint violation
(`|q| − 1`, `‖RᵀR − I‖`) of the run.

## Results

| Example | n | rank | max Re λ | min ζ | e₀ pos. (m) | e₀ att. (rad) | t_s | peak ‖u − u_eq‖∞ | drift | converged |
|---|---|---|---|---|---|---|---|---|---|---|
| `spherical_pendulum_inverted` | 4 | 4 | -3.17 | 1.00 | — | 0.5 | 1.39 s | 9.68 N·m | 7.0e-15 | yes |
| `double_spherical_pendulum_inverted` | 8 | 8 | -2.55 | 1.00 | — | 0.312 | 1.60 s | 6.33 N·m | 9.9e-15 | yes |
| `n_link_spherical_chain_inverted` | 12 | 12 | -2.17 | 1.00 | — | 0.305 | 1.88 s | 7.72 N·m | 1.5e-14 | yes |
| `elastic_spherical_pendulum_inverted` | 6 | 6 | -0.559 | 0.09 | 0.1 | 0.342 | 4.38 s | 0.527 N, 5.29 N·m | 7.1e-15 | yes |
| `cart_planar_pendulum_inverted` | 4 | 4 | -1.51 | 0.79 | 0.3 | 0.175 | 3.01 s | 21.4 N | 1.1e-16 | yes |
| `rigid_pendulum_inverted` | 6 | 6 | -3.16 | 1.00 | — | 0.436 | 1.22 s | 3.8 N·m | 4.1e-14 | yes |
| `double_rigid_pendulum_inverted` | 12 | 12 | -2.97 | 1.00 | — | 0.315 | 1.50 s | 4.32 N·m | 4.1e-14 | yes |
| `n_link_rigid_chain_inverted` | 18 | 18 | -2.42 | 1.00 | — | 0.306 | 1.57 s | 5.87 N·m | 5.1e-14 | yes |
| `free_rigid_body_attitude` | 6 | 6 | -1.66 | 0.78 | — | 1.92 | 2.27 s | 4.96 N·m | 1.1e-13 | yes |
| `spacecraft_reaction_wheels_slew` † | 9 | 6 | -0.371 | 0.72 | — | 1.05 | 10.15 s | 0.191 N·m | 7.0e-14 | yes |
| `spacecraft_reaction_wheels_s1_slew` † | 9 | 6 | -0.371 | 0.72 | — | 1.05 | 10.15 s | 0.191 N·m | 7.0e-14 | yes |
| `spacecraft_solar_arrays_slew` † | 13 | 10 | -0.281 | 0.73 | — | 1.44 | 13.52 s | 0.191 N·m | 1.3e-13 | yes |
| `biped_3d_balance` | 12 | 12 | -0.882 | 0.21 | — | 0.145 | 1.87 s | 20.4 N·m | 6.9e-14 | yes |
| `rigid_body_pose` | 12 | 12 | -0.924 | 0.73 | 1.22 | 1.57 | 3.93 s | 3.16 N, 4.97 N·m | 5.6e-14 | yes |
| `quadrotor_hover` | 12 | 12 | -0.924 | 0.71 | 1.87 | 1.05 | 3.78 s | 1.58 N, 4.76 N·m | 9.5e-14 | yes |
| `quadrotor_offset_com_hover` | 12 | 12 | -0.923 | 0.68 | 1.22 | 0.785 | 3.64 s | 1.53 N, 3.23 N·m | 1.3e-13 | yes |
| `quadrotor_suspended_load_hover` | 16 | 16 | -0.596 | 0.36 | 0.781 | 0.92 | 4.49 s | 0.761 N, 2.33 N·m | 7.3e-14 | yes |
| `quadrotor_elastic_load_hover` | 18 | 18 | -0.596 | 0.04 | 0.473 | 0.406 | 4.23 s | 2.51 N, 0.752 N·m | 9.7e-14 | yes |
| `quadrotor_flexible_cable_hover` | 24 | 24 | -0.166 | 0.01 | 0.374 | 0.203 | 21.70 s | 0.12 N, 0.808 N·m | 1.4e-13 | yes |
| `multi_quadrotor_point_mass_hover` | 36 | 36 | -0.45 | 0.48 | 0.412 | 0.139 | 6.79 s | 0.221 N, 0.276 N·m | 1.5e-13 | yes |
| `multi_quadrotor_rigid_load_hover` | 42 | 42 | -0.438 | 0.17 | 0.412 | 0.281 | 6.81 s | 0.143 N, 0.321 N·m | 1.3e-13 | yes |
| `two_rigid_bodies_spherical_joint_hold` | 18 | 18 | -1.03 | 0.74 | 0.412 | 0.554 | 3.64 s | 0.967 N, 1.11 N·m | 9.5e-14 | yes |

† Not stabilizable in the strict sense: uncontrollable modes with eigenvalue 0 (conserved quantities). LQR on the controllable subspace converges from initial states on the invariant subspace (for the spacecraft, zero total angular momentum).

**Not stabilizable** (uncontrollable modes listed; LQR on the controllable subspace):

| Example | n | rank | uncontrollable eigenvalues | max Re λ (controlled) | e₀ att. (rad) | t_s | converged |
|---|---|---|---|---|---|---|---|
| `cart_spherical_pendulum_inverted_in_plane` | 6 | 4 | -3.5, +3.5 | -1.51 | 0.174 | 3.01 s | yes (in-plane start only) |
| `cart_spherical_pendulum_out_of_plane` | 6 | 4 | -3.5, +3.5 | -1.51 | 0.174 | — | diverged |

### Notes per example

- **Planar cart-pole.** `cart_planar_pendulum_inverted` is the in-plane run of the 3-D
  cart-pole on `cart_planar_pendulum` (R × S¹), with the same start (0.3 m, 10° tilt),
  weights, `R = 0.1` and horizon. Unlike the 3-D cart-pole, whose *linearization* is not
  stabilizable (the sideways swing is uncontrollable to first order, with eigenvalues
  ±√(g/l)), the planar cart-pole is fully controllable (rank 4 of 4) and the library's
  `lqr` applies directly. Its closed-loop eigenvalues are those of the controllable part
  of the 3-D design (max Re λ = −1.51, min ζ = 0.79) and it settles in the same
  t_s = 3.01 s, final ν = 2.2e-4; the peak rail force is 21.4 N (21.3 N in 3-D, where the
  S² tangent coordinate of the 10° tilt is sin 10° rather than 10°).
- **Cart-pole (3-D): the linearization is not stabilizable.** The cart force on a rail
  along e₁ exerts no torque about the rail axis, so at the upright equilibrium the
  pendulum's swing out of the rail's vertical plane is uncontrollable to first order, with
  eigenvalues ±√(g/l) = ±3.5 s⁻¹ (PBH and the controllable-subspace rank 4 of 6 agree).
  Every smooth (C¹) state feedback therefore leaves the eigenvalue +3.5 in the closed-loop
  linearization and the upright equilibrium unstable (Lyapunov's indirect method); LQR, or
  any controller designed on the linearization, cannot stabilize it. The force does reach
  the sideways swing through nonlinear coupling (once the pendulum is tilted both in and
  out of the plane, the rail force changes the rate about the vertical and so the
  sideways angle), so this does not rule out nonsmooth or time-varying nonlinear
  controllers; none is tried here. `cart_spherical_pendulum_inverted_in_plane` designs LQR
  on the controllable (in-plane) subspace and starts in the plane. It converges because
  the out-of-plane component starts at exactly zero and stays zero by symmetry (the run
  also converges over 12 s); any out-of-plane perturbation grows under this controller.
  `cart_spherical_pendulum_out_of_plane` adds a 1°
  out-of-plane tilt to the same run: it diverges (the simulation stops at t ≈ 3.9 s).
- **Spacecraft with reaction wheels (†).** The wheel angles are cyclic (their columns of
  `A` and `B` are zero) and are left out of the design; with the wheel torques internal,
  the three components of the total angular momentum are conserved, so three modes
  (eigenvalue 0) are uncontrollable. LQR is designed on the 6-D controllable subspace.
  The rest-to-rest start (zero total momentum) is essential, not a convenience: internal
  torques cannot change the total angular momentum, so a spacecraft that starts with
  momentum cannot come to rest at `R = I` with wheel torques alone. With `R = 1e4` the
  60° slew settles in 10.2 s with a peak wheel torque of 0.19 N·m (no saturation is
  modelled).
- **Spacecraft with S¹ wheels (†).** `spacecraft_reaction_wheels_s1_slew` is the slew
  above on `spacecraft_reaction_wheels_s1`, with the same weights, `R = 1e4` and horizon.
  The wheel directions qᵢ are cyclic like θᵢ (their angle errors γᵢ are dropped), and the
  reduced linearization is the R³ one with its coordinates reordered, so the result is
  the same: t_s = 10.15 s, final ν = 6.7e-5, peak wheel torque 0.191 N·m, identical in
  every reported digit to `spacecraft_reaction_wheels_slew`.
- **Spacecraft with solar arrays (†).** `spacecraft_solar_arrays_slew`: the same 60° slew
  on `spacecraft_solar_arrays`, with the three wheel torques and the two array-drive
  torques, ending with both panel normals turned 40° about b₁ from b₃ (the commanded
  sun-pointing angle; the panels start at 0°). The wheel angles are dropped as above; the
  panel angles are kept and controlled through the drive torques. The internal torques
  conserve the total angular momentum: rank 10 of 13, three uncontrollable modes at
  eigenvalue 0, so the system is stabilizable only on the zero-momentum subspace, which
  the rest-to-rest start lies on. It converges: t_s = 13.5 s, final ν = 2.0e-4, the
  panels at 40.0° (error below 1e-4 rad). The drive torques are weighted 10× the wheel
  torques and the panel rates 100: peak wheel torque 0.191 N·m, peak drive torque
  0.008 N·m (the table shows the larger, per unit). The panels' rate relative to the bus
  peaks at 0.23 rad/s, mostly the bus turning under them while the weak drives hold them.
- **Lightly damped modes.** The elastic cable of `quadrotor_elastic_load_hover` has a
  closed-loop cable-spring mode near 28 rad/s with damping ratio 0.015 under uniform
  weights; weighting the stretch rate 10 gives 0.039, which is used. Heavier weights damp
  it more (0.12, 0.36) but need thrust deviations of 12–28 N against a hover thrust of
  17.7 N; the chosen weights and a smaller perturbation keep the peak thrust deviation at
  2.5 N (14 %) and the peak moment at 0.75 N·m. The three-link cable of
  `quadrotor_flexible_cable_hover` keeps two modes near 25 rad/s lightly damped (ratio
  0.007 with the moments weighted R = 10, which keeps the peak moment at 0.81 N·m);
  the error rings, and the 40 s horizon is long enough for the tail-window criterion
  to be met honestly. The other small damping ratios (elastic pendulum 0.09, plate team
  0.17, biped 0.21) are left as LQR gives them. No input saturation is modelled in any
  example: the peaks are what the linear gains ask for.
- **Biped.** A static posture of the swing-phase model: stance leg vertical, swing leg
  held 20° forward; `u_eq` is the hip torque 8.39 N·m about e₂ that holds the swing leg
  (and its reaction on the stance leg, cancelled by the ankle torque). The stance leg
  starts 5° and the swing leg 6.6° from the equilibrium. Contact is not modelled: the
  stance foot is a spherical joint able to apply any ankle torque, and the swing foot may
  pass through the ground. This is balance of the model, not walking.
- **Elastic pendulum upright** sits at the compressed spring length `l = l₀ − m g / k`.
- **Quadrotor with an offset centre of mass** hovers level with thrust `m g` and the
  constant moment `M = m g ρ × e₃ = (−0.196, −0.392, 0)` N·m from `u_eq`.
- **Teams of quadrotors**: symmetric equilibria with the cables at 0.5 rad (point mass)
  and 0.3 rad (rigid plate) from the vertical; each quadrotor's thrust direction balances
  its weight and the cable tension `m_L g / (n cos α)` (`team_equilibrium` in `run.py`).
- **Two rigid bodies**: body 1 held level in mid-air, body 2 hanging below the joint,
  force `F = (m₁ + m₂) g e₃` on body 1 and the moment that the offset joint load needs.

## Time-varying LQR tracking (stretch, `tracking.py`)

The Riccati differential equation is integrated backward with scipy along the reference,
with the linearization module evaluated at each reference point
(`reduced_linearization(..., a=a_ref, basis=...)`) and three corrections for a moving S²
reference (derived in the module docstring): the time-varying basis (`−TᵀMṪ`), the normal
component of the velocity perturbation (`q·δω = −(q × ω)·ξ`) and the row projection of
the equation's normal residual (`(q·E_r) Bᵀq̂B`). Errors: position (m) and attitude
`‖log(R_refᵀR)‖` (rad) separately.

| Example | reference | e₀ pos. / att. | RMS pos. error (2nd half) | final pos. / att. error | peak ‖u − u_ref‖∞ |
|---|---|---|---|---|---|
| `quadrotor_figure_eight` | figure-eight by flatness, 8 s period, ≤ 1.67 m/s | 0.52 m / 0.35 rad | 5.2e-05 m | 7.2e-08 m / 4.1e-09 rad | 0.92 N (f), 1.14 N·m (M) |
| `suspended_load_circle` | simulated nominal on a 1 m circle, 5 s period, ≤ 0.95 m/s | 0.41 m / 0.35 rad | 9.3e-05 m | 1.4e-06 m / 2.2e-06 rad | 0.77 N (f), 1.26 N·m (M) |

**Flatness reference (figure-eight).** Analytic derivatives of the flat outputs (jerk,
snap) and closed-form `Ω_d = vee(R_dᵀṘ_d)`, `Ω̇_d = vee(skew(R_dᵀR̈_d))`. Independent
checks: `max ‖Ṙ_d − R_dΩ̂_d‖` with `Ṙ_d` from 5-point differences of `R_d`:
8.3e-13; `max |∫Ω̇_d dt − ΔΩ_d|` (trapezoidal, h = 1e-4):
1.3e-09; for `max |Ω_d|` = 0.315 rad/s and
`max |Ω̇_d|` = 0.460 rad/s².

**Reduced A(t) along the suspended-load reference** (`validate_reduction`): central
differences of the nonlinear error dynamics at four reference samples
(`max |ω_q|` = 0.22 rad/s):

| Reduced A(t) | max relative difference from finite differences |
|---|---|
| equilibrium reduction only | 3.1e-02 |
| + basis rate (−TᵀMṪ) | 3.1e-02 |
| + normal velocity coupling | 1.4e-02 |
| + curvature term (all) | 4.0e-07 |

At this speed the basis-rate term (`−TᵀMṪ`) changes nothing visible (3.1e-2 with and
without it: the basis turns slowly); the improvement, from 3.1e-2 to 4.0e-7, comes from the
normal-velocity coupling (to 1.4e-2) and the curvature term (to 4.0e-7). The load reference is
a trajectory of the nonlinear model itself (the hover LQR driven by a set point on a 1 m
circle, recorded with its inputs), so it is feasible by construction; it is not the
differential-flatness construction of Sreenath, Lee and Kumar (2013), which is not
implemented here.
