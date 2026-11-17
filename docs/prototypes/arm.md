# Six-joint arm

```{note}
Prototype: built with geomech, not part of `geomech.systems` and not covered by the test
suite. geomech derives the unconstrained equations of the links and rotors; the hinges and
gears are constraints that the prototype code adds and solves, and the controller is
prototype code too. The parameters are illustrative.
```

## Animation

A pick-and-place sequence of 7.5 s: five moves of 1.5 s between joint-space waypoints,
each a smooth quintic, carrying a 1 kg payload. Drag to orbit; the full report (the
unconstrained equations, standard form, parameters) opens in a
<a href="../_static/prototypes/arm.html">separate page</a>.

<iframe class="model-report" src="../_static/prototypes/arm.html?view" title="Six-joint arm, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Model

Six links, each a rigid body with attitude $R_k \in SO(3)$, and six motor rotors, each a
scalar spin angle $\psi_k$. Link $k$ starts at the end of link $k-1$ (the first at the
origin) and its centre of mass is at the joint plus $\tfrac12 R_k l_k$. Rotor $k$ sits on
link $k-1$ (on the ground for the first) and spins about the joint axis $a_k$, so its
angular velocity is $\Omega_{k-1} + a_k \dot\psi_k$. The Lagrangian sums, over the links
and rotors,

$$
\tfrac12 m_k \|\dot x_k\|^2 + \tfrac12 \Omega_k \cdot J_k \Omega_k - m_k g\, e_3 \cdot x_k
+ \tfrac12 (\Omega_{k-1} + a_k \dot\psi_k) \cdot J_{r,k} (\Omega_{k-1} + a_k \dot\psi_k),
$$

plus the translational energy and weight of the rotors, and the motor torques do the
virtual work $\delta W = \sum_k \tau_k\, \delta\psi_k$. geomech derives the equations on
$SO(3)^6 \times \mathbb{R}^6$, 24 coordinates.

The prototype then adds 18 constraint rows: two per hinge, $(R_{k-1} b) \cdot (R_k a_k) = 0$
for the two directions $b$ perpendicular to the axis, and one per gear,
$\dot\psi_k = n_k\, a_k \cdot (\Omega_k - \Omega_{k-1})$, with $n_k = 100$ for the first
three joints and 50 for the last three. That leaves 6 degrees of freedom. The constraint
forces are solved together with the equations (a KKT system) inside geomech's Lie-group
RK4 integrator, and the state is projected back onto the constraints every 10 ms. The
controller is PD with gravity compensation, applied through the gears.

## Checks

`check_eom` of the unconstrained equations: $1.3 \times 10^{-7}$ (12 equation blocks).
From the stored trajectories of the prototype's runs:

- Pick and place: joint tracking error at most 0.58°, motor torque at most 0.91 N·m,
  constraint residual after each projection at most $5 \times 10^{-14}$.
- A separate 4 s passive swing (no motor torque) keeps its energy, 52.8 J, to
  $2.6 \times 10^{-6}$ J.
