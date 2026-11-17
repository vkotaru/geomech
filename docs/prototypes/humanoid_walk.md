# Walking humanoid

```{note}
Prototype: built with geomech, not part of `geomech.systems` and not covered by the test
suite. The gait planner, the controller, the switching between stance models and the
foot-strike impacts are prototype code; geomech derives the three models they use. The
masses are illustrative.
```

## Animation

Six steps on flat ground, 0.30 m long and 0.54 s each: the pelvis moves 1.80 m forward in
3.25 s. Drag to orbit; the full report (equations of the floating model, standard form,
parameters) opens in a <a href="../_static/prototypes/humanoid_walk.html">separate page</a>.

<iframe class="model-report" src="../_static/prototypes/humanoid_walk.html?view" title="Walking humanoid, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Model

Twelve rigid bodies, every one on $SO(3)$: the ten of the
[floating humanoid](humanoid.md) and two flat feet on ankle ball joints, 56.7 kg in all.
geomech derives three models from the same Lagrangian (a sum over the bodies of
$\tfrac12 m_b \|\dot x_b\|^2 + \tfrac12 \Omega_b \cdot J_b \Omega_b - m_b g\, e_3 \cdot x_b$):

- **Floating**, for the foot strikes: pelvis position $x \in \mathbb{R}^3$ and all twelve
  attitudes, $\mathbb{R}^3 \times SO(3)^{12}$, 39 degrees of freedom.
- **Left stance** and **right stance**: the stance foot is flat and fixed, so it is part
  of the ground, and the stance ankle is a ball joint at the origin. The positions are
  built inward from the ankle, so $x$ is no longer a variable: $SO(3)^{11}$, 33 degrees of
  freedom. The inputs are the ankle moment from the ground and the ten joint moments,
  33 components: the model is fully actuated.

Each joint moment enters the virtual work on the child and, with the opposite sign, on
the parent, as for the floating humanoid. In the prototype, computed-torque control
tracks a planned gait (a linear inverted pendulum for the pelvis), each stance runs until
the swing foot reaches the ground, and the strike is a plastic flat-foot impact solved on
the floating model; the next stance model starts from the velocities after the impact.

## Checks

From the run behind the animation:

- `check_eom`: $1.3 \times 10^{-7}$ (floating model), $1.7 \times 10^{-7}$ (left stance),
  $8.1 \times 10^{-8}$ (right stance).
- Every stance phase: the centre of pressure stays inside the sole, the friction ratio
  stays at most 0.19, the vertical ground force at least 677 N, and the change of energy
  equals the work of the joint and ankle moments to $1.1 \times 10^{-4}$ J.
- Every foot strike: the new foot's velocity constraint holds to $3 \times 10^{-16}$, the
  linear and angular momentum balances with the ground impulse to $3 \times 10^{-14}$, and
  0.60 J of kinetic energy is lost.
