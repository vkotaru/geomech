# Floating humanoid

```{note}
Prototype: built with geomech, not part of `geomech.systems` and not covered by the test
suite. The joint controller is prototype code; the parameters are illustrative.
```

## Animation

Eight seconds in zero gravity. PD joint moments windmill both arms at 0.75 revolutions
per second (with knees and hips held tucked); since the moments are internal, the whole
body turns the other way. Drag to orbit; the full report (equations, standard form,
parameters) opens in a <a href="../_static/prototypes/humanoid_float.html">separate page</a>.

<iframe class="model-report" src="../_static/prototypes/humanoid_float.html?view" title="Floating humanoid, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Model

Ten rigid bodies, every one on $SO(3)$: pelvis, torso, two upper arms, two forearms, two
thighs and two shanks. The pelvis is the floating base, with centre of mass
$x \in \mathbb{R}^3$ and attitude $R_p$; the torso hangs on a waist ball joint, the arms
on shoulder and elbow ball joints, the legs on hip and knee ball joints. The configuration
manifold is $\mathbb{R}^3 \times SO(3)^{10}$, 33 degrees of freedom.

Every position is built outward from the pelvis: a joint point is the parent's joint point
plus $R_{\text{parent}}$ times a constant offset, and a body's centre of mass is its joint
point plus $R_b \rho_b$. The Lagrangian is the sum over the ten bodies

$$
L = \sum_b \tfrac12 m_b \|\dot x_b\|^2 + \tfrac12 \Omega_b \cdot J_b \Omega_b
  - m_b g\, e_3 \cdot x_b ,
$$

with $\dot R_b = R_b \widehat{\Omega}_b$. Each of the nine joints carries an inertial-frame
moment $M_c$ that acts on the child and, with the opposite sign, on the parent:

$$
\delta W = F \cdot \delta x + \sum_{\text{joints}}
  \eta_c \cdot R_c^T M_c - \eta_p \cdot R_p^T M_c ,
$$

with $\eta_b$ the variation of $R_b$ ($\delta R_b = R_b \widehat{\eta}_b$). There are
30 joint-moment components and a force $F$ on the pelvis, which is zero here.

## Checks

From the run behind the animation (gravity $g = 0$):

- `check_eom`: $2.8 \times 10^{-7}$ (11 equation blocks).
- The total linear and angular momentum stay zero: at most $3.8 \times 10^{-11}$ kg·m/s
  and $4.8 \times 10^{-12}$ kg·m²/s over the run. The joint moments are internal and
  there is no gravity, so both are conserved.
- After 8 s the pelvis has pitched about 100°.
