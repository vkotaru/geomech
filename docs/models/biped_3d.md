# 3-D biped

The swing phase of a three-dimensional compass-gait biped. The stance leg is pinned at
the stance foot by a spherical joint, the swing leg hangs from the hip by a second
spherical joint, and a point mass sits at the hip. The configuration is the pair of leg
attitudes $(R_1, R_2) \in SO(3) \times SO(3)$, with body angular velocities
$\Omega_1, \Omega_2$ and $\dot R_k = R_k \widehat{\Omega_k}$; no Euler angles, so no
singular attitude. Only this continuous single-support phase is modelled: there are no
foot impacts, no switching of the stance foot and no ground contact, so it is not a
walking model. Not taken from the literature; the equations are verified numerically
(`check_eom`). Library system `biped_3d`.

## Animation

One planar step, unforced, 0.775 s, in the $e_1$–$e_3$ plane, from a start that is not
the library's initial state (that one is not planar). The stance leg starts
0.25 rad behind the vertical and vaults forward over the stance foot; the swing foot
starts on the ground 0.5 m behind and swings forward. The hip stays at least 0.96 m
up and moves from 0.25 m behind the stance foot to 0.28 m in front of it. There is no
contact model: mid-swing the swing foot scuffs up to about 9 mm below the ground plane
(at about 0.14 s), then rises up to 0.13 m. The run ends at touch-down, with the swing
foot back on the ground 0.57 m in front of the stance foot; the impact and the switch
of stance leg are not modelled. Drag to orbit; the full report (equations, standard
form, parameters) opens in a <a href="../model-reports/biped_3d.html">separate page</a>.

<iframe class="model-report" src="../model-reports/biped_3d.html?view" title="3-D biped, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 280" role="img" aria-labelledby="bp-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="bp-title">3-D biped: stance leg R1 and swing leg R2, leg length l, centre-of-mass distance a, leg mass m and inertia J, hip mass m_H, angular velocities Omega1 and Omega2, hip torque M_H, ankle torque M_A and gravity</title>
  <defs>
    <marker id="bp-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- ground (drawn only, no contact is modelled) -->
    <line x1="40" y1="250" x2="340" y2="250" stroke-dasharray="5 4" opacity="0.6"/>
    <!-- inertial axes -->
    <line x1="50" y1="90" x2="50" y2="50" marker-end="url(#bp-arr)"/>
    <line x1="50" y1="90" x2="90" y2="90" marker-end="url(#bp-arr)"/>
    <!-- stance leg (foot to hip) and swing leg (hip to foot) -->
    <line x1="150" y1="250" x2="191.4" y2="95.4" stroke-width="3"/>
    <line x1="191.4" y1="95.4" x2="246.1" y2="245.8" stroke-width="3"/>
  </g>
  <!-- angular velocities of the legs -->
  <path d="M 171.1,204.7 A 50 50 0 0 1 188.3,217.9" fill="none" stroke="#3f8fa0" stroke-width="1.8" marker-end="url(#bp-arr)"/>
  <path d="M 221.4,147.4 A 60 60 0 0 0 240.5,129.8" fill="none" stroke="#3f8fa0" stroke-width="1.8" marker-end="url(#bp-arr)"/>
  <!-- hip torque, ankle torque, gravity -->
  <path d="M 170.7,87.9 A 22 22 0 0 1 183.9,74.7" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#bp-arr)"/>
  <path d="M 126.4,245.8 A 24 24 0 0 1 134.6,231.6" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#bp-arr)"/>
  <line x1="330" y1="40" x2="330" y2="90" stroke="#e08a2c" stroke-width="2" marker-end="url(#bp-arr)"/>
  <circle cx="150" cy="250" r="4" fill="currentColor"/>
  <circle cx="246.1" cy="245.8" r="4" fill="currentColor"/>
  <circle cx="170.7" cy="172.7" r="3" fill="currentColor"/>
  <circle cx="218.8" cy="170.6" r="3" fill="currentColor"/>
  <circle cx="191.4" cy="95.4" r="10" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="146" y="158" fill="#3f8fa0">R<tspan font-size="12" dy="4">1</tspan></text>
    <text x="242" y="206" fill="#3f8fa0">R<tspan font-size="12" dy="4">2</tspan></text>
    <text x="196" y="228" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">1</tspan></text>
    <text x="250" y="128" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">2</tspan></text>
    <text x="160" y="134">l</text>
    <text x="138" y="212">a</text>
    <text x="120" y="182">m, J</text>
    <text x="236" y="176">m, J</text>
    <text x="208" y="94">m<tspan font-size="12" dy="4">H</tspan></text>
    <text x="134" y="72" fill="#e08a2c">M<tspan font-size="12" dy="4">H</tspan></text>
    <text x="96" y="238" fill="#e08a2c">M<tspan font-size="12" dy="4">A</tspan></text>
    <text x="338" y="80" fill="#e08a2c">g</text>
    <text x="44" y="43">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="95" y="95">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>The stance foot is fixed at the origin. Each leg has length <span class="math notranslate nohighlight">\(l\)</span> along its
body axis <span class="math notranslate nohighlight">\(e_3\)</span>, which points from the foot to the hip, so the hip is at
<span class="math notranslate nohighlight">\(l R_1 e_3\)</span> and the swing foot at <span class="math notranslate nohighlight">\(l R_1 e_3 - l R_2 e_3\)</span>. Each leg has mass
<span class="math notranslate nohighlight">\(m\)</span> with its centre of mass (small dots) a distance <span class="math notranslate nohighlight">\(a\)</span> from its foot and
inertia <span class="math notranslate nohighlight">\(J\)</span> about it; the hip carries a point mass <span class="math notranslate nohighlight">\(m_H\)</span>. The hip
torque <span class="math notranslate nohighlight">\(M_H\)</span> acts on the swing leg and <span class="math notranslate nohighlight">\(-M_H\)</span> on the
stance leg; the ankle torque <span class="math notranslate nohighlight">\(M_A\)</span> acts on the stance leg from the
ground. Both are inertial-frame moments. Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span>. The dashed ground is
not a constraint.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $R_1$ | attitude of the stance leg, body frame to inertial frame | $SO(3)$ |
| $R_2$ | attitude of the swing leg, body frame to inertial frame | $SO(3)$ |
| $\Omega_1$ | body angular velocity of the stance leg, $\dot R_1 = R_1 \widehat{\Omega_1}$ | $\mathbb{R}^3$ |
| $\Omega_2$ | body angular velocity of the swing leg, $\dot R_2 = R_2 \widehat{\Omega_2}$ | $\mathbb{R}^3$ |
| $\eta_1$ | variation of $R_1$: $\delta R_1 = R_1 \widehat{\eta_1}$ | $\mathbb{R}^3$ |
| $\eta_2$ | variation of $R_2$: $\delta R_2 = R_2 \widehat{\eta_2}$ | $\mathbb{R}^3$ |
| $M_H$ | hip torque, on the swing leg, inertial frame | $\mathbb{R}^3$ |
| $M_A$ | ankle torque, on the stance leg from the ground, inertial frame | $\mathbb{R}^3$ |
| $J$ | inertia of each leg about its centre of mass, body frame | constant, symmetric |
| $l, a$ | leg length, foot to centre of mass of a leg | constants |
| $m, m_H, g, e_3$ | leg mass, hip mass, gravity, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m_H \left(l R_1 \widehat{\Omega_1} e_3 \cdot l R_1 \widehat{\Omega_1} e_3\right) \\
& + \frac{1}{2} m \left(a R_1 \widehat{\Omega_1} e_3 \cdot a R_1 \widehat{\Omega_1} e_3\right) \\
& + \frac{1}{2} m \left(\left(l R_1 \widehat{\Omega_1} e_3 -\left(l -a\right) R_2 \widehat{\Omega_2} e_3\right) \cdot \left(l R_1 \widehat{\Omega_1} e_3 -\left(l -a\right) R_2 \widehat{\Omega_2} e_3\right)\right) \\
& + \frac{1}{2} \left(\Omega_1 \cdot J \Omega_1\right) \\
& + \frac{1}{2} \left(\Omega_2 \cdot J \Omega_2\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& m_H g \left(e_3 \cdot l R_1 e_3\right) \\
& + m g \left(e_3 \cdot a R_1 e_3\right) \\
& + m g \left(e_3 \cdot \left(l R_1 e_3 -\left(l -a\right) R_2 e_3\right)\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m_H \left(l R_1 \widehat{\Omega_1} e_3 \cdot l R_1 \widehat{\Omega_1} e_3\right) \\
& + \frac{1}{2} m \left(a R_1 \widehat{\Omega_1} e_3 \cdot a R_1 \widehat{\Omega_1} e_3\right) \\
& + \frac{1}{2} m \left(\left(l R_1 \widehat{\Omega_1} e_3 -\left(l -a\right) R_2 \widehat{\Omega_2} e_3\right) \cdot \left(l R_1 \widehat{\Omega_1} e_3 -\left(l -a\right) R_2 \widehat{\Omega_2} e_3\right)\right) \\
& + \frac{1}{2} \left(\Omega_1 \cdot J \Omega_1\right) \\
& + \frac{1}{2} \left(\Omega_2 \cdot J \Omega_2\right) \\
& -m_H g \left(e_3 \cdot l R_1 e_3\right) \\
& -m g \left(e_3 \cdot a R_1 e_3\right) \\
& -m g \left(e_3 \cdot \left(l R_1 e_3 -\left(l -a\right) R_2 e_3\right)\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \eta_2 \cdot R_2^T M_H \\
& -\eta_1 \cdot R_1^T M_H \\
& + \eta_1 \cdot R_1^T M_A
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); each is written in the body frame of its leg.

<!-- model:eom -->
<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_1\)</span>, the equation of <span class="math notranslate nohighlight">\(R_1\)</span> (23 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_H l l \left(e_3 \times \left(\dot{\Omega_1} \times e_3\right)\right) \\
& -m_H l l \left(\Omega_1 \times \left(e_3 \times \left(\Omega_1 \times e_3\right)\right)\right) \\
& -m a a \left(e_3 \times \left(\dot{\Omega_1} \times e_3\right)\right) \\
& -m a a \left(\Omega_1 \times \left(e_3 \times \left(\Omega_1 \times e_3\right)\right)\right) \\
& -m l l \left(\left(\Omega_1 \times e_3\right) \times \left(R_1^T R_2 \left(\Omega_2 \times e_3\right)\right)\right) \\
& + m l a \left(\left(\Omega_1 \times e_3\right) \times \left(R_1^T R_2 \left(\Omega_2 \times e_3\right)\right)\right) \\
& -m l l \left(e_3 \times \left(\dot{\Omega_1} \times e_3\right)\right) \\
& -m l l \left(e_3 \times \left(\Omega_1 \times \left(R_1^T R_2 \left(\Omega_2 \times e_3\right)\right)\right)\right) \\
& + m l l \left(e_3 \times \left(R_1^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times e_3\right)\right)\right)\right) \\
& + m l l \left(e_3 \times \left(R_1^T R_2 \left(\dot{\Omega_2} \times e_3\right)\right)\right) \\
& + m l a \left(e_3 \times \left(\Omega_1 \times \left(R_1^T R_2 \left(\Omega_2 \times e_3\right)\right)\right)\right) \\
& -m l a \left(e_3 \times \left(R_1^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times e_3\right)\right)\right)\right) \\
& -m l a \left(e_3 \times \left(R_1^T R_2 \left(\dot{\Omega_2} \times e_3\right)\right)\right) \\
& -m l l \left(\Omega_1 \times \left(e_3 \times \left(\Omega_1 \times e_3\right)\right)\right) \\
& + m l l \left(\Omega_1 \times \left(e_3 \times \left(R_1^T R_2 \left(\Omega_2 \times e_3\right)\right)\right)\right) \\
& -m l a \left(\Omega_1 \times \left(e_3 \times \left(R_1^T R_2 \left(\Omega_2 \times e_3\right)\right)\right)\right) \\
& -J \dot{\Omega_1} \\
& -\Omega_1 \times \left(J \Omega_1\right) \\
& -l m_H g \left(e_3 \times \left(R_1^T e_3\right)\right) \\
& -a m g \left(e_3 \times \left(R_1^T e_3\right)\right) \\
& -l m g \left(e_3 \times \left(R_1^T e_3\right)\right) \\
& -R_1^T M_H \\
& + R_1^T M_A
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_2\)</span>, the equation of <span class="math notranslate nohighlight">\(R_2\)</span> (21 terms)</summary>

$$
\begin{aligned}
0 = {}& -m l l \left(\left(\Omega_2 \times e_3\right) \times \left(R_2^T R_1 \left(\Omega_1 \times e_3\right)\right)\right) \\
& + m a l \left(\left(\Omega_2 \times e_3\right) \times \left(R_2^T R_1 \left(\Omega_1 \times e_3\right)\right)\right) \\
& -m l l \left(e_3 \times \left(\Omega_2 \times \left(R_2^T R_1 \left(\Omega_1 \times e_3\right)\right)\right)\right) \\
& + m l l \left(e_3 \times \left(R_2^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times e_3\right)\right)\right)\right) \\
& + m l l \left(e_3 \times \left(R_2^T R_1 \left(\dot{\Omega_1} \times e_3\right)\right)\right) \\
& -m l l \left(e_3 \times \left(\dot{\Omega_2} \times e_3\right)\right) \\
& + 2 m l a \left(e_3 \times \left(\dot{\Omega_2} \times e_3\right)\right) \\
& + m a l \left(e_3 \times \left(\Omega_2 \times \left(R_2^T R_1 \left(\Omega_1 \times e_3\right)\right)\right)\right) \\
& -m a l \left(e_3 \times \left(R_2^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times e_3\right)\right)\right)\right) \\
& -m a l \left(e_3 \times \left(R_2^T R_1 \left(\dot{\Omega_1} \times e_3\right)\right)\right) \\
& -m a a \left(e_3 \times \left(\dot{\Omega_2} \times e_3\right)\right) \\
& + m l l \left(\Omega_2 \times \left(e_3 \times \left(R_2^T R_1 \left(\Omega_1 \times e_3\right)\right)\right)\right) \\
& -m l l \left(\Omega_2 \times \left(e_3 \times \left(\Omega_2 \times e_3\right)\right)\right) \\
& + 2 m l a \left(\Omega_2 \times \left(e_3 \times \left(\Omega_2 \times e_3\right)\right)\right) \\
& -m a l \left(\Omega_2 \times \left(e_3 \times \left(R_2^T R_1 \left(\Omega_1 \times e_3\right)\right)\right)\right) \\
& -m a a \left(\Omega_2 \times \left(e_3 \times \left(\Omega_2 \times e_3\right)\right)\right) \\
& -J \dot{\Omega_2} \\
& -\Omega_2 \times \left(J \Omega_2\right) \\
& + m g l \left(e_3 \times \left(R_2^T e_3\right)\right) \\
& -a m g \left(e_3 \times \left(R_2^T e_3\right)\right) \\
& + R_2^T M_H
\end{aligned}
$$

</details>
<!-- /model:eom -->

## Code

```python
from geomech.systems import get

# Hip xH = l R1 e3, leg centres of mass x1 = a R1 e3 and x2 = xH - (l - a) R2 e3;
# L = 1/2 mH xH'.xH' + 1/2 m x1'.x1' + 1/2 m x2'.x2' + 1/2 Om1.J Om1 + 1/2 Om2.J Om2
#     - g e3.(mH xH + m x1 + m x2), and dW = eta2.(R2^T MH) - eta1.(R1^T MH) + eta1.(R1^T MA).
system = get("biped_3d")
eom = system.derive()
```

The same model is `geomech.systems.get("biped_3d")`, with default parameters
($m_H$ = 10 kg, $m$ = 5 kg, $l$ = 1 m, $a$ = 0.5 m), a default initial state (not the
step animated above) and the 3-D scene used above; its Lagrangian is `system.L` and its virtual work `system.dW`.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 0.775 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-10}$, with $H(0) = 154$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
