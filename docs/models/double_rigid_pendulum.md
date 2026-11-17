# Double rigid pendulum

Two rigid links under gravity: link 1 hangs from a fixed pivot by a spherical joint, and
link 2 hangs from the point $R_1 l_1$ of link 1 by a second spherical joint. Each link
can swing and spin in every direction. The configuration is the pair of attitudes
$(R_1, R_2) \in SO(3) \times SO(3)$, with body angular velocities $\Omega_1, \Omega_2$ and
$\dot R_k = R_k \widehat{\Omega_k}$; no Euler angles, so no singular attitude. Not taken
from the literature; the equations are verified numerically (`check_eom`). Library
system `double_rigid_pendulum`.

## Animation

Ten seconds, unforced, from the library's initial state: link 1 tilted about 44° and
spinning at 4 rad/s about its long axis, link 2 tilted about 67° and spinning at
−5 rad/s about its own. Both links swing round the vertical in irregular loops: link 1
stays within 51° of the downward vertical and link 2 within 72°, and both keep spinning
about their long axes. Drag to orbit; the full report (equations, standard form,
parameters) opens in a
<a href="../model-reports/double_rigid_pendulum.html">separate page</a>.

<iframe class="model-report" src="../model-reports/double_rigid_pendulum.html?view" title="Double rigid pendulum, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 300" role="img" aria-labelledby="drp-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="drp-title">Double rigid pendulum: two links with attitudes R1, R2, centre-of-mass offsets rho1, rho2, joint offset l1, angular velocities Omega1, Omega2, gravity and input moments M1, M2</title>
  <defs>
    <marker id="drp-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <line x1="155.0" y1="42.0" x2="215.0" y2="42.0" stroke-width="2"/>
    <path d="M157,42.0 l-8,-10 M169,42.0 l-8,-10 M181,42.0 l-8,-10 M193,42.0 l-8,-10 M205,42.0 l-8,-10 M217,42.0 l-8,-10"/>
    <line x1="40" y1="282" x2="40" y2="244" marker-end="url(#drp-arr)"/>
    <line x1="40" y1="282" x2="78" y2="282" marker-end="url(#drp-arr)"/>
    <line x1="185.0" y1="42.0" x2="185.0" y2="272.0" stroke-dasharray="3 4" opacity="0.45"/>
    <path d="M178.9,45.5 L237.9,147.7 L250.1,140.7 L191.1,38.5 Z" stroke-width="1.8"/>
    <line x1="171.1" y1="50.0" x2="200.6" y2="101.1" stroke-width="1.6" marker-end="url(#drp-arr)"/>
    <line x1="155.6" y1="59.0" x2="214.6" y2="161.2" marker-end="url(#drp-arr)"/>
    <path d="M237.7,141.1 L186.0,247.2 L198.6,253.3 L250.3,147.3 Z" stroke-width="1.8"/>
    <line x1="229.6" y1="137.2" x2="203.8" y2="190.2" stroke-width="1.6" marker-end="url(#drp-arr)"/>
  </g>
  <path d="M219.4,125.7 A 13 5 -30 0 0 242.2,117.1" fill="none" stroke="#3f8fa0" stroke-width="1.8" marker-end="url(#drp-arr)"/>
  <path d="M200.1,40.2 Q212.5,37.7 208.9,51.3" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#drp-arr)"/>
  <path d="M193.9,219.5 A 13 5 26 0 0 213.7,233.7" fill="none" stroke="#3f8fa0" stroke-width="1.8" marker-end="url(#drp-arr)"/>
  <path d="M254.0,155.7 Q263.0,164.6 249.6,169.2" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#drp-arr)"/>
  <line x1="340" y1="40" x2="340" y2="90" stroke="#e08a2c" stroke-width="2" marker-end="url(#drp-arr)"/>
  <circle cx="185.0" cy="42.0" r="4" fill="currentColor"/>
  <circle cx="214.5" cy="93.1" r="3" fill="currentColor"/>
  <circle cx="244.0" cy="144.2" r="4" fill="currentColor"/>
  <circle cx="218.1" cy="197.2" r="3" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="143" y="81">R<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan>ρ<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan></text>
    <text x="153" y="138">R<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan>l<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan></text>
    <text x="249" y="115" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan></text>
    <text x="225" y="82">m<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan>, J<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan></text>
    <text x="215" y="36" fill="#e08a2c">M<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan></text>
    <text x="177" y="162">R<tspan font-size="12" dy="4">2</tspan><tspan dy="-4"></tspan>ρ<tspan font-size="12" dy="4">2</tspan><tspan dy="-4"></tspan></text>
    <text x="224" y="240" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">2</tspan><tspan dy="-4"></tspan></text>
    <text x="238" y="203">m<tspan font-size="12" dy="4">2</tspan><tspan dy="-4"></tspan>, J<tspan font-size="12" dy="4">2</tspan><tspan dy="-4"></tspan></text>
    <text x="267" y="167" fill="#e08a2c">M<tspan font-size="12" dy="4">2</tspan><tspan dy="-4"></tspan></text>
    <text x="348" y="80" fill="#e08a2c">g</text>
    <text x="30" y="238">e<tspan font-size="12" dy="4">3</tspan><tspan dy="-4"></tspan></text>
    <text x="84" y="287">e<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan></text>
  </g>
</svg>
<figcaption>Link 1 hangs from the fixed pivot; link 2 hangs from the joint at
<span class="math notranslate nohighlight">\(R_1 l_1\)</span>. The centres of mass (small dots) are at <span class="math notranslate nohighlight">\(R_1 \rho_1\)</span> and
<span class="math notranslate nohighlight">\(R_1 l_1 + R_2 \rho_2\)</span>. <span class="math notranslate nohighlight">\(\Omega_1\)</span>, <span class="math notranslate nohighlight">\(\Omega_2\)</span> are drawn about the long axes, where most
of the default spin is. Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span>; the inputs <span class="math notranslate nohighlight">\(M_1\)</span>, <span class="math notranslate nohighlight">\(M_2\)</span> are
body-frame moments on each link.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $R_1$ | attitude of link 1, body frame to inertial frame | $SO(3)$ |
| $R_2$ | attitude of link 2, body frame to inertial frame | $SO(3)$ |
| $\Omega_1$ | body angular velocity of link 1, $\dot R_1 = R_1 \widehat{\Omega_1}$ | $\mathbb{R}^3$ |
| $\Omega_2$ | body angular velocity of link 2, $\dot R_2 = R_2 \widehat{\Omega_2}$ | $\mathbb{R}^3$ |
| $\eta_1$ | variation of $R_1$: $\delta R_1 = R_1 \widehat{\eta_1}$ | $\mathbb{R}^3$ |
| $\eta_2$ | variation of $R_2$: $\delta R_2 = R_2 \widehat{\eta_2}$ | $\mathbb{R}^3$ |
| $M_1$ | input moment on link 1, body frame | $\mathbb{R}^3$ |
| $M_2$ | input moment on link 2, body frame | $\mathbb{R}^3$ |
| $J_1, J_2$ | inertias about the centres of mass, body frames | constant, symmetric |
| $\rho_1, \rho_2$ | upper joint to centre of mass of each link, body frame | constants |
| $l_1$ | pivot to the second joint, body frame of link 1 | constant |
| $m_1, m_2, g, e_3$ | masses, gravity, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} \left(\Omega_1 \cdot J_1 \Omega_1\right) \\
& + \frac{1}{2} \left(\Omega_2 \cdot J_2 \Omega_2\right) \\
& + \frac{1}{2} \left(R_1 \widehat{\Omega_1} \rho_1 \cdot R_1 \widehat{\Omega_1} \rho_1\right) m_1 \\
& + \frac{1}{2} \left(\left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} \rho_2\right) \cdot \left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} \rho_2\right)\right) m_2
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& \left(m_1 g \left(R_1 \rho_1 \cdot e_3\right) + m_2 g \left(\left(R_1 l_1 + R_2 \rho_2\right) \cdot e_3\right)\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} \left(\Omega_1 \cdot J_1 \Omega_1\right) \\
& + \frac{1}{2} \left(\Omega_2 \cdot J_2 \Omega_2\right) \\
& + \frac{1}{2} \left(R_1 \widehat{\Omega_1} \rho_1 \cdot R_1 \widehat{\Omega_1} \rho_1\right) m_1 \\
& + \frac{1}{2} \left(\left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} \rho_2\right) \cdot \left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} \rho_2\right)\right) m_2 \\
& -\left(m_1 g \left(R_1 \rho_1 \cdot e_3\right) + m_2 g \left(\left(R_1 l_1 + R_2 \rho_2\right) \cdot e_3\right)\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \eta_1 \cdot M_1 \\
& + \eta_2 \cdot M_2
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); each is written in the body frame of its link.

<!-- model:eom -->
Coefficient of $\eta_1$, the equation of $R_1$:

$$
\begin{aligned}
0 = {}& -J_1 \dot{\Omega_1} \\
& -\Omega_1 \times \left(J_1 \Omega_1\right) \\
& -m_1 \left(\rho_1 \times \left(\dot{\Omega_1} \times \rho_1\right)\right) \\
& -m_1 \left(\Omega_1 \times \left(\rho_1 \times \left(\Omega_1 \times \rho_1\right)\right)\right) \\
& + m_2 \left(\left(\Omega_1 \times l_1\right) \times \left(R_1^T R_2 \left(\Omega_2 \times \rho_2\right)\right)\right) \\
& -m_2 \left(l_1 \times \left(\dot{\Omega_1} \times l_1\right)\right) \\
& + m_2 \left(l_1 \times \left(\Omega_1 \times \left(R_1^T R_2 \left(\Omega_2 \times \rho_2\right)\right)\right)\right) \\
& -m_2 \left(l_1 \times \left(R_1^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times \rho_2\right)\right)\right)\right) \\
& -m_2 \left(l_1 \times \left(R_1^T R_2 \left(\dot{\Omega_2} \times \rho_2\right)\right)\right) \\
& -m_2 \left(\Omega_1 \times \left(l_1 \times \left(\Omega_1 \times l_1\right)\right)\right) \\
& -m_2 \left(\Omega_1 \times \left(l_1 \times \left(R_1^T R_2 \left(\Omega_2 \times \rho_2\right)\right)\right)\right) \\
& -m_1 g \left(\rho_1 \times \left(R_1^T e_3\right)\right) \\
& -m_2 g \left(l_1 \times \left(R_1^T e_3\right)\right) \\
& + M_1
\end{aligned}
$$

Coefficient of $\eta_2$, the equation of $R_2$:

$$
\begin{aligned}
0 = {}& -J_2 \dot{\Omega_2} \\
& -\Omega_2 \times \left(J_2 \Omega_2\right) \\
& + m_2 \left(\left(\Omega_2 \times \rho_2\right) \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right) \\
& + m_2 \left(\rho_2 \times \left(\Omega_2 \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_2 \left(\rho_2 \times \left(R_2^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_2 \left(\rho_2 \times \left(R_2^T R_1 \left(\dot{\Omega_1} \times l_1\right)\right)\right) \\
& -m_2 \left(\rho_2 \times \left(\dot{\Omega_2} \times \rho_2\right)\right) \\
& -m_2 \left(\Omega_2 \times \left(\rho_2 \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_2 \left(\Omega_2 \times \left(\rho_2 \times \left(\Omega_2 \times \rho_2\right)\right)\right) \\
& -m_2 g \left(\rho_2 \times \left(R_2^T e_3\right)\right) \\
& + M_2
\end{aligned}
$$
<!-- /model:eom -->

## Code

```python
from geomech import SO3, Dot, Matrix, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
J1 = Matrix("J1", attr=["Constant", "SymmetricMatrix"])
J2 = Matrix("J2", attr=["Constant", "SymmetricMatrix"])
rho1, rho2 = Vector("\\rho_1", attr=["Constant"]), Vector("\\rho_2", attr=["Constant"])
l1 = Vector("l_1", attr=["Constant"])
m1, m2, g = getScalars("m1 m2 g", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
M1, M2 = Vector("M1"), Vector("M2")
R1, R2 = SO3("R1"), SO3("R2")
Om1, Om2 = R1.get_tangent_vector(), R2.get_tangent_vector()
eta1, eta2 = R1.get_variation_vector(), R2.get_variation_vector()
x1, x2 = R1 * rho1, R1 * l1 + R2 * rho2  # centres of mass
v1, v2 = x1.t_diff(), x2.t_diff()

KE = (
    Dot(Om1, J1 * Om1) * half
    + Dot(Om2, J2 * Om2) * half
    + Dot(v1, v1) * m1 * half
    + Dot(v2, v2) * m2 * half
)
PE = m1 * g * Dot(x1, e3) + m2 * g * Dot(x2, e3)
L = KE - PE
dW = Dot(eta1, M1) + Dot(eta2, M2)
eom = compute_eom(L, dW, SystemVariables(matrices=[R1, R2]))
```

The same model is `geomech.systems.get("double_rigid_pendulum")`, with default
parameters (two slender links), a default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-8}$, with $H(0) = -10.2$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
