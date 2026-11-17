# Two rigid bodies, spherical joint

Two rigid bodies flying freely under gravity, joined by a spherical joint. Body 1 has
its centre of mass at $x$ and the joint at $x + R_1 \rho_1$; body 2 has its centre of
mass at $x + R_1 \rho_1 + R_2 \rho_2$. The configuration is
$(x, R_1, R_2) \in \mathbb{R}^3 \times SO(3) \times SO(3)$, with body angular velocities
$\Omega_1, \Omega_2$ and $\dot R_k = R_k \widehat{\Omega_k}$; no Euler angles, so no
singular attitude. A force $F$ acts at the centre of mass of body 1 and a body-frame
moment on each body. Not taken from the literature; the equations are verified
numerically (`check_eom`). Library system `two_rigid_bodies_spherical_joint`.

## Animation

One second, unforced, from a start that is not the library's initial state: body 1
starts 1 m up, thrown upward at 4.9 m/s with a small sideways velocity, with both bodies
turning at three times the library's initial rates (the library's initial state throws
body 1 upward at 2 m/s); they tumble at up to 8.4 rad/s. The pair is in free
fall: $x$ peaks 2.37 m up at about 0.55 s and is back at 1.29 m after 1 s. Body 1
turns through about 1.1 revolutions and body 2 about 0.8, tumbling about the joint
relative to each other. Drag to orbit; the full report (equations, standard form,
parameters) opens in a
<a href="../model-reports/two_rigid_bodies_spherical_joint.html">separate page</a>.

<iframe class="model-report" src="../model-reports/two_rigid_bodies_spherical_joint.html?view" title="Two rigid bodies, spherical joint, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 280" role="img" aria-labelledby="trb-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="trb-title">Two rigid bodies joined by a spherical joint: position x, offsets rho1, rho2, attitudes, angular velocities Omega1, Omega2, force F, moments M1, M2 and gravity</title>
  <defs>
    <marker id="trb-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- inertial axes -->
    <line x1="50" y1="240" x2="50" y2="200" marker-end="url(#trb-arr)"/>
    <line x1="50" y1="240" x2="90" y2="240" marker-end="url(#trb-arr)"/>
    <!-- body 1 and body 2 -->
    <path d="M72.1,159.7 L222.5,105.3 L227.9,120.3 L77.5,174.7 Z" stroke-width="1.8"/>
    <path d="M230.8,108.6 L310.6,214.8 L299.4,223.2 L219.6,117 Z" stroke-width="1.8"/>
    <!-- R1 rho1: centre of mass of body 1 to the joint; R2 rho2: joint to centre of mass of body 2 -->
    <line x1="150" y1="140" x2="221" y2="114.3" stroke-width="2" marker-end="url(#trb-arr)"/>
    <line x1="225.2" y1="112.8" x2="263.4" y2="163.8" stroke-width="2" marker-end="url(#trb-arr)"/>
  </g>
  <!-- position and velocity of body 1, angular velocities -->
  <line x1="50" y1="240" x2="148" y2="142" stroke="#3f8fa0" stroke-width="2.4" marker-end="url(#trb-arr)"/>
  <line x1="150" y1="140" x2="150" y2="60" stroke="#3f8fa0" stroke-width="2" marker-end="url(#trb-arr)"/>
  <line x1="150" y1="140" x2="110" y2="70" stroke="#3f8fa0" stroke-width="3" marker-end="url(#trb-arr)"/>
  <line x1="265.2" y1="166.2" x2="335" y2="140" stroke="#3f8fa0" stroke-width="3" marker-end="url(#trb-arr)"/>
  <!-- force at the centre of mass of body 1, moments, gravity -->
  <line x1="150" y1="140" x2="185" y2="200" stroke="#e08a2c" stroke-width="2" marker-end="url(#trb-arr)"/>
  <path d="M 66,148 A 16 16 0 0 1 92,140" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#trb-arr)"/>
  <path d="M 322,206 A 18 18 0 0 1 320,232" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#trb-arr)"/>
  <line x1="360" y1="40" x2="360" y2="90" stroke="#e08a2c" stroke-width="2" marker-end="url(#trb-arr)"/>
  <circle cx="50" cy="240" r="3" fill="currentColor"/>
  <circle cx="150" cy="140" r="3.5" fill="currentColor"/>
  <circle cx="225.2" cy="112.8" r="4.5" fill="currentColor"/>
  <circle cx="265.2" cy="166.2" r="3.5" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="80" y="186" fill="#3f8fa0">x</text>
    <text x="156" y="62" fill="#3f8fa0">ẋ</text>
    <text x="84" y="66" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">1</tspan></text>
    <text x="326" y="126" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">2</tspan></text>
    <text x="172" y="98">R<tspan font-size="12" dy="4">1</tspan><tspan dy="-4">ρ</tspan><tspan font-size="12" dy="4">1</tspan></text>
    <text x="200" y="160">R<tspan font-size="12" dy="4">2</tspan><tspan dy="-4">ρ</tspan><tspan font-size="12" dy="4">2</tspan></text>
    <text x="112" y="196">m<tspan font-size="12" dy="4">1</tspan><tspan dy="-4">, J</tspan><tspan font-size="12" dy="4">1</tspan></text>
    <text x="226" y="214">m<tspan font-size="12" dy="4">2</tspan><tspan dy="-4">, J</tspan><tspan font-size="12" dy="4">2</tspan></text>
    <text x="190" y="214" fill="#e08a2c">F</text>
    <text x="44" y="124" fill="#e08a2c">M<tspan font-size="12" dy="4">1</tspan></text>
    <text x="332" y="240" fill="#e08a2c">M<tspan font-size="12" dy="4">2</tspan></text>
    <text x="366" y="80" fill="#e08a2c">g</text>
    <text x="44" y="193">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="95" y="245">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>Body 1 (left) has its centre of mass at <span class="math notranslate nohighlight">\(x\)</span>; the joint (large dot) is at
<span class="math notranslate nohighlight">\(x + R_1 \rho_1\)</span>, and the centre of mass of body 2 at
<span class="math notranslate nohighlight">\(x + R_1 \rho_1 + R_2 \rho_2\)</span>. <span class="math notranslate nohighlight">\(\dot x\)</span> and <span class="math notranslate nohighlight">\(F\)</span> are
inertial-frame vectors at the centre of mass of body 1; <span class="math notranslate nohighlight">\(\Omega_1\)</span>, <span class="math notranslate nohighlight">\(\Omega_2\)</span> and the moments
<span class="math notranslate nohighlight">\(M_1\)</span>, <span class="math notranslate nohighlight">\(M_2\)</span> are body-frame vectors. Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span>.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $x$ | position of the centre of mass of body 1, inertial frame | $\mathbb{R}^3$ |
| $\dot x$ | its velocity, inertial frame | $\mathbb{R}^3$ |
| $\delta x$ | variation of $x$ | $\mathbb{R}^3$ |
| $R_1$ | attitude of body 1, body frame to inertial frame | $SO(3)$ |
| $R_2$ | attitude of body 2, body frame to inertial frame | $SO(3)$ |
| $\Omega_1$ | body angular velocity of body 1, $\dot R_1 = R_1 \widehat{\Omega_1}$ | $\mathbb{R}^3$ |
| $\Omega_2$ | body angular velocity of body 2, $\dot R_2 = R_2 \widehat{\Omega_2}$ | $\mathbb{R}^3$ |
| $\eta_1$ | variation of $R_1$: $\delta R_1 = R_1 \widehat{\eta_1}$ | $\mathbb{R}^3$ |
| $\eta_2$ | variation of $R_2$: $\delta R_2 = R_2 \widehat{\eta_2}$ | $\mathbb{R}^3$ |
| $F$ | input force at the centre of mass of body 1, inertial frame | $\mathbb{R}^3$ |
| $M_1$ | input moment on body 1, body frame | $\mathbb{R}^3$ |
| $M_2$ | input moment on body 2, body frame | $\mathbb{R}^3$ |
| $J_1, J_2$ | inertias about the centres of mass, body frames | constant, symmetric |
| $\rho_1$ | centre of mass of body 1 to the joint, body frame of body 1 | constant |
| $\rho_2$ | joint to the centre of mass of body 2, body frame of body 2 | constant |
| $m_1, m_2, g, e_3$ | masses, gravity, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m_1 \left(\dot{x} \cdot \dot{x}\right) \\
& + \frac{1}{2} m_2 \left(\left(\dot{x} + R_1 \widehat{\Omega_1} \rho_1 + R_2 \widehat{\Omega_2} \rho_2\right) \cdot \left(\dot{x} + R_1 \widehat{\Omega_1} \rho_1 + R_2 \widehat{\Omega_2} \rho_2\right)\right) \\
& + \frac{1}{2} \left(\Omega_1 \cdot J_1 \Omega_1\right) \\
& + \frac{1}{2} \left(\Omega_2 \cdot J_2 \Omega_2\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& m_1 g \left(e_3 \cdot x\right) \\
& + m_2 g \left(e_3 \cdot \left(x + R_1 \rho_1 + R_2 \rho_2\right)\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m_1 \left(\dot{x} \cdot \dot{x}\right) \\
& + \frac{1}{2} m_2 \left(\left(\dot{x} + R_1 \widehat{\Omega_1} \rho_1 + R_2 \widehat{\Omega_2} \rho_2\right) \cdot \left(\dot{x} + R_1 \widehat{\Omega_1} \rho_1 + R_2 \widehat{\Omega_2} \rho_2\right)\right) \\
& + \frac{1}{2} \left(\Omega_1 \cdot J_1 \Omega_1\right) \\
& + \frac{1}{2} \left(\Omega_2 \cdot J_2 \Omega_2\right) \\
& -m_1 g \left(e_3 \cdot x\right) \\
& -m_2 g \left(e_3 \cdot \left(x + R_1 \rho_1 + R_2 \rho_2\right)\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \delta x \cdot F \\
& + \eta_1 \cdot M_1 \\
& + \eta_2 \cdot M_2
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); the equation of $x$ is written in the inertial
frame and the equations of $R_1$ and $R_2$ in the body frames.

<!-- model:eom -->
Coefficient of $\delta x$, the equation of $x$:

$$
\begin{aligned}
0 = {}& -m_1 \ddot{x} \\
& -m_2 \ddot{x} \\
& -m_2 R_1 \left(\Omega_1 \times \left(\Omega_1 \times \rho_1\right)\right) \\
& -m_2 R_1 \left(\dot{\Omega_1} \times \rho_1\right) \\
& -m_2 R_2 \left(\Omega_2 \times \left(\Omega_2 \times \rho_2\right)\right) \\
& -m_2 R_2 \left(\dot{\Omega_2} \times \rho_2\right) \\
& -m_1 g e_3 \\
& -m_2 g e_3 \\
& + F
\end{aligned}
$$

Coefficient of $\eta_1$, the equation of $R_1$:

$$
\begin{aligned}
0 = {}& m_2 \left(\left(\Omega_1 \times \rho_1\right) \times \left(R_1^T \dot{x}\right)\right) \\
& + m_2 \left(\left(\Omega_1 \times \rho_1\right) \times \left(R_1^T R_2 \left(\Omega_2 \times \rho_2\right)\right)\right) \\
& + m_2 \left(\rho_1 \times \left(\Omega_1 \times \left(R_1^T \dot{x}\right)\right)\right) \\
& -m_2 \left(\rho_1 \times \left(R_1^T \ddot{x}\right)\right) \\
& -m_2 \left(\rho_1 \times \left(\dot{\Omega_1} \times \rho_1\right)\right) \\
& + m_2 \left(\rho_1 \times \left(\Omega_1 \times \left(R_1^T R_2 \left(\Omega_2 \times \rho_2\right)\right)\right)\right) \\
& -m_2 \left(\rho_1 \times \left(R_1^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times \rho_2\right)\right)\right)\right) \\
& -m_2 \left(\rho_1 \times \left(R_1^T R_2 \left(\dot{\Omega_2} \times \rho_2\right)\right)\right) \\
& -m_2 \left(\Omega_1 \times \left(\rho_1 \times \left(R_1^T \dot{x}\right)\right)\right) \\
& -m_2 \left(\Omega_1 \times \left(\rho_1 \times \left(\Omega_1 \times \rho_1\right)\right)\right) \\
& -m_2 \left(\Omega_1 \times \left(\rho_1 \times \left(R_1^T R_2 \left(\Omega_2 \times \rho_2\right)\right)\right)\right) \\
& -J_1 \dot{\Omega_1} \\
& -\Omega_1 \times \left(J_1 \Omega_1\right) \\
& -m_2 g \left(\rho_1 \times \left(R_1^T e_3\right)\right) \\
& + M_1
\end{aligned}
$$

Coefficient of $\eta_2$, the equation of $R_2$:

$$
\begin{aligned}
0 = {}& m_2 \left(\left(\Omega_2 \times \rho_2\right) \times \left(R_2^T \dot{x}\right)\right) \\
& + m_2 \left(\left(\Omega_2 \times \rho_2\right) \times \left(R_2^T R_1 \left(\Omega_1 \times \rho_1\right)\right)\right) \\
& + m_2 \left(\rho_2 \times \left(\Omega_2 \times \left(R_2^T \dot{x}\right)\right)\right) \\
& -m_2 \left(\rho_2 \times \left(R_2^T \ddot{x}\right)\right) \\
& + m_2 \left(\rho_2 \times \left(\Omega_2 \times \left(R_2^T R_1 \left(\Omega_1 \times \rho_1\right)\right)\right)\right) \\
& -m_2 \left(\rho_2 \times \left(R_2^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times \rho_1\right)\right)\right)\right) \\
& -m_2 \left(\rho_2 \times \left(R_2^T R_1 \left(\dot{\Omega_1} \times \rho_1\right)\right)\right) \\
& -m_2 \left(\rho_2 \times \left(\dot{\Omega_2} \times \rho_2\right)\right) \\
& -m_2 \left(\Omega_2 \times \left(\rho_2 \times \left(R_2^T \dot{x}\right)\right)\right) \\
& -m_2 \left(\Omega_2 \times \left(\rho_2 \times \left(R_2^T R_1 \left(\Omega_1 \times \rho_1\right)\right)\right)\right) \\
& -m_2 \left(\Omega_2 \times \left(\rho_2 \times \left(\Omega_2 \times \rho_2\right)\right)\right) \\
& -J_2 \dot{\Omega_2} \\
& -\Omega_2 \times \left(J_2 \Omega_2\right) \\
& -m_2 g \left(\rho_2 \times \left(R_2^T e_3\right)\right) \\
& + M_2
\end{aligned}
$$
<!-- /model:eom -->

## Code

```python
from geomech import SO3, Dot, Matrix, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
m1, m2, g = getScalars("m1 m2 g", attr=["Constant"])
J1 = Matrix("J1", attr=["Constant", "SymmetricMatrix"])
J2 = Matrix("J2", attr=["Constant", "SymmetricMatrix"])
rho1, rho2 = Vector("\\rho_1", attr=["Constant"]), Vector("\\rho_2", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
x, R1, R2 = Vector("x"), SO3("R1"), SO3("R2")
F, M1, M2 = Vector("F"), Vector("M1"), Vector("M2")
Om1, Om2 = R1.get_tangent_vector(), R2.get_tangent_vector()
eta1, eta2 = R1.get_variation_vector(), R2.get_variation_vector()
x2 = x + R1 * rho1 + R2 * rho2  # centre of mass of body 2
v1, v2 = x.t_diff(), x2.t_diff()

L = (
    half * m1 * Dot(v1, v1)
    + half * m2 * Dot(v2, v2)
    + half * Dot(Om1, J1 * Om1)
    + half * Dot(Om2, J2 * Om2)
    - m1 * g * Dot(e3, x)
    - m2 * g * Dot(e3, x2)
)
dW = Dot(x.delta(), F) + Dot(eta1, M1) + Dot(eta2, M2)
eom = compute_eom(L, dW, SystemVariables(vectors=[x], matrices=[R1, R2]))
```

The same model is `geomech.systems.get("two_rigid_bodies_spherical_joint")`, with
default parameters, a default initial state (not the toss animated above) and the 3-D
scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 1 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-7}$, with $H(0) = 39.1$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
