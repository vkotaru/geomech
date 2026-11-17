# Free rigid body

A rigid body turning freely about its centre of mass, driven by a moment $M$ in the body
frame. The configuration is the attitude $R \in SO(3)$, which maps body axes to inertial
axes, with body angular velocity $\Omega_{R}$ and $\dot R = R \widehat{\Omega_{R}}$. No
Euler angles are used, so the equations have no singular attitude. Lee et al. §6.6.1,
eq. 6.16. Library system `free_rigid_body`.

## Animation

Ten seconds, unforced, from the library's initial state, $\Omega_{R} = (0.3, -1.2, 0.8)$
rad/s. The body spins mostly about $b_2$, its axis of largest inertia: the $b_2$ component of
$\Omega_{R}$ stays between $-1.25$ and $-1.13$ rad/s, while the other two components
oscillate with a period of about 10 s. In the 10 s it turns through about 2.3 revolutions, and the inertial
angular momentum $R J \Omega_{R}$ stays constant. Drag to orbit; the full report
(equations, standard form, parameters) opens in a
<a href="../model-reports/free_rigid_body.html">separate page</a>.

<iframe class="model-report" src="../model-reports/free_rigid_body.html?view" title="Free rigid body, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 260" role="img" aria-labelledby="frb-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="frb-title">Free rigid body: a box with body axes, attitude R, angular velocity Omega_R, inertia J and input moment M</title>
  <defs>
    <marker id="frb-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- box, oblique view -->
    <path d="M130,110 L250,110 L250,170 L130,170 Z" stroke-width="1.8"/>
    <path d="M130,110 L170,80 L290,80 L250,110 M250,170 L290,140 L290,80" stroke-width="1.8"/>
    <!-- body axes from the centre of mass -->
    <line x1="210" y1="125" x2="340" y2="125" marker-end="url(#frb-arr)"/>
    <line x1="210" y1="125" x2="282" y2="71" marker-end="url(#frb-arr)"/>
    <line x1="210" y1="125" x2="210" y2="40" marker-end="url(#frb-arr)"/>
    <!-- inertial axes -->
    <line x1="40" y1="230" x2="40" y2="190" marker-end="url(#frb-arr)"/>
    <line x1="40" y1="230" x2="80" y2="230" marker-end="url(#frb-arr)"/>
  </g>
  <!-- attitude R: inertial frame to body frame -->
  <path d="M 52,222 Q 70,170 122,168" fill="none" stroke="#3f8fa0" stroke-width="1.6" stroke-dasharray="4 3" marker-end="url(#frb-arr)"/>
  <!-- angular velocity -->
  <line x1="210" y1="125" x2="140" y2="50" stroke="#3f8fa0" stroke-width="3" marker-end="url(#frb-arr)"/>
  <!-- input moment about b3 -->
  <path d="M 192,62 A 22 8 0 0 0 230,60" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#frb-arr)"/>
  <circle cx="210" cy="125" r="3.5" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="90" y="202" fill="#3f8fa0">R</text>
    <text x="106" y="46" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">R</tspan></text>
    <text x="238" y="64" fill="#e08a2c">M</text>
    <text x="180" y="152">J</text>
    <text x="346" y="130">b<tspan font-size="12" dy="4">1</tspan></text>
    <text x="288" y="68">b<tspan font-size="12" dy="4">2</tspan></text>
    <text x="216" y="38">b<tspan font-size="12" dy="4">3</tspan></text>
    <text x="34" y="183">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="85" y="235">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>The body axes <span class="math notranslate nohighlight">\(b_1\)</span>, <span class="math notranslate nohighlight">\(b_2\)</span>, <span class="math notranslate nohighlight">\(b_3\)</span> are fixed in the body at its centre of mass; the
attitude <span class="math notranslate nohighlight">\(R\)</span> maps them to the inertial axes (<span class="math notranslate nohighlight">\(e_1\)</span> and <span class="math notranslate nohighlight">\(e_3\)</span> drawn). The angular velocity
<span class="math notranslate nohighlight">\(\Omega_R\)</span> and the moment <span class="math notranslate nohighlight">\(M\)</span> are body-frame vectors. There is no
gravity and no translation.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $R$ | attitude, body frame to inertial frame | $SO(3)$ |
| $\Omega_{R}$ | body angular velocity, $\dot R = R \widehat{\Omega_{R}}$ | $\mathbb{R}^3$ |
| $\eta_{R}$ | variation of $R$: $\delta R = R \widehat{\eta_{R}}$ | $\mathbb{R}^3$ |
| $M$ | input moment, body frame | $\mathbb{R}^3$ |
| $J$ | inertia about the centre of mass, body frame | constant, symmetric |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& 0
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \eta_{R} \cdot M
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); it is written in the body frame.

<!-- model:eom -->
Coefficient of $\eta_{R}$, the equation of $R$:

$$
\begin{aligned}
0 = {}& -J \dot{\Omega_{R}} \\
& -\Omega_{R} \times \left(J \Omega_{R}\right) \\
& + M
\end{aligned}
$$
<!-- /model:eom -->

Multiplied by $-1$, this is Lee et al. eq. 6.16 with the input added:
$J \dot\Omega_R + \Omega_R \times J\Omega_R = M$.

## Code

```python
from geomech import SO3, Dot, Matrix, Scalar, SystemVariables, Vector, compute_eom

half = Scalar("0.5", value=0.5, attr=["Constant"])
J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
M = Vector("M")
R = SO3("R")
Om, eta = R.get_tangent_vector(), R.get_variation_vector()

L = Dot(Om, J * Om) * half
dW = Dot(eta, M)
eom = compute_eom(L, dW, SystemVariables(matrices=[R]))
```

The same model is `geomech.systems.get("free_rigid_body")`, with default parameters, a
default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-7}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-12}$, with $H(0) = 0.666$ |
| Closed form, Lee et al. eq. 6.16 | equal numerically at random states (`tests/test_systems_book.py`) |
<!-- /model:checks -->
