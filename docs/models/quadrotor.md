# Quadrotor

A quadrotor as a rigid body: position $x \in \mathbb{R}^3$ of its centre of mass and
attitude $R \in SO(3)$, so the configuration space is $\mathbb{R}^3 \times SO(3)$. The
rotors give a total thrust $f$ along the body axis $R e_3$, through the centre of mass, and
a moment $M$ in the body frame. Model from the quadrotor paper of Lee, Leok and McClamroch
(2010), with $e_3$ pointing up here (the paper's points down); the equations are verified
numerically (`check_eom`). Library system `quadrotor`.

## Animation

Six seconds, open loop: there is no controller. The start is not the library's initial
state (which is tilted and turning): the vehicle starts level and not rotating, at the
library's initial position and velocity, $x = (0, 0.5, 1)$ m and
$\dot x = (0.3, -0.2, 0.5)$ m/s, with constant thrust $f = m g = 19.6$ N and $M = 0$.
Thrust and weight cancel, so the attitude stays exactly level and the vehicle coasts at
its initial velocity: 2.2 m horizontally, and up from 1 m to 4 m. Drag to orbit; the full
report (equations, standard form, parameters) opens in a
<a href="../model-reports/quadrotor.html">separate page</a>.

<iframe class="model-report" src="../model-reports/quadrotor.html?view" title="Quadrotor, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 300" role="img" aria-labelledby="qr-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="qr-title">Quadrotor: arms along b1 and b2 with four rotors, position x, attitude R, angular velocity Omega_R, thrust f along R e3, moment M and gravity</title>
  <defs>
    <marker id="qr-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- arms along b1 and b2, rotors at their ends (the rotor on +b1 is marked) -->
    <line x1="95" y1="122" x2="285" y2="178" stroke-width="2.4"/>
    <line x1="130" y1="182" x2="250" y2="118" stroke-width="2.4"/>
    <ellipse cx="285" cy="178" rx="32" ry="10" stroke-width="1.6" fill="currentColor" fill-opacity="0.15"/>
    <ellipse cx="250" cy="118" rx="32" ry="10" stroke-width="1.6"/>
    <ellipse cx="95" cy="122" rx="32" ry="10" stroke-width="1.6"/>
    <ellipse cx="130" cy="182" rx="32" ry="10" stroke-width="1.6"/>
    <!-- inertial axes -->
    <line x1="340" y1="270" x2="340" y2="236" marker-end="url(#qr-arr)"/>
    <line x1="340" y1="270" x2="372" y2="270" marker-end="url(#qr-arr)"/>
  </g>
  <rect x="176" y="143" width="28" height="14" fill="currentColor"/>
  <!-- position of the thrust point -->
  <line x1="340" y1="270" x2="205" y2="162" stroke="#3f8fa0" stroke-width="2.2" marker-end="url(#qr-arr)"/>
  <!-- body angular velocity, here about the thrust axis -->
  <path d="M168,78 A22,7 0 0 0 212,78" fill="none" stroke="#3f8fa0" stroke-width="2" marker-end="url(#qr-arr)"/>
  <!-- thrust along R e3, moment, gravity -->
  <line x1="190" y1="143" x2="190" y2="40" stroke="#e08a2c" stroke-width="2.4" marker-end="url(#qr-arr)"/>
  <path d="M172,105 A18,6 0 0 0 208,105" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#qr-arr)"/>
  <line x1="190" y1="160" x2="190" y2="225" stroke="#e08a2c" stroke-width="2" marker-end="url(#qr-arr)"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="198" y="50" fill="#e08a2c">f</text>
    <text x="146" y="102" fill="#e08a2c">M</text>
    <text x="222" y="72" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">R</tspan></text>
    <text x="262" y="152" fill="#3f8fa0">R</text>
    <text x="104" y="160">m, J</text>
    <text x="262" y="240" fill="#3f8fa0">x</text>
    <text x="310" y="210">b<tspan font-size="12" dy="4">1</tspan></text>
    <text x="286" y="104">b<tspan font-size="12" dy="4">2</tspan></text>
    <text x="148" y="222" fill="#e08a2c">m g</text>
    <text x="346" y="242">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="356" y="290">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>The body axes <span class="math notranslate nohighlight">\(b_1\)</span>, <span class="math notranslate nohighlight">\(b_2\)</span> run along the arms and <span class="math notranslate nohighlight">\(b_3 = R e_3\)</span> is
the thrust axis; the rotor on <span class="math notranslate nohighlight">\(+b_1\)</span> is shaded. <span class="math notranslate nohighlight">\(x\)</span> is the centre of mass, where the
thrust <span class="math notranslate nohighlight">\(f R e_3\)</span> and the weight act. <span class="math notranslate nohighlight">\(\Omega_R\)</span> and <span class="math notranslate nohighlight">\(M\)</span> are drawn about <span class="math notranslate nohighlight">\(b_3\)</span> for
illustration; in general both point anywhere. Gravity acts along
<span class="math notranslate nohighlight">\(-e_3\)</span>.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $x$ | position of the centre of mass | $\mathbb{R}^3$ |
| $\dot x$ | velocity of the centre of mass | $\mathbb{R}^3$ |
| $\delta x$ | variation of $x$ | $\mathbb{R}^3$ |
| $R$ | attitude, body frame to inertial frame | $SO(3)$ |
| $\Omega_{R}$ | body angular velocity, $\dot R = R \widehat{\Omega_{R}}$ | $\mathbb{R}^3$ |
| $\eta_{R}$ | variation of $R$: $\delta R = R \widehat{\eta_{R}}$ | $\mathbb{R}^3$ |
| $f$ | total thrust, along $R e_3$ | $\mathbb{R}$ |
| $M$ | moment of the rotors, body frame | $\mathbb{R}^3$ |
| $m, J, g, e_3$ | mass, inertia about the centre of mass (body frame), gravity, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m \left(\dot{x} \cdot \dot{x}\right) \\
& + \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& m g \left(e_3 \cdot x\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m \left(\dot{x} \cdot \dot{x}\right) \\
& + \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right) \\
& -m g \left(e_3 \cdot x\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \delta x \cdot f R e_3 \\
& + \eta_{R} \cdot M
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); the equation of $x$ is in the inertial frame and the
equation of $R$ in the body frame.

<!-- model:eom -->
Coefficient of $\delta x$, the equation of $x$:

$$
\begin{aligned}
0 = {}& -m \ddot{x} \\
& -m g e_3 \\
& + f R e_3
\end{aligned}
$$

Coefficient of $\eta_{R}$, the equation of $R$:

$$
\begin{aligned}
0 = {}& -J \dot{\Omega_{R}} \\
& -\Omega_{R} \times \left(J \Omega_{R}\right) \\
& + M
\end{aligned}
$$
<!-- /model:eom -->

## Code

```python
from geomech import SO3, Dot, Matrix, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
m, g = getScalars("m g", attr=["Constant"])
J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
e3 = Vector("e3", attr=["Constant"])
x, R = Vector("x"), SO3("R")
Om, eta = R.get_tangent_vector(), R.get_variation_vector()
f, M = Scalar("f"), Vector("M")
v = x.t_diff()

L = half * m * Dot(v, v) + half * Dot(Om, J * Om) - m * g * Dot(e3, x)
dW = Dot(x.delta(), (R * e3) * f) + Dot(eta, M)
eom = compute_eom(L, dW, SystemVariables(vectors=[x], matrices=[R]))
```

The same model is `geomech.systems.get("quadrotor")`, with default parameters
($m = 2$ kg), a default initial state (tilted 13° and turning at
$\Omega_R = (0.3, -0.2, 0.5)$ rad/s, not the level start animated above) and the 3-D scene
used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 6 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-8}$, with $H(0) = 20$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
