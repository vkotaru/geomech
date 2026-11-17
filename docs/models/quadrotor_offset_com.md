# Quadrotor, offset centre of mass

The [quadrotor](quadrotor.md) with its centre of mass off the thrust axis. $x$ is the centre
of the rotor plane, where the thrust acts, and the centre of mass is at $x + R\rho$, with
$\rho$ constant in the body frame. The configuration space is still
$\mathbb{R}^3 \times SO(3)$, but the offset couples translation and rotation, through the
kinetic energy and through gravity. Not taken from the literature; the equations are
verified numerically (`check_eom`). Library system `quadrotor_offset_com`.

## Animation

Two seconds, open loop, with the same start and inputs as the [quadrotor](quadrotor.md),
a start that is not the library's initial state (which is tilted and turning): level and
not rotating, at the library's initial position and velocity, with $f = m g$ and $M = 0$.
Here that does not hover. The default offset
is $\rho = (0.02, -0.01, -0.08)$ m, so the thrust through $x$ and the weight through the
centre of mass form a couple, and nothing cancels it. The vehicle tips 32° in the first
0.25 s, is more than 90° off vertical at 0.42 s and then tumbles, at up to 23 rad/s, while
it falls from 1 m to −13.6 m. This is the effect the model captures: a hover input that
ignores the offset is not an equilibrium. The run stops at 2 s because the tumble only
continues. Drag to orbit; the full report (equations, standard form, parameters) opens in
a
<a href="../model-reports/quadrotor_offset_com.html">separate page</a>.

<iframe class="model-report" src="../model-reports/quadrotor_offset_com.html?view" title="Quadrotor, offset centre of mass, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 300" role="img" aria-labelledby="qo-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="qo-title">Quadrotor with an offset centre of mass: thrust point x, offset rho to the centre of mass x_c, attitude R, angular velocity Omega_R, thrust f along R e3, moment M and gravity at x_c</title>
  <defs>
    <marker id="qo-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
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
    <line x1="340" y1="270" x2="340" y2="236" marker-end="url(#qo-arr)"/>
    <line x1="340" y1="270" x2="372" y2="270" marker-end="url(#qo-arr)"/>
  </g>
  <rect x="176" y="143" width="28" height="14" fill="currentColor"/>
  <!-- offset of the centre of mass, constant in the body frame -->
  <line x1="190" y1="150" x2="203" y2="186" stroke="currentColor" stroke-width="1.6" marker-end="url(#qo-arr)"/>
  <circle cx="204" cy="190" r="4.5" fill="currentColor"/>
  <!-- position of the thrust point -->
  <line x1="340" y1="270" x2="208" y2="164" stroke="#3f8fa0" stroke-width="2.2" marker-end="url(#qo-arr)"/>
  <!-- body angular velocity, here about the thrust axis -->
  <path d="M168,78 A22,7 0 0 0 212,78" fill="none" stroke="#3f8fa0" stroke-width="2" marker-end="url(#qo-arr)"/>
  <!-- thrust along R e3, moment, gravity -->
  <line x1="190" y1="143" x2="190" y2="40" stroke="#e08a2c" stroke-width="2.4" marker-end="url(#qo-arr)"/>
  <path d="M172,105 A18,6 0 0 0 208,105" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#qo-arr)"/>
  <line x1="204" y1="196" x2="204" y2="252" stroke="#e08a2c" stroke-width="2" marker-end="url(#qo-arr)"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="198" y="50" fill="#e08a2c">f</text>
    <text x="146" y="102" fill="#e08a2c">M</text>
    <text x="222" y="72" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">R</tspan></text>
    <text x="262" y="152" fill="#3f8fa0">R</text>
    <text x="104" y="160">m, J</text>
    <text x="262" y="240" fill="#3f8fa0">x</text>
    <text x="310" y="210">b<tspan font-size="12" dy="4">1</tspan></text>
    <text x="286" y="104">b<tspan font-size="12" dy="4">2</tspan></text>
    <text x="172" y="184">ρ</text>
    <text x="214" y="200">x<tspan font-size="12" dy="4">c</tspan></text>
    <text x="160" y="246" fill="#e08a2c">m g</text>
    <text x="346" y="242">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="356" y="290">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>As for the <a href="quadrotor.html">quadrotor</a>, but the thrust
<span class="math notranslate nohighlight">\(f R e_3\)</span> acts at <span class="math notranslate nohighlight">\(x\)</span>, the centre of the rotor plane, and the weight at the centre
of mass <span class="math notranslate nohighlight">\(x_c = x + R \rho\)</span>. With <span class="math notranslate nohighlight">\(\rho\)</span> not along <span class="math notranslate nohighlight">\(b_3\)</span> the two
forces do not act along one line.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $x$ | position of the centre of the rotor plane, where the thrust acts | $\mathbb{R}^3$ |
| $\dot x$ | velocity of $x$ | $\mathbb{R}^3$ |
| $\delta x$ | variation of $x$ | $\mathbb{R}^3$ |
| $R$ | attitude, body frame to inertial frame | $SO(3)$ |
| $\Omega_{R}$ | body angular velocity, $\dot R = R \widehat{\Omega_{R}}$ | $\mathbb{R}^3$ |
| $\eta_{R}$ | variation of $R$: $\delta R = R \widehat{\eta_{R}}$ | $\mathbb{R}^3$ |
| $f$ | total thrust, along $R e_3$, acting at $x$ | $\mathbb{R}$ |
| $M$ | moment of the rotors, body frame | $\mathbb{R}^3$ |
| $\rho$ | offset of the centre of mass from $x$, body frame | constant |
| $m, J, g, e_3$ | mass, inertia about the centre of mass (body frame), gravity, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m \left(\left(\dot{x} + R \widehat{\Omega_{R}} \rho\right) \cdot \left(\dot{x} + R \widehat{\Omega_{R}} \rho\right)\right) \\
& + \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& m g \left(e_3 \cdot \left(x + R \rho\right)\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m \left(\left(\dot{x} + R \widehat{\Omega_{R}} \rho\right) \cdot \left(\dot{x} + R \widehat{\Omega_{R}} \rho\right)\right) \\
& + \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right) \\
& -m g \left(e_3 \cdot \left(x + R \rho\right)\right)
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
& -m R \left(\Omega_{R} \times \left(\Omega_{R} \times \rho\right)\right) \\
& -m R \left(\dot{\Omega_{R}} \times \rho\right) \\
& -m g e_3 \\
& + f R e_3
\end{aligned}
$$

Coefficient of $\eta_{R}$, the equation of $R$:

$$
\begin{aligned}
0 = {}& m \left(\left(\Omega_{R} \times \rho\right) \times \left(R^T \dot{x}\right)\right) \\
& + m \left(\rho \times \left(\Omega_{R} \times \left(R^T \dot{x}\right)\right)\right) \\
& -m \left(\rho \times \left(R^T \ddot{x}\right)\right) \\
& -m \left(\rho \times \left(\dot{\Omega_{R}} \times \rho\right)\right) \\
& -m \left(\Omega_{R} \times \left(\rho \times \left(R^T \dot{x}\right)\right)\right) \\
& -m \left(\Omega_{R} \times \left(\rho \times \left(\Omega_{R} \times \rho\right)\right)\right) \\
& -J \dot{\Omega_{R}} \\
& -\Omega_{R} \times \left(J \Omega_{R}\right) \\
& -m g \left(\rho \times \left(R^T e_3\right)\right) \\
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
rho = Vector("\\rho", attr=["Constant"])
x, R = Vector("x"), SO3("R")
Om, eta = R.get_tangent_vector(), R.get_variation_vector()
f, M = Scalar("f"), Vector("M")
xc = x + R * rho  # centre of mass
vc = xc.t_diff()

L = half * m * Dot(vc, vc) + half * Dot(Om, J * Om) - m * g * Dot(e3, xc)
dW = Dot(x.delta(), (R * e3) * f) + Dot(eta, M)
eom = compute_eom(L, dW, SystemVariables(vectors=[x], matrices=[R]))
```

The same model is `geomech.systems.get("quadrotor_offset_com")`, with default parameters,
a default initial state (tilted 13° and turning, not the level start animated above) and
the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-7}$ |
| Energy $H = T + V$, 2 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-9}$, with $H(0) = 18.4$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
