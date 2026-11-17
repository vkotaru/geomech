# Rigid pendulum

A rigid body hung from a fixed pivot by a spherical joint, free to swing and spin in every
direction under gravity; its centre of mass is at $R \rho$ from the pivot. The
configuration is the attitude $R \in SO(3)$, with body angular velocity $\Omega_{R}$ and
$\dot R = R \widehat{\Omega_{R}}$. No Euler angles are used, so the equations have no
singular attitude. Lee et al. §6.6.2, eq. 6.20. Library system `rigid_pendulum`.

## Animation

Ten seconds, unforced, from the library's initial state: the link tilted about 44° from
the hanging position, spinning at 4 rad/s about its long axis with a small swinging
rate. The link swings between 10° and 44° from the downward vertical and goes round the
vertical about six times, in a rosette, while it keeps spinning about its long axis.
Drag to orbit; the full report (equations, standard form, parameters) opens in a
<a href="../model-reports/rigid_pendulum.html">separate page</a>.

<iframe class="model-report" src="../model-reports/rigid_pendulum.html?view" title="Rigid pendulum, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 300" role="img" aria-labelledby="rp-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="rp-title">Rigid pendulum: pivot O, centre of mass at R rho, mass m, inertia J, angular velocity Omega_R, gravity and input moment M</title>
  <defs>
    <marker id="rp-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <line x1="155.0" y1="42.0" x2="215.0" y2="42.0" stroke-width="2"/>
    <path d="M157,42.0 l-8,-10 M169,42.0 l-8,-10 M181,42.0 l-8,-10 M193,42.0 l-8,-10 M205,42.0 l-8,-10 M217,42.0 l-8,-10"/>
    <line x1="40" y1="282" x2="40" y2="244" marker-end="url(#rp-arr)"/>
    <line x1="40" y1="282" x2="78" y2="282" marker-end="url(#rp-arr)"/>
    <line x1="185.0" y1="42.0" x2="185.0" y2="272.0" stroke-dasharray="3 4" opacity="0.45"/>
    <path d="M178.5,44.6 L242.2,202.2 L255.2,197.0 L191.5,39.4 Z" stroke-width="1.8"/>
    <line x1="170.2" y1="48.0" x2="202.0" y2="126.8" stroke-width="1.6" marker-end="url(#rp-arr)"/>
  </g>
  <path d="M222.3,166.3 A 13 5 -22 0 0 246.0,161.0" fill="none" stroke="#3f8fa0" stroke-width="1.8" marker-end="url(#rp-arr)"/>
  <path d="M200.2,42.3 Q212.9,41.5 207.3,54.5" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#rp-arr)"/>
  <line x1="340" y1="40" x2="340" y2="90" stroke="#e08a2c" stroke-width="2" marker-end="url(#rp-arr)"/>
  <circle cx="185.0" cy="42.0" r="4" fill="currentColor"/>
  <circle cx="216.8" cy="120.8" r="3" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="143" y="91">R ρ</text>
    <text x="254" y="160" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">R</tspan><tspan dy="-4"></tspan></text>
    <text x="228" y="107">m, J</text>
    <text x="216" y="40" fill="#e08a2c">M</text>
    <text x="150" y="62">O</text>
    <text x="348" y="80" fill="#e08a2c">g</text>
    <text x="30" y="238">e<tspan font-size="12" dy="4">3</tspan><tspan dy="-4"></tspan></text>
    <text x="84" y="287">e<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan></text>
  </g>
</svg>
<figcaption>The pivot O is fixed. The centre of mass is at <span class="math notranslate nohighlight">\(R \rho\)</span> from the pivot,
with <span class="math notranslate nohighlight">\(\rho\)</span> a constant body-frame vector (by default along body <span class="math notranslate nohighlight">\(-e_3\)</span>);
<span class="math notranslate nohighlight">\(J\)</span> is the inertia about the pivot. <span class="math notranslate nohighlight">\(\Omega_R\)</span> is drawn about the long
axis, where most of the default spin is. Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span>; the input <span class="math notranslate nohighlight">\(M\)</span> is a
body-frame moment.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $R$ | attitude, body frame to inertial frame | $SO(3)$ |
| $\Omega_{R}$ | body angular velocity, $\dot R = R \widehat{\Omega_{R}}$ | $\mathbb{R}^3$ |
| $\eta_{R}$ | variation of $R$: $\delta R = R \widehat{\eta_{R}}$ | $\mathbb{R}^3$ |
| $M$ | input moment, body frame | $\mathbb{R}^3$ |
| $J$ | inertia about the pivot, body frame | constant, symmetric |
| $\rho$ | pivot to centre of mass, body frame | constant |
| $m, g, e_3$ | mass, gravity, vertical unit vector | constants |

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
V = {}& m g \left(R \rho \cdot e_3\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right) \\
& -m g \left(R \rho \cdot e_3\right)
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
& -m g \left(\rho \times \left(R^T e_3\right)\right) \\
& + M
\end{aligned}
$$
<!-- /model:eom -->

Multiplied by $-1$, this is Lee et al. eq. 6.20 with the book's $\rho_0 = -\rho$ (the
book measures from the centre of mass to the pivot) and the input added:
$J \dot\Omega_R + \Omega_R \times J\Omega_R + m g\, \rho \times R^T e_3 = M$.

## Code

```python
from geomech import SO3, Dot, Matrix, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
rho = Vector("\\rho", attr=["Constant"])
m, g = getScalars("m g", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
M = Vector("M")
R = SO3("R")
Om, eta = R.get_tangent_vector(), R.get_variation_vector()

L = Dot(Om, J * Om) * half - m * g * Dot(R * rho, e3)
dW = Dot(eta, M)
eom = compute_eom(L, dW, SystemVariables(matrices=[R]))
```

The same model is `geomech.systems.get("rigid_pendulum")`, with default parameters (a
1 m slender link), a default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-7}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-9}$, with $H(0) = -3.45$ |
| Closed form, Lee et al. eq. 6.20 | equal numerically at random states (`tests/test_systems_book.py`) |
<!-- /model:checks -->
