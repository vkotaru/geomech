# Heavy top

A rigid body spinning about a fixed pivot under gravity, its centre of mass a distance
$l$ from the pivot along the body axis $e_3$. The configuration is the attitude
$R \in SO(3)$, with body angular velocity $\Omega_{R}$ and $\dot R = R \widehat{\Omega_{R}}$.
No Euler angles are used, so the equations have no singular attitude. Lee et al. §6.6.2,
eq. 6.20. Library system `heavy_top`.

## Animation

Ten seconds, unforced, from the library's initial state: the axis tilted 0.3 rad
(17°) from the vertical, spinning at 20 rad/s about it. The default inertia is
axisymmetric, so the spin component of $\Omega_{R}$ stays at 20 rad/s. At this spin rate
the nutation is large: the axis nutates between 17° and 71° from the
vertical with a period of about 1.5 s, and precesses about 5.6 times round the vertical.
Drag to orbit; the full report (equations, standard form, parameters) opens in a
<a href="../model-reports/heavy_top.html">separate page</a>.

<iframe class="model-report" src="../model-reports/heavy_top.html?view" title="Heavy top, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 300" role="img" aria-labelledby="ht-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="ht-title">Heavy top: pivot O, body axis R e3, distance l, mass m, inertia J, spin Omega_R, gravity and input moment M</title>
  <defs>
    <marker id="ht-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- pivot support -->
    <path d="M170,250 L158,268 L182,268 Z"/>
    <line x1="146" y1="268" x2="194" y2="268" stroke-width="2"/>
    <!-- vertical through the pivot -->
    <line x1="170" y1="250" x2="170" y2="60" stroke-dasharray="3 4" opacity="0.5"/>
    <!-- inertial axes -->
    <line x1="50" y1="285" x2="50" y2="245" marker-end="url(#ht-arr)"/>
    <line x1="50" y1="285" x2="90" y2="285" marker-end="url(#ht-arr)"/>
    <!-- rod from the pivot to the centre of mass -->
    <line x1="170" y1="250" x2="240" y2="129" stroke-width="2.2"/>
    <!-- the top: a disc perpendicular to its axis -->
    <ellipse cx="240" cy="129" rx="60" ry="16" transform="rotate(30 240 129)" stroke-width="1.8"/>
  </g>
  <!-- body axis R e3 beyond the centre of mass -->
  <line x1="240" y1="129" x2="275" y2="68" stroke="#3f8fa0" stroke-width="3" marker-end="url(#ht-arr)"/>
  <!-- spin: angular velocity, mostly about the body axis -->
  <path d="M 246.4,81 A 18 6 30 0 0 277.6,99" fill="none" stroke="#3f8fa0" stroke-width="1.8" marker-end="url(#ht-arr)"/>
  <!-- gravity at the centre of mass -->
  <line x1="240" y1="129" x2="240" y2="215" stroke="#e08a2c" stroke-width="2" marker-end="url(#ht-arr)"/>
  <!-- input moment -->
  <path d="M 144,244 A 26 26 0 0 1 158,227" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#ht-arr)"/>
  <circle cx="170" cy="250" r="4" fill="currentColor"/>
  <circle cx="240" cy="129" r="3.5" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="188" y="258">O</text>
    <text x="280" y="62" fill="#3f8fa0">R e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="290" y="104" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">R</tspan></text>
    <text x="216" y="200">l</text>
    <text x="306" y="170">m</text>
    <text x="198" y="162">J</text>
    <text x="250" y="212" fill="#e08a2c">m g</text>
    <text x="118" y="232" fill="#e08a2c">M</text>
    <text x="44" y="238">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="95" y="290">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>The pivot O is fixed. The body axis, <span class="math notranslate nohighlight">\(R e_3\)</span> in the inertial frame,
carries the centre of mass at <span class="math notranslate nohighlight">\(l R e_3\)</span> from the pivot; <span class="math notranslate nohighlight">\(J\)</span> is the inertia
about the pivot. The top spins mostly about its axis (<span class="math notranslate nohighlight">\(\Omega_R\)</span>, body frame).
Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span>; the input <span class="math notranslate nohighlight">\(M\)</span> is a body-frame moment.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $R$ | attitude, body frame to inertial frame | $SO(3)$ |
| $\Omega_{R}$ | body angular velocity, $\dot R = R \widehat{\Omega_{R}}$ | $\mathbb{R}^3$ |
| $\eta_{R}$ | variation of $R$: $\delta R = R \widehat{\eta_{R}}$ | $\mathbb{R}^3$ |
| $M$ | input moment, body frame | $\mathbb{R}^3$ |
| $J$ | inertia about the pivot, body frame | constant, symmetric |
| $m, l, g, e_3$ | mass, pivot to centre of mass along body $e_3$, gravity, vertical unit vector | constants |

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
V = {}& m g l \left(R e_3 \cdot e_3\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right) \\
& -m g l \left(R e_3 \cdot e_3\right)
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
& -m g l \left(e_3 \times \left(R^T e_3\right)\right) \\
& + M
\end{aligned}
$$
<!-- /model:eom -->

Multiplied by $-1$, this is Lee et al. eq. 6.20 with the book's pivot offset
$\rho_0 = -l e_3$ and the input added:
$J \dot\Omega_R + \Omega_R \times J\Omega_R + m g l\, e_3 \times R^T e_3 = M$.

## Code

```python
from geomech import SO3, Dot, Matrix, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
m, g, l = getScalars("m g l", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
M = Vector("M")
R = SO3("R")
Om, eta = R.get_tangent_vector(), R.get_variation_vector()

L = Dot(Om, J * Om) * half - m * g * l * Dot(R * e3, e3)
dW = Dot(eta, M)
eom = compute_eom(L, dW, SystemVariables(matrices=[R]))
```

The same model is `geomech.systems.get("heavy_top")`, with default parameters (an
axisymmetric top), a default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-7}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-7}$, with $H(0) = 12.8$ |
| Closed form, Lee et al. eq. 6.20 | equal numerically at random states (`tests/test_systems_book.py`) |
<!-- /model:checks -->
