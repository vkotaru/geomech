# Point mass

A particle of mass $m$ moving freely in space under an applied force $F$. The
configuration is its position $x \in \mathbb{R}^3$, a vector space, so the equation of
motion is Newton's law $m \ddot x = F$ with no geometry to handle. It is the simplest
system in the library and a baseline for the others. Lee et al. §3.2. Library system
`point_mass`.

## Animation

Ten seconds, unforced, from the library's initial state: with $F = 0$ the mass moves on
a straight line at its initial velocity $\dot x = (1, 0.5, -0.4)$. Drag to orbit;
the full report (equations, standard form, parameters) opens in a
<a href="../model-reports/point_mass.html">separate page</a>.

<iframe class="model-report" src="../model-reports/point_mass.html?view" title="Point mass, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 260" role="img" aria-labelledby="pm-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="pm-title">Point mass: origin O, position x, mass m, velocity ẋ and applied force F</title>
  <defs>
    <marker id="pm-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- inertial axes -->
    <line x1="40" y1="240" x2="40" y2="205" marker-end="url(#pm-arr)"/>
    <line x1="40" y1="240" x2="75" y2="240" marker-end="url(#pm-arr)"/>
    <line x1="40" y1="240" x2="62" y2="222" marker-end="url(#pm-arr)"/>
  </g>
  <!-- x: position from the origin -->
  <line x1="110" y1="200" x2="257.6" y2="98.5" stroke="#3f8fa0" stroke-width="3" marker-end="url(#pm-arr)"/>
  <!-- velocity -->
  <line x1="283.4" y1="83.3" x2="330" y2="60" stroke="#3f8fa0" stroke-width="2" marker-end="url(#pm-arr)"/>
  <!-- applied force -->
  <line x1="276.7" y1="103.4" x2="300" y2="150" stroke="#e08a2c" stroke-width="2" marker-end="url(#pm-arr)"/>
  <circle cx="110" cy="200" r="3.5" fill="currentColor"/>
  <circle cx="270" cy="90" r="13" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="90" y="216">O</text>
    <text x="178" y="134" fill="#3f8fa0">x</text>
    <text x="240" y="74">m</text>
    <text x="336" y="62" fill="#3f8fa0">ẋ</text>
    <text x="306" y="156" fill="#e08a2c">F</text>
    <text x="47" y="206">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="69" y="218">e<tspan font-size="12" dy="4">2</tspan></text>
    <text x="80" y="245">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>The mass is at <span class="math notranslate nohighlight">\(x\)</span>, measured from a fixed origin O, and moves with
velocity <span class="math notranslate nohighlight">\(\dot x\)</span>. The applied force <span class="math notranslate nohighlight">\(F\)</span> is the only force: there is no gravity in
this model.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $x$ | position of the mass | $\mathbb{R}^3$ |
| $\dot x$ | velocity | $\mathbb{R}^3$ |
| $\delta x$ | variation of $x$ | $\mathbb{R}^3$ |
| $F$ | applied force | $\mathbb{R}^3$ |
| $m$ | mass | constant |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m \left(\dot{x} \cdot \dot{x}\right)
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
L = {}& \frac{1}{2} m \left(\dot{x} \cdot \dot{x}\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \delta x \cdot F
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); here it is one vector equation in $\mathbb{R}^3$.

<!-- model:eom -->
Coefficient of $\delta x$, the equation of $x$:

$$
\begin{aligned}
0 = {}& -m \ddot{x} \\
& + F
\end{aligned}
$$
<!-- /model:eom -->

## Code

```python
from geomech import Dot, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
(m,) = getScalars("m", attr=["Constant"])
x, F = Vector("x"), Vector("F")
v = x.t_diff()

L = half * m * Dot(v, v)
dW = Dot(x.delta(), F)
eom = compute_eom(L, dW, SystemVariables(vectors=[x]))
```

The same model is `geomech.systems.get("point_mass")`, with default parameters, a
default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-7}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-14}$, with $H(0) = 1.06$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
