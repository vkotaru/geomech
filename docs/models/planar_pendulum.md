# Planar pendulum

A point mass $m$ on a massless rod of length $l$, swinging in a vertical plane about a
fixed pivot. The configuration is the direction of the rod, a unit vector $q$ on the
circle $S^1$ normal to the axis $e_2$, with a scalar angular velocity $\omega_q$: no
angle coordinate, so nothing to wrap at $\pm\pi$. Lee et al. §4.6.1. Library system
`planar_pendulum`.

## Animation

Ten seconds, unforced, from the library's initial state: 1 rad off the vertical with
$\omega_q = 1.5$ rad/s. The bob swings back and forth in the plane, up to $65°$ from
hanging on each side. Drag to orbit; the full report (equations, standard form,
parameters) opens in a
<a href="../model-reports/planar_pendulum.html">separate page</a>.

<iframe class="model-report" src="../model-reports/planar_pendulum.html?view" title="Planar pendulum, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 300" role="img" aria-labelledby="pp-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="pp-title">Planar pendulum: pivot, rod direction q, length l, bob m, angular velocity ω_q about e2, gravity and input moment τ</title>
  <defs>
    <marker id="pp-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- support -->
    <line x1="150" y1="50" x2="210" y2="50" stroke-width="2"/>
    <path d="M152,50 l-8,-10 M164,50 l-8,-10 M176,50 l-8,-10 M188,50 l-8,-10 M200,50 l-8,-10 M212,50 l-8,-10"/>
    <!-- circle of reachable positions -->
    <path d="M 15.4,145 A 190 190 0 0 0 344.6,145" stroke-dasharray="4 5" opacity="0.5"/>
    <line x1="180" y1="50" x2="180" y2="250" stroke-dasharray="3 4" opacity="0.5"/>
    <!-- inertial axes; e2 points into the page -->
    <circle cx="60" cy="270" r="7"/>
    <path d="M55,265 L65,275 M65,265 L55,275"/>
    <line x1="60" y1="263" x2="60" y2="225" marker-end="url(#pp-arr)"/>
    <line x1="67" y1="270" x2="105" y2="270" marker-end="url(#pp-arr)"/>
    <!-- rod -->
    <line x1="180" y1="50" x2="289" y2="206" stroke-width="2.2"/>
  </g>
  <!-- q: unit vector along the rod -->
  <line x1="180" y1="50" x2="220" y2="107" stroke="#3f8fa0" stroke-width="3" marker-end="url(#pp-arr)"/>
  <!-- omega: rotation about e2, clockwise in this view -->
  <path d="M 271.9,127.1 A 120 120 0 0 1 221.0,162.8" fill="none" stroke="#3f8fa0" stroke-width="2" marker-end="url(#pp-arr)"/>
  <!-- gravity on the bob -->
  <line x1="289" y1="206" x2="289" y2="270" stroke="#e08a2c" stroke-width="2" marker-end="url(#pp-arr)"/>
  <!-- input moment about e2 -->
  <path d="M 170.4,76.3 A 28 28 0 0 1 153.7,59.6" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#pp-arr)"/>
  <circle cx="180" cy="50" r="4" fill="currentColor"/>
  <circle cx="289" cy="206" r="13" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="214" y="64">O</text>
    <text x="222" y="96" fill="#3f8fa0">q</text>
    <text x="252" y="202">l</text>
    <text x="307" y="214">m</text>
    <text x="297" y="262" fill="#e08a2c">m g</text>
    <text x="140" y="86" fill="#e08a2c">τ</text>
    <text x="192" y="180" fill="#3f8fa0">ω<tspan font-size="12" dy="4">q</tspan></text>
    <text x="67" y="238">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="110" y="275">e<tspan font-size="12" dy="4">1</tspan></text>
    <text x="28" y="276">e<tspan font-size="12" dy="4">2</tspan></text>
  </g>
</svg>
<figcaption>The bob is at <span class="math notranslate nohighlight">\(l q\)</span>, on the circle of radius <span class="math notranslate nohighlight">\(l\)</span> about the pivot O
(dashed), in the plane normal to <span class="math notranslate nohighlight">\(e_2\)</span>. The axis <span class="math notranslate nohighlight">\(e_2\)</span> points into the page, so
<span class="math notranslate nohighlight">\(\dot q = \omega_q\, e_2 \times q\)</span> turns <span class="math notranslate nohighlight">\(q\)</span> clockwise in this picture when
<span class="math notranslate nohighlight">\(\omega_q &gt; 0\)</span>. Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span>; the input <span class="math notranslate nohighlight">\(\tau\)</span> is a moment about
<span class="math notranslate nohighlight">\(e_2\)</span>.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $q$ | direction of the rod, pivot to bob, $q \cdot e_2 = 0$ | $S^1$ |
| $\omega_q$ | angular velocity about $e_2$, $\dot q = \omega_q\, e_2 \times q$ | $\mathbb{R}$ |
| $\gamma_q$ | variation of $q$: $\delta q = \gamma_q\, e_2 \times q$ | $\mathbb{R}$ |
| $\tau$ | input moment about $e_2$ | $\mathbb{R}$ |
| $m, l, g$ | mass, rod length, gravity | constants |
| $e_2, e_3$ | axis of the plane of motion, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m l l \omega_{q} \omega_{q}
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& m g l \left(q \cdot e_3\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m l l \omega_{q} \omega_{q} \\
& -m g l \left(q \cdot e_3\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \gamma_{q} \tau
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); on $S^1$ it is a scalar equation, the coefficient
of the scalar variation $\gamma_q$.

<!-- model:eom -->
Coefficient of $\gamma_{q}$, the equation of $q$:

$$
\begin{aligned}
0 = {}& -m l l \dot{\omega_{q}} \\
& -m g l \left(e_2 \times q \cdot e_3\right) \\
& + \tau
\end{aligned}
$$
<!-- /model:eom -->

Multiplied by $-1$, this is Lee et al. eq. 4.33 with the input added,
$m l^2 \dot\omega_q + m g l\, e_1^T q = \tau$, where the book's $e_1^T q$ (in its planar
frame) is $(e_2 \times q) \cdot e_3$ here.

## Code

```python
from geomech import S1, Dot, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
m, g, l = getScalars("m g l", attr=["Constant"])
e2 = Vector("e2", attr=["Constant", "UnitNorm"])
e3 = Vector("e3", attr=["Constant"])
q, tau = S1("q", axis=e2), Scalar("\\tau")
omega, gamma = q.get_tangent_vector(), q.get_variation_vector()

L = half * m * l * l * omega * omega - m * g * l * Dot(q, e3)
dW = gamma * tau
eom = compute_eom(L, dW, SystemVariables(vectors=[q]))
```

The same model is `geomech.systems.get("planar_pendulum")`, with default parameters
($e_2 = (0, 1, 0)$), a default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-7}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-10}$, with $H(0) = -4.18$ |
| Closed form, Lee et al. eq. 4.33 | equal numerically at random states (`tests/test_systems_book.py`) |
<!-- /model:checks -->
