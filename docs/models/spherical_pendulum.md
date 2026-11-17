# Spherical pendulum

A point mass $m$ on a massless rod of length $l$, free to swing in every direction about
a fixed pivot. The configuration is the direction of the rod, a unit vector $q \in S^2$:
no angles, so no singularity at the bottom or the top. Lee et al. §5.6.1. Library system
`spherical_pendulum`.

## Animation

Ten seconds, unforced, from a start that is not the library's initial state: 1 rad off
the vertical, with $\omega_q$ out of the vertical plane, so the bob traces a
three-dimensional rosette. The library's initial state ($\omega_q$ almost along $e_2$)
swings in one vertical plane instead. Drag to orbit;
the full report (equations, standard form, parameters) opens in a
<a href="../model-reports/spherical_pendulum.html">separate page</a>.

<iframe class="model-report" src="../model-reports/spherical_pendulum.html?view" title="Spherical pendulum, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 300" role="img" aria-labelledby="sp-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="sp-title">Spherical pendulum: pivot, rod direction q, length l, bob m, angular velocity ω_q, gravity and input moment M</title>
  <defs>
    <marker id="sp-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- support -->
    <line x1="150" y1="50" x2="210" y2="50" stroke-width="2"/>
    <path d="M152,50 l-8,-10 M164,50 l-8,-10 M176,50 l-8,-10 M188,50 l-8,-10 M200,50 l-8,-10 M212,50 l-8,-10"/>
    <!-- sphere of reachable positions -->
    <path d="M 15.4,145 A 190 190 0 0 0 344.6,145" stroke-dasharray="4 5" opacity="0.5"/>
    <line x1="180" y1="50" x2="180" y2="250" stroke-dasharray="3 4" opacity="0.5"/>
    <!-- inertial axes -->
    <line x1="60" y1="270" x2="60" y2="225" marker-end="url(#sp-arr)"/>
    <line x1="60" y1="270" x2="105" y2="270" marker-end="url(#sp-arr)"/>
    <!-- rod -->
    <line x1="180" y1="50" x2="289" y2="206" stroke-width="2.2"/>
  </g>
  <!-- q: unit vector along the rod -->
  <line x1="180" y1="50" x2="220" y2="107" stroke="#3f8fa0" stroke-width="3" marker-end="url(#sp-arr)"/>
  <!-- gravity on the bob -->
  <line x1="289" y1="206" x2="289" y2="270" stroke="#e08a2c" stroke-width="2" marker-end="url(#sp-arr)"/>
  <!-- input moment about the pivot -->
  <path d="M 153.7,59.6 A 28 28 0 0 0 170.4,76.3" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#sp-arr)"/>
  <!-- omega: perpendicular to q (out of the page here) -->
  <circle cx="320" cy="175" r="8" fill="none" stroke="#3f8fa0" stroke-width="1.8"/>
  <circle cx="320" cy="175" r="2.2" fill="#3f8fa0"/>
  <circle cx="180" cy="50" r="4" fill="currentColor"/>
  <circle cx="289" cy="206" r="13" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="214" y="64">O</text>
    <text x="222" y="96" fill="#3f8fa0">q</text>
    <text x="250" y="140">l</text>
    <text x="307" y="214">m</text>
    <text x="297" y="262" fill="#e08a2c">m g</text>
    <text x="140" y="86" fill="#e08a2c">M</text>
    <text x="333" y="171" fill="#3f8fa0">ω<tspan font-size="12" dy="4">q</tspan></text>
    <text x="54" y="218">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="110" y="275">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>The bob is at <span class="math notranslate nohighlight">\(l q\)</span>, on the sphere of radius <span class="math notranslate nohighlight">\(l\)</span> about the pivot O
(dashed). The angular velocity <span class="math notranslate nohighlight">\(\omega_q\)</span> is perpendicular to <span class="math notranslate nohighlight">\(q\)</span>, out of the page in this
picture, and <span class="math notranslate nohighlight">\(\dot q = \omega_q \times q\)</span>. Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span>; the input <span class="math notranslate nohighlight">\(M\)</span> is a moment
about the pivot.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $q$ | direction of the rod, pivot to bob | $S^2$ |
| $\omega_q$ | angular velocity of the rod, $\dot q = \omega_q \times q$, $q \cdot \omega_q = 0$ | $T_q S^2$ |
| $\xi_q$ | variation of $q$: $\delta q = \xi_q \times q$ | $T_q S^2$ |
| $M$ | input moment about the pivot | $\mathbb{R}^3$ |
| $m, l, g, e_3$ | mass, rod length, gravity, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m l l \left(\omega_{q} \cdot \omega_{q}\right)
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
L = {}& \frac{1}{2} m l l \left(\omega_{q} \cdot \omega_{q}\right) \\
& -m g l \left(q \cdot e_3\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \xi_{q} \cdot M
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); it holds on the tangent plane at $q$.

<!-- model:eom -->
Coefficient of $\xi_{q}$, the equation of $q$:

$$
\begin{aligned}
0 = {}& -m l l \dot{\omega_{q}} \\
& -m g l \left(q \times e_3\right) \\
& + M
\end{aligned}
$$
<!-- /model:eom -->

Multiplied by $-1$, this is Lee et al. eq. 5.38 with the input added:
$m l^2 \dot\omega_q + m g l\, q \times e_3 = M$.

## Code

```python
from geomech import S2, Dot, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
m, g, l = getScalars("m g l", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
q, M = S2("q"), Vector("M")
omega, xi = q.get_tangent_vector(), q.get_variation_vector()

L = half * m * l * l * Dot(omega, omega) - m * g * l * Dot(q, e3)
dW = Dot(xi, M)
eom = compute_eom(L, dW, SystemVariables(vectors=[q]))
```

The same model is `geomech.systems.get("spherical_pendulum")`, with default parameters,
a default initial state (1 rad off the vertical with $\omega_q$ almost along $e_2$, a
swing in the $e_1$–$e_3$ plane; not the start animated above) and the 3-D scene used
above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-10}$, with $H(0) = -2.58$ |
| Closed form, Lee et al. eq. 5.38 | equal numerically at random states (`tests/test_systems_book.py`) |
<!-- /model:checks -->
