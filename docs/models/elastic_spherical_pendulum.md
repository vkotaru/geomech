# Elastic spherical pendulum

A point mass $m$ on a linear spring of stiffness $k$ and rest length $l_0$, hinged at a
fixed pivot and free to swing in every direction. The configuration is the spring length
$l \in \mathbb{R}$ and its direction $q \in S^2$, so the manifold is
$\mathbb{R} \times S^2$: a scalar and a unit vector, with no angles. Not taken from the
literature; the equations are verified numerically (`check_eom`). Library system
`elastic_spherical_pendulum`.

## Animation

Ten seconds, unforced, from a start that is not the library's initial state: 1 rad off
the vertical, the spring stretched to $l = 1.1$
($l_0 = 1$) and lengthening at $\dot l = 0.2$, with $\omega_q$ out of the vertical
plane. The bob traces a three-dimensional path, between $23°$ and $65°$ from hanging,
while the spring length oscillates between $1.06$ and $1.61$. The library's initial
state (the same, but with $\omega_q$ almost along $e_2$) swings in one vertical plane
instead. Drag to orbit; the full report
(equations, standard form, parameters) opens in a
<a href="../model-reports/elastic_spherical_pendulum.html">separate page</a>.

<iframe class="model-report" src="../model-reports/elastic_spherical_pendulum.html?view" title="Elastic spherical pendulum, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 300" role="img" aria-labelledby="esp-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="esp-title">Elastic spherical pendulum: pivot, spring direction q, length l and its rate, stiffness k, bob m, angular velocity ω_q, gravity, input force f_l along the spring and input moment M_q</title>
  <defs>
    <marker id="esp-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- support -->
    <line x1="150" y1="50" x2="210" y2="50" stroke-width="2"/>
    <path d="M152,50 l-8,-10 M164,50 l-8,-10 M176,50 l-8,-10 M188,50 l-8,-10 M200,50 l-8,-10 M212,50 l-8,-10"/>
    <line x1="180" y1="50" x2="180" y2="250" stroke-dasharray="3 4" opacity="0.5"/>
    <!-- inertial axes -->
    <line x1="60" y1="270" x2="60" y2="225" marker-end="url(#esp-arr)"/>
    <line x1="60" y1="270" x2="105" y2="270" marker-end="url(#esp-arr)"/>
    <!-- spring -->
    <path d="M180.0,50.0 L223.0,111.5 L233.8,113.0 L229.3,134.5 L251.0,137.6 L246.5,159.1 L268.2,162.2 L263.7,183.7 L274.5,185.3 L281.4,195.1" stroke-width="1.8"/>
  </g>
  <!-- l: length of the spring, pivot to bob -->
  <line x1="158.7" y1="64.9" x2="267.7" y2="220.9" stroke="#3f8fa0" stroke-width="1.6" marker-start="url(#esp-arr)" marker-end="url(#esp-arr)"/>
  <!-- q: unit vector along the spring -->
  <line x1="180" y1="50" x2="220" y2="107" stroke="#3f8fa0" stroke-width="3" marker-end="url(#esp-arr)"/>
  <!-- rate of change of l -->
  <line x1="307.0" y1="193.4" x2="324.2" y2="218.0" stroke="#3f8fa0" stroke-width="2" marker-end="url(#esp-arr)"/>
  <!-- omega: perpendicular to q (out of the page here) -->
  <circle cx="310" cy="150" r="8" fill="none" stroke="#3f8fa0" stroke-width="1.8"/>
  <circle cx="310" cy="150" r="2.2" fill="#3f8fa0"/>
  <!-- gravity on the bob -->
  <line x1="289" y1="221" x2="289" y2="275" stroke="#e08a2c" stroke-width="2" marker-end="url(#esp-arr)"/>
  <!-- input force along the spring -->
  <line x1="297.6" y1="218.3" x2="318" y2="247.5" stroke="#e08a2c" stroke-width="2" marker-end="url(#esp-arr)"/>
  <!-- input moment about the pivot -->
  <path d="M 153.7,59.6 A 28 28 0 0 0 170.4,76.3" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#esp-arr)"/>
  <circle cx="180" cy="50" r="4" fill="currentColor"/>
  <circle cx="289" cy="206" r="13" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="214" y="64">O</text>
    <text x="222" y="96" fill="#3f8fa0">q</text>
    <text x="198" y="158" fill="#3f8fa0">l</text>
    <text x="262" y="140">k</text>
    <text x="300" y="180">m</text>
    <text x="330" y="214" fill="#3f8fa0">l̇</text>
    <text x="324" y="156" fill="#3f8fa0">ω<tspan font-size="12" dy="4">q</tspan></text>
    <text x="236" y="272" fill="#e08a2c">m g</text>
    <text x="324" y="256" fill="#e08a2c">f<tspan font-size="12" dy="4">l</tspan></text>
    <text x="132" y="86" fill="#e08a2c">M<tspan font-size="12" dy="4">q</tspan></text>
    <text x="54" y="218">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="110" y="275">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>The bob is at <span class="math notranslate nohighlight">\(l q\)</span> from the pivot O. The spring pulls it back toward
the rest length <span class="math notranslate nohighlight">\(l_0\)</span> (not drawn) with stiffness <span class="math notranslate nohighlight">\(k\)</span>. The angular velocity
<span class="math notranslate nohighlight">\(\omega_q\)</span> is perpendicular to <span class="math notranslate nohighlight">\(q\)</span>, out of the page in this picture, and
<span class="math notranslate nohighlight">\(\dot q = \omega_q \times q\)</span>. Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span>; the input <span class="math notranslate nohighlight">\(f_l\)</span> is a
force along the spring and <span class="math notranslate nohighlight">\(M_q\)</span> a moment about the pivot.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $l$ | length of the spring, pivot to bob | $\mathbb{R}$ |
| $\dot l$ | rate of change of the length | $\mathbb{R}$ |
| $\delta l$ | variation of $l$ | $\mathbb{R}$ |
| $q$ | direction of the spring, pivot to bob | $S^2$ |
| $\omega_q$ | angular velocity of the spring, $\dot q = \omega_q \times q$, $q \cdot \omega_q = 0$ | $T_q S^2$ |
| $\xi_q$ | variation of $q$: $\delta q = \xi_q \times q$ | $T_q S^2$ |
| $f_l$ | input force along the spring | $\mathbb{R}$ |
| $M_q$ | input moment about the pivot | $\mathbb{R}^3$ |
| $m, k, l_0, g, e_3$ | mass, stiffness, rest length, gravity, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m \left(\dot{l} \dot{l} + l l \left(\omega_{q} \cdot \omega_{q}\right)\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& \left(m g l \left(q \cdot e_3\right) + \frac{1}{2} k \left(l -l_0\right) \left(l -l_0\right)\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m \left(\dot{l} \dot{l} + l l \left(\omega_{q} \cdot \omega_{q}\right)\right) \\
& -\left(m g l \left(q \cdot e_3\right) + \frac{1}{2} k \left(l -l_0\right) \left(l -l_0\right)\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \delta l f_l \\
& + \xi_{q} \cdot M_q
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); the equation of $l$ is scalar, and the equation of
$q$ holds on the tangent plane at $q$.

<!-- model:eom -->
Coefficient of $\delta l$, the equation of $l$:

$$
\begin{aligned}
0 = {}& -m \ddot{l} \\
& + m l \left(\omega_{q} \cdot \omega_{q}\right) \\
& -m g \left(q \cdot e_3\right) \\
& -k l \\
& + k l_0 \\
& + f_l
\end{aligned}
$$

Coefficient of $\xi_{q}$, the equation of $q$:

$$
\begin{aligned}
0 = {}& -m l l \dot{\omega_{q}} \\
& -2 m l \dot{l} \omega_{q} \\
& -m g l \left(q \times e_3\right) \\
& + M_q
\end{aligned}
$$
<!-- /model:eom -->

## Code

```python
from geomech import S2, Dot, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
m, g, k, l0 = getScalars("m g k l0", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
l, q = Scalar("l"), S2("q")
omega, xi = q.get_tangent_vector(), q.get_variation_vector()
f_l, M_q = Scalar("f_l"), Vector("M_q")

ld = l.t_diff()
KE = half * m * (ld * ld + l * l * Dot(omega, omega))
PE = m * g * l * Dot(q, e3) + half * k * (l - l0) * (l - l0)
L = KE - PE
dW = l.delta() * f_l + Dot(xi, M_q)
eom = compute_eom(L, dW, SystemVariables(scalars=[l], vectors=[q]))
```

The same model is `geomech.systems.get("elastic_spherical_pendulum")`, with default
parameters, a default initial state (not the start animated above) and the 3-D scene
used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-9}$, with $H(0) = -2.32$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
