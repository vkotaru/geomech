# Particle on a rotating sphere

A particle of mass $m$ that slides without friction on a sphere of radius $L$, while the
sphere turns about the vertical axis $e_3$ at a constant rate $\Omega$. The model is
written in the frame that turns with the sphere, so the configuration is the particle's
direction $q \in S^2$ in that frame, and the rotation enters the Lagrangian as a
Coriolis coupling and a centrifugal term. There is no gravity. On this page $L$ is also
the sphere radius, as in Lee et al. §5.6.2, besides the Lagrangian, and $\Omega$ is a
scalar rotation rate. Lee et al. §5.6.2. Library system `particle_on_rotating_sphere`.

## Animation

Ten seconds, unforced, from the library's initial state, shown in the frame that turns
with the sphere ($\Omega = 2$ rad/s). The Coriolis and centrifugal terms bend the
particle's path: it crosses the equator again and again and stays within $52°$ of it in
latitude. Drag to orbit; the full report (equations, standard form, parameters) opens in
a <a href="../model-reports/particle_on_rotating_sphere.html">separate page</a>.

<iframe class="model-report" src="../model-reports/particle_on_rotating_sphere.html?view" title="Particle on a rotating sphere, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 300" role="img" aria-labelledby="prs-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="prs-title">Particle on a rotating sphere: sphere of radius L turning at rate Ω about e3, particle direction q, mass m, angular velocity ω_q, input moment M</title>
  <defs>
    <marker id="prs-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- sphere and equator -->
    <circle cx="190" cy="165" r="110" stroke-width="1.6"/>
    <path d="M 80,165 A 110 26 0 0 0 300,165" stroke-dasharray="4 5" opacity="0.6"/>
    <path d="M 80,165 A 110 26 0 0 1 300,165" stroke-dasharray="2 6" opacity="0.35"/>
    <!-- rotation axis -->
    <line x1="190" y1="55" x2="190" y2="275" stroke-dasharray="3 4" opacity="0.5"/>
    <line x1="190" y1="55" x2="190" y2="25" marker-end="url(#prs-arr)"/>
    <!-- sphere rotation about e3 -->
    <path d="M 172,40 A 18 6 0 0 0 208,40" marker-end="url(#prs-arr)"/>
    <!-- radius to the particle -->
    <line x1="190" y1="165" x2="258.2" y2="110"/>
  </g>
  <!-- q: unit vector from the centre to the particle -->
  <line x1="190" y1="165" x2="225" y2="136.7" stroke="#3f8fa0" stroke-width="3" marker-end="url(#prs-arr)"/>
  <!-- omega: tangent to the sphere at q -->
  <line x1="265.1" y1="118.6" x2="283" y2="141" stroke="#3f8fa0" stroke-width="2" marker-end="url(#prs-arr)"/>
  <!-- input moment about the centre -->
  <path d="M 167.4,156.8 A 24 24 0 0 1 181.8,142.4" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#prs-arr)"/>
  <circle cx="190" cy="165" r="3.5" fill="currentColor"/>
  <circle cx="258.2" cy="110" r="8" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="170" y="24">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="214" y="46">Ω</text>
    <text x="216" y="168" fill="#3f8fa0">q</text>
    <text x="210" y="128">L</text>
    <text x="240" y="96">m</text>
    <text x="262" y="156" fill="#3f8fa0">ω<tspan font-size="12" dy="4">q</tspan></text>
    <text x="150" y="128" fill="#e08a2c">M</text>
  </g>
</svg>
<figcaption>The particle is at <span class="math notranslate nohighlight">\(L q\)</span> on the sphere, with <span class="math notranslate nohighlight">\(q\)</span> taken in the frame
that turns with the sphere at the rate <span class="math notranslate nohighlight">\(\Omega\)</span> about <span class="math notranslate nohighlight">\(e_3\)</span>. Its angular velocity <span class="math notranslate nohighlight">\(\omega_q\)</span> in
that frame is perpendicular to <span class="math notranslate nohighlight">\(q\)</span>, with <span class="math notranslate nohighlight">\(\dot q = \omega_q \times q\)</span>. The input
<span class="math notranslate nohighlight">\(M\)</span> is a moment about the centre.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $q$ | direction of the particle from the centre, in the sphere's frame | $S^2$ |
| $\omega_q$ | angular velocity in the sphere's frame, $\dot q = \omega_q \times q$, $q \cdot \omega_q = 0$ | $T_q S^2$ |
| $\xi_q$ | variation of $q$: $\delta q = \xi_q \times q$ | $T_q S^2$ |
| $M$ | input moment about the centre | $\mathbb{R}^3$ |
| $m, L, \Omega$ | mass, sphere radius, rotation rate | constants |
| $e_3, C$ | rotation axis, $C = \mathrm{diag}(1, 1, 0)$ | constants |

## Energy, Lagrangian and virtual work

The Lagrangian is a single term with a velocity in it, so all of it appears as $T$
below and $V = 0$. The energy that the run conserves is the Jacobi integral
$\frac{1}{2} m L^2 (\omega_q \cdot \omega_q - \Omega^2\, q \cdot C q)$, not $T$.

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m L L \left(\omega_{q} \cdot \omega_{q} + 2 \left(q \cdot \widehat{\Omega e_3} \widehat{q} \omega_{q}\right) + \Omega \Omega \left(q \cdot C q\right)\right)
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
L = {}& \frac{1}{2} m L L \left(\omega_{q} \cdot \omega_{q} + 2 \left(q \cdot \widehat{\Omega e_3} \widehat{q} \omega_{q}\right) + \Omega \Omega \left(q \cdot C q\right)\right)
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
0 = {}& -m L L \dot{\omega_{q}} \\
& + m L L \Omega \left(q \times \left(e_3 \times \left(q \times \omega_{q}\right)\right)\right) \\
& -m L L \Omega \left(q \times \left(\omega_{q} \times \left(e_3 \times q\right)\right)\right) \\
& -m L L \Omega \left(\left(\left(\omega_{q} \times q\right) \times e_3\right) \times q\right) \\
& -m L L \Omega \left(\left(q \times e_3\right) \times \left(\omega_{q} \times q\right)\right) \\
& + 2 m L L \Omega \left(\omega_{q} \times \left(q \times \left(e_3 \times q\right)\right)\right) \\
& + m L L \Omega \Omega \left(q \times \left(C q\right)\right) \\
& + M
\end{aligned}
$$
<!-- /model:eom -->

Lee et al. eq. 5.43 gives the same dynamics in $(q, \dot q)$ form; with
$\dot q = \omega_q \times q$, the equation above is $-q \times (\text{eq. 5.43}) + M$ on
the tangent plane at $q$.

## Code

```python
from geomech import S2, Dot, Hat, Matrix, MVMul, Scalar, SystemVariables, Vector
from geomech import compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
m, L_, Om = getScalars("m L \\Omega", attr=["Constant"])  # L_ is the radius L
e3 = Vector("e3", attr=["Constant"])
C = Matrix("C", attr=["Constant", "SymmetricMatrix"])  # diag(1, 1, 0)
q, M = S2("q"), Vector("M")
omega, xi = q.get_tangent_vector(), q.get_variation_vector()

coupling = Dot(q, MVMul(Hat(e3 * Om), MVMul(Hat(q), omega))) * 2
L = half * m * L_ * L_ * (Dot(omega, omega) + coupling + Om * Om * Dot(q, C * q))
dW = Dot(xi, M)
eom = compute_eom(L, dW, SystemVariables(vectors=[q]))
```

The same model is `geomech.systems.get("particle_on_rotating_sphere")`, with default
parameters, a default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-7}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-11}$, with $H(0) = -0.753$ |
| Closed form, Lee et al. eq. 5.43 | equal numerically at random states (`tests/test_systems_book.py`) |
<!-- /model:checks -->
