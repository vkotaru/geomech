# Double spherical pendulum

Two point masses on massless rods, the second hinged to the first by a spherical joint,
free to swing in every direction. The configuration is a pair of unit vectors
$(q_1, q_2) \in S^2 \times S^2$, the directions of the two links: no angles, so no
singular configurations. Lee et al. §5.6.5. Library system `double_spherical_pendulum`.

## Animation

Ten seconds, unforced, from the library's initial state: both links tilted to the same
side, slightly out of the $e_1$–$e_3$ plane, and swinging mostly in that plane. The
motion does not stay planar: both masses trace three-dimensional paths, with link 1
within $47°$ of hanging and link 2 within $59°$.
Drag to orbit; the full report (equations, standard form, parameters) opens in a
<a href="../model-reports/double_spherical_pendulum.html">separate page</a>.

<iframe class="model-report" src="../model-reports/double_spherical_pendulum.html?view" title="Double spherical pendulum, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 400 320" role="img" aria-labelledby="dsp-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="dsp-title">Double spherical pendulum: link directions q1 and q2, lengths l1 and l2, masses m1 and m2, angular velocities ω1 and ω2, input moments M1 and M2, gravity</title>
  <defs>
    <marker id="dsp-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- support -->
    <line x1="130" y1="50" x2="190" y2="50" stroke-width="2"/>
    <path d="M132,50 l-8,-10 M144,50 l-8,-10 M156,50 l-8,-10 M168,50 l-8,-10 M180,50 l-8,-10 M192,50 l-8,-10"/>
    <line x1="160" y1="50" x2="160" y2="190" stroke-dasharray="3 4" opacity="0.5"/>
    <!-- inertial axes -->
    <line x1="50" y1="290" x2="50" y2="250" marker-end="url(#dsp-arr)"/>
    <line x1="50" y1="290" x2="90" y2="290" marker-end="url(#dsp-arr)"/>
    <!-- links -->
    <line x1="160" y1="50" x2="225" y2="162.6" stroke-width="2.2"/>
    <line x1="225" y1="162.6" x2="315.6" y2="204.9" stroke-width="2.2"/>
  </g>
  <!-- q1, q2: unit vectors along the links -->
  <line x1="160" y1="50" x2="187.5" y2="97.6" stroke="#3f8fa0" stroke-width="3" marker-end="url(#dsp-arr)"/>
  <line x1="225" y1="162.6" x2="265.8" y2="181.6" stroke="#3f8fa0" stroke-width="3" marker-end="url(#dsp-arr)"/>
  <!-- omega1, omega2: perpendicular to q1, q2 (out of the page here) -->
  <circle cx="232" cy="128" r="7" fill="none" stroke="#3f8fa0" stroke-width="1.8"/>
  <circle cx="232" cy="128" r="2" fill="#3f8fa0"/>
  <circle cx="275" cy="165" r="7" fill="none" stroke="#3f8fa0" stroke-width="1.8"/>
  <circle cx="275" cy="165" r="2" fill="#3f8fa0"/>
  <!-- input moments on link 1 and link 2 -->
  <path d="M 133.7,59.6 A 28 28 0 0 0 150.4,76.3" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#dsp-arr)"/>
  <path d="M 216.1,187.0 A 26 26 0 0 0 241.7,182.5" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#dsp-arr)"/>
  <!-- gravity -->
  <line x1="360" y1="60" x2="360" y2="120" stroke="#e08a2c" stroke-width="2" marker-end="url(#dsp-arr)"/>
  <circle cx="160" cy="50" r="4" fill="currentColor"/>
  <circle cx="225" cy="162.6" r="11" fill="currentColor"/>
  <circle cx="315.6" cy="204.9" r="11" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="196" y="92" fill="#3f8fa0">q<tspan font-size="12" dy="4">1</tspan></text>
    <text x="244" y="204" fill="#3f8fa0">q<tspan font-size="12" dy="4">2</tspan></text>
    <text x="246" y="128" fill="#3f8fa0">ω<tspan font-size="12" dy="4">1</tspan></text>
    <text x="288" y="160" fill="#3f8fa0">ω<tspan font-size="12" dy="4">2</tspan></text>
    <text x="112" y="90" fill="#e08a2c">M<tspan font-size="12" dy="4">1</tspan></text>
    <text x="218" y="212" fill="#e08a2c">M<tspan font-size="12" dy="4">2</tspan></text>
    <text x="170" y="126">l<tspan font-size="12" dy="4">1</tspan></text>
    <text x="279" y="222">l<tspan font-size="12" dy="4">2</tspan></text>
    <text x="190" y="170">m<tspan font-size="12" dy="4">1</tspan></text>
    <text x="336" y="210">m<tspan font-size="12" dy="4">2</tspan></text>
    <text x="368" y="98" fill="#e08a2c">g</text>
    <text x="57" y="262">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="95" y="295">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>The first mass is at <span class="math notranslate nohighlight">\(l_1 q_1\)</span> from the pivot, the second at
<span class="math notranslate nohighlight">\(l_1 q_1 + l_2 q_2\)</span>. Each angular velocity <span class="math notranslate nohighlight">\(\omega_i\)</span> is perpendicular to its
<span class="math notranslate nohighlight">\(q_i\)</span>, out of the page in this picture, and <span class="math notranslate nohighlight">\(\dot q_i = \omega_i \times q_i\)</span>. Gravity acts along
<span class="math notranslate nohighlight">\(-e_3\)</span>; <span class="math notranslate nohighlight">\(M_1\)</span> and <span class="math notranslate nohighlight">\(M_2\)</span> are moments on link 1 and link 2, with virtual work
<span class="math notranslate nohighlight">\(\xi_1 \cdot M_1 + \xi_2 \cdot M_2\)</span>.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $q_1$ | direction of link 1, pivot to $m_1$ | $S^2$ |
| $q_2$ | direction of link 2, $m_1$ to $m_2$ | $S^2$ |
| $\omega_1$ | angular velocity of link 1, $\dot q_1 = \omega_1 \times q_1$, $q_1 \cdot \omega_1 = 0$ | $T_{q_1} S^2$ |
| $\omega_2$ | angular velocity of link 2, $\dot q_2 = \omega_2 \times q_2$, $q_2 \cdot \omega_2 = 0$ | $T_{q_2} S^2$ |
| $\xi_1$ | variation of $q_1$: $\delta q_1 = \xi_1 \times q_1$ | $T_{q_1} S^2$ |
| $\xi_2$ | variation of $q_2$: $\delta q_2 = \xi_2 \times q_2$ | $T_{q_2} S^2$ |
| $M_1$ | input moment on link 1 | $\mathbb{R}^3$ |
| $M_2$ | input moment on link 2 | $\mathbb{R}^3$ |
| $m_1, m_2, l_1, l_2, g, e_3$ | masses, link lengths, gravity, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m_1 \left(l_1 \left(\omega_1 \times q_1\right) \cdot l_1 \left(\omega_1 \times q_1\right)\right) \\
& + \frac{1}{2} m_2 \left(\left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right) \cdot \left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right)\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& \left(m_1 g \left(l_1 q_1 \cdot e_3\right) + m_2 g \left(\left(l_1 q_1 + l_2 q_2\right) \cdot e_3\right)\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m_1 \left(l_1 \left(\omega_1 \times q_1\right) \cdot l_1 \left(\omega_1 \times q_1\right)\right) \\
& + \frac{1}{2} m_2 \left(\left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right) \cdot \left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right)\right) \\
& -\left(m_1 g \left(l_1 q_1 \cdot e_3\right) + m_2 g \left(\left(l_1 q_1 + l_2 q_2\right) \cdot e_3\right)\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \xi_1 \cdot M_1 \\
& + \xi_2 \cdot M_2
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); the equation of $q_i$ holds on the tangent plane at
$q_i$.

<!-- model:eom -->
Coefficient of $\xi_1$, the equation of $q_1$:

$$
\begin{aligned}
0 = {}& -m_1 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -m_2 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -m_2 l_1 l_2 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_2 l_1 l_2 \left(q_1 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_2 l_1 l_2 \left(q_1 \times q_2\right) \\
& + 2 m_2 l_1 l_2 \left(\omega_1 \times \left(q_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_2 l_1 l_2 \left(q_1 \times \left(\omega_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -l_1 m_1 g \left(q_1 \times e_3\right) \\
& -l_1 m_2 g \left(q_1 \times e_3\right) \\
& + M_1
\end{aligned}
$$

Coefficient of $\xi_2$, the equation of $q_2$:

$$
\begin{aligned}
0 = {}& -m_2 l_2 l_1 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_2 l_2 l_1 \left(q_2 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_2 l_2 l_1 \left(q_2 \times q_1\right) \\
& -m_2 l_2 l_2 \left(q_2 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + 2 m_2 l_2 l_1 \left(\omega_2 \times \left(q_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_2 l_2 l_1 \left(q_2 \times \left(\omega_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -l_2 m_2 g \left(q_2 \times e_3\right) \\
& + M_2
\end{aligned}
$$
<!-- /model:eom -->

Multiplied by $-1$, these are the rows of Lee et al. eq. 5.71 with the inputs added,
equal on the tangent planes at $q_1$ and $q_2$ though the terms are grouped differently
(`tests/test_systems_book.py` compares them numerically).

## Code

```python
from geomech import S2, Dot, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
m1, m2 = getScalars("m1 m2", attr=["Constant"])
l1, l2 = getScalars("l1 l2", attr=["Constant"])
(g,) = getScalars("g", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
q1, q2 = S2("q1"), S2("q2")
M1, M2 = Vector("M1"), Vector("M2")

x1 = q1 * l1
x2 = x1 + q2 * l2
KE = sum(half * m * Dot(x.t_diff(), x.t_diff()) for m, x in [(m1, x1), (m2, x2)])
PE = sum(m * g * Dot(x, e3) for m, x in [(m1, x1), (m2, x2)])
L = KE - PE
dW = sum(Dot(q.get_variation_vector(), M) for q, M in [(q1, M1), (q2, M2)])
eom = compute_eom(L, dW, SystemVariables(vectors=[q1, q2]))
```

The same model is `geomech.systems.get("double_spherical_pendulum")`, with default
parameters, a default initial state and the 3-D scene used above; it is also
`get("n_link_spherical_chain", n_links=2)`.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-10}$, with $H(0) = -17.4$ |
| Closed form, Lee et al. eq. 5.71 | equal numerically at random states (`tests/test_systems_book.py`) |
<!-- /model:checks -->
