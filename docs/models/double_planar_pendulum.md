# Double planar pendulum

Two point masses on massless rods, the second hinged to the first, swinging in a
vertical plane. The configuration is a pair of unit vectors $(q_1, q_2)$ on the circle
$S^1$ normal to the axis $e_2$, each with a scalar angular velocity, so the configuration
manifold is $S^1 \times S^1$ and the equations need no angle coordinates. Lee et al.
§4.6.4. Library system `double_planar_pendulum`.

## Animation

Ten seconds, unforced, from the library's initial state: link 1 at 0.8 rad and link 2
at 1.6 rad from hanging, with $\omega_1 = 0.5$ and $\omega_2 = -1.2$ rad/s. The motion is
irregular: link 1 stays within $66°$ of hanging, while link 2 swings up to $108°$, above
the horizontal. Drag to orbit; the full report (equations, standard form, parameters)
opens in a <a href="../model-reports/double_planar_pendulum.html">separate page</a>.

<iframe class="model-report" src="../model-reports/double_planar_pendulum.html?view" title="Double planar pendulum, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 400 320" role="img" aria-labelledby="dpp-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="dpp-title">Double planar pendulum: link directions q1 and q2, lengths l1 and l2, masses m1 and m2, angular velocities ω1 and ω2 about e2, input moments τ1 and τ2, gravity</title>
  <defs>
    <marker id="dpp-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- support -->
    <line x1="130" y1="50" x2="190" y2="50" stroke-width="2"/>
    <path d="M132,50 l-8,-10 M144,50 l-8,-10 M156,50 l-8,-10 M168,50 l-8,-10 M180,50 l-8,-10 M192,50 l-8,-10"/>
    <line x1="160" y1="50" x2="160" y2="190" stroke-dasharray="3 4" opacity="0.5"/>
    <!-- inertial axes; e2 points into the page -->
    <circle cx="50" cy="290" r="7"/>
    <path d="M45,285 L55,295 M55,285 L45,295"/>
    <line x1="50" y1="283" x2="50" y2="250" marker-end="url(#dpp-arr)"/>
    <line x1="57" y1="290" x2="92" y2="290" marker-end="url(#dpp-arr)"/>
    <!-- links -->
    <line x1="160" y1="50" x2="225" y2="162.6" stroke-width="2.2"/>
    <line x1="225" y1="162.6" x2="315.6" y2="204.9" stroke-width="2.2"/>
  </g>
  <!-- q1, q2: unit vectors along the links -->
  <line x1="160" y1="50" x2="187.5" y2="97.6" stroke="#3f8fa0" stroke-width="3" marker-end="url(#dpp-arr)"/>
  <line x1="225" y1="162.6" x2="265.8" y2="181.6" stroke="#3f8fa0" stroke-width="3" marker-end="url(#dpp-arr)"/>
  <!-- omega1, omega2: rotation of each link about e2, clockwise in this view -->
  <path d="M 190.9,58.3 A 32 32 0 0 1 182.6,72.6" fill="none" stroke="#3f8fa0" stroke-width="2" marker-end="url(#dpp-arr)"/>
  <path d="M 234.6,136.3 A 28 28 0 0 1 251.3,153.0" fill="none" stroke="#3f8fa0" stroke-width="2" marker-end="url(#dpp-arr)"/>
  <!-- input moments about e2 on link 1 and link 2 -->
  <path d="M 150.4,76.3 A 28 28 0 0 1 133.7,59.6" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#dpp-arr)"/>
  <path d="M 241.7,182.5 A 26 26 0 0 1 216.1,187.0" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#dpp-arr)"/>
  <!-- gravity -->
  <line x1="360" y1="60" x2="360" y2="120" stroke="#e08a2c" stroke-width="2" marker-end="url(#dpp-arr)"/>
  <circle cx="160" cy="50" r="4" fill="currentColor"/>
  <circle cx="225" cy="162.6" r="11" fill="currentColor"/>
  <circle cx="315.6" cy="204.9" r="11" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="196" y="92" fill="#3f8fa0">q<tspan font-size="12" dy="4">1</tspan></text>
    <text x="270" y="166" fill="#3f8fa0">q<tspan font-size="12" dy="4">2</tspan></text>
    <text x="200" y="60" fill="#3f8fa0">ω<tspan font-size="12" dy="4">1</tspan></text>
    <text x="258" y="136" fill="#3f8fa0">ω<tspan font-size="12" dy="4">2</tspan></text>
    <text x="112" y="90" fill="#e08a2c">τ<tspan font-size="12" dy="4">1</tspan></text>
    <text x="218" y="212" fill="#e08a2c">τ<tspan font-size="12" dy="4">2</tspan></text>
    <text x="170" y="126">l<tspan font-size="12" dy="4">1</tspan></text>
    <text x="279" y="222">l<tspan font-size="12" dy="4">2</tspan></text>
    <text x="190" y="170">m<tspan font-size="12" dy="4">1</tspan></text>
    <text x="336" y="210">m<tspan font-size="12" dy="4">2</tspan></text>
    <text x="368" y="98" fill="#e08a2c">g</text>
    <text x="57" y="262">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="97" y="295">e<tspan font-size="12" dy="4">1</tspan></text>
    <text x="18" y="296">e<tspan font-size="12" dy="4">2</tspan></text>
  </g>
</svg>
<figcaption>The first mass is at <span class="math notranslate nohighlight">\(l_1 q_1\)</span> from the pivot, the second at
<span class="math notranslate nohighlight">\(l_1 q_1 + l_2 q_2\)</span>. Both links swing in the plane normal to <span class="math notranslate nohighlight">\(e_2\)</span>,
which points into the page: <span class="math notranslate nohighlight">\(\dot q_i = \omega_i\, e_2 \times q_i\)</span> turns link <span class="math notranslate nohighlight">\(i\)</span> clockwise in this
picture when <span class="math notranslate nohighlight">\(\omega_i &gt; 0\)</span>. Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span>; <span class="math notranslate nohighlight">\(\tau_1\)</span> and <span class="math notranslate nohighlight">\(\tau_2\)</span> are moments about <span class="math notranslate nohighlight">\(e_2\)</span> on link 1
and link 2.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $q_1$ | direction of link 1, pivot to $m_1$, $q_1 \cdot e_2 = 0$ | $S^1$ |
| $q_2$ | direction of link 2, $m_1$ to $m_2$, $q_2 \cdot e_2 = 0$ | $S^1$ |
| $\omega_1$ | angular velocity of link 1 about $e_2$, $\dot q_1 = \omega_1\, e_2 \times q_1$ | $\mathbb{R}$ |
| $\omega_2$ | angular velocity of link 2 about $e_2$, $\dot q_2 = \omega_2\, e_2 \times q_2$ | $\mathbb{R}$ |
| $\gamma_1$ | variation of $q_1$: $\delta q_1 = \gamma_1\, e_2 \times q_1$ | $\mathbb{R}$ |
| $\gamma_2$ | variation of $q_2$: $\delta q_2 = \gamma_2\, e_2 \times q_2$ | $\mathbb{R}$ |
| $\tau_1$ | input moment about $e_2$ on link 1 | $\mathbb{R}$ |
| $\tau_2$ | input moment about $e_2$ on link 2 | $\mathbb{R}$ |
| $m_1, m_2, l_1, l_2, g$ | masses, link lengths, gravity | constants |
| $e_2, e_3$ | axis of the plane of motion, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m_1 \left(l_1 \omega_1 \left(e_2 \times q_1\right) \cdot l_1 \omega_1 \left(e_2 \times q_1\right)\right) \\
& + \frac{1}{2} m_2 \left(\left(l_1 \omega_1 \left(e_2 \times q_1\right) + l_2 \omega_2 \left(e_2 \times q_2\right)\right) \cdot \left(l_1 \omega_1 \left(e_2 \times q_1\right) + l_2 \omega_2 \left(e_2 \times q_2\right)\right)\right)
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
L = {}& \frac{1}{2} m_1 \left(l_1 \omega_1 \left(e_2 \times q_1\right) \cdot l_1 \omega_1 \left(e_2 \times q_1\right)\right) \\
& + \frac{1}{2} m_2 \left(\left(l_1 \omega_1 \left(e_2 \times q_1\right) + l_2 \omega_2 \left(e_2 \times q_2\right)\right) \cdot \left(l_1 \omega_1 \left(e_2 \times q_1\right) + l_2 \omega_2 \left(e_2 \times q_2\right)\right)\right) \\
& -\left(m_1 g \left(l_1 q_1 \cdot e_3\right) + m_2 g \left(\left(l_1 q_1 + l_2 q_2\right) \cdot e_3\right)\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \gamma_1 \tau_1 \\
& + \gamma_2 \tau_2
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); on $S^1$ each is a scalar equation, the coefficient
of the scalar variation $\gamma_i$.

<!-- model:eom -->
Coefficient of $\gamma_1$, the equation of $q_1$:

$$
\begin{aligned}
0 = {}& -m_1 l_1 \dot{\omega_1} l_1 \\
& -m_2 l_1 \dot{\omega_1} l_1 \\
& -m_2 l_1 \dot{\omega_2} l_2 \left(q_1 \cdot q_2\right) \\
& -m_2 l_1 \omega_2 l_2 \omega_2 \left(q_1 \cdot e_2 \times q_2\right) \\
& -m_1 g l_1 \left(e_2 \times q_1 \cdot e_3\right) \\
& -m_2 g l_1 \left(e_2 \times q_1 \cdot e_3\right) \\
& + \tau_1
\end{aligned}
$$

Coefficient of $\gamma_2$, the equation of $q_2$:

$$
\begin{aligned}
0 = {}& -m_2 l_2 \dot{\omega_1} l_1 \left(q_2 \cdot q_1\right) \\
& + m_2 l_2 \omega_1 l_1 \omega_1 \left(q_1 \cdot e_2 \times q_2\right) \\
& -m_2 l_2 \dot{\omega_2} l_2 \\
& -m_2 g l_2 \left(e_2 \times q_2 \cdot e_3\right) \\
& + \tau_2
\end{aligned}
$$
<!-- /model:eom -->

Multiplied by $-1$, these are the rows of Lee et al. eq. 4.62 with the inputs added,
except for the velocity terms: the book prints $\omega_2(\omega_2 - \omega_1)$ and
$\omega_1(\omega_1 - \omega_2)$, which leave out the kinetic part of
$q_i^T S\, \partial L / \partial q_i$ in its general form eq. 4.12, and the derived
$\omega_2^2$ and $\omega_1^2$ are those of the textbook double pendulum in angles (both
comparisons are in `tests/test_systems_book.py`).

## Code

```python
from geomech import S1, Dot, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
m1, m2, l1, l2, g = getScalars("m1 m2 l1 l2 g", attr=["Constant"])
e2 = Vector("e2", attr=["Constant", "UnitNorm"])
e3 = Vector("e3", attr=["Constant"])
q1, q2 = S1("q1", axis=e2), S1("q2", axis=e2)
tau1, tau2 = Scalar("\\tau_1"), Scalar("\\tau_2")

x1 = q1 * l1
x2 = x1 + q2 * l2
v1, v2 = x1.t_diff(), x2.t_diff()
KE = half * m1 * Dot(v1, v1) + half * m2 * Dot(v2, v2)
PE = m1 * g * Dot(x1, e3) + m2 * g * Dot(x2, e3)
L = KE - PE
dW = q1.get_variation_vector() * tau1 + q2.get_variation_vector() * tau2
eom = compute_eom(L, dW, SystemVariables(vectors=[q1, q2]))
```

The same model is `geomech.systems.get("double_planar_pendulum")`, with default
parameters, a default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-8}$, with $H(0) = -11.9$ |
| Closed form, Lee et al. eq. 4.62 | equal numerically at random states (`tests/test_systems_book.py`) |
<!-- /model:checks -->
