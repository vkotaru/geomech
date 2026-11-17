# Quadrotors carrying a point mass

$n$ quadrotors carry one point-mass load $m_L$, each on its own taut massless cable. Cable
$i$ has length $l_i$ and direction $q_i \in S^2$ from quadrotor $i$ to the load, so
quadrotor $i$ is at $x_L - l_i q_i$. The configuration manifold is
$\mathbb{R}^3 \times (SO(3) \times S^2)^n$. The equations below are derived once for
every $n$ with `geomech.indexed`, then for $n = 5$ with the library system, which the
animation uses. Not taken from the literature; the equations are verified numerically
(`check_eom`). Library system `multi_quadrotor_point_mass`.

## Animation

Six seconds, five quadrotors, open loop, with no controller. The start is not the
library's initial state: each quadrotor starts at rest, tilted 2.0° outward so that its
constant thrust (15.7 N) balances its weight and the pull of its cable: with the cables at
the library's initial directions (29° off the vertical, leaning in) and the load at rest,
this is a static equilibrium. The load and cables start with the library's initial
velocities, so the system moves away from it. With $M_i = 0$ the
attitude equations do not involve the load, and each vehicle keeps its tilt. The cables
swing between 6° and 56° from the vertical, the vehicles come within 0.38 m of each other
(centre to centre; nothing models a collision), and the load drifts 2.3 m sideways while
it rises from 0.5 m to 2.1 m.
Drag to orbit; the full report (equations, standard form, parameters) opens in a
<a href="../model-reports/multi_quadrotor_point_mass.html">separate page</a>.

<iframe class="model-report" src="../model-reports/multi_quadrotor_point_mass.html?view" title="Quadrotors carrying a point mass, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 310" role="img" aria-labelledby="mpm-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="mpm-title">Quadrotors carrying a point mass: load m_L at x_L; cable i with direction q_i from quadrotor i to the load, length l_i; quadrotor i with attitude R_i, thrust f_i and moment M_i</title>
  <defs>
    <marker id="mpm-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <line x1="30" y1="290" x2="30" y2="250" marker-end="url(#mpm-arr)"/>
    <line x1="30" y1="290" x2="70" y2="290" marker-end="url(#mpm-arr)"/>
  </g>
  <g fill="none" stroke="currentColor" stroke-width="1.6">
    <line x1="80" y1="95" x2="190" y2="232"/>
    <line x1="310" y1="95" x2="190" y2="232"/>
    <line x1="210" y1="52" x2="190" y2="232"/>
  </g>
  <!-- quadrotor at (80,95) -->
  <g fill="none" stroke="currentColor" stroke-width="1.6">
    <line x1="42" y1="95" x2="118" y2="95"/>
    <line x1="60" y1="105" x2="100" y2="85"/>
    <ellipse cx="42" cy="95" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="118" cy="95" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="100" cy="85" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="60" cy="105" rx="13" ry="4" stroke-width="1.2"/>
  </g>
  <rect x="71" y="90" width="18" height="10" rx="2" fill="currentColor"/>
  <line x1="80.0" y1="87.0" x2="80.0" y2="47.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#mpm-arr)"/>
  <path d="M 67,65.0 A 13 4 0 1 0 93,65.0" fill="none" stroke="#e08a2c" stroke-width="1.6" marker-end="url(#mpm-arr)"/>
  <!-- quadrotor at (310,95) -->
  <g fill="none" stroke="currentColor" stroke-width="1.6">
    <line x1="272" y1="95" x2="348" y2="95"/>
    <line x1="290" y1="105" x2="330" y2="85"/>
    <ellipse cx="272" cy="95" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="348" cy="95" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="330" cy="85" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="290" cy="105" rx="13" ry="4" stroke-width="1.2"/>
  </g>
  <rect x="301" y="90" width="18" height="10" rx="2" fill="currentColor"/>
  <line x1="310.0" y1="87.0" x2="310.0" y2="47.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#mpm-arr)"/>
  <!-- quadrotor at (210,52) -->
  <g fill="none" stroke="currentColor" stroke-width="1.6">
    <line x1="172" y1="52" x2="248" y2="52"/>
    <line x1="190" y1="62" x2="230" y2="42"/>
    <ellipse cx="172" cy="52" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="248" cy="52" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="230" cy="42" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="190" cy="62" rx="13" ry="4" stroke-width="1.2"/>
  </g>
  <rect x="201" y="47" width="18" height="10" rx="2" fill="currentColor"/>
  <line x1="210.0" y1="44.0" x2="210.0" y2="8.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#mpm-arr)"/>
  <line x1="80.0" y1="95.0" x2="108.2" y2="130.1" stroke="#3f8fa0" stroke-width="3" marker-end="url(#mpm-arr)"/>
  <line x1="310.0" y1="95.0" x2="280.3" y2="128.9" stroke="#3f8fa0" stroke-width="3" marker-end="url(#mpm-arr)"/>
  <line x1="210.0" y1="52.0" x2="205.0" y2="96.7" stroke="#3f8fa0" stroke-width="3" marker-end="url(#mpm-arr)"/>
  <circle cx="190" cy="232" r="12" fill="currentColor"/>
  <line x1="190.0" y1="244.0" x2="190.0" y2="290.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#mpm-arr)"/>
  <circle cx="160" cy="156" r="7" fill="none" stroke="#3f8fa0" stroke-width="1.6"/>
  <circle cx="160" cy="156" r="2" fill="#3f8fa0"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="60" y="52" fill="#e08a2c">f<tspan font-size="12" dy="4">1</tspan></text>
    <text x="98" y="68" fill="#e08a2c">M<tspan font-size="12" dy="4">1</tspan></text>
    <text x="316" y="58" fill="#e08a2c">f<tspan font-size="12" dy="4">2</tspan></text>
    <text x="216" y="20" fill="#e08a2c">f<tspan font-size="12" dy="4">3</tspan></text>
    <text x="8" y="130" fill="#3f8fa0">R<tspan font-size="12" dy="4">1</tspan><tspan dy="-4">, Ω</tspan><tspan font-size="12" dy="4">1</tspan></text>
    <text x="8" y="154">m<tspan font-size="12" dy="4">1</tspan><tspan dy="-4">, J</tspan><tspan font-size="12" dy="4">1</tspan></text>
    <text x="171" y="161" fill="#3f8fa0">ω<tspan font-size="12" dy="4">1</tspan></text>
    <text x="24" y="242">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="76" y="295">e<tspan font-size="12" dy="4">1</tspan></text>
    <text x="118" y="126" fill="#3f8fa0">q<tspan font-size="12" dy="4">1</tspan></text>
    <text x="110" y="182">l<tspan font-size="12" dy="4">1</tspan></text>
    <text x="250" y="142" fill="#3f8fa0">q<tspan font-size="12" dy="4">2</tspan></text>
    <text x="217" y="100" fill="#3f8fa0">q<tspan font-size="12" dy="4">3</tspan></text>
    <text x="208" y="238">m<tspan font-size="12" dy="4">L</tspan></text>
    <text x="158" y="238" fill="#3f8fa0">x<tspan font-size="12" dy="4">L</tspan></text>
    <text x="198" y="286" fill="#e08a2c">m<tspan font-size="12" dy="4">L</tspan><tspan dy="-4"> g</tspan></text>
  </g>
</svg>
<figcaption>Drawn with three quadrotors (five on this page). Quadrotor 1 is labelled in
full; the others are the same with their own index.
Each <span class="math notranslate nohighlight">\(q_i\)</span> points from quadrotor <span class="math notranslate nohighlight">\(i\)</span> to the load, and quadrotor
<span class="math notranslate nohighlight">\(i\)</span> is at <span class="math notranslate nohighlight">\(x_L - l_i q_i\)</span>. <span class="math notranslate nohighlight">\(\omega_1 \perp q_1\)</span> is out of the page here. Each thrust <span class="math notranslate nohighlight">\(f_i\)</span> acts along
<span class="math notranslate nohighlight">\(R_i e_3\)</span>; each vehicle also takes a body-frame moment <span class="math notranslate nohighlight">\(M_i\)</span>
(shown for the first). Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span> on every mass (shown on the load only).</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $x_L$ | position of the load | $\mathbb{R}^3$ |
| $\dot x_L$ | velocity of the load | $\mathbb{R}^3$ |
| $\delta x_L$ | variation of $x_L$ | $\mathbb{R}^3$ |
| $q_i$ | direction of cable $i$, quadrotor $i$ to load | $S^2$ |
| $\omega_i$ | angular velocity of cable $i$, $\dot q_i = \omega_i \times q_i$, $q_i \cdot \omega_i = 0$ | $T_{q_i} S^2$ |
| $\xi_i$ | variation of $q_i$: $\delta q_i = \xi_i \times q_i$ | $T_{q_i} S^2$ |
| $R_i$ | attitude of quadrotor $i$, body frame to inertial frame | $SO(3)$ |
| $\Omega_i$ | body angular velocity of quadrotor $i$, $\dot R_i = R_i \widehat{\Omega_i}$ | $\mathbb{R}^3$ |
| $\eta_i$ | variation of $R_i$: $\delta R_i = R_i \widehat{\eta_i}$ | $\mathbb{R}^3$ |
| $f_i$ | thrust of quadrotor $i$ along $R_i e_3$ | $\mathbb{R}$ |
| $M_i$ | moment of quadrotor $i$, body frame | $\mathbb{R}^3$ |
| $m_L, m_i, J_i, l_i, g, e_3$ | load mass, quadrotor masses and inertias, cable lengths, gravity, vertical unit vector | constants |

## Lagrangian and equations for every $n$

<!-- model:indexed -->
Lagrangian for every $n$, as `geomech.indexed.systems.multi_quadrotor_point_mass()` holds it (the velocities expanded, so the inner sums are the inertia couplings):

$$
\begin{aligned}
L = {}& \frac{1}{2} m_L \left(\dot{x_L} \cdot \dot{x_L}\right) \\
& -m_L g \left(e_3 \cdot x_L\right) \\
& + \sum_{j=1}^{n} \frac{1}{2} m_{j} \left(\left(\dot{x_L} -l_{j} \left(\omega_{j} \times q_{j}\right)\right) \cdot \left(\dot{x_L} -l_{j} \left(\omega_{j} \times q_{j}\right)\right)\right) \\
& + \sum_{j=1}^{n} \frac{1}{2} \left(\Omega_{j} \cdot J_{j} \Omega_{j}\right) \\
& -\sum_{j=1}^{n} m_{j} g \left(e_3 \cdot \left(x_L -l_{j} q_{j}\right)\right)
\end{aligned}
$$

Virtual work:

$$
\begin{aligned}
\delta W = {}& \sum_{j=1}^{n} \left(\delta x_L -l_{j} \left(\xi_{j} \times q_{j}\right)\right) \cdot f_{j} R_{j} e_3 \\
& + \sum_{j=1}^{n} \eta_{j} \cdot M_{j}
\end{aligned}
$$

Equations from `compute_indexed_eom`, derived once for every $n$; $i$ is the body, $1 \le i \le n$:

Coefficient of $\delta x_L$, the equation of $x_L$:

$$
\begin{aligned}
0 = {}& -m_L \ddot{x_L} \\
& -m_L g e_3 \\
& -\sum_{j=1}^{n} m_{j} \ddot{x_L} \\
& + \sum_{j=1}^{n} m_{j} l_{j} \left(\dot{\omega_{j}} \times q_{j}\right) \\
& -\sum_{j=1}^{n} \left(\omega_{j} \cdot \omega_{j}\right) m_{j} l_{j} q_{j} \\
& -\sum_{j=1}^{n} m_{j} g e_3 \\
& + \sum_{j=1}^{n} f_{j} R_{j} e_3
\end{aligned}
$$

Coefficient of $\xi_{i}$, the equation of $q_i$:

$$
\begin{aligned}
0 = {}& m_{i} l_{i} \left(\left(\omega_{i} \times q_{i}\right) \times \dot{x_L}\right) \\
& + m_{i} l_{i} \left(q_{i} \times \ddot{x_L}\right) \\
& -m_{i} l_{i} l_{i} \left(q_{i} \times \left(\dot{\omega_{i}} \times q_{i}\right)\right) \\
& -2 m_{i} l_{i} \left(\omega_{i} \times \left(q_{i} \times \dot{x_L}\right)\right) \\
& + m_{i} l_{i} \left(q_{i} \times \left(\omega_{i} \times \dot{x_L}\right)\right) \\
& + m_{i} g l_{i} \left(q_{i} \times e_3\right) \\
& -l_{i} f_{i} \left(q_{i} \times \left(R_{i} e_3\right)\right)
\end{aligned}
$$

Coefficient of $\eta_{i}$, the equation of $R_i$:

$$
\begin{aligned}
0 = {}& -J_{i} \dot{\Omega_{i}} \\
& -\Omega_{i} \times \left(J_{i} \Omega_{i}\right) \\
& + M_{i}
\end{aligned}
$$
<!-- /model:indexed -->

## Energy, Lagrangian and virtual work at $n = 5$

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m_L \left(\dot{x_L} \cdot \dot{x_L}\right) \\
& + \frac{1}{2} m_1 \left(\left(\dot{x_L} -l_1 \left(\omega_1 \times q_1\right)\right) \cdot \left(\dot{x_L} -l_1 \left(\omega_1 \times q_1\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_1 \cdot J_1 \Omega_1\right) \\
& + \frac{1}{2} m_2 \left(\left(\dot{x_L} -l_2 \left(\omega_2 \times q_2\right)\right) \cdot \left(\dot{x_L} -l_2 \left(\omega_2 \times q_2\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_2 \cdot J_2 \Omega_2\right) \\
& + \frac{1}{2} m_3 \left(\left(\dot{x_L} -l_3 \left(\omega_3 \times q_3\right)\right) \cdot \left(\dot{x_L} -l_3 \left(\omega_3 \times q_3\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_3 \cdot J_3 \Omega_3\right) \\
& + \frac{1}{2} m_4 \left(\left(\dot{x_L} -l_4 \left(\omega_4 \times q_4\right)\right) \cdot \left(\dot{x_L} -l_4 \left(\omega_4 \times q_4\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_4 \cdot J_4 \Omega_4\right) \\
& + \frac{1}{2} m_5 \left(\left(\dot{x_L} -l_5 \left(\omega_5 \times q_5\right)\right) \cdot \left(\dot{x_L} -l_5 \left(\omega_5 \times q_5\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_5 \cdot J_5 \Omega_5\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& m_L g \left(e_3 \cdot x_L\right) \\
& + m_1 g \left(e_3 \cdot \left(x_L -l_1 q_1\right)\right) \\
& + m_2 g \left(e_3 \cdot \left(x_L -l_2 q_2\right)\right) \\
& + m_3 g \left(e_3 \cdot \left(x_L -l_3 q_3\right)\right) \\
& + m_4 g \left(e_3 \cdot \left(x_L -l_4 q_4\right)\right) \\
& + m_5 g \left(e_3 \cdot \left(x_L -l_5 q_5\right)\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m_L \left(\dot{x_L} \cdot \dot{x_L}\right) \\
& -m_L g \left(e_3 \cdot x_L\right) \\
& + \frac{1}{2} m_1 \left(\left(\dot{x_L} -l_1 \left(\omega_1 \times q_1\right)\right) \cdot \left(\dot{x_L} -l_1 \left(\omega_1 \times q_1\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_1 \cdot J_1 \Omega_1\right) \\
& -m_1 g \left(e_3 \cdot \left(x_L -l_1 q_1\right)\right) \\
& + \frac{1}{2} m_2 \left(\left(\dot{x_L} -l_2 \left(\omega_2 \times q_2\right)\right) \cdot \left(\dot{x_L} -l_2 \left(\omega_2 \times q_2\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_2 \cdot J_2 \Omega_2\right) \\
& -m_2 g \left(e_3 \cdot \left(x_L -l_2 q_2\right)\right) \\
& + \frac{1}{2} m_3 \left(\left(\dot{x_L} -l_3 \left(\omega_3 \times q_3\right)\right) \cdot \left(\dot{x_L} -l_3 \left(\omega_3 \times q_3\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_3 \cdot J_3 \Omega_3\right) \\
& -m_3 g \left(e_3 \cdot \left(x_L -l_3 q_3\right)\right) \\
& + \frac{1}{2} m_4 \left(\left(\dot{x_L} -l_4 \left(\omega_4 \times q_4\right)\right) \cdot \left(\dot{x_L} -l_4 \left(\omega_4 \times q_4\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_4 \cdot J_4 \Omega_4\right) \\
& -m_4 g \left(e_3 \cdot \left(x_L -l_4 q_4\right)\right) \\
& + \frac{1}{2} m_5 \left(\left(\dot{x_L} -l_5 \left(\omega_5 \times q_5\right)\right) \cdot \left(\dot{x_L} -l_5 \left(\omega_5 \times q_5\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_5 \cdot J_5 \Omega_5\right) \\
& -m_5 g \left(e_3 \cdot \left(x_L -l_5 q_5\right)\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \left(\delta x_L -l_1 \left(\xi_1 \times q_1\right)\right) \cdot f_1 R_1 e_3 \\
& + \eta_1 \cdot M_1 \\
& + \left(\delta x_L -l_2 \left(\xi_2 \times q_2\right)\right) \cdot f_2 R_2 e_3 \\
& + \eta_2 \cdot M_2 \\
& + \left(\delta x_L -l_3 \left(\xi_3 \times q_3\right)\right) \cdot f_3 R_3 e_3 \\
& + \eta_3 \cdot M_3 \\
& + \left(\delta x_L -l_4 \left(\xi_4 \times q_4\right)\right) \cdot f_4 R_4 e_3 \\
& + \eta_4 \cdot M_4 \\
& + \left(\delta x_L -l_5 \left(\xi_5 \times q_5\right)\right) \cdot f_5 R_5 e_3 \\
& + \eta_5 \cdot M_5
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion at $n = 5$

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); each $R_i$ equation is in the body frame and each
$q_i$ equation holds on the tangent plane at $q_i$.

<!-- model:eom -->
<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\delta x_L\)</span>, the equation of <span class="math notranslate nohighlight">\(x_L\)</span> (27 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_L \ddot{x_L} \\
& -m_L g e_3 \\
& -m_1 \ddot{x_L} \\
& + m_1 l_1 \left(\dot{\omega_1} \times q_1\right) \\
& -\left(\omega_1 \cdot \omega_1\right) m_1 l_1 q_1 \\
& -m_1 g e_3 \\
& -m_2 \ddot{x_L} \\
& + m_2 l_2 \left(\dot{\omega_2} \times q_2\right) \\
& -\left(\omega_2 \cdot \omega_2\right) m_2 l_2 q_2 \\
& -m_2 g e_3 \\
& -m_3 \ddot{x_L} \\
& + m_3 l_3 \left(\dot{\omega_3} \times q_3\right) \\
& -\left(\omega_3 \cdot \omega_3\right) m_3 l_3 q_3 \\
& -m_3 g e_3 \\
& -m_4 \ddot{x_L} \\
& + m_4 l_4 \left(\dot{\omega_4} \times q_4\right) \\
& -\left(\omega_4 \cdot \omega_4\right) m_4 l_4 q_4 \\
& -m_4 g e_3 \\
& -m_5 \ddot{x_L} \\
& + m_5 l_5 \left(\dot{\omega_5} \times q_5\right) \\
& -\left(\omega_5 \cdot \omega_5\right) m_5 l_5 q_5 \\
& -m_5 g e_3 \\
& + f_1 R_1 e_3 \\
& + f_2 R_2 e_3 \\
& + f_3 R_3 e_3 \\
& + f_4 R_4 e_3 \\
& + f_5 R_5 e_3
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_1\)</span>, the equation of <span class="math notranslate nohighlight">\(q_1\)</span> (7 terms)</summary>

$$
\begin{aligned}
0 = {}& m_1 l_1 \left(\left(\omega_1 \times q_1\right) \times \dot{x_L}\right) \\
& + m_1 l_1 \left(q_1 \times \ddot{x_L}\right) \\
& -m_1 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -2 m_1 l_1 \left(\omega_1 \times \left(q_1 \times \dot{x_L}\right)\right) \\
& + m_1 l_1 \left(q_1 \times \left(\omega_1 \times \dot{x_L}\right)\right) \\
& + m_1 g l_1 \left(q_1 \times e_3\right) \\
& -l_1 f_1 \left(q_1 \times \left(R_1 e_3\right)\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_2\)</span>, the equation of <span class="math notranslate nohighlight">\(q_2\)</span> (7 terms)</summary>

$$
\begin{aligned}
0 = {}& m_2 l_2 \left(\left(\omega_2 \times q_2\right) \times \dot{x_L}\right) \\
& + m_2 l_2 \left(q_2 \times \ddot{x_L}\right) \\
& -m_2 l_2 l_2 \left(q_2 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& -2 m_2 l_2 \left(\omega_2 \times \left(q_2 \times \dot{x_L}\right)\right) \\
& + m_2 l_2 \left(q_2 \times \left(\omega_2 \times \dot{x_L}\right)\right) \\
& + m_2 g l_2 \left(q_2 \times e_3\right) \\
& -l_2 f_2 \left(q_2 \times \left(R_2 e_3\right)\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_3\)</span>, the equation of <span class="math notranslate nohighlight">\(q_3\)</span> (7 terms)</summary>

$$
\begin{aligned}
0 = {}& m_3 l_3 \left(\left(\omega_3 \times q_3\right) \times \dot{x_L}\right) \\
& + m_3 l_3 \left(q_3 \times \ddot{x_L}\right) \\
& -m_3 l_3 l_3 \left(q_3 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& -2 m_3 l_3 \left(\omega_3 \times \left(q_3 \times \dot{x_L}\right)\right) \\
& + m_3 l_3 \left(q_3 \times \left(\omega_3 \times \dot{x_L}\right)\right) \\
& + m_3 g l_3 \left(q_3 \times e_3\right) \\
& -l_3 f_3 \left(q_3 \times \left(R_3 e_3\right)\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_4\)</span>, the equation of <span class="math notranslate nohighlight">\(q_4\)</span> (7 terms)</summary>

$$
\begin{aligned}
0 = {}& m_4 l_4 \left(\left(\omega_4 \times q_4\right) \times \dot{x_L}\right) \\
& + m_4 l_4 \left(q_4 \times \ddot{x_L}\right) \\
& -m_4 l_4 l_4 \left(q_4 \times \left(\dot{\omega_4} \times q_4\right)\right) \\
& -2 m_4 l_4 \left(\omega_4 \times \left(q_4 \times \dot{x_L}\right)\right) \\
& + m_4 l_4 \left(q_4 \times \left(\omega_4 \times \dot{x_L}\right)\right) \\
& + m_4 g l_4 \left(q_4 \times e_3\right) \\
& -l_4 f_4 \left(q_4 \times \left(R_4 e_3\right)\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_5\)</span>, the equation of <span class="math notranslate nohighlight">\(q_5\)</span> (7 terms)</summary>

$$
\begin{aligned}
0 = {}& m_5 l_5 \left(\left(\omega_5 \times q_5\right) \times \dot{x_L}\right) \\
& + m_5 l_5 \left(q_5 \times \ddot{x_L}\right) \\
& -m_5 l_5 l_5 \left(q_5 \times \left(\dot{\omega_5} \times q_5\right)\right) \\
& -2 m_5 l_5 \left(\omega_5 \times \left(q_5 \times \dot{x_L}\right)\right) \\
& + m_5 l_5 \left(q_5 \times \left(\omega_5 \times \dot{x_L}\right)\right) \\
& + m_5 g l_5 \left(q_5 \times e_3\right) \\
& -l_5 f_5 \left(q_5 \times \left(R_5 e_3\right)\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_1\)</span>, the equation of <span class="math notranslate nohighlight">\(R_1\)</span> (3 terms)</summary>

$$
\begin{aligned}
0 = {}& -J_1 \dot{\Omega_1} \\
& -\Omega_1 \times \left(J_1 \Omega_1\right) \\
& + M_1
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_2\)</span>, the equation of <span class="math notranslate nohighlight">\(R_2\)</span> (3 terms)</summary>

$$
\begin{aligned}
0 = {}& -J_2 \dot{\Omega_2} \\
& -\Omega_2 \times \left(J_2 \Omega_2\right) \\
& + M_2
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_3\)</span>, the equation of <span class="math notranslate nohighlight">\(R_3\)</span> (3 terms)</summary>

$$
\begin{aligned}
0 = {}& -J_3 \dot{\Omega_3} \\
& -\Omega_3 \times \left(J_3 \Omega_3\right) \\
& + M_3
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_4\)</span>, the equation of <span class="math notranslate nohighlight">\(R_4\)</span> (3 terms)</summary>

$$
\begin{aligned}
0 = {}& -J_4 \dot{\Omega_4} \\
& -\Omega_4 \times \left(J_4 \Omega_4\right) \\
& + M_4
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_5\)</span>, the equation of <span class="math notranslate nohighlight">\(R_5\)</span> (3 terms)</summary>

$$
\begin{aligned}
0 = {}& -J_5 \dot{\Omega_5} \\
& -\Omega_5 \times \left(J_5 \Omega_5\right) \\
& + M_5
\end{aligned}
$$

</details>
<!-- /model:eom -->

## Code

The model for every $n$, as `geomech.indexed.systems.multi_quadrotor_point_mass()` writes
it, and the library system for $n = 5$:

```python
from geomech import S2, SO3, Cross, Dot, Matrix, Scalar, Vector
from geomech.core.base.constants import HALF
from geomech.indexed import Family, IndexedVariables, Sum, compute_indexed_eom, indices, n
from geomech.systems import get

(i,) = indices("i")
mL, g = Scalar("m_L", attr=["Constant"]), Scalar("g", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
xL = Vector("x_L")
m, l = Family(Scalar, "m", ["Constant"]), Family(Scalar, "l", ["Constant"])
J = Family(Matrix, "J", attr=["Constant", "SymmetricMatrix"])
R, q = Family(SO3, "R"), Family(S2, "q")  # quadrotor attitudes, cable directions
f, M = Family(Scalar, "f"), Family(Vector, "M")
vL = xL.t_diff()
xQ = xL - q[i] * l[i]  # quadrotor i
Om = R[i].get_tangent_vector()

L = HALF * mL * Dot(vL, vL) - mL * g * Dot(e3, xL)
L = L + Sum(
    i,
    1,
    n,
    HALF * m[i] * Dot(xQ.t_diff(), xQ.t_diff())
    + HALF * Dot(Om, J[i] * Om)
    - m[i] * g * Dot(e3, xQ),
)
dxQ = xL.delta() - Cross(q[i].get_variation_vector(), q[i]) * l[i]
dW = Sum(i, 1, n, Dot(dxQ, (R[i] * e3) * f[i]) + Dot(R[i].get_variation_vector(), M[i]))
variables = IndexedVariables(families=[q, R], vectors=[xL])
eom = compute_indexed_eom(L, dW, variables)  # every n

system = get("multi_quadrotor_point_mass", n_quads=5)  # n from 1 to 20
eom_5 = system.derive()
```

The same model at $n = 5$ is `geomech.systems.get("multi_quadrotor_point_mass", n_quads=5)`,
with default parameters, a default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-5}$ |
| Energy $H = T + V$, 6 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-9}$, with $H(0) = 105$ |
| Equations for every $n$ above, at $n = 5$: `check_eom` against the Euler–Lagrange evaluation of the generic $L$, $\delta W$ at $n = 5$ | below $10^{-5}$ |
| Generic $L$ at $n = 5$ against the library system's $L$, at random states | relative difference below $10^{-14}$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
