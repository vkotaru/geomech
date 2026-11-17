# Quadrotors carrying a rigid body

$n$ quadrotors carry one rigid-body load, each on its own taut massless cable. The load
has centre of mass $x_L$, attitude $R_L$, mass $m_L$ and inertia $J_L$; cable $i$ is
attached at $x_L + R_L \rho_i$, has length $l_i$ and direction $q_i \in S^2$ from
quadrotor $i$ to the attachment point. The configuration manifold is
$\mathbb{R}^3 \times SO(3) \times (SO(3) \times S^2)^n$. The equations below are
derived once for every $n$ with `geomech.indexed`, then for $n = 5$ with the library
system, which the animation uses. Not taken from the literature; the equations are
verified numerically (`check_eom`). Library system `multi_quadrotor_rigid_load`.

## Animation

Six seconds, five quadrotors, open loop, with no controller. The start is not the
library's initial state: each quadrotor starts at rest, tilted 1.7° outward so that its
constant thrust (16.3 N) balances its weight and the pull of its cable: with the load level
and at rest and the cables at the library's initial directions (17° off the vertical),
this is a static equilibrium. The load and cables start with the library's initial
velocities, including a rotation of the load, so the system moves away from it. With $M_i = 0$ the attitude equations do not involve the load, and each vehicle
keeps its tilt. The load rocks up to 37° from level, the cables stay between 3° and 46°
from the vertical, the vehicles come within 0.49 m of each other (centre to centre), and
the load drifts 2.2 m sideways while it rises from 0.5 m to 2.2 m. Drag to orbit; the full report (equations,
standard form, parameters) opens in a
<a href="../model-reports/multi_quadrotor_rigid_load.html">separate page</a>.

<iframe class="model-report" src="../model-reports/multi_quadrotor_rigid_load.html?view" title="Quadrotors carrying a rigid body, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 320" role="img" aria-labelledby="mrl-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="mrl-title">Quadrotors carrying a rigid body: load at x_L with attitude R_L; cable i from quadrotor i to the attachment point x_L + R_L ρ_i, direction q_i, length l_i; quadrotor i with attitude R_i, thrust f_i and moment M_i</title>
  <defs>
    <marker id="mrl-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <line x1="30" y1="300" x2="30" y2="260" marker-end="url(#mrl-arr)"/>
    <line x1="30" y1="300" x2="70" y2="300" marker-end="url(#mrl-arr)"/>
  </g>
  <!-- the load: a flat box -->
  <path d="M 80,222 L 230,192 L 310,230 L 160,262 Z M 80,222 L 80,234 L 160,274 L 310,242 L 310,230 M 160,262 L 160,274" fill="none" stroke="currentColor" stroke-width="1.6"/>
  <g fill="none" stroke="currentColor" stroke-width="1.6">
    <line x1="80" y1="95" x2="125" y2="226"/>
    <line x1="310" y1="95" x2="265" y2="232"/>
    <line x1="210" y1="48" x2="205" y2="203"/>
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
  <line x1="80.0" y1="87.0" x2="80.0" y2="47.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#mrl-arr)"/>
  <path d="M 67,65.0 A 13 4 0 1 0 93,65.0" fill="none" stroke="#e08a2c" stroke-width="1.6" marker-end="url(#mrl-arr)"/>
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
  <line x1="310.0" y1="87.0" x2="310.0" y2="47.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#mrl-arr)"/>
  <!-- quadrotor at (210,48) -->
  <g fill="none" stroke="currentColor" stroke-width="1.6">
    <line x1="172" y1="48" x2="248" y2="48"/>
    <line x1="190" y1="58" x2="230" y2="38"/>
    <ellipse cx="172" cy="48" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="248" cy="48" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="230" cy="38" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="190" cy="58" rx="13" ry="4" stroke-width="1.2"/>
  </g>
  <rect x="201" y="43" width="18" height="10" rx="2" fill="currentColor"/>
  <line x1="210.0" y1="40.0" x2="210.0" y2="4.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#mrl-arr)"/>
  <line x1="80.0" y1="95.0" x2="94.6" y2="137.6" stroke="#3f8fa0" stroke-width="3" marker-end="url(#mrl-arr)"/>
  <line x1="310.0" y1="95.0" x2="296.0" y2="137.8" stroke="#3f8fa0" stroke-width="3" marker-end="url(#mrl-arr)"/>
  <line x1="210.0" y1="48.0" x2="208.5" y2="93.0" stroke="#3f8fa0" stroke-width="3" marker-end="url(#mrl-arr)"/>
  <line x1="195.0" y1="226.5" x2="127.0" y2="226.0" stroke="currentColor" stroke-width="1.4" marker-end="url(#mrl-arr)"/>
  <circle cx="125" cy="226" r="3.5" fill="currentColor"/>
  <circle cx="265" cy="232" r="3.5" fill="currentColor"/>
  <circle cx="205" cy="203" r="3.5" fill="currentColor"/>
  <circle cx="195" cy="226.5" r="4" fill="#3f8fa0"/>
  <line x1="195.0" y1="232.5" x2="195.0" y2="300.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#mrl-arr)"/>
  <circle cx="130" cy="173" r="7" fill="none" stroke="#3f8fa0" stroke-width="1.6"/>
  <circle cx="130" cy="173" r="2" fill="#3f8fa0"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="60" y="52" fill="#e08a2c">f<tspan font-size="12" dy="4">1</tspan></text>
    <text x="98" y="68" fill="#e08a2c">M<tspan font-size="12" dy="4">1</tspan></text>
    <text x="316" y="58" fill="#e08a2c">f<tspan font-size="12" dy="4">2</tspan></text>
    <text x="216" y="16" fill="#e08a2c">f<tspan font-size="12" dy="4">3</tspan></text>
    <text x="8" y="130" fill="#3f8fa0">R<tspan font-size="12" dy="4">1</tspan><tspan dy="-4">, Ω</tspan><tspan font-size="12" dy="4">1</tspan></text>
    <text x="8" y="154">m<tspan font-size="12" dy="4">1</tspan><tspan dy="-4">, J</tspan><tspan font-size="12" dy="4">1</tspan></text>
    <text x="141" y="178" fill="#3f8fa0">ω<tspan font-size="12" dy="4">1</tspan></text>
    <text x="24" y="252">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="76" y="305">e<tspan font-size="12" dy="4">1</tspan></text>
    <text x="106" y="134" fill="#3f8fa0">q<tspan font-size="12" dy="4">1</tspan></text>
    <text x="76" y="178">l<tspan font-size="12" dy="4">1</tspan></text>
    <text x="304" y="140" fill="#3f8fa0">q<tspan font-size="12" dy="4">2</tspan></text>
    <text x="217" y="96" fill="#3f8fa0">q<tspan font-size="12" dy="4">3</tspan></text>
    <text x="148" y="246">ρ<tspan font-size="12" dy="4">1</tspan></text>
    <text x="203" y="300" fill="#e08a2c">m<tspan font-size="12" dy="4">L</tspan><tspan dy="-4"> g</tspan></text>
    <text x="246" y="280" fill="#3f8fa0">x<tspan font-size="12" dy="4">L</tspan><tspan dy="-4">, R</tspan><tspan font-size="12" dy="4">L</tspan><tspan dy="-4">, Ω</tspan><tspan font-size="12" dy="4">R</tspan><tspan font-size="9" dy="2">L</tspan></text>
  </g>
</svg>
<figcaption>Drawn with three quadrotors (five on this page). Quadrotor 1 is labelled in
full; the others are the same with their own index.
Cable <span class="math notranslate nohighlight">\(i\)</span> is attached at <span class="math notranslate nohighlight">\(x_L + R_L \rho_i\)</span>
(<span class="math notranslate nohighlight">\(\rho_i\)</span> constant in the load frame, shown for the first), and quadrotor
<span class="math notranslate nohighlight">\(i\)</span> is at <span class="math notranslate nohighlight">\(x_L + R_L \rho_i - l_i q_i\)</span>. <span class="math notranslate nohighlight">\(\omega_1 \perp q_1\)</span> is out of the page here. Each
thrust <span class="math notranslate nohighlight">\(f_i\)</span> acts along <span class="math notranslate nohighlight">\(R_i e_3\)</span>; each vehicle also takes a
body-frame moment <span class="math notranslate nohighlight">\(M_i\)</span> (shown for the first). Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span> on
every body (shown on the load only).</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $x_L$ | centre of mass of the load | $\mathbb{R}^3$ |
| $\dot x_L$ | velocity of the centre of mass of the load | $\mathbb{R}^3$ |
| $\delta x_L$ | variation of $x_L$ | $\mathbb{R}^3$ |
| $R_L$ | attitude of the load, body frame to inertial frame | $SO(3)$ |
| $\Omega_{R_L}$ | body angular velocity of the load, $\dot R_L = R_L \widehat{\Omega_{R_L}}$ | $\mathbb{R}^3$ |
| $\eta_{R_L}$ | variation of $R_L$: $\delta R_L = R_L \widehat{\eta_{R_L}}$ | $\mathbb{R}^3$ |
| $q_i$ | direction of cable $i$, quadrotor $i$ to its attachment point | $S^2$ |
| $\omega_i$ | angular velocity of cable $i$, $\dot q_i = \omega_i \times q_i$, $q_i \cdot \omega_i = 0$ | $T_{q_i} S^2$ |
| $\xi_i$ | variation of $q_i$: $\delta q_i = \xi_i \times q_i$ | $T_{q_i} S^2$ |
| $R_i$ | attitude of quadrotor $i$, body frame to inertial frame | $SO(3)$ |
| $\Omega_i$ | body angular velocity of quadrotor $i$, $\dot R_i = R_i \widehat{\Omega_i}$ | $\mathbb{R}^3$ |
| $\eta_i$ | variation of $R_i$: $\delta R_i = R_i \widehat{\eta_i}$ | $\mathbb{R}^3$ |
| $f_i$ | thrust of quadrotor $i$ along $R_i e_3$ | $\mathbb{R}$ |
| $M_i$ | moment of quadrotor $i$, body frame | $\mathbb{R}^3$ |
| $m_L, J_L, \rho_i$ | load mass and inertia, attachment points in the load frame | constants |
| $m_i, J_i, l_i, g, e_3$ | quadrotor masses and inertias, cable lengths, gravity, vertical unit vector | constants |

## Lagrangian and equations for every $n$

<!-- model:indexed -->
Lagrangian for every $n$, as `geomech.indexed.systems.multi_quadrotor_rigid_load()` holds it (the velocities expanded, so the inner sums are the inertia couplings):

$$
\begin{aligned}
L = {}& \frac{1}{2} m_L \left(\dot{x_L} \cdot \dot{x_L}\right) \\
& + \frac{1}{2} \left(\Omega_{R_L} \cdot J_L \Omega_{R_L}\right) \\
& -m_L g \left(e_3 \cdot x_L\right) \\
& + \sum_{j=1}^{n} \frac{1}{2} m_{j} \left(\left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_{j} -l_{j} \left(\omega_{j} \times q_{j}\right)\right) \cdot \left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_{j} -l_{j} \left(\omega_{j} \times q_{j}\right)\right)\right) \\
& + \sum_{j=1}^{n} \frac{1}{2} \left(\Omega_{j} \cdot J_{j} \Omega_{j}\right) \\
& -\sum_{j=1}^{n} m_{j} g \left(e_3 \cdot \left(x_L + R_L \rho_{j} -l_{j} q_{j}\right)\right)
\end{aligned}
$$

Virtual work:

$$
\begin{aligned}
\delta W = {}& \sum_{j=1}^{n} \left(\delta x_L + R_L \left(\eta_{R_L} \times \rho_{j}\right) -l_{j} \left(\xi_{j} \times q_{j}\right)\right) \cdot f_{j} R_{j} e_3 \\
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
& -\sum_{j=1}^{n} m_{j} R_L \left(\Omega_{R_L} \times \left(\Omega_{R_L} \times \rho_{j}\right)\right) \\
& -\sum_{j=1}^{n} m_{j} R_L \left(\dot{\Omega_{R_L}} \times \rho_{j}\right) \\
& + \sum_{j=1}^{n} m_{j} l_{j} \left(\dot{\omega_{j}} \times q_{j}\right) \\
& -\sum_{j=1}^{n} \left(\omega_{j} \cdot \omega_{j}\right) m_{j} l_{j} q_{j} \\
& -\sum_{j=1}^{n} m_{j} g e_3 \\
& + \sum_{j=1}^{n} f_{j} R_{j} e_3
\end{aligned}
$$

Coefficient of $\eta_{R_L}$, the equation of $R_L$:

$$
\begin{aligned}
0 = {}& -J_L \dot{\Omega_{R_L}} \\
& -\Omega_{R_L} \times \left(J_L \Omega_{R_L}\right) \\
& + \sum_{j=1}^{n} m_{j} \left(\left(\Omega_{R_L} \times \rho_{j}\right) \times \left(R_L^T \dot{x_L}\right)\right) \\
& -\sum_{j=1}^{n} m_{j} l_{j} \left(\left(\Omega_{R_L} \times \rho_{j}\right) \times \left(R_L^T \left(\omega_{j} \times q_{j}\right)\right)\right) \\
& + \sum_{j=1}^{n} m_{j} \left(\rho_{j} \times \left(\Omega_{R_L} \times \left(R_L^T \dot{x_L}\right)\right)\right) \\
& -\sum_{j=1}^{n} m_{j} \left(\rho_{j} \times \left(R_L^T \ddot{x_L}\right)\right) \\
& -\sum_{j=1}^{n} m_{j} \left(\rho_{j} \times \left(\dot{\Omega_{R_L}} \times \rho_{j}\right)\right) \\
& -\sum_{j=1}^{n} m_{j} l_{j} \left(\rho_{j} \times \left(\Omega_{R_L} \times \left(R_L^T \left(\omega_{j} \times q_{j}\right)\right)\right)\right) \\
& + \sum_{j=1}^{n} m_{j} l_{j} \left(\rho_{j} \times \left(R_L^T \left(\dot{\omega_{j}} \times q_{j}\right)\right)\right) \\
& -\sum_{j=1}^{n} \left(\omega_{j} \cdot \omega_{j}\right) m_{j} l_{j} \left(\rho_{j} \times \left(R_L^T q_{j}\right)\right) \\
& -\sum_{j=1}^{n} m_{j} \left(\Omega_{R_L} \times \left(\rho_{j} \times \left(R_L^T \dot{x_L}\right)\right)\right) \\
& -\sum_{j=1}^{n} m_{j} \left(\Omega_{R_L} \times \left(\rho_{j} \times \left(\Omega_{R_L} \times \rho_{j}\right)\right)\right) \\
& + \sum_{j=1}^{n} m_{j} l_{j} \left(\Omega_{R_L} \times \left(\rho_{j} \times \left(R_L^T \left(\omega_{j} \times q_{j}\right)\right)\right)\right) \\
& -\sum_{j=1}^{n} m_{j} g \left(\rho_{j} \times \left(R_L^T e_3\right)\right) \\
& + \sum_{j=1}^{n} f_{j} \left(\rho_{j} \times \left(R_L^T R_{j} e_3\right)\right)
\end{aligned}
$$

Coefficient of $\xi_{i}$, the equation of $q_i$:

$$
\begin{aligned}
0 = {}& m_{i} l_{i} \left(\left(\omega_{i} \times q_{i}\right) \times \dot{x_L}\right) \\
& + m_{i} l_{i} \left(q_{i} \times \ddot{x_L}\right) \\
& + m_{i} l_{i} \left(\left(\omega_{i} \times q_{i}\right) \times \left(R_L \left(\Omega_{R_L} \times \rho_{i}\right)\right)\right) \\
& + m_{i} l_{i} \left(q_{i} \times \left(R_L \left(\Omega_{R_L} \times \left(\Omega_{R_L} \times \rho_{i}\right)\right)\right)\right) \\
& + m_{i} l_{i} \left(q_{i} \times \left(R_L \left(\dot{\Omega_{R_L}} \times \rho_{i}\right)\right)\right) \\
& -m_{i} l_{i} l_{i} \left(q_{i} \times \left(\dot{\omega_{i}} \times q_{i}\right)\right) \\
& -2 m_{i} l_{i} \left(\omega_{i} \times \left(q_{i} \times \dot{x_L}\right)\right) \\
& -2 m_{i} l_{i} \left(\omega_{i} \times \left(q_{i} \times \left(R_L \left(\Omega_{R_L} \times \rho_{i}\right)\right)\right)\right) \\
& + m_{i} l_{i} \left(q_{i} \times \left(\omega_{i} \times \dot{x_L}\right)\right) \\
& + m_{i} l_{i} \left(q_{i} \times \left(\omega_{i} \times \left(R_L \left(\Omega_{R_L} \times \rho_{i}\right)\right)\right)\right) \\
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
& + \frac{1}{2} \left(\Omega_{R_L} \cdot J_L \Omega_{R_L}\right) \\
& + \frac{1}{2} m_1 \left(\left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_1 -l_1 \left(\omega_1 \times q_1\right)\right) \cdot \left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_1 -l_1 \left(\omega_1 \times q_1\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_1 \cdot J_1 \Omega_1\right) \\
& + \frac{1}{2} m_2 \left(\left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_2 -l_2 \left(\omega_2 \times q_2\right)\right) \cdot \left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_2 -l_2 \left(\omega_2 \times q_2\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_2 \cdot J_2 \Omega_2\right) \\
& + \frac{1}{2} m_3 \left(\left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_3 -l_3 \left(\omega_3 \times q_3\right)\right) \cdot \left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_3 -l_3 \left(\omega_3 \times q_3\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_3 \cdot J_3 \Omega_3\right) \\
& + \frac{1}{2} m_4 \left(\left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_4 -l_4 \left(\omega_4 \times q_4\right)\right) \cdot \left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_4 -l_4 \left(\omega_4 \times q_4\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_4 \cdot J_4 \Omega_4\right) \\
& + \frac{1}{2} m_5 \left(\left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_5 -l_5 \left(\omega_5 \times q_5\right)\right) \cdot \left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_5 -l_5 \left(\omega_5 \times q_5\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_5 \cdot J_5 \Omega_5\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& m_L g \left(e_3 \cdot x_L\right) \\
& + m_1 g \left(e_3 \cdot \left(x_L + R_L \rho_1 -l_1 q_1\right)\right) \\
& + m_2 g \left(e_3 \cdot \left(x_L + R_L \rho_2 -l_2 q_2\right)\right) \\
& + m_3 g \left(e_3 \cdot \left(x_L + R_L \rho_3 -l_3 q_3\right)\right) \\
& + m_4 g \left(e_3 \cdot \left(x_L + R_L \rho_4 -l_4 q_4\right)\right) \\
& + m_5 g \left(e_3 \cdot \left(x_L + R_L \rho_5 -l_5 q_5\right)\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m_L \left(\dot{x_L} \cdot \dot{x_L}\right) \\
& + \frac{1}{2} \left(\Omega_{R_L} \cdot J_L \Omega_{R_L}\right) \\
& -m_L g \left(e_3 \cdot x_L\right) \\
& + \frac{1}{2} m_1 \left(\left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_1 -l_1 \left(\omega_1 \times q_1\right)\right) \cdot \left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_1 -l_1 \left(\omega_1 \times q_1\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_1 \cdot J_1 \Omega_1\right) \\
& -m_1 g \left(e_3 \cdot \left(x_L + R_L \rho_1 -l_1 q_1\right)\right) \\
& + \frac{1}{2} m_2 \left(\left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_2 -l_2 \left(\omega_2 \times q_2\right)\right) \cdot \left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_2 -l_2 \left(\omega_2 \times q_2\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_2 \cdot J_2 \Omega_2\right) \\
& -m_2 g \left(e_3 \cdot \left(x_L + R_L \rho_2 -l_2 q_2\right)\right) \\
& + \frac{1}{2} m_3 \left(\left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_3 -l_3 \left(\omega_3 \times q_3\right)\right) \cdot \left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_3 -l_3 \left(\omega_3 \times q_3\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_3 \cdot J_3 \Omega_3\right) \\
& -m_3 g \left(e_3 \cdot \left(x_L + R_L \rho_3 -l_3 q_3\right)\right) \\
& + \frac{1}{2} m_4 \left(\left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_4 -l_4 \left(\omega_4 \times q_4\right)\right) \cdot \left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_4 -l_4 \left(\omega_4 \times q_4\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_4 \cdot J_4 \Omega_4\right) \\
& -m_4 g \left(e_3 \cdot \left(x_L + R_L \rho_4 -l_4 q_4\right)\right) \\
& + \frac{1}{2} m_5 \left(\left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_5 -l_5 \left(\omega_5 \times q_5\right)\right) \cdot \left(\dot{x_L} + R_L \widehat{\Omega_{R_L}} \rho_5 -l_5 \left(\omega_5 \times q_5\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_5 \cdot J_5 \Omega_5\right) \\
& -m_5 g \left(e_3 \cdot \left(x_L + R_L \rho_5 -l_5 q_5\right)\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \left(\delta x_L + R_L \left(\eta_{R_L} \times \rho_1\right) -l_1 \left(\xi_1 \times q_1\right)\right) \cdot f_1 R_1 e_3 \\
& + \eta_1 \cdot M_1 \\
& + \left(\delta x_L + R_L \left(\eta_{R_L} \times \rho_2\right) -l_2 \left(\xi_2 \times q_2\right)\right) \cdot f_2 R_2 e_3 \\
& + \eta_2 \cdot M_2 \\
& + \left(\delta x_L + R_L \left(\eta_{R_L} \times \rho_3\right) -l_3 \left(\xi_3 \times q_3\right)\right) \cdot f_3 R_3 e_3 \\
& + \eta_3 \cdot M_3 \\
& + \left(\delta x_L + R_L \left(\eta_{R_L} \times \rho_4\right) -l_4 \left(\xi_4 \times q_4\right)\right) \cdot f_4 R_4 e_3 \\
& + \eta_4 \cdot M_4 \\
& + \left(\delta x_L + R_L \left(\eta_{R_L} \times \rho_5\right) -l_5 \left(\xi_5 \times q_5\right)\right) \cdot f_5 R_5 e_3 \\
& + \eta_5 \cdot M_5
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion at $n = 5$

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); the $R_L$ and $R_i$ equations are in the body frame
and each $q_i$ equation holds on the tangent plane at $q_i$.

<!-- model:eom -->
<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\delta x_L\)</span>, the equation of <span class="math notranslate nohighlight">\(x_L\)</span> (37 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_L \ddot{x_L} \\
& -m_L g e_3 \\
& -m_1 \ddot{x_L} \\
& -m_1 R_L \left(\Omega_{R_L} \times \left(\Omega_{R_L} \times \rho_1\right)\right) \\
& -m_1 R_L \left(\dot{\Omega_{R_L}} \times \rho_1\right) \\
& + m_1 l_1 \left(\dot{\omega_1} \times q_1\right) \\
& -\left(\omega_1 \cdot \omega_1\right) m_1 l_1 q_1 \\
& -m_1 g e_3 \\
& -m_2 \ddot{x_L} \\
& -m_2 R_L \left(\Omega_{R_L} \times \left(\Omega_{R_L} \times \rho_2\right)\right) \\
& -m_2 R_L \left(\dot{\Omega_{R_L}} \times \rho_2\right) \\
& + m_2 l_2 \left(\dot{\omega_2} \times q_2\right) \\
& -\left(\omega_2 \cdot \omega_2\right) m_2 l_2 q_2 \\
& -m_2 g e_3 \\
& -m_3 \ddot{x_L} \\
& -m_3 R_L \left(\Omega_{R_L} \times \left(\Omega_{R_L} \times \rho_3\right)\right) \\
& -m_3 R_L \left(\dot{\Omega_{R_L}} \times \rho_3\right) \\
& + m_3 l_3 \left(\dot{\omega_3} \times q_3\right) \\
& -\left(\omega_3 \cdot \omega_3\right) m_3 l_3 q_3 \\
& -m_3 g e_3 \\
& -m_4 \ddot{x_L} \\
& -m_4 R_L \left(\Omega_{R_L} \times \left(\Omega_{R_L} \times \rho_4\right)\right) \\
& -m_4 R_L \left(\dot{\Omega_{R_L}} \times \rho_4\right) \\
& + m_4 l_4 \left(\dot{\omega_4} \times q_4\right) \\
& -\left(\omega_4 \cdot \omega_4\right) m_4 l_4 q_4 \\
& -m_4 g e_3 \\
& -m_5 \ddot{x_L} \\
& -m_5 R_L \left(\Omega_{R_L} \times \left(\Omega_{R_L} \times \rho_5\right)\right) \\
& -m_5 R_L \left(\dot{\Omega_{R_L}} \times \rho_5\right) \\
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
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_1\)</span>, the equation of <span class="math notranslate nohighlight">\(q_1\)</span> (12 terms)</summary>

$$
\begin{aligned}
0 = {}& m_1 l_1 \left(\left(\omega_1 \times q_1\right) \times \dot{x_L}\right) \\
& + m_1 l_1 \left(q_1 \times \ddot{x_L}\right) \\
& + m_1 l_1 \left(\left(\omega_1 \times q_1\right) \times \left(R_L \left(\Omega_{R_L} \times \rho_1\right)\right)\right) \\
& + m_1 l_1 \left(q_1 \times \left(R_L \left(\Omega_{R_L} \times \left(\Omega_{R_L} \times \rho_1\right)\right)\right)\right) \\
& + m_1 l_1 \left(q_1 \times \left(R_L \left(\dot{\Omega_{R_L}} \times \rho_1\right)\right)\right) \\
& -m_1 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -2 m_1 l_1 \left(\omega_1 \times \left(q_1 \times \dot{x_L}\right)\right) \\
& -2 m_1 l_1 \left(\omega_1 \times \left(q_1 \times \left(R_L \left(\Omega_{R_L} \times \rho_1\right)\right)\right)\right) \\
& + m_1 l_1 \left(q_1 \times \left(\omega_1 \times \dot{x_L}\right)\right) \\
& + m_1 l_1 \left(q_1 \times \left(\omega_1 \times \left(R_L \left(\Omega_{R_L} \times \rho_1\right)\right)\right)\right) \\
& + m_1 g l_1 \left(q_1 \times e_3\right) \\
& -l_1 f_1 \left(q_1 \times \left(R_1 e_3\right)\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_2\)</span>, the equation of <span class="math notranslate nohighlight">\(q_2\)</span> (12 terms)</summary>

$$
\begin{aligned}
0 = {}& m_2 l_2 \left(\left(\omega_2 \times q_2\right) \times \dot{x_L}\right) \\
& + m_2 l_2 \left(q_2 \times \ddot{x_L}\right) \\
& + m_2 l_2 \left(\left(\omega_2 \times q_2\right) \times \left(R_L \left(\Omega_{R_L} \times \rho_2\right)\right)\right) \\
& + m_2 l_2 \left(q_2 \times \left(R_L \left(\Omega_{R_L} \times \left(\Omega_{R_L} \times \rho_2\right)\right)\right)\right) \\
& + m_2 l_2 \left(q_2 \times \left(R_L \left(\dot{\Omega_{R_L}} \times \rho_2\right)\right)\right) \\
& -m_2 l_2 l_2 \left(q_2 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& -2 m_2 l_2 \left(\omega_2 \times \left(q_2 \times \dot{x_L}\right)\right) \\
& -2 m_2 l_2 \left(\omega_2 \times \left(q_2 \times \left(R_L \left(\Omega_{R_L} \times \rho_2\right)\right)\right)\right) \\
& + m_2 l_2 \left(q_2 \times \left(\omega_2 \times \dot{x_L}\right)\right) \\
& + m_2 l_2 \left(q_2 \times \left(\omega_2 \times \left(R_L \left(\Omega_{R_L} \times \rho_2\right)\right)\right)\right) \\
& + m_2 g l_2 \left(q_2 \times e_3\right) \\
& -l_2 f_2 \left(q_2 \times \left(R_2 e_3\right)\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_3\)</span>, the equation of <span class="math notranslate nohighlight">\(q_3\)</span> (12 terms)</summary>

$$
\begin{aligned}
0 = {}& m_3 l_3 \left(\left(\omega_3 \times q_3\right) \times \dot{x_L}\right) \\
& + m_3 l_3 \left(q_3 \times \ddot{x_L}\right) \\
& + m_3 l_3 \left(\left(\omega_3 \times q_3\right) \times \left(R_L \left(\Omega_{R_L} \times \rho_3\right)\right)\right) \\
& + m_3 l_3 \left(q_3 \times \left(R_L \left(\Omega_{R_L} \times \left(\Omega_{R_L} \times \rho_3\right)\right)\right)\right) \\
& + m_3 l_3 \left(q_3 \times \left(R_L \left(\dot{\Omega_{R_L}} \times \rho_3\right)\right)\right) \\
& -m_3 l_3 l_3 \left(q_3 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& -2 m_3 l_3 \left(\omega_3 \times \left(q_3 \times \dot{x_L}\right)\right) \\
& -2 m_3 l_3 \left(\omega_3 \times \left(q_3 \times \left(R_L \left(\Omega_{R_L} \times \rho_3\right)\right)\right)\right) \\
& + m_3 l_3 \left(q_3 \times \left(\omega_3 \times \dot{x_L}\right)\right) \\
& + m_3 l_3 \left(q_3 \times \left(\omega_3 \times \left(R_L \left(\Omega_{R_L} \times \rho_3\right)\right)\right)\right) \\
& + m_3 g l_3 \left(q_3 \times e_3\right) \\
& -l_3 f_3 \left(q_3 \times \left(R_3 e_3\right)\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_4\)</span>, the equation of <span class="math notranslate nohighlight">\(q_4\)</span> (12 terms)</summary>

$$
\begin{aligned}
0 = {}& m_4 l_4 \left(\left(\omega_4 \times q_4\right) \times \dot{x_L}\right) \\
& + m_4 l_4 \left(q_4 \times \ddot{x_L}\right) \\
& + m_4 l_4 \left(\left(\omega_4 \times q_4\right) \times \left(R_L \left(\Omega_{R_L} \times \rho_4\right)\right)\right) \\
& + m_4 l_4 \left(q_4 \times \left(R_L \left(\Omega_{R_L} \times \left(\Omega_{R_L} \times \rho_4\right)\right)\right)\right) \\
& + m_4 l_4 \left(q_4 \times \left(R_L \left(\dot{\Omega_{R_L}} \times \rho_4\right)\right)\right) \\
& -m_4 l_4 l_4 \left(q_4 \times \left(\dot{\omega_4} \times q_4\right)\right) \\
& -2 m_4 l_4 \left(\omega_4 \times \left(q_4 \times \dot{x_L}\right)\right) \\
& -2 m_4 l_4 \left(\omega_4 \times \left(q_4 \times \left(R_L \left(\Omega_{R_L} \times \rho_4\right)\right)\right)\right) \\
& + m_4 l_4 \left(q_4 \times \left(\omega_4 \times \dot{x_L}\right)\right) \\
& + m_4 l_4 \left(q_4 \times \left(\omega_4 \times \left(R_L \left(\Omega_{R_L} \times \rho_4\right)\right)\right)\right) \\
& + m_4 g l_4 \left(q_4 \times e_3\right) \\
& -l_4 f_4 \left(q_4 \times \left(R_4 e_3\right)\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_5\)</span>, the equation of <span class="math notranslate nohighlight">\(q_5\)</span> (12 terms)</summary>

$$
\begin{aligned}
0 = {}& m_5 l_5 \left(\left(\omega_5 \times q_5\right) \times \dot{x_L}\right) \\
& + m_5 l_5 \left(q_5 \times \ddot{x_L}\right) \\
& + m_5 l_5 \left(\left(\omega_5 \times q_5\right) \times \left(R_L \left(\Omega_{R_L} \times \rho_5\right)\right)\right) \\
& + m_5 l_5 \left(q_5 \times \left(R_L \left(\Omega_{R_L} \times \left(\Omega_{R_L} \times \rho_5\right)\right)\right)\right) \\
& + m_5 l_5 \left(q_5 \times \left(R_L \left(\dot{\Omega_{R_L}} \times \rho_5\right)\right)\right) \\
& -m_5 l_5 l_5 \left(q_5 \times \left(\dot{\omega_5} \times q_5\right)\right) \\
& -2 m_5 l_5 \left(\omega_5 \times \left(q_5 \times \dot{x_L}\right)\right) \\
& -2 m_5 l_5 \left(\omega_5 \times \left(q_5 \times \left(R_L \left(\Omega_{R_L} \times \rho_5\right)\right)\right)\right) \\
& + m_5 l_5 \left(q_5 \times \left(\omega_5 \times \dot{x_L}\right)\right) \\
& + m_5 l_5 \left(q_5 \times \left(\omega_5 \times \left(R_L \left(\Omega_{R_L} \times \rho_5\right)\right)\right)\right) \\
& + m_5 g l_5 \left(q_5 \times e_3\right) \\
& -l_5 f_5 \left(q_5 \times \left(R_5 e_3\right)\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_{R_L}\)</span>, the equation of <span class="math notranslate nohighlight">\(R_L\)</span> (67 terms)</summary>

$$
\begin{aligned}
0 = {}& -J_L \dot{\Omega_{R_L}} \\
& -\Omega_{R_L} \times \left(J_L \Omega_{R_L}\right) \\
& + m_1 \left(\left(\Omega_{R_L} \times \rho_1\right) \times \left(R_L^T \dot{x_L}\right)\right) \\
& -m_1 l_1 \left(\left(\Omega_{R_L} \times \rho_1\right) \times \left(R_L^T \left(\omega_1 \times q_1\right)\right)\right) \\
& + m_1 \left(\rho_1 \times \left(\Omega_{R_L} \times \left(R_L^T \dot{x_L}\right)\right)\right) \\
& -m_1 \left(\rho_1 \times \left(R_L^T \ddot{x_L}\right)\right) \\
& -m_1 \left(\rho_1 \times \left(\dot{\Omega_{R_L}} \times \rho_1\right)\right) \\
& -m_1 l_1 \left(\rho_1 \times \left(\Omega_{R_L} \times \left(R_L^T \left(\omega_1 \times q_1\right)\right)\right)\right) \\
& + m_1 l_1 \left(\rho_1 \times \left(R_L^T \left(\dot{\omega_1} \times q_1\right)\right)\right) \\
& -\left(\omega_1 \cdot \omega_1\right) m_1 l_1 \left(\rho_1 \times \left(R_L^T q_1\right)\right) \\
& -m_1 \left(\Omega_{R_L} \times \left(\rho_1 \times \left(R_L^T \dot{x_L}\right)\right)\right) \\
& -m_1 \left(\Omega_{R_L} \times \left(\rho_1 \times \left(\Omega_{R_L} \times \rho_1\right)\right)\right) \\
& + m_1 l_1 \left(\Omega_{R_L} \times \left(\rho_1 \times \left(R_L^T \left(\omega_1 \times q_1\right)\right)\right)\right) \\
& -m_1 g \left(\rho_1 \times \left(R_L^T e_3\right)\right) \\
& + m_2 \left(\left(\Omega_{R_L} \times \rho_2\right) \times \left(R_L^T \dot{x_L}\right)\right) \\
& -m_2 l_2 \left(\left(\Omega_{R_L} \times \rho_2\right) \times \left(R_L^T \left(\omega_2 \times q_2\right)\right)\right) \\
& + m_2 \left(\rho_2 \times \left(\Omega_{R_L} \times \left(R_L^T \dot{x_L}\right)\right)\right) \\
& -m_2 \left(\rho_2 \times \left(R_L^T \ddot{x_L}\right)\right) \\
& -m_2 \left(\rho_2 \times \left(\dot{\Omega_{R_L}} \times \rho_2\right)\right) \\
& -m_2 l_2 \left(\rho_2 \times \left(\Omega_{R_L} \times \left(R_L^T \left(\omega_2 \times q_2\right)\right)\right)\right) \\
& + m_2 l_2 \left(\rho_2 \times \left(R_L^T \left(\dot{\omega_2} \times q_2\right)\right)\right) \\
& -\left(\omega_2 \cdot \omega_2\right) m_2 l_2 \left(\rho_2 \times \left(R_L^T q_2\right)\right) \\
& -m_2 \left(\Omega_{R_L} \times \left(\rho_2 \times \left(R_L^T \dot{x_L}\right)\right)\right) \\
& -m_2 \left(\Omega_{R_L} \times \left(\rho_2 \times \left(\Omega_{R_L} \times \rho_2\right)\right)\right) \\
& + m_2 l_2 \left(\Omega_{R_L} \times \left(\rho_2 \times \left(R_L^T \left(\omega_2 \times q_2\right)\right)\right)\right) \\
& -m_2 g \left(\rho_2 \times \left(R_L^T e_3\right)\right) \\
& + m_3 \left(\left(\Omega_{R_L} \times \rho_3\right) \times \left(R_L^T \dot{x_L}\right)\right) \\
& -m_3 l_3 \left(\left(\Omega_{R_L} \times \rho_3\right) \times \left(R_L^T \left(\omega_3 \times q_3\right)\right)\right) \\
& + m_3 \left(\rho_3 \times \left(\Omega_{R_L} \times \left(R_L^T \dot{x_L}\right)\right)\right) \\
& -m_3 \left(\rho_3 \times \left(R_L^T \ddot{x_L}\right)\right) \\
& -m_3 \left(\rho_3 \times \left(\dot{\Omega_{R_L}} \times \rho_3\right)\right) \\
& -m_3 l_3 \left(\rho_3 \times \left(\Omega_{R_L} \times \left(R_L^T \left(\omega_3 \times q_3\right)\right)\right)\right) \\
& + m_3 l_3 \left(\rho_3 \times \left(R_L^T \left(\dot{\omega_3} \times q_3\right)\right)\right) \\
& -\left(\omega_3 \cdot \omega_3\right) m_3 l_3 \left(\rho_3 \times \left(R_L^T q_3\right)\right) \\
& -m_3 \left(\Omega_{R_L} \times \left(\rho_3 \times \left(R_L^T \dot{x_L}\right)\right)\right) \\
& -m_3 \left(\Omega_{R_L} \times \left(\rho_3 \times \left(\Omega_{R_L} \times \rho_3\right)\right)\right) \\
& + m_3 l_3 \left(\Omega_{R_L} \times \left(\rho_3 \times \left(R_L^T \left(\omega_3 \times q_3\right)\right)\right)\right) \\
& -m_3 g \left(\rho_3 \times \left(R_L^T e_3\right)\right) \\
& + m_4 \left(\left(\Omega_{R_L} \times \rho_4\right) \times \left(R_L^T \dot{x_L}\right)\right) \\
& -m_4 l_4 \left(\left(\Omega_{R_L} \times \rho_4\right) \times \left(R_L^T \left(\omega_4 \times q_4\right)\right)\right) \\
& + m_4 \left(\rho_4 \times \left(\Omega_{R_L} \times \left(R_L^T \dot{x_L}\right)\right)\right) \\
& -m_4 \left(\rho_4 \times \left(R_L^T \ddot{x_L}\right)\right) \\
& -m_4 \left(\rho_4 \times \left(\dot{\Omega_{R_L}} \times \rho_4\right)\right) \\
& -m_4 l_4 \left(\rho_4 \times \left(\Omega_{R_L} \times \left(R_L^T \left(\omega_4 \times q_4\right)\right)\right)\right) \\
& + m_4 l_4 \left(\rho_4 \times \left(R_L^T \left(\dot{\omega_4} \times q_4\right)\right)\right) \\
& -\left(\omega_4 \cdot \omega_4\right) m_4 l_4 \left(\rho_4 \times \left(R_L^T q_4\right)\right) \\
& -m_4 \left(\Omega_{R_L} \times \left(\rho_4 \times \left(R_L^T \dot{x_L}\right)\right)\right) \\
& -m_4 \left(\Omega_{R_L} \times \left(\rho_4 \times \left(\Omega_{R_L} \times \rho_4\right)\right)\right) \\
& + m_4 l_4 \left(\Omega_{R_L} \times \left(\rho_4 \times \left(R_L^T \left(\omega_4 \times q_4\right)\right)\right)\right) \\
& -m_4 g \left(\rho_4 \times \left(R_L^T e_3\right)\right) \\
& + m_5 \left(\left(\Omega_{R_L} \times \rho_5\right) \times \left(R_L^T \dot{x_L}\right)\right) \\
& -m_5 l_5 \left(\left(\Omega_{R_L} \times \rho_5\right) \times \left(R_L^T \left(\omega_5 \times q_5\right)\right)\right) \\
& + m_5 \left(\rho_5 \times \left(\Omega_{R_L} \times \left(R_L^T \dot{x_L}\right)\right)\right) \\
& -m_5 \left(\rho_5 \times \left(R_L^T \ddot{x_L}\right)\right) \\
& -m_5 \left(\rho_5 \times \left(\dot{\Omega_{R_L}} \times \rho_5\right)\right) \\
& -m_5 l_5 \left(\rho_5 \times \left(\Omega_{R_L} \times \left(R_L^T \left(\omega_5 \times q_5\right)\right)\right)\right) \\
& + m_5 l_5 \left(\rho_5 \times \left(R_L^T \left(\dot{\omega_5} \times q_5\right)\right)\right) \\
& -\left(\omega_5 \cdot \omega_5\right) m_5 l_5 \left(\rho_5 \times \left(R_L^T q_5\right)\right) \\
& -m_5 \left(\Omega_{R_L} \times \left(\rho_5 \times \left(R_L^T \dot{x_L}\right)\right)\right) \\
& -m_5 \left(\Omega_{R_L} \times \left(\rho_5 \times \left(\Omega_{R_L} \times \rho_5\right)\right)\right) \\
& + m_5 l_5 \left(\Omega_{R_L} \times \left(\rho_5 \times \left(R_L^T \left(\omega_5 \times q_5\right)\right)\right)\right) \\
& -m_5 g \left(\rho_5 \times \left(R_L^T e_3\right)\right) \\
& + f_1 \left(\rho_1 \times \left(R_L^T R_1 e_3\right)\right) \\
& + f_2 \left(\rho_2 \times \left(R_L^T R_2 e_3\right)\right) \\
& + f_3 \left(\rho_3 \times \left(R_L^T R_3 e_3\right)\right) \\
& + f_4 \left(\rho_4 \times \left(R_L^T R_4 e_3\right)\right) \\
& + f_5 \left(\rho_5 \times \left(R_L^T R_5 e_3\right)\right)
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

The model for every $n$, as `geomech.indexed.systems.multi_quadrotor_rigid_load()` writes
it, and the library system for $n = 5$:

```python
from geomech import S2, SO3, Cross, Dot, Matrix, Scalar, Vector
from geomech.core.base.constants import HALF
from geomech.indexed import Family, IndexedVariables, Sum, compute_indexed_eom, indices, n
from geomech.systems import get

(i,) = indices("i")
mL, g = Scalar("m_L", attr=["Constant"]), Scalar("g", attr=["Constant"])
JL = Matrix("J_L", attr=["Constant", "SymmetricMatrix"])
e3 = Vector("e3", attr=["Constant"])
xL, RL = Vector("x_L"), SO3("R_L")
OmL, etaL = RL.get_tangent_vector(), RL.get_variation_vector()
rho = Family(Vector, "rho", ["Constant"])  # attachment points, load frame
m, l = Family(Scalar, "m", ["Constant"]), Family(Scalar, "l", ["Constant"])
J = Family(Matrix, "J", attr=["Constant", "SymmetricMatrix"])
R, q = Family(SO3, "R"), Family(S2, "q")
f, M = Family(Scalar, "f"), Family(Vector, "M")
vL = xL.t_diff()
xQ = xL + RL * rho[i] - q[i] * l[i]  # quadrotor i
Om = R[i].get_tangent_vector()

L = HALF * mL * Dot(vL, vL) + HALF * Dot(OmL, JL * OmL) - mL * g * Dot(e3, xL)
L = L + Sum(
    i,
    1,
    n,
    HALF * m[i] * Dot(xQ.t_diff(), xQ.t_diff())
    + HALF * Dot(Om, J[i] * Om)
    - m[i] * g * Dot(e3, xQ),
)
dxQ = xL.delta() + RL * Cross(etaL, rho[i]) - Cross(q[i].get_variation_vector(), q[i]) * l[i]
dW = Sum(i, 1, n, Dot(dxQ, (R[i] * e3) * f[i]) + Dot(R[i].get_variation_vector(), M[i]))
variables = IndexedVariables(families=[q, R], vectors=[xL], matrices=[RL])
eom = compute_indexed_eom(L, dW, variables)  # every n

system = get("multi_quadrotor_rigid_load", n_quads=5)  # n from 1 to 20
eom_5 = system.derive()
```

The same model at $n = 5$ is `geomech.systems.get("multi_quadrotor_rigid_load", n_quads=5)`,
with default parameters, a default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-5}$ |
| Energy $H = T + V$, 6 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-8}$, with $H(0) = 112$ |
| Equations for every $n$ above, at $n = 5$: `check_eom` against the Euler–Lagrange evaluation of the generic $L$, $\delta W$ at $n = 5$ | below $10^{-5}$ |
| Generic $L$ at $n = 5$ against the library system's $L$, at random states | relative difference below $10^{-14}$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
