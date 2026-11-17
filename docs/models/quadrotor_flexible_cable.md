# Quadrotor with a flexible cable

A quadrotor carrying a load on a flexible cable, modelled as $n$ rigid massless links with
a point mass at the far end of each; the last mass is the load. Link $i$ has direction
$q_i \in S^2$, pointing away from the quadrotor, and length $l_i$, so mass $i$ is at
$x_Q + \sum_{j \le i} l_j q_j$. The configuration manifold is
$\mathbb{R}^3 \times SO(3) \times (S^2)^n$. The equations below are derived once for
every $n$ with `geomech.indexed`, then for $n = 5$ with the library system, which the
animation uses. Not taken from the literature; the equations are verified numerically
(`check_eom`). Library system `quadrotor_flexible_cable`.

## Animation

Six seconds with constant thrust equal to the total weight and zero moments: open loop,
with no controller. The start is not the library's initial state: the quadrotor starts
level and not rotating, while the five links keep the library's initial state, 23° to 48°
off the vertical. With
$M = 0$ and the thrust through the centre of mass, the attitude equation does not involve
the load, so a vehicle that starts level and at rest stays exactly level. The cable swings
as a whole, each link up to 48° from the vertical, and bends little: the first and last
links stay within 26° of each other. Nothing holds position: the quadrotor drifts 1.2 m
sideways and rises from 2 m to 3.7 m. Drag to orbit; the full report (equations, standard
form, parameters) opens in a
<a href="../model-reports/quadrotor_flexible_cable.html">separate page</a>.

<iframe class="model-report" src="../model-reports/quadrotor_flexible_cable.html?view" title="Quadrotor with a flexible cable, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 320" role="img" aria-labelledby="qfc-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="qfc-title">Quadrotor with a flexible cable: quadrotor at x_Q with attitude R, thrust f and moment M; three links with directions q_1, q_2, q_3 from the quadrotor end, lengths l_1, l_2, l_3 and point masses m_1, m_2 and the load m_3</title>
  <defs>
    <marker id="qfc-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <line x1="40" y1="300" x2="40" y2="260" marker-end="url(#qfc-arr)"/>
    <line x1="40" y1="300" x2="80" y2="300" marker-end="url(#qfc-arr)"/>
  </g>
  <g fill="none" stroke="currentColor" stroke-width="1.6">
    <line x1="120" y1="66" x2="150" y2="136"/>
    <line x1="150" y1="136" x2="195" y2="201"/>
    <line x1="195" y1="201" x2="255" y2="256"/>
  </g>
  <!-- quadrotor at (120,66) -->
  <g fill="none" stroke="currentColor" stroke-width="1.6">
    <line x1="70" y1="66" x2="170" y2="66"/>
    <line x1="92" y1="80" x2="148" y2="52"/>
    <ellipse cx="70" cy="66" rx="17" ry="5" stroke-width="1.2"/>
    <ellipse cx="170" cy="66" rx="17" ry="5" stroke-width="1.2"/>
    <ellipse cx="148" cy="52" rx="17" ry="5" stroke-width="1.2"/>
    <ellipse cx="92" cy="80" rx="17" ry="5" stroke-width="1.2"/>
  </g>
  <rect x="111" y="61" width="18" height="10" rx="2" fill="currentColor"/>
  <line x1="120.0" y1="58.0" x2="120.0" y2="14.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#qfc-arr)"/>
  <path d="M 107,33.8 A 13 5 0 1 0 133,33.8" fill="none" stroke="#e08a2c" stroke-width="1.6" marker-end="url(#qfc-arr)"/>
  <line x1="120.0" y1="66.0" x2="133.8" y2="98.2" stroke="#3f8fa0" stroke-width="3" marker-end="url(#qfc-arr)"/>
  <line x1="150.0" y1="136.0" x2="169.9" y2="164.8" stroke="#3f8fa0" stroke-width="3" marker-end="url(#qfc-arr)"/>
  <line x1="195.0" y1="201.0" x2="220.8" y2="224.7" stroke="#3f8fa0" stroke-width="3" marker-end="url(#qfc-arr)"/>
  <line x1="255.0" y1="267.0" x2="255.0" y2="300.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#qfc-arr)"/>
  <circle cx="160" cy="118" r="7" fill="none" stroke="#3f8fa0" stroke-width="1.6"/>
  <circle cx="160" cy="118" r="2" fill="#3f8fa0"/>
  <circle cx="150" cy="136" r="5" fill="currentColor"/>
  <circle cx="195" cy="201" r="5" fill="currentColor"/>
  <circle cx="255" cy="256" r="11" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="127" y="22" fill="#e08a2c">f</text>
    <text x="140" y="36" fill="#e08a2c">M</text>
    <text x="196" y="72" fill="#3f8fa0">R, Ω<tspan font-size="12" dy="4">R</tspan></text>
    <text x="32" y="72" fill="#3f8fa0">x<tspan font-size="12" dy="4">Q</tspan></text>
    <text x="48" y="106">m<tspan font-size="12" dy="4">Q</tspan><tspan dy="-4">, J</tspan></text>
    <text x="145" y="100" fill="#3f8fa0">q<tspan font-size="12" dy="4">1</tspan></text>
    <text x="171" y="123" fill="#3f8fa0">ω<tspan font-size="12" dy="4">1</tspan></text>
    <text x="162" y="142">m<tspan font-size="12" dy="4">1</tspan></text>
    <text x="183" y="162" fill="#3f8fa0">q<tspan font-size="12" dy="4">2</tspan></text>
    <text x="207" y="206">m<tspan font-size="12" dy="4">2</tspan></text>
    <text x="236" y="223" fill="#3f8fa0">q<tspan font-size="12" dy="4">3</tspan></text>
    <text x="272" y="261">m<tspan font-size="12" dy="4">3</tspan></text>
    <text x="263" y="302" fill="#e08a2c">m<tspan font-size="12" dy="4">3</tspan><tspan dy="-4"> g</tspan></text>
    <text x="110" y="114">l<tspan font-size="12" dy="4">1</tspan></text>
    <text x="150" y="186">l<tspan font-size="12" dy="4">2</tspan></text>
    <text x="204" y="248">l<tspan font-size="12" dy="4">3</tspan></text>
    <text x="34" y="252">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="86" y="305">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>Each <span class="math notranslate nohighlight">\(q_i\)</span> points from the quadrotor end of link <span class="math notranslate nohighlight">\(i\)</span> to its far end,
where the mass <span class="math notranslate nohighlight">\(m_i\)</span> sits; the last mass is the load (<span class="math notranslate nohighlight">\(m_3\)</span> in this
drawing of three links; <span class="math notranslate nohighlight">\(m_5\)</span> on this page). Each link has an angular
velocity <span class="math notranslate nohighlight">\(\omega_i \perp q_i\)</span> (shown for the first, out of the page
here). The thrust <span class="math notranslate nohighlight">\(f\)</span> acts along <span class="math notranslate nohighlight">\(R e_3\)</span> at <span class="math notranslate nohighlight">\(x_Q\)</span>. Gravity acts
along <span class="math notranslate nohighlight">\(-e_3\)</span> on every mass (shown on the load only).</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $x_Q$ | position of the quadrotor | $\mathbb{R}^3$ |
| $\dot x_Q$ | velocity of the quadrotor | $\mathbb{R}^3$ |
| $\delta x_Q$ | variation of $x_Q$ | $\mathbb{R}^3$ |
| $q_i$ | direction of link $i$, away from the quadrotor | $S^2$ |
| $\omega_i$ | angular velocity of link $i$, $\dot q_i = \omega_i \times q_i$, $q_i \cdot \omega_i = 0$ | $T_{q_i} S^2$ |
| $\xi_i$ | variation of $q_i$: $\delta q_i = \xi_i \times q_i$ | $T_{q_i} S^2$ |
| $R$ | attitude of the quadrotor, body frame to inertial frame | $SO(3)$ |
| $\Omega_{R}$ | body angular velocity, $\dot R = R \widehat{\Omega_{R}}$ | $\mathbb{R}^3$ |
| $\eta_{R}$ | variation of $R$: $\delta R = R \widehat{\eta_{R}}$ | $\mathbb{R}^3$ |
| $f$ | thrust along $R e_3$ | $\mathbb{R}$ |
| $M$ | moment, body frame | $\mathbb{R}^3$ |
| $m_Q, J, m_i, l_i, g, e_3$ | quadrotor mass and inertia, link masses (load $m_n$; $m_5$ on this page), link lengths, gravity, vertical unit vector | constants |

## Lagrangian and equations for every $n$

<!-- model:indexed -->
Lagrangian for every $n$, as `geomech.indexed.systems.quadrotor_flexible_cable()` holds it (the velocities expanded, so the inner sums are the inertia couplings):

$$
\begin{aligned}
L = {}& \frac{1}{2} m_Q \left(\dot{x_Q} \cdot \dot{x_Q}\right) \\
& + \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right) \\
& -m_Q g \left(e_3 \cdot x_Q\right) \\
& + \sum_{j=1}^{n} \frac{1}{2} m_{j} \left(\dot{x_Q} \cdot \dot{x_Q}\right) \\
& + \sum_{j=1}^{n} \sum_{k=j}^{n} \frac{1}{2} m_{k} \left(\dot{x_Q} \cdot l_{j} \left(\omega_{j} \times q_{j}\right)\right) \\
& + \sum_{j=1}^{n} \sum_{k=j}^{n} \frac{1}{2} m_{k} \left(l_{j} \left(\omega_{j} \times q_{j}\right) \cdot \dot{x_Q}\right) \\
& + \sum_{j=1}^{n} \sum_{k=1}^{n} \sum_{r=\max(j, k)}^{n} \frac{1}{2} m_{r} \left(l_{j} \left(\omega_{j} \times q_{j}\right) \cdot l_{k} \left(\omega_{k} \times q_{k}\right)\right) \\
& -\sum_{j=1}^{n} m_{j} g \left(x_Q \cdot e_3\right) \\
& -\sum_{j=1}^{n} \sum_{k=j}^{n} m_{k} g \left(l_{j} q_{j} \cdot e_3\right)
\end{aligned}
$$

Virtual work:

$$
\begin{aligned}
\delta W = {}& \delta x_Q \cdot f R e_3 \\
& + \eta_{R} \cdot M
\end{aligned}
$$

Equations from `compute_indexed_eom`, derived once for every $n$; $i$ is the body, $1 \le i \le n$:

Coefficient of $\delta x_Q$, the equation of $x_Q$:

$$
\begin{aligned}
0 = {}& -m_Q \ddot{x_Q} \\
& -m_Q g e_3 \\
& + f R e_3 \\
& -\sum_{j=1}^{n} m_{j} \ddot{x_Q} \\
& -\sum_{j=1}^{n} m_{j} g e_3 \\
& -\sum_{j=1}^{n} \sum_{k=j}^{n} m_{k} l_{j} \left(\dot{\omega_{j}} \times q_{j}\right) \\
& + \sum_{j=1}^{n} \sum_{k=j}^{n} m_{k} l_{j} \left(\omega_{j} \cdot \omega_{j}\right) q_{j}
\end{aligned}
$$

Coefficient of $\eta_{R}$, the equation of $R$:

$$
\begin{aligned}
0 = {}& -J \dot{\Omega_{R}} \\
& -\Omega_{R} \times \left(J \Omega_{R}\right) \\
& + M
\end{aligned}
$$

Coefficient of $\xi_{i}$, the equation of $q_i$:

$$
\begin{aligned}
0 = {}& -\sum_{j=i}^{n} m_{j} l_{i} \left(\left(\omega_{i} \times q_{i}\right) \times \dot{x_Q}\right) \\
& -\sum_{j=i}^{n} m_{j} l_{i} \left(q_{i} \times \ddot{x_Q}\right) \\
& + \sum_{j=i}^{n} 2 m_{j} l_{i} \left(\omega_{i} \times \left(q_{i} \times \dot{x_Q}\right)\right) \\
& -\sum_{j=i}^{n} m_{j} l_{i} \left(q_{i} \times \left(\omega_{i} \times \dot{x_Q}\right)\right) \\
& -\sum_{j=i}^{n} l_{i} m_{j} g \left(q_{i} \times e_3\right) \\
& -\sum_{j=1}^{n} \sum_{k=\max(i, j)}^{n} m_{k} l_{i} l_{j} \left(\left(\omega_{i} \times q_{i}\right) \times \left(\omega_{j} \times q_{j}\right)\right) \\
& -\sum_{j=1}^{n} \sum_{k=\max(i, j)}^{n} m_{k} l_{i} l_{j} \left(q_{i} \times \left(\dot{\omega_{j}} \times q_{j}\right)\right) \\
& + \sum_{j=1}^{n} \sum_{k=\max(i, j)}^{n} m_{k} l_{i} l_{j} \left(\omega_{j} \cdot \omega_{j}\right) \left(q_{i} \times q_{j}\right) \\
& + \sum_{j=1}^{n} \sum_{k=\max(i, j)}^{n} 2 m_{k} l_{i} l_{j} \left(\omega_{i} \times \left(q_{i} \times \left(\omega_{j} \times q_{j}\right)\right)\right) \\
& -\sum_{j=1}^{n} \sum_{k=\max(i, j)}^{n} m_{k} l_{i} l_{j} \left(q_{i} \times \left(\omega_{i} \times \left(\omega_{j} \times q_{j}\right)\right)\right)
\end{aligned}
$$
<!-- /model:indexed -->

## Energy, Lagrangian and virtual work at $n = 5$

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m_Q \left(\dot{x_Q} \cdot \dot{x_Q}\right) \\
& + \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right) \\
& + \frac{1}{2} m_1 \left(\left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right)\right) \cdot \left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right)\right)\right) \\
& + \frac{1}{2} m_2 \left(\left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right) \cdot \left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right)\right) \\
& + \frac{1}{2} m_3 \left(\left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right) \cdot \left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right)\right) \\
& + \frac{1}{2} m_4 \left(\left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right)\right) \cdot \left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right)\right)\right) \\
& + \frac{1}{2} m_5 \left(\left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right) + l_5 \left(\omega_5 \times q_5\right)\right) \cdot \left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right) + l_5 \left(\omega_5 \times q_5\right)\right)\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& m_Q g \left(e_3 \cdot x_Q\right) \\
& + m_1 g \left(e_3 \cdot \left(x_Q + l_1 q_1\right)\right) \\
& + m_2 g \left(e_3 \cdot \left(x_Q + l_1 q_1 + l_2 q_2\right)\right) \\
& + m_3 g \left(e_3 \cdot \left(x_Q + l_1 q_1 + l_2 q_2 + l_3 q_3\right)\right) \\
& + m_4 g \left(e_3 \cdot \left(x_Q + l_1 q_1 + l_2 q_2 + l_3 q_3 + l_4 q_4\right)\right) \\
& + m_5 g \left(e_3 \cdot \left(x_Q + l_1 q_1 + l_2 q_2 + l_3 q_3 + l_4 q_4 + l_5 q_5\right)\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m_Q \left(\dot{x_Q} \cdot \dot{x_Q}\right) \\
& + \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right) \\
& -m_Q g \left(e_3 \cdot x_Q\right) \\
& + \frac{1}{2} m_1 \left(\left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right)\right) \cdot \left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_1 g \left(e_3 \cdot \left(x_Q + l_1 q_1\right)\right) \\
& + \frac{1}{2} m_2 \left(\left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right) \cdot \left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_2 g \left(e_3 \cdot \left(x_Q + l_1 q_1 + l_2 q_2\right)\right) \\
& + \frac{1}{2} m_3 \left(\left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right) \cdot \left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_3 g \left(e_3 \cdot \left(x_Q + l_1 q_1 + l_2 q_2 + l_3 q_3\right)\right) \\
& + \frac{1}{2} m_4 \left(\left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right)\right) \cdot \left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_4 g \left(e_3 \cdot \left(x_Q + l_1 q_1 + l_2 q_2 + l_3 q_3 + l_4 q_4\right)\right) \\
& + \frac{1}{2} m_5 \left(\left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right) + l_5 \left(\omega_5 \times q_5\right)\right) \cdot \left(\dot{x_Q} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right) + l_5 \left(\omega_5 \times q_5\right)\right)\right) \\
& -m_5 g \left(e_3 \cdot \left(x_Q + l_1 q_1 + l_2 q_2 + l_3 q_3 + l_4 q_4 + l_5 q_5\right)\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \delta x_Q \cdot f R e_3 \\
& + \eta_{R} \cdot M
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion at $n = 5$

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); the $R$ equation is in the body frame and each $q_i$
equation holds on the tangent plane at $q_i$.

<!-- model:eom -->
<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\delta x_Q\)</span>, the equation of <span class="math notranslate nohighlight">\(x_Q\)</span> (43 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_Q \ddot{x_Q} \\
& -m_Q g e_3 \\
& -m_1 \ddot{x_Q} \\
& -m_1 l_1 \left(\dot{\omega_1} \times q_1\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_1 l_1 q_1 \\
& -m_1 g e_3 \\
& -m_2 \ddot{x_Q} \\
& -m_2 l_1 \left(\dot{\omega_1} \times q_1\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_2 l_1 q_1 \\
& -m_2 l_2 \left(\dot{\omega_2} \times q_2\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_2 l_2 q_2 \\
& -m_2 g e_3 \\
& -m_3 \ddot{x_Q} \\
& -m_3 l_1 \left(\dot{\omega_1} \times q_1\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_3 l_1 q_1 \\
& -m_3 l_2 \left(\dot{\omega_2} \times q_2\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_3 l_2 q_2 \\
& -m_3 l_3 \left(\dot{\omega_3} \times q_3\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_3 l_3 q_3 \\
& -m_3 g e_3 \\
& -m_4 \ddot{x_Q} \\
& -m_4 l_1 \left(\dot{\omega_1} \times q_1\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_4 l_1 q_1 \\
& -m_4 l_2 \left(\dot{\omega_2} \times q_2\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_4 l_2 q_2 \\
& -m_4 l_3 \left(\dot{\omega_3} \times q_3\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_4 l_3 q_3 \\
& -m_4 l_4 \left(\dot{\omega_4} \times q_4\right) \\
& + \left(\omega_4 \cdot \omega_4\right) m_4 l_4 q_4 \\
& -m_4 g e_3 \\
& -m_5 \ddot{x_Q} \\
& -m_5 l_1 \left(\dot{\omega_1} \times q_1\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_5 l_1 q_1 \\
& -m_5 l_2 \left(\dot{\omega_2} \times q_2\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_5 l_2 q_2 \\
& -m_5 l_3 \left(\dot{\omega_3} \times q_3\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_5 l_3 q_3 \\
& -m_5 l_4 \left(\dot{\omega_4} \times q_4\right) \\
& + \left(\omega_4 \cdot \omega_4\right) m_5 l_4 q_4 \\
& -m_5 l_5 \left(\dot{\omega_5} \times q_5\right) \\
& + \left(\omega_5 \cdot \omega_5\right) m_5 l_5 q_5 \\
& -m_5 g e_3 \\
& + f R e_3
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_1\)</span>, the equation of <span class="math notranslate nohighlight">\(q_1\)</span> (80 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_1 l_1 \left(\left(\omega_1 \times q_1\right) \times \dot{x_Q}\right) \\
& -m_1 l_1 \left(q_1 \times \ddot{x_Q}\right) \\
& -m_1 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + 2 m_1 l_1 \left(\omega_1 \times \left(q_1 \times \dot{x_Q}\right)\right) \\
& -m_1 l_1 \left(q_1 \times \left(\omega_1 \times \dot{x_Q}\right)\right) \\
& -l_1 m_1 g \left(q_1 \times e_3\right) \\
& -m_2 l_1 \left(\left(\omega_1 \times q_1\right) \times \dot{x_Q}\right) \\
& -m_2 l_1 \left(q_1 \times \ddot{x_Q}\right) \\
& -m_2 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -m_2 l_1 l_2 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_2 l_1 l_2 \left(q_1 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_2 l_1 l_2 \left(q_1 \times q_2\right) \\
& + 2 m_2 l_1 \left(\omega_1 \times \left(q_1 \times \dot{x_Q}\right)\right) \\
& + 2 m_2 l_1 l_2 \left(\omega_1 \times \left(q_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_2 l_1 \left(q_1 \times \left(\omega_1 \times \dot{x_Q}\right)\right) \\
& -m_2 l_1 l_2 \left(q_1 \times \left(\omega_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -l_1 m_2 g \left(q_1 \times e_3\right) \\
& -m_3 l_1 \left(\left(\omega_1 \times q_1\right) \times \dot{x_Q}\right) \\
& -m_3 l_1 \left(q_1 \times \ddot{x_Q}\right) \\
& -m_3 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -m_3 l_1 l_2 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_3 l_1 l_2 \left(q_1 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_3 l_1 l_2 \left(q_1 \times q_2\right) \\
& -m_3 l_1 l_3 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_3 l_1 l_3 \left(q_1 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_3 l_1 l_3 \left(q_1 \times q_3\right) \\
& + 2 m_3 l_1 \left(\omega_1 \times \left(q_1 \times \dot{x_Q}\right)\right) \\
& + 2 m_3 l_1 l_2 \left(\omega_1 \times \left(q_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_3 l_1 l_3 \left(\omega_1 \times \left(q_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_3 l_1 \left(q_1 \times \left(\omega_1 \times \dot{x_Q}\right)\right) \\
& -m_3 l_1 l_2 \left(q_1 \times \left(\omega_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_3 l_1 l_3 \left(q_1 \times \left(\omega_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -l_1 m_3 g \left(q_1 \times e_3\right) \\
& -m_4 l_1 \left(\left(\omega_1 \times q_1\right) \times \dot{x_Q}\right) \\
& -m_4 l_1 \left(q_1 \times \ddot{x_Q}\right) \\
& -m_4 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -m_4 l_1 l_2 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_4 l_1 l_2 \left(q_1 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_4 l_1 l_2 \left(q_1 \times q_2\right) \\
& -m_4 l_1 l_3 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_4 l_1 l_3 \left(q_1 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_4 l_1 l_3 \left(q_1 \times q_3\right) \\
& -m_4 l_1 l_4 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_4 \times q_4\right)\right) \\
& -m_4 l_1 l_4 \left(q_1 \times \left(\dot{\omega_4} \times q_4\right)\right) \\
& + \left(\omega_4 \cdot \omega_4\right) m_4 l_1 l_4 \left(q_1 \times q_4\right) \\
& + 2 m_4 l_1 \left(\omega_1 \times \left(q_1 \times \dot{x_Q}\right)\right) \\
& + 2 m_4 l_1 l_2 \left(\omega_1 \times \left(q_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_4 l_1 l_3 \left(\omega_1 \times \left(q_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& + 2 m_4 l_1 l_4 \left(\omega_1 \times \left(q_1 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_4 l_1 \left(q_1 \times \left(\omega_1 \times \dot{x_Q}\right)\right) \\
& -m_4 l_1 l_2 \left(q_1 \times \left(\omega_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_4 l_1 l_3 \left(q_1 \times \left(\omega_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_4 l_1 l_4 \left(q_1 \times \left(\omega_1 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -l_1 m_4 g \left(q_1 \times e_3\right) \\
& -m_5 l_1 \left(\left(\omega_1 \times q_1\right) \times \dot{x_Q}\right) \\
& -m_5 l_1 \left(q_1 \times \ddot{x_Q}\right) \\
& -m_5 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -m_5 l_1 l_2 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_5 l_1 l_2 \left(q_1 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_5 l_1 l_2 \left(q_1 \times q_2\right) \\
& -m_5 l_1 l_3 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_5 l_1 l_3 \left(q_1 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_5 l_1 l_3 \left(q_1 \times q_3\right) \\
& -m_5 l_1 l_4 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_4 \times q_4\right)\right) \\
& -m_5 l_1 l_4 \left(q_1 \times \left(\dot{\omega_4} \times q_4\right)\right) \\
& + \left(\omega_4 \cdot \omega_4\right) m_5 l_1 l_4 \left(q_1 \times q_4\right) \\
& -m_5 l_1 l_5 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_5 \times q_5\right)\right) \\
& -m_5 l_1 l_5 \left(q_1 \times \left(\dot{\omega_5} \times q_5\right)\right) \\
& + \left(\omega_5 \cdot \omega_5\right) m_5 l_1 l_5 \left(q_1 \times q_5\right) \\
& + 2 m_5 l_1 \left(\omega_1 \times \left(q_1 \times \dot{x_Q}\right)\right) \\
& + 2 m_5 l_1 l_2 \left(\omega_1 \times \left(q_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_5 l_1 l_3 \left(\omega_1 \times \left(q_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& + 2 m_5 l_1 l_4 \left(\omega_1 \times \left(q_1 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& + 2 m_5 l_1 l_5 \left(\omega_1 \times \left(q_1 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -m_5 l_1 \left(q_1 \times \left(\omega_1 \times \dot{x_Q}\right)\right) \\
& -m_5 l_1 l_2 \left(q_1 \times \left(\omega_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_5 l_1 l_3 \left(q_1 \times \left(\omega_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_5 l_1 l_4 \left(q_1 \times \left(\omega_1 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_5 l_1 l_5 \left(q_1 \times \left(\omega_1 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -l_1 m_5 g \left(q_1 \times e_3\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_2\)</span>, the equation of <span class="math notranslate nohighlight">\(q_2\)</span> (74 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_2 l_2 \left(\left(\omega_2 \times q_2\right) \times \dot{x_Q}\right) \\
& -m_2 l_2 \left(q_2 \times \ddot{x_Q}\right) \\
& -m_2 l_2 l_1 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_2 l_2 l_1 \left(q_2 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_2 l_2 l_1 \left(q_2 \times q_1\right) \\
& -m_2 l_2 l_2 \left(q_2 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + 2 m_2 l_2 \left(\omega_2 \times \left(q_2 \times \dot{x_Q}\right)\right) \\
& + 2 m_2 l_2 l_1 \left(\omega_2 \times \left(q_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_2 l_2 \left(q_2 \times \left(\omega_2 \times \dot{x_Q}\right)\right) \\
& -m_2 l_2 l_1 \left(q_2 \times \left(\omega_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -l_2 m_2 g \left(q_2 \times e_3\right) \\
& -m_3 l_2 \left(\left(\omega_2 \times q_2\right) \times \dot{x_Q}\right) \\
& -m_3 l_2 \left(q_2 \times \ddot{x_Q}\right) \\
& -m_3 l_2 l_1 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_3 l_2 l_1 \left(q_2 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_3 l_2 l_1 \left(q_2 \times q_1\right) \\
& -m_3 l_2 l_2 \left(q_2 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& -m_3 l_2 l_3 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_3 l_2 l_3 \left(q_2 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_3 l_2 l_3 \left(q_2 \times q_3\right) \\
& + 2 m_3 l_2 \left(\omega_2 \times \left(q_2 \times \dot{x_Q}\right)\right) \\
& + 2 m_3 l_2 l_1 \left(\omega_2 \times \left(q_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_3 l_2 l_3 \left(\omega_2 \times \left(q_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_3 l_2 \left(q_2 \times \left(\omega_2 \times \dot{x_Q}\right)\right) \\
& -m_3 l_2 l_1 \left(q_2 \times \left(\omega_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_3 l_2 l_3 \left(q_2 \times \left(\omega_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -l_2 m_3 g \left(q_2 \times e_3\right) \\
& -m_4 l_2 \left(\left(\omega_2 \times q_2\right) \times \dot{x_Q}\right) \\
& -m_4 l_2 \left(q_2 \times \ddot{x_Q}\right) \\
& -m_4 l_2 l_1 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_4 l_2 l_1 \left(q_2 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_4 l_2 l_1 \left(q_2 \times q_1\right) \\
& -m_4 l_2 l_2 \left(q_2 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& -m_4 l_2 l_3 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_4 l_2 l_3 \left(q_2 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_4 l_2 l_3 \left(q_2 \times q_3\right) \\
& -m_4 l_2 l_4 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_4 \times q_4\right)\right) \\
& -m_4 l_2 l_4 \left(q_2 \times \left(\dot{\omega_4} \times q_4\right)\right) \\
& + \left(\omega_4 \cdot \omega_4\right) m_4 l_2 l_4 \left(q_2 \times q_4\right) \\
& + 2 m_4 l_2 \left(\omega_2 \times \left(q_2 \times \dot{x_Q}\right)\right) \\
& + 2 m_4 l_2 l_1 \left(\omega_2 \times \left(q_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_4 l_2 l_3 \left(\omega_2 \times \left(q_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& + 2 m_4 l_2 l_4 \left(\omega_2 \times \left(q_2 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_4 l_2 \left(q_2 \times \left(\omega_2 \times \dot{x_Q}\right)\right) \\
& -m_4 l_2 l_1 \left(q_2 \times \left(\omega_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_4 l_2 l_3 \left(q_2 \times \left(\omega_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_4 l_2 l_4 \left(q_2 \times \left(\omega_2 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -l_2 m_4 g \left(q_2 \times e_3\right) \\
& -m_5 l_2 \left(\left(\omega_2 \times q_2\right) \times \dot{x_Q}\right) \\
& -m_5 l_2 \left(q_2 \times \ddot{x_Q}\right) \\
& -m_5 l_2 l_1 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_5 l_2 l_1 \left(q_2 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_5 l_2 l_1 \left(q_2 \times q_1\right) \\
& -m_5 l_2 l_2 \left(q_2 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& -m_5 l_2 l_3 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_5 l_2 l_3 \left(q_2 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_5 l_2 l_3 \left(q_2 \times q_3\right) \\
& -m_5 l_2 l_4 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_4 \times q_4\right)\right) \\
& -m_5 l_2 l_4 \left(q_2 \times \left(\dot{\omega_4} \times q_4\right)\right) \\
& + \left(\omega_4 \cdot \omega_4\right) m_5 l_2 l_4 \left(q_2 \times q_4\right) \\
& -m_5 l_2 l_5 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_5 \times q_5\right)\right) \\
& -m_5 l_2 l_5 \left(q_2 \times \left(\dot{\omega_5} \times q_5\right)\right) \\
& + \left(\omega_5 \cdot \omega_5\right) m_5 l_2 l_5 \left(q_2 \times q_5\right) \\
& + 2 m_5 l_2 \left(\omega_2 \times \left(q_2 \times \dot{x_Q}\right)\right) \\
& + 2 m_5 l_2 l_1 \left(\omega_2 \times \left(q_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_5 l_2 l_3 \left(\omega_2 \times \left(q_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& + 2 m_5 l_2 l_4 \left(\omega_2 \times \left(q_2 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& + 2 m_5 l_2 l_5 \left(\omega_2 \times \left(q_2 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -m_5 l_2 \left(q_2 \times \left(\omega_2 \times \dot{x_Q}\right)\right) \\
& -m_5 l_2 l_1 \left(q_2 \times \left(\omega_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_5 l_2 l_3 \left(q_2 \times \left(\omega_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_5 l_2 l_4 \left(q_2 \times \left(\omega_2 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_5 l_2 l_5 \left(q_2 \times \left(\omega_2 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -l_2 m_5 g \left(q_2 \times e_3\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_3\)</span>, the equation of <span class="math notranslate nohighlight">\(q_3\)</span> (63 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_3 l_3 \left(\left(\omega_3 \times q_3\right) \times \dot{x_Q}\right) \\
& -m_3 l_3 \left(q_3 \times \ddot{x_Q}\right) \\
& -m_3 l_3 l_1 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_3 l_3 l_1 \left(q_3 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_3 l_3 l_1 \left(q_3 \times q_1\right) \\
& -m_3 l_3 l_2 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_3 l_3 l_2 \left(q_3 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_3 l_3 l_2 \left(q_3 \times q_2\right) \\
& -m_3 l_3 l_3 \left(q_3 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + 2 m_3 l_3 \left(\omega_3 \times \left(q_3 \times \dot{x_Q}\right)\right) \\
& + 2 m_3 l_3 l_1 \left(\omega_3 \times \left(q_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_3 l_3 l_2 \left(\omega_3 \times \left(q_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_3 l_3 \left(q_3 \times \left(\omega_3 \times \dot{x_Q}\right)\right) \\
& -m_3 l_3 l_1 \left(q_3 \times \left(\omega_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_3 l_3 l_2 \left(q_3 \times \left(\omega_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -l_3 m_3 g \left(q_3 \times e_3\right) \\
& -m_4 l_3 \left(\left(\omega_3 \times q_3\right) \times \dot{x_Q}\right) \\
& -m_4 l_3 \left(q_3 \times \ddot{x_Q}\right) \\
& -m_4 l_3 l_1 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_4 l_3 l_1 \left(q_3 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_4 l_3 l_1 \left(q_3 \times q_1\right) \\
& -m_4 l_3 l_2 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_4 l_3 l_2 \left(q_3 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_4 l_3 l_2 \left(q_3 \times q_2\right) \\
& -m_4 l_3 l_3 \left(q_3 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& -m_4 l_3 l_4 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_4 \times q_4\right)\right) \\
& -m_4 l_3 l_4 \left(q_3 \times \left(\dot{\omega_4} \times q_4\right)\right) \\
& + \left(\omega_4 \cdot \omega_4\right) m_4 l_3 l_4 \left(q_3 \times q_4\right) \\
& + 2 m_4 l_3 \left(\omega_3 \times \left(q_3 \times \dot{x_Q}\right)\right) \\
& + 2 m_4 l_3 l_1 \left(\omega_3 \times \left(q_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_4 l_3 l_2 \left(\omega_3 \times \left(q_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_4 l_3 l_4 \left(\omega_3 \times \left(q_3 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_4 l_3 \left(q_3 \times \left(\omega_3 \times \dot{x_Q}\right)\right) \\
& -m_4 l_3 l_1 \left(q_3 \times \left(\omega_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_4 l_3 l_2 \left(q_3 \times \left(\omega_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_4 l_3 l_4 \left(q_3 \times \left(\omega_3 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -l_3 m_4 g \left(q_3 \times e_3\right) \\
& -m_5 l_3 \left(\left(\omega_3 \times q_3\right) \times \dot{x_Q}\right) \\
& -m_5 l_3 \left(q_3 \times \ddot{x_Q}\right) \\
& -m_5 l_3 l_1 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_5 l_3 l_1 \left(q_3 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_5 l_3 l_1 \left(q_3 \times q_1\right) \\
& -m_5 l_3 l_2 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_5 l_3 l_2 \left(q_3 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_5 l_3 l_2 \left(q_3 \times q_2\right) \\
& -m_5 l_3 l_3 \left(q_3 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& -m_5 l_3 l_4 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_4 \times q_4\right)\right) \\
& -m_5 l_3 l_4 \left(q_3 \times \left(\dot{\omega_4} \times q_4\right)\right) \\
& + \left(\omega_4 \cdot \omega_4\right) m_5 l_3 l_4 \left(q_3 \times q_4\right) \\
& -m_5 l_3 l_5 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_5 \times q_5\right)\right) \\
& -m_5 l_3 l_5 \left(q_3 \times \left(\dot{\omega_5} \times q_5\right)\right) \\
& + \left(\omega_5 \cdot \omega_5\right) m_5 l_3 l_5 \left(q_3 \times q_5\right) \\
& + 2 m_5 l_3 \left(\omega_3 \times \left(q_3 \times \dot{x_Q}\right)\right) \\
& + 2 m_5 l_3 l_1 \left(\omega_3 \times \left(q_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_5 l_3 l_2 \left(\omega_3 \times \left(q_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_5 l_3 l_4 \left(\omega_3 \times \left(q_3 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& + 2 m_5 l_3 l_5 \left(\omega_3 \times \left(q_3 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -m_5 l_3 \left(q_3 \times \left(\omega_3 \times \dot{x_Q}\right)\right) \\
& -m_5 l_3 l_1 \left(q_3 \times \left(\omega_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_5 l_3 l_2 \left(q_3 \times \left(\omega_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_5 l_3 l_4 \left(q_3 \times \left(\omega_3 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_5 l_3 l_5 \left(q_3 \times \left(\omega_3 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -l_3 m_5 g \left(q_3 \times e_3\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_4\)</span>, the equation of <span class="math notranslate nohighlight">\(q_4\)</span> (47 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_4 l_4 \left(\left(\omega_4 \times q_4\right) \times \dot{x_Q}\right) \\
& -m_4 l_4 \left(q_4 \times \ddot{x_Q}\right) \\
& -m_4 l_4 l_1 \left(\left(\omega_4 \times q_4\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_4 l_4 l_1 \left(q_4 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_4 l_4 l_1 \left(q_4 \times q_1\right) \\
& -m_4 l_4 l_2 \left(\left(\omega_4 \times q_4\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_4 l_4 l_2 \left(q_4 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_4 l_4 l_2 \left(q_4 \times q_2\right) \\
& -m_4 l_4 l_3 \left(\left(\omega_4 \times q_4\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_4 l_4 l_3 \left(q_4 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_4 l_4 l_3 \left(q_4 \times q_3\right) \\
& -m_4 l_4 l_4 \left(q_4 \times \left(\dot{\omega_4} \times q_4\right)\right) \\
& + 2 m_4 l_4 \left(\omega_4 \times \left(q_4 \times \dot{x_Q}\right)\right) \\
& + 2 m_4 l_4 l_1 \left(\omega_4 \times \left(q_4 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_4 l_4 l_2 \left(\omega_4 \times \left(q_4 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_4 l_4 l_3 \left(\omega_4 \times \left(q_4 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_4 l_4 \left(q_4 \times \left(\omega_4 \times \dot{x_Q}\right)\right) \\
& -m_4 l_4 l_1 \left(q_4 \times \left(\omega_4 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_4 l_4 l_2 \left(q_4 \times \left(\omega_4 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_4 l_4 l_3 \left(q_4 \times \left(\omega_4 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -l_4 m_4 g \left(q_4 \times e_3\right) \\
& -m_5 l_4 \left(\left(\omega_4 \times q_4\right) \times \dot{x_Q}\right) \\
& -m_5 l_4 \left(q_4 \times \ddot{x_Q}\right) \\
& -m_5 l_4 l_1 \left(\left(\omega_4 \times q_4\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_5 l_4 l_1 \left(q_4 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_5 l_4 l_1 \left(q_4 \times q_1\right) \\
& -m_5 l_4 l_2 \left(\left(\omega_4 \times q_4\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_5 l_4 l_2 \left(q_4 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_5 l_4 l_2 \left(q_4 \times q_2\right) \\
& -m_5 l_4 l_3 \left(\left(\omega_4 \times q_4\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_5 l_4 l_3 \left(q_4 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_5 l_4 l_3 \left(q_4 \times q_3\right) \\
& -m_5 l_4 l_4 \left(q_4 \times \left(\dot{\omega_4} \times q_4\right)\right) \\
& -m_5 l_4 l_5 \left(\left(\omega_4 \times q_4\right) \times \left(\omega_5 \times q_5\right)\right) \\
& -m_5 l_4 l_5 \left(q_4 \times \left(\dot{\omega_5} \times q_5\right)\right) \\
& + \left(\omega_5 \cdot \omega_5\right) m_5 l_4 l_5 \left(q_4 \times q_5\right) \\
& + 2 m_5 l_4 \left(\omega_4 \times \left(q_4 \times \dot{x_Q}\right)\right) \\
& + 2 m_5 l_4 l_1 \left(\omega_4 \times \left(q_4 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_5 l_4 l_2 \left(\omega_4 \times \left(q_4 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_5 l_4 l_3 \left(\omega_4 \times \left(q_4 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& + 2 m_5 l_4 l_5 \left(\omega_4 \times \left(q_4 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -m_5 l_4 \left(q_4 \times \left(\omega_4 \times \dot{x_Q}\right)\right) \\
& -m_5 l_4 l_1 \left(q_4 \times \left(\omega_4 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_5 l_4 l_2 \left(q_4 \times \left(\omega_4 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_5 l_4 l_3 \left(q_4 \times \left(\omega_4 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_5 l_4 l_5 \left(q_4 \times \left(\omega_4 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -l_4 m_5 g \left(q_4 \times e_3\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_5\)</span>, the equation of <span class="math notranslate nohighlight">\(q_5\)</span> (26 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_5 l_5 \left(\left(\omega_5 \times q_5\right) \times \dot{x_Q}\right) \\
& -m_5 l_5 \left(q_5 \times \ddot{x_Q}\right) \\
& -m_5 l_5 l_1 \left(\left(\omega_5 \times q_5\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_5 l_5 l_1 \left(q_5 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_5 l_5 l_1 \left(q_5 \times q_1\right) \\
& -m_5 l_5 l_2 \left(\left(\omega_5 \times q_5\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_5 l_5 l_2 \left(q_5 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_5 l_5 l_2 \left(q_5 \times q_2\right) \\
& -m_5 l_5 l_3 \left(\left(\omega_5 \times q_5\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_5 l_5 l_3 \left(q_5 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_5 l_5 l_3 \left(q_5 \times q_3\right) \\
& -m_5 l_5 l_4 \left(\left(\omega_5 \times q_5\right) \times \left(\omega_4 \times q_4\right)\right) \\
& -m_5 l_5 l_4 \left(q_5 \times \left(\dot{\omega_4} \times q_4\right)\right) \\
& + \left(\omega_4 \cdot \omega_4\right) m_5 l_5 l_4 \left(q_5 \times q_4\right) \\
& -m_5 l_5 l_5 \left(q_5 \times \left(\dot{\omega_5} \times q_5\right)\right) \\
& + 2 m_5 l_5 \left(\omega_5 \times \left(q_5 \times \dot{x_Q}\right)\right) \\
& + 2 m_5 l_5 l_1 \left(\omega_5 \times \left(q_5 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_5 l_5 l_2 \left(\omega_5 \times \left(q_5 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_5 l_5 l_3 \left(\omega_5 \times \left(q_5 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& + 2 m_5 l_5 l_4 \left(\omega_5 \times \left(q_5 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_5 l_5 \left(q_5 \times \left(\omega_5 \times \dot{x_Q}\right)\right) \\
& -m_5 l_5 l_1 \left(q_5 \times \left(\omega_5 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_5 l_5 l_2 \left(q_5 \times \left(\omega_5 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_5 l_5 l_3 \left(q_5 \times \left(\omega_5 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_5 l_5 l_4 \left(q_5 \times \left(\omega_5 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -l_5 m_5 g \left(q_5 \times e_3\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_{R}\)</span>, the equation of <span class="math notranslate nohighlight">\(R\)</span> (3 terms)</summary>

$$
\begin{aligned}
0 = {}& -J \dot{\Omega_{R}} \\
& -\Omega_{R} \times \left(J \Omega_{R}\right) \\
& + M
\end{aligned}
$$

</details>
<!-- /model:eom -->

## Code

The model for every $n$, as `geomech.indexed.systems.quadrotor_flexible_cable()` writes
it, and the library system for $n = 5$:

```python
from geomech import S2, SO3, Dot, Matrix, Scalar, Vector
from geomech.core.base.constants import HALF
from geomech.indexed import Family, IndexedVariables, Sum, compute_indexed_eom, indices, n
from geomech.systems import get

i, j = indices("i j")
mQ, g, f = Scalar("m_Q", attr=["Constant"]), Scalar("g", attr=["Constant"]), Scalar("f")
J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
e3 = Vector("e3", attr=["Constant"])
xQ, R, M = Vector("x_Q"), SO3("R"), Vector("M")
q = Family(S2, "q")  # q_1, ..., q_n
m, l = Family(Scalar, "m", ["Constant"]), Family(Scalar, "l", ["Constant"])
Om, eta = R.get_tangent_vector(), R.get_variation_vector()
vQ = xQ.t_diff()


def x(k):  # mass k at x_Q + l_1 q_1 + ... + l_k q_k
    return xQ + Sum(j, 1, k, q[j] * l[j])


L = HALF * mQ * Dot(vQ, vQ) + HALF * Dot(Om, J * Om) - mQ * g * Dot(e3, xQ)
L = L + Sum(i, 1, n, HALF * m[i] * x(i).t_diff().dot(x(i).t_diff()))
L = L - Sum(i, 1, n, m[i] * g * x(i).dot(e3))
dW = Dot(xQ.delta(), (R * e3) * f) + Dot(eta, M)
variables = IndexedVariables(families=[q], vectors=[xQ], matrices=[R])
eom = compute_indexed_eom(L, dW, variables)  # every n

system = get("quadrotor_flexible_cable", n_links=5)  # n from 1 to 20
eom_5 = system.derive()
```

The same model at $n = 5$ is `geomech.systems.get("quadrotor_flexible_cable", n_links=5)`,
with default parameters, a default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-5}$ |
| Energy $H = T + V$, 6 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-9}$, with $H(0) = 36.1$ |
| Equations for every $n$ above, at $n = 5$: `check_eom` against the Euler–Lagrange evaluation of the generic $L$, $\delta W$ at $n = 5$ | below $10^{-5}$ |
| Generic $L$ at $n = 5$ against the library system's $L$, at random states | relative difference below $10^{-14}$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
