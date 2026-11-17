# Two quadrotors carrying a hose

Two quadrotors carry a flexible hose, one at each end. The hose is $n$ massless links of
lengths $l_i$ and directions $q_i \in S^2$, joined by spherical joints, with point masses
$m_k$ at the joints $x_k = x_0 + \sum_{i \le k} l_i q_i$, $k = 0, \dots, n$. The
quadrotors (mass $m_{Qk}$, inertia $J_k$, attitude $R_k$) sit at the end joints $x_0$ and
$x_n$, with their centres of mass on the joint. The configuration manifold is
$\mathbb{R}^3 \times (S^2)^n \times SO(3)^2$. `geomech.indexed` has no hose model yet, so
this page derives the default $n = 3$ only. Model of Kotaru and Sreenath (2020); the
equations are verified numerically (`check_eom`). Library system `multi_quadrotor_hose`.

## Animation

Six seconds with constant thrust equal to the total weight (each vehicle carries an equal
share) and zero moments: open loop, with no controller. The start is not the library's
initial state: the quadrotors start level and not rotating, while the hose keeps the
library's initial state, sagging between ends 1.6 m apart. With $M_0 = M_3 = 0$ and the thrust at the centre of mass, the attitude equations do not involve
the hose, so both vehicles stay exactly level. With the thrusts vertical, nothing balances
the inward pull of the hose, and it still folds: the ends swing towards each other, down
to 0.27 m apart at 2.1 s (0.52 m at 6 s), and the middle link turns over. The system
drifts 1.7 m sideways. Drag to orbit; the full report (equations, standard form,
parameters) opens in a
<a href="../model-reports/multi_quadrotor_hose.html">separate page</a>.

<iframe class="model-report" src="../model-reports/multi_quadrotor_hose.html?view" title="Two quadrotors carrying a hose, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 280" role="img" aria-labelledby="mqh-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="mqh-title">Two quadrotors carrying a hose: joints x_0 to x_3 joined by links of directions q_1, q_2, q_3 and lengths l_1, l_2, l_3; point masses at the joints; quadrotors 0 and 3 at the ends with attitudes R_0, R_3, thrusts f_0, f_3 and moments M_0, M_3</title>
  <defs>
    <marker id="mqh-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <line x1="30" y1="266" x2="30" y2="226" marker-end="url(#mqh-arr)"/>
    <line x1="30" y1="266" x2="70" y2="266" marker-end="url(#mqh-arr)"/>
  </g>
  <g fill="none" stroke="currentColor" stroke-width="2">
    <line x1="70" y1="100" x2="130" y2="180"/>
    <line x1="130" y1="180" x2="230" y2="190"/>
    <line x1="230" y1="190" x2="300" y2="110"/>
  </g>
  <!-- quadrotor at (70,100) -->
  <g fill="none" stroke="currentColor" stroke-width="1.6">
    <line x1="32" y1="100" x2="108" y2="100"/>
    <line x1="50" y1="110" x2="90" y2="90"/>
    <ellipse cx="32" cy="100" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="108" cy="100" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="90" cy="90" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="50" cy="110" rx="13" ry="4" stroke-width="1.2"/>
  </g>
  <rect x="61" y="95" width="18" height="10" rx="2" fill="currentColor"/>
  <!-- quadrotor at (300,110) -->
  <g fill="none" stroke="currentColor" stroke-width="1.6">
    <line x1="262" y1="110" x2="338" y2="110"/>
    <line x1="320" y1="120" x2="280" y2="100"/>
    <ellipse cx="262" cy="110" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="338" cy="110" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="280" cy="100" rx="13" ry="4" stroke-width="1.2"/>
    <ellipse cx="320" cy="120" rx="13" ry="4" stroke-width="1.2"/>
  </g>
  <rect x="291" y="105" width="18" height="10" rx="2" fill="currentColor"/>
  <line x1="70.0" y1="92.0" x2="70.0" y2="48.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#mqh-arr)"/>
  <path d="M 57,67.8 A 13 4 0 1 0 83,67.8" fill="none" stroke="#e08a2c" stroke-width="1.6" marker-end="url(#mqh-arr)"/>
  <line x1="300.0" y1="102.0" x2="300.0" y2="58.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#mqh-arr)"/>
  <path d="M 287,77.8 A 13 4 0 1 0 313,77.8" fill="none" stroke="#e08a2c" stroke-width="1.6" marker-end="url(#mqh-arr)"/>
  <line x1="70.0" y1="100.0" x2="94.0" y2="132.0" stroke="#3f8fa0" stroke-width="3" marker-end="url(#mqh-arr)"/>
  <line x1="130.0" y1="180.0" x2="169.8" y2="184.0" stroke="#3f8fa0" stroke-width="3" marker-end="url(#mqh-arr)"/>
  <line x1="230.0" y1="190.0" x2="256.3" y2="159.9" stroke="#3f8fa0" stroke-width="3" marker-end="url(#mqh-arr)"/>
  <circle cx="130" cy="180" r="5" fill="currentColor"/>
  <line x1="130.0" y1="186.0" x2="130.0" y2="224.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#mqh-arr)"/>
  <circle cx="230" cy="190" r="5" fill="currentColor"/>
  <line x1="230.0" y1="196.0" x2="230.0" y2="234.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#mqh-arr)"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="76" y="52" fill="#e08a2c">f<tspan font-size="12" dy="4">0</tspan></text>
    <text x="88" y="74" fill="#e08a2c">M<tspan font-size="12" dy="4">0</tspan></text>
    <text x="306" y="62" fill="#e08a2c">f<tspan font-size="12" dy="4">3</tspan></text>
    <text x="318" y="84" fill="#e08a2c">M<tspan font-size="12" dy="4">3</tspan></text>
    <text x="14" y="84" fill="#3f8fa0">x<tspan font-size="12" dy="4">0</tspan></text>
    <text x="8" y="140" fill="#3f8fa0">R<tspan font-size="12" dy="4">0</tspan><tspan dy="-4">, Ω</tspan><tspan font-size="12" dy="4">0</tspan></text>
    <text x="8" y="186">m<tspan font-size="12" dy="4">0</tspan><tspan dy="-4">, m</tspan><tspan font-size="12" dy="4">Q0</tspan><tspan dy="-4">, J</tspan><tspan font-size="12" dy="4">0</tspan></text>
    <text x="292" y="145" fill="#3f8fa0">R<tspan font-size="12" dy="4">3</tspan><tspan dy="-4">, Ω</tspan><tspan font-size="12" dy="4">3</tspan></text>
    <text x="262" y="204">m<tspan font-size="12" dy="4">3</tspan><tspan dy="-4">, m</tspan><tspan font-size="12" dy="4">Q3</tspan><tspan dy="-4">, J</tspan><tspan font-size="12" dy="4">3</tspan></text>
    <text x="104" y="128" fill="#3f8fa0">q<tspan font-size="12" dy="4">1</tspan></text>
    <text x="80" y="160">l<tspan font-size="12" dy="4">1</tspan></text>
    <text x="160" y="208" fill="#3f8fa0">q<tspan font-size="12" dy="4">2</tspan></text>
    <text x="196" y="174">l<tspan font-size="12" dy="4">2</tspan></text>
    <text x="268" y="176" fill="#3f8fa0">q<tspan font-size="12" dy="4">3</tspan></text>
    <text x="250" y="140">l<tspan font-size="12" dy="4">3</tspan></text>
    <text x="88" y="220" fill="#e08a2c">m<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"> g</tspan></text>
    <text x="236" y="230" fill="#e08a2c">m<tspan font-size="12" dy="4">2</tspan><tspan dy="-4"> g</tspan></text>
    <text x="24" y="218">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="76" y="271">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>Each <span class="math notranslate nohighlight">\(q_i\)</span> points from joint <span class="math notranslate nohighlight">\(i - 1\)</span> to joint <span class="math notranslate nohighlight">\(i\)</span>, from quadrotor
0 towards quadrotor 3. The end joints carry the hose masses <span class="math notranslate nohighlight">\(m_0\)</span>, <span class="math notranslate nohighlight">\(m_3\)</span> and the
quadrotors <span class="math notranslate nohighlight">\(m_{Q0}\)</span>, <span class="math notranslate nohighlight">\(m_{Q3}\)</span>. Each thrust <span class="math notranslate nohighlight">\(f_k\)</span>
acts along <span class="math notranslate nohighlight">\(R_k e_3\)</span> at its joint; <span class="math notranslate nohighlight">\(M_k\)</span> are the body-frame
moments. Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span> on every mass (shown on the inner joints only).</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $x_0$ | first joint, where quadrotor 0 sits | $\mathbb{R}^3$ |
| $\dot x_0$ | velocity of the first joint | $\mathbb{R}^3$ |
| $\delta x_0$ | variation of $x_0$ | $\mathbb{R}^3$ |
| $q_i$ | direction of link $i$, joint $i-1$ to joint $i$ | $S^2$ |
| $\omega_i$ | angular velocity of link $i$, $\dot q_i = \omega_i \times q_i$, $q_i \cdot \omega_i = 0$ | $T_{q_i} S^2$ |
| $\xi_i$ | variation of $q_i$: $\delta q_i = \xi_i \times q_i$ | $T_{q_i} S^2$ |
| $R_0$ | attitude of quadrotor 0, body frame to inertial frame | $SO(3)$ |
| $R_3$ | attitude of quadrotor 3, body frame to inertial frame | $SO(3)$ |
| $\Omega_0$ | body angular velocity of quadrotor 0, $\dot R_0 = R_0 \widehat{\Omega_0}$ | $\mathbb{R}^3$ |
| $\Omega_3$ | body angular velocity of quadrotor 3, $\dot R_3 = R_3 \widehat{\Omega_3}$ | $\mathbb{R}^3$ |
| $\eta_0$ | variation of $R_0$: $\delta R_0 = R_0 \widehat{\eta_0}$ | $\mathbb{R}^3$ |
| $\eta_3$ | variation of $R_3$: $\delta R_3 = R_3 \widehat{\eta_3}$ | $\mathbb{R}^3$ |
| $f_0$ | thrust of quadrotor 0 along $R_0 e_3$ | $\mathbb{R}$ |
| $f_3$ | thrust of quadrotor 3 along $R_3 e_3$ | $\mathbb{R}$ |
| $M_0$ | moment of quadrotor 0, body frame | $\mathbb{R}^3$ |
| $M_3$ | moment of quadrotor 3, body frame | $\mathbb{R}^3$ |
| $m_0, \dots, m_3, l_i$ | hose masses at the joints, link lengths | constants |
| $m_{Q0}, m_{Q3}, J_0, J_3, g, e_3$ | quadrotor masses and inertias, gravity, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m_0 \left(\dot{x_0} \cdot \dot{x_0}\right) \\
& + \frac{1}{2} m_1 \left(\left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right)\right) \cdot \left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right)\right)\right) \\
& + \frac{1}{2} m_2 \left(\left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right) \cdot \left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right)\right) \\
& + \frac{1}{2} m_3 \left(\left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right) \cdot \left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right)\right) \\
& + \frac{1}{2} m_{Q0} \left(\dot{x_0} \cdot \dot{x_0}\right) \\
& + \frac{1}{2} \left(\Omega_0 \cdot J_0 \Omega_0\right) \\
& + \frac{1}{2} m_{Q3} \left(\left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right) \cdot \left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_3 \cdot J_3 \Omega_3\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& m_0 g \left(e_3 \cdot x_0\right) \\
& + m_1 g \left(e_3 \cdot \left(x_0 + l_1 q_1\right)\right) \\
& + m_2 g \left(e_3 \cdot \left(x_0 + l_1 q_1 + l_2 q_2\right)\right) \\
& + m_3 g \left(e_3 \cdot \left(x_0 + l_1 q_1 + l_2 q_2 + l_3 q_3\right)\right) \\
& + m_{Q0} g \left(e_3 \cdot x_0\right) \\
& + m_{Q3} g \left(e_3 \cdot \left(x_0 + l_1 q_1 + l_2 q_2 + l_3 q_3\right)\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m_0 \left(\dot{x_0} \cdot \dot{x_0}\right) \\
& -m_0 g \left(e_3 \cdot x_0\right) \\
& + \frac{1}{2} m_1 \left(\left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right)\right) \cdot \left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_1 g \left(e_3 \cdot \left(x_0 + l_1 q_1\right)\right) \\
& + \frac{1}{2} m_2 \left(\left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right) \cdot \left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_2 g \left(e_3 \cdot \left(x_0 + l_1 q_1 + l_2 q_2\right)\right) \\
& + \frac{1}{2} m_3 \left(\left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right) \cdot \left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_3 g \left(e_3 \cdot \left(x_0 + l_1 q_1 + l_2 q_2 + l_3 q_3\right)\right) \\
& + \frac{1}{2} m_{Q0} \left(\dot{x_0} \cdot \dot{x_0}\right) \\
& + \frac{1}{2} \left(\Omega_0 \cdot J_0 \Omega_0\right) \\
& -m_{Q0} g \left(e_3 \cdot x_0\right) \\
& + \frac{1}{2} m_{Q3} \left(\left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right) \cdot \left(\dot{x_0} + l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right)\right) \\
& + \frac{1}{2} \left(\Omega_3 \cdot J_3 \Omega_3\right) \\
& -m_{Q3} g \left(e_3 \cdot \left(x_0 + l_1 q_1 + l_2 q_2 + l_3 q_3\right)\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \delta x_0 \cdot f_0 R_0 e_3 \\
& + \eta_0 \cdot M_0 \\
& + \left(\delta x_0 + l_1 \left(\xi_1 \times q_1\right) + l_2 \left(\xi_2 \times q_2\right) + l_3 \left(\xi_3 \times q_3\right)\right) \cdot f_3 R_3 e_3 \\
& + \eta_3 \cdot M_3
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); each $R_k$ equation is in the body frame and each
$q_i$ equation holds on the tangent plane at $q_i$.

<!-- model:eom -->
<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\delta x_0\)</span>, the equation of <span class="math notranslate nohighlight">\(x_0\)</span> (32 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_0 \ddot{x_0} \\
& -m_0 g e_3 \\
& -m_1 \ddot{x_0} \\
& -m_1 l_1 \left(\dot{\omega_1} \times q_1\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_1 l_1 q_1 \\
& -m_1 g e_3 \\
& -m_2 \ddot{x_0} \\
& -m_2 l_1 \left(\dot{\omega_1} \times q_1\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_2 l_1 q_1 \\
& -m_2 l_2 \left(\dot{\omega_2} \times q_2\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_2 l_2 q_2 \\
& -m_2 g e_3 \\
& -m_3 \ddot{x_0} \\
& -m_3 l_1 \left(\dot{\omega_1} \times q_1\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_3 l_1 q_1 \\
& -m_3 l_2 \left(\dot{\omega_2} \times q_2\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_3 l_2 q_2 \\
& -m_3 l_3 \left(\dot{\omega_3} \times q_3\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_3 l_3 q_3 \\
& -m_3 g e_3 \\
& -m_{Q0} \ddot{x_0} \\
& -m_{Q0} g e_3 \\
& -m_{Q3} \ddot{x_0} \\
& -m_{Q3} l_1 \left(\dot{\omega_1} \times q_1\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_{Q3} l_1 q_1 \\
& -m_{Q3} l_2 \left(\dot{\omega_2} \times q_2\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_{Q3} l_2 q_2 \\
& -m_{Q3} l_3 \left(\dot{\omega_3} \times q_3\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_{Q3} l_3 q_3 \\
& -m_{Q3} g e_3 \\
& + f_0 R_0 e_3 \\
& + f_3 R_3 e_3
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_1\)</span>, the equation of <span class="math notranslate nohighlight">\(q_1\)</span> (50 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_1 l_1 \left(\left(\omega_1 \times q_1\right) \times \dot{x_0}\right) \\
& -m_1 l_1 \left(q_1 \times \ddot{x_0}\right) \\
& -m_1 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + 2 m_1 l_1 \left(\omega_1 \times \left(q_1 \times \dot{x_0}\right)\right) \\
& -m_1 l_1 \left(q_1 \times \left(\omega_1 \times \dot{x_0}\right)\right) \\
& -l_1 m_1 g \left(q_1 \times e_3\right) \\
& -m_2 l_1 \left(\left(\omega_1 \times q_1\right) \times \dot{x_0}\right) \\
& -m_2 l_1 \left(q_1 \times \ddot{x_0}\right) \\
& -m_2 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -m_2 l_1 l_2 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_2 l_1 l_2 \left(q_1 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_2 l_1 l_2 \left(q_1 \times q_2\right) \\
& + 2 m_2 l_1 \left(\omega_1 \times \left(q_1 \times \dot{x_0}\right)\right) \\
& + 2 m_2 l_1 l_2 \left(\omega_1 \times \left(q_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_2 l_1 \left(q_1 \times \left(\omega_1 \times \dot{x_0}\right)\right) \\
& -m_2 l_1 l_2 \left(q_1 \times \left(\omega_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -l_1 m_2 g \left(q_1 \times e_3\right) \\
& -m_3 l_1 \left(\left(\omega_1 \times q_1\right) \times \dot{x_0}\right) \\
& -m_3 l_1 \left(q_1 \times \ddot{x_0}\right) \\
& -m_3 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -m_3 l_1 l_2 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_3 l_1 l_2 \left(q_1 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_3 l_1 l_2 \left(q_1 \times q_2\right) \\
& -m_3 l_1 l_3 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_3 l_1 l_3 \left(q_1 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_3 l_1 l_3 \left(q_1 \times q_3\right) \\
& + 2 m_3 l_1 \left(\omega_1 \times \left(q_1 \times \dot{x_0}\right)\right) \\
& + 2 m_3 l_1 l_2 \left(\omega_1 \times \left(q_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_3 l_1 l_3 \left(\omega_1 \times \left(q_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_3 l_1 \left(q_1 \times \left(\omega_1 \times \dot{x_0}\right)\right) \\
& -m_3 l_1 l_2 \left(q_1 \times \left(\omega_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_3 l_1 l_3 \left(q_1 \times \left(\omega_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -l_1 m_3 g \left(q_1 \times e_3\right) \\
& -m_{Q3} l_1 \left(\left(\omega_1 \times q_1\right) \times \dot{x_0}\right) \\
& -m_{Q3} l_1 \left(q_1 \times \ddot{x_0}\right) \\
& -m_{Q3} l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -m_{Q3} l_1 l_2 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_{Q3} l_1 l_2 \left(q_1 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_{Q3} l_1 l_2 \left(q_1 \times q_2\right) \\
& -m_{Q3} l_1 l_3 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_{Q3} l_1 l_3 \left(q_1 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_{Q3} l_1 l_3 \left(q_1 \times q_3\right) \\
& + 2 m_{Q3} l_1 \left(\omega_1 \times \left(q_1 \times \dot{x_0}\right)\right) \\
& + 2 m_{Q3} l_1 l_2 \left(\omega_1 \times \left(q_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_{Q3} l_1 l_3 \left(\omega_1 \times \left(q_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_{Q3} l_1 \left(q_1 \times \left(\omega_1 \times \dot{x_0}\right)\right) \\
& -m_{Q3} l_1 l_2 \left(q_1 \times \left(\omega_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_{Q3} l_1 l_3 \left(q_1 \times \left(\omega_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -l_1 m_{Q3} g \left(q_1 \times e_3\right) \\
& + f_3 l_1 \left(q_1 \times \left(R_3 e_3\right)\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_2\)</span>, the equation of <span class="math notranslate nohighlight">\(q_2\)</span> (44 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_2 l_2 \left(\left(\omega_2 \times q_2\right) \times \dot{x_0}\right) \\
& -m_2 l_2 \left(q_2 \times \ddot{x_0}\right) \\
& -m_2 l_2 l_1 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_2 l_2 l_1 \left(q_2 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_2 l_2 l_1 \left(q_2 \times q_1\right) \\
& -m_2 l_2 l_2 \left(q_2 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + 2 m_2 l_2 \left(\omega_2 \times \left(q_2 \times \dot{x_0}\right)\right) \\
& + 2 m_2 l_2 l_1 \left(\omega_2 \times \left(q_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_2 l_2 \left(q_2 \times \left(\omega_2 \times \dot{x_0}\right)\right) \\
& -m_2 l_2 l_1 \left(q_2 \times \left(\omega_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -l_2 m_2 g \left(q_2 \times e_3\right) \\
& -m_3 l_2 \left(\left(\omega_2 \times q_2\right) \times \dot{x_0}\right) \\
& -m_3 l_2 \left(q_2 \times \ddot{x_0}\right) \\
& -m_3 l_2 l_1 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_3 l_2 l_1 \left(q_2 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_3 l_2 l_1 \left(q_2 \times q_1\right) \\
& -m_3 l_2 l_2 \left(q_2 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& -m_3 l_2 l_3 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_3 l_2 l_3 \left(q_2 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_3 l_2 l_3 \left(q_2 \times q_3\right) \\
& + 2 m_3 l_2 \left(\omega_2 \times \left(q_2 \times \dot{x_0}\right)\right) \\
& + 2 m_3 l_2 l_1 \left(\omega_2 \times \left(q_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_3 l_2 l_3 \left(\omega_2 \times \left(q_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_3 l_2 \left(q_2 \times \left(\omega_2 \times \dot{x_0}\right)\right) \\
& -m_3 l_2 l_1 \left(q_2 \times \left(\omega_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_3 l_2 l_3 \left(q_2 \times \left(\omega_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -l_2 m_3 g \left(q_2 \times e_3\right) \\
& -m_{Q3} l_2 \left(\left(\omega_2 \times q_2\right) \times \dot{x_0}\right) \\
& -m_{Q3} l_2 \left(q_2 \times \ddot{x_0}\right) \\
& -m_{Q3} l_2 l_1 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_{Q3} l_2 l_1 \left(q_2 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_{Q3} l_2 l_1 \left(q_2 \times q_1\right) \\
& -m_{Q3} l_2 l_2 \left(q_2 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& -m_{Q3} l_2 l_3 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_{Q3} l_2 l_3 \left(q_2 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_{Q3} l_2 l_3 \left(q_2 \times q_3\right) \\
& + 2 m_{Q3} l_2 \left(\omega_2 \times \left(q_2 \times \dot{x_0}\right)\right) \\
& + 2 m_{Q3} l_2 l_1 \left(\omega_2 \times \left(q_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_{Q3} l_2 l_3 \left(\omega_2 \times \left(q_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_{Q3} l_2 \left(q_2 \times \left(\omega_2 \times \dot{x_0}\right)\right) \\
& -m_{Q3} l_2 l_1 \left(q_2 \times \left(\omega_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_{Q3} l_2 l_3 \left(q_2 \times \left(\omega_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -l_2 m_{Q3} g \left(q_2 \times e_3\right) \\
& + f_3 l_2 \left(q_2 \times \left(R_3 e_3\right)\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_3\)</span>, the equation of <span class="math notranslate nohighlight">\(q_3\)</span> (33 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_3 l_3 \left(\left(\omega_3 \times q_3\right) \times \dot{x_0}\right) \\
& -m_3 l_3 \left(q_3 \times \ddot{x_0}\right) \\
& -m_3 l_3 l_1 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_3 l_3 l_1 \left(q_3 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_3 l_3 l_1 \left(q_3 \times q_1\right) \\
& -m_3 l_3 l_2 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_3 l_3 l_2 \left(q_3 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_3 l_3 l_2 \left(q_3 \times q_2\right) \\
& -m_3 l_3 l_3 \left(q_3 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + 2 m_3 l_3 \left(\omega_3 \times \left(q_3 \times \dot{x_0}\right)\right) \\
& + 2 m_3 l_3 l_1 \left(\omega_3 \times \left(q_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_3 l_3 l_2 \left(\omega_3 \times \left(q_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_3 l_3 \left(q_3 \times \left(\omega_3 \times \dot{x_0}\right)\right) \\
& -m_3 l_3 l_1 \left(q_3 \times \left(\omega_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_3 l_3 l_2 \left(q_3 \times \left(\omega_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -l_3 m_3 g \left(q_3 \times e_3\right) \\
& -m_{Q3} l_3 \left(\left(\omega_3 \times q_3\right) \times \dot{x_0}\right) \\
& -m_{Q3} l_3 \left(q_3 \times \ddot{x_0}\right) \\
& -m_{Q3} l_3 l_1 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_{Q3} l_3 l_1 \left(q_3 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_{Q3} l_3 l_1 \left(q_3 \times q_1\right) \\
& -m_{Q3} l_3 l_2 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_{Q3} l_3 l_2 \left(q_3 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_{Q3} l_3 l_2 \left(q_3 \times q_2\right) \\
& -m_{Q3} l_3 l_3 \left(q_3 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + 2 m_{Q3} l_3 \left(\omega_3 \times \left(q_3 \times \dot{x_0}\right)\right) \\
& + 2 m_{Q3} l_3 l_1 \left(\omega_3 \times \left(q_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_{Q3} l_3 l_2 \left(\omega_3 \times \left(q_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_{Q3} l_3 \left(q_3 \times \left(\omega_3 \times \dot{x_0}\right)\right) \\
& -m_{Q3} l_3 l_1 \left(q_3 \times \left(\omega_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_{Q3} l_3 l_2 \left(q_3 \times \left(\omega_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -l_3 m_{Q3} g \left(q_3 \times e_3\right) \\
& + f_3 l_3 \left(q_3 \times \left(R_3 e_3\right)\right)
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_0\)</span>, the equation of <span class="math notranslate nohighlight">\(R_0\)</span> (3 terms)</summary>

$$
\begin{aligned}
0 = {}& -J_0 \dot{\Omega_0} \\
& -\Omega_0 \times \left(J_0 \Omega_0\right) \\
& + M_0
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
<!-- /model:eom -->

## Code

```python
from geomech.systems import get

# Joints x_k = x_0 + l_1 q_1 + ... + l_k q_k, k = 0..n_links, with point masses m_k;
# quadrotors at x_0 and x_n.
# L = Σ_k ( ½ m_k ẋ_k·ẋ_k − m_k g e3·x_k )
#     + Σ_{k=0,n} ( ½ m_Qk ẋ_k·ẋ_k + ½ Ω_k·J_kΩ_k − m_Qk g e3·x_k )
# δW = Σ_{k=0,n} ( f_k (R_k e3)·δx_k + M_k·η_k )
system = get("multi_quadrotor_hose")  # n_links=3
eom = system.derive()
```

The same model is `geomech.systems.get("multi_quadrotor_hose")`, with default parameters,
a default initial state and the 3-D scene used above.
`geomech.systems.get("multi_quadrotor_hose", n_links=n)` builds the hose with $n$ links (1
to 20); `geomech.systems.aerial.hose_system(n_links, attach)` puts quadrotors at any set
of joints.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-5}$ |
| Energy $H = T + V$, 6 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-9}$, with $H(0) = 65.8$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
