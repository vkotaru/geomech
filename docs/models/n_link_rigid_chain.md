# Chain of rigid bodies

A chain of $n$ rigid links under gravity: link 1 hangs from a fixed pivot by a spherical
joint, and each further link hangs from the lower end of the one above by a spherical
joint. The configuration is $(R_1, \dots, R_n) \in SO(3)^n$, with body angular velocities
$\Omega_k$ and $\dot R_k = R_k \widehat{\Omega_k}$; no Euler angles, so no singular
attitude. The inputs are torques $\tau_k$ at the joints, in the inertial frame. The
equations below are derived once for every $n$ with `geomech.indexed`, then for $n = 5$
with the library system, which the animation uses. Not taken from the literature; the
equations are verified numerically (`check_eom`). Library system `n_link_rigid_chain`.

## Animation

Ten seconds, unforced, five links, from the library's initial state: the links tilted
from 32° (link 1) to 72° (link 5) off the downward vertical, spinning at 2 rad/s about
their long axes in alternating directions, with small swinging rates. The chain swings
round the vertical in irregular loops; link 1 stays within 38° of the downward vertical
and link 5 within 80°. Each link is 0.8 times as long and as heavy as the one above. Drag to orbit; the full report (equations, standard form,
parameters) opens in a
<a href="../model-reports/n_link_rigid_chain.html">separate page</a>.

<iframe class="model-report" src="../model-reports/n_link_rigid_chain.html?view" title="Chain of rigid bodies, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 300" role="img" aria-labelledby="nrc-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="nrc-title">Chain of three rigid links: offsets R1 rho1 and R1 l1, angular velocities Omega1 to Omega3, masses and inertias, joint torques tau1 to tau3 and gravity</title>
  <defs>
    <marker id="nrc-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <line x1="155.0" y1="42.0" x2="215.0" y2="42.0" stroke-width="2"/>
    <path d="M157,42.0 l-8,-10 M169,42.0 l-8,-10 M181,42.0 l-8,-10 M193,42.0 l-8,-10 M205,42.0 l-8,-10 M217,42.0 l-8,-10"/>
    <line x1="40" y1="282" x2="40" y2="244" marker-end="url(#nrc-arr)"/>
    <line x1="40" y1="282" x2="78" y2="282" marker-end="url(#nrc-arr)"/>
    <line x1="185.0" y1="42.0" x2="185.0" y2="272.0" stroke-dasharray="3 4" opacity="0.45"/>
    <path d="M179.2,45.9 L226.2,115.6 L237.8,107.7 L190.8,38.1 Z" stroke-width="1.8"/>
    <line x1="171.7" y1="50.9" x2="195.2" y2="85.8" stroke-width="1.6" marker-end="url(#nrc-arr)"/>
    <line x1="156.8" y1="61.0" x2="203.8" y2="130.7" marker-end="url(#nrc-arr)"/>
    <path d="M225.9,108.1 L183.9,180.9 L196.0,187.9 L238.0,115.1 Z" stroke-width="1.8"/>
    <path d="M184.2,188.3 L231.1,257.9 L242.7,250.1 L195.8,180.5 Z" stroke-width="1.8"/>
  </g>
  <path d="M210.8,101.6 A 13 5 -34 0 0 232.9,91.5" fill="none" stroke="#3f8fa0" stroke-width="1.8" marker-end="url(#nrc-arr)"/>
  <path d="M200.0,39.1 Q212.1,35.8 209.4,49.6" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#nrc-arr)"/>
  <path d="M189.7,160.9 A 13 5 30 0 0 208.4,176.4" fill="none" stroke="#3f8fa0" stroke-width="1.8" marker-end="url(#nrc-arr)"/>
  <path d="M241.1,123.8 Q249.5,133.3 235.8,137.0" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#nrc-arr)"/>
  <path d="M215.7,244.0 A 13 5 -34 0 0 237.9,233.9" fill="none" stroke="#3f8fa0" stroke-width="1.8" marker-end="url(#nrc-arr)"/>
  <path d="M181.7,197.2 Q174.0,207.2 187.9,209.9" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#nrc-arr)"/>
  <line x1="340" y1="40" x2="340" y2="90" stroke="#e08a2c" stroke-width="2" marker-end="url(#nrc-arr)"/>
  <circle cx="185.0" cy="42.0" r="4" fill="currentColor"/>
  <circle cx="208.5" cy="76.8" r="3" fill="currentColor"/>
  <circle cx="232.0" cy="111.6" r="4" fill="currentColor"/>
  <circle cx="211.0" cy="148.0" r="3" fill="currentColor"/>
  <circle cx="190.0" cy="184.4" r="4" fill="currentColor"/>
  <circle cx="213.5" cy="219.2" r="3" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="141" y="76">R<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan>ρ<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan></text>
    <text x="146" y="119">R<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan>l<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan></text>
    <text x="239" y="89" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan></text>
    <text x="220" y="67">m<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan>, J<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan></text>
    <text x="214" y="34" fill="#e08a2c">τ<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan></text>
    <text x="219" y="184" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">2</tspan><tspan dy="-4"></tspan></text>
    <text x="230" y="157">m<tspan font-size="12" dy="4">2</tspan><tspan dy="-4"></tspan>, J<tspan font-size="12" dy="4">2</tspan><tspan dy="-4"></tspan></text>
    <text x="254" y="136" fill="#e08a2c">τ<tspan font-size="12" dy="4">2</tspan><tspan dy="-4"></tspan></text>
    <text x="244" y="231" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">3</tspan><tspan dy="-4"></tspan></text>
    <text x="225" y="210">m<tspan font-size="12" dy="4">3</tspan><tspan dy="-4"></tspan>, J<tspan font-size="12" dy="4">3</tspan><tspan dy="-4"></tspan></text>
    <text x="156" y="210" fill="#e08a2c">τ<tspan font-size="12" dy="4">3</tspan><tspan dy="-4"></tspan></text>
    <text x="348" y="80" fill="#e08a2c">g</text>
    <text x="30" y="238">e<tspan font-size="12" dy="4">3</tspan><tspan dy="-4"></tspan></text>
    <text x="84" y="287">e<tspan font-size="12" dy="4">1</tspan><tspan dy="-4"></tspan></text>
  </g>
</svg>
<figcaption>Drawn with three links. Link <span class="math notranslate nohighlight">\(k\)</span> hangs from its upper joint (dots);
its centre of mass is at <span class="math notranslate nohighlight">\(R_k \rho_k\)</span> and its lower joint at
<span class="math notranslate nohighlight">\(R_k l_k\)</span> from that joint (drawn for link 1). The joint torque
<span class="math notranslate nohighlight">\(\tau_k\)</span> acts on link <span class="math notranslate nohighlight">\(k\)</span> and <span class="math notranslate nohighlight">\(-\tau_k\)</span> on link
<span class="math notranslate nohighlight">\(k - 1\)</span> (<span class="math notranslate nohighlight">\(\tau_1\)</span> between the ground and link 1). Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span>.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $R_k$ | attitude of link $k$, body frame to inertial frame | $SO(3)$ each |
| $\Omega_k$ | body angular velocity of link $k$, $\dot R_k = R_k \widehat{\Omega_k}$ | $\mathbb{R}^3$ each |
| $\eta_k$ | variation: $\delta R_k = R_k \widehat{\eta_k}$ | $\mathbb{R}^3$ each |
| $\tau_k$ | input torque at joint $k$, inertial frame | $\mathbb{R}^3$ each |
| $M_k$ | body moment on link $k$ in the generic model; $M_k = R_k^T(\tau_k - \tau_{k+1})$, $\tau_{n+1} = 0$ | $\mathbb{R}^3$ each |
| $J_k$ | inertia of link $k$ about its centre of mass, body frame | constant, symmetric |
| $\rho_k$ | upper joint to centre of mass of link $k$, body frame | constant |
| $l_k$ | upper joint to lower joint of link $k$, body frame ($k < n$) | constant |
| $m_k, g, e_3$ | masses, gravity, vertical unit vector | constants |

## Lagrangian and equations for every $n$

<!-- model:indexed -->
Lagrangian for every $n$, as `geomech.indexed.systems.rigid_chain()` holds it (the velocities expanded, so the inner sums are the inertia couplings):

$$
\begin{aligned}
L = {}& \sum_{j=1}^{n} \frac{1}{2} \left(\Omega_{j} \cdot J_{j} \Omega_{j}\right) \\
& + \sum_{j=1}^{n} \sum_{k=1}^{n} \sum_{r=\max(j+1, k+1)}^{n} \frac{1}{2} m_{r} \left(R_{j} \widehat{\Omega_{j}} l_{j} \cdot R_{k} \widehat{\Omega_{k}} l_{k}\right) \\
& + \sum_{j=1}^{n} \sum_{k=1}^{j-1} \frac{1}{2} m_{j} \left(R_{k} \widehat{\Omega_{k}} l_{k} \cdot R_{j} \widehat{\Omega_{j}} \rho_{j}\right) \\
& + \sum_{j=1}^{n} \sum_{k=1}^{j-1} \frac{1}{2} m_{j} \left(R_{j} \widehat{\Omega_{j}} \rho_{j} \cdot R_{k} \widehat{\Omega_{k}} l_{k}\right) \\
& + \sum_{j=1}^{n} \frac{1}{2} m_{j} \left(R_{j} \widehat{\Omega_{j}} \rho_{j} \cdot R_{j} \widehat{\Omega_{j}} \rho_{j}\right) \\
& -\sum_{j=1}^{n} \sum_{k=j+1}^{n} m_{k} g \left(R_{j} l_{j} \cdot e_3\right) \\
& -\sum_{j=1}^{n} m_{j} g \left(R_{j} \rho_{j} \cdot e_3\right)
\end{aligned}
$$

Virtual work:

$$
\begin{aligned}
\delta W = {}& \sum_{j=1}^{n} \eta_{j} \cdot M_{j}
\end{aligned}
$$

Equations from `compute_indexed_eom`, derived once for every $n$; $i$ is the body, $1 \le i \le n$:

Coefficient of $\eta_{i}$, the equation of $R_i$:

$$
\begin{aligned}
0 = {}& -J_{i} \dot{\Omega_{i}} \\
& -\Omega_{i} \times \left(J_{i} \Omega_{i}\right) \\
& -m_{i} \left(\rho_{i} \times \left(\dot{\Omega_{i}} \times \rho_{i}\right)\right) \\
& -m_{i} \left(\Omega_{i} \times \left(\rho_{i} \times \left(\Omega_{i} \times \rho_{i}\right)\right)\right) \\
& -m_{i} g \left(\rho_{i} \times \left(R_{i}^T e_3\right)\right) \\
& + M_{i} \\
& + \sum_{j=i+1}^{n} m_{j} \left(\left(\Omega_{i} \times l_{i}\right) \times \left(R_{i}^T R_{j} \left(\Omega_{j} \times \rho_{j}\right)\right)\right) \\
& + \sum_{j=i+1}^{n} m_{j} \left(l_{i} \times \left(\Omega_{i} \times \left(R_{i}^T R_{j} \left(\Omega_{j} \times \rho_{j}\right)\right)\right)\right) \\
& -\sum_{j=i+1}^{n} m_{j} \left(l_{i} \times \left(R_{i}^T R_{j} \left(\Omega_{j} \times \left(\Omega_{j} \times \rho_{j}\right)\right)\right)\right) \\
& -\sum_{j=i+1}^{n} m_{j} \left(l_{i} \times \left(R_{i}^T R_{j} \left(\dot{\Omega_{j}} \times \rho_{j}\right)\right)\right) \\
& -\sum_{j=i+1}^{n} m_{j} \left(\Omega_{i} \times \left(l_{i} \times \left(R_{i}^T R_{j} \left(\Omega_{j} \times \rho_{j}\right)\right)\right)\right) \\
& -\sum_{j=i+1}^{n} m_{j} g \left(l_{i} \times \left(R_{i}^T e_3\right)\right) \\
& + \sum_{j=1}^{i-1} m_{i} \left(\left(\Omega_{i} \times \rho_{i}\right) \times \left(R_{i}^T R_{j} \left(\Omega_{j} \times l_{j}\right)\right)\right) \\
& + \sum_{j=1}^{i-1} m_{i} \left(\rho_{i} \times \left(\Omega_{i} \times \left(R_{i}^T R_{j} \left(\Omega_{j} \times l_{j}\right)\right)\right)\right) \\
& -\sum_{j=1}^{i-1} m_{i} \left(\rho_{i} \times \left(R_{i}^T R_{j} \left(\Omega_{j} \times \left(\Omega_{j} \times l_{j}\right)\right)\right)\right) \\
& -\sum_{j=1}^{i-1} m_{i} \left(\rho_{i} \times \left(R_{i}^T R_{j} \left(\dot{\Omega_{j}} \times l_{j}\right)\right)\right) \\
& -\sum_{j=1}^{i-1} m_{i} \left(\Omega_{i} \times \left(\rho_{i} \times \left(R_{i}^T R_{j} \left(\Omega_{j} \times l_{j}\right)\right)\right)\right) \\
& + \sum_{j=1}^{n} \sum_{k=\max(i+1, j+1)}^{n} m_{k} \left(\left(\Omega_{i} \times l_{i}\right) \times \left(R_{i}^T R_{j} \left(\Omega_{j} \times l_{j}\right)\right)\right) \\
& + \sum_{j=1}^{n} \sum_{k=\max(i+1, j+1)}^{n} m_{k} \left(l_{i} \times \left(\Omega_{i} \times \left(R_{i}^T R_{j} \left(\Omega_{j} \times l_{j}\right)\right)\right)\right) \\
& -\sum_{j=1}^{n} \sum_{k=\max(i+1, j+1)}^{n} m_{k} \left(l_{i} \times \left(R_{i}^T R_{j} \left(\Omega_{j} \times \left(\Omega_{j} \times l_{j}\right)\right)\right)\right) \\
& -\sum_{j=1}^{n} \sum_{k=\max(i+1, j+1)}^{n} m_{k} \left(l_{i} \times \left(R_{i}^T R_{j} \left(\dot{\Omega_{j}} \times l_{j}\right)\right)\right) \\
& -\sum_{j=1}^{n} \sum_{k=\max(i+1, j+1)}^{n} m_{k} \left(\Omega_{i} \times \left(l_{i} \times \left(R_{i}^T R_{j} \left(\Omega_{j} \times l_{j}\right)\right)\right)\right)
\end{aligned}
$$
<!-- /model:indexed -->

## Energy, Lagrangian and virtual work at $n = 5$

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} \left(\Omega_1 \cdot J_1 \Omega_1\right) \\
& + \frac{1}{2} \left(\Omega_2 \cdot J_2 \Omega_2\right) \\
& + \frac{1}{2} \left(\Omega_3 \cdot J_3 \Omega_3\right) \\
& + \frac{1}{2} \left(\Omega_4 \cdot J_4 \Omega_4\right) \\
& + \frac{1}{2} \left(\Omega_5 \cdot J_5 \Omega_5\right) \\
& + \frac{1}{2} \left(R_1 \widehat{\Omega_1} \rho_1 \cdot R_1 \widehat{\Omega_1} \rho_1\right) m_1 \\
& + \frac{1}{2} \left(\left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} \rho_2\right) \cdot \left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} \rho_2\right)\right) m_2 \\
& + \frac{1}{2} \left(\left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} l_2 + R_3 \widehat{\Omega_3} \rho_3\right) \cdot \left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} l_2 + R_3 \widehat{\Omega_3} \rho_3\right)\right) m_3 \\
& + \frac{1}{2} \left(\left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} l_2 + R_3 \widehat{\Omega_3} l_3 + R_4 \widehat{\Omega_4} \rho_4\right) \cdot \left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} l_2 + R_3 \widehat{\Omega_3} l_3 + R_4 \widehat{\Omega_4} \rho_4\right)\right) m_4 \\
& + \frac{1}{2} \left(\left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} l_2 + R_3 \widehat{\Omega_3} l_3 + R_4 \widehat{\Omega_4} l_4 + R_5 \widehat{\Omega_5} \rho_5\right) \cdot \left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} l_2 + R_3 \widehat{\Omega_3} l_3 + R_4 \widehat{\Omega_4} l_4 + R_5 \widehat{\Omega_5} \rho_5\right)\right) m_5
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& \left(m_1 g \left(R_1 \rho_1 \cdot e_3\right) + m_2 g \left(\left(R_1 l_1 + R_2 \rho_2\right) \cdot e_3\right) + m_3 g \left(\left(R_1 l_1 + R_2 l_2 + R_3 \rho_3\right) \cdot e_3\right) + m_4 g \left(\left(R_1 l_1 + R_2 l_2 + R_3 l_3 + R_4 \rho_4\right) \cdot e_3\right) + m_5 g \left(\left(R_1 l_1 + R_2 l_2 + R_3 l_3 + R_4 l_4 + R_5 \rho_5\right) \cdot e_3\right)\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} \left(\Omega_1 \cdot J_1 \Omega_1\right) \\
& + \frac{1}{2} \left(\Omega_2 \cdot J_2 \Omega_2\right) \\
& + \frac{1}{2} \left(\Omega_3 \cdot J_3 \Omega_3\right) \\
& + \frac{1}{2} \left(\Omega_4 \cdot J_4 \Omega_4\right) \\
& + \frac{1}{2} \left(\Omega_5 \cdot J_5 \Omega_5\right) \\
& + \frac{1}{2} \left(R_1 \widehat{\Omega_1} \rho_1 \cdot R_1 \widehat{\Omega_1} \rho_1\right) m_1 \\
& + \frac{1}{2} \left(\left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} \rho_2\right) \cdot \left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} \rho_2\right)\right) m_2 \\
& + \frac{1}{2} \left(\left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} l_2 + R_3 \widehat{\Omega_3} \rho_3\right) \cdot \left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} l_2 + R_3 \widehat{\Omega_3} \rho_3\right)\right) m_3 \\
& + \frac{1}{2} \left(\left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} l_2 + R_3 \widehat{\Omega_3} l_3 + R_4 \widehat{\Omega_4} \rho_4\right) \cdot \left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} l_2 + R_3 \widehat{\Omega_3} l_3 + R_4 \widehat{\Omega_4} \rho_4\right)\right) m_4 \\
& + \frac{1}{2} \left(\left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} l_2 + R_3 \widehat{\Omega_3} l_3 + R_4 \widehat{\Omega_4} l_4 + R_5 \widehat{\Omega_5} \rho_5\right) \cdot \left(R_1 \widehat{\Omega_1} l_1 + R_2 \widehat{\Omega_2} l_2 + R_3 \widehat{\Omega_3} l_3 + R_4 \widehat{\Omega_4} l_4 + R_5 \widehat{\Omega_5} \rho_5\right)\right) m_5 \\
& -\left(m_1 g \left(R_1 \rho_1 \cdot e_3\right) + m_2 g \left(\left(R_1 l_1 + R_2 \rho_2\right) \cdot e_3\right) + m_3 g \left(\left(R_1 l_1 + R_2 l_2 + R_3 \rho_3\right) \cdot e_3\right) + m_4 g \left(\left(R_1 l_1 + R_2 l_2 + R_3 l_3 + R_4 \rho_4\right) \cdot e_3\right) + m_5 g \left(\left(R_1 l_1 + R_2 l_2 + R_3 l_3 + R_4 l_4 + R_5 \rho_5\right) \cdot e_3\right)\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \eta_1 \cdot R_1^T \tau_1 \\
& + \eta_2 \cdot R_2^T \tau_2 \\
& -\eta_1 \cdot R_1^T \tau_2 \\
& + \eta_3 \cdot R_3^T \tau_3 \\
& -\eta_2 \cdot R_2^T \tau_3 \\
& + \eta_4 \cdot R_4^T \tau_4 \\
& -\eta_3 \cdot R_3^T \tau_4 \\
& + \eta_5 \cdot R_5^T \tau_5 \\
& -\eta_4 \cdot R_4^T \tau_5
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion at $n = 5$

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); each is written in the body frame of its link.

<!-- model:eom -->
<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_1\)</span>, the equation of <span class="math notranslate nohighlight">\(R_1\)</span> (69 terms)</summary>

$$
\begin{aligned}
0 = {}& -J_1 \dot{\Omega_1} \\
& -\Omega_1 \times \left(J_1 \Omega_1\right) \\
& -m_1 \left(\rho_1 \times \left(\dot{\Omega_1} \times \rho_1\right)\right) \\
& -m_1 \left(\Omega_1 \times \left(\rho_1 \times \left(\Omega_1 \times \rho_1\right)\right)\right) \\
& + m_2 \left(\left(\Omega_1 \times l_1\right) \times \left(R_1^T R_2 \left(\Omega_2 \times \rho_2\right)\right)\right) \\
& -m_2 \left(l_1 \times \left(\dot{\Omega_1} \times l_1\right)\right) \\
& + m_2 \left(l_1 \times \left(\Omega_1 \times \left(R_1^T R_2 \left(\Omega_2 \times \rho_2\right)\right)\right)\right) \\
& -m_2 \left(l_1 \times \left(R_1^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times \rho_2\right)\right)\right)\right) \\
& -m_2 \left(l_1 \times \left(R_1^T R_2 \left(\dot{\Omega_2} \times \rho_2\right)\right)\right) \\
& -m_2 \left(\Omega_1 \times \left(l_1 \times \left(\Omega_1 \times l_1\right)\right)\right) \\
& -m_2 \left(\Omega_1 \times \left(l_1 \times \left(R_1^T R_2 \left(\Omega_2 \times \rho_2\right)\right)\right)\right) \\
& + m_3 \left(\left(\Omega_1 \times l_1\right) \times \left(R_1^T R_2 \left(\Omega_2 \times l_2\right)\right)\right) \\
& + m_3 \left(\left(\Omega_1 \times l_1\right) \times \left(R_1^T R_3 \left(\Omega_3 \times \rho_3\right)\right)\right) \\
& -m_3 \left(l_1 \times \left(\dot{\Omega_1} \times l_1\right)\right) \\
& + m_3 \left(l_1 \times \left(\Omega_1 \times \left(R_1^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_3 \left(l_1 \times \left(R_1^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_3 \left(l_1 \times \left(R_1^T R_2 \left(\dot{\Omega_2} \times l_2\right)\right)\right) \\
& + m_3 \left(l_1 \times \left(\Omega_1 \times \left(R_1^T R_3 \left(\Omega_3 \times \rho_3\right)\right)\right)\right) \\
& -m_3 \left(l_1 \times \left(R_1^T R_3 \left(\Omega_3 \times \left(\Omega_3 \times \rho_3\right)\right)\right)\right) \\
& -m_3 \left(l_1 \times \left(R_1^T R_3 \left(\dot{\Omega_3} \times \rho_3\right)\right)\right) \\
& -m_3 \left(\Omega_1 \times \left(l_1 \times \left(\Omega_1 \times l_1\right)\right)\right) \\
& -m_3 \left(\Omega_1 \times \left(l_1 \times \left(R_1^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_3 \left(\Omega_1 \times \left(l_1 \times \left(R_1^T R_3 \left(\Omega_3 \times \rho_3\right)\right)\right)\right) \\
& + m_4 \left(\left(\Omega_1 \times l_1\right) \times \left(R_1^T R_2 \left(\Omega_2 \times l_2\right)\right)\right) \\
& + m_4 \left(\left(\Omega_1 \times l_1\right) \times \left(R_1^T R_3 \left(\Omega_3 \times l_3\right)\right)\right) \\
& + m_4 \left(\left(\Omega_1 \times l_1\right) \times \left(R_1^T R_4 \left(\Omega_4 \times \rho_4\right)\right)\right) \\
& -m_4 \left(l_1 \times \left(\dot{\Omega_1} \times l_1\right)\right) \\
& + m_4 \left(l_1 \times \left(\Omega_1 \times \left(R_1^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_4 \left(l_1 \times \left(R_1^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_4 \left(l_1 \times \left(R_1^T R_2 \left(\dot{\Omega_2} \times l_2\right)\right)\right) \\
& + m_4 \left(l_1 \times \left(\Omega_1 \times \left(R_1^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_4 \left(l_1 \times \left(R_1^T R_3 \left(\Omega_3 \times \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_4 \left(l_1 \times \left(R_1^T R_3 \left(\dot{\Omega_3} \times l_3\right)\right)\right) \\
& + m_4 \left(l_1 \times \left(\Omega_1 \times \left(R_1^T R_4 \left(\Omega_4 \times \rho_4\right)\right)\right)\right) \\
& -m_4 \left(l_1 \times \left(R_1^T R_4 \left(\Omega_4 \times \left(\Omega_4 \times \rho_4\right)\right)\right)\right) \\
& -m_4 \left(l_1 \times \left(R_1^T R_4 \left(\dot{\Omega_4} \times \rho_4\right)\right)\right) \\
& -m_4 \left(\Omega_1 \times \left(l_1 \times \left(\Omega_1 \times l_1\right)\right)\right) \\
& -m_4 \left(\Omega_1 \times \left(l_1 \times \left(R_1^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_4 \left(\Omega_1 \times \left(l_1 \times \left(R_1^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_4 \left(\Omega_1 \times \left(l_1 \times \left(R_1^T R_4 \left(\Omega_4 \times \rho_4\right)\right)\right)\right) \\
& + m_5 \left(\left(\Omega_1 \times l_1\right) \times \left(R_1^T R_2 \left(\Omega_2 \times l_2\right)\right)\right) \\
& + m_5 \left(\left(\Omega_1 \times l_1\right) \times \left(R_1^T R_3 \left(\Omega_3 \times l_3\right)\right)\right) \\
& + m_5 \left(\left(\Omega_1 \times l_1\right) \times \left(R_1^T R_4 \left(\Omega_4 \times l_4\right)\right)\right) \\
& + m_5 \left(\left(\Omega_1 \times l_1\right) \times \left(R_1^T R_5 \left(\Omega_5 \times \rho_5\right)\right)\right) \\
& -m_5 \left(l_1 \times \left(\dot{\Omega_1} \times l_1\right)\right) \\
& + m_5 \left(l_1 \times \left(\Omega_1 \times \left(R_1^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_5 \left(l_1 \times \left(R_1^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_5 \left(l_1 \times \left(R_1^T R_2 \left(\dot{\Omega_2} \times l_2\right)\right)\right) \\
& + m_5 \left(l_1 \times \left(\Omega_1 \times \left(R_1^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_5 \left(l_1 \times \left(R_1^T R_3 \left(\Omega_3 \times \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_5 \left(l_1 \times \left(R_1^T R_3 \left(\dot{\Omega_3} \times l_3\right)\right)\right) \\
& + m_5 \left(l_1 \times \left(\Omega_1 \times \left(R_1^T R_4 \left(\Omega_4 \times l_4\right)\right)\right)\right) \\
& -m_5 \left(l_1 \times \left(R_1^T R_4 \left(\Omega_4 \times \left(\Omega_4 \times l_4\right)\right)\right)\right) \\
& -m_5 \left(l_1 \times \left(R_1^T R_4 \left(\dot{\Omega_4} \times l_4\right)\right)\right) \\
& + m_5 \left(l_1 \times \left(\Omega_1 \times \left(R_1^T R_5 \left(\Omega_5 \times \rho_5\right)\right)\right)\right) \\
& -m_5 \left(l_1 \times \left(R_1^T R_5 \left(\Omega_5 \times \left(\Omega_5 \times \rho_5\right)\right)\right)\right) \\
& -m_5 \left(l_1 \times \left(R_1^T R_5 \left(\dot{\Omega_5} \times \rho_5\right)\right)\right) \\
& -m_5 \left(\Omega_1 \times \left(l_1 \times \left(\Omega_1 \times l_1\right)\right)\right) \\
& -m_5 \left(\Omega_1 \times \left(l_1 \times \left(R_1^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_5 \left(\Omega_1 \times \left(l_1 \times \left(R_1^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_5 \left(\Omega_1 \times \left(l_1 \times \left(R_1^T R_4 \left(\Omega_4 \times l_4\right)\right)\right)\right) \\
& -m_5 \left(\Omega_1 \times \left(l_1 \times \left(R_1^T R_5 \left(\Omega_5 \times \rho_5\right)\right)\right)\right) \\
& -m_1 g \left(\rho_1 \times \left(R_1^T e_3\right)\right) \\
& -m_2 g \left(l_1 \times \left(R_1^T e_3\right)\right) \\
& -m_3 g \left(l_1 \times \left(R_1^T e_3\right)\right) \\
& -m_4 g \left(l_1 \times \left(R_1^T e_3\right)\right) \\
& -m_5 g \left(l_1 \times \left(R_1^T e_3\right)\right) \\
& + R_1^T \tau_1 \\
& -R_1^T \tau_2
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_2\)</span>, the equation of <span class="math notranslate nohighlight">\(R_2\)</span> (66 terms)</summary>

$$
\begin{aligned}
0 = {}& -J_2 \dot{\Omega_2} \\
& -\Omega_2 \times \left(J_2 \Omega_2\right) \\
& + m_2 \left(\left(\Omega_2 \times \rho_2\right) \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right) \\
& + m_2 \left(\rho_2 \times \left(\Omega_2 \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_2 \left(\rho_2 \times \left(R_2^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_2 \left(\rho_2 \times \left(R_2^T R_1 \left(\dot{\Omega_1} \times l_1\right)\right)\right) \\
& -m_2 \left(\rho_2 \times \left(\dot{\Omega_2} \times \rho_2\right)\right) \\
& -m_2 \left(\Omega_2 \times \left(\rho_2 \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_2 \left(\Omega_2 \times \left(\rho_2 \times \left(\Omega_2 \times \rho_2\right)\right)\right) \\
& + m_3 \left(\left(\Omega_2 \times l_2\right) \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right) \\
& + m_3 \left(\left(\Omega_2 \times l_2\right) \times \left(R_2^T R_3 \left(\Omega_3 \times \rho_3\right)\right)\right) \\
& + m_3 \left(l_2 \times \left(\Omega_2 \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_3 \left(l_2 \times \left(R_2^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_3 \left(l_2 \times \left(R_2^T R_1 \left(\dot{\Omega_1} \times l_1\right)\right)\right) \\
& -m_3 \left(l_2 \times \left(\dot{\Omega_2} \times l_2\right)\right) \\
& + m_3 \left(l_2 \times \left(\Omega_2 \times \left(R_2^T R_3 \left(\Omega_3 \times \rho_3\right)\right)\right)\right) \\
& -m_3 \left(l_2 \times \left(R_2^T R_3 \left(\Omega_3 \times \left(\Omega_3 \times \rho_3\right)\right)\right)\right) \\
& -m_3 \left(l_2 \times \left(R_2^T R_3 \left(\dot{\Omega_3} \times \rho_3\right)\right)\right) \\
& -m_3 \left(\Omega_2 \times \left(l_2 \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_3 \left(\Omega_2 \times \left(l_2 \times \left(\Omega_2 \times l_2\right)\right)\right) \\
& -m_3 \left(\Omega_2 \times \left(l_2 \times \left(R_2^T R_3 \left(\Omega_3 \times \rho_3\right)\right)\right)\right) \\
& + m_4 \left(\left(\Omega_2 \times l_2\right) \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right) \\
& + m_4 \left(\left(\Omega_2 \times l_2\right) \times \left(R_2^T R_3 \left(\Omega_3 \times l_3\right)\right)\right) \\
& + m_4 \left(\left(\Omega_2 \times l_2\right) \times \left(R_2^T R_4 \left(\Omega_4 \times \rho_4\right)\right)\right) \\
& + m_4 \left(l_2 \times \left(\Omega_2 \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_4 \left(l_2 \times \left(R_2^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_4 \left(l_2 \times \left(R_2^T R_1 \left(\dot{\Omega_1} \times l_1\right)\right)\right) \\
& -m_4 \left(l_2 \times \left(\dot{\Omega_2} \times l_2\right)\right) \\
& + m_4 \left(l_2 \times \left(\Omega_2 \times \left(R_2^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_4 \left(l_2 \times \left(R_2^T R_3 \left(\Omega_3 \times \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_4 \left(l_2 \times \left(R_2^T R_3 \left(\dot{\Omega_3} \times l_3\right)\right)\right) \\
& + m_4 \left(l_2 \times \left(\Omega_2 \times \left(R_2^T R_4 \left(\Omega_4 \times \rho_4\right)\right)\right)\right) \\
& -m_4 \left(l_2 \times \left(R_2^T R_4 \left(\Omega_4 \times \left(\Omega_4 \times \rho_4\right)\right)\right)\right) \\
& -m_4 \left(l_2 \times \left(R_2^T R_4 \left(\dot{\Omega_4} \times \rho_4\right)\right)\right) \\
& -m_4 \left(\Omega_2 \times \left(l_2 \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_4 \left(\Omega_2 \times \left(l_2 \times \left(\Omega_2 \times l_2\right)\right)\right) \\
& -m_4 \left(\Omega_2 \times \left(l_2 \times \left(R_2^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_4 \left(\Omega_2 \times \left(l_2 \times \left(R_2^T R_4 \left(\Omega_4 \times \rho_4\right)\right)\right)\right) \\
& + m_5 \left(\left(\Omega_2 \times l_2\right) \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right) \\
& + m_5 \left(\left(\Omega_2 \times l_2\right) \times \left(R_2^T R_3 \left(\Omega_3 \times l_3\right)\right)\right) \\
& + m_5 \left(\left(\Omega_2 \times l_2\right) \times \left(R_2^T R_4 \left(\Omega_4 \times l_4\right)\right)\right) \\
& + m_5 \left(\left(\Omega_2 \times l_2\right) \times \left(R_2^T R_5 \left(\Omega_5 \times \rho_5\right)\right)\right) \\
& + m_5 \left(l_2 \times \left(\Omega_2 \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_5 \left(l_2 \times \left(R_2^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_5 \left(l_2 \times \left(R_2^T R_1 \left(\dot{\Omega_1} \times l_1\right)\right)\right) \\
& -m_5 \left(l_2 \times \left(\dot{\Omega_2} \times l_2\right)\right) \\
& + m_5 \left(l_2 \times \left(\Omega_2 \times \left(R_2^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_5 \left(l_2 \times \left(R_2^T R_3 \left(\Omega_3 \times \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_5 \left(l_2 \times \left(R_2^T R_3 \left(\dot{\Omega_3} \times l_3\right)\right)\right) \\
& + m_5 \left(l_2 \times \left(\Omega_2 \times \left(R_2^T R_4 \left(\Omega_4 \times l_4\right)\right)\right)\right) \\
& -m_5 \left(l_2 \times \left(R_2^T R_4 \left(\Omega_4 \times \left(\Omega_4 \times l_4\right)\right)\right)\right) \\
& -m_5 \left(l_2 \times \left(R_2^T R_4 \left(\dot{\Omega_4} \times l_4\right)\right)\right) \\
& + m_5 \left(l_2 \times \left(\Omega_2 \times \left(R_2^T R_5 \left(\Omega_5 \times \rho_5\right)\right)\right)\right) \\
& -m_5 \left(l_2 \times \left(R_2^T R_5 \left(\Omega_5 \times \left(\Omega_5 \times \rho_5\right)\right)\right)\right) \\
& -m_5 \left(l_2 \times \left(R_2^T R_5 \left(\dot{\Omega_5} \times \rho_5\right)\right)\right) \\
& -m_5 \left(\Omega_2 \times \left(l_2 \times \left(R_2^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_5 \left(\Omega_2 \times \left(l_2 \times \left(\Omega_2 \times l_2\right)\right)\right) \\
& -m_5 \left(\Omega_2 \times \left(l_2 \times \left(R_2^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_5 \left(\Omega_2 \times \left(l_2 \times \left(R_2^T R_4 \left(\Omega_4 \times l_4\right)\right)\right)\right) \\
& -m_5 \left(\Omega_2 \times \left(l_2 \times \left(R_2^T R_5 \left(\Omega_5 \times \rho_5\right)\right)\right)\right) \\
& -m_2 g \left(\rho_2 \times \left(R_2^T e_3\right)\right) \\
& -m_3 g \left(l_2 \times \left(R_2^T e_3\right)\right) \\
& -m_4 g \left(l_2 \times \left(R_2^T e_3\right)\right) \\
& -m_5 g \left(l_2 \times \left(R_2^T e_3\right)\right) \\
& + R_2^T \tau_2 \\
& -R_2^T \tau_3
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_3\)</span>, the equation of <span class="math notranslate nohighlight">\(R_3\)</span> (58 terms)</summary>

$$
\begin{aligned}
0 = {}& -J_3 \dot{\Omega_3} \\
& -\Omega_3 \times \left(J_3 \Omega_3\right) \\
& + m_3 \left(\left(\Omega_3 \times \rho_3\right) \times \left(R_3^T R_1 \left(\Omega_1 \times l_1\right)\right)\right) \\
& + m_3 \left(\left(\Omega_3 \times \rho_3\right) \times \left(R_3^T R_2 \left(\Omega_2 \times l_2\right)\right)\right) \\
& + m_3 \left(\rho_3 \times \left(\Omega_3 \times \left(R_3^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_3 \left(\rho_3 \times \left(R_3^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_3 \left(\rho_3 \times \left(R_3^T R_1 \left(\dot{\Omega_1} \times l_1\right)\right)\right) \\
& + m_3 \left(\rho_3 \times \left(\Omega_3 \times \left(R_3^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_3 \left(\rho_3 \times \left(R_3^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_3 \left(\rho_3 \times \left(R_3^T R_2 \left(\dot{\Omega_2} \times l_2\right)\right)\right) \\
& -m_3 \left(\rho_3 \times \left(\dot{\Omega_3} \times \rho_3\right)\right) \\
& -m_3 \left(\Omega_3 \times \left(\rho_3 \times \left(R_3^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_3 \left(\Omega_3 \times \left(\rho_3 \times \left(R_3^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_3 \left(\Omega_3 \times \left(\rho_3 \times \left(\Omega_3 \times \rho_3\right)\right)\right) \\
& + m_4 \left(\left(\Omega_3 \times l_3\right) \times \left(R_3^T R_1 \left(\Omega_1 \times l_1\right)\right)\right) \\
& + m_4 \left(\left(\Omega_3 \times l_3\right) \times \left(R_3^T R_2 \left(\Omega_2 \times l_2\right)\right)\right) \\
& + m_4 \left(\left(\Omega_3 \times l_3\right) \times \left(R_3^T R_4 \left(\Omega_4 \times \rho_4\right)\right)\right) \\
& + m_4 \left(l_3 \times \left(\Omega_3 \times \left(R_3^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_4 \left(l_3 \times \left(R_3^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_4 \left(l_3 \times \left(R_3^T R_1 \left(\dot{\Omega_1} \times l_1\right)\right)\right) \\
& + m_4 \left(l_3 \times \left(\Omega_3 \times \left(R_3^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_4 \left(l_3 \times \left(R_3^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_4 \left(l_3 \times \left(R_3^T R_2 \left(\dot{\Omega_2} \times l_2\right)\right)\right) \\
& -m_4 \left(l_3 \times \left(\dot{\Omega_3} \times l_3\right)\right) \\
& + m_4 \left(l_3 \times \left(\Omega_3 \times \left(R_3^T R_4 \left(\Omega_4 \times \rho_4\right)\right)\right)\right) \\
& -m_4 \left(l_3 \times \left(R_3^T R_4 \left(\Omega_4 \times \left(\Omega_4 \times \rho_4\right)\right)\right)\right) \\
& -m_4 \left(l_3 \times \left(R_3^T R_4 \left(\dot{\Omega_4} \times \rho_4\right)\right)\right) \\
& -m_4 \left(\Omega_3 \times \left(l_3 \times \left(R_3^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_4 \left(\Omega_3 \times \left(l_3 \times \left(R_3^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_4 \left(\Omega_3 \times \left(l_3 \times \left(\Omega_3 \times l_3\right)\right)\right) \\
& -m_4 \left(\Omega_3 \times \left(l_3 \times \left(R_3^T R_4 \left(\Omega_4 \times \rho_4\right)\right)\right)\right) \\
& + m_5 \left(\left(\Omega_3 \times l_3\right) \times \left(R_3^T R_1 \left(\Omega_1 \times l_1\right)\right)\right) \\
& + m_5 \left(\left(\Omega_3 \times l_3\right) \times \left(R_3^T R_2 \left(\Omega_2 \times l_2\right)\right)\right) \\
& + m_5 \left(\left(\Omega_3 \times l_3\right) \times \left(R_3^T R_4 \left(\Omega_4 \times l_4\right)\right)\right) \\
& + m_5 \left(\left(\Omega_3 \times l_3\right) \times \left(R_3^T R_5 \left(\Omega_5 \times \rho_5\right)\right)\right) \\
& + m_5 \left(l_3 \times \left(\Omega_3 \times \left(R_3^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_5 \left(l_3 \times \left(R_3^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_5 \left(l_3 \times \left(R_3^T R_1 \left(\dot{\Omega_1} \times l_1\right)\right)\right) \\
& + m_5 \left(l_3 \times \left(\Omega_3 \times \left(R_3^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_5 \left(l_3 \times \left(R_3^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_5 \left(l_3 \times \left(R_3^T R_2 \left(\dot{\Omega_2} \times l_2\right)\right)\right) \\
& -m_5 \left(l_3 \times \left(\dot{\Omega_3} \times l_3\right)\right) \\
& + m_5 \left(l_3 \times \left(\Omega_3 \times \left(R_3^T R_4 \left(\Omega_4 \times l_4\right)\right)\right)\right) \\
& -m_5 \left(l_3 \times \left(R_3^T R_4 \left(\Omega_4 \times \left(\Omega_4 \times l_4\right)\right)\right)\right) \\
& -m_5 \left(l_3 \times \left(R_3^T R_4 \left(\dot{\Omega_4} \times l_4\right)\right)\right) \\
& + m_5 \left(l_3 \times \left(\Omega_3 \times \left(R_3^T R_5 \left(\Omega_5 \times \rho_5\right)\right)\right)\right) \\
& -m_5 \left(l_3 \times \left(R_3^T R_5 \left(\Omega_5 \times \left(\Omega_5 \times \rho_5\right)\right)\right)\right) \\
& -m_5 \left(l_3 \times \left(R_3^T R_5 \left(\dot{\Omega_5} \times \rho_5\right)\right)\right) \\
& -m_5 \left(\Omega_3 \times \left(l_3 \times \left(R_3^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_5 \left(\Omega_3 \times \left(l_3 \times \left(R_3^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_5 \left(\Omega_3 \times \left(l_3 \times \left(\Omega_3 \times l_3\right)\right)\right) \\
& -m_5 \left(\Omega_3 \times \left(l_3 \times \left(R_3^T R_4 \left(\Omega_4 \times l_4\right)\right)\right)\right) \\
& -m_5 \left(\Omega_3 \times \left(l_3 \times \left(R_3^T R_5 \left(\Omega_5 \times \rho_5\right)\right)\right)\right) \\
& -m_3 g \left(\rho_3 \times \left(R_3^T e_3\right)\right) \\
& -m_4 g \left(l_3 \times \left(R_3^T e_3\right)\right) \\
& -m_5 g \left(l_3 \times \left(R_3^T e_3\right)\right) \\
& + R_3^T \tau_3 \\
& -R_3^T \tau_4
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_4\)</span>, the equation of <span class="math notranslate nohighlight">\(R_4\)</span> (45 terms)</summary>

$$
\begin{aligned}
0 = {}& -J_4 \dot{\Omega_4} \\
& -\Omega_4 \times \left(J_4 \Omega_4\right) \\
& + m_4 \left(\left(\Omega_4 \times \rho_4\right) \times \left(R_4^T R_1 \left(\Omega_1 \times l_1\right)\right)\right) \\
& + m_4 \left(\left(\Omega_4 \times \rho_4\right) \times \left(R_4^T R_2 \left(\Omega_2 \times l_2\right)\right)\right) \\
& + m_4 \left(\left(\Omega_4 \times \rho_4\right) \times \left(R_4^T R_3 \left(\Omega_3 \times l_3\right)\right)\right) \\
& + m_4 \left(\rho_4 \times \left(\Omega_4 \times \left(R_4^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_4 \left(\rho_4 \times \left(R_4^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_4 \left(\rho_4 \times \left(R_4^T R_1 \left(\dot{\Omega_1} \times l_1\right)\right)\right) \\
& + m_4 \left(\rho_4 \times \left(\Omega_4 \times \left(R_4^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_4 \left(\rho_4 \times \left(R_4^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_4 \left(\rho_4 \times \left(R_4^T R_2 \left(\dot{\Omega_2} \times l_2\right)\right)\right) \\
& + m_4 \left(\rho_4 \times \left(\Omega_4 \times \left(R_4^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_4 \left(\rho_4 \times \left(R_4^T R_3 \left(\Omega_3 \times \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_4 \left(\rho_4 \times \left(R_4^T R_3 \left(\dot{\Omega_3} \times l_3\right)\right)\right) \\
& -m_4 \left(\rho_4 \times \left(\dot{\Omega_4} \times \rho_4\right)\right) \\
& -m_4 \left(\Omega_4 \times \left(\rho_4 \times \left(R_4^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_4 \left(\Omega_4 \times \left(\rho_4 \times \left(R_4^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_4 \left(\Omega_4 \times \left(\rho_4 \times \left(R_4^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_4 \left(\Omega_4 \times \left(\rho_4 \times \left(\Omega_4 \times \rho_4\right)\right)\right) \\
& + m_5 \left(\left(\Omega_4 \times l_4\right) \times \left(R_4^T R_1 \left(\Omega_1 \times l_1\right)\right)\right) \\
& + m_5 \left(\left(\Omega_4 \times l_4\right) \times \left(R_4^T R_2 \left(\Omega_2 \times l_2\right)\right)\right) \\
& + m_5 \left(\left(\Omega_4 \times l_4\right) \times \left(R_4^T R_3 \left(\Omega_3 \times l_3\right)\right)\right) \\
& + m_5 \left(\left(\Omega_4 \times l_4\right) \times \left(R_4^T R_5 \left(\Omega_5 \times \rho_5\right)\right)\right) \\
& + m_5 \left(l_4 \times \left(\Omega_4 \times \left(R_4^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_5 \left(l_4 \times \left(R_4^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_5 \left(l_4 \times \left(R_4^T R_1 \left(\dot{\Omega_1} \times l_1\right)\right)\right) \\
& + m_5 \left(l_4 \times \left(\Omega_4 \times \left(R_4^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_5 \left(l_4 \times \left(R_4^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_5 \left(l_4 \times \left(R_4^T R_2 \left(\dot{\Omega_2} \times l_2\right)\right)\right) \\
& + m_5 \left(l_4 \times \left(\Omega_4 \times \left(R_4^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_5 \left(l_4 \times \left(R_4^T R_3 \left(\Omega_3 \times \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_5 \left(l_4 \times \left(R_4^T R_3 \left(\dot{\Omega_3} \times l_3\right)\right)\right) \\
& -m_5 \left(l_4 \times \left(\dot{\Omega_4} \times l_4\right)\right) \\
& + m_5 \left(l_4 \times \left(\Omega_4 \times \left(R_4^T R_5 \left(\Omega_5 \times \rho_5\right)\right)\right)\right) \\
& -m_5 \left(l_4 \times \left(R_4^T R_5 \left(\Omega_5 \times \left(\Omega_5 \times \rho_5\right)\right)\right)\right) \\
& -m_5 \left(l_4 \times \left(R_4^T R_5 \left(\dot{\Omega_5} \times \rho_5\right)\right)\right) \\
& -m_5 \left(\Omega_4 \times \left(l_4 \times \left(R_4^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_5 \left(\Omega_4 \times \left(l_4 \times \left(R_4^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_5 \left(\Omega_4 \times \left(l_4 \times \left(R_4^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_5 \left(\Omega_4 \times \left(l_4 \times \left(\Omega_4 \times l_4\right)\right)\right) \\
& -m_5 \left(\Omega_4 \times \left(l_4 \times \left(R_4^T R_5 \left(\Omega_5 \times \rho_5\right)\right)\right)\right) \\
& -m_4 g \left(\rho_4 \times \left(R_4^T e_3\right)\right) \\
& -m_5 g \left(l_4 \times \left(R_4^T e_3\right)\right) \\
& + R_4^T \tau_4 \\
& -R_4^T \tau_5
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_5\)</span>, the equation of <span class="math notranslate nohighlight">\(R_5\)</span> (26 terms)</summary>

$$
\begin{aligned}
0 = {}& -J_5 \dot{\Omega_5} \\
& -\Omega_5 \times \left(J_5 \Omega_5\right) \\
& + m_5 \left(\left(\Omega_5 \times \rho_5\right) \times \left(R_5^T R_1 \left(\Omega_1 \times l_1\right)\right)\right) \\
& + m_5 \left(\left(\Omega_5 \times \rho_5\right) \times \left(R_5^T R_2 \left(\Omega_2 \times l_2\right)\right)\right) \\
& + m_5 \left(\left(\Omega_5 \times \rho_5\right) \times \left(R_5^T R_3 \left(\Omega_3 \times l_3\right)\right)\right) \\
& + m_5 \left(\left(\Omega_5 \times \rho_5\right) \times \left(R_5^T R_4 \left(\Omega_4 \times l_4\right)\right)\right) \\
& + m_5 \left(\rho_5 \times \left(\Omega_5 \times \left(R_5^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_5 \left(\rho_5 \times \left(R_5^T R_1 \left(\Omega_1 \times \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_5 \left(\rho_5 \times \left(R_5^T R_1 \left(\dot{\Omega_1} \times l_1\right)\right)\right) \\
& + m_5 \left(\rho_5 \times \left(\Omega_5 \times \left(R_5^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_5 \left(\rho_5 \times \left(R_5^T R_2 \left(\Omega_2 \times \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_5 \left(\rho_5 \times \left(R_5^T R_2 \left(\dot{\Omega_2} \times l_2\right)\right)\right) \\
& + m_5 \left(\rho_5 \times \left(\Omega_5 \times \left(R_5^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_5 \left(\rho_5 \times \left(R_5^T R_3 \left(\Omega_3 \times \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_5 \left(\rho_5 \times \left(R_5^T R_3 \left(\dot{\Omega_3} \times l_3\right)\right)\right) \\
& + m_5 \left(\rho_5 \times \left(\Omega_5 \times \left(R_5^T R_4 \left(\Omega_4 \times l_4\right)\right)\right)\right) \\
& -m_5 \left(\rho_5 \times \left(R_5^T R_4 \left(\Omega_4 \times \left(\Omega_4 \times l_4\right)\right)\right)\right) \\
& -m_5 \left(\rho_5 \times \left(R_5^T R_4 \left(\dot{\Omega_4} \times l_4\right)\right)\right) \\
& -m_5 \left(\rho_5 \times \left(\dot{\Omega_5} \times \rho_5\right)\right) \\
& -m_5 \left(\Omega_5 \times \left(\rho_5 \times \left(R_5^T R_1 \left(\Omega_1 \times l_1\right)\right)\right)\right) \\
& -m_5 \left(\Omega_5 \times \left(\rho_5 \times \left(R_5^T R_2 \left(\Omega_2 \times l_2\right)\right)\right)\right) \\
& -m_5 \left(\Omega_5 \times \left(\rho_5 \times \left(R_5^T R_3 \left(\Omega_3 \times l_3\right)\right)\right)\right) \\
& -m_5 \left(\Omega_5 \times \left(\rho_5 \times \left(R_5^T R_4 \left(\Omega_4 \times l_4\right)\right)\right)\right) \\
& -m_5 \left(\Omega_5 \times \left(\rho_5 \times \left(\Omega_5 \times \rho_5\right)\right)\right) \\
& -m_5 g \left(\rho_5 \times \left(R_5^T e_3\right)\right) \\
& + R_5^T \tau_5
\end{aligned}
$$

</details>
<!-- /model:eom -->

## Code

The model for every $n$, as `geomech.indexed.systems.rigid_chain()` writes it, and the
library system for $n = 5$. The generic model takes a body moment $M_k$ on each link;
the library's joint torques give $M_k = R_k^T(\tau_k - \tau_{k+1})$ ($\tau_{n+1} = 0$),
and the Lagrangians are the same.

```python
from geomech import SO3, Dot, Matrix, Scalar, Vector
from geomech.core.base.constants import HALF
from geomech.indexed import Family, IndexedVariables, Sum, compute_indexed_eom, indices, n
from geomech.systems import get

i, j = indices("i j")
g = Scalar("g", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
R = Family(SO3, "R")  # R_1, ..., R_n
J = Family(Matrix, "J", attr=["Constant", "SymmetricMatrix"])
m = Family(Scalar, "m", ["Constant"])
rho, l = Family(Vector, "rho", ["Constant"]), Family(Vector, "l", ["Constant"])
M = Family(Vector, "M")


def x(k):  # centre of mass of link k: R_1 l_1 + ... + R_{k-1} l_{k-1} + R_k rho_k
    return Sum(j, 1, k - 1, R[j] * l[j]) + R[k] * rho[k]


def omega(k):
    return R[k].get_tangent_vector()


T = Sum(i, 1, n, HALF * Dot(omega(i), J[i] * omega(i)))
T = T + Sum(i, 1, n, HALF * m[i] * x(i).t_diff().dot(x(i).t_diff()))
V = Sum(i, 1, n, m[i] * g * x(i).dot(e3))
L = T - V
dW = Sum(i, 1, n, Dot(R[i].get_variation_vector(), M[i]))
eom = compute_indexed_eom(L, dW, IndexedVariables(families=[R]))  # every n

system = get("n_link_rigid_chain", n_links=5)  # n from 1 to 20; joint torques tau_k
eom_5 = system.derive()
```

The same model at $n = 5$ is `geomech.systems.get("n_link_rigid_chain", n_links=5)`, with
default parameters, a default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-7}$, with $H(0) = -40.4$ |
| Equations for every $n$ above, at $n = 5$: `check_eom` against the Euler–Lagrange evaluation of the generic $L$, $\delta W$ at $n = 5$ | below $10^{-6}$ |
| Generic $L$ at $n = 5$ against the library system's $L$, at random states | relative difference below $10^{-14}$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
