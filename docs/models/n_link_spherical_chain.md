# Chain of spherical pendulums

A chain of $n$ point masses on massless rods, each rod hinged to the previous mass by a
spherical joint and the first to a fixed pivot. The configuration is one unit vector per
link, $(q_1, \dots, q_n) \in (S^2)^n$, so the equations hold for any $n$ without
angle coordinates. The equations below are derived once for every $n$ with
`geomech.indexed`, then for $n = 5$ with the library system, which the animation uses.
Lee et al. §5.3 (general form eq. 5.14); for $n = 2$ the model is the
[double spherical pendulum](double_spherical_pendulum.md) of §5.6.5. Library system
`n_link_spherical_chain`.

## Animation

Ten seconds, unforced, five links, from the library's initial state: the links tilted
to the same side by increasing angles, from 40° (link 1) to 107° (link 5) off hanging,
slightly out of the $e_1$–$e_3$ plane. The chain swings in three dimensions; link 1
stays within 66° of hanging and link 5 within 114°. Drag to orbit; the full report (equations, standard form, parameters)
opens in a <a href="../model-reports/n_link_spherical_chain.html">separate page</a>.

<iframe class="model-report" src="../model-reports/n_link_spherical_chain.html?view" title="Chain of spherical pendulums, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 400 300" role="img" aria-labelledby="nsc-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="nsc-title">Chain of three spherical pendulums: link directions q1, q2, q3, lengths l1, l2, l3, masses m1, m2, m3, angular velocities ω1, ω2, ω3, input moments M1, M2, M3, gravity</title>
  <defs>
    <marker id="nsc-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- support -->
    <line x1="80" y1="45" x2="140" y2="45" stroke-width="2"/>
    <path d="M82,45 l-8,-10 M94,45 l-8,-10 M106,45 l-8,-10 M118,45 l-8,-10 M130,45 l-8,-10 M142,45 l-8,-10"/>
    <!-- inertial axes -->
    <line x1="40" y1="270" x2="40" y2="230" marker-end="url(#nsc-arr)"/>
    <line x1="40" y1="270" x2="80" y2="270" marker-end="url(#nsc-arr)"/>
    <!-- links -->
    <line x1="110.0" y1="45.0" x2="150.4" y2="155.9" stroke-width="2.2"/>
    <line x1="150.4" y1="155.9" x2="229.3" y2="229.5" stroke-width="2.2"/>
    <line x1="229.3" y1="229.5" x2="322.5" y2="259.8" stroke-width="2.2"/>
  </g>
  <!-- q_i: unit vectors along the links; M_i: input moments on the links -->
  <line x1="110.0" y1="45.0" x2="124.4" y2="84.5" stroke="#3f8fa0" stroke-width="3" marker-end="url(#nsc-arr)"/>
  <path d="M 88.1,43.1 A 22 22 0 0 0 97.4,63.0" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#nsc-arr)"/>
  <line x1="150.4" y1="155.9" x2="181.1" y2="184.5" stroke="#3f8fa0" stroke-width="3" marker-end="url(#nsc-arr)"/>
  <path d="M 130.0,164.1 A 22 22 0 0 0 147.3,177.7" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#nsc-arr)"/>
  <line x1="229.3" y1="229.5" x2="269.3" y2="242.5" stroke="#3f8fa0" stroke-width="3" marker-end="url(#nsc-arr)"/>
  <path d="M 214.3,245.6 A 22 22 0 0 0 235.8,250.6" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#nsc-arr)"/>
  <!-- gravity -->
  <line x1="360" y1="60" x2="360" y2="120" stroke="#e08a2c" stroke-width="2" marker-end="url(#nsc-arr)"/>
  <!-- omega_i: perpendicular to q_i (out of the page here); masses -->
  <circle cx="110" cy="45" r="4" fill="currentColor"/>
  <circle cx="163.5" cy="115.9" r="7" fill="none" stroke="#3f8fa0" stroke-width="1.8"/>
  <circle cx="163.5" cy="115.9" r="2" fill="#3f8fa0"/>
  <circle cx="150.4" cy="155.9" r="10" fill="currentColor"/>
  <circle cx="225.0" cy="189.9" r="7" fill="none" stroke="#3f8fa0" stroke-width="1.8"/>
  <circle cx="225.0" cy="189.9" r="2" fill="#3f8fa0"/>
  <circle cx="229.3" cy="229.5" r="10" fill="currentColor"/>
  <circle cx="304.5" cy="226.6" r="7" fill="none" stroke="#3f8fa0" stroke-width="1.8"/>
  <circle cx="304.5" cy="226.6" r="2" fill="#3f8fa0"/>
  <circle cx="322.5" cy="259.8" r="10" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="138.6" y="68.2" text-anchor="middle" fill="#3f8fa0">q<tspan font-size="12" dy="4">1</tspan></text>
    <text x="176.5" y="121.9" fill="#3f8fa0">ω<tspan font-size="12" dy="4">1</tspan></text>
    <text x="113.3" y="112.6" text-anchor="middle">l<tspan font-size="12" dy="4">1</tspan></text>
    <text x="73.7" y="67.9" text-anchor="middle" fill="#e08a2c">M<tspan font-size="12" dy="4">1</tspan></text>
    <text x="120.4" y="149.9" text-anchor="middle">m<tspan font-size="12" dy="4">1</tspan></text>
    <text x="183.7" y="164.3" text-anchor="middle" fill="#3f8fa0">q<tspan font-size="12" dy="4">2</tspan></text>
    <text x="238.0" y="195.9" fill="#3f8fa0">ω<tspan font-size="12" dy="4">2</tspan></text>
    <text x="177.6" y="211.9" text-anchor="middle">l<tspan font-size="12" dy="4">2</tspan></text>
    <text x="125.7" y="193.4" text-anchor="middle" fill="#e08a2c">M<tspan font-size="12" dy="4">2</tspan></text>
    <text x="199.3" y="231.5" text-anchor="middle">m<tspan font-size="12" dy="4">2</tspan></text>
    <text x="260.6" y="223.6" text-anchor="middle" fill="#3f8fa0">q<tspan font-size="12" dy="4">3</tspan></text>
    <text x="317.5" y="232.6" fill="#3f8fa0">ω<tspan font-size="12" dy="4">3</tspan></text>
    <text x="270.4" y="267.8" text-anchor="middle">l<tspan font-size="12" dy="4">3</tspan></text>
    <text x="220.3" y="274.5" text-anchor="middle" fill="#e08a2c">M<tspan font-size="12" dy="4">3</tspan></text>
    <text x="338.5" y="287.8" text-anchor="middle">m<tspan font-size="12" dy="4">3</tspan></text>
    <text x="368" y="98" fill="#e08a2c">g</text>
    <text x="47" y="242">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="85" y="275">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>Drawn with three links. Mass <span class="math notranslate nohighlight">\(i\)</span> is at <span class="math notranslate nohighlight">\(x_i = l_1 q_1 + \dots + l_i q_i\)</span> from
the pivot. Each angular velocity <span class="math notranslate nohighlight">\(\omega_i\)</span> is perpendicular to its <span class="math notranslate nohighlight">\(q_i\)</span>, out of the page
in this picture, and <span class="math notranslate nohighlight">\(\dot q_i = \omega_i \times q_i\)</span>. Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span>; <span class="math notranslate nohighlight">\(M_i\)</span> is a moment
on link <span class="math notranslate nohighlight">\(i\)</span>, with virtual work <span class="math notranslate nohighlight">\(\sum_i \xi_i \cdot M_i\)</span>.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $q_i$ | direction of link $i$, from mass $i-1$ (or the pivot) to mass $i$ | $S^2$ each |
| $\omega_i$ | angular velocity of link $i$, $\dot q_i = \omega_i \times q_i$, $q_i \cdot \omega_i = 0$ | $T_{q_i} S^2$ |
| $\xi_i$ | variation of $q_i$: $\delta q_i = \xi_i \times q_i$ | $T_{q_i} S^2$ |
| $M_i$ | input moment on link $i$ | $\mathbb{R}^3$ each |
| $m_i, l_i$ | mass at the end of link $i$, length of link $i$ | constants |
| $g, e_3$ | gravity, vertical unit vector | constants |

## Lagrangian and equations for every $n$

<!-- model:indexed -->
Lagrangian for every $n$, as `geomech.indexed.systems.spherical_chain()` holds it (the velocities expanded, so the inner sums are the inertia couplings):

$$
\begin{aligned}
L = {}& \sum_{j=1}^{n} \sum_{k=1}^{n} \sum_{r=\max(j, k)}^{n} \frac{1}{2} m_{r} \left(l_{j} \left(\omega_{j} \times q_{j}\right) \cdot l_{k} \left(\omega_{k} \times q_{k}\right)\right) \\
& -\sum_{j=1}^{n} \sum_{k=j}^{n} m_{k} g \left(l_{j} q_{j} \cdot e_3\right)
\end{aligned}
$$

Virtual work:

$$
\begin{aligned}
\delta W = {}& \sum_{j=1}^{n} \xi_{j} \cdot M_{j}
\end{aligned}
$$

Equations from `compute_indexed_eom`, derived once for every $n$; $i$ is the body, $1 \le i \le n$:

Coefficient of $\xi_{i}$, the equation of $q_i$:

$$
\begin{aligned}
0 = {}& M_{i} \\
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
T = {}& \frac{1}{2} m_1 \left(l_1 \left(\omega_1 \times q_1\right) \cdot l_1 \left(\omega_1 \times q_1\right)\right) \\
& + \frac{1}{2} m_2 \left(\left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right) \cdot \left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right)\right) \\
& + \frac{1}{2} m_3 \left(\left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right) \cdot \left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right)\right) \\
& + \frac{1}{2} m_4 \left(\left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right)\right) \cdot \left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right)\right)\right) \\
& + \frac{1}{2} m_5 \left(\left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right) + l_5 \left(\omega_5 \times q_5\right)\right) \cdot \left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right) + l_5 \left(\omega_5 \times q_5\right)\right)\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& \left(m_1 g \left(l_1 q_1 \cdot e_3\right) + m_2 g \left(\left(l_1 q_1 + l_2 q_2\right) \cdot e_3\right) + m_3 g \left(\left(l_1 q_1 + l_2 q_2 + l_3 q_3\right) \cdot e_3\right) + m_4 g \left(\left(l_1 q_1 + l_2 q_2 + l_3 q_3 + l_4 q_4\right) \cdot e_3\right) + m_5 g \left(\left(l_1 q_1 + l_2 q_2 + l_3 q_3 + l_4 q_4 + l_5 q_5\right) \cdot e_3\right)\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m_1 \left(l_1 \left(\omega_1 \times q_1\right) \cdot l_1 \left(\omega_1 \times q_1\right)\right) \\
& + \frac{1}{2} m_2 \left(\left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right) \cdot \left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right)\right)\right) \\
& + \frac{1}{2} m_3 \left(\left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right) \cdot \left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right)\right)\right) \\
& + \frac{1}{2} m_4 \left(\left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right)\right) \cdot \left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right)\right)\right) \\
& + \frac{1}{2} m_5 \left(\left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right) + l_5 \left(\omega_5 \times q_5\right)\right) \cdot \left(l_1 \left(\omega_1 \times q_1\right) + l_2 \left(\omega_2 \times q_2\right) + l_3 \left(\omega_3 \times q_3\right) + l_4 \left(\omega_4 \times q_4\right) + l_5 \left(\omega_5 \times q_5\right)\right)\right) \\
& -\left(m_1 g \left(l_1 q_1 \cdot e_3\right) + m_2 g \left(\left(l_1 q_1 + l_2 q_2\right) \cdot e_3\right) + m_3 g \left(\left(l_1 q_1 + l_2 q_2 + l_3 q_3\right) \cdot e_3\right) + m_4 g \left(\left(l_1 q_1 + l_2 q_2 + l_3 q_3 + l_4 q_4\right) \cdot e_3\right) + m_5 g \left(\left(l_1 q_1 + l_2 q_2 + l_3 q_3 + l_4 q_4 + l_5 q_5\right) \cdot e_3\right)\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \xi_1 \cdot M_1 \\
& + \xi_2 \cdot M_2 \\
& + \xi_3 \cdot M_3 \\
& + \xi_4 \cdot M_4 \\
& + \xi_5 \cdot M_5
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion at $n = 5$

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); the equation of $q_i$ holds on the tangent plane at
$q_i$.

<!-- model:eom -->
<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_1\)</span>, the equation of <span class="math notranslate nohighlight">\(q_1\)</span> (61 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_1 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -m_2 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -m_2 l_1 l_2 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_2 l_1 l_2 \left(q_1 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_2 l_1 l_2 \left(q_1 \times q_2\right) \\
& + 2 m_2 l_1 l_2 \left(\omega_1 \times \left(q_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_2 l_1 l_2 \left(q_1 \times \left(\omega_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_3 l_1 l_1 \left(q_1 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& -m_3 l_1 l_2 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_3 l_1 l_2 \left(q_1 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_3 l_1 l_2 \left(q_1 \times q_2\right) \\
& -m_3 l_1 l_3 \left(\left(\omega_1 \times q_1\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_3 l_1 l_3 \left(q_1 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_3 l_1 l_3 \left(q_1 \times q_3\right) \\
& + 2 m_3 l_1 l_2 \left(\omega_1 \times \left(q_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_3 l_1 l_3 \left(\omega_1 \times \left(q_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_3 l_1 l_2 \left(q_1 \times \left(\omega_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_3 l_1 l_3 \left(q_1 \times \left(\omega_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
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
& + 2 m_4 l_1 l_2 \left(\omega_1 \times \left(q_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_4 l_1 l_3 \left(\omega_1 \times \left(q_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& + 2 m_4 l_1 l_4 \left(\omega_1 \times \left(q_1 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_4 l_1 l_2 \left(q_1 \times \left(\omega_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_4 l_1 l_3 \left(q_1 \times \left(\omega_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_4 l_1 l_4 \left(q_1 \times \left(\omega_1 \times \left(\omega_4 \times q_4\right)\right)\right) \\
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
& + 2 m_5 l_1 l_2 \left(\omega_1 \times \left(q_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_5 l_1 l_3 \left(\omega_1 \times \left(q_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& + 2 m_5 l_1 l_4 \left(\omega_1 \times \left(q_1 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& + 2 m_5 l_1 l_5 \left(\omega_1 \times \left(q_1 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -m_5 l_1 l_2 \left(q_1 \times \left(\omega_1 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_5 l_1 l_3 \left(q_1 \times \left(\omega_1 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_5 l_1 l_4 \left(q_1 \times \left(\omega_1 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_5 l_1 l_5 \left(q_1 \times \left(\omega_1 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -l_1 m_1 g \left(q_1 \times e_3\right) \\
& -l_1 m_2 g \left(q_1 \times e_3\right) \\
& -l_1 m_3 g \left(q_1 \times e_3\right) \\
& -l_1 m_4 g \left(q_1 \times e_3\right) \\
& -l_1 m_5 g \left(q_1 \times e_3\right) \\
& + M_1
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_2\)</span>, the equation of <span class="math notranslate nohighlight">\(q_2\)</span> (59 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_2 l_2 l_1 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_2 l_2 l_1 \left(q_2 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_2 l_2 l_1 \left(q_2 \times q_1\right) \\
& -m_2 l_2 l_2 \left(q_2 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + 2 m_2 l_2 l_1 \left(\omega_2 \times \left(q_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_2 l_2 l_1 \left(q_2 \times \left(\omega_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_3 l_2 l_1 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_3 l_2 l_1 \left(q_2 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_3 l_2 l_1 \left(q_2 \times q_1\right) \\
& -m_3 l_2 l_2 \left(q_2 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& -m_3 l_2 l_3 \left(\left(\omega_2 \times q_2\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_3 l_2 l_3 \left(q_2 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_3 l_2 l_3 \left(q_2 \times q_3\right) \\
& + 2 m_3 l_2 l_1 \left(\omega_2 \times \left(q_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_3 l_2 l_3 \left(\omega_2 \times \left(q_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_3 l_2 l_1 \left(q_2 \times \left(\omega_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_3 l_2 l_3 \left(q_2 \times \left(\omega_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
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
& + 2 m_4 l_2 l_1 \left(\omega_2 \times \left(q_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_4 l_2 l_3 \left(\omega_2 \times \left(q_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& + 2 m_4 l_2 l_4 \left(\omega_2 \times \left(q_2 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_4 l_2 l_1 \left(q_2 \times \left(\omega_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_4 l_2 l_3 \left(q_2 \times \left(\omega_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_4 l_2 l_4 \left(q_2 \times \left(\omega_2 \times \left(\omega_4 \times q_4\right)\right)\right) \\
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
& + 2 m_5 l_2 l_1 \left(\omega_2 \times \left(q_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_5 l_2 l_3 \left(\omega_2 \times \left(q_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& + 2 m_5 l_2 l_4 \left(\omega_2 \times \left(q_2 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& + 2 m_5 l_2 l_5 \left(\omega_2 \times \left(q_2 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -m_5 l_2 l_1 \left(q_2 \times \left(\omega_2 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_5 l_2 l_3 \left(q_2 \times \left(\omega_2 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_5 l_2 l_4 \left(q_2 \times \left(\omega_2 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_5 l_2 l_5 \left(q_2 \times \left(\omega_2 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -l_2 m_2 g \left(q_2 \times e_3\right) \\
& -l_2 m_3 g \left(q_2 \times e_3\right) \\
& -l_2 m_4 g \left(q_2 \times e_3\right) \\
& -l_2 m_5 g \left(q_2 \times e_3\right) \\
& + M_2
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_3\)</span>, the equation of <span class="math notranslate nohighlight">\(q_3\)</span> (52 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_3 l_3 l_1 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_3 l_3 l_1 \left(q_3 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_3 l_3 l_1 \left(q_3 \times q_1\right) \\
& -m_3 l_3 l_2 \left(\left(\omega_3 \times q_3\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_3 l_3 l_2 \left(q_3 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_3 l_3 l_2 \left(q_3 \times q_2\right) \\
& -m_3 l_3 l_3 \left(q_3 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + 2 m_3 l_3 l_1 \left(\omega_3 \times \left(q_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_3 l_3 l_2 \left(\omega_3 \times \left(q_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_3 l_3 l_1 \left(q_3 \times \left(\omega_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_3 l_3 l_2 \left(q_3 \times \left(\omega_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
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
& + 2 m_4 l_3 l_1 \left(\omega_3 \times \left(q_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_4 l_3 l_2 \left(\omega_3 \times \left(q_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_4 l_3 l_4 \left(\omega_3 \times \left(q_3 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_4 l_3 l_1 \left(q_3 \times \left(\omega_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_4 l_3 l_2 \left(q_3 \times \left(\omega_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_4 l_3 l_4 \left(q_3 \times \left(\omega_3 \times \left(\omega_4 \times q_4\right)\right)\right) \\
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
& + 2 m_5 l_3 l_1 \left(\omega_3 \times \left(q_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_5 l_3 l_2 \left(\omega_3 \times \left(q_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_5 l_3 l_4 \left(\omega_3 \times \left(q_3 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& + 2 m_5 l_3 l_5 \left(\omega_3 \times \left(q_3 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -m_5 l_3 l_1 \left(q_3 \times \left(\omega_3 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_5 l_3 l_2 \left(q_3 \times \left(\omega_3 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_5 l_3 l_4 \left(q_3 \times \left(\omega_3 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_5 l_3 l_5 \left(q_3 \times \left(\omega_3 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -l_3 m_3 g \left(q_3 \times e_3\right) \\
& -l_3 m_4 g \left(q_3 \times e_3\right) \\
& -l_3 m_5 g \left(q_3 \times e_3\right) \\
& + M_3
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_4\)</span>, the equation of <span class="math notranslate nohighlight">\(q_4\)</span> (40 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_4 l_4 l_1 \left(\left(\omega_4 \times q_4\right) \times \left(\omega_1 \times q_1\right)\right) \\
& -m_4 l_4 l_1 \left(q_4 \times \left(\dot{\omega_1} \times q_1\right)\right) \\
& + \left(\omega_1 \cdot \omega_1\right) m_4 l_4 l_1 \left(q_4 \times q_1\right) \\
& -m_4 l_4 l_2 \left(\left(\omega_4 \times q_4\right) \times \left(\omega_2 \times q_2\right)\right) \\
& -m_4 l_4 l_2 \left(q_4 \times \left(\dot{\omega_2} \times q_2\right)\right) \\
& + \left(\omega_2 \cdot \omega_2\right) m_4 l_4 l_2 \left(q_4 \times q_2\right) \\
& -m_4 l_4 l_3 \left(\left(\omega_4 \times q_4\right) \times \left(\omega_3 \times q_3\right)\right) \\
& -m_4 l_4 l_3 \left(q_4 \times \left(\dot{\omega_3} \times q_3\right)\right) \\
& + \left(\omega_3 \cdot \omega_3\right) m_4 l_4 l_3 \left(q_4 \times q_3\right) \\
& -m_4 l_4 l_4 \left(q_4 \times \left(\dot{\omega_4} \times q_4\right)\right) \\
& + 2 m_4 l_4 l_1 \left(\omega_4 \times \left(q_4 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_4 l_4 l_2 \left(\omega_4 \times \left(q_4 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_4 l_4 l_3 \left(\omega_4 \times \left(q_4 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_4 l_4 l_1 \left(q_4 \times \left(\omega_4 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_4 l_4 l_2 \left(q_4 \times \left(\omega_4 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_4 l_4 l_3 \left(q_4 \times \left(\omega_4 \times \left(\omega_3 \times q_3\right)\right)\right) \\
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
& + 2 m_5 l_4 l_1 \left(\omega_4 \times \left(q_4 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_5 l_4 l_2 \left(\omega_4 \times \left(q_4 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_5 l_4 l_3 \left(\omega_4 \times \left(q_4 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& + 2 m_5 l_4 l_5 \left(\omega_4 \times \left(q_4 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -m_5 l_4 l_1 \left(q_4 \times \left(\omega_4 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_5 l_4 l_2 \left(q_4 \times \left(\omega_4 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_5 l_4 l_3 \left(q_4 \times \left(\omega_4 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_5 l_4 l_5 \left(q_4 \times \left(\omega_4 \times \left(\omega_5 \times q_5\right)\right)\right) \\
& -l_4 m_4 g \left(q_4 \times e_3\right) \\
& -l_4 m_5 g \left(q_4 \times e_3\right) \\
& + M_4
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\xi_5\)</span>, the equation of <span class="math notranslate nohighlight">\(q_5\)</span> (23 terms)</summary>

$$
\begin{aligned}
0 = {}& -m_5 l_5 l_1 \left(\left(\omega_5 \times q_5\right) \times \left(\omega_1 \times q_1\right)\right) \\
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
& + 2 m_5 l_5 l_1 \left(\omega_5 \times \left(q_5 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& + 2 m_5 l_5 l_2 \left(\omega_5 \times \left(q_5 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& + 2 m_5 l_5 l_3 \left(\omega_5 \times \left(q_5 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& + 2 m_5 l_5 l_4 \left(\omega_5 \times \left(q_5 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -m_5 l_5 l_1 \left(q_5 \times \left(\omega_5 \times \left(\omega_1 \times q_1\right)\right)\right) \\
& -m_5 l_5 l_2 \left(q_5 \times \left(\omega_5 \times \left(\omega_2 \times q_2\right)\right)\right) \\
& -m_5 l_5 l_3 \left(q_5 \times \left(\omega_5 \times \left(\omega_3 \times q_3\right)\right)\right) \\
& -m_5 l_5 l_4 \left(q_5 \times \left(\omega_5 \times \left(\omega_4 \times q_4\right)\right)\right) \\
& -l_5 m_5 g \left(q_5 \times e_3\right) \\
& + M_5
\end{aligned}
$$

</details>
<!-- /model:eom -->

## Code

The model for every $n$, as `geomech.indexed.systems.spherical_chain()` writes it, and the
library system for $n = 5$:

```python
from geomech import S2, Dot, Scalar, Vector
from geomech.core.base.constants import HALF
from geomech.indexed import Family, IndexedVariables, Sum, compute_indexed_eom, indices, n
from geomech.systems import get

i, j = indices("i j")
g = Scalar("g", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
q = Family(S2, "q")  # q_1, ..., q_n
m, l = Family(Scalar, "m", ["Constant"]), Family(Scalar, "l", ["Constant"])
M = Family(Vector, "M")


def x(k):  # mass k at x_k = l_1 q_1 + ... + l_k q_k
    return Sum(j, 1, k, q[j] * l[j])


def v(k):
    return x(k).t_diff()


L = Sum(i, 1, n, HALF * m[i] * v(i).dot(v(i))) - Sum(i, 1, n, m[i] * g * x(i).dot(e3))
dW = Sum(i, 1, n, Dot(q[i].get_variation_vector(), M[i]))
eom = compute_indexed_eom(L, dW, IndexedVariables(families=[q]))  # every n

system = get("n_link_spherical_chain", n_links=5)  # n from 1 to 20
eom_5 = system.derive()
```

`eom.instantiate(5)` gives the equations at $n = 5$ in the form `compute_eom` returns.
The same model at $n = 5$ is `geomech.systems.get("n_link_spherical_chain", n_links=5)`,
with default parameters, a default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-5}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-7}$, with $H(0) = -52.8$ |
| Equations for every $n$ above, at $n = 5$: `check_eom` against the Euler–Lagrange evaluation of the generic $L$, $\delta W$ at $n = 5$ | below $10^{-6}$ |
| Generic $L$ at $n = 5$ against the library system's $L$, at random states | relative difference below $10^{-14}$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
