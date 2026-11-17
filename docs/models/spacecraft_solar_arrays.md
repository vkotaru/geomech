# Spacecraft with solar arrays

The spacecraft of [Spacecraft with reaction wheels on S¹](spacecraft_reaction_wheels_s1.md)
with two solar-array panels. Each panel is a rigid plate that a drive turns relative to the
bus about a bus-fixed axis $a_k$; its orientation is its normal $n_k \in S^1$. The
configuration space is $SO(3) \times (S^1)^3 \times (S^1)^2$. The spacecraft floats free
with zero linear momentum, so its centre of mass is a fixed point of the bus frame and a
panel's centre of mass moves only as the bus turns. Not taken from the literature; the
equations are verified numerically (`check_eom`) and against an independent model of the
bodies (`tests/test_systems_physics.py`). Library system `spacecraft_solar_arrays`.

## Animation

Ten seconds, unforced, from the library's initial state: the bus and wheels start as in the
reaction-wheel models (bus at $\Omega_R = (0.05, -0.1, 0.08)$ rad/s, wheels at 30, −20 and
10 rad/s), both panels face along $b_3$, and the drives turn at $\omega_{n_1} = 0.5$ and
$\omega_{n_2} = -0.4$ rad/s. With no drive torque the panels keep turning, but their rates
are not constant: through the coupling with the bus, $\omega_{n_1}$ moves between 0.39 and
0.52 rad/s and $\omega_{n_2}$ between −0.42 and −0.30 rad/s. The bus ends about 80° away from
its initial attitude; the wheel rates stay within 0.12 rad/s of their start. Drag to orbit; the full report (equations, standard form,
parameters) opens in a
<a href="../model-reports/spacecraft_solar_arrays.html">separate page</a>.

<iframe class="model-report" src="../model-reports/spacecraft_solar_arrays.html?view" title="Spacecraft with solar arrays, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 400 250" role="img" aria-labelledby="ssa-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="ssa-title">Spacecraft with solar arrays: bus with attitude R and angular velocity Omega_R, two panels on drive axes a1 and a2 at distances d1 and d2, panel normals n1 and n2, drive torques mu1 and mu2</title>
  <defs>
    <marker id="ssa-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- bus: front face, top face, right face -->
    <rect x="170" y="125" width="60" height="60" stroke-width="1.6"/>
    <path d="M170,125 L194,113 L254,113 L230,125 M254,113 L254,173 L230,185" stroke-width="1.6"/>
    <!-- drive axis b1 through the centre of the bus (dashed inside it) -->
    <line x1="170" y1="149" x2="242" y2="149" stroke-dasharray="3 4" opacity="0.6"/>
    <line x1="242" y1="149" x2="385" y2="149" marker-end="url(#ssa-arr)"/>
    <line x1="170" y1="149" x2="30" y2="149" marker-end="url(#ssa-arr)"/>
    <!-- panels, turned about the drive axis -->
    <path d="M239,133 L347,133 L357,165 L249,165 z" fill="currentColor" fill-opacity="0.08"/>
    <path d="M67,133 L158,133 L168,165 L77,165 z" fill="currentColor" fill-opacity="0.08"/>
    <!-- distances from the centre of the bus to the panel centres -->
    <path d="M117.5,212 L298,212 M117.5,206 L117.5,218 M212,206 L212,218 M298,206 L298,218" opacity="0.7"/>
  </g>
  <!-- panel normals -->
  <line x1="298" y1="149" x2="329" y2="113" stroke="#3f8fa0" stroke-width="2.2" marker-end="url(#ssa-arr)"/>
  <line x1="117.5" y1="149" x2="148.5" y2="113.5" stroke="#3f8fa0" stroke-width="2.2" marker-end="url(#ssa-arr)"/>
  <!-- body angular velocity -->
  <line x1="212" y1="119" x2="230" y2="55" stroke="#3f8fa0" stroke-width="2.2" marker-end="url(#ssa-arr)"/>
  <circle cx="212" cy="149" r="3" fill="currentColor"/>
  <circle cx="298" cy="149" r="2.5" fill="currentColor"/>
  <circle cx="117.5" cy="149" r="2.5" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="366" y="136">a<tspan font-size="12" dy="4">1</tspan></text>
    <text x="20" y="136">a<tspan font-size="12" dy="4">2</tspan></text>
    <text x="338" y="104" fill="#3f8fa0">n<tspan font-size="12" dy="4">1</tspan></text>
    <text x="128" y="100" fill="#3f8fa0">n<tspan font-size="12" dy="4">2</tspan></text>
    <text x="242" y="64" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">R</tspan></text>
    <text x="176" y="180" fill="#3f8fa0">R</text>
    <text x="330" y="190" fill="#e08a2c">μ<tspan font-size="12" dy="4">1</tspan></text>
    <text x="70" y="190" fill="#e08a2c">μ<tspan font-size="12" dy="4">2</tspan></text>
    <text x="250" y="234">d<tspan font-size="12" dy="4">1</tspan></text>
    <text x="158" y="234">d<tspan font-size="12" dy="4">2</tspan></text>
  </g>
</svg>
<figcaption>Panel <span class="math notranslate nohighlight">\(k\)</span> turns about its drive axis <span class="math notranslate nohighlight">\(a_k\)</span> (<span class="math notranslate nohighlight">\(a_1 = b_1\)</span>,
<span class="math notranslate nohighlight">\(a_2 = -b_1\)</span>), with its centre of mass on that axis at <span class="math notranslate nohighlight">\(d_k a_k\)</span>
from the centre of mass of the spacecraft; <span class="math notranslate nohighlight">\(n_k\)</span> is its normal, written in the
body frame, and the drive applies <span class="math notranslate nohighlight">\(\mu_k\)</span> to the panel and <span class="math notranslate nohighlight">\(-\mu_k\)</span> to
the bus. <span class="math notranslate nohighlight">\(R\)</span> and <span class="math notranslate nohighlight">\(\Omega_R\)</span> are the attitude and body angular velocity of the
bus. The three reaction wheels (<span class="math notranslate nohighlight">\(q_i\)</span>, <span class="math notranslate nohighlight">\(\tau_i\)</span>) are those of the
<a href="spacecraft_reaction_wheels_s1.html"><span class="math notranslate nohighlight">\(S^1\)</span> page</a> and are left out of this drawing. Not to
scale.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $R$ | attitude of the bus, body frame to inertial frame | $SO(3)$ |
| $\Omega_{R}$ | body angular velocity of the bus, $\dot R = R \widehat{\Omega_{R}}$ | $\mathbb{R}^3$ |
| $\eta_{R}$ | variation of $R$: $\delta R = R \widehat{\eta_{R}}$ | $\mathbb{R}^3$ |
| $q_i$ | wheel-fixed direction of wheel $i$, orthogonal to $b_i$, in the body frame, $i = 1, 2, 3$ | $S^1$ |
| $\omega_{q_i}$ | spin rate of wheel $i$ relative to the bus, $\dot q_i = \omega_{q_i}\, b_i \times q_i$ | $\mathbb{R}$ |
| $\gamma_{q_i}$ | variation of $q_i$: $\delta q_i = \gamma_{q_i}\, b_i \times q_i$ | $\mathbb{R}$ |
| $\tau_i$ | motor torque on wheel $i$ about $b_i$ (the bus gets $-\tau_i b_i$) | $\mathbb{R}$ |
| $n_k$ | normal of panel $k$ in the body frame, orthogonal to $a_k$, $k = 1, 2$ | $S^1$ |
| $\omega_{n_k}$ | drive rate of panel $k$ relative to the bus, $\dot n_k = \omega_{n_k}\, a_k \times n_k$ | $\mathbb{R}$ |
| $\gamma_{n_k}$ | variation of $n_k$: $\delta n_k = \gamma_{n_k}\, a_k \times n_k$ | $\mathbb{R}$ |
| $\mu_k$ | drive torque on panel $k$ about $a_k$ (the bus gets $-\mu_k a_k$) | $\mathbb{R}$ |
| $J_b, I_i, b_i$ | bus inertia (wheels locked, minus their axial inertias), wheel axial inertias, wheel axes | constants |
| $a_k, d_k$ | drive axis of panel $k$ (unit, body frame), distance of its centre of mass along $a_k$ | constants |
| $m_{p1}, m_{p2}$ | panel masses | constants |
| $I_{a,k}, I_{n,k}, I_{w,k}$ | principal moments of inertia of panel $k$ about its centre of mass: about $a_k$, $n_k$ and $a_k \times n_k$ | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} \left(\Omega_{R} \cdot J_b \Omega_{R}\right) \\
& + \frac{1}{2} I_1 \left(\omega_{q_1} + b_1 \cdot \Omega_{R}\right) \left(\omega_{q_1} + b_1 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_2 \left(\omega_{q_2} + b_2 \cdot \Omega_{R}\right) \left(\omega_{q_2} + b_2 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_3 \left(\omega_{q_3} + b_3 \cdot \Omega_{R}\right) \left(\omega_{q_3} + b_3 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} m_{p1} d_1 d_1 \left(\Omega_{R} \times a_1 \cdot \Omega_{R} \times a_1\right) \\
& + \frac{1}{2} I_{a,1} \left(\omega_{n_1} + a_1 \cdot \Omega_{R}\right) \left(\omega_{n_1} + a_1 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_{n,1} \left(n_1 \cdot \Omega_{R}\right) \left(n_1 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_{w,1} \left(a_1 \times n_1 \cdot \Omega_{R}\right) \left(a_1 \times n_1 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} m_{p2} d_2 d_2 \left(\Omega_{R} \times a_2 \cdot \Omega_{R} \times a_2\right) \\
& + \frac{1}{2} I_{a,2} \left(\omega_{n_2} + a_2 \cdot \Omega_{R}\right) \left(\omega_{n_2} + a_2 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_{n,2} \left(n_2 \cdot \Omega_{R}\right) \left(n_2 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_{w,2} \left(a_2 \times n_2 \cdot \Omega_{R}\right) \left(a_2 \times n_2 \cdot \Omega_{R}\right)
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
L = {}& \frac{1}{2} \left(\Omega_{R} \cdot J_b \Omega_{R}\right) \\
& + \frac{1}{2} I_1 \left(\omega_{q_1} + b_1 \cdot \Omega_{R}\right) \left(\omega_{q_1} + b_1 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_2 \left(\omega_{q_2} + b_2 \cdot \Omega_{R}\right) \left(\omega_{q_2} + b_2 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_3 \left(\omega_{q_3} + b_3 \cdot \Omega_{R}\right) \left(\omega_{q_3} + b_3 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} m_{p1} d_1 d_1 \left(\Omega_{R} \times a_1 \cdot \Omega_{R} \times a_1\right) \\
& + \frac{1}{2} I_{a,1} \left(\omega_{n_1} + a_1 \cdot \Omega_{R}\right) \left(\omega_{n_1} + a_1 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_{n,1} \left(n_1 \cdot \Omega_{R}\right) \left(n_1 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_{w,1} \left(a_1 \times n_1 \cdot \Omega_{R}\right) \left(a_1 \times n_1 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} m_{p2} d_2 d_2 \left(\Omega_{R} \times a_2 \cdot \Omega_{R} \times a_2\right) \\
& + \frac{1}{2} I_{a,2} \left(\omega_{n_2} + a_2 \cdot \Omega_{R}\right) \left(\omega_{n_2} + a_2 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_{n,2} \left(n_2 \cdot \Omega_{R}\right) \left(n_2 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_{w,2} \left(a_2 \times n_2 \cdot \Omega_{R}\right) \left(a_2 \times n_2 \cdot \Omega_{R}\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \tau_1 \left(\gamma_{q_1} + b_1 \cdot \eta_{R}\right) \\
& -\tau_1 \left(b_1 \cdot \eta_{R}\right) \\
& + \tau_2 \left(\gamma_{q_2} + b_2 \cdot \eta_{R}\right) \\
& -\tau_2 \left(b_2 \cdot \eta_{R}\right) \\
& + \tau_3 \left(\gamma_{q_3} + b_3 \cdot \eta_{R}\right) \\
& -\tau_3 \left(b_3 \cdot \eta_{R}\right) \\
& + \mu_1 \left(\gamma_{n_1} + a_1 \cdot \eta_{R}\right) \\
& -\mu_1 \left(a_1 \cdot \eta_{R}\right) \\
& + \mu_2 \left(\gamma_{n_2} + a_2 \cdot \eta_{R}\right) \\
& -\mu_2 \left(a_2 \cdot \eta_{R}\right)
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); the equation of $R$ is in the body frame, and the
equations of $q_i$ and $n_k$ are the scalar coefficients of $\gamma_{q_i}$ and $\gamma_{n_k}$.

<!-- model:eom -->
<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\gamma_{q_1}\)</span>, the equation of <span class="math notranslate nohighlight">\(q_1\)</span> (3 terms)</summary>

$$
\begin{aligned}
0 = {}& -I_1 \dot{\omega_{q_1}} \\
& -I_1 \left(b_1 \cdot \dot{\Omega_{R}}\right) \\
& + \tau_1
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\gamma_{q_2}\)</span>, the equation of <span class="math notranslate nohighlight">\(q_2\)</span> (3 terms)</summary>

$$
\begin{aligned}
0 = {}& -I_2 \dot{\omega_{q_2}} \\
& -I_2 \left(b_2 \cdot \dot{\Omega_{R}}\right) \\
& + \tau_2
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\gamma_{q_3}\)</span>, the equation of <span class="math notranslate nohighlight">\(q_3\)</span> (3 terms)</summary>

$$
\begin{aligned}
0 = {}& -I_3 \dot{\omega_{q_3}} \\
& -I_3 \left(b_3 \cdot \dot{\Omega_{R}}\right) \\
& + \tau_3
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\gamma_{n_1}\)</span>, the equation of <span class="math notranslate nohighlight">\(n_1\)</span> (5 terms)</summary>

$$
\begin{aligned}
0 = {}& -I_{a,1} \dot{\omega_{n_1}} \\
& -I_{a,1} \left(a_1 \cdot \dot{\Omega_{R}}\right) \\
& + I_{n,1} \left(a_1 \times n_1 \cdot \Omega_{R}\right) \left(n_1 \cdot \Omega_{R}\right) \\
& -I_{w,1} \left(n_1 \cdot \Omega_{R}\right) \left(a_1 \times n_1 \cdot \Omega_{R}\right) \\
& + \mu_1
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\gamma_{n_2}\)</span>, the equation of <span class="math notranslate nohighlight">\(n_2\)</span> (5 terms)</summary>

$$
\begin{aligned}
0 = {}& -I_{a,2} \dot{\omega_{n_2}} \\
& -I_{a,2} \left(a_2 \cdot \dot{\Omega_{R}}\right) \\
& + I_{n,2} \left(a_2 \times n_2 \cdot \Omega_{R}\right) \left(n_2 \cdot \Omega_{R}\right) \\
& -I_{w,2} \left(n_2 \cdot \Omega_{R}\right) \left(a_2 \times n_2 \cdot \Omega_{R}\right) \\
& + \mu_2
\end{aligned}
$$

</details>

<details>
<summary>Coefficient of <span class="math notranslate nohighlight">\(\eta_{R}\)</span>, the equation of <span class="math notranslate nohighlight">\(R\)</span> (42 terms)</summary>

$$
\begin{aligned}
0 = {}& -J_b \dot{\Omega_{R}} \\
& -\Omega_{R} \times \left(J_b \Omega_{R}\right) \\
& -I_1 \dot{\omega_{q_1}} b_1 \\
& -I_1 \left(b_1 \cdot \dot{\Omega_{R}}\right) b_1 \\
& -I_1 \omega_{q_1} \left(\Omega_{R} \times b_1\right) \\
& -I_1 \left(b_1 \cdot \Omega_{R}\right) \left(\Omega_{R} \times b_1\right) \\
& -I_2 \dot{\omega_{q_2}} b_2 \\
& -I_2 \left(b_2 \cdot \dot{\Omega_{R}}\right) b_2 \\
& -I_2 \omega_{q_2} \left(\Omega_{R} \times b_2\right) \\
& -I_2 \left(b_2 \cdot \Omega_{R}\right) \left(\Omega_{R} \times b_2\right) \\
& -I_3 \dot{\omega_{q_3}} b_3 \\
& -I_3 \left(b_3 \cdot \dot{\Omega_{R}}\right) b_3 \\
& -I_3 \omega_{q_3} \left(\Omega_{R} \times b_3\right) \\
& -I_3 \left(b_3 \cdot \Omega_{R}\right) \left(\Omega_{R} \times b_3\right) \\
& -m_{p1} d_1 d_1 \left(a_1 \times \left(\dot{\Omega_{R}} \times a_1\right)\right) \\
& -m_{p1} d_1 d_1 \left(\Omega_{R} \times \left(a_1 \times \left(\Omega_{R} \times a_1\right)\right)\right) \\
& -I_{a,1} \dot{\omega_{n_1}} a_1 \\
& -I_{a,1} \left(a_1 \cdot \dot{\Omega_{R}}\right) a_1 \\
& -I_{a,1} \omega_{n_1} \left(\Omega_{R} \times a_1\right) \\
& -I_{a,1} \left(a_1 \cdot \Omega_{R}\right) \left(\Omega_{R} \times a_1\right) \\
& -I_{n,1} \omega_{n_1} \left(n_1 \cdot \Omega_{R}\right) \left(a_1 \times n_1\right) \\
& -I_{n,1} \omega_{n_1} \left(a_1 \times n_1 \cdot \Omega_{R}\right) n_1 \\
& -I_{n,1} \left(n_1 \cdot \dot{\Omega_{R}}\right) n_1 \\
& -I_{n,1} \left(n_1 \cdot \Omega_{R}\right) \left(\Omega_{R} \times n_1\right) \\
& + I_{w,1} \omega_{n_1} \left(a_1 \times n_1 \cdot \Omega_{R}\right) n_1 \\
& + I_{w,1} \omega_{n_1} \left(n_1 \cdot \Omega_{R}\right) \left(a_1 \times n_1\right) \\
& -I_{w,1} \left(a_1 \times n_1 \cdot \dot{\Omega_{R}}\right) \left(a_1 \times n_1\right) \\
& -I_{w,1} \left(a_1 \times n_1 \cdot \Omega_{R}\right) \left(\Omega_{R} \times \left(a_1 \times n_1\right)\right) \\
& -m_{p2} d_2 d_2 \left(a_2 \times \left(\dot{\Omega_{R}} \times a_2\right)\right) \\
& -m_{p2} d_2 d_2 \left(\Omega_{R} \times \left(a_2 \times \left(\Omega_{R} \times a_2\right)\right)\right) \\
& -I_{a,2} \dot{\omega_{n_2}} a_2 \\
& -I_{a,2} \left(a_2 \cdot \dot{\Omega_{R}}\right) a_2 \\
& -I_{a,2} \omega_{n_2} \left(\Omega_{R} \times a_2\right) \\
& -I_{a,2} \left(a_2 \cdot \Omega_{R}\right) \left(\Omega_{R} \times a_2\right) \\
& -I_{n,2} \omega_{n_2} \left(n_2 \cdot \Omega_{R}\right) \left(a_2 \times n_2\right) \\
& -I_{n,2} \omega_{n_2} \left(a_2 \times n_2 \cdot \Omega_{R}\right) n_2 \\
& -I_{n,2} \left(n_2 \cdot \dot{\Omega_{R}}\right) n_2 \\
& -I_{n,2} \left(n_2 \cdot \Omega_{R}\right) \left(\Omega_{R} \times n_2\right) \\
& + I_{w,2} \omega_{n_2} \left(a_2 \times n_2 \cdot \Omega_{R}\right) n_2 \\
& + I_{w,2} \omega_{n_2} \left(n_2 \cdot \Omega_{R}\right) \left(a_2 \times n_2\right) \\
& -I_{w,2} \left(a_2 \times n_2 \cdot \dot{\Omega_{R}}\right) \left(a_2 \times n_2\right) \\
& -I_{w,2} \left(a_2 \times n_2 \cdot \Omega_{R}\right) \left(\Omega_{R} \times \left(a_2 \times n_2\right)\right)
\end{aligned}
$$

</details>
<!-- /model:eom -->

## Code

```python
from geomech.systems import get

# L: the bus and wheels of spacecraft_reaction_wheels_s1, plus for each panel k
#   ½ m_pk d_k² |Ω × a_k|² + ½ Ia_k (ω_{n_k} + a_k·Ω)² + ½ In_k (n_k·Ω)² + ½ Iw_k ((a_k × n_k)·Ω)²
# δW: each wheel motor and each drive acts equal and opposite on the bus.
system = get("spacecraft_solar_arrays")
eom = system.derive()
```

The library builder writes $L$ and $\delta W$ in a loop over the wheels and the panels.
The same model is `geomech.systems.get("spacecraft_solar_arrays")`, with default
parameters (1 kg panels of 0.9 m × 0.35 m at $d_k = 0.72$ m), a default initial state and
the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-11}$, with $H(0) = 14.1$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
