# Spacecraft with reaction wheels

A rigid spacecraft (the bus) carrying three reaction wheels that spin about the body axes
$b_1, b_2, b_3$. The configuration is the attitude of the bus, $R \in SO(3)$, and the three
wheel angles relative to the bus, $\theta_i \in \mathbb{R}$: the space is
$SO(3) \times \mathbb{R}^3$. The wheel motors are internal, so they turn the bus only by
exchanging angular momentum with the wheels. Lee et al. §9.10, eqs. 9.61–9.62 (a special
case: three symmetric wheels on the body axes, no translation). Library system
`spacecraft_reaction_wheels`.

## Animation

Ten seconds from the library's initial state, with a constant motor torque
$\tau_1 = 0.02$ N·m on wheel 1 and none on the others. The bus starts at
$\Omega_R = (0.05, -0.1, 0.08)$ rad/s and the wheels at 30, −20 and 10 rad/s relative to it.
Wheel 1 speeds up steadily, from 30 to 40.1 rad/s; wheels 2 and 3 change by less than
0.2 rad/s. The bus takes the reaction: its rate about $b_1$ goes from 0.05 to −0.03 rad/s
(it does not fall monotonically, because the three axes are coupled), and the bus ends 24.5°
from its initial attitude, against 64° in an unforced run. The wheels are drawn but not turned: a
scene cannot rotate a body by $\theta_i$, and a symmetric wheel looks the same at every angle
(the [S¹ version](spacecraft_reaction_wheels_s1.md) marks the rim). Drag to orbit; the full
report (equations, standard form, parameters) opens in a
<a href="../model-reports/spacecraft_reaction_wheels.html">separate page</a>.

<iframe class="model-report" src="../model-reports/spacecraft_reaction_wheels.html?view" title="Spacecraft with reaction wheels, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 300" role="img" aria-labelledby="srw-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="srw-title">Spacecraft with reaction wheels: bus with body axes b1, b2, b3, attitude R, angular velocity Omega_R, and a wheel on each axis with angle theta and motor torque tau</title>
  <defs>
    <marker id="srw-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- bus: front face, top face, right face -->
    <rect x="130" y="130" width="100" height="80" stroke-width="1.6"/>
    <path d="M130,130 L174,108 L274,108 L230,130 M274,108 L274,188 L230,210" stroke-width="1.6"/>
    <!-- body axes from the centre of the bus (dashed inside it) -->
    <path d="M202,159 L202,119 M202,159 L252,159 M202,159 L180,170" stroke-dasharray="3 4" opacity="0.6"/>
    <line x1="202" y1="119" x2="202" y2="28" marker-end="url(#srw-arr)"/>
    <line x1="252" y1="159" x2="360" y2="159" marker-end="url(#srw-arr)"/>
    <line x1="180" y1="170" x2="70" y2="225" marker-end="url(#srw-arr)"/>
    <!-- wheels, one on each axis -->
    <ellipse cx="202" cy="72" rx="26" ry="10" stroke-width="2"/>
    <ellipse cx="305" cy="159" rx="10" ry="24" stroke-width="2"/>
    <ellipse cx="100" cy="210" rx="20" ry="24" stroke-width="2"/>
  </g>
  <!-- body angular velocity -->
  <line x1="202" y1="159" x2="120" y2="70" stroke="#3f8fa0" stroke-width="2.2" marker-end="url(#srw-arr)"/>
  <circle cx="202" cy="159" r="3" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="212" y="36">b<tspan font-size="12" dy="4">3</tspan></text>
    <text x="336" y="146">b<tspan font-size="12" dy="4">2</tspan></text>
    <text x="40" y="245">b<tspan font-size="12" dy="4">1</tspan></text>
    <text x="80" y="74" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">R</tspan></text>
    <text x="282" y="214" fill="#3f8fa0">R</text>
    <text x="118" y="252">I<tspan font-size="12" dy="4">1</tspan></text>
    <text x="326" y="190">I<tspan font-size="12" dy="4">2</tspan></text>
    <text x="152" y="78">I<tspan font-size="12" dy="4">3</tspan></text>
    <text x="195" y="198">J<tspan font-size="12" dy="4">b</tspan></text>
    <text x="242" y="80"><tspan fill="#3f8fa0">θ<tspan font-size="12" dy="4">3</tspan></tspan><tspan dy="-4">, </tspan><tspan fill="#e08a2c">τ<tspan font-size="12" dy="4">3</tspan></tspan></text>
    <text x="290" y="116"><tspan fill="#3f8fa0">θ<tspan font-size="12" dy="4">2</tspan></tspan><tspan dy="-4">, </tspan><tspan fill="#e08a2c">τ<tspan font-size="12" dy="4">2</tspan></tspan></text>
    <text x="26" y="190"><tspan fill="#3f8fa0">θ<tspan font-size="12" dy="4">1</tspan></tspan><tspan dy="-4">, </tspan><tspan fill="#e08a2c">τ<tspan font-size="12" dy="4">1</tspan></tspan></text>
  </g>
</svg>
<figcaption>The bus with its body axes <span class="math notranslate nohighlight">\(b_1\)</span>, <span class="math notranslate nohighlight">\(b_2\)</span>, <span class="math notranslate nohighlight">\(b_3\)</span> from its centre of mass;
<span class="math notranslate nohighlight">\(R\)</span> takes body coordinates to inertial ones and <span class="math notranslate nohighlight">\(\Omega_R\)</span> is the angular
velocity in the body frame. Wheel <span class="math notranslate nohighlight">\(i\)</span> sits on <span class="math notranslate nohighlight">\(b_i\)</span> and turns by
<span class="math notranslate nohighlight">\(\theta_i\)</span> relative to the bus; its motor applies <span class="math notranslate nohighlight">\(\tau_i\)</span> to the wheel
about <span class="math notranslate nohighlight">\(b_i\)</span> and <span class="math notranslate nohighlight">\(-\tau_i\)</span> to the bus. The animation also draws two
solar panels fixed to the bus; they are drawn with the bus and are not separate bodies in
this model.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $R$ | attitude of the bus, body frame to inertial frame | $SO(3)$ |
| $\Omega_{R}$ | body angular velocity of the bus, $\dot R = R \widehat{\Omega_{R}}$ | $\mathbb{R}^3$ |
| $\eta_{R}$ | variation of $R$: $\delta R = R \widehat{\eta_{R}}$ | $\mathbb{R}^3$ |
| $\theta_i$ | angle of wheel $i$ relative to the bus, $i = 1, 2, 3$ | $\mathbb{R}$ |
| $\dot\theta_i$ | spin rate of wheel $i$ relative to the bus | $\mathbb{R}$ |
| $\delta\theta_i$ | variation of $\theta_i$ | $\mathbb{R}$ |
| $\tau_i$ | motor torque on wheel $i$ about $b_i$ (the bus gets $-\tau_i b_i$) | $\mathbb{R}$ |
| $J_b$ | inertia of the spacecraft about its centre of mass with the wheels locked, minus the wheels' axial inertias | constant |
| $I_i, b_i$ | axial inertia of wheel $i$, its unit spin axis in the body frame | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} \left(\Omega_{R} \cdot J_b \Omega_{R}\right) \\
& + \frac{1}{2} I_1 \left(\dot{\theta_1} + b_1 \cdot \Omega_{R}\right) \left(\dot{\theta_1} + b_1 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_2 \left(\dot{\theta_2} + b_2 \cdot \Omega_{R}\right) \left(\dot{\theta_2} + b_2 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_3 \left(\dot{\theta_3} + b_3 \cdot \Omega_{R}\right) \left(\dot{\theta_3} + b_3 \cdot \Omega_{R}\right)
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
& + \frac{1}{2} I_1 \left(\dot{\theta_1} + b_1 \cdot \Omega_{R}\right) \left(\dot{\theta_1} + b_1 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_2 \left(\dot{\theta_2} + b_2 \cdot \Omega_{R}\right) \left(\dot{\theta_2} + b_2 \cdot \Omega_{R}\right) \\
& + \frac{1}{2} I_3 \left(\dot{\theta_3} + b_3 \cdot \Omega_{R}\right) \left(\dot{\theta_3} + b_3 \cdot \Omega_{R}\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \tau_1 \left(\delta \theta_1 + b_1 \cdot \eta_{R}\right) \\
& -\tau_1 \left(b_1 \cdot \eta_{R}\right) \\
& + \tau_2 \left(\delta \theta_2 + b_2 \cdot \eta_{R}\right) \\
& -\tau_2 \left(b_2 \cdot \eta_{R}\right) \\
& + \tau_3 \left(\delta \theta_3 + b_3 \cdot \eta_{R}\right) \\
& -\tau_3 \left(b_3 \cdot \eta_{R}\right)
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); the equation of $R$ is in the body frame.

<!-- model:eom -->
Coefficient of $\delta \theta_1$, the equation of $\theta_1$:

$$
\begin{aligned}
0 = {}& -I_1 \ddot{\theta_1} \\
& -I_1 \left(b_1 \cdot \dot{\Omega_{R}}\right) \\
& + \tau_1
\end{aligned}
$$

Coefficient of $\delta \theta_2$, the equation of $\theta_2$:

$$
\begin{aligned}
0 = {}& -I_2 \ddot{\theta_2} \\
& -I_2 \left(b_2 \cdot \dot{\Omega_{R}}\right) \\
& + \tau_2
\end{aligned}
$$

Coefficient of $\delta \theta_3$, the equation of $\theta_3$:

$$
\begin{aligned}
0 = {}& -I_3 \ddot{\theta_3} \\
& -I_3 \left(b_3 \cdot \dot{\Omega_{R}}\right) \\
& + \tau_3
\end{aligned}
$$

Coefficient of $\eta_{R}$, the equation of $R$:

$$
\begin{aligned}
0 = {}& -J_b \dot{\Omega_{R}} \\
& -\Omega_{R} \times \left(J_b \Omega_{R}\right) \\
& -I_1 \ddot{\theta_1} b_1 \\
& -I_1 \left(b_1 \cdot \dot{\Omega_{R}}\right) b_1 \\
& -I_1 \dot{\theta_1} \left(\Omega_{R} \times b_1\right) \\
& -I_1 \left(b_1 \cdot \Omega_{R}\right) \left(\Omega_{R} \times b_1\right) \\
& -I_2 \ddot{\theta_2} b_2 \\
& -I_2 \left(b_2 \cdot \dot{\Omega_{R}}\right) b_2 \\
& -I_2 \dot{\theta_2} \left(\Omega_{R} \times b_2\right) \\
& -I_2 \left(b_2 \cdot \Omega_{R}\right) \left(\Omega_{R} \times b_2\right) \\
& -I_3 \ddot{\theta_3} b_3 \\
& -I_3 \left(b_3 \cdot \dot{\Omega_{R}}\right) b_3 \\
& -I_3 \dot{\theta_3} \left(\Omega_{R} \times b_3\right) \\
& -I_3 \left(b_3 \cdot \Omega_{R}\right) \left(\Omega_{R} \times b_3\right)
\end{aligned}
$$
<!-- /model:eom -->

Multiplied by $-1$, with $J' = J_b + \sum_i I_i b_i b_i^T$ (the inertia with the wheels
locked), $\alpha_i = I_i$, $s_i = b_i$, $\omega_i = \dot\theta_i$ and no external moment,
these are Lee et al. eqs. 9.61–9.62:
$J' \dot\Omega_R + \sum_i I_i \ddot\theta_i b_i + \Omega_R \times (J' \Omega_R + \sum_i I_i \dot\theta_i b_i) = 0$
and $I_i\, b_i \cdot \dot\Omega_R + I_i \ddot\theta_i = \tau_i$.

## Code

```python
from geomech import SO3, Dot, Matrix, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
R = SO3("R")
Om, eta = R.get_tangent_vector(), R.get_variation_vector()
Jb = Matrix("J_b", attr=["Constant", "SymmetricMatrix"])
Is = getScalars("I_1 I_2 I_3", attr=["Constant"])
bs = [Vector(f"b{i}", attr=["Constant", "UnitNorm"]) for i in (1, 2, 3)]
thetas = [Scalar(f"\\theta_{i}") for i in (1, 2, 3)]
taus = [Scalar(f"\\tau_{i}") for i in (1, 2, 3)]

L = half * Dot(Om, Jb * Om)
dW = None
for I_i, b, th, tau in zip(Is, bs, thetas, taus):
    spin = th.t_diff() + Dot(b, Om)  # absolute spin rate of the wheel about b
    L = L + half * I_i * spin * spin
    term = tau * (th.delta() + Dot(b, eta)) - tau * Dot(b, eta)  # on the wheel, on the bus
    dW = term if dW is None else dW + term
eom = compute_eom(L, dW, SystemVariables(scalars=thetas, matrices=[R]))
```

The same model is `geomech.systems.get("spacecraft_reaction_wheels")`, with default
parameters, a default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-12}$, with $H(0) = 14.1$ |
| Closed form, Lee et al. eqs. 9.61–9.62 | equal numerically at random states (`tests/test_systems_book.py`) |
<!-- /model:checks -->
