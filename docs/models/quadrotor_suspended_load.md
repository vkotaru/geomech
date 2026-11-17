# Quadrotor with a suspended load

A quadrotor carrying a point-mass load $m_L$ on a cable of length $l$. The cable is taut
and massless, so the quadrotor sits at $x_Q = x_L - l q$, where $q \in S^2$ points from
the quadrotor to the load. The configuration manifold is
$\mathbb{R}^3 \times SO(3) \times S^2$: load position $x_L$, attitude $R$, cable direction
$q$. The thrust $f$ acts along the body axis $R e_3$ at the quadrotor. Model from
Sreenath, Lee and Kumar (2013); the equations are verified numerically (`check_eom`).
Library system `quadrotor_suspended_load`.

## Animation

Six seconds with constant thrust equal to the total weight and zero moments: open loop,
with no controller. The start is not the library's initial state: the quadrotor starts
level and not rotating, while the load keeps the library's initial state, 29° off the
vertical and swinging.
With $M = 0$ and the thrust through the centre of mass, the attitude equation does not
involve the load, so a vehicle that starts level and at rest stays exactly level. The load
swings through the vertical, up to 33° to either side, mostly in one vertical plane.
Nothing holds position: the system keeps the momentum of its initial velocity and drifts
5.1 m sideways while the load rises from 0.5 m to 4.6 m. Drag to orbit; the full report
(equations, standard form, parameters) opens in a
<a href="../model-reports/quadrotor_suspended_load.html">separate page</a>.

<iframe class="model-report" src="../model-reports/quadrotor_suspended_load.html?view" title="Quadrotor with a suspended load, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 310" role="img" aria-labelledby="qsl-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="qsl-title">Quadrotor with a suspended load: quadrotor at x_Q with attitude R, thrust f and moment M; cable direction q from the quadrotor to the load, length l; load m_L at x_L</title>
  <defs>
    <marker id="qsl-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <line x1="40" y1="290" x2="40" y2="250" marker-end="url(#qsl-arr)"/>
    <line x1="40" y1="290" x2="80" y2="290" marker-end="url(#qsl-arr)"/>
  </g>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <line x1="150" y1="80" x2="240" y2="230" stroke-width="1.6"/>
  </g>
  <!-- quadrotor at (150,80) -->
  <g fill="none" stroke="currentColor" stroke-width="1.6">
    <line x1="100" y1="80" x2="200" y2="80"/>
    <line x1="122" y1="94" x2="178" y2="66"/>
    <ellipse cx="100" cy="80" rx="17" ry="5" stroke-width="1.2"/>
    <ellipse cx="200" cy="80" rx="17" ry="5" stroke-width="1.2"/>
    <ellipse cx="178" cy="66" rx="17" ry="5" stroke-width="1.2"/>
    <ellipse cx="122" cy="94" rx="17" ry="5" stroke-width="1.2"/>
  </g>
  <rect x="141" y="75" width="18" height="10" rx="2" fill="currentColor"/>
  <line x1="150.0" y1="72.0" x2="150.0" y2="26.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#qsl-arr)"/>
  <path d="M 137,46.7 A 13 5 0 1 0 163,46.7" fill="none" stroke="#e08a2c" stroke-width="1.6" marker-end="url(#qsl-arr)"/>
  <!-- q: from the quadrotor to the load -->
  <line x1="150.0" y1="80.0" x2="175.7" y2="122.9" stroke="#3f8fa0" stroke-width="3" marker-end="url(#qsl-arr)"/>
  <line x1="240.0" y1="242.0" x2="240.0" y2="290.0" stroke="#e08a2c" stroke-width="2" marker-end="url(#qsl-arr)"/>
  <circle cx="225" cy="160" r="7" fill="none" stroke="#3f8fa0" stroke-width="1.6"/>
  <circle cx="225" cy="160" r="2" fill="#3f8fa0"/>
  <circle cx="240" cy="230" r="12" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="158" y="32" fill="#e08a2c">f</text>
    <text x="171" y="47" fill="#e08a2c">M</text>
    <text x="226" y="86" fill="#3f8fa0">R, Ω<tspan font-size="12" dy="4">R</tspan></text>
    <text x="106" y="122">m<tspan font-size="12" dy="4">Q</tspan><tspan dy="-4">, J</tspan></text>
    <text x="64" y="86">x<tspan font-size="12" dy="4">Q</tspan></text>
    <text x="186" y="118" fill="#3f8fa0">q</text>
    <text x="237" y="165" fill="#3f8fa0">ω<tspan font-size="12" dy="4">q</tspan></text>
    <text x="258" y="236">m<tspan font-size="12" dy="4">L</tspan></text>
    <text x="208" y="236" fill="#3f8fa0">x<tspan font-size="12" dy="4">L</tspan></text>
    <text x="248" y="284" fill="#e08a2c">m<tspan font-size="12" dy="4">L</tspan><tspan dy="-4"> g</tspan></text>
    <text x="34" y="242">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="86" y="295">e<tspan font-size="12" dy="4">1</tspan></text>
    <text x="175" y="170">l</text>
  </g>
</svg>
<figcaption>The quadrotor is at <span class="math notranslate nohighlight">\(x_Q = x_L - l q\)</span>; <span class="math notranslate nohighlight">\(q\)</span> points from the
quadrotor to the load, so a hanging cable has <span class="math notranslate nohighlight">\(q = -e_3\)</span>. The cable's angular velocity
<span class="math notranslate nohighlight">\(\omega_q\)</span> is perpendicular to <span class="math notranslate nohighlight">\(q\)</span> (out of the page here) and <span class="math notranslate nohighlight">\(\dot q = \omega_q \times q\)</span>. The thrust
<span class="math notranslate nohighlight">\(f\)</span> acts along <span class="math notranslate nohighlight">\(R e_3\)</span>; <span class="math notranslate nohighlight">\(M\)</span> is the moment in the body frame. Gravity acts
along <span class="math notranslate nohighlight">\(-e_3\)</span> on both masses (shown on the load only).</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $x_L$ | position of the load | $\mathbb{R}^3$ |
| $\dot x_L$ | velocity of the load | $\mathbb{R}^3$ |
| $\delta x_L$ | variation of $x_L$ | $\mathbb{R}^3$ |
| $q$ | cable direction, quadrotor to load | $S^2$ |
| $\omega_{q}$ | angular velocity of the cable, $\dot q = \omega_{q} \times q$, $q \cdot \omega_{q} = 0$ | $T_q S^2$ |
| $\xi_{q}$ | variation of $q$: $\delta q = \xi_{q} \times q$ | $T_q S^2$ |
| $R$ | attitude of the quadrotor, body frame to inertial frame | $SO(3)$ |
| $\Omega_{R}$ | body angular velocity, $\dot R = R \widehat{\Omega_{R}}$ | $\mathbb{R}^3$ |
| $\eta_{R}$ | variation of $R$: $\delta R = R \widehat{\eta_{R}}$ | $\mathbb{R}^3$ |
| $f$ | thrust along $R e_3$ | $\mathbb{R}$ |
| $M$ | moment, body frame | $\mathbb{R}^3$ |
| $m_Q, J, m_L, l, g, e_3$ | quadrotor mass and inertia, load mass, cable length, gravity, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m_Q \left(\left(\dot{x_L} -l \left(\omega_{q} \times q\right)\right) \cdot \left(\dot{x_L} -l \left(\omega_{q} \times q\right)\right)\right) \\
& + \frac{1}{2} m_L \left(\dot{x_L} \cdot \dot{x_L}\right) \\
& + \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& m_Q g \left(e_3 \cdot \left(x_L -l q\right)\right) \\
& + m_L g \left(e_3 \cdot x_L\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m_Q \left(\left(\dot{x_L} -l \left(\omega_{q} \times q\right)\right) \cdot \left(\dot{x_L} -l \left(\omega_{q} \times q\right)\right)\right) \\
& + \frac{1}{2} m_L \left(\dot{x_L} \cdot \dot{x_L}\right) \\
& + \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right) \\
& -m_Q g \left(e_3 \cdot \left(x_L -l q\right)\right) \\
& -m_L g \left(e_3 \cdot x_L\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \left(\delta x_L -l \left(\xi_{q} \times q\right)\right) \cdot f R e_3 \\
& + \eta_{R} \cdot M
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); the $R$ equation is in the body frame and the $q$
equation holds on the tangent plane at $q$.

<!-- model:eom -->
Coefficient of $\delta x_L$, the equation of $x_L$:

$$
\begin{aligned}
0 = {}& -m_Q \ddot{x_L} \\
& + m_Q l \left(\dot{\omega_{q}} \times q\right) \\
& -\left(\omega_{q} \cdot \omega_{q}\right) m_Q l q \\
& -m_L \ddot{x_L} \\
& -m_Q g e_3 \\
& -m_L g e_3 \\
& + f R e_3
\end{aligned}
$$

Coefficient of $\xi_{q}$, the equation of $q$:

$$
\begin{aligned}
0 = {}& m_Q l \left(\left(\omega_{q} \times q\right) \times \dot{x_L}\right) \\
& + m_Q l \left(q \times \ddot{x_L}\right) \\
& -m_Q l l \left(q \times \left(\dot{\omega_{q}} \times q\right)\right) \\
& -2 m_Q l \left(\omega_{q} \times \left(q \times \dot{x_L}\right)\right) \\
& + m_Q l \left(q \times \left(\omega_{q} \times \dot{x_L}\right)\right) \\
& + m_Q g l \left(q \times e_3\right) \\
& -l f \left(q \times \left(R e_3\right)\right)
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
<!-- /model:eom -->

## Code

```python
from geomech import S2, SO3, Cross, Dot, Matrix, Scalar, SystemVariables, Vector, compute_eom
from geomech import getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
mQ, mL, g, l = getScalars("m_Q m_L g l", attr=["Constant"])
J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
e3 = Vector("e3", attr=["Constant"])
xL, q, R = Vector("x_L"), S2("q"), SO3("R")
xi, Om, eta = q.get_variation_vector(), R.get_tangent_vector(), R.get_variation_vector()
f, M = Scalar("f"), Vector("M")

xQ = xL - q * l  # the quadrotor hangs above the load
vQ, vL = xQ.t_diff(), xL.t_diff()
L = (
    half * mQ * Dot(vQ, vQ)
    + half * mL * Dot(vL, vL)
    + half * Dot(Om, J * Om)
    - mQ * g * Dot(e3, xQ)
    - mL * g * Dot(e3, xL)
)
dW = Dot(xL.delta() - Cross(xi, q) * l, (R * e3) * f) + Dot(eta, M)  # thrust acts at x_Q
eom = compute_eom(L, dW, SystemVariables(vectors=[xL, q], matrices=[R]))
```

The same model is `geomech.systems.get("quadrotor_suspended_load")`, with default
parameters, a default initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 6 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-9}$, with $H(0) = 23.1$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
