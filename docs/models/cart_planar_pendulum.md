# Cart with a planar pendulum

The planar cart-pole: a cart of mass $M$ on a straight rail along $e_1$ carrying a
pendulum (point mass $m$, rod length $l$) that swings in the vertical plane of the rail.
The configuration is the cart position $s \in \mathbb{R}$ and the rod direction $q$ on
the circle $S^1$ normal to $e_2$, so the manifold is $\mathbb{R} \times S^1$ and the
upright and hanging positions need no special treatment. Not taken from the literature;
the equations are verified numerically (`check_eom`). Library system
`cart_planar_pendulum`.

## Animation

Ten seconds, unforced, from the library's initial state: the rod $0.6$ rad ($34°$) from
upright, $\omega_q = 0.5$ rad/s, and the cart moving at $\dot s = 0.3$ m/s. The pendulum
falls, swings through the hanging position to about $34°$ from upright on the other
side and back, a full swing every $2.9$ s, without going over the top. The cart moves back and forth against the
swing while drifting along the rail, to $s = 3.8$ m at the end. Drag to orbit; the full
report (equations, standard form, parameters) opens in a
<a href="../model-reports/cart_planar_pendulum.html">separate page</a>.

<iframe class="model-report" src="../model-reports/cart_planar_pendulum.html?view" title="Cart with a planar pendulum, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 400 300" role="img" aria-labelledby="cpp-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="cpp-title">Cart with a planar pendulum: cart position s, velocity ṡ, cart mass M, force f, rod direction q, angular velocity ω_q about e2, length l, bob m, gravity</title>
  <defs>
    <marker id="cpp-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- rail -->
    <line x1="20" y1="230" x2="380" y2="230" stroke-width="2"/>
    <!-- cart -->
    <rect x="150" y="196" width="80" height="30" stroke-width="1.6"/>
    <circle cx="165" cy="226" r="4"/>
    <circle cx="215" cy="226" r="4"/>
    <!-- upright reference -->
    <line x1="190" y1="196" x2="190" y2="66" stroke-dasharray="3 4" opacity="0.5"/>
    <!-- cart position extent -->
    <line x1="60" y1="234" x2="60" y2="262" stroke-dasharray="3 4" opacity="0.6"/>
    <line x1="190" y1="234" x2="190" y2="262" stroke-dasharray="3 4" opacity="0.6"/>
    <!-- inertial axes; e2 points into the page -->
    <circle cx="40" cy="80" r="7"/>
    <path d="M35,75 L45,85 M45,75 L35,85"/>
    <line x1="40" y1="73" x2="40" y2="38" marker-end="url(#cpp-arr)"/>
    <line x1="47" y1="80" x2="82" y2="80" marker-end="url(#cpp-arr)"/>
    <!-- rod -->
    <line x1="190" y1="196" x2="255" y2="83.4" stroke-width="2.2"/>
  </g>
  <!-- s: cart position along the rail -->
  <line x1="60" y1="252" x2="190" y2="252" stroke="#3f8fa0" stroke-width="2.4" marker-end="url(#cpp-arr)"/>
  <!-- cart velocity -->
  <line x1="236" y1="211" x2="276" y2="211" stroke="#3f8fa0" stroke-width="2" marker-end="url(#cpp-arr)"/>
  <!-- q: unit vector along the rod -->
  <line x1="190" y1="196" x2="215" y2="152.7" stroke="#3f8fa0" stroke-width="3" marker-end="url(#cpp-arr)"/>
  <!-- omega: rotation about e2, clockwise in this view -->
  <path d="M 161.8,185.7 A 30 30 0 0 1 184.8,166.5" fill="none" stroke="#3f8fa0" stroke-width="2" marker-end="url(#cpp-arr)"/>
  <!-- force on the cart -->
  <line x1="100" y1="211" x2="146" y2="211" stroke="#e08a2c" stroke-width="2" marker-end="url(#cpp-arr)"/>
  <!-- gravity on the bob -->
  <line x1="255" y1="97" x2="255" y2="150" stroke="#e08a2c" stroke-width="2" marker-end="url(#cpp-arr)"/>
  <circle cx="60" cy="230" r="3.5" fill="currentColor"/>
  <circle cx="190" cy="196" r="3.5" fill="currentColor"/>
  <circle cx="255" cy="83.4" r="12" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="190" y="217" text-anchor="middle">M</text>
    <text x="42" y="220">O</text>
    <text x="120" y="274" fill="#3f8fa0">s</text>
    <text x="282" y="216" fill="#3f8fa0">ṡ</text>
    <text x="116" y="200" fill="#e08a2c">f</text>
    <text x="220" y="178" fill="#3f8fa0">q</text>
    <text x="141" y="168" fill="#3f8fa0">ω<tspan font-size="12" dy="4">q</tspan></text>
    <text x="214" y="110">l</text>
    <text x="275" y="80">m</text>
    <text x="262" y="145" fill="#e08a2c">m g</text>
    <text x="47" y="44">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="87" y="85">e<tspan font-size="12" dy="4">1</tspan></text>
    <text x="10" y="86">e<tspan font-size="12" dy="4">2</tspan></text>
  </g>
</svg>
<figcaption>The cart is at <span class="math notranslate nohighlight">\(s\, e_1\)</span> along the rail, measured from the origin O, and
the bob at <span class="math notranslate nohighlight">\(s\, e_1 + l q\)</span>. The rod direction <span class="math notranslate nohighlight">\(q\)</span> lies in the plane normal
to <span class="math notranslate nohighlight">\(e_2\)</span>, which points into the page, so <span class="math notranslate nohighlight">\(\dot q = \omega_q\, e_2 \times q\)</span> turns the rod
clockwise in this picture when <span class="math notranslate nohighlight">\(\omega_q &gt; 0\)</span>. Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span>; the force
<span class="math notranslate nohighlight">\(f\)</span> acts on the cart along the rail, and the pendulum is unactuated.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $s$ | cart position along the rail | $\mathbb{R}$ |
| $\dot s$ | cart velocity | $\mathbb{R}$ |
| $\delta s$ | variation of $s$ | $\mathbb{R}$ |
| $q$ | direction of the rod, cart to bob, $q \cdot e_2 = 0$ | $S^1$ |
| $\omega_q$ | angular velocity of the rod about $e_2$, $\dot q = \omega_q\, e_2 \times q$ | $\mathbb{R}$ |
| $\gamma_q$ | variation of $q$: $\delta q = \gamma_q\, e_2 \times q$ | $\mathbb{R}$ |
| $f$ | force on the cart along the rail | $\mathbb{R}$ |
| $M, m, l, g$ | cart mass, bob mass, rod length, gravity | constants |
| $e_1, e_2, e_3$ | rail direction, axis of the plane of motion, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} M \dot{s} \dot{s} \\
& + \frac{1}{2} m \left(\left(\dot{s} e_1 + l \omega_{q} \left(e_2 \times q\right)\right) \cdot \left(\dot{s} e_1 + l \omega_{q} \left(e_2 \times q\right)\right)\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& \left(\left(M + m\right) g s \left(e_1 \cdot e_3\right) + m g l \left(q \cdot e_3\right)\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} M \dot{s} \dot{s} \\
& + \frac{1}{2} m \left(\left(\dot{s} e_1 + l \omega_{q} \left(e_2 \times q\right)\right) \cdot \left(\dot{s} e_1 + l \omega_{q} \left(e_2 \times q\right)\right)\right) \\
& -\left(\left(M + m\right) g s \left(e_1 \cdot e_3\right) + m g l \left(q \cdot e_3\right)\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \delta s f
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); both are scalar equations, the coefficients of
$\delta s$ and of the scalar variation $\gamma_q$.

<!-- model:eom -->
Coefficient of $\delta s$, the equation of $s$:

$$
\begin{aligned}
0 = {}& -M \ddot{s} \\
& -m \ddot{s} \\
& -m \dot{\omega_{q}} l \left(e_1 \cdot e_2 \times q\right) \\
& + m \omega_{q} l \omega_{q} \left(e_1 \cdot q\right) \\
& -M g \left(e_1 \cdot e_3\right) \\
& -m g \left(e_1 \cdot e_3\right) \\
& + f
\end{aligned}
$$

Coefficient of $\gamma_{q}$, the equation of $q$:

$$
\begin{aligned}
0 = {}& -m l \ddot{s} \left(e_2 \times q \cdot e_1\right) \\
& -m l \dot{\omega_{q}} l \\
& -m g l \left(e_2 \times q \cdot e_3\right)
\end{aligned}
$$
<!-- /model:eom -->

## Code

```python
from geomech import S1, Dot, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
M, m, g, l = getScalars("M m g l", attr=["Constant"])
e1 = Vector("e1", attr=["Constant", "UnitNorm"])
e2 = Vector("e2", attr=["Constant", "UnitNorm"])
e3 = Vector("e3", attr=["Constant"])
s, q, f = Scalar("s"), S1("q", axis=e2), Scalar("f")

cart = e1 * s
bob = cart + q * l
sd, vb = s.t_diff(), bob.t_diff()
PE = (M + m) * g * s * Dot(e1, e3) + m * g * l * Dot(q, e3)
L = half * M * sd * sd + half * m * Dot(vb, vb) - PE
dW = s.delta() * f
eom = compute_eom(L, dW, SystemVariables(scalars=[s], vectors=[q]))
```

The same model is `geomech.systems.get("cart_planar_pendulum")`, with default
parameters (a horizontal rail, $e_1 = (1, 0, 0)$; the rail may be inclined), a default
initial state and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-9}$, with $H(0) = 3.44$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
