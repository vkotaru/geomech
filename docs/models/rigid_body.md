# Rigid body in space

A rigid body moving freely in space under gravity, with six degrees of freedom, driven by
a force $F$ at its centre of mass and a moment $M$ in the body frame. The configuration
is the position $x \in \mathbb{R}^3$ of the centre of mass and the attitude
$R \in SO(3)$, with body angular velocity $\Omega_{R}$ and $\dot R = R \widehat{\Omega_{R}}$;
SE(3) is modelled as the product $\mathbb{R}^3 \times SO(3)$. No Euler angles are used, so
the equations have no singular attitude. Lee et al. §7.6.1, eqs. 7.39–7.40. Library
system `rigid_body`.

## Animation

Ten seconds, forced, from a start that is not the library's initial state: the
library's start, but drifting at 0.1 m/s along $e_1$ instead of thrown upward at
$\dot x = (0.5, 0, 3)$ m/s. A constant force $F = m g e_3$ at the centre of mass cancels
gravity, and $M = 0$. The body starts 1 m up, tumbling with
$\Omega_{R} = (0.3, -1.2, 0.8)$ rad/s. The centre of mass moves 1 m along
$e_1$ at constant height. The rotation is that of a free rigid body: about 2.4
revolutions in the 10 s, with $\lvert \Omega_{R} \rvert$ at most 1.52 rad/s. Drag to
orbit; the full report (equations, standard form, parameters) opens in a
<a href="../model-reports/rigid_body.html">separate page</a>.

<iframe class="model-report" src="../model-reports/rigid_body.html?view" title="Rigid body in space, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 280" role="img" aria-labelledby="rb-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="rb-title">Rigid body in space: position x, velocity, attitude, angular velocity Omega_R, mass m, inertia J, force F, moment M and gravity</title>
  <defs>
    <marker id="rb-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- box, oblique view -->
    <path d="M170,110 L290,110 L290,170 L170,170 Z" stroke-width="1.8"/>
    <path d="M170,110 L210,80 L330,80 L290,110 M290,170 L330,140 L330,80" stroke-width="1.8"/>
    <!-- inertial axes -->
    <line x1="50" y1="240" x2="50" y2="195" marker-end="url(#rb-arr)"/>
    <line x1="50" y1="240" x2="95" y2="240" marker-end="url(#rb-arr)"/>
  </g>
  <!-- position of the centre of mass -->
  <line x1="50" y1="240" x2="247" y2="127" stroke="#3f8fa0" stroke-width="2.4" marker-end="url(#rb-arr)"/>
  <!-- velocity -->
  <line x1="250" y1="125" x2="345" y2="45" stroke="#3f8fa0" stroke-width="2" marker-end="url(#rb-arr)"/>
  <!-- angular velocity -->
  <line x1="250" y1="125" x2="190" y2="40" stroke="#3f8fa0" stroke-width="3" marker-end="url(#rb-arr)"/>
  <!-- force at the centre of mass, gravity, moment -->
  <line x1="250" y1="125" x2="340" y2="190" stroke="#e08a2c" stroke-width="2" marker-end="url(#rb-arr)"/>
  <line x1="250" y1="125" x2="250" y2="235" stroke="#e08a2c" stroke-width="2" marker-end="url(#rb-arr)"/>
  <path d="M 150,150 A 22 22 0 0 1 150,112" fill="none" stroke="#e08a2c" stroke-width="2" marker-end="url(#rb-arr)"/>
  <circle cx="50" cy="240" r="3" fill="currentColor"/>
  <circle cx="250" cy="125" r="3.5" fill="currentColor"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="156" y="200" fill="#3f8fa0">x</text>
    <text x="350" y="48" fill="#3f8fa0">ẋ</text>
    <text x="156" y="38" fill="#3f8fa0">Ω<tspan font-size="12" dy="4">R</tspan></text>
    <text x="182" y="135">m, J</text>
    <text x="346" y="202" fill="#e08a2c">F</text>
    <text x="258" y="232" fill="#e08a2c">m g</text>
    <text x="114" y="136" fill="#e08a2c">M</text>
    <text x="44" y="188">e<tspan font-size="12" dy="4">3</tspan></text>
    <text x="100" y="245">e<tspan font-size="12" dy="4">1</tspan></text>
  </g>
</svg>
<figcaption>The centre of mass is at <span class="math notranslate nohighlight">\(x\)</span> in the inertial frame and moves with the
inertial velocity <span class="math notranslate nohighlight">\(\dot x\)</span>; the attitude <span class="math notranslate nohighlight">\(R\)</span> maps body axes to inertial axes.
<span class="math notranslate nohighlight">\(\Omega_R\)</span> and the moment <span class="math notranslate nohighlight">\(M\)</span> are body-frame vectors; the force <span class="math notranslate nohighlight">\(F\)</span>
acts at the centre of mass in the inertial frame. Gravity acts along <span class="math notranslate nohighlight">\(-e_3\)</span>.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $x$ | position of the centre of mass, inertial frame | $\mathbb{R}^3$ |
| $\dot x$ | velocity of the centre of mass, inertial frame | $\mathbb{R}^3$ |
| $\delta x$ | variation of $x$ | $\mathbb{R}^3$ |
| $R$ | attitude, body frame to inertial frame | $SO(3)$ |
| $\Omega_{R}$ | body angular velocity, $\dot R = R \widehat{\Omega_{R}}$ | $\mathbb{R}^3$ |
| $\eta_{R}$ | variation of $R$: $\delta R = R \widehat{\eta_{R}}$ | $\mathbb{R}^3$ |
| $F$ | input force at the centre of mass, inertial frame | $\mathbb{R}^3$ |
| $M$ | input moment, body frame | $\mathbb{R}^3$ |
| $J$ | inertia about the centre of mass, body frame | constant, symmetric |
| $m, g, e_3$ | mass, gravity, vertical unit vector | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m \left(\dot{x} \cdot \dot{x}\right) \\
& + \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right)
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& m g \left(e_3 \cdot x\right)
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m \left(\dot{x} \cdot \dot{x}\right) \\
& + \frac{1}{2} \left(\Omega_{R} \cdot J \Omega_{R}\right) \\
& -m g \left(e_3 \cdot x\right)
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \delta x \cdot F \\
& + \eta_{R} \cdot M
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); the equation of $x$ is written in the inertial
frame and the equation of $R$ in the body frame.

<!-- model:eom -->
Coefficient of $\delta x$, the equation of $x$:

$$
\begin{aligned}
0 = {}& -m \ddot{x} \\
& -m g e_3 \\
& + F
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

Multiplied by $-1$, these are Lee et al. eqs. 7.39–7.40 with the inputs added:
$m \ddot x + m g e_3 = F$ and $J \dot\Omega_R + \Omega_R \times J\Omega_R = M$.

## Code

```python
from geomech import SO3, Dot, Matrix, Scalar, SystemVariables, Vector, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
m, g = getScalars("m g", attr=["Constant"])
J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
e3 = Vector("e3", attr=["Constant"])
x, R, F, M = Vector("x"), SO3("R"), Vector("F"), Vector("M")
Om, eta = R.get_tangent_vector(), R.get_variation_vector()
v = x.t_diff()

L = half * m * Dot(v, v) + half * Dot(Om, J * Om) - m * g * Dot(e3, x)
dW = Dot(x.delta(), F) + Dot(eta, M)
eom = compute_eom(L, dW, SystemVariables(vectors=[x], matrices=[R]))
```

The same model is `geomech.systems.get("rigid_body")`, with default parameters, a
default initial state (thrown upward at $\dot x = (0.5, 0, 3)$ m/s, not the drift
animated above) and the 3-D scene used above.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-6}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-7}$, with $H(0) = 19.9$ |
| Closed form, Lee et al. eqs. 7.39–7.40 | equal numerically at random states (`tests/test_systems_book.py`) |
<!-- /model:checks -->
