# Spring–mass

A mass $m$ on a linear spring of stiffness $k$, moving along a line, with an applied
force $f$. The configuration is the displacement $x \in \mathbb{R}$ from the spring's
rest position: a scalar configuration variable, the harmonic oscillator
$m \ddot x + k x = f$. Not taken from the literature; the equations are verified
numerically (`check_eom`). Library system `spring_mass`.

## Animation

Ten seconds, unforced, from the library's initial state $x = 0.4$, $\dot x = 0.3$: the
mass oscillates along $e_1$ with amplitude $0.44$ and period $2\pi\sqrt{m/k} = 3.85$ s
(default $m = 1.5$, $k = 4$). The wall and the direction $e_1$ are drawing aids only;
the model is one-dimensional. Drag to orbit; the full report (equations, standard form,
parameters) opens in a <a href="../model-reports/spring_mass.html">separate page</a>.

<iframe class="model-report" src="../model-reports/spring_mass.html?view" title="Spring–mass, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Diagram

<figure class="model-diagram">
<svg viewBox="0 0 380 260" role="img" aria-labelledby="sm-title" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%;height:auto;color:inherit">
  <title id="sm-title">Spring–mass: wall, spring of stiffness k, mass m, displacement x, velocity ẋ and applied force f</title>
  <defs>
    <marker id="sm-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.2">
    <!-- wall and floor -->
    <line x1="40" y1="100" x2="40" y2="200" stroke-width="2"/>
    <path d="M40,105 l-10,10 M40,122 l-10,10 M40,139 l-10,10 M40,156 l-10,10 M40,173 l-10,10 M40,190 l-10,10"/>
    <line x1="40" y1="200" x2="350" y2="200" stroke-width="2"/>
    <path d="M60,200 l-10,10 M85,200 l-10,10 M110,200 l-10,10 M135,200 l-10,10 M160,200 l-10,10 M185,200 l-10,10 M210,200 l-10,10 M235,200 l-10,10 M260,200 l-10,10 M285,200 l-10,10 M310,200 l-10,10 M335,200 l-10,10"/>
    <!-- spring -->
    <path d="M40,170 L60,170 l7.5,-10 l15,20 l15,-20 l15,20 l15,-20 l15,20 l15,-20 l15,20 l15,-20 l7.5,10 L210,170"/>
    <!-- mass -->
    <rect x="210" y="140" width="60" height="60" stroke-width="1.6"/>
    <!-- rest position and displacement extent -->
    <line x1="150" y1="114" x2="150" y2="236" stroke-dasharray="3 4" opacity="0.6"/>
    <line x1="210" y1="204" x2="210" y2="236" stroke-dasharray="3 4" opacity="0.6"/>
  </g>
  <!-- x: displacement from the rest position -->
  <line x1="150" y1="228" x2="210" y2="228" stroke="#3f8fa0" stroke-width="2.4" marker-end="url(#sm-arr)"/>
  <!-- velocity -->
  <line x1="220" y1="122" x2="262" y2="122" stroke="#3f8fa0" stroke-width="2" marker-end="url(#sm-arr)"/>
  <!-- applied force -->
  <line x1="270" y1="170" x2="330" y2="170" stroke="#e08a2c" stroke-width="2" marker-end="url(#sm-arr)"/>
  <g fill="currentColor" font-family="serif" font-style="italic" font-size="17">
    <text x="150" y="104" text-anchor="middle" font-size="15">x = 0</text>
    <text x="118" y="146">k</text>
    <text x="240" y="177" text-anchor="middle">m</text>
    <text x="176" y="251" fill="#3f8fa0">x</text>
    <text x="270" y="127" fill="#3f8fa0">ẋ</text>
    <text x="337" y="175" fill="#e08a2c">f</text>
  </g>
</svg>
<figcaption>The displacement <span class="math notranslate nohighlight">\(x\)</span> is measured from the position where the spring is
at rest (dashed), positive to the right; <span class="math notranslate nohighlight">\(\dot x\)</span> is the velocity and <span class="math notranslate nohighlight">\(f\)</span> the
applied force along the line.</figcaption>
</figure>

| Symbol | Meaning | Space |
|---|---|---|
| $x$ | displacement of the mass from the spring's rest position | $\mathbb{R}$ |
| $\dot x$ | velocity | $\mathbb{R}$ |
| $\delta x$ | variation of $x$ | $\mathbb{R}$ |
| $f$ | applied force along the line | $\mathbb{R}$ |
| $m, k$ | mass, spring stiffness | constants |

## Energy, Lagrangian and virtual work

<!-- model:energy -->
Kinetic energy (the terms of $L$ with a velocity):

$$
\begin{aligned}
T = {}& \frac{1}{2} m \dot{x} \dot{x}
\end{aligned}
$$

Potential energy (the other terms, with the sign changed):

$$
\begin{aligned}
V = {}& \frac{1}{2} k x x
\end{aligned}
$$

Lagrangian, $L = T - V$, as the library holds it:

$$
\begin{aligned}
L = {}& \frac{1}{2} m \dot{x} \dot{x} \\
& -\frac{1}{2} k x x
\end{aligned}
$$

Virtual work of the inputs:

$$
\begin{aligned}
\delta W = {}& \delta x f
\end{aligned}
$$
<!-- /model:energy -->

## Equations of motion

Derived by `compute_eom`. Each equation is $E = 0$ with
$E = -(\text{Euler–Lagrange}) + \text{forcing}$
([conventions](../conventions.md)); here it is one scalar equation.

<!-- model:eom -->
Coefficient of $\delta x$, the equation of $x$:

$$
\begin{aligned}
0 = {}& -m \ddot{x} \\
& -k x \\
& + f
\end{aligned}
$$
<!-- /model:eom -->

## Code

```python
from geomech import Scalar, SystemVariables, compute_eom, getScalars

half = Scalar("0.5", value=0.5, attr=["Constant"])
m, k = getScalars("m k", attr=["Constant"])
x, f = Scalar("x"), Scalar("f")
v = x.t_diff()

L = half * m * v * v - half * k * x * x
dW = x.delta() * f
eom = compute_eom(L, dW, SystemVariables(scalars=[x]))
```

The same model is `geomech.systems.get("spring_mass")`, with default parameters and a
default initial state.

## Checks

<!-- model:checks -->
| Check | Result |
|---|---|
| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random states | below $10^{-7}$ |
| Energy $H = T + V$, 10 s unforced from the animation's initial state, RKMK4 with the default step: $\max_t \lvert H(t) - H(0) \rvert$ | below $10^{-13}$, with $H(0) = 0.388$ |
| Closed form | none compared; the equations are verified numerically (`check_eom`) |
<!-- /model:checks -->
