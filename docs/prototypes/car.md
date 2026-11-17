# Off-road car

```{note}
Prototype: built with geomech, not part of `geomech.systems` and not covered by the test
suite. The terrain, tyre, suspension and drive forces and the driver are prototype code,
passed to the plant as inputs. The parameters are illustrative.
```

## Animation

Nine seconds over rolling hills, a washboard and four rocks, with four-wheel drive
speeding up to 6 m/s and a pure-pursuit driver steering the front wheels. Wheels leave
the ground often. Drag to orbit; the full report (equations, standard form, parameters)
opens in a <a href="../_static/prototypes/car_offroad.html">separate page</a>.

<iframe class="model-report" src="../_static/prototypes/car_offroad.html?view" title="Off-road car, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

The same model on a 15° banked circular track of radius 25 m at 14 m/s
(<a href="../_static/prototypes/car_circle.html">full report</a>):

<iframe class="model-report" src="../_static/prototypes/car_circle.html?view" title="Car on a banked turn, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

Wheel spin and steering are not drawn.

## Model

The chassis has centre of mass $x \in \mathbb{R}^3$ and attitude $R \in SO(3)$. Each wheel
$i$ has a suspension travel $s_i$ along the chassis axis $R e_3$ and a spin angle $\phi_i$
about the axle, the chassis axis $e_2$. The configuration manifold is
$\mathbb{R}^3 \times SO(3) \times \mathbb{R}^4 \times \mathbb{R}^4$, 14 degrees of freedom.
The hub of wheel $i$ is at $x_i = x + R r_i + s_i R e_3$. Each wheel is symmetric about its
axle, so its rotational kinetic energy needs only $\Omega + \dot\phi_i e_2$, and

$$
L = \tfrac12 M \|\dot x\|^2 + \tfrac12 \Omega \cdot J \Omega - M g\, e_3 \cdot x
  + \sum_i \tfrac12 m_w \|\dot x_i\|^2
  + \tfrac12 (\Omega + \dot\phi_i e_2) \cdot J_w (\Omega + \dot\phi_i e_2)
  - m_w g\, e_3 \cdot x_i .
$$

Everything else enters through the virtual work: per wheel a world-frame force $F_i$ at
the hub, a body-frame moment $M_i$ on the chassis, a suspension force $S_i$ and a spin
torque $\tau_i$,

$$
\delta W = \sum_i F_i \cdot \delta x_i + \eta \cdot M_i + S_i\, \delta s_i + \tau_i\, \delta\phi_i ,
$$

32 input components in all. A numeric callback computes them at each step: the terrain
height and normal, a tyre with linear slip and a friction circle ($\mu = 0.9$) that can
leave the ground, a spring-damper suspension with bump stops, and the drive torque.

## Checks

`check_eom`: $2.1 \times 10^{-7}$ (10 equation blocks). From the stored trajectories:

- Off-road: at least one wheel is off the ground in 42 % of the stored states.
- Banked turn, over the last 3 s: the car holds a radius of 24.94 m at 13.85 m/s. The
  ratio of the tyres' lateral to normal force is 0.4268; the banked-curve formula at the
  same speed and radius gives 0.4266. The inner wheels carry down to 0.66 of their static
  load.
