# Hand

```{note}
Prototype: built with geomech, not part of `geomech.systems` and not covered by the test
suite. The joint controller is prototype code; the parameters are illustrative.
```

## Animation

Three seconds. PD joint torques close the four fingers and the thumb from almost straight
into a cupped pose, reached smoothly by 1.6 s, while a constant upward force equal to the
weight of the whole hand holds up the palm. Drag to orbit; the full report (equations,
standard form, parameters) opens in a
<a href="../_static/prototypes/hand.html">separate page</a>.

<iframe class="model-report" src="../_static/prototypes/hand.html?view" title="Hand, interactive report" loading="lazy" style="width:100%;height:560px;border:1px solid var(--color-background-border);border-radius:6px"></iframe>

## Model

A palm of mass $M$ on $\mathbb{R}^3$ (it translates and does not rotate), four fingers of
three links and a thumb of two links. Every link is a thin rod on $S^1$: its direction
$q$ turns in the flexion plane of its digit, about the plane's unit normal $a$, with
scalar angular velocity $\omega$ ($\dot q = \omega\, a \times q$). The thumb's plane is
tilted 40° from the fingers'. The configuration manifold is
$\mathbb{R}^3 \times (S^1)^{14}$, 17 degrees of freedom.

The first joint of a digit is fixed to the palm at $x + b$; link $k$ has its centre at the
joint plus $\tfrac12 l_k q_k$ and ends at the joint plus $l_k q_k$, where the next link
starts. With $x_k$ the centre of link $k$,

$$
L = \tfrac12 M \|\dot x\|^2 - M g\, e_3 \cdot x
  + \sum_k \tfrac12 m_k \|\dot x_k\|^2 + \tfrac12 I_k \omega_k^2 - m_k g\, e_3 \cdot x_k ,
$$

and the virtual work of a force $F$ on the palm and a torque $\tau_k$ at each joint,
acting between link $k$ and the link before it, is

$$
\delta W = F \cdot \delta x + \sum_k \tau_k \,(\gamma_k - \gamma_{k-1}),
$$

with $\gamma_k$ the variation of $q_k$ and $\gamma_0 = 0$ for the first link of a digit.

## Checks

From the run behind the animation: `check_eom` $1.9 \times 10^{-7}$ (15 equation blocks).
