---
title: 'geomech: Coordinate-free equations of motion on manifolds for robotics'
tags:
  - Python
  - geometric mechanics
  - equations of motion
  - symbolic computation
  - robotics
  - Lie groups
authors:
  - name: Prasanth Kotaru
    orcid: 0000-0001-6608-7469
    affiliation: 1
affiliations:
  - name: Independent Researcher
    index: 1
date: 15 November 2026
bibliography: paper.bib
---

# Summary

The configuration spaces of many robotic and mechanical systems are not Euclidean.
For instance, the attitude of a drone evolves on the rotation group $SO(3)$, and the
direction of a cable or a pendulum evolves on the unit sphere $S^2$. Equations of
motion written directly on these manifolds are compact and free of the singularities
of local coordinates such as Euler angles, and geometric controllers are designed from
them [@lee2010geometric; @lee2018global].

`geomech` is a Python package that derives these equations symbolically. The user
writes the Lagrangian and the virtual work of the applied forces (the work done by the
inputs along an admissible variation of the configuration) in terms of configuration
variables on $\mathbb{R}$, $\mathbb{R}^3$, $S^1$, $S^2$, $SO(3)$ and their products.
The package returns one equation per variable, expressed with rotation matrices, unit
vectors, cross products and hat maps (the map from a vector to the skew-symmetric
matrix of its cross product). For a model of $n$ bodies, such as a chain or a cable of
$n$ links, it can derive the equation of the $i$-th body with $i$ and $n$ symbolic. Each derived equation can be linearized, exported as NumPy [@harris2020numpy] or CasADi code, simulated on
the manifold and used for LQR design.

# Statement of need

Coordinate-free models on $SO(3)$ and $S^2$ are the basis of geometric controllers for
quadrotors [@lee2010geometric] and for quadrotors with cable-suspended loads
[@sreenath2013geometric], and Lee, Leok and McClamroch [@lee2018global] develop the
variational calculus behind them. Deriving such a model by hand is laborious. One
takes the variation of the action, substitutes the constrained variations
($\delta q = \xi \times q$ on $S^2$ and $\delta R = R\hat\eta$ on $SO(3)$), integrates
by parts, and simplifies the result with vector identities. For multi-body systems on
product manifolds, this derivation spans several pages of algebra, and errors are difficult to detect. `geomech` automates this derivation for researchers and
students in robotics and control who model systems such as aerial manipulators,
cable-suspended loads, pendulums and chains of rigid bodies.

# Example

The spherical pendulum shows the workflow. A point mass $m$ hangs from a massless rod
of length $l$ in the direction $q \in S^2$, under gravity and an applied moment $\tau$:

```python
from geomech import *

m, g, l = getScalars("m g l", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
q, tau = S2("q"), Vector("tau")
xi = q.get_variation_vector()          # delta q = xi x q

x = l * q                              # position of the mass
v = x.t_diff()                         # l (omega x q)
L = 0.5 * m * Dot(v, v) - m * g * Dot(x, e3)
dW = Dot(xi, tau)

X = SystemVariables(vectors=[q])
eom = compute_eom(L, dW, X)
check_eom(L, dW, X, eom, inputs=[tau])  # about 3e-8
```

Projected onto the tangent plane at $q$, the derived equation is
$m l^2 \dot\omega + m g l\, q \times e_3 = (I - q q^\top)\tau$, which reduces to Lee et
al. [@lee2018global, Eq. (5.38)] when $\tau = 0$. The same steps apply to the 30
systems of the package's `systems` library, from this pendulum to quadrotors carrying a
rigid load or a flexible hose.

# State of the field

SymPy's mechanics module and PyDy [@meurer2017sympy; @gede2013constrained] derive
symbolic equations by Kane's or Lagrange's method. They cover far more mechanics than
`geomech`, including constraints, many joint types and code generation, but they
formulate the dynamics in generalized coordinates and speeds (Euler angles, quaternions
or body-axis angular velocity components), so the result is a set of scalar equations
rather than the vector equations on $S^2$ and $SO(3)$ used in geometric control.
Commercial tools such as MapleSim and MotionGenesis also derive symbolic multibody
equations in generalized coordinates. Drake [@drake], Pinocchio
[@carpentier2019pinocchio] and MuJoCo [@todorov2012mujoco] evaluate the dynamics of
kinematic trees efficiently, with contact and derivatives, and treat
rotational joints with unit quaternions and tangent-space velocities. Drake and Pinocchio can also evaluate their algorithms with a symbolic
scalar type, which yields the mass matrix and bias terms as symbolic expressions in
generalized coordinates. None of these tools returns the equations in the vector and
matrix form of Lee et al., written with $q \in S^2$, $R \in SO(3)$, $\omega$ and
$\Omega$. CasADi [@andersson2019casadi] provides automatic differentiation and optimal
control, and `geomech` uses it as a code generation target.

| Tool | Symbolic | Vector form | Symbolic $n$ | Linearization | Codegen | Contact |
|---|---|---|---|---|---|---|
| SymPy, PyDy | yes | no | no | yes | yes | no |
| Drake | partial | no | no | yes | yes | yes |
| Pinocchio | partial | no | no | yes | yes | yes |
| MuJoCo | no | no | no | yes | no | yes |
| Bittner and Sreenath | yes | yes | no | no | no | no |
| `geomech` | yes | yes | yes | yes | yes | no |

: Capabilities of tools that derive or evaluate equations of motion, as documented in
November 2026: symbolic equations, equations in the vector form on $S^2$ and $SO(3)$,
symbolic number of bodies, linearization, code generation and contact. *Partial*: the equations are available as symbolic expressions in
generalized coordinates, obtained by evaluating the numerical algorithms with a
symbolic scalar type. \label{tab:compare}

Bittner and Sreenath [@bittner2020symbolic] were the first to automate the derivation
of coordinate-free equations of motion on manifolds, in a Scala package
[@dynamics_on_manifolds_repo] that derives the equations on products of $S^2$ and
$SO(3)$ directly from the Lagrangian and prints them in LaTeX. `geomech` is written in Python, so the derived equations can be passed directly to
numerical tools such as NumPy and CasADi. Its expression tree treats vectors and
matrices as whole objects ($q$, $\omega$, $R$, $\Omega$) and carries the manifold
structure in its leaves, so identities that hold only on the manifold
($\omega \cdot q = 0$, $R^\top R = I$) are rules that state their hypotheses. \autoref{tab:compare} compares the tools discussed in
this section.


`geomech` is limited to $\mathbb{R}$, $\mathbb{R}^3$, $S^1$, $S^2$, $SO(3)$ and their
products. It does not handle contact, constraints other than the manifold structure or
higher-order Lagrangians, nor the Hamiltonian and discrete forms of Lee et al.

# Software design

Every expression in `geomech` is an immutable tree. The leaves are scalars, vectors,
matrices and the manifold points `S1`, `S2` and `SO3`, and Python operators build typed
nodes such as `Dot`, `Cross`, `Hat` and `Variation`. Manifold points carry their own
kinematics: for $q \in S^2$, `q.t_diff()` returns $\omega \times q$ and `q.delta()`
returns $\xi \times q$, so a derivation never differentiates a manifold point directly.

The function `compute_eom` follows the steps of a derivation by hand. It forms the
variation $\delta L + \delta W$, substitutes the variations of the angular velocities
on each manifold, simplifies, integrates by parts, and extracts the coefficient of each
variation vector. Each simplification is a named rule with a stated hypothesis (for
example, $\omega \cdot q = 0$ for a tangent vector at $q$), and every rule applied can
be logged. When a term cannot be resolved, the pipeline raises an error instead of
dropping it, so a forgotten degree of freedom is reported.

The function `check_eom` verifies a derived equation numerically: it forms the
Euler–Lagrange residual of the Lagrangian by finite differences and compares it with
the equation at random states. For parameters of order one, a correct equation agrees
to about $10^{-7}$, and a wrong sign or a missing term gives an error of at least
$10^{-2}$. Each rewrite rule is also tested on random expressions by comparing values
before and after the rewrite.

For each system in the `systems` library, the Lagrangian and virtual work are compared
with an independent NumPy model of its bodies, the equations with the closed-form
equations of Lee et al. where the book gives them, and the printed equations with
stored reference outputs.

A model of $n$ bodies is written with finite sums, such as $x_i = \sum_{j \le i} l_j
q_j$ for a chain. A sum is kept outside the expression tree as a list of summands with
their bound indices and order constraints, so every summand passes through
`compute_eom` unchanged. The coefficient of
the variation of the $j$-th body under a bound index $j$ becomes a term of the equation
of body $i$ by the substitution $j := i$, which turns the constraints on $j$ into bounds
on the remaining sums. For the $n$-link spherical chain this yields the inertia
coefficient $l_i l_j \sum_{k \ge \max(i,j)} m_k$, and the generic equation agrees with
Lee et al. [@lee2018global, Eq. (5.17)] for every $n$ tested. For every such model the
tests expand the generic equation for $n = 1, \dots, 6$ and check it with `check_eom`.

Downstream modules use only the derived equations. The `linearize` function computes the
variation-based linearization of Wu and Sreenath [@wu2015variation] in state-space form;
`codegen` writes standalone NumPy or CasADi modules; `simulation` integrates them on the
Lie group with a fourth-order Runge–Kutta–Munthe-Kaas method [@munthekaas1999high] and
designs LQR controllers about equilibria; and `viz` renders a 3-D animation through the
same generated code. The repository's control examples design LQR controllers for 24
regulation problems, 19 of them with a stabilizable linearization, and for two
tracking problems.

# Research impact statement

`geomech` derives models across robotics and mechanics. In aerial transport, it
reproduces the models of a quadrotor carrying a load on an elastic or flexible cable
and of several quadrotors carrying a flexible hose, which the author previously derived
by hand [@kotaru2017elastic; @kotaru2018flexible; @kotaru2020hose]; the flexible cable
is also derived for a symbolic number of links. Other examples include the
three-dimensional rigid pendulum [@siravuru2020nonlinear], double rigid pendulums and
chains of rigid links, spacecraft with reaction wheels and solar arrays
[@lee2018global], and the swing phase of a three-dimensional biped. The package is
available on PyPI with documentation, control examples and an agent skill for deriving
new models.

# AI usage disclosure

Much of the code, tests, automated code reviews and documentation written in 2026 was
produced with Anthropic's Claude models through Claude Code, under the author's direction. This paper was drafted with Claude and revised by the author. The original 2021
implementation was written without AI assistance. The author made the design
decisions, reviewed the changes, and is responsible for the content. 

# Acknowledgements

The author thanks Brian Bittner, Koushil Sreenath and the Hybrid Robotics Lab for their
earlier work on the symbolic computation of dynamics on manifolds and for their Scala
source code.

# References
