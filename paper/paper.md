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
date: 14 November 2026
bibliography: paper.bib
---

# Summary

Many robots and mechanical systems have configurations that are not lists of
independent numbers. The attitude of a drone is a rotation, an element of the
rotation group $SO(3)$; the direction of a taut cable or a pendulum is a unit vector,
a point on the sphere $S^2$. `geomech` is a Python package that derives the equations
of motion of such systems symbolically and without coordinates. The user writes a
Lagrangian (kinetic minus potential energy) and the virtual work of the applied forces
in terms of configuration variables on $\mathbb{R}$, $\mathbb{R}^3$, $S^1$, $S^2$, $SO(3)$ and
their products ($SE(3)$ is modelled as $\mathbb{R}^3 \times SO(3)$), and `compute_eom`
returns one vector equation per variable, expressed with rotation matrices, unit vectors,
cross products and hat maps, valid globally and free of the singularities of Euler
angles.

Around this core the package provides: an independent numerical check, `check_eom`, that tests every derived
equation against a finite-difference Euler–Lagrange residual computed without the
library's own rewrite rules; the standard form $M a + f + G u = 0$; the variation-based
linearization $M\,\delta\dot x = A\,\delta x + B\,\delta u$ along an arbitrary trajectory
[@wu2015variation]; generation of standalone NumPy [@harris2020numpy] or CasADi
[@andersson2019casadi] modules; simulation with a fourth-order Runge–Kutta–Munthe-Kaas
integrator on the Lie group [@munthekaas1999high] and LQR design from the
linearization; reports in a single HTML file with a 3-D animation, the equations and the
verification numbers; a library of 29 ready-made systems, from a spherical pendulum, a
spacecraft with reaction wheels and a 3-D biped (single-support phase) to teams of
quadrotors carrying a rigid load; and an agent skill that walks a coding
assistant through the whole procedure. NumPy is the only required dependency; SciPy and CasADi are optional.

# Statement of need

Coordinate-free models on $SO(3)$ and $S^2$ underpin geometric controllers for
quadrotors [@lee2010geometric] and quadrotors with cable-suspended loads
[@sreenath2013geometric], and the systematic treatment of Lee, Leok and McClamroch
[@lee2018global] gives the variational calculus behind them. Euler angles and other
minimal coordinates are singular somewhere on $SO(3)$, so equations written in them hold
only locally; working on the manifold avoids this but moves the burden to the
derivation. By hand, one takes the variation of the action, substitutes the constrained
variations ($\delta q = \xi \times q$ on $S^2$, $\delta R = R\hat\eta$ on $SO(3)$),
integrates by parts, and simplifies with cross-product and orthogonality identities. For
a multi-body system on a product manifold this runs to pages of algebra, and the errors
are hard to see: a dropped factor in a gyroscopic coupling still yields a plausible
equation.

`geomech` is for robotics and dynamics researchers and students who model aerial
manipulators, cable-suspended loads, pendulums and chains of rigid bodies on these
manifolds, and who want equations they can read, publish, linearize for control design
and turn into simulation or optimization code, with a numeric check that the equations
are right.

# State of the field

SymPy's mechanics module and PyDy [@meurer2017sympy; @gede2013constrained] derive
symbolic equations by Kane's or Lagrange's method in generalized coordinates, with
orientations parameterized by coordinates the user chooses (Euler angles, quaternions);
they cover far more mechanics (constraints, many joint types, code generation) than
`geomech`, but do not produce equations on $S^2$ and $SO(3)$. Drake
[@drake], Pinocchio [@carpentier2019pinocchio] and MuJoCo [@todorov2012mujoco] are
numerical: they evaluate the dynamics of kinematic trees quickly and at scale, with
contact and standard model formats (all three) and analytical derivatives (Pinocchio),
none of which `geomech` offers; they do not return the closed-form, coordinate-free equation
from which geometric controllers are designed. CasADi [@andersson2019casadi] is a
symbolic framework for automatic differentiation and optimal control with no notion of
mechanics; `geomech` uses it as a code generation target.

The method `geomech` implements is that of Lee, Leok and McClamroch [@lee2018global].
Bittner and Sreenath [@bittner2020symbolic] first automated it, in Scala; their paper
describes the derivation of the equations, and their public implementation
[@dynamics_on_manifolds_repo] also prints them in LaTeX and contains a routine that
varies the equations about a desired trajectory. `geomech` is a Python implementation
of the same approach that adds a numeric check of every equation, the variation-based
linearization of Wu and Sreenath [@wu2015variation] in the state-space form used for
LQR design, code generation, Lie-group simulation and reporting. It is limited
to $\mathbb{R}$, $\mathbb{R}^3$, $S^1$, $S^2$, $SO(3)$ and their products; it does not handle
contact, constraints other than the manifold structure, or higher-order Lagrangians.

# Software design

Every expression is an immutable tree. Leaves are scalars, vectors, matrices and the
manifold points `S2` and `SO3`; Python operators build typed nodes (n-ary sums, binary
products chosen by operand type, `Dot`, `Cross`, `Hat`, `Vee`, `Transpose`, and the
calculus nodes `TimeDerivative` and `Variation`), and a mistyped operand raises an error
naming both operands. Manifold points are eager: `q.t_diff()` returns $\omega \times q$
and `q.delta()` returns $\xi \times q$, so a tree never contains the derivative of a
manifold point.

`compute_eom` runs a fixed pipeline: the variation $\delta S = \delta L + \delta W$; the
manifold rules that replace the variation of each tangent vector; simplification to a
fixed point; integration by parts against the time derivative of each variation vector,
after `collect` brings that derivative to the left of a dot product; expansion; and
extraction of the coefficient of each variation vector by tree truncation. Extraction
raises a typed error instead of returning zero for a term it cannot resolve, and the
pipeline raises if the variation of an undeclared symbol is left over, so a forgotten
degree of freedom is reported rather than dropped. An optional trace object keeps every
stage, so each can be rerun by hand.

Every rewrite is a named rule (for example `vector.dot_tangent_orthogonal`,
$\omega \cdot q = 0$ for a tangent vector at $q$), listed with the hypothesis it needs
in its module's docstring and logged on demand. Tests check that the documented and
implemented rule names agree and that every rewrite of random trees logs a rule.

Correctness rests on two numeric checks, not on the rules alone. The check `check_eom`
evaluates the Lagrangian numerically and forms the Euler–Lagrange residual by central
differences (on $S^2$ through the constrained Lagrangian and compared on the tangent
plane; on $SO(3)$ in the body-frame form [@lee2018global]), then compares it with the
derived equation at random states and accelerations; a correct equation scores about
$10^{-7}$ and a wrong sign or dropped term $2 \times 10^{-2}$ or more for parameters of
order one. Every pinned system of the test suite passes it; those taken from Lee et al.
are also pinned as strings and compared numerically with the textbook equation.
Separately, a fuzzer applies every simplification and manifold rewrite, the product
rules of $\delta$ and $d/dt$, and the coefficient extraction to random well-typed trees
(200 per campaign in its quick mode, 3000 by default), evaluated as truncated Taylor
series in time and variation, and compares input and output.

The $S^2$ variation rule shows why both are needed. The variation of the angular
velocity, from Lee et al. eq. (5.13), contains a term $q(q \cdot \dot\xi)$ that is easy
to drop as "normal to the sphere"; on the constraint set it equals $\omega \times \xi$,
giving $\delta\omega = \dot\xi - 2\,\omega \times \xi$. An earlier version of `geomech`
used $\dot\xi - \omega \times \xi$. The two agree whenever $\partial L/\partial\omega$ is
orthogonal to $q$, as for every pendulum example, so string tests pinned against the
textbook passed; `check_eom` exposed that gyroscopic and coupling terms otherwise came
out at half strength. Such systems are now in the test suite.

The `systems` library makes these checks reusable. Each of its 29 systems carries its
Lagrangian, virtual work, default parameters, an initial state and its source, and the
test suite checks each one with `check_eom`; its Lagrangian and virtual work against
an independent numeric model of its bodies (positions, kinetic energy by finite
differences, heights, and the work of forces at their points of application), which
`check_eom` alone cannot do because it only shows that the equations follow from the stated
Lagrangian; the closed forms of Lee et al. where the book gives one; and golden
reference outputs (equations stored as strings and as values at fixed states) that fail
on any change of value.

Downstream modules consume only the symbolic output. `codegen` maps symbols to a
documented state layout and writes modules that need only NumPy (or CasADi), and the
tests compare the generated functions with direct evaluation of the tree on every
pinned system.
`simulation` integrates generated modules on the manifold and designs LQR gains from
the generated linearization, which the repository's control examples apply to 24
equilibria of the library systems (22 stabilizable) and two tracking problems; `viz` evaluates scene geometry through the same generated
code, so the animation has no kinematics of its own.

# Research impact statement

The first use of `geomech` is reproducing the coordinate-free models the author
derived by hand in earlier work on aerial transport: a quadrotor with a load on an
elastic cable [@kotaru2017elastic] and on a flexible cable modelled as a chain of links
[@kotaru2018flexible], which the library now derives and verifies as
`quadrotor_elastic_load` and `quadrotor_flexible_cable`. Deriving and checking such
models was the bottleneck in that work, and the multi-quadrotor flexible-hose model of
[@kotaru2020hose] is built from the same pieces (rigid bodies and chains of links on
$S^2$). The library's systems and the control examples in the repository are the
current evidence that the package does what it claims; they are verification, not
research results.

# AI usage disclosure

Substantial parts of the code, tests, code reviews and documentation written in 2026
were produced with Anthropic's Claude models (Opus and Sonnet families), used through
Claude Code, under the author's direction; commit trailers in the repository record
this. This paper was drafted with Claude and revised by the author. The original 2021
implementation was written without AI assistance.

The author made the design decisions, reviewed the changes, and is responsible for the
content. Correctness of the derived equations does not rest on review alone: the
equation of every pinned system is checked by the test suite against the numerical
Euler–Lagrange check `check_eom`, the textbook systems also against the equations of Lee et al.
[@lee2018global], and the rewrite rules are checked by numeric fuzzing.

# Acknowledgements

This work received no external funding.

# References
