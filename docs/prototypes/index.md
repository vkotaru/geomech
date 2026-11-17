# Prototypes

Larger systems built with geomech: a floating and a walking humanoid, a hand, a six-joint
arm with geared motors, and an off-road car. Each one writes its Lagrangian and virtual
work by hand, derives the equations with `compute_eom`, checks them with `check_eom`,
generates a plant with `geomech.codegen`, simulates it with `geomech.simulation` and draws
it with `geomech.viz`.

```{note}
These are prototypes. They are not part of `geomech.systems`, the test suite does not
cover them, and their parameters are illustrative. Controllers, contact, constraints and
impacts are written in the prototype code, outside geomech; each page says which parts.
The animations are checked-in reports of one run each.
```

| Prototype | Configuration manifold | DOF | Scenario |
| --- | --- | --- | --- |
| [Floating humanoid](humanoid.md) | $\mathbb{R}^3 \times SO(3)^{10}$ | 33 | windmilling arms in zero gravity |
| [Walking humanoid](humanoid_walk.md) | $SO(3)^{11}$ per stance, $\mathbb{R}^3 \times SO(3)^{12}$ at impacts | 33 / 39 | six steps on flat ground |
| [Hand](hand.md) | $\mathbb{R}^3 \times (S^1)^{14}$ | 17 | closing into a cupped pose |
| [Six-joint arm](arm.md) | $SO(3)^6 \times \mathbb{R}^6$, 18 constraints | 6 | pick and place |
| [Off-road car](car.md) | $\mathbb{R}^3 \times SO(3) \times \mathbb{R}^4 \times \mathbb{R}^4$ | 14 | hills and rocks; a banked turn |

```{toctree}
:maxdepth: 1
:hidden:

humanoid
humanoid_walk
hand
arm
car
```
