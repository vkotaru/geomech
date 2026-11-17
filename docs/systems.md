# Systems library

`geomech.systems` holds known mechanical systems as ready-made objects: the Lagrangian,
the virtual work, the configuration variables and inputs, default parameter values, a
consistent initial state, the source the model comes from and, for most systems, a 3-D
scene. Start from one of them instead of retyping a Lagrangian; every system is verified
by the test suite (below).

```python
from geomech import check_eom, to_standard_form
from geomech.codegen import load_module, module_source
from geomech.simulation import simulate
from geomech.systems import get, list_systems
from geomech.viz import Report

sys = get("double_spherical_pendulum")     # or geomech.systems.double_spherical_pendulum()
sys.L, sys.dW, sys.variables, sys.inputs   # what compute_eom / to_standard_form take
sys.params                                 # default parameters, keyed by printed symbol
sys.x0                                     # a state on the manifold, StateLayout order
sys.reference                              # Reference(source=..., equation=..., note=...)

eom = sys.derive()                         # compute_eom(sys.L, sys.dW, sys.variables)
err = check_eom(sys.L, sys.dW, sys.variables, eom, inputs=sys.inputs)
sf = to_standard_form(eom, sys.variables, sys.inputs)
layout = sys.layout()                      # StateLayout with every parameter
plant = load_module(module_source(sf, layout, lagrangian=sys.L), "plant")
traj = simulate(plant, sys.x0, 2.0, sys.params)
Report("Double spherical pendulum", sys.as_tuple(), eom=eom, standard_form=sf,
       scene=sys.scene(layout), trajectory=traj).write("double_spherical_pendulum.html")
```

`get(name, **kwargs)` passes keyword arguments to the builder; the parametric systems
take a size (`get("n_link_spherical_chain", n_links=5)`,
`get("multi_quadrotor_point_mass", n_quads=4)`). Builders derive nothing: building a
system is cheap, and `sys.derive()` derives it each time it is called. `System` is a
frozen dataclass; `sys.symbols` holds the leaves and point expressions of the model
(`sys.symbols["xQ"]`, the quadrotor position of a load system) for building
controllers or scenes of your own. `sys.scene(layout)` is `None` for a system without a
scene (`spring_mass`).

## The systems

`list_systems()` returns one `SystemInfo` per system (name, configuration, degrees of
freedom at the defaults, reference, builder arguments, one-line summary), in the order
below. A test checks that this table is `geomech.systems.systems_table()`.

<!-- systems-table -->
| Name | Configuration | dof | Source | Closed form compared |
|---|---|---|---|---|
| `point_mass` | R³ | 3 | Lee et al. §3.2 | — |
| `spring_mass` | R | 1 | elementary | — |
| `planar_pendulum` | S¹ | 1 | Lee et al. §4.6.1; general form eq. 4.12 | Lee et al. eq. 4.33 |
| `double_planar_pendulum` | S¹ × S¹ | 2 | Lee et al. §4.6.4 | Lee et al. eq. 4.62 |
| `spherical_pendulum` | S² | 2 | Lee et al. §5.6.1; general form eq. 5.14 | Lee et al. eq. 5.38 |
| `particle_on_rotating_sphere` | S² | 2 | Lee et al. §5.6.2 | Lee et al. eq. 5.43 |
| `double_spherical_pendulum` | S² × S² | 4 | Lee et al. §5.6.5 | Lee et al. eq. 5.71 |
| `n_link_spherical_chain(n_links=…)` | (S²)ⁿ | 6 | Lee et al. §5.3, eq. 5.14; §5.6.5 for n = 2 | — |
| `elastic_spherical_pendulum` | R × S² | 3 | — | — |
| `cart_spherical_pendulum` | R × S² | 3 | — | — |
| `cart_planar_pendulum` | R × S¹ | 2 | elementary | — |
| `free_rigid_body` | SO(3) | 3 | Lee et al. §6.6.1 | Lee et al. eq. 6.16 |
| `heavy_top` | SO(3) | 3 | Lee et al. §6.6.2 | Lee et al. eq. 6.20 |
| `rigid_pendulum` | SO(3) | 3 | Lee et al. §6.6.2; general form eq. 6.8 | Lee et al. eq. 6.20 |
| `double_rigid_pendulum` | SO(3) × SO(3) | 6 | — | — |
| `n_link_rigid_chain(n_links=…)` | (SO(3))ⁿ | 9 | — | — |
| `two_rigid_bodies_spherical_joint` | R³ × SO(3) × SO(3) | 9 | — | — |
| `rigid_body` | R³ × SO(3) | 6 | Lee et al. §7.6.1 | Lee et al. eqs. 7.39–7.40 |
| `spacecraft_reaction_wheels` | SO(3) × R³ | 6 | Lee et al. §9.10 | Lee et al. eqs. 9.61–9.62 |
| `spacecraft_reaction_wheels_s1` | SO(3) × (S¹)³ | 6 | Lee et al. §9.10 | Lee et al. eqs. 9.61–9.62 |
| `spacecraft_solar_arrays` | SO(3) × (S¹)³ × (S¹)² | 8 | — | — |
| `biped_3d` | SO(3) × SO(3) | 6 | — | — |
| `quadrotor` | R³ × SO(3) | 6 | Lee, Leok, McClamroch (2010) | — |
| `quadrotor_offset_com` | R³ × SO(3) | 6 | — | — |
| `quadrotor_suspended_load` | R³ × SO(3) × S² | 8 | Sreenath, Lee, Kumar (2013) | — |
| `quadrotor_elastic_load` | R³ × SO(3) × R × S² | 9 | — | — |
| `quadrotor_flexible_cable(n_links=…)` | R³ × SO(3) × (S²)ⁿ | 12 | — | — |
| `multi_quadrotor_point_mass(n_quads=…)` | R³ × (SO(3) × S²)ⁿ | 18 | — | — |
| `multi_quadrotor_rigid_load(n_quads=…)` | R³ × SO(3) × (SO(3) × S²)ⁿ | 21 | — | — |
| `multi_quadrotor_hose(n_links=…)` | R³ × (S²)ⁿ × SO(3)² | 15 | — | — |
<!-- /systems-table -->

Each system has a model page, listed under [Models](models/index.md), that animates it,
draws it and shows its derived equations: [`point_mass`](models/point_mass.md),
[`spring_mass`](models/spring_mass.md), [`planar_pendulum`](models/planar_pendulum.md),
[`double_planar_pendulum`](models/double_planar_pendulum.md),
[`spherical_pendulum`](models/spherical_pendulum.md),
[`particle_on_rotating_sphere`](models/particle_on_rotating_sphere.md),
[`double_spherical_pendulum`](models/double_spherical_pendulum.md),
[`n_link_spherical_chain`](models/n_link_spherical_chain.md),
[`elastic_spherical_pendulum`](models/elastic_spherical_pendulum.md),
[`cart_spherical_pendulum`](models/cart_spherical_pendulum.md),
[`cart_planar_pendulum`](models/cart_planar_pendulum.md),
[`free_rigid_body`](models/free_rigid_body.md), [`heavy_top`](models/heavy_top.md),
[`rigid_pendulum`](models/rigid_pendulum.md),
[`double_rigid_pendulum`](models/double_rigid_pendulum.md),
[`n_link_rigid_chain`](models/n_link_rigid_chain.md),
[`two_rigid_bodies_spherical_joint`](models/two_rigid_bodies_spherical_joint.md),
[`rigid_body`](models/rigid_body.md),
[`spacecraft_reaction_wheels`](models/spacecraft_reaction_wheels.md),
[`spacecraft_reaction_wheels_s1`](models/spacecraft_reaction_wheels_s1.md),
[`spacecraft_solar_arrays`](models/spacecraft_solar_arrays.md), [`biped_3d`](models/biped_3d.md),
[`quadrotor`](models/quadrotor.md), [`quadrotor_offset_com`](models/quadrotor_offset_com.md),
[`quadrotor_suspended_load`](models/quadrotor_suspended_load.md),
[`quadrotor_elastic_load`](models/quadrotor_elastic_load.md),
[`quadrotor_flexible_cable`](models/quadrotor_flexible_cable.md),
[`multi_quadrotor_point_mass`](models/multi_quadrotor_point_mass.md),
[`multi_quadrotor_rigid_load`](models/multi_quadrotor_rigid_load.md),
[`multi_quadrotor_hose`](models/multi_quadrotor_hose.md).

Each builder's docstring, in the [API reference](api/systems.rst), states the Lagrangian,
the virtual work and the inputs. The conventions are those of the rest of the package
([Conventions](conventions.md)): gravity acts along $-e_3$ ($V = m g\, e_3 \cdot x$), a pendulum or
cable direction $q$ points from the pivot (or quadrotor) to the mass, so $q = -e_3$
hangs down, and $\Omega$ is the body-frame angular velocity. The planar pendulums
(`planar_pendulum`, `double_planar_pendulum`) swing in the vertical plane normal to the
$S^1$ axis $e_2 = (0, 1, 0)$, a parameter (`e2`) like `e3`.

`cart_planar_pendulum` is `cart_spherical_pendulum` with the pendulum held in the rail's
vertical plane: the same parameters ($M$, $m$, $l$, the rail direction `e1`), and $q$ on the
circle normal to $e_2$ (the horizontal normal to that plane, $e_1 \cdot e_2 = 0$). Near upright,
$q = (\sin\varphi, 0, \cos\varphi)$ and $\omega_q = \dot\varphi$. `tests/test_systems_book.py` compares its equations
with the textbook cart-pole in the angle $\varphi$, and `tests/test_s1_manifold.py` checks that
they agree with those of `cart_spherical_pendulum` at in-plane states ($q \cdot e_2 = 0$,
$\omega_q \parallel e_2$) and that the two models simulate alike from an in-plane start.

`double_planar_pendulum` is compared with Lee et al. eq. 4.62 with its velocity terms
corrected: the book prints $m_2 L_1 L_2\, \omega_2(\omega_2 - \omega_1)\, q_1^T S q_2$ and $m_2 L_1 L_2\, \omega_1(\omega_1 - \omega_2)\, q_2^T S q_1$,
which is $\tfrac{d}{dt} \partial \tilde L/\partial \omega_i$ alone; the term $q_i^T S\, \partial \tilde L/\partial q_i$ of eq. 4.12 adds
$m_2 L_1 L_2\, \omega_1 \omega_2\, q_i^T S q_j$ and leaves $\omega_2^2$ and $\omega_1^2$, the textbook double pendulum
(`tests/test_systems_book.py` checks the corrected form against both eq. 4.62 and the
equations in the angles $\theta_i$). `spacecraft_reaction_wheels_s1` is
`spacecraft_reaction_wheels` with the wheel angles on $S^1$; its equations print the same
after renaming $\dot\theta_i \to \omega_i$ and its trajectories agree (`tests/test_s1_manifold.py`).

**Sources.** "Lee et al." is the book *Global Formulations of Lagrangian and Hamiltonian
Dynamics on Manifolds* (`lee2018global` in
[`docs/references.bib`](https://github.com/vkotaru/geomech/blob/main/docs/references.bib)).
"Source" is where the model is treated (`Reference.section`): background reading, not a
comparison. "Closed form compared" (`Reference.equation`) names exactly the equations a
test compares the derived equations with; `n_link_spherical_chain` is compared with
eq. 5.71 at `n_links=2` only. For the quadrotor and the quadrotor with a suspended load the
source is the paper the model is taken from (`lee2010geometric`, `sreenath2013geometric`;
$e_3$ points up here, down in the first); their equations are not compared with the papers.
A "—" in the last column means no closed form is compared: the equations are verified by
the numerical check ([Numeric verification](numeric-verification.md)) and the Lagrangian by
the independent model below.

**Multibody systems.** `spacecraft_reaction_wheels` writes the wheel angles $\theta_i$ as
scalars in $\mathbb{R}$; `spacecraft_reaction_wheels_s1` puts them on $S^1$ as the book does (the angle
is cyclic, so the equations are the same); a test also checks that the inertial angular momentum is
conserved without input. `spacecraft_solar_arrays` adds two solar-array panels to the $S^1$
spacecraft, each a rigid plate with its own mass and inertia turned by a drive about the
bus-fixed axis $\pm b_1$, its normal an $S^1$ point. It is a free-floating model: with no
external force and zero linear momentum, the composite centre of mass is at rest and is a
fixed point of the bus frame (each panel's centre of mass lies on its drive axis), so a
panel's translation is the one the bus rotation induces. Tests check the Lagrangian
against a model of the separate bodies, the conservation of the inertial angular
momentum under wheel and drive torques, and that with massless panels and zero drive
torques its bus and wheel equations are those of `spacecraft_reaction_wheels_s1`. `biped_3d` is only the continuous single-support (swing) phase
of a compass-gait walker: impacts, the switch of the stance foot and the ground contact
are not modelled, so it is not a walking model.

What the scenes can and cannot show: a scene draws expressions of the configuration
variables and parameters, and a shape's size and a cylinder's axis are numbers. The
wheels sit at $0.3\, b_i$ (they follow the parameters) with their discs along the body axes
$e_1$, $e_2$, $e_3$, the default $b_i$; they do not turn with $\theta_i$, since there is no sine or cosine
to build $R \exp(\theta_i \widehat{b_i})$, and an axisymmetric wheel looks the same when it spins anyway.
The solar-array panels do turn: the panel attitude $R P_k$ is built from the panel normal
$n_k$ and the drive axis $a_k$, $P_k = a_k b_1^T + (n_k \times a_k)\, b_2^T + n_k b_3^T$, a rotation when $b_1$, $b_2$,
$b_3$ are the body axes (the default).
The biped's legs and the rigid chain's links are placed from the parameters but drawn with
their default lengths. The biped's ground plane is drawn only.

**Rigid chain.** `n_link_rigid_chain(n_links=2)` has the Lagrangian of
`double_rigid_pendulum` term for term (same symbols); its inputs are joint torques
instead of per-link moments, and a test checks that both give the same equations when the
moments are the joint torques' net moments. Lee et al. §10.2 treat a chain of point
masses on massless links (here `n_link_spherical_chain`), not a chain of rigid bodies;
the rigid chain is checked by `check_eom` and the independent model.

## How the systems are verified

The tests `tests/test_systems*.py` check, for every system at its defaults and every
parametric system also at sizes 1 and 2:

- `compute_eom` succeeds and `check_eom` is at most 1e-5 at three seeds;
- `to_standard_form` gives a `G` block for every input, and `StateLayout` plus
  `generate_module` write an importable module;
- `simulate` runs 1 s from `x0` without input and conserves the energy to 1e-8;
- the scene renders through `Report`;
- the equations of the systems with a closed form in Lee et al. (eqs. 4.33, 4.62, 5.38,
  5.43, 5.71, 6.16, 6.20, 7.39–7.40, 9.61–9.62) agree with it at random admissible states;
- the Lagrangian and the virtual work agree (to 1e-8) with an independent model of the
  system's bodies, written in NumPy from its physical description (positions of the point
  masses and centres of mass, inertias about the centres of mass, heights, springs and
  the points where the forces act), at random admissible states; parametric systems at
  sizes 1, 2 and 3 (`tests/test_systems_physics.py`);
- the default parameters are physical: $e_3 = (0, 0, 1)$, positive masses and lengths, and
  every inertia, shifted to the centre of mass where it is given about a pivot,
  symmetric positive definite with principal moments satisfying the triangle inequality;
- the derived equations and standard-form blocks, the default `params` and `x0` match
  `tests/golden/<name>.json`.

The golden files store each equation's printed string and its value at three fixed
states. A change of value fails as a change of the equations. A change of the printed
form with every value unchanged (for example after a change to the simplifier) also
fails, with a message to check the numbers and regenerate the files:

```bash
python scripts/regenerate_golden.py --check   # report changes, write nothing
python scripts/regenerate_golden.py           # rewrite tests/golden/
```

Every parametric system derives in under 2 s at its default size (a test; the suite prints
the derivation time of every system in its summary). Larger sizes take longer: about
0.3 s for `n_link_rigid_chain(n_links=5)` and 4 s for `n_links=10` on a laptop.

`System` is immutable: `params`, `symbols` and `options` are read-only mappings and every
array in them is a read-only copy owned by the system. To change a parameter, copy the
mapping: `p = dict(sys.params); p["m"] = 2.0`. The sizes `n_links` and `n_quads` run from
1 to 20 (`geomech.systems.MAX_SIZE`); other values raise `ValueError`.

## Not supported

Systems that need a configuration variable on $S^1$ whose sine or cosine enters the
Lagrangian (planar pendulums written with an angle, the Furuta pendulum), nonholonomic
constraints (rolling), closed kinematic chains, contact, or flexible bodies other than
chains of rigid links are not in the library: the engine derives systems on $\mathbb{R}$, $\mathbb{R}^3$, $S^2$,
$SO(3)$ and their products only.
