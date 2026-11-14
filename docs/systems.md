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
| `spherical_pendulum` | S² | 2 | Lee et al. §5.6.1; general form eq. 5.14 | Lee et al. eq. 5.38 |
| `particle_on_rotating_sphere` | S² | 2 | Lee et al. §5.6.2 | Lee et al. eq. 5.43 |
| `double_spherical_pendulum` | S² × S² | 4 | Lee et al. §5.6.5 | Lee et al. eq. 5.71 |
| `n_link_spherical_chain(n_links=…)` | (S²)ⁿ | 6 | Lee et al. §5.3, eq. 5.14; §5.6.5 for n = 2 | — |
| `elastic_spherical_pendulum` | R × S² | 3 | — | — |
| `cart_spherical_pendulum` | R × S² | 3 | — | — |
| `free_rigid_body` | SO(3) | 3 | Lee et al. §6.6.1 | Lee et al. eq. 6.16 |
| `heavy_top` | SO(3) | 3 | Lee et al. §6.6.2 | Lee et al. eq. 6.20 |
| `rigid_pendulum` | SO(3) | 3 | Lee et al. §6.6.2; general form eq. 6.8 | Lee et al. eq. 6.20 |
| `double_rigid_pendulum` | SO(3) × SO(3) | 6 | — | — |
| `two_rigid_bodies_spherical_joint` | R³ × SO(3) × SO(3) | 9 | — | — |
| `rigid_body` | R³ × SO(3) | 6 | Lee et al. §7.6.1 | Lee et al. eqs. 7.39–7.40 |
| `quadrotor` | R³ × SO(3) | 6 | Lee, Leok, McClamroch (2010) | — |
| `quadrotor_offset_com` | R³ × SO(3) | 6 | — | — |
| `quadrotor_suspended_load` | R³ × SO(3) × S² | 8 | Sreenath, Lee, Kumar (2013) | — |
| `quadrotor_elastic_load` | R³ × SO(3) × R × S² | 9 | — | — |
| `quadrotor_flexible_cable(n_links=…)` | R³ × SO(3) × (S²)ⁿ | 12 | — | — |
| `multi_quadrotor_point_mass(n_quads=…)` | R³ × (SO(3) × S²)ⁿ | 18 | — | — |
| `multi_quadrotor_rigid_load(n_quads=…)` | R³ × SO(3) × (SO(3) × S²)ⁿ | 21 | — | — |
<!-- /systems-table -->

Each builder's docstring, in the [API reference](api/systems.rst), states the Lagrangian,
the virtual work and the inputs. The conventions are those of the rest of the package
([Conventions](conventions.md)): gravity acts along −e₃ (`V = m g e₃·x`), a pendulum or
cable direction `q` points from the pivot (or quadrotor) to the mass, so `q = −e₃`
hangs down, and Ω is the body-frame angular velocity.

**Sources.** "Lee et al." is the book *Global Formulations of Lagrangian and Hamiltonian
Dynamics on Manifolds* (`lee2018global` in
[`docs/references.bib`](https://github.com/vkotaru/geomech/blob/main/docs/references.bib)).
"Source" is where the model is treated (`Reference.section`): background reading, not a
comparison. "Closed form compared" (`Reference.equation`) names exactly the equations a
test compares the derived equations with; `n_link_spherical_chain` is compared with
eq. 5.71 at `n_links=2` only. For the quadrotor and the quadrotor with a suspended load the
source is the paper the model is taken from (`lee2010geometric`, `sreenath2013geometric`;
e₃ points up here, down in the first); their equations are not compared with the papers.
A "—" in the last column means no closed form is compared: the equations are verified by
the numeric oracle ([Numeric verification](numeric-verification.md)) and the Lagrangian by
the independent model below.

## How the systems are verified

The tests `tests/test_systems*.py` check, for every system at its defaults and every
parametric system also at sizes 1 and 2:

- `compute_eom` succeeds and `check_eom` is at most 1e-5 at three seeds;
- `to_standard_form` gives a `G` block for every input, and `StateLayout` plus
  `generate_module` write an importable module;
- `simulate` runs 1 s from `x0` without input and conserves the energy to 1e-8;
- the scene renders through `Report`;
- the equations of the systems with a closed form in Lee et al. (eqs. 5.38, 5.43, 5.71,
  6.16, 6.20, 7.39–7.40) agree with it at random admissible states;
- the Lagrangian and the virtual work agree (to 1e-8) with an independent model of the
  system's bodies, written in NumPy from its physical description (positions of the point
  masses and centres of mass, inertias about the centres of mass, heights, springs and
  the points where the forces act), at random admissible states; parametric systems at
  sizes 1, 2 and 3 (`tests/test_systems_physics.py`);
- the default parameters are physical: `e₃ = (0, 0, 1)`, positive masses and lengths, and
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
the derivation time of every system in its summary).

`System` is immutable: `params`, `symbols` and `options` are read-only mappings and every
array in them is a read-only copy owned by the system. To change a parameter, copy the
mapping: `p = dict(sys.params); p["m"] = 2.0`. The sizes `n_links` and `n_quads` run from
1 to 20 (`geomech.systems.MAX_SIZE`); other values raise `ValueError`.

## Not supported

Systems that need a configuration variable on S¹ (planar pendulums written with an angle,
the Furuta pendulum), nonholonomic constraints (rolling), closed kinematic chains,
contact, or flexible bodies other than chains of rigid links are not in the library:
the engine derives systems on R, R³, S², SO(3) and their products only.
