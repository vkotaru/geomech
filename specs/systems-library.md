# Spec: verified systems library (`geomech.systems`)

Roadmap feature 7. Branch `feature/systems-library` from `main`.

## Goal

Turn the example systems into reproducible, citable objects: one function per system that
returns everything needed to derive, verify, linearize, generate, simulate and visualize
it, with the source it comes from, and a stored reference output that the test suite
re-checks. A user (or an agent following the skill) can then start from a known-good
system instead of retyping a Lagrangian.

## API

```python
from geomech.systems import get, list_systems, System

sys = get("double_rigid_pendulum")          # or geomech.systems.double_rigid_pendulum()
sys.L, sys.dW, sys.variables, sys.inputs     # what compute_eom / to_standard_form need
sys.params                                   # default numeric parameter values (dict, str -> float/array)
sys.x0                                       # a consistent initial state (packed with StateLayout)
sys.reference                                # Reference(source=..., equation=..., note=...)
sys.scene(layout)                            # optional geomech.viz Scene of the bodies
list_systems()                               # names, manifolds, dof, reference, in a stable order
```

- `System` is a frozen dataclass. Builder functions take only keyword arguments with
  defaults (for example `n_links=3`, `n_quads=3`) and are cheap (no derivation at build).
- Every symbol name is unique within a system and follows the conventions in
  `docs/conventions.md` (constants flagged `Constant`, inputs not constant, parameters with
  `params=` only when prescribed in time).
- `reference.source` is a citation key that exists in a new `docs/references.bib` (or a
  short text reference), and `reference.equation` names the equation or section when the
  source gives the equations (for example Lee et al. eq. 5.14).

## Systems

Only systems the library derives today, on R, R³, S², SO(3) and their products. Systems
needing S¹ variables or nonholonomic constraints are out of scope (listed in the docs as
not supported).

| Name | Configuration | Source |
|---|---|---|
| `point_mass` | R³ | Lee et al. ch. 3 |
| `spring_mass` | R | elementary |
| `spherical_pendulum` | S² | Lee et al. §5.3, eq. 5.14 |
| `particle_on_rotating_sphere` | S² | Lee et al. §5.6.2, eq. 5.43 |
| `double_spherical_pendulum` | S² × S² | Lee et al. §5.6 |
| `n_link_spherical_chain(n_links)` | (S²)ⁿ | Lee et al. ch. 5 |
| `elastic_spherical_pendulum` | R × S² | |
| `cart_spherical_pendulum` (3-D cart-pole) | R × S² | |
| `free_rigid_body` | SO(3) | Lee et al. §6, eq. 6.16 |
| `heavy_top` | SO(3) | Lee et al. §6 |
| `rigid_pendulum` | SO(3) | Lee et al. §6, eq. 6.8 |
| `double_rigid_pendulum` | SO(3) × SO(3) | |
| `two_rigid_bodies_spherical_joint` | R³ × SO(3) × SO(3) | |
| `rigid_body` (6-DOF) | R³ × SO(3) | Lee et al. ch. 6 |
| `quadrotor` | R³ × SO(3) | Lee, Leok, McClamroch 2010 (geometric control) |
| `quadrotor_offset_com` | R³ × SO(3) | |
| `quadrotor_suspended_load` | R³ × SO(3) × S² | Sreenath et al. 2013 |
| `quadrotor_elastic_load` | R³ × SO(3) × R × S² | |
| `quadrotor_flexible_cable(n_links)` | R³ × SO(3) × (S²)ⁿ | |
| `multi_quadrotor_point_mass(n_quads)` | R³ × (SO(3) × S²)ⁿ | |
| `multi_quadrotor_rigid_load(n_quads)` | R³ × SO(3) × (SO(3) × S²)ⁿ | |

Sources are filled in from the literature only when the library's equations can be
compared with them; otherwise the entry says what the system is and that it is verified by
`check_eom` alone. Do not invent citations: if a source cannot be confirmed, leave it out.

## Verification (definition of done)

1. Every system: `compute_eom` succeeds; `check_eom ≤ 1e-5` at three seeds;
   `to_standard_form` with the system's inputs succeeds; `StateLayout` + `generate_module`
   succeed; `simulate` for 1 s from `x0` runs and, for unforced runs, conserves energy to
   1e-8; `sys.scene` (when present) renders through `Report` without error.
2. Where the source gives the equation (Lee et al. book equations), a test compares the
   library's equation with the source's closed form numerically at random admissible states
   (tangent plane for S²), with the hand derivation in the docstring.
3. Reference outputs: `tests/golden/<name>.json` stores, per equation, the printed equation
   string and its value at fixed seeded states (and the standard-form blocks). A test
   re-derives every system and fails on any change of value; a change of printed string
   alone fails with a message that tells the author to regenerate the files with
   `python scripts/regenerate_golden.py` after checking the numeric equality.
4. Parametric systems are tested at two sizes (for example n = 2 and 3) and derive in under
   10 s each at the larger size; report the derivation time per system.
5. `list_systems()` output is stable and documented; a docs page "Systems library"
   (Guides) is generated from it (table of name, configuration, dof, source), and each
   system's docstring appears in the API reference.
6. The existing test builders (`tests/vbl_systems.py`, `tests/phase1_systems.py`, …) that
   duplicate a library system are switched to import it, without changing any test
   assertion; builders that exist only for adversarial tests stay in `tests/`.
7. The skill (`skills/geomech/SKILL.md` and its byte-identical copy under
   `.claude/skills/geomech/`) mentions `geomech.systems.get(...)` as the starting point
   for known systems.
8. ruff, mypy, the full suite, and the docs build with `-W` pass.

## Non-goals

S¹ variables, constraints (rolling, closed chains), contact, flexible bodies beyond link
chains, generating text for any document outside this repository.
