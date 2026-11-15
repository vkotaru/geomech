---
name: geomech
description: Derive, verify, linearize, generate code for, simulate and visualize the equations of motion of mechanical systems on R, R³, S¹, S², SO(3) and their products (rigid bodies, planar and spherical pendulums, multi-link chains, reaction wheels, aerial vehicles with cable-suspended loads) with the geomech Python package. Use when asked to "derive the equations of motion of X", to check or linearize such equations, to generate NumPy/CasADi code for them, to simulate them, or to show them as a 3-D report.
---

# geomech: equations of motion on manifolds

Use this procedure whenever a task asks for the dynamics of a mechanical system whose
configuration lives on R (scalars), R³ (positions), S¹ (angles without wrap-around:
planar pendulums, wheels, joints about a fixed axis), S² (unit vectors: pendulums, cable
directions) or SO(3) (attitudes), or a product of these. SE(3) is R³ × SO(3).
Do not derive by hand and do not use Euler angles: geomech derives coordinate-free
equations (Lee et al., *Global Formulations of Lagrangian and Hamiltonian Dynamics on
Manifolds*).

## Setup

The skill needs the `geomech` package (Python 3.10 or later; NumPy is its only required
dependency). Check with `python -c "import geomech"`; if it fails, install it into the
project's environment:

```bash
pip install geomech              # derivation, checks, linearization, code generation
pip install "geomech[sim]"       # adds SciPy, needed by lqr and method="solve_ivp"
pip install "geomech[casadi]"    # adds CasADi for the CasADi code target
```

Run every snippet below with that environment's Python, and state the installed version
(`geomech.__version__`) with the result.

## Known systems: start from the library

If the system is one of `geomech.systems.list_systems()` (point mass, spring, planar and
double planar pendulums on S¹, spherical,
double spherical and n-link spherical pendulums, particle on a rotating sphere, elastic
pendulum, 3-D cart-pole, free rigid body, heavy top, rigid and double rigid pendulums,
n-link rigid chain,
two bodies with a spherical joint, 6-DOF rigid body, spacecraft with reaction wheels
(wheel angles on R³ or on (S¹)³),
3-D compass-gait biped (swing phase), quadrotor, quadrotor with an offset
centre of mass, with a suspended, elastic or flexible-cable load, several quadrotors
carrying a point mass or a rigid body), start from `sys = geomech.systems.get(name)`
(`get("multi_quadrotor_point_mass", n_quads=4)` for a size) instead of retyping the
Lagrangian: `sys.L`, `sys.dW`, `sys.variables`, `sys.inputs` go into the steps below,
`sys.params` and `sys.x0` into `simulate`, `sys.layout()` is the `StateLayout`, and
`sys.scene(layout)` the `Scene`. Still run step 4 (`check_eom`) and report it, and cite
`sys.reference`. Write the model by hand only for a system that is not in the list.

## Models of n bodies

When the user asks for a chain, cable or team of *n* bodies with n left general, use
`geomech.indexed` (`Family`, `Sum`, `IndexedVariables`, `compute_indexed_eom`; the five
models in `geomech.indexed.systems`) instead of fixing n. Verify the generic equation
with `eom.instantiate(k)` and `check_eom` on `model.unroll(k)` for several k, and say
which k were checked. See `docs/indexed.md`.

## Procedure

Follow the steps in order. Each call is exact; `reference.md` lists the API and the
error messages.

1. **Define symbols and variables.** Constants need `attr=["Constant"]` (case-sensitive);
   inertias also `"SymmetricMatrix"`. Configuration variables: `Scalar("s")`,
   `Vector("x")`, `S1("q", axis=e)`, `S2("q")`, `SO3("R")`. Velocities: `x.t_diff()`,
   `q.get_tangent_vector()` (ω, `q̇ = ω × q` on S²; a **scalar** ω with `q̇ = ω e × q` on
   S¹), `R.get_tangent_vector()` (Ω, body frame). Variations: `x.delta()`,
   `q.get_variation_vector()` (ξ on S²; the scalar γ on S¹, `δq = γ e × q`),
   `R.get_variation_vector()` (η). List every configuration variable in
   `SystemVariables(scalars=[…], vectors=[…], matrices=[…])` (S¹ and S² points go in
   `vectors`, SO(3) in `matrices`).
   **S¹** is the unit circle in the plane normal to a constant unit axis `e`, which must
   be built as `e = Vector("e", attr=["Constant", "UnitNorm"])` (an angle θ about `e`;
   the axis is a parameter, `p["e"]`, of the generated code). Several S¹ points may share
   one axis (a planar chain) or have their own (three reaction wheels on `b1`, `b2`, `b3`:
   `S1("q_1", axis=b1)`, …). Each S¹ point gives one scalar equation. A moment `τ` about
   the axis acting on an S¹ point enters as `γ * τ`.
2. **Write L and δW.** `L = T − V` with `Dot`, `Cross`, `*`; a half is
   `Scalar("0.5", value=0.5, attr=["Constant"])`. Forces enter as `Dot(variation, force)`.
   A force at a point enters as `Dot(δ(point), force)`, where `δ(point)` is the
   variation of the point's position expression term by term: `x` → `x.delta()`, `q` →
   `Cross(ξ, q)` (`δq = ξ × q`), `R` → `R * Cross(η, ·)` (`δR = R η̂`). A force `F` at
   `x + l q` gives `Dot(x.delta() + Cross(xi, q) * l, F)`; at a body point `x + R ρ` it
   gives `Dot(x.delta() + R * Cross(eta, rho), F)`. A moment `M` on `R` is
   `Dot(eta, M)`.
3. **Derive:** `eom = compute_eom(L, dW, variables)`.
4. **Verify (mandatory):** `err = check_eom(L, dW, variables, eom, inputs=inputs)`.
   Report the number. A correct equation gives about 1e-7; never present an equation
   whose error is above 1e-5. If it fails, report the failure and the error text; do not
   edit the equation by hand.
5. **Standard form:** `sf = to_standard_form(eom, variables, inputs)` gives
   `M a + f + G u = 0` per equation.
6. Optionally **linearize:** `lin = linearize(eom, variables, inputs)`
   (`M δẋ = A δx + B δu`).
7. **Generate code:** `layout = StateLayout(variables, inputs, collect_params([L, dW],
   variables, inputs))`, then `plant = load_module(module_source(sf, layout,
   lagrangian=L), "name")` (or `generate_module(sf, layout, "file.py", lagrangian=L)`).
8. **Simulate:** `traj = simulate(plant, layout.pack(...), t_final, p)`; check
   `traj.energy()` drift and `traj.max_drift()`. The `pack` keys are the configuration
   names and `x_dot`, `omega_<name>`, `Omega_<name>` (`omega_q1` for `S2("q1")` or
   `S1("q1", …)`, a scalar for S¹, `Omega_R2` for `SO3("R2")`); `layout.slot_names()`
   lists them. An S¹ point must start on its circle (`|q| = 1`, `q·e = 0`), and `p` must
   hold its axis. Energy is conserved
   only without inputs: with nonzero inputs (`u=`, recorded in `traj.u`) `E(t) − E(0)`
   is the work the inputs did, not an integrator error; call it "energy change", never
   "drift". For a forced run, also simulate once from the same `x0` over the same
   horizon with `u=None` and report that run's energy drift as the integrator check.
9. **Report:** build a `Scene(layout)` (rigid bodies as boxes with a non-square
   cross-section, point masses, rods, traces) and write `Report(...).write("x.html")`.
   Do not pass parameters to the `Scene`: the report draws the scene and fills the
   parameter table with `traj.p`, the parameters that were simulated, so the picture
   cannot show a different system. A `params=` (on `Scene` or `Report`) that disagrees
   with `traj.p` raises `ValueError` naming the parameter and both values. A forced run
   is labelled "Energy change" in the summary automatically. The view is framed on the
   whole trajectory: if an open-loop run drifts far away, the bodies look small at
   t = 0, so for the report simulate a shorter horizon or close the loop with a
   controller (`lqr`).

## Conventions to state with every result

- Each equation is `E = 0` with **E = −(Euler–Lagrange) + forcing** (a point mass gives
  `−m ẍ + F = 0`).
- In the standard form **M = −(physical inertia)** (`M[ẍ] = −m I`, `M[Ω̇] = −J`).
- **S² equations hold on the tangent plane** of their point: only `(I − qqᵀ) E = 0` is
  an equation; the S² mass block can be singular, which is expected.
- **An S¹ equation is a scalar**, `E = −d/dt(∂L/∂ω) + (e × q)·∂L/∂q + forcing` (Lee et al.
  eq. 4.12); the book's `S x` is printed `Cross(e,x)`, so its `pᵀ S r` is
  `Dot(p,Cross(e,r))`.
- `print_eom(eom)` prints each equation as LaTeX source (`\int \xi_{q} \cdot (… -m l …)
  dt`), like `to_latex(expr)`; `str(expr)` (and `print(expr)`) is the plain form, where
  products are written with `*` (`m*l*l*(-1)`). Say which one you show.
- Ω is the body-frame angular velocity; `Ṙ = R Ω̂`. Gravity is along −e₃ with
  `V = m g e₃·x`.

## Failure modes: report the error text, do not work around it

| Error text (excerpt) | Cause | What to do |
|---|---|---|
| `compute_eom: δS still contains the variation symbols '\delta{y}' of the symbols 'y', which are neither entries of SystemVariables nor declared in params` | stray variation: `y` (or a constant such as `m` without `Constant`) is not declared | list it in `SystemVariables`, declare a prescribed time-varying symbol with `params=[y]`, or mark it `attr=["Constant"]` |
| `attr: unknown attribute 'constant' for symbol 'm'` | attribute typo (case-sensitive) | `Constant`, `SymmetricMatrix`, … as listed in the message |
| `to_standard_form: input 'g' has the name of the constant 'g' in the equations` | input/constant clash | rename one; an input is a non-constant leaf |
| an `AlgebraicError` naming a calculus operator (`\frac{d}{dt}(\frac{d}{dt}(x))` in `L`) | acceleration-dependent Lagrangian | geomech supports first-order Lagrangians only; rewrite `L` in velocities |
| `NonLinearError: … nonlinear — …`, `AlgebraicError: … no left matrix coefficient`, `collect: the target … occurs inside the matrix factor` | unsupported extraction shape (a variation entering nonlinearly, `M·Hat(v)`, a velocity inside a matrix factor) | restructure the term as the message suggests (e.g. `(a·Ω)(c·b)`); otherwise report it |
| `check_eom` error above 1e-5 | a wrong equation | stop and report; `docs/troubleshooting.md` has the localisation recipe |

Quote the error text in the answer; do not silence or work around an error.

## Presentation

Always present the report (`Report(...).write(...)`, or `report._repr_html_()` in a
notebook) or, in a terminal, the equations (`print_eom(eom)` for LaTeX source, or
`print(eq)` for each `(variation, eq)` of `eom.values()` in the `*` form) with the
`check_eom` number, the energy drift (or, for a forced run, the energy change, that it
is the work of the inputs, and the drift of the same run with `u=None`) and the
conventions above. Never draw a diagram by hand.

## Worked example: spherical pendulum (Lee et al. Eq. 5.14)

```python
import numpy as np
from geomech import S2, Dot, Scalar, SystemVariables, Vector, check_eom, compute_eom
from geomech import getScalars, print_eom, to_standard_form
from geomech.codegen import StateLayout, collect_params, load_module, module_source
from geomech.core.base.expressions import ZeroVector
from geomech.simulation import simulate
from geomech.viz import PointMass, Report, Rod, Scene, Trace

half = Scalar("0.5", value=0.5, attr=["Constant"])
m, g, l = getScalars("m g l", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
q, f = S2("q"), Vector("f")
omega, xi = q.get_tangent_vector(), q.get_variation_vector()
L = half * m * l * l * Dot(omega, omega) - m * g * l * Dot(q, e3)
dW = Dot(xi, f)
variables, inputs = SystemVariables(vectors=[q]), [f]

eom = compute_eom(L, dW, variables)
err = check_eom(L, dW, variables, eom, inputs=inputs)
assert err < 1e-5, err
print_eom(eom)
sf = to_standard_form(eom, variables, inputs)

layout = StateLayout(variables, inputs, collect_params([L, dW], variables, inputs))
plant = load_module(module_source(sf, layout, lagrangian=L), "spherical_pendulum")
p = {"e3": np.array([0.0, 0.0, 1.0]), "g": 9.81, "l": 1.0, "m": 1.0}
q0 = np.array([np.sin(1.0), 0.0, -np.cos(1.0)])
w0 = np.array([0.0, 1.5, 0.0]) - (np.array([0.0, 1.5, 0.0]) @ q0) * q0
traj = simulate(plant, layout.pack(q=q0, omega_q=w0), 2.0, p)
E = traj.energy()
print(f"check_eom {err:.1e}, energy drift {np.max(np.abs(E - E[0])) / abs(E[0]):.1e}")

bob = q * l
scene = Scene(layout)                  # drawn with traj.p, the simulated parameters
scene.add(Rod(ZeroVector, bob))
scene.add(PointMass(bob, radius=0.08))
scene.add(Trace(bob))
Report("Spherical pendulum", system=(L, dW, variables, inputs), eom=eom, standard_form=sf,
       scene=scene, trajectory=traj).write("spherical_pendulum.html")
```

## Worked example: double rigid pendulum (SO(3) × SO(3))

Give each link a box with a non-square cross-section and an initial Ω with a component
along the link so the spin about the link's own axis shows.

```python
import numpy as np
from geomech import SO3, Dot, Matrix, Scalar, SystemVariables, Vector, check_eom
from geomech import compute_eom, getScalars, to_standard_form
from geomech.codegen import StateLayout, collect_params, load_module, module_source
from geomech.core.base.expressions import ZeroVector
from geomech.numeric import expm_so3
from geomech.simulation import simulate
from geomech.viz import Box, Report, RigidBody, Scene

half = Scalar("0.5", value=0.5, attr=["Constant"])
J1 = Matrix("J1", attr=["Constant", "SymmetricMatrix"])
J2 = Matrix("J2", attr=["Constant", "SymmetricMatrix"])
rho1, rho2 = Vector("\\rho_1", attr=["Constant"]), Vector("\\rho_2", attr=["Constant"])
l1 = Vector("l_1", attr=["Constant"])
m1, m2, g = getScalars("m1 m2 g", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
M1, M2 = Vector("M1"), Vector("M2")
R1, R2 = SO3("R1"), SO3("R2")
Om1, Om2 = R1.get_tangent_vector(), R2.get_tangent_vector()
x1, x2 = R1 * rho1, R1 * l1 + R2 * rho2
v1, v2 = x1.t_diff(), x2.t_diff()
L = (half * Dot(Om1, J1 * Om1) + half * Dot(Om2, J2 * Om2) + half * m1 * Dot(v1, v1)
     + half * m2 * Dot(v2, v2) - m1 * g * Dot(x1, e3) - m2 * g * Dot(x2, e3))
dW = Dot(R1.get_variation_vector(), M1) + Dot(R2.get_variation_vector(), M2)
variables, inputs = SystemVariables(matrices=[R1, R2]), [M1, M2]

eom = compute_eom(L, dW, variables)
err = check_eom(L, dW, variables, eom, inputs=inputs)
assert err < 1e-5, err
sf = to_standard_form(eom, variables, inputs)
layout = StateLayout(variables, inputs, collect_params([L, dW], variables, inputs))
plant = load_module(module_source(sf, layout, lagrangian=L), "double_pendulum")

p = {"J1": np.diag([0.084, 0.085, 0.0014]), "J2": np.diag([0.043, 0.044, 0.0011]),
     "\\rho_1": np.array([0, 0, -0.5]), "\\rho_2": np.array([0, 0, -0.4]),
     "l_1": np.array([0, 0, -1.0]), "m1": 1.0, "m2": 0.8, "g": 9.81,
     "e3": np.array([0, 0, 1.0])}
x0 = layout.pack(R1=expm_so3([0.7, 0.3, 0.0]), R2=expm_so3([1.1, -0.4, 0.3]),
                 Omega_R1=[0.4, -0.6, 4.0], Omega_R2=[-0.8, 0.5, -5.0])
traj = simulate(plant, x0, 1.0, p)

scene = Scene(layout)
scene.add(RigidBody(R1, origin=ZeroVector, shape=Box(0.12, 0.05, 1.0), offset=[0, 0, -0.5]))
scene.add(RigidBody(R2, origin=R1 * l1, shape=Box(0.12, 0.05, 0.8), offset=[0, 0, -0.4]))
report = Report("Double rigid pendulum", system=(L, dW, variables, inputs), eom=eom,
                standard_form=sf, scene=scene, trajectory=traj)
report.write("double_pendulum.html")
print(f"check_eom {report.check_error:.1e}, energy drift {report.energy_drift():.1e}")
```

## Worked example: quadrotor with a cable-suspended load (R³ × S² × SO(3))

```python
import numpy as np
from geomech import S2, SO3, Cross, Dot, Matrix, Scalar, SystemVariables, Vector
from geomech import check_eom, compute_eom, getScalars, to_standard_form
from geomech.codegen import StateLayout, collect_params, load_module, module_source
from geomech.simulation import simulate
from geomech.viz import Box, Cylinder, PointMass, Report, RigidBody, Rod, Scene

half = Scalar("0.5", value=0.5, attr=["Constant"])
mQ, mL, g, l = getScalars("m_Q m_L g l", attr=["Constant"])
J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
e3 = Vector("e3", attr=["Constant"])
xL, q, R = Vector("x_L"), S2("q"), SO3("R")
Om, eta, xi = R.get_tangent_vector(), R.get_variation_vector(), q.get_variation_vector()
f, M = Scalar("f"), Vector("M")
xQ = xL - q * l                       # q points from the quadrotor to the load
vQ, vL = xQ.t_diff(), xL.t_diff()
L = (half * mQ * Dot(vQ, vQ) + half * mL * Dot(vL, vL) + half * Dot(Om, J * Om)
     - mQ * g * Dot(e3, xQ) - mL * g * Dot(e3, xL))
dW = Dot(xL.delta() - Cross(xi, q) * l, (R * e3) * f) + Dot(eta, M)
variables, inputs = SystemVariables(vectors=[xL, q], matrices=[R]), [f, M]

eom = compute_eom(L, dW, variables)
err = check_eom(L, dW, variables, eom, inputs=inputs)
assert err < 1e-5, err
sf = to_standard_form(eom, variables, inputs)
layout = StateLayout(variables, inputs, collect_params([L, dW], variables, inputs))
plant = load_module(module_source(sf, layout, lagrangian=L), "quadrotor_load")

p = {"J": np.diag([0.02, 0.025, 0.04]), "e3": np.array([0, 0, 1.0]), "g": 9.81,
     "l": 1.0, "m_L": 0.3, "m_Q": 1.5}
x0 = layout.pack(x_L=np.zeros(3), q=[0, 0, -1], R=np.eye(3), Omega_R=[0, 0, 1.0])
hover = np.r_[(p["m_Q"] + p["m_L"]) * p["g"], 0, 0, 0]       # u = [f, M]
traj = simulate(plant, x0, 1.0, p, u=hover)
free = simulate(plant, x0, 1.0, p)          # u=None: the integrator check of a forced run
E = free.energy()
print(f"unforced energy drift {np.max(np.abs(E - E[0])) / abs(E[0]):.1e}")

scene = Scene(layout)
scene.add(RigidBody(R, origin=xQ, shape=Box(0.16, 0.12, 0.06)))
for off in [(0.25, 0, 0.03), (0, 0.25, 0.03), (-0.25, 0, 0.03), (0, -0.25, 0.03)]:
    scene.add(RigidBody(R, xQ, Cylinder(0.09, 0.012), offset=off, axes=False, color="rotor"))
scene.add(Rod(xQ, xL))
scene.add(PointMass(xL, radius=0.07))
report = Report("Quadrotor with a load", system=(L, dW, variables, inputs), eom=eom,
                standard_form=sf, scene=scene, trajectory=traj)
report.write("quadrotor_load.html")
assert report.forced()     # hover thrust: the summary shows "Energy change", the inputs' work
```

See `reference.md` for the API, the state layout and the error messages, and
`examples/report_*.py` in the repository for complete scripts (the quadrotor one adds an
LQR controller about hover).
