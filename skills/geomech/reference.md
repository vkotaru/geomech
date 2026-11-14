# geomech reference for agents

The calls used by `SKILL.md`, with their arguments, and the facts an answer must state.
The full documentation is in `docs/` of the repository (`conventions.md`, `codegen.md`,
`simulation.md`, `visualization.md`, `troubleshooting.md`).

## Symbols

| Kind | Build | Notes |
|---|---|---|
| constant scalar | `m, g = getScalars("m g", attr=["Constant"])` | names are the identity: two leaves with one name are one symbol |
| number | `Scalar("0.5", value=0.5, attr=["Constant"])`, or `2 * x` | a named scalar with a value stays a symbol |
| constant vector | `Vector("e3", attr=["Constant"])` | `e3 = (0, 0, 1)` in `p`; gravity potential `m g e3·x` |
| inertia | `Matrix("J", attr=["Constant", "SymmetricMatrix"])` | |
| configuration | `Scalar("s")`, `Vector("x")`, `S2("q")`, `SO3("R")` | list in `SystemVariables` |
| velocity | `x.t_diff()`, `q.get_tangent_vector()` (ω), `R.get_tangent_vector()` (Ω) | `q̇ = ω × q`, `Ṙ = R Ω̂` |
| variation | `x.delta()`, `q.get_variation_vector()` (ξ), `R.get_variation_vector()` (η) | `δq = ξ × q`, `δR = R η̂` |
| input | `Scalar("f")`, `Vector("M")` (not constant) | passed as `inputs=[f, M]`, order matters |

Operators: `Dot(a, b)`, `Cross(a, b)`, `Hat(v)`, `Vee(M)`, `A * v`, `s * v`, `v.T()`.
`compute_eom(L, dW, variables, params=[l])` declares prescribed time-varying symbols.

## Pipeline

| Call | Returns |
|---|---|
| `compute_eom(L, dW, variables)` | `{str(variation): (variation, E)}`, `E = 0` |
| `check_eom(L, dW, variables, eom, inputs=inputs)` | worst relative error; ~1e-7 when correct |
| `to_standard_form(eom, variables, inputs)` | `{key: StandardFormEquation}` with `.M[accel]`, `.f`, `.G[input]` |
| `linearize(eom, variables, inputs)` | `LinearizedSystem` (`M δẋ = A δx + B δu`) |
| `print_eom(eom)`, `to_latex(expr, names=latex_names(eom))` | text and LaTeX |

## Code generation and simulation

```python
import numpy as np
from geomech import Dot, Scalar, SystemVariables, Vector, compute_eom, getScalars
from geomech import to_standard_form
from geomech.codegen import StateLayout, collect_params, load_module, module_source
from geomech.simulation import simulate

half = Scalar("0.5", value=0.5, attr=["Constant"])
m, k = getScalars("m k", attr=["Constant"])
x, F = Vector("x"), Vector("F")
L = half * m * Dot(x.t_diff(), x.t_diff()) - half * k * Dot(x, x)
dW = Dot(x.delta(), F)
V = SystemVariables(vectors=[x])
sf = to_standard_form(compute_eom(L, dW, V), V, [F])
layout = StateLayout(V, [F], collect_params([L, dW], V, [F]))   # params: k, m
plant = load_module(module_source(sf, layout, lagrangian=L), "oscillator")
p = {"m": 2.0, "k": 8.0}
traj = simulate(plant, layout.pack(x=[1.0, 0, 0], x_dot=[0, 1.0, 0]), 3.0, p)
assert np.allclose(traj["x"][-1][0], np.cos(2.0 * 3.0), atol=1e-6)   # ω = sqrt(k/m) = 2
print(traj.energy()[[0, -1]], traj.max_drift())
```

- `state` = configuration slots (variable order: scalars, vectors, matrices; `R`
  row-major), then velocity slots. `layout.pack(x=..., q=..., R=..., x_dot=...,
  omega_q=..., Omega_R=...)` builds it; missing velocities are zero;
  `layout.unpack(traj.state)` splits it.
- `p` holds exactly the parameters of the layout, keyed by printed name
  (`"\\rho_1"`, `"J1"`); a missing key raises `KeyError`, an extra one `ValueError`.
- `simulate(module, x0, t_final, p, u=None, method="rkmk4", dt=1e-3, output_every=1)`;
  `u` is `None`, a constant array (`[f, M]` stacked) or `u(t, state)`.
  `x0` must be on the manifold (`|q| = 1`, `q·ω = 0`, `R ∈ SO(3)`).
- `lqr(lin_module, x_eq, u_eq, p, Q, R, module=plant)` (needs scipy) returns a
  controller usable as `u`.

## Visualization (`geomech.viz`)

| Element | Call |
|---|---|
| rigid body | `RigidBody(R, origin, shape, offset=None, axes=True, color=None)`; shape centred at `origin + R offset`, rotated by `R` |
| shapes | `Box(lx, ly, lz)` (use `lx ≠ ly` for a link along e₃), `Cylinder(radius, length, axis="z")`, `Sphere(radius)` |
| point mass | `PointMass(position, radius=0.05)` |
| cable, massless link | `Rod(a, b, radius=0.008)` |
| path | `Trace(point)` |
| colours | `"link1"`…`"link4"`, `"load"`, `"cable"`, `"rotor"`, `"trace"` or `"#rrggbb"` |

Expressions may use configuration variables and parameters only (`ValueError` naming
the leaf otherwise). An S² point has no attitude: draw it as a rod plus a point mass.
Several `RigidBody` entries may share one `R` (a quadrotor body, arms and rotors).

```python
import numpy as np
from geomech import SO3, Dot, Matrix, Scalar, SystemVariables, Vector
from geomech.codegen import StateLayout
from geomech.core.base.expressions import ZeroVector
from geomech.viz import Box, RigidBody, Scene

R = SO3("R")
J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
half = Scalar("0.5", value=0.5, attr=["Constant"])
Om = R.get_tangent_vector()
layout = StateLayout(SystemVariables(matrices=[R]), [], [J])
scene = Scene(layout, params={"J": np.diag([1.0, 2.0, 3.0])})
scene.add(RigidBody(R, origin=ZeroVector, shape=Box(0.3, 0.1, 1.0)))
states = np.r_[np.eye(3).reshape(9), [0.0, 0.0, 1.0]][None, :]
frames = scene.frames(states)                  # or scene.frames(traj)
print(frames.bodies[0]["rotation"][0], frames.bodies[0]["center"][0])
try:
    scene.add(RigidBody(R, origin=Om, shape=Box(0.1, 0.1, 0.1)))
except ValueError as err:
    print(err)                                 # names the velocity '\Omega_{R}'
```

`Report(title, system=(L, dW, variables, inputs), eom=None, standard_form=None,
scene=None, trajectory=None, params=None, checks=True, units=None, source=None,
derive_time=None, fps=60)` derives what is missing, runs `check_eom` and writes one
HTML page with `report.write(path)`: summary (derivation time, `check_eom` error, energy
and constraint drift, integrator), the 3-D animation (three.js and MathJax from pinned
CDN URLs, the only external resources), L, the equations one term per line, the
standard form (block overview, each M, f and G block), the parameters and the source.

## Facts to state

- `E = −(Euler–Lagrange) + forcing`; the equation is `E = 0`.
- `M = −(physical inertia)` in the standard form; the linearization's M is positive.
- S² equations hold on the tangent plane only; the S² mass block may be singular.
- Products print with `*`; `str` is the identity of a symbol.
