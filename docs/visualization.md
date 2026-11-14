# Visualization and reports

`geomech.viz` draws a simulated system in 3-D and writes one self-contained HTML report
with the animation, the derived equations, the standard form and the verification
numbers. It opens in a browser, in a notebook (`report` as the last expression of a
cell) or as a Claude artifact.

```python
from geomech.core.base.expressions import ZeroVector
from geomech.viz import Box, PointMass, Report, RigidBody, Rod, Scene, Trace

scene = Scene(layout, params=p)                       # the layout of the simulated module
scene.add(RigidBody(R1, origin=ZeroVector, shape=Box(0.12, 0.05, 1.0), offset=[0, 0, -0.5]))
scene.add(RigidBody(R2, origin=R1 * l1, shape=Box(0.12, 0.05, 0.8), offset=[0, 0, -0.4]))
scene.add(Trace(R1 * l1 + R2 * (rho2 * 2)))           # path of the tip of link 2

report = Report("Double rigid pendulum", system=(L, dW, variables, inputs), eom=eom,
                standard_form=sf, scene=scene, trajectory=traj, params=p)
report.write("double_rigid_pendulum.html")
```

`import geomech` does not import `geomech.viz`, and `geomech.viz` needs NumPy only (no
scipy). The page needs a browser with WebGL and network access to two CDNs.

## The scene

Geometry is declared once with geomech expressions and evaluated with
`geomech.codegen.lambdify` on the `StateLayout` of the module that was simulated, the
same code path as the generated equations: the visualizer has no kinematics of its own.

| Element | Call | Drawn as |
|---|---|---|
| rigid body | `RigidBody(R, origin, shape, offset=None, axes=True, color=None, name=None)` | the shape, centred at `origin + R offset` and rotated by `R`, with its body frame |
| point mass | `PointMass(position, radius=0.05)` | a sphere |
| cable, massless link | `Rod(a, b, radius=0.008)` | a cylinder from `a` to `b` |
| path of a point | `Trace(point)` | a line up to the current time |

- `R` is an `SO3` configuration variable (or a matrix expression of the configuration
  variables and parameters that is a rotation; `frames` raises `ValueError` otherwise).
  `origin`, `position`, `a`, `b` and `point` are vector expressions, or `None` / three
  numbers for a fixed point. `offset` is three numbers or a constant vector expression
  (a parameter such as `\rho_1`), in the body frame.
- Shapes, in the scene's units: `Box(lx, ly, lz)` with edges along the body axes,
  `Cylinder(radius, length, axis="z")`, `Sphere(radius)`. No meshes.
- Expressions may use the **configuration variables and parameters** of the layout only.
  A velocity, an input or a symbol the layout does not know raises `ValueError` naming
  it when the element is added, e.g. `PointMass 'R1*\Omega_{R1}': the position uses the
  velocity '\Omega_{R1}'; scene geometry may use configuration variables and parameters
  only`.
- Several `RigidBody` entries may share one `R`: a quadrotor is a box, two arms and four
  rotor cylinders with the same attitude and different offsets
  (`examples/report_quadrotor_load.py`); `axes=False` on all but one.
- **S² points have no attitude.** A pendulum on S² is `Rod(origin, origin + l q)` plus a
  `PointMass`. A rigid body hanging from an S² point cannot be drawn: the configuration
  does not say how it is turned about the cable; model it on SO(3) instead.
- Colours: `"link1"` … `"link4"`, `"load"`, `"cable"`, `"rotor"`, `"trace"` (theme
  colours, different in light and dark mode) or `"#rrggbb"`. Bodies cycle through the
  link colours by default.

`scene.frames(traj)` evaluates the scene at every stored state (or `indices=` of them)
and returns a `SceneFrames` with `t` and, per body, `"position"` (the origin, `(N, 3)`),
`"rotation"` (`(N, 3, 3)`) and `"center"` (`origin + R offset`); per point mass its
position; per rod its end points `(N, 2, 3)`; per trace the point. A bare configuration
slot (`R1`, `x_L`) is read from the state directly (vectorised); every other expression
is a generated function called once per state. The test suite compares the frames with
an independent NumPy evaluation of the same geometry to 1e-12.

### Seeing the orientation

A link drawn as a line hides its rotation about its own axis. The report therefore
draws each rigid body as a solid with lighting and shadows, its **+e₁ face tinted red**
and its +e₂ face slightly green, dark edges, and its body frame as arrows (e₁ red, e₂
green, e₃ blue) that stick out of the shape on every side. Give a link a non-square
cross-section (`Box(0.12, 0.05, 1.0)` for a link along e₃) and a spin about the link,
i.e. a component of Ω along e₃, is unmistakable: the double pendulum example starts with
`Ω₁ = (0.4, −0.6, 4.0)` and turns its links several times about their axes while they
swing.

## The report

`Report(title, system=(L, dW, variables, inputs), eom=None, standard_form=None,
scene=None, trajectory=None, params=None, checks=True, *, units=None, source=None,
derive_time=None, fps=60, description="")`

- `eom` and `standard_form` are derived (and the derivation timed) when not given; pass
  `derive_time` when you derived them yourself.
- `checks=True` runs `check_eom(L, dW, variables, eom, inputs=inputs)`; an error above
  `1e-5` is shown in red as a failure. Never present equations that fail it.
- `units={"m": "kg"}` fills the units column of the parameter table; `source` is shown
  verbatim at the end (the examples pass `inspect.getsource(model)`).
- `fps` is the number of frames per second of simulated time kept in the animation; the
  states are picked from the trajectory (the first stored state at or after each
  `k / fps`, and the last), never interpolated. A 10 s run at `dt = 1e-3` stores 10 001
  states; the page keeps 601.
- `report.write(path)` writes the page; `report.html()` returns it; in Jupyter
  `_repr_html_` shows it in an `<iframe srcdoc>` so that its MathJax does not clash with
  the notebook's.

Sections, in order:

1. **Summary**: derivation time, `check_eom` error, energy drift
   `max |E − E₀| / |E₀|` of `trajectory.energy()` (for a forced run this is the energy
   change, not an error), constraint drift, integrator and step, simulated time.
2. **Animation**: play/pause, a time scrubber, speed (0.1× to 4×), reset camera, toggles
   for body axes and traces, a theme switch. Drag to orbit, scroll to zoom, right-drag to
   pan.
3. **Lagrangian** and virtual work, one term per line.
4. **Equations of motion**, one per configuration variable, one term per line.
5. **Standard form**: the block overview `[M][a] + [f] + [G][u] = 0`, then each inertia
   block `M_ij`, each bias `f_i` (one term per line) and each input block `G_ik`.
6. **Parameters**: name, value, units.
7. **System definition**: the source passed in.

### Rendering

- **External resources**: three.js r147
  (`https://cdn.jsdelivr.net/npm/three@0.147.0/build/three.min.js` and
  `.../examples/js/controls/OrbitControls.js`) and MathJax 3.2.2
  (`https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-svg-full.min.js`).
  Nothing else is fetched; the trajectory is inlined as JSON (rounded to 5 decimals).
  Without network the equations stay readable as TeX source and the viewer says that
  three.js could not be loaded.
- **Frames**: the inertial frame has e₃ up and gravity along −e₃. three.js is y-up, so
  every inertial object is placed in one group rotated by −90° about x (inertial
  `(x, y, z)` is drawn at three.js `(x, z, −y)`); nothing else is converted. A ground
  grid lies in the inertial xy plane (at z = 0 when the scene is just above it,
  otherwise just below the lowest point), with the inertial axes at the origin, a
  perspective camera with orbit controls, and lights fixed in the world.
- **Theme**: colours are CSS variables, switched by `prefers-color-scheme` and by the
  `data-theme="light"`/`"dark"` attribute of the root element; the canvas background,
  grid and body colours follow them.
- **Layout**: responsive down to 400 px; the 3-D view stays square (at most 720 px);
  each equation scrolls horizontally in its own box.
- **Size**: about 175 kB for the 10 s double pendulum, 260 kB for the quadrotor with a
  load (8 bodies): far below the 2 MB bound the tests enforce.

### Readable symbols

The engine prints a leaf by its name, so `J1 \dot{\Omega_{R1}}` would be shown. The report
builds a map from leaf names to LaTeX with `latex_names(L, dW, eom, sf, …)` and prints
with `to_latex(expr, names=...)`, which applies it to leaf names only, never to an
assembled formula (a pattern over a whole formula can match across symbol boundaries
and corrupt it).

| Leaf name | Printed |
|---|---|
| letters (or one `\command`) followed by digits: `J1`, `e3`, `m12`, `R1` | `J_1`, `e_3`, `m_{12}`, `R_1` |
| tangent or variation vector of a point named letters + digits: `\Omega_{R1}`, `\eta_{R2}` | `\Omega_1`, `\eta_2` |
| anything else: `m_L`, `\rho_1`, `q`, `\omega_{q}`, numbers | unchanged |

The map is one-to-one: when two leaves would print alike (`J1` next to a leaf named
`J_1`, or `\Omega_{R1}` and `\Omega_{Q1}`), neither is renamed. `pretty_name(name,
point=None)` applies the rule to one name.

## Examples

`examples/report_double_rigid_pendulum.py`, `examples/report_spherical_pendulum.py` and
`examples/report_quadrotor_load.py` (LQR about hover; needs scipy) each derive, verify,
simulate 10 s and write their report to the current directory (or to the path given as
the first argument). CI runs them and discards the output.

| Example | Generation time | Size |
|---|---|---|
| double rigid pendulum | 5.3 s (the simulation: 4.8 s) | 175 kB |
| spherical pendulum | 2.5 s | 85 kB |
| quadrotor with a load (LQR) | 5.7 s | 259 kB |
