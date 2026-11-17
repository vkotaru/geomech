# Visualization

`geomech.viz` draws a simulated system in 3-D and writes one self-contained HTML report
with the animation, the derived equations, the standard form and the verification
numbers. It opens in a browser, in a notebook (`report` as the last expression of a
cell) or as a Claude artifact.

```python
from geomech.core.base.expressions import ZeroVector
from geomech.viz import Box, PointMass, Report, RigidBody, Rod, Scene, Trace

scene = Scene(layout)                                 # the layout of the simulated module
scene.add(RigidBody(R1, origin=ZeroVector, shape=Box(0.12, 0.05, 1.0), offset=[0, 0, -0.5]))
scene.add(RigidBody(R2, origin=R1 * l1, shape=Box(0.12, 0.05, 0.8), offset=[0, 0, -0.4]))
scene.add(Trace(R1 * l1 + R2 * (rho2 * 2)))           # path of the tip of link 2

report = Report("Double rigid pendulum", system=(L, dW, variables, inputs), eom=eom,
                standard_form=sf, scene=scene, trajectory=traj)     # drawn with traj.p
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
| point mass | `PointMass(position, radius=None)` | a sphere |
| cable, massless link | `Rod(a, b, radius=None)` | a cylinder from `a` to `b` |
| path of a point | `Trace(point)` | a line up to the current time |

- `R` is an `SO3` configuration variable (or a matrix expression of the configuration
  variables and parameters that is a rotation; `frames` raises `ValueError` otherwise).
  `origin`, `position`, `a`, `b` and `point` are vector expressions, or `None` / three
  numbers for a fixed point. `offset` is three numbers or a constant vector expression
  (a parameter such as `\rho_1`), in the body frame.
- Shapes, in the scene's units: `Box(lx, ly, lz)` with edges along the body axes,
  `Cylinder(radius, length, axis="z")`, `Sphere(radius)`. No meshes.
- **Scale.** Nothing in the view has an absolute size: the camera frames everything that
  is drawn (shapes, body-axis arrows, point and rod radii) with no minimum, the arrow
  thickness is a fraction of the arrow length, positions keep six significant figures of
  the scene size, and `radius=None` (the default) makes a point mass 6 % and a rod 1 % of
  the scene's half-size (at most 35 % and 6 % of the largest body). A 1 mm and a 1 km
  copy of a scene give the same picture.
- Expressions may use the **configuration variables and parameters** of the layout only.
  A velocity, an input or a symbol the layout does not know raises `ValueError` naming
  it when the element is added, e.g. `PointMass 'R1*\Omega_{R1}': the position uses the
  velocity '\Omega_{R1}'; scene geometry may use configuration variables and parameters
  only`.
- Several `RigidBody` entries may share one `R`: a quadrotor is a box, two arms and four
  rotor cylinders with the same attitude and different offsets
  (`examples/report_quadrotor_load.py`); `axes=False` on all but one.
- **$S^2$ points have no attitude.** A pendulum on $S^2$ is `Rod(origin, origin + l q)` plus a
  `PointMass`. A rigid body hanging from an $S^2$ point cannot be drawn: the configuration
  does not say how it is turned about the cable; model it on $SO(3)$ instead.
- Colours: `"link1"` … `"link4"`, `"load"`, `"cable"`, `"rotor"`, `"trace"` (theme
  colours, different in light and dark mode) or `"#rrggbb"`. Bodies cycle through the
  link colours by default.

### Parameters: one source of truth

`Scene(layout)` takes no parameters. `scene.frames(traj)` of a `Trajectory`, and the
`Report`, evaluate the scene with `traj.p`, the parameters that were simulated, and the
report's parameter table shows the same `traj.p`; so the picture, the table and the
equations always describe the same system. `Scene(layout, params=p)` and
`Report(..., params=p)` are still accepted, but a value that differs from `traj.p`
raises `ValueError` naming the parameter and both values (for example a pendulum drawn
with `l = 2.5` next to a simulation with `l = 1`). `params` is needed only to evaluate
bare state arrays, `scene.frames(states)`.

### Evaluation

`scene.frames(traj)` evaluates the scene at every stored state (or `indices=` of them)
and returns a `SceneFrames` with `t` and, per body, `"position"` (the origin, `(N, 3)`),
`"rotation"` (`(N, 3, 3)`) and `"center"` (`origin + R offset`); per point mass its
position; per rod its end points `(N, 2, 3)`; per trace the point. A bare configuration
slot (`R1`, `x_L`) is read from the state directly (vectorised); every other expression
is a generated function called once per state. The test suite compares the frames with
an independent NumPy evaluation of the same geometry to 1e-12.

A state with a NaN or infinite entry (a diverged simulation) raises `ValueError` in
`frames` and in `Report`, naming the first such time and state entry, e.g.
`Scene.frames: the state at t = 0.05 (stored state 5) is not finite: R[1, 1] = inf; …`.
Shorten `t_span` to end before it, or reduce `dt`.

### Seeing the orientation

A link drawn as a line hides its rotation about its own axis. The report therefore
draws each rigid body as a solid with lighting and shadows, its **$+e_1$ face tinted red**
and its $+e_2$ face tinted green, dark edges, and its body frame as arrows ($e_1$ red, $e_2$
green, $e_3$ blue) that stick out of the shape on every side. A cylinder carries a tinted
60° sector of its side and end caps on body $e_1$ (on $e_2$ for a cylinder along $e_1$), and a
sphere a cap on $e_1$ and a smaller one on $e_2$, so the spin of a wheel, rotor or ball shows
with its axes hidden (`axes=False` or the "Body axes" toggle). Each tint is the first mix
of the body colour towards the axis colour that differs from the body colour by at least
$\Delta E = 25$ (CIE76); when the axis colour itself is too close (a red `load` body and the red
$e_1$), its lightness is shifted away from the body's instead, so the tint shows on every
palette and custom colour. Give a link a non-square
cross-section (`Box(0.12, 0.05, 1.0)` for a link along $e_3$) and a spin about the link,
i.e. a component of $\Omega$ along $e_3$, is unmistakable: the double pendulum example starts with
$\Omega_1 = (0.4, -0.6, 4.0)$ and turns its links several times about their axes while they
swing.

The body arrows start at the shape's centre. A link hinged at the origin and hanging
along $-e_3$ therefore points its $e_3$ arrow up through the pivot, along the inertial $z$ arrow
of the origin triad, so two blue arrows overlap there; orbit the camera, or switch off
"Body axes", to separate them (the corner gizmo always shows the inertial frame).

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
- With a trajectory, `params` defaults to (and must agree with) `trajectory.p`; see
  [Parameters](#parameters-one-source-of-truth).
- `fps` is the number of frames per second of simulated time kept in the animation; the
  states are picked from the trajectory (the first stored state at or after each
  `k / fps`, and the last), never interpolated. A 10 s run at `dt = 1e-3` stores 10 001
  states; the page keeps 601.
- `report.write(path)` writes the page; `report.html()` returns it; in Jupyter
  `_repr_html_` shows it in an `<iframe srcdoc>` so that its MathJax does not clash with
  the notebook's.
- A written report opened with `?view` in its URL (`report.html?view`) shows only the 3-D
  view and its playback controls: no title, summary, equations or tables. The model pages
  of these docs embed their reports this way, in an `<iframe>`.

Sections, in order:

1. **Summary**: derivation time, `check_eom` error, energy drift
   $\max \lvert E - E_0 \rvert / \lvert E_0 \rvert$ of `trajectory.energy()`, constraint drift, integrator and
   step, simulated time. When the run is forced (`trajectory.u`, the inputs at the stored
   states, is nonzero somewhere) the card reads **Energy change** with the note that it
   is the work of the inputs, not an integrator error (`report.forced()`). A trajectory
   without recorded inputs (`trajectory.u is None`, e.g. built by hand) also reads
   **Energy change**, with the note "inputs not recorded": "drift" is claimed only for a
   run whose inputs are known to be zero.
2. **Animation**: play/pause, a time scrubber, speed (0.1× to 4×), reset camera, toggles
   for body axes and traces, a theme switch. Drag to orbit, scroll to zoom, right-drag to
   pan. Playback shows every kept frame, the last one included, before it starts over; a
   one-state trajectory has Play disabled.
3. **Lagrangian** and virtual work, one term per line.
4. **Equations of motion**, one per configuration variable, one term per line.
5. **Standard form**: the block overview $[M][a] + [f] + [G][u] = 0$, then each inertia
   block $M_{ij}$, each bias $f_i$ (one term per line) and each input block $G_{ik}$.
6. **Parameters**: name, value, units.
7. **System definition**: the source passed in.

### Rendering

- **External resources**: three.js r147
  (`https://cdn.jsdelivr.net/npm/three@0.147.0/build/three.min.js` and
  `.../examples/js/controls/OrbitControls.js`) and MathJax 3.2.2
  (`https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-svg-full.min.js`).
  Nothing else is fetched; the trajectory is inlined as JSON (positions to six
  significant figures of the scene size, rotation entries to 5 decimals).
  Without network the equations stay readable as TeX source and the viewer says that
  three.js could not be loaded.
- **Frames**: the inertial frame has $e_3$ up and gravity along $-e_3$. three.js is y-up, so
  every inertial object is placed in one group rotated by $-90^\circ$ about $x$ (inertial
  $(x, y, z)$ is drawn at three.js $(x, z, -y)$); nothing else is converted. A ground
  grid lies in the inertial $xy$ plane (at $z = 0$ when the scene is just above it,
  otherwise just below the lowest point), with the inertial axes at the origin, a
  perspective camera with orbit controls, and lights fixed in the world. A small
  inertial triad in the bottom-left corner, drawn in its own viewport and turned with the
  camera, is always visible; when the origin is more than three scene radii away, a
  second inertial triad also stands on the ground behind the scene.
- **Equations and user text**: MathJax typesets the equation containers only, never the
  title, description or units, and runs without the TeX packages that emit links, CSS or
  HTML attributes or load other packages (`html`, `action`, `bbox`, `unicode`,
  `require`, `autoload`, `setoptions`). Symbol names are user text: a name such as
  `\href{javascript:…}{x}` renders as an error, never as a link. (The `ui/safe`
  extension is not used because loading it would fetch a fourth URL.)
- **Theme**: colours are CSS variables, switched by `prefers-color-scheme` and by the
  `data-theme="light"`/`"dark"` attribute of the root element; the canvas background,
  grid and body colours follow them.
- **Layout**: responsive down to 400 px; the 3-D view stays square (at most 720 px);
  each equation scrolls horizontally in its own box, which shows a shadow on each side
  that has more to scroll to.
- **Framing**: the view is framed on the whole trajectory (every drawn element at every
  kept frame). An open-loop run that drifts far away therefore shows small bodies at
  $t = 0$; for the report, simulate a shorter horizon or close the loop with a controller.
- **Ready hook**: the root element carries `data-gm-viewer` (`"ready"` once the scene
  is built, `"no-three"` or `"no-webgl"` when it cannot be drawn, `"none"` without a
  scene), next to `data-gm-errors`, `data-gm-mathjax` and `data-gm-canvas-colors`; the
  browser tests wait for it. They are opt-in (`GEOMECH_BROWSER_TESTS=1`, see
  `tests/browser.py`) and skip without Chrome, the CDNs or WebGL.
- **Size**: about 180 kB for the 10 s double pendulum, 290 kB for the quadrotor with a
  load (8 bodies): far below the 2 MB bound the tests enforce.

### Readable symbols

The engine prints a leaf by its name, so `J1 \dot{\Omega_{R1}}` would be shown. The report
builds a map from leaf names to LaTeX with `latex_names(L, dW, eom, sf, …)` and prints
with `to_latex(expr, names=...)`, which applies it to leaf names only, never to an
assembled formula (a pattern over a whole formula can match across symbol boundaries
and corrupt it).

| Leaf name | Printed |
|---|---|
| letters, or one symbol command (a Greek letter, `\ell`, …), followed by digits: `J1`, `e3`, `m12`, `\theta2` | `J_1`, `e_3`, `m_{12}`, `\theta_2` |
| tangent or variation vector of a point named letters + digits: `\Omega_{R1}`, `\eta_{R2}` | `\Omega_1`, `\eta_2` |
| anything else: `m_L`, `\rho_1`, `q`, `\omega_{q}`, numbers, and a command that takes an argument (`\bar1`, `\hat1`, `\vec1`, `\mathbf1`: `\bar_1` would be invalid TeX) | unchanged |

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
| double rigid pendulum | 7.2 s (mostly the simulation) | 181 kB |
| spherical pendulum | 2.7 s | 91 kB |
| quadrotor with a load (LQR) | 6.3 s | 287 kB |
