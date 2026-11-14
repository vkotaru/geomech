# Spec: 3-D visualizer and report (`geomech.viz`) and the agent skill

Branch `feature/viz-skill`, on top of `feature/codegen-simulation`.

## Goal

1. Show a system the way an engineer pictures it: rigid bodies as solid shapes whose
   orientation is visible (body-frame axes, a non-symmetric cross-section, so a spin about
   the link's own axis, i.e. yaw of the body, is obvious), point masses as spheres, cables
   and links as rods, with the trajectory produced by `geomech.simulation`.
2. Put that animation, the derived equations, the standard form and the verification
   numbers into one self-contained HTML report that opens in a browser, in a notebook, or
   as a Claude artifact.
3. Ship a skill so that an agent asked "derive the equations of motion of X" follows one
   procedure: define, derive, verify, standard form, generate, simulate, report.

## Scene description

Geometry is declared once, with geomech expressions, so it is evaluated through the same
code path as everything else (`geomech.codegen.lambdify` on configuration variables and
parameters; no hand-written kinematics in the visualizer).

```python
from geomech.viz import Scene, RigidBody, PointMass, Rod, Box, Cylinder, Sphere

scene = Scene(layout, params=p)
scene.add(RigidBody(R1, origin=Zero3, shape=Box(0.06, 0.03, 0.5), offset=rho1, color="link1"))
scene.add(RigidBody(R2, origin=R1 * l1, shape=Box(0.06, 0.03, 0.5), offset=rho2))
scene.add(PointMass(x_L, radius=0.05))          # an R³ or S²-derived position expression
scene.add(Rod(x_Q, x_L))                        # a cable between two position expressions
scene.add(Trace(R1 * l1 + R2 * (rho2 * 2)))     # path of a point
```

- `RigidBody(R, origin, shape, offset=None, axes=True)`: `R` is an `SO3` configuration
  variable; `origin` a vector expression (the body's reference point in the inertial
  frame); the shape is centred at `origin + R·offset` and rotated by `R`. `axes=True`
  draws the body frame (e₁ red, e₂ green, e₃ blue by convention, length from the shape).
- `S2` variables have no attitude; a pendulum on S² is drawn as `Rod(origin, origin + l q)`
  plus a `PointMass`. A rigid body hanging from an S² point is not representable without
  a body attitude; the docs say so.
- Shapes: `Box(lx, ly, lz)`, `Cylinder(radius, length, axis="z")`, `Sphere(radius)`.
  Sizes are numbers in the scene's units; no meshes.
- Expressions may use configuration variables and parameters only; anything else raises
  `ValueError` naming the leaf (reuse `StateLayout` / codegen validation).
- `scene.frames(traj)` returns, per stored state, each body's position (3) and rotation
  (3×3), each point's position, each rod's endpoints. Evaluating the scene is vectorised
  over states where cheap, otherwise a loop over the generated callables.

## Report

```python
from geomech.viz import Report
report = Report(title, system=(L, dW, variables, inputs), eom=eom, standard_form=sf,
                scene=scene, trajectory=traj, params=p, checks=True)
report.write("double_pendulum.html")      # self-contained HTML file
report._repr_html_()                      # Jupyter display
```

Sections, in order: summary strip (derive time, `check_eom` error, energy drift if the
trajectory is unforced or the drift of whatever `traj.energy()` returns, integrator and
step); 3-D view with play/pause, scrubber, speed, camera reset, toggles for body axes and
traces; Lagrangian; equations of motion one term per line; standard form (block matrix
overview, each inertia block, each f, G); parameter table (name, value, units if given);
the system definition as source (passed in by the caller as a string, optional).

Rendering rules:

- Three.js r147 (`three@0.147.0/build/three.min.js` and
  `examples/js/controls/OrbitControls.js` from cdn.jsdelivr.net, pinned) and MathJax 3.2.2
  (`tex-svg-full`, cdnjs). These are the only external resources, and they are allowed by
  the Claude artifact CSP. No other network access; trajectory data is inlined as JSON.
- Inertial frame: z up, gravity along −e₃ in the docs' convention; a ground grid and the
  inertial axes at the origin; perspective camera with orbit controls; lights fixed in the
  world.
- Theme: light and dark via `prefers-color-scheme` and the `data-theme` attribute, as in
  the artifact contract (token CSS variables; canvas background from a token).
- Responsive down to 400 px; the 3-D view stays square; equations scroll horizontally in
  their own containers.
- LaTeX: `to_latex` output with readable subscripts (`J1` → `J_1`, `\Omega_{R1}` →
  `\Omega_1`) via a documented, tested mapping applied to symbol names only (never by
  regex over whole formulas, which can corrupt them).
- Size: a 10 s run at the default `output_every` must stay under 2 MB of HTML.

## Skill

`skills/geomech/SKILL.md` (Agent Skills format: YAML front matter `name`, `description`)
plus `skills/geomech/reference.md`, and a copy for Claude Code at
`.claude/skills/geomech/` in the repo (or a symlink if the docs build tolerates it).

Content:
- When to use: deriving, checking, linearizing, generating or simulating dynamics of
  rigid bodies, pendulums, aerial vehicles with loads, on R, R³, S², SO(3) and products.
- The procedure, as numbered steps with the exact calls: define symbols and
  `SystemVariables`, write L and δW, `compute_eom`, `check_eom` (mandatory; report the
  number; never present an equation that fails it), `to_standard_form`, optionally
  `linearize`, `StateLayout` + `generate_module`, `simulate`, `Report`.
- Conventions the agent must state when presenting results: E = −(Euler–Lagrange) +
  forcing; M = −physical inertia; S² equations hold on the tangent plane; products
  print with `*`.
- Failure modes and what to do (from `docs/troubleshooting.md`): stray variation,
  attribute typo, input/constant clash, acceleration-dependent Lagrangian, unsupported
  extraction shapes. Report the error text instead of working around it.
- Presentation: always the report or, in a terminal, the equations plus the numbers; never
  a hand-drawn diagram.
- Worked examples: spherical pendulum, double rigid pendulum, quadrotor with a cable load
  (each a short, runnable snippet; a test runs every snippet in the skill).

## Verification (definition of done)

1. For every body in the demo scenes, the position and rotation in `scene.frames(traj)`
   equal an independent NumPy evaluation of the same geometry from the trajectory state
   (1e-12).
2. Rotations in the frames are orthogonal with det +1 (1e-12); a spin about the link axis
   is visible in the data (the body e₁ axis rotates about the link axis in a test case
   with an initial Ω along the link).
3. The generated HTML: contains no external URLs other than the pinned three/MathJax
   ones; parses as HTML (html5lib or the stdlib parser without errors); its inline JSON
   round-trips; size bound holds; MathJax source strings are balanced.
4. A headless browser smoke test if Chrome or Chromium is available on the machine
   (skipped otherwise): load the page, wait for MathJax, assert no console errors and
   that the canvas is non-blank.
5. Examples: `examples/report_double_rigid_pendulum.py` (boxes, initial spin about the
   link axis), `examples/report_quadrotor_load.py` (quadrotor as a box with four rotor
   cylinders, cable, load), `examples/report_spherical_pendulum.py`; each writes its
   report; CI runs them (outputs not committed).
6. Skill: every code snippet in the skill files runs in the test suite; the skill's
   front matter validates.
7. Docs page `docs/visualization.md`, nav entry, README pointer, CHANGELOG; ruff, mypy
   clean; `import geomech` does not import the visualizer's heavy pieces (nothing heavy
   exists, but the module must not import scipy).

## Non-goals

Meshes or CAD import, collision, real-time interaction with the dynamics (dragging bodies),
video export, an MCP server (next step).
