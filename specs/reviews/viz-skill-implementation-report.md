# Implementation report: 3-D visualizer, report and agent skill

Branch `feature/viz-skill` (from 3faf73c), spec `specs/viz-and-skill.md`. Nothing is
committed.

## What was built

| Piece | Files |
|---|---|
| Scene | `geomech/viz/scene.py`: `Scene`, `RigidBody`, `PointMass`, `Rod`, `Trace`, `Box`, `Cylinder`, `Sphere`, `SceneFrames` |
| Report | `geomech/viz/report.py`: `Report` (`write`, `html`, `_repr_html_`, `stats`, `energy_drift`, `frame_indices`, `scene_data`), `EXTERNAL_URLS` |
| HTML/JS template | `geomech/viz/template.html`, package data (`[tool.setuptools.package-data] "geomech.viz" = ["template.html"]`), read with `importlib.resources` |
| LaTeX names | `geomech/viz/latex.py`: `pretty_name`, `latex_names`; `to_latex(expr, names=None)` in `geomech/utils/printing.py` |
| Examples | `examples/report_double_rigid_pendulum.py`, `examples/report_spherical_pendulum.py`, `examples/report_quadrotor_load.py` |
| Skill | `skills/geomech/SKILL.md`, `skills/geomech/reference.md`, identical copy in `.claude/skills/geomech/` |
| Tests | `tests/test_viz.py` (38), `tests/test_viz_skill.py` (9) |
| Docs | `docs/visualization.md` (nav "Visualization and Reports"), README pointer and package tree, CHANGELOG Unreleased/Added |
| CI | ruff format/check include `examples/report_*.py`; a step in the 3.12 test job runs the three examples into a temp dir and deletes it |
| Evidence | `specs/reviews/viz-double-rigid-pendulum.png` (headless Chrome, 1000×2600) |

## Public API

```python
from geomech.viz import Scene, RigidBody, PointMass, Rod, Trace, Box, Cylinder, Sphere
from geomech.viz import Report, latex_names, pretty_name

Scene(layout, params=None).add(element) -> element
Scene.frames(traj_or_states, indices=None) -> SceneFrames   # t, bodies[{position, rotation, center}], points, rods (N,2,3), traces
RigidBody(R, origin, shape, offset=None, axes=True, color=None, name=None)
PointMass(position, radius=0.05, color=None, name=None)
Rod(a, b, radius=0.008, color=None, name=None)
Trace(point, color=None, name=None)
Box(lx, ly, lz); Cylinder(radius, length, axis="z"); Sphere(radius)
Report(title, system=(L, dW, variables, inputs), eom=None, standard_form=None, scene=None,
       trajectory=None, params=None, checks=True, *, units=None, source=None,
       derive_time=None, fps=60.0, description="")
Report.write(path) -> Path; Report.html() -> str; Report._repr_html_() -> str
to_latex(expr, names=None)          # backwards compatible
```

Expressions are validated when added: every leaf must be a configuration variable or a
parameter of the layout (or a literal); velocities, inputs and unknown symbols raise
`ValueError` naming the leaf. Evaluation: a bare configuration slot is sliced from the
state (vectorised); everything else goes through `geomech.codegen.lambdify` and is
called per state.

## Decisions not fixed by the spec

- **Animation sampling.** The page keeps the first stored state at or after each
  multiple of `1/fps` (default 60 frames per simulated second) plus the last; no
  interpolation. 10 s → 601 frames. Values rounded to 5 decimals. Rotations are stored
  once per distinct `R` expression (the quadrotor's 7 shapes share one track).
- **World rotation.** All inertial objects sit in one three.js group rotated −90° about
  x (inertial `(x, y, z)` → three `(x, z, −y)`); camera, target and grid use the same
  map. Lights are in the world.
- **Orientation cues.** Solid `MeshStandardMaterial` bodies with dark edges, the +e₁
  face of a box tinted towards red (+e₂ slightly towards green), body-frame arrows whose
  length is half the shape along that axis plus a margin `max(0.25·extent, 4·smallest
  dimension)` so all three stick out, drawn at the shape centre; directional light with
  PCF soft shadows on a shadow-catcher plane, hemisphere fill light.
- **Ground.** At z = 0 when the lowest point of the scene is between 0 and 2 scene radii
  above it, otherwise just below the lowest point; grid step a 1-2-5 number near R/4.
- **One body = one shape.** The API has one shape per `RigidBody`; multi-shape bodies
  (quadrotor) are several `RigidBody` entries with the same `R`/origin and different
  offsets, `axes=False` on all but one.
- **`origin`/`offset`.** `origin` also accepts `None` or 3 numbers (fixed point);
  `offset` accepts 3 numbers or a constant vector expression (configuration-dependent
  offsets raise). `R` may be any matrix expression; `frames` raises if it is not a
  rotation to 1e-6.
- **Colours.** Palette names (`link1`–`link4`, `load`, `cable`, `rotor`, `trace`)
  resolve to theme CSS variables `--viz-<name>` at runtime, so bodies follow
  light/dark; `#rrggbb` also accepted.
- **Summary.** "Energy drift" is always `max |E − E₀| / |E₀|` of `traj.energy()` (the
  spec's "or whatever traj.energy() returns"); for the forced quadrotor the description
  says so. `check_eom` errors above `1e-5` are shown red with "FAILS". Derivation time is
  measured when `Report` derives; otherwise `derive_time=` (examples pass it).
- **Standard form.** Blocks named `\mathbf{M}_{ij}`, `\mathbf{f}_i`, `\mathbf{G}_{ik}`
  (bold, to avoid a clash with an input called `M`); the overview is three `bmatrix`
  products; every block printed with one addend per line.
- **LaTeX map.** Rule set: letters/one command + digits → subscript; the tangent and
  variation vectors of a point named letters + digits take its index
  (`\Omega_{R1}` → `\Omega_1`); everything else unchanged; any collision drops all
  colliding entries (one-to-one). `to_latex` carries the map in a `ContextVar` so the
  recursive handlers inherit it; default output byte-identical.
- **Jupyter.** `_repr_html_` wraps the page in `<iframe srcdoc>` so MathJax configs don't
  clash.
- **Smoke-test hooks.** The page records `data-gm-errors` (window errors, unhandled
  rejections, `console.error`), `data-gm-mathjax="done"` and `data-gm-canvas-colors`
  (distinct colours on a 31×31 sample of the canvas after the first render) on `<html>`,
  so `--dump-dom` is enough to assert "no console errors" and "canvas non-blank".
  The renderer uses `preserveDrawingBuffer: true` for this.
- **Theme switch.** A small "Theme: auto/light/dark" button sets `data-theme` (not
  required by the spec, but it exercises the attribute path).
- **`.claude/skills`.** `.gitignore` ignored `.claude/`; changed to `.claude/*` +
  `!.claude/skills/` so the copy can be committed while `settings.local.json` stays
  ignored. Copy, not symlink; a test asserts byte equality.
- **Skill snippets** run with `python -W error -c` in a temp cwd with the repo on
  `PYTHONPATH`; simulations in them are 1–3 s to keep the suite fast.

## Verification

| Spec item | Test / evidence |
|---|---|
| 1. frames vs NumPy (1e-12) | `test_double_pendulum_frames_match_numpy`, `test_quadrotor_load_frames_match_numpy`, `test_spherical_pendulum_frames_match_numpy` |
| 2. orthogonal, det +1; spin visible | `test_example_rotations_are_orthogonal_and_spin_is_visible` (example: Ω₁·e₃ integrates to > 1 turn); `test_spin_about_link_axis` (hanging pendulum with Ω along the link: `R e₁ = (cos wt, sin wt, 0)` to 1e-9, link axis fixed) |
| 3. HTML: URLs, parses, JSON, size, balanced TeX | `test_only_pinned_external_urls`, `test_html_parses` (stdlib parser + tag balance), `test_inline_json_round_trips`, `test_size_bound`, `test_mathjax_sources_balanced`, plus section order, one-term-per-line, failing check shown red |
| 4. headless browser | `test_headless_browser_smoke` (Chrome `--headless=new --disable-gpu --dump-dom`): no errors, MathJax done, >20 `mjx-container`, canvas present, ≥ 5 canvas colours; skipped without Chrome or without network to the three CDN URLs. It ran (not skipped) here. |
| 5. examples | the three scripts; CI step runs them |
| 6. skill | `tests/test_viz_skill.py`: front matter, identical copy, every fenced python block (5) runs |
| 7. docs, ruff, mypy, imports | `docs/visualization.md`, `test_imports_stay_light` (`import geomech` does not import `geomech.viz`; `geomech.viz` does not import scipy) |

Commands run (all clean):

- `ruff format --check` and `ruff check` on `geomech/ tests/ scripts/ examples/generated/
  examples/inverted_rigid_pendulum_lqr.py examples/report_*.py`
- `mypy geomech --ignore-missing-imports`: no issues in 50 files
- `python -m pytest tests/ -q -W error`: 3346 passed, 2 xfailed (pre-existing), 75 s; 47 of them new (38 in `test_viz.py`, 9 in `test_viz_skill.py`). The viz, skill and printing tests were rerun after the last edits (99 passed)
- `mkdocs build --strict`: built (after adding annotations griffe asked for)
- `python scripts/regenerate_generated.py --check`: 18 modules up to date

## Sizes and times (this machine)

| Example | Total run | of which simulation | `Report.html()` | HTML size |
|---|---|---|---|---|
| double rigid pendulum (10 s, 10 001 states) | 5.3 s | ≈4.8 s | 0.24 s | 175 kB |
| spherical pendulum (10 s) | 2.5 s | | | 85 kB |
| quadrotor with load, LQR (10 s, 8 shapes) | 5.7 s | | | 259 kB |

All far below the 2 MB bound. Derivation times shown in the reports: 0.02 s (double
pendulum), 0.002 s (spherical), 0.02 s (quadrotor).

## Screenshot

`specs/reviews/viz-double-rigid-pendulum.png`. Checked by eye: both links are solid
boxes with lighting, edges and shadows, the tinted +e₁ face and the three body-frame
arrows are visible on each link, the inertial axes with x/y/z labels and the ground grid
are present, the tip trace is drawn, and every equation (L, δW, both EOMs, the standard
form) is typeset by MathJax with readable subscripts (`J_1`, `\Omega_1`, `R_1^T R_2`,
`m_1`, `e_3`). A narrow (420 px) dark-mode render was also checked: layout reflows, the
summary cards wrap two per row, the view stays square and the theme colours apply to
canvas, grid and text.

## Not done / caveats

- The smoke test uses the page's own hooks (`data-gm-*`) rather than Chrome's console
  log; an error thrown before the first inline script (impossible in practice) would be
  missed.
- WebGL in headless Chrome here ran on a GPU path despite `--disable-gpu`; on a machine
  where WebGL is unavailable the page shows "WebGL is not available" and the smoke test
  would fail on the canvas-colour assertion rather than skip.
- `examples/report_quadrotor_load.py` needs scipy (`lqr`); CI installs it via `[dev]`.
- No pixel-level check of the screenshot in the tests (no PIL dependency); the canvas
  colour count stands in for "non-blank".
