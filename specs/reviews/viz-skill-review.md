# Adversarial review: 3-D visualizer, report and agent skill

Branch `feature/viz-skill`, range `3faf73c..309d5ad`, against `specs/viz-and-skill.md` and
`specs/reviews/viz-skill-implementation-report.md`. Tests: `tests/test_viz_adversarial.py`
(20 passing checks, 11 strict xfails `F<n>`). Screenshots: `viz-review-1-dark-far-from-origin.png`,
`viz-review-2-millimetre-scene.png`. Nothing committed; no existing file modified.

**Verdict: GO WITH CHANGES.** The rendered geometry is right: for the three example scenes
and three of my own, every three.js matrix read back from the live page matches an
independent NumPy pose to the 5-decimal rounding. The z-up world map is right, handedness
and row/column order are right, there is no transpose, and the axes are coloured e₁ red,
e₂ green, e₃ blue. Box spin is visible. Fix F1 (the scene can silently draw a different
system from the one simulated) and F5 (`javascript:` links through MathJax) before
reports are shared. F2–F4 are the places where orientation is still hidden.

## How the page was checked

A probe is injected into the generated HTML. A pre-script wraps `THREE.Scene` and
`THREE.PerspectiveCamera` to capture them, and a post-script pauses playback, drives the
real scrubber (`input` events → the page's own `setFrame`), calls `updateMatrixWorld` and
dumps every body group's `matrixWorld`, the axis-arrow directions and colours, the
point/rod/trace world positions, the trace draw ranges, the time label, the camera and the
inertial triads. `--dump-dom` from headless Chrome returns the result. For playback,
`requestAnimationFrame` is replaced by a queue, and synthetic 60 Hz timestamps are pumped
through the page's own `tick()`, which makes the playback check deterministic. Real-time
headless rAF is throttled and cannot be used for timing.

## Findings, ranked

Ranked by impact: misrepresents the motion > hides orientation > security > crash > skill
errors > cosmetic.

### F1: the scene can draw a different system from the one simulated (misrepresents the motion)

- **Reproducer:** `test_scene_params_must_match_trajectory`. Simulate the spherical
  pendulum with `l = 1`, build `Scene(layout, params=dict(p, l=2.5))`, then `Report(...,
  scene=scene, trajectory=traj)`.
- **Observed:** the bob is drawn at |x| = 2.500. The parameter table (from
  `trajectory.p`) says `l = 1`. There is no warning.
- **Expected:** `Scene.frames(traj)` evaluates with `traj.p` when it is given, or `Report`
  raises `ValueError` when `scene.params` differs from `trajectory.p` on any key the
  geometry uses.
- **Where:** `geomech/viz/scene.py:381-387` (the params are stored independently),
  `scene.py:516` (`p = self.params`) and `geomech/viz/report.py:223-225`.
- **Why it ranks first:** parameter sweeps and copy-paste make this easy to hit. The
  picture then contradicts the equations and the table beside it. A related point:
  `Scene(layout)` without params raises `KeyError` for parameters the geometry does not use
  (`e3`, `g`).

### F2: very small and very large scenes hide the body (hides orientation)

- **Reproducers:**
  - `test_millimetre_scene_is_framed` and screenshot 2: a 1 mm link with an
    `R`-dependent offset.
  - `test_kilometre_body_axes_are_visible`.
- **Observed:**
  - `radius = max(…, 0.5)` (`report.py:365`) frames a 1 mm body in a 1 m view, so it covers
    0.1 % of the canvas, about 0.7 px, and only the inertial axes are visible.
  - Positions are rounded to 1e-5 (`report.py:44`), which is 1 % of a 1 mm part.
  - Body-axis arrows have an absolute radius of `clamp(0.02·len, 0.004, 0.03)`
    (`template.html:370`). On a 1 km body they are 3 cm wide over 750 m, ratio 4e-5, which
    is sub-pixel.
  - The default `PointMass` radius (0.05) and `Rod` radius (0.008) are also absolute
    (`scene.py:192, 209`).
- **Expected:**
  - Framing from the scene's own size, with no metre floor.
  - Rounding relative to the scene radius.
  - Arrow and default radii as fractions of the scene or shape size.
- **Holds:** 1 m and 1 km framing (`test_metre_and_kilometre_scene_framing`).

### F3: spin of a Cylinder or Sphere body is invisible with its axes hidden (hides orientation)

- **Reproducer:** `test_cylinder_and_sphere_spin_is_visible_without_axes`. The probe checks
  that the drawn body, without its axes group, is a single colour and that its vertex set
  is unchanged by a quarter turn about body e₃.
- **Observed:** the cylinder (40 segments; `EdgesGeometry(geo, 30)` keeps only the two cap
  circles) and the sphere are single-material and quarter-turn invariant. With
  `axes=False`, which the docs recommend for all but one shape sharing an `R`, or with the
  "Body axes" toggle off, a spinning wheel, rotor, top or ball looks static.
- **Where:** `template.html:353-361` (only boxes get tinted faces) and `template.html:366`.
- **Expected:** a seam or tinted sector on cylinders and a marked meridian or band on
  spheres, so that "unmistakable spin" holds for every shape, not only boxes. Boxes hold:
  `test_box_faces_and_axis_colours` (≥ 3 face colours, not quarter-turn invariant).

### F4: far from the origin there is no inertial frame in view (hides orientation)

- **Reproducer:** `test_inertial_frame_in_view_far_from_origin` and screenshot 1. A free
  body is placed at (50, −30, 20).
- **Observed:** the only inertial triad sits at the origin (`template.html:326-336`), about
  62 units from the camera, beyond `far = 100 R ≈ 60` (`template.html:409`). It is not in
  the frustum. The grid does not show which way x points, so the body's yaw relative to
  the inertial frame cannot be read. The near-origin case holds
  (`test_inertial_frame_in_view_near_origin`).
- **Expected:** an inertial triad near the scene (for example at the grid corner) or a
  fixed view gizmo.

### F5: `\href{javascript:…}` in symbol names, title, description and units becomes clickable links (security)

- **Reproducer:** `test_no_javascript_links_after_typesetting`. A constant named
  `\href{javascript:alert(1)}{click}`, plus the same text in `title`, `description` and
  `units`.
- **Observed:** 20 `javascript:` hrefs in the DOM after MathJax. One click runs script in
  the page origin.
- **Where:** `template.html:19-30`. MathJax is configured without `ui/safe` and typesets the
  whole body, including `<h1>`, the subtitle and the units cells (`report.py:543,
  584-587`).
- **Expected:**
  - Load `ui/safe` (`loader: {load: ["ui/safe"]}`), or drop the `html` TeX package.
  - Restrict typesetting to `.math` (`options: {ignoreHtmlClass: "…", processHtmlClass:
    "math"}`).
- **What held** (`test_hostile_strings_are_escaped`, `test_hostile_strings_do_not_execute`):
  - HTML escaping: `</title><script>` in the title, `<img onerror>` in the description and
    `</code></pre><script>` in the source are all inert.
  - JSON escaping: `</script>` in a body name becomes `</script>` and
    round-trips.

### F6: non-finite states pass the rotation check, then the report crashes (crash)

- **Reproducer:** `test_non_finite_state_is_rejected_by_frames`. One NaN in an `R` entry.
- **Observed:**
  - `Scene.frames` returns NaN rotations. `err > tol` and `|det − 1| > tol` are False for
    NaN (`scene.py:521-527`), and under `-W error` `np.linalg.det` raises a
    `RuntimeWarning` instead.
  - With NaN in the second half of a trajectory, `Report.html()` fails in `json.dumps`:
    `ValueError: Out of range float values are not JSON compliant: nan`
    (`report.py:65`). The error names no state and no time.
- **Expected:** a `ValueError` naming the first non-finite state and its time, or a report
  truncated at the divergence that says so. A diverged simulation is the case where the
  picture is needed most.

### F7: skill and repository-hygiene errors (skill errors)

- **Print convention.** `SKILL.md:53` says "Printed products use `*`", and `SKILL.md:73`
  says to present `print_eom(eom)` in a terminal. `print_eom` prints LaTeX (`\int \xi_{q}
  \cdot (… -m l …) dt`); `*` appears only in `str(expr)`. An agent would state a convention
  that the output it shows contradicts (`test_skill_print_convention_matches_print_eom`).
- **Scene parameters.** The skill never says that the `Scene` must be built with the same
  `p` that was simulated. This is the trap behind F1. All three worked examples happen to
  reuse `p` (`SKILL.md:112, 167, 213`).
- **Forced energy.** The quadrotor worked example is forced (hover thrust). The summary
  card still reads "Energy drift" with no "forced" marker: the report cannot know. The
  skill tells the agent to check drift "(unforced)" (`SKILL.md:42`) but not to label the
  forced number as energy change rather than error, so it can be misreported as an
  integrator problem.
- **`.gitignore`.** `.gitignore:17-18` (`.claude/*` + `!.claude/skills/`) un-ignores
  every future `.claude/skills/<x>`, not just the geomech copy. A personal skill placed
  there becomes committable to the public repo
  (`test_gitignore_unignores_only_the_geomech_skill`). `settings.local.json`,
  `.claude/agents` and `.claude/commands` stay ignored (`test_gitignore_keeps_settings_ignored`).
  Narrow it to `!.claude/skills/geomech/`.

### F8: the name map produces invalid TeX for argument-taking commands (cosmetic)

- **Observed:** `pretty_name("\\bar1") == "\\bar_1"`, and likewise `\hat1`, `\vec1` and
  `\mathbf1`. MathJax renders "Missing superscript or subscript argument" (6 `merror`s in a
  test page), where the unmapped `\bar1` was valid.
- **Where:** `latex.py:34` (`\\[A-Za-z]+` accepts any command).
- **Reproducer:** `test_latex_map_keeps_argument_commands_valid`.

### F9: playback and framing details (cosmetic)

- **Last frame.** The last stored state is never shown while playing. The clock wraps when
  `simTime > t[N−1]` before `frameAt` can return `N−1` (`template.html:524`). Pumped at
  60 Hz the frames run 599 → 0, and frame 600 (t = 10.00 s) is never displayed
  (`test_playback_shows_last_frame`).
- **Lag and single-state pages.** There is a one-frame lag: at wall 1.000 s the page shows
  t = 0.98 s. With `N = 1`, pressing Play loops `requestAnimationFrame` forever
  (`template.html:533`).
- **Arrow clipping.** Framing (`report.py:365`) ignores the axis-arrow lengths: in
  screenshot 1 the e₁ arrow is cut at the canvas edge.
- **Tints that vanish.** Face tints depend on the body colour. Measured as ΔE (CIE76):
  - A body coloured `load` has a +e₁ tint of ΔE 3.5 in light theme and 4.0 in dark: it
    vanishes.
  - `link3` has a +e₂ tint of ΔE 4.1.
  - In every case one of the two tinted faces stays at ΔE ≥ 19, so box spin about any
    axis stays visible except for a load-coloured cube turning about e₂ with its axes off.

## What held (with counts)

- **Geometry fidelity, 1.** `test_rendered_poses_match_numpy` covers 5 scenes: the double
  rigid pendulum, quadrotor with load (5 shapes sharing one `R`), a free box spinning at
  9 rad/s about e₃ (not its long axis) far from the origin, two chained S² pendulums, and
  the new skill-built SO(3)×S² system. Seven frames each, including 0, 1 and the last two.
  - Rotation `Pᵀ·M₃ = R` and position = `origin + R·offset`: worst 5.0e-6.
  - Body-arrow directions `= R e_j`: worst 5.6e-6.
  - Points, rods and traces: worst 4.9e-6. Trace `drawRange = k + 1`.
  - Time label equals `traj.t` to 0.005 s.
  - The world matrix equals `[[1,0,0],[0,0,1],[0,−1,0]]` exactly (proper, z-up), and
    `camera.up` is inertial e₃.
  - The inertial arrows point along e₁, e₂, e₃ in the axis-token colours. Each body's
    arrows use the same colours.

  A first exploratory pass over the three full 10 s examples (601 frames, six frames each)
  gave the same worst errors.
- **Spin in the data.** The free box's long axis sweeps more than two turns
  (`test_spin_about_non_long_axis_is_in_the_data`); across the full 3 s run it sweeps
  1545°.
- **Time, 3.** Playback is real time × speed for speeds 1 and 2, within one frame
  (`test_playback_is_real_time_times_speed`). Thinning keeps the first and last state with
  no gap above 1/fps + dt. Frames are picked, never interpolated, as documented.
- **LaTeX, 4.**
  - Collisions drop both entries: `J1` with `J_1`, and `\omega1` with `\omega_{q1}`.
  - Unicode, `^`, `*` and already-subscripted names are unchanged.
  - The map is one-to-one on the new system.
  - The default `to_latex(expr)` is byte-identical to `3faf73c` for 50 pinned systems
    (388 strings), checked against a `git archive` of the base
    (`test_to_latex_default_byte_identical_to_base`).
- **CSP and resources, 5.** Under a strict artifact-style CSP (only the three pinned
  scripts, `connect-src 'none'`) there are no violations and no errors. MathJax finishes,
  the canvas is drawn (34 colours), and the only resources fetched are the three pinned
  URLs (`test_page_under_artifact_csp`). The dark theme renders: canvas, grid and body
  tokens switch (screenshot 1).
- **Degenerate scenes, 6.** An empty scene, a traces-only scene, a one-state trajectory,
  and Cylinder and Sphere bodies all load with 0 errors (`test_degenerate_scenes_load_without_errors`).
  A 120 s run is 0.78 MB.
- **Skill, 7.**
  - Front matter is valid for Agent Skills: `name` matches the directory, `[a-z0-9-]`,
    ≤ 64 chars; `description` is 461 chars with no angle brackets. The `.claude/skills`
    copy is byte-identical.
  - Following SKILL.md alone, the new system (an S² pendulum on a swinging, yawing rigid
    body) derived, passed `check_eom` at 2.1e-8 with energy drift 6e-11, and produced a
    correct report on the first try (`test_skill_procedure_on_new_system`).
  - The failure-mode texts in the skill match the real messages: stray variation,
    attribute typo, input/constant clash, acceleration-dependent `L`.
  - `check_eom` is mandatory and asserted in every worked example, and `Report` runs it by
    default; `checks=False` shows "not run" in red.
- **Regression, 8.**
  - `pytest tests/ -W error`: 3367 passed, 13 xfailed (3346 + 2 before this file).
  - `ruff format --check` and `ruff check` on the CI paths: clean.
  - `mypy geomech`: clean.
  - `mkdocs build --strict`: builds.
  - The three examples ran with `-W error`, `check_eom` at 5.0e-8, 4.0e-8 and 2.3e-8,
    sizes 175, 259 and 85 kB.

## Fix log

Fixes on `feature/viz-skill` on top of `49765b3` (not committed). Every strict xfail but
one was flipped by deleting its marker only; no assertion was changed.

| Finding | Fix (file:line) | Test |
|---|---|---|
| F1 | `Scene(layout)` needs no params; `Scene.params_for` uses `traj.p` (`geomech/viz/scene.py:609`), `check_same_params` raises `ValueError` naming the parameter and both values (`scene.py:365`), called by `Report` for `params=` and `Scene(params=)` (`geomech/viz/report.py:244-246`). Parameter offsets are evaluated lazily with those params (`_Offset`, `scene.py:329`). Docs: `docs/visualization.md` "Parameters: one source of truth"; skill step 9. | `test_scene_params_must_match_trajectory` (flipped); `test_viz_fixes.py`: `test_scene_without_params_uses_trajectory_p`, `test_conflicting_params_name_the_parameter_and_both_values[scene,report]`, `test_same_params_given_twice_are_accepted`, `test_bare_states_need_params`, `test_parameter_offset_is_evaluated_with_trajectory_p` |
| F2 | `scene_data` (`report.py:382`): bounds widened by each element's reach (shape half-diagonal, axis-arrow length, radii; `report.py:420`), radius from the actual extent with no floor (`report.py:430`); positions to `_POS_SIG = 6` significant figures of the radius (`report.py:47, 83`); default `PointMass`/`Rod` radii `None` → 6 % / 1 % of the scene half-size, capped at 35 % / 6 % of the largest body (`scene.py:196, 218`, `report.py:401`); arrow radius `0.02 × arrow length`, no clamp (`template.html:457`). | `test_millimetre_scene_is_framed`, `test_kilometre_body_axes_are_visible` (flipped); `test_scene_data_is_scale_invariant`, `test_millimetre_positions_are_not_rounded_away`, `test_default_point_and_rod_radii_follow_the_scene`, `test_framing_includes_the_body_axis_arrows` |
| F3 | `marker()` (`template.html:410`): a cylinder gets a tinted 60° sector of side and caps on body e₁ (e₂ for `axis="x"`), a sphere caps on e₁ and e₂ (`template.html:436, 444`); independent of the axes toggle. | `test_cylinder_and_sphere_spin_is_visible_without_axes` (flipped); `test_marker_direction_follows_body_axis_with_axes_hidden` reads each marker's world direction back at two frames with the axes toggled off and compares it with `R e_j` (atol 2e-3, body turned > 30° between the frames) |
| F4 | Corner gizmo: its own scene created before the main one (`template.html:245`), rendered in a bottom-left scissor viewport, turned with the camera (`renderGizmo`, `template.html:630`); a second inertial triad on the ground behind the scene when the origin is more than 3 radii away (`template.html:355`). | `test_inertial_frame_in_view_far_from_origin` (flipped); `test_corner_gizmo_is_an_inertial_triad` |
| F5 | MathJax without `html`, `action`, `bbox`, `unicode`, `require`, `autoload`, `setoptions` (`template.html:27`), typesetting only `.math` (`template.html:30`, `<body class="gm-nomath">` at `:190`). `ui/safe` was not used: it is not in `tex-svg-full` and loading it would fetch a fourth URL, breaking the artifact-CSP test. | `test_no_javascript_links_after_typesetting` (flipped: 20 → 0 links); `test_tex_injection_strings_stay_inert` (11 strings: `\href`, `\require{html}\href`, `\text{\href}`, `\class`, `\cssId`, `\style`, `\bbox`, `\unicode`, `\toggle`, `\texttip`, `\setOptions` in names, title, description and units: 0 links, 0 attacker classes/ids/styles, nothing typeset in header or units, equations still typeset); `test_review_hostile_report_has_no_links` |
| F6 | `check_finite_states` (`scene.py:393`) in `Scene.frames` (`scene.py:640`) and `Report.__init__` (`report.py:247`): `ValueError` naming the first non-finite time, stored-state index and entry (`R[1, 1] = inf`); `np.isfinite` only, so no RuntimeWarning under `-W error`. | `test_non_finite_state_is_rejected_by_frames` (flipped); `test_non_finite_state_names_time_and_slot` |
| F7 | `SKILL.md:42-52, 60-62` and the presentation section: `print_eom` prints LaTeX, `str` the `*` form; the scene is drawn with `traj.p`; a forced run's energy change is the work of the inputs. `Trajectory.u` (inputs at the stored states, `geomech/simulation/integrators.py:69, 257`) lets `Report.forced()` (`report.py:356`) label the card "Energy change … the work of the inputs, not an integrator error" (`report.py:332`). `.gitignore:17-20` un-ignores `.claude/skills/geomech/` only; `.claude/settings.local.json`, `settings.json`, `agents/`, `commands/`, `skills/private/` and `skills/x.md` checked ignored with `git check-ignore`. `.claude/skills/geomech` is byte-identical to `skills/geomech`. | `test_gitignore_unignores_only_the_geomech_skill` (flipped), `test_gitignore_keeps_settings_ignored`; `test_forced_run_labels_energy_change`; the skill snippets run (`test_viz_skill.py`) |
| F8 | `latex.py:51-67`: only letters or a symbol command (Greek letters, `\ell`, `\hbar`, `\partial`, …: `SYMBOL_COMMANDS`) followed by digits are subscripted; `\bar1`, `\hat1`, `\vec1`, `\mathbf1`, unknown commands stay unchanged. | `test_latex_map_keeps_argument_commands_valid` (flipped); `test_every_mapped_name_is_valid_tex` (balanced braces, `_`/`^` only after letters or known symbol commands); `test_every_mapped_name_typesets_in_mathjax` (MathJax 3.2.2 in headless Chrome: each mapped name typesets with 0 errors) |
| F9 | Past the end, the clock first shows the last stored state, then wraps (`template.html:616`); Play and the scrubber are disabled for one state (`template.html:554`); framing includes arrows (F2); tints via `contrastTint` (`template.html:373-405`): first mix towards the axis colour with ΔE ≥ 25, else the axis hue shifted in lightness away from the body. The one-frame display lag of playback is unchanged. | `test_playback_shows_last_frame` (flipped); `test_one_state_trajectory_does_not_loop`; `test_orientation_tints_stand_out_on_every_palette_colour[light,dark]` (box faces and markers on all 8 palette colours: worst ΔE ≥ 24.5, was 3.5) |

### Not flipped: `test_skill_print_convention_matches_print_eom`

Its assertions require both that SKILL.md keeps the sentence "Printed products use `*`"
and that `print_eom(eom)` output contains `*`. The requested fix is in the skill (say
what `print_eom` prints versus `str`), which removes that sentence; `print_eom` still
prints LaTeX, by design (`geomech/utils/printing.py:39`). Making the test pass would mean
changing `print_eom`'s output, which the brief did not ask for and which other tests and
the docs rely on. The assertion encodes the wrong fix, so the strict xfail is kept as is
(it still xfails) and is reported here instead of being edited.

### Checks

- Screenshots, each looked at: `viz-double-rigid-pendulum.png` (regenerated, full report
  at t = 1.5 s: corner gizmo, links, trace, arrows framed), `viz-fix-millimetre-scene.png`
  and `viz-fix-kilometre-scene.png` (the same picture, apart from the title, of a 1 mm and a
  1 km link), `viz-fix-far-from-origin-dark.png` (body at (50, −30, 20): local triad
  under the body and the corner gizmo), `viz-fix-markers-frame0.png` and
  `viz-fix-markers-frame6.png` (rotor, x-cylinder and sphere with axes off; the sector
  and caps move between the frames), `viz-fix-quadrotor-forced.png` ("Energy change"
  card, rotor sectors).
- New sanity check: the 6 new browser tests (markers, gizmo, injection, one-state, tints
  ×2) fail against the `49765b3` template and pass with the new one.
- `pytest tests/ -W error`: **3398 passed, 3 xfailed** (was 3367 passed, 13 xfailed:
  10 xfails flipped, 21 new tests in `tests/test_viz_fixes.py`; the 3 xfails are the two
  older ones outside this review and the print-convention test above).
- `ruff format --check` and `ruff check` on `geomech/ tests/ scripts/
  examples/report_*.py`: clean.
- `mypy geomech`: 2 errors, both "Library stubs not installed for scipy" in
  `simulation/integrators.py:433` and `simulation/control.py:197`; identical on the
  unmodified `49765b3` tree in this venv (environment, not this change).
- `mkdocs build --strict`: builds.
- Examples with `-W error`: double rigid pendulum `check_eom` 5.0e-8, 181 kB, 7.2 s;
  spherical pendulum 2.3e-8, 91 kB, 2.7 s; quadrotor with load 4.0e-8, 287 kB, 6.3 s
  (summary reads "Energy change"). Sizes grew 3-11 % from the extra position digits for
  sub-metre scenes; the 2 MB bound holds.

### Arbitration (orchestrator)

`test_skill_print_convention_matches_print_eom` asserted the skill's old wording ("Printed
products use `*`") together with `*` in `print_eom` output; the agreed fix corrected the
skill instead (print_eom prints LaTeX source, `str` of an equation gives the `*` form). The
test now asserts both statements in the skill and checks each against the code; its xfail
marker is removed.
