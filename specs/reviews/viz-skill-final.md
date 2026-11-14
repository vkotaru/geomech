# Final review: 3-D visualizer, report and agent skill

Branch `feature/viz-skill`, range `a4961bb..18b72fe` (spec, implementation, adversarial
review, fixes), for a PR into `main`. Read: `specs/viz-and-skill.md`, the implementation
report, the adversarial review with its fix log and arbitration, and the full diff of
`geomech/ skills/ .claude/ docs/ README.md CHANGELOG.md .github/ .gitignore pyproject.toml
examples/`. Screenshots: `viz-final-light-420.png`, `viz-final-dark-420.png`. Nothing
committed. No source, test or skill file was changed by this review.

## Verdict: GO WITH CHANGES

The library code, docs, skill and examples are ready, and every gate is green on this
machine. One change is required before the PR, and it is in the tests, not the library:
the headless-browser tests will run in CI, because GitHub's `ubuntu-latest` image has
Google Chrome and network access. They depend on live CDNs and on software WebGL, and
when WebGL is missing they fail instead of skipping. On a protected branch with three
required test jobs, that makes them a source of blocked merges.

## Required changes

**R1. Browser tests must skip, not fail, when the browser cannot draw.**
- **Where:** `tests/test_viz.py:518` and `tests/test_viz_adversarial.py:430`. The second
  is the `browser` fixture, which `tests/test_viz_fixes.py` reuses.
- **Evidence:** I ran Chrome with `--disable-gpu --disable-software-rasterizer`
  (no WebGL). The page degrades correctly ("WebGL is not available in this browser…"),
  but `three.js` calls `console.error` 7 times, so the root reads
  `data-gm-errors="7"`. `test_headless_browser_smoke` then fails on
  `'data-gm-errors="0"' in root`, and the probe-based tests fail too.
- **Today:** Chrome 154 still falls back to SwiftShader under `--disable-gpu` (checked:
  `ANGLE (… SwiftShader driver)`), so CI probably passes. Chrome has deprecated that
  automatic fallback, though.
- **Flakes:** in 11 runs of the browser test files, two runs failed, each on a
  different test:
  - `test_page_under_artifact_csp` under Python 3.10, which depends on the network and
    wall-clock time;
  - a probe test in `tests/test_viz_adversarial.py`, whose injected probe threw
    `TypeError: Cannot read properties of undefined (reading 'children')` at `world`
    (`tests/test_viz_adversarial.py:283`,
    `sc.children.filter(…isGroup)[0]` came back undefined).

  The probe failure is probably a race between the captured `THREE.Scene` and the corner
  gizmo's scene, which is now created first. The rest of each run passed.
- **Change:**
  1. In both commands, add `"--enable-unsafe-swiftshader"` after `"--disable-gpu"`.
  2. After the run, skip when the page could not draw: in `test_viz.py`, after
     `dom = proc.stdout`, add
     `if "WebGL is not available" in dom: pytest.skip("no WebGL in this headless Chrome")`.
     Do the same in the fixture, on `proc.stdout`, before the probe is parsed.
  3. Choose one: run the browser tests in the 3.12 job only (an env flag such as
     `GEOMECH_BROWSER_TESTS=1`, set by that job), or keep them in all three jobs and
     accept the network flake on required checks.

## Recommended (not blocking)

1. **The test pinned to a commit outside the branch.**
   - `tests/test_viz_adversarial.py:670,674,677` archive `3faf73c`. That commit is a
     pre-rebase base: it is not an ancestor of `HEAD`, it exists only in this clone's
     object store, and CI's shallow checkout skips the test.
   - Pin `a4961bb` (the merge base on `main`) instead, so the test compares against
     history the branch actually has.
2. **"Unforced" stated without evidence.**
   - `geomech/viz/report.py:345` labels the card "unforced" whenever `trajectory.u` is
     `None`. That happens when a callable `u` raised in `_stored_inputs`, or when a
     `Trajectory` was built by hand. In those cases the inputs are unknown, not zero.
   - Suggest `note="max |E − E₀| / |E₀|; inputs not recorded"` when `u is None`.
3. **`_stored_inputs` can crash `simulate` after the integration has finished.**
   - `geomech/simulation/integrators.py:261` catches only
     `ValueError, TypeError, ArithmeticError`. A controller that raises `KeyError` or
     `IndexError` at a stored state fails `simulate` after the whole integration.
   - It also calls the controller a second time at every stored state. That is wrong for
     a stateful or stochastic controller.
   - Either document both behaviours, or catch `Exception`.
4. **`reference.md:111-113` lists the `Report` signature without the `*`.**
   - The signature shows `units=None, source=None, derive_time=None, fps=60`, but these
     are keyword-only.
   - Write `checks=True, *, units=None, …`, as `docs/visualization.md` does.
5. **Skill gaps found by the agent walkthrough (below):**
   - the `layout.pack` keys for indexed variables;
   - the virtual work of a force at a point that depends on S²;
   - how to check the integrator on a forced run.
6. **Visual polish, from the screenshots (below):**
   - the grid is faint in the light theme;
   - the corner-gizmo labels are about 8 px at 420 px wide;
   - the e₃ arrow of a body hinged at the origin lies along the inertial z arrow, so two
     blue arrowheads stack.
7. **`specs/follow-ups.md` records none of the deferred items listed under Open items.**
   Add them.
8. **CHANGELOG.** The "Visualizer and skill fixes after the adversarial review" bullet
   lists fixes to features that were never released, under "Added". Fold them into the
   feature bullets, or keep them as they are on purpose.

## Screenshots

`examples/report_double_rigid_pendulum.py` (run with `-W error`), opened in headless
Chrome 154 (`--headless=new --disable-gpu`, 420×2200 px, 8 s virtual time). The light
screenshot uses the default colour scheme. The dark one uses `--force-dark-mode`, which
makes `prefers-color-scheme: dark` match; the theme button stays on "auto". The animation
was playing: t = 0.57 s (light) and 0.68 s (dark).

- **Bodies are solid.**
  - Link 1 is a blue box and link 2 a tan/orange box. Both are lit and shaded, with dark
    edges.
  - Each box casts a shadow on the ground: soft grey in the light theme, black in the
    dark one. The arrows cast shadows too.
- **Orientation cues.**
  - On link 2, a narrow red-tinted strip shows along one long face: the +e₁ face. At this
    camera angle link 1's tinted face is turned away.
  - Body arrows: link 1's e₂ (green) sticks out clearly. Its e₁ (red) is foreshortened
    to a red dot at the box edge. Its e₃ (blue) points up along the link, collinear with
    the inertial z arrow at the pivot, so two blue arrowheads stack.
  - Link 2's e₁ and e₂ arrows show at its centre. Its e₃ arrow pokes out near the joint.
  - One still frame cannot show the spin itself. The tinted face and the arrows that it
    relies on are present, and the adversarial review checked the poses in the page
    against NumPy.
- **Inertial frame.**
  - The origin triad is labelled x, y, z in red, green and blue. It is large and
    readable.
  - The corner gizmo is in the bottom-left, about 78 px, and turns with the camera. Its
    arrows are clear, but its x/y/z letters are about 8 px. They are readable when zoomed
    in, and only just legible at native size. The y label touches the green arrow.
- **Ground.**
  - Dark theme: a high-contrast grid.
  - Light theme: the grid is very faint (`--grid-minor` `#c7ccd4` on `--canvas-bg`
    `#e9edf2`). The ground plane is hard to read, and only the shadows anchor the bodies.
- **Tip trace.** A thin violet line, visible in both themes.
- **Layout at 420 px.**
  - The summary cards wrap two per row, with nothing clipped. The check_eom card is
    green ("4.99e-08, passes (tolerance 1e-05)"). The energy drift card reads "4.12e-11,
    unforced".
  - The view is square and fills the width.
  - The controls wrap onto two rows (Pause/scrubber/time; Speed/Reset camera/Body axes;
    Traces/Theme), with no overlap.
- **Equations.** MathJax typeset every equation, with readable subscripts (J₁, Ω₁, m₁,
  R₁ᵀR₂, e₃).
  - The δW, R₁ and R₂ equation blocks fit the width.
  - The two longest Lagrangian addends (the m₂ kinetic term and the potential term) run
    past the right edge of their box. That box scrolls horizontally, by design. The
    screenshot hides scrollbars, so this is not a defect, but a reader at 420 px gets no
    visual hint that it scrolls.
- **Problems:** nothing overlaps, and nothing outside the equation scrollers is clipped.

## Agent walkthrough (only `SKILL.md` and `reference.md`)

**System:** two quadrotors carrying a slender rod. The configuration space is
R³ × S² × S² × S² × SO(3) × SO(3):
- rod centre `x` and rod direction `r ∈ S²`, with rotational kinetic energy ½ J_r ω_r·ω_r;
- cable directions `q1, q2 ∈ S²`, from each quadrotor to its rod end;
- quadrotor attitudes `R1, R2`;
- inputs: thrusts `f1, f2` along `R_i e₃`, and moments `M1, M2`.

No such system exists in the repository. The script is
a temporary directory; it is not committed.

| Step (skill) | Result |
|---|---|
| 1-2 define, L, δW | Done from the step list and the quadrotor example. δW of a thrust at `x_Qi = x ± a r − l_i q_i` was written as `Dot(x.delta() ± Cross(ξ_r, r)·a − Cross(ξ_i, q_i)·l_i, R_i e₃ f_i)` |
| 3 `compute_eom` | 0.07 s |
| 4 `check_eom` | **1.56e-7**, passes |
| 5 standard form, 7 codegen | First try; params `J1 J2 J_r a e3 g l1 l2 m_Q m_R` |
| 8 simulate (3 s, hover thrust, swung cable, both quads yawing) | Constraint drift 2.9e-13; energy change 2.25e-2 (forced) |
| 9 Scene + Report | Two quadrotors (box + 4 rotor cylinders each), the rod as a `Rod` between its ends, two cables, centre mass, end trace. Page: `data-gm-errors="0"`, MathJax done, 22 canvas colours, 143 kB. The summary reads "Energy change … the work of the inputs, not an integrator error". The rotors' red sectors and both quadrotors' axes show |

Everything worked with `python -W error`, first try. Where the skill was not enough:

1. **`layout.pack` keys for indexed variables.**
   - `reference.md` shows `omega_q=`, `Omega_R=`, `x_dot=`. That `S2("q1")` takes
     `omega_q1` and `SO3("R1")` takes `Omega_R1` has to be inferred from the double
     pendulum example.
   - Add one line: "velocity keys are `x_dot`, `omega_<name>`, `Omega_<name>`".
2. **Virtual work of a force at a point that depends on S².**
   - Step 2 says only "Forces enter as `Dot(variation, force)`". The rule I needed is
     δ(point) = `x.delta()` + Σ `Cross(ξ, q)·length` (with `δq = ξ × q` from the
     reference). I got it by pattern-matching the quadrotor example's sign.
   - One sentence would fix it: "the variation of a point is the variation of its
     position expression, with `δq = Cross(ξ, q)`, `δR = R Hat(η)`".
3. **No integrator check for a forced run.**
   - The skill correctly says a forced energy change is not drift. It then leaves the
     agent with no integrator check at all.
   - Suggest: "for a forced run, also simulate once with `u=None` and report that drift".
4. **The view is framed on the whole trajectory.** An open-loop run that drifts makes the
   bodies small at t = 0. This is documented as the framing rule; the skill could say to
   use short horizons or a controller for the report.
5. **Minor:** the `Report` signature in `reference.md` has no `*` (Recommended 4).

The skill misled nowhere. Its conventions (E = −EL + forcing, M = −inertia, S² tangent
plane, `print_eom` LaTeX versus `str` `*`) match the code.

## Gates

| Gate | Result |
|---|---|
| `python -m pytest tests/ -q -W error` | **3399 passed, 2 xfailed, 0 skipped, 114 s** (Chrome and network present, so the browser tests ran). The 2 xfails predate this branch |
| Viz tests on Python 3.10.20 (`-W error`, numpy/scipy/pytest) | First run: 99 passed, 1 failed (`test_page_under_artifact_csp`). Six reruns of the browser files: five 53/53, one with 1 probe failure, `reading 'children'` (see R1) |
| `ruff format --check` / `ruff check` on `geomech/ tests/ scripts/ examples/generated/ examples/inverted_rigid_pendulum_lqr.py examples/report_*.py` (ruff 0.17.4) | 161 files formatted; all checks passed |
| `mypy geomech --ignore-missing-imports` | no issues in 50 files |
| `mkdocs build --strict` | exit 0 (only the upstream mkdocs-2 notice) |
| `scripts/regenerate_generated.py --check` | 18 modules up to date |
| Examples, `python -W error` | double rigid pendulum: check_eom 4.99e-8, drift 4.1e-11, 181 kB, 5.7 s. Spherical: 2.32e-8, drift 7.7e-13, 91 kB, 2.7 s. Quadrotor + load (LQR): 3.98e-8, 287 kB, 6.4 s |
| CI file | Triggers: push and pull_request on `main` only. The lint job's ruff steps include `examples/report_*.py`. The 3.12 job runs all three report examples into a temp dir (scipy comes from `[dev]`). The browser tests skip without Chrome or without CDN reachability, but not without WebGL (R1) |

## `.gitignore` and the skill copies

- `git check-ignore` results:
  - ignored: `.claude/settings.local.json`, `.claude/settings.json`,
    `.claude/agents/x.md`, `.claude/commands/y.md`, `.claude/skills/private/SKILL.md`,
    `.claude/skills/x.md`;
  - not ignored: only `.claude/skills/geomech/*`.
- `git ls-files .claude` lists only `.claude/skills/geomech/SKILL.md` and
  `reference.md`.
- `skills/geomech/` and `.claude/skills/geomech/` are byte-identical (`cmp`,
  `diff -rq`).

## Docs and CHANGELOG against behaviour

- **The `traj.p` rule.** The docs, the skill (step 9), the CHANGELOG and the docstrings
  all say what the code does. `Scene.params_for` and `Report.__init__` check with
  `check_same_params` and use `traj.p`.
- **Energy label.** "Energy change" when `traj.u` is nonzero somewhere (`Report.forced`)
  matches the docs, the skill and `docs/simulation.md`. The `u is None` case is stated
  inaccurately (Recommended 2).
- **`Trajectory.u`.** Zeros for `u=None`, `None` when a callable fails: that is the code,
  and the docs say so.
- **MathJax restrictions.** The excluded packages and `processHtmlClass: "math"` in the
  template match `docs/visualization.md` and the CHANGELOG. The note on why `ui/safe` is
  not used is accurate.
- **Scale rules.** The defaults are 6 % and 1 % of the half-size, capped at 35 % and 6 %
  of the largest body. Positions keep six significant figures; arrow radius is 0.02 × the
  arrow length; there is no framing floor. All of this matches `report.py` and the
  template.
- **Sizes and times** in the docs table match my runs (181/91/287 kB).

## Pre-push PII scan

The remote `vkotaru/geomech` is **public**. Scanned:
- `git log a4961bb..HEAD -p`: 4 commits, all patches;
- the commit messages;
- the author and committer fields;
- the names and text chunks of the 11 committed PNGs and the 2 new ones.

| Check | Result |
|---|---|
| Author/committer | all `vkotaru <17078723+vkotaru@users.noreply.github.com>` (noreply) |
| User or real names, other emails | none |
| Home, scratch or tmp paths (`/home/`, `/tmp/`, scratchpad, `.local`) | none |
| Host, network or IP names | none |
| Names of the author's other local repositories and workspace | none |
| Local notes and memory files of the development setup | none added. `paper_context.md` appears only as unchanged context in the `.gitignore` hunk, and it is already in the public `.gitignore` |
| Tokens/keys | none (every "token" hit is a CSS theme token) |
| Scala references | none (every hit is `Scalar`) |
| Employer/affiliation | none |
| PNG files | descriptive names only; no text chunks with paths |

Hits that are not PII:
- Commit hashes from before the rebase (`3faf73c`, `49765b3`, `309d5ad`) appear in the
  review documents and in one test (Recommended 1).
- `Claude-Session:` trailers are in the commit messages, as the attribution rule
  requires.

## Open items (deferred), and whether `specs/follow-ups.md` records them

`specs/follow-ups.md` has no entry for this feature; none of the items below is recorded.

| Item | Source | In follow-ups.md |
|---|---|---|
| One-frame display lag during playback (unchanged) | fix log, F9 | no |
| MCP server | spec, non-goals ("next step") | no |
| `ui/safe` not loaded (a fourth URL would break the artifact CSP); TeX packages are excluded instead | fix log, F5 | no |
| The smoke test reads page hooks and misses errors thrown before the first inline script | implementation report | no |
| No WebGL: browser tests fail instead of skip; CDN dependence of the browser tests | implementation report; R1 here | no |
| No pixel-level check of screenshots (canvas colour count only) | implementation report | no |
| Forced-run integrator check; `u=None` label (Recommended 2-3) | this review | no |
| Rigid body hanging from an S² point cannot be drawn | spec, docs | documented in docs (a model limitation), not a follow-up |
| The arbitration rewrote the reviewer's strict-xfail test `test_skill_print_convention_matches_print_eom` | review, arbitration | recorded in the review |

That rewrite was a test from this branch whose assertion encoded the wrong fix, and the
change is documented; it does not count against "do not edit passing tests".

## Fix log

Applied on `feature/viz-skill` on top of `18b72fe`; nothing committed. No assertion of an
existing test was changed except the pinned base hash (Recommended 1); the
`assert "error" not in out` of the old fixture moved into the shared helper unchanged.

| Item | Fix |
|---|---|
| **R1** | New `tests/browser.py`, used by `test_viz.py`, `test_viz_adversarial.py` (the `browser` fixture) and `test_viz_fixes.py`. Flags: `--headless=new --disable-gpu --enable-unsafe-swiftshader --no-sandbox`. Opt-in with `GEOMECH_BROWSER_TESTS=1` (set only in the 3.12 CI job; the local command is in the module docstring and `docs/visualization.md`). Skips without Chrome, without network to the CDNs, when three.js or MathJax does not load (after one retry; MathJax loading is detected from `MathJax.version`), when Chrome produces no DOM, and when the page reports no WebGL. |
| R1, ready hook | The page now sets `data-gm-viewer` on the root (`ready`, `no-three`, `no-webgl`, `none`). Every scene probe (`_POST`, `_PUMP`, `_MARKERS`, `_GIZMO`, `_TINTS`, `_PRESS`) is wrapped in `when_ready`, which waits for the hook and reports `unavailable` instead of reading an undefined scene. The `children` TypeError was the probe reading `window.__gmScene` after the main script had returned early (WebGL context creation failed or three.js did not load), so no scene was ever built; it now skips with the reason. A plain "WebGL is not available" substring check would always match (the text is in the inline script), hence the attribute. |
| R1, CSP test | `test_page_under_artifact_csp` goes through the same helper: one retry when a CDN script did not load or MathJax loaded but did not finish typesetting without errors, then a skip naming the reason. Verified by hand: a page with an unreachable three.js URL skips ("three.js did not load from its CDN (2 attempts)"); with `--disable-software-rasterizer` both the smoke load and a probe skip ("no WebGL"). |
| R1, 10 runs | `GEOMECH_BROWSER_TESTS=1 pytest tests/test_viz.py tests/test_viz_adversarial.py tests/test_viz_fixes.py -W error`, 10 consecutive runs: **10/10**, each 91 passed, 0 skipped, 0 failed (49–53 s). An earlier attempt hit `Disk quota exceeded` on this machine's `/tmp` (filled by an unrelated process), which made Chrome's GPU process abort; those runs skipped (no DOM) rather than failed. The 10/10 loop ran with `TMPDIR` on the home disk. |
| Rec. 1 | `BASE = "a4961bb"` (the merge base on `main`); the test first checks `git cat-file -e a4961bb^{commit}` and skips with "base commit a4961bb is not in this clone (a shallow checkout?)". |
| Rec. 2 | `report.py`: "Energy drift … unforced" only when `traj.u` is recorded and all zero; `traj.u is None` gives "Energy change … inputs not recorded"; nonzero inputs keep "Energy change … the work of the inputs". `Report.forced()` docstring says which case is which. |
| Rec. 3 | `simulate`: RKMK4 records each stored step's first-stage input (the singularity check's `u(t0, x0)` is reused for step 0), plus one call at the final state; the controller is never called twice at one state (`4·n_steps + 1` calls). `_stored_inputs` is gone. An exception from a callable `u` propagates as soon as it is raised, unchanged in type, with a note `simulate: raised by u(t, state) at t = …` (`add_note`, or `__notes__` on 3.10). With `solve_ivp` the stored states are interpolated, so `u` is evaluated once at each after the run (documented; errors propagate the same way). |
| Rec. 4 | `reference.md` shows `checks=True, *, units=None, source=None, derive_time=None, fps=60, description=""`. |
| Walkthrough 1 | New `StateLayout.slot_names()`; SKILL step 8, `reference.md` and `docs/codegen.md` state the keys `x_dot`, `omega_<name>`, `Omega_<name>` (`omega_q1`, `Omega_R2`). |
| Walkthrough 2 | SKILL step 2 and `reference.md`: δ(point) term by term, `x → x.delta()`, `q → Cross(ξ, q)`, `R → R * Cross(η, ·)`, with the S² and body-point examples; both checked with `check_eom` (2.4e-8, 2.8e-8). |
| Walkthrough 3 | SKILL step 8 and Presentation: for a forced run, also simulate with `u=None` and report that drift; the quadrotor worked example does it (6.8e-14). |
| Walkthrough 4 | Documented (the smaller change): SKILL step 9, `reference.md`, `docs/visualization.md` ("Framing"): the view is framed on the whole trajectory; use a shorter horizon or a controller. A `frame=` option is in `specs/follow-ups.md`. |
| Visual polish | Light grid `#7b8492`/`#aab2bd` (was `#9aa2ad`/`#c7ccd4`). Corner gizmo at least 96 px (24 % of the view, at most 140), labels 0.6 in a 3.3-wide view drawn at 52/64 of the sprite: about 14 px at 420 px, clear of the arrow tips. Equation boxes show a scroll shadow on each side with more to scroll to (`--scroll-hint`, CSS only) and a thin scrollbar. The body-e₃/inertial-z overlap at a pivot is documented in `docs/visualization.md` (an arrow inside the box would be invisible, and a triad off the body's centre would misplace the frame); listed in follow-ups. Screenshot: `viz-final-fix-light-420.png` (darker grid, readable gizmo letters, shadow on the right of the Lagrangian box). |
| Rec. 7 | `specs/follow-ups.md` "Visualizer and skill", items 30–39: MCP server, playback lag, `ui/safe`, smoke-test blind spot, CDN/WebGL dependence, no pixel check, trajectory framing, pivot arrow overlap, `solve_ivp` post-run `u` evaluation, rigid body on an S² point. |
| Rec. 8 | CHANGELOG: the "fixes after the adversarial review" bullet is folded into the `Scene`, `Report`, `latex_names`, skill and `simulate` bullets; `slot_names` added. |
| New tests | `tests/test_viz_final.py` (8): controller call count and recorded values (stateful controller, `output_every` 1 and 7), `KeyError` at t = 0.025 propagates at once with the time note, `IndexError` with `solve_ivp`, `solve_ivp` records `u`, unrecorded inputs labelled "Energy change … inputs not recorded", `slot_names`, ready-hook states in the page. |

### Gates after the fixes

| Gate | Result |
|---|---|
| `pytest tests/ -q -W error` | 3381 passed, 26 skipped (the opt-in browser tests), 2 xfailed |
| `GEOMECH_BROWSER_TESTS=1 pytest tests/ -q -W error` | 3407 passed, 2 xfailed, 0 skipped |
| Browser files ×10 | 10/10 (91 passed each) |
| `ruff format --check` / `ruff check` (CI paths) | 163 files formatted; all checks passed |
| `mypy geomech --ignore-missing-imports` | no issues in 50 files |
| `scripts/regenerate_generated.py --check` | 18 modules up to date |
| `mkdocs build --strict` | exit 0 |
| Examples, `python -W error` | double pendulum 4.99e-8, drift 4.12e-11, 182 kB; spherical 2.32e-8, drift 7.74e-13, 92 kB; quadrotor + load 3.98e-8, 289 kB |
| Skill copies | `skills/geomech` and `.claude/skills/geomech` byte-identical |

Not run here: the suite on Python 3.10/3.11 (CI covers them; the 3.10 path of the
exception note is `__notes__`).
