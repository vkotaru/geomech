# Phase 1, item 1.6 — numeric oracle in the package (`geomech.numeric`)

Branch `fix/phase1-equations`, uncommitted. No file under `geomech/core/` or
`geomech/dynamics/` was touched; no existing test function was edited.

## What was built

| File | Change |
|---|---|
| `geomech/numeric.py` (new) | The evaluator and finite-difference helpers from `tests/numeric_helpers.py`, with docstrings, plus `evaluate`, `check_eom`, `assert_eom`. |
| `geomech/__init__.py` | Exports `Evaluator`, `evaluate`, `check_eom`, `assert_eom`, `finite_difference_check`. `numeric` imports only `core.base` and `core.operations` at runtime (`SystemVariables` under `TYPE_CHECKING`), so there is no import cycle and the package imports cleanly from any cwd. |
| `tests/numeric_helpers.py` | Thin re-export of `geomech.numeric` with an explicit `__all__`. `tests/test_vbl.py`, `tests/test_vbl_adversarial.py`, the audit harness and `engine_fuzz.py` import it unchanged. |
| `tests/test_numeric.py` (new) | 26 tests: SO(3) primitives, `evaluate` on known trees, oracle ≈ 0 on four pinned systems, prescribed-parameter and no-forcing cases, and the *detection* tests (below). |
| `tests/test_eom_oracle.py` (new) | 21 systems, each asserting `check_eom(...) < 1e-5`. |
| `docs/numeric-verification.md` (new), `mkdocs.yml` | "Verifying an equation numerically": key convention, 15-line example, what is and is not verified, how to use it on a new system. `mkdocs build --strict` passes. |
| `CHANGELOG.md` | Entry under Unreleased / Added. |

## API

```python
evaluate(expr, env: dict[str, np.ndarray | float] | None = None, *, seed=0) -> np.ndarray | float
```
Leaves are keyed by `str(leaf)`: `"m"`, `"x"`, `str(TimeDerivative(x))` for `ẋ`,
`str(q.get_tangent_vector())` for `ω`, `str(R.get_tangent_vector())` for `Ω`. Unlisted
leaves are drawn once from `default_rng(seed)`; `env` is copied, not mutated; a `Variation`
must be given. Zero/identity constants are resolved by type.

```python
check_eom(lagrangian, inf_work, variables, eom, *, inputs=(), params=(), trials=5, seed=0, h=1e-5) -> float
assert_eom(..., rtol=1e-5) -> float      # AssertionError naming the variation vector and the error
```
Port of `check_system` from `specs/reviews/scripts/euler_lagrange_check.py`, restructured
into short functions (`_Coordinates`, `_constant_env`, `_forcing`, `_euler_lagrange`,
`_relative_error`, `_eom_errors`). Sign convention documented in the docstring and the docs
page: `E = −(Euler–Lagrange) + forcing`. Residuals: R/R³ `g = d/dt ∂L/∂v − ∂L/∂q`; S² via
`L̃(q, q̇) = L(q, q × q̇)`, `E_true = −q × g + f`, compared on the tangent plane; SO(3) via
`g = d/dt ∂L/∂Ω − (∂L/∂Ω) × Ω − vee(A − Aᵀ)`, `A = Rᵀ ∂L/∂R`. `params` are prescribed
time-varying symbols (the harness's `extra_params`); `inf_work` may be `None` or the zero
constant. The error is `max |P (E_lib − E_true)| / max(1, |P E_true|)`.

## One change to the evaluator's behaviour

`Evaluator.leaf` now honours `UnitNorm` (unit vector), `OrthogonalMatrix` (rotation) and
`SkewSymmetricMatrix` when drawing a leaf, in addition to `SymmetricMatrix`, `S2`, `SO3`
and `Scalar.value` as before. Without it the R³ × R prismatic-slide system
(`tests/test_scalar_adversarial.py::TestR3xRCoupledMassMatrix`, `e1` declared `UnitNorm`)
scored 0.96: the engine correctly uses `e1·e1 = 1`, the evaluator drew a non-unit `e1`.
No existing numeric test uses such a leaf, so the random streams of the existing suite are
unchanged (verified: `test_vbl.py` and `test_vbl_adversarial.py` pass as before).

## Evidence that the oracle is not a tautology

- `test_check_eom_detects_a_flipped_term` flips the sign of *every* additive term of every
  equation of the spherical pendulum, free rigid body, double rigid pendulum and quadrotor
  with load, and asserts the error exceeds 1e-2 (observed: 0.1 to 1.4). The one term skipped
  is `ω × (q × ẋ_L) = q (ω·ẋ_L)` in the load equation, which is normal to `q` and therefore
  invisible on the tangent plane by design; the test detects such terms numerically
  (projection at random on-manifold states), not by index, so it survives engine changes
  that reorder terms.
- `test_check_eom_detects_a_wrong_coriolis_term`: hand-written `−JΩ̇ − Ω×JΩ + M` passes,
  `−JΩ̇ + Ω×JΩ + M` scores > 0.1.
- `test_check_eom_with_a_prescribed_parameter`: pendulum with cable length `l(t)` passes
  with `params=[l]` (2e-7) and fails without it (1.8), which is the documented contract.

## Systems pinned in `tests/test_eom_oracle.py` (all pass, worst error per system)

| System | Manifold | Error |
|---|---|---|
| point_mass_r3, point_mass_gravity, free_particle (δW = 0) | R³ | ≤ 7e-7 |
| spring_mass_r, spring_mass_damper, double_spring_mass, varying_inertia, squared_inertia | R, R² | ≤ 9e-7 |
| free_rigid_body, rigid_pendulum_pivot, rigid_pendulum_com | SO(3) | ≤ 8e-7 |
| spherical_pendulum_s2, spherical_pendulum_qdot | S² | ≤ 1.6e-6 |
| double_rigid_pendulum | SO(3)² | 5.8e-6 |
| quadrotor, point_mass_and_rigid_body | R³ × SO(3) | ≤ 2.9e-6 |
| elastic_spherical_pendulum, elastic_spherical_pendulum_qdot | R × S² | ≤ 1.7e-6 |
| prismatic_slide | R³ × R | 5.8e-7 |
| sliding_mass_on_rotating_body | SO(3) × R | 3.5e-6 |
| quadrotor_with_load | R³ × S² × SO(3) | 1.6e-6 |

Every system reconstructed from `tests/test_manifold_pipelines.py`,
`tests/test_standard_form.py` and `tests/test_scalar_adversarial.py` is included; none
had to be excluded, and no engine bug surfaced on these systems (as the plan expected: the
known bugs are masked for them).

## Counts

- `python -m pytest tests/ -q`: **1375 passed, 5 xfailed** (baseline 1322 + 5; +53 new tests:
  26 in `test_numeric.py`, 21 parametrized cases + 0 others in `test_eom_oracle.py`, and 6
  parametrized cases of the two oracle tests in `test_numeric.py`). The new files run in
  under 2 s.
- `ruff format --check` and `ruff check` on `geomech/ tests/`: clean.
- `mkdocs build --strict`: passes (the only output is a pre-existing plugin deprecation
  banner).

## Not done / notes for the next task

- The plan says every pinned pipeline test "gains a `check_eom` assertion". Editing the
  existing test functions was out of scope, so the gate lives in `tests/test_eom_oracle.py`,
  which re-derives the same Lagrangians. When the engine fixes (1.1–1.5) land, add the
  audit's six wrong rows and four crash rows to `SYSTEMS` there.
- `check_eom` does not verify `to_standard_form` or `linearize`; the latter is covered by
  `finite_difference_check`, the former is not covered numerically yet (the audit did it
  ad hoc). A `check_standard_form` would be a 20-line addition on top of `evaluate`.
- `evaluate` ignores a `value=` on `Vector`/`Matrix` leaves (only `Scalar.value` is used),
  as the evaluator always did. Honouring it would be a one-line change but alters the
  random stream of existing tests that construct such leaves, so it was left alone.
