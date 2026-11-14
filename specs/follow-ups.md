# Open follow-ups after the 2026-11-13 fix plan

Source: the closure table in `specs/reviews/phase5b-final.md` (Part B). The fix plan
(`specs/fix-plan.md`) is complete: 167 items done, 16 deferred with a record, and the 23
below open (the 24th, an input named like a constant leaf, is fixed; so is the later
finding that `to_standard_form` accepted a `Constant`, duplicate or configuration-variable
input that `linearize` rejects, now one shared `check_inputs`). None of them changes a
derived equation.

## Behaviour

1. `TimeIntegral(c).is_constant` is `True` for a constant `c` (engine audit F17).
2. `Evaluator`/`evaluate` ignore `value=` on `Vector`/`Matrix` leaves (draws random values).

## Features

3. `check_standard_form`: numeric check of `to_standard_form` (`Σ M a + f + G u = E`).
4. `check_rewrite(before, after)`: rule-level numeric check, optionally on every rule
   firing under `GEOMECH_DEBUG`.
5. Promote the fuzzer from `specs/reviews/scripts/engine_fuzz.py` to `scripts/` or
   `geomech.numeric`, since `docs/pipeline.md` and `docs/troubleshooting.md` point users to it.

## Documentation

6. `StandardFormEquation` docstring: for a scalar row with a vector acceleration, `M[a]`
   means `M[a]ᵀ a`.
7. README pointer to `check_eom` / `docs/numeric-verification.md`.
8. The accompanying paper's S² variation paragraph could add "and from ω × ξ ∥ q" to the
   justification of `qqᵀξ̇ = ω × ξ`.

## Cosmetic code

9. `(d/dt J)ᵀ` of a symmetric `J` is not simplified to `d/dt J`.
10. Magic numbers: `full_simplify(max_iter=10)`, `len(name) <= 10`, `gap=3`, `[-2000:]`
    in `geomech/utils/printing.py`.
11. The "no left coefficient" error says "3-tensor" also for a matrix query.
12. Export `TermResolution` and `DEFAULT_H` from `geomech`.
13. The `inf_work` time-derivative hint is phrased for R³ (`Dot(x.delta(), f)`) on S² too.
14. Mixed error-message prefixes (`attr must be a list…` has no `attr:` prefix).
15. `_is_row`/`_q` one-line shims in `geomech/core/operations/multiplication.py`.
16. `display_latex` keeps a sum's outer parentheses.
17. `__rtruediv__` quotes its operand with `!r`, unshortened.
18. `print_tree` uses `├─` for the last child (no `└─`).

## Test hygiene

19. `HALF` is defined in seven test modules; some system builders are duplicated.
20. Five tests are still named `test_pinnedequation`.
21. `tests/test_validation.py::TestStrayManifoldVariation` covers more than its name says.
22. Two test docstrings (`TestRotatingSphere` "using ω ⊥ q"; the SO(3)×SO(3) orientation
    note) were not updated after phase 1.
23. `row_spring` and `row_quadratic_form` (`tests/test_row_forms.py`) are checked
    numerically only; pin their strings.

## Code generation and simulation

From the final review of the input-validation, code-generation and simulation branch
(`specs/reviews/codegen-simulation-final.md`, "Open items"). Its R1 to R5 and the state
helper (`StateLayout.pack`/`unpack`) are fixed, and the codegen energy now imports the
public `perturbation_variables`/`partition_terms` of `geomech.dynamics.vbl`. Open:

24. `energy_function`/`momenta` do not support prescribed parameters
    (`compute_eom(..., params=[l])`): the variation of `L` carries `δl`, which is no
    velocity or configuration perturbation, and the term partition rejects it, as in
    `linearize`.
25. `TimeDerivative(Transpose(x))` is rejected by the code generator as "a TimeDerivative
    of the non-leaf expression". Conforming (no engine output holds it), but
    `d/dt(xᵀ) = ẋᵀ` could be generated.
26. Simulation stepper overhead: 35–55 µs per RKMK4 stage outside the generated functions
    (preallocated stage arrays, one fused `mass`/`bias`/`input_matrix` call per stage).
    `Trajectory.energy` calls `energy` once per stored state and is not vectorised.
27. `lqr` on an S² block: an anisotropic `Q` weights directions of `tangent_basis(q_eq)`,
    which turns by 90° under a tiny change of `q_eq` where the two smallest components of
    `q_eq` tie (near `±e₃`, for example). Documented in `docs/simulation.md` ("Anisotropic
    weights on an S² block depend on the basis"); a `bases=` argument to `lqr` taking a
    user basis per S² variable would remove the dependence.
28. The test suite grew from about 17 s to about 60 s (the 10 s RKMK4 energy runs of
    `tests/test_simulation_conservation.py`, slowest 5.3 s). Option: a registered `slow`
    marker on those rows, deselected locally with `-m "not slow"` and run in CI.
29. `StateLayout.pack` covers `state` only; `u` is still stacked by hand in the layout's
    input order (`layout.input_slices`). A `pack_inputs(**inputs)` would complete it.
