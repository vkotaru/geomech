# Open follow-ups after the 2026-11-13 fix plan

Source: the closure table in `specs/reviews/phase5b-final.md` (Part B). The fix plan
(`specs/fix-plan.md`) is complete: 167 items done, 16 deferred with a record, and the 24
below open. None of them changes a derived equation; the first one can give a silently
incomplete linearization and should be fixed first.

## Behaviour

1. **An input with the name of a constant leaf is accepted silently.**
   `linearize(eom, V, [f, Scalar("g")])` on `… − m g e₃·x` (with `g` a `Constant` in L)
   returns no `B[δg]` column and no error; `to_standard_form` gives `G[g] = −m e₃`.
   Reject an input whose printed name equals a constant leaf in the equations.
2. `TimeIntegral(c).is_constant` is `True` for a constant `c` (engine audit F17).
3. `Evaluator`/`evaluate` ignore `value=` on `Vector`/`Matrix` leaves (draws random values).

## Features

4. `check_standard_form`: numeric check of `to_standard_form` (`Σ M a + f + G u = E`).
5. `check_rewrite(before, after)`: rule-level numeric check, optionally on every rule
   firing under `GEOMECH_DEBUG`.
6. Promote the fuzzer from `specs/reviews/scripts/engine_fuzz.py` to `scripts/` or
   `geomech.numeric`, since `docs/pipeline.md` and `docs/troubleshooting.md` point users to it.

## Documentation

7. `StandardFormEquation` docstring: for a scalar row with a vector acceleration, `M[a]`
   means `M[a]ᵀ a`.
8. README pointer to `check_eom` / `docs/numeric-verification.md`.
9. The accompanying paper's S² variation paragraph could add "and from ω × ξ ∥ q" to the
   justification of `qqᵀξ̇ = ω × ξ`.

## Cosmetic code

10. `(d/dt J)ᵀ` of a symmetric `J` is not simplified to `d/dt J`.
11. Magic numbers: `full_simplify(max_iter=10)`, `len(name) <= 10`, `gap=3`, `[-2000:]`
    in `geomech/utils/printing.py`.
12. The "no left coefficient" error says "3-tensor" also for a matrix query.
13. Export `TermResolution` and `DEFAULT_H` from `geomech`.
14. The `inf_work` time-derivative hint is phrased for R³ (`Dot(x.delta(), f)`) on S² too.
15. Mixed error-message prefixes (`attr must be a list…` has no `attr:` prefix).
16. `_is_row`/`_q` one-line shims in `geomech/core/operations/multiplication.py`.
17. `display_latex` keeps a sum's outer parentheses.
18. `__rtruediv__` quotes its operand with `!r`, unshortened.
19. `print_tree` uses `├─` for the last child (no `└─`).

## Test hygiene

20. `HALF` is defined in seven test modules; some system builders are duplicated.
21. Five tests are still named `test_pinnedequation`.
22. `tests/test_validation.py::TestStrayManifoldVariation` covers more than its name says.
23. Two test docstrings (`TestRotatingSphere` "using ω ⊥ q"; the SO(3)×SO(3) orientation
    note) were not updated after phase 1.
24. `row_spring` and `row_quadratic_form` (`tests/test_row_forms.py`) are checked
    numerically only; pin their strings.
