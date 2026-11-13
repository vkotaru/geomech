# Implementation report: variation-based linearization

Spec: `specs/vbl.md`. Branch `feature/vbl`. Nothing committed; everything is in the
working tree.

## What was built

| File | Content |
|---|---|
| `geomech/dynamics/vbl.py` | `linearize(eom, variables, inputs=None)` and the frozen dataclass `LinearizedSystem` (`rows`, `states`, `dstates`, `inputs`, `M`, `A`, `B`, `block_shape`, `__str__`). Module docstring carries the perturbation table, kinematic rows, sign convention and completeness rule. |
| `geomech/dynamics/__init__.py`, `geomech/__init__.py` | Export `linearize`, `LinearizedSystem` next to `to_standard_form`. |
| `tests/numeric_helpers.py` | `hat`, `vee`, `rand_rot`, `expm_so3` (Rodrigues), `Evaluator` (promoted from the review script, with `env` injection, `TimeDerivative` chains over leaves, `Variation` values required from `env`, zero/identity constants resolved by type), `sample_state`, `sample_direction`, `perturbed_state`, `linearized_rhs`, `finite_difference_check`. The review script `specs/reviews/scripts/double_rigid_pendulum_numeric_check.py` is untouched. |
| `tests/vbl_systems.py` | The nine systems as `System(L, dW, variables, inputs, symbols)`. |
| `tests/test_vbl.py` | 100 tests: systems 1–9, cross-cutting block properties, completeness errors, input validation, numeric-helper self-tests. |
| `docs/vbl.md`, `mkdocs.yml` | User-facing page (perturbation variables, rows, sign convention, worked free-rigid-body example, how to read the blocks) and nav entry. |
| `CHANGELOG.md` | Entry under Unreleased / Added. |
| `README.md` | First line made accurate (R, R³, S², SO(3) and products; SE(3) as R³ × SO(3)); one pipeline step added for `to_standard_form` / `linearize`. |

Final counts: `python -m pytest tests/ -q` → **1234 passed, 5 xfailed** (baseline 1134 + 100
new; the 5 strict xfails are the pre-existing ones). `ruff format --check` and
`ruff check` clean. No existing test file was edited.

## Algorithm

For each dynamic row, `δE = expand(full_simplify(E.delta()))`. The additive terms are
assigned to exactly one perturbation symbol by testing `_has_algebraic(term, sym)` in the
priority order *d/dt of velocity perturbations → velocity perturbations → configuration
perturbations → input variations*, so that a term in `d/dt(δΩ)` is never tested against
`δΩ` (which would raise inside `_has_algebraic`). Each bucket is combined and passed to
`extract_linear_coeff`, then `full_simplify`; the coefficient is checked to be free of every
perturbation symbol (including `ξ̇`, `η̇`). Mass blocks are negated before simplification
(`M = −coeff`); a scalar-row × vector-column coefficient is wrapped in `Transpose` to make
the 1×3 block. Kinematic rows are written directly from the formulas in the spec.

`apply_manifold_rules` is not called anywhere in the VBL path.

## Hand derivations (recorded in the test docstrings)

1. Point mass: E = −m ẍ + F → M[δẍ] = m I, B[δF] = I, A = 0. Exact strings.
2. Spring-mass: E = −m s̈ − k s + f → M = m, A[δs] = −k, B = 1. Exact strings.
3. Free rigid body (Lee et al. 6.16): δE = −J δΩ̇ + (Hat(JΩ) − Hat(Ω) J) δΩ + δM.
   Exact strings (`J`, `(Hat(JΩ)+Hat(Ω)J(-1))`, `I`) and numeric equality.
4. Rigid pendulum about the pivot (Lee et al. 6.8): δ(Rᵀe₃) = Hat(Rᵀe₃) η gives
   A[η] = −m g Hat(ρ) Hat(Rᵀe₃); exact string `Hat(ρ)Hat((R)'e3)(-1)mg` and numeric equality.
5. Spherical pendulum (Lee et al. 5.14): δ(q × e₃) = Hat(e₃) Hat(q) ξ gives
   M = m l² I, A[ξ] = −m g l Hat(e₃) Hat(q), B = I; kinematic row
   `q(q)'Hat(ω)`, `(I+q(q)'(-1))`. Exact strings; the hanging-equilibrium check
   (A[ξ] ξ = −m g l ξ at q = −e₃) and a numeric check of the kinematic row against
   q qᵀ ω̂ ξ + (I − q qᵀ)(ξ × ω + w) with d/dt(ξ·q) = 0.
6. Quadrotor: δ(f R e₃) = R e₃ δf − f R Hat(e₃) η → A[η] = `RHat(e3)(-1)f`, B[δf] = `Re3`
   (3×1), M = m I; rotational row as system 3. Exact strings and numeric equality.
7. Double rigid pendulum: exact kinematic rows; the mass blocks are the negatives of the
   standard-form blocks pinned in `TestDoubleRigidPendulumMassMatrixExact`
   (`(J1+Hat(ρ₁)Hat(ρ₁)m1(-1)+Hat(l₁)Hat(l₁)m2(-1))`, `Hat(l₁)(R1)'R2Hat(ρ₂)m2(-1)`, …);
   finite differences for everything else.
8. Elastic spherical pendulum: block shapes (1×3, 3×1) and the 1×3 block
   A[δl][δω] = (2 m l ω)ᵀ from δ(m l ω·ω); finite differences.
9. Quadrotor with cable-suspended load: `compute_eom` derives it without any engine change,
   so no `xfail` was needed. B[δf] = `Re3` (load row) and `Cross(q,Re3)l(-1)` (cable row)
   are pinned; finite differences on all three dynamic rows.

Finite-difference check (every system): 5 random base states on the manifold, random
directions tangent to the constraints (S²: q(ε) = exp(ε ξ̂) q, ω(ε) = exp(ε ξ̂)(ω + ε w),
ω̇(ε) = exp(ε ξ̂)(ω̇ + ε a) with ξ, w, a ⊥ q, inducing δω = ξ × ω + w; SO(3):
R exp(ε η̂), Ω + ε δΩ; R, R³, inputs: linear), ε = 1e-5, central differences, error
measured relative to max(1, |rhs|) with tolerance 1e-6. Observed worst errors are
1e-10 to 1e-12. `test_finite_difference_check_detects_a_wrong_block` flips the sign of one
block and asserts that the check fails, so the check is evidence and not a tautology.

## Engine changes

**None.** No file under `geomech/core/` was modified.

The spec anticipated a `simplify.py` rule for the `Hat((R)'e3)(-1)(-1)` chain. The probe
that produced that string printed the *raw* output of `extract_linear_coeff`; the existing
`pull` + `_eliminate` pass already folds `SMMul(SMMul(X, -1), -1)` to `X`, and `linearize`
always runs `full_simplify` on each block. `test_no_unfolded_numeric_factors` asserts that no
block of any of the nine systems prints a `(-1)(-1)` or `(-1.0)(-1` chain.

## Design decisions not fixed by the spec

- `LinearizedSystem` carries two extra fields, `row_types` and `col_types`
  (`dict[str, ExprType]`), so `block_shape` works from labels alone; `block_shape` accepts
  the column symbol or its `str`.
- The 1×3 zero block is `Transpose(ZeroVector)` (`is_zero` propagates through `Transpose`);
  the other zero blocks are `Zero`, `ZeroVector`, `ZeroMatrix`.
- In dynamic rows the `M` columns belonging to configuration dstates (`ξ̇`, `η̇`, and
  `d/dt(δs)` / `d/dt(δx)`, whose printed names coincide with the velocity states) are always
  zero: a term in `δṡ` is a state term and goes to `A`. The finite-difference helper still
  supplies values for these columns (from the kinematic rows) so a non-zero entry there
  would be caught.
- Validation: `variables` go through `compute_eom`'s `_check_variables`; the `eom` keys must
  equal the variation-vector names of `variables` (`ValueError` otherwise) and each value
  must be a `(variation_vector, Expr)` pair; inputs must be non-constant `Scalar`/`Vector`
  leaves that are not configuration variables; a non-`SO3` matrix configuration variable
  raises `TypeError` (the spec's perturbation table has no entry for it).
- A dynamic equation that is identically zero gives an all-zero row instead of raising.
- Terms with a perturbation symbol only inside a calculus node raise `AlgebraicError` with
  the row label, the offending term and the symbol (wrapping the `_has_algebraic` error).

## Things not done / not verifiable here

- `mkdocs` is not installed in the venv, so the docs build with the new page and nav entry
  was not run; the page uses the same arithmatex `\(...\)` / `$$` syntax as the existing
  pages.
- Nothing from the spec was impossible; no part of the spec was skipped.
