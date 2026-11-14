# Phase 2 implementation report: error instead of silence

Branch `fix/phase2-validation`, stacked on `fix/phase1-equations` at `865699f`. Nothing is
committed. Spec: `specs/fix-plan.md`, phase 2 table (2.1–2.12) plus the error-export
follow-up from `specs/reviews/phase1-final.md`. Round 2 applied the coordinator's four
decisions (section "Decisions applied").

Baseline on this machine: 1552 passed, 5 xfailed. Final: **1630 passed, 0 failed,
5 xfailed** (82 new tests). `ruff format --check` and `ruff check` pass; `mkdocs build
--strict` passes; `engine_fuzz.py --quick` reports 0 wrong-value cases and no new
exception class; all 6 example scripts run.

## Per item

| # | Where (file:line) | Exception and message | Tests (`tests/test_validation.py` unless noted) |
|---|---|---|---|
| 2.1 | `geomech/core/base/types.py:15-28` (`VALID_ATTRS`), `:45` (`from_attr_list(attrs, name=None)`), `:77-84` (`"TangentVector"` accepted as a marker; `case _` raises); `geomech/core/base/expressions.py:161,269,421` pass the leaf name | `ValueError: unknown attribute 'constant' for symbol 'm'; valid attributes are ['Constant', 'Identity', 'Manifold', 'Ones', 'OrthogonalMatrix', 'SkewSymmetricMatrix', 'SymmetricMatrix', 'TangentVector', 'UnitNorm', 'Zero']` | `TestUnknownAttr` (6) |
| 2.2 | `geomech/dynamics/eom.py:120-162` (`_check_params`), `:165-199` (`_stray_variation_symbols(expr, allowed, param_names)`), `:205-235` (`compute_eom(..., params=())`, docstring), `:281-293` (check after extraction over δS and every extracted equation) | `ValueError: compute_eom: δS still contains the variation symbols '\delta{y}' of the symbols 'y', which are neither entries of SystemVariables nor declared in params; for each symbol, list it in SystemVariables, declare it in params=..., or mark it 'Constant'` (manifold symbols print as `'\xi_{q}', '\delta{\omega_{q}}' of the symbols 'q'`). `params` validation: `TypeError: compute_eom: params entry '\omega_{q}' is a tangent or variation vector, not a prescribed parameter`; `TypeError: … params entry 'p' is a manifold point, not a prescribed parameter; list it in SystemVariables`; `TypeError: … params entry <repr> must be a plain Scalar, Vector or Matrix leaf, got Mul`; `ValueError: … params entry 'g' is constant; a prescribed parameter varies in time (drop it from params or remove the 'Constant' attribute)`; `ValueError: … params entry 'x' is also a configuration variable (SystemVariables.vectors); a symbol is either prescribed or a degree of freedom`; `ValueError: … params entry 'l' is listed more than once` | `TestStrayManifoldVariation` (10: unlisted S2/SO3 point, `S2("q")` listed as `Vector("q")`, unlisted `y` in L (dbg-A24), `δp·F` in δW (dyn-F5), unflagged `e3` (dyn-F5), `l(t)` with and without `params`, force in δW unaffected, listed points pass), `TestParamsValidation` (8, incl. `to_standard_form` on a params equation and `linearize` keeping its `AlgebraicError` on the `δl` term) |
| 2.3 | `geomech/dynamics/eom.py:23-103` (`_check_variables`); `geomech/dynamics/vbl.py:360-373` (`_check_eom` compares `len(eom)` with the expected list) | tangent/variation vector: `TypeError: SystemVariables.vectors entry '\omega_{q}' is a tangent or variation vector, not a configuration variable; list the manifold point 'q' instead`; calculus node: `TypeError: SystemVariables.vectors entry '\frac{d}{dt}(x)' is a TimeDerivative node, not a configuration variable; list the underlying variable 'x' instead`; non-SO3: `TypeError: SystemVariables.matrices entry 'J' is a general Matrix; only SO3 rotation matrices are supported as matrix configuration variables`; duplicate: `ValueError: SystemVariables: the name 'x' is listed more than once (scalars and vectors); configuration variables must be distinct symbols`; name/variation clash: `ValueError: SystemVariables.matrices entry 'R' has the variation vector '\eta_{R}', which is also the name of the entry listed in vectors; …` and `… entry '\xi_{q}' has the printed name of the variation vector of 'q'; …`; `_check_eom`: `ValueError: linearize: eom keys […] do not match the variation vectors […] of the given SystemVariables` | `TestSystemVariablesValidation` (12) |
| 2.4 | verified, no code: `expressions.py:195` (Scalar), `:308` (Vector), `:338` (S2), `:455` (Matrix), `:463-466` (SkewSymmMatrix), `:474` (SO3) all copy `attr`; `TS2`/`TSO3` call `super().__init__(s)` without an attr and append to their own list | — | `TestAttrListCopied` (6, parametrised) |
| 2.5 | `geomech/core/transformations/vector_rules.py:43,45,111` (`a.S2 == b` instead of `is`; `==` on leaves is type + printed name) | no exception; `Dot(ω_q, q₂) → 0`, `q₂ × (ω_q × q₂) → ω_q` for a second `S2("q")` | `TestTangentParentByName` (4, incl. a pipeline test that the two-object equation equals the one-object equation) |
| 2.6 | `geomech/dynamics/eom.py:105-117` (`_check_scalar_input`), `:236-239` (`None`/`0`/`0.0` → `Zero`, then both arguments checked) | `TypeError: compute_eom: inf_work must be a scalar expression, got a VECTOR expression 'f'`; `TypeError: compute_eom: inf_work must be a scalar expression, got str '0'`; same text with `lagrangian` | `TestInfWork` (4) |
| 2.7 | `geomech/core/transformations/manifold_rules.py:143-147` | `TypeError: apply_manifold_rules: cannot rebuild a _Opaque node ('Opaque(\delta{\omega_{q}})') with substituted children; add the class to manifold_rules._rebuild` | `TestRebuildUnknownNode` (3) |
| 2.8 | `geomech/core/operations/calculus.py:9-21` (`_reject_manifold_leaf`), `:38` (`Variation.__init__`), `:87` (`TimeDerivative.__init__`); `geomech/core/operations/geometry.py:123-127` (`Transpose()` without operand) | `ValueError: TimeDerivative(q) is not represented in the tree for the manifold point 'q'; use q.t_diff(), which returns the kinematic expression ω × q` (and `Variation(R) … use R.delta() … R·Hat(η)`); `TypeError: Transpose() needs an operand: Transpose(expr)` | `TestCalculusOfManifoldLeaf` (5); `tests/test_geometry.py::TestTransposeConstruction::test_no_operand_raises` (replaces five childless-Transpose tests, see below) |
| 2.9 | `geomech/core/math/collect.py:33-40` (non-scalar argument), `:105-113` (both-sides check now precedes the `l == vec` / `r == vec` shortcuts) | `TypeError: collect: expected a scalar expression, got Vector 'v' of type VECTOR`; `NonLinearError: collect: 'eta' appears on both sides of 'Dot(eta(a)'b,eta)'; the expression is not linear in the target` | `TestCollectGuards` (4); the phase1-final "one ok case" (`Dot(MVMul(VVMul(η,·),·), η)`) is the pinned reproducer |
| 2.10 | `geomech/core/transformations/simplify.py:60-72` (`for … else: warnings.warn(…)`) | `RuntimeWarning: full_simplify: no fixed point after 1 passes; result may be unsimplified: <repr of the result>` | `TestFullSimplifyCapWarning` (3): a genuine two-pass expression `Dot(Hat(q)(ω×q), q)` warns at `max_iter=1` and reaches `0` without a warning at the default |
| 2.11 | `geomech/core/base/expressions.py:518-530` (`Number`: `Zero` for `value == 0`, `Ones` for `value == 1`), `:142,234,399` (operator wrapping goes through `Number`), `geomech/core/operations/multiplication.py:24` (`_wrap_numeric` → `Number`), `geomech/core/transformations/simplify.py:453-463` (`_simplify_mul`: a folded 0 returns the `Zero` leaf; other folds via `Number`), `:397`, `:476-480` | no exception; `Number(0).is_zero`, `full_simplify(Number(0) * x) → 0`, `full_simplify(Mul(Number(0), Number(3))) is Zero` | `tests/test_simplify_numbers.py` (14) |
| 2.12 | `examples/rigid_pendulum.py:32` (`J2 = Matrix("J2", …)`) | — | `TestExampleSymbolNames` (1); the example runs and prints `J1`/`J2` on the two rows |
| exports | `geomech/__init__.py:112-117`, `geomech/utils/__init__.py:1-6` | `AlgebraicError`, `NonLinearError`, `ExpressionMismatchError`, `SizeMismatchError` importable from `geomech` | `TestErrorExports` (2) |

`str(Number(0))` is still `(0)`; only the flags changed, so no pinned string moved. The
full suite emits no `RuntimeWarning` from 2.10 (no pipeline is near the cap of 10).

## Decisions applied (round 2)

1. **2.2 strict, with `params=`.** `compute_eom(lagrangian, inf_work, variables,
   params=())` declares prescribed time-varying parameters, mirroring
   `numeric.check_eom(params=...)`. `_check_params` (`eom.py:120-162`) requires distinct,
   non-constant plain `Scalar`/`Vector`/`Matrix` leaves that are neither configuration
   variables (by printed name, across all groups) nor tangent/variation vectors nor
   manifold points. After extraction, every surviving `Variation` node whose printed form
   is not a listed variation vector and whose operand is not a declared param, and every
   ξ/η/δω/δΩ of an unlisted manifold point, raises the `ValueError` quoted in the table.
   `Variation(p)` of a declared param is dropped as before (its variation is not free);
   its rate terms are derived as before. `to_standard_form` on such an equation works
   (`M` carries `l`, `f` carries `l̇`); `linearize` still raises its clear
   `AlgebraicError` on the `δl` term (`… contains no perturbation symbol …`), unchanged.
   The four prescribed-parameter tests now pass `params=[l]` / `params=[p]` with their
   assertions untouched: `tests/test_numeric.py:204`,
   `tests/test_phase1_adversarial.py:186`, `tests/test_scalar_adversarial.py:556,847`.
   Documented in the `compute_eom` docstring, `docs/getting-started.md` (paragraph after
   the first derivation), `docs/numeric-verification.md` (the "Unlisted time-varying
   parameters" bullet) and CHANGELOG (Added: `params`; Changed: the strict check with the
   old symptom).
2. **2.8 `Transpose()` → `TypeError`** (`geometry.py:123-127`). `tests/test_geometry.py`:
   the four listed tests (`test_none`, `test_expr_none`, `test_is_constant_none`,
   `test_is_zero_none`) plus a fifth of the same pattern that the list missed
   (`TestTransposeProperties::test_size_none`, `Transpose().size is None`) are replaced by
   `TestTransposeConstruction::test_no_operand_raises`, whose docstring records why
   (dbg-F8: a childless node had `type None` and leaked into the engine).
3. **`tests/test_vbl_validation.py:60`**: `"Symmetric"` → `"SymmetricMatrix"`. The test's
   assertion is unchanged; it was the one red test of round 1.
4. **`examples/sanity_check.py`**: the two calls that the engine rejects by design,
   `xi.delta()` (line 83) and `eta.delta()` (line 92), are removed with a comment; the
   script still exercises every other constructor and now exits 0. Its failure predated
   this branch.

## Checks run (final)

- `ruff format --check geomech/ tests/ examples/` and `ruff check geomech/ tests/`: pass.
- `python -m pytest tests/ -q`: 1630 passed, 0 failed, 5 xfailed. New files:
  `tests/test_validation.py` (68 tests, one class per item in the style of
  `tests/test_vbl_validation.py`, plus `TestParamsValidation`), `tests/test_simplify_numbers.py`
  (14). Existing tests touched: the four `params=[…]` additions (decision 1), the
  five-for-one replacement in `tests/test_geometry.py` (decision 2), the one-token typo
  fix in `tests/test_vbl_validation.py` (decision 3).
- `mkdocs build --strict`: pass (the docs extras were installed into the development
  virtual environment with `pip install -e ".[docs]"`).
- `python specs/reviews/scripts/engine_fuzz.py --quick`: 0 wrong-value cases, 42
  findings as at baseline, same exception classes (`AlgebraicError`, `NonLinearError`,
  one pre-existing `ExpressionMismatchError` in `simplify[adversarial]` from the
  row-vector encoding, phase 3.1/3.2). Expected shifts from 2.9: `collect[target occurs
  several times]` ok 14 → 10 / `NonLinearError` 23 → 27, and `collect[target inside matrix
  factor]` ok 1 → 0 / `NonLinearError` 11 → 12 (the "one ok case" of the phase 1 final
  review now raises). The generator itself built `TimeDerivative(S2)` directly, which 2.8
  rejects, so `engine_fuzz.py:627-636` (`wrap_calc`) now uses the eager `.t_diff()`/
  `.delta()` for manifold leaves and `:1103` skips a bare manifold leaf in the
  delta/t_diff-rule campaign; this reshuffles the per-row case counts of the `extract`
  campaigns (same random stream, different consumption) with no new category.
- Examples: `action_integral.py`, `nary_check.py`, `point_mass.py`, `quick_check.py`,
  `rigid_pendulum.py`, `sanity_check.py` all exit 0.
- CHANGELOG: one bullet per item under Unreleased / Changed with the old symptom, the
  `params` bullet under Added, plus the export and example bullets.

## Side effect to know about

The first `pip install -e ".[docs]"` resolved to a `pip` outside the development virtual
environment and installed `geomech` plus the docs extras there.
`geomech` was uninstalled from there again; `mkdocs`/`mkdocstrings`/`mkdocs-material` and
friends may now be present in that environment if they were not before.
