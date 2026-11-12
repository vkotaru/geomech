# Review: scalar configuration variables (commits 16b9f57, e06e4c0, 31a0c45)

Spec: `specs/scalar-config-variables.md`.
Tests: `tests/test_scalar_adversarial.py` (46 pass, 20 strict-xfail; each xfail is one
finding below). Full suite after adding the file: 1082 passed, 20 xfailed; ruff clean.

Every expected value was derived by hand from the Lagrangian before running the pipeline.
Convention used throughout (same as the S2/SO3 pipeline tests): the pipeline returns the
coefficient of the variation after IBP, i.e. `-(d/dt ∂L/∂q̇ - ∂L/∂q) + f = 0`, and
`to_standard_form` gives `Σ M[a] a + f + Σ G[u] u = 0` with a *negative* mass coefficient.
The scalar path follows that convention consistently (`m(-1.0)` for spring-mass,
`I(-1)mll` for the S2 pendulum); no sign discrepancy was found.

## Findings, by severity

### F1 — WRONG MATH: IBP drops `d/dt` of a non-constant scalar factor hoisted outside `Dot(ξ̇, ·)`

**Where.** `geomech/core/math/ibp.py`, `_apply_ibp`, `case Mul()`: when the target is not a
direct factor it recurses with `Mul(_apply_ibp(left), _apply_ibp(right))`. The `Dot` branch
then differentiates only the Dot's right operand. Any scalar factor that `pull` hoisted out
of the Dot (which it always does for `Dot(ξ̇, SVMul(ω, s))`) is treated as a constant.

The branch is pre-existing (it is byte-identical in `9eea193`), but before this feature every
scalar that could multiply a kinetic-energy `Dot` was a `Constant`, so it was latent. The spec's
own target system makes it live.

**Expected (hand derivation, elastic spherical pendulum, spec row 4).**
`δKE = m l̇ δl̇ + m l (ω·ω) δl + m l² ω·ξ̇`; IBP of the last term is
`−ξ·d/dt(m l² ω) = −ξ·(2 m l l̇ ω + m l² ω̇)`. ξ equation:

```
−m l² ω̇ − 2 m l l̇ ω − m g l (q × e3) + f_q = 0
```

**Code produced** (no `l̇` anywhere):

```
(\frac{d}{dt}(\omega_{q})(-1)mll+Cross(q,e3)mgl(-1)+f_q)
```

i.e. the fixed-length spherical pendulum. `to_standard_form` therefore reports
`f = Cross(q,e3)mgl(-1)` instead of `−2 m l l̇ ω − m g l q×e3`. The δl equation is correct
(see F-OK below), so the error is confined to the vector equations of every mixed system.

**Same bug, three other exhibits.**

1. R³ × R with `KE = ½ m |ẋ + l̇ e1|²`. Hand: δx: `−m ẍ − m l̈ e1 + F = 0`;
   δl: `−m e1·ẍ − m l̈ + f_l = 0`. Code: the δl equation is exact
   (`M[ẍ] = e1m(-1)`, `M[l̈] = m(-1)`), but the δx equation is
   `(\frac{d}{dt}(\frac{d}{dt}(x))(-1)m+F)`: the `−m l̈ e1` column is gone, and the mass
   matrix is not symmetric. The term lost is `m l̇ Dot(δẋ, e1)` → should be
   `−δx·(m l̈ e1)`, code gives `m l̇ Dot(−δx, ė1) = 0`.
2. Elastic pendulum in the `x = l q`, `KE = ½ m ẋ·ẋ` form: ξ equation has
   `ω l̇ l m (-1.0) + ω l l̇ m` (the two non-IBP terms, which should cancel) and no
   `−2 m l l̇ ω`. Net Coriolis coefficient 0 instead of `−2 m l l̇`.
3. SO(3) × R sliding mass `x = R ρ l`: η equation lacks `J̇_eff Ω = −2 m l l̇ ρ×(ρ×Ω)`;
   the two `l̇` terms present cancel.
4. Pure S2 pendulum with a *prescribed* non-constant `l(t)` (no `scalars` at all) loses
   `2 m l l̇ ω` — demonstrates the bug does not need the new code path to trigger.

**Minimal repro.**

```python
a = Scalar("a"); xi, om = Vector("xi"), Vector("om")
integrate_by_parts(Mul(a, Dot(TimeDerivative(xi), om)), [TimeDerivative(xi)])
# → a(-1)Dot(xi, d/dt(om))          expected: −ȧ Dot(xi, om) − a Dot(xi, d/dt(om))
```

**Tests.** `TestElasticSphericalPendulum::test_xi_eom_has_coriolis_term`,
`::test_xi_f_has_coriolis_term`, `TestElasticSphericalPendulumQdotForm::test_xi_f_exact_structure`,
`TestR3xRCoupledMassMatrix::test_x_equation_has_l_ddot`, `::test_mass_matrix_symmetric`,
`TestIntegrateByPartsScalarFactors::test_vector_target_with_nonconstant_prefactor`,
`::test_s2_pendulum_with_time_varying_length`,
`TestStandardFormHygiene::test_eta_equation_has_inertia_rate_term`.

**Fix direction.** In `_apply_ibp`, for `Mul` whose chain contains a `Dot(target, rhs)` factor,
IBP the whole product: `Mul(prefactor, Dot(target, rhs)) → Dot(−∫target, d/dt(prefactor · rhs))`,
i.e. push the prefactor back inside before differentiating (or differentiate the product
`prefactor * rhs` with the product rule). The new scalar path (`_remove_factor`) already does
this correctly for scalar targets (`TestIntegrateByPartsScalarFactors::test_scalar_target_with_nonconstant_prefactor` passes);
the vector path needs the same treatment.

### F2 — SILENT WRONG ANSWER: `extract_linear_coeff` returns zero for `Hat`, `Transpose`, `Vee`

`_scalar_from_matrix` has no `Hat` case and `_scalar_from_vector` has no `Transpose`/`Vee`
case; all fall to `case _: return Zero*`. These operators are linear, so the coefficient of
`l` in `Hat(a l)` is `Hat(a)`, never zero.

```python
extract_linear_coeff(Hat(SVMul(a, l)), l)              # → 0      expected Hat(a)
extract_linear_coeff(MVMul(Hat(SVMul(a, l)), w), l)    # → 0w     expected Hat(a)w
extract_linear_coeff(Transpose(SVMul(a, l)), l)        # → 0v     expected (a)'
extract_linear_coeff(Vee(SMMul(M, l)), l)              # → 0v     expected Vee(M)
```

Through `compute_eom` this is hard to reach because `_eliminate` rewrites `Hat(v) w` to
`Cross(v, w)` before extraction (verified: the SO(3)×R system written with `Hat(ρ l)` gives
the same δl equation as the `Cross` form). It is reachable through the public API and through
any expression where a `Hat`/`Transpose` survives simplification. The spec says unsupported
patterns must raise a typed error naming the pattern; these return a wrong value instead.
Also note `0w` is `MVMul(ZeroMatrix, w)` left unsimplified.

Tests: `TestScalarQueryInsideVectorOps::test_silent_zero_branches[Hat|Hat*vec|Transpose|Vee]`.

### F3 — SILENT WRONG ANSWER: higher time derivatives vanish

1. `to_standard_form` on `m x⃛ + k x`: `has(ẍ)` is true for the jerk term (it recurses into
   `TimeDerivative`), so the term is put in the acceleration bucket; `extract_linear_coeff`
   treats `TimeDerivative` as opaque and returns `Zero`. Result: `M[ẍ] = 0`, `f = kx`; the jerk
   term is gone with no error.
2. `compute_eom(L = ½ m ẍ², f δx)` returns the EOM `f`: `δL = m ẍ δẍ`, `δẍ` is not an IBP
   target, the term survives IBP and is then opaque to extraction. The whole Lagrangian is
   discarded silently. (Correct answer is `−m d⁴x/dt⁴ + f = 0`; a typed error would satisfy
   the spec.)

The `has`/`_has_algebraic` disagreement in (1) is the structural cause: the partition in
`standard_form._extract_single` uses the recursive `has`, extraction uses the opaque
`_has_algebraic`. Any term where they disagree is dropped. Suggested guard: after extracting
`M`, verify `accel_expr − Σ M[a]·a` simplifies to zero, or make the partition use
`_has_algebraic`.

Tests: `TestNonLinearAcceleration::test_jerk_term_is_not_silently_dropped`,
`::test_higher_derivative_lagrangian_raises`.

### F4 — WRONG ANSWER on name aliasing: `_has_algebraic` compares `str()` without a type check

`Mul.__str__` is bare concatenation, so `Mul(k, x)` prints as `kx` and `Mul(x, One)` prints as
`x1`. `Expr.__eq__`/`has` check `type()` first and are safe; `_has_algebraic` and the
`str(expr) == str(query)` short-circuits in `_scalar_from_scalar`/`_vector_from_vector` do not:

```python
extract_linear_coeff(Mul(k, x), Scalar("kx"))   # → 1   expected 0
extract_linear_coeff(Mul(x, One), Scalar("x1")) # → 1   expected 0
```

`One` is produced by `_scalar_from_scalar` itself, so `Mul(·, One)` does occur inside the
extraction. In the pipeline the queries are `\delta{x}`, `\frac{d}{dt}(…)` and input names,
which makes a collision unlikely but not impossible (an input named like a product of two
symbols). The pipeline tests with `x`/`x1` and `l`/`l0` pass (`TestNameCollisions::test_x_and_x1_coupled_exact`).
Fix: add `type(expr) is type(query)` to the string comparison, or compare with `==`.

Tests: `TestNameCollisions::test_product_string_does_not_alias_a_leaf`,
`::test_x_times_one_does_not_alias_x1`.

### F5 — MISSING ERROR: a `Constant` scalar in `SystemVariables.scalars`

`Scalar("xc", attr=["Constant"]).delta()` is `Zero`, so `compute_eom` adds the entry
`eom["0"] = (Zero, Zero)` and carries on. Should raise. Test:
`TestInputValidation::test_constant_scalar_in_scalars_raises`.

### F6 — QUALITY: numeric coefficients and like-term combination

1. `-1` vs `-1.0`. The same arithmetic prints as `(-1)` in some terms and `(-1.0)` in
   others, even inside one equation: spring-mass-damper
   `m ẍ(-1.0) + k x(-1.0) + f + c ẋ(-1)`; elastic δl: `m l̈(-1) … k l(-1.0)`. Cause: `½ m ẋ δẋ`
   is produced twice by the product rule as `Mul(Mul(m, δv), v)` and `Mul(Mul(m, v), δv)`,
   which `_commutative_equal` does not identify (only the top-level `Mul` is commuted), so
   they are IBP'd separately and summed later as `0.5 + 0.5 = 1.0` float. The S2 path
   combines first and drops the `1.0`. Every exact-string test is therefore pinned to an
   artefact of evaluation order.
   Test: `TestSpringMassDamperExact::test_numeric_coefficients_have_one_type`.
2. Like terms with a buried numeric do not combine. `½ (m0 + m1 x) ẋ²` gives
   `f = 0.5 m1 ẋ ẋ + m1 ẋ ẋ (-1.0)` (should be `m1 ẋ ẋ (-0.5)`); `_split_coeff` only looks at
   the outermost `Mul`. Same cause leaves `ω l̇ l m (-1.0) + ω l l̇ m` uncancelled in the
   elastic `qdot` form (factor order `l̇ l` vs `l l̇`).
   Test: `TestVaryingInertia::test_f_combines_like_terms`.
3. `M` keeps identically-zero coefficients: SO(3)×R δl equation reports `M[Ω̇] = 0v`
   (from `ρ·(Ω̇×ρ)`, which is only recognised as zero after extraction).
   Test: `TestStandardFormHygiene::test_zero_coefficients_are_pruned_from_M`.
4. `Dot(R a, R b)` is only stripped inside `collect`, not in `simplify`, so
   `M[l̈] = m Dot(Rρ, Rρ)(-1.0)` instead of `m(-1)`; `(ω×q)·(ω×q)` is not reduced to `ω·ω`
   in the elastic `qdot` δl equation. Pre-existing, not scalar-specific.
5. Row-vector form `MVMul(Transpose(a l), M)` crashes `extract_linear_coeff` with
   `ExpressionMismatchError: MVMul cannot be performed between MATRIX and MATRIX`
   (untyped, misleading). Rare.
6. `expand` for `Vee` tests `isinstance(inner, VAdd)` but `Vee`'s operand is a matrix
   (`MAdd`); linearity never fires. Pre-existing, unrelated to scalars.

## Verified correct (passing tests, exact strings)

- Spring-mass, spring-mass-damper, double spring-mass: EOM, `M`, `f`, `G` exact and equal to
  the hand derivation (`TestSpringMassExact`, `TestSpringMassDamperExact`,
  `TestDoubleSpringMassExact`).
- Sign convention consistent with the S2 pipeline (`test_sign_convention_matches_s2_pendulum`).
- Coupled scalar accelerations `L = ½ m (ẋ1 + ẋ2)²`: both `M` keys present and `−m`.
- Varying inertia `½ (m0 + m1 x) ẋ²` and `½ m x² ẋ²`: mathematically correct
  (`M = −(m0 + m1 x)`, `M = −m x²`); only the term combination is poor.
- Elastic spherical pendulum δl equation exact:
  `−m l̈ + m l ω·ω − m g q·e3 − k l + k l0 + f_l`; `M[ω̇] = I(-1)mll`, `G[f_q] = I`; no
  cross-contamination between δl and ξ.
- Scalar inside `Cross`, `Dot`, `MVMul`, `VVMul`: coefficients exact
  (`Hat(b)(-1)a`, `Dot(Hat(b)(-1)a,c)`, `Ma`, `(a)'b`); quadratic occurrences raise
  `NonLinearError`.
- Scalar IBP target inside `Dot` after `pull`, with and without a non-constant prefactor,
  is correct (`_remove_factor` path keeps `ȧ`).
- `ẍ²`, `m ẍ ẍ`, `ẍ² + m x`, `SVMul(a, ẍ ẍ)` all raise `NonLinearError`; matrix query raises
  `NotImplementedError`.
- `c ẋ ẍ + m ẍ + k x` partitions to `M = (c ẋ + m)`, `f = k x`; `ẋ` never matches `ẍ`.
- `x`/`x1`, `l`/`l0` do not collide in `has` or in the EOM dict keys.

## Verdict

**Do not merge as-is.** The spec's own fourth system (R × S² elastic pendulum) comes out
wrong: the ξ equation is missing `−2 m l l̇ ω` (F1), and every mixed scalar/vector system
loses the scalar-rate coupling in its vector equation (R³×R mass matrix asymmetric,
SO(3)×R missing `J̇ Ω`). The three pure-scalar systems are exact.

Must fix before merge:

1. F1 — IBP of `Mul(prefactor, Dot(ξ̇, rhs))` must differentiate the prefactor. Add the
   elastic pendulum ξ equation and the R³×R cross term to the pipeline tests as exact
   strings once fixed.
2. F3 — make the acceleration partition and extraction agree (or verify the residual is
   zero) so that no term can be dropped silently; raise on `δẍ`.
3. F2 — give `_scalar_from_matrix`/`_scalar_from_vector` `Hat`/`Transpose`/`Vee` cases (or
   raise `NotImplementedError` naming the node) instead of returning zero.

Should fix (small): F4 type check in `_has_algebraic`; F5 validation of `scalars`.
Can defer: F6 (quality), though item 1 means the exact-string tests pin an evaluation-order
artefact and will churn whenever simplification changes.

## Fix log

Scope: F1–F5 fixed; F6 untouched (deferred). Suite after the fixes:
`1123 passed, 7 xfailed, 1 failed` (the one failure is a pre-existing test that pins F1, see
"Pre-existing test that pins the bug" below); ruff format + check clean.
Files touched: `geomech/core/math/ibp.py`, `geomech/core/math/extract.py`,
`geomech/core/math/collect.py`, `geomech/dynamics/standard_form.py`, `geomech/dynamics/eom.py`,
`tests/test_scalar_adversarial.py` (13 xfail markers removed, section 9 added),
`tests/test_ibp_prefactor.py` (new), this file.

### F1 — IBP under a non-constant scalar prefactor (`ibp.py`)

**Bug.** `_apply_ibp`'s `Mul` branch recursed `Mul(_apply_ibp(l), _apply_ibp(r))`, so the `Dot`
branch differentiated only the Dot's right operand; a scalar prefactor hoisted by `pull` was
treated as constant.

**Fix.** The `Mul` branch now splits the chain into `(prefactor, Dot(target, rhs))`
(`_split_dot_factor`). The `a·Dot(−x, ṙ)` term is produced by replacing the Dot *in place*
(`_replace_in_chain`), so for a constant prefactor the output tree is byte-identical to before
and every pinned string in the SO(3)/S² pipelines is unchanged. For a non-constant prefactor
the product rule's other half, `−Dot(x, ȧ r)`, is added, with `ȧ` formed by `_t_diff_product`:
`d/dt(∏ fᵢ^{nᵢ}) = Σ nᵢ fᵢ^{nᵢ−1} ḟᵢ ∏_{j≠i} fⱼ`, so `d/dt(m l l) = 2 m l l̇` is one term and the
Coriolis term comes out as a single addend `ω m l l̇ (−2)`. Replacing in place also avoids
double-counting when the Dot sits in the middle of a chain (`Mul(Mul(a, Dot), b)`), which a
naive "recurse then add" produced.

**Verified by hand (elastic spherical pendulum, spec row 4).** ξ equation now
`(ω̇(−1)mll + ω m l l̇(−2) + Cross(q,e3)mgl(−1) + f_q)`, i.e.
`−m l² ω̇ − 2 m l l̇ ω − m g l (q×e₃) + f_q = 0`; `M[ω̇] = I(−1)mll` unchanged, `f` has the
two expected addends. R³×R: δx equation `(ẍ(−1)m + e1 m l̈(−1) + F)`, `M_x[l̈] = e1m(−1)`
equals `M_l[ẍ]` (mass matrix symmetric). Prescribed `l(t)` S² pendulum gains `2 m l l̇ ω`.

**Tests.** `TestElasticSphericalPendulum::test_xi_eom_has_coriolis_term`,
`TestR3xRCoupledMassMatrix::test_x_equation_has_l_ddot`,
`TestIntegrateByPartsScalarFactors::test_vector_target_with_nonconstant_prefactor`,
`::test_s2_pendulum_with_time_varying_length` (markers removed);
`TestFixF1PipelineExact` (exact strings for the three systems above);
`tests/test_ibp_prefactor.py` (nine pure IBP identities with the hand derivation in each
docstring, including constant prefactor unchanged, repeated factor → single `(−2)` term,
constant rhs, Dot mid-chain not double counted).

**Pre-existing test that pins the bug.** `tests/test_ibp.py::TestIbpRecursion::test_mul_recurses`
asserts `_apply_ibp(Mul(s, Dot(η̇, w)), η̇)` is a `Mul` with left `s` where `s = Scalar("s")` is
*not* constant. Correct math is `−ṡ Dot(η, w) − s Dot(η, ẇ)` (two terms; `ṡ ≠ 0`), so the
assertion can only hold with the bug. Not edited; left failing for arbitration. If `s` were
declared `attr=["Constant"]` the test would pass unchanged.

**Side finding exposed by F3 (fixed in `collect.py`).** With the F3 rule in place, the SO(3)×SO(3)
double rigid pendulum raised: the coupling term `m₂ Dot(R₁(η̇₁×l₁), R₂(Ω₂×ρ₂))` survived IBP
because `_collect_dot`'s `(MVMul, MVMul)` case returned unchanged when the two rotations
differ, and extraction then dropped it silently (the old M had no `Ω̇₂` block in the η₁
equation). The case now peels the matrix off the side holding the target,
`Dot(M f(v), w) = Dot(f(v), Mᵀ w)`, and recurses. Result, checked by hand
(`∂KE_c/∂Ω₁ = −m₂ l̂₁ R₁ᵀR₂ ρ̂₂ Ω₂`): `M₁[Ω̇₂] = Hat(l_1)(R1)'R2Hat(ρ_2)m2`,
`M₂[Ω̇₁] = Hat(ρ_2)(R2)'R1Hat(l_1)m2 = M₁[Ω̇₂]ᵀ`. The existing double-pendulum tests
(substring checks) still pass. Test: `TestCollectTwoRotations::test_cross_inertia_blocks_exact`.

### F2 — `Hat`, `Transpose`, `Vee`, outer product for a scalar query (`extract.py`)

**Bug.** `_scalar_from_matrix` had no `Hat`/`Transpose`/`VVMul` case and `_scalar_from_vector`
no `Transpose`/`Vee` case; all fell to `case _: return Zero*`.

**Fix.** Linearity cases added with one-line math docstrings: `Hat(c l) = l Hat(c)`,
`(c l)ᵀ = l cᵀ`, `Vee(C l) = l Vee(C)`, `(a l) bᵀ = l (a bᵀ)`. The `MVMul` case dispatches
each factor on its *type* (`_scalar_from`), so the row-vector form `(a l)ᵀ M` now gives
`(a)'M` instead of crashing (review F6.5, incidentally). Every remaining `case _` fallback
in `extract.py` goes through `_zero_or_raise`: zero only if the node is free of the query,
otherwise `AlgebraicError` naming the node — except `_vector_from_matrix`, whose silent
`ZeroMatrix` for a vector query in `Hat(v)·M` is pinned by three tests in
`tests/test_extract.py::TestNotImplementedErrors` (out of this feature's scope; left as is).

**Tests.** `TestScalarQueryInsideVectorOps::test_silent_zero_branches[Hat|Hat*vec|Transpose|Vee]`
(markers removed); `TestFixF2LinearUnaryOps` (row-vector × matrix, outer product both sides,
matrix transpose, quadratic `Hat(a l)(w l)` raises `NonLinearError`).

### F3 — higher derivatives must raise, never vanish (`extract.py`, `standard_form.py`)

**Rule adopted.** A query that occurs only *inside* a calculus operator which is not itself
the query (ẍ inside d/dt(ẍ), δx inside d/dt(d/dt(δx))) raises `AlgebraicError` with the
message "… appears only inside the calculus operator … — a derivative of higher order than
supported". Implemented in `_has_algebraic` itself, so every extraction path inherits it;
`standard_form._extract_single` partitions with `_has_algebraic` instead of the recursive
`has`, so the partition and the extraction agree and no term can be bucketed as
acceleration/input and then extracted as zero. Lower derivatives stay opaque and silent
(ẋ contains no ẍ; `c ẋ + m ẍ` still partitions to `M = m`, `f = c ẋ`).

**Tests.** `TestNonLinearAcceleration::test_jerk_term_is_not_silently_dropped`,
`::test_higher_derivative_lagrangian_raises` (markers removed);
`TestFixF3HigherDerivativeRule` (message match in extraction and in the partition, δẍ through
`compute_eom`, lower derivative unaffected).

### F4 — typed identity in `_has_algebraic` (`extract.py`)

**Bug.** `str(expr) == str(query)` with no type check, so `Mul(k, x)` ("kx") matched
`Scalar("kx")` and `Mul(x, One)` ("x1") matched `Scalar("x1")`.

**Fix.** `_same(expr, query)`: `type(expr) is type(query)` (the guard `Expr.__eq__` uses), then
leaves compare by `name`, non-leaves by rendered string. Used by `_has_algebraic` and by the
`One`/`IdentityMatrix` short-circuits in `_scalar_from_scalar`, `_vector_from_vector`,
`_matrix_from_matrix`. No structural equality attempted.

**Tests.** `TestNameCollisions::test_product_string_does_not_alias_a_leaf`,
`::test_x_times_one_does_not_alias_x1` (markers removed); `TestFixF4TypedIdentity`.

### F5 — constant configuration variable (`eom.py`)

**Fix.** `_check_variables` runs first in `compute_eom`: each entry of `scalars`/`vectors`/
`matrices` must have the matching `ExprType` (`TypeError` otherwise) and a non-vanishing
variation — `is_constant`, `is_numeric`, `is_zero` or `delta().is_zero` raise `ValueError`
naming the group and the variable. Applied to all three groups; no existing test passes a
constant.

**Tests.** `TestInputValidation::test_constant_scalar_in_scalars_raises` (marker removed);
`TestFixF5Validation` (constant vector, constant matrix, wrong type in either direction).

### xfail markers left in place (5, after arbitration)

All five are F6 (deferred like-term / numeric-type quality items):

- `TestSpringMassDamperExact::test_numeric_coefficients_have_one_type`,
  `TestVaryingInertia::test_f_combines_like_terms`,
  `TestStandardFormHygiene::test_zero_coefficients_are_pruned_from_M`.
- `TestElasticSphericalPendulumQdotForm::test_xi_f_exact_structure` — the IBP term
  `ω m l l̇(−2)` is present since F1; the non-IBP pair `ω l̇ l m(−1.0) + ω l l̇ m` is not
  combined (F6.2), so `f` has three addends instead of two.
- `TestStandardFormHygiene::test_eta_equation_has_inertia_rate_term` — the IBP term
  `−2 m l l̇ ρ×(ρ×Ω)` is present since F1; the two cancelling non-IBP `l̇` terms are not
  combined (F6.2), so `f` has three `l̇` addends, not one.

Arbitration by the session owner's agent after the fix step:

- `TestElasticSphericalPendulum::test_xi_f_has_coriolis_term` — the original assertion
  asked for `isinstance(f, Add)`; a vector sum is a `VAdd` in this library. Assertion
  corrected to `VAdd`, marker removed; the test now passes and pins the F1 fix.
- `TestR3xRCoupledMassMatrix::test_mass_matrix_symmetric` — deleted. Its third assertion
  compared a row vector with a column vector by string and could never hold; the
  representation-correct check lives in
  `TestFixF1PipelineExact::test_r3_x_r_cross_inertia_exact`.

### Final review follow-ups (see `scalar-config-variables-final.md`)

- Finding 1: `collect` flip loop guarded (`(_, MVMul())` only flips when the target is
  inside the `MVMul`; `(MVMul(), _)` returns unchanged when the target is nested in an
  unexposable form). Tests: `tests/test_collect_termination.py`.
- Finding 2: double rigid pendulum mass-matrix blocks pinned exactly in
  `tests/test_manifold_pipelines.py::TestDoubleRigidPendulumMassMatrixExact`.
- Finding 3: `to_standard_form` raises `NonLinearError` when a mass-matrix coefficient still
  contains an acceleration. Test in `tests/test_collect_termination.py`.
- Finding 4: `_has_algebraic` error text no longer claims "higher order derivative" for
  integrals and variations.
- Findings 5, 6, 7, 9 (two product-rule helpers, private cross-module import, query-free
  guard on unary branches, docstring signature): deferred to the structural-equality /
  extraction feature; not user-visible.
- Pre-existing `tests/test_ibp.py::TestIbpRecursion::test_mul_recurses` pins the old,
  wrong IBP output for a non-constant prefactor. Left untouched pending the owner's decision.
