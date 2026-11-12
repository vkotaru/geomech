# Final review: scalar configuration variables (branch `feature/scalar-implementations`)

Spec: `specs/scalar-config-variables.md`. Adversarial review + fix log:
`specs/reviews/scalar-config-variables-review.md`.
Scope reviewed: `git diff origin/geomech...HEAD -- geomech/ tests/`, the uncommitted diff to
`geomech/core/math/{ibp,collect,extract}.py`, `geomech/dynamics/{eom,standard_form}.py`,
`tests/test_scalar_adversarial.py`, and the untracked `tests/test_ibp_prefactor.py`.

State at review: `1124 passed, 6 xfailed, 1 failed`; `ruff format --check` and `ruff check`
clean. The one failure (`tests/test_ibp.py::TestIbpRecursion::test_mul_recurses`) is confirmed
below as a pre-existing test that pins the F1 bug; it is the owner's call.

## Verdict: GO WITH CHANGES

The mathematics is right. Every identity I re-derived independently and ran through the code
came out correct, including the one that changes a previously "verified" system (the SO(3)×SO(3)
double rigid pendulum), which I checked numerically against a from-scratch hand derivation of
the full η₁ and η₂ equations. None of the findings below is a wrong equation. The changes
required before merge are small: one robustness fix in `collect.py` (a new path into a
pre-existing infinite recursion), the CHANGELOG entry, and two documentation/test-hygiene
corrections. The rest are recommendations.

## A. Mathematical correctness (re-derived)

### A1. `ibp.py` — IBP under a non-constant scalar prefactor (F1) — CORRECT

Identity: ∫ a·Dot(ẋ, r) dt = −∫ [ȧ·Dot(x, r) + a·Dot(x, ṙ)] dt.

Cases I constructed and ran (`_apply_ibp` and the full `integrate_by_parts`), with `a, b, l`
non-constant and `c, d` constant, `D = Dot(ξ̇, r)`:

| Input | Output (after `full_simplify`) | Correct |
|---|---|---|
| `Mul(Mul(D,a),b)` (Dot first) | `−ab Dot(ξ,ṙ) − ȧb Dot(ξ,r) − aḃ Dot(ξ,r)` | yes |
| `Mul(Mul(a,b),D)` (Dot last) | same three terms | yes |
| `Mul(Mul(a,D),b)` (Dot middle) | same three terms, once each | yes |
| `Mul(a,Mul(D,b))`, `Mul(a,Mul(b,D))` | same | yes |
| `Mul(Mul(Mul(a,b),l),D)` | `−abl Dot(ξ,ṙ) − ȧbl… − aḃl… − abl̇…` (all three product-rule terms) | yes |
| `c l l D` | `−cll Dot(ξ,ṙ) + c l l̇ Dot(ξ,r)(−2)` | yes |
| `a l l D` | `… + ȧll(−1) + a l l̇(−2)` | yes |
| `Mul(c,D)`, `Mul(Mul(c,d),D)` (constant) | raw `_apply_ibp` output `cDot(ξ(-1),ṙ)` / `cdDot(ξ(-1),ṙ)`: byte-identical to the pre-fix form; `Mul(c,d).is_constant` is True | yes |
| `Dot(w,w) · D` (non-leaf scalar factor) | `… + Dot(ẇ,w)Dot(ξ,r)(−2)` | yes |
| `(c + a) · D` (Add factor) | `−c Dot(ξ,ṙ) − a Dot(ξ,ṙ) − ȧ Dot(ξ,r)` | yes |

The in-place replacement (`_replace_in_chain`, identity-based) is why the constant-prefactor
tree is unchanged; `_split_dot_factor` returns the actual node so identity holds. Single-rotation
SO(3) and S² pipeline tests pin exact strings (`test_manifold_pipelines.py:244-489`) and are
green, so no pre-fix output changed.

Quadratic-in-target chains (`Mul(Dot(ξ̇,r), Dot(ξ̇,w))`) produce garbage, but a variation is linear
by construction and this was equally wrong before; not a regression.

### A2. `collect.py` — `Dot(M f, w) = Dot(f, Mᵀ w)` and the double rigid pendulum — CORRECT

The identity (M f)·w = fᵀMᵀw = f·(Mᵀw) holds for any matrix M, not just orthogonal, so applying
it whenever `l` is an `MVMul` containing the target is valid.

**Hand derivation (Lee et al. §6 two-body chain).** x₂ = R₁l₁ + R₂ρ₂, ẋ₂ = R₁(Ω₁×l₁) + R₂(Ω₂×ρ₂).
Cross kinetic energy KE_c = m₂ (R₁(Ω₁×l₁))·(R₂(Ω₂×ρ₂)) = m₂ (−l̂₁Ω₁)ᵀ R₁ᵀR₂ (−ρ̂₂Ω₂)
= −m₂ Ω₁ᵀ l̂₁ R₁ᵀR₂ ρ̂₂ Ω₂ (using l̂₁ᵀ = −l̂₁). Hence
∂KE_c/∂Ω₁ = −m₂ l̂₁R₁ᵀR₂ρ̂₂ Ω₂ and ∂KE_c/∂Ω₂ = −m₂ ρ̂₂R₂ᵀR₁l̂₁ Ω₁. In the library convention
(coefficient of η after IBP, i.e. −ṗ + …), the Ω̇₂ block of the η₁ equation is
**M₁[Ω̇₂] = +m₂ l̂₁R₁ᵀR₂ρ̂₂** and M₂[Ω̇₁] = +m₂ ρ̂₂R₂ᵀR₁l̂₁ = M₁[Ω̇₂]ᵀ. Sign is consistent with the
diagonal blocks M₁[Ω̇₁] = −(J₁ − m₁ρ̂₁ρ̂₁ − m₂l̂₁l̂₁) (the whole mass matrix is negated). The fixer's
claim is correct.

**Numeric check.** I wrote an evaluator over the expression tree (random R₁,R₂ ∈ SO(3), random
J, ρ, l, Ω, Ω̇, inputs) and compared the library's *full* η₁ and η₂ equations with my hand
expressions
η₁: −ṗ₁ + p₁×Ω₁ + m₂ (a × R₁ᵀR₂b) − m₁g (ρ₁ × R₁ᵀe₃) − m₂g (l₁ × R₁ᵀe₃) + M₁,
p₁ = J_eff,1 Ω₁ − m₂ l̂₁R₁ᵀR₂ρ̂₂ Ω₂, ṗ₁ including d/dt(R₁ᵀ) = −Ω̂₁R₁ᵀ and d/dt(R₂) = R₂Ω̂₂,
a = Ω₁×l₁, b = Ω₂×ρ₂ (η₂ symmetric). Relative error over 5 random trials: ≤ 1e-15 for both
equations and for all four M blocks, including M₁[Ω̇₂] = M₂[Ω̇₁]ᵀ.
Script: `specs/reviews/scripts/double_rigid_pendulum_numeric_check.py` (run from the repo root with the venv active).

Consequence worth stating plainly: **before this branch the double rigid pendulum was wrong**
(the whole −η₁·d/dt(∂KE_c/∂Ω₁) contribution was dropped, not just the Ω̇₂ block), and
`TestDoubleRigidPendulumPipeline` did not catch it because it only checks substrings. See C2.

### A3. `extract.py` — F2/F3/F4 — CORRECT, with two low-severity edges

New branches and the identity each implements:

- `_scalar_from_vector` `Transpose`: (c l)ᵀ = l cᵀ. `Vee`: Vee(C l) = l Vee(C).
- `_scalar_from_matrix` `Hat`: Hat(c l) = l Hat(c). `Transpose`: (C l)ᵀ = l Cᵀ.
  `VVMul` (outer): (a l)bᵀ = l (abᵀ); a(b l)ᵀ = l (abᵀ); both-sides → `NonLinearError`.
- `_scalar_from_vector` `MVMul` now dispatches each factor by type (`_scalar_from`), so the
  row-vector form (a l)ᵀM gives `(a)'M` and (a l)ᵀ(M l) raises `NonLinearError`.
- `_zero_or_raise` on every `case _` except `_vector_from_matrix` (documented, pinned by
  `tests/test_extract.py::TestNotImplementedErrors`, out of scope).

I ran 20 cases (all of the above, nested `Transpose(Transpose(·))`, `M·Vee(Hat(a l))`,
`TimeIntegral(l)`, l under `d/dt` inside a `Cross`): every result was the correct coefficient or
a typed error naming the node. No silent zero or dropped term escapes for a scalar query.

`_same` type guard — cannot miss a legitimate match: `TS2` vs fresh `TS2` of the same name
(with or without the `S2=` back-reference) → True; `TSO3` likewise; `Variation(x)`,
`d/dt(d/dt(x))` → True by rendered string; `Number(2)` vs `Scalar("(2)")` → True (same class,
same name); `Scalar("x")` vs `Scalar("x", attr=Constant)` → True. It correctly rejects
`TS2("x")` vs `Vector("x")` and `SO3("R")` vs `Matrix("R")`, which is exactly what `Expr.__eq__`
does, so extraction and `has` now agree.

`_has_algebraic` raise rule: lower derivative stays opaque (`_has_algebraic(ẋ, ẍ)` → False);
`d/dt(ẍ)`, `d/dt(δx)`, `TimeIntegral(l)` raise; a sum `Dot(ξ,w) + Dot(ξ̇,w)` raises on the
second addend. Correct.

Edges (not blocking, see findings 4 and 7): the message always says "a derivative of higher
order than supported" even when the operator is `∫` or the query is an *input* whose rate
appears; the three new unary branches have no "query-free → zero" guard, so
`extract_linear_coeff(Transpose(a), l)` returns `(0v)'` (is_zero True, but `full_simplify`
does not fold it to `0v`).

### A4. `eom.py` `_check_variables` (F5) and `standard_form.py` partition (F3) — CORRECT

`_check_variables`: S2, SO3, plain Vector, non-SO3 Matrix accepted; `ZeroVector`, numeric
`Scalar("2", value=2)`, `Scalar("c", attr=Constant)` → `ValueError` naming group and variable;
`Vector` in `scalars`, a `str` in `scalars` → `TypeError` naming group and variable. Both error
types are typed and name the offending pattern.

Partition/extraction agreement: both sides now use `_has_algebraic`, and the `Add` branches of
extraction filter addends with the same predicate, so a term classified as acceleration cannot
be extracted as zero. Verified: `m ẍ + c ẋ + k x + f`, damping-times-acceleration, `ẍ` inside
`Dot(a ẍ, b)`, jerk (raises), `x²` (fine). One pre-existing hole remains (finding 3): a term
bilinear in two accelerations `c ẍ₁ ẍ₂` gives `M[ẍ₁] = c ẍ₂, M[ẍ₂] = c ẍ₁`, so Σ M[a]·a
double-counts the term. The spec says nonlinear acceleration must never return a partial
coefficient; this case is not covered by the `ẍ²` guard.

### Pre-existing failing test

`tests/test_ibp.py::TestIbpRecursion::test_mul_recurses` constructs `s = Scalar("s")` (not
constant) and asserts `_apply_ibp(Mul(s, Dot(η̇, w)), η̇)` is a `Mul` with left `s`. Correct
math is −ṡ Dot(η,w) − s Dot(η,ẇ): two addends. The fixer's diagnosis is right. I verified that
with `s = Scalar("s", attr=["Constant"])` the test passes unchanged (output `sDot(eta(-1),ẇ)`,
a `Mul` with left `s`). Recommendation to the owner: declare `s` constant in the fixture and
keep the test; the non-constant case is now covered by `tests/test_ibp_prefactor.py`.

## Findings, by severity

### 1. MEDIUM — new `collect` branch opens a path into a pre-existing infinite recursion

`geomech/core/math/collect.py:92-93`. The peel `Dot(M f(v), w) → collect(Dot(f(v), Mᵀ w), v)`
recurses on whatever `f(v)` is. If `f(v)` is not a `Cross`/`MVMul` (e.g. `SVMul`, `VAdd`), the
result lands in the `(MVMul, _)` ↔ `(_, MVMul)` flip pair (lines 99-119), which flips forever
when the non-MVMul side holds the target and the MVMul side does not.

```python
from geomech import Dot, MVMul, SVMul, VAdd, Matrix, Vector, Scalar
from geomech.core.math.collect import collect
v, a, b, A, B, s = Vector("v"), Vector("a"), Vector("b"), Matrix("A"), Matrix("B"), Scalar("s")
collect(Dot(MVMul(A, SVMul(v, s)), MVMul(B, b)), v)   # RecursionError (returned unchanged before)
collect(Dot(MVMul(A, VAdd(v, a)), MVMul(B, b)), v)    # RecursionError (returned unchanged before)
collect(Dot(SVMul(v, s), MVMul(B, b)), v)             # RecursionError — PRE-EXISTING
```

Not reachable through `integrate_by_parts`/`compute_eom` (expand distributes the `VAdd`, pull
hoists the scalar out of the Dot first), so no pipeline result is affected; reachable through
the public `collect`. Fix either the root (`case (MVMul(), _)` line 115: flip only if
`r.has(vec)`, otherwise `return expr`) or guard the new peel to `isinstance(l.right, (Cross,
MVMul))`. Add the three lines above as tests (expect "returned unchanged" or a typed error).

### 2. MEDIUM — the double rigid pendulum test does not pin the equation it now gets right

`tests/test_manifold_pipelines.py:497-560` checks only that substrings (`J1`, `\rho_1`, `m2`,
`e3`, …) appear. That is why a wrong equation passed for the life of the test. The project
standard (the contributor notes "Testing conventions") requires the actual derived equation. Now that the
system is correct, pin the four M blocks exactly (I verified all four numerically):

```
M₁[Ω̇₁] = (J1(-1.0)+Hat(\rho_1)Hat(\rho_1)m1+Hat(l_1)Hat(l_1)m2)
M₁[Ω̇₂] = Hat(l_1)(R1)'R2Hat(\rho_2)m2
M₂[Ω̇₁] = Hat(\rho_2)(R2)'R1Hat(l_1)m2
M₂[Ω̇₂] = (J2(-1.0)+Hat(\rho_2)Hat(\rho_2)m2)
```

`TestCollectTwoRotations::test_cross_inertia_blocks_exact` already pins the two off-diagonal
blocks with the hand derivation in its docstring; the diagonal blocks and `G` belong in the
pipeline test. The f terms are long; pinning them is optional, but a docstring stating the
hand form (A2 above) would let a future reader check them.

### 3. LOW — bilinear acceleration terms are double-counted, not rejected

`geomech/dynamics/standard_form.py:254-261`. Repro:

```python
x1, x2 = Scalar("x1"), Scalar("x2"); a1, a2 = x1.t_diff().t_diff(), x2.t_diff().t_diff()
to_standard_form({k: (Variation(x1), Mul(Mul(c, a1), a2))}, SystemVariables(scalars=[x1, x2]), [])
# M = {ẍ1: c ẍ2, ẍ2: c ẍ1}   →  Σ M[a]·a = 2 c ẍ1 ẍ2
```

Pre-existing and physically unreachable from a second-order Lagrangian, but it is exactly the
"partial coefficient" the spec forbids. One guard after `M[str(a)] = coeff`: if
`any(_has_algebraic(coeff, b) for b in accel_vars)` raise `NonLinearError`.

### 4. LOW — `_has_algebraic` error text is wrong for two of the cases it covers

`geomech/core/math/extract.py:69-73`. The message says "a derivative of higher order than
supported" for any query found only inside `d/dt`, `δ` or `∫`. For `TimeIntegral(l)` and for an
input whose rate appears (`to_standard_form(… m ẍ + ḟ …, inputs=[f])` → raises with this text)
the diagnosis is misleading. Suggested: "'{query}' occurs only inside the calculus operator
'{expr}', which extract_linear_coeff treats as opaque; no linear coefficient exists".

### 5. LOW (maintainability) — two product-rule implementations with different normal forms

`geomech/core/math/ibp.py:87-112` (`_t_diff_product`, merges repeated factors → single
`(-2)` term) versus the scalar-target path at line 141 (`remainder.t_diff()`, the generic
`_BinaryMixin._apply_rule`, which does not merge). Same integrand, different output:

```
c l l δẋ  → δx(-1)c l̇ l + δx(-1)c l l̇          (two addends)
c l l Dot(ξ̇,r) → … + c l l̇ Dot(ξ,r)(-2)      (one addend)
```

Both are correct, but a stranger will ask why. Either route the scalar path through
`_t_diff_product` too (one-line change, and it would improve the F6.2 exhibits) or state in the
docstring why the vector path needs the merge and the scalar path does not. Relatedly,
`_remove_factor` (line 26) and `_split_dot_factor` (line 50) are two walkers over the same Mul
chain, one testing `str(x) == "1"` and the other `x is One`; pick one idiom.

### 6. LOW (maintainability) — private cross-module import

`geomech/dynamics/standard_form.py:123` imports `_has_algebraic` from `extract.py`. It is now
load-bearing for the partition contract (F3) and is used by two modules; make it public
(`has_algebraic`) and document it in the `extract.py` module docstring, or re-export it from
`geomech.core.math`.

### 7. LOW — new unary branches lack the query-free guard the binary branches have

`extract.py:262-268, 310-316`. `Hat`/`Transpose`/`Vee` always recurse, so a query-free operand
yields `Hat(0v)` / `(0v)'` / `Vee(0)` rather than the typed zero. `full_simplify` folds
`Hat(0v)` to `0` but leaves `(0v)'`. Harmless inside the pipeline (callers filter by
`_has_algebraic`) and visible only through the public API; a one-line
`if not _has_algebraic(expr.expr, query): return ZeroVector` per branch matches the rest of
the file.

### 8. LOW (process/documentation) — fix log is stale and one committed test was edited

- The fix log (`scalar-config-variables-review.md`, "xfail markers left in place (7)") says
  `TestElasticSphericalPendulum::test_xi_f_has_coriolis_term` keeps its marker. The working
  tree has the marker removed and the assertion changed from `isinstance(f, Add)` to
  `isinstance(f, VAdd)` with an "Arbitration" note (`tests/test_scalar_adversarial.py:421-428`),
  so the suite has 6 xfails, not 7. The edit is justified (the original assertion could never
  hold for a vector `f`), but that test was committed in 218f6b4 and the spec's definition of
  done says "No existing test modified"; the owner should see this explicitly.
- `TestR3xRCoupledMassMatrix::test_mass_matrix_symmetric` (line 506) has the same property (an
  assertion that can never hold because of the column/row representation) but was left as a
  strict xfail with reason "Representation: …". That is not a deferred quality item; it is a
  dead test. Apply the same arbitration (replace the third assertion with the representation-
  correct `str(M_x[l̈]) == str(M_l[ẍ])`, already in section 9) or delete it.

### 9. NIT — docstring signature

`extract.py:11` says `extract_linear_coeff(query, expr)`; the signature is `(expr, query)`.

## B. Maintainability summary

- Every new rule has a one-line math docstring (`_t_diff_product`, the `Hat`/`Transpose`/`Vee`/
  `VVMul` branches, the `collect` peel, `_apply_ibp`). Names are honest. No dead code found
  (`collect.py:96 return expr` is live for the both-sides case).
- Error messages name the offending node and query (`_zero_or_raise`, `_has_algebraic`,
  `_nonlinear`, `_check_variables`).
- No new public functions; `_check_variables` has type hints, the other new private helpers do
  not (consistent with the surrounding file).
- Duplication: findings 5 and 6.
- CHANGELOG: the `[Unreleased]` section does not mention the IBP fix, the `collect` change,
  the new `AlgebraicError` behaviour of `extract_linear_coeff`/`to_standard_form`, or the
  `compute_eom` validation. Recommended text below.

## C. Tests

- `tests/test_ibp_prefactor.py`: nine identities, each with the hand derivation in the
  docstring and an exact expected string; `test_single_pass_keeps_tree_shape` pins the in-place
  replacement, which is an implementation detail but the docstring says so. Not tautological:
  I reproduced each expected string from the identity before running the code.
- `tests/test_scalar_adversarial.py` section 9 (`TestFixF1PipelineExact`,
  `TestCollectTwoRotations`, `TestFixF2LinearUnaryOps`, `TestFixF3HigherDerivativeRule`,
  `TestFixF4TypedIdentity`, `TestFixF5Validation`): exact strings with derivations; the
  double-pendulum off-diagonal blocks match my hand derivation and numeric check.
- The 6 remaining xfails: `test_numeric_coefficients_have_one_type`,
  `test_f_combines_like_terms`, `test_zero_coefficients_are_pruned_from_M`,
  `test_xi_f_exact_structure`, `test_eta_equation_has_inertia_rate_term` are genuine F6
  quality items and labelled as such. `test_mass_matrix_symmetric` is not (finding 8).
- Weak existing test: finding 2.

## Recommended CHANGELOG entry (`[Unreleased]`)

```markdown
### Added
- Scalar configuration variables (`SystemVariables(scalars=[...])`) participate in
  `compute_eom` and `to_standard_form`. Spring-mass, spring-mass-damper, coupled spring-mass
  and the elastic spherical pendulum (R × S²) derive end to end.
- `extract_linear_coeff` replaces `extract_coeff` with a documented 3×3 dispatch grid
  (query type × expression type) and typed errors `NonLinearError` and `AlgebraicError`;
  scalar queries are supported inside `Hat`, `Transpose`, `Vee`, outer products and
  row-vector products.
- Integration by parts for scalar targets inside `Mul` chains.
- LaTeX rendering (`to_latex`, `display_eom`, `display_standard_form`) and example
  notebooks for the rigid pendulum, spherical pendulum, quadrotor and quadrotor with payload.

### Changed
- `Expr.__eq__` now returns `False` across different node types before comparing strings.
- `extract_linear_coeff` raises `AlgebraicError` instead of returning zero when the query
  occurs only inside a calculus operator (e.g. a jerk term, or `δẍ` from a Lagrangian in `ẍ`)
  or inside an unsupported node. `to_standard_form` partitions acceleration and input terms
  with the same test, so such terms can no longer be dropped silently.
- `compute_eom` validates `SystemVariables`: a constant, numeric or zero entry raises
  `ValueError`; an entry of the wrong expression type raises `TypeError`.

### Fixed
- Integration by parts dropped `d/dt` of a non-constant scalar factor multiplying
  `Dot(ξ̇, ·)`, so every system with a time-varying scalar in its kinetic energy (elastic
  pendulum, prismatic joint, prescribed `l(t)`) lost its Coriolis / rate-of-inertia terms
  (`−2 m l l̇ ω`, `−m l̈ e₁`, `J̇ Ω`). The prefactor is now differentiated with the product rule.
- `collect` could not move the variation out of `Dot(R₁ a, R₂ b)` with two different rotations,
  which dropped the cross-inertia coupling of multi-body SO(3) systems. The SO(3)×SO(3) double
  rigid pendulum now has `M₁[Ω̇₂] = m₂ l̂₁ R₁ᵀ R₂ ρ̂₂` and `M₂[Ω̇₁]` its transpose.
- `extract_linear_coeff` compared nodes by rendered string only, so `Mul(k, x)` could be
  mistaken for a leaf named `kx`.
```

## Required before merge

1. Finding 1 (guard the `collect` recursion) with the three-line regression test.
2. CHANGELOG entry as above.
3. Finding 8: update the fix log to the actual xfail set, and resolve
   `test_mass_matrix_symmetric` (arbitrate or delete). Owner to decide on
   `test_ibp.py::test_mul_recurses` (recommend: make `s` constant).

Findings 2–7 and 9 are recommended; 2 and 3 are cheap and close real gaps.
