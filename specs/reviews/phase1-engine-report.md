# Phase 1 engine report: items 1.1–1.5 of `specs/fix-plan.md`

Branch `fix/phase1-equations`, on top of 636e963 (`geomech.numeric`). Nothing committed.

| | Baseline (636e963) | After 1.1–1.5 |
|---|---|---|
| `python -m pytest tests/ -q` | 1375 passed, 5 xfailed | **1441 passed, 5 xfailed** |
| `ruff format --check` + `ruff check` | clean | clean |
| `mkdocs build --strict` | clean | clean |
| `engine_fuzz.py --quick` wrong-value findings | 5 (`collect`) + 1302/5000 at full size (`manifold_rules`) | **0** |
| Oracle systems (`tests/test_eom_oracle.py`) | 21 | 36 |

Line numbers below refer to the working tree after the change.

---

## 1.1 S² variation rule `δω = ξ̇ − 2 ω×ξ`

### What changed

- `geomech/core/transformations/manifold_rules.py:43` `_NEG1` → `_NEG2 = Scalar("(-2)", value=-2)`;
  `:86` the substituted node is `VAdd(TimeDerivative(ξ), SVMul(Cross(ω, ξ), _NEG2))`.
  Module docstring (`:1-40`) rewritten: cites Lee et al. eq. (5.13)
  `δω = −S(ω)ξ + (I − qqᵀ)ξ̇`, shows `q·ξ̇ = −ξ·q̇ = q·(ω×ξ)` from differentiating `ξ·q = 0`,
  hence `qqᵀξ̇ = ω×ξ` (because `ω×ξ ∥ q`), states that δω has the normal component
  `q·δω = −q·ξ̇` and names the symptom (gyroscopic terms at half strength).
- `geomech/core/base/expressions.py:312-331` `TSO3.delta(substitute=True)` and
  `TS2.delta(substitute=True)` deleted; both classes now inherit `Vector.delta()` →
  `Variation(self)`. grep: the `substitute=True` path was referenced only by
  `tests/test_manifold_rules.py::TestTangentDelta::test_ts2_delta_substitute` and
  `::test_tso3_delta_substitute`, which tested nothing but that path; both deleted. The
  `*_no_substitute` tests stay and pass.
- `docs/manifold-implementation-notes.md` rewritten: "needs verification" removed, the
  derivation above, the history/symptom paragraph, and the organisation of the substitution
  (the single source of truth is now `manifold_rules`). Rule lists updated for 1.2/1.3.
- `tests/test_manifold_rules.py::TestS2ManifoldRules::test_variation_of_omega_expands`
  (authorised): expects `(-2)`, docstring cites eq. (5.13). Module docstring line 4 updated
  to the −2 formula.
- `specs/reviews/scripts/engine_fuzz.py:1127-1129`: the `manifold_rules[tangential delta]`
  campaign, whose reference was the old projection formula, removed; the `[full delta]`
  campaign (exact jets) is the reference. Comment at `:322-326` updated.
- Paper: the accompanying paper (background and toolbox sections) already carries
  `\dot{\xi} - 2\,\omega \times \xi`; no change needed.

### Pinned tests (`tests/test_s2_variation_fix.py`, 9 tests, each `str(E)` + `assert_eom`)

All hand derivations use, for an S² variable with p = ∂L/∂ω, the coefficient of ξ in
`p·δω + (∂L/∂q)·(ξ×q) + δW` after IBP, i.e. on the tangent plane
`E ≡ −ṗ + 2 ω×p + q×∂L/∂q + f (mod q)`.

| | Lagrangian | hand result | engine string (= hand result) | `check_eom` |
|---|---|---|---|---|
| (a) | `½ml²ω·ω + c·ω − mgl q·e₃` | p = ml²ω + c → `−ml²ω̇ + 2ω×c − mgl q×e₃ + f` | `(d/dt(ω)(-1)mll+Cross(ω,c)(2)+Cross(q,e3)mgl(-1)+f)` | 6.8e-7 |
| (b) | `½ω·Aω − mgl q·e₃` | p = Aω → `−Aω̇ + 2ω×Aω − mgl q×e₃ + f` | `(A d/dt(ω)(-1.0)+Cross(ω,Aω)(2.0)+Cross(q,e3)mgl(-1)+f)` | 5.2e-7 |
| (c) | Lee et al. §5.6.2, ω-form and q̇-form | see below | 7-term strings pinned in `TestRotatingSphere` | 5.3e-7 / 4.1e-7 |
| (d) | `½ml²(|ω+W₀|² − (W₀·q)²) − mgl q·e₃` | p = ml²(ω+W₀) → `−ml²ω̇ + 2ml²ω×W₀ − ml²(W₀·q) q×W₀ − mgl q×e₃ + f` | `(d/dt(ω)(-1)mll+Cross(ω,W0)mll(2.0)+Cross(q,W0)mllDot(W0,q)(-1.0)+Cross(q,e3)mgl(-1)+f)` | 3.4e-7 |
| (e) | `½mẋ·ẋ + ½ml²ω·ω + ẋ·ω` | x: `−mẍ − ω̇ + F`; q: p = ml²ω + ẋ → `−ml²ω̇ − ẍ + 2ω×ẋ + f` | both pinned | 5.9e-7 |
| (f1) | `½mω·ω + ½Ω·JΩ + kω·Ω` | ξ: `−mω̇ − kΩ̇ + 2kω×Ω + f`; η: p_Ω = JΩ + kω, `E_η = −ṗ_Ω + p_Ω×Ω + M = −JΩ̇ − Ω×JΩ − kω̇ − kΩ×ω + M` | both pinned | 7.5e-7 |
| (f2) | `½mω·ω + ½Ω·JΩ + kω·(Rᵀe₃)` | ξ: p = mω + kRᵀe₃, ṗ = mω̇ − kΩ×Rᵀe₃ → `−mω̇ + kΩ×Rᵀe₃ + 2kω×Rᵀe₃ + f`; η: δ(Rᵀe₃) = (Rᵀe₃)×η → `−JΩ̇ − Ω×JΩ − k(Rᵀe₃)×ω + M` | both pinned | 2.7e-6 |

(c) Rotating sphere, `L̃ = ½mL²{|ω|² + 2qᵀS(Ωe₃)S(q)ω + Ω²qᵀCq}`, a = Ωe₃:
p = mL²(ω + q×(a×q)); ṗ = mL²(ω̇ − (ω×q)(q·a) − q((ω×q)·a)); 2ω×p = 2mL²(ω×a − (ω×q)(q·a));
q×∂L/∂q = mL²(−(q×ω)(q·a) + Ω² q×Cq). The `(ω×q)(q·a)` terms cancel:
`E ≡ −mL²ω̇ + 2mL²Ω ω×e₃ + mL²Ω² q×Cq + f (mod q)`. This equals `−q×(5.43) + f` on the tangent
plane since `q×q̈ = ω̇` and `q×PS(a)q̇ = −(ω×q)(q·a) ≡ −ω×a (mod q)`. The test evaluates the
engine's string with `evaluate` at 5 random tangent states (q unit, ω ⊥ q, ω̇ ⊥ q,
e₃ = (0,0,1), C = diag(1,1,0)) and compares the tangent components with both the compact
hand form and the book's (5.43) (`q̈ = ω̇×q + ω×(ω×q)`), tolerance 1e-9; the ω-form and the
q̇-form are also compared with each other. Note: here ∂L/∂ω ⊥ q, so (c) is insensitive to
the normal component of δω; it is a cross-check against the book, not a witness of the bug.
The six other rows are witnesses (the dynamics audit's rel. errors 0.05–1.2 before the fix).

All eight builders live in `tests/phase1_systems.py` and are in the oracle `SYSTEMS` list.

## 1.2 `collect`: peel the matrix, delete the false identity

### What changed (`geomech/core/math/collect.py`)

- `:70-75` top of `_collect_dot`: target on both sides → `NonLinearError` (fuzz finding 5;
  was `RecursionError`).
- `:84-104` `(MVMul(), Cross())`: the two correct triple-product rotations kept
  (`r.left == vec`, `r.right == vec`); the `(M·b)×c` rewrite deleted; `l.has(vec)` →
  `collect(Dot(l.right, MVMul(Transpose(l.left), r)), vec)`; target deeper in the Cross →
  rotate the triple product directly (flipping would flip straight back).
- `:127-138` `(MVMul(), _)`: `l.has(vec)` → same peel, replacing the one-level
  `inner.left == vec` / `inner.right == vec` special cases and the "return unchanged" dead end.
- Termination argument in the `_collect_dot` docstring (`:46-52`) and inline: every step
  either flips once into a shape that is handled without flipping back, or removes one
  matrix layer from the side holding the target; a side holding the target is never wrapped
  in a new `MVMul`. `tests/test_collect_termination.py` (4 tests) still passes unchanged.
- `tests/test_collect.py::test_mvmul_cross_vec_in_mvmul` (authorised): expects
  `Dot(v, MVMul(Transpose(M), Cross(w, u)))`, string `Dot(v,(M)'Cross(w,u))`, and numeric
  equality with `evaluate` under a shared env at three random points.

### Pinned tests (`tests/test_collect_fix.py`, 9 tests, each `str(E)` + `assert_eom`)

SO(3) hand form: p = ∂L/∂Ω, `E_η = −ṗ + p×Ω + (δR-terms) + M`, with
`δ(Rᵀe₃) = (Rᵀe₃)×η`, `d/dt(Rᵀe₃) = (Rᵀe₃)×Ω`, r := Rᵀe₃.

| Lagrangian | hand result | engine | `check_eom` |
|---|---|---|---|
| `½Ω·JΩ + (JΩ)·(r×b)` | p = JΩ + J(r×b); `−JΩ̇ − Ω×JΩ + J((Ω×r)×b) − Ω×J(r×b) − r×(b×JΩ) + M` | identical, pinned | 1.7e-6 (was 0.7) |
| `½Ω·JΩ + Ω·J(r×b)` (twin) | same Lagrangian | byte-identical string; numeric equality asserted | 7.2e-7 |
| `½Ω·JΩ + (R(JΩ))·e₃` | p = JΩ + Jr; `−JΩ̇ − Ω×JΩ + (JΩ)×r + J(Ω×r) − Ω×Jr + M` | identical, pinned | 1.0e-6 (crashed) |
| `((RJ)Ω)·e₃` grouping | same | byte-identical; numeric equality asserted | 6.5e-7 (crashed) |
| S² `½q̇·Aq̇ − mg q·e₃` | p = q×A(ω×q); `−(ω×q)×A(ω×q) − q×A(ω̇×q) − q×A(ω×(ω×q)) + 2ω×(q×A(ω×q)) − q×(ω×A(ω×q)) − mg q×e₃ + f` | identical term by term | 7.2e-7 (crashed) |
| SO(3)×SO(3) `kΩ₁·(R₁ᵀR₂Ω₂)` | Q = R₁ᵀR₂; E₁, E₂ derived in the docstring; net cross terms `+kΩ₁×QΩ₂` / `+kΩ₂×QᵀΩ₁` | pinned (engine keeps the identically-zero `Q(Ω×Ω)` terms, phase-5 item) | 1.6e-6 (crashed) |
| `½(Ω×ρ)·J(Ω×ρ)` | p = ρ×J(Ω×ρ); `−ρ×J(Ω̇×ρ) − Ω×(ρ×J(Ω×ρ)) + M` | identical | 6.0e-7 (crashed) |

Plus `TestTargetOnBothSides` (2 tests): the fuzz report's two `RecursionError` reproducers
raise `NonLinearError`. All seven builders are in the oracle `SYSTEMS` list.

## 1.3 `Hat(Vee(X)) → X` only for skew-symmetric X

- `geomech/core/transformations/simplify.py:168` the rule fires only if
  `_is_skew_symmetric(inner.expr)` (`:246-248`: a `Hat`, or `is_skew_symmetric`).
- `geomech/core/operations/geometry.py:82-85` `Hat.is_skew_symmetric` → `True`;
  `:145-152` `Transpose.is_symmetric` / `is_skew_symmetric` from the operand.
  `geomech/core/operations/addition.py:86-93` `MAdd.is_symmetric` / `is_skew_symmetric` (all
  nodes). `geomech/core/operations/multiplication.py:113-120` `SMMul.is_symmetric` /
  `is_skew_symmetric` (matrix operand). Leaves already read `ExprFlags`; operation nodes had
  no such property (the base `Expr._flag` returns `False`), so properties were added.
- **Latent bug found and fixed**: `geomech/core/base/expressions.py:462-470`
  `SkewSymmMatrix` appended `"SkewSymmetricMatrix"` to `attr` *after* `super().__init__`
  had built the flags, so no `SkewSymmMatrix` leaf ever had `is_skew_symmetric`. The attr
  is now appended before the call (same pattern as `SO3`). Without this, the authorised
  test change could not pass.
- `tests/test_simplify.py::test_hat_vee_cancel` (authorised) uses `SkewSymmMatrix`. New
  tests: `Hat(Vee(I))` not `I` and evaluates to zero; `Hat(Vee(J))` (J symmetric) unchanged
  and evaluates to a skew matrix ≠ J; plain `Matrix` unchanged; `Hat(Vee(Hat(v))) → Hat(v)`;
  propagation through `SMMul`/`MAdd`/`Transpose` (positive and negative cases). Note: with
  this package's `vee` (reads the (2,1),(0,2),(1,0) entries) `Hat(Vee(J))` for symmetric J
  is not zero, so the test asserts "skew and ≠ J", not zero.

## 1.4 `extract_linear_coeff` holes

`geomech/core/math/extract.py`:

- `:1-22` module docstring grid corrected: vector query / matrix expression raises; matrix
  query / scalar or vector expression raises; left-coefficient semantics stated.
  `paper_context.md:49-50` table rows corrected (`M·Hat(v)` → `AlgebraicError`; matrix query
  partial), and its code sketch says `NonLinearError`.
- `:153-157` dispatch fallthrough: `NotImplementedError` → `AlgebraicError` carrying the
  expression ("a MATRIX query has no left coefficient in a SCALAR expression | …").
- `:433-437` `Cross` branch: `left == right → ZeroMatrix` kept first (the existing
  `test_extract.py` pins `Cross(v, v) → 0`), then the both-sides guard, then the shortcuts.
- `:450-452` `MVMul` branch: both-sides guard → `NonLinearError`.
- `:603` `_vector_from_vector_svmul_scalar` fallthrough: `NotImplementedError` →
  `_zero_or_raise` (extra; same family).
- `:610-614` `_no_left_coefficient` helper; `:668-676` `_vector_from_matrix`: explicit
  `Hat() | Transpose()` case raises when the query is inside, `VVMul` raises only when the
  query is inside, fallthrough is `_zero_or_raise`.
- `:715-726` `_matrix_from_matrix` `SMMul`: both-sides → `NonLinearError`, scalar factor
  containing the query → `AlgebraicError`, fallthrough `_zero_or_raise`.

`tests/test_extract_fix.py` (27 tests): the fuzz reproducers (`Hat(x)`, `MMMul(M, Hat(q))`,
`Transpose(Hat(q))`, `Transpose(SMMul(VVMul(q, cᵀ), s))`, outer product → `AlgebraicError`;
`MVMul(Hat(q), q)`, `MVMul(q aᵀ, q)`, `MVMul(Hat(q), a×q)`, `MVMul(SMMul(M, q·a), q)`,
`Cross(Cross(a,q),q)`, `Cross(q, Cross(a,q))` → `NonLinearError`; `SMMul(Z, Dot(a, M b))`
query M → `AlgebraicError`; matrix query in scalar/vector expr → `AlgebraicError`), plus
seven linear shapes whose coefficient is verified numerically (`expr(q) == coeff·q`).

## 1.5 Zero outer product keeps its type

- `geomech/core/transformations/simplify.py:142-145` `VVMul` with a zero factor →
  `_zero_for(expr)`: `ZeroMatrix` for the outer product, `Zero` for the inner product.
- Tests (`tests/test_simplify.py::TestZeroOuterProductType`, 3): types of both products and
  `full_simplify(MAdd(VVMul(ZeroVector, Transpose(a)), M)) == M` (used to raise
  `ExpressionMismatchError`).

---

## Pinned equation strings in the existing suite

None changed. `git diff` of `tests/test_manifold_pipelines.py`, `tests/test_standard_form.py`,
`tests/test_vbl.py`, `tests/vbl_systems.py`, `tests/test_eom.py`, `tests/test_ibp.py`,
`tests/test_vbl_adversarial.py` is empty and the suite passes, matching the dynamics audit's
prediction (every pinned S² system has `∂L/∂ω ⊥ q`).

One unpinned printed equation of a shipped system did change (review F4):
`vbl_systems.quadrotor_with_load`, load equation, the term parallel to `q`,
`Cross(ω, Cross(q, ẋ_L)) m_Q l`, changed coefficient from −1 to −2. It is
`q (ω · ẋ_L)`, invisible on the tangent plane, so the oracle, `test_vbl.py` and the
notebook are unaffected and the dynamics do not change.

## Existing tests changed

Authorised by the plan (3): `test_manifold_rules.py::test_variation_of_omega_expands`,
`test_collect.py::test_mvmul_cross_vec_in_mvmul`, `test_simplify.py::test_hat_vee_cancel`.

Required by the specified behaviour but not listed in the plan (4) — **please confirm**:

- `tests/test_extract.py::TestNotImplementedErrors::test_vector_from_matrix_mmmul_left`,
  `::test_vector_from_matrix_mmmul_right_nonleaf`, `::test_vector_from_matrix_smmul_mat_side`
  asserted the silent zero (`MMMul(0, M)`, `SMMul(0, s)` for `Hat(v)` inside) that 1.4 replaces
  with `AlgebraicError`; they now expect the error (docstrings say why).
- `tests/test_scalar_adversarial.py::TestNonLinearAcceleration::test_matrix_query_raises_typed`
  asserted `NotImplementedError`, which 1.4 replaces with `AlgebraicError`.

Deleted (2): `test_ts2_delta_substitute`, `test_tso3_delta_substitute` (tested only the deleted
code path). Docstring-only edits: `tests/test_manifold_rules.py` module docstring and the hand
derivation header of `tests/test_scalar_adversarial.py::TestElasticSphericalPendulum` (stated
the old δω; its result is unaffected because `ω·(ω×ξ) = 0`).

## Test counts

Baseline 1375 → 1441 passed (+66 net): +9 `test_s2_variation_fix.py`, +9 `test_collect_fix.py`,
+27 `test_extract_fix.py`, +15 oracle systems (21 → 36), +8 in `test_simplify.py`, −2 deleted
in `test_manifold_rules.py`. 5 xfailed unchanged.

## Fuzz (`python specs/reviews/scripts/engine_fuzz.py --quick`, seed 20261113)

**0 wrong-value findings.** Before the harness update the only wrong-value cluster was the
`manifold_rules[tangential delta]` campaign (5/200), whose reference was the pre-fix
projection formula; that campaign was removed (the exact `[full delta]` campaign is
200/200 value-preserving, up from 3698/5000 at full size). Remaining 38 distinct findings,
all informational per the fuzz report's own classification: 14 `flag-dependent`
(adversarial mode, rules firing on declared `UnitNorm`/`Orthogonal`/`Symmetric` flags, fuzz
finding 12), 3 `documented-error-on-legal-input` (2 row-vector `ExpressionMismatchError`
cases of fuzz finding 6, phase 3.2; 1 is the new `NonLinearError` for a target on both sides
of a `Dot`, the intended replacement of the `RecursionError`), 21 `unsupported-linear-term`
(`AlgebraicError` on `Transpose`/`Vee` around the query, matrix queries, and the vector-in-
matrix cell that now raises by design; fuzz finding 13). `collect[target occurs once]`
167/167, `simplify`/`full_simplify` 200/200, `extract-nonlinear` 200/200 raised
`NonLinearError`.

## Not done / out of scope

- Phase 5 items that the new pinned strings make visible: `Q(Ω×Ω)` is kept inside an `MVMul`
  (vector rules do not descend into `MVMul`, eng-F8), and the three `k Ω₁×QΩ₂`-type terms in
  the SO(3)×SO(3) equation are not merged (commutative like-term limitation).
- `Transpose(Transpose(v))` / `Vee(Hat(v))` around a query still raise in
  `_vector_from_vector` (documented; `full_simplify` removes them before the pipeline extracts).
- The accompanying paper needed no edit (formula already −2).
