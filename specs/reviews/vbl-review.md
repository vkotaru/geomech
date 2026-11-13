# Adversarial review: variation-based linearization (`geomech.dynamics.vbl`)

Reviewed: commit `a0497d2` on `feature/vbl` against `specs/vbl.md`. Nothing in the diff was
modified; the review adds `tests/test_vbl_adversarial.py` (51 passing tests that pin the
attack surface, 7 strict `xfail` tests that pin the findings) and this document.

Suite after the review: `python -m pytest tests/ -q` → **1285 passed, 12 xfailed**
(1234 + 51 new passing; 5 pre-existing + 7 new strict xfails). `ruff format --check` and
`ruff check` clean.

## Independent derivation (done before reading `vbl.py`)

- **S² kinematic row.** With q̇ = ω × q, δq = ξ × q, ω ⊥ q, ξ ⊥ q: equating δ(q̇) and
  d/dt(δq) and using Jacobi gives (δω − ξ̇ − ξ × ω) × q = 0, so δω − ξ̇ − ξ × ω = c q.
  Dotting with q and using q·δω = q·(ξ × ω) (from δ(q·ω) = 0) and q·ξ̇ = q·(ω × ξ)
  (from d/dt(q·ξ) = 0) gives c = q·(ξ × ω); the full relation is δω = ξ̇ + 2 ξ × ω on the
  constraint set. Projecting: ξ̇ = q qᵀ ω̂ ξ + (I − q qᵀ) δω, **identical to the spec**. The
  row preserves ξ ⊥ q for any δω (d/dt(q·ξ) = q·(ξ×ω) + q·(ω×ξ) = 0). Note the EOM
  manifold rule δω = ξ̇ − ω × ξ (`manifold_rules.py`) drops the component along q; that is
  harmless for the EOM pipeline (it is only ever dotted with tangent vectors) and VBL does
  not use it.
- **SO(3) kinematic row.** δΩ = η̇ + Ω × η ⇒ η̇ = −Ω̂ η + δΩ. Matches.
- **Sign convention.** E is the post-IBP coefficient (negative Euler–Lagrange), so
  δE = −M δv̇ + A δx + B δu with M = −coeff(δv̇). Free rigid body: M = J,
  A[δΩ] = Hat(JΩ) − Hat(Ω) J, B = I. Rigid pendulum: gravity moment in the body frame is
  ρ × Rᵀ(−m g e₃), so E = −JΩ̇ − Ω×JΩ − m g ρ×Rᵀe₃ + M and A[η] = −m g Hat(ρ) Hat(Rᵀe₃)
  (δ(Rᵀe₃) = Hat(Rᵀe₃) η). Spherical pendulum: A[ξ] = −m g l Hat(e₃) Hat(q); at q = −e₃
  this is −m g l on the tangent plane (pendulum stiffness, q = −e₃ is the minimum of
  U = m g l q·e₃). Quadrotor: A[η] = −f R Hat(e₃), B[δf] = R e₃. Elastic spherical
  pendulum (R × S²): full block table in `TestElasticSphericalPendulumAnalytic`.
  **Every value pinned in `tests/test_vbl.py` agrees with these derivations.** No pinned
  test asserts a wrong block.
- Both trajectory families in `TestKinematicRowsIndependent` (exact q̇, Ω, d/dt δq, d/dt δR;
  only ∂_ε is differenced) reproduce the two kinematic rows to 1e-7.

## Findings (ordered by severity)

### F1 — Silent wrong output: a term of δE carrying two perturbation symbols is assigned to one of them and the other contribution is dropped

- **What.** `expand` does not distribute `VVMul` (outer product `v wᵀ`, inner product
  `vᵀ w`) over `VAdd`. After δ, a factor such as `(q + ω)` becomes `(ξ × q + δω)`, so one
  additive term of `expand(δE)` contains both ξ and δω. `_partition_terms` assigns the term
  to the first symbol of the priority list that is found (δω), and the `VAdd` branch of
  `extract_linear_coeff` keeps only the addends containing the query, so the ξ-contribution
  disappears. `_check_coefficient_free` cannot see it (the dropped addend never reaches the
  coefficient). Result: `A[ξ] = 0` with no error.
- **Why it is wrong.** The spec's completeness rule requires *exactly one* perturbation
  symbol per term and an `AlgebraicError` otherwise; the implementation checks for *at
  least one*. The output is mathematically wrong and the finite-difference check fails on
  it.
- **Minimal repro.** `E = MVMul(VVMul(q + ω, Transpose(e3)), e3)` (= (e₃·e₃)(q + ω)) for an
  S² variable q gives `A[ξ] = 0`; the correct block is −(e₃·e₃) Hat(q). Scalar-row variant:
  `E_s = VVMul(Transpose(q + ω), e3)` gives `A[δs][ξ] = 0` instead of (q × e₃)ᵀ.
  `compute_eom` itself emits such outer products: with a position-dependent inertia
  `L = ½ m ẋ·(p pᵀ) ẋ`, `p = x − l q`, the derived E contains `(x + q l(−1))((ẋ + …))ᵀ ẋ`
  terms; there `linearize` happens to raise on an unrelated `Transpose(MAdd)` first, which
  is luck rather than protection.
- **Tests.** `TestF1OuterProductOfSumsDropsATerm::{test_vector_row_left_factor,
  test_vector_row_right_factor, test_scalar_row_inner_product}` (xfail strict). They accept
  either the correct block or an `AlgebraicError`.
- **Fix.** In `_partition_terms`, after finding the owner, assert that no *other*
  perturbation symbol occurs in the term (`_has_algebraic(term, other)` for every other
  symbol, treating the calculus-boundary error as "present") and raise `AlgebraicError`
  naming the term; that implements the spec's rule and turns the silent zero into an
  error. Independently, teach `expand` to distribute `VVMul` over `VAdd` in both factors
  (and through `Transpose(VAdd)`), which makes these cases linearize correctly.

### F2 — A tangent vector (Ω, ω) is accepted as an input and the A[δΩ] block is silently moved to B

- **What.** `_check_inputs` rejects only configuration variables. `Ω = R.get_tangent_vector()`
  (a `TSO3`, hence a `Vector` leaf) passes, so δΩ is both a state and an input column. In the
  dynamic loop the owner is the state symbol, but `elif sym in input_perts` compares by
  `__eq__` (= `str`), is true, and the textbook block Hat(JΩ) − Hat(Ω) J is written to
  `B[δΩ]` while `A[δΩ]` stays zero. Same for `ω` on S².
- **Tests.** `TestF2TangentVectorAcceptedAsInput::test_so3_tangent_vector_as_input_raises`,
  `::test_s2_tangent_vector_as_input_raises` (xfail strict);
  `::test_so3_tangent_vector_as_input_currently_moves_the_block` documents the present
  behaviour.
- **Fix.** In `_check_inputs`, reject any input whose `delta()` prints as one of the state
  or dstate symbols (equivalently reject `TS2`/`TSO3` leaves and any `u` with
  `str(u.delta())` in the state set), and use identity on the symbol list (`is`/index)
  rather than `in` when deciding the destination table.

### F3 — Duplicate inputs accepted; `B` is no longer dense over `inputs`

- `inputs=[M, M]` yields `lin.inputs == [δM, δM]` but every `B[row]` has a single column, so
  `len(B[row]) != len(lin.inputs)` and iteration over `inputs` double-counts in any consumer.
- **Test.** `TestF3DuplicateInputs::test_duplicate_input_raises` (xfail strict).
- **Fix.** Raise `ValueError` on repeated `str(u)` in `_check_inputs`.

### F4 — Equation type not checked against the variation vector; wrongly shaped block emitted

- `{str(η): (η, Dot(Ω, Ω))}` is accepted. `row_type` is taken from η (vector) but the
  variation of the scalar equation produces a 3×1 `Vector` coefficient that is stored in a
  slot `block_shape` reports as 3×3. Spec: "Dense blocks of correct shape".
- **Tests.** `TestF4EquationTypeMismatch::test_scalar_equation_for_vector_variation_raises`
  (xfail strict); `::test_mismatch_currently_yields_a_wrongly_shaped_block` documents it.
- **Fix.** In `_check_eom`, require `eq.type == var_vec.type` (and `eq.type` in
  {SCALAR, VECTOR}).

No other correctness defects were found. F1 is the only finding that affects the output of
a system a user could plausibly build from a Lagrangian; F2–F4 are input-validation gaps
that produce wrong structures instead of errors.

## What was attacked and held

Every item below is pinned by a passing test in `tests/test_vbl_adversarial.py`.

- **Kinematic rows against the truth** (not against the spec formula): S² and SO(3) rows
  reproduce ξ̇ and η̇ of exact perturbed trajectory families with ω ≠ 0, Ω ≠ 0
  (`TestKinematicRowsIndependent`, 8 cases).
- **R × S² with scalar and vector inputs** (elastic spherical pendulum): all 20 blocks of the
  two dynamic rows match the analytic linearization with an independent evaluator,
  including the 1×3, 3×1 and 1×1 zero blocks and shapes
  (`TestElasticSphericalPendulumAnalytic`).
- **q-dependent mass matrix** (q̇-form spherical pendulum): M[δω̇] = m l² (I − q qᵀ), and
  A[ξ] carries the variation of the mass term, δ(−m l² q × (ω̇ × q)), exactly
  (`TestSphericalPendulumQdotForm`).
- **δω is not treated as tangent**: E = c q × (ω × q) linearizes to c (I − q qᵀ) δω + … and
  the row sum is c δω on tangent directions; the S² rules in `vector_rules.py` do not
  recurse into `Variation`, so no block depends on a false δω ⊥ q
  (`TestVariationOfOmegaIsNotTangent`).
- **Terms vanishing on the constraint set but with non-zero differential**
  (`q (q·ω)`): correct (`dot_q_omega_times_q`).
- **Damping / velocity-dependent forces** go to A[δṡ] = −c and A[δẋ] = −c I, not to M.
- **Nonlinear inputs and velocities** (`F (F·F)`, `Ω (Ω·Ω)`): B and A blocks include the
  rank-one `2 F Fᵀ`, `2 Ω Ωᵀ` terms.
- **Quadrotor with load**: thrust blocks `B[x_L][δf] = R e₃`, `B[ξ][δf] = −l q × R e₃`, the
  cross-inertia `M[x_L][δω̇] = m_Q l Hat(q)` and the attitude couplings `A[x_L][η] = −f R
  Hat(e₃)`, `A[ξ][η] = l f Hat(q) R Hat(e₃)` match hand values.
- **Sums inside products**: 13 vector-row and 2 scalar-row equations with `VAdd`/`MAdd`/`Add`
  nested under `SVMul`, `MVMul`, `MMMul`, `SMMul`, `Hat`, `Cross`, `Dot`, `Transpose(Hat(·))`
  are fully distributed by `expand` and pass finite differences. Only `VVMul` is missed
  (F1).
- **Error paths**: stray η or `Variation(Ω)` in E → `ValueError` from the engine (δ of a
  variation is unrepresentable); `d/dt(F)` with F an input → `AlgebraicError` naming the
  calculus node; `Transpose(MAdd)`, `Transpose(VVMul)`, `Vee(MAdd)`, `Vee(Hat(VAdd))` →
  `AlgebraicError` (opaque to extraction, but never silent); variation vector as input →
  `ValueError`; input named like a configuration variable → `ValueError`; eom keys from a
  different `SystemVariables`, non-dict eom, non-pair values → errors; non-SO3 matrix
  configuration → `TypeError`.
- **Exports**: `geomech.linearize`, `geomech.LinearizedSystem`, `geomech.dynamics.__all__`.
- **Docs page**: the printed free-rigid-body example reproduces the actual `str(lin)`
  output; the `to_standard_form` M = −J claim is true; the quadrotor and elastic-pendulum
  block strings quoted in "Reading the blocks" are the actual ones.

## Is the finite-difference evidence trustworthy?

**Yes, for the dynamic rows**, with two caveats.

Audited points of `tests/numeric_helpers.py`:

- The perturbation really stays on the manifold: q(ε) = exp(ε ξ̂) q, ω(ε) = exp(ε ξ̂)(ω + ε w)
  with ξ, w ⊥ q keeps |q| = 1 and q·ω = 0 exactly; R(ε) = R exp(ε η̂) matches δR = R Hat(η);
  the induced directions δω = ξ × ω + w, δω̇ = ξ × ω̇ + a are the ones fed to the blocks and
  satisfy q·δω = q·(ξ × ω) (checked in `TestKinematicRowsIndependent`).
- `Evaluator` keys every leaf and every `TimeDerivative` chain by its printed name, draws
  constants once on the base evaluator and copies the environment to the ±ε evaluators, so
  E(±ε) and the blocks see identical parameters; `Variation` symbols must be supplied
  (never random); zero constants are resolved by type; operation nodes are not `Vector`
  subclasses, so a surviving `d/dt(non-leaf)` raises instead of becoming a random leaf.
- A wrong block cannot pass: both sides are O(1) for every system, the tolerance is
  1e-6 relative to max(1, |rhs|), observed errors are 1e-10–1e-12, and
  `test_finite_difference_check_detects_a_wrong_block` plus the F1 repros (the check fails
  on the dropped-term blocks) show it is discriminating. The same check fails if the S²
  coupling δω = ξ × ω + w is replaced by w alone.

Caveats:

1. **The kinematic rows are not finite-differenced.** `_dstate_values` evaluates ξ̇ and η̇
   from the very formulas `linearize` writes, so a wrong kinematic row would pass
   `finite_difference_check` and also `test_kinematic_row_numeric` in `tests/test_vbl.py`
   (same formula on both sides). The independent trajectory tests in this review close that
   gap.
2. Directions are tangent to the constraints only, so the blocks are verified as linear maps
   on the constraint set. That is the correct notion (δE depends on ξ only through ξ × q and
   on δω only through its constrained values), but a block that differs off the tangent
   space would not be detected; none does, by construction of the extraction.

## Spec compliance

| Spec item | Status |
|---|---|
| `linearize(eom, variables, inputs=None)`; `_check_variables` reused | Yes |
| `LinearizedSystem` frozen dataclass with `rows`, `states`, `dstates`, `inputs`, `M`, `A`, `B`, `block_shape`, `__str__` | Yes (plus `row_types`, `col_types`) |
| Row order: kinematic rows (`"kin:" + str(var)`) then eom keys | Yes |
| States: configuration perturbations then velocity perturbations, `SystemVariables` order | Yes |
| Dense blocks, shapes by row × column type, zero constants of right shape | Yes for valid input; the 1×3 zero is `Transpose(ZeroVector)` (no library constant exists); violated under F3/F4 |
| Kinematic rows S², SO(3), R, R³ as specified | Yes, verified independently |
| Sign convention δE = −M δv̇ + A δx + B δu | Yes |
| Completeness: exactly one perturbation symbol per term, `AlgebraicError` otherwise | **Partially**: "none" and "inside calculus node" raise; "more than one" is not checked (F1) |
| `apply_manifold_rules` not called | Yes |
| Blocks through `full_simplify`, no `(-1)(-1)` chains | Yes (no engine change needed) |
| Systems 1–9 with hand derivations in docstrings; FD check for all | Yes; derivations verified here |
| Scalar inputs | Yes (quadrotor, elastic pendulum) |
| `tests/numeric_helpers.py` with `hat`, `vee`, `rand_rot`, `expm_so3`, `Evaluator`; review script untouched | Yes |
| Exports in `geomech/dynamics/__init__.py` and `geomech/__init__.py` | Yes |
| `docs/vbl.md` content and `mkdocs.yml` nav | Yes; claims checked, none false |
| CHANGELOG entry; README first line | Yes |
| No existing test modified; suite green; ruff clean | Yes |
| Implementation report | Yes; accurate (it correctly reports no engine change and no xfail for system 9) |

## Code-quality notes (lower severity)

- **String keys everywhere.** `M[row][str(col)]`, `_check_eom`, `_check_inputs` and the
  `sym in input_perts` test all rely on `str()` equality (`Expr.__eq__`). Two different
  `SystemVariables` with the same printed names are indistinguishable (an eom from
  `SO3("R")` is accepted with another `SO3("R")`), and F2 is a direct consequence. Use
  identity (`is`, or an index into `priority`) when routing a symbol to M/A/B.
- **Frozen dataclass with mutable dict fields.** `LinearizedSystem` is `frozen=True` but
  `M`, `A`, `B`, `rows`, `states` are plain mutable containers shared by reference;
  `tests/test_vbl.py` already builds a "bad" system by copying and mutating them. Consider
  `MappingProxyType`/tuples or document that the object is not immutable.
- **Priority-order partitioning is fragile.** Correctness of `_partition_terms` depends on
  (i) `expand` fully distributing products over sums and (ii) every term being linear. Only
  (ii) is asserted (via `_check_coefficient_free`), and only when the second symbol
  survives into the coefficient. An explicit exactly-one check (F1 fix) removes the
  dependence on (i).
- **`_check_inputs` rejects `S2` but not `TS2`/`TSO3`/`SO3`-derived leaves** (F2). The
  intent "a plain `Scalar`/`Vector` leaf" should be an allow-list on exact types.
- **Equation type unchecked** (F4); also `eq` may be any `Expr`, including a `Matrix`,
  which would reach `_zero_block` and raise an opaque `TypeError` from the `match`.
- **`finite_difference_check` populates constants by evaluating `eq`** before sampling the
  direction; a block leaf absent from `eq` (impossible today) would be drawn independently
  on the base evaluator and never appear in E(±ε), silently decoupling the two sides. A
  defensive check that every leaf of every block was seen in some `eq` would make the
  helper robust against future changes.
- **Error type for stray perturbation symbols in E** is the engine's `ValueError`
  (δ of a variation), not `AlgebraicError`. Acceptable, but the docstring of `linearize`
  lists only `TypeError`, `ValueError` for invalid *arguments*; a stray η in an equation
  is a different failure and could be caught and re-raised with the row label.
- `tests/test_vbl.py::TestSphericalPendulumS2::test_kinematic_row_numeric` compares the row
  with the same formula it was written from (not independent); keep it, but the
  trajectory tests in `tests/test_vbl_adversarial.py` are the actual evidence.
- Transpose of sums and outer products, and `Vee`, are opaque to `extract_linear_coeff`
  (raise rather than extract). Not a VBL bug, but a limitation worth a note in
  `docs/vbl.md`: equations containing `Transpose(MAdd)` or `Vee` cannot be linearized yet.

## Fix log

Applied on `feature/vbl` after this review. No test in this file's companion
`tests/test_vbl_adversarial.py` was weakened; the seven `xfail(strict=True)` markers on the
F1–F4 tests were removed once the tests passed. No other existing test file was edited.

### F1 — engine: `expand` distributes `VVMul` over `VAdd`

- `geomech/core/transformations/expand.py`: new `case VVMul()` handled by `_expand_vvmul`.
  Outer product `v wᵀ`: a `VAdd` in `v` (after expanding `v`) or a `VAdd` inside
  `Transpose(w)` (after expanding the inner of the transpose) is distributed to an `MAdd`
  of outer products. Inner product `vᵀ w`: a `VAdd` inside `Transpose(v)` or in `w` is
  distributed to an `Add` of inner products. Each addend of a `Transpose(VAdd)` becomes
  `Transpose(addend)`, so every produced `VVMul` keeps the (plain, Transpose) shape the
  constructor requires. A `Transpose` factor whose inner does not expand to a `VAdd` is kept
  as the original node. `Transpose` is still not an `expand` case of its own, so
  `Transpose(MAdd)`, `Transpose(VVMul)` and `Transpose(Hat(VAdd))` are unchanged and keep
  raising in `extract_linear_coeff` (`TestErrorPathsHeld::test_unsupported_structures_raise_not_drop`
  and `transpose_hat_vadd` still hold).
- Termination: `expand` is a single bottom-up pass; the new rule only recurses into the
  operands of the node it is handling and builds the result from already-expanded
  subterms without calling `expand` on its output, so it inherits the structural
  recursion's termination. `tests/test_expand_vvmul.py` checks the rule, the opaque
  `Transpose` cases, and that a second pass is a fixed point.
- The three F1 repros now linearize to the correct blocks (`−(e₃·e₃) Hat(q)`,
  `e₃ (q × e₃)ᵀ`, `(q × e₃)ᵀ`), not to an `AlgebraicError`.

### F1 — defence in depth: exactly one perturbation symbol per term

- `geomech/dynamics/vbl.py::_partition_terms`: after the owner is found, every later
  symbol in the priority list is checked with the new `_occurs_outside(term, sym, owner)`,
  which walks the term skipping subtrees identical to the owner (so `δω` inside the owner
  `d/dt(δω)` is not a second symbol) and treats a foreign calculus node containing the
  symbol as an occurrence. A hit raises `AlgebraicError` naming the row, the term and both
  symbols. Buckets are now a list indexed by position in `priority`, not a dict keyed by
  `Expr` (whose hash is the printed name). Covered by `tests/test_vbl_validation.py`.

### F2 — tangent vectors as inputs

- `vbl.py::_check_inputs` now receives the state and d/dt-of-state symbols and raises
  `ValueError` for `TS2`/`TSO3` leaves and for any input whose `str(u.delta())` is already
  a state symbol. (`ξ`/`η` as inputs, previously caught by the configuration-name check for
  S² only, are caught by the same rule.) The dynamic-row loop routes each bucket to
  `M`/`A`/`B` by its index range in `priority`, not by `sym in input_perts`.

### F3 — duplicate inputs

- `vbl.py::_check_inputs` raises `ValueError` on a repeated `str(u)`.

### F4 — equation type vs variation vector

- `vbl.py::_check_eom` raises `TypeError` when `eq.type` is not scalar or vector, or
  differs from `var_vec.type`. Validation order is now `_check_variables` →
  `_perturbation_variables` (so a non-SO3 matrix still fails with the "SO3" message) →
  `_check_eom` → `_check_inputs`.

### Quality notes

- Addressed: identity-based routing of blocks; `linearize` docstring lists the new error
  conditions; `expand` module docstring documents the rule and why `Transpose` stays opaque.
- Deliberately left: a stray variation symbol in `E` still surfaces as the engine's
  `ValueError` ("variation vector" / "variation of a variation"). Two held tests in
  `tests/test_vbl_adversarial.py` pin exactly that `ValueError` and message, and
  `AlgebraicError` does not subclass `ValueError`, so re-raising would break them. The
  `linearize` docstring now says so. The frozen-dataclass-with-mutable-fields note and the
  `finite_difference_check` leaf-coverage note are unchanged.
- The two non-xfail tests that "document the present (wrong) behaviour"
  (`TestF2TangentVectorAcceptedAsInput::test_so3_tangent_vector_as_input_currently_moves_the_block`,
  `TestF4EquationTypeMismatch::test_mismatch_currently_yields_a_wrongly_shaped_block`)
  assert the bug and therefore fail once F2 and F4 are fixed; they were left for a human
  to remove, per the fixer's brief.
