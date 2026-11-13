# Final review: variation-based linearization (branch `feature/vbl`)

Spec: `specs/vbl.md`. Adversarial review + fix log: `specs/reviews/vbl-review.md`.
Implementer's report: `specs/reviews/vbl-implementation-report.md`.
Scope reviewed: `git diff db86568 HEAD` (commits `a0497d2`, `293fab6`, `1d2261c`):
`geomech/dynamics/vbl.py`, `geomech/core/transformations/expand.py`, exports,
`tests/{test_vbl,test_vbl_adversarial,test_vbl_validation,test_expand_vvmul}.py`,
`tests/{vbl_systems,numeric_helpers}.py`, `docs/vbl.md`, `mkdocs.yml`, `CHANGELOG.md`,
`README.md`.

State at review (HEAD `1d2261c`): `python -m pytest tests/ -q` → **1322 passed, 5 xfailed**
(the 5 are the pre-existing strict xfails in `tests/test_scalar_adversarial.py`);
`ruff format --check` and `ruff check` clean; `mkdocs build --strict` succeeds (the
implementation report says mkdocs was not installed; it is now, and the page builds).
Same counts after the edits made in this review (docs and CHANGELOG only).

## Verdict: GO WITH CHANGES

The mathematics is right. I re-derived both kinematic rows and the sign convention from
scratch, checked two pinned systems (rigid pendulum on SO(3), spherical pendulum on S²)
against my own closed forms with my own evaluator, and finite-differenced the quadrotor-
with-load translational row with code that shares nothing with `tests/numeric_helpers.py`.
Every block agrees. The `expand` change is correct for both legal `VVMul` shapes, keeps
`expand` a single bottom-up pass with a fixed point, and leaves the printed output of every
pinned `compute_eom`/`to_standard_form` system byte-identical (verified by running the
pre-feature tree side by side). No test was weakened. The code is readable and raises
instead of returning silent zeros.

The "changes" are small and none touches code logic: (1) the documentation and CHANGELOG
edits I applied in the working tree (see section D; they must be committed with the
feature), and (2) two comment-only clean-ups in test files that I was not permitted to
make. One mathematical subtlety that the docs did not mention (S² blocks are unique only on
the constraint set) is now documented; it is not a defect.

## A. Findings

Severity: **required** = before merge; **recommended** = soon, no behaviour change;
**optional** = judgement call.

### A1 (required, comment-only) Test comments still narrate the xfail history

- `tests/test_vbl_adversarial.py:15-16` ("Section 8 holds the tests that pinned findings
  F1-F4 …; they were strict xfails until the fixes landed and now run as ordinary tests.")
  and `tests/test_vbl_adversarial.py:639` ("# 8. Findings F1-F4 (fixed; formerly strict
  xfail)").
- Why: the brief for this review asks that no reference to xfail or to the fixing process
  remain in test comments; in a year the history is noise and `specs/reviews/vbl-review.md`
  already records it.
- Change: line 15-16 → "Section 8 covers outer and inner products with a sum as a factor,
  tangent vectors offered as inputs, duplicate inputs and equation/variation type
  mismatches." Line 639 → "# 8. Outer products of sums, tangent-vector inputs, duplicate
  inputs, type mismatches". The class docstrings in section 8 (lines 644-649, 693-697,
  715-717, 727-730) still say "F1: …", "Must raise"; keep the math, drop the finding labels.

### A2 (required) Commit the documentation and CHANGELOG edits made in this review

- `docs/vbl.md` and `CHANGELOG.md` are modified in the working tree (section D). They add
  the limitations the adversarial review asked for, the constraint-set note (A3), fix the
  `M[δẋ]` indexing in the R/R³ kinematic bullet, and move the `expand` entry to *Changed*.
  Review them and commit with the feature; nothing else in the tree is modified.

### A3 (recommended, documented) S² blocks are unique only on the constraint set

- `geomech/dynamics/vbl.py:272` takes δ of the equation *as simplified by `compute_eom`*,
  which has already used ω·q = 0 (e.g. ω × (ω × q) → −|ω|² q). The variation of the
  simplified and unsimplified forms differ by terms ∝ (q·δω − ξ·(ω × q)), which vanish on
  the constraint set δ(q·ω) = 0 but not as raw matrices.
- Numbers (quadrotor with load, row δx_L, 5 random states): my hand blocks from
  q̈ = ω̇ × q + ω × (ω × q) differ from the library's `A[ξ]`, `A[δω]` by up to **3.47**
  entrywise, yet `(A_lib − A_hand)·(ξ, ξ×ω + w)` on tangent directions is **1.6e-15**, and
  my own central difference of E agrees with my hand linearization to **1.6e-9** relative.
  So both are correct VBL; the printed matrices are not canonical off the constraint set.
- Why it matters: a user feeding the 3-dimensional (ξ, δω) blocks to an LQR without
  projecting gets a system whose off-constraint part is an artefact of the simplification
  order. This is inherent to VBL (Wu & Sreenath work on the tangent space), not a bug.
- Change: documented in `docs/vbl.md` ("Blocks on S² are defined on the constraint set").
  Optional follow-up: a `project_tangent(lin)` helper or a note in `LinearizedSystem`'s
  docstring. No code change required.

### A4 (recommended) `vbl.py` imports private names from `extract.py`

- `geomech/dynamics/vbl.py:98`: `from geomech.core.math.extract import _has_algebraic,
  _same, extract_linear_coeff`.
- Why: a rename inside `extract.py` breaks `linearize` with an `ImportError` at import
  time of `geomech` itself (`geomech/__init__.py` imports `linearize`). Cross-module use of
  underscore names also hides that these two are now part of the engine's contract.
- Change: expose them as `has_algebraic` / `same_node` (keep the underscore aliases for one
  release) or move them to `geomech/core/math/_tree.py` with a docstring saying `vbl.py`
  depends on them.

### A5 (recommended) `_negate` builds nodes three different ways

- `geomech/dynamics/vbl.py:513-522`: scalar → `Mul(expr, -1)`, vector → `expr * -1`
  (operator dispatch), matrix → `SMMul(expr, -1)`.
- Change: `SVMul(expr, -1)` for the vector case so all three are explicit constructors;
  or all three via the operator. Behaviour is identical today.

### A6 (recommended) Review-process labels in test docstrings

- `tests/test_vbl_validation.py:1-8` ("Defence-in-depth checks added … after the
  adversarial review of VBL … (F1 defence …)"), `tests/test_expand_vvmul.py:68` ("the F1
  shape from the VBL review").
- Change: name the mathematics, not the finding: "an outer product whose left factor is a
  sum, under an enclosing matrix–vector product" and "`_partition_terms` rejects a term
  carrying two perturbation symbols even when `expand` could not split it". The tests
  themselves are exact and fine.

### A7 (optional) Trivial tidy-ups in `vbl.py`

- `vbl.py:260`: `priority + [d for d in dstates[: len(cfg_perts)]]` → `priority +
  dstates[: len(cfg_perts)]`.
- `vbl.py:141-149`: dataclass fields are typed `list`/`dict`; the docstring carries the
  element types. `list[str]`, `list[Expr]`, `dict[str, dict[str, Expr]]` would let a reader
  skip the docstring (the codebase's `StandardFormEquation` uses the same plain style, so
  this is consistency, not a defect).
- `LinearizedSystem` is `frozen=True` but holds mutable dicts (noted in the adversarial
  review and deliberately left). A one-line docstring sentence "fields are not deep-frozen"
  would settle it.

### A8 (optional) A constant equation yields a silent all-zero dynamic row

- `vbl.py:273-274`: if δE is identically zero the row is all zeros. This is the documented
  and tested behaviour (`test_zero_equation_gives_zero_dynamic_row`) and consistent with the
  completeness rule (there is no *term* to complain about). A constant non-zero `E` (an
  inconsistent equation) also lands here. Consider a `warnings.warn` when `E` has no
  dependence on any configuration, velocity or input; not required.

### A9 (optional) `expand` on ill-typed `VVMul` operands

- `geomech/core/transformations/expand.py:81-95`: if a `VAdd` factor contains a
  `Transpose` addend (a row vector added to a column vector), distribution produces
  `VVMul(Transpose, Transpose)` and the constructor raises `ExpressionMismatchError` from
  inside `expand`. The input is already ill-typed, so this is garbage-in; mention only.

## B. What I verified independently

### B1 Kinematic rows (re-derived, then checked numerically against exact curves)

- S²: with q̇ = ω×q, δq = ξ×q and Jacobi, (δω − ξ̇ − ξ×ω) × q = 0; (I − qqᵀ)(ξ×ω) = 0
  because ξ×ω ∥ q; q·ξ̇ = −q̇·ξ = qᵀ Hat(ω) ξ from d/dt(q·ξ) = 0. Hence
  ξ̇ = q qᵀ Hat(ω) ξ + (I − q qᵀ) δω. Matches `vbl.py:246-248` and the printed blocks
  `q(q)'Hat(ω)`, `(I+q(q)'(-1))`.
- SO(3): δṘ = R η̂ Ω̂ + R δΩ̂ and d/dt(δR) = R Ω̂ η̂ + R η̂̇ give η̂̇ = η̂Ω̂ − Ω̂η̂ + δΩ̂ =
  (η×Ω)^ + δΩ̂, so η̇ = −Hat(Ω) η + δΩ. Matches `vbl.py:250-251`.
- Numeric: my own trajectory families q(t, ε) = exp(t Hat(a+εb)) exp(ε Hat(c)) q₀ (ω := q×q̇)
  and R(t, ε) = exp(t Hat(a+εb)) R₀ exp(ε Hat(c)) (Ω = Rᵀ(a+εb), η = vee(Rᵀ∂_εR)), with
  nested central differences h = 1e-5: |A[ξ]ξ + A[δω]δω − ξ̇| = **2.7e-7**,
  |A[η]η + A[δΩ]δΩ − η̇| = **2.8e-7** (terms are O(1); differencing noise).

### B2 Sign convention

E is the post-IBP coefficient (negative Euler–Lagrange); the pinned free rigid body
E = −JΩ̇ − Ω×JΩ + M gives δE = −J δΩ̇ + (Hat(JΩ) − Hat(Ω)J) δΩ + δM, so M = −coeff(δΩ̇) = J
(positive inertia), A[δΩ] = Hat(JΩ) − Hat(Ω) J, B = I. `to_standard_form` reports M = −J
for the same system (checked), confirming the sign flip is VBL's and deliberate.

### B3 Pinned systems against my own closed forms (own evaluator, 5 random states each)

- Rigid pendulum about the pivot (SO(3), Lee et al. Eq. 6.8): gravity moment
  −m g ρ × Rᵀe₃, δ(Rᵀe₃) = Hat(Rᵀe₃) η ⇒ A[η] = −m g Hat(ρ) Hat(Rᵀe₃); A[δΩ] =
  Hat(JΩ) − Hat(Ω) J; M = J; kinematic A[η] = −Hat(Ω). Max |library − hand| = **1.8e-15**.
- Spherical pendulum (S², Lee et al. Eq. 5.14): (ξ×q)×e₃ = Hat(e₃) Hat(q) ξ ⇒
  A[ξ] = −m g l Hat(e₃) Hat(q); M = m l² I; kinematic blocks as in B1. Max difference
  **0.0** (exact). Hanging equilibrium q = −e₃, ξ = (0.3, −0.4, 0): A[ξ]ξ =
  (−2.678, 3.571, 0) = −m g l ξ with m g l = 8.927. Restoring, as it must be.

### B4 Quadrotor with cable-suspended load, row δx_L (own derivation, own evaluator, own FD)

E_x = −(m_Q+m_L) ẍ_L + m_Q l q̈ − (m_Q+m_L) g e₃ + f R e₃ with q̈ = ω̇×q + ω×(ω×q).
Hand blocks: M[δẍ_L] = (m_Q+m_L) I, M[δω̇] = m_Q l Hat(q), A[η] = −f R Hat(e₃),
B[δf] = R e₃, A[ξ] = −m_Q l (Hat(ω̇) + Hat(ω)²) Hat(q), A[δω] = −m_Q l (Hat(ω×q) + Hat(ω) Hat(q)).
- Library M[δω̇], M[δẍ_L], A[η], B[δf] equal the hand matrices to machine precision.
- A[ξ], A[δω]: see A3 — differ as matrices (≤ 3.47) but agree on tangent directions
  (1.6e-15).
- Own central difference of E_x along q(ε) = exp(ε ξ̂) q, ω(ε) = exp(ε ξ̂)(ω + ε w),
  R(ε) = R exp(ε η̂), ε = 1e-6, versus my hand linearization: **1.6e-9** relative, 5 trials.
  The library's row passes the same check through `finite_difference_check` (observed
  1e-10 to 1e-12).

### B5 Engine change (`expand` VVMul rule)

- Shapes: the `VVMul` constructor admits exactly (plain, Transpose) → matrix and
  (Transpose, plain) → scalar; `_expand_vvmul` branches on `expr.type` and handles both.
  Each addend of `Transpose(VAdd)` becomes `Transpose(addend)`, so every produced `VVMul`
  is well-typed. `VAdd`/`MAdd` flatten nested sums (`_flatten_nodes`), so the result is
  flat for `_flatten_addends`.
- Invariants: children are expanded before distribution (bottom-up); the rule never calls
  `expand` on its own output (single pass); on a second pass no factor expands to a `VAdd`,
  so `expand(expand(x)) == expand(x)` (pinned by `TestIdempotent`). A `Transpose` factor
  whose inner is not a sum keeps the original node, so no other rule's input changes
  (pinned by `test_transpose_factor_without_a_sum_keeps_its_original_node`).
- Effect on existing output: I extracted the tree at `db86568` (`git archive`) and ran the
  same script against both trees, printing `compute_eom` and `to_standard_form` for the
  nine spec systems plus one extra Lagrangian with an explicit outer product over sums
  (½ m ẋ·(p pᵀ) ẋ, p = x − l q). **The nine systems are byte-identical before and after.**
  Only the outer-product Lagrangian changes, from `(x + q l(-1))((x + q l(-1)))' …` to the
  distributed sum; mathematically identical. Hence the CHANGELOG entry belongs under
  *Changed* with that scope stated (done, section D).

### B6 Completeness and robustness probes

- `Expr.has` recurses on composite nodes (`mixins.py:14`), so `_has_algebraic` raises for a
  perturbation two levels inside a calculus node (`d/dt(c δΩ)`), as the completeness rule
  needs. `_occurs_outside` sees through `Hat`, `Transpose`, `Vee`, `VVMul`, `SVMul`
  (probed directly).
- Undeclared non-constant symbol → `AlgebraicError` naming the term and suggesting the fix
  ("Is an input missing from `inputs`, or a symbol missing the 'Constant' attribute?").
- `Hat(ω) ω` (identically zero) linearizes to zero blocks; `Hat(q) ω (q·e₃)` linearizes
  (nonlinear in the state, linear in the perturbation).

### B7 Documentation claims

Ran the free-rigid-body snippet in `docs/vbl.md`: `print(lin)` output is identical to the
quoted block, line for line. `B[δf] = Re3`, `A[η] = RHat(e3)(-1)f` (quadrotor) and
`A[δω] = (ω m l (2))'` with `block_shape` = (1, 3) (elastic pendulum) are the actual
strings. `mkdocs build --strict` passes with the new nav entry.

### B8 Test hygiene

`git diff 293fab6 HEAD -- tests/test_vbl_adversarial.py`: seven `xfail(strict=True)`
decorators removed, the two "documents the present (wrong) behaviour" tests deleted
(`test_so3_tangent_vector_as_input_currently_moves_the_block`,
`test_mismatch_currently_yields_a_wrongly_shaped_block`), two comments updated. No
assertion changed. No `pytest.skip`/`skipif` in any new test file. Pinned tests assert
exact strings (`"Hat(\\rho)Hat((R)'e3)(-1)mg"`, `"q(q)'Hat(\\omega_{q})"`, …) or numeric
equality to a hand formula; docstrings carry the derivation. No mention of the Scala
implementations in any new Python or docs file (word-boundary grep over `geomech/`,
`tests/`, `docs/vbl.md`, `CHANGELOG.md`; the pre-existing citation in `README.md:7` and
`docs/manifold-implementation-notes.md` are outside this diff).

## C. Spec compliance (delta to the adversarial review's table)

All rows of the adversarial review's compliance table hold at HEAD; the one "Partially"
(exactly-one rule) is now "Yes" (`_partition_terms` + `_occurs_outside`, pinned in
`tests/test_vbl_validation.py` and the three F1 repros). `apply_manifold_rules` is not
called in the VBL path. No existing test was modified (`git diff db86568 HEAD --stat`
touches no pre-existing test file).

## D. Edits made in this review (working tree, uncommitted)

- `docs/vbl.md`: R/R³ bullet now indexes `M[d/dt δx]` and says the two symbols print
  identically; new subsection "Blocks on S² are defined on the constraint set" (A3, with
  the quadrotor-with-load example); new "Limitations" section (Transpose of sums / outer
  products and `Vee` opaque to extraction with a rewrite hint; stray perturbation symbols
  surface as the engine's `ValueError`; zero variation gives a zero row).
- `CHANGELOG.md`: the `expand` VVMul rule moved from *Added* to *Changed* with its exact
  scope (only Lagrangians containing an outer/inner product with a sum as a factor; none of
  the pinned systems affected). The separate "`linearize` validates more strictly …
  previously only at least one" bullet, which described a pre-release intermediate state,
  folded into the VBL entry as the behaviour a user sees.
- Suite, ruff and `mkdocs build --strict` re-run after the edits: 1322 passed, 5 xfailed;
  clean; builds.

## E. Final counts

`python -m pytest tests/ -q` → 1322 passed, 5 xfailed (1.3 s). `ruff format --check`:
65 files already formatted. `ruff check`: all checks passed. `mkdocs build --strict`: OK.
