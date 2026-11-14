# Fix plan for the 2026-11-13 audits

Sources: `specs/reviews/engine-fuzz-report.md` (73 535 fuzz cases),
`specs/reviews/engine-audit.md` (read-through with reproducers),
`specs/reviews/dynamics-audit.md` (independent Euler–Lagrange harness),
`specs/reviews/debuggability-audit.md` (future-debugger walkthrough).
Baseline: `feature/vbl` at 54d97b9, 1322 passed / 5 xfailed.

Findings are referenced as `fuzz-N`, `eng-FN`, `dyn-FN`, `dbg-FN`. Three independent
methods (fuzzing, reading, EL harness) converged on the same two silent-wrong-equation
bugs, which is strong evidence that the list below is the complete set of wrong-value
defects reachable today.

Each phase is one branch stacked on the previous one, one PR each, suite + ruff + mkdocs
strict green at every commit. Existing passing tests that must change are listed per
phase; per project rules they are changed only with the author's consent, which is given
by approving this plan.

---

## Phase 0 — Record the audits (no code)

Commit the four reports and `specs/reviews/scripts/engine_fuzz.py`. Move the dynamics
audit's Euler–Lagrange harness (its Appendix A) into
`specs/reviews/scripts/euler_lagrange_check.py` so the report's "scratch directory"
references point into the repo.

## Phase 1 — Wrong equations (mathematics)

Everything here changes derived output. Each item ships with a pinned regression test
whose expected equation is derived from Lee et al. or by hand in the docstring, plus a
numeric check against the independent Euler–Lagrange harness.

### 1.1 S² variation rule `δω = ξ̇ − 2 ω×ξ`   (dyn-F1, fuzz-2)

- `geomech/core/transformations/manifold_rules.py:26,66-70`: coefficient −1 → −2, with
  the derivation in the module docstring: Lee et al. eq. (5.13),
  `δω = −S(ω)ξ + (I − qqᵀ)ξ̇`, and `qqᵀξ̇ = ω×ξ` on the constraint set `ξ ⊥ q`.
- `geomech/core/base/expressions.py:343-351` `TS2.delta(substitute=True)`: delete (never
  called; second source of truth) or route through `manifold_rules`.
- `docs/manifold-implementation-notes.md`: replace "needs verification" with the
  derivation and the symptom (gyroscopic terms come out at half strength).
- Paper: the accompanying paper and `03_toolbox.tex:116` state the wrong
  formula; fix both.
- Existing test to change: `tests/test_manifold_rules.py::TestS2ManifoldRules::test_variation_of_omega_expands`
  (pins the −1).
- New pinned pipeline tests, each with a numeric EL check: spherical pendulum on a
  rotating base (Lee et al. eq. 5.43, `∂L/∂ω` not parallel to `ω`), gyroscopic `c·ω`
  term (hand derivation: `E = −m l² ω̇ + 2 ω×c − m g l q×e₃`), anisotropic `½ ω·Aω`,
  S²×SO(3) coupling `k ω·Ω`. Also a test that the ω-form and q̇-form of the same
  Lagrangian give numerically identical equations.

### 1.2 `collect` false identity and dead end   (eng-F1, eng-F6, dyn-F2, dyn-F3, fuzz-1)

- `geomech/core/math/collect.py:64-77` and `:99-117`: when the target is anywhere inside
  `MVMul.right`, peel the matrix with its transpose and recurse,
  `Dot(M·f(v), r) → collect(Dot(f(v), Mᵀ r), v)`; delete the `(M·b)×c` branch.
  Terminates because each peel removes one matrix layer.
- Existing test to change: `tests/test_collect.py:204-213 test_mvmul_cross_vec_in_mvmul`
  (pins the false identity); new expectation `Dot(v, Mᵀ(b×c))` plus numeric equality.
- New pinned pipeline tests with EL check: SO(3) body with `(JΩ)·((Rᵀe₃)×b)`; inertial
  angular momentum `(R(JΩ))·e₃`; S² anisotropic in q̇-form `½ q̇·A q̇`; SO(3)×SO(3)
  relative-rate coupling `Ω₁·(R₁ᵀR₂Ω₂)`; `½(Ω×ρ)·J(Ω×ρ)` (currently crashes).

### 1.3 `Hat(Vee(X)) → X` only for skew-symmetric X   (eng-F2, fuzz-3)

- `simplify.py:161-167`: fire only if `X` is a `Hat`, a `SkewSymmMatrix`, or
  `is_skew_symmetric`. Add `is_skew_symmetric` to `Hat` (always), and to
  `SMMul`/`MAdd`/`Transpose` of skew operands.
- Existing test to change: `tests/test_simplify.py:376-380 test_hat_vee_cancel` (uses a
  plain `Matrix`); use `SkewSymmMatrix` and add a test that `Hat(Vee(I))` is not `I`.

### 1.4 `extract_linear_coeff` holes   (eng-F3, eng-F4, dbg-F7, fuzz-7, fuzz-8, fuzz-9)

- `extract.py:650-651` and `:689-693`: replace the silent `ZeroMatrix` with
  `_zero_or_raise`; explicit `Hat`/`Transpose`-of-query cases raise `AlgebraicError`
  ("linear but has no left matrix coefficient"). Correct the dispatch table in
  `paper_context.md` and the docstring claim at `extract.py:84`.
- `extract.py:426-450`: move the both-sides guard to the top of the `Cross` branch and
  add it to `MVMul` → `NonLinearError`.
- `extract.py:148`: `NotImplementedError` → `AlgebraicError` with the expression.
- Tests: the fuzz harness's quadratic trees as unit tests; `Hat(x)` query raises.

### 1.5 Zero outer product has the wrong type   (fuzz-4)

- `simplify.py:141-144`: return `ZeroMatrix` for a matrix-typed `VVMul` with a zero
  factor, `Zero` only for the inner product.

### 1.6 Numeric oracle in the package   (dbg-F13; roadmap feature 1)

- New `geomech/numeric.py`: `Evaluator`, `hat`, `vee`, `rand_rot`, `expm_so3`,
  `sample_state` moved from `tests/numeric_helpers.py` (which becomes a re-export), plus
  `evaluate(expr, env)` and `check_eom(L, dW, variables, eom, trials=5)` built from the
  dynamics audit's harness: finite-difference Euler–Lagrange residual on R, R³, S², SO(3)
  and products, returning the worst relative error.
- Every pinned pipeline test in `tests/test_manifold_pipelines.py` and
  `tests/vbl_systems.py` gains a `check_eom` assertion. This is the gate for 1.1–1.5 and
  for every future engine change.
- `docs/`: one page "Verifying an equation numerically".

Done when: suite green, `engine_fuzz.py` reports 0 wrong-value cases, the dynamics
audit's six wrong rows and four crash rows all pass `check_eom` at ≤1e-5, paper formula
fixed.

## Phase 2 — Error instead of silence (validation)

No change to any correct output; every item turns a silent wrong or ill-posed result into
an exception that names the symbol. Message style: the existing `linearize` messages.

| # | Item | Findings |
|---|---|---|
| 2.1 | Unknown `attr` string → `ValueError` listing valid attributes; accept `"TangentVector"` | eng-F14, dbg-F1 |
| 2.2 | After extraction, any `Variation` symbol left in δS that is not a listed variable → `ValueError` ("add the variable or mark it Constant") | dbg-F2, dyn-F5, eng-F17 |
| 2.3 | `SystemVariables`: reject duplicate names (across groups and across Scalar/Vector with the same name), tangent/variation vectors, non-SO3 matrices; `linearize._check_eom` compares list lengths | dyn-F4, dyn nits, dbg-F8 |
| 2.4 | `S2`/`SO3`/`Matrix` constructors copy the `attr` list | eng-F15 |
| 2.5 | S² tangent rules compare the parent by name, not `is` | dbg-C6, eng-F17 |
| 2.6 | `compute_eom(L, 0/None, …)` accepted as `Zero`; vector-typed `dW` → `TypeError` | dyn-F6 |
| 2.7 | `manifold_rules._rebuild` unknown node → `TypeError`, not silent passthrough | dbg-F7 |
| 2.8 | `TimeDerivative`/`Variation` of an `S2`/`SO3` leaf → `ValueError`; `Transpose()` without child → `TypeError` | dbg-F8 |
| 2.9 | `collect`: target on both sides → `NonLinearError` instead of `RecursionError`; non-scalar argument → `TypeError` | fuzz-5, dbg-A18 |
| 2.10 | `full_simplify` emits a `warning` when the iteration cap is hit | dbg-F12 |
| 2.11 | `Number(0)` flagged `Zero`, `Number(1)` flagged `Ones`; `_simplify_mul` returns `Zero` on a zero product | eng-F11 |
| 2.12 | `examples/rigid_pendulum.py:32` `J2 = Matrix("J1")` → `"J2"` | dbg-A27 |

Tests: one test per item asserting the exception type and that the message names the
symbol (style of `tests/test_vbl_validation.py`).

## Phase 3 — Crashes on legal input, operators, printing

| # | Item | Findings |
|---|---|---|
| 3.1 | `Transpose.__mul__` branches on `self.type`: `Rᵀ v`, `Rᵀ M`, `Rᵀ s` work; `Transpose(v) * s` returns `SVMul(Transpose(v), s)` not a double transpose | eng-F5, eng-F17, dyn-F7 |
| 3.2 | `Transpose(Transpose(x))` collapses to `x` at construction, so the row-vector encoding survives `expand`/`simplify` | fuzz-6 |
| 3.3 | Numeric operands on both sides: `2*x`, `x+1`, `1+x`, `x-1`, `-x` via `__radd__`, `__rmul__`, `__rsub__`, `__neg__` and wrapping in `__add__`/`__sub__` | eng-F7, dbg-A29–A32 |
| 3.4 | Constructor errors name the operands (`Dot cannot combine SCALAR 'm' and VECTOR 'x'`); unary ops say "expected VECTOR, got SCALAR 'm'"; remove the bare `raise ExpressionMismatchError` and the three empty `NotImplementedError`s; delete unreachable `UndefinedCaseError` | dbg-F5, dbg-A6, dbg-A16 |
| 3.5 | `to_latex` is precedence-aware: parenthesise `Cross`/`VAdd`/`MVMul` operands of `Cross`, `Transpose`, `MVMul`; print `(b·c)\,a`; coefficient before symbols. `print_eom`, `eom_to_latex`, `render_eom` use `to_latex`, not `str` | dyn-F8, dbg-F6 |
| 3.6 | `repr` prefixes non-scalar leaves (`v:x`, `M:J`, `S2:q`) and marks constants; `print_tree` defaults to `style="indent"` and shows orthogonal/identity/skew flags and the TS2/TSO3 parent | dbg-F6 |

Tests: round-trip tests for 3.5 (every ambiguous pair in the dynamics audit's table
prints differently) and operator tests for 3.1–3.3. The paper's listings are regenerated
from `to_latex` after 3.5 and diffed by eye.

## Phase 4 — Debuggability

| # | Item | Findings |
|---|---|---|
| 4.1 | `compute_eom(..., trace=EOMTrace())` keeps δL, δS, after manifold rules, after simplify, after IBP, after expand, raw and final per-variable equations; same for `linearize` (per-row δE, buckets, raw blocks). Export `apply_manifold_rules` | dbg-F3 |
| 4.2 | Engine errors re-raised with pipeline context: stage, variable, term ("a d/dt(η) survived integration by parts in term …"); `_check_coefficient_free` keeps the inner message | dbg-F4, dbg-C3 |
| 4.3 | `logging.getLogger("geomech")`, one debug line per stage; `GEOMECH_DEBUG=1` enables it plus the invariant checker | dbg-F11 |
| 4.4 | `geomech/utils/invariants.py: check_tree(expr)` (operand types, no nested same-class Add, VVMul has one Transpose child, no calculus node over a manifold leaf, tangent parents set); run in debug mode at pipeline entry and unconditionally in tests | dbg-F8 |
| 4.5 | Named rewrite rules through one `rule(name, before, after)` helper in `vector_rules`, `simplify`, `pull`, `collect`, `expand`; `GEOMECH_TRACE=1` logs every firing; names listed in each module docstring next to their hypotheses | dbg-F9, dbg-C1 |
| 4.6 | Rename transformation modules so the package no longer shadows them (`_expand.py` …, or export under different names) | dbg-F14 |
| 4.7 | Docs: "Pipeline stage by stage" (one order, matching the code, replacing the three inconsistent ones in `eom.py`, `README.md`, `getting-started.md`); "Troubleshooting" (the silent cases of 2.x with their old symptoms, the two AlgebraicError causes, `print_tree`); "Conventions" (E = −Euler–Lagrange, M = −inertia, VBL sign convention, zero constants, `str` identity). Fix stale text: `variables.py:12`, `extract.py:11-13`, `manifold-implementation-notes.md` | dbg-F10 |

## Phase 5 — Simplifier completeness and structure

Lower priority: nothing here produces a wrong equation. Worth doing before the paper's
examples are generated, because unreduced terms make listings longer.

| # | Item | Findings |
|---|---|---|
| 5.1 | `vector_rules` and `pull` recurse through every node type (`Dot`, `Hat`, `Vee`, `Transpose`, `MVMul`, `MMMul`, `SMMul`, `MAdd`, calculus) using one `_rebuild` helper | eng-F8, dbg-C1, fuzz-11 |
| 5.2 | `expand` distributes `Vee` over `MAdd` | fuzz-10 |
| 5.3 | `IdentityMatrix` flagged symmetric; `Transpose(skew) → −skew`; `MMMul` reassociation mirror of the `MVMul` rule | eng-F9 |
| 5.4 | Like-term merging flattens `Mul` chains before splitting the numeric coefficient | eng-F10 |
| 5.5 | One `core/base/constants.py` (`NEG_ONE`, `HALF`, `negate`, `zero_like`); dedupe `_flatten_addends`, `_flatten_mul`, `_zero_for`, `_combine`; promote `has_algebraic`, `same_node`, `check_variables`, `wrap_numeric` to public names; delete dead `"TangentVector"` parsing gap, `has_nested_add`/`is_leaf`, no-op lines; one vocabulary (variation vector / perturbation) | dbg-F15, vbl-final A3–A9 |
| 5.6 | Split `_eliminate`, `linearize`, `to_latex`, `_collect_dot`, `pull` by node class (also what 4.5 needs) | dbg-F15 |
| 5.7 | Type annotations on mixins/bases; `mypy geomech` in CI (9 errors today, 102 pyright) | dbg-F16 |

## Decisions (answered by the author on 2026-11-13)

All phases 0–5 approved. Answers, as given: D1 (a) named scalars with a value stay
symbolic; D2 (a) separator in the product family plus a type guard in `has()`; D3 and D4
"compare by name; drop immutability claim" (D3 (b) now, freezing in phase 5; D4 (a)); D5
keyword argument `trace=` (the recommendation, implemented in phase 4).

The options as originally presented:

| # | Question | Options | Recommendation |
|---|---|---|---|
| D1 | `Scalar("m", value=2)` is folded as the literal 2 and loses its name (eng-F12) | (a) fold only anonymous numbers, keep named symbols symbolic; (b) document that `value` makes a literal | (a) |
| D2 | Equality is `type + str`; `Mul(ab,c) == Mul(a,bc)`, `Variation(Scalar x) == Variation(Vector x)` (eng-F13) | (a) print the `Mul` family with a separator and give `has()` the type guard; (b) full structural equality | (a) now; (b) is a larger refactor and not needed for correctness |
| D3 | Tree is not actually immutable (eng-F16) | (a) `frozen=True` + tuple children; (b) drop the claim from the contributor notes | (b) now, (a) when 5.6 touches the classes anyway |
| D4 | Two `S2("q")` objects compare equal but behave differently in the tangent rules (dbg-C6) | (a) compare by name everywhere (2.5); (b) intern manifold leaves by name | (a) |
| D5 | Trace API shape (4.1) | `compute_eom(..., trace=obj)` vs `compute_eom_trace(...) -> EOMTrace` | keyword argument; one function to remember |

## Order and size

| Phase | Content | Size |
|---|---|---|
| 0 | commit audits | minutes |
| 1 | wrong equations + numeric oracle | largest; ~15 new pinned systems with EL checks; 3 existing tests change |
| 2 | validation | medium; 12 small changes, 12 tests |
| 3 | crashes, operators, printing | medium; `to_latex` rewrite is the big item |
| 4 | debuggability | medium; trace object, logging, invariants, docs |
| 5 | completeness, structure, typing | large but mechanical; can be split |

Phase 1 is the only one that changes equations and the only one that touches the paper.
Phases 2–4 are independent of each other and could run as parallel loops after phase 1.
Each phase uses the same loop as the previous features: spec (this document's section) →
implement → adversarial review with tests → fix → final review.
