# Phase 5a implementation report: simplifier completeness

Branch `fix/phase5a-simplifier`, from `fix/phase4-debuggability` at 47cce2d. Not committed.
Scope: fix plan items 5.1–5.4, decision D1 (option a), and the simplifier follow-ups that
the phase 1–4 final reviews assigned to phase 5. Items 5.5–5.7 (constants module,
splitting by node class, typing) are not part of this phase.

## Gate results

| check | result |
|---|---|
| suite, `-W error` | 2171 passed, 2 xfailed (base: 2120 passed, 5 xfailed) |
| suite, `GEOMECH_DEBUG=1 -W error` | 2171 passed, 2 xfailed |
| `ruff format --check`, `ruff check` (geomech, tests, specs/reviews/scripts) | clean |
| `mkdocs build --strict` | builds |
| `tests/test_rule_names.py`, `tests/test_rule_coverage.py` | pass (new rules are in the module docstrings) |
| `engine_fuzz.py --cases 3000` | 0 wrong-value findings; 249 distinct findings vs 202 at base, the difference is all `flag-dependent` in the adversarial campaign (142 vs 95: rules that fire on declared flags now reach more nodes); every other category is identical |
| examples (`examples/*.py`) | all run |
| notebooks | re-executed in place, `execution` cell metadata stripped, kernel `language_info.version` kept as before; 3 of 4 change output text only (numbers last, `(-1.0)` → `(-1)`, parenthesised negation in `repr`) |
| numeric equivalence (`pinned_snapshot.py compare` against 47cce2d) | 655 changed outputs; 649 identical in value to ≤ 7.3e-16 relative; 6 S² `linearize` blocks differ as matrices but every linearization row is identical along admissible directions; 0 rows differ; `check_eom` unchanged on all 197 calls (largest ratio 1.000000006); details in `phase5a-string-changes.md` |

## Timing

`compute_eom` on the quadrotor with load (`tests/vbl_systems.quadrotor_with_load`), median
of 20 runs after one warm-up, three repetitions each:

| tree | median |
|---|---|
| 47cce2d | 10.8–11.0 ms |
| phase 5a | 12.0 ms |

A slowdown of 1.1×, well inside the 2× budget. The full suite runs in the same 9 s.
Equation lengths of that system: `δx_L` 180 → 180, `ξ` 285 → 285 (same terms, numbers
last), `η` 69 → 67.

## Items

### 5.1 `vector_rules` and `pull` recurse through every node

- `geomech/core/transformations/utils.py:32` `rebuild(expr, children, caller)` and `:57`
  `map_children(expr, fn, caller)`: the one structural helper. N-ary sums, the six binary
  products, `Dot`, `Cross`, and the unary `Hat`, `Vee`, `Transpose`, `Variation`,
  `TimeDerivative`, `TimeIntegral`; any other class raises `TypeError` naming the caller.
  `map_children` returns the same object when no child changed.
- `geomech/core/transformations/manifold_rules.py:119` now uses it (its private
  `_rebuild` is gone; the `TypeError` for an unknown node is kept, message
  "apply_manifold_rules: cannot rebuild a … node").
- `geomech/core/transformations/vector_identities.py:152` `vector_rules` = `map_children`
  + `_at_node` (`:158`): bottom up over every node class. Every rewrite still goes
  through `rule()`; the rule set is unchanged. Docstring updated (the "recurses only
  through" sentence is replaced).
- `geomech/core/transformations/pulling.py:94` `pull` = `map_children` + `_at_node`
  (`:100`). New rules, all through `rule()`:
  `pull.row_mvmul_scaled` (`(rᵀ s)(A t) = (rᵀ A)(s t)`, `:84`), `pull.vvmul_scaled`
  (scaled column of an inner product, scaled factors of an outer product, `:65`),
  `pull.transpose_scaled` (`(v s)ᵀ = vᵀ s`, `(M s)ᵀ = Mᵀ s`, `:190`), `pull.hat_scaled`
  (`:196`), `pull.vee_scaled`. `pull.vvmul_scaled_row` now also hoists a scaled column on
  the right in the same step. Calculus nodes are recursed into, with no rule (the
  operand of `d/dt` is not constant).
- Effect: every reproducer of engine audit F8 and fuzz finding 11 reduces (pinned in
  `tests/test_phase5a_simplifier.py::test_rules_fire_below_every_node`); `Q(Ω×Ω)` in the
  SO(3)×SO(3) coupling is gone (phase 1 follow-up 5.1).

### 5.2 `expand` distributes `Vee` over `MAdd`

- `geomech/core/transformations/expansion.py:223` rule `expand.vee_over_madd`; the "Vee
  is not distributed" paragraph of the docstring is replaced by the rule line.

### 5.3 identity, skew transpose, `MMMul` associativity

- `geomech/core/base/expressions.py:627` `IdentityMatrix` carries `SymmetricMatrix`, so
  `Iᵀ → I` by `simplify.transpose_symmetric`; `_is_identity` (`simplification.py:412`)
  also accepts `Transpose(I)`.
- `simplification.py:343` `simplify.transpose_skew`: `Sᵀ = −S` for a skew-symmetric
  operand that is not a `Hat` (kept as `simplify.transpose_hat_skew`) or a scaled
  matrix (handled by `simplify.transpose_scaled`).
- `simplification.py:273` `simplify.mmmul_reassociate` (`(A·B)·C → A·(B·C)`) and `:282`
  `simplify.mmmul_combine` (`A·(B·C) → (A·B)·C` when `A·B` simplifies to a non-product),
  the mirror of the `MVMul` pair: `A Rᵀ R → A`, `Rᵀ (R A) → A`.

### 5.4 like terms and numbers

- `simplification.py:524` `_split_term`: an addend becomes (product of all numeric
  literals of its flattened `Mul` chain, including the scalars of nested `SVMul`/`SMMul`;
  the other scalar factors; the vector/matrix base). `:597` `_same_factors`: multiset
  equality with `_commutative_equal` per factor (which now also compares `Mul` chains as
  flattened multisets). `:635` `_collect_like_terms` merges equal (base, factors) pairs
  and keeps an unmerged term exactly as it was (no churn for terms that did not merge).
  Vector and matrix terms with symbolic coefficients cancel: `a k − a k → 0v`
  (phase 3 follow-up F4).
- `(vᵀ s) w` and `s (vᵀ w)` merge: `pull` brings the first to `(vᵀ w) s` and the
  multiset comparison equates the two (phase 3 follow-up F3/R1).
- `simplification.py:683` `_simplify_mul`: the numeric literals of a product chain are
  folded into one number placed last (`(-1)mll → mll(-1)`); a chain already in that form
  is returned as built. `:434` `_exact`: a float with an integer value becomes an `int`
  (`(-1.0) → (-1)`, `(2.0) → (2)`), used by `mul_fold`, `add_fold` and the like-term
  rebuild. Non-integral floats are kept (`x(0.30000000000000004)`).
- The three `k Ω₁ × QΩ₂` terms of the SO(3)×SO(3) equation now merge (phase 1 follow-up
  5.4): `E₁` drops from 8 to 5 terms, `E₂` from 7 to 5 (counting the torque).
- No new rule names: these are the existing `simplify.add_fold`/`vadd_fold`/`madd_fold`/
  `mul_fold` rules, whose docstring entries now state the like-term definition.

### D1: named scalars with a value stay symbolic

- `geomech/core/base/expressions.py:669` `is_number_literal(expr)`: a scalar leaf whose
  printed name parses to its value (`Number(2)` → `(2)`, `Zero`/`One`,
  `Scalar("0.5", value=0.5)`, `Scalar("(-1)", value=-1)`).
- `simplification.py:407` `_is_one` and `:425` `_is_numeric_leaf` use it, so
  `Scalar("m", value=2)` is never folded, dropped (value 1) or summed with another value;
  `printing.py:505` (`to_latex`) and `printing.py:172` (`repr` minus sign) as well.
- `geomech.numeric.Evaluator`/`evaluate` are unchanged and still use `value=` for named
  scalars, so `check_eom` of `stiff_point_mass(k_value=…)` and the phase-1 D1 test pass.
- Changed test: `tests/test_phase1_adversarial.py::TestOracle::test_scalar_with_value_is_folded_by_engine_and_used_by_evaluator`
  pinned the old folding (`"m" not in …`); it now asserts `"m" in …` (docstring says
  why). Engine audit F12's four reproducers are pinned in
  `tests/test_phase5a_simplifier.py::test_named_scalar_with_value_is_a_symbol`.

### Follow-ups assigned to phase 5

| follow-up | where | done |
|---|---|---|
| unmerged `kΩ₁×QΩ₂` (phase 1) | 5.4 above | merged |
| `Q(Ω×Ω)` unreduced (phase 1) | 5.1 above | reduced |
| rule log lines whose two sides print the same (phase 4 A1) | `geomech/utils/trace.py:121` (`str` shown when the `repr`s agree), `:132` `_cut_pair` (a side over 300 characters is cut around the first difference) | never logged; `tests/test_phase5a_simplifier.py::test_no_logged_rule_line_has_two_equal_sides` |
| nested negation in `repr` (phase 4 A2) | `geomech/utils/printing.py:126` `_paren_negative` | `-(-hat(v:ρ*))`, `hat(v:ρ*)·(-(-1))` |
| `collect.flip` with 6 call sites (phase 4 item 9) | `geomech/core/math/collecting.py:167, 179, 227, 238, 242, 247` | `collect.flip_target_right`, `flip_cross_mvmul`, `flip_mvmul_pair`, `flip_mvmul_left`, `flip_mvmul_right`, `flip_cross_left` (docstring lists each with its call-site hypothesis) |
| `ibp.dot_prefactor` hides `_t_diff_product` (phase 4 item 9) | `geomech/core/math/ibp.py` (`_apply_ibp`, the `a_rate` line) | `ibp.product_rule` logs `d/dt(a) ==> Σ …` before `ibp.dot_prefactor` |
| `Vee` over `MAdd` (phase 4 item 2) | 5.2 | done |
| vector rules do not descend (phase 4 item 1, troubleshooting row) | 5.1; `docs/troubleshooting.md:17` rewritten (only `q × (ω̇ × q)` is still kept, by design) | done |

## Rule names added

`pull.row_mvmul_scaled`, `pull.vvmul_scaled`, `pull.transpose_scaled`,
`pull.hat_scaled`, `pull.vee_scaled`, `expand.vee_over_madd`,
`simplify.transpose_skew`, `simplify.mmmul_reassociate`, `simplify.mmmul_combine`,
`ibp.product_rule`, and `collect.flip` split into `collect.flip_target_right`,
`collect.flip_cross_mvmul`, `collect.flip_mvmul_pair`, `collect.flip_mvmul_left`,
`collect.flip_mvmul_right`, `collect.flip_cross_left`. All are in their module
docstrings with the hypothesis; `test_rule_names.py` and `test_rule_coverage.py` pass.

## Numeric equivalence (the procedure required for changing pinned strings)

1. `specs/reviews/scripts/pinned_snapshot.py` was written first and run on a worktree of
   47cce2d (`snapshot`). It records every `compute_eom`, `to_standard_form`, `linearize`
   and `integrate_by_parts` call of the whole suite (not only the listed files) plus a
   chain for each builder system of `tests/vbl_systems.py`, `tests/phase1_systems.py`,
   `tests/test_eom_oracle.py`: 509 records. A self-compare of the snapshot reports 0
   changes (the fixed states are deterministic). Re-create the baseline with
   `git worktree add <dir> 47cce2d`, copy the script there, and run `snapshot`; the
   baseline file (2.9 MB) is not in the repository.
2. `compare` with the phase-5a tree: see the gate table and
   `phase5a-string-changes.md`. One deviation from "every linearize block numerically
   identical": six blocks of three S² builder systems differ as matrices. They differ by
   terms that vanish along admissible perturbations (on S² the blocks are not unique; the
   old engine produced them from an equation written differently), and the comparison of
   the linearization rows along admissible directions, which is what
   `finite_difference_check` tests, is identical to 1e-15. No test pins these blocks. I
   judged this not to be a bug; if the author wants raw-block identity on S², the only
   way is to stop applying tangency identities inside nested nodes, which undoes 5.1 on
   S².
3. Pinned strings were updated mechanically from the replayed outputs (a test passes only
   with the replayed string, and the replayed string's value was checked in step 2).
   The full table is in `phase5a-string-changes.md` (93 hunks in 11 test files).

## Counts

| | count |
|---|---|
| pinned-string hunks changed in existing tests | 93 (11 files) |
| strict xfails that now pass (markers removed, docstring says "Fixed in phase 5a") | 3, all in `tests/test_scalar_adversarial.py`: `test_numeric_coefficients_have_one_type`, `test_f_combines_like_terms`, `test_xi_f_exact_structure` |
| structural asserts adjusted to numbers-last | 3: `tests/test_ibp.py` (2), `tests/test_simplify.py::test_chain_mul_numerics` |
| semantic test changes | 3: the D1 test above; `test_hat_vee_simplification_preserves_value[Transpose]` compares with `full_simplify(X)` since `Kᵀ` now prints `K(-1)`; `TestErrorPathsHeld::test_unsupported_structures_raise_not_drop` split, two cases moved to a new test asserting the linearization is correct (`finite_difference_check`) |
| new tests | 48 in `tests/test_phase5a_simplifier.py` (every F8/F9/F10/F12 reproducer, fuzz 10/11, the follow-ups, each with a numeric check) |
| docs | `pipeline.md` (stage 3, 4, 5 text, the regenerated `summary(chars=70)` block, the rule-firing list: 42 lines), `conventions.md` (numbers, like terms, D1), `troubleshooting.md` (row 17 and two quoted messages that changed form), `vbl.md` (one quoted block) |
| CHANGELOG | Added (phase 5a simplifier items, `is_number_literal`, rule names, snapshot script) and Changed (printed forms, D1, S² linearize blocks, rule-trace lines) |

## Notes for the reviewer

- The numbers-last order is a global change of printed form (`(-1)mll` → `mll(-1)`). It
  is what makes like-term merging and the pull changes print deterministically; without
  it, the same equation printed `(-1)mgl` or `mgl(-1)` depending on the hoisting path.
- Symbolic coefficients are not factored: `x a + x b` stays, so the quadrotor-with-load
  translational equation keeps `ẍ_L m_Q(-1) + ẍ_L m_L(-1)`.
- `q × (ω̇ × q)` is still kept (tangency of `ω̇` is not structurally provable; unchanged
  design).
- `vector.transpose_transpose` is now unreachable (the constructor collapses a double
  transpose) but kept; deleting dead code is item 5.5.
