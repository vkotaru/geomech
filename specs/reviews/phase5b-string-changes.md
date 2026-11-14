# Phase 5b: changed pinned strings

Every printed form that phase 5b changed, with the check that licensed it. Base:
`fix/phase5a-simplifier` at fc0d386. The procedure is the one of phase 5a
(`specs/reviews/scripts/pinned_snapshot.py`, `phase5a-string-changes.md`): a snapshot
records every `compute_eom`, `to_standard_form`, `linearize` and `integrate_by_parts`
call of the suite plus a chain for each builder system, and `compare` replays them with
the working tree, comparing `str` and the value of every output at five fixed states,
the linearization rows along admissible directions, and `check_eom`.

Two snapshots were taken:

| snapshot | tree | records | outputs |
|---|---|---|---|
| `base.json` | fc0d386 (worktree) | 562 | 4272 (+152 expected errors) |
| `step4.json` | phase 5b after steps 1–4 (structure, splitting, freeze, typing), before D2 | 563 (one more call, from `tests/test_phase5b_structure.py`) | 4273 (+152) |

## Gates per step

| step | compare against | strings changed | value difference | `check_eom` worse |
|---|---|---|---|---|
| 1. 5.5 constants, dedupe, public names, dead code, hoisted imports | fc0d386 | 0 | — | 0 / 216 |
| 2a. 5.6 split by node class; like terms grouped in one pass (`term_key`) | fc0d386 | 0 | — | 0 / 216 |
| 2b. inner products `aᵀb = bᵀa` in like terms | fc0d386 | 0 | — | 0 / 216 |
| 2c. literal spelling (`(2.0)` → `(2)`, `Scalar("0.5")` → `(0.5)`) | fc0d386 | 1 (below) | 0 | 0 / 216 |
| 2d. constant δ/d/dt are `ZeroVector`/`ZeroMatrix` | fc0d386 | 1 (the same one) | 0 | 0 / 216 |
| 3. D3 freeze | fc0d386 | 1 (the same one) | 0 | 0 / 216 |
| 4. 5.7 typing | fc0d386 | 1 (the same one) | 0 | 0 / 216 |
| 5. D2 separator | `step4.json` | 1594 | 0 (all 1594 identical values) | 0 / 217 |
| final | fc0d386 | 1593 | 0 | 0 / 216 |

The inner-product merge (2b) changed no recorded output: no pipeline output of the suite
or of a builder system holds both `aᵀb` and `bᵀa` in one sum. It is pinned by
`tests/test_phase5b_structure.py::test_inner_products_merge_either_way_round` with a
numeric check.

## Records that are not a pure separator insertion

The D2 gate asserts, for every changed output, `new.replace("*", "") == old` (no baseline
string contains a `*`: 0 of 4273). Against `step4.json`: **0 of 1594 records fail it**, and
every one of them has value difference 0 at all five states. Against fc0d386 one record
fails it, and it is the step-2c change, not D2:

| record | old (fc0d386) | new | value |
|---|---|---|---|
| `integrate_by_parts` output of `tests/test_trace.py::test_each_stage_reproduces_from_the_previous_one` (rigid pendulum, stage `ibp`) | `...J\Omega_{R})0.5+Dot(\Omega_{R},JCross(\Omega_{R},\eta_{R}))0.5+...` | `...J*\Omega_{R})*(0.5)+Dot(\Omega_{R},J*Cross(\Omega_{R},\eta_{R}))*(0.5)+...` | identical (0 relative difference) |

`Scalar("0.5", value=0.5)` (the `HALF` of the test systems) is now written `(0.5)`, the
spelling `Number(0.5)` gives it (phase 5a advisory 2). It is a raw IBP stage, not pinned by
any test.

## Changed outputs by kind (D2, against `step4.json`)

1594 outputs in 331 of the 563 recorded calls. By output label of the compare script
(blocks of the builder chains are prefixed `sf:` and `lin:`, those of direct calls carry
the block name only):

| label | outputs |
|---|---|
| `E` (`compute_eom` equations) | 282 |
| `M`, `A`, `B` (`linearize`, direct calls) and `M`, `f`, `G` (`to_standard_form`, direct calls) | `A` 382, `M` 232, `f` 106, `B` 8, `G` 6 |
| `lin:M`, `lin:A`, `lin:B` (builder chains) | 69, 288, 3 |
| `sf:M`, `sf:f`, `sf:G` (builder chains) | 117, 84, 3 |
| `ibp` (`integrate_by_parts`) | 14 |

The 2679 unchanged outputs hold no product (zero blocks, identities, sums of leaves).

## Pinned strings in the tests (D2)

Updated mechanically, then checked: a pytest plugin recorded every failing string
comparison whose two sides differ only by `*`; a token-level rewriter replaced the old
literal by the new one, only inside the failing test function, its class body or module
level, never in a call argument (a symbol name such as `Scalar("m1")`) and never as the
left operand of `in` unless the substring was matched in the new output; a second pass
did the same for f-strings and implicitly concatenated literals (placeholders aligned by
regex, separators inserted only in literal text). 13 rounds; the remaining cases were
edited by hand (listed below). Afterwards every changed test line was checked:
**318 changed lines in 25 files, every one equal to its old line with `*`s inserted.**

| file | hunks | lines |
|---|---|---|
| `tests/test_addition.py` | 3 | 3 |
| `tests/test_collect.py` | 1 | 1 |
| `tests/test_collect_fix.py` | 6 | 17 |
| `tests/test_collect_termination.py` | 3 | 3 |
| `tests/test_error_context.py` | 4 | 5 |
| `tests/test_extract.py` | 3 | 3 |
| `tests/test_ibp_prefactor.py` | 11 | 12 |
| `tests/test_manifold_pipelines.py` | 21 | 30 |
| `tests/test_multiplication.py` | 27 | 27 |
| `tests/test_operators.py` | 14 | 20 |
| `tests/test_phase1_adversarial.py` | 1 | 1 |
| `tests/test_phase2_adversarial.py` | 7 | 10 |
| `tests/test_phase3_adversarial.py` | 5 | 5 |
| `tests/test_phase3_fixes.py` | 3 | 3 |
| `tests/test_phase5a_adversarial.py` | 10 | 12 |
| `tests/test_phase5a_simplifier.py` | 12 | 19 |
| `tests/test_pull.py` | 15 | 15 |
| `tests/test_repr_debug.py` | 1 | 2 |
| `tests/test_row_forms.py` | 2 | 8 |
| `tests/test_rule_names.py` | 1 | 1 |
| `tests/test_s2_variation_fix.py` | 11 | 27 |
| `tests/test_scalar_adversarial.py` | 42 | 58 |
| `tests/test_validation.py` | 2 | 3 |
| `tests/test_vbl.py` | 21 | 29 |
| `tests/test_vbl_adversarial.py` | 3 | 4 |
| **total** | **229** | **318** |

After the check one line of `tests/test_s2_variation_fix.py` (the `s2_gyroscopic`
expected string, 110 characters with separators) was split into two concatenated
literals for the 99-character limit; the value is unchanged.

Edited by hand (each still a separator insertion): `tests/test_operators.py`
(`"x" + text` → `"x*" + text`, `text + "x"` → `text + "*x"`), `tests/test_phase2_adversarial.py`
(`endswith("+f*(2))")`), `tests/test_phase5a_adversarial.py` (`"m*g*l" in s`),
`tests/test_scalar_adversarial.py` (the f-string `f"(m*x*{dx}*{dx}*(-1)+m*x*x*{DDX}*(-1)+f)"`,
`"m*x*x*(-1)"`, the parametrize id `"Hat(a)*w"`), `tests/test_validation.py`
(`str(om) + "*(-1)"`).

## Test changes that are not pinned strings

| file | change | why |
|---|---|---|
| `tests/test_vectors.py::test_delta_constant`, `tests/test_matrices.py::test_delta_constant` | `"0"` → `"0v"`, `"O"` → `"0"` | step 2d: the variation of a constant vector/matrix is `ZeroVector`/`ZeroMatrix` (phase 4 follow-up 5); a zero either way |
| `tests/test_addition.py`, `tests/test_calculus.py`, `tests/test_geometry.py` (`nodes == [a]`) | `== (a,)` | D3: operands are a tuple |
| `tests/test_invariants.py`, `tests/test_phase4_adversarial.py` | `node.nodes[i] = x`, `node.nodes = [...]`, `bad.flags.is_numeric = False` → `object.__setattr__` (helper `_replace_child`) | D3: these tests build malformed trees for `check_tree`; that now needs the explicit bypass. Without the change one parametrized case of `test_check_tree_names_the_invariant_and_path` would have been *skipped* (its builder raised), not failed |
| `tests/test_utils.py` | tests of `is_leaf`/`has_nested_add` replaced by tests of `addends`, `sum_of`, `flatten_mul` | 5.5: the two functions are deleted |

## Docs, notebooks, README

| file | change |
|---|---|
| `docs/pipeline.md` | `summary(chars=70)` block regenerated (4 lines: `·0.5` → `·(0.5)`, step 2c; `repr` has no separator change); like-term examples with separators and `aᵀb s − bᵀa s → 0`. The rule list (42 lines, the same counts) reproduces unchanged |
| `docs/conventions.md` | products print `*`; zero leaves; numbers (literal spelling, bare `Number(2.0)`); like terms (inner products commute, `term_key`); the `str` identity section (calculus operand class, `has()`, separators, frozen nodes) |
| `docs/troubleshooting.md` | 7 quoted messages regenerated from the code (each a separator insertion; checked by `tests/test_phase4_adversarial.py` and `tests/test_troubleshooting_doc.py`); one History row for the D2 collision |
| `docs/vbl.md` | the `print(lin)` block (2 lines) and three inline blocks |
| `README.md`, `docs/getting-started.md` | the `str` note (products print `*`); README's `_eliminate` sketch shows the per-class handlers |
| `notebooks/quadrotor.ipynb`, `notebooks/quadrotor_payload.ipynb` | re-executed; one output line each, `+ Re3f` → `+ R e3 f` (`display_standard_form` now writes the separator of an input key as a space instead of printing `str`); `execution` cell metadata stripped, `language_info.version` kept as it was. The other two notebooks re-executed byte-identical |
