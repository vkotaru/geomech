# Phase 5b implementation report: structure, freeze, typing, product separator

Branch `fix/phase5b-structure`, from `fix/phase5a-simplifier` at fc0d386. Committed as f852c11 (implementation) and
d102af4 (review fixes F1–F5).
Scope: fix plan items 5.5–5.7, decisions D2 (a) and D3 (a), the phase-5 follow-ups of
`phase4-final.md` and `phase5a-final.md`, and the phase-5a advisory items. The string
changes are in `phase5b-string-changes.md`.

## Gate results

| check | result |
|---|---|
| suite, `-W error` | 2307 passed, 2 xfailed (base: 2238 passed, 2 xfailed) |
| suite, `GEOMECH_DEBUG=1 -W error` | 2307 passed, 2 xfailed |
| suite on Python 3.10 (uv venv, numpy + pytest only) | 2307 passed, 2 xfailed |
| `ruff format --check`, `ruff check` (geomech, tests, specs/reviews/scripts) | clean |
| `mypy geomech --ignore-missing-imports` | `Success: no issues found in 34 source files` (base: 15 errors in 6 files) |
| pyright 1.1.414, default settings | 5 errors (base fc0d386: 109; audit F16 counted 102 earlier) |
| `mkdocs build --strict` | builds (the only output is MkDocs' own "MkDocs 2" banner) |
| `tests/test_rule_names.py`, `tests/test_rule_coverage.py` | pass; rule names unchanged except the deleted `vector.transpose_transpose` |
| `engine_fuzz.py --cases 2000` | path line `fuzzing geomech from <repo>/geomech` (the working tree); 0 wrong-value findings; 178 distinct findings vs 183 at fc0d386: `flag-dependent` 99 vs 104, `unsupported-linear-term` 74 vs 74, `documented-error-on-legal-input` 5 vs 5. **Corrected (phase 5b review F4):** the difference was not deduplication by printed form. Tuple operands broke the script's shrinker (`nodes[:i] + [repl] + ...` raised a swallowed `TypeError`), and the shrinker draws from the campaign RNG, so after the first finding the adversarial campaigns generated different trees. With the shrinker fixed (measured on the fix tree, 2026-11-13): `--cases 2000` gives 183 distinct findings (104/74/5) and a per-campaign table identical to fc0d386, and the findings are identical to fc0d386's with `*` removed; `--cases 3000 --seed 55` gives 257 vs 258 at fc0d386, the one difference being `extract[vector query, matrix expr]` 89 → 88: findings there are deduplicated by the first 60 characters of the error message, and the inserted `*` moves the character that told two outer-product messages apart (`'b*D…` vs `'b*(…`) past the cut. Two runs with the same seed give identical output (timings aside) |
| examples (`examples/*.py`) | all six run |
| notebooks | re-executed in place; two output lines change (`+ Re3f` → `+ R e3 f`, see D2); `execution` metadata stripped, `language_info.version` kept |
| snapshot gates | per step in `phase5b-string-changes.md`: steps 1, 2a, 2b, 3, 4 change 0 strings; 2c/2d change 1 (literal spelling, value identical); D2 changes 1594 outputs, every one the old string with `*` inserted, all values identical, `check_eom` 0 worse |
| pickle / copy / deepcopy of a frozen tree | round trip equal (`==`, `hash`, `str`), operands stay tuples, the copy is still frozen, a `TS2` keeps its point; a tree pickled with list operands unpickles with tuples (`tests/test_phase5b_structure.py`) |

## Timing

Median of 20 runs (`linearize` 10) after one warm-up, three repetitions, same machine:

| system | fc0d386 `compute_eom` | 5b `compute_eom` | fc0d386 `linearize` | 5b `linearize` |
|---|---|---|---|---|
| quadrotor with load (`tests/vbl_systems.quadrotor_with_load`) | 12.0–12.1 ms | 12.5–12.8 ms | 7.0–7.1 ms | 7.3–7.4 ms |
| double rigid pendulum (`tests/vbl_systems.double_rigid_pendulum`) | 20.5–20.7 ms | 19.7–20.1 ms | 21.2–21.4 ms | 21.2–21.9 ms |

Within ±5%. The one-pass like-term grouping did not buy the speed-up phase 5a hoped for:
`term_key` builds a string per addend, which costs about what the pairwise comparison of
short sums cost. It removes the quadratic case for long sums (200 addends merge in one
pass, `test_like_terms_group_in_one_pass_and_keep_order`).

## Step 1: 5.5 structure (0 strings changed)

- **`geomech/core/base/constants.py`** (new): `NEG_ONE`, `NEG_TWO`, `HALF` (each a
  `Number`, printed `(-1)`, `(-2)`, `(0.5)`), `negate(expr)` (`expr * NEG_ONE`, the same
  tree as `-expr` for every type, row vectors included), `zero_like(expr)`,
  `zero_of(expr_type)`. Exported from `geomech.core.base`. Replaced:
  - `simplification.py` inline `Scalar("(-1)", value=-1, attr=["Constant"])` (twice,
    `transpose_hat_skew`/`transpose_skew`) → `NEG_ONE`;
  - `vector_identities.py` `_neg1` local → `NEG_ONE`; `_ZERO`/`_ONE` lambdas → `Zero`/`One`;
  - `manifold_rules.py` `_NEG2` → `NEG_TWO`;
  - `ibp.py` `x.t_integrate() * (-1)` (three sites) → `negate(...)`;
  - `vbl.py` `_negate` → `negate`;
  - `simplification.py` `_zero_for`, `standard_form.py` `_zero_for` and the six
    `match self.type: ... Zero/ZeroVector/ZeroMatrix` blocks of `mixins.py` → `zero_like`.
- **Duplicates** → `geomech/core/transformations/utils.py`: `addends` (was
  `standard_form._flatten_addends`, `vbl._flatten_addends`, `eom._terms`), `sum_of` (was
  `standard_form._combine`, `vbl._combine`), `flatten_mul` (was `simplification._flatten_mul`,
  `ibp._flatten_mul`, `eom._mul_factors`). `printing._scalar_factors` stays (printing
  cannot import the transformations at load time, see below).
- **Public names**: `has_algebraic`, `same_node` (`core/math/extract.py`),
  `check_variables` (`dynamics/eom.py`), `wrap_numeric` (`core/operations/multiplication.py`,
  `_wrap_numeric` kept as an alias because `tests/test_simplify_numbers.py` imports it; the
  engine itself no longer calls it). Importers updated (`vbl.py`, `standard_form.py`).
- **Dead code**, proved by a coverage run of the suite (`coverage run --source=geomech`)
  and a grep of the repository, a downstream application and the accompanying paper:
  - `vector.transpose_transpose` (`vector_identities.py`, lines 201/210 never executed:
    `Transpose.__new__` returns `x` for `Transpose(Transpose(x))`), rule and docstring
    line deleted. `simplify.transpose_transpose` is *kept*: it is reachable when the
    operand of a `Transpose` simplifies to a `Transpose` (`(I·Rᵀ)ᵀ`).
  - `is_leaf`, `has_nested_add` (`transformations/utils.py`): used only by
    `tests/test_utils.py` and re-exported from `geomech` and
    `geomech.core.transformations`. Deleted with their exports (CHANGELOG "Removed");
    `tests/test_utils.py` now tests `addends`, `sum_of`, `flatten_mul`.
  - The "`TangentVector` parsing gap" and `Matrix.__init__` no-op lines were already gone
    at fc0d386 (`VALID_ATTRS` accepts the marker). `ibp._remove_factor` tested
    `str(x) == "1"` (a user `Scalar("1")` would have matched); it now tests `x is One`.
- **Function-local imports hoisted** where no cycle exists: `transformations/utils.py`
  (`calculus_node`, `rebuild`), `simplification.py` (`pull`, `vector_rules`, `expand`),
  `multiplication.py` (`Number`, `_is_number`, `_operand_str`), `calculus.py` (`TS2`,
  `TSO3`), `extract.py` (`One`, `IdentityMatrix`), `invariants.py` (all node classes),
  `trace.py` (`print_tree`), `eom.py`/`vbl.py`/`standard_form.py` (`check_tree`),
  `printing.py` (`numpy`, `re`, `shutil`). Each module also imports alone in a fresh
  interpreter. **Kept, with the reason in a comment:** `printing.py`'s imports of the node
  classes (`geomech.core.base.expressions` imports `geomech.utils.errors`, whose package
  `__init__` loads `printing` before any node class exists), `expressions.py` → operations
  (operations import expressions), `mixins._apply_rule` → `addition`,
  `multiplication._is_row` → `geometry` (geometry imports multiplication).
- **Vocabulary**: "variation vector" on the EOM side and "perturbation" on the VBL side
  were already consistent in docstrings and the glossary (`docs/conventions.md`); the new
  code follows it. Error messages were left as they are (tests pin some).

## Step 2: 5.6 split by node class, like terms

Rule names unchanged; `test_rule_names.py`/`test_rule_coverage.py` pass.

| function | now |
|---|---|
| `simplification._eliminate` (154 lines) | `_eliminate` (`simplification.py:184`) dispatches on `type(expr)` through `_ELIMINATE` (`:711`) to `_elim_add`, `_elim_vadd`, `_elim_madd`, `_elim_mul`, `_elim_svmul`, `_elim_smmul`, `_elim_mvmul`, `_elim_mmmul`, `_elim_vvmul`, `_elim_dot`, `_elim_cross`, `_elim_hat`, `_elim_vee`, `_elim_transpose`, `_elim_variation`, `_elim_time_derivative`, `_elim_time_integral`. The three sum folds stay separate functions because `test_rule_names.py` reads the literal rule names |
| `pulling._at_node` (the `pull` match) | `_PULL` table (`pulling.py:223`): `_pull_dot`, `_pull_svmul`, `_pull_cross`, `_pull_mvmul` (row case → `_pull_row_mvmul`), `_pull_smmul`, `_pull_mmmul`, `_pull_vvmul`, `_pull_transpose`, `_pull_hat`, `_pull_vee` |
| `printing.to_latex` (119 lines) | `to_latex` (`printing.py:641`) + `_latex_handlers()` (`:667`, built on first use: the node classes cannot be imported at load time) + `_latex_leaf`, `_latex_sum`, `_latex_scaled`, `_latex_product`, `_latex_dot`, `_latex_cross`, `_latex_transpose`, `_latex_time_derivative`, `_latex_variation` |
| `collecting._collect_dot` (110 lines) | preamble + one `match` dispatching to `_collect_mvmul_cross`, `_collect_mvmul_pair`, `_collect_mvmul_other`, `_collect_triple_product` (`collecting.py:123`–`:252`); case order unchanged, so the termination argument in the docstring still holds |
| `vbl.linearize` (128 lines) | validation and setup in `linearize`; `_Blocks` (`vbl.py:336`, rows and the three block dicts with `add_row`), `_kinematic_row` (`:358`), `_dynamic_row` (`:378`) |
| `eom._extraction_error` (≈ 90 lines) | the cause selection stays in `_extraction_error` (`eom.py:656`); the three messages are `_higher_order_error`, `_ibp_left_error`, `_collect_error` |
| phase 4 A3 (context in a private attribute set on the exception) | `_separate_variations` raises `_ExtractionFailure(engine, vec)` (`eom.py:550`); `compute_eom` and `separate_variations` turn it into the engine's error type with `failure.retold(dS, ..., lagrangian, variables)`, `from` the engine error. The signature `_separate_variations(dS, vecs, record)` is kept (a test spies on it) |

Like terms (`simplification.py`):

- `term_key(expr)` (`:562`): a canonical string, injective on trees up to the swap of `Dot`
  operands, the order of a flattened `Mul` chain and, new, the swap of an inner product
  `VVMul(Transpose(a), b)` of two columns. A leaf is `Class#len#name` (self-delimiting).
  `_commutative_equal` (kept; tests import it) is `term_key(a) == term_key(b)`.
- `_collect_like_terms` (`:633`) groups addends by `_like_key(base, factors)` in one pass.
  Groups keep the order of their first term and their coefficients are summed in addend
  order, so the result is identical to the pairwise loop (0 strings changed in step 2a).
- Inner products merge (phase 5a follow-up 1): `aᵀb s − bᵀa s → 0`,
  `aᵀb s + bᵀa s → (a)'*b*s*(2)`. No recorded output changed (no output holds both).

Phase 5a advisory items (step 2c/2d, one recorded output changed, value identical):

- `_canonical_literal` (`simplification.py:448`, replaces `_exact_literal`): `_eliminate`
  puts every numeric literal leaf in the `Number` spelling with `_exact` applied, so a bare
  `Number(2.0)` simplifies to `(2)` and `Scalar("0.5", value=0.5)` to `(0.5)`. `Zero`/`One`
  keep their names; a named scalar with a value is untouched (D1).
- Phase 4 follow-up 5: `Vector.delta/t_diff` and `Matrix.delta/t_diff` of a `Constant`
  leaf return `ZeroVector`/`ZeroMatrix` (`expressions.py:352` `_zero_leaf`; a size other
  than 3/3×3 keeps a same-sized zero leaf). Before: a vector `0` and a matrix `O`.
  `tests/test_vectors.py`/`tests/test_matrices.py::test_delta_constant` pinned the old
  names and were updated (`"0v"`, `"0"`).
- `tests/test_scalar_adversarial.py:328` `test_f_is_mathematically_correct_but_uncombined`:
  `git blame` shows it was written on 2026-11-11 (commit 218f6b4e, the scalar-variables
  adversarial review), with its assertion updated in phase 5a. It is **not a duplicate** of
  `:334` `test_f_combines_like_terms`: `:328` pins the exact string, `:334` only that the
  result is not a sum. Left as is (the name is stale: it now asserts the combined form).

## Step 3: D3 freeze (0 strings changed)

- Every node class is `@dataclass(eq=False, repr=False, frozen=True)`, including the
  undecorated subclasses `S2`, `SO3`, `TS2`, `TSO3`, `SkewSymmMatrix` (otherwise a
  subclass could still set non-field attributes). Operands are `nodes: tuple = ()`, set by
  `_set_nodes` (`expressions.py:43`); leaf constructors set their fields with
  `_set_fields`. `ExprFlags` and `ManifoldInfo` are frozen; `ExprFlags.from_attr_list`
  builds the flags from a table instead of setting attributes; `Scalar` uses
  `dataclasses.replace(flags, is_numeric=True)`.
- `Expr.__setstate__` turns list operands of an old pickle into a tuple (the snapshot
  baseline was pickled at fc0d386 and replays correctly).
- Places that mutated a node: the node constructors (now `object.__setattr__` through
  the two helpers), `TS2`/`TSO3` (appended `"TangentVector"` to `self.attr`; now a new
  list), `VVMul._result_type`, `pinned_snapshot.py::_refresh` (sets `attr`/`flags` of
  unpickled identity leaves; now `object.__setattr__`), and tests that build malformed
  trees for `check_tree` (`tests/test_invariants.py`, `tests/test_phase4_adversarial.py`,
  now an explicit `object.__setattr__` bypass) or compared `nodes == [a]` (3 tests, now
  `(a,)`). A grep for attribute assignment on nodes, `.nodes.append`, `flags.is_*=` finds
  nothing else in `geomech/`, `examples/`, `notebooks/`, `docs/`.
- Not frozen, deliberately: `attr` stays a list (mutating it changes nothing: the flags
  were built at construction; `tests/test_matrices.py` relies on appending to it), and a
  leaf's `value` array is a numpy array (its entries can be written).

## Step 4: 5.7 typing (0 strings changed)

- Mixins: `_BaseMixin(_NodeBase)` where `_NodeBase` is `Expr` under `TYPE_CHECKING` and
  `object` at run time, plus `nodes: tuple[Any, ...]`; this removes every `mixins.py`
  error. `Expr` declares `type` (`TYPE_CHECKING` property) and annotates `delta`,
  `t_diff`, `t_integrate`, `get_variation_vector`, `T`, the flag properties and `has`.
- `str | None` where a default is `None` (leaf `name`, `EOMTrace.variation_vectors`/
  `ibp_targets`, `StandardFormEquation.f`); `flags` defaults to `ExprFlags()`; `S2`/`SO3`
  declare `manifold_info`, `TS2.S2`/`TSO3.SO3` are declared fields (typed `Any`: the class
  name and the attribute name coincide).
- Public functions annotated: `compute_eom`, `separate_variations`, `to_standard_form`,
  `linearize`, `check_group_is_sequence` (dynamics); `hat`, `vee`, `rand_rot`,
  `evaluate`, `sample_state`, `sample_direction`, `perturbed_state`, `linearized_rhs`
  (numeric; `finite_difference_check`, `check_eom`, `assert_eom` already were); every
  public function of `utils/trace.py`, `utils/printing.py`, `utils/invariants.py`,
  `utils/errors.py`. Annotating made mypy check those bodies; two local variable renames
  (`numeric.sample_direction`, `invariants.check_tree`) were needed.
- CI: the lint job installs `ruff mypy numpy` and runs `mypy geomech
  --ignore-missing-imports` after ruff (`.github/workflows/tests.yml`); `mypy` added to the
  `dev` extras (`pyproject.toml`).
- pyright, before → after: 109 → 5 (four `.value` on an `Expr` after an
  `is_number_literal` check in `simplification.py`, one optional access in
  `invariants.py`). Not required to be 0.

## Step 5: D2 separator (every pinned string)

- `PRODUCT_SEPARATOR = "*"` (`multiplication.py:49`): `Mul`, `SVMul`, `SMMul`, `MVMul`,
  `MMMul`, `VVMul` print `left*right`. `*` occurs in no other printed form (no baseline
  output contains one); `repr` keeps `·` between factors (and `*` as its constant
  marker), so the two forms stay distinct (`str(x*m) == "x*m"`, `repr == "v:x·m*"`).
- `Expr.__eq__` (`expressions.py:242`) compares class, `str` and `_identity_extra()`;
  `_CalcUnaryMixin._identity_extra` (`mixins.py:185`) is the chain of classes down to the
  first non-calculus operand, so `Variation(Scalar x) != Variation(Vector x)` and
  `TimeDerivative(Scalar m) != TimeDerivative(Vector m)`. `hash` is still `hash(str)`
  (equal objects have equal hashes).
- `Expr.has` (leaves) is `self == elem` (the type guard of `__eq__`); `_BaseMixin.has`
  already used `==`. `same_node` falls back to `==` instead of `str` for nameless nodes.
  The engine audit F13 reproducer (a scalar `x` next to a vector `x`) now extracts `k` and
  `F` instead of raising.
- `display_standard_form` printed an input key with `str` inside LaTeX (`+ Re3f`); with the
  separator that became `R*e3*f`, so `_latex_str_key` now turns the separator into a space
  (`+ R e3 f`). This is the only notebook output change.
- Strings regenerated as in phase 5a: snapshot at the step-4 state, apply D2, compare.
  1594 outputs changed, every one `new.replace("*", "") == old`, all values identical,
  `check_eom` 0 worse. Tests updated mechanically (318 lines in 25 files, each line its old
  self with `*`s inserted), docs and notebooks regenerated from the code; details and
  counts in `phase5b-string-changes.md`.

## New tests

`tests/test_phase5b_structure.py`, 34 tests: the constants and `negate`/`zero_like`, the
constant-variation zeros, the public names, the inner-product merge (numeric), `term_key`
symmetries and asymmetries, one-pass grouping of a 200-term sum (numeric), literal
spelling, frozen nodes (six mutation kinds, each `FrozenInstanceError` with `str`/`hash`
unchanged), tuple operands and frozen flags, pickle/copy/deepcopy, unpickling list
operands, shared literals, `check_eom` of a pinned system, `_remove_factor` with `One`,
the D2 separator and collisions, `repr` vs `str`, calculus identity, the `has()` guard,
and the F13 extraction reproducer.

## CHANGELOG

`[Unreleased]`: "Added" (constants module, shared helpers, public names, typing and mypy
in CI), "Changed" (printed form of products, equality, frozen nodes, one-pass like terms
and symmetric inner products, literal spelling, constant-variation zeros,
`display_standard_form`, the splits), "Removed" (`is_leaf`, `has_nested_add`,
`vector.transpose_transpose`).

## Stale lines in the contributor notes (not edited)

| line | now | proposed |
|---|---|---|
| venv setup note ("Recreate the venv before running anything") | the development virtual environment works (Python 3.12) | Say the existing virtual environment works and use `python -m pip`; keep the recreate commands as the fallback |
| `pip install -e ".[dev]"  # editable install + pytest, ruff, jupyterlab` | the extras include mypy | `# editable install + pytest, ruff, mypy, jupyterlab` |
| `python -m pytest tests/ -q  # full suite (~1000 tests, ...)` | 2307 tests | `# full suite (~2300 tests, numpy is the only runtime dep)` |
| `ruff format geomech/ tests/ && ruff check geomech/ tests/  # CI runs ...` and "CI (`.github/workflows/tests.yml`) lints first" | CI also runs mypy | add `mypy geomech --ignore-missing-imports` to the command block; "CI lints first (`ruff format --check`, `ruff check`, `mypy geomech --ignore-missing-imports`), then tests on 3.10–3.12" |
| "Everything is an immutable **expression tree**." | now enforced | "Everything is an immutable **expression tree**: node classes are frozen dataclasses, operands are a tuple (`expr.nodes`), flags are frozen; build a new node instead of changing one (`object.__setattr__` is the only way past it, used by constructors and by tests that build malformed trees)." |
| layer 1 `core/base/` | has `constants.py` | add: "`constants.py` holds `NEG_ONE`, `NEG_TWO`, `HALF`, `negate(expr)`, `zero_like(expr)`; rewrites build −1, −2 and typed zeros from them." |
| layer 2: "`mixins.py` holds the operator dispatch (`__add__`, `__mul__`, …) that picks the right node class from `(self.type, other.type)` and raises `ExpressionMismatchError` otherwise." | stale since phase 3: the operators live on `Expr` (`core/base/expressions.py`, `_add`/`_sub`/`_mul`) | "`Expr` (`core/base/expressions.py`) holds the operator dispatch (`__add__`, `__mul__`, … via `_add`/`_sub`/`_mul`), which picks the node class from `(a.type, b.type)`; `mixins.py` holds what node families share (`has`, `left`/`right`/`expr`, the product rule of `delta`/`t_diff`, the calculus identity)." |
| layer 3 names only functions | modules have other names | "`core/transformations/` — modules `expansion` (`expand`), `pulling` (`pull`), `vector_identities` (`vector_rules`), `simplification` (`simplify`, `full_simplify`; `_eliminate` dispatches per node class through `_ELIMINATE`, like terms are keyed by `term_key`), `manifold_rules`, and `utils` (`rebuild`, `map_children`, `addends`, `sum_of`, `flatten_mul`)." |
| layer 4 | modules | "`core/math/` — `extract.py` (`extract_linear_coeff`, `has_algebraic`, `same_node`), `collecting.py` (`collect`), `ibp.py` (`integrate_by_parts`)." |
| "`Expr.__eq__`/`__hash__` compare `str()` after a type check. Commutative rearrangements are therefore *not* equal; an earlier code review lists this as a known limitation." | products print `*`; calculus nodes compare their operand class; `has()` uses `==` | "`Expr.__eq__` compares class and `str()` (and, for a calculus node, the operand classes); `__hash__` is `hash(str)`. Products print `*` between factors, so two different trees never print alike (D2); `has()` uses the same `==`. Commutative rearrangements are still *not* equal (`a*b` vs `b*a`); like-term merging uses `term_key`, which knows `Dot`, `Mul` and inner products commute." |
| Testing conventions | pinned strings carry the separator | optional: "Pinned strings print products with `*` (`m*l*l*(-1)`). A change that reprints them goes through `specs/reviews/scripts/pinned_snapshot.py` (snapshot before, compare after, numeric equality)." |

## Notes for the reviewer

- The accompanying paper and a downstream application were not touched. Any listing in the
  paper produced by `str` now carries `*`; LaTeX produced by `to_latex` is unchanged.
- The `*` collides with nothing in `str`, but it is also `repr`'s constant marker; the two
  forms are never mixed in one string.
- `term_key` is a string per addend; a cached key per node (possible now that nodes are
  frozen) would make like-term grouping cheaper. Not done: no measurable gain at today's
  sizes.
- `vector_rules._at_node` (80 lines) was not split; the task list did not include it.
- Remaining pyright errors (5) and the `Any`-typed `TS2.S2`/`TSO3.SO3` are listed above.
