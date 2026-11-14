# Phase 5b adversarial review

Branch `fix/phase5b-structure`, commit f852c11, range `fc0d386..f852c11`. Base for every
differential run: a detached worktree of fc0d386 (removed afterwards) sharing the repo
venv (Python 3.12.13), with `geomech` imported from the worktree (checked per run).
New tests: `tests/test_phase5b_adversarial.py` (20 pass, 7 strict xfail). Nothing
committed; no file under `geomech/` and no existing test was changed.

## Verdict: GO WITH CHANGES

No silent behavioural change in the engine was found. The value-level claims hold on
every differential run below. The changes needed before merge are small and in the tests,
tooling and docs: two negative assertions that the separator made vacuous (F1), the fuzz
shrinker broken by tuple operands (F4), and a decision on `*` in symbol names (F2), which
reopens the F13 collision that D2 was meant to close.

## Findings, ranked

| # | rank | finding | where | reproducer |
|---|---|---|---|---|
| F2 | equality hole (medium) | Leaf names may contain `*`, so a product and a lumped symbol print alike and compare equal again: `Mul(Scalar("m*l"), g) == Mul(m, Mul(l, g))` (both `m*l*g`). This is engine-audit F13 moved from names like `ab` to names with `*`. In a path that uses `==`, the value is wrong: `_t_diff_product(Mul(f, g))` with `f = Dot(x*m*l, y)` and `g = Dot(x*L, y)` (`L = Scalar("m*l")`) returns `2 f ḟ` instead of `ḟ g + f ġ`. Like-term merging is not fooled (`term_key` writes leaves as `Class#len#name`), but `full_simplify(m*l*g − L*g)` then prints `(m*l*Dot(x,y)+m*l*Dot(x,y)*(-1))`, which reads as zero and is not. Names with `*` occur in practice: `tests/test_latex_precedence.py:112` uses `Vector("x^*")`. Suggested fix: reject the separator in leaf names in `_attr_list`/the leaf constructors (or escape it in `__str__`), and correct the claim in `docs/conventions.md` ("two different products never print alike") and the `extract.same_node` docstring | `geomech/core/operations/multiplication.py:59` (`PRODUCT_SEPARATOR`), `geomech/core/math/ibp.py:108` (`f == g`) | `test_f2_separator_in_a_name_does_not_collide_with_a_product`, `test_f2_t_diff_product_keeps_a_lumped_symbol_apart` |
| F1 | lost assertion coverage | Two negative substring assertions were not updated and can no longer fail. `tests/test_manifold_pipelines.py:485` `assert "(R)'R" not in str(self.eq.f)`: RᵀR now prints `(R)'*R`. `tests/test_vbl.py:651-653` `"(-1)(-1)"`, `"(-1.0)(-1"`, `"(-1)(-1.0"` not in the block strings: chained numbers now print `(-1)*(-1)`. The string-update tooling only rewrote *failing* comparisons, so it could not see these. The intended checks still hold with the separator (passing tests `test_f1_rigid_pendulum_f_has_no_RtR_with_the_separator` and `test_f1_vbl_blocks_have_no_unfolded_numbers_with_the_separator[8 systems]`). Fix: insert `*` in the four literals. No other vacuous negative was found (every `not in`, `!=`, `count`, `re.*` on a string literal in `tests/` was checked by hand) | the two test files | `test_f1_manifold_pipelines_negative_check_uses_the_separator`, `test_f1_vbl_negative_check_uses_the_separator` |
| F4 | lost coverage (tooling) | `specs/reviews/scripts/engine_fuzz.py:853` `shrink_candidates` builds `nodes[:i] + [repl] + nodes[i + 1 :]`. `nodes` is now a tuple, so this raises `TypeError`, which the surrounding `except Exception: pass` swallows. The fresh-leaf shrink candidates are silently gone (`Dot(v*a*b, v)` gives 3 candidates instead of 5). Detection is unaffected, but reproducers shrink less. The shrinker also consumes the campaign RNG, so a run with the same seed generates different trees after the first finding. With `--seed 55` the first divergence is case 1 (`vector_rules[adversarial]`): same output and same status, different RNG state afterwards. That, and not dedup by printed form, is why the adversarial counts differ between fc0d386 and 5b (`vector_rules[adversarial]` ok 913 → 912, `simplify` 846 → 849, `full_simplify` 846 → 848). The implementation report's attribution of 183 → 178 findings (and my 258 → 264 at `--cases 3000 --seed 55`) should be restated. Fix: `list(nodes[:i]) + [repl] + list(nodes[i + 1 :])` | the script | `test_f4_engine_fuzz_shrinker_offers_fresh_leaf_candidates` |
| F3 | equality hole (low, not a regression) | `Transpose` preserves its operand's type, as the calculus nodes do, but `_identity_extra` was added only to `_CalcUnaryMixin`. So `Transpose(Vector("x")) == Transpose(Matrix("x"))` (types VECTOR and MATRIX, equal hash), and `Variation(Transpose(·))` inherits this, because the class chain stops at the first non-calculus node. No engine path was found where this changes a value: `term_key` keeps the two apart, and configuration names are unique across types. Fix: give `Transpose` the same `_identity_extra` | `geomech/core/operations/geometry.py:127`, `mixins.py:185` | `test_f3_transpose_compares_the_operand_class` |
| F5 | cosmetic | `_latex_str_key` turns every `*` of a key into a space, so `display_standard_form` renders an input named `u^*` as `u^ `, which is invalid LaTeX. Fix: render the key from the input expression with `to_latex` rather than from its `str` | `geomech/utils/printing.py:861` | `test_f5_latex_key_keeps_a_star_in_a_name` |

Advisory (no test):

- `_ELIMINATE`, `_PULL` and `_latex_handlers()` dispatch on the exact `type(expr)`, where
  the old `match` accepted subclasses. No subclass of a node class exists in `geomech/`,
  `tests/`, `examples/` or a downstream application. A future subclass would fall through to
  leaf handling (`_canonical_literal`, `str`) without an error.
- CI installs `mypy` unpinned (`pip install ruff mypy numpy`). It passes today with mypy
  2.4.1 and `--python-version` 3.10, 3.11 and 3.12. Pinning the version avoids a
  surprise lint failure.
- The `extract.same_node` docstring cites `Scalar("k*x")` as the leaf the type guard keeps
  apart from `Mul(k, x)`. That is true one level down, but a parent node collides (F2).

## What held (with counts)

1. **Pure-separator claim (tests).** I compared every changed test-line pair independently
   (`git diff -U0`, pairing within hunks, checking that the old line equals the new one with
   some `*` removed, and that each inserted `*` lies inside a string token). Results:
   - **305** line pairs are pure insertions inside strings.
   - **12** reflowed hunks are pure: their joined literals are equal with `*` removed and
     their code is equal.
   - The remaining changes are the documented ones: `nodes == (a,)` (3), the
     `object.__setattr__` bypass for malformed trees (`test_invariants.py` 11 lines,
     `test_phase4_adversarial.py` 2), the zero-leaf names `0v`/`0` (2), and `test_utils.py`
     (13 tests of deleted helpers replaced by 8 tests of `addends`/`sum_of`/`flatten_mul`).
   - Assertion and test counts per file are unchanged except `test_utils.py` (13 → 12
     asserts, deleted functions) and the new `test_phase5b_structure.py`.
   - No comparison operator, expected count or tolerance changed.
   - The only lost coverage is F1.
2. **Equality.**
   - Pairs built two ways are still equal: operator vs constructor for every product,
     `-x` vs `x*-1` for each type, `v.T()*w` vs `VVMul`, and calculus of equal operands.
   - The F13 pairs (`Mul(ab,c)`/`Mul(a,bc)`, the `MMMul` analogue, `Variation`/`d/dt` of
     a scalar vs a vector `x`) are unequal, and a dict keyed on `Variation(x)` misses
     `Variation(X)`.
   - The `str(` lookups in the engine (`manifold_rules` tangent map, `eom` reserved
     names, stray and allowed sets, `eom[str(vec)]`, `standard_form` M/G keys, the `vbl`
     blocks, the `numeric` env) are all keyed on leaves or calculus-of-leaf nodes, and
     both sides are built with `str`. None can miss because of the separator.
   - Configuration variables with one name across types are rejected
     (`SystemVariables: the name 'x' is listed more than once`).
3. **Like terms.**
   - 3000 seeded sums (scalar, vector and matrix; inner and outer products, aᵀMb,
     Cross, MVMul, MMMul, Dot, Mul chains with numeric and symbolic coefficients, swapped
     variants) were run through old and new `_collect_like_terms`, `simplify` and
     `full_simplify`. There was **0 value difference** against the input, and 0 between
     old and new.
   - 270 `_collect_like_terms` and 191 `simplify` outputs differ in string. All of them
     contain a transpose, and each is a merge: the same number of addends or fewer.
     Sums without a transpose: 1516 identical out of 1516.
   - The merge fires only for `VVMul(Transpose(column), column)` of scalar type. Outer
     products are not merged, aᵀMb vs bᵀMa is not merged (symmetric J included,
     conservatively), and (Ma)ᵀb + bᵀ(Ma) is merged. All pinned in
     `test_inner_product_merge_is_value_preserving_and_restricted`, plus 300 extra seeded
     sums in `test_like_term_grouping_random_sums_keep_their_value`.
4. **Frozen nodes.**
   - The only `object.__setattr__` calls in the library are in `_set_fields`/`_set_nodes`
     (constructors) and `VVMul._result_type` (constructor).
   - There is no `cached_property`, `lru_cache`, `__post_init__` or misused
     `dataclasses.replace` (only `replace(flags, is_numeric=True)`).
   - `__setstate__` converts list operands to a tuple.
   - Pickle and deepcopy round trips keep `==`, `hash`, `str`, tuple operands, frozenness
     and `VVMul.type`.
   - Shared `NEG_ONE`/`HALF` constants cause no aliasing.
   - Timing, median of 15, two alternating rounds, against fc0d386:
     - quadrotor with load: 15.6/17.3 ms vs 17.1/15.9 ms
     - double pendulum: 24.1/26.3 ms vs 26.9/27.1 ms
     - quadrotor: 3.7/3.8 ms vs 3.4/3.7 ms
     - `linearize` is within noise.
   - Tuple operands cost nothing measurable.
5. **Function splits.**
   - On the engine_fuzz generator, 5000 trees (seed 55, depth 1–4), `expand`, `pull`,
     `vector_rules`, `simplify`, `full_simplify` and `to_latex` gave **5000/5000 identical**
     outputs (strings with `*` removed, LaTeX byte for byte).
   - On the same trees, `collect` gave 1932/1932 identical outputs, with 20 identical
     errors.
   - Every output's value matched old at fresh jets.
   - The `full_simplify` rule-name sequences under `set_trace` were 5000/5000 identical.
   - On the 67 builder systems (`pinned_snapshot.builder_systems`), the `compute_eom`,
     `to_standard_form` and `linearize` traces had **201/201 identical rule-name
     sequences** (9271 rule lines).
   - Line content differs only by the documented `(0.5)` spelling and the
     `0v`/`0`-for-`O` zero leaves (including where they shift `_cut_pair` truncation).
     EOM LaTeX is identical for all 67 systems, and outputs are equal with `*` removed.
   - Error messages of the 5 failing builder calls are separator insertions.
6. **mypy in CI.**
   - `.github/workflows/tests.yml` runs `mypy geomech --ignore-missing-imports` in the
     lint job, on Python 3.12 only, before the 3.10–3.12 test matrix.
   - Locally it gives `Success: no issues found in 34 source files` with
     `--python-version` 3.10, 3.11 and 3.12.
   - The suite passes on a real 3.10.20 interpreter (uv venv, numpy 2.2.6): 2267 passed,
     2 xfailed, `-W error`.
7. **Regression.**
   - Suite `-W error`: 2267 passed, 2 xfailed. With `GEOMECH_DEBUG=1`: the same.
   - With the new file: 2287 passed, 9 xfailed.
   - `ruff format --check` and `ruff check` on `geomech`/`tests`: clean.
   - `engine_fuzz.py --cases 3000 --seed 55 --no-json`: 0 wrong-value findings, old and
     new. The counts differ only in the adversarial campaigns (F4).
   - The 6 examples run on both trees with identical output (with `*` removed).
   - The 4 notebooks re-execute with outputs identical to the committed ones.
8. **Docs.**
   - README, `pipeline.md`, `conventions.md`, `troubleshooting.md` and `vbl.md` show the
     `*` form throughout. No juxtaposed product remains in a printed-output example.
   - `tests/test_troubleshooting_doc.py` passes (the quoted messages match).
   - The one inaccurate statement is the "never print alike" claim (F2).

## Scripts (temporary directory, not committed)

`rev5b/sepcheck.py` (test-line pairs), `rev5b/diffeng.py` + `cmpeng.py`/`cmpsum.py` (5000
trees, 3000 sums, old vs new), `rev5b/sysdiff.py` + `cmpsys.py` (67 systems, rule traces),
`rev5b/advtrace.py` (F4 RNG divergence), `rev5b/bench.py` (timing).

## Fix log

Fixes on top of 3b45974 (uncommitted). Every strict xfail of
`tests/test_phase5b_adversarial.py` was flipped by deleting its marker; no assertion of an
xfail was changed.

- **F1** (orchestrator, before this log): the four literals in
  `tests/test_manifold_pipelines.py` and `tests/test_vbl.py` now contain `*`
  (`"(R)'*R"`, `"(-1)*(-1)"`, `"(-1.0)*(-1"`, `"(-1)*(-1.0"`). Markers removed; both pass.
- **F2**: `PRODUCT_SEPARATOR` moved to `geomech/core/base/expressions.py` (re-imported by
  `multiplication.py`), and `_check_leaf_name` rejects a name containing it in
  `Scalar.__init__`, `Vector.__init__` and `Matrix.__init__`, which every leaf class
  (`S2`, `SO3`, `TS2`, `TSO3`, `SkewSymmMatrix`) and `Number(str)` go through:
  `ValueError: symbol name 'm*l' contains '*', which is reserved as the product
  separator; ...`. Grep of `geomech/`, `tests/`, `examples/`, `notebooks/`, `docs/`, README
  and a downstream application for a leaf name with `*`: only `tests/test_latex_precedence.py:112`
  (`Vector("x^*")`, phase 3, 2026-11-13), renamed to `Vector("x^\\ast")` with the expected
  LaTeX `\left(x^\ast\right)^T` (the test still checks that a superscripted name is not
  given a second superscript), plus the review's own tests. No pre-2026-11-11 test or
  example uses `*` in a name. `docs/conventions.md`, the `Expr.__eq__` and
  `extract.same_node` docstrings now state the guarantee precisely: every `*` in a printed
  form is a product boundary, so two products print alike only when they multiply factors
  that print alike in the same order (regroupings `(a*b)*c`/`a*(b*c)` print alike and are
  equal). The two F2 xfails now pass through their `except ValueError: return` branch.
  `test_f2_like_terms_keep_a_lumped_symbol_apart` (a plain test) built `Scalar("m*l")`,
  which the constructor now rejects: it builds that leaf around the constructor
  (`__setstate__` of a `Scalar("ml")` with the name replaced) so it keeps checking that
  `term_key` does not depend on the reservation; its assertions are unchanged.
- **F3**: `Transpose._identity_extra` returns the operand's class plus the operand's own
  extra, and `_CalcUnaryMixin._identity_extra` appends the extra of the first non-calculus
  operand, so `Transpose(Vector("x")) != Transpose(Matrix("x"))` and
  `Variation(Transpose(·))`/`TimeDerivative(Transpose(·))` differ too; equal trees keep an
  equal hash. `Hat` and `Vee` are unchanged: their operand has one fixed type (VECTOR,
  MATRIX), so their printed form already fixes the operand's type. They do not tell a
  `Vector("x")` operand from a `TS2("x")` one, but neither does any binary node
  (`Dot(Vector x, w) == Dot(TS2 x, w)`); that is the leaf-name rule of
  `docs/conventions.md` (one printed name, one symbol), not a type hole.
- **F4**: `shrink_candidates` builds `(*nodes[:i], repl, *nodes[i + 1 :])` (and
  `(*nodes[:i], *nodes[i + 1 :])` for the addend drop, which worked on tuples but is now
  written the same way). No other list concatenation on node operands in
  `specs/reviews/scripts/`. `engine_fuzz.py --cases 3000 --seed 55 --no-json` twice: identical
  output except the timing lines; 257 distinct findings, 0 wrong-value; adversarial
  campaigns `vector_rules` ok 913, `simplify` 846, `full_simplify` 846, the same as fc0d386.
  fc0d386 (its own script) gives 258: the one difference is
  `extract[vector query, matrix expr]` 89 → 88, because findings there are deduplicated by
  the first 60 characters of the error message and the inserted `*` pushes the
  distinguishing character of two outer-product messages past the cut. `--cases 2000`
  (default seed): 183 findings (104/74/5), per-campaign table identical to fc0d386 and
  findings identical with `*` removed (f852c11 with the broken shrinker: 178). The
  implementation report's row is corrected with these numbers.
- **F5**: `_latex_str_key` turns a `*` into a space only where it can be the separator:
  not at the end of the key and not after `^`, `_` or `{` (no factor prints empty or ends
  with those). `u^*` stays `u^*`; `R*e3*f` → `R e3 f`. With F2 an engine-built key has no
  other `*`. `display_standard_form` on the quadrotor and quadrotor-with-load systems: the
  input keys are `f` and `M` (unchanged), and the coefficient `R e3 f` comes from
  `to_latex`.
- **Advisory, mypy**: pinned to the venv's 2.4.1 in the `dev` extras and in CI
  (`pip install ruff "mypy==2.4.1" numpy`).
- **Advisory, dispatch**: exact-class lookup kept (an MRO walk would cost a walk per leaf, or
  a cache). `_PULL`, `_ELIMINATE` and `_latex_handlers()` carry a comment that a subclass of
  a node class is not dispatched, and `test_dispatch_tables_match_the_exact_class_only`
  (`tests/test_phase5b_structure.py`) pins it: `simplify`/`pull` return a `Dot` subclass
  unchanged, `to_latex` prints its `str`.
- New tests in `tests/test_phase5b_structure.py`: the name rule for every leaf class and `Number`,
  `ml`/`x^\ast` accepted, `Transpose` operand class (incl. under `Variation` and
  `TimeDerivative`, and equal trees with equal hashes), the dispatch test. CHANGELOG
  (Unreleased, Changed) records the name rule, the `Transpose` equality and the mypy pin.

Verification:

- Pinned snapshot against f852c11 (`pinned_snapshot.py`, base taken in a worktree of
  f852c11): 563 records, 4273 outputs unchanged, **0 strings changed**, 152 expected errors
  kept, 0 failures, 217 `check_eom` compared, 0 lower.
- Suite `-W error`: 2307 passed, 2 xfailed; with `GEOMECH_DEBUG=1`: the same.
- `ruff format --check`, `ruff check` on `geomech/ tests/ specs/reviews/scripts/`: clean.
  `mypy geomech --ignore-missing-imports` (3.10, 3.11, 3.12): no issues in 34 files.
- `mkdocs build --strict`: builds.
- The six examples: output byte-identical to f852c11. The four notebooks re-executed in a
  copy: every cell output identical to the committed one.
