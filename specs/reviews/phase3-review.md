# Phase 3 adversarial review: operators, printing

Branch `fix/phase3-operators-printing`, commit `475c596`, range `1e10349..475c596`.
New tests: `tests/test_phase3_adversarial.py` (31 pass, 12 strict xfail). Nothing committed.

## Verdict: GO WITH CHANGES

I found no silent wrong equation. Every operator-built system I derived passes
`check_eom`. Outside explicit double transposes, `str` held across 4500 random trees and 6
transforms, and no two random trees with equal `str` or equal LaTeX evaluate differently.
The findings below need fixing before merge:

- F1 is a small `str` change on an existing operator path, and the CHANGELOG describes it
  wrongly.
- F2 is a dead end on an input that the new docs and error messages recommend.
- F3 is LaTeX that reads as a subtraction.
- F6 and F7 are notebook hygiene issues.
- F4 and F5 are small.

## Findings, ranked

### F1: `s * vᵀ` builds a different tree (equality change)

- **Reproducer:** `str(x * a.T())`, and likewise `Number(2) * a.T()`.
  - At `1e10349`: `(a)'x`, which is `SVMul(Transpose(a), x)`. This is one transpose and a
    correct row.
  - Now: `(ax)'`, which is `Transpose(SVMul(a, x))`.
  - Consequence: `x * a.T() == SVMul(Transpose(a), x)` was True and is now False.
- **Cause:** `geomech/core/base/expressions.py:106` (`_scale`) rewrites both orders.
- **Plan scope:** only the `Transpose(v) * s` path built the double transpose that the plan
  set out to fix.
- **CHANGELOG is wrong.** It says `s * Transpose(v)` built the double transpose. It did not.
- **Why the stage dump missed it:** no pipeline path builds `s * vᵀ` through the operator.
- **Second effect (expected: cancel; observed: two terms).** The two row forms now coexist
  and do not merge: `full_simplify((x*a.T())*b - x*(a.T()*b))` → `((ax)'b+x(a)'b(-1))`.
- **Expected:** either
  - keep `SVMul(Transpose(v), s)` for a scalar on the left. `is_row_vector` already
    classifies it as a row, so it multiplies fine, and the old pinned tests only cover the
    transpose-on-left order; or
  - choose one canonical row form, enforce it in `simplify`, and correct the CHANGELOG.

### F2: row times matrix (`vᵀ M`) is accepted but the engine cannot use it (crash on legal input)

- **Where it is recommended:**
  - `README.md` and `getting-started.md` advertise `x.T() * y`.
  - The new MVMul error message says "write … `v.T() * M` for the row vᵀ M".
- **Reproducer:** with `K` symmetric and constant, `L = m*Dot(v,v)/2 - q.T()*K*q/2` and
  `dW = Dot(q.delta(), f)`, `compute_eom` raises
  `AlgebraicError: _vector_from_vector: unsupported node Transpose contains '\delta{q}' in '((\delta{q})'K)''`.
  The same happens for:
  - `x.T()*A*x` and `x.T()*A.T()*x`;
  - `(x-e3).T()*K*(x-e3)` and `(x.T()-e3.T())*K*(x-e3)`;
  - `dW = -(dx.T()*A*v)`.
- **Forms that work:** `x.T()*(K*x)` and `(A*x).T()*x` pass the oracle (err ≈ 4e-9).
- **`expand` does not distribute the new row shapes either** (`expand.py:56-96`). Examples:
  - `(a.T()*(K+J))*b` → `(a)'(K+J)b`
  - `((a+b).T()*K)*c`
  - `VVMul(SVMul(Transpose(a+b), x), c)`
  - `(a.T()+b.T())*c`
  - `d/dt((a+b)ᵀ) c`

  They all come back unchanged, and then `extract_linear_coeff` (`extract.py:412`) raises.
- **Before phase 3:** these raised `ExpressionMismatchError` at construction, so this is
  not a lost capability. But the failure has moved from a clear constructor error to a
  deep internal one, on a path the docs now recommend.
- **Expected:** one canonical row form. Rewrite `MVMul(vᵀ, M)` → `Transpose(MVMul(Mᵀ, v))`
  and `SVMul(vᵀ, s)` → `Transpose(SVMul(v, s))`, in `simplify`/`expand` or at extraction,
  so that the existing `Transpose(column)` paths apply. Until that is done, stop
  recommending `v.T() * M` in the error message.

### F3: a negated transpose inside a product prints as a subtraction (misleading LaTeX)

| Input | `to_latex` | Reads as |
|---|---|---|
| `x * (-a.T() * b)` | `x -a^T b` | x − aᵀb (a well-typed, different scalar) |
| `A * (-K).T()` | `A -K^T` | A − Kᵀ (a well-typed, different matrix) |
| `c * -a.T()` | `c -a^T` | |
| `-((-K).T())` | `--K^T` | |
| `VVMul(-a.T(), b) * -1` | `--a^T b` | |

- **Cause:** `printing.py:647-654`. A `Transpose` of an `SVMul`/`SMMul` returns
  `_product_latex(...)`, which can start with `-`. Neither `_needs_parens_in_product`
  (`:536`), `_scalar_factor_latex` (`:496`) nor the SVMul/SMMul base (`:616-620`) wraps an
  operand whose LaTeX starts with `-`.
- **Scope of the fuzz:** the 2000 random trees from `engine_fuzz.Gen` never contain a
  negated transpose. My operator-built generator hit it in 4 of 2000 trees for `--`, and
  more often for the subtraction reading.
- **Expected:** wrap any non-leading operand or base whose LaTeX starts with `-`, or fold
  the sign of a transposed scalar multiple into the enclosing coefficient.

### F4: `sum([a, b, c])` over vectors raises (crash on legal input, low)

- **Observed:** `ExpressionMismatchError: Add cannot combine SCALAR '(0)' and VECTOR 'a'`.
  `__radd__` wraps the start value `0` as `Number(0)`.
- **Scalars:** `sum([x, y, z])` works but keeps the `(0)` term, giving `((0)+x+y+z)`.
- **Expected:** `__radd__` returns `self` when `other == 0` (int). Otherwise the message
  should suggest `sum(..., ZeroVector)`.

### F5: two error messages quote an operand the user did not write (misleading)

- **The first addend instead of the operand.** `(Dot(a,a) + Dot(b,a)*x) + a` reports
  `Add cannot combine SCALAR 'Dot(a,a)' and VECTOR 'a'`. The cause is `addition.py:51`,
  which quotes `nodes[0]` of the flattened sum instead of the left operand. With a long
  Lagrangian this points at an arbitrary first term.
- **The internal form of `a - 1`.** It reports
  `VAdd cannot combine VECTOR 'a' and SCALAR '(-1)'`, the internal `x + (-1)` form
  (`expressions.py:173`). Expected: `Sub cannot combine VECTOR 'a' and SCALAR '(1)'`.
- **Cosmetic, same area:**
  - `Dot([big], x)` puts an untruncated `repr` of the list in a `TypeError` (1076 characters).
  - `_raise_constant_integral` does not truncate long names.

### F6: `print_tree` prints memory addresses (nondeterministic output)

- **Where:** `printing.py:407`, `of=S2 'q' @0x…`.
- **Effect:**
  - `notebooks/rigid_pendulum.ipynb` and `spherical_pendulum.ipynb` change on every run. I
    re-executed both; after masking `@0x…` they are byte-identical to the committed outputs.
  - The address is also misleading. The engine compares tangent parents with `==` by name
    (`vector_rules.py:43,45,111`; phase 2 made `S2("q") == S2("q")`). Two different
    addresses therefore suggest "different points" for objects the engine treats as the
    same point.
- **Proposed deterministic alternative:** print `of=S2 'q'`. Within one `print_tree` call,
  append `#1`, `#2` (order of first appearance) only when two *distinct objects* with the
  same name occur. That keeps the debugging signal without nondeterminism.

### F7: the quadrotor notebooks keep the pre-phase-3 LaTeX (stale committed output)

- **What the report says:** `quadrotor.ipynb` and `quadrotor_payload.ipynb` were left
  unexecuted.
- **What they still show:** the ambiguous forms that 3.5 removed, for example
  `-\Omega_{R} \times J \Omega_{R}` and `\dot{x} -l \omega_{q} \times q`. Executed now they
  differ in 5 and 6 cells.
- **Expected:** regenerate them with the same nbconvert command.

### Cosmetic (no test)

- Division brings float coefficients into like-term merging:
  `full_simplify(x*0.1*3 - x*0.3)` → `x(5.55e-17)`. This is pre-existing folding, made
  easier to reach by `/`.
- `simplify(Transpose(SVMul(Transpose(a), x)))` stays `((a)'x)'`, where the target is
  `ax`. `to_latex` already prints `x a` correctly.
- A named scalar with a value prints its value, and the value now folds:
  `to_latex(g*m)` with `g=9.81`, `m=2` gives `19.62`. Before it was `9.81 2`.
- `tree_str` (topdown) labels `S*V`/`M*V` use `*` as multiplication, next to leaves where
  `*` marks a constant. `repr` switched to `·`; the topdown labels did not.
- Pre-existing, outside the phase:
  - `ZeroMatrix` and `Zero` both print `0`, so `TimeDerivative(Zero) == TimeDerivative(ZeroMatrix)`.
  - `VAdd(aᵀ, b)` mixing a row and a column is accepted. Its orientation now depends on
    which addend comes first (`is_row_vector`).
  - `dataclasses.replace(Transpose(a))` raises `TypeError`, as it did at base. The engine
    does not use `replace`.

## What held

- **`str` invariance.**
  - Base vs new on 4500 `engine_fuzz.Gen` trees (seeds 31-33): the tree and `expand`,
    `pull`, `vector_rules`, `simplify`, `full_simplify` and `apply_manifold_rules`.
  - 59 of 4500 records differ: 57 contain an explicit `Transpose(Transpose(·))` in the
    input (the planned collapse of 3.2), and 2 raised `ExpressionMismatchError` at base and
    now succeed.
  - 0 unexplained differences.
- **Operator vs constructor:** 12 pairs build trees with equal `str` and compare equal.
  - `x*-1`/`-x`, `x-y`/`x+y*-1`, `Rᵀ*v`/`MVMul(Transpose(R), v)`, `x/2`/`x*0.5`
  - `x-1`/`x+(-1)`, `np.float64(2)*x`/`2.0*x`, `x*np.int64(2)`/`x*2`, `a*x`/`x*a`
  - `Number(2)*a`, `-R`, `-Dot`, `aᵀ*s`

  By design these differ: `x-1` vs `x-Number(1)`, and `2*x` (`(2)x`) vs `x*2`.
- **Numeric folding:** `x-1+1`, `1-1+x`, `x*2 - x*2.0`, `x/2 + x/2`, `-(-x)` all reduce.
- **Rewrite rules keyed on old structure:**
  - Every `isinstance(…, Transpose)` left in the engine is a matrix-transpose rule or the
    `Transpose(VAdd)` split, not a row test: `simplify.py:147,150,196`, `expand.py:64`,
    `vector_rules.py:195`, `numeric.py:211`. Row tests are all `is_row_vector`.
  - The `.expr` accesses on VVMul operands in `extract.py` were all migrated.
- **`Transpose.__new__`:**
  - `copy`, `deepcopy` and `pickle` round-trip.
  - `Transpose(Transpose(row))` returns the same row object.
  - row × matrix types as a row (evaluates to shape 1×3).
  - The `_rebuild` paths in `simplify`/`expand`/`extract` only rely on the collapsed value,
    which is mathematically the same. The `extract` rewrite `left.expr → Transpose(left)`
    is the identical object for a `Transpose` node.
- **Operators:**
  - numpy arrays raise a clear `TypeError`, on either side.
  - `x/0.0` and `x/np.float64(0)` raise `ZeroDivisionError`.
  - `bool` and `np.bool_` are rejected on every operator.
  - complex, `Fraction` and `Decimal` are rejected with a clear message.
  - Matrix and vector division by numbers works.
  - `-Zero` simplifies to `0`.
  - `x/nan` and `x*inf` are accepted (cosmetic).
- **Operator-built systems:** 10 systems, all through `compute_eom` and `check_eom` with
  err < 1e-8, including:
  - `(k*q).T()*q/2`, `k*q.T()*q/2`, `q.T()*(A*q)/2`, `(A*q).T()*q`;
  - `dW = f.T()*dq`, `2*f.T()*dq/2`, `(-c v)ᵀ dq`, `dqᵀ(-c v)`;
  - `1 - …` and `-(-KE)`.
- **Error messages:**
  - Operands longer than 80 characters are truncated in every `ExpressionMismatchError`
    and operator `TypeError` (messages are 51-264 characters for a 720-character sum).
  - No message path calls `str` on something that raises.
  - `ExpressionMismatchError` keeps base `Exception` and the 3-argument positional
    signature with the old text.
  - `UndefinedCaseError` removal: no user outside `specs/` (including a downstream application).
- **LaTeX fuzz:** 4000 trees, 2000 from `engine_fuzz.Gen` and 2000 from an operator-built
  generator with rows, negations, numeric coefficients, `/` and `TimeIntegral`.
  - 0 unbalanced `{}`/`()`; `\left`/`\right` counts all match.
  - 0 `^T^T` or `^T ^`.
  - 2787 unique strings compiled with `pdflatex -halt-on-error` in chunks of 250: 0 failures.
  - 56 groups of structurally different trees share a LaTeX string. All 55 evaluable groups
    agree numerically (shared jet environment), so there is 0 equal-LaTeX/different-value.
  - `--` appears 4 times (F3).
  - Equal-`str`/different-value: 0 introduced by phase 3. The 1 found is the pre-existing
    `Zero`/`ZeroMatrix` name clash.
- **`repr` and `print_tree`:** no raise on 6000 random trees or on any leaf type (`TS2`
  without a parent, `SkewSymmMatrix`, `IdentityMatrix`, numerics, valued vectors and
  matrices). The `·`/`*` legend is documented in `README.md`, `getting-started.md` and
  `repr_str`.
- **Regression:**
  - `engine_fuzz.py --cases 3000 --seed 31`: 189 findings, against 202 at base with the same
    seed. All 13 removed findings are the old "documented error on legal input" rows in
    simplify/full_simplify plus extract rows that shifted with the random stream. 0
    wrong-value, the same as base.
  - The 6 examples exit 0 under `-W error`.
  - The 4 notebooks execute with 0 error outputs.
  - The full suite with `-W error`: 1878 passed, 17 xfailed (1847 + 5 before my file).
  - `ruff format --check` and `ruff check` are clean.

## Fix log

Fixes on top of `d4349ac`; nothing committed. Two of the decided fixes (F1 tree, F6) are
**blocked**. In each case a phase-3 test pins the exact behaviour the decision changes,
and the project rule forbids editing existing tests. Both xfails stay in place. Ready
patches are described under "Blocked".

### Done

- **F2: a row times a matrix works in the pipeline.**
  - **Deviation from the decision.** The decision was to canonicalise at construction:
    `vᵀM → Transpose(MVMul(Mᵀ, v))`. That breaks three pinned tests:
    - `tests/test_geometry.py::TestTransposeOperators::test_mul_matrix`, a pre-phase-3 test,
      asserts `isinstance(Transpose(x) * M, MVMul)`.
    - `tests/test_operators.py::test_row_times_matrix_is_a_row` asserts the same.
    - The F2 xfail `test_f2_expand_row_times_matrix_sum` pins the row-form string
      `(a)'Kb+(a)'Jb`.

    The decided fallback (reject `v.T() * M`) breaks the same tests and also
    `test_error_messages.py:127,137`. Neither option is possible, so I kept the row form at
    construction and canonicalised where the engine needs a column.
  - `geomech/core/operations/geometry.py:221` adds `column_of(row)`. It builds the column
    structurally: `(vᵀ)ᵀ = v`, `(rᵀM)ᵀ = Mᵀr`, `(rᵀs)ᵀ = r s`, a sum of rows addend by
    addend, and d/dt, δ, ∫ through. For a `Transpose` node it returns the same object as
    `Transpose(row)`.
  - `geomech/core/math/extract.py:210,331,398,404,405,540,601` use `column_of(left)` in
    place of `Transpose(left)`.
  - `geomech/core/transformations/expand.py:78-107` (`_row_addends`, `_expand_row_mvmul`,
    `_expand_row_svmul`, plus the `VVMul` operands) distributes every row shape listed in
    F2 and keeps it in row form.
    - Exception: `d/dt((a+b)ᵀ) c` stays as it is, because `expand` does not distribute d/dt
      over a column sum either. Extraction now handles it.
  - `geomech/core/transformations/pull.py:16` hoists the scalar of a scaled row out of an
    inner product: `(aᵀs) b → aᵀb s`. This is restricted to `SVMul(row, s)`. A general
    `VVMul` hoist (including `(c s)ᵀ b`) changed 10 pinned systems, so it is not used.
  - The following now derive and pass `assert_eom`, and each equation matches the
    hand-derived Euler–Lagrange equation numerically:
    - `L = ½mẋ·ẋ − ½xᵀKx` (written `x.T()*K*x/2`);
    - `x.T()*A*x`;
    - `(x−e3).T()*K*(x−e3)`;
    - `dW = −(dx.T()*A*v)`.

    The tests are in `tests/test_row_forms.py` (14). The four systems are also
    `row_spring`, `row_quadratic_form`, `row_shifted_spring` and `row_damping_work` in the
    `tests/test_eom_oracle.py` SYSTEMS list.
  - Not supported: a row inside the kinetic energy (`v.T()*A*v`). This is pre-existing:
    `(Aᵀv)ᵀv` and `v.T()*(A*v)` also fail at HEAD, because IBP/collect handle only `Dot`.
  - README and `getting-started.md` gained the line `x.T() * K * x`. The MVMul error
    message is unchanged, and it is now accurate.
- **F3.** `geomech/utils/printing.py:562` adds `_signed_factor`. It is used in
  `_product_latex` (`:521-524`), in the right operand of a product, `Dot` and `Cross`, and
  for a transpose atom.
  - Outputs: `x \left(-a^T b\right)`, `A \left(-K^T\right)`, `c \left(-a^T\right)`,
    `-\left(-K^T\right)` and `-\left(-a^T b\right)`.
  - On 4000 random operator-built trees: 0 `--`, 0 `\times -`, 0 `\cdot -`.
- **F4.** `geomech/core/base/expressions.py:176`: `0 + x` is `x`.
  - `sum([a,b,c])` gives `(a+b+c)`, and `sum([x,y])` gives `(x+y)`.
  - `x + 0` is unchanged.
- **F5.**
  - `addition.py:47` quotes the first operand as written: `Add cannot combine SCALAR
    '(Dot(a,a)+Dot(b,a)x)' and VECTOR 'a'`.
  - `expressions.py:184`: `Sub cannot combine VECTOR 'a' and SCALAR '(1)'`.
  - Cosmetic: `errors.py:12` (`_short_repr`) truncates non-expression operands in
    `TypeError`s. `Dot(list(range(500)), a)` used to give a message of 1000 or more
    characters; it is now under 200. `_raise_constant_integral` truncates the name.
- **F7.** `quadrotor.ipynb` and `quadrotor_payload.ipynb` were re-executed and stored
  with `nbconvert --execute --ExecutePreprocessor.record_timing=False`, with
  `language_info.version` kept at `3.12.3`.
  - The changes are the LaTeX of 3.5 (14 and 16 lines).
  - A second execution is byte-identical.
  - 0 error outputs.
- **Cosmetic.**
  - `simplify.py:294` (`_cancelled`): a float coefficient sum at most `1e-12 ×` the
    largest summand is 0. `full_simplify(x*0.1*3 − x*0.3)` → `0`, and `x*1e-20 + y` is
    kept.
  - `simplify.py:207`: `((a)'x)'` → `ax`.
  - The topdown labels are `S·V`, `S·M`, `M·V`, `M·M`, `V·V`.
  - None of these changes a pinned string.
- **CHANGELOG.** The `s * Transpose(v)` sentence is corrected. It now says that the scalar
  on the left already built the correct row `(v)'s` and now builds `(vs)'`. New entries
  cover F2–F5 and the cosmetic items.

### Blocked (decision conflicts with a pinned phase-3 test)

- **F1 (restore `x * a.T()` → `SVMul(Transpose(a), x)`).**
  - The blocker is `tests/test_operators.py:124`
    `test_scalar_times_row_matches_row_times_scalar`. It asserts
    `str(s*aᵀ) == str(aᵀ*s)`, while `:120` pins `str(aᵀ*s) == "(as)'"`. Together with the
    xfail `str(x*a.T()) == "(a)'x"`, these cannot all hold.
  - The fix is ready as a one-line change: `expressions.py:120`, `_scale(b, a)` →
    `SVMul(b, a)`. The patch (`f1_code.patch`) was kept in a temporary directory; the `_mul`
    docstring needs the matching sentence.
  - I tried it in place:
    - Both F1 xfails pass.
    - The 61-system dump is byte-identical.
    - The fuzz `str` dump has 0 differences.
    - The only failure is the test at `:124`.
  - The restricted `pull` hoist above makes `test_f1_row_forms_cancel` pass once the tree
    is restored.
- **F6 (no address in `print_tree`, `#k` labels).**
  - The blocker is `tests/test_repr_debug.py:137` `test_tangent_vector_shows_its_parent`.
    It asserts `tree_str(om) == f"TS2(..., of=S2 'q' @{id(q):#x})"` and
    `of=SO3 'R' @{id(R):#x}`.
  - A ready implementation (`f6_code.patch`) was kept in a temporary directory. It adds
    `_point_labels`, which numbers distinct points with the same name, covering both
    `S2`/`SO3` leaves and tangent parents, in order of first appearance.
  - Applied, the results were:
    - `test_f6_tree_is_deterministic` passes.
    - `test_two_parents_with_the_same_name_are_told_apart` still passes
      (`of=S2 'q' #1` / `S2('q' #2, …)`).
    - The only failure is the test at `:137`.
  - Still to do after that: regenerate `rigid_pendulum.ipynb`/`spherical_pendulum.ipynb`,
    which makes `test_f6_notebooks_have_no_addresses` pass, and update the CHANGELOG
    `print_tree` entry.
  - Currently, re-executing these two notebooks changes only the `@0x…` addresses (3 and
    1 lines). After masking the addresses they are byte-identical to the stored output,
    so they were left as stored.

### Checks

- **`str` byte-identity against `1e10349`.** I ran the 61 argument-free systems of
  `vbl_systems`/`phase1_systems`/`test_eom_oracle` through `L`, `dW`, `δL+δW`,
  `apply_manifold_rules`, `full_simplify`, `compute_eom`, `to_standard_form` and
  `linearize`. Output: **235 396 bytes, `cmp` equal**. The 4 new oracle systems are
  excluded because they do not build at base.
  - Fuzz `str` dump: 4500 `engine_fuzz.Gen` trees (seeds 31–33) × 6 transforms against
    HEAD, **0 differences**.
  - Both worktrees were removed and `git worktree prune` was run.
- **Suite with `-W error`:** 1927 passed, 9 xfailed. The xfails are 5 pre-existing plus
  F1 ×2 and F6 ×2.
  - Flipped: F2 ×2, F3 ×2, F4, F5 ×2, F7, with no assertion changed.
  - New tests: `test_row_forms.py` (14), `test_phase3_fixes.py` (19), and +4 oracle
    systems.
- **Lint:**
  - `ruff format --check` and `ruff check` on `geomech/ tests/`: clean.
  - `ruff check examples/`: the same pre-existing errors as base. Example scripts are
    exempt.
- **Docs:** `mkdocs build --strict` exits 0.
- **Fuzz:**
  - `engine_fuzz.py --quick`: the summary table is identical to HEAD (36 distinct
    findings).
  - `--cases 3000 --seed 31`: identical to HEAD (189 findings, 0 wrong-value).
- **Examples:** all 6 exit 0 under `-W error`.
- **Notebooks:**
  - All four execute with 0 error outputs.
  - The quadrotor notebooks are byte-identical on a second run. The pendulum notebooks
    differ only in addresses (F6 blocked).
  - No stored output contains `/home/`, `/tmp/` or a user name.

### Arbitration of the blocked items (orchestrator)

F1 and F6 were blocked only by two tests the phase 3 implementer wrote this week
(`tests/test_operators.py` scalar-times-row and `tests/test_repr_debug.py` tangent-parent
address); both pinned exactly the behaviour the review found wrong. The fixer's patches
were applied, those two tests now assert the corrected behaviour (the pre-phase-3 tree
for `s * vᵀ`, a deterministic `of=SO3 'R'` with no address), the four F1/F6 xfail markers
were removed with assertions unchanged, and the two pendulum notebooks were regenerated
with execution timestamps stripped. Suite: 1931 passed, 5 xfailed (pre-existing).
