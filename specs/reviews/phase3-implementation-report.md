# Phase 3 implementation report: crashes on legal input, operators, printing

Branch `fix/phase3-operators-printing`, created from `1e10349` (`fix/phase2-validation`).
Nothing is committed. Plan: `specs/fix-plan.md`, phase 3, items 3.1 to 3.6.

## Summary

| Check | Result |
|---|---|
| Suite | **1847 passed, 5 xfailed**. Base was 1699 passed, 5 xfailed. The 148 new tests are in four new files. No existing test was edited. |
| `str` byte-identity | **Identical.** I ran all 61 argument-free systems of `tests/vbl_systems.py`, `tests/phase1_systems.py` and `tests/test_eom_oracle.py` through `str(L)`, `str(dW)`, `str(δL+δW)`, `apply_manifold_rules`, `full_simplify`, `compute_eom`, `to_standard_form` (fields via `str`) and `linearize` (every attribute via `str`), in a worktree at `1e10349` and in this tree. Output: 235 396 bytes each, `cmp` equal. The narrower dump without the stages is 166 636 bytes and also equal. The same 3 pre-existing typed errors occur in both trees. The worktree is removed (`git worktree prune`). |
| `ruff format --check` / `ruff check` on `geomech/ tests/` | clean (80 files) |
| `ruff check examples/` | 111 errors (101 F405, 5 F403, 3 I001, 2 F841), all pre-existing. The count is identical at base. the contributor notes exempts the example scripts. |
| `mkdocs build --strict` | exit 0. Only the existing INFO lines remain (`manifold-implementation-notes.md` not in nav; `reference/` link). |
| `engine_fuzz.py --quick` | runs. expand/pull/vector_rules/simplify/full_simplify: 200/200 ok, 0 findings. |
| `engine_fuzz.py` (3000 cases, default seed), base vs new | "documented-error-on-legal-input" in `expand`/`simplify`/`full_simplify` (incl. adversarial) went from **9 findings, 17 `ExpressionMismatchError`s to 0**. 0 wrong-value in both. Base had 221 distinct findings, new has 202. The other rows match except the extract rows, which drew a slightly different random stream (cases 798→757 etc.). |
| Examples (6) | all exit 0. Output differs from base only in `print_eom` lines (`action_integral`, `point_mass`, `rigid_pendulum`, now `to_latex`) and in the error message of `sanity_check` (`Add cannot combine SCALAR 'a' and VECTOR 'x'`). |
| Notebooks (4) | all execute, with 0 error outputs. Differences are listed below. |
| LaTeX compile | `pdflatex -halt-on-error` on `eom_to_latex` for the double rigid pendulum and the quadrotor with load, including their standard-form `M`/`f`/`G` blocks: rc 0, PDF produced, no `^T^T`. `tests/test_latex_precedence.py::test_eom_document_compiles` repeats this; it is skipped when `pdflatex` is not installed. |

New test files and counts: `tests/test_operators.py` (46), `tests/test_error_messages.py` (34),
`tests/test_latex_precedence.py` (39), `tests/test_repr_debug.py` (29).

## 3.1 `Transpose(...) * x`

- `geomech/core/base/expressions.py:94` `_mul(a, b)` is now the only `*` dispatch, and it
  is reached from `Expr.__mul__`/`__rmul__` (`:182`, `:185`). `Transpose.__mul__` is
  deleted, so a matrix-typed transpose dispatches like any matrix: `Rᵀ * v → MVMul`,
  `Rᵀ * M → MMMul`, `Rᵀ * s → SMMul`.
- `expressions.py:106` `_scale`: a scalar times the row `Transpose(v)`, on either side,
  builds **`Transpose(SVMul(v, s))`**. That is one transpose of a column and no double
  transpose. `str` is `(vs)'`, the old one was `((v)'s)'`.
  - **Deviation from the plan text.** The plan says `SVMul(Transpose(v), s)`. But
    `tests/test_geometry.py:605-616` and `:624-628` (`TestTransposeOperators.test_mul_scalar`,
    `test_mul_numeric`) pin `isinstance(result, Transpose)` and
    `isinstance(result.expr, SVMul)`. I tried to edit those two tests and the edit was
    refused, because the instructions forbid editing existing tests. `Transpose(SVMul(v, s))`
    is the same row (`(v s)ᵀ = s vᵀ`) and passes both the old tests and the
    no-double-transpose requirement. `to_latex` prints it as `s v^T`.
- Row vectors that are not a `Transpose` node (`vᵀK`, `s·vᵀ`) multiply correctly:
  `VVMul` and `MVMul` decide orientation with `is_row_vector` (see 3.2).

## 3.2 `Transpose(Transpose(x))` collapses at construction

- `geomech/core/operations/geometry.py:133` `Transpose.__new__` returns `expr.expr` when
  `expr` is a `Transpose`. The `__init__` guard after it skips re-initialisation.
- Collapsing at construction changed **no** pinned `str` (byte-identity check above), so
  the `simplify` fallback was not needed.
- `geometry.py:192` `is_row_vector(expr)` gives the structural orientation: a row is
  `Transpose(column)`, `MVMul(row, M)`, `SVMul(row, s)`, a `VAdd` whose first addend is a
  row, or d/dt/δ/∫ of a row. It replaces the test "the operand is a `Transpose` node" in
  - `multiplication.py:152-177` (`MVMul`),
  - `multiplication.py:222-251` (`VVMul`),
  - the old `MatrixExpr.__mul__` check.
- Fuzz finding 6 and follow-up F-4 (`multiplication.py:148,207` at base): all four
  reproducers now succeed and evaluate correctly. They are tested in
  `tests/test_operators.py::TestRowOrientationSurvivesTransformations`:
  - `expand((a + (bᵀK)ᵀ) cᵀ)`
  - `expand(cᵀ(a + (bᵀK)ᵀ))`
  - `expand((a + bᵀᵀ) cᵀ)`
  - `simplify`/`full_simplify(((aᵀK)ᵀ)ᵀ Q)` → `(a)'KQ`
- `geomech/core/math/extract.py:210, 331, 398, 404, 405, 540, 601, 605`: `left.expr` (the
  column of a row) became `Transpose(left)`. With the collapse this returns the identical
  object for a `Transpose` node, which the byte-identity check confirms. A row that is not
  a `Transpose` node no longer gives an `AttributeError`; it gives a coefficient, or the
  existing typed `AlgebraicError`.
  - Known limitation: extracting `a` from `(aᵀK) b` raises
    `AlgebraicError: _vector_from_vector: unsupported node Transpose ...`. That shape could
    not be built before phase 3. The error is typed, not silent.
- `specs/reviews/scripts/engine_fuzz.py`:
  - The `UndefinedCaseError` import is removed (3.4).
  - `_col(v)` at line ~1196 now uses `is_row_vector`. The old "strip `Transpose` layers"
    turned the scaled row `Transpose(SVMul(Transpose(c), s))` into the row
    `SVMul(Transpose(c), s)`. Base accepted that as the left of an outer product and the
    new constructor rightly refuses it.

## 3.3 Numeric operands

All of these are in `geomech/core/base/expressions.py`:
- `Expr.__add__/__radd__/__sub__/__rsub__/__mul__/__rmul__/__neg__/__truediv__/__rtruediv__`
  (`:167-210`) apply to leaves and operation nodes alike.
- The duplicate implementations on `ScalarExpr`, `VectorExpr`, `MatrixExpr` and
  `_BaseMixin` are removed (`mixins.py` lost 79 lines).
- The trees are the old ones:
  - `a - b` → `Add(a, Mul(b, -1))` / `VAdd(a, SVMul(b, -1))` / `MAdd(a, SMMul(b, -1))`.
  - `-x` is `x * -1`, so `str(-x) == str(x * -1)`, for example `x(-1)`. A test checks this
    for scalar, vector, matrix, `Dot` and `Cross`.
- New trees:
  - `2*x → (2)x`; `2*v → v(2)`, the same as `v*2`.
  - `x+1 → (x+(1))`, `1+x → ((1)+x)`.
  - `x-1 → (x+(-1))` (`:173`, not `x + (1)(-1)`).
  - `1-x → ((1)+x(-1))`.
  - `x/2 → x(0.5)`. Division by an expression raises `TypeError`, and division by zero
    raises `ZeroDivisionError`.
- numpy (F-3):
  - `_is_number` (`:32`) accepts `np.integer`/`np.floating`.
  - `Number` (`:637`) converts them to the equal Python `int`/`float`, so
    `str(Number(np.int64(2))) == "(2)"`.
  - `Expr.__array_ufunc__ = None` (`:163`) makes `np.float64(2) * x` defer to `__rmul__`.
  - `bool` and `np.bool_` are rejected. The message keeps "got bool True", which
    `tests/test_phase2_adversarial.py:712` pins.
- The `Add`/`VAdd`/`MAdd`/`Mul`/`SVMul`/`SMMul` constructors wrap numbers:
  `_node_operands(..., numeric=True)` (`expressions.py:132`), called from `addition.py:48`
  and `multiplication.py:73,95,119`.
- `Expr.T()` (`expressions.py:249`) moved up from `VectorExpr`, so `x.T().T()` and `R.T()`
  work.

## 3.4 Error messages

`geomech/utils/errors.py:16-55`: `ExpressionMismatchError(operation, ltype, rtype, lhs=None, rhs=None, *, msg=None)`
and `ExpressionMismatchError.unary(op, expected, expr)`. Without operands it keeps the old
"cannot be performed between" text. Operands longer than 80 characters are shortened with
`...`. Class names and bases are unchanged.

Exact new messages:

| Input | Message |
|---|---|
| `Dot(m, x)` (A3) | `ExpressionMismatchError: Dot cannot combine SCALAR 'm' and VECTOR 'x'` |
| `m + x` (A4) | `ExpressionMismatchError: Add cannot combine SCALAR 'm' and VECTOR 'x'` |
| `x * J` (A5) | `ExpressionMismatchError: MVMul: a column vector times a matrix is not defined: 'x' * 'J'; write M * v for the column M v, or v.T() * M for the row vᵀ M` |
| `J * x.T()` (A6, was bare) | `ExpressionMismatchError: MVMul: a matrix times a row vector is not defined: 'J' * '(x)''; a matrix multiplies a column vector (M * v), and a row vector multiplies a matrix from the left (v.T() * M)` |
| `x.T() * x.T()` (A7) | `ExpressionMismatchError: VVMul: two row vectors cannot be multiplied: '(x)'' * '(x)''; write a.T() * b for the inner product aᵀb or a * b.T() for the outer product a bᵀ` |
| `x * y` | `ExpressionMismatchError: VVMul: two column vectors cannot be multiplied: 'x' * 'y'; write Dot(a, b) or a.T() * b for the inner product, a * b.T() for the outer product, or Cross(a, b)` |
| `x - m` (A8) | `ExpressionMismatchError: Sub cannot combine VECTOR 'x' and SCALAR 'm'` |
| `Hat(m)` (A9) | `ExpressionMismatchError: Hat expects a VECTOR, got SCALAR 'm'` |
| `Vee(x)` | `ExpressionMismatchError: Vee expects a MATRIX, got VECTOR 'x'` |
| `VAdd(m, x)` | `ExpressionMismatchError: VAdd expects VECTOR operands, got SCALAR 'm'` |
| `J + x` | `ExpressionMismatchError: MAdd cannot combine MATRIX 'J' and VECTOR 'x'` |
| `Cross`/`Mul`/`SVMul`/`SMMul`/`MVMul`/`MMMul`/`VVMul` type errors | `<Node> cannot combine <TYPE> '<lhs>' and <TYPE> '<rhs>'` |
| `Scalar('m', attr=['Constant']).t_integrate()` (A16) | `AlgebraicError: t_integrate: cannot integrate the constant scalar 'm' in time: ∫ m dt = m·t, and the tree has no time symbol. Only terms that are time derivatives (d/dt(x)) integrate to an expression` (`vector`/`matrix` for the other two leaves) |
| `m * "x"` (A29) | `TypeError: cannot apply '*' to Scalar 'm' and str 'x'; operands must be geomech expressions or numbers (int, float)` |
| `x * True` | `TypeError: cannot apply '*' to Vector 'x': got bool True; a bool is not a number here, write 1 or 0` |
| `Add(m, None)` | `TypeError: Add: operands must be geomech expressions or numbers, got NoneType None` |
| `Dot(x, 2)` | `TypeError: Dot: operands must be geomech expressions, got int 2` |
| `x / m` (A33) | `TypeError: cannot divide Vector 'x' by the expression Scalar 'm'; only division by a number is supported. For a symbolic divisor declare its inverse as a symbol, e.g. Scalar('1/m', attr=['Constant'])` |
| `2 / x` | `TypeError: cannot divide 2 by the expression Vector 'x'; only division by a number ...` |
| `x / 0` | `ZeroDivisionError: cannot divide Vector 'x' by zero` |

Where each message is raised:
- `geometry.py:36,62,84,111` (`Dot`, `Cross`, `Hat`, `Vee`)
- `multiplication.py:77,101,125,155,169,177,200,235,243,251`
- `addition.py:48-55`
- `expressions.py:44` (`_as_operand`), `:79` (`_sub`), `:192-210` (division), `:301`
  (`_raise_constant_integral`, called from `:360,436,549`).

Removed items:
- `UndefinedCaseError` is deleted from `errors.py`. A grep showed its only raises were the
  `case _` arms of the old type dispatch, which is unreachable because every node type is
  SCALAR/VECTOR/MATRIX. It was not exported from `geomech`; the only outside user was
  `engine_fuzz.py`, which is updated.
- No bare `raise ExpressionMismatchError` or empty `NotImplementedError` is left. A test
  (`test_error_messages.py::TestErrorClasses::test_no_bare_raise_left`) greps the package.
  The only hits it allows are the abstract `Expr.__str__` and `get_tangent_vector`, which
  `tests/test_manifolds.py` expects to raise `NotImplementedError`.

## 3.5 Precedence-aware `to_latex`

`geomech/utils/printing.py:445-682`:
- helpers `_wrap`, `_numeric_value`, `_format_number`, `_coefficient`, `_scalar_factors`,
  `_product_latex`, `_needs_parens_in_cross`, `_needs_parens_in_product`,
  `_is_transpose_atom`;
- `to_latex` at `:564`.

Rules:
- `Cross` operands that are `Cross`/`MVMul`/`SVMul`/`VVMul` are wrapped. A `VAdd` is
  always wrapped.
- `MVMul`/`MMMul`/`VVMul` operands that are `Cross`/`SVMul`/`SMMul` are wrapped.
  - **Deviation:** a nested `MVMul`/`MMMul` is *not* wrapped, because the matrix product is
    associative: `J R a` reads correctly either way, and wrapping would add a pair to every
    `R1^T R2 (…)`. The F8 cases involving `MVMul` all involve a `Cross`, and they all print
    differently now.
- `Transpose` puts `^T` directly only on an atom: a leaf without `^`, `\dot`/`\ddot` of a
  leaf, or `\widehat{}`. Any other operand becomes `\left(…\right)^T`, so `^T^T` is
  impossible. `Transpose(SVMul(v, s))` prints `s v^T`.
- A `Dot` factor is wrapped whenever anything else multiplies it: `\left(b \cdot c\right) a`.
- The numeric factors of a `Mul` chain, or of a chain of `SVMul`/`SMMul` scalars, fold into
  one coefficient printed first. `0.5 → \frac{1}{2}` and `2.0 → 2`.
- A scalar multiple of a `Cross` is wrapped (`m \left(a \times b\right)`), unless the
  multiplier is only a sign.
- `TimeDerivative`/`Variation` of a composite prints as `\frac{d}{dt}\left(…\right)` or
  `\delta\left(…\right)`.

`print_eom` (`:25`), `eom_to_latex` (`:787`) and `render_eom` (through `eom_to_latex`) now
use `_eom_line` (`:19`): `\int <to_latex(ξ)> \cdot \left(<to_latex(E)>\right) \, dt = 0`.
The old output was `\int{ξ \cdot \Big(str(E)\Big)}dt=0`.

Two related fixes:
- `_strip_outer_parens` (`:692`) now strips only a matching outer pair. The old version
  would turn `\left(a\right) \left(b\right)` into unbalanced LaTeX.
- `display_standard_form` (`:767`) strips the outer parentheses of `f`, which is one more
  summand.

### LaTeX before/after (dynamics audit F8 table and extras)

| Expression | Before (`1e10349`) | After |
|---|---|---|
| `Cross(a, Cross(b, c))` | `a \times b \times c` | `a \times \left(b \times c\right)` |
| `Cross(Cross(a, b), c)` | `a \times b \times c` | `\left(a \times b\right) \times c` |
| `MVMul(Transpose(R), Cross(a, b))` | `R^T a \times b` | `R^T \left(a \times b\right)` |
| `Cross(MVMul(Transpose(R), a), b)` | `R^T a \times b` | `\left(R^T a\right) \times b` |
| `Cross(MVMul(R, a), b)` | `R a \times b` | `\left(R a\right) \times b` |
| `MVMul(R, Cross(a, b))` | `R a \times b` | `R \left(a \times b\right)` |
| `Transpose(Cross(a, b))` | `a \times b^T` | `\left(a \times b\right)^T` |
| `Transpose(MVMul(R, a))` | `R a^T` | `\left(R a\right)^T` |
| `SVMul(a, Dot(b, c))` | `b \cdot c a` | `\left(b \cdot c\right) a` |
| `Cross(a, SVMul(b, -1))` | `a \times -b` | `a \times \left(-b\right)` |
| `SVMul(Cross(a, b), -1)` | `-a \times b` | `-a \times b` |
| `Mul(Mul(s, s), -0.5)` | `-s s \frac{1}{2}` | `-\frac{1}{2} s s` |
| `Mul(s, 2.5)` | `s 2.5` | `2.5 s` |
| `Transpose(v) * s` | `s v^T^T` (does not compile) | `s v^T` |
| `Mul(Dot(a,b), Dot(b,c))` | `a \cdot b b \cdot c` | `\left(a \cdot b\right) \left(b \cdot c\right)` |
| double rigid pendulum, η_R1 row (excerpt) | `-m1 \Omega_{R1} \times \rho_1 \times \Omega_{R1} \times \rho_1` | `-m1 \left(\Omega_{R1} \times \left(\rho_1 \times \left(\Omega_{R1} \times \rho_1\right)\right)\right)` |
| same row | `… l_1 \times R1^T R2 \Omega_{R2} \times \rho_2` | `… \left(\Omega_{R1} \times l_1\right) \times \left(R1^T R2 \left(\Omega_{R2} \times \rho_2\right)\right)` |

`tests/test_latex_precedence.py` covers the following:
- `test_ambiguous_pairs_print_differently`: 8 pairs, every row of F8, where `str` differs
  but the old LaTeX was identical or misleading.
- `test_exact_text`: 23 exact strings.
- `print_eom`/`eom_to_latex` use `to_latex`, with no `Cross(`/`Dot(`/`Hat(` left.
- The real double-pendulum row.
- `pdflatex` compiles both documents.

**Existing tests that pin old LaTeX text: none.** I edited no existing test.
- `tests/test_printing.py::TestEomToLatex` asserts only `\documentclass`, `\begin{equation}`,
  `\end{document}` and `"a" in latex`. All of these still hold.
- `tests/test_vbl.py:669` only checks `isinstance(to_latex(block), str)`.
- `tests/test_eom.py` calls `print_eom` as a smoke test.

## 3.6 `repr` and `print_tree`

`geomech/utils/printing.py`:
- `_leaf_prefix` (`:87`): `v:` for Vector/TS2/TSO3, `M:` for Matrix/SkewSymmMatrix,
  `S2:`, `SO3:`; scalars get none. `_compact_leaf` reuses it.
- `_repr_leaf` (`:105`) adds a trailing `*` for constants. Numeric leaves print their name
  bare (`(2)`, `0.5`).
- `repr_str` (`:121`):
  - **Products are now joined with `·` instead of `*`.** Otherwise the constant marker
    collides with multiplication (`m**g*`). Example: `m*·g*·<v:e3*, v:x>`.
  - The leaf dot-accent becomes `v:ẋ`, with the prefix kept.
  - A transposed composite gets parentheses: `(SO3:R·v:x)'`.
- `print_tree` (`:40`) defaults to `style="indent"`.
  - `tree_str` keeps `style="topdown"` because `tests/test_printing.py::TestTreeStrTopdown`
    pins it (`"m*" in tree_str(m)`, `"v:x" in tree_str(x)`).
  - `_LEAF_FLAGS` (`:366`) adds `is_identity`, `is_skew_symmetric` and `is_orthogonal`.
  - `_leaf_detail` (`:380`) appends `of=S2 'q' @0x…` / `of=SO3 'R' @0x…` for `TS2`/`TSO3`.
    Two `S2("q")` therefore render differently (audit C6).

`str` is untouched. The one repr-dependent test (`tests/test_validation.py:606`) compares
`repr(out)` with the warning that also uses `repr`, so it still passes.

Docs: `README.md` and `docs/getting-started.md` "Operations" gained the new operator forms,
the error-message style and a **Printing** block (`to_latex`, `print_eom`, `repr` legend,
`print_tree`). No doc page showed repr or print_tree output.

Notebooks: `notebooks/rigid_pendulum.ipynb` and `notebooks/spherical_pendulum.ipynb` show
`print_tree(eqn)` output, which is now in the indent style with the new repr. I regenerated
their stored outputs with `nbconvert --execute --inplace --ExecutePreprocessor.record_timing=False`
and restored the `language_info.version` line. Their print_tree cells now contain memory
addresses (`@0x…`), so those lines change on every run. `quadrotor.ipynb` and
`quadrotor_payload.ipynb` only import `print_tree`, so their stored outputs are left as they
were.

### Notebook output differences (executed at `1e10349` vs executed now)

All four notebooks execute with exit 0 in both trees and produce no error outputs. Every
difference is a rendering change from 3.5/3.6; no equation changed.

| Notebook | Cells differing | Nature |
|---|---|---|
| `quadrotor` | 4 of 7 | Display of L: `\frac{1}{2} m \left(\dot{x} \cdot \dot{x}\right)` was `m \dot{x} \cdot \dot{x} \frac{1}{2}`. EOM: `-\Omega_{R} \times \left(J \Omega_{R}\right)` was `-\Omega_{R} \times J \Omega_{R}`. Standard form: `f` is no longer wrapped in an extra `+ \left(…\right)`. |
| `quadrotor_payload` | 5 of 8 | Same as above, plus `-2 m_L l \left(\omega_{q} \times \left(q \times \dot{x}\right)\right)`, which was `-m_L l 2 \omega_{q} \times q \times \dot{x}`. The old form reads as `(ω×q)×ẋ`, a different vector. |
| `rigid_pendulum` | 5 of 8 | LaTeX as above (`-m g \left(\rho \times \left(R^T e3\right)\right)`, was `-m g \rho \times R^T e3`). `print_tree` cell: indent style with flags and `of=SO3 'R' @0x…`, and the `expr:` line uses the new repr. |
| `spherical_pendulum` | 4 of 7 | LaTeX `-m g l \left(q \times e3\right)`. `print_tree` cell as above. |

## Notes for the reviewer

- **Unchanged behaviour.**
  - Every operator that existed before builds the same tree; `_mul`/`_sub` copy the old
    dispatch arm by arm.
  - A tree that the old `isinstance(x, Transpose)` test classified as a row but that is a
    column (`Transpose(row)`) is now refused as the left of an inner product. No pipeline
    path builds one: the suite and the 61-system stage dump are unchanged.
- **`x - 1` vs `x - y`.** `x - 1` builds `x + (-1)`, and `x - y` builds `x + y(-1)` as
  before. Under the existing `str` equality, `x - 1 != x + Number(1) * -1`.
- **Division.** `x / 3` gives the float `0.3333333333333333` as a coefficient. Exact
  rationals are not attempted.
- **Not done (out of scope).** `Dot` operands are not parenthesised (`a \cdot b \times c`
  is the standard triple-product notation). `VAdd` can still mix a row and a column; no
  constructor check was added for that.
