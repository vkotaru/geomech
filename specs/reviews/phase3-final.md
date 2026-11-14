# Phase 3 final review: operators, printing

Branch `fix/phase3-operators-printing`, range `1e10349..743a29b` (implementation `475c596`,
review `d4349ac`, fixes `743a29b`). Read-only on `geomech/` and `tests/`; nothing committed.

## Verdict: GO WITH CHANGES

Every check in the brief passes. No pinned equation changed: `str` output of 61 systems
through the whole pipeline is byte-identical to `1e10349`. The fuzz has no wrong-value
finding, and the four new row-form systems agree with a hand derivation. Three
documentation errors remain, and all three are things a user would copy:

1. The CHANGELOG describes the `s * vᵀ` tree that the arbitration reverted.
2. The README and getting-started operator line shows `x + 1` with `x` a vector, which
   raises.
3. Nothing a user reads says that a row in a kinetic energy is unsupported. The error
   message does not point to `Dot(v, A * v)` either.

All three are text changes. No code change is required to merge.

## Required changes

### R1. The CHANGELOG still describes the reverted F1 tree (`CHANGELOG.md:216-219`)

**Current text:**

> The other order, `s * Transpose(v)`, was not affected: it built the correct row
> `SVMul(Transpose(v), s)` (`(v)'s`). It now builds the same `Transpose(SVMul(v, s))`
> (`(vs)'`) as `Transpose(v) * s`, so its `str` and `==` changed; no pipeline builds it.

**What HEAD does:**
- `str(x * a.T())` is `(a)'x`. That is the old tree, restored by the arbitration and
  pinned by `test_f1_scalar_times_row_keeps_old_tree`.
- `str(a.T() * x)` is `(ax)'`.

**Replace the text with:**

> The other order, `s * Transpose(v)`, was not affected and is unchanged: it builds the
> row `SVMul(Transpose(v), s)` (`(v)'s`). The two orders therefore build different trees
> (`(v)'s` and `(vs)'`); `simplify` rewrites `(s vᵀ) w` to `s (vᵀ w)`, but not yet
> `(vᵀ s) w`.

The last clause is accurate:
- `full_simplify((s*a.T())*b - s*(a.T()*b))` gives `0`.
- `full_simplify((a.T()*s)*b - s*(a.T()*b))` gives `((as)'b+s(a)'b(-1))`.

### R2. `x + 1` with a vector `x` raises (`README.md:200`, `docs/getting-started.md:91`)

**The problem:** in that code block `x` is a vector (`x + y  # vector addition`, `M * x`).
`x + 1` raises:

```
ExpressionMismatchError: VAdd cannot combine VECTOR 'x' and SCALAR '(1)'
```

**Replace the line in both files with:**

```python
2 * x, a + 1, 1 - a, -x, x / 2   # numbers on either side (a number adds only to a scalar); / divides by a number only
```

If 99 characters is a concern, split the comment over two lines.

**Optionally**, `CHANGELOG.md:10` could say "`x + 1`, `1 + x`, `x - 1`, `1 - x` for a
scalar `x`".

### R3. A row in a kinetic energy is undocumented where users hit it (`README.md:199`, `docs/getting-started.md:90`)

**Where it is documented now:** only `CHANGELOG.md:28-29` mentions the limitation.

**What a user sees:** the operator line advertises `x.T() * K * x` with no qualification.
A user who writes the kinetic energy the same way, `v.T() * A * v / 2` (or
`v.T() * (A * v) / 2`), gets:

```
AlgebraicError: extract_linear_coeff: '\delta{x}' appears only inside the calculus operator
'\frac{d}{dt}(\delta{x})' (a higher order derivative, a variation or an integral of the
query); such a term has no linear coefficient
```

That message does not name the cause, and it does not suggest `Dot(v, A * v)`. With
`Dot(v, A * v) / 2`, the same system derives correctly: `check_eom` = 6.7e-9.

**Add this line directly after the `x.T() * K * x` line in both files:**

```python
                 # (in a potential or in dW; in a kinetic energy write Dot(v, K * v):
                 #  integration by parts handles Dot only)
```

The actionable error message belongs to item 4.2 of the plan (re-raise engine errors with
pipeline context); see follow-up F1. It does not block this merge, but 4.2 must cover
this case.

## Verified (commands and counts)

All commands use `source venv/bin/activate`. The venv is Python 3.12.13. `nbconvert`
7.17.2 was already present.

### 1. Suite, lint, docs, fuzz, examples, notebooks

| Check | Command | Result |
|---|---|---|
| Suite | `python -m pytest tests/ -q -W error -p no:cacheprovider` | **1931 passed, 5 xfailed** (the 5 are pre-existing), exit 0 |
| Lint | `ruff format --check geomech/ tests/ && ruff check geomech/ tests/` | 83 files formatted; all checks passed |
| Docs | `mkdocs build --strict -d <tmp>/site` | exit 0. Only the two existing INFO lines remain: `manifold-implementation-notes.md` is not in the nav, and the `reference/` link |
| Fuzz | `python -W error specs/reviews/scripts/engine_fuzz.py --cases 2000 --seed 5150 --no-json` | exit 0; **143 distinct findings**: 63 flag-dependent, 75 unsupported-linear-term, 5 documented-error-on-legal-input (all in `collect`, same counts as base); **0 wrong-value** |
| Examples | each of the 6 scripts with `python -W error examples/<f>.py` | all exit 0; no `/home` paths or `0x` addresses in the output |
| Notebooks | each of the 4 notebooks executed twice with `jupyter nbconvert --execute --ExecutePreprocessor.record_timing=False --output-dir <tmp>/runN` | details below |

**Fuzz at base**, same seed, script from the `1e10349` worktree: 158 findings. The extra
findings are:
- 9 documented-error-on-legal-input, including the `ExpressionMismatchError` rows in
  `expand` (1), `simplify` (1) and `full_simplify` (2/3) that phase 3 removes;
- 86 unsupported-linear-term;
- 0 wrong-value.

**Notebooks:**
- Run 1 and run 2 are byte-identical (`cmp`) for all four notebooks.
- Every cell's outputs equal the committed outputs (0 differing cells: 13/16/16/14 cells).
- There are 0 error outputs.
- The committed files contain no `0x…` address, no `/home`, no `/tmp` and no user name.
- There is no per-cell `execution` timing metadata.

### 2. `str` invariance against `1e10349`

**Method:**
- I added the worktree `git worktree add --detach <tmp>/wt_base 1e10349`. It is
  removed and pruned; `git worktree list` now shows only the main tree.
- I used my own dump script, which imports `geomech` and `tests.*` from the given root and
  asserts the import path. It ran every argument-free system builder in
  `tests/vbl_systems.py`, `tests/phase1_systems.py` and `tests/test_eom_oracle.py`.
- For each system it recorded:
  - `str(L)`, `str(dW)` and `str(δL+δW)`;
  - `apply_manifold_rules` and then `full_simplify`;
  - `compute_eom` (each key, variation vector and equation);
  - `to_standard_form(..., inputs)` and `linearize(..., inputs)`, with every field
    recursively passed through `str`.

**Result:**
- Base: 61 systems. HEAD: 65 systems. The only systems HEAD adds are `row_damping_work`,
  `row_quadratic_form`, `row_shifted_spring` and `row_spring`.
- For the 61 common systems the dumps are **199 909 bytes each and `cmp`-identical**.
- The same 3 pre-existing typed errors occur in both trees: a linearize
  `AlgebraicError` (term without a perturbation symbol) and 2 `collect` "target inside the
  matrix factor" errors.

### 3. Row-form systems

**Hand derivation (my own, not one of the four pinned systems).** K is a constant
non-symmetric matrix.
- L = ½ m ẋ·ẋ − ½ xᵀKx and δW = δx·f.
- ∂/∂x(½ xᵀKx) = ½(K + Kᵀ)x.
- The Euler–Lagrange equation in `compute_eom`'s sign convention is
  E = −m ẍ − ½(K + Kᵀ)x + f.

**Derived by HEAD** from `L = m*Dot(v,v)/2 - x.T()*K*x/2`:

```
(\frac{d}{dt}(\frac{d}{dt}(x))(-1)m+Kx(-0.5)+(K)'x(-0.5)+f)
latex: \left(-m \ddot{x} -\frac{1}{2} K x -\frac{1}{2} K^T x + f\right)
```

- The difference from the hand formula at a random point is 2.2e-16.
- `check_eom` = 6.9e-9.
- `to_standard_form`: `M = -m I`, `f = (Kx(-0.5)+(K)'x(-0.5))`, `G = I`.

**The four committed systems at HEAD:**
- `row_spring` (K symmetric): `…+Kx(-1.0)+f`
- `row_quadratic_form`: `…+Ax(-1)+(A)'x(-1)+f`, which is (A + Aᵀ)x
- `row_shifted_spring`: `…+Kx(-2)+Ke3(2)+f`
- `row_damping_work`: `…+A\frac{d}{dt}(x)(-1)+f`

All four match the docstrings of `tests/test_eom_oracle.py:186-216`.

`tests/test_row_forms.py::test_row_times_matrix_derives` compares each equation
numerically with the hand-derived formula (atol 1e-12) and with `assert_eom`. **No `str`
is pinned for these four systems**; see advisory A2.

**The kinetic-energy limitation** is neither documented nor actionable. See R3 and F1.

### 4. One example per row of the phase 3 table (output from HEAD)

```
== 3.1 Transpose(...) * x
  R.T() * v    -> MVMul  str=(R)'v  latex=R^T v
  R.T() * J    -> MMMul  str=(R)'J  latex=R^T J
  R.T() * s    -> SMMul  str=(R)'s  latex=s R^T
  a.T() * s    -> Transpose str=(as)'  latex=s a^T
  s * a.T()    -> SVMul  str=(a)'s  latex=s a^T
== 3.2 double transpose collapses
  Transpose(Transpose(a)) is a: True | ((aᵀJ)ᵀ)ᵀ: (a)'J
  simplify(((aᵀJ)ᵀ)ᵀ J): (a)'JJ | expand((a+(bᵀJ)ᵀ) cᵀ): (a(c)'+((b)'J)'(c)')
== 3.3 numbers on both sides
  2*s -> (2)s   s+1 -> (s+(1))   1+s -> ((1)+s)   s-1 -> (s+(-1))   1-s -> ((1)+s(-1))
  -s -> s(-1)   -a -> a(-1)      a/2 -> a(0.5)    sum([a, b, c]) -> (a+b+c)
== 3.4 error messages
  Dot(m, x)      ExpressionMismatchError: Dot cannot combine SCALAR 'm' and VECTOR 'x'
  m + x          ExpressionMismatchError: Add cannot combine SCALAR 'm' and VECTOR 'x'
  x - m          ExpressionMismatchError: Sub cannot combine VECTOR 'x' and SCALAR 'm'
  Hat(m)         ExpressionMismatchError: Hat expects a VECTOR, got SCALAR 'm'
  Vee(x)         ExpressionMismatchError: Vee expects a MATRIX, got VECTOR 'x'
  J * x.T()      ExpressionMismatchError: MVMul: a matrix times a row vector is not defined: 'J' * '(x)''; a matrix multiplies a column vector (M * v), and a row vector multiplies a matrix from the left (v.T() * M)
  x * J          ExpressionMismatchError: MVMul: a column vector times a matrix is not defined: 'x' * 'J'; write M * v for the column M v, or v.T() * M for the row vᵀ M
  x * a          ExpressionMismatchError: VVMul: two column vectors cannot be multiplied: 'x' * 'a'; write Dot(a, b) or a.T() * b for the inner product, a * b.T() for the outer product, or Cross(a, b)
  x.T() * a.T()  ExpressionMismatchError: VVMul: two row vectors cannot be multiplied: '(x)'' * '(a)''; write a.T() * b for the inner product aᵀb or a * b.T() for the outer product a bᵀ
  m * 'x'        TypeError: cannot apply '*' to Scalar 'm' and str 'x'; operands must be geomech expressions or numbers (int, float)
  x / m          TypeError: cannot divide Vector 'x' by the expression Scalar 'm'; only division by a number is supported. For a symbolic divisor declare its inverse as a symbol, e.g. Scalar('1/m', attr=['Constant'])
  a - 1          ExpressionMismatchError: Sub cannot combine VECTOR 'a' and SCALAR '(1)'
  Scalar('k', attr=['Constant']).t_integrate()  AlgebraicError: t_integrate: cannot integrate the constant scalar 'k' in time: ∫ k dt = k·t, and the tree has no time symbol. Only terms that are time derivatives (d/dt(x)) integrate to an expression
== 3.5 LaTeX of the dynamics audit F8 pairs
  Cross(a,Cross(b,c))   a \times \left(b \times c\right)   | Cross(Cross(a,b),c)  \left(a \times b\right) \times c
  (R)'Cross(a,b)        R^T \left(a \times b\right)        | Cross((R)'a,b)       \left(R^T a\right) \times b
  RCross(a,b)           R \left(a \times b\right)          | Cross(Ra,b)          \left(R a\right) \times b
  (Cross(a,b))'         \left(a \times b\right)^T          | (Ra)'                \left(R a\right)^T
  aDot(b,c)             \left(b \cdot c\right) a           | Cross(a,b(-1))       a \times \left(-b\right)
  ss(-0.5)              -\frac{1}{2} s s                   | s(2.5)               2.5 s
== 3.6 repr / print_tree
  repr: m·g*·<v:e3*, v:x> | cross(v:\omega_{q}, S2:q) | (SO3:R·v:x)' | M:J M:K*
Dot
  ├─ L: TSO3('\eta_{R}', of=SO3 'R')
  ├─ R: MVMul
    ├─ L: Matrix('I3', is_constant, is_identity)
    ├─ R: Cross
      ├─ L: TSO3('\Omega_{R}', of=SO3 'R')
      ├─ R: MVMul
        ├─ L: Transpose
          ├─ SO3('R', is_orthogonal, is_manifold)
        ├─ R: Vector('e3', is_constant)

expr: <v:\eta_{R}, M:I3*·cross(v:\Omega_{R}, SO3:R'·v:e3*)>
(two distinct S2('q') objects, tree_str(..., style="indent")):
Cross
  ├─ L: TS2('\omega_{q}', of=S2 'q' #1)
  ├─ R: S2('q' #2, is_unit_norm, is_manifold)
```

**Status of each row:**
- **3.1 to 3.6 are all done.** The two deviations from the plan text are justified:
  - 3.1: `vᵀ * s` is `Transpose(SVMul(v, s))`, not `SVMul(Transpose(v), s)`. That is one
    transpose either way, and pre-phase-3 tests pin the tree.
  - 3.5: a nested `MVMul`/`MMMul` is not wrapped, because the product is associative.
- **Debris removed:** no bare `raise ExpressionMismatchError`, no empty
  `NotImplementedError` except the two abstract ones that tests expect, and
  `UndefinedCaseError` is gone.
- **Real double rigid pendulum output** (`examples/rigid_pendulum.py`): every nested cross
  product is parenthesised, for example
  `m2 \left(\left(\Omega_{R1} \times l1\right) \times \left(R1^T R2 \left(\Omega_{R2} \times \rho2\right)\right)\right)`.

### 5. Operator and printing sections of README and getting-started

Both files carry identical blocks. I ran every line with `a, b` scalars, `x, y` vectors,
`M, K` matrices and `R` an `SO3`:
- Every snippet behaves as its comment says, except `x + 1` (R2).
- `to_latex`, `print_eom` and `print_tree` are exported from `geomech`.
- The documented `repr` example `m*·g*·<v:e3*, v:x>` is reproduced exactly when `m` and
  `g` are constants.
- `print_tree` defaults to indent style, and `style="topdown"` gives the sketch with the
  `m*`/`v:e3*` leaves.

## Advisory (no merge block)

- **A1.** `row_spring` gives `Kx(-1.0)`, a float from 0.5 + 0.5 when `Kᵀx` merges into
  `Kx`. `to_latex` prints `-K x`, so only `str` shows it.
- **A2.** The four row systems are checked only numerically. The contributor notes ask pipeline tests
  to assert the derived equation itself. The numeric comparison with the hand-derived
  formula meets the intent, but pinning `str` for `row_spring`/`row_quadratic_form`
  (above) would catch a regression to an equivalent but uglier form.
- **A3.** `language_info.version` differs between notebooks: `3.12.3` for the quadrotor
  notebooks, `3.12.13` for the pendulum notebooks, which the arbitration step regenerated.
  This is harmless but churns on the next run.
- **A4.** Notebook display cells print sums with their outer parentheses
  (`x_L = \left(x -l q\right)`). `display_latex` (`printing.py:736`) does not apply
  `_strip_outer_parens`, while `display_standard_form` now does.
- **A5.** `__rtruediv__` (`geomech/core/base/expressions.py:221`) quotes `{other!r}`
  without `_short_repr`, unlike every other operator message.
- **A6.** `_CANCEL_RTOL = 1e-12` (`simplify.py:291`) drops float like-terms that cancel to
  roundoff. This is sound, but it is a convention to list in the 4.7 "Conventions" page.
- **A7.** The plan's "regenerate the paper's listings from `to_latex` and diff by eye" is
  moot. The accompanying paper's sources contain no `to_latex` output (the `lstlisting` blocks
  are Python input; the equations are hand-typed). Neither phase report says so.
- **A8.** `print_tree` uses `├─` for every child, including the last one. This is
  cosmetic.

## Follow-ups for phases 4–5

- **F1 (4.2, must cover).** A row in the kinetic energy fails with the non-actionable
  message at `geomech/core/math/extract.py:73-77`.
  - **Cause:** IBP rewrites only `Dot(ẋ, f)` (`geomech/core/math/ibp.py:132`), so a
    `VVMul` inner product holding `d/dt(δx)` survives.
  - **4.2 should raise:** "a d/dt(δx) survived integration by parts in term '…': write
    inner products in the kinetic energy as Dot(v, A * v), not v.T() * A * v".
  - **Better (5.x):** have `collect` (`geomech/core/math/collect.py:49`) rewrite an inner
    product `VVMul(row, col)` to `Dot(column_of(row), col)` before IBP. That lifts the
    limitation.
- **F2 (4.4 `check_tree`).** Rows are still accepted where only columns make sense:
  - `Dot(a.T(), b)`, `geometry.py:27`;
  - `Cross(a.T(), b)`, `:53`;
  - `Hat(a.T())`, `:79`;
  - a `VAdd` that mixes a row and a column (`a.T() + b`, `addition.py:103`), whose
    orientation `is_row_vector` takes from the first addend (`geometry.py:213`).

  All of these are pre-existing; the invariant checker should flag them.
- **F3 (5.1).** `pull` hoists the scalar of a scaled row only for `SVMul(row, s)` on the
  left of a `VVMul`, and it does not recurse into a `VVMul` otherwise
  (`geomech/core/transformations/pull.py:16`).
  - **Consequence:** `(vᵀ s) w` and `s (vᵀ w)` do not merge (see R1).
  - **Fix:** a general hoist, once the 10 pinned systems it touches are re-derived.
- **F4 (5.4).** Vector like terms with a symbolic coefficient do not cancel. This is
  pre-existing and unchanged:
  - `full_simplify(a*x - a*x)` gives `(xa+xa(-1))`;
  - `full_simplify(a*(x+y) - a*x)` gives `(xa+ya+xa(-1))`.

  The location is `geomech/core/transformations/simplify.py:431` (`_collect_like_terms`).
- **F5 (5.5).** `_wrap_numeric` (`geomech/core/operations/multiplication.py:21`) now
  duplicates `_node_operands(numeric=True)`. `_is_row`/`_q` (`:30`, `:36`) are one-line
  shims around `is_row_vector`/`_operand_str`. Fold them into the constants and helpers
  module.
- **F6 (5.6).** `to_latex` (`geomech/utils/printing.py:614`) is now about 120 lines of
  `isinstance` chains plus about 10 helpers. Split it by node class along with
  `_eliminate`/`pull`.
- **F7 (4.7 docs).** The "Conventions" page should state:
  - the two row forms (`(v)'s` from `s * vᵀ` and `(vs)'` from `vᵀ * s`);
  - that `0 + x` is `x` but `x + 0` keeps `(0)` until `simplify`;
  - the float-cancellation tolerance (A6).
