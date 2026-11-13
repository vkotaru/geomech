# Engine audit: symbolic core of `geomech` (branch `feature/vbl`, HEAD `54d97b9`)

Adversarial, read-only review of `geomech/core/base`, `geomech/core/operations`,
`geomech/core/transformations`, `geomech/core/math`, `geomech/utils/errors.py` and
`geomech/utils/printing.py`, against the intended semantics in
`docs/manifold-implementation-notes.md` and `paper_context.md`. Nothing under `geomech/`
or `tests/` was modified. State at review: `python -m pytest tests/ -q` → 1322 passed,
5 xfailed.

Every finding below was reproduced with the public API (snippets are verbatim from the
probe scripts; numeric checks use `tests/numeric_helpers.Evaluator`). Severity order:
**silent-wrong-result** > **crash-on-legal-input** > **rule-does-not-fire** > **cosmetic**.

## Summary

| # | Severity | Where | One line |
|---|---|---|---|
| F1 | silent-wrong-result | `core/math/collect.py:72-74` | `Dot(M·v, b×c) → Dot(v, (M·b)×c)` is a false identity; it reaches `compute_eom` and produces a wrong SO(3) equation. A test pins the wrong rule. |
| F2 | silent-wrong-result | `core/transformations/simplify.py:165-166` | `Hat(Vee(X)) → X` fires for any `X`, not only skew-symmetric ones (`Hat(Vee(I)) → I`). A test pins it with a plain `Matrix`. |
| F3 | silent-wrong-result (low reach) | `core/math/extract.py:650-651` | `_vector_from_matrix` returns `ZeroMatrix` for any unhandled node containing the query: `extract_linear_coeff(MMMul(M, Hat(q)), q) → M·0`. Contradicts `paper_context.md`. |
| F4 | silent-wrong-result (low reach) | `core/math/extract.py:430-447` | `_vector_from_vector` accepts nonlinear `Cross`/`MVMul` terms and returns a "coefficient" that still contains the query (`MVMul(Hat(q), q) → Hat(q)`), unlike `SVMul`/`Dot` which raise `NonLinearError`. |
| F5 | crash-on-legal-input | `core/operations/geometry.py:153-165` | `Transpose(R) * v` (and `* M`, `* s`) always raises: `Transpose.__mul__` assumes its operand is a vector. |
| F6 | crash-on-legal-input | `core/math/collect.py:64-75` | `(MVMul, Cross)` branch dead-ends when the target is deeper in `MVMul.right`; `L = ½(Ω×ρ)·J(Ω×ρ)` fails in `compute_eom` with `AlgebraicError`, while the `m·I` version works. |
| F7 | crash-on-legal-input (API) | `core/base/expressions.py:120-150`, `operations/mixins.py:28-62` | `2*x`, `x+1`, `x-1`, `-x` raise (`TypeError`/`AttributeError`); only `x*2` wraps numerics. |
| F8 | rule-does-not-fire | `vector_rules.py`, `pull.py` | `vector_rules` and `pull` do not descend into `Dot`, `Hat`, `Vee`, `Transpose`, `MVMul`, `MMMul`, `SMMul`, `MAdd`, calculus nodes: `Dot(v, q×(ω×q))`, `Hat(a×a)`, `J·(a×a)`, `d/dt(a×a)`, `Hat(v·s)` survive `full_simplify`. |
| F9 | rule-does-not-fire | `simplify.py:177-193`, `expressions.py:527` | `Transpose(I)` never simplifies (`I` is not flagged symmetric and `Transpose` carries no flags); `Transpose(SkewSymmMatrix)`, `MMMul(MMMul(A, Rᵀ), R)` likewise. |
| F10 | rule-does-not-fire | `simplify.py:312-332, 431-465` | Like-term collection misses buried numerics: `(x*2)*y + (x*3)*y` is not merged; `SVMul(v, Mul(s,2)) + SVMul(v, s)` is not merged. |
| F11 | rule-does-not-fire | `expressions.py:537-543`, `simplify.py:440-442` | `Number(0)` is not flagged `Zero` (`Number(0)*x` stays `(0)x`); folding two numerics to 0 also produces an unflagged `(0)`. |
| F12 | API hazard | `simplify.py:234-245` | `Scalar("m", value=1)` is treated as the literal `1` and dropped from products; `Scalar("m", value=2) + Scalar("l", value=3)` becomes `(5)`; the symbol name is discarded. Yet `.is_constant` is `False`. |
| F13 | cosmetic / latent | `expressions.py:20-23, 99`, `multiplication.py` `__str__` | Equality is `type + str`, and `Mul`/`MMMul`/`MVMul`/`SVMul` print children with no separator: `Mul(ab, c) == Mul(a, bc)`. Calculus nodes compare equal across inner leaf types: `Variation(Scalar x) == Variation(Vector x)`. Leaf `has()` has no type check. |
| F14 | API hazard | `core/base/types.py:27-54` | Unknown `attr` strings are silently ignored (`"constant"`, `"Orthogonal"`, even the internally used `"TangentVector"`). |
| F15 | API hazard | `expressions.py:357-362, 491-496` | `S2(...)`/`SO3(...)` append to the caller's `attr` list; a later `Vector("e4", attr=shared)` becomes unit-norm + manifold. |
| F16 | design note | all node classes | "Immutable expression tree" is not enforced: dataclasses are not frozen, `nodes` is a mutable list, `flags` is mutable; hashes change after mutation. |
| F17 | latent | `mixins.py:266-289` | `TimeIntegral(c).is_constant` is `True` for constant `c`. `Transpose(v)*s` builds `Transpose(SVMul(Transpose(v), s))` (LaTeX `s v^T^T`, invalid double superscript). The `Transpose` node does not distinguish row from column. |

Everything else that was targeted (listed in "Checked and clean") held up, including two
numeric fuzz campaigns: 6000 random trees through `expand`/`pull`/`simplify`/`full_simplify`
(only F2 fires) and 10 000 linear-by-construction trees through `extract_linear_coeff`
(zero mismatches for vector and scalar queries).

---

## F1 — `collect` rewrites `Dot(M·vec, b×c)` with a false identity (silent wrong EOM)

**File**: `geomech/core/math/collect.py:64-75`, branch `case (MVMul(), Cross())`, lines 72-74:

```python
if l.right == vec:
    # Dot(M*vec, Cross(b, c)) = Dot(vec, Cross(M*b, c))
    return Dot(vec, Cross(MVMul(l.left, r.left), r.right))
```

The claimed identity `(M v)·(b×c) = v·((M b)×c)` is false for general `M` (it holds only
for `M = sI`). The correct rewrite is `(M v)·(b×c) = v·(Mᵀ(b×c))`, which is exactly what the
generic branch `case (MVMul(), _)` at lines 100-102 already does; the specialised
`(MVMul, Cross)` branch shadows it.

**Reproducer** (unit level):

```python
M, v, b, c = Matrix("M"), Vector("v"), Vector("b"), Vector("c")
expr = Dot(MVMul(M, v), Cross(b, c))
out = collect(expr, v)            # Dot(v,Cross(Mb,c))
ev = Evaluator(); ev(expr), ev(out)   # 0.0816 vs -0.2626  (not equal)
```

Same mismatch with `J` flagged `SymmetricMatrix`. The reference rewrite
`Dot(v, MVMul(Transpose(M), Cross(b, c)))` evaluates equal.

**Reproducer** (pipeline): a rigid body with a gyroscopic coupling linear in `Ω`,
`L = ½ Ω·JΩ + (JΩ)·(q×e₃) + ½ ω·ω`, `J` constant symmetric, `q ∈ S²`, `R ∈ SO(3)`:

```python
eom = compute_eom(L, Zero, SystemVariables(vectors=[q], matrices=[R]))
# engine  η-row: -J Ω̇ - Ω×JΩ - (J(ω×q))×e3 - Ω×J(q×e3)
```

Hand derivation (Euler–Poincaré with `p = ∂L/∂Ω = JΩ + Jᵀ(q×e₃)`, engine convention
`-(ṗ + Ω×p) = 0`): the third term must be `-Jᵀ((ω×q)×e₃)`, not `-(J(ω×q))×e₃`. Numeric
check at a random state satisfying the constraints:

```
engine  : [-3.0936  3.4438  1.2309]
expected: [-3.0779  5.6230  0.7128]   match: False
```

With the branch replaced at runtime by `Dot(vec, MVMul(Transpose(l.left), r))` the engine
output matches `expected` exactly, and a 4000-case numeric fuzz of `collect` (random
`Dot(l, r)` over `Cross`/`MVMul`/`SVMul` with `R, R2 ∈ SO(3)`, symmetric `J`, generic `M`)
goes from 5 mismatching shape classes (all of which route through this branch after
flips) to **0 mismatches**.

**Why no test caught it**: `tests/test_collect.py:204-213`
(`test_mvmul_cross_vec_in_mvmul`) asserts the wrong shape (`result.right.left` is an
`MVMul`), i.e. it pins the bug. The pinned systems in `test_manifold_pipelines.py` never
have `M·ξ̇` on one side of a `Dot` and a `Cross` on the other (for `L = ½Ω·JΩ` the `Dot` is
`Dot(δΩ, JΩ)`, which takes other branches).

**Suggested fix**: replace lines 72-75 by

```python
if l.has(vec):      # covers l.right == vec and deeper (see F6)
    return collect(Dot(l.right, MVMul(Transpose(l.left), r)), vec)
return expr
```

keeping the two `r.left == vec` / `r.right == vec` cases above it (those are correct:
`a·(b×c) = b·(c×a) = c·(a×b)`), and fix the test to expect
`Dot(v, MVMul(Transpose(M), Cross(w, u)))`. `_eliminate` already turns `Transpose(J)` into
`J` for symmetric `J`.

## F2 — `Hat(Vee(X)) → X` without the skew-symmetry hypothesis

**File**: `geomech/core/transformations/simplify.py:161-167`:

```python
case Hat():
    inner = _eliminate(expr.expr)
    ...
    if isinstance(inner, Vee):
        return inner.expr          # Hat(Vee(X)) = X only if X is skew-symmetric
```

`Vee(Hat(v)) → v` (lines 173-174) is always valid; the converse is not. Reproduced:

```python
full_simplify(Hat(Vee(IdentityMatrix)))              # → I     (numerically Hat(Vee(I)) = 0)
full_simplify(Hat(Vee(Matrix("J", attr=["SymmetricMatrix"]))))   # → J
full_simplify(MVMul(Transpose(R), MVMul(Hat(Vee(I)), e3)))       # → (R)'e3, should be 0v
```

In the 6000-tree simplification fuzz (two seeds) **every** numeric mismatch of
`_eliminate`/`simplify`/`full_simplify` (48 + 73 cases) contained `Hat(Vee(X))` with a
non-skew `X` (the one exception, `Hat(Vee(J) + …)`, expands to it). `tests/test_simplify.py:376-380`
(`test_hat_vee_cancel`) pins the rule with a plain `Matrix("M")`.

Reach: the engine never creates a `Vee` node itself (`expand`/`_eliminate`/`extract`
only pass it through), so this is user-built input only, but `Vee` is the documented way
to write SO(3) equations.

**Suggested fix**: fire only when `inner.expr` is a `Hat(...)`, a `SkewSymmMatrix`, or
has `is_skew_symmetric`; add `is_skew_symmetric` properties to `Hat` (always `True`) and
to `SMMul`/`MAdd`/`Transpose` of skew operands. Change the test to use `SkewSymmMatrix`.

## F3 — `_vector_from_matrix` silently returns `ZeroMatrix` for unhandled nodes

**File**: `geomech/core/math/extract.py:604-651`; the fallthrough at 650-651 is
`case _: return ZeroMatrix` with no `_zero_or_raise` guard, unlike every other
extractor in the module.

```python
q, M = Vector("q"), Matrix("M")
extract_linear_coeff(Hat(q), q)             # → 0   (Matrix)
extract_linear_coeff(MMMul(M, Hat(q)), q)   # → M*0
extract_linear_coeff(Transpose(Hat(q)), q)  # → 0
extract_linear_coeff(MAdd(Hat(q), M), q)    # → 0
```

`paper_context.md` lists `MMMul(M, Hat(v))` (matrix expr, vector query) as a supported
row of the dispatch grid; it returns a wrong zero. Mathematically `v ↦ M·Hat(v)` is linear
but has no 3×3 matrix coefficient `C` with `C v = M Hat(v)`, so the honest answer is
`AlgebraicError` (as the module already does for outer products at line 629-634).

Reach: `compute_eom` only extracts from the scalar `δS`; `vbl.py` and
`standard_form.py` extract from scalar/vector rows, so this path is not hit by the
pipeline today. It is a public-API hazard and a documentation error.

**Suggested fix**: `case _: return _zero_or_raise("_vector_from_matrix", expr, query, ZeroMatrix)`
plus explicit `Hat`/`Transpose` cases that raise `AlgebraicError`; correct the table in
`paper_context.md`.

## F4 — nonlinear `Cross`/`MVMul` terms are not rejected by `_vector_from_vector`

**File**: `geomech/core/math/extract.py:426-450`. In the `Cross` branch the
`right == query` / `left == query` shortcuts (430, 432) run *before* the linearity guard
(434-435); in the `MVMul` branch there is no guard at all (445: `if vec == query: return mat`
even if `mat` contains the query; 447: `MMMul(mat, …)` without checking `mat`).

```python
extract_linear_coeff(Cross(Cross(a, q), q), q)             # → hat(cross(a, q))     (contains q)
extract_linear_coeff(Cross(q, Cross(a, q)), q)             # → -hat(cross(a, q))
extract_linear_coeff(MVMul(Hat(q), q), q)                  # → hat(q)
extract_linear_coeff(MVMul(SMMul(M, Dot(q, a)), q), q)     # → M*<q, a>
extract_linear_coeff(MVMul(Hat(q), Cross(a, q)), q)        # → hat(q)*hat(a)
# for contrast, these raise NonLinearError as designed:
extract_linear_coeff(SVMul(q, Dot(q, a)), q); extract_linear_coeff(Dot(q, q), q)
```

Reach: variations enter `δS` linearly by construction, so this needs a user-supplied
`inf_work` that is nonlinear in `ξ`/`η`; it then produces an EOM whose "coefficient"
contains the variation vector, with no error from `compute_eom`.

**Suggested fix**: move the `_has_algebraic(left) and _has_algebraic(right)` guard to the
top of the `Cross` branch and add the same guard to the `MVMul` branch
(`_has_algebraic(mat, query) and _has_algebraic(vec, query)` → `NonLinearError`).

## F5 — `Transpose(matrix) * anything` raises

**File**: `geomech/core/operations/geometry.py:153-165`, `Transpose.__mul__` dispatches
on `other.type` only and assumes `self` is a row vector:

```python
Transpose(SO3("R")) * Vector("v")   # ExpressionMismatchError: VVMul ... MATRIX and VECTOR
Transpose(SO3("R")) * Matrix("M")   # ExpressionMismatchError: MVMul ... MATRIX and MATRIX
Transpose(SO3("R")) * Scalar("s")   # ExpressionMismatchError: SVMul ... MATRIX and SCALAR
MVMul(Transpose(SO3("R")), Vector("v"))   # fine
```

`Rᵀ v` is the most common thing one writes for SO(3) systems. The test suite only builds
it through `MVMul(...)` directly.

**Suggested fix**: in `Transpose.__mul__`, branch on `self.type`: for `MATRIX` delegate to
`_BaseMixin.__mul__` (which already produces `SMMul`/`MVMul`/`MMMul`); keep the row-vector
cases for `VECTOR`. Also the `SCALAR` case returns `Transpose(SVMul(Transpose(v), s))`
(a transposed row, i.e. a column again); return `SVMul(self, other)` instead (see F17).

## F6 — `collect` dead-ends on `Dot(M·f(target), b×c)`, crashing a legal Lagrangian

**File**: `geomech/core/math/collect.py:64-75`. When the IBP target is inside
`MVMul.right` but is not `l.right` itself, the `(MVMul(), Cross())` branch returns `expr`
unchanged (line 75), so the target is never moved to the left, IBP does not fire, and
extraction later raises.

```python
R = SO3("R"); Om = R.get_tangent_vector(); rho = Vector("rho", attr=["Constant"])
J = Matrix("J", attr=["SymmetricMatrix", "Constant"])
L = Dot(Cross(Om, rho), MVMul(J, Cross(Om, rho))) * 0.5
compute_eom(L, Zero, SystemVariables(matrices=[R]))
# AlgebraicError: extract_linear_coeff: '\eta_{R}' appears only inside the calculus operator
#                 '\frac{d}{dt}(\eta_{R})' ...
# whereas L2 = Dot(Cross(Om, rho), Cross(Om, rho)) * (m*0.5) derives fine.
```

The offending `Dot` is `Dot(J·(η̇×ρ), Ω×ρ)` (the second product-rule term of `δL`). The
generic `(MVMul(), _)` branch (lines 99-120) knows how to peel the matrix; the
`(MVMul, Cross)` branch intercepts it first. The fix in F1 (peel with `Transpose(l.left)`
whenever `l.has(vec)` and recurse) resolves this case too. In the collect fuzz, 1844/4000
random `Dot`s were left with the target not on the left; many are legitimately
uncollectable shapes (`SVMul` wrapping that `pull` would remove first), but
`Dot(a, Cross(J·M·t, R·M·b))` and `Dot(J·a, J·M·t·s)` are examples that should collect.

## F7 — mixed numeric arithmetic only works on the right of `*`

**Files**: `geomech/core/base/expressions.py:120-150` (`ScalarExpr.__add__/__sub__/__mul__`),
`geomech/core/operations/mixins.py:28-62`. `__mul__` wraps `int`/`float`; `__add__`,
`__sub__` do not; there is no `__radd__`, `__rmul__`, `__neg__`.

```python
x = Scalar("x")
2 * x     # TypeError: unsupported operand type(s) for *: 'int' and 'Scalar'
x + 1     # AttributeError: 'int' object has no attribute 'type'
x - 1     # AttributeError
-x        # TypeError
x * 2     # Mul(x, (2))  fine
```

**Suggested fix**: `_wrap_numeric` in `__add__`/`__sub__`; add `__radd__ = __add__`,
`__rmul__ = __mul__`, `__neg__ = lambda self: self * -1` on the three `*Expr` bases and
`_BaseMixin`.

## F8 — `vector_rules` and `pull` do not recurse through most node types

**Files**: `geomech/core/transformations/vector_rules.py` (`vector_rules`, lines ≈122-201:
recursion only through `Add`, `Mul`, `VAdd`, `SVMul`, and the generic `Cross` case; the
`Dot` case inspects its children but never rewrites them; `Hat`, `Vee`, `Transpose`, `MVMul`,
`MMMul`, `SMMul`, `MAdd`, `Variation`, `TimeDerivative`, `TimeIntegral` fall to `case _`),
`geomech/core/transformations/pull.py` (no cases for `Hat`, `Vee`, `Transpose`, `VVMul`).
Reproduced with `full_simplify` (10 passes, fixed point reached immediately):

```
Dot(v, Cross(q, Cross(ω_q, q)))   → unchanged   (should be Dot(v, ω_q))
Dot(Cross(q, Cross(ω_q, q)), v)   → unchanged
Hat(Cross(a, a))                  → unchanged   (should be 0)
MVMul(J, Cross(a, a))             → unchanged   (should be 0v)
TimeDerivative(Cross(a, a))       → unchanged
Hat(SVMul(v, s))                  → unchanged   (pull: should be SMMul(Hat(v), s))
Hat(Cross(SVMul(a, s), a))        → unchanged
Transpose(Cross(a, a))            → unchanged
SMMul(J, Dot(q, q))               → unchanged   (but SVMul(v, Dot(q, q)) → v works)
```

Not wrong results, but the pipeline's derived equations can carry terms that are
identically zero or unreduced whenever a `Cross` ends up under `Hat`/`MVMul`/`Dot`/`d/dt`,
and the asymmetry (`SVMul` yes, `SMMul` no) is surprising. `Vee(Hat(Cross(a,a)))` happens
to reduce (to `0v`) only because `_eliminate` cancels `Vee(Hat)` first.

**Suggested fix**: make `vector_rules` and `pull` structurally recursive over all node
types (the same `_rebuild` helper used by `manifold_rules._rebuild` would do), keeping the
rewrite rules at each node.

## F9 — `Transpose(IdentityMatrix)` and friends never simplify

**File**: `geomech/core/transformations/simplify.py:177-193` (`Transpose` case) and
`geomech/core/base/expressions.py:527` (`IdentityMatrix = Matrix("I", attr=["Constant","Identity"])`,
no `SymmetricMatrix`). `_is_identity` reads `flags`, which `Transpose` nodes do not have.

```
full_simplify(MVMul(Transpose(IdentityMatrix), v))    → (I)'v
full_simplify(Transpose(SkewSymmMatrix("S")))         → (S)'      (should be -S)
full_simplify(MMMul(MMMul(A, Transpose(R)), R))       → A(R)'R    (MVMul has associativity, MMMul does not)
full_simplify(MMMul(Transpose(R), MMMul(R, A)))       → (R)'RA
```

`Transpose(I)` is produced by `_vector_from_scalar` (`MVMul(Transpose(coeff), right)`)
whenever the query sits under a wrapper on one side of a `Dot`, e.g.
`full_simplify(extract_linear_coeff(Dot(SVMul(q, s), w), q))` → `(I)'ws` instead of `ws`,
and by `vbl.py` (`coeff = Transpose(coeff)` for scalar rows with vector symbols). **Fix**: flag `IdentityMatrix` symmetric (and treat
`_is_identity(Transpose(X))` as `_is_identity(X)`); add `Transpose(S) → SMMul(S, -1)` for
`is_skew_symmetric`; add the `MMMul` reassociation mirror of the `MVMul` rule.

## F10 — like-term collection misses buried numeric factors

**File**: `geomech/core/transformations/simplify.py:312-332` (`_split_coeff` looks only at
a top-level numeric child) and `431-465` (`_simplify_mul` refolds a chain only when it
holds ≥ 2 numerics).

```
full_simplify((x*2)*y + (x*3)*y)                 → (x(2)y+x(3)y)   expected x y (5)
full_simplify(x*(y*2) + x*(y*3))                 → (xy(2)+xy(3))
full_simplify(SVMul(v, Mul(s,2)) + SVMul(v, s))  → (vs(2)+vs)
full_simplify(x*2 + x*3)                         → x(5)            (works)
```

**Fix**: in `_split_coeff`, flatten the `Mul` chain (and the scalar of `SVMul`/`SMMul`),
pull out all numeric leaves, and rebuild the remainder; or always normalise `Mul` chains
to "symbols first, single numeric last" in `_simplify_mul` (drop the `len(numerics) >= 2`
condition).

## F11 — `Number(0)` is not a zero

**File**: `geomech/core/base/expressions.py:537-543` (`Number` sets only `Constant`),
`simplify.py:440-442` (numeric product folded to `0` without the `Zero` flag).

```
Number(0).is_zero                          → False
full_simplify(Number(0) * x)               → (0)x
full_simplify(SVMul(v, Number(0)))         → v(0)
full_simplify(Mul(Number(0), Number(3)))   → (0)   with .is_zero False
Number(1).is_ones                          → False  (harmless: _is_one uses .value)
```

Inside an `Add` the term is dropped by `_make_term(coeff=0)`, so only standalone products
are affected. **Fix**: in `Number`/`_wrap_numeric`, add `"Zero"` when `value == 0` and
`"Ones"` when `value == 1`; in `_simplify_mul` return `Zero` when `combined == 0`.

## F12 — `Scalar(name, value=…)` is folded as a literal and loses its name

**File**: `simplify.py:234-245` (`_is_one`, `_is_numeric_leaf` key on `.value`), used by
`_simplify_mul`, `_simplify_add`, `_eliminate`.

```
full_simplify(Scalar("m", value=1) * x)                        → x
full_simplify(Scalar("m", value=2) + Scalar("l", value=3))     → (5)
full_simplify(Scalar("m", value=2) * Scalar("l", value=3))     → (6)
full_simplify(Scalar("g", value=9.81) + Scalar("g", value=9.81)) → (19.62)
Scalar("m", value=2).is_constant → False   but .t_diff() → 0 and .delta() → 0
```

The docs use `value=` only for `Scalar("0.5", value=0.5)` and `Vector("e3", value=np.array(...))`,
so the folding is probably intended for anonymous numbers only. A user attaching a numeric
value to a named parameter (the `Evaluator` in `tests/numeric_helpers.py` reads
`node.value` exactly for that) silently loses the symbol. **Fix**: fold only when the name
is the printed number (i.e. created by `Number`/`_wrap_numeric`), or document that
`value` makes the scalar a literal; make `is_constant` true for numeric leaves either way.

## F13 — string equality: separator-free printing and type-blind calculus equality

**Files**: `geomech/core/base/expressions.py:20-23` (`__eq__`: same class and same
`str`), `:99` (`has` on a leaf compares `str` only), `core/operations/multiplication.py`
(`Mul`, `SVMul`, `SMMul`, `MVMul`, `MMMul`, `VVMul.__str__` are `str(left) + str(right)`).

```
Mul(Scalar("ab"), Scalar("c")) == Mul(Scalar("a"), Scalar("bc"))   → True   (both "abc")
MMMul(Matrix("AB"), Matrix("C")) == MMMul(Matrix("A"), Matrix("BC")) → True
Variation(Scalar("x")) == Variation(Vector("x"))                   → True   (both "\delta{x}")
TimeDerivative(Scalar("m")) == TimeDerivative(Vector("m"))         → True
Scalar("x").has(Vector("x")); Dot(Vector("m"), w).has(Scalar("m")) → True, True
S2("q") == Vector("q") → False (type-guarded);   Add(a, b) == Add(b, a) → False (documented)
```

Consequence that was reproduced: a scalar configuration variable `x` alongside a vector
named `x` makes extraction fail:

```python
expr = Add(Mul(Variation(Scalar("x")), k), Dot(Variation(Vector("x")), F))
extract_linear_coeff(expr, Variation(Scalar("x")))
# AlgebraicError: _scalar_from_vector: unsupported node Variation contains '\delta{x}' in '\delta{x}'
```

The multiplication collision needs names that concatenate ambiguously (`m`, `l`, `ml`), so
it is unlikely in practice but it is a correctness hole in the one primitive every rewrite
relies on. **Fix**: print `Mul`-family nodes with a separator (`*`) or, better, make
`__eq__`/`__hash__` structural (class, child equality, leaf name + class); give `Expr.has`
the same type guard `__eq__` has. `_same` in `extract.py` inherits the calculus-node
blindness because it falls back to `str` for nameless nodes.

## F14 — unknown `attr` strings are silently ignored

**File**: `geomech/core/base/types.py:27-54` (`match a:` with no default).

```python
Scalar("m", attr=["constant"]).is_constant   → False    (also "Const", "Orthogonal", "unitnorm")
Scalar("m", attr=["TangentVector"])           → no flag  (the string TS2/TSO3 append themselves)
compute_eom(Dot(ẋ, ẋ) * Scalar("m", attr=["constant"]) * 0.5, Zero, SystemVariables(vectors=[x]))
# → (d/dt(ẋ)*-m + ẋ*-ṁ)   vs   d/dt(ẋ)*-m with "Constant"
```

The result is visible (`ṁ` appears) rather than silent, but it is the kind of typo that
costs an afternoon. **Fix**: `case _: raise ValueError(f"unknown attribute {a!r}")` and add
`"TangentVector"` to the accepted set.

## F15 — `S2`/`SO3` constructors mutate the caller's `attr` list

**File**: `geomech/core/base/expressions.py:357-362` (`S2.__init__`: `attr.append(...)` on
the list passed in), `:491-496` (`SO3`). `Matrix.__init__` keeps a reference too.

```python
shared = ["Constant"]
e3 = Vector("e3", attr=shared); q = S2("q", attr=shared)
shared                       → ['Constant', 'Manifold', 'UnitNorm']
e4 = Vector("e4", attr=shared)
e4.is_unit_norm, e4.is_manifold   → True, True      (and e4.delta() → 0 because Constant)
```

**Fix**: `attr = list(attr) if attr else []` before appending.

## F16 — the tree is not actually immutable

All node dataclasses are `@dataclass(eq=False, repr=False)` without `frozen=True`; `nodes`
is a plain list; `flags` is a mutable dataclass.

```python
e = Add(a, b); h = hash(e); e.nodes.append(Scalar("c")); e → (a+b+c), hash changed
a.name = "zz"; e → (zz+b+c)
v = Vector("v"); v.flags.is_constant = True; v.delta() → 0
```

Nothing in the engine relies on this today (no memoisation keyed by node), so this is a
design note: either freeze (`tuple` nodes, `frozen=True`, immutable flags) or drop the
"immutable" claim from the contributor notes.

## F17 — minor latent issues

- `mixins.py:266-289`: `_CalcUnaryMixin` inherits `_UnaryMixin.is_constant`, so
  `TimeIntegral(Scalar("c", attr=["Constant"])).is_constant → True`. No rewrite currently
  uses it on an integral, but `_CalcUnaryMixin.t_integrate` would return `Zero` for it.
- `geometry.py:159`: `Transpose(v) * s → Transpose(SVMul(Transpose(v), s))`; `repr` is
  `v'*s'`, `to_latex` gives `s v^T^T` (double superscript, will not compile).
- The `Transpose` node type is used as the row/column marker (`VVMul.__init__`,
  `MVMul.__init__`, `Transpose.__mul__` all test `isinstance(x, Transpose)`), so
  `Transpose(Transpose(v))`, a column, is treated as a row until `vector_rules` cancels
  it. Fuzzing with that shape produced an `ExpressionMismatchError` from
  `_scalar_from_matrix` (`VVMul(Transpose(...), Transpose(w))`); it needs an ill-typed input
  to trigger, so it is listed here rather than as a finding.
- `_simplify_add` prints float sums verbatim: `Number(0.1) + Number(0.2) → (0.30000000000000004)`.
- `vector_rules.py:35-45, 97-115`, `_is_tangent_orthogonal` / `_is_tangent_to_s2` use instance
  identity (`a.S2 is b`) while everything else uses name equality; a second `S2("q")`
  instance (equal under `==`) does not get `Dot(ω, q) → 0` or the BAC-CAB rules
  (`full_simplify(Dot(ω_other, q)) → <ω_q, q>`). Missed simplification only; the EOM in
  the probed case came out identical.
- `compute_eom(L_built_with_S2("q"), Zero, SystemVariables(vectors=[Vector("q")]))` returns
  `{'\\delta{q}': '0v'}` silently (the variation `Variation(Vector q)` never appears in
  `δS`). A check that every configuration variable's variation vector occurs in `δS`
  would catch this class of user error.

---

## Checked and clean

The following were targeted explicitly and no defect was found; the absence of findings
there is backed by the probes described.

**Numeric soundness of the rewrite passes** (`expand`, `pull`, `vector_rules`, `_eliminate`,
`simplify`, `full_simplify`): 2 × 3000 random trees (depth ≤ 4; scalar, vector and matrix
typed) over leaves `q ∈ S²` with its `ω`, `ξ` (values drawn ⊥ `q`), `R ∈ SO(3)` with `Ω`, `η`,
symmetric constant `J`, generic `A`, skew `K`, `I`, zero constants, constant `e3`, numbers
`2, -1, 0.5, 1, 0, -2.0`, and `d/dt` of leaves. Each pass's output was evaluated against the
input with all leaves pre-seeded. Result: 0 exceptions, 0 non-convergence (fixed point of
`simplify` always reached in ≤ 10 passes; no oscillating rule pair was found), and every
mismatch was F2. In particular the S² rules (`Dot(q,q)=1`, `Dot(ω,q)=0`,
`q×(a×q)=a`, `q×(q×a)=-a`, `(q×a)×q=a`, `(a×q)×q=-a`, `a×(a×q)=-q(a·a)`) only fire with the
unit-norm `S2` instance and a structurally tangent `a`, and are numerically correct under
those hypotheses; `RᵀR = I`/`RRᵀ = I` require the `OrthogonalMatrix` flag on the exact leaf;
`Transpose(J) → J` requires `SymmetricMatrix`; `Transpose(Hat(v)) → -Hat(v)`,
`Hat(v)·w → v×w`, `(AB)v → A(Bv)`, `Transpose(AB) → BᵀAᵀ`, `Transpose(sM) → s Mᵀ` are
unconditional and correct.

**Flag propagation**: `is_zero`/`is_constant` on `_BinaryMixin` (`or`/`and`), `_NaryMixin`
(`all`/`all`), `_UnaryMixin` and `Transpose` are sound (confirmed by the fuzz: no zero was
ever manufactured). `is_unit_norm`, `is_orthogonal`, `is_symmetric`, `is_skew_symmetric`,
`is_identity` are leaf-only (`Expr._flag` reads `self.flags`, which operation nodes lack),
so `Transpose(R)`, `SMMul(R, s)`, `Cross(a, q)`, `MMMul(R1, R2)` never claim a property they
might not have. Conservative, never unsound (but see F9 for the cost).

**`extract_linear_coeff`**: 3000 vector-typed + 3000 scalar-typed random trees containing a
vector query exactly once (shapes: `Cross` either side, `MVMul` with `R`/`J`/`M`, `SVMul`,
`VAdd`, `Hat(·)·v`, `Transpose(Hat(·))·v`, outer-product matrices, `SVMul(vec, Dot(·,·))`,
`Add`, `Mul`, `Dot`, `VVMul` inner products) and 4000 scalar/vector/matrix trees containing
a scalar query exactly once (`Mul`, `Add`, `Dot`, `VVMul`, `SVMul`, `Cross`, `MVMul`, `VAdd`,
`Vee`, `SMMul`, `MMMul`, `Hat`, `Transpose`, `MAdd`, outer products). Check: `C·q` (resp.
`c·q`, `c·z`) equals `expr(2q) − expr(q)`. **0 mismatches, 0 spurious raises.** So the
Cross/Hat sign conventions (`l×q → Hat(l)`, `q×r → -Hat(r)`, `Hat(q)v → -Hat(v)`,
`Hat(q)ᵀv → Hat(v)`), the `SVMul`-scalar redistribution (`Mul`, `Dot`, `VVMul`, `Add`), the
outer-product rule and the row-vector (`Transpose(v)·M`) paths are all right. Calculus nodes
are opaque and a query found only inside one raises `AlgebraicError` rather than returning
zero (`_has_algebraic`, verified). `Dot`/`SVMul`/`SMMul`/`Mul`/`VVMul`/`MMMul` all raise
`NonLinearError` when the query is on both sides (only F4's two branches do not).

**`collect`**: all rewrites except F1 were verified numerically (after patching F1 at
runtime, 4000/4000 exact): Dot flip, `a·(vec×c) = vec·(c×a)`, `a·(b×vec) = vec·(a×b)`,
`M(vec×b)·w = vec·(b×Mᵀw)`, `M(b×vec)·w = vec·(Mᵀw×b)`, `(M vec)·w = vec·(Mᵀw)`, nested
`f(vec)` inside `Cross`, and `Dot(R a, R b) → Dot(a, b)` which correctly requires
`isinstance(SO3)` on both sides with the same name (a `Matrix` merely flagged
`OrthogonalMatrix`, or two different `SO3`s, are not stripped). The documented
"flip-only-if-other-side-has-vec" guard prevents the `(Cross, MVMul)`/`(MVMul, Cross)`
ping-pong; no non-termination was observed in 8000 random calls.

**Integration by parts** (`ibp.py`): `∫ a·Dot(ẋ, f) = -∫ Dot(x, ȧ f + a ḟ)` is implemented
for non-constant prefactors (`Mul(s, Dot(ξ̇, f))`, prefactor chains with repeated factors
`s·s → 2 s ṡ`, prefactors containing `Dot(q, e3)` of an S² point so that `q̇ = ω×q` is used);
signs are right; `TimeDerivative(ξ).t_integrate() → ξ`; scalar-chain IBP removes the exact
factor and differentiates the remainder; `TimeDerivative(x).delta()` and
`Variation(x).t_diff()` build the identical tree `d/dt(δx)` so the IBP target matches. An
uncollected `Dot(f, ξ̇)` is left alone and surfaces as an `AlgebraicError` in extraction
(never silently dropped). A nonlinear `Mul(Dot(ξ̇,f), Dot(ξ̇,f))` produces `d/dt(d/dt(η))`
and likewise errors later.

**`manifold_rules`**: `δω = ξ̇ − ω×ξ` was re-derived from `q̇ = ω×q`, `δq = ξ×q` and the
Jacobi identity; the component along `q` vanishes because `ξ̇·q = q·(ω×ξ)` (so no `λq`
term), i.e. the formula is correct under `ω ⊥ q`, `ξ ⊥ q`, which `TS2` guarantees. `δΩ = η̇ + Ω×η`
matches `Ṙ = RΩ̂`, `δR = Rη̂`. Substitution is keyed by the printed tangent-vector name, so
`q` and `q1` (`\omega_{q}` vs `\omega_{q1}`) cannot collide; a second application is a no-op
because the result contains no `Variation`; `Variation(ξ)`/`Variation(η)` are rejected at
construction; `Variation(ω)` for an `S2` absent from `SystemVariables` is left untouched
(and then visibly survives into any coefficient it shares with a listed variation vector).
User-built `Variation(SVMul(ω, s))` is not expanded (the pipeline never produces it since
`.delta()` pushes `δ` to the leaves).

**Leaves**: `S2.t_diff()`/`.delta()` and `SO3.t_diff()`/`.delta()` are eager and follow Lee
et al.'s conventions; `__hash__` is consistent with `__eq__` (`Scalar("x")` and `Vector("x")`
hash equal but compare unequal, which is allowed); `Vector.delta()`/`t_diff()` of
`Constant` leaves return flagged zeros; `Scalar("x") == Vector("x")` is `False`.

**Printing**: `Dot`, `Cross`, `Hat`, `Vee`, `Transpose`, `Add` and the calculus nodes print
with delimiters, so the only separator-free forms are the multiplication family (F13).
`to_latex`/`repr_str` handle `-1` and `0.5` folding consistently.

---

## Tests that currently pin wrong behaviour

- `tests/test_collect.py:204-213` `test_mvmul_cross_vec_in_mvmul` — asserts F1's false
  identity shape.
- `tests/test_simplify.py:376-380` `test_hat_vee_cancel` — asserts `Hat(Vee(M)) == M` for a
  generic `Matrix("M")` (F2).

Per the project's testing convention, fixing F1/F2 requires changing these two tests to
the mathematically correct expectation (with a numeric cross-check, as in
`tests/numeric_helpers.py`), not weakening them.

## Suggested priority

1. F1 (wrong EOM, reachable from any `Dot(M·δ-derivative, Cross)` term) together with F6
   (same branch), then re-run the pinned pipeline tests and add a gyroscopic-coupling
   pipeline test with a numeric check.
2. F2 and F4 (make the hypotheses explicit; cheap).
3. F5, F7 (operator ergonomics; one-line fixes).
4. F3, F11, F14, F15 (error instead of silence).
5. F8, F9, F10 (simplifier completeness) and F12/F13/F16 (semantics worth deciding
   before the paper is written).
