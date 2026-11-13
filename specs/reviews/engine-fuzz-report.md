# Engine fuzz report: numerical property testing of the symbolic core

Branch `feature/vbl` at `54d97b9` (tests: 1322 passed, 5 xfailed). No source file under
`geomech/` or `tests/` was modified; every finding below was reproduced against the
unmodified engine with the snippet shown.

## How to rerun

```bash
cd pygeomech && source venv/bin/activate
python specs/reviews/scripts/engine_fuzz.py                 # 3000 cases per campaign, ~30 s
python specs/reviews/scripts/engine_fuzz.py --cases 5000    # the run reported here, 44 s
python specs/reviews/scripts/engine_fuzz.py --quick          # 200 cases, ~2 s smoke test
```

The script prints a per-transform table and the minimised findings, and writes
`engine_fuzz_findings.json` next to itself (`--json PATH` to redirect, `--no-json` to skip,
`--seed N` to change the seed; the run below used the default seed `20261113`).

## What was fuzzed and how

**Numeric model.** Every leaf is a curve `x(t, ε)` stored as a truncated bivariate Taylor
series (a "jet", order 5 in `t`, 2 in `ε`). Jet products are truncated Cauchy products, so
`TimeDerivative` and `Variation` nodes are evaluated *exactly* as coefficient shifts; no
finite differences are involved. Manifold leaves are sampled as curves that satisfy their
constraints to all orders:

| leaf | construction | consequence |
|---|---|---|
| `S2` point `q` | `p/|p|` for a random jet `p` | `|q| = 1` to all orders |
| `TS2` ω, ξ of `q` | `ω := q × ∂ₜq`, `ξ := q × ∂_ε q` | `q̇ = ω×q`, `δq = ξ×q`, `ω⊥q`, `ξ⊥q` exactly |
| `SO3` `R` | Cayley `R₀(I−A)⁻¹(I+A)`, `A = Hat(d)`, `d` nilpotent | `RᵀR = I` to all orders |
| `TSO3` Ω, η of `R` | `Ω := Vee(Rᵀ∂ₜR)`, `η := Vee(Rᵀ∂_εR)` | `Ṙ = R·Hat(Ω)`, `δR = R·Hat(η)` exactly |
| `UnitNorm` vector `u` | normalised jet | `|u| = 1` |
| `SymmetricMatrix` / `SkewSymmetricMatrix` / `OrthogonalMatrix` | `A+Aᵀ`, `A−Aᵀ`, Cayley | property holds |
| `Constant` / numeric / `Zero` / `Identity` | only the `(0,0)` coefficient | derivative is zero |

In **adversarial** mode the flagged leaves violate their declared properties (ω, ξ random and
not ⊥ q; `|u| ≠ 1`; `R` a random matrix), which shows exactly which rewrite rules rely on
which declared property. The evaluator is self-tested against `tests/numeric_helpers.Evaluator`
on 300 calculus-free trees and against the manifold identities above before every run.

**Generator.** Typed random trees over all node classes (`Add/VAdd/MAdd`, `Mul`, `SVMul`,
`SMMul`, `MVMul` incl. the row-vector form `MVMul(Transpose(v), M)`, `MMMul`, `VVMul` inner
and outer, `Dot`, `Cross`, `Hat`, `Vee`, `Transpose`, `TimeDerivative`, `Variation`), the
leaves above, numeric constants, and S²/SO(3) "motifs" that exercise every rule in
`vector_rules.py` and `simplify.py` both under their hypotheses (`Cross(q, Cross(ω, q))`,
`Dot(q, ξ)`, `RᵀR`, `Hat(Vee(Hat(v)))`, …) and outside them (`Cross(q, Cross(a, q))` with a
generic `a`, `Cross(q, Cross(ω_p, q))` with the tangent vector of a *different* point,
`Cross(u, Cross(a, u))` with a unit-norm vector that is not an `S2`, `MᵀM`, `Hat(Vee(M))`).
Depth 1–4. Trees whose nested `d/dt` depth exceeds the exact margin of the jets are skipped
(column `unsup`). `TimeIntegral` is excluded (not evaluable pointwise).

**Oracle.** `value(before) ≈ value(after)` at `(t, ε) = (0, 0)` with relative tolerance
`1e-8` (scale `max(1, |before|, |after|)`). On a mismatch the tree is shrunk greedily (replace
by a same-typed descendant, drop an addend, replace a child by a fresh leaf) while the
mismatch persists. Each case runs under a 10 s alarm; any exception that is not one of the
types in `geomech/utils/errors.py` is recorded as a crash, and a documented error raised *by
a transformation on a legal tree* is recorded separately.

**`extract_linear_coeff`.** Trees linear in a chosen query are built by construction (query
leaves: `ξ`, `η`, a plain `Vector`, a plain `Scalar`, a plain `Matrix`, `Variation(b)`,
`Variation(s)`; the generic sub-trees are regenerated until they are free of the query). The
check is `expr(query = v) − expr(query = 0) == coeff·v` with the documented convention per
cell (`c·s`, `qᵀc`, `C·q`). For the two cells without a numeric convention (matrix query;
vector query with matrix expression) the check is that a term depending on the query is never
returned as a zero coefficient. A separate campaign builds trees of degree 2 in the query and
expects `NonLinearError`; when something is returned instead, the tree is verified to be
genuinely quadratic numerically (`expr(2v) − expr(0) ≠ 2(expr(v) − expr(0))`).

### Case counts (run of 2026-11-13, seed 20261113, 73,535 cases, 43.6 s)

| campaign | cases | value-preserving / ok | skipped | documented errors raised by the transform |
|---|---:|---:|---:|---|
| `expand` | 5000 | 4999 | 0 | ExpressionMismatchError: 1 (finding 6) |
| `pull` | 5000 | 5000 | 0 | – |
| `vector_rules` | 5000 | 5000 | 0 | – |
| `simplify` | 5000 | 4951 | 0 | ExpressionMismatchError: 6 (finding 6) |
| `full_simplify` | 5000 | 4954 | 0 | ExpressionMismatchError: 7 (finding 6) |
| `full_simplify` reaches a fixed point | 5000 | 4993 (7 = the exceptions above) | 0 | – |
| `expand` / `pull` [adversarial] | 1666 / 1666 | 1665 / 1666 | 0 | 1 / – |
| `vector_rules` [adversarial] | 1666 | 1551 | 0 | – (115 flag-dependent, finding 12) |
| `simplify` / `full_simplify` [adversarial] | 1666 / 1666 | 1440 / 1433 | 0 | 3 / 5 |
| `Variation(e)` vs `e.delta()` | 3400 | 3400 | 0 | – |
| `TimeDerivative(e)` vs `e.t_diff()` | 3400 | 3394 | 6 | – |
| `apply_manifold_rules` [full δω] | 5000 | 3698 | 0 | – (finding 2) |
| `apply_manifold_rules` [tangential δω] | 5000 | 5000 | 0 | – |
| `collect` [target occurs once] | 4069 | 3826 | 0 | – (finding 1) |
| `collect` [target occurs several times] | 931 | 796 | 0 | – (findings 1, 5) |
| `ibp._t_diff_product` | 3405 | 3401 | 4 | – |
| `extract` [scalar query, scalar expr] | 628 | 628 | 0 | – |
| `extract` [scalar query, vector expr] | 566 | 566 | 0 | – |
| `extract` [scalar query, matrix expr] | 213 | 213 | 0 | – |
| `extract` [vector query, scalar expr] | 1281 | 1240 | 0 | AlgebraicError: 41 (finding 13) |
| `extract` [vector query, vector expr] | 1167 | 1013 | 0 | AlgebraicError: 154 (finding 13) |
| `extract` [vector query, matrix expr] | 435 | 44 | 0 | AlgebraicError: 97 (findings 8, 13) |
| `extract` [matrix query, matrix expr] | 111 | 70 | 0 | AlgebraicError: 33 (findings 8, 13) |
| `extract` [matrix query, scalar/vector expr] | 308 / 291 | 0 / 0 | 0 | NotImplementedError (finding 9) |
| `extract` quadratic trees, scalar query | 694 + 291 | 985 raised `NonLinearError` | 0 | – |
| `extract` quadratic trees, matrix query | 1014 | 1014 raised `NonLinearError` | 0 | – |
| `extract` quadratic trees, vector query, scalar expr | 1209 | 1209 raised `NonLinearError` | 0 | – |
| `extract` quadratic trees, vector query, vector expr | 1792 | 1573 | 0 | AlgebraicError: 50; 169 missed (finding 7) |

No timeout, no `RecursionError` outside `collect`, and no `AttributeError`/`TypeError`/
`KeyError` occurred anywhere. 798 raw findings were recorded; after shrinking they reduce to
the 13 root causes below (the in-process confirmation runs mentioned under findings 1 and 7
show that the two largest clusters each have a single cause).

## Findings

Severity legend: **wrong-value** (a value-preserving transformation changes the value),
**crash** (an exception on a legal tree), **rule-fires-outside-hypothesis**, **cosmetic**.

### 1. `collect` rewrites `Dot(M·v, b×c)` as `Dot(v, (M·b)×c)` — wrong-value, reaches `compute_eom`

`geomech/core/math/collect.py:72-74` (the `l.right == vec` branch of `case (MVMul(), Cross())`):

```python
if l.right == vec:
    # Dot(M*vec, Cross(b, c)) = Dot(vec, Cross(M*b, c))
    return Dot(vec, Cross(MVMul(l.left, r.left), r.right))
```

The identity is false: `(Mv)·(b×c) = v·Mᵀ(b×c)`, and `Mᵀ(b×c) ≠ (Mb)×c` for any `M ≠ I`.

```python
from geomech import *
M, v, b, c = Matrix("M"), Vector("v"), Vector("b"), Vector("c")
collect(Dot(MVMul(M, v), Cross(b, c)), v)      # -> Dot(v,Cross(Mb,c))
# jets, random state:  before = 6.558   after = -4.418
```

All 243 `collect` mismatches on single-occurrence trees, and all 131 non-crashing mismatches
on multi-occurrence trees, are this branch, reached directly or after the flips at lines 62,
95, 119, 124, 143. Confirmation: with this one branch replaced in-process by
`Dot(vec, MVMul(Transpose(l.left), r))`, a 2000-case rerun of the campaign gives 1630/1631
single-occurrence cases value-preserving (the remaining case raised a documented error) and
no value mismatch at all on the 369 multi-occurrence cases (only the recursion crashes of
finding 5 remain).

Pipeline reachability: `integrate_by_parts` calls `collect` on every term, so any Lagrangian
whose variation contains `Dot(J·ξ̇, b×c)` is affected. Example, S² with a gyroscopic coupling:

```python
q = S2("q"); w = q.get_tangent_vector(); xi = q.get_variation_vector()
J = Matrix("J", attr=["SymmetricMatrix", "Constant"]); e3 = Vector("e_3", attr=["Constant"])
integrate_by_parts(Dot(MVMul(J, TimeDerivative(xi)), Cross(e3, q)), [TimeDerivative(xi)])
# -> (-1)Dot(xi, Cross(J e_3, Cross(omega, q)))      correct: -Dot(xi, J (e_3 x (omega x q)))
compute_eom(Dot(MVMul(J, w), Cross(e3, q)), Zero, SystemVariables(vectors=[q]))
# -> xi: Cross(J e_3, Cross(omega,q))(-1) + ...       (first term wrong)
```

Numerically (tangential components, which are the only ones that matter for a coefficient of
ξ ⊥ q): `P·X_pkg = [-0.142, -0.172, -0.149]` vs `P·X_true = [0.471, -0.048, -0.533]`.
The existing suite does not hit this branch because its Lagrangians only produce
`Dot(ξ̇, J·Ω)`-shaped terms (`MVMul` on the *right* of the `Dot`), which go through the
correct `(Cross(), MVMul())`/`(MVMul(), _)` paths.

### 2. The S² substitution `δω → ξ̇ − ω×ξ` is only the tangential part of δω — wrong-value when `∂L/∂ω` has a component along `q`

`geomech/core/transformations/manifold_rules.py:65-70` (and the same formula in
`TS2.delta(substitute=True)`, `geomech/core/base/expressions.py:349`, and the module
docstring). `docs/manifold-implementation-notes.md` lines 3-10 and 72-77 list this formula as
"needs verification"; this is the verification.

With `q̇ = ω×q`, `δq = ξ×q`, `ω⊥q`, `ξ⊥q` (the hypotheses the engine enforces through
`vector_rules`), the angular velocity of the varied curve is `ω_ε = q_ε × q̇_ε`, and its
variation is

    δω = ξ̇ − 2 ω×ξ,            with  ξ̇ − ω×ξ = (I − qqᵀ) δω.

The docstring derivation ("set δ(q̇) = d/dt(δq) and solve for δω") determines δω only modulo
`q`; the component along `q` follows from `ω·q ≡ 0`: `q·δω = −ω·δq = −q·(ω×ξ)`, and `ω×ξ`
is parallel to `q`. Two independent numerical checks (exact jets, and plain central finite
differences on `q(t,ε) = p(t,ε)/|p(t,ε)|` with `ω := q×q̇`, `ξ := q×∂_ε q`) both give:

```
delta omega           : [0.0915 0.2125 0.7391]
xi_dot - omega x xi   : [0.1004 0.2014 0.7410]   |err| = 1.1e-2
xi_dot - 2 omega x xi : [0.0915 0.2125 0.7391]   |err| = 1.4e-9
difference is along q : |(diff) x q| < 1.1e-9
```

(`python specs/reviews/scripts/engine_fuzz.py` prints the jet version as its first
self-test line.) The SO(3) formula `δΩ = η̇ + Ω×η` is exact (checked the same two ways).

Consequence. The pipeline forms `Dot(∂L/∂ω, δω)`; the missing term `−ω×ξ ∥ q` is
annihilated whenever `∂L/∂ω ⊥ q`, which is the case for every Lagrangian that depends on ω
only through `q̇ = ω×q` (all the pendulum systems in the test suite: `|ω|²`, `ẋ·(ω×q)`, …).
It is **not** annihilated for a Lagrangian written directly in ω with a generic inertia:

```python
q = S2("q"); w = q.get_tangent_vector(); J = Matrix("J", attr=["SymmetricMatrix", "Constant"])
eom = compute_eom(Mul(Dot(w, MVMul(J, w)), 0.5), Zero, SystemVariables(vectors=[q]))
# xi: (J d/dt(omega) (-1.0) + Cross(omega, J omega))
```

versus the Euler–Lagrange equation of the same Lagrangian written as
`L(q, q̇) = ½ (q×q̇)·J(q×q̇)`, evaluated with exact jets in the Lee et al. form
`q × (d/dt ∂L/∂q̇ − ∂L/∂q) = 0` (the gyroscopic term carries a factor 2):

```
tangential part of package EOM : [ 1.0446 -0.6697 -0.3814]
tangential part of true EOM    : [ 1.0417 -0.6731 -0.3774]   (= -(q x (d/dt dL/dqdot - dL/dq)))
```

Campaign evidence: 1302/5000 mismatches with the full δω, 0/5000 when `Variation(ω)` is
evaluated as the tangential projection `(I − qqᵀ)δω`. Either the substitution should be
`ξ̇ − 2ω×ξ` (then `vector_rules` must not assume `δω ⊥ q`), or the documentation should state
that the kinematic substitution is valid only for Lagrangians with `∂L/∂ω ⊥ q`.

### 3. `simplify` rewrites `Hat(Vee(M))` as `M` for any matrix — rule-fires-outside-hypothesis (wrong-value)

`geomech/core/transformations/simplify.py:165-166` (`_eliminate`, `case Hat()`):

```python
if isinstance(inner, Vee):
    return inner.expr
```

`Hat(Vee(M)) = M` holds only for skew-symmetric `M` (`Vee` reads three entries). No flag is
checked; `Vee` accepts any matrix.

```python
M = Matrix("M")
simplify(Hat(Vee(M)))      # -> M        (value: skew part of M   vs   M)
simplify(Hat(Vee(I)))      # -> I        (value: 0 vs I)
simplify(Hat(Vee(Matrix("J", attr=["SymmetricMatrix", "Constant"]))))   # -> J  (value 0 vs J)
full_simplify(MVMul(Hat(MVMul(I, Vee(M))), Vector("z")))               # -> Mz  (value: skew(M) z)
```

44 of the 47 `simplify`/`full_simplify` mismatches are this rule (the other 3 are finding 4).
The companion rule `Vee(Hat(v)) → v` (line 173-174) is always valid and never failed. The rule
is safe for `Hat(Vee(X))` when `X` is a `Hat(...)`, a `SkewSymmMatrix` leaf, or `MMMul(Rᵀ, Ṙ)`
with `R` an `SO3`; it should be restricted to those.

### 4. `simplify` turns an outer product with a zero factor into the *scalar* `Zero` — wrong type, crashes downstream

`geomech/core/transformations/simplify.py:141-144` (`_eliminate`, `case VVMul()`):

```python
if _is_zero(l) or _is_zero(r):
    return Zero          # scalar, also for the outer product v wᵀ (type MATRIX)
```

```python
a = Vector("a")
simplify(VVMul(ZeroVector, Transpose(a))).type   # -> ExprType.SCALAR  (should be a 3x3 zero)
MAdd(simplify(VVMul(ZeroVector, Transpose(a))), Matrix("M"))
# -> ExpressionMismatchError: MAdd cannot be performed between MATRIX and SCALAR
```

Use `ZeroMatrix` when `expr.type == ExprType.MATRIX`.

### 5. `collect` recurses forever when the target occurs on both sides of a `Dot` — crash (`RecursionError`)

`geomech/core/math/collect.py:127-129` and `139-145`: `(Cross(), _) if l.has(vec)` flips the
`Dot`, and the `(_, Cross())` branch flips it back whenever the *other* side also contains the
target but not as a direct cross-product operand.

```python
a, b, v = Vector("a"), Vector("b"), Vector("v")
collect(Dot(Cross(a, v), Cross(b, v)), v)                                  # RecursionError
collect(Dot(a, Cross(Cross(b, VAdd(v, a)), Cross(a, VAdd(v, b)))), v)      # RecursionError
```

4 crashes in 931 multi-occurrence cases, none in 4069 single-occurrence cases. The pipeline
always calls `collect` with a target that occurs linearly, so this is reachable only through
the public API (`collect`, or `integrate_by_parts` on a user-built expression), but it is the
only non-documented exception the engine raised in 73k cases. A depth guard or a check that
`l.has(vec) and r.has(vec)` returns `expr` unchanged would close it.

### 6. Row vectors are encoded as "the operand is literally a `Transpose` node", and `expand`/`simplify` break that invariant on legal trees — crash (documented error on legal input)

Reproducers (8 distinct cases in the run; both patterns below):

```python
a, b, c, K, Q = Vector("a"), Vector("b"), Vector("c"), Matrix("K"), Matrix("Q")
col = Transpose(MVMul(Transpose(b), K))      # (bᵀK)ᵀ: a legal column vector (type VECTOR)
expand(VVMul(VAdd(a, col), Transpose(c)))    # ExpressionMismatchError: VVMul ... VECTOR and VECTOR
expand(VVMul(Transpose(c), VAdd(a, col)))    # same
expand(VVMul(VAdd(a, Transpose(Transpose(b))), Transpose(c)))   # same

row = MVMul(Transpose(a), K)                 # aᵀK: a legal row vector
simplify(MVMul(Transpose(Transpose(row)), Q))   # ExpressionMismatchError: MVMul ... VECTOR and MATRIX
```

Diagnosis. `_expand_vvmul` (`expand.py:88,95`) distributes over the `VAdd` and builds
`VVMul(Transpose(...), Transpose(c))`, which `VVMul.__init__`
(`multiplication.py:198`) rejects because *both* operands are `Transpose` instances even
though one is a column. In the second case `_eliminate`'s `Transpose` rule
(`simplify.py:179-180`) strips `(row)ᵀᵀ` to `row = MVMul(Transpose(a), K)`, and the enclosing
`MVMul` is rebuilt as `MVMul(row, Q)`, which `MVMul.__init__` (`multiplication.py:139`)
rejects because the left operand is no longer a `Transpose` node. Everything the mixins let a
user build here is type-correct; the transformations lose the row/column orientation.

### 7. `extract_linear_coeff` misses the nonlinearity `MVMul(f(q), g(q))` and returns a coefficient that still contains the query — wrong-value

`geomech/core/math/extract.py:442-450` (`_vector_from_vector`, `case MVMul()`): unlike every
other binary case in the module, there is no "query in both children → `NonLinearError`"
check; the vector side wins and the matrix side is returned verbatim.

```python
v, a = Vector("v"), Vector("a")
extract_linear_coeff(MVMul(Hat(v), v), v)                       # -> Hat(v)      (expr is v x v = 0)
extract_linear_coeff(MVMul(VVMul(v, Transpose(a)), v), v)       # -> v(a)'       (expr = v (a.v), quadratic)
extract_linear_coeff(MVMul(Hat(v), Cross(a, v)), v)             # -> Hat(v)Hat(a)
```

169 of 1792 quadratic vector-in-vector trees were accepted; all 169 returned coefficients
contain the query. Confirmation: adding the two-sided check to that one branch in-process
makes the same campaign raise `NonLinearError` for 1804/1804 quadratic trees with no new
false positives on the 5000 linear trees. In the pipeline the query is a variation vector,
which occurs linearly by construction, so this affects direct users of the API (and any future
code relying on `NonLinearError` as a guard).

### 8. `extract_linear_coeff` silently returns a zero coefficient for terms that depend on the query (two cells) — wrong-value (silent drop)

(a) `extract.py:650-651` (`_vector_from_matrix`, fallthrough `case _: return ZeroMatrix`),
unlike the other extractors which use `_zero_or_raise`:

```python
v, M, c = Vector("v"), Matrix("M"), Vector("c")
extract_linear_coeff(Hat(v), v)                                   # -> 0
extract_linear_coeff(MMMul(Hat(v), M), v)                         # -> 0M
extract_linear_coeff(Transpose(SMMul(VVMul(v, Transpose(c)), Scalar("s"))), v)   # -> 0
```

(b) `extract.py:689-693` (`_matrix_from_matrix`, `case SMMul()`): a scalar factor containing
the matrix query is ignored:

```python
Z, a, b, M = Matrix("Z"), Vector("a"), Vector("b"), Matrix("M")
extract_linear_coeff(SMMul(Z, Dot(a, MVMul(M, b))), M)            # -> 0
```

34 distinct cases (26 in cell (a), 8 in cell (b)). The coefficient of a vector query in a
matrix expression is a 3-tensor, so these cells cannot return a correct answer; they should
raise `AlgebraicError` as the `VVMul` and `SMMul` branches of `_vector_from_matrix` already do.

### 9. `extract_linear_coeff` raises `NotImplementedError` for a matrix query in a scalar or vector expression — cosmetic

`extract.py:148` (dispatch) and `extract.py:598` (`_vector_from_vector_svmul_scalar`
fallthrough). `NotImplementedError` is not one of the package's error types; the grid in the
module docstring marks these cells "n/a", so `AlgebraicError`/`UndefinedCaseError` would be
the documented way to refuse. 599 cases, all raised.

```python
extract_linear_coeff(Dot(Vector("a"), MVMul(Matrix("M"), Vector("b"))), Matrix("M"))  # NotImplementedError
```

### 10. `expand` never distributes `Vee` over a matrix sum — cosmetic

`expand.py:157-161` tests `isinstance(inner, VAdd)` on the (matrix) operand of `Vee`, which
can never be true; the `Hat` rule above it correctly tests `VAdd`.

```python
expand(Vee(MAdd(Hat(Vector("a")), Hat(Vector("b")))))   # -> Vee((Hat(a)+Hat(b)))   (Hat(VAdd) does split)
```

### 11. `vector_rules` does not recurse into the operands of a `Dot` — cosmetic (missed simplification)

`vector_rules.py:138-167`: every `Dot` branch returns either a constant or `expr` itself, so
S²/zero rules inside a `Dot` never fire, even under `full_simplify`:

```python
q = S2("q"); w = q.get_tangent_vector(); a, b, c = Vector("a"), Vector("b"), Vector("c")
full_simplify(Cross(q, Cross(w, q)))              # -> omega_q
full_simplify(Dot(a, Cross(q, Cross(w, q))))      # -> Dot(a,Cross(q,Cross(omega_q,q)))   unchanged
full_simplify(Dot(a, Cross(b, Cross(c, c))))      # -> unchanged (Cross(c,c) is 0)
```

### 12. Rules that fire purely on declared flags (adversarial mode) — informational, not bugs

With the declared properties violated numerically, these and only these rewrites changed the
value (115 `vector_rules`, 223 `simplify`, 228 `full_simplify` cases out of 1666 each, the
remaining 3/5 non-preserving cases being the exceptions of finding 6; the rest of the engine
is flag-independent). They are correct under their hypotheses and no rule
fired on an *unflagged* leaf in any of the 25,000 consistent-mode cases:

| rule | source | relies on |
|---|---|---|
| `Dot(u,u) → 1` | `vector_rules.py:157-160` | `UnitNorm` flag (also on a plain `Vector`) |
| `Dot(ω,q) → 0`, `Dot(ξ,q) → 0` | `vector_rules.py:163-164`, `_is_tangent_orthogonal` | `TS2.S2 is q` (identity, not name) |
| `Cross(q, Cross(a, q)) → a` and the four BAC-CAB variants | `_try_s2_cross_reduction` | `S2` leaf + `a` a `TS2` of the same point or a `Cross(·, q)` |
| `RᵀR → I`, `RRᵀ → I` | `simplify.py:132-137` | `OrthogonalMatrix` flag (`SO3` or plain `Matrix`) |
| `Jᵀ → J`, `(JM)ᵀ → MᵀJ` | `simplify.py:182-183` | `SymmetricMatrix` flag |
| `Hat(Vee(M)) → M` | `simplify.py:165-166` | nothing — see finding 3 |

Negative checks that passed in all consistent-mode cases: `Cross(q, Cross(a, q))` with a plain
`a`, with `TimeDerivative(ω)`, with `Variation(ω)`, or with the tangent vector of a different
`S2` point is left alone; `Cross(u, Cross(a, u))` with a unit-norm non-`S2` vector is left
alone; `MᵀM` for a plain matrix is left alone; `Dot(p, ξ_q)` is left alone.

### 13. Linear terms that `extract_linear_coeff` refuses with `AlgebraicError` — informational

Documented behaviour, listed so the counts above are interpretable (375 cases):

- `Transpose(Transpose(v))` and `Vee(Hat(v))` around the query (`_vector_from_vector` has no
  `Transpose`/`Vee` case). `full_simplify` removes both shapes before the pipeline extracts,
  so this is only visible through direct calls.
- `Transpose(M)`, `Hat(M·v)`, outer products for a matrix query (`_matrix_from_matrix` is
  explicitly partial).
- Outer products and `SMMul(M, Dot(q, w))` for a vector query in a matrix expression (the
  error messages say why; these have no `C·q` factorisation).

## No discrepancies found (with case counts)

- `expand`: 4999/4999 evaluable trees value-preserving (5000 consistent + 1666 adversarial).
- `pull`: 6666/6666.
- `vector_rules` in consistent mode: 5000/5000, including every S² motif under its
  hypotheses and the negative motifs listed under finding 12.
- `simplify` / `full_simplify` in consistent mode: all 47 mismatches are findings 3 and 4;
  zero elimination, numeric folding, like-term collection, `Hat(v)·w → v×w`,
  `(AB)v → A(Bv)`, `A(Bv) → (AB)v` when `AB` collapses, `Transpose` of
  `Hat`/`SMMul`/`MMMul`, `Vee(Hat(v)) → v` never changed a value. `full_simplify(max_iter=10)` reached a fixed point in all 4993
  non-raising cases (`simplify(result) == result`); no oscillation or growth was seen.
- `Variation(e)` vs `e.delta()` (product/chain rules through every node type, incl. the eager
  `S2.delta() = ξ×q`, `SO3.delta() = R·Hat(η)`, `Transpose`, `VVMul`, constants): 3400/3400.
- `TimeDerivative(e)` vs `e.t_diff()` (incl. `S2.t_diff() = ω×q`, `SO3.t_diff() = R·Hat(Ω)`,
  `d/dt` of `Variation`): 3394/3394 evaluable.
- `apply_manifold_rules`: the SO(3) substitution `δΩ → η̇ + Ω×η` is exact in all cases; the S²
  substitution is exact for the tangential part in 5000/5000 cases (finding 2 is the normal
  part only).
- `ibp._t_diff_product` (merged-power product rule, repeated factors, constants, `Dot`
  factors): 3401/3401.
- `extract_linear_coeff`, scalar query: 1407/1407 (scalar, vector and matrix expressions);
  vector query in scalar expression: 1240/1240 of the accepted trees; vector query in vector
  expression: 1013/1013 of the accepted trees (`Cross`, `MVMul`, `SVMul`, `VAdd`, `Hat`,
  `Transpose(Hat)`, `Transpose(SMMul)`, `Transpose(MMMul)`, outer-product matrices, `SVMul`
  with `Dot`/`VVMul`/`Mul`/`Add` scalars). No coefficient returned for a linear tree contained
  the query. No linear-by-construction tree was rejected as nonlinear.
- Nonlinearity detection: `Mul`, `Dot`, `VVMul`, `SVMul`, `SMMul`, `Cross`, `MMMul` raised
  `NonLinearError` for 4788/4788 quadratic trees; only `MVMul` (finding 7) misses.
- Exceptions: no `AttributeError`, `TypeError`, `KeyError` or timeout anywhere; the only
  non-documented exception is the `RecursionError` of finding 5, and the only documented
  errors raised by a transformation on a legal tree are the row-vector cases of finding 6.

## Harness limitations

- `TimeIntegral` and `integrate_by_parts` as a whole are not checked (not pointwise
  identities); the product rule it uses (`_t_diff_product`) and the `collect` step it calls
  are.
- Trees with nested `d/dt` depth > 4 are skipped (10 cases); second-order variations
  (`Variation` inside `Variation`) are not generated because the constructor rejects the
  direct form and the pipeline never builds the indirect form.
- Equality in the engine is string-based, so the fuzzer cannot distinguish two different leaves
  with the same name; all leaves used have distinct names.
