# Troubleshooting

Three kinds of trouble: the equation is wrong (or looks odd) and nothing was raised; a
typed error was raised; or you want to see what the engine did. The tools are listed at
the end.

## The equation looks wrong and nothing was raised

These can happen with the current version. Bugs that used to give a plausible wrong
equation silently, and are fixed, are listed separately under
[History](#history-silent-bugs-fixed-in-earlier-versions) so that an old notebook or an
old result can be recognised; they are not suspects for a wrong equation today.

| Symptom | Cause | Fix |
|---|---|---|
| A `δm`/`ṁ` term, or a `ValueError` "δS still contains the variation symbols '\delta{m}'" | `m` is not `Constant` | `Scalar("m", attr=["Constant"])`; attribute names are case-sensitive |
| `q × (ω̇ × q)` stays in an S² equation | `ω̇ ⊥ q` holds on the sphere but is never assumed: the tangency test sees `ω`, `ξ` and `Cross(·, q)` only, not a `d/dt` (`vector_rules` reaches every node, so `ω·q` and `q × (ω × q)` are reduced wherever they sit) | correct but unreduced; `to_standard_form` reads it as the projection `(I − qqᵀ) ω̇` |
| Partially simplified output, with a `RuntimeWarning` "no confirmed fixed point after N passes" | `full_simplify` stopped at `max_iter` | raise `max_iter`, or look at the expression in the warning |
| Two different inertias print the same and are merged | Two leaves with the same name are the same symbol (`str` is the identity): `J2 = Matrix("J1")` | give every symbol its own name; `print_tree` shows the names |
| On S², `to_standard_form`'s 3×3 solve leaves the sphere | `E` is meaningful on the tangent plane only | solve `(I − qqᵀ)(M ω̇ + f + G u) = 0` with `ω̇ ⊥ q` ([Conventions](conventions.md)) |
| `M` is negative | `E = −(Euler–Lagrange)`, so `M = −(physical inertia)` | by convention ([Conventions](conventions.md)) |
| A wrong sign, a term missing, an extra term, and no row here matches | an engine rule that rewrites a term wrongly | [the recipe below](#the-equation-is-wrong-and-no-row-matches) |

## The equation is wrong and no row matches

First confirm it numerically: `check_eom(L, dW, variables, eom)` compares the equation
with a finite-difference Euler–Lagrange residual computed without the library's rules
([Verifying an Equation Numerically](numeric-verification.md)); a relative error far
above `1e-5` means a wrong term. Then localise it in three steps.

1. **Bisect with `summary()`.** Keep only the terms of `L` (and `dW`) that the wrong
   term can come from, run `compute_eom(..., trace=tr)` with `tr = EOMTrace()`, and read
   `tr.summary()` against the formulas of [Pipeline](pipeline.md): the first stage that
   disagrees (a substitution that is not `ξ̇ − 2 ω × ξ`, a symbol that disappears, a
   node count that drops sharply) is the stage to suspect. `tr.ibp[t]["collect"]`,
   `["integrate"]` and `["simplify"]` split integration by parts per target.
2. **Run the engine fuzzer.** From the repository root, `python
   scripts/engine_fuzz.py --quick` applies every rewrite (expand, pull,
   vector_rules, simplify, collect, the manifold rules, the product rules) to thousands
   of random trees and compares input and output numerically. A rule that changes a
   value is reported as a `[wrong-value]` finding with a shrunk counterexample naming
   the transformation (`manifold_rules[full delta]`, `collect[target occurs once]`,
   `vector_rules: z_v Dot(q, ω_p) → z_v 0`); a correct engine has none. No hand
   derivation is needed.
3. **Rerun the suspect stage under `set_trace`.** `geomech.set_trace(True)`, then call
   the stage on the previous stage of the trace by hand, for example
   `full_simplify(tr["manifold_rules"])` or `collect(tr["ibp.prepare"], target)`. Each
   line is `rule name: before ==> after`; the module docstring lists every rule with
   its identity and the hypothesis it needs, so a line that is false on sight, or a
   rule that fired without its hypothesis, is the bug.

## History: silent bugs fixed in earlier versions

These used to give a plausible wrong equation (or a wrong print) without an error. They
are fixed: each now raises or derives correctly. They are listed so that an old notebook
or an old result can be recognised; for a wrong equation with the current version, use
[the recipe above](#the-equation-is-wrong-and-no-row-matches) instead.

| Old symptom | Old cause | Now |
|---|---|---|
| A `δm` term although `attr=["constant"]` was given | a typo in `attr` was ignored | `ValueError: attr: unknown attribute` |
| A variable had no equation, its terms were missing from the others | it was not listed in `SystemVariables` and its terms were dropped | `ValueError` naming `δy` |
| `−ḟ` where `f` was expected | `δW` written with the variation of the velocity, `Dot(x.t_diff().delta(), f)` | `ValueError` |
| `ω·q` or `q × (ω × q)` stayed with two `S2("q")` objects | the tangent rules compared objects | compared by name |
| Gyroscopic terms on S² at half strength (`ω × c` instead of `2 ω × c`) | S² rule `δω = ξ̇ − ω × ξ` | `ξ̇ − 2 ω × ξ`; re-derive (pendulum-type systems are unchanged) |
| A wrong SO(3) equation with `(M v)·(b × c)` terms; missing cross-inertia `M₁[Ω̇₂]` of a double rigid pendulum | `collect` used a false identity and could not separate `Dot(R₁a, R₂b)` | re-derive |
| Coriolis / rate-of-inertia terms (`−2 m l l̇ ω`, `J̇ Ω`) missing | integration by parts dropped `d/dt` of a non-constant scalar prefactor | re-derive |
| `Hat(Vee(I))` simplified to `I` | `Hat(Vee(X)) = X` was applied to any `X` | skew-symmetric `X` only |
| A coefficient of `0` for a term that holds the variable (`M·Hat(v)`) | `extract_linear_coeff` returned a silent zero | `AlgebraicError` "no left matrix coefficient" |
| One equation for two variables; duplicate rows in `linearize` | the same printed name listed twice | `ValueError` "is listed more than once" |
| The equation `(ω m + f)` | a tangent vector listed as a configuration variable | `TypeError` |
| An entry like `x + y` or `q.t_diff()` in `SystemVariables` gave an unrelated stray-symbol error | composite entries were accepted | `TypeError` "is a VAdd expression, not a single symbol" |
| `(0)x` in a result, `is_zero` false | `Number(0) * x` was not flagged zero | flagged zero |
| Two different products printed alike and compared equal (`Mul(ab, c)` and `Mul(a, bc)` both `abc`; `Variation(Scalar("x")) == Variation(Vector("x"))`) | products printed their factors with no separator, and `==` compared the printed form only | products print `*` (`ab*c`, `a*bc`), and a calculus node also compares its operand's class (phase 5b) |
| `x(5.55e-17)` | float roundoff of a cancelled like term | sums within `1e-12` of the largest coefficient are dropped ([Conventions](conventions.md)) |
| A substitution silently not applied under an unusual node | `apply_manifold_rules` returned the old node when it could not rebuild one | `TypeError` |
| Partially simplified output without a word | `full_simplify` stopped at `max_iter` silently | `RuntimeWarning` "no confirmed fixed point after N passes" |
| `x * True` behaved like `x * 1` | `bool` was taken as a number | `TypeError` |
| Two leaves built from one `attr` list shared flags | the list was aliased | copied |
| A parameter named like `\omega_{q}` merged with ω | reserved printed names were accepted | `ValueError` "has the printed name of the tangent vector" |
| `to_standard_form` gave `G[g] = −m e₃` and `linearize` a zero `B[δg]` column for an input `Scalar("g")` while `L` held the constant `g` | inputs were matched to leaves by printed name only, so a constant parameter counted as an input in one function and not in the other | `ValueError` "has the name of the constant 'g' in the equations" |
| `Transpose(v) * s` printed as `s v^T^T` and was treated as a row | it built the double transpose `Transpose(SVMul(Transpose(v), s))`, a column | the row `Transpose(SVMul(v, s))` |
| `ExpressionMismatchError` from `expand`/`simplify` on a legal tree with rows (`(a + (bᵀK)ᵀ) cᵀ`), or a row turned into a column | `Transpose(Transpose(x))` was kept, and orientation was read from "is a `Transpose` node", so the row encoding was lost under `expand` | `Transpose(Transpose(x))` is `x` at construction; orientation is `is_row_vector` |
| A correct equation printed as a wrong one: `a × b × c` (which product first?), `(b·c) a` printed as `b·c a`, a coefficient after its symbols | `to_latex` ignored precedence; `print_eom`/`render_eom` used `str` | `to_latex` parenthesises `Cross`/`VAdd`/`MVMul` operands and prints `(b·c)\,a`; the printers use `to_latex` |

## The `AlgebraicError`s about a calculus operator

All three end with the engine text "appears only inside the calculus operator … such a term
has no linear coefficient", but they have unrelated causes; the first sentence tells them
apart. `compute_eom` decides between the first and the third by the order of the
surviving derivative and by looking for an acceleration in `L`, before it says anything
about `collect`.

**1. A `d/dt` of a variation vector survived integration by parts** (from
`compute_eom`). The derivative is of first order and `L` holds no acceleration, but
`collect` could not bring its term into the shape `Dot(d/dt(η), f)` (for a scalar
variable: a direct factor `d/dt(δs) g` of a product), so the IBP pass did not integrate
it, and the coefficient of `η` cannot be read off:

```
AlgebraicError: compute_eom: while extracting the equation for '\eta_{R}': the time
derivative '\frac{d}{dt}(\eta_{R})' of the variation vector '\eta_{R}' survived integration
by parts in the term 'Dot(w,Vee(Hat(\frac{d}{dt}(\eta_{R}))*A))'
of δS; collect could not expose it as the left operand of a Dot (Dot(\frac{d}{dt}(\eta_{R}), f)),
so the term was not integrated by parts and has no coefficient. Inspect the term with
print_tree(term, style='indent'), or run compute_eom(..., trace=EOMTrace()) and look at the
stage 'ibp.collect[\frac{d}{dt}(\eta_{R})]'. (engine: extract_linear_coeff: '\eta_{R}' appears
only inside the calculus operator '\frac{d}{dt}(\eta_{R})' (a higher order derivative, a
variation or an integral of the query); such a term has no linear coefficient)
```

Fix: rewrite the term so the velocity sits in a vector operand (`Dot(w, Vee(Hat(Ω) A))`
above: express the `Vee` explicitly, or use an equivalent `Dot`/`Cross` form). A row
form `v.T() * A * w` in a kinetic energy is handled (`collect` turns an inner product
holding the target into a `Dot`); if one still survives, the message adds "write it as
Dot(v, A * w)".

**2. A jerk** (from `to_standard_form`): an acceleration appears only under a further
time derivative, so the equation is not second order:

```
AlgebraicError: to_standard_form: in equation '\delta{xs}', term
'\frac{d}{dt}(\frac{d}{dt}(\frac{d}{dt}(xs)))': the acceleration
'\frac{d}{dt}(\frac{d}{dt}(xs))' occurs only inside
'\frac{d}{dt}(\frac{d}{dt}(\frac{d}{dt}(xs)))', a third or higher time derivative of the
configuration (a jerk): the equation is not second order, so it has no standard form
M a + f + G u = 0 (engine: extract_linear_coeff: '\frac{d}{dt}(\frac{d}{dt}(xs))' appears
only inside the calculus operator '\frac{d}{dt}(\frac{d}{dt}(\frac{d}{dt}(xs)))' (a higher
order derivative, a variation or an integral of the query); such a term has no linear
coefficient)
```

Fix: a jerk comes from a virtual work that holds a third time derivative of the
configuration; `compute_eom` output is fine, only the standard form does not exist.

**3. A Lagrangian in accelerations** (from `compute_eom`): `L` holds `ẍ` (or `ω̇`, `Ω̇`),
so δL holds `d/dt(d/dt(δx))`. Integration by parts runs once per variation vector,
against its first derivative, so the second derivative survives. For `L = ½ ẍ·ẍ`:

```
AlgebraicError: compute_eom: while extracting the equation for '\delta{x}': the Lagrangian
depends on accelerations (L contains '\frac{d}{dt}(\frac{d}{dt}(x))'); higher-order
Lagrangians are not supported. Integration by parts runs once per variation vector, against
its first time derivative '\frac{d}{dt}(\delta{x})', so the term
'Dot(\frac{d}{dt}(\frac{d}{dt}(\delta{x})),\frac{d}{dt}(\frac{d}{dt}(x)))' of δS keeps
'\frac{d}{dt}(\frac{d}{dt}(\delta{x}))' and has no coefficient. Write L in the configuration
variables and their velocities only; run compute_eom(..., trace=EOMTrace()) and look at the
stage 'dL' to see where the acceleration enters. (engine: extract_linear_coeff: '\delta{x}'
appears only inside the calculus operator '\frac{d}{dt}(\frac{d}{dt}(\delta{x}))' (a higher
order derivative, a variation or an integral of the query); such a term has no linear
coefficient)
```

Fix: higher-order Lagrangians are not supported; write `L` in the configuration
variables and their velocities.

## Typed errors

Exact messages, produced by the code as of this version (the quoted symbols are the
ones of the example). Every message starts with the function or node that raised it.

### Building expressions

| Mistake | Message |
|---|---|
| `Dot(m, x)` | `ExpressionMismatchError: Dot cannot combine SCALAR 'm' and VECTOR 'x'` |
| `m + x` | `ExpressionMismatchError: Add cannot combine SCALAR 'm' and VECTOR 'x'` |
| `x - m` | `ExpressionMismatchError: Sub cannot combine VECTOR 'x' and SCALAR 'm'` |
| `x * J` | `ExpressionMismatchError: MVMul: a column vector times a matrix is not defined: 'x' * 'J'; write M * v for the column M v, or v.T() * M for the row vᵀ M` |
| `J * x.T()` | `ExpressionMismatchError: MVMul: a matrix times a row vector is not defined: 'J' * '(x)''; a matrix multiplies a column vector (M * v), and a row vector multiplies a matrix from the left (v.T() * M)` |
| `x.T() * y.T()` | `ExpressionMismatchError: VVMul: two row vectors cannot be multiplied: '(x)'' * '(y)''; write a.T() * b for the inner product aᵀb or a * b.T() for the outer product a bᵀ` |
| `x * y` | `ExpressionMismatchError: VVMul: two column vectors cannot be multiplied: 'x' * 'y'; write Dot(a, b) or a.T() * b for the inner product, a * b.T() for the outer product, or Cross(a, b)` |
| `Hat(m)` | `ExpressionMismatchError: Hat expects a VECTOR, got SCALAR 'm'` |
| `Vee(x)` | `ExpressionMismatchError: Vee expects a MATRIX, got VECTOR 'x'` |
| `Vector("a", size=(2,)) + x` | `SizeMismatchError: VAdd size mismatch: expected (2,), got (3,)` |
| `m * "x"` | `TypeError: cannot apply '*' to Scalar 'm' and str 'x'; operands must be geomech expressions or numbers (int, float)` |
| `m * True` | `TypeError: cannot apply '*' to Scalar 'm': got bool True; a bool is not a number here, write 1 or 0` |
| `x / m` | `TypeError: cannot divide Vector 'x' by the expression Scalar 'm'; only division by a number is supported. For a symbolic divisor declare its inverse as a symbol, e.g. Scalar('1/m', attr=['Constant'])` |
| `1 / m` | `TypeError: cannot divide 1 by the expression Scalar 'm'; only division by a number is supported. For a symbolic divisor declare its inverse as a symbol, e.g. Scalar('1/m', attr=['Constant'])` |
| `x / 0` | `ZeroDivisionError: cannot divide Vector 'x' by zero` |
| `Number([1])` | `TypeError: Number expects int, float, or str, got list` |
| `Scalar("m", attr=["constant"])` | `ValueError: attr: unknown attribute 'constant' for symbol 'm'; valid attributes are ['Constant', 'Identity', 'Manifold', 'Ones', 'OrthogonalMatrix', 'SkewSymmetricMatrix', 'SymmetricMatrix', 'TangentVector', 'UnitNorm', 'Zero']` |
| `Scalar("m", attr="Constant")` | `TypeError: attr must be a list of strings, got str 'Constant' for symbol 'm'; write attr=['Constant']` |
| `Transpose()` | `TypeError: Transpose: Transpose() needs an operand; write Transpose(expr)` |
| `Add()` | `TypeError: Add() needs at least one operand; use the zero leaf of the type for an empty sum` |
| `x.delta().delta()` | `ValueError: Cannot take variation of a variation: δ(δ(x)). Second-order variations vanish in Hamilton's principle.` |
| `xi.delta()` | `ValueError: Cannot take variation of variation vector \xi_{q}. The variation vector is already an infinitesimal perturbation.` |
| `TimeDerivative(q)` | `ValueError: TimeDerivative: TimeDerivative(q) is not represented in the tree for the manifold point 'q'; use q.t_diff(), which returns the kinematic expression ω × q` |
| `Variation(R)` | `ValueError: Variation: Variation(R) is not represented in the tree for the manifold point 'R'; use R.delta(), which returns the kinematic expression R·Hat(η)` |
| `m.t_integrate()`, `m` constant | `AlgebraicError: t_integrate: cannot integrate the constant scalar 'm' in time: ∫ m dt = m·t, and the tree has no time symbol. Only terms that are time derivatives (d/dt(x)) integrate to an expression` |

### `compute_eom`

| Mistake | Message |
|---|---|
| `SystemVariables(scalars=[x])`, `x` a Vector | `TypeError: compute_eom: SystemVariables.scalars entry 'x' has type ExprType.VECTOR, expected ExprType.SCALAR` |
| `SystemVariables(vectors=x)` | `TypeError: SystemVariables.vectors must be a list of configuration variables, got Vector 'x'; write vectors=[x]` |
| a `Constant` variable | `ValueError: compute_eom: SystemVariables.vectors entry 'e3' is constant: its variation vanishes, so it cannot be a configuration variable` |
| `ω` listed as a variable | `TypeError: compute_eom: SystemVariables.vectors entry '\omega_{q}' is a tangent or variation vector, not a configuration variable; list the manifold point 'q' instead` |
| `ẋ` listed | `TypeError: compute_eom: SystemVariables.vectors entry '\frac{d}{dt}(x)' is a TimeDerivative node, not a configuration variable; list the underlying variable 'x' instead` |
| `x + y` listed | `TypeError: compute_eom: SystemVariables.vectors entry '(x+y)' is a VAdd expression, not a single symbol; a configuration variable is one Scalar, Vector, S2 or SO3 leaf (list the symbols it is built from instead)` |
| a plain `Matrix` in `matrices` | `TypeError: compute_eom: SystemVariables.matrices entry 'A' is a general Matrix; only SO3 rotation matrices are supported as matrix configuration variables` |
| a name listed twice | `ValueError: compute_eom: SystemVariables: the name 'x' is listed more than once (twice in SystemVariables.vectors); configuration variables must be distinct symbols` |
| a variable named like `ω_q` | `ValueError: compute_eom: SystemVariables.vectors entry '\omega_{q}' has the printed name of the tangent vector of 'q'; configuration variables must be distinct symbols (rename it)` |
| `SystemVariables()` | `ValueError: compute_eom: no configuration variables listed in SystemVariables; a system with no degree of freedom has no equation of motion` — only for an `L` without variables (a constant): the check runs after the stray-symbol check, so `L = ½ m ẋ·ẋ` with `SystemVariables()` gives the stray-symbol message of the row "`y` used but not listed" for `'\delta{x}'` |
| a vector Lagrangian | `TypeError: compute_eom: lagrangian must be a scalar expression, got a VECTOR expression '\frac{d}{dt}(x)'` |
| `inf_work=True` | `TypeError: compute_eom: inf_work must be a scalar expression, None or 0, got bool True` |
| `y` used but not listed (or `mm` not `Constant`) | `ValueError: compute_eom: δS still contains the variation symbols '\delta{y}' of the symbols 'y', which are neither entries of SystemVariables nor declared in params; for each symbol, list it in SystemVariables, declare it in params=..., or mark it 'Constant'` |
| `Dot(x.t_diff().delta(), f)` in δW | `ValueError: compute_eom: inf_work contains '\frac{d}{dt}(\delta{x})', the time derivative of the variation vector '\delta{x}'; virtual work must be linear in the variations themselves, not their time derivatives (write Dot(x.delta(), f), not Dot(x.t_diff().delta(), f))` |
| a hand-built `Variation(x.t_diff())` | `ValueError: compute_eom: δS still contains the variation symbols '\delta{\frac{d}{dt}(x)}' of the symbols 'x', which are neither entries of SystemVariables nor declared in params; '\delta{\frac{d}{dt}(x)}' is the variation of a time derivative, which the engine never builds: write the variation of the configuration variable itself, x.delta() (the variation of the velocity is x.t_diff().delta(), built as d/dt of the variation, and is not allowed in inf_work)` |
| `params=l` | `TypeError: compute_eom: params must be a sequence of leaves, got Scalar 'l'; write params=[l]` |
| a `Constant` param | `ValueError: compute_eom: params entry 'm' is constant; a prescribed parameter varies in time (drop it from params or remove the 'Constant' attribute)` |
| a param that is a variable | `ValueError: compute_eom: params entry 'x' is also a configuration variable (SystemVariables.vectors); a symbol is either prescribed or a degree of freedom` |
| a param named like a constant of `L` (`params=[Scalar("l")]`, `l` `Constant` in `L`) | `ValueError: compute_eom: params entry 'l' has the name of the constant 'l' in the Lagrangian or the virtual work; a symbol is either a constant parameter or a prescribed parameter (rename one, or drop attr=['Constant'] if 'l' is a prescribed parameter)` |
| a velocity inside a matrix factor, `a·((Ω cᵀ) b)` | `AlgebraicError: compute_eom: while integrating by parts (stage ibp, targets \frac{d}{dt}(\eta_{R})): collect: the target '\frac{d}{dt}(\eta_{R})' occurs inside the matrix factor '\frac{d}{dt}(\eta_{R})*(c2)'' of the term 'Dot(\frac{d}{dt}(\eta_{R})*(c2)'*c3,c1)'; transposition exposes a target only in the vector operand ((M f(v))·w = f(v)·(Mᵀ w)), so this term cannot be integrated by parts` — write it as `(a·Ω)(c·b)` |
| an acceleration in `L` (`½ ẍ·ẍ`) | see [the AlgebraicErrors about a calculus operator](#the-algebraicerrors-about-a-calculus-operator) |
| a `d/dt(η)` survived IBP | see [the AlgebraicErrors about a calculus operator](#the-algebraicerrors-about-a-calculus-operator) |
| a term quadratic in a variation | `NonLinearError: compute_eom: while extracting the equation for '\delta{x}' from δS (extraction, after the stage 'expand'): _vector_from_scalar: nonlinear — '\delta{x}' appears in both sides of VVMul in '(\delta{x})'*\delta{x}'` |

### `to_standard_form`

`to_standard_form` and `linearize` validate `inputs` with the same function
(`geomech.dynamics.eom.check_inputs`): every input row of either table is raised by both,
with the caller's name in front.

| Mistake | Message |
|---|---|
| a jerk | see [the AlgebraicErrors about a calculus operator](#the-algebraicerrors-about-a-calculus-operator) |
| a term bilinear in two accelerations | `NonLinearError: to_standard_form: in equation '\delta{x1}', the coefficient of '\frac{d}{dt}(\frac{d}{dt}(x1))' still contains the acceleration '\frac{d}{dt}(\frac{d}{dt}(x2))' in '\frac{d}{dt}(\frac{d}{dt}(x2))*m'; the equation is not linear in the accelerations` |
| an input named like a constant of `L` (`inputs=[Scalar("g")]`, `g` `Constant` in `L`) | `ValueError: to_standard_form: input 'g' has the name of the constant 'g' in the equations; a symbol is either a constant parameter or an input (rename one, or drop attr=['Constant'] if 'g' is an input)` |
| a `Constant` input (`inputs=[g]`, `g` `Constant` in `L`) | `ValueError: to_standard_form: input 'g' is constant: its variation vanishes, so it cannot be an input` |
| an input listed twice (`inputs=[f, f]`) | `ValueError: to_standard_form: input 'f' is listed more than once (an entry that prints like 'f' is the same symbol)` |
| a set (or dict) of inputs, `inputs={f, M}` | `TypeError: to_standard_form: inputs must be an ordered sequence (list or tuple), got set; a set has no order, so the columns of B would depend on the hash seed` |
| an input printed like a velocity (`Vector(str(x.t_diff()))`) | `ValueError: to_standard_form: input '\frac{d}{dt}(x)' prints like the velocity '\frac{d}{dt}(x)' of the configuration variable 'x'; leaves are matched by printed name, so the two would be confused (rename the input)` |
| an input printed like an operation node of the equations (`Scalar(str(l.t_diff()))`, `l` a param) | `ValueError: to_standard_form: input '\frac{d}{dt}(l)' has the printed name of the TimeDerivative node '\frac{d}{dt}(l)' in the equations, not of a leaf; leaves are matched by printed name, so the two would be confused (rename the input)` |
| a param of `compute_eom` whose rate is in the equations given as an input (`inputs=[f, l]`) | `AlgebraicError: to_standard_form: input 'l' occurs in the equations inside the calculus operator '\frac{d}{dt}(l)'; an input enters the equations as itself, not through its rate (is 'l' a prescribed parameter of compute_eom(..., params=[...])? It can be an input only if its rate does not occur)` |
| an input of another size than its leaf (`Vector("f", size=(6,))`) | `ValueError: to_standard_form: input 'f' (Vector of size (6,)) has the name of the Vector of size (3,) 'f' in the equations, a different symbol; leaves are matched by printed name, so the two would be confused (rename one, or build the input with the class, attr and size it has in the equations)` |
| an input named like a named numeric leaf (`Scalar("k")`, `L` holds `Scalar("k", value=2.0)`) | `ValueError: to_standard_form: input 'k' has the name of the constant 'k' in the equations; a symbol is either a constant parameter or an input (rename one, or drop value= from the leaf 'k' if it is an input)` |
| not the `compute_eom` dict, `to_standard_form([], V)` | `TypeError: to_standard_form: eom must be the dict returned by compute_eom, got <class 'list'>` |

### `linearize`

| Mistake | Message |
|---|---|
| an input missing from `inputs` (or a symbol missing `Constant`) | `AlgebraicError: linearize: in row '\delta{x}', the term '\delta{f}' of the variation contains no perturbation symbol (configuration, velocity or input variation); the linearization is incomplete. It carries '\delta{f}' (of 'f'), which is neither a configuration, velocity nor input variation. If the symbol is an input, list it in `inputs`; if it is a constant, give it attr=['Constant']; if it is a prescribed parameter (compute_eom(..., params=[...])), linearize does not support prescribed parameters yet (pass params=[...] to linearize to have it say so).` |
| equations with a prescribed param, `linearize(..., params=[l])` | `AlgebraicError: linearize: in row '\xi_{q}', the term '\frac{d}{dt}(\omega_{q})*m*\delta{l}*l*(-2)' of the variation contains no perturbation symbol (configuration, velocity or input variation); the linearization is incomplete. It carries the variation of the prescribed parameter(s) 'l' declared in params: linearize does not support prescribed parameters yet (the variation of a prescribed parameter vanishes along the reference, but its rate does not in general).` |
| other `SystemVariables` than for `compute_eom` | `ValueError: linearize: eom keys ['\\delta{x}'] do not match the variation vectors ['\\xi_{q}'] of the given SystemVariables` |
| an input that is a variable | `ValueError: linearize: input 'x' is also a configuration variable` |
| an input named like a constant, `linearize(eom, V, [f, Scalar("g")])`, `g` `Constant` in `L` | `ValueError: linearize: input 'g' has the name of the constant 'g' in the equations; a symbol is either a constant parameter or an input (rename one, or drop attr=['Constant'] if 'g' is an input)` |
| an input of another class than its leaf (`Scalar("F")`, the equations hold the constant `Vector("F")`) | `ValueError: linearize: input 'F' (Scalar) has the name of the constant Vector 'F' in the equations, a different symbol; leaves are matched by printed name, so the two would be confused (rename one, or build the input with the class and attr it has in the equations)` |
| `ω` as an input | `ValueError: linearize: input '\omega_{q}' is a tangent or variation vector of a configuration variable; its variation is part of the perturbation state, not an input` |
| `ξ`, `Ω` or `η` as an input | `ValueError: linearize: input '\eta_{R}' is a tangent or variation vector of a configuration variable; its variation is part of the perturbation state, not an input` |
| a plain `Vector` printed like `Ω` as an input | `ValueError: linearize: the variation '\delta{\Omega_{R}}' of input '\Omega_{R}' is already a perturbation state (is '\Omega_{R}' a tangent vector of a configuration variable?)` |
| an input that is not a plain leaf (`ẋ`, `F + u`, a `Matrix`, an `S2`/`SO3` point) | `TypeError: linearize: input '\frac{d}{dt}(x)' must be a Scalar or Vector leaf, got TimeDerivative` |
| one input instead of a list, `linearize(eom, V, f)` | `TypeError: linearize: inputs must be a sequence of Scalar/Vector leaves, got Vector 'f'; write inputs=[f]` |
| not the `compute_eom` dict | `TypeError: linearize: eom must be the dict returned by compute_eom, got <class 'list'>` |

### Code generation and simulation

Generate the standard-form and the linearization module of a system from one
`StateLayout`, built with the `SystemVariables` and the `inputs` (same order) passed to
`to_standard_form` and `linearize` ([Code generation](codegen.md#the-generated-module)).
The first three rows are the checks that keep a controller from acting on the wrong
inputs or states.

| Mistake | Message |
|---|---|
| a linearization with other inputs than the layout (`linearize(eom, V, [T, f])`, `StateLayout(V, [f, T])`) | `ValueError: generate_module: the linearization's inputs ['\delta{T}', '\delta{f}'] are not the variations of the layout's inputs ['f', 'T'] in that order (expected ['\delta{f}', '\delta{T}']); build the StateLayout with the SystemVariables and the inputs (same order) passed to linearize and to_standard_form, and generate both modules from that one StateLayout` |
| a linearization with another variable order than the layout (`SystemVariables(vectors=[y, x])` for `linearize`, `[x, y]` for the layout) | `ValueError: generate_module: the linearization's configuration variables ['y', 'x'] (rows ['kin:y', 'kin:x', '\delta{x}', '\delta{y}']) are not the layout's variables ['x', 'y'] in that order (rows ['kin:x', 'kin:y', '\delta{x}', '\delta{y}']); build the StateLayout with the SystemVariables and the inputs (same order) passed to linearize and to_standard_form, and generate both modules from that one StateLayout` |
| `lqr(lin_module, ..., module=plant)` with modules of two layouts | `ValueError: lqr: module does not have the layout of lin_module: lin_module has the variables ['x'] and inputs ['f', 'T'], module has the variables ['x'] and inputs ['T', 'f']; build one StateLayout with the SystemVariables and the inputs (same order) passed to to_standard_form and linearize, and generate both modules from it` |
| a variation vector as a parameter (`collect_params` leaves them out) | `ValueError: StateLayout: parameter '\xi_{q}' is a tangent or variation vector of a configuration variable; it is part of the state or of the perturbation state, not a parameter (a variation in δW belongs to the derivation, not to the generated code)` |
| something without symbols in `collect_params` | `TypeError: collect_params: cannot find symbols in str 'dW'; pass the expressions of the system ([L, dW], or [L] plus the equations), the dict of compute_eom or of to_standard_form, a LinearizedSystem, or an iterable of these` |

### Engine functions called directly

| Call | Message |
|---|---|
| `extract_linear_coeff(a*a, a)` | `NonLinearError: _scalar_from_scalar: nonlinear — 'a' appears in both sides of Mul in 'a*a'` |
| `extract_linear_coeff(x * x.T(), x)` | `AlgebraicError: Cannot extract vector query 'x' from outer product 'x*(x)'': rank-1 matrix does not factor as C * q` |
| `extract_linear_coeff(Dot(ẋ, f), x)` | `AlgebraicError: extract_linear_coeff: 'x' appears only inside the calculus operator '\frac{d}{dt}(x)' (a higher order derivative, a variation or an integral of the query); such a term has no linear coefficient` |
| `extract_linear_coeff(Hat(x), x)` | `AlgebraicError: _vector_from_matrix: 'Hat(x)' is linear in 'x' but has no left matrix coefficient (the map is a 3-tensor, not C * query)` |
| `collect(x, x)` | `TypeError: collect: expected a scalar expression, got Vector 'x' of type VECTOR` |
| `full_simplify(expr, max_iter=1)` at the cap | `RuntimeWarning: full_simplify: no confirmed fixed point after 1 passes (the last pass still changed the expression, so no fixed point after 1 passes was confirmed by an unchanged pass); result may be unsimplified: <repr>` |

### `check_tree`

`check_tree(expr)` raises `InvariantError` naming the invariant, the path from the root
(`L`/`R` for binary nodes, the index for sums, `x` for a unary operand) and the subtree:

| Tree | Message |
|---|---|
| `Add([Add(m, g), k])` (a list is not flattened) | `InvariantError: check_tree: Add nested directly in a Add at root.0: Add '(m+g)'` |
| `Dot(x.T(), y)` | `InvariantError: check_tree: Dot operand is a row vector; Dot takes columns (use the column, or VVMul for rowᵀ·column) at root.L: Transpose '(x)''` |
| `x.T() + y` | `InvariantError: check_tree: VAdd mixes row and column vectors at root: VAdd '((x)'+y)'` |
| `TS2("w")` | `InvariantError: check_tree: TS2 leaf has no S2 point (.S2 is None) at root: TS2 'w'` |
| `Scalar("(-1)")` | `InvariantError: check_tree: Scalar is printed like a number but has no value; it would be treated as a symbol (use Number(...) for a numeric constant) at root: Scalar '(-1)'` |

The full list of invariants is in the docstring of `geomech/utils/invariants.py`.

## Tools

- **`compute_eom(..., trace=EOMTrace())`** keeps every stage (also after an exception:
  `tr.failed_at` then names the stage that raised, and `tr.summary()` starts with
  `FAILED at <stage>: <exception>`); `tr.summary()`, `tr.show("simplify")`. A trace is
  cleared when a pipeline starts, so reusing one object keeps only the last run.
  `linearize(..., trace=LinearizationTrace())` keeps δE, the term buckets, the raw and
  the final blocks per row. See [Pipeline](pipeline.md#seeing-the-stages).
- **`GEOMECH_DEBUG=1`** (environment variable, read at import): one `DEBUG` line per
  pipeline stage on the `geomech` logger (stage, node count, start of the `repr`) and
  `check_tree` on the inputs and outputs of `compute_eom`, `to_standard_form` and
  `linearize`. In a notebook, `logging.basicConfig(level=logging.DEBUG)` gives the same
  lines without the checks (call it before `set_trace(True)`; `geomech` adds its own
  `stderr` handler only when no handler is configured, so each line prints once).
- **`GEOMECH_TRACE=1`** or **`geomech.set_trace(True)`**: every firing of a named rewrite
  rule, `name: before ==> after`, on the `geomech.rules` logger. The rule names, with the
  hypothesis each needs, are in the docstrings of the transformation modules
  ([Pipeline](pipeline.md#which-rule-fired)).
- **`print_tree(expr)`** (default `style="indent"`): one node per line with class names,
  `L`/`R`/index labels, every flag of each leaf and, for a tangent or variation vector, its
  manifold point (`of=S2 'q'`); two different points with one name are numbered `#1`,
  `#2`. `repr(expr)` marks constants with `*` and the leaf kind (`v:x`, `M:J`, `S2:q`).
- **`check_tree(expr)`**: the structural invariants of a hand-built or mutated tree.
- **`check_eom(L, dW, variables, eom)`** / **`assert_eom`**: compare a derived equation
  with a finite-difference Euler–Lagrange residual computed without the library's rules
  ([Verifying an Equation Numerically](numeric-verification.md)). A wrong term shows up as
  a relative error far above `1e-5`.
- **`scripts/engine_fuzz.py --quick`**: checks every rewrite numerically
  on random trees; a `[wrong-value]` finding names the transformation and a minimal
  counterexample (step 2 of [the recipe](#the-equation-is-wrong-and-no-row-matches)).
- The modules are importable for breakpoints and patches under names distinct from their
  functions: `import geomech.core.transformations.simplification as m` gives the module
  (`expansion`, `pulling`, `vector_identities`, `simplification`,
  `geomech.core.math.collecting`), while `from geomech.core.transformations import
  simplify` gives the function.
