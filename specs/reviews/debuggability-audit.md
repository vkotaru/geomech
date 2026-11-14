# Debuggability audit: `geomech` (branch `feature/vbl`, HEAD `54d97b9`)

Question asked: in 6-12 months the author feeds a new mechanical system to `compute_eom` /
`to_standard_form` / `linearize`, gets either an exception or an equation that is wrong,
and must find the cause fast. How much does the code help?

Method: read-only. Every claim below was produced by running the public API in the
`venv/` (Python 3.12.13, suite at `1322 passed, 5 xfailed`); exact outputs are pasted.
Line numbers refer to the working tree at `54d97b9`. Nothing under `geomech/` or
`tests/` was modified.

Short verdict: the *validation* layer added with scalars and VBL (`_check_variables`,
`_check_inputs`, `_check_eom`, `_partition_terms`, `_zero_or_raise`) is in good shape and
its messages are the best in the package. The *engine* underneath it is still opaque: a
wrong equation gives the debugger no way to see intermediate stages without importing
non-exported functions, no way to know which rule rewrote what, no invariant checks, and
several places where user mistakes produce a plausible wrong equation silently (a typo in
an `attr` string, a configuration variable missing from `SystemVariables`, a second
`S2('q')` object). Those silent paths are the expensive ones; they are ranked first in
section (ii).

---

## (i) Scenario walkthroughs

### Scenario (a): the user gets an exception

Thirty-three deliberate mistakes through the public API. Exact messages, with the frames
inside `geomech/` from the traceback:

| # | Mistake | Exception and message | Names the symbol? | Names the stage? | Usable? |
|---|---|---|---|---|---|
| A1 | `SystemVariables(scalars=[x])` with `x` a Vector | `TypeError: SystemVariables.scalars entry 'x' has type ExprType.VECTOR, expected ExprType.SCALAR` (`eom.py:35`) | yes | implicitly | **good** |
| A2 | constant vector as config var | `ValueError: SystemVariables.vectors entry 'e3' is constant: its variation vanishes, so it cannot be a configuration variable` (`eom.py:40`) | yes | yes | **good** |
| A3 | `Dot(m, x)` | `ExpressionMismatchError: Dot cannot be performed between SCALAR and VECTOR` (`geometry.py:34`) | **no** | n/a | weak: in a 40-line Lagrangian you do not know *which* Dot |
| A4 | `m + x` | `ExpressionMismatchError: Add cannot be performed between SCALAR and VECTOR` (`addition.py:44`) | **no** | n/a | weak |
| A5 | `x * J` | `ExpressionMismatchError: MVMul cannot be performed between VECTOR and MATRIX` (`multiplication.py:139`) | **no** | n/a | weak, and "MVMul" is an internal class name the user never typed |
| A6 | `J * x.T()` | `ExpressionMismatchError: operation cannot be performed` (`expressions.py:425`, bare `raise ExpressionMismatchError`) | **no** | **no** | **useless** |
| A7 | `x.T() * x.T()` | `ExpressionMismatchError: VVMul cannot be performed between VECTOR and VECTOR` | no | n/a | misleading: both *are* vectors; the problem is the two transposes |
| A8 | `x - m` | `ExpressionMismatchError: Sub cannot be performed between VECTOR and SCALAR` | no | n/a | weak |
| A9 | `Hat(m)` | `ExpressionMismatchError: Hat cannot be performed between VECTOR and SCALAR` | no | n/a | wording wrong: reads as "between" two operands, but Hat is unary (it is "expected VECTOR, got SCALAR") |
| A10 | `Vector('a', size=(2,)) + x` | `SizeMismatchError: VAdd size mismatch: expected (2,), got (3,)` | no | n/a | ok |
| A11 | `x.delta().delta()` | `ValueError: Cannot take variation of a variation: δ(δ(x)). Second-order variations vanish in Hamilton's principle.` | yes | n/a | **good** |
| A12 | `xi.delta()` | `ValueError: Cannot take variation of variation vector \xi_{q}. ...` | yes | n/a | **good** |
| A13 | `extract_linear_coeff(a*a, a)` | `NonLinearError: _scalar_from_scalar: nonlinear — 'a' appears in both sides of Mul in 'aa'` | yes | function name | good (but `'aa'` is the ambiguous `str`, see F8) |
| A14 | extract `x` from `x xᵀ` | `AlgebraicError: Cannot extract vector query 'x' from outer product 'x(x)'': rank-1 matrix does not factor as C * q` | yes | no | good |
| A15 | extract `x` from `Dot(ẋ, f)` | `AlgebraicError: extract_linear_coeff: 'x' appears only inside the calculus operator '\frac{d}{dt}(x)' (a higher order derivative, a variation or an integral of the query); such a term has no linear coefficient` | yes | yes | good in isolation; **inside the pipeline it does not say which term of δS or which equation** (see C3 below) |
| A16 | `m.t_integrate()` with `m` constant | `NotImplementedError:` (empty message, `expressions.py:185`; same at `:298`, `:466`) | **no** | **no** | **useless** |
| A17 | `Number([1])` | `TypeError: Number expects int, float, or str, got list` | n/a | n/a | good |
| A18 | `collect(x, x)` on a vector | returns `x` unchanged, no error (`collect.py:36-37`) | - | - | silent passthrough on a type the function is documented not to accept |
| A19 | `linearize(eom, vars)` without `inputs` | `AlgebraicError: linearize: in row '\delta{x}', the term '\delta{f}' of the variation contains no perturbation symbol (configuration, velocity or input variation); the linearization is incomplete. Is an input missing from \`inputs\`, or a symbol missing the 'Constant' attribute?` | yes | yes (row + term) | **exemplary** |
| A20 | wrong `SystemVariables` to `linearize` | `ValueError: linearize: eom keys ['\\delta{x}'] do not match the variation vectors ['\\xi_{q}'] of the given SystemVariables` | yes | yes | **good** |
| A21 | input is a config var | `ValueError: linearize: input 'x' is also a configuration variable` | yes | yes | good |
| A22 | input is `ω` | `ValueError: linearize: input '\omega_{q}' is a tangent or variation vector of a configuration variable; its variation is part of the perturbation state, not an input` | yes | yes | good |
| A23 | `m = Scalar('m')` (forgot `Constant`) | **no error**; EOM gains `\frac{d}{dt}(x)\frac{d}{dt}(m)(-1)` | - | - | silent, plausible-looking wrong equation |
| A24 | `y` used in `L` but not in `SystemVariables` | **no error**; equation for `x` is correct, the `δy` terms are **dropped**, no equation for `y` | - | - | **silent loss of a degree of freedom** |
| A25 | `dW = Dot(ẋ.delta(), f)` instead of `x.delta()` | **no error**; produces `... + \frac{d}{dt}(f)(-1)` (IBP applied to the virtual work) | - | - | silent |
| A26 | `to_standard_form` on an equation with a jerk | `AlgebraicError: extract_linear_coeff: '\frac{d}{dt}(\frac{d}{dt}(x))' appears only inside the calculus operator '\frac{d}{dt}(\frac{d}{dt}(\frac{d}{dt}(x)))' ...` | yes | partially (`to_standard_form` is only in the traceback) | ok |
| A27 | `Matrix('J1')` twice (this exact bug is in `examples/rigid_pendulum.py:32`: `J2 = Matrix("J1", ...)`) | no error; the two objects compare equal, the double pendulum silently uses one inertia twice | - | - | silent |
| A28 | `SystemVariables([x])` positional | same as A1 (good) | yes | | good |
| A29 | `m * "x"` | `AttributeError: 'str' object has no attribute 'type'` (`expressions.py:143`) | no | no | raw Python error from inside the package |
| A30 | `m + 1` | `AttributeError: 'int' object has no attribute 'type'` (`addition.py:43`); note `m * 1` works, `m + 1` does not | no | no | inconsistent numeric wrapping |
| A31 | `1 + m` | `TypeError: unsupported operand type(s) for +: 'int' and 'Scalar'` (no `__radd__`) | | | |
| A32 | `-x` | `TypeError: bad operand type for unary -: 'Vector'` (no `__neg__`) | | | users will write `-m*g*...` |
| A33 | `x / m` | `TypeError: unsupported operand type(s) for /` (no `__truediv__`) | | | |

Conclusions for (a):

- Everything raised by `dynamics/` (A1, A2, A11, A12, A19-A22) says *which symbol* and
  *which rule* was violated. Keep this style; it is the model for the rest.
- Everything raised while *building the tree* (A3-A9) names only the two `ExprType`s. In a
  real Lagrangian that is the same message for every one of a dozen `Dot`s. The fix is
  one line per constructor: include `str(l)` and `str(r)` (F5).
- Two bare raises (A6, A16) carry no message at all.
- The silent cases (A23, A24, A25, A27) are the real risk: each yields a *plausible*
  equation. A24 in particular loses a DOF without any signal and is a ten-line check (F2).

### Scenario (b): the derived equation is wrong; can I look inside?

Reconstructing `compute_eom` for the rigid pendulum on SO(3) by hand (`README.md` example)
from the public API:

```python
from geomech import *                                                   # public
from geomech.core.transformations.manifold_rules import apply_manifold_rules  # NOT exported from geomech
from geomech.core.math.ibp import _apply_ibp                            # private (single IBP pass)

dL = L.delta()                       # 1
dS = dL + dW                         # 2
dS = apply_manifold_rules(dS, vars_) # 3   not in geomech.__init__, not in core.transformations.__init__
dS = full_simplify(dS)               # 4
dS = integrate_by_parts(dS, [eta.t_diff()])   # 5   (internally: expand, full_simplify, then per target collect -> _apply_ibp -> full_simplify)
dS = expand(dS)                      # 6
eq = extract_linear_coeff(dS, eta)   # 7
eq = full_simplify(eq)               # 8
assert str(eq) == str(compute_eom(L, dW, vars_)[str(eta)][1])   # True
```

Count: one function that is public-named but not exported (`apply_manifold_rules`), one
private function if you want to see the IBP pass *without* its surrounding
`expand`/`full_simplify`/`collect` (`_apply_ibp`), and the knowledge that
`integrate_by_parts` re-expands and re-simplifies internally. There is no `trace=`,
`verbose=`, logger, env var, or stages object. `compute_eom` (`eom.py:46-97`) computes
eight intermediate values and keeps none of them.

What each stage looks like (the 400-char `str` of stage 4 is the thing a debugger has to
read; the `repr` is somewhat better):

```
--- 3 apply_manifold_rules ---  [Add, 444 chars]
str : ((Dot((\frac{d}{dt}(\eta_{R})+Cross(\Omega_{R},\eta_{R})),J\Omega_{R})+Dot(\Omega_{R},J(\frac{d}{dt}(\eta_{R})+...
repr: ((<(d/dt(\eta_{R}) + cross(\Omega_{R}, \eta_{R})), J*\Omega_{R}> + <\Omega_{R}, J*(d/dt(\eta_{R}) + cross(\Omega_{R}, \eta_{R}))>)*(0.5) + ...
--- 5b _apply_ibp (private) ---  [Add, 461 chars]
repr: (<-\eta_{R}, J*d/dt(\Omega_{R})>*(0.5) + <cross(\Omega_{R}, \eta_{R}), J*\Omega_{R}>*(0.5) + <-\eta_{R}, J'*d/dt(\Omega_{R})>*(0.5) + ...
--- 8 full_simplify(eq) ---  [VAdd, 207 chars]
repr: (-J*d/dt(\Omega_{R}) -cross(\Omega_{R}, J*\Omega_{R}) + cross(\rho, cross(d/dt(\Omega_{R}), \rho))*-m + cross(\Omega_{R}, cross(\rho, cross(\Omega_{R}, \rho)))*-m + cross(\rho, R'*e_3)*-m*g + M)
```

Printing. There are three printers and none is right for debugging:

- `str()` (what `print_eom` and `__eq__` use) mixes LaTeX (`\frac{d}{dt}`, `\eta_{R}`)
  with Python (`Dot(`, `Cross(`) and juxtaposition for every product. It is not LaTeX and
  not readable. `print_eom` output for the rigid pendulum is
  `\int{\eta_{R} \cdot \Big((J\frac{d}{dt}(\Omega_{R})(-1.0)+Cross(\Omega_{R},J\Omega_{R})(-1)+...)\Big)}dt=0`.
- `repr()` (`printing.py:76-147`) is the most readable but shows **no class, no type, no
  flags**: `repr(Scalar('m')) == repr(Vector('m')) == 'm'`;
  `repr(Scalar('k', attr=['Constant'])) == repr(Scalar('k')) == 'k'`;
  `repr(Mul(Mul(a,-1),-1)) == '--a'`.
- `print_tree` (`printing.py:37-57`) defaults to `style="topdown"`, which for the
  six-term rigid pendulum equation is **230 columns wide** and wraps in any terminal. The
  `style="indent"` output is excellent (class names, `L`/`R` labels, leaf flags) but is
  not the default and is not mentioned in any doc:
  ```
  VAdd
    ├─ 0: SVMul
      ├─ L: MVMul
        ├─ L: Matrix('J', is_constant, is_symmetric)
        ├─ R: TimeDerivative
          ├─ TSO3('\Omega_{R}')
      ├─ R: Scalar('(-1.0)', value=-1.0, is_constant, is_numeric)
  ```
  Even there, `_leaf_detail` (`printing.py:338-348`) omits `is_orthogonal`,
  `is_identity` and `is_skew_symmetric`, so `SO3('R')` prints as `SO3('R', is_manifold)`
  with the orthogonality that drives the `RᵀR = I` rule invisible, and `TSO3('\Omega_{R}')`
  does not show *which* `SO3` it belongs to (the `is`-identity link that drives the S²
  tangent rules, see C6).

`str` collisions measured (`__eq__` is `type(a) is type(b) and str(a) == str(b)`,
`expressions.py:20-23`; leaf `has()` is `str(elem) == str(self)` with **no** type check,
`expressions.py:99-100`):

| pair | `str` equal | `==` | note |
|---|---|---|---|
| `Mul(a,b)` vs `Scalar('ab')` | yes (`'ab'`) | no (type) | but `Scalar('ab').has(Mul(a,b))` is **True** |
| `Scalar('m')` vs `Vector('m')` | yes | no | `Vector('m').has(Scalar('m'))` is **True**; same hash |
| `Zero` vs `ZeroMatrix` | yes (`'0'`) | no | `Evaluator` in tests has a special case for exactly this |
| `ZeroVector` vs `Vector('c', Constant).delta()` | **no** (`'0v'` vs `'0'`) | no | two different printed zero vectors in one tree |
| `SVMul(x,-1)` vs `SVMul(x, Scalar('(-1)'))` | yes | yes | a named scalar `(-1)` is indistinguishable from the number |
| `Scalar('k', Constant)` vs `Scalar('k')` | yes | **yes** | flags are not part of identity; whichever object a rule sees first wins |
| `Matrix('J1')` vs a second `Matrix('J1')` | yes | yes | `examples/rigid_pendulum.py:32` |

### Scenario (c): a transformation rule misfires

There is no rule-level trace and rules have no names. `vector_rules` (`vector_rules.py:122-201`)
and `_eliminate` (`simplify.py:69-222`) are single `match` statements with comments; the
only way to find which branch fired is a debugger breakpoint. Monkeypatching to trace is
possible but obstructed: `geomech/core/transformations/__init__.py:1-5` imports the
functions `expand`, `pull`, `simplify`, `vector_rules` into the package namespace, which
**shadows the submodules of the same name**, so
`import geomech.core.transformations.simplify as S` binds `S` to the *function*
(`type(T.simplify).__name__ == 'function'`); you need `sys.modules[...]` to reach the
module. A 15-line wrapper around `pull`/`vector_rules`/`_eliminate` then reports only
"module X changed the string", never the rule.

Three concrete misfires found while probing, each a mismatch between the module docstring
and the code:

C1. **The S² rules do not fire below `Dot`, `MVMul`, `Hat`, `MMMul`, `MAdd`, `Transpose`
or any calculus node.** `vector_rules` recurses only through `Add`, `Mul`, `VAdd`, `SVMul`
and `Cross` (`vector_rules.py:125-135, 185-190`); the `Dot` branch (`:138-167`) inspects
its children but never calls `vector_rules` on them. `full_simplify` iterates to a fixed
point but never reaches the nested node. The docstring (`vector_rules.py:1-17`) lists
"cross(q, cross(ω, q)) = ω" and "Linearity through Add and Mul" and does not state the
recursion limit.
```
full_simplify(Cross(q, Cross(om, q)))            -> \omega_{q}                       # top level: fires
full_simplify(Dot(a, Cross(q, Cross(om, q))))    -> Dot(a,Cross(q,Cross(\omega_{q},q)))   # unchanged
full_simplify(MVMul(J, Cross(q, Cross(om, q))))  -> JCross(q,Cross(\omega_{q},q))         # unchanged
full_simplify(Hat(Cross(q, Cross(om, q))))       -> Hat(Cross(q,Cross(\omega_{q},q)))     # unchanged
full_simplify(Hat(SVMul(a, Dot(q, q))))          -> Hat(aDot(q,q))                       # unit norm not applied
```
The pipeline mostly avoids this because `expand` lifts `Cross` to the top of each term,
but a user calling `full_simplify` on their own expression, or an EOM whose mass block
is `MVMul(J, Cross(q, Cross(ω̇, q)))`, gets an unreduced but correct expression and no
way to know why the rule "listed in the docstring" did not apply.

C2. **`extract_linear_coeff` returns a silent zero for two node shapes.** Every other
fallback goes through `_zero_or_raise` (`extract.py:81-90`), but `_vector_from_matrix`
ends in `case _: return ZeroMatrix` (`extract.py:650-651`) and `_matrix_from_matrix`'s
`SMMul` branch ignores a query inside the scalar factor (`extract.py:689-693`):
```
extract_linear_coeff(Hat(x), x)                                   -> 0     # Hat is linear in x; 0 is wrong
extract_linear_coeff(Transpose(Hat(x)), x)                        -> 0
extract_linear_coeff(SMMul(A, VVMul(Transpose(a), MVMul(J, w))), J) -> 0
```
These shapes are reachable from `to_standard_form` on a matrix-valued equation and from
any direct call; the module docstring says "Never returns zero for a term that contains
the query" (`extract.py:84`), which is false for these two branches.

C3. **When `collect` cannot expose the IBP target, the failure surfaces two stages later
in `extract_linear_coeff`, naming neither the term nor the stage.** With
`L = Dot(w, Vee(MMMul(Hat(Ω), A)))` (a shape `collect.py:40-149` has no rule for):
```
AlgebraicError: extract_linear_coeff: '\eta_{R}' appears only inside the calculus operator
'\frac{d}{dt}(\eta_{R})' (a higher order derivative, a variation or an integral of the
query); such a term has no linear coefficient
```
True but unhelpful: the real cause is that `integrate_by_parts` left a `d/dt(η)` behind
because `collect` returned the `Dot(w, Vee(...))` term unchanged (`collect.py:149`,
`case _: return expr`). The message should say "a d/dt(η) survived integration by parts
in term '...'; collect could not expose it" and name the term. The same message is what
`to_standard_form` produces for a jerk (A26), so the two unrelated causes are not
distinguishable from the text.

C6. **S² tangent rules key on object identity, everything else keys on the printed
name.** `_is_tangent_orthogonal` and `_is_tangent_to_s2` use `a.S2 is b`
(`vector_rules.py:41, 43, 109`), while `__eq__`, `manifold_rules` (`tangent_map[str(omega)]`,
`manifold_rules.py:42, 49`) and `_check_variables` use strings:
```
q1, q2 = S2('q'), S2('q')
q1 == q2                                          -> True
full_simplify(Dot(q1.get_tangent_vector(), q1))   -> 0
full_simplify(Dot(q1.get_tangent_vector(), q2))   -> Dot(\omega_{q},q)        # rule silently skipped
full_simplify(Cross(q2, Cross(q1.get_tangent_vector(), q2))) -> unchanged
```
A user who builds `q` in a helper and again in the main script gets an equation with
unreduced `ω·q` and `q×(ω×q)` terms and no message. Nothing in the docs says "construct
each manifold point exactly once".

---

## (ii) Findings, ranked by minutes of future debugging saved per hour of work now

Each item: where, what happens today, the concrete proposal (message text or API shape),
and an estimate of the work.

### F1. Unknown `attr` strings are silently ignored (10 min; saves hours per incident)

`ExprFlags.from_attr_list` (`types.py:27-54`) is a `match` with no `case _`. A typo
(`attr=['constant']`, `attr=['Const']`, `attr=['Symmetric']`) produces a symbol that is
quietly non-constant, so `δm` and `ṁ` terms appear (A23 above) and the equation is
plausibly wrong. Proposal:

```python
_VALID_ATTRS = {"Constant", "Zero", "Ones", "Identity", "UnitNorm", "SymmetricMatrix",
                "SkewSymmetricMatrix", "OrthogonalMatrix", "Manifold", "TangentVector"}
case _:
    raise ValueError(f"unknown attribute '{a}' for symbol; valid attributes are {sorted(_VALID_ATTRS)}")
```
`"TangentVector"` must be added to the set because `TS2`/`TSO3` append it
(`expressions.py:320, 341`) although no rule reads it.

### F2. A variation symbol left over after extraction is dropped without a word (30 min)

A24: a configuration variable used in `L` but absent from `SystemVariables` loses its
equation and its terms, silently. Likewise A25 (`δẋ` in `dW`) and a `Variation(y)` of any
non-listed symbol. `separate_variations` (`eom.py:100-109`) extracts per listed vector and
ignores everything else. Proposal, after the loop in `separate_variations` (or at the end
of `compute_eom`):

```python
stray = _variation_symbols(inf_action_integral) - {str(v) for v in variation_vectors}
if stray:
    raise ValueError(
        f"compute_eom: δS still contains the variation symbols {sorted(stray)} which are "
        f"not variations of any entry of SystemVariables; add the corresponding variable "
        f"or mark the symbol 'Constant'")
```
where `_variation_symbols` walks the tree collecting `Variation` leaves and `TS2`/`TSO3`
leaves whose name equals their parent's `variation_vector_name`. The VBL module already
has the exact analogue (`_partition_terms`, `vbl.py:465-470`) and its message is the one
to copy.

### F3. A stages object for the pipeline, and export `apply_manifold_rules` (1-2 h)

`compute_eom` (`eom.py:46-97`) should keep what it computes. Proposed shape, zero cost on
the happy path:

```python
@dataclass
class EOMTrace:
    dL: Expr; dS: Expr; after_manifold_rules: Expr; after_simplify: Expr
    variation_vectors: list; after_ibp: Expr; after_expand: Expr
    raw_equations: dict      # before the final full_simplify, per variation vector
    equations: dict

def compute_eom(lagrangian, inf_work, variables, *, trace: EOMTrace | None = None): ...
# usage:  tr = EOMTrace(); eqs = compute_eom(L, dW, vars_, trace=tr); print_tree(tr.after_ibp, style="indent")
```
Alternatively `compute_eom_trace(...) -> EOMTrace` with `compute_eom` returning
`.equations`. Either way, add `apply_manifold_rules` to `geomech/__init__.py` and to
`core/transformations/__init__.py`; it is the one stage a user cannot currently reproduce
without a deep import. Do the same for `linearize` (`dE` per row, the buckets from
`_partition_terms`, the raw coefficient before `_negate`).

### F4. Re-raise engine errors with pipeline context (30 min)

Wrap the two call sites where engine errors escape with no stage/variable/term:

- `separate_variations` (`eom.py:106-109`): catch `AlgebraicError | NonLinearError`,
  re-raise `type(err)(f"compute_eom: while extracting the equation for '{vec}' from δS: {err}. "
  f"A 'd/dt({vec})' that survived means integration by parts could not reach that term; "
  f"inspect the term with print_tree(..., style='indent')") from err`.
- `_extract_single` (`standard_form.py:122-177`): same with
  `f"to_standard_form: equation '{key}', term '{term}': {err}"` (the `key` is available in
  `to_standard_form`, pass it down).
- `_has_algebraic` (`extract.py:67-74`) should include the enclosing term: callers that
  own a term (`_partition_terms`, `_extract_single`, every `Add`/`VAdd` branch in
  `extract.py`) can pass it; the message today shows only the calculus node.
- `_check_coefficient_free` (`vbl.py:504-505`) catches `AlgebraicError` and discards its
  message (`leftover = True`); include `err` in the re-raise like `_partition_terms`
  (`vbl.py:455-460`) already does.

### F5. Constructor errors must name the operands (20 min)

`ExpressionMismatchError.__init__` (`errors.py:4-9`) takes only an operation name and two
types. Add optional `lhs`/`rhs` expressions and print them:

```
Dot cannot be performed between SCALAR 'm' and VECTOR 'x'
Hat expects a VECTOR, got SCALAR 'm'
Matrix * Transpose(Vector) is not defined: 'J' * '(x)''; use J * x or Transpose(x) * J
```
Specific sites: `geometry.py:34, 59, 80, 101, 165`; `multiplication.py:63, 87, 111, 139,
141, 163, 198, 200`; `addition.py:38, 44`; `mixins.py:51-59, 94`; the **bare**
`raise ExpressionMismatchError` at `expressions.py:425` (A6) and the three bare
`raise NotImplementedError` for `t_integrate` of constants at `expressions.py:185, 298, 466`
(A16; should be `AlgebraicError("cannot integrate the constant 'm'")`). `UndefinedCaseError`
(`errors.py:46-53`) is unreachable (the type is always one of three) and carries no
context; delete it or give it a message. Also wrap ints in `__add__`/`__sub__` the way
`__mul__` does (`m * 1` works, `m + 1` raises `AttributeError`, A30), add `__radd__`,
`__rmul__`, `__neg__`.

### F6. Make the debug printers show what the engine keys on (30 min)

- `print_tree` default `style="indent"` (`printing.py:37`); the top-down renderer is
  unusable past three terms.
- `_leaf_detail` (`printing.py:338-348`): add `is_orthogonal`, `is_identity`,
  `is_skew_symmetric`, and for `TS2`/`TSO3` the parent: `TSO3('\Omega_{R}', of=SO3 'R' @0x..)`
  so that C6 (two parents with the same name) becomes visible.
- `repr_str` leaf branch (`printing.py:85-93`): prefix non-scalars the way `_compact_leaf`
  already does (`v:x`, `M:J`, `S2:q`, `SO3:R`) and mark constants with `*`. Then
  `repr(Scalar('m')) != repr(Vector('m'))` and `print(eq)` in a notebook shows which
  symbols the engine treats as constant.
- `print_eom`/`eom_to_latex` (`printing.py:19-29, 624-640`) use `str()`, which is neither
  LaTeX nor readable; switch them to `to_latex` (already correct and used by
  `display_eom`).

### F7. Replace the two silent zeros in `extract.py` with `_zero_or_raise` (10 min)

`extract.py:650-651` (`_vector_from_matrix`, default case) and `extract.py:689-693`
(`_matrix_from_matrix`, `SMMul` with the query in the scalar). C2 above shows
`extract_linear_coeff(Hat(x), x) -> 0`. Both should raise
`AlgebraicError(f"_vector_from_matrix: cannot factor '{expr}' as C * '{query}' (Hat/Transpose of the query is linear but not a left matrix product)")`.
`manifold_rules._rebuild` (`manifold_rules.py:126`) has the same pattern: an unknown node
class returns the **original** node with the **old** children, silently discarding the
substitution; make it `raise TypeError(f"apply_manifold_rules: cannot rebuild {type(expr).__name__}")`.

### F8. An invariant checker and three constructor guards (1 h)

Nothing checks tree invariants after construction. Today you can build
`TimeDerivative(S2('q'))`, `Variation(S2('q'))`, `Variation(SO3('R'))` (the contributor notes'
invariant "the tree never contains `TimeDerivative(q)` for a manifold point" is a
convention, not a check; `calculus.py:69-70` has no guard) and `Transpose()` with no
child (`geometry.py:118-122`, `type` is then `None`). Proposal:

```python
# geomech/utils/invariants.py
def check_tree(expr) -> None:
    """Raise InvariantError naming the offending subtree.  Checks: operand types of every node
    match its class; Add/VAdd/MAdd children are not the same Add class; VVMul has exactly one
    Transpose child; Transpose has a child; no TimeDerivative/Variation whose child is S2/SO3;
    TS2/TSO3 parent is not None; numeric Scalars have value is not None iff is_numeric."""
```
Call it at the top of `compute_eom`, `to_standard_form`, `linearize` under
`if os.environ.get("GEOMECH_DEBUG")` (see F11), and unconditionally in tests. Also add the
`S2`/`SO3` guard in `TimeDerivative.__init__`/`Variation.__init__` (same style as the
existing `Variation(Variation)` guard, `calculus.py:23-42`), and reject `Transpose()`.
`_check_variables` (`eom.py:21-43`) should also reject duplicates
(`SystemVariables(vectors=[x, x])` yields one equation silently, C4g) and should check
that each listed variable actually occurs in `L` or `dW` (an all-zero equation today).

### F9. Name the rules and add a trace hook (2 h)

In `vector_rules.py:122-201`, `simplify.py:69-222`, `pull.py:13-109`, `collect.py:40-149`
and `expand.py:98-174`, each `case` that rewrites should go through one helper:

```python
# geomech/utils/trace.py
_TRACE = os.environ.get("GEOMECH_TRACE") == "1"   # or geomech.set_trace(True)
def rule(name, before, after):
    if _TRACE and after is not before:
        logger.debug("%-28s %s  ==>  %s", name, repr(before), repr(after))
    return after
# vector_rules.py
case Dot() if ...: return rule("dot.tangent_orthogonal", expr, _ZERO())
```
Rule names then double as documentation and as the key for a per-rule test. The names
should be listed in the module docstring next to their hypothesis (today the docstring of
`vector_rules.py:1-17` lists identities without the conditions the code checks: the S²
rules require `TS2` with the *same parent object* or a `Cross(*, q)` child, C1/C6; the
unit-norm rule requires both sides `is_unit_norm` **and** `==`).

### F10. Documentation that would have answered scenarios (b) and (c) (2-3 h)

Missing pages (`mkdocs.yml` nav: Home, Getting Started, Examples, VBL, API):

- **Pipeline, stage by stage**: the real order is δ → `apply_manifold_rules` →
  `full_simplify` → `integrate_by_parts` (which itself does expand → full_simplify →
  per target: `collect` → IBP → full_simplify) → `expand` → `extract_linear_coeff` →
  `full_simplify`. Three docs give three different orders: `eom.py:1-8` ("variation →
  simplify → IBP → expand → extract", no manifold rules), `README.md` "How it works"
  ("variation → expand → integrate by parts → simplify → extract"), `docs/getting-started.md`
  ("Variation, Expand, Integrate by parts, Simplify, Extract"). Pick the code's order and
  put it in one place with the F3 trace object as the way to see each stage.
- **Troubleshooting**: the silent cases A23/A24/A25/A27/C6 with their symptom
  ("a `ṁ` term appears" → missing `Constant`; "`ω·q` does not vanish" → two `S2('q')`
  objects; "a variable has no equation" → not in `SystemVariables`), the two
  `AlgebraicError` texts with their causes (jerk vs. IBP could not reach the term), and
  `print_tree(..., style="indent")`.
- **Glossary and conventions**: `E = 0` is the *negative* Euler-Lagrange expression;
  `to_standard_form` gives `M = −(inertia)`; VBL's `δE = −M δv̇ + A δx + B δu`
  (`vbl.py:48-66` has it, `docs/vbl.md` "Sign convention" has it, nothing EOM-side does);
  `SVMul(vec, scalar)` is normalized vector-first; `(-1)` is a `Scalar` named `"(-1)"`;
  `Zero`/`ZeroVector`/`ZeroMatrix` print `0`/`0v`/`0`; `str()` is the identity.
- Fix stale text: `variables.py:12` says scalars are "not yet used in EOM pipeline";
  `docs/manifold-implementation-notes.md:3-9` says the S² variation formula "is not
  confirmed correct" (it is pinned by `tests/test_manifold_pipelines.py` and re-derived in
  `specs/reviews/vbl-review.md`); that page is not in the nav and should be rewritten as
  the pipeline page or removed. `extract.py:11` documents the signature as
  `extract_linear_coeff(query, expr)`; it is `(expr, query)`. `extract.py:13` "Assumes
  query appears linearly" contradicts the `NonLinearError` the code raises.

### F11. Logging and a `GEOMECH_DEBUG` switch (30 min)

There is no `logging` use anywhere in `geomech/`. A module logger
(`logging.getLogger("geomech")`) with one `debug` line per pipeline stage in
`compute_eom`, `integrate_by_parts` (per target) and `linearize` (per row), printing
`repr` and the node count, costs nothing when disabled and is the cheapest way to get
scenario (b) without the API change in F3. Wire `GEOMECH_DEBUG=1` to set the level and to
enable `check_tree` (F8) and the rule trace (F9).

### F12. Warn when `full_simplify` hits its cap (5 min)

`simplify.py:56-60` stops after `max_iter=10` passes with no signal. One `compute_eom` of
the rigid pendulum spends 9 `simplify` passes across its `full_simplify` calls; whether
any single call is near the cap is invisible from outside. Add
`warnings.warn(f"full_simplify: no fixed point after {max_iter} passes; result may be unsimplified: {expr!r}")`.

### F13. Make `Evaluator` importable and add an EOM-level numeric oracle (1-2 h)

`tests/numeric_helpers.py` (`Evaluator`, `finite_difference_check`) cannot be imported
from outside the repo root (`ModuleNotFoundError: No module named 'tests'` from any other
cwd; `pyproject.toml` packages only `geomech*`). The reviewer of the double pendulum
re-implemented an evaluator in `specs/reviews/scripts/double_rigid_pendulum_numeric_check.py`
because of this. `finite_difference_check` verifies a *linearization*; there is no entry
point that verifies an *EOM*. Proposal: move `Evaluator`, `hat`, `vee`, `rand_rot`,
`expm_so3`, `sample_state` to `geomech/numeric.py` (numpy is already the runtime
dependency) and add

```python
def evaluate(expr, env: dict[str, np.ndarray], seed=0) -> np.ndarray
def check_eom(eom, variables, reference: Callable[[dict], dict], trials=5, rtol=1e-8) -> float
    """Evaluate each derived equation at random states and compare with a user-supplied
    closed form (e.g. the textbook equation written with numpy); returns the worst error."""
```
Keep `tests/numeric_helpers.py` as a thin re-export so the suite is unchanged.

### F14. Un-shadow the transformation submodules (15 min)

`core/transformations/__init__.py:1-5` binds `expand`, `pull`, `simplify`,
`vector_rules` to the functions, hiding the modules. Either rename the modules
(`_expand.py`, ...) or keep the modules importable by not re-exporting the function under
the module's own name (`from .simplify import simplify as simplify` still shadows;
renaming is the only clean fix). Without this, F9's trace and any monkeypatch need
`sys.modules`.

### F15. Duplicates, naming, dead code, magic numbers (1-2 h, maintenance value)

- Duplicated helpers: `_flatten_addends` (`standard_form.py:97`, `vbl.py:424`);
  `_flatten_mul` (`ibp.py:78`, `simplify.py:468`); `_zero_for` (`standard_form.py:104`,
  `simplify.py:248`); `_combine` (`standard_form.py:111`, `vbl.py:430`); the constant −1
  is built five ways (`manifold_rules.py:26 _NEG1`, `vector_rules.py:61 _neg1`,
  `simplify.py:186` inline, `_wrap_numeric(-1)` through `* -1`, `vbl.py:513 _negate`);
  `_ZERO()`/`_ONE()` lambdas (`vector_rules.py:118-119`) duplicate the `Zero`/`One`
  singletons (`expressions.py:523-524`) for no reason since identity is by string.
  One `geomech/core/base/constants.py` with `NEG_ONE`, `HALF`, `negate(expr)` and
  `zero_like(expr)` removes all of them.
- Private cross-module imports: `vbl.py:98` (`_has_algebraic`, `_same`),
  `standard_form.py:13` (`_has_algebraic`), `vbl.py:105` (`_check_variables`),
  `geometry.py:154` (`_wrap_numeric`). Promote `has_algebraic`, `same_node`,
  `check_variables`, `wrap_numeric` to public names in their modules.
- Two sources of truth for the manifold variation formulas: `manifold_rules.py:57-82`
  and the never-called `TS2.delta(substitute=True)` / `TSO3.delta(substitute=True)`
  (`expressions.py:322-331, 343-351`). Delete the latter or make the former call them.
- Dead or inert: `"TangentVector"` attr appended but never parsed; `SkewSymmMatrix`
  (`expressions.py:482`) used only by the printer; `is_skew_symmetric` never consulted
  (the transpose rule tests `isinstance(inner, Hat)`, `simplify.py:185`);
  `has_nested_add`/`is_leaf` exported but unused (`expand.py:4` says so);
  `Matrix.__init__` lines 477-478 are a no-op; `UndefinedCaseError` unreachable.
- Same concept, several names: "variation vector" (`eom.py`), "perturbation"
  (`vbl.py`, 29 uses), `var_vec`/`variation_vec`/`variation_vectors`; `inf_work` /
  `dW` / `deltaW`; `targets` (ibp) / `variation_vector_dots` (eom) / `dvel_perts` /
  `dstates` (vbl). Pick "variation vector" for the EOM side, "perturbation" for VBL, and
  say so in the glossary (F10).
- Long functions with nested `match`: `_eliminate` 154 lines (`simplify.py:69`),
  `linearize` 128 (`vbl.py:181`), `to_latex` 119 (`printing.py:420`), `_collect_dot` 110
  (`collect.py:40`), `pull` 97 (7 nested `match`, `pull.py:13`), `vector_rules` 80
  (`vector_rules.py:122`). Splitting `_eliminate` per node class (`_elim_mvmul`, ...)
  is also what makes F9's rule names natural.
- Magic numbers: `max_iter=10` (`simplify.py:45`), `len(inner) <= 10`
  (`printing.py:140`), `gap=3` (`printing.py:214`), `[-2000:]` (`printing.py:672`).
- `ibp._remove_factor` tests `str(left_removed) == "1"` (`ibp.py:38, 44`) instead of
  `is_one`/`is One`; a user `Scalar('1')` would be treated as the number one.
- `examples/rigid_pendulum.py:32`: `J2 = Matrix("J1", ...)` is a bug (two inertias, one
  name) and a good troubleshooting example once F6 makes it visible.

### F16. Type checking (reported, not fixed)

`pyright geomech` (pyright 1.1.414, default settings): **102 errors, 0 warnings** —
`reportIncompatibleMethodOverride` 47, `reportAttributeAccessIssue` 33,
`reportAssignmentType` 7, `reportOptionalMemberAccess` 7, `reportCallIssue` 7,
`reportArgumentType` 1; by file `mixins.py` 31, `expressions.py` 24, `geometry.py` 13,
`multiplication.py` 12, `addition.py` 6, `calculus.py` 6, `simplify.py` 5, `vbl.py` 3,
`standard_form.py` 2. `mypy geomech --ignore-missing-imports`: **9 errors** (8
`[assignment]`: dataclass fields declared `str`/`ExprFlags` defaulting to `None`,
`expressions.py:155, 158, 262, 266, 435, 439`, `standard_form.py:38`,
`manifold_rules.py:49`; 1 `[attr-defined]`, `addition.py:41`). Only 18 of 181 functions
carry any annotation. Root causes are structural and cheap to fix: the mixins access
`self.nodes`/`self.type` without declaring them (a `Protocol` or class-level annotations
on `_BaseMixin` removes all 31 in `mixins.py`); `Scalar.delta` returning `Zero` where
`Expr.delta` is declared to return `Variation` (annotate `-> Expr`); `name: str = None`
should be `str | None` or non-optional. Adding `mypy --ignore-missing-imports geomech`
to CI after those fixes would catch the `type`-less `Transpose()` and the `str`-on-int
cases (A29/A30) at edit time.

### F17. `collect` and the S² identity link (20 min)

`vector_rules.py:41, 43, 109`: compare parents with `==` (string identity, consistent with
the rest of the engine) or, better, have `S2.__init__`/`SO3.__init__` intern instances by
name so a second `S2('q')` returns the first (then C6 cannot happen). Document whichever is
chosen. `collect` (`collect.py:21-37`) should raise `TypeError` for a non-scalar top-level
argument instead of returning it (A18).

---

## (iii) Already good: keep these

- **Typed, contextual engine errors**: `AlgebraicError` and `NonLinearError` carry the
  query, the node type and the offending expression (`extract.py:67-74, 87-89, 108-110,
  630-634, 643-647`); `_zero_or_raise` (`extract.py:81-90`) is the right pattern and only
  two branches miss it (F7). No bare `except`, no `assert` for validation, no swallowed
  exceptions except the one message loss in `vbl.py:504-505`.
- **The `dynamics/` validation layer**: `_check_variables` (`eom.py:21-43`),
  `_check_inputs` / `_check_eom` / `_perturbation_variables` (`vbl.py:316-416`) and
  `_partition_terms` (`vbl.py:438-479`) name the symbol, the row, the term and the likely
  user mistake ("Is an input missing from `inputs`, or a symbol missing the 'Constant'
  attribute?"). A19-A22 are the messages every other raise should look like.
- **The completeness rule in VBL** (every term must carry exactly one perturbation symbol,
  `vbl.py:68-76, 438-479`): it converts silent zeros into errors. F2 asks for the same
  rule on the EOM side.
- **`Variation` guards** (`calculus.py:23-42`): `δ(δx)` and `δ(ξ)` are rejected with a
  sentence of physics.
- **`print_tree(style="indent")`** (`printing.py:304-383`): class names, `L`/`R`/index
  labels, leaf flags and values. It only needs to be the default and to show three more
  flags.
- **Tests as documentation**: system tests are classes named after the system
  (`TestRigidPendulumPivotPipeline`, `TestSphericalPendulumS2`, `TestQuadrotorWithLoad`),
  pin the **full equation string** and cite Lee et al. equation numbers in docstrings
  (`tests/test_manifold_pipelines.py`, `tests/vbl_systems.py`, `tests/test_vbl.py`);
  `tests/test_vbl_validation.py` tests the *text* of error messages. `grep -l "class
  Test<System>" tests/` finds the test for a system in one step. The strict `xfail`
  convention for known gaps (`tests/test_scalar_adversarial.py`) is worth keeping.
- **`finite_difference_check`** (`tests/numeric_helpers.py:334-363`) and the independent
  evaluator in `tests/test_vbl_adversarial.py` are a real numeric oracle for VBL; F13 only
  asks to make them reachable and to add an EOM-level counterpart.
- **Module docstrings with derivations**: `vbl.py:1-77` (perturbation table, kinematic
  rows, sign convention, worked free-rigid-body example, completeness rule),
  `manifold_rules.py:1-15`, `expand.py:1-22`, `extract.py:1-14` (dispatch grid),
  `ibp.py:1-9`. The gaps are the hypotheses the code checks (F9, F10), not the math.
- **An earlier code review, `CHANGELOG.md`, `specs/reviews/*`**: the known-limitation list
  (string equality, commutative rearrangements) and the spec/adversarial-review/final-review
  cycle are exactly what a returning author needs to re-orient. Do not let them go stale.
- **Immutability and `match`/`case` rewrites**: every transformation is a pure function
  on an immutable tree, so inserting a trace or a stages object (F3, F9, F11) is a local
  change with no risk to results.
