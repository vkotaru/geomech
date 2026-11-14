# Changelog

All notable changes to `geomech` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Debugging aids (phase 4; `docs/pipeline.md`, `docs/troubleshooting.md`):
  - `compute_eom(..., trace=EOMTrace())` (keyword only) keeps every stage — δL, δW, δS,
    after the manifold rules, after `full_simplify`, per integration-by-parts target the
    result of `collect`, of the IBP pass and of `full_simplify`, after `expand`, and the
    raw and final equation per variation vector — filled as the pipeline runs, so it also
    holds the stages before an exception. `linearize(..., trace=LinearizationTrace())`
    keeps per row δE, the term buckets, the raw coefficients (before `full_simplify` and
    before the mass blocks are negated) and the final blocks. Both have `summary()` (stage
    name, node count, first 200 characters of the `repr`) and `show(stage)`
    (`print_tree(style="indent")`). `EOMTrace`, `LinearizationTrace` and
    `apply_manifold_rules` are exported from `geomech`.
  - The `geomech` logger: one `DEBUG` line per stage of `compute_eom`,
    `integrate_by_parts` (per target), `to_standard_form` (per block) and `linearize` (per
    row). `GEOMECH_DEBUG=1` (read at import) sets the level, attaches a `stderr` handler if
    none is configured, and runs `check_tree` on the inputs and outputs of the three
    pipelines. Nothing is computed for the log when the level is above `DEBUG`. The
    handler is added only when neither the `geomech` loggers nor the root logger have
    one, so under `logging.basicConfig()` each line is printed once.
  - `geomech.check_tree(expr)` and `InvariantError`: operand types per node class, no
    nested same-class sums, exactly one row operand in a `VVMul`, columns only in `Dot`,
    `Cross` and `Hat`, no `VAdd` mixing rows and columns, no calculus node directly over an
    `S2`/`SO3` leaf, tangent vectors with their point set and named after it, numeric
    scalars consistent with their value. The message names the invariant, the subtree and
    its path. The test suite runs it on every tree that goes into or comes out of a
    successful `compute_eom`, `to_standard_form` or `linearize` call (`tests/conftest.py`).
  - Named rewrite rules: every rewriting case of `vector_rules`, `simplify` (`_eliminate`),
    `pull`, `expand`, `collect`, integration by parts and the manifold rules returns through
    `rule(name, before, after)` (`geomech/utils/trace.py`) with a dotted name
    (`vector.dot_tangent_orthogonal`, `simplify.hat_vee_skew`, `collect.peel_matrix`, …).
    Each module docstring lists its rule names with the hypothesis each needs;
    `tests/test_rule_names.py` checks that the docstring and the code list the same names.
    `geomech.set_trace(True)` or `GEOMECH_TRACE=1` logs every firing on `geomech.rules` as
    `name: before ==> after`. A `simplify.*_fold` line is logged only when the fold itself
    changed the node, not for a rewrite made by a child rule. The scalar triple product
    of `collect` has one name per rotation (`collect.triple_product_left`, `_right`,
    `_left_nested`, `_right_nested`, and the `_mvmul_*` variants).
    `tests/test_rule_coverage.py` checks that every rewrite of random trees logs a rule.
  - A trace is cleared when a pipeline starts (reusing one object keeps the last run
    only); after an exception `failed_at` names the stage that raised and `summary()`
    starts with `FAILED at <stage>: <exception>`.
  - `troubleshooting.md` has a recipe for a wrong equation that matches no row (bisect
    with `summary()`, run `specs/reviews/scripts/engine_fuzz.py --quick`, rerun the
    suspect stage under `set_trace`), and lists the silent bugs fixed in earlier versions
    in a separate history section.
  - Docs: "The Pipeline, Stage by Stage", "Conventions" and "Troubleshooting" (every
    silent case fixed in phases 1–3 with its old symptom, the exact text of every typed
    error a user can hit, the three `AlgebraicError` causes told apart); the manifold
    derivation notes are in the nav as "Manifold Rules: Derivation".
- A kinetic energy written with rows derives: `collect` rewrites an inner product that
  holds the integration-by-parts target, `VVMul(rᵀ, b)`, to `Dot(column_of(rᵀ), b)`, so
  `½ ẋᵀAẋ` and `½ ΩᵀJΩ` go through (the latter prints the same equation as `½ Ω·JΩ`). Two
  new oracle systems (`row_kinetic_energy`, `row_kinetic_energy_so3` in
  `tests/test_eom_oracle.py`). No pinned equation changes: before, any such inner product
  made extraction fail.
- `linearize(..., params=[...])`: the prescribed parameters given to `compute_eom`. They
  are still not supported; a term carrying the variation of a declared param raises
  `AlgebraicError` saying so instead of suggesting a missing input or `Constant`.
- Numbers on either side of the arithmetic operators: `2 * x`, `x + 1`, `1 + x`, `x - 1`,
  `1 - x`, unary `-x` (the same tree as `x * -1`, so its `str` matches the existing
  `x(-1)` output) and `x / 2` (division by a number only; dividing by an expression raises
  `TypeError` suggesting an inverse symbol such as `Scalar('1/m')`). numpy integer and
  floating scalars (`np.int64(2)`, `np.float64(0.5)`) are wrapped like `int`/`float`, also
  on the left (`Expr.__array_ufunc__ = None` makes numpy defer); `bool` is still rejected.
  `x.T()` is available on every expression, not only on vectors.
- `geomech.core.operations.geometry.is_row_vector(expr)`: the row/column orientation of a
  vector expression, decided structurally (`Transpose(column)`, `MVMul(row, M)`, a scalar
  multiple or sum of rows, d/dt/δ/∫ of a row), and `column_of(row)`, the column whose
  transpose is a row, built structurally (`(aᵀK)ᵀ = Kᵀa`, `(aᵀs)ᵀ = a s`).
- A row vector times a matrix goes through the pipeline: `L = ½ m ẋ·ẋ − ½ xᵀK x` written
  `x.T() * K * x`, `(x - e3).T() * K * (x - e3)` and `δW = -(dx.T() * A * v)` derive and
  pass the finite-difference oracle (`tests/test_row_forms.py`, four new systems in
  `tests/test_eom_oracle.py`). `expand` distributes the row shapes and keeps them in row
  form (`(aᵀ(K + J)) b → aᵀK b + aᵀJ b`, `((a + b)ᵀK) c`, `(aᵀ + bᵀ) c`, `((a + b)ᵀ s) c`),
  `extract_linear_coeff` takes the column of a row with `column_of`, and `pull` hoists the
  scalar of a scaled row out of an inner product (`(aᵀs) b → aᵀb s`). A row in a kinetic
  energy (`v.T() * A * v`) is handled by `collect` (see the row kinetic energy entry above).
- `0 + x` is `x` for every expression type, so `sum([a, b, c])` works for vectors and
  matrices and leaves no `(0)` term in a scalar sum.
- `compute_eom(L, dW, variables, params=[...])`: prescribed time-varying parameters (a
  cable length `l(t)`, a moving attachment point `b(t)`) are declared explicitly, mirroring
  `check_eom(params=...)`. Entries must be distinct, non-constant plain
  `Scalar`/`Vector`/`Matrix` leaves that are not configuration variables or tangent vectors
  (`TypeError`/`ValueError` naming the entry otherwise). A declared symbol contributes its
  rate terms (`l̇`) and its variation is not free (`δl` is dropped); an undeclared one is an
  error (see Changed). `to_standard_form` accepts such equations; `linearize` does not yet
  (it raises `AlgebraicError` on `δl`).
- Numeric oracle in the package, `geomech.numeric` (`docs/numeric-verification.md`):
  `evaluate(expr, env)` evaluates a tree with leaf values keyed by `str(leaf)` and draws
  the rest deterministically; `check_eom(L, dW, variables, eom)` returns the worst relative
  error of the output of `compute_eom` against a finite-difference Euler–Lagrange residual
  computed without the library's variational rules (R, R³, S² on the tangent plane, SO(3)
  in body-frame form, products of these, prescribed time-varying parameters via `params`);
  `assert_eom` raises naming the offending variable. `Evaluator`, `hat`, `vee`, `rand_rot`,
  `expm_so3` and `finite_difference_check` moved here from the test helpers, which now
  re-export them; the evaluator now honours `UnitNorm`, `OrthogonalMatrix` and
  `SkewSymmetricMatrix` when drawing a leaf. Every pinned system (`tests/test_eom_oracle.py`,
  36 systems) passes at 1e-5 over 10 seeds (worst correct score 2.3e-7), and flipping any
  single term of any of them is detected (worst score 2.1e-2). The default step is
  `h=1e-4` (`numeric.DEFAULT_H`): `d/dt ∂L/∂v` nests two central differences, so the
  roundoff is `≈ 1e-16 |L| / h²`, which at `1e-5` was `≈ 1e-6` and let a correct equation
  score up to `1.4e-5`. The error is relative to the size of the whole equation, so a term
  that is small next to a large one (a stiff `k = 1e7` next to gravity) can be wrong
  without moving it; `assert_eom` now checks that every additive term is resolved at
  `rtol` (`2 |P t| / max(1, |E_true|) ≥ rtol` at some trial) and raises naming the term,
  its size and the resolution (`require_resolved=False` to opt out). The grey zone below
  finite-difference noise remains and is documented (`docs/numeric-verification.md`).
- Variation-based linearization: `linearize(eom, variables, inputs)` turns the output of
  `compute_eom` into the symbolic blocks of `M δẋ = A δx + B δu` along an arbitrary
  trajectory, with manifold-respecting perturbation variables (δs, δx, ξ, η and the
  variations of the velocities) and kinematic rows for R, R³, S² and SO(3)
  (`LinearizedSystem`, `docs/vbl.md`). Verified against hand derivations for the point mass,
  spring-mass, free rigid body, rigid pendulum, spherical pendulum and quadrotor, and by
  finite differences for those and the double rigid pendulum, the elastic spherical pendulum
  and the quadrotor with a cable-suspended load. Every additive term of δE must carry
  exactly one perturbation symbol; a term with none, with two, or with one hidden inside a
  calculus node raises `AlgebraicError` (no silent zeros). Inputs must be distinct,
  non-constant `Scalar`/`Vector` leaves that are neither configuration variables nor
  tangent vectors (`Ω`, `ω`), and each equation must have the type of its variation vector;
  violations raise `ValueError`/`TypeError`.
- Scalar configuration variables (`SystemVariables(scalars=[...])`) participate in
  `compute_eom` and `to_standard_form`. Spring-mass, spring-mass-damper, coupled spring-mass
  and the elastic spherical pendulum (R × S²) derive end to end.
- `extract_linear_coeff` replaces `extract_coeff` with a documented 3×3 dispatch grid
  (query type × expression type) and typed errors `NonLinearError` and `AlgebraicError`;
  scalar queries are supported inside `Hat`, `Transpose`, `Vee`, outer products and
  row-vector products.
- Integration by parts for scalar targets inside `Mul` chains.
- LaTeX rendering (`to_latex`, `display_eom`, `display_standard_form`) and example
  notebooks for the rigid pendulum, spherical pendulum, quadrotor and quadrotor with payload.

### Changed
- Engine errors reaching the user through a pipeline are re-raised with the same type,
  the original chained as `__cause__`, and the pipeline context: `compute_eom` names the
  variation vector being extracted and, when a time derivative of a variation vector
  survived integration by parts, the term it survived in and that `collect` could not
  expose it (pointing to `print_tree` and the `ibp.collect[...]` trace stage); a `collect`
  error during integration by parts is prefixed with the stage; `to_standard_form` names
  the equation and the term, and a jerk says the equation is not second order (the two
  causes of the "appears only inside the calculus operator" error are now told apart);
  `linearize` keeps the inner message when a coefficient still contains a perturbation,
  and its "no perturbation symbol" hint names the stray variation and lists every remedy.
- The transformation modules no longer share their names with their functions:
  `expand.py`, `pull.py`, `simplify.py`, `vector_rules.py` and `core/math/collect.py` are
  now `expansion.py`, `pulling.py`, `simplification.py`, `vector_identities.py` and
  `collecting.py`. The packages still export the *functions* under the old names, so
  `from geomech.core.transformations import simplify; simplify(expr)` is unchanged (a
  plain function: picklable, with its own signature), and
  `import geomech.core.transformations.simplification as m` gives the module. Before,
  the package re-exported the function under the module's name and hid the module.
  Imports through the old module paths now raise `ModuleNotFoundError`, for example
  `from geomech.core.transformations.simplify import simplify` and
  `from geomech.core.math.collect import collect`; use the new module names, or import
  the functions from `geomech`, `geomech.core.transformations` or `geomech.core.math`,
  which are unchanged.
- A Lagrangian in accelerations (`½ ẍ·ẍ`) is reported as such by `compute_eom`
  ("the Lagrangian depends on accelerations …; higher-order Lagrangians are not
  supported"), instead of as a term that `collect` could not expose with a trace stage
  that did not exist. A surviving first-order derivative already in the IBP shape points
  to `ibp.integrate[t]`; a scalar one is described as a factor of a product, not as the
  operand of a `Dot`.
- The public `integrate_by_parts` logs its stages on the `geomech` logger at `DEBUG`, as
  its docstring said.
- `full_simplify`'s cap warning says "no confirmed fixed point after N passes" (the last
  allowed pass changed the expression, so no further pass confirmed it).
- A composite `SystemVariables` entry (`x + y`, `m * z`, `q.t_diff()`) raises `TypeError`
  "is a VAdd expression, not a single symbol"; before, it got past validation and failed
  later with a stray-variation message that did not name the entry.
- `render_eom` reports through the `geomech` logger (`WARNING` with the end of the
  `pdflatex` log on failure, `INFO` with the path on success) instead of `print`.
- `is_row_vector`/`column_of` import their node classes once at module level instead of on
  every call.
- `to_latex` is precedence-aware: composite operands of a cross product, of a matrix
  product and of a transpose are wrapped in `\left( \right)` (`a \times \left(b \times
  c\right)`, `R^T \left(a \times b\right)`, `\left(R a\right)^T`), a dot product
  multiplied by anything is wrapped (`\left(b \cdot c\right) a`), the numeric factors of
  a product are folded into one coefficient printed first (`-\frac{1}{2} s s`, `2.5 s`),
  a scalar multiple of a cross product is wrapped (`m \left(a \times b\right)`), the
  time derivative or variation of a composite prints `\frac{d}{dt}\left(…\right)` /
  `\delta\left(…\right)`, and a transposed scalar multiple prints `s v^T`. Every
  ambiguous pair of the dynamics audit (F8) now prints differently, and no output contains
  a double superscript (`^T^T`). `print_eom`, `eom_to_latex` and `render_eom` typeset the
  equations with `to_latex` instead of embedding `str(eqn)` in LaTeX math.
- `repr` shows the leaf type and constness: `v:x`, `M:J`, `S2:q`, `SO3:R`, a trailing `*`
  for a constant (`m*`, `M:J*`), numbers bare; products are joined with `·` so `*` is never
  a multiplication (`m*·g*·<v:e3*, v:x>`). `str`, the equality key, is unchanged.
- `print_tree` defaults to `style="indent"`; each leaf lists every flag that is set
  (now including `is_identity`, `is_skew_symmetric`, `is_orthogonal`) and a `TS2`/`TSO3`
  names its manifold point (`of=SO3 'R'`); when two distinct points share a name within
  one call they are labelled `#1`, `#2`. The output is deterministic (no addresses). `tree_str` keeps `style="topdown"` as its default.
- `to_latex` wraps a factor whose LaTeX starts with `-` (a negated transpose) unless it
  leads an unsigned product: `x \left(-a^T b\right)`, `A \left(-K^T\right)`,
  `-\left(-K^T\right)`, `a \cdot \left(-b\right)`; before, these printed `x -a^T b` (read
  as a subtraction) and `--K^T`.
- Like terms whose float coefficients cancel to roundoff (`|sum| ≤ 1e-12 ·` the largest
  coefficient) are dropped: `full_simplify(x*0.1*3 - x*0.3)` is `0`, was `x(5.55e-17)`.
  `simplify` rewrites a transposed scaled row `((a)'x)'` to the column `ax`.
- The `tree_str`/`print_tree(style="topdown")` labels of the products are `S·V`, `S·M`,
  `M·V`, `M·M`, `V·V` (were `S*V`, …), since `*` marks a constant leaf.
- Tree-building errors name the operands: `Dot cannot combine SCALAR 'm' and VECTOR 'x'`,
  `Hat expects a VECTOR, got SCALAR 'm'`, `Sub cannot combine VECTOR 'x' and SCALAR 'm'`;
  `J * x.T()` (was the bare `operation cannot be performed`), `x * J`, `x * y` and
  `x.T() * y.T()` explain the row/column problem and what to write instead. The operands
  are quoted as written: `(p + q) + a` names the sum `p + q`, not its first addend, and
  `a - 1` reports `Sub cannot combine VECTOR 'a' and SCALAR '(1)'` (was the internal
  `VAdd … '(-1)'`). Non-expression operands in a `TypeError` are shortened past 80
  characters like expressions. A non-expression
  operand (`m * "x"`, `Add(m, None)`) raises `TypeError` naming it instead of an
  `AttributeError` from inside the package. `ExpressionMismatchError` takes the operands as
  optional `lhs`/`rhs` arguments; the class and its base are unchanged.
- `t_integrate()` of a constant `Scalar`/`Vector`/`Matrix` raises `AlgebraicError` naming the
  symbol (was an empty `NotImplementedError`).
- Unknown `attr=[...]` strings raise `ValueError` naming the symbol and listing the valid
  attributes (`Constant`, `Zero`, `Ones`, `Identity`, `UnitNorm`, `SymmetricMatrix`,
  `SkewSymmetricMatrix`, `OrthogonalMatrix`, `Manifold`, `TangentVector`). Before, a typo
  such as `attr=["constant"]` was ignored and the symbol was silently non-constant, so `ṁ`
  and `δm` terms appeared in the equations. A bare string (`attr="Constant"`) raises
  `TypeError: attr must be a list of strings, got str 'Constant'`; before, it was iterated
  character by character and reported as `unknown attribute 'C'`.
- `compute_eom` rejects a time derivative of a listed variation vector in `inf_work`
  (`Dot(x.t_diff().delta(), f)`, `Dot(ξ.t_diff(), f)`, `δ(q̇)·f` on S², `η̇·M` on SO(3))
  with a `ValueError` naming the term: virtual work is linear in the variations
  themselves. Before, integration by parts moved the derivative onto the force and the
  equation silently carried `−ḟ` instead of `f`. A velocity-dependent force
  (`δx·(−c ẋ)`) is unaffected; δL is not checked.
- A `SystemVariables` entry or a `params` entry printed like the tangent vector (`ω_q`,
  `Ω_R`) of a listed `S2`/`SO3` point raises `ValueError` naming both (params also may not
  reuse a variation vector name). Before, two different leaves with one printed name were
  merged by printing, `check_eom` and `linearize`, and the equation was wrong.
- `params` given as a bare leaf and a `SystemVariables` group given as a bare leaf
  (`vectors=x`) raise `TypeError` saying a list is expected; `compute_eom` with no
  configuration variables raises `ValueError` (an empty `SystemVariables()` can still be
  built and passed to `check_eom`). `Number(True)`, `x * True` and `compute_eom(L, False,
  …)` raise `TypeError` (bool was taken as the integer 1/0). A hand-built
  `Variation(x.t_diff())` is reported under `x` with a hint to write `x.delta()`.
- `compute_eom` raises `ValueError` when a variation symbol survives extraction that is
  neither the variation of a listed configuration variable nor of a symbol declared in
  `params` (`δy` of an unlisted `y`, `δe3` of an `e3` missing `Constant`, `ξ`/`δω` of an
  `S2("q")` listed as `Vector("q")`, `δp` from `δW = … + δp·F`), naming the symbols and
  saying to list them in `SystemVariables`, declare them in `params=...`, or mark them
  `Constant`. Before, those terms were dropped from every equation without a message, so a
  forgotten degree of freedom lost its equation silently and a prescribed parameter was
  indistinguishable from a typo. Prescribed time-varying parameters are now declared with
  `compute_eom(..., params=[l])` (see Added).
- `SystemVariables` entries are validated by `compute_eom` and `linearize`: a printed name
  listed twice (across groups, or a `Scalar` and a `Vector` of the same name), a tangent or
  variation vector (`ω`, `ξ`, `Ω`, `η`) or a `TimeDerivative`/`Variation` node listed as a
  configuration variable, and a non-`SO3` matrix raise `TypeError`/`ValueError` naming the
  entry; `linearize._check_eom` also compares the number of keys. Before, a duplicate merged
  two equations into one and `linearize` emitted duplicate rows and states, a tangent vector
  gave the meaningless equation `(ω m + f)`, and a plain matrix failed later with an
  `AlgebraicError` from `extract_linear_coeff`/`collect` that did not name the entry.
- The S² tangent rules (`ω·q = 0`, `q × (ω × q) = ω`) compare the tangent vector's parent by
  name, as every other rule does, so a second `S2("q")` object is the same point. Before,
  the rules were skipped for it and `ω·q`, `q×(ω×q)` terms stayed in the equation.
- `compute_eom(L, None, …)` and `compute_eom(L, 0, …)` mean δW = 0; a non-scalar or
  non-expression `lagrangian`/`inf_work` raises `TypeError` naming its type. Before, `None`
  and `0` crashed with `AttributeError` and a vector `dW` with `ExpressionMismatchError: Add`.
- `apply_manifold_rules` raises `TypeError` naming the node class when it cannot rebuild a
  node whose children were substituted. Before, the original node was returned and the
  substitution silently discarded.
- `TimeDerivative(q)` and `Variation(q)` for an `S2`/`SO3` leaf raise `ValueError` pointing
  to `q.t_diff()`/`q.delta()`. Before, the node was built and carried none of the
  kinematics.
- `collect` raises `TypeError` for a non-scalar argument, and a target on both sides of a
  `Dot` raises `NonLinearError` before the flip shortcut. Before, a non-scalar argument was
  returned unchanged and `Dot(f(v), v)` with `v` inside `f` was flipped and accepted as
  linear.
- `full_simplify` emits a `RuntimeWarning` carrying the expression when it stops at
  `max_iter` passes without reaching a fixed point. Before, the partially simplified
  expression was returned silently.
- `Number(0)` carries the `Zero` flag and `Number(1)` the `Ones` flag (also for the numbers
  wrapped by `2 * x`-style operators), and a product of numerics that folds to 0 is the
  `Zero` leaf. Before, `Number(0) * x` simplified to `(0)x` with `is_zero == False`.
- `AlgebraicError`, `NonLinearError`, `ExpressionMismatchError` and `SizeMismatchError` are
  exported from `geomech` (and `geomech.utils`) so user code can catch them.
- `examples/rigid_pendulum.py`: the second inertia matrix of the double rigid pendulum is
  named `J2` (it was a second `Matrix("J1")`, so both links printed the same symbol).
- `expand` distributes `VVMul` over `VAdd`: outer products `(a + b) cᵀ`, `a (b + c)ᵀ` and
  inner products `(a + b)ᵀ c`, `aᵀ (b + c)` are split into sums of products, as every other
  product already was. `Transpose` itself is still not expanded. The printed output of
  `compute_eom` and `to_standard_form` changes only for Lagrangians that contain an outer or
  inner product with a sum as a factor (for example `ẋ·(p pᵀ) ẋ` with `p = x − l q`); the
  result is the distributed, mathematically identical expression. None of the pinned
  systems is affected.
- `Expr.__eq__` now returns `False` across different node types before comparing strings.
- `extract_linear_coeff` raises `AlgebraicError` instead of returning zero when the query
  occurs only inside a calculus operator (e.g. a jerk term, or `δẍ` from a Lagrangian in `ẍ`)
  or inside an unsupported node. `to_standard_form` partitions acceleration and input terms
  with the same test, so such terms can no longer be dropped silently, and raises
  `NonLinearError` when a term is bilinear in two acceleration variables.
- `compute_eom` validates `SystemVariables`: a constant, numeric or zero entry raises
  `ValueError`; an entry of the wrong expression type raises `TypeError`.

### Removed
- `UndefinedCaseError` (unreachable: every node type is SCALAR, VECTOR or MATRIX; it
  carried no context) is removed from `geomech.utils.errors`. It was not exported from
  `geomech`.

### Fixed
- `Transpose(R) * v`, `Transpose(R) * M` and `Transpose(R) * s` raised
  `ExpressionMismatchError` (the transpose was always treated as a row vector); they now
  build `MVMul`, `MMMul` and `SMMul`. `Transpose(v) * s` built the double transpose
  `Transpose(SVMul(Transpose(v), s))`, a column that the engine treated as a row and that
  `to_latex` printed as `s v^T^T`; it is now the row `Transpose(SVMul(v, s))`. The other
  order, `s * Transpose(v)`, was not affected and is unchanged: it builds the row
  `SVMul(Transpose(v), s)` (`(v)'s`). The two orders therefore build different trees
  (`(v)'s` and `(vs)'`); `simplify` rewrites `(s vᵀ) w` to `s (vᵀ w)`, but not yet
  `(vᵀ s) w`.
- `Transpose(Transpose(x))` is `x` at construction, and the row/column orientation is
  decided by `is_row_vector` instead of "the operand is a `Transpose` node". Legal trees
  that `expand`/`simplify` used to break with `ExpressionMismatchError` (engine fuzz
  report, finding 6: `(a + (bᵀK)ᵀ) cᵀ`, `cᵀ(a + (bᵀK)ᵀ)`, `((aᵀK)ᵀ)ᵀ Q`) now go through;
  `extract_linear_coeff` takes the column of a row with `column_of(row)` (the same object
  as `Transpose(row)` for a `Transpose` node), so it also accepts rows that are not a
  `Transpose` node. No pinned equation changes: `str` of every system
  through `compute_eom`, `to_standard_form` and `linearize` is byte-identical.
- The S² variation rule substituted `δω = ξ̇ − ω × ξ`, the tangential projection of the true
  `δω = ξ̇ − 2 ω × ξ` (Lee et al. eq. 5.13 with `qqᵀξ̇ = ω × ξ` on the constraint set). Every
  Lagrangian whose `∂L/∂ω` has a component along `q` (gyroscopic `c·ω`, anisotropic
  `½ ω·Aω`, rotating frame `|ω + W₀|²`, couplings `ẋ·ω`, `ω·Ω`, `ω·(Rᵀe₃)`) came out with its
  gyroscopic terms at half strength; pendulum-type systems (`∂L/∂ω ⊥ q`) were unaffected
  and their pinned equations are unchanged. The one shipped system whose printed output
  changes is the quadrotor with a cable-suspended load (`tests/vbl_systems.py`, notebook):
  in the load equation the term parallel to `q`, `ω × (q × ẋ_L) m_Q l`, changed coefficient
  from −1 to −2; it is invisible on the tangent plane and does not affect the dynamics.
  `TS2.delta(substitute=True)` / `TSO3.delta(substitute=True)`, an unused second copy of the
  formulas, were removed.
- `collect` rewrote `(M v)·(b × c)` as `v·((M b) × c)`, a false identity, which produced a
  wrong SO(3) equation for `½ Ω·JΩ + (JΩ)·((Rᵀe₃) × b)`; and it dead-ended when the target
  sat more than one level under a matrix inside a `Dot`, so `(R(JΩ))·e₃`, `½ q̇·A q̇` on S²,
  `Ω₁·(R₁ᵀR₂Ω₂)` and `½ (Ω×ρ)·J(Ω×ρ)` raised `AlgebraicError`. The matrix is now peeled off
  the side holding the target, `Dot(M f(v), w) → Dot(f(v), Mᵀ w)`, but only when the target
  is in the vector operand `f(v)`: a target inside the matrix factor (`Hat(v)`, `v cᵀ`,
  `M Hat(v)`) cannot be exposed by transposition and raises `AlgebraicError` naming the term
  and the target (the first version of the peel moved such a matrix to the other side and
  looped until `RecursionError`, e.g. for `½ Ω·JΩ + a·((Ω cᵀ) b)`). A target on both sides of
  a `Dot` raises `NonLinearError` instead of `RecursionError`.
- `simplify` rewrote `Hat(Vee(X)) → X` for any matrix, so `Hat(Vee(I))` became `I` (it is
  zero). The rule now fires only for skew-symmetric `X` (`Hat(·)`, `SkewSymmMatrix`, and
  `SMMul`/`MAdd`/`Transpose` of such, which now report `is_skew_symmetric`). `SkewSymmMatrix`
  leaves did not carry their flag at all (it was appended after the flags were built).
- `extract_linear_coeff` returned a silent zero for a vector query inside `Hat`, `Transpose`
  or an `MMMul` chain in a matrix expression (`M·Hat(v) → M·0`) and for a matrix query in the
  scalar factor of an `SMMul`; these raise `AlgebraicError` ("no left matrix coefficient").
  `Cross(Cross(a, q), q)` and `MVMul(Hat(q), q)` returned a coefficient that still contained
  the query; they raise `NonLinearError`. A matrix query in a scalar or vector expression
  raises `AlgebraicError` instead of `NotImplementedError`.
- `simplify` turned an outer product with a zero factor into the scalar `Zero`, so
  `MAdd(VVMul(0v, aᵀ), M)` raised `ExpressionMismatchError`; it is now `ZeroMatrix`.
- Integration by parts dropped `d/dt` of a non-constant scalar factor multiplying
  `Dot(ξ̇, ·)`, so every system with a time-varying scalar in its kinetic energy (elastic
  pendulum, prismatic joint, prescribed `l(t)`) lost its Coriolis / rate-of-inertia terms
  (`−2 m l l̇ ω`, `−m l̈ e₁`, `J̇ Ω`). The prefactor is now differentiated with the product rule.
- `collect` could not move the variation out of `Dot(R₁ a, R₂ b)` with two different rotations,
  which dropped the cross-inertia coupling of multi-body SO(3) systems. The SO(3)×SO(3) double
  rigid pendulum now has `M₁[Ω̇₂] = m₂ l̂₁ R₁ᵀ R₂ ρ̂₂` and `M₂[Ω̇₁]` its transpose.
- `collect` no longer recurses forever on `Dot(M a, N b)` when the target is nested in the
  vector operand in a form the rules cannot expose (`SVMul`, `VAdd`); the expression is
  returned unchanged. (A target inside the matrix factor raises, see above.)
- `SkewSymmMatrix`, `SO3`, `S2`, `Matrix`, `Vector` and `Scalar` copy the `attr` list they are
  given; `SkewSymmMatrix`, `SO3` and `S2` appended their own flags to the caller's list, so
  two leaves built from one list shared each other's flags (`Matrix("J", attr=shared)` was
  skew-symmetric after `SkewSymmMatrix("K", attr=shared)`).
- `Add()`, `VAdd()` and `MAdd()` with no operands raise `TypeError`; an empty `MAdd` was
  vacuously symmetric and skew-symmetric (`all([])`). Nothing in the package built one.
- `extract_linear_coeff` compared nodes by rendered string only, so `Mul(k, x)` could be
  mistaken for a leaf named `kx`.

## [0.1.0] — 2026-03

Initial release: expression trees, S² and SO(3) manifold support, `compute_eom`,
`to_standard_form`, printing utilities.
