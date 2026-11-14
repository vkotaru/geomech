# Changelog

All notable changes to `geomech` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Upgrading from 0.1.0

The 2026-11-13 fix plan (`specs/fix-plan.md`) changed derived equations and the public
surface. In short:

- **Equations on S² change value** when `∂L/∂ω` has a component along `q`: the variation
  rule is now `δω = ξ̇ − 2 ω×ξ` (Lee et al. eq. 5.13). Gyroscopic and coupling terms
  were at half strength before. A false `collect` identity that corrupted
  `Dot(M v, b×c)` terms is also fixed. Re-derive any such system.
- **Printed form**: products print with `*` (`m*l*l*(-1)`), numbers come last, and
  equations are shorter after the simplifier changes. Any string you compared against
  changes; values do not.
- **Equality**: `Mul(ab, c) != Mul(a, bc)`; calculus nodes and `Transpose` compare their
  operand's class; symbol names may not contain `*`.
- **Immutable nodes**: node classes are frozen dataclasses with tuple operands.
- **Module renames**: `core.transformations.{expand, pull, simplify, vector_rules}` and
  `core.math.collect` are now `expansion`, `pulling`, `simplification`,
  `vector_identities`, `collecting`. The functions keep their import paths from
  `geomech`, `geomech.core.transformations` and `geomech.core.math`.
- **New errors on input that used to pass silently**: unknown `attr` strings, a
  configuration variable missing from `SystemVariables` (declare prescribed time-varying
  symbols with `compute_eom(..., params=[...])`), duplicates, a time derivative of a
  variation in the virtual work, a calculus node over a manifold leaf, and
  `extract_linear_coeff` raising instead of returning zero. `docs/troubleshooting.md`
  lists every message.
- **Warnings**: `full_simplify` emits a `RuntimeWarning` at its iteration cap, which is
  an exception under `-W error`.
- **New tools**: `geomech.numeric` (`check_eom`, `assert_eom`, `evaluate`),
  `linearize`, `compute_eom(..., trace=EOMTrace())`, `GEOMECH_DEBUG`, `GEOMECH_TRACE`,
  `check_tree`. See `docs/pipeline.md`.

Open follow-ups are listed in `specs/follow-ups.md`.

### Added
- Simulation, `geomech.simulation` (roadmap feature 6, `specs/simulation.md`,
  `docs/simulation.md`); it consumes only `geomech.codegen` modules and imports scipy
  (extra `geomech[sim]`) only for `solve_ivp` and `lqr`:
  - `accelerations(module, state, u, p)`: solves `mass @ a + bias + G u = 0` with each S²
    row block replaced by `P_q · row + q qᵀ` on its own `ω̇` column (right-hand side
    `−P_q (bias + G u)`), so `q · ω̇ = 0` and the tangent-plane equation holds whether the
    derived S² mass block is singular (`Hat(q) Hat(q) m l²`) or full rank (`−m l² I`).
    A singular system (condition number > 1e12) raises `LinAlgError` with the condition
    number and the state.
  - `simulate(module, x0, t_span, p, u=None, method="rkmk4", dt=1e-3)` → `Trajectory`
    (`t`, `state`, `constraint_drift`, `energy()`, slot access `traj["q"]`). `"rkmk4"`:
    fixed-step Runge–Kutta–Munthe-Kaas of order 4 (S² rotated from the left, SO(3) from
    the right, `dexp⁻¹` to the `ad²` term, classical RK4 on vector slots); constraint
    drift at round-off. `"solve_ivp"`: DOP853 on the embedded state, outputs projected
    (`q/|q|`, polar factor of `R`), drift before projection reported. `x0` is checked
    (unit `q`, `ω ⊥ q`, `R ∈ SO(3)`) with `ValueError` naming the slot; `u` is `None`, a
    constant or a callable `u(t, state)`.
  - `lqr(lin_module, x_eq, u_eq, p, Q, R)` → `LQRController` (gain, Riccati solution,
    reduced `A`, `B`, `error(state)`, callable as `u(t, state)`), with SO(3) error
    `vee(log(R_eqᵀ R))` and S² perturbations reduced to a 2-D tangent basis
    (`tangent_basis`). `geomech.numeric.log_so3` (principal logarithm, accurate near 0
    and π).
  - `examples/inverted_rigid_pendulum_lqr.py`: Lagrangian → equations → generated code →
    LQR → closed-loop simulation of an inverted rigid pendulum on SO(3).
  - Tests `tests/test_simulation*.py`: energy conservation over 10 s for six systems and
    the momentum `R J Ω` of the free rigid body with both integrators, the spherical
    pendulum (Lee et al. Eq. 5.14) and Euler's equations integrated in NumPy (1e-8 over
    2 s), the order of RKMK4, constraint drift, the S² solve on every S² system of the
    oracle gate, the LQR demos, `x0` checks and a 1 s run of every VBL system.
- Code generation, `geomech.codegen` (roadmap feature 5, `specs/codegen.md`,
  `docs/codegen.md`):
  - `StateLayout(variables, inputs, params)`: the flat `state` (configuration slots, then
    velocity slots; R row-major, S² `q` and `ω`), `a` (accelerations), `u` (inputs) and
    `p` (a dict keyed by `str(symbol)`), with `slices`, `to_dict()` and conversions to and
    from `Evaluator` envs; `collect_params` lists the parameters of a set of expressions.
  - `to_source(expr, layout, target)` and `lambdify(...)` (`f(state, u, p, a=None)`), and
    `generate_module(system, layout, path, target)` for the dict of `compute_eom`
    (`residual`), of `to_standard_form` (`mass`, `bias`, `input_matrix`) and a
    `LinearizedSystem` (`M_lin`, `A_lin`, `B_lin`), each with `kinematics`, `LAYOUT` and
    the `_hat`/`_vee` helpers. Generated modules are ruff-clean and need only `numpy`.
    Leaves are checked: a stray symbol, a `TimeDerivative` of a non-leaf, a missing or an
    extra parameter raise naming the symbol.
  - Targets `"numpy"` (default) and `"casadi"` (optional extra `geomech[casadi]`, pinned
    to `casadi==3.8.2`; functions accept CasADi symbols, so `ca.jacobian` applies).
  - `energy_function(L, layout)` (and `energy_expr`, `momenta`): the total energy
    `Σ (∂L/∂v)·v − L`, with `∂L/∂v` extracted from the first variation of `L`;
    `generate_module(..., lagrangian=L)` adds `energy(state, p)`.
  - `examples/generated/`: the standard-form and linearized modules of the nine systems of
    `tests/vbl_systems.py`, rewritten by `scripts/regenerate_generated.py`; CI regenerates
    them and fails on a diff. Optional extra `geomech[sim]` (scipy) for integrating them.
  - Tests `tests/test_codegen*.py`: every block of every oracle, phase-1 and VBL system
    against `evaluate` at 20 admissible states (1e-12), CasADi against NumPy, each node
    class of the translation table, the layout and its errors, `∂L/∂v` and `E` against
    finite differences, `dE/dt = 0` along the unforced dynamics, and a 10 s free rigid
    body run with `solve_ivp` conserving the generated energy to 1e-8.
- Structure and typing (phase 5b of `specs/fix-plan.md`, items 5.5–5.7, decisions D2 (a)
  and D3 (a)):
  - `geomech.core.base.constants`: the literals `NEG_ONE`, `NEG_TWO`, `HALF` (printed
    `(-1)`, `(-2)`, `(0.5)`, the same as `Number(...)`), `negate(expr)` (the tree of
    `-expr`) and `zero_like(expr)`/`zero_of(expr_type)`; every rewrite builds its −1, −2
    and typed zeros from them.
  - `geomech.core.transformations.utils`: `addends`, `sum_of`, `flatten_mul`, shared by the
    simplifier, integration by parts, `compute_eom`, `to_standard_form` and `linearize`.
  - Public names for helpers other modules used privately: `has_algebraic`, `same_node`
    (`geomech.core.math.extract`), `check_variables` (`geomech.dynamics.eom`),
    `wrap_numeric` (`geomech.core.operations.multiplication`; `_wrap_numeric` stays as an
    alias), `term_key` (`geomech.core.transformations.simplification`).
  - Type annotations: the operation mixins declare `nodes` (and are `Expr` for the type
    checkers), `Expr` declares `type`, `str | None` where a default is `None`, and the
    public functions of `geomech.dynamics`, `geomech.numeric` and `geomech.utils` are
    annotated. `mypy geomech --ignore-missing-imports` passes and runs in CI (lint job);
    `mypy` is in the `dev` extras, pinned to 2.4.1 there and in CI. pyright: 109 errors at the start of the phase, 5 now.
- Simplifier completeness (phase 5a of `specs/fix-plan.md`); no derived equation changes
  value, most print shorter:
  - `vector_rules` and `pull` recurse through every node class with one structural helper
    (`geomech.core.transformations.utils.rebuild`/`map_children`, shared with
    `apply_manifold_rules`), so the geometric identities fire below a `Dot` operand, an
    `MVMul`, `MMMul`, `SMMul`, `MAdd`, `VVMul`, `Hat`, `Vee`, `Transpose` or a calculus node:
    `Dot(v, q×(ω×q)) → v·ω`, `J(a×a) → 0v`, `Q(Ω×Ω) → 0`. New `pull` rules
    `pull.row_mvmul_scaled`, `pull.vvmul_scaled`, `pull.transpose_scaled`,
    `pull.hat_scaled`, `pull.vee_scaled` (`Hat(v s) = Hat(v) s`, `(v s)ᵀ = vᵀ s`, ...).
  - `expand.vee_over_madd`: `Vee(A + B) = Vee(A) + Vee(B)`.
  - `simplify.transpose_skew` (`Sᵀ = −S` for a skew-symmetric leaf, sum or transpose),
    `simplify.mmmul_reassociate` and `simplify.mmmul_combine` (the `MMMul` mirror of the
    `MVMul` associativity: `A Rᵀ R → A`, `Rᵀ(R A) → A`); `IdentityMatrix` is flagged
    symmetric and `Iᵀ` is the identity.
  - Like terms are compared on flattened products: an addend is split into the product
    of all its numeric literals, its other scalar factors (a multiset) and its vector or
    matrix base. `x(2)y + xy(3) → xy(5)`, `v a − v a → 0v`, `(vᵀ s) w + s (vᵀ w)` merge;
    the three `k Ω₁ × QΩ₂` terms of the SO(3)×SO(3) coupling now sum to one.
  - `is_number_literal(expr)` (`geomech.core.base.expressions`): a scalar leaf printed as
    its own value (`Number(2)`, `Zero`, `Scalar("0.5", value=0.5)`).
  - `ibp.product_rule` names the rate `ȧ` that `ibp.dot_prefactor` builds, and the six
    `collect.flip` call sites have their own names (`collect.flip_target_right`,
    `collect.flip_cross_mvmul`, `collect.flip_mvmul_pair`, `collect.flip_mvmul_left`,
    `collect.flip_mvmul_right`, `collect.flip_cross_left`).
  - `specs/reviews/scripts/pinned_snapshot.py`: records every pipeline call of the test
    suite and every builder system at one commit and replays them later, comparing each
    output's value at fixed manifold states (and each linearization row along admissible
    directions) and `check_eom`; the gate for changes that reprint pinned equations.
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
- `to_standard_form` validates `inputs` like `linearize`: both call the new
  `geomech.dynamics.eom.check_inputs`, so they accept and reject the same inputs with the
  same message after the caller's name. Before, `to_standard_form` only ran the name-clash
  check: `inputs=[Scalar("g", attr=["Constant"])]` returned a `G[g]` block, a duplicate
  input was merged, `ξ`/`η` or a `Matrix`, an `S2` point or `F + u` returned an empty or
  wrong `G` entry, and a configuration variable or `ω`/`Ω` failed later with an engine
  error. Each now raises the `linearize` error (`TypeError` for a non-leaf, `ValueError`
  otherwise). A generator of inputs is no longer consumed before the extraction (it gave
  an empty `G`), and one leaf passed instead of a list (`inputs=f`) raises `TypeError`
  "inputs must be a sequence of Scalar/Vector leaves" in both functions.
- An input printed like a different leaf of the equations raises `ValueError` naming both
  (formerly item 1 of `specs/follow-ups.md`). With `L = ½ m ẋ·ẋ − m g e₃·x` and `g` `Constant`,
  `linearize(eom, V, [f, Scalar("g")])` used to return a zero `B[δg]` column (the
  constant has no variation) while `to_standard_form(eom, V, inputs=[Scalar("g")])`
  matched it by name and returned `G[g] = −m e₃`; both now raise "input 'g' has the name
  of the constant 'g' in the equations; a symbol is either a constant parameter or an
  input". The same holds for a different class (`Scalar("F")` against `Vector("F")`),
  other flags or another size, and for `compute_eom(..., params=[...])` against a leaf of
  `L` or `δW` (a declared param named like a constant had no effect). An input that is
  the same symbol as its leaf in the equations is unchanged. The check is
  `geomech.dynamics.eom.check_symbol_clash`.
- Products print a `*` between their factors (phase 5b, decision D2 (a)): `Mul`
  (`m*l`), `SVMul` (`x*m`), `SMMul` (`J*(2)`), `MVMul` (`J*\Omega`), `MMMul` (`R*J`),
  `VVMul` (`(a)'*b`). Every pinned equation string changes form; each changed output was
  checked to be the old string with separators inserted and numerically identical
  (`specs/reviews/phase5b-string-changes.md`). `repr` is unchanged (`·` between factors).
- Equality (D2 (a)): two products whose factors concatenate alike are no longer equal
  (`Mul(ab, c) != Mul(a, bc)`); a calculus node also compares the class of its operand
  (`Variation(Scalar("x")) != Variation(Vector("x"))`); `expr.has(x)` uses the same `==`,
  so `Scalar("m")` does not occur in `Dot(Vector("m"), w)`. A scalar and a vector variable
  with one name now extract separately.
- The product separator is reserved (phase 5b review F2): a symbol name containing `*`
  (`Scalar("m*l")`, `Vector("x^*")`, any leaf class or `Number("a*b")`) raises
  `ValueError`, since `Mul(Scalar("m*l"), g)` and `Mul(m, Mul(l, g))` would both print
  `m*l*g` and compare equal. Write `ml` for a lumped constant and `x^\ast` for a star.
  `Transpose` compares the class of its operand like the calculus nodes
  (`Transpose(Vector("x")) != Transpose(Matrix("x"))`, review F3).
- Nodes are immutable (D3 (a)): every node class is a frozen dataclass, operands are a
  tuple (`expr.nodes == (a, b)`), `ExprFlags` and `ManifoldInfo` are frozen. Assigning to a
  node raises `FrozenInstanceError`; pickling, `copy` and `deepcopy` work, and a tree
  pickled before (operands in a list) unpickles with tuple operands.
- Like terms are grouped in one pass by a canonical key (`term_key`) instead of pairwise,
  and an inner product of two columns is symmetric: `aᵀb s − bᵀa s → 0`.
- Number literals are normalised by `simplify`: a bare `Number(2.0)` becomes `(2)`, and
  `Scalar("0.5", value=0.5)` is written `(0.5)` like `Number(0.5)` (`x0.5` → `x*(0.5)`).
- The variation and the derivative of a `Constant` vector or matrix are `ZeroVector`
  (`0v`) and `ZeroMatrix` (`0`) instead of a vector `0` and a matrix `O`.
- `display_standard_form` writes an input key such as `R*e3*f` as `R e3 f`.
- Long functions are split by node class: `_eliminate` (`_ELIMINATE` table), `pull`
  (`_PULL`), `to_latex` (`_latex_handlers()`), `collect`'s `_collect_dot`, `linearize`
  (`_kinematic_row`, `_dynamic_row`) and `compute_eom`'s extraction message
  (`_higher_order_error`, `_ibp_left_error`, `_collect_error`); `compute_eom` passes the
  extraction context with a small internal exception instead of an attribute set on the
  error. Rule names are unchanged.
- Pinned equation strings change form (phase 5a; every changed output was checked to be
  numerically identical, ≤ 1e-15 relative, at fixed states, see
  `specs/reviews/phase5a-string-changes.md`):
  - the numbers of a product are folded into one placed last, `(-1)mll → mll(-1)`,
    `Mul(-1, Dot(…)) → Mul(Dot(…), -1)`;
  - a float with an integer value is an `int`: `J\frac{d}{dt}(\Omega_{R})(-1.0)` is
    now `…(-1)`, `(2.0)` is `(2)`; one equation no longer mixes `(-1)` and `(-1.0)`;
  - terms that are identically zero on the manifold or like terms that were kept apart
    are gone (the double spherical pendulum's two equations shrink from 14 and 12 terms
    to 10 and 8).
- `Scalar("m", value=2)` stays a symbol in every simplification and in `to_latex`
  (decision D1): before, `m + l` with values was folded to `(5)` and the names lost. Only
  `geomech.numeric` uses the value.
- `linearize` on S²: blocks of an equation that changed form may differ as matrices from
  the old ones by terms that vanish along admissible perturbations (`ξ ⊥ q`,
  `q·δω = −ω·(ξ × q)`); every row `−M δẋ + A δx + B δu` is unchanged along such
  directions (the blocks were never unique on S²). `MVMul(Transpose(Hat(a) + Hat(b)), w)`
  and `Vee(Hat(a + b))` in an equation now linearize instead of raising.
- Rule trace lines: a firing whose two sides have the same `repr` but a different `str`
  shows the `str` forms, a side longer than 300 characters is cut around the first
  difference, and a line whose two sides would print the same is not logged. `repr`
  parenthesises a negated operand (`-(-hat(v:ρ*))`, `hat(v:ρ*)·(-(-1))`).
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
- `is_leaf` and `has_nested_add` (unused by the engine; they were exported from `geomech`
  and `geomech.core.transformations`) and the unreachable rule `vector.transpose_transpose`
  (the `Transpose` constructor already collapses a double transpose). Phase 5b.
- `UndefinedCaseError` (unreachable: every node type is SCALAR, VECTOR or MATRIX; it
  carried no context) is removed from `geomech.utils.errors`. It was not exported from
  `geomech`.

### Fixed
- A rewrite under `TimeDerivative`/`Variation` that reduces the operand to a bare
  `S2`/`SO3` leaf (`d/dt(q·(q·q))`, `d/dt(q·1)`, `δ(I q)`, `d/dt(RᵀR q)`) raised
  `ValueError` from `simplify`/`full_simplify`; the node is now rebuilt in eager form,
  `q.t_diff()` or `q.delta()` (`transformations.utils.calculus_node`).
- The scalar factor of a vector or matrix term is normalised like a scalar chain:
  `v·(2.0)` simplifies to `v(2)`, the same string as `v·(2)`.
- `Number(10**400) * x` raised `OverflowError` in `simplify`; an `int` literal is
  compared exactly.
- `repr` of a named scalar with a value prints it as a constant symbol (`m*`), not as a
  bare number.
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
