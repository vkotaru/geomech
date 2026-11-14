# Changelog

All notable changes to `geomech` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
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

### Fixed
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
