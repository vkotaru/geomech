# Changelog

All notable changes to `geomech` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
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
- Integration by parts dropped `d/dt` of a non-constant scalar factor multiplying
  `Dot(ξ̇, ·)`, so every system with a time-varying scalar in its kinetic energy (elastic
  pendulum, prismatic joint, prescribed `l(t)`) lost its Coriolis / rate-of-inertia terms
  (`−2 m l l̇ ω`, `−m l̈ e₁`, `J̇ Ω`). The prefactor is now differentiated with the product rule.
- `collect` could not move the variation out of `Dot(R₁ a, R₂ b)` with two different rotations,
  which dropped the cross-inertia coupling of multi-body SO(3) systems. The SO(3)×SO(3) double
  rigid pendulum now has `M₁[Ω̇₂] = m₂ l̂₁ R₁ᵀ R₂ ρ̂₂` and `M₂[Ω̇₁]` its transpose.
- `collect` no longer recurses forever on `Dot(M a, N b)` when the target is nested in a form
  the rules cannot expose; the expression is returned unchanged.
- `extract_linear_coeff` compared nodes by rendered string only, so `Mul(k, x)` could be
  mistaken for a leaf named `kx`.

## [0.1.0] — 2026-03

Initial release: expression trees, S² and SO(3) manifold support, `compute_eom`,
`to_standard_form`, printing utilities.
