# Spec: variation-based linearization (VBL)

Branch: `feature/vbl` (on top of `feature/scalar-implementations`). Roadmap feature 4.
Primary references: Wu & Sreenath, "Variation-Based Linearization of Nonlinear Systems
Evolving on SO(3) and S²", IEEE Access 3, 2015; Lee et al. Appendix B (linearization on
TS², TSO(3)), §5.3 (S² pendulum, Eq. 5.14), §6 (rigid body, Eqs. 6.8, 6.16).

## Goal

Given the equations of motion produced by `compute_eom`, produce their linearization along
an arbitrary trajectory in the first-order form

    M δẋ = A δx + B δu

where the perturbation state δx is built from manifold-respecting perturbation variables
(ξ on S², η on SO(3), δx on R³, δs on R, and the variations of the corresponding
velocities), not from local coordinates. The result is a set of symbolic block matrices
that are functions of the trajectory (q, v, v̇, u). This is the formulation of Wu &
Sreenath 2015 and of the quadrotor-load papers from the Hybrid Robotics lab; it is what
LQR and variation-based controllers consume.

## Mathematics

### Perturbation variables

For each configuration variable the perturbation state consists of a configuration
perturbation and a velocity perturbation. The velocity perturbation of every variable is
the `Variation` of its velocity symbol and is **not** substituted by the EOM-pipeline
manifold rules (`apply_manifold_rules` must not be called in VBL).

| Configuration | Config. perturbation | Velocity perturbation | d/dt of each |
|---|---|---|---|
| `s: Scalar` (R) | δs = `Variation(s)` | δṡ = `TimeDerivative(Variation(s))` | δṡ, δs̈ |
| `x: Vector` (R³) | δx = `Variation(x)` | δẋ = `TimeDerivative(Variation(x))` | δẋ, δẍ |
| `q: S2` | ξ = `q.get_variation_vector()` (TS2, ξ ⊥ q) | δω = `Variation(ω)`, ω = `q.get_tangent_vector()` | ξ̇, δω̇ |
| `R: SO3` | η = `R.get_variation_vector()` (TSO3) | δΩ = `Variation(Ω)`, Ω = `R.get_tangent_vector()` | η̇, δΩ̇ |

The engine already canonicalises δ(v̇) as `TimeDerivative(Variation(v))` (δ and d/dt
commute on vector-space quantities), and `S2.delta()`, `SO3.delta()` are eager:
δq = ξ × q, δR = R·Hat(η). δ(constant) = 0. δu = `Variation(u)` for each input
(scalar or vector).

### Kinematic rows (one per configuration variable)

Derived from commutativity of δ and d/dt applied to the manifold kinematics.

- R and R³: d/dt(δs) = δṡ and d/dt(δx) = δẋ. These are identities between the same
  symbols; the row is included so that every configuration variable has a kinematic row
  and the block structure is uniform.
- S² (q̇ = ω × q, ω ⊥ q, δq = ξ × q, ξ ⊥ q). Equating δ(q̇) = δω × q + ω × (ξ × q) with
  d/dt(δq) = ξ̇ × q + ξ × (ω × q) and using Jacobi gives (δω − ξ̇ − ξ × ω) × q = 0, so
  δω − ξ̇ − ξ × ω = c q. The tangency constraint d/dt(ξ·q) = 0 fixes qᵀξ̇ = qᵀ ω̂ ξ, and
  (I − qqᵀ)(ω × ξ) = 0 because ω × ξ ∥ q. Hence

      ξ̇ = q qᵀ ω̂ ξ + (I − q qᵀ) δω.

  In the M δẋ = A δx + B δu form: M[ξ̇] = I, A[ξ] = q qᵀ Hat(ω), A[δω] = I − q qᵀ.
- SO(3) (Ṙ = R·Hat(Ω), δR = R·Hat(η)). The EOM manifold rule δΩ = η̇ + Ω × η
  (`manifold_rules.py` docstring) rearranged:

      η̇ = −Hat(Ω) η + δΩ.

  M[η̇] = I, A[η] = −Hat(Ω), A[δΩ] = I.

### Dynamic rows (one per equation from `compute_eom`)

`compute_eom` returns, for each variation vector, an expression E(q, v, v̇, u) with E = 0.
By the library convention E is the coefficient of the variation in the post-IBP integrand,
i.e. the negative of the Euler–Lagrange expression; `to_standard_form` therefore reports
M = −(physical inertia). VBL takes the first variation δE (which is linear in the
perturbation variables) and defines the blocks so that the physical inertia comes out
positive:

    δE = −M δv̇ + A δx + B δu

    M[row][δv̇_j] = −coeff(δv̇_j, δE)
    A[row][δx_j] =  coeff(δx_j, δE)       (δx_j ranges over config. and velocity perturbations)
    B[row][δu_j] =  coeff(δu_j, δE)

Check (free rigid body, Lee et al. Eq. 6.16): E = −JΩ̇ − Ω × JΩ + M_in. δE =
−J δΩ̇ − δΩ × JΩ − Ω × J δΩ + δM = −J δΩ̇ + (Hat(JΩ) − Hat(Ω) J) δΩ + δM. So
M = J, A[δΩ] = Hat(JΩ) − Hat(Ω) J, B = I, which is the textbook linearization.

Identities the extraction must handle (all already supported by `extract_linear_coeff`,
verified by probe): coeff(η, η × ρ) = −Hat(ρ); coeff(η, ρ × η) = Hat(ρ);
coeff(η, Hat(η) ρ) = −Hat(ρ); coeff(δΩ, δΩ × JΩ) = −Hat(JΩ); coeff(δΩ, Ω × J δΩ) =
Hat(Ω) J; δ(Rᵀ e₃) = −η × Rᵀe₃ = Hat(Rᵀe₃) η; coeff of a vector perturbation in a scalar
equation is a row vector `Transpose(w)`.

### Completeness

After extraction, every additive term of `expand(δE)` must contain exactly one
perturbation symbol (configuration, velocity, d/dt of velocity, or input variation). A term
with none, or with a perturbation symbol only inside a non-canonical calculus node, means
the linearization is incomplete and `linearize` must raise `AlgebraicError` naming the
term. No silent zeros.

### Simplification of blocks

Each block is passed through `full_simplify` (as `to_standard_form` does). The probe
showed `Hat((R)'e3)(-1)(-1)` surviving simplification; numeric factors inside a `Mul`
chain must be folded so that blocks print without `(-1)(-1)`. If that needs a new rule in
`simplify.py`, add it with its own unit tests, and do not alter any existing test's
expectation.

## Public API

```python
from geomech.dynamics import linearize, LinearizedSystem

lin = linearize(eom, variables, inputs=None)
```

- `eom`: the dict returned by `compute_eom` (`{name: (variation_vec, eq)}`).
- `variables`: the `SystemVariables` used for `compute_eom`. Validated with the same
  `_check_variables` as `compute_eom`.
- `inputs`: list of `Scalar` / `Vector` inputs (same meaning as in `to_standard_form`).
  Scalar inputs (quadrotor thrust) must work.

`LinearizedSystem` is a frozen dataclass:

- `rows`: ordered row labels. Kinematic rows first, one per configuration variable in
  `SystemVariables` order (scalars, vectors, matrices), labelled `"kin:" + str(var)`;
  then dynamic rows labelled by the `eom` keys.
- `states`: ordered list of perturbation symbols for δx: all configuration perturbations
  in `SystemVariables` order, then all velocity perturbations in the same order.
- `dstates`: `[s.t_diff() for s in states]`.
- `inputs`: `[u.delta() for u in inputs]`.
- `M`, `A`, `B`: `dict[str, dict[str, Expr]]`, `M[row][str(dstate)]`,
  `A[row][str(state)]`, `B[row][str(input_variation)]`. Dense: every (row, col) present.
  Block shape follows row type × column type: vector row × vector col → 3×3 Matrix;
  scalar row × vector col → 1×3 (`Transpose(Vector)`); vector row × scalar col → 3×1
  Vector; scalar × scalar → Scalar. Zero blocks use the existing zero constants of the
  right shape.
- `block_shape(row, col) -> tuple[int, int]` helper.
- `__str__` lists the non-zero blocks row by row. No LaTeX beyond `to_latex(entry)`,
  which already works on any Expr.

Also export `linearize` and `LinearizedSystem` from `geomech/__init__.py` next to
`to_standard_form`.

## Systems that must be reproduced

Each system is defined once in `tests/vbl_systems.py` (functions returning
`(L, dW, variables, inputs)` plus the symbols needed by the tests) so the adversarial
reviewer can reuse them. Expected blocks below are hand derivations from the library's
pinned EOM in `tests/test_manifold_pipelines.py`; the implementer re-derives each before
writing the test and records the derivation in the test docstring.

1. **Point mass R³** (`m ẍ = F`): states [δx, δẋ]; M = blkdiag(I, m I);
   A = [[0, I], [0, 0]]; B = [[0], [I]]. Exact string/structure test.
2. **Spring-mass R** (`m s̈ + k s = 0`, pinned in `TestSpringMassPipeline`): states
   [δs, δṡ]; M = diag(1, m); A = [[0, 1], [−k, 0]]. Exact.
3. **Free rigid body SO(3)** (Lee et al. Eq. 6.16): states [η, δΩ];
   M = blkdiag(I, J); A = [[−Hat(Ω), I], [0, Hat(JΩ) − Hat(Ω) J]]; B = [[0], [I]].
   Exact up to printed form; numeric equality mandatory.
4. **Rigid pendulum SO(3) about the pivot** (Lee et al. Eq. 6.8, pinned in
   `TestRigidPendulumPivotPipeline`): the gravity block follows from
   δ(ρ × Rᵀe₃) = Hat(ρ) Hat(Rᵀe₃) η; the sign comes from the pinned f. Numeric equality
   against the hand formula at random states.
5. **Spherical pendulum S², ω-form** (Lee et al. Eq. 5.14, pinned in
   `TestSphericalPendulumPipeline`: E = −m l² ω̇ − m g l q × e₃ + f): states [ξ, δω];
   M = blkdiag(I, m l² I); A = [[q qᵀ Hat(ω), I − q qᵀ], [−m g l Hat(e₃) Hat(q), 0]];
   B = [[0], [I]]. Sanity check to include in the docstring: at q = −e₃, ω = 0 the ξ
   block reduces to m l² ξ̈ = −m g l ξ on the tangent plane, the pendulum frequency.
6. **Quadrotor R³ × SO(3)** (Lee, Leok, McClamroch CDC 2010; Wu & Sreenath 2015):
   L = ½ m ẋ·ẋ + ½ Ω·JΩ − m g e₃·x, dW = f (R e₃)·δx + M·η, inputs [f (Scalar), M].
   Translational row: M = m I, A[η] = −f R Hat(e₃) (from δ(f R e₃) = δf R e₃ +
   f R Hat(η) e₃ = δf R e₃ − f R Hat(e₃) η), B[δf] = R e₃. Rotational row as system 3.
   Numeric equality against the hand formula.
7. **Double rigid pendulum SO(3) × SO(3)** (system of
   `TestDoubleRigidPendulumMassMatrixExact`): finite-difference check only, plus exact
   kinematic rows.
8. **Elastic spherical pendulum S² × R** (from `tests/test_scalar_adversarial.py`):
   finite-difference check only.
9. **Quadrotor with cable-suspended point-mass load R³ × S² × SO(3)** (Sreenath, Lee,
   Kumar CDC 2013): load position x_L, cable direction q (from load to quadrotor, length
   l), x_Q = x_L − l q; L = ½ m_Q ẋ_Q·ẋ_Q + ½ m_L ẋ_L·ẋ_L + ½ Ω·JΩ − m_Q g e₃·x_Q −
   m_L g e₃·x_L; dW = f (R e₃)·(δx_L − l ξ × q) + M·η. Finite-difference check. If
   `compute_eom` itself fails on this system, do not patch the engine: record the failure
   in the implementation report and mark the test `xfail(strict=True)` with the reason.

### Finite-difference check (mandatory for every system above)

Promote the numeric evaluator from
`specs/reviews/scripts/double_rigid_pendulum_numeric_check.py` (`hat`, `vee`,
`rand_rot`, `Evaluator`) into `tests/numeric_helpers.py`, adding `expm_so3(v)`
(Rodrigues) and the perturbation below. Leave the review script untouched.

For a base state sampled on the manifold (|q| = 1, q·ω = 0, R ∈ SO(3)) and a random
direction (δx, δv, δv̇, δu) tangent to the constraints, define the perturbed state

- R, R³: s(ε) = s + ε δs, x(ε) = x + ε δx, velocities and accelerations likewise;
- SO(3): R(ε) = R·expm(ε Hat(η)), Ω(ε) = Ω + ε δΩ;
- S²: q(ε) = expm(ε Hat(ξ)) q with ξ ⊥ q, ω(ε) = expm(ε Hat(ξ)) (ω + ε w) with w ⊥ q, so
  that |q(ε)| = 1 and q(ε)·ω(ε) = 0 hold exactly; the induced velocity direction is
  δω = ξ × ω + w;
- inputs: u(ε) = u + ε δu.

Then with ε = 1e-5 and central differences, for every dynamic row

    (E(ε) − E(−ε)) / (2ε)  ≈  −M δv̇ + A δx + B δu

to relative tolerance 1e-6, at 5 random base states and directions, where the right-hand
side is the numeric evaluation of the symbolic blocks summed over columns. This check is
independent of the symbolic rewrite rules and is the primary evidence of correctness.

## Non-goals

- Reduction to minimal coordinates at an equilibrium (Lee et al. Appendix B basis σ);
  the output stays in the full perturbation vectors as in Wu & Sreenath.
- LQR, controller synthesis, discrete-time VBL, code generation, simulation.
- Constraints between bodies, closed chains.
- Structural equality or canonical form of expressions (roadmap feature 2); numeric
  equality in tests is the accepted substitute for exact-string tests where the printed
  form is long.
- Any change to `compute_eom`, `to_standard_form`, or existing tests.

## Definition of done

- `geomech/dynamics/vbl.py` with `linearize` and `LinearizedSystem`; module docstring
  carries the mathematics above (kinematic rows, sign convention, completeness rule).
- Exports in `geomech/dynamics/__init__.py` and `geomech/__init__.py`.
- `tests/numeric_helpers.py`, `tests/vbl_systems.py`, `tests/test_vbl.py` covering
  systems 1–9, the completeness error path, input validation, scalar inputs, and the
  kinematic rows exactly.
- `docs/vbl.md` (user-facing: what VBL is, the perturbation variables, the sign
  convention, a worked free-rigid-body example, how to read the blocks) and a nav entry
  in `mkdocs.yml`.
- `CHANGELOG.md` entry under Unreleased / Added. README's first line made accurate.
- Full suite green, no existing test modified, `ruff format --check` and `ruff check`
  clean.
- `specs/reviews/vbl-implementation-report.md`: what was built, what was derived by
  hand, anything that could not be done and why, engine changes made (if any).
