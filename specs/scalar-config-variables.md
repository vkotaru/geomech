# Spec: scalar configuration variables in the EOM pipeline

Branch: `feature/scalar-implementations` (commits 16b9f57, e06e4c0, 31a0c45 on top of `geomech`).
Status: implemented, awaiting adversarial review.

## Goal

`SystemVariables(scalars=[x, ...])` participates fully in `compute_eom` and
`to_standard_form`, so systems with scalar degrees of freedom (spring-mass, cable length,
elastic pendulum, varying inertia) derive without `NotImplementedError`.

## Mathematics

Let `x ∈ R` be a configuration variable with Lagrangian `L(x, ẋ)` and virtual work
`δW = f δx`. Hamilton's principle with the boundary terms dropped:

```
∫ (δL + δW) dt = 0
δL = (∂L/∂x) δx + (∂L/∂ẋ) δẋ
∫ (∂L/∂ẋ) δẋ dt = −∫ d/dt(∂L/∂ẋ) δx dt        (integration by parts)
⇒ coefficient of δx:  ∂L/∂x − d/dt(∂L/∂ẋ) + f = 0
```

Equivalent to the Euler–Lagrange form `d/dt(∂L/∂ẋ) − ∂L/∂x = f` (Lee et al. ch. 3). The
library's EOM is the coefficient of the variation in the action integrand after IBP, set to
zero. The overall sign of each equation therefore differs from Euler–Lagrange by −1; this is
a convention, not an error, but it must be the same convention as for vector and manifold
variables, and `to_standard_form` must report `M`, `f`, `G` consistently with
`Σ M[a] a + f + Σ G[u] u = 0`.

Mixed systems: for `(x, q) ∈ R × S²` (elastic spherical pendulum, Lee et al. §5.x), the
scalar variation `δx` and the S² variation `ξ` are independent, so the pipeline must
produce one equation per variable with no cross-contamination; the standard form must
produce the correct cross-coupling in `M` and `f`.

## Public API

- `SystemVariables.scalars: list[Scalar]` — now used (previously accepted but raised).
- `compute_eom(L, dW, variables)` returns `{str(Variation(x)): (Variation(x), eq)}` for each
  scalar `x`, in addition to vector/matrix entries.
- `to_standard_form(eom, variables, inputs)` treats `TimeDerivative(TimeDerivative(x))` as
  the acceleration variable for scalar `x`.
- `extract_linear_coeff(expr, query)` supports scalar queries in scalar, vector and matrix
  expressions (3×3 dispatch grid documented in `extract.py`).
- `integrate_by_parts` handles scalar targets inside `Mul` chains via `_remove_factor`.

## Systems that must be reproduced exactly

| System | L | δW | Expected (Euler–Lagrange form) |
|---|---|---|---|
| Spring-mass | ½ m ẋ² − ½ k x² | f δx | m ẍ + k x = f |
| Spring-mass-damper | ½ m ẋ² − ½ k x² | (f − c ẋ) δx | m ẍ + c ẋ + k x = f |
| Double spring-mass | ½ m₁ ẋ₁² + ½ m₂ ẋ₂² − ½ k₁ x₁² − ½ k₂ (x₂ − x₁)² | f₁ δx₁ + f₂ δx₂ | m₁ ẍ₁ + k₁ x₁ − k₂ (x₂ − x₁) = f₁ ; m₂ ẍ₂ + k₂ (x₂ − x₁) = f₂ |
| Elastic spherical pendulum (R × S²) | ½ m (l̇² + l² ω·ω) − m g l (q·e₃) − ½ k (l − l₀)² | f_l δl + ξ·f_q | m l̈ − m l ω·ω + m g q·e₃ + k (l − l₀) = f_l ; m l² ω̇ + 2 m l l̇ ω + m g l (q × e₃) ... = ... (derive by hand, including the sign of the gravity term under the S² convention ξ × q) |

## Non-goals

- Scalar variables appearing nonlinearly in the *acceleration* (e.g. `ẍ²`) — must raise
  `NonLinearError`, never return a partial coefficient.
- Constraints between scalars and manifold variables.
- Codegen / numeric evaluation (separate feature).

## Definition of done

- All four systems above reproduce the expected equations exactly (string match on the
  simplified EOM or on `M`, `f`, `G`), with the equation derived by hand in the test
  docstring.
- Every unsupported pattern raises a typed error with a message naming the pattern.
- No existing test modified. Full suite + ruff green.
