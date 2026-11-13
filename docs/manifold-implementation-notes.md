# Manifold Implementation Notes

## S² variation rule: δω = ξ̇ − 2 ω × ξ

`manifold_rules.py` substitutes `Variation(ω)` for an S² tangent vector ω by

```
δω = ξ̇ − 2 ω × ξ
```

### Derivation

The kinematics are `q̇ = ω × q` with `ω ⊥ q` (equivalently `ω = q × q̇`) and the
variation is `δq = ξ × q` with `ξ ⊥ q`. Lee et al. eq. (5.13) give

```
δω = −S(ω) ξ + (I − q qᵀ) ξ̇ = ξ̇ − ω × ξ − q (q · ξ̇).
```

The last term does not vanish. Differentiating the constraint `ξ · q = 0` along the curve,

```
q · ξ̇ = −ξ · q̇ = −ξ · (ω × q) = q · (ω × ξ),
```

and since ω and ξ both lie in the tangent plane, `ω × ξ` is parallel to q, so
`q (q · ξ̇) = q (q · (ω × ξ)) = ω × ξ`. On the constraint set therefore

```
δω = ξ̇ − 2 ω × ξ.
```

δω is **not** tangent to the sphere: `q · δω = −q · ξ̇ = −q · (ω × ξ)`, which is what keeps
`δ(ω · q) = 0`.

### History and symptom

Until the 2026-11 audits the rule was `δω = ξ̇ − ω × ξ`, which is only the tangential
projection `(I − q qᵀ) δω`. The two agree whenever `∂L/∂ω ⊥ q`, which is the case for every
Lagrangian that depends on ω only through `q̇ = ω × q` (all the pendulum examples), so the
pinned pendulum equations never changed. The projection is wrong as soon as `∂L/∂ω` has a
component along q: with `L = ½ m l² ω·ω + c·ω − m g l q·e₃` the true equation is
`−m l² ω̇ + 2 ω × c − m g l q × e₃ + f` and the old rule printed `ω × c`, i.e. **gyroscopic
terms came out at half strength** (likewise for `½ ω·Aω`, a rotating frame
`|ω + W₀|²`, R³×S² couplings `ẋ·ω` and S²×SO(3) couplings `ω·Ω`, `ω·(Rᵀe₃)`). The
regression tests are in `tests/test_s2_variation_fix.py`; each pinned equation is also
checked against the finite-difference Euler–Lagrange oracle in `geomech.numeric`.

### SO(3) variation rule

For `Ṙ = R·Hat(Ω)` and `δR = R·Hat(η)`, equating the mixed derivatives `d/dt(δR) = δ(Ṙ)`
gives `δΩ = η̇ + Ω × η` exactly (Lee et al. eq. 6.9); it has no hidden normal component.

---

## How the S² substitution is organised

- `S2.t_diff()` eagerly returns `Cross(ω, q)` — the tree never holds `TimeDerivative(q)`.
- `S2.delta()` eagerly returns `Cross(ξ, q)` — the tree never holds `Variation(q)`.
- `manifold_rules` substitutes `Variation(ω)` → `ξ̇ − 2 ω × ξ` and `Variation(Ω)` →
  `η̇ + Ω × η`. This is the only place these formulas live; `TS2`/`TSO3` carry no
  substitution of their own.

An equivalent formulation works at the level of q (substituting `δ(q̇)`, `q̇` and `δq`
directly and never forming δω). The two are related by the Jacobi identity; the ω-level
formulation needs the extra simplification rules below (BAC-CAB, skew-symmetry) because it
produces more intermediate cross products.

### S² simplification rules (`vector_rules.py`)
- BAC-CAB identity: `cross(q, cross(a, q)) = a` when `‖q‖ = 1` and `a` is tangent
- Variants: `cross(q, cross(q, a)) = −a`, `cross(cross(q, a), q) = a`, etc.
- `cross(a, cross(a, q)) = −q · dot(a, a)` when `a` is tangent
- Recursion into Cross children, VAdd, SVMul

Tangency is only assumed for `TS2` leaves (ω, ξ) and for `Cross(·, q)`; it is never
assumed for `d/dt(ω)`, `d/dt(ξ)` or the substituted δω, which is correct because δω is
not tangent.

### Simplification rules (`simplify.py`)
- `Hat(v) * w → Cross(v, w)` in MVMul
- `MMMul(A, B) * v → A * (B * v)` matrix associativity in MVMul
- `Transpose(Hat(v)) → −Hat(v)` (skew-symmetry)
- `Transpose(SMMul(M, s)) → SMMul(Transpose(M), s)`
- `Transpose(MMMul(A, B)) → MMMul(Transpose(B), Transpose(A))`
- `Hat(Vee(X)) → X` only when X is skew-symmetric (a `Hat`, a `SkewSymmMatrix`, or a
  scalar multiple / sum / transpose of such); `Vee(Hat(v)) → v` always

### Expand rules (`expand.py`)
- `Hat(a + b) → Hat(a) + Hat(b)` (linearity)
- `Vee(a + b) → Vee(a) + Vee(b)` (linearity)

### Collect rules (`collect.py`)
- `Dot(R*a, R*b) → Dot(a, b)` for SO3 rotation stripping
- `Dot(Cross(vec, a), b)` — flip and use the scalar triple product
- `Dot(M*f(vec), w) → Dot(f(vec), Mᵀ*w)` — the matrix is peeled off the side holding the
  target, one layer per step, until the target is a direct operand

---

## Files

### Source
- `geomech/core/base/expressions.py` — SO3/TS2/TSO3 types
- `geomech/core/transformations/manifold_rules.py` — S2 + SO3 substitutions
- `geomech/core/transformations/vector_rules.py` — S2 cross product identities
- `geomech/core/transformations/simplify.py` — Hat/Transpose/associativity
- `geomech/core/transformations/expand.py` — Hat/Vee linearity
- `geomech/core/math/collect.py` — SO3 rotation stripping, matrix peeling

### Tests
- `tests/test_manifold_rules.py` — S2/SO3 manifold substitutions
- `tests/test_manifold_pipelines.py` — end-to-end pipelines against Lee et al.
- `tests/test_s2_variation_fix.py` — systems with `∂L/∂ω` not orthogonal to q
- `tests/test_collect_fix.py` — velocities nested under matrices inside a `Dot`
- `tests/test_eom_oracle.py` — every pinned system against the numeric Euler–Lagrange
  oracle
