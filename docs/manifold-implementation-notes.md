# Manifold rules

The kinematic substitutions of pipeline stage 2 ([Pipeline](pipeline.md)), why the S²
rule has a factor 2, and how the manifold points are kept out of the tree.

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

### Rules that rely on the manifold structure

The rewrite rules are listed, with the hypothesis each one checks, in the docstrings of
their modules; the ones that depend on S² or SO(3) are:

- `vector.dot_tangent_orthogonal` (`ω·q = 0`), `vector.dot_unit_norm` (`q·q = 1`),
  `vector.s2_tangent_recovery` (`q × (a × q) = a`), `vector.s2_tangent_negated`
  (`q × (q × a) = −a`) and `vector.s2_double_cross` (`a × (a × q) = −q (a·a)`) in
  `vector_identities.py`. Tangency is only assumed for a `TS2` leaf of the same point (compared
  by name) and for `Cross(·, q)`; never for `d/dt(ω)`, `d/dt(ξ)` or the substituted δω,
  which is correct because δω is not tangent. `vector_rules` recurses through `Add`,
  `Mul`, `VAdd`, `SVMul` and `Cross` only.
- `simplify.orthogonal_transpose` (`RᵀR = RRᵀ = I` for an `SO3` point),
  `simplify.mvmul_hat_to_cross` (`Hat(v) w = v × w`), `simplify.transpose_hat_skew`
  (`Hat(v)ᵀ = −Hat(v)`), `simplify.hat_vee_skew` (`Hat(Vee(X)) = X` for skew-symmetric `X`
  only) and `simplify.vee_hat` in `simplification.py`.
- `expand.hat_over_vadd` (`Hat(a + b) = Hat(a) + Hat(b)`) in `expansion.py`.
- `collect.strip_rotation` (`Dot(R a, R b) = Dot(a, b)`), `collect.peel_matrix`
  (`Dot(M f(v), w) = Dot(f(v), Mᵀ w)`, one matrix layer per step, only when the target is
  in the vector operand) and `collect.triple_product_*` (one name per rotation) in `collecting.py`.

`geomech.set_trace(True)` logs each firing ([Pipeline](pipeline.md#which-rule-fired)).

---

## Files

### Source
- `geomech/core/base/expressions.py` — SO3/TS2/TSO3 types
- `geomech/core/transformations/manifold_rules.py` — S2 + SO3 substitutions
- `geomech/core/transformations/vector_identities.py` — S2 cross product identities
- `geomech/core/transformations/simplification.py` — Hat/Transpose/associativity
- `geomech/core/transformations/expansion.py` — Hat/Vee linearity
- `geomech/core/math/collecting.py` — SO3 rotation stripping, matrix peeling

### Tests
- `tests/test_manifold_rules.py` — S2/SO3 manifold substitutions
- `tests/test_manifold_pipelines.py` — end-to-end pipelines against Lee et al.
- `tests/test_s2_variation_fix.py` — systems with `∂L/∂ω` not orthogonal to q
- `tests/test_collect_fix.py` — velocities nested under matrices inside a `Dot`
- `tests/test_eom_oracle.py` — every pinned system against the numeric Euler–Lagrange
  oracle
