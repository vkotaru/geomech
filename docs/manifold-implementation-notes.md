# Manifold rules

The kinematic substitutions of pipeline stage 2 ([Pipeline](pipeline.md)), why the S²
rule has a factor 2, how S¹ is embedded in R³, and how the manifold points are kept out
of the tree.

## S¹ variation rule: δω = γ̇

S¹ is embedded as the unit circle in the plane normal to a constant unit axis `e`:
`S1("q", axis=e)` with `q·e = 0`, `|q| = 1`. The engine is typed in three dimensions
(`Cross`, `Hat`), and in that plane `e × ·` acts exactly as the rotation
`S = [[0, −1], [1, 0]]` of Lee et al. ch. 4: `e × (e × q) = −q`, `|e × q| = 1` and
`(e × a)·b = −a·(e × b)`, the skew-symmetry of `S`. So `q̇ = ω e × q` is eq. 4.1 and
`δq = γ e × q` is eq. 4.10, with scalar `ω` and `γ`.

### Derivation

With `ω = (e × q)·q̇` and `ė = 0`,

```
δq̇ = d/dt(δq) = γ̇ e × q + γ e × q̇ = γ̇ e × q − γ ω q,
δω = (e × δq)·q̇ + (e × q)·δq̇
   = γ (e × (e × q))·q̇ + γ̇ |e × q|² − γ ω (e × q)·q
   = −γ q·q̇ + γ̇ = γ̇,
```

because `q·q̇ = 0` and `(e × q)·q = 0`. This is eq. 4.11. The circle is abelian, so,
unlike S² and SO(3), `δω` has no term in `ω`, and being a scalar it has no normal
component either. `manifold_rules.py` substitutes `Variation(ω) → TimeDerivative(γ)`
(`manifold.s1_variation_omega`); integration by parts then treats `γ̇` like the rate of a
scalar configuration variable, and each S¹ point gives the scalar equation

```
E = −d/dt(∂L/∂ω) + (e × q)·∂L/∂q + forcing,
```

eq. 4.12 with every term moved to one side and the library's sign.

### Identities

The rules of `vector_identities.py` that use the S¹ structure, each checking that `q` is
an `S1` point whose axis equals `e` by name (the axis is a constant `UnitNorm` leaf, which
the `S1` constructor requires):

- `vector.s1_axis_orthogonal`: `e·q = q·e = 0`;
- `vector.s1_axis_double_cross`: `e × (e × q) = −q` and its three mirrors (it reduces
  `d/dt(e × q) = ω e × (e × q)` to `−ω q`);
- `vector.s1_axis_recovery`: `q × (e × q) = e` and its three mirrors;
- `vector.s1_rotated_dot`: `(e × a)·(e × b) = a·b` for two points of one circle
  (Lagrange's identity with `|e| = 1`, `a·e = b·e = 0`); `a = b` gives `|e × q| = 1`, so a
  kinetic energy written as `½ m l² q̇·q̇` reduces to `½ m l² ω²`;
- `vector.s1_planar_triple`: the triple product `e·(a × b)` of two points of one circle,
  the book's `aᵀ S b`, is written in one form, `Dot(p, Cross(e, r))` with
  `str(p) < str(r)`. `(e × a)·b = −a·(e × b)` holds for any vectors; it is applied to S¹
  points only, so that the two coupling terms of the double planar pendulum (from
  `d/dt ∂L̃/∂ω₁` and from `q₁ᵀS ∂L̃/∂q₁`) are recognised as like terms and cancel.

`scripts/engine_fuzz.py` samples each S¹ point as a curve on its circle to all orders
(`q = P p/|P p|`, `P = I − e eᵀ`, `ω = (e × q)·q̇`, `γ = (e × q)·∂q/∂ε`), so every rewrite
above, and the manifold rule, is checked numerically on random trees; a point of another
circle, a generic vector and a unit vector that is not the axis are near misses the rules
must leave alone.

### Errata

The (q, ω) form of the double planar pendulum, Lee et al. eq. 4.62, prints the velocity
terms `m₂L₁L₂ ω₂(ω₂ − ω₁) q₁ᵀSq₂` and `m₂L₁L₂ ω₁(ω₁ − ω₂) q₂ᵀSq₁`. That is
`d/dt ∂L̃/∂ωᵢ` alone; eq. 4.12 adds `qᵢᵀS ∂L̃/∂qᵢ`, whose kinetic part is
`m₂L₁L₂ ω₁ω₂ qᵢᵀSqⱼ`, so the equations have `ω₂²` and `ω₁²`, those of the textbook double
pendulum in angles. `tests/test_systems_book.py` compares the derived equations with the
corrected eq. 4.62 and with the angle form, and checks that the printed form differs.

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

(S¹ is organised the same way: `S1.t_diff()` eagerly returns `(e × q) ω`, `S1.delta()`
returns `(e × q) γ`, and the only S¹ substitution is `Variation(ω) → γ̇`.)

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
- `geomech/core/base/expressions.py` — S1/TS1, S2/TS2, SO3/TSO3 types
- `geomech/core/transformations/manifold_rules.py` — S1, S2 and SO3 substitutions
- `geomech/core/transformations/vector_identities.py` — S1 and S2 cross product identities
- `geomech/core/transformations/simplification.py` — Hat/Transpose/associativity
- `geomech/core/transformations/expansion.py` — Hat/Vee linearity
- `geomech/core/math/collecting.py` — SO3 rotation stripping, matrix peeling

### Tests
- `tests/test_manifold_rules.py` — S2/SO3 manifold substitutions
- `tests/test_s1_manifold.py` — S1: leaf, rule, identities, Lee et al. §4.6 equations,
  standard form, VBL, codegen, simulation, LQR, and the spacecraft on SO(3) × (S¹)³
- `tests/test_manifold_pipelines.py` — end-to-end pipelines against Lee et al.
- `tests/test_s2_variation_fix.py` — systems with `∂L/∂ω` not orthogonal to q
- `tests/test_collect_fix.py` — velocities nested under matrices inside a `Dot`
- `tests/test_eom_oracle.py` — every pinned system against the numeric Euler–Lagrange
  oracle
