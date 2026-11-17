# Manifold rules

The kinematic substitutions of pipeline stage 2 ([Pipeline](pipeline.md)), why the $S^2$
rule has a factor 2, how $S^1$ is embedded in $\mathbb{R}^3$, and how the manifold points
are kept out of the tree.

## S¹ variation rule: δω = γ̇

$S^1$ is embedded as the unit circle in the plane normal to a constant unit axis $e$:
`S1("q", axis=e)` with $q \cdot e = 0$, $\lvert q \rvert = 1$. The engine is typed in three
dimensions (`Cross`, `Hat`), and in that plane $e \times \cdot$ acts exactly as the rotation
$S = \begin{bmatrix} 0 & -1 \\ 1 & 0 \end{bmatrix}$ of Lee et al. ch. 4:
$e \times (e \times q) = -q$, $\lvert e \times q \rvert = 1$ and
$(e \times a) \cdot b = -a \cdot (e \times b)$, the skew-symmetry of $S$. So
$\dot q = \omega\, e \times q$ is eq. 4.1 and $\delta q = \gamma\, e \times q$ is eq. 4.10,
with scalar $\omega$ and $\gamma$.

### Derivation

With $\omega = (e \times q) \cdot \dot q$ and $\dot e = 0$,

$$
\begin{aligned}
\delta\dot q &= \frac{d}{dt}(\delta q) = \dot\gamma\, e \times q + \gamma\, e \times \dot q = \dot\gamma\, e \times q - \gamma \omega q, \\
\delta\omega &= (e \times \delta q) \cdot \dot q + (e \times q) \cdot \delta\dot q \\
&= \gamma\, (e \times (e \times q)) \cdot \dot q + \dot\gamma\, \lvert e \times q \rvert^2 - \gamma \omega\, (e \times q) \cdot q \\
&= -\gamma\, q \cdot \dot q + \dot\gamma = \dot\gamma,
\end{aligned}
$$

because $q \cdot \dot q = 0$ and $(e \times q) \cdot q = 0$. This is eq. 4.11. The circle is
abelian, so, unlike $S^2$ and $SO(3)$, $\delta\omega$ has no term in $\omega$, and being a
scalar it has no normal component either. `manifold_rules.py` substitutes
$\delta\omega \to \dot\gamma$ (the `Variation` node of $\omega$ by the `TimeDerivative` node
of $\gamma$, `manifold.s1_variation_omega`); integration by parts then treats $\dot\gamma$
like the rate of a scalar configuration variable, and each $S^1$ point gives the scalar
equation

$$
E = -\frac{d}{dt}\left(\frac{\partial L}{\partial \omega}\right) + (e \times q) \cdot \frac{\partial L}{\partial q} + \text{forcing},
$$

eq. 4.12 with every term moved to one side and the library's sign.

### Identities

The rules of `vector_identities.py` that use the $S^1$ structure, each checking that $q$ is
an `S1` point whose axis equals $e$ by name (the axis is a constant `UnitNorm` leaf, which
the `S1` constructor requires):

- `vector.s1_axis_orthogonal`: $e \cdot q = q \cdot e = 0$;
- `vector.s1_axis_double_cross`: $e \times (e \times q) = -q$ and its three mirrors (it
  reduces $\frac{d}{dt}(e \times q) = \omega\, e \times (e \times q)$ to $-\omega q$);
- `vector.s1_axis_recovery`: $q \times (e \times q) = e$ and its three mirrors;
- `vector.s1_rotated_dot`: $(e \times a) \cdot (e \times b) = a \cdot b$ for two points of
  one circle (Lagrange's identity with $\lvert e \rvert = 1$, $a \cdot e = b \cdot e = 0$);
  $a = b$ gives $\lvert e \times q \rvert = 1$, so a kinetic energy written as
  $\tfrac12 m l^2\, \dot q \cdot \dot q$ reduces to $\tfrac12 m l^2 \omega^2$;
- `vector.s1_planar_triple`: the triple product $e \cdot (a \times b)$ of two points of one
  circle, the book's $a^T S b$, is written in one form, `Dot(p, Cross(e, r))` with
  `str(p) < str(r)`. $(e \times a) \cdot b = -a \cdot (e \times b)$ holds for any vectors;
  it is applied to $S^1$ points only, so that the two coupling terms of the double planar
  pendulum (from $\frac{d}{dt} \partial\tilde L/\partial\omega_1$ and from
  $q_1^T S\, \partial\tilde L/\partial q_1$) are recognised as like terms and cancel.

`scripts/engine_fuzz.py` samples each $S^1$ point as a curve on its circle to all orders
($q = P p / \lVert P p \rVert$, $P = I - e e^T$, $\omega = (e \times q) \cdot \dot q$,
$\gamma = (e \times q) \cdot \partial q / \partial \varepsilon$), so every rewrite above, and
the manifold rule, is checked numerically on random trees; a point of another circle, a
generic vector and a unit vector that is not the axis are near misses the rules must leave
alone.

### Errata

The $(q, \omega)$ form of the double planar pendulum, Lee et al. eq. 4.62, prints the
velocity terms $m_2 L_1 L_2\, \omega_2 (\omega_2 - \omega_1)\, q_1^T S q_2$ and
$m_2 L_1 L_2\, \omega_1 (\omega_1 - \omega_2)\, q_2^T S q_1$. That is
$\frac{d}{dt} \partial\tilde L/\partial\omega_i$ alone; eq. 4.12 adds
$q_i^T S\, \partial\tilde L/\partial q_i$, whose kinetic part is
$m_2 L_1 L_2\, \omega_1 \omega_2\, q_i^T S q_j$, so the equations have $\omega_2^2$ and
$\omega_1^2$, those of the textbook double pendulum in angles.
`tests/test_systems_book.py` compares the derived equations with the corrected eq. 4.62 and
with the angle form, and checks that the printed form differs.

## S² variation rule: δω = ξ̇ − 2 ω × ξ

`manifold_rules.py` substitutes the `Variation` of an $S^2$ tangent vector $\omega$ by

$$
\delta\omega = \dot\xi - 2\, \omega \times \xi
$$

### Derivation

The kinematics are $\dot q = \omega \times q$ with $\omega \perp q$ (equivalently
$\omega = q \times \dot q$) and the variation is $\delta q = \xi \times q$ with $\xi \perp q$.
Lee et al. eq. (5.13) give

$$
\delta\omega = -\widehat{\omega}\, \xi + (I - q q^T)\, \dot\xi = \dot\xi - \omega \times \xi - q\, (q \cdot \dot\xi).
$$

The last term does not vanish. Differentiating the constraint $\xi \cdot q = 0$ along the
curve,

$$
q \cdot \dot\xi = -\xi \cdot \dot q = -\xi \cdot (\omega \times q) = q \cdot (\omega \times \xi),
$$

and since $\omega$ and $\xi$ both lie in the tangent plane, $\omega \times \xi$ is parallel
to $q$, so $q\, (q \cdot \dot\xi) = q\, (q \cdot (\omega \times \xi)) = \omega \times \xi$.
On the constraint set therefore

$$
\delta\omega = \dot\xi - 2\, \omega \times \xi.
$$

$\delta\omega$ is **not** tangent to the sphere:
$q \cdot \delta\omega = -q \cdot \dot\xi = -q \cdot (\omega \times \xi)$, which is what keeps
$\delta(\omega \cdot q) = 0$.

### History and symptom

Until the 2026-11 audits the rule was $\delta\omega = \dot\xi - \omega \times \xi$, which is
only the tangential projection $(I - q q^T)\, \delta\omega$. The two agree whenever
$\partial L/\partial\omega \perp q$, which is the case for every Lagrangian that depends on
$\omega$ only through $\dot q = \omega \times q$ (all the pendulum examples), so the pinned
pendulum equations never changed. The projection is wrong as soon as
$\partial L/\partial\omega$ has a component along $q$: with
$L = \tfrac12 m l^2\, \omega \cdot \omega + c \cdot \omega - m g l\, q \cdot e_3$ the true
equation is $-m l^2 \dot\omega + 2\, \omega \times c - m g l\, q \times e_3 + f$ and the old
rule printed $\omega \times c$, i.e. **gyroscopic terms came out at half strength**
(likewise for $\tfrac12 \omega \cdot A\omega$, a rotating frame
$\lvert \omega + W_0 \rvert^2$, $\mathbb{R}^3 \times S^2$ couplings $\dot x \cdot \omega$ and
$S^2 \times SO(3)$ couplings $\omega \cdot \Omega$, $\omega \cdot (R^T e_3)$). The
regression tests are in `tests/test_s2_variation_fix.py`; each pinned equation is also
checked with the finite-difference Euler–Lagrange check in `geomech.numeric`.

### SO(3) variation rule

For $\dot R = R \widehat{\Omega}$ and $\delta R = R \widehat{\eta}$, equating the mixed
derivatives $\frac{d}{dt}(\delta R) = \delta(\dot R)$ gives
$\delta\Omega = \dot\eta + \Omega \times \eta$ exactly (Lee et al. eq. 6.9); it has no
hidden normal component.

---

## How the S² substitution is organised

($S^1$ is organised the same way: `S1.t_diff()` eagerly returns $(e \times q)\, \omega$,
`S1.delta()` returns $(e \times q)\, \gamma$, and the only $S^1$ substitution is
$\delta\omega \to \dot\gamma$.)

- `S2.t_diff()` eagerly returns the `Cross` node $\omega \times q$ — the tree never holds `TimeDerivative(q)`.
- `S2.delta()` eagerly returns the `Cross` node $\xi \times q$ — the tree never holds `Variation(q)`.
- `manifold_rules` substitutes $\delta\omega \to \dot\xi - 2\, \omega \times \xi$ and
  $\delta\Omega \to \dot\eta + \Omega \times \eta$. This is the only place these formulas
  live; `TS2`/`TSO3` carry no substitution of their own.

An equivalent formulation works at the level of $q$ (substituting $\delta(\dot q)$,
$\dot q$ and $\delta q$ directly and never forming $\delta\omega$). The two are related by
the Jacobi identity; the $\omega$-level formulation needs the extra simplification rules
below (BAC-CAB, skew-symmetry) because it produces more intermediate cross products.

### Rules that rely on the manifold structure

The rewrite rules are listed, with the hypothesis each one checks, in the docstrings of
their modules; the ones that depend on $S^2$ or $SO(3)$ are:

- `vector.dot_tangent_orthogonal` ($\omega \cdot q = 0$), `vector.dot_unit_norm`
  ($q \cdot q = 1$), `vector.s2_tangent_recovery` ($q \times (a \times q) = a$),
  `vector.s2_tangent_negated` ($q \times (q \times a) = -a$) and `vector.s2_double_cross`
  ($a \times (a \times q) = -q\, (a \cdot a)$) in `vector_identities.py`. Tangency is only
  assumed for a `TS2` leaf of the same point (compared by name) and for `Cross(…, q)`;
  never for $\dot\omega$, $\dot\xi$ or the substituted $\delta\omega$, which is correct
  because $\delta\omega$ is not tangent. `vector_rules` recurses through `Add`, `Mul`,
  `VAdd`, `SVMul` and `Cross` only.
- `simplify.orthogonal_transpose` ($R^T R = R R^T = I$ for an `SO3` point),
  `simplify.mvmul_hat_to_cross` ($\widehat{v}\, w = v \times w$),
  `simplify.transpose_hat_skew` ($\widehat{v}^T = -\widehat{v}$), `simplify.hat_vee_skew`
  (`Hat(Vee(X)) = X` for skew-symmetric $X$ only) and `simplify.vee_hat` in
  `simplification.py`.
- `expand.hat_over_vadd` ($\widehat{a + b} = \widehat{a} + \widehat{b}$) in `expansion.py`.
- `collect.strip_rotation` (`Dot(R a, R b) = Dot(a, b)`), `collect.peel_matrix`
  ($(M f(v)) \cdot w = f(v) \cdot (M^T w)$, one matrix layer per step, only when the target
  is in the vector operand) and `collect.triple_product_*` (one name per rotation) in
  `collecting.py`.

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
  standard form, VBL, codegen, simulation, LQR, and the spacecraft on $SO(3) \times (S^1)^3$
- `tests/test_manifold_pipelines.py` — end-to-end pipelines against Lee et al.
- `tests/test_s2_variation_fix.py` — systems with $\partial L/\partial\omega$ not orthogonal
  to $q$
- `tests/test_collect_fix.py` — velocities nested under matrices inside a `Dot`
- `tests/test_eom_oracle.py` — every pinned system against the numerical
  Euler–Lagrange check
