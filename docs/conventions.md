# Conventions

The choices below are fixed by the code; each one has surprised someone reading an
equation. Equation numbers refer to Lee et al., *Global Formulations of Lagrangian and
Hamiltonian Dynamics on Manifolds*.

## Signs

**$E = -(\text{Euler–Lagrange}) + \text{forcing}$.** `compute_eom` returns, for each
variation vector, the coefficient $E$ of that variation in the integrand of $\delta S$
after integration by parts, and the equation is $E = 0$. For $L(x, \dot x)$ and
$\delta W = \delta x \cdot f$ this is

$$
E = -\left(\frac{d}{dt}\frac{\partial L}{\partial \dot x} - \frac{\partial L}{\partial x}\right) + f,
$$

so the point mass $L = \tfrac12 m\, \dot x \cdot \dot x$ with $\delta W = \delta x \cdot F$
gives $E = -m \ddot x + F$, and the free rigid body gives
$E = -J \dot\Omega - \Omega \times J\Omega + M$ (eq. 6.16 with every term moved to one side).

**$M = -(\text{physical inertia})$ in `to_standard_form`.** The blocks split the same $E$:
$E = M a + f + G u$, so $M[\ddot x] = -m I$ for the point mass, $M[\dot\Omega] = -J$ for a
rigid body, and $G[F] = I$ for a force entering as $\delta x \cdot F$.

**VBL.** `linearize` takes $\delta E$ and defines its blocks so that the physical inertia is
positive: $\delta E = -M \delta\dot v + A \delta x + B \delta u$, i.e.
$M[\delta\dot v] = -\operatorname{coeff}(\delta\dot v, \delta E)$, and the result reads
$M \delta\dot x = A \delta x + B \delta u$ ([Linearization](vbl.md#sign-convention)).

## S¹ points live on a circle normal to an axis

`S1("q", axis=e)` is the unit circle $\{q : \lvert q \rvert = 1,\ q \cdot e = 0\}$ in
$\mathbb{R}^3$, with $e$ a constant `Vector` leaf flagged `UnitNorm` (required:
`Vector("e", attr=["Constant", "UnitNorm"])`). In that plane $e \times \cdot$ is the
rotation $S$ by $\pi/2$ of Lee et al. ch. 4, so $\dot q = \omega\, e \times q$,
$\delta q = \gamma\, e \times q$ (scalar $\omega$, $\gamma$) and $p^T S r$ is printed
`Dot(p,Cross(e,r))`. Each $S^1$ point has one scalar equation
$E = -\frac{d}{dt}\frac{\partial L}{\partial \omega} + (e \times q) \cdot \frac{\partial L}{\partial q} + \text{forcing}$
(eq. 4.12 in the library's sign), keyed by `\gamma_{q}`, with no projection: unlike $S^2$,
it has no normal component. The axis is a parameter of the generated code (`p[str(e)]`).
Write the Lagrangian in $\omega$ (the modified Lagrangian of the book) or in `q.t_diff()`;
the $S^1$ identities reduce $\dot q \cdot \dot q$ to $\omega^2$.

## S² equations live on the tangent plane

For a point $q \in S^2$ the variation is $\delta q = \xi \times q$ with $\xi \perp q$, so
only the part of $E$ orthogonal to $q$ is an equation: $(I - q q^T) E = 0$. $E$ as derived
in general has a component along $q$ (from the forcing, from $\omega \times c$, from
$A \dot\omega$), and `to_standard_form` and `linearize` return the terms as derived, not
projected. To integrate, solve $(I - q q^T)(M \dot\omega + f + \sum G u) = 0$ together with
$\dot\omega \perp q$; the unprojected $3 \times 3$ solve gives $q \cdot \dot\omega \neq 0$
and a trajectory that leaves the sphere. `numeric.check_eom` compares on the tangent plane
for the same reason.

## Rows and columns

There is one `VECTOR` type; a row is recognised by its structure (`is_row_vector`):
`Transpose(column)`, a row times a matrix $v^T M$ (`MVMul(row, M)`), a scalar multiple of a
row, a sum whose first addend is a row, and $d/dt$, $\delta$, $\int$ of a row.

- `s * v.T()` builds `SVMul(Transpose(v), s)`, printed `(v)'*s`; `v.T() * s` builds
  `Transpose(SVMul(v, s))`, printed `(v*s)'`. Both are the same row; `simplify` rewrites
  $(s v^T) w$ to $s (v^T w)$, but $(v^T s) w$ and $s (v^T w)$ are not merged yet.
- `v.T() * w` is the inner product (`VVMul`, a scalar), `v * w.T()` the outer product (a
  matrix). Integration by parts works on `Dot`; `collect` turns an inner product holding
  the velocity of a variation vector into a `Dot`, so `v.T() * A * v / 2` in a kinetic
  energy derives (tested for $\tfrac12 \dot x^T A \dot x$ and $\tfrac12 \Omega^T J \Omega$,
  which gives the same string as $\tfrac12 \Omega \cdot J\Omega$). `Dot(v, A * v)` is still
  the form the rest of the engine knows best: prefer it when there is a choice.
- `Dot`, `Cross` and `Hat` take columns; `Dot(a.T(), b)`, `Cross(a, b.T())`, `Hat(a.T())`
  and a sum mixing a row and a column (`a.T() + b`) are still built (existing code
  constructs them) but `check_tree` rejects them.

## Zeros, numbers and folding

- **Zero leaves.** `Zero` prints `0`, `ZeroVector` prints `0v`, `ZeroMatrix` prints `0`;
  `Number(0)` prints `(0)` and is flagged zero. `Transpose(ZeroVector)` (`(0v)'`) is the
  zero row that `linearize` uses for an empty scalar-row × vector-column block. The
  variation and the derivative of a `Constant` symbol are the zero leaves of its type
  (`Zero`, `ZeroVector`, `ZeroMatrix`; a vector or matrix of a size other than 3 or
  $3 \times 3$ gets a zero leaf of its own size). `geomech.core.base.constants` has
  `zero_like(expr)` and `negate(expr)` for code that builds trees.
- **`0 + x` is `x`** for every type, so `sum([a, b, c])` works for vectors and matrices.
  **`x + 0` keeps the zero** until `simplify` for a scalar (`m + 0` prints `(m+(0))`), and
  is an error for a vector or a matrix (`0` is a scalar: `VAdd cannot combine VECTOR 'x'
  and SCALAR '(0)'`); write `x + ZeroVector` or leave the zero out.
- **Products.** Every product node prints its two factors joined by `*`: `Mul` (`m*l`),
  `SVMul` (`x*m`), `SMMul` (`J*(2)`), `MVMul` (`J*\Omega`), `MMMul` (`R*J`), `VVMul`
  (`(a)'*b`, `a*(b)'`). Before phase 5b the factors were juxtaposed, so
  `Mul(ab, c)` and `Mul(a, bc)` both printed `abc` and compared equal. `repr` keeps `·` between factors and `*` as the constant marker.
- **Numbers.** `Number(2.5)` prints `(2.5)`, `Number(-1)` `(-1)`; `x * -1` and `-x` are
  `SVMul(x, (-1))`, printed `x*(-1)`. A scalar factor is always normalised to the right of
  `SVMul`/`SMMul` (`m * x` is `SVMul(x, m)`, printed `x*m`). `simplify` multiplies the
  numeric factors of a product chain into one number placed last
  (`Mul(-0.5, Mul(m, 2)) → Mul(m, -1)`, `(-1)*m*l*l → m*l*l*(-1)`); a float with an integer
  value becomes an `int` ($0.5 \cdot 2 \cdot (-1)$ prints `(-1)`, not `(-1.0)`, and a bare
  `Number(2.0)` simplifies to `(2)`), and every number literal is written the way
  `Number` writes it (`Scalar("0.5", value=0.5)` simplifies to `(0.5)`), so every equation
  prints its coefficients the same way. Numeric addends are summed.
- **Like terms.** Two addends merge when, after flattening their products, they have the
  same vector or matrix base and the same symbolic factors as a multiset; their numbers
  are added: `x*(2)*y + x*y*(3) → x*y*(5)`, $v a - v a \to 0$,
  $(v^T w) s + s (v^T w) \to 2\,(v^T w)\, s$ ($(v^T s) w$ is first pulled to
  $(v^T w) s$). Addends with different symbolic coefficients are not factored ($x a + x b$
  stays). Scalar products, `Dot` and the inner product of two columns are treated as
  commutative ($a^T b\, s - b^T a\, s$ is $0$); bases and factors built with MMMul, MVMul,
  Cross, SVMul, SMMul or an outer product are compared operand by operand, so
  $A B v s - B A v s$ and $(a \times b) m + (b \times a) m$ are not merged. Grouping is one
  pass, keyed by `geomech.core.transformations.simplification.term_key`.
- **Float cancellation.** Like terms are merged by their numeric coefficients; a sum
  whose magnitude is at most `1e-12` times the largest coefficient that went into it is
  treated as an exact zero (`_CANCEL_RTOL` in `simplification.py`): `x*0.1*3 - x*0.3` is `0`, not
  `x(5.55e-17)`.
- **A named number.** `Scalar("(-1)")` without `value` prints exactly like `Number(-1)`
  but is a symbol; `check_tree` rejects it. `Scalar("m", value=2)` is a *symbol* with a
  value: every simplification keeps it by name (`m + m` is `m*(2)`, never `(4)`; with
  `l = Scalar("l", value=3)`, $m \cdot l$ stays `m*l`, never `(6)`), `to_latex` prints `m`, and only
  `geomech.numeric` (`Evaluator`, `evaluate`, `check_eom`) uses the 2. A number literal is a scalar leaf printed as its own value:
  `Number(2)`, `Zero`, `One`, `Scalar("0.5", value=0.5)` and the numbers the operators
  wrap (`2 * x`); `geomech.core.base.expressions.is_number_literal` tells them apart.

## `str` is the identity

`a == b` is `type(a) is type(b) and str(a) == str(b)`, plus, for a calculus node
(`Variation`, `TimeDerivative`, `TimeIntegral`), the classes of the nodes down to its first
non-calculus operand; hashing uses `str`. `expr.has(x)` uses the same `==`. Hence:

- two leaves with the same printed name are the same symbol, whatever their flags
  (`Scalar("k", attr=["Constant"]) == Scalar("k")`), so `J2 = Matrix("J1")` silently
  merges two inertias, and two `S2("q")` objects are the same point;
- leaves of different classes differ even with one name: `Scalar("x") != Vector("x")`,
  `Variation(Scalar("x")) != Variation(Vector("x"))`,
  `Transpose(Vector("x")) != Transpose(Matrix("x"))`, and `Scalar("m")` does not occur in
  `Dot(Vector("m"), w)`;
- products print a `*` between their factors, and `*` is reserved for that: a symbol
  name may not contain it (`Scalar("m*l")` raises `ValueError`; write `ml`, or `x^\ast`
  for a superscript star). Every `*` in a printed form is therefore a product boundary,
  so two products print alike only when they multiply factors that print alike, in the
  same order: `Mul(ab, c)` is `ab*c` and `Mul(a, bc)` is `a*bc`, but the regroupings
  `(a*b)*c` and `a*(b*c)` both print `a*b*c` and are equal (they have one value);
- commutative rearrangements are different trees (`a*b` vs `b*a`), so like-term merging
  uses a structural key (`term_key`) that knows `Dot`, `Mul` and inner products commute,
  and a few equivalent terms are not merged;
- nodes are frozen dataclasses with tuple operands: a tree cannot change after it is
  built (assignment raises `FrozenInstanceError`), so a node can be shared between trees
  and used as a dictionary key; pickling and `copy.deepcopy` work;
- the pinned tests compare `str` of every equation byte for byte. `repr` is the readable
  debugging form (leaf kinds `v:x`, `M:J`, `S2:q`, `SO3:R`, a trailing `*` for a constant,
  `·` for products) and is not an identity.

## Glossary

- **configuration variable**: A leaf listed in `SystemVariables`: a `Scalar`, a `Vector`, an `S1` point, an `S2`
  point or an `SO3` point. Hamilton's principle varies it; each one gets an equation.
- **tangent vector**: The velocity coordinate of a manifold point: the scalar $\omega_q$ for `S1("q", axis=e)`
  ($\dot q = \omega\, e \times q$, a `TS1` leaf), $\omega_q$ for `S2("q")` ($\dot q = \omega \times q$,
  $\omega \perp q$) and $\Omega_R$ for `SO3("R")` ($\dot R = R \widehat{\Omega}$, body frame). Obtained with
  `q.get_tangent_vector()`; a `TS1`/`TS2`/`TSO3` leaf that remembers its point.
- **variation vector**: The free variation of a configuration variable, whose coefficient is the equation:
  $\delta x$ for a `Scalar`/`Vector`, the scalar $\gamma_q$ for `S1` ($\delta q = \gamma\, e \times q$), $\xi_q$ for
  `S2` ($\delta q = \xi \times q$, $\xi \perp q$) and $\eta_R$ for
  `SO3` ($\delta R = R \widehat{\eta}$). The keys of the `compute_eom` result are the printed variation
  vectors. Used on the EOM side.
- **perturbation**: A variable of the linearization (VBL side): the configuration perturbations ($\delta x$,
  $\xi$, $\eta$), the velocity perturbations ($\delta\dot x$, $\delta\omega$, $\delta\Omega$), their derivatives and the
  input variations $\delta u$. They are the columns of $M$, $A$, $B$.
- **params**: Prescribed time-varying symbols declared with `compute_eom(..., params=[l])`: their
  rate ($\dot l$) appears in the equations, their variation is not free and is dropped.
  `to_standard_form` accepts them; `linearize` does not yet.
- **inputs**: The forces and torques passed to `to_standard_form`/`linearize` (`inputs=[f, M]`):
  $G$ (standard form) and $B$ (linearization) are their coefficients. For `compute_eom`
  they are ordinary non-constant symbols in $\delta W$. An input must be the same symbol as its
  leaf in the equations (class, `attr`, size): `Scalar("g")` against a constant
  `Scalar("g", attr=["Constant"])` raises `ValueError`, as does a param of `compute_eom`.
  Both functions, and `StateLayout`, accept and reject the same inputs (`check_inputs`):
  an ordered sequence (a list, tuple or generator) of distinct non-constant
  `Scalar`/`Vector` leaves that are neither a configuration variable, a tangent or
  variation vector ($\omega$, $\xi$, $\Omega$, $\eta$), nor printed like a velocity, an acceleration or
  any other node of the equations; `inputs=None` means none. A `set` or a `dict` raises
  `TypeError`: the order of `inputs` is the order of the columns of $B$ and of the
  generated $u$ vector, and a set has no order (it would depend on the hash seed). A size
  may be spelled `3`, `[3]` or `(3,)`: leaves store it as a tuple. A param of
  `compute_eom` may also be an input only while its rate does not occur in the equations
  (`AlgebraicError` otherwise).
- **constant**: A leaf with `attr=["Constant"]` (or a number): its variation and derivative vanish.
  A symbol that is neither a configuration variable, a param nor a constant makes
  `compute_eom` raise ([Troubleshooting](troubleshooting.md)).
