# Conventions

The choices below are fixed by the code; each one has surprised someone reading an
equation. Equation numbers refer to Lee et al., *Global Formulations of Lagrangian and
Hamiltonian Dynamics on Manifolds*.

## Signs

**`E = −(Euler–Lagrange) + forcing`.** `compute_eom` returns, for each variation vector,
the coefficient `E` of that variation in the integrand of δS after integration by parts,
and the equation is `E = 0`. For `L(x, ẋ)` and `δW = δx·f` this is

    E = −(d/dt ∂L/∂ẋ − ∂L/∂x) + f,

so the point mass `L = ½ m ẋ·ẋ` with `δW = δx·F` gives `E = −m ẍ + F`, and the free rigid
body gives `E = −J Ω̇ − Ω × JΩ + M` (Eq. 6.16 with every term moved to one side).

**`M = −(physical inertia)` in `to_standard_form`.** The blocks split the same `E`:
`E = M a + f + G u`, so `M[ẍ] = −m I` for the point mass, `M[Ω̇] = −J` for a rigid body, and
`G[F] = I` for a force entering as `δx·F`.

**VBL.** `linearize` takes `δE` and defines its blocks so that the physical inertia is
positive: `δE = −M δv̇ + A δx + B δu`, i.e. `M[δv̇] = −coeff(δv̇, δE)`, and the result reads
`M δẋ = A δx + B δu` ([Variation-Based Linearization](vbl.md#sign-convention)).

## S² equations live on the tangent plane

For a point `q ∈ S²` the variation is `δq = ξ × q` with `ξ ⊥ q`, so only the part of `E`
orthogonal to `q` is an equation: `(I − q qᵀ) E = 0`. `E` as derived in general has a
component along `q` (from the forcing, from `ω × c`, from `A ω̇`), and `to_standard_form` and
`linearize` return the terms as derived, not projected. To integrate, solve
`(I − qqᵀ)(M ω̇ + f + Σ G u) = 0` together with `ω̇ ⊥ q`; the unprojected 3×3 solve gives
`q·ω̇ ≠ 0` and a trajectory that leaves the sphere. `numeric.check_eom` compares on the
tangent plane for the same reason.

## Rows and columns

There is one `VECTOR` type; a row is recognised by its structure (`is_row_vector`):
`Transpose(column)`, a row times a matrix `vᵀM` (`MVMul(row, M)`), a scalar multiple of a
row, a sum whose first addend is a row, and d/dt, δ, ∫ of a row.

- `s * v.T()` builds `SVMul(vᵀ, s)`, printed `(v)'s`; `v.T() * s` builds
  `Transpose(SVMul(v, s))`, printed `(vs)'`. Both are the same row; `simplify` rewrites
  `(s vᵀ) w` to `s (vᵀ w)`, but `(vᵀ s) w` and `s (vᵀ w)` are not merged yet.
- `v.T() * w` is the inner product (`VVMul`, a scalar), `v * w.T()` the outer product (a
  matrix). Integration by parts works on `Dot`; `collect` turns an inner product holding
  the velocity of a variation vector into a `Dot`, so `v.T() * A * v / 2` in a kinetic
  energy derives (tested for `½ ẋᵀAẋ` and `½ ΩᵀJΩ`, which gives the same string as
  `½ Ω·JΩ`). `Dot(v, A * v)` is still the form the rest of the engine knows best:
  prefer it when there is a choice.
- `Dot`, `Cross` and `Hat` take columns; `Dot(a.T(), b)`, `Cross(a, b.T())`, `Hat(a.T())`
  and a sum mixing a row and a column (`a.T() + b`) are still built (existing code
  constructs them) but `check_tree` rejects them.

## Zeros, numbers and folding

- **Zero leaves.** `Zero` prints `0`, `ZeroVector` prints `0v`, `ZeroMatrix` prints `0`;
  `Number(0)` prints `(0)` and is flagged zero. `Transpose(ZeroVector)` (`(0v)'`) is the
  zero row that `linearize` uses for an empty scalar-row × vector-column block. The
  variation and the derivative of a `Constant` symbol are zero leaves of its type, but
  not always the shared constants: for a constant vector both print `0` (a zero vector
  leaf, not `ZeroVector`'s `0v`), for a constant matrix `δ` prints `O` and `d/dt` prints
  `0`. All of them are flagged zero, and `simplify` removes them; only the printed form
  differs.
- **`0 + x` is `x`** for every type, so `sum([a, b, c])` works for vectors and matrices.
  **`x + 0` keeps the zero** until `simplify` for a scalar (`m + 0` prints `(m+(0))`), and
  is an error for a vector or a matrix (`0` is a scalar: `VAdd cannot combine VECTOR 'x'
  and SCALAR '(0)'`); write `x + ZeroVector` or leave the zero out.
- **Numbers.** `Number(2.5)` prints `(2.5)`, `Number(-1)` `(-1)`; `x * -1` and `-x` are
  `SVMul(x, (-1))`, printed `x(-1)`. A scalar factor is always normalised to the right of
  `SVMul`/`SMMul` (`m * x` is `SVMul(x, m)`, printed `xm`). `simplify` multiplies the
  numeric factors of a product chain into one number placed last
  (`Mul(-0.5, Mul(m, 2)) → Mul(m, -1)`, `(-1)mll → mll(-1)`); a float with an integer
  value becomes an `int` (`0.5 · 2 · (-1)` prints `(-1)`, not `(-1.0)`), so every
  equation prints integer coefficients the same way. Numeric addends are summed.
- **Like terms.** Two addends merge when, after flattening their products, they have the
  same vector or matrix base and the same symbolic factors as a multiset; their numbers
  are added: `x(2)y + xy(3) → xy(5)`, `v a − v a → 0`, `(vᵀw)s + s(vᵀw) → (vᵀw)s(2)`
  (`(vᵀ s) w` is first pulled to `(vᵀ w) s`). Addends with different symbolic
  coefficients are not factored (`x a + x b` stays). Only scalar products and Dot are treated as commutative: bases and
factors built with MMMul, MVMul, VVMul, Cross, SVMul or SMMul are compared operand by
operand, so A B v s − B A v s, (a×b) m + (b×a) m and aᵀ b s − bᵀ a s are not merged (the
last is equal but not recognised).
- **Float cancellation.** Like terms are merged by their numeric coefficients; a sum
  whose magnitude is at most `1e-12` times the largest coefficient that went into it is
  treated as an exact zero (`_CANCEL_RTOL` in `simplification.py`): `x*0.1*3 − x*0.3` is `0`, not
  `x(5.55e-17)`.
- **A named number.** `Scalar("(-1)")` without `value` prints exactly like `Number(-1)`
  but is a symbol; `check_tree` rejects it. `Scalar("m", value=2)` is a *symbol* with a
  value: every simplification keeps it by name (`m + m` is `m(2)`, never `(4)`; with
  `l = Scalar("l", value=3)`, `m·l` stays `ml`, never `(6)`), `to_latex` prints `m`, and only
  `geomech.numeric` (`Evaluator`, `evaluate`, `check_eom`) uses the 2 (decision D1 in
  `specs/fix-plan.md`). A number literal is a scalar leaf printed as its own value:
  `Number(2)`, `Zero`, `One`, `Scalar("0.5", value=0.5)` and the numbers the operators
  wrap (`2 * x`); `geomech.core.base.expressions.is_number_literal` tells them apart.

## `str` is the identity

`a == b` is `type(a) is type(b) and str(a) == str(b)`, and hashing uses `str`. Hence:

- two leaves with the same printed name are the same symbol, whatever their flags
  (`Scalar("k", attr=["Constant"]) == Scalar("k")`), so `J2 = Matrix("J1")` silently
  merges two inertias, and two `S2("q")` objects are the same point;
- commutative rearrangements are different trees (`a*b` vs `b*a`), so like-term merging
  uses a structural comparison that knows `Dot` and `Mul` commute, and a few equivalent
  terms are not merged;
- the pinned tests compare `str` of every equation byte for byte. `repr` is the readable
  debugging form (leaf kinds `v:x`, `M:J`, `S2:q`, `SO3:R`, a trailing `*` for a constant,
  `·` for products) and is not an identity.

## Glossary

- **configuration variable**: A leaf listed in `SystemVariables`: a `Scalar`, a `Vector`, an `S2` point or an `SO3`
  point. Hamilton's principle varies it; each one gets an equation.
- **tangent vector**: The velocity coordinate of a manifold point: `ω_q` for `S2("q")` (`q̇ = ω × q`,
  `ω ⊥ q`) and `Ω_R` for `SO3("R")` (`Ṙ = R Hat(Ω)`, body frame). Obtained with
  `q.get_tangent_vector()`; a `TS2`/`TSO3` leaf that remembers its point.
- **variation vector**: The free variation of a configuration variable, whose coefficient is the equation:
  `δx` for a `Scalar`/`Vector`, `ξ_q` for `S2` (`δq = ξ × q`, `ξ ⊥ q`) and `η_R` for
  `SO3` (`δR = R Hat(η)`). The keys of the `compute_eom` result are the printed variation
  vectors. Used on the EOM side.
- **perturbation**: A variable of the linearization (VBL side): the configuration perturbations (`δx`,
  `ξ`, `η`), the velocity perturbations (`δẋ`, `δω`, `δΩ`), their derivatives and the
  input variations `δu`. They are the columns of `M`, `A`, `B`.
- **params**: Prescribed time-varying symbols declared with `compute_eom(..., params=[l])`: their
  rate (`l̇`) appears in the equations, their variation is not free and is dropped.
  `to_standard_form` accepts them; `linearize` does not yet.
- **inputs**: The forces and torques passed to `to_standard_form`/`linearize` (`inputs=[f, M]`):
  `G` (standard form) and `B` (linearization) are their coefficients. For `compute_eom`
  they are ordinary non-constant symbols in `δW`.
- **constant**: A leaf with `attr=["Constant"]` (or a number): its variation and derivative vanish.
  A symbol that is neither a configuration variable, a param nor a constant makes
  `compute_eom` raise ([Troubleshooting](troubleshooting.md)).
