# Linearization

`linearize` turns the equations of motion produced by `compute_eom` into their
linearization along an arbitrary trajectory,

$$
\mathbf{M}\,\delta\dot{x} = \mathbf{A}\,\delta x + \mathbf{B}\,\delta u ,
$$

where the perturbation state $\delta x$ is built from manifold-respecting
perturbation variables rather than local coordinates. This is the variation-based
linearization of Wu & Sreenath (2015) and of Lee et al., Appendix B; it is the model
that LQR and variation-based geometric controllers consume. The result is a set of
symbolic block matrices that are functions of the trajectory $(q, v, \dot v, u)$.

```python
from geomech import compute_eom, linearize

eom = compute_eom(L, dW, variables)
lin = linearize(eom, variables, inputs=[M])
print(lin)
```

## Perturbation variables

Each configuration variable contributes a configuration perturbation and a velocity
perturbation. The velocity perturbation is the `Variation` of the velocity symbol; unlike
in `compute_eom`, it is **not** rewritten by the manifold rules, because it is a state of
the linearized system.

| Configuration | Configuration perturbation | Velocity perturbation |
|---|---|---|
| `s: Scalar` ($\mathbb{R}$) | $\delta s$ = `Variation(s)` | $\delta\dot s$ = `TimeDerivative(Variation(s))` |
| `x: Vector` ($\mathbb{R}^3$) | $\delta x$ = `Variation(x)` | $\delta\dot x$ = `TimeDerivative(Variation(x))` |
| `q: S1` | $\gamma$ = `q.get_variation_vector()` (scalar) | $\delta\omega$ = `Variation(ω)`, $\omega$ = `q.get_tangent_vector()` (scalar) |
| `q: S2` | $\xi$ = `q.get_variation_vector()`, $\xi \perp q$ | $\delta\omega$ = `Variation(ω)`, $\omega$ = `q.get_tangent_vector()` |
| `R: SO3` | $\eta$ = `R.get_variation_vector()` | $\delta\Omega$ = `Variation(Ω)`, $\Omega$ = `R.get_tangent_vector()` |

The engine produces $\delta q = \gamma\, e \times q$ (S¹), $\delta q = \xi \times q$ (S²) and
$\delta R = R\,\hat\eta$ eagerly,
$\delta(\text{constant}) = 0$, and canonicalises $\delta(\dot v)$ as
`TimeDerivative(Variation(v))`. Each input $u$ (a `Scalar` or `Vector`) contributes
$\delta u$ = `Variation(u)`.

`LinearizedSystem.states` lists the configuration perturbations in `SystemVariables`
order (scalars, vectors, matrices) followed by the velocity perturbations in the same
order; `dstates` is `[s.t_diff() for s in states]` and `inputs` is
`[u.delta() for u in inputs]`.

## Rows

The system has one **kinematic row** per configuration variable, labelled
`"kin:" + str(var)`, followed by one **dynamic row** per equation of motion, labelled by
the `compute_eom` key (the variation vector).

The kinematic rows come from commuting $\delta$ and $d/dt$ on the manifold
kinematics:

- $\mathbb{R}$, $\mathbb{R}^3$: $\tfrac{d}{dt}\delta x = \delta\dot x$, so
  $\mathbf{M}[\tfrac{d}{dt}\delta x] = I$, $\mathbf{A}[\delta\dot x] = I$ (the two symbols print
  identically, `TimeDerivative(Variation(x))`; the row is an identity included so that every
  configuration variable has a kinematic row).
- $S^1$ ($\dot q = \omega\, e\times q$, $\delta q = \gamma\, e \times q$, scalar $\omega$,
  $\gamma$): $\delta\omega = \dot\gamma$ (Lee et al. eq. 4.11), so $\dot\gamma = \delta\omega$,
  $\mathbf{M}[\dot\gamma] = 1$, $\mathbf{A}[\delta\omega] = 1$: one scalar tangent
  coordinate, like $\mathbb{R}$.
- $S^2$ ($\dot q = \omega\times q$, $\omega\perp q$):
  $\dot\xi = q q^\top \hat\omega\, \xi + (I - q q^\top)\,\delta\omega$, so
  $\mathbf{M}[\dot\xi] = I$, $\mathbf{A}[\xi] = q q^\top \hat\omega$,
  $\mathbf{A}[\delta\omega] = I - q q^\top$.
- $SO(3)$ ($\dot R = R\hat\Omega$): $\dot\eta = -\hat\Omega\,\eta + \delta\Omega$, so
  $\mathbf{M}[\dot\eta] = I$, $\mathbf{A}[\eta] = -\hat\Omega$,
  $\mathbf{A}[\delta\Omega] = I$.

## Sign convention

`compute_eom` returns, for each variation vector, an expression $E(q, v, \dot v, u)$
with $E = 0$. By the library convention $E$ is the coefficient of the variation in
the integrand after integration by parts, i.e. the *negative* of the Euler–Lagrange
expression (which is why `to_standard_form` reports $M = -J$ for a rigid body). VBL takes
the first variation $\delta E$, which is linear in the perturbation variables, and
defines the blocks so that the physical inertia is positive:

$$
\delta E = -\mathbf{M}\,\delta\dot v + \mathbf{A}\,\delta x + \mathbf{B}\,\delta u,
\qquad
\mathbf{M}[\delta\dot v_j] = -\operatorname{coeff}(\delta\dot v_j, \delta E),\quad
\mathbf{A}[\delta x_j] = \operatorname{coeff}(\delta x_j, \delta E),\quad
\mathbf{B}[\delta u_j] = \operatorname{coeff}(\delta u_j, \delta E).
$$

Every additive term of $\delta E$ must contain exactly one perturbation symbol. If a
term contains none (an input that was not declared in `inputs`, a symbol missing the
`Constant` attribute) or contains a perturbation only inside a calculus node that the
extraction cannot see through, `linearize` raises `AlgebraicError` naming the term.
Nothing is silently set to zero.

## Worked example: free rigid body

```python
from geomech import SO3, Dot, Matrix, Scalar, SystemVariables, Vector, compute_eom, linearize

J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
M = Vector("M")
R = SO3("R")
Om, eta = R.get_tangent_vector(), R.get_variation_vector()
half = Scalar("0.5", value=0.5, attr=["Constant"])

eom = compute_eom(Dot(Om, J * Om) * half, Dot(eta, M), SystemVariables(matrices=[R]))
lin = linearize(eom, SystemVariables(matrices=[R]), inputs=[M])
print(lin)
```

```
LinearizedSystem: M d/dt(dx) = A dx + B du
  states : [\eta_{R}, \delta{\Omega_{R}}]
  inputs : [\delta{M}]
  row kin:R
    M[\frac{d}{dt}(\eta_{R})] = I
    A[\eta_{R}] = Hat(\Omega_{R})*(-1)
    A[\delta{\Omega_{R}}] = I
  row \eta_{R}
    M[\frac{d}{dt}(\delta{\Omega_{R}})] = J
    A[\delta{\Omega_{R}}] = (Hat(J*\Omega_{R})+Hat(\Omega_{R})*J*(-1))
    B[\delta{M}] = I
```

The equation of motion is $E = -J\dot\Omega - \Omega\times J\Omega + M$
(Lee et al. Eq. 6.16). Its variation is

$$
\delta E = -J\,\delta\dot\Omega - \delta\Omega\times J\Omega - \Omega\times J\,\delta\Omega + \delta M
        = -J\,\delta\dot\Omega + \big(\widehat{J\Omega} - \hat\Omega J\big)\,\delta\Omega + \delta M,
$$

so the blocks read, with $\delta x = (\eta, \delta\Omega)$,

$$
\begin{bmatrix} I & 0 \\ 0 & J \end{bmatrix}
\begin{bmatrix} \dot\eta \\ \delta\dot\Omega \end{bmatrix}
=
\begin{bmatrix} -\hat\Omega & I \\ 0 & \widehat{J\Omega} - \hat\Omega J \end{bmatrix}
\begin{bmatrix} \eta \\ \delta\Omega \end{bmatrix}
+
\begin{bmatrix} 0 \\ I \end{bmatrix} \delta M ,
$$

which is the textbook linearization of the rigid body on $TSO(3)$.

## Reading the blocks

`M`, `A` and `B` are nested dicts indexed by row label and by the printed column symbol:
`lin.M[row][str(dstate)]`, `lin.A[row][str(state)]`, `lin.B[row][str(input_variation)]`.
They are dense: every (row, column) pair is present, and zero blocks are the library's
zero constants. The shape of a block follows from the row and column types,
`lin.block_shape(row, col)`:

| row \ column | vector perturbation | scalar perturbation |
|---|---|---|
| vector equation | 3×3 `Matrix` expression | 3×1 `Vector` expression |
| scalar equation | 1×3 `Transpose(Vector)` | `Scalar` expression |

For example, in the quadrotor with thrust `f: Scalar` and moment `M: Vector`, the
translational row has `B[δf] = R*e3` (3×1) and `A[η] = R*Hat(e3)*f*(-1)` (3×3), and the
elastic spherical pendulum's scalar row has the 1×3 block `A[δω] = (ω*m*l*(2))'`.

`str(lin)` lists the non-zero blocks row by row; `to_latex(block)` renders any block.

### Blocks on $S^2$ are defined on the constraint set

For an $S^2$ variable the perturbation $(\xi, \delta\omega)$ is constrained:
$\xi \perp q$ and, from $\delta(q\cdot\omega) = 0$,
$q\cdot\delta\omega = \xi\cdot(\omega\times q)$. `compute_eom` simplifies with
$q\cdot\omega = 0$ before `linearize` takes the variation (for example
$\omega\times(\omega\times q) \to -|\omega|^2 q$), and the variations of the two forms
differ by terms proportional to $q\cdot\delta\omega - \xi\cdot(\omega\times q)$, which vanish
on the constraint set. The $\mathbf{A}$ blocks against $\xi$ and $\delta\omega$ are therefore
unique as linear maps on the constraint set, and agree there with any hand derivation, but
two correct derivations can print different $3\times 3$ matrices. For the quadrotor with a
cable-suspended load, `A[δx_L][ξ]` prints as
$-m_Q l\,\hat{\dot\omega}\hat q + m_Q l\,|\omega|^2 \hat q$ while the derivation from
$\ddot q = \dot\omega\times q + \omega\times(\omega\times q)$ gives
$-m_Q l\,(\hat{\dot\omega} + \hat\omega^2)\hat q$; both act identically on tangent
$(\xi, \delta\omega)$. A consumer that needs a minimal (two-dimensional) representation
must project onto the tangent space (Lee et al., Appendix B).

## What `linearize` does not do

- It does not reduce the perturbation to minimal coordinates at an equilibrium (Lee et al.
  Appendix B basis $\sigma$); the output stays in the full perturbation vectors, as in
  Wu & Sreenath.
- It does not synthesise controllers, discretise, simulate or generate code itself:
  [Code generation](codegen.md) turns a `LinearizedSystem` into `M_lin`, `A_lin`, `B_lin`
  functions, and [Simulation](simulation.md) designs an LQR gain from them
  (`geomech.simulation.lqr`) and simulates the closed loop.
- It does not call the EOM manifold rules on the velocity variations, and it does not
  change `compute_eom` or `to_standard_form`.

## Limitations

- The transpose of a sum or of an outer product (`Transpose(MAdd)`, `Transpose(VVMul)`,
  `Transpose(Hat(a + b))`) and `Vee` are opaque to `extract_linear_coeff`. `expand`
  distributes `VVMul` over `VAdd` factors (including a `VAdd` inside the transposed factor)
  but does not expand `Transpose` itself, so an equation whose variation carries a
  perturbation inside one of these nodes raises `AlgebraicError` instead of linearizing.
  Rewrite the equation without the transpose (for example `Transpose(Hat(a)) = -Hat(a)`,
  `Transpose(a bᵀ) = b aᵀ`) before calling `linearize`.
- A stray perturbation symbol in an equation (`η`, `Variation(Ω)`, `ξ̇`) is not an input
  `linearize` can interpret; it surfaces as the engine's `ValueError` (the variation of a
  variation is not representable) or as the `AlgebraicError` described above.
- An equation whose variation is identically zero produces an all-zero dynamic row; it is
  not reported as an error.
- Prescribed time-varying parameters (`compute_eom(..., params=[l])`) are not supported:
  the variation of the equation contains `δl`, which is neither a state nor an input
  perturbation, and `linearize` raises `AlgebraicError` naming the term. Listing `l` in
  `inputs` does not help, because `l̇` brings in `d/dt(δl)`.
