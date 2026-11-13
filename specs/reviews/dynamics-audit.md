# Adversarial audit: `geomech.dynamics` (branch `feature/vbl`, HEAD `54d97b9`)

Scope: `compute_eom` / `SystemVariables` / `to_standard_form` / `linearize` and their
interaction with `manifold_rules`, `ibp`, `collect`, `extract`. Read-only audit: nothing under
`geomech/` or `tests/` was modified; this file is the only addition. Baseline at HEAD:
`python -m pytest tests/ -q` → 1322 passed, 5 xfailed; `ruff` clean.

## Method

Every equation produced by `compute_eom` was checked numerically against an **independent
derivation that does not use the library's variational rules**. The harness
(`specs/reviews/scripts/euler_lagrange_check.py`, pasted in Appendix A) evaluates the *numeric* Lagrangian with the
test `Evaluator` and forms the Euler–Lagrange residual by central finite differences:

- R, R³: `g = d/dt ∂L/∂v − ∂L/∂q`;
- S²: constrained Euler–Lagrange in R³ with `L̃(q, q̇) = L(q, q × q̇)` and the multiplier
  removed by projection, `E_true = −q × g + f` (only tangent components compared);
- SO(3): `g = d/dt ∂L/∂Ω − (∂L/∂Ω) × Ω − vee(A − Aᵀ)`, `A = Rᵀ ∂L/∂R`, gradients taken
  on the nine entries of `R`.

`d/dt` is a directional derivative along the sampled `(q̇, q̈)` so no integration is needed:
the library's `E(q, v, v̇)` and `−g + forcing` must agree for *random* accelerations. The
harness reproduces all nine pinned systems in `tests/vbl_systems.py` to 1e-6 or better
(finite-difference noise), so a mismatch of 1e-1 is a wrong equation, not noise.
`to_standard_form` was checked by evaluating `Σ M[a]·a + f + Σ G[u]·u` against `E` at random
points; `linearize` with `tests/numeric_helpers.finite_difference_check`.

The probe scripts were run from a scratch directory; the harness is kept at `specs/reviews/scripts/euler_lagrange_check.py`.

## Findings, by severity

### F1 — SILENT WRONG EOM: the S² rule `δω = ξ̇ − ω × ξ` drops the component along q

**Severity: high (silent wrong output on legal input).**

`geomech/core/transformations/manifold_rules.py:26,66-70` substitutes
`Variation(ω) → ξ̇ − ω × ξ` (also `expressions.py:348-349`, `TS2.delta(substitute=True)`).
The correct variation for `q̇ = ω × q`, `δq = ξ × q`, `ξ ⊥ q` is (Lee et al. §5.3.3, p. 216)

    δω = −S(ω) ξ + (I − q qᵀ) ξ̇  =  ξ̇ − ω × ξ − q (q·ξ̇)  =  ξ̇ − 2 ω × ξ

(using `q·ξ̇ = −q̇·ξ = q·(ξ × ω)` and `ω × ξ ∥ q`). The library keeps `−ω × ξ` but
omits the equal-sized projection term `−q(q·ξ̇)`, so the substituted `δω` violates
`δ(q·ω) = 0` (`δ_lib(q·ω) = ξ·(q × ω) ≠ 0`). The error is invisible whenever `∂L/∂ω ∥ ω`
(every pinned system: isotropic `½ m l² ω·ω`, and every "q̇-form" Lagrangian, where ω enters
only through `ω × q`). It is wrong as soon as `∂L/∂ω` has a component along `q`:

| Lagrangian (S² variable q, `Dot(ξ, f)` forcing) | rel. error of E vs truth | after fix |
|---|---|---|
| anisotropic `½ ω·Aω − m g l q·e₃`, A symmetric constant | **2.3e-1** | 1.4e-7 |
| pendulum in a rotating frame, `½ m l² (|ω+W₀|² − (W₀·q)²) − m g l q·e₃` | **1.6e-1** | 2.2e-7 |
| gyroscopic term `½ m l² ω·ω + c·ω − m g l q·e₃` | **2.9e-1** | 1.5e-7 |
| R³×S² coupling `½ m ẋ·ẋ + ½ m l² ω·ω + ẋ·ω` | **1.8e-1** | 6.4e-7 |
| S²×SO(3) coupling `½ m ω·ω + ½ Ω·JΩ + k ω·Ω` | **1.2e0** | 6.9e-7 |
| S²×SO(3) coupling `… + k ω·(Rᵀe₃)` | **4.5e-2** | 5.6e-6 |

Hand check of the third row: `∂L/∂ω = m l² ω + c`; truth `E = −m l² ω̇ + 2 ω × c − m g l q × e₃`
(tangent part); library prints `… + Cross(ω, c) + …` — the gyroscopic term is **half** its true
value. Same physics written in q̇-form (`½ m l² |(ω+W₀)×q|²`) is correct (2.0e-7), i.e. the
derived equations depend on how the user spells the Lagrangian.

Reproducer (public API):

```python
from geomech import *
HALF = Scalar("0.5", value=0.5, attr=["Constant"])
m, g, l = getScalars("m g l", attr=["Constant"]); e3 = Vector("e3", attr=["Constant"])
c = Vector("c", attr=["Constant"]); q, f = S2("q"), Vector("f")
om, xi = q.get_tangent_vector(), q.get_variation_vector()
L = HALF*m*l*l*Dot(om, om) + Dot(c, om) - m*g*l*Dot(q, e3)
eom = compute_eom(L, Dot(xi, f), SystemVariables(vectors=[q]))
print(eom[str(xi)][1])
# (d/dt(ω)(-1)mll + Cross(ω,c) + Cross(q,e3)mgl(-1) + f)      <- expected Cross(ω,c)(2)
```

Independent check: `L̃(q, q̇) = ½ m l² |q̇|² + q̇·(c × q) − m g l q·e₃`,
`g = d/dt(m l² q̇ + c × q) − q̇ × c = m l² q̈ + 2 c × q̇`, `−q × g` projected gives
`2(c·q)(ω × q) = P(2 ω × c)`; the library's `P(ω × c) = (c·q)(ω × q)`.

**Fix.** `manifold_rules.py:26`: `_NEG1 → Scalar("(-2)", value=-2, …)` (and the comment on
line 66, `expressions.py:348-349`, `docs/manifold-implementation-notes.md:7,67,74` — the latter
already says "needs to be verified"). Validated by monkeypatching `_NEG1 = -2` at runtime
(runtime monkeypatch of `_NEG1`): all six wrong rows above drop to ≤ 6e-6, all pinned S²
systems (spherical, elastic, quadrotor-with-load, double spherical pendulum) are unchanged, and
the full suite with the patch applied via a pytest plugin gives **1321 passed, 1 failed**: only
`tests/test_manifold_rules.py::TestS2ManifoldRules::test_variation_of_omega_expands`, which pins
the `(-1)` coefficient, i.e. pins the bug. Note that with the fix the printed q̇-form equations
may change by terms parallel to `q` (e.g. `k (ω×q)·Ω` gains a `2` on `ω × (q × Ω)`); those are
not errors (S² equations are unique only on the tangent space, cf. `docs/vbl.md` "constraint
set" note), and the numeric check is unchanged (4.6e-7 before and after).

### F2 — SILENT WRONG EOM: `collect` rewrites `Dot(M v, b × c)` as `Dot(v, (M b) × c)`

**Severity: high (silent wrong output), narrow trigger.**

`geomech/core/math/collect.py:64-77`, case `(MVMul(), Cross())`, branch `l.right == vec`
(line 72-74) uses the identity `(M v)·(b × c) = v·((M b) × c)`, which is false for every
`M` except `I` (`Mᵀ(b × c) ≠ (M b) × c`, including rotations: `Rᵀ(b×c) = (Rᵀb)×(Rᵀc)`).

```python
from geomech import *
from geomech.core.math.collect import collect
J = Matrix("J", attr=["Constant", "SymmetricMatrix"]); v, x, y = Vector("v"), Vector("x"), Vector("y")
print(collect(Dot(MVMul(J, v), Cross(x, y)), v))   # Dot(v,Cross(Jx,y))  -- numeric lhs-rhs = -0.207
```

Pipeline reproducer (SO(3), angular momentum dotted with a configuration-dependent axis):

```python
HALF = Scalar("0.5", value=0.5, attr=["Constant"])
J = Matrix("J", attr=["Constant", "SymmetricMatrix"]); b, e3 = Vector("b", attr=["Constant"]), Vector("e3", attr=["Constant"])
R = SO3("R"); Om, eta = R.get_tangent_vector(), R.get_variation_vector(); M = Vector("M")
L1 = HALF*Dot(Om, J*Om) + Dot(J*Om, Cross(MVMul(Transpose(R), e3), b))      # J on the velocity
L2 = HALF*Dot(Om, J*Om) + Dot(Om, J*Cross(MVMul(Transpose(R), e3), b))      # identical physics
eom1 = compute_eom(L1, Dot(eta, M), SystemVariables(matrices=[R]))  # prints ...+Cross(JCross(Ω,(R)'e3),b)+...
eom2 = compute_eom(L2, Dot(eta, M), SystemVariables(matrices=[R]))  # prints ...+JCross(Cross(Ω,(R)'e3),b)+...
```

Harness: `eom1` rel. error **7.0e-1**, `eom2` 3.2e-7. Same with `k (RΩ)·(a × b)` (R the
configuration): **9.8e-1** vs 1.7e-7 for `k Ω·(Rᵀ(a × b))`. The trigger needs the velocity
*directly* under a matrix inside a `Dot` whose other side is a `Cross` with configuration
dependence; `Dot(JΩ, a×b)` with constants hides the bug because the wrong IBP term is
`d/dt(const) = 0`.

**Fix.** Replace line 74 by `return Dot(vec, MVMul(Transpose(l.left), r))` (the `(MVMul(), _)`
case at line 100-102 already does exactly this). Validated by monkeypatch
(runtime monkeypatch of `collect`): 3.4e-7, pinned SO(3) systems unchanged.

### F3 — CRASH on legal input: velocity nested under two matrices, or under a matrix inside a `Dot` with a `Cross`

**Severity: medium (crash, `AlgebraicError` from extraction, not silent).**

`collect.py:99-117` (`(MVMul(), _)`) and `:64-77` (`(MVMul(), Cross())`) give up
(`return expr`) when the target is inside the `MVMul` but not directly its right operand or
inside a `Cross`. The term survives IBP with `ξ̇`/`η̇` inside and `extract_linear_coeff`
raises *"'η_R' appears only inside the calculus operator d/dt(η_R)"*. All of these are
ordinary Lagrangians:

| Lagrangian | current | with fix |
|---|---|---|
| `½ Ω·JΩ + (R(JΩ))·e₃` or `((R J)Ω)·e₃` (inertial angular momentum) | AlgebraicError | 2e-7, 6e-7 |
| S² anisotropic q̇-form `½ q̇·A q̇ − m g q·e₃`, `q̇ = ω × q` | AlgebraicError | 4.8e-7 |
| S² `½ m q̇·q̇ + k (A q̇)·(q × e₃)` | AlgebraicError | 1.2e-6 |
| SO(3)×SO(3) `½ Ω₁·J₁Ω₁ + ½ Ω₂·J₂Ω₂ + k Ω₁·(R₁ᵀ(R₂ Ω₂))` (relative-rate coupling) | AlgebraicError | 7.2e-7 |

The second row is the same system as F1's "anisotropic A" written the "safe" way, so a user
who avoids F1 by using q̇-form runs into F3.

**Fix (validated together with F2 in `p12_fix_collect.py`).** In both cases, when `l.has(vec)`
and `l.right` is not `vec` (and, in `(MVMul(), _)`, not a `Cross`), peel the matrix instead of
returning: `return collect(Dot(l.right, MVMul(Transpose(l.left), r)), vec)`. Each peel removes
one matrix layer from the side containing `vec`, so it terminates (the loop that
`specs/reviews/scalar-config-variables-final.md` Finding 1 guarded against came from *flipping*
sides, not from peeling). With the patch the six pinned SO(3)/R³ systems are byte-identical.

### F4 — MISSING VALIDATION: duplicate configuration variables are accepted and silently merged

**Severity: medium (ill-posed system accepted; `linearize` emits duplicate rows/states).**

`eom.py:21-44 _check_variables` checks type and constancy only. Observed:

```python
q = S2("q"); ...
compute_eom(L, dW, SystemVariables(vectors=[q, q]))          # one key, no error
compute_eom(L + ½ ω'·ω', dW, SystemVariables(vectors=[q, S2("q")]))   # two objects, same name:
#   one equation with both d/dt(ω) terms merged: (d/dt(ω)(-1)mll + ... + d/dt(ω)(-1) + f)
linearize(eom, SystemVariables(vectors=[q, q]), [f])
#   rows ['kin:q', 'kin:q', 'ξ_q'], states ['ξ_q', 'ξ_q', 'δω_q', 'δω_q']  (dict rows overwritten)
```

`vbl._check_eom` (`vbl.py:355-372`) compares *sets* of keys, so duplicates pass there too.
A `Scalar("x")` and a `Vector("x")` together (both print `\delta{x}`) crash with an
`AlgebraicError` from `_scalar_from_vector`. **Fix:** in `_check_variables`, reject repeated
`str(var)` and repeated `str(var.get_variation_vector())` across all three groups
(`ValueError`), and have `_check_eom` compare the key *list* length as well.

### F5 — MISSING VALIDATION: `Variation` of a non-configuration symbol is silently dropped

**Severity: low–medium (silent, but the result is right for prescribed `p(t)`).**

A `Variation(p)` that survives in δS (from `dW = Dot(x.delta(), F) + Dot(p.delta(), F)`, or from
a parameter the user forgot to flag `Constant`, e.g. `e3 = Vector("e3")`) is dropped by
`extract_linear_coeff` because it contains no variation vector; `compute_eom` returns the same
equation as if the term were absent. For `L = ½ m ẋ·ẋ − m g e₃·x` with unflagged `e3` the
output is identical to the flagged case, so a user cannot tell whether `e3` was treated as a
configuration or as prescribed. Suggest: after `separate_variations`, walk δS for `Variation`
nodes whose inner is not a configuration variable / tangent vector and raise (or warn)
naming the symbol. (`linearize` already does the analogous check for its own perturbations.)

### F6 — CRASH on `inf_work = 0` / `None`

**Severity: low.** `eom.py:64 dS = dL + inf_work` → `AttributeError: 'int' object has no
attribute 'type'` for `compute_eom(L, 0, vars)` and for `None`; `Zero` works. Suggest
`inf_work = Zero if inf_work in (None, 0) else inf_work` and a type check (`ScalarExpr`) with a
clear `TypeError` (a vector-typed `dW` currently surfaces as `ExpressionMismatchError: Add`).

### F7 — CRASH on legal input: `Transpose(R) * v` raises `ExpressionMismatchError`

**Severity: low (usability; forces `MVMul(Transpose(R), v)`).**
`geometry.py:153-163 Transpose.__mul__` dispatches on `other.type` only; for a matrix-typed
`Transpose` times a vector it builds `VVMul` and the constructor raises
*"VVMul cannot be performed between MATRIX and VECTOR"*. `Rᵀ e₃` is the single most common
body-frame expression. Fix: branch on `self.type` first (`MATRIX × VECTOR → MVMul`,
`MATRIX × MATRIX → MMMul`). Related nit: `Scalar` has no `__neg__`, so `-m*g*Dot(e3, x)`
raises `TypeError`.

### F8 — MISLEADING OUTPUT: `to_latex` drops parentheses that change the meaning

**Severity: low–medium (a human copying the LaTeX gets a different equation).**
`printing.py:508-511` (`Cross`), `:492-497` (`MVMul`/`MMMul`), `:520-521` (`Transpose`),
`:470-480` (`SVMul` with a `Dot` scalar) never parenthesise composite operands:

| expression | `to_latex` | reads as |
|---|---|---|
| `Cross(a, Cross(b, c))` **and** `Cross(Cross(a, b), c)` | `a \times b \times c` (identical) | ambiguous; e.g. `Ω × (Ω × ρ)` reads as `(Ω×Ω)×ρ = 0` |
| `MVMul(Transpose(R), Cross(a, b))` **and** `Cross(MVMul(R,a), b)` | `R^T a \times b` / `R a \times b` | `(Rᵀa) × b ≠ Rᵀ(a × b)` |
| `Transpose(Cross(a, b))`, `Transpose(MVMul(R, a))` | `a \times b^T`, `R a^T` | transposes the last factor only |
| `SVMul(a, Dot(b, c))` | `b \cdot c a` | `b·(c a)` |
| `Cross(a, SVMul(b, -1))` | `a \times -b` | — |
| `Mul(Mul(s, s), -0.5)`, `Mul(s, 2.5)` | `-s s \frac{1}{2}`, `s 2.5` | coefficient after the symbols |

Real output, double rigid pendulum `η_{R1}` row: `-m1 \Omega_{R1} \times \rho_1 \times \Omega_{R1}
\times \rho_1 … -m2 l_1 \times R1^T R2 \Omega_{R2} \times \rho_2`. Suggest a precedence-aware
`to_latex` (wrap `Cross`, `VAdd`, `MVMul` operands of `Cross`/`Transpose`/`MVMul` in
`\left(\right)`, print `(b·c)\,a`). Separately, `print_eom` (`:19-28`), `eom_to_latex`/
`render_eom` (`:624-640`) emit `str(eqn)` inside LaTeX math, so the PDF shows
`Cross(ΩR, JΩR)(−1)+Cross(ρ,(R)′e3)mg(−1)` rather than typeset maths; they should call
`to_latex`. (`render_eom` itself works with the system `pdflatex`, `open_pdf=False`.)
`repr()` output (`cross(a, b)`, `<a, b>`, `R'*cross(a, b)`) is unambiguous and fine.

### F9 — `linearize` and non-constant parameters (documented limitation, not a bug)

A time-varying parameter (`m = Scalar("m")`, no `Constant`) makes `compute_eom` produce the
correct `ṁ` terms (verified, 8.9e-7) but `linearize` raises `AlgebraicError … contains no
perturbation symbol` on the `δm` term; listing `m` as an input instead raises on
`d/dt(δm)`. Both messages are clear. Worth one sentence in `docs/vbl.md` ("parameters must be
`Constant`; VBL does not linearise with respect to prescribed time-varying parameters").

### Nits (documentation / consistency)

- `variables.py:12` docstring: "scalars: … (not yet used in EOM pipeline)" is stale.
- `standard_form.py:25` `Σ M[a_j] * a_j`: for a scalar row with a vector acceleration `M[a]`
  is a 3×1 `Vector` meaning `M[a]ᵀ a` (e.g. telescoping arm `M[Ω̇] = Cross(rho,e1)m(-1)`);
  `linearize` transposes the analogous block, `to_standard_form` does not. Document the shape.
- `_has_algebraic`/`_same` compare type and name only, so a `Scalar("g")` input and a
  `Scalar("g", attr=["Constant"])` parameter are the same symbol: `to_standard_form(…, inputs=
  [Scalar("g")])` on `L = … − m g e₃·x` returns `G[g] = (e3 m (-1) + e3)`. No symbol table
  exists; at least `linearize._check_inputs` could reject an input whose name occurs in the
  equations as a constant leaf.
- `SystemVariables.vectors` accepts a tangent vector (`ω = q.get_tangent_vector()`) as a
  configuration variable and returns the meaningless equation `(ω m + f)`; a plain `Matrix`
  configuration reaches `NotImplementedError` in `extract_linear_coeff`. Both could be
  rejected in `_check_variables` (`TS2`/`TSO3` leaves; non-`SO3` matrices).
- `J(t)` symmetric but non-constant prints `(d/dt(J))'Ω` (the transpose of a derivative of a
  symmetric matrix is not simplified). Cosmetic.

## Areas probed that came back clean

All numbers are worst relative errors of the independent Euler–Lagrange check unless stated.

- **Pinned systems** (`tests/vbl_systems.py`): point mass R³ 2.8e-7, spring–mass R 8.7e-8,
  free rigid body 1.4e-7, rigid pendulum about pivot 7.3e-7, spherical pendulum 2.3e-6,
  quadrotor R³×SO(3) 2.2e-7, double rigid pendulum SO(3)×SO(3) 7.6e-7, elastic spherical
  pendulum R×S² 2.5e-7, quadrotor with cable load R³×S²×SO(3) 5.1e-7.
- **Multi-body products.** Double spherical pendulum S²×S² in position form and in ω-form
  (Lee et al. §5.6.5) 5.5e-6 / 3.7e-6, also with the colliding names `q`, `q1` (5.5e-6:
  `ξ_{q}` vs `ξ_{q1}` are routed correctly). Load with a two-link cable R³×S²×S² 4.8e-6.
  SO(3)×SO(3) with potential `k (R₁e₃)·(R₂e₃)` and velocity coupling `k (R₁Ω₁)·(R₂Ω₂)` 7.2e-7.
  S²×SO(3) with `k q̇·(R e₃)` 9.6e-7. Cross-terms between bodies, off-diagonal mass blocks
  (`M[ξ₁][ω̇₂] = Hat(q₁)Hat(q₂) m₂ l₁ l₂`, `M[x_L][ω̇] = −m_Q l Hat(q)`) all consistent.
- **Non-constant parameters.** `m(t)` in R³ (`ṁ ẋ` term present) 8.9e-7; `J(t)` symmetric on
  SO(3) (`½ J̇Ω + ½ J̇ᵀΩ` present) 6.0e-7; `K(t)` non-symmetric 1.9e-7; `K` constant
  non-symmetric (`½(K+Kᵀ)` structure) 1.8e-7; `l(t), m(t)` spherical pendulum 1.6e-7;
  `ρ(t)` on SO(3) 1.6e-6; python float/int constants mixed (`Dot(v,v)*m*0.5`, `*1`, `2.0*0.25`)
  4.6e-7. IBP of `⟨η̇, J w⟩` with non-constant J produces the `J̇` term (no constancy assumed
  without the flag).
- **Virtual work.** Moment `Dot(ξ, f)` vs force at the bob `Dot(δq·l, F)` (→ `l q × F`) 1.4e-6;
  configuration-dependent `Dot(δq, q × B)` 2.3e-6; Rayleigh-like `−k ξ·ω` 4.3e-8; Lorentz
  `Dot(δx, ẋ × B)` 3.9e-7; `Dot(F, δx)` order, `Dot(δx, F+B) k`, input appearing in both L and
  dW (`u x·e₃` + `u e₃·δx` → `e3 u (2)`) 4.1e-7. Nonlinear dW (`Dot(δx,F) Dot(δx,B)`, `δs δs`)
  raises `NonLinearError` as it should; `dW = Zero` fine.
- **to_standard_form.** `Σ M a + f + G u = E` to ≤ 4e-16 on: all pinned systems, double
  spherical pendulum, R×SO(3) telescoping arm (`x = R(ρ + s e₁)`, scalar row with vector
  acceleration), `m(t)` point mass (`ṁ ẋ` lands in `f`), R³×S² with a `(ẋ·q)²` term
  (`M[ẍ] = −m I − 2k q qᵀ` rank-one block), input under `Cross` (`G[F] = Hat(x)`), equation
  with no acceleration (`M = {}`). Nonlinear inputs `u u`, `F (F·F)` raise `NonLinearError`.
  `f` never contained an acceleration (`_has_algebraic` partition). Sign convention
  `M = −(physical inertia)` matches the docstring `M a + f + G u = 0`.
- **linearize** (finite-difference check ≤ 6e-10): double spherical pendulum (two kinematic
  rows, block routing by index), R³×S²×S² two-link cable, input in the Lagrangian (scalar
  `u x·e₃`, vector `F·x`), vector input under `Cross` with S² (`ξ·(q × f)`, `ξ·(f × ω)`, force at
  the bob), elastic pendulum in q̇-form (R×S²), SO(3)×SO(3) coupled potential, S²×SO(3)
  velocity coupling, S² damping, `u s²` in L with scalar config. `δω` appearing both as
  `Variation(ω)` and inside `d/dt(Variation(ω))` is partitioned correctly in every S² case.
- **Degenerate inputs.** L independent of a listed variable → `0v` equation, empty `M`, zero
  VBL row (no crash); empty `SystemVariables` → `{}`; zero Lagrangian → `E = F`; constant
  vector as configuration → `ValueError`; `S2` in `matrices` → `TypeError`; configuration `x`
  used as a direction in another body's potential → consistent equations (`−m g q` / `−m g q×x`).
- **Name collisions.** `q`/`q1`, `R1`/`R2`, `m1`/`m2`/`l1` all distinct through `str`-keyed
  dicts; `Mul(m, l)` (`"ml"`) never matches a leaf `Scalar("ml")` (type guard in `_same`).
- **Printing.** `print_eom`, `eom_to_latex`, `render_eom` (pdflatex), `str(LinearizedSystem)`,
  `repr(StandardFormEquation)`, `tree_str` run without error on all systems; only the
  parenthesisation issues of F8 were found.

## Appendix A — independent Euler–Lagrange harness (`specs/reviews/scripts/euler_lagrange_check.py`)

```python
"""E_lib(q, v, v̇) must equal −g + forcing, g the finite-difference Euler–Lagrange residual."""
import numpy as np
from geomech import S2, SO3, Scalar, Vector, Matrix, TimeDerivative
from tests.numeric_helpers import Evaluator, hat, vee, rand_rot

def cross(a, b): return np.cross(np.asarray(a).reshape(3), np.asarray(b).reshape(3)).reshape(3, 1)
def _dkey(sym): return str(TimeDerivative(sym))

class Coords:
    """Numeric coordinates (state), their d/dt (vel) and the library env built from them."""
    def __init__(self, variables, rng, extra_params=()):
        self.vars = list(variables.scalars) + list(variables.vectors) + list(variables.matrices)
        self.state, self.vel, self.extra = {}, {}, {}
        for v in self.vars:
            if isinstance(v, S2):
                q = rng.normal(size=(3, 1)); q /= np.linalg.norm(q); P = np.eye(3) - q @ q.T
                w, wd = P @ rng.normal(size=(3, 1)), P @ rng.normal(size=(3, 1))
                qd, qdd = cross(w, q), cross(wd, q) + cross(w, cross(w, q))
                self.state[("q", v)], self.state[("qd", v)] = q, qd
                self.vel[("q", v)], self.vel[("qd", v)] = qd, qdd; self.extra[v] = (w, wd)
            elif isinstance(v, SO3):
                R, Om, Omd = rand_rot(rng), rng.normal(size=(3, 1)), rng.normal(size=(3, 1))
                self.state[("R", v)], self.state[("Om", v)] = R, Om
                self.vel[("R", v)], self.vel[("Om", v)] = R @ hat(Om), Omd
            elif isinstance(v, Scalar):
                s, sd, sdd = (float(rng.normal()) for _ in range(3))
                self.state[("s", v)], self.state[("sd", v)] = s, sd; self.vel[("s", v)], self.vel[("sd", v)] = sd, sdd
            else:
                x, xd, xdd = (rng.normal(size=(3, 1)) for _ in range(3))
                self.state[("x", v)], self.state[("xd", v)] = x, xd; self.vel[("x", v)], self.vel[("xd", v)] = xd, xdd
        for p in extra_params:   # prescribed p(t): value, ṗ, p̈
            shape = () if isinstance(p, Scalar) else ((3, 3) if isinstance(p, Matrix) else (3, 1))
            a, b, c = (rng.normal(size=shape) if shape else float(rng.normal()) for _ in range(3))
            if isinstance(p, Matrix) and p.is_symmetric: a, b, c = (m + m.T for m in (a, b, c))
            self.state[("p", p)], self.state[("pd", p)] = a, b; self.vel[("p", p)], self.vel[("pd", p)] = b, c

    def env(self, state):
        env = {}
        for v in self.vars:
            if isinstance(v, S2):
                q, qd = state[("q", v)], state[("qd", v)]; env[str(v)] = q
                env[str(v.get_tangent_vector())] = cross(q, qd)          # ω = q × q̇
            elif isinstance(v, SO3):
                env[str(v)], env[str(v.get_tangent_vector())] = state[("R", v)], state[("Om", v)]
            elif isinstance(v, Scalar):
                env[str(v)], env[_dkey(v)] = state[("s", v)], state[("sd", v)]
            else:
                env[str(v)], env[_dkey(v)] = state[("x", v)], state[("xd", v)]
        for (kind, p), val in state.items():
            if kind == "p": env[str(p)] = val
            if kind == "pd": env[_dkey(p)] = val
        return env

    def full_env(self):   # base point including accelerations, for evaluating E
        env = self.env(self.state)
        for v in self.vars:
            if isinstance(v, S2): env[_dkey(v.get_tangent_vector())] = self.extra[v][1]
            elif isinstance(v, SO3): env[_dkey(v.get_tangent_vector())] = self.vel[("Om", v)]
            elif isinstance(v, Scalar): env[_dkey(v.t_diff())] = self.vel[("sd", v)]
            else: env[_dkey(v.t_diff())] = self.vel[("xd", v)]
        for (kind, p), val in self.vel.items():
            if kind == "pd": env[_dkey(p.t_diff())] = val
        return env

def _shift(state, key, idx, h):
    s = dict(state); val = state[key]
    if np.isscalar(val): s[key] = val + h
    else: val = np.array(val, dtype=float, copy=True); val.flat[idx] += h; s[key] = val
    return s

def _grad(fun, state, key, h=1e-5):
    val = state[key]
    if np.isscalar(val): return (fun(_shift(state, key, 0, h)) - fun(_shift(state, key, 0, -h))) / (2 * h)
    g = np.zeros_like(np.asarray(val, dtype=float))
    for i in range(g.size): g.flat[i] = (fun(_shift(state, key, i, h)) - fun(_shift(state, key, i, -h))) / (2 * h)
    return g

def _along(state, vel, h): return {k: (v + h * vel[k] if k in vel else v) for k, v in state.items()}
def _ddt(fun, state, vel, h=1e-5): return (fun(_along(state, vel, h)) - fun(_along(state, vel, -h))) / (2 * h)

def check_system(L, dW, variables, inputs, eom, trials=3, seed=0, extra_params=(), verbose=True):
    rng = np.random.default_rng(seed); worst = 0.0
    for t in range(trials):
        co = Coords(variables, rng, extra_params); ev = Evaluator(env=co.full_env(), rng=rng)
        for u in inputs: ev(u)
        for _, (_, eq) in eom.items(): ev(eq)       # draw every constant leaf once
        ev(L)
        var_keys = {w: (str(w.delta()) if not isinstance(w, (S2, SO3)) else str(w.get_variation_vector())) for w in co.vars}
        zero = {k: (0.0 if isinstance(w, Scalar) else np.zeros((3, 1))) for w, k in var_keys.items()}
        if hasattr(dW, "type"):
            e = dict(ev.env); e.update(zero); Evaluator(env=e, rng=rng)(dW)
            ev.env.update({k: v for k, v in e.items() if k not in ev.env and k not in zero})
        consts = dict(ev.env)
        def Lnum(state):
            e = dict(consts); e.update(co.env(state)); return float(Evaluator(env=e, rng=rng)(L))
        for v in co.vars:
            if isinstance(v, S2):
                dLdqd = lambda st: _grad(Lnum, st, ("qd", v))
                g = _ddt(dLdqd, co.state, co.vel) - _grad(Lnum, co.state, ("q", v))
                q = co.state[("q", v)]; Etrue = -cross(q, g); proj = np.eye(3) - q @ q.T
            elif isinstance(v, SO3):
                dLdOm = lambda st: _grad(Lnum, st, ("Om", v))
                R, Om = co.state[("R", v)], co.state[("Om", v)]; A = R.T @ _grad(Lnum, co.state, ("R", v))
                g = _ddt(dLdOm, co.state, co.vel) - cross(dLdOm(co.state), Om) - vee(A - A.T)
                Etrue = -g; proj = np.eye(3)
            elif isinstance(v, Scalar):
                dLdsd = lambda st: _grad(Lnum, st, ("sd", v))
                Etrue = -(_ddt(dLdsd, co.state, co.vel) - _grad(Lnum, co.state, ("s", v))); proj = 1.0
            else:
                dLdxd = lambda st: _grad(Lnum, st, ("xd", v))
                Etrue = -(_ddt(dLdxd, co.state, co.vel) - _grad(Lnum, co.state, ("x", v))); proj = np.eye(3)
            key = var_keys[v]
            if hasattr(dW, "type"):       # forcing = dW with this variation set to basis vectors
                basis = [1.0] if isinstance(v, Scalar) else [np.eye(3)[:, [k]] for k in range(3)]
                f = []
                for b in basis:
                    e = dict(consts); e.update(zero); e[key] = b; f.append(float(Evaluator(env=e, rng=rng)(dW)))
                Etrue = Etrue + (f[0] if isinstance(v, Scalar) else np.array(f).reshape(3, 1))
            Elib = Evaluator(env=dict(consts), rng=rng)(eom[key][1])
            if np.isscalar(proj):
                err, scale = abs(float(Elib) - float(Etrue)), max(1.0, abs(float(Etrue)))
            else:
                d = proj @ (np.asarray(Elib, dtype=float).reshape(3, 1) - Etrue)
                err, scale = float(np.max(np.abs(d))), max(1.0, float(np.max(np.abs(proj @ Etrue))))
            worst = max(worst, err / scale)
            if verbose: print(f"  trial {t} {key}: rel err {err / scale:.2e}")
    return worst
```

Usage: `eom = compute_eom(L, dW, variables); check_system(L, dW, variables, inputs, eom)`.
