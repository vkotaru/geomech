# Numeric verification

`geomech.numeric` evaluates expression trees with numpy and checks the output of
`compute_eom` against a finite-difference Euler–Lagrange residual that does not use the
library's variational rules. It is the regression gate for the engine and the quickest way
to gain confidence in a derivation for a new system.

```python
from geomech import evaluate, check_eom, assert_eom, Evaluator
```

## Leaf values: the key convention

Every leaf of a tree is identified by its printed name, `str(leaf)`. `evaluate` takes a
dictionary from those names to numpy values and returns a float (scalar expression), a
`(3, 1)` array (vector) or a `(3, 3)` array (matrix):

| Leaf | Key | Value |
|---|---|---|
| `Scalar("m")` | `"m"` | `float` |
| `Vector("x")` | `"x"` | `(3, 1)` array |
| `Matrix("J")` | `"J"` | `(3, 3)` array |
| velocity `x.t_diff()` | `str(TimeDerivative(x))` | `(3, 1)` array |
| acceleration `x.t_diff().t_diff()` | `str(TimeDerivative(TimeDerivative(x)))` | `(3, 1)` array |
| `q = S2("q")` | `"q"` | unit `(3, 1)` array |
| `omega = q.get_tangent_vector()` | `str(omega)` | `(3, 1)` array $\perp q$ |
| `R = SO3("R")` | `"R"` | rotation matrix |
| `Omega = R.get_tangent_vector()` | `str(Omega)` | `(3, 1)` array |

Velocities and tangent vectors are leaves of their own: the tree never holds
`TimeDerivative(q)` for a manifold point (`q.t_diff()` *is* `Cross(omega, q)`), and $\dot x$ is a
value you supply, not something computed from $x$. The zero and identity constants are
resolved by type and need no key.

Any leaf missing from the dictionary is drawn once from `numpy.random.default_rng(seed)`,
respecting the leaf's declared structure (`UnitNorm` vectors have unit norm,
`SymmetricMatrix` leaves are symmetric, `SO3` leaves are rotations, a `Scalar` with a
`value` takes that value). A `Variation` symbol has no natural random value and must be
given explicitly.

```python
import numpy as np
from geomech import Vector, Cross, Dot, evaluate

x, e3 = Vector("x"), Vector("e3", attr=["Constant"])
v = x.t_diff()
env = {"x": np.array([[1.0], [0.0], [0.0]]),
       str(v): np.array([[0.0], [2.0], [0.0]]),
       "e3": np.array([[0.0], [0.0], [1.0]])}
evaluate(Dot(Cross(x, v), e3), env)   # 2.0
```

`Evaluator(env, rng)` is the underlying class; it stores every value it draws in `env`, so
successive calls on the same instance see the same numbers.

## Checking an equation of motion

```python
from geomech import (S2, Vector, Dot, Scalar, SystemVariables, getScalars,
                     compute_eom, check_eom)

m, g, l = getScalars("m g l", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
q, M = S2("q"), Vector("M")      # M: a moment about the pivot
omega, xi = q.get_tangent_vector(), q.get_variation_vector()
half = Scalar("0.5", value=0.5, attr=["Constant"])

L = half * m * l * l * Dot(omega, omega) - m * g * l * Dot(q, e3)
dW = Dot(xi, M)
variables = SystemVariables(vectors=[q])

eom = compute_eom(L, dW, variables)
err = check_eom(L, dW, variables, eom, inputs=[M])
print(err)          # ~1e-7: the equation −m l² ω̇ − m g l (q × e₃) + M = 0 is right
assert err < 1e-5
```

`check_eom` returns the worst relative error over `trials` random trajectory points and
over all configuration variables; `assert_eom(..., rtol=1e-5)` raises an `AssertionError`
naming the variation vector whose equation is wrong. For a system whose parameters are of
order one a correct equation scores at finite-difference noise, about `1e-7` (worst
`2.3e-7` measured over 36 pinned systems and 10 seeds), and a wrong sign or a dropped term
scores `0.02` or more (worst over every term of the same systems). The gap is four orders
of magnitude for such systems; it is **not** unconditional, see "Resolution" below.

### What it computes

`compute_eom` returns, for each configuration variable, the coefficient $E$ of its
variation vector in $\delta\!\int L\,dt + \int \delta W\,dt$ after integration by parts,
so the equation of motion is $E = 0$ and $E = -(\text{Euler–Lagrange}) + \text{forcing}$
(for a point mass, $E = -m\ddot x + f$). `check_eom` evaluates the *numeric* Lagrangian
with `Evaluator` and forms the Euler–Lagrange residual $g$ by central differences:

- $\mathbb{R}$, $\mathbb{R}^3$: $g = \frac{d}{dt}\frac{\partial L}{\partial v} -
  \frac{\partial L}{\partial q}$, $E_{\text{true}} = -g + f$;
- $S^2$: the constrained Lagrangian $\tilde L(q, \dot q) = L(q, q \times \dot q)$ on
  $\lvert q \rvert = 1$ with $\dot q = \omega \times q$, $E_{\text{true}} = -q \times g + f$,
  compared **on the tangent plane only** (Lee et al., Chapter 5);
- $SO(3)$: the body-frame form $g = \frac{d}{dt}\frac{\partial L}{\partial \Omega} -
  \frac{\partial L}{\partial \Omega} \times \Omega - (A - A^T)^\vee$ with
  $A = R^T \frac{\partial L}{\partial R}$, $E_{\text{true}} = -g + f$
  (Lee et al., Chapter 6).

The forcing $f$ is `inf_work` evaluated with the variable's variation set to each basis
vector and the other variations zero. $d/dt$ is a directional derivative along the sampled
velocity and acceleration, so nothing is integrated: the library's $E(q, v, \dot v, u)$ and
$E_{\text{true}}$ must agree at random accelerations, which pins the mass matrix, the
velocity terms and the forces at once.

### Resolution: what a term must be for the check to see it

The error is measured relative to the size of the whole equation,
$\max(1, \lvert P\,E_{\text{true}} \rvert)$. A term $t$ that is wrong (sign, missing factor) moves
the error by $2\lvert P\,t \rvert / \max(1, \lvert P\,E_{\text{true}} \rvert)$, so it is **resolved** at `rtol`
only if that quantity exceeds `rtol`. With parameters of order one every term is resolved
by a wide margin. With a large parameter in the equation it is not: a point mass on a
spring with `Scalar("k", value=1e7)` has $\lvert E_{\text{true}} \rvert \approx 10^7$, the gravity
term $m g\, e_3$ is of order one, a flipped gravity term moves the error by
$\approx 4 \cdot 10^{-7}$, and `check_eom` returns `~1e-9` for the correct and
`~1e-7` for the wrong equation: both pass at `1e-5`. The finite-difference noise scales
with $\lvert L \rvert$ as well (it is `6e-3` to `2e-2` in absolute terms for that system at the
default step), so this is a limit of the finite-difference check, not of the
normalisation alone.

`assert_eom` therefore checks the resolution of every additive term of each equation
(its tangent part on $S^2$) and raises, naming the term, its size and the resolution
$\mathrm{rtol} \cdot \max(1, \lvert E_{\text{true}} \rvert)$, when a term is unresolved at every trial. Terms parallel to
$q$ or zero on the constraint set have nothing to resolve and are ignored. If you need
to accept such an equation on the relative error alone, pass `require_resolved=False`;
the better fix is to keep the parameters symbolic (random, order one) or to scale them so
the terms are of comparable size. `check_eom` reports the relative error only.

### What it does not verify

- **The normal component on $S^2$.** The coefficient of $\xi \perp q$ is only defined on
  the tangent plane, so a term such as $\omega \times (q \times \dot x) = q\,(\omega\cdot\dot x)$
  is invisible to the check. This is correct behaviour, not a gap: such a term has no
  physical content. The same applies to the equation itself: on $S^2$ the equation
  $E = 0$ returned by `compute_eom` is only meaningful on the tangent plane,
  $(I - q q^T)E = 0$. $E$ generally has a component along $q$ (from the forcing,
  from $\omega \times c$, from $A\dot\omega$, …), so anyone solving $E = 0$ for
  $\dot\omega$, including from the blocks of `to_standard_form` or `linearize`, must
  project: solve $(I - q q^T)(M\dot\omega + f + G u) = 0$ with $\dot\omega \perp q$.
  The unprojected $3 \times 3$ solve gives $q \cdot \dot\omega \neq 0$ and a trajectory
  that leaves the sphere.
- **Anything below finite-difference noise.** $\frac{d}{dt}\frac{\partial L}{\partial v}$ nests two central differences
  with the same step, so the noise is about $h^2 + 10^{-16} \lvert L \rvert / h^2$: `~1e-8` for a
  Lagrangian of order one at the default `h=1e-4` (the minimum of the two terms; at
  `h=1e-5` the roundoff term is `~1e-6`, within a factor 10 of the gate, and a correct
  equation can fail it). A term whose tangent part is below the noise cannot be verified
  by any normalisation; this grey zone always exists, and it widens with $\lvert L \rvert$ (see
  "Resolution").
- **Unlisted time-varying parameters.** A non-constant leaf that is neither a configuration
  variable nor listed in `params` is frozen during the check, so a correct equation that
  contains its rate looks wrong. List it: `check_eom(..., params=[l])`. `compute_eom`
  takes the same `params=[l]`; without it, the surviving $\delta l$ is reported as an error
  rather than dropped, so the two calls declare the prescribed symbols identically.
- **The symbolic form.** Two equations that agree numerically are accepted even if one is
  not simplified; `check_eom` checks values, the pinned string tests check form.
- **`to_standard_form` and `linearize`.** The standard-form split is a rearrangement of
  the same expression; linearizations are checked separately by
  `geomech.numeric.finite_difference_check`.

## Using it on a new system

1. Build `L`, `dW` and `SystemVariables` as usual and call `compute_eom`.
2. Call `check_eom(L, dW, variables, eom, inputs=[...])`. List every prescribed
   time-varying symbol in `params`.
3. If the error is `~1e-7`, the equation is right. If it is `~1`, use `assert_eom` to
   learn which variable's equation is wrong, then compare `print_eom(eom)` with the
   textbook form; a wrong equation is an engine bug worth reporting with the Lagrangian
   that triggers it. If `assert_eom` reports an unresolved term, the check could not
   see that term: use parameters of order one.
4. Pin the system in a test with `assert_eom(...)` (or `assert check_eom(...) < 1e-5`
   when the parameters are of order one) next to the string assertion on the equation,
   so the form and the value are both guarded.
