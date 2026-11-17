# Installation and first system

## Installation

Requires Python >= 3.10.

```bash
pip install geomech
```

Optional extras: `geomech[sim]` (SciPy, for `solve_ivp` and `lqr`) and `geomech[casadi]`
(the CasADi code generation target). From a clone, for development:
`pip install -e ".[dev]"`.

## Your first derivation

Let's derive Newton's equation of motion for a 3D point mass.

```python
from geomech import *
import numpy as np

# Constants
m, g = getScalars('m g', attr=['Constant'])
e3 = Vector('e3', attr=['Constant'], value=np.array([0., 0., 1]))

# State variables
x, f = getVectors(['x', 'f'])
v = x.t_diff()

# Energies
KE = m * Dot(v, v) * 0.5
PE = m * x.dot(g * e3)
L  = KE - PE

# Virtual work
dW = Dot(x.delta(), f)

# Derive EOM
eqs = compute_eom(L, dW, SystemVariables(vectors=[x]))
print_eom(eqs)
```

This produces the equation $m \ddot{x} = f - mge_3$, output as a LaTeX integral equation.

Every symbol in `L` and `dW` is one of three things: a configuration variable (listed in
`SystemVariables`), a constant (`attr=["Constant"]`), or a prescribed time-varying
parameter declared with `compute_eom(L, dW, variables, params=[l])` — a cable length
$l(t)$, say, whose rate $\dot l$ appears in the equations but whose variation is not free. A
symbol that is none of these and whose variation survives (a forgotten `y`, an `e3`
missing `Constant`) raises `ValueError` naming it; the terms are never dropped silently.

## How it works

The `compute_eom` pipeline, in the order the code runs it:

$\delta S = \delta L + \delta W$ → manifold rules ($\delta\omega$, $\delta\Omega$) → `full_simplify` → integration by parts
(`expand`, `full_simplify`, then per variation vector: `collect` → IBP → `full_simplify`) →
`expand` → extract the coefficient of each variation vector → `full_simplify`.

1. **Variation** -- $\delta L$ by the product rule; manifold points are eager
   ($\delta R = R\widehat{\eta}$ for $SO(3)$, $\delta q = \xi \times q$ for $S^2$),
   and $\delta S = \delta L + \delta W$
2. **Manifold rules** -- $\delta\omega = \dot\xi - 2\,\omega \times \xi$ on $S^2$,
   $\delta\Omega = \dot\eta + \Omega \times \eta$ on $SO(3)$
3. **Simplify** -- expand once, then pull scalars, apply vector identities (orthogonality,
   unit norm, tangency) and eliminate zeros until nothing changes
4. **Integrate by parts** -- for each variation vector, `collect` rotates every term so its
   time derivative is the left operand of a `Dot`, then moves the derivative off it
5. **Expand** -- distributes products over sums
6. **Extract** -- the coefficient of each variation vector, simplified, is its equation

[Deriving equations](pipeline.md) describes each stage and how to see it
(`compute_eom(..., trace=EOMTrace())`); [Troubleshooting](troubleshooting.md) lists the
error messages and the debugging switches; [Conventions](conventions.md) the sign and
printing conventions.

## Expression types

| Type | Description | Constructor |
|------|-------------|-------------|
| `Scalar` | Symbolic scalar | `Scalar('m', attr=['Constant'])` |
| `Vector` | Column vector (default 3x1) | `Vector('x')` |
| `Matrix` | Square matrix (default 3x3) | `Matrix('J', attr=['Constant'])` |
| `S1` | Unit vector on the circle normal to the constant unit axis `e` (`attr=['Constant', 'UnitNorm']`); scalar $\omega$, $\gamma$ | `S1('q', axis=e)` |
| `S2` | Unit vector on the 2-sphere | `S2('q')` |
| `SO3` | Rotation matrix on $SO(3)$ | `SO3('R')` |

Attributes are passed as a list of strings, `attr=['Constant']`, and the names are
case-sensitive (`'constant'` raises `ValueError` listing the valid names).

Batch constructors: `getScalars('a b c')`, `getVectors(['x', 'y'])`, `getMatrices('M N')`.

## Operations

**Arithmetic** -- natural Python operators build expression trees:

```python
a * b          # scalar multiplication
a * x          # scalar-vector multiplication
M * x          # matrix-vector multiplication
x + y          # vector addition
Dot(x, y)      # dot product
Cross(x, y)    # cross product
Hat(x)         # hat map (R^3 -> so(3))
Vee(M)         # vee map (so(3) -> R^3)
Transpose(x)   # transpose (x.T() is the same; Transpose(Transpose(x)) is x)
R.T() * x      # a transposed matrix multiplies like a matrix: R^T x
x.T() * y      # row times column: the inner product x^T y; x * y.T() is the outer product
x.T() * K * x  # a row times a matrix is a row: the quadratic form x^T K x
# (also in a kinetic energy, v.T() * K * v: collect turns it into Dot(v, K * v) for integration by parts)
2 * x, a + 1, 1 - a, -x, x / 2   # numbers on either side (a number adds only to a scalar); / divides by a number only
```

A mistyped operand raises `ExpressionMismatchError` naming both operands, e.g.
`Dot cannot combine SCALAR 'm' and VECTOR 'x'`.

**Calculus**:

```python
x.t_diff()       # time derivative d/dt
x.t_integrate()  # time integral
x.delta()        # variation (delta)
```

**Transformations**:

```python
expand(expr)        # distribute products over sums
simplify(expr)      # pull scalars, apply identities, eliminate zeros
full_simplify(expr) # expand + simplify until convergence
```

**Printing** -- `str(expr)` is the internal form and the equality key (products print
`*` between their factors: `m*g*Dot(e3,x)`); for reading:

```python
to_latex(expr)      # LaTeX, parenthesised where juxtaposition would be ambiguous
print_eom(eqs)      # one LaTeX integral equation per variable (via to_latex)
repr(expr)          # compact: v:x vector, M:J matrix, S2:q, SO3:R; trailing * = constant;
                    # products joined by ·, e.g. m*·g*·<v:e3*, v:x>
print_tree(expr)    # one node per line with every leaf's flags (style="topdown" for a sketch)
```
