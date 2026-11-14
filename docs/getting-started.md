# Getting Started

## Installation

Requires Python >= 3.10.

```bash
pip install -e .
```

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

This produces the equation \(m \ddot{x} = f - mge_3\), output as a LaTeX integral equation.

Every symbol in `L` and `dW` is one of three things: a configuration variable (listed in
`SystemVariables`), a constant (`attr=["Constant"]`), or a prescribed time-varying
parameter declared with `compute_eom(L, dW, variables, params=[l])` — a cable length
`l(t)`, say, whose rate `l̇` appears in the equations but whose variation is not free. A
symbol that is none of these and whose variation survives (a forgotten `y`, an `e3`
missing `Constant`) raises `ValueError` naming it; the terms are never dropped silently.

## How it works

The `compute_eom` pipeline:

1. **Variation** -- takes \(\delta L\) by applying the product rule and manifold-specific variation formulas (\(\delta R = R\hat{\eta}\) for \(SO(3)\), \(\delta q = \xi \times q\) for \(S^2\))
2. **Expand** -- distributes products over sums
3. **Integrate by parts** -- moves time derivatives off variation vectors
4. **Simplify** -- pulls scalars, applies vector identities (orthogonality, unit norm), eliminates zeros
5. **Extract** -- collects coefficients of each variation vector to produce the final equations

## Expression types

| Type | Description | Constructor |
|------|-------------|-------------|
| `Scalar` | Symbolic scalar | `Scalar('m', attr=['Constant'])` |
| `Vector` | Column vector (default 3x1) | `Vector('x')` |
| `Matrix` | Square matrix (default 3x3) | `Matrix('J', attr=['Constant'])` |
| `S2` | Unit vector on the 2-sphere | `S2('q')` |
| `SO3` | Rotation matrix on \(SO(3)\) | `SO3('R')` |

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
Transpose(x)   # transpose
```

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
