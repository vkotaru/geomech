# geomech

A Python library for the symbolic derivation of coordinate-free equations of motion, and their variation-based linearization, on $R$, $R^3$, $S^1$, $S^2$, $SO(3)$ and their products (e.g. $SE(3)$ as $R^3 \times SO(3)$), including models with a symbolic number of bodies.

`geomech` implements the variational method of Lee, Leok and McClamroch,
[Global Formulations of Lagrangian and Hamiltonian Dynamics on Manifolds](https://link.springer.com/book/10.1007/978-3-319-56953-6).
Bittner and Sreenath were the first to automate the derivation of coordinate-free
equations of motion on manifolds, in
[Symbolic Computation of Dynamics on Smooth Manifolds](https://link.springer.com/chapter/10.1007/978-3-030-43089-4_22)
and its [Scala package](https://github.com/HybridRobotics/dynamics_on_manifolds).

Documentation: https://geomech.readthedocs.io

## What it does

Given a Lagrangian and virtual work, `geomech` automatically derives
globally valid, singularity-free equations of motion on smooth manifolds
(SO(3), S2, S1) using Hamilton's principle and variational calculus --
no Euler angles or quaternions required.

The pipeline:

1. Define symbolic scalars, vectors, and matrices
2. Build kinetic/potential energy and virtual work expressions
3. Call `compute_eom()` -- the library takes the variation of the action, applies the
   manifold kinematics, simplifies, integrates by parts and extracts the equations
   ([pipeline, stage by stage](https://geomech.readthedocs.io/en/latest/pipeline.html))
4. Optionally call `to_standard_form()` for `M(q) a + f + G u = 0`, or
   `linearize()` for the variation-based linearization `M δẋ = A δx + B δu`
   along a trajectory
5. Optionally generate NumPy (or CasADi) code for the result with `geomech.codegen`
   (`generate_module`, `lambdify`, `energy_function`; [code generation](https://geomech.readthedocs.io/en/latest/codegen.html));
   `examples/generated/` holds the modules of nine systems
6. Optionally simulate a generated model on its manifold with `geomech.simulation`
   (`simulate` with a Lie-group RK4 integrator or `solve_ivp`, energy checks, and `lqr`
   from the linearization; [simulation](https://geomech.readthedocs.io/en/latest/simulation.html),
   `examples/inverted_rigid_pendulum_lqr.py`)
7. Optionally draw the simulated system in 3-D and write one self-contained HTML report
   (animation, equations, standard form, verification numbers) with `geomech.viz`
   ([visualization](https://geomech.readthedocs.io/en/latest/visualization.html), `examples/report_*.py`); an agent skill that
   follows this whole procedure is in `skills/geomech/`
8. For a model of n bodies (a chain, a cable of n links, a team of n quadrotors), write
   the Lagrangian with `Sum` over a `Family` of variables and call
   `geomech.indexed.compute_indexed_eom`: it derives the equation of the i-th body with
   i and n symbolic, e.g. the inertia of an n-link chain as
   `Σ_j Σ_{k ≥ max(i, j)} m_k l_i l_j (...)` ([sums over n bodies](https://geomech.readthedocs.io/en/latest/indexed.html))

## How it works

Expressions are **symbolic trees**. Python operators (`*`, `+`, etc.) build
the tree; transformations rewrite it using `match`/`case` pattern matching.

**Node arity:**

| Arity  | Nodes | Example |
| ------ | ----- | ------- |
| N-ary  | `Add`, `VAdd`, `MAdd` | `a + b + c` → one node, three children |
| Binary | `Mul`, `Dot`, `Cross`, `SVMul`, `MVMul`, ... | `Dot(x, y)` → two children |
| Unary  | `Hat`, `Vee`, `Transpose`, `TimeDerivative`, `Variation` | `Hat(x)` → one child |

**Sample trees:**

`M * (x + y + Cross(m, n))` -- matrix-vector multiply with N-ary addition:

```
          MVMul                   ← Binary (matrix × vector → vector)
          ╱   ╲
         M    VAdd                ← N-ary (3 children)
            ╱  |  ╲
           x   y  Cross           ← Binary (vector × vector)
                  ╱   ╲
                 m     n
```

`KE = m * Dot(ẋ, ẋ) * 0.5` -- kinetic energy:

```
              Mul                 ← Binary (scalar × scalar)
             ╱   ╲
          Mul     0.5
         ╱   ╲
        m    Dot                  ← Binary (vector · vector → scalar)
            ╱   ╲
         d/dt   d/dt             ← Unary (time derivative)
          |      |
          x      x               ← Vector leaves
```

**Tree rewrites** -- one function per node class, picked by the class of the node
(`match`/`case` inside a handler where a rule depends on the operands' classes):

```python
def _elim_dot(expr):                         # binary
    l, r = _eliminate(expr.left), _eliminate(expr.right)
    if _is_zero(l) or _is_zero(r):
        return rule("simplify.dot_zero", expr, Zero)
    return Dot(l, r)


def _elim_hat(expr):                         # unary
    inner = _eliminate(expr.expr)
    if isinstance(inner, Vee) and _is_skew_symmetric(inner.expr):
        return rule("simplify.hat_vee_skew", expr, inner.expr)  # Hat(Vee(X)) = X
    return Hat(inner)


_ELIMINATE = {Dot: _elim_dot, Hat: _elim_hat, Add: _elim_add, ...}  # a leaf has none


def _eliminate(expr):
    handler = _ELIMINATE.get(type(expr))
    return _canonical_literal(expr) if handler is None else handler(expr)
```

Nodes are frozen: a rewrite always builds new nodes, and a subtree can be shared by
several trees.

Every rewrite is logged through `rule(name, before, after)`; the names, with the hypothesis
each rule needs, are listed in the module docstrings, and `GEOMECH_TRACE=1` logs every
firing.

The full pipeline (`compute_eom`), each stage a recursive tree rewrite:

δS = δL + δW → manifold rules (δω, δΩ) → `full_simplify` → integration by parts
(`expand`, `full_simplify`, then per variation vector: `collect` → IBP → `full_simplify`) →
`expand` → extract the coefficient of each variation vector → `full_simplify`.

[Deriving equations](https://geomech.readthedocs.io/en/latest/pipeline.html) has one paragraph per stage;
`compute_eom(..., trace=EOMTrace())` keeps every stage, `GEOMECH_DEBUG=1` logs them, and
`GEOMECH_TRACE=1` logs every named rewrite rule that fires
([Troubleshooting](https://geomech.readthedocs.io/en/latest/troubleshooting.html)).

## Installation

Requires Python >= 3.10.

```bash
pip install geomech
pip install "geomech[casadi]"   # optional: the CasADi code generation target
pip install "geomech[sim]"      # optional: scipy, for solve_ivp and lqr in geomech.simulation
```

From a clone, for development: `pip install -e ".[dev]"`.

## Quick start

**Point mass** -- derives `m * ddot_x = f - m*g*e3`:

```python
from geomech import *
import numpy as np

m, g = getScalars('m g', attr=['Constant'])
e3 = Vector('e3', attr=['Constant'], value=np.array([0., 0., 1]))
x, f = getVectors(['x', 'f'])

v = x.t_diff()

KE = m * Dot(v, v) * 0.5
PE = m * x.dot(g * e3)
L  = KE - PE

dW = Dot(x.delta(), f)

eqs = compute_eom(L, dW, SystemVariables(vectors=[x]))
print_eom(eqs)
```

**Spherical pendulum** on S2:

```python
m, g, l = getScalars('m g l', attr=['Constant'])
e3 = Vector('e3', attr=['Constant'], value=np.array([0., 0., 1]))

q = S2('q')             # attitude on the 2-sphere
f = Vector('f')

x = l * q
v = x.t_diff()

KE = m * Dot(v, v) * 0.5
PE = m * x.dot(g * e3)
L  = KE - PE
dW = Dot(q.delta(), f)

eqs = compute_eom(L, dW, SystemVariables(vectors=[q]))
print_eom(eqs)
```

**Rigid body** on SO(3):

```python
J = Matrix('J', attr=['Constant', 'SymmetricMatrix'])
rho = Vector('\\rho', attr=['Constant'])
m, g = getScalars('m g', attr=['Constant'])
e3 = Vector('e_3', attr=['Constant'], value=np.array([0., 0., 1]))

R = SO3('R')             # rotation matrix on SO(3)
Om = R.get_tangent_vector()
eta = R.get_variation_vector()
M = Vector('M')          # external torque

x = R * rho
v = x.t_diff()

KE = Dot(Om, J * Om) * 0.5 + Dot(v, v) * m * 0.5
PE = m * g * Dot(x, e3)
L  = KE - PE
dW = Dot(eta, M)

eqs = compute_eom(L, dW, SystemVariables(matrices=[R]))
print_eom(eqs)
```

More examples in [`examples/`](https://github.com/vkotaru/geomech/blob/main/examples/).

## Using geomech with AI coding agents

The repository includes an [Agent Skills](https://agentskills.io/specification) skill,
[`skills/geomech`](skills/geomech), that walks an agent through deriving, verifying,
linearizing, simulating and visualizing a system with this package. The same folder
works in Claude Code, Codex and Gemini CLI. Copy it (or symlink it) into the directory
your agent reads:

| Agent | In one project | For all your projects |
| ----- | -------------- | --------------------- |
| Claude Code | `.claude/skills/geomech/` | `~/.claude/skills/geomech/` |
| Codex | `.agents/skills/geomech/` | `~/.agents/skills/geomech/` |
| Gemini CLI | `.agents/skills/geomech/` (or `.gemini/skills/`) | `~/.agents/skills/geomech/` (or `~/.gemini/skills/`) |

For example, from the root of your project:

```bash
git clone --depth 1 https://github.com/vkotaru/geomech /tmp/geomech
mkdir -p .agents/skills && cp -r /tmp/geomech/skills/geomech .agents/skills/
pip install geomech
```

Then ask the agent, for instance, to "derive the equations of motion of a double
spherical pendulum and check them". The skill's `Setup` section installs the package if
it is missing.

Agents working inside this repository are pointed to the skill by
[`AGENTS.md`](AGENTS.md) (`CLAUDE.md` and `GEMINI.md` import it).

## Expression types

| Type     | Description                 | Constructor                      |
| -------- | --------------------------- | -------------------------------- |
| `Scalar` | Symbolic scalar             | `Scalar('m', attr=['Constant'])` |
| `Vector` | Column vector (default 3x1) | `Vector('x')`                    |
| `Matrix` | Square matrix (default 3x3) | `Matrix('J', attr=['Constant'])` |
| `S1`     | Unit vector on the circle normal to a constant unit axis | `S1('q', axis=e)` |
| `S2`     | Unit vector on the 2-sphere | `S2('q')`                        |
| `SO3`    | Rotation matrix on SO(3)    | `SO3('R')`                       |

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
x.t_diff()       # time derivative
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

## Package structure

```
geomech/
    core/
        base/           # Expr types, flags, scalars/vectors/matrices, manifolds
        operations/     # Add, Mul, Dot, Cross, Hat, Variation, TimeDerivative, ...
        transformations/# expand, simplify, pull, vector_rules
        math/           # extract_linear_coeff, collect, integrate_by_parts
    dynamics/           # compute_eom, to_standard_form, linearize, SystemVariables
    indexed/            # Family, Sum, compute_indexed_eom: models with n bodies, n symbolic
    systems/            # library of 30 verified systems (get, list_systems)
    codegen/            # StateLayout, to_source, lambdify, generate_module, energy_function
    simulation/         # accelerations, simulate (RKMK4, solve_ivp), Trajectory, lqr
    viz/                # Scene (3-D bodies from expressions), Report (HTML), latex_names
    numeric.py          # evaluate, check_eom, finite_difference_check
    utils/              # print_eom, print_tree, render_eom
```

## Running tests

```bash
pip install -e ".[dev]"
pytest
```

## Citation

If you use geomech in your research, please cite this package and the book of Lee, Leok and McClamroch whose method it implements, [Global Formulations of Lagrangian and Hamiltonian Dynamics on Manifolds](https://link.springer.com/book/10.1007/978-3-319-56953-6). The first automated derivation of these equations is described in [Symbolic Computation of Dynamics on Smooth Manifolds](https://link.springer.com/chapter/10.1007/978-3-030-43089-4_22).

```bibtex
@software{kotaru2026geomech,
  title={geomech: Coordinate-free equations of motion on manifolds for robotics},
  version={1.0.0},
  author={Kotaru, Prasanth},
  url={https://github.com/vkotaru/geomech},
  year={2026}
}
```

## License

BSD 3-Clause. See [LICENSE](LICENSE).
