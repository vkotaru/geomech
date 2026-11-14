# Code generation

`geomech.codegen` turns the symbolic output of the pipeline into plain numeric functions:
the blocks of `to_standard_form`, the equations of `compute_eom`, the blocks of
`linearize` and the total energy of a Lagrangian. Every generated function goes through
one code path, and the test suite checks it against `geomech.numeric.evaluate` on every
pinned system.

```python
from geomech.codegen import StateLayout, to_source, lambdify, generate_module, energy_function

layout = StateLayout(variables, inputs=[f, M], params=[m, J, g, e3])
src = to_source(expr, layout)                 # one expression -> Python source
fn = lambdify(expr, layout)                   # -> fn(state, u, p, a=None)
generate_module(system, layout, "model.py")   # writes a module, returns its source
E = energy_function(L, layout)                # -> E(state, p)
```

Two targets:

| Target | Dependency | Use |
|---|---|---|
| `"numpy"` (default) | `numpy` only | simulation (`scipy.integrate.solve_ivp`), tests, examples |
| `"casadi"` | extra `geomech[casadi]` | optimal control, automatic derivatives; the generated module imports `casadi`, `geomech` does not |

A generated module needs only `numpy` (or `casadi` and `numpy`), not `geomech`.

## The state layout

`StateLayout(variables, inputs=(), params=())` is the single source of truth for how the
symbols map to arrays, and is embedded in every generated module as `LAYOUT` (plain data,
`StateLayout.to_dict()`).

| Variable | Configuration slot | Velocity slot | Acceleration |
|---|---|---|---|
| `s: Scalar` (\(\mathbb{R}\)) | `s` (1) | `ṡ` (1) | `s̈` |
| `x: Vector` (\(\mathbb{R}^3\)) | `x` (3) | `ẋ` (3) | `ẍ` |
| `q: S2` | `q` (3, unit) | `ω` (3, \(\perp q\)) | `ω̇` |
| `R: SO3` | `R` (9, row-major) | `Ω` (3) | `Ω̇` |

- `state` holds every configuration slot in variable order (the `SystemVariables` order:
  scalars, vectors, matrices), then every velocity slot in the same order.
  `layout.slices[str(symbol)]` is the slice of a configuration or velocity symbol;
  `layout.accel_slices` and `layout.input_slices` those of `a` and `u`.
- `a` stacks the accelerations in variable order; it is an argument only of the functions
  whose expressions contain an acceleration (`compute_eom` equations, linearization blocks).
- `u` stacks the inputs in the order given (a scalar takes 1 entry, a vector 3).
- `p` is a `dict` keyed by `str(symbol)` of every parameter: a float for a scalar, any
  array of 3 or 3×3 entries for a vector or matrix. A missing key raises `KeyError`, an
  extra one `ValueError`, both naming the key. Numeric literals (`0.5`, `(-1)`) and the
  zero and identity constants are inlined, never parameters. A *named* scalar with a value
  (`Scalar("m", value=2)`) is a symbol, so it is a parameter.
- `collect_params(exprs, variables, inputs)` lists the symbols of the expressions that are
  not covered by the variables and inputs (sorted by name): the `params` to pass.
- `layout.state_from_env(env)`, `accel_from_env`, `inputs_from_env`, `params_from_env`
  convert an `Evaluator` env (printed-name keys, `(3, 1)` columns) to the flat arrays, and
  `env_from_arrays(state, a, u, p)` converts back.
- Each slot also has a Python identifier for the generated code (`\omega_{q}` → `omega_q`,
  `ω̇` → `omega_q_dot`, `ẋ` → `x_dot`, `ẍ` → `x_ddot`).

Generation checks every leaf: it must be a configuration, velocity, acceleration, input or
parameter of the layout, a numeric literal, or a zero or identity constant. A stray leaf,
a `TimeDerivative` of a non-leaf, a `Variation` or a time integral raises `ValueError`
naming it. Two layout symbols that print alike raise at construction.

## The generated module

For a standard-form system (`to_standard_form`):

```python
LAYOUT = {...}                              # StateLayout.to_dict()
def mass(state, p)         -> (n, n_a)      # block matrix of M[a_j] per equation
def bias(state, p)         -> (n,)          # stacked f
def input_matrix(state, p) -> (n, n_u)      # block matrix of G[u_k]
def kinematics(state)      -> (n_config,)   # ṡ, ẋ, ω × q, R hat(Ω) row-major
def energy(state, p)                        # with generate_module(..., lagrangian=L)
def layout()                                # a copy of LAYOUT
```

so that the equations read `mass @ a + bias + input_matrix @ u = 0`. Equations (block
rows) are in the order of the configuration variables; a scalar equation takes one row and
a vector equation three. For `compute_eom` output the module has
`residual(state, a, u, p)` (the stacked equations `E`) and `kinematics`. For a
`LinearizedSystem` it has `M_lin`, `A_lin`, `B_lin`, each `(state, a, u, p)`, with the rows
of `LinearizedSystem.rows` and the columns of `dstates`, `states` and `inputs`; the slices
are recorded in `LINEARIZATION`.

`examples/generated/` holds the modules of the nine systems of `tests/vbl_systems.py`
(`<name>.py` from the standard form with `energy`, `<name>_lin.py` from `linearize`).
`python scripts/regenerate_generated.py` rewrites them; CI fails when they drift.

## Worked example: the spherical pendulum

Lee et al. Eq. 5.14, \(L = \tfrac12 m l^2\,\omega\cdot\omega - m g l\, q\cdot e_3\) and
\(\delta W = \xi\cdot f\):

```python
import numpy as np
from geomech import S2, Dot, Scalar, SystemVariables, Vector, compute_eom, getScalars
from geomech import to_standard_form
from geomech.codegen import StateLayout, generate_module, load_module

m, g, l = getScalars("m g l", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
q, f = S2("q"), Vector("f")
omega, xi = q.get_tangent_vector(), q.get_variation_vector()
half = Scalar("0.5", value=0.5, attr=["Constant"])
L = half * m * l * l * Dot(omega, omega) - m * g * l * Dot(q, e3)
variables = SystemVariables(vectors=[q])

eom = compute_eom(L, Dot(xi, f), variables)
sf = to_standard_form(eom, variables, inputs=[f])
layout = StateLayout(variables, inputs=[f], params=[e3, g, l, m])
src = generate_module(sf, layout, "spherical_pendulum.py", lagrangian=L)
```

`state` is `[q (3), ω (3)]`, `a` is `[ω̇ (3)]`, `u` is `[f (3)]`. The generated functions
(from `examples/generated/spherical_pendulum_s2.py`):

```python
def mass(state, p):
    """Mass matrix M(state): blocks M[a_j] per equation (M = -physical inertia)."""
    _check_params(p, _PARAMS)
    l = p["l"]
    m = p["m"]
    out = np.zeros((3, 3))
    out[0:3, 0:3] = np.eye(3) * (m * l * l * -1)
    return out


def bias(state, p):
    """Bias f(state): the terms of each equation free of accelerations and inputs."""
    _check_params(p, _PARAMS)
    state = _arr(state, 6, "state")
    q = state[0:3]
    e3 = _vec(p["e3"])
    g = p["g"]
    l = p["l"]
    m = p["m"]
    out = np.zeros(3)
    out[0:3] = np.cross(q, e3) * (m * g * l * -1)
    return out


def energy(state, p):
    """Total energy E = sum((dL/dv) . v) - L at the state."""
    ...
    _t0 = m * l * l * (omega_q @ omega_q)
    return _t0 * 0.5 + m * g * l * (q @ e3)
```

Used directly:

```python
mod = load_module("spherical_pendulum.py")   # or: import spherical_pendulum as mod
p = {"e3": [0, 0, 1], "g": 9.81, "l": 1.0, "m": 1.0}
q0 = np.array([np.sin(0.3), 0.0, np.cos(0.3)])
state = np.concatenate([q0, [0.0, 0.5, 0.0]])          # ω ⊥ q
M, f = mod.mass(state, p), mod.bias(state, p)
mod.energy(state, p)                                     # ½ m l² |ω|² + m g l q·e₃
```

## Translation and conventions

Each node class maps to one construct; vectors are `(3,)` arrays and matrices `(3, 3)`
(CasADi: `3×1` columns and `3×3`).

| Node | NumPy | CasADi | Note |
|---|---|---|---|
| `Add`/`VAdd`/`MAdd` | `a + b + …` | same | |
| `Mul`/`SVMul`/`SMMul` | `a * s` | same | scalar times anything |
| `MVMul` | `M @ v` | `ca.mtimes(M, v)` | a row `vᵀM` is generated as `Mᵀv` |
| `MMMul` | `A @ B` | `ca.mtimes(A, B)` | |
| `VVMul` | inner `a @ b`, outer `np.outer(a, b)` | `ca.dot`, `ca.mtimes(a, b.T)` | by `is_row_vector` |
| `Dot` | `a @ b` | `ca.dot(a, b)` | |
| `Cross` | `np.cross(a, b)` | `ca.cross(a, b)` | |
| `Hat` / `Vee` | `_hat(v)` / `_vee(M)` | same | helpers emitted once per module |
| `Transpose` | `M.T`; of a vector, the vector | same | |
| `TimeDerivative` of a leaf | its velocity or acceleration slot | same | anything else raises |
| zero / identity | `0.0`, `np.zeros(3)`, `np.zeros((3, 3))`, `np.eye(3)` | `ca.DM.zeros`, `ca.DM.eye(3)` | |

- **Rows.** There is one vector type; a row (`is_row_vector`) is generated as the same
  3-vector as its column, and the structure only decides between the inner and the outer
  product of a `VVMul` and the side of a row `MVMul`. A scalar equation's block for a
  vector acceleration (or input, or linearization column) fills a 1×3 block of the matrix.
- **Common subexpressions.** A subtree that occurs more than once (keyed by `str`, the
  engine's identity) is computed once into a local `_tN`; long expressions are split into
  locals so that every line is short. Nothing else is simplified.
- **Signs.** Those of [Conventions](conventions.md): `M` is −(physical inertia) in the
  standard form (`mass` of a point mass is `−m I`), and the linearization blocks have the
  positive inertia (`M δẋ = A δx + B δu`).
- **Values.** The CasADi functions take NumPy arrays or CasADi expressions (`SX`, `MX`,
  `DM`) for `state`, `a`, `u` and the values of `p` and return CasADi objects, so
  `ca.Function` and `ca.jacobian` apply to them directly.

### The S² mass block

An S² equation is meaningful only on the tangent plane, `(I − qqᵀ) E = 0`
([Conventions](conventions.md#s2-equations-live-on-the-tangent-plane)), and the module
returns the blocks as derived, not projected. Depending on the Lagrangian the S² mass
block is singular (rank 2: `Hat(q) Hat(q) m_Q l²` in the quadrotor with a suspended load)
or full rank (`−m l² I` for the spherical pendulum, where the engine already used
`q · ω̇ = 0`). Either way the unprojected 3×3 solve is not the equation of motion: the
generated module does **not** solve for the accelerations. Replace each S² block row by
`P_q · row` (with `P_q = I − qqᵀ`) plus `q qᵀ` on its own `ω̇` column and its right-hand
side by `−P_q (bias + G u)`; the system is then square and full rank and its solution
satisfies `q · ω̇ = 0`. `geomech.simulation.accelerations` does exactly this
([Simulation](simulation.md#accelerations-and-the-s2-solve)).

## Total energy

`energy_function(L, layout)` returns `E(state, p) = Σ (∂L/∂v)·v − L` over the velocities
`ṡ`, `ẋ`, `ω`, `Ω`, generated through the same code path (`energy_expr` gives the
expression, `momenta` the `∂L/∂v`). `∂L/∂v` is derived with the engine: the first
variation `δL` is expanded, its terms are split by perturbation symbol as `linearize`
does, and `extract_linear_coeff` gives the coefficient of each velocity perturbation
(`d/dt(δs)`, `d/dt(δx)`, `δω`, `δΩ`). `Σ (∂L/∂v)·v` is the derivative of `L(q, λv)` at
`λ = 1`, which does not depend on the velocity coordinates as long as they are linear in
the rates (`q̇ = ω × q`, `Ṙ = R Ω̂`), so `E = L₂ − L₀` is the energy of Lee et al. The
tests compare `∂L/∂v` with a finite difference of the generated Lagrangian, `E` with a
finite difference in `λ`, and check `dE/dt = 0` along the unforced dynamics of the nine
systems of `tests/vbl_systems.py`.
