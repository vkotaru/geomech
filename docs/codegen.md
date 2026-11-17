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
| `s: Scalar` ($\mathbb{R}$) | $s$ (1) | $\dot s$ (1) | $\ddot s$ |
| `x: Vector` ($\mathbb{R}^3$) | $x$ (3) | $\dot x$ (3) | $\ddot x$ |
| `q: S1` | $q$ (3, unit, $\perp e$) | $\omega_q$ (1) | $\dot{\omega_q}$ (1) |
| `q: S2` | $q$ (3, unit) | $\omega_q$ (3, $\perp q$) | $\dot{\omega_q}$ |
| `R: SO3` | $R$ (9, row-major) | $\Omega_R$ (3) | $\dot{\Omega_R}$ |

- `state` holds every configuration slot in variable order (the `SystemVariables` order:
  scalars, vectors, matrices), then every velocity slot in the same order.
  `layout.slices[str(symbol)]` is the slice of a configuration or velocity symbol;
  `layout.accel_slices` and `layout.input_slices` those of `a` and `u`.
- `a` stacks the accelerations in variable order; it is an argument only of the functions
  whose expressions contain an acceleration (`compute_eom` equations, linearization blocks).
- `u` stacks the inputs in the order given (a scalar takes 1 entry, a vector 3).
- **Required arguments.** A generated function converts (and so requires) only the arrays
  its expressions use: `a` when they contain an acceleration, `u` when they contain an
  input. Passing `None` for a required array raises `ValueError`
  (`u must have 1 entries, got None (u is required here)`), and so does a wrong size.
  A function whose expressions have no input accepts `u=None`, which is always the case
  for a system without inputs; likewise `a=None` (the default of `lambdify`) without
  accelerations. The generated code never reads `None` as zero input. That shorthand
  belongs to `geomech.simulation` (`accelerations(module, state, u=None, p)`,
  `simulate(..., u=None)`), whose functions document `None` as zero input; `mass`,
  `bias`, `input_matrix` and `energy` take no `u` at all.
- `p` is a `dict` keyed by `str(symbol)` of every parameter: for a scalar a Python or
  NumPy scalar, or a 0-d or single-element array or sequence (converted to a float; more
  entries raise `ValueError` naming the parameter); any array of 3 or $3 \times 3$ entries for a
  vector or matrix. A missing key raises `KeyError`, an extra one `ValueError`, both
  naming the key. Numeric literals (`0.5`, `(-1)`) and the
  zero and identity constants are inlined, never parameters. A *named* scalar with a value
  (`Scalar("m", value=2)`) is a symbol, so it is a parameter.
- `collect_params(system, variables, inputs)` lists the symbols that the variables and
  inputs do not cover (sorted by name): the `params` to pass. Pass what the system is made
  of: `[L, dW]` (the usual call), `[L]` plus the equations, the dict of `compute_eom` or
  of `to_standard_form`, a `LinearizedSystem`, one expression, or a list of these. The
  tangent and variation vectors of the variables ($\omega_q$, $\xi_q$, $\Omega_R$, $\eta_R$, which `dW` and
  `compute_eom` hold) are left out, and `StateLayout` rejects them as parameters.
  Anything it cannot walk (a number, a string, an object without symbols) raises
  `TypeError` instead of giving `[]`.
- `layout.pack(**slots)` builds a flat `state` from named slot values and
  `layout.unpack(state)` splits one (or a `(N, n_state)` array such as
  `Trajectory.state`) into a dict. The names are the slot identifiers below (`x_L`, `q`,
  `R`, `x_L_dot`, `omega_q`, `Omega_R`); `pack` also takes a dict keyed by those names or
  by the printed symbols. Every configuration slot is required, a velocity left out is
  zero, `R` is a $3 \times 3$ matrix or its 9 row-major entries. `pack` does not project onto the
  manifold; `simulate` and `lqr` check the state. Velocity names follow the variable's
  name, indices included: `x_dot`, `omega_q1` for `S2("q1")`, `Omega_R2` for
  `SO3("R2")`; `layout.slot_names()` lists every key in state order.

  ```python
  x_eq = layout.pack(x_L=np.zeros(3), q=[0, 0, -1], R=np.eye(3))   # at rest
  layout.unpack(traj.state)["q"]                                    # (N, 3)
  layout.slot_names()          # ['x_L', 'q', 'R', 'x_L_dot', 'omega_q', 'Omega_R']
  ```
- `layout.state_from_env(env)`, `accel_from_env`, `inputs_from_env`, `params_from_env`
  convert an `Evaluator` env (printed-name keys, `(3, 1)` columns) to the flat arrays, and
  `env_from_arrays(state, a, u, p)` converts back.
- Each slot also has a Python identifier for the generated code (`\omega_{q}` → `omega_q`,
  $\dot{\omega_q}$ → `omega_q_dot`, $\dot x$ → `x_dot`, $\ddot x$ → `x_ddot`).

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
def kinematics(state)      -> (n_config,)   # ṡ, ẋ, ω e × q, ω × q, R hat(Ω) row-major
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

**One layout for everything.** Generate the standard-form and the linearization module
from the same `StateLayout`, built with the same `SystemVariables` and the same `inputs`
list (same order) as `to_standard_form` and `linearize`. The columns of `input_matrix`
follow the layout's inputs and those of `B_lin` follow `lin.inputs`; with a different
order the two would pair a column with a different entry of `u`, and a controller
designed on `B_lin` would act on the wrong channels. `generate_module` therefore raises
`ValueError` when the variables (rows and states) or inputs of a `LinearizedSystem` are
not those of the layout in the same order, naming both orders, and
`lqr(..., module=...)` checks that the two modules have the same `LAYOUT` (see
[Troubleshooting](troubleshooting.md#code-generation-and-simulation)).

`examples/generated/` holds the modules of the nine systems of `tests/vbl_systems.py`
(`<name>.py` from the standard form with `energy`, `<name>_lin.py` from `linearize`).
`python scripts/regenerate_generated.py` rewrites them (and deletes the files of systems
no longer listed); `--check` writes nothing and exits 1 when a module is out of date or
missing, when a stale file is left over, or when a generated file is not tracked by git.
CI runs `--check`.

## Worked example: the spherical pendulum

Lee et al. eq. 5.14, $L = \tfrac12 m l^2\,\omega\cdot\omega - m g l\, q\cdot e_3$ and
$\delta W = \xi\cdot M$, with $M$ a moment about the pivot:

```python
import numpy as np
from geomech import S2, Dot, Scalar, SystemVariables, Vector, compute_eom, getScalars
from geomech import to_standard_form
from geomech.codegen import StateLayout, generate_module, load_module

m, g, l = getScalars("m g l", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
q, M = S2("q"), Vector("M")
omega, xi = q.get_tangent_vector(), q.get_variation_vector()
half = Scalar("0.5", value=0.5, attr=["Constant"])
L = half * m * l * l * Dot(omega, omega) - m * g * l * Dot(q, e3)
variables = SystemVariables(vectors=[q])

eom = compute_eom(L, Dot(xi, M), variables)
sf = to_standard_form(eom, variables, inputs=[M])
layout = StateLayout(variables, inputs=[M], params=[e3, g, l, m])
src = generate_module(sf, layout, "spherical_pendulum.py", lagrangian=L)
```

`state` is $[q, \omega_q]$ (3 + 3 entries), `a` is $[\dot{\omega_q}]$ (3), `u` is $[M]$ (3). The generated functions
(from `examples/generated/spherical_pendulum_s2.py`):

```python
def mass(state, p):
    """Mass matrix M(state): blocks M[a_j] per equation (M = -physical inertia)."""
    _check_params(p, _PARAMS)
    l = _scalar(p["l"], "l")
    m = _scalar(p["m"], "m")
    out = np.zeros((3, 3))
    out[0:3, 0:3] = np.eye(3) * (m * l * l * -1)
    return out


def bias(state, p):
    """Bias f(state): the terms of each equation free of accelerations and inputs."""
    _check_params(p, _PARAMS)
    state = _arr(state, 6, "state")
    q = state[0:3]
    e3 = _vec(p["e3"])
    g = _scalar(p["g"], "g")
    l = _scalar(p["l"], "l")
    m = _scalar(p["m"], "m")
    out = np.zeros(3)
    out[0:3] = _cross(q, e3) * (m * g * l * -1)
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
M_a, f = mod.mass(state, p), mod.bias(state, p)       # mass matrix and bias
mod.energy(state, p)                                     # ½ m l² |ω|² + m g l q·e₃
```

## Translation and conventions

Each node class maps to one construct; vectors are `(3,)` arrays and matrices `(3, 3)`
(CasADi: $3 \times 1$ columns and $3 \times 3$).

| Node | NumPy | CasADi | Note |
|---|---|---|---|
| `Add`/`VAdd`/`MAdd` | `a + b + …` | same | |
| `Mul`/`SVMul`/`SMMul` | `a * s` | same | scalar times anything |
| `MVMul` | `M @ v` | `ca.mtimes(M, v)` | a row $v^T M$ is generated as $M^T v$ |
| `MMMul` | `A @ B` | `ca.mtimes(A, B)` | |
| `VVMul` | inner `a @ b`, outer `np.outer(a, b)` | `ca.dot`, `ca.mtimes(a, b.T)` | by `is_row_vector` |
| `Dot` | `a @ b` | `ca.dot(a, b)` | |
| `Cross` | `_cross(a, b)` | `ca.cross(a, b)` | NumPy: component-wise helper, same values as `np.cross` at a fraction of its per-call cost |
| `Hat` / `Vee` | `_hat(v)` / `_vee(M)` | same | helpers emitted once per module |
| `Transpose` | `M.T`; of a vector, the vector | same | |
| `TimeDerivative` of a leaf | its velocity or acceleration slot | same | anything else raises |
| zero / identity | `0.0`, `np.zeros(3)`, `np.zeros((3, 3))`, `np.eye(3)` | `ca.DM.zeros`, `ca.DM.eye(3)` | |
| `±inf` / `nan` literal | `np.inf`, `-np.inf`, `np.nan` | `ca.inf`, `-ca.inf`, `np.nan` | inlined, never a parameter |

- **Rows.** There is one vector type; a row (`is_row_vector`) is generated as the same
  3-vector as its column, and the structure only decides between the inner and the outer
  product of a `VVMul` and the side of a row `MVMul`. A scalar equation's block for a
  vector acceleration (or input, or linearization column) fills a $1 \times 3$ block of the matrix.
- **Common subexpressions.** A subtree that occurs more than once is computed once into a
  local `_tN`; long expressions are split into locals so that every line is short.
  Nothing else is simplified. Subtrees are compared structurally (a leaf by its `str`,
  the engine's identity; an operation by its class and operands), with the identity of
  each node object computed once per function, so generation time is linear in the
  number of node objects, also for trees that reuse a subtree many times.
- **Signs.** Those of [Conventions](conventions.md): $M$ is $-(\text{physical inertia})$ in the
  standard form (`mass` of a point mass is $-m I$), and the linearization blocks have the
  positive inertia ($M\, \delta\dot x = A\, \delta x + B\, \delta u$).
- **Values.** The CasADi functions take NumPy arrays or CasADi expressions (`SX`, `MX`,
  `DM`) for `state`, `a`, `u` and the values of `p` and return CasADi objects, so
  `ca.Function` and `ca.jacobian` apply to them directly.

### The S² mass block

An $S^2$ equation is meaningful only on the tangent plane, $(I - q q^T) E = 0$
([Conventions](conventions.md#s2-equations-live-on-the-tangent-plane)), and the module
returns the blocks as derived, not projected. Depending on the Lagrangian the $S^2$ mass
block is singular (rank 2: $\widehat{q}\,\widehat{q}\, m_Q l^2$ in the quadrotor with a suspended load)
or full rank ($-m l^2 I$ for the spherical pendulum, where the engine already used
$q \cdot \dot{\omega_q} = 0$). Either way the unprojected $3 \times 3$ solve is not the equation of motion: the
generated module does **not** solve for the accelerations. Replace each $S^2$ block row by
$P_q \cdot \text{row}$ (with $P_q = I - q q^T$) plus a multiple $s\, q q^T$ on its own $\dot{\omega_q}$ column ($s$ =
the tangent inertia, see [Simulation](simulation.md#accelerations-and-the-s2-solve)) and
its right-hand
side by $-P_q (\text{bias} + G u)$; the system is then square and full rank and its solution
satisfies $q \cdot \dot{\omega_q} = 0$. `geomech.simulation.accelerations` does exactly this
([Simulation](simulation.md#accelerations-and-the-s2-solve)).

## Total energy

`energy_function(L, layout)` returns `E(state, p)` $= \sum (\partial L/\partial v) \cdot v - L$ over the velocities
$\dot s$, $\dot x$, $\omega_q$, $\Omega_R$, generated through the same code path (`energy_expr` gives the
expression, `momenta` the $\partial L/\partial v$). $\partial L/\partial v$ is derived with the engine: the first
variation $\delta L$ is expanded, its terms are split by perturbation symbol as `linearize`
does, and `extract_linear_coeff` gives the coefficient of each velocity perturbation
($\tfrac{d}{dt}\delta s$, $\tfrac{d}{dt}\delta x$, $\delta\omega_q$, $\delta\Omega_R$). $\sum (\partial L/\partial v) \cdot v$ is the derivative of $L(q, \lambda v)$ at
$\lambda = 1$, which does not depend on the velocity coordinates as long as they are linear in
the rates ($\dot q = \omega_q \times q$, $\dot R = R \widehat{\Omega_R}$), so $E = L_2 - L_0$ is the energy of Lee et al. The
tests compare $\partial L/\partial v$ with a finite difference of the generated Lagrangian, $E$ with a
finite difference in $\lambda$, and check $dE/dt = 0$ along the unforced dynamics of the nine
systems of `tests/vbl_systems.py`.
