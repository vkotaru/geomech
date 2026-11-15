# Spec: code generation (`geomech.codegen`)

Roadmap feature 5. Branch `feature/codegen` from `geomech`.

## Goal

Turn the symbolic output of the pipeline into plain numeric functions that simulation,
control design and other tools can call, through one
code path that is verified against `geomech.numeric.evaluate`.

Two targets:

| Target | Dependency | Use |
|---|---|---|
| NumPy | none beyond the runtime `numpy` | default; simulation with `scipy.integrate.solve_ivp`, tests, examples |
| CasADi | optional extra `geomech[casadi]` | optimal control, automatic derivatives; the generated code imports `casadi` only when this target is chosen |

## Public API

```python
from geomech.codegen import to_source, lambdify, generate_module, StateLayout

layout = StateLayout(variables, inputs=[f, M], params=[m, J, g, e3])
src: str = to_source(expr, layout, target="numpy")          # one expression → Python source
fn = lambdify(expr, layout, target="numpy")                  # → callable(state, u, p)
generate_module(system, layout, path, target="numpy")        # writes a .py file
```

`system` is any of: the `dict` returned by `compute_eom`, the `dict[str,
StandardFormEquation]` from `to_standard_form`, or a `LinearizedSystem`. The module
contains one function per block, plus `kinematics` and `layout()`; see "Generated module".

### `StateLayout`

The single source of truth for how symbols map to arrays. It is built from
`SystemVariables`, `inputs` and `params`, in that order of variables as listed, and is
stored in the generated module so that callers never guess.

| Variable | Configuration slot | Velocity slot | Acceleration symbol |
|---|---|---|---|
| `s: Scalar` (R) | `s` (1) | `ṡ` (1) | `s̈` |
| `x: Vector` (R³) | `x` (3) | `ẋ` (3) | `ẍ` |
| `q: S2` | `q` (3, unit) | `ω` (3, ⊥ q) | `ω̇` |
| `R: SO3` | `R` (9, row-major) | `Ω` (3) | `Ω̇` |

- `state` is a flat `np.ndarray` = all configuration slots in variable order, then all
  velocity slots in the same order. `layout.slices` gives the slice of each symbol.
- `u` is a flat array of the inputs in the given order (scalars 1, vectors 3).
- `p` is a `dict[str, float | np.ndarray]` keyed by `str(symbol)` for every parameter
  (named constant). A missing or extra key raises `KeyError`/`ValueError` naming it.
  Numeric literals are inlined, never parameters.
- Every leaf of the expression must be a configuration, velocity, acceleration (only for
  `compute_eom` output), input, or parameter of the layout; anything else raises
  `ValueError` at generation time naming the leaf. This is the codegen analogue of the
  stray-symbol check.

## Translation rules (NumPy target)

Each node class maps to one NumPy construct; vectors are shape `(3,)`, matrices `(3, 3)`,
rows are produced only inside inner products.

| Node | NumPy | Note |
|---|---|---|
| `Add`/`VAdd`/`MAdd` | `a + b + …` | |
| `Mul`/`SVMul`/`SMMul` | `a * s` | scalar times anything |
| `MVMul` | `M @ v` | a row operand `vᵀ M` becomes `(Mᵀ v)` |
| `MMMul` | `A @ B` | |
| `VVMul` | inner `a @ b`, outer `np.outer(a, b)` | decided by `is_row_vector` |
| `Dot` | `a @ b` | |
| `Cross` | `_cross(a, b)` | component-wise helper (fix F5 of `specs/reviews/codegen-review.md`) |
| `Hat` / `Vee` | `_hat(v)` / `_vee(M)` | helpers emitted once per module |
| `Transpose` | `M.T`; of a vector, the column (row form handled by the parent) | |
| `TimeDerivative` of a leaf | the corresponding velocity or acceleration slot | anything else raises |
| zero / identity constants | `0.0`, `np.zeros(3)`, `np.zeros((3,3))`, `np.eye(3)` | typed |

The CasADi target uses the same table with `ca.mtimes`, `ca.cross`, `ca.dot`, `ca.SX`
helpers. Common subexpressions: NumPy target emits one local variable per repeated
subtree (keyed by `str`, which is the engine's identity); CasADi does its own.

## Generated module

For a standard-form system with configuration variables v₁…vₙ:

```python
LAYOUT = {...}                         # StateLayout as plain data
def mass(state, p) -> np.ndarray       # block matrix of M[aⱼ] per row, assembled
def bias(state, p) -> np.ndarray       # stacked f
def input_matrix(state, p) -> np.ndarray   # block matrix of G[uₖ]
def kinematics(state) -> np.ndarray    # configuration rates: ṡ, ẋ, ω×q, R·hat(Ω) (flattened)
```

so that the equations read `mass @ a + bias + input_matrix @ u = 0` with `a` the stacked
accelerations. For a `LinearizedSystem`: `M_lin(state, a, u, p)`, `A_lin(...)`,
`B_lin(...)` with the row/column order of `LinearizedSystem.rows/states/inputs`.

Sign and shape conventions are those of `docs/conventions.md` (`M` = −physical inertia;
a scalar row with a vector acceleration contributes a 1×3 block). The S² mass block is
singular (rank 2) by construction (correction: only for some systems, e.g. the quadrotor
with a load; the spherical pendulum's is `−m l² I`, full rank); the module does **not**
solve for accelerations.
Solving with the tangency constraint q·ω̇ = 0 belongs to the simulation feature.

## Verification (definition of done)

1. For every system in `tests/vbl_systems.py`, `tests/phase1_systems.py` and the numerical
   check systems: each generated NumPy block equals `geomech.numeric.evaluate` of the symbolic
   block at 20 random admissible states (`sample_state`) to 1e-12 relative.
2. CasADi target equals NumPy target on the same states (skipped if `casadi` is absent;
   CI installs it in one job).
3. `generate_module` output for the nine VBL systems is checked into
   `examples/generated/`, ruff-clean, importable, and regenerated in CI; any drift fails
   the build.
4. Layout errors: missing parameter, extra parameter, stray leaf, `TimeDerivative` of a
   non-leaf, each raises with the symbol named.
5. A free rigid body simulated for 10 s with `scipy.integrate.solve_ivp` from the
   generated `mass`/`bias`/`kinematics` conserves energy to 1e-8 (smoke test only; full
   simulation is the next feature).
6. Docs page `docs/codegen.md` with the layout table, one worked example and the
   conventions; CHANGELOG entry; mypy and ruff clean; no new runtime dependency.

## Non-goals

C/C++ or JAX output, symbolic simplification of generated code beyond CSE, solving
for accelerations, integrators, controllers.
