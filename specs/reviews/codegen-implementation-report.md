# Code generation (`specs/codegen.md`): implementation report

Branch `feature/codegen-simulation`, on top of 19cab76. Nothing is committed.

## What was built

`geomech/codegen/` (1656 lines over 8 files, dispatch per node class like the rest of the engine):

| File | Contents |
|---|---|
| `layout.py` | `StateLayout`, `Slot`, `collect_params`, `is_literal_leaf`, `base_symbol` |
| `emitter.py` | `Emitter`: tree → statements through a target. It resolves leaves, does CSE by `str`, splits long lines and handles row/column and precedence |
| `numpy_target.py` | `NumpyTarget`: the translation table, helpers (`_arr`, `_vec`, `_mat`, `_hat`, `_vee`, `_check_params`) and block/stack assembly |
| `casadi_target.py` | `CasadiTarget(NumpyTarget)`: the same table with `ca.mtimes`/`ca.dot`/`ca.cross`/`vertcat`/`horzcat`, and `require_casadi` |
| `source.py` | `function_source` (the one code path for every generated function), `to_source`, `lambdify`, `pyliteral`, `get_target` |
| `module.py` | `generate_module`, `module_source`, `load_module` |
| `energy.py` | `momenta`, `energy_expr`, `energy_source`, `energy_function` |

The other deliverables:

- `examples/generated/`: 18 modules, `<name>.py` (standard form plus `energy`) and `<name>_lin.py` (linearization), for the nine systems of `tests/vbl_systems.py`.
- `scripts/regenerate_generated.py`: rewrites the modules; `--check` exits 1 on drift.
- `docs/codegen.md`: layout table, generated module, worked spherical-pendulum example, translation table and conventions, the S² mass block note and the energy derivation. It has an mkdocs nav entry ("Code Generation").
- README pointer (pipeline step 5, install extras, package tree) and a CHANGELOG entry under Unreleased/Added.
- `pyproject.toml`: extras `casadi = ["casadi==3.8.2"]` and `sim = ["scipy"]`, and `scipy` added to `dev`.
- CI (`.github/workflows/tests.yml`):
  - The lint job also lints `scripts/` and `examples/generated/`.
  - The lint job runs `python scripts/regenerate_generated.py && git diff --exit-code examples/generated/`.
  - The 3.12 test job installs `.[casadi]`.
- `scripts/gen_ref_pages.py` was reformatted (one blank line) so that `ruff format --check scripts/` passes. It was the only pre-existing file that failed.

## API

```python
from geomech.codegen import StateLayout, collect_params, to_source, lambdify, generate_module
from geomech.codegen import energy_function, energy_expr, momenta, load_module, module_source

layout = StateLayout(variables, inputs=[f, M], params=[m, J, g, e3])
params = collect_params(exprs, variables, inputs)       # the leaves the layout must list
src = to_source(expr, layout, target="numpy")           # "_PARAMS = (...)" + def f(state, u, p, a=None)
fn = lambdify(expr, layout, target="numpy")             # fn(state, u, p, a=None); fn.source
src = generate_module(system, layout, path=None, target="numpy", lagrangian=None, description="")
E = energy_function(L, layout)                          # E(state, p); E.expr, E.source
mod = load_module(src_or_path)                          # exec a generated module (tests, scripts)
```

`system` can be a `compute_eom` dict, a `to_standard_form` dict or a `LinearizedSystem`. The generated module provides:

| System | Functions |
|---|---|
| Every module | `LAYOUT`, `layout()` (a deep copy), `_PARAMS`, the helpers |
| Standard form | `mass(state, p)`, `bias(state, p)`, `input_matrix(state, p)`, `kinematics(state)`, and `energy(state, p)` when `lagrangian=L` |
| `compute_eom` | `residual(state, a, u, p)` (the stacked `E`) and `kinematics` |
| `LinearizedSystem` | `M_lin`, `A_lin`, `B_lin` of `(state, a, u, p)`, `kinematics`, and `LINEARIZATION` (the slice of each row, dstate, state and input) |

Every generated function behaves the same way:

- It checks `p` first. A missing key raises `KeyError("missing parameter 'm' ...")` and an extra key raises `ValueError("unknown parameter 'k' ...")`.
- It checks the sizes of `state`, `a` and `u` (`ValueError("state must have 6 entries, got 5")`).
- It returns `float`, `(3,)` or `(3, 3)` for a single expression. Block functions return `(n_rows, n_cols)`.

## Layout decisions the spec left open

- **Variable order** is the `SystemVariables` group order (scalars, vectors, matrices), the order `linearize` and `numeric` already use. Equations (block rows) follow the same order, keyed by variation vector.
- **`a` is the last, optional argument of `lambdify`** (`f(state, u, p, a=None)`). The spec gives `callable(state, u, p)`, but `compute_eom` equations and linearization blocks contain accelerations. The `compute_eom` and linearization module functions take `(state, a, u, p)`, as the spec states for `M_lin`.
- **Parameters**:
  - A parameter can be any `Scalar`/`Vector`/`Matrix` leaf that is not a configuration variable or input, including a non-constant leaf such as the frozen `R` of `adv_virtual_work_under_rotation`. It can also be a `TimeDerivative` chain of a leaf (a prescribed rate `l̇`, keyed `"\frac{d}{dt}(l)"`).
  - A named scalar with a value (`Scalar("m", value=2)`) is a parameter, following decision D1: the name is the symbol and the value only matters for `evaluate`.
  - Literals and zero/identity constants are rejected as parameters with a message saying they are inlined.
  - Vectors may be passed as any 3-entry array and matrices as any 9-entry array.
- **Rows as 3-vectors**: rows and columns are generated as the same 3-vector. This works because the engine never relies on orientation outside `VVMul` (inner or outer product) and the row `MVMul` (`vᵀM` → `M.T @ v`). A scalar row with a vector column fills a 1×3 slice (CasADi transposes it to `1×3`).
- **Identifiers**: backslashes and non-identifier characters are stripped (`\omega_{q}` → `omega_q`, `\rho_1` → `rho_1`). Velocities and accelerations are named `x_dot`/`x_ddot` for R and R³, and `omega_q`/`omega_q_dot` or `Omega_R`/`Omega_R_dot` for S² and SO(3). Python keywords and the argument names (`state`, `a`, `u`, `p`, `np`, `ca`, `out`) get a trailing `_`, and collisions get `_2`, `_3`.
- **Slot kinds per function**: `mass`/`bias`/`input_matrix`/`energy` may only use configuration, velocity and parameter symbols. An acceleration or input in a standard-form block, or an input in `L`, raises `ValueError` at generation time. `kinematics` uses only state symbols.
- **CSE and line length** for both targets:
  - A subtree that occurs twice (by `(class, str)`) goes into a local `_tN`.
  - Operands longer than 30 characters and expressions longer than 72 also go into locals.
  - Long sums are accumulated in pieces.
  - Parentheses follow `*`/`@`/`+` precedence.

  Every line is under 99 characters, and `ruff format` leaves the generated files unchanged. `LAYOUT` is written as fully expanded collections with magic trailing commas, which is also ruff-stable. Output is deterministic: it was identical under `PYTHONHASHSEED=123` and `999`.
- **CasADi** functions accept NumPy arrays or `SX`/`MX`/`DM` (state, `a`, `u` and the values of `p`) and return CasADi objects. `R` is unpacked as `ca.reshape(state[i:i+9], 3, 3).T` to keep the row-major convention. CasADi gets the same CSE as NumPy; the spec said "CasADi does its own", and the extra CSE is harmless and keeps the source short.
- **`generate_module(path=None)`** returns the source without writing. `module_source` is the same function without a path.

## Energy: derivation and verification

`E = Σ (∂L/∂v)·v − L`, with `v ∈ {ṡ, ẋ, ω, Ω}`.

`∂L/∂v` comes from the engine, not from a separate differentiator:

1. Take `expand(full_simplify(L.delta()))`.
2. Partition its terms by perturbation symbol with `linearize`'s own `_partition_terms`, using the velocity perturbations (`d/dt(δs)`, `d/dt(δx)`, `δω`, `δΩ`, from `_perturbation_variables`) before the configuration ones. Every term must hold exactly one perturbation symbol.
3. `extract_linear_coeff` gives the coefficient of each velocity perturbation.
4. `full_simplify` the coefficient.

`Σ (∂L/∂v)·v` is `d/dλ L(q, λv)|₁`, which is invariant under any velocity coordinates that are linear in the rates (`q̇ = ω × q` with `ω ⊥ q`, `Ṙ = R Ω̂`). E is therefore `L₂ − L₀`, the energy of Lee et al. For S², `∂L/∂ω` may have a component along `q`. That component drops out of `(∂L/∂ω)·ω` because `ω ⊥ q`.

Verification (`tests/test_codegen_energy.py`):

1. `momenta` against a central difference (`h = 1e-5`) of the generated Lagrangian in every velocity entry, for all 64 systems at 5 states each. S² is compared after projecting with `I − qqᵀ`. Tolerance is `1e-7 · max(1, |∂L/∂v|)`.
2. `E` against `d/dλ L(λv) − L` by finite difference, for all 64 systems.
3. Closed forms:
   - point mass: `str(E) == str(m(ẋ·ẋ)·0.5)`;
   - spherical pendulum: `½ml²|ω|² + mgl q·e₃` (to 1e-12);
   - free rigid body: `½ΩᵀJΩ` (to 1e-12).
4. `dE/dt = 0` along the unforced dynamics of the nine VBL systems. Accelerations come from the generated `mass`/`bias`, with S² rows solved on the tangent plane as in `specs/simulation.md`. Masses are positive, inertias positive definite, and the check uses a central difference along `[kinematics, a]` (≤ 1e-7 relative).
5. Spec item 5: a free rigid body integrated for 10 s with `solve_ivp` (DOP853, rtol = atol = 1e-12) from the generated `mass`/`bias`/`kinematics` keeps the generated energy to **1.2e-13 relative** (the gate is 1e-8). The run takes 1046 RHS evaluations and 0.027 s.

## Tests

There are 484 new tests in `tests/test_codegen*.py`, all passing. The full suite is 2881 passed and 2 xfailed under `-W error` in about 17 s; it was 2397 before.

| File | Tests | What it checks |
|---|---|---|
| `test_codegen_numeric.py` | 194 | Spec item 1, for all **64 systems** (the 42 oracle-gate builders: 9 VBL, 12 pinned elsewhere, 15 phase-1 S²/collect and 6 row systems; plus the 22 phase-1 adversarial ones): every `M`/`f`/`G` block, every `compute_eom` equation, the kinematics, and every `M_lin`/`A_lin`/`B_lin` block (63 systems; `adv_virtual_work_under_rotation` does not linearize) against `Evaluator` at 20 `sample_state` points, `|gen − ref| ≤ 1e-12·max(1, |ref|)`. Also `lambdify` checks. |
| `test_codegen_casadi.py` | 65 | Spec item 2: every function of the three modules of the 64 systems, CasADi against NumPy at the same 20 states (1e-12). SX arguments, `ca.Function` and `ca.jacobian` of `bias` against a finite difference. Skipped with an install message when casadi is absent. |
| `test_codegen_translation.py` | 55 | Every row of the translation table on hand-built trees, both targets: row `MVMul` with a non-symmetric matrix, inner and outer `VVMul`, `Vee`, `Hat`, transposes, zeros, identity, number literals, velocities, accelerations, inputs. Also CSE, the line-length bound and result shapes. |
| `test_codegen_layout.py` | 25 | Slot order and sizes, identifiers, R row-major, `to_dict`, `collect_params`. Spec item 4 errors: missing parameter, extra parameter, stray leaf (expression and block), `TimeDerivative` of a non-leaf, of an unknown leaf, a `Variation`, a kind restriction, wrong array sizes, and the `StateLayout` validation messages. |
| `test_codegen_energy.py` | 139 | The energy checks above. |
| `test_codegen_generated.py` | 6 | Spec item 3: regeneration produces no diff. All 18 files run in a `python -I` subprocess with `geomech` blocked by a meta-path finder, and their values match the in-process modules. `ruff format --check` and `ruff check` pass. `import geomech, geomech.codegen` loads neither scipy nor casadi. The CasADi target with casadi blocked raises `ImportError` mentioning `geomech[casadi]`. |

The translation tests exist because of a mutation check. Dropping the transpose of a row `MVMul` went unnoticed by the 64-system comparison: no derived block contains a row `MVMul`, an inner `VVMul` or a `Vee`, and every derived block contains only outer products, `Hat` and transposed columns. The translation tests catch that mutation in 5 tests.

Other checks run:

- `ruff format --check` and `ruff check` on `geomech/ tests/ scripts/ examples/generated/` pass.
- `mypy geomech --ignore-missing-imports` passes (42 files).
- `mkdocs build --strict` passes; griffe asked for annotations on two signatures, which were added.
- `python -c "import geomech"` does not import scipy or casadi.
- The codegen tests with casadi shadowed by a module raising `ModuleNotFoundError`: 394 passed, 26 skipped.

## Timing (nine VBL systems, this machine, NumPy 2.5, Python 3.12)

| System | Derive (s) | Generate sf + lin (ms) | Lines | `mass`+`bias`+`G`+`kinematics` (µs/call) | `M_lin`+`A_lin`+`B_lin` (µs/call) |
|---|---|---|---|---|---|
| point_mass_r3 | 0.001 | 0.8 | 327 | 4.0 | 5.2 |
| spring_mass_r | 0.002 | 1.2 | 347 | 1.9 | 1.4 |
| free_rigid_body | 0.003 | 1.1 | 341 | 17.7 | 10.0 |
| rigid_pendulum_pivot | 0.004 | 1.5 | 407 | 30.4 | 15.5 |
| spherical_pendulum_s2 | 0.002 | 1.4 | 389 | 25.2 | 14.0 |
| quadrotor | 0.005 | 1.9 | 500 | 22.9 | 18.9 |
| double_rigid_pendulum | 0.046 | 9.0 | 728 | 319.7 | 275.6 |
| elastic_spherical_pendulum | 0.008 | 3.4 | 530 | 28.8 | 31.5 |
| quadrotor_with_load | 0.023 | 6.1 | 692 | 121.7 | 120.9 |

"Derive" is `compute_eom` + `to_standard_form` + `linearize` + layout. "Lines" counts both modules. Regenerating all 18 files takes about 0.25 s. CasADi source generation takes 0.2–0.9 ms per module. Call cost is Python overhead on 3-vectors. The double pendulum is the outlier: its unsimplified cross-coupling terms need many small numpy calls.

## Notes and things not done

- **The S² mass block is not always rank 2.** The spec says it is "singular (rank 2) by construction". That holds for the quadrotor with load (`Hat(q)Hat(q) m_Q l²`), but for the spherical pendulum the engine already simplified with `q·ω̇ = 0` and returns `−m l² I`, which is full rank. `docs/codegen.md` says so. In both cases the tangent-plane solve of `specs/simulation.md` is needed and works: `tests/codegen_systems.tangent_accelerations` implements it for the `dE/dt` test, and the simulation feature can reuse it.
- **`energy` uses private `linearize` helpers** (`geomech.dynamics.vbl._perturbation_variables`, `_partition_terms`). They are imported, not copied, so the energy derivation uses exactly the term partition of the linearization. Making them public is a small follow-up.
- **The spec's `layout()`** is provided as a function next to the `LAYOUT` constant.
- **Prescribed parameters** (`compute_eom(..., params=[l])`) are supported through `params=[l, l.t_diff()]`. `energy_function` does not support them, because their variation is not a perturbation symbol and the partition raises `AlgebraicError`.
- **The CasADi CI job** is the 3.12 job of the test matrix. `casadi==3.8.2` is pinned exactly; it is the current PyPI release and was installed in the venv for these tests.
- **Not done** (non-goals of the spec): C/C++ or JAX output, simplification beyond CSE, solving for accelerations, integrators.
