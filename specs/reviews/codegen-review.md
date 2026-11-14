# Code generation adversarial review

Branch `feature/codegen-simulation`, range `19cab76..cc1f58f` ("Add code generation
(geomech.codegen)"). The later simulation commit changes only the `generate_module`
docstring under `geomech/codegen/`, so the checkout was reviewed in place. Repo venv
(Python 3.12.13, numpy 2.5.4, casadi 3.8.2). The CI drift checks ran in a detached worktree
of cc1f58f, which was removed afterwards. New tests are in `tests/test_codegen_adversarial.py`:
53 pass and 6 are strict xfail. The full suite passes under `-W error`, and ruff is clean.
Nothing was committed, and no existing file was changed.

## Verdict: GO WITH CHANGES

No wrong number was found on any input the contract allows. 12 000 fuzzed trees agree with
`numeric.evaluate` on both targets. The layout, the kinematics, the energy and the CI
regeneration all hold. The changes needed before merge are in the argument helpers of the
generated code:

- **F1:** an omitted `a` or `u` of size 1 silently becomes NaN.
- **F2:** a scalar parameter passed as a 1-element list is repeated as a list instead of being multiplied.

Both are short fixes in `_arr` and `unpack`. F3 to F6 are small and can follow.

## Findings, ranked

| # | rank | finding | where | reproducer |
|---|---|---|---|---|
| F1 | wrong numbers (silent NaN) | `_arr(None, n, name)` runs `np.asarray(None, dtype=float)`, which is `array(nan)` of size 1. `lambdify` functions default to `a=None`, so `lambdify(s̈ + f, layout)(state, u, p)` on a system with one scalar variable returns `nan` instead of raising. The same happens when a module's `residual`/`M_lin`/`A_lin`/`B_lin` gets `u=None` (or `a=None`) and the input (or acceleration) has size 1: `spring_mass_r` gives `residual(state, [3.0], None, p) == [nan]`. When the size is larger, the error is misleading: `"u must have 6 entries, got 1"`. `geomech.simulation.accelerations` accepts `u=None` as "zero input", so users will expect the generated functions to do the same. **Expected:** a `TypeError`/`ValueError` saying the argument is required. The CasADi `_arr` behaves the same way (`ca.DM(nan)`) | `geomech/codegen/numpy_target.py:28-32`, `casadi_target.py:179-184`, `source.py:111` (`a=None`) | `test_omitted_acceleration_raises` |
| F2 | wrong numbers (input outside the documented types) | Scalar parameters are bound raw (`m = p["m"]`), while vectors and matrices go through `_vec`/`_mat`. A scalar given as `[1.5]` is list-multiplied: `lambdify(m*3 + s)(…, {"m": [1.5]})` returns `[2.5, 2.5, 2.5]` (expected `5.5`), and `m*2` gives `[1.5, 1.5]`. A `(1,)` array gives a `(1,)` result from `lambdify`/`energy`, and `ValueError: setting an array element with a sequence` from `mass`/`bias`. Vectors accept any 3-entry sequence, so this asymmetry is easy to hit, for example by slicing a parameter vector from an optimizer. **Fix:** add a `_scalar(x)` helper: `float(np.asarray(x).reshape(()))` for NumPy, and pass through 1×1 CasADi objects | `numpy_target.py:91-94` | `test_scalar_parameter_as_one_element_sequence` |
| F3 | invalid code | `Number(float("inf"))` is emitted as `inf`, which raises `NameError` at call time; generation succeeds. `Number(nan)` fails `is_number_literal` because `nan != nan`, so `collect_params` makes it a *parameter* named `"(nan)"`. **Fix:** emit `float("inf")`/`np.inf`, or reject non-finite literals at generation with the node named | `numpy_target.py:62-67` (`_number`), `layout.py:94` | `test_infinite_literal_generates_valid_code` |
| F4 | tooling | `scripts/regenerate_generated.py` only rewrites and compares the files of `SYSTEMS`. When a system is removed from the list, its stale `.py` stays in place, and `--check` exits 0. The CI lint step (`regenerate && git diff --exit-code`) is also blind to *untracked* output: adding a system without committing its files passes that step. Both cases were simulated in the worktree. The test job still catches both through `test_regeneration_produces_no_diff`, which compares the file set, so CI as a whole is safe. The script's own `--check` and the lint step are not. **Fix:** in `main`, report (and in write mode delete) `*.py` files in `OUT` that are not in `sources()`, and use `git status --porcelain examples/generated/` in CI | `scripts/regenerate_generated.py:82-98`, `.github/workflows/tests.yml` (lint job) | `test_check_mode_reports_stale_files` |
| F5 | performance | CSE is effective: no function of the 18 checked-in modules repeats a call or binary-op subexpression (AST scan). The 320 µs of the double rigid pendulum come from `np.cross`, about 9.4 µs per 3-vector call in numpy 2.5. `bias` makes 28 of them (283 µs). A component-wise `_cross(a, b)` helper gives identical values in 47 µs (6×), and cuts `energy` from 91 to 24 µs. `mass` (16 µs) and `kinematics` (5 µs) are already cheap | `numpy_target.py:115-116` | `test_generated_numpy_code_avoids_np_cross` |
| F6 | cosmetic / packaging | The optional extra pins `casadi==3.8.2` exactly. Wheels exist for 3.10 (cp310) and 3.11+ (cp311-abi3), so CI installs it. But an exact pin in a library extra forces a downgrade or a resolver conflict for users who already have casadi. **Fix:** `casadi>=3.6,<4` (only `SX/MX/DM`, `mtimes`, `dot`, `cross`, `vertcat`/`horzcat`, `reshape` are used) | `pyproject.toml:28-30` | `test_casadi_extra_is_not_an_exact_pin` |

Advisory (no test):

- **Generation cost is quadratic in tree depth and exponential on shared DAGs.** The
  emitter computes `str(node)` (uncached, O(subtree)) twice per node (`emitter.py:116,146`).
  A chain of depth 800 (2401 nodes) takes 0.77 s, with the time growing about 4.5× per
  doubling. A flat sum of 800 terms takes 5 ms. A tree that reuses each level twice
  (`VAdd(Cross(l, a), SVMul(l, s))`) generates in 0.58 s at depth 16 and 8.2 s at depth 20.
  The *generated code* stays linear at 30 lines and 192 µs, so CSE works; the cost comes
  from `str` expanding the sharing. That is an engine property (`str` is the identity).
  Every derived system generates in 9 ms or less, so this is not a problem today.
- **Orientation of row results.** A row-vector expression is returned as `(3,)` (NumPy) or
  `3×1` (CasADi), while `evaluate` returns `(1, 3)`. This is documented ("rows as
  3-vectors"). `TimeDerivative(Transpose(x))` is rejected as "a TimeDerivative of a non-leaf",
  although it equals `(ẋ)ᵀ`. This conforms to the spec; no engine output contains it.
- `Number(-0.0)` carries the Zero flag and is written `0.0`. Only the sign of zero differs.
- With the literal `1e20` in the fuzz pool, 1 to 3 of 3000 trees disagree with `evaluate`.
  In every case, perturbing the inputs by 1e-15 (relative) moves the reference itself by
  ~10⁴ against a discrepancy of 0.7, so this is ill-conditioning (cancellation of 1e20
  terms) and not translation. Without `1e20`, 0 of 12 000 disagree. Literal text is
  bit-exact (`repr` round trip, checked for 1e±20, 5e-324, 1e300, 2**70, negatives).

## What held

1. **Translation (12 000 trees).** The generator of `specs/reviews/scripts/engine_fuzz.py`
   was used with variation symbols mapped to plain vectors and `wrap_calc` replaced by
   `t_diff()`, so whole-subtree product rules occur. Leaves cover R (×2), R³ (×2), S² (×2),
   SO(3), two inputs, constant/symmetric/skew/orthogonal/unit-norm parameters, prescribed
   rates `d/dt(M)` (as params), and the literals −2.5, −3, ±1e-20, 0.5, 0, 1. Every node
   class occurred, with between ~350 and ~3 800 nodes per class per 3000 trees. Rows,
   outer products, nested `Transpose`, `Vee`, `Hat`, row `MVMul` chains and
   `Transpose(MVMul(Transpose(v), M))` were all present. NumPy and CasADi matched
   `evaluate` to 1e-9 relative in all 11 820 generatable trees over seeds 1–4. The 180 rejected
   trees contained a `Variation` (the S² δq motif) and were rejected by design. The
   hand-built orientation cases (TT-outer, `(aᵀM)N`, `((aᵀM)N)ᵀ`, `aᵀMb`, `(aᵀM)ᵀbᵀ`,
   `aᵀs·b`, `(aᵀ+bᵀM)a`, `((MN)ᵀ)ᵀ`) are also correct.
2. **CSE and identifiers.** A same-`str` subtree with a different value was not found
   (printing puts parentheses around sums and negative numbers, and `Cross`/`Dot`/`Hat`
   print explicitly). Rows and columns of one symbol print differently (`(a)'`). 41
   awkward names were tested: keywords, builtins (`sum`, `max`, `float`, `range`),
   reserved names (`np`, `ca`, `out`, `state`, `p`, `u`, `a`), helpers and temporaries
   (`_t0`, `_hat`, `_vec`, `_arr`, `_check_params`, `_PARAMS`), module names (`LAYOUT`,
   `copy`, `mass`), LaTeX (`x^\ast`, `\omega_{q}`), unicode (`α`, `ω̇`), quotes,
   trailing backslashes, `1x` vs `s_1x`, and `m` vs `m `. All produce unique, valid
   identifiers that never start with `_`. Leading underscores are stripped, so a name
   cannot shadow a helper or a `_tN` temporary. Both targets give correct values.
3. **StateLayout.** R×R×S²×R³×S²×SO(3)×SO(3) layout:
   - Order and sizes are as listed (29 configuration, 46 state and 17 acceleration entries).
   - `state_from_env` and `env_from_arrays` round-trip exactly.
   - `kinematics` equals a central difference of the exact flow (`R exp(hΩ̂)`, `exp(hω̂) q`) to 6e-11, row-major, on both targets.
   - The `LAYOUT` slices equal the layout's.
   - 1×3 and 3×1 mass blocks between R and R³, S² and SO(3) match `evaluate`, and CasADi transposes the 1×3 block correctly.
   - A decoupled system (a variable absent from L, no inputs) gives zero rows and `(n, 0)` input matrices.
4. **Energy.** Five systems absent from the suite were checked against hand-derived closed forms:
   - SO(3) body with a **non-symmetric** J;
   - **gyroscopic** particle `c ẋ·(b×x)`;
   - **S²×SO(3)** body with a pendulum hung at ρ;
   - variable inertia `½(1+s²)ΩᵀJΩ` on R×SO(3);
   - two coupled spheres with a linear `ω₁·b` term.

   `E` agrees to 1e-12. `dE/dt` (central difference along `[kinematics, accelerations]`)
   equals the input power (`F·ẋ`, `M·Ω`, `f·ω`, `f_s ṡ`) with random nonzero inputs, to
   ≤ 9e-10 relative.
5. **Generated modules.** The 18 sources are byte-identical under `PYTHONHASHSEED` 0, 1,
   999, 4242 and random. The CI regeneration step runs in a numpy-only venv (as in the lint
   job), and it fails on a hand edit of a committed file and on a codegen change without
   regeneration. Importing without geomech is already covered by
   `test_generated_files_run_without_geomech`.
6. **Performance.** Generation of the nine systems takes 1–9 ms. There is no textual CSE
   miss in any checked-in function, and shared trees produce linear code.

## Fix log

Fixed on `feature/codegen-simulation` on top of 83f80d6 (uncommitted). The six strict
xfails of `tests/test_codegen_adversarial.py` were flipped by deleting their markers; no
assertion was changed. New tests are in `tests/test_codegen_arguments.py` (59 tests).
Repo venv (Python 3.12.13, numpy 2.5.4, casadi 3.8.2, scipy 1.18.1).

| # | fix | where |
|---|---|---|
| F1 | `_arr(None, n, name)` raises `ValueError("u must have 1 entries, got None (u is required here)")` in both targets. A function converts only the arrays its expressions use, so `u=None` (or `a=None`) is accepted exactly when the function does not depend on it: always for a system with no inputs. The generated code never reads `None` as zero. That shorthand stays in `geomech.simulation` (`accelerations`, `simulate`, `lqr`), which documents `u=None` as zero input and converts it with `Model.as_input` before calling the module; `mass`/`bias`/`input_matrix`/`energy` take no `u`. The message keeps the old `"<name> must have n entries"` prefix, so `test_codegen_layout.py::test_acceleration_needs_a` passes unchanged. Checked for `residual`, `M_lin`, `A_lin` and `B_lin` of all nine VBL systems on both targets: `None` raises iff the function's source converts that array, and otherwise gives the same value as a zero array. | `geomech/codegen/numpy_target.py:38-40`, `casadi_target.py:25-27` |
| F2 | `_scalar(x, name)`: Python/NumPy scalars and 0-d or single-element arrays/sequences give a `float`; more (or zero) entries raise `ValueError("parameter 'm_1' is a scalar, got 2 entries")`. The CasADi version passes a 1-entry `SX/MX/DM` through as `1×1`. Every scalar parameter is unpacked through it (`m = _scalar(p["m"], "m")`). | `numpy_target.py:47-51`, `:121-126`; `casadi_target.py:35-44` |
| F3 | `_number` writes `np.inf`, `-np.inf`, `np.nan`; the CasADi target maps `±np.inf` to `±ca.inf` (CasADi has no `nan` constant, so `np.nan` stays). `is_number_literal` treats `nan` as equal to its printed `(nan)`, so `Number(nan)` is inlined and `collect_params` leaves it out. | `numpy_target.py:87-97`, `casadi_target.py:106-107`, `geomech/core/base/expressions.py:784` |
| F4 | `regenerate_generated.py` deletes stale `*.py` in write mode. `--check` reports `out of date`, `missing`, `stale (no longer generated)` and, inside a git work tree, `not tracked by git` (`git ls-files`; skipped outside a repo or without git), and exits 1 on any of them. The CI lint step is now `python scripts/regenerate_generated.py --check`. A fresh CI checkout holds only committed files, so a forgotten new output shows there as `missing` and a removed system as `stale`. | `scripts/regenerate_generated.py:83-143`, `.github/workflows/tests.yml:33` |
| F5 | NumPy `Cross` is written `_cross(a, b)`, a module-level helper of three multiply-subtracts. It equals `np.cross` bit for bit: 1500 random pairs at scales 1, 1e±200, shape `(3,)`. The helper is only ever called on 3-vectors. `examples/generated/` was regenerated: 63 `np.cross` calls became `_cross`, and some locals moved because `_cross(` is shorter than the line-split threshold. Every function of all nine standard-form modules returns bit-identical values old vs. new on 20 perturbed states each, and 0.5 s simulations of the nine systems give identical trajectories (max diff 0.0). The CasADi target keeps `ca.cross`. | `numpy_target.py:70-73`, `:148` |
| F6 | `casadi>=3.6,<4` | `pyproject.toml:29` |
| advisory | The emitter's CSE key is structural: `(class, str)` for a leaf and `(class, operand keys)` for an operation, interned to an int and memoised per emission pass in a dict keyed by `id(node)`. The dict also holds the node, so an id cannot be reused. `str` is never taken of an operation node (`_time_derivative` checks `base_symbol` before the `str`-keyed slot lookup). The 18 checked-in modules are byte-identical to what `(class, str)` keys produce (checked by monkeypatching the key). | `geomech/codegen/emitter.py:108-112` (state), `:153-166`, `:248-250` |

Timings, same machine, best of 3–5 runs (scratch script, double rigid pendulum from
`tests/simulation_systems.SETUPS`):

| measurement | before | after |
|---|---|---|
| `double_rigid_pendulum.bias` | 272 µs | 45 µs |
| `double_rigid_pendulum.energy` | 88 µs | 23 µs |
| `double_rigid_pendulum.mass` (no cross product) | 15.5 µs | 16.3 µs |
| `simulate(..., 1.0 s, rkmk4, dt=1e-3)` | 1.36 s | 0.45 s (final state identical) |
| `to_source`, depth-20 shared tree `VAdd(Cross(e, a), SVMul(e, s))` | 7991 ms | 0.4 ms |
| `to_source`, chain of depth 800 (alternating `VAdd(Cross(e, a), x)` / `SVMul(e, s)`) | 730 ms | 4.5 ms |
| full suite, `pytest tests/ -q -W error` | 82.4 s (3192 passed, 16 xfailed) | 57.1 s (3257 passed, 10 xfailed) |

The suite is 25 s faster, mostly from the simulation tests through the cheaper `bias`.

Other changes:

- `tests/test_codegen_translation.py::test_common_subexpressions_computed_once` counted
  `"np.cross("`. It now counts `"_cross("`, with the same assertion (one cross product).
  This text change follows from F5. The docstring of
  `tests/test_codegen_generated.py` now names `--check`.
- `docs/codegen.md`: required arguments and `None`, scalar parameter forms, the
  `_cross`/non-finite rows of the translation table, structural CSE, `--check`.
  `specs/codegen.md` translation row for `Cross`. `geomech/codegen/source.py`: the
  `to_source` docstring lists the helpers and the `None` rule. CHANGELOG "Fixed" entry.
  The casadi pin was corrected under "Added".
- Checks run: ruff format/check on `geomech/ tests/ scripts/ examples/generated/`
  (clean), `mypy geomech --ignore-missing-imports` (clean), full suite under `-W error`,
  `mkdocs build --strict` (clean), and `scripts/regenerate_generated.py --check` (exit 0).
  `test_generated_source_is_independent_of_hash_seed` (seeds 0, 1, 4242) passes, and the sources hash identically under PYTHONHASHSEED 0, 1, 999, 4242 and random, so the
  generated modules are byte-stable.
