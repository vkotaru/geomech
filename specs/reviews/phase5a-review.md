# Phase 5a adversarial review: simplifier completeness and D1

Branch `fix/phase5a-simplifier`, commit cbd3653, range `47cce2d..cbd3653`. Reviewer tests:
`tests/test_phase5a_adversarial.py` (50 passing, 7 strict xfails). Nothing committed.

## Verdict

**GO WITH CHANGES.** I found no wrong equation, no wrong standard-form block and no
linearization row that differs on admissible perturbations. Every rule firing I captured
holds along curves on the constraint manifold, including the new firings inside
`TimeDerivative`/`Variation`. The findings below are one robustness regression (F1, a
crash outside the pipeline), two edge-case or cosmetic regressions (F2–F4) and one
tooling bug from an earlier phase (F5) that weakens one of the implementer's gate numbers.
None of them blocks the phase. F1 should be fixed, or recorded as a known limitation,
before 5.5/5.6 build on `rebuild`.

## Findings (ranked)

Severity order: wrong equation > wrong linearization on admissible perturbations > rule
firing outside its hypothesis > equality change > cosmetic. **No finding is in the first
two classes.**

### F1 — simplifying a calculus node whose operand reduces to a manifold leaf raises (robustness regression)

`vector_rules` now recurses into `TimeDerivative`/`Variation` (5.1). When the operand
reduces to a bare `S2` leaf, `map_children` rebuilds `TimeDerivative(q)` or
`Variation(q)`, and the constructor rejects it with `ValueError`
(`geomech/core/transformations/utils.py:50`, reached from `vector_identities.py:154`).

```python
q = S2("q")
full_simplify(TimeDerivative(SVMul(q, Dot(q, q))))          # 47cce2d: d/dt(q Dot(q,q)); 5a: ValueError
full_simplify(Variation(SVMul(q, Dot(q, q))))               # 47cce2d: δ(q Dot(q,q));   5a: ValueError
full_simplify(TimeDerivative(VAdd(q, Cross(q, Cross(q, q))))) # same
```

The rewrite itself is sound (`q·q = 1` holds identically, so it may be applied under
d/dt). Only the rebuild is wrong: it should give `q.t_diff()`/`q.delta()` (the eager
kinematic tree) rather than construct the forbidden node. This crash class already
existed (`d/dt(q·1)`, `d/dt(RᵀR q)` and `δ(I q)` raise at 47cce2d too, via `_eliminate`).
5a makes more inputs reach it. The default fuzz finds 2 such inputs in 8000 trees
(1 at base). `compute_eom` and `linearize` cannot reach it, because they take `.delta()`
and `.t_diff()` eagerly. Strict xfail:
`test_calculus_node_over_reduced_manifold_leaf_does_not_raise` (3 cases).

### F2 — integer-valued floats are normalised only in scalar chains (equality inconsistency)

`_exact` is applied in `mul_fold`, `add_fold` and the like-term rebuild. It is not
applied to a vector/matrix term with one numeric factor, or to a bare `Number` leaf:

```python
full_simplify(x * Number(2.0))   # x(2)      == full_simplify(x * Number(2))
full_simplify(v * Number(2.0))   # v(2.0)    != full_simplify(v * Number(2))  -> v(2)
full_simplify(Number(2.0))       # (2.0)
```

Values are unaffected, and like-term merging still sums such coefficients. Under the
`str` identity, though, the same term compares equal when it is a scalar and unequal
when it is a vector. That contradicts the CHANGELOG claim that "a float with an integer
value becomes an int". Fix: normalise in the `SVMul`/`SMMul` cases of `_eliminate`
(`simplification.py:211`, `:220`), or document the limit. Strict xfail:
`test_integer_valued_float_normalised_for_vectors_and_matrices`.

### F3 — `is_number_literal` raises `OverflowError` for an int beyond the float range

`expressions.py:685` evaluates `float(text) == float(expr.value)`. `Number(10**400) * x`
was folded at 47cce2d. It now raises `OverflowError` from `simplify`. This is an edge
case. Fix: compare `int` values exactly, or also catch `OverflowError`. Strict xfail:
`test_huge_integer_literal_folds`.

### F4 — `repr` of a named constant with a value drops the `*` marker (cosmetic)

`printing.py:120` (`_repr_leaf`) prints every `is_numeric` leaf by bare name.
`repr(Scalar("m", value=2, attr=["Constant"]))` is `m` and the valueless one is `m*`.
Since D1 these symbols survive into equations, so `print_eom`/`repr` output now shows
them unmarked (`v:x·m·(2)`). Use `is_number_literal` there as `_numeric_value` does.
Strict xfail: `test_repr_of_named_constant_with_value_keeps_marker`.

### F5 — `engine_fuzz.py` does not fuzz its own tree (tooling bug from an earlier phase)

`specs/reviews/scripts/engine_fuzz.py:55` puts `dirname×3(__file__)` first on `sys.path`.
That is `specs/`, not the repository root, so `import geomech` resolves to the editable
install. Run from a worktree of 47cce2d, it printed "fuzzing geomech from
…/pygeomech/geomech" and fuzzed **the phase-5a tree**. My first base run had to be
redone with `PYTHONPATH=<worktree>`. This does not affect the engine. It does affect the
reproducibility of the report's "249 distinct findings vs 202 at base" line: with seed 77
and a correct base I measured 237 vs 186. Strict xfail:
`test_engine_fuzz_imports_its_own_tree`.

### Notes (no test)

- `tests/test_phase1_adversarial.py:568` is still named
  `test_scalar_with_value_is_folded_by_engine_and_used_by_evaluator`, but since D1 it
  asserts the opposite (`"m" in …`). The name is stale.
- D1 is a semantic change that the fix plan lists under "not changed without a yes". I
  reviewed it as approved. It is consistent everywhere I looked (below). The visible
  consequence: `Scalar("a", value=2)·q̇₁·q̇₂ − Scalar("b", value=2)·q̇₂·q̇₁` used to
  cancel. Now it keeps both terms, and `to_standard_form` gets an extra (numerically
  zero) mass block.

## What held

### 1. Independent numeric equivalence (own harness, not `pinned_snapshot.py`)

I used my own git worktree of 47cce2d (removed afterwards) and my own harness.
Constants are drawn from a generator seeded by the leaf's printed name, so values do not
depend on how a tree is written. States come from `sample_state`, and S² equations are
compared on the tangent plane. Coverage: 115 systems, made of 67 builder systems
(`test_eom_oracle.SYSTEMS`, `phase1_systems.ADVERSARIAL_SYSTEMS`, the two outer-product
reproducers, `stiff_point_mass(k=3)`) and **48 new Lagrangians**:

- S² with non-symmetric, skew, `Iᵀ` and `Hat·Hat` matrices
- tangency identities under `MVMul`/`Hat`
- row forms
- named values and float-integer literals
- a time-varying parameter `l(t)`
- double and skew-coupled S²×S²
- SO(3) with `RᵀRJ`, `RᵀKR`, `Vee(ARᵀ − RAᵀ)`, `Ω̂JΩ̂`, `Ω̂ᵀJΩ̂`, row products, outer products, scaled hats and scaled rows
- SO(3)×SO(3) and SE(3) with an offset centre of mass
- R, R², R×S², R×SO(3), R³ with skew and row forms, R³×S²×SO(3)

| check | result |
|---|---|
| systems derived by both trees | 108 (5 raise in both with the same error class; 2 that raised at base now derive: `x_so3_scaled_hat`, `x_so3_row_scaled`, `check_eom` 9.6e-9 and 1.3e-8) |
| equations compared / strings changed | 143 / 126 |
| equation values, 4 states each, worst relative difference | 1.6e-15 (S² projected; the unprojected values are also identical) |
| `to_standard_form` blocks compared, worst relative difference | 1868, 5.1e-16 |
| linearization blocks compared / differing as matrices | 8648 / 24 (= the 6 S² blocks × 4 states the implementer reported; no other system) |
| linearization rows along admissible directions (3 per state) | 3384, worst relative difference 2.0e-15 |
| `check_eom`, new tree, 110 systems | worst 1.9e-7, same as base |
| `finite_difference_check`, new tree, 108 linearizations | worst 7.8e-10 |

### 2. The six changed S² linearize blocks

For `s2_anisotropic_qdot`, `adv_omega_and_qdot_mixed` and `adv_omega_dot_S_qdot` I took
the A blocks at 4 states and 3 admissible directions each (ξ ⊥ q,
δω = ξ × ω + w with w ⊥ q):

| system | max ‖(A_new − A_old)·p‖ admissible | max ‖D_δω (I − qqᵀ)‖ | with δω + 0.7 q (not admissible) |
|---|---|---|---|
| s2_anisotropic_qdot | 1.5e-14 (blocks of size up to 7) | 2.0e-15 | 3.1 |
| adv_omega_and_qdot_mixed | 3.6e-15 | 1.3e-15 | 1.1 |
| adv_omega_dot_S_qdot | 4.0e-15 | 2.4e-15 | 2.9 |

The implementer's claim is confirmed. The difference vanishes on admissible
perturbations (≤ 1e-12), the δω column difference has the form `u qᵀ`, and it shows up
only when the perturbation leaves the tangent bundle. None of my 21 S²-containing
systems has a block that differs, not even as a raw matrix.

### 3. Each new rule in isolation (per-firing soundness)

I wrapped `rule()` in `vector_identities`, `pulling`, `simplification` and `expansion`,
with tracing switched on so that the `*_fold` rules fire too. Each `(before, after)`
pair was then compared as **full (t, ε) Taylor jets**, using `engine_fuzz.JetEvaluator`
with every coefficient compared, not just the value at t = 0. An identity that held only
at an instant would therefore fail. The inputs were 1500 random fuzz trees plus 1500
targeted trees: each new rule's pattern, wrapped in `d/dt`, `d²/dt²`, `δ` and `d/dt δ`,
and placed under `Dot`, `MVMul`, `Hat`, `Transpose`, `Vee` or `SVMul`.

- Consistent mode: **23 630 firings over 66 rule names, 0 mismatches.** The only
  residual is 9.8e-8 on `vector.dot_unit_norm`, which is the jet-normalisation
  truncation. The phase-5a rules were all exercised (`pull.transpose_scaled` 167,
  `pull.hat_scaled` 60, `pull.vee_scaled` 45, `pull.row_mvmul_scaled` 44,
  `pull.vvmul_scaled` 59, `expand.vee_over_madd` 42, `simplify.transpose_skew` 94,
  `simplify.mmmul_reassociate` 33, `simplify.mmmul_combine` 52, the four folds
  1,500–2,500 each; counts from the 400-case run).
- Adversarial mode (flags declared but violated): mismatches occur only in rules that
  depend on a flag or constraint. Those are tangency (`dot_tangent_orthogonal`,
  `s2_tangent_recovery`/`negated`, `s2_double_cross`), `dot_unit_norm`,
  `orthogonal_transpose`, `transpose_symmetric`, `transpose_skew` and `hat_vee_skew`.
  Associativity, identity and product rules mismatch only when their `after` embeds one
  of those nested rewrites. I checked every worst case: `(Q)'Q…`, `R(R)'…` and `(K)'`
  always appear inside. No `pull.*`, `expand.*`, fold or like-term rewrite mismatched on
  its own.
- Hypotheses are checked in the new contexts. `transpose_skew` does not fire on `KK`,
  `â b̂`, `J + K`, `JK` or `d/dt(K)`; their `is_skew_symmetric` is false, which is
  conservative. Tangency is by name and does not look inside `d/dt(ω)`. Under a calculus
  node, a tangency or unit-norm rewrite is applied to an operand for which the identity
  holds **identically along the curve or variation**, so the rewrite commutes with
  d/dt and δ. This is confirmed exactly by the jet comparison.
- The engine fuzz (`--cases 3000 --seed 77`) shows 0 wrong-value findings in consistent
  mode on both trees. Distinct findings went from 186 (base, run correctly; see F5) to
  237. The +51 are all `flag-dependent` in the adversarial campaign (`vector_rules` 11→29,
  `simplify` 39→55, `full_simplify` 40→57), plus the one extra F1 input.

### 4. Like-term merging

The multiset key compares bases with order-sensitive equality for `MMMul`, `MVMul`,
`VVMul`, `Cross`, `SVMul` and `SMMul`. Only the flattened `Mul` factors and `Dot` are
treated as commutative. I tried 7 counterexamples, all held (not merged; value kept):

- `AB s − BA s`
- `(a×b)m + (b×a)m`
- `(a×b)·c m + m (b×a)·c`
- `aᵀAb m − m bᵀAa` with A non-symmetric
- `A a m − Aᵀ a m`
- `â b̂ m − b̂ â m`
- `a bᵀ m − b aᵀ m`

The valid merges are made: `x(2)y + xy(3) → xy(5)`, `Dot(a,b)m + m Dot(b,a) →
Dot(a,b)m(2)`, and the three `k m (a×ω)·b` copies in a pipeline merge to `km(-4)`.

### 5. D1

`Scalar("m", value=2)`, `Scalar("l", value=1)` (not dropped as a one) and
`Scalar("z", value=0)` (not dropped as a zero) stay by name through `expand`, `pull`,
`vector_rules`, `simplify`, `full_simplify`, `compute_eom`, `to_standard_form`
(`M = Imll(-1)`), `linearize` (blocks contain `mgl`), `to_latex` (`-k x`) and `repr`.
`evaluate` uses the values: `px + rx` with p = r = 2 gives 6.0 at x = 1.5, and
`check_eom`/`finite_difference_check` pass. Two named scalars with equal values do not
merge (`px + rx`, `p + r(-1)`), the same name does (`px − xp → 0`), and a named value
does not merge with a literal (`px + x(2)`). At 47cce2d these gave `x(4)` and `0`.

### 6. Printing

`x·(-1.0)` and `x·(-1)` simplify to the same `x(-1)`. A symbol `Scalar("(-1)")` with no
value is not a literal: it is neither folded nor dropped, and `repr` prints it as a name,
not as a minus sign. An equality that used to fail now holds: `simplify(x·2.0) ==
simplify(x·2)`. I found none that used to hold and now fails, apart from the F2
asymmetry. Pipelines compare by `str` only at fixed-point checks, and those converge
(fuzz convergence 2997/2997 non-raising cases).

### 7. Regression

| check | result |
|---|---|
| suite `-W error` (with the reviewer file) | 2221 passed, 9 xfailed (2171 + 2 without it, as reported) |
| suite `GEOMECH_DEBUG=1 -W error` | 2221 passed, 9 xfailed (with the reviewer file) |
| suite `GEOMECH_TRACE=1` | 2171 passed, 2 xfailed |
| `engine_fuzz.py --cases 3000 --seed 77` | 0 wrong values in consistent mode; 237 distinct findings, all classes above |
| examples (6) | all run |
| notebooks (4), executed into a temporary directory | all run; outputs identical to the committed ones |
| timing, median `compute_eom` (ms), base → 5a | quadrotor+load 10.9 → 11.8; double rigid pendulum 17.6 → 20.1; double spherical pendulum 22.2 → 17.6; R³×S²×SO(3) 13.5 → 14.4; SO(3)×SO(3) 8.0 → 8.6; worst 1.14×. `linearize`: worst 1.2×, two systems faster |
| `ruff format --check`, `ruff check` on the new test file | clean |

### Reviewer artefacts

The harness lived in a temporary directory and is not part of the repository:
`run.py <tree> <out.pkl>`, `compare.py`, `s2blocks.py`, `fuzz_rules.py <tree> <cases> <seed>`,
`my_systems.py`. The 47cce2d worktree was removed after the runs.

Note (orchestrator): D1 was approved by the author on 2026-11-13 (option "Stay symbolic"),
recorded in the Decisions section of `specs/fix-plan.md`.

## Fix log

Fixes on `fix/phase5a-simplifier` on top of cb4697b (not committed). The 7 strict xfails
were flipped by removing their markers; no assertion was changed.

- **F1.** New `calculus_node(cls, operand)` in `geomech/core/transformations/utils.py`:
  a `TimeDerivative`/`Variation` whose operand is a bare `S2`/`SO3` leaf is built in eager
  form (`q.t_diff()`/`R.t_diff()`, `q.delta()`/`R.delta()`), anything else as `cls(operand)`.
  `rebuild` uses it, and so do the two places that construct the nodes directly:
  `simplification._eliminate` (the pre-existing `d/dt(q·1)`, `δ(I q)`, `d/dt(RᵀR q)` path)
  and `expansion.expand`. New test `test_calculus_node_over_reduced_manifold_leaf_is_eager`
  (8 cases: `d/dt` and `δ` of `q·1`, `q (q·q)`, `RᵀR q`/`I q`, `R·1`/`I R`) runs every pass
  and checks that the result equals `full_simplify` of the eager form.
- **F2.** `_exact_literal` normalises the numeric factor of `SVMul`/`SMMul` in
  `_eliminate`: `v·(2.0)` → `v(2)`, `A·(-1.0)` → `A(-1)`. Bare `Number(2.0)` leaves are
  left as they are (out of scope; no pipeline produces one).
- **F3.** `is_number_literal` compares an `int` value with `int(text)` exactly and falls
  back to `float` only otherwise, catching `OverflowError`. `Number(10**400)·x·Number(10**400)`
  folds to `x(10**800)` exactly (new test).
- **F4.** `_repr_leaf` prints a leaf bare only when `is_number_literal`; a named scalar
  with a value gets the `*` (also without the `Constant` attribute, since a value makes it
  constant in time: `delta()`/`t_diff()` return zero). New test pins `m*`, `(2)`, `0.5`,
  `v:v·m*`.
- **F5.** `engine_fuzz.py` inserts `os.path.abspath(os.path.join(__file__, *[os.pardir] * 4))`,
  the repository root. `double_rigid_pendulum_numeric_check.py` had no `sys.path` insert at
  all (same bug: it imported the installed `geomech`); it now inserts the same root.
  `pinned_snapshot.py` and `euler_lagrange_check.py` were already correct (`../../..` from
  the script directory). Verified from a git worktree of cbd3653 with the fixed scripts
  copied in, run from a temporary directory: `engine_fuzz.py --quick` printed "fuzzing geomech from
  <worktree>/geomech", and all four scripts resolve `geomech` to the worktree.
- **Note.** `tests/test_phase1_adversarial.py` test renamed to
  `test_scalar_with_value_stays_symbolic_in_engine_and_is_used_by_evaluator`; body
  unchanged. Older reports that cite the old name were left as written.
- CHANGELOG: four entries under Unreleased / Fixed.

| check | result |
|---|---|
| `pinned_snapshot.py compare` against a snapshot taken at cbd3653 | 510 records, 3727 outputs unchanged, **0 strings changed**, 0 failures, `check_eom` 198 compared, 0 worse |
| suite `-W error` | 2238 passed, 2 xfailed (the 2 pre-existing) |
| suite `GEOMECH_DEBUG=1 -W error` | 2238 passed, 2 xfailed |
| `ruff format --check`, `ruff check` on `geomech/ tests/ specs/reviews/scripts/` | clean |
| `mkdocs build --strict` | builds |
| `engine_fuzz.py --quick` | 0 findings in consistent mode; adversarial findings are flag-dependent only |
| `engine_fuzz.py --cases 3000 --seed 77` | 0 consistent-mode findings; 0 manifold-leaf `ValueError`; 209 distinct findings (115 flag-dependent, 89 unsupported-linear-term, 5 documented-error-on-legal-input) |
| examples (6) | all run |
