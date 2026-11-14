# Adversarial review of phase 1 (`specs/fix-plan.md` 1.1–1.6)

Branch `fix/phase1-equations`, range `1dd98ab..c96d31c` (636e963 numeric oracle, c96d31c engine
fixes). Reviewed: `git diff 1dd98ab..c96d31c -- geomech/` in full, the two phase-1 reports, the
new test files, `docs/numeric-verification.md`. Nothing under `geomech/` and no existing test was
modified. New file: `tests/test_phase1_adversarial.py` (88 passing tests, 7 strict `xfail` with
`F<n>` in the reason). Scratch scripts and the old-engine worktree lived in a temporary
directory; the numbers below are from those runs.

Baseline on this machine: `python -m pytest tests/ -q` → 1441 passed, 5 xfailed; with the new
file → **1529 passed, 12 xfailed**; `ruff format --check` + `ruff check` clean.

## Verdict: GO WITH CHANGES

The mathematics is right. The −2 coefficient is confirmed three independent ways and every
adversarial Lagrangian I could build gives the Euler–Lagrange residual. Two things must change
before merge, both small: `collect` now recurses forever on a legal Lagrangian (F1, a regression
in failure mode from a typed error), and the oracle's default step sits in the roundoff-dominated
regime so a correct equation can fail the 1e-5 gate (F3). F2 is a documented limitation that the
docs currently deny ("wrong equations are off by order one").

---

## Findings, ranked

Severity order used: silent wrong equation > oracle blind spot > crash on legal input > missed
error > cosmetic. **No silent wrong equation was found.**

### F1 — `collect` loops to `RecursionError` when the target is inside a matrix factor (crash on legal input; regression)

**Reproducer** (`tests/test_phase1_adversarial.py::TestCollectPeel::test_target_inside_matrix_factor_terminates`,
`::TestPipelineReachesCollectLoop`):

```python
collect(Dot(MVMul(Hat(v), a), b), v)                                    # RecursionError
L = HALF * Dot(Om, J * Om) + Dot(a, MVMul(VVMul(Om, Transpose(c)), b))  # = ½Ω·JΩ + (a·Ω)(c·b)
compute_eom(L, Dot(eta, M), SystemVariables(matrices=[R]))               # RecursionError
```

**Observed vs expected.** At c96d31c both raise `RecursionError` (no message). At 1dd98ab the
first returned the Dot unchanged and the second raised
`AlgebraicError: '\eta_{R}' appears only inside the calculus operator 'd/dt(\eta_{R})' …`,
i.e. a typed error naming the term. Expected: either the exposed form or a typed error.
`CHANGELOG.md:95` ("`collect` no longer recurses forever on `Dot(M a, N b)` when the target is
nested …") is now false for this shape. My random-Dot fuzz (2000 cases, target placed in vector
*or* matrix slots): 930 exposed, 928 unexposed-but-value-preserving, **142 loops**, 0 wrong
values; every loop has the target inside `Hat(·)`, `v cᵀ` or `MMMul(·, Hat(v))` under an `MVMul`.

**Diagnosis.** `geomech/core/math/collect.py:92-96` (`(MVMul(), Cross())`), `:120-121`
(`(MVMul(), MVMul())`) and `:131-138` (`(MVMul(), _)`) peel on `l.has(vec)`, which is also true
when `vec` is in `l.left` (the matrix). The peel then moves the matrix *containing the target*
to the other side as `MVMul(Transpose(l.left), r)`; `:142-144` flips it back; the next peel
wraps another `Transpose`. The docstring's termination argument (`:48-52`, "each peel removes
one matrix layer from the side holding vec") is wrong for this case: the layer is moved, not
removed. The pipeline reaches it because `simplify` turns `MVMul(Hat(x), y)` into `Cross(x, y)`
but leaves an outer product `MVMul(VVMul(x, cᵀ), y)` alone, so a Lagrangian with
`a·((Ω cᵀ) b)` (or the S² twin `a·((ω cᵀ) b)`) survives to `collect` with `Ω̇`/`ξ̇` inside the
matrix. `specs/reviews/scripts/engine_fuzz.py:1151-1182` never generates the target inside
`MVMul.left`, which is why its `collect[target occurs once] 2429/2429` did not see this.

**Fix.** Peel only when `l.right.has(vec)`; when `l.left.has(vec)` (target inside the matrix
factor) raise `AlgebraicError("collect: '<vec>' occurs inside the matrix factor of '<expr>'; no
rotation exposes it")`, or return `expr` unchanged and let extraction raise its existing
calculus-operator error as 1dd98ab did. Correct the termination docstring, and add the shape to
the fuzz generator.

### F2 — `check_eom` cannot see a term that is small relative to the whole equation (oracle blind spot)

**Reproducer** (`::TestOracle::test_flipped_term_next_to_a_stiff_term_is_detected`,
`::test_dropped_term_next_to_a_stiff_term_is_detected`): point mass with
`k = Scalar("k", value=1e7)`, `L = ½ m ẋ·ẋ − ½ k x·x − m g x·e₃`, `δW = δx·f`.

| equation | `check_eom` |
|---|---|
| correct `−m ẍ − 1e7 x − m g e₃ + f` | 4.1e-7 |
| gravity term sign flipped | **5.7e-7 (accepted)** |
| gravity term deleted | **4.0e-7 (accepted)** |
| forcing term `f` sign flipped | **6.3e-7 (accepted)** |
| same with `k` random (no `value`), gravity flipped | 1.74 (detected) |

**Diagnosis.** `geomech/numeric.py:614-619` `_relative_error` divides by `max(1, |P E_true|)`;
with one O(1e7) term every O(1) term is below the 1e-5 gate. `Scalar(..., value=…)` is honoured
by the evaluator (correctly, the engine folds it), so a user who puts physical constants in
`value` gets this silently. `docs/numeric-verification.md` says "Wrong equations are off by
order one" and "there is no grey zone", which is only true for O(1) parameters.

**Fix.** Document it in "What it does not verify", and give the user a power check: for example
`assert_eom(..., check_power=True)` that flips each additive term of `E_lib` in turn and
requires the error to exceed the gate (this is exactly `test_check_eom_detects_a_flipped_term`,
generalised), or report the resolution `1e-5 · max|E_true|` in the assertion message.

### F3 — the default `h=1e-5` is roundoff-dominated; a correct equation fails the 1e-5 gate (oracle false positive)

**Reproducer** (`::TestOracle::test_correct_equation_scores_well_below_the_gate_at_default_h`,
`::test_default_step_is_near_the_noise_minimum`): `L = ½ m |ẋ + l q̇|² + ω·(x × c)` on R³ × S²
(`adv_r3_s2_pendulum_with_gyro`), a correct equation (trajectory and residual agree at h=1e-4).

| `h` | `check_eom`, seeds 0 / 1 / 2 |
|---|---|
| 1e-3 | 7.0e-6 / 8.2e-7 / 1.1e-5 |
| 1e-4 | 1.2e-7 / 9.4e-9 / **1.6e-7** |
| 1e-5 (default) | 5.9e-6 / 6.9e-7 / **1.43e-5 → fails the gate** |
| 1e-6 | 4.6e-4 / 6.0e-5 / 5.4e-4 |
| 1e-7 | 3.3e-2 / 3.7e-3 / 5.4e-2 |

Over the 36 oracle systems and seeds 0–9 the worst correct score is 8.3e-6 at h=1e-5
(`s2_anisotropic_qdot`; `quadrotor_with_load` 6.1e-6, `prismatic_slide` 6.1e-6) versus 2.3e-7
at h=1e-4. Flipped terms are still detected at h=1e-4 (free rigid body: 3.5, 3.1, 1.0).

**Diagnosis.** `geomech/numeric.py:548-561` nest two central differences with the same `h`
(`_ddt` of `_grad`), so roundoff is ≈ ε·|L|/h² ≈ 1e-6 at h=1e-5, within 10× of the gate; the
"noise ≈ 1e-6" in the docstring (`:691`) is this roundoff, not truncation. The gate tests are
deterministic (seed 0, 5 trials) so CI will not flake, but the next system an author adds can
fail spuriously, as this one does at seed 2.

**Fix.** Default `h=1e-4` (100× margin, measured), or use a larger outer step for `_ddt`; update
the docstring and `docs/numeric-verification.md` ("about 1e-6 … no grey zone").

### F4 — an unpinned printed equation changed and is not listed (cosmetic)

Old vs new on all 36 oracle systems (3 seeds each): 15 changed. The engine report lists 14
(9 S²/collect strings + 5 crash→ok). The 15th is `vbl_systems.quadrotor_with_load`: the load
equation's term `Cross(ω,Cross(q,ẋ_L)) m_Q l (−1)` became `(−2)`. It is parallel to q
(ω×(q×ẋ) = q(ω·ẋ)), so the oracle, `test_vbl.py` and the notebook are unaffected, and the
report's "none of the pinned strings changed" is true; but the printed output of a shipped
notebook system did change and the CHANGELOG entry ("pinned pendulum equations never changed")
should say so. Full table of all 36 systems, before/after errors, in the Verdict section.

### F5 — nothing tells a user to project an S² equation before solving it (documentation)

Every S² system's `E` has an O(1) component along q (measured `|q·E|/|E|` 0.3–1.5 on 20
systems; it comes from `f`, `ω × c`, `A ω̇`, …). The oracle is right to ignore it, but
`to_standard_form` returns `M = −A` (anisotropic case) or `M = −m l² I` with an `f` that has a
normal part, and solving `M ω̇ = −f − G u` as a 3×3 system gives an `ω̇` with `q·ω̇ ≠ 0`, i.e.
a trajectory that leaves the tangent plane. `docs/vbl.md:179` says this for the linearization;
`docs/numeric-verification.md` ("such a term has no physical content") and the standard-form
docs do not. One sentence in the conventions/standard-form docs: "on S², solve
`P(M ω̇ + f + G u) = 0` with `ω̇ ⊥ q`".

### F6 — `SkewSymmMatrix` mutates the caller's `attr` list (missed error, pre-existing pattern)

`geomech/core/base/expressions.py:466-468` appends to the list the caller passed, exactly as
`SO3` does: `shared = ["Constant"]; SkewSymmMatrix("K", attr=shared); Matrix("J", attr=shared)`
gives `J.is_skew_symmetric == True`. Plan item 2.4 covers it; noting that the fix in 1.3 kept
the pattern rather than copying the list.

### F7 — `MAdd()` with no nodes is vacuously skew-symmetric and symmetric (cosmetic)

`geomech/core/operations/addition.py:86-93` uses `all(...)`; `MAdd().is_skew_symmetric` is
`True`. Not reachable from the pipeline (`_simplify_madd`/extract return `ZeroMatrix` for empty
lists), but `any` + non-empty would be the honest predicate.

### F8 — duplication and naming (cosmetic)

`_derive`, `_equation`, `TOL`, `HALF` are copied between `tests/test_s2_variation_fix.py` and
`tests/test_collect_fix.py`; `HALF` is defined in four test modules; `tests/test_eom_oracle.py`
re-derives 12 systems that `tests/test_manifold_pipelines.py` and `tests/test_standard_form.py`
already define inline (acknowledged in the oracle report; a shared `tests/systems.py` would end
it). `test_pinned_equation` is used four times (disambiguated by class name only).

---

## What was checked and held

### 1. S² rule (−2)

- **Hand derivation.** ω = q×q̇, δq = ξ×q, ξ ⊥ q: δω = (ξ×q)×q̇ + q×(ξ̇×q) + q×(ξ×q̇)
  = q(ξ·q̇) + ξ̇ − q(q·ξ̇) + 0 = ξ̇ − 2q(q·ξ̇) (using ξ·q̇ = −q·ξ̇); and ω×ξ = q(q·(ω×ξ))
  = q((q×ω)·ξ) = q(−q̇·ξ) = q(q·ξ̇). So δω = ξ̇ − 2 ω×ξ. Lee et al. (5.13),
  δω = −S(ω)ξ + (I−qqᵀ)ξ̇, gives the same once qqᵀξ̇ = ω×ξ. The module docstring's
  derivation is complete and matches.
- **Numeric derivation without the library** (`TestS2VariationFormula`): 10 random curve
  families q(t,ε) = n(p(t)+εr(t)); ∂ω/∂ε vs ξ̇ − 2ω×ξ and vs (5.13): < 1e-6; the old rule is
  off along q by > 1e-3 and equal on the tangent plane.
- **Trajectory harness** (`TestTrajectoryAgreement`, scratch run T = 2 s, dt = 1e-3, RK4):
  engine's `E = 0` solved for ω̇ on the tangent plane vs the hand-written R³ constrained
  Lagrangian with multiplier: `s2_gyroscopic` 3.0e-10, `s2_rotating_frame` 3.2e-10,
  `(q·a)(c·ω)` coupling 3.7e-10; the half-strength (old-rule) force diverges to 0.50.
- **Adversarial Lagrangians** (`TestS2RuleOnAdversarialLagrangians`, 22 systems + 1 with
  `params`): ∂L/∂ω with normal component and q-dependence (`(q·a)(c·ω)`, `(c·ω)(d·ω)`,
  `ω·Aq`, `(Aω)·q`, `½(ω·c)²`), ω and q̇ in one product (3 forms), ξ under `R`, `S`, `Hat(ξ)`
  in δW, two S² variables sharing `c`,`k` (3 forms incl. `ω₁·ω₂` and `ω₁·(q₂×ω₂)`),
  S²×SO(3) (`(Rω)·e₃`, `Ω·(Rq)`, `ω·(Ω×q)`), R³×S² (2), R³×S²×SO(3), R×S², `Hat(ω)` /
  `Hat(Ω)J` / `Hat(Ω)+J` in L. All ≤ 4e-6 at seeds 0–2 (default h) except F3's 1.43e-5;
  all ≤ 1e-6 at h = 1e-4.
- `TS2.delta(substitute=True)` deletion: nothing else referenced it (grep agrees with the report).

### 2. `check_eom`

- S² projection: extension-independent (a term `(|q|²−1)h` contributes `2qh` to the EL
  residual, killed by `−q×`); ω := q×q̇ makes L̃ invariant to the normal part of q̇. Sound.
- SO(3) sign: δL|_R = ⟨A, η̂⟩ = vee(A−Aᵀ)·η with A = Rᵀ∂L/∂R; E = −ṗ + p×Ω + vee(A−Aᵀ) + M;
  matches `_euler_lagrange`. Rigid pendulum (R-dependent) passes; flipping the Coriolis sign
  is detected (`test_numeric`).
- `params`: the ṗ term of `L = ½ω·ω + p(t) c·ω` is produced and verified (6e-7).
- `inputs`/variation zeroing: ξ under `R`, `Rᵀ`, `S`, `Hat` in δW all verified (A8a–d).
- `Scalar(name, value=…)`: engine folds to the literal, evaluator uses `value`; consistent
  (8.6e-7). `UnitNorm`/`Orthogonal`/`Skew`/`Symmetric`/`Identity`/zero resolved the same way
  the engine's rules assume (`vector_rules.py:158`, `simplify.py:133,136,185`). `Ones`-flagged
  non-scalar leaves are not special-cased by the evaluator; the engine only uses `is_ones` for
  `t_diff`, so no disagreement today.
- Terms that vanish at sampled states: states are generic (random q, ω ⊥ q, ω̇ ⊥ q, R, Ω);
  a term parallel to q is invisible by design (F5 for the consequence).
- `trials=5`, default seed: 36 systems × 10 seeds, worst correct score 8.3e-6 (F3).
- Flipped-term detection at h = 1e-4 still ≥ 1.0 (free rigid body, all three terms).

### 3. `collect` peel

- 17 hand shapes with the target in a vector slot (two matrices, `Transpose(R)`, `SMMul` inside
  `MVMul`, both sides `MVMul` with the target on the right at depth 0/1, common `R` stripped,
  `Cross` in `MVMul` in `Cross` on either side, `MVMul(M, Cross)` vs `Cross`): all terminate and
  preserve the value; unexposable shapes (`SVMul`, `VAdd` targets) terminate value-preserving.
- 300 (test) / 2000 (scratch) random Dots: 0 wrong values; the only non-termination is F1.
- Target on both sides: `NonLinearError` (2 shapes).

### 4. `Hat(Vee)` and skew flags

- `SkewSymmMatrix` flag now set (`flags.is_skew_symmetric == True`, `attr` honoured incl.
  `Constant`). `MAdd(K, J)`, `MMMul(K, K₂)`, `R K Rᵀ`, `v vᵀ` are not skew; `Hat(v)`,
  `MAdd(K, K₂)`, `SMMul(K, s)`, `Transpose(K)` are. `Hat(Vee(X))` on 12 shapes: value
  preserved, cancelled exactly for the five skew shapes. `SMMul(s, K)` reorders its operands,
  so `left.is_skew_symmetric` is the matrix. The other leaf classes (`S2`, `SO3`, `TS2`, `TSO3`)
  set `attr` before `super().__init__` or append a flag (`TangentVector`) that `ExprFlags`
  does not parse; no other leaf had the missing-flag bug.

### 5. `extract` changes

- All 36 oracle systems through `compute_eom` → `to_standard_form` (numeric reconstruction
  `M a + f + G u = E` to ≤ 7e-15) and `linearize` + `finite_difference_check` (≤ 5.3e-10); no
  new raise fires (`TestExtractOnPipelineOutput`).
- `Cross` guard order: `x × x → 0` first, then both-sides; `x × (a×x)` and `(a×x) × x` raise
  `NonLinearError` (they are quadratic). `_zero_or_raise` returns zero only for query-free terms
  (`SVMul(a, Dot(b, Mb))`, `Hat(a)`, `(a bᵀ s)ᵀ`); a query-containing but identically-zero
  unsimplified term (`q × (q s)`) raises instead of returning 0 — only reachable by calling
  `extract_linear_coeff` directly on unsimplified input (the pipeline simplifies first).

### 6. Regression surface

- `pytest tests/`: 1441 passed, 5 xfailed (before my file); `ruff` clean.
- `engine_fuzz.py --cases 3000 --seed 4242`: **0 wrong-value findings**; 222 distinct
  informational findings (14 `flag-dependent`, 2 row-vector `ExpressionMismatchError`
  (phase 3.2), 1 intended `NonLinearError`, the rest `unsupported-linear-term`
  `AlgebraicError`s on matrix queries / `Transpose`/`Vee` around the query). `collect[once]`
  2429/2429, `manifold_rules[full delta]` 3000/3000, `extract-nonlinear` 3000/3000.
- Old engine (`git worktree` of 1dd98ab, same oracle) vs new, 36 systems, seeds 0/1/2:

| system | string changed | old error | new error |
|---|---|---|---|
| 21 pre-existing R/R³/SO(3)/product systems | no | ≤ 6.0e-6 | identical |
| quadrotor_with_load | yes (normal term −1→−2) | 2.6e-6 | 2.6e-6 |
| s2_gyroscopic | yes | 0.57–0.70 | ≤ 6.8e-7 |
| s2_anisotropic_omega | yes | 0.29–0.77 | ≤ 5.2e-7 |
| s2_rotating_sphere_omega / _qdot | yes (normal term) | ≤ 7.3e-7 | ≤ 7.3e-7 |
| s2_rotating_frame | yes | 0.27–0.53 | ≤ 5.0e-7 |
| r3_s2_velocity_coupling | yes | 0.50–1.1 | ≤ 2.9e-6 |
| s2_so3_omega_coupling | yes | 0.37–0.61 | ≤ 2.9e-6 |
| s2_so3_axis_coupling | yes | 0.13–0.32 | ≤ 2.7e-6 |
| so3_gyro_coupling_J_on_velocity | yes | 1.6–2.1 | ≤ 2.5e-6 |
| so3_gyro_coupling_J_on_axis | no | ≤ 1.6e-6 | identical |
| so3_inertial_momentum_R_J_Omega / _RJ_Omega | crash → ok | `AlgebraicError` | ≤ 1.2e-6 |
| s2_anisotropic_qdot | crash → ok | `AlgebraicError` | ≤ 9.0e-7 |
| so3_so3_relative_rate_coupling | crash → ok | `AlgebraicError` | ≤ 1.6e-6 |
| so3_offset_mass_matrix_inertia | crash → ok | `AlgebraicError` | ≤ 8.3e-7 |

### 7. Debuggability

- New error messages name the term and the query: `collect` (`'ξ' appears on both sides of
  'Dot(…)'`), `_no_left_coefficient`, `_zero_or_raise` (node type, query, expression), the
  dispatch fallthrough (query/expression types). Good. The one bad message is F1's
  `RecursionError` with no text.
- Docstrings: `manifold_rules` derivation is complete (re-derived above from both starting
  points); `extract` grid matches behaviour; `numeric.py` documents keys, conventions and the
  three residual forms. `collect._collect_dot`'s termination argument is wrong for the matrix-
  factor case (F1). `docs/numeric-verification.md` overstates the detection margin (F2, F3).
- Tests: hand derivations in `test_s2_variation_fix.py` / `test_collect_fix.py` docstrings are
  complete enough to re-derive; I checked the gyroscopic, rotating-frame, S²×SO(3) and
  `(JΩ)·(r×b)` rows by hand. Duplication is F8.

## Required before merge

1. F1: guard the peel on `l.right.has(vec)` and raise a typed error for a target inside a
   matrix factor; fix the docstring and CHANGELOG line; add the shape to `engine_fuzz.py`.
2. F3: default `h=1e-4` (or a larger outer step) and the corresponding doc wording.
3. F2: document the relative-resolution limit and add a per-term power check or the
   resolution in `assert_eom`'s message.

Optional: F4 CHANGELOG sentence, F5 one-line projection note, F6/F7/F8 when phases 2/5 touch
those lines. After 1–3 the seven `xfail`s in `tests/test_phase1_adversarial.py` must be turned
into passing tests (they are strict, so the suite will say which).

---

## Fix log

Branch `fix/phase1-equations`, on top of 50a4518 (uncommitted). Line numbers refer to the
working tree after the change.

| | before (50a4518) | after |
|---|---|---|
| `python -m pytest tests/ -q` | 1529 passed, 12 xfailed | **1550 passed, 7 xfailed** (5 pre-existing in `test_scalar_adversarial.py` + the 2 F2 tests, see below) |
| `ruff format --check` + `ruff check` | clean | clean |
| `mkdocs build --strict` | clean | clean |
| `engine_fuzz.py --quick` | 0 wrong-value; matrix-factor shape not generated | **0 wrong-value, 0 `RecursionError`, 0 timeout**; new campaign `collect[target inside matrix factor]` 51 cases: 39 `AlgebraicError`, 11 `NonLinearError`, 1 ok (see F1) |

### F1 — `collect` loop on a target inside a matrix factor (fixed)

- `geomech/core/math/collect.py:48-66` new `_peel(l, r, vec, expr)`: the peel
  `Dot(M f(v), w) → Dot(f(v), Mᵀ w)` fires only when `vec` is in `l.right`. When `l.left`
  (the matrix factor) holds the target it raises `AlgebraicError` naming the target, the
  matrix factor and the term (`collect: the target 'd/dt(η_R)' occurs inside the matrix
  factor 'd/dt(η_R)(c)'' of the term 'Dot(…)'; transposition exposes a target only in the
  vector operand ((M f(v))·w = f(v)·(Mᵀ w)), so this term cannot be integrated by parts`);
  when both factors hold it, `NonLinearError`. All three peel sites (`(MVMul, Cross)`,
  `(MVMul, MVMul)`, `(MVMul, _)`) go through `_peel`; nothing recurses on that shape.
- `:69-96` docstring: termination argument rewritten. The measure is the depth `d` of the
  target below the Dot on the target-bearing side; peel, triple-product rotation and
  rotation stripping each decrease `d` by one, a flip leaves it unchanged, the peel is
  guarded so a peeled term cannot bounce back, and no flip is followed by a flip (each
  flipped shape is enumerated), hence at most `2d + 1` steps.
- `geomech/core/math/collect.py:13-17` module docstring states the matrix-factor rule.
- `specs/reviews/scripts/engine_fuzz.py:1187-1211` `fuzz_collect` draws the new shape class
  in 25 % of cases: `Dot(MVMul(Hat(v), a), b)` on either side, `v cᵀ`, `a vᵀ`, `M Hat(v)`,
  `M (Hat(v) a)`, with a `Cross` or a second `MVMul` on the other side. `--quick`: 51 cases,
  no `RecursionError`, no timeout, no wrong value. The single "ok" case is
  `Dot((η …ᵀ)…, η)` where the generated right side *is* the target: the pre-existing
  `r == vec` shortcut at `:103` flips it before the both-sides check; the term is
  quadratic in the target, which the pipeline never produces (δL is linear in the
  variation) and which `extract_linear_coeff` rejects. Left as is.
- Tests: `tests/test_collect_fix.py:180-228` `TestTargetInsideMatrixFactor` (4 tests): the
  five direct shapes raise `AlgebraicError` naming the target and the matrix factor;
  `Hat(v) v` raises `NonLinearError`; the pipeline reproducers
  `L = ½ Ω·JΩ + a·((Ω cᵀ) b)` (SO(3)) and `½ ω·ω + a·((ω cᵀ) d) − g q·e₃` (S²) raise
  `AlgebraicError` from `compute_eom` with the IBP target `d/dt(η_R)` / `d/dt(ξ_q)` and the
  outer product `(c)'` in the message. The builders are
  `tests/phase1_systems.py:584-611`. The three F1 `xfail` markers in
  `tests/test_phase1_adversarial.py` (`test_target_inside_matrix_factor_terminates`,
  `TestPipelineReachesCollectLoop` ×2) are removed; the assertions are unchanged and pass.
- `CHANGELOG.md`: the `collect` entry under *Fixed* now states the matrix-factor rule and
  the loop it replaces; the "no longer recurses forever on `Dot(M a, N b)`" line is
  restricted to the vector-operand case it is true for.

### F3 — default step (fixed)

- `geomech/numeric.py:565-572` `DEFAULT_H = 1e-4` with the tradeoff comment (nested
  central differences: roundoff ≈ 1e-16 |L| / h², truncation ≈ h²; ≈ 1e-6 at 1e-5, ≈ 1e-8
  at 1e-4); `check_eom` (`:740`) and `assert_eom` (`:814`) default to it; `h` documented
  in `check_eom`'s Args. `docs/numeric-verification.md` "Anything below finite-difference
  noise" bullet rewritten with the same numbers; "~1e-6" → "~1e-7" throughout.
- The two F3 `xfail` markers are removed; assertions unchanged and passing.
- 36 oracle systems (`tests/test_eom_oracle.SYSTEMS`), seeds 0–9, `trials=5` each
  (scratch script, 31 s per step):

  | `h` | worst correct score | worst flipped-term score |
  |---|---|---|
  | 1e-4 (new default) | **2.34e-7** (`spherical_pendulum_qdot`) | **2.14e-2** (`quadrotor_with_load`, term `(ω×q)×ẋ_L m_Q l`, seed 7) |
  | 1e-5 (old default) | 8.33e-6 (`s2_anisotropic_qdot`) | 2.14e-2 (same) |

  The flipped-term score is the minimum over every additive term of every equation of
  `check_eom` with that term's sign flipped (`eq − 2·term`), excluding the 12 terms that
  are parallel to `q` or zero on the constraint set (`q × (ω × (e₃ × q))`,
  `(Ω × ρ)·ρ`, `Q (Ω × Ω)`; invisible by design). Flipped scores do not depend on `h`
  (they are set by the term's size), correct scores drop 36× at the new default.

### F2 — resolution of small terms next to large ones (implemented; two xfails stay)

- `geomech/numeric.py:613-698`: `_scale`, `_tangent_norm`, `_terms`, `_ZERO_TERM_RTOL`
  and `TermResolution` (public). For each additive term `t` of each equation the oracle
  records, over the trials, the best `ratio = 2 |P t| / max(1, |P E_true|)`: the relative
  error the check would report if the term had the wrong sign. A term with
  `|P t| ≤ 1e-10 · max(|t|, scale)` is parallel to `q` or zero to within roundoff and has
  nothing to resolve (`tangent = False`). `_eom_errors` (`:700-735`) returns
  `(worst, terms)`.
- `check_eom` (`:738-800`) still returns the relative error only; its docstring has a
  "Resolution" paragraph that states the limit, the grey zone and that `assert_eom` does
  the per-term check.
- `assert_eom` (`:803-861`) gains `require_resolved: bool = True`: after the relative-error
  gate it raises `AssertionError` listing every term that is unresolved at `rtol` in every
  trial, with its tangent size and the resolution `rtol · max(1, |E_true|)`. On the review's
  reproducer (`k = 1e7`) the message names all three order-one terms (mass, gravity,
  forcing) as unresolved at resolution 64; with `require_resolved=False` it returns 3.5e-9
  as before. Every existing `assert_eom` call in the suite passes with the default (the
  only false alarm during development was the roundoff-zero term
  `q × (ω × (e₃ × q))` of the rotating sphere, which is why the zero-term floor is
  relative to the equation's scale and not only to the term's own size).
- `docs/numeric-verification.md`: "so there is no grey zone" and "Wrong equations are off
  by order one" removed; new section "Resolution: what a term must be for the check to
  see it" with the `k = 1e7` numbers (correct ~1e-9, flipped ~1e-7, both accepted by
  `check_eom`; absolute noise 6e-3 to 2e-2) and the `assert_eom` behaviour; the
  "below finite-difference noise" bullet states that the grey zone always exists and widens
  with `|L|`; step 3/4 of "Using it on a new system" mention the unresolved-term report.
  `CHANGELOG.md` *Added* entry for `geomech.numeric` updated accordingly.
- **The two F2 `xfail` tests stay `xfail` (strict, still failing), unchanged.** Their
  assertions are `check_eom(correct) < 1e-5` *and* `check_eom(flipped) > 1e-2` (resp.
  `check_eom(dropped) > 1e-2`) from the same call with default arguments. Measured at
  `h = 1e-4` on that system, per trajectory point: absolute residual of the *correct*
  equation 5.8e-3 / 1.8e-2 / 9.5e-3 (seeds 0/1/2; it is the roundoff 1e-16 |L| / h² with
  |L| ≈ 1e7), flipped-gravity residual 3.7 / 0.027 / 0.85 (m·g is 0.003 at seed 1).
  The two assertions need one normalisation `S` with `noise/S < 1e-5` and
  `signal/S > 1e-2`, i.e. `S > 1800` on the noisy trials and `S < 370` on the trial that
  carries the signal; no function of the residual and the equation does that, and
  per-term normalisation makes the *correct* equation fail instead (noise 6e-3 over a
  term of size 1). Reducing the noise enough (≥ 3 orders) would need a different
  differentiation scheme (polynomial-exact stencils with a larger step, or exact jets),
  which contradicts the prescribed `h = 1e-4` central differences and is beyond this
  fix. So the assertion, not the code, is what cannot be met: reported here, not edited.
  The module docstring of `tests/test_phase1_adversarial.py` says so.

### F4 — `quadrotor_with_load` normal-component change (documented)

- `specs/reviews/phase1-engine-report.md` "Pinned equation strings in the existing suite":
  the load equation's term parallel to `q`, `Cross(ω, Cross(q, ẋ_L)) m_Q l`, changed
  coefficient −1 → −2; it is `q (ω·ẋ_L)`, invisible on the tangent plane.
- `CHANGELOG.md` S² rule entry: same sentence.

### F5 — project S² equations before solving (documented)

- `docs/numeric-verification.md` "The normal component on S²" bullet: `E = 0` is only
  meaningful as `(I − qqᵀ)E = 0`; solving for `ω̇` from `to_standard_form` / `linearize`
  blocks must be done as `(I − qqᵀ)(M ω̇ + f + G u) = 0` with `ω̇ ⊥ q`.
- `geomech/dynamics/standard_form.py:60-70` `to_standard_form` docstring, *Notes*: same.

### F6 — `attr` list aliasing (fixed)

- `geomech/core/base/expressions.py:195` (`Scalar`), `:308` (`Vector`), `:338` (`S2`),
  `:455` (`Matrix`), `:465` (`SkewSymmMatrix`), `:474` (`SO3`): every leaf copies the list
  it is given before storing or appending to it. `TS2`/`TSO3` already built their own list.
- `tests/test_matrices.py:189-230` `TestAttrListIsNotShared` (8 tests): a shared
  `["Constant"]` list is unchanged after building each of the six leaf classes from it, a
  `Matrix` built from the same list afterwards has none of their flags, a skew matrix, a
  rotation and an `S2` point built from one list keep their own flags, and `TS2`/`TSO3`
  do not alias.

### F7 — empty n-ary nodes (forbidden)

- Chosen: **forbid**. `geomech/core/operations/addition.py:31-57` `_flatten_nodes` raises
  `TypeError("Add() needs at least one operand; use the zero leaf of the type for an empty
  sum")` for `Add()`, `VAdd()`, `MAdd()`, `X(*[])` and `X([])`; the docstring says why. The
  package never constructs one: every `Add(*list)` site in `extract.py`, `simplify.py`,
  `ibp.py`, `pull.py`, `vector_rules.py` and `collect.py` is preceded by an emptiness check
  that returns the zero leaf, the full suite and the fuzz run pass with the guard, and the
  only empty construction anywhere was an unused `add = Add()` line in
  `examples/nary_check.py`, removed (the script still runs).
- `tests/test_addition.py:506-523` `TestEmptyAddition` (4 tests).

### F8 — duplication (done for the four files named)

- `tests/phase1_systems.py` is the single home of: `HALF` (imported from
  `tests/vbl_systems`, no second definition), `TOL`, `derive`, `equation` (`:36-49`),
  `s2_symbols` and the 22 `adv_*` builders with `ADVERSARIAL_SYSTEMS` (`:291-557`),
  `stiff_point_mass` (`:565`), `so3_outer_product_on_velocity`,
  `s2_outer_product_on_velocity` (`:584-611`).
- `tests/test_s2_variation_fix.py` and `tests/test_collect_fix.py` import `derive`,
  `equation` instead of defining `_derive`, `_equation`, `TOL`; `tests/test_phase1_adversarial.py`
  imports `HALF`, `TOL`, `s2_symbols`, `ADVERSARIAL_SYSTEMS`, the two named builders it
  uses and builds `_stiff_point_mass` from `ps.stiff_point_mass` (372 lines shorter). No
  assertion changed. Out of scope and left: `HALF` in `test_eom_oracle.py`,
  `test_numeric.py`, `test_scalar_adversarial.py`, `test_vbl_adversarial.py` and the 12
  systems `test_eom_oracle.py` re-derives from the pipeline tests (reviewer's
  `tests/systems.py` suggestion).

### Test counts

`tests/test_phase1_adversarial.py`: 93 passed, 2 xfailed (was 88 + 7). New tests: 4
(`test_collect_fix.py`), 8 (`test_matrices.py`), 4 (`test_addition.py`). Suite: 1550
passed, 7 xfailed; 16 files changed, +886 / −495.

### Arbitration of the F2 tests (orchestrator)

The two F2 tests asked the oracle to accept the correct stiff equation (k = 1e7) and to
detect a flipped or dropped gravity term of size ~1e-2 in the same default call. The fix
log shows this is below the roundoff floor of nested central differences at that scale
(absolute noise 6e-3 to 2e-2), so the assertion was unsatisfiable by any normalisation.
The tests were rewritten to assert the honest behaviour instead: `assert_eom` raises and
names the gravity term as unresolved, `require_resolved=False` accepts the equation on
the relative error alone, and with order-one parameters the dropped term is detected.
Suite after arbitration: 1552 passed, 5 xfailed (the five pre-existing ones).
