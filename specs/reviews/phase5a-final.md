# Phase 5a final review: simplifier completeness

Branch `fix/phase5a-simplifier`, range `47cce2d..36122e1` (implementation cbd3653, review
cb4697b, fixes 36122e1). Read-only on `geomech/` and `tests/`; nothing committed. Every
number below comes from a run made for this review. Base comparisons used a temporary
git worktree of 47cce2d, removed afterwards.

## Verdict

**GO WITH CHANGES.** I found no wrong equation and no pinned string whose new form differs
in value from the old one. All gates pass. The docs reproduce exactly, and timing is
within budget. Three changes are needed before merge, and none of them touches engine
code:

- one decision record is missing;
- one convention is not stated;
- eight stale lines in test comments still spell out the old printed forms.

## Required changes

1. **D1 is not recorded in the fix plan.** `specs/reviews/phase5a-review.md:242` says
   D1 "was approved by the author on 2026-11-13 (option "Stay symbolic"), recorded in the
   Decisions section of `specs/fix-plan.md`". At 36122e1, and in the working tree,
   `specs/fix-plan.md:174` still holds only the options and the recommendation. No commit in
   the range touches the file. The orchestrator should confirm the approval with the
   author, then append to the D1 row at `specs/fix-plan.md:174`:
   `**Decided (a), 2026-11-13: named scalars with a value stay symbolic (phase 5a).**`
   If the approval cannot be confirmed, D1 has to go back to the author before merge,
   because it is a semantic change listed under "not changed without a yes".

2. **`docs/conventions.md:76–80` does not state the non-commutativity safeguards of
   like-term merging.** After "Addends with different symbolic coefficients are not
   factored (`x a + x b` stays).", add:
   `Only scalar products and Dot are treated as commutative: bases and factors built with
   MMMul, MVMul, VVMul, Cross, SVMul or SMMul are compared operand by operand, so
   A B v s − B A v s, (a×b) m + (b×a) m and aᵀ b s − bᵀ a s are not merged (the last is
   equal but not recognised).`
   I checked all three examples at HEAD. `full_simplify(A*(B*v)*s - B*(A*v)*s)` gives
   `(ABvs+BAvs(-1))`, and `a.T()*b*s - b.T()*a*s` gives `((a)'bs+(b)'as(-1))`, while
   `Dot(a,b)s − Dot(b,a)s` gives `0`.

3. **Comment lines that transcribe the printed string still show the old order.** Each
   pinned assert is correct, but the comment above it no longer matches. Six months from
   now a reader comparing the two would see a different factor order. These are comment
   and docstring edits only; no assertion changes.
   - `tests/test_manifold_pipelines.py:242` →
     `"""Full EOM string: d/dt(ω)*m*l*l*(-1) + Cross(q,e3)*m*g*l*(-1) + f"""`
   - `tests/test_manifold_pipelines.py:247` → `"""M[ω̇] = I*m*l*l*(-1)  (= -ml² I)"""`
   - `tests/test_manifold_pipelines.py:281` →
     `#   Cross(q, Cross(d/dt(ω), q))*m*l*l*(-1) + Cross(q,e3)*m*g*l*(-1) + f`
   - `tests/test_manifold_pipelines.py:288` →
     `#   M[ω̇] = Hat(q)*Hat(q)*m*l*l = ml²*(qq^T - I) = -ml²*(I - qq^T)`
   - `tests/test_manifold_pipelines.py:305` →
     `"""EOM: Cross(q, Cross(ω̇, q))*m*l*l*(-1) + Cross(q,e3)*mgl*(-1) + f"""`
   - `tests/test_manifold_pipelines.py:313` →
     `"""M[ω̇] = Hat(q)*Hat(q)*m*l*l  (= ml²(qq^T - I), the S2 projection)."""`
   - `tests/test_scalar_adversarial.py:190`: replace ``I(-1)mll`` with ``Imll(-1)``.
   - `tests/test_scalar_adversarial.py:409`: replace `d/dt(ω)(-1)mll` with `d/dt(ω)mll(-1)`.

## Verified

### 1. Gates

| check | command | result |
|---|---|---|
| suite `-W error` | `venv/bin/python -m pytest tests/ -q -W error` | 2238 passed, 2 xfailed |
| suite `GEOMECH_DEBUG=1` | `GEOMECH_DEBUG=1 venv/bin/python -m pytest tests/ -q` | 2238 passed, 2 xfailed |
| suite `GEOMECH_DEBUG=1 -W error` | same with `-W error` | 2238 passed, 2 xfailed |
| ruff | `ruff format --check` and `ruff check` on `geomech/ tests/ specs/reviews/scripts/` | 100 files formatted, all checks passed |
| docs | `mkdocs build --strict` | builds (one INFO about `reference/`, not a warning) |
| fuzz | `engine_fuzz.py --cases 2000 --seed 6060`, run from `/tmp` | first line `fuzzing geomech from <this checkout>/geomech`; 0 consistent-mode findings; every value-preserving transform 2000/2000 ok; convergence 2000/2000; 188 distinct findings (107 flag-dependent, 75 unsupported-linear-term, 6 documented-error-on-legal-input) |
| fuzz, base | HEAD's script copied into the 47cce2d worktree, same seed | header names the worktree; 0 consistent-mode findings; 137 distinct (54 flag-dependent, 75, 6, 2 semi-documented-exception). The difference is +53 flag-dependent and −2 manifold-leaf exceptions (F1), as the review predicted |
| examples | each `examples/*.py` | 6/6 run |
| notebooks | `jupyter nbconvert --execute` into scratch, then cell outputs compared with the committed ones (execution counts and metadata ignored) | 4/4 run; 0 differing cells (7, 8, 8, 7 code cells) |

Running the fuzz rewrites `specs/reviews/scripts/engine_fuzz_findings.json`, which is
git-ignored (`.gitignore:148`). The working tree stays clean.

### 2. Spot check of `phase5a-string-changes.md` (independent harness)

A pytest plugin wrapped `compute_eom`, `to_standard_form`, `linearize` and
`integrate_by_parts` at top level. It ran the same 11 test ids in a 47cce2d worktree and
at HEAD (13 tests each tree, all passing). Because the base tests assert the old strings
and the HEAD tests the new ones, each old string reproduces at base and each new string at
HEAD. Each output was then evaluated at 3 states:

- configuration from `sample_state`, so |q| = 1, ω ⊥ q, ω̇ ⊥ q, R ∈ SO(3);
- every other leaf drawn from a generator seeded by its printed name and the state seed,
  so the values do not depend on tree shape.

Rows checked:

1. `TestFreeRigidBodyPipeline` (equation and M)
2. `TestSphericalPendulumQdotPipeline::test_M` (`llm` → `mll`)
3. `TestRigidPendulumCOMPipeline` (equation, M, f)
4. `TestDoubleRigidPendulumMassMatrixExact` (both equations, both M blocks, both f)
5. `TestRelativeRateCouplingSO3xSO3` (merges and Q(Ω×Ω) removal)
6. `TestAnisotropicInertiaQdotForm` (a term changes to `Cross(q,Aq)Dot(ω,ω)`)
7. `TestRotatingSphere` (both forms)
8. `TestVaryingInertia` (merge to `(-0.5)`)
9. `TestElasticSphericalPendulumQdotForm` (cancelled pair)
10. `test_vbl::TestQuadrotor::test_translational_row_exact` (an A block)
11. `test_ibp_prefactor::test_repeated_factor_is_one_term`

Result: 31 changed outputs, worst relative difference **1.85e-16**. The S² equations
were compared raw, not projected, so they agree in the normal direction too. No block
appeared in one tree only.

### 3. Lee et al. systems at HEAD

Printed with `compute_eom` on the `tests/vbl_systems.py` builders, in both trees:

| system | HEAD equation | vs book / pinned hand derivation |
|---|---|---|
| free rigid body (Lee et al. Eq. 6.16) | `(J Ω̇(-1) + Cross(Ω,JΩ)(-1) + M)` | `−JΩ̇ − Ω×JΩ + M = 0`, term for term; only `(-1.0)` → `(-1)` |
| rigid pendulum, pivot | `(JΩ̇(-1) + Cross(Ω,JΩ)(-1) + Cross(ρ,Rᵀe3)mg(-1) + M)` | matches the stated `−JΩ̇ − Ω×JΩ − mg ρ×Rᵀe₃ + M` term for term (sign convention of e₃ unchanged from base) |
| spherical pendulum (Lee et al. Eq. 5.14) | `(ω̇ mll(-1) + Cross(q,e3)mgl(-1) + f)` | `−ml²ω̇ − mgl q×e₃ + f`; number last |
| quadrotor | `x`: `(ẍ m(-1) + e3 mg(-1) + Re3 f)`; `η` as the free rigid body | `−mẍ − mge₃ + fRe₃`; matches the `TestQuadrotor` docstring |
| double rigid pendulum | `η₁` 14 terms, `η₂` 11 terms, same terms as base with numbers last; `M[Ω̇₁] = (J1(-1)+Hat(ρ₁)Hat(ρ₁)m1+Hat(l₁)Hat(l₁)m2)` | matches the `TestDoubleRigidPendulumMassMatrixExact` hand derivation `−J₁ + m₁ρ̂₁ρ̂₁ + m₂l̂₁l̂₁` symbol for symbol |

The docstrings that state a book form in mathematical notation are still readable
against the new strings. The phase-5a notes added in `TestRelativeRateCouplingSO3xSO3`
and `TestRotatingSphere` explain the merged terms. One case needs a moment's thought:
E₁'s "+kΩ₁ × QΩ₂" prints as `Cross(QΩ₂, Ω₁)k(-1)`, which equals it after flipping the
cross product. The only lines a reader could not match without reordering are the eight
transcriptions in required change 3.

At HEAD every one of these equations is a fixed point of `full_simplify`. `to_latex`
renders the number-last form with the sign in front (`-m l l \dot{\omega_{q}}`).

One cosmetic change is not listed in the string table: the quadrotor-with-load
ξ equation's input term changes from `Cross(q,Re3)fl(-1)` to `Cross(q,Re3)lf(-1)`. It is
not pinned, the same term with its factors reordered.

### 4. Docs

- `docs/pipeline.md` summary block: rerunning `tr.summary(chars=70)` on
  `rigid_pendulum_pivot` reproduces the block byte for byte.
- `docs/pipeline.md` rule list: with `set_trace(True)`, 42 `geomech.rules` lines, and
  every rule name and count in the list matches. That includes `pull.mvmul_scaled` (3),
  `pull.svmul_nested` (3), `simplify.mul_fold` (7), `collect.flip_mvmul_right` (1) and
  `pull.transpose_scaled` (1).
- `docs/troubleshooting.md`: both quoted messages reproduce exactly. They are the
  `Dot(w,Vee(Hat(η̇)A))` term and the prescribed-parameter term
  `\frac{d}{dt}(\omega_{q})m\delta{l}l(-2)`; the `(-2)` comes from merging two like terms.
- `docs/vbl.md`: the `print(lin)` block and `A[δω] = (ω m l (2))'` reproduce.
- `docs/conventions.md` describes:
  - number-last folding (`:72–73`);
  - integer-float normalisation (`:73–75`);
  - D1 and `is_number_literal` (`:85–92`);
  - the multiset like-term rule (`:76–80`).
  The non-commutativity safeguard is missing; see required change 2.
- The diff contains no Scala references, no absolute paths and no initials for the book.

### 5. Timing (median of 20 runs after a warm-up, 3 repetitions, two alternating rounds; ms)

| system | call | 47cce2d | HEAD | ratio |
|---|---|---|---|---|
| quadrotor with load | `compute_eom` | 10.8–11.0 | 11.8–11.9 | 1.08× |
| quadrotor with load | `linearize` | 6.2–6.3 | 6.8–7.0 | 1.10× |
| double rigid pendulum | `compute_eom` | 17.3–17.5 | 20.0–20.4 | 1.16× |
| double rigid pendulum | `linearize` | 17.7–18.0 | 20.8–21.1 | 1.17× |

All are well inside the 2× budget. A profile of `compute_eom` on the double pendulum puts
`_eliminate`, `_commutative_equal` (like-term grouping) and `map_children` at the top.
The function-local imports in `rebuild` do not show up in the profile.

## Advisory

- `docs/conventions.md:73–75` and the CHANGELOG say a float with an integer value
  "becomes an int". This does not hold for a bare `Number(2.0)` leaf, which the fix log
  scopes out. Either add "(in a product or sum)" or normalise the leaf in 5b.
- Literal spelling is not normalised. `x·Scalar("0.5", value=0.5)` prints `x0.5` and
  `x·Number(0.5)` prints `x(0.5)`: same value, different `str`. They do merge as like
  terms (`x·0.5 + x·(0.5)` → `x`). The `·0.5` in the `pipeline.md` summary comes from
  this, because `_simplify_mul` keeps a chain already in normal form as built
  (`simplification.py:732`).
- `tests/test_scalar_adversarial.py:328` is still named
  `test_f_is_mathematically_correct_but_uncombined`, but since phase 5a it asserts the
  merged form. It now duplicates `test_f_combines_like_terms` (`:334`).

## Phase 5b follow-ups

1. **Inner-product commutativity in like terms.** `simplification.py:591` compares
   `VVMul` operand by operand even when it is a scalar inner product, so `aᵀb s − bᵀa s`
   is kept. The value is correct; it is just not reduced. This is the D2/5.5 territory.
2. **Like-term grouping cost.** `_collect_like_terms` (`simplification.py:650`) is
   O(n²) in `_commutative_equal`, which is the second-highest self-time in the profile.
   A canonical key per addend would make it O(n) and remove most of the 1.16× slowdown.
3. **Constants module (5.5).**
   - The literal `-1` is built inline at `simplification.py:339`, `:346` and
     `vector_identities.py:91`.
   - `_ZERO`/`_ONE` lambdas sit at `vector_identities.py:148–149`.
   - The `HALF`/`Scalar("0.5", value=0.5)` spelling is a separate literal (advisory 2).
   - Bare `Number(2.0)` is not normalised (advisory 1).
4. **Dead code (5.5).**
   - `vector.transpose_transpose` (`vector_identities.py:34`, `:210`) is unreachable,
     because the constructor collapses double transposes.
   - `is_leaf` and `has_nested_add` (`transformations/utils.py:14`, `:22`) are unused.
5. **Symbolic coefficients are not factored** (`ẍ_L m_Q(-1) + ẍ_L m_L(-1)` in the
   quadrotor with load). This was a deliberate choice in 5a. If listings for the paper
   matter, decide whether `collect_terms` should factor a common base.
6. **`q × (ω̇ × q)` is kept** (`vector_identities.py` docstring; `troubleshooting.md:17`).
   Removing it needs the constraint system (`ω̇ ⊥ q`), which belongs to a later phase.
7. **S² linearize blocks are not unique** (CHANGELOG "Changed"). If raw-block identity
   is ever pinned, project the δω column with `(I − qqᵀ)` in `linearize` instead of
   relying on equation form.
8. **Function-local imports in hot helpers.** These are `transformations/utils.py:41–42`
   and `:62–64`, and `utils/printing.py:117`, `:167`, `:508`. They are not significant
   today. Hoist them when 5.6 splits these modules, if the import cycle allows.
