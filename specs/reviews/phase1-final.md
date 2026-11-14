# Final review of phase 1 (`specs/fix-plan.md` 1.1–1.6)

Branch `fix/phase1-equations`, range `1dd98ab..8e28578` (636e963 oracle, c96d31c engine
fixes, 50a4518 adversarial review, 8e28578 fixes F1–F8). Read: the contributor notes, the plan, the
three phase-1 reports including the fix log and the arbitration note, and
`git diff 1dd98ab..8e28578 -- geomech/ docs/ CHANGELOG.md README.md` in full, plus the
changed existing tests and the new test files. Nothing under `geomech/` or `tests/` was
modified; nothing was committed. Scratch scripts are in a temporary directory.

## Verdict: GO WITH CHANGES

The mathematics is right and independently confirmed: the S² rule re-derives to
`δω = ξ̇ − 2 ω × ξ` from first principles and from Lee et al. (5.13) (checked against the
book text, p. 216), all ten dynamics-audit rows rebuilt from the public API pass
`check_eom` at ≤ 2.7e-6 with two seeds and two steps, every hand derivation I checked in the
pinned tests is complete and reproduces the pinned string, and the engine diff has no
remaining hole I could find. Four small things must change before merge, none in
`geomech/`: the fuzz gate script crashes at the seed the plan's "Done when" is run with, a
shipped notebook still teaches the old formula, the paper fix is sitting uncommitted in the
accompanying paper's working tree, and the adversarial test file's docstring describes a state (xfail) that
no longer exists.

## Required changes

1. **`specs/reviews/scripts/engine_fuzz.py:1195` crashes at seed 777 / 2000 cases**, so
   the plan's gate ("`engine_fuzz.py` reports 0 wrong-value cases") cannot be run as asked:
   `python specs/reviews/scripts/engine_fuzz.py --cases 2000 --seed 777` dies in the new
   `collect[target inside matrix factor]` campaign with
   `ExpressionMismatchError: VVMul cannot be performed between VECTOR and VECTOR`
   (`multiplication.py:207`) because `gen.vector(d)` can return a `Transpose` (the row-vector
   encoding) and `VVMul(Transpose(..), Transpose(vec))` is refused. `--quick` (default seed)
   happens not to draw that shape. Fix in the builder list, exact text:

   ```python
           def _col(v):  # VVMul(v, wᵀ) needs a column vector on the left
               while isinstance(v, Transpose):
                   v = v.expr
               return v

           matrix_factor_builders = [
               ...
               lambda: Dot(gen.vector(d), MVMul(VVMul(_col(gen.vector(d)), Transpose(vec)), gen.vector(d))),
   ```

   With exactly this patch applied to a scratch copy, `--cases 2000 --seed 777` completes:
   **0 wrong-value findings, 0 `RecursionError`, 0 timeouts** (details under "Verified").

2. **`notebooks/spherical_pendulum.ipynb:184`** still reads
   `"3. Apply manifold rules: δ(ω) = ξ̇ − ω × ξ\n",` — the old rule, in a shipped notebook that
   mirrors the example. Replace with `"3. Apply manifold rules: δ(ω) = ξ̇ − 2 ω × ξ\n",`
   (the notebook's derived equation is a pendulum, `∂L/∂ω ⊥ q`, so its stored outputs are
   unaffected; no other notebook or example states the formula — grepped).

3. **The paper fix is not committed.** The accompanying paper's background and toolbox
   sections do read `\dot{\xi} - 2\,\omega \times \xi`, and the surrounding
   sentence is mathematically right (see "Verified", item 4), but `git status` in
   the paper's repository shows both files **modified and uncommitted** (together with
   three other files). The engine report's "paper already
   carries −2; no change needed" was read off the dirty working tree. "Done when … paper
   formula fixed" needs the two files committed and pushed to the paper's remote (the
   other three modified files are outside this review; commit them or not as you intend).

4. **`tests/test_phase1_adversarial.py` wording** (three lines, no assertion changes):
   - `:4-11` still opens with "The confirmed findings were strict ``xfail`` tests …" although
     no `xfail` remains (`grep xfail` finds only this sentence). Proposed replacement for
     lines 4–11:

     > All findings of the review are now ordinary passing tests (the history is in
     > `specs/reviews/phase1-review.md`). F1 (`collect` loop) and F3 (default step) were
     > fixed in the engine and the oracle. F2 is a resolution limit of the oracle, not a
     > bug: the relative error is normalised by `max(1, |E_true|)`, so next to a stiff
     > `k = 1e7` term a flipped order-one gravity term moves it by ~5e-7 at any step, and
     > the absolute roundoff of the residual (`1e-16 |L| / h²`, 6e-3 to 2e-2 at `h = 1e-4`)
     > is comparable to the term itself; `assert_eom` therefore reports such terms as
     > unresolved instead of passing silently, and the F2 tests assert that.
   - `:595` "while m g is of order 1e-2" is wrong for the default seeds (`m`, `g` are unit
     normals; `2 m g ≈ 3.7` at seed 0, `0.003` at seed 1). Proposed: "while the gravity
     term is of order one (its tangent size is 1.8 at seed 0)".
   - `:625` comment `# "no grey zone" (docs) needs a margin below the gate` cites a docs
     sentence that F2 removed. Proposed: `# well below the 1e-5 gate; the docs promise ~1e-7
     for order-one parameters`.

## Verified

Commands run from `pygeomech/` with `source venv/bin/activate` (Python 3.12.13):

| check | command | result |
|---|---|---|
| suite | `python -m pytest tests/ -q` | **1552 passed, 5 xfailed** in 5.9 s (the 5 are the pre-existing ones in `test_scalar_adversarial.py`) |
| lint | `ruff format --check geomech/ tests/`; `ruff check geomech/ tests/` | 73 files already formatted; all checks passed |
| docs | `mkdocs build --strict` | built in 0.9 s; only the pre-existing INFO that `manifold-implementation-notes.md` is not in `nav` |
| fuzz, as asked | `python specs/reviews/scripts/engine_fuzz.py --cases 2000 --seed 777` | **crashes** (required change 1) |
| fuzz, default seed | `… --quick` | 0 wrong-value; `collect[target inside matrix factor]` 51 cases: 39 `AlgebraicError`, 11 `NonLinearError`, 1 ok; `manifold_rules[full delta]` 200/200 |
| fuzz, patched copy | scratch copy with the one-line fix of change 1, `--cases 2000 --seed 777` | **0 wrong-value** (162 distinct findings: 86 `unsupported-linear-term`, 69 `flag-dependent`, 7 `documented-error-on-legal-input`; 0 `RecursionError`, 0 timeout); `collect[target occurs once]` 1192/1192, `[inside matrix factor]` 501 cases → 428 `AlgebraicError` + 72 `NonLinearError` + 1 ok, `manifold_rules[full delta]` 2000/2000, `extract-nonlinear` 2000/2000 |
| import | `cd /tmp && python -I -c "from geomech import *"` | works from a foreign cwd; `Evaluator`, `evaluate`, `check_eom`, `assert_eom`, `finite_difference_check` exported; `geomech` has no `__all__` (pre-existing) |
| docs example | the 15-line `check_eom` example and the `evaluate` example of `docs/numeric-verification.md`, run verbatim | 2.3e-8 and `2.0` as the page says |

### 1. The S² rule, re-derived

Take ω := q × q̇ (so ω ⊥ q and q̇ = ω × q) and δq = ξ × q with ξ ⊥ q along the whole
family of curves. Vary ω directly: δω = δq × q̇ + q × δq̇ = (ξ × q) × q̇ + q × (ξ̇ × q) +
q × (ξ × q̇). With |q| = 1, q·q̇ = 0 and q·ξ = 0 the three BAC-CAB expansions are
(ξ × q) × q̇ = q (ξ·q̇), q × (ξ̇ × q) = ξ̇ − q (q·ξ̇), q × (ξ × q̇) = 0, so
δω = ξ̇ + q (ξ·q̇ − q·ξ̇). Because ξ·q ≡ 0 along the curve, ξ̇·q + ξ·q̇ = 0, hence
δω = ξ̇ − 2 q (q·ξ̇). Finally q·ξ̇ = −ξ·q̇ = −ξ·(ω × q) = q·(ω × ξ), and ω × ξ is parallel to
q because both factors are tangent, so q (q·ξ̇) = ω × ξ and δω = ξ̇ − 2 ω × ξ. The term
the old rule dropped, −q (q·ξ̇) = −ω × ξ, is the normal component that keeps δ(q·ω) = 0;
it is exactly the intermediate line of Lee et al. (5.13), `δω = (q̇·γ)q + γ̇ − (q·γ̇)q`
(book p. 215–216, checked in the PDF), after the same substitution. The module docstring
of `manifold_rules.py` and `docs/manifold-implementation-notes.md` carry this derivation
correctly, and the identities they state (`q·δω = −q·ξ̇`, the old rule being
`(I − qqᵀ)δω`) check out.

### 2. The dynamics audit's rows, rebuilt from the public API

A script in a temporary directory builds the six F1 rows and the four F3 rows of
`specs/reviews/dynamics-audit.md` (plus the `((RJ)Ω)·e₃` grouping and `½(Ω×ρ)·J(Ω×ρ)`)
without importing any test builder, runs `compute_eom` and `check_eom` with seeds 0 and 1
at `h = 1e-4` (default) and `h = 1e-5`, and `assert_eom` with the resolution check:

| row | h=1e-4, seeds 0/1 | h=1e-5, seeds 0/1 |
|---|---|---|
| F1 anisotropic `½ ω·Aω` | 2.6e-8 / 8.7e-9 | 5.2e-7 / 3.4e-7 |
| F1 rotating frame | 8.3e-9 / 1.0e-8 | 3.4e-7 / 4.8e-7 |
| F1 gyroscopic `c·ω` | 5.9e-9 / 7.0e-9 | 6.8e-7 / 3.2e-7 |
| F1 R³×S² `ẋ·ω` | 8.9e-9 / 9.3e-9 | 5.9e-7 / 1.7e-6 |
| F1 S²×SO(3) `k ω·Ω` | 5.3e-9 / 3.0e-8 | 7.5e-7 / 1.4e-6 |
| F1 S²×SO(3) `k ω·(Rᵀe₃)` | 3.8e-8 / 1.1e-8 | 2.7e-6 / 2.4e-6 |
| F3 `(R(JΩ))·e₃` | 6.9e-9 / 4.5e-9 | 1.0e-6 / 6.0e-7 |
| F3 `((RJ)Ω)·e₃` | 9.0e-9 / 4.5e-9 | 6.5e-7 / 7.6e-7 |
| F3 S² `½ q̇·A q̇` | 4.0e-8 / 7.8e-8 | 7.2e-7 / 5.3e-7 |
| F3 S² `½ m q̇·q̇ + k (A q̇)·(q × e₃)` (not pinned by the engine report) | 2.4e-8 / 4.0e-8 | 1.9e-6 / 1.5e-6 |
| F3 SO(3)×SO(3) `k Ω₁·(R₁ᵀ(R₂Ω₂))` | 2.2e-8 / 2.2e-8 | 1.6e-6 / 1.2e-6 |
| `½ (Ω×ρ)·J(Ω×ρ)` | 3.0e-9 / 1.2e-8 | 6.0e-7 / 8.3e-7 |

All ≤ 1e-5; worst 2.7e-6 (at the old step), 7.8e-8 at the default; every row "resolved".
No crash. The printed gyroscopic equation is
`d/dt(ω)(-1)mll + Cross(ω,c)(2) + Cross(q,e3)mgl(-1) + f`, i.e. the audit's expected `(2)`.

Additionally, over the 36 systems of `tests/test_eom_oracle.SYSTEMS` and seeds 0–9 the
worst correct score is **2.34e-7** (`spherical_pendulum_qdot`, seed 3), matching the
CHANGELOG, and `assert_eom` with the resolution check passes on all 36. Flipping every
additive term of every equation (seed 0): 233 visible terms all score ≥ **0.15**; the 12
terms the resolution logic classifies as invisible are the ones parallel to `q` or
identically zero (`ω×(q×ẋ_L)`, `(Ω×ρ)·ρ`, `Q(Ω×Ω)`, …), as documented.

### 3. Hand derivations in the pinned tests

Checked line by line against the pinned strings, using `E ≡ −ṗ + 2 ω×p + q×∂L/∂q + f (mod q)`
on S² and `E = −ṗ + p×Ω + (δR-terms) + M` on SO(3):

- `test_s2_variation_fix.py`: gyroscopic, anisotropic ω-form, rotating frame, R³×S²,
  both S²×SO(3) couplings (6 of 7) — every intermediate (`ṗ`, `2ω×p`, `q×∂L/∂q`, the
  `d/dt(Rᵀe₃) = (Rᵀe₃)×Ω` and `δ(Rᵀe₃) = (Rᵀe₃)×η` steps) is stated and correct, and the
  stated result is the pinned string term for term. The rotating sphere (`TestRotatingSphere`)
  is the heaviest: I re-derived `p`, `ṗ`, `2ω×p`, `q×∂L/∂q` from
  `L̃ = ½mL²{|ω|² + 2qᵀS(Ωe₃)S(q)ω + Ω²qᵀCq}` (the book's form, p. 240), confirmed the
  `(ω×q)(q·a)` cancellation and that `−q×(5.43) + f` equals the compact form on the
  tangent plane (`q×q̈ = ω̇`, `q×PS(a)q̇ ≡ −ω×a`); the docstring is complete and the
  numeric `_book_543` matches eq. (5.43) as printed in the book. One subtlety a reader
  should know: `q×∂L/∂q` drops a `(q×a)(ω·q)` term that is zero only on the constraint set;
  the docstring's result is right but does not say that this step uses `ω ⊥ q`. Minor.
- `test_collect_fix.py`: gyroscopic coupling `(JΩ)·(r×b)`, inertial angular momentum,
  anisotropic q̇-form, SO(3)×SO(3) relative rate, offset-mass inertia — all reproduce the
  pinned strings. The SO(3)×SO(3) docstring lists the three `k`-cross terms in the engine's
  orientation (`a×b` vs `−b×a`) rather than the orientation its own `ṗ₁`/`p₁×Ω₁` lines
  give; the net `+kΩ₁×QΩ₂` it states is right, so it is reproducible, but a reader has to
  notice the swap.

None is incomplete; all can be reproduced in six months from the docstring plus the
module header's formula.

### 4. Error messages, docstrings, exports, paper sentence

- Every new or changed message names the term and the target/query (run, not read):
  `collect` both-sides (`'v' appears on both sides of 'Dot(Cross(a,v),Cross(b,v))'`),
  `_peel` (`the target 'v' occurs inside the matrix factor 'Hat(v)' of the term
  'Dot(Hat(v)a,b)'; transposition exposes …`), both-factors `NonLinearError`,
  `_no_left_coefficient` (`'Hat(x)' is linear in 'x' but has no left matrix coefficient`),
  `_zero_or_raise` (node type, query, expression), the dispatch fallthrough (`a MATRIX query
  has no left coefficient in a SCALAR expression | query='A', expr='ms'`), `Add()`
  `TypeError`, `assert_eom` (variable, term, tangent size, resolution). One wording slip:
  `extract.py:613` says "the map is a 3-tensor" also for the matrix-query case at `:722`,
  where the comment correctly says 4-tensor (A1).
- Docstrings: `evaluate`/`Evaluator` give the key convention (`str(leaf)`,
  `str(TimeDerivative(x))`, tangent vectors by `str(q.get_tangent_vector())`, zero/identity
  by type, `Variation` must be supplied); `check_eom` gives the sign convention
  `E = −(Euler–Lagrange) + forcing`, the three residual forms with book chapters, the
  `params` contract and the resolution limit; `assert_eom`, `TermResolution`, `DEFAULT_H`
  are documented; `to_standard_form` gained the S² projection note. `finite_difference_check`
  (moved, now public) documents the measure but leaves the sign convention of the rows to
  `linearized_rhs` (A3).
- `geomech/__init__.py` exports the five oracle names; `TermResolution` and `DEFAULT_H`
  are reachable via `geomech.numeric` only (fine). `AlgebraicError`/`NonLinearError` are not
  exported from `geomech` (pre-existing; now more visible, see follow-ups).
- Paper: both lines read `−2`. The sentence "the second form uses `qqᵀξ̇ = ω × ξ`, which
  follows from differentiating `ξᵀq = 0` along (eq:s2_kin). Note that δω has a component
  along q, since ω must remain orthogonal to the varied point q + δq; dropping it … gives
  wrong equations whenever ∂L/∂ω is not tangent to the sphere" is mathematically right
  (`δ(ω·q) = 0` forces `q·δω = −ω·δq = −q·(ω×ξ) ≠ 0`). It omits the second ingredient of
  the projection step, `ω × ξ ∥ q` (A2).

### 5. The F2 arbitration

Honest. The reviewer's two tests demanded, in one default `check_eom` call on the
`k = 1e7` point mass, `correct < 1e-5` **and** `flipped/dropped gravity > 1e-2`. I measured
the flipped score at `h = 1e-4`, `1e-3`, `1e-2`: **5.73e-7 at every step** — the blindness
is the normalisation by `|E_true| ≈ 1e7`, independent of roundoff, so no step change
rescues the assertion; and at the prescribed step the absolute roundoff (correct equation,
seed 0: 3.5e-9 relative ≈ 3e-2 absolute) is of the order of the gravity term itself, so a
per-term normalisation cannot separate correct from flipped either, as the fix log says.
The rewritten tests assert the behaviour the docs now declare: `assert_eom` raises
"below the resolution" naming `e3` (verified: it names the mass, gravity and forcing terms
with tangent sizes 1.65/1.83/1.80 against resolution 64), `require_resolved=False` accepts
on the relative error, and order-one parameters detect the dropped term. The review file
records that the tests were rewritten and why. Docs wording is accurate: "Resolution"
section, "the grey zone always exists", the old "no grey zone"/"off by order one" sentences
are gone. The only inaccuracy is the "m g is of order 1e-2" phrase in the test docstring
(required change 4).

### 6. Half-done / stale

- Required: fuzz builder (1), notebook (2), paper commit (3), adversarial docstring (4).
- `README.md:16-22` and `:88-89`, `docs/getting-started.md:44-50`, `geomech/dynamics/eom.py:51`
  give three pipeline orders, none with `collect`, two without manifold rules: pre-existing,
  plan item 4.7, not made worse here.
- `paper_context.md` is gitignored (`.gitignore:7`), so the corrected dispatch table the
  engine report cites is not in the repository; the contributor notes tell readers to
  consult it.
  Not a phase-1 defect, but the report's claim should be read with that in mind.
- CHANGELOG: complete and accurate for every item I could test (S² rule and the
  quadrotor-with-load normal term, `collect` including the matrix-factor rule and the
  `Dot(M a, N b)` restriction, `Hat(Vee)` and the `SkewSymmMatrix` flag, extraction holes,
  zero outer product, `attr` copying, empty `Add()`, oracle with the 36-system and
  `h = 1e-4` numbers, which I reproduced). It does not mention the notebook (because the
  notebook was not updated; after change 2 it needs no mention).
- Test names: seven `test_pinnedequation` (sic, no underscore) across two files,
  disambiguated by class only; `HALF` still defined in five test modules
  (`tests/test_eom_oracle.py:29`, `test_vbl_adversarial.py:56`, `test_scalar_adversarial.py:53`,
  `test_numeric.py:41`, `vbl_systems.py:26`) — the fix log's declared out-of-scope F8 remainder.
- `tests/test_phase1_adversarial.py` mentions xfail only in the sentence of change 4.

## Advisory

- A1. `geomech/core/math/extract.py:610-614` `_no_left_coefficient`: make the tensor rank
  follow the query type, e.g. `f"(the map is a {'4' if query.type is ExprType.MATRIX else '3'}-tensor, not C * query)"`,
  or drop the parenthesis and keep "has no left matrix coefficient".
- A2. Accompanying paper, background section: "which follows from differentiating
  `ξᵀq = 0` along (eq:s2_kin)" → "which follows from differentiating `ξᵀq = 0` along
  (eq:s2_kin) and from `ω × ξ ∥ q` for tangent `ω`, `ξ`".
- A3. `geomech/numeric.py:407-414` `finite_difference_check` docstring: state the row
  convention it checks (`(E(ε) − E(−ε))/2ε` against `−M δv̇ + A δx + B δu`, as
  `linearized_rhs` says) since the function is now part of the public API.
- A4. `geomech/numeric.py` `TermResolution` and `DEFAULT_H`: consider exporting from
  `geomech` next to `assert_eom`, since the assertion message refers to the resolution and
  users will want to tune `h`.
- A5. `tests/test_s2_variation_fix.py` `TestRotatingSphere` docstring: add "(using ω ⊥ q)"
  after `q × ∂L/∂q = mL²(−(q×ω)(q·a) + Ω² q×Cq)`; `tests/test_collect_fix.py:120-131`
  SO(3)×SO(3) docstring: write the three k-cross terms in the orientation of its own `ṗ₁`
  line (`+kΩ₁×QΩ₂ + k(QΩ₂)×Ω₁ + kΩ₁×QΩ₂`) and then say the engine prints the last two as
  `−kΩ₁×QΩ₂ − k(QΩ₂)×Ω₁`.
- A6. Rename the seven `test_pinnedequation` to `test_pinned_equation`.
- A7. `docs/manifold-implementation-notes.md` is still outside `mkdocs.yml` `nav`
  (strict build reports it); add it under a "Design notes" entry or keep it deliberately
  internal and say so in the contributor notes.
- A8. README has no pointer to the oracle; one line under "The pipeline" ("`check_eom`
  verifies a derived equation against a finite-difference Euler–Lagrange residual, see
  docs/numeric-verification.md") would make the new gate discoverable.

## Follow-ups surfaced by this branch for phases 2–5

| phase | item | where |
|---|---|---|
| 2.9 | the `r == vec` shortcut runs before the both-sides check, so `Dot(MVMul(VVMul(η,·),·), η)` (quadratic) is flipped and "ok" instead of `NonLinearError` (the one "ok" case of the matrix-factor fuzz campaign) | `geomech/core/math/collect.py:103` |
| 2 / 3 | `AlgebraicError`/`NonLinearError` are not importable from `geomech`; `compute_eom` now raises them on legal-looking input (`a·((Ω cᵀ) b)`), so users need `geomech.utils.errors` | `geomech/__init__.py:88-110`, `geomech/utils/__init__.py` |
| 3.1 / 3.2 | fuzz still reports row-vector `ExpressionMismatchError` on legal input, e.g. `(z_v)'(… + ((ξ)'0)' + …)` | `geomech/core/operations/multiplication.py:207` (seed 777 findings 1–2) |
| 3.4 | three `NotImplementedError`s remain in `extract.py` (unreachable: `Dot` operands are always vectors; `_scalar_from` on a typeless expr) | `geomech/core/math/extract.py:108, 376, 385` |
| 4.2 | the new `collect` `AlgebraicError` reaches the user from `compute_eom` without stage/variable context ("in term … while integrating by parts for ξ_q") | `geomech/core/math/collect.py:55-65`, `geomech/dynamics/eom.py:66-80` |
| 4.7 | three pipeline descriptions, none naming `collect`; `variables.py` says scalars are "not yet used in EOM pipeline" | `README.md:16-22, 88-89`, `docs/getting-started.md:44-50`, `geomech/dynamics/eom.py:51`, `geomech/dynamics/variables.py:12` |
| 4.7 | `paper_context.md` is gitignored yet cited by the contributor notes as a design document; decide whether it is tracked | `.gitignore:7` |
| 4.7 | `manifold-implementation-notes.md` not in nav | `mkdocs.yml:55-61` |
| 5.1 | vector rules do not descend into `MVMul`, so `Q(Ω×Ω)` and `Q(Ω̇)`-type terms are kept unreduced | `tests/test_collect_fix.py:131-137` (pinned strings show it), `geomech/core/transformations/vector_rules.py` |
| 5.4 / D2 | the three `kΩ₁×QΩ₂` terms of the SO(3)×SO(3) equation are not merged (commutative like-term limit); `_peel` relies on `has()` which is string containment | `tests/test_collect_fix.py:124-131`, `geomech/core/base/expressions.py:99-100` |
| 5.3 | `IdentityMatrix` is not flagged symmetric; `Transpose(skew)` could fold to `−skew` now that `is_skew_symmetric` propagates | `geomech/core/base/expressions.py` (identity leaf), `geomech/core/transformations/simplify.py:168-172` |
| 5.5 | `HALF` in five test modules; `test_eom_oracle.py` re-derives 12 pipeline systems; `_derive`-style helpers now live in `tests/phase1_systems.py` but `vbl_systems`/`test_eom_oracle` have their own | `tests/test_eom_oracle.py:29-170`, `tests/vbl_systems.py:26` |
| 5 / D1 | `Evaluator` ignores `value=` on `Vector`/`Matrix` leaves (only `Scalar.value` is honoured) | `geomech/numeric.py:154-167` |
| feature | `check_standard_form` / numeric check of `to_standard_form` blocks (only `linearize` has `finite_difference_check`) | `geomech/numeric.py` (oracle report "Not done") |
| 3.5 | `str(E)` is what the pinned tests compare; the printed strings in this branch (`(R1)'R2`, `(-1.0)` vs `(-1)`) will change when `to_latex`/printing is reworked, so phase 3 must regenerate 16 pinned strings | `tests/test_s2_variation_fix.py`, `tests/test_collect_fix.py` |
