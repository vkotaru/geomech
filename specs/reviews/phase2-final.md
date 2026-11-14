# Phase 2 final review: error instead of silence

Branch `fix/phase2-validation`, range `865699f..feeb476` (implementation c18fb55, adversarial
review ca4baf6, fixes feeb476). Read: the contributor notes, the phase 2 table of `specs/fix-plan.md`,
`specs/reviews/phase2-implementation-report.md`, `specs/reviews/phase2-review.md` with its fix
log, and `git diff 865699f..feeb476 -- geomech/ docs/ examples/ CHANGELOG.md` in full. No
file under `geomech/` or `tests/` was modified, and nothing was committed.

## Verdict: GO WITH CHANGES

All twelve rows are done (2.4 by verification). Every misuse I tried raises a typed error that
names the symbol. Every pinned system, example and notebook gives byte-identical output to
`865699f`. The only two output changes are intended: `rigid_pendulum.py` now prints `J2`
(2.12), and `sanity_check.py` now exits 0. The δW time-derivative check runs on `inf_work`
only. The four rewritten reviewer tests pin the fixed behaviour and lose no coverage.

The one required change is documentation. Brief item 4 asks for the behaviour of
`to_standard_form` and `linearize` on a system with `params` to be documented, and it is not
documented anywhere. `linearize` cannot handle such a system at all, and the hint in its error
message points to the wrong remedy.

## Required changes

1. **Document `to_standard_form`/`linearize` with `params`.** I verified the behaviour on a
   spherical pendulum with `l(t)` declared in `params=[l]`:
   - `to_standard_form` works. `M = {ω̇: −m l l I}`, `f` carries `l̇`
     (`ω m l l̇ (−2) − m g l q×e3`), and `G = {f: I}`.
   - `linearize(eom, V, inputs=[f])` raises:
     `AlgebraicError: linearize: in row '\xi_{q}', the term '\frac{d}{dt}(\omega_{q})(-1)m\delta{l}l' of the variation contains no perturbation symbol (configuration, velocity or input variation); the linearization is incomplete. Is an input missing from `inputs`, or a symbol missing the 'Constant' attribute?`
     Both remedies in that hint are wrong for a declared param. Following the first one
     (`inputs=[f, l]`) gives a second `AlgebraicError`: `δl` appears only inside
     `d/dt(δl)`.

   Add the following text:
   - `geomech/dynamics/eom.py:352` (end of the `compute_eom` docstring), after
     "``numeric.check_eom`` takes the same ``params``.": add "``to_standard_form`` accepts the
     result (a param's rate stays in ``f``, the param itself in ``M``); ``linearize`` does
     not support prescribed parameters yet and raises ``AlgebraicError`` on the ``δl``
     term."
   - `docs/vbl.md`, section "Limitations" (after the bullet that ends "or as the
     `AlgebraicError` described above."): add "- Prescribed time-varying parameters
     (`compute_eom(..., params=[l])`) are not supported: the variation of the equation
     contains `δl`, which is neither a state nor an input perturbation, and `linearize` raises
     `AlgebraicError` naming the term. Listing `l` in `inputs` does not help, because `l̇`
     brings in `d/dt(δl)`."
   - `CHANGELOG.md:16`, Added/`params` bullet: change "an undeclared one is an error (see
     Changed)." to "an undeclared one is an error (see Changed). `to_standard_form` accepts
     such equations; `linearize` does not yet (it raises `AlgebraicError` on `δl`)."

   I am not asking for a code change here. Whether `linearize` should take `params=` and drop
   `δl`/`d/dt(δl)` is a follow-up (F-1 below).

## Verified

### 1. Suite, lint, docs, fuzz, examples, notebooks

| Check | Command | Result |
|---|---|---|
| Suite | `python -m pytest tests/ -q -W error::RuntimeWarning` | **1699 passed, 5 xfailed**, 0 warnings promoted (6.2 s) |
| Format | `ruff format --check geomech/ tests/` | 76 files already formatted |
| Lint | `ruff check geomech/ tests/` | All checks passed |
| Docs | `mkdocs build --strict -d <tmp>` | exit 0 (INFO only: `manifold-implementation-notes.md` not in nav; `reference/` link) |
| Fuzz | `engine_fuzz.py --cases 2000 --seed 2024` at HEAD and at `865699f` | 149 findings in both: 9 documented-error-on-legal-input, 62 flag-dependent, 78 unsupported-linear-term. **0 wrong-value** in both. No new exception class. |
| Examples | all 6 `examples/*.py` at HEAD and `865699f`, stdout diffed | all exit 0 at HEAD. `action_integral`, `nary_check`, `point_mass`, `quick_check`: identical. `rigid_pendulum`: only `J1`→`J2` on the second link (2.12). `sanity_check`: exited 1 at base (`xi.delta()`), 0 now. |
| Notebooks | `jupyter nbconvert --execute --to notebook --stdout` for all 4 notebooks, in both trees | all exit 0, 0 error outputs; cell outputs (text/latex/plain) **identical** to `865699f` for `quadrotor`, `quadrotor_payload`, `rigid_pendulum`, `spherical_pendulum`. A probe notebook confirmed that the base kernel imported the base tree. |
| README/docs snippets | every ```` ```python ```` block of `README.md` and `docs/*.md`, executed in sequence per file, in both trees | identical results. The 5 failures in each tree are `NameError`s from fragments (`a`, `expr`, `L` undefined), the same at base. None comes from the new validation. |

Fuzz row shifts versus base, all accounted for:
- `collect[target occurs several times]`: ok 136 → 114, `NonLinearError` 164 → 188. This is
  2.9.
- `collect[target inside matrix factor]`: ok 2 → 0. This is 2.9.
- `collect[target occurs once]` 1233 → 1231 cases and `t_diff-rule` unsupported 2 → 3: same
  random stream, consumed differently because `wrap_calc` now uses `.t_diff()`/`.delta()` on
  manifold leaves (`engine_fuzz.py:629-634`).

### 2. Each row, one misuse through the public API

Each message below is pasted verbatim. The last column says whether a user six months from
now could fix their input from the message alone.

| # | Misuse | Message | Fixable? |
|---|---|---|---|
| 2.1 | `Scalar("m", attr=["constant"])` | `ValueError: attr: unknown attribute 'constant' for symbol 'm'; valid attributes are ['Constant', 'Identity', 'Manifold', 'Ones', 'OrthogonalMatrix', 'SkewSymmetricMatrix', 'SymmetricMatrix', 'TangentVector', 'UnitNorm', 'Zero']` | yes |
| 2.1 | `Scalar("m", attr="Constant")` | `TypeError: attr must be a list of strings, got str 'Constant' for symbol 'm'; write attr=['Constant']` | yes |
| 2.2 | unlisted `y` in L | `ValueError: compute_eom: δS still contains the variation symbols '\delta{y}' of the symbols 'y', which are neither entries of SystemVariables nor declared in params; for each symbol, list it in SystemVariables, declare it in params=..., or mark it 'Constant'` | yes |
| 2.2 | `e_3` without `Constant` | same text with `'\delta{e_3}'`/`'e_3'` | yes |
| 2.2 | `S2` point listed as `Vector("q")` | `… the variation symbols '\delta{\omega_{q}}', '\xi_{q}' of the symbols 'q' …` | yes; it names `q`, though not the S2/Vector mix-up |
| 2.2 | `l(t)` not declared | `… '\delta{l}' of the symbols 'l' …` (with `params=[l]` it derives `ξ_q`) | yes |
| 2.3 | `vectors=[x, Vector("x")]` | `ValueError: compute_eom: SystemVariables: the name 'x' is listed more than once (twice in SystemVariables.vectors); configuration variables must be distinct symbols` | yes |
| 2.3 | `scalars=[Scalar("x")], vectors=[x]` | `… (in SystemVariables.scalars and vectors) …` | yes |
| 2.3 | `vectors=[ω_q]` | `TypeError: compute_eom: SystemVariables.vectors entry '\omega_{q}' is a tangent or variation vector, not a configuration variable; list the manifold point 'q' instead` | yes |
| 2.3 | `matrices=[Matrix("J")]` | `TypeError: compute_eom: SystemVariables.matrices entry 'J' is a general Matrix; only SO3 rotation matrices are supported as matrix configuration variables` | yes |
| 2.3 | `SystemVariables(vectors=x)` | `TypeError: SystemVariables.vectors must be a list of configuration variables, got Vector 'x'; write vectors=[x]` | yes |
| 2.3 | `linearize(eom, SystemVariables(vectors=[x, x]))` | `ValueError: linearize: SystemVariables: the name 'x' is listed more than once (twice in SystemVariables.vectors); …` | yes |
| 2.3 | `linearize` keys mismatch | `ValueError: linearize: eom keys ['\\delta{x}'] do not match the variation vectors ['\\delta{z}'] of the given SystemVariables` | yes |
| 2.4 | `a=["Constant"]` passed to `Matrix`, `S2`, `SO3` | `a` is still `['Constant']` afterwards | n/a |
| 2.5 | `Dot(ω_q, S2("q"))` with two `S2("q")` objects | simplifies to `0` | n/a |
| 2.6 | `dW = f` (vector) | `TypeError: compute_eom: inf_work must be a scalar expression, got a VECTOR expression 'f'` | yes |
| 2.6 | `dW = "0"` / `False` | `… got str '0'` / `TypeError: compute_eom: inf_work must be a scalar expression, None or 0, got bool False`; `None` and `0` derive | yes |
| 2.7 | `_Opaque` node with a substituted child (the test's class, `tests/test_validation.py:501`) | `TypeError: apply_manifold_rules: cannot rebuild a _Opaque node (…) with substituted children; add the class to manifold_rules._rebuild` | yes (developer-facing) |
| 2.8 | `TimeDerivative(S2("q"))` | `ValueError: TimeDerivative: TimeDerivative(q) is not represented in the tree for the manifold point 'q'; use q.t_diff(), which returns the kinematic expression ω × q` | yes |
| 2.8 | `Variation(SO3("R"))` | `ValueError: Variation: Variation(R) is not represented in the tree for the manifold point 'R'; use R.delta(), which returns the kinematic expression R·Hat(η)` | yes |
| 2.8 | `Transpose()` | `TypeError: Transpose: Transpose() needs an operand; write Transpose(expr)` | yes |
| 2.9 | `collect(x, eta)` | `TypeError: collect: expected a scalar expression, got Vector 'x' of type VECTOR` | yes |
| 2.9 | `collect(Dot(Cross(eta, x), eta), eta)` | `NonLinearError: collect: 'eta' appears on both sides of 'Dot(Cross(eta,x),eta)'; the expression is not linear in the target` | yes |
| 2.10 | `full_simplify(Dot(Hat(q)·(ω×q), q), max_iter=1)` | `RuntimeWarning: full_simplify: no fixed point after 1 passes; result may be unsimplified: <cross(q, cross(\omega_{q}, q)), q>` | yes |
| 2.11 | `Number(0).is_zero`, `full_simplify(Number(0)*x)`, `x*0` | `True`, `0v`, `0v`; `Number(True)` → `TypeError: Number expects int, float, or str, got bool True; write Number(1) or Number(0)` | n/a / yes |
| 2.12 | `examples/rigid_pendulum.py` | second link prints `J2` | n/a |
| F1 | `Dot(x.t_diff().delta(), f)` | `ValueError: compute_eom: inf_work contains '\frac{d}{dt}(\delta{x})', the time derivative of the variation vector '\delta{x}'; virtual work must be linear in the variations themselves, not their time derivatives (write Dot(x.delta(), f), not Dot(x.t_diff().delta(), f))` | yes |
| F1 | `Dot(q.t_diff().delta(), f)` on S² | `… inf_work contains '\frac{d}{dt}(\xi_{q})' after the manifold kinematic substitutions, the time derivative of the variation vector '\xi_{q}'; …` | yes, but see A-2 |
| F2 | `params=[Vector("\omega_{q}")]` | `ValueError: compute_eom: params entry '\omega_{q}' has the printed name of the tangent vector of 'q'; rename the parameter (…)` | yes |
| F3 | `params=l` / empty `SystemVariables()` | `TypeError: compute_eom: params must be a sequence of leaves, got Scalar 'l'; write params=[l]` / `ValueError: compute_eom: no configuration variables listed in SystemVariables; …` | yes |
| F4 | hand-built `Variation(x.t_diff())` in δW | names `x` and points to `x.delta()` | yes (long, but correct) |

### 3. The four rewritten reviewer tests (`tests/test_phase2_adversarial.py`, `ca4baf6..feeb476`)

- `test_message_names_the_derivative_and_the_fix` replaces
  `test_observed_symptom_is_a_derivative_of_the_force`. It now asserts that the bad δW raises,
  that the message quotes `d/dt(δx)` and the "not Dot(x.t_diff().delta(), f)" hint, and,
  still, that the correct δW ends in `+f)`. The only assertion removed was the `−ḟ` symptom,
  which can no longer be produced. No coverage lost.
- `test_renamed_param_passes_the_oracle` replaces
  `test_observed_symptom_param_case_fails_the_oracle`. It keeps the oracle check on the renamed
  param (`< TOL`). The clashing case it used to evaluate now raises and is covered by
  `test_param_named_like_the_tangent_vector_raises` in the same class. No coverage lost.
- `test_string_attr_does_raise` still loops over every entry of `VALID_ATTRS` and now expects
  `TypeError` "attr must be a list of strings, got str". It asserts the fixed behaviour, and
  the property "never silently accepted" is kept.
- `test_bool_literals_are_rejected` replaces `test_bool_literals_are_ints`. It covers
  `Number(True)`, `Number(False)`, `Scalar * True` and `compute_eom(L, False, …)`. Every path
  the old test touched is still exercised. No coverage lost.

The 13 strict xfails had their markers removed with their assertions unchanged. I checked
this with `git diff ca4baf6..feeb476`.

### 4. `params` documentation

The `compute_eom` docstring (`eom.py:338-353`), `docs/getting-started.md:42-47`,
`docs/numeric-verification.md:153-157` and CHANGELOG (Added `params`, Changed stray check)
agree on four points:
- a symbol is a variable, a constant, or a declared param;
- an entry must be a distinct non-constant plain leaf;
- the rate of a param is kept and its variation is dropped;
- an undeclared param is an error.

`check_eom` takes the same `params`. What is missing is the `to_standard_form`/`linearize`
behaviour (required change 1).

### 5. The δW time-derivative check never sees δL

By reading: `_check_virtual_work(inf_work, variables, allowed)` (`eom.py:297`) is called once
(`eom.py:385`) with the user's `inf_work`, before `dL = lagrangian.delta()` is formed. It
works on its own copy (`full_simplify(apply_manifold_rules(inf_work, …))`), and δS is built
from the untouched `inf_work`.

By running: I wrapped `_check_virtual_work` with a spy and derived the pinned builders.

| System | Argument the check received | L passed? | `check_eom` error |
|---|---|---|---|
| spherical pendulum (`vbl_systems.spherical_pendulum_s2`) | `Dot(ξ_q, f)` | no | 2.3e-8 |
| free rigid body | `Dot(η_R, M)` | no | 4.8e-9 |
| quadrotor with load | `(Dot((δx_L + Cross(ξ_q,q) l (−1)), R e3 f) + Dot(η_R, M))` | no | 4.0e-8 |
| spring–mass–damper (`test_eom_oracle.spring_mass_damper`) | `(δx f + c ẋ δx (−1))` | no | 7.2e-10 |
| elastic spherical pendulum (with scalar `l`) | `(δl f_l + Dot(ξ_q, f_q))` | no | 2.2e-8 |

Each equation has the expected form. For example, the rigid body gives
`J Ω̇ (−1) − Ω×JΩ + M`, which is Lee et al. (6.5), and the spherical pendulum gives
`−m l² ω̇ − m g l q×e3 + f`. These ξ̇/η̇-bearing Lagrangians all derive, so δL is not checked.

I also ran four false-positive probes on δW, and all derive: `Dot(δx, ḟ)`, `Dot(δx, d/dt(l f))`
with `params=[l]`, and `Dot(δx, m ẍ)`. `Dot(ẏ.delta(), f)` with an unlisted `y` is caught by
the stray check, which names `y`.

### 6. Byte-identical pinned output against `865699f`

I made a temporary git worktree at `865699f` and have since removed it with
`git worktree remove` and `git worktree prune`. In both trees I ran every argument-free system
builder in `tests/vbl_systems.py`, `tests/phase1_systems.py` and `tests/test_eom_oracle.py`
(61 distinct systems; the union of the 85 list entries the review counted) through
`compute_eom → to_standard_form (repr) → linearize (str + every attribute)`. The output
matched exactly (`cmp`, 206 253 bytes each). The same three pre-existing typed errors appear
in both trees:
- one `linearize` `AlgebraicError` (`adv_virtual_work_under_rotation`);
- two `collect` `AlgebraicError`s (outer product on velocity, S² and SO(3)).

Examples and notebooks: see section 1.

## Advisory (not blocking)

- **A-1 `_check_variables` accepts composite nodes.** Entries such as `q.t_diff()` (a
  `Cross`), `x + y` and `m*z` get past `eom.py:35-131`, because `get_variation_vector()` falls
  back to `delta()`. They are not silent: the stray check raises later. The message, though,
  says "δx of the symbols 'x' … list it", not "entry 'Cross(ω,q)' is not a leaf". A check
  that the entry is a `Scalar`/`Vector`/`S2`/`SO3` leaf, like the one `_check_params` already
  makes (`eom.py:187`), would say what is wrong.
- **A-2 F1 hint is R³-shaped.** `eom.py:315-317` always suggests
  `Dot(x.delta(), f), not Dot(x.t_diff().delta(), f)`. For S²/SO(3), where the user wrote
  `ξ̇·f` or `δ(q̇)·f`, a hint built from the actual variation vector (`Dot(ξ_q, f)`) would be
  clearer.
- **A-3 `_check_eom` length comparison is unreachable** (`vbl.py:368`). `linearize` runs
  `_check_variables` first, and that already rejects duplicates. This is harmless belt and
  braces.
- **A-4 Message prefixes are still mixed.** The F5 message (`attr must be a list…`,
  `expressions.py:26`) and the `SystemVariables(...)` constructor message (`variables.py:16`)
  carry no `attr:`/`SystemVariables:` prefix, while the 2.1 message does (`attr: unknown
  attribute`).
- **A-5 `docs/numeric-verification.md:166`.** Step 1 ("call `compute_eom`") could say "with
  the same `params`". Step 2 already says it for `check_eom`.
- **A-6 Existing tests were edited.** The four `params=[…]` additions, the five-for-one
  `Transpose()` replacement and the `"Symmetric"` → `"SymmetricMatrix"` typo were coordinator
  decisions recorded in the implementation report. I checked each one: the assertions are
  unchanged, apart from the deleted `Transpose()`-is-None tests, which pinned the defect.
- **A-7 Side effect from the implementer.** The implementer reported that the docs extras may
  have been installed into an environment outside the development virtual environment.
  This is worth checking, but it does not affect the branch.

## Follow-ups this branch surfaced (phases 3–5)

- **F-1 (phase 4/5) `linearize(…, params=[…])`.** `vbl.py:181` has no `params` argument, and
  the hint at `vbl.py:469` ("Is an input missing … or a symbol missing the 'Constant'
  attribute?") is wrong for a declared param. Either accept `params` and drop
  `δl`/`d/dt(δl)` (zero along a prescribed reference), or name params in the hint.
- **F-2 (phase 4) Composite configuration-variable entries.** See A-1 (`eom.py:35`). This
  belongs with the phase 4.4 invariant checker that the review mentions for reserved printed
  names anywhere in `L`.
- **F-3 (phase 3.3) Numeric operands other than `int`/`float`.** `numpy.int64` is neither, so
  `x * np.int64(2)` still fails with `AttributeError` (`expressions.py:157,249,414` check
  `(int, float)` only). `bool` is now rejected (`expressions.py:540`).
- **F-4 (phase 3.4) `ExpressionMismatchError` from the row-vector encoding.** Fuzz findings 1
  and 2 (`simplify`/`full_simplify`: `MVMul cannot be performed between VECTOR and MATRIX`,
  `multiplication.py:148,207`) are unchanged from base. They are phase 3.1/3.2.
- **F-5 (cosmetic) Test class name.** `tests/test_validation.py:126` `TestStrayManifoldVariation`
  also covers the R³/scalar strays, so the name undersells it. The fixer deliberately did not
  rename it.
- **F-6 (phase 4) `full_simplify` cap warning.** The warning fires when the last allowed pass
  reached the fixed point but no confirming pass ran (`simplify.py:56-74`). This is
  documented. If phase 4 adds a "converged" flag to the debug output, it can say "unconfirmed"
  instead.
