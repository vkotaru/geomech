# Phase 5b final review, and closure of the fix plan

Branch `fix/phase5b-structure`, HEAD d102af4. Phase 5b range `fc0d386..d102af4`; whole plan
`54d97b9..d102af4` (phases 0–5b, stacked on `feature/vbl`). Working tree clean at review
time. Everything below was rerun at HEAD with the repo venv (Python 3.12.13).

## Part A — phase 5b verdict: GO

No required change. Every gate passes; the five review findings (F1–F5) and the two
advisories (mypy pin, exact-class dispatch) are fixed and pinned by tests.

| gate | command | result |
|---|---|---|
| suite, warnings as errors | `python -m pytest tests/ -q -W error` | 2307 passed, 2 xfailed (both pre-plan strict xfails in `tests/test_scalar_adversarial.py:785,793`) |
| suite, debug mode | `GEOMECH_DEBUG=1 python -m pytest tests/ -q` | 2307 passed, 2 xfailed |
| ruff | `ruff format --check` / `ruff check` on `geomech/ tests/` and `specs/reviews/scripts/` | clean (99 + 4 files) |
| mypy | `mypy geomech --ignore-missing-imports --python-version {3.10,3.11,3.12}` (mypy 2.4.1) | `Success: no issues found in 34 source files` ×3 (one `annotation-unchecked` note, `mixins.py:200`) |
| docs | `mkdocs build --strict` | builds; one INFO (`index.md` link `reference/`) |
| fuzz | `engine_fuzz.py --cases 2000 --seed 1234` | path line `fuzzing geomech from <repo>/geomech` (this checkout); **0 wrong-value**; every value-preserving transform 2000/2000 ok, convergence 2000/2000; 188 distinct findings: 114 flag-dependent, 69 unsupported-linear-term, 5 documented-error-on-legal-input (all `collect` `NonLinearError`/`AlgebraicError` on trees that are nonlinear or have the target in a matrix factor) |
| examples | `python examples/*.py` | all six exit 0 |
| notebooks | re-executed in a copy with `nbconvert --execute` | four notebooks, every code-cell output identical to the committed one |
| CI file | `.github/workflows/tests.yml` | lint job: `ruff format --check`, `ruff check`, then `mypy geomech --ignore-missing-imports` (ruff, `mypy==2.4.1`, numpy installed); test job `needs: lint`, matrix 3.10/3.11/3.12, `pip install -e ".[dev]"`, `pytest tests/ -v`. Correct |
| docs/README vs `*` and name rule | grep of README, `docs/`, examples, notebooks, a downstream application | no juxtaposed product in printed output; no leaf name containing `*` outside `docs/conventions.md:119` (the rule itself) and the guard docstring; `conventions.md:105-133` states the separator, the name rule, the `Transpose` operand class and frozen nodes |

Spot checks of the fix log at HEAD: `Transpose(Vector("x")) != Transpose(Matrix("x"))`
(also under `Variation`); `_latex_str_key("u^*") == "u^*"`, `"R*e3*f" → "R e3 f"`;
pickle and `deepcopy` keep `==`, `hash`, `str`, tuple operands; `S2`, `SO3`, `Vector`,
`Matrix`, `Number`, `SkewSymmMatrix` all reject `"a*b"`; nodes and `ExprFlags` raise
`FrozenInstanceError` on assignment.

Advisory (not blocking):

- `specs/reviews/phase5b-implementation-report.md:3` still says "Not committed", and its
  gate table (`:12-14`) gives the pre-fix count 2267; the fix log
  (`phase5b-review.md:211`) has the current 2307.
- The paper fix of plan item 1.1 (background and toolbox sections of the accompanying
  paper, now `ξ̇ − 2 ω × ξ`) exists only as an uncommitted change in the paper's
  working tree.

## Part B — closure of the whole plan

Rule used: **DONE** = fixed at HEAD, with the commit and a check I ran today; **DEFERRED** =
not done, with an explicit deferral or a user-facing statement of the limitation in
`docs/`, `CHANGELOG.md` or a later phase report; **OPEN** = not done and recorded only in
the review that surfaced it. Commits: phase 0 `1dd98ab`; phase 1 `636e963`, `c96d31c`
(+ `8e28578`, `865699f`); phase 2 `c18fb55` (+ `feeb476`, `1e10349`); phase 3 `475c596`
(+ `743a29b`, `57232da`); phase 4 `bffba73` (+ `86e7f26`, `47cce2d`); phase 5a `cbd3653`
(+ `36122e1`, `fc0d386`); phase 5b `f852c11` (+ `d102af4`).

### Dynamics-audit rows and engine F1/F2 reproducers, rerun at HEAD

Rebuilt from the public API without the test builders; `check_eom` seeds 0/1 at the
default `h`, and `assert_eom` (with the resolution check) on each. No crash.

| row | seeds 0 / 1 |
|---|---|
| F1 anisotropic `½ ω·Aω` | 2.6e-8 / 8.7e-9 |
| F1 rotating frame | 8.3e-9 / 1.0e-8 |
| F1 gyroscopic `c·ω` (prints `Cross(\omega_{q},c)*(2)`) | 5.9e-9 / 7.0e-9 |
| F1 R³×S² `ẋ·ω` | 8.9e-9 / 9.3e-9 |
| F1 S²×SO(3) `k ω·Ω` | 5.3e-9 / 3.0e-8 |
| F1 S²×SO(3) `k ω·(Rᵀe₃)` | 3.8e-8 / 1.1e-8 |
| F3 `(R(JΩ))·e₃`, `((RJ)Ω)·e₃` | 6.9e-9 / 4.5e-9; 9.0e-9 / 4.5e-9 |
| F3 S² `½ q̇·A q̇` | 4.0e-8 / 7.8e-8 |
| F3 S² `½ m q̇·q̇ + k (A q̇)·(q × e₃)` | 2.4e-8 / 4.0e-8 |
| F3 SO(3)×SO(3) `k Ω₁·(R₁ᵀ(R₂Ω₂))` | 2.2e-8 / 2.2e-8 |
| engine F6 `½(Ω×ρ)·J(Ω×ρ)` | 3.0e-9 / 1.2e-8 |
| engine F1 pipeline `½ Ω·JΩ + (JΩ)·(q×e₃) + ½ ω·ω` | 1.7e-8 / 6.6e-9 |
| dynamics F2 `(JΩ)·((Rᵀe₃)×b)`, `k (RΩ)·(a×b)` | 1.4e-8 / 7.1e-9; 1.5e-8 / 4.2e-9 |

Engine F1 unit: `collect(Dot(M v, b×c), v)` → `Dot(v,(M)'*Cross(b,c))`, value equal to
4e-15 over 5 random draws (general and symmetric `M`). Engine F2: `full_simplify(Hat(Vee(I)))`
stays `Hat(Vee(I))`, symmetric `J` stays `Hat(Vee(J))`, `Rᵀ(Hat(Vee(I)) e₃)` evaluates to 0,
skew `S` → `S`. Worst score overall 7.8e-8 (gate ≤ 1e-5).

### Fix plan rows (phases 0–5, decisions)

| item | status | proof at HEAD |
|---|---|---|
| 0 audits + harness | DONE 1dd98ab | four reports and `scripts/euler_lagrange_check.py` tracked |
| 1.1 S² rule | DONE c96d31c | six F1 rows above; `TS2.delta(substitute=True)` gone (`TypeError`); paper fixed but uncommitted in the accompanying paper |
| 1.2 collect | DONE c96d31c | engine F1 unit/pipeline and F3 rows above |
| 1.3 Hat(Vee) | DONE c96d31c | engine F2 above |
| 1.4 extract holes | DONE c96d31c | `Hat(x)`, `Transpose(Hat(x))`, `M·Hat(x)` → `AlgebraicError`; `Hat(x)·x`, `(Mx)×x` → `NonLinearError`; matrix query in scalar → `AlgebraicError`; `paper_context.md` table corrected (file is gitignored) |
| 1.5 zero outer product | DONE c96d31c | `simplify(VVMul(0v, aᵀ)).type == MATRIX` |
| 1.6 oracle | DONE 636e963 (deviation) | `geomech.numeric.check_eom` used above; the gate is `tests/test_eom_oracle.py` (re-derived systems), not a `check_eom` line in `test_manifold_pipelines.py`/`vbl_systems.py` (0 occurrences); deviation recorded in `phase1-oracle-report.md:94-99` |
| 1 done-criteria | DONE | fuzz 0 wrong-value; rows ≤ 1e-5 |
| 2.1–2.12 | DONE c18fb55 | each reproduced: `attr=["constant"]` ValueError; stray `δy`, `δm` ValueError; duplicate `q`/`x` ValueError, tangent vector and plain `Matrix` TypeError; shared attr list untouched; second `S2("q")` gives `Dot(ω,q2) → 0`; `dW` 0/None accepted, vector `dW` TypeError; `_Opaque` rebuild TypeError (`tests/test_validation.py:502`); `TimeDerivative(q)` ValueError, `Transpose()` TypeError; collect both-sides NonLinearError, non-scalar TypeError; cap RuntimeWarning; `Number(0).is_zero`; `examples/rigid_pendulum.py:32` is `J2` |
| 3.1–3.6 | DONE 475c596 | `Rᵀ*v/M/s` build; `vᵀ*s → (v*s)'`; `Transpose(Transpose(x)) is x`; `2*x`, `m+1`, `1+m`, `m-1`, `-x`, `x/2`; `Dot cannot combine SCALAR 'm' and VECTOR 'x'`, `Hat expects a VECTOR, got SCALAR 'm'`; `a \times \left(b \times c\right)` vs `\left(a \times b\right) \times c`, `\left(b \cdot c\right) a`; `repr` `v:x`, `M:J`, `S2:q`, `k*`, indent `print_tree` with flags |
| 4.1–4.7 | DONE bffba73 | `EOMTrace().summary()` lists 13 stages; engine-C3 term names the term and stage; jerk message distinct; `GEOMECH_DEBUG=1` suite green; `set_trace` logs `simplify.hat_vee_skew: …`; renamed modules import; `pipeline.md`, `conventions.md`, `troubleshooting.md` build strict |
| 5.1–5.4 | DONE cbd3653 | `Dot(v,q×(ω×q)) → Dot(v,ω)`, `Hat(a×a) → 0`, `J(a×a) → 0v`; `Vee(A+B)` distributed; `Iᵀ → I`, `Sᵀ → S*(-1)`, `(ARᵀ)R → A`; `(m*2)*y+(m*3)*y → m*y*(5)` |
| 5.5–5.7 | DONE f852c11, d102af4 | `constants` module; `is_leaf`/`has_nested_add` gone; splits present; mypy clean ×3 |
| D1 | DONE cbd3653 | `m(=2)+l(=3) → (m+l)` |
| D2 | DONE f852c11, d102af4 | `Mul(ab,c) != Mul(a,bc)`; `Variation(Scalar x) != Variation(Vector x)`; `Dot(Vector m, w).has(Scalar m)` is False; `*` reserved |
| D3 | DONE f852c11 | frozen nodes and flags, tuple operands |
| D4 | DONE c18fb55 | = 2.5 |
| D5 | DONE bffba73 | `trace=` keyword |

### Original audit findings

Many findings overlap (one fix closes `eng-F1`, `dyn-F2`, `fuzz-1`); each is counted in its own audit.

| finding | status | proof / where recorded |
|---|---|---|
| eng-F1, F2, F3, F4, F6 | DONE c96d31c | reproducers above |
| eng-F5, F7 | DONE 475c596 | 3.1, 3.3 |
| eng-F8, F9, F10, F12 | DONE cbd3653 | 5.1, 5.3, 5.4, D1 |
| eng-F11, F14, F15 | DONE c18fb55 | 2.11, 2.1, 2.4 |
| eng-F13, F16 | DONE f852c11, d102af4 | D2, D3 |
| eng-F17 `vᵀ*s` double transpose; Transpose as row marker | DONE 475c596 | `(v*s)'`, `to_latex` `s v^T`; `is_row_vector` |
| eng-F17 S² rules by identity; S² point listed as `Vector("q")` | DONE c18fb55 | `Dot(ω,q2) → 0`; stray `δω_q` ValueError |
| eng-F17 float sum `(0.30000000000000004)` | DEFERRED | `phase5a-implementation-report.md:97` (non-integral floats kept) |
| **eng-F17 `TimeIntegral(c).is_constant` is True** | **OPEN** | still `True` at HEAD; recorded only in `engine-audit.md:444` |
| dyn-F1, F2, F3 | DONE c96d31c | rows above |
| dyn-F4, F5, F6 | DONE c18fb55 | 2.3, 2.2, 2.6 |
| dyn-F7, F8 | DONE 475c596 | 3.1, 3.5 |
| dyn-F9 linearize with time-varying params | DEFERRED | `docs/vbl.md:204-207` |
| dyn nit `variables.py` docstring | DONE | `variables.py:26` "scalar configuration variables" |
| dyn nit tangent vector / plain `Matrix` as configuration | DONE c18fb55 | TypeError for both |
| **dyn nit `to_standard_form` M-block shape for a scalar row with a vector acceleration** | **OPEN** | `StandardFormEquation` docstring (`standard_form.py:26-39`) does not say `M[a]` means `M[a]ᵀ a` |
| **dyn nit input with the name of a constant leaf** | **OPEN** | `linearize(eom, V, [f, Scalar("g")])` on `… − m g e₃·x` returns no `B[δg]` column and no error; `to_standard_form` gives `G[g] = −m e₃`; only the general name rule is documented (`conventions.md:111`) |
| **dyn nit `(d/dt(J))ᵀ` of a symmetric J not simplified** | **OPEN** | `full_simplify(Transpose(J.t_diff()))` → `(\frac{d}{dt}(J))'` |
| dbg A3–A9, A16, A29–A33 | DONE 475c596 | messages above; `m*"x"` TypeError; division TypeError with hint |
| dbg A18, A23, A24, A25, A27 | DONE c18fb55 | collect TypeError; stray `δm`, `δy` ValueError; `ẋ.delta()` in `dW` ValueError; `J2` |
| dbg A26, C3 | DONE bffba73 | jerk vs IBP messages differ |
| dbg C1 | DONE cbd3653 | 5.1 |
| dbg C2 | DONE c96d31c | 1.4 |
| dbg C6 | DONE c18fb55 | 2.5 |
| dbg F1, F2, F12, F17 | DONE c18fb55 | 2.1, 2.2, 2.10, 2.5/2.9 |
| dbg F3, F4, F9, F10, F11, F14 | DONE bffba73 | 4.1, 4.2, 4.5, 4.7, 4.3, 4.6 |
| dbg F5, F6 | DONE 475c596 | 3.4, 3.5/3.6 |
| dbg F7 | DONE c96d31c | 1.4 |
| dbg F8 | DONE c18fb55, bffba73 | 2.8 guards, 4.4 `check_tree` |
| dbg F13 | DONE 636e963 | `geomech.numeric` |
| dbg F15 duplicates, constants, public names, dead code, vocabulary, `_remove_factor` | DONE f852c11 | constants module; `UndefinedCaseError` gone; `x is One` |
| dbg F15 `vector_rules._at_node` (80 lines) not split | DEFERRED | `phase5b-implementation-report.md:249` |
| **dbg F15 magic numbers** | **OPEN** | `full_simplify(max_iter=10)`, `printing.py:214` `len(name) <= 10`, `:278` `gap=3`, `:966` `[-2000:]` unchanged; not in plan 5.5, not recorded later |
| dbg F16 | DONE f852c11 | 5.7 |
| fuzz 1, 2, 3, 4, 7, 8, 9 | DONE c96d31c | 1.1–1.5 |
| fuzz 5 | DONE c18fb55 | 2.9 |
| fuzz 6 | DONE 475c596 | 3.2 |
| fuzz 10, 11 | DONE cbd3653 | 5.2, 5.1 |
| fuzz 12 flag-dependent rules | DEFERRED (by design) | `engine-fuzz-report.md:373` "informational, not bugs"; 114 at seed 1234 |
| fuzz 13 linear terms refused with AlgebraicError | DEFERRED (by design) | `docs/troubleshooting.md` AlgebraicError causes; 69 at seed 1234 |

### Follow-ups and advisories of the phase final reviews

Informational advisories with no action (phase 2 A-3, A-6, A-7; phase 3 A3, A7) are not counted.

| source | item | status | proof / where |
|---|---|---|---|
| P1 fu | `r == vec` shortcut before both-sides check | DONE c18fb55 | fuzz-5 reproducer NonLinearError |
| P1 fu | errors importable from `geomech` | DONE c18fb55 | `geomech.AlgebraicError` etc. |
| P1 fu | row-vector ExpressionMismatchError | DONE 475c596 | no such error in today's fuzz |
| P1 fu | three `NotImplementedError` in `extract.py` (`:111, :377, :386`) | DEFERRED | unreachable, carry messages; 3.4 removed only the empty ones (`phase3-implementation-report.md:152-155`) |
| P1 fu | collect error context | DONE bffba73 | "a `collect` error during integration by parts is prefixed with the stage" |
| P1 fu | three pipeline descriptions; `variables.py` | DONE bffba73 | `docs/pipeline.md` |
| P1 fu | `paper_context.md` gitignored but cited | DEFERRED | author decision, `phase4-final.md:186` |
| P1 fu | manifold notes not in nav | DONE | `mkdocs.yml:64` |
| P1 fu | `Q(Ω×Ω)` unreduced; `kΩ₁×QΩ₂` unmerged; `_peel` by string `has()` | DONE cbd3653, f852c11 | 5a report table; `has()` uses `==` |
| P1 fu | Identity symmetric, `Transpose(skew)` | DONE cbd3653 | 5.3 |
| **P1 fu** | **`HALF` defined in 7 test modules; duplicated system builders in tests** | **OPEN** | `tests/test_numeric.py`, `test_scalar_adversarial.py`, `vbl_systems.py`, `test_eom_oracle.py`, `test_validation.py`, `test_vbl_adversarial.py`, `test_phase2_adversarial.py`; 5b deduplicated the library only |
| **P1 fu** | **`Evaluator` ignores `value=` on `Vector`/`Matrix` leaves** | **OPEN** | `evaluate(Vector("u", value=[1,2,3], attr=["Constant"]))` returns a random vector |
| **P1 fu** | **`check_standard_form` (numeric check of `to_standard_form`)** | **OPEN** (feature) | not in `geomech.numeric` |
| P1 fu | regenerate pinned strings after 3.5 | DONE 475c596 | |
| P1 A3 | `finite_difference_check` docstring convention | DONE | `numeric.py:419` |
| P1 A7 | notes in nav | DONE | `mkdocs.yml:64` |
| **P1 A1** | **`_no_left_coefficient` tensor rank** | **OPEN** | a matrix query (`SMMul(A, aᵀJw)` for `J`) still says "3-tensor" |
| **P1 A2** | **paper sentence "and from ω × ξ ∥ q"** | **OPEN** | accompanying paper, background section, unchanged |
| **P1 A4** | **export `TermResolution`, `DEFAULT_H` from `geomech`** | **OPEN** | `hasattr(geomech, …)` False |
| **P1 A5** | **test docstrings (`TestRotatingSphere` "using ω ⊥ q"; SO(3)×SO(3) orientation)** | **OPEN** | not changed |
| **P1 A6** | **rename `test_pinnedequation`** | **OPEN** | 5 remain |
| **P1 A8** | **README pointer to `check_eom`** | **OPEN** | README has no mention |
| P2 F-1 | `linearize(params=)` support | DEFERRED | `docs/vbl.md:204-207` |
| P2 F-2 / A-1 | composite `SystemVariables` entries | DONE | CHANGELOG.md:274 TypeError |
| P2 F-3 | numpy scalar operands | DONE 475c596 | `x*np.int64(2)` → `x*(2)` |
| P2 F-4 | row-vector errors | DONE 475c596 | |
| **P2 F-5** | **`TestStrayManifoldVariation` name** | **OPEN** | `tests/test_validation.py:126` |
| P2 F-6 | cap warning wording | DONE bffba73 | "no confirmed fixed point" |
| **P2 A-2** | **`inf_work` hint is R³-shaped** | **OPEN** | an S² `ξ̇·f` still gets "write Dot(x.delta(), f)" (`eom.py:331`) |
| **P2 A-4** | **message prefixes mixed** | **OPEN** | `attr must be a list…` has no `attr:` prefix |
| P2 A-5 | numeric-verification "same params" | DONE | `docs/numeric-verification.md:156` |
| P3 F1 | row in kinetic energy | DONE bffba73, 475c596 | collect rewrites `VVMul(rᵀ,b)`; CHANGELOG.md:103 |
| P3 F2 | rows in `Dot`/`Cross`/`Hat`, mixed `VAdd` | DONE bffba73 (in `check_tree`); construction guard DEFERRED | `docs/conventions.md:51-53` |
| P3 F3, F4 | general hoist; `a*x − a*x` | DONE cbd3653 | `full_simplify(a*x − a*x)` → `0v` |
| **P3 F5** | **`_is_row`/`_q` one-line shims in `multiplication.py:35,42`** | **OPEN** | `wrap_numeric` made public; shims remain |
| P3 F6 | split `to_latex` | DONE f852c11 | `_latex_handlers()` |
| P3 F7 | conventions: row forms, `0 + x`, tolerance | DONE | `conventions.md:39-43, 64-65, 94` |
| P3 A1, A6 | `Kx(-1.0)`; tolerance documented | DONE cbd3653 | ints normalised; `conventions.md:94` |
| **P3 A2** | **pin `str` of `row_spring`/`row_quadratic_form`** | **OPEN** | `tests/test_row_forms.py` checks them numerically only |
| **P3 A4** | **`display_latex` keeps outer parentheses** | **OPEN** | does not call `_strip_outer_parens` |
| **P3 A5** | **`__rtruediv__` quotes `{other!r}` unshortened** | **OPEN** | `expressions.py:260` |
| **P3 A8** | **`print_tree` uses `├─` for the last child** | **OPEN** | no `└─` in `printing.py` |
| P4 fu 1, 2 | vector rules descend; Vee over MAdd | DONE cbd3653 | |
| P4 fu 3 | construction-time row rejection | DEFERRED | `docs/conventions.md:51-53` |
| P4 fu 4 | `linearize` params | DEFERRED | `docs/vbl.md:204` |
| P4 fu 5, 6 | `zero_like` for constant variations; split `_extraction_error`, `_eliminate` | DONE f852c11 | `0v`/`0`; `_higher_order_error` … |
| **P4 fu 7** | **`check_rewrite(before, after)` rule-level numeric check** | **OPEN** | not in `geomech` |
| P4 fu 8, 9 | rule-log readability; `collect.flip`/`ibp.product_rule` names | DONE cbd3653 | |
| P4 A1, A2, A3 | log lines; nested negation `repr`; exception context | DONE cbd3653, f852c11 | `_ExtractionFailure` |
| **P4 A4** | **fuzzer lives under `specs/reviews/scripts/` though user docs point to it** | **OPEN** | `docs/pipeline.md:252`, `troubleshooting.md:38,289` |
| P4 A5 | `paper_context.md` gitignored | DEFERRED | author decision |
| P5a fu 1–4, 8 | inner-product merge; one-pass grouping; constants; dead code; local imports | DONE f852c11 | 5b report steps 1–2 |
| P5a fu 5 | symbolic coefficients not factored | DEFERRED | `docs/conventions.md:86` |
| P5a fu 6 | `q × (ω̇ × q)` kept | DEFERRED | `docs/troubleshooting.md:17` |
| P5a fu 7 | S² linearize blocks not unique | DEFERRED | CHANGELOG.md:232-236 |
| P5a A1, A2 | bare `Number(2.0)`; literal spelling | DONE f852c11 | `_canonical_literal` |
| P5a A3 | stale test name `test_f_is_mathematically_correct_but_uncombined` | DEFERRED | `phase5b-implementation-report.md:128-132` |

### Counts

| group | DONE | DEFERRED | OPEN |
|---|---|---|---|
| fix plan rows, done-criteria, decisions (45) | 45 | 0 | 0 |
| engine audit (22, F17 split) | 20 | 1 | 1 |
| dynamics audit (14, nits split) | 10 | 1 | 3 |
| debuggability audit (42, F15 split) | 40 | 1 | 1 |
| engine fuzz report (13) | 11 | 2 | 0 |
| phase final follow-ups and advisories (71) | 41 | 11 | 19 |
| **total (207)** | **167** | **16** | **24** |

All 24 OPEN items are cosmetic, documentation, test-hygiene or small features. None
changes a derived equation; the only behavioural one is the dyn nit on an input named like
a constant leaf (silently missing `B[δg]` column in `linearize`), which is a user-error
case. This file now records them, so a later pass can treat it as their register.

## Part C — contributor notes and pre-push notes

### Stale lines in the contributor notes

| line | replace with |
|---|---|
| 23 | The development virtual environment at `venv/` (Python 3.12) works: `source venv/bin/activate`; use `python -m pip`. If the venv breaks, recreate it: |
| 26 | `rm -rf venv && uv venv --seed --python 3.12 venv && source venv/bin/activate   # only if the venv is broken (the current venv was made by uv; python3.12 -m venv venv also works)` |
| 27 | `python -m pip install -e ".[dev]"       # editable install + pytest, ruff, mypy (pinned 2.4.1), jupyterlab` |
| 28 | `python -m pytest tests/ -q              # full suite (~2300 tests, numpy is the only runtime dep)` |
| 31 | `ruff format geomech/ tests/ && ruff check geomech/ tests/ && mypy geomech --ignore-missing-imports   # what CI's lint job runs (with --check for format)` |
| new after 31 | `python specs/reviews/scripts/engine_fuzz.py --quick                 # numeric fuzz of every rewrite; must report 0 wrong-value` |
| 32 | `python -m pip install -e ".[docs]" && mkdocs build --strict          # docs; API pages generated by scripts/gen_ref_pages.py (mkdocs serve to browse)` |
| 35 | CI (`.github/workflows/tests.yml`) lints first (`ruff format --check`, `ruff check`, `mypy geomech --ignore-missing-imports`, on 3.12), then tests on 3.10–3.12. Target branch for PRs is `geomech`, not `main`. Docs deploy from `main`/`geomech` via `mkdocs gh-deploy`. |
| 41 | Everything is an immutable **expression tree**: node classes are frozen dataclasses, operands are a tuple (`expr.nodes`), `ExprFlags` is frozen; build a new node instead of changing one (`object.__setattr__` is the only way past it, used by constructors and by tests that build malformed trees). Python operators on nodes build the tree; every algorithm is a recursive rewrite, dispatched by node class. Layers, in import order (see `geomech/__init__.py`): |
| 43 (append) | Unknown `attr` strings raise `ValueError`; `TimeDerivative(q)`/`Variation(q)` of a manifold leaf raise `ValueError`. `constants.py` holds `NEG_ONE`, `NEG_TWO`, `HALF`, `negate(expr)`, `zero_like(expr)`; rewrites build −1, −2 and typed zeros from them. A symbol name may not contain `*` (the product separator). |
| 44 (last sentence) | `Expr` (`core/base/expressions.py`) holds the operator dispatch (`__add__`, `__mul__`, … via `_add`/`_sub`/`_mul`), which picks the node class from `(a.type, b.type)` and raises `ExpressionMismatchError` naming the operands; `mixins.py` holds what node families share (`has`, `left`/`right`/`expr`, the product rule of `delta`/`t_diff`, the calculus identity). |
| 45 | **`core/transformations/`** — modules `expansion` (`expand`), `pulling` (`pull`), `vector_identities` (`vector_rules`, incl. S² tangent rules like `q×(ω×q)=ω`), `simplification` (`simplify` = pull + vector_rules + `_eliminate`, which dispatches per node class through `_ELIMINATE`; `full_simplify` = expand once, then simplify to a fixed point, warning at `max_iter`; like terms keyed by `term_key`), `manifold_rules` (substitutes `Variation(ω)` → `ξ̇ − 2ω×ξ` and `Variation(Ω)` → `η̇ + Ω×η`; derivation in the module docstring), and `utils` (`rebuild`, `map_children`, `addends`, `sum_of`, `flatten_mul`). Every rewrite returns through `rule(name, before, after)`; names are listed in each module docstring. |
| 46 | **`core/math/`** — `extract.py` (`extract_linear_coeff(expr, query)`: coefficient by *tree truncation*, 3×3 dispatch grid, raises `NonLinearError`/`AlgebraicError`, never returns zero for a term containing the query; also `has_algebraic`, `same_node`), `collecting.py` (`collect`: rotates `Dot`/triple products so the target lands on the left of a `Dot`, peels a matrix `Dot(M f(v), w) → Dot(f(v), Mᵀw)`, rewrites a row inner product to `Dot`, strips `Dot(R a, R b)→Dot(a,b)`), `ibp.py` (`integrate_by_parts`). |
| 47 | **`dynamics/`** — `compute_eom(L, dW, SystemVariables, params=(), trace=None)` validates the variables, then δL + δW → manifold rules → full_simplify → IBP against each variation vector's derivative → expand → `extract_linear_coeff` per variation vector, and raises if a variation symbol is left over (`docs/pipeline.md`). Returns `{name: (variation_vec, eq)}`. `to_standard_form` splits each equation into `M·a + f + G·u = 0` (`StandardFormEquation`); `linearize` (`vbl.py`) gives `M δẋ = A δx + B δu`. `geomech/numeric.py` holds the numeric oracle (`evaluate`, `check_eom`, `assert_eom`, `finite_difference_check`). |
| 48 | **`utils/`** — printing (`print_eom`, `print_tree`, `to_latex`, `display_*`, `render_eom`), `trace.py` (`EOMTrace`, `LinearizationTrace`, `rule`, `set_trace`, the `geomech` logger, `GEOMECH_DEBUG`/`GEOMECH_TRACE`), `invariants.py` (`check_tree`), and the error types. |
| 52 | `Expr.__eq__` compares class and `str()` (plus, for a calculus node or `Transpose`, the operand classes); `__hash__` is `hash(str)`; `has()` uses the same `==`. Products print `*` between factors and `*` is reserved in names, so two different products never print alike except for regroupings of one product. Commutative rearrangements are still *not* equal (`a*b` vs `b*a`); like-term merging uses `term_key`, which knows `Dot`, `Mul` and inner products commute. See `docs/conventions.md`. |
| 55 | `docs/manifold-implementation-notes.md` (in the docs nav) and the local, gitignored `paper_context.md` record the reasoning behind the S²/SO(3) substitution rules and the extraction algorithm; `docs/pipeline.md`, `docs/conventions.md` and `docs/troubleshooting.md` describe the pipeline, the identity rules and every typed error. Read them before changing those modules. |
| 59 (append) | Every pinned system also passes `check_eom` in `tests/test_eom_oracle.py` (≤ 1e-5); add new systems there too. Pinned strings print products with `*` (`m*l*l*(-1)`); a change that reprints them goes through `specs/reviews/scripts/pinned_snapshot.py` (snapshot before, compare after, numeric equality). Run the full suite, ruff and mypy before every commit. |
| 15 (project roadmap and earlier code review) | add: "The engine weaknesses in the earlier code review are superseded by `specs/fix-plan.md` and `specs/reviews/` (2026-11-13 audits and phases 1–5b)." |

### What a user must know before this stack is pushed

| change | in CHANGELOG |
|---|---|
| Printed form: products print `*` (`m*l`, `J*\Omega`, `(a)'*b`); every pinned string changes | yes, `CHANGELOG.md:186-190` |
| Equality: `Mul(ab,c) != Mul(a,bc)`; calculus nodes and `Transpose` compare operand class; `has()` uses `==` | yes, `:191-195`, `:200-201` |
| Name rule: a symbol name containing `*` raises `ValueError` | yes, `:196-199` |
| Frozen nodes: assignment raises `FrozenInstanceError`, `nodes` is a tuple, flags frozen; old pickles load | yes, `:202-205` (not stated: `attr` stays a mutable list) |
| Module renames: `transformations.{expand,pull,simplify,vector_rules}` and `core.math.collect` modules are `expansion`, `pulling`, `simplification`, `vector_identities`, `collecting`; old module paths raise `ModuleNotFoundError` | yes, `:251-263` |
| Removed: `is_leaf`, `has_nested_add`, `UndefinedCaseError`, `TS2/TSO3.delta(substitute=True)` | yes, `:403-409`, `:448-449` |
| Derived equations change value on S² when `∂L/∂ω` has a component along `q` | yes, under Fixed `:439-449` |
| New errors on input that used to pass: unknown `attr` (`:321`), stray variation symbols incl. a parameter missing `Constant` (`:344`), `SystemVariables` duplicates, tangent vectors, plain matrices, composite entries (`:353`, `:274`), `d/dt` of a variation in `dW` (`:328`), `TimeDerivative(q)` of a manifold leaf (`:370`), `collect` on non-scalar (`:373`), `extract_linear_coeff` raising instead of returning zero (`:464`), `t_integrate` of a constant (`:319`), `bool` operands, division by an expression (`:114`, `:341`) | yes, each |
| `full_simplify` warns (`RuntimeWarning`) at the iteration cap; becomes an exception under `-W error` | yes, `:377` (the `-W error` consequence is not spelled out) |
| `repr`, `print_tree` default style, `tree_str` labels, error-message texts changed | yes, `:292-307`, `:308` |
| Constant-symbol variations are `0v`/`0`, `display_standard_form` key spacing, literal spelling `(0.5)` | yes, `:208-212` |

CHANGELOG states every breaking change. Advisory: the `[Unreleased]` section is about 400
lines with the breaking items mixed into "Changed" and "Fixed"; a short "Upgrading from
0.1.0" list at its top (the rows above) would help, and `pyproject.toml` is still
`0.2.0.dev0`.

Other pre-push notes:

- The remote `vkotaru/pygeomech` is public. Commit author/committer on all 25 commits is the
  GitHub noreply address; no home path, host name or credential appears in the added
  lines. References in `specs/` to the author's local setup have been generalised
  (contributor notes, temporary directories, the accompanying paper, a downstream
  application).
- The paper fix for item 1.1 sits uncommitted in the accompanying paper.
- `paper_context.md` is cited as a design document but gitignored, so it will not travel
  with the push.
