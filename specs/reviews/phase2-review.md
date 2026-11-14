# Adversarial review of phase 2 (`specs/fix-plan.md` 2.1–2.12, "error instead of silence")

Branch `fix/phase2-validation`, range `865699f..c18fb55`. Reviewed: the contributor notes, the phase 2
table, `specs/reviews/phase2-implementation-report.md`, `git diff 865699f..c18fb55 --
geomech/ examples/ docs/` in full, `tests/test_validation.py`, `tests/test_simplify_numbers.py`.
Nothing under `geomech/` and no existing test was modified. New file:
`tests/test_phase2_adversarial.py` (56 passing, 13 strict `xfail` with `F<n>` in the reason).
Probe scripts, the 85-system dumps and the `865699f` worktree lived in a temporary directory.

Baseline on this machine: `python -m pytest tests/ -q -W error::RuntimeWarning` → 1630 passed,
5 xfailed (matches the report; no `RuntimeWarning` from 2.10 anywhere in the suite). With the
new file → **1686 passed, 18 xfailed**; `ruff format --check` + `ruff check` clean.

## Verdict: GO WITH CHANGES

Every guard in the table does what the report says, no legitimate system is rejected (85
pinned systems derive byte-identically to `865699f` through `compute_eom`,
`to_standard_form` and `linearize`), and the fuzz shifts are exactly the two predicted rows.
Two silent cases remain that the phase's own charter covers, and both are ten-line name
checks in `geomech/dynamics/eom.py`:

- **F1** dbg-A25 (`dW = Dot(ẋ.delta(), f)`) still derives `−ḟ` instead of `+f` without a
  message, on R, R³, S² and SO(3). The debuggability audit lists A25 next to A24 as "the
  real risk" and 2.2 cites that finding (dbg-F2), but the stray-variation walk only looks
  for *unlisted* variation symbols, and `d/dt(δx)` of a *listed* `x` passes it.
- **F2** a symbol printed like a listed variable's *tangent* vector (`Vector("\omega_{q}")`
  next to `S2("q")`) passes `_check_variables` and `_check_params`, which only guard the
  *variation* vector names; as a param it produces an equation the oracle rejects (error 1.97).

F3–F5 are message quality. Nothing else moved.

---

## Findings, ranked

Severity order: silent misuse still accepted > legitimate system now rejected > wrong
message > cosmetic. **No legitimate system is rejected.**

### F1 — a time derivative of a variation vector in δW is accepted silently (silent misuse; dbg-A25)

**Reproducer** (`tests/test_phase2_adversarial.py::TestF1DerivativeOfVariationInVirtualWork`,
five strict xfails plus the pinned symptom):

```python
L, dW, x, f, *_ = point_mass()                       # L = ½ m ẋ·ẋ − m g e₃·x
compute_eom(L, Dot(x.t_diff().delta(), f), SystemVariables(vectors=[x]))
# {'δx': (d/dt(ẋ)·(−m) + e3·(−m g) + d/dt(f)·(−1))}        ← −ḟ, no message
compute_eom(L_s2, Dot(xi.t_diff(), f), SystemVariables(vectors=[q]))     # … + d/dt(f)(−1)
compute_eom(L_s2, Dot(q.t_diff().delta(), f), SystemVariables(vectors=[q]))
# … + Cross(Cross(ω,q),f)(−1) + Cross(q,ḟ)(−1) + Cross(ω,Cross(q,f))(2) + Cross(q,Cross(ω,f))(−1)
compute_eom(L_so3, Dot(eta.t_diff(), M), SystemVariables(matrices=[R]))  # … + d/dt(M)(−1)
compute_eom(L_r, x.t_diff().delta() * f, SystemVariables(scalars=[x]))   # … + d/dt(f)(−1)
```

**Observed vs expected.** At c18fb55 and at 865699f all five return a plausible equation
whose force term is the negated time derivative of the force (or, for `δ(q̇)·f`, four
cross-product terms). `check_eom` cannot be used as a backstop: it raises
`KeyError: "no value for the perturbation symbol 'd/dt(δx)' in env"`. Expected: `ValueError`
naming `inf_work`, the derivative and the variation vector, e.g. `compute_eom: inf_work
contains 'd/dt(δx)', the time derivative of the variation vector 'δx'; virtual work is linear
in the variation vectors themselves (write Dot(x.delta(), f), not Dot(x.t_diff().delta(), f))`.

**Diagnosis.** `geomech/dynamics/eom.py:165-199` (`_stray_variation_symbols`) walks
`TimeDerivative(Variation(x))`, finds `Variation(x)` with `str(node) in allowed` and returns.
Nothing inspects `inf_work` on its own: δW is added to δL at `:245` and the IBP at `:262-269`
then happily moves the derivative onto `f`. The check belongs before the IBP and after the
manifold rules have been applied to `inf_work` alone (the `δ(q̇)·f` case only exposes
`d/dt(ξ)` once `δω → ξ̇ − 2ω×ξ` has fired): any `TimeDerivative` node in δW whose subtree
contains a listed variation vector (a `Variation` in `allowed`, or a `TS2`/`TSO3` leaf named
like a listed point's variation vector) is an error. A velocity-dependent *force*
(`Dot(x.delta(), −c ẋ)`, damping) must still pass; it is pinned in the same class.

### F2 — a symbol printed like a listed variable's tangent vector passes the name checks (silent misuse; name collision)

**Reproducer** (`::TestF2TangentVectorNameClash`, three strict xfails plus the pinned symptom):

```python
L, q, om, xi, f = spherical_pendulum()
w = Vector(q.manifold_info.tangent_vector_name)               # a plain Vector "\omega_{q}"
eom = compute_eom(L + Dot(w, om), Dot(xi, f), SystemVariables(vectors=[q]), params=[w])
# ξ: d/dt(ω_q)(−1) m l l + Cross(q,e3) m g l (−1) + d/dt(ω_q)(−1) + Cross(ω_q, ω_q)(2) + f
check_eom(..., params=[w])  → 1.97          # the same L with w = Vector("w") → 7.7e-9
compute_eom(L + ½ ẇ·ẇ, Dot(xi, f) + Dot(w.delta(), f), SystemVariables(vectors=[q, w]))
# accepted; keys 'ξ_q' and 'δ(ω_q)', the latter is also the printed form of Variation(ω_q)
```

**Observed vs expected.** `_check_params` accepts `w` (it is not a configuration variable
by name), `_check_variables` accepts it in `vectors` (the clash set holds only `ξ_q`/`η_R`).
The tree keeps two different leaves (`Vector` vs `TS2`) with one printed name, so every
str-keyed consumer merges them: the printed equation shows `Cross(ω_q, ω_q)` where
`Cross(ω_q, w)` is meant, `check_eom` evaluates both to one value, and `linearize` would
bucket them together. Expected: `ValueError` from both checks, in the style of the existing
variation-vector clash: `compute_eom: params entry '\omega_{q}' has the printed name of the
tangent vector of 'q'; …`.

**Diagnosis.** `geomech/dynamics/eom.py:76-79` builds `seen_variations` from
`get_variation_vector()` only; `:120-162` (`_check_params`) compares against
`config_names` only. Add the tangent-vector names (`get_tangent_vector()` for `S2`/`SO3`
entries; `str(s.t_diff())` is not a leaf name and needs no entry) to both clash sets. The
general version of this (any leaf anywhere in `L` whose printed name is a reserved symbol)
is what phase 4.4's invariant checker is for; the two entry points above are enough for
phase 2 because they are the only places a user can *declare* such a symbol.

### F3 — input-shape mistakes surface as Python's own message (wrong message)

`::TestF3InputShape` (three strict xfails):

| call | observed | expected |
|---|---|---|
| `compute_eom(..., params=l)` (bare leaf) | `TypeError: 'Scalar' object is not iterable` | `TypeError: compute_eom: params must be a sequence of leaves, got Scalar 'l'; write params=[l]` |
| `SystemVariables(vectors=x)` | `TypeError: 'Vector' object is not iterable` (from `_check_variables`'s `for var in group`) | message naming `SystemVariables.vectors` |
| `compute_eom(m*m, None, SystemVariables())` | returns `{}` | `ValueError: compute_eom: no configuration variables listed` (with a real `L` the stray check fires correctly, pinned) |

**Diagnosis.** `eom.py:134` iterates `params` without a type check; `:44` iterates each
group likewise; `separate_variations` with an empty list returns `{}`.

### F4 — a hand-built `Variation(TimeDerivative(x))` is reported as a stray `d/dt(x)` (wrong message)

`::TestF4VariationOfDerivativeMessage`: `compute_eom(L, dW + Dot(Variation(x.t_diff()), f),
SystemVariables(vectors=[x]))` raises (good) with `… variation symbols 'δ(d/dt(x))' of the
symbols 'd/dt(x)' … list it in SystemVariables …`. Listing `d/dt(x)` is exactly what
`_check_variables` rejects two lines earlier. Reachable only by hand (the engine builds
`d/dt(δx)`), so low. **Diagnosis:** `eom.py:180-183` reports `str(inner)` without peeling
calculus nodes; the message should name `x` and point to `x.t_diff().delta()`.

### F5 — `attr="Constant"` raises about the attribute `'C'` (wrong message)

`::TestF5AttrAsString`: `Scalar("m", attr="Constant")` → `ValueError: unknown attribute 'C'
for symbol 'm'; valid attributes are […]`. It never silently accepts (no single character is
valid, pinned), but the cause — a string where a list is expected — is not named.
**Diagnosis:** `expressions.py:195` does `list(attr)` on a str; a one-line
`isinstance(attr, str)` check in `ExprFlags.from_attr_list` (`types.py:45`) would say
`attr must be a list of attribute names, got the string 'Constant'; write attr=['Constant']`.

### Cosmetic (no test)

- `_check_variables` duplicate message reads `listed more than once (vectors and vectors)`
  when both entries are in the same group (`eom.py:81-84`).
- `bool` is an `int`: `Number(True)` prints `(True)` with the `Ones` flag, `x * True` prints
  `x(True)`, `compute_eom(L, False, …)` means δW = 0 (`expressions.py:526`, `eom.py:236`).
  Pre-existing shape, now with flags; pinned as a note in `TestCapWarningAndNumbers`.
- 2.10 warns when the *last* pass produced the fixed point but no confirming pass ran
  (`max_iter=2` on a two-pass expression returns `0` and warns "no fixed point after 2
  passes"). Correct per the loop's contract and the message says "may be"; pinned.
- Message prefixes: `compute_eom:`, `collect:`, `apply_manifold_rules:`, `full_simplify:`
  follow `linearize:`; `SystemVariables.<group> entry …` (pre-existing), `TimeDerivative(q)
  is not represented …`, `Transpose() needs an operand`, `unknown attribute …` carry none.
- `TestStrayManifoldVariation` in `tests/test_validation.py` also holds the R³/scalar cases
  (`y`, `e3`, `l`); the name undersells it.
- CHANGELOG, 2.3 bullet: "a plain matrix reached `NotImplementedError`". At `865699f` both
  shapes I tried (`Matrix("J")` in `matrices` with and without a `J̇` term in `L`) raise
  `AlgebraicError` from `extract`/`collect`, not `NotImplementedError`. If the claim comes
  from a `linearize` path, say so; otherwise it reads as the symptom of a different case.

---

## What was checked and held

### 2.2 — stray variations and `params` (23 passing tests)

Legitimate systems that derive and pass `check_eom(params=…)` at `< 1e-5`: an input in
both `L` and δW (`params=[f]`; without it the error names `δf`, pinned); a vector and a
scalar param that appear **only under `TimeDerivative`** (`δ(ḃ)` is built as `d/dt(δb)`, so
the walk sees `Variation(b)` with `b` declared — no false positive); a prescribed inertia
`J(t)` (`J̇` in the equation, `to_standard_form` carries `J` in `M`); one `l` shared by two
S² bodies; a param under `Hat` and inside an outer product `b cᵀ` (both go through
`_rebuild`); a pendulum hung from a prescribed point `x = b + l q` through `compute_eom`,
`check_eom`, `to_standard_form`, and `linearize`'s documented `AlgebraicError` on `δb`; a
non-constant force only in δW without `params`; quadrotor thrust `δx·(R e₃ u)` through the
whole pipeline; `params` as tuple, list, generator and set give one printed equation. All
nine `vbl_systems` run `compute_eom → to_standard_form → linearize`.

Misuse that raises with the symbol and the three remedies: unlisted `Vector` (dbg-A24),
unlisted `Scalar`, an `S2` point that appears only statically (`ξ_p`), an `SO3` point that
appears only through `R ρ` (`η_R`, `δΩ_R`), `S2("q")` in `L` with `Vector("q")` listed and
the reverse, three strays at once (all named, sorted), a param that is a variable, a
constant param. Cancelling stray terms (`Dot(x,y) − Dot(y,x)`) are not reported.

Also seen, not findings: `Scalar("x")` in `L` with `Vector("x")` listed raises
`AlgebraicError: _vector_from_scalar: unsupported node Variation` — not silent; the message
is the str-equality limitation of D2 (open decision, out of phase 2). `Dot(b.delta(), f)` in
δW with `b` in `params` is dropped — consistent with the declaration.

### 2.3 — `SystemVariables` (5 passing)

`SO3` in `vectors`, `S2` in `matrices`/`scalars`, `None`/`str`/`int` entries, an
`OrthogonalMatrix`-flagged `Matrix`, two distinct objects with one name, `Vector("q")` +
`S2("q")`: all raise naming the entry; `linearize` runs the same checks. Report's 12 tests
re-read and hold.

### 2.1 — `attr` (3 passing)

Every attr string in the repo (`examples/`, `docs/`, `README.md`, `tests/`, `geomech/`,
`notebooks/*.ipynb`): `Constant` ×423, `SymmetricMatrix` ×46, `Zero` ×35, `UnitNorm` ×9,
`Ones` ×9, `TangentVector`, `SkewSymmetricMatrix`, `OrthogonalMatrix`, `Identity` once each —
all in `VALID_ATTRS`. The six lowercase/abbreviated strings found (`constant`, `unitnorm`,
`skew`, `Symmetric`, `Orthogonal`, `Const`) are the negative cases of
`tests/test_validation.py` itself. Case sensitivity: the error lists the valid names
verbatim and CHANGELOG spells them; `docs/getting-started.md` shows only `Constant`
examples and does not say the names are case-sensitive (a half-sentence in the attributes
table at `:63` would close it). `None`, `1`, `""`, `" Constant"`, `"CONSTANT"` all raise.
Tuples and sets are accepted (iterated), duplicates are harmless.

### 2.5 — parent by name (5 passing)

Two `S2("q")` objects mixed in one `L` (`ω` from one, `q` from the other, `q̇·q̇` across
both) give the same printed equation as one object **and pass `check_eom`**;
`S2("q")`/`S2("q1")` never reduce; `S2("q", attr=["Constant"])` is the same point as
`S2("q")` (by name, as decided in D4); a plain `Vector("q")` and a parentless `TS2` are not
reduced. No `is` comparison on a manifold parent is left in `geomech/` (grep).

### 2.8 — calculus guards (2 passing)

`TimeDerivative(ω)`, `Variation(ω)`, `TimeDerivative(ξ)`, nested `d/dt`, `δ(d/dt x)`,
`d/dt(δx)`, `(l q).t_diff()`, `(A q).t_diff()`, `Hat(q).t_diff()`, `q.T()`, `q.t_integrate()`
all build; `TimeDerivative(S2)`/`Variation(SO3)` point to `.t_diff()`/`.delta()`.

### 2.9 — `collect` order (3 passing; fuzz; 85-system dump)

`engine_fuzz.py --quick` and `--cases 2000 --seed 99`, run at `c18fb55` and in a `865699f`
worktree: the only rows that move are `collect[target occurs several times]` (ok 125 → 109,
`NonLinearError` 155 → 171 at 2000 cases; 14 → 10 / 23 → 27 quick) and `collect[target
inside matrix factor]` (ok 2 → 0; 1 → 0 quick); `collect[target occurs once]` is 1222/1222
in both; 0 wrong-value cases in both; no new exception class. The one behavioural change on
a single `Dot` is `Dot(v, v)` (base: returned as "canonical"; now `NonLinearError`), which
is quadratic and not pipeline-reachable; `collect`'s only caller is `ibp.py:173`.

All 85 pinned systems (36 `test_eom_oracle.SYSTEMS`, 9 `vbl_systems`, 40 `phase1_systems`)
dumped through `compute_eom → to_standard_form → linearize` at both commits: **byte-identical**
(including the three pre-existing typed errors: `s2/so3_outer_product_on_velocity` raise
`AlgebraicError` from `collect`, `adv_virtual_work_under_rotation` from `linearize`).

### 2.10 / 2.11 (6 passing)

Warning at `max_iter=1` (unsimplified) and at `max_iter=2` (fixed point reached but not
confirmed, see cosmetic), silent at 3; the full suite under `-W error::RuntimeWarning`
emits none. `Number(0)` in like-term merging (`a + 0·a + a → 2a`, `0 + a + 0 → a`,
`2 + (−2) + a → a`); `Number(1)` as a genuine coefficient (`1·a + a → 2a`, `v·1 + v → 2v`,
`1 + 1 → (2)`, `1·1 → (1)`, `extract_linear_coeff(Dot(v,w)·1, v) → w(1)`); `-0.0`/`0.0`
fold to the Zero leaf, `1.0` drops; `1e-300` is not zero; `Scalar("m", value=0)` keeps its
name (D1 open).

### 2.6 / 2.7 / 2.12 / exports

`None`, `0`, `0.0`, `False`, `Number(0)` mean δW = 0; `1`, `"0"`, `ZeroVector`, `Dot(f, f)`
(scalar, no variation: accepted, adds nothing) behave as documented; the `_Opaque` test
covers 2.7 and every node class the engine defines is in `_rebuild`'s three tuples;
`examples/rigid_pendulum.py` prints `J1`/`J2`; the four error classes are the engine's.

### Debuggability

No new `except`, no `assert` used for validation in the diff (grep). `compute_eom`'s
docstring and `docs/getting-started.md` state the three-way contract (variable / constant /
`params`). CHANGELOG "before" claims verified at `865699f`: `None`/`0` →
`AttributeError`, vector δW → `ExpressionMismatchError: Add`, tangent vector listed →
`ω m`-style equation, `Number(0)*x → (0)x` with `is_zero False`, attr typo ignored,
duplicate collapsed to one key, unlisted `y` dropped — all as written, except the
`NotImplementedError` line (cosmetic above). Test names in `tests/test_validation.py` read
as sentences.

---

## Required before merge

1. **F1**: reject a `TimeDerivative` of a listed variation vector in `inf_work` (after
   manifold rules on δW, before IBP), message naming the derivative and the vector. The five
   xfails in `TestF1DerivativeOfVariationInVirtualWork` must pass; the damping test next
   to them must keep passing.
2. **F2**: add the tangent-vector names of listed `S2`/`SO3` entries to the clash sets of
   `_check_variables` and `_check_params`. The three xfails in `TestF2TangentVectorNameClash`
   must pass.

Recommended, same commit: F3 (two `isinstance` checks and an empty-variables check), F5
(one `isinstance(attr, str)` branch), F4 (peel calculus nodes when naming the stray), the
`(vectors and vectors)` wording, and the CHANGELOG `NotImplementedError` line. Flip the
corresponding xfails to passing tests rather than deleting them.

---

## Fix log

Working tree on `fix/phase2-validation` at `ca4baf6` (not committed). All 13 strict xfails in
`tests/test_phase2_adversarial.py` had their marker removed with the assertions unchanged
and pass. The pipeline that builds δS is untouched: the new δW check runs on a separate
copy of `inf_work`, so all 85 pinned systems (`test_eom_oracle.SYSTEMS`, `vbl_systems`,
`phase1_systems`) dump byte-identically through `compute_eom → linearize` against
`ca4baf6` (same three pre-existing typed errors).

### F1 — derivative of a variation vector in δW (`geomech/dynamics/eom.py`)

- `:267` `_derivatives_of_variations(expr, allowed)`: outermost `TimeDerivative` whose
  subtree holds a listed variation vector (`Variation` or `TS2`/`TSO3` leaf with its printed
  name in `allowed`).
- `:297` `_check_virtual_work`: `full_simplify(apply_manifold_rules(inf_work, variables))`,
  then the walk above. Called at `:385`, after the variation vectors are gathered and before
  δL is formed; δL is never passed to it (ξ̇/η̇ are legitimate there).
- Message (`:314`), e.g. R³:
  `compute_eom: inf_work contains '\frac{d}{dt}(\delta{x})', the time derivative of the
  variation vector '\delta{x}'; virtual work must be linear in the variations themselves,
  not their time derivatives (write Dot(x.delta(), f), not Dot(x.t_diff().delta(), f))`.
  For `δ(q̇)·f` on S², where the derivative only appears once `δω → ξ̇ − …` has fired, the
  term is followed by ` after the manifold kinematic substitutions`; more than one bad term
  adds ` (and N more)`.
- Damping `δx·(−c ẋ)` and every pinned system pass.

### F2 — tangent-vector name clash (`geomech/dynamics/eom.py`)

- `:25` `_tangent_names`; `_check_variables` (`:107`):
  `compute_eom: SystemVariables.vectors entry '\omega_{q}' has the printed name of the
  tangent vector of 'q'; configuration variables must be distinct symbols (rename it)`
  (order-independent: the tangent names of all listed points are collected first).
- `_check_params` (`:165`, `:204`) also reserves the variation-vector names:
  `compute_eom: params entry '\omega_{q}' has the printed name of the tangent vector of 'q';
  rename the parameter (two different symbols with one printed name are merged when the
  equation is printed or evaluated)`.

### F3 — input shape

- `geomech/dynamics/variables.py:8` `check_group_is_sequence`, run in
  `SystemVariables.__post_init__` (`:39`) and again in `_check_variables` (prefixed):
  `SystemVariables.vectors must be a list of configuration variables, got Vector 'x';
  write vectors=[x]` (list or tuple accepted).
- `eom.py:158`: `compute_eom: params must be a sequence of leaves, got Scalar 'l'; write
  params=[l]` (list, tuple, set, generator still accepted).
- Empty `SystemVariables()`: decided to keep it constructible (`test_eom.py` pins the empty
  defaults, `check_eom(..., SystemVariables(), {})` returns 0 in `test_numeric.py`) and to
  reject it in `compute_eom` *after* the stray check (`eom.py:434`), so a real Lagrangian
  still gets the more useful stray message naming `x` (pinned by
  `test_empty_variables_with_a_real_lagrangian_raises_stray`):
  `compute_eom: no configuration variables listed in SystemVariables; a system with no
  degree of freedom has no equation of motion`. Documented in the `SystemVariables` and
  `compute_eom` docstrings.

### F4 — hand-built `Variation(x.t_diff())` (`eom.py:244`, `:426`)

`_stray_variation_symbols` peels `TimeDerivative` under a `Variation` and records the leaf;
the message keeps the generic part (pinned by `test_it_does_raise`) and, for such symbols,
replaces the "list it in SystemVariables" remedy with:
`…; '\delta{\frac{d}{dt}(x)}' is the variation of a time derivative, which the engine never
builds: write the variation of the configuration variable itself, x.delta() (the variation
of the velocity is x.t_diff().delta(), built as d/dt of the variation, and is not allowed
in inf_work)`. The generic remedy is still appended when another, ordinary stray is present.

### F5 — `attr` as a string

`geomech/core/base/expressions.py:16` `_attr_list` (used by all six leaf constructors) and
`geomech/core/base/types.py:52` (`ExprFlags.from_attr_list`, for direct callers):
`TypeError: attr must be a list of strings, got str 'Constant' for symbol 'm'; write
attr=['Constant']`.

### Cosmetic

- Duplicate wording (`eom.py:97`): `… listed more than once (twice in
  SystemVariables.vectors)` / `(in SystemVariables.scalars and vectors)`.
- bool: `Number(True)` (`expressions.py:540`) → `TypeError: Number expects int, float, or
  str, got bool True; write Number(1) or Number(0)` (also reached by `x * True`);
  `compute_eom(L, False, …)` (`eom.py:356`) → `TypeError: compute_eom: inf_work must be a
  scalar expression, None or 0, got bool False`.
- 2.10 warning: **no code change.** The loop already breaks on the first pass that changes
  nothing, so reaching the cap means the last pass *did* change the expression; the
  `max_iter=2` case changes `…` to `0` on pass 2 and warns correctly by that criterion.
  Not warning there would need an extra confirming pass beyond `max_iter`, which breaks the
  cap's meaning and the pinned `test_warning_fires_only_when_the_loop_exhausts`. The
  docstring (`simplify.py:56`) now says this explicitly.
- Prefixes: `_check_variables` takes `caller` (`compute_eom` by default, `linearize` from
  `vbl.py:216`) and every message starts `compute_eom: …`/`linearize: …`;
  `TimeDerivative: TimeDerivative(q) is not represented …` (`calculus.py:19`, likewise
  `Variation:`); `Transpose: Transpose() needs an operand; write Transpose(expr)`
  (`geometry.py:126`); `attr: unknown attribute 'constant' for symbol 'm'; …`
  (`types.py:88`). The bare-`SystemVariables(vectors=x)` message from the constructor
  carries the class name as its prefix.
- CHANGELOG: the 2.3 line now reads "a plain matrix failed later with an `AlgebraicError`
  from `extract_linear_coeff`/`collect` that did not name the entry"; Changed gained
  bullets for F1–F5 and bool.
- `docs/getting-started.md`: after the expression-types table, "Attributes are passed as a
  list of strings, `attr=['Constant']`, and the names are case-sensitive (`'constant'`
  raises `ValueError` listing the valid names)."
- Not done: renaming `TestStrayManifoldVariation` in `tests/test_validation.py` (existing
  test file, cosmetic only).

### Pinned tests changed (they pinned the defects being fixed)

Four non-xfail tests in `tests/test_phase2_adversarial.py` asserted the old behaviour and
could not survive the fixes; each was rewritten to pin the new behaviour, docstring noting
the old one:
`test_observed_symptom_is_a_derivative_of_the_force` → `test_message_names_the_derivative_and_the_fix`
(`:125`, the bad δW now raises; the correct δW still ends in `+f`);
`test_observed_symptom_param_case_fails_the_oracle` → `test_renamed_param_passes_the_oracle`
(`:183`, the clashing param now raises; the renamed one passes `check_eom`);
`test_string_attr_does_raise` (`:281`, now expects `TypeError` "attr must be a list of
strings, got str" instead of `ValueError` "unknown attribute");
`test_bool_literals_are_ints` → `test_bool_literals_are_rejected` (`:705`).

### Final counts

- `ruff format geomech/ tests/` (76 files unchanged) + `ruff check` — clean.
- `python -m pytest tests/ -q -W error::RuntimeWarning` — **1699 passed, 5 xfailed** (the
  five pre-existing xfails; 1686 + the 13 flipped).
- `mkdocs build --strict` — builds (INFO only).
- `specs/reviews/scripts/engine_fuzz.py --quick` — 0 wrong-value cases; rows as at
  `c18fb55` (`collect[target occurs several times]` 10 ok / 27 `NonLinearError`,
  `collect[target inside matrix factor]` 0 ok).
- `examples/*.py` — all six run.
