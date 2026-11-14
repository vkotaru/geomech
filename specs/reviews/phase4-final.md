# Phase 4 final review: debuggability

Branch `fix/phase4-debuggability`, range `57232da..86e7f26` (implementation `bffba73`,
adversarial review `6239da7`, fixes `86e7f26`). Read: `specs/fix-plan.md` phase 4,
`phase4-implementation-report.md`, `phase4-review.md` (hunts, fix log, arbitration), and
`git diff 57232da..86e7f26 -- geomech/ docs/ README.md CHANGELOG.md` in full. Nothing in
`geomech/` or `tests/` was changed; nothing is committed.

## Verdict: GO WITH CHANGES

The code is sound. Output is byte-identical to `57232da` except for the two expected error
prefixes. Fuzz findings are identical, there is no runtime cost, and every review finding
was fixed or arbitrated. A fresh blind hunt located a new kind of bug in about 40 s with
the trace alone. The four required changes are small text fixes. One of them matters for
trust: run from a second checkout, step 2 of the new "no row matches" recipe silently
tests a different tree and reports a false all-clear (R1).

## Required changes

**R1. `specs/reviews/scripts/engine_fuzz.py` fuzzes whichever `geomech` is installed, not
the checkout it sits in.** `REPO` is derived from `_geomech.__file__` (`engine_fuzz.py:1800`),
and `geomech` is imported before anything puts the script's repository on `sys.path`. In a
worktree or a copy that shares the main venv (this project's normal review workflow),
`python specs/reviews/scripts/engine_fuzz.py --quick`, the exact command in
`docs/troubleshooting.md:37-38` and `docs/pipeline.md:245`, fuzzes the editable install. It
reported 0 `[wrong-value]` on my bugged copy. Fix: before the first `geomech` import, add

```python
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
```

(the repository root, three levels above the script), and print `geomech.__file__` in the
header line. Alternatively, have the docs say to run it with `PYTHONPATH=.`. The script also
needs `tests/numeric_helpers.py` for its self-test (`:1756`), so a copy without `tests/`
crashes. That one is fine, because it is loud.

**R2. `docs/pipeline.md:191-197`: orphaned paragraph.** The F4 note was inserted between the
code block and its continuation, so a paragraph now starts with a lowercase "or run with
`GEOMECH_DEBUG=1` …" that follows the note about handlers. Move the line
"or run with `GEOMECH_DEBUG=1`, which also runs `check_tree` on the inputs and outputs of
`compute_eom`, `to_standard_form` and `linearize`." up so that it directly follows the
```` ``` ```` that closes the `logging` block. Put the paragraph
"`logging.basicConfig(level=logging.DEBUG)` works as well. Configure logging before …" after it.

**R3. `README.md:84-85`: the `_eliminate` snippet still shows the pre-F3 fold.** It reads
`return rule("simplify.add_fold", expr, _simplify_add([_eliminate(n) for n in nodes]))`,
which is the form F3 removed, because it logs the fold for a child's rewrite. Replace it with:

```python
        case Add(nodes=nodes):               # n-ary
            children = [_eliminate(n) for n in nodes]
            out = _simplify_add(children)
            if _trace._TRACE_RULES:          # a fold logs only its own work
                rule("simplify.add_fold", _folded_from(Add, children, expr), out)
            return out
```

and change `README.md:88` "Every rewrite returns through `rule(name, before, after)`" to
"Every rewrite is logged through `rule(name, before, after)`".

**R4. CHANGELOG: the AlgebraicError count and the import path that breaks.**
- `CHANGELOG.md:56`: "the two `AlgebraicError` causes told apart" → "the three
  `AlgebraicError` causes told apart". `troubleshooting.md:84-90` has three: surviving
  d/dt, jerk, and a Lagrangian in accelerations.
- `CHANGELOG.md:159-160`: the example given is a private name. The public form that breaks
  is the one ten test files used at base. Replace "Code that imported from the old module
  paths (`from geomech.core.transformations.simplify import _eliminate`) must use the new
  ones." with "Imports through the old *module* paths now raise `ModuleNotFoundError`
  (`from geomech.core.transformations.simplify import simplify`, `from
  geomech.core.math.collect import collect`, …); import the functions from `geomech`,
  `geomech.core.transformations` or `geomech.core.math`, which are unchanged, or use the
  new module names."

## Verified

| Check | Command | Result |
|---|---|---|
| Suite | `python -m pytest tests/ -q -W error` | **2120 passed, 5 xfailed**, 8.8 s |
| Suite, debug | `GEOMECH_DEBUG=1 … -W error` | 2120 passed, 5 xfailed, 9.7 s |
| Suite, rule trace | `GEOMECH_TRACE=1 … -W error` | 2120 passed, 5 xfailed, 11.1 s |
| Lint | `ruff format --check geomech/ tests/ && ruff check geomech/ tests/` | 94 files formatted, all checks passed |
| Docs | `mkdocs build --strict` | exit 0; only INFO is the existing `reference/` link; the API pages are generated for `expansion`, `pulling`, `simplification`, `vector_identities`, `collecting`, `utils/trace`, `utils/invariants` |
| Fuzz | `engine_fuzz.py --cases 2000 --seed 8080 -W error`, head and base | exit 0 both; **159 findings each, the same list** (88 flag-dependent, 66 unsupported-linear-term, 5 documented-error), **0 wrong-value**; texts identical apart from `where:` |
| Examples (6) | `python -W error examples/X.py`, head vs base | all exit 0, stdout+stderr byte-identical |
| Notebooks (4) | `jupyter nbconvert --execute` into scratch, compare code-cell outputs to committed | 18/28/12/9 outputs, all identical |
| `str` dump, 65 systems | my own script: every argument-free builder of `tests/vbl_systems.py`, `tests/phase1_systems.py`, `tests/test_eom_oracle.py` (minus the 2 new row systems): `L`, `dW`, δS, `apply_manifold_rules`, `full_simplify`, `integrate_by_parts`, `expand`, `compute_eom`, `to_standard_form` (M, f, G), `linearize` (rows, states, dstates, inputs, M, A, B) | 347 278 B at base; **2 differing lines**, the expected `compute_eom: while integrating by parts (stage ibp, …):` prefix on `s2_/so3_outer_product_on_velocity` |
| Timing | `compute_eom`, median of 40, 3 alternating rounds, tracing off | quadrotor with load: base 11.31 / 10.92 / 10.87 ms, head 10.77 / 10.86 / 11.20 ms (noise); double rigid pendulum: base 18.48 / 18.11 / 17.78 ms, head 17.61 / 17.56 / 17.93 ms |
| Logging | `GEOMECH_DEBUG=1` with and without `basicConfig(level=DEBUG)` first | 23 lines both ways, no duplicates |

The base worktree (`57232da`, under a temporary directory) and all scratch copies were removed;
`git worktree list` shows only the main checkout.

### Docs read as a newcomer (item 3)

- **Snippets run.** `pipeline.md`: the `EOMTrace` block, `show("simplify")`,
  `tr.ibp[...]["collect"]`, `show("ibp.collect[…]")`, `tr["ibp"] is tr.after_ibp`, and trace
  reuse (13 stages after a second run, keys of the second run only). The
  `summary(chars=70)` block for the rigid pendulum reproduces exactly. The rule log is
  **38 lines with exactly the per-rule counts listed** (`manifold.so3_variation_Omega` 2,
  `simplify.mul_fold` 4, …). `troubleshooting.md` recipe step 3
  (`collect(tr["ibp.prepare"], t)`, `full_simplify(tr["manifold_rules"])` under `set_trace`)
  works. `check_eom(L, dW, variables, eom)` matches the real signature. A failed run
  (`½ ẍ·ẍ`) gives `failed_at = raw[\delta{X}]` and a summary that starts with
  `FAILED at …`; `SystemVariables()` with a real `L` gives `failed_at = stray_check`, as the
  doc row now says. `test_troubleshooting_doc.py` covers the typed-error rows.
- **`conventions.md` claims checked:** `(v)'s` vs `(vs)'`, zero leaves `0`/`0v`/`0`/`(0)`/`(0v)'`,
  constant vector δ and d/dt printing `0`, constant matrix δ `O` and d/dt `0`, `0 + x = x`,
  `m + 0 → (m+(0))`, `x + 0` raising for a vector, `x(-1)`/`xm`,
  `Mul(-0.5, Mul(m,2)) → m(-1.0)`, `x*0.1*3 − x*0.3 → 0` (prints `0v` for a vector `x`),
  `Scalar("k", Constant) == Scalar("k")`, `check_tree` rejecting `Scalar("(-1)")`,
  `Scalar("m", value=2)` folding (`m2·3 → (6)`). Point mass: `E = −mẍ + F`, `M = I(-1)m`,
  `G = I`.
- **Consistency.** README, `getting-started.md`, `index.md`, the `compute_eom` docstring
  and `pipeline.md` state the same stage order as `eom.py`. All of them use the new module
  names; no doc names `vector_rules.py`, `simplify.py`, `expand.py`, `pull.py` or
  `collect.py`. The two inconsistencies are R2 and R4.

### Module rename (item 4)

- Old module paths, `from geomech.core.transformations.{simplify,expand,pull,vector_rules} import …`
  and `from geomech.core.math.collect import …`, now raise `ModuleNotFoundError` (R4 asks
  the CHANGELOG to say so plainly).
- Package-level imports are unchanged and work: from `geomech`,
  `geomech.core.transformations` (incl. `apply_manifold_rules`) and `geomech.core.math`.
  They are plain functions, `pickle.dumps(simplify)` works, and `geomech.simplify is
  geomech.core.transformations.simplify`.
- Grep across README, docs, examples, notebooks, a downstream application and the accompanying paper
  (grep only): **no old module path in use.** Examples and notebooks import from `geomech`,
  `geomech.utils.printing` and `geomech.core.operations.multiplication`, all unchanged. The
  paper names functions only (`apply_manifold_rules()`, `full_simplify()`,
  `integrate_by_parts()`, `separate_variations()`, `expand`, `simplify`) plus the helper
  `_is_tangent_to_s2` (`03_toolbox.tex:176`), which still exists in `vector_identities.py:124`.
  The only stale path is `specs/fix-plan.md:55` (`collect.py`), a historical spec; it is fine.
- The 10 test files that imported old module paths had only their import lines changed
  (checked line by line); `test_rule_names.py::test_submodule_is_the_module` was rewritten
  by author arbitration.

## Blind bug hunt (item 2)

**Setup.** A separate agent copied the tree to a temporary directory, injected one silent bug of
its choice (excluded: the three bugs of the review), wrote `scenario.py`, and sealed the
answer in a file outside the copy. I did not see the bug until I had found it.

**Scenario.** An elastic spherical pendulum on R × S², `L = ½ m |d/dt(l q)|² − m g l q·e₃ −
½ k (l − l₀)²`. `check_eom` = **0.607** (correct code: 1.0e-7).

**Steps (wall clock from the first command):**
1. Ran the scenario. Reading the ξ equation against the physics, the Coriolis terms
   `−l l̇ m ω − l̇ l m ω + l l̇ m ω` sum to `−m l l̇ ω`; Lee et al. needs `−2 m l l̇ ω`.
2. Recipe step 2, the fuzzer: **0 wrong-value. A false negative (R1).** It had fuzzed the
   real tree through the editable install.
3. Recipe step 1, bisect: `L = ½ m l² ω·ω` alone, then `tr.summary()` and
   `tr.ibp[t][step]`. Every stage up to `ibp.collect[ξ̇]` was fine. `ibp.integrate[ξ̇]`
   produced `Dot(ξ(-1), ω l l̇ m (1))`, where the product rule on the prefactor `l·l·m`
   should give `2 l l̇ m`. **Stage localised.**
4. Recipe step 3, `set_trace(True)` on `_apply_ibp(Dot(ξ̇, w)·pre, ξ̇)` for three
   prefactors. The log line `ibp.dot_prefactor: … ==> … ·l·l̇·m*·(1)` named the rule, and
   `a·l` and `l·a·m` came out right, so the bug needs a *repeated* factor.
5. Read `_t_diff_product` (`ibp.py:102`, called by `ibp.dot_prefactor`): `Mul(term, count - 1)`
   where its own docstring says `nᵢ`. **Found in ≈ 40 s, 4 tool steps.** The sealed answer
   matched exactly.
6. Rerun with `PYTHONPATH=<copy>` and `tests/numeric_helpers.py` added: the fuzzer then
   reported 4 `[wrong-value] in ibp._t_diff_product` with minimal counterexamples
   (`d/dt((2) s₂ s₂)`, rel_err 0.5). Total ≈ 70 s.

**Verdict on the tooling.** The per-target IBP sub-stages and the rule log did the
localisation with no code reading. The docstring formula next to the rule made the
mistake visible on sight. A hand check of one term was still needed to know which stage
was wrong, as in the review's three hunts. The fuzzer would have done it without that
check, but only after R1.

## Advisory (not blocking)

- **A1.** `rule()` logs a line whose two sides have the same `repr` when only `str`
  differs (`trace.py:117`), e.g. `pull.mvmul_scaled: -hat(Ω)·J*·Ω ==> -hat(Ω)·J*·Ω` in the
  rigid-pendulum log. To a reader the rule fired and changed nothing. Print `str` as well
  when the `repr`s agree, or log only when the `repr`s differ.
- **A2.** `repr` of nested negation is ambiguous: `--hat(v:\rho*)`,
  `hat(ρ)·-(-1)` (`simplify.transpose_scaled` line in the same log).
- **A3.** `compute_eom` hands context from `_separate_variations` back to itself through
  a private attribute set on the exception (`eom.py:555`, read at `:488`). It works, but a
  small exception subclass or a return value would be clearer.
- **A4.** `pipeline.md:241-246` and `troubleshooting.md:37,288` send users to a script under
  `specs/reviews/scripts/`. If the fuzzer is now part of the debugging story, move it to
  `scripts/` or into `geomech.numeric` (F7's `check_rewrite`).
- **A5.** `paper_context.md` is still gitignored (`.gitignore:7`) while the contributor notes cites
  it (unchanged from the implementation report; author decision).

## Phase 5 follow-ups surfaced

1. `vector_rules` and `pull` still do not descend into `Dot`, `MVMul`, `Hat`, `Transpose` or
   calculus nodes (`geomech/core/transformations/vector_identities.py:3`, the "recurses
   through" docstring; `pulling.py:45`). This is plan item 5.1; troubleshooting row
   `troubleshooting.md:17`.
2. `Vee` is not distributed over `MAdd` (5.2); the dead `vee_over_vadd` rule was deleted.
3. Rows in `Dot`/`Cross`/`Hat` and mixed-orientation `VAdd` are rejected by `check_tree`
   only, not at construction (`geomech/utils/invariants.py`; `conventions.md:51-53`).
   This needs a cheaper `is_row_vector` and a decision on `tests/test_geometry.py`
   `test_add_vector`/`test_sub_vector` and `tests/test_latex_precedence.py:50`.
4. `linearize` accepts `params` only to improve the error message; supporting prescribed
   parameters is still open (`geomech/dynamics/vbl.py:213-217`).
5. Constant-symbol variations print `0`/`O` instead of `0v`/`0` (`conventions.md:60-64`),
   which plan item 5.5's `zero_like` should cover.
6. Split `_extraction_error` (`geomech/dynamics/eom.py`, ≈ 90 lines with three message
   branches) and the `_eliminate` match (`simplification.py`) by node class (5.6). A3 fits
   here.
7. A rule-level numeric spot check: the fuzzer's jet evaluator exposed as
   `check_rewrite(before, after)` and optionally run per `rule()` firing under
   `GEOMECH_DEBUG` (review F7). Hunt (iv) shows it would replace the remaining hand check.
8. A1/A2 (rule-log readability).
9. Rule names still cover families: `collect.flip` has 6 call sites
   (`geomech/core/math/collecting.py`), and `ibp.dot_prefactor` hides
   `_t_diff_product`, where the hunt bug was. A `stacklevel` in the log line, or a
   sub-name for the product rule, would finish the localisation.
