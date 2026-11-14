# The Pipeline, Stage by Stage

`compute_eom(L, dW, variables)` runs the stages below, in this order (this is the order
of the code in `geomech/dynamics/eom.py`; the README, Getting Started and the
`compute_eom` docstring state the same order):

δS = δL + δW → manifold rules → `full_simplify` → integration by parts (`expand`,
`full_simplify`, then per variation vector: `collect` → IBP → `full_simplify`) →
`expand` → extract the coefficient of each variation vector → `full_simplify`.

Every stage is a pure function on an expression tree, so each one can be rerun by hand
from the previous one. The names in brackets are the stage names of
[`EOMTrace`](#seeing-the-stages), and the rule names (`vector.dot_tangent_orthogonal`, …)
are the ones logged by the [rule trace](#which-rule-fired).

## Validation (before any stage)

`compute_eom` first checks its inputs and raises a typed error naming the entry:
`SystemVariables` entries (type, constant, duplicate names, tangent vectors and composite
expressions listed as variables), `params` entries, `lagrangian`/`inf_work` types, and
`inf_work` holding a time derivative of a variation vector (virtual work must be linear
in the variations themselves). See [Troubleshooting](troubleshooting.md) for the
messages. With `GEOMECH_DEBUG=1`, `check_tree` also runs on `L` and `dW`.

## 1. Variation `[dL]`, `[dW]`, `[dS]`

`L.delta()` applies the product rule through the tree. Manifold leaves are eager:
`q.delta()` for an `S2` point is `ξ × q`, `R.delta()` for an `SO3` point is `R Hat(η)`,
and the variation of a velocity is built as the time derivative of the variation
(`x.t_diff().delta()` is `d/dt(δx)`). Tangent vectors (`ω`, `Ω`) are leaves at this
point, so their variations stay as `δω`, `δΩ`. δS is `δL + δW` (`dW=None` or `0` means
δW = 0).

## 2. Manifold rules `[manifold_rules]`

`apply_manifold_rules(dS, variables)` replaces the variation of each listed tangent
vector by its kinematic expansion (rules `manifold.s2_variation_omega`,
`manifold.so3_variation_Omega`):

- S²: `δω → ξ̇ − 2 ω × ξ` (Lee et al. eq. 5.13 on the constraint set ξ ⊥ q; the
  derivation, and why the factor is 2, is in
  [Manifold Rules: Derivation](manifold-implementation-notes.md));
- SO(3): `δΩ → η̇ + Ω × η` (Lee et al. eq. 6.9).

A tangent vector is matched by its printed name, so only the points listed in
`SystemVariables` are expanded.

## 3. Simplify `[simplify]`

`full_simplify` runs `expand` once, then repeats one `simplify` pass until the printed
tree stops changing (at most `max_iter=10` passes; hitting the cap emits the
`RuntimeWarning` "no confirmed fixed point after N passes"). One pass is

1. `pull`: hoist scalar factors out of products (`pull.dot_scaled`,
   `pull.mvmul_scaled`, …);
2. `vector_rules`: geometric identities, top down — orthogonality
   (`vector.dot_cross_orthogonal`), unit norm (`vector.dot_unit_norm`), tangency on S²
   (`vector.dot_tangent_orthogonal`, `vector.s2_tangent_recovery`, …). It recurses only
   through `Add`, `Mul`, `VAdd`, `SVMul` and `Cross`: a rule does not fire below a `Dot`
   operand, an `MVMul`, a `Hat` or a calculus node;
3. `_eliminate`: bottom up, zeros, identities and numbers (`simplify.mul_fold`,
   `simplify.mvmul_hat_to_cross`, `simplify.orthogonal_transpose`,
   `simplify.transpose_symmetric`, …).

The rules and the hypothesis each needs are listed in the docstrings of
`geomech.core.transformations.pulling`, `.vector_identities` and `.simplification`
(the modules of the functions `pull`, `vector_rules` and `simplify`).

## 4. Integration by parts `[ibp.prepare]`, `[ibp.collect[t]]`, `[ibp.integrate[t]]`, `[ibp.simplify[t]]`, `[ibp]`

`integrate_by_parts(dS, targets)` takes the time derivative of each variation vector as
a target (`d/dt(δx)`, `d/dt(ξ_q)`, `d/dt(η_R)`, in the order scalars, vectors,
matrices). It first runs `expand` and `full_simplify` (`[ibp.prepare]`), then for each
target in turn:

1. `collect(expr, target)` `[ibp.collect[t]]` rewrites every term so that the target is
   the *left* operand of a `Dot`: commutativity (`collect.flip`), the cyclic triple
   product (`collect.triple_product_left`, `_right`, `_left_nested`, … : one name per
   rotation), peeling a matrix off the side holding the target,
   `Dot(M f, w) = Dot(f, Mᵀ w)` (`collect.peel_matrix`), stripping a common rotation
   `Dot(R a, R b) = Dot(a, b)` (`collect.strip_rotation`), and turning an inner product
   written with rows, `vᵀ A w`, into a `Dot` (`collect.inner_product_to_dot`). A target
   inside the *matrix* factor of a product (`Hat(Ω)`, `Ω cᵀ`) cannot be exposed and
   raises `AlgebraicError` naming the term;
2. one IBP pass `[ibp.integrate[t]]`: `a Dot(ẋ, f) → a Dot(−x, ḟ)` (`ibp.dot`), with the
   product rule when the scalar prefactor `a` is not constant (`ibp.dot_prefactor`), and
   `ẋ g → −x ġ` for a scalar target in a product (`ibp.scalar_factor`). Boundary terms
   are dropped (Hamilton's principle);
3. `full_simplify` `[ibp.simplify[t]]`.

A target that `collect` could not bring into `Dot(target, f)` survives this stage; the
extraction in stage 6 then fails, and `compute_eom` names the term (see
[Troubleshooting](troubleshooting.md)).

## 5. Expand `[expand]`

`expand` distributes every product over sums (`expand.dot_over_vadd`,
`expand.mvmul_over_sum`, …) so that each additive term of δS carries exactly one
variation vector as a direct factor.

## 6. Extraction `[raw[ξ]]`, `[equation[ξ]]`

For each variation vector, `extract_linear_coeff(dS, vec)` returns its coefficient by
truncating the tree at the variation vector (the coefficient is always a left factor;
`geomech/core/math/extract.py` has the dispatch grid). A term where the variation
vector appears twice raises `NonLinearError`; one where it appears only inside a
calculus operator (a surviving `d/dt(η)`) raises `AlgebraicError`. The coefficient is
then simplified with `full_simplify`. The result is `{str(vec): (vec, E)}` with `E = 0`
the equation of motion (sign convention in [Conventions](conventions.md)).

Finally, `compute_eom` checks that no variation symbol of an unlisted symbol is left in
δS (a forgotten degree of freedom, a symbol missing `Constant`) and raises `ValueError`
naming it.

## After `compute_eom`

- `to_standard_form(eom, variables, inputs)` expands each equation, partitions its terms
  into acceleration, input and remaining terms, and extracts `M`, `G` and `f`
  (`E = M a + f + G u`).
- `linearize(eom, variables, inputs)` takes the variation δE of each equation, expands
  it, partitions its terms by the perturbation symbol each carries and extracts the
  blocks of `M δẋ = A δx + B δu` ([Variation-Based Linearization](vbl.md)).

## Seeing the stages

Pass an `EOMTrace` (keyword only); it is filled as the pipeline runs, so after an
exception it holds every stage up to the failing one:

```python
from geomech import EOMTrace, compute_eom

tr = EOMTrace()
eom = compute_eom(L, dW, variables, trace=tr)
tr.summary()                 # one line per stage: name, node count, start of the repr
tr.show("simplify")          # print_tree(..., style="indent") of one stage
tr.after_ibp                 # the same stages as attributes
tr.ibp["\\frac{d}{dt}(\\eta_{R})"]["collect"]   # per target: collect / integrate / simplify
tr.raw_equations, tr.equations
```

Two names look alike: the *stage* `"ibp"` (`tr["ibp"]`, `tr.show("ibp")`) is δS after
integration by parts, the attribute `tr.after_ibp`; the *attribute* `tr.ibp` is the dict
of per-target steps. To print one step, give its stage name:
`tr.show("ibp.collect[\\frac{d}{dt}(\\eta_{R})]")`.

The pipeline clears the trace when it starts, so passing the same `EOMTrace` to a second
`compute_eom` (re-running a notebook cell) keeps only the second run. When the run
raises, the trace keeps the stages before the failure, `tr.failed_at` names the stage
that raised (`raw[\eta_{R}]`, `ibp.collect[...]`, `validation`, `stray_check` …) and
`tr.error` the exception class; `summary()` then starts with
`FAILED at <stage>: <exception> (the stages below ran before it)`.

For the rigid pendulum (`L = ½ Ω·JΩ − m g e₃·Rρ`, `dW = η·M`), `tr.summary(chars=70)`
prints

```
dL                                        30 nodes  ((<δ(v:\Omega_{R}), M:J*·v:\Omega_{R}> + <v:\Omega_{R}, M:J*·δ(v:\O...
dW                                         3 nodes  <v:\eta_{R}, v:M>
dS                                        33 nodes  ((<δ(v:\Omega_{R}), M:J*·v:\Omega_{R}> + <v:\Omega_{R}, M:J*·δ(v:\O...
manifold_rules                            41 nodes  ((<(d/dt(v:\eta_{R}) + cross(v:\Omega_{R}, v:\eta_{R})), M:J*·v:\Om...
simplify                                  51 nodes  (<d/dt(v:\eta_{R}), M:J*·v:\Omega_{R}>·(0.5) + <cross(v:\Omega_{R},...
ibp.prepare                               51 nodes  (<d/dt(v:\eta_{R}), M:J*·v:\Omega_{R}>·(0.5) + <cross(v:\Omega_{R},...
ibp.collect[\frac{d}{dt}(\eta_{R})]       52 nodes  (<d/dt(v:\eta_{R}), M:J*·v:\Omega_{R}>·(0.5) + <cross(v:\Omega_{R},...
ibp.integrate[\frac{d}{dt}(\eta_{R})]     56 nodes  (<-v:\eta_{R}, M:J*·d/dt(v:\Omega_{R})>·(0.5) + <cross(v:\Omega_{R}...
ibp.simplify[\frac{d}{dt}(\eta_{R})]      43 nodes  (-<v:\eta_{R}, M:J*·d/dt(v:\Omega_{R})> + <cross(v:\Omega_{R}, v:\e...
ibp                                       43 nodes  (-<v:\eta_{R}, M:J*·d/dt(v:\Omega_{R})> + <cross(v:\Omega_{R}, v:\e...
expand                                    43 nodes  (-<v:\eta_{R}, M:J*·d/dt(v:\Omega_{R})> + <cross(v:\Omega_{R}, v:\e...
raw[\eta_{R}]                             41 nodes  (-M:J*·d/dt(v:\Omega_{R}) + (hat(v:\Omega_{R}))'·M:J*·v:\Omega_{R}·...
equation[\eta_{R}]                        27 nodes  (-M:J*·d/dt(v:\Omega_{R}) -cross(v:\Omega_{R}, M:J*·v:\Omega_{R}) +...
```

Every stage can be reproduced from the previous one with exported functions
(`apply_manifold_rules`, `full_simplify`, `integrate_by_parts`, `collect`, `expand`,
`extract_linear_coeff`); `tests/test_trace.py` does exactly that.

`linearize(..., trace=LinearizationTrace())` keeps, per dynamic row, `dE` (the expanded
variation of the equation), `buckets` (its terms grouped by the perturbation symbol each
carries), `raw_blocks` (the coefficient straight from `extract_linear_coeff`, before
`full_simplify` and before the mass blocks are negated) and `blocks` (the final entries
of `M`, `A`, `B`).

Without a trace object, the same stages are logged on the `geomech` logger at `DEBUG`
(one line per stage with its node count and the first 200 characters of the `repr`):

```python
import logging
logging.basicConfig()
logging.getLogger("geomech").setLevel(logging.DEBUG)
```

or run with `GEOMECH_DEBUG=1`, which also runs `check_tree` on the inputs and outputs of
`compute_eom`, `to_standard_form` and `linearize`.

`logging.basicConfig(level=logging.DEBUG)` works as well. Configure logging before
`set_trace(True)` (or before importing with `GEOMECH_DEBUG=1`): `geomech` attaches its own
`stderr` handler only when neither its loggers nor the root logger have one, so each
line is printed once.

## Which rule fired

`geomech.set_trace(True)` (or `GEOMECH_TRACE=1`) logs every firing of a named rewrite
rule on the `geomech.rules` logger as `name: before ==> after`:

```
geomech.rules: manifold.so3_variation_Omega: δ(v:\Omega_{R}) ==> (d/dt(v:\eta_{R}) + cross(v:\Omega_{R}, v:\eta_{R}))
geomech.rules: vector.dot_tangent_orthogonal: <v:\omega_{q}, S2:q> ==> 0
```

For the rigid pendulum above, 38 lines are logged. The rules that fire, with their line
counts, are `manifold.so3_variation_Omega` (2), `expand.dot_over_vadd` (2),
`expand.mul_over_add`, `expand.mvmul_over_sum`, `pull.dot_scaled` (2),
`pull.mvmul_scaled` (2), `pull.svmul_nested` (2), `pull.smmul_nested`, `collect.flip`,
`collect.peel_matrix`, `ibp.dot` (2), and in `_eliminate` `simplify.mul_fold` (4),
`simplify.mvmul_hat_to_cross` (3), `simplify.mvmul_reassociate` (3),
`simplify.transpose_hat_skew` (3), `simplify.transpose_symmetric` (2),
`simplify.transpose_product` (2), `simplify.transpose_scaled`, `simplify.add_fold` (2),
`simplify.vadd_fold`. A `*_fold` line means the fold itself changed the node (zeros
dropped, numbers combined, like terms merged); a rewrite made by a child rule is logged
once, by that rule, not again by every enclosing sum or product. The same rewrite can
still be logged once per pass of `full_simplify` and once per copy of the subterm. A
rule that did not fire is often the answer to "why is this term not reduced": check the
hypothesis listed next to its name in the module docstring.

| Module | Rule prefix | Where it runs |
|---|---|---|
| `geomech.core.transformations.manifold_rules` | `manifold.` | stage 2 |
| `geomech.core.transformations.expansion` | `expand.` | stages 3–5 (inside `full_simplify` and on its own) |
| `geomech.core.transformations.pulling` | `pull.` | every `simplify` pass |
| `geomech.core.transformations.vector_identities` | `vector.` | every `simplify` pass |
| `geomech.core.transformations.simplification` | `simplify.` | every `simplify` pass (`_eliminate`) |
| `geomech.core.math.collecting` | `collect.` | stage 4, per target |
| `geomech.core.math.ibp` | `ibp.` | stage 4, per target |

`tests/test_rule_names.py` checks that the names in each docstring and the names in the
code are the same set, and `tests/test_rule_coverage.py` that every rewrite of random
trees logs at least one named rule (no rewrite bypasses `rule()`).

## When the equation is wrong

The trace shows *where* each change happens, not which change is false. To find a
false rewrite without deriving the equation by hand, run the engine fuzzer from the
repository root:

```
python specs/reviews/scripts/engine_fuzz.py --quick
```

It applies every rewrite to thousands of random trees and compares input and output
numerically; a `[wrong-value]` finding names the transformation and a minimal
counterexample. [Troubleshooting](troubleshooting.md#the-equation-is-wrong-and-no-row-matches)
has the full recipe: bisect with `summary()`, run the fuzzer, rerun the suspect stage
under `set_trace`.
