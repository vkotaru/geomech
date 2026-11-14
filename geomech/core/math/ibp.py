"""Integration by parts for scalar expressions.

    integrate_by_parts(expr, targets) — expand → full_simplify, then for each target in
    order: collect → one IBP pass → full_simplify

IBP identity (boundary terms dropped):
    ∫ Dot(ẋ, f) dt  =  -∫ Dot(x, ḟ) dt

So:  Dot(target, rhs)  →  Dot(-target.t_integrate(), rhs.t_diff())

Only a ``Dot`` whose *left* operand is the target, or a scalar ``Mul`` chain with the target
as a direct factor, is integrated; ``collect`` is what brings each term into that shape.
A target left anywhere else survives, and extraction then fails on it (``compute_eom``
names the term).

Rules (name, identity, hypothesis):

  ibp.dot             a Dot(ẋ, f) → a Dot(-x, ḟ)       target is the left operand; a is
                                                       constant (or absent)
  ibp.dot_prefactor   a Dot(ẋ, f) → a Dot(-x, ḟ) - Dot(x, ȧ f)
                                                       a is a non-constant scalar factor
  ibp.scalar_factor   ẋ g → -x ġ                       target is a factor of a Mul chain
  ibp.product_rule    d/dt(∏ fᵢ^nᵢ) → Σ nᵢ fᵢ^(nᵢ-1) ḟᵢ ∏_{j≠i} fⱼ^nⱼ
                                                       the rate ȧ of ibp.dot_prefactor;
                                                       constant factors are skipped and
                                                       repeated factors merged
"""

from __future__ import annotations

from geomech.core.base.expressions import One, Zero
from geomech.core.math.collecting import collect
from geomech.core.operations.addition import Add
from geomech.core.operations.calculus import TimeDerivative
from geomech.core.operations.geometry import Dot
from geomech.core.operations.multiplication import Mul, SVMul
from geomech.core.transformations.expansion import expand
from geomech.core.transformations.simplification import full_simplify
from geomech.utils.trace import rule, stage_recorder

# ---------------------------------------------------------------------------
# Single IBP pass
# ---------------------------------------------------------------------------


def _remove_factor(expr, target):
    """Remove *target* from a Mul chain, returning the remaining product.

    If target is not a direct factor, returns None.
    Example: _remove_factor(Mul(Mul(a, target), b), target) → Mul(a, b)
    """
    if expr == target:
        return One
    if isinstance(expr, Mul):
        # Try removing from left
        left_removed = _remove_factor(expr.left, target)
        if left_removed is not None:
            if str(left_removed) == "1":
                return expr.right
            return Mul(left_removed, expr.right)
        # Try removing from right
        right_removed = _remove_factor(expr.right, target)
        if right_removed is not None:
            if str(right_removed) == "1":
                return expr.left
            return Mul(expr.left, right_removed)
    return None


def _split_dot_factor(expr, target):
    """Split a Mul chain into ``(prefactor, Dot(target, rhs))``; None if no such Dot factor.

    Example: Mul(Mul(a, Dot(target, r)), b) → (Mul(a, b), Dot(target, r))
    """
    if isinstance(expr, Dot):
        return (One, expr) if expr.left == target else None
    if isinstance(expr, Mul):
        found = _split_dot_factor(expr.left, target)
        if found is not None:
            pre, dot = found
            return (expr.right if pre is One else Mul(pre, expr.right)), dot
        found = _split_dot_factor(expr.right, target)
        if found is not None:
            pre, dot = found
            return (expr.left if pre is One else Mul(expr.left, pre)), dot
    return None


def _replace_in_chain(expr, old, new):
    """Return the Mul chain *expr* with the factor *old* (by identity) replaced by *new*."""
    if expr is old:
        return new
    if isinstance(expr, Mul):
        return Mul(_replace_in_chain(expr.left, old, new), _replace_in_chain(expr.right, old, new))
    return expr


def _flatten_mul(expr, factors):
    """Flatten a nested Mul chain into its list of factors (left to right)."""
    if isinstance(expr, Mul):
        _flatten_mul(expr.left, factors)
        _flatten_mul(expr.right, factors)
    else:
        factors.append(expr)


def _t_diff_product(expr):
    """d/dt of a scalar product with repeated factors merged:
    d/dt(∏ᵢ fᵢ^{nᵢ}) = Σᵢ nᵢ fᵢ^{nᵢ-1} ḟᵢ ∏_{j≠i} fⱼ^{nⱼ}  (e.g. d/dt(m l l) = 2 m l l̇).
    """
    factors = []
    _flatten_mul(expr, factors)
    terms = []
    seen = []
    for f in factors:
        if f.is_constant or any(f == g for g in seen):
            continue
        seen.append(f)
        count = sum(1 for g in factors if g == f)
        last = max(j for j, g in enumerate(factors) if g == f)
        new = [f.t_diff() if j == last else g for j, g in enumerate(factors)]
        term = new[0]
        for g in new[1:]:
            term = Mul(term, g)
        if count > 1:
            term = Mul(term, count)
        terms.append(term)
    if not terms:
        return Zero
    if len(terms) == 1:
        return terms[0]
    return Add(*terms)


def _apply_ibp(expr, target):
    """Apply integration by parts to *expr* with respect to *target*.

    Recurses through Add and Mul.  At each Dot where the left side
    equals *target*, applies the IBP identity:
      ∫ Dot(ẋ, f) dt = -∫ Dot(x, ḟ) dt

    For scalar Mul chains where *target* is a direct factor:
      ∫ ẋ * f dt = -∫ x * ḟ dt

    For a Mul chain with a non-constant scalar prefactor *a* on a Dot:
      ∫ a Dot(ẋ, f) dt = -∫ Dot(x, d/dt(a f)) dt = -∫ Dot(x, ȧ f + a ḟ) dt
    """
    match expr:
        case Add(nodes=nodes):
            return Add(*[_apply_ibp(n, target) for n in nodes])

        case Dot():
            if expr.left == target:
                return rule(
                    "ibp.dot", expr, Dot(expr.left.t_integrate() * (-1), expr.right.t_diff())
                )
            return expr

        case Mul():
            # Scalar IBP: find target as a direct factor in the Mul chain
            remainder = _remove_factor(expr, target)
            if remainder is not None:
                return rule(
                    "ibp.scalar_factor",
                    expr,
                    Mul(target.t_integrate() * (-1), remainder.t_diff()),
                )
            found = _split_dot_factor(expr, target)
            if found is None:
                # No Dot(ẋ, ·) factor in this chain: recurse to reach deeper nodes
                return Mul(_apply_ibp(expr.left, target), _apply_ibp(expr.right, target))
            # a Dot(ẋ, f) → a Dot(-x, ḟ), replacing the Dot in place (tree shape preserved)
            prefactor, dot = found
            out = _replace_in_chain(expr, dot, _apply_ibp(dot, target))
            if prefactor.is_constant:
                return out
            # Non-constant prefactor a: the product rule also needs -Dot(x, ȧ f)
            a_rate = rule(
                "ibp.product_rule", TimeDerivative(prefactor), _t_diff_product(prefactor)
            )
            a_dot_f = SVMul(dot.right, a_rate)
            return rule(
                "ibp.dot_prefactor", expr, Add(out, Dot(target.t_integrate() * (-1), a_dot_f))
            )

        case _:
            return expr


# ---------------------------------------------------------------------------
# Full IBP pipeline
# ---------------------------------------------------------------------------


def integrate_by_parts(expr, targets):
    """Expand, simplify, then apply collect + IBP for each target.

    *targets* is a list of variation vector time-derivatives (e.g. dot_eta)
    that should be integrated by parts to remove the time derivative.

    With the ``geomech`` logger at ``DEBUG`` each stage (``ibp.prepare``, then per target
    ``ibp.collect[t]``, ``ibp.integrate[t]``, ``ibp.simplify[t]``) is logged as
    ``integrate_by_parts: <stage> [<n> nodes] <repr>``; otherwise nothing is recorded.
    """
    return _integrate_by_parts(expr, targets, stage_recorder(None, "integrate_by_parts"))


def _integrate_by_parts(expr, targets, record=None):
    """``integrate_by_parts`` with an optional ``record(stage, expr)`` hook.

    The stages are ``ibp.prepare`` (after expand + full_simplify) and, per target,
    ``ibp.collect[t]``, ``ibp.integrate[t]`` and ``ibp.simplify[t]``. ``compute_eom``
    passes its stage recorder (trace object and/or the ``geomech`` debug logger).
    """
    expr = expand(expr)
    expr = full_simplify(expr)
    if record is not None:
        record("ibp.prepare", expr)
    for target in targets:
        expr = collect(expr, target)
        if record is not None:
            record(f"ibp.collect[{target}]", expr)
        expr = _apply_ibp(expr, target)
        if record is not None:
            record(f"ibp.integrate[{target}]", expr)
        expr = full_simplify(expr)
        if record is not None:
            record(f"ibp.simplify[{target}]", expr)
    return expr
