from __future__ import annotations

from dataclasses import dataclass

from geomech.core.base.expressions import (
    MANIFOLD_POINTS,
    S1,
    S2,
    TANGENT_LEAVES,
    Expr,
    _set_nodes,
    manifold_point_of,
)
from geomech.core.operations.mixins import _CalcUnaryMixin


def _reject_manifold_leaf(expr, op, method, s1_form, s2_form, so3_form):
    """``TimeDerivative(q)``/``Variation(q)`` for a manifold point is not a tree node.

    The kinematics of ``S1``/``S2``/``SO3`` leaves are eager: ``q.t_diff()`` returns
    ``ω × q`` (``(e × q) ω`` on S¹) and ``R.delta()`` returns ``R·Hat(η)``. A calculus node
    wrapping the leaf itself would carry none of that structure, so it is rejected at
    construction.
    """
    if isinstance(expr, MANIFOLD_POINTS):
        if isinstance(expr, S1):
            form = s1_form
        else:
            form = s2_form if isinstance(expr, S2) else so3_form
        raise ValueError(
            f"{op}: {op}({expr}) is not represented in the tree for the manifold point"
            f" '{expr}'; use {expr}.{method}(), which returns the kinematic expression {form}"
        )


# ---------------------------------------------------------------------------
# Variation  (δ operator as a tree node)
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False, frozen=True)
class Variation(_CalcUnaryMixin, Expr):
    """Variation operator δ{expr}.  Preserves the type of its inner expression."""

    nodes: tuple = ()

    def __init__(self, expr):
        _reject_manifold_leaf(expr, "Variation", "delta", "(e × q) γ", "ξ × q", "R·Hat(η)")
        # δ(δ(...)) is a second-order variation — vanishes in Hamilton's principle
        if isinstance(expr, Variation):
            raise ValueError(
                f"Cannot take variation of a variation: δ(δ({expr.expr})). "
                f"Second-order variations vanish in Hamilton's principle."
            )
        # δ(γ), δ(ξ) or δ(η) where γ/ξ/η are already variation vectors is nonsensical
        point = manifold_point_of(expr) if isinstance(expr, TANGENT_LEAVES) else None
        if point is not None:
            var_name = point.manifold_info.variation_vector_name
            if expr.name == var_name:
                raise ValueError(
                    f"Cannot take variation of variation vector {expr}. "
                    f"The variation vector is already an infinitesimal perturbation."
                )
        _set_nodes(self, (expr,))

    def __str__(self):
        return "\\delta{" + str(self.expr) + "}"

    def t_diff(self):
        """d/dt(δx) — put TimeDerivative on the outside.

        δ and d/dt commute, but the expression tree produces
        TimeDerivative(Variation(x)) when taking the variation of d/dt(x).
        This override ensures the IBP target matches that structure.
        """
        return TimeDerivative(self)


# ---------------------------------------------------------------------------
# TimeDerivative  (d/dt operator as a tree node)
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False, frozen=True)
class TimeDerivative(_CalcUnaryMixin, Expr):
    """Time derivative d/dt{expr}.  Preserves the type of its inner expression."""

    nodes: tuple = ()

    def __init__(self, expr):
        _reject_manifold_leaf(expr, "TimeDerivative", "t_diff", "(e × q) ω", "ω × q", "R·Hat(Ω)")
        _set_nodes(self, (expr,))

    def __str__(self):
        return "\\frac{d}{dt}(" + str(self.expr) + ")"

    def t_integrate(self):
        """∫(d/dt(x)) dt = x — cancellation."""
        return self.expr


# ---------------------------------------------------------------------------
# TimeIntegral  (∫ dt operator as a tree node)
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False, frozen=True)
class TimeIntegral(_CalcUnaryMixin, Expr):
    """Time integral ∫{expr}dt.  Preserves the type of its inner expression."""

    nodes: tuple = ()

    def __init__(self, expr):
        _set_nodes(self, (expr,))

    def __str__(self):
        return "\\int{" + str(self.expr) + "}dt"

    def t_diff(self):
        """d/dt(∫x dt) = x — cancellation."""
        return self.expr
