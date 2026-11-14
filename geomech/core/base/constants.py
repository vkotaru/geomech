"""Shared numeric literals and the typed zero / negation helpers.

Every rewrite that needs ``-1``, ``-2`` or ``1/2`` uses the leaves below instead of
building its own: ``NEG_ONE`` prints ``(-1)``, exactly as ``Number(-1)`` and the ``-1`` the
operators wrap (``x * -1``, ``-x``), so a rule's output and a user's ``-x`` are the same
tree. The leaves are ordinary ``Number`` literals (``is_number_literal`` is true for each),
shared by every tree that uses them; that is safe because nodes are never mutated.

``negate(expr)`` is ``expr * (-1)`` for an expression of any type (``Mul``, ``SVMul`` or
``SMMul``; a row vector stays a row), and ``zero_like(expr)`` is the zero leaf of the type
of *expr* (``Zero``, ``ZeroVector``, ``ZeroMatrix``); ``zero_of(expr_type)`` the same
for a bare ``ExprType``.
"""

from __future__ import annotations

from geomech.core.base.expressions import Expr, Number, Zero, ZeroMatrix, ZeroVector
from geomech.core.base.types import ExprType

#: The literal -1, printed ``(-1)``.
NEG_ONE = Number(-1)
#: The literal -2, printed ``(-2)`` (the S² variation rule ``ξ̇ − 2 ω × ξ``).
NEG_TWO = Number(-2)
#: The literal 1/2, printed ``(0.5)``.
HALF = Number(0.5)


def negate(expr: Expr) -> Expr:
    """``expr · (-1)``: ``Mul`` for a scalar, ``SVMul`` for a vector, ``SMMul`` for a matrix.

    The same tree as ``-expr`` (``Expr.__neg__``): a scaled row vector stays a row
    (``(v s)ᵀ``). Raises ``TypeError`` for an object without an expression type.
    """
    if getattr(expr, "type", None) not in (ExprType.SCALAR, ExprType.VECTOR, ExprType.MATRIX):
        raise TypeError(f"negate: cannot negate {type(expr).__name__} {expr!r}")
    return expr * NEG_ONE


def zero_like(expr: Expr) -> Expr:
    """The zero leaf of the type of *expr*: ``ZeroVector``, ``ZeroMatrix``, else ``Zero``.

    An object without an expression type (``None``) gets the scalar ``Zero``.
    """
    return zero_of(getattr(expr, "type", None))


def zero_of(expr_type: ExprType | None) -> Expr:
    """The zero leaf of *expr_type* (``Zero`` for ``SCALAR`` and for ``None``)."""
    match expr_type:
        case ExprType.VECTOR:
            return ZeroVector
        case ExprType.MATRIX:
            return ZeroMatrix
        case _:
            return Zero


__all__ = ["HALF", "NEG_ONE", "NEG_TWO", "negate", "zero_like", "zero_of"]
