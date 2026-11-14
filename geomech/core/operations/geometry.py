from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from geomech.core.base.expressions import (
    Expr,
    MatrixExpr,
    ScalarExpr,
    VectorExpr,
    _node_operands,
    _set_nodes,
)
from geomech.core.base.types import ExprType
from geomech.core.operations.addition import VAdd
from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
from geomech.core.operations.mixins import _BaseMixin, _BinaryMixin, _UnaryMixin
from geomech.core.operations.multiplication import MVMul, SVMul
from geomech.utils.errors import ExpressionMismatchError, SizeMismatchError

# ---------------------------------------------------------------------------
# Dot  (scalar result from two vectors)
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False, frozen=True)
class Dot(_BinaryMixin, ScalarExpr):
    """Dot product of two vectors → scalar."""

    nodes: tuple = ()

    def __init__(self, l, r):
        l, r = _node_operands(type(self).__name__, l, r)
        if l.type == ExprType.VECTOR and r.type == ExprType.VECTOR:
            ls = getattr(l, "size", None)
            rs = getattr(r, "size", None)
            if ls is not None and rs is not None and ls != rs:
                raise SizeMismatchError("Dot", ls, rs)
            _set_nodes(self, (l, r))
        else:
            raise ExpressionMismatchError("Dot", l.type, r.type, l, r)

    def __str__(self):
        return "Dot(" + str(self.left) + "," + str(self.right) + ")"


# ---------------------------------------------------------------------------
# Cross  (vector result from two vectors)
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False, frozen=True)
class Cross(_BinaryMixin, VectorExpr):
    """Cross product of two 3-vectors → vector."""

    nodes: tuple = ()

    def __init__(self, l, r):
        l, r = _node_operands(type(self).__name__, l, r)
        if l.type == ExprType.VECTOR and r.type == ExprType.VECTOR:
            ls = getattr(l, "size", None)
            rs = getattr(r, "size", None)
            if ls is not None and rs is not None and ls != rs:
                raise SizeMismatchError("Cross", ls, rs)
            _set_nodes(self, (l, r))
        else:
            raise ExpressionMismatchError("Cross", l.type, r.type, l, r)

    def __str__(self):
        return "Cross(" + str(self.left) + "," + str(self.right) + ")"


# ---------------------------------------------------------------------------
# Hat  (vector → skew-symmetric matrix)
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False, frozen=True)
class Hat(_UnaryMixin, MatrixExpr):
    """Hat map: R^3 → so(3)."""

    nodes: tuple = ()

    def __init__(self, expr):
        (expr,) = _node_operands("Hat", expr)
        if expr.type == ExprType.VECTOR:
            _set_nodes(self, (expr,))
        else:
            raise ExpressionMismatchError.unary("Hat", ExprType.VECTOR, expr)

    @property
    def is_skew_symmetric(self):
        """Hat(v) is skew-symmetric for every v."""
        return True

    def __str__(self):
        return "Hat(" + str(self.expr) + ")"


# ---------------------------------------------------------------------------
# Vee  (skew-symmetric matrix → vector)
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False, frozen=True)
class Vee(_UnaryMixin, VectorExpr):
    """Vee map: so(3) → R^3."""

    nodes: tuple = ()

    def __init__(self, expr):
        (expr,) = _node_operands("Vee", expr)
        if expr.type == ExprType.MATRIX:
            _set_nodes(self, (expr,))
        else:
            raise ExpressionMismatchError.unary("Vee", ExprType.MATRIX, expr)

    def __str__(self):
        return "Vee(" + str(self.expr) + ")"


# ---------------------------------------------------------------------------
# Transpose  (type-preserving unary)
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False, frozen=True)
class Transpose(_BaseMixin, Expr):
    """Transpose operator.  Preserves the type of its inner expression.

    ``Transpose(Transpose(x))`` is ``x``: the constructor returns the inner operand
    instead of building a double transpose, so a transposed row vector is a column
    again and never reaches the engine as a ``Transpose`` node.
    """

    nodes: tuple = ()

    def __new__(cls, expr=None) -> Any:
        if isinstance(expr, Transpose):
            return expr.expr
        return super().__new__(cls)

    def __init__(self, expr=None):
        if expr is None:
            # a childless Transpose had type None and leaked into the engine
            raise TypeError("Transpose: Transpose() needs an operand; write Transpose(expr)")
        _node_operands("Transpose", expr)
        if isinstance(expr, Transpose):
            # __new__ returned expr.expr, which is already built
            return
        _set_nodes(self, (expr,))

    @property
    def expr(self) -> Any:
        # None only for a malformed node without an operand, which check_tree reports
        return self.nodes[0] if self.nodes else None

    @property
    def type(self) -> ExprType:
        return self.expr.type if self.expr else None  # type: ignore[return-value]

    @property
    def is_constant(self):
        return self.expr.is_constant if self.expr else False

    @property
    def is_zero(self):
        return self.expr.is_zero if self.expr else False

    @property
    def is_symmetric(self):
        return self.expr.is_symmetric if self.expr else False

    @property
    def is_skew_symmetric(self):
        """The transpose of a skew-symmetric matrix (−X) is skew-symmetric."""
        return self.expr.is_skew_symmetric if self.expr else False

    @property
    def size(self):
        return getattr(self.expr, "size", None)

    def __str__(self):
        return "(" + str(self.expr) + ")'"

    def _identity_extra(self) -> tuple:
        """The operand's class (and what it adds): ``Transpose(Vector x)`` (a row) and
        ``Transpose(Matrix x)`` both print ``(x)'`` but differ. ``Hat`` and ``Vee`` need no
        such guard: their operand has one fixed type (``VECTOR``/``MATRIX``)."""
        expr = self.expr
        if expr is None:
            return ()
        return (type(expr), *expr._identity_extra())

    def delta(self):
        return Transpose(self.expr.delta())

    def t_diff(self):
        return Transpose(self.expr.t_diff())


# ---------------------------------------------------------------------------
# Row vectors
# ---------------------------------------------------------------------------


def is_row_vector(expr):
    """True if the vector-typed *expr* is a row vector.

    The tree has one ``VECTOR`` type for columns and rows; orientation is structural.
    A row is ``Transpose(column)``, a row times a matrix (``MVMul(vᵀ, M)``), a scalar
    multiple of a row, a sum whose first addend is a row, or a calculus operator
    (d/dt, δ, ∫) of a row. Every other vector expression is a column.
    """
    if getattr(expr, "type", None) != ExprType.VECTOR:
        return False
    match expr:
        case Transpose():
            return not is_row_vector(expr.expr)
        case MVMul():
            return expr.left.type == ExprType.VECTOR
        case SVMul():
            return is_row_vector(expr.left)
        case VAdd():
            return is_row_vector(expr.nodes[0])
        case TimeDerivative() | Variation() | TimeIntegral():
            return is_row_vector(expr.expr)
        case _:
            return False


def column_of(row):
    """The column vector whose transpose is the row vector *row*, built structurally.

    ``(vᵀ)ᵀ = v``, ``(rᵀ M)ᵀ = Mᵀ r``, ``(r s)ᵀ = rᵀ s``, a sum of rows maps addend by
    addend, and d/dt, δ, ∫ commute with the transpose. A row of any other shape falls back
    to ``Transpose(row)``. For a ``Transpose`` node this is the same object as
    ``Transpose(row)``; for the other row shapes (``aᵀK``, ``s aᵀ``) it gives a column the
    engine can work with, where ``Transpose(MVMul(aᵀ, K))`` would be opaque.
    """
    match row:
        case Transpose():
            return row.expr
        case MVMul() if row.left.type == ExprType.VECTOR:
            return MVMul(Transpose(row.right), column_of(row.left))
        case SVMul() if is_row_vector(row.left):
            return SVMul(column_of(row.left), row.right)
        case VAdd() if is_row_vector(row):
            return VAdd(*[column_of(n) for n in row.nodes])
        case TimeDerivative() if is_row_vector(row.expr):
            return TimeDerivative(column_of(row.expr))
        case Variation() if is_row_vector(row.expr):
            return Variation(column_of(row.expr))
        case TimeIntegral() if is_row_vector(row.expr):
            return TimeIntegral(column_of(row.expr))
        case _:
            return Transpose(row)
