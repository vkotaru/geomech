from __future__ import annotations

from dataclasses import dataclass, field

from geomech.core.base.expressions import (
    Expr,
    MatrixExpr,
    ScalarExpr,
    VectorExpr,
    _node_operands,
)
from geomech.core.base.types import ExprType
from geomech.core.operations.mixins import _BinaryMixin
from geomech.utils.errors import ExpressionMismatchError, SizeMismatchError

# ---------------------------------------------------------------------------
# Binary multiplication helpers
# ---------------------------------------------------------------------------


def _wrap_numeric(val):
    """Wrap a number (int, float, numpy scalar) as a ``Number``; pass anything else through."""
    from geomech.core.base.expressions import Number, _is_number

    if _is_number(val):
        return Number(val)
    return val


def _is_row(expr):
    from geomech.core.operations.geometry import is_row_vector

    return is_row_vector(expr)


def _q(expr):
    from geomech.utils.errors import _operand_str

    return _operand_str(expr)


def _check_mul_sizes(op_name: str, l, r):
    """Validate matrix/vector dimension compatibility for multiplication.

    Rules (when both operands carry a known size):
      MVMul: M(m,n) * v(p,) → n must equal p
      MMMul: M(m,n) * N(p,q) → n must equal p
      VVMul (inner): v(m,)^T * w(p,) → m must equal p
    """
    ls = getattr(l, "size", None)
    rs = getattr(r, "size", None)
    if ls is None or rs is None:
        return
    # columns of left must match rows/length of right
    l_cols = ls[-1]  # last dim: cols for matrix, length for vector
    r_rows = rs[0]  # first dim: rows for matrix, length for vector
    if l_cols != r_rows:
        raise SizeMismatchError(op_name, "left cols=" + str(l_cols), "right rows=" + str(r_rows))


# ---------------------------------------------------------------------------
# Scalar * Scalar
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False)
class Mul(_BinaryMixin, ScalarExpr):
    """Scalar multiplication."""

    nodes: list = field(default_factory=list)

    def __init__(self, l, r):
        l, r = _node_operands(type(self).__name__, l, r, numeric=True)
        if l.type == ExprType.SCALAR and r.type == ExprType.SCALAR:
            self.nodes = [l, r]
        else:
            raise ExpressionMismatchError("Mul", l.type, r.type, l, r)

    def __str__(self):
        return str(self.left) + str(self.right)


# ---------------------------------------------------------------------------
# Scalar * Vector
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False)
class SVMul(_BinaryMixin, VectorExpr):
    """Scalar-Vector multiplication.  Normalized: left=vector, right=scalar."""

    nodes: list = field(default_factory=list)

    def __init__(self, l, r):
        l, r = _node_operands(type(self).__name__, l, r, numeric=True)
        if l.type == ExprType.VECTOR and r.type == ExprType.SCALAR:
            self.nodes = [l, r]
        elif l.type == ExprType.SCALAR and r.type == ExprType.VECTOR:
            self.nodes = [r, l]
        else:
            raise ExpressionMismatchError("SVMul", l.type, r.type, l, r)

    def __str__(self):
        return str(self.left) + str(self.right)


# ---------------------------------------------------------------------------
# Scalar * Matrix
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False)
class SMMul(_BinaryMixin, MatrixExpr):
    """Scalar-Matrix multiplication.  Normalized: left=matrix, right=scalar."""

    nodes: list = field(default_factory=list)

    def __init__(self, l, r):
        l, r = _node_operands(type(self).__name__, l, r, numeric=True)
        if l.type == ExprType.MATRIX and r.type == ExprType.SCALAR:
            self.nodes = [l, r]
        elif l.type == ExprType.SCALAR and r.type == ExprType.MATRIX:
            self.nodes = [r, l]
        else:
            raise ExpressionMismatchError("SMMul", l.type, r.type, l, r)

    @property
    def is_symmetric(self):
        return self.left.is_symmetric

    @property
    def is_skew_symmetric(self):
        """A scalar multiple of a skew-symmetric matrix is skew-symmetric."""
        return self.left.is_skew_symmetric

    def __str__(self):
        return str(self.left) + str(self.right)


# ---------------------------------------------------------------------------
# Matrix * Vector
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False)
class MVMul(_BinaryMixin, VectorExpr):
    """Matrix-Vector multiplication."""

    nodes: list = field(default_factory=list)

    def __init__(self, l, r):
        l, r = _node_operands("MVMul", l, r)
        if l.type == ExprType.MATRIX and r.type == ExprType.VECTOR:
            if _is_row(r):
                raise ExpressionMismatchError(
                    msg=(
                        f"MVMul: a matrix times a row vector is not defined: {_q(l)} * {_q(r)};"
                        " a matrix multiplies a column vector (M * v), and a row vector"
                        " multiplies a matrix from the left (v.T() * M)"
                    )
                )
            _check_mul_sizes("MVMul", l, r)
            self.nodes = [l, r]
        elif l.type == ExprType.VECTOR and r.type == ExprType.MATRIX:
            # row vector times matrix: vᵀ M (a row)
            if _is_row(l):
                self.nodes = [l, r]
            else:
                raise ExpressionMismatchError(
                    msg=(
                        f"MVMul: a column vector times a matrix is not defined: {_q(l)} *"
                        f" {_q(r)}; write M * v for the column M v, or v.T() * M for the"
                        " row vᵀ M"
                    )
                )
        else:
            raise ExpressionMismatchError("MVMul", l.type, r.type, l, r)

    def __str__(self):
        return str(self.left) + str(self.right)


# ---------------------------------------------------------------------------
# Matrix * Matrix
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False)
class MMMul(_BinaryMixin, MatrixExpr):
    """Matrix-Matrix multiplication."""

    nodes: list = field(default_factory=list)

    def __init__(self, l, r):
        l, r = _node_operands("MMMul", l, r)
        if l.type == ExprType.MATRIX and r.type == ExprType.MATRIX:
            _check_mul_sizes("MMMul", l, r)
            self.nodes = [l, r]
        else:
            raise ExpressionMismatchError("MMMul", l.type, r.type, l, r)

    def __str__(self):
        return str(self.left) + str(self.right)


# ---------------------------------------------------------------------------
# Vector * Vector
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False)
class VVMul(_BinaryMixin, Expr):
    """Vector-Vector multiplication.
    Result type depends on Transpose:
      Transpose(v) * w  → scalar
      v * Transpose(w)  → matrix
    """

    nodes: list = field(default_factory=list)

    def __init__(self, l, r):
        l, r = _node_operands("VVMul", l, r)
        if l.type == ExprType.VECTOR and r.type == ExprType.VECTOR:
            l_row, r_row = _is_row(l), _is_row(r)
            if l_row and not r_row:
                # v^T * w → scalar: sizes must match
                _check_mul_sizes("VVMul", l, r)
                self.nodes = [l, r]
                self._result_type = ExprType.SCALAR
            elif not l_row and r_row:
                # v * w^T → matrix (outer product): any sizes valid
                self.nodes = [l, r]
                self._result_type = ExprType.MATRIX
            elif l_row:
                raise ExpressionMismatchError(
                    msg=(
                        f"VVMul: two row vectors cannot be multiplied: {_q(l)} * {_q(r)};"
                        " write a.T() * b for the inner product aᵀb or a * b.T() for the"
                        " outer product a bᵀ"
                    )
                )
            else:
                raise ExpressionMismatchError(
                    msg=(
                        f"VVMul: two column vectors cannot be multiplied: {_q(l)} * {_q(r)};"
                        " write Dot(a, b) or a.T() * b for the inner product, a * b.T() for"
                        " the outer product, or Cross(a, b)"
                    )
                )
        else:
            raise ExpressionMismatchError("VVMul", l.type, r.type, l, r)

    @property
    def type(self):
        return self._result_type

    def __str__(self):
        return str(self.left) + str(self.right)
