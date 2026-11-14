from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from geomech.core.base.types import ExprFlags, ExprType, ManifoldInfo
from geomech.utils.errors import (
    AlgebraicError,
    ExpressionMismatchError,
    _operand_str,
    _short_repr,
)

# ---------------------------------------------------------------------------
# Base
# ---------------------------------------------------------------------------


def _attr_list(attr, name):
    """Copy *attr* into a fresh list (never alias the caller's), rejecting a bare string.

    ``attr="Constant"`` would otherwise be iterated character by character and fail with
    ``unknown attribute 'C'``, which hides the actual mistake.
    """
    if attr is None:
        return None
    if isinstance(attr, str):
        where = f" for symbol '{name}'" if name is not None else ""
        raise TypeError(
            f"attr must be a list of strings, got str {attr!r}{where}; write attr=[{attr!r}]"
        )
    return list(attr)


def _is_number(val):
    """True for a Python or numpy int/float; ``bool`` (and ``numpy.bool_``) is not a number."""
    if isinstance(val, (bool, np.bool_)):
        return False
    return isinstance(val, (int, float, np.integer, np.floating))


def _describe(expr):
    """``Scalar 'm'``: class name and quoted ``str`` of an operand, for error messages."""
    return f"{type(expr).__name__} {_operand_str(expr)}"


def _as_operand(val, op, other):
    """Return the operand *val* of ``other <op> val`` as an expression.

    Expressions pass through and numbers (``int``, ``float`` and the numpy scalar types)
    become ``Number(val)``. Anything else, ``bool`` included, raises ``TypeError`` naming
    the expression on the other side.
    """
    if isinstance(val, Expr):
        return val
    if _is_number(val):
        return Number(val)
    if isinstance(val, (bool, np.bool_)):
        raise TypeError(
            f"cannot apply '{op}' to {_describe(other)}: got bool {val}; a bool is not a"
            " number here, write 1 or 0"
        )
    raise TypeError(
        f"cannot apply '{op}' to {_describe(other)} and {type(val).__name__} {_short_repr(val)};"
        " operands must be geomech expressions or numbers (int, float)"
    )


def _add(a, b):
    """``a + b`` for two expressions, dispatched on the type of *a*."""
    from geomech.core.operations import Add, MAdd, VAdd

    match a.type:
        case ExprType.SCALAR:
            return Add(a, b)
        case ExprType.VECTOR:
            return VAdd(a, b)
        case _:
            return MAdd(a, b)


def _sub(a, b):
    """``a - b`` = ``a + b·(-1)`` for two expressions of the same type."""
    from geomech.core.operations import Add, MAdd, Mul, SMMul, SVMul, VAdd

    if a.type != b.type:
        raise ExpressionMismatchError("Sub", a.type, b.type, a, b)
    match a.type:
        case ExprType.SCALAR:
            return Add(a, Mul(b, -1))
        case ExprType.VECTOR:
            return VAdd(a, SVMul(b, -1))
        case _:
            return MAdd(a, SMMul(b, -1))


def _mul(a, b):
    """``a * b`` for two expressions: the node class is picked from the operand types.

    A scalar factor is normalised to the right of ``SVMul``/``SMMul``. A scalar on the
    left of a row ``vᵀ`` keeps the pre-phase-3 form ``SVMul(vᵀ, s)`` (printed ``(v)'s``),
    while ``vᵀ * s`` is ``(v s)ᵀ``; both are one transpose of a column. Vector times
    vector is ``VVMul`` (inner product for a row times a column, outer product for a
    column times a row); a row vector times a matrix is the row ``MVMul(vᵀ, M)``. A
    transposed matrix dispatches like any matrix (``Rᵀ v`` is ``MVMul``).
    """
    from geomech.core.operations import MMMul, Mul, MVMul, SMMul, SVMul, Transpose, VVMul
    from geomech.core.operations.geometry import is_row_vector

    def _scale(v, s):
        if isinstance(v, Transpose) and not is_row_vector(v.expr):
            return Transpose(SVMul(v.expr, s))
        return SVMul(v, s)

    match (a.type, b.type):
        case (ExprType.SCALAR, ExprType.SCALAR):
            return Mul(a, b)
        case (ExprType.SCALAR, ExprType.VECTOR):
            return SVMul(b, a)
        case (ExprType.SCALAR, ExprType.MATRIX):
            return SMMul(b, a)
        case (ExprType.VECTOR, ExprType.SCALAR):
            return _scale(a, b)
        case (ExprType.VECTOR, ExprType.VECTOR):
            return VVMul(a, b)
        case (ExprType.VECTOR, ExprType.MATRIX):
            return MVMul(a, b)
        case (ExprType.MATRIX, ExprType.SCALAR):
            return SMMul(a, b)
        case (ExprType.MATRIX, ExprType.VECTOR):
            return MVMul(a, b)
        case _:
            return MMMul(a, b)


def _node_operands(op, *vals, numeric=False):
    """Validate the operands of the node constructor *op*.

    Expressions pass through; with ``numeric=True`` numbers are wrapped in ``Number``
    (the scalar slot of ``Mul``/``SVMul``/``SMMul``). Anything else raises ``TypeError``
    naming the node, instead of an ``AttributeError`` from inside the constructor.
    """
    out = []
    for v in vals:
        if isinstance(v, Expr):
            out.append(v)
        elif numeric and _is_number(v):
            out.append(Number(v))
        else:
            kinds = "expressions or numbers" if numeric else "expressions"
            raise TypeError(
                f"{op}: operands must be geomech {kinds}, got {type(v).__name__} {_short_repr(v)}"
            )
    return out


@dataclass(eq=False, repr=False)
class Expr:
    """Base expression node.

    Python operators build the tree: ``+``, ``-``, ``*`` and unary ``-`` accept another
    expression or a number on either side (``2 * x``, ``x + 1``, ``1 - x``, ``-x``);
    ``x / 2`` divides by a number only.
    """

    # numpy scalars (np.float64(2) * x) defer to __rmul__ & co. instead of broadcasting
    __array_ufunc__ = None

    # ------ arithmetic (one dispatch for leaves and operation nodes) ------

    def __add__(self, other):
        return _add(self, _as_operand(other, "+", self))

    def __radd__(self, other):
        if _is_number(other) and other == 0:
            # 0 + x is x, for every type: sum([a, b, c]) starts from int 0
            return self
        return _add(_as_operand(other, "+", self), self)

    def __sub__(self, other):
        if _is_number(other):
            if self.type != ExprType.SCALAR:
                # name the operands as written: 'a - 1', not the internal 'a + (-1)'
                raise ExpressionMismatchError(
                    "Sub", self.type, ExprType.SCALAR, self, Number(other)
                )
            # x - 1 is x + (-1), not x + (1)(-1)
            return _add(self, Number(-other))
        return _sub(self, _as_operand(other, "-", self))

    def __rsub__(self, other):
        return _sub(_as_operand(other, "-", self), self)

    def __mul__(self, other):
        return _mul(self, _as_operand(other, "*", self))

    def __rmul__(self, other):
        return _mul(_as_operand(other, "*", self), self)

    def __neg__(self):
        # the same tree as ``self * -1``, so str(-x) == str(x * -1)
        return self * -1

    def __truediv__(self, other):
        if isinstance(other, Expr):
            raise TypeError(
                f"cannot divide {_describe(self)} by the expression {_describe(other)};"
                " only division by a number is supported. For a symbolic divisor declare"
                " its inverse as a symbol, e.g. Scalar('1/m', attr=['Constant'])"
            )
        if not _is_number(other):
            _as_operand(other, "/", self)  # raises TypeError naming the operand
        if other == 0:
            raise ZeroDivisionError(f"cannot divide {_describe(self)} by zero")
        return self * (1 / other)

    def __rtruediv__(self, other):
        raise TypeError(
            f"cannot divide {other!r} by the expression {_describe(self)}; only division"
            " by a number is supported. For a symbolic divisor declare its inverse as a"
            " symbol, e.g. Scalar('1/m', attr=['Constant'])"
        )

    def __eq__(self, other):
        if type(self) is not type(other):
            return False
        return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        raise NotImplementedError

    def __repr__(self):
        from geomech.utils.printing import repr_str

        return repr_str(self)

    def delta(self):
        from geomech.core.operations.calculus import Variation

        return Variation(self)

    def t_diff(self):
        from geomech.core.operations.calculus import TimeDerivative

        return TimeDerivative(self)

    def t_integrate(self):
        from geomech.core.operations.calculus import TimeIntegral

        return TimeIntegral(self)

    def get_variation_vector(self):
        return self.delta()

    def get_tangent_vector(self):
        raise NotImplementedError

    def T(self):
        """``Transpose(self)``; ``x.T().T()`` is ``x``."""
        from geomech.core.operations import Transpose

        return Transpose(self)

    # --- flag properties (read from flags if present, else False) ---

    def _flag(self, name):
        flags = getattr(self, "flags", None)
        return getattr(flags, name, False) if flags else False

    @property
    def is_constant(self):
        return self._flag("is_constant")

    @property
    def is_zero(self):
        return self._flag("is_zero")

    @property
    def is_ones(self):
        return self._flag("is_ones")

    @property
    def is_unit_norm(self):
        return self._flag("is_unit_norm")

    @property
    def is_numeric(self):
        return self._flag("is_numeric")

    @property
    def is_symmetric(self):
        return self._flag("is_symmetric")

    @property
    def is_skew_symmetric(self):
        return self._flag("is_skew_symmetric")

    @property
    def is_orthogonal(self):
        return self._flag("is_orthogonal")

    @property
    def is_manifold(self):
        return self._flag("is_manifold")

    def has(self, elem):
        return str(elem) == str(self)


def _raise_constant_integral(expr):
    """``t_integrate`` of a constant leaf: its antiderivative ``c·t`` needs a time symbol."""
    name = _operand_str(expr)[1:-1]
    raise AlgebraicError(
        f"t_integrate: cannot integrate the constant {type(expr).__name__.lower()}"
        f" '{name}' in time: ∫ {name} dt = {name}·t, and the tree has no time symbol. Only"
        " terms that are time derivatives (d/dt(x)) integrate to an expression"
    )


# ---------------------------------------------------------------------------
# Scalar
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False)
class ScalarExpr(Expr):
    """Base for scalar-typed expressions."""

    @property
    def type(self):
        return ExprType.SCALAR

    @property
    def size(self):
        return (1,)


@dataclass(eq=False, repr=False)
class Scalar(ScalarExpr):
    name: str = None
    value: Any = None
    attr: list[str] | None = None
    flags: ExprFlags = field(init=False, repr=False, default=None)

    def __post_init__(self):
        self.flags = ExprFlags.from_attr_list(self.attr, self.name)
        if self.value is not None:
            self.flags.is_numeric = True

    def __str__(self):
        return self.name

    def delta(self):
        if self.is_constant or self.is_numeric:
            return Zero
        else:
            from geomech.core.operations import Variation

            return Variation(self)

    def t_diff(self):
        if self.is_constant or self.is_numeric:
            return Zero
        from geomech.core.operations.calculus import TimeDerivative

        return TimeDerivative(self)

    def t_integrate(self):
        if self.is_constant:
            _raise_constant_integral(self)
        from geomech.core.operations.calculus import TimeIntegral

        return TimeIntegral(self)

    # support Scalar(s='a') via custom __init__
    def __init__(self, s=None, *, name=None, value=None, attr=None, size=None):
        # 's' is the legacy positional arg; 'name' is the dataclass field
        self.name = s if s is not None else name
        self.value = value
        self.attr = _attr_list(attr, self.name)  # never alias the caller's list
        self.__post_init__()


# ---------------------------------------------------------------------------
# Vector
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False)
class VectorExpr(Expr):
    """Base for vector-typed expressions."""

    @property
    def type(self):
        return ExprType.VECTOR

    def dot(self, other):
        from geomech.core.operations import Dot

        return Dot(self, other)

    def cross(self, other):
        from geomech.core.operations import Cross

        return Cross(self, other)


@dataclass(eq=False, repr=False)
class Vector(VectorExpr):
    name: str = None
    size: tuple = (3,)
    value: Any = None
    attr: list[str] | None = None
    flags: ExprFlags = field(init=False, repr=False, default=None)

    def __post_init__(self):
        self.flags = ExprFlags.from_attr_list(self.attr, self.name)
        if self.value is None:
            self.value = np.empty(self.size, dtype="object")
        else:
            self.size = self.value.shape

    def __str__(self):
        return self.name

    def delta(self):
        if self.is_ones or self.is_zero or self.is_constant:
            return Vector("0", attr=["Constant", "Zero"])
        else:
            from geomech.core.operations import Variation

            return Variation(self)

    def t_diff(self):
        if self.is_constant:
            return Vector(s="0", size=self.size, attr=["Constant", "Zero"])
        from geomech.core.operations.calculus import TimeDerivative

        return TimeDerivative(self)

    def get_variation_vector(self):
        return self.delta()

    def t_integrate(self):
        if self.is_constant:
            _raise_constant_integral(self)
        from geomech.core.operations.calculus import TimeIntegral

        return TimeIntegral(self)

    # support Vector(s='x') via custom __init__
    def __init__(self, s=None, *, name=None, size=(3,), value=None, attr=None):
        self.name = s if s is not None else name
        self.size = size
        self.value = value
        self.attr = _attr_list(attr, self.name)  # never alias the caller's list
        self.__post_init__()


class TSO3(Vector):
    """Tangent space of SO3 manifold."""

    def __init__(self, s, *, SO3=None):
        super().__init__(s)
        self.SO3 = SO3
        if self.attr is None:
            self.attr = []
        self.attr.append("TangentVector")


class TS2(Vector):
    """Tangent space of S2 manifold."""

    def __init__(self, s, *, S2=None):
        super().__init__(s)
        self.S2 = S2
        if self.attr is None:
            self.attr = []
        self.attr.append("TangentVector")


class S2(Vector):
    """S2 manifold (unit sphere in R3, ||q|| = 1)."""

    def __init__(self, s=None, *, size=(3,), value=None, attr=None):
        attr = _attr_list(attr, s) or []  # copy: the caller's list is not mutated
        attr.append("Manifold")
        attr.append("UnitNorm")
        super().__init__(s, size=size, value=value, attr=attr)
        self.manifold_info = ManifoldInfo(
            tangent_vector_name="\\omega_{" + self.name + "}",
            variation_vector_name="\\xi_{" + self.name + "}",
        )

    def delta(self):
        from geomech.core.operations import Cross

        return Cross(self.get_variation_vector(), self)

    def get_tangent_vector(self):
        return TS2(self.manifold_info.tangent_vector_name, S2=self)

    def get_variation_vector(self):
        return TS2(self.manifold_info.variation_vector_name, S2=self)

    def t_diff(self):
        from geomech.core.operations import Cross

        return Cross(self.get_tangent_vector(), self)


# ---------------------------------------------------------------------------
# Matrix
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False)
class MatrixExpr(Expr):
    """Base for matrix-typed expressions."""

    @property
    def type(self):
        return ExprType.MATRIX


@dataclass(eq=False, repr=False)
class Matrix(MatrixExpr):
    name: str = None
    size: tuple = (3, 3)
    value: Any = None
    attr: list[str] | None = None
    flags: ExprFlags = field(init=False, repr=False, default=None)

    def __post_init__(self):
        self.flags = ExprFlags.from_attr_list(self.attr, self.name)
        if self.value is None:
            self.value = np.empty(self.size, dtype="object")

    def __str__(self):
        return self.name

    def delta(self):
        if self.is_ones or self.is_zero or self.is_constant:
            return Matrix("O", attr=["Constant", "Zero"])
        else:
            from geomech.core.operations import Variation

            return Variation(self)

    def t_diff(self):
        if self.is_constant:
            return Matrix(s="0", size=self.size, attr=["Constant", "Zero"])
        from geomech.core.operations.calculus import TimeDerivative

        return TimeDerivative(self)

    def t_integrate(self):
        if self.is_constant:
            _raise_constant_integral(self)
        from geomech.core.operations.calculus import TimeIntegral

        return TimeIntegral(self)

    # support Matrix(s='M') via custom __init__
    def __init__(self, s=None, *, name=None, size=(3, 3), value=None, attr=None):
        self.name = s if s is not None else name
        self.size = size
        self.value = value
        self.attr = _attr_list(attr, self.name) or []  # never alias the caller's list
        self.__post_init__()


class SkewSymmMatrix(Matrix):
    """Skew-symmetric matrix leaf (Xᵀ = −X); carries the ``SkewSymmetricMatrix`` flag."""

    def __init__(self, s=None, *, size=(3, 3), value=None, attr=None):
        # copy (the caller's list is not mutated) and append before super().__init__ so
        # that __post_init__ sees the flag when it builds ExprFlags
        attr = _attr_list(attr, s) or []
        attr.append("SkewSymmetricMatrix")
        super().__init__(s, size=size, value=value, attr=attr)


class SO3(Matrix):
    """SO(3) rotation matrix manifold (R^T R = I, det(R) = 1)."""

    def __init__(self, s=None, *, size=(3, 3), value=None, attr=None):
        attr = _attr_list(attr, s) or []  # copy: the caller's list is not mutated
        attr.append("Manifold")
        attr.append("OrthogonalMatrix")
        super().__init__(s, size=size, value=value, attr=attr)
        self.manifold_info = ManifoldInfo(
            tangent_vector_name="\\Omega_{" + self.name + "}",
            variation_vector_name="\\eta_{" + self.name + "}",
        )

    def delta(self):
        from geomech.core.operations import Hat, MMMul

        return MMMul(self, Hat(self.get_variation_vector()))

    def get_tangent_vector(self):
        return TSO3(self.manifold_info.tangent_vector_name, SO3=self)

    def get_variation_vector(self):
        return TSO3(self.manifold_info.variation_vector_name, SO3=self)

    def t_diff(self):
        from geomech.core.operations import Hat, MMMul

        return MMMul(self, Hat(self.get_tangent_vector()))


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

Zero = Scalar("0", value=0, attr=["Constant", "Zero"])
One = Scalar("1", value=1, attr=["Constant", "Ones"])
ZeroVector = Vector(s="0v", attr=["Constant", "Zero"])
ZeroMatrix = Matrix("0", attr=["Constant", "Zero"])
IdentityMatrix = Matrix("I", attr=["Constant", "Identity"])
O = ZeroMatrix
I = IdentityMatrix


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def Number(value):
    """Anonymous numeric constant ``(value)``.

    ``value`` is an ``int`` or ``float`` (numpy integer and floating scalars are converted
    to the equal Python number) or a ``str`` name; ``bool`` is rejected.

    ``Number(0)`` carries the ``Zero`` flag and ``Number(1)`` the ``Ones`` flag, so the
    zero/identity rules of ``simplify`` see them like ``Zero``/``One``.
    """
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(
            f"Number expects int, float, or str, got bool {value}; write Number(1) or Number(0)"
        )
    # numpy scalars become the equal Python number, so that str() is the same as for
    # the literal (np.int64(2) -> "(2)", np.float64(0.5) -> "(0.5)")
    if isinstance(value, np.integer):
        value = int(value)
    elif isinstance(value, np.floating):
        value = float(value)
    if isinstance(value, (int, float)):
        attr = ["Constant"]
        if value == 0:
            attr.append("Zero")
        elif value == 1:
            attr.append("Ones")
        return Scalar("(" + str(value) + ")", value=value, attr=attr)
    elif isinstance(value, str):
        return Scalar(s=value, attr=["Constant"])
    else:
        raise TypeError(f"Number expects int, float, or str, got {type(value).__name__}")


def getScalars(x, attr=None):
    if isinstance(x, list):
        variables = x
    elif isinstance(x, str):
        variables = x.split()
    else:
        return None
    return tuple(Scalar(v, attr=attr) for v in variables)


def getVectors(x, attr=None):
    if isinstance(x, list):
        variables = x
    elif isinstance(x, str):
        variables = x.split()
    else:
        return None
    return tuple(Vector(v, attr=attr) for v in variables)


def getMatrices(x):
    if isinstance(x, list):
        vars_ = x
    elif isinstance(x, str):
        vars_ = x.split()
    else:
        return None
    return tuple(Matrix(v) for v in vars_)
