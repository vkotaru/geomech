from __future__ import annotations

from dataclasses import dataclass, field

from geomech.core.base.expressions import (
    MatrixExpr,
    ScalarExpr,
    VectorExpr,
    _node_operands,
)
from geomech.core.base.types import ExprType
from geomech.core.operations.mixins import _NaryMixin
from geomech.utils.errors import ExpressionMismatchError, SizeMismatchError, _operand_str

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _check_sizes(nodes, op_name: str):
    """Validate that all operands with a known size are compatible."""
    ref_size = None
    for n in nodes:
        s = getattr(n, "size", None)
        if s is not None:
            if ref_size is None:
                ref_size = s
            elif s != ref_size:
                raise SizeMismatchError(op_name, ref_size, s)


def _flatten_nodes(args, expected_type: type, expr_type: ExprType) -> list:
    """Flatten nested same-type adds and validate types.

    Numbers are wrapped as ``Number`` (``Add(x, 1)``); any other non-expression operand
    raises ``TypeError``. An operand of the wrong type raises ``ExpressionMismatchError``
    naming it and the first operand as written: for ``(p + q) + a`` that is the sum
    ``p + q``, not its first addend ``p``.

    An addition with no operands is refused (``TypeError``): the zero of each type is a
    leaf (``Zero``, ``ZeroVector``, ``ZeroMatrix``), and an empty sum would make the
    ``all(...)``-based predicates (``is_constant``, ``is_zero``, ``MAdd.is_symmetric``,
    ``MAdd.is_skew_symmetric``) vacuously true.
    """
    op = expected_type.__name__
    nodes = []
    written = []  # the first operand as the caller passed it, before flattening

    def _take(a):
        (a,) = _node_operands(op, a, numeric=True)
        if a.type != expr_type:
            if written:
                first = written[0]
                raise ExpressionMismatchError(op, first.type, a.type, first, a)
            raise ExpressionMismatchError(
                msg=f"{op} expects {expr_type.name} operands, got {a.type.name} {_operand_str(a)}"
            )
        return a

    for arg in args:
        if isinstance(arg, (list, tuple)):
            for a in arg:
                nodes.append(_take(a))
                written[:] = written or nodes[-1:]
        elif isinstance(arg, expected_type):
            nodes.extend(arg.nodes)
            written[:] = written or [arg]
        else:
            nodes.append(_take(arg))
            written[:] = written or nodes[-1:]
    if not nodes:
        raise TypeError(
            f"{op}() needs at least one operand; use the zero leaf of the type for an empty sum"
        )
    # Size validation for vector/matrix additions
    if expr_type in (ExprType.VECTOR, ExprType.MATRIX):
        _check_sizes(nodes, op)
    return nodes


# ---------------------------------------------------------------------------
# N-ary addition
# ---------------------------------------------------------------------------


@dataclass(eq=False, repr=False)
class Add(_NaryMixin, ScalarExpr):
    """Scalar addition (n-ary)."""

    nodes: list = field(default_factory=list)

    def __init__(self, *args):
        self.nodes = _flatten_nodes(args, Add, ExprType.SCALAR)


@dataclass(eq=False, repr=False)
class VAdd(_NaryMixin, VectorExpr):
    """Vector addition (n-ary)."""

    nodes: list = field(default_factory=list)

    def __init__(self, *args):
        self.nodes = _flatten_nodes(args, VAdd, ExprType.VECTOR)


@dataclass(eq=False, repr=False)
class MAdd(_NaryMixin, MatrixExpr):
    """Matrix addition (n-ary)."""

    nodes: list = field(default_factory=list)

    def __init__(self, *args):
        self.nodes = _flatten_nodes(args, MAdd, ExprType.MATRIX)

    @property
    def is_symmetric(self):
        return all(n.is_symmetric for n in self.nodes)

    @property
    def is_skew_symmetric(self):
        """A sum of skew-symmetric matrices is skew-symmetric."""
        return all(n.is_skew_symmetric for n in self.nodes)
