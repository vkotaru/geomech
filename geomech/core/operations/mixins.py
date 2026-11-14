from __future__ import annotations

from typing import TYPE_CHECKING, Any

from geomech.core.base.constants import zero_like
from geomech.core.base.types import ExprType

if TYPE_CHECKING:
    # For the type checkers the mixins are expressions (every class that uses one is an
    # Expr subclass), so ``self.type``, ``self.flags`` and the operators are known. At run
    # time they stay plain classes ahead of the Expr base in the MRO.
    from geomech.core.base.expressions import Expr as _NodeBase
else:
    _NodeBase = object

# ---------------------------------------------------------------------------
# Base — shared by all operation nodes
# ---------------------------------------------------------------------------


class _BaseMixin(_NodeBase):
    """Shared base for all operation nodes.

    ``nodes`` is the tuple of operands, set once by the node's constructor (the node
    classes are frozen dataclasses).
    """

    nodes: tuple[Any, ...]

    def has(self, elem):
        if self == elem:
            return True
        return any(n.has(elem) for n in self.nodes)

    @property
    def arity(self):
        return len(self.nodes)

    def __len__(self):
        return len(self.nodes)

    # Arithmetic operators (+, -, *, /, unary -) are defined once on ``Expr``
    # (``geomech/core/base/expressions.py``) for leaves and operation nodes alike.


# ---------------------------------------------------------------------------
# N-ary (addition)
# ---------------------------------------------------------------------------


class _NaryMixin(_BaseMixin):
    """Shared properties for n-ary addition operations."""

    @property
    def is_constant(self):
        return all(n.is_constant for n in self.nodes)

    @property
    def is_zero(self):
        return all(n.is_zero for n in self.nodes)

    def __str__(self):
        return "(" + "+".join(str(n) for n in self.nodes) + ")"

    def delta(self):
        if self.is_constant:
            return zero_like(self)
        cls: Any = type(self)
        return cls(*[n.delta() for n in self.nodes])

    def t_diff(self):
        if self.is_constant:
            return zero_like(self)
        cls: Any = type(self)
        return cls(*[n.t_diff() for n in self.nodes])


# ---------------------------------------------------------------------------
# Binary (multiplication)
# ---------------------------------------------------------------------------


class _BinaryMixin(_BaseMixin):
    """Shared properties for binary operations using unified nodes."""

    @property
    def left(self):
        return self.nodes[0]

    @property
    def right(self):
        return self.nodes[1]

    @property
    def is_constant(self):
        return self.left.is_constant and self.right.is_constant

    @property
    def is_zero(self):
        return self.left.is_zero or self.right.is_zero

    def _apply_rule(self, op_name):
        """Apply linearity rule (delta or t_diff) with product-rule logic.

        op_name: 'delta' or 't_diff'
        """
        if self.is_constant:
            return zero_like(self)
        op = lambda node: getattr(node, op_name)()
        cls: Any = type(self)
        if self.left.is_constant:
            return cls(self.left, op(self.right))
        if self.right.is_constant:
            return cls(op(self.left), self.right)
        # Product rule: op(l*r) = op(l)*r + l*op(r)
        term1 = cls(op(self.left), self.right)
        term2 = cls(self.left, op(self.right))
        from geomech.core.operations.addition import Add, MAdd, VAdd

        match self.type:
            case ExprType.SCALAR:
                return Add(term1, term2)
            case ExprType.VECTOR:
                return VAdd(term1, term2)
            case ExprType.MATRIX:
                return MAdd(term1, term2)

    def delta(self):
        return self._apply_rule("delta")

    def t_diff(self):
        return self._apply_rule("t_diff")


# ---------------------------------------------------------------------------
# Unary (fixed output type — Hat, Vee, etc.)
# ---------------------------------------------------------------------------


class _UnaryMixin(_BaseMixin):
    """Shared properties for unary operation nodes with fixed output type."""

    @property
    def expr(self):
        return self.nodes[0]

    @property
    def is_constant(self):
        return self.expr.is_constant

    @property
    def is_zero(self):
        return self.expr.is_zero

    def delta(self):
        if self.is_constant:
            return zero_like(self)
        cls: Any = type(self)
        return cls(self.expr.delta())

    def t_diff(self):
        if self.is_constant:
            return zero_like(self)
        cls: Any = type(self)
        return cls(self.expr.t_diff())


# ---------------------------------------------------------------------------
# Calculus unary (Variation, TimeDerivative, TimeIntegral, Delta)
# ---------------------------------------------------------------------------


class _CalcUnaryMixin(_UnaryMixin):
    """Type-preserving unary — delegates type to inner expression.

    Pushes t_integrate() inside by default so that e.g.
    Variation(expr).t_integrate() → Variation(expr.t_integrate()).
    Subclasses (TimeDerivative, TimeIntegral) override for cancellation.
    """

    @property
    def type(self):
        return self.expr.type

    def _identity_extra(self) -> tuple:
        """The classes from this node down to the first operand that is not a calculus
        node, and what that operand adds itself (a ``Transpose`` adds its operand's class):
        ``Variation(Scalar x)`` and ``Variation(Vector x)`` print alike but differ."""
        chain = []
        node = self.expr
        while isinstance(node, _CalcUnaryMixin):
            chain.append(type(node))
            node = node.expr
        chain.append(type(node))
        return (*chain, *node._identity_extra())

    def t_integrate(self):
        if self.is_constant:
            return zero_like(self)
        cls: Any = type(self)
        return cls(self.expr.t_integrate())
