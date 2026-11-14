from __future__ import annotations

from geomech.core.base.types import ExprType

# ---------------------------------------------------------------------------
# Base — shared by all operation nodes
# ---------------------------------------------------------------------------


class _BaseMixin:
    """Shared base for all operation nodes."""

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
        from geomech.core.base.expressions import Zero, ZeroMatrix, ZeroVector

        if self.is_constant:
            match self.type:
                case ExprType.SCALAR:
                    return Zero
                case ExprType.VECTOR:
                    return ZeroVector
                case ExprType.MATRIX:
                    return ZeroMatrix
        return type(self)(*[n.delta() for n in self.nodes])

    def t_diff(self):
        from geomech.core.base.expressions import Zero, ZeroMatrix, ZeroVector

        if self.is_constant:
            match self.type:
                case ExprType.SCALAR:
                    return Zero
                case ExprType.VECTOR:
                    return ZeroVector
                case ExprType.MATRIX:
                    return ZeroMatrix
        return type(self)(*[n.t_diff() for n in self.nodes])


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
        from geomech.core.base.expressions import Zero, ZeroMatrix, ZeroVector

        if self.is_constant:
            match self.type:
                case ExprType.SCALAR:
                    return Zero
                case ExprType.VECTOR:
                    return ZeroVector
                case ExprType.MATRIX:
                    return ZeroMatrix
        op = lambda node: getattr(node, op_name)()
        if self.left.is_constant:
            return type(self)(self.left, op(self.right))
        if self.right.is_constant:
            return type(self)(op(self.left), self.right)
        # Product rule: op(l*r) = op(l)*r + l*op(r)
        term1 = type(self)(op(self.left), self.right)
        term2 = type(self)(self.left, op(self.right))
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
        from geomech.core.base.expressions import Zero, ZeroMatrix, ZeroVector

        if self.is_constant:
            match self.type:
                case ExprType.SCALAR:
                    return Zero
                case ExprType.VECTOR:
                    return ZeroVector
                case ExprType.MATRIX:
                    return ZeroMatrix
        return type(self)(self.expr.delta())

    def t_diff(self):
        from geomech.core.base.expressions import Zero, ZeroMatrix, ZeroVector

        if self.is_constant:
            match self.type:
                case ExprType.SCALAR:
                    return Zero
                case ExprType.VECTOR:
                    return ZeroVector
                case ExprType.MATRIX:
                    return ZeroMatrix
        return type(self)(self.expr.t_diff())


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

    def t_integrate(self):
        from geomech.core.base.expressions import Zero, ZeroMatrix, ZeroVector

        if self.is_constant:
            match self.type:
                case ExprType.SCALAR:
                    return Zero
                case ExprType.VECTOR:
                    return ZeroVector
                case ExprType.MATRIX:
                    return ZeroMatrix
        return type(self)(self.expr.t_integrate())
