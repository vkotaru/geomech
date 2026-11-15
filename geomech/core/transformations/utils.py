"""Structural helpers shared by the tree rewrites.

``rebuild(expr, children)`` makes a node of the same class as *expr* with new children, and
``map_children(expr, fn)`` applies *fn* to every child and rebuilds only when a child
changed. ``vector_rules``, ``pull`` and ``apply_manifold_rules`` recurse through every node
class with these two functions, so a rule written for one node fires wherever that node
sits in the tree (under a ``Dot``, ``Hat``, ``MVMul``, a calculus operator, ...).

``addends(expr)``, ``sum_of(terms, expr_type)`` and ``flatten_mul(expr)`` take sums and
product chains apart and put them back together; the simplifier, integration by parts,
``to_standard_form`` and ``linearize`` share them.
"""

from __future__ import annotations

from geomech.core.base.constants import zero_of
from geomech.core.base.expressions import MANIFOLD_POINTS, Expr
from geomech.core.base.types import ExprType
from geomech.core.operations.addition import Add, MAdd, VAdd
from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
from geomech.core.operations.geometry import Cross, Dot, Hat, Transpose, Vee
from geomech.core.operations.multiplication import MMMul, Mul, MVMul, SMMul, SVMul, VVMul


def addends(expr: Expr) -> list[Expr]:
    """The addends of a sum (``Add``, ``VAdd``, ``MAdd``) as a list; ``[expr]`` otherwise."""
    if isinstance(expr, (Add, VAdd, MAdd)):
        return list(expr.nodes)
    return [expr]


def sum_of(terms: list[Expr], expr_type: ExprType | None) -> Expr:
    """``terms`` added up as a sum of *expr_type*: the inverse of ``addends``.

    No term gives the zero of the type (``zero_of``), one term is returned as it is,
    several make an ``Add``, ``VAdd`` or ``MAdd``. A missing type counts as scalar.
    """
    if not terms:
        return zero_of(expr_type)
    if len(terms) == 1:
        return terms[0]
    match expr_type:
        case ExprType.VECTOR:
            return VAdd(*terms)
        case ExprType.MATRIX:
            return MAdd(*terms)
        case _:
            return Add(*terms)


def flatten_mul(expr: Expr) -> list[Expr]:
    """The factors of a scalar ``Mul`` chain, left to right (``[expr]`` for a non-``Mul``)."""
    factors: list[Expr] = []
    stack = [expr]
    while stack:
        node = stack.pop()
        if isinstance(node, Mul):
            stack.append(node.right)
            stack.append(node.left)
        else:
            factors.append(node)
    return factors


def calculus_node(cls, operand):
    """``cls(operand)`` for a calculus class, in eager form over a manifold point.

    A rewrite of the operand of ``TimeDerivative``/``Variation`` can reduce it to a bare
    ``S1``/``S2``/``SO3`` leaf (``d/dt(q·(q·q))`` → ``d/dt(q)``). The tree never holds that node
    (the constructor rejects it), so the leaf's own kinematics are returned instead:
    ``q.t_diff()``/``R.t_diff()`` for ``TimeDerivative`` and ``q.delta()``/``R.delta()`` for
    ``Variation``. Any other operand gives ``cls(operand)``.
    """
    if isinstance(operand, MANIFOLD_POINTS):
        if cls is TimeDerivative:
            return operand.t_diff()
        if cls is Variation:
            return operand.delta()
    return cls(operand)


def rebuild(expr, children, caller="rebuild"):
    """A node of the class of *expr* with *children* in place of its operands.

    N-ary sums take all children, the binary products and ``Dot``/``Cross`` two, the unary
    operators (``Hat``, ``Vee``, ``Transpose`` and the calculus nodes) one. A
    ``TimeDerivative``/``Variation`` whose new operand is a bare manifold leaf is built in
    eager form (see ``calculus_node``). Any other class
    raises ``TypeError`` naming *caller*: returning the old node would silently drop the
    rewrite of a child.
    """
    cls = type(expr)
    if isinstance(expr, (Add, VAdd, MAdd)):
        return cls(*children)
    if isinstance(expr, (Mul, SVMul, SMMul, MVMul, MMMul, VVMul, Dot, Cross)):
        return cls(children[0], children[1])
    if isinstance(expr, (Variation, TimeDerivative)):
        return calculus_node(cls, children[0])
    if isinstance(expr, (Hat, Vee, Transpose, TimeIntegral)):
        return cls(children[0])
    raise TypeError(
        f"{caller}: cannot rebuild a {cls.__name__} node ('{expr}') with rewritten"
        " children; add the class to geomech.core.transformations.utils.rebuild"
    )


def map_children(expr, fn, caller="map_children"):
    """*expr* with *fn* applied to each child; the same object when no child changed.

    A leaf (no children) is returned as is.
    """
    nodes = getattr(expr, "nodes", None)
    if not nodes:
        return expr
    new = [fn(n) for n in nodes]
    if all(a is b for a, b in zip(new, nodes)):
        return expr
    return rebuild(expr, new, caller)
