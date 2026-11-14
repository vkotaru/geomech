"""Algebraic identities for vector/matrix geometry operations (``vector_rules``).

One bottom-up pass over every node class (``map_children`` in ``transformations.utils``):
the operands of a node are rewritten first, then the rules below are tried at the node, so
an identity fires wherever its node sits, under a ``Dot``, ``Hat``, ``MVMul``, ``Transpose``
or a calculus operator included. Equality of two operands is ``==`` (type + printed name),
so two ``S2("q")`` objects are the same point.

Rules (name, identity, hypothesis the code checks):

  vector.dot_zero                  Dot(x, 0) = Dot(0, x) = 0
                                   an operand is flagged zero
  vector.dot_cross_orthogonal      Dot(x, Cross(x, y)) = Dot(x, Cross(y, x)) = 0,
                                   Dot(z, Cross(y, y)) = 0, and the mirror with the Cross
                                   on the left; x equal to a Cross operand, or the two
                                   Cross operands equal
  vector.dot_unit_norm             Dot(q, q) = 1
                                   both operands flagged unit norm (S2 points) and equal
  vector.dot_tangent_orthogonal    Dot(w, q) = Dot(q, w) = 0
                                   w a TS2 leaf (ω or ξ) whose S2 point equals q by name
  vector.cross_self                Cross(x, x) = 0
                                   operands equal
  vector.cross_zero                Cross(x, 0) = Cross(0, x) = 0
                                   an operand is flagged zero
  vector.cross_scaled_self         Cross(a, a s) = Cross(a s, a) = 0
                                   the SVMul's vector equals the other operand
  vector.s2_tangent_recovery       Cross(q, Cross(a, q)) = a and Cross(Cross(q, a), q) = a
                                   q an S2 point (unit norm), a tangent to q: a TS2 leaf
                                   of q (by name) or a Cross with q as an operand
  vector.s2_tangent_negated        Cross(q, Cross(q, a)) = -a and Cross(Cross(a, q), q) = -a
                                   same hypothesis as vector.s2_tangent_recovery
  vector.s2_double_cross           Cross(a, Cross(a, q)) = -q (a·a)
                                   q an S2 point, a tangent to q (as above)
  vector.hat_zero                  Hat(0) = 0 (matrix)
                                   the operand is flagged zero

Before 5.1 the pass stopped at ``Dot`` operands, ``Hat``, ``MVMul`` and the calculus
nodes, so ``Dot(v, q×(ω×q))``, ``Hat(a×a)`` or ``J(a×a)`` survived ``full_simplify``.

The tangency test does not look inside ``TimeDerivative``/``Variation``: ``ω̇ ⊥ q`` holds
on S² but needs a constraint system to prove, so ``Cross(q, Cross(ω̇, q))`` is kept.
"""

from geomech.core.base.constants import NEG_ONE
from geomech.core.base.expressions import S2, TS2, One, Zero, ZeroMatrix, ZeroVector
from geomech.core.operations.geometry import Cross, Dot, Hat
from geomech.core.operations.multiplication import Mul, SVMul
from geomech.core.transformations.utils import map_children
from geomech.utils.trace import rule


def _is_orthogonal_dot_cross(vec, cross_expr):
    """Check if dot(vec, cross(a, b)) is zero by orthogonality.

    Returns True when:
      - vec == a or vec == b  (x perpendicular to x x y)
      - a == b                (x x x = 0)
    """
    return vec == cross_expr.left or vec == cross_expr.right or cross_expr.left == cross_expr.right


def _is_tangent_orthogonal(a, b):
    """Check if dot(a, b) = 0 due to manifold tangent-point orthogonality.

    On S2: tangent vectors (ω, ξ) are orthogonal to the manifold point q.
    Returns True when one is a TS2 and the other is its parent S2. The parent is
    compared by name (``==`` on leaves is type + printed name), as everywhere else in
    the engine, so two ``S2("q")`` objects are the same point.
    """
    if isinstance(a, TS2) and isinstance(b, S2) and a.S2 == b:
        return True
    if isinstance(b, TS2) and isinstance(a, S2) and b.S2 == a:
        return True
    return False


def _try_s2_cross_reduction(expr):
    """Simplify Cross expressions using S2 unit-norm + tangent orthogonality.

    BAC-CAB: a × (b × c) = b(a·c) - c(a·b)

    When q is S2 (‖q‖=1) and a is known tangent (a ⊥ q):
      q × (a × q) = a(q·q) - q(q·a) = a

    Only fires when tangency is structurally provable (direct TS2 or
    Cross(*, q)). For other cases (e.g. d/dt(ω)), the constraint
    system (see GitHub issue #3) is needed.
    """
    l, r = expr.left, expr.right

    # Pattern: Cross(q, Cross(a, q)) where q is S2, a is tangent to q
    if isinstance(l, S2) and isinstance(r, Cross) and r.right == l:
        a = r.left
        if _is_tangent_to_s2(a, l):
            return rule("vector.s2_tangent_recovery", expr, a)

    # Pattern: Cross(q, Cross(q, a)) = -a (when a is tangent)
    if isinstance(l, S2) and isinstance(r, Cross) and r.left == l:
        a = r.right
        if _is_tangent_to_s2(a, l):
            return rule("vector.s2_tangent_negated", expr, SVMul(a, NEG_ONE))

    # Pattern: Cross(Cross(q, a), q) = a (when a is tangent)
    if isinstance(l, Cross) and isinstance(r, S2) and l.left == r:
        a = l.right
        if _is_tangent_to_s2(a, r):
            return rule("vector.s2_tangent_recovery", expr, a)

    # Pattern: Cross(Cross(a, q), q) = -a (when a is tangent)
    if isinstance(l, Cross) and isinstance(r, S2) and l.right == r:
        a = l.left
        if _is_tangent_to_s2(a, r):
            return rule("vector.s2_tangent_negated", expr, SVMul(a, NEG_ONE))

    # Pattern: Cross(a, Cross(a, q)) = -q*(a·a) when a is tangent to q
    if isinstance(r, Cross) and r.right is not None and isinstance(r.right, S2) and l == r.left:
        q = r.right
        a = l
        if _is_tangent_to_s2(a, q):
            return rule("vector.s2_double_cross", expr, SVMul(q, Mul(Dot(a, a), NEG_ONE)))

    return None


def _is_tangent_to_s2(expr, s2):
    """Check if expr is a tangent vector to the given S2 manifold.

    Returns True for TS2 vectors whose parent is the given S2 (compared by name),
    or for Cross(anything, q) which is perpendicular to q by definition.

    Does NOT recurse into TimeDerivative or Variation — while ω̇ ⊥ q
    is provable for S2, this requires a constraint system to verify
    properly (see GitHub issue #3). Without constraints, we only
    apply BAC-CAB to vectors we can structurally confirm as tangent.
    """
    # Direct TS2
    if isinstance(expr, TS2) and expr.S2 == s2:
        return True
    # Cross(a, q) is tangent to q (perpendicular to q by definition)
    if isinstance(expr, Cross):
        if expr.right == s2 or expr.left == s2:
            return True
    return False


def vector_rules(expr):
    """Apply the geometric identities listed in the module docstring, bottom up, once."""
    expr = map_children(expr, vector_rules, "vector_rules")
    return _at_node(expr)


def _at_node(expr):
    """The identities of the module docstring at the root of *expr* (operands done)."""
    match expr:
        # --- Dot product rules ---
        case Dot():
            # zero absorption: dot(0, x) or dot(x, 0)
            if expr.left.is_zero or expr.right.is_zero:
                return rule("vector.dot_zero", expr, Zero)

            match (expr.left, expr.right):
                # dot(x, cross(a, b)) — cross on right
                case (_, Cross()) if _is_orthogonal_dot_cross(expr.left, expr.right):
                    return rule("vector.dot_cross_orthogonal", expr, Zero)

                # dot(cross(a, b), x) — cross on left
                case (Cross(), _) if _is_orthogonal_dot_cross(expr.right, expr.left):
                    return rule("vector.dot_cross_orthogonal", expr, Zero)

                # dot(q, q) = 1 when unit norm
                case _ if (
                    expr.left.is_unit_norm and expr.right.is_unit_norm and expr.left == expr.right
                ):
                    return rule("vector.dot_unit_norm", expr, One)

                # dot(ω, q) = 0 — tangent vector orthogonal to manifold point
                case _ if _is_tangent_orthogonal(expr.left, expr.right):
                    return rule("vector.dot_tangent_orthogonal", expr, Zero)

                case _:
                    return expr

        # --- Cross product rules ---
        case Cross() if expr.left == expr.right:
            return rule("vector.cross_self", expr, ZeroVector)

        case Cross() if expr.left.is_zero or expr.right.is_zero:
            return rule("vector.cross_zero", expr, ZeroVector)

        # cross(a, s*a) = s*(a × a) = 0  and  cross(s*a, a) = 0
        case Cross() if isinstance(expr.right, SVMul) and expr.right.left == expr.left:
            return rule("vector.cross_scaled_self", expr, ZeroVector)

        case Cross() if isinstance(expr.left, SVMul) and expr.left.left == expr.right:
            return rule("vector.cross_scaled_self", expr, ZeroVector)

        # S2 BAC-CAB identity with tangent orthogonality
        case Cross():
            reduced = _try_s2_cross_reduction(expr)
            return expr if reduced is None else reduced

        # --- Hat rules ---
        case Hat() if expr.expr.is_zero:
            return rule("vector.hat_zero", expr, ZeroMatrix)

        case _:
            return expr
