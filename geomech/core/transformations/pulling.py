"""Pull scalar factors out of vector/matrix operations.

Bottom-up traversal over every node class (``map_children`` in ``transformations.utils``):
the operands are pulled first, then the rules below are tried at the node. O(n), one pass;
``full_simplify`` repeats it to a fixed point.

Rules (name, identity; none needs a hypothesis beyond the operand classes shown; ``s``,
``t`` are scalars, ``SVMul(v, s)`` is printed ``vs``, a row is any row vector, see
``is_row_vector``):

  pull.dot_scaled            Dot(a s, b t) = s (t Dot(a, b)); Dot(a s, b) = s Dot(a, b);
                             Dot(a, b t) = t Dot(a, b)
  pull.svmul_nested          (v s) t = v (s t)
  pull.cross_scaled          Cross(a s, b t) = Cross(a, b) (s t); Cross(a s, b) =
                             Cross(a, b) s; Cross(a, b t) = Cross(a, b) t
  pull.mvmul_scaled          (A s)(v t) = (A v)(s t); (A s) v = (A v) s; A (v t) = (A v) t
  pull.row_mvmul_scaled      (rᵀ s)(A t) = (rᵀ A)(s t); (rᵀ s) A = (rᵀ A) s;
                             rᵀ (A t) = (rᵀ A) t        a row times a matrix
  pull.smmul_nested          (M s) t = M (s t)
  pull.mmmul_scaled          (A s)(B t) = (A B)(s t); (A s) B = (A B) s; A (B t) = (A B) t
  pull.vvmul_scaled_row      (rᵀ s) b = (rᵀ b) s        inner product, the left factor is a
                                                        scaled row ``SVMul(row, s)``
  pull.vvmul_scaled          rᵀ (b t) = (rᵀ b) t (inner product, scaled column);
                             (a s) rᵀ = (a rᵀ) s and a (rᵀ t) = (a rᵀ) t (outer product)
  pull.transpose_scaled      (v s)ᵀ = vᵀ s; (M s)ᵀ = Mᵀ s
  pull.hat_scaled            Hat(v s) = Hat(v) s
  pull.vee_scaled            Vee(M s) = Vee(M) s

An inner product ``(aᵀ s) b`` (a scaled row ``SVMul(aᵀ, s)``, which ``a.T() - b.T()``
builds for its second addend, or ``(a s)ᵀ``, which becomes ``aᵀ s`` first) becomes
``aᵀb s``, so it merges with ``s (aᵀ b)``. The calculus nodes have no rule: ``d/dt(v s)``
is not ``d/dt(v) s`` unless ``s`` is constant, so their operands are only recursed into.
"""

from geomech.core.base.types import ExprType
from geomech.core.operations.geometry import Cross, Dot, Hat, Transpose, Vee, is_row_vector
from geomech.core.operations.multiplication import MMMul, Mul, MVMul, SMMul, SVMul, VVMul
from geomech.core.transformations.utils import map_children
from geomech.utils.trace import rule


def _scaled_vector(v):
    """``(vector, scalar)`` of ``SVMul(vector, scalar)``, else ``(v, None)``."""
    if isinstance(v, SVMul):
        return v.left, v.right
    return v, None


def _scaled_matrix(m):
    """``(matrix, scalar)`` of ``SMMul(matrix, scalar)``, else ``(m, None)``."""
    if isinstance(m, SMMul):
        return m.left, m.right
    return m, None


def _times(s, t):
    """The scalar ``s t`` with ``None`` meaning "no factor"."""
    if s is None:
        return t
    if t is None:
        return s
    return Mul(s, t)


def _pull_vvmul(expr):
    """Hoist the scalars of scaled factors out of an inner or outer product."""
    l, r = expr.left, expr.right
    if expr.type == ExprType.SCALAR:
        # inner product rowᵀ · column
        if isinstance(l, SVMul) and is_row_vector(l.left):
            r0, t = _scaled_vector(r)
            return rule("pull.vvmul_scaled_row", expr, Mul(VVMul(l.left, r0), _times(l.right, t)))
        if isinstance(r, SVMul):
            return rule("pull.vvmul_scaled", expr, Mul(VVMul(l, r.left), r.right))
        return expr
    # outer product column · rowᵀ
    l0, s = _scaled_vector(l)
    r0, t = _scaled_vector(r) if isinstance(r, SVMul) and is_row_vector(r.left) else (r, None)
    if s is None and t is None:
        return expr
    return rule("pull.vvmul_scaled", expr, SMMul(VVMul(l0, r0), _times(s, t)))


def _pull_row_mvmul(expr):
    """``rᵀ A`` with a scaled row or a scaled matrix: hoist the scalars (a scaled row)."""
    l, r = expr.left, expr.right
    l0, s = _scaled_vector(l)
    r0, t = _scaled_matrix(r)
    if s is None and t is None:
        return expr
    return rule("pull.row_mvmul_scaled", expr, SVMul(MVMul(l0, r0), _times(s, t)))


def pull(expr):
    """Hoist scalar factors out of products, bottom up (rules in the module docstring)."""
    expr = map_children(expr, pull, "pull")
    return _at_node(expr)


def _at_node(expr):
    """The rules of the module docstring at the root of *expr* (operands already pulled).

    One function per node class, picked from ``_PULL`` by the class of *expr*; sums,
    scalar products, calculus nodes and leaves have no rule.
    """
    handler = _PULL.get(type(expr))
    return expr if handler is None else handler(expr)


def _pull_dot(expr):
    """``Dot(vec, vec)``: the scalars of scaled operands come out."""
    l, r = expr.left, expr.right
    match (l, r):
        case (SVMul(), SVMul()):
            return rule("pull.dot_scaled", expr, Mul(l.right, Mul(r.right, Dot(l.left, r.left))))
        case (SVMul(), _):
            return rule("pull.dot_scaled", expr, Mul(l.right, Dot(l.left, r)))
        case (_, SVMul()):
            return rule("pull.dot_scaled", expr, Mul(r.right, Dot(l, r.left)))
        case _:
            return expr


def _pull_svmul(expr):
    """``(v s) t = v (s t)``."""
    if isinstance(expr.left, SVMul):
        l = expr.left
        return rule("pull.svmul_nested", expr, SVMul(l.left, Mul(l.right, expr.right)))
    return expr


def _pull_cross(expr):
    """``Cross(vec, vec)``: the scalars of scaled operands come out."""
    l, r = expr.left, expr.right
    match (l, r):
        case (SVMul(), SVMul()):
            return rule(
                "pull.cross_scaled",
                expr,
                SVMul(Cross(l.left, r.left), Mul(l.right, r.right)),
            )
        case (SVMul(), _):
            return rule("pull.cross_scaled", expr, SVMul(Cross(l.left, r), l.right))
        case (_, SVMul()):
            return rule("pull.cross_scaled", expr, SVMul(Cross(l, r.left), r.right))
        case _:
            return expr


def _pull_mvmul(expr):
    """A matrix times a vector (or a row times a matrix): scaled operands come out."""
    if expr.left.type == ExprType.VECTOR:
        return _pull_row_mvmul(expr)
    l, r = expr.left, expr.right
    match (l, r):
        case (SMMul(), SVMul()):
            return rule(
                "pull.mvmul_scaled",
                expr,
                SVMul(MVMul(l.left, r.left), Mul(l.right, r.right)),
            )
        case (SMMul(), _):
            return rule("pull.mvmul_scaled", expr, SVMul(MVMul(l.left, r), l.right))
        case (_, SVMul()):
            return rule("pull.mvmul_scaled", expr, SVMul(MVMul(l, r.left), r.right))
        case _:
            return expr


def _pull_smmul(expr):
    """``(M s) t = M (s t)``."""
    if isinstance(expr.left, SMMul):
        l = expr.left
        return rule("pull.smmul_nested", expr, SMMul(l.left, Mul(l.right, expr.right)))
    return expr


def _pull_mmmul(expr):
    """A matrix product: the scalars of scaled operands come out."""
    l, r = expr.left, expr.right
    match (l, r):
        case (SMMul(), SMMul()):
            return rule(
                "pull.mmmul_scaled",
                expr,
                SMMul(MMMul(l.left, r.left), Mul(l.right, r.right)),
            )
        case (SMMul(), _):
            return rule("pull.mmmul_scaled", expr, SMMul(MMMul(l.left, r), l.right))
        case (_, SMMul()):
            return rule("pull.mmmul_scaled", expr, SMMul(MMMul(l, r.left), r.right))
        case _:
            return expr


def _pull_transpose(expr):
    """``(v s)ᵀ = vᵀ s``, ``(M s)ᵀ = Mᵀ s``."""
    inner = expr.expr
    if isinstance(inner, SVMul):
        return rule("pull.transpose_scaled", expr, SVMul(Transpose(inner.left), inner.right))
    if isinstance(inner, SMMul):
        return rule("pull.transpose_scaled", expr, SMMul(Transpose(inner.left), inner.right))
    return expr


def _pull_hat(expr):
    """``Hat(v s) = Hat(v) s``."""
    inner = expr.expr
    if isinstance(inner, SVMul):
        return rule("pull.hat_scaled", expr, SMMul(Hat(inner.left), inner.right))
    return expr


def _pull_vee(expr):
    """``Vee(M s) = Vee(M) s``."""
    inner = expr.expr
    if isinstance(inner, SMMul):
        return rule("pull.vee_scaled", expr, SVMul(Vee(inner.left), inner.right))
    return expr


#: ``_at_node``'s handler for each node class with a pull rule.
#: Looked up by the exact class (``type(expr)``), not along the MRO: a subclass of a
#: node class is not dispatched and is returned unchanged, like a leaf. No such subclass exists;
#: one added later needs its own entry here (pinned in test_phase5b_structure.py).
_PULL = {
    Dot: _pull_dot,
    SVMul: _pull_svmul,
    Cross: _pull_cross,
    MVMul: _pull_mvmul,
    SMMul: _pull_smmul,
    MMMul: _pull_mmmul,
    VVMul: _pull_vvmul,
    Transpose: _pull_transpose,
    Hat: _pull_hat,
    Vee: _pull_vee,
}
