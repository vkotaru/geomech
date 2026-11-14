"""Pull scalar factors out of vector/matrix operations.

Bottom-up traversal: recurse into children first, then handle the
current node.  O(n) single pass — no has_nested_scalars predicate needed.

An inner product ``(aᵀ s) b`` (a scaled row ``SVMul(aᵀ, s)``, which ``a.T() - b.T()``
builds for its second addend) becomes ``aᵀb s``, so it merges with ``s (aᵀ b)``.
"""

from geomech.core.operations.addition import Add, MAdd, VAdd
from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
from geomech.core.operations.geometry import Cross, Dot, is_row_vector
from geomech.core.operations.multiplication import MMMul, Mul, MVMul, SMMul, SVMul, VVMul


def _pull_vvmul(expr):
    """``(rᵀ s) b`` → ``(rᵀ b) s``: hoist the scalar of a scaled row out of an inner product.

    Only the row form ``SVMul(row, s)`` on the left is rewritten; every other ``VVMul`` is
    returned unchanged (in particular ``(c s)ᵀ b``, which the pinned pipelines build).
    """
    left = pull(expr.left)
    if isinstance(left, SVMul) and is_row_vector(left.left):
        return Mul(VVMul(left.left, expr.right), left.right)
    return expr


def pull(expr):
    match expr:
        # --- N-ary (additions) ---
        case Add(nodes=nodes):
            return Add(*[pull(n) for n in nodes])
        case VAdd(nodes=nodes):
            return VAdd(*[pull(n) for n in nodes])
        case MAdd(nodes=nodes):
            return MAdd(*[pull(n) for n in nodes])

        # --- Scalar * Scalar ---
        case Mul():
            return Mul(pull(expr.left), pull(expr.right))

        # --- Dot(vec, vec) → scalar ---
        case Dot():
            l, r = pull(expr.left), pull(expr.right)
            match (l, r):
                case (SVMul(), SVMul()):
                    return Mul(l.right, Mul(r.right, Dot(l.left, r.left)))
                case (SVMul(), _):
                    return Mul(l.right, Dot(l.left, r))
                case (_, SVMul()):
                    return Mul(r.right, Dot(l, r.left))
                case _:
                    return Dot(l, r)

        # --- Scalar * Vector ---
        case SVMul():
            l, r = pull(expr.left), pull(expr.right)
            match l:
                case SVMul():
                    return SVMul(l.left, Mul(l.right, r))
                case _:
                    return SVMul(l, r)

        # --- Cross(vec, vec) → vector ---
        case Cross():
            l, r = pull(expr.left), pull(expr.right)
            match (l, r):
                case (SVMul(), SVMul()):
                    return SVMul(Cross(l.left, r.left), Mul(l.right, r.right))
                case (SVMul(), _):
                    return SVMul(Cross(l.left, r), l.right)
                case (_, SVMul()):
                    return SVMul(Cross(l, r.left), r.right)
                case _:
                    return Cross(l, r)

        # --- Matrix * Vector → vector ---
        case MVMul():
            l, r = pull(expr.left), pull(expr.right)
            match (l, r):
                case (SMMul(), SVMul()):
                    return SVMul(MVMul(l.left, r.left), Mul(l.right, r.right))
                case (SMMul(), _):
                    return SVMul(MVMul(l.left, r), l.right)
                case (_, SVMul()):
                    return SVMul(MVMul(l, r.left), r.right)
                case _:
                    return MVMul(l, r)

        # --- Scalar * Matrix ---
        case SMMul():
            l, r = pull(expr.left), pull(expr.right)
            match l:
                case SMMul():
                    return SMMul(l.left, Mul(l.right, r))
                case _:
                    return SMMul(l, r)

        # --- Matrix * Matrix → matrix ---
        case MMMul():
            l, r = pull(expr.left), pull(expr.right)
            match (l, r):
                case (SMMul(), SMMul()):
                    return SMMul(MMMul(l.left, r.left), Mul(l.right, r.right))
                case (SMMul(), _):
                    return SMMul(MMMul(l.left, r), l.right)
                case (_, SMMul()):
                    return SMMul(MMMul(l, r.left), r.right)
                case _:
                    return MMMul(l, r)

        # --- Vector * Vector (inner or outer product) ---
        case VVMul():
            return _pull_vvmul(expr)

        # --- Calculus unary: recurse into child ---
        case Variation():
            return Variation(pull(expr.expr))

        case TimeDerivative():
            return TimeDerivative(pull(expr.expr))

        case TimeIntegral():
            return TimeIntegral(pull(expr.expr))

        # --- Leaf / unhandled ---
        case _:
            return expr
