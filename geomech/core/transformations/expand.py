"""Distribute multiplication over addition (expand products).

Bottom-up traversal: expand children first, then distribute at the
current node.  Single pass — no has_nested_add predicate needed.

Distribution rules implemented:
  Mul  over Add   → Add     (scalar)
  Dot  over VAdd  → Add     (scalar from vectors)
  Cross over VAdd → VAdd    (vector)
  SVMul over VAdd/Add → VAdd (vector, left=vector right=scalar)
  MVMul over MAdd/VAdd → VAdd (vector)
  SMMul over MAdd/Add → MAdd (matrix, left=matrix right=scalar)
  MMMul over MAdd → MAdd    (matrix)
  VVMul over VAdd → MAdd    (outer product  v wᵀ; a VAdd in v or inside Transpose(w))
  VVMul over VAdd → Add     (inner product  vᵀ w; a VAdd inside Transpose(v) or in w)
  Hat/Vee — recurse into child

``Transpose`` is not expanded on its own (Transpose(MAdd) and Transpose(VVMul)
are left intact); the VVMul rule is the only place that looks through a
``Transpose(VAdd)`` factor, because a bare ``VAdd`` of transposed vectors is not
a valid ``VVMul`` operand.
"""

from geomech.core.base.types import ExprType
from geomech.core.operations.addition import Add, MAdd, VAdd
from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
from geomech.core.operations.geometry import Cross, Dot, Hat, Transpose, Vee
from geomech.core.operations.multiplication import MMMul, Mul, MVMul, SMMul, SVMul, VVMul


def _distribute(l, r, op_cls, sum_cls, l_sum_cls, r_sum_cls):
    """Distribute *op_cls* over addition on either or both sides.

    l, r        — already-expanded left and right children
    op_cls      — constructor for the binary operation (Mul, Dot, …)
    sum_cls     — constructor for the result addition (Add, VAdd, MAdd)
    l_sum_cls   — addition type to check on the left  (VAdd, MAdd, Add)
    r_sum_cls   — addition type to check on the right (VAdd, MAdd, Add)

    If both sides are sums  → cartesian product of terms.
    If one side is a sum    → distribute over that side.
    Otherwise               → reconstruct with expanded children.
    """
    l_is_sum = isinstance(l, l_sum_cls)
    r_is_sum = isinstance(r, r_sum_cls)

    if l_is_sum and r_is_sum:
        return sum_cls(*[op_cls(nl, nr) for nl in l.nodes for nr in r.nodes])
    if l_is_sum:
        return sum_cls(*[op_cls(nl, r) for nl in l.nodes])
    if r_is_sum:
        return sum_cls(*[op_cls(l, nr) for nr in r.nodes])
    return op_cls(l, r)


def _transposed_addends(node):
    """Addends of a ``Transpose(VAdd)`` factor as transposed vectors, else ``None``.

    The inner expression is expanded first so that ``Transpose(Cross(a, b + c))``
    is also split; a ``Transpose`` whose inner does not expand to a ``VAdd`` is
    returned untouched by the caller (the original node, not the expanded one,
    so that no other rule's input changes).
    """
    if not isinstance(node, Transpose):
        return None
    inner = expand(node.expr)
    if isinstance(inner, VAdd):
        return [Transpose(n) for n in inner.nodes]
    return None


def _expand_vvmul(expr):
    """Distribute an outer product ``v wᵀ`` (→ MAdd) or an inner product ``vᵀ w``
    (→ Add) over a ``VAdd`` in either factor.

    Each addend of a ``Transpose(VAdd)`` factor becomes ``Transpose(addend)`` so
    every produced ``VVMul`` has the same (plain, Transpose) shape as the input.
    The plain factor is expanded as usual; the ``Transpose`` factor is only split
    (never replaced by its expanded inner) so that a non-sum inner is untouched.
    """
    if expr.type == ExprType.MATRIX:
        # outer product: left is a plain vector, right is Transpose(w)
        l = expand(expr.left)
        ls = list(l.nodes) if isinstance(l, VAdd) else [l]
        rs = _transposed_addends(expr.right) or [expr.right]
        if len(ls) == 1 and len(rs) == 1:
            return VVMul(l, expr.right)
        return MAdd(*[VVMul(nl, nr) for nl in ls for nr in rs])
    # inner product: left is Transpose(v), right is a plain vector
    r = expand(expr.right)
    ls = _transposed_addends(expr.left) or [expr.left]
    rs = list(r.nodes) if isinstance(r, VAdd) else [r]
    if len(ls) == 1 and len(rs) == 1:
        return VVMul(expr.left, r)
    return Add(*[VVMul(nl, nr) for nl in ls for nr in rs])


def expand(expr):
    match expr:
        # --- N-ary additions: expand each term ---
        case Add(nodes=nodes):
            return Add(*[expand(n) for n in nodes])

        case VAdd(nodes=nodes):
            return VAdd(*[expand(n) for n in nodes])

        case MAdd(nodes=nodes):
            return MAdd(*[expand(n) for n in nodes])

        # --- Scalar * Scalar ---
        case Mul():
            l, r = expand(expr.left), expand(expr.right)
            return _distribute(l, r, Mul, Add, Add, Add)

        # --- Dot product (vector x vector → scalar) ---
        case Dot():
            l, r = expand(expr.left), expand(expr.right)
            return _distribute(l, r, Dot, Add, VAdd, VAdd)

        # --- Cross product (vector x vector → vector) ---
        case Cross():
            l, r = expand(expr.left), expand(expr.right)
            return _distribute(l, r, Cross, VAdd, VAdd, VAdd)

        # --- Scalar * Vector  (left=vector, right=scalar) ---
        case SVMul():
            l, r = expand(expr.left), expand(expr.right)
            return _distribute(l, r, SVMul, VAdd, VAdd, Add)

        # --- Matrix * Vector ---
        case MVMul():
            l, r = expand(expr.left), expand(expr.right)
            return _distribute(l, r, MVMul, VAdd, MAdd, VAdd)

        # --- Scalar * Matrix  (left=matrix, right=scalar) ---
        case SMMul():
            l, r = expand(expr.left), expand(expr.right)
            return _distribute(l, r, SMMul, MAdd, MAdd, Add)

        # --- Matrix * Matrix ---
        case MMMul():
            l, r = expand(expr.left), expand(expr.right)
            return _distribute(l, r, MMMul, MAdd, MAdd, MAdd)

        # --- Vector * Vector (outer product → matrix, inner product → scalar) ---
        case VVMul():
            return _expand_vvmul(expr)

        # --- Unary ops: recurse into child ---
        # Hat and Vee are linear: Hat(a + b) = Hat(a) + Hat(b)
        case Hat():
            inner = expand(expr.expr)
            if isinstance(inner, VAdd):
                return MAdd(*[Hat(n) for n in inner.nodes])
            return Hat(inner)

        case Vee():
            inner = expand(expr.expr)
            if isinstance(inner, VAdd):
                return VAdd(*[Vee(n) for n in inner.nodes])
            return Vee(inner)

        case Variation():
            return Variation(expand(expr.expr))

        case TimeDerivative():
            return TimeDerivative(expand(expr.expr))

        case TimeIntegral():
            return TimeIntegral(expand(expr.expr))

        # --- Leaves and everything else ---
        case _:
            return expr
