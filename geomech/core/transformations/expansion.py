"""Distribute multiplication over addition (expand products).

Bottom-up traversal: expand children first, then distribute at the
current node.  Single pass.

Distribution rules (name, identity; none needs a hypothesis: they are bilinearity of
each product and linearity of Hat). A rule fires only when an operand is a sum.

  expand.mul_over_add          Mul over Add → Add (scalar)
  expand.dot_over_vadd         Dot over VAdd → Add (scalar from vectors)
  expand.cross_over_vadd       Cross over VAdd → VAdd (vector)
  expand.svmul_over_sum        SVMul over VAdd/Add → VAdd (vector, left=vector right=scalar)
  expand.row_svmul_over_sum    a scaled row rᵀ s over a sum of rows or an Add scalar
  expand.mvmul_over_sum        MVMul over MAdd/VAdd → VAdd (vector)
  expand.row_mvmul_over_sum    a row times a matrix, rᵀ M, over a sum of rows on the left
                               or an MAdd on the right
  expand.smmul_over_sum        SMMul over MAdd/Add → MAdd (matrix, left=matrix right=scalar)
  expand.mmmul_over_madd       MMMul over MAdd → MAdd (matrix)
  expand.vvmul_outer_over_vadd VVMul over VAdd → MAdd (outer product v wᵀ; a VAdd in v or
                               inside Transpose(w))
  expand.vvmul_inner_over_vadd VVMul over VAdd → Add (inner product vᵀ w; a VAdd inside
                               Transpose(v) or in w)
  expand.hat_over_vadd         Hat(a + b) = Hat(a) + Hat(b)
  expand.vee_over_madd         Vee(A + B) = Vee(A) + Vee(B)

``Transpose`` is not expanded on its own (Transpose(MAdd) and Transpose(VVMul)
are left intact); the VVMul rule is the only place that looks through a
``Transpose(VAdd)`` factor, because a bare ``VAdd`` of transposed vectors is not
a valid ``VVMul`` operand.

Row vectors that are not a ``Transpose`` node (``aᵀK``, ``s aᵀ``, ``aᵀ + bᵀ``) are
distributed in place and keep their row form: ``(aᵀ(K + J)) b`` expands to
``aᵀK b + aᵀJ b``. Inside such a row, a ``Transpose(VAdd)`` factor is split the same way
as in a ``VVMul``.
"""

from geomech.core.base.types import ExprType
from geomech.core.operations.addition import Add, MAdd, VAdd
from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
from geomech.core.operations.geometry import Cross, Dot, Hat, Transpose, Vee, is_row_vector
from geomech.core.operations.multiplication import MMMul, Mul, MVMul, SMMul, SVMul, VVMul
from geomech.core.transformations.utils import calculus_node
from geomech.utils.trace import rule


def _distribute(l, r, op_cls, sum_cls, l_sum_cls, r_sum_cls, name: str = "", expr=None):
    """Distribute *op_cls* over addition on either or both sides.

    l, r        — already-expanded left and right children
    op_cls      — constructor for the binary operation (Mul, Dot, …)
    sum_cls     — constructor for the result addition (Add, VAdd, MAdd)
    l_sum_cls   — addition type to check on the left  (VAdd, MAdd, Add)
    r_sum_cls   — addition type to check on the right (VAdd, MAdd, Add)
    name, expr  — the rule name and the original node, for ``rule`` (tracing only)

    If both sides are sums  → cartesian product of terms.
    If one side is a sum    → distribute over that side.
    Otherwise               → reconstruct with expanded children.
    """
    l_is_sum = isinstance(l, l_sum_cls)
    r_is_sum = isinstance(r, r_sum_cls)

    if l_is_sum and r_is_sum:
        return rule(name, expr, sum_cls(*[op_cls(nl, nr) for nl in l.nodes for nr in r.nodes]))
    if l_is_sum:
        return rule(name, expr, sum_cls(*[op_cls(nl, r) for nl in l.nodes]))
    if r_is_sum:
        return rule(name, expr, sum_cls(*[op_cls(l, nr) for nr in r.nodes]))
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


def _row_addends(row):
    """The addends of the expanded row vector *row*, each a row vector.

    A ``Transpose`` node is only split (``_transposed_addends``), as in ``VVMul``; any
    other row shape is expanded and a resulting ``VAdd`` of rows is split.
    """
    if isinstance(row, Transpose):
        return _transposed_addends(row) or [row]
    out = expand(row)
    return list(out.nodes) if isinstance(out, VAdd) else [out]


def _expand_row_mvmul(expr):
    """Distribute the row ``rᵀ M`` over a sum of rows and over an ``MAdd`` matrix."""
    ls = _row_addends(expr.left)
    r = expand(expr.right)
    rs = list(r.nodes) if isinstance(r, MAdd) else [r]
    if len(ls) == 1 and len(rs) == 1:
        return MVMul(ls[0], rs[0])
    return rule(
        "expand.row_mvmul_over_sum", expr, VAdd(*[MVMul(nl, nr) for nl in ls for nr in rs])
    )


def _expand_row_svmul(expr):
    """Distribute the row ``rᵀ s`` over a sum of rows and over an ``Add`` scalar."""
    ls = _row_addends(expr.left)
    r = expand(expr.right)
    rs = list(r.nodes) if isinstance(r, Add) else [r]
    if len(ls) == 1 and len(rs) == 1:
        return SVMul(ls[0], rs[0])
    return rule(
        "expand.row_svmul_over_sum", expr, VAdd(*[SVMul(nl, nr) for nl in ls for nr in rs])
    )


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
        rs = _row_addends(expr.right)
        if len(ls) == 1 and len(rs) == 1:
            return VVMul(l, rs[0])
        return rule(
            "expand.vvmul_outer_over_vadd", expr, MAdd(*[VVMul(nl, nr) for nl in ls for nr in rs])
        )
    # inner product: left is Transpose(v), right is a plain vector
    r = expand(expr.right)
    ls = _row_addends(expr.left)
    rs = list(r.nodes) if isinstance(r, VAdd) else [r]
    if len(ls) == 1 and len(rs) == 1:
        return VVMul(ls[0], r)
    return rule(
        "expand.vvmul_inner_over_vadd", expr, Add(*[VVMul(nl, nr) for nl in ls for nr in rs])
    )


def expand(expr):
    """Distribute every product over sums, bottom up (rules in the module docstring)."""
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
            return _distribute(l, r, Mul, Add, Add, Add, "expand.mul_over_add", expr)

        # --- Dot product (vector x vector → scalar) ---
        case Dot():
            l, r = expand(expr.left), expand(expr.right)
            return _distribute(l, r, Dot, Add, VAdd, VAdd, "expand.dot_over_vadd", expr)

        # --- Cross product (vector x vector → vector) ---
        case Cross():
            l, r = expand(expr.left), expand(expr.right)
            return _distribute(l, r, Cross, VAdd, VAdd, VAdd, "expand.cross_over_vadd", expr)

        # --- Scalar * Vector  (left=vector, right=scalar) ---
        case SVMul() if is_row_vector(expr.left):
            return _expand_row_svmul(expr)

        case SVMul():
            l, r = expand(expr.left), expand(expr.right)
            return _distribute(l, r, SVMul, VAdd, VAdd, Add, "expand.svmul_over_sum", expr)

        # --- Matrix * Vector ---
        case MVMul() if expr.left.type == ExprType.VECTOR:
            return _expand_row_mvmul(expr)

        case MVMul():
            l, r = expand(expr.left), expand(expr.right)
            return _distribute(l, r, MVMul, VAdd, MAdd, VAdd, "expand.mvmul_over_sum", expr)

        # --- Scalar * Matrix  (left=matrix, right=scalar) ---
        case SMMul():
            l, r = expand(expr.left), expand(expr.right)
            return _distribute(l, r, SMMul, MAdd, MAdd, Add, "expand.smmul_over_sum", expr)

        # --- Matrix * Matrix ---
        case MMMul():
            l, r = expand(expr.left), expand(expr.right)
            return _distribute(l, r, MMMul, MAdd, MAdd, MAdd, "expand.mmmul_over_madd", expr)

        # --- Vector * Vector (outer product → matrix, inner product → scalar) ---
        case VVMul():
            return _expand_vvmul(expr)

        # --- Unary ops: recurse into child ---
        # Hat is linear: Hat(a + b) = Hat(a) + Hat(b)
        case Hat():
            inner = expand(expr.expr)
            if isinstance(inner, VAdd):
                return rule("expand.hat_over_vadd", expr, MAdd(*[Hat(n) for n in inner.nodes]))
            return Hat(inner)

        # Vee is linear: Vee(A + B) = Vee(A) + Vee(B)
        case Vee():
            inner = expand(expr.expr)
            if isinstance(inner, MAdd):
                return rule("expand.vee_over_madd", expr, VAdd(*[Vee(n) for n in inner.nodes]))
            return Vee(inner)

        case Variation():
            return calculus_node(Variation, expand(expr.expr))

        case TimeDerivative():
            return calculus_node(TimeDerivative, expand(expr.expr))

        case TimeIntegral():
            return TimeIntegral(expand(expr.expr))

        # --- Leaves and everything else ---
        case _:
            return expr
