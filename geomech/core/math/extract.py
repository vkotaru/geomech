"""Coefficient extraction: extract a query from an expression.

Dispatch grid (query_type × expr_type → result_type):

    query \\ expr  |  Scalar  |  Vector  |  Matrix
    --------------|---------|---------|--------
    Scalar        | Scalar  | Vector  | Matrix
    Vector        | Vector  | Matrix  | (raises)
    Matrix        | (raises)| (raises)| Matrix (partial)

extract_linear_coeff(expr, query) dispatches to the appropriate function.

The coefficient is always a *left* factor: ``c·s``, ``cᵀ·v``, ``C·v`` or ``C·N``. A term
that is linear in the query but has no such factor (a vector query inside a matrix
expression, ``M·Hat(v)`` or ``v wᵀ``, whose coefficient would be a 3-tensor; a scalar
factor ``⟨a, M b⟩`` of a matrix query) raises ``AlgebraicError``; a term where the query
occurs in both children of a product raises ``NonLinearError``. Nothing that contains the
query is ever returned as a zero coefficient.
"""

from __future__ import annotations

from geomech.core.base.expressions import Zero, ZeroMatrix, ZeroVector
from geomech.core.base.types import ExprType
from geomech.core.operations.addition import Add, MAdd, VAdd
from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
from geomech.core.operations.geometry import Cross, Dot, Hat, Transpose, Vee, column_of
from geomech.core.operations.multiplication import (
    MMMul,
    Mul,
    MVMul,
    SMMul,
    SVMul,
    VVMul,
)
from geomech.utils.errors import AlgebraicError, NonLinearError

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _same(expr, query) -> bool:
    """Node identity: same type and, for leaves, same name; otherwise the same rendered tree.

    The type guard is the one ``Expr.__eq__`` uses, so ``Mul(k, x)`` (rendered ``kx``)
    never matches the leaf ``Scalar("kx")``.
    """
    if type(expr) is not type(query):
        return False
    e_name, q_name = getattr(expr, "name", None), getattr(query, "name", None)
    if e_name is not None or q_name is not None:
        return e_name == q_name
    return str(expr) == str(query)


def _has_algebraic(expr, query) -> bool:
    """Check if query appears algebraically in expr (not inside calculus ops).

    Calculus operators (d/dt, δ, ∫) are opaque boundaries — d/dt(η) does
    not algebraically contain η. The pipeline always expands
    d/dt(A+B) → d/dt(A) + d/dt(B) before extraction, so calculus ops
    only wrap leaf nodes by the time this is called.

    Raises ``AlgebraicError`` if the query is found *only inside* a calculus
    operator that is not itself the query (e.g. query ẍ inside d/dt(ẍ), or
    δx inside d/dt(δx)): such a term has no linear coefficient and must
    never be silently dropped.
    """
    if _same(expr, query):
        return True
    if isinstance(expr, (TimeDerivative, TimeIntegral, Variation)):
        if expr.has(query):
            raise AlgebraicError(
                f"extract_linear_coeff: '{query}' appears only inside the calculus"
                f" operator '{expr}' (a higher order derivative, a variation or an"
                " integral of the query); such a term has no linear coefficient"
            )
        return False
    nodes = getattr(expr, "nodes", None)
    if nodes:
        return any(_has_algebraic(n, query) for n in nodes)
    return False


def _zero_or_raise(func, expr, query, zero):
    """Fallback for an unhandled node: *zero* if it is free of *query*, else a typed error.

    Never returns zero for a term that contains the query.
    """
    if _has_algebraic(expr, query):
        raise AlgebraicError(
            f"{func}: unsupported node {type(expr).__name__} contains '{query}' in '{expr}'"
        )
    return zero


def _scalar_from(expr, query):
    """Dispatch a scalar query on *expr*'s type (scalar/vector/matrix expr)."""
    match getattr(expr, "type", None):
        case ExprType.SCALAR:
            return _scalar_from_scalar(expr, query)
        case ExprType.VECTOR:
            return _scalar_from_vector(expr, query)
        case ExprType.MATRIX:
            return _scalar_from_matrix(expr, query)
        case _:
            raise NotImplementedError(f"_scalar_from: expr type | query='{query}', expr='{expr}'")


def _nonlinear(func, node_type, expr, query):
    """Raise NonLinearError for nonlinear (query in both children)."""
    raise NonLinearError(
        f"{func}: nonlinear — '{query}' appears in both sides of {node_type} in '{expr}'"
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def extract_linear_coeff(expr, query):
    """Extract the coefficient of *query* from *expr*.

    Dispatches based on the types of query and expr:
      scalar query + scalar expr → scalar coefficient
      scalar query + vector expr → vector coefficient
      scalar query + matrix expr → matrix coefficient
      vector query + scalar expr → vector coefficient
      vector query + vector expr → matrix coefficient
      vector query + matrix expr → matrix coefficient
    """
    q_type = getattr(query, "type", None)
    e_type = getattr(expr, "type", None)

    match (e_type, q_type):
        case (ExprType.SCALAR, ExprType.SCALAR):
            return _scalar_from_scalar(expr, query)
        case (ExprType.VECTOR, ExprType.SCALAR):
            return _scalar_from_vector(expr, query)
        case (ExprType.MATRIX, ExprType.SCALAR):
            return _scalar_from_matrix(expr, query)
        case (ExprType.SCALAR, ExprType.VECTOR):
            return _vector_from_scalar(expr, query)
        case (ExprType.VECTOR, ExprType.VECTOR):
            return _vector_from_vector(expr, query)
        case (ExprType.MATRIX, ExprType.VECTOR):
            return _vector_from_matrix(expr, query)
        case (ExprType.MATRIX, ExprType.MATRIX):
            return _matrix_from_matrix(expr, query)
        case _:
            raise AlgebraicError(
                f"extract_linear_coeff: a {q_type.name if q_type else q_type} query has no"
                f" left coefficient in a {e_type.name if e_type else e_type} expression"
                f" | query='{query}', expr='{expr}'"
            )


# ---------------------------------------------------------------------------
# Scalar query extraction
# ---------------------------------------------------------------------------


def _scalar_from_scalar(expr, query):
    """Extract scalar query from scalar expr → scalar coefficient.

    expr = coeff * query + ...  →  returns coeff (scalar).
    """
    if _same(expr, query):
        from geomech.core.base.expressions import One

        return One

    match expr:
        case Add(nodes=nodes):
            extracted = [_scalar_from_scalar(n, query) for n in nodes if _has_algebraic(n, query)]
            if not extracted:
                return Zero
            if len(extracted) == 1:
                return extracted[0]
            return Add(*extracted)

        case Mul():
            left, right = expr.left, expr.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear("_scalar_from_scalar", "Mul", expr, query)
            if _has_algebraic(left, query):
                return Mul(_scalar_from_scalar(left, query), right)
            if _has_algebraic(right, query):
                return Mul(left, _scalar_from_scalar(right, query))
            return Zero

        case Dot():
            left, right = expr.left, expr.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear("_scalar_from_scalar", "Dot", expr, query)
            if _has_algebraic(left, query):
                return Dot(_scalar_from_vector(left, query), right)
            if _has_algebraic(right, query):
                return Dot(_scalar_from_vector(right, query), left)
            return Zero

        case VVMul():
            left, right = expr.left, expr.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear("_scalar_from_scalar", "VVMul", expr, query)
            if _has_algebraic(left, query):
                # left is Transpose(v), query inside v → Transpose(coeff) * right
                return VVMul(Transpose(_scalar_from_vector(column_of(left), query)), right)
            if _has_algebraic(right, query):
                # right is plain vector → left * coeff
                return VVMul(left, _scalar_from_vector(right, query))
            return Zero

        case _:
            return _zero_or_raise("_scalar_from_scalar", expr, query, Zero)


def _scalar_from_vector(expr, query):
    """Extract scalar query from vector expr → vector coefficient.

    expr = coeff * query + ...  →  returns coeff (vector).
    """
    match expr:
        case VAdd(nodes=nodes):
            extracted = [_scalar_from_vector(n, query) for n in nodes if _has_algebraic(n, query)]
            if not extracted:
                return ZeroVector
            if len(extracted) == 1:
                return extracted[0]
            return VAdd(*extracted)

        case SVMul():
            vec, scl = expr.left, expr.right
            if _has_algebraic(vec, query) and _has_algebraic(scl, query):
                _nonlinear("_scalar_from_vector", "SVMul", expr, query)
            if _has_algebraic(vec, query):
                return SVMul(_scalar_from_vector(vec, query), scl)
            if _has_algebraic(scl, query):
                if scl == query:
                    return vec
                return SVMul(vec, _scalar_from_scalar(scl, query))
            return ZeroVector

        case MVMul():
            # M v or (row) vᵀ M: both factors are linear, so coefficient = M c or cᵀ M
            mat, vec = expr.left, expr.right
            if _has_algebraic(mat, query) and _has_algebraic(vec, query):
                _nonlinear("_scalar_from_vector", "MVMul", expr, query)
            if _has_algebraic(mat, query):
                return MVMul(_scalar_from(mat, query), vec)
            if _has_algebraic(vec, query):
                return MVMul(mat, _scalar_from(vec, query))
            return ZeroVector

        case Cross():
            left, right = expr.left, expr.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear("_scalar_from_vector", "Cross", expr, query)
            if _has_algebraic(left, query):
                return MVMul(SMMul(Hat(right), -1), _scalar_from_vector(left, query))
            if _has_algebraic(right, query):
                return MVMul(Hat(left), _scalar_from_vector(right, query))
            return ZeroVector

        case Transpose():
            # (c l)ᵀ = l cᵀ: coefficient of l is cᵀ
            return Transpose(_scalar_from_vector(expr.expr, query))

        case Vee():
            # Vee(C l) = l Vee(C): coefficient of l is Vee(C)
            return Vee(_scalar_from_matrix(expr.expr, query))

        case _:
            return _zero_or_raise("_scalar_from_vector", expr, query, ZeroVector)


def _scalar_from_matrix(expr, query):
    """Extract scalar query from matrix expr → matrix coefficient.

    expr = coeff * query + ...  →  returns coeff (matrix).
    """
    match expr:
        case MAdd(nodes=nodes):
            extracted = [_scalar_from_matrix(n, query) for n in nodes if _has_algebraic(n, query)]
            if not extracted:
                return ZeroMatrix
            if len(extracted) == 1:
                return extracted[0]
            return MAdd(*extracted)

        case SMMul():
            mat, scl = expr.left, expr.right
            if _has_algebraic(mat, query) and _has_algebraic(scl, query):
                _nonlinear("_scalar_from_matrix", "SMMul", expr, query)
            if _has_algebraic(mat, query):
                return SMMul(_scalar_from_matrix(mat, query), scl)
            if _has_algebraic(scl, query):
                if scl == query:
                    return mat
                return SMMul(mat, _scalar_from_scalar(scl, query))
            return ZeroMatrix

        case MMMul():
            left, right = expr.left, expr.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear("_scalar_from_matrix", "MMMul", expr, query)
            if _has_algebraic(left, query):
                return MMMul(_scalar_from_matrix(left, query), right)
            if _has_algebraic(right, query):
                return MMMul(left, _scalar_from_matrix(right, query))
            return ZeroMatrix

        case Hat():
            # Hat(c l) = l Hat(c): coefficient of l is Hat(c)
            return Hat(_scalar_from_vector(expr.expr, query))

        case Transpose():
            # (C l)ᵀ = l Cᵀ: coefficient of l is Cᵀ
            return Transpose(_scalar_from_matrix(expr.expr, query))

        case VVMul():
            # Outer product (a l) bᵀ = l (a bᵀ): coefficient of l is a bᵀ (or a cᵀ when l is in b)
            left, right = expr.left, expr.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear("_scalar_from_matrix", "VVMul", expr, query)
            if _has_algebraic(left, query):
                return VVMul(_scalar_from_vector(left, query), right)
            if _has_algebraic(right, query):
                return VVMul(left, Transpose(_scalar_from_vector(column_of(right), query)))
            return ZeroMatrix

        case _:
            return _zero_or_raise("_scalar_from_matrix", expr, query, ZeroMatrix)


# ---------------------------------------------------------------------------
# Vector query extraction
# ---------------------------------------------------------------------------


def _vector_from_scalar(expr, query):
    """Extract vector query from scalar expr → vector coefficient.

    expr = query^T * coeff + ...  →  returns coeff (vector).
    """
    match expr:
        case Add(nodes=nodes):
            extracted = [_vector_from_scalar(n, query) for n in nodes if _has_algebraic(n, query)]
            if not extracted:
                return ZeroVector
            if len(extracted) == 1:
                return extracted[0]
            return VAdd(*extracted)

        case Mul():
            left, right = expr.left, expr.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear("_vector_from_scalar", "Mul", expr, query)
            if _has_algebraic(left, query):
                return SVMul(_vector_from_scalar(left, query), right)
            if _has_algebraic(right, query):
                return SVMul(_vector_from_scalar(right, query), left)
            return ZeroVector

        case Dot():
            left, right = expr.left, expr.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear("_vector_from_scalar", "Dot", expr, query)
            if _has_algebraic(left, query):
                if left == query:
                    return right
                if left.type == ExprType.VECTOR:
                    return MVMul(Transpose(_vector_from_vector(left, query)), right)
                raise NotImplementedError(
                    f"_vector_from_scalar: Dot.left type={left.type}"
                    f" | query='{query}', expr='{expr}'"
                )
            if _has_algebraic(right, query):
                if right == query:
                    return left
                if right.type == ExprType.VECTOR:
                    return MVMul(Transpose(_vector_from_vector(right, query)), left)
                raise NotImplementedError(
                    f"_vector_from_scalar: Dot.right type={right.type}"
                    f" | query='{query}', expr='{expr}'"
                )
            return ZeroVector

        case VVMul():
            # VVMul(Transpose(v), w) = v^T * w — scalar inner product
            left, right = expr.left, expr.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear("_vector_from_scalar", "VVMul", expr, query)
            if _has_algebraic(left, query):
                # left = Transpose(v), query in v → coeff^T * right
                inner = column_of(left)
                if inner == query:
                    return right
                return MVMul(Transpose(_vector_from_vector(inner, query)), right)
            if _has_algebraic(right, query):
                if right == query:
                    return column_of(left)
                return MVMul(Transpose(_vector_from_vector(right, query)), column_of(left))
            return ZeroVector

        case _:
            return _zero_or_raise("_vector_from_scalar", expr, query, ZeroVector)


def _vector_from_vector(expr, query):
    """Extract vector query from vector expr → matrix coefficient.

    expr = coeff * query + ...  →  returns coeff (matrix).
    """
    from geomech.core.base.expressions import IdentityMatrix

    if _same(expr, query):
        return IdentityMatrix

    match expr:
        case VAdd(nodes=nodes):
            extracted = [_vector_from_vector(n, query) for n in nodes if _has_algebraic(n, query)]
            if not extracted:
                return ZeroMatrix
            if len(extracted) == 1:
                return extracted[0]
            return MAdd(*extracted)

        case Cross():
            left, right = expr.left, expr.right
            if left == right:
                return ZeroMatrix  # x × x = 0 identically
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                # before the shortcuts: Cross(q, Cross(a, q)) is quadratic in q
                _nonlinear("_vector_from_vector", "Cross", expr, query)
            if right == query:
                return Hat(left)
            if left == query:
                return SMMul(Hat(right), -1)
            if _has_algebraic(left, query):
                return MMMul(SMMul(Hat(right), -1), _vector_from_vector(left, query))
            if _has_algebraic(right, query):
                return MMMul(Hat(left), _vector_from_vector(right, query))
            return ZeroMatrix

        case MVMul():
            mat, vec = expr.left, expr.right
            if _has_algebraic(mat, query) and _has_algebraic(vec, query):
                # MVMul(Hat(q), q), MVMul(q aᵀ, q), ... are quadratic in q
                _nonlinear("_vector_from_vector", "MVMul", expr, query)
            if _has_algebraic(vec, query):
                if vec == query:
                    return mat
                return MMMul(mat, _vector_from_vector(vec, query))
            if _has_algebraic(mat, query):
                return _vector_from_mvmul_mat(mat, vec, query)
            return ZeroMatrix

        case SVMul():
            vec, scl = expr.left, expr.right
            if _has_algebraic(vec, query) and _has_algebraic(scl, query):
                _nonlinear("_vector_from_vector", "SVMul", expr, query)
            if _has_algebraic(vec, query):
                return SMMul(_vector_from_vector(vec, query), scl)
            if _has_algebraic(scl, query):
                return _vector_from_vector_svmul_scalar(vec, scl, query)
            return ZeroMatrix

        case _:
            return _zero_or_raise("_vector_from_vector", expr, query, ZeroMatrix)


def _vector_from_mvmul_mat(mat, vec, query):
    """Handle _vector_from_vector(MVMul(mat, vec), query) where *mat* contains query."""
    match mat:
        case MAdd(nodes=nodes):
            extracted = [
                _vector_from_vector(MVMul(n, vec), query)
                for n in nodes
                if _has_algebraic(n, query)
            ]
            if not extracted:
                return ZeroMatrix
            if len(extracted) == 1:
                return extracted[0]
            return MAdd(*extracted)

        case MMMul():
            left, right = mat.left, mat.right
            if _has_algebraic(left, query):
                return _vector_from_vector(MVMul(left, MVMul(right, vec)), query)
            if _has_algebraic(right, query):
                return MMMul(left, _vector_from_vector(MVMul(right, vec), query))
            return ZeroMatrix

        case SMMul():
            inner_mat, scl = mat.left, mat.right
            if _has_algebraic(inner_mat, query):
                return _vector_from_vector(MVMul(inner_mat, SVMul(vec, scl)), query)
            if _has_algebraic(scl, query):
                # SMMul(M, scl) * vec = M * (scl * vec) → query in SVMul scalar
                return _vector_from_vector(MVMul(inner_mat, SVMul(vec, scl)), query)
            return ZeroMatrix

        case Hat():
            if mat.expr == query:
                return SMMul(Hat(vec), -1)
            if _has_algebraic(mat.expr, query):
                return MMMul(SMMul(Hat(vec), -1), _vector_from_vector(mat.expr, query))
            return ZeroMatrix

        case Transpose():
            inner = mat.expr
            match inner:
                case Hat():
                    # Tr(Hat(y)) * vec = -Hat(y) * vec → rewrite
                    return _vector_from_vector(MVMul(Hat(inner.expr), SVMul(vec, -1)), query)
                case SMMul():
                    # Tr(SMMul(M, s)) * vec = M^T * (s * vec)
                    return _vector_from_vector(
                        MVMul(Transpose(inner.left), SVMul(vec, inner.right)),
                        query,
                    )
                case MMMul():
                    # Tr(MMMul(A, B)) * vec = B^T * A^T * vec
                    left, right = inner.left, inner.right
                    return _vector_from_vector(
                        MVMul(Transpose(right), MVMul(Transpose(left), vec)),
                        query,
                    )
                case _:
                    return _zero_or_raise("_vector_from_mvmul_mat", mat, query, ZeroMatrix)

        case VVMul():
            # VVMul(v, Transpose(w)) is outer product (matrix).
            # (v * w^T) * vec = v * (w^T * vec) = SVMul(v, Dot(w, vec))
            outer_left, outer_right = mat.left, mat.right
            return _vector_from_vector(SVMul(outer_left, Dot(column_of(outer_right), vec)), query)

        case _:
            return _zero_or_raise("_vector_from_mvmul_mat", mat, query, ZeroMatrix)


def _vector_from_vector_svmul_scalar(vec, scl, query):
    """Handle _vector_from_vector(SVMul(vec, scl), query) where *scl* contains query.

    Redistributes the scalar structure so extraction can proceed:
      SVMul(vec, Add(a, b))  → extract from SVMul(vec, a) + SVMul(vec, b)
      SVMul(vec, Mul(a, b))  → extract from SVMul(SVMul(vec, other), containing)
      SVMul(vec, Dot(a, b))  → extract from MVMul(VVMul(vec, Transpose(other)), containing)
    """
    full_expr = SVMul(vec, scl)
    match scl:
        case Add(nodes=nodes):
            extracted = [
                _vector_from_vector(SVMul(vec, n), query)
                for n in nodes
                if _has_algebraic(n, query)
            ]
            if not extracted:
                return ZeroMatrix
            if len(extracted) == 1:
                return extracted[0]
            return MAdd(*extracted)

        case Mul():
            left, right = scl.left, scl.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear("_vector_from_vector", "Mul(SVMul scalar)", full_expr, query)
            if _has_algebraic(left, query):
                return _vector_from_vector(SVMul(SVMul(vec, right), left), query)
            if _has_algebraic(right, query):
                return _vector_from_vector(SVMul(SVMul(vec, left), right), query)
            return ZeroMatrix

        case Dot():
            left, right = scl.left, scl.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear("_vector_from_vector", "Dot(SVMul scalar)", full_expr, query)
            if _has_algebraic(left, query):
                return _vector_from_vector(MVMul(VVMul(vec, Transpose(right)), left), query)
            if _has_algebraic(right, query):
                return _vector_from_vector(MVMul(VVMul(vec, Transpose(left)), right), query)
            return ZeroMatrix

        case VVMul():
            # VVMul(Transpose(a), b) = a^T * b — same structure as Dot
            left, right = scl.left, scl.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear(
                    "_vector_from_vector",
                    "VVMul(SVMul scalar)",
                    full_expr,
                    query,
                )
            # left = Transpose(a), query in a
            if _has_algebraic(left, query):
                return _vector_from_vector(
                    MVMul(VVMul(vec, Transpose(right)), column_of(left)), query
                )
            # right = b, query in b
            if _has_algebraic(right, query):
                return _vector_from_vector(MVMul(VVMul(vec, left), right), query)
            return ZeroMatrix

        case _:
            return _zero_or_raise("_vector_from_vector", full_expr, query, ZeroMatrix)


def _no_left_coefficient(func, expr, query):
    raise AlgebraicError(
        f"{func}: '{expr}' is linear in '{query}' but has no left matrix coefficient"
        " (the map is a 3-tensor, not C * query)"
    )


def _vector_from_matrix(expr, query):
    """Extract vector query from matrix expr → matrix coefficient.

    A matrix that depends linearly on a vector (``Hat(v)``, ``v wᵀ``, products and sums
    of these) has no 3×3 left coefficient ``C`` with ``C v = expr``; every such term raises
    ``AlgebraicError``. Only terms free of the query reduce, to ``ZeroMatrix``.
    """
    match expr:
        case MAdd(nodes=nodes):
            extracted = [_vector_from_matrix(n, query) for n in nodes if _has_algebraic(n, query)]
            if not extracted:
                return ZeroMatrix
            if len(extracted) == 1:
                return extracted[0]
            return MAdd(*extracted)

        case MMMul():
            left, right = expr.left, expr.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear("_vector_from_matrix", "MMMul", expr, query)
            if _has_algebraic(left, query):
                if left == query:
                    return right
                return MMMul(_vector_from_matrix(left, query), right)
            if _has_algebraic(right, query):
                if right == query:
                    return left
                return MMMul(left, _vector_from_matrix(right, query))
            return ZeroMatrix

        case VVMul():
            if _has_algebraic(expr, query):
                raise AlgebraicError(
                    f"Cannot extract vector query '{query}' from outer"
                    f" product '{expr}': rank-1 matrix does not"
                    " factor as C * q"
                )
            return ZeroMatrix

        case SMMul():
            mat, scl = expr.left, expr.right
            if _has_algebraic(mat, query) and _has_algebraic(scl, query):
                _nonlinear("_vector_from_matrix", "SMMul", expr, query)
            if _has_algebraic(mat, query):
                return SMMul(_vector_from_matrix(mat, query), scl)
            if _has_algebraic(scl, query):
                raise AlgebraicError(
                    f"Cannot extract vector query '{query}' from"
                    f" '{expr}': scalar factor depends on query"
                    " (q^T·w), does not factor as C * q"
                )
            return ZeroMatrix

        case Hat() | Transpose():
            if _has_algebraic(expr, query):
                _no_left_coefficient("_vector_from_matrix", expr, query)
            return ZeroMatrix

        case _:
            return _zero_or_raise("_vector_from_matrix", expr, query, ZeroMatrix)


# ---------------------------------------------------------------------------
# Matrix query extraction (partial)
# ---------------------------------------------------------------------------


def _matrix_from_matrix(expr, query):
    """Extract matrix query from matrix expr → matrix coefficient."""
    from geomech.core.base.expressions import IdentityMatrix

    if _same(expr, query):
        return IdentityMatrix

    match expr:
        case MAdd(nodes=nodes):
            extracted = [_matrix_from_matrix(n, query) for n in nodes if _has_algebraic(n, query)]
            if not extracted:
                return ZeroMatrix
            if len(extracted) == 1:
                return extracted[0]
            return MAdd(*extracted)

        case MMMul():
            left, right = expr.left, expr.right
            if _has_algebraic(left, query) and _has_algebraic(right, query):
                _nonlinear("_matrix_from_matrix", "MMMul", expr, query)
            if _has_algebraic(right, query):
                if right == query:
                    return left
                return MMMul(left, _matrix_from_matrix(right, query))
            if _has_algebraic(left, query):
                if left == query:
                    return right
                return MMMul(_matrix_from_matrix(left, query), right)
            return ZeroMatrix

        case SMMul():
            mat, scl = expr.left, expr.right
            if _has_algebraic(mat, query) and _has_algebraic(scl, query):
                _nonlinear("_matrix_from_matrix", "SMMul", expr, query)
            if _has_algebraic(mat, query):
                return SMMul(_matrix_from_matrix(mat, query), scl)
            if _has_algebraic(scl, query):
                # Z * <a, M b> is linear in M but its coefficient is a 4-tensor
                _no_left_coefficient("_matrix_from_matrix", expr, query)
            return ZeroMatrix

        case _:
            return _zero_or_raise("_matrix_from_matrix", expr, query, ZeroMatrix)
