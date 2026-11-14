"""Simplification passes.

simplify(expr)      — one pass: ``pull`` (hoist scalars), ``vector_rules`` (geometric
                      identities), then ``_eliminate`` (bottom-up zeros, identities and
                      numeric folding, the rules below).
full_simplify(expr) — ``expand`` once, then ``simplify`` until the printed tree stops
                      changing (at most ``max_iter`` passes; hitting the cap warns).

Rules of ``_eliminate`` (name, identity, hypothesis the code checks). "Zero" means the
operand is flagged zero (``Zero``, ``ZeroVector``, ``ZeroMatrix``, a ``Constant``
symbol's variation...), "one" a numeric leaf with value 1, "identity" a leaf flagged
``Identity``.

  simplify.add_fold                Add: zeros dropped, numeric leaves summed, like terms
                                   merged; float sums within 1e-12 of the largest
                                   magnitude are cancelled to 0 (``_CANCEL_RTOL``)
  simplify.vadd_fold               VAdd: zeros dropped, like terms merged
  simplify.madd_fold               MAdd: zeros dropped, like terms merged
  simplify.mul_fold                Mul: a zero factor gives 0, a factor one is dropped,
                                   numeric factors of a Mul chain are multiplied together
  simplify.svmul_zero              v·0 = 0·s = 0 (vector)
  simplify.svmul_one               v·1 = v
  simplify.smmul_zero              M·0 = 0·s = 0 (matrix)
  simplify.smmul_one               M·1 = M
  simplify.mvmul_zero              0·v = M·0 = 0
  simplify.mvmul_identity          I·v = v                    left factor flagged identity
  simplify.mvmul_hat_to_cross      Hat(a)·v = Cross(a, v)
  simplify.mvmul_reassociate       (A·B)·v = A·(B·v)
  simplify.mvmul_combine           A·(B·v) = (A·B)·v          only when A·B simplifies to a
                                                              non-product (e.g. RᵀR = I)
  simplify.mmmul_zero              0·B = A·0 = 0
  simplify.mmmul_identity          I·B = B, A·I = A           factor flagged identity
  simplify.orthogonal_transpose    RᵀR = I and RRᵀ = I        R flagged orthogonal (SO3),
                                                              both factors equal by name
  simplify.vvmul_zero              0ᵀw = vᵀ0 = 0 (scalar), v0ᵀ = 0 (matrix)
  simplify.dot_zero                Dot(0, w) = Dot(v, 0) = 0
  simplify.cross_zero              Cross(0, w) = Cross(v, 0) = 0
  simplify.hat_zero                Hat(0) = 0
  simplify.hat_vee_skew            Hat(Vee(X)) = X            X skew-symmetric: a Hat, a
                                                              leaf flagged
                                                              SkewSymmetricMatrix, or a
                                                              sum/scaling/transpose of such
  simplify.vee_zero                Vee(0) = 0
  simplify.vee_hat                 Vee(Hat(v)) = v
  simplify.transpose_transpose     (Xᵀ)ᵀ = X
  simplify.transpose_symmetric     Jᵀ = J                     J flagged symmetric
  simplify.transpose_hat_skew      Hat(v)ᵀ = -Hat(v)
  simplify.transpose_scaled        (M s)ᵀ = Mᵀ s
  simplify.transpose_scaled_row    (aᵀ s)ᵀ = a s
  simplify.transpose_product       (A·B)ᵀ = Bᵀ·Aᵀ
  simplify.variation_zero          δ(0) = 0
  simplify.time_derivative_zero    d/dt(0) = 0
  simplify.time_derivative_integral   d/dt(∫x dt) = x
  simplify.time_integral_zero      ∫0 dt = 0
  simplify.time_integral_derivative   ∫d/dt(x) dt = x

The rules of the other two stages of a pass are listed in ``pull`` and ``vector_rules``.
"""

from __future__ import annotations

import warnings

from geomech.core.base.expressions import (
    IdentityMatrix,
    Number,
    Scalar,
    Zero,
    ZeroMatrix,
    ZeroVector,
)
from geomech.core.base.types import ExprType
from geomech.core.operations.addition import Add, MAdd, VAdd
from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
from geomech.core.operations.geometry import Cross, Dot, Hat, Transpose, Vee
from geomech.core.operations.multiplication import (
    MMMul,
    Mul,
    MVMul,
    SMMul,
    SVMul,
    VVMul,
)
from geomech.utils import trace as _trace
from geomech.utils.trace import rule

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def simplify(expr):
    """Pull scalars, apply vector identities, then eliminate zeros and
    combine numeric constants — single call does everything."""
    from geomech.core.transformations.pulling import pull
    from geomech.core.transformations.vector_identities import vector_rules

    expr = pull(expr)
    expr = vector_rules(expr)
    return _eliminate(expr)


def full_simplify(expr, max_iter=10):
    """Expand once, then loop simplify until convergence.

    Parameters
    ----------
    max_iter : int
        Safety cap on the number of passes (default 10). Hitting it emits a
        ``RuntimeWarning`` ("no confirmed fixed point after N passes") naming the
        expression: the result may be unsimplified. A pass that changes nothing ends the
        loop, so the warning fires only when the last allowed pass still changed the
        expression; that last result may already be the fixed point, which only one more
        pass could confirm (hence "confirmed").
    """
    from geomech.core.transformations.expansion import expand

    expr = expand(expr)
    for _ in range(max_iter):
        prev = expr
        expr = simplify(expr)
        if str(expr) == str(prev):
            break
    else:
        why = (
            "max_iter=0: no simplify pass ran"
            if max_iter <= 0
            else f"the last pass still changed the expression, so no fixed point after"
            f" {max_iter} passes was confirmed by an unchanged pass"
        )
        warnings.warn(
            f"full_simplify: no confirmed fixed point after {max_iter} passes ({why});"
            f" result may be unsimplified: {expr!r}",
            RuntimeWarning,
            stacklevel=2,
        )
    return expr


# ---------------------------------------------------------------------------
# Bottom-up zero / identity / constant elimination
# ---------------------------------------------------------------------------


def _folded_from(cls, children, expr):
    """The node a fold rule starts from: *expr* rebuilt from its eliminated *children*.

    The ``*_fold`` rules log ``before ==> after`` against this node, so a fold is logged
    only when the fold itself changed something, not when a child rule did (the child
    logs its own line). Called only while rule tracing is on.
    """
    try:
        return cls(*children)
    except Exception:  # noqa: BLE001 - a debugging aid must not change the result
        return expr


def _eliminate(expr):
    """Bottom-up pass: zeros, identities, numeric constant combining."""
    match expr:
        # ---- n-ary addition ----
        case Add(nodes=nodes):
            children = [_eliminate(n) for n in nodes]
            out = _simplify_add(children)
            if _trace._TRACE_RULES:
                rule("simplify.add_fold", _folded_from(Add, children, expr), out)
            return out
        case VAdd(nodes=nodes):
            children = [_eliminate(n) for n in nodes]
            out = _simplify_vadd(children)
            if _trace._TRACE_RULES:
                rule("simplify.vadd_fold", _folded_from(VAdd, children, expr), out)
            return out
        case MAdd(nodes=nodes):
            children = [_eliminate(n) for n in nodes]
            out = _simplify_madd(children)
            if _trace._TRACE_RULES:
                rule("simplify.madd_fold", _folded_from(MAdd, children, expr), out)
            return out

        # ---- scalar * scalar ----
        case Mul():
            children = [_eliminate(expr.left), _eliminate(expr.right)]
            out = _simplify_mul(*children)
            if _trace._TRACE_RULES:
                rule("simplify.mul_fold", _folded_from(Mul, children, expr), out)
            return out

        # ---- scalar * vector ----
        case SVMul():
            vec, scl = _eliminate(expr.left), _eliminate(expr.right)
            if _is_zero(scl) or _is_zero(vec):
                return rule("simplify.svmul_zero", expr, ZeroVector)
            if _is_one(scl):
                return rule("simplify.svmul_one", expr, vec)
            return SVMul(vec, scl)

        # ---- scalar * matrix ----
        case SMMul():
            mat, scl = _eliminate(expr.left), _eliminate(expr.right)
            if _is_zero(scl) or _is_zero(mat):
                return rule("simplify.smmul_zero", expr, ZeroMatrix)
            if _is_one(scl):
                return rule("simplify.smmul_one", expr, mat)
            return SMMul(mat, scl)

        # ---- matrix * vector ----
        case MVMul():
            mat, vec = _eliminate(expr.left), _eliminate(expr.right)
            if _is_zero(mat) or _is_zero(vec):
                return rule("simplify.mvmul_zero", expr, ZeroVector)
            if _is_identity(mat):
                return rule("simplify.mvmul_identity", expr, vec)
            # Hat(v) * w = Cross(v, w)
            if isinstance(mat, Hat):
                return rule("simplify.mvmul_hat_to_cross", expr, Cross(mat.expr, vec))
            # (A * B) * v → A * (B * v)  (matrix associativity)
            if isinstance(mat, MMMul):
                return rule(
                    "simplify.mvmul_reassociate",
                    expr,
                    _eliminate(MVMul(mat.left, MVMul(mat.right, vec))),
                )
            # A * (B * v) → (A*B) * v  when A*B simplifies (e.g. R'*R = I)
            if isinstance(vec, MVMul):
                combined = _eliminate(MMMul(mat, vec.left))
                if not isinstance(combined, MMMul):
                    # A*B simplified (e.g. to I), so use the result
                    return rule(
                        "simplify.mvmul_combine", expr, _eliminate(MVMul(combined, vec.right))
                    )
            return MVMul(mat, vec)

        # ---- matrix * matrix ----
        case MMMul():
            l, r = _eliminate(expr.left), _eliminate(expr.right)
            if _is_zero(l) or _is_zero(r):
                return rule("simplify.mmmul_zero", expr, ZeroMatrix)
            if _is_identity(l):
                return rule("simplify.mmmul_identity", expr, r)
            if _is_identity(r):
                return rule("simplify.mmmul_identity", expr, l)
            # R^T * R = I  when R is orthogonal (SO3)
            if isinstance(l, Transpose) and l.expr.is_orthogonal and l.expr == r:
                return rule("simplify.orthogonal_transpose", expr, IdentityMatrix)
            # R * R^T = I  when R is orthogonal (SO3)
            if isinstance(r, Transpose) and r.expr.is_orthogonal and r.expr == l:
                return rule("simplify.orthogonal_transpose", expr, IdentityMatrix)
            return MMMul(l, r)

        # ---- vector * vector ----
        case VVMul():
            l, r = _eliminate(expr.left), _eliminate(expr.right)
            if _is_zero(l) or _is_zero(r):
                # inner product vᵀw → scalar 0; outer product v wᵀ → 3×3 zero
                return rule("simplify.vvmul_zero", expr, _zero_for(expr))
            return VVMul(l, r)

        # ---- binary geometry ops ----
        case Dot():
            l, r = _eliminate(expr.left), _eliminate(expr.right)
            if _is_zero(l) or _is_zero(r):
                return rule("simplify.dot_zero", expr, Zero)
            return Dot(l, r)

        case Cross():
            l, r = _eliminate(expr.left), _eliminate(expr.right)
            if _is_zero(l) or _is_zero(r):
                return rule("simplify.cross_zero", expr, ZeroVector)
            return Cross(l, r)

        # ---- unary geometry ops ----
        case Hat():
            inner = _eliminate(expr.expr)
            if _is_zero(inner):
                return rule("simplify.hat_zero", expr, ZeroMatrix)
            # Hat(Vee(X)) = X only for skew-symmetric X (Vee reads three entries of X,
            # so for a general matrix Hat(Vee(X)) is the skew part of X, e.g. Hat(Vee(I)) = 0).
            if isinstance(inner, Vee) and _is_skew_symmetric(inner.expr):
                return rule("simplify.hat_vee_skew", expr, inner.expr)
            return Hat(inner)

        case Vee():
            inner = _eliminate(expr.expr)
            if _is_zero(inner):
                return rule("simplify.vee_zero", expr, ZeroVector)
            if isinstance(inner, Hat):
                return rule("simplify.vee_hat", expr, inner.expr)
            return Vee(inner)

        case Transpose():
            inner = _eliminate(expr.expr)
            if isinstance(inner, Transpose):
                return rule("simplify.transpose_transpose", expr, inner.expr)
            # Symmetric matrix: J^T = J
            if inner.is_symmetric:
                return rule("simplify.transpose_symmetric", expr, inner)
            # Hat is skew-symmetric: Hat(v)^T = -Hat(v)
            if isinstance(inner, Hat):
                return rule(
                    "simplify.transpose_hat_skew",
                    expr,
                    SMMul(inner, Scalar("(-1)", value=-1, attr=["Constant"])),
                )
            # (s*M)^T = s*M^T
            if isinstance(inner, SMMul):
                return rule(
                    "simplify.transpose_scaled",
                    expr,
                    SMMul(_eliminate(Transpose(inner.left)), inner.right),
                )
            # a transposed scaled row is the scaled column: (aᵀ s)ᵀ = a s
            if isinstance(inner, SVMul) and isinstance(inner.left, Transpose):
                return rule(
                    "simplify.transpose_scaled_row", expr, SVMul(inner.left.expr, inner.right)
                )
            # (A*B)^T = B^T * A^T
            if isinstance(inner, MMMul):
                return rule(
                    "simplify.transpose_product",
                    expr,
                    _eliminate(MMMul(Transpose(inner.right), Transpose(inner.left))),
                )
            return Transpose(inner)

        # ---- calculus ops ----
        case Variation():
            inner = _eliminate(expr.expr)
            if _is_zero(inner):
                return rule("simplify.variation_zero", expr, _zero_for(inner))
            return Variation(inner)

        case TimeDerivative():
            inner = _eliminate(expr.expr)
            if _is_zero(inner):
                return rule("simplify.time_derivative_zero", expr, _zero_for(inner))
            # d/dt(∫x dt) = x
            if isinstance(inner, TimeIntegral):
                return rule("simplify.time_derivative_integral", expr, inner.expr)
            return TimeDerivative(inner)

        case TimeIntegral():
            inner = _eliminate(expr.expr)
            if _is_zero(inner):
                return rule("simplify.time_integral_zero", expr, _zero_for(inner))
            # ∫(d/dt(x)) dt = x
            if isinstance(inner, TimeDerivative):
                return rule("simplify.time_integral_derivative", expr, inner.expr)
            return TimeIntegral(inner)

        # ---- leaves ----
        case _:
            return expr


# ---------------------------------------------------------------------------
# Predicates
# ---------------------------------------------------------------------------


def _is_zero(expr) -> bool:
    return getattr(expr, "is_zero", False)


def _is_one(expr) -> bool:
    return getattr(expr, "value", None) == 1


def _is_identity(expr) -> bool:
    flags = getattr(expr, "flags", None)
    return flags is not None and flags.is_identity


def _is_skew_symmetric(expr) -> bool:
    """True for Hat(·), a SkewSymmMatrix leaf, and SMMul/MAdd/Transpose of such."""
    return isinstance(expr, Hat) or getattr(expr, "is_skew_symmetric", False)


def _is_numeric_leaf(expr) -> bool:
    """True for leaf Scalars with a numeric value (not operation nodes)."""
    return isinstance(expr, Scalar) and not hasattr(expr, "nodes") and expr.value is not None


def _zero_for(expr):
    """Return the appropriate typed zero constant."""
    match getattr(expr, "type", None):
        case ExprType.VECTOR:
            return ZeroVector
        case ExprType.MATRIX:
            return ZeroMatrix
        case _:
            return Zero


# ---------------------------------------------------------------------------
# Per-node helpers
# ---------------------------------------------------------------------------


# Relative size below which a sum of float coefficients is roundoff of an exact zero:
# x*0.1*3 - x*0.3 sums 0.30000000000000004 - 0.3 = 5.55e-17, which is not a term.
_CANCEL_RTOL = 1e-12


def _cancelled(total, scale):
    """*total* with float cancellation residue (``|total| <= 1e-12·scale``) set to 0."""
    if isinstance(total, float) and abs(total) <= _CANCEL_RTOL * scale:
        return 0
    return total


def _simplify_add(nodes):
    """Scalar addition: filter zeros, combine numerics, collect like terms."""
    numeric_total = 0
    scale = 0
    kept = []
    for n in nodes:
        if _is_zero(n):
            continue
        if _is_numeric_leaf(n):
            numeric_total += n.value
            scale = max(scale, abs(n.value))
        else:
            kept.append(n)
    numeric_total = _cancelled(numeric_total, scale)
    if numeric_total != 0:
        kept.append(Number(numeric_total))
    kept = _collect_like_terms(kept, ExprType.SCALAR)
    if len(kept) == 0:
        return Zero
    if len(kept) == 1:
        return kept[0]
    return Add(*kept)


def _simplify_vadd(nodes):
    """Vector addition: filter zeros, collect like terms."""
    kept = [n for n in nodes if not _is_zero(n)]
    kept = _collect_like_terms(kept, ExprType.VECTOR)
    if len(kept) == 0:
        return ZeroVector
    if len(kept) == 1:
        return kept[0]
    return VAdd(*kept)


def _simplify_madd(nodes):
    """Matrix addition: filter zeros, collect like terms."""
    kept = [n for n in nodes if not _is_zero(n)]
    kept = _collect_like_terms(kept, ExprType.MATRIX)
    if len(kept) == 0:
        return ZeroMatrix
    if len(kept) == 1:
        return kept[0]
    return MAdd(*kept)


# ---------------------------------------------------------------------------
# Like-term collection
# ---------------------------------------------------------------------------


def _split_coeff(expr):
    """Split expr into (numeric_coefficient, term).

    For Mul(term, Number) or Mul(Number, term), returns (number, term).
    For SVMul(vec, Number) or SMMul(mat, Number), returns (number, vec/mat).
    Otherwise returns (1, expr).
    """
    if isinstance(expr, Mul):
        if _is_numeric_leaf(expr.left):
            return expr.left.value, expr.right
        if _is_numeric_leaf(expr.right):
            return expr.right.value, expr.left
    if isinstance(expr, SVMul):
        if _is_numeric_leaf(expr.right):
            return expr.right.value, expr.left
    if isinstance(expr, SMMul):
        if _is_numeric_leaf(expr.right):
            return expr.right.value, expr.left
    if _is_numeric_leaf(expr):
        return expr.value, None
    return 1, expr


def _commutative_equal(a, b):
    """Structural equality that respects commutativity of Dot and Mul."""
    if a is b:
        return True
    if type(a) is not type(b):
        return False

    # Commutative binary ops: check both orderings
    if isinstance(a, Dot):
        return (_commutative_equal(a.left, b.left) and _commutative_equal(a.right, b.right)) or (
            _commutative_equal(a.left, b.right) and _commutative_equal(a.right, b.left)
        )
    if isinstance(a, Mul):
        return (_commutative_equal(a.left, b.left) and _commutative_equal(a.right, b.right)) or (
            _commutative_equal(a.left, b.right) and _commutative_equal(a.right, b.left)
        )

    # Non-commutative binary ops: order matters
    if isinstance(a, (SVMul, SMMul, MVMul, MMMul, VVMul, Cross)):
        return _commutative_equal(a.left, b.left) and _commutative_equal(a.right, b.right)

    # N-ary: same length and all children match in order
    a_nodes = getattr(a, "nodes", None)
    b_nodes = getattr(b, "nodes", None)
    if a_nodes is not None and b_nodes is not None:
        if len(a_nodes) != len(b_nodes):
            return False
        return all(_commutative_equal(an, bn) for an, bn in zip(a_nodes, b_nodes))

    # Unary ops
    a_inner = getattr(a, "expr", None)
    b_inner = getattr(b, "expr", None)
    if a_inner is not None and b_inner is not None:
        return _commutative_equal(a_inner, b_inner)

    # Leaves: fall back to string comparison
    return str(a) == str(b)


def _make_term(coeff, term, expr_type):
    """Reconstruct coeff * term for the given type."""
    if term is None:
        # Pure numeric
        return Number(coeff)
    if coeff == 0:
        return None
    if coeff == 1:
        return term
    scalar = Number(coeff)
    if expr_type == ExprType.VECTOR:
        return SVMul(term, scalar)
    if expr_type == ExprType.MATRIX:
        return SMMul(term, scalar)
    return Mul(term, scalar)


def _collect_like_terms(nodes, expr_type):
    """Group terms by commutative equality and sum their coefficients."""
    if len(nodes) <= 1:
        return nodes

    # Split each node into (coeff, term)
    pairs = [_split_coeff(n) for n in nodes]

    # Group by commutative equality
    merged_coeffs = list(range(len(pairs)))  # index of canonical representative
    coeffs = [c for c, _ in pairs]
    scales = [abs(c) for c in coeffs]
    terms = [t for _, t in pairs]

    for i in range(len(pairs)):
        if merged_coeffs[i] != i:
            continue  # already merged into an earlier group
        for j in range(i + 1, len(pairs)):
            if merged_coeffs[j] != j:
                continue  # already merged
            ti = terms[i]
            tj = terms[j]
            # Both pure numeric (term is None) — already combined in _simplify_add
            if ti is None and tj is None:
                continue
            if ti is not None and tj is not None and _commutative_equal(ti, tj):
                coeffs[i] += coeffs[j]
                scales[i] = max(scales[i], abs(coeffs[j]))
                coeffs[j] = 0
                merged_coeffs[j] = i

    # Reconstruct
    result = []
    for i in range(len(pairs)):
        if merged_coeffs[i] != i:
            continue
        made = _make_term(_cancelled(coeffs[i], scales[i]), terms[i], expr_type)
        if made is not None:
            result.append(made)
    return result


def _simplify_mul(l, r):
    """Scalar multiplication: zeros, identity, numeric combining."""
    if _is_zero(l) or _is_zero(r):
        return Zero
    if _is_one(l):
        return r
    if _is_one(r):
        return l
    # both numeric leaves → combine (a zero product is the Zero leaf, not an unflagged "(0)")
    if _is_numeric_leaf(l) and _is_numeric_leaf(r):
        combined = l.value * r.value
        if combined == 0:
            return Zero
        return Number(combined)
    # Fold numeric constants through nested Muls:
    # e.g. Mul(-0.5, Mul(m, 2)) → Mul(m, -1.0)
    result = Mul(l, r)
    factors = []
    _flatten_mul(result, factors)
    numerics = [f for f in factors if _is_numeric_leaf(f)]
    if len(numerics) >= 2:
        others = [f for f in factors if not _is_numeric_leaf(f)]
        combined = 1
        for n in numerics:
            combined *= n.value
        if combined == 0:
            return Zero
        if combined != 1:
            others.append(Number(combined))
        if len(others) == 0:
            return Number(combined)
        # Rebuild left-associative chain
        out = others[0]
        for o in others[1:]:
            out = Mul(out, o)
        return out
    return result


def _flatten_mul(expr, factors):
    """Flatten nested Mul into a list of factors."""
    if isinstance(expr, Mul):
        _flatten_mul(expr.left, factors)
        _flatten_mul(expr.right, factors)
    else:
        factors.append(expr)
