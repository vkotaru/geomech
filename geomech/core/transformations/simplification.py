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

Like terms (all three folds): an addend is split into a number (the product of every
numeric literal in its flattened ``Mul`` chain and, for vectors and matrices, in the
scalars of its ``SVMul``/``SMMul``), its other scalar factors and its vector/matrix base;
two addends merge when their bases are equal and their factors are equal as multisets
(``x(2)y + xy(3) → xy(5)``, ``v a + v a(-1) → 0``, ``(vᵀw)s + s(vᵀw) → (vᵀw)s(2)``).
"Equal" is ``term_key``: the operands of a ``Dot`` and of an inner product ``aᵀb`` commute
(``aᵀb s − bᵀa s → 0``) and so do the factors of a ``Mul``; every other product keeps its
order. Grouping is one pass over the addends, keyed by ``term_key``.
Numeric literals are anonymous numbers only: a named ``Scalar("m", value=2)`` is a symbol
and is never folded (decision D1); its value is used by ``geomech.numeric`` only.
  simplify.mul_fold                Mul: a zero factor gives 0, a factor one is dropped,
                                   numeric factors of a Mul chain are multiplied into one
                                   number placed last; a float with an integer value
                                   becomes an int (``(-1.0) → (-1)``)

Every numeric literal leaf is put in the spelling ``Number`` gives it (``(2.0)`` → ``(2)``,
``Scalar("0.5", value=0.5)`` → ``(0.5)``), so equal numbers print, and compare, the same.
  simplify.svmul_zero              v·0 = 0·s = 0 (vector)
  simplify.svmul_one               v·1 = v
  simplify.smmul_zero              M·0 = 0·s = 0 (matrix)
  simplify.smmul_one               M·1 = M
  simplify.mvmul_zero              0·v = M·0 = 0
  simplify.mvmul_identity          I·v = v                    left factor flagged identity
                                                              (or its transpose)
  simplify.mvmul_hat_to_cross      Hat(a)·v = Cross(a, v)
  simplify.mvmul_reassociate       (A·B)·v = A·(B·v)
  simplify.mvmul_combine           A·(B·v) = (A·B)·v          only when A·B simplifies to a
                                                              non-product (e.g. RᵀR = I)
  simplify.mmmul_zero              0·B = A·0 = 0
  simplify.mmmul_identity          I·B = B, A·I = A           factor flagged identity (or
                                                              its transpose)
  simplify.mmmul_reassociate       (A·B)·C = A·(B·C)
  simplify.mmmul_combine           A·(B·C) = (A·B)·C          only when A·B simplifies to a
                                                              non-product (e.g. RᵀR = I)
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
  simplify.transpose_symmetric     Jᵀ = J                     J flagged symmetric (the
                                                              identity ``I`` is)
  simplify.transpose_hat_skew      Hat(v)ᵀ = -Hat(v)
  simplify.transpose_skew          Sᵀ = -S                    S skew-symmetric (a leaf flagged
                                                              SkewSymmetricMatrix, a sum or
                                                              transpose of such), not a Hat
                                                              or a scaled matrix
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

from geomech.core.base.constants import NEG_ONE, zero_like
from geomech.core.base.expressions import (
    IdentityMatrix,
    Number,
    Zero,
    ZeroMatrix,
    ZeroVector,
    is_number_literal,
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
from geomech.core.transformations.expansion import expand
from geomech.core.transformations.pulling import pull
from geomech.core.transformations.utils import calculus_node, flatten_mul
from geomech.core.transformations.vector_identities import vector_rules
from geomech.utils import trace as _trace
from geomech.utils.trace import rule

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def simplify(expr):
    """Pull scalars, apply vector identities, then eliminate zeros and
    combine numeric constants — single call does everything."""
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
    """Bottom-up pass: zeros, identities, numeric constant combining.

    One function per node class (``_elim_add``, ``_elim_mvmul``, ...), picked from
    ``_ELIMINATE`` by the class of *expr*; a leaf is returned as it is, except a numeric
    literal, which is put in its normal form (``_canonical_literal``).
    """
    handler = _ELIMINATE.get(type(expr))
    return _canonical_literal(expr) if handler is None else handler(expr)


def _elim_add(expr):
    children = [_eliminate(n) for n in expr.nodes]
    out = _simplify_add(children)
    if _trace._TRACE_RULES:
        rule("simplify.add_fold", _folded_from(Add, children, expr), out)
    return out


def _elim_vadd(expr):
    children = [_eliminate(n) for n in expr.nodes]
    out = _simplify_vadd(children)
    if _trace._TRACE_RULES:
        rule("simplify.vadd_fold", _folded_from(VAdd, children, expr), out)
    return out


def _elim_madd(expr):
    children = [_eliminate(n) for n in expr.nodes]
    out = _simplify_madd(children)
    if _trace._TRACE_RULES:
        rule("simplify.madd_fold", _folded_from(MAdd, children, expr), out)
    return out


def _elim_mul(expr):
    children = [_eliminate(expr.left), _eliminate(expr.right)]
    out = _simplify_mul(*children)
    if _trace._TRACE_RULES:
        rule("simplify.mul_fold", _folded_from(Mul, children, expr), out)
    return out


def _elim_svmul(expr):
    vec, scl = _eliminate(expr.left), _eliminate(expr.right)
    if _is_zero(scl) or _is_zero(vec):
        return rule("simplify.svmul_zero", expr, ZeroVector)
    if _is_one(scl):
        return rule("simplify.svmul_one", expr, vec)
    return SVMul(vec, _canonical_literal(scl))


def _elim_smmul(expr):
    mat, scl = _eliminate(expr.left), _eliminate(expr.right)
    if _is_zero(scl) or _is_zero(mat):
        return rule("simplify.smmul_zero", expr, ZeroMatrix)
    if _is_one(scl):
        return rule("simplify.smmul_one", expr, mat)
    return SMMul(mat, _canonical_literal(scl))


def _elim_mvmul(expr):
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
            return rule("simplify.mvmul_combine", expr, _eliminate(MVMul(combined, vec.right)))
    return MVMul(mat, vec)


def _elim_mmmul(expr):
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
    # (A * B) * C → A * (B * C)  (matrix associativity, as for MVMul)
    if isinstance(l, MMMul):
        return rule(
            "simplify.mmmul_reassociate",
            expr,
            _eliminate(MMMul(l.left, MMMul(l.right, r))),
        )
    # A * (B * C) → (A*B) * C  when A*B simplifies (e.g. R'*R = I)
    if isinstance(r, MMMul):
        combined = _eliminate(MMMul(l, r.left))
        if not isinstance(combined, MMMul):
            return rule("simplify.mmmul_combine", expr, _eliminate(MMMul(combined, r.right)))
    return MMMul(l, r)


def _elim_vvmul(expr):
    l, r = _eliminate(expr.left), _eliminate(expr.right)
    if _is_zero(l) or _is_zero(r):
        # inner product vᵀw → scalar 0; outer product v wᵀ → 3×3 zero
        return rule("simplify.vvmul_zero", expr, zero_like(expr))
    return VVMul(l, r)


def _elim_dot(expr):
    l, r = _eliminate(expr.left), _eliminate(expr.right)
    if _is_zero(l) or _is_zero(r):
        return rule("simplify.dot_zero", expr, Zero)
    return Dot(l, r)


def _elim_cross(expr):
    l, r = _eliminate(expr.left), _eliminate(expr.right)
    if _is_zero(l) or _is_zero(r):
        return rule("simplify.cross_zero", expr, ZeroVector)
    return Cross(l, r)


def _elim_hat(expr):
    inner = _eliminate(expr.expr)
    if _is_zero(inner):
        return rule("simplify.hat_zero", expr, ZeroMatrix)
    # Hat(Vee(X)) = X only for skew-symmetric X (Vee reads three entries of X,
    # so for a general matrix Hat(Vee(X)) is the skew part of X, e.g. Hat(Vee(I)) = 0).
    if isinstance(inner, Vee) and _is_skew_symmetric(inner.expr):
        return rule("simplify.hat_vee_skew", expr, inner.expr)
    return Hat(inner)


def _elim_vee(expr):
    inner = _eliminate(expr.expr)
    if _is_zero(inner):
        return rule("simplify.vee_zero", expr, ZeroVector)
    if isinstance(inner, Hat):
        return rule("simplify.vee_hat", expr, inner.expr)
    return Vee(inner)


def _elim_transpose(expr):
    inner = _eliminate(expr.expr)
    if isinstance(inner, Transpose):
        return rule("simplify.transpose_transpose", expr, inner.expr)
    # Symmetric matrix: J^T = J
    if inner.is_symmetric:
        return rule("simplify.transpose_symmetric", expr, inner)
    # Hat is skew-symmetric: Hat(v)^T = -Hat(v)
    if isinstance(inner, Hat):
        return rule("simplify.transpose_hat_skew", expr, SMMul(inner, NEG_ONE))
    # any other skew-symmetric matrix: S^T = -S
    if _is_skew_symmetric(inner) and not isinstance(inner, SMMul):
        return rule("simplify.transpose_skew", expr, SMMul(inner, NEG_ONE))
    # (s*M)^T = s*M^T
    if isinstance(inner, SMMul):
        return rule(
            "simplify.transpose_scaled",
            expr,
            SMMul(_eliminate(Transpose(inner.left)), inner.right),
        )
    # a transposed scaled row is the scaled column: (aᵀ s)ᵀ = a s
    if isinstance(inner, SVMul) and isinstance(inner.left, Transpose):
        return rule("simplify.transpose_scaled_row", expr, SVMul(inner.left.expr, inner.right))
    # (A*B)^T = B^T * A^T
    if isinstance(inner, MMMul):
        return rule(
            "simplify.transpose_product",
            expr,
            _eliminate(MMMul(Transpose(inner.right), Transpose(inner.left))),
        )
    return Transpose(inner)


def _elim_variation(expr):
    inner = _eliminate(expr.expr)
    if _is_zero(inner):
        return rule("simplify.variation_zero", expr, zero_like(inner))
    return calculus_node(Variation, inner)


def _elim_time_derivative(expr):
    inner = _eliminate(expr.expr)
    if _is_zero(inner):
        return rule("simplify.time_derivative_zero", expr, zero_like(inner))
    # d/dt(∫x dt) = x
    if isinstance(inner, TimeIntegral):
        return rule("simplify.time_derivative_integral", expr, inner.expr)
    return calculus_node(TimeDerivative, inner)


def _elim_time_integral(expr):
    inner = _eliminate(expr.expr)
    if _is_zero(inner):
        return rule("simplify.time_integral_zero", expr, zero_like(inner))
    # ∫(d/dt(x)) dt = x
    if isinstance(inner, TimeDerivative):
        return rule("simplify.time_integral_derivative", expr, inner.expr)
    return TimeIntegral(inner)


# ---------------------------------------------------------------------------
# Predicates
# ---------------------------------------------------------------------------


def _is_zero(expr) -> bool:
    return getattr(expr, "is_zero", False)


def _is_one(expr) -> bool:
    """A numeric literal equal to 1 (a named ``Scalar("m", value=1)`` is a symbol, D1)."""
    return is_number_literal(expr) and expr.value == 1


def _is_identity(expr) -> bool:
    """A leaf flagged ``Identity``, or the transpose of one (``Iᵀ = I``)."""
    if isinstance(expr, Transpose):
        return _is_identity(expr.expr)
    flags = getattr(expr, "flags", None)
    return flags is not None and flags.is_identity


def _is_skew_symmetric(expr) -> bool:
    """True for Hat(·), a SkewSymmMatrix leaf, and SMMul/MAdd/Transpose of such."""
    return isinstance(expr, Hat) or getattr(expr, "is_skew_symmetric", False)


def _is_numeric_leaf(expr) -> bool:
    """True for an anonymous numeric literal (``Number``, ``Zero``, ``2 * x``'s ``(2)``).

    A named scalar with a value (``Scalar("m", value=2)``) is not one: it stays a symbol
    in every simplification (fix plan decision D1).
    """
    return is_number_literal(expr)


def _exact(value):
    """*value* as an ``int`` when it is a float with an exact integer value (``-1.0 → -1``).

    ``0.5 * 2 * -1`` folds to ``-1.0``; printing ``(-1)`` instead of ``(-1.0)`` makes the
    term equal to the ``(-1)`` terms it should merge with. Non-integral floats are kept.
    """
    if isinstance(value, float) and value.is_integer() and abs(value) < 2.0**53:
        return int(value)
    return value


def _canonical_literal(expr):
    """A numeric literal in the one spelling ``Number`` gives it; anything else unchanged.

    The value goes through ``_exact`` (``(2.0)`` → ``(2)``) and the name becomes the
    ``Number`` form, so the literal one half prints ``(0.5)`` whether it was built as
    ``Number(0.5)``, ``0.5 * x`` or ``Scalar("0.5", value=0.5)``. ``Zero`` and ``One``
    (printed ``0`` and ``1``) keep their names. A named scalar with a value
    (``Scalar("m", value=2)``) is not a literal and is never touched (decision D1).
    """
    if not is_number_literal(expr) or expr.name in ("0", "1"):
        return expr
    canonical = Number(_exact(expr.value))
    return expr if canonical.name == expr.name else canonical


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
    numeric_total = _exact(_cancelled(numeric_total, scale))
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


def _split_term(expr, expr_type):
    """Split an addend into ``(coefficient, base, factors)``.

    *coefficient* is the product of every numeric literal of the term, *factors* the
    other scalar factors (a flattened ``Mul`` chain, in order) and *base* the vector or
    matrix the scalars multiply (``None`` for a scalar term). ``x(2)y`` gives
    ``(2, None, [x, y])``, ``v(s(-1))`` gives ``(-1, v, [s])``. Two addends are like terms
    when their bases are equal and their factors are equal as multisets.
    """
    scalars = []
    base = None
    if expr_type == ExprType.SCALAR:
        scalars.append(expr)
    else:
        node = expr
        cls = SVMul if expr_type == ExprType.VECTOR else SMMul
        while isinstance(node, cls):
            scalars.append(node.right)
            node = node.left
        base = node
    factors = [f for sc in scalars for f in flatten_mul(sc)]
    coeff = 1
    rest = []
    for f in factors:
        if _is_numeric_leaf(f):
            coeff *= f.value
        else:
            rest.append(f)
    return _exact(coeff), base, rest


def term_key(expr) -> str:
    """A canonical string for *expr* up to the commutativity the like-term merge knows.

    Two expressions have the same key exactly when they are the same tree except that the
    operands of a ``Dot`` may be swapped, the factors of a ``Mul`` chain may be reordered,
    and an inner product ``aᵀb`` of two columns (``VVMul(Transpose(a), b)``) may be
    written ``bᵀa`` (at any depth). Every other node compares its operands in order, and a leaf
    by class and printed name. The key is injective on such classes: a leaf is written
    ``Class#len#name`` (the length makes it self-delimiting), a node
    ``Class(key,key,...)`` with the operand keys sorted for ``Dot`` and the flattened
    ``Mul`` factors, and an inner product ``VVMul.inner(key,key)`` with the keys of its
    two columns sorted.
    """
    nodes = getattr(expr, "nodes", None)
    cls = type(expr).__name__
    if not nodes:
        name = str(expr)
        return f"{cls}#{len(name)}#{name}"
    if isinstance(expr, Mul):
        keys = sorted(term_key(f) for f in flatten_mul(expr))
    elif (
        isinstance(expr, VVMul)
        and expr.type == ExprType.SCALAR
        and isinstance(expr.left, Transpose)
    ):
        # the inner product aᵀb of two columns is symmetric: aᵀb = bᵀa
        keys = sorted((term_key(expr.left.expr), term_key(expr.right)))
        return "VVMul.inner(" + ",".join(keys) + ")"
    elif isinstance(expr, Dot):
        keys = sorted(term_key(n) for n in nodes)
    else:
        keys = [term_key(n) for n in nodes]
    return cls + "(" + ",".join(keys) + ")"


def _commutative_equal(a, b):
    """Equality up to the commutativity of ``Dot``, ``Mul`` and inner products (``term_key``)."""
    return a is b or term_key(a) == term_key(b)


def _like_key(base, factors):
    """The like-term key of an addend split as (base, factors): equal keys merge."""
    return (
        None if base is None else term_key(base),
        tuple(sorted(term_key(f) for f in factors)),
    )


def _chain(factors):
    """Left-associated ``Mul`` chain of *factors* (at least one)."""
    out = factors[0]
    for f in factors[1:]:
        out = Mul(out, f)
    return out


def _make_term(coeff, base, factors, expr_type):
    """Rebuild ``coeff · factors · base`` (the inverse of ``_split_term``); ``None`` for 0."""
    if coeff == 0:
        return None
    scalars = list(factors)
    if coeff != 1:
        scalars.append(Number(coeff))
    if expr_type == ExprType.SCALAR:
        return _chain(scalars) if scalars else Number(1)
    if not scalars:
        return base
    cls = SVMul if expr_type == ExprType.VECTOR else SMMul
    return cls(base, _chain(scalars))


def _collect_like_terms(nodes, expr_type):
    """Merge like terms (same base, same symbolic factors) by summing their coefficients.

    A term that merges with nothing is kept as it is; a merged group is rebuilt from its
    first term with the summed coefficient (dropped when the sum is 0). Terms are grouped
    by ``_like_key`` in one pass (linear in the number of addends); the groups keep the
    order of their first terms and their coefficients are summed in addend order.
    """
    if len(nodes) <= 1:
        return nodes

    split = [_split_term(n, expr_type) for n in nodes]
    first = {}  # like key → index of the group's first term
    group = []  # index of each term's group (its first term)
    coeffs = [c for c, _, _ in split]
    scales = [abs(c) for c in coeffs]
    merged = [False] * len(nodes)

    for j, (_, base, factors) in enumerate(split):
        i = first.setdefault(_like_key(base, factors), j)
        group.append(i)
        if i != j:
            coeffs[i] += coeffs[j]
            scales[i] = max(scales[i], abs(coeffs[j]))
            merged[i] = True

    result = []
    for i in range(len(nodes)):
        if group[i] != i:
            continue
        if not merged[i]:
            result.append(nodes[i])
            continue
        _, base, factors = split[i]
        made = _make_term(_exact(_cancelled(coeffs[i], scales[i])), base, factors, expr_type)
        if made is not None:
            result.append(made)
    return result


def _simplify_mul(l, r):
    """Scalar multiplication: zeros, identity, numeric combining.

    The chain ``l·r`` is flattened; its numeric literals are multiplied into one number
    placed last (``(-0.5)·(m·2) → m(-1)``), the other factors keep their order.
    """
    if _is_zero(l) or _is_zero(r):
        return Zero
    if _is_one(l):
        return r
    if _is_one(r):
        return l
    factors = flatten_mul(Mul(l, r))
    numerics = [f for f in factors if _is_numeric_leaf(f)]
    if not numerics:
        return Mul(l, r)
    others = [f for f in factors if not _is_numeric_leaf(f)]
    combined = 1
    for n in numerics:
        combined *= n.value
    combined = _exact(combined)
    if combined == 0:
        return Zero
    if not others:
        return Number(combined)
    if combined == 1:
        return _chain(others)
    if (
        len(numerics) == 1
        and factors[-1] is numerics[0]
        and type(numerics[0].value) is type(combined)
        and numerics[0].value == combined
    ):
        return Mul(l, r)  # already in normal form: keep the tree as built
    return Mul(_chain(others), Number(combined))


#: ``_eliminate``'s handler for each node class (a leaf has none and is kept).
#: Looked up by the exact class (``type(expr)``), not along the MRO: a subclass of a
#: node class is not dispatched and is kept as it is, like a leaf. No such subclass exists;
#: one added later needs its own entry here (pinned in test_phase5b_structure.py).
_ELIMINATE = {
    Add: _elim_add,
    VAdd: _elim_vadd,
    MAdd: _elim_madd,
    Mul: _elim_mul,
    SVMul: _elim_svmul,
    SMMul: _elim_smmul,
    MVMul: _elim_mvmul,
    MMMul: _elim_mmmul,
    VVMul: _elim_vvmul,
    Dot: _elim_dot,
    Cross: _elim_cross,
    Hat: _elim_hat,
    Vee: _elim_vee,
    Transpose: _elim_transpose,
    Variation: _elim_variation,
    TimeDerivative: _elim_time_derivative,
    TimeIntegral: _elim_time_integral,
}
