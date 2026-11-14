"""Collect: rearrange a scalar expression so that a target vector
appears on the left side of Dot.

This is a preprocessing step for integration by parts (ibp), which
only matches ``Dot(target, ...)`` on the left.

Rules (name, identity, hypothesis the code checks):

  collect.inner_product_to_dot   VVMul(rᵀ, b) = Dot(r, b), with r = column_of(rᵀ)
                                 (``(vᵀ)ᵀ = v``, ``(vᵀM)ᵀ = Mᵀv``, ``(vᵀs)ᵀ = v s``);
                                 only for an inner product that contains the target, so
                                 ``v.T() * A * v`` in a kinetic energy can be integrated
                                 by parts
  collect.flip_target_right      Dot(a, v) = Dot(v, a)    the target is the right operand
  collect.flip_cross_mvmul       Dot(b×c, M x) = Dot(M x, b×c)
                                                         the target is on either side; the
                                                         MVMul goes left for the rules below
  collect.flip_mvmul_pair        Dot(M x, N f(v)) = Dot(N f(v), M x)
                                                         the target is in the right MVMul
  collect.flip_mvmul_left        Dot(M x, f(v)) = Dot(f(v), M x)
                                                         the target is in the non-MVMul side
  collect.flip_mvmul_right       Dot(a, M f(v)) = Dot(M f(v), a)
                                                         the target is in the right MVMul
  collect.flip_cross_left        Dot(f(v)×c, a) = Dot(a, f(v)×c)
                                                         the target is in the left Cross
                                 all six are commutativity of Dot (no hypothesis); the
                                 names say which call site flipped
  collect.triple_product_left    a·(v×c) = v·(c×a)       the target v is the left operand
                                                         of the Cross
  collect.triple_product_right   a·(b×v) = v·(a×b)       the target is the right operand
  collect.triple_product_left_nested
                                 a·(f(v)×c) = f(v)·(c×a)  the target is inside the left
                                                         operand (collect then recurses)
  collect.triple_product_right_nested
                                 a·(b×f(v)) = f(v)·(a×b)  inside the right operand
  collect.triple_product_mvmul_left
  collect.triple_product_mvmul_right
  collect.triple_product_mvmul_left_nested
  collect.triple_product_mvmul_right_nested
                                 the same four rotations when the other operand of the
                                 Dot is a matrix product, a = M x (the cases of
                                 Dot(MVMul, Cross)); all eight are the cyclic scalar
                                 triple product, so none needs a hypothesis
  collect.peel_matrix            Dot(M f, w) = Dot(f, Mᵀ w)
                                 the target is in the vector operand f, not in M
  collect.strip_rotation         Dot(R a, R b) = Dot(a, b)
                                 both matrices are SO3 points equal by name

A target that occurs on both sides of a Dot is not linear; ``collect`` raises
``NonLinearError`` instead of looking for a rotation that cannot exist.  A target that
sits inside the *matrix* factor of an ``MVMul`` (``Hat(v)``, ``v cᵀ``, ``M Hat(v)``) cannot
be exposed by transposition, which only moves the matrix: ``collect`` raises
``AlgebraicError`` naming the term and the target.
"""

from __future__ import annotations

from geomech.core.base.expressions import SO3
from geomech.core.base.types import ExprType
from geomech.core.operations.addition import Add
from geomech.core.operations.geometry import Cross, Dot, Transpose, column_of
from geomech.core.operations.multiplication import Mul, MVMul, VVMul
from geomech.utils.errors import AlgebraicError, NonLinearError
from geomech.utils.trace import rule


def collect(expr, vec):
    """Collect *expr* with respect to *vec*.

    Rearranges so *vec* appears on the left of every Dot that contains it. Only a scalar
    expression can hold a ``Dot``; anything else raises ``TypeError``.
    """
    if getattr(expr, "type", None) != ExprType.SCALAR:
        kind = getattr(getattr(expr, "type", None), "name", None)
        raise TypeError(
            f"collect: expected a scalar expression, got {type(expr).__name__} '{expr}'"
            f" of type {kind}"
        )
    match expr:
        case Add(nodes=nodes):
            return Add(*[collect(n, vec) for n in nodes])

        case Mul():
            return Mul(collect(expr.left, vec), collect(expr.right, vec))

        case Dot():
            return _collect_dot(expr, vec)

        case VVMul() if expr.has(vec):
            # an inner product rᵀb holding the target: as Dot(r, b) the Dot rules apply
            return _collect_dot(
                rule("collect.inner_product_to_dot", expr, Dot(column_of(expr.left), expr.right)),
                vec,
            )

        case _:
            return expr


def _peel(l, r, vec, expr):
    """``Dot(M f(vec), w) → Dot(f(vec), Mᵀ w)`` for ``l = MVMul(M, f(vec))``, then recurse.

    Only fires when *vec* is in the vector operand ``l.right``. A target inside the
    matrix factor ``l.left`` cannot be exposed by transposition (it would only move the
    matrix to the other side, and the flip rule would move it back), so that case raises.
    """
    if l.left.has(vec):
        if l.right.has(vec):
            raise NonLinearError(
                f"collect: '{vec}' appears in both factors of '{l}' in '{expr}'; the term is"
                " not linear in the target"
            )
        raise AlgebraicError(
            f"collect: the target '{vec}' occurs inside the matrix factor '{l.left}' of the"
            f" term '{expr}'; transposition exposes a target only in the vector operand"
            " ((M f(v))·w = f(v)·(Mᵀ w)), so this term cannot be integrated by parts"
        )
    return collect(
        rule("collect.peel_matrix", expr, Dot(l.right, MVMul(Transpose(l.left), r))), vec
    )


def _collect_dot(expr, vec):
    """Rearrange a Dot so *vec* is on the left.

    Termination.  After the both-sides check exactly one side holds *vec*; let ``d`` be
    the number of nodes on the path from that side's root down to *vec*. Every rule is one
    of:

    - terminal: returns a Dot with *vec* on the left, returns *expr* unchanged, or raises;
    - a flip ``Dot(l, r) → Dot(r, l)``, which leaves ``d`` unchanged;
    - a peel ``Dot(M f(vec), w) → Dot(f(vec), Mᵀ w)``, a triple-product rotation
      ``Dot(a, Cross(f(vec), c)) → Dot(f(vec), Cross(c, a))`` (or its mirror) or the
      stripping of a common rotation ``Dot(R a, R b) → Dot(a, b)``, each of which removes
      the node directly above *vec*'s branch and so decreases ``d`` by one. The side that
      receives a new ``MVMul``/``Cross`` layer never holds *vec*.

    The peel is guarded (``_peel``): it fires only when *vec* is in the vector operand of
    the ``MVMul``. A target inside the matrix factor raises instead of being moved to the
    other side, so a peeled term cannot bounce back through the flip rules. A flip is never
    followed by another flip: ``(Cross, MVMul)`` flips into ``(MVMul, Cross)``, which
    rotates, peels or terminates; ``(MVMul, MVMul)`` with *vec* on the right flips into the
    same case with *vec* on the left, which strips, peels or terminates; ``(MVMul, other)``
    with *vec* in ``other`` (not a ``Cross`` or ``MVMul``) flips into ``(other, MVMul)`` with
    *vec* on the left, which only the terminal default matches; ``(other, MVMul)`` and
    ``(Cross, other)`` with *vec* in the named side flip into cases that do not flip. Hence
    at most ``2d + 1`` steps.
    """
    l, r = expr.left, expr.right

    # vec on both sides: the Dot is at least quadratic in vec; no rotation can
    # expose a single linear occurrence (and flipping would recurse forever).  This
    # runs before the shortcuts below so that Dot(vec, f(vec)) and Dot(f(vec), vec)
    # are reported, not flipped into a "canonical" form.
    if l.has(vec) and r.has(vec):
        raise NonLinearError(
            f"collect: '{vec}' appears on both sides of '{expr}'; the expression is not"
            " linear in the target"
        )

    # Already canonical
    if l == vec:
        return expr

    # Simple flip
    if r == vec:
        return rule("collect.flip_target_right", expr, Dot(r, l))

    # Neither side contains vec — nothing to do
    if not l.has(vec) and not r.has(vec):
        return expr

    # --- structural patterns (exactly one side has vec) ---

    match (l, r):
        # MVMul + Cross: kinetic energy patterns
        case (Cross(), MVMul()):
            # Flip so MVMul is on the left, then recurse
            return collect(rule("collect.flip_cross_mvmul", expr, Dot(r, l)), vec)

        case (MVMul(), Cross()):
            # Dot(M*a, Cross(b, c)) — use triple product to put vec on left
            if r.left == vec:
                # Dot(M*a, Cross(vec, c)) = Dot(vec, Cross(c, M*a))
                return rule("collect.triple_product_mvmul_left", expr, Dot(vec, Cross(r.right, l)))
            if r.right == vec:
                # Dot(M*a, Cross(b, vec)) = Dot(vec, Cross(M*a, b))
                return rule("collect.triple_product_mvmul_right", expr, Dot(vec, Cross(l, r.left)))
            if l.has(vec):
                # vec is inside the MVMul (directly as l.right or deeper): peel the
                # matrix, Dot(M*f(vec), b×c) = Dot(f(vec), Mᵀ(b×c)), and recurse.  Note
                # that (M v)·(b×c) is NOT v·((M b)×c); the only valid move is Mᵀ.
                return _peel(l, r, vec, expr)
            # vec is nested deeper in the Cross: rotate the triple product so the
            # factor holding vec comes out (flipping here would flip straight back).
            if r.left.has(vec):
                # Dot(M*a, Cross(f(vec), c)) = Dot(f(vec), Cross(c, M*a))
                return collect(
                    rule(
                        "collect.triple_product_mvmul_left_nested",
                        expr,
                        Dot(r.left, Cross(r.right, l)),
                    ),
                    vec,
                )
            # Dot(M*a, Cross(b, f(vec))) = Dot(f(vec), Cross(M*a, b))
            return collect(
                rule(
                    "collect.triple_product_mvmul_right_nested",
                    expr,
                    Dot(r.right, Cross(l, r.left)),
                ),
                vec,
            )

        # Both sides are MVMul with the same rotation matrix prefix
        # Dot(R*a, R*b) = Dot(a, b) for orthogonal R (SO3)
        case (MVMul(), MVMul()):
            if isinstance(l.left, SO3) and isinstance(r.left, SO3) and l.left == r.left:
                # Strip the common rotation: Dot(R*a, R*b) → Dot(a, b)
                return collect(rule("collect.strip_rotation", expr, Dot(l.right, r.right)), vec)
            if l.has(vec):
                # Dot(M*vec, w) → Dot(vec, Mᵀ*w); deeper, peel the matrix off and let the
                # Cross/MVMul rules expose vec (e.g. Dot(R1*(vec×a), R2*b) with R1 ≠ R2).
                return _peel(l, r, vec, expr)
            # vec is in the right MVMul: flip so it is on the left
            return collect(rule("collect.flip_mvmul_pair", expr, Dot(r, l)), vec)

        # MVMul on one side only
        case (MVMul(), _):
            if l.has(vec):
                # Dot(M*vec, w) → Dot(vec, Mᵀ*w); deeper, peel one matrix layer and let
                # the Cross/MVMul rules expose vec (f(vec) = vec×b gives
                # Dot(vec, Cross(b, Mᵀ*w))).  If f(vec) is a shape no rule can expose
                # (SVMul, VAdd, ...) the recursion stops with the peeled Dot unchanged.
                return _peel(l, r, vec, expr)
            # vec not in MVMul — the other side holds it
            return collect(rule("collect.flip_mvmul_left", expr, Dot(r, l)), vec)

        case (_, MVMul()) if r.has(vec):
            # vec is inside the MVMul: flip so MVMul is on the left, recurse
            return collect(rule("collect.flip_mvmul_right", expr, Dot(r, l)), vec)

        # Cross on one or both sides — scalar triple product
        case (Cross(), _) if l.has(vec):
            # Flip so vec's Cross is on the right, then use triple product
            return collect(rule("collect.flip_cross_left", expr, Dot(r, l)), vec)

        case (_, Cross()):
            if r.left == vec:
                # Dot(a, Cross(vec, c)) = Dot(vec, Cross(c, a))
                return rule("collect.triple_product_left", expr, Dot(vec, Cross(r.right, l)))
            if r.right == vec:
                # Dot(a, Cross(b, vec)) = Dot(vec, Cross(a, b))
                return rule("collect.triple_product_right", expr, Dot(vec, Cross(l, r.left)))
            # vec is nested deeper in the Cross (e.g. Cross(f(vec), b))
            if r.left.has(vec):
                # Dot(a, Cross(f(vec), c)) — use triple product
                # = Dot(f(vec), Cross(c, a))
                return collect(
                    rule(
                        "collect.triple_product_left_nested", expr, Dot(r.left, Cross(r.right, l))
                    ),
                    vec,
                )
            if r.right.has(vec):
                # Dot(a, Cross(b, f(vec))) = Dot(f(vec), Cross(a, b))
                return collect(
                    rule(
                        "collect.triple_product_right_nested", expr, Dot(r.right, Cross(l, r.left))
                    ),
                    vec,
                )
            return expr

        case _:
            return expr
