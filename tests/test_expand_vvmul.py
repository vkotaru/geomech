"""``expand`` distributes ``VVMul`` (outer product ``v wᵀ``, inner product ``vᵀ w``) over
``VAdd`` in either factor, including a ``VAdd`` inside the ``Transpose`` factor.

``Transpose`` itself is still not expanded: ``Transpose(MAdd)``, ``Transpose(VVMul)`` and
``Transpose(Hat(VAdd))`` must come out of ``expand`` unchanged (they are opaque to
``extract_linear_coeff`` and raise there rather than being silently mis-handled).
"""

from geomech.core.base.expressions import Vector, getVectors
from geomech.core.base.types import ExprType
from geomech.core.operations.addition import Add, MAdd, VAdd
from geomech.core.operations.geometry import Cross, Hat, Transpose
from geomech.core.operations.multiplication import MMMul, MVMul, VVMul
from geomech.core.transformations.expand import expand

a, b, c, d = getVectors("a b c d")


def _addends(expr, cls):
    assert isinstance(expr, cls), f"expected {cls.__name__}, got {type(expr).__name__}: {expr}"
    return [str(n) for n in expr.nodes]


class TestOuterProduct:
    def test_left_vadd(self):
        """(a + b) cᵀ → a cᵀ + b cᵀ."""
        out = expand(VVMul(a + b, Transpose(c)))
        assert _addends(out, MAdd) == [str(VVMul(a, Transpose(c))), str(VVMul(b, Transpose(c)))]
        assert out.type == ExprType.MATRIX

    def test_right_transpose_vadd(self):
        """a (b + c)ᵀ → a bᵀ + a cᵀ."""
        out = expand(VVMul(a, Transpose(b + c)))
        assert _addends(out, MAdd) == [str(VVMul(a, Transpose(b))), str(VVMul(a, Transpose(c)))]

    def test_both_sides(self):
        """(a + b)(c + d)ᵀ → four outer products in cartesian order."""
        out = expand(VVMul(a + b, Transpose(c + d)))
        assert _addends(out, MAdd) == [
            str(VVMul(a, Transpose(c))),
            str(VVMul(a, Transpose(d))),
            str(VVMul(b, Transpose(c))),
            str(VVMul(b, Transpose(d))),
        ]

    def test_nested_sum_in_left_factor(self):
        """Cross(a, b + c) dᵀ: the left factor expands to a VAdd first, then distributes."""
        out = expand(VVMul(Cross(a, b + c), Transpose(d)))
        assert _addends(out, MAdd) == [
            str(VVMul(Cross(a, b), Transpose(d))),
            str(VVMul(Cross(a, c), Transpose(d))),
        ]

    def test_nested_sum_inside_transpose(self):
        """a Cross(b, c + d)ᵀ: the inner of the Transpose is expanded and split."""
        out = expand(VVMul(a, Transpose(Cross(b, c + d))))
        assert _addends(out, MAdd) == [
            str(VVMul(a, Transpose(Cross(b, c)))),
            str(VVMul(a, Transpose(Cross(b, d)))),
        ]

    def test_no_sum_is_unchanged(self):
        E = VVMul(a, Transpose(b))
        assert str(expand(E)) == str(E)
        assert isinstance(expand(E), VVMul)

    def test_distributes_under_an_enclosing_product(self):
        """((a + b) cᵀ) d → (a cᵀ) d + (b cᵀ) d: an outer product whose left factor is a sum,
        under an enclosing matrix-vector product."""
        out = expand(MVMul(VVMul(a + b, Transpose(c)), d))
        assert _addends(out, VAdd) == [
            str(MVMul(VVMul(a, Transpose(c)), d)),
            str(MVMul(VVMul(b, Transpose(c)), d)),
        ]

    def test_outer_sum_then_matrix_product(self):
        """((a + b) cᵀ) Hat(d) → (a cᵀ) Hat(d) + (b cᵀ) Hat(d) via the MMMul rule."""
        out = expand(MMMul(VVMul(a + b, Transpose(c)), Hat(d)))
        assert _addends(out, MAdd) == [
            str(MMMul(VVMul(a, Transpose(c)), Hat(d))),
            str(MMMul(VVMul(b, Transpose(c)), Hat(d))),
        ]


class TestInnerProduct:
    def test_left_transpose_vadd(self):
        """(a + b)ᵀ c → aᵀ c + bᵀ c (scalar)."""
        out = expand(VVMul(Transpose(a + b), c))
        assert _addends(out, Add) == [str(VVMul(Transpose(a), c)), str(VVMul(Transpose(b), c))]
        assert out.type == ExprType.SCALAR

    def test_right_vadd(self):
        """aᵀ (b + c) → aᵀ b + aᵀ c."""
        out = expand(VVMul(Transpose(a), b + c))
        assert _addends(out, Add) == [str(VVMul(Transpose(a), b)), str(VVMul(Transpose(a), c))]

    def test_both_sides(self):
        out = expand(VVMul(Transpose(a + b), c + d))
        assert _addends(out, Add) == [
            str(VVMul(Transpose(a), c)),
            str(VVMul(Transpose(a), d)),
            str(VVMul(Transpose(b), c)),
            str(VVMul(Transpose(b), d)),
        ]

    def test_no_sum_is_unchanged(self):
        E = VVMul(Transpose(a), b)
        assert str(expand(E)) == str(E)
        assert isinstance(expand(E), VVMul)


class TestTransposeStaysOpaque:
    """The rule looks through ``Transpose(VAdd)`` only as a VVMul factor."""

    def test_transpose_of_hat_of_sum_is_not_expanded(self):
        E = MVMul(Transpose(Hat(a + b)), c)
        assert str(expand(E)) == str(E)

    def test_transpose_of_madd_is_not_expanded(self):
        E = MVMul(Transpose(Hat(a) + Hat(b)), c)
        assert str(expand(E)) == str(E)

    def test_transpose_of_outer_product_is_not_expanded(self):
        E = MVMul(Transpose(VVMul(a, Transpose(b))), c)
        assert str(expand(E)) == str(E)

    def test_transpose_factor_without_a_sum_keeps_its_original_node(self):
        """A Transpose(Cross(b, c)) factor has no sum: the very same node is kept."""
        t = Transpose(Cross(b, c))
        out = expand(VVMul(a, t))
        assert out.right is t


class TestIdempotent:
    def test_second_pass_is_a_fixed_point(self):
        E = MVMul(VVMul(Cross(a, b + c), Transpose(c + d)), a + b)
        once = expand(E)
        assert str(expand(once)) == str(once)

    def test_terminates_on_deep_nesting(self):
        v = Vector("v")
        E = VVMul(v, Transpose(v))
        for w in (a, b, c, d):
            E = VVMul(MVMul(E, w + v), Transpose(w + v))
        out = expand(E)
        assert isinstance(out, MAdd)
        assert str(expand(out)) == str(out)
