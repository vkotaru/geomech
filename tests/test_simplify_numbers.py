"""Phase 2 item 2.11: anonymous numerics carry the Zero / Ones flags.

Before, ``Number(0)`` was only ``Constant``: ``Number(0) * x`` simplified to ``(0)x`` with
``is_zero == False`` and a product of numerics that folded to 0 was an unflagged ``(0)``.
"""

from geomech import (
    Matrix,
    Mul,
    Number,
    Scalar,
    SMMul,
    SVMul,
    Vector,
    Zero,
    full_simplify,
)
from geomech.core.operations.multiplication import _wrap_numeric


class TestNumberFlags:
    def test_zero_is_flagged(self):
        for z in (Number(0), Number(0.0), Number(-0.0)):
            assert z.is_zero and z.is_constant and z.is_numeric
            assert not z.is_ones

    def test_one_is_flagged(self):
        for o in (Number(1), Number(1.0)):
            assert o.is_ones and o.is_constant and o.is_numeric
            assert not o.is_zero

    def test_other_values_carry_neither(self):
        for n in (Number(2), Number(-1), Number(0.5)):
            assert n.is_constant and n.is_numeric
            assert not n.is_zero and not n.is_ones

    def test_printed_name_unchanged(self):
        assert str(Number(0)) == "(0)"
        assert str(Number(1)) == "(1)"
        assert str(Number(2.5)) == "(2.5)"

    def test_named_numeric_scalar_is_not_touched(self):
        """D1 (fix plan) is open: Scalar("m", value=0) keeps its own flags."""
        m = Scalar("m", value=0)
        assert m.is_numeric and not m.is_zero and not m.is_constant

    def test_operator_wrapping_uses_the_same_flags(self):
        assert _wrap_numeric(0).is_zero
        assert _wrap_numeric(1).is_ones
        x, v, M = Scalar("x"), Vector("v"), Matrix("M")
        assert (x * 0).right.is_zero
        assert (v * 0).right.is_zero
        assert (M * 1).right.is_ones


class TestZeroProducts:
    def test_number_zero_times_scalar(self):
        x = Scalar("x")
        out = full_simplify(Number(0) * x)
        assert out.is_zero and str(out) == "0"

    def test_number_zero_times_vector(self):
        v = Vector("v")
        out = full_simplify(SVMul(v, Number(0)))
        assert out.is_zero and str(out) == "0v"
        assert full_simplify(v * 0).is_zero

    def test_number_zero_times_matrix(self):
        M = Matrix("M")
        out = full_simplify(SMMul(M, Number(0)))
        assert out.is_zero and str(out) == "0"
        assert full_simplify(M * 0).is_zero

    def test_folded_numeric_product_is_the_zero_leaf(self):
        out = full_simplify(Mul(Number(0), Number(3)))
        assert out is Zero
        out = full_simplify(Mul(Number(3), Number(0)))
        assert out is Zero

    def test_nested_fold_to_zero(self):
        a = Scalar("a")
        assert full_simplify(Mul(Mul(a, Number(2)), Number(0))).is_zero
        assert full_simplify(Mul(Mul(a, Number(2)), Number(0))) is Zero

    def test_nested_fold_to_one_drops_the_factor(self):
        a = Scalar("a")
        assert str(full_simplify(Mul(Mul(a, Number(2)), Number(0.5)))) == "a"

    def test_folded_numeric_sum_is_flagged(self):
        out = full_simplify(Number(2) + Number(-1))
        assert out.is_ones and str(out) == "(1)"
        assert full_simplify(Number(2) + Number(-2)).is_zero

    def test_ones_factor_is_dropped(self):
        a, v = Scalar("a"), Vector("v")
        assert str(full_simplify(Mul(a, Number(1)))) == "a"
        assert str(full_simplify(SVMul(v, Number(1)))) == "v"
