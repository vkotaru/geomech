"""Fixes from the phase 3 review (``specs/reviews/phase3-review.md``) beyond its xfail tests:
``sum()`` over expressions, float cancellation, the transposed scaled row, the hoist of a
scaled row out of an inner product, LaTeX of negated operands, and error messages."""

import re

import pytest

from geomech import (
    AlgebraicError,
    Cross,
    Dot,
    ExpressionMismatchError,
    Matrix,
    Number,
    Scalar,
    SVMul,
    Transpose,
    Vector,
    full_simplify,
    simplify,
    to_latex,
    tree_str,
)

x, y = Scalar("x"), Scalar("y")
a, b, c = Vector("a"), Vector("b"), Vector("c")
K = Matrix("K")


class TestSum:
    def test_zero_plus_is_identity(self):
        for e in (x, a, K):
            assert 0 + e is e
            assert 0.0 + e is e

    def test_sum_of_each_type(self):
        assert str(sum([x, y])) == "(x+y)"
        assert str(sum([a, b, c])) == "(a+b+c)"
        assert str(sum([K, K])) == "(K+K)"

    def test_other_numbers_still_add(self):
        assert str(1 + x) == "((1)+x)"
        with pytest.raises(ExpressionMismatchError, match="Add cannot combine SCALAR"):
            1 + a


class TestSimplify:
    def test_float_cancellation_is_zero(self):
        assert str(full_simplify(x * 0.1 * 3 - x * 0.3)) == "0"
        assert str(full_simplify(a * 0.1 * 3 - a * 0.3)) == "0v"
        assert str(full_simplify(Number(0.1) * 3 - 0.3 + x)) == "x"

    def test_genuine_small_coefficient_kept(self):
        assert str(full_simplify(x * 1e-20 + y)) == "(x(1e-20)+y)"

    def test_transposed_scaled_row_is_scaled_column(self):
        assert str(simplify(Transpose(SVMul(Transpose(a), x)))) == "ax"

    def test_scaled_row_hoisted_from_inner_product(self):
        assert str(full_simplify((a.T() - b.T()) * c - a.T() * c + b.T() * c)) == "0"
        assert str(full_simplify(SVMul(Transpose(a), x) * b - x * (a.T() * b))) == "0"


class TestLatexSigns:
    @pytest.mark.parametrize(
        "build, text",
        [
            (lambda: x * (-a.T() * b), "x \\left(-a^T b\\right)"),
            (lambda: K * (-K).T(), "K \\left(-K^T\\right)"),
            (lambda: -((-K).T()), "-\\left(-K^T\\right)"),
            (lambda: Dot(a, -b), "a \\cdot \\left(-b\\right)"),
            (lambda: Cross(a, -b), "a \\times \\left(-b\\right)"),
            (lambda: -a.T() * b, "-a^T b"),
            (lambda: x * -a, "-x a"),
        ],
    )
    def test_negative_operand_wrapped(self, build, text):
        s = to_latex(build())
        assert s == text
        assert "--" not in s and not re.search(r"[^\s(]\s-", s)


class TestMessages:
    def test_sum_operand_quoted_as_written(self):
        with pytest.raises(ExpressionMismatchError) as ei:
            (Dot(a, a) + Dot(b, a) * x) + a
        assert str(ei.value) == ("Add cannot combine SCALAR '(Dot(a,a)+Dot(b,a)x)' and VECTOR 'a'")

    def test_vector_minus_number(self):
        with pytest.raises(ExpressionMismatchError) as ei:
            a - 1
        assert str(ei.value) == "Sub cannot combine VECTOR 'a' and SCALAR '(1)'"

    def test_long_non_expression_operand_shortened(self):
        with pytest.raises(TypeError) as ei:
            Dot(list(range(500)), a)
        assert len(str(ei.value)) < 200 and "..." in str(ei.value)

    def test_long_constant_name_shortened(self):
        m = Scalar("m" * 200, attr=["Constant"])
        with pytest.raises(AlgebraicError) as ei:
            m.t_integrate()
        assert "m" * 100 not in str(ei.value) and "..." in str(ei.value)


def test_topdown_product_labels_use_a_dot():
    text = tree_str(K * a)
    assert "M·V" in text and "M*V" not in text
