"""Operator dispatch (fix plan phase 3, items 3.1-3.3).

3.1  ``Transpose(R) * x`` dispatches on the type of the transpose: a transposed matrix
     multiplies like a matrix; a row vector times a scalar is the row ``(v s)ᵀ``, one
     transpose of a column, never the double transpose ``((v)' s)'``.
3.2  ``Transpose(Transpose(x))`` is ``x`` at construction, and the row/column orientation
     survives ``expand``/``simplify`` (engine fuzz report, finding 6).
3.3  Numbers on either side of ``+``, ``-``, ``*``; unary ``-``; ``/`` by a number.

Each structural assertion is backed by a numeric evaluation of the built tree against the
numpy expression it stands for.
"""

from __future__ import annotations

import numpy as np
import pytest

from geomech import (
    SO3,
    Add,
    Cross,
    Dot,
    ExpressionMismatchError,
    MAdd,
    Matrix,
    MMMul,
    Mul,
    MVMul,
    Number,
    Scalar,
    SMMul,
    SVMul,
    Transpose,
    VAdd,
    Vector,
    VVMul,
)
from geomech.core.operations.geometry import is_row_vector
from geomech.core.transformations import expand, full_simplify, simplify
from geomech.numeric import evaluate

rng = np.random.default_rng(3)
A_VAL = rng.normal(size=(3, 1))
B_VAL = rng.normal(size=(3, 1))
C_VAL = rng.normal(size=(3, 1))
K_VAL = rng.normal(size=(3, 3))
Q_VAL = rng.normal(size=(3, 3))
R_VAL = np.linalg.qr(rng.normal(size=(3, 3)))[0]
S_VAL = 1.7
ENV = {"a": A_VAL, "b": B_VAL, "c": C_VAL, "K": K_VAL, "Q": Q_VAL, "R": R_VAL, "s": S_VAL}


@pytest.fixture
def sym():
    return dict(
        a=Vector("a"),
        b=Vector("b"),
        c=Vector("c"),
        K=Matrix("K"),
        Q=Matrix("Q"),
        R=SO3("R"),
        s=Scalar("s"),
    )


def num(expr):
    return evaluate(expr, ENV)


# ---------------------------------------------------------------------------
# 3.1  Transpose(...) * x
# ---------------------------------------------------------------------------


class TestTransposedMatrixProducts:
    """Rᵀ v, Rᵀ M, Rᵀ s used to raise (engine audit F5, dynamics audit F7)."""

    def test_transposed_matrix_times_vector_is_mvmul(self, sym):
        out = Transpose(sym["R"]) * sym["a"]
        assert isinstance(out, MVMul)
        assert str(out) == str(MVMul(Transpose(sym["R"]), sym["a"])) == "(R)'a"
        np.testing.assert_allclose(num(out), R_VAL.T @ A_VAL)

    def test_transposed_matrix_times_matrix_is_mmmul(self, sym):
        out = Transpose(sym["R"]) * sym["K"]
        assert isinstance(out, MMMul)
        assert str(out) == "(R)'K"
        np.testing.assert_allclose(num(out), R_VAL.T @ K_VAL)

    def test_transposed_matrix_times_scalar_is_smmul(self, sym):
        out = Transpose(sym["R"]) * sym["s"]
        assert isinstance(out, SMMul)
        assert str(out) == "(R)'s"
        np.testing.assert_allclose(num(out), S_VAL * R_VAL.T)

    def test_transposed_matrix_times_number(self, sym):
        out = Transpose(sym["K"]) * 2
        assert isinstance(out, SMMul)
        np.testing.assert_allclose(num(out), 2 * K_VAL.T)

    def test_matrix_times_transposed_matrix(self, sym):
        out = sym["K"] * Transpose(sym["R"])
        assert isinstance(out, MMMul)
        np.testing.assert_allclose(num(out), K_VAL @ R_VAL.T)

    def test_body_frame_gravity_reads_as_written(self, sym):
        """``Rᵀ e₃`` (the most common SO(3) expression) through ``*``."""
        e3 = Vector("e3", attr=["Constant"])
        assert str(Transpose(sym["R"]) * e3) == str(MVMul(Transpose(sym["R"]), e3))


class TestRowVectorProducts:
    def test_row_times_scalar_is_one_transpose_of_a_column(self, sym):
        a, s = sym["a"], sym["s"]
        out = Transpose(a) * s
        assert isinstance(out, Transpose)
        assert isinstance(out.expr, SVMul)
        assert not isinstance(out.expr.left, Transpose)  # no ((a)'s)'
        assert str(out) == "(as)'"
        assert is_row_vector(out)
        np.testing.assert_allclose(num(out), S_VAL * A_VAL.T)

    def test_scalar_times_row_keeps_the_pre_phase3_tree(self, sym):
        """A scalar on the left of a row builds ``SVMul(Transpose(a), s)`` as it did before
        phase 3 (printed ``(a)'s``), so existing trees keep comparing equal (phase 3
        review F1). It is the same row as ``Transpose(a) * s`` numerically."""
        a, s = sym["a"], sym["s"]
        assert str(s * Transpose(a)) == "(a)'s"
        np.testing.assert_allclose(num(s * Transpose(a)), num(Transpose(a) * s))
        np.testing.assert_allclose(num(2 * Transpose(a)), num(Transpose(a) * 2))

    def test_scaled_row_is_still_a_row(self, sym):
        a, b, s = sym["a"], sym["b"], sym["s"]
        out = (Transpose(a) * s) * b
        assert isinstance(out, VVMul)
        assert out.type.name == "SCALAR"
        assert num(out) == pytest.approx(S_VAL * (A_VAL.T @ B_VAL).item())

    def test_row_times_column_is_inner_product(self, sym):
        out = Transpose(sym["a"]) * sym["b"]
        assert isinstance(out, VVMul) and out.type.name == "SCALAR"
        assert num(out) == pytest.approx((A_VAL.T @ B_VAL).item())

    def test_column_times_row_is_outer_product(self, sym):
        out = sym["a"] * Transpose(sym["b"])
        assert isinstance(out, VVMul) and out.type.name == "MATRIX"
        np.testing.assert_allclose(num(out), A_VAL @ B_VAL.T)

    def test_row_times_matrix_is_a_row(self, sym):
        out = Transpose(sym["a"]) * sym["K"]
        assert isinstance(out, MVMul)
        assert is_row_vector(out)
        np.testing.assert_allclose(num(out), A_VAL.T @ K_VAL)

    def test_row_matrix_chain(self, sym):
        """aᵀK Q and aᵀK b: a row built from MVMul (not a Transpose node) keeps working."""
        row = Transpose(sym["a"]) * sym["K"]
        np.testing.assert_allclose(num(row * sym["Q"]), A_VAL.T @ K_VAL @ Q_VAL)
        assert num(row * sym["b"]) == pytest.approx((A_VAL.T @ K_VAL @ B_VAL).item())

    def test_negated_row(self, sym):
        out = -Transpose(sym["a"])
        assert is_row_vector(out)
        np.testing.assert_allclose(num(out), -A_VAL.T)


class TestIsRowVector:
    def test_orientation_of_each_shape(self, sym):
        a, K, s = sym["a"], sym["K"], sym["s"]
        row = MVMul(Transpose(a), K)
        assert not is_row_vector(a)
        assert is_row_vector(Transpose(a))
        assert is_row_vector(row)
        assert not is_row_vector(Transpose(row))  # (aᵀK)ᵀ = Kᵀa is a column
        assert is_row_vector(SVMul(Transpose(a), s))
        assert not is_row_vector(MVMul(K, a))
        assert not is_row_vector(Cross(a, sym["b"]))
        assert is_row_vector(Transpose(a).t_diff())
        assert not is_row_vector(s)


# ---------------------------------------------------------------------------
# 3.2  Transpose(Transpose(x)) collapses at construction
# ---------------------------------------------------------------------------


class TestDoubleTransposeCollapses:
    def test_vector(self, sym):
        a = sym["a"]
        assert Transpose(Transpose(a)) is a

    def test_matrix(self, sym):
        assert Transpose(Transpose(sym["R"])) is sym["R"]

    def test_composite(self, sym):
        inner = MVMul(sym["K"], sym["a"])
        assert Transpose(Transpose(inner)) is inner
        assert Transpose(Transpose(Transpose(inner))) == Transpose(inner)

    def test_method_form(self, sym):
        a = sym["a"]
        assert a.T().T() is a

    def test_delta_and_t_diff_of_a_transpose(self, sym):
        a = sym["a"]
        assert str(Transpose(a).t_diff()) == "(\\frac{d}{dt}(a))'"
        assert str(Transpose(a).delta()) == "(\\delta{a})'"


class TestRowOrientationSurvivesTransformations:
    """Engine fuzz report finding 6: legal trees whose rows/columns ``expand``/``simplify``
    used to confuse (``ExpressionMismatchError``)."""

    def test_expand_outer_product_with_transposed_row_addend(self, sym):
        a, b, c, K = sym["a"], sym["b"], sym["c"], sym["K"]
        col = Transpose(MVMul(Transpose(b), K))  # (bᵀK)ᵀ, a column
        expr = VVMul(VAdd(a, col), Transpose(c))
        out = expand(expr)
        assert isinstance(out, MAdd)
        np.testing.assert_allclose(num(out), (A_VAL + K_VAL.T @ B_VAL) @ C_VAL.T)

    def test_expand_inner_product_with_transposed_row_addend(self, sym):
        a, b, c, K = sym["a"], sym["b"], sym["c"], sym["K"]
        col = Transpose(MVMul(Transpose(b), K))
        expr = VVMul(Transpose(c), VAdd(a, col))
        out = expand(expr)
        assert isinstance(out, Add)
        assert num(out) == pytest.approx((C_VAL.T @ (A_VAL + K_VAL.T @ B_VAL)).item())

    def test_expand_double_transpose_addend(self, sym):
        a, b, c = sym["a"], sym["b"], sym["c"]
        expr = VVMul(VAdd(a, Transpose(Transpose(b))), Transpose(c))
        assert str(expand(expr)) == "(a(c)'+b(c)')"

    @pytest.mark.parametrize("transform", [simplify, full_simplify])
    def test_simplify_row_times_matrix(self, sym, transform):
        a, K, Q = sym["a"], sym["K"], sym["Q"]
        row = MVMul(Transpose(a), K)
        expr = MVMul(Transpose(Transpose(row)), Q)
        out = transform(expr)
        assert str(out) == "(a)'KQ"
        np.testing.assert_allclose(num(out), A_VAL.T @ K_VAL @ Q_VAL)


# ---------------------------------------------------------------------------
# 3.3  Numbers on either side
# ---------------------------------------------------------------------------


class TestNumericOperands:
    def setup_method(self):
        self.x = Scalar("x")
        self.v = Vector("v")
        self.M = Matrix("M")

    def test_rmul(self):
        out = 2 * self.x
        assert isinstance(out, Mul)
        assert str(out) == "(2)x"
        assert str(2 * self.v) == str(self.v * 2) == "v(2)"
        assert str(2 * self.M) == str(self.M * 2) == "M(2)"

    def test_add_and_radd(self):
        assert str(self.x + 1) == "(x+(1))"
        assert str(1 + self.x) == "((1)+x)"
        assert isinstance(1.5 + self.x, Add)

    def test_sub_and_rsub(self):
        assert str(self.x - 1) == "(x+(-1))"
        assert str(1 - self.x) == "((1)+x(-1))"
        assert str(self.x - 2.5) == "(x+(-2.5))"

    def test_neg_builds_the_same_tree_as_times_minus_one(self):
        for e in (self.x, self.v, self.M, Dot(self.v, self.v), Cross(self.v, self.v)):
            assert str(-e) == str(e * -1)
            assert type(-e) is type(e * -1)
        assert str(-self.x) == "x(-1)"
        assert str(-self.v) == "v(-1)"

    def test_neg_matches_existing_pinned_sub(self):
        """``a - b`` has always built ``a + b(-1)``; ``a + -b`` is the same tree."""
        y = Scalar("y")
        assert str(self.x - y) == str(self.x + -y)
        w = Vector("w")
        assert str(self.v - w) == str(self.v + -w)

    def test_mixed_sign_lagrangian(self):
        """``-m*g*Dot(e3, x)`` used to raise TypeError (dynamics audit F7 nit)."""
        m, g = Scalar("m", attr=["Constant"]), Scalar("g", attr=["Constant"])
        e3, x = Vector("e3", attr=["Constant"]), Vector("x")
        assert str(-m * g * Dot(e3, x)) == "m(-1)gDot(e3,x)"

    def test_div_by_number(self):
        assert str(self.x / 2) == str(self.x * 0.5) == "x(0.5)"
        assert str(self.v / 4) == "v(0.25)"
        assert str(self.M / 2) == "M(0.5)"

    def test_div_numeric_value(self):
        env = {"v": A_VAL}
        np.testing.assert_allclose(evaluate(self.v / 4, env), A_VAL / 4)

    def test_div_by_expression_raises(self):
        with pytest.raises(TypeError, match="cannot divide Scalar 'x' by the expression"):
            self.x / Scalar("m")
        with pytest.raises(TypeError, match="only division by a number is supported"):
            2 / self.x

    def test_div_by_zero(self):
        with pytest.raises(ZeroDivisionError, match="cannot divide Scalar 'x' by zero"):
            self.x / 0

    @pytest.mark.parametrize(
        "value, text",
        [(np.int64(2), "(2)"), (np.int32(3), "(3)"), (np.float64(0.5), "(0.5)")],
    )
    def test_numpy_scalars_wrap_like_python_numbers(self, value, text):
        """Follow-up F-3 of the phase 2 review: ``x * np.int64(2)`` raised AttributeError."""
        assert str(self.x * value) == "x" + text
        assert str(value * self.x) == text + "x"
        assert str(self.x + value) == "(x+" + text + ")"
        assert str(Number(value)) == text

    def test_numpy_scalar_times_vector(self):
        assert str(np.float64(2.0) * self.v) == str(self.v * 2.0)

    def test_numpy_int_value_is_a_python_int(self):
        assert type(Number(np.int64(2)).value) is int
        assert type(Number(np.float64(2.5)).value) is float

    @pytest.mark.parametrize("value", [True, False, np.bool_(True)])
    def test_bool_rejected(self, value):
        with pytest.raises(TypeError, match="got bool"):
            self.x * value
        with pytest.raises(TypeError, match="got bool"):
            value + self.x
        with pytest.raises(TypeError, match="got bool"):
            Number(value)

    def test_constructors_wrap_numbers(self):
        assert str(Add(self.x, 1)) == "(x+(1))"
        assert str(Mul(self.x, np.int64(3))) == "x(3)"

    def test_vector_plus_number_names_the_operands(self):
        with pytest.raises(ExpressionMismatchError, match="VAdd cannot combine VECTOR 'v'"):
            self.v + 1

    def test_inplace_add(self):
        y = self.x
        y += 1
        assert str(y) == "(x+(1))"
