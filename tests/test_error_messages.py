"""Constructor errors name their operands (fix plan 3.4; debuggability audit A3-A9, A16,
A29-A33, F5).

Before phase 3 every tree-building error named only the two ``ExprType``s
(``Dot cannot be performed between SCALAR and VECTOR``), which is the same message for
every ``Dot`` of a long Lagrangian; ``J * x.T()`` raised a bare
``operation cannot be performed``; ``m.t_integrate()`` raised an empty
``NotImplementedError``; ``m * "x"`` and ``m + 1`` raised ``AttributeError`` from inside the
package. The error classes and their bases are unchanged, so ``pytest.raises`` on them
still matches.
"""

from __future__ import annotations

import re

import pytest

import geomech.utils.errors as errors
from geomech import (
    S2,
    SO3,
    Add,
    AlgebraicError,
    Cross,
    Dot,
    ExpressionMismatchError,
    Hat,
    MAdd,
    Matrix,
    MMMul,
    Mul,
    MVMul,
    Scalar,
    SMMul,
    SVMul,
    Transpose,
    VAdd,
    Vector,
    Vee,
    VVMul,
)


@pytest.fixture
def s():
    return dict(
        m=Scalar("m"),
        x=Vector("x"),
        y=Vector("y"),
        J=Matrix("J"),
        mc=Scalar("m", attr=["Constant"]),
    )


def message(fn):
    with pytest.raises(Exception) as info:
        fn()
    return type(info.value), str(info.value)


class TestBinaryConstructorsNameTheOperands:
    @pytest.mark.parametrize(
        "build, text",
        [
            (lambda s: Dot(s["m"], s["x"]), "Dot cannot combine SCALAR 'm' and VECTOR 'x'"),
            (lambda s: Cross(s["x"], s["m"]), "Cross cannot combine VECTOR 'x' and SCALAR 'm'"),
            (lambda s: Mul(s["m"], s["x"]), "Mul cannot combine SCALAR 'm' and VECTOR 'x'"),
            (lambda s: SVMul(s["x"], s["y"]), "SVMul cannot combine VECTOR 'x' and VECTOR 'y'"),
            (lambda s: SMMul(s["J"], s["x"]), "SMMul cannot combine MATRIX 'J' and VECTOR 'x'"),
            (lambda s: MVMul(s["J"], s["m"]), "MVMul cannot combine MATRIX 'J' and SCALAR 'm'"),
            (lambda s: MMMul(s["J"], s["x"]), "MMMul cannot combine MATRIX 'J' and VECTOR 'x'"),
            (lambda s: VVMul(s["J"], s["x"]), "VVMul cannot combine MATRIX 'J' and VECTOR 'x'"),
        ],
    )
    def test_message(self, s, build, text):
        kind, msg = message(lambda: build(s))
        assert kind is ExpressionMismatchError
        assert msg == text

    def test_operand_in_a_long_lagrangian_is_identifiable(self, s):
        """A3: the message tells which ``Dot`` of a Lagrangian is wrong."""
        q = S2("q")
        om = q.get_tangent_vector()
        kind, msg = message(lambda: Dot(Cross(om, q), s["m"]))
        assert msg == "Dot cannot combine VECTOR 'Cross(\\omega_{q},q)' and SCALAR 'm'"

    def test_long_operands_are_shortened(self):
        long = Vector("x" * 200)
        _, msg = message(lambda: Dot(long, Scalar("m")))
        assert "..." in msg
        assert len(msg) < 140


class TestOperatorMessages:
    def test_add_scalar_and_vector(self, s):
        """A4"""
        assert message(lambda: s["m"] + s["x"]) == (
            ExpressionMismatchError,
            "Add cannot combine SCALAR 'm' and VECTOR 'x'",
        )

    def test_sub_vector_and_scalar(self, s):
        """A8"""
        assert message(lambda: s["x"] - s["m"]) == (
            ExpressionMismatchError,
            "Sub cannot combine VECTOR 'x' and SCALAR 'm'",
        )

    def test_matrix_plus_vector(self, s):
        assert message(lambda: s["J"] + s["x"])[1] == (
            "MAdd cannot combine MATRIX 'J' and VECTOR 'x'"
        )

    def test_sum_whose_first_operand_has_the_wrong_type(self, s):
        assert message(lambda: VAdd(s["m"], s["x"]))[1] == (
            "VAdd expects VECTOR operands, got SCALAR 'm'"
        )
        assert message(lambda: MAdd(s["x"]))[1] == "MAdd expects MATRIX operands, got VECTOR 'x'"

    def test_column_times_matrix(self, s):
        """A5: names the operands and says what to write."""
        kind, msg = message(lambda: s["x"] * s["J"])
        assert kind is ExpressionMismatchError
        assert msg == (
            "MVMul: a column vector times a matrix is not defined: 'x' * 'J'; write M * v for"
            " the column M v, or v.T() * M for the row vᵀ M"
        )

    def test_matrix_times_row(self, s):
        """A6: was the bare ``operation cannot be performed``."""
        kind, msg = message(lambda: s["J"] * s["x"].T())
        assert kind is ExpressionMismatchError
        assert msg == (
            "MVMul: a matrix times a row vector is not defined: 'J' * '(x)''; a matrix"
            " multiplies a column vector (M * v), and a row vector multiplies a matrix from"
            " the left (v.T() * M)"
        )
        # the same through the constructor
        assert message(lambda: MVMul(s["J"], Transpose(s["x"])))[1] == msg

    def test_two_rows(self, s):
        """A7: said "VECTOR and VECTOR", which reads as if two vectors could not multiply."""
        kind, msg = message(lambda: s["x"].T() * s["y"].T())
        assert kind is ExpressionMismatchError
        assert msg.startswith("VVMul: two row vectors cannot be multiplied: '(x)'' * '(y)''")
        assert "a.T() * b for the inner product" in msg

    def test_two_columns(self, s):
        _, msg = message(lambda: s["x"] * s["y"])
        assert msg.startswith("VVMul: two column vectors cannot be multiplied: 'x' * 'y'")
        assert "Dot(a, b)" in msg and "Cross(a, b)" in msg


class TestUnaryConstructors:
    def test_hat(self, s):
        """A9: was "Hat cannot be performed between VECTOR and SCALAR"."""
        assert message(lambda: Hat(s["m"])) == (
            ExpressionMismatchError,
            "Hat expects a VECTOR, got SCALAR 'm'",
        )

    def test_vee(self, s):
        assert message(lambda: Vee(s["x"]))[1] == "Vee expects a MATRIX, got VECTOR 'x'"


class TestNonExpressionOperands:
    def test_string_operand(self, s):
        """A29: was ``AttributeError: 'str' object has no attribute 'type'``."""
        assert message(lambda: s["m"] * "x") == (
            TypeError,
            "cannot apply '*' to Scalar 'm' and str 'x'; operands must be geomech"
            " expressions or numbers (int, float)",
        )

    def test_none_in_a_sum(self, s):
        kind, msg = message(lambda: Add(s["m"], None))
        assert kind is TypeError
        assert msg == "Add: operands must be geomech expressions or numbers, got NoneType None"

    def test_number_in_dot(self, s):
        assert message(lambda: Dot(s["x"], 2)) == (
            TypeError,
            "Dot: operands must be geomech expressions, got int 2",
        )

    def test_list_operand(self, s):
        kind, msg = message(lambda: s["x"] + [1, 2, 3])
        assert kind is TypeError
        assert "Vector 'x' and list [1, 2, 3]" in msg

    def test_division_by_expression(self, s):
        """A33"""
        kind, msg = message(lambda: s["x"] / s["m"])
        assert kind is TypeError
        assert msg.startswith("cannot divide Vector 'x' by the expression Scalar 'm'")


class TestConstantIntegral:
    """A16: three empty ``NotImplementedError``s became ``AlgebraicError`` with the symbol."""

    @pytest.mark.parametrize(
        "leaf, kind",
        [
            (Scalar("m", attr=["Constant"]), "scalar"),
            (Vector("e3", attr=["Constant"]), "vector"),
            (Matrix("J", attr=["Constant"]), "matrix"),
        ],
    )
    def test_message(self, leaf, kind):
        with pytest.raises(AlgebraicError) as info:
            leaf.t_integrate()
        name = str(leaf)
        assert str(info.value) == (
            f"t_integrate: cannot integrate the constant {kind} '{name}' in time:"
            f" ∫ {name} dt = {name}·t, and the tree has no time symbol. Only terms that are"
            " time derivatives (d/dt(x)) integrate to an expression"
        )

    def test_non_constant_still_integrates(self):
        x = Scalar("x")
        assert str(x.t_integrate()) == "\\int{x}dt"
        assert x.t_diff().t_integrate() is x


class TestErrorClasses:
    def test_classes_and_bases_unchanged(self):
        assert issubclass(ExpressionMismatchError, Exception)
        assert issubclass(AlgebraicError, Exception)
        assert not issubclass(AlgebraicError, NotImplementedError)

    def test_old_type_only_form_still_available(self):
        from geomech.core.base.types import ExprType

        err = ExpressionMismatchError("Foo", ExprType.SCALAR, ExprType.VECTOR)
        assert str(err) == "Foo cannot be performed between SCALAR and VECTOR"

    def test_undefined_case_error_is_gone(self):
        """It was unreachable (every node type is SCALAR, VECTOR or MATRIX) and carried no
        context."""
        assert not hasattr(errors, "UndefinedCaseError")

    def test_no_bare_raise_left(self):
        """Every ``raise ExpressionMismatchError`` in the package passes arguments."""
        import pathlib

        import geomech

        root = pathlib.Path(geomech.__file__).parent
        bare = []
        for path in root.rglob("*.py"):
            for i, line in enumerate(path.read_text().splitlines(), 1):
                if re.search(r"raise (ExpressionMismatchError|NotImplementedError)\s*$", line):
                    bare.append(f"{path.relative_to(root)}:{i}")
        # Expr.__str__ / get_tangent_vector are abstract hooks, not error paths
        assert bare == [
            "core/base/expressions.py:" + str(_line("def __str__")),
            "core/base/expressions.py:" + str(_line("def get_tangent_vector")),
        ], bare


def _line(marker):
    """Line number of the ``raise`` right after the first *marker* in expressions.py."""
    import pathlib

    import geomech

    text = (pathlib.Path(geomech.__file__).parent / "core/base/expressions.py").read_text()
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if marker in line)
    return next(i for i in range(start, len(lines)) if "raise" in lines[i]) + 1


def test_so3_transpose_dispatch_errors_name_operands():
    R = SO3("R")
    with pytest.raises(ExpressionMismatchError, match="Cross cannot combine MATRIX '\\(R\\)''"):
        Cross(Transpose(R), Vector("x"))
