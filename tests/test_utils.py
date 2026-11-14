"""The shared structural helpers of ``geomech.core.transformations.utils``.

``is_leaf`` and ``has_nested_add`` were unused by the engine and are deleted (fix plan
5.5); their tests are replaced by tests of the helpers that took over the duplicated code
(``addends``, ``sum_of``, ``flatten_mul``).
"""

from geomech.core.base.expressions import (
    Matrix,
    Scalar,
    Zero,
    ZeroMatrix,
    ZeroVector,
    getVectors,
)
from geomech.core.base.types import ExprType
from geomech.core.operations.addition import Add, MAdd, VAdd
from geomech.core.operations.geometry import Dot
from geomech.core.operations.multiplication import Mul
from geomech.core.transformations.utils import addends, flatten_mul, sum_of


class TestAddends:
    def test_scalar_sum(self):
        a, b, c = Scalar("a"), Scalar("b"), Scalar("c")
        assert addends(Add(a, b, c)) == [a, b, c]

    def test_vector_and_matrix_sums(self):
        x, y = getVectors(["x", "y"])
        A, B = Matrix("A"), Matrix("B")
        assert addends(VAdd(x, y)) == [x, y]
        assert addends(MAdd(A, B)) == [A, B]

    def test_non_sum_is_one_addend(self):
        a, b = Scalar("a"), Scalar("b")
        assert addends(Mul(a, b)) == [Mul(a, b)]


class TestSumOf:
    def test_round_trip(self):
        a, b = Scalar("a"), Scalar("b")
        x, y = getVectors(["x", "y"])
        assert str(sum_of(addends(Add(a, b)), ExprType.SCALAR)) == str(Add(a, b))
        assert str(sum_of(addends(VAdd(x, y)), ExprType.VECTOR)) == str(VAdd(x, y))

    def test_empty_is_typed_zero(self):
        assert sum_of([], ExprType.SCALAR) is Zero
        assert sum_of([], ExprType.VECTOR) is ZeroVector
        assert sum_of([], ExprType.MATRIX) is ZeroMatrix

    def test_single_term_is_returned_as_is(self):
        x, y = getVectors(["x", "y"])
        d = Dot(x, y)
        assert sum_of([d], ExprType.SCALAR) is d


class TestFlattenMul:
    def test_left_to_right(self):
        a, b, c, d = (Scalar(n) for n in "abcd")
        assert flatten_mul(Mul(Mul(a, b), Mul(c, d))) == [a, b, c, d]

    def test_non_mul(self):
        a = Scalar("a")
        assert flatten_mul(a) == [a]
