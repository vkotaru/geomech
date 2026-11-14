import pytest

from geomech.core.base.expressions import (
    S2,
    SO3,
    TS2,
    TSO3,
    I,
    IdentityMatrix,
    Matrix,
    O,
    Scalar,
    SkewSymmMatrix,
    Vector,
    ZeroMatrix,
    getMatrices,
)
from geomech.core.base.types import ExprType
from geomech.core.operations.addition import MAdd
from geomech.core.operations.calculus import Variation
from geomech.core.operations.multiplication import MMMul, MVMul, SMMul


class TestMatrixCreation:
    def test_basic_matrix(self):
        M = Matrix("M")
        assert str(M) == "M"
        assert M.type == ExprType.MATRIX
        assert M.size == (3, 3)
        assert not M.is_constant

    def test_constant_matrix(self):
        J = Matrix("J", attr=["Constant"])
        assert J.is_constant

    def test_symmetric_matrix(self):
        J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
        assert J.is_symmetric
        assert J.is_constant

    def test_zero_matrix(self):
        assert ZeroMatrix.is_zero
        assert ZeroMatrix.is_constant
        assert str(ZeroMatrix) == "0"

    def test_identity_matrix(self):
        assert IdentityMatrix.is_constant
        assert str(IdentityMatrix) == "I"

    def test_aliases(self):
        assert O is ZeroMatrix
        assert I is IdentityMatrix

    def test_get_matrices(self):
        M, N = getMatrices("M N")
        assert str(M) == "M"
        assert str(N) == "N"

    def test_get_matrices_list(self):
        A, B = getMatrices(["A", "B"])
        assert str(A) == "A"
        assert str(B) == "B"


class TestMatrixArithmetic:
    def test_addition(self):
        M, N = getMatrices("M N")
        result = M + N
        assert isinstance(result, MAdd)
        assert result.type == ExprType.MATRIX

    def test_addition_size_mismatch_fails(self):
        A = Matrix("A", size=(3, 3))
        B = Matrix("B", size=(4, 4))
        with pytest.raises(Exception):
            A + B

    def test_scalar_mul(self):
        M = Matrix("M")
        a = Scalar("a")
        result = M * a
        assert isinstance(result, SMMul)

    def test_scalar_left_mul(self):
        M = Matrix("M")
        a = Scalar("a")
        result = a * M
        assert isinstance(result, SMMul)

    def test_scalar_mul_commutativity(self):
        M = Matrix("M")
        a = Scalar("a")
        left = a * M  # ScalarExpr.__mul__
        right = M * a  # MatrixExpr.__mul__
        assert left == right

    def test_int_mul(self):
        M = Matrix("M")
        result = M * 2
        assert isinstance(result, SMMul)

    def test_vector_mul(self):
        M = Matrix("M")
        x = Vector("x")
        result = M * x
        assert isinstance(result, MVMul)
        assert result.type == ExprType.VECTOR

    def test_matrix_mul(self):
        M, N = getMatrices("M N")
        result = M * N
        assert isinstance(result, MMMul)
        assert result.type == ExprType.MATRIX


class TestMatrixOperations:
    def test_delta_variable(self):
        M = Matrix("M")
        d = M.delta()
        assert isinstance(d, Variation)

    def test_delta_constant(self):
        J = Matrix("J", attr=["Constant"])
        d = J.delta()
        # phase 5b: the zero leaf of the type, ZeroMatrix (printed 0; was a matrix "O")
        assert str(d) == "0"
        assert d.is_constant

    def test_diff_variable(self):
        M = Matrix("M")
        dM = M.t_diff()
        assert str(dM) == "\\frac{d}{dt}(M)"

    def test_diff_constant(self):
        J = Matrix("J", attr=["Constant"])
        dJ = J.t_diff()
        assert dJ.is_constant
        assert dJ.is_zero


class TestMatrixIntegrate:
    def test_integrate_undoes_diff(self):
        M = Matrix("M")
        dM = M.t_diff()
        result = dM.t_integrate()
        assert str(result) == "M"

    def test_integrate_adds_prefix(self):
        M = Matrix("M")
        result = M.t_integrate()
        assert str(result) == "\\int{M}dt"


class TestSkewSymmMatrix:
    def test_creation_with_name(self):
        S = SkewSymmMatrix("S")
        assert str(S) == "S"
        assert "SkewSymmetricMatrix" in S.attr


class TestSO3Manifold:
    def test_so3_creation(self):
        R = SO3("R")
        assert R.is_manifold
        assert R.type == ExprType.MATRIX

    def test_so3_variation_vector(self):
        R = SO3("R")
        eta = R.get_variation_vector()
        assert isinstance(eta, Vector)
        assert "\\eta" in str(eta)

    def test_so3_tangent_vector(self):
        R = SO3("R")
        Om = R.get_tangent_vector()
        assert isinstance(Om, TSO3)
        assert "\\Omega" in str(Om)

    def test_so3_delta(self):
        R = SO3("R")
        dR = R.delta()
        assert isinstance(dR, MMMul)

    def test_so3_diff(self):
        R = SO3("R")
        dR = R.t_diff()
        assert isinstance(dR, MMMul)


class TestAttrListIsNotShared:
    """A leaf copies the caller's ``attr`` list; appending its own flags must not leak into
    another leaf built from the same list (review F6, fix plan 2.4)."""

    @pytest.mark.parametrize(
        "make",
        [
            lambda attr: SkewSymmMatrix("K", attr=attr),
            lambda attr: SO3("R", attr=attr),
            lambda attr: S2("q", attr=attr),
            lambda attr: Matrix("M", attr=attr),
            lambda attr: Vector("v", attr=attr),
            lambda attr: Scalar("s", attr=attr),
        ],
        ids=["SkewSymmMatrix", "SO3", "S2", "Matrix", "Vector", "Scalar"],
    )
    def test_shared_attr_list_is_not_mutated(self, make):
        shared = ["Constant"]
        leaf = make(shared)
        assert shared == ["Constant"]
        assert leaf.is_constant
        J = Matrix("J", attr=shared)
        assert not J.is_skew_symmetric and not J.flags.is_orthogonal and not J.is_unit_norm
        leaf.attr.append("SymmetricMatrix")  # mutating the leaf's own list does not leak
        assert shared == ["Constant"]

    def test_skew_and_rotation_from_one_list_keep_their_own_flags(self):
        shared = ["Constant"]
        K, R, q = SkewSymmMatrix("K", attr=shared), SO3("R", attr=shared), S2("q", attr=shared)
        assert K.is_skew_symmetric and not K.flags.is_orthogonal
        assert R.flags.is_orthogonal and not R.is_skew_symmetric
        assert q.is_unit_norm and shared == ["Constant"]

    def test_tangent_vectors_do_not_alias(self):
        a, b = TS2("a"), TSO3("b")
        assert a.attr is not b.attr and a.attr == ["TangentVector"]
