"""``extract_linear_coeff`` never returns a zero or a query-containing coefficient for a
term that depends on the query (fix plan 1.4; engine audit F3/F4, fuzz findings 7, 8, 9).

- A vector query inside a matrix expression (``Hat(v)``, ``M Hat(v)``, ``(v wᵀ s)ᵀ``) is
  linear in the query but has no 3×3 left coefficient: ``AlgebraicError``, not ``0``.
- ``Cross``/``MVMul`` terms with the query in both children are quadratic:
  ``NonLinearError``, not a "coefficient" that still contains the query.
- A matrix query in a scalar/vector expression, or under a scalar factor of an ``SMMul``,
  has no left coefficient either: ``AlgebraicError`` (previously ``NotImplementedError``).
"""

import numpy as np
import pytest

from geomech import (
    Cross,
    Dot,
    Hat,
    Matrix,
    MMMul,
    MVMul,
    Scalar,
    SMMul,
    SVMul,
    Transpose,
    Vector,
    VVMul,
    extract_linear_coeff,
)
from geomech.numeric import evaluate
from geomech.utils.errors import AlgebraicError, NonLinearError


@pytest.fixture
def syms():
    return dict(
        q=Vector("q"),
        a=Vector("a"),
        b=Vector("b"),
        c=Vector("c"),
        M=Matrix("M"),
        Z=Matrix("Z"),
        s=Scalar("s"),
    )


class TestVectorQueryInMatrixExpressionRaises:
    def test_hat_of_query(self, syms):
        with pytest.raises(AlgebraicError, match="no left matrix coefficient"):
            extract_linear_coeff(Hat(syms["q"]), syms["q"])

    def test_matrix_times_hat_of_query(self, syms):
        """paper_context used to list M·Hat(v) as supported; it returned M·0."""
        with pytest.raises(AlgebraicError, match="no left matrix coefficient"):
            extract_linear_coeff(MMMul(syms["M"], Hat(syms["q"])), syms["q"])

    def test_hat_times_matrix(self, syms):
        with pytest.raises(AlgebraicError):
            extract_linear_coeff(MMMul(Hat(syms["q"]), syms["M"]), syms["q"])

    def test_transpose_of_hat(self, syms):
        with pytest.raises(AlgebraicError, match="no left matrix coefficient"):
            extract_linear_coeff(Transpose(Hat(syms["q"])), syms["q"])

    def test_transpose_of_scaled_outer_product(self, syms):
        q, c, s = syms["q"], syms["c"], syms["s"]
        with pytest.raises(AlgebraicError):
            extract_linear_coeff(Transpose(SMMul(VVMul(q, Transpose(c)), s)), q)

    def test_outer_product(self, syms):
        with pytest.raises(AlgebraicError, match="outer"):
            extract_linear_coeff(VVMul(syms["q"], Transpose(syms["a"])), syms["q"])

    def test_query_free_matrix_terms_are_zero(self, syms):
        """Terms without the query still reduce to the zero matrix."""
        M, a, q = syms["M"], syms["a"], syms["q"]
        assert extract_linear_coeff(Hat(a), q).is_zero
        assert extract_linear_coeff(Transpose(M), q).is_zero
        assert extract_linear_coeff(MMMul(M, Hat(a)), q).is_zero


class TestQuadraticVectorTermsRaise:
    def test_mvmul_hat_q_q(self, syms):
        """MVMul(Hat(q), q) = q × q is quadratic; it used to return Hat(q)."""
        q = syms["q"]
        with pytest.raises(NonLinearError, match="both sides"):
            extract_linear_coeff(MVMul(Hat(q), q), q)

    def test_mvmul_outer_product_q(self, syms):
        q, a = syms["q"], syms["a"]
        with pytest.raises(NonLinearError):
            extract_linear_coeff(MVMul(VVMul(q, Transpose(a)), q), q)

    def test_mvmul_hat_q_cross(self, syms):
        q, a = syms["q"], syms["a"]
        with pytest.raises(NonLinearError):
            extract_linear_coeff(MVMul(Hat(q), Cross(a, q)), q)

    def test_mvmul_scaled_matrix_dot_q(self, syms):
        q, a, M = syms["q"], syms["a"], syms["M"]
        with pytest.raises(NonLinearError):
            extract_linear_coeff(MVMul(SMMul(M, Dot(q, a)), q), q)

    def test_cross_cross_a_q_q(self, syms):
        """Cross(Cross(a, q), q) used to return Hat(Cross(a, q)) (contains q)."""
        q, a = syms["q"], syms["a"]
        with pytest.raises(NonLinearError, match="both sides"):
            extract_linear_coeff(Cross(Cross(a, q), q), q)

    def test_cross_q_cross_a_q(self, syms):
        q, a = syms["q"], syms["a"]
        with pytest.raises(NonLinearError):
            extract_linear_coeff(Cross(q, Cross(a, q)), q)

    def test_cross_q_q_is_zero_matrix(self, syms):
        """Cross(q, q) is identically zero: still the zero matrix, not an error."""
        q = syms["q"]
        assert extract_linear_coeff(Cross(q, q), q).is_zero


class TestLinearVectorTermsStillExtract:
    """The guards must not reject linear terms; check the coefficient numerically."""

    @staticmethod
    def _check(expr, q, trials=3):
        rng = np.random.default_rng(0)
        coeff = extract_linear_coeff(expr, q)
        assert not coeff.has(q)
        for _ in range(trials):
            env = {"a": rng.normal(size=(3, 1)), "b": rng.normal(size=(3, 1))}
            env["M"] = rng.normal(size=(3, 3))
            env["s"] = float(rng.normal())
            qv = rng.normal(size=(3, 1))
            lhs = evaluate(expr, {**env, "q": qv})
            rhs = evaluate(coeff, env) @ qv
            assert np.allclose(lhs, rhs)

    def test_cross_a_q(self, syms):
        self._check(Cross(syms["a"], syms["q"]), syms["q"])

    def test_cross_q_b(self, syms):
        self._check(Cross(syms["q"], syms["b"]), syms["q"])

    def test_mvmul_M_q(self, syms):
        self._check(MVMul(syms["M"], syms["q"]), syms["q"])

    def test_mvmul_hat_a_cross_b_q(self, syms):
        a, b, q = syms["a"], syms["b"], syms["q"]
        self._check(MVMul(Hat(a), Cross(b, q)), q)

    def test_mvmul_hat_q_b(self, syms):
        self._check(MVMul(Hat(syms["q"]), syms["b"]), syms["q"])

    def test_cross_M_q_b(self, syms):
        M, q, b = syms["M"], syms["q"], syms["b"]
        self._check(Cross(MVMul(M, q), b), q)

    def test_svmul_q_s(self, syms):
        self._check(SVMul(syms["q"], syms["s"]), syms["q"])


class TestMatrixQueryRaisesAlgebraicError:
    def test_matrix_query_in_scalar_expr(self, syms):
        a, b, M = syms["a"], syms["b"], syms["M"]
        with pytest.raises(AlgebraicError, match="no left coefficient"):
            extract_linear_coeff(Dot(a, MVMul(M, b)), M)

    def test_matrix_query_in_vector_expr(self, syms):
        b, M = syms["b"], syms["M"]
        with pytest.raises(AlgebraicError, match="no left coefficient"):
            extract_linear_coeff(MVMul(M, b), M)

    def test_matrix_query_in_scalar_factor_of_smmul(self, syms):
        """SMMul(Z, <a, M b>) is linear in M but its coefficient is a 4-tensor; used to
        return 0 (fuzz finding 8b)."""
        Z, a, b, M = syms["Z"], syms["a"], syms["b"], syms["M"]
        with pytest.raises(AlgebraicError, match="no left matrix coefficient"):
            extract_linear_coeff(SMMul(Z, Dot(a, MVMul(M, b))), M)

    def test_matrix_query_in_both_smmul_children(self, syms):
        a, b, M = syms["a"], syms["b"], syms["M"]
        with pytest.raises(NonLinearError):
            extract_linear_coeff(SMMul(M, Dot(a, MVMul(M, b))), M)

    def test_matrix_query_smmul_matrix_side_still_works(self, syms):
        Z, M, s = syms["Z"], syms["M"], syms["s"]
        result = extract_linear_coeff(SMMul(MMMul(Z, M), s), M)
        assert str(result) == str(SMMul(Z, s))

    def test_no_notimplemented_from_dispatch(self, syms):
        a, M = syms["a"], syms["M"]
        try:
            extract_linear_coeff(Dot(a, a), M)
        except NotImplementedError:  # pragma: no cover - the old behaviour
            pytest.fail("extract_linear_coeff raised NotImplementedError")
        except AlgebraicError:
            pass
