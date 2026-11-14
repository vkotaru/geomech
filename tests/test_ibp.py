"""Tests for geomech.core.math.ibp — integration by parts."""

import pytest

from geomech.core.base.expressions import Scalar, Vector
from geomech.core.math.ibp import _apply_ibp, integrate_by_parts
from geomech.core.operations.addition import Add
from geomech.core.operations.calculus import TimeDerivative
from geomech.core.operations.geometry import Dot
from geomech.core.operations.multiplication import Mul


@pytest.fixture
def vectors():
    v = Vector("v")
    w = Vector("w")
    u = Vector("u")
    return v, w, u


@pytest.fixture
def dot_vectors():
    """Time-derivative vectors and their integrals."""
    eta = Vector("eta")
    dot_eta = TimeDerivative(eta)
    return dot_eta, eta


# ===========================================================================
# _apply_ibp — single pass
# ===========================================================================


class TestIbpPassthrough:
    """Expressions that don't match the IBP pattern return unchanged."""

    def test_scalar_leaf(self, vectors):
        """Scalar leaf passes through."""
        v, _, _ = vectors
        s = Scalar("s")
        assert _apply_ibp(s, v) is s

    def test_dot_left_not_target(self, vectors):
        """Dot(w, u) where left != target returns unchanged."""
        v, w, u = vectors
        expr = Dot(w, u)
        result = _apply_ibp(expr, v)
        assert str(result) == str(expr)

    def test_vector_leaf(self, vectors):
        """Vector leaf passes through."""
        v, w, _ = vectors
        assert _apply_ibp(w, v) is w


class TestIbpDot:
    """IBP on Dot where left == target."""

    def test_dot_left_is_target(self, dot_vectors):
        """Dot(d/dt(eta), w) → Dot(-eta, d/dt(w))."""
        dot_eta, eta = dot_vectors
        w = Vector("w")
        expr = Dot(dot_eta, w)
        result = _apply_ibp(expr, dot_eta)
        # Result should be Dot(-integrate(dot_eta), diff(w))
        assert isinstance(result, Dot)
        # Left side: TimeDerivative(eta).t_integrate() * (-1) → SVMul(eta, -1)
        # Right side: w.t_diff() → TimeDerivative(w)
        assert result.left.has(eta)
        assert isinstance(result.right, TimeDerivative)


class TestIbpRecursion:
    """IBP recurses through Add and Mul."""

    def test_add_recurses(self, dot_vectors):
        """Add(Dot(target, w), Dot(target, u)) applies IBP to both terms."""
        dot_eta, eta = dot_vectors
        w = Vector("w")
        u = Vector("u")
        expr = Add(Dot(dot_eta, w), Dot(dot_eta, u))
        result = _apply_ibp(expr, dot_eta)
        assert isinstance(result, Add)
        assert len(result.nodes) == 2
        # Both terms should have been transformed
        assert isinstance(result.nodes[0], Dot)
        assert isinstance(result.nodes[1], Dot)

    def test_mul_recurses_constant_prefactor(self, dot_vectors):
        """Mul(c, Dot(target, w)) with constant c: IBP acts inside the Mul only.

        ∫ c ⟨η̇, w⟩ dt = −∫ c ⟨η, ẇ⟩ dt, so c stays as the untouched left factor.
        """
        dot_eta, eta = dot_vectors
        w = Vector("w")
        c = Scalar("c", attr=["Constant"])
        expr = Mul(c, Dot(dot_eta, w))
        result = _apply_ibp(expr, dot_eta)
        assert isinstance(result, Mul)
        assert str(result.left) == "c"
        assert result.right == Dot(eta * (-1), TimeDerivative(w))

    def test_mul_time_varying_prefactor(self, dot_vectors):
        """Mul(s, Dot(target, w)) with time-varying s: the product rule applies.

        ∫ s ⟨η̇, w⟩ dt = −∫ [ s ⟨η, ẇ⟩ + ⟨η, ṡ w⟩ ] dt. Treating s as constant
        drops the ṡ term, which loses every Coriolis / rate-of-inertia term in
        systems with a scalar configuration variable.
        """
        dot_eta, eta = dot_vectors
        w = Vector("w")
        s = Scalar("s")
        expr = Mul(s, Dot(dot_eta, w))
        result = _apply_ibp(expr, dot_eta)
        assert isinstance(result, Add)
        assert len(result.nodes) == 2
        assert result.nodes[0] == Mul(s, Dot(eta * (-1), TimeDerivative(w)))
        assert result.nodes[1] == Dot(eta * (-1), w * TimeDerivative(s))

    def test_add_mixed_terms(self, dot_vectors):
        """Add with some terms matching target and some not."""
        dot_eta, _ = dot_vectors
        w = Vector("w")
        u = Vector("u")
        # Dot(dot_eta, w) + Dot(u, w)  — only first term gets transformed
        expr = Add(Dot(dot_eta, w), Dot(u, w))
        result = _apply_ibp(expr, dot_eta)
        assert isinstance(result, Add)
        # Second term should be unchanged
        assert str(result.nodes[1].left) == "u"


# ===========================================================================
# integrate_by_parts — full pipeline
# ===========================================================================


class TestIntegrateByParts:
    """Full IBP pipeline: expand → simplify → collect → _apply_ibp → simplify."""

    def test_simple_dot(self, dot_vectors):
        """Simple Dot(d/dt(eta), w) through the full pipeline.

        IBP gives Dot(-eta, d/dt(w)), simplify pulls -1 out → Mul(-1, Dot(eta, d/dt(w))).
        Since phase 5a a numeric factor is placed last: Mul(Dot(eta, d/dt(w)), -1).
        """
        dot_eta, eta = dot_vectors
        w = Vector("w")
        expr = Dot(dot_eta, w)
        result = integrate_by_parts(expr, [dot_eta])
        # full_simplify pulls the -1 scalar out of the Dot
        assert isinstance(result, Mul)
        assert result.right.value == -1
        assert isinstance(result.left, Dot)
        assert str(result.left.left) == "eta"
        assert isinstance(result.left.right, TimeDerivative)

    def test_flipped_dot_gets_collected_then_ibpd(self, dot_vectors):
        """Dot(w, d/dt(eta)) — collect flips it, then IBP applies."""
        dot_eta, eta = dot_vectors
        w = Vector("w")
        expr = Dot(w, dot_eta)
        result = integrate_by_parts(expr, [dot_eta])
        # collect flips to Dot(dot_eta, w), then IBP → Mul(Dot(eta, d/dt(w)), -1)
        # (a numeric factor is placed last since phase 5a)
        assert isinstance(result, Mul)
        assert result.right.value == -1
        assert isinstance(result.left, Dot)
