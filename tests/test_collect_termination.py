"""Regression tests from the final review of the scalar-configuration-variable feature.

1. ``collect`` must terminate on ``Dot(M a, N b)`` when the target is nested in a
   form its rules cannot expose (``SVMul``, ``VAdd``): return the expression
   unchanged instead of flipping the two sides forever.
2. ``to_standard_form`` must reject an equation that is bilinear in two
   acceleration variables rather than report a partial mass matrix.
"""

import pytest

from geomech import (
    Dot,
    Mul,
    MVMul,
    Scalar,
    SVMul,
    SystemVariables,
    VAdd,
    Variation,
    getMatrices,
    getVectors,
    to_standard_form,
)
from geomech.core.math.collecting import collect
from geomech.utils.errors import NonLinearError


class TestCollectTerminates:
    def setup_method(self):
        self.A, self.B = getMatrices("A B")
        self.v, self.a, self.b = getVectors(["v", "a", "b"])
        self.s = Scalar("s")

    def test_scalar_times_target_inside_left_mvmul(self):
        """Dot(A (v s), B b) = Dot(v s, Aᵀ B b): the matrix is peeled off the side
        holding v, then no rule can expose v inside the SVMul, so it stops there."""
        e = Dot(MVMul(self.A, SVMul(self.v, self.s)), MVMul(self.B, self.b))
        assert str(collect(e, self.v)) == "Dot(v*s,(A)'*B*b)"

    def test_sum_containing_target_inside_left_mvmul(self):
        """Dot(A (v + a), B b) = Dot(v + a, Aᵀ B b); expand would distribute first
        in the pipeline, so the raw VAdd is left as is after the peel."""
        e = Dot(MVMul(self.A, VAdd(self.v, self.a)), MVMul(self.B, self.b))
        assert str(collect(e, self.v)) == "Dot((v+a),(A)'*B*b)"

    def test_scalar_times_target_against_mvmul(self):
        """Dot(v s, B b): pull would hoist s first; raw input is unchanged."""
        e = Dot(SVMul(self.v, self.s), MVMul(self.B, self.b))
        assert str(collect(e, self.v)) == str(e)

    def test_target_inside_mvmul_still_moves_left(self):
        """Dot(A v, w) → Dot(v, Aᵀ w): the useful case still works."""
        e = Dot(MVMul(self.A, self.v), self.b)
        assert str(collect(e, self.v)) == "Dot(v,(A)'*b)"


class TestStandardFormRejectsBilinearAcceleration:
    def test_bilinear_acceleration_raises(self):
        """c ẍ₁ ẍ₂ has no linear coefficient in either acceleration: must raise,
        never report M = {ẍ₁: c ẍ₂, ẍ₂: c ẍ₁} (which would double-count)."""
        x1, x2, c = Scalar("x1"), Scalar("x2"), Scalar("c", attr=["Constant"])
        a1, a2 = x1.t_diff().t_diff(), x2.t_diff().t_diff()
        eom = {str(Variation(x1)): (Variation(x1), Mul(Mul(c, a1), a2))}
        with pytest.raises(NonLinearError):
            to_standard_form(eom, SystemVariables(scalars=[x1, x2]), [])
