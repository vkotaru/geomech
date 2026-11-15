"""Pinned equations for the ``collect`` fix (fix plan 1.2).

``collect`` moves the integration-by-parts target to the left of every ``Dot``. Two
defects were fixed:

- the rewrite ``(M v)·(b × c) → v·((M b) × c)`` is a false identity (engine audit F1,
  dynamics audit F2); the only valid move is ``(M v)·w = v·(Mᵀ w)``;
- a target nested deeper than one level under a matrix (``Dot(M f(v), w)``) was left in
  place, so the term survived integration by parts and extraction raised
  ``AlgebraicError`` on a legal Lagrangian (engine audit F6, dynamics audit F3).

Both are handled by peeling the matrix off the side that holds the target,
``Dot(M f(v), w) → Dot(f(v), Mᵀ w)``, one layer per step.

Hand derivations use the body-frame form on SO(3): with p = ∂L/∂Ω, δΩ = η̇ + Ω × η and
δR = R η̂, the coefficient of η is E = −ṗ + p × Ω + (q-terms) + M, where the q-terms come from
δ(Rᵀe₃) = (Rᵀe₃) × η and d/dt(Rᵀe₃) = (Rᵀe₃) × Ω.  Every pinned equation is also checked
with the finite-difference Euler-Lagrange check (``geomech.numeric.assert_eom``).
"""

import numpy as np
import pytest

from geomech import Cross, Dot, Hat, Matrix, MMMul, MVMul, Transpose, Vector, VVMul, compute_eom
from geomech.core.math.collecting import collect
from geomech.numeric import evaluate
from geomech.utils.errors import AlgebraicError, NonLinearError
from tests import phase1_systems as ps
from tests.phase1_systems import derive, equation


def _agree_numerically(eq1, eq2, trials=5):
    for seed in range(trials):
        a, b = evaluate(eq1, seed=seed), evaluate(eq2, seed=seed)
        assert np.max(np.abs(a - b)) < 1e-9


class TestGyroscopicCouplingSO3:
    """L = ½ Ω·JΩ + (JΩ)·(r × b), r = Rᵀe₃ (dynamics audit F2).

    p = JΩ + J(r × b);  ṗ = JΩ̇ + J((r × Ω) × b);  p × Ω = −Ω × JΩ − Ω × J(r × b);
    δR-term: (JΩ)·((r × η) × b) = −η·(r × (b × JΩ)).
    E = −JΩ̇ − Ω × JΩ + J((Ω × r) × b) − Ω × J(r × b) − r × (b × JΩ) + M.
    Before the fix the third term was (J(Ω × r)) × b (the false identity), rel. error 0.7.
    """

    EXPECTED = (
        "(J*\\frac{d}{dt}(\\Omega_{R})*(-1)+Cross(\\Omega_{R},J*\\Omega_{R})*(-1)"
        "+J*Cross(Cross(\\Omega_{R},(R)'*e3),b)+Cross(\\Omega_{R},J*Cross((R)'*e3,b))*(-1)"
        "+Cross((R)'*e3,Cross(b,J*\\Omega_{R}))*(-1)+M)"
    )

    def test_J_on_velocity(self):
        s, eom, _ = derive(ps.so3_gyro_coupling_J_on_velocity)
        assert str(equation(eom, s.symbols["R"])) == self.EXPECTED

    def test_J_on_axis_is_identical(self):
        """½ Ω·JΩ + Ω·J((Rᵀe₃) × b) is the same Lagrangian (J symmetric); same equation."""
        s1, eom1, _ = derive(ps.so3_gyro_coupling_J_on_velocity)
        s2, eom2, _ = derive(ps.so3_gyro_coupling_J_on_axis)
        eq1, eq2 = equation(eom1, s1.symbols["R"]), equation(eom2, s2.symbols["R"])
        assert str(eq2) == self.EXPECTED
        _agree_numerically(eq1, eq2)


class TestInertialAngularMomentumSO3:
    """L = ½ Ω·JΩ + (R(JΩ))·e₃ = ½ Ω·JΩ + Ω·(J r), r = Rᵀe₃ (dynamics audit F3, row 1).

    p = JΩ + J r;  ṗ = JΩ̇ + J(r × Ω);  p × Ω = −Ω × JΩ − Ω × J r;
    δR-term: (JΩ)·(r × η) = η·((JΩ) × r).
    E = −JΩ̇ − Ω × JΩ + (JΩ) × r + J(Ω × r) − Ω × J r + M.
    Both groupings crashed with AlgebraicError before the fix.
    """

    EXPECTED = (
        "(J*\\frac{d}{dt}(\\Omega_{R})*(-1)+Cross(\\Omega_{R},J*\\Omega_{R})*(-1)"
        "+Cross(J*\\Omega_{R},(R)'*e3)+J*Cross(\\Omega_{R},(R)'*e3)+Cross(\\Omega_{R},J*(R)'*e3)*(-1)"
        "+M)"
    )

    def test_R_times_J_Omega(self):
        s, eom, _ = derive(ps.so3_inertial_momentum_R_J_Omega)
        assert str(equation(eom, s.symbols["R"])) == self.EXPECTED

    def test_RJ_times_Omega(self):
        s1, eom1, _ = derive(ps.so3_inertial_momentum_R_J_Omega)
        s2, eom2, _ = derive(ps.so3_inertial_momentum_RJ_Omega)
        eq1, eq2 = equation(eom1, s1.symbols["R"]), equation(eom2, s2.symbols["R"])
        assert str(eq2) == self.EXPECTED
        _agree_numerically(eq1, eq2)


class TestAnisotropicInertiaQdotForm:
    def test_pinnedequation(self):
        """L = ½ q̇·A q̇ − m g q·e₃ with q̇ = ω × q (dynamics audit F3, row 2).

        In ω-form, L = ½ (ω × q)·A(ω × q) − m g q·e₃, so p = ∂L/∂ω = q × A(ω × q),
        ṗ = (ω × q) × A(ω × q) + q × A(ω̇ × q) + q × A(ω × (ω × q)),
        2 ω × p = 2 ω × (q × A(ω × q)), q × ∂L/∂q = −q × (ω × A(ω × q)) − m g q × e₃.
        E = −(ω × q) × A(ω × q) − q × A(ω̇ × q) − q × A(ω × (ω × q))
            + 2 ω × (q × A(ω × q)) − q × (ω × A(ω × q)) − m g q × e₃ + f.
        Crashed with AlgebraicError before the fix (target two levels under A).
        The new printed form is the same expression simplified further: the third term
        prints as q × (A q)(ω·ω), since ω × (ω × q) = −q (ω·ω) for ω ⊥ q (phase 5a).
        """
        s, eom, _ = derive(ps.s2_anisotropic_qdot)
        expected = (
            "(Cross(Cross(\\omega_{q},q),A*Cross(\\omega_{q},q))*(-1)"
            "+Cross(q,A*Cross(\\frac{d}{dt}(\\omega_{q}),q))*(-1)"
            "+Cross(q,A*q)*Dot(\\omega_{q},\\omega_{q})"
            "+Cross(\\omega_{q},Cross(q,A*Cross(\\omega_{q},q)))*(2)"
            "+Cross(q,Cross(\\omega_{q},A*Cross(\\omega_{q},q)))*(-1)+Cross(q,e3)*m*g*(-1)+f)"
        )
        assert str(equation(eom, s.symbols["q"])) == expected


class TestRelativeRateCouplingSO3xSO3:
    def test_pinned_equations(self):
        """L = ½ Ω₁·J₁Ω₁ + ½ Ω₂·J₂Ω₂ + k Ω₁·(QΩ₂), Q = R₁ᵀR₂ (dynamics audit F3, row 4).

        Q̇ = −Ω̂₁Q + QΩ̂₂; δQ = −η̂₁Q + Qη̂₂.
        η₁: p₁ = J₁Ω₁ + kQΩ₂, ṗ₁ = J₁Ω̇₁ − kΩ₁ × QΩ₂ + kQ(Ω₂ × Ω₂) + kQΩ̇₂,
            p₁ × Ω₁ = −Ω₁ × J₁Ω₁ + k(QΩ₂) × Ω₁, δR₁-term k η₁·(Ω₁ × QΩ₂):
            E₁ = −J₁Ω̇₁ − Ω₁ × J₁Ω₁ + kΩ₁ × QΩ₂ − kQ(Ω₂ × Ω₂) − kQΩ̇₂ − kΩ₁ × QΩ₂
                 − k(QΩ₂) × Ω₁ + M₁   (the three k-cross terms sum to +kΩ₁ × QΩ₂).
        η₂: p₂ = J₂Ω₂ + kQᵀΩ₁, ṗ₂ = J₂Ω̇₂ − kΩ₂ × QᵀΩ₁ + kQᵀ(Ω₁ × Ω₁) + kQᵀΩ̇₁,
            p₂ × Ω₂ = −Ω₂ × J₂Ω₂ − kΩ₂ × QᵀΩ₁, δR₂-term k η₂·(Ω₂ × QᵀΩ₁):
            E₂ = −J₂Ω̇₂ − Ω₂ × J₂Ω₂ + 2kΩ₂ × QᵀΩ₁ − kQᵀ(Ω₁ × Ω₁) − kQᵀΩ̇₁ − kΩ₂ × QᵀΩ₁ + M₂.
        The engine keeps the identically-zero Q(Ω × Ω) terms (no rule descends into MVMul).
        The new printed form is the same expression simplified further (phase 5a): the
        Q(Ω × Ω) terms vanish and the k-cross terms are merged.
        Crashed with AlgebraicError before the fix.
        """
        s, eom, _ = derive(ps.so3_so3_relative_rate_coupling)
        assert str(equation(eom, s.symbols["R1"])) == (
            "(J1*\\frac{d}{dt}(\\Omega_{R1})*(-1)+Cross(\\Omega_{R1},J1*\\Omega_{R1})*(-1)"
            "+(R1)'*R2*\\frac{d}{dt}(\\Omega_{R2})*k*(-1)"
            "+Cross((R1)'*R2*\\Omega_{R2},\\Omega_{R1})*k*(-1)+M1)"
        )
        assert str(equation(eom, s.symbols["R2"])) == (
            "(J2*\\frac{d}{dt}(\\Omega_{R2})*(-1)+Cross(\\Omega_{R2},J2*\\Omega_{R2})*(-1)"
            "+Cross(\\Omega_{R2},(R2)'*R1*\\Omega_{R1})*k+(R2)'*R1*\\frac{d}{dt}(\\Omega_{R1})*k*(-1)"
            "+M2)"
        )


class TestOffsetMassMatrixInertia:
    def test_pinnedequation(self):
        """L = ½ (Ω × ρ)·J(Ω × ρ) (engine audit F6).

        p = ρ × J(Ω × ρ) (∂(Ω × ρ)/∂Ω = −ρ̂, (−ρ̂)ᵀ = ρ̂), ṗ = ρ × J(Ω̇ × ρ),
        p × Ω = −Ω × (ρ × J(Ω × ρ)).  E = −ρ × J(Ω̇ × ρ) − Ω × (ρ × J(Ω × ρ)) + M.
        Crashed with AlgebraicError before the fix (the m·I version worked).
        """
        s, eom, _ = derive(ps.so3_offset_mass_matrix_inertia)
        assert str(equation(eom, s.symbols["R"])) == (
            "(Cross(\\rho,J*Cross(\\frac{d}{dt}(\\Omega_{R}),\\rho))*(-1)"
            "+Cross(\\Omega_{R},Cross(\\rho,J*Cross(\\Omega_{R},\\rho)))*(-1)+M)"
        )


class TestTargetOnBothSides:
    """fuzz finding 5: a target on both sides of a Dot looped until RecursionError."""

    def test_cross_both_sides_raises(self):
        a, b, v = Vector("a"), Vector("b"), Vector("v")
        with pytest.raises(NonLinearError, match="both sides"):
            collect(Dot(Cross(a, v), Cross(b, v)), v)

    def test_nested_both_sides_raises(self):
        from geomech import VAdd

        a, b, v = Vector("a"), Vector("b"), Vector("v")
        expr = Dot(a, Cross(Cross(b, VAdd(v, a)), Cross(a, VAdd(v, b))))
        with pytest.raises(NonLinearError):
            collect(expr, v)


class TestTargetInsideMatrixFactor:
    """Review F1: a target inside the matrix factor of an ``MVMul`` looped to
    ``RecursionError`` (the peel moved the matrix across, the flip moved it back and each
    round added a ``Transpose``). Transposition can only expose a target in the vector
    operand, ``(M f(v))·w = f(v)·(Mᵀ w)``, so ``collect`` now raises ``AlgebraicError``
    naming the term and the target, and the pipeline reports it instead of hanging.
    """

    def test_collect_raises_naming_term_and_target(self):
        v, a, b, c, M = Vector("v"), Vector("a"), Vector("b"), Vector("c"), Matrix("M")
        for expr, factor in (
            (Dot(MVMul(Hat(v), a), b), Hat(v)),
            (Dot(b, MVMul(Hat(v), a)), Hat(v)),  # flipped first, then named as Dot(Hat(v)a,b)
            (Dot(MVMul(VVMul(v, Transpose(c)), a), b), VVMul(v, Transpose(c))),
            (Dot(MVMul(MMMul(M, Hat(v)), a), b), MMMul(M, Hat(v))),
            (Dot(MVMul(Hat(v), a), MVMul(M, b)), Hat(v)),
        ):
            with pytest.raises(AlgebraicError, match="matrix factor") as info:
                collect(expr, v)
            msg = str(info.value)
            assert "the target 'v'" in msg and f"matrix factor '{factor}'" in msg
            assert "the term 'Dot(" in msg

    def test_target_in_both_factors_is_nonlinear(self):
        v, b = Vector("v"), Vector("b")
        with pytest.raises(NonLinearError, match="both factors"):
            collect(Dot(MVMul(Hat(v), v), b), v)

    def test_pipeline_so3_outer_product_on_the_velocity(self):
        """L = ½ Ω·JΩ + a·((Ω cᵀ) b): the error names the IBP target η̇ and the term."""
        s = ps.so3_outer_product_on_velocity()
        with pytest.raises(AlgebraicError, match="matrix factor") as info:
            compute_eom(s.L, s.dW, s.variables)
        msg = str(info.value)
        eta_dot = str(s.symbols["R"].get_variation_vector().t_diff())
        assert f"the target '{eta_dot}'" in msg
        assert "(c)'" in msg and "cannot be integrated by parts" in msg

    def test_pipeline_s2_outer_product_on_the_velocity(self):
        s = ps.s2_outer_product_on_velocity()
        with pytest.raises(AlgebraicError, match="matrix factor") as info:
            compute_eom(s.L, s.dW, s.variables)
        xi_dot = str(s.symbols["q"].get_variation_vector().t_diff())
        assert f"the target '{xi_dot}'" in str(info.value)
