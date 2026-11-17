"""Pinned equations for the S² variation rule δω = ξ̇ − 2 ω × ξ (fix plan 1.1).

Lee et al. eq. (5.13) give δω = −S(ω)ξ + (I − qqᵀ)ξ̇ for q̇ = ω × q, δq = ξ × q, ω ⊥ q,
ξ ⊥ q. Differentiating ξ·q = 0 gives q·ξ̇ = q·(ω × ξ), and ω × ξ ∥ q, so qqᵀξ̇ = ω × ξ and
δω = ξ̇ − 2 ω × ξ. The old rule ξ̇ − ω × ξ was the tangential projection of δω only; it was
invisible while ∂L/∂ω ⊥ q and halved every gyroscopic term otherwise.

With p = ∂L/∂ω the pipeline produces, for every S² variable, the coefficient of ξ in
p·δω + (∂L/∂q)·(ξ × q) + δW after integration by parts, i.e. on the tangent plane

    E ≡ −ṗ + 2 ω × p + q × ∂L/∂q + f        (mod q).

Each test below pins the printed equation, derives it by hand from this formula in its
docstring and checks it with the finite-difference Euler-Lagrange check
(``geomech.numeric.assert_eom``), which uses none of the library's variational rules.
"""

import numpy as np

from geomech.numeric import evaluate, hat
from tests import phase1_systems as ps
from tests.phase1_systems import derive, equation


class TestGyroscopicTerm:
    def test_pinnedequation(self):
        """L = ½ m l² ω·ω + c·ω − m g l q·e₃, δW = ξ·f (dynamics audit F1, row 3).

        p = m l² ω + c, ṗ = m l² ω̇, 2 ω × p ≡ 2 ω × c (mod q, since ω × ω = 0),
        q × ∂L/∂q = −m g l q × e₃.  Hence E = −m l² ω̇ + 2 ω × c − m g l q × e₃ + f.
        The printed string is exactly this (before the fix the middle term was ω × c).
        """
        s, eom, _ = derive(ps.s2_gyroscopic)
        expected = (
            "(\\frac{d}{dt}(\\omega_{q})*m*l*l*(-1)+Cross(\\omega_{q},c)*(2)"
            "+Cross(q,e3)*m*g*l*(-1)+f)"
        )
        assert str(equation(eom, s.symbols["q"])) == expected


class TestAnisotropicInertiaOmegaForm:
    def test_pinnedequation(self):
        """L = ½ ω·Aω − m g l q·e₃, A constant symmetric (dynamics audit F1, row 1).

        p = Aω, ṗ = Aω̇, 2 ω × p = 2 ω × Aω, q × ∂L/∂q = −m g l q × e₃.
        E = −Aω̇ + 2 ω × Aω − m g l q × e₃ + f.
        """
        s, eom, _ = derive(ps.s2_anisotropic_omega)
        expected = (
            "(A*\\frac{d}{dt}(\\omega_{q})*(-1)+Cross(\\omega_{q},A*\\omega_{q})*(2)"
            "+Cross(q,e3)*m*g*l*(-1)+f)"
        )
        assert str(equation(eom, s.symbols["q"])) == expected


class TestRotatingSphere:
    """Lee et al. §5.6.2: particle on a sphere rotating about e₃ at rate Ω.

    Book equation (5.43), with a = Ω e₃, P = I − qqᵀ, C = diag(1, 1, 0):

        m L² q̈ + m L² |q̇|² q + 2 m L² P S(a) q̇ − m L² Ω² P C q = 0.

    Hand derivation of the ω-form.  With p = ∂L/∂ω = m L² (ω + q × (a × q)) one has
    ṗ = m L² (ω̇ − (ω × q)(q·a) − q((ω × q)·a)), 2 ω × p = 2 m L² (ω × a − (ω × q)(q·a)) and
    q × ∂L/∂q = m L² (−(q × ω)(q·a) + Ω² q × Cq); the (ω × q)(q·a) terms cancel, leaving

        E ≡ −m L² ω̇ + 2 m L² Ω ω × e₃ + m L² Ω² q × Cq + f      (mod q),

    which is −q × (5.43) + f on the tangent plane because q × q̈ = ω̇ and
    q × P S(a) q̇ = −(ω × q)(q·a) ≡ −ω × a (mod q).  Here ∂L/∂ω ⊥ q, so this system is
    insensitive to the normal component of δω; it pins the printed form and cross-checks the
    pipeline against the book in both the ω-form and the (q, q̇)-form of the Lagrangian.
    """

    OMEGA_FORM = (
        "(\\frac{d}{dt}(\\omega_{q})*m*L*L*(-1)+Cross(q,Cross(e3,Cross(q,\\omega_{q})))*m*L*L*\\Omega"
        "+Cross(q,Cross(\\omega_{q},Cross(e3,q)))*m*L*L*\\Omega*(-1)"
        "+Cross(Cross(Cross(\\omega_{q},q),e3),q)*m*L*L*\\Omega*(-1)"
        "+Cross(Cross(q,e3),Cross(\\omega_{q},q))*m*L*L*\\Omega*(-1)"
        "+Cross(\\omega_{q},Cross(q,Cross(e3,q)))*m*L*L*\\Omega*(2)+Cross(q,C*q)*m*L*L*\\Omega*\\Omega+M)"
    )
    QDOT_FORM = (
        "(Cross(q,Cross(\\frac{d}{dt}(\\omega_{q}),q))*m*L*L*(-1)"
        "+Cross(q,Cross(e3,Cross(\\omega_{q},q)))*m*L*L*\\Omega*(-1)"
        "+Cross(Cross(\\omega_{q},q),Cross(q,e3))*m*L*L*\\Omega"
        "+Cross(q,Cross(Cross(\\omega_{q},q),e3))*m*L*L*\\Omega"
        "+Cross(\\omega_{q},Cross(q,Cross(e3,q)))*m*L*L*\\Omega*(2)"
        "+Cross(q,Cross(\\omega_{q},Cross(e3,q)))*m*L*L*\\Omega*(-1)+Cross(q,C*q)*m*L*L*\\Omega*\\Omega"
        "+M)"
    )

    @staticmethod
    def _random_state(rng):
        q = rng.normal(size=(3, 1))
        q /= np.linalg.norm(q)
        P = np.eye(3) - q @ q.T
        return {
            "q": q,
            "\\omega_{q}": P @ rng.normal(size=(3, 1)),
            "\\frac{d}{dt}(\\omega_{q})": P @ rng.normal(size=(3, 1)),
            "e3": np.array([[0.0], [0.0], [1.0]]),
            "C": np.diag([1.0, 1.0, 0.0]),
            "\\Omega": float(rng.normal()),
            "m": float(rng.uniform(0.5, 2.0)),
            "L": float(rng.uniform(0.5, 2.0)),
            "M": rng.normal(size=(3, 1)),
        }

    @staticmethod
    def _book_543(env):
        """−q × (left-hand side of Lee et al. 5.43) + M, with q̈ from ω̇ and ω."""
        q, w, wd = env["q"], env["\\omega_{q}"], env["\\frac{d}{dt}(\\omega_{q})"]
        m, L, Om, C = env["m"], env["L"], env["\\Omega"], env["C"]
        a = Om * env["e3"]
        P = np.eye(3) - q @ q.T
        qd = hat(w) @ q
        qdd = hat(wd) @ q + hat(w) @ (hat(w) @ q)
        g = (
            m * L**2 * qdd
            + m * L**2 * float((qd.T @ qd)[0, 0]) * q
            + 2 * m * L**2 * P @ hat(a) @ qd
            - m * L**2 * Om**2 * P @ C @ q
        )
        return -hat(q) @ g + env["M"]

    @staticmethod
    def _hand_form(env):
        """−m L² ω̇ + 2 m L² Ω ω × e₃ + m L² Ω² q × Cq + M (tangent part)."""
        q, w, wd = env["q"], env["\\omega_{q}"], env["\\frac{d}{dt}(\\omega_{q})"]
        m, L, Om, C = env["m"], env["L"], env["\\Omega"], env["C"]
        return (
            -m * L**2 * wd
            + 2 * m * L**2 * Om * hat(w) @ env["e3"]
            + m * L**2 * Om**2 * hat(q) @ (C @ q)
            + env["M"]
        )

    def _check_against_book(self, eq):
        rng = np.random.default_rng(1)
        for _ in range(5):
            env = self._random_state(rng)
            q = env["q"]
            P = np.eye(3) - q @ q.T
            e_lib = evaluate(eq, env)
            for e_ref in (self._book_543(env), self._hand_form(env)):
                assert np.max(np.abs(P @ (e_lib - e_ref))) < 1e-9

    def test_omega_form_pinned_and_matches_book(self):
        s, eom, _ = derive(ps.s2_rotating_sphere_omega)
        eq = equation(eom, s.symbols["q"])
        assert str(eq) == self.OMEGA_FORM
        self._check_against_book(eq)

    def test_qdot_form_pinned_and_matches_book(self):
        s, eom, _ = derive(ps.s2_rotating_sphere_qdot)
        eq = equation(eom, s.symbols["q"])
        assert str(eq) == self.QDOT_FORM
        self._check_against_book(eq)

    def test_both_forms_agree(self):
        """The ω-form and the q̇-form of the same Lagrangian give the same tangent equation."""
        s1, eom1, _ = derive(ps.s2_rotating_sphere_omega)
        s2, eom2, _ = derive(ps.s2_rotating_sphere_qdot)
        eq1, eq2 = equation(eom1, s1.symbols["q"]), equation(eom2, s2.symbols["q"])
        rng = np.random.default_rng(2)
        for _ in range(5):
            env = self._random_state(rng)
            P = np.eye(3) - env["q"] @ env["q"].T
            assert np.max(np.abs(P @ (evaluate(eq1, env) - evaluate(eq2, env)))) < 1e-9


class TestRotatingFrame:
    def test_pinnedequation(self):
        """L = ½ m l² (|ω + W₀|² − (W₀·q)²) − m g l q·e₃ (dynamics audit F1, row 2).

        p = m l² (ω + W₀), ṗ = m l² ω̇, 2 ω × p ≡ 2 m l² ω × W₀,
        q × ∂L/∂q = −m l² (W₀·q) q × W₀ − m g l q × e₃.
        E = −m l² ω̇ + 2 m l² ω × W₀ − m l² (W₀·q) q × W₀ − m g l q × e₃ + f.
        """
        s, eom, _ = derive(ps.s2_rotating_frame)
        expected = (
            "(\\frac{d}{dt}(\\omega_{q})*m*l*l*(-1)+Cross(\\omega_{q},W0)*m*l*l*(2)"
            "+Cross(q,W0)*m*l*l*Dot(W0,q)*(-1)+Cross(q,e3)*m*g*l*(-1)+f)"
        )
        assert str(equation(eom, s.symbols["q"])) == expected


class TestR3xS2VelocityCoupling:
    def test_pinned_equations(self):
        """L = ½ m ẋ·ẋ + ½ m l² ω·ω + ẋ·ω (dynamics audit F1, row 4).

        x: p_x = m ẋ + ω,  E_x = −ṗ_x + F = −m ẍ − ω̇ + F.
        q: p = m l² ω + ẋ,  E = −m l² ω̇ − ẍ + 2 ω × ẋ + f.
        """
        s, eom, _ = derive(ps.r3_s2_velocity_coupling)
        assert str(eom[str(s.symbols["x"].delta())][1]) == (
            "(\\frac{d}{dt}(\\frac{d}{dt}(x))*m*(-1)+\\frac{d}{dt}(\\omega_{q})*(-1)+F)"
        )
        assert str(equation(eom, s.symbols["q"])) == (
            "(\\frac{d}{dt}(\\omega_{q})*m*l*l*(-1)+\\frac{d}{dt}(\\frac{d}{dt}(x))*(-1)"
            "+Cross(\\omega_{q},\\frac{d}{dt}(x))*(2)+f)"
        )


class TestS2xSO3Couplings:
    def test_omega_dot_Omega(self):
        """L = ½ m ω·ω + ½ Ω·JΩ + k ω·Ω (dynamics audit F1, row 5).

        ξ: p = m ω + k Ω,  E = −m ω̇ − k Ω̇ + 2 k ω × Ω + f.
        η: p_Ω = JΩ + k ω, δΩ = η̇ + Ω × η, so E_η = −ṗ_Ω + p_Ω × Ω + M
           = −JΩ̇ − Ω × JΩ − k ω̇ − k Ω × ω + M.
        """
        s, eom, _ = derive(ps.s2_so3_omega_coupling)
        assert str(equation(eom, s.symbols["q"])) == (
            "(\\frac{d}{dt}(\\omega_{q})*m*(-1)+\\frac{d}{dt}(\\Omega_{R})*k*(-1)"
            "+Cross(\\omega_{q},\\Omega_{R})*k*(2)+f)"
        )
        assert str(equation(eom, s.symbols["R"])) == (
            "(J*\\frac{d}{dt}(\\Omega_{R})*(-1)+Cross(\\Omega_{R},J*\\Omega_{R})*(-1)"
            "+\\frac{d}{dt}(\\omega_{q})*k*(-1)+Cross(\\Omega_{R},\\omega_{q})*k*(-1)+M)"
        )

    def test_omega_dot_body_axis(self):
        """L = ½ m ω·ω + ½ Ω·JΩ + k ω·(Rᵀe₃) (dynamics audit F1, row 6).

        ξ: p = m ω + k Rᵀe₃, ṗ = m ω̇ − k Ω × (Rᵀe₃) (Ṙᵀ = −Ω̂ Rᵀ),
           E = −m ω̇ + k Ω × (Rᵀe₃) + 2 k ω × (Rᵀe₃) + f.
        η: δ(Rᵀe₃) = (Rᵀe₃) × η, so k ω·((Rᵀe₃) × η) = −k η·((Rᵀe₃) × ω) and
           E_η = −JΩ̇ − Ω × JΩ − k (Rᵀe₃) × ω + M.
        """
        s, eom, _ = derive(ps.s2_so3_axis_coupling)
        assert str(equation(eom, s.symbols["q"])) == (
            "(\\frac{d}{dt}(\\omega_{q})*m*(-1)+Cross(\\Omega_{R},(R)'*e3)*k"
            "+Cross(\\omega_{q},(R)'*e3)*k*(2)+f)"
        )
        assert str(equation(eom, s.symbols["R"])) == (
            "(J*\\frac{d}{dt}(\\Omega_{R})*(-1)+Cross(\\Omega_{R},J*\\Omega_{R})*(-1)"
            "+Cross((R)'*e3,\\omega_{q})*k*(-1)+M)"
        )
