"""Variation-based linearization (``geomech.dynamics.vbl``).

Every expected block below is derived by hand from the equation of motion pinned in
``tests/test_manifold_pipelines.py`` / ``tests/test_scalar_adversarial.py`` (quoted as ``E`` in
each docstring). Convention (``vbl`` module docstring): ``E`` is the negative of the
Euler–Lagrange expression, and

    δE = −M δv̇ + A δx + B δu,   M[δv̇] = −coeff(δv̇, δE),  A[δx] = coeff(δx, δE),  B = coeff(δu, δE)

so the physical inertia appears positive in ``M``. Besides the exact or numeric block
checks, every system is verified by the finite-difference check of
``tests/numeric_helpers.py``: ``(E(ε) − E(−ε)) / 2ε`` along a constraint-respecting
perturbation of a random trajectory point must equal ``−M δv̇ + A δx + B δu``.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pytest

from geomech import (
    SO3,
    Dot,
    Matrix,
    Scalar,
    SVMul,
    SystemVariables,
    TimeDerivative,
    Variation,
    Vector,
    compute_eom,
    linearize,
    to_latex,
)
from geomech.core.base.types import ExprType
from geomech.dynamics import LinearizedSystem
from geomech.utils.errors import AlgebraicError
from tests import vbl_systems
from tests.numeric_helpers import (
    Evaluator,
    expm_so3,
    finite_difference_check,
    hat,
    linearized_rhs,
    sample_direction,
    sample_state,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def _derive(name):
    """(system, eom, lin) for the named constructor of ``tests/vbl_systems.py``."""
    system = getattr(vbl_systems, name)()
    eom = compute_eom(system.L, system.dW, system.variables)
    lin = linearize(eom, system.variables, system.inputs)
    return system, eom, lin


def _zero_blocks_except(blocks, nonzero):
    """Every block not listed in *nonzero* must be a zero constant."""
    for col, block in blocks.items():
        if col not in nonzero:
            assert block.is_zero, f"expected zero block at column {col}, got {block}"


def _all_block_strings(lin):
    for row in lin.rows:
        for table in (lin.M, lin.A, lin.B):
            for block in table[row].values():
                yield str(block)


# ---------------------------------------------------------------------------
# 1. Point mass in R³
# ---------------------------------------------------------------------------


class TestPointMassR3:
    """E = −m ẍ + F, so δE = −m δẍ + δF.

    States [δx, δẋ]. Kinematic row d/dt(δx) = δẋ: M[δẋ] = I, A[δẋ] = I.
    Dynamic row: M[δẍ] = −(−m I) = m I, A = 0, B[δF] = I.
    """

    def setup_method(self):
        self.s, self.eom, self.lin = _derive("point_mass_r3")
        x = self.s.symbols["x"]
        self.dx = x.delta()
        self.dv = self.dx.t_diff()
        self.da = self.dv.t_diff()
        self.dF = self.s.symbols["F"].delta()

    def test_pinned_eom(self):
        assert str(self.eom[str(self.dx)][1]) == "(\\frac{d}{dt}(\\frac{d}{dt}(x))*m*(-1)+F)"

    def test_structure(self):
        assert self.lin.rows == ["kin:x", str(self.dx)]
        assert [str(s) for s in self.lin.states] == [str(self.dx), str(self.dv)]
        assert [str(d) for d in self.lin.dstates] == [str(self.dv), str(self.da)]
        assert [str(u) for u in self.lin.inputs] == [str(self.dF)]

    def test_kinematic_row_exact(self):
        row = "kin:x"
        assert str(self.lin.M[row][str(self.dv)]) == "I"
        assert str(self.lin.A[row][str(self.dv)]) == "I"
        _zero_blocks_except(self.lin.M[row], {str(self.dv)})
        _zero_blocks_except(self.lin.A[row], {str(self.dv)})
        _zero_blocks_except(self.lin.B[row], set())

    def test_dynamic_row_exact(self):
        row = str(self.dx)
        assert str(self.lin.M[row][str(self.da)]) == "I*m"
        assert str(self.lin.B[row][str(self.dF)]) == "I"
        _zero_blocks_except(self.lin.M[row], {str(self.da)})
        _zero_blocks_except(self.lin.A[row], set())

    def test_zero_blocks_have_matrix_shape(self):
        for row in self.lin.rows:
            for col, block in self.lin.A[row].items():
                assert self.lin.block_shape(row, col) == (3, 3)
                if block.is_zero:
                    assert block.type == ExprType.MATRIX

    def test_finite_difference(self):
        finite_difference_check(self.lin, self.eom, self.s.variables, self.s.inputs)


# ---------------------------------------------------------------------------
# 2. Spring-mass on R
# ---------------------------------------------------------------------------


class TestSpringMassR:
    """E = −m s̈ − k s + f (pinned in ``TestSpringMassPipeline``), δE = −m δs̈ − k δs + δf.

    States [δs, δṡ]. Kinematic row: M[δṡ] = 1, A[δṡ] = 1.
    Dynamic row: M[δs̈] = m, A[δs] = −k, A[δṡ] = 0, B[δf] = 1.
    """

    def setup_method(self):
        self.s, self.eom, self.lin = _derive("spring_mass_r")
        self.dx = self.s.symbols["x"].delta()
        self.dv = self.dx.t_diff()
        self.da = self.dv.t_diff()
        self.df = self.s.symbols["f"].delta()

    def test_pinned_eom(self):
        expected = "(m*\\frac{d}{dt}(\\frac{d}{dt}(x))*(-1)+k*x*(-1)+f)"
        assert str(self.eom[str(self.dx)][1]) == expected

    def test_kinematic_row_exact(self):
        row = "kin:x"
        assert str(self.lin.M[row][str(self.dv)]) == "1"
        assert str(self.lin.A[row][str(self.dv)]) == "1"
        _zero_blocks_except(self.lin.M[row], {str(self.dv)})
        _zero_blocks_except(self.lin.A[row], {str(self.dv)})
        assert str(self.lin.B[row][str(self.df)]) == "0"

    def test_dynamic_row_exact(self):
        row = str(self.dx)
        assert str(self.lin.M[row][str(self.da)]) == "m"
        assert str(self.lin.A[row][str(self.dx)]) == "k*(-1)"
        assert str(self.lin.B[row][str(self.df)]) == "1"
        _zero_blocks_except(self.lin.M[row], {str(self.da)})
        _zero_blocks_except(self.lin.A[row], {str(self.dx)})

    def test_block_shapes_are_scalar(self):
        for row in self.lin.rows:
            for col in list(self.lin.M[row]) + list(self.lin.A[row]) + list(self.lin.B[row]):
                assert self.lin.block_shape(row, col) == (1, 1)

    def test_finite_difference(self):
        finite_difference_check(self.lin, self.eom, self.s.variables, self.s.inputs)


# ---------------------------------------------------------------------------
# 3. Free rigid body on SO(3)
# ---------------------------------------------------------------------------


class TestFreeRigidBody:
    """Lee et al. Eq. 6.16. E = −J Ω̇ − Ω × JΩ + M (pinned in ``TestFreeRigidBodyPipeline``).

    δE = −J δΩ̇ − δΩ × JΩ − Ω × J δΩ + δM = −J δΩ̇ + (Hat(JΩ) − Hat(Ω) J) δΩ + δM,
    using a × b = Hat(a) b and δΩ × JΩ = −Hat(JΩ) δΩ.

    States [η, δΩ]. Kinematic row η̇ = −Hat(Ω) η + δΩ: M[η̇] = I, A[η] = −Hat(Ω), A[δΩ] = I.
    Dynamic row: M[δΩ̇] = J, A[η] = 0, A[δΩ] = Hat(JΩ) − Hat(Ω) J, B[δM] = I.
    """

    def setup_method(self):
        self.s, self.eom, self.lin = _derive("free_rigid_body")
        R = self.s.symbols["R"]
        self.Om, self.eta = R.get_tangent_vector(), R.get_variation_vector()
        self.dOm = Variation(self.Om)
        self.ddOm = self.dOm.t_diff()
        self.dM = self.s.symbols["M"].delta()

    def test_pinned_eom(self):
        expected = "(J*\\frac{d}{dt}(\\Omega_{R})*(-1)+Cross(\\Omega_{R},J*\\Omega_{R})*(-1)+M)"
        assert str(self.eom[str(self.eta)][1]) == expected

    def test_structure(self):
        assert self.lin.rows == ["kin:R", str(self.eta)]
        assert [str(s) for s in self.lin.states] == [str(self.eta), str(self.dOm)]
        assert [str(d) for d in self.lin.dstates] == [str(self.eta.t_diff()), str(self.ddOm)]

    def test_kinematic_row_exact(self):
        row = "kin:R"
        assert str(self.lin.M[row][str(self.eta.t_diff())]) == "I"
        assert str(self.lin.A[row][str(self.eta)]) == "Hat(\\Omega_{R})*(-1)"
        assert str(self.lin.A[row][str(self.dOm)]) == "I"
        _zero_blocks_except(self.lin.M[row], {str(self.eta.t_diff())})
        _zero_blocks_except(self.lin.B[row], set())

    def test_dynamic_row_exact(self):
        row = str(self.eta)
        assert str(self.lin.M[row][str(self.ddOm)]) == "J"
        assert str(self.lin.A[row][str(self.dOm)]) == (
            "(Hat(J*\\Omega_{R})+Hat(\\Omega_{R})*J*(-1))"
        )
        assert str(self.lin.B[row][str(self.dM)]) == "I"
        _zero_blocks_except(self.lin.M[row], {str(self.ddOm)})
        _zero_blocks_except(self.lin.A[row], {str(self.dOm)})

    def test_dynamic_row_numeric(self):
        row = str(self.eta)
        for seed in range(3):
            ev = Evaluator(rng=np.random.default_rng(seed))
            A_dOm = ev(self.lin.A[row][str(self.dOm)])
            J, W = ev.env["J"], ev.env[str(self.Om)]
            expected = hat(J @ W) - hat(W) @ J
            assert np.allclose(A_dOm, expected, atol=1e-12)
            assert np.allclose(ev(self.lin.M[row][str(self.ddOm)]), J)
            assert np.allclose(ev(self.lin.A["kin:R"][str(self.eta)]), -hat(W))

    def test_finite_difference(self):
        finite_difference_check(self.lin, self.eom, self.s.variables, self.s.inputs)


# ---------------------------------------------------------------------------
# 4. Rigid pendulum on SO(3), inertia about the pivot
# ---------------------------------------------------------------------------


class TestRigidPendulumPivot:
    """Lee et al. Eq. 6.8. E = −J Ω̇ − Ω × JΩ − m g ρ × Rᵀe₃ + M
    (pinned in ``TestRigidPendulumPivotPipeline``).

    δ(Rᵀe₃) = (R Hat(η))ᵀ e₃ = −Hat(η) Rᵀe₃ = Hat(Rᵀe₃) η, so
    δ(ρ × Rᵀe₃) = Hat(ρ) Hat(Rᵀe₃) η and the gravity block is
    A[η] = −m g Hat(ρ) Hat(Rᵀe₃). The remaining blocks are those of the free rigid body.
    """

    def setup_method(self):
        self.s, self.eom, self.lin = _derive("rigid_pendulum_pivot")
        R = self.s.symbols["R"]
        self.Om, self.eta = R.get_tangent_vector(), R.get_variation_vector()
        self.dOm = Variation(self.Om)
        self.ddOm = self.dOm.t_diff()

    def test_pinned_eom(self):
        expected = (
            "(J*\\frac{d}{dt}(\\Omega_{R})*(-1)+Cross(\\Omega_{R},J*\\Omega_{R})*(-1)"
            "+Cross(\\rho,(R)'*e3)*m*g*(-1)+M)"
        )
        assert str(self.eom[str(self.eta)][1]) == expected

    def test_gravity_block_exact(self):
        row = str(self.eta)
        assert str(self.lin.A[row][str(self.eta)]) == "Hat(\\rho)*Hat((R)'*e3)*m*g*(-1)"
        assert str(self.lin.M[row][str(self.ddOm)]) == "J"
        assert str(self.lin.A[row][str(self.dOm)]) == (
            "(Hat(J*\\Omega_{R})+Hat(\\Omega_{R})*J*(-1))"
        )

    def test_blocks_numeric(self):
        row = str(self.eta)
        for seed in range(3):
            ev = Evaluator(rng=np.random.default_rng(seed))
            A_eta = ev(self.lin.A[row][str(self.eta)])
            A_dOm = ev(self.lin.A[row][str(self.dOm)])
            M_block = ev(self.lin.M[row][str(self.ddOm)])
            R, rho, e3 = ev.env["R"], ev.env["\\rho"], ev.env["e3"]
            m, g = ev.env["m"], ev.env["g"]
            J, W = ev.env["J"], ev.env[str(self.Om)]
            assert np.allclose(A_eta, -m * g * hat(rho) @ hat(R.T @ e3), atol=1e-12)
            assert np.allclose(A_dOm, hat(J @ W) - hat(W) @ J)
            assert np.allclose(M_block, J)

    def test_finite_difference(self):
        finite_difference_check(self.lin, self.eom, self.s.variables, self.s.inputs)


# ---------------------------------------------------------------------------
# 5. Spherical pendulum on S², ω-form
# ---------------------------------------------------------------------------


class TestSphericalPendulumS2:
    """Lee et al. Eq. 5.14. E = −m l² ω̇ − m g l q × e₃ + M (``TestSphericalPendulumPipeline``).

    δq = ξ × q, so δ(q × e₃) = (ξ × q) × e₃ = −Hat(e₃) Hat(ξ) q = Hat(e₃) Hat(q) ξ and
    δE = −m l² δω̇ − m g l Hat(e₃) Hat(q) ξ + δM.

    States [ξ, δω]. Kinematic row ξ̇ = q qᵀ Hat(ω) ξ + (I − q qᵀ) δω.
    Dynamic row: M[δω̇] = m l² I, A[ξ] = −m g l Hat(e₃) Hat(q), A[δω] = 0, B[δM] = I.

    Sanity check at q = −e₃, ω = 0: Hat(e₃) Hat(−e₃) ξ = −e₃ × (e₃ × ξ) = ξ for ξ ⊥ e₃,
    and ξ̇ = δω on the tangent plane, so m l² ξ̈ = −m g l ξ: the pendulum frequency √(g/l).
    """

    def setup_method(self):
        self.s, self.eom, self.lin = _derive("spherical_pendulum_s2")
        q = self.s.symbols["q"]
        self.q = q
        self.om, self.xi = q.get_tangent_vector(), q.get_variation_vector()
        self.dom = Variation(self.om)
        self.ddom = self.dom.t_diff()
        self.dM = self.s.symbols["M"].delta()

    def test_pinned_eom(self):
        expected = "(\\frac{d}{dt}(\\omega_{q})*m*l*l*(-1)+Cross(q,e3)*m*g*l*(-1)+M)"
        assert str(self.eom[str(self.xi)][1]) == expected

    def test_kinematic_row_exact(self):
        row = "kin:q"
        assert str(self.lin.M[row][str(self.xi.t_diff())]) == "I"
        assert str(self.lin.A[row][str(self.xi)]) == "q*(q)'*Hat(\\omega_{q})"
        assert str(self.lin.A[row][str(self.dom)]) == "(I+q*(q)'*(-1))"
        _zero_blocks_except(self.lin.M[row], {str(self.xi.t_diff())})
        _zero_blocks_except(self.lin.B[row], set())

    def test_dynamic_row_exact(self):
        row = str(self.xi)
        assert str(self.lin.M[row][str(self.ddom)]) == "I*m*l*l"
        assert str(self.lin.A[row][str(self.xi)]) == "Hat(e3)*Hat(q)*m*g*l*(-1)"
        assert str(self.lin.B[row][str(self.dM)]) == "I"
        _zero_blocks_except(self.lin.M[row], {str(self.ddom)})
        _zero_blocks_except(self.lin.A[row], {str(self.xi)})

    def test_hanging_equilibrium_restoring_force(self):
        """At q = −e₃ the ξ block acts as −m g l on tangent vectors (pendulum stiffness)."""
        row = str(self.xi)
        e3 = np.array([[0.0], [0.0], [1.0]])
        env = {"q": -e3, "e3": e3, "m": 2.0, "g": 9.81, "l": 0.5}
        ev = Evaluator(env=env)
        A_xi = ev(self.lin.A[row][str(self.xi)])
        xi = np.array([[0.3], [-0.7], [0.0]])  # tangent at q = −e₃
        assert np.allclose(A_xi @ xi, -2.0 * 9.81 * 0.5 * xi)
        A_kin_dom = ev(self.lin.A["kin:q"][str(self.dom)])
        assert np.allclose(A_kin_dom @ xi, xi)  # (I − q qᵀ) is the identity on the tangent plane

    def test_kinematic_row_numeric(self):
        """ξ̇ from the kinematic row equals d/dε of ξ along the perturbed curve.

        With q(ε) = exp(ε ξ̂) q and ω(ε) = exp(ε ξ̂)(ω + ε w) the tangency of ξ is preserved;
        the row must reproduce q qᵀ ω̂ ξ + (I − q qᵀ)(ξ × ω + w).
        """
        rng = np.random.default_rng(3)
        env = sample_state(self.s.variables, self.s.inputs, rng)
        ev = Evaluator(env=env, rng=rng)
        raw, pert = sample_direction(self.s.variables, self.s.inputs, ev.env, rng)
        q, w = ev.env["q"], ev.env[str(self.om)]
        xi, dom = pert[str(self.xi)], pert[str(self.dom)]
        lhs = (
            ev(self.lin.A["kin:q"][str(self.xi)]) @ xi
            + ev(self.lin.A["kin:q"][str(self.dom)]) @ dom
        )
        P = q @ q.T
        assert np.allclose(lhs, P @ hat(w) @ xi + (np.eye(3) - P) @ dom)
        assert abs((q.T @ lhs).item() - (q.T @ hat(w) @ xi).item()) < 1e-12  # d/dt(ξ·q) = 0

    def test_finite_difference(self):
        finite_difference_check(self.lin, self.eom, self.s.variables, self.s.inputs)


# ---------------------------------------------------------------------------
# 6. Quadrotor on R³ × SO(3), scalar thrust input
# ---------------------------------------------------------------------------


class TestQuadrotor:
    """Translational E = −m ẍ − m g e₃ + f R e₃; rotational E as the free rigid body.

    δ(f R e₃) = δf R e₃ + f R Hat(η) e₃ = R e₃ δf − f R Hat(e₃) η, so the translational row
    has M[δẍ] = m I, A[η] = −f R Hat(e₃), B[δf] = R e₃ (a 3×1 block for the scalar input)
    and zero blocks for δx, δẋ, δΩ and δM. States [δx, η, δẋ, δΩ]; inputs [δf, δM].
    """

    def setup_method(self):
        self.s, self.eom, self.lin = _derive("quadrotor")
        x, R = self.s.symbols["x"], self.s.symbols["R"]
        self.dx = x.delta()
        self.Om, self.eta = R.get_tangent_vector(), R.get_variation_vector()
        self.dOm = Variation(self.Om)
        self.df = self.s.symbols["f"].delta()
        self.dM = self.s.symbols["M"].delta()

    def test_pinned_eom(self):
        assert str(self.eom[str(self.dx)][1]) == (
            "(\\frac{d}{dt}(\\frac{d}{dt}(x))*m*(-1)+e3*m*g*(-1)+R*e3*f)"
        )

    def test_structure(self):
        assert self.lin.rows == ["kin:x", "kin:R", str(self.dx), str(self.eta)]
        assert [str(s) for s in self.lin.states] == [
            str(self.dx),
            str(self.eta),
            str(self.dx.t_diff()),
            str(self.dOm),
        ]
        assert [str(u) for u in self.lin.inputs] == [str(self.df), str(self.dM)]

    def test_translational_row_exact(self):
        row = str(self.dx)
        assert str(self.lin.M[row][str(self.dx.t_diff().t_diff())]) == "I*m"
        assert str(self.lin.A[row][str(self.eta)]) == "R*Hat(e3)*f*(-1)"
        assert str(self.lin.B[row][str(self.df)]) == "R*e3"
        _zero_blocks_except(self.lin.M[row], {str(self.dx.t_diff().t_diff())})
        _zero_blocks_except(self.lin.A[row], {str(self.eta)})
        _zero_blocks_except(self.lin.B[row], {str(self.df)})

    def test_rotational_row_exact(self):
        row = str(self.eta)
        assert str(self.lin.M[row][str(self.dOm.t_diff())]) == "J"
        assert str(self.lin.A[row][str(self.dOm)]) == (
            "(Hat(J*\\Omega_{R})+Hat(\\Omega_{R})*J*(-1))"
        )
        assert str(self.lin.B[row][str(self.dM)]) == "I"
        _zero_blocks_except(self.lin.A[row], {str(self.dOm)})
        _zero_blocks_except(self.lin.B[row], {str(self.dM)})

    def test_block_shapes(self):
        row = str(self.dx)
        assert self.lin.block_shape(row, self.df) == (3, 1)
        assert self.lin.block_shape(row, self.dM) == (3, 3)
        assert self.lin.block_shape(row, self.eta) == (3, 3)
        assert self.lin.B[row][str(self.dM)].type == ExprType.MATRIX
        assert self.lin.B[str(self.eta)][str(self.df)].type == ExprType.VECTOR
        assert self.lin.B[str(self.eta)][str(self.df)].is_zero

    def test_translational_row_numeric(self):
        row = str(self.dx)
        for seed in range(3):
            ev = Evaluator(rng=np.random.default_rng(seed))
            A_eta = ev(self.lin.A[row][str(self.eta)])
            B_f = ev(self.lin.B[row][str(self.df)])
            M_block = ev(self.lin.M[row][str(self.dx.t_diff().t_diff())])
            R, e3, f, m = ev.env["R"], ev.env["e3"], ev.env["f"], ev.env["m"]
            assert np.allclose(A_eta, -f * R @ hat(e3), atol=1e-12)
            assert np.allclose(B_f, R @ e3)
            assert np.allclose(M_block, m * np.eye(3))

    def test_finite_difference(self):
        finite_difference_check(self.lin, self.eom, self.s.variables, self.s.inputs)


# ---------------------------------------------------------------------------
# 7. Double rigid pendulum on SO(3) × SO(3)
# ---------------------------------------------------------------------------


class TestDoubleRigidPendulum:
    """System of ``TestDoubleRigidPendulumMassMatrixExact``; finite differences plus exact
    kinematic rows (one SO(3) row per body, η̇ᵢ = −Hat(Ωᵢ) ηᵢ + δΩᵢ, no cross coupling)."""

    def setup_method(self):
        self.s, self.eom, self.lin = _derive("double_rigid_pendulum")
        self.R1, self.R2 = self.s.symbols["R1"], self.s.symbols["R2"]

    def test_structure(self):
        eta1, eta2 = self.R1.get_variation_vector(), self.R2.get_variation_vector()
        assert self.lin.rows == ["kin:R1", "kin:R2", str(eta1), str(eta2)]
        assert [str(s) for s in self.lin.states] == [
            str(eta1),
            str(eta2),
            str(Variation(self.R1.get_tangent_vector())),
            str(Variation(self.R2.get_tangent_vector())),
        ]

    def test_kinematic_rows_exact(self):
        for R in (self.R1, self.R2):
            row = "kin:" + str(R)
            eta, Om = R.get_variation_vector(), R.get_tangent_vector()
            assert str(self.lin.M[row][str(eta.t_diff())]) == "I"
            assert str(self.lin.A[row][str(eta)]) == f"Hat(\\Omega_{{{R}}})*(-1)"
            assert str(self.lin.A[row][str(Variation(Om))]) == "I"
            _zero_blocks_except(self.lin.M[row], {str(eta.t_diff())})
            _zero_blocks_except(self.lin.A[row], {str(eta), str(Variation(Om))})
            _zero_blocks_except(self.lin.B[row], set())

    def test_mass_blocks_match_standard_form_inertia(self):
        """M[ηᵢ][δΩ̇ⱼ] is the negative of the pinned standard-form block (M = −coeff)."""
        eta1, eta2 = self.R1.get_variation_vector(), self.R2.get_variation_vector()
        dOm1 = str(Variation(self.R1.get_tangent_vector()).t_diff())
        dOm2 = str(Variation(self.R2.get_tangent_vector()).t_diff())
        M = self.lin.M
        assert (
            str(M[str(eta1)][dOm1])
            == "(J1+Hat(\\rho_1)*Hat(\\rho_1)*m1*(-1)+Hat(l_1)*Hat(l_1)*m2*(-1))"
        )
        assert str(M[str(eta1)][dOm2]) == "Hat(l_1)*(R1)'*R2*Hat(\\rho_2)*m2*(-1)"
        assert str(M[str(eta2)][dOm1]) == "Hat(\\rho_2)*(R2)'*R1*Hat(l_1)*m2*(-1)"
        assert str(M[str(eta2)][dOm2]) == "(J2+Hat(\\rho_2)*Hat(\\rho_2)*m2*(-1))"

    def test_finite_difference(self):
        finite_difference_check(self.lin, self.eom, self.s.variables, self.s.inputs)


# ---------------------------------------------------------------------------
# 8. Elastic spherical pendulum on R × S²
# ---------------------------------------------------------------------------


class TestElasticSphericalPendulum:
    """System of ``tests/test_scalar_adversarial.py``; finite differences plus block shapes.

    The scalar row δl has 1×3 blocks against the vector perturbations (ξ, δω, δM_q) and the
    vector row ξ has 3×1 blocks against the scalar perturbations (δl, δl̇, δf_l).
    Hand check of one 1×3 block: E_l contains m l ω·ω, whose variation in δω is
    2 m l ωᵀ δω, so A[δl][δω] = (2 m l ω)ᵀ.
    """

    def setup_method(self):
        self.s, self.eom, self.lin = _derive("elastic_spherical_pendulum")
        self.l, self.q = self.s.symbols["l"], self.s.symbols["q"]
        self.dl = self.l.delta()
        self.xi = self.q.get_variation_vector()
        self.dom = Variation(self.q.get_tangent_vector())

    def test_structure(self):
        assert self.lin.rows == ["kin:l", "kin:q", str(self.dl), str(self.xi)]
        assert [str(s) for s in self.lin.states] == [
            str(self.dl),
            str(self.xi),
            str(self.dl.t_diff()),
            str(self.dom),
        ]

    def test_block_shapes_and_types(self):
        scalar_row, vector_row = str(self.dl), str(self.xi)
        assert self.lin.block_shape(scalar_row, self.xi) == (1, 3)
        assert self.lin.block_shape(scalar_row, self.dl) == (1, 1)
        assert self.lin.block_shape(vector_row, self.dl) == (3, 1)
        assert self.lin.block_shape(vector_row, self.dom) == (3, 3)
        for col, block in self.lin.A[scalar_row].items():
            shape = self.lin.block_shape(scalar_row, col)
            ev = Evaluator(rng=np.random.default_rng(0))
            assert np.atleast_2d(np.asarray(ev(block))).shape == shape, (col, block)
        for col, block in self.lin.A[vector_row].items():
            shape = self.lin.block_shape(vector_row, col)
            ev = Evaluator(rng=np.random.default_rng(0))
            assert np.atleast_2d(np.asarray(ev(block))).shape == shape, (col, block)

    def test_row_block_for_omega_exact(self):
        assert str(self.lin.A[str(self.dl)][str(self.dom)]) == "(\\omega_{q}*m*l*(2))'"

    def test_finite_difference(self):
        finite_difference_check(self.lin, self.eom, self.s.variables, self.s.inputs)


# ---------------------------------------------------------------------------
# 9. Quadrotor with cable-suspended load on R³ × S² × SO(3)
# ---------------------------------------------------------------------------


class TestQuadrotorWithLoad:
    """Sreenath, Lee, Kumar CDC 2013; finite differences on all three dynamic rows."""

    def setup_method(self):
        self.s, self.eom, self.lin = _derive("quadrotor_with_load")

    def test_structure(self):
        xL, q, R = self.s.symbols["xL"], self.s.symbols["q"], self.s.symbols["R"]
        assert self.lin.rows == [
            "kin:x_L",
            "kin:q",
            "kin:R",
            str(xL.delta()),
            str(q.get_variation_vector()),
            str(R.get_variation_vector()),
        ]
        assert len(self.lin.states) == 6
        assert [str(u) for u in self.lin.inputs] == ["\\delta{f}", "\\delta{M}"]

    def test_thrust_input_blocks(self):
        """B[δf] is R e₃ in the load row and −l q × R e₃ in the cable row.

        Both follow from δW = f (R e₃)·(δx_L − l ξ × q).
        """
        xL, q = self.s.symbols["xL"], self.s.symbols["q"]
        B = self.lin.B
        assert str(B[str(xL.delta())]["\\delta{f}"]) == "R*e3"
        assert str(B[str(q.get_variation_vector())]["\\delta{f}"]) == "Cross(q,R*e3)*l*(-1)"

    def test_finite_difference(self):
        finite_difference_check(self.lin, self.eom, self.s.variables, self.s.inputs)


# ---------------------------------------------------------------------------
# Cross-cutting properties
# ---------------------------------------------------------------------------

ALL_SYSTEMS = [
    "point_mass_r3",
    "spring_mass_r",
    "free_rigid_body",
    "rigid_pendulum_pivot",
    "spherical_pendulum_s2",
    "quadrotor",
    "double_rigid_pendulum",
    "elastic_spherical_pendulum",
    "quadrotor_with_load",
]


@pytest.mark.parametrize("name", ALL_SYSTEMS)
def test_blocks_are_dense_and_typed(name):
    s, eom, lin = _derive(name)
    assert isinstance(lin, LinearizedSystem)
    cols_M = [str(d) for d in lin.dstates]
    cols_A = [str(x) for x in lin.states]
    cols_B = [str(u) for u in lin.inputs]
    for row in lin.rows:
        assert list(lin.M[row]) == cols_M
        assert list(lin.A[row]) == cols_A
        assert list(lin.B[row]) == cols_B
        for table in (lin.M, lin.A, lin.B):
            for col, block in table[row].items():
                shape = lin.block_shape(row, col)
                if shape == (3, 3):
                    assert block.type == ExprType.MATRIX, (row, col, block)
                elif shape == (1, 1):
                    assert block.type == ExprType.SCALAR, (row, col, block)
                else:
                    assert block.type == ExprType.VECTOR, (row, col, block)


@pytest.mark.parametrize("name", ALL_SYSTEMS)
def test_no_unfolded_numeric_factors(name):
    """Blocks print without a `(-1)(-1)` chain: numeric factors are folded by simplify."""
    _, _, lin = _derive(name)
    for text in _all_block_strings(lin):
        assert "(-1)*(-1)" not in text, text
        assert "(-1.0)*(-1" not in text, text
        assert "(-1)*(-1.0" not in text, text


@pytest.mark.parametrize("name", ALL_SYSTEMS)
def test_blocks_are_free_of_perturbation_symbols(name):
    _, _, lin = _derive(name)
    symbols = [str(x) for x in lin.states + lin.dstates + lin.inputs]
    for text in _all_block_strings(lin):
        for sym in symbols:
            assert sym not in text, (sym, text)


@pytest.mark.parametrize("name", ALL_SYSTEMS)
def test_str_and_latex(name):
    _, _, lin = _derive(name)
    text = str(lin)
    for row in lin.rows:
        assert f"row {row}" in text
    for row in lin.rows:
        for table in (lin.M, lin.A, lin.B):
            for block in table[row].values():
                assert isinstance(to_latex(block), str)
                if not block.is_zero:
                    assert str(block) in text


def test_finite_difference_check_detects_a_wrong_block():
    """The numeric check is only evidence if it fails on a wrong linearization."""
    s, eom, lin = _derive("free_rigid_body")
    row = str(s.symbols["R"].get_variation_vector())
    dOm = str(Variation(s.symbols["R"].get_tangent_vector()))
    bad_A = {r: dict(blocks) for r, blocks in lin.A.items()}
    bad_A[row][dOm] = lin.A[row][dOm] * -1
    bad = LinearizedSystem(
        rows=lin.rows,
        states=lin.states,
        dstates=lin.dstates,
        inputs=lin.inputs,
        M=lin.M,
        A=bad_A,
        B=lin.B,
        row_types=lin.row_types,
        col_types=lin.col_types,
    )
    with pytest.raises(AssertionError):
        finite_difference_check(bad, eom, s.variables, s.inputs)


def test_linearized_rhs_matches_hand_formula_free_rigid_body():
    """−M δΩ̇ + A δx + B δu evaluates to −J δΩ̇ + (Hat(JΩ) − Hat(Ω)J) δΩ + δM."""
    s, eom, lin = _derive("free_rigid_body")
    R = s.symbols["R"]
    row = str(R.get_variation_vector())
    rng = np.random.default_rng(11)
    env = sample_state(s.variables, s.inputs, rng)
    ev = Evaluator(env=env, rng=rng)
    ev(eom[row][1])
    raw, pert = sample_direction(s.variables, s.inputs, ev.env, rng)
    dOm, ddOm = Variation(R.get_tangent_vector()), Variation(R.get_tangent_vector()).t_diff()
    dvals = {str(ddOm): pert[str(ddOm)], str(R.get_variation_vector().t_diff()): np.zeros((3, 1))}
    rhs = linearized_rhs(lin, row, ev, pert, dvals)
    J, W = ev.env["J"], ev.env[str(R.get_tangent_vector())]
    expected = (
        -J @ pert[str(ddOm)] + (hat(J @ W) - hat(W) @ J) @ pert[str(dOm)] + pert["\\delta{M}"]
    )
    assert np.allclose(rhs, expected.reshape(-1))


# ---------------------------------------------------------------------------
# Completeness and validation
# ---------------------------------------------------------------------------


class TestCompleteness:
    def test_missing_input_raises_naming_the_term(self):
        """Without inputs=[M] the term δM has no perturbation symbol: no silent zero."""
        s, eom, _ = _derive("free_rigid_body")
        with pytest.raises(AlgebraicError) as err:
            linearize(eom, s.variables)
        assert "\\delta{M}" in str(err.value)
        assert "no perturbation symbol" in str(err.value)

    def test_undeclared_non_constant_symbol_raises(self):
        """A vector in E that is neither constant, configuration nor input is reported."""
        s, eom, _ = _derive("spherical_pendulum_s2")
        q = s.symbols["q"]
        xi = q.get_variation_vector()
        d = Vector("d")  # disturbance without the Constant attribute
        bad = {str(xi): (xi, eom[str(xi)][1] + d)}
        with pytest.raises(AlgebraicError) as err:
            linearize(bad, s.variables, s.inputs)
        assert "\\delta{d}" in str(err.value)

    def test_perturbation_inside_non_canonical_calculus_node_raises(self):
        """A d/dt wrapping a product of δΩ leaves δΩ inside a calculus node."""
        R = SO3("R")
        Om, eta = R.get_tangent_vector(), R.get_variation_vector()
        c = Scalar("c", attr=["Constant"])
        E = TimeDerivative(SVMul(Om, c))
        variables = SystemVariables(matrices=[R])
        with pytest.raises(AlgebraicError) as err:
            linearize({str(eta): (eta, E)}, variables)
        assert "calculus operator" in str(err.value)
        assert str(Variation(Om)) in str(err.value)

    def test_zero_equation_gives_zero_dynamic_row(self):
        """A vanishing equation linearizes to an all-zero row rather than raising."""
        R = SO3("R")
        eta = R.get_variation_vector()
        variables = SystemVariables(matrices=[R])
        lin = linearize({str(eta): (eta, Vector("0", attr=["Constant", "Zero"]))}, variables)
        for table in (lin.M, lin.A, lin.B):
            for block in table[str(eta)].values():
                assert block.is_zero


class TestValidation:
    def test_constant_input_rejected(self):
        s, eom, _ = _derive("free_rigid_body")
        with pytest.raises(ValueError, match="constant"):
            linearize(eom, s.variables, [Vector("F0", attr=["Constant"])])

    def test_non_leaf_input_rejected(self):
        s, eom, _ = _derive("free_rigid_body")
        M = s.symbols["M"]
        with pytest.raises(TypeError, match="Scalar or Vector leaf"):
            linearize(eom, s.variables, [Dot(M, M)])
        with pytest.raises(TypeError, match="Scalar or Vector leaf"):
            linearize(eom, s.variables, [Matrix("K")])

    def test_configuration_variable_as_input_rejected(self):
        s, eom, _ = _derive("point_mass_r3")
        with pytest.raises(ValueError, match="configuration variable"):
            linearize(eom, s.variables, [s.symbols["x"]])

    def test_eom_keys_must_match_variables(self):
        s, eom, _ = _derive("free_rigid_body")
        other = SystemVariables(matrices=[SO3("Q")])
        with pytest.raises(ValueError, match="do not match"):
            linearize(eom, other, s.inputs)

    def test_eom_must_be_dict_of_pairs(self):
        s, eom, _ = _derive("free_rigid_body")
        with pytest.raises(TypeError):
            linearize([eom], s.variables, s.inputs)
        key = next(iter(eom))
        with pytest.raises(TypeError, match="pair"):
            linearize({key: eom[key][1]}, s.variables, s.inputs)

    def test_constant_configuration_variable_rejected(self):
        """Same validation as compute_eom (``_check_variables``)."""
        s, eom, _ = _derive("free_rigid_body")
        with pytest.raises(ValueError, match="constant"):
            linearize(eom, SystemVariables(matrices=[Matrix("K", attr=["Constant"])]))

    def test_general_matrix_configuration_rejected(self):
        K = Matrix("K")
        variables = SystemVariables(matrices=[K])
        with pytest.raises(TypeError, match="SO3"):
            linearize({str(K.delta()): (K.delta(), K)}, variables)

    def test_block_shape_unknown_keys(self):
        _, _, lin = _derive("free_rigid_body")
        with pytest.raises(KeyError):
            lin.block_shape("nope", lin.states[0])
        with pytest.raises(KeyError):
            lin.block_shape(lin.rows[0], "nope")

    def test_inputs_default_to_empty(self):
        """A system without inputs linearizes with empty B rows."""
        R = SO3("R")
        J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
        Om, eta = R.get_tangent_vector(), R.get_variation_vector()
        half = Scalar("0.5", value=0.5, attr=["Constant"])
        variables = SystemVariables(matrices=[R])
        eom = compute_eom(
            Dot(Om, J * Om) * half, Dot(eta, Vector("0", attr=["Constant", "Zero"])), variables
        )
        lin = linearize(eom, variables)
        assert lin.inputs == []
        assert all(lin.B[row] == {} for row in lin.rows)
        assert str(lin.M[str(eta)][str(Variation(Om).t_diff())]) == "J"


# ---------------------------------------------------------------------------
# Numeric helpers
# ---------------------------------------------------------------------------


class TestNumericHelpers:
    def test_expm_so3_is_a_rotation_matching_hat(self):
        rng = np.random.default_rng(1)
        v = rng.normal(size=(3, 1))
        E = expm_so3(v)
        assert np.allclose(E.T @ E, np.eye(3))
        assert abs(np.linalg.det(E) - 1) < 1e-12
        # first-order: exp(ε v̂) ≈ I + ε v̂
        eps = 1e-6
        assert np.allclose((expm_so3(eps * v) - np.eye(3)) / eps, hat(v), atol=1e-5)
        # exact: rotation about v by |v| leaves v fixed
        assert np.allclose(E @ v, v)

    def test_s2_perturbation_stays_on_the_manifold(self):
        s = vbl_systems.spherical_pendulum_s2()
        rng = np.random.default_rng(2)
        env = sample_state(s.variables, s.inputs, rng)
        raw, pert = sample_direction(s.variables, s.inputs, env, rng)
        from tests.numeric_helpers import perturbed_state

        q, w = s.symbols["q"], s.symbols["q"].get_tangent_vector()
        for eps in (1e-3, -1e-3, 0.3):
            env_e = perturbed_state(s.variables, s.inputs, env, raw, eps)
            assert abs(np.linalg.norm(env_e["q"]) - 1) < 1e-12
            assert abs((env_e["q"].T @ env_e[str(w)]).item()) < 1e-12
            assert abs((env_e["q"].T @ env_e[str(w.t_diff())]).item()) < 1e-12
        # induced δq = ξ × q and δω = ξ × ω + w to first order
        eps = 1e-6
        env_p = perturbed_state(s.variables, s.inputs, env, raw, eps)
        env_m = perturbed_state(s.variables, s.inputs, env, raw, -eps)
        xi = pert[str(q.get_variation_vector())]
        dq = (env_p["q"] - env_m["q"]) / (2 * eps)
        assert np.allclose(
            dq, np.cross(xi.reshape(3), env["q"].reshape(3)).reshape(3, 1), atol=1e-8
        )
        dw = (env_p[str(w)] - env_m[str(w)]) / (2 * eps)
        assert np.allclose(dw, pert[str(Variation(w))], atol=1e-8)

    def test_evaluator_requires_values_for_perturbations(self):
        x = Vector("x")
        ev = Evaluator()
        with pytest.raises(KeyError):
            ev(Variation(x))
        with pytest.raises(KeyError):
            ev(TimeDerivative(Variation(x)))
        ev = Evaluator(env={str(Variation(x)): np.ones((3, 1))})
        assert np.allclose(ev(SVMul(Variation(x), 2)), 2 * np.ones((3, 1)))

    def test_evaluator_zero_constants_by_type(self):
        from geomech import Zero, ZeroMatrix, ZeroVector

        ev = Evaluator()
        assert ev(Zero) == 0.0
        assert ev(ZeroVector).shape == (3, 1)
        assert ev(ZeroMatrix).shape == (3, 3)
        assert ev(Vector("0", attr=["Constant", "Zero"])).shape == (3, 1)
