"""Systems used by the variation-based linearization tests.

Each function returns ``System(L, dW, variables, inputs, symbols)`` where ``symbols`` is a
dict of the leaves the tests need (configuration variables, parameters, inputs). The
definitions mirror those pinned in ``tests/test_manifold_pipelines.py`` and
``tests/test_scalar_adversarial.py`` so that the linearizations can be derived by hand from
the pinned equations of motion.
"""

from __future__ import annotations

from typing import NamedTuple

from geomech import (
    S2,
    SO3,
    Cross,
    Dot,
    Matrix,
    Scalar,
    SystemVariables,
    Vector,
    getScalars,
)

HALF = Scalar("0.5", value=0.5, attr=["Constant"])


class System(NamedTuple):
    L: object
    dW: object
    variables: SystemVariables
    inputs: list
    symbols: dict


def point_mass_r3():
    """Point mass in R³: m ẍ = F.  L = ½ m ẋ·ẋ, δW = δx·F."""
    (m,) = getScalars("m", attr=["Constant"])
    x, F = Vector("x"), Vector("F")
    v = x.t_diff()
    L = HALF * m * Dot(v, v)
    dW = Dot(x.delta(), F)
    return System(L, dW, SystemVariables(vectors=[x]), [F], dict(m=m, x=x, F=F))


def spring_mass_r():
    """Spring-mass on R: m s̈ + k s = f (``TestSpringMassPipeline``)."""
    m, k = getScalars("m k", attr=["Constant"])
    x, f = Scalar("x"), Scalar("f")
    v = x.t_diff()
    L = HALF * m * v * v - HALF * k * x * x
    dW = x.delta() * f
    return System(L, dW, SystemVariables(scalars=[x]), [f], dict(m=m, k=k, x=x, f=f))


def free_rigid_body():
    """Free rigid body on SO(3), Lee et al. Eq. 6.16 (``TestFreeRigidBodyPipeline``)."""
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    M = Vector("M")
    R = SO3("R")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = Dot(Om, J * Om) * HALF
    dW = Dot(eta, M)
    return System(L, dW, SystemVariables(matrices=[R]), [M], dict(J=J, M=M, R=R))


def rigid_pendulum_pivot():
    """Rigid pendulum on SO(3), inertia about the pivot, Lee et al. Eq. 6.8
    (``TestRigidPendulumPivotPipeline``)."""
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    rho = Vector("\\rho", attr=["Constant"])
    m, g = getScalars("m g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    M = Vector("M")
    R = SO3("R")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = Dot(Om, J * Om) * HALF - m * g * Dot(R * rho, e3)
    dW = Dot(eta, M)
    symbols = dict(J=J, rho=rho, m=m, g=g, e3=e3, M=M, R=R)
    return System(L, dW, SystemVariables(matrices=[R]), [M], symbols)


def spherical_pendulum_s2():
    """Spherical pendulum on S², ω-form, Lee et al. Eq. 5.14
    (``TestSphericalPendulumPipeline``)."""
    m, g, l = getScalars("m g l", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    q, f = S2("q"), Vector("f")
    omega, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * m * l * l * Dot(omega, omega) - m * g * l * Dot(q, e3)
    dW = Dot(xi, f)
    symbols = dict(m=m, g=g, l=l, e3=e3, q=q, f=f)
    return System(L, dW, SystemVariables(vectors=[q]), [f], symbols)


def quadrotor():
    """Quadrotor on R³ × SO(3) (Lee, Leok, McClamroch CDC 2010; Wu & Sreenath 2015).

    L = ½ m ẋ·ẋ + ½ Ω·JΩ − m g e₃·x,  δW = f (R e₃)·δx + M·η,  inputs [f (Scalar), M].
    """
    m, g = getScalars("m g", attr=["Constant"])
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    e3 = Vector("e3", attr=["Constant"])
    x, R = Vector("x"), SO3("R")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    f, M = Scalar("f"), Vector("M")
    v = x.t_diff()
    L = HALF * m * Dot(v, v) + HALF * Dot(Om, J * Om) - m * g * Dot(e3, x)
    dW = Dot(x.delta(), (R * e3) * f) + Dot(eta, M)
    symbols = dict(m=m, g=g, J=J, e3=e3, x=x, R=R, f=f, M=M)
    return System(L, dW, SystemVariables(vectors=[x], matrices=[R]), [f, M], symbols)


def double_rigid_pendulum():
    """Double rigid pendulum on SO(3) × SO(3) (``TestDoubleRigidPendulumMassMatrixExact``)."""
    J1 = Matrix("J1", attr=["Constant", "SymmetricMatrix"])
    J2 = Matrix("J2", attr=["Constant", "SymmetricMatrix"])
    rho1 = Vector("\\rho_1", attr=["Constant"])
    rho2 = Vector("\\rho_2", attr=["Constant"])
    l1 = Vector("l_1", attr=["Constant"])
    m1, m2, g = getScalars("m1 m2 g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    M1, M2 = Vector("M1"), Vector("M2")
    R1, R2 = SO3("R1"), SO3("R2")
    Om1, Om2 = R1.get_tangent_vector(), R2.get_tangent_vector()
    eta1, eta2 = R1.get_variation_vector(), R2.get_variation_vector()
    x1 = R1 * rho1
    x2 = R1 * l1 + R2 * rho2
    v1, v2 = x1.t_diff(), x2.t_diff()
    KE = (
        Dot(Om1, J1 * Om1) * HALF
        + Dot(Om2, J2 * Om2) * HALF
        + Dot(v1, v1) * m1 * HALF
        + Dot(v2, v2) * m2 * HALF
    )
    PE = m1 * g * Dot(x1, e3) + m2 * g * Dot(x2, e3)
    dW = Dot(eta1, M1) + Dot(eta2, M2)
    symbols = dict(J1=J1, J2=J2, rho1=rho1, rho2=rho2, l1=l1, m1=m1, m2=m2, g=g, R1=R1, R2=R2)
    return System(KE - PE, dW, SystemVariables(matrices=[R1, R2]), [M1, M2], symbols)


def elastic_spherical_pendulum():
    """Elastic spherical pendulum on R × S² (``tests/test_scalar_adversarial.py``).

    KE = ½ m (l̇² + l² ω·ω), PE = m g l (q·e₃) + ½ k (l − l₀)², δW = f_l δl + ξ·f_q.
    """
    m, g, k, l0 = getScalars("m g k l0", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    l, q = Scalar("l"), S2("q")
    omega, xi = q.get_tangent_vector(), q.get_variation_vector()
    f_l, f_q = Scalar("f_l"), Vector("f_q")
    ld = l.t_diff()
    KE = HALF * m * (ld * ld + l * l * Dot(omega, omega))
    PE = m * g * l * Dot(q, e3) + HALF * k * (l - l0) * (l - l0)
    dW = l.delta() * f_l + Dot(xi, f_q)
    symbols = dict(m=m, g=g, k=k, l0=l0, e3=e3, l=l, q=q, f_l=f_l, f_q=f_q)
    return System(KE - PE, dW, SystemVariables(scalars=[l], vectors=[q]), [f_l, f_q], symbols)


def quadrotor_with_load():
    """Quadrotor with a cable-suspended point-mass load on R³ × S² × SO(3)
    (Sreenath, Lee, Kumar CDC 2013).

    Load position x_L, cable direction q from load to quadrotor, cable length l,
    x_Q = x_L − l q.  L = ½ m_Q ẋ_Q·ẋ_Q + ½ m_L ẋ_L·ẋ_L + ½ Ω·JΩ − m_Q g e₃·x_Q − m_L g e₃·x_L,
    δW = f (R e₃)·(δx_L − l ξ × q) + M·η.
    """
    mQ, mL, g, l = getScalars("m_Q m_L g l", attr=["Constant"])
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    e3 = Vector("e3", attr=["Constant"])
    xL, q, R = Vector("x_L"), S2("q"), SO3("R")
    xi = q.get_variation_vector()
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    f, M = Scalar("f"), Vector("M")
    xQ = xL - q * l
    vQ, vL = xQ.t_diff(), xL.t_diff()
    L = (
        HALF * mQ * Dot(vQ, vQ)
        + HALF * mL * Dot(vL, vL)
        + HALF * Dot(Om, J * Om)
        - mQ * g * Dot(e3, xQ)
        - mL * g * Dot(e3, xL)
    )
    dW = Dot(xL.delta() - Cross(xi, q) * l, (R * e3) * f) + Dot(eta, M)
    symbols = dict(mQ=mQ, mL=mL, g=g, l=l, J=J, e3=e3, xL=xL, q=q, R=R, f=f, M=M)
    variables = SystemVariables(vectors=[xL, q], matrices=[R])
    return System(L, dW, variables, [f, M], symbols)
