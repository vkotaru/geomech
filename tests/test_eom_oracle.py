"""Every pinned system passes the finite-difference Euler-Lagrange oracle.

The systems are the nine of ``tests/vbl_systems.py`` plus the ones pinned in
``tests/test_manifold_pipelines.py``, ``tests/test_standard_form.py`` and
``tests/test_scalar_adversarial.py`` that are not already there, rebuilt here from the
same Lagrangians, and the systems of ``tests/phase1_systems.py`` pinned by the phase-1
engine fixes. ``check_eom`` compares each equation of ``compute_eom`` with a residual
derived without the library's variational rules (see ``docs/numeric-verification.md``), so
this file is the regression gate for every change to the engine.
"""

import pytest

from geomech import (
    S2,
    SO3,
    Dot,
    Matrix,
    Scalar,
    SystemVariables,
    Vector,
    check_eom,
    compute_eom,
    getScalars,
)
from tests import phase1_systems, vbl_systems
from tests.vbl_systems import System

HALF = Scalar("0.5", value=0.5, attr=["Constant"])
TOL = 1e-5


# ---------------------------------------------------------------------------
# Systems pinned elsewhere that vbl_systems does not have
# ---------------------------------------------------------------------------


def spherical_pendulum_qdot():
    """Spherical pendulum, L(q, q̇) form with x = l q (``TestSphericalPendulumQdotPipeline``)."""
    m, g, l = getScalars("m g l", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    q, f = S2("q"), Vector("f")
    v = (l * q).t_diff()
    L = m * Dot(v, v) * HALF - m * g * l * Dot(q, e3)
    dW = Dot(q.get_variation_vector(), f)
    return System(L, dW, SystemVariables(vectors=[q]), [f], {})


def rigid_pendulum_com():
    """Rigid pendulum with J about the centre of mass (``TestRigidPendulumCOMPipeline``)."""
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    rho = Vector("\\rho", attr=["Constant"])
    m, g = getScalars("m g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    R, M = SO3("R"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    x = R * rho
    v = x.t_diff()
    L = Dot(Om, J * Om) * HALF + Dot(v, v) * m * HALF - m * g * Dot(x, e3)
    return System(L, Dot(eta, M), SystemVariables(matrices=[R]), [M], {})


def point_mass_and_rigid_body():
    """Uncoupled R³ × SO(3) system (``TestMixedR3SO3Pipeline``)."""
    (m,) = getScalars("m", attr=["Constant"])
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    x, R, F, M = Vector("x"), SO3("R"), Vector("F"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    v = x.t_diff()
    L = m * Dot(v, v) * HALF + Dot(Om, J * Om) * HALF
    dW = Dot(x.delta(), F) + Dot(eta, M)
    return System(L, dW, SystemVariables(vectors=[x], matrices=[R]), [F, M], {})


def point_mass_gravity():
    """Point mass under gravity (``TestPointMassStandardForm``)."""
    m, g = getScalars("m g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    x, f = Vector("x"), Vector("f")
    v = x.t_diff()
    L = m * Dot(v, v) * HALF - m * Dot(x, g * e3)
    return System(L, Dot(x.delta(), f), SystemVariables(vectors=[x]), [f], {})


def free_particle():
    """Point mass with no forcing, δW the zero constant (``TestNoInputs``)."""
    (m,) = getScalars("m", attr=["Constant"])
    x = Vector("x")
    v = x.t_diff()
    L = m * Dot(v, v) * HALF
    dW = Scalar("0", value=0, attr=["Constant", "Zero"])
    return System(L, dW, SystemVariables(vectors=[x]), [], {})


def spring_mass_damper():
    """Spring-mass-damper, damping through δW (``TestSpringMassDamperPipeline``)."""
    m, k, c = getScalars("m k c", attr=["Constant"])
    x, f = Scalar("x"), Scalar("f")
    v = x.t_diff()
    L = HALF * m * v * v - HALF * k * x * x
    dW = x.delta() * f - c * v * x.delta()
    return System(L, dW, SystemVariables(scalars=[x]), [f], {})


def double_spring_mass():
    """Two coupled spring-masses on R² (``TestDoubleSpringMassPipeline``)."""
    m1, m2, k1, k2 = getScalars("m1 m2 k1 k2", attr=["Constant"])
    x1, x2, f1, f2 = Scalar("x1"), Scalar("x2"), Scalar("f1"), Scalar("f2")
    v1, v2 = x1.t_diff(), x2.t_diff()
    KE = HALF * m1 * v1 * v1 + HALF * m2 * v2 * v2
    PE = HALF * k1 * x1 * x1 + HALF * k2 * (x2 - x1) * (x2 - x1)
    dW = x1.delta() * f1 + x2.delta() * f2
    return System(KE - PE, dW, SystemVariables(scalars=[x1, x2]), [f1, f2], {})


def elastic_spherical_pendulum_qdot():
    """Elastic spherical pendulum, L(q, q̇) form (``TestElasticSphericalPendulumQdotForm``)."""
    m, g, k, l0 = getScalars("m g k l0", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    l, q = Scalar("l"), S2("q")
    f_l, f_q = Scalar("f_l"), Vector("f_q")
    v = (l * q).t_diff()
    KE = HALF * m * Dot(v, v)
    PE = m * g * l * Dot(q, e3) + HALF * k * (l - l0) * (l - l0)
    dW = l.delta() * f_l + Dot(q.get_variation_vector(), f_q)
    return System(KE - PE, dW, SystemVariables(scalars=[l], vectors=[q]), [f_l, f_q], {})


def prismatic_slide():
    """Point mass on a prismatic slide, R³ × R with an inertial cross term
    (``TestR3xRCoupledMassMatrix``)."""
    (m,) = getScalars("m", attr=["Constant"])
    e1 = Vector("e1", attr=["Constant", "UnitNorm"])
    x, l, F, f_l = Vector("x"), Scalar("l"), Vector("F"), Scalar("f_l")
    v = x.t_diff() + e1 * l.t_diff()
    L = HALF * m * Dot(v, v)
    dW = Dot(x.delta(), F) + l.delta() * f_l
    return System(L, dW, SystemVariables(scalars=[l], vectors=[x]), [F, f_l], {})


def varying_inertia():
    """L = ½ (m₀ + m₁ x) ẋ² (``TestVaryingInertia``)."""
    m0, m1 = getScalars("m0 m1", attr=["Constant"])
    x, f = Scalar("x"), Scalar("f")
    v = x.t_diff()
    L = HALF * (m0 + m1 * x) * v * v
    return System(L, x.delta() * f, SystemVariables(scalars=[x]), [f], {})


def squared_inertia():
    """L = ½ m x² ẋ² (``TestVaryingInertia.test_x_squared_inertia_exact``)."""
    (m,) = getScalars("m", attr=["Constant"])
    x, f = Scalar("x"), Scalar("f")
    v = x.t_diff()
    return System(HALF * m * x * x * v * v, x.delta() * f, SystemVariables(scalars=[x]), [f], {})


def sliding_mass_on_rotating_body():
    """Sliding mass on a rotating body, SO(3) × R (``TestStandardFormHygiene``)."""
    (m,) = getScalars("m", attr=["Constant"])
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    rho = Vector("rho", attr=["Constant", "UnitNorm"])
    R, l = SO3("R"), Scalar("l")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    f_l, M = Scalar("f_l"), Vector("M")
    v = (R * (rho * l)).t_diff()
    L = HALF * m * Dot(v, v) + HALF * Dot(Om, J * Om)
    dW = l.delta() * f_l + Dot(eta, M)
    return System(L, dW, SystemVariables(scalars=[l], matrices=[R]), [f_l, M], {})


# ---------------------------------------------------------------------------
# Row vectors built with the operators (phase 3 review, F2): a row times a matrix
# ---------------------------------------------------------------------------


def _row_symbols():
    (m,) = getScalars("m", attr=["Constant"])
    K = Matrix("K", attr=["Constant", "SymmetricMatrix"])
    A = Matrix("A", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    x, f = Vector("x"), Vector("f")
    return m, K, A, e3, x, f


def row_spring():
    """L = ½ m ẋ·ẋ − ½ xᵀK x (K symmetric): m ẍ + K x = f."""
    m, K, A, e3, x, f = _row_symbols()
    v = x.t_diff()
    L = m * Dot(v, v) / 2 - x.T() * K * x / 2
    return System(L, Dot(x.delta(), f), SystemVariables(vectors=[x]), [f], dict(m=m, K=K, x=x))


def row_quadratic_form():
    """L = ½ m ẋ·ẋ − xᵀA x (A not symmetric): m ẍ + (A + Aᵀ) x = f."""
    m, K, A, e3, x, f = _row_symbols()
    v = x.t_diff()
    L = m * Dot(v, v) / 2 - x.T() * A * x
    return System(L, Dot(x.delta(), f), SystemVariables(vectors=[x]), [f], dict(m=m, A=A, x=x))


def row_shifted_spring():
    """L = ½ m ẋ·ẋ − (x − e₃)ᵀK(x − e₃) (K symmetric): m ẍ + 2K(x − e₃) = f."""
    m, K, A, e3, x, f = _row_symbols()
    v = x.t_diff()
    L = m * Dot(v, v) / 2 - (x - e3).T() * K * (x - e3)
    return System(L, Dot(x.delta(), f), SystemVariables(vectors=[x]), [f], dict(m=m, K=K, x=x))


def row_damping_work():
    """L = ½ m ẋ·ẋ, δW = −δxᵀA ẋ + δx·f: m ẍ + A ẋ = f."""
    m, K, A, e3, x, f = _row_symbols()
    v = x.t_diff()
    L = m * Dot(v, v) / 2
    dW = -(x.delta().T() * A * v) + Dot(x.delta(), f)
    return System(L, dW, SystemVariables(vectors=[x]), [f], dict(m=m, A=A, x=x))


def row_kinetic_energy():
    """L = ½ ẋᵀA ẋ − ½ xᵀK x (A not symmetric, K symmetric): ½(A + Aᵀ) ẍ + K x = f.

    The kinetic energy is written with rows; ``collect`` turns the inner products that
    hold d/dt(δx) into ``Dot`` so that integration by parts reaches them (phase 4.2).
    """
    m, K, A, e3, x, f = _row_symbols()
    v = x.t_diff()
    L = v.T() * A * v / 2 - x.T() * K * x / 2
    return System(L, Dot(x.delta(), f), SystemVariables(vectors=[x]), [f], dict(A=A, K=K, x=x))


def row_kinetic_energy_so3():
    """Free rigid body with a row kinetic energy, L = ½ Ωᵀ J Ω: Lee et al. Eq. 6.16."""
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    R, M = SO3("R"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = Om.T() * J * Om / 2
    return System(L, Dot(eta, M), SystemVariables(matrices=[R]), [M], dict(J=J, R=R, M=M))


SYSTEMS = [
    vbl_systems.point_mass_r3,
    vbl_systems.spring_mass_r,
    vbl_systems.free_rigid_body,
    vbl_systems.rigid_pendulum_pivot,
    vbl_systems.spherical_pendulum_s2,
    vbl_systems.quadrotor,
    vbl_systems.double_rigid_pendulum,
    vbl_systems.elastic_spherical_pendulum,
    vbl_systems.quadrotor_with_load,
    spherical_pendulum_qdot,
    rigid_pendulum_com,
    point_mass_and_rigid_body,
    point_mass_gravity,
    free_particle,
    spring_mass_damper,
    double_spring_mass,
    elastic_spherical_pendulum_qdot,
    prismatic_slide,
    varying_inertia,
    squared_inertia,
    sliding_mass_on_rotating_body,
    *phase1_systems.S2_VARIATION_SYSTEMS,
    *phase1_systems.COLLECT_SYSTEMS,
    row_spring,
    row_quadratic_form,
    row_shifted_spring,
    row_damping_work,
    row_kinetic_energy,
    row_kinetic_energy_so3,
]


@pytest.mark.parametrize("build", SYSTEMS, ids=lambda f: f.__name__)
def test_eom_matches_euler_lagrange_residual(build):
    s = build()
    eom = compute_eom(s.L, s.dW, s.variables)
    err = check_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs)
    assert err < TOL, f"{build.__name__}: worst relative error {err:.2e}"
