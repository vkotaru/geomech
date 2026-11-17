"""Systems pinned by the phase-1 engine fixes (S² variation rule, ``collect``) and the
adversarial systems of the phase-1 review, with the helpers the phase-1 tests share.

Each builder returns ``System(L, dW, variables, inputs, symbols)`` like ``tests/vbl_systems``.
They are shared by the pinned-equation tests (``tests/test_s2_variation_fix.py``,
``tests/test_collect_fix.py``), the review tests (``tests/test_phase1_adversarial.py``) and
the numerical check gate (``tests/test_eom_oracle.py``); ``derive`` and ``equation`` are the
"derive and check numerically" helpers of the pinned-equation tests.

Every S² Lagrangian in the pinned sections has ``∂L/∂ω`` with a component along ``q`` or a
velocity nested under a matrix inside a ``Dot``, which is exactly what the pre-fix engine
got wrong.
"""

from __future__ import annotations

from geomech import (
    S2,
    SO3,
    Cross,
    Dot,
    Hat,
    Matrix,
    MVMul,
    Scalar,
    SystemVariables,
    Transpose,
    Vector,
    VVMul,
    compute_eom,
    getScalars,
    systems,
)
from geomech.numeric import assert_eom
from tests.vbl_systems import HALF, System

TOL = 1e-5


def derive(build):
    """``(system, eom, error)``: derive the system and check it with ``check_eom``."""
    s = build()
    eom = compute_eom(s.L, s.dW, s.variables)
    err = assert_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs, rtol=TOL)
    return s, eom, err


def equation(eom, var):
    """The equation of a manifold variable (keyed by its variation vector)."""
    return eom[str(var.get_variation_vector())][1]


# ---------------------------------------------------------------------------
# 1.1  S² variation rule  δω = ξ̇ − 2 ω × ξ
# ---------------------------------------------------------------------------


def s2_gyroscopic():
    """Spherical pendulum with a gyroscopic term: L = ½ m l² ω·ω + c·ω − m g l q·e₃."""
    m, g, l = getScalars("m g l", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    c = Vector("c", attr=["Constant"])
    q, f = S2("q"), Vector("f")
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * m * l * l * Dot(om, om) + Dot(c, om) - m * g * l * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], dict(q=q, c=c, f=f))


def s2_anisotropic_omega():
    """Anisotropic inertia on S² in ω-form: L = ½ ω·Aω − m g l q·e₃, A constant symmetric."""
    m, g, l = getScalars("m g l", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    A = Matrix("A", attr=["Constant", "SymmetricMatrix"])
    q, f = S2("q"), Vector("f")
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * Dot(om, A * om) - m * g * l * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], dict(q=q, A=A, f=f))


def _rotating_sphere_symbols():
    m, L_, Om = getScalars("m L \\Omega", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    C = Matrix("C", attr=["Constant", "SymmetricMatrix"])
    q, M = S2("q"), Vector("M")
    return m, L_, Om, e3, C, q, M


def s2_rotating_sphere_omega():
    """Particle on a sphere rotating about e₃ at rate Ω, Lee et al. §5.6.2, ω-form.

    L̃ = ½ m L² { |ω|² + 2 qᵀ S(Ωe₃) S(q) ω + Ω² qᵀ C q },  C = diag(1, 1, 0) (a constant
    symmetric matrix here).
    """
    s = systems.particle_on_rotating_sphere()  # the same Lagrangian and symbols
    symbols = {k: s.symbols[k] for k in ("q", "C", "Om", "e3", "M")}
    return System(s.L, s.dW, s.variables, list(s.inputs), symbols)


def s2_rotating_sphere_qdot():
    """Same Lagrangian as ``s2_rotating_sphere_omega`` in (q, q̇) form with q̇ = ω × q.

    L = ½ m L² { |q̇|² + 2 qᵀ Sᵀ(Ωe₃) q̇ + Ω² qᵀ C q }.
    """
    m, L_, Om, e3, C, q, M = _rotating_sphere_symbols()
    xi = q.get_variation_vector()
    qd = q.t_diff()
    S_a = Hat(e3 * Om)
    coupling = Dot(q, MVMul(Transpose(S_a), qd)) * 2
    L = HALF * m * L_ * L_ * (Dot(qd, qd) + coupling + Om * Om * Dot(q, C * q))
    return System(
        L, Dot(xi, M), SystemVariables(vectors=[q]), [M], dict(q=q, C=C, Om=Om, e3=e3, M=M)
    )


def s2_rotating_frame():
    """Spherical pendulum in a frame rotating at W₀ (dynamics audit F1, row 2).

    L = ½ m l² ( |ω + W₀|² − (W₀·q)² ) − m g l q·e₃.
    """
    m, g, l = getScalars("m g l", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    W0 = Vector("W0", attr=["Constant"])
    q, f = S2("q"), Vector("f")
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    w = om + W0
    L = HALF * m * l * l * (Dot(w, w) - Dot(W0, q) * Dot(W0, q)) - m * g * l * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], dict(q=q, W0=W0, f=f))


def r3_s2_velocity_coupling():
    """R³ × S² with a velocity coupling: L = ½ m ẋ·ẋ + ½ m l² ω·ω + ẋ·ω."""
    m, l = getScalars("m l", attr=["Constant"])
    x, F = Vector("x"), Vector("F")
    q, f = S2("q"), Vector("f")
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    xd = x.t_diff()
    L = HALF * m * Dot(xd, xd) + HALF * m * l * l * Dot(om, om) + Dot(xd, om)
    dW = Dot(x.delta(), F) + Dot(xi, f)
    return System(L, dW, SystemVariables(vectors=[x, q]), [F, f], dict(x=x, q=q, F=F, f=f))


def _s2_so3_symbols():
    m, k = getScalars("m k", attr=["Constant"])
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    e3 = Vector("e3", attr=["Constant"])
    q, f = S2("q"), Vector("f")
    R, M = SO3("R"), Vector("M")
    return m, k, J, e3, q, f, R, M


def s2_so3_omega_coupling():
    """S² × SO(3) with L = ½ m ω·ω + ½ Ω·JΩ + k ω·Ω."""
    m, k, J, e3, q, f, R, M = _s2_so3_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = HALF * m * Dot(om, om) + HALF * Dot(Om, J * Om) + k * Dot(om, Om)
    dW = Dot(xi, f) + Dot(eta, M)
    return System(
        L, dW, SystemVariables(vectors=[q], matrices=[R]), [f, M], dict(q=q, R=R, f=f, M=M)
    )


def s2_so3_axis_coupling():
    """S² × SO(3) with L = ½ m ω·ω + ½ Ω·JΩ + k ω·(Rᵀe₃)."""
    m, k, J, e3, q, f, R, M = _s2_so3_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = HALF * m * Dot(om, om) + HALF * Dot(Om, J * Om) + k * Dot(om, MVMul(Transpose(R), e3))
    dW = Dot(xi, f) + Dot(eta, M)
    return System(
        L, dW, SystemVariables(vectors=[q], matrices=[R]), [f, M], dict(q=q, R=R, f=f, M=M)
    )


# ---------------------------------------------------------------------------
# 1.2  collect: peel a matrix off the side holding the target
# ---------------------------------------------------------------------------


def _so3_body_symbols():
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    b = Vector("b", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    R, M = SO3("R"), Vector("M")
    return J, b, e3, R, M


def so3_gyro_coupling_J_on_velocity():
    """SO(3) body: L = ½ Ω·JΩ + (JΩ)·((Rᵀe₃) × b)   (J on the velocity side of the Dot)."""
    J, b, e3, R, M = _so3_body_symbols()
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = HALF * Dot(Om, J * Om) + Dot(J * Om, Cross(MVMul(Transpose(R), e3), b))
    return System(L, Dot(eta, M), SystemVariables(matrices=[R]), [M], dict(R=R, J=J, M=M))


def so3_gyro_coupling_J_on_axis():
    """Same physics as ``so3_gyro_coupling_J_on_velocity``: L = ½ Ω·JΩ + Ω·J((Rᵀe₃) × b)."""
    J, b, e3, R, M = _so3_body_symbols()
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = HALF * Dot(Om, J * Om) + Dot(Om, J * Cross(MVMul(Transpose(R), e3), b))
    return System(L, Dot(eta, M), SystemVariables(matrices=[R]), [M], dict(R=R, J=J, M=M))


def so3_inertial_momentum_R_J_Omega():
    """SO(3) body with the inertial angular momentum: L = ½ Ω·JΩ + (R(JΩ))·e₃."""
    J, b, e3, R, M = _so3_body_symbols()
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = HALF * Dot(Om, J * Om) + Dot(MVMul(R, MVMul(J, Om)), e3)
    return System(L, Dot(eta, M), SystemVariables(matrices=[R]), [M], dict(R=R, J=J, M=M))


def so3_inertial_momentum_RJ_Omega():
    """Same as ``so3_inertial_momentum_R_J_Omega`` with the product grouped as ((RJ)Ω)·e₃."""
    J, b, e3, R, M = _so3_body_symbols()
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = HALF * Dot(Om, J * Om) + Dot(MVMul(R * J, Om), e3)
    return System(L, Dot(eta, M), SystemVariables(matrices=[R]), [M], dict(R=R, J=J, M=M))


def s2_anisotropic_qdot():
    """Anisotropic inertia on S² in q̇-form: L = ½ q̇·A q̇ − m g q·e₃ with q̇ = ω × q."""
    m, g = getScalars("m g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    A = Matrix("A", attr=["Constant", "SymmetricMatrix"])
    q, f = S2("q"), Vector("f")
    xi = q.get_variation_vector()
    qd = q.t_diff()
    L = HALF * Dot(qd, A * qd) - m * g * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], dict(q=q, A=A, f=f))


def so3_so3_relative_rate_coupling():
    """SO(3) × SO(3): L = ½ Ω₁·J₁Ω₁ + ½ Ω₂·J₂Ω₂ + k Ω₁·(R₁ᵀ(R₂Ω₂))."""
    (k,) = getScalars("k", attr=["Constant"])
    J1 = Matrix("J1", attr=["Constant", "SymmetricMatrix"])
    J2 = Matrix("J2", attr=["Constant", "SymmetricMatrix"])
    R1, R2 = SO3("R1"), SO3("R2")
    M1, M2 = Vector("M1"), Vector("M2")
    Om1, eta1 = R1.get_tangent_vector(), R1.get_variation_vector()
    Om2, eta2 = R2.get_tangent_vector(), R2.get_variation_vector()
    rel = MVMul(Transpose(R1), MVMul(R2, Om2))
    L = HALF * Dot(Om1, J1 * Om1) + HALF * Dot(Om2, J2 * Om2) + k * Dot(Om1, rel)
    dW = Dot(eta1, M1) + Dot(eta2, M2)
    return System(
        L, dW, SystemVariables(matrices=[R1, R2]), [M1, M2], dict(R1=R1, R2=R2, M1=M1, M2=M2)
    )


def so3_offset_mass_matrix_inertia():
    """SO(3) with the velocity under a Cross under a matrix: L = ½ (Ω×ρ)·J(Ω×ρ)."""
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    rho = Vector("\\rho", attr=["Constant"])
    R, M = SO3("R"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    v = Cross(Om, rho)
    L = HALF * Dot(v, J * v)
    return System(L, Dot(eta, M), SystemVariables(matrices=[R]), [M], dict(R=R, J=J, M=M))


S2_VARIATION_SYSTEMS = [
    s2_gyroscopic,
    s2_anisotropic_omega,
    s2_rotating_sphere_omega,
    s2_rotating_sphere_qdot,
    s2_rotating_frame,
    r3_s2_velocity_coupling,
    s2_so3_omega_coupling,
    s2_so3_axis_coupling,
]

COLLECT_SYSTEMS = [
    so3_gyro_coupling_J_on_velocity,
    so3_gyro_coupling_J_on_axis,
    so3_inertial_momentum_R_J_Omega,
    so3_inertial_momentum_RJ_Omega,
    s2_anisotropic_qdot,
    so3_so3_relative_rate_coupling,
    so3_offset_mass_matrix_inertia,
]


# ---------------------------------------------------------------------------
# Adversarial systems of the phase-1 review (``tests/test_phase1_adversarial.py``)
# ---------------------------------------------------------------------------


def s2_symbols():
    """``(m, g, l, k, e3, a, c, d, A, S, J, q, f)`` shared by the adversarial S² systems."""
    m, g, l, k = (Scalar(n, attr=["Constant"]) for n in ("m", "g", "l", "k"))
    e3, a, c, d = (Vector(n, attr=["Constant"]) for n in ("e3", "a", "c", "d"))
    A = Matrix("A", attr=["Constant"])
    S = Matrix("S", attr=["Constant", "SymmetricMatrix"])
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    q, f = S2("q"), Vector("f")
    return m, g, l, k, e3, a, c, d, A, S, J, q, f


def adv_q_dependent_gyro():
    """∂L/∂ω = m l² ω + (q·a) c: normal component and explicit q-dependence."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * m * l * l * Dot(om, om) + Dot(q, a) * Dot(c, om) - m * g * l * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], dict(q=q))


def adv_bilinear_omega():
    """L = ½ ω·ω + (c·ω)(d·ω) − g q·e₃."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * Dot(om, om) + Dot(c, om) * Dot(d, om) - g * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], dict(q=q))


def adv_omega_dot_Aq():
    """∂L/∂ω = A q with A a general (non-symmetric) constant matrix."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * Dot(om, om) + Dot(om, A * q) - g * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], dict(q=q))


def adv_Aomega_dot_q():
    """Same physics as ``adv_omega_dot_Aq`` with the matrix on the ω side of the Dot."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * Dot(om, om) + Dot(A * om, q) - g * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], dict(q=q))


def adv_omega_and_qdot_mixed():
    """ω and q̇ = ω × q in one product: L = ½ q̇·q̇ + c·ω + (q̇·d)(ω·e₃) − g q·e₃."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi, qd = q.get_tangent_vector(), q.get_variation_vector(), q.t_diff()
    L = HALF * Dot(qd, qd) + Dot(c, om) + Dot(qd, d) * Dot(om, e3) - g * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], dict(q=q))


def adv_omega_dot_qdot_cross():
    """L = ½ q̇·q̇ + ω·(q̇ × c) − g q·e₃."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi, qd = q.get_tangent_vector(), q.get_variation_vector(), q.t_diff()
    L = HALF * Dot(qd, qd) + Dot(om, Cross(qd, c)) - g * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], dict(q=q))


def adv_omega_dot_S_qdot():
    """L = ½ q̇·q̇ + ω·(S q̇) − g q·e₃ (velocity under a matrix, ω outside)."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi, qd = q.get_tangent_vector(), q.get_variation_vector(), q.t_diff()
    L = HALF * Dot(qd, qd) + Dot(om, S * qd) - g * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], dict(q=q))


def adv_virtual_work_under_rotation():
    """δW = (R ξ)·f with R a frozen SO(3) parameter: E gets Rᵀ f."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    R = SO3("R")
    L = HALF * m * l * l * Dot(om, om) - m * g * l * Dot(q, e3)
    return System(L, Dot(MVMul(R, xi), f), SystemVariables(vectors=[q]), [f], dict(q=q))


def adv_virtual_work_hat_xi():
    """δW = (ξ̂ q)·f = (ξ × q)·f, a force on the point x = q."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * Dot(om, om) - g * Dot(q, e3)
    return System(L, Dot(MVMul(Hat(xi), q), f), SystemVariables(vectors=[q]), [f], dict(q=q))


def adv_virtual_work_symmetric_matrix():
    """δW = (S ξ)·f."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * Dot(om, om) - g * Dot(q, e3)
    return System(L, Dot(MVMul(S, xi), f), SystemVariables(vectors=[q]), [f], dict(q=q))


def adv_two_spheres_shared_parameter():
    """Two S² variables sharing c and k:

    L = ½|ω₁|² + ½|ω₂|² + c·(ω₁+ω₂) + k ω₁·ω₂ − g (q₁+q₂)·e₃.
    """
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    q1, q2, f1, f2 = S2("q1"), S2("q2"), Vector("f1"), Vector("f2")
    w1, x1 = q1.get_tangent_vector(), q1.get_variation_vector()
    w2, x2 = q2.get_tangent_vector(), q2.get_variation_vector()
    L = (
        HALF * Dot(w1, w1)
        + HALF * Dot(w2, w2)
        + Dot(c, w1 + w2)
        + k * Dot(w1, w2)
        - g * Dot(q1 + q2, e3)
    )
    dW = Dot(x1, f1) + Dot(x2, f2)
    return System(L, dW, SystemVariables(vectors=[q1, q2]), [f1, f2], dict(q1=q1, q2=q2))


def adv_two_spheres_mixed_forms():
    """L = ½|ω₁|² + ½|ω₂|² + k q̇₁·q̇₂ + (q₁·q₂)(c·ω₁)."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    q1, q2, f1, f2 = S2("q1"), S2("q2"), Vector("f1"), Vector("f2")
    w1, x1 = q1.get_tangent_vector(), q1.get_variation_vector()
    w2, x2 = q2.get_tangent_vector(), q2.get_variation_vector()
    L = (
        HALF * Dot(w1, w1)
        + HALF * Dot(w2, w2)
        + k * Dot(q1.t_diff(), q2.t_diff())
        + Dot(q1, q2) * Dot(c, w1)
    )
    dW = Dot(x1, f1) + Dot(x2, f2)
    return System(L, dW, SystemVariables(vectors=[q1, q2]), [f1, f2], dict(q1=q1, q2=q2))


def adv_two_spheres_cross_coupling():
    """L = ½|ω₁|² + ½|ω₂|² + k ω₁·(q₂ × ω₂)."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    q1, q2, f1, f2 = S2("q1"), S2("q2"), Vector("f1"), Vector("f2")
    w1, x1 = q1.get_tangent_vector(), q1.get_variation_vector()
    w2, x2 = q2.get_tangent_vector(), q2.get_variation_vector()
    L = HALF * Dot(w1, w1) + HALF * Dot(w2, w2) + k * Dot(w1, Cross(q2, w2))
    dW = Dot(x1, f1) + Dot(x2, f2)
    return System(L, dW, SystemVariables(vectors=[q1, q2]), [f1, f2], dict(q1=q1, q2=q2))


def adv_s2_so3_R_omega():
    """S² × SO(3): L = ½ ω·ω + ½ Ω·JΩ + k (R ω)·e₃ (R on the ω side of the Dot)."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    R, M = SO3("R"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = HALF * Dot(om, om) + HALF * Dot(Om, J * Om) + k * Dot(MVMul(R, om), e3)
    dW = Dot(xi, f) + Dot(eta, M)
    return System(L, dW, SystemVariables(vectors=[q], matrices=[R]), [f, M], dict(q=q, R=R))


def adv_s2_so3_Omega_Rq():
    """S² × SO(3): L = ½ ω·ω + ½ Ω·JΩ + k Ω·(R q)."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    R, M = SO3("R"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = HALF * Dot(om, om) + HALF * Dot(Om, J * Om) + k * Dot(Om, MVMul(R, q))
    dW = Dot(xi, f) + Dot(eta, M)
    return System(L, dW, SystemVariables(vectors=[q], matrices=[R]), [f, M], dict(q=q, R=R))


def adv_s2_so3_omega_Omega_cross_q():
    """S² × SO(3): L = ½ ω·ω + ½ Ω·JΩ + k ω·(Ω × q)."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    R, M = SO3("R"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = HALF * Dot(om, om) + HALF * Dot(Om, J * Om) + k * Dot(om, Cross(Om, q))
    dW = Dot(xi, f) + Dot(eta, M)
    return System(L, dW, SystemVariables(vectors=[q], matrices=[R]), [f, M], dict(q=q, R=R))


def adv_r3_s2_q_dependent_velocity_coupling():
    """R³ × S²: L = ½ m ẋ·ẋ + ½ ω·ω + (q·a)(ẋ·ω)."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    x, F = Vector("x"), Vector("F")
    xd = x.t_diff()
    L = HALF * m * Dot(xd, xd) + HALF * Dot(om, om) + Dot(q, a) * Dot(xd, om)
    dW = Dot(x.delta(), F) + Dot(xi, f)
    return System(L, dW, SystemVariables(vectors=[x, q]), [F, f], dict(x=x, q=q))


def adv_r3_s2_pendulum_with_gyro():
    """R³ × S²: L = ½ m |ẋ + l q̇|² + ω·(x × c)."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi, qd = q.get_tangent_vector(), q.get_variation_vector(), q.t_diff()
    x, F = Vector("x"), Vector("F")
    xd = x.t_diff()
    L = HALF * m * Dot(xd + l * qd, xd + l * qd) + Dot(om, Cross(x, c))
    dW = Dot(x.delta(), F) + Dot(xi, f)
    return System(L, dW, SystemVariables(vectors=[x, q]), [F, f], dict(x=x, q=q))


def adv_r3_s2_so3_triple():
    """R³ × S² × SO(3): L = ½ m ẋ·ẋ + ½ ω·ω + ½ Ω·JΩ + (Rᵀẋ)·ω + k Ω·ω."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    x, F = Vector("x"), Vector("F")
    xd = x.t_diff()
    R, M = SO3("R"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = (
        HALF * m * Dot(xd, xd)
        + HALF * Dot(om, om)
        + HALF * Dot(Om, J * Om)
        + Dot(MVMul(Transpose(R), xd), om)
        + k * Dot(Om, om)
    )
    dW = Dot(x.delta(), F) + Dot(xi, f) + Dot(eta, M)
    return System(
        L, dW, SystemVariables(vectors=[x, q], matrices=[R]), [F, f, M], dict(x=x, q=q, R=R)
    )


def adv_r_s2_scalar_length():
    """R × S²: L = ½ ṡ² + ½ s² ω·ω + s c·ω − g s q·e₃ (a scalar variable scales the gyro term)."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    s, fs = Scalar("s"), Scalar("f_s")
    sd = s.t_diff()
    L = HALF * sd * sd + HALF * s * s * Dot(om, om) + s * Dot(c, om) - g * s * Dot(q, e3)
    dW = s.delta() * fs + Dot(xi, f)
    return System(L, dW, SystemVariables(scalars=[s], vectors=[q]), [fs, f], dict(s=s, q=q))


def adv_hat_omega_in_lagrangian_s2():
    """L = ½ ω·ω + a·(ω̂ c) − g q·e₃: the tangent vector inside a Hat."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * Dot(om, om) + Dot(a, MVMul(Hat(om), c)) - g * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], dict(q=q))


def adv_hat_Omega_times_matrix_so3():
    """SO(3): L = ½ Ω·JΩ + a·((Ω̂ J) c) and + a·((Ω̂ + J) c) (Hat of the velocity in products)."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    R, M = SO3("R"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = HALF * Dot(Om, J * Om) + Dot(a, MVMul(Hat(Om) * J, c)) + Dot(a, MVMul(Hat(Om) + J, c))
    return System(L, Dot(eta, M), SystemVariables(matrices=[R]), [M], dict(R=R))


ADVERSARIAL_SYSTEMS = [
    adv_q_dependent_gyro,
    adv_bilinear_omega,
    adv_omega_dot_Aq,
    adv_Aomega_dot_q,
    adv_omega_and_qdot_mixed,
    adv_omega_dot_qdot_cross,
    adv_omega_dot_S_qdot,
    adv_virtual_work_under_rotation,
    adv_virtual_work_hat_xi,
    adv_virtual_work_symmetric_matrix,
    adv_two_spheres_shared_parameter,
    adv_two_spheres_mixed_forms,
    adv_two_spheres_cross_coupling,
    adv_s2_so3_R_omega,
    adv_s2_so3_Omega_Rq,
    adv_s2_so3_omega_Omega_cross_q,
    adv_r3_s2_q_dependent_velocity_coupling,
    adv_r3_s2_pendulum_with_gyro,
    adv_r3_s2_so3_triple,
    adv_r_s2_scalar_length,
    adv_hat_omega_in_lagrangian_s2,
    adv_hat_Omega_times_matrix_so3,
]


# ---------------------------------------------------------------------------
# Reproducers that are not numerical-check systems
# ---------------------------------------------------------------------------


def stiff_point_mass(k_value=None):
    """Point mass on a spring: L = ½ m ẋ·ẋ − ½ k x·x − m g x·e₃, δW = δx·f.

    With ``k_value`` the stiffness is a numeric constant (``Scalar("k", value=k_value)``)
    that the engine folds and the evaluator honours; ``k = 1e7`` is the review's F2
    reproducer (a term of order one next to one of order 1e7).
    """
    m, g = Scalar("m", attr=["Constant"]), Scalar("g", attr=["Constant"])
    k = (
        Scalar("k", value=k_value, attr=["Constant"])
        if k_value
        else Scalar("k", attr=["Constant"])
    )
    e3, x, f = Vector("e3", attr=["Constant"]), Vector("x"), Vector("f")
    v = x.t_diff()
    L = HALF * m * Dot(v, v) - HALF * k * Dot(x, x) - m * g * Dot(x, e3)
    return System(L, Dot(x.delta(), f), SystemVariables(vectors=[x]), [f], dict(x=x, k=k))


def so3_outer_product_on_velocity():
    """SO(3): L = ½ Ω·JΩ + a·((Ω cᵀ) b), i.e. (a·Ω)(c·b) written with an outer product.

    Legal, but the outer product survives ``simplify`` and puts the IBP target inside the
    matrix factor of an ``MVMul``; ``compute_eom`` raises ``AlgebraicError`` from
    ``collect`` (review F1: it used to recurse forever).
    """
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    a, b, c = (Vector(n, attr=["Constant"]) for n in "abc")
    R, M = SO3("R"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = HALF * Dot(Om, J * Om) + Dot(a, MVMul(VVMul(Om, Transpose(c)), b))
    return System(L, Dot(eta, M), SystemVariables(matrices=[R]), [M], dict(R=R))


def s2_outer_product_on_velocity():
    """S² twin of ``so3_outer_product_on_velocity``: L = ½ ω·ω + a·((ω cᵀ) d) − g q·e₃."""
    m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * Dot(om, om) + Dot(a, MVMul(VVMul(om, Transpose(c)), d)) - g * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], dict(q=q))
