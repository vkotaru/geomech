"""Multibody systems on SO(3) and its products: a spacecraft with reaction wheels (wheel
angles on R³ or on (S¹)³) and the swing phase of a three-dimensional compass-gait biped.

Ω is the body-frame angular velocity (``Ṙ = R Ω̂``), η the variation (``δR = R η̂``); a
body-frame moment ``M`` on ``R`` enters as ``η·M`` and an inertial-frame moment ``M`` as
``η·(RᵀM)``. Gravity acts along −e₃. Each docstring states the Lagrangian and the virtual
work.
"""

from __future__ import annotations

import numpy as np

from geomech.core.base import SO3, I, Matrix, Scalar, Vector, ZeroVector, getScalars
from geomech.core.operations import Dot, MVMul, Transpose
from geomech.dynamics.variables import SystemVariables
from geomech.systems._base import E3, HALF, Reference, frozen_array, make_system, register, rot

LEE = "lee2018global"

#: Bus inertia of the spacecraft about its centre of mass (kg m²), including the wheels'
#: masses and transverse inertias; products of inertia are not zero.
J_BUS = frozen_array([[1.2, 0.02, 0.01], [0.02, 1.0, 0.03], [0.01, 0.03, 0.9]])
#: Inertia of a compass-gait leg (5 kg, 1 m slender rod) about its centre of mass (kg m²).
J_LEG = frozen_array(np.diag([0.42, 0.42, 0.005]))


def _inertia(name):
    return Matrix(name, attr=["Constant", "SymmetricMatrix"])


def _spacecraft_scene(R, bs, markers=()):
    """The bus (a box with two solar panels) and a cylinder per wheel centred at ``0.3 bᵢ``
    with its disc normal along the body axis e₁, e₂ or e₃; *markers* are body-frame points
    drawn as small spheres (the rim mark of an S¹ wheel)."""

    def scene(layout, s):
        from geomech.viz import Box, Cylinder, PointMass, RigidBody, Scene

        sc = Scene(layout)
        sc.add(RigidBody(R, ZeroVector, Box(0.5, 0.4, 0.35), color="link1", name="bus"))
        for side in (1, -1):
            sc.add(
                RigidBody(
                    R,
                    ZeroVector,
                    Box(0.9, 0.35, 0.015),
                    offset=[side * 0.72, 0, 0],
                    axes=False,
                    color="link2",
                    name=f"panel {'+' if side > 0 else '-'}",
                )
            )
        for i, (axis, b) in enumerate(zip("xyz", bs), 1):
            sc.add(
                RigidBody(
                    R,
                    ZeroVector,
                    Cylinder(0.09, 0.04, axis),
                    offset=b * 0.3,
                    axes=False,
                    color="rotor",
                    name=f"wheel {i}",
                )
            )
        for point in markers:
            sc.add(PointMass(MVMul(R, point), radius=0.02))
        return sc

    return scene


@register("SO(3) × R³")
def spacecraft_reaction_wheels():
    """Rigid spacecraft with three reaction wheels spinning about the body axes b₁, b₂, b₃,
    on SO(3) × R³: attitude R, body angular velocity Ω, wheel angles θᵢ relative to the bus.

    Each wheel is axisymmetric about bᵢ with its centre of mass on the axis; ``Iᵢ`` is its
    axial moment of inertia and its absolute spin rate about bᵢ is ``θ̇ᵢ + bᵢ·Ω``. ``J_b``
    is the inertia of the spacecraft about its centre of mass with the wheels locked,
    minus the wheels' axial inertias: the bus, plus each wheel's mass at its position and
    its transverse inertia. Then
    ``L = ½ Ω·J_bΩ + Σ ½ Iᵢ (θ̇ᵢ + bᵢ·Ω)²`` (no potential).

    The motor of wheel i applies ``+τᵢ bᵢ`` to the wheel and ``−τᵢ bᵢ`` to the bus; the
    wheel's rotation variation is ``δθᵢ bᵢ + η`` and the bus's ``η``, so
    ``δW = Σ τᵢ (δθᵢ + bᵢ·η) − Σ τᵢ bᵢ·η`` (which is ``Σ τᵢ δθᵢ``), inputs
    ``[τ₁, τ₂, τ₃]``. Without input the energy and the inertial angular momentum
    ``R (J_bΩ + Σ Iᵢ (θ̇ᵢ + bᵢ·Ω) bᵢ)`` are conserved.

    Lee et al. §9.10 treat n wheels with general axes on SE(3) × (S¹)ⁿ; with symmetric
    wheels whose axes pass through their centres of mass, no translation and ``sᵢ = bᵢ``,
    their eqs. 9.61–9.62 are the equations of this system with ``αᵢ = Iᵢ``, ``ωᵢ = θ̇ᵢ``
    and ``J′ = J_b + Σ Iᵢ bᵢbᵢᵀ``. The book writes ``J′ = J_b + Σ (mᵢ S(ρᵢ)ᵀS(ρᵢ) + Jᵢ)``,
    using ``Qᵢ Jᵢ Qᵢᵀ = Jᵢ``, which holds for a wheel inertia ``Jᵢ = diag(αᵢ, βᵢ, βᵢ)`` only
    when ``sᵢ = e₁``. For a general axis ``Qᵢ Jᵢ Qᵢᵀ = αᵢ sᵢsᵢᵀ + βᵢ (I − sᵢsᵢᵀ)``; this
    library uses that form, with the transverse part ``βᵢ (I − sᵢsᵢᵀ)`` and the mass term
    in ``J_b``. Here θᵢ ∈ R (the angle is cyclic: L does not depend on it).

    Scene: the bus is a box with two solar panels; each wheel is a cylinder centred at
    ``0.3 bᵢ`` (it follows the parameters), its disc normal along the body axis e₁, e₂ or
    e₃ (the default bᵢ; a Cylinder's axis is a fixed body axis, so a non-default bᵢ moves
    the wheel but does not tilt its disc). The wheels do not turn with θᵢ: a scene cannot
    express ``R expm(θᵢ b̂ᵢ)`` (no sine or cosine), and an axisymmetric wheel looks the same
    when it spins anyway.
    """
    R = SO3("R")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    Jb = _inertia("J_b")
    Is = getScalars("I_1 I_2 I_3", attr=["Constant"])
    bs = [Vector(f"b{i}", attr=["Constant", "UnitNorm"]) for i in (1, 2, 3)]
    thetas = [Scalar(f"\\theta_{i}") for i in (1, 2, 3)]
    taus = [Scalar(f"\\tau_{i}") for i in (1, 2, 3)]
    L = HALF * Dot(Om, Jb * Om)
    dW = None
    for I_i, b, th, tau in zip(Is, bs, thetas, taus):
        spin = th.t_diff() + Dot(b, Om)  # absolute spin rate of the wheel about b
        L = L + HALF * I_i * spin * spin
        term = tau * (th.delta() + Dot(b, eta)) - tau * Dot(b, eta)
        dW = term if dW is None else dW + term
    scene = _spacecraft_scene(R, bs)

    return make_system(
        "spacecraft_reaction_wheels",
        "SO(3) × R³",
        L,
        dW,
        SystemVariables(scalars=thetas, matrices=[R]),
        taus,
        {
            "J_b": J_BUS,
            "I_1": 0.02,
            "I_2": 0.02,
            "I_3": 0.02,
            "b1": np.array([1.0, 0.0, 0.0]),
            "b2": np.array([0.0, 1.0, 0.0]),
            "b3": np.array([0.0, 0.0, 1.0]),
        },
        Reference(
            LEE,
            equation="eqs. 9.61–9.62",
            section="§9.10",
            note="Three axisymmetric wheels on the body axes; the special case of eqs."
            " 9.61–9.62 with J = J_b + Σ Iᵢ bᵢbᵢᵀ.",
        ),
        x0={
            "R": rot([0.1, 0.2, 0.3]),
            "Omega_R": [0.05, -0.1, 0.08],
            "theta_1": 0.0,
            "theta_2": 0.0,
            "theta_3": 0.0,
            "theta_1_dot": 30.0,
            "theta_2_dot": -20.0,
            "theta_3_dot": 10.0,
        },
        symbols=dict(R=R, J_b=Jb, I=tuple(Is), b=tuple(bs), theta=tuple(thetas), tau=tuple(taus)),
        scene=scene,
    )


@register("SO(3) × (S¹)³")
def spacecraft_reaction_wheels_s1():
    """The spacecraft of ``spacecraft_reaction_wheels`` with each wheel angle on S¹, on
    SO(3) × (S¹)³: attitude R, body angular velocity Ω, wheel directions qᵢ.

    ``qᵢ`` is a wheel-fixed direction orthogonal to the wheel axis bᵢ, written in the body
    frame: an ``S1`` point with axis ``bᵢ``, so ``q̇ᵢ = ωᵢ bᵢ × qᵢ`` with ``ωᵢ`` the spin rate
    of the wheel relative to the bus (``θ̇ᵢ`` of the R³ model) and ``δqᵢ = γᵢ bᵢ × qᵢ``.
    ``L = ½ Ω·J_bΩ + Σ ½ Iᵢ (ωᵢ + bᵢ·Ω)²`` (no potential) and
    ``δW = Σ τᵢ (γᵢ + bᵢ·η) − Σ τᵢ bᵢ·η``, inputs ``[τ₁, τ₂, τ₃]``. The equations are those of
    ``spacecraft_reaction_wheels`` with ``θ̇ᵢ → ωᵢ``, ``θ̈ᵢ → ω̇ᵢ`` and ``δθᵢ → γᵢ``
    (``tests/test_s1_systems.py`` checks this numerically, and that both simulate to the
    same trajectory). This is the configuration space of Lee et al. §9.10, SE(3) × (S¹)ⁿ,
    without the translation; ``qᵢ`` is a global coordinate of the wheel angle (no 2π
    wrap-around).

    Scene: the bus and wheels of ``spacecraft_reaction_wheels``, with a rim mark at
    ``0.3 bᵢ + 0.09 qᵢ`` on each wheel, which turns with it.
    """
    from geomech.core.base import S1  # not at import time: see geomech.systems.circles

    R = SO3("R")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    Jb = _inertia("J_b")
    Is = getScalars("I_1 I_2 I_3", attr=["Constant"])
    bs = [Vector(f"b{i}", attr=["Constant", "UnitNorm"]) for i in (1, 2, 3)]
    qs = [S1(f"q_{i}", axis=b) for i, b in zip((1, 2, 3), bs)]
    taus = [Scalar(f"\\tau_{i}") for i in (1, 2, 3)]
    L = HALF * Dot(Om, Jb * Om)
    dW = None
    for I_i, b, q, tau in zip(Is, bs, qs, taus):
        spin = q.get_tangent_vector() + Dot(b, Om)  # absolute spin rate of the wheel about b
        L = L + HALF * I_i * spin * spin
        term = tau * (q.get_variation_vector() + Dot(b, eta)) - tau * Dot(b, eta)
        dW = term if dW is None else dW + term
    marks = [b * 0.3 + q * 0.09 for b, q in zip(bs, qs)]

    return make_system(
        "spacecraft_reaction_wheels_s1",
        "SO(3) × (S¹)³",
        L,
        dW,
        SystemVariables(vectors=qs, matrices=[R]),
        taus,
        {
            "J_b": J_BUS,
            "I_1": 0.02,
            "I_2": 0.02,
            "I_3": 0.02,
            "b1": np.array([1.0, 0.0, 0.0]),
            "b2": np.array([0.0, 1.0, 0.0]),
            "b3": np.array([0.0, 0.0, 1.0]),
        },
        Reference(
            LEE,
            equation="eqs. 9.61–9.62",
            section="§9.10",
            note="spacecraft_reaction_wheels with the wheel angles on S¹ (the configuration"
            " space of §9.10 without translation); same equations.",
        ),
        x0={
            "R": rot([0.1, 0.2, 0.3]),
            "Omega_R": [0.05, -0.1, 0.08],
            "q_1": [0.0, 1.0, 0.0],
            "q_2": [0.0, 0.0, 1.0],
            "q_3": [1.0, 0.0, 0.0],
            "omega_q_1": 30.0,
            "omega_q_2": -20.0,
            "omega_q_3": 10.0,
        },
        symbols=dict(R=R, J_b=Jb, I=tuple(Is), b=tuple(bs), q=tuple(qs), tau=tuple(taus)),
        scene=_spacecraft_scene(R, bs, marks),
    )


@register("SO(3) × SO(3)")
def biped_3d():
    """Swing phase of a three-dimensional compass-gait biped on SO(3) × SO(3).

    The stance leg (attitude R₁) is pinned at the stance foot, at the origin, by a
    spherical joint; the swing leg (attitude R₂) hangs from the hip by a spherical joint.
    Each leg is a rigid body of mass m and length l along its body axis e₃, which points
    from the foot to the hip, with its centre of mass a distance a from the foot and
    inertia J about its centre of mass; a point mass m_H sits at the hip. Hip
    ``x_H = R₁e₃ l``, stance-leg centre of mass ``x₁ = R₁e₃ a``, swing-leg centre of mass
    ``x₂ = x_H − R₂e₃ (l − a)``, swing foot ``x_H − R₂e₃ l``.
    ``L = ½ m_H ẋ_H·ẋ_H + ½ m ẋ₁·ẋ₁ + ½ m ẋ₂·ẋ₂ + ½ Ω₁·JΩ₁ + ½ Ω₂·JΩ₂
    − g e₃·(m_H x_H + m x₁ + m x₂)``.

    Inputs ``[M_H, M_A]``, both inertial-frame moments: the hip torque ``M_H`` acts on the
    swing leg and ``−M_H`` on the stance leg; the ankle torque ``M_A`` acts on the stance
    leg from the ground (zero for a point foot):
    ``δW = η₂·(R₂ᵀM_H) − η₁·(R₁ᵀM_H) + η₁·(R₁ᵀM_A)``.

    Only the continuous single-support phase is modelled: foot impacts, the switching of
    the stance foot, ground contact forces and the ground constraint (the swing foot may
    pass through the ground, the stance foot never lifts off) are not, so this is not a
    walking (hybrid) model. Defaults: m_H = 10 kg, m = 5 kg, l = 1 m, a = 0.5 m, J of a
    slender 1 m rod.

    Scene: the legs are cylinders placed from the parameters (``l`` sets where they sit),
    but a Cylinder's length is a number, so they are drawn with the default ``l = 1 m``;
    a run with another ``l`` places the legs correctly and draws them at the wrong
    length. The ground plane is drawn only.
    """
    R1, R2 = SO3("R1"), SO3("R2")
    Om1, eta1 = R1.get_tangent_vector(), R1.get_variation_vector()
    Om2, eta2 = R2.get_tangent_vector(), R2.get_variation_vector()
    mH, m, g, l, a = getScalars("m_H m g l a", attr=["Constant"])
    J = _inertia("J")
    e3 = Vector("e3", attr=["Constant"])
    MH, MA = Vector("M_H"), Vector("M_A")
    xH = R1 * e3 * l
    x1 = R1 * e3 * a
    x2 = xH - R2 * e3 * (l - a)
    vH, v1, v2 = xH.t_diff(), x1.t_diff(), x2.t_diff()
    L = (
        HALF * mH * Dot(vH, vH)
        + HALF * m * Dot(v1, v1)
        + HALF * m * Dot(v2, v2)
        + HALF * Dot(Om1, J * Om1)
        + HALF * Dot(Om2, J * Om2)
        - mH * g * Dot(e3, xH)
        - m * g * Dot(e3, x1)
        - m * g * Dot(e3, x2)
    )
    dW = (
        Dot(eta2, MVMul(Transpose(R2), MH))
        - Dot(eta1, MVMul(Transpose(R1), MH))
        + Dot(eta1, MVMul(Transpose(R1), MA))
    )
    swing_foot = xH - R2 * e3 * l
    params = {"m_H": 10.0, "m": 5.0, "g": 9.81, "l": 1.0, "a": 0.5, "J": J_LEG, "e3": E3}
    leg = params["l"]  # a Cylinder's length is a number: the drawing uses the default l

    def scene(layout, s):
        from geomech.viz import Box, Cylinder, PointMass, RigidBody, Scene, Trace

        sc = Scene(layout)
        sc.add(
            RigidBody(
                R1,
                ZeroVector,
                Cylinder(0.035, leg, "z"),
                offset=e3 * l * HALF,
                color="link1",
                name="stance leg",
            )
        )
        sc.add(
            RigidBody(
                R2,
                xH,
                Cylinder(0.035, leg, "z"),
                offset=e3 * l * (-0.5),
                color="link2",
                name="swing leg",
            )
        )
        sc.add(
            RigidBody(
                R1,
                ZeroVector,
                Box(0.18, 0.09, 0.03),
                axes=False,
                color="cable",
                name="stance foot",
            )
        )
        sc.add(
            RigidBody(
                R2, swing_foot, Box(0.18, 0.09, 0.03), axes=False, color="cable", name="swing foot"
            )
        )
        sc.add(PointMass(xH, radius=0.09, color="load", name="hip"))
        # the ground plane z = 0, drawn only (no contact is modelled)
        sc.add(
            RigidBody(
                I,
                ZeroVector,
                Box(2.4, 2.4, 0.01),
                offset=[0.0, 0.0, -0.005],
                axes=False,
                color="#c8c8c8",
                name="ground",
            )
        )
        sc.add(Trace(swing_foot, name="swing foot path"))
        return sc

    return make_system(
        "biped_3d",
        "SO(3) × SO(3)",
        L,
        dW,
        SystemVariables(matrices=[R1, R2]),
        [MH, MA],
        params,
        Reference(
            "check_eom",
            note="Single-support phase only; impacts and contact are not modelled.",
        ),
        x0={
            "R1": rot([0.0, 0.2, 0.0]),
            "R2": rot([0.0, -0.3, 0.05]),
            "Omega_R1": [0.1, -0.8, 0.05],
            "Omega_R2": [-0.2, 1.5, 0.1],
        },
        symbols=dict(R1=R1, R2=R2, xH=xH, x1=x1, x2=x2, swing_foot=swing_foot, M_H=MH, M_A=MA),
        scene=scene,
    )


__all__ = ["spacecraft_reaction_wheels", "biped_3d"]
