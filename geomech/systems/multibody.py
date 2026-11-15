"""Multibody systems on SO(3) and its products: a spacecraft with reaction wheels (wheel
angles on R³ or on (S¹)³), the same spacecraft with two solar-array drives, and the swing
phase of a three-dimensional compass-gait biped.

Ω is the body-frame angular velocity (``Ṙ = R Ω̂``), η the variation (``δR = R η̂``); a
body-frame moment ``M`` on ``R`` enters as ``η·M`` and an inertial-frame moment ``M`` as
``η·(RᵀM)``. Gravity acts along −e₃. Each docstring states the Lagrangian and the virtual
work.
"""

from __future__ import annotations

import numpy as np

from geomech.core.base import SO3, I, Matrix, Scalar, Vector, ZeroVector, getScalars
from geomech.core.operations import Cross, Dot, MVMul, Transpose
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


#: Solar-array panel: a thin plate of 0.9 m along its drive axis and 0.35 m across, its
#: centre 0.72 m from the bus centre (the scene's panel geometry).
PANEL_SIZE = (0.9, 0.35, 0.015)
PANEL_OFFSET = 0.72


def _spacecraft_scene(R, bs, markers=(), panels=None):
    """The bus (a box with two solar panels) and a cylinder per wheel centred at ``0.3 bᵢ``
    with its disc normal along the body axis e₁, e₂ or e₃; *markers* are body-frame points
    drawn as small spheres (the rim mark of an S¹ wheel).

    *panels* (``None``: two panels fixed to the bus at ``±0.72 e₁``) is a list of
    ``(P, offset)``: the panel's attitude ``P`` relative to the bus (a matrix expression)
    and the centre of its plate in the panel frame (a constant vector expression); the
    panel is drawn with attitude ``R P``."""

    def scene(layout, s):
        from geomech.viz import Box, Cylinder, PointMass, RigidBody, Scene

        sc = Scene(layout)
        sc.add(RigidBody(R, ZeroVector, Box(0.5, 0.4, 0.35), color="link1", name="bus"))
        if panels is None:
            fixed = [(R, [side * PANEL_OFFSET, 0, 0], side) for side in (1, -1)]
        else:
            fixed = [(R * P, off, k) for k, (P, off) in enumerate(panels, 1)]
        for att, off, tag in fixed:
            name = f"panel {'+' if tag == 1 else '-'}" if panels is None else f"panel {tag}"
            sc.add(
                RigidBody(
                    att,
                    ZeroVector,
                    Box(*PANEL_SIZE),
                    offset=off,
                    axes=False,
                    color="link2",
                    name=name,
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


@register("SO(3) × (S¹)³ × (S¹)²")
def spacecraft_solar_arrays():
    """The spacecraft of ``spacecraft_reaction_wheels_s1`` with two solar-array panels, each
    a rigid plate turned relative to the bus by a solar-array drive, on
    SO(3) × (S¹)³ × (S¹)²: attitude R, body angular velocity Ω, wheel directions qᵢ and
    panel normals nₖ.

    Panel k turns about the bus-fixed unit axis ``aₖ`` (its drive axis, along the panel's
    length; ``a₁ = b₁``, ``a₂ = −b₁``, the panel geometry of the scene). ``nₖ`` is the
    panel's normal written in the body frame: an ``S1`` point with axis ``aₖ``, so
    ``ṅₖ = ϖₖ aₖ × nₖ`` with ``ϖₖ`` (printed ``ω_{n_k}``) the drive rate relative to the
    bus, and ``δnₖ = γₖ aₖ × nₖ``. The panel's centre of mass lies on its drive axis at
    ``ρₖ = dₖ aₖ``, so it does not move in the bus frame when the panel turns; its mass is
    ``m_pk`` and its principal moments of inertia about its centre of mass are ``Ia_k``
    about ``aₖ``, ``In_k`` about ``nₖ`` and ``Iw_k`` about the in-plane width direction
    ``aₖ × nₖ``.

    Kinetic energy of each body, from first principles. Model assumption: free floating,
    no external force, total linear momentum zero. The composite centre of mass then stays
    at the inertial origin, and because every body's centre of mass is fixed in the bus
    frame (the wheels' on their axes, the panels' on their drive axes) it is a fixed point
    of the bus frame, taken as the body-frame origin. ``ρₖ`` is measured from it and
    ``J_b`` (the bus with the wheels locked, minus the wheels' axial inertias, as in
    ``spacecraft_reaction_wheels``) is about it; for the default symmetric panels it is
    the bus centre. Then the velocity of panel k's centre of mass is ``R (Ω × ρₖ)``
    (translation induced by the bus rotation), its angular velocity in the body frame is
    ``Ω + ϖₖ aₖ``, and its inertia in the body frame is
    ``Ia_k aₖaₖᵀ + In_k nₖnₖᵀ + Iw_k (aₖ × nₖ)(aₖ × nₖ)ᵀ``; since ``nₖ·aₖ = 0``,

    ``L = ½ Ω·J_bΩ + Σᵢ ½ Iᵢ (ωᵢ + bᵢ·Ω)²
    + Σₖ [½ m_pk dₖ² |Ω × aₖ|² + ½ Ia_k (ϖₖ + aₖ·Ω)² + ½ In_k (nₖ·Ω)²
    + ½ Iw_k ((aₖ × nₖ)·Ω)²]`` (no potential).

    Each motor (wheel torque τᵢ, drive torque μₖ) acts equal and opposite on the bus:
    ``δW = Σᵢ τᵢ (γᵢ + bᵢ·η) − Σᵢ τᵢ bᵢ·η + Σₖ μₖ (γₖ + aₖ·η) − Σₖ μₖ aₖ·η``, inputs
    ``[τ₁, τ₂, τ₃, μ₁, μ₂]``. Without input (or with internal torques only) the inertial
    angular momentum
    ``R (J_bΩ + Σᵢ Iᵢ(ωᵢ + bᵢ·Ω) bᵢ + Σₖ [m_pk dₖ² (Ω − (aₖ·Ω) aₖ) + Jₖ(Ω + ϖₖ aₖ)])``
    is conserved, with ``Jₖ`` the panel inertia above. With ``m_pk = Ia_k = In_k = Iw_k =
    0`` and ``μₖ = 0`` the bus and wheel equations are those of
    ``spacecraft_reaction_wheels_s1``. Verified by ``check_eom`` and by an independent
    model of the bodies (``tests/test_systems_physics.py``); no closed form is compared.

    Defaults: the bus and wheels of ``spacecraft_reaction_wheels``; each panel a uniform
    thin plate of 1 kg, 0.9 m along ``aₖ`` and 0.35 m across, centred at ``dₖ = 0.72`` m,
    so ``Ia = m w²/12 = 0.0102``, ``Iw = m ℓ²/12 = 0.0675`` and ``In = Ia + Iw = 0.0777``
    (kg m²).

    Scene: the bus and wheels of ``spacecraft_reaction_wheels_s1`` (with their rim marks);
    panel k is drawn as the plate with attitude ``R Pₖ``,
    ``Pₖ = aₖb₁ᵀ + (nₖ × aₖ)b₂ᵀ + nₖb₃ᵀ`` (a rotation when b₁, b₂, b₃ are an orthonormal
    right-handed frame, as by default), centred at ``R dₖ aₖ``; it turns with ``nₖ``.
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
    axes = [Vector(f"a{k}", attr=["Constant", "UnitNorm"]) for k in (1, 2)]
    ns = [S1(f"n_{k}", axis=a) for k, a in zip((1, 2), axes)]
    mus = [Scalar(f"\\mu_{k}") for k in (1, 2)]
    masses, dists, Ias, Ins, Iws = [], [], [], [], []
    for k, a, n, mu in zip((1, 2), axes, ns, mus):
        m, d, Ia, In, Iw = getScalars(f"m_p{k} d_{k} Ia_{k} In_{k} Iw_{k}", attr=["Constant"])
        v = Cross(Om, a)  # the velocity of the centre of mass in the body frame, over d
        spin = n.get_tangent_vector() + Dot(a, Om)  # the panel's absolute rate about a
        w = Cross(a, n)  # the in-plane width direction
        L = (
            L
            + HALF * m * d * d * Dot(v, v)
            + HALF * Ia * spin * spin
            + HALF * In * Dot(n, Om) * Dot(n, Om)
            + HALF * Iw * Dot(w, Om) * Dot(w, Om)
        )
        dW = dW + mu * (n.get_variation_vector() + Dot(a, eta)) - mu * Dot(a, eta)
        for lst, x in zip((masses, dists, Ias, Ins, Iws), (m, d, Ia, In, Iw)):
            lst.append(x)
    marks = [b * 0.3 + q * 0.09 for b, q in zip(bs, qs)]
    b1, b2, b3 = bs
    panels = [
        (a * Transpose(b1) + Cross(n, a) * Transpose(b2) + n * Transpose(b3), b1 * d)
        for a, n, d in zip(axes, ns, dists)
    ]
    m_p = 1.0
    length, width, _ = PANEL_SIZE
    J_a, J_w = m_p * width**2 / 12, m_p * length**2 / 12
    params = {
        "J_b": J_BUS,
        "I_1": 0.02,
        "I_2": 0.02,
        "I_3": 0.02,
        "b1": np.array([1.0, 0.0, 0.0]),
        "b2": np.array([0.0, 1.0, 0.0]),
        "b3": np.array([0.0, 0.0, 1.0]),
        "a1": np.array([1.0, 0.0, 0.0]),
        "a2": np.array([-1.0, 0.0, 0.0]),
    }
    for k in (1, 2):
        params.update(
            {
                f"m_p{k}": m_p,
                f"d_{k}": PANEL_OFFSET,
                f"Ia_{k}": J_a,
                f"In_{k}": J_a + J_w,
                f"Iw_{k}": J_w,
            }
        )

    return make_system(
        "spacecraft_solar_arrays",
        "SO(3) × (S¹)³ × (S¹)²",
        L,
        dW,
        SystemVariables(vectors=qs + ns, matrices=[R]),
        taus + mus,
        params,
        Reference(
            "check_eom",
            note="spacecraft_reaction_wheels_s1 with two rigid solar-array panels on drives"
            " about ±b₁; free floating with zero linear momentum, so the composite centre of"
            " mass is a fixed point of the bus frame and the panel translation is induced by"
            " the bus rotation.",
        ),
        x0={
            "R": rot([0.1, 0.2, 0.3]),
            "Omega_R": [0.05, -0.1, 0.08],
            "q_1": [0.0, 1.0, 0.0],
            "q_2": [0.0, 0.0, 1.0],
            "q_3": [1.0, 0.0, 0.0],
            "n_1": [0.0, 0.0, 1.0],
            "n_2": [0.0, 0.0, 1.0],
            "omega_q_1": 30.0,
            "omega_q_2": -20.0,
            "omega_q_3": 10.0,
            "omega_n_1": 0.5,
            "omega_n_2": -0.4,
        },
        symbols=dict(
            R=R,
            J_b=Jb,
            I=tuple(Is),
            b=tuple(bs),
            q=tuple(qs),
            tau=tuple(taus),
            a=tuple(axes),
            n=tuple(ns),
            mu=tuple(mus),
            m_p=tuple(masses),
            d=tuple(dists),
            Ia=tuple(Ias),
            In=tuple(Ins),
            Iw=tuple(Iws),
        ),
        scene=_spacecraft_scene(R, bs, marks, panels),
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


__all__ = [
    "spacecraft_reaction_wheels",
    "spacecraft_reaction_wheels_s1",
    "spacecraft_solar_arrays",
    "biped_3d",
]
