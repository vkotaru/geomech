"""Systems on SO(3) and on products with SO(3): rigid bodies, tops and rigid pendulums.

Ω is the body-frame angular velocity (``Ṙ = R Ω̂``), η the variation (``δR = R η̂``), and a
moment ``M`` on ``R`` enters as ``η·M``. Gravity acts along −e₃. Each docstring states the
Lagrangian and the virtual work.
"""

from __future__ import annotations

import numpy as np

from geomech.core.base import SO3, Matrix, Vector, ZeroVector, getScalars
from geomech.core.operations import Dot, MVMul, Transpose
from geomech.dynamics.variables import SystemVariables
from geomech.systems._base import (
    E3,
    HALF,
    J_A,
    J_B,
    Reference,
    check_count,
    frozen_array,
    make_system,
    register,
    rot,
)

LEE = "lee2018global"

#: Inertias of two slender links about their centres of mass (kg m²), long axis e₃.
J_LINK1 = frozen_array(np.diag([0.084, 0.085, 0.0014]))
J_LINK2 = frozen_array(np.diag([0.043, 0.044, 0.0011]))


def parallel_axis(J_c, m, r) -> np.ndarray:
    """The inertia about a point at ``-r`` from the centre of mass (Steiner):
    ``J_c + m (|r|² I − r rᵀ)``, with *r* from that point to the centre of mass."""
    r = np.asarray(r, dtype=float)
    return frozen_array(np.asarray(J_c) + m * ((r @ r) * np.eye(3) - np.outer(r, r)))


def _inertia(name):
    return Matrix(name, attr=["Constant", "SymmetricMatrix"])


@register("SO(3)")
def free_rigid_body():
    """Free rigid body on SO(3) driven by a body-frame moment: ``J Ω̇ + Ω × JΩ = M``.

    ``L = ½ Ω·JΩ``, ``δW = η·M``, inputs ``[M]``. Lee et al. §6.6.1, eq. 6.16.
    """
    J = _inertia("J")
    M = Vector("M")
    R = SO3("R")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = Dot(Om, J * Om) * HALF
    dW = Dot(eta, M)

    def scene(layout, s):
        from geomech.viz import Box, RigidBody, Scene

        sc = Scene(layout)
        sc.add(RigidBody(R, ZeroVector, Box(0.6, 0.4, 0.2), name="body"))
        return sc

    return make_system(
        "free_rigid_body",
        "SO(3)",
        L,
        dW,
        SystemVariables(matrices=[R]),
        [M],
        {"J": J_B},
        Reference(
            LEE,
            equation="eq. 6.16",
            section="§6.6.1",
            note="Euler's equations; the inertial angular momentum R J Ω is conserved without"
            " input.",
        ),
        x0=dict(R=rot([0.1, 0.2, 0.3]), Omega_R=[0.3, -1.2, 0.8]),
        symbols=dict(J=J, M=M, R=R),
        scene=scene,
    )


@register("SO(3)")
def heavy_top():
    """Heavy top on SO(3): a rigid body spinning about a fixed point, its centre of mass a
    distance l along the body axis e₃ from the pivot.

    ``L = ½ Ω·JΩ − m g l (R e₃)·e₃`` (J about the pivot), ``δW = η·M``, inputs ``[M]``.
    Lee et al. §6.6.2, eq. 6.20 with the book's pivot offset ``ρ₀ = −l e₃``:
    ``J Ω̇ + Ω × JΩ + m g l e₃ × Rᵀe₃ = M``. The default J is axisymmetric (a Lagrange top)
    and x0 spins fast about the body axis.
    """
    J = _inertia("J")
    m, g, l = getScalars("m g l", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    M = Vector("M")
    R = SO3("R")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = Dot(Om, J * Om) * HALF - m * g * l * Dot(R * e3, e3)
    dW = Dot(eta, M)

    def scene(layout, s):
        from geomech.viz import Cylinder, RigidBody, Rod, Scene, Trace

        sc = Scene(layout)
        sc.add(Rod(ZeroVector, R * e3 * l, radius=0.01))
        sc.add(RigidBody(R, ZeroVector, Cylinder(0.25, 0.06, "z"), offset=e3 * l, name="top"))
        sc.add(Trace(R * e3 * l))
        return sc

    return make_system(
        "heavy_top",
        "SO(3)",
        L,
        dW,
        SystemVariables(matrices=[R]),
        [M],
        # about the pivot: a top with J_c = diag(0.03, 0.03, 0.05) about its centre of mass
        {
            "J": parallel_axis(np.diag([0.03, 0.03, 0.05]), 1.0, [0.0, 0.0, 0.3]),
            "m": 1.0,
            "g": 9.81,
            "l": 0.3,
            "e3": E3,
        },
        Reference(
            LEE,
            equation="eq. 6.20",
            section="§6.6.2",
            note="The three-dimensional pendulum of §6.6.2 with the centre"
            " of mass on the body axis above the pivot.",
        ),
        x0=dict(R=rot([0.3, 0.0, 0.0]), Omega_R=[0.0, 0.0, 20.0]),
        symbols=dict(J=J, m=m, g=g, l=l, e3=e3, M=M, R=R),
        scene=scene,
    )


@register("SO(3)")
def rigid_pendulum():
    """Rigid (three-dimensional) pendulum on SO(3) hinged at its pivot, centre of mass at
    ``R ρ`` from the pivot.

    ``L = ½ Ω·JΩ − m g (Rρ)·e₃`` (J about the pivot), ``δW = η·M``, inputs ``[M]``.
    Lee et al. eq. 6.8; §6.6.2 eq. 6.20 with the book's ``ρ₀ = −ρ``:
    ``J Ω̇ + Ω × JΩ + m g ρ × Rᵀe₃ = M``.
    """
    J = _inertia("J")
    rho = Vector("\\rho", attr=["Constant"])
    m, g = getScalars("m g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    M = Vector("M")
    R = SO3("R")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    L = Dot(Om, J * Om) * HALF - m * g * Dot(R * rho, e3)
    dW = Dot(eta, M)

    def scene(layout, s):
        from geomech.viz import Box, RigidBody, Scene, Trace

        sc = Scene(layout)
        sc.add(RigidBody(R, ZeroVector, Box(0.12, 0.05, 1.0), offset=[0, 0, -0.5]))
        sc.add(Trace(R * rho))
        return sc

    return make_system(
        "rigid_pendulum",
        "SO(3)",
        L,
        dW,
        SystemVariables(matrices=[R]),
        [M],
        {
            "J": parallel_axis(J_LINK1, 1.0, [0.0, 0.0, -0.5]),  # link 1 about its end
            "\\rho": np.array([0.0, 0.0, -0.5]),
            "e3": E3,
            "g": 9.81,
            "m": 1.0,
        },
        Reference(
            LEE,
            equation="eq. 6.20",
            section="§6.6.2; general form eq. 6.8",
            note="The three-dimensional pendulum of §6.6.2.",
        ),
        x0=dict(R=rot([0.7, 0.3, 0.0]), Omega_R=[0.4, -0.6, 4.0]),
        symbols=dict(J=J, rho=rho, m=m, g=g, e3=e3, M=M, R=R),
        scene=scene,
    )


@register("SO(3) × SO(3)")
def double_rigid_pendulum():
    """Double rigid pendulum on SO(3) × SO(3): link 1 hinged at the origin, link 2 hinged
    at the point ``R₁ l₁`` of link 1 by a spherical joint.

    Centres of mass ``x₁ = R₁ρ₁``, ``x₂ = R₁l₁ + R₂ρ₂``;
    ``L = ½ Ω₁·J₁Ω₁ + ½ Ω₂·J₂Ω₂ + ½ m₁ ẋ₁·ẋ₁ + ½ m₂ ẋ₂·ẋ₂ − m₁ g x₁·e₃ − m₂ g x₂·e₃``
    (Jᵢ about the centre of mass), ``δW = η₁·M₁ + η₂·M₂``, inputs ``[M1, M2]``.
    """
    J1, J2 = _inertia("J1"), _inertia("J2")
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

    def scene(layout, s):
        from geomech.viz import Box, RigidBody, Scene, Trace

        sc = Scene(layout)
        sc.add(RigidBody(R1, ZeroVector, Box(0.12, 0.05, 1.0), offset=[0, 0, -0.5]))
        sc.add(RigidBody(R2, R1 * l1, Box(0.12, 0.05, 0.8), offset=[0, 0, -0.4]))
        sc.add(Trace(R1 * l1 + R2 * rho2 * 2))
        return sc

    return make_system(
        "double_rigid_pendulum",
        "SO(3) × SO(3)",
        KE - PE,
        dW,
        SystemVariables(matrices=[R1, R2]),
        [M1, M2],
        {
            "J1": J_LINK1,
            "J2": J_LINK2,
            "\\rho_1": np.array([0.0, 0.0, -0.5]),
            "\\rho_2": np.array([0.0, 0.0, -0.4]),
            "l_1": np.array([0.0, 0.0, -1.0]),
            "m1": 1.0,
            "m2": 0.8,
            "g": 9.81,
            "e3": E3,
        },
        Reference("check_eom", note="Two slender links joined by a spherical joint."),
        x0=dict(
            R1=rot([0.7, 0.3, 0.0]),
            R2=rot([1.1, -0.4, 0.3]),
            Omega_R1=[0.4, -0.6, 4.0],
            Omega_R2=[-0.8, 0.5, -5.0],
        ),
        symbols=dict(
            J1=J1,
            J2=J2,
            rho1=rho1,
            rho2=rho2,
            l1=l1,
            m1=m1,
            m2=m2,
            g=g,
            R1=R1,
            R2=R2,
            M1=M1,
            M2=M2,
        ),
        scene=scene,
    )


#: Cross-section of the drawn chain links (m); also the box used for their default inertia.
LINK_WIDTH = (0.06, 0.03)


def chain_link_defaults(k: int) -> dict:
    """Default parameters of link *k* (1-based) of ``n_link_rigid_chain``: a uniform box
    of length ``L = 0.8^(k−1)`` m, mass ``m = 0.8^(k−1)`` kg and cross-section
    ``LINK_WIDTH``, hanging along body −e₃ with its centre of mass half-way, and the box's
    inertia about its centre of mass."""
    L = m = 0.8 ** (k - 1)
    wx, wy = LINK_WIDTH
    J = np.diag([m * (wy**2 + L**2), m * (wx**2 + L**2), m * (wx**2 + wy**2)]) / 12
    return {
        f"m{k}": m,
        f"J{k}": J,
        f"\\rho_{k}": np.array([0.0, 0.0, -L / 2]),
        f"l_{k}": np.array([0.0, 0.0, -L]),
    }


@register("(SO(3))ⁿ")
def n_link_rigid_chain(*, n_links: int = 3):
    """Chain of n rigid links on SO(3)ⁿ: link 1 hinged at the origin by a spherical joint,
    link k + 1 hinged at the end of link k by a spherical joint.

    Link k has attitude Rₖ, mass mₖ, inertia Jₖ about its centre of mass, its centre of
    mass at ``ρₖ`` from its upper joint and its lower joint at ``lₖ`` from its upper joint
    (both constant body-frame vectors; by default along body −e₃, so the chain hangs at
    ``Rₖ = I``). Joints ``p₁ = 0``, ``pₖ₊₁ = pₖ + Rₖlₖ``; centres of mass
    ``xₖ = pₖ + Rₖρₖ``. ``L = Σ ½ Ωₖ·JₖΩₖ + Σ ½ mₖ ẋₖ·ẋₖ − Σ mₖ g xₖ·e₃``; ``lₙ`` of the last
    link does not enter L and is not a parameter.

    Inputs ``[τ₁, …, τₙ]``: inertial-frame torques at the joints, τₖ acting on link k and
    −τₖ on link k − 1 (τ₁ between the ground and link 1):
    ``δW = Σₖ ηₖ·(RₖᵀTₖ)`` with ``Tₖ = τₖ − τₖ₊₁`` (``τₙ₊₁ = 0``), written term by term as
    ``Σₖ (ηₖ·(Rₖᵀτₖ) − ηₖ₋₁·(Rₖ₋₁ᵀτₖ))``.

    With ``n_links=2`` the Lagrangian is that of ``double_rigid_pendulum`` term for term
    (same symbols ``J1, J2, \\rho_1, \\rho_2, l_1, m1, m2, g, e3``); the virtual work
    differs only in the inputs (body moments ``M₁, M₂`` there, joint torques here, related
    by ``M₁ = R₁ᵀ(τ₁ − τ₂)``, ``M₂ = R₂ᵀτ₂``). The default links are uniform boxes
    (``chain_link_defaults``), each 0.8 times the one above. The scene draws each link as
    a box from its upper joint (placed with ``ρₖ``, drawn with the default length) and a
    sphere at each joint. ``n_links`` is an int from 1 to 20 (``MAX_SIZE``); any other
    value raises ``ValueError``.
    """
    n = check_count(n_links, "n_link_rigid_chain: n_links")
    ks = range(1, n + 1)
    Rs = [SO3(f"R{k}") for k in ks]
    Js = [_inertia(f"J{k}") for k in ks]
    rhos = [Vector(f"\\rho_{k}", attr=["Constant"]) for k in ks]
    ls = [Vector(f"l_{k}", attr=["Constant"]) for k in range(1, n)]
    ms = getScalars(" ".join(f"m{k}" for k in ks), attr=["Constant"])
    (g,) = getScalars("g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    taus = [Vector(f"\\tau_{k}") for k in ks]
    joints, xs, p = [], [], None
    for k in range(n):
        joints.append(ZeroVector if p is None else p)
        x = Rs[k] * rhos[k] if p is None else p + Rs[k] * rhos[k]
        xs.append(x)
        if k < n - 1:
            p = Rs[k] * ls[k] if p is None else p + Rs[k] * ls[k]
    terms = [
        Dot(R.get_tangent_vector(), J * R.get_tangent_vector()) * HALF for R, J in zip(Rs, Js)
    ]
    terms += [Dot(x.t_diff(), x.t_diff()) * m * HALF for x, m in zip(xs, ms)]
    KE = terms[0]
    for t in terms[1:]:
        KE = KE + t
    PE = ms[0] * g * Dot(xs[0], e3)
    for m, x in zip(ms[1:], xs[1:]):
        PE = PE + m * g * Dot(x, e3)
    dW = Dot(Rs[0].get_variation_vector(), MVMul(Transpose(Rs[0]), taus[0]))
    for k in range(1, n):
        dW = (
            dW
            + Dot(Rs[k].get_variation_vector(), MVMul(Transpose(Rs[k]), taus[k]))
            - Dot(Rs[k - 1].get_variation_vector(), MVMul(Transpose(Rs[k - 1]), taus[k]))
        )
    params: dict = {"g": 9.81, "e3": E3}
    for k in ks:
        params.update(chain_link_defaults(k))
    params.pop(f"l_{n}")
    lengths = [-chain_link_defaults(k)[f"l_{k}"][2] for k in ks]
    tip = joints[-1] + Rs[-1] * rhos[-1] * 2

    def scene(layout, s):
        from geomech.viz import Box, PointMass, RigidBody, Scene, Trace

        sc = Scene(layout)
        for k, (R, joint, rho, L_k) in enumerate(zip(Rs, joints, rhos, lengths), 1):
            sc.add(RigidBody(R, joint, Box(*LINK_WIDTH, L_k), offset=rho, name=f"link {k}"))
            sc.add(PointMass(joint, radius=0.04, color="cable", name=f"joint {k}"))
        sc.add(Trace(tip, name="tip"))
        return sc

    x0: dict = {}
    for k in ks:
        sign = 1 if k % 2 else -1
        x0[f"R{k}"] = rot([0.4 + 0.15 * k, 0.1 * k * sign, 0.05 * k])
        x0[f"Omega_R{k}"] = [0.3, -0.2 * k, 2.0 * sign]
    return make_system(
        "n_link_rigid_chain",
        "(SO(3))ⁿ",
        KE - PE,
        dW,
        SystemVariables(matrices=Rs),
        taus,
        params,
        Reference("check_eom", note="Chain of rigid links with spherical joints."),
        x0=x0,
        symbols=dict(R=tuple(Rs), joints=tuple(joints), x=tuple(xs), tau=tuple(taus), tip=tip),
        scene=scene,
        options={"n_links": n},
    )


@register("R³ × SO(3) × SO(3)")
def two_rigid_bodies_spherical_joint():
    """Two free-flying rigid bodies joined by a spherical joint, on R³ × SO(3) × SO(3).

    Body 1 has its centre of mass at x; the joint is at ``x + R₁ρ₁`` and body 2 has its
    centre of mass at ``x₂ = x + R₁ρ₁ + R₂ρ₂``.
    ``L = ½ m₁ ẋ·ẋ + ½ m₂ ẋ₂·ẋ₂ + ½ Ω₁·J₁Ω₁ + ½ Ω₂·J₂Ω₂ − m₁ g e₃·x − m₂ g e₃·x₂``,
    ``δW = δx·F + η₁·M₁ + η₂·M₂`` (F at the centre of mass of body 1), inputs
    ``[F, M1, M2]``.
    """
    m1, m2, g = getScalars("m1 m2 g", attr=["Constant"])
    J1, J2 = _inertia("J1"), _inertia("J2")
    rho1 = Vector("\\rho_1", attr=["Constant"])
    rho2 = Vector("\\rho_2", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    x, R1, R2 = Vector("x"), SO3("R1"), SO3("R2")
    F, M1, M2 = Vector("F"), Vector("M1"), Vector("M2")
    Om1, Om2 = R1.get_tangent_vector(), R2.get_tangent_vector()
    eta1, eta2 = R1.get_variation_vector(), R2.get_variation_vector()
    x2 = x + R1 * rho1 + R2 * rho2
    v1, v2 = x.t_diff(), x2.t_diff()
    L = (
        HALF * m1 * Dot(v1, v1)
        + HALF * m2 * Dot(v2, v2)
        + HALF * Dot(Om1, J1 * Om1)
        + HALF * Dot(Om2, J2 * Om2)
        - m1 * g * Dot(e3, x)
        - m2 * g * Dot(e3, x2)
    )
    dW = Dot(x.delta(), F) + Dot(eta1, M1) + Dot(eta2, M2)

    def scene(layout, s):
        from geomech.viz import Box, RigidBody, Scene, Trace

        sc = Scene(layout)
        sc.add(RigidBody(R1, x, Box(0.6, 0.15, 0.1), name="body 1"))
        sc.add(
            RigidBody(R2, x + R1 * rho1, Box(0.5, 0.12, 0.08), offset=[0.25, 0, 0], name="body 2")
        )
        sc.add(Trace(x))
        return sc

    return make_system(
        "two_rigid_bodies_spherical_joint",
        "R³ × SO(3) × SO(3)",
        L,
        dW,
        SystemVariables(vectors=[x], matrices=[R1, R2]),
        [F, M1, M2],
        {
            "m1": 1.0,
            "m2": 0.6,
            "g": 9.81,
            "J1": np.diag([0.004, 0.03, 0.032]),
            "J2": np.diag([0.002, 0.013, 0.014]),
            "\\rho_1": np.array([0.3, 0.0, 0.0]),
            "\\rho_2": np.array([0.25, 0.0, 0.0]),
            "e3": E3,
        },
        Reference("check_eom", note="Free-flying two-body chain."),
        x0=dict(
            x=[0.0, 0.0, 1.0],
            R1=rot([0.1, 0.2, 0.0]),
            R2=rot([0.0, -0.4, 0.9]),
            x_dot=[0.2, 0.0, 2.0],
            Omega_R1=[0.5, -1.0, 2.0],
            Omega_R2=[-1.0, 0.8, -1.5],
        ),
        symbols=dict(x=x, x2=x2, R1=R1, R2=R2, F=F, M1=M1, M2=M2),
        scene=scene,
    )


@register("R³ × SO(3)")
def rigid_body():
    """Rigid body in R³ × SO(3) (six degrees of freedom) under gravity, driven by a force
    at the centre of mass and a body-frame moment.

    ``L = ½ m ẋ·ẋ + ½ Ω·JΩ − m g e₃·x``, ``δW = δx·F + η·M``, inputs ``[F, M]``.
    Lee et al. §7.6.1, eqs. 7.39–7.40 (``J Ω̇ + Ω × JΩ = M``, ``m ẍ + m g e₃ = F``).
    """
    m, g = getScalars("m g", attr=["Constant"])
    J = _inertia("J")
    e3 = Vector("e3", attr=["Constant"])
    x, R, F, M = Vector("x"), SO3("R"), Vector("F"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    v = x.t_diff()
    L = HALF * m * Dot(v, v) + HALF * Dot(Om, J * Om) - m * g * Dot(e3, x)
    dW = Dot(x.delta(), F) + Dot(eta, M)

    def scene(layout, s):
        from geomech.viz import Box, RigidBody, Scene, Trace

        sc = Scene(layout)
        sc.add(RigidBody(R, x, Box(0.6, 0.4, 0.2), name="body"))
        sc.add(Trace(x))
        return sc

    return make_system(
        "rigid_body",
        "R³ × SO(3)",
        L,
        dW,
        SystemVariables(vectors=[x], matrices=[R]),
        [F, M],
        {"m": 2.0, "g": 9.81, "J": J_A * 10, "e3": E3},
        Reference(
            LEE,
            equation="eqs. 7.39–7.40",
            section="§7.6.1",
            note="SE(3) modelled as the product R³ × SO(3), with the inertial velocity ẋ.",
        ),
        x0=dict(
            x=[0.0, 0.0, 1.0],
            R=rot([0.2, -0.1, 0.4]),
            x_dot=[0.5, 0.0, 3.0],
            Omega_R=[0.3, -1.2, 0.8],
        ),
        symbols=dict(m=m, g=g, J=J, e3=e3, x=x, R=R, F=F, M=M),
        scene=scene,
    )


__all__ = [
    "free_rigid_body",
    "heavy_top",
    "rigid_pendulum",
    "double_rigid_pendulum",
    "n_link_rigid_chain",
    "two_rigid_bodies_spherical_joint",
    "rigid_body",
]
