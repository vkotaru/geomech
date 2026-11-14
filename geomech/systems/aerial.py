"""Quadrotors and quadrotors carrying loads, on products of R³, SO(3), S² and R.

A quadrotor is a rigid body (attitude R, body-frame Ω) whose four rotors give a thrust
``f`` along the body axis ``R e₃`` and a body-frame moment ``M``: ``δW = f (R e₃)·δx_Q + η·M``
with ``δx_Q`` the variation of the point where the thrust acts. Gravity acts along −e₃
(``e₃`` points up here). A cable direction ``q`` points from the quadrotor to what hangs
from it, so a hanging cable has ``q = −e₃``. Each docstring states the Lagrangian.
"""

from __future__ import annotations

import numpy as np

from geomech.core.base import S2, SO3, Matrix, Scalar, Vector, getScalars
from geomech.core.operations import Cross, Dot
from geomech.dynamics.variables import SystemVariables
from geomech.systems._base import (
    E3,
    HALF,
    Reference,
    check_count,
    frozen_array,
    make_system,
    register,
    rot,
    tangent,
    unit,
)

J_QUAD = frozen_array([[0.02, 0.0005, 0.0], [0.0005, 0.025, 0.001], [0.0, 0.001, 0.04]])
ARM = 0.25


def _inertia(name):
    return Matrix(name, attr=["Constant", "SymmetricMatrix"])


def add_quadrotor(scene, R, x, label="", first=True):
    """A quadrotor drawn as a box with two arms and four rotor cylinders, all with the
    attitude *R* about the point *x* (the rotor on +e₁ is marked)."""
    from geomech.viz import Box, Cylinder, RigidBody

    tag = f" {label}" if label else ""
    scene.add(RigidBody(R, x, Box(0.16, 0.12, 0.06), color="link1", name=f"body{tag}"))
    scene.add(
        RigidBody(R, x, Box(2 * ARM, 0.025, 0.02), axes=False, color="cable", name=f"arm 1{tag}")
    )
    scene.add(
        RigidBody(R, x, Box(0.025, 2 * ARM, 0.02), axes=False, color="cable", name=f"arm 2{tag}")
    )
    offsets = [(ARM, 0, 0.025), (0, ARM, 0.025), (-ARM, 0, 0.025), (0, -ARM, 0.025)]
    for k, off in enumerate(offsets):
        color = "link2" if k == 0 else "rotor"
        scene.add(
            RigidBody(
                R,
                x,
                Cylinder(0.09, 0.012, "z"),
                offset=off,
                axes=False,
                color=color,
                name=f"rotor {k + 1}{tag}",
            )
        )


def _hover_x0(prefix_R="R"):
    return {prefix_R: rot([0.2, -0.1, 0.4]), f"Omega_{prefix_R}": [0.3, -0.2, 0.5]}


# ---------------------------------------------------------------------------
# A single quadrotor
# ---------------------------------------------------------------------------


@register("R³ × SO(3)")
def quadrotor():
    """Quadrotor on R³ × SO(3), thrust through the centre of mass.

    ``L = ½ m ẋ·ẋ + ½ Ω·JΩ − m g e₃·x``, ``δW = f (R e₃)·δx + M·η``, inputs ``[f, M]``
    (f a scalar). The model of Lee, Leok and McClamroch (2010), with e₃ pointing up here.
    """
    m, g = getScalars("m g", attr=["Constant"])
    J = _inertia("J")
    e3 = Vector("e3", attr=["Constant"])
    x, R = Vector("x"), SO3("R")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    f, M = Scalar("f"), Vector("M")
    v = x.t_diff()
    L = HALF * m * Dot(v, v) + HALF * Dot(Om, J * Om) - m * g * Dot(e3, x)
    dW = Dot(x.delta(), (R * e3) * f) + Dot(eta, M)

    def scene(layout, s):
        from geomech.viz import Scene, Trace

        sc = Scene(layout)
        add_quadrotor(sc, R, x)
        sc.add(Trace(x))
        return sc

    return make_system(
        "quadrotor",
        "R³ × SO(3)",
        L,
        dW,
        SystemVariables(vectors=[x], matrices=[R]),
        [f, M],
        {"J": J_QUAD, "e3": E3, "g": 9.81, "m": 2.0},
        Reference(
            "lee2010geometric",
            note="Rigid body with thrust along R e₃; the paper's e₃ points down.",
        ),
        x0=dict(x=[0.0, 0.5, 1.0], x_dot=[0.3, -0.2, 0.5], **_hover_x0()),
        symbols=dict(m=m, g=g, J=J, e3=e3, x=x, R=R, f=f, M=M),
        scene=scene,
    )


@register("R³ × SO(3)")
def quadrotor_offset_com():
    """Quadrotor whose centre of mass is offset from the rotor plane centre, on R³ × SO(3).

    x is the centre of the rotor plane, where the thrust acts; the centre of mass is at
    ``x_c = x + Rρ`` (ρ constant in the body frame). ``L = ½ m ẋ_c·ẋ_c + ½ Ω·JΩ − m g e₃·x_c``
    (J about the centre of mass), ``δW = f (R e₃)·δx + M·η``, inputs ``[f, M]``. The offset
    couples translation and rotation through the inertia and through gravity.
    """
    m, g = getScalars("m g", attr=["Constant"])
    J = _inertia("J")
    e3 = Vector("e3", attr=["Constant"])
    rho = Vector("\\rho", attr=["Constant"])
    x, R = Vector("x"), SO3("R")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    f, M = Scalar("f"), Vector("M")
    xc = x + R * rho
    vc = xc.t_diff()
    L = HALF * m * Dot(vc, vc) + HALF * Dot(Om, J * Om) - m * g * Dot(e3, xc)
    dW = Dot(x.delta(), (R * e3) * f) + Dot(eta, M)

    def scene(layout, s):
        from geomech.viz import PointMass, Scene, Trace

        sc = Scene(layout)
        add_quadrotor(sc, R, x)
        sc.add(PointMass(xc, radius=0.03, color="load", name="centre of mass"))
        sc.add(Trace(x))
        return sc

    return make_system(
        "quadrotor_offset_com",
        "R³ × SO(3)",
        L,
        dW,
        SystemVariables(vectors=[x], matrices=[R]),
        [f, M],
        {"J": J_QUAD, "e3": E3, "g": 9.81, "m": 2.0, "\\rho": np.array([0.02, -0.01, -0.08])},
        Reference("check_eom", note="Centre of mass offset from the thrust point."),
        x0=dict(x=[0.0, 0.5, 1.0], x_dot=[0.3, -0.2, 0.5], **_hover_x0()),
        symbols=dict(m=m, g=g, J=J, e3=e3, rho=rho, x=x, x_c=xc, R=R, f=f, M=M),
        scene=scene,
    )


# ---------------------------------------------------------------------------
# A quadrotor with a load
# ---------------------------------------------------------------------------


def _load_scene(R, xQ, xL, rods=()):
    def build(layout, s):
        from geomech.viz import PointMass, Rod, Scene, Trace

        sc = Scene(layout)
        add_quadrotor(sc, R, xQ)
        pts = [xQ, *rods, xL]
        for a, b in zip(pts, pts[1:]):
            sc.add(Rod(a, b, radius=0.006))
        for p in rods:
            sc.add(PointMass(p, radius=0.02, color="cable"))
        sc.add(PointMass(xL, radius=0.07, color="load", name="load"))
        sc.add(Trace(xL))
        return sc

    return build


@register("R³ × SO(3) × S²")
def quadrotor_suspended_load():
    """Quadrotor with a point-mass load on a rigid massless cable, on R³ × SO(3) × S².

    Load at x_L, cable direction q from the quadrotor to the load, cable length l, so
    ``x_Q = x_L − l q``. ``L = ½ m_Q ẋ_Q·ẋ_Q + ½ m_L ẋ_L·ẋ_L + ½ Ω·JΩ − m_Q g e₃·x_Q
    − m_L g e₃·x_L``, ``δW = f (R e₃)·(δx_L − l ξ × q) + M·η``, inputs ``[f, M]``. The
    model of Sreenath, Lee and Kumar (2013).
    """
    mQ, mL, g, l = getScalars("m_Q m_L g l", attr=["Constant"])
    J = _inertia("J")
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
    q0 = unit([np.sin(0.5), 0.0, -np.cos(0.5)])
    return make_system(
        "quadrotor_suspended_load",
        "R³ × SO(3) × S²",
        L,
        dW,
        SystemVariables(vectors=[xL, q], matrices=[R]),
        [f, M],
        {"J": J_QUAD, "e3": E3, "g": 9.81, "l": 1.0, "m_L": 0.3, "m_Q": 1.5},
        Reference("sreenath2013geometric", note="Point-mass load on a taut massless cable."),
        x0=dict(
            x_L=[0.0, 0.2, 0.5],
            q=q0,
            x_L_dot=[0.2, -0.1, 0.3],
            omega_q=tangent(q0, [0.0, 1.0, 0.0]),
            **_hover_x0(),
        ),
        symbols=dict(mQ=mQ, mL=mL, g=g, l=l, J=J, e3=e3, xL=xL, q=q, R=R, f=f, M=M, xQ=xQ),
        scene=_load_scene(R, xQ, xL),
    )


@register("R³ × SO(3) × R × S²")
def quadrotor_elastic_load():
    """Quadrotor with a point-mass load on an elastic cable, on R³ × SO(3) × R × S².

    As ``quadrotor_suspended_load`` with the cable length l a configuration variable and
    a linear spring of stiffness k and rest length l₀: ``x_Q = x_L − l q``,
    ``L = ½ m_Q ẋ_Q·ẋ_Q + ½ m_L ẋ_L·ẋ_L + ½ Ω·JΩ − m_Q g e₃·x_Q − m_L g e₃·x_L
    − ½ k (l − l₀)²``, ``δW = f (R e₃)·(δx_L − δl q − l ξ × q) + M·η``, inputs ``[f, M]``.
    The cable stays taut only while the spring is stretched (not enforced).
    """
    mQ, mL, g, k, l0 = getScalars("m_Q m_L g k l0", attr=["Constant"])
    J = _inertia("J")
    e3 = Vector("e3", attr=["Constant"])
    l = Scalar("l")
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
        - HALF * k * (l - l0) * (l - l0)
    )
    dW = Dot(xL.delta() - q * l.delta() - Cross(xi, q) * l, (R * e3) * f) + Dot(eta, M)
    q0 = unit([np.sin(0.5), 0.0, -np.cos(0.5)])
    return make_system(
        "quadrotor_elastic_load",
        "R³ × SO(3) × R × S²",
        L,
        dW,
        SystemVariables(scalars=[l], vectors=[xL, q], matrices=[R]),
        [f, M],
        {"J": J_QUAD, "e3": E3, "g": 9.81, "k": 200.0, "l0": 1.0, "m_L": 0.3, "m_Q": 1.5},
        Reference("check_eom", note="Point-mass load on an elastic cable."),
        x0=dict(
            l=1.05,
            x_L=[0.0, 0.2, 0.5],
            q=q0,
            l_dot=0.1,
            x_L_dot=[0.2, -0.1, 0.3],
            omega_q=tangent(q0, [0.0, 1.0, 0.0]),
            **_hover_x0(),
        ),
        symbols=dict(xL=xL, q=q, R=R, l=l, f=f, M=M, xQ=xQ),
        scene=_load_scene(R, xQ, xL),
    )


@register("R³ × SO(3) × (S²)ⁿ")
def quadrotor_flexible_cable(*, n_links: int = 3):
    """Quadrotor carrying a load on a flexible cable modelled as n rigid links, on
    R³ × SO(3) × (S²)ⁿ.

    Quadrotor at x_Q; link i has direction qᵢ (from the quadrotor end), length lᵢ and a
    point mass mᵢ at its far end, ``xᵢ = x_Q + Σ_{j≤i} lⱼ qⱼ``; mₙ is the load.
    ``L = ½ m_Q ẋ_Q·ẋ_Q + ½ Ω·JΩ + Σ ½ mᵢ ẋᵢ·ẋᵢ − m_Q g e₃·x_Q − Σ mᵢ g e₃·xᵢ``,
    ``δW = f (R e₃)·δx_Q + M·η``, inputs ``[f, M]``.
    ``n_links`` is an int from 1 to 20 (``MAX_SIZE``); any other value raises
    ``ValueError``.
    """
    n = check_count(n_links, "quadrotor_flexible_cable: n_links")
    mQ, g = getScalars("m_Q g", attr=["Constant"])
    masses = getScalars(" ".join(f"m{i}" for i in range(1, n + 1)), attr=["Constant"])
    lengths = getScalars(" ".join(f"l{i}" for i in range(1, n + 1)), attr=["Constant"])
    J = _inertia("J")
    e3 = Vector("e3", attr=["Constant"])
    xQ, R = Vector("x_Q"), SO3("R")
    qs = [S2(f"q{i}") for i in range(1, n + 1)]
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    f, M = Scalar("f"), Vector("M")
    points, xp = [], xQ
    for qi, li in zip(qs, lengths):
        xp = xp + qi * li
        points.append(xp)
    vQ = xQ.t_diff()
    L = HALF * mQ * Dot(vQ, vQ) + HALF * Dot(Om, J * Om) - mQ * g * Dot(e3, xQ)
    for mi, p in zip(masses, points):
        vp = p.t_diff()
        L = L + HALF * mi * Dot(vp, vp) - mi * g * Dot(e3, p)
    dW = Dot(xQ.delta(), (R * e3) * f) + Dot(eta, M)
    params = {"J": J_QUAD, "e3": E3, "g": 9.81, "m_Q": 1.5}
    x0: dict = dict(x_Q=[0.0, 0.0, 2.0], x_Q_dot=[0.2, -0.1, 0.3], **_hover_x0())
    for i in range(1, n + 1):
        params[f"m{i}"] = 0.05 if i < n else 0.3
        params[f"l{i}"] = 1.0 / n
        q_val = unit([np.sin(0.3 + 0.1 * i), 0.05 * i, -np.cos(0.3 + 0.1 * i)])
        x0[f"q{i}"] = q_val
        x0[f"omega_q{i}"] = tangent(q_val, [0.1 * i, 0.8 - 0.2 * i, 0.0])
    return make_system(
        "quadrotor_flexible_cable",
        "R³ × SO(3) × (S²)ⁿ",
        L,
        dW,
        SystemVariables(vectors=[xQ, *qs], matrices=[R]),
        [f, M],
        params,
        Reference("check_eom", note="Lumped-mass cable."),
        x0=x0,
        symbols=dict(xQ=xQ, R=R, qs=qs, points=points, f=f, M=M),
        scene=_load_scene(R, xQ, points[-1], rods=points[:-1]),
        options={"n_links": n},
    )


# ---------------------------------------------------------------------------
# Several quadrotors carrying one load
# ---------------------------------------------------------------------------


def _team(n):
    """Per-quadrotor symbols (mᵢ, Jᵢ, lᵢ, Rᵢ, qᵢ, fᵢ, Mᵢ) for i = 1..n."""
    ms = getScalars(" ".join(f"m{i}" for i in range(1, n + 1)), attr=["Constant"])
    ls = getScalars(" ".join(f"l{i}" for i in range(1, n + 1)), attr=["Constant"])
    Js = [_inertia(f"J{i}") for i in range(1, n + 1)]
    Rs = [SO3(f"R{i}") for i in range(1, n + 1)]
    qs = [S2(f"q{i}") for i in range(1, n + 1)]
    fs = [Scalar(f"f{i}") for i in range(1, n + 1)]
    Ms = [Vector(f"M{i}") for i in range(1, n + 1)]
    return ms, ls, Js, Rs, qs, fs, Ms


def _team_values(n, params, x0, tilt=0.5):
    for i in range(1, n + 1):
        th = 2 * np.pi * (i - 1) / n
        params[f"m{i}"] = 1.5
        params[f"l{i}"] = 1.0
        params[f"J{i}"] = J_QUAD
        qi = unit([-np.sin(tilt) * np.cos(th), -np.sin(tilt) * np.sin(th), -np.cos(tilt)])
        x0[f"q{i}"] = qi
        x0[f"R{i}"] = rot([0.05 * i, -0.03 * i, 0.1 * i])
        x0[f"omega_q{i}"] = tangent(qi, [0.1, -0.2 * i, 0.1])
        x0[f"Omega_R{i}"] = [0.1 * i, -0.1, 0.2]


def _team_scene(Rs, xQs, attach, xL, R_L=None):
    def build(layout, s):
        from geomech.viz import Box, PointMass, RigidBody, Rod, Scene, Trace

        sc = Scene(layout)
        for i, (R, xQ, a) in enumerate(zip(Rs, xQs, attach), 1):
            add_quadrotor(sc, R, xQ, label=str(i))
            sc.add(Rod(xQ, a, radius=0.006))
        if R_L is None:
            sc.add(PointMass(xL, radius=0.08, color="load", name="load"))
        else:
            sc.add(RigidBody(R_L, xL, Box(0.8, 0.6, 0.1), color="load", name="load"))
        sc.add(Trace(xL))
        return sc

    return build


@register("R³ × (SO(3) × S²)ⁿ")
def multi_quadrotor_point_mass(*, n_quads: int = 3):
    """n quadrotors carrying one point-mass load on rigid massless cables, on
    R³ × (SO(3) × S²)ⁿ.

    Load at x_L; quadrotor i at ``x_Qi = x_L − lᵢ qᵢ`` (qᵢ from the quadrotor to the load).
    ``L = ½ m_L ẋ_L·ẋ_L − m_L g e₃·x_L + Σ ( ½ mᵢ ẋ_Qi·ẋ_Qi + ½ Ωᵢ·JᵢΩᵢ − mᵢ g e₃·x_Qi )``,
    ``δW = Σ ( fᵢ (Rᵢe₃)·(δx_L − lᵢ ξᵢ × qᵢ) + Mᵢ·ηᵢ )``, inputs ``[f1, M1, …, fn, Mn]``.
    ``n_quads`` is an int from 1 to 20 (``MAX_SIZE``); any other value raises
    ``ValueError``.
    """
    n = check_count(n_quads, "multi_quadrotor_point_mass: n_quads")
    mL, g = getScalars("m_L g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    xL = Vector("x_L")
    ms, ls, Js, Rs, qs, fs, Ms = _team(n)
    vL = xL.t_diff()
    L = HALF * mL * Dot(vL, vL) - mL * g * Dot(e3, xL)
    dW = None
    xQs = []
    for m, l, J, R, q, f, M in zip(ms, ls, Js, Rs, qs, fs, Ms):
        xQ = xL - q * l
        xQs.append(xQ)
        vQ, Om = xQ.t_diff(), R.get_tangent_vector()
        L = L + HALF * m * Dot(vQ, vQ) + HALF * Dot(Om, J * Om) - m * g * Dot(e3, xQ)
        dxQ = xL.delta() - Cross(q.get_variation_vector(), q) * l
        term = Dot(dxQ, (R * e3) * f) + Dot(R.get_variation_vector(), M)
        dW = term if dW is None else dW + term
    params = {"m_L": 0.5, "g": 9.81, "e3": E3}
    x0: dict = dict(x_L=[0.0, 0.0, 0.5], x_L_dot=[0.2, -0.1, 0.3])
    _team_values(n, params, x0)
    inputs = [u for pair in zip(fs, Ms) for u in pair]
    return make_system(
        "multi_quadrotor_point_mass",
        "R³ × (SO(3) × S²)ⁿ",
        L,
        dW,
        SystemVariables(vectors=[xL, *qs], matrices=Rs),
        inputs,
        params,
        Reference("check_eom", note="Point-mass load carried by n quadrotors."),
        x0=x0,
        symbols=dict(xL=xL, qs=qs, Rs=Rs, xQs=xQs),
        scene=_team_scene(Rs, xQs, [xL] * n, xL),
        options={"n_quads": n},
    )


@register("R³ × SO(3) × (SO(3) × S²)ⁿ")
def multi_quadrotor_rigid_load(*, n_quads: int = 3):
    """n quadrotors carrying one rigid-body load on rigid massless cables, on
    R³ × SO(3) × (SO(3) × S²)ⁿ.

    Load: centre of mass x_L, attitude R_L, mass m_L, inertia J_L; cable i is attached at
    ``x_L + R_L ρᵢ`` and quadrotor i is at ``x_Qi = x_L + R_L ρᵢ − lᵢ qᵢ``.
    ``L = ½ m_L ẋ_L·ẋ_L + ½ Ω_L·J_LΩ_L − m_L g e₃·x_L
    + Σ ( ½ mᵢ ẋ_Qi·ẋ_Qi + ½ Ωᵢ·JᵢΩᵢ − mᵢ g e₃·x_Qi )``,
    ``δW = Σ ( fᵢ (Rᵢe₃)·(δx_L + R_L(η_L × ρᵢ) − lᵢ ξᵢ × qᵢ) + Mᵢ·ηᵢ )``, inputs
    ``[f1, M1, …, fn, Mn]``.
    ``n_quads`` is an int from 1 to 20 (``MAX_SIZE``); any other value raises
    ``ValueError``.
    """
    n = check_count(n_quads, "multi_quadrotor_rigid_load: n_quads")
    mL, g = getScalars("m_L g", attr=["Constant"])
    JL = _inertia("J_L")
    e3 = Vector("e3", attr=["Constant"])
    xL, RL = Vector("x_L"), SO3("R_L")
    OmL, etaL = RL.get_tangent_vector(), RL.get_variation_vector()
    rhos = [Vector(f"\\rho_{i}", attr=["Constant"]) for i in range(1, n + 1)]
    ms, ls, Js, Rs, qs, fs, Ms = _team(n)
    vL = xL.t_diff()
    L = HALF * mL * Dot(vL, vL) + HALF * Dot(OmL, JL * OmL) - mL * g * Dot(e3, xL)
    dW = None
    xQs, attach = [], []
    for m, l, J, R, q, f, M, rho in zip(ms, ls, Js, Rs, qs, fs, Ms, rhos):
        a = xL + RL * rho
        xQ = a - q * l
        attach.append(a)
        xQs.append(xQ)
        vQ, Om = xQ.t_diff(), R.get_tangent_vector()
        L = L + HALF * m * Dot(vQ, vQ) + HALF * Dot(Om, J * Om) - m * g * Dot(e3, xQ)
        dxQ = xL.delta() + RL * Cross(etaL, rho) - Cross(q.get_variation_vector(), q) * l
        term = Dot(dxQ, (R * e3) * f) + Dot(R.get_variation_vector(), M)
        dW = term if dW is None else dW + term
    params = {"m_L": 0.8, "g": 9.81, "e3": E3, "J_L": np.diag([0.03, 0.05, 0.08])}
    for i in range(1, n + 1):
        th = 2 * np.pi * (i - 1) / n
        params[f"\\rho_{i}"] = np.array([0.4 * np.cos(th), 0.4 * np.sin(th), 0.0])
    x0: dict = dict(
        x_L=[0.0, 0.0, 0.5],
        R_L=rot([0.1, -0.05, 0.2]),
        x_L_dot=[0.2, -0.1, 0.3],
        Omega_R_L=[0.2, 0.1, -0.3],
    )
    _team_values(n, params, x0, tilt=0.3)
    inputs = [u for pair in zip(fs, Ms) for u in pair]
    return make_system(
        "multi_quadrotor_rigid_load",
        "R³ × SO(3) × (SO(3) × S²)ⁿ",
        L,
        dW,
        SystemVariables(vectors=[xL, *qs], matrices=[RL, *Rs]),
        inputs,
        params,
        Reference("check_eom", note="Rigid-body load carried by n quadrotors."),
        x0=x0,
        symbols=dict(xL=xL, RL=RL, qs=qs, Rs=Rs, xQs=xQs, attach=attach),
        scene=_team_scene(Rs, xQs, attach, xL, R_L=RL),
        options={"n_quads": n},
    )


__all__ = [
    "quadrotor",
    "quadrotor_offset_com",
    "quadrotor_suspended_load",
    "quadrotor_elastic_load",
    "quadrotor_flexible_cable",
    "multi_quadrotor_point_mass",
    "multi_quadrotor_rigid_load",
]
