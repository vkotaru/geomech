"""Systems on S¹ and its products: planar pendulums.

Each pendulum swings in the vertical plane normal to the unit axis ``e2`` (default
``(0, 1, 0)``), so its link direction ``q`` lies on the circle S¹ = {q : |q| = 1, q·e₂ = 0}
with ``q̇ = ω e₂ × q`` and a scalar angular velocity ω (Lee et al. ch. 4, with the planar
rotation S written as ``e₂ × ·``). Gravity acts along −e₃ (``V = m g e₃·x``); ``q = −e₃``
hangs down. Each docstring states the Lagrangian and the virtual work.
"""

from __future__ import annotations

import numpy as np

from geomech.core.base import Scalar, Vector, getScalars
from geomech.core.operations import Dot
from geomech.dynamics.variables import SystemVariables
from geomech.systems._base import E3, HALF, Reference, make_system, register
from geomech.systems.spheres import _chain_scene

LEE = "lee2018global"

#: The plane normal of the planar pendulums (the S¹ axis), e₂.
E2 = np.array([0.0, 1.0, 0.0])


def _circle(angle):
    """The point of the circle normal to e₂ at *angle* from hanging (``−e₃``)."""
    return np.array([np.sin(angle), 0.0, -np.cos(angle)])


def _axis():
    return Vector("e2", attr=["Constant", "UnitNorm"])


def _s1(name, axis):
    """``S1(name, axis=axis)``. ``S1`` is imported here, not at module level, so that
    ``geomech.systems`` imports against an older ``geomech.core`` (``tests/
    test_viz_adversarial.py`` loads it next to the core of an earlier commit)."""
    from geomech.core.base import S1

    return S1(name, axis=axis)


@register("S¹")
def planar_pendulum():
    """Planar pendulum on S¹ in (q, ω) form.

    ``L = ½ m l² ω² − m g l q·e₃``, ``δW = γ τ``, inputs ``[τ]`` (a moment about the pivot
    axis e₂), with ``q̇ = ω e₂ × q`` and ``δq = γ e₂ × q``. Lee et al. §4.6.1, eq. 4.33
    (``m l² ω̇ + m g l e₁ᵀq = 0`` in the book's planar frame), from the general form
    eq. 4.12.
    """
    m, g, l = getScalars("m g l", attr=["Constant"])
    e2, e3 = _axis(), Vector("e3", attr=["Constant"])
    q, tau = _s1("q", e2), Scalar("\\tau")
    omega, gamma = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * m * l * l * omega * omega - m * g * l * Dot(q, e3)
    dW = gamma * tau
    return make_system(
        "planar_pendulum",
        "S¹",
        L,
        dW,
        SystemVariables(vectors=[q]),
        [tau],
        {"e2": E2, "e3": E3, "g": 9.81, "l": 1.0, "m": 1.0},
        Reference(
            LEE,
            equation="eq. 4.33",
            section="§4.6.1; general form eq. 4.12",
            note="Point mass on a massless rod of length l swinging in the plane normal to e₂.",
        ),
        x0=dict(q=_circle(1.0), omega_q=1.5),
        symbols=dict(m=m, g=g, l=l, e2=e2, e3=e3, q=q, tau=tau, bob=q * l),
        scene=_chain_scene(None, ["bob"]),
    )


@register("S¹ × S¹")
def double_planar_pendulum():
    """Double planar pendulum on (S¹)²: two point masses on massless links in the plane
    normal to e₂.

    ``x₁ = l₁ q₁``, ``x₂ = l₁ q₁ + l₂ q₂``; ``L = ½ m₁ ẋ₁·ẋ₁ + ½ m₂ ẋ₂·ẋ₂ − m₁ g e₃·x₁
    − m₂ g e₃·x₂`` (``ẋ`` through ``q̇ᵢ = ωᵢ e₂ × qᵢ``), ``δW = γ₁ τ₁ + γ₂ τ₂`` (moments
    about e₂ on each link), inputs ``[τ1, τ2]``. Lee et al. §4.6.4, eq. 4.62; the velocity
    terms of eq. 4.62 as printed (``ω₂(ω₂ − ω₁)``, ``ω₁(ω₁ − ω₂)``) omit ``qᵢᵀ S ∂L̃/∂qᵢ``
    of the kinetic energy, and the equations here have ``ω₂²`` and ``ω₁²``, those of
    the textbook double pendulum (see ``tests/test_s1_systems.py``).
    """
    m1, m2, l1, l2, g = getScalars("m1 m2 l1 l2 g", attr=["Constant"])
    e2, e3 = _axis(), Vector("e3", attr=["Constant"])
    q1, q2 = _s1("q1", e2), _s1("q2", e2)
    tau1, tau2 = Scalar("\\tau_1"), Scalar("\\tau_2")
    x1 = q1 * l1
    x2 = x1 + q2 * l2
    v1, v2 = x1.t_diff(), x2.t_diff()
    KE = HALF * m1 * Dot(v1, v1) + HALF * m2 * Dot(v2, v2)
    PE = m1 * g * Dot(x1, e3) + m2 * g * Dot(x2, e3)
    dW = q1.get_variation_vector() * tau1 + q2.get_variation_vector() * tau2
    return make_system(
        "double_planar_pendulum",
        "S¹ × S¹",
        KE - PE,
        dW,
        SystemVariables(vectors=[q1, q2]),
        [tau1, tau2],
        {"e2": E2, "e3": E3, "g": 9.81, "m1": 1.0, "m2": 0.8, "l1": 1.0, "l2": 0.7},
        Reference(
            LEE,
            equation="eq. 4.62",
            section="§4.6.4",
            note="Lagrangian of §4.6.4.1; compared with the (q, ω) form eq. 4.62, whose"
            " velocity terms are corrected to ω₂² and ω₁².",
        ),
        x0=dict(q1=_circle(0.8), q2=_circle(1.6), omega_q1=0.5, omega_q2=-1.2),
        symbols=dict(m1=m1, m2=m2, l1=l1, l2=l2, g=g, e2=e2, e3=e3, q1=q1, q2=q2, x1=x1, x2=x2),
        scene=_chain_scene(None, ["x1", "x2"]),
    )
