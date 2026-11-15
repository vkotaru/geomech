"""Systems on S¹ and its products for the pipeline tests (``tests/test_s1_manifold.py``) and
the Euler–Lagrange oracle (``tests/test_eom_oracle.py``).

Each builder returns a ``tests.vbl_systems.System``. The planar frame of Lee et al. ch. 4 is
embedded as the plane normal to a constant unit axis ``n``: the book's ``S`` is ``n × ·``.
"""

from __future__ import annotations

from geomech import S1, S2, Dot, Scalar, SystemVariables, Vector, getScalars, systems
from tests.vbl_systems import HALF, System, _from_library


def _axis(name="n"):
    return Vector(name, attr=["Constant", "UnitNorm"])


def planar_pendulum():
    """``geomech.systems.planar_pendulum``: Lee et al. eq. 4.33."""
    return _from_library(systems.planar_pendulum(), "m g l e2 e3 q tau")


def planar_pendulum_qdot():
    """The planar pendulum written with ``q̇``: ``L = ½ m l² q̇·q̇ − m g l e₃·q`` with
    ``q̇ = ω n × q``; the S¹ identities reduce ``q̇·q̇`` to ``ω²``. Lee et al. eq. 4.33."""
    m, g, l = getScalars("m g l", attr=["Constant"])
    n, e3 = _axis(), Vector("e3", attr=["Constant"])
    q, tau = S1("q", axis=n), Scalar("\\tau")
    qd = q.t_diff()
    L = HALF * m * l * l * Dot(qd, qd) - m * g * l * Dot(e3, q)
    dW = q.get_variation_vector() * tau
    return System(L, dW, SystemVariables(vectors=[q]), [tau], dict(m=m, g=g, l=l, n=n, q=q))


def double_planar_pendulum():
    """``geomech.systems.double_planar_pendulum``: Lee et al. eq. 4.62."""
    return _from_library(systems.double_planar_pendulum(), "m1 m2 l1 l2 g e2 e3 q1 q2")


def rotating_hoop():
    """Particle on a hoop rotating at the rate Ω about a diameter, in the hoop frame (Lee
    et al. §4.6.2.2): ``L̃ = ½ m L² (ω² + Ω² (e₁ᵀq)²)``, no potential, ``δW = γ τ``;
    eq. 4.40: ``m L² ω̇ + m L² Ω² qᵀe₂ e₁ᵀq = 0`` with ``e₂ = S e₁ = n × e₁``."""
    m, L_, Om = getScalars("m L \\Omega", attr=["Constant"])
    n, e1 = _axis(), Vector("e1", attr=["Constant", "UnitNorm"])
    q, tau = S1("q", axis=n), Scalar("\\tau")
    w = q.get_tangent_vector()
    L = HALF * m * L_ * L_ * (w * w + Om * Om * Dot(e1, q) * Dot(e1, q))
    dW = q.get_variation_vector() * tau
    return System(L, dW, SystemVariables(vectors=[q]), [tau], dict(m=m, L=L_, Om=Om, e1=e1, q=q))


def elastic_planar_pendulums():
    """Two elastically connected planar pendulums (Lee et al. §4.6.3.2):
    ``L̃ = ½ m L² ω₁² + ½ m L² ω₂² − κ (1 − q₁ᵀq₂)``, ``δW = γ₁τ₁ + γ₂τ₂``; eqs. 4.49–4.50:
    ``m L² ω̇₁ + κ q₁ᵀSq₂ = 0``, ``m L² ω̇₂ + κ q₂ᵀSq₁ = 0``."""
    m, L_, k = getScalars("m L \\kappa", attr=["Constant"])
    n = _axis()
    q1, q2 = S1("q_1", axis=n), S1("q_2", axis=n)
    t1, t2 = Scalar("\\tau_1"), Scalar("\\tau_2")
    w1, w2 = q1.get_tangent_vector(), q2.get_tangent_vector()
    L = HALF * m * L_ * L_ * w1 * w1 + HALF * m * L_ * L_ * w2 * w2 - k * (1 - Dot(q1, q2))
    dW = q1.get_variation_vector() * t1 + q2.get_variation_vector() * t2
    variables = SystemVariables(vectors=[q1, q2])
    return System(L, dW, variables, [t1, t2], dict(m=m, L=L_, k=k, n=n, q1=q1, q2=q2))


def cart_planar_pendulum():
    """Cart-pole on R × S¹: cart at ``s e₁``, bob at ``s e₁ + l q`` with ``q`` on the circle
    normal to ``n``; ``L = ½ M ṡ² + ½ m |ẋ_b|² − m g l e₃·q``, ``δW = δs f + γ τ``. Mixes a
    scalar configuration variable with an S¹ point (``check_eom`` only)."""
    M, m, g, l = getScalars("M m g l", attr=["Constant"])
    n, e1, e3 = _axis(), Vector("e1", attr=["Constant"]), Vector("e3", attr=["Constant"])
    s, q = Scalar("s"), S1("q", axis=n)
    f, tau = Scalar("f"), Scalar("\\tau")
    xb = e1 * s + q * l
    vb = xb.t_diff()
    sd = s.t_diff()
    L = HALF * M * sd * sd + HALF * m * Dot(vb, vb) - m * g * l * Dot(e3, q)
    dW = s.delta() * f + q.get_variation_vector() * tau
    variables = SystemVariables(scalars=[s], vectors=[q])
    return System(L, dW, variables, [f, tau], dict(s=s, q=q, n=n))


def cart_planar_pendulum_library():
    """``geomech.systems.cart_planar_pendulum``: the planar cart-pole on R × S¹ with the
    rail force only (the builder above also actuates the pendulum)."""
    return _from_library(systems.cart_planar_pendulum(), "M m g l e1 e2 e3 s q f")


def gyroscopic_s1():
    """An S¹ point whose momentum depends on the configuration: ``L = ½ J ω² + c ω a·q −
    k e₃·q``, ``δW = γ τ``. ``∂L/∂ω = J ω + c a·q``, so integration by parts differentiates
    the prefactor ``a·q``; the two ``c`` terms then cancel (``c ω a·q dt`` is a closed
    1-form on the circle, which exerts no force), leaving ``E = −J ω̇ − k e₃·(n × q) + τ``."""
    J, c, k = getScalars("J c k", attr=["Constant"])
    n, a, e3 = _axis(), Vector("a", attr=["Constant"]), Vector("e3", attr=["Constant"])
    q, tau = S1("q", axis=n), Scalar("\\tau")
    w = q.get_tangent_vector()
    L = HALF * J * w * w + c * w * Dot(a, q) - k * Dot(e3, q)
    dW = q.get_variation_vector() * tau
    return System(L, dW, SystemVariables(vectors=[q]), [tau], dict(q=q, n=n, a=a))


def s1_hanging_from_s2():
    """A planar pendulum hung from the bob of a spherical pendulum, on S² × S¹: ``x₁ = l₁ p``,
    ``x₂ = x₁ + l₂ q`` with ``p ∈ S²`` and ``q`` on the circle normal to ``n``;
    ``L = ½ m₁|ẋ₁|² + ½ m₂|ẋ₂|² − g e₃·(m₁x₁ + m₂x₂)``, ``δW = ξ·f + γ τ`` (``check_eom``
    only: the S² and S¹ substitutions in one δS)."""
    m1, m2, l1, l2, g = getScalars("m_1 m_2 l_1 l_2 g", attr=["Constant"])
    n, e3 = _axis(), Vector("e3", attr=["Constant"])
    p, q = S2("p"), S1("q", axis=n)
    f, tau = Vector("f"), Scalar("\\tau")
    x1 = p * l1
    x2 = x1 + q * l2
    v1, v2 = x1.t_diff(), x2.t_diff()
    L = HALF * m1 * Dot(v1, v1) + HALF * m2 * Dot(v2, v2) - g * Dot(e3, x1 * m1 + x2 * m2)
    dW = Dot(p.get_variation_vector(), f) + q.get_variation_vector() * tau
    variables = SystemVariables(vectors=[p, q])
    return System(L, dW, variables, [f, tau], dict(p=p, q=q, n=n))


def spacecraft_reaction_wheels_s1():
    """``geomech.systems.spacecraft_reaction_wheels_s1``: SO(3) × (S¹)³."""
    return _from_library(systems.spacecraft_reaction_wheels_s1(), "R J_b b q tau")


def spacecraft_solar_arrays():
    """``geomech.systems.spacecraft_solar_arrays``: SO(3) × (S¹)³ × (S¹)²."""
    return _from_library(systems.spacecraft_solar_arrays(), "R J_b b q tau a n mu")


SYSTEMS = [
    planar_pendulum,
    planar_pendulum_qdot,
    double_planar_pendulum,
    rotating_hoop,
    elastic_planar_pendulums,
    cart_planar_pendulum,
    cart_planar_pendulum_library,
    gyroscopic_s1,
    s1_hanging_from_s2,
    spacecraft_reaction_wheels_s1,
    spacecraft_solar_arrays,
]
