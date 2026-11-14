"""Systems on R, R³, S² and their products: point masses, springs and pendulums.

Gravity acts along −e₃ (``V = m g e₃·x``); a pendulum's ``q`` points from its pivot to its
mass, so ``q = −e₃`` hangs down. Each docstring states the Lagrangian and the virtual work.
"""

from __future__ import annotations

import numpy as np

from geomech.core.base import S2, Matrix, Scalar, Vector, ZeroVector, getScalars
from geomech.core.operations import Dot, Hat, MVMul
from geomech.dynamics.variables import SystemVariables
from geomech.systems._base import (
    E3,
    HALF,
    Reference,
    check_count,
    make_system,
    register,
    tangent,
    unit,
)

LEE = "lee2018global"

#: An initial cable or link direction 1 rad from hanging (``q = −e₃`` is down).
Q0 = unit([np.sin(1.0), 0.0, -np.cos(1.0)])
W0 = tangent(Q0, [0.3, 1.5, -0.2])


# ---------------------------------------------------------------------------
# Scenes
# ---------------------------------------------------------------------------


def _particle_scene(point_key):
    def build(layout, s):
        from geomech.viz import PointMass, Scene, Trace

        scene = Scene(layout)
        scene.add(PointMass(s[point_key], radius=0.06))
        scene.add(Trace(s[point_key]))
        return scene

    return build


def _chain_scene(origin_key, point_keys):
    def build(layout, s):
        from geomech.viz import PointMass, Rod, Scene, Trace

        scene = Scene(layout)
        prev = s[origin_key] if origin_key else ZeroVector
        for k in point_keys:
            scene.add(Rod(prev, s[k], radius=0.01))
            scene.add(PointMass(s[k], radius=0.06))
            prev = s[k]
        scene.add(Trace(s[point_keys[-1]]))
        return scene

    return build


# ---------------------------------------------------------------------------
# R and R³
# ---------------------------------------------------------------------------


@register("R³")
def point_mass():
    """Point mass in R³ driven by a force: ``m ẍ = F``.

    ``L = ½ m ẋ·ẋ``, ``δW = δx·F``, inputs ``[F]``. Lee et al. ch. 3 (Lagrangian dynamics
    on Rⁿ, §3.2).
    """
    (m,) = getScalars("m", attr=["Constant"])
    x, F = Vector("x"), Vector("F")
    v = x.t_diff()
    L = HALF * m * Dot(v, v)
    dW = Dot(x.delta(), F)
    return make_system(
        "point_mass",
        "R³",
        L,
        dW,
        SystemVariables(vectors=[x]),
        [F],
        {"m": 1.5},
        Reference(
            LEE, section="§3.2", note="Euler–Lagrange equations on R³ with an external force."
        ),
        x0=dict(x=[0.1, -0.2, 0.3], x_dot=[1.0, 0.5, -0.4]),
        symbols=dict(m=m, x=x, F=F),
        scene=_particle_scene("x"),
    )


@register("R")
def spring_mass():
    """Mass on a linear spring, on R: ``m ẍ + k x = f``.

    ``L = ½ m ẋ² − ½ k x²``, ``δW = δx f``, inputs ``[f]``.
    """
    m, k = getScalars("m k", attr=["Constant"])
    x, f = Scalar("x"), Scalar("f")
    v = x.t_diff()
    L = HALF * m * v * v - HALF * k * x * x
    dW = x.delta() * f
    return make_system(
        "spring_mass",
        "R",
        L,
        dW,
        SystemVariables(scalars=[x]),
        [f],
        {"m": 1.5, "k": 4.0},
        Reference("elementary", note="The harmonic oscillator."),
        x0=dict(x=0.4, x_dot=0.3),
        symbols=dict(m=m, k=k, x=x, f=f),
    )


# ---------------------------------------------------------------------------
# S²
# ---------------------------------------------------------------------------


@register("S²")
def spherical_pendulum():
    """Spherical pendulum on S² in (q, ω) form.

    ``L = ½ m l² ω·ω − m g l q·e₃``, ``δW = ξ·f``, inputs ``[f]`` (a moment about the
    pivot). Lee et al. §5.6.1, eq. 5.38 (``m l² ω̇ + m g l q × e₃ = 0``), from the general
    form eq. 5.14.
    """
    m, g, l = getScalars("m g l", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    q, f = S2("q"), Vector("f")
    omega, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * m * l * l * Dot(omega, omega) - m * g * l * Dot(q, e3)
    dW = Dot(xi, f)
    return make_system(
        "spherical_pendulum",
        "S²",
        L,
        dW,
        SystemVariables(vectors=[q]),
        [f],
        {"e3": E3, "g": 9.81, "l": 1.0, "m": 1.0},
        Reference(
            LEE,
            equation="eq. 5.38",
            section="§5.6.1; general form eq. 5.14",
            note="Point mass on a massless rod of length l.",
        ),
        x0=dict(q=Q0, omega_q=W0),
        symbols=dict(m=m, g=g, l=l, e3=e3, q=q, f=f, bob=q * l),
        scene=_chain_scene(None, ["bob"]),
    )


@register("S²")
def particle_on_rotating_sphere():
    """Particle on a sphere of radius L rotating about e₃ at the constant rate Ω, in the
    sphere-fixed frame, (q, ω) form.

    ``L̃ = ½ m L² ( ω·ω + 2 qᵀ S(Ωe₃) S(q) ω + Ω² qᵀ C q )`` with ``C = diag(1, 1, 0)`` (a
    constant symmetric matrix parameter), no gravity, ``δW = ξ·f``, inputs ``[f]``.
    Lee et al. §5.6.2, eq. 5.43 (the same dynamics in (q, q̇) form). The conserved energy is
    the Jacobi integral ``½ m L² (|ω|² − Ω² qᵀCq)``.
    """
    m, L_, Om = getScalars("m L \\Omega", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    C = Matrix("C", attr=["Constant", "SymmetricMatrix"])
    q, f = S2("q"), Vector("f")
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    S_a = Hat(e3 * Om)  # S(Ω e₃)
    coupling = Dot(q, MVMul(S_a, MVMul(Hat(q), om))) * 2
    L = HALF * m * L_ * L_ * (Dot(om, om) + coupling + Om * Om * Dot(q, C * q))
    q0 = unit([0.6, 0.2, 0.7])
    return make_system(
        "particle_on_rotating_sphere",
        "S²",
        L,
        Dot(xi, f),
        SystemVariables(vectors=[q]),
        [f],
        {"m": 1.0, "L": 1.0, "\\Omega": 2.0, "e3": E3, "C": np.diag([1.0, 1.0, 0.0])},
        Reference(
            LEE,
            equation="eq. 5.43",
            section="§5.6.2",
            note="Modified Lagrangian of §5.6.2.2; compared with eq. 5.43 through q̇ = ω × q.",
        ),
        x0=dict(q=q0, omega_q=tangent(q0, [0.5, -0.4, 0.3])),
        symbols=dict(q=q, C=C, Om=Om, e3=e3, f=f, m=m, L=L_, particle=q * L_),
        scene=_particle_scene("particle"),
    )


def _chain(n, name, configuration, reference):
    masses = getScalars(" ".join(f"m{i}" for i in range(1, n + 1)), attr=["Constant"])
    lengths = getScalars(" ".join(f"l{i}" for i in range(1, n + 1)), attr=["Constant"])
    (g,) = getScalars("g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    qs = [S2(f"q{i}") for i in range(1, n + 1)]
    fs = [Vector(f"f{i}") for i in range(1, n + 1)]
    points, x = [], None
    for q, l in zip(qs, lengths):
        x = q * l if x is None else x + q * l
        points.append(x)
    KE = sum(HALF * m * Dot(p.t_diff(), p.t_diff()) for m, p in zip(masses, points))
    PE = sum(m * g * Dot(p, e3) for m, p in zip(masses, points))
    dW = sum(Dot(q.get_variation_vector(), f) for q, f in zip(qs, fs))
    params = {"g": 9.81, "e3": E3}
    x0 = {}
    for i in range(1, n + 1):
        params[f"m{i}"] = 1.0 - 0.2 * (i - 1) / n
        params[f"l{i}"] = 1.0 - 0.3 * (i - 1) / n
        qi = unit([np.sin(0.4 + 0.3 * i), 0.1 * i, -np.cos(0.4 + 0.3 * i)])
        x0[f"q{i}"] = qi
        x0[f"omega_q{i}"] = tangent(qi, [0.2 * i, 1.0 - 0.3 * i, 0.1])
    symbols = {f"x{i}": p for i, p in enumerate(points, 1)}
    symbols.update(g=g, e3=e3, qs=qs, inputs=fs)
    return make_system(
        name,
        configuration,
        KE - PE,
        dW,
        SystemVariables(vectors=qs),
        fs,
        params,
        reference,
        x0=x0,
        symbols=symbols,
        scene=_chain_scene(None, list(symbols)[:n]),
        options={"n_links": n},
    )


@register("S² × S²")
def double_spherical_pendulum():
    """Double spherical pendulum on S² × S²: two point masses on massless links.

    ``x₁ = l₁ q₁``, ``x₂ = l₁ q₁ + l₂ q₂``; ``L = ½ m₁ ẋ₁·ẋ₁ + ½ m₂ ẋ₂·ẋ₂ − m₁ g e₃·x₁
    − m₂ g e₃·x₂`` (``ẋ`` through ``q̇ᵢ = ωᵢ × qᵢ``), ``δW = ξ₁·f₁ + ξ₂·f₂``, inputs
    ``[f1, f2]``. Lee et al. §5.6.5, eq. 5.71. Same model as
    ``n_link_spherical_chain(n_links=2)``.
    """
    return _chain(
        2,
        "double_spherical_pendulum",
        "S² × S²",
        Reference(
            LEE,
            equation="eq. 5.71",
            section="§5.6.5",
            note="Lagrangian of §5.6.5.1; compared with the (q, ω) form eq. 5.71.",
        ),
    )


@register("(S²)ⁿ")
def n_link_spherical_chain(*, n_links: int = 3):
    """Chain of n spherical pendulums on (S²)ⁿ: point masses mᵢ at the ends of massless
    links of length lᵢ, the first link hinged at the origin.

    ``xᵢ = Σ_{j≤i} lⱼ qⱼ``; ``L = Σ ½ mᵢ ẋᵢ·ẋᵢ − Σ mᵢ g e₃·xᵢ``, ``δW = Σ ξᵢ·fᵢ``, inputs
    ``[f1, …, fn]``. Lee et al. ch. 5 (eq. 5.14 with the Lagrangian of eq. 5.22 for
    n links); verified by ``check_eom``, and for n = 2 against eq. 5.71.
    ``n_links`` is an int from 1 to 20 (``MAX_SIZE``); any other value raises
    ``ValueError``.
    """
    n = check_count(n_links, "n_link_spherical_chain: n_links")
    return _chain(
        n,
        "n_link_spherical_chain",
        "(S²)ⁿ",
        Reference(
            LEE,
            equation="eq. 5.71" if n == 2 else "",
            section="§5.3, eq. 5.14; §5.6.5 for n = 2",
            note="n-link generalisation of the double spherical pendulum of §5.6.5; compared"
            " with eq. 5.71 at n = 2.",
        ),
    )


@register("R × S²")
def elastic_spherical_pendulum():
    """Spherical pendulum on an elastic rod, on R × S²: length l and direction q.

    ``L = ½ m (l̇² + l² ω·ω) − m g l q·e₃ − ½ k (l − l₀)²``, ``δW = δl f_l + ξ·f_q``,
    inputs ``[f_l, f_q]``.
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
    return make_system(
        "elastic_spherical_pendulum",
        "R × S²",
        KE - PE,
        dW,
        SystemVariables(scalars=[l], vectors=[q]),
        [f_l, f_q],
        {"e3": E3, "g": 9.81, "k": 40.0, "l0": 1.0, "m": 1.0},
        Reference(
            "check_eom",
            note="Point mass on a linear spring of rest length l0 hinged at the origin.",
        ),
        x0=dict(l=1.1, q=Q0, l_dot=0.2, omega_q=W0),
        symbols=dict(m=m, g=g, k=k, l0=l0, e3=e3, l=l, q=q, f_l=f_l, f_q=f_q, bob=q * l),
        scene=_chain_scene(None, ["bob"]),
    )


@register("R × S²")
def cart_spherical_pendulum():
    """Three-dimensional cart-pole on R × S²: a cart of mass M on a rail along the unit
    vector e₁ carrying a spherical pendulum (point mass m, rod length l).

    Cart at ``s e₁``, bob at ``x_b = s e₁ + l q``;
    ``L = ½ M ṡ² + ½ m ẋ_b·ẋ_b − (M + m) g s e₁·e₃ − m g l q·e₃``, ``δW = δs f``, inputs
    ``[f]`` (the force on the cart; the pendulum is unactuated). The rail may be inclined
    (``e₁·e₃ ≠ 0``); the default ``e₁ = (1, 0, 0)`` is horizontal.
    """
    M, m, g, l = getScalars("M m g l", attr=["Constant"])
    e1 = Vector("e1", attr=["Constant", "UnitNorm"])
    e3 = Vector("e3", attr=["Constant"])
    s, q, f = Scalar("s"), S2("q"), Scalar("f")
    cart = e1 * s
    bob = cart + q * l
    sd = s.t_diff()
    vb = bob.t_diff()
    PE = (M + m) * g * s * Dot(e1, e3) + m * g * l * Dot(q, e3)
    L = HALF * M * sd * sd + HALF * m * Dot(vb, vb) - PE
    dW = s.delta() * f
    q0 = unit([np.sin(0.6), 0.2, np.cos(0.6)])  # near upright
    return make_system(
        "cart_spherical_pendulum",
        "R × S²",
        L,
        dW,
        SystemVariables(scalars=[s], vectors=[q]),
        [f],
        {"M": 2.0, "m": 0.5, "g": 9.81, "l": 0.8, "e1": np.array([1.0, 0.0, 0.0]), "e3": E3},
        Reference("check_eom", note="Cart on a straight rail with a spherical pendulum."),
        x0=dict(s=0.0, q=q0, s_dot=0.3, omega_q=tangent(q0, [0.4, -0.6, 0.0])),
        symbols=dict(cart=cart, bob=bob, s=s, q=q, f=f),
        scene=_cart_scene,
    )


def _cart_scene(layout, s):
    from geomech.viz import PointMass, Rod, Scene, Trace

    scene = Scene(layout)
    scene.add(PointMass(s["cart"], radius=0.1, color="link1", name="cart"))
    scene.add(Rod(s["cart"], s["bob"], radius=0.01))
    scene.add(PointMass(s["bob"], radius=0.06))
    scene.add(Trace(s["bob"]))
    return scene


__all__ = [
    "point_mass",
    "spring_mass",
    "spherical_pendulum",
    "particle_on_rotating_sphere",
    "double_spherical_pendulum",
    "n_link_spherical_chain",
    "elastic_spherical_pendulum",
    "cart_spherical_pendulum",
]
