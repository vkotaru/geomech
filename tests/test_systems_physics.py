"""Every Lagrangian and virtual work of ``geomech.systems`` agrees with an independent model
of the system's bodies.

For each system (parametric ones at sizes 1, 2 and 3) the model below is written in NumPy
from the system's physical description, without the library's expressions:

- the positions of the point masses and centres of mass as functions of the
  configuration; the kinetic energy is ``Σ ½ m |ẋ|²`` with ``ẋ`` by central differences
  along the trajectory through the state (``q(t) = exp(t ω̂) q`` on S²,
  ``q(t) = exp(t ω ê) q`` on S¹, ``R(t) = R exp(t Ω̂)``, ``x(t) = x + t ẋ``), plus
  ``½ Ωᵀ J_c Ω`` per rigid body with ``J_c`` about its centre of mass (a pivot inertia is
  shifted to the centre of mass with the parallel-axis theorem);
- the potential: ``Σ m g (height)`` with gravity along −e₃, plus the springs;
- the virtual work: each force times the central difference of its application point
  along a random admissible variation (``δq = ξ × q``, ``δq = γ e × q``, ``δR = R η̂``),
  each body moment ``M·η`` and each moment about a pivot ``f·ξ`` (the work of a moment on
  the rotation ξ).

At random admissible states and the default parameters, ``L`` and ``δW`` evaluated with
``geomech.numeric.Evaluator`` agree with the model to 1e-8 relative to ``max(1, |model|)``.
"""

from __future__ import annotations

import numpy as np
import pytest

from geomech.core.base import S1, S2, SO3, Scalar
from geomech.numeric import Evaluator, expm_so3, sample_state
from geomech.systems import get, names

H = 1e-6
TOL = 1e-8
SEEDS = (0, 1, 2)
UP = np.array([0.0, 0.0, 1.0])  # the physical vertical, independent of the param "e3"
SIZED = {
    "n_link_spherical_chain": "n_links",
    "n_link_rigid_chain": "n_links",
    "quadrotor_flexible_cable": "n_links",
    "multi_quadrotor_point_mass": "n_quads",
    "multi_quadrotor_rigid_load": "n_quads",
}


def col(v):
    return np.asarray(v, dtype=float).reshape(3)


def kinematic_keys(variables):
    """``(variable, configuration key, velocity key)`` per configuration variable."""
    out = []
    for v in variables.scalars:
        out.append((v, str(v), str(v.t_diff())))
    for v in variables.vectors:
        vel = v.get_tangent_vector() if isinstance(v, (S1, S2)) else v.t_diff()
        out.append((v, str(v), str(vel)))
    for v in variables.matrices:
        out.append((v, str(v), str(v.get_tangent_vector())))
    return out


def flow(keys, env, t, direction=None):
    """The configuration at time *t* along the velocity (or along *direction*)."""
    cfg = {}
    for var, kq, kv in keys:
        d = env[kv] if direction is None else direction[kq]
        if isinstance(var, Scalar):
            cfg[kq] = float(env[kq]) + t * float(d)
        elif isinstance(var, S1):  # the rotation about the axis by the angle t·ω (or t·γ)
            cfg[kq] = expm_so3(t * float(d) * col(env[str(var.axis)])) @ col(env[kq])
        elif isinstance(var, S2):
            cfg[kq] = expm_so3(t * col(d)) @ col(env[kq])
        elif isinstance(var, SO3):
            cfg[kq] = np.asarray(env[kq]) @ expm_so3(t * col(d))
        else:
            cfg[kq] = col(env[kq]) + t * col(d)
    return cfg


class Model:
    """The bodies of one system: ``points(c, t)`` → [(mass, position)], ``bodies`` →
    [(J about the centre of mass, Ω key)], ``springs(c)`` → potential, ``forces(c, u)`` →
    [(force, point function)], ``moments(u)`` → [(moment, configuration key)],
    ``scalar_forces(u)`` → [(generalized force, scalar key)]."""

    def __init__(self, s, name):
        self.p = s.params
        n = s.options.get("n_links", s.options.get("n_quads", 0))
        self.points = lambda c, t: []
        self.bodies: list = []
        self.springs = lambda c: 0.0
        self.forces = lambda c, u: []
        self.moments = lambda u: []
        self.scalar_forces = lambda u: []
        # kinetic energy and virtual work that need the attitudes themselves:
        # extra_T(cp, cm) and extra_W(c0, vp, vm, u), with cp/cm the configurations at ±H
        # along the velocity and vp/vm along the variation
        self.extra_T = lambda cp, cm: 0.0
        self.extra_W = lambda c0, vp, vm, u: 0.0
        self.g = self.p.get("g", 0.0)
        getattr(self, name)(n)

    def v(self, k):
        val = self.p[k]
        return float(val) if np.isscalar(val) else np.asarray(val, dtype=float)

    def body_moment(self, *pairs):
        self.moments = lambda u: [(u[m], k) for m, k in pairs]

    # ---- R, R³, S¹, S²
    def point_mass(self, n):
        self.points = lambda c, t: [(self.v("m"), c["x"])]
        self.forces = lambda c, u: [(u["F"], lambda c: c["x"])]

    def spring_mass(self, n):
        self.points = lambda c, t: [(self.v("m"), c["x"] * np.array([1.0, 0, 0]))]
        self.springs = lambda c: 0.5 * self.v("k") * c["x"] ** 2
        self.scalar_forces = lambda u: [(u["f"], "x")]

    def planar_pendulum(self, n):
        # τ is the moment about the pivot axis: its work on the rotation γ e₂ is τ γ
        self.points = lambda c, t: [(self.v("m"), self.v("l") * c["q"])]
        self.scalar_forces = lambda u: [(u["\\tau"], "q")]

    def double_planar_pendulum(self, n):
        self.points = self._chain(2)
        self.scalar_forces = lambda u: [(u["\\tau_1"], "q1"), (u["\\tau_2"], "q2")]

    def spherical_pendulum(self, n):
        self.points = lambda c, t: [(self.v("m"), self.v("l") * c["q"])]
        self.body_moment(("f", "q"))

    def particle_on_rotating_sphere(self, n):
        # sphere-fixed coordinates q; inertial position Rz(Ω t) L q; no gravity
        def pts(c, t):
            Rz = expm_so3(self.v("\\Omega") * t * UP)
            return [(self.v("m"), self.v("L") * Rz @ c["q"])]

        self.points = pts
        self.g = 0.0
        self.body_moment(("f", "q"))

    def _chain(self, n, origin=None):
        def pts(c, t):
            out, x = [], (np.zeros(3) if origin is None else c[origin])
            for i in range(1, n + 1):
                x = x + self.v(f"l{i}") * c[f"q{i}"]
                out.append((self.v(f"m{i}"), x))
            return out

        return pts

    def double_spherical_pendulum(self, n):
        self.n_link_spherical_chain(2)

    def n_link_spherical_chain(self, n):
        self.points = self._chain(n)
        self.body_moment(*[(f"f{i}", f"q{i}") for i in range(1, n + 1)])

    def elastic_spherical_pendulum(self, n):
        self.points = lambda c, t: [(self.v("m"), c["l"] * c["q"])]
        self.springs = lambda c: 0.5 * self.v("k") * (c["l"] - self.v("l0")) ** 2
        # f_l pushes the bob along the rod; f_q is a moment about the pivot
        self.forces = lambda c, u: [(u["f_l"] * c["q"], lambda c: c["l"] * c["q"])]
        self.body_moment(("f_q", "q"))

    def cart_spherical_pendulum(self, n):
        e1 = self.v("e1")
        self.points = lambda c, t: [
            (self.v("M"), c["s"] * e1),
            (self.v("m"), c["s"] * e1 + self.v("l") * c["q"]),
        ]
        self.forces = lambda c, u: [(u["f"] * e1, lambda c: c["s"] * e1)]

    def cart_planar_pendulum(self, n):
        # the same bodies; q stays in the plane normal to e₂ along the S¹ flow
        self.cart_spherical_pendulum(n)

    # ---- rigid bodies
    def free_rigid_body(self, n):
        self.bodies = [(self.v("J"), "\\Omega_{R}")]
        self.body_moment(("M", "R"))

    def _pivoted(self, r_body):
        """A body hinged at the origin, centre of mass at R r_body, J given about the pivot."""
        m, J = self.v("m"), self.v("J")
        J_c = J - m * ((r_body @ r_body) * np.eye(3) - np.outer(r_body, r_body))
        self.points = lambda c, t: [(m, c["R"] @ r_body)]
        self.bodies = [(J_c, "\\Omega_{R}")]
        self.body_moment(("M", "R"))

    def heavy_top(self, n):
        self._pivoted(self.v("l") * UP)

    def rigid_pendulum(self, n):
        self._pivoted(self.v("\\rho"))

    def double_rigid_pendulum(self, n):
        self.points = lambda c, t: [
            (self.v("m1"), c["R1"] @ self.v("\\rho_1")),
            (self.v("m2"), c["R1"] @ self.v("l_1") + c["R2"] @ self.v("\\rho_2")),
        ]
        self.bodies = [(self.v("J1"), "\\Omega_{R1}"), (self.v("J2"), "\\Omega_{R2}")]
        self.body_moment(("M1", "R1"), ("M2", "R2"))

    def n_link_rigid_chain(self, n):
        """Joints p₁ = 0, pₖ₊₁ = pₖ + Rₖlₖ; centres of mass pₖ + Rₖρₖ; Jₖ about the centre
        of mass; joint torque τₖ acts on link k and −τₖ on link k − 1 (inertial frame),
        its work taken on finite-difference rotation variations of the attitudes."""

        def pts(c, t):
            out, p = [], np.zeros(3)
            for k in range(1, n + 1):
                out.append((self.v(f"m{k}"), p + c[f"R{k}"] @ self.v(f"\\rho_{k}")))
                if k < n:
                    p = p + c[f"R{k}"] @ self.v(f"l_{k}")
            return out

        self.points = pts
        self.bodies = [(self.v(f"J{k}"), f"\\Omega_{{R{k}}}") for k in range(1, n + 1)]

        def spatial_variation(Rp, Rm):
            Om = ((Rp - Rm) / (2 * H)) @ (0.5 * (Rp + Rm)).T
            return np.array([Om[2, 1], Om[0, 2], Om[1, 0]])

        def extra_W(c0, vp, vm, u):
            w = {k: spatial_variation(vp[f"R{k}"], vm[f"R{k}"]) for k in range(1, n + 1)}
            W = 0.0
            for k in range(1, n + 1):
                tau = u[f"\\tau_{k}"]
                W += tau @ w[k] - (tau @ w[k - 1] if k > 1 else 0.0)
            return W

        self.extra_W = extra_W

    def two_rigid_bodies_spherical_joint(self, n):
        self.points = lambda c, t: [
            (self.v("m1"), c["x"]),
            (self.v("m2"), c["x"] + c["R1"] @ self.v("\\rho_1") + c["R2"] @ self.v("\\rho_2")),
        ]
        self.bodies = [(self.v("J1"), "\\Omega_{R1}"), (self.v("J2"), "\\Omega_{R2}")]
        self.forces = lambda c, u: [(u["F"], lambda c: c["x"])]
        self.body_moment(("M1", "R1"), ("M2", "R2"))

    def rigid_body(self, n):
        self.points = lambda c, t: [(self.v("m"), c["x"])]
        self.bodies = [(self.v("J"), "\\Omega_{R}")]
        self.forces = lambda c, u: [(u["F"], lambda c: c["x"])]
        self.body_moment(("M", "R"))

    # ---- multibody
    def spacecraft_reaction_wheels(self, n):
        """The bus and three wheels as separate rigid bodies. Each wheel i: mass m_w at
        ρᵢ = 0.3 bᵢ, attitude ``R expm(θᵢ b̂ᵢ)``, inertia ``Iᵢ bᵢbᵢᵀ + J_t (I − bᵢbᵢᵀ)`` in
        its own frame. ``J_b`` is documented as the bus plus the wheels' masses and
        transverse inertias, so the bus alone has
        ``J_bus = J_b − Σ (J_t (I − bᵢbᵢᵀ) + m_w (|ρᵢ|² I − ρᵢρᵢᵀ))``. ``m_w`` and ``J_t``
        are chosen here only; a correct L does not depend on them. Angular velocities and
        rotation variations come from central differences of the attitudes. The motor
        applies +τᵢ bᵢ to wheel i and −τᵢ bᵢ to the bus."""
        m_w, J_t = 0.7, 0.013
        bs = [self.v(f"b{i}") for i in (1, 2, 3)]
        Is = [self.v(f"I_{i}") for i in (1, 2, 3)]
        rhos = [0.3 * b for b in bs]
        J_bus = self.v("J_b") - sum(
            J_t * (np.eye(3) - np.outer(b, b)) + m_w * ((r @ r) * np.eye(3) - np.outer(r, r))
            for b, r in zip(bs, rhos)
        )

        def wheel(c, i):
            return c["R"] @ expm_so3(self.wheel_angle(c, i) * bs[i - 1])

        def body_rate(Rp, Rm):  # body-frame angular velocity (or variation) by differences
            Rd = (Rp - Rm) / (2 * H)
            W = 0.5 * (Rp + Rm)
            Om = W.T @ Rd
            return np.array([Om[2, 1], Om[0, 2], Om[1, 0]])

        def extra_T(cp, cm):
            T = 0.5 * body_rate(cp["R"], cm["R"]) @ J_bus @ body_rate(cp["R"], cm["R"])
            for i, (b, I_i, r) in enumerate(zip(bs, Is, rhos), 1):
                w = body_rate(wheel(cp, i), wheel(cm, i))
                Jw = I_i * np.outer(b, b) + J_t * (np.eye(3) - np.outer(b, b))
                v = (cp["R"] @ r - cm["R"] @ r) / (2 * H)
                T += 0.5 * w @ Jw @ w + 0.5 * m_w * v @ v
            return T

        def extra_W(c0, vp, vm, u):
            eta = body_rate(vp["R"], vm["R"])
            W = 0.0
            for i, b in enumerate(bs, 1):
                tau = u[f"\\tau_{i}"]
                W += tau * b @ body_rate(wheel(vp, i), wheel(vm, i))  # +τ b on the wheel
                W -= tau * b @ eta  # −τ b on the bus
            return W

        self.extra_T, self.extra_W = extra_T, extra_W

    def wheel_angle(self, c, i):
        """θᵢ of the R³ model, or the angle about bᵢ of the S¹ point qᵢ from a fixed
        reference direction rᵢ ⊥ bᵢ (any rᵢ: the wheel is axisymmetric about bᵢ)."""
        if f"\\theta_{i}" in c:
            return c[f"\\theta_{i}"]
        b, q = self.v(f"b{i}"), col(c[f"q_{i}"])
        r = np.cross(b, [0.3, -0.5, 0.8])
        r /= np.linalg.norm(r)
        return np.arctan2(b @ np.cross(r, q), r @ q)

    def spacecraft_reaction_wheels_s1(self, n):
        """The bodies of ``spacecraft_reaction_wheels``; wheel i turned by the angle of
        its S¹ point ``q_i`` about bᵢ (``wheel_angle``)."""
        self.spacecraft_reaction_wheels(n)

    def spacecraft_solar_arrays(self, n):
        """The bodies of ``spacecraft_reaction_wheels_s1`` plus two panels. Panel k: attitude
        ``R expm(φₖ âₖ)`` with φₖ the angle about aₖ of its S¹ point ``n_k`` from a fixed
        reference normal ``n⁰ₖ ⊥ aₖ``; inertia ``Ia aₖaₖᵀ + In n⁰ₖn⁰ₖᵀ + Iw wₖwₖᵀ``
        (``wₖ = aₖ × n⁰ₖ``) in its own frame; mass ``m_pk`` at ``R dₖ aₖ`` (measured, like
        ``J_b``, from the composite centre of mass, which is at rest). The drive applies
        +μₖ aₖ to the panel and −μₖ aₖ to the bus."""
        self.spacecraft_reaction_wheels(n)
        wheels_T, wheels_W = self.extra_T, self.extra_W
        axes = [self.v(f"a{k}") for k in (1, 2)]

        def frame(k):
            a = axes[k - 1]
            n0 = np.cross(a, [0.3, -0.5, 0.8])
            n0 /= np.linalg.norm(n0)
            return a, n0

        def panel(c, k):
            a, n0 = frame(k)
            n = col(c[f"n_{k}"])
            phi = np.arctan2(a @ np.cross(n0, n), n0 @ n)
            return c["R"] @ expm_so3(phi * a)

        def J0(k):
            a, n0 = frame(k)
            w = np.cross(a, n0)
            return (
                self.v(f"Ia_{k}") * np.outer(a, a)
                + self.v(f"In_{k}") * np.outer(n0, n0)
                + self.v(f"Iw_{k}") * np.outer(w, w)
            )

        def body_rate(Rp, Rm):
            Om = (0.5 * (Rp + Rm)).T @ ((Rp - Rm) / (2 * H))
            return np.array([Om[2, 1], Om[0, 2], Om[1, 0]])

        self.points = lambda c, t: [
            (self.v(f"m_p{k}"), c["R"] @ (self.v(f"d_{k}") * axes[k - 1])) for k in (1, 2)
        ]

        def extra_T(cp, cm):
            T = wheels_T(cp, cm)
            for k in (1, 2):
                w = body_rate(panel(cp, k), panel(cm, k))
                T += 0.5 * w @ J0(k) @ w
            return T

        def extra_W(c0, vp, vm, u):
            W = wheels_W(c0, vp, vm, u)
            eta = body_rate(vp["R"], vm["R"])
            for k in (1, 2):
                a, mu = axes[k - 1], u[f"\\mu_{k}"]
                W += mu * a @ body_rate(panel(vp, k), panel(vm, k))  # +μ a on the panel
                W -= mu * a @ eta  # −μ a on the bus
            return W

        self.extra_T, self.extra_W = extra_T, extra_W

    def biped_3d(self, n):
        """Hip point mass at the top of the stance leg, legs of length l along their body
        axis e₃ (foot to hip), centres of mass a from the foot, J about the centre of mass;
        the hip torque acts +M_H on the swing leg and −M_H on the stance leg, the ankle
        torque +M_A on the stance leg, all inertial-frame moments."""
        l, a = self.v("l"), self.v("a")

        def hip(c):
            return c["R1"] @ np.array([0.0, 0.0, l])

        self.points = lambda c, t: [
            (self.v("m_H"), hip(c)),
            (self.v("m"), c["R1"] @ np.array([0.0, 0.0, a])),
            (self.v("m"), hip(c) - c["R2"] @ np.array([0.0, 0.0, l - a])),
        ]
        self.bodies = [(self.v("J"), "\\Omega_{R1}"), (self.v("J"), "\\Omega_{R2}")]

        def spatial_variation(Rp, Rm):  # inertial rotation variation by differences
            Om = ((Rp - Rm) / (2 * H)) @ (0.5 * (Rp + Rm)).T
            return np.array([Om[2, 1], Om[0, 2], Om[1, 0]])

        def extra_W(c0, vp, vm, u):
            w1 = spatial_variation(vp["R1"], vm["R1"])
            w2 = spatial_variation(vp["R2"], vm["R2"])
            return u["M_H"] @ w2 - u["M_H"] @ w1 + u["M_A"] @ w1

        self.extra_W = extra_W

    # ---- aerial
    def _thrust(self, R, point):
        return lambda c, u: [(u["f"] * c[R] @ UP, point)]

    def quadrotor(self, n):
        self.points = lambda c, t: [(self.v("m"), c["x"])]
        self.bodies = [(self.v("J"), "\\Omega_{R}")]
        self.forces = self._thrust("R", lambda c: c["x"])
        self.body_moment(("M", "R"))

    def quadrotor_offset_com(self, n):
        self.points = lambda c, t: [(self.v("m"), c["x"] + c["R"] @ self.v("\\rho"))]
        self.bodies = [(self.v("J"), "\\Omega_{R}")]
        self.forces = self._thrust("R", lambda c: c["x"])  # thrust at the rotor centre
        self.body_moment(("M", "R"))

    def _quad_load(self, length):
        def xQ(c):
            return c["x_L"] - length(c) * c["q"]

        self.points = lambda c, t: [(self.v("m_Q"), xQ(c)), (self.v("m_L"), c["x_L"])]
        self.bodies = [(self.v("J"), "\\Omega_{R}")]
        self.forces = self._thrust("R", xQ)
        self.body_moment(("M", "R"))

    def quadrotor_suspended_load(self, n):
        self._quad_load(lambda c: self.v("l"))

    def quadrotor_elastic_load(self, n):
        self._quad_load(lambda c: c["l"])
        self.springs = lambda c: 0.5 * self.v("k") * (c["l"] - self.v("l0")) ** 2

    def quadrotor_flexible_cable(self, n):
        chain = self._chain(n, origin="x_Q")
        self.points = lambda c, t: [(self.v("m_Q"), c["x_Q"]), *chain(c, t)]
        self.bodies = [(self.v("J"), "\\Omega_{R}")]
        self.forces = self._thrust("R", lambda c: c["x_Q"])
        self.body_moment(("M", "R"))

    def _team(self, n, rigid):
        def attach(c, i):
            return c["x_L"] + (c["R_L"] @ self.v(f"\\rho_{i}") if rigid else 0.0)

        def xQ(c, i):
            return attach(c, i) - self.v(f"l{i}") * c[f"q{i}"]

        self.points = lambda c, t: (
            [(self.v("m_L"), c["x_L"])] + [(self.v(f"m{i}"), xQ(c, i)) for i in range(1, n + 1)]
        )
        self.bodies = [(self.v(f"J{i}"), f"\\Omega_{{R{i}}}") for i in range(1, n + 1)]
        if rigid:
            self.bodies.append((self.v("J_L"), "\\Omega_{R_L}"))
        self.forces = lambda c, u: [
            (u[f"f{i}"] * c[f"R{i}"] @ UP, (lambda c, i=i: xQ(c, i))) for i in range(1, n + 1)
        ]
        self.body_moment(*[(f"M{i}", f"R{i}") for i in range(1, n + 1)])

    def multi_quadrotor_point_mass(self, n):
        self._team(n, rigid=False)

    def multi_quadrotor_rigid_load(self, n):
        self._team(n, rigid=True)


def errors(s, name, seed):
    """Relative errors of L and δW against the model at one random admissible state."""
    rng = np.random.default_rng(seed)
    keys = kinematic_keys(s.variables)
    # an S¹ point is drawn in the plane of its default axis; for the other systems the base
    # env changes no random draw
    env = sample_state(s.variables, list(s.inputs), rng, env=dict(s.params))
    for k, val in s.params.items():
        env[k] = float(val) if np.isscalar(val) else np.asarray(val, float).reshape(np.shape(val))
        if np.ndim(val) == 1:
            env[k] = env[k].reshape(3, 1)
    mdl = Model(s, name)
    L_lib = Evaluator(env=dict(env), rng=rng)(s.L)
    c0, cp, cm = flow(keys, env, 0.0), flow(keys, env, H), flow(keys, env, -H)
    T = sum(0.5 * col(env[k]) @ J @ col(env[k]) for J, k in mdl.bodies)
    for (m, xp), (_, xm) in zip(mdl.points(cp, H), mdl.points(cm, -H)):
        v = (xp - xm) / (2 * H)
        T += 0.5 * m * v @ v
    T += mdl.extra_T(cp, cm)
    V = mdl.springs(c0) + sum(m * mdl.g * (x @ UP) for m, x in mdl.points(c0, 0.0))
    L_true = T - V

    u = {str(k): float(env[str(k)]) if isinstance(k, Scalar) else col(env[str(k)])
         for k in s.inputs}  # fmt: skip
    direction, venv = {}, dict(env)
    for var, kq, _ in keys:
        if isinstance(var, Scalar):
            d = float(rng.normal())
            venv[str(var.delta())] = d
        elif isinstance(var, S1):  # δq = γ e × q
            d = float(rng.normal())
            venv[str(var.get_variation_vector())] = d
        else:
            d = rng.normal(size=3)
            if isinstance(var, S2):
                d = d - (d @ col(env[kq])) * col(env[kq])
            key = var.get_variation_vector() if isinstance(var, (S2, SO3)) else var.delta()
            venv[str(key)] = d.reshape(3, 1)
        direction[kq] = d
    dW_lib = Evaluator(env=venv, rng=rng)(s.dW)
    vp, vm = flow(keys, env, H, direction), flow(keys, env, -H, direction)
    dW = sum(F @ ((pt(vp) - pt(vm)) / (2 * H)) for F, pt in mdl.forces(c0, u))
    dW += sum(M @ col(direction[k]) for M, k in mdl.moments(u))
    dW += sum(f * direction[k] for f, k in mdl.scalar_forces(u))
    dW += mdl.extra_W(c0, vp, vm, u)
    return (
        abs(L_lib - L_true) / max(1.0, abs(L_true)),
        abs(dW_lib - dW) / max(1.0, abs(dW)),
    )


def _cases():
    out = []
    for name in names():
        if name in SIZED:
            out += [(name, {SIZED[name]: n}) for n in (1, 2, 3)]
        else:
            out.append((name, {}))
    return out


CASES = _cases()


@pytest.mark.parametrize(
    "name, kwargs", CASES, ids=[n + "".join(f"-{v}" for v in kw.values()) for n, kw in CASES]
)
def test_lagrangian_and_virtual_work_match_a_model_of_the_bodies(name, kwargs):
    s = get(name, **kwargs)
    for seed in SEEDS:
        eL, eW = errors(s, name, seed)
        assert eL < TOL, f"{name} {kwargs}: L differs from the model by {eL:.1e} (seed {seed})"
        assert eW < TOL, f"{name} {kwargs}: δW differs from the model by {eW:.1e} (seed {seed})"


def test_every_system_has_a_model():
    assert {n for n, _ in CASES} == set(names())


def test_inclined_cart_rail():
    """The cart-pole with a rail inclined 30° (``e₁·e₃ ≠ 0``): the cart's weight enters."""
    import dataclasses
    from types import MappingProxyType

    s = get("cart_spherical_pendulum")
    p = dict(s.params)
    p["e1"] = np.array([np.cos(np.pi / 6), 0.0, np.sin(np.pi / 6)])
    tilted = dataclasses.replace(s, params=MappingProxyType(p))
    for seed in SEEDS:
        eL, eW = errors(tilted, "cart_spherical_pendulum", seed)
        assert eL < TOL and eW < TOL, (eL, eW)


def test_the_check_detects_a_missing_term():
    """Sanity check of the check: the heavy top with the free rigid body's Lagrangian (no
    gravity; same symbols J, R) fails, and so does a virtual work with the moment
    reversed."""
    import dataclasses

    top = get("heavy_top")
    no_gravity = dataclasses.replace(top, L=get("free_rigid_body").L)
    assert errors(no_gravity, "heavy_top", 0)[0] > 1e-3
    reversed_moment = dataclasses.replace(top, dW=-top.dW)
    assert errors(reversed_moment, "heavy_top", 0)[1] > 1e-3


def _replaced(s, **changes):
    import dataclasses

    return dataclasses.replace(s, **changes)


def test_the_check_detects_a_dropped_hip_mass():
    """The biped's Lagrangian without the hip point mass fails the model check."""
    from geomech import Dot, Scalar, Vector
    from geomech.core.base import HALF

    s = get("biped_3d")
    mH, g = Scalar("m_H", attr=["Constant"]), Scalar("g", attr=["Constant"])
    e3, xH = Vector("e3", attr=["Constant"]), s.symbols["xH"]
    vH = xH.t_diff()
    no_hip = s.L - HALF * mH * Dot(vH, vH) + mH * g * Dot(e3, xH)
    assert errors(_replaced(s, L=no_hip), "biped_3d", 0)[0] > 1e-3
    assert errors(s, "biped_3d", 0)[0] < TOL


def test_the_check_detects_a_flipped_wheel_torque():
    """The spacecraft's virtual work with the sign of τ₁ on its wheel flipped fails."""
    s = get("spacecraft_reaction_wheels")
    tau1, th1 = s.symbols["tau"][0], s.symbols["theta"][0]
    flipped = s.dW - tau1 * th1.delta() * 2
    assert errors(_replaced(s, dW=flipped), "spacecraft_reaction_wheels", 0)[1] > 1e-3
    # and the bus reaction is needed too: the wheel torque alone, without −τ b on the bus
    b1, eta = s.symbols["b"][0], s.symbols["R"].get_variation_vector()
    from geomech import Dot

    no_reaction = s.dW + tau1 * Dot(b1, eta)
    assert errors(_replaced(s, dW=no_reaction), "spacecraft_reaction_wheels", 0)[1] > 1e-3


def test_the_check_detects_a_missing_joint_reaction():
    """The rigid chain's virtual work without the reaction −τ₂ on link 1 fails, and so
    does its Lagrangian without the rotational energy of link 2."""
    from geomech import Dot, MVMul, Transpose
    from geomech.core.base import HALF

    s = get("n_link_rigid_chain", n_links=2)
    R1, R2 = s.symbols["R"]
    tau2 = s.symbols["tau"][1]
    no_reaction = s.dW + Dot(R1.get_variation_vector(), MVMul(Transpose(R1), tau2))
    assert errors(_replaced(s, dW=no_reaction), "n_link_rigid_chain", 0)[1] > 1e-3
    Om2 = R2.get_tangent_vector()
    from geomech import Matrix

    J2 = Matrix("J2", attr=["Constant", "SymmetricMatrix"])
    no_spin = s.L - Dot(Om2, J2 * Om2) * HALF
    assert errors(_replaced(s, L=no_spin), "n_link_rigid_chain", 0)[0] > 1e-3
    assert max(errors(s, "n_link_rigid_chain", 0)) < TOL


def test_the_check_detects_a_missing_panel_term():
    """The solar-array spacecraft without the translational energy of panel 1, or without
    the bus reaction of drive 1, fails the model check."""
    from geomech import Cross, Dot
    from geomech.core.base import HALF

    s = get("spacecraft_solar_arrays")
    assert max(errors(s, "spacecraft_solar_arrays", 0)) < TOL
    Om = s.symbols["R"].get_tangent_vector()
    a1, m1, d1 = s.symbols["a"][0], s.symbols["m_p"][0], s.symbols["d"][0]
    v = Cross(Om, a1)
    no_translation = s.L - HALF * m1 * d1 * d1 * Dot(v, v)
    assert errors(_replaced(s, L=no_translation), "spacecraft_solar_arrays", 0)[0] > 1e-3
    mu1, eta = s.symbols["mu"][0], s.symbols["R"].get_variation_vector()
    no_reaction = s.dW + mu1 * Dot(a1, eta)
    assert errors(_replaced(s, dW=no_reaction), "spacecraft_solar_arrays", 0)[1] > 1e-3


def test_the_check_detects_a_wrong_planar_cart_pole():
    """The planar cart-pole's Lagrangian without the cart–pendulum coupling ``m ṡ e₁·(l q̇)``
    fails the model check, and so does the rail force applied to the bob instead of the
    cart."""
    from geomech import Dot

    s = get("cart_planar_pendulum")
    assert max(errors(s, "cart_planar_pendulum", 0)) < TOL
    sym = s.symbols
    m, l, e1, q, f = sym["m"], sym["l"], sym["e1"], sym["q"], sym["f"]
    no_coupling = s.L - m * l * sym["s"].t_diff() * Dot(e1, q.t_diff())
    assert errors(_replaced(s, L=no_coupling), "cart_planar_pendulum", 0)[0] > 1e-3
    on_bob = s.dW + f * l * Dot(e1, q.delta())
    assert errors(_replaced(s, dW=on_bob), "cart_planar_pendulum", 0)[1] > 1e-3
