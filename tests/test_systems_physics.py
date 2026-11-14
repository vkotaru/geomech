"""Every Lagrangian and virtual work of ``geomech.systems`` agrees with an independent model
of the system's bodies.

For each system (parametric ones at sizes 1, 2 and 3) the model below is written in NumPy
from the system's physical description, without the library's expressions:

- the positions of the point masses and centres of mass as functions of the
  configuration; the kinetic energy is ``Σ ½ m |ẋ|²`` with ``ẋ`` by central differences
  along the trajectory through the state (``q(t) = exp(t ω̂) q``, ``R(t) = R exp(t Ω̂)``,
  ``x(t) = x + t ẋ``), plus ``½ Ωᵀ J_c Ω`` per rigid body with ``J_c`` about its centre of
  mass (a pivot inertia is shifted to the centre of mass with the parallel-axis theorem);
- the potential: ``Σ m g (height)`` with gravity along −e₃, plus the springs;
- the virtual work: each force times the central difference of its application point
  along a random admissible variation (``δq = ξ × q``, ``δR = R η̂``), each body moment
  ``M·η`` and each moment about a pivot ``f·ξ`` (the work of a moment on the rotation ξ).

At random admissible states and the default parameters, ``L`` and ``δW`` evaluated with
``geomech.numeric.Evaluator`` agree with the model to 1e-8 relative to ``max(1, |model|)``.
"""

from __future__ import annotations

import numpy as np
import pytest

from geomech.core.base import S2, SO3, Scalar
from geomech.numeric import Evaluator, expm_so3, sample_state
from geomech.systems import get, names

H = 1e-6
TOL = 1e-8
SEEDS = (0, 1, 2)
UP = np.array([0.0, 0.0, 1.0])  # the physical vertical, independent of the param "e3"
SIZED = {
    "n_link_spherical_chain": "n_links",
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
        vel = v.get_tangent_vector() if isinstance(v, S2) else v.t_diff()
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
        self.g = self.p.get("g", 0.0)
        getattr(self, name)(n)

    def v(self, k):
        val = self.p[k]
        return float(val) if np.isscalar(val) else np.asarray(val, dtype=float)

    def body_moment(self, *pairs):
        self.moments = lambda u: [(u[m], k) for m, k in pairs]

    # ---- R, R³, S²
    def point_mass(self, n):
        self.points = lambda c, t: [(self.v("m"), c["x"])]
        self.forces = lambda c, u: [(u["F"], lambda c: c["x"])]

    def spring_mass(self, n):
        self.points = lambda c, t: [(self.v("m"), c["x"] * np.array([1.0, 0, 0]))]
        self.springs = lambda c: 0.5 * self.v("k") * c["x"] ** 2
        self.scalar_forces = lambda u: [(u["f"], "x")]

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
    env = sample_state(s.variables, list(s.inputs), rng)
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
    V = mdl.springs(c0) + sum(m * mdl.g * (x @ UP) for m, x in mdl.points(c0, 0.0))
    L_true = T - V

    u = {str(k): float(env[str(k)]) if isinstance(k, Scalar) else col(env[str(k)])
         for k in s.inputs}  # fmt: skip
    direction, venv = {}, dict(env)
    for var, kq, _ in keys:
        if isinstance(var, Scalar):
            d = float(rng.normal())
            venv[str(var.delta())] = d
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
