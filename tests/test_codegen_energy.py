"""``energy_function``: E = Σ (∂L/∂v)·v − L, derived with the engine, checked numerically.

1. ``momenta`` (∂L/∂v from the first variation of L) equals a central difference of the
   generated Lagrangian in each velocity slot, for every system (on S² compared on the
   tangent plane: the velocity ω is tangent, so only that part of ∂L/∂ω enters E).
2. The energy equals ``d/dλ L(q, λv)|₁ − L`` by a finite difference in λ.
3. Closed forms for a point mass, the spherical pendulum and the free rigid body.
4. ``dE/dt = 0`` along the unforced dynamics of the nine VBL systems (accelerations from
   the generated standard form, S² rows solved on the tangent plane).
5. Spec verification item 5: a free rigid body integrated for 10 s with
   ``scipy.integrate.solve_ivp`` from the generated ``mass``/``bias``/``kinematics``
   conserves the generated energy to 1e-8.
"""

from __future__ import annotations

import numpy as np
import pytest

from geomech import Dot, Mul
from geomech.codegen import energy_expr, energy_function, lambdify, momenta
from tests import vbl_systems
from tests.codegen_systems import (
    ALL_SYSTEMS,
    VBL_SYSTEMS,
    derived,
    modules,
    sample_points,
    tangent_accelerations,
)

_ids = dict(ids=lambda f: f.__name__)
H = 1e-5


def _lagrangian(d):
    return lambdify(d.system.L, d.layout)


@pytest.mark.parametrize("build", ALL_SYSTEMS, **_ids)
def test_momenta_match_finite_difference(build):
    d = derived(build)
    L = _lagrangian(d)
    mom = {
        k: lambdify(p, d.layout) for k, (_, p) in momenta(d.system.L, d.system.variables).items()
    }
    for pt in sample_points(d, n=5):
        for cfg, vel in zip(d.layout.configuration, d.layout.velocity):
            fd = np.zeros(vel.size)
            for i in range(vel.size):
                e = np.zeros(d.layout.n_state)
                e[vel.start + i] = H
                fd[i] = (L(pt.state + e, pt.u, pt.p) - L(pt.state - e, pt.u, pt.p)) / (2 * H)
            got = np.atleast_1d(mom[vel.key](pt.state, pt.u, pt.p))
            if cfg.manifold == "S2":
                q = pt.state[cfg.slice]
                P = np.eye(3) - np.outer(q, q)
                got, fd = P @ got, P @ fd
            scale = max(1.0, float(np.max(np.abs(fd))))
            assert np.max(np.abs(got - fd)) < 1e-7 * scale, f"dL/d{vel.key}: {got} vs {fd}"


@pytest.mark.parametrize("build", ALL_SYSTEMS, **_ids)
def test_energy_is_velocity_scaling_derivative(build):
    d = derived(build)
    L, E = _lagrangian(d), energy_function(d.system.L, d.layout)
    vel = slice(d.layout.n_config, d.layout.n_state)
    for pt in sample_points(d, n=5):

        def scaled(lam):
            st = pt.state.copy()
            st[vel] *= lam
            return L(st, pt.u, pt.p)

        ref = (scaled(1 + H) - scaled(1 - H)) / (2 * H) - scaled(1.0)
        assert abs(E(pt.state, pt.p) - ref) < 1e-7 * max(1.0, abs(ref))


def test_energy_closed_forms():
    pm = derived(vbl_systems.point_mass_r3)
    x = pm.system.symbols["x"]
    m = pm.system.symbols["m"]
    v = x.t_diff()
    assert str(energy_expr(pm.system.L, pm.system.variables)) == str(Mul(m, Dot(v, v)) * 0.5)

    for build, closed in (
        (vbl_systems.spherical_pendulum_s2, _spherical_pendulum_energy),
        (vbl_systems.free_rigid_body, _rigid_body_energy),
    ):
        d = derived(build)
        E = energy_function(d.system.L, d.layout)
        for pt in sample_points(d, n=5):
            ref = closed(d.layout, pt)
            assert abs(E(pt.state, pt.p) - ref) < 1e-12 * max(1.0, abs(ref))


def _spherical_pendulum_energy(layout, pt):
    q, w = pt.state[0:3], pt.state[3:6]
    m, g, l, e3 = pt.p["m"], pt.p["g"], pt.p["l"], pt.p["e3"]
    return 0.5 * m * l * l * w @ w + m * g * l * q @ e3


def _rigid_body_energy(layout, pt):
    W = pt.state[9:12]
    return 0.5 * W @ pt.p["J"] @ W


def _physical(d, pt):
    """Parameters with positive masses and positive-definite symmetric inertias."""
    p = dict(pt.p)
    for key, val in p.items():
        if np.shape(val) == (3, 3):
            p[key] = np.diag([1.0, 2.0, 3.0]) + 0.1 * (val + val.T)
        elif np.shape(val) == ():
            p[key] = 1.0 + abs(val)
    return p


@pytest.mark.parametrize("build", VBL_SYSTEMS, **_ids)
def test_energy_conserved_along_unforced_dynamics(build):
    d = derived(build)
    sf_mod = modules(build)[0]
    for pt in sample_points(d, n=5):
        p, u = _physical(d, pt), np.zeros(d.layout.n_inputs)
        a = tangent_accelerations(sf_mod, d.layout, pt.state, u, p)
        rate = np.concatenate([sf_mod.kinematics(pt.state), a])
        E = sf_mod.energy
        dE = (E(pt.state + H * rate, p) - E(pt.state - H * rate, p)) / (2 * H)
        scale = max(1.0, float(np.max(np.abs(rate))) * abs(E(pt.state, p)))
        assert abs(dE) < 1e-7 * scale, f"dE/dt = {dE}"


def test_free_rigid_body_energy_smoke():
    """10 s of a free rigid body with solve_ivp: |E(t) − E(0)| / E(0) < 1e-8."""
    integrate = pytest.importorskip(
        "scipy.integrate", reason="the simulation smoke test needs scipy (geomech[sim])"
    )
    mod = modules(vbl_systems.free_rigid_body)[0]
    p = {"J": np.array([[1.0, 0.1, 0.0], [0.1, 2.0, 0.2], [0.0, 0.2, 3.0]])}
    state0 = np.concatenate([np.eye(3).reshape(9), [0.3, -1.2, 0.8]])

    def rhs(t, x):
        a = np.linalg.solve(mod.mass(x, p), -mod.bias(x, p))
        return np.concatenate([mod.kinematics(x), a])

    sol = integrate.solve_ivp(rhs, (0.0, 10.0), state0, method="DOP853", rtol=1e-12, atol=1e-12)
    assert sol.success
    E = np.array([mod.energy(x, p) for x in sol.y.T])
    assert np.max(np.abs(E - E[0])) / E[0] < 1e-8
