"""Physical checks of ``simulate``: conserved quantities, closed forms, order, drift.

1. The energy/momentum table of ``specs/simulation.md``: over 10 s from the states of
   ``tests/simulation_systems.py``, with ``rkmk4`` (dt = 1e-3) and with ``solve_ivp``
   (DOP853, rtol = 1e-10), the generated energy is conserved to the stated relative
   tolerance, and for the free rigid body also the inertial angular momentum ``R J Ω``.
   Inputs are zero (the quadrotors fall freely).
2. Closed forms integrated in NumPy without geomech: the spherical pendulum (Lee et al.
   Eq. 5.14 in R³) and the free rigid body (Euler's equations with ``Ṙ = R Ω̂``). Both
   integrators agree with them to 1e-8 over 2 s.
3. RKMK4 has order 4: halving dt divides the error at t = 2 s by about 16, on SO(3) (where
   the ``dexp⁻¹`` correction matters) and on S².
4. RKMK4 keeps ``|q| − 1`` and ``‖RᵀR − I‖`` at round-off (< 1e-12) over 10 s without
   projection.
"""

from __future__ import annotations

import functools

import numpy as np
import pytest

from geomech.numeric import hat
from geomech.simulation import simulate
from tests.simulation_systems import SETUPS, standard_module

T_END = 10.0
DT = 1e-3

#: (system, relative tolerance on E over 10 s), from the table of specs/simulation.md.
ENERGY_TABLE = [
    ("spherical_pendulum_s2", 1e-8),
    ("free_rigid_body", 1e-9),
    ("rigid_pendulum_pivot", 1e-8),
    ("double_rigid_pendulum", 1e-7),
    ("elastic_spherical_pendulum", 1e-7),
    ("quadrotor_with_load", 1e-7),
]
METHODS = ["rkmk4", "solve_ivp"]


@functools.cache
def run(name: str, method: str):
    p, x0 = SETUPS[name]
    return simulate(standard_module(name), x0, T_END, p, method=method, dt=DT, rtol=1e-10)


def relative_energy_drift(traj) -> float:
    E = traj.energy()
    return float(np.max(np.abs(E - E[0])) / abs(E[0]))


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("name, tol", ENERGY_TABLE)
def test_energy_conserved(name, tol, method):
    traj = run(name, method)
    assert traj.t[-1] == T_END
    E = traj.energy()
    assert abs(E[0]) > 0.3  # the relative drift is meaningful
    drift = relative_energy_drift(traj)
    assert drift < tol, f"{name} ({method}): relative energy drift {drift:.2e} >= {tol:.0e}"


@pytest.mark.parametrize("method", METHODS)
def test_free_rigid_body_inertial_momentum_conserved(method):
    traj = run("free_rigid_body", method)
    J = SETUPS["free_rigid_body"][0]["J"]
    pi = np.einsum("nij,jk,nk->ni", traj["R"], J, traj["Omega_R"])
    assert np.linalg.norm(pi[0]) > 1.0
    err = float(np.max(np.abs(pi - pi[0])))
    assert err < 1e-9, f"|R J Ω − π₀| = {err:.2e}"


def test_quadrotors_fall_freely_with_zero_inputs():
    """The table runs have u = 0: the centre of mass accelerates at −g e₃."""
    traj = run("quadrotor_with_load", "rkmk4")
    p = SETUPS["quadrotor_with_load"][0]
    xc = (p["m_L"] * traj["x_L"] + p["m_Q"] * (traj["x_L"] - p["l"] * traj["q"])) / (
        p["m_L"] + p["m_Q"]
    )
    t = traj.t
    v0 = (xc[1] - xc[0]) / (t[1] - t[0]) + 0.5 * p["g"] * (t[1] - t[0]) * np.array([0, 0, 1])
    ref = xc[0] + np.outer(t, v0) - 0.5 * p["g"] * np.outer(t**2, [0, 0, 1])
    np.testing.assert_allclose(xc, ref, atol=1e-6)


@pytest.mark.parametrize("name", [n for n, _ in ENERGY_TABLE])
def test_rkmk4_constraint_drift_at_round_off(name):
    traj = run(name, "rkmk4")
    for key, values in traj.constraint_drift.items():
        if "." in key:  # q·ω: a velocity constraint, kept by the dynamics to O(h⁴)
            assert values.max() < 1e-9, (key, values.max())
        else:
            assert values.max() < 1e-12, (key, values.max())


# ---------------------------------------------------------------------------
# Closed forms written in NumPy (no geomech)
# ---------------------------------------------------------------------------


def _reference(rhs, y0, t):
    integrate = pytest.importorskip("scipy.integrate")
    sol = integrate.solve_ivp(rhs, (t[0], t[-1]), y0, "DOP853", t_eval=t, rtol=1e-13, atol=1e-14)
    assert sol.success
    return sol.y.T


def _sim(name, method, t_end=2.0):
    p, x0 = SETUPS[name]
    return simulate(standard_module(name), x0, t_end, p, method=method, rtol=1e-12, atol=1e-13)


@pytest.mark.parametrize("method", METHODS)
def test_spherical_pendulum_matches_lee_eq_5_14(method):
    """Lee et al. Eq. 5.14, ``m l² ω̇ = m g l e₃ × q`` for the potential ``m g l e₃·q`` of
    this Lagrangian, with ``q̇ = ω × q``, integrated as an ODE in R³ × R³."""
    p, x0 = SETUPS["spherical_pendulum_s2"]
    g, l, e3 = p["g"], p["l"], p["e3"]

    def rhs(t, y):
        q, w = y[:3], y[3:]
        return np.r_[np.cross(w, q), (g / l) * np.cross(e3, q)]

    traj = _sim("spherical_pendulum_s2", method)
    ref = _reference(rhs, x0, traj.t)
    assert np.ptp(ref[:, 2]) > 0.3  # it swings
    err = float(np.max(np.abs(traj.state - ref)))
    assert err < 1e-8, f"max deviation {err:.2e}"


@pytest.mark.parametrize("method", METHODS)
def test_free_rigid_body_matches_euler_equations(method):
    """``J Ω̇ = J Ω × Ω``, ``Ṙ = R Ω̂``, integrated as an ODE in R⁹ × R³."""
    p, x0 = SETUPS["free_rigid_body"]
    J = p["J"]

    def rhs(t, y):
        R, W = y[:9].reshape(3, 3), y[9:]
        return np.r_[(R @ hat(W)).reshape(9), np.linalg.solve(J, np.cross(J @ W, W))]

    traj = _sim("free_rigid_body", method)
    ref = _reference(rhs, x0, traj.t)
    assert np.ptp(ref[:, 9:], axis=0).max() > 0.1  # Ω changes in the body frame
    err = float(np.max(np.abs(traj.state - ref)))
    assert err < 1e-8, f"max deviation {err:.2e}"


# ---------------------------------------------------------------------------
# Order of RKMK4
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["free_rigid_body", "spherical_pendulum_s2"])
def test_rkmk4_is_fourth_order(name):
    p, x0 = SETUPS[name]
    x0 = x0.copy()
    if name == "free_rigid_body":
        x0[9:] = [2.0, -3.0, 1.5]  # fast tumbling: the brackets of dexp⁻¹ matter
    mod = standard_module(name)
    ref = simulate(mod, x0, 2.0, p, method="solve_ivp", rtol=1e-13, atol=1e-13).state[-1]
    errs = [
        float(np.max(np.abs(simulate(mod, x0, 2.0, p, dt=dt).state[-1] - ref)))
        for dt in (0.04, 0.02, 0.01)
    ]
    ratios = [errs[0] / errs[1], errs[1] / errs[2]]
    assert all(13.0 < r < 19.0 for r in ratios), (errs, ratios)
