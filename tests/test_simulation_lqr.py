"""``geomech.simulation.lqr``: gains from the generated linearization, checked in closed loop.

1. Rigid pendulum about its pivot (SO(3)), inverted equilibrium: from 20° off about a
   random axis, at rest, ``u = −K δx`` brings ``‖η‖ = ‖vee(log(R_eqᵀ R))‖`` below 1e-3
   within 5 s (several axes).
2. Inverted spherical pendulum (q near +e₃, the maximum of the potential ``m g l e₃·q``)
   with the tangent torque input: it converges to the upright point.
3. The S² tangent-basis reduction: the 3-D linearization of an S² point has an
   uncontrollable mode along ``q`` (PBH rank test at λ = 0), the 2-D one does not; the
   basis is orthonormal and tangent; for the quadrotor with a load the full ``M_lin`` is
   singular (``Hat(q) Hat(q)``) and the reduced one is not, and the gain stabilises it.
4. The state error: ``η`` is exact for ``R = R_eq exp(η̂)``, ``ξ = q_eq × q``.
"""

from __future__ import annotations

import numpy as np
import pytest

from geomech.numeric import expm_so3
from geomech.simulation import lqr, simulate, tangent_basis
from tests.simulation_systems import E3, J_A, SETUPS, linearized_module, standard_module, unit

pytest.importorskip("scipy.linalg", reason="lqr needs scipy (geomech[sim])")

R_UP = np.diag([1.0, -1.0, -1.0])  # half turn about e1: R ρ points up for ρ = −0.5 e3
P_PENDULUM = {"J": J_A, "\\rho": np.array([0.0, 0.0, -0.5]), "e3": E3, "g": 9.81, "m": 1.0}
P_SPHERICAL = SETUPS["spherical_pendulum_s2"][0]


def pendulum_controller():
    x_eq = np.r_[R_UP.reshape(9), np.zeros(3)]
    Q = np.r_[[100.0] * 3, [10.0] * 3]
    return lqr(linearized_module("rigid_pendulum_pivot"), x_eq, np.zeros(3), P_PENDULUM, Q, 1.0)


def test_inverted_pendulum_open_loop_is_unstable():
    ctrl = pendulum_controller()
    assert np.max(np.linalg.eigvals(ctrl.A).real) > 1.0
    assert np.max(ctrl.closed_loop_eigenvalues.real) < 0
    assert ctrl.K.shape == (3, 6)
    assert ctrl.labels == [f"eta_R[{i}]" for i in range(3)] + [
        f"delta_Omega_R[{i}]" for i in range(3)
    ]


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_inverted_rigid_pendulum_stabilised(seed):
    ctrl = pendulum_controller()
    axis = unit(np.random.default_rng(seed).normal(size=3))
    x0 = np.r_[(R_UP @ expm_so3(np.deg2rad(20.0) * axis)).reshape(9), np.zeros(3)]
    mod = standard_module("rigid_pendulum_pivot")
    assert np.linalg.norm(ctrl.error(x0)[:3]) == pytest.approx(np.deg2rad(20.0))
    traj = simulate(mod, x0, 5.0, P_PENDULUM, u=ctrl)
    eta = np.array([np.linalg.norm(ctrl.error(x)[:3]) for x in traj.state])
    assert eta[-1] < 1e-3, eta[-1]
    assert np.linalg.norm(traj.state[-1][9:]) < 1e-2
    # without control the pendulum falls: the tilt of R ρ from e₃ more than doubles in 1 s
    fall = simulate(mod, x0, 1.0, P_PENDULUM)
    assert _tilt(fall.state[-1]) > 2.0 * _tilt(x0) > 0.1


def _tilt(x):
    up = x[:9].reshape(3, 3) @ P_PENDULUM["\\rho"]
    return float(np.arccos(np.clip(up @ E3 / np.linalg.norm(up), -1.0, 1.0)))


def spherical_controller():
    x_eq = np.r_[E3, np.zeros(3)]
    Q = np.r_[[10.0] * 2, [1.0] * 2]
    return lqr(linearized_module("spherical_pendulum_s2"), x_eq, np.zeros(3), P_SPHERICAL, Q, 1.0)


def test_inverted_spherical_pendulum_stabilised():
    ctrl = spherical_controller()
    assert ctrl.K.shape == (3, 4)
    assert ctrl.labels == ["xi_q[0]", "xi_q[1]", "delta_omega_q[0]", "delta_omega_q[1]"]
    # the input along q has no effect on a tangent equation: its row of K is zero
    np.testing.assert_allclose(E3 @ ctrl.K, 0.0, atol=1e-12)
    q0 = expm_so3(np.deg2rad(15.0) * unit([1.0, 2.0, 0.0])) @ E3
    x0 = np.r_[q0, np.zeros(3)]
    mod = standard_module("spherical_pendulum_s2")
    traj = simulate(mod, x0, 5.0, P_SPHERICAL, u=ctrl)
    assert np.linalg.norm(traj["q"][-1] - E3) < 1e-4
    assert np.linalg.norm(traj["omega_q"][-1]) < 1e-3
    assert traj.max_drift() < 1e-12
    fall = simulate(mod, x0, 2.0, P_SPHERICAL)
    assert np.linalg.norm(fall["q"][-1] - E3) > 0.2


def _uncontrollable_at_zero(A, B) -> bool:
    n = A.shape[0]
    return np.linalg.matrix_rank(np.hstack([A, B]), tol=1e-9) < n


def test_s2_reduction_removes_the_uncontrollable_normal_mode():
    lin = linearized_module("spherical_pendulum_s2")
    x_eq, a0, u0 = np.r_[E3, np.zeros(3)], np.zeros(3), np.zeros(3)
    M = lin.M_lin(x_eq, a0, u0, P_SPHERICAL)
    A3 = np.linalg.solve(M, lin.A_lin(x_eq, a0, u0, P_SPHERICAL))
    B3 = np.linalg.solve(M, lin.B_lin(x_eq, a0, u0, P_SPHERICAL))
    assert _uncontrollable_at_zero(A3, B3)  # ξ·q and δω·q: a zero mode no input reaches
    ctrl = spherical_controller()
    assert not _uncontrollable_at_zero(ctrl.A, ctrl.B)
    np.testing.assert_allclose(ctrl.T.T @ ctrl.T, np.eye(4), atol=1e-15)
    np.testing.assert_allclose(E3 @ ctrl.T[:3], 0.0, atol=1e-15)
    np.testing.assert_allclose(E3 @ ctrl.T[3:], 0.0, atol=1e-15)


@pytest.mark.parametrize("q", [E3, -E3, unit([1.0, 1.0, 1.0]), unit([0.3, -2.0, 0.1])])
def test_tangent_basis(q):
    B = tangent_basis(q)
    np.testing.assert_allclose(B.T @ B, np.eye(2), atol=1e-15)
    np.testing.assert_allclose(q @ B, 0.0, atol=1e-15)
    np.testing.assert_allclose(np.cross(B[:, 0], B[:, 1]), q, atol=1e-15)


def test_quadrotor_with_load_singular_s2_block_reduced_and_stabilised():
    p = SETUPS["quadrotor_with_load"][0]
    x_eq = np.r_[np.zeros(3), -E3, np.eye(3).reshape(9), np.zeros(9)]  # quadrotor above load
    u_eq = np.r_[(p["m_Q"] + p["m_L"]) * p["g"], 0.0, 0.0, 0.0]
    lin = linearized_module("quadrotor_with_load")
    M = lin.M_lin(x_eq, np.zeros(9), u_eq, p)
    assert np.linalg.matrix_rank(M) == 17  # Hat(q) Hat(q) m l² on the ω̇ columns
    ctrl = lqr(lin, x_eq, u_eq, p, Q=100.0, R=1.0)
    assert ctrl.K.shape == (4, 16)
    assert np.max(ctrl.closed_loop_eigenvalues.real) < 0
    x0 = x_eq.copy()
    x0[:3] = [0.3, -0.2, 0.1]
    x0[3:6] = expm_so3([0.2, 0.1, 0.0]) @ -E3
    x0[6:15] = expm_so3([0.1, -0.15, 0.2]).reshape(9)
    traj = simulate(standard_module("quadrotor_with_load"), x0, 5.0, p, u=ctrl)
    e0, e1 = np.abs(ctrl.error(x0)).max(), np.abs(ctrl.error(traj.state[-1])).max()
    assert e1 < 0.1 * e0, (e0, e1)  # slowest closed-loop mode ≈ −1/s


def test_state_error_coordinates():
    ctrl = pendulum_controller()
    eta = np.array([0.1, -0.2, 0.05])
    W = np.array([0.3, 0.0, -0.1])
    x = np.r_[(R_UP @ expm_so3(eta)).reshape(9), W]
    np.testing.assert_allclose(ctrl.error(x), np.r_[eta, W], atol=1e-14)
    np.testing.assert_allclose(ctrl(0.0, x), -ctrl.K @ np.r_[eta, W], atol=1e-12)
    sph = spherical_controller()
    q = unit([0.02, -0.01, 1.0])
    B = sph.bases["q"]
    np.testing.assert_allclose(sph.error(np.r_[q, 0, 0, 0])[:2], B.T @ np.cross(E3, q))


def test_lqr_argument_errors():
    lin = linearized_module("spherical_pendulum_s2")
    x_eq = np.r_[E3, np.zeros(3)]
    with pytest.raises(ValueError, match=r"Q must be a scalar, \(4,\) or \(4, 4\), got \(6, 6\)"):
        lqr(lin, x_eq, np.zeros(3), P_SPHERICAL, np.eye(6), 1.0)
    with pytest.raises(ValueError, match="x_eq: the S2 slot 'q' must have unit norm"):
        lqr(lin, 2 * x_eq, np.zeros(3), P_SPHERICAL, 1.0, 1.0)
    with pytest.raises(ValueError, match="u_eq must have 3 entries, got 1"):
        lqr(lin, x_eq, 0.0, P_SPHERICAL, 1.0, 1.0)
    with pytest.raises(TypeError, match="no M_lin\\(\\); generate it from linearize"):
        lqr(standard_module("spherical_pendulum_s2"), x_eq, np.zeros(3), P_SPHERICAL, 1.0, 1.0)
