"""Adversarial review of ``geomech.simulation`` (``specs/reviews/simulation-review.md``).

Every reference here is written directly in NumPy, never through geomech:

- S² accelerations from a constrained Lagrangian in R³ with one multiplier per sphere
  (``M q̈ = F + Σ λᵢ ∂cᵢ``, ``qᵢ · q̈ᵢ = −|q̇ᵢ|²``), mapped to ``ω̇ = q × q̈``;
- rigid bodies from Euler's equations with ``Ṙ = R Ω̂`` in R¹²;
- LQR gains from a central-difference linearization of those reference dynamics in the
  coordinates ``R = R_eq exp(η̂)``, ``q = exp((Bξ)^) q_eq``, ``ω = B ν``.

Confirmed findings are strict xfails whose reason starts with ``F<n>:``; the other tests
record what held.
"""

from __future__ import annotations

import numpy as np
import pytest

from geomech import S2, Dot, Scalar, SystemVariables, Vector, getScalars
from geomech.simulation import accelerations, lqr, simulate
from geomech.simulation import integrators as _integrators
from tests.codegen_systems import modules
from tests.simulation_systems import (
    J_A,
    Q0,
    SETUPS,
    W0,
    linearized_module,
    standard_module,
    tangent,
    unit,
)
from tests.vbl_systems import HALF, System

scipy_linalg = pytest.importorskip("scipy.linalg", reason="needs scipy (geomech[sim])")
scipy_integrate = pytest.importorskip("scipy.integrate", reason="needs scipy (geomech[sim])")

E3 = np.array([0.0, 0.0, 1.0])


# ---------------------------------------------------------------------------
# Systems built here (not in the shared helpers)
# ---------------------------------------------------------------------------


def double_spherical_pendulum():
    """Two point masses on massless links, x₁ = l₁ q₁, x₂ = l₁ q₁ + l₂ q₂ (Lee et al. §5.4),
    potential with e₃ up, δW = ξ₁·f₁ + ξ₂·f₂."""
    m1, m2, l1, l2, g = getScalars("m1 m2 l1 l2 g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    q1, q2 = S2("q1"), S2("q2")
    xi1, xi2 = q1.get_variation_vector(), q2.get_variation_vector()
    f1, f2 = Vector("f1"), Vector("f2")
    v1 = q1.t_diff() * l1
    v2 = q1.t_diff() * l1 + q2.t_diff() * l2
    L = (
        HALF * m1 * Dot(v1, v1)
        + HALF * m2 * Dot(v2, v2)
        - m1 * g * l1 * Dot(q1, e3)
        - m2 * g * Dot(q1 * l1 + q2 * l2, e3)
    )
    dW = Dot(xi1, f1) + Dot(xi2, f2)
    return System(L, dW, SystemVariables(vectors=[q1, q2]), [f1, f2], {})


def gyro_spherical_pendulum():
    """Spherical pendulum with a gyroscopic term b e₃·ω (ω = q × q̇, a magnetic-type force)."""
    m, g, l, b = getScalars("m g l b", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    q, f = S2("q"), Vector("f")
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * m * l * l * Dot(om, om) + b * Dot(e3, om) - m * g * l * Dot(q, e3)
    return System(L, Dot(xi, f), SystemVariables(vectors=[q]), [f], {})


def free_sphere_no_inputs():
    """Free particle on a sphere with no inputs (δW = 0)."""
    m, l = getScalars("m l", attr=["Constant"])
    q = S2("q")
    om = q.get_tangent_vector()
    zero = Scalar("0", value=0.0, attr=["Constant"])
    return System(HALF * m * l * l * Dot(om, om), zero, SystemVariables(vectors=[q]), [], {})


P_DSP = {"m1": 1.3, "m2": 0.7, "l1": 1.1, "l2": 0.6, "g": 9.81, "e3": E3}
P_GYRO = {"m": 0.9, "g": 9.81, "l": 1.2, "b": 1.5, "e3": E3}


# ---------------------------------------------------------------------------
# Independent NumPy references
# ---------------------------------------------------------------------------


def hat(v):
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


def expm(v):
    th = np.linalg.norm(v)
    K = hat(v)
    if th < 1e-12:
        return np.eye(3) + K
    return np.eye(3) + np.sin(th) / th * K + (1 - np.cos(th)) / th**2 * (K @ K)


def kkt(M, F, C, rhs):
    """Solve ``M q̈ = F + Cᵀ λ``, ``C q̈ = rhs`` for ``q̈``."""
    n, k = M.shape[0], C.shape[0]
    K = np.zeros((n + k, n + k))
    K[:n, :n], K[:n, n:], K[n:, :n] = M, -C.T, C
    return np.linalg.solve(K, np.r_[F, rhs])[:n]


def dsp_qdd(q1, q2, qd1, qd2, f1, f2, p):
    m1, m2, l1, l2, g = p["m1"], p["m2"], p["l1"], p["l2"], p["g"]
    eye = np.eye(3)
    M = np.block(
        [[(m1 + m2) * l1**2 * eye, m2 * l1 * l2 * eye], [m2 * l1 * l2 * eye, m2 * l2**2 * eye]]
    )
    # generalized force of δW = ξ·f with δq = ξ × q is f × q
    F = np.r_[-(m1 + m2) * g * l1 * E3 + np.cross(f1, q1), -m2 * g * l2 * E3 + np.cross(f2, q2)]
    C = np.zeros((2, 6))
    C[0, :3], C[1, 3:] = q1, q2
    return kkt(M, F, C, np.array([-qd1 @ qd1, -qd2 @ qd2]))


def dsp_accel(x, u, p):
    q1, q2, w1, w2 = x[0:3], x[3:6], x[6:9], x[9:12]
    qdd = dsp_qdd(q1, q2, np.cross(w1, q1), np.cross(w2, q2), u[:3], u[3:], p)
    return np.r_[np.cross(q1, qdd[:3]), np.cross(q2, qdd[3:])]  # ω̇ = q × q̈


def dsp_accel_scaled(x, p):
    """The double spherical pendulum, zero inputs, solved with each sphere's equations
    divided by their own inertia (well conditioned for any mass ratio)."""
    m1, m2, l1, l2, g = p["m1"], p["m2"], p["l1"], p["l2"], p["g"]
    q1, q2, w1, w2 = x[0:3], x[3:6], x[6:9], x[9:12]
    qd1, qd2 = np.cross(w1, q1), np.cross(w2, q2)
    H1, H2 = -hat(q1), -hat(q2)  # q̈ = H ω̇ + ω × q̇
    c1, c2 = np.cross(w1, qd1), np.cross(w2, qd2)
    P1, P2 = np.eye(3) - np.outer(q1, q1), np.eye(3) - np.outer(q2, q2)
    s1, s2 = (m1 + m2) * l1**2, m2 * l2**2
    A, b = np.zeros((6, 6)), np.zeros(6)
    A[:3, :3] = P1 @ (s1 * H1) / s1 + np.outer(q1, q1)
    A[:3, 3:] = P1 @ (m2 * l1 * l2 * H2) / s1
    b[:3] = -P1 @ (s1 * c1 + m2 * l1 * l2 * c2 + (m1 + m2) * g * l1 * E3) / s1
    A[3:, :3] = P2 @ (m2 * l1 * l2 * H1) / s2
    A[3:, 3:] = P2 @ (s2 * H2) / s2 + np.outer(q2, q2)
    b[3:] = -P2 @ (m2 * l1 * l2 * c1 + s2 * c2 + m2 * g * l2 * E3) / s2
    return np.linalg.solve(A, b)


def gyro_qdd(q, qd, f, p):
    m, g, l, b = p["m"], p["g"], p["l"], p["b"]
    # d/dt ∂L/∂q̇ − ∂L/∂q for b (e₃ × q)·q̇ gives 2 b e₃ × q̇
    F = -2 * b * np.cross(E3, qd) - m * g * l * E3 + np.cross(f, q)
    return kkt(m * l * l * np.eye(3), F, q[None, :], np.array([-qd @ qd]))


def qwl_accel(x, u, p):
    """Quadrotor with load in (x_L, q) ∈ R⁶ with |q| = 1, plus Euler's equations."""
    mQ, mL, l, g, J = p["m_Q"], p["m_L"], p["l"], p["g"], p["J"]
    q, R = x[3:6], x[6:15].reshape(3, 3)
    w, Om = x[18:21], x[21:24]
    f, Mt = u[0], u[1:4]
    eye = np.eye(3)
    M = np.block([[(mQ + mL) * eye, -mQ * l * eye], [-mQ * l * eye, mQ * l * l * eye]])
    thrust = f * R @ E3
    F = np.r_[-(mQ + mL) * g * E3 + thrust, mQ * g * l * E3 - l * thrust]
    qd = np.cross(w, q)
    C = np.zeros((1, 6))
    C[0, 3:] = q
    sol = kkt(M, F, C, np.array([-qd @ qd]))
    Omd = np.linalg.solve(J, np.cross(J @ Om, Om) + Mt)
    return np.r_[sol[:3], np.cross(q, sol[3:]), Omd]


def rigid_body_rhs(t, y, J, torque):
    R, Om = y[:9].reshape(3, 3), y[9:]
    return np.r_[(R @ hat(Om)).reshape(9), np.linalg.solve(J, np.cross(J @ Om, Om) + torque(t))]


def pendulum_accel(R, Om, M, p):
    """Rigid pendulum about its pivot, L = ½ Ω·JΩ − m g e₃·(R ρ): body torque M."""
    J, rho, m, g = p["J"], p["\\rho"], p["m"], p["g"]
    tau = m * g * np.cross(R.T @ E3, rho)
    return np.linalg.solve(J, np.cross(J @ Om, Om) + tau + M)


def fd_jacobians(f, n, u0, eps=1e-6):
    A = np.column_stack([(f(e, u0) - f(-e, u0)) / (2 * eps) for e in eps * np.eye(n)])
    B = np.column_stack(
        [
            (f(np.zeros(n), u0 + e) - f(np.zeros(n), u0 - e)) / (2 * eps)
            for e in eps * np.eye(u0.size)
        ]
    )
    return A, B


def ref_lqr(A, B, Q, R):
    return np.linalg.solve(R, B.T @ scipy_linalg.solve_continuous_are(A, B, Q, R))


def rel(a, b):
    return float(np.max(np.abs(a - b)) / max(1.0, np.max(np.abs(b))))


def random_tangent_state(rng, n_spheres):
    qs = [unit(rng.normal(size=3)) for _ in range(n_spheres)]
    ws = [tangent(q, rng.normal(size=3)) for q in qs]
    return np.concatenate(qs + ws)


# ---------------------------------------------------------------------------
# (2) the S² projected solve against a constrained-Lagrangian solve with multipliers
# ---------------------------------------------------------------------------


def test_s2_solve_double_spherical_pendulum_matches_multiplier_solve():
    """Both S² rows couple to the other sphere's ω̇; inputs on both."""
    mod = modules(double_spherical_pendulum)[0]
    rng = np.random.default_rng(1)
    for _ in range(10):
        x, u = random_tangent_state(rng, 2), rng.normal(size=6)
        assert rel(accelerations(mod, x, u, P_DSP), dsp_accel(x, u, P_DSP)) < 1e-13


def test_s2_solve_quadrotor_with_load_matches_multiplier_solve():
    """The S² row couples to ẍ_L and its mass block is rank 2 (Hat(q)Hat(q) m_Q l²)."""
    mod = standard_module("quadrotor_with_load")
    p = SETUPS["quadrotor_with_load"][0]
    rng = np.random.default_rng(2)
    for _ in range(10):
        q = unit(rng.normal(size=3))
        x = np.r_[
            rng.normal(size=3),
            q,
            expm(rng.normal(size=3)).reshape(9),
            rng.normal(size=3),
            tangent(q, rng.normal(size=3)),
            rng.normal(size=3),
        ]
        u = np.r_[20.0, 0.0, 0.0, 0.0] + rng.normal(size=4)
        assert rel(accelerations(mod, x, u, p), qwl_accel(x, u, p)) < 1e-13


def test_s2_solve_gyroscopic_term_matches_multiplier_solve():
    mod = modules(gyro_spherical_pendulum)[0]
    rng = np.random.default_rng(3)
    for _ in range(10):
        x, f = random_tangent_state(rng, 1), rng.normal(size=3)
        q, w = x[:3], x[3:]
        ref = np.cross(q, gyro_qdd(q, np.cross(w, q), f, P_GYRO))
        assert rel(accelerations(mod, x, f, P_GYRO), ref) < 1e-13


def test_s2_solve_mass_ratio_1e6_no_false_alarm():
    """A light second link (m₂/m₁ = 1e-6) is accepted and solved to 1e-10."""
    mod = modules(double_spherical_pendulum)[0]
    p = {**P_DSP, "m1": 1.0, "m2": 1e-6, "l1": 1.0, "l2": 1.0}
    x = random_tangent_state(np.random.default_rng(4), 2)
    assert rel(accelerations(mod, x, None, p), dsp_accel_scaled(x, p)) < 1e-10


def test_s2_solve_light_link_accurate():
    mod = modules(double_spherical_pendulum)[0]
    p = {**P_DSP, "m1": 1.0, "m2": 1e-9, "l1": 1.0, "l2": 1.0}
    q2 = unit([0.3, -0.5, -0.8])
    x = np.r_[Q0, q2, W0, tangent(q2, [1.0, 0.2, 0.3])]
    assert rel(accelerations(mod, x, None, p), dsp_accel_scaled(x, p)) < 1e-12


def test_s2_solve_tether_not_rejected():
    mod = standard_module("spherical_pendulum_s2")
    p = {"m": 1e3, "l": 1e5, "g": 9.81, "e3": E3}
    x = np.r_[Q0, 1e-3 * W0]
    a = accelerations(mod, x, None, p)
    ref = np.cross(Q0, gyro_qdd(Q0, np.cross(1e-3 * W0, Q0), np.zeros(3), {**p, "b": 0.0}))
    assert rel(a, ref) < 1e-12


# ---------------------------------------------------------------------------
# (1) RKMK4: order 4 with time-varying inputs, exact flows, the right-action sign
# ---------------------------------------------------------------------------


def _ref_ivp(rhs, y0, T):
    sol = scipy_integrate.solve_ivp(rhs, (0, T), y0, method="DOP853", rtol=1e-13, atol=1e-14)
    return sol.y[:, -1]


def _order_case(name):
    """(module, x0, p, u(t), reference final state in the geomech layout, T, dts)."""
    T = 2.0
    if name == "so3_time_varying_torque":
        mod, p = standard_module("free_rigid_body"), {"J": J_A}
        x0 = SETUPS["free_rigid_body"][1].copy()
        x0[9:] = [2.0, -3.0, 1.5]

        def torque(t):
            return np.array([np.sin(3 * t), np.cos(2 * t), 0.5 * t])

        ref = _ref_ivp(lambda t, y: rigid_body_rhs(t, y, J_A, torque), x0, T)
        return mod, x0, p, torque, ref, T, [0.02, 0.01, 0.005]
    if name == "s2_gyroscopic":
        mod, p = modules(gyro_spherical_pendulum)[0], P_GYRO
        x0 = np.r_[Q0, W0]

        def f(t):
            return np.array([np.cos(t), 0.5 * np.sin(2 * t), 0.2])

        def rhs(t, y):
            return np.r_[y[3:], gyro_qdd(y[:3], y[3:], f(t), p)]

        y = _ref_ivp(rhs, np.r_[Q0, np.cross(W0, Q0)], T)
        return mod, x0, p, f, np.r_[y[:3], np.cross(y[:3], y[3:])], T, [0.02, 0.01, 0.005]
    if name == "s2xs2_double_spherical":
        mod, p = modules(double_spherical_pendulum)[0], P_DSP
        q2 = unit([0.3, -0.5, -0.8])
        w2 = tangent(q2, [1.0, 0.2, 0.3])
        x0 = np.r_[Q0, q2, W0, w2]

        def f(t):
            return np.r_[0.3 * np.sin(t), 0.0, 0.1, 0.0, 0.2 * np.cos(t), 0.0]

        def rhs(t, y):
            u = f(t)
            return np.r_[y[6:], dsp_qdd(y[0:3], y[3:6], y[6:9], y[9:12], u[:3], u[3:], p)]

        y = _ref_ivp(rhs, np.r_[Q0, q2, np.cross(W0, Q0), np.cross(w2, q2)], T)
        ref = np.r_[y[:6], np.cross(y[0:3], y[6:9]), np.cross(y[3:6], y[9:12])]
        return mod, x0, p, f, ref, T, [0.01, 0.005, 0.0025]
    assert name == "r3xs2xso3_quadrotor_with_load"
    mod = standard_module("quadrotor_with_load")
    p, x0 = SETUPS["quadrotor_with_load"]

    def f(t):
        return np.r_[20 + 3 * np.sin(2 * t), 0.1 * np.cos(t), -0.2, 0.05 * t]

    def rhs(t, y):
        q, qd, R, Om = y[3:6], y[18:21], y[6:15].reshape(3, 3), y[21:24]
        w = np.cross(q, qd)
        a = qwl_accel(np.r_[y[:18], w, Om], f(t), p)
        qdd = np.cross(a[3:6], q) + np.cross(w, qd)
        return np.r_[y[15:18], qd, (R @ hat(Om)).reshape(9), a[:3], qdd, a[6:]]

    y0 = x0.copy()
    y0[18:21] = np.cross(x0[18:21], x0[3:6])
    y = _ref_ivp(rhs, y0, T)
    y[18:21] = np.cross(y[3:6], y[18:21])
    return mod, x0, p, f, y, T, [0.02, 0.01, 0.005]


ORDER_CASES = [
    "so3_time_varying_torque",
    "s2_gyroscopic",
    "s2xs2_double_spherical",
    "r3xs2xso3_quadrotor_with_load",
]


def _errors(case):
    mod, x0, p, f, ref, T, dts = case
    return [
        float(np.max(np.abs(simulate(mod, x0, T, p, u=lambda t, x: f(t), dt=h).state[-1] - ref)))
        for h in dts
    ]


@pytest.mark.parametrize("name", ORDER_CASES)
def test_rkmk4_fourth_order_with_time_varying_input(name):
    errs = _errors(_order_case(name))
    ratios = [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]
    assert errs[-1] < 1e-4, errs
    assert 13.0 < ratios[-1] < 19.0, ratios


def test_rkmk4_right_action_sign_is_load_bearing(monkeypatch):
    """Flipping the odd dexp⁻¹ term for the right action on SO(3) drops the order to 2,
    so the sign in ``_LieRK.dexpinv`` is the right one and the order test can see it.
    (On S² the odd term is parallel to q at leading order and its sign is invisible.)"""
    cross = _integrators._cross

    def flipped(U, G, h):
        e, s, r = U
        ge, gs, gr = G
        if s.size:
            us = cross(s, gs)
            gs = gs - 0.5 * us + cross(s, us) / 12.0
        if r.size:
            ur = cross(r, gr)
            gr = gr - 0.5 * ur + cross(r, ur) / 12.0
        return h * ge, h * gs, h * gr

    case = _order_case("so3_time_varying_torque")
    monkeypatch.setattr(_integrators._LieRK, "dexpinv", staticmethod(flipped))
    errs = _errors(case)
    assert 3.5 < errs[-2] / errs[-1] < 4.5, errs


def test_rkmk4_s2_is_left_rotation_about_omega_exactly():
    """g = 0: ω is constant and q(t) = exp(t ω̂) q₀ (q̇ = ω × q); RKMK4 is exact."""
    mod = standard_module("spherical_pendulum_s2")
    p = {"m": 1.0, "l": 1.0, "g": 0.0, "e3": E3}
    w0 = tangent(Q0, [0.3, 2.5, -0.2])
    traj = simulate(mod, np.r_[Q0, w0], 3.0, p, dt=0.05)
    for t, x in zip(traj.t[::10], traj.state[::10]):
        np.testing.assert_allclose(x[:3], expm(t * w0) @ Q0, atol=1e-13)
        np.testing.assert_allclose(x[3:], w0, atol=1e-13)


def test_rkmk4_so3_is_right_multiplication_exactly():
    """J = I, no torque: Ω is constant and R(t) = R₀ exp(t Ω̂) (Ṙ = R Ω̂); RKMK4 is exact."""
    mod = standard_module("free_rigid_body")
    R0, Om = expm([0.1, 0.2, 0.3]), np.array([0.7, -2.0, 1.1])
    traj = simulate(mod, np.r_[R0.reshape(9), Om], 3.0, {"J": np.eye(3)}, dt=0.05)
    for t, x in zip(traj.t[::10], traj.state[::10]):
        np.testing.assert_allclose(x[:9].reshape(3, 3), R0 @ expm(t * Om), atol=1e-13)


# ---------------------------------------------------------------------------
# (3) solve_ivp output projection and restarting from a trajectory
# ---------------------------------------------------------------------------


def test_solve_ivp_final_state_is_a_valid_x0():
    p, x0 = SETUPS["spherical_pendulum_s2"]
    mod = standard_module("spherical_pendulum_s2")
    traj = simulate(mod, x0, 5.0, p, method="solve_ivp", rtol=1e-6, atol=1e-8, dt=1e-2)
    simulate(mod, traj.state[-1], 0.1, p, method="solve_ivp", rtol=1e-6)


def test_rkmk4_final_state_is_a_valid_x0():
    p, x0 = SETUPS["spherical_pendulum_s2"]
    mod = standard_module("spherical_pendulum_s2")
    traj = simulate(mod, x0, 20.0, p, dt=0.02)
    simulate(mod, traj.state[-1], 0.1, p)


def test_solve_ivp_drift_is_reported_before_projection():
    p, x0 = SETUPS["rigid_pendulum_pivot"]
    mod = standard_module("rigid_pendulum_pivot")
    traj = simulate(mod, x0, 2.0, p, method="solve_ivp", rtol=1e-6, atol=1e-8, dt=1e-2)
    assert traj.max_drift() > 1e-8  # the raw states left SO(3)
    R = traj["R"]
    assert np.max(np.abs(np.einsum("nki,nkj->nij", R, R) - np.eye(3))) < 1e-13


# ---------------------------------------------------------------------------
# (4) LQR against an independent finite-difference linearization
# ---------------------------------------------------------------------------

R_UP = np.diag([1.0, -1.0, -1.0])
P_PEND = {"J": J_A, "\\rho": np.array([0.0, 0.0, -0.5]), "e3": E3, "g": 9.81, "m": 1.0}


def test_lqr_inverted_rigid_pendulum_matches_independent_gain():
    def f(z, M):  # z = (η, Ω); η̇ = Ω to first order about Ω = 0
        R = R_UP @ expm(z[:3])
        return np.r_[z[3:], pendulum_accel(R, z[3:], M, P_PEND)]

    A, B = fd_jacobians(f, 6, np.zeros(3))
    Q = np.r_[[100.0] * 3, [10.0] * 3]
    Kref = ref_lqr(A, B, np.diag(Q), np.eye(3))
    x_eq = np.r_[R_UP.reshape(9), np.zeros(3)]
    ctrl = lqr(
        linearized_module("rigid_pendulum_pivot"),
        x_eq,
        np.zeros(3),
        P_PEND,
        Q,
        1.0,
        module=standard_module("rigid_pendulum_pivot"),
    )
    assert rel(ctrl.A, A) < 1e-8
    assert rel(ctrl.B, B) < 1e-8
    assert rel(ctrl.K, Kref) < 1e-8


def _sphere_fd(q_eq, basis):
    p = {**SETUPS["spherical_pendulum_s2"][0], "b": 0.0}

    def f(z, u):
        q, w = expm(basis @ z[:2]) @ q_eq, basis @ z[2:]
        qd = np.cross(w, q)
        wd = np.cross(q, gyro_qdd(q, qd, u, p))
        return np.r_[basis.T @ np.cross(q_eq, qd), basis.T @ wd]

    return fd_jacobians(f, 4, np.zeros(3))


@pytest.mark.parametrize("sign", [1.0, -1.0])
def test_lqr_spherical_pendulum_matches_independent_gain_at_poles(sign):
    """At q = ±e₃ (where a fixed basis could degenerate) the gain matches one built in
    the explicit basis (e₁, ±e₂), which has b₁ × b₂ = q."""
    q_eq = sign * E3
    basis = np.c_[[1.0, 0.0, 0.0], [0.0, sign, 0.0]]
    A, B = _sphere_fd(q_eq, basis)
    Q = np.r_[10.0, 10.0, 1.0, 1.0]
    Kref = ref_lqr(A, B, np.diag(Q), np.eye(3))
    p = SETUPS["spherical_pendulum_s2"][0]
    ctrl = lqr(
        linearized_module("spherical_pendulum_s2"),
        np.r_[q_eq, 0, 0, 0],
        np.zeros(3),
        p,
        Q,
        1.0,
        module=standard_module("spherical_pendulum_s2"),
    )
    np.testing.assert_allclose(ctrl.bases["q"], basis, atol=1e-15)
    assert rel(ctrl.K, Kref) < 1e-8


def test_lqr_quadrotor_with_load_matches_independent_gain():
    """Coupled S² row (ẍ_L in the q equation) and a rank-2 S² mass block, 16 coordinates."""
    p = SETUPS["quadrotor_with_load"][0]
    q_eq, basis = -E3, np.c_[[1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]
    u_eq = np.r_[(p["m_Q"] + p["m_L"]) * p["g"], 0.0, 0.0, 0.0]

    def f(z, u):
        q, R, w = expm(basis @ z[3:5]) @ q_eq, expm(z[5:8]), basis @ z[11:13]
        x = np.r_[z[0:3], q, R.reshape(9), z[8:11], w, z[13:16]]
        a = qwl_accel(x, u, p)
        qd = np.cross(w, q)
        return np.r_[
            z[8:11], basis.T @ np.cross(q_eq, qd), z[13:16], a[:3], basis.T @ a[3:6], a[6:]
        ]

    A, B = fd_jacobians(f, 16, u_eq)
    Kref = ref_lqr(A, B, 100.0 * np.eye(16), np.eye(4))
    x_eq = np.r_[np.zeros(3), q_eq, np.eye(3).reshape(9), np.zeros(9)]
    ctrl = lqr(
        linearized_module("quadrotor_with_load"),
        x_eq,
        u_eq,
        p,
        Q=100.0,
        R=1.0,
        module=standard_module("quadrotor_with_load"),
    )
    np.testing.assert_allclose(ctrl.bases["q"], basis, atol=1e-15)
    assert rel(ctrl.A, A) < 1e-8
    assert rel(ctrl.K, Kref) < 1e-7


@pytest.mark.filterwarnings("ignore:lqr. .x_eq, u_eq. is not checked:UserWarning")
def test_lqr_feedback_independent_of_tangent_basis_for_isotropic_weights():
    """Two equilibria 1e-9 apart where ``tangent_basis`` picks different axes (b₁ jumps by
    90°): with Q isotropic on each S² block the feedback u(state) is the same."""
    p = SETUPS["spherical_pendulum_s2"][0]
    lin = linearized_module("spherical_pendulum_s2")
    qa, qb = unit([1.0, 1.0 + 1e-9, 5.0]), unit([1.0 + 1e-9, 1.0, 5.0])
    ca = lqr(lin, np.r_[qa, 0, 0, 0], np.zeros(3), p, [10.0, 10.0, 1.0, 1.0], 1.0)
    cb = lqr(lin, np.r_[qb, 0, 0, 0], np.zeros(3), p, [10.0, 10.0, 1.0, 1.0], 1.0)
    assert abs(ca.bases["q"][:, 0] @ cb.bases["q"][:, 0]) < 0.8  # the bases differ
    x = np.r_[expm([0.05, -0.02, 0.0]) @ qa, tangent(qa, [0.1, 0.2, 0.0])]
    x[3:] = tangent(x[:3], x[3:])
    np.testing.assert_allclose(ca(0.0, x), cb(0.0, x), atol=1e-6)


def test_lqr_rejects_a_non_equilibrium():
    x_bad = np.r_[(R_UP @ expm(np.deg2rad(30.0) * np.array([1.0, 0, 0]))).reshape(9), np.zeros(3)]
    with pytest.raises(ValueError, match="equilibrium"):
        lqr(
            linearized_module("rigid_pendulum_pivot"),
            x_bad,
            np.zeros(3),
            P_PEND,
            1.0,
            1.0,
            module=standard_module("rigid_pendulum_pivot"),  # F2 fix: the check needs it
        )


# ---------------------------------------------------------------------------
# (5) inputs and x0
# ---------------------------------------------------------------------------


def test_input_lengths_are_checked():
    p, x0 = SETUPS["quadrotor_with_load"]
    mod = standard_module("quadrotor_with_load")
    with pytest.raises(ValueError, match="u must have 4 entries, got 3"):
        simulate(mod, x0, 0.01, p, u=np.zeros(3))
    with pytest.raises(ValueError, match=r"u\(t, state\) must have 4 entries, got 1"):
        simulate(mod, x0, 0.01, p, u=lambda t, x: 1.0)
    sf = modules(free_sphere_no_inputs)[0]
    assert sf.LAYOUT["n_inputs"] == 0
    x = np.r_[E3, 1.0, 0.0, 0.0]
    simulate(sf, x, 0.01, {"m": 1.0, "l": 1.0}, u=[])
    with pytest.raises(ValueError, match="u must have 0 entries, got 1"):
        simulate(sf, x, 0.01, {"m": 1.0, "l": 1.0}, u=[1.0])


def test_callable_input_sees_the_stage_state():
    """u(t, state) is evaluated at the RKMK4 stage states: a feedback u = −k ω on the
    gyroscopic pendulum dissipates energy monotonically."""
    mod = modules(gyro_spherical_pendulum)[0]
    traj = simulate(mod, np.r_[Q0, W0], 2.0, P_GYRO, u=lambda t, x: -2.0 * x[3:6], dt=1e-2)
    m, g, l = P_GYRO["m"], P_GYRO["g"], P_GYRO["l"]
    E = 0.5 * m * l * l * np.sum(traj.state[:, 3:] ** 2, axis=1) + m * g * l * traj.state[:, 2]
    assert np.all(np.diff(E) < 0)


def test_callable_input_returning_none_is_rejected():
    p, x0 = SETUPS["quadrotor_with_load"]

    def controller(t, x):
        np.r_[20.0, 0.0, 0.0, 0.0]  # forgot to return

    with pytest.raises((TypeError, ValueError)):
        simulate(standard_module("quadrotor_with_load"), x0, 0.01, p, u=controller)


def test_nan_input_is_reported_as_an_input_error():
    p, x0 = SETUPS["quadrotor_with_load"]
    with pytest.raises(ValueError, match=r"\bu\b") as exc:
        simulate(standard_module("quadrotor_with_load"), x0, 0.01, p, u=[20.0, 0.0, 0.0, np.nan])
    assert not isinstance(exc.value, np.linalg.LinAlgError)


def test_infinite_t_span_is_a_value_error():
    p, x0 = SETUPS["spherical_pendulum_s2"]
    with pytest.raises(ValueError):
        simulate(standard_module("spherical_pendulum_s2"), x0, np.inf, p)


def test_x0_tolerance_is_1e9_on_the_norm():
    p, x0 = SETUPS["quadrotor_with_load"]
    mod = standard_module("quadrotor_with_load")
    ok, bad = x0.copy(), x0.copy()
    ok[3:6] *= 1 + 5e-10
    bad[3:6] *= 1 + 2e-9
    simulate(mod, ok, 0.01, p)
    with pytest.raises(ValueError, match="S2 slot 'q' must have unit norm"):
        simulate(mod, bad, 0.01, p)
