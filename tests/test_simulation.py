"""``geomech.simulation``: accelerations, state checks, the ``simulate`` API, ``log_so3``.

- ``accelerations`` on every S² system of the oracle gate: ``ω̇`` is tangent
  (``|q·ω̇| < 1e-12``) and the tangent-plane equation ``P_q (M a + f + G u) = 0`` holds to
  1e-10, whether the derived S² mass block is singular (quadrotor with load) or full rank
  (spherical pendulum). A singular system raises ``LinAlgError`` with the condition
  number and the state.
- ``x0`` is checked (unit ``q``, ``ω ⊥ q``, ``R`` orthogonal with ``det R = +1``): each
  violation raises ``ValueError`` naming the slot.
- ``u`` as ``None``, a constant or a callable; ``Trajectory`` slot access and energy.
- ``simulate`` runs 1 s on each of the nine systems of ``tests/vbl_systems.py`` with both
  methods.
- ``import geomech.simulation`` does not import scipy.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from geomech.numeric import expm_so3, log_so3
from geomech.simulation import accelerations, check_state, simulate
from tests.codegen_systems import ALL_SYSTEMS, derived, modules, sample_points
from tests.simulation_systems import SETUPS, linearized_module, standard_module, unit

ROOT = Path(__file__).resolve().parents[1]
_ids = dict(ids=lambda f: f.__name__)

S2_SYSTEMS = [
    b for b in ALL_SYSTEMS if any(s.manifold == "S2" for s in derived(b).layout.configuration)
]


def _physical(p):
    """Positive scalars and positive-definite symmetric matrices (random values kept)."""
    out = {}
    for key, val in p.items():
        if np.shape(val) == (3, 3):
            out[key] = np.diag([1.0, 2.0, 3.0]) + 0.1 * (val + val.T)
        elif np.shape(val) == ():
            out[key] = 1.0 + abs(val)
        else:
            out[key] = val
    return out


# ---------------------------------------------------------------------------
# accelerations
# ---------------------------------------------------------------------------


def test_there_are_s2_systems_of_both_kinds():
    assert len(S2_SYSTEMS) >= 10
    names = {b.__name__ for b in S2_SYSTEMS}
    assert {"spherical_pendulum_s2", "quadrotor_with_load"} <= names


#: Kinetic energy ½ m |ẋ + l q̇|²: the velocity ẋ = −l q̇ costs nothing, so the mass
#: matrix is singular for every parameter value and ``accelerations`` must raise.
DEGENERATE = {"adv_r3_s2_pendulum_with_gyro"}


@pytest.mark.parametrize("build", S2_SYSTEMS, **_ids)
def test_accelerations_tangent_and_equation_on_tangent_plane(build):
    d = derived(build)
    mod = modules(build)[0]
    rng = np.random.default_rng(3)
    for pt in sample_points(d, n=5):
        p = _physical(pt.p)
        u = rng.normal(size=d.layout.n_inputs)
        if build.__name__ in DEGENERATE:
            with pytest.raises(np.linalg.LinAlgError, match="condition number"):
                accelerations(mod, pt.state, u, p)
            continue
        a = accelerations(mod, pt.state, u, p)
        residual = mod.mass(pt.state, p) @ a + mod.bias(pt.state, p)
        residual = residual + mod.input_matrix(pt.state, p) @ u
        scale = max(1.0, float(np.max(np.abs(mod.bias(pt.state, p)))), float(np.abs(a).max()))
        for cfg, acc in zip(d.layout.configuration, d.layout.acceleration):
            if cfg.manifold != "S2":
                assert np.max(np.abs(residual[acc.slice])) < 1e-10 * scale
                continue
            q = pt.state[cfg.slice]
            assert abs(q @ a[acc.slice]) < 1e-12 * max(1.0, np.abs(a[acc.slice]).max())
            P = np.eye(3) - np.outer(q, q)
            assert np.max(np.abs(P @ residual[acc.slice])) < 1e-10 * scale


def test_s2_mass_block_is_singular_or_full_rank_and_both_solve():
    """The derived S² mass block is −m l² I for the spherical pendulum (rank 3) and
    Hat(q) Hat(q) m l² for the quadrotor with load (rank 2); the projected solve handles
    both, while the plain 3×3 solve gives q·ω̇ ≠ 0 or fails."""
    p, x = SETUPS["spherical_pendulum_s2"]
    mod = standard_module("spherical_pendulum_s2")
    assert np.linalg.matrix_rank(mod.mass(x, p)) == 3
    p2, x2 = SETUPS["quadrotor_with_load"]
    mod2 = standard_module("quadrotor_with_load")
    M2 = mod2.mass(x2, p2)
    assert np.linalg.matrix_rank(M2[3:6, 3:6]) == 2
    u2 = np.array([20.0, 0.1, -0.2, 0.3])
    a2 = accelerations(mod2, x2, u2, p2)
    assert abs(x2[3:6] @ a2[3:6]) < 1e-12
    # without the projection the full system is singular
    assert np.linalg.matrix_rank(M2) < M2.shape[0]


def test_accelerations_point_mass_and_rigid_body_closed_form():
    p, x = SETUPS["point_mass_r3"]
    mod = standard_module("point_mass_r3")
    F = np.array([1.0, -2.0, 0.5])
    np.testing.assert_allclose(accelerations(mod, x, F, p), F / p["m"], rtol=1e-14)
    p, x = SETUPS["free_rigid_body"]
    mod = standard_module("free_rigid_body")
    W, J, Mt = x[9:], p["J"], np.array([0.1, 0.2, -0.3])
    ref = np.linalg.solve(J, np.cross(J @ W, W) + Mt)
    np.testing.assert_allclose(accelerations(mod, x, Mt, p), ref, rtol=1e-13, atol=1e-14)


def test_accelerations_singular_system_raises_with_condition_number():
    mod = standard_module("point_mass_r3")
    x = SETUPS["point_mass_r3"][1]
    with pytest.raises(np.linalg.LinAlgError, match=r"singular \(condition number .*state"):
        accelerations(mod, x, None, {"m": 0.0})


def test_simulate_singular_system_raises_before_integrating():
    p, x = SETUPS["free_rigid_body"]
    J = np.diag([1.0, 1.0, 0.0])
    with pytest.raises(np.linalg.LinAlgError, match="condition number"):
        simulate(standard_module("free_rigid_body"), x, 1.0, {"J": J})


def test_accelerations_size_errors():
    p, x = SETUPS["spherical_pendulum_s2"]
    mod = standard_module("spherical_pendulum_s2")
    with pytest.raises(ValueError, match="state must have 6 entries, got 5"):
        accelerations(mod, x[:5], None, p)
    with pytest.raises(ValueError, match="u must have 3 entries, got 2"):
        accelerations(mod, x, [1.0, 2.0], p)


def test_accelerations_needs_a_standard_form_module():
    lin = linearized_module("spherical_pendulum_s2")
    with pytest.raises(TypeError, match="no mass\\(\\); generate it from to_standard_form"):
        accelerations(lin, SETUPS["spherical_pendulum_s2"][1], None, {})
    with pytest.raises(TypeError, match="generated by geomech.codegen"):
        accelerations(object(), np.zeros(6), None, {})


# ---------------------------------------------------------------------------
# x0 checks
# ---------------------------------------------------------------------------


def _bad(name, edit):
    p, x = SETUPS[name]
    x = x.copy()
    edit(x)
    return standard_module(name), x, p


@pytest.mark.parametrize(
    "name, edit, message",
    [
        (
            "spherical_pendulum_s2",
            lambda x: x.__setitem__(slice(0, 3), x[0:3] * 1.01),
            r"x0: the S2 slot 'q' must have unit norm, got \|q\| = 1\.01",
        ),
        (
            "spherical_pendulum_s2",
            lambda x: x.__setitem__(slice(3, 6), x[3:6] + 1e-6 * x[0:3]),
            r"x0: the velocity slot '\\omega_\{q\}' must be orthogonal to 'q'",
        ),
        (
            "elastic_spherical_pendulum",
            lambda x: x.__setitem__(slice(1, 4), x[1:4] * 0.5),
            r"the S2 slot 'q' must have unit norm",
        ),
        (
            "free_rigid_body",
            lambda x: x.__setitem__(slice(0, 9), x[0:9] + 1e-6),
            r"x0: the SO3 slot 'R' must be orthogonal",
        ),
        (
            "free_rigid_body",
            lambda x: x.__setitem__(slice(0, 9), -x[0:9]),
            r"x0: the SO3 slot 'R' must have determinant \+1, got -1",
        ),
        (
            "double_rigid_pendulum",
            lambda x: x.__setitem__(slice(9, 18), np.eye(3).reshape(9) * 2),
            r"the SO3 slot 'R2' must be orthogonal",
        ),
        (
            "quadrotor_with_load",
            lambda x: x.__setitem__(18, x[18] + 0.1),
            r"the velocity slot '\\omega_\{q\}' must be orthogonal to 'q'",
        ),
        ("point_mass_r3", lambda x: x.__setitem__(0, np.nan), r"x0: non-finite entries \[0\]"),
    ],
)
def test_x0_validation_names_the_slot(name, edit, message):
    mod, x, p = _bad(name, edit)
    with pytest.raises(ValueError, match=message):
        simulate(mod, x, 0.01, p)


def test_x0_size_and_tolerance():
    p, x = SETUPS["spherical_pendulum_s2"]
    mod = standard_module("spherical_pendulum_s2")
    with pytest.raises(ValueError, match="x0 must have 6 entries, got 7"):
        simulate(mod, np.r_[x, 0.0], 0.01, p)
    y = x.copy()
    y[:3] *= 1 + 5e-10  # inside the 1e-9 tolerance
    check_state(mod, y)
    with pytest.raises(ValueError, match="unit norm"):
        check_state(mod, y, tol=1e-12)


# ---------------------------------------------------------------------------
# simulate API
# ---------------------------------------------------------------------------


def test_simulate_argument_errors():
    p, x = SETUPS["point_mass_r3"]
    mod = standard_module("point_mass_r3")
    with pytest.raises(ValueError, match="method must be one of"):
        simulate(mod, x, 1.0, p, method="euler")
    with pytest.raises(ValueError, match="tf > t0"):
        simulate(mod, x, (1.0, 0.5), p)
    with pytest.raises(ValueError, match="dt must be positive"):
        simulate(mod, x, 1.0, p, dt=0.0)
    with pytest.raises(ValueError, match="u must have 3 entries, got 2"):
        simulate(mod, x, 1.0, p, u=[1.0, 2.0])
    with pytest.raises(ValueError, match=r"u\(t, state\) must have 3 entries, got 1"):
        simulate(mod, x, 1.0, p, u=lambda t, s: 1.0)
    with pytest.raises(KeyError, match="missing parameter 'm'"):
        simulate(mod, x, 1.0, {})


@pytest.mark.parametrize("method", ["rkmk4", "solve_ivp"])
def test_point_mass_inputs_constant_and_callable(method):
    """ẍ = F/m: constant F gives the parabola, F(t) = c t the cubic (both exact for RK4)."""
    p, x0 = SETUPS["point_mass_r3"]
    mod = standard_module("point_mass_r3")
    F = np.array([1.0, -2.0, 0.5])
    tr = simulate(mod, x0, (0.5, 1.5), p, u=F, method=method, dt=0.01)
    assert tr.t[0] == 0.5 and tr.t[-1] == 1.5 and len(tr.t) == 101
    tau = tr.t - 0.5
    ref = x0[:3] + np.outer(tau, x0[3:]) + 0.5 * np.outer(tau**2, F / p["m"])
    np.testing.assert_allclose(tr["x"], ref, atol=1e-12)
    tr = simulate(mod, x0, 1.0, p, u=lambda t, s: F * t, method=method, dt=0.01)
    ref = x0[:3] + np.outer(tr.t, x0[3:]) + np.outer(tr.t**3 / 6, F / p["m"])
    np.testing.assert_allclose(tr["x"], ref, atol=1e-12)
    np.testing.assert_allclose(tr["x_dot"][-1], x0[3:] + 0.5 * F / p["m"], atol=1e-12)
    assert tr.constraint_drift == {} and tr.max_drift() == 0.0


def test_dt_is_shortened_to_divide_the_span():
    p, x0 = SETUPS["spring_mass_r"]
    tr = simulate(standard_module("spring_mass_r"), x0, 0.25, p, dt=0.1)
    assert len(tr.t) == 4 and tr.t[-1] == 0.25
    assert tr.info["dt"] == pytest.approx(0.25 / 3)


def test_trajectory_slots_and_energy():
    p, x0 = SETUPS["quadrotor_with_load"]
    mod = standard_module("quadrotor_with_load")
    tr = simulate(mod, x0, 0.05, p)
    n = len(tr.t)
    assert tr["x_L"].shape == (n, 3) and tr["q"].shape == (n, 3)
    assert tr["R"].shape == (n, 3, 3) and tr["\\Omega_{R}"].shape == (n, 3)
    np.testing.assert_array_equal(tr["R"][0], x0[6:15].reshape(3, 3))
    with pytest.raises(KeyError, match="no state slot 'omega'"):
        tr["omega"]
    E = tr.energy()
    np.testing.assert_allclose(tr.energy(lambda s, p: mod.energy(s, p) + 1.0), E + 1.0)
    assert set(tr.constraint_drift) == {"|q|-1", "q.omega_q", "R^T R-I"}
    p1, x1 = SETUPS["spring_mass_r"]
    tr1 = simulate(standard_module("spring_mass_r"), x1, 0.01, p1)
    assert tr1["x"].shape == (len(tr1.t),)
    tr1.module = type("NoEnergy", (), {"LAYOUT": tr1.module.LAYOUT})
    with pytest.raises(TypeError, match="no energy\\(\\)"):
        tr1.energy()


@pytest.mark.parametrize("name", list(SETUPS))
@pytest.mark.parametrize("method", ["rkmk4", "solve_ivp"])
def test_simulate_one_second_every_vbl_system(name, method):
    p, x0 = SETUPS[name]
    tr = simulate(standard_module(name), x0, 1.0, p, method=method)
    assert tr.state.shape == (1001, x0.size) and np.all(np.isfinite(tr.state))
    assert tr.max_drift() < (1e-12 if method == "rkmk4" else 1e-8)
    check_state(tr.module, tr.state[-1], tol=1e-8)


def test_solve_ivp_output_is_projected_and_drift_reported():
    p, x0 = SETUPS["spherical_pendulum_s2"]
    mod = standard_module("spherical_pendulum_s2")
    tr = simulate(mod, x0, 2.0, p, method="solve_ivp", rtol=1e-6)
    assert tr.constraint_drift["|q|-1"].max() > 1e-9  # the raw solution left the sphere
    np.testing.assert_allclose(np.linalg.norm(tr["q"], axis=1), 1.0, atol=1e-15)
    assert tr.info["nfev"] > 0 and tr.info["ivp_method"] == "DOP853"


def test_solve_ivp_without_scipy_has_install_hint(monkeypatch):
    p, x0 = SETUPS["spring_mass_r"]
    monkeypatch.setitem(sys.modules, "scipy.integrate", None)
    with pytest.raises(ImportError, match=r"geomech\[sim\]"):
        simulate(standard_module("spring_mass_r"), x0, 0.1, p, method="solve_ivp")


def test_import_simulation_does_not_import_scipy():
    code = (
        "import sys, geomech, geomech.simulation;"
        " bad = [m for m in ('scipy', 'casadi') if m in sys.modules];"
        " assert not bad, bad"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=ROOT)
    assert proc.returncode == 0, proc.stderr


# ---------------------------------------------------------------------------
# log_so3
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "angle", [0.0, 1e-9, 1e-5, 1e-3, 0.5, 2.0, 3.0, np.pi - 1e-4, np.pi - 1e-8]
)
def test_log_so3_inverts_expm(angle):
    rng = np.random.default_rng(int(angle * 1e3) % 97)
    for _ in range(20):
        v = angle * unit(rng.normal(size=3))
        w = log_so3(expm_so3(v))
        assert w.shape == (3,)
        np.testing.assert_allclose(w, v, atol=1e-7 if angle > 3.1 else 1e-12)
        np.testing.assert_allclose(expm_so3(w), expm_so3(v), atol=1e-13)


def test_log_so3_at_pi():
    for axis in np.eye(3).tolist() + [unit([1.0, 2.0, -2.0]).tolist()]:
        R = expm_so3(np.pi * np.array(axis))
        w = log_so3(R)
        assert np.linalg.norm(w) == pytest.approx(np.pi)
        np.testing.assert_allclose(expm_so3(w), R, atol=1e-12)
