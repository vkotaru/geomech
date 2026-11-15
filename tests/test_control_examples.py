"""The closed-loop LQR examples of ``examples/control/`` (a fast version of
``python examples/control/run.py --all``).

For every example: the equilibrium residual is at round-off, the controllability claim
holds (full rank, except the documented cases), the closed loop of the controllable part
is Hurwitz, and the nonlinear simulation reduces the normalized error (the square root of
the LQR cost-to-go, relative to its initial value) at least tenfold and stays converged
over the last 10 % of the horizon. The simulations use a coarser step than the examples'
1 ms (``dt = min(5 ms, 2 / |λ|max)``, inside the RKMK4 stability region of the fastest
closed-loop eigenvalue) over the example's own horizon. ``tracking.py`` is smoke-tested on
short horizons, and the README table is checked against the stored results.
"""

from __future__ import annotations

import functools
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("scipy")

HERE = Path(__file__).resolve().parents[1] / "examples" / "control"


def _load(name, key):
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    if key in sys.modules:
        return sys.modules[key]
    spec = importlib.util.spec_from_file_location(key, HERE / f"{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[key] = mod
    spec.loader.exec_module(mod)
    return mod


# run.py imports the helper as ``control_helper``: register it under that name first
helper = _load("control_helper", "control_helper")
examples = _load("run", "control_examples_run")
NAMES = [ex.name for ex in examples.EXAMPLES]

#: name → (rank deficiency, uncontrollable eigenvalues as a function of the parameters)
NOT_CONTROLLABLE = {
    # the pendulum's swing out of the rail's vertical plane: ±√(g/l)
    "cart_spherical_pendulum_inverted_in_plane": lambda p: [
        np.sqrt(p["g"] / p["l"]),
        -np.sqrt(p["g"] / p["l"]),
    ],
    "cart_spherical_pendulum_out_of_plane": lambda p: [
        np.sqrt(p["g"] / p["l"]),
        -np.sqrt(p["g"] / p["l"]),
    ],
    # the three components of the conserved angular momentum
    "spacecraft_reaction_wheels_slew": lambda p: [0.0, 0.0, 0.0],
    "spacecraft_reaction_wheels_s1_slew": lambda p: [0.0, 0.0, 0.0],
    "spacecraft_solar_arrays_slew": lambda p: [0.0, 0.0, 0.0],
}
EXPECTED_TO_FAIL = {"cart_spherical_pendulum_out_of_plane"}


@functools.cache
def designed(name):
    ex = examples.BY_NAME[name]
    sys_, layout, plant, plant_lin, _ = helper.build(ex.system, ex.kwargs)
    p = dict(sys_.params)
    x_eq = layout.pack(**ex.equilibrium(sys_, layout))
    u_eq, residual = helper.equilibrium_input(plant, x_eq, p)
    d = helper.design(ex, sys_, layout, plant, plant_lin, x_eq, u_eq)
    return ex, sys_, residual, d


@pytest.mark.parametrize("name", NAMES)
def test_equilibrium_and_controllability(name):
    ex, sys_, residual, d = designed(name)
    assert residual < 1e-10, f"{name}: equilibrium residual {residual:.2e}"
    if name in NOT_CONTROLLABLE:
        expected = np.sort(np.array(NOT_CONTROLLABLE[name](sys_.params), dtype=float))
        got = np.sort(np.array([z.real for z in d.uncontrollable]))
        assert d.rank == d.n_reduced - len(expected)
        assert np.allclose(got, expected, atol=1e-6), (got, expected)
        assert np.allclose([z.imag for z in d.uncontrollable], 0, atol=1e-6)
    else:
        assert d.rank == d.n_reduced and not d.uncontrollable
        assert d.used_library_lqr or ex.drop
    assert d.cl_max_real < 0, f"{name}: closed loop not Hurwitz ({d.cl_max_real})"


@pytest.mark.parametrize("name", NAMES)
def test_nonlinear_closed_loop_reduces_the_error(name):
    ex, sys_, residual, d = designed(name)
    dt = min(5e-3, 2.0 / d.cl_fastest)
    horizon = ex.horizon
    m = helper.run(ex, horizon=horizon, dt=dt, save=False)
    if name in EXPECTED_TO_FAIL:
        assert not m["converged"]
        return
    assert m["final_normalized_error"] <= 0.1, (
        f"{name}: normalized error {m['final_normalized_error']:.3g} after {horizon:.1f} s"
    )
    assert m["converged"], f"{name}: not settled over the last 10 % of the horizon"
    assert m["constraint_drift"] < 1e-10


def test_settling_definition():
    t = np.linspace(0.0, 10.0, 1001)
    nu = np.exp(-t)
    ts, ok = helper.settling(t, nu)
    assert ok and abs(ts - np.log(50)) < 0.02
    # re-crossing the band late: settled only after the last exit, not converged if it
    # is still above the band in the final 10 %
    nu2 = nu.copy()
    nu2[950] = 0.05
    ts2, ok2 = helper.settling(t, nu2)
    assert ts2 == t[951] and not ok2
    nu3 = nu.copy()
    nu3[600] = 0.05
    ts3, ok3 = helper.settling(t, nu3)
    assert ts3 == t[601] and ok3


def test_readme_table_is_the_stored_results():
    """``results/table.md`` is what ``run.py`` produces from the stored JSON, and the
    README contains it."""
    import json

    stored = [
        json.loads((HERE / "results" / f"{ex.name}.json").read_text()) for ex in examples.EXAMPLES
    ]
    table = examples.table(stored)
    assert (HERE / "results" / "table.md").read_text() == table
    readme = (HERE / "README.md").read_text()
    assert table in readme
    for ex in examples.EXAMPLES:
        assert f"`{ex.name}`" in readme


tracking = _load("tracking", "control_examples_tracking")


def test_tracking_figure_eight_smoke():
    """Short horizon, coarse step: the flat reference is kinematically consistent
    (independent finite-difference and integral checks) and TV-LQR reduces the position
    error tenfold."""
    m = tracking.quadrotor_figure_eight(horizon=4.0, dt=5e-3, save=False)
    checks = m["flatness_checks"]
    assert checks["kinematic_error_Rdot"] < 1e-8
    assert checks["integral_error_Omega"] < 1e-6
    assert m["final_position_error"] < 0.1 * m["initial_position_error"]


def test_tracking_suspended_load_smoke():
    """The corrected time-varying reduction matches finite differences of the nonlinear
    error dynamics along the reference much better than the equilibrium reduction, and
    TV-LQR reduces the load position error tenfold."""
    ref, layout, plant, plant_lin, p = tracking.load_reference(6.0, dt=5e-3)
    n = len(ref.t)
    errs = tracking.validate_reduction(plant, plant_lin, ref, p, [n // 2, 3 * n // 4])
    full = errs["+ curvature term (all)"]
    assert full < 1e-5 and full < 1e-2 * errs["equilibrium reduction only"]
    m = tracking.suspended_load_circle(horizon=6.0, dt=5e-3, save=False, validate=False)
    assert m["final_position_error"] < 0.1 * m["initial_position_error"]


def test_s1_spacecraft_slew_matches_the_r3_slew():
    """The S¹ wheel model's slew has the R³ model's closed loop: same eigenvalues, and the
    same settling time and normalized error, both in the stored results and in a run."""
    import json

    a, b = (
        json.loads((HERE / "results" / f"{n}.json").read_text())
        for n in ("spacecraft_reaction_wheels_slew", "spacecraft_reaction_wheels_s1_slew")
    )
    for key in ("settling_time", "closed_loop_max_real", "closed_loop_min_damping"):
        assert a[key] == pytest.approx(b[key], rel=1e-9, abs=1e-12), key
    assert b["final_normalized_error"] == pytest.approx(a["final_normalized_error"], rel=1e-6)
    runs = [
        helper.run(examples.BY_NAME[n], horizon=12.0, dt=5e-3, save=False)
        for n in ("spacecraft_reaction_wheels_slew", "spacecraft_reaction_wheels_s1_slew")
    ]
    assert runs[0]["settling_time"] == runs[1]["settling_time"]
    assert runs[1]["final_normalized_error"] == pytest.approx(
        runs[0]["final_normalized_error"], rel=1e-6
    )


def test_solar_arrays_slew_ends_at_the_commanded_panel_angle():
    """The panels end at the commanded sun-pointing angle about b₁ and the bus at R = I,
    with the wheel torques below 0.2 N·m and the drive torques below 0.02 N·m."""
    import json

    name = "spacecraft_solar_arrays_slew"
    ex = examples.BY_NAME[name]
    sys_, layout, plant, _, _ = helper.build(ex.system, ex.kwargs)
    from geomech.simulation import simulate

    _, _, _, d = designed(name)
    x0 = ex.perturbation(sys_, layout, layout.pack(**ex.equilibrium(sys_, layout)))
    traj = simulate(plant, x0, ex.horizon, dict(sys_.params), u=d.ctrl, dt=5e-3)
    b1 = np.asarray(sys_.params["b1"])
    for k in (1, 2):
        n = traj[f"n_{k}"][-1]
        angle = np.degrees(np.arctan2(np.cross(helper.UP, n) @ b1, n @ helper.UP))
        assert abs(angle - examples.SUN_ANGLE) < 0.05, (k, angle)
    assert np.allclose(traj["R"][-1], np.eye(3), atol=1e-3)
    m = json.loads((HERE / "results" / f"{name}.json").read_text())
    peaks = {k: v["peak_deviation"] for k, v in m["peak_inputs"].items()}
    assert max(peaks[f"\\tau_{i}"] for i in (1, 2, 3)) < 0.2
    assert max(peaks[f"\\mu_{k}"] for k in (1, 2)) < 0.02
