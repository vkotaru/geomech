"""Regression tests for the final review of the visualizer and skill
(``specs/reviews/viz-skill-final.md``): recorded inputs, controller errors, the energy
label without recorded inputs, ``StateLayout.slot_names`` and the page's ready hook."""

from __future__ import annotations

import copy

import numpy as np
import pytest

from geomech import S2, SO3, Scalar, SystemVariables, Vector
from geomech.codegen import StateLayout
from geomech.simulation import simulate
from geomech.viz import PointMass, Report, Scene
from tests.codegen_systems import derived
from tests.simulation_systems import BUILDERS, SETUPS, standard_module


def _pendulum():
    p, x0 = SETUPS["spherical_pendulum_s2"]
    return p, x0, standard_module("spherical_pendulum_s2")


@pytest.mark.parametrize("output_every", [1, 7])
def test_controller_is_called_once_per_stage_and_recorded(output_every):
    p, x0, mod = _pendulum()
    calls = []

    def ctrl(t, x):
        calls.append(t)
        return np.array([0.0, 0.5 * len(calls), 0.0])  # stateful: depends on the call count

    traj = simulate(mod, x0, 0.05, p, u=ctrl, dt=1e-2, output_every=output_every)
    n = traj.info["n_steps"]
    assert len(calls) == 4 * n + 1  # four stages per step, plus the final state
    # every stored row is a value the controller returned during the run
    returned = 0.5 * np.arange(1, len(calls) + 1)
    assert np.isin(traj.u[:, 1], returned).all()
    assert traj.u[0, 1] == 0.5  # the first call, shared by the singularity check and K1
    assert traj.u[-1, 1] == returned[-1]


def test_controller_error_propagates_at_once_with_the_time():
    p, x0, mod = _pendulum()
    calls = []

    def ctrl(t, x):
        calls.append(t)
        if t > 0.02:
            raise KeyError("gain table has no entry")
        return np.zeros(3)

    with pytest.raises(KeyError, match="gain table") as exc:
        simulate(mod, x0, 1.0, p, u=ctrl, dt=1e-2)
    assert max(calls) < 0.04  # raised during the run, not after all 100 steps
    notes = getattr(exc.value, "__notes__", [])
    assert any("u(t, state) at t = 0.025" in n for n in notes), notes


def test_controller_error_with_solve_ivp_keeps_its_type():
    pytest.importorskip("scipy")
    p, x0, mod = _pendulum()

    def ctrl(t, x):
        raise IndexError("no such channel")

    with pytest.raises(IndexError, match="no such channel"):
        simulate(mod, x0, 0.05, p, u=ctrl, method="solve_ivp")


def test_solve_ivp_records_inputs():
    pytest.importorskip("scipy")
    p, x0, mod = _pendulum()
    traj = simulate(mod, x0, 0.05, p, u=lambda t, x: np.array([0.0, t, 0.0]),
                    method="solve_ivp", dt=1e-2)  # fmt: skip
    np.testing.assert_allclose(traj.u[:, 1], traj.t)


def test_unrecorded_inputs_are_not_called_unforced():
    d = derived(BUILDERS["spherical_pendulum_s2"])
    p, x0, mod = _pendulum()
    traj = simulate(mod, x0, 0.2, p, output_every=10)
    blind = copy.copy(traj)
    blind.u = None
    scene = Scene(d.layout)
    scene.add(PointMass(d.system.symbols["q"] * d.system.symbols["l"]))
    sysd = tuple(d.system)[:4]
    rep = Report("x", system=sysd, eom=d.eom, standard_form=d.sf, scene=scene,
                 trajectory=blind)  # fmt: skip
    assert not rep.forced()
    page = rep.html()
    assert "Energy change" in page and "Energy drift" not in page
    assert "inputs not recorded" in page and "unforced" not in page
    known = Report("x", system=sysd, eom=d.eom, standard_form=d.sf, scene=scene,
                   trajectory=traj).html()  # fmt: skip
    assert "Energy drift" in known and "unforced" in known


def test_slot_names_follow_variable_names():
    variables = SystemVariables(
        scalars=[Scalar("s")], vectors=[Vector("x"), S2("q1")], matrices=[SO3("R2")]
    )
    layout = StateLayout(variables, [], [])
    names = layout.slot_names()
    assert names == ["s", "x", "q1", "R2", "s_dot", "x_dot", "omega_q1", "Omega_R2"]
    state = layout.pack(s=1.0, x=[1, 2, 3], q1=[0, 0, 1], R2=np.eye(3), omega_q1=[1, 0, 0])
    assert list(layout.unpack(state)) == names


def test_page_has_the_ready_hook():
    d = derived(BUILDERS["spherical_pendulum_s2"])
    p, x0, mod = _pendulum()
    traj = simulate(mod, x0, 0.1, p, output_every=10)
    scene = Scene(d.layout)
    scene.add(PointMass(d.system.symbols["q"] * d.system.symbols["l"]))
    page = Report("x", system=tuple(d.system)[:4], eom=d.eom, standard_form=d.sf,
                  scene=scene, trajectory=traj).html()  # fmt: skip
    for state in ("ready", "no-webgl", "no-three", "none"):
        assert f'"{state}"' in page
    assert "data-gm-viewer" in page
