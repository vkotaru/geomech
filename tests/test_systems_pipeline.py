"""Every system of ``geomech.systems`` runs through the whole pipeline (spec
``specs/systems-library.md``, verification items 1 and 4).

For each system at its defaults, and each parametric system also at sizes 1 and 2:
``compute_eom`` succeeds and passes ``check_eom`` (≤ 1e-5) at three seeds;
``to_standard_form`` with the system's inputs succeeds; ``StateLayout`` and
``generate_module`` write an importable module; ``simulate`` runs 1 s from ``x0`` with no
input and conserves the energy to 1e-8 (relative to ``max(1, |E₀|)``); the scene, when the
system has one, renders through ``Report``. Parametric systems derive in under 10 s at the
larger size. The derivation times are printed in the pytest summary
(``TIMINGS``, see ``tests/conftest.py``).
"""

from __future__ import annotations

import functools
import time

import numpy as np
import pytest

from geomech import check_eom, compute_eom, to_standard_form
from geomech.codegen import generate_module, load_module
from geomech.simulation import simulate
from geomech.viz import Report
from tests.test_systems import CASES, IDS, PARAMETRIC

TOL = 1e-5
SEEDS = (0, 1, 2)
ENERGY_TOL = 1e-8
T_FINAL = 1.0
MAX_DERIVE_SECONDS = 10.0  # the spec
DOCS_DERIVE_SECONDS = 2.0  # the bound docs/systems.md states

#: ``{case id: (dof, derivation seconds)}``, filled as the systems are derived.
TIMINGS: dict[str, tuple[int, float]] = {}


def _id(name, kwargs):
    return IDS[CASES.index((name, kwargs))]


@functools.cache
def _derived(name, items):
    from geomech.systems import get

    s = get(name, **dict(items))
    t0 = time.perf_counter()
    eom = compute_eom(s.L, s.dW, s.variables)
    TIMINGS[_id(name, dict(items))] = (s.dof, time.perf_counter() - t0)
    sf = to_standard_form(eom, s.variables, list(s.inputs))
    return s, eom, sf


def derived(name, kwargs):
    return _derived(name, tuple(sorted(kwargs.items())))


@functools.cache
def _module(name, items, tmp):
    s, eom, sf = _derived(name, items)
    layout = s.layout()
    path = generate_module(sf, layout, f"{tmp}/{_id(name, dict(items))}.py", lagrangian=s.L)
    return load_module(path, f"systems_{name}"), layout


@pytest.fixture(scope="module")
def outdir(tmp_path_factory):
    return str(tmp_path_factory.mktemp("systems"))


@pytest.mark.parametrize("name, kwargs", CASES, ids=IDS)
def test_derive_and_check(name, kwargs):
    s, eom, _ = derived(name, kwargs)
    assert len(eom) == len(s.variables.scalars + s.variables.vectors + s.variables.matrices)
    for seed in SEEDS:
        err = check_eom(s.L, s.dW, s.variables, eom, inputs=list(s.inputs), seed=seed)
        assert err < TOL, f"{name} {kwargs}: check_eom {err:.2e} at seed {seed}"


@pytest.mark.parametrize("name, kwargs", CASES, ids=IDS)
def test_standard_form_has_every_input(name, kwargs):
    s, eom, sf = derived(name, kwargs)
    assert set(sf) == set(eom)
    used = {u for eq in sf.values() for u in eq.G}
    assert used == {str(u) for u in s.inputs}, f"{name}: inputs without a G block"


@pytest.mark.parametrize("name, kwargs", CASES, ids=IDS)
def test_generate_simulate_conserve_energy(name, kwargs, outdir):
    s, _, _ = derived(name, kwargs)
    module, layout = _module(name, tuple(sorted(kwargs.items())), outdir)
    assert module.LAYOUT["n_state"] == s.x0.size == layout.n_state
    traj = simulate(module, s.x0, T_FINAL, s.params)
    assert traj.t[-1] == T_FINAL and np.all(np.isfinite(traj.state))
    E = traj.energy()
    drift = float(np.max(np.abs(E - E[0]))) / max(1.0, abs(float(E[0])))
    assert drift < ENERGY_TOL, f"{name}: relative energy drift {drift:.2e}"
    assert traj.max_drift() < 1e-10


@pytest.mark.parametrize("name, kwargs", CASES, ids=IDS)
def test_scene_renders_through_report(name, kwargs, outdir):
    s, eom, sf = derived(name, kwargs)
    if not s.has_scene:
        assert s.scene() is None
        pytest.skip(f"{name} has no scene")
    module, layout = _module(name, tuple(sorted(kwargs.items())), outdir)
    traj = simulate(module, s.x0, 0.2, s.params, dt=2e-3)
    scene = s.scene(layout)
    assert scene.bodies or scene.points
    report = Report(name, s.as_tuple(), eom=eom, standard_form=sf, scene=scene, trajectory=traj)
    html = report.html()
    assert "<html" in html.lower() and report.check_error < TOL


def test_quadrotors_are_drawn_as_solid_bodies():
    from geomech.viz import Box, Cylinder

    for name in ("quadrotor", "quadrotor_suspended_load", "multi_quadrotor_rigid_load"):
        s, _, _ = derived(name, {})
        scene = s.scene()
        shapes = [type(b.shape) for b in scene.bodies]
        assert shapes.count(Box) >= 3 and shapes.count(Cylinder) >= 4, (name, shapes)


@pytest.mark.parametrize("name, param", list(PARAMETRIC.items()))
def test_parametric_systems_derive_in_under_ten_seconds(name, param):
    """The larger tested size (the default, 3) derives in under 10 s (the spec), and in
    under 2 s (the claim of ``docs/systems.md``)."""
    from geomech.systems import get

    s = get(name)
    assert s.options[param] == 3
    t0 = time.perf_counter()
    compute_eom(s.L, s.dW, s.variables)
    elapsed = time.perf_counter() - t0
    assert elapsed < min(MAX_DERIVE_SECONDS, DOCS_DERIVE_SECONDS), f"{name}: {elapsed:.1f} s"
