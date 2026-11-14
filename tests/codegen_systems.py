"""Systems and helpers shared by the code generation tests (``tests/test_codegen*.py``).

``ALL_SYSTEMS`` is every builder of the numeric oracle gate (``tests/test_eom_oracle.py``,
which includes the nine of ``tests/vbl_systems.py`` and the pinned phase-1 systems) plus the
adversarial systems of ``tests/phase1_systems.py``. ``derived(build)`` runs the pipeline
once per builder and caches the result.
"""

from __future__ import annotations

import functools
from typing import NamedTuple

import numpy as np

from geomech import compute_eom, linearize, to_standard_form
from geomech.codegen import StateLayout, collect_params, load_module, module_source
from geomech.numeric import Evaluator, sample_state
from geomech.utils.errors import AlgebraicError
from tests import phase1_systems, vbl_systems
from tests.test_eom_oracle import SYSTEMS as ORACLE_SYSTEMS

VBL_SYSTEMS = [
    vbl_systems.point_mass_r3,
    vbl_systems.spring_mass_r,
    vbl_systems.free_rigid_body,
    vbl_systems.rigid_pendulum_pivot,
    vbl_systems.spherical_pendulum_s2,
    vbl_systems.quadrotor,
    vbl_systems.double_rigid_pendulum,
    vbl_systems.elastic_spherical_pendulum,
    vbl_systems.quadrotor_with_load,
]

ALL_SYSTEMS = list(dict.fromkeys([*ORACLE_SYSTEMS, *phase1_systems.ADVERSARIAL_SYSTEMS]))

N_STATES = 20
RTOL = 1e-12


class Derived(NamedTuple):
    system: vbl_systems.System
    eom: dict
    sf: dict
    lin: object  # LinearizedSystem, or None when linearize does not apply
    layout: StateLayout


@functools.cache
def derived(build) -> Derived:
    s = build()
    eom = compute_eom(s.L, s.dW, s.variables)
    sf = to_standard_form(eom, s.variables, s.inputs)
    try:
        lin = linearize(eom, s.variables, s.inputs)
    except AlgebraicError:
        lin = None  # a non-constant symbol that is not a configuration variable (frozen R)
    exprs = [s.L] + [eq for _, eq in eom.values()]
    layout = StateLayout(s.variables, s.inputs, collect_params(exprs, s.variables, s.inputs))
    return Derived(s, eom, sf, lin, layout)


@functools.cache
def modules(build, target="numpy"):
    """``(standard form, eom, linearized or None)`` generated modules of *build*."""
    d = derived(build)
    sf = load_module(module_source(d.sf, d.layout, target, lagrangian=d.system.L))
    eom = load_module(module_source(d.eom, d.layout, target))
    lin = None if d.lin is None else load_module(module_source(d.lin, d.layout, target))
    return sf, eom, lin


class Point(NamedTuple):
    env: dict
    state: np.ndarray
    a: np.ndarray
    u: np.ndarray
    p: dict


def sample_points(d: Derived, n=N_STATES, seed=0) -> list[Point]:
    """*n* admissible states (``sample_state``) with random parameter values."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        env = sample_state(d.system.variables, d.system.inputs, rng)
        ev = Evaluator(env=env, rng=rng)
        for sym in d.layout.params:
            ev(sym)
        lay = d.layout
        state, a = lay.state_from_env(ev.env), lay.accel_from_env(ev.env)
        u, p = lay.inputs_from_env(ev.env), lay.params_from_env(ev.env)
        out.append(Point(dict(ev.env), state, a, u, p))
    return out


def reference(expr, env):
    """``Evaluator`` value of *expr* at *env*, as a flat float array."""
    val = Evaluator(env=dict(env))(expr)
    return np.atleast_1d(np.asarray(val, dtype=float)).reshape(-1)


def assert_close(got, ref, what, rtol=RTOL):
    got = np.atleast_1d(np.asarray(got, dtype=float)).reshape(-1)
    ref = np.atleast_1d(np.asarray(ref, dtype=float)).reshape(-1)
    assert got.shape == ref.shape, f"{what}: shape {got.shape} vs {ref.shape}"
    scale = max(1.0, float(np.max(np.abs(ref)))) if ref.size else 1.0
    err = float(np.max(np.abs(got - ref))) if ref.size else 0.0
    assert err <= rtol * scale, f"{what}: generated {got} vs evaluate {ref} (err {err:.2e})"
