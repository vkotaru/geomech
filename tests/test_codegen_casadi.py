"""The CasADi target equals the NumPy target (spec verification item 2).

Skipped when ``casadi`` is not installed (``pip install geomech[casadi]``); CI installs it
in one job. Every block function of the standard-form, ``compute_eom`` and linearized
modules of every system is compared on the 20 states of the NumPy check.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.codegen_systems import ALL_SYSTEMS, derived, modules, sample_points

ca = pytest.importorskip(
    "casadi", reason="the CasADi target needs the optional dependency: pip install casadi"
)

_ids = dict(ids=lambda f: f.__name__)


def _np(x):
    return np.array(ca.DM(x), dtype=float) if isinstance(x, (ca.DM, ca.SX)) else np.asarray(x)


def _same(got, ref, what):
    got = _np(got).reshape(np.shape(ref))
    scale = max(1.0, float(np.max(np.abs(ref)))) if np.size(ref) else 1.0
    err = float(np.max(np.abs(got - ref))) if np.size(ref) else 0.0
    assert err <= 1e-12 * scale, f"{what}: casadi {got} vs numpy {ref}"


@pytest.mark.parametrize("build", ALL_SYSTEMS, **_ids)
def test_casadi_equals_numpy(build):
    d = derived(build)
    np_mods, ca_mods = modules(build, "numpy"), modules(build, "casadi")
    for pt in sample_points(d):
        calls = {
            "mass": (pt.state, pt.p),
            "bias": (pt.state, pt.p),
            "input_matrix": (pt.state, pt.p),
            "kinematics": (pt.state,),
            "energy": (pt.state, pt.p),
            "residual": (pt.state, pt.a, pt.u, pt.p),
            "M_lin": (pt.state, pt.a, pt.u, pt.p),
            "A_lin": (pt.state, pt.a, pt.u, pt.p),
            "B_lin": (pt.state, pt.a, pt.u, pt.p),
        }
        for np_mod, ca_mod in zip(np_mods, ca_mods):
            if np_mod is None:
                continue
            for name, args in calls.items():
                if hasattr(np_mod, name):
                    ref = getattr(np_mod, name)(*args)
                    _same(getattr(ca_mod, name)(*args), ref, f"{build.__name__}.{name}")


def test_symbolic_arguments_and_jacobian():
    """The CasADi module accepts SX symbols; its Jacobian matches a finite difference."""
    from tests.vbl_systems import quadrotor_with_load

    d = derived(quadrotor_with_load)
    np_mod, ca_mod = (
        modules(quadrotor_with_load, "numpy")[0],
        modules(quadrotor_with_load, "casadi")[0],
    )
    pt = sample_points(d, n=1)[0]
    x = ca.SX.sym("x", d.layout.n_state)
    bias = ca.Function("bias", [x], [ca_mod.bias(x, pt.p)])
    jac = ca.Function("dbias", [x], [ca.jacobian(ca_mod.bias(x, pt.p), x)])
    np.testing.assert_allclose(_np(bias(pt.state)).ravel(), np_mod.bias(pt.state, pt.p))
    h, J = 1e-6, _np(jac(pt.state))
    for i in range(d.layout.n_state):
        e = np.zeros(d.layout.n_state)
        e[i] = h
        fd = (np_mod.bias(pt.state + e, pt.p) - np_mod.bias(pt.state - e, pt.p)) / (2 * h)
        np.testing.assert_allclose(J[:, i], fd, atol=1e-6)
    p_sym = {k: ca.SX.sym(k.replace("\\", ""), *np.shape(v) or (1,)) for k, v in pt.p.items()}
    assert isinstance(ca_mod.mass(x, p_sym), ca.SX)
