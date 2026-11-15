"""Generated NumPy code equals ``geomech.numeric.evaluate`` of the symbolic blocks.

Spec ``specs/codegen.md``, verification item 1: for every system of
``tests/vbl_systems.py``, ``tests/phase1_systems.py`` and the numerical check gate, each block of
the generated module (standard form ``M``/``f``/``G``, the ``compute_eom`` equations, the
linearization ``M``/``A``/``B`` and the kinematics) is compared with the ``Evaluator`` value
of the symbolic block at 20 random admissible states (``sample_state``), to 1e-12 relative
(``|generated − evaluate| ≤ 1e-12 · max(1, |evaluate|)``).
"""

from __future__ import annotations

import numpy as np
import pytest

from geomech.codegen import lambdify
from tests.codegen_systems import (
    ALL_SYSTEMS,
    VBL_SYSTEMS,
    assert_close,
    derived,
    modules,
    reference,
    sample_points,
)

_ids = dict(ids=lambda f: f.__name__)


def _spans(sizes):
    out, pos = [], 0
    for n in sizes:
        out.append(slice(pos, pos + n))
        pos += n
    return out


def _eq_keys(d):
    """The equation keys in layout order with their row slices."""
    keys = []
    for var in d.layout.config_variables:
        vv = var.get_variation_vector() if hasattr(var, "manifold_info") else var.delta()
        keys.append(str(vv))
    rows = _spans(1 if d.eom[k][1].type.name == "SCALAR" else 3 for k in keys)
    return list(zip(keys, rows))


def test_system_count():
    assert len(VBL_SYSTEMS) == 9
    assert set(VBL_SYSTEMS) <= set(ALL_SYSTEMS)
    assert len(ALL_SYSTEMS) == 64


@pytest.mark.parametrize("build", ALL_SYSTEMS, **_ids)
def test_standard_form_blocks(build):
    d = derived(build)
    mod, _, _ = modules(build)
    acc = {s.key: s.slice for s in d.layout.acceleration}
    inp = {s.key: s.slice for s in d.layout.input_slots}
    for pt in sample_points(d):
        M, f = mod.mass(pt.state, pt.p), mod.bias(pt.state, pt.p)
        G = mod.input_matrix(pt.state, pt.p)
        assert M.shape == (len(f), d.layout.n_accel) and G.shape == (len(f), d.layout.n_inputs)
        for key, rows in _eq_keys(d):
            eq = d.sf[key]
            assert_close(f[rows], reference(eq.f, pt.env), f"f[{key}]")
            for col, block in eq.M.items():
                assert_close(M[rows, acc[col]], reference(block, pt.env), f"M[{key}][{col}]")
            for col, block in eq.G.items():
                assert_close(G[rows, inp[col]], reference(block, pt.env), f"G[{key}][{col}]")
            missing = set(acc) - set(eq.M)
            for col in missing:
                assert not np.any(M[rows, acc[col]]), f"M[{key}][{col}] should be zero"


@pytest.mark.parametrize("build", ALL_SYSTEMS, **_ids)
def test_eom_residual_and_kinematics(build):
    d = derived(build)
    _, mod, _ = modules(build)
    for pt in sample_points(d):
        E = mod.residual(pt.state, pt.a, pt.u, pt.p)
        for key, rows in _eq_keys(d):
            assert_close(E[rows], reference(d.eom[key][1], pt.env), f"E[{key}]")
        kin = mod.kinematics(pt.state)
        for var, slot in zip(d.layout.config_variables, d.layout.configuration):
            assert_close(kin[slot.slice], reference(var.t_diff(), pt.env), f"d/dt {var}")


@pytest.mark.parametrize("build", [b for b in ALL_SYSTEMS if derived(b).lin is not None], **_ids)
def test_linearized_blocks(build):
    d = derived(build)
    _, _, mod = modules(build)
    lin = d.lin
    rows = dict(zip(lin.rows, _spans(lin.block_shape(r, lin.states[0])[0] for r in lin.rows)))

    def cols(symbols):
        return dict(
            zip(map(str, symbols), _spans(lin.block_shape(lin.rows[0], s)[1] for s in symbols))
        )

    for pt in sample_points(d):
        args = (pt.state, pt.a, pt.u, pt.p)
        for name, blocks, symbols in (
            ("M_lin", lin.M, lin.dstates),
            ("A_lin", lin.A, lin.states),
            ("B_lin", lin.B, lin.inputs),
        ):
            got, spans = getattr(mod, name)(*args), cols(symbols)
            for row, row_blocks in blocks.items():
                for col, block in row_blocks.items():
                    ref = reference(block, pt.env)
                    assert_close(got[rows[row], spans[col]], ref, f"{name}[{row}][{col}]")


def test_lambdify_matches_module():
    d = derived(VBL_SYSTEMS[-1])
    pt = sample_points(d, n=1)[0]
    for key, (_, eq) in d.eom.items():
        fn = lambdify(eq, d.layout)
        assert "def f(state, u, p, a=None):" in fn.source
        assert_close(fn(pt.state, pt.u, pt.p, pt.a), reference(eq, pt.env), key)


def test_lambdify_shapes():
    d = derived(VBL_SYSTEMS[3])  # rigid pendulum about the pivot
    pt = sample_points(d, n=1)[0]
    J = d.system.symbols["J"]
    R = d.system.symbols["R"]
    assert np.shape(lambdify(J, d.layout)(pt.state, pt.u, pt.p)) == (3, 3)
    assert np.shape(lambdify(R * d.system.symbols["rho"], d.layout)(pt.state, pt.u, pt.p)) == (3,)
    R_val = lambdify(R, d.layout)(pt.state, pt.u, pt.p)
    np.testing.assert_allclose(R_val, pt.env[str(R)])
    np.testing.assert_allclose(pt.state[0:9], np.asarray(pt.env[str(R)]).reshape(9))
