"""``geomech.simulation.reduced_linearization``: the matrices ``lqr`` designs with."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy")

from geomech.simulation import ReducedLinearization, lqr, model_of, reduced_linearization
from geomech.simulation.control import tangent_basis
from tests.simulation_systems import E3, SETUPS, linearized_module, standard_module

#: systems with an equilibrium at the configuration of their SETUPS state, zero velocity
#: and zero input (no gravity, or hanging), and one with S²
CASES = ["free_rigid_body", "spring_mass_r", "spherical_pendulum_s2"]


def _equilibrium(name):
    p, x0 = SETUPS[name]
    mod = standard_module(name)
    model = model_of(mod)
    x = np.array(x0, dtype=float)
    for var in model.variables:
        x[var.vel] = 0.0
        if var.manifold == "S2":
            x[var.cfg] = -E3  # hanging
        if var.manifold in ("R", "R3"):
            x[var.cfg] = 0.0
    return p, x, np.zeros(model.n_inputs)


@pytest.mark.parametrize("name", CASES)
def test_matches_lqr(name):
    p, x, u = _equilibrium(name)
    lin = linearized_module(name)
    red = reduced_linearization(lin, x, u, p)
    assert isinstance(red, ReducedLinearization)
    ctrl = lqr(lin, x, u, p, Q=1.0, R=1.0, module=standard_module(name))
    assert np.array_equal(red.A, ctrl.A) and np.array_equal(red.B, ctrl.B)
    assert np.array_equal(red.T, ctrl.T) and red.labels == ctrl.labels
    assert red.bases.keys() == ctrl.bases.keys()
    for k in red.bases:
        assert np.array_equal(red.bases[k], ctrl.bases[k])
    n = red.T.shape[1]
    assert red.A.shape == (n, n) and red.B.shape == (n, model_of(lin).n_inputs)
    assert len(red.labels) == n
    Mr = red.T.T @ red.M_full @ red.T
    assert np.allclose(Mr @ red.A, red.T.T @ red.A_full @ red.T)


def test_labels_and_s2_reduction():
    p, x, u = _equilibrium("spherical_pendulum_s2")
    red = reduced_linearization(linearized_module("spherical_pendulum_s2"), x, u, p)
    assert red.labels == ["xi_q[0]", "xi_q[1]", "delta_omega_q[0]", "delta_omega_q[1]"]
    B = red.bases["q"]
    assert np.allclose(B, tangent_basis(-E3)) and np.allclose(B.T @ B, np.eye(2))
    # the hanging spherical pendulum: ẍ_tangent = −(g/l) ξ, so A = [[0, I], [−g/l I, 0]]
    gl = p["g"] / p["l"]
    assert np.allclose(
        red.A, np.block([[np.zeros((2, 2)), np.eye(2)], [-gl * np.eye(2), np.zeros((2, 2))]])
    )


def test_custom_basis():
    p, x, u = _equilibrium("spherical_pendulum_s2")
    lin = linearized_module("spherical_pendulum_s2")

    def rotated(q):
        B = tangent_basis(q)
        c, s = np.cos(0.3), np.sin(0.3)
        return B @ np.array([[c, -s], [s, c]])

    a = reduced_linearization(lin, x, u, p)
    b = reduced_linearization(lin, x, u, p, basis=rotated)
    assert np.allclose(b.bases["q"], rotated(x[0:3]))
    # the same system in rotated coordinates: z_b = Gᵀ z_a with G the 2×2 rotation
    c, s_ = np.cos(0.3), np.sin(0.3)
    G = np.array([[c, -s_], [s_, c]])
    Gf = np.kron(np.eye(2), G)
    assert np.allclose(b.A, Gf.T @ a.A @ Gf) and np.allclose(b.B, Gf.T @ a.B)


def test_rejects_a_standard_form_module():
    p, x, u = _equilibrium("free_rigid_body")
    with pytest.raises(TypeError, match="M_lin"):
        reduced_linearization(standard_module("free_rigid_body"), x, u, p)
