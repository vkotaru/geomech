"""Checks specific to the multibody systems of ``geomech.systems``: the spacecraft with
reaction wheels conserves its inertial angular momentum, and the default parameters and
scenes of both systems are what their docstrings state."""

from __future__ import annotations

import functools

import numpy as np

from geomech import to_standard_form
from geomech.codegen import load_module, module_source
from geomech.simulation import simulate
from geomech.systems import get


@functools.cache
def plant(name):
    s = get(name)
    eom = s.derive()
    sf = to_standard_form(eom, s.variables, list(s.inputs))
    return s, load_module(module_source(sf, s.layout(), lagrangian=s.L), f"mb_{name}")


def test_spacecraft_conserves_energy_and_inertial_angular_momentum():
    """Unforced: E and ``π = R (J_bΩ + Σ Iᵢ (θ̇ᵢ + bᵢ·Ω) bᵢ)`` are constant over 5 s, while
    the bus and the wheels exchange momentum (Ω and θ̇ change)."""
    s, module = plant("spacecraft_reaction_wheels")
    traj = simulate(module, s.x0, 5.0, s.params)
    E = traj.energy()
    assert np.max(np.abs(E - E[0])) / abs(E[0]) < 1e-9
    p = s.params
    Om = traj["Omega_R"]
    pi = np.einsum("nij,nj->ni", traj["R"], Om @ np.asarray(p["J_b"]).T)
    for i in (1, 2, 3):
        b = np.asarray(p[f"b{i}"])
        spin = traj[f"theta_{i}_dot"] + Om @ b
        pi += p[f"I_{i}"] * spin[:, None] * np.einsum("nij,j->ni", traj["R"], b)
    assert np.linalg.norm(pi[0]) > 0.1
    assert float(np.max(np.abs(pi - pi[0]))) < 1e-10
    assert float(np.max(np.abs(Om - Om[0]))) > 1e-3  # the motion is not trivial


def test_spacecraft_wheel_torque_spins_the_bus_the_other_way():
    """A positive torque on wheel 3 from rest spins the wheel up (+b₃) and the bus down
    (−b₃): the total angular momentum stays zero."""
    s, module = plant("spacecraft_reaction_wheels")
    x0 = s.layout().pack(R=np.eye(3), theta_1=0.0, theta_2=0.0, theta_3=0.0)
    traj = simulate(module, x0, 0.5, s.params, u=np.array([0.0, 0.0, 0.01]))
    assert traj["theta_3_dot"][-1] > 0 and traj["Omega_R"][-1, 2] < 0


def test_spacecraft_and_biped_defaults_are_physical():
    sc = get("spacecraft_reaction_wheels").params
    B = np.array([sc[f"b{i}"] for i in (1, 2, 3)])
    assert np.allclose(B @ B.T, np.eye(3))  # orthonormal wheel axes
    assert all(sc[f"I_{i}"] > 0 for i in (1, 2, 3))
    Jp = np.asarray(sc["J_b"]) + sum(sc[f"I_{i}"] * np.outer(b, b) for i, b in enumerate(B, 1))
    assert np.all(np.linalg.eigvalsh(Jp) > 0)
    bp = get("biped_3d").params
    assert 0 < bp["a"] < bp["l"] and bp["m_H"] > 0 and bp["m"] > 0


def test_scenes_draw_solid_bodies():
    from geomech.viz import Box, Cylinder

    s = get("spacecraft_reaction_wheels")
    shapes = [type(b.shape) for b in s.scene().bodies]
    assert shapes.count(Box) == 3 and shapes.count(Cylinder) == 3  # bus, 2 panels, 3 wheels
    s = get("biped_3d")
    sc = s.scene()
    shapes = [type(b.shape) for b in sc.bodies]
    assert shapes.count(Cylinder) == 2 and shapes.count(Box) == 3  # legs; feet and ground
    assert len(sc.points) == 1  # the hip


def test_api_reference_documents_the_multibody_module():
    from pathlib import Path

    page = (Path(__file__).resolve().parents[1] / "docs" / "api" / "systems.rst").read_text()
    assert "geomech.systems.multibody" in page and "spacecraft_reaction_wheels" in page


def test_rigid_chain_with_two_links_is_the_double_rigid_pendulum():
    """``n_link_rigid_chain(n_links=2)`` has the Lagrangian of ``double_rigid_pendulum``
    term for term, and its virtual work and equations agree numerically when the body
    moments of the double pendulum are the joint torques' net moments
    ``M₁ = R₁ᵀ(τ₁ − τ₂)``, ``M₂ = R₂ᵀτ₂``."""
    from geomech.numeric import Evaluator, sample_state

    chain, double = get("n_link_rigid_chain", n_links=2), get("double_rigid_pendulum")
    assert str(chain.L) == str(double.L) and repr(chain.L) == repr(double.L)
    eom_c, eom_d = chain.derive(), double.derive()
    assert list(eom_c) == list(eom_d)
    R1, R2 = chain.symbols["R"]
    for seed in range(3):
        rng = np.random.default_rng(seed)
        env = sample_state(chain.variables, list(chain.inputs), rng)
        ev = Evaluator(env=env, rng=rng)
        tau1, tau2 = ev(chain.symbols["tau"][0]), ev(chain.symbols["tau"][1])
        env["M1"] = ev(R1).T @ (tau1 - tau2)
        env["M2"] = ev(R2).T @ tau2
        for key in eom_c:
            got, ref = ev(eom_c[key][1]), ev(eom_d[key][1])
            assert np.allclose(got, ref, rtol=1e-12, atol=1e-12), key
        for k, eta in (("eta1", R1.get_variation_vector()), ("eta2", R2.get_variation_vector())):
            env[str(eta)] = rng.normal(size=(3, 1))
        assert np.isclose(ev(chain.dW), ev(double.dW), rtol=1e-12, atol=1e-12)


def test_rigid_chain_defaults_are_uniform_boxes():
    """Each default link is a uniform box drawn with its own size: ρₖ = lₖ / 2 along −e₃
    and the inertia is the box's about its centre of mass."""
    from geomech.systems.rigid import LINK_WIDTH

    s = get("n_link_rigid_chain", n_links=4)
    p = s.params
    wx, wy = LINK_WIDTH
    for k in range(1, 5):
        rho = np.asarray(p[f"\\rho_{k}"])
        L = -2 * rho[2]
        assert np.allclose(rho[:2], 0) and L > 0
        if k < 4:
            assert np.allclose(p[f"l_{k}"], 2 * rho)
        m = p[f"m{k}"]
        J = np.diag([m * (wy**2 + L**2), m * (wx**2 + L**2), m * (wx**2 + wy**2)]) / 12
        assert np.allclose(p[f"J{k}"], J)
    assert "l_4" not in p
    lengths = [b.shape.lz for b in s.scene().bodies]
    assert np.allclose(lengths, [-2 * p[f"\\rho_{k}"][2] for k in range(1, 5)])


def test_scenes_follow_the_parameters_they_can():
    """The wheels sit at 0.3 bᵢ (the default bᵢ = eᵢ is also their disc axis) and the
    biped's legs are drawn with the default leg length."""
    from geomech.viz import Cylinder

    s = get("spacecraft_reaction_wheels")
    sc = s.scene()
    wheels = [b for b in sc.bodies if isinstance(b.shape, Cylinder)]
    for i, (w, axis) in enumerate(zip(wheels, "xyz"), 1):
        assert w.shape.axis == axis
        assert np.allclose(sc.body_offset(w, dict(s.params)), 0.3 * s.params[f"b{i}"])
    b = get("biped_3d")
    legs = [x for x in b.scene().bodies if isinstance(x.shape, Cylinder)]
    assert [leg.shape.length for leg in legs] == [b.params["l"]] * 2
