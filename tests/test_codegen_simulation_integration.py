"""Code generation and simulation together: the final review of the branch.

- R1: a linearization generated with another input or variable order than its
  ``StateLayout`` raises at ``generate_module`` (before it was accepted, and LQR on it
  diverged); ``lqr`` matches ``LINEARIZATION`` to ``LAYOUT`` by symbol and checks that the
  plant module has the same layout.
- The reviewer's user story from the public API only: the quadrotor with a
  cable-suspended load (R³ × S² × SO(3)), LQR about hover, the closed loop from 10° off.
- R2/R3: ``collect_params`` leaves out tangent and variation vectors, accepts the three
  system types and raises ``TypeError`` for what it cannot walk; ``StateLayout`` rejects a
  variation vector as a parameter.
- ``StateLayout.pack``/``unpack``.
- R5: ``examples/inverted_rigid_pendulum_lqr.py`` runs under ``-W error``.
"""

from __future__ import annotations

import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

from geomech import (
    S2,
    SO3,
    Cross,
    Dot,
    Matrix,
    Scalar,
    SystemVariables,
    Vector,
    compute_eom,
    getScalars,
    linearize,
    to_standard_form,
)
from geomech.codegen import StateLayout, collect_params, generate_module, lambdify, load_module
from geomech.numeric import expm_so3

ROOT = Path(__file__).resolve().parent.parent
HALF = Scalar("0.5", value=0.5, attr=["Constant"])
E3 = np.array([0.0, 0.0, 1.0])
P = {
    "m_Q": 1.5,
    "m_L": 0.3,
    "g": 9.81,
    "l": 1.0,
    "J": np.diag([0.02, 0.02, 0.04]),
    "e3": E3,
}
PARAM_NAMES = ["J", "e3", "g", "l", "m_L", "m_Q"]


class Quadrotor:
    """The quadrotor with a cable-suspended load as a user writes it from the docs.

    ``x_Q = x_L − l q``; ``L = ½ m_Q|ẋ_Q|² + ½ m_L|ẋ_L|² + ½ Ω·JΩ − m_Q g e₃·x_Q −
    m_L g e₃·x_L``; ``δW = δx_Q · (R e₃ f) + η · M``.
    """

    def __init__(self):
        mQ, mL, g, l = getScalars("m_Q m_L g l", attr=["Constant"])
        J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
        e3 = Vector("e3", attr=["Constant"])
        self.xL, self.q, self.R = Vector("x_L"), S2("q"), SO3("R")
        self.f, self.M = Scalar("f"), Vector("M")
        Om, eta = self.R.get_tangent_vector(), self.R.get_variation_vector()
        xQ = self.xL - self.q * l
        vQ, vL = xQ.t_diff(), self.xL.t_diff()
        self.L = (
            HALF * mQ * Dot(vQ, vQ)
            + HALF * mL * Dot(vL, vL)
            + HALF * Dot(Om, J * Om)
            - mQ * g * Dot(e3, xQ)
            - mL * g * Dot(e3, self.xL)
        )
        self.dW = Dot(xQ.delta(), (self.R * e3) * self.f) + Dot(eta, self.M)
        self.V = SystemVariables(vectors=[self.xL, self.q], matrices=[self.R])
        self.inputs = [self.f, self.M]
        self.eom = compute_eom(self.L, self.dW, self.V)
        self.sf = to_standard_form(self.eom, self.V, inputs=self.inputs)
        self.lin = linearize(self.eom, self.V, inputs=self.inputs)
        self.layout = StateLayout(
            self.V, self.inputs, collect_params([self.L, self.dW], self.V, self.inputs)
        )
        self.plant = load_module(generate_module(self.sf, self.layout, lagrangian=self.L))
        self.lin_module = load_module(generate_module(self.lin, self.layout))
        self.x_eq = self.layout.pack(x_L=np.zeros(3), q=-E3, R=np.eye(3))  # hover, at rest
        self.u_eq = np.r_[(P["m_Q"] + P["m_L"]) * P["g"], 0.0, 0.0, 0.0]


@pytest.fixture(scope="module")
def quad():
    return Quadrotor()


# ---------------------------------------------------------------------------
# R1: one layout for both modules
# ---------------------------------------------------------------------------


def test_swapped_inputs_raise_at_generate_module(quad):
    """The reviewer's case: ``StateLayout(V, [f, M])`` with ``linearize(eom, V, [M, f])``.
    Before the fix the module was generated, ``lqr`` reported max Re λ = −0.60 and the
    closed loop diverged (``|x_L| ≈ 35 m`` after 3 s); now nothing is generated."""
    lin = linearize(quad.eom, quad.V, inputs=[quad.M, quad.f])
    msg = (
        "generate_module: the linearization's inputs ['\\delta{M}', '\\delta{f}'] are not"
        " the variations of the layout's inputs ['f', 'M'] in that order"
    )
    with pytest.raises(ValueError, match=re.escape(msg)) as info:
        generate_module(lin, quad.layout)
    assert "generate both modules from that one StateLayout" in str(info.value)


def test_swapped_variables_raise_at_generate_module(quad):
    """``linearize`` with ``SystemVariables(vectors=[q, x_L])`` against a layout of
    ``[x_L, q]``: before the fix ``lqr`` failed with a misleading singular ``M_lin``."""
    V2 = SystemVariables(vectors=[quad.q, quad.xL], matrices=[quad.R])
    lin = linearize(quad.eom, V2, inputs=quad.inputs)
    msg = (
        "the linearization's configuration variables ['q', 'x_L', 'R'] (rows ['kin:q',"
        " 'kin:x_L', 'kin:R', '\\delta{x_L}', '\\xi_{q}', '\\eta_{R}']) are not the"
        " layout's variables ['x_L', 'q', 'R'] in that order"
    )
    with pytest.raises(ValueError, match=re.escape(msg)):
        generate_module(lin, quad.layout)


def test_same_order_generates(quad):
    data = quad.lin_module.LINEARIZATION
    assert [r["label"] for r in data["rows"]][:3] == ["kin:x_L", "kin:q", "kin:R"]
    assert [u["symbol"] for u in data["inputs"]] == ["\\delta{f}", "\\delta{M}"]


def test_lqr_rejects_a_plant_with_another_input_order(quad):
    pytest.importorskip("scipy")
    from geomech.simulation import lqr

    layout2 = StateLayout(quad.V, [quad.M, quad.f], quad.layout.params)
    sf2 = to_standard_form(quad.eom, quad.V, inputs=[quad.M, quad.f])
    plant2 = load_module(generate_module(sf2, layout2))
    msg = (
        "lqr: module does not have the layout of lin_module: lin_module has the variables"
        " ['x_L', 'q', 'R'] and inputs ['f', 'M'], module has the variables"
        " ['x_L', 'q', 'R'] and inputs ['M', 'f']"
    )
    with pytest.raises(ValueError, match=re.escape(msg)):
        lqr(quad.lin_module, quad.x_eq, quad.u_eq, P, 1.0, 1.0, module=plant2)


def _tampered(module, key, transform):
    """A copy of the generated *module* whose ``LINEARIZATION[key]`` is transformed."""
    import copy
    import types

    clone = types.ModuleType("tampered")
    clone.__dict__.update(module.__dict__)
    clone.LINEARIZATION = copy.deepcopy(module.LINEARIZATION)
    clone.LINEARIZATION[key] = transform(clone.LINEARIZATION[key])
    return clone


def test_lqr_matches_linearization_to_layout_by_symbol(quad):
    """Defence in depth for a module generated before the check (or edited): the reduction
    compares the symbols of ``LINEARIZATION`` with ``LAYOUT``, not only the positions."""
    pytest.importorskip("scipy")
    from geomech.simulation import lqr

    swapped = _tampered(quad.lin_module, "inputs", lambda e: e[::-1])
    with pytest.raises(ValueError, match=re.escape("inputs ['\\delta{M}', '\\delta{f}']")):
        lqr(swapped, quad.x_eq, quad.u_eq, P, 1.0, 1.0, module=quad.plant)

    def relabel(rows):
        rows[0]["label"], rows[1]["label"] = rows[1]["label"], rows[0]["label"]
        return rows

    relabelled = _tampered(quad.lin_module, "rows", relabel)
    msg = "LINEARIZATION of lin_module has the variables ['q', 'x_L', 'R']"
    with pytest.raises(ValueError, match=re.escape(msg)):
        lqr(relabelled, quad.x_eq, quad.u_eq, P, 1.0, 1.0, module=quad.plant)


# ---------------------------------------------------------------------------
# The user story: LQR about hover, closed loop from 10° off
# ---------------------------------------------------------------------------


def test_user_story_quadrotor_with_load_lqr(quad):
    pytest.importorskip("scipy")
    from geomech.simulation import lqr, simulate

    assert quad.layout.param_names == PARAM_NAMES
    ctrl = lqr(
        quad.lin_module,
        quad.x_eq,
        quad.u_eq,
        P,
        Q=np.eye(16),
        R=np.eye(4),
        module=quad.plant,
    )
    assert ctrl.K.shape == (4, 16)
    assert len(ctrl.labels) == 16  # 3 + 2 + 3 + 3 + 2 + 3
    assert np.max(ctrl.closed_loop_eigenvalues.real) == pytest.approx(-0.596, abs=0.01)

    q0 = expm_so3(np.deg2rad(10.0) * np.array([1.0, 0.0, 0.0])) @ -E3
    x0 = quad.layout.pack(x_L=np.zeros(3), q=q0, R=np.eye(3))
    traj = simulate(quad.plant, x0, 4.0, P, u=ctrl)
    end = quad.layout.unpack(traj.state[-1])
    angle0 = np.degrees(np.arccos(np.clip(q0 @ -E3, -1, 1)))
    angle = np.degrees(np.arccos(np.clip(end["q"] @ -E3, -1, 1)))
    e0, e1 = np.abs(ctrl.error(x0)).max(), np.abs(ctrl.error(traj.state[-1])).max()
    # slowest closed-loop mode −0.60/s: e^(−2.4) ≈ 0.09 over 4 s; the swapped-input
    # controller of R1 reached |x_L| = 35 m after 3 s
    assert e1 < 0.2 * e0, (e0, e1)
    assert angle < 0.2 * angle0, (angle0, angle)
    assert np.linalg.norm(end["x_L"]) < 0.1
    assert traj.max_drift() < 1e-10


def test_user_story_unforced_energy(quad):
    from geomech.simulation import simulate

    q0 = expm_so3(np.array([0.3, -0.2, 0.0])) @ -E3
    x0 = quad.layout.pack(
        x_L=[0.1, 0.0, 0.2],
        q=q0,
        R=expm_so3(np.array([0.1, 0.2, -0.1])),
        x_L_dot=[0.2, -0.1, 0.3],
        omega_q=np.cross(q0, [0.5, 0.2, 0.1]),
        Omega_R=[0.3, -0.4, 0.5],
    )
    traj = simulate(quad.plant, x0, 1.0, P, output_every=100)
    E = traj.energy()
    assert np.max(np.abs(E - E[0])) < 1e-9 * max(1.0, abs(E[0]))
    assert set(traj.constraint_drift) == {"|q|-1", "q.omega_q", "R^T R-I"}


# ---------------------------------------------------------------------------
# R2/R3: collect_params and parameters
# ---------------------------------------------------------------------------


def test_collect_params_leaves_out_variation_vectors(quad):
    """``[L, dW]``, the natural call: ``ξ`` and ``η`` of ``δW`` are not parameters."""
    params = collect_params([quad.L, quad.dW], quad.V, quad.inputs)
    assert [str(s) for s in params] == PARAM_NAMES


def test_collect_params_accepts_the_system_types(quad):
    expected = set(PARAM_NAMES)

    def names(x):
        return {str(s) for s in collect_params(x, quad.V, quad.inputs)}

    assert names(quad.eom) == expected  # the variation vector keys are left out
    assert names(quad.eom.values()) == expected
    assert names(quad.sf) == expected
    assert names(quad.sf.values()) == expected  # gave [] before the fix
    assert names(quad.lin) == expected
    assert names(quad.L) == expected
    assert names([quad.L, quad.eom, quad.lin]) == expected
    assert names([]) == set()


@pytest.mark.parametrize(
    "bad, what",
    [(3, "int"), ("L", "str"), (None, "NoneType"), (object(), "object"), ([1.0], "float")],
)
def test_collect_params_raises_for_what_it_cannot_walk(quad, bad, what):
    with pytest.raises(TypeError, match=f"collect_params: cannot find symbols in {what}"):
        collect_params(bad, quad.V, quad.inputs)
    with pytest.raises(TypeError, match="pass the expressions of the system"):
        collect_params([quad.L, bad], quad.V, quad.inputs)


@pytest.mark.parametrize("which", ["xi", "eta"])
def test_layout_rejects_a_variation_vector_as_param(quad, which):
    sym = quad.q.get_variation_vector() if which == "xi" else quad.R.get_variation_vector()
    msg = f"StateLayout: parameter '{sym}' is a tangent or variation vector of a configuration"
    with pytest.raises(ValueError, match=re.escape(msg)):
        StateLayout(quad.V, quad.inputs, [*quad.layout.params, sym])


def test_stray_variation_vector_in_an_expression_is_named(quad):
    xi = quad.q.get_variation_vector()
    with pytest.raises(ValueError, match=r"a tangent or variation vector that is not a velocity"):
        lambdify(Dot(xi, Cross(quad.q, xi)), quad.layout)


# ---------------------------------------------------------------------------
# StateLayout.pack / unpack
# ---------------------------------------------------------------------------


def test_pack_matches_the_layout_slices(quad):
    lay = quad.layout
    R = expm_so3(np.array([0.1, 0.2, 0.3]))
    vals = dict(
        x_L=[1.0, 2.0, 3.0],
        q=-E3,
        R=R,
        x_L_dot=[4.0, 5.0, 6.0],
        omega_q=[0.1, 0.2, 0.0],
        Omega_R=[7.0, 8.0, 9.0],
    )
    x = lay.pack(**vals)
    s = lay.slices
    np.testing.assert_array_equal(x[s["x_L"]], vals["x_L"])
    np.testing.assert_array_equal(x[s["q"]], -E3)
    np.testing.assert_array_equal(x[s["R"]], R.reshape(9))  # row-major
    np.testing.assert_array_equal(x[s["\\omega_{q}"]], vals["omega_q"])
    np.testing.assert_array_equal(x[s["\\Omega_{R}"]], vals["Omega_R"])
    # printed names as dict keys, R as 9 entries
    by_key = {"x_L": vals["x_L"], "q": -E3, "R": R.reshape(9), "\\Omega_{R}": vals["Omega_R"]}
    y = lay.pack(by_key, x_L_dot=vals["x_L_dot"], omega_q=vals["omega_q"])
    np.testing.assert_array_equal(x, y)
    out = lay.unpack(x)
    assert set(out) == set(vals)
    for k, v in vals.items():
        np.testing.assert_array_equal(out[k], np.asarray(v, float))
    assert out["R"].shape == (3, 3)


def test_pack_defaults_velocities_to_zero(quad):
    x = quad.layout.pack(x_L=np.zeros(3), q=-E3, R=np.eye(3))
    np.testing.assert_array_equal(x, np.r_[np.zeros(3), -E3, np.eye(3).reshape(9), np.zeros(9)])


def test_unpack_a_trajectory_array(quad):
    X = np.stack([quad.x_eq, quad.x_eq])
    out = quad.layout.unpack(X)
    assert out["q"].shape == (2, 3) and out["R"].shape == (2, 3, 3)
    with pytest.raises(ValueError, match="StateLayout.unpack: state must have 24 entries"):
        quad.layout.unpack(np.zeros(5))


def test_pack_errors(quad):
    lay = quad.layout
    with pytest.raises(KeyError, match=r"no value for the configuration slot\(s\) \['R'\]"):
        lay.pack(x_L=np.zeros(3), q=-E3)
    with pytest.raises(KeyError, match="'omega' is not a configuration or velocity slot"):
        lay.pack(x_L=np.zeros(3), q=-E3, R=np.eye(3), omega=np.zeros(3))
    with pytest.raises(ValueError, match="'q' needs 3 entries"):
        lay.pack(x_L=np.zeros(3), q=[0.0, 1.0], R=np.eye(3))
    with pytest.raises(ValueError, match="'q' given twice"):
        lay.pack({"q": -E3}, x_L=np.zeros(3), q=-E3, R=np.eye(3))


def test_pack_scalar_slots():
    s = Scalar("s")
    lay = StateLayout(SystemVariables(scalars=[s]))
    np.testing.assert_array_equal(lay.pack(s=0.5, s_dot=-1.0), [0.5, -1.0])
    assert lay.unpack([0.5, -1.0]) == {"s": 0.5, "s_dot": -1.0}


# ---------------------------------------------------------------------------
# R5: the LQR example runs
# ---------------------------------------------------------------------------


def test_lqr_example_runs_under_warnings_as_errors():
    pytest.importorskip("scipy")
    script = ROOT / "examples" / "inverted_rigid_pendulum_lqr.py"
    start = time.perf_counter()
    run = subprocess.run(
        [sys.executable, "-W", "error", str(script)],
        capture_output=True,
        text=True,
        timeout=60,
        cwd=ROOT,
    )
    elapsed = time.perf_counter() - start
    assert run.returncode == 0, run.stderr
    assert elapsed < 10.0, elapsed
    final = re.search(r"t = 5\.0 s\s+\|eta\| = ([0-9.e+-]+) rad", run.stdout)
    assert final is not None, run.stdout
    assert float(final.group(1)) < 1e-6  # docs/simulation.md: 4.846e-08 rad
    assert "max |RᵀR − I| over the run" in run.stdout
