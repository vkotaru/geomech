"""``geomech.viz``: scene geometry against NumPy, the LaTeX name map, the HTML report and
a headless-browser smoke test (opt-in, see ``tests/browser.py``)."""

from __future__ import annotations

import html.parser
import importlib.util
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

from geomech import SO3, Dot, Matrix, Scalar, Vector, to_latex
from geomech.core.base.expressions import ZeroVector
from geomech.simulation import simulate
from geomech.viz import (
    Box,
    Cylinder,
    PointMass,
    Report,
    RigidBody,
    Rod,
    Scene,
    Sphere,
    Trace,
    latex_names,
    pretty_name,
)
from geomech.viz.report import EXTERNAL_URLS
from tests import browser
from tests.codegen_systems import derived
from tests.simulation_systems import BUILDERS, SETUPS, standard_module

ROOT = Path(__file__).resolve().parents[1]
TOL = 1e-12


def _load_example(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "examples" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _short_traj(name, t=0.5):
    p, x0 = SETUPS[name]
    return simulate(standard_module(name), x0, t, p, dt=1e-3, output_every=10)


# ---------------------------------------------------------------------------
# LaTeX names
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, point, expected",
    [
        ("J1", None, "J_1"),
        ("e3", None, "e_3"),
        ("m12", None, "m_{12}"),
        ("l0", None, "l_0"),
        ("\\rho1", None, "\\rho_1"),
        ("m_L", None, "m_L"),
        ("\\rho_1", None, "\\rho_1"),
        ("q", None, "q"),
        ("0.5", None, "0.5"),
        ("(2)", None, "(2)"),
        ("x1y", None, "x1y"),
        ("\\Omega_{R1}", "R1", "\\Omega_1"),
        ("\\eta_{R2}", "R2", "\\eta_2"),
        ("\\omega_{q}", "q", "\\omega_{q}"),
        ("\\Omega_{R1}", None, "\\Omega_{R1}"),
    ],
)
def test_pretty_name(name, point, expected):
    assert pretty_name(name, point) == expected


def test_latex_names_of_double_pendulum():
    d = derived(BUILDERS["double_rigid_pendulum"])
    names = latex_names(d.system.L, d.eom, d.sf)
    assert names["J1"] == "J_1" and names["m2"] == "m_2" and names["e3"] == "e_3"
    assert names["\\Omega_{R1}"] == "\\Omega_1" and names["\\eta_{R2}"] == "\\eta_2"
    assert names["R1"] == "R_1"
    assert "\\rho_1" not in names and "l_1" not in names
    tex = to_latex(d.eom["\\eta_{R1}"][1], names=names)
    assert "\\dot{\\Omega_1}" in tex and "J_1" in tex and "R_1^T R_2" in tex
    assert "J1" not in tex and "{R1}" not in tex


def test_latex_names_one_to_one():
    # 'J1' would become 'J_1', which another leaf already is: neither is renamed
    a, b = Matrix("J1"), Matrix("J_1")
    x = Vector("x")
    assert "J1" not in latex_names(a * x + b * x)
    # two SO3 points whose vectors would both print \Omega_1
    R, Q = SO3("R1"), SO3("Q1")
    names = latex_names(Dot(R * R.get_tangent_vector(), Q * Q.get_tangent_vector()))
    assert "\\Omega_{R1}" not in names and "\\Omega_{Q1}" not in names
    assert names["R1"] == "R_1" and names["Q1"] == "Q_1"


def test_to_latex_names_default_unchanged_and_leaf_only():
    J1, m1 = Matrix("J1"), Scalar("m1")
    Om = SO3("R1").get_tangent_vector()
    expr = J1 * Om.t_diff() * m1
    assert to_latex(expr) == to_latex(expr, names=None) == to_latex(expr, names={})
    # the map is applied to leaf names only: a formula containing the text 'J1' elsewhere
    # (here, inside another leaf's name) is not touched
    other = Matrix("J1x")
    assert to_latex(other * Om, names={"J1": "J_1"}) == to_latex(other * Om)
    out = to_latex(expr, names={"J1": "J_1", "\\Omega_{R1}": "\\Omega_1"})
    assert "J_1 \\dot{\\Omega_1}" in out
    # nested calls do not leak the map
    assert "J1" in to_latex(expr)


# ---------------------------------------------------------------------------
# Scene geometry against NumPy
# ---------------------------------------------------------------------------


def _rot(state, sl):
    return state[:, sl].reshape(-1, 3, 3)


def test_double_pendulum_frames_match_numpy():
    d = derived(BUILDERS["double_rigid_pendulum"])
    s = d.system.symbols
    p, _ = SETUPS["double_rigid_pendulum"]
    traj = _short_traj("double_rigid_pendulum")
    scene = Scene(d.layout, params=p)
    scene.add(RigidBody(s["R1"], ZeroVector, Box(0.1, 0.04, 1.0), offset=[0, 0, -0.5]))
    scene.add(RigidBody(s["R2"], s["R1"] * s["l1"], Box(0.1, 0.04, 0.8), offset=s["rho2"]))
    scene.add(Trace(s["R1"] * s["l1"] + s["R2"] * (s["rho2"] * 2)))
    fr = scene.frames(traj)
    X = traj.state
    R1, R2 = _rot(X, slice(0, 9)), _rot(X, slice(9, 18))
    l1, rho2 = p["l_1"], p["\\rho_2"]
    assert len(fr) == len(traj.t) and np.array_equal(fr.t, traj.t)
    np.testing.assert_allclose(fr.bodies[0]["rotation"], R1, atol=TOL, rtol=0)
    np.testing.assert_allclose(fr.bodies[1]["rotation"], R2, atol=TOL, rtol=0)
    np.testing.assert_allclose(fr.bodies[0]["position"], 0.0, atol=TOL)
    np.testing.assert_allclose(fr.bodies[0]["center"], R1 @ [0, 0, -0.5], atol=TOL, rtol=0)
    np.testing.assert_allclose(fr.bodies[1]["position"], R1 @ l1, atol=TOL, rtol=0)
    np.testing.assert_allclose(fr.bodies[1]["center"], R1 @ l1 + R2 @ rho2, atol=TOL, rtol=0)
    np.testing.assert_allclose(fr.traces[0], R1 @ l1 + 2 * (R2 @ rho2), atol=TOL, rtol=0)


def test_quadrotor_load_frames_match_numpy():
    d = derived(BUILDERS["quadrotor_with_load"])
    s = d.system.symbols
    p, _ = SETUPS["quadrotor_with_load"]
    traj = _short_traj("quadrotor_with_load")
    xL, q, R, l = s["xL"], s["q"], s["R"], s["l"]
    xQ = xL - q * l
    offsets = [(0.25, 0, 0.02), (0, 0.25, 0.02), (-0.25, 0, 0.02), (0, -0.25, 0.02)]
    scene = Scene(d.layout, params=p)
    scene.add(RigidBody(R, xQ, Box(0.16, 0.12, 0.05)))
    for off in offsets:
        scene.add(RigidBody(R, xQ, Cylinder(0.08, 0.01), offset=off, axes=False))
    scene.add(Rod(xQ, xL))
    scene.add(PointMass(xL))
    fr = scene.frames(traj)
    X = traj.state
    x_L, q_n, R_n = X[:, 0:3], X[:, 3:6], _rot(X, slice(6, 15))
    x_Q = x_L - p["l"] * q_n
    for i, body in enumerate(fr.bodies):
        off = np.zeros(3) if i == 0 else np.array(offsets[i - 1], float)
        np.testing.assert_allclose(body["rotation"], R_n, atol=TOL, rtol=0)
        np.testing.assert_allclose(body["position"], x_Q, atol=TOL, rtol=0)
        np.testing.assert_allclose(body["center"], x_Q + R_n @ off, atol=TOL, rtol=0)
    np.testing.assert_allclose(fr.rods[0][:, 0], x_Q, atol=TOL, rtol=0)
    np.testing.assert_allclose(fr.rods[0][:, 1], x_L, atol=TOL, rtol=0)
    np.testing.assert_allclose(fr.points[0], x_L, atol=TOL, rtol=0)


def test_spherical_pendulum_frames_match_numpy():
    d = derived(BUILDERS["spherical_pendulum_s2"])
    s = d.system.symbols
    p, _ = SETUPS["spherical_pendulum_s2"]
    traj = _short_traj("spherical_pendulum_s2")
    bob = s["q"] * s["l"]
    scene = Scene(d.layout, params=p)
    scene.add(Rod(ZeroVector, bob))
    scene.add(PointMass(bob))
    fr = scene.frames(traj)
    q = traj.state[:, 0:3]
    np.testing.assert_allclose(fr.points[0], p["l"] * q, atol=TOL, rtol=0)
    np.testing.assert_allclose(fr.rods[0][:, 0], 0.0, atol=TOL)
    np.testing.assert_allclose(fr.rods[0][:, 1], p["l"] * q, atol=TOL, rtol=0)


def test_example_rotations_are_orthogonal_and_spin_is_visible(dp_report):
    report, _ = dp_report
    fr = report.scene.frames(report.trajectory)
    for body in fr.bodies:
        R = body["rotation"]
        err = np.abs(np.einsum("nki,nkj->nij", R, R) - np.eye(3)).max()
        assert err < TOL, err
        np.testing.assert_allclose(np.linalg.det(R), 1.0, atol=TOL)
    # spin about each link's own axis: the body-frame Ω along e3 integrates to several turns
    Om1 = report.trajectory.state[:, 18:21]
    turns = np.trapezoid(Om1[:, 2], report.trajectory.t) / (2 * np.pi)
    assert abs(turns) > 1.0, turns


def test_spin_about_link_axis():
    """A rigid pendulum hanging straight (ρ along e3) with Ω along the link spins about
    it at a constant rate: the body e1 axis turns about the link axis, R(t) = Rz(w t)."""
    d = derived(BUILDERS["rigid_pendulum_pivot"])
    p = {"J": np.diag([0.3, 0.3, 0.02]), "\\rho": np.array([0, 0, -0.5]),
         "e3": np.array([0.0, 0, 1]), "g": 9.81, "m": 1.0}  # fmt: skip
    w = 2.0
    x0 = d.layout.pack(R=np.eye(3), Omega_R=[0, 0, w])
    traj = simulate(standard_module("rigid_pendulum_pivot"), x0, 2.0, p, output_every=50)
    scene = Scene(d.layout, params=p)
    scene.add(RigidBody(d.system.symbols["R"], ZeroVector, Box(0.1, 0.04, 1.0), [0, 0, -0.5]))
    R = scene.frames(traj).bodies[0]["rotation"]
    e1, e3 = R[:, :, 0], R[:, :, 2]
    np.testing.assert_allclose(e3, np.tile([0, 0, 1.0], (len(R), 1)), atol=1e-9)
    expected = np.stack([np.cos(w * traj.t), np.sin(w * traj.t), 0 * traj.t], axis=1)
    np.testing.assert_allclose(e1, expected, atol=1e-9)
    angle = np.unwrap(np.arctan2(e1[:, 1], e1[:, 0]))
    assert angle[-1] - angle[0] > np.pi  # more than half a turn: visible


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def dp_setup():
    d = derived(BUILDERS["double_rigid_pendulum"])
    p, _ = SETUPS["double_rigid_pendulum"]
    return d, p


def test_scene_rejects_velocity_input_and_stray_leaves(dp_setup):
    d, p = dp_setup
    s = d.system.symbols
    scene = Scene(d.layout, params=p)
    Om1 = s["R1"].get_tangent_vector()
    with pytest.raises(ValueError, match=r"velocity '\\Omega_\{R1\}'"):
        scene.add(PointMass(s["R1"] * Om1))
    with pytest.raises(ValueError, match="neither a configuration variable nor a parameter"):
        scene.add(Trace(s["R1"] * Vector("d")))
    with pytest.raises(ValueError, match="'M1'"):
        scene.add(Rod(ZeroVector, Vector("M1")))
    with pytest.raises(ValueError, match="must be a vector"):
        scene.add(PointMass(s["m1"]))
    with pytest.raises(ValueError, match="must be a matrix"):
        scene.add(RigidBody(s["l1"], ZeroVector, Sphere(0.1)))
    with pytest.raises(ValueError, match="must be constant"):
        scene.add(RigidBody(s["R1"], ZeroVector, Sphere(0.1), offset=s["R2"] * s["l1"]))
    with pytest.raises(TypeError, match="Box, Cylinder or Sphere"):
        RigidBody(s["R1"], ZeroVector, "box")
    with pytest.raises(ValueError, match="positive"):
        Box(0.1, -1, 1)
    with pytest.raises(ValueError, match="color"):
        PointMass(ZeroVector, color="blue")
    assert not scene.bodies and not scene.points and not scene.rods and not scene.traces


def test_scene_rejects_non_rotation_and_wrong_layout(dp_setup):
    d, p = dp_setup
    s = d.system.symbols
    scene = Scene(d.layout, params=p)
    scene.add(RigidBody(s["R1"] * s["J1"], ZeroVector, Sphere(0.1), name="bad"))
    _, x0 = SETUPS["double_rigid_pendulum"]
    with pytest.raises(ValueError, match="is not a rotation"):
        scene.frames(x0)
    with pytest.raises(ValueError, match="24 entries"):
        Scene(d.layout, p).frames(np.zeros(5))
    other = _short_traj("rigid_pendulum_pivot", 0.01)
    with pytest.raises(ValueError):
        Scene(d.layout, p).frames(other)


# ---------------------------------------------------------------------------
# The HTML report
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def dp_report(tmp_path_factory):
    """The double-pendulum example: 10 s at dt = 1e-3, default output_every."""
    out = tmp_path_factory.mktemp("viz") / "double_rigid_pendulum.html"
    t0 = time.perf_counter()
    report = _load_example("report_double_rigid_pendulum").main(str(out))
    report.generation_time = time.perf_counter() - t0
    return report, out


class _Parser(html.parser.HTMLParser):
    VOID = {"meta", "link", "br", "img", "input", "hr", "source", "area", "col", "wbr"}

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack, self.errors, self.urls = [], [], []
        self.math, self._math_depth, self._buf = [], None, []
        self.scripts, self._script = {}, None
        self.canvas_host = False

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        for key in ("src", "href"):
            if a.get(key):
                self.urls.append(a[key])
        if tag in self.VOID:
            return
        self.stack.append(tag)
        if "math" in (a.get("class") or "").split() and self._math_depth is None:
            self._math_depth, self._buf = len(self.stack), []
        if tag == "script":
            self._script = a.get("id") or f"script{len(self.scripts)}"
            self.scripts[self._script] = ""
        if a.get("id") == "gm-viewer":
            self.canvas_host = True

    def handle_endtag(self, tag):
        if tag in self.VOID:
            return
        if not self.stack or self.stack[-1] != tag:
            self.errors.append(f"unexpected </{tag}> (open: {self.stack[-3:]})")
            return
        if self._math_depth == len(self.stack):
            self.math.append("".join(self._buf))
            self._math_depth = None
        self.stack.pop()
        if tag == "script":
            self._script = None

    def handle_data(self, data):
        if self._math_depth is not None:
            self._buf.append(data)
        if self._script is not None:
            self.scripts[self._script] += data


@pytest.fixture(scope="module")
def parsed(dp_report):
    _, out = dp_report
    text = out.read_text(encoding="utf-8")
    p = _Parser()
    p.feed(text)
    p.close()
    return text, p


def test_html_parses(parsed):
    text, p = parsed
    assert text.startswith("<!doctype html>")
    assert not p.errors, p.errors
    assert p.stack == [], p.stack
    assert p.canvas_host


def test_only_pinned_external_urls(parsed):
    text, p = parsed
    found = set(re.findall(r"""(?:https?:)?//[A-Za-z0-9.-]+\.[a-z]{2,}[^\s"'<>)]*""", text))
    assert found == set(EXTERNAL_URLS), found - set(EXTERNAL_URLS)
    assert set(p.urls) == set(EXTERNAL_URLS)
    assert "three@0.147.0" in EXTERNAL_URLS[0] and "mathjax/3.2.2" in EXTERNAL_URLS[2]
    for pattern in ("@import", "url(", "fetch(", "XMLHttpRequest", "import("):
        assert pattern not in text


def test_inline_json_round_trips(parsed, dp_report):
    report, _ = dp_report
    _, p = parsed
    data = json.loads(p.scripts["gm-data"])
    assert data["title"] == report.title
    assert data["scene"] == json.loads(json.dumps(report.scene_data()))
    sc = data["scene"]
    n = len(sc["t"])
    assert 590 <= n <= 610  # 10 s at 60 frames per second of simulated time
    assert all(len(b["pos"]) == 3 * n for b in sc["bodies"])
    assert all(len(r) == 9 * n for r in sc["rotations"])
    # the inlined rotations are the frames' rotations to the stored precision
    idx = report.frame_indices()
    R1 = report.scene.frames(report.trajectory, idx).bodies[0]["rotation"]
    np.testing.assert_allclose(np.reshape(sc["rotations"][0], (n, 3, 3)), R1, atol=6e-6)


def test_size_bound(dp_report):
    report, out = dp_report
    assert len(report.trajectory.t) == 10001  # 10 s at the default output_every
    assert out.stat().st_size < 2_000_000, out.stat().st_size


def _balanced(tex):
    depth = 0
    for i, ch in enumerate(tex):
        if ch == "{" and (i == 0 or tex[i - 1] != "\\"):
            depth += 1
        elif ch == "}" and (i == 0 or tex[i - 1] != "\\"):
            depth -= 1
            if depth < 0:
                return False
    if depth:
        return False
    if len(re.findall(r"\\left\b", tex)) != len(re.findall(r"\\right\b", tex)):
        return False
    begins = re.findall(r"\\begin\{(\w+)\}", tex)
    ends = re.findall(r"\\end\{(\w+)\}", tex)
    return sorted(begins) == sorted(ends)


def test_mathjax_sources_balanced(parsed):
    _, p = parsed
    assert len(p.math) > 20
    for tex in p.math:
        assert tex.startswith(("\\[", "\\(")) and tex.endswith(("\\]", "\\)")), tex[:80]
        assert _balanced(tex[2:-2]), tex[:200]
        assert "J1" not in tex and "\\Omega_{R1}" not in tex  # readable names


def test_section_order(parsed):
    text, _ = parsed
    ids = ["summary", "view", "lagrangian", "eom", "standard-form", "parameters", "source"]
    pos = [text.index(f'id="{i}"') for i in ids]
    assert pos == sorted(pos)


def test_eom_one_term_per_line_and_summary(dp_report, parsed):
    report, _ = dp_report
    text, _ = parsed
    eq = report.eom["\\eta_{R1}"][1]
    n_terms = len(eq.nodes)
    block = text[text.index('id="eom"') : text.index('id="standard-form"')]
    first = block[block.index("\\begin{aligned}") : block.index("\\end{aligned}")]
    assert first.count("\\\\") == n_terms - 1
    assert report.check_error < 1e-5 and "passes" in text
    assert report.energy_drift() < 1e-8


def test_repr_html_is_iframe(dp_report):
    report, _ = dp_report
    out = report._repr_html_()
    assert out.startswith("<iframe srcdoc=") and "&lt;!doctype html&gt;" in out


def test_report_without_scene_derives_and_checks():
    d = derived(BUILDERS["spherical_pendulum_s2"])
    s = d.system
    r = Report("x", system=(s.L, s.dW, s.variables, s.inputs))
    assert r.derive_time is not None and r.check_error < 1e-5
    page = r.html()
    assert 'id="view"' not in page and 'id="standard-form"' in page
    with pytest.raises(ValueError, match="both scene and trajectory"):
        Report("x", system=(s.L, s.dW, s.variables, s.inputs), scene=Scene(d.layout))


def test_failing_check_is_shown_as_failure():
    d = derived(BUILDERS["spherical_pendulum_s2"])
    s = d.system
    key, (var, eq) = next(iter(d.eom.items()))
    wrong = {key: (var, eq + eq)}  # every term doubled: a wrong equation
    r = Report("x", system=(s.L, s.dW, s.variables, s.inputs), eom=wrong)
    assert r.check_error > 1e-3
    assert "FAILS" in r.html() and 'class="stat bad"' in r.html()


def test_template_is_package_data():
    from importlib import resources

    assert resources.files("geomech.viz").joinpath("template.html").is_file()
    assert '"geomech.viz" = ["template.html"]' in (ROOT / "pyproject.toml").read_text()


def test_imports_stay_light():
    code = (
        "import sys, geomech; assert 'geomech.viz' not in sys.modules;"
        " import geomech.viz; assert 'scipy' not in sys.modules, 'scipy imported'"
    )
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    subprocess.run([sys.executable, "-c", code], check=True, env=env)


# ---------------------------------------------------------------------------
# Headless browser
# ---------------------------------------------------------------------------


def test_headless_browser_smoke(dp_report):
    """Opt-in (``GEOMECH_BROWSER_TESTS=1``); see ``tests/browser.py``."""
    _, out = dp_report
    dom = browser.load(out.read_text(encoding="utf-8"))
    root = re.search(r"<html[^>]*>", dom).group(0)
    assert 'data-gm-errors="0"' in root, root
    assert 'data-gm-mathjax="done"' in root, root
    assert dom.count("<mjx-container") > 20
    assert "<canvas" in dom
    colors = re.search(r'data-gm-canvas-colors="(\d+)"', root)
    assert colors and int(colors.group(1)) >= 5, root  # the scene is drawn, not blank
