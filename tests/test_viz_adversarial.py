# ruff: noqa: E501  (the injected JavaScript probes are kept on readable single lines)
"""Adversarial review of ``geomech.viz``.

Confirmed findings are strict xfails (``reason="F<n>: …"``); everything else here held.
The browser tests inject a probe into the generated page, read the three.js object
matrices back from headless Chrome and compare them with poses recomputed with NumPy
from the trajectory (``tests/browser.py``: opt-in with ``GEOMECH_BROWSER_TESTS=1``;
they skip without Chrome, the CDNs or WebGL).
"""

from __future__ import annotations

import copy
import io
import json
import os
import re
import shutil
import subprocess
import sys
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
import pytest

from geomech import (
    S2,
    SO3,
    Dot,
    Matrix,
    Scalar,
    SystemVariables,
    Vector,
    check_eom,
    compute_eom,
    getScalars,
    print_eom,
    to_standard_form,
)
from geomech.codegen import StateLayout, collect_params, load_module, module_source
from geomech.core.base.expressions import ZeroVector
from geomech.numeric import expm_so3
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
from tests import browser as browser_mod
from tests.browser import when_ready
from tests.codegen_systems import derived
from tests.simulation_systems import BUILDERS, SETUPS, standard_module

ROOT = Path(__file__).resolve().parents[1]
ROUND = 2e-5  # the page stores 5 decimals

# inertial (x, y, z) is drawn at three.js (x, z, -y): three = P @ inertial
P = np.array([[1.0, 0, 0], [0, 0, 1.0], [0, -1.0, 0]])


# ---------------------------------------------------------------------------
# Systems
# ---------------------------------------------------------------------------


def _sys(d):
    return tuple(d.system)[:4]


def _slice(layout, var):
    return layout.slot(var).slice


@pytest.fixture(scope="module")
def free_body():
    """A box long along e₁ with a large spin about body e₃ (not the long axis), placed far
    from the inertial origin."""
    d = derived(BUILDERS["free_rigid_body"])
    R = d.system.symbols["R"]
    p = {"J": np.diag([0.0042, 0.084, 0.087])}
    x0 = d.layout.pack(R=expm_so3([0.3, 0.2, 0.1]), Omega_R=[0.05, 0.3, 9.0])
    traj = simulate(standard_module("free_rigid_body"), x0, 1.5, p, output_every=5)
    origin = np.array([50.0, -30.0, 20.0])
    scene = Scene(d.layout, params=p)
    scene.add(RigidBody(R, origin=origin, shape=Box(1.0, 0.2, 0.05), offset=[0.1, 0, 0]))
    rep = Report("Free body", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene,
                 trajectory=traj, params=p)  # fmt: skip
    X = traj.state[rep.frame_indices()]
    Rn = X[:, _slice(d.layout, R)].reshape(-1, 3, 3)
    exp = dict(bodies=[(Rn, origin + Rn @ np.array([0.1, 0, 0]))], points=[], rods=[], traces=[])
    return rep, exp


@pytest.fixture(scope="module")
def double_s2():
    """Two chained S² pendulums, derived here (not one of the pinned systems)."""
    half = Scalar("0.5", value=0.5, attr=["Constant"])
    m1, m2, l1, l2, g = getScalars("m1 m2 l1 l2 g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    q1, q2 = S2("q1"), S2("q2")
    x1 = q1 * l1
    x2 = x1 + q2 * l2
    v1, v2 = x1.t_diff(), x2.t_diff()
    L = (half * m1 * Dot(v1, v1) + half * m2 * Dot(v2, v2) - m1 * g * Dot(x1, e3)
         - m2 * g * Dot(x2, e3))  # fmt: skip
    f1, f2 = Vector("f1"), Vector("f2")
    dW = Dot(q1.get_variation_vector(), f1) + Dot(q2.get_variation_vector(), f2)
    V, inputs = SystemVariables(vectors=[q1, q2]), [f1, f2]
    eom = compute_eom(L, dW, V)
    assert check_eom(L, dW, V, eom, inputs=inputs) < 1e-5
    sf = to_standard_form(eom, V, inputs)
    layout = StateLayout(V, inputs, collect_params([L, dW], V, inputs))
    plant = load_module(module_source(sf, layout, lagrangian=L), "adv_double_s2")
    p = {"m1": 1.0, "m2": 0.7, "l1": 1.0, "l2": 0.8, "g": 9.81, "e3": np.array([0, 0, 1.0])}
    a = np.array([np.sin(1.0), 0, -np.cos(1.0)])
    b = np.array([0, np.sin(0.5), -np.cos(0.5)])
    w1 = np.array([0, 1.0, 0.0]) - a[1] * a
    traj = simulate(plant, layout.pack(q1=a, q2=b, omega_q1=w1), 2.0, p, output_every=5)
    scene = Scene(layout, params=p)
    scene.add(Rod(ZeroVector, x1))
    scene.add(Rod(x1, x2))
    scene.add(PointMass(x1, radius=0.06))
    scene.add(PointMass(x2, radius=0.06))
    scene.add(Trace(x2))
    rep = Report("Double spherical pendulum", system=(L, dW, V, inputs), eom=eom,
                 standard_form=sf, scene=scene, trajectory=traj, params=p)  # fmt: skip
    X = traj.state[rep.frame_indices()]
    X1 = p["l1"] * X[:, _slice(layout, q1)]
    X2 = X1 + p["l2"] * X[:, _slice(layout, q2)]
    exp = dict(bodies=[], points=[X1, X2], rods=[(0 * X1, X1), (X1, X2)], traces=[X2])
    return rep, exp


@pytest.fixture(scope="module")
def dp_short():
    d = derived(BUILDERS["double_rigid_pendulum"])
    s = d.system.symbols
    p, x0 = SETUPS["double_rigid_pendulum"]
    traj = simulate(standard_module("double_rigid_pendulum"), x0, 1.0, p, output_every=5)
    scene = Scene(d.layout, params=p)
    scene.add(RigidBody(s["R1"], ZeroVector, Box(0.12, 0.05, 1.0), offset=[0, 0, -0.5]))
    scene.add(RigidBody(s["R2"], s["R1"] * s["l1"], Box(0.12, 0.05, 0.8), offset=s["rho2"]))
    scene.add(Trace(s["R1"] * s["l1"] + s["R2"] * (s["rho2"] * 2)))
    rep = Report("dp", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene,
                 trajectory=traj, params=p)  # fmt: skip
    X = traj.state[rep.frame_indices()]
    R1 = X[:, _slice(d.layout, s["R1"])].reshape(-1, 3, 3)
    R2 = X[:, _slice(d.layout, s["R2"])].reshape(-1, 3, 3)
    l1, rho2 = p["l_1"], p["\\rho_2"]
    exp = dict(
        bodies=[(R1, R1 @ [0, 0, -0.5]), (R2, R1 @ l1 + R2 @ rho2)],
        points=[], rods=[], traces=[R1 @ l1 + 2 * (R2 @ rho2)],
    )  # fmt: skip
    return rep, exp


@pytest.fixture(scope="module")
def quad_short():
    d = derived(BUILDERS["quadrotor_with_load"])
    s = d.system.symbols
    p, x0 = SETUPS["quadrotor_with_load"]
    traj = simulate(standard_module("quadrotor_with_load"), x0, 1.0, p, output_every=5)
    xL, q, R, l = s["xL"], s["q"], s["R"], s["l"]
    xQ = xL - q * l
    offs = [(0, 0, 0), (0.25, 0, 0.02), (0, 0.25, 0.02), (-0.25, 0, 0.02), (0, -0.25, 0.02)]
    scene = Scene(d.layout, params=p)
    scene.add(RigidBody(R, xQ, Box(0.16, 0.12, 0.05)))
    for off in offs[1:]:
        scene.add(RigidBody(R, xQ, Cylinder(0.08, 0.01), offset=off, axes=False))
    scene.add(Rod(xQ, xL))
    scene.add(PointMass(xL))
    scene.add(Trace(xL))
    rep = Report("quad", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene,
                 trajectory=traj, params=p)  # fmt: skip
    X = traj.state[rep.frame_indices()]
    x_L = X[:, _slice(d.layout, xL)]
    Rn = X[:, _slice(d.layout, R)].reshape(-1, 3, 3)
    x_Q = x_L - p["l"] * X[:, _slice(d.layout, q)]
    exp = dict(
        bodies=[(Rn, x_Q + Rn @ np.array(o, float)) for o in offs],
        points=[x_L], rods=[(x_Q, x_L)], traces=[x_L],
    )  # fmt: skip
    return rep, exp


def _body_with_pendulum():
    """New system, built by following skills/geomech/SKILL.md only: an S² pendulum
    hanging from a rigid body that swings about a fixed pivot (SO(3) × S²)."""
    half = Scalar("0.5", value=0.5, attr=["Constant"])
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    rho, d = Vector("\\rho", attr=["Constant"]), Vector("d", attr=["Constant"])
    mB, m, g, l = getScalars("m_B m g l", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    R, q = SO3("R"), S2("q")
    Om = R.get_tangent_vector()
    M, f = Vector("M"), Vector("f")
    xB, xp = R * rho, R * d + q * l
    vB, vp = xB.t_diff(), xp.t_diff()
    L = (half * Dot(Om, J * Om) + half * mB * Dot(vB, vB) + half * m * Dot(vp, vp)
         - mB * g * Dot(xB, e3) - m * g * Dot(xp, e3))  # fmt: skip
    dW = Dot(R.get_variation_vector(), M) + Dot(q.get_variation_vector(), f)
    return L, dW, SystemVariables(vectors=[q], matrices=[R]), [M, f], dict(R=R, q=q, d=d, l=l)


@pytest.fixture(scope="module")
def body_pendulum():
    L, dW, V, inputs, s = _body_with_pendulum()
    eom = compute_eom(L, dW, V)
    err = check_eom(L, dW, V, eom, inputs=inputs)
    sf = to_standard_form(eom, V, inputs)
    layout = StateLayout(V, inputs, collect_params([L, dW], V, inputs))
    plant = load_module(module_source(sf, layout, lagrangian=L), "adv_body_pendulum")
    p = {"J": np.diag([0.09, 0.08, 0.004]), "\\rho": np.array([0.0, 0.0, -0.4]),
         "d": np.array([0.05, 0.0, -0.8]), "m_B": 1.0, "m": 0.5, "g": 9.81, "l": 0.6,
         "e3": np.array([0.0, 0.0, 1.0])}  # fmt: skip
    q0 = np.array([np.sin(0.6), 0.0, -np.cos(0.6)])
    w0 = np.array([0.0, 2.0, 0.5]) - (np.array([0.0, 2.0, 0.5]) @ q0) * q0
    x0 = layout.pack(R=expm_so3([0.4, -0.2, 0.0]), q=q0, Omega_R=[0.2, 0.1, 6.0], omega_q=w0)
    traj = simulate(plant, x0, 2.0, p, output_every=5)
    R, q, d, l = s["R"], s["q"], s["d"], s["l"]
    scene = Scene(layout, params=p)
    scene.add(RigidBody(R, ZeroVector, Box(0.12, 0.05, 0.9), offset=[0, 0, -0.45]))
    scene.add(Rod(R * d, R * d + q * l))
    scene.add(PointMass(R * d + q * l, radius=0.06))
    scene.add(Trace(R * d + q * l))
    rep = Report("Body with a pendulum", system=(L, dW, V, inputs), eom=eom, standard_form=sf,
                 scene=scene, trajectory=traj, params=p)  # fmt: skip
    X = traj.state[rep.frame_indices()]
    Rn = X[:, _slice(layout, R)].reshape(-1, 3, 3)
    hinge = Rn @ p["d"]
    bob = hinge + p["l"] * X[:, _slice(layout, q)]
    exp = dict(bodies=[(Rn, Rn @ [0, 0, -0.45])], points=[bob], rods=[(hinge, bob)],
               traces=[bob])  # fmt: skip
    return rep, exp, err, eom


# ---------------------------------------------------------------------------
# Headless Chrome with a probe injected into the page
# ---------------------------------------------------------------------------

_MAIN = '<script>\n(function () {\n  "use strict";'

_PRE = r"""<script>
(function () {
  var O = THREE.Scene;
  THREE.Scene = function () { var s = new O(); window.__gmScene = s; return s; };
  THREE.Scene.prototype = O.prototype;
  var C = THREE.PerspectiveCamera;
  THREE.PerspectiveCamera = function (a, b, c, d) { var s = new C(a, b, c, d); window.__gmCam = s; return s; };
  THREE.PerspectiveCamera.prototype = C.prototype;
})();
</script>
"""

_RAF = """<script>
window.__rafQ = [];
window.requestAnimationFrame = function (cb) { window.__rafQ.push(cb); return window.__rafQ.length; };
</script>
"""

_POST = r"""<script>
(function () {
  var FR = __FRAMES__, out = {};
  function m4(m) { return Array.prototype.slice.call(m.elements); }
  function hex(c) { return c.getHexString(); }
  try {
    var sc = window.__gmScene, cam = window.__gmCam;
    var world = sc.children.filter(function (c) { return c.isGroup; })[0];
    var play = document.getElementById("gm-play");
    if (play.textContent === "Pause") play.click();
    var scrub = document.getElementById("gm-scrub");
    sc.updateMatrixWorld(true); cam.updateMatrixWorld(true);
    out.world = m4(world.matrixWorld);
    out.cameraUp = cam.up.toArray();
    var tokens = ["--axis-x", "--axis-y", "--axis-z"].map(function (n) {
      return getComputedStyle(document.documentElement).getPropertyValue(n).trim().slice(1).toLowerCase(); });
    function arrows(grp) {
      return grp.children.filter(function (c) { return c.isGroup && c.children.length === 2; }).map(function (a) {
        var shaft = a.children[0], head = a.children[1];
        return { dir: new THREE.Vector3(0, 1, 0).transformDirection(a.matrixWorld).toArray(),
                 color: hex(shaft.material.color),
                 length: shaft.geometry.parameters.height + head.geometry.parameters.height,
                 radius: shaft.geometry.parameters.radiusTop };
      });
    }
    function isTriad(g) {
      var a = arrows(g);
      return a.length >= 3 && a[0].color === tokens[0] && a[1].color === tokens[1] && a[2].color === tokens[2];
    }
    var groups = world.children.filter(function (c) { return c.isGroup && c.matrixAutoUpdate === false; });
    var pts = world.children.filter(function (c) { return c.isMesh && c.geometry.type === "SphereGeometry"; });
    var rods = world.children.filter(function (c) { return c.isMesh && c.geometry.type === "CylinderGeometry"; });
    var tg = world.children.filter(function (c) { return c.isGroup && c.children.length && c.children[0].isLine && !c.children[0].isLineSegments; })[0];
    // inertial triads (not attached to a body) and whether one is inside the view frustum
    var triads = [];
    world.traverse(function (o) {
      if (o.isGroup && groups.indexOf(o) < 0 && (!o.parent || groups.indexOf(o.parent) < 0) && isTriad(o)) {
        var v = new THREE.Vector3().setFromMatrixPosition(o.matrixWorld).project(cam);
        triads.push({ ndc: v.toArray(), visible: Math.abs(v.x) <= 1 && Math.abs(v.y) <= 1 && v.z >= -1 && v.z <= 1,
                      axes: arrows(o) });
      }
    });
    out.triads = triads;
    // spin visibility: is a body (axes hidden) unchanged by a quarter turn about body e3?
    out.bodyStatic = groups.map(function (g) {
      var kids = g.children.filter(function (c) { return !c.isGroup; });
      var colors = {}, kinds = [];
      kids.forEach(function (c) {
        kinds.push(c.geometry.type);
        if (c.isMesh) (Array.isArray(c.material) ? c.material : [c.material]).forEach(function (m) { colors[hex(m.color)] = 1; });
      });
      var sym = true;
      kids.forEach(function (c) {
        var pos = c.geometry.attributes.position, set = {};
        function key(x, y, z) { return [x, y, z].map(function (v) { return Math.round(v * 1e4); }).join(","); }
        for (var i = 0; i < pos.count; i++) set[key(pos.getX(i), pos.getY(i), pos.getZ(i))] = 1;
        for (var j = 0; j < pos.count; j++) if (!set[key(-pos.getY(j), pos.getX(j), pos.getZ(j))]) { sym = false; break; }
      });
      return { kinds: kinds, nColors: Object.keys(colors).length, quarterTurnInvariant: sym,
               axes: (g.children.filter(function (c) { return c.isGroup; })[0] || null) !== null };
    });
    out.frames = FR.map(function (k) {
      scrub.value = String(k);
      scrub.dispatchEvent(new Event("input"));
      sc.updateMatrixWorld(true);
      return {
        k: k, label: document.getElementById("gm-time").textContent,
        bodies: groups.map(function (g) {
          var ax = g.children.filter(function (c) { return c.isGroup; })[0];
          return { m: m4(g.matrixWorld), axes: ax ? arrows(ax) : null };
        }),
        points: pts.map(function (p) { return new THREE.Vector3().setFromMatrixPosition(p.matrixWorld).toArray(); }),
        rods: rods.map(function (r) { return [r.localToWorld(new THREE.Vector3(0, 0, 0)).toArray(), r.localToWorld(new THREE.Vector3(0, 1, 0)).toArray()]; }),
        traces: tg ? tg.children.map(function (l) {
          var a = l.geometry.attributes.position.array, n = l.geometry.drawRange.count, i = Math.min(n, a.length / 3) - 1;
          return { count: n, last: new THREE.Vector3(a[3 * i], a[3 * i + 1], a[3 * i + 2]).applyMatrix4(l.matrixWorld).toArray() };
        }) : []
      };
    });
  } catch (e) { out.error = String(e) + "\n" + (e.stack || ""); }
  var pre = document.createElement("pre"); pre.id = "gm-probe"; pre.textContent = JSON.stringify(out);
  document.body.appendChild(pre);
})();
</script>
"""

_PUMP = """<script>(function () {
  document.getElementById("gm-speed").value = "__SPEED__";
  var s = document.getElementById("gm-scrub"), lab = document.getElementById("gm-time");
  var rec = [], hz = 60, T = __T__, seen = {};
  for (var i = 0; i <= Math.round(T * hz); i++) {
    var q = window.__rafQ; window.__rafQ = [];
    q.forEach(function (cb) { cb(1000 + (1000 * i) / hz); });
    seen[+s.value] = 1;
    rec.push([i / hz, +s.value, parseFloat(lab.textContent.split("=")[1])]);
  }
  var pre = document.createElement("pre"); pre.id = "gm-probe";
  pre.textContent = JSON.stringify({ rec: rec, max: +s.max, sawLast: !!seen[+s.max] });
  document.body.appendChild(pre);
})();</script>"""

_COLLECT = r"""<script>
(function () {
  var viol = [], t0 = Date.now();
  document.addEventListener("securitypolicyviolation", function (e) { viol.push(e.violatedDirective + " " + e.blockedURI); });
  function done() {
    var out = { errors: window.__gmErrors, viol: viol,
      resources: performance.getEntriesByType("resource").map(function (r) { return r.name; }),
      mathjax: document.documentElement.getAttribute("data-gm-mathjax"),
      colors: document.documentElement.getAttribute("data-gm-canvas-colors"),
      jsLinks: Array.prototype.filter.call(document.querySelectorAll("*"), function (el) {
        var h = el.getAttribute("href") || el.getAttributeNS("http://www.w3.org/1999/xlink", "href") || "";
        return /^\s*javascript:/i.test(h); }).length,
      injected: document.querySelectorAll("img, iframe, object, embed").length,
      title: document.title };
    var pre = document.createElement("pre"); pre.id = "gm-probe"; pre.textContent = JSON.stringify(out);
    document.body.appendChild(pre);
  }
  (function wait() {
    if (document.documentElement.getAttribute("data-gm-mathjax") === "done" || Date.now() - t0 > 7000) setTimeout(done, 300);
    else setTimeout(wait, 100);
  })();
})();
</script>
"""

_CSP = (
    '<meta http-equiv="Content-Security-Policy" content="default-src \'none\'; script-src'
    " 'unsafe-inline' " + " ".join(EXTERNAL_URLS) + "; style-src 'unsafe-inline';"
    " img-src data: blob:; font-src data:; connect-src 'none'\">"
)


@pytest.fixture(scope="module")
def browser():
    """Run a page in headless Chrome and return its probe output (``tests/browser.py``:
    opt-in with ``GEOMECH_BROWSER_TESTS=1``, skips without Chrome, CDNs or WebGL)."""
    browser_mod.require()

    def run(page: str, budget: int = 8000) -> dict:
        return browser_mod.probe(page, budget)

    return run


def _before_main(page: str, snippet: str) -> str:
    assert page.count(_MAIN) == 1
    return page.replace(_MAIN, snippet + _MAIN)


def _at_end(page: str, snippet: str) -> str:
    i = page.rindex("</body>")
    return page[:i] + snippet + page[i:]


def _probe(browser, page: str, frames) -> dict:
    page = _before_main(page, _PRE)
    post = when_ready(_POST.replace("__FRAMES__", json.dumps(list(frames))))
    return browser(_at_end(page, post))


def _mat(elems):
    M = np.asarray(elems, float).reshape(4, 4).T  # three.js stores column-major
    return M[:3, :3], M[:3, 3]


def _inertial(v):
    return P.T @ np.asarray(v, float)


# ---------------------------------------------------------------------------
# 1. Geometry fidelity: what the page renders against NumPy (held)
# ---------------------------------------------------------------------------

SCENES = ["dp_short", "quad_short", "free_body", "double_s2", "body_pendulum"]


@pytest.mark.parametrize("name", SCENES)
def test_rendered_poses_match_numpy(name, request, browser):
    rep, exp = request.getfixturevalue(name)[:2]
    t = np.asarray(rep.trajectory.t)[rep.frame_indices()]
    N = len(t)
    ks = sorted({0, 1, 7, N // 3, N // 2, N - 2, N - 1})
    out = _probe(browser, rep.html(), ks)
    assert 'data-gm-errors="0"' in out["root"]
    W, _ = _mat(out["world"])
    np.testing.assert_allclose(W, P, atol=1e-12)  # z-up world map, proper rotation
    assert out["cameraUp"] == [0, 1, 0]  # three.js y = inertial e3
    for tri in out["triads"][:1]:  # the inertial axes: e1 red along x, e2 green, e3 blue
        for j, a in enumerate(tri["axes"][:3]):
            np.testing.assert_allclose(_inertial(a["dir"]), np.eye(3)[j], atol=1e-9)
    assert len(out["frames"][0]["bodies"]) == len(exp["bodies"])
    for fr in out["frames"]:
        k = fr["k"]
        assert abs(float(fr["label"].split("=")[1].split()[0]) - t[k]) <= 0.0051
        for (R, c), b in zip(exp["bodies"], fr["bodies"]):
            M, pos = _mat(b["m"])
            np.testing.assert_allclose(P.T @ M, R[k], atol=ROUND)  # no transpose, handedness
            np.testing.assert_allclose(_inertial(pos), c[k], atol=ROUND)  # origin + R offset
            for j, a in enumerate(b["axes"] or []):  # body e_j along R e_j, in the e_j colour
                np.testing.assert_allclose(_inertial(a["dir"]), R[k][:, j], atol=ROUND)
                assert a["color"] == out["triads"][0]["axes"][j]["color"]
        for x, pp in zip(exp["points"], fr["points"]):
            np.testing.assert_allclose(_inertial(pp), x[k], atol=ROUND)
        for (a, b), rr in zip(exp["rods"], fr["rods"]):
            np.testing.assert_allclose(_inertial(rr[0]), a[k], atol=ROUND)
            np.testing.assert_allclose(_inertial(rr[1]), b[k], atol=ROUND)
        for x, tr in zip(exp["traces"], fr["traces"]):
            assert tr["count"] == k + 1  # drawn up to the frame shown
            np.testing.assert_allclose(_inertial(tr["last"]), x[k], atol=ROUND)


def test_box_faces_and_axis_colours(dp_short, browser):
    out = _probe(browser, dp_short[0].html(), [0])
    tri = out["triads"][0]["axes"]
    for b in out["frames"][0]["bodies"]:
        assert [a["color"] for a in b["axes"]] == [a["color"] for a in tri[:3]]
    assert all(s["nColors"] >= 3 for s in out["bodyStatic"])  # +e1, +e2 tints on boxes
    assert not any(s["quarterTurnInvariant"] for s in out["bodyStatic"])


def test_spin_about_non_long_axis_is_in_the_data(free_body):
    rep, exp = free_body
    R = exp["bodies"][0][0]
    e1 = R[:, :, 0]
    swept = np.ptp(np.unwrap(np.arctan2(e1[:, 1], e1[:, 0])))
    assert swept > 4 * np.pi  # the long axis sweeps several turns about the spin axis


# ---------------------------------------------------------------------------
# 2. Visibility of orientation
# ---------------------------------------------------------------------------


def _cyl_sphere_page():
    d = derived(BUILDERS["rigid_pendulum_pivot"])
    p, x0 = SETUPS["rigid_pendulum_pivot"]
    R = d.system.symbols["R"]
    traj = simulate(standard_module("rigid_pendulum_pivot"), x0, 0.2, p, output_every=10)
    scene = Scene(d.layout, params=p)
    scene.add(RigidBody(R, ZeroVector, Cylinder(0.05, 1.0), offset=[0, 0, -0.5], axes=False))
    scene.add(RigidBody(R, [1.0, 0, 0], Sphere(0.2), axes=False))
    return Report("cs", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene,
                  trajectory=traj, params=p).html()  # fmt: skip


def test_cylinder_and_sphere_spin_is_visible_without_axes(browser):
    out = _probe(browser, _cyl_sphere_page(), [0])
    for s in out["bodyStatic"]:
        assert s["nColors"] > 1 or not s["quarterTurnInvariant"], s


def test_inertial_frame_in_view_far_from_origin(free_body, browser):
    out = _probe(browser, free_body[0].html(), [0])
    assert any(t["visible"] for t in out["triads"]), out["triads"]


def test_inertial_frame_in_view_near_origin(dp_short, browser):
    out = _probe(browser, dp_short[0].html(), [0])
    assert any(t["visible"] for t in out["triads"])


def _scaled_report(L):
    d = derived(BUILDERS["rigid_pendulum_pivot"])
    p, x0 = SETUPS["rigid_pendulum_pivot"]
    R = d.system.symbols["R"]
    traj = simulate(standard_module("rigid_pendulum_pivot"), x0, 0.2, p, output_every=10)
    scene = Scene(d.layout, params=p)
    scene.add(RigidBody(R, ZeroVector, Box(0.12 * L, 0.05 * L, L), offset=[0, 0, -L / 2]))
    return Report("s", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene,
                  trajectory=traj, params=p)  # fmt: skip


def test_millimetre_scene_is_framed():
    sd = _scaled_report(1e-3).scene_data()
    assert sd["radius"] < 0.01, sd["radius"]


def test_metre_and_kilometre_scene_framing():
    for L in (1.0, 1000.0):
        sd = _scaled_report(L).scene_data()
        assert 0.25 * L < sd["radius"] < 2 * L


def test_kilometre_body_axes_are_visible(browser):
    out = _probe(browser, _scaled_report(1000.0).html(), [0])
    for a in out["frames"][0]["bodies"][0]["axes"]:
        assert a["radius"] / a["length"] > 0.005, a


# ---------------------------------------------------------------------------
# 3. Time
# ---------------------------------------------------------------------------


def _pump(browser, page, speed, T):
    page = _before_main(page, _RAF)
    pump = when_ready(_PUMP.replace("__SPEED__", speed).replace("__T__", str(T)))
    return browser(_at_end(page, pump))


@pytest.mark.parametrize("speed", ["1", "2"])
def test_playback_is_real_time_times_speed(double_s2, browser, speed):
    rep, _ = double_s2
    T_end = float(rep.trajectory.t[-1] - rep.trajectory.t[0])
    out = _pump(browser, rep.html(), speed, T_end / float(speed) - 0.1)
    for wall, _k, shown in out["rec"][2:]:
        assert abs(shown - float(speed) * wall) <= float(speed) * 2.5 / 60 + 0.006


def test_playback_shows_last_frame(double_s2, browser):
    rep, _ = double_s2
    T_end = float(rep.trajectory.t[-1] - rep.trajectory.t[0])
    assert _pump(browser, rep.html(), "1", T_end + 0.5)["sawLast"]


def test_frame_thinning_keeps_ends_and_spacing(dp_short):
    rep, _ = dp_short
    t = np.asarray(rep.trajectory.t)
    idx = rep.frame_indices()
    assert idx[0] == 0 and idx[-1] == len(t) - 1
    gaps = np.diff(t[idx])
    assert gaps.max() <= 1 / 60 + (t[1] - t[0]) + 1e-12  # nothing dropped beyond 1/fps


# ---------------------------------------------------------------------------
# 1b. Parameters: the scene must draw the simulated system
# ---------------------------------------------------------------------------


def test_scene_params_must_match_trajectory():
    d = derived(BUILDERS["spherical_pendulum_s2"])
    p, x0 = SETUPS["spherical_pendulum_s2"]
    traj = simulate(standard_module("spherical_pendulum_s2"), x0, 0.2, p, output_every=10)
    q, l = d.system.symbols["q"], d.system.symbols["l"]
    scene = Scene(d.layout, params=dict(p, l=2.5))
    scene.add(PointMass(q * l))
    try:
        sd = Report("x", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene,
                    trajectory=traj).scene_data()  # fmt: skip
    except ValueError:
        return
    bob = np.asarray(sd["points"][0]["pos"]).reshape(-1, 3)
    np.testing.assert_allclose(np.linalg.norm(bob, axis=1), traj.p["l"], atol=1e-4)


# ---------------------------------------------------------------------------
# 4. LaTeX names
# ---------------------------------------------------------------------------


def test_latex_map_collisions_and_unicode():
    assert latex_names(Matrix("J1") * Vector("x") + Matrix("J_1") * Vector("x")) == {}
    q, w = S2("q1"), Scalar("\\omega1")
    names = latex_names(Dot(q.get_tangent_vector(), q.get_tangent_vector()) * w)
    assert "\\omega1" not in names and "\\omega_{q1}" not in names  # both would be \omega_1
    for name in ("θ1", "x^d", "x^d1", "a*1", "x_1", "\\rho_1"):
        assert pretty_name(name) == name
    targets = list(latex_names(_body_with_pendulum()[:4]).values())
    assert len(targets) == len(set(targets))


def test_latex_map_keeps_argument_commands_valid():
    for name in ("\\bar1", "\\hat1", "\\vec1", "\\mathbf1", "\\tilde2"):
        assert not re.match(r"^\\[A-Za-z]+_", pretty_name(name)), pretty_name(name)


BASE = "5acff46"  # the merge base of the visualizer branch on main


def test_to_latex_default_byte_identical_to_base(tmp_path):
    """``to_latex(expr)`` without names prints every pinned system exactly as at 5acff46,
    the base of the branch on ``main`` (skipped when that commit is not in the clone)."""
    base = tmp_path / "base"
    base.mkdir()
    try:
        subprocess.run(["git", "cat-file", "-e", f"{BASE}^{{commit}}"], cwd=ROOT,
                       capture_output=True, check=True)  # fmt: skip
        arch = subprocess.run(["git", "archive", BASE, "geomech"], cwd=ROOT,
                              capture_output=True, check=True).stdout  # fmt: skip
    except (OSError, subprocess.CalledProcessError):
        pytest.skip(f"base commit {BASE} is not in this clone (a shallow checkout?)")
    subprocess.run(["tar", "-x", "-C", str(base)], input=arch, check=True)
    code = (
        "import inspect, json, sys\n"
        f"sys.path.insert(1, {str(ROOT)!r})\n"
        "from geomech import compute_eom, to_standard_form, to_latex\n"
        "import tests.vbl_systems as vs, tests.phase1_systems as p1\n"
        "out = {}\n"
        "for mod in (vs, p1):\n"
        "    for n, fn in inspect.getmembers(mod, inspect.isfunction):\n"
        "        if fn.__module__ != mod.__name__ or n.startswith('_'): continue\n"
        "        try: s = fn()\n"
        "        except TypeError: continue\n"
        "        if not (isinstance(s, tuple) and len(s) >= 4): continue\n"
        "        L, dW, V, u = s[:4]\n"
        "        try: eom = compute_eom(L, dW, V); sf = to_standard_form(eom, V, u)\n"
        "        except Exception as e: out[n] = type(e).__name__; continue\n"
        "        xs = [L] + [e for _, e in eom.values()]\n"
        "        for q in sf.values(): xs += [*q.M.values(), q.f, *q.G.values()]\n"
        "        out[n] = [to_latex(x) for x in xs]\n"
        "print(json.dumps(out, sort_keys=True))\n"
    )

    def dump(first):
        env = dict(os.environ, PYTHONPATH=str(first))
        return subprocess.run([sys.executable, "-c", code], cwd=first, env=env, check=True,
                              capture_output=True, text=True).stdout  # fmt: skip

    old, new = dump(base), dump(ROOT)
    assert sum(isinstance(v, list) for v in json.loads(new).values()) >= 40
    assert old == new


# ---------------------------------------------------------------------------
# 5. HTML and security
# ---------------------------------------------------------------------------

EVIL_TEX = r"\href{javascript:alert(1)}{click}"


def _hostile_report():
    half = Scalar("0.5", value=0.5, attr=["Constant"])
    m = Scalar(EVIL_TEX, attr=["Constant"])
    g, l = getScalars("g l", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    q, f = S2("q"), Vector("f")
    w = q.get_tangent_vector()
    L = half * m * l * l * Dot(w, w) - m * g * l * Dot(q, e3)
    dW = Dot(q.get_variation_vector(), f)
    V = SystemVariables(vectors=[q])
    layout = StateLayout(V, [f], collect_params([L, dW], V, [f]))
    p = {EVIL_TEX: 1.0, "g": 9.81, "l": 1.0, "e3": np.array([0, 0, 1.0])}
    plant = load_module(module_source(to_standard_form(compute_eom(L, dW, V), V, [f]), layout,
                                      lagrangian=L), "adv_hostile")  # fmt: skip
    traj = simulate(plant, layout.pack(q=[0, 0, -1.0]), 0.05, p, output_every=10)
    scene = Scene(layout, params=p)
    scene.add(PointMass(q * l, name='</script><script>document.title="PWNED-json"</script>'))
    return Report(
        '</title><script>document.title="PWNED-title"</script>\\(\\href{javascript:alert(2)}{t}\\)',
        system=(L, dW, V, [f]), scene=scene, trajectory=traj, params=p,
        description='<img src=x onerror="document.title=1"> \\(\\href{javascript:alert(3)}{d}\\)',
        units={"g": "\\(\\href{javascript:alert(4)}{u}\\)"},
        source='</code></pre><script>document.title="PWNED-src"</script>',
    )  # fmt: skip


def test_hostile_strings_are_escaped():
    page = _hostile_report().html()
    assert page.count("<script") == 6  # the template's own scripts only
    assert "<img" not in page
    data = re.search(r'<script type="application/json" id="gm-data">(.*?)</script>', page, re.S)
    assert json.loads(data.group(1))["scene"]["points"][0]["name"].startswith("</script>")
    assert "\\u003c/script\\u003e" in data.group(1)


def test_hostile_strings_do_not_execute(browser):
    rep = _hostile_report()
    out = browser(_at_end(rep.html(), _COLLECT), budget=10000)
    assert out["errors"] == [] and out["injected"] == 0
    assert out["title"] == rep.title  # no injected script or handler rewrote it


def test_no_javascript_links_after_typesetting(browser):
    out = browser(_at_end(_hostile_report().html(), _COLLECT), budget=10000)
    assert out["jsLinks"] == 0, out["jsLinks"]


def test_page_under_artifact_csp(quad_short, browser):
    page = quad_short[0].html().replace('<meta charset="utf-8">', '<meta charset="utf-8">' + _CSP)
    out = browser(_at_end(page, _COLLECT), budget=10000)
    assert out["viol"] == [] and out["errors"] == []
    assert set(out["resources"]) <= set(EXTERNAL_URLS)
    assert out["mathjax"] == "done" and int(out["colors"]) >= 5


# ---------------------------------------------------------------------------
# 6. Robustness
# ---------------------------------------------------------------------------


def _rp():
    d = derived(BUILDERS["rigid_pendulum_pivot"])
    p, x0 = SETUPS["rigid_pendulum_pivot"]
    traj = simulate(standard_module("rigid_pendulum_pivot"), x0, 0.2, p, output_every=10)
    return d, p, traj, d.system.symbols["R"]


def _robust_pages():
    d, p, traj, R = _rp()

    def rep(scene, tr=traj):
        return Report("r", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene,
                      trajectory=tr, params=p).html()  # fmt: skip

    one = copy.copy(traj)
    one.t, one.state = traj.t[:1], traj.state[:1]
    s1 = Scene(d.layout, params=p)
    s1.add(RigidBody(R, ZeroVector, Box(0.1, 0.04, 1.0)))
    s2 = Scene(d.layout, params=p)
    s2.add(Trace(R * d.system.symbols["rho"]))
    return {"empty": rep(Scene(d.layout, params=p)), "only_traces": rep(s2),
            "one_state": rep(s1, one)}  # fmt: skip


def test_degenerate_scenes_load_without_errors(browser):
    for name, page in _robust_pages().items():
        out = browser(_at_end(page, _COLLECT), budget=10000)
        assert out["errors"] == [], name
        assert out["mathjax"] == "done", name


def test_non_finite_state_is_rejected_by_frames():
    d, p, traj, R = _rp()
    scene = Scene(d.layout, params=p)
    scene.add(RigidBody(R, ZeroVector, Box(0.1, 0.04, 1.0)))
    st = traj.state.copy()
    st[3, 0] = np.nan
    with pytest.raises(ValueError, match="(?i)finite|nan"):
        scene.frames(st)


# ---------------------------------------------------------------------------
# 7. Skill
# ---------------------------------------------------------------------------

SKILL = ROOT / "skills" / "geomech"


def test_skill_front_matter_agent_skills_format():
    text = (SKILL / "SKILL.md").read_text(encoding="utf-8")
    fm = re.match(r"^---\n(.*?)\n---\n", text, re.S).group(1)
    fields = dict(line.split(":", 1) for line in fm.splitlines())
    name, desc = fields["name"].strip(), fields["description"].strip()
    assert re.fullmatch(r"[a-z0-9]+(-[a-z0-9]+)*", name) and len(name) <= 64
    assert name == SKILL.name
    assert 0 < len(desc) <= 1024 and "<" not in desc and ">" not in desc
    for f in ("SKILL.md", "reference.md"):
        assert (SKILL / f).read_bytes() == (
            ROOT / ".claude" / "skills" / "geomech" / f
        ).read_bytes()


def test_skill_procedure_on_new_system(body_pendulum):
    """SO(3) body + S² pendulum built from the skill alone: check_eom passes, energy is
    conserved, the report's own check agrees and the scene evaluates."""
    rep, exp, err, eom = body_pendulum
    assert err < 1e-5
    assert rep.check_error == pytest.approx(err)
    assert rep.energy_drift() < 1e-8
    assert set(eom) == {"\\xi_{q}", "\\eta_{R}"}
    fr = rep.scene.frames(rep.trajectory, indices=rep.frame_indices())
    np.testing.assert_allclose(fr.points[0], exp["points"][0], atol=1e-12)


def test_gitignore_keeps_settings_ignored():
    def ignored(path):
        r = subprocess.run(["git", "check-ignore", "-q", path], cwd=ROOT)
        return r.returncode == 0

    if shutil.which("git") is None or not (ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    assert ignored(".claude/settings.local.json")
    assert ignored(".claude/agents/x.md") and ignored(".claude/commands/x.md")
    assert not ignored(".claude/skills/geomech/SKILL.md")


def test_gitignore_unignores_only_the_geomech_skill():
    if shutil.which("git") is None or not (ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    r = subprocess.run(["git", "check-ignore", "-q", ".claude/skills/private/SKILL.md"], cwd=ROOT)
    assert r.returncode == 0


def test_skill_print_convention_matches_print_eom(body_pendulum):
    """F7 (fixed in the skill): print_eom prints LaTeX source, str/print of an equation
    gives the '*' form; the skill must say both, and both must be true."""
    text = (SKILL / "SKILL.md").read_text(encoding="utf-8")
    assert "`print_eom(eom)` prints each equation as LaTeX source" in text
    assert "`*` form" in text
    buf = io.StringIO()
    with redirect_stdout(buf):
        print_eom(body_pendulum[3])
    assert "\\int" in buf.getvalue()  # LaTeX source
    eq = next(iter(body_pendulum[3].values()))[1]
    assert "*" in str(eq)  # the '*' form
