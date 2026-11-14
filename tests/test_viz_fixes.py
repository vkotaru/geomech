# ruff: noqa: E501  (the injected JavaScript probes are kept on readable single lines)
"""Regression tests for the fixes of ``specs/reviews/viz-skill-review.md`` (F1–F9).

The browser tests reuse the headless-Chrome harness of ``tests/browser.py`` through the
``browser`` fixture of ``tests/test_viz_adversarial.py`` (opt-in with
``GEOMECH_BROWSER_TESTS=1``; they skip without Chrome, the CDNs or WebGL).
"""

from __future__ import annotations

import copy
import json
import re

import numpy as np
import pytest

from geomech import (
    S2,
    Dot,
    Scalar,
    SystemVariables,
    Vector,
    compute_eom,
    getScalars,
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
    latex_names,
    pretty_name,
)
from geomech.viz.latex import SYMBOL_COMMANDS
from geomech.viz.report import MATHJAX_URL
from geomech.viz.scene import PALETTE
from tests.browser import when_ready
from tests.codegen_systems import derived
from tests.simulation_systems import BUILDERS, SETUPS, standard_module
from tests.test_viz_adversarial import (  # noqa: F401  (browser is a fixture)
    _COLLECT,
    _PRE,
    _at_end,
    _before_main,
    _hostile_report,
    _probe,
    _rp,
    browser,
)


@pytest.fixture(scope="module")
def chrome(browser):  # noqa: F811  (the imported fixture, requested by name)
    """The headless-Chrome page runner of ``tests/test_viz_adversarial.py``."""
    return browser


P = np.array([[1.0, 0, 0], [0, 0, 1.0], [0, -1.0, 0]])


def _sys(d):
    return tuple(d.system)[:4]


def _pendulum(t_final=0.2, every=10):
    d = derived(BUILDERS["spherical_pendulum_s2"])
    p, x0 = SETUPS["spherical_pendulum_s2"]
    traj = simulate(standard_module("spherical_pendulum_s2"), x0, t_final, p, output_every=every)
    return d, p, traj


# ---------------------------------------------------------------------------
# F1: one source of truth for the parameters
# ---------------------------------------------------------------------------


def test_scene_without_params_uses_trajectory_p():
    d, p, traj = _pendulum()
    q, l = d.system.symbols["q"], d.system.symbols["l"]
    scene = Scene(d.layout)  # no params
    scene.add(PointMass(q * l))
    rep = Report("x", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene, trajectory=traj)
    bob = np.asarray(rep.scene_data()["points"][0]["pos"]).reshape(-1, 3)
    np.testing.assert_allclose(np.linalg.norm(bob, axis=1), p["l"], atol=1e-5)
    assert rep.params == traj.p


@pytest.mark.parametrize("where", ["scene", "report"])
def test_conflicting_params_name_the_parameter_and_both_values(where):
    d, p, traj = _pendulum()
    q, l = d.system.symbols["q"], d.system.symbols["l"]
    other = dict(p, l=2.5)
    scene = Scene(d.layout, params=other if where == "scene" else p)
    scene.add(PointMass(q * l))
    kw = {"params": other} if where == "report" else {}
    with pytest.raises(ValueError) as exc:
        Report("x", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene,
               trajectory=traj, **kw)  # fmt: skip
    msg = str(exc.value)
    assert "'l'" in msg and "2.5" in msg and f"{p['l']:g}" in msg and "trajectory.p" in msg
    if where == "scene":
        with pytest.raises(ValueError, match="Scene.params=.*'l' is 2.5"):
            scene.frames(traj)


def test_same_params_given_twice_are_accepted():
    d, p, traj = _pendulum()
    scene = Scene(d.layout, params=copy.deepcopy(p))
    scene.add(PointMass(d.system.symbols["q"] * d.system.symbols["l"]))
    Report("x", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene, trajectory=traj,
           params=dict(p)).scene_data()  # fmt: skip


def test_bare_states_need_params():
    d, p, traj = _pendulum()
    scene = Scene(d.layout)
    scene.add(PointMass(d.system.symbols["q"] * d.system.symbols["l"]))
    with pytest.raises(ValueError, match="pass a Trajectory"):
        scene.frames(traj.state)
    with_p = Scene(d.layout, params=p)
    with_p.add(PointMass(d.system.symbols["q"] * d.system.symbols["l"]))
    np.testing.assert_allclose(with_p.frames(traj.state).points[0], scene.frames(traj).points[0])


def test_parameter_offset_is_evaluated_with_trajectory_p():
    d = derived(BUILDERS["double_rigid_pendulum"])
    s = d.system.symbols
    p, x0 = SETUPS["double_rigid_pendulum"]
    traj = simulate(standard_module("double_rigid_pendulum"), x0, 0.05, p, output_every=10)
    scene = Scene(d.layout)  # offset rho2 is a parameter: evaluated lazily
    scene.add(RigidBody(s["R2"], s["R1"] * s["l1"], Box(0.1, 0.04, 0.8), offset=s["rho2"]))
    fr = scene.frames(traj)
    R1 = traj.state[:, 0:9].reshape(-1, 3, 3)
    R2 = traj.state[:, 9:18].reshape(-1, 3, 3)
    np.testing.assert_allclose(fr.bodies[0]["center"], R1 @ p["l_1"] + R2 @ p["\\rho_2"],
                               atol=1e-12)  # fmt: skip


# ---------------------------------------------------------------------------
# F2: sizes and rounding relative to the scene
# ---------------------------------------------------------------------------


def _scaled(L, extra=()):
    d = derived(BUILDERS["rigid_pendulum_pivot"])
    p, x0 = SETUPS["rigid_pendulum_pivot"]
    R = d.system.symbols["R"]
    traj = simulate(standard_module("rigid_pendulum_pivot"), x0, 0.2, p, output_every=10)
    scene = Scene(d.layout)
    scene.add(RigidBody(R, ZeroVector, Box(0.12 * L, 0.05 * L, L), offset=[0, 0, -L / 2]))
    k = Scalar(repr(2.0 * L), value=2.0 * L, attr=["Constant"])  # a number: rho is 0.5 m
    tip = R * d.system.symbols["rho"] * k
    for kind in extra:
        scene.add(PointMass(tip) if kind == "point" else Rod(ZeroVector, tip))
    return Report("s", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene,
                  trajectory=traj)  # fmt: skip


def test_scene_data_is_scale_invariant():
    """A 1 mm and a 1 km copy of one scene give the same picture: every length is the
    1 m one times the scale, to the stored precision, and the frame is the same."""
    base = _scaled(1.0).scene_data()
    for L in (1e-3, 1e3):
        sd = _scaled(L).scene_data()
        assert sd["radius"] == pytest.approx(L * base["radius"], rel=1e-5)
        np.testing.assert_allclose(sd["center"], L * np.asarray(base["center"]), rtol=1e-5,
                                   atol=1e-6 * L)  # fmt: skip
        np.testing.assert_allclose(sd["bodies"][0]["axisLengths"],
                                   L * np.asarray(base["bodies"][0]["axisLengths"]), rtol=1e-5)  # fmt: skip
        # positions keep six significant figures of the scene size, not 1e-5 absolute
        pos = np.asarray(sd["bodies"][0]["pos"])
        np.testing.assert_allclose(pos, L * np.asarray(base["bodies"][0]["pos"]),
                                   atol=2e-5 * L * base["radius"])  # fmt: skip


def test_millimetre_positions_are_not_rounded_away():
    rep = _scaled(1e-3)
    sd = rep.scene_data()
    fr = rep.scene.frames(rep.trajectory, indices=rep.frame_indices())
    exact = fr.bodies[0]["center"].reshape(-1)
    err = np.max(np.abs(np.asarray(sd["bodies"][0]["pos"]) - exact))
    assert err <= 1e-5 * sd["radius"], err  # was 1e-5 m: 1 % of the part


def test_default_point_and_rod_radii_follow_the_scene():
    r = {}
    for L in (1e-3, 1.0, 1e3):
        sd = _scaled(L, extra=("point", "rod")).scene_data()
        r[L] = (sd["points"][0]["radius"], sd["rods"][0]["radius"])
        assert 0 < r[L][1] < r[L][0] < 0.5 * L
    for L in (1e-3, 1e3):
        np.testing.assert_allclose(r[L], L * np.asarray(r[1.0]), rtol=1e-5)


def test_framing_includes_the_body_axis_arrows():
    sd = _scaled(1.0).scene_data()
    body = sd["bodies"][0]
    c, R = np.asarray(sd["center"]), sd["radius"]
    rot = np.asarray(sd["rotations"][0]).reshape(-1, 3, 3)
    pos = np.asarray(body["pos"]).reshape(-1, 3)
    for k in range(len(pos)):
        for j, length in enumerate(body["axisLengths"]):
            tip = pos[k] + length * rot[k][:, j]
            assert np.all(np.abs(tip - c) <= R * (1 + 1e-6)), (k, j)


# ---------------------------------------------------------------------------
# F3: orientation markers on bodies of revolution
# ---------------------------------------------------------------------------

_MARKERS = r"""<script>
(function () {
  var out = {};
  try {
    var sc = window.__gmScene, world = sc.children.filter(function (c) { return c.isGroup; })[0];
    var play = document.getElementById("gm-play"); if (play.textContent === "Pause") play.click();
    var axes = document.getElementById("gm-axes"); axes.checked = false; axes.dispatchEvent(new Event("change"));
    var scrub = document.getElementById("gm-scrub");
    var groups = world.children.filter(function (c) { return c.isGroup && c.matrixAutoUpdate === false; });
    out.frames = __FRAMES__.map(function (k) {
      scrub.value = String(k); scrub.dispatchEvent(new Event("input")); sc.updateMatrixWorld(true);
      return groups.map(function (g) {
        var c = new THREE.Vector3().setFromMatrixPosition(g.matrixWorld);
        var ax = g.children.filter(function (o) { return o.isGroup; })[0];
        return { axesVisible: ax ? ax.visible : null,
                 markers: g.children.filter(function (o) { return o.userData.gmMarker !== undefined; }).map(function (m) {
                   m.geometry.computeBoundingBox();
                   var w = m.geometry.boundingBox.getCenter(new THREE.Vector3()).applyMatrix4(m.matrixWorld).sub(c).normalize();
                   return { dir: m.userData.gmMarker, world: w.toArray(), visible: m.visible,
                            color: m.material.color.getHexString(), base: g.children[0].material.color.getHexString() };
                 }) };
      });
    });
  } catch (e) { out.error = String(e) + "\n" + (e.stack || ""); }
  var pre = document.createElement("pre"); pre.id = "gm-probe"; pre.textContent = JSON.stringify(out);
  document.body.appendChild(pre);
})();
</script>
"""


def _marker_report():
    d = derived(BUILDERS["free_rigid_body"])
    R = d.system.symbols["R"]
    p = {"J": np.diag([0.0042, 0.084, 0.087])}
    x0 = d.layout.pack(R=expm_so3([0.3, 0.2, 0.1]), Omega_R=[0.4, 0.3, 9.0])
    traj = simulate(standard_module("free_rigid_body"), x0, 0.5, p, output_every=5)
    scene = Scene(d.layout)
    scene.add(RigidBody(R, [0, 0, 0.5], Cylinder(0.3, 0.15), axes=False, color="rotor"))
    scene.add(RigidBody(R, [0.9, 0, 0.5], Cylinder(0.1, 0.6, axis="x"), axes=True))
    scene.add(RigidBody(R, [0, 0.9, 0.5], Cylinder(0.1, 0.6, axis="y"), axes=False))
    scene.add(RigidBody(R, [-0.9, 0, 0.5], Sphere(0.25), axes=False, color="load"))
    rep = Report("m", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene, trajectory=traj)
    Rn = traj.state[rep.frame_indices()][:, d.layout.slot(R).slice].reshape(-1, 3, 3)
    return rep, Rn


def test_marker_direction_follows_body_axis_with_axes_hidden(chrome):
    rep, Rn = _marker_report()
    ks = [0, len(Rn) // 2]
    page = _before_main(rep.html(), _PRE)
    out = chrome(_at_end(page, when_ready(_MARKERS.replace("__FRAMES__", json.dumps(ks)))))
    expected_dirs = [[0], [1], [0], [0, 1]]  # cylinder z, cylinder x (on e2), cylinder y, sphere
    turned = 0.0
    for k, frame in zip(ks, out["frames"]):
        for body, dirs in zip(frame, expected_dirs):
            assert body["axesVisible"] in (False, None)  # the axes toggle is off
            assert [m["dir"] for m in body["markers"]] == dirs
            for m in body["markers"]:
                assert m["visible"] and m["color"] != m["base"]
                world = P.T @ np.asarray(m["world"])  # inertial
                np.testing.assert_allclose(world, Rn[k][:, m["dir"]], atol=2e-3)
    turned = np.degrees(np.arccos(np.clip(Rn[ks[0]][:, 0] @ Rn[ks[1]][:, 0], -1, 1)))
    assert turned > 30  # the marker moved visibly between the two frames


# ---------------------------------------------------------------------------
# F4: the corner gizmo
# ---------------------------------------------------------------------------

_ALL_SCENES = r"""<script>
(function () {
  var O = THREE.Scene; window.__gmScenes = [];
  THREE.Scene = function () { var s = new O(); window.__gmScenes.push(s); return s; };
  THREE.Scene.prototype = O.prototype;
})();
</script>
"""

_GIZMO = r"""<script>
(function () {
  var out = {};
  try {
    var sc = window.__gmScenes;
    out.n = sc.length;
    var g = sc[0].children.filter(function (c) { return c.isGroup; })[0];
    var tri = g.children[0];
    out.arrows = tri.children.filter(function (c) { return c.isGroup; }).map(function (a) {
      a.updateMatrixWorld(true);
      return { dir: new THREE.Vector3(0, 1, 0).transformDirection(a.matrixWorld).toArray(),
               color: a.children[0].material.color.getHexString() };
    });
    out.tokens = ["--axis-x", "--axis-y", "--axis-z"].map(function (n) {
      return getComputedStyle(document.documentElement).getPropertyValue(n).trim().slice(1).toLowerCase(); });
  } catch (e) { out.error = String(e); }
  var pre = document.createElement("pre"); pre.id = "gm-probe"; pre.textContent = JSON.stringify(out);
  document.body.appendChild(pre);
})();
</script>
"""


def test_corner_gizmo_is_an_inertial_triad(chrome):
    d, p, traj = _pendulum()
    scene = Scene(d.layout)
    scene.add(PointMass(d.system.symbols["q"] * d.system.symbols["l"]))
    page = Report("g", system=_sys(d), scene=scene, trajectory=traj).html()
    out = chrome(_at_end(_before_main(page, _ALL_SCENES), when_ready(_GIZMO)))
    assert out["n"] == 2  # the gizmo scene and the main scene
    for j, a in enumerate(out["arrows"]):
        np.testing.assert_allclose(P.T @ np.asarray(a["dir"]), np.eye(3)[j], atol=1e-9)
        assert a["color"] == out["tokens"][j]
    assert len(out["arrows"]) == 3


# ---------------------------------------------------------------------------
# F5: TeX injection
# ---------------------------------------------------------------------------

INJECTIONS = [
    r"\href{javascript:alert(1)}{click}",
    r"\require{html}\href{javascript:alert(5)}{r}",
    r"\text{\href{javascript:alert(6)}{t}}",
    r"\class{gm-evil}{x}",
    r"\cssId{gm-evil-id}{x}",
    r"\style{position:fixed;top:0;left:0;width:100%;height:100%;background:red}{x}",
    r"\bbox[position:fixed;top:0;background:red]{x}",
    r"\unicode[gm-evil-font]{x41}",
    r"\toggle{a}{b}\endtoggle",
    r"\texttip{x}{gm-evil-tip}",
    r"\setOptions{}{}",
]


def _injected_page():
    half = Scalar("0.5", value=0.5, attr=["Constant"])
    g, l = getScalars("g l", attr=["Constant"])
    ms = [Scalar(s, attr=["Constant"]) for s in INJECTIONS]
    e3 = Vector("e3", attr=["Constant"])
    q, f = S2("q"), Vector("f")
    w = q.get_tangent_vector()
    mass = ms[0]
    for m in ms[1:]:
        mass = mass + m
    L = half * mass * l * l * Dot(w, w) - g * l * Dot(q, e3)
    dW = Dot(q.get_variation_vector(), f)
    V = SystemVariables(vectors=[q])
    layout = StateLayout(V, [f], collect_params([L, dW], V, [f]))
    p = {s: 1.0 for s in INJECTIONS} | {"g": 9.81, "l": 1.0, "e3": np.array([0, 0, 1.0])}
    plant = load_module(module_source(to_standard_form(compute_eom(L, dW, V), V, [f]), layout,
                                      lagrangian=L), "fix_injected")  # fmt: skip
    traj = simulate(plant, layout.pack(q=[0, 0, -1.0]), 0.05, p, output_every=10)
    scene = Scene(layout)
    scene.add(PointMass(q * l))
    text = " ".join(f"\\({s}\\)" for s in INJECTIONS)
    return Report(text, system=(L, dW, V, [f]), scene=scene, trajectory=traj, description=text,
                  units={s: f"\\({s}\\)" for s in INJECTIONS}).html()  # fmt: skip


_EVIL = r"""<script>
(function () {
  function done() {
    var all = Array.prototype.slice.call(document.querySelectorAll("*"));
    var out = {
      mathjax: document.documentElement.getAttribute("data-gm-mathjax"),
      jsLinks: all.filter(function (el) {
        var h = el.getAttribute("href") || el.getAttributeNS("http://www.w3.org/1999/xlink", "href") || "";
        return /^\s*javascript:/i.test(h); }).length,
      links: document.querySelectorAll("a, [href]:not(use):not(link)").length,
      evil: all.filter(function (el) {
        var s = (el.getAttribute("class") || "") + " " + (el.id || "") + " " + (el.getAttribute("style") || "") + " " + (el.getAttribute("title") || "");
        return /gm-evil|position:\s*fixed/.test(s); }).length,
      typesetInHeader: document.querySelectorAll("header mjx-container, td:last-child mjx-container").length,
      typeset: document.querySelectorAll(".math mjx-container").length
    };
    var pre = document.createElement("pre"); pre.id = "gm-probe"; pre.textContent = JSON.stringify(out);
    document.body.appendChild(pre);
  }
  var t0 = Date.now();
  (function wait() {
    if (document.documentElement.getAttribute("data-gm-mathjax") === "done" || Date.now() - t0 > 8000) setTimeout(done, 300);
    else setTimeout(wait, 100);
  })();
})();
</script>
"""


def test_tex_injection_strings_stay_inert(chrome):
    out = chrome(_at_end(_injected_page(), _EVIL), budget=12000)
    assert out["mathjax"] == "done"
    assert out["jsLinks"] == 0 and out["links"] == 0, out
    assert out["evil"] == 0, out
    assert out["typesetInHeader"] == 0  # title, description and units are not typeset
    assert out["typeset"] > 10  # the equations still are


def test_review_hostile_report_has_no_links(chrome):
    out = chrome(_at_end(_hostile_report().html(), _COLLECT), budget=10000)
    assert out["jsLinks"] == 0 and out["errors"] == [] and out["mathjax"] == "done"


# ---------------------------------------------------------------------------
# F6: non-finite states
# ---------------------------------------------------------------------------


def test_non_finite_state_names_time_and_slot():
    d, p, traj, R = _rp()
    scene = Scene(d.layout)
    scene.add(RigidBody(R, ZeroVector, Box(0.1, 0.04, 1.0)))
    bad = copy.copy(traj)
    bad.state = traj.state.copy()
    bad.state[5:, 4] = np.inf
    bad.state[7, 0] = np.nan
    t5 = f"{traj.t[5]:.6g}"
    with pytest.raises(ValueError, match=rf"t = {t5} \(stored state 5\).*R\[1, 1\] = inf"):
        scene.frames(bad)
    with pytest.raises(ValueError, match=rf"Report: the state at t = {t5}.*R\[1, 1\]"):
        Report("x", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene, trajectory=bad)
    with pytest.raises(ValueError, match=r"stored state 7\).*R\[0, 0\] = nan"):
        scene.frames(bad, indices=[0, 7])


# ---------------------------------------------------------------------------
# F7: forced runs and the energy figure
# ---------------------------------------------------------------------------


def test_forced_run_labels_energy_change():
    d, p, traj = _pendulum()
    scene = Scene(d.layout)
    scene.add(PointMass(d.system.symbols["q"] * d.system.symbols["l"]))
    assert traj.u.shape == (len(traj.t), 3) and not traj.u.any()
    free = Report("x", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene, trajectory=traj)
    assert not free.forced()
    assert "Energy drift" in free.html() and "Energy change" not in free.html()
    _, x0 = SETUPS["spherical_pendulum_s2"]
    pushed = simulate(standard_module("spherical_pendulum_s2"), x0, 0.2, p,
                      u=lambda t, x: np.array([0.0, 2.0 * t, 0.0]), output_every=10)  # fmt: skip
    np.testing.assert_allclose(pushed.u[:, 1], 2.0 * pushed.t)
    forced = Report("x", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene,
                    trajectory=pushed)  # fmt: skip
    assert forced.forced()
    page = forced.html()
    assert "Energy change" in page and "Energy drift" not in page
    assert "the work of the inputs, not an integrator error" in page


# ---------------------------------------------------------------------------
# F8: every mapped name is valid TeX
# ---------------------------------------------------------------------------

NAMES = [
    "J1", "e3", "m12", "R1", "x2", "\\theta1", "\\Omega12", "\\varphi3", "\\ell2", "\\bar1",
    "\\hat1", "\\vec1", "\\mathbf1", "\\tilde2", "\\dot3", "\\mathrm4", "\\operatorname5",
    "\\sqrt2", "\\frac12", "\\left1", "\\begin1", "\\foo1", "\\Omega_{R1}", "\\rho_1",
]  # fmt: skip
_KNOWN = SYMBOL_COMMANDS | {"Omega", "omega", "eta", "xi"}


def _valid(tex: str) -> bool:
    if tex.count("{") != tex.count("}"):
        return False
    # a subscript may only follow a letter, a closing brace or an argument-free command
    for cmd in re.findall(r"\\([A-Za-z]+)(?=[_^])", tex):
        if cmd not in _KNOWN:
            return False
    return True


def test_every_mapped_name_is_valid_tex():
    for name in NAMES:
        out = pretty_name(name)
        assert out == name or _valid(out), (name, out)
    assert pretty_name("\\theta1") == "\\theta_1" and pretty_name("\\bar1") == "\\bar1"
    S = S2("q1")
    names = latex_names(Dot(S.get_tangent_vector(), S.get_tangent_vector()))
    assert all(_valid(v) for v in names.values())


_TYPESET = r"""<!doctype html><html><head><meta charset="utf-8"><script>
window.MathJax = { startup: { pageReady: function () {
  return MathJax.startup.defaultPageReady().then(function () {
    var out = Array.prototype.map.call(document.querySelectorAll("p"), function (p) {
      return { tex: p.getAttribute("data-tex"), typeset: p.querySelectorAll("mjx-container").length,
               errors: p.querySelectorAll("[data-mjx-error], mjx-merror").length }; });
    var pre = document.createElement("pre"); pre.id = "gm-probe"; pre.textContent = JSON.stringify({ rows: out });
    document.body.appendChild(pre);
  }); } } };
</script><script src="__URL__"></script></head><body>__BODY__</body></html>"""


def test_every_mapped_name_typesets_in_mathjax(chrome):
    import html as html_mod

    mapped = [pretty_name(n) for n in NAMES if pretty_name(n) != n]
    body = "".join(
        f'<p data-tex="{html_mod.escape(t)}">\\({html_mod.escape(t, quote=False)} x\\)</p>'
        for t in mapped
    )
    page = _TYPESET.replace("__URL__", MATHJAX_URL).replace("__BODY__", body)
    rows = chrome(page, budget=10000)["rows"]
    assert len(rows) == len(mapped) >= 8
    for row in rows:
        assert row["typeset"] == 1 and row["errors"] == 0, row


# ---------------------------------------------------------------------------
# F9: playback and tints
# ---------------------------------------------------------------------------

_PLAY = r"""<script>
window.__rafQ = []; window.__rafN = 0;
window.requestAnimationFrame = function (cb) { window.__rafN++; window.__rafQ.push(cb); return window.__rafN; };
</script>
"""

_PRESS = r"""<script>(function () {
  var play = document.getElementById("gm-play");
  for (var i = 0; i < 5; i++) { var q = window.__rafQ; window.__rafQ = []; q.forEach(function (cb) { cb(1000 + 16 * i); }); }
  var before = window.__rafN;
  play.click();
  for (var j = 0; j < 20; j++) { var q2 = window.__rafQ; window.__rafQ = []; q2.forEach(function (cb) { cb(2000 + 16 * j); }); }
  var pre = document.createElement("pre"); pre.id = "gm-probe";
  pre.textContent = JSON.stringify({ disabled: play.disabled, label: play.textContent, pending: window.__rafQ.length,
                                     requested: window.__rafN - before });
  document.body.appendChild(pre);
})();</script>"""


def test_one_state_trajectory_does_not_loop(chrome):
    d, p, traj, R = _rp()
    one = copy.copy(traj)
    one.t, one.state, one.u = traj.t[:1], traj.state[:1], traj.u[:1]
    scene = Scene(d.layout)
    scene.add(RigidBody(R, ZeroVector, Box(0.1, 0.04, 1.0)))
    page = Report("r", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene,
                  trajectory=one).html()  # fmt: skip
    out = chrome(_at_end(_before_main(page, _PLAY), when_ready(_PRESS)))
    assert out["disabled"] is True and out["label"] == "Play"
    assert out["pending"] == 0 and out["requested"] <= 1  # no requestAnimationFrame loop


_TINTS = r"""<script>
(function () {
  var out = {};
  try {
    var sc = window.__gmScene, world = sc.children.filter(function (c) { return c.isGroup; })[0];
    var groups = world.children.filter(function (c) { return c.isGroup && c.matrixAutoUpdate === false; });
    out.bodies = groups.map(function (g) {
      var mesh = g.children[0];
      var mats = Array.isArray(mesh.material) ? mesh.material : [mesh.material];
      return { faces: mats.map(function (m) { return m.color.getHexString(); }),
               markers: g.children.filter(function (o) { return o.userData.gmMarker !== undefined; }).map(function (m) { return m.material.color.getHexString(); }) };
    });
  } catch (e) { out.error = String(e); }
  var pre = document.createElement("pre"); pre.id = "gm-probe"; pre.textContent = JSON.stringify(out);
  document.body.appendChild(pre);
})();
</script>
"""


def _lab(hexcolor: str) -> np.ndarray:
    c = np.array([int(hexcolor[i : i + 2], 16) / 255 for i in (0, 2, 4)])
    lin = np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    M = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])
    xyz = M @ lin / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > 216 / 24389, np.cbrt(xyz), (24389 / 27 * xyz + 16) / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def _de(a: str, b: str) -> float:
    return float(np.linalg.norm(_lab(a) - _lab(b)))


@pytest.mark.parametrize("theme", ["light", "dark"])
def test_orientation_tints_stand_out_on_every_palette_colour(chrome, theme):
    d, p, traj, R = _rp()
    scene = Scene(d.layout)
    for i, color in enumerate(PALETTE):
        scene.add(RigidBody(R, [i, 0, 0], Box(0.3, 0.2, 0.1), color=color, axes=False))
        scene.add(RigidBody(R, [i, 2, 0], Cylinder(0.1, 0.2), color=color, axes=False))
        scene.add(RigidBody(R, [i, 4, 0], Sphere(0.1), color=color, axes=False))
    page = Report("t", system=_sys(d), eom=d.eom, standard_form=d.sf, scene=scene,
                  trajectory=traj).html()  # fmt: skip
    page = page.replace('<html lang="en">', f'<html lang="en" data-theme="{theme}">')
    out = chrome(_at_end(_before_main(page, _PRE), when_ready(_TINTS)))
    assert "error" not in out
    worst = 1e9
    for b in out["bodies"]:
        base = b["faces"][1] if len(b["faces"]) > 1 else b["faces"][0]
        tints = ([b["faces"][0], b["faces"][2]] if len(b["faces"]) == 6 else []) + b["markers"]
        assert tints, b
        for t in tints:
            worst = min(worst, _de(t, base))
    assert worst >= 24.5, worst  # CIE76; the review measured ΔE 3.5 for a load-coloured box
