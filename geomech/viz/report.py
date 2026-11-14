"""A self-contained HTML report: 3-D animation, equations, standard form and checks.

.. code-block:: python

    report = Report("Double rigid pendulum", system=(L, dW, variables, inputs), eom=eom,
                    standard_form=sf, scene=scene, trajectory=traj, params=p)
    report.write("double_pendulum.html")

The page loads three.js r147 and MathJax 3.2.2 from pinned CDN URLs (the only external
resources) and inlines everything else, the trajectory included, as JSON.
"""

from __future__ import annotations

import html
import json
import re
import time
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import TYPE_CHECKING, Any, Mapping, Sequence

import numpy as np

from geomech.viz.latex import latex_names

if TYPE_CHECKING:
    from geomech.simulation import Trajectory
    from geomech.viz.scene import Scene

__all__ = ["Report", "EXTERNAL_URLS", "CHECK_EOM_TOL"]

THREE_URL = "https://cdn.jsdelivr.net/npm/three@0.147.0/build/three.min.js"
ORBIT_URL = "https://cdn.jsdelivr.net/npm/three@0.147.0/examples/js/controls/OrbitControls.js"
MATHJAX_URL = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-svg-full.min.js"
#: The only external resources a report loads.
EXTERNAL_URLS = (THREE_URL, ORBIT_URL, MATHJAX_URL)

#: ``check_eom`` errors above this are shown as a failure (``docs/numeric-verification.md``).
CHECK_EOM_TOL = 1e-5

#: Digits kept in the inlined trajectory (positions in scene units, rotation entries).
_DIGITS = 5

_MARKER = re.compile(r"__GM_([A-Z_]+)__")


def _template() -> str:
    return resources.files("geomech.viz").joinpath("template.html").read_text(encoding="utf-8")


def _fill(template: str, values: Mapping[str, str]) -> str:
    """Replace the ``__GM_NAME__`` markers of *template* in one pass (a value is never
    scanned again, so user text that contains a marker is inserted verbatim)."""

    def sub(m):
        return values[m.group(1)]

    return _MARKER.sub(sub, template)


def _json_script(obj) -> str:
    """JSON safe inside ``<script type="application/json">``."""
    text = json.dumps(obj, separators=(",", ":"), allow_nan=False)
    return text.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")


def _flat(arr) -> list:
    a = np.round(np.asarray(arr, dtype=float).reshape(-1), _DIGITS)
    a[a == 0] = 0.0  # no "-0.0"
    return a.tolist()


# ---------------------------------------------------------------------------
# LaTeX assembly
# ---------------------------------------------------------------------------


def _terms(expr) -> list:
    from geomech.core.operations.addition import Add, MAdd, VAdd

    if isinstance(expr, (Add, VAdd, MAdd)):
        return list(expr.nodes)
    return [expr]


def _term_lines(expr, names) -> list[str]:
    """One LaTeX string per addend, each with its sign (``+ a``, ``-b``) except the first."""
    from geomech.utils.printing import _strip_outer_parens, to_latex

    lines = []
    for i, term in enumerate(_terms(expr)):
        s = to_latex(term, names=names)
        if len(_terms(expr)) == 1:
            s = _strip_outer_parens(s)
        if i > 0 and not s.startswith("-"):
            s = "+ " + s
        lines.append(s)
    return lines


def _aligned(lhs: str, lines: Sequence[str], rhs: str = "") -> str:
    r"""``lhs = term₁ \\ + term₂ \\ …`` (or ``term₁ … = rhs``), one addend per line."""
    if not lines:
        lines = ["0"]
    rows = []
    for i, line in enumerate(lines):
        lead = f"{lhs} = {{}}& " if (i == 0 and lhs) else "& "
        if i == 0 and not lhs:
            lead = "& "
        rows.append(lead + line)
    if rhs:
        rows[-1] += f" = {rhs}"
    return "\\begin{aligned}\n" + " \\\\\n".join(rows) + "\n\\end{aligned}"


def _display(tex: str) -> str:
    return f'<div class="eq math">\\[{html.escape(tex, quote=False)}\\]</div>'


def _inline(tex: str) -> str:
    return f'<span class="math">\\({html.escape(tex, quote=False)}\\)</span>'


def _bmatrix(rows: Sequence[Sequence[str]]) -> str:
    body = " \\\\ ".join(" & ".join(r) for r in rows)
    return f"\\begin{{bmatrix}} {body} \\end{{bmatrix}}"


def _number(x: float) -> str:
    return f"{float(x):.6g}"


def _value_tex(val) -> str:
    arr = np.asarray(val, dtype=float)
    if arr.size == 1:
        return _number(arr.reshape(-1)[0])
    if arr.size == 3:
        return _bmatrix([[_number(v)] for v in arr.reshape(-1)])
    if arr.size == 9:
        return _bmatrix([[_number(v) for v in row] for row in arr.reshape(3, 3)])
    return html.escape(np.array2string(arr, precision=6))


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------


@dataclass
class _Stat:
    label: str
    value: str
    status: str = ""  # "", "ok" or "bad"
    note: str = ""


class Report:
    """The derivation, verification and simulation of one system as an HTML page.

    Args:
        title: the page title.
        system: ``(L, dW, variables, inputs)``, as passed to ``compute_eom`` and
            ``to_standard_form``.
        eom: the ``compute_eom`` result; derived (and timed) here when ``None``.
        standard_form: the ``to_standard_form`` result; computed here when ``None``.
        scene: a ``Scene``; with *trajectory*, the 3-D view is drawn.
        trajectory: a ``Trajectory`` from ``geomech.simulation.simulate``.
        params: the parameter values (``p``); defaults to ``trajectory.p``.
        checks: run ``check_eom`` and show its error (``True`` by default).
        units: ``{parameter name: unit}`` for the parameter table.
        source: the system definition as Python source, shown at the end.
        derive_time: the time ``compute_eom`` took, in seconds, when *eom* is given.
        fps: frames per second of simulated time kept in the animation (states are
            picked from the trajectory, never interpolated).
        description: a sentence under the title.

    ``write(path)`` writes the page; in Jupyter the report displays itself in an
    ``<iframe>`` (``_repr_html_``).
    """

    def __init__(
        self,
        title: str,
        system: tuple,
        eom: dict | None = None,
        standard_form: dict | None = None,
        scene: Scene | None = None,
        trajectory: Trajectory | None = None,
        params: dict | None = None,
        checks: bool = True,
        *,
        units: Mapping[str, str] | None = None,
        source: str | None = None,
        derive_time: float | None = None,
        fps: float = 60.0,
        description: str = "",
    ):
        from geomech import compute_eom, to_standard_form

        if not isinstance(system, (tuple, list)) or len(system) != 4:
            raise TypeError("Report: system must be (L, dW, variables, inputs)")
        if not (fps > 0 and np.isfinite(fps)):
            raise ValueError(f"Report: fps must be positive, got {fps}")
        if (scene is None) != (trajectory is None):
            raise ValueError("Report: pass both scene and trajectory, or neither")
        self.title = str(title)
        self.description = str(description)
        self.L, self.dW, self.variables, inputs = system
        self.inputs = list(inputs or [])
        if eom is None:
            t0 = time.perf_counter()
            eom = compute_eom(self.L, self.dW, self.variables)
            derive_time = time.perf_counter() - t0
        self.eom = eom
        self.derive_time = derive_time
        if standard_form is None:
            standard_form = to_standard_form(eom, self.variables, self.inputs)
        self.standard_form = standard_form
        self.scene: Any = scene
        self.trajectory: Any = trajectory
        if params is None and trajectory is not None:
            params = trajectory.p
        self.params = dict(params or {})
        self.units = dict(units or {})
        self.source = source
        self.fps = float(fps)
        self.check_error: float | None = None
        self.check_time: float | None = None
        if checks:
            from geomech.numeric import check_eom

            t0 = time.perf_counter()
            self.check_error = float(
                check_eom(self.L, self.dW, self.variables, self.eom, inputs=self.inputs)
            )
            self.check_time = time.perf_counter() - t0
        self.names = latex_names(
            self.L,
            self.dW,
            self.eom,
            self.standard_form,
            list(self.variables.scalars),
            list(self.variables.vectors),
            list(self.variables.matrices),
            self.inputs,
        )
        self._html: str | None = None

    # ------------------------------------------------------------------ pieces

    def _tex(self, expr) -> str:
        from geomech.utils.printing import to_latex

        return to_latex(expr, names=self.names)

    def _config_variables(self) -> list:
        v = self.variables
        return [*v.scalars, *v.vectors, *v.matrices]

    def _symbol_tex(self, key: str) -> str:
        """LaTeX of an acceleration or input key of the standard form."""
        sym = self._lookup(key)
        return self._tex(sym) if sym is not None else html.escape(key)

    def _lookup(self, key: str):
        from geomech.codegen.layout import _kinematic_symbols

        for var in self._config_variables():
            _, vel, acc = _kinematic_symbols(var)
            for sym in (var, vel, acc):
                if str(sym) == key:
                    return sym
        for u in self.inputs:
            if str(u) == key:
                return u
        return None

    def stats(self) -> list[_Stat]:
        """The numbers of the summary strip."""
        out = []
        if self.derive_time is not None:
            out.append(_Stat("Derivation", f"{self.derive_time:.3g} s", note="compute_eom"))
        else:
            out.append(_Stat("Derivation", "—", note="compute_eom (not timed)"))
        if self.check_error is None:
            out.append(_Stat("check_eom", "not run", "bad"))
        else:
            ok = self.check_error <= CHECK_EOM_TOL
            out.append(
                _Stat(
                    "check_eom error",
                    f"{self.check_error:.2e}",
                    "ok" if ok else "bad",
                    note=f"{'passes' if ok else 'FAILS'} (tolerance {CHECK_EOM_TOL:g})",
                )
            )
        traj = self.trajectory
        if traj is not None:
            drift = self.energy_drift()
            if drift is None:
                out.append(_Stat("Energy drift", "—", note="module has no energy()"))
            else:
                out.append(_Stat("Energy drift", f"{drift:.2e}", note="max |E − E₀| / |E₀|"))
            out.append(
                _Stat("Constraint drift", f"{traj.max_drift():.1e}", note="|q|−1, ‖RᵀR − I‖")
            )
            dt = traj.info.get("dt")
            step = f"dt = {dt:.3g} s" if dt is not None else ""
            out.append(_Stat("Integrator", traj.method, note=step))
            t = np.asarray(traj.t)
            out.append(_Stat("Simulated", f"{t[-1] - t[0]:.3g} s", note=f"{len(t)} stored states"))
        return out

    def energy_drift(self) -> float | None:
        """``max |E(t) − E(0)| / |E(0)|`` of ``trajectory.energy()`` (``None`` without an
        energy function); the absolute drift when ``E(0) = 0``."""
        try:
            E = np.asarray(self.trajectory.energy(), dtype=float)
        except TypeError:
            return None
        scale = abs(E[0]) if abs(E[0]) > 0 else 1.0
        return float(np.max(np.abs(E - E[0])) / scale)

    def frame_indices(self) -> np.ndarray:
        """The stored states kept in the animation: the first state at or after each
        multiple of ``1/fps`` of simulated time, and the last state."""
        t = np.asarray(self.trajectory.t, dtype=float)
        if len(t) < 2:
            return np.arange(len(t))
        targets = np.arange(t[0], t[-1], 1.0 / self.fps)
        idx = np.searchsorted(t, targets - 1e-12)
        return np.unique(np.r_[idx[idx < len(t)], len(t) - 1])

    def scene_data(self) -> dict:
        """The animation data inlined in the page (rounded to 5 decimals)."""
        scene = self.scene
        idx = self.frame_indices()
        fr = scene.frames(self.trajectory, indices=idx)
        rot_tracks: list[list] = []
        rot_keys: dict[str, int] = {}
        bodies = []
        extent = 0.0
        for body, ev in zip(scene.bodies, fr.bodies):
            key = str(body.R)
            if key not in rot_keys:
                rot_keys[key] = len(rot_tracks)
                rot_tracks.append(_flat(ev["rotation"]))
            shape = body.shape
            extent = max(extent, shape.extent)
            bodies.append(
                {
                    "name": body.name,
                    "shape": shape.to_dict(),
                    "color": body.color,
                    "axes": bool(body.axes),
                    "axisLengths": _axis_lengths(shape),
                    "rot": rot_keys[key],
                    "pos": _flat(ev["center"]),
                }
            )
        lo, hi = fr.bounds()
        center = (lo + hi) / 2
        radius = float(max(np.max(hi - lo) / 2 + extent / 2, 0.5))
        zlow = float(lo[2] - extent / 2)
        ground = 0.0 if 0.0 <= zlow <= 2 * radius else zlow - 0.05 * radius
        return {
            "t": _flat(fr.t - fr.t[0]),
            "t0": float(fr.t[0]),
            "rotations": rot_tracks,
            "bodies": bodies,
            "points": [
                {"name": pm.name, "radius": pm.radius, "color": pm.color, "pos": _flat(v)}
                for pm, v in zip(scene.points, fr.points)
            ],
            "rods": [
                {"name": r.name, "radius": r.radius, "color": r.color, "ends": _flat(v)}
                for r, v in zip(scene.rods, fr.rods)
            ],
            "traces": [
                {"name": tr.name, "color": tr.color, "pos": _flat(v)}
                for tr, v in zip(scene.traces, fr.traces)
            ],
            "center": [float(c) for c in center],
            "radius": radius,
            "ground": float(ground),
        }

    # ------------------------------------------------------------------ sections

    def _section_summary(self) -> str:
        cards = []
        for s in self.stats():
            cls = f"stat {s.status}".strip()
            note = f'<div class="stat-note">{html.escape(s.note)}</div>' if s.note else ""
            cards.append(
                f'<div class="{cls}"><div class="stat-label">{html.escape(s.label)}</div>'
                f'<div class="stat-value">{html.escape(s.value)}</div>{note}</div>'
            )
        return f'<section id="summary" class="summary">{"".join(cards)}</section>'

    def _section_view(self) -> str:
        if self.scene is None:
            return ""
        return """<section id="view">
<h2>Animation</h2>
<div class="viewer-wrap">
  <div id="gm-viewer" class="viewer" aria-label="3-D view of the trajectory">
    <div id="gm-viewer-msg" class="viewer-msg">Loading three.js…</div>
  </div>
  <div class="controls" role="group" aria-label="Playback">
    <button id="gm-play" type="button">Pause</button>
    <input id="gm-scrub" type="range" min="0" max="0" value="0" step="1"
           aria-label="Time">
    <span id="gm-time" class="time">t = 0.00 s</span>
  </div>
  <div class="controls" role="group" aria-label="View options">
    <label>Speed <select id="gm-speed">
      <option value="0.1">0.1×</option><option value="0.25">0.25×</option>
      <option value="0.5">0.5×</option><option value="1" selected>1×</option>
      <option value="2">2×</option><option value="4">4×</option></select></label>
    <button id="gm-reset" type="button">Reset camera</button>
    <label><input id="gm-axes" type="checkbox" checked> Body axes</label>
    <label><input id="gm-traces" type="checkbox" checked> Traces</label>
    <button id="gm-theme" type="button" title="Switch theme">Theme: auto</button>
  </div>
  <p class="legend">Inertial frame: e₁ <span class="sw x"></span> e₂ <span class="sw y"></span>
  e₃ <span class="sw z"></span> (e₃ up, gravity along −e₃). Body frames use the same
  colours. Drag to orbit, scroll to zoom, right-drag to pan.</p>
</div>
</section>"""

    def _section_lagrangian(self) -> str:
        tex = _aligned("L", _term_lines(self.L, self.names))
        parts = [f'<section id="lagrangian"><h2>Lagrangian</h2>{_display(tex)}']
        if self.dW is not None and hasattr(self.dW, "type"):
            dw = _aligned("\\delta W", _term_lines(self.dW, self.names))
            parts.append(f"<p>Virtual work:</p>{_display(dw)}")
        parts.append("</section>")
        return "".join(parts)

    def _variation_label(self, key: str) -> str:
        for var in self._config_variables():
            if str(var.get_variation_vector()) == key:
                return (
                    f"Equation of {_inline(self._tex(var))}, coefficient of "
                    f"{_inline(self._tex(var.get_variation_vector()))}"
                )
        return html.escape(key)

    def _section_eom(self) -> str:
        parts = ['<section id="eom"><h2>Equations of motion</h2>']
        parts.append(
            '<p class="note">E = −(Euler–Lagrange) + forcing; each equation is E = 0, one term'
            " per line. S² equations hold on the tangent plane of their point.</p>"
        )
        for key, (_, eq) in self.eom.items():
            parts.append(f"<h3>{self._variation_label(key)}</h3>")
            parts.append(_display(_aligned("0", _term_lines(eq, self.names))))
        parts.append("</section>")
        return "".join(parts)

    def _section_standard_form(self) -> str:
        sf = self.standard_form
        rows = list(sf.keys())
        acc_keys: list[str] = []
        inp_keys: list[str] = []
        for eq in sf.values():
            for k in eq.M:
                if k not in acc_keys:
                    acc_keys.append(k)
            for k in eq.G:
                if k not in inp_keys:
                    inp_keys.append(k)
        # columns in variable / input order where known
        order = {str(acc): i for i, acc in enumerate(self._accelerations())}
        acc_keys.sort(key=lambda k: order.get(k, len(order)))
        iorder = {str(u): i for i, u in enumerate(self.inputs)}
        inp_keys.sort(key=lambda k: iorder.get(k, len(iorder)))
        n = len(rows)
        sep = "" if max(n, len(acc_keys), len(inp_keys)) < 10 else ","

        def sub(*ij):
            return "_{" + sep.join(str(i + 1) for i in ij) + "}"

        M_over = [
            [f"\\mathbf{{M}}{sub(i, j)}" if a in sf[r].M else "0" for j, a in enumerate(acc_keys)]
            for i, r in enumerate(rows)
        ]
        overview = _bmatrix(M_over) + " " + _bmatrix([[self._symbol_tex(a)] for a in acc_keys])
        overview += " + " + _bmatrix([[f"\\mathbf{{f}}{sub(i)}"] for i in range(n)])
        if inp_keys:
            G_over = [
                [
                    f"\\mathbf{{G}}{sub(i, k)}" if u in sf[r].G else "0"
                    for k, u in enumerate(inp_keys)
                ]
                for i, r in enumerate(rows)
            ]
            overview += " + " + _bmatrix(G_over) + " "
            overview += _bmatrix([[self._symbol_tex(u)] for u in inp_keys])
        overview += " = 0"
        parts = ['<section id="standard-form"><h2>Standard form</h2>']
        parts.append(
            '<p class="note">M a + f + G u = 0 with M = −(physical inertia). Rows are the'
            " equations above, in the same order; blocks that are not shown are zero.</p>"
        )
        parts.append(_display(overview))
        parts.append("<h3>Inertia blocks</h3>")
        for i, r in enumerate(rows):
            for j, a in enumerate(acc_keys):
                if a in sf[r].M:
                    lhs = f"\\mathbf{{M}}{sub(i, j)}"
                    parts.append(_display(_aligned(lhs, _term_lines(sf[r].M[a], self.names))))
        parts.append("<h3>Bias terms</h3>")
        for i, r in enumerate(rows):
            lhs = f"\\mathbf{{f}}{sub(i)}"
            parts.append(_display(_aligned(lhs, _term_lines(sf[r].f, self.names))))
        if inp_keys:
            parts.append("<h3>Input blocks</h3>")
            for i, r in enumerate(rows):
                for k, u in enumerate(inp_keys):
                    if u in sf[r].G:
                        lhs = f"\\mathbf{{G}}{sub(i, k)}"
                        parts.append(_display(_aligned(lhs, _term_lines(sf[r].G[u], self.names))))
        parts.append("</section>")
        return "".join(parts)

    def _accelerations(self) -> list:
        from geomech.codegen.layout import _kinematic_symbols

        return [_kinematic_symbols(v)[2] for v in self._config_variables()]

    def _section_params(self) -> str:
        if not self.params:
            return ""
        rows = []
        for name in sorted(self.params):
            pretty = self.names.get(name, name)
            unit = html.escape(self.units.get(name, ""))
            rows.append(
                f"<tr><td>{_inline(pretty)}</td><td>{_inline(_value_tex(self.params[name]))}"
                f"</td><td>{unit}</td></tr>"
            )
        return (
            '<section id="parameters"><h2>Parameters</h2><div class="table-wrap"><table>'
            "<thead><tr><th>Name</th><th>Value</th><th>Units</th></tr></thead><tbody>"
            + "".join(rows)
            + "</tbody></table></div></section>"
        )

    def _section_source(self) -> str:
        if not self.source:
            return ""
        return (
            '<section id="source"><h2>System definition</h2>'
            f"<pre><code>{html.escape(self.source.strip())}</code></pre></section>"
        )

    # ------------------------------------------------------------------ output

    def html(self) -> str:
        """The page as one string (cached)."""
        if self._html is not None:
            return self._html
        body = "\n".join(
            s
            for s in (
                self._section_summary(),
                self._section_view(),
                self._section_lagrangian(),
                self._section_eom(),
                self._section_standard_form(),
                self._section_params(),
                self._section_source(),
            )
            if s
        )
        data: dict[str, Any] = {"title": self.title}
        if self.scene is not None:
            data["scene"] = self.scene_data()
        subtitle = html.escape(self.description) if self.description else ""
        self._html = _fill(
            _template(),
            {
                "TITLE": html.escape(self.title),
                "SUBTITLE": subtitle,
                "BODY": body,
                "DATA": _json_script(data),
                "THREE_URL": THREE_URL,
                "ORBIT_URL": ORBIT_URL,
                "MATHJAX_URL": MATHJAX_URL,
            },
        )
        return self._html

    def write(self, path) -> Path:
        """Write the page to *path* (UTF-8) and return the path."""
        path = Path(path)
        path.write_text(self.html(), encoding="utf-8")
        return path

    def _repr_html_(self) -> str:
        """Jupyter: the page in an ``<iframe srcdoc>``, so that its scripts and MathJax
        configuration do not interfere with the notebook's."""
        doc = html.escape(self.html(), quote=True)
        return (
            f'<iframe srcdoc="{doc}" style="width:100%;height:900px;border:0"'
            f' title="{html.escape(self.title)}"></iframe>'
        )


def _axis_lengths(shape) -> list[float]:
    """Body-axis arrow lengths: half the shape along each axis plus a margin, so every
    arrow sticks out of the shape."""
    from geomech.viz.scene import Box, Cylinder

    if isinstance(shape, Box):
        half = [shape.lx / 2, shape.ly / 2, shape.lz / 2]
        small = sorted([shape.lx, shape.ly, shape.lz])[0]
    elif isinstance(shape, Cylinder):
        half = [shape.radius] * 3
        half["xyz".index(shape.axis)] = shape.length / 2
        small = min(2 * shape.radius, shape.length)
    else:
        half = [shape.radius] * 3
        small = 2 * shape.radius
    margin = max(0.25 * shape.extent, 4 * small)
    return [round(h + margin, 6) for h in half]
